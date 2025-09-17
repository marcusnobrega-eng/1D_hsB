%% ========================================================================
%  Recharge–Release Test with k,f Calibration via fmincon (min RMSE)
%  ------------------------------------------------------------------------
%  Author: Marcus Nobrega Gomes Junior, PhD  (extended with calibration)
%  Affiliation: University of Arizona, Dept. of Hydrology & Atmospheric Sciences
%  Date: Sep/2025
%
%  PURPOSE
%  --------
%  (A) Set up a recharge–release experiment and compare hsB recession to
%      a linear-reservoir target Q(t) = a * exp(-b t).
%  (B) Calibrate two hsB parameters — hydraulic parameter k [m/s] and
%      drainable porosity f [-] — using fmincon to MINIMIZE the RMSE
%      (in mm/day) between hsB recession and the linear reservoir.
%
%  HOW CALIBRATION WORKS
%  ---------------------
%   * Decision vars: p = [k, f].
%   * Objective: RMSE_mmday( p ) from compare_recession_LR_days(…).
%   * Constraints: simple bounds lb <= p <= ub (no nonlinear constraints).
%   * We recommend *sensible* bounds (see CAL.lb/ub below) to avoid
%     unphysical regions and help the optimizer.
%
%  TUNING fmincon (rationales inline)
%  ----------------------------------
%   - Algorithm:
%       'interior-point' (robust for bound-constrained smooth problems),
%       'sqp' (good for small problems), or 'active-set' (legacy).
%   - Tolerances:
%       * 'OptimalityTolerance' ~ 1e-6 to 1e-8 for steady results.
%       * 'StepTolerance' tight enough (e.g., 1e-10…1e-8) to avoid
%         premature stopping on flat plateaus.
%       * 'MaxFunctionEvaluations' set high enough (100–200+) since each
%         eval runs the hsB solver.
%   - Scaling:
%       Use 'TypicalX' to tell fmincon the *scale* of k and f. This improves
%       line searches and step sizes (e.g., TypicalX = [1e-4, 0.05]).
%   - Parallel:
%       'UseParallel' can speed finite-difference gradients **if** you have
%       Parallel Computing Toolbox. Start a pool before running if desired.
%
%% ========================================================================

clear variables; close all; clc;
addpath('hsB_catchment')

if exist('fmincon','file') ~= 2
    error('This script requires Optimization Toolbox (fmincon).');
end

%% ---------------- 0) Inputs (xH, H, metrics) ----------------
% Requires xH, H, and metrics (for hillslope length) in 'input_data.mat'
%  - xH: centers of distance-to-stream bins [m]
%  - H : width per bin [m]
%  - metrics.max_hillslope_length_m
% The file should also include (if you use it) hillslope angle, etc.
% Alternatively, you can run a pre-processing algorithm available here
% This script must produce: xH [1×n], H [1×n], Forcing (struct array), metrics (struct)
%  - xH: centers of uniform distance-to-stream bins (m)
%  - H : hillslope width per bin (m) such that area ≈ sum(H)*dxH
%  - metrics.max_hillslope_length_m: max flow distance to stream (m)
% run('hsB_PreProcessing.m')
% save('input_data','metrics','xH','H','Forcing','t_forcing_days','t0_datetime');
load('input_data.mat');   % provides: xH, H, metrics

% --- FORCE DOUBLE on inputs that might have been saved as single ---
xH     = double(xH);
H      = double(H);
metrics.max_hillslope_length_m = double(metrics.max_hillslope_length_m);
if isfield(metrics,'hillslope_angle_deg')
    metrics.hillslope_angle_deg = double(metrics.hillslope_angle_deg);
end

%% ---------------- User settings: scenario + LR params ----------
% Constant recharge "on", then "off" to produce recession
recharge_rate_mmday = (11.6518/283941248)*1000*86400;    % [mm/day]
N_recharge_days     = 30;   % days with recharge on
N_release_days      = 60;   % days with recharge off (recession window)
dt_forcing_hours    = 24;   % forcing stamp interval (hrs)

% Linear reservoir target (used for the RMSE objective)
%   Q_LR(t_days) = a_LR * exp(-b_LR * t_days)
a_LR  = 11.6518;     % [m^3/s] amplitude at recession start
b_LR  = 0.353;       % [1/day]  decay rate

%% ---------------- 1) Grid + ghost-padded width ----------------
xH  = xH(:).';  H = H(:).';
if numel(xH) < 2, error('xH must have at least two bins.'); end

dxH_vec = diff(xH);
dxH = median(dxH_vec);
if any(abs(dxH_vec - dxH) > 1e-6*dxH)
    warning('xH not perfectly uniform; using median(diff(xH)).');
end

L = metrics.max_hillslope_length_m;   % domain length [m]

% Clustered + uniform grid (simple defaults)
Nx_target = numel(xH);
M  = 3;  r = 1.15;  d0 = 5;                  % near-outlet clustering
dx_geo  = d0 * r.^(0:M-1);
len_geo = sum(dx_geo);
if len_geo >= L, error('Clustered block (%.2f m) >= L=%.2f m.', len_geo, L); end
Nx = Nx_target;  if Nx <= M, error('Nx must be > M.'); end
dx_uni = (L - len_geo) / (Nx - M);
dx_all = [dx_geo, dx_uni * ones(1, Nx-M)];
xe = [0; cumsum(dx_all(:))];
x  = 0.5 * (xe(1:end-1) + xe(2:end));

% Ghost-pad and map widths smoothly
[xHg, Hg, dxH_check, ~] = add_ghosts_to_width_function(xH, H);
if abs(dxH_check - dxH) > 1e-9*dxH, dxH = dxH_check; end
w = interp1(xHg(:), Hg(:), x(:), 'pchip', 'extrap');  % [Nx×1]

% Area closure
A_hist = sum(H) * dxH;
A_grid = sum(w(:) .* diff(xe));
if A_grid <= 0, error('A_grid <= 0.'); end
w = w * (A_hist / A_grid);
fprintf('Area closure OK. A_hist=%.3e, A_grid(adj)=%.3e\n', A_hist, sum(w(:).*diff(xe)));

% --- Hillslope area (m^2) for mm/day conversion ---
dx   = diff(xe(:));          % [Nx×1] cell widths [m]
A_hs = sum(w(:) .* dx);      % total hillslope area [m^2]

%% ---------------- 2) Forcing: constant recharge, then zero ----
t_total_days = N_recharge_days + N_release_days;
dt_days      = dt_forcing_hours/24;
t_days       = (0:dt_days:t_total_days).';              % [nF×1]
N_mmday      = zeros(size(t_days));
N_mmday(t_days < N_recharge_days) = recharge_rate_mmday; % recharge ON then OFF

% Convert to SI (m/s)
N_mps = (N_mmday/1000) / 86400;

% Package forcing for hsB solver
ForcingNative = struct();
ForcingNative.t_days = t_days;    % [nF×1]
ForcingNative.N_mps  = N_mps;     % [nF×1] (uniform in space)

%% ---------------- 3) Physical parameters ----------------------
% Initial (pre-calibration) guesses
params = struct();
params.k    = 6e-4;                               % [m/s]  initial guess
params.f    = 0.05;                               % [-]    initial guess
params.iota = metrics.hillslope_angle_deg*pi/180; % [rad]
params.S0   = zeros(numel(x),1);                  % start dry (will fill)
params.D    = 0.4;                                % [m]

%% ---------------- 4) Solver options ---------------------------
opts = struct();
opts.theta       = 1.0;
opts.omega       = 0.5;
opts.picard_max  = 20;
opts.picard_tol  = 1e-8;
opts.safeguard   = true;
opts.mass_check  = true;
opts.plot_live   = false;
opts.verbose     = false;

opts.substep_mode    = 'adaptive';
opts.max_sub_dt      = 24*3600;
opts.min_sub_dt      = 300;
opts.cfl_adv         = 0.1;
opts.cfl_diff        = 0.1;
opts.picard_easy     = 4;
opts.picard_hard     = 10;
opts.frac_change_eta = 0.4;
opts.fixed_nsub      = 4;
opts.no_recharge_outlet = false;

% No final figs from the solver during optimization; we’ll plot once at end
PLOT.final_figs = false; PLOT.stream_stamps = false;
opts.stream_plot_stamps = PLOT.stream_stamps; do_plot = PLOT.final_figs;
sty = struct('fs',13,'axlw',2,'plotlw',2,'tickdir','in','ticklength',[0.02 0.02]);

%% ---------------- 5) CALIBRATION SETTINGS (SIMPLE) --------------
CAL.do_calibrate = true;

% Work in log10-space for better scaling:
% q = [log10(k), log10(f)]
CAL.lb_lin = [1e-7;  1e-3];           % physical lower bounds
CAL.ub_lin = [1e-2;  0.5];            % physical upper bounds

CAL.q0   = log10([params.k; params.f]);
CAL.qlb  = log10(CAL.lb_lin);
CAL.qub  = log10(CAL.ub_lin);

% Simple, robust fmincon options (no fragile/advanced knobs)
options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...            % good for small, smooth problems
    'Display','iter', ...
    'MaxFunctionEvaluations',300, ... % bump if needed
    'OptimalityTolerance',1e-6, ...
    'UseParallel',false);             % set true if you have a pool

% Notes on choosing these:
%  * Start with wider bounds; if the optimum sticks to a bound, widen or
%    revisit your prior. If exploring too wide, the solver may waste steps.
%  * If the surface is noisy, loosen OptTol/StepTol or use 'sqp'.
%  * If gradients are tricky, keep 'forward' differencing; 'central' is
%    more accurate but doubles model runs.
%  * TypicalX should be near expected magnitude of the optimum.

%% ---------------- 6) Run base model once (optional sanity) ----
[Qout0, ~, ~, ~, ~, diag0] = hsB_solver(x, w, ForcingNative, params, opts, do_plot, sty);
fprintf('\nBase run complete: forcing steps=%d, Nx=%d, internal steps=%d\n', ...
    numel(ForcingNative.t_days), numel(x), diag0.n_internal_steps);

% Compute initial metrics to see where we start
m0 = compare_recession_LR_days( ...
    ForcingNative.t_days, Qout0, ForcingNative.N_mps, ...
    a_LR, b_LR, false, A_hs, ...
    recharge_rate_mmday, N_recharge_days, N_release_days);
fprintf('Initial RMSE (mm/day): %.4g | NSE: %.3f | PBIAS: %.2f%%\n', ...
    m0.RMSE_mmday, m0.NSE, m0.PBIAS_percent);

%% ---------------- 7) fmincon calibration ----------------------
if CAL.do_calibrate
    fprintf('\n=== Starting fmincon calibration of [k, f] in log10-space (simple) ===\n');

    obj_q = @(q) obj_rmse_q( ...
        q, x, w, ForcingNative, params, opts, sty, ...
        A_hs, recharge_rate_mmday, N_recharge_days, N_release_days, ...
        a_LR, b_LR);

    A = []; b = []; Aeq = []; beq = []; nonlcon = [];

    [q_best, fval, exitflag, output] = fmincon( ...
        obj_q, CAL.q0(:), A, b, Aeq, beq, CAL.qlb(:), CAL.qub(:), nonlcon, options);

    fprintf('\n=== fmincon done (exitflag=%d) ===\n', exitflag); disp(output);

    % Map back to linear space
    p_best   = 10.^q_best(:);
    params.k = p_best(1);
    params.f = p_best(2);
    fprintf('Best params: k = %.6g [m/s], f = %.6g [-]\n', params.k, params.f);
    fprintf('Best objective (RMSE mm/day): %.6g\n', fval);
else
    fprintf('\nCalibration disabled; using initial params.\n');
end

%% ---------------- 8) Final run & plot with metrics -------------
[Qout, Qsurf, Qtotal, h, S, diag] = hsB_solver(x, w, ForcingNative, params, opts, false, sty);
fprintf('\nFinal run complete: forcing steps=%d, Nx=%d, internal steps=%d\n', ...
    numel(ForcingNative.t_days), numel(x), diag.n_internal_steps);

metrics_out = compare_recession_LR_days( ...
    ForcingNative.t_days, ...      % days
    Qout, ...                      % m^3/s
    ForcingNative.N_mps, ...       % m/s
    a_LR, b_LR, ...                % LR params (a: m^3/s, b: 1/day)
    true, ...                      % do_plot (final figure)
    A_hs, ...                      % m^2  (area)
    recharge_rate_mmday, ...       % mm/day
    N_recharge_days, ...           % days with recharge
    N_release_days, ...            % dry days
    params.k, params.f ...         % <<< show these on the plot
);

disp(metrics_out);

%% ========================================================================
%  Helper functions
% ========================================================================

function [xHg, Hg, dxH, Lphys] = add_ghosts_to_width_function(xH, H)
    xH = xH(:).'; H = H(:).';
    if numel(xH) < 2, error('xH must have ≥2 bins.'); end
    d = diff(xH);
    dxH = median(d);
    if any(abs(d - dxH) > 1e-6*dxH)
        warning('add_ghosts_to_width_function: xH non-uniform; using median dxH.');
    end
    Lphys = xH(end) + dxH/2;
    xLghost = -dxH/2;  xRghost = Lphys + dxH/2;
    HLghost = H(1);    HRghost = H(end);
    xHg = [xLghost, xH, xRghost];
    Hg  = [HLghost, H,  HRghost];
end

function metrics = compare_recession_LR_days( ...
    t_days, Q_hsB, N_mps, a, b, do_plot, A_hs, ...
    rec_mmday, n_recharge_days, n_dry_days, k_val, f_val)
% Compare hsB recession to Q(t)=a*exp(-b*t) and PLOT IN mm/day.
% Adds a LaTeX metrics box with recharge rate, recharge/dry days, RMSE, NSE, PBIAS.
    
    if nargin < 11, k_val = []; end
    if nargin < 12, f_val = []; end

    if nargin < 10
        error('Pass A_hs, rec_mmday, n_recharge_days, n_dry_days.');
    end
    if nargin < 6 || isempty(do_plot), do_plot = true; end

    % ----- Find start of recession (first zero recharge) -----
    k0 = find(N_mps <= 0, 1, 'first');
    if isempty(k0) || k0 == numel(t_days)
        error('No recession window found.');
    end

    % Time since recession start (days)
    t_rec_days = t_days(k0:end) - t_days(k0);

    % Observed hsB discharge (flip sign -> positive) [m^3/s]
    Q_obs_m3s = -Q_hsB(k0:end);

    % Linear reservoir [m^3/s]  (t in days)
    Q_sim_m3s = a * exp(-b * t_rec_days);

    % ----- Convert to mm/day using area -----
    m3s_to_mmday = (1000 * 86400) / A_hs;  % m^3/s -> mm/day
    Q_obs = Q_obs_m3s * m3s_to_mmday;
    Q_sim = Q_sim_m3s * m3s_to_mmday;

    % ----- Metrics -----
    rmse  = sqrt(mean((Q_sim - Q_obs).^2));
    denom = sum((Q_obs - mean(Q_obs)).^2);
    NSE   = NaN;  if denom > 0, NSE = 1 - sum((Q_sim - Q_obs).^2)/denom; end
    pbias = 100 * (sum(Q_sim - Q_obs) / sum(Q_obs));

    metrics = struct( ...
        'RMSE_mmday', rmse, 'NSE', NSE, 'PBIAS_percent', pbias, ...
        't_rec_days', t_rec_days, 'Q_obs_mmday', Q_obs, 'Q_sim_mmday', Q_sim, ...
        'a_m3s', a, 'b_dayinv', b, 'A_hs_m2', A_hs, ...
        'recharge_mmday', rec_mmday, 'n_recharge_days', n_recharge_days, 'n_dry_days', n_dry_days, ...
        'k0', k0, ...
        'k_mps', k_val, 'f', f_val );

    if ~do_plot, return; end

    % ----- Plot -----
    figure('Color','w','Position',[120 120 900 420]);
    plot(t_rec_days, Q_obs, '-',  'LineWidth', 2); hold on;
    plot(t_rec_days, Q_sim, '--', 'LineWidth', 2);
    grid on; box on;

    xlabel('Time since recession start (days)');
    ylabel('Discharge, $Q$ (mm/day)', 'Interpreter','latex');

    % Use \textsuperscript in title for units
    title(sprintf(['Recession: hsB vs. Q(t) = a exp(-b t)', ...
        '   (a = %.3g m$^{\\textsuperscript{3}}$/s, b = %.3g 1/day)'], a, b), ...
        'Interpreter','latex');

    legend({'hsB baseflow','Linear reservoir'}, 'Location','northeast');

    % LaTeX metrics box with one metric per line (use cell array -> one line/cell)
    ymax = max([Q_obs; Q_sim]); if ~isfinite(ymax) || ymax<=0, ymax = 1; end
    
    % ---------- Metrics box (LaTeX; one metric per line) ----------
    lines = {
        sprintf('\\textbf{Recharge:} %.3g\\,mm/day', rec_mmday)
        sprintf('\\textbf{Days on:} %d', n_recharge_days)
        sprintf('\\textbf{Dry days:} %d', n_dry_days)
    };
    
    % Optionally show calibrated/used parameters
    if ~isempty(k_val)
        lines{end+1} = sprintf('\\textbf{k:} %.3g\\,m/s', k_val);
    end
    if ~isempty(f_val)
        lines{end+1} = sprintf('\\textbf{f:} %.3f', f_val);
    end
    
    % Then the fit metrics
    lines = [lines; {
        sprintf('\\textbf{RMSE:} %.4g\\,mm/day', rmse)
        sprintf('\\textbf{NSE:} %.3f', NSE)
        sprintf('\\textbf{PBIAS:} %.2f\\,\\%%', pbias)
    }];

    % Add textbox annotation (normalized figure coords)
    annotation('textbox', [0.24 0.54 0.26 0.32], ...
        'String', lines, ...
        'Interpreter','latex', ...
        'FontSize', 11, ...
        'EdgeColor','k', ...
        'LineWidth', 1.0, ...
        'BackgroundColor','w', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FitBoxToText','on');
end

function rmse = obj_rmse_mmday(p, x, w, ForcingNative, base_params, opts, sty, ...
                               A_hs, rec_mmday, n_recharge_days, n_dry_days, ...
                               a_LR, b_LR)
% Objective wrapper for fmincon: return RMSE in mm/day for given [k, f].
% Safe-guards: if solver fails or returns non-finite, penalize heavily.

    % Enforce bounds softly (fmincon also enforces them)
    k = p(1); f = p(2);
    if ~isfinite(k) || ~isfinite(f) || k <= 0 || f <= 0
        rmse = 1e12; return;
    end

    % Set params and run model (no plotting inside objective)
    params = base_params;
    params.k = k;
    params.f = f;

    try
        do_plot = false;
        [Qout, ~, ~, ~, ~, diag] = hsB_solver(x, w, ForcingNative, params, opts, do_plot, sty); %#ok<ASGLU>

        % Compute metrics vs LR (no plot)
        M = compare_recession_LR_days( ...
            ForcingNative.t_days, Qout, ForcingNative.N_mps, ...
            a_LR, b_LR, false, A_hs, rec_mmday, n_recharge_days, n_dry_days);

        rmse = M.RMSE_mmday;

        % Guard against NaN/Inf
        if ~isfinite(rmse), rmse = 1e12; end

    catch ME
        % Penalize failures so optimizer avoids this region
        warning('obj_rmse_mmday: failure at k=%.3g, f=%.3g | %s', k, f, ME.message);
        rmse = 1e12;
    end
end

function rmse = obj_rmse_q(q, x, w, ForcingNative, base_params, opts, sty, ...
                           A_hs, rec_mmday, n_recharge_days, n_dry_days, ...
                           a_LR, b_LR)
% q = [log10(k), log10(f)]  ->  p = [k, f] in linear space

    q   = double(q(:));
    k   = 10.^q(1);
    f   = 10.^q(2);

    % Hard guards (fmincon enforces bounds, but be defensive)
    if ~isfinite(k) || ~isfinite(f) || k<=0 || f<=0
        rmse = 1e12; rmse = double(rmse); return;
    end

    % Evaluate the original (linear-space) objective
    rmse = obj_rmse_mmday([k; f], x, w, ForcingNative, base_params, opts, sty, ...
                          A_hs, rec_mmday, n_recharge_days, n_dry_days, ...
                          a_LR, b_LR);
    rmse = double(rmse);
    if ~isfinite(rmse), rmse = 1e12; end
end
