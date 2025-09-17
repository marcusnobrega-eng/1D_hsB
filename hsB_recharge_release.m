%% ========================================================================
%  Recharge–Release Test Script for hsB Model with Linear Reservoir Comparison
%  ------------------------------------------------------------------------
%  Author: Marcus Nobrega Gomes Junior, PhD
%  Affiliation: University of Arizona, Dept. of Hydrology & Atmospheric Sciences
%  Date: Sep/2025
%
%  PURPOSE
%  --------
%  This script sets up and runs a simple "recharge–release" scenario to test
%  the hillslope-storage Boussinesq (hsB) model. The goal is to validate the
%  simulated baseflow recession against an analytical linear reservoir model:
%
%       Q(t) = a * exp(-b * t)
%
%  where:
%     - Q(t)  = discharge [m³/s]
%     - a     = discharge amplitude at start of recession [m³/s]
%     - b     = exponential decay rate [1/day]
%     - t     = time since recharge stops [days]
%
%  WORKFLOW
%  --------
%   (1) Load preprocessed hillslope geometry and forcing (input_data.mat).
%   (2) Build a clustered grid with ghost-padded widths.
%   (3) Apply a constant recharge (mm/day) for a fixed number of days.
%   (4) Switch recharge to zero for a given number of "dry" days to trigger
%       a baseflow recession.
%   (5) Run the hsB solver with these inputs to obtain simulated outflow.
%   (6) Convert model discharges (m³/s) to equivalent depth fluxes (mm/day).
%   (7) Compare simulated recession to the analytical reservoir model.
%   (8) Compute performance metrics (RMSE, NSE, PBIAS).
%   (9) Plot hsB and linear reservoir hydrographs, with annotated scenario
%       parameters and fit metrics.
%
%  USER INPUTS
%  -----------
%   - recharge_rate_mmday : Recharge intensity [mm/day] during the "on" phase
%   - N_recharge_days     : Number of recharge days (constant recharge applied)
%   - N_release_days      : Number of dry days (recharge = 0, recession only)
%   - a_LR                : Linear reservoir amplitude [m³/s]
%   - b_LR                : Linear reservoir decay constant [1/day]
%   - dt_forcing_hours    : Time step between forcing stamps [hours]
%
%  INTERNALS
%  ---------
%   - hsB_solver: Solves storage/discharge dynamics on a 1D hillslope grid
%   - add_ghosts_to_width_function: Adds ghost bins for flux consistency
%   - compare_recession_LR_days: Post-processing function that converts model
%       outputs to mm/day, applies the linear reservoir, and computes metrics
%
%  OUTPUTS
%  -------
%   - Time series of hsB discharge (mm/day) and linear reservoir discharge
%   - Recession plot with scenario information and performance metrics
%   - Metrics reported: RMSE [mm/day], NSE [-], PBIAS [%]
%
%  METRICS DEFINITIONS
%  --------------------
%   - RMSE  (Root Mean Square Error): Magnitude of deviation [mm/day]
%   - NSE   (Nash–Sutcliffe Efficiency): Fit skill relative to mean of obs
%   - PBIAS (% Bias): Relative bias between observed and simulated volumes
%
%  NOTES
%  -----
%   * The observed "Q_obs" from hsB is taken as -Qout (to enforce positive flows).
%   * Hillslope discharge is normalized by the total contributing area to
%     express flows in depth-equivalent units (mm/day).
%   * Parameters a and b must be tuned manually to approximate the observed
%     recession slope.
%
%  EXAMPLE
%  -------
%     recharge_rate_mmday = 5;
%     N_recharge_days     = 30;
%     N_release_days      = 60;
%     a_LR = 0.1; b_LR = 0.05;   % test linear reservoir params
%     run('hsB_recharge_release_driver.m')
%
%  REFERENCES
%  ----------
%   - Troch, P. A., et al. (2003). "Catchment-scale hydrological modeling and
%     the hillslope-storage Boussinesq (hsB) model."
%   - McDonnell, J. J., et al. (2007). "Moving beyond heterogeneity and process
%     complexity: A new vision for watershed hydrology."
%   - Gupta, H. V., et al. (2009). "Decomposition of the mean squared error
%     and NSE for improved hydrologic model evaluation."
%
%% ========================================================================

clear variables; close all; clc;
addpath('hsB_catchment')

%% ---------------- 0) Inputs (xH, H, metrics) ----------------
% Requires xH, H, and metrics (for hillslope length) in 'input_data.mat'
%  - xH: centers of distance-to-stream bins [m]
%  - H : width per bin [m]
%  - metrics.max_hillslope_length_m
% Alternatively, you can run a pre-processing algorithm available here
% This script must produce: xH [1×n], H [1×n], Forcing (struct array), metrics (struct)
%  - xH: centers of uniform distance-to-stream bins (m)
%  - H : hillslope width per bin (m) such that area ≈ sum(H)*dxH
%  - metrics.max_hillslope_length_m: max flow distance to stream (m)
% run('hsB_PreProcessing.m')
% save('input_data','metrics','xH','H','Forcing','t_forcing_days','t0_datetime');

load('input_data.mat');   % provides: xH, H, metrics

%% ---------------- User settings: scenario + LR params ----------
recharge_rate_mmday = (11.6518/283941248)*1000*86400;    % constant recharge during filling [mm/day]
N_recharge_days     = 30;     % days with recharge on
N_release_days      = 60;     % days with recharge off (recession window)
dt_forcing_hours    = 24;      % forcing stamp interval

% Linear reservoir parameters to compare against hsB recession:
% Q_LR(t_days) = a_LR * exp(-b_LR * t_days)
a_LR  = 11.6518;        % [m^3/s]  amplitude at recession start (t=0 of recession)
b_LR  = 0.353;       % [1/day]  decay rate (days^-1). Tune manually.

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
params = struct();
params.k    = 6e-4;                               % [m/s]
params.f    = 0.05;                                % [-]
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
opts.verbose     = true;

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

% No final figs from the solver; we’ll do the comparison plot
PLOT.final_figs = false; PLOT.stream_stamps = false;
opts.stream_plot_stamps = PLOT.stream_stamps; do_plot = PLOT.final_figs;
sty = struct('fs',13,'axlw',2,'plotlw',2,'tickdir','in','ticklength',[0.02 0.02]);

%% ---------------- 5) Run hsB solver ---------------------------
[Qout, Qsurf, Qtotal, h, S, diag] = hsB_solver(x, w, ForcingNative, params, opts, do_plot, sty);

fprintf('\nRun complete: forcing steps=%d, Nx=%d, internal steps=%d\n', ...
    numel(ForcingNative.t_days), numel(x), diag.n_internal_steps);

%% ---------------- 6) Compare recession to Q=a*exp(-b*t_days) --
metrics_out = compare_recession_LR_days( ...
    ForcingNative.t_days, ...      % days
    Qout, ...                      % m^3/s
    ForcingNative.N_mps, ...       % m/s
    a_LR, b_LR, ...                % LR params (a: m^3/s, b: 1/day)
    true, ...                      % do_plot
    A_hs, ...                      % m^2  (area)
    recharge_rate_mmday, ...       % mm/day
    N_recharge_days, ...           % days with recharge
    N_release_days ...             % dry days
);

disp(metrics_out);

% -------------------------------------------------------------------------
% Helper: ghost-padding for (xH,H)
% -------------------------------------------------------------------------
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

function metrics = compare_recession_LR_days(t_days, Q_hsB, N_mps, a, b, do_plot, A_hs, rec_mmday, n_recharge_days, n_dry_days)
% Compare hsB recession to Q(t)=a*exp(-b*t) and PLOT IN mm/day.
% Adds a LaTeX metrics box with recharge rate, recharge/dry days, RMSE, NSE, PBIAS.

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
        'k0', k0 );

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

    % LaTeX metrics box with one metric per line
    ymax = max([Q_obs; Q_sim]); if ~isfinite(ymax) || ymax<=0, ymax = 1; end
    xtext = 0.02 * max(t_rec_days);
    ytext = 0.98 * ymax;

% ---------- Metrics box (LaTeX; one metric per line) ----------
lines = {
    sprintf('\\textbf{Recharge:} %.3g\\,mm/day', rec_mmday)
    sprintf('\\textbf{Days on:} %d', n_recharge_days)
    sprintf('\\textbf{Dry days:} %d', n_dry_days)
    sprintf('\\textbf{RMSE:} %.3g\\,mm/day', rmse)
    sprintf('\\textbf{NSE:} %.3f', NSE)
    sprintf('\\textbf{PBIAS:} %.2f\\,\\%%', pbias)
};

% Add textbox annotation
annotation('textbox', [0.25 0.55 0.25 0.3], ...
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
