function [Qout, Qsurf, Qtotal, h, S, diag] = hsB_solver(x, w, N_in, varargin)
% hsB_solver  Wrapper around the HSB core integrator with two call styles:
%
%   PART A — WRAPPER (dispatcher) + PLOTTING
%   ---------------------------------------
%   (A) Legacy, fixed-step path (backward compatible):
%       [Qout,Qsurf,Qtotal,h,S,diag] = hsB_solver(x, w, N, dt, Nt, params, opts, do_plot, sty)
%           x       [Nx×1] cell-center coordinates [m], strictly increasing
%           w       [Nx×1] planform width at centers [m] (>0)
%           N       [Nt×1] or [Nt×Nx] recharge [m s^-1] on the model time grid
%           dt      scalar  model time step [s]
%           Nt      scalar  number of model steps
%           params  struct  .k [m/s], .f [-], .iota [rad], .S0 [Nx×1], .D [m]
%           opts    struct  usual HSB options (theta, omega, etc.)
%           do_plot logical (default false) publication-style figures
%           sty     struct  plotting style (optional; sensible defaults)
%
%       -> Runs hsb_core ONCE on the fixed model grid (as before).
%          Output length is Nt (legacy behavior).
%
%   (B) Native forcing with internal substepping (new; preferred):
%       [Qout,Qsurf,Qtotal,h,S,diag] = hsB_solver(x, w, Forcing, params, opts, do_plot, sty)
%           Forcing is a struct with native-resolution series:
%              .t_days  [nF×1] increasing times in DAYS since a reference epoch
%              .N_mps   [nF×1] or [nF×Nx] recharge rate in m/s at those times
%           params, opts, do_plot, sty as above.
%
%       -> Integrates interval-by-interval between forcing time stamps using
%          INTERNAL SUBSTEPS chosen by a light CFL/accuracy policy:
%              opts.substep_mode   = 'adaptive'|'fixed' (default 'adaptive')
%              opts.max_sub_dt     = 6*3600  [s] upper cap per substep
%              opts.min_sub_dt     = 300     [s] lower cap per substep
%              opts.cfl_adv        = 0.6     target advective CFL
%              opts.cfl_diff       = 0.3     target diffusive CFL (accuracy)
%              opts.frac_change_eta= 0.2     max allowed fractional ΔS per substep
%              opts.fixed_nsub     = 4       if substep_mode=='fixed'
%          Output is reported at the NATIVE forcing stamps (length nF).
%
%   PART B — SUBFUNCTIONS
%   ---------------------
%     hsB_solver_native(...)     : runs the native-forcing path (internal substeps)
%     plot_hsb_results(...)      : shared plotting helper (coarse time axis)
%     build_edges_from_centers() : small grid helper (duplicated locally)
%     get_opt()                  : safe option getter with defaults
%
%   REQUIREMENTS
%   ------------
%   - The core time integrator `hsb_core.m` must be on your MATLAB path.
%     This wrapper does NOT modify `hsb_core`. It simply chooses when/how
%     to call it (once in legacy mode; many small substeps in native mode).
%
%   SIGN CONVENTION
%   ----------------
%   - Inside hsb_core: Qout is negative for outward flow (seepage).
%   - This wrapper returns:
%         Qout    [same sign as core]  (negative = outward)
%         Qsurf   [≥0] saturation-excess runoff (outward-equivalent)
%         Qtotal  [outward-positive]   = -Qout + Qsurf
%
%   MASS BALANCE DIAGNOSTICS
%   ------------------------
%   - Legacy: diag.mass_residual is the per-step residual from hsb_core.
%   - Native: the wrapper accumulates volumes per forcing interval:
%         diag.forcing_Rvol, diag.forcing_Outvol, diag.forcing_dV
%     and counts diag.n_internal_steps.
%
%   PLOTTING
%   --------
%   - In both modes, if do_plot==true, the wrapper plots:
%       Fig 1: semi-log discharge (subsurface + surface) + inverted recharge bars
%       Fig 2: h(x,t) image
%       Fig 3: mass-balance residual (legacy only; optional)
%
%   Author: UA-Hydrology (2025). Heavily commented for clarity.

% ----------------------------- DISPATCH -----------------------------
% Decide which call signature we are using based on the 3rd argument.
if isstruct(N_in) && isfield(N_in,'t_days') && isfield(N_in,'N_mps')
    % =================== NATIVE-FORCING CODE PATH ===================
    %
    % Signature: hsB_solver(x, w, Forcing, params, opts, do_plot, sty)
    %
    Forcing = N_in;
    % Pull the remaining args with defaults.
    params  = varargin{1};
    opts    = varargin{2};
    do_plot = false;  if numel(varargin)>=3 && ~isempty(varargin{3}), do_plot = varargin{3}; end
    sty     = struct(); if numel(varargin)>=4 && ~isempty(varargin{4}), sty = varargin{4}; end

    % Delegate to the native substepper
    [Qout, Qsurf, Qtotal, h, S, diag] = hsB_solver_native(x, w, Forcing, params, opts, do_plot, sty);
    return
else
    % ======================= LEGACY CODE PATH =======================
    %
    % Signature: hsB_solver(x, w, N, dt, Nt, params, opts, do_plot, sty)
    %
    N       = N_in;
    dt      = varargin{1};
    Nt      = varargin{2};
    params  = varargin{3};
    opts    = varargin{4};
    do_plot = false;  if numel(varargin)>=5 && ~isempty(varargin{5}), do_plot = varargin{5}; end
    sty     = struct(); if numel(varargin)>=6 && ~isempty(varargin{6}), sty = varargin{6}; end

    % Run the core ONCE with the fixed model time step/grid
    [Qout, Qsurf, Qtotal, h, S, diag] = hsb_core(x, w, N, dt, Nt, params, opts);

    % ---- PLOTTING (legacy) ----
    if do_plot
        % Build plotting time axis in days for Nt fixed steps
        t_days = (0:Nt-1)' * (dt/86400);

        % Area-weight recharge to mm/day for bars
        [xe, dx] = build_edges_from_centers(x); %#ok<ASGLU>
        Ahs = sum(w(:).*dx(:));
        if size(N,2) == numel(w)
            aw = w(:).*dx(:);
            N_area = (N * aw) / Ahs;  % [Nt×1] m/s
        else
            N_area = N;               % [Nt×1] m/s
        end
        N_mmday = N_area * 86400 * 1000;

        % Subsurface and surface discharges to mm/day
        Q_sub_mmday  = (-Qout / Ahs) * 86400 * 1000;   % outward-positive magnitude
        Q_surf_mmday = (Qsurf / Ahs) * 86400 * 1000;

        % Fig 1 + Fig 2 (+ Fig 3 mass residuals if available)
        plot_hsb_results(t_days, x, h, Q_sub_mmday, Q_surf_mmday, N_mmday, diag, opts, sty, 'legacy');
    end
end
end

% =========================================================================
% ====================== PART B — SUBFUNCTIONS ============================
% =========================================================================

function [Qout_F, Qsurf_F, Qtotal_F, h_F, S_F, diag] = hsB_solver_native(x, w, Forcing, params, opts, do_plot, sty)
% hsB_solver_native
% Integrate the HSB model over each forcing interval [t_k, t_{k+1}] using
% INTERNAL solver substeps chosen for stability/accuracy. Report outputs on
% the COARSE (forcing) grid only.
%
% INPUT
%   x,w        : grid and widths (as in legacy)
%   Forcing    : struct with native-resolution time series
%                   .t_days  [nF×1], increasing times in DAYS
%                   .N_mps   [nF×1] or [nF×Nx], recharge rate in m/s
%   params,opts: as usual; opts may include substep policy knobs (see header)
%   do_plot,sty: plotting toggle and style
%
% OUTPUT
%   Qout_F, Qsurf_F, Qtotal_F : [nF×1] coarse series at forcing stamps
%   h_F, S_F                  : [nF×Nx] states at forcing stamps
%   diag                      : diagnostics including:
%         .n_internal_steps
%         .forcing_Rvol, .forcing_Outvol, .forcing_dV (per forcing interval)
%         .dx, .Ahs, .Qsurf (coarse), .overflow_steps, .max_Qsurf

% --- enforce column vectors for grid + widths ---
x = x(:);
w = w(:);
Nx = numel(x);

% ------------------- Read & validate forcing -------------------
tF_days = Forcing.t_days(:);
N_F     = Forcing.N_mps;
nF      = numel(tF_days);

if size(N_F,1) ~= nF
    error('Forcing.N_mps must have nF rows to match Forcing.t_days.');
end
if isvector(N_F), N_F = N_F(:); end  % [nF×1] acceptable; [nF×Nx] also fine

% ------------------- Prepare coarse storage --------------------
Nx = numel(x);
Qout_F   = zeros(nF,1);   % subsurface at stamps (same sign conv as core)
Qsurf_F  = zeros(nF,1);   % surface (≥0) at stamps
Qtotal_F = zeros(nF,1);   % outward-positive at stamps
h_F      = zeros(nF,Nx);  % head at stamps
S_F      = zeros(nF,Nx);  % storage per unit-x at stamps

% ------------------- Rolling initial state ---------------------
params_work = params;
if ~isfield(params_work,'S0') || isempty(params_work.S0)
    params_work.S0 = zeros(Nx,1);
else
    params_work.S0 = params_work.S0(:);   % ensure column
end
S_last = params_work.S0;                   % Nx×1

% ------------------- Diagnostics & static metrics ---------------
[xe, dx]  = build_edges_from_centers(x);
Ahs       = sum(w(:).*dx(:));
n_internal_steps = 0;
forcing_Rvol   = zeros(max(nF-1,0),1);
forcing_Outvol = zeros(max(nF-1,0),1);
forcing_dV     = zeros(max(nF-1,0),1);

% ----- Coarse plotting axes and recharge (native stamps) -----
t_days = tF_days - tF_days(1);      % [nF×1] days since first stamp

% Area-weighted recharge (coarse) for bars in mm/day
if size(N_F,2) == numel(w)
    aw = w(:).*dx(:);               % [Nx×1] area weights
    N_area = (N_F * aw) / Ahs;      % [nF×1] m/s
else
    N_area = N_F;                   % [nF×1] m/s
end
N_mmday_full = N_area * 86400 * 1000;   % [nF×1] mm/day

% --- Verbosity for substeps (we’ll print our own concise line) ---
log_substeps      = get_opt(opts,'verbose',true);  % reuse existing verbose flag
opts_core         = opts;                          % copy user opts
opts_core.verbose = false;                         % silence hsb_core header prints

% ------------------- Substep policy parameters ------------------
mode       = get_opt(opts,'substep_mode','adaptive');  % 'adaptive'|'fixed'
dt_max     = get_opt(opts,'max_sub_dt',6*3600);        % [s]
dt_min     = get_opt(opts,'min_sub_dt',300);           % [s]
cfl_adv    = get_opt(opts,'cfl_adv',0.6);
cfl_diff   = get_opt(opts,'cfl_diff',0.3);
eta_max    = get_opt(opts,'frac_change_eta',0.2);
fixed_nsub = get_opt(opts,'fixed_nsub',4);

% CFL estimates (static, conservative)
k  = params.k; f = params.f; iota = params.iota;
c_adv  = (k/f)*sin(iota);                % [m/s]
dx_min = min(dx);
dt_cfl_adv  = (cfl_adv  * dx_min) / max(c_adv, eps);

% Diffusive accuracy cap: nu_eff ~ (k/f)*cos(iota)*(f*w_max*h_typ)
wmax  = max(w);
D     = inf; if isfield(params,'D') && isfinite(params.D), D = params.D; end
h_typ = min(0.3*D, 0.2);                 % mild 0.2 m if D large
S_star    = f * wmax * max(h_typ, 1e-4); % [m^2]
nu_eff    = (k/f)*cos(iota) * S_star;    % [m^2/s]
dt_cfl_diff = (cfl_diff * dx_min^2) / max(nu_eff, eps);

% ------------------- Record at the first stamp ------------------
h0 = S_last ./ max(f.*w, eps);   % Nx×1 elementwise
h_F(1,:) = h0.';                  % 1×Nx row
S_F(1,:) = S_last.';
Qout_F(1)= 0;       % not known yet; can be filled after 1st interval if desired
Qsurf_F(1)=0; Qtotal_F(1)= -Qout_F(1)+Qsurf_F(1);

% --- Stream (live) plot only at forcing stamps? ---
stream_on = do_plot && get_opt(opts,'stream_plot_stamps',false);

% --- ensure core doesn't live-plot every substep ---
opts_core = opts;
opts_core.plot_live = false;   % we will do our own quick plot at forcing stamps

% ------------------- March interval-by-interval -----------------
for kF = 1:(nF-1)
    t0 = tF_days(kF);
    t1 = tF_days(kF+1);
    Tspan = (t1 - t0) * 86400;   % [s] duration of current forcing interval
    if Tspan <= 0
        warning('Non-increasing forcing time at index k=%d. Skipping.', kF);
        continue
    end

    % Pick a starting substep size
    if strcmpi(mode,'fixed')
        dt_try = max(dt_min, min(dt_max, Tspan / max(fixed_nsub,1)));
    else
        dt_try = max(dt_min, min([dt_max, dt_cfl_adv, dt_cfl_diff, Tspan/4]));
    end

    % Recharge is piecewise-constant over [t0, t1) equal to value at t0
    Nk = N_F(kF,:);      % row: 1×1 or 1×Nx

    % Interval accumulators for coarse mass balance
    Rvol_k  = 0; Out_k = 0; dV_k = 0;

    % Integrate with internal substeps; always land on t1 exactly
    t_done = 0;
    while t_done < Tspan - 1e-9
        dt = min(dt_try, Tspan - t_done);

        % Build tiny forcing block for this substep
        if size(Nk,2)==1
            N_block = Nk(1);             % scalar → [1×1] (hsb_core will expand)
        else
            N_block = Nk;                 % [1×Nx]
        end

        % Run one substep with hsb_core
        params_work.S0 = S_last;          % rolling initial state
        [Qout_sub, Qsurf_sub, Qtotal_sub, h_sub, S_sub, diag_sub] = ...
              hsb_core(x, w, N_block, dt, 1, params_work, opts_core);

        % Fractional change monitor for adaptivity
        S_new = S_sub(end,:).';
        frac_change = max(abs(S_new - S_last) ./ max(abs(S_last), 1e-6));

        % If too large change, shrink dt and retry this chunk
        if strcmpi(mode,'adaptive') && (frac_change > eta_max) && (dt > dt_min + 1e-12)
            dt_try = max(dt_min, 0.5*dt);
            % (do not advance time; retry)
            continue
        end

        % Accept substep
        n_internal_steps = n_internal_steps + 1;
        t_done = t_done + dt;
        S_last = S_new;

        % Accumulate interval volumes (using core’s strict accounting)
        if isfield(diag_sub,'Rvol_step'), Rvol_k = Rvol_k + sum(diag_sub.Rvol_step); end
        if isfield(diag_sub,'Out_step'),  Out_k  = Out_k  + sum(diag_sub.Out_step);  end
        if isfield(diag_sub,'dV_step'),   dV_k   = dV_k   + sum(diag_sub.dV_step);   end

        % -------- Pretty substep log (t, dt, Nx, N, Storage, Out) --------
        if log_substeps
            % Area-weighted recharge for this substep (m/s)
            if size(Nk,2) == numel(w)
                aw     = w(:) .* dx(:);            % [Nx×1]
                N_area = (Nk * aw) / Ahs;          % scalar (row*col)
            else
                N_area = Nk(1);                    % uniform recharge (scalar)
            end
        
            % Convert to requested units
            N_mmday   = N_area * 86400 * 1000;                 % mm/day
            V_m3      = sum(S_last .* dx);                     % m^3
            Stor_mm   = (V_m3 / Ahs) * 1000;                   % mm (equiv depth)
            Out_mmday = (Qtotal_sub(1) / Ahs) * 86400 * 1000;  % mm/day (outward +)
        
            % Current absolute time at end of this substep (days since first stamp)
            t_cur_days = (t0 - tF_days(1)) + (t_done / 86400);
        
            fprintf('HSB solver: t=%.3f d | dt=%gs | Nx=%d | N=%.2f mm/d | Storage=%.1f mm | Out=%.2f mm/d\n', ...
                    t_cur_days, dt, numel(x), N_mmday, Stor_mm, Out_mmday);
        end

        % If the substep was very “easy”, try gently increasing dt
        if strcmpi(mode,'adaptive') && (frac_change < 0.03)
            dt_try = min([dt_max, 1.5*dt_try, Tspan - t_done]);
        end
    end

    % End of interval: record the state at t1
    h_end = S_last ./ max(f.*w, eps); % Nx×1 elementwise
    h_F(kF+1,:) = h_end.';            % 1×Nx row
    S_F(kF+1,:) = S_last.';

    % Also record coarse volumes for this interval
    forcing_Rvol(kF)   = Rvol_k;
    forcing_Outvol(kF) = Out_k;
    forcing_dV(kF)     = dV_k;

    % For plotting/reporting at stamp k+1, get an outlet flux consistent
    % with the *current* state. Easiest robust approach: a 1-second core call.
    params_tmp = params; params_tmp.S0 = S_last;
    N_here = N_F(kF+1,:); if isvector(N_here), N_here = N_here(:).'; end
    [Qout_chk, Qsurf_chk, Qtotal_chk, ~, ~, ~] = ...
         hsb_core(x, w, N_here, 1.0, 1, params_tmp, opts_core);  % muted core header
    Qout_F(kF+1)   = Qout_chk(1);           % (negative = outward)
    Qsurf_F(kF+1)  = Qsurf_chk(1);
    Qtotal_F(kF+1) = -Qout_chk(1) + Qsurf_chk(1);

    % ----------------- OPTIONAL: quick plot only at forcing stamps -----------------
    if do_plot && get_opt(opts,'stream_plot_stamps',false)
        % latest state at the new forcing stamp
        h_now = h_F(kF+1,:).';     % [Nx×1]
        % aquifer cap (for y-lim), if available
        Dloc = inf; 
        if isfield(params,'D') && ~isempty(params.D) && isfinite(params.D)
            Dloc = params.D;
        end
        % time in days since first stamp (for the title)
        t_now_days = (tF_days(kF+1) - tF_days(1));
        % draw the light 2-panel plot (h and w)
        quick_plot_head_width(x, w, h_now, t_now_days, Dloc);
        drawnow limitrate;
    end

end

% ------------- Assemble diagnostics -------------
diag = struct();
diag.dx               = dx;
diag.Ahs              = Ahs;
diag.n_internal_steps = n_internal_steps;
diag.forcing_Rvol     = forcing_Rvol;
diag.forcing_Outvol   = forcing_Outvol;
diag.forcing_dV       = forcing_dV;
diag.Qsurf            = Qsurf_F;                 % coarse series
diag.overflow_steps   = sum(Qsurf_F > 0);
diag.max_Qsurf        = max(Qsurf_F, [], 'omitnan');

% ------------- Optional plotting -------------
if do_plot
    % Coarse time axis starting at zero (days)
    t_days = tF_days - tF_days(1);

    % Area-weight recharge to mm/day for bars (coarse)
    if size(N_F,2) == numel(w)
        aw = w(:).*dx(:);
        N_area = (N_F * aw) / Ahs;   % [nF×1] m/s
    else
        N_area = N_F;                % [nF×1] m/s
    end
    N_mmday = N_area * 86400 * 1000;

    % Subsurface/surface to mm/day (coarse stamps)
    Q_sub_mmday  = (-Qout_F / Ahs) * 86400 * 1000;
    Q_surf_mmday = (Qsurf_F / Ahs) * 86400 * 1000;

    plot_hsb_results(t_days, x, h_F, Q_sub_mmday, Q_surf_mmday, N_mmday, diag, opts, sty, 'native');
end
end

% -------------------------------------------------------------------------
function plot_hsb_results(t_days, x, h_xt, Q_sub_mmday, Q_surf_mmday, N_mmday, diag, opts, sty, mode_tag)
% plot_hsb_results  Publication-style figures used by both paths.
%
% INPUT
%   t_days       : [Nt|nF × 1] time in days (monotone increasing)
%   x            : [Nx × 1] grid centers [m]
%   h_xt         : [Nt|nF × Nx] head above bedrock [m]
%   Q_sub_mmday  : [Nt|nF × 1] subsurface (outward-positive magnitude) [mm/day]
%   Q_surf_mmday : [Nt|nF × 1] surface runoff (outward-positive) [mm/day]
%   N_mmday      : [Nt|nF × 1] recharge (area-weighted) [mm/day]
%   diag         : diagnostics (contains diag.dx, diag.Ahs; mass_residual if legacy)
%   opts         : options (mass_check for legacy residual plot)
%   sty          : style (fields are optional; sensible defaults applied)
%   mode_tag     : 'legacy' or 'native' (affects residual figure)

% ---- style defaults ----
if ~exist('sty','var') || isempty(sty), sty = struct(); end
if ~isfield(sty,'axlw'),        sty.axlw = 2.5; end
if ~isfield(sty,'plotlw'),      sty.plotlw = 2.5; end
if ~isfield(sty,'tickdir'),     sty.tickdir = 'in'; end
if ~isfield(sty,'ticklength'),  sty.ticklength = [0.02 0.02]; end
if ~isfield(sty,'fs'),          sty.fs = 13; end
if ~isfield(sty,'c') || ~isfield(sty.c,'teal'), sty.c.teal = [40 141 141]/255; end
if ~isfield(sty,'c') || ~isfield(sty.c,'red'),  sty.c.red  = [159 0 0]/255; end
if ~isfield(sty,'c') || ~isfield(sty.c,'gray'), sty.c.gray = [72 72 72]/255; end
if ~isfield(sty,'c') || ~isfield(sty.c,'blue'), sty.c.blue = [0 114 178]/255; end
if ~isfield(sty,'cmap_name'),   sty.cmap_name = 'turbo'; end
if ~isfield(sty,'cmap_levels'), sty.cmap_levels = 21; end
if ~isfield(sty,'cmap_clim'),   sty.cmap_clim = []; end
if ~isfield(sty,'n_profiles'),  sty.n_profiles = 5; end

Ntplot = numel(t_days);
Nx     = numel(x);

% ---------------- Fig 1: semi-log discharges + inverted recharge bars ---------------
figure('Color','w','Position',[120 90 1100 450]);
tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Left tile: discharges on log y-axis
ax1 = nexttile(tl);
yyaxis(ax1,'left');

yfloor = 1e-6; % mm/day floor so log scale behaves with zeros
Q1 = max(Q_sub_mmday,  yfloor);
Q2 = max(Q_surf_mmday, yfloor);
p1 = plot(ax1, t_days, Q1, 'LineWidth', sty.plotlw, 'Color', sty.c.teal); hold(ax1,'on');
p2 = plot(ax1, t_days, Q2, 'LineWidth', sty.plotlw, 'Color', sty.c.red, 'LineStyle','--'); hold(ax1,'off');
set(ax1,'YScale','log');
grid(ax1,'on'); box(ax1,'on');
set(ax1,'TickDir',sty.tickdir,'TickLength',sty.ticklength, ...
    'LineWidth',sty.axlw,'FontSize',sty.fs, ...
    'GridColor',sty.c.gray,'GridAlpha',0.25, ...
    'TickLabelInterpreter','latex');
xlabel(ax1,'\textbf{Time (days)}','Interpreter','latex');
ylabel(ax1,'\textbf{Discharge (mm d$^{-1}$)}','Interpreter','latex');

% --- Right axis: recharge bars (put bars first so lines sit on top) ---
yyaxis(ax1,'right');
p3 = bar(ax1, t_days, N_mmday, ...
    'FaceColor', sty.c.blue, ...
    'EdgeColor', 'none', ...
    'BarWidth', 1, ...
    'FaceAlpha', 0.35, ...              % make bars translucent
    'DisplayName', 'Recharge N');       % legend label for bars
set(ax1,'YDir','reverse');              % downward bars = positive recharge
ylabel(ax1,'\textbf{$N$ (mm d$^{-1}$)}','Interpreter','latex');
ax1.YAxis(2).Color = sty.c.blue;

% --- Left axis back on top: ensure line y-axis is black and lines above bars ---
yyaxis(ax1,'left');
ax1.YAxis(1).Color = 'k';
uistack(p1,'top');                       % keep subsurface line above bars
uistack(p2,'top');                       % keep surface line above bars

% Legend with explicit handles (correct labels)
legend(ax1, [p1 p2 p3], ...
    {'$\mathrm{Subsurface}\ Q_{\mathrm{out}}$', ...
     '$\mathrm{Surface}\ Q_{\mathrm{surf}}$', ...
     '$\mathrm{Recharge}\ N$'}, ...
    'Interpreter','latex','Location','best');

% Right tile: water-table profiles (sampled across time)
ax2 = nexttile(tl);
idx = unique(round(linspace(1, Ntplot, max(2,round(sty.n_profiles)))));
if exist(sty.cmap_name,'file')==2 || exist(sty.cmap_name,'builtin')==5
    cmapP = feval(sty.cmap_name, numel(idx));
else
    cmapP = lines(numel(idx));
end
hold(ax2,'on');
for k = 1:numel(idx)
    plot(ax2, x, h_xt(idx(k),:), 'LineWidth', sty.plotlw, 'Color', cmapP(k,:));
end
hold(ax2,'off');
grid(ax2,'on'); box(ax2,'on');
set(ax2,'TickDir',sty.tickdir,'TickLength',sty.ticklength, ...
    'LineWidth',sty.axlw,'FontSize',sty.fs, ...
    'GridColor',sty.c.gray,'GridAlpha',0.25, ...
    'TickLabelInterpreter','latex');
xlabel(ax2,'\textbf{x (m)}','Interpreter','latex');
ylabel(ax2,'\textbf{h (m above bedrock)}','Interpreter','latex');
title(ax2,'\textbf{Water table profiles}','Interpreter','latex');
ylim(ax2,[0, max(h_xt(:))*1.05]);
leg = arrayfun(@(tt) sprintf('$t=%.2f$ d', t_days(idx(tt))), 1:numel(idx), 'UniformOutput', false);
legend(ax2, leg, 'Interpreter','latex','Location','best');

% ---------------- Fig 2: h(x,t) image ---------------------------------------------
figure('Color','w','Position',[120 560 880 420]);
ax3 = axes();
imagesc(ax3, x, t_days, h_xt); set(ax3,'YDir','normal');
grid(ax3,'on'); box(ax3,'on');
set(ax3,'TickDir',sty.tickdir,'TickLength',sty.ticklength, ...
    'LineWidth',sty.axlw,'FontSize',sty.fs, ...
    'GridColor',sty.c.gray,'GridAlpha',0.15, ...
    'TickLabelInterpreter','latex');
xlabel(ax3,'\textbf{x (m)}','Interpreter','latex');
ylabel(ax3,'\textbf{time (days)}','Interpreter','latex');
title(ax3,'\textbf{$h(x,t)$}','Interpreter','latex');
if ~(exist(sty.cmap_name,'file')==2 || exist(sty.cmap_name,'builtin')==5)
    warning('Colormap "%s" not found. Using "parula".', sty.cmap_name);
    sty.cmap_name = 'parula';
end
cmap = feval(sty.cmap_name, max(2,round(sty.cmap_levels)));
colormap(ax3, cmap);
if ~isempty(sty.cmap_clim) && numel(sty.cmap_clim)==2
    caxis(ax3, sort(sty.cmap_clim));
end
cb = colorbar(ax3);
set(cb,'LineWidth',sty.axlw,'TickDirection',sty.tickdir,'TickLabelInterpreter','latex');

% ---------------- Fig 3: mass-balance residual (legacy only) -----------------------
is_legacy = strcmpi(mode_tag,'legacy');
if is_legacy && isfield(opts,'mass_check') && opts.mass_check && isfield(diag,'mass_residual') ...
        && ~isempty(diag.mass_residual)
    figure('Color','w','Position',[1020 560 520 420]);
    ax4 = axes();
    plot(ax4, t_days, diag.mass_residual, 'LineWidth', sty.plotlw, 'Color', sty.c.teal);
    hold(ax4,'on'); yline(ax4,0,'-','LineWidth',sty.axlw,'Color',sty.c.gray); hold(ax4,'off');
    grid(ax4,'on'); box(ax4,'on');
    set(ax4,'TickDir',sty.tickdir,'TickLength',sty.ticklength, ...
        'LineWidth',sty.axlw,'FontSize',sty.fs, ...
        'GridColor',sty.c.gray,'GridAlpha',0.25, ...
        'TickLabelInterpreter','latex');
    xlabel(ax4,'\textbf{time (days)}','Interpreter','latex');
    ylabel(ax4,'\textbf{MB residual (m$^3$)}','Interpreter','latex');
    title(ax4,'\textbf{Mass-balance residual per step}','Interpreter','latex');
end
end

% -------------------------------------------------------------------------
function [xe, dx] = build_edges_from_centers(x)
% build_edges_from_centers  Deduce face locations (xe) and cell widths (dx)
% from strictly increasing centers x (1D). The left boundary is set at 0,
% the right boundary by linear extrapolation of the last gap.
Nx = numel(x);
xe = zeros(Nx+1,1);
xe(1)     = 0.0;                                % outlet face
if Nx >= 2
    xe(2:Nx)  = 0.5 * (x(1:Nx-1) + x(2:Nx));    % interior midpoints
    xe(Nx+1)  = x(end) + 0.5*(x(end) - x(end-1)); % extrapolated right face
else
    xe(2)     = x(1);                            % degenerate 1-cell case
end
dx = diff(xe);
end

% -------------------------------------------------------------------------
function v = get_opt(s, f, dflt)
% get_opt  Safe options getter with default
v = dflt;
if nargin<3, dflt = []; end
if ~isempty(s) && isstruct(s) && isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
end
end

% -------------------------------------------------------------------------
% ============================ CORE INTEGRATOR =============================
% -------------------------------------------------------------------------
function [Qout, Qsurf, Qtotal, h_store, S_store, diag] = hsb_core(x, w, N, dt, Nt, params, opts)
% hsb_core  Semi-implicit FV + Picard with stabilized fluxes and overflow sink.
%
% INPUTS (same as legacy path):
%   x      [Nx×1] strictly increasing centers
%   w      [Nx×1] positive widths
%   N      [Nt×1] or [Nt×Nx] recharge (m/s) on the model time grid
%   dt     scalar  time step (s)
%   Nt     scalar  number of steps
%   params struct  .k [m/s], .f [-], .iota [rad], .S0 [Nx×1], .D [m] (cap)
%   opts   struct  .theta, .omega, .picard_max, .picard_tol, .safeguard,
%                  .verbose, .mass_check, .plot_live, .no_recharge_outlet
%
% OUTPUTS:
%   Qout    [Nt×1] subsurface discharge at outlet face (negative = outward)
%   Qsurf   [Nt×1] saturation-excess (outward-equiv) (≥0)
%   Qtotal  [Nt×1] outward-positive = -Qout + Qsurf
%   h_store [Nt×Nx] head (m)
%   S_store [Nt×Nx] storage per unit-x (= f w h) (m^2)
%   diag    struct  dx, xe, Ahs, mass_residual, per-step volumes, etc.

% ---- sanity & defaults ----
x = x(:); w = w(:);
Nx = numel(x);
assert(all(diff(x)>0), 'x must be strictly increasing');
assert(numel(w)==Nx && all(w>0), 'w must be length Nx and positive');

if isvector(N), N = N(:); end
if size(N,1)~=Nt, error('N must have Nt rows'); end
if size(N,2)==1
    N = repmat(N,1,Nx);          % broadcast uniform recharge
elseif size(N,2)~=Nx
    error('N must be Nt×1 or Nt×Nx');
end

% Optionally zero recharge at outlet micro-cell
no_recharge_outlet = get_opt(opts,'no_recharge_outlet',false);
if no_recharge_outlet, N(:,1) = 0; end

% ---- Parameters ----
k    = params.k;
f    = params.f;
iota = params.iota;
if isfield(params,'S0') && ~isempty(params.S0)
    S = params.S0(:);
else
    S = zeros(Nx,1);
end

% Finite aquifer thickness cap (0 <= h <= D)
if isfield(params,'D') && isfinite(params.D) && params.D > 0
    D    = params.D;
    Smax = f .* w * D;           % [Nx×1]
else
    D    = inf;
    Smax = inf(Nx,1);
end

% ---- Options ----
theta      = get_opt(opts,'theta',1.0);      % fully implicit
omega0     = get_opt(opts,'omega',0.6);      % base Picard damping
picard_max = get_opt(opts,'picard_max',20);
picard_tol = get_opt(opts,'picard_tol',1e-8);
verbose    = get_opt(opts,'verbose',true);
safeguard  = get_opt(opts,'safeguard',true);
mass_check = get_opt(opts,'mass_check',true);
plot_live  = get_opt(opts,'plot_live',false);

% ---- Grid metrics ----
[xe, dx] = build_edges_from_centers(x);

% ---- Allocate outputs ----
Qout    = zeros(Nt,1);
Qsurf   = zeros(Nt,1);          % saturation-excess overland flow (m^3 s^-1)
h_store = zeros(Nt,Nx);
S_store = zeros(Nt,Nx);
mass_residual = nan(Nt,1);

% Per-step volumes for strict balance
Rvol_step = zeros(Nt,1);        % recharge volume per step   [m^3]
Out_step  = zeros(Nt,1);        % outward-positive outflow   [m^3]
dV_step   = zeros(Nt,1);        % storage change per step    [m^3]

% ---- Helpers and constants ----
cosI = cos(iota);
sinI = sin(iota);
Ahs  = sum(w .* dx);             % plan area (for diagnostics)

if verbose
    fprintf('HSB solver: Nx=%d, Nt=%d, dt=%.0fs, theta=%.2f\n', Nx, Nt, dt, theta);
end

% ====================== TIME LOOP ======================
for n = 1:Nt
    S_old = S;
    Nn    = N(n,:).';            % recharge at step n (column)

    % Picard initialization + per-step damping reset
    omega     = omega0;
    incr_prev = Inf;
    S_iter    = S_old;

    % -------------------- Picard iterations --------------------
    for m = 1:picard_max
        % state used to evaluate fluxes at n+theta
        S_theta = (1-theta)*S_old + theta*S_iter;

        % fluxes at faces using stabilized scheme
        h_cent = S_theta ./ max(f.*w, eps);
        Qface  = compute_face_fluxes(x, xe, w, h_cent, S_theta, k, f, cosI, sinI);

        % conservative update for storage
        divQ   = -(Qface(2:end) - Qface(1:end-1)) ./ dx;
        RHS    = divQ + Nn .* w;                    % m^3 s^-1 per unit-x
        S_new  = S_old + dt * RHS;

        % safeguards: positivity + optional head cap h<=D
        if safeguard, S_new = max(S_new, 0); end
        if isfinite(D), S_new = min(S_new, Smax); end

        % relaxed Picard update
        S_relaxed = (1-omega)*S_iter + omega*S_new;

        % simple adaptive damping
        incr_now = max(abs(S_relaxed - S_iter));  % L-infinity increment
        if incr_now > 1.05*incr_prev
            omega = max(0.2, 0.5*omega);          % damp more (not below 0.2)
        end

        % convergence test on relaxed iterate
        if incr_now < picard_tol
            S_iter = S_relaxed;
            break
        end

        % continue iterating
        S_iter    = S_relaxed;
        incr_prev = incr_now;
    end
    % ------------------ end Picard iterations ------------------

    % accept step
    S = S_iter;

    % store outputs
    h_now = S ./ max(f.*w, eps);
    h_store(n,:) = h_now.';
    S_store(n,:) = S.';

    % θ-state consistent with the time discretization
    S_theta = (1-theta)*S_old + theta*S;
    h_theta = S_theta ./ max(f.*w, eps);

    % Face fluxes and divergence at θ-state
    Qface_theta = compute_face_fluxes(x, xe, w, h_theta, S_theta, k, f, cosI, sinI);
    Qout(n)     = Qface_theta(1);   % NOTE: negative = outward
    divQ_theta  = -(Qface_theta(2:end) - Qface_theta(1:end-1)) ./ dx;

    % Provisional unconstrained update at θ-state (for overflow accounting)
    S_unclipped = S_old + dt * (divQ_theta + Nn .* w);  % could violate caps
    if safeguard, S_unclipped = max(S_unclipped, 0); end
    S_clip      = min(S_unclipped, Smax);
    excess      = max(S_unclipped - S_clip, 0);         % [m^2] per cell
    Qsurf(n)    = sum(excess .* dx) / dt;               % [m^3 s^-1]

    % --- Step volumes (strict accounting with actual dt)
    Rvol_step(n) = sum((Nn .* w) .* dx) * dt;           % recharge m^3 in step
    Out_step(n)  = (-Qout(n) + Qsurf(n)) * dt;          % outward-positive m^3
    dV_step(n)   = sum((S - S_old) .* dx);              % storage change m^3

    % --- Mass-balance residual at this step (should be ~0)
    if mass_check
        R_excess      = excess / dt;                    % [m^2 s^-1] sink per unit-x
        res_cell      = (S - S_old)/dt - (divQ_theta + Nn .* w - R_excess);
        mass_residual(n) = sum(res_cell .* dx);         % [m^3]
    end

    % (optional) ultra-light live plot — generally disable in substepping
    if plot_live && mod(n, max(1,round(Nt/100)))==0
        % quick_plot_step(x, w, h_now, n, D, dt); % keep commented by default
        drawnow;
    end
end
% ==================== END TIME LOOP ====================

% outward-positive total discharge
Qtotal = -Qout + Qsurf;

% Diagnostics out
diag = struct( ...
    'mass_residual',   mass_residual, ...
    'dx',              dx, ...
    'xe',              xe, ...
    'Ahs',             Ahs, ...
    'Qsurf',           Qsurf, ...
    'Rvol_step',       Rvol_step, ...    % recharge m^3 per step
    'Out_step',        Out_step,  ...    % outward m^3 per step
    'dV_step',         dV_step,   ...    % storage change per step
    'overflow_mask',   Qsurf > 0, ...
    'overflow_steps',  nnz(Qsurf > 0), ...
    'max_Qsurf',       max(Qsurf, [], 'omitnan') ...
);
end

% -------------------------------------------------------------------------
function Qface = compute_face_fluxes(xc, xe, w, h_c, S_c, k, f, cosI, sinI)
% Monotone face flux = centered diffusion (harmonic S) + upwinded drift
% Q = Qdiff + Qadv, with small Lax–Friedrichs viscosity on the drift.
Nx = numel(xc);
Qface = zeros(Nx+1,1);

Kd_over_f = (k/f);

% Stronger minimum "film" to keep conductance from vanishing when w is tiny
h_floor   = 2e-3;          % was 1e-4 → now 2 mm equivalent head
Sfilm_abs = 1e-5;          % absolute floor per face [m^2] (tune 1e-6–1e-4)

S_floor_L   = max(f * w(1) * h_floor, Sfilm_abs);
S_floor_int = max(f * 0.5*(w(1:end-1)+w(2:end)) * h_floor, Sfilm_abs);

% ---- outlet face (x=0): Dirichlet h(0)=0
dxL   = max(xc(1) - xe(1), eps);
gradL = (h_c(1) - 0) / dxL;
Sf    = harmonic_mean(max(S_c(1),S_floor_L), max(S_c(1),S_floor_L));
Qdiff = - Kd_over_f * cosI * Sf * gradL;

Sup   = S_c(1);                                   % upwind for drift (downslope)
Qadv  = - Kd_over_f * sinI * Sup;

lf_fac = 2.0;                                     % ↑ a bit more dissipation
a      = lf_fac * abs(Kd_over_f * sinI);          % LF viscosity speed
Qface(1) = Qdiff + Qadv + 0.5*a*(0 - Sup);        % boundary has no right state

% ---- interior faces i=2..Nx
for i = 2:Nx
    dxloc = max(xc(i) - xc(i-1), eps);
    grad  = (h_c(i) - h_c(i-1)) / dxloc;

    % centered diffusion with harmonic S and small floor
    Sf    = harmonic_mean(max(S_c(i),S_floor_int(i-1)), ...
                          max(S_c(i-1),S_floor_int(i-1)));
    Qdiff = - Kd_over_f * cosI * Sf * grad;

    % downslope drift: left cell is upwind
    Sup   = S_c(i-1);
    Sdn   = S_c(i);
    Qadv  = - Kd_over_f * sinI * Sup;

    % local Lax–Friedrichs stabilization on the drift
    lf_fac = 2.0;                                     % more Rusanov dissipation
    a      = lf_fac * abs(Kd_over_f * sinI);
    Qnum  = 0.5*a*(Sup - Sdn);

    Qface(i) = Qdiff + Qadv + Qnum;
end

% ---- divide face: no-flow
Qface(Nx+1) = 0.0;
end

% -------------------------------------------------------------------------
function Sh = harmonic_mean(a,b)
a = max(a,0); b = max(b,0);
Sh = 2*a.*b ./ max(a+b, eps);
end

function quick_plot_head_width(x, w, h, t_now_days, D)
% quick_plot_head_width  Minimal live view: h(x) and w(x) side-by-side.
% - x [Nx×1], w [Nx×1], h [Nx×1]
% - t_now_days: time (days) to show in the title
% - D: aquifer thickness cap (for y-limit); use inf if not capped

% --- minimalist style (same vibe as earlier quick plot) ---
axlw   = 2.0;   plw = 2.0;   fs = 12;
tickdir = 'in'; ticklen = [0.02 0.02];
c_teal = [40 141 141]/255;   c_gray = [72 72 72]/255;

if ~ishandle(9999), figure(9999); clf; end
fig = figure(9999);
set(fig,'Color','w','Position',[360 360 880 360]);

tl = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

% ---- Panel 1: head h(x) ----
ax1 = nexttile(tl);
plot(ax1, x, h, '-', 'LineWidth', plw, 'Color', c_teal);
grid(ax1,'on'); box(ax1,'on');
set(ax1,'TickDir',tickdir,'TickLength',ticklen, ...
    'LineWidth',axlw,'FontSize',fs, ...
    'GridColor',c_gray,'GridAlpha',0.25, ...
    'TickLabelInterpreter','latex');
xlabel(ax1,'\textbf{x (m)}','Interpreter','latex');
ylabel(ax1,'\textbf{h (m)}','Interpreter','latex');
if isfinite(D)
    ylim(ax1,[0, D]);
else
    ylim(ax1,[0, max(h)*1.05 + eps]);
end
title(ax1,'\textbf{Head}','Interpreter','latex');

% ---- Panel 2: width w(x) ----
ax2 = nexttile(tl);
plot(ax2, x, w, '-', 'LineWidth', plw, 'Color', c_gray);
grid(ax2,'on'); box(ax2,'on');
set(ax2,'TickDir',tickdir,'TickLength',ticklen, ...
    'LineWidth',axlw,'FontSize',fs, ...
    'GridColor',c_gray,'GridAlpha',0.25, ...
    'TickLabelInterpreter','latex');
xlabel(ax2,'\textbf{x (m)}','Interpreter','latex');
ylabel(ax2,'\textbf{w (m)}','Interpreter','latex');
title(ax2,'\textbf{Width}','Interpreter','latex');

% ---- Title across both ----
sgtitle(tl, sprintf('\\textbf{t = %.2f d (forcing stamp)}', t_now_days), 'Interpreter','latex');
end