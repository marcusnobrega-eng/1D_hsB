% -------------------------------------------------------------------------
% Driver script for the Hillslope-Storage Boussinesq (HSB) model
% -------------------------------------------------------------------------
% Author: Marcus Nobrega Gomes Junior, PhD
% University of Arizona, Department of Hydrology and Atmospheric Sciences
% Sep/2025
%
% Purpose:
%   Build inputs (grid, widths, recharge, physical parameters, solver opts),
%   then run the HSB solver and optionally plot results.
%
% What this script does:
%   1) Loads preprocessed catchment inputs (xH,H,Forcing,metrics).
%   2) Builds a 1-D hillslope grid (clustered near the outlet).
%   3) Pads the histogram-based width function with **ghost bins** on both
%      ends (purely geometric; no mass is added). This gives consistent
%      ghost nodes for centered-difference fluxes **without touching the solver**.
%   4) Maps widths from (ghost-padded) bins to the custom grid cells.
%   5) (Part 2) Sets time/recharge, parameters, options, runs the solver.
% -------------------------------------------------------------------------

clear variables; close all; clc;

addpath('hsB_catchment')

%% ---------------- 0) Preprocess Inputs ----------------
% This script must produce: xH [1×n], H [1×n], Forcing (struct array), metrics (struct)
%  - xH: centers of uniform distance-to-stream bins (m)
%  - H : hillslope width per bin (m) such that area ≈ sum(H)*dxH
%  - metrics.max_hillslope_length_m: max flow distance to stream (m)
% run('hsB_PreProcessing.m')
% save('input_data','metrics','xH','H','Forcing','t_forcing_days','t0_datetime');

%% ---------------- 1) Grid + Ghost-padded width function ----------------
clear all; close all; clc;
load('input_data.mat')

% --- 1a) Infer the histogram bin size dxH from xH (must be ~uniform)
xH  = xH(:).';                   % row
H   = H(:).';                    % row

% H(1:end) = H(1);               % Assuming fixed hillslope width
if numel(xH) < 2
    error('xH must have at least two bins.');
end
dxH_vec = diff(xH);
dxH     = median(dxH_vec);       % robust to roundoff
if any(abs(dxH_vec - dxH) > 1e-6*dxH)
    warning('xH is not perfectly uniform. Proceeding with dxH = median(diff(xH)).');
end

% Physical domain length represented by bins: [0, Lphys]
% If xH are bin centers, the physical length is the last center + half a bin.
Lphys_bins = xH(end) + dxH/2;     % m (from histogram)
Lphys_dem  = metrics.max_hillslope_length_m;   % m (from DEM)
% Choose what you want to model as domain length:
L = Lphys_dem;                     % recommended: use DEM-based length

% --- 1b) Build a **clustered** grid near the outlet (x=0) that ends at x=L
% Define the total number of *model* cells you want.
% A common choice is to keep #cells ≈ #bins, but you can change it.
Nx_target = numel(xH);

% Cluster M small cells near outlet using a geometric sequence, then uniform
M   = 3;                 % number of tiny near-outlet cells
r   = 1.15;              % geometric growth factor (>1)
d0  = 5;                 % smallest cell width [m] at the outlet

dx_geo = d0 * r.^(0:M-1);          % geometric small cells
len_geo = sum(dx_geo);             % total length of clustered block
if len_geo >= L
    error('Clustered block (%.2f m) exceeds or equals the domain length L=%.2f m.', len_geo, L);
end

Nx = Nx_target;                    % total model cells
if Nx <= M
    error('Nx (=#cells) must be > M (clustered cells). Set Nx_target > M.');
end

dx_uni = (L - len_geo) / (Nx - M); % uniform width for the rest
dx_all = [dx_geo, dx_uni * ones(1, Nx-M)];  % [1×Nx]

% Edges/centers of the **model grid** (no ghosts are added to x)
xe = [0; cumsum(dx_all(:))];       % [Nx+1 × 1], runs 0 → L
x  = 0.5 * (xe(1:end-1) + xe(2:end));  % [Nx × 1], cell centers

% --- 1c) Ghost-pad (xH,H) to form (xHg,Hg) with one ghost on each side
% Purpose: provide consistent "outside" widths for centered-difference
% fluxes at boundaries *without* special-casing the solver.
[xHg, Hg, dxH_check, Lphys_bins_check] = add_ghosts_to_width_function(xH, H);

% Sanity: ensure we matched dxH
if abs(dxH_check - dxH) > 1e-9*dxH, dxH = dxH_check; end
Lphys_bins = Lphys_bins_check;   %#ok<NASGU>

% --- 1d) Map widths from **ghost-padded bins** to model cells
% Bins (with ghosts) conceptually have edges at: ..., -dxH, 0, dxH, 2dxH, ...
% We want the **first physical model cell** (x≈0+) to use **H(1)**, i.e.,
% the first *physical* bin, not the left-ghost width. To do that cleanly:
%   idx = floor(x/dxH) + 2
% gives:
%   x ∈ [0, dxH)     → idx = 2 → Hg(2) = H(1)  (1st physical bin)
%   x ∈ [dxH, 2dxH)  → idx = 3 → Hg(3) = H(2)  (2nd physical bin), etc.
% idx = floor(x / dxH) + 2;                          % [Nx×1], 1-based
% idx = max(1, min(numel(Hg), idx));                 % clamp safely
% 
% w = Hg(idx);                                       % [Nx×1] width per model cell

% --- 1d) Map widths from **ghost-padded bins** to model cells (SMOOTH)
% Use a shape-preserving cubic so w(x) is C^1 and monotone where data are.
w = interp1(xHg(:), Hg(:), x(:), 'pchip', 'extrap');   % [Nx×1]

% --- 1e) Optional: exact **area closure** between histogram and model grid
% Area implied by histogram (physical part only): A_hist = sum(H)*dxH
% Area on the model grid:                         A_grid = sum(w .* dx)
A_hist = sum(H) * dxH;
A_grid = sum(w(:) .* diff(xe));
if A_grid <= 0
    error('Computed A_grid <= 0. Check grid construction.');
end
scale_area = A_hist / A_grid;
w = w * scale_area;                   % enforce ∑ w_j dx_j = ∑ H_k dxH

% Quick diagnostic print
fprintf('Area closure: A_hist=%.3e m^2, A_grid=%.3e m^2, scale=%.6f\n', A_hist, A_grid, scale_area);

%% ---------------- 2) Forcing at native resolution (no interpolation) ----------------
% Use the forcing exactly as provided: times in days and rates in mm/day.
tF_days = t_forcing_days(:);                 % [nF×1] days since t0_datetime
prec_mm_day = [Forcing.Prec].';             % [nF×1] precipitation rate (mm/day)

% If your input were totals per interval, convert to rate here (optional):
forcing_is_rate = true;
if ~forcing_is_rate
    dF = diff(tF_days);
    if isempty(dF), error('Need ≥2 forcing timestamps to infer interval.'); end
    d_last = median(dF);
    step_days = [dF; d_last];
    prec_mm_day = prec_mm_day ./ step_days;
end

% Convert to SI (m/s), still at native resolution; no time interpolation:
N_forcing_mps = (prec_mm_day / 1000) / 86400;    % [nF×1] m/s

% Package forcing to pass to the solver wrapper:
ForcingNative = struct();
ForcingNative.t_days = tF_days;                  % increasing, length nF
ForcingNative.N_mps  = N_forcing_mps;            % [nF×1] or later [nF×Nx]



%% ---------------- 3) Physical parameters ----------------
Nx = numel(x);                          % number of model cells (NOT #bins)
params.k     = 5e-3;                    % saturated hydraulic conductivity [m/s]
params.f     = 0.1;                     % drainable porosity [-]
params.iota  = metrics.hillslope_angle_deg * pi/180; % bedrock slope angle [rad]
params.S0    = zeros(Nx,1);             % initial storage per unit-x [m^2]
params.D     = 1.0;                     % aquifer thickness [m] (finite)

%% ---------------- 4) Solver options ----------------
opts.theta      = 1.0;        % 1.0 = fully implicit (stable & robust)
opts.omega      = 0.5;        % Picard relaxation (0.4–0.7 typical)
opts.picard_max = 20;         % max Picard iterations per step
opts.picard_tol = 1e-8;       % L∞ tolerance on Picard increment (m^2)
opts.safeguard  = true;       % enforce S>=0
opts.mass_check = true;       % compute mass-balance residuals
opts.plot_live  = false;       % light live h(x) plot while integrating
opts.verbose    = true;       % print progress each step block

% --- Internal substepping policy (used by hsB_solver wrapper) ---
opts.substep_mode   = 'adaptive';   % 'adaptive' or 'fixed'
opts.max_sub_dt     = 24*3600;      % [s] upper cap per substep (e.g., 6 h)
opts.min_sub_dt     = 300;          % [s] floor (e.g., 5 min)
opts.cfl_adv        = 0.1;          % target CFL for drift
opts.cfl_diff       = 0.1;          % target CFL for diffusion (accuracy)
opts.picard_easy    = 4;            % if iter <= this, may increase dt
opts.picard_hard    = 10;           % if iter >  this, cut dt and retry
opts.frac_change_eta= 0.4;          % max allowed fractional change per step
opts.fixed_nsub     = 4;            % used if substep_mode == 'fixed'


% NOTE: With the ghost-node approach, your first cell is fully physical.
% If you previously used a "micro seepage cell", you might have zeroed
% recharge at the outlet. That is no longer required; treat all cells equally:
opts.no_recharge_outlet = false;

%% ---------------- 5) Plot controls (simple & explicit) ----------------
% What to draw?
PLOT.final_figs    = true;   % publication-style figures AFTER the run
PLOT.stream_stamps = false;  % live QUICK head/width figure, updated ONLY at forcing stamps

% Hand these two flags to the solver:
opts.stream_plot_stamps = PLOT.stream_stamps;   % quick head/width during run
do_plot = PLOT.final_figs;                      % final publication plots at end

% Shared styling (used by both final figs and quick plot)
sty = struct();
sty.fs     = 13;     % font size
sty.axlw   = 2.0;    % axes line width
sty.plotlw = 2.0;    % curve line width
sty.tickdir    = 'in';
sty.ticklength = [0.02 0.02];

% Colors (final figures)
sty.c = struct( ...
    'teal', [40 141 141]/255, ...  % subsurface discharge
    'red',  [159 0 0]/255,   ...   % surface discharge
    'gray', [72 72 72]/255,  ...   % grids/axes
    'blue', [0 114 178]/255  ...   % recharge bars
);

% Colormap for h(x,t) in the final figure
sty.cmap_name   = 'turbo';   % e.g., 'turbo' | 'parula'
sty.cmap_levels = 21;
sty.cmap_clim   = [];        % [] auto; or [cmin cmax]
sty.n_profiles  = 5;         % number of water-table profiles overlaid

% Quick-plot layout (used ONLY when PLOT.stream_stamps == true)
% If your solver’s quick plot reads these fields, they’ll be applied; otherwise they’re harmless.
sty.quick_xlim  = [0, L];           % x-range for quick head/width (m); [] to auto
sty.quick_hmax  = [];               % y-limit for head (m); [] to auto
sty.quick_title = 't = %.2f d (forcing stamp)';  % title format

%% ---------------- 6) Run the model (native forcing; internal substepping) ----------------
[Qout, Qsurf, Qtotal, h, S, diag] = hsB_solver(x, w, ForcingNative, params, opts, do_plot, sty);

% ================= Figure: Flow Duration Curves (baseflow, surface, recharge) =================
% Uses coarse (forcing-stamp) outputs: Qout, Qsurf, ForcingNative.N_mps
% and geometry (w, x/xe or diag.dx). Plots in mm/day with Weibull exceedance.

% ---- Plan area & cell widths ----
if isfield(diag,'Ahs');      Ahs = diag.Ahs; else
    % fallback if diag.Ahs missing
    if exist('xe','var')==1
        dx_local = diff(xe(:));
    else
        % rebuild edges from centers
        xe_tmp = zeros(numel(x)+1,1);
        xe_tmp(1)       = 0.0;
        xe_tmp(2:end-1) = 0.5*(x(1:end-1)+x(2:end));
        xe_tmp(end)     = x(end) + 0.5*(x(end)-x(end-1));
        dx_local = diff(xe_tmp);
    end
    Ahs = sum(w(:).*dx_local(:));
end
if isfield(diag,'dx')
    dx_local = diag.dx(:);
elseif exist('dx_local','var')~=1
    dx_local = diff(xe(:));
end

% ---- Discharges at forcing stamps → mm/day ----
Q_base_mmday = (-Qout   / Ahs) * 86400 * 1000;   % baseflow (subsurface), outward-positive magnitude
Q_surf_mmday = ( Qsurf  / Ahs) * 86400 * 1000;   % surface runoff

% ---- Area-averaged recharge at forcing stamps → mm/day ----
N_F = ForcingNative.N_mps;   % [nF×1] or [nF×Nx] in m/s
if size(N_F,2) == numel(w)
    aw = w(:).*dx_local(:);                 % [Nx×1] cell areas per unit length
    N_area = (N_F * aw) / Ahs;              % [nF×1] m/s
else
    N_area = N_F(:);                         % [nF×1] m/s (already uniform)
end
N_mmday = N_area * 86400 * 1000;             % [nF×1] mm/day

% ---- Weibull plotting positions (exceedance %) ----
weibull = @(n) 100*(( (1:n)' - 0.44 ) ./ (n + 0.12));

% floor so semilogy works when zeros occur
yfloor = 1e-6;  % mm/day

% Baseflow FDC
y_base = max(Q_base_mmday(:), yfloor);
y_base = sort(y_base,'descend');
p_base = weibull(numel(y_base));

% Surface FDC
y_surf = max(Q_surf_mmday(:), yfloor);
y_surf = sort(y_surf,'descend');
p_surf = weibull(numel(y_surf));

% Recharge FDC (only positive recharge)
maskN  = N_mmday(:) > 0;
if any(maskN)
    y_rech = max(N_mmday(maskN), yfloor);
    y_rech = sort(y_rech,'descend');
    p_rech = weibull(numel(y_rech));
else
    y_rech = []; p_rech = [];
    warning('FDC:Recharge','No positive recharge values found for the recharge FDC.');
end

% ---- Plot (semi-log y) ----
figF = figure('Color','w','Position',[120 120 760 520]);
axF  = axes('Parent',figF);

semilogy(axF, p_base, y_base, '-',  'LineWidth', sty.plotlw, 'Color', sty.c.teal); hold(axF,'on');
semilogy(axF, p_surf, y_surf, '--', 'LineWidth', sty.plotlw, 'Color', sty.c.red);
if ~isempty(y_rech)
    semilogy(axF, p_rech, y_rech, ':',  'LineWidth', sty.plotlw, 'Color', sty.c.blue);
end
hold(axF,'off');

grid(axF,'on'); box(axF,'on');
set(axF,'TickDir',sty.tickdir,'TickLength',sty.ticklength, ...
        'LineWidth',sty.axlw,'FontSize',sty.fs, ...
        'GridColor',sty.c.gray,'GridAlpha',0.25, ...
        'TickLabelInterpreter','latex');

xlabel(axF,'\textbf{Exceedance probability (\%)}','Interpreter','latex');
ylabel(axF,'\textbf{Rate (mm d$^{-1}$)}','Interpreter','latex');
title(axF,'\textbf{Flow Duration Curves}','Interpreter','latex');

xlim(axF,[0 100]);
ymaxF = 1.2*max([y_base(:); y_surf(:); (y_rech(:))], [], 'omitnan');
if isempty(ymaxF) || ~isfinite(ymaxF), ymaxF = 1; end
ylim(axF,[yfloor, ymaxF]);

if ~isempty(y_rech)
    legend(axF, {'$\mathrm{Baseflow}\ (-Q_{\mathrm{out}})$', ...
                 '$\mathrm{Surface\ flow}\ (Q_{\mathrm{surf}})$', ...
                 '$\mathrm{Recharge}\ (N>0)$'}, ...
           'Interpreter','latex','Location','southwest');
else
    legend(axF, {'$\mathrm{Baseflow}\ (-Q_{\mathrm{out}})$', ...
                 '$\mathrm{Surface\ flow}\ (Q_{\mathrm{surf}})$'}, ...
           'Interpreter','latex','Location','southwest');
end

%% ---------------- MASS BALANCE (forcing-interval accounting) ----------------
% Geometry
dx    = diff(xe(:));                    % [Nx×1] cell widths [m]
Acell = w(:) .* dx;                     % [Nx×1] cell areas [m^2]
Ahs   = sum(Acell);                     % total hillslope area [m^2]

% Storage volume at forcing stamps: V_k = sum_j S(k,j)*dx_j  [m^3]
V_forcing = S * dx;                     % [nF×1], S is [nF×Nx]
nF = numel(ForcingNative.t_days);

% Volumes accumulated by the wrapper over each forcing interval k=1..nF-1
Rvol_k   = diag.forcing_Rvol(:);        % [nF-1×1] recharge volume in interval k
Outvol_k = diag.forcing_Outvol(:);      % [nF-1×1] outward-positive discharge volume
dV_k     = diag.forcing_dV(:);          % [nF-1×1] storage change actually taken in substeps

% Independent coarse check using V at stamps:
dV_from_V = diff(V_forcing);            % [nF-1×1]

% Step-wise mass-balance residuals at forcing resolution:
% res_k = (V_{k+1}-V_k) - (R_k - Out_k)
res_k = dV_from_V - (Rvol_k - Outvol_k);

% Diagnostics
res_abs_max   = max(abs(res_k));
res_rms       = sqrt(mean(res_k.^2));
inflow_total  = sum(Rvol_k);
outflow_total = sum(Outvol_k);
net_storage   = V_forcing(end) - V_forcing(1);

fprintf('\n--- Mass balance (forcing intervals) ---\n');
fprintf('Intervals: %d\n', nF-1);
fprintf('Total inflow   R = %.3e m^3\n', inflow_total);
fprintf('Total outflow  O = %.3e m^3\n', outflow_total);
fprintf('Storage change ΔV = %.3e m^3\n', net_storage);
fprintf('Closure (R - O - ΔV) = %.3e m^3\n', inflow_total - outflow_total - net_storage);
fprintf('Residual per interval:  max = %.3e m^3,  RMS = %.3e m^3\n', res_abs_max, res_rms);

% (Optional) attach to diag for saving/plotting
diag2 = struct();
diag2.dx                 = dx;
diag2.Acell              = Acell;
diag2.Ahs                = Ahs;

% time series at forcing stamps
diag2.t_forcing_days     = ForcingNative.t_days(:);
diag2.V_forcing          = V_forcing(:);      % [nF×1] storage volume at each stamp

% per-interval (k=1..nF-1) volumes & checks
diag2.Rvol_interval      = Rvol_k(:);         % recharge volume per forcing interval
diag2.Outvol_interval    = Outvol_k(:);       % outward-positive discharge volume
diag2.dV_from_V_interval = dV_from_V(:);      % ΔV computed from V_forcing
diag2.residual_interval  = res_k(:);          % mass-balance residual per interval

% summary stats
diag2.res_abs_max        = res_abs_max;
diag2.res_rms            = res_rms;
%% ---------------- 7) Quick textual summary ----------------
nF = numel(ForcingNative.t_days);
fprintf('\nRun complete: forcing steps=%d, Nx=%d, internal steps=%d\n', ...
    nF, numel(x), diag.n_internal_steps);
if isfield(diag,'Ahs')
    fprintf('Hillslope area used by solver: %.3e m^2\n', diag.Ahs);
end
if isfield(diag,'overflow_steps') && diag.overflow_steps > 0
    fprintf('Saturation-excess in %d/%d forcing stamps (max Q_surf = %.3g m^3/s)\n', ...
        diag.overflow_steps, nF, diag.max_Qsurf);
end

% -------------------------------------------------------------------------
% Helper: ghost-padding for (xH,H)
% -------------------------------------------------------------------------
function [xHg, Hg, dxH, Lphys] = add_ghosts_to_width_function(xH, H)
% ADD_GHOSTS_TO_WIDTH_FUNCTION  Pad (xH,H) with one ghost bin at each end.
% Assumptions:
%   - xH are bin centers with (approximately) uniform spacing dxH.
%   - The physical domain spans [0, Lphys], with Lphys = xH(end) + dxH/2.
%   - Ghost centers sit at: -dxH/2  and  Lphys + dxH/2.
%   - Ghost widths reuse boundary widths: H(1) at left, H(end) at right.
%
% Returns:
%   xHg  : [1×(n+2)] centers including ghosts
%   Hg   : [1×(n+2)] widths including ghosts
%   dxH  : inferred bin thickness from xH
%   Lphys: physical domain length implied by (xH,H)

    xH = xH(:).'; H = H(:).';                 % row vectors
    if numel(xH) < 2, error('xH must have ≥ 2 bins.'); end

    d = diff(xH);
    dxH = median(d);
    if any(abs(d - dxH) > 1e-6*dxH)
        warning('add_ghosts_to_width_function: xH not perfectly uniform; using median dxH.');
    end

    % physical domain (0 … Lphys) if xH are centers
    Lphys = xH(end) + dxH/2;

    % ghost centers & widths
    xLghost = -dxH/2;
    xRghost =  Lphys + dxH/2;
    HLghost = H(1);
    HRghost = H(end);

    xHg = [xLghost, xH, xRghost];
    Hg  = [HLghost, H,  HRghost];
end
