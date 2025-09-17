% -------------------------------------------------------------------------
% Driver script for the Hillslope-Storage Boussinesq (HSB) model
% -------------------------------------------------------------------------
% Author: Marcus Nobrega Gomes Junior, PhD
%
% University of Arizona
% Department of Hydrology and Atmospheric Sciences
% Sep/2025
%
% Purpose:
%   Build inputs (grid, widths, recharge, physical parameters, solver opts),
%   then run the HSB solver and optionally plot results.
%
% What this script does:
%   1) Defines a 1-D hillslope grid (clustered near the outlet) and a width
%      profile w(x).
%   2) Sets the simulation clock (dt, Nt) and recharge N(t) (uniform in
%      space by default, but can be extended to N(t,x)).
%   3) Sets physical parameters (k, f, iota, S0, D).
%   4) Sets solver options (theta, Picard damping & tolerances, etc.).
%   5) Optionally defines plotting style (line widths, colors, colormap).
%   6) Calls:  [Qout, Qsurf, Qtotal, h, S, diag] = hsb_solver(...)
%
% Key notes:
%   • Units: Recharge N must be [m s^-1]. Use (mm/day) / 1000 / 86400.
%   • If you want space–time recharge, build N as [Nt×Nx] instead of [Nt×1].
%   • Set opts.no_recharge_outlet=true to zero recharge in the first
%     micro-cell at the seepage face to avoid spurious spikes.
%   • The solver returns:
%       Qout   [Nt×1]  subsurface outlet flux (negative = out of domain)
%       Qsurf  [Nt×1]  saturation-excess overland flow (≥0)
%       Qtotal [Nt×1]  outward-positive total discharge ( = -Qout + Qsurf )
%       h      [Nt×Nx] head above bedrock [m]
%       S      [Nt×Nx] storage per unit-x (= f w h) [m^2]
%       diag   struct  diagnostics (mass balance residuals, dx, xe, area, …)
%
% Requirements:
%   • The function file containing: hsb_solver(x,w,N,dt,Nt,params,opts,do_plot,sty)
%     must be on the MATLAB path. That function handles the model time-stepping
%     and (optionally) all figures based on do_plot and sty.
% -------------------------------------------------------------------------

clear variables; close all; clc;

%% ---------------- 0) Preprocess Inputs ----------------
run('hsB_PreProcessing.m')


%% ---------------- 1) Grid (clustered near outlet) ----------------
L   = metrics.max_hillslope_length_m;             % hillslope length [m]
Nx  = length(xH);               % number of cell centers

M   = 3;                % how many clustered cells near outlet
r   = 1.15;             % geometric growth factor for cluster
d0  = 5;                % first clustered cell width [m]

dx_geo = d0 * r.^(0:M-1);                 % clustered widths
len_geo = sum(dx_geo);
dx_uni  = (L - len_geo) / (Nx - M);       % remaining cells uniform
dx_all  = [dx_geo, dx_uni*ones(1, Nx-M)];

xe = [0; cumsum(dx_all(:))];              % cell edges  (0…L)
x  = 0.5 * (xe(1:end-1) + xe(2:end));     % cell centers

%% ---------------- 2) Width profile w(x) ----------------
w = xH;                                      % positive width profile pm]


%% ---------------- 3) Time & recharge ----------------
dt_hours = 1;                 % time step [hour]
dt       = dt_hours*3600;     % [s]

T_days   = 730;               % total duration = 2 years
Nt       = round(T_days*86400/dt);      % number of steps
tvec     = (0:Nt-1)'*dt/86400;          % time vector [days] (handy for checks)

% --- Recharge (time-only by default). For space–time, create Nt×Nx. ---
N_mm_per_day_peak = 1;        % peak recharge [mm/day]
dur_hours         = 365*24;   % recharge duration [h]

N0         = zeros(Nt,1);
N0(1:dur_hours) = N_mm_per_day_peak;    % simple block/pulse
N          = (N0/1000)/86400;           % convert mm/day -> m/s, Nt×1

% Example: to switch to space–time recharge, do something like:
% N_spacetime = repmat(N,1,Nx);               % start from uniform
% N_spacetime(:, round(Nx/2):end) = 1.5*N_spacetime(:, round(Nx/2):end);
% N = N_spacetime;                            % now Nt×Nx

% Alternatively, you can simply input a .mat file or a .txt file

%% ---------------- 4) Physical parameters ----------------
params.k     = 1e-3;          % saturated hydraulic conductivity [m/s]
params.f     = 0.05;          % drainable porosity [-]
params.iota  = 1*pi/180;      % bedrock slope angle [rad]
params.S0    = zeros(Nx,1);   % initial storage per unit-x [m^2] (= f*w*h)
params.D     = 3.0;           % aquifer thickness [m] (finite => overflow possible)

%% ---------------- 5) Solver options ----------------
opts.theta      = 1.0;        % 1.0 = fully implicit (stable & robust)
opts.omega      = 0.5;        % Picard relaxation (0.4–0.7 typical)
opts.picard_max = 40;         % max Picard iterations per step
opts.picard_tol = 1e-10;      % L∞ tolerance on Picard increment (m^2)
opts.safeguard  = true;       % enforce S>=0
opts.mass_check = true;       % compute mass-balance residuals
opts.plot_live  = true;       % light live h(x) plot while integrating
opts.verbose    = true;       % print progress each step block
opts.no_recharge_outlet = true; % set true to force N(:,1)=0

%% ---------------- 6) Plotting control & style ----------------
do_plot = true;   % set false to skip all post-run figures

% Optional style overrides (pass empty struct to use function defaults)
sty = struct();
sty.axlw        = 2.5;                 % axes/tick thickness
sty.plotlw      = 2.5;                 % line thickness
sty.tickdir     = 'in';                % ticks inside
sty.ticklength  = [0.02 0.02];
sty.fs          = 13;                  % base font size

% Colors
sty.c.teal      = [40 141 141]/255;    % for discharge line
sty.c.red       = [159 0 0]/255;       % not used here, available
sty.c.gray      = [72 72 72]/255;      % grid/axes
sty.c.blue      = [0 114 178]/255;     % recharge bars

% Colormap for h(x,t)
sty.cmap_name   = 'turbo';             % e.g., 'turbo' | 'parula'
sty.cmap_levels = 21;                  % discrete levels (>=2)
sty.cmap_clim   = [];                  % [] = auto; or [cmin cmax]

% How many water-table profiles to overlay in the “profiles” plot
sty.n_profiles  = 5;

%% ---------------- 7) Run the model ----------------
[Qout, Qsurf, Qtotal, h, S, diag] = hsB_solver(x, w, N, dt, Nt, params, opts, do_plot, sty);

% (Optional) quick textual summary
fprintf('\nRun complete: Nt=%d, Nx=%d, dt=%.0fs, area=%.2e m^2\n', Nt, numel(x), dt, diag.Ahs);
if isfield(diag,'overflow_steps') && diag.overflow_steps > 0
    fprintf('Saturation-excess occurred in %d/%d steps (max Q_surf = %.3g m^3/s)\n', ...
        diag.overflow_steps, Nt, diag.max_Qsurf);
end
if isfield(diag,'mass_residual')
    fprintf('Mean |MB residual| per step: %.3e m^3\n', mean(abs(diag.mass_residual),'omitnan'));
end

% All outputs (Qout, Qsurf, Qtotal, h, S, diag) remain in the workspace.
% ---------------- 8) Access outputs in workspace ----------------
% Examples:
%   max_outflow_mmday = max((-Qout/diag.Ahs)*86400*1000);
%   any_overflow = any(Qsurf>0);


