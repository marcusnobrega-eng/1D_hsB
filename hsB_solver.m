function [Qout, Qsurf, Qtotal, h, S, diag] = hsB_solver(x, w, N, dt, Nt, params, opts, do_plot, sty)
%HSB_SOLVER  Run the Hillslope-Storage Boussinesq (HSB) model with given inputs.
%   [Qout,Qsurf,Qtotal,h,S,diag] = hsb_solver(x,w,N,dt,Nt,params,opts,do_plot,sty)
%   runs the semi-implicit finite-volume HSB solver on the provided grid,
%   width profile, and recharge. Optionally plots results if do_plot==true.
%
%   INPUTS
%     x      [Nx×1]  strictly increasing cell-center coordinates [m]
%     w      [Nx×1]  planform width at centers [m] (>0)
%     N      [Nt×1] or [Nt×Nx] recharge [m s^-1] (time-only or space–time)
%     dt     scalar  time step [s]
%     Nt     scalar  number of time steps
%     params struct  fields:
%            .k     [m s^-1] saturated hydraulic conductivity
%            .f     [-]      drainable porosity
%            .iota  [rad]    bedrock slope angle
%            .S0    [Nx×1]   initial storage per unit-x (= f w h) [m^2]
%            .D     [m]      aquifer thickness (cap). If finite => overflow possible
%     opts   struct  optional:
%            .theta, .omega, .picard_max, .picard_tol, .safeguard,
%            .verbose, .mass_check, .plot_live, .no_recharge_outlet
%     do_plot logical  if true, produce publication-style plots (default: false)
%     sty    struct    plotting style (optional). If omitted, sensible defaults used.
%
%   OUTPUTS
%     Qout    [Nt×1]  subsurface discharge at outlet face (x=0) [m^3 s^-1]
%                     NOTE: model-internal sign: negative = outward
%     Qsurf   [Nt×1]  saturation-excess overland flow (outlet-equiv) [m^3 s^-1] (≥0)
%     Qtotal  [Nt×1]  outward-positive total discharge = -Qout + Qsurf [m^3 s^-1]
%     h       [Nt×Nx] head above bedrock [m]
%     S       [Nt×Nx] storage per unit-x (= f w h) [m^2]
%     diag    struct  diagnostics: mass_residual, dx, xe, Ahs, Qsurf, etc.
%
%   NOTES
%     - Set opts.no_recharge_outlet=true to zero recharge in the first micro-cell.
%     - Live plotting inside the solver occurs only if opts.plot_live && do_plot.

% ------------------ defaults for plotting args ------------------
if nargin < 8 || isempty(do_plot)
    do_plot = false; 
end
if nargin < 9 || isempty(sty)
    sty.axlw        = 2.5;                 % axes & tick thickness
    sty.plotlw      = 2.5;                 % line thickness
    sty.tickdir     = 'in';                % ticks inside
    sty.ticklength  = [0.02 0.02];
    sty.fs          = 13;                  % base font size
    sty.c.teal      = [40 141 141]/255;    % #288D8D
    sty.c.red       = [159 0 0]/255;       % #9F0000
    sty.c.gray      = [72 72 72]/255;      % dark gray
    sty.c.blue      = [0 114 178]/255;     % nice blue for recharge bars
    sty.cmap_name   = 'turbo';
    sty.cmap_levels = 21;
    sty.cmap_clim   = [];
    sty.n_profiles  = 5;
end

% --------------- run core solver (no plotting here) ---------------
opts_local = opts;
opts_local.plot_live = isfield(opts,'plot_live') && opts.plot_live && do_plot;

[Qout, Qsurf, Qtotal, h, S, diag] = hsb_core(x, w, N, dt, Nt, params, opts_local);

% --------------- optional plotting ---------------
if do_plot
    tvec = (0:Nt-1)'*dt/86400;  % days for plotting

    % ----- helper for area-weighted recharge (mm/day) -----
    Ahs      = diag.Ahs;                        % hillslope plan area [m^2]
    Q_mmday  = (-Qout / Ahs) * 86400 * 1000;    % outward-positive, mm/day

    if size(N,2) == numel(w)
        aw     = w(:).*diag.dx(:);              % area weights [m^2 per m in x]
        N_area = (N * aw) / Ahs;               % [Nt×1] m/s
    else
        N_area = N;                             % [Nt×1] m/s
    end
    N_mmday = N_area * 86400 * 1000;           % mm/day

    % ================= Figure 1: Q (left) + recharge bars (right) ================
    figure('Color','w','Position',[120 90 1100 450]);
    tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % ---- Left tile: discharge line (left y-axis) + recharge bars (right y-axis)
    ax1 = nexttile(tl);

    yyaxis(ax1,'left');
    plot(ax1, tvec, Q_mmday, 'LineWidth', sty.plotlw, 'Color', sty.c.teal);
    grid(ax1,'on'); box(ax1,'on');
    set(ax1,'TickDir',sty.tickdir,'TickLength',sty.ticklength, ...
            'LineWidth',sty.axlw,'FontSize',sty.fs, ...
            'GridColor',sty.c.gray,'GridAlpha',0.25, ...
            'TickLabelInterpreter','latex');
    xlabel(ax1,'\textbf{Time (days)}','Interpreter','latex');
    ylabel(ax1,'\textbf{$Q_{\mathrm{out}}$ (mm d$^{-1}$)}','Interpreter','latex');
    qmax = max(Q_mmday(:)); ylim(ax1,[0, 1.5*max(qmax, eps)]);

    yyaxis(ax1,'right');
    bar(ax1, tvec, N_mmday, 'FaceColor', sty.c.blue, 'EdgeColor','none', 'BarWidth', 1);
    set(ax1,'YDir','reverse');
    ylabel(ax1,'\textbf{$N$ (mm d$^{-1}$)}','Interpreter','latex');
    nmax = max(N_mmday(:)); ylim(ax1,[0, 4*max(nmax, eps)]);
    ax1.YAxis(1).Color = sty.c.teal;
    ax1.YAxis(2).Color = sty.c.blue;
    title(ax1,'\textbf{Outlet hydrograph \& recharge}','Interpreter','latex');

    % ================= Right tile: water-table profiles ==========================
    ax2 = nexttile(tl);
    npr  = max(2, round(sty.n_profiles));
    idx  = unique(round(linspace(1, Nt, npr)));

    if exist(sty.cmap_name,'file')==2 || exist(sty.cmap_name,'builtin')==5
        cmapP = feval(sty.cmap_name, numel(idx));
    else
        cmapP = lines(numel(idx));
    end
    hold(ax2,'on');
    for k = 1:numel(idx)
        plot(ax2, x, h(idx(k),:), 'LineWidth', sty.plotlw, 'Color', cmapP(k,:));
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
    ylim(ax2,[0, max(h(:))*1.05]);
    legLabels = arrayfun(@(tt) sprintf('$t=%.2f$ d', tt), tvec(idx), 'UniformOutput', false);
    legend(ax2, legLabels, 'Interpreter','latex','Location','best');

    % ================= Figure 2: h(x,t) image ===================================
    fig2 = figure('Color','w','Position',[120 560 880 420]);
    ax3  = axes('Parent',fig2);
    imagesc(ax3, x, tvec, h); set(ax3,'YDir','normal');
    grid(ax3,'on'); box(ax3,'on');
    set(ax3,'TickDir',sty.tickdir,'TickLength',sty.ticklength, ...
            'LineWidth',sty.axlw,'FontSize',sty.fs, ...
            'GridColor',sty.c.gray,'GridAlpha',0.15, ...
            'TickLabelInterpreter','latex');
    xlabel(ax3,'\textbf{x (m)}','Interpreter','latex');
    ylabel(ax3,'\textbf{time (days)}','Interpreter','latex');
    title(ax3,'\textbf{$h(x,t)$}','Interpreter','latex');

    % discrete colormap + optional fixed limits
    if ~(exist(sty.cmap_name,'file')==2 || exist(sty.cmap_name,'builtin')==5)
        warning('Colormap "%s" not found. Using "parula".', sty.cmap_name);
        sty.cmap_name = 'parula';
    end
    cmap = feval(sty.cmap_name, max(2, round(sty.cmap_levels)));
    colormap(ax3, cmap);
    if ~isempty(sty.cmap_clim) && numel(sty.cmap_clim)==2
        caxis(ax3, sort(sty.cmap_clim));
    end
    cb = colorbar(ax3);
    set(cb,'LineWidth',sty.axlw,'TickDirection',sty.tickdir,'TickLabelInterpreter','latex');

    % ================= Figure 3: mass-balance residual ==========================
    if isfield(opts,'mass_check') && opts.mass_check
        fig3 = figure('Color','w','Position',[1020 560 520 420]);
        ax4  = axes('Parent',fig3);
        plot(ax4, tvec, diag.mass_residual, 'LineWidth', sty.plotlw, 'Color', sty.c.teal);
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
end
% =============================== CORE SOLVER ===============================
function [Qout, Qsurf, Qtotal, h_store, S_store, diag] = hsb_core(x, w, N, dt, Nt, params, opts)
% HSB core: semi-implicit FV + Picard, stabilized fluxes, overflow sink.

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
theta      = get_opt(opts,'theta',1.0);      % fully implicit by default
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
Qsurf   = zeros(Nt,1);    % saturation-excess overland flow (m^3 s^-1)
h_store = zeros(Nt,Nx);
S_store = zeros(Nt,Nx);
mass_residual = nan(Nt,1);

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
    Nn    = N(n,:).';            % recharge at step n

    % Picard initialization + per-step damping reset
    omega     = omega0;
    incr_prev = Inf;
    S_iter    = S_old;

    % -------------------- Picard iterations --------------------
    for m = 1:picard_max
        % state used to evaluate fluxes at n+theta
        S_theta = (1-theta)*S_old + theta*S_iter;

        % fluxes at faces using stabilized scheme
        h_cent = S_theta ./ max(f*w, eps);
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
    h_now = S ./ max(f*w, eps);
    h_store(n,:) = h_now.';
    S_store(n,:) = S.';

    % outlet discharge from accepted state
    Qface    = compute_face_fluxes(x, xe, w, h_now, S, k, f, cosI, sinI);
    Qout(n)  = Qface(1);     % NOTE: negative = outward

    % ---- Compute saturation-excess overflow consistently (using accepted operator)
    divQ2   = -(Qface(2:end) - Qface(1:end-1)) ./ dx;
    rawS    = S_old + dt * (divQ2 + Nn .* w);      % [Nx×1], could exceed caps
    if safeguard, rawSpos = max(rawS, 0); else, rawSpos = rawS; end
    Sclip   = min(rawSpos, Smax);
    excess  = max(rawSpos - Sclip, 0);             % per-cell overflow [m^2]
    Qsurf(n)= sum(excess .* dx) / dt;              % outlet-equivalent overland flow [m^3 s^-1]

    % ---- Mass-balance residual (should be ~0) with explicit overflow sink
    if mass_check
        R_excess = excess / dt;                    % [m^2 s^-1] per unit-x
        res_cell = (S - S_old)/dt - (divQ2 + Nn .* w - R_excess);
        mass_residual(n) = sum(res_cell .* dx);    % [m^3]
    end

    % light live plot
    if plot_live && mod(n, max(1,round(Nt/100)))==0
        quick_plot_step(x, w, h_now, n, D, dt);
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
    'overflow_mask',   Qsurf > 0, ...
    'overflow_steps',  nnz(Qsurf > 0), ...
    'max_Qsurf',       max(Qsurf, [], 'omitnan') ...
);

% Informative warning if any saturation-excess occurred
if any(Qsurf > 0)
    warning('HSB:Overflow', ...
        'Saturation-excess overland flow occurred in %d/%d steps (max Q_surf = %.3g m^3 s^{-1}).', ...
        diag.overflow_steps, Nt, diag.max_Qsurf);
end
end

% =============================== HELPERS ===============================
function [xe, dx] = build_edges_from_centers(x)
Nx = numel(x);
xe = zeros(Nx+1,1);
xe(1)     = 0.0;                                 % outlet face
xe(2:Nx)  = 0.5*(x(1:Nx-1)+x(2:Nx));             % interior midpoints
xe(Nx+1)  = x(end) + 0.5*(x(end)-x(end-1));      % extrapolate right edge
dx = diff(xe);
end

function Qface = compute_face_fluxes(xc, xe, w, h_c, S_c, k, f, cosI, sinI)
% Monotone face flux = centered diffusion (harmonic S) + upwinded drift
% Q = Qdiff + Qadv, with small Lax–Friedrichs viscosity on the drift.
Nx = numel(xc);
Qface = zeros(Nx+1,1);

Kd_over_f = (k/f);

% Small numerical film near the seepage face to prevent S_f→0
h_floor     = 1e-4;                                % 0.1–1 mm
S_floor_L   = f * w(1) * h_floor;
S_floor_int = f * 0.5*(w(1:end-1)+w(2:end)) * h_floor;

% ---- outlet face (x=0): Dirichlet h(0)=0
dxL   = max(xc(1) - xe(1), eps);
gradL = (h_c(1) - 0) / dxL;
Sf    = harmonic_mean(max(S_c(1),S_floor_L), max(S_c(1),S_floor_L));
Qdiff = - Kd_over_f * cosI * Sf * gradL;

Sup   = S_c(1);                                   % upwind for drift (downslope)
Qadv  = - Kd_over_f * sinI * Sup;

a     = abs(Kd_over_f * sinI);                    % LF viscosity speed
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
    a     = abs(Kd_over_f * sinI);
    Qnum  = 0.5*a*(Sup - Sdn);

    Qface(i) = Qdiff + Qadv + Qnum;
end

% ---- divide face: no-flow
Qface(Nx+1) = 0.0;
end

function Sh = harmonic_mean(a,b)
a = max(a,0); b = max(b,0);
Sh = 2*a.*b ./ max(a+b, eps);
end

function quick_plot_step(x, w, h, n, D, dt)
% Styled live view (consistent with driver’s figures)
sty.axlw        = 2.5;  sty.plotlw = 2.5;
sty.tickdir     = 'in'; sty.ticklength = [0.02 0.02];
sty.fs          = 13;
sty.c.teal      = [40, 141, 141]/255;
sty.c.gray      = [72, 72, 72]/255;

if ~ishandle(9999), figure(9999); clf; end
fig = figure(9999); set(fig,'Color','w','Position',[80 80 880 360]);
tl = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

ax1 = nexttile(tl);
plot(ax1, x, h, '-', 'LineWidth', sty.plotlw, 'Color', sty.c.teal);
grid(ax1,'on'); box(ax1,'on');
set(ax1,'TickDir',sty.tickdir,'TickLength',sty.ticklength, ...
        'LineWidth',sty.axlw,'FontSize',sty.fs, ...
        'GridColor',sty.c.gray,'GridAlpha',0.25, ...
        'TickLabelInterpreter','latex');
xlabel(ax1,'\textbf{x (m)}','Interpreter','latex');
ylabel(ax1,'\textbf{h (m)}','Interpreter','latex');
title(ax1, sprintf('\\textbf{h(x), step %d}', n), 'Interpreter','latex');
if isfinite(D), ylim(ax1,[0, D]); else, ylim(ax1,[0, max(h)*1.05]); end

ax2 = nexttile(tl);
plot(ax2, x, w, '-', 'LineWidth', sty.plotlw, 'Color', sty.c.gray);
grid(ax2,'on'); box(ax2,'on');
set(ax2,'TickDir',sty.tickdir,'TickLength',sty.ticklength, ...
        'LineWidth',sty.axlw,'FontSize',sty.fs, ...
        'GridColor',sty.c.gray,'GridAlpha',0.25, ...
        'TickLabelInterpreter','latex');
xlabel(ax2,'\textbf{x (m)}','Interpreter','latex');
ylabel(ax2,'\textbf{w (m)}','Interpreter','latex');
title(ax2,'\textbf{Width profile}','Interpreter','latex');

sgtitle(tl, sprintf('\\textbf{$t = %.2f$ days}', n*dt/86400), 'Interpreter','latex');
end

function v = get_opt(s, f, d)
if nargin<3, d = []; end
if ~isempty(s) && isstruct(s) && isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end

