# ⛰️ hsB: Hillslope-Storage Boussinesq Model

A MATLAB implementation of the 1D Hillslope-Storage Boussinesq (HSB) model for simulating subsurface flow along a hillslope with optional saturation-excess overland flow when the water table reaches a finite aquifer thickness.

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)
![](https://img.shields.io/github/issues/marcusnobrega-eng/DRAIN-LID)
![](https://img.shields.io/github/forks/marcusnobrega-eng/DRAIN-LID)
![](https://img.shields.io/github/last-commit/marcusnobrega-eng/DRAIN-LID)

<img width="1563" height="430" alt="image" src="https://github.com/user-attachments/assets/5330e1dd-c076-46d3-a01d-f529d7a7066d" />


---

## 📘 Project Summary

hsB solves a storage form of the Boussinesq equation on a 1D hillslope using a finite-volume, semi-implicit (θ-scheme) Picard solver. 
It supports non-uniform grids, arbitrary planform width w(x), time- and space-varying recharge N(t) or N(x,t), and a finite aquifer thickness DDD that activates saturation-excess overflow.


## Designed for

🏔️ Catchment/hillslope subsurface flow and hydrograph generation

🧮 Mass-conservative scheme with stabilization and diagnostics

🔁 Quick experimentation in MATLAB (no special toolboxes required)


### 🚀 Quick Start

Clone or download this repository.

Open MATLAB in the repo folder.

Edit run_hsB.m to set grid, width, recharge, parameters, and options.

Run:

[Qout, Qsurf, Qtotal, h, S, diag] = hsB_solver(x, w, N, dt, Nt, params, opts, do_plot);
---

## 🔍 Model Capabilities

| Feature                        | Description                                                                   |
| ------------------------------ | ----------------------------------------------------------------------------- |
| 📦 Finite-volume in space      | Control volumes on non-uniform centers; conservative fluxes                   |
| ⏱️ Semi-implicit time stepping | θ-scheme (default θ=1, backward Euler) with Picard relaxation                 |
| 🌊 Flux splitting              | Diffusion (harmonic-mean storage) + slope-driven drift with LF stabilization  |
| 🧭 Boundaries                  | Outlet seepage face: Dirichlet $h(0)=0$; Divide: no-flow                      |
| 📈 Recharge forcing            | $N(t)$ (Nt×1) or $N(x,t)$ (Nt×Nx); **enter as mm/day in driver**              |
| 📈 Arbitrary width             | Any positive $w(x)$, including convergent/divergent planforms                 |
| 🧱 Finite thickness $D$        | Enforces $0\le h\le D$; exceedance yields saturation-excess $Q_{\text{surf}}$ |
| ✅ Mass balance                 | Per-step residual (m³), diagnostics, optional live plot                       |
| 🖼️ Visualization              | Publication-style plots (thick frames, LaTeX labels, discrete colormap)       |


## 📥 Inputs (set in run_hsB.m)
•	Grid

o	L [m]: hillslope length

o	x [Nx×1]: cell centers (non-uniform OK; script can cluster near outlet)

o	w [Nx×1, m]: planform width at centers (positive)

•	Time & Recharge

o	dt [s], Nt [–]; tvec built for plotting

o	N: mm/day in driver; converted to m/s internally. Shape Nt×1 (uniform) or Nt×Nx (space–time).

•	Physical Parameters (params)

o	k [m/s], f [–], iota [rad], S0 [Nx×1, m² per unit-x = f·w·h], D [m]

•	Solver Options (opts)

o	theta, omega, picard_max, picard_tol, safeguard, mass_check, plot_live, no_recharge_outlet

•	Plotting

o	do_plot = true/false; optional style struct sty (colors, line widths, colormap levels)

## 📤 Outputs
Both run_hsB and hsB_solver return:

•	Qout [Nt×1] — subsurface discharge at outlet face (model sign: negative = outward) [m³/s]

•	Qsurf [Nt×1] — saturation-excess overland flow (outlet-equivalent) [m³/s]

•	Qtotal [Nt×1] — outward-positive total discharge: Qtotal = -Qout + Qsurf [m³/s]

•	h [Nt×Nx] — water table height above bedrock [m]

•	S [Nt×Nx] — storage per unit-x (S = f·w·h) [m²]

•	diag — diagnostics (dx, xe, Ahs, mass_residual, Qsurf, overflow stats, …)

Areal hydrograph conversion to mm/day:


## 📚 Documentation
A full user manual is included in the repository, covering:

📖 Mathematical formulation and boundary condition types

💻 Numerical methods (discretization, solver)


## 👤 Developer
**Marcus N. Gomes Jr.**  
Postdoctoral Researcher II, University of Arizona  
📧 Email: [marcusnobrega.engcivil@gmail.com](mailto:marcusnobrega.engcivil@gmail.com)  
🌐 Website: [marcusnobrega-eng.github.io/profile](https://marcusnobrega-eng.github.io/profile)  
📄 CV: [Download PDF](https://marcusnobrega-eng.github.io/profile/files/CV___Marcus_N__Gomes_Jr_.pdf)  
🧪 ORCID: [0000-0002-8250-8195](https://orcid.org/0000-0002-8250-8195)  
🐙 GitHub: [@marcusnobrega-eng](https://github.com/marcusnobrega-eng)
