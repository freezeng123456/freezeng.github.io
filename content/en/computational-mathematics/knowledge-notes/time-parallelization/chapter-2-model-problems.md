---
title: "Chapter 2: Model Problems"
description: Heat, advection-diffusion, Burgers, and wave equations as a dissipation-to-propagation spectrum
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-2-model-problems
tags:
  - parallel-in-time
  - PDE
---

The four model equations form a controlled progression from strong dissipation to persistent propagation. This progression exposes which components of a coarse temporal model are required for convergence.

## 2.1 Heat equation

$$
u_t-\nu\Delta u=f.
$$

A Fourier mode $e^{ikx}$ decays as $e^{-\nu k^2t}$. High-frequency error is removed first, so coarse grids and coarse time steps can often approximate the long-time dynamics. The heat equation is therefore a favorable reference problem for Parareal, MGRiT, and STMG.

Three effects must be separated in an implementation:

- **Discretization plateau:** further PinT iterations are not useful after the iteration error falls below the fine discretization error.
- **Stability versus accuracy:** backward Euler is strongly damping but first order; the trapezoidal rule has different amplitude and phase properties.
- **Spatial coarsening:** temporal coarsening alone does not remove every spatial high-frequency mode.

## 2.2 Advection-diffusion equation

$$
u_t+a u_x-\nu u_{xx}=f.
$$

The Fourier symbol contains both phase, $e^{-iakt}$, and damping, $e^{-\nu k^2t}$. As the Péclet number increases or $\nu$ decreases, phase accuracy becomes more important than additional numerical damping. A stable coarse propagator may still be ineffective if it travels at the wrong numerical speed.

For a mode with wavenumber $k$, the ratio

$$
\chi_k=\frac{|ak|}{\nu k^2}=\frac{|a|}{\nu|k|}
$$

measures the competition between propagation and dissipation. The low-frequency modes have the largest $\chi_k$ and are also the most persistent in time.

### Recomputed solution

The formal Python experiment uses a source-free Dirichlet problem with $\Delta t=\Delta x=10^{-3}$, $T=3$, $\nu=5\times10^{-4}$, and advection speed $a=1$. The initial condition is $\sin^2(8\pi(1-x)^2)$. Its final $L^\infty$ norm is $7.022\times10^{-99}$.

![[assets/pint/model-advection-diffusion.png]]

The near-zero final norm is specific to this dissipative, outflow-dominated Dirichlet experiment. It should not be transferred to periodic or forced variants of the equation.

## 2.3 Viscous Burgers equation

$$
u_t+u u_x-\nu u_{xx}=f.
$$

The propagation speed now depends on the solution. Fine and coarse propagators must therefore reconcile nonlinearity, phase, and diffusion at every PinT iteration. When viscosity is small, steep gradients couple spatial discretization error, nonlinear-solver error, and temporal iteration error.

Parareal retains the update

$$
U_{n+1}^{k+1}
=G(U_n^{k+1})+F(U_n^k)-G(U_n^k),
$$

but $F-G$ is no longer a fixed linear correction.

### Recomputed solution

The formal experiment uses $\Delta t=\Delta x=1/400$, $T=3$, and $\nu=5\times10^{-4}$ with homogeneous Dirichlet boundaries. The solution reaches a maximum value of $1.045940$ and has final $L^\infty$ norm $0.325871$.

![[assets/pint/model-burgers.png]]

## 2.4 Second-order wave equation

$$
u_{tt}-c^2\Delta u=f.
$$

After conversion to a first-order system, the ideal propagation eigenvalues lie on or near the unit circle. Coarse phase error accumulates because the dynamics do not remove it. Useful approaches often employ characteristic variables, waveform relaxation, exponential propagation, diagonalization, or a phase-aware coarse model.

### Recomputed solution

The source-free formal experiment uses the trapezoidal rule with $\Delta t=\Delta x=1/400$, $T=3$, and spatial coefficient $0.2$. It starts from the same oscillatory displacement as the previous two experiments and zero initial velocity. The final displacement has $L^\infty$ norm $0.948217$.

![[assets/pint/model-wave.png]]

The persistent amplitude contrasts with the nearly extinguished advection-diffusion solution. This comparison illustrates the temporal nonlocality that makes wave propagation difficult for dissipative coarse-grid corrections.

## Numerical summary

| Formal experiment   | Grid and horizon                   |                      Reported metric |
| ------------------- | ---------------------------------- | -----------------------------------: |
| advection-diffusion | $\Delta x=\Delta t=10^{-3}$, $T=3$ | final $L^\infty=7.022\times10^{-99}$ |
| viscous Burgers     | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.325871$ |
| wave equation       | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.948217$ |

These experiments visualize the model dynamics. They are not parallel speedup measurements.
