---
title: "Chapter 2: Model Problems"
description: A paragraph-level study of temporal locality and long-range propagation in four model PDEs
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-2-model-problems
tags:
  - parallel-in-time
  - PDE
---

> [!note] Reading scope
> This chapter follows Section 2 of the paper (pp. 388–396). Sections 2.1–2.5 explain the models, parameters, boundary conditions, and Figures 2.1–2.4 in source order. Section 2.6 presents three recomputed experiments from this site. Claims from the paper and results produced by our Python project are identified separately.

## 2.1 Why begin with four model problems?

After spatial discretization, the paper writes a linear PDE as

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol g(t),
\qquad \boldsymbol u(0)=\boldsymbol u_0, \tag{2.1}
$$

and a nonlinear PDE as

$$
\boldsymbol u'(t)=\boldsymbol f(\boldsymbol u(t),t),
\qquad \boldsymbol u(0)=\boldsymbol u_0. \tag{2.2}
$$

The matrix $A$ and vector function $\boldsymbol f$ result from the spatial discretization. All examples use $\Omega=(0,1)$. A one-dimensional domain isolates the temporal mechanism: how much past information the solution forgets, how long fine scales survive, and whether a boundary allows information to leave. Additional spatial dimensions change the cost of each spatial solve but do not remove these mechanisms.

The four equations form a controlled dynamical spectrum. The heat equation is strongly diffusive. Advection–diffusion combines decay with directed transport. Burgers' equation adds a solution-dependent velocity and shocks. The wave equation preserves and reflects detailed structures for long times. This spectrum serves as a common stress test throughout the paper.

![The four model equations span temporal locality to long-range memory](assets/diagrams/pint/en/model-memory-spectrum.svg)

> [!important] Meaning of “favorable for PinT”
> The phrase means that a coarse temporal representation and local corrections may communicate useful information quickly. It does not establish wall-clock speedup, nor does it guarantee robustness for every discretization, boundary condition, and parameter regime.

## 2.2 Heat equation: diffusion creates temporal locality

The first problem is

$$
u_t-\nu u_{xx}=f(x,t),
\qquad (x,t)\in\Omega\times(0,T), \tag{2.3}
$$

with $\nu>0$. Figure 2.1 uses $\nu=0.1$ and the localized source

$$
f(x,t)=e^{-\sigma(x-\frac12)^2}
\Bigl(
e^{-\sigma(t-0.1)^2}+e^{-\sigma(t-0.6)^2}
+e^{-\sigma(t-1.35)^2}+e^{-\sigma(t-1.85)^2}
\Bigr),
\qquad \sigma=200. \tag{2.4}
$$

The source is concentrated around $x=1/2$ and fires near four distinct times. These well-separated space–time events make the fading memory of the solution easy to observe.

### Homogeneous Dirichlet boundaries

With $u(0,t)=u(1,t)=0$, heat can leave the domain. Each pulse in Figure 2.1(a) creates a local hot region that rapidly spreads and decays. A late time interval retains little detailed information about events far in the past; only a small collection of slowly decaying components remains.

This behavior is the temporal locality emphasized in the paper. After a long horizon is partitioned, useful work on later intervals may begin from an inaccurate initial value because diffusion will have suppressed much of its high-frequency error before corrected interface data arrives. This setting favors Nievergelt-style precomputation and parabolic coarse-grid correction.

### Homogeneous Neumann boundaries

With $u_x(0,t)=u_x(1,t)=0$, no heat crosses the boundary. The local structure in Figure 2.1(b) is still smoothed, while the spatial mean remains. When the source is removed and the initial condition is oscillatory, Figure 2.1(c) shows the oscillations decaying toward a spatial constant.

The Neumann problem therefore contains two temporal scales. High frequencies are easy to eliminate locally; the constant mode must be communicated accurately across the complete horizon. A PinT method may handle the first component rapidly while still being limited by the second.

### Periodic boundaries

Periodic boundaries also retain the constant spatial mode. Figure 2.1(d) resembles the Neumann case: oscillatory details disappear quickly, but the mean persists. The PDE class supplies only the first prediction of PinT behavior. Boundary conditions determine which low-frequency information can escape.

### PinT implication

A Fourier mode $e^{ikx}$ of the heat equation decays as $e^{-\nu k^2t}$. Larger $|k|$ means faster forgetting. Parareal, MGRiT, and STMG can exploit this scale separation: coarse levels carry slowly varying components, while smoothing or fine propagation resolves local high-frequency error. Dirichlet boundaries are especially favorable; Neumann and periodic problems demand accurate treatment of near-zero spatial frequencies.

## 2.3 Advection–diffusion: persistent transport as diffusion weakens

The second model is

$$
u_t+u_x-\nu u_{xx}=f(x,t),
\qquad (x,t)\in\Omega\times(0,T). \tag{2.5}
$$

The advection speed is one. Figure 2.2 compares $\nu=1$, $10^{-2}$, and $5\times10^{-4}$. The same signal spreads rapidly in the diffusion-dominated regime and follows a distinct characteristic path when diffusion is small.

### Information exits under Dirichlet boundaries

Figures 2.2(a–c) use homogeneous Dirichlet boundaries and the source in (2.4). At $\nu=1$, diffusion largely hides the directed motion. At $\nu=10^{-2}$, right-leaning signal bands are visible. At $\nu=5\times10^{-4}$, the signal is transported almost along a characteristic and leaves through the right boundary.

Figure 2.2(d) removes the source and uses an oscillatory initial condition. Small diffusion preserves high frequencies longer, yet outflow eventually clears the domain. The paper observes that candidate solutions on sufficiently late intervals can then be precomputed and combined or corrected after the incoming state becomes available. This opportunity comes from eventual outflow, not from small $\nu$ by itself.

### Periodicity exposes the hyperbolic limit

Figures 2.2(e–h) use periodic boundaries. Diffusion still removes detail rapidly at $\nu=1$. At $\nu=10^{-2}$, slanted bands remain recognizable after wrapping around the domain. At $\nu=5\times10^{-4}$, sharp structures traverse the periodic boundary repeatedly, so events from the distant past keep influencing the future.

For wavenumber $k$, the continuous propagator contains phase $e^{-ikt}$ and damping $e^{-\nu k^2t}$. The ratio

$$
\chi_k=\frac{|k|}{\nu k^2}=\frac{1}{\nu|k|}
$$

measures transport relative to diffusion. Small $\nu$ moves more modes toward pure advection. A stable coarse propagator may still be ineffective when its phase speed places the correction at the wrong location.

> [!warning] Designing a meaningful stress test
> A Dirichlet outflow problem can overstate robustness near the transport limit. Periodic boundaries, small viscosity, and a fine-scale initial condition reveal phase error and long-range memory much more clearly.

## 2.4 Viscous Burgers equation: nonlinearity amplifies the challenge

The third model is

$$
u_t+u\,u_x-\nu u_{xx}=f(x,t),
\qquad (x,t)\in\Omega\times(0,T). \tag{2.6}
$$

The transport velocity is now the solution itself. High-amplitude regions move faster, causing different parts of the profile to catch up and form steep fronts. The paper retains the viscosity values and boundary settings used for advection–diffusion so that the linear and nonlinear results can be compared directly.

### Deformation, diffusion, and outflow with Dirichlet boundaries

Figures 2.3(a–d) show strong diffusion at large $\nu$. As viscosity decreases, source pulses deform during transport and develop sharper gradients. Homogeneous Dirichlet boundaries still allow these structures to leave, so very late intervals may gradually lose detailed dependence on the initial events.

### Persistent shocks with periodic boundaries

In Figures 2.3(e–h), periodicity recirculates the signal. Small viscosity both weakens damping and permits the nonlinearity to generate new high spatial frequencies. A phase error now changes the local transport velocity and shock position as well. The difference between fine and coarse propagators depends on the current state and cannot be represented completely by one fixed linear correction.

The plots expose a three-way coupling relevant to PinT: the spatial mesh controls shock resolution, the nonlinear solver controls each local propagation, and the time-parallel iteration reconciles states across intervals. Periodicity with small viscosity is therefore a stricter test than its linear advection–diffusion counterpart.

## 2.5 Second-order wave equation: fine information crosses the whole horizon

The final model is

$$
u_{tt}-c^2u_{xx}=f(x,t),
\qquad (x,t)\in\Omega\times(0,T), \tag{2.7}
$$

with $c^2=0.2$. The paper uses both the localized source (2.4) and an oscillatory initial displacement, always with zero initial velocity. Figure 2.4 compares Dirichlet, Neumann, and periodic boundaries.

### All three boundary conditions preserve long-range influence

Dirichlet boundaries reflect waves with one sign relation, Neumann boundaries produce another reflection relation, and periodic boundaries return a wave through the opposite side. Their detailed behavior differs, yet every case contains fine propagation paths that cross the full time interval. A wavefront created by an early source or initial condition can still affect the state after several reflections and interactions.

The second-order equation transports information along two characteristic directions. In a first-order formulation, ideal propagation modes lie on or near the unit circle. The dynamics do not damp phase error, so a small coarse-grid wave-speed error becomes a visible positional error over time.

This is the central difficulty of hyperbolic PinT. Coarse-grid methods based on dissipation and temporal locality cannot rely on an automatic forgetting mechanism. Effective algorithms must use finite propagation speed, characteristics, exponential propagation, or global time algebra. SWR, ParaExp, and ParaDiag in Chapter 3 embody these alternatives.

## 2.6 Site supplement: three recomputed model solutions

The following Python experiments contrast diffusion with outflow against persistent propagation. They use their own initial data, grids, and source-free settings and are not pixel-level reproductions of Figures 2.2–2.4.

### Advection–diffusion solution

The experiment uses homogeneous Dirichlet boundaries, $\Delta t=\Delta x=10^{-3}$, $T=3$, $\nu=5\times10^{-4}$, unit advection speed, and initial data $\sin^2(8\pi(1-x)^2)$. The final $L^\infty$ norm is $7.022\times10^{-99}$.

![Recomputed space–time solution of the advection–diffusion equation](assets/pint/model-advection-diffusion.svg)

This near-zero final state is specific to the present diffusive outflow setting. A periodic problem would not forget at the same rate.

### Viscous Burgers solution

The experiment uses $\Delta t=\Delta x=1/400$, $T=3$, $\nu=5\times10^{-4}$, and homogeneous Dirichlet boundaries. The maximum over the space–time domain is $1.045940$, and the final $L^\infty$ norm is $0.325871$.

![Recomputed space–time solution of the viscous Burgers equation](assets/pint/model-burgers.svg)

### Wave solution

The source-free experiment applies the trapezoidal rule with $\Delta t=\Delta x=1/400$, $T=3$, and $c^2=0.2$. It uses the same oscillatory initial displacement as the previous experiments and zero initial velocity. The final displacement has $L^\infty$ norm $0.948217$.

![Recomputed space–time solution of the wave equation](assets/pint/model-wave.svg)

### Numerical summary

| Supplemental experiment | Grid and horizon                   |                      Reported metric |
| ----------------------- | ---------------------------------- | -----------------------------------: |
| advection–diffusion     | $\Delta x=\Delta t=10^{-3}$, $T=3$ | final $L^\infty=7.022\times10^{-99}$ |
| viscous Burgers         | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.325871$ |
| wave equation           | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.948217$ |

These values characterize the discrete dynamics. They do not measure PinT convergence or hardware speedup.

## 2.7 Comparing the four models

| Model               | Dominant mechanism                           | Information retained over long times              | Primary PinT requirement                                            |
| ------------------- | -------------------------------------------- | ------------------------------------------------- | ------------------------------------------------------------------- |
| heat                | diffusion                                    | constant and a few low-frequency modes            | communicate slow modes accurately on coarse levels                  |
| advection–diffusion | directed transport and diffusion             | phase and fine scales at small viscosity          | control phase mismatch between fine and coarse propagation          |
| Burgers             | state-dependent transport, shocks, diffusion | shock position, amplitude, and nonlinear velocity | coordinate spatial, nonlinear-solver, and temporal iteration errors |
| wave                | bidirectional transport and reflection       | long-time amplitude, phase, and propagation paths | use characteristics and finite speed, or a global algebraic method  |

Boundary conditions decide whether information can leave. Viscosity sets the lifetime of fine scales. Nonlinearity determines whether fine scales are regenerated. These three questions recur in the analysis of every method that follows.

## 2.8 Source coverage audit

| Source location                   | This page | Material covered                                                                                                               |
| --------------------------------- | --------- | ------------------------------------------------------------------------------------------------------------------------------ |
| opening of Section 2, pp. 388–389 | 2.1       | semidiscrete systems (2.1)–(2.2), one-dimensional domain, purpose of the four-model comparison                                 |
| Section 2.1, pp. 389–391          | 2.2       | heat equation and source (2.3)–(2.4), three boundary types, all panels of Figure 2.1, temporal locality                        |
| Section 2.2, pp. 391–393          | 2.3       | equation (2.5), three viscosities, two boundary settings, all eight panels of Figure 2.2, periodic small-viscosity stress test |
| Section 2.3, pp. 393–395          | 2.4       | Burgers equation (2.6), nonlinear deformation and shocks, all eight panels of Figure 2.3, periodic small-viscosity difficulty  |
| Section 2.4, pp. 395–396          | 2.5       | wave equation (2.7), source and initial data, three boundary types, Figure 2.4, long-range bidirectional propagation           |

## Summary

The four models turn the parabolic–hyperbolic distinction into a difference in temporal memory. Heat rapidly forgets fine scales. Advection–diffusion preserves phase as viscosity falls. Burgers' equation generates new fine scales through nonlinearity. The wave equation transports amplitude and phase throughout the time domain. Boundary conditions can shift the observed behavior substantially. Chapter 3 studies methods that retain effectiveness under long-range propagation; Chapter 4 focuses on methods that exploit dissipation.

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 2, pp. 388–396.
