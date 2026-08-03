---
title: "Chapter 2: Model Problems Linking the Parabolic and Hyperbolic World"
description: A source-aligned reading of the heat, advection–diffusion, Burgers, and wave equations
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-2-model-problems
tags:
  - parallel-in-time
  - PDE
---

> [!note] Reading scope and numbering
> This page follows Section 2 of the paper (printed pp. 388–396). Its main headings retain the source numbering exactly: the Section 2 introduction, 2.1 Heat equation, 2.2 Advection–diffusion equation, 2.3 Burgers' equation, and 2.4 Second-order wave equation. Interpretive material and Python results produced for this site appear in explicitly labeled supplements and do not occupy source section numbers.

> [!info] Source figures
> Figures 2.1–2.4 are extracted directly from the paper. Their graphics, axes, panel labels, and order are unchanged. The text below identifies the setup, visible behavior, and PinT conclusion of every panel. The paper is distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); full attribution appears at the end.

## Source-to-page map

| Source location                     | Page section                     | Equations and figures        | Main question                                                   |
| ----------------------------------- | -------------------------------- | ---------------------------- | --------------------------------------------------------------- |
| Section 2 introduction, pp. 388–389 | Section 2 introduction           | (2.1)–(2.2)                  | Why four PDEs link parabolic and hyperbolic dynamics            |
| Section 2.1, pp. 389–391            | 2.1 Heat equation                | (2.3)–(2.4), Figure 2.1(a–d) | How diffusion and boundaries create temporal locality           |
| Section 2.2, pp. 391–393            | 2.2 Advection–diffusion equation | (2.5), Figure 2.2(a–h)       | How far fine scales travel as diffusion weakens                 |
| Section 2.3, pp. 393–395            | 2.3 Burgers' equation            | (2.6), Figure 2.3(a–h)       | How nonlinearity generates and transports shocks                |
| Section 2.4, pp. 395–396            | 2.4 Second-order wave equation   | (2.7), Figure 2.4(a–d)       | Why multidirectional propagation and reflection preserve memory |

## Section 2 introduction: why these four models?

The paper uses PDEs to place parabolic and hyperbolic dynamics in one sequence. Linear advection–diffusion contains both diffusion and transport, so varying the diffusion parameter moves the problem from a parabolic-dominated regime toward a hyperbolic limit. Burgers' equation adds nonlinearity. The heat and second-order wave equations provide cleaner parabolic and hyperbolic reference problems.

Many PinT methods are formulated and analyzed after spatial semidiscretization. The paper first introduces the linear system

$$
\begin{aligned}
\boldsymbol u'(t)&=A\boldsymbol u(t)+\boldsymbol g(t),
&&t\in(0,T],\\
\boldsymbol u(0)&=\boldsymbol u_0,
\end{aligned}
\tag{2.1}
$$

where $A\in\mathbb R^{N_x\times N_x}$ is the matrix produced by spatial semidiscretization. Domain-decomposition methods are the stated exception because they often operate on continuous space–time subproblems rather than an ODE system.

The nonlinear form is

$$
\begin{aligned}
\boldsymbol u'(t)&=\boldsymbol f(\boldsymbol u(t),t),
&&t\in(0,T],\\
\boldsymbol u(0)&=\boldsymbol u_0.
\end{aligned}
\tag{2.2}
$$

Here $\boldsymbol f:\mathbb R^{N_x}\times\mathbb R\to\mathbb R^{N_x}$ is nonlinear in its first argument. The Burgers-type example in the paper is
$\boldsymbol f(\boldsymbol u(t),t)=A\boldsymbol u(t)+B\boldsymbol u^2(t)+\boldsymbol g(t)$.

Every model is posed on the one-dimensional unit interval $\Omega=(0,1)$. The authors choose one dimension to avoid unnecessary notation and state that the applicability and convergence properties of PinT methods generally do not depend on spatial dimension. Dimension changes the cost of a spatial solve, while the mechanisms studied here remain: decay of past information, survival of fine scales, and the escape or return of information through a boundary.

## 2.1 Heat equation

> [!abstract] Source location
> Section 2.1, printed pp. 389–391; equations (2.3)–(2.4); Figure 2.1(a–d).

The parabolic model is

$$
\partial_tu(x,t)=\partial_{xx}u(x,t)+g(x,t),
\qquad (x,t)\in\Omega\times(0,T], \tag{2.3}
$$

with initial data $u(x,0)=u_0(x)$. The text first introduces homogeneous Dirichlet or Neumann boundaries. The numerical figure also includes periodic boundaries to test whether heat can leave the domain.

Figure 2.1(a–c) uses zero initial data and four localized heating events:

$$
g(x,t)=10\sum_{j=1}^{4}
\exp\!\left(
-\sigma\left[(t-t_j)^2+(x-0.5)^2\right]
\right), \tag{2.4}
$$

where

$$
(t_1,t_2,t_3,t_4)=(0.1,0.6,1.35,1.85),
\qquad \sigma=200.
$$

The source is centered at $x=0.5$. Its four separated time centers produce four bright bands whose later behavior depends on diffusion and the boundary condition.

![Source Figure 2.1: heat-equation solutions for three boundary conditions and oscillatory initial data](assets/papers/time-parallelization/source-figures/figure-2-1.svg)

### Figure 2.1(a): homogeneous Dirichlet boundaries

Panel (a) uses $u(0,t)=u(1,t)=0$, zero initial data, and source (2.4). Each heating event forms a bright localized region that spreads and decays rapidly. Heat can leave through the boundary, so early pulses leave little detailed influence on later intervals.

The paper makes this PinT intuition concrete. The contribution of the fourth source on $t\in(1.7,2.2)$ can be computed without first resolving the complete earlier solution. A living-room analogy follows: a temperature prediction a week or month ahead depends chiefly on whether the heater will be on and the windows closed, while the detailed earlier temperature field matters little.

### Figure 2.1(b): homogeneous Neumann boundaries

Panel (b) uses $u_x(0,t)=u_x(1,t)=0$. Diffusion still removes spatial variation from each pulse, but its spatial mean cannot escape. The solution on $t\in(1.7,2.2)$ therefore retains contributions from the first three heating events.

The paper compares this case with a perfectly insulated room. Its later temperature depends on how often and how long the heater ran because the added heat remains. A real room leaks heat slowly, and the authors note that Robin boundary conditions would model that situation more realistically.

### Figure 2.1(c): periodic boundaries

Panel (c) retains zero initial data and source (2.4) but uses periodic boundaries. Its long-time behavior resembles the Neumann case. Spatial variation is smoothed rapidly, while the constant component accumulates. The solution on $t\in(1.7,2.2)$ still contains information from the first three sources because heat cannot leave the periodic domain.

### Figure 2.1(d): oscillatory data with periodic boundaries

Panel (d) has a different experiment. The source is zero, the boundary remains periodic, and

$$
u_0(x)=\sin^2\!\left(8\pi(1-x)^2\right).
$$

Fine spatial stripes survive only in a narrow region near $t=0$. The plot then becomes nearly uniform. This uniform field has a precise interpretation: periodic heat flow damps every nonconstant Fourier mode rapidly and preserves the spatial mean of the initial condition. The paper further observes that the surviving constant is approximately the same as the constant produced by the first two sources in panels (b) and (c).

Panel (d) visualizes frequency selection. Diffusion quickly removes the exact oscillatory profile, while the constant mode spans the full horizon. The nearly featureless upper part of the panel records persistent low-frequency information.

### PinT conclusion of Section 2.1

The paper describes the Dirichlet heat problem as highly local in time and cites Gander, Ohlberger, and Rave (2024). It also compares temporal locality with spatial locality in solvation models from computational chemistry, citing Ciaramella and Gander (2017, 2018a, 2018b).

Neumann and periodic problems remain candidates for PinT if low-frequency components, such as the constant mode, can be communicated accurately over long times. A coarse grid is the mechanism proposed in the paper. Figure 2.1 therefore separates two tasks: diffusion removes high-frequency error locally, while slow components require global temporal communication.

## 2.2 Advection–diffusion equation

> [!abstract] Source location
> Section 2.2, printed pp. 391–393; equation (2.5); Figure 2.2(a–h).

The model on the unit interval is

$$
\partial_tu(x,t)+\partial_xu(x,t)
-\nu\partial_{xx}u(x,t)
=g(x,t),
\qquad (x,t)\in\Omega\times(0,T], \tag{2.5}
$$

with $u(x,0)=u_0(x)$ and $\nu>0$. The advection speed is one, and $\nu$ sets the diffusion strength. The paper compares homogeneous Dirichlet and periodic boundaries. A footnote states that Neumann conditions would add no new qualitative behavior, so they are omitted.

![Source Figure 2.2: eight advection–diffusion solutions with Dirichlet and periodic boundaries](assets/papers/time-parallelization/source-figures/figure-2-2.svg)

The top row, panels (a–d), uses zero Dirichlet boundaries. The bottom row, panels (e–h), uses periodic boundaries. Panels (a–c) and (e–g) use $u_0=0$, source (2.4), and

$$
\nu=1,\qquad 10^{-2},\qquad 5\times10^{-4},
$$

from left to right. Panels (d) and (h) use zero source, the oscillatory data from Figure 2.1(d), and $\nu=5\times10^{-4}$.

### Figures 2.2(a) and (e): $\nu=1$

Diffusion dominates panel (a), so the four heating events look much like the horizontal bands in Figure 2.1(a). Dirichlet boundaries allow them to decay and leave. Panel (e) is also diffusion-dominated, but periodic boundaries retain the constant component. The four inputs build a layered long-time background. The contrast comes chiefly from the boundary treatment of low frequencies.

### Figures 2.2(b) and (f): $\nu=10^{-2}$

The bright bands in panel (b) lean to the right, showing that transport now controls signal location. Each band leaves after reaching $x=1$. Panel (f) returns the signal from $x=1$ to $x=0$, so trajectories from earlier pulses remain visible much longer. Diffusion still broadens each band during propagation.

### Figures 2.2(c) and (g): $\nu=5\times10^{-4}$

Panel (c) contains narrower trajectories that are closer to characteristic translation. Fine scales survive longer, although Dirichlet outflow eventually removes every signal. Panel (g) recirculates the same sharp trajectories across the periodic boundary. Events far in the past still determine detailed features of the current state.

### Figures 2.2(d) and (h): one initial condition, two outcomes

Both panels use $g=0$, $\nu=5\times10^{-4}$, and the oscillatory initial data. Fine stripes move from left to right in panel (d) and have left the Dirichlet domain by approximately $t=1$, after which the solution is nearly zero. In panel (h), the stripes repeatedly wrap around the periodic domain and retain substantial high-frequency content through $t=3$. This pair isolates the effect of outflow from the effect of weak diffusion.

### PinT conclusion of Section 2.2

For Dirichlet boundaries, every component eventually diffuses or leaves, for both large and small $\nu$. The paper therefore states that the solution on $t\in(1.25,2.5)$ can be computed before the complete earlier solution is available.

The conclusion changes with $\nu$ under periodic boundaries. Panel (e) retains only coarse, low-frequency information, so a coarse grid may suffice. Panels (f) and (g) carry progressively finer information over long times as diffusion weakens. Panel (h) shows that fine-scale initial data also persists. A later interval can no longer be determined independently of the earlier state. The difficulty is strongest as $\nu\to0$ in the hyperbolic limit.

The paper emphasizes the test configuration. Periodicity and small diffusion must appear together to expose the long-range transport challenge. A Dirichlet outflow test hides it.

> [!note] Interpretive supplement: phase and damping
> For a continuous Fourier mode $e^{ikx}$, propagation contains the phase factor $e^{-ikt}$ and damping $e^{-\nu k^2t}$. Small $\nu$ places more modes near pure advection. A stable coarse propagator may still fail when a phase-speed error places a correction at the wrong spatial location. This frequency-domain formula is supplied for interpretation; it is not written in Section 2.2 of the paper.

## 2.3 Burgers' equation

> [!abstract] Source location
> Section 2.3, printed pp. 393–395; equation (2.6); Figure 2.3(a–h).

To compare PinT methods in a nonlinear setting, the paper uses

$$
\begin{aligned}
\partial_tu(x,t)-\nu\partial_{xx}u(x,t)
+\frac12\partial_x\!\left(u^2(x,t)\right)
&=g(x,t),
&& (x,t)\in\Omega\times(0,T],\\
u(x,0)&=u_0(x),
&&x\in\Omega,
\end{aligned}
\tag{2.6}
$$

where $\nu>0$. The boundary conditions, source, initial data, and three diffusion values match Figure 2.2, enabling a direct comparison between linear and nonlinear transport.

![Source Figure 2.3: eight Burgers solutions with Dirichlet and periodic boundaries](assets/papers/time-parallelization/source-figures/figure-2-3.svg)

### Figures 2.3(a) and (e): $\nu=1$

Diffusion dominates both panels. They resemble the Dirichlet and periodic heat cases, respectively. Four inputs remain separated and decay in panel (a), while panel (e) retains an accumulated low-frequency background. Nonlinear transport has not yet produced prominent steep fronts.

### Figures 2.3(b) and (f): $\nu=10^{-2}$

The pulses in panel (b) move right and deform asymmetrically. Regions of larger amplitude travel faster, so the leading and trailing edges are no longer nearly parallel as in linear advection–diffusion. Panel (f) recirculates these deformed structures, leaving persistent backgrounds and slanted trajectories from earlier inputs.

### Figures 2.3(c) and (g): $\nu=5\times10^{-4}$

With weaker diffusion, a smooth source produces steep edges: viscous shocks containing high spatial frequencies. These fronts eventually leave the Dirichlet domain in panel (c). Panel (g) returns them through the periodic boundary, carrying fine structure far in space and time. The source text points back to Figure 2.2(b,c) as the linear transport comparison. Shape change and shock generation are the additional phenomena in Figure 2.3.

### Figures 2.3(d) and (h): oscillatory data forms sharper fronts

Both panels use $g=0$, $\nu=5\times10^{-4}$, and oscillatory initial data. The initial profile already contains high-frequency components. Nonlinear evolution compresses gradients further and forms sharper shocks. Dirichlet outflow gradually clears panel (d), leaving a weak tail at late times. Periodicity recirculates multiple deformed fronts across the entire horizon in panel (h).

### PinT conclusion of Section 2.3

Every component eventually diffuses or leaves in the Dirichlet case. As in Section 2.2, the paper states that the later interval $t\in(1.25,2.5)$ may be computed before the full earlier solution is known.

Small diffusion and periodicity amplify the difficulty found in linear advection–diffusion. Nonlinearity continually generates high-frequency shocks. A successful PinT algorithm must communicate their locations and shapes over long distances in space and time, which is difficult on a coarse grid. Panel (h) shows that initial data alone can sustain this fine structure.

As $\nu\to0$, the problem approaches a hyperbolic limit with natural shocks. Later states retain dependence on the full frequency content of earlier states, removing the precomputation opportunity. The paper again identifies periodic boundaries and small diffusion as the key stress-test conditions. The wave equation in the next section exhibits long-range detailed propagation without those two additional conditions.

## 2.4 Second-order wave equation

> [!abstract] Source location
> Section 2.4, printed pp. 395–396; equation (2.7); Figure 2.4(a–d). Figure 2.4 floats ahead of the Section 2.4 heading in the paper layout, but its experiment and discussion belong to Section 2.4.

The hyperbolic model is

$$
\begin{aligned}
\partial_{tt}u(x,t)&=c^2\partial_{xx}u(x,t)+g(x,t),
&& (x,t)\in(0,1)\times(0,T],\\
u(x,0)&=u_0(x),
&&x\in(0,1),\\
\partial_tu(x,0)&=0,
&&x\in(0,1),
\end{aligned}
\tag{2.7}
$$

where $c>0$. Figure 2.4 uses $c^2=0.2$.

![Source Figure 2.4: wave-equation solutions for three boundary conditions and oscillatory initial data](assets/papers/time-parallelization/source-figures/figure-2-4.svg)

### Figure 2.4(a): Dirichlet boundaries and localized sources

Panel (a) uses zero initial displacement, zero initial velocity, and source (2.4). Every localized event launches waves in both directions, followed by reflections at the boundaries. Propagation paths from early events still influence later states after several reflections, producing V-shaped and overlapping structures across the horizon.

### Figure 2.4(b): Neumann boundaries and localized sources

Panel (b) uses the same source and zero initial data. Neumann reflection changes how waves combine at the boundary but does not remove long-range propagation. The influence of all four events spreads along bidirectional characteristics and overlaps at later times.

### Figure 2.4(c): periodic boundaries and localized sources

Panel (c) returns waves through the opposite side of the domain. Its large-scale bands resemble the Neumann case, while smaller ripples record each source event. Periodicity reconnects the paths, but detailed propagation remains.

### Figure 2.4(d): oscillatory data propagates all frequencies

Panel (d) uses zero source, periodic boundaries, zero initial velocity, and

$$
u_0(x)=\sin^2\!\left(8\pi(1-x)^2\right).
$$

Multiple frequencies travel in both directions and repeatedly superpose. Dense interference patterns remain throughout $0<t<3$. The paper states that the same detailed dependence on every initial frequency would also occur with Dirichlet or Neumann boundaries. Panel (d) therefore represents the general long-term memory of hyperbolic problems, and the same dependence persists under all three boundary types.

### PinT conclusion of Section 2.4

First-order advection has a preferred direction. Periodic boundaries are needed in the advection–diffusion and Burgers examples to recirculate signals and expose the small-$\nu$ difficulty fully. Waves travel in multiple directions and reflect under Dirichlet and Neumann boundaries, so all three boundary types retain detailed long-time information.

This propagation makes PinT more difficult than in parabolic problems. Phase or wave-speed error in a coarse temporal model is not damped automatically. Effective methods need a mechanism that communicates detailed information over long times. Chapter 3 therefore studies SWR, parallel IDC, ParaExp, and ParaDiag.

## Site interpretation: a unified view of the four models

The diagram and table below synthesize the four source sections.

![The four model equations span temporal locality to long-range memory](assets/diagrams/pint/en/model-memory-spectrum.svg)

| Model               | Dominant mechanism                         | Information retained over long times                                                                                | PinT implication in the paper                                |
| ------------------- | ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| heat                | diffusion                                  | rapid forgetting with Dirichlet boundaries; constants and other low frequencies with Neumann or periodic boundaries | communicate slow modes accurately on coarse levels           |
| advection–diffusion | directed transport and diffusion           | phase and fine scales under periodic boundaries and small $\nu$                                                     | move progressively finer information across long times       |
| Burgers             | nonlinear transport, shocks, and diffusion | high-frequency fronts that persist and are regenerated under periodic boundaries and small $\nu$                    | track shock location and shape in a coarse representation    |
| wave                | bidirectional propagation and reflection   | amplitude, phase, and propagation paths under every boundary type                                                   | use a temporal mechanism designed for long-range propagation |

Boundary conditions decide whether information can leave. Diffusion sets the lifetime of fine scales. Nonlinearity determines whether high frequencies are regenerated. These three questions test whether an experiment exposes the actual challenge faced by a PinT method.

## Site numerical supplement: three recomputed solutions

The following Python experiments contrast diffusion plus outflow with persistent propagation. They use independent initial data, grids, and source-free settings. They do not reproduce the source-driven experiments in Figures 2.2–2.4. Their numerical values are reported separately from the paper figures.

### Advection–diffusion solution

The experiment uses homogeneous Dirichlet boundaries, $\Delta t=\Delta x=10^{-3}$, $T=3$, $\nu=5\times10^{-4}$, unit advection speed, and initial data $\sin^2(8\pi(1-x)^2)$. The final $L^\infty$ norm is $7.022\times10^{-99}$.

![Recomputed space–time solution of the advection–diffusion equation](assets/pint/model-advection-diffusion.svg)

The near-zero final state results from the present combination of diffusion and Dirichlet outflow. A periodic problem retains recirculating signals.

### Viscous Burgers solution

The experiment uses $\Delta t=\Delta x=1/400$, $T=3$, $\nu=5\times10^{-4}$, and homogeneous Dirichlet boundaries. The maximum over space and time is $1.045940$, and the final $L^\infty$ norm is $0.325871$.

![Recomputed space–time solution of the viscous Burgers equation](assets/pint/model-burgers.svg)

### Wave solution

The source-free experiment uses the trapezoidal rule with $\Delta t=\Delta x=1/400$, $T=3$, and $c^2=0.2$. It uses the same oscillatory initial displacement as the previous experiments and zero initial velocity. The final displacement has $L^\infty$ norm $0.948217$.

![Recomputed space–time solution of the wave equation](assets/pint/model-wave.svg)

### Numerical summary

| Supplemental experiment | Grid and horizon                   |                      Reported metric |
| ----------------------- | ---------------------------------- | -----------------------------------: |
| advection–diffusion     | $\Delta x=\Delta t=10^{-3}$, $T=3$ | final $L^\infty=7.022\times10^{-99}$ |
| viscous Burgers         | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.325871$ |
| wave equation           | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.948217$ |

These quantities characterize the discrete dynamics. They do not measure PinT convergence or hardware speedup.

## Completeness audit

| Item checked                                                                                                                                                    | Result             |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------ |
| The two semidiscrete systems, matrix dimensions, nonlinear example, domain-decomposition exception, and one-dimensional rationale in the Section 2 introduction | covered            |
| Figure 2.1(a–d): boundaries, source/data, specific time interval, room analogy, Robin comment, and low-frequency conclusion                                     | covered            |
| Figure 2.2(a–h): two boundary types, three $\nu$ values, two data configurations, precomputable interval, and Neumann footnote                                  | covered            |
| Figure 2.3(a–h): linear comparison, deformation, shocks, high-frequency generation, precomputable interval, and limiting difficulty                             | covered            |
| Figure 2.4(a–d): three boundary types, bidirectional propagation, reflections, all frequencies, and cross-boundary conclusion                                   | covered            |
| Notation, domains, initial conditions, and numbering of equations (2.1)–(2.7)                                                                                   | checked            |
| Source claims, site frequency interpretation, synthesis diagram, and Python experiments                                                                         | labeled separately |

## Summary

Section 2 lengthens temporal memory by varying boundaries, diffusion, and nonlinearity. The Dirichlet heat equation forgets old information rapidly. Neumann and periodic heat problems preserve constants and other low-frequency components. Directed transport extends memory; small diffusion and periodic boundaries carry fine scales across long times. Burgers nonlinearity continues to generate shocks and high-frequency fronts. The wave equation propagates and reflects detailed information under every common boundary type.

This progression motivates the organization of the paper. Chapter 3 treats methods capable of long-range propagation. Chapter 4 studies methods that primarily exploit parabolic dissipation and temporal locality.

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 2, pp. 388–396.
