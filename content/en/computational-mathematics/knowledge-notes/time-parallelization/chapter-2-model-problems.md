---
title: "Chapter 2: Model Problems Linking the Parabolic and Hyperbolic World"
description: A paragraph-by-paragraph reading of the heat, advection–diffusion, Burgers, and wave equations, explaining how boundary conditions, diffusion strength, and nonlinearity push solutions from temporal locality toward long-range memory, and what this means fundamentally for time parallelization
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-2-model-problems
tags:
  - parallel-in-time
  - PDE
---

> [!note] Reading scope and numbering
> This page follows Section 2 of the paper (printed pp. 388–396). The main headings retain the source numbering exactly: the Section 2 introduction, 2.1 Heat equation, 2.2 Advection–diffusion equation, 2.3 Burgers' equation, and 2.4 Second-order wave equation. Interpretation, comparison diagrams, and Python experiments added for this site are collected under "interpretive supplements," "Insight," or "numerical supplements" and do not occupy the source section numbers.

> [!info] Source figures
> Figures 2.1–2.4 are extracted directly from the paper, with their graphics, axes, panel labels, and order unchanged. The text explains, item by item, the parameters of each panel, the visible phenomena, and the PinT judgment the paper draws from it. The paper is distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); full attribution appears at the end.

## Source-to-page map

| Source location                     | Page section                     | Equations and figures        | Main question                                                   |
| ----------------------------------- | -------------------------------- | ---------------------------- | --------------------------------------------------------------- |
| Section 2 introduction, pp. 388–389 | Section 2 introduction           | (2.1)–(2.2)                  | Why four PDEs link parabolic and hyperbolic dynamics            |
| Section 2.1, pp. 389–391            | 2.1 Heat equation                | (2.3)–(2.4), Figure 2.1(a–d) | Diffusion, boundary conditions, and temporal locality           |
| Section 2.2, pp. 391–393            | 2.2 Advection–diffusion equation | (2.5), Figure 2.2(a–h)       | How far fine scales travel once viscosity drops                 |
| Section 2.3, pp. 393–395            | 2.3 Burgers' equation            | (2.6), Figure 2.3(a–h)       | How nonlinearity generates and transports high-frequency shocks |
| Section 2.4, pp. 395–396            | 2.4 Second-order wave equation   | (2.7), Figure 2.4(a–d)       | Why multidirectional propagation and reflection preserve memory |

## Section 2 introduction: why these four models?

The paper uses a set of PDEs to place parabolic and hyperbolic dynamics in one chain of observation, so as to show that parallel-in-time (PinT) methods must behave fundamentally differently on the two classes of problems. This requires a set of models that both transition continuously and provide pure reference points at each end. The four equations the paper selects form exactly such a chain: the heat equation gives a fairly pure parabolic reference, the second-order wave equation gives a fairly pure hyperbolic reference, the linear advection–diffusion equation slides continuously between the two ends through a single diffusion parameter, and Burgers' equation then adds nonlinearity on top of the advection–diffusion skeleton.

The linear advection–diffusion equation is the hinge of this chain, because it contains both a parabolic component and a hyperbolic component in the same equation: the second-order diffusion term $-\nu\partial_{xx}u$ is parabolic, while the first-order advection term $\partial_xu$ is hyperbolic. Lowering the diffusion parameter $\nu$ weakens the parabolic component and highlights the hyperbolic one, so the problem moves smoothly from diffusion-dominated toward advection- (hyperbolic-) dominated. It is precisely with this single adjustable parameter that the paper displays, continuously, the key PinT question of whether the solution retains early fine detail over long times.

> [!tip] Insight: why advection–diffusion can serve as a bridge
> The two components have opposite physical characters, and that is exactly why the equation can connect the two worlds. The diffusion term is isotropic and directionless; it smooths energy across all wavenumbers and decays them as $e^{-\nu k^2t}$, with higher wavenumbers decaying faster, so a parabolic solution rapidly "forgets" high-frequency detail. The advection term is first order and directional; it translates a signal at finite speed without decay, so the hyperbolic component carries detail intact to distant locations. The parameter $\nu$ measures their relative strength: when $\nu$ is large, decay overwhelms translation and the solution behaves like the heat equation; as $\nu\to0$ the equation degenerates into pure advection $\partial_tu+\partial_xu=0$, becoming hyperbolic, and detail can travel arbitrarily far. With this in mind, one can anticipate, in every later figure, how far fine scales are able to travel.

Many time-parallel methods take the semidiscrete ODE as their object of description and analysis. The paper first writes down the linear system

$$
\begin{aligned}
\boldsymbol u'(t)&=A\boldsymbol u(t)+\boldsymbol g(t),
&&t\in(0,T],\\
\boldsymbol u(0)&=\boldsymbol u_0,
\end{aligned}
\tag{2.1}
$$

where $A\in\mathbb R^{N_x\times N_x}$ comes from the spatial semidiscretization of the PDE. The reason for first collapsing the PDE into such a family of ordinary differential equations in time is that most time-parallel algorithms (for example the ParaExp, ParaDiag, and parallel IDC methods discussed later) are designed and analyzed for ODE systems: time is the direction to be parallelized, while space has already been discretized into a finite-dimensional vector $\boldsymbol u(t)\in\mathbb R^{N_x}$. The paper also points out one exception—domain-decomposition methods usually operate directly on continuous space–time subproblems and do not pass through this ODE formulation, so the way they are analyzed later is also different.

The nonlinear problem is written as

$$
\begin{aligned}
\boldsymbol u'(t)&=\boldsymbol f(\boldsymbol u(t),t),
&&t\in(0,T],\\
\boldsymbol u(0)&=\boldsymbol u_0.
\end{aligned}
\tag{2.2}
$$

Here $\boldsymbol f:\mathbb R^{N_x}\times\mathbb R\to\mathbb R^{N_x}$ is nonlinear in its first argument. The Burgers-type example given in the paper is
$\boldsymbol f(\boldsymbol u(t),t)=A\boldsymbol u(t)+B\boldsymbol u^2(t)+\boldsymbol g(t)$,
where the linear part $A\boldsymbol u$ comes from the semidiscretization of the diffusion term and the quadratic term $B\boldsymbol u^2$ comes from the nonlinear advection $\tfrac12\partial_x(u^2)$. The point of writing it this way is that the nonlinearity is captured in this single quadratic term, making it easy to compare term by term with the linear system (2.1) and to see clearly "what extra appears once nonlinearity is added."

All models are posed on the one-dimensional unit interval $\Omega=(0,1)$. The authors adopt a one-dimensional setting to reduce notation and explicitly state that the applicability and convergence properties of PinT methods generally do not depend on the spatial dimension. The spatial dimension significantly changes the cost of a single spatial solve, but the three mechanisms this chapter focuses on are independent of dimension and remain: the decay of old information over time, the preservation of fine scales, and the escape or return of information through the boundary. It is these three mechanisms that decide "whether the late solution must depend on the early solution," and hence how many parallelizable degrees of freedom remain in the time direction.

## 2.1 Heat equation

> [!abstract] Source location
> Section 2.1, printed pp. 389–391; equations (2.3)–(2.4); Figure 2.1(a–d).

The parabolic model in the paper is

$$
\partial_tu(x,t)=\partial_{xx}u(x,t)+g(x,t),
\qquad (x,t)\in\Omega\times(0,T], \tag{2.3}
$$

with initial data $u(x,0)=u_0(x)$. The text first introduces homogeneous Dirichlet or Neumann boundaries, while the numerical figure additionally includes periodic boundaries, to contrast how the single question of "whether heat can leave the domain" differs across the three boundary types. The heat equation is a pure parabolic reference because it has only a diffusion term and no transport term: every spatial frequency decays monotonically, and the equation has no mechanism to carry a signal intact.

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

The heat source is concentrated at $x=0.5$ in space and lands at four well-separated centers in time. The large value $\sigma=200$ makes each Gaussian bump narrow in space and time, so each heating event forms a clear, isolated bright band in the figure whose later shape is determined jointly by diffusion and the boundary. The advantage of this design is that each pulse is an "event" that can be tracked on its own, making it easy to judge by eye whether it is quickly forgotten, retained over the long term, or carried elsewhere.

![Source Figure 2.1: heat-equation solutions for three boundary conditions and oscillatory initial data](assets/papers/time-parallelization/source-figures/figure-2-1.svg)

### Figure 2.1(a): homogeneous Dirichlet boundaries

Panel (a) uses $u(0,t)=u(1,t)=0$, zero initial data, and source (2.4). The four bright bands are well separated, and the temperature peak formed by each heating event quickly spreads and decays. There are two reasons: diffusion rapidly flattens the local peak, and the Dirichlet boundary conducts heat out of the domain, so together they leave the early pulses with almost no trace in later intervals.

The paper draws from this a key parallel intuition: when computing the solution produced by the fourth heat source on $t\in(1.7,2.2)$, one can first ignore the complete solution at earlier times. The authors stress that this is a canonical example of being able to "compute a later segment first, despite causality"—causality says the future is determined by the past, but once the influence of the past is rapidly erased by diffusion and outflow, the late state depends so weakly on early detail that the late interval can be solved almost independently. The authors use a winter room as an analogy: predicting the room temperature a week or a month ahead depends chiefly on whether the heater is on and the windows are closed at that time, while the fine temperature distribution at some earlier moment matters almost not at all.

### Figure 2.1(b): homogeneous Neumann boundaries

Panel (b) uses $u_x(0,t)=u_x(1,t)=0$. Diffusion still flattens the local structure formed by each heating event, but the zero-flux boundary does not allow heat to flow out through the endpoints, so the total heat within the domain (the spatial mean) cannot dissipate. Therefore on $t\in(1.7,2.2)$ the solution still contains the contributions accumulated from the first three heat sources. The paper points out that this is exactly an expression of causality: here the late state does remember what happened early.

The paper compares this with a perfectly insulated room: one must know how often and how long the heater ran in the past, because every unit of input heat stays in the room forever. The authors then add a realistic correction—a real room always dissipates heat slowly and can never be perfectly insulated, and this "slow leakage" situation is better modeled by Robin boundary conditions, which lie between Dirichlet (free dissipation) and Neumann (no dissipation at all).

### Figure 2.1(c): periodic boundaries

Panel (c) still uses zero initial data and source (2.4), with the boundary changed to a periodic condition. Its long-time behavior resembles the Neumann case: local spatial variation is quickly smoothed by diffusion, while the constant component of the heat keeps accumulating; the state on $t\in(1.7,2.2)$ is still influenced by the first three sources. The only difference is the mechanism—the periodic boundary lets heat that leaves one end return from the other, so the total heat likewise cannot leave the domain. In other words, Neumann conserves total heat by "sealing the endpoints," while periodicity conserves it by "joining the two ends," and the two have the same net effect on the low-frequency component.

### Figure 2.1(d): oscillatory data with periodic boundaries

Panel (d) uses a different experimental setup than the first three figures. Here the external source is turned off, i.e. $g=0$, the periodic boundary is kept, and

$$
u_0(x)=\sin^2\!\left(8\pi(1-x)^2\right)
$$

is used. This initial condition deliberately contains many fine stripes of high spatial frequency. In the figure these very fine stripes survive only near $t=0$, after which the field quickly becomes an almost uniform block of color. This block has a definite meaning: the periodic heat equation damps all nonconstant Fourier modes within a very short time and preserves only the spatial mean of the initial condition. The source further notes that this surviving constant is roughly comparable to the constant accumulated by the first two heat sources in panels (b) and (c)—low-frequency components from different sources end up at the same order of magnitude.

Panel (d) therefore clearly displays a process of "frequency selection": the exact oscillatory shape (high frequency) is completely erased within a very short time, while the constant mode (the lowest frequency) survives across the entire time interval. The latter part of the figure looks structureless, but it in fact records the result that "only low-frequency information can be retained over the long term."

> [!tip] Insight: why Dirichlet forgets while Neumann/periodic keeps only the constant
> The difference between the three figures can be explained in one sentence in spectral language. Expand the solution in the eigenfunctions of the Laplacian; each mode evolves as $e^{\lambda t}$, with $\lambda$ the corresponding eigenvalue. Under homogeneous Dirichlet conditions all eigenvalues are strictly negative (in one dimension about $-(k\pi)^2$), so **every** mode decays exponentially and the solution "completely forgets" the initial data and early sources—this is the origin of temporal locality. Neumann and periodic boundaries, by contrast, admit a constant mode with a zero eigenvalue ($\lambda=0$) that neither decays nor grows, so uniquely this lowest-frequency component survives over the long term while the remaining high-frequency modes still vanish rapidly as $e^{-\nu k^2t}$. This explains why the later text says PinT remains feasible so long as one "effectively carries low-frequency components such as the constant far into time"—the degrees of freedom that truly require global-in-time communication are very few. This spectral view is added by this site to tie the three panels together; the source states it through a physical analogy and does not write out the eigenvalues.

### PinT conclusion of Section 2.1

The paper calls the Dirichlet heat equation a highly temporally local case: despite causality, the solution of the heat equation and more general parabolic problems under Dirichlet boundaries is "completely local in time," so time parallelization is in fact quite easy, and the authors cite the related results of Gander, Ohlberger, and Rave (2024) for this. The authors also compare this "temporal locality" with the "spatial locality" in solvation models from computational chemistry, where distant detail has weak influence on the local solution; see Ciaramella and Gander (2017, 2018a, 2018b).

The Neumann and periodic problems still leave room for time parallelization, but the premise has changed: the algorithm must be able to carry low-frequency components such as the constant effectively far into time, and the typical mechanism the paper offers is a coarse grid (a coarse time level). Figure 2.1 therefore divides the PinT task into two kinds—high-frequency error is eliminated automatically and rapidly by diffusion and needs no long-range communication; the low-frequency component does not decay and must rely on communication across the full time scale to be computed correctly. Grasping this dichotomy explains why later algorithms generally adopt a structure of "fine levels for local advancement + coarse levels to carry low frequencies."

## 2.2 Advection–diffusion equation

> [!abstract] Source location
> Section 2.2, printed pp. 391–393; equation (2.5); Figure 2.2(a–h).

The paper considers, on the unit interval,

$$
\partial_tu(x,t)+\partial_xu(x,t)
-\nu\partial_{xx}u(x,t)
=g(x,t),
\qquad (x,t)\in\Omega\times(0,T], \tag{2.5}
$$

with $u(x,0)=u_0(x)$ and $\nu>0$. The transport speed is $1$ (the coefficient of the advection term $\partial_xu$), and $\nu$ controls the diffusion strength. Compared with the heat equation, the additional first-order term $\partial_xu$ brings direction: a signal is translated from left to right at speed $1$, and this is precisely the hyperbolic component. The paper compares homogeneous Dirichlet and periodic boundaries, and states in a footnote that Neumann conditions bring no new qualitative phenomenon, so no separate set of figures is included—its low-frequency retention behavior is similar to the periodic case and provides no new information.

![Source Figure 2.2: eight advection–diffusion solutions with Dirichlet and periodic boundaries](assets/papers/time-parallelization/source-figures/figure-2-2.svg)

The top row of Figure 2.2, panels (a–d), uses zero Dirichlet boundaries; the bottom row, panels (e–h), uses periodic boundaries. Panels (a–c) and (e–g) take $u_0=0$ and source (2.4), with the diffusion parameter decreasing from left to right as

$$
\nu=1,\qquad 10^{-2},\qquad 5\times10^{-4}.
$$

Panels (d) and (h) turn off the external source and use the same oscillatory initial data as Figure 2.1(d), with $\nu=5\times10^{-4}$. The purpose of this arrangement is to cross-compare the same source/data across two dimensions—"outflow boundary vs. return boundary" and "large diffusion vs. small diffusion"—so as to separate the roles of the boundary and of diffusion.

### Figures 2.2(a) and (e): $\nu=1$

In panel (a) diffusion dominates and the first-order transport is almost drowned out by strong diffusion; the four heating events appear as nearly horizontal bright bands close to those in Figure 2.1(a), and the Dirichlet boundary lets them gradually decay and leave the domain. Panel (e) is likewise diffusion-dominated, but the periodic boundary preserves the constant component, so the four inputs form a long-time background that accumulates and layers successively. The difference between the two comes almost entirely from how the boundary treats low-frequency information, not from advection—at large $\nu$, advection–diffusion is qualitatively a heat equation with a slight drift.

### Figures 2.2(b) and (f): $\nu=10^{-2}$

The bright bands in panel (b) clearly lean to the right, showing that transport begins to dominate signal location: each band translates along the characteristic $x=x_0+t$ and, upon reaching $x=1$, crosses the boundary and leaves the domain. Panel (f) allows the signal to return from $x=1$ to $x=0$ (periodic recirculation), so the slanted trajectories left by early pulses persist across a longer time and reappear repeatedly. Diffusion is still visible here—the stripes gradually broaden and fade during propagation, indicating that high frequencies are slowly lost—but far less immediately than the flattening at $\nu=1$.

### Figures 2.2(c) and (g): $\nu=5\times10^{-4}$

In panel (c) diffusion is further weakened, the trajectories are narrower and closer to pure translation along characteristics, and fine scales are preserved longer; but Dirichlet outflow still eventually clears out each signal. Panel (g) lets the same sharp trajectories circle the periodic boundary many times, so information produced far earlier than the current moment can still determine the fine structure of the current solution—this is precisely the direct picture of "the late strongly depends on the early." Comparing (c) and (g) shows that small diffusion is responsible for "preserving detail," while the boundary type is responsible for "how far the detail travels."

### Figures 2.2(d) and (h): two fates of the same oscillatory data

Both figures take $g=0$, $\nu=5\times10^{-4}$, and the oscillatory initial data, with the only difference being the boundary. The fine stripes in panel (d) translate from left to right and all cross $x=1$ to leave the Dirichlet domain by roughly $t=1$; thereafter the solution is nearly zero, because no signal is fed in and none returns. In panel (h) the same stripes keep circling back and still retain substantial high-frequency information at $t=3$. This pair of panels cleanly separates two roles: the outflow boundary decides "whether detail can be drained," while small diffusion decides "how complete the detail is before being drained."

### PinT conclusion of Section 2.2

For the Dirichlet problem, regardless of the size of $\nu$, all components are eventually flattened by diffusion or leave through outflow. The paper therefore points out that one can compute the later interval $t\in(1.25,2.5)$ first, before the exact early solution is available—the late state's dependence on early detail vanishes naturally over time, and temporal locality still holds.

The conclusion for the periodic problem changes significantly with $\nu$. Panel (e) (large diffusion) retains only low-frequency components such as the constant over the long term, so a coarse grid may still be effective, consistent with the Neumann/periodic heat equation. Panels (f) and (g) retain progressively finer information as diffusion weakens, and a PinT algorithm must have a mechanism to carry these fine scales effectively far into time; panel (h) shows more directly that high-frequency initial data can travel very far. Here a later interval cannot be determined in advance independently of the early state, and the difficulty is most apparent in the hyperbolic limit $\nu\to0$. The authors stress that this is fundamentally different from the heat equation, and that it **appears only when periodic boundaries are combined with small diffusion**.

The paper specifically warns about the test setup: periodic boundaries and small diffusion must appear together to fully expose the challenge that transport-dominated problems pose for PinT. Using only Dirichlet outflow examples masks long-range propagation and thus overestimates an algorithm's performance on hyperbolic-type problems—this is an easy trap to fall into when evaluating PinT methods.

> [!note] Interpretive supplement: phase and damping
> For a continuous Fourier mode $e^{ikx}$, propagation contains both a phase factor $e^{-ikt}$ and diffusive damping $e^{-\nu k^2t}$. Small $\nu$ places more modes close to pure transport (weak damping, strong phase translation). Even a numerically stable coarse propagator will, if its phase speed is biased, place correction information at the wrong spatial location and thus destroy convergence. This frequency explanation is intended to interpret the source figures; Section 2.2 of the paper does not write out this formula.

## 2.3 Burgers' equation

> [!abstract] Source location
> Section 2.3, printed pp. 393–395; equation (2.6); Figure 2.3(a–h).

To compare how different PinT methods perform on nonlinear problems, the paper uses

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

where $\nu>0$. The key change here is to replace the linear transport term $\partial_xu$ of advection–diffusion with the nonlinear flux $\tfrac12\partial_x(u^2)$, whose equivalent propagation speed becomes $u$ itself, i.e. the local amplitude of the solution. The boundary conditions, source, initial data, and three diffusion values correspond one-to-one with Figure 2.2, and this deliberately parallel design lets linear and nonlinear transport be compared panel by panel, highlighting exactly what nonlinearity adds.

![Source Figure 2.3: eight Burgers solutions with Dirichlet and periodic boundaries](assets/papers/time-parallelization/source-figures/figure-2-3.svg)

### Figures 2.3(a) and (e): $\nu=1$

Diffusion dominates and nonlinear transport has not yet formed a prominent sharp front, so the two figures resemble the Dirichlet and periodic heat-equation forms, respectively. The four inputs in panel (a) remain separated and decay, while panel (e) retains an accumulated low-frequency background. Here strong diffusion rapidly flattens amplitude differences, $u$ varies gently, and nonlinearity plays almost no role.

### Figures 2.3(b) and (f): $\nu=10^{-2}$

The pulses in panel (b) propagate to the right and also deform asymmetrically: regions of higher value move forward faster (because the speed is $u$), so the leading edge of a bright band steepens while the trailing edge stretches, and the front and back edges are no longer nearly parallel as in linear advection–diffusion. Panel (f) repeatedly transports these already deformed structures through the periodic boundary, and the early inputs form persistent slanted trajectories and backgrounds. Here, for the first time, one can see by eye the "shape change" brought by nonlinearity, not merely a change of location.

### Figures 2.3(c) and (g): $\nu=5\times10^{-4}$

Once diffusion is further weakened, nonlinearity keeps amplifying amplitude differences: even if the source itself is smooth, the solution spontaneously forms very steep edges, i.e. viscous shocks containing high spatial frequencies. In panel (c) these steep fronts eventually flow out of the Dirichlet domain and are cleared. Panel (g) lets the shocks circle back through the domain, and fine scales travel far in space and time. In this passage the paper points back to Figure 2.2(b,c) as the linear-transport reference, making clear that the key new phenomenon Figure 2.3 adds relative to advection–diffusion is **shape change and shock generation**, not merely longer propagation.

> [!tip] Insight: why a smooth source can create high-frequency shocks
> In linear advection–diffusion, once a frequency has been weakened by diffusion it is never regenerated, and high frequencies only monotonically drain away; but the nonlinear flux $\tfrac12\partial_x(u^2)$ of Burgers makes the propagation speed vary with $u$, so crests catch up with troughs, gradients are continually compressed, and **new high frequencies are generated continually**. When $\nu$ is very small, diffusion cannot flatten these gradients in time, and compression and dissipation balance at a very steep front, forming a viscous shock of width about $\mathcal O(\nu)$. This explains why, at small $\nu$, ever-sharper structures appear in the solution even when the initial data and source are smooth—and it is especially thorny for PinT, because a coarse grid simply cannot resolve, let alone carry, such $\mathcal O(\nu)$-scale fronts. This mechanistic explanation is added by this site; the source summarizes the phenomenon with the phrase "shock waves."

### Figures 2.3(d) and (h): oscillatory data forms sharper fronts

Both figures take $g=0$, $\nu=5\times10^{-4}$, and the oscillatory initial data. The initial data already contains high-frequency components, and nonlinear evolution further compresses gradients and forms shocks sharper than in the source-driven case. Panel (d) is gradually cleared under Dirichlet outflow, leaving only a very weak tail at late times; panel (h) recirculates continually through the periodic boundary, and multiple deformed fronts span the entire time interval. This is consistent with the (d)/(h) comparison in advection–diffusion, except that nonlinearity sharpens the fronts further and enlarges the difficulty.

### PinT conclusion of Section 2.3

In the Dirichlet case, all components of the solution are eventually cleared by diffusion or outflow. As in Section 2.2, the paper points out that one can compute the later interval $t\in(1.25,2.5)$ first, without waiting for the exact early solution.

Periodic boundaries combined with small diffusion further amplify the difficulty found in linear advection–diffusion: nonlinearity continually generates high-frequency shocks, and a successful PinT algorithm must carry not only these fronts but also their **locations** effectively far into space and time, while a coarse grid can barely represent, let alone transport, sharp shocks. Panel (h) shows that even without an external source, the fine structure produced by the initial data alone is retained over the long term and sharpened by nonlinearity.

As $\nu\to0$, the problem approaches a hyperbolic limit with natural shocks, where all frequency components travel far, the later solution depends on the full frequency content of the early solution, the precomputation opportunity vanishes, and PinT is harder than in the linear case. The paper again reminds the reader that periodic boundaries and small diffusion are the key test conditions for exposing this kind of difficulty; and it previews that the wave equation in the next section needs neither of these two additional conditions and exhibits long-term fine propagation under any boundary—because the advection term of advection–diffusion/Burgers is first order and directional and needs periodic recirculation to reveal long-range memory, whereas the wave equation propagates bidirectionally by itself.

## 2.4 Second-order wave equation

> [!abstract] Source location
> Section 2.4, printed pp. 395–396; equation (2.7); Figure 2.4(a–d). In the paper layout, Figure 2.4 floats ahead of the Section 2.4 heading, but the figure's content and discussion still belong to Section 2.4.

The hyperbolic model in the paper is

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

where $c>0$, and Figure 2.4 uses $c^2=0.2$. Note that (2.7) is second order in time and therefore needs two initial conditions: the initial displacement $u_0$ and the initial velocity $\partial_tu(x,0)=0$ (zero initial velocity). The most essential difference from the previous parabolic/semiparabolic models is that the second-order wave equation has no diffusion term at all—every frequency does not decay but only oscillates, and signals propagate simultaneously in two directions (left- and right-going waves), reflecting when they meet a boundary. This is why it serves as a pure hyperbolic reference.

![Source Figure 2.4: wave-equation solutions for three boundary conditions and oscillatory initial data](assets/papers/time-parallelization/source-figures/figure-2-4.svg)

### Figure 2.4(a): Dirichlet boundaries and localized sources

Panel (a) takes zero initial displacement, zero initial velocity, and source (2.4). Each localized excitation emits waves in both the left and right directions at once, and the waves are reflected back into the domain after reaching the boundary. Since there is no diffusion to weaken them, the propagation paths left by the early sources keep influencing later states even after many reflections, so the figure shows V-shaped trajectories and overlapping interference structures that span long times. The paper concludes from this that the solution of the wave equation depends, for a very long time and in a complex and detailed way, on the source at every space–time location.

### Figure 2.4(b): Neumann boundaries and localized sources

Panel (b) uses the same source and zero initial data, changing only the boundary to Neumann. The phase combination of the reflections changes with the boundary type (Dirichlet reflection flips the sign, Neumann reflection does not), but this does not remove long-range propagation: the influence of the four excitations still spreads along bidirectional characteristics and superposes repeatedly at later times. In other words, changing the boundary changes only the details of the interference pattern, not the fundamental property that "information does not dissipate over the long term."

### Figure 2.4(c): periodic boundaries and localized sources

Panel (c) lets waves return from one end to the other. Its global color bands resemble the Neumann panel, while the local ripples still record, one by one, the bidirectional propagation and reflection of each excitation. Periodicity changes only how the paths connect (joining the two ends rather than reflecting at the endpoints); the core phenomenon that "detail is preserved over the long term" is unchanged. The three boundaries give three different interference patterns yet yield the same PinT conclusion.

### Figure 2.4(d): oscillatory data and long-term propagation of all frequencies

Panel (d) turns off the external source and uses periodic boundaries, zero initial velocity, and

$$
u_0(x)=\sin^2\!\left(8\pi(1-x)^2\right).
$$

The multiple frequencies in the initial data travel in both directions and superpose repeatedly, and dense interference textures are retained throughout the interval $0<t<3$, with no frequency weakened. The source states explicitly that this detailed dependence on **all** initial frequencies also exists under Dirichlet and Neumann boundaries and is not unique to the periodic boundary. Panel (d) therefore exhibits the general long-term memory of hyperbolic problems—all three boundary types preserve this dependence on the full spectrum.

> [!tip] Insight: first-order transport needs "recirculation," while second-order waves are naturally "bidirectional"
> The division of labor between this section and the previous two can be understood uniformly through the dispersion relation. The advection term of advection/Burgers is first order, and its sign fixes the propagation direction (rightward only), so on an open interval the signal will sooner or later flow out and must be sent back by a periodic boundary before long-range memory can be seen at small $\nu$. The time derivative of the wave equation is second order; substituting a plane wave $e^{i(kx-\omega t)}$ gives $\omega=\pm ck$, whose two signs correspond to left-going and right-going families of waves, with $|\text{amplitude}|\equiv1$ containing no decay. So under any boundary, signals propagate bidirectionally, reflect, and superpose losslessly, and long-range fine memory is intrinsic and does not need the pair of additional conditions "periodicity + small diffusion" to be artificially triggered. This is the mechanism-level explanation for what the paper calls the fact that "the difficulty of hyperbolic problems already appears under any boundary"; this derivation is added by this site, while the source states it through the phenomenon and directionality.

### PinT conclusion of Section 2.4

Here the paper contracts the logic of the three sections into a single comparison: first-order advection is directional, and the advection–diffusion and Burgers equations need periodic boundaries to keep signals recirculating, so as to fully reveal the difficulty at small $\nu$; the wave equation, however, propagates in multiple directions and reflects under both Dirichlet and Neumann boundaries, so all three boundary types retain detailed long-term information, without any additional condition.

It is exactly this "all-frequency, multidirectional, decay-free, long-time" propagation that makes time parallelization of hyperbolic problems harder than parabolic ones: once a coarse temporal model has an error in phase or propagation speed, the dynamics do not automatically dissipate the error as parabolic problems do, and the error accumulates along the wave. An effective method must therefore be specifically designed with a mechanism to convey long-range detail. Chapter 3 of the paper accordingly discusses SWR, parallel IDC, ParaExp, and ParaDiag for hyperbolic problems.

## Site interpretation: a unified view of the four models

The diagram and comparison table below are a synthesis by this site, used to connect the conclusions of the four source sections.

![The four models transition from temporal locality to long-range memory](assets/diagrams/pint/en/model-memory-spectrum.svg)

| Model               | Dominant mechanism                         | Information retained over long times                                                                             | PinT implication given in the paper                                |
| ------------------- | ------------------------------------------ | ---------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| heat                | diffusion                                  | rapid forgetting under Dirichlet; constants and other low frequencies retained under Neumann/periodic boundaries | coarse levels must convey slow modes effectively                   |
| advection–diffusion | directed transport and diffusion           | phase and fine scales retained under small $\nu$ and periodic boundaries                                         | must carry progressively finer information far into time           |
| Burgers             | nonlinear transport, shocks, and diffusion | high-frequency fronts retained and regenerated over the long term under small $\nu$ and periodic boundaries      | a coarse representation must track shock location and shape        |
| wave                | bidirectional propagation and reflection   | amplitude, phase, and propagation paths retained under every boundary type                                       | need a time-parallel mechanism designed for long-range propagation |

The whole chain can be summarized with three questions: the boundary conditions decide whether information can leave the domain; the diffusion parameter decides how fast fine scales decay; and nonlinearity decides whether high-frequency structure is regenerated. The heat equation varies only on the first question, advection–diffusion adds the second, Burgers adds the third, and the wave equation, through its intrinsic bidirectional propagation, makes the "hardest setting" of all three questions appear simultaneously under all boundaries. When analyzing each subsequent algorithm, one can use these three questions to judge whether its test setup truly touches the difficulty—being especially wary of an artificially easy "Dirichlet outflow + large diffusion only" setting that masks long-range propagation.

## Site numerical supplement: three recomputed model solutions

The following three Python experiments are used to contrast "diffusion plus outflow" with "persistent propagation." They use their own initial data, grids, and source-free settings, and do not reproduce the source-driven experiments of Figures 2.2–2.4. The numerical values are reported separately from the paper figures.

### Advection–diffusion solution

The experiment uses homogeneous Dirichlet boundaries, $\Delta t=\Delta x=10^{-3}$, $T=3$, $\nu=5\times10^{-4}$, advection speed $1$, and initial data $\sin^2(8\pi(1-x)^2)$. The final $L^\infty$ norm is $7.022\times10^{-99}$.

![Recomputed space–time solution of the advection–diffusion equation](assets/pint/model-advection-diffusion.svg)

This near-zero final state results from the present combination of diffusion and Dirichlet outflow. A periodic problem would retain recirculating signals over the long term.

### Viscous Burgers solution

The experiment uses $\Delta t=\Delta x=1/400$, $T=3$, $\nu=5\times10^{-4}$, and homogeneous Dirichlet boundaries. The maximum over the whole space–time domain is $1.045940$, and the final $L^\infty$ norm is $0.325871$.

![Recomputed space–time solution of the viscous Burgers equation](assets/pint/model-burgers.svg)

### Wave solution

The source-free experiment uses the trapezoidal rule with $\Delta t=\Delta x=1/400$, $T=3$, and $c^2=0.2$. The initial displacement reuses the oscillatory function of the previous two experiments, with zero initial velocity. The final displacement has $L^\infty$ norm $0.948217$.

![Recomputed space–time solution of the wave equation](assets/pint/model-wave.svg)

### Numerical summary

| Supplemental experiment | Grid and horizon                   |                      Reported metric |
| ----------------------- | ---------------------------------- | -----------------------------------: |
| advection–diffusion     | $\Delta x=\Delta t=10^{-3}$, $T=3$ | final $L^\infty=7.022\times10^{-99}$ |
| viscous Burgers         | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.325871$ |
| wave equation           | $\Delta x=\Delta t=1/400$, $T=3$   |            final $L^\infty=0.948217$ |

These metrics describe the dynamics of the discrete solution; they do not measure time-parallel iteration speed or hardware speedup.

## Completeness audit

| Item checked                                                                                                                                                    | Result             |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------ |
| The two semidiscrete systems, matrix dimensions, nonlinear example, domain-decomposition exception, and one-dimensional rationale in the Section 2 introduction | covered            |
| Advection–diffusion as a bridge (containing both parabolic/hyperbolic components) and the argument that $\nu$ tunes between the two ends                        | covered            |
| Figure 2.1(a–d): boundaries, source/data, time interval, causality discussion, room analogy, Robin comment, and low-frequency conclusion                        | covered            |
| Figure 2.2(a–h): two boundary types, three $\nu$ values, two initial-data experiments, precomputable interval, and Neumann footnote                             | covered            |
| Figure 2.3(a–h): linear reference, deformation, shock generation, high frequencies, precomputable interval, and limiting difficulty                             | covered            |
| Figure 2.4(a–d): three boundary types, bidirectional propagation, reflection, all frequencies, and cross-boundary conclusion                                    | covered            |
| The difficulty of the $\nu\to0$ hyperbolic limit and the emphasis on the "periodic + small diffusion" test condition                                            | covered            |
| Notation, domains, initial conditions, and numbering of equations (2.1)–(2.7)                                                                                   | checked            |
| Source results, site frequency/spectral interpretation and insights, synthesis diagram, and Python experiments                                                  | labeled separately |

## Summary

Chapter 2 uses three knobs—boundary conditions, diffusion strength, and nonlinearity—to progressively lengthen the temporal memory of the solution. The Dirichlet heat equation forgets old information quickly; the Neumann and periodic heat equations retain only low-frequency components such as the constant. Advection adds directed transport, and small diffusion together with periodic boundaries lets fine scales span very long times. Burgers nonlinearity continues to generate and sharpen shocks and high-frequency fronts. The wave equation has multidirectional propagation and reflection under all common boundaries, so long-term memory becomes a universal phenomenon and does not depend on the artificial condition of "periodic + small diffusion."

This progression from "temporal locality" to "long-range memory" explains the paper's subsequent classification: the methods of Chapter 3 are designed specifically to handle long-range propagation and are therefore often equally effective on parabolic problems; the methods of Chapter 4 mainly exploit the dissipation and temporal locality of parabolic problems, and fail once moved to hyperbolic problems.

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 2, pp. 388–396.
