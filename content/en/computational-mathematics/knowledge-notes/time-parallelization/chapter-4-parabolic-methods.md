---
title: "Chapter 4: PinT Methods Designed for Parabolic Problems"
description: Complete analysis of Parareal, PFASST, MGRiT, diagonalized Parareal, and space–time multigrid, including the design motivation of the Section 4 introduction and the 4.1 historical development
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-parabolic-methods
tags:
  - parallel-in-time
  - parabolic-PDE
---

> [!note] Reading scope
> This chapter covers Section 4 of the paper (pp. 443–481). The main text strictly retains the original hierarchy of the Section 4 introduction, 4.1 historical development, 4.2 Parareal, 4.3 PFASST, 4.4 MGRiT, the two diagonalized Parareal variants in 4.5, and 4.6 STMG. Python reproductions, parameter comparisons, and coverage audits are uniformly marked as site supplements and do not take paper section numbers. Theoretical contraction factors, measured iteration curves, and wall-clock performance carry distinct labels.

## Source-to-page map

This page retains the chapter-level synthesis and site reproductions. The equation-by-equation, theorem-by-theorem, and figure-by-figure arguments are split into the following pages:

| Paper section      | Source pages | Close-reading page                                                                                                                      | Coverage                                                 |
| ------------------ | ------------ | --------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| Section 4, 4.1–4.2 | pp. 443–452  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-1-parareal\|Historical development and Parareal]]         | (4.1)–(4.9), Theorems 4.1–4.4, Figures 4.1–4.5           |
| 4.3–4.4            | pp. 452–460  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-2-pfasst-mgrit\|PFASST and MGRiT]]                        | (4.10)–(4.13), Theorems 4.5–4.6, Figures 4.6–4.11        |
| 4.5                | pp. 460–472  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-3-diagonalized-parareal\|Diagonalization-based Parareal]] | (4.14)–(4.29), Theorems 4.7–4.8, Figures 4.12–4.17       |
| 4.6                | pp. 472–481  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-4-stmg\|Space–time multigrid]]                            | (4.30)–(4.44), Theorem 4.9, Figures 4.18–4.22, Table 4.1 |

## Section 4 introduction and 4.1 Historical development

### What parabolic methods exploit

Chapter 2 already explained intuitively that the solution of a parabolic problem is **local** in time: apart from a few very low-frequency components, information is quickly smoothed away by diffusion. The solution of a hyperbolic problem, by contrast, is highly **nonlocal** in time, and that nonlocality spans every frequency component from the lowest to the highest. The methods of Chapter 3 were designed precisely for hyperbolic problems, treating all frequencies nonlocally at once. They often perform even better on parabolic problems, because over a long time window the only components that truly require nonlocal treatment are the few very low-frequency modes that remain highly nonlocal in time.

The difficulty is that most Chapter 3 methods were **designed for linear problems** and are most effective in the linear case; each has a weakness against nonlinearity:

- SWR (especially optimized SWR) needs transmission/Robin parameters to be determined, and converges poorly without reasonable parameters;
- ParaExp is mainly applicable to linear problems;
- the nonlinear versions of ParaDiag I and II use Newton iteration as the outer solver, and nonlinearity slows Newton's convergence, possibly even diverging over longer time windows.

This chapter turns instead to methods **designed specifically for parabolic problems** that exploit their temporal locality directly: Parareal, PFASST, MGRiT, two forms of diagonalized Parareal, and space–time multigrid (STMG). They are equally effective for linear and nonlinear problems, and their convergence mechanism is **completely different** from that of the Chapter 3 methods; for exactly this reason, transferring them directly to hyperbolic problems usually converges very slowly or even diverges.

### Why dissipation makes cheap coarse propagation work

These methods share a single premise: dissipation. Diffusion damps high-frequency error at a rate proportional to the square of the frequency, so a later state is no longer sensitive to fine information from long ago—no matter how badly the coarse and fine propagators disagree on a high mode, that mode is flattened by its own strong dissipation and the error does not accumulate over time. As a result, a cheap coarse propagator that only needs to carry the **slow modes** correctly is enough to transmit the long-range causal information; the fine propagation merely has to fill in the local details concurrently within each large interval.

Parareal's per-mode convergence factor $\rho_\ell=\dfrac{|R_g-R_f^J|}{1-|R_g|}$ captures this mechanism precisely: the numerator measures the difference between the coarse and fine propagators, and the denominator $1-|R_g|$ is the coarse propagator's **dissipation margin**. Parabolic modes lie on the negative real spectrum, $|R_g|$ is clearly less than $1$, and the denominator holds the difference firmly in check; once one enters the hyperbolic limit and the modes approach the imaginary axis, $|R_g|\to1$, the denominator collapses, a coarse/fine phase mismatch of the same size is sharply amplified, and the factor crosses $1$. This also explains why the same coarse-grid correction mechanism fails on hyperbolic problems: without dissipation, the phase and propagation errors at all frequencies do not vanish on their own, while a coarse propagator accurate enough to represent those errors is often so expensive that it loses its cost advantage.

The first three methods (Parareal, PFASST, MGRiT) construct their coarse/fine hierarchy mainly in **time** and usually keep a single spatial grid; STMG instead coarsens space and time together and uses a block Jacobi smoother that can be parallelized across time.

![Coarse/fine hierarchy and parallel work in parabolic PinT methods](assets/diagrams/pint/en/parabolic-multilevel-map.svg)

> [!tip] Insight
> One can read $1-|R_g|$ as a "memory-forgetting rate." Parabolic problems forget quickly, so a cheap coarse model only needs to remember a few slow modes; hyperbolic problems forget almost nothing, so any phase information discarded by the coarse propagator persists over the long run and turns back against the iteration. This single thread uniformly explains the degradation of all four methods in this chapter under small viscosity near the hyperbolic limit, and it foreshadows the phenomenon seen throughout the Chapter 4 figures in which the convergence factors uniformly approach $1$ as $\nu$ decreases (Gander, Wu, and Zhou, 2025).

### 4.1 Historical development: from multiple shooting to space–time multigrid

**The Parareal line.** The first method to appear in this chapter is the Parareal scheme of Lions, Maday, and Turinici (2001) (already mentioned in Section 3.4 when discussing the nonlinear ParaExp variant). Although it was proposed independently in a virtual-control setting, its roots trace back to the earlier **multiple shooting** technique for evolution problems (Bellen and Zennaro, 1989; Chartier and Philippe, 1993). Saha, Stadel, and Tremaine (1997) had already given the same algorithm in solar-system simulations, replacing expensive propagation with a "coarse model" (rather than a coarse grid), and noted its connection to waveform relaxation. An even earlier precursor is Nievergelt (1964), though the method there was **not yet iterative**—it precomputes several candidate trajectories all at once and then stitches them together. In the two-plus decades since Parareal was proposed, its convergence has become well understood (Gander and Vandewalle, 2007; Gander and Hairer, 2008, 2014; Gander and Lunet, 2024), and it has become the template for deriving more efficient PinT methods.

**Variants grown from Parareal.** The literature contains many adaptations for different problems or purposes: the parallel implicit time-integrator algorithm PITA (Farhat and Chandesris, 2003; Farhat et al., 2006; Cortial and Farhat, 2009), the parallel full approximation scheme in space and time PFASST (Minion, 2011; Emmett and Minion, 2012; Minion et al., 2015), multigrid reduction in time MGRiT (Falgout et al., 2014; Dobrev et al., 2017; Hessenthaler et al., 2020), and combinations of Parareal with ParaDiag (Wu, 2018; Gander and Wu, 2020). Parareal-based PinT methods share one feature: they use two (or more) levels of grid in the **temporal** discretization while using only a single grid in the **spatial** discretization.

**The space–time multigrid line.** The idea of using multigrid in space and time simultaneously goes back to Hackbusch's (1984) parabolic multigrid, with an elegant analysis by Lubich and Ostermann (1987) in the form of multigrid waveform relaxation. Along this route, coarsening in the time direction was long difficult to make effective; important improvements for strongly advective problems appear in Vandewalle and Van de Velde (1994), Horton and Vandewalle (1995), Janssen and Vandewalle (1996), and Van Lent and Vandewalle (2002). Temporal coarsening was truly unlocked by Gander and Neumüller (2016): using only standard multigrid components, they introduced and analyzed a **block Jacobi smoother** along time as the crucial new ingredient—this is STMG. It is one of the most powerful PinT algorithms currently available for parabolic problems, with excellent strong and weak scalability (see also Neumüller and Smears, 2019). This chapter presents it last and shows that it is equally effective for nonlinear parabolic problems.

> [!tip] Insight
> The dividing line between the two historical threads is precisely "intrusiveness": Parareal → PITA → PFASST/MGRiT wrap existing coarse and fine propagators in a two-level temporal grid and stay as nonintrusive as possible; the space–time multigrid line of Hackbusch → Gander and Neumüller (2016) instead requires access to the full temporal discretization, the smoother, and the grid transfers—more intrusive, but in return delivering the theoretically most thorough and most scalable temporal coarsening. Understanding this tradeoff lets you predict which family to try first on an existing codebase.

## 4.2 Parareal

### Algorithm and parallel structure

Partition $[0,T]$ into $N$ large intervals. Let $\mathcal F$ be the accurate, expensive fine propagator and $\mathcal G$ the cheap coarse propagator. Parareal updates

$$
U_{n+1}^{k+1}
=\mathcal G(U_n^{k+1})
+\mathcal F(U_n^k)-\mathcal G(U_n^k). \tag{4.1}
$$

In iteration $k$, all $\mathcal F(U_n^k)$ can be computed simultaneously; the new $\mathcal G(U_n^{k+1})$ still propagate sequentially in $n$. The coarse term establishes the latest causal prediction, and the difference $\mathcal F(U_n^k)-\mathcal G(U_n^k)$ repairs the previous iteration's fine/coarse mismatch across intervals.

Parareal can be viewed as an approximate Newton method for multiple shooting, or as coarse-propagator preconditioned iteration of the lower-triangular all-at-once system. It is a nonintrusive algorithm: existing fine and coarse time integrators can be plugged in as long as they expose an interface that "advances one large interval from a given initial value."

### Linear modal analysis

For a spatial eigenmode $\lambda_\ell$, let one fine step of the fine propagator have stability function $R_f(\lambda_\ell\Delta t)$, with $J$ fine steps in each large interval, and let the coarse stability function be $R_g(\lambda_\ell\Delta T)$. Theorem 4.1 writes the error iteration as a strictly lower-triangular Toeplitz matrix. Strict lower-triangularity yields a finite-step property: in exact arithmetic the sequential fine solution is reached in at most $N$ iterations, and the first $k$ large intervals are already exact after iteration $k$.

Finite-step termination does not automatically bring speedup. If reaching the tolerance requires nearly $N$ iterations, each of which contains one sequential coarse propagation, the temporal-parallel gain is small.

Theorem 4.2 provides two complementary estimates. When the number of time intervals is small, a combinatorial upper bound reflects superlinear convergence that accelerates with $k$; over long time horizons, the more explanatory quantity is the per-mode linear factor

$$
\rho_\ell
=\frac{|R_g-R_f^J|}{1-|R_g|}. \tag{4.5b}
$$

The numerator measures the fine/coarse propagation difference, and the denominator measures the coarse propagator's dissipation margin. Figure 4.2 shows that the short-horizon superlinear phase and the long-horizon nearly linear phase can coexist.

### Nonlinear finite-step property

When the coarse propagator has Lipschitz stability and a local error of order $p$, Theorem 4.3 gives finite-step and superlinear error bounds for nonlinear Parareal. The core mechanism is still the interval-by-interval establishment of exactness. Nonlinearity makes the constants depend on the regularity of the solution and the Lipschitz constants of the propagators, so the observed iteration count can be highly sensitive to the time window and the dynamical parameters.

### How the time integrator sets the limiting factor

Theorem 4.4 studies parabolic modes on the negative real axis. If the fine propagation uses an L-stable method, the coarse propagation uses backward Euler, and each large interval contains enough fine steps, the worst long-time factor is about $0.3$. If the fine method is only A-stable, high-frequency modes are not sufficiently damped, and the required $J$ grows with the most dangerous frequency. Figure 4.3 compares combinations such as backward Euler, the trapezoidal rule, and Radau IIA, showing that once the coarse method is fixed, the fine method and the fine/coarse step ratio still change the convergence region significantly. A suitable Radau IIA combination can reduce the theoretical worst factor to about $0.068$.

### From the heat equation toward the hyperbolic limit

Figures 4.4–4.5 use periodic ADE and Burgers equations with $T=4$, $\Delta T=0.1$, $\Delta x=1/128$, $J=32$, backward Euler coarse propagation, and a second-order L-stable SDIRK fine propagation. As viscosity decreases, the fine/coarse phase mismatch of ADE grows; the Burgers equation additionally changes the local propagation speed and shock position. Both equations slow markedly, and Burgers shows approximate divergence near $\nu\le10^{-3}$.

The wave equation is more demanding. Unless the coarse propagation is almost as accurate as the fine propagation, Parareal can hardly control the phase error; a coarse propagator that reaches such accuracy usually also loses its cost advantage. Linear advection can be improved with a semi-Lagrangian or phase-optimized coarse propagator, but a general nonlinear hyperbolic coarse model remains an open problem.

## 4.3 PFASST

### From SDC to a temporal multilevel hierarchy

PFASST combines Parareal's large-interval concurrency with the high-order collocation solve of SDC. On each time step, choose $M_f$ fine collocation nodes and write the collocation equation as

$$
\boldsymbol U_f
=\boldsymbol U_{0,f}
+\Delta t\,Q_f\boldsymbol f(\boldsymbol U_f). \tag{4.10}
$$

Directly solving the dense collocation system is expensive. SDC uses an easily solved lower-triangular approximation based on implicit Euler as a preconditioning sweep, and each sweep eliminates part of the collocation residual. PFASST then introduces a coarse level with $M_c$ nodes, pipelining the fine and coarse sweeps across different time steps.

Restriction and prolongation between fine and coarse nodes use Lagrange interpolation. The coarse level uses the full approximation scheme (FAS), carrying the fine-level residual into the coarse collocation equation as a $\tau$ correction. One PFASST iteration comprises a parallel fine sweep, fine-to-coarse transfer, a sequential or pipelined coarse sweep, and coarse-to-fine correction. This structure can be understood either from Parareal or as multigrid on the collocation equations.

> [!note] Two equivalent formulations of PFASST
> This page follows the FAS/$\tau$-correction formulation, which is the **common way of naming** the full approximation scheme. The **precise characterization** given in the paper's §4.3 is the block-iteration form (Gander et al., 2023b): the two-level PFASST update is composed of three block operators $B_{01}$, $B_{10}$, $B_{00}$, written schematically as
>
> $$
> \boldsymbol U_{n+1}^{k+1}
> =B_{01}\,\boldsymbol U_{n+1}^{k}
> +B_{10}\,\boldsymbol U_{n}^{k+1}
> +B_{00}\,\boldsymbol U_{n}^{k}+\cdots,
> $$
>
> where each block is composed of the fine and coarse collocation matrices, their implicit-Euler approximation $\tilde M$, and the inter-node transfer matrices $\mathcal T_{f\to c}$, $\mathcal T_{c\to f}$ (for example, structures of the form $B_{10}=\mathcal T_{c\to f}\tilde M_c^{-1}\mathcal T_{f\to c}$). The two formulations point to the same algorithm: the block-iteration form is more convenient for convergence analysis, while FAS/$\tau$ explains its multigrid origin.

### Numerical observation and limitations

Figure 4.6 uses periodic heat and ADE equations with $T=3$, $\Delta x=1/128$, $\Delta t=1/64$, source parameter $\sigma=1000$, three Radau IIA nodes on the fine level, and two nodes on the coarse level. The heat equation converges rapidly; once ADE viscosity decreases, high-frequency propagation across time steps becomes increasingly difficult for the coarse collocation level to represent, and convergence slows accordingly.

PFASST suits scenarios that require high-order temporal accuracy and where the per-step collocation solve is expensive. Actual performance depends on the collocation nodes, the number of SDC sweeps, the coarse-level cost, node transfers, pipeline fill, and the allocation of spatial parallel resources. Reporting only the collocation order is not enough to judge parallel efficiency.

## 4.4 MGRiT

### F points, C points, and FCF relaxation

MGRiT marks one C point every $J$ fine points on the time grid, with the rest being F points. F relaxation advances fine propagation concurrently between neighboring C points; C relaxation updates the coarse points; coarse-grid correction handles long-distance propagation. A two-level FCF iteration can be written as an overlapping update similar to Parareal: first F relaxation, then C relaxation and a second F relaxation.

One FCF iteration uses roughly two sets of fine propagation and therefore usually costs more than one Parareal iteration. The extra CF segment provides overlap, so in theory each iteration can bring two large intervals into the exact region; at most about $\lceil N/2\rceil$ iterations reach the sequential fine solution. The more general $F(CF)^\nu$ trades more overlap for stronger contraction.

### Convergence factor and fair-cost comparison

Theorem 4.5 gives the long-time modal factor for two-level FCF:

$$
\rho_{\mathrm{MGRiT},\ell}
=|R_f^J|
\frac{|R_g-R_f^J|}{1-|R_g|}.
$$

It has one more $|R_f^J|$ than the Parareal factor, so dissipative modes contract additionally. Each additional CF segment multiplies by a similar fine-propagation factor while adding another set of fine-solve cost.

A fair comparison places one FCF iteration and two Parareal iterations at comparable fine-solve workloads. Theorem 4.6 gives representative worst factors for L-stable fine methods. With backward Euler coarse propagation, Parareal and one FCF are about $0.2984$ and $0.1115$; a second-order Lobatto IIIC combination is about $0.0817$ and $0.0197$. One FCF has the smaller per-iteration factor, but it can still be slightly worse than the square of the Parareal factor. Figure 4.8 further displays this equal-cost relationship through convergence regions in the complex plane.

### ADE and Burgers experiments

Figures 4.9–4.10 use $T=5$, $J=20$, $\Delta T=1/8$, $\Delta x=1/160$, backward Euler coarse propagation, and SDIRK22 fine propagation. The per-mode factors of the heat equation are far below 1. As ADE viscosity decreases from $0.1$ to $0.01$ and $0.002$, the dangerous modes gradually approach the high-phase, low-dissipation region. At $\nu=0.002$, the linear-stage factors of Parareal and FCF are about $1.4211$ and $1.2812$ respectively, and both grow transiently.

MGRiT's finite-step decay still holds; it only states that the strictly lower-triangular error is eventually cleared and cannot remove the amplification in the first few dozen iterations. If the actual tolerance requires stopping during the transient-growth phase, the method has not formed a scalable solver.

Figure 4.11 obtains the same main pattern for the nonlinear Burgers equation. With adequate diffusion, FCF outperforms Parareal per iteration; once balanced by the number of fine propagations, one FCF is often comparable to two Parareal iterations. When the coarse propagation cannot accurately represent the nonlinear shock, both degrade.

## 4.5 Diagonalization-based Parareal

The paper discusses two different locations for diagonalization. The first rewrites the sequential coarse-grid correction over the $N$ large intervals as an FFT-diagonalizable head–tail system; the second, inside each large interval, constructs a low-cost coarse propagator from an $\alpha$-circulant system over the $J$ fine time points.

### 4.5.1 Diagonalization-based CGC

The coarse correction of standard Parareal is carried out sequentially in time. The paper adds the condition

$$
U_0^{k+1}=u_0+\alpha U_N^{k+1}
$$

and modifies the correction right-hand side accordingly so that the original initial-value problem remains the fixed point of the iteration. The linear backward Euler case produces an $\alpha$-circulant all-at-once matrix that can be reduced by FFT to $N$ independent shifted coarse spatial problems.

As $\alpha\to0$ the standard sequential coarse correction is recovered, but the roundoff error of the transform grows. Theorem 4.7 shows that if the convergence factor of standard Parareal is $\rho$, then choosing

$$
\alpha\le \frac{\rho}{1+\rho}
$$

preserves convergence of the same order. In Figure 4.12, the heat equation measures $\rho\approx0.22$ with a threshold of about $0.18$; ADE at $\nu=0.1$ measures $\rho\approx0.39$ with a threshold of about $0.28$, consistent with the theory.

The nonlinear version applies quasi-Newton linearization to the all-at-once coarse equation and then builds a block $\alpha$-circulant system from an averaged Jacobian. The Burgers experiment in Figure 4.13 shows the same $\alpha$-threshold phenomenon.

Directly inserting the head–tail condition into MGRiT changes its fixed point and causes divergence. The paper gives a consistent alternative condition in which the head–tail coupling uses the difference between the two most recent iterates. At small $\alpha$, it retains the same convergence rate as the original MGRiT; this condition can also be used in Parareal.

### 4.5.2 Diagonalization-based coarse propagator

The second method has the fine and coarse propagation use the same time integrator and the same small step. The fine propagation executes $J$ small steps sequentially; the coarse propagation places these $J$ unknown time points into an $\alpha$-circulant all-at-once system and solves them in parallel. At $\alpha=0$ the coarse propagation equals the fine propagation, and Parareal converges in one iteration but without a parallel cost advantage; at $\alpha>0$ a diagonalizable approximation results, and the coarse propagation can be about $J$ times cheaper.

Theorem 4.8 gives a clear distinction by equation type. A parabolic problem with negative real spectrum has the factor

$$
\rho=\alpha,
$$

and it is independent of the number $N$ of large intervals. A wave problem with purely imaginary spectrum satisfies the upper bound

$$
\rho\le \frac{2\alpha N}{1+\alpha}.
$$

Figure 4.14 verifies the $\alpha$ slope for the heat equation. The wave experiment in Figure 4.15 shows that at $\alpha=0.01$, increasing $N$ slows convergence; the actual curve still has a superlinear phase, and increasing $N$ from 24 to 960 costs only about two more iterations to reach discretization error. The linear upper bound is closest to reality in the region of small $\alpha N$.

Figures 4.16–4.17 apply the method to the Burgers equation. The error factor is approximately linear in $\alpha$, and at $\alpha=10^{-3}$ it is fairly robust with respect to the number of time intervals; decreasing $\alpha$ also mitigates the degradation caused by small viscosity. This coarse propagator works for both parabolic and hyperbolic problems and is the more broadly applicable of the two diagonalized Parareal variants.

## 4.6 Space–time multigrid

### All-at-once system and block Jacobi smoothing

For the linear system $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$, the $\theta$ method forms the all-at-once system

$$
K\boldsymbol U
=\left(B\otimes I_x-\Delta t\,\widetilde B\otimes A\right)\boldsymbol U
=\boldsymbol b. \tag{4.31}
$$

STMG uses damped block Jacobi:

$$
\boldsymbol U^{j+1}
=\boldsymbol U^j
+\eta\left(I_t\otimes r_1\right)^{-1}
(\boldsymbol b-K\boldsymbol U^j),
\qquad r_1=I_x-\theta\Delta t A. \tag{4.32}
$$

$I_t\otimes r_1$ is block diagonal in time, so the spatial blocks at all time points can be solved in parallel. The damping $\eta$ determines the rate at which high-frequency error is reduced.

A two-level cycle first performs $s_1$ presmoothing steps, restricts the residual simultaneously to the coarse spatial and coarse temporal grids, solves the coarse-grid equation, prolongs the correction, and then performs $s_2$ postsmoothing steps. Recursive application yields the full STMG. It coarsens $\Delta X=2\Delta x$ and $\Delta T=2\Delta t$ together, handing the long-wave error in both space and time to the coarse level.

Early parabolic multigrid often used pointwise Gauss–Seidel advancing along time. That smoother is essentially sequential forward substitution: although spatial coarsening is effective, the time direction does not scale. Block Jacobi changes the location of parallelism and is the key to modern STMG.

### Local Fourier analysis

For the one-dimensional heat equation, centered spatial differences, and backward Euler, local Fourier analysis decomposes the error into space–time frequencies and derives the smoothing symbol. Theorem 4.9 gives the optimal damping

$$
\eta=\frac12.
$$

The smoothing factor for high temporal frequencies is at most $1/\sqrt2$, so temporal coarsening can be done safely. In the paper's normalized heat equation, if $\Delta t/\Delta x^2\ge1/\sqrt2$, high spatial frequencies are also pushed under the same upper bound, permitting spatial coarsening at the same time; after restoring the diffusion coefficient, the corresponding nondimensional ratio is $\nu\Delta t/\Delta x^2$. The ADE symbol contains an imaginary part, and on this basis the paper still recommends $\eta=1/2$ as a sensible starting point for backward Euler.

### Damping, smoothing count, and the time integrator

Figure 4.19 scans $\eta$ for two-level backward-Euler STMG. Both the heat equation and ADE at $\nu=0.01$ perform well near $1/2$. Figure 4.20 fixes $\eta=1/2$ and compares one versus three block Jacobi smoothing steps. More smoothing makes each cycle more expensive while significantly reducing the number of cycles; ADE is slower than the heat equation, but as the number of smoothing steps increases its sensitivity to viscosity weakens and a superlinear phase appears.

After Figure 4.21 switches to the trapezoidal rule, the heat equation can hardly converge stably even with ten smoothing steps; ADE can converge, more smoothing improves the speed, and the better damping in the sampling is about $0.8$. This shows that the effectiveness of STMG depends on the stability function of the time integrator. The theoretical damping for backward Euler cannot be transplanted directly to the trapezoidal rule.

### Large-scale scaling results

Table 4.1 summarizes modern STMG data on the three-dimensional heat equation. In weak scaling, the number of cores increases from 1 to 262,144, the time steps from 2 to 524,288, and the degrees of freedom from 59,768 to 15,667,822,592; the iteration count stays at 7, and the total time changes from about 28.8 seconds to 30.0 seconds. The sequential-time-marching column, parallel only in space, grows with the problem size to about 4,988,060 seconds.

In strong scaling, a problem with 512 time steps and 15,300,608 degrees of freedom drops from 7,635.2 seconds on 1 core to 30.0 seconds on 256 cores; a larger 524,288-time-step problem drops from 15,205.9 seconds on 512 cores to 30.0 seconds on 262,144 cores. These data come from the three-dimensional implementation cited by the paper and demonstrate the weak- and strong-scaling potential of STMG for parabolic problems.

### Nonlinear FAS-STMG

The nonlinear semi-discrete system $\boldsymbol u'=\boldsymbol f(\boldsymbol u)$, after the $\theta$ method, forms

$$
K(\boldsymbol U)
=(B\otimes I_x)\boldsymbol U
-\Delta t(\widetilde B\otimes I_x)\boldsymbol f(\boldsymbol U)
=\boldsymbol b. \tag{4.42}
$$

Nonlinear block Jacobi solves the local correction in parallel on each time block, optionally using Newton internally. Local Fourier analysis no longer applies directly, and the damping needs to be tested per problem. The coarse level uses FAS: restrict the current solution and residual, solve the nonlinear coarse problem with a $\tau$ consistency correction, prolong the coarse-solution difference, and then postsmooth.

Figure 4.22 uses two block Jacobi smoothing steps for the Burgers equation. With adequate diffusion, nonlinear STMG converges rapidly; once viscosity decreases, it deteriorates markedly. The better damping in the experiment is $\eta=1/4$, reflecting the difference between the nonlinear case and the linear backward-Euler theory.

Taking Section 4.6 as a whole, STMG is one of the most powerful PinT methods currently available for parabolic problems and has demonstrated large-scale scaling. It requires access to the full temporal discretization, the smoother, and the grid transfers, making it more intrusive than Parareal. Its robustness for hyperbolic problems and across different time integrators still requires further study.

## Site numerical supplement: Python reproduction results

### Parareal baselines and Figure 4.5

| Problem             | Formal parameters                     | Iterations |   Final maximum error |
| ------------------- | ------------------------------------- | ---------: | --------------------: |
| advection–diffusion | $N_x=128$, $N=40$, $J=32$, $\nu=0.02$ |         39 | $6.106\times10^{-15}$ |
| Burgers             | $N_x=128$, $N=40$, $J=32$, $\nu=1$    |         16 | $3.544\times10^{-12}$ |

![Parareal convergence baseline on the advection–diffusion equation](assets/pint/parareal-ade-baseline.svg)

![Parareal convergence baseline on the viscous Burgers equation](assets/pint/parareal-burgers-baseline.svg)

With $T=4$, $\Delta T=0.1$, $\Delta x=1/128$, and $J=32$ fixed, the iterations required to reduce the maximum error below $10^{-10}$ are:

| Equation            | $\nu=1$ | $\nu=0.1$ | $\nu=0.02$ |
| ------------------- | ------: | --------: | ---------: |
| advection–diffusion |      14 |        24 |         35 |
| Burgers             |      14 |        21 |         25 |

![Parareal convergence for ADE and Burgers as diffusion weakens](assets/pint/parareal-figure-4-5.svg)

These counts measure convergence toward the sequential fine solution. GPU performance appears in [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view#gpu-acceleration-and-profiling|Chapter 5]].

### MGRiT baseline and Figures 4.9–4.10

The near-hyperbolic baseline uses $N_x=160$, $N=40$, $J=20$, $T=5$, $\nu=0.002$. At the end of the reported window, Parareal has a maximum error of $2.895\times10^2$; two-level MGRiT eventually drops to $5.551\times10^{-16}$ through the finite-step property.

![Baseline comparison of Parareal and MGRiT on near-hyperbolic ADE](assets/pint/mgrit-baseline.svg)

Comparing fine-solve cost as "one FCF versus two Parareal iterations," the per-mode long-time factors are:

| Problem          | Parareal factor | FCF-MGRiT factor |
| ---------------- | --------------: | ---------------: |
| heat             |          0.2824 |           0.0835 |
| ADE, $\nu=0.1$   |          0.4453 |           0.2719 |
| ADE, $\nu=0.01$  |          1.0501 |           0.9021 |
| ADE, $\nu=0.002$ |          1.4211 |           1.2812 |

![Per-mode long-time convergence factors for heat and ADE at different viscosities](assets/pint/mgrit-figure-4-9.svg)

![Equal-fine-solve-cost curves of Parareal and FCF-MGRiT](assets/pint/mgrit-figure-4-10.svg)

> [!note] A numerical difference from Figure 4.9
> The original paper's figure annotates the maximum Parareal factor at $\nu=0.01$ as $0.9986$; the Python conversion project, computing directly from the paper's equation (4.5b), the printed parameters, and the stability function in the upstream `MGRiT_Heat_ADE.m`, obtains $1.0501$. The MGRiT factor $0.9021$ for the same case and the other two pairs of factors all agree closely with the original figure. This page retains the reproduced value and records the difference explicitly, rather than overwriting the computed result with the figure's annotation. Both numbers lie near $1$ and support the same qualitative conclusion: long-time convergence is already close to a critical state at this viscosity.

### STMG damping and smoothing validation

The trapezoidal-rule baseline uses $N_x=N_t=255$, $\nu=10^{-3}$, and three pre- and postsmoothing steps. In the sampling, the minimum error after 15 cycles occurs at $\eta\approx0.98$.

![Baseline scan of the STMG damping parameter with the trapezoidal rule](assets/pint/stmg-baseline.svg)

The paper-grid validation for backward Euler gives:

| Problem         | Best sampled $\eta$ after 15 iterations |
| --------------- | --------------------------------------: |
| heat            |                                   0.500 |
| ADE, $\nu=0.01$ |                                   0.372 |

![Damping-parameter scan for backward-Euler STMG](assets/pint/stmg-figure-4-19.svg)

The discrete-sampling minimum for ADE is $0.372$, whereas the paper's theory and plot take $1/2$ as the robust choice. The two measure different things—the 15-cycle error on a particular finite grid versus a high-frequency smoothing bound—so their conclusions are at different levels.

With $\eta=0.5$ fixed, three pre- and postsmoothing steps reduce the error in fewer cycles than one smoothing step:

![STMG convergence comparison of one versus three pre- and postsmoothing steps](assets/pint/stmg-figure-4-20.svg)

A real performance decision should compare

$$
\text{work to tolerance}
=(\text{cost per cycle})\times(\text{cycle count})
+\text{communication and memory cost}.
$$

## Site method comparison: parameters and failure modes

| Method                | Key parameters                                           | Tradeoff controlled                                        | Typical failure mode                                                 |
| --------------------- | -------------------------------------------------------- | ---------------------------------------------------------- | -------------------------------------------------------------------- |
| Parareal              | interval count $N$, fine/coarse step ratio, $\mathcal G$ | concurrency, coarse cost, and propagation accuracy         | iteration count approaches $N$; phase error grows                    |
| PFASST                | collocation nodes, SDC sweeps, coarse nodes              | high-order accuracy, pipeline depth, and coarse-level cost | coarse collocation level distorts propagation under weak diffusion   |
| MGRiT                 | coarsening factor, F/FCF/$F(CF)^\nu$                     | overlapping contraction versus fine-solve work             | equal-cost advantage disappears; transient amplification             |
| diagonalized Parareal | $\alpha$, location of diagonalization                    | sequential coarse-correction fraction versus roundoff      | too-small $\alpha$ is ill-conditioned; too-large changes convergence |
| STMG                  | $\eta$, smoothing count, coarsening, integrator          | high-frequency smoothing versus cycle cost                 | integrator mismatch; hyperbolic/low-viscosity degradation            |

## Source coverage audit

| Source location                 | This page | Material covered                                                                                                                                                                                    |
| ------------------------------- | --------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Sections 4 and 4.1, pp. 443–444 | 4.1       | parabolic temporal locality, the nonlinear limitations of the Chapter 3 methods, the mechanism by which dissipation supports coarse correction, and the two historical threads of Parareal and STMG |
| Section 4.2, pp. 444–452        | 4.2       | Parareal update, Theorems 4.1–4.4, Figures 4.1–4.5, nonlinear and hyperbolic degradation                                                                                                            |
| Section 4.3, pp. 452–455        | 4.3       | collocation equations, SDC, FAS transfer, block-iteration (B01/B10/B00) characterization, PFASST iteration and Figure 4.6                                                                           |
| Section 4.4, pp. 455–460        | 4.4       | FCF structure, Theorems 4.5–4.6, Figures 4.7–4.11, cost-fair comparison, and Burgers                                                                                                                |
| Section 4.5, pp. 460–472        | 4.5       | two diagonalization locations, Theorems 4.7–4.8, Figures 4.12–4.17, nonlinearity and the consistent MGRiT condition                                                                                 |
| Section 4.6, pp. 472–481        | 4.6       | all-at-once STMG, block Jacobi, Theorem 4.9, Figures 4.18–4.22, Table 4.1, nonlinear FAS                                                                                                            |

## Summary

This group of methods turns the diffusion-induced temporal locality into an algorithmic advantage. Parareal corrects parallel fine propagation with cheap coarse propagation; PFASST superimposes SDC and FAS on collocation nodes; MGRiT strengthens coarse correction with a temporal hierarchy and overlapping relaxation; the diagonalized variants reduce sequential coarse propagation; STMG treats spatial and temporal scales together. They can be highly efficient on parabolic problems, with performance still depending on the spectral matching of coarse and fine propagation, the time integrator, cost-fair comparisons, and full hardware overhead. As viscosity falls, phase and long-range memory again become dominant, and standard coarse-grid mechanisms progressively lose their advantage.

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4, pp. 443–481.
