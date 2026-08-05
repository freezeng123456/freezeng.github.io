---
title: Parallel-in-Time Algorithms
description: Thirteen papers that rewrite sequential time stepping as a diagonalisable or preconditionable algebraic problem
lang: en
translation: computational-mathematics/paper-notes/parallel-in-time
tags:
  - computational-mathematics
  - paper-notes
  - parallel-in-time
---

This topic holds 13 papers written between 2015 and 2025, mostly with Shu-Lin Wu, and summarised in a 2025 _Acta Numerica_ survey (paper 85). They all answer one question: **time stepping is inherently sequential, so how can that sequential structure be rewritten as an algebraic problem that admits concurrency?**

![Replace the sequential recurrence by a diagonalisable time matrix](assets/diagrams/tao-zhou-papers/en/pint-diagonalization.svg)

> [!info] Relation to the survey close reading on this site
> The section-by-section close reading of paper 85 is a separate topic here: [[en/computational-mathematics/knowledge-notes/time-parallelization/index|Time Parallelization for Hyperbolic and Parabolic Problems]], including every original figure and reproducible Python experiments. This page does not repeat it; it explains what each of the twelve original papers underneath the survey actually solved.

## Two technical routes

### The iterative route: coarse and fine propagators

The parareal iteration

$$
\boldsymbol u_{n+1}^{k+1}=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^{k})
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k})
$$

parallelises the expensive fine propagation $\mathcal F$ at the price of a coarse correction that still advances sequentially. For a linear problem $\boldsymbol u'+A\boldsymbol u=g$ the error decouples by mode, and the contraction is governed by

$$
\varrho_l(J,z)=\frac{\bigl|R_g(z)-R_f^{J}(z/J)\bigr|}{1-|R_g(z)|},
\qquad z=\Delta T\lambda,\ \lambda\in\sigma(A),
$$

where $J=\Delta T/\Delta t$ is the coarsening ratio and $R_g,R_f$ are the stability functions of the two propagators. The numerator is the coarse-fine discrepancy and the denominator is the dissipation margin of the coarse propagator.

### The direct route: diagonalising the all-at-once system

Stacking the unknowns at all time levels, a one-step discretisation gives a block lower-triangular Toeplitz system. If the temporal matrix is diagonalisable, $B=VDV^{-1}$, then multiplying by $V^{-1}$ decouples it into **independent complex-shifted spatial problems** that can be solved fully concurrently before multiplying back by $V$. This route has no outer iteration; the cost moves into the conditioning of $V$.

## Three close readings

| Close reading                                                                                                                      | Papers             | Technical core                                    |
| ---------------------------------------------------------------------------------------------------------------------------------- | ------------------ | ------------------------------------------------- |
| [[en/computational-mathematics/paper-notes/parallel-in-time/parareal-convergence\|Convergence analysis for parareal]]              | 12, 20, 30, 77     | contraction factors and fine-propagator stability |
| [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique\|The diagonalisation technique]]             | 31, 39, 46, 53     | making the temporal matrix diagonalisable         |
| [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners\|All-at-once preconditioners and spectra]] | 59, 65, 71, 84, 85 | the spectrum of the preconditioned system         |

## Core idea of each paper

### Convergence analysis for parareal

- **12 (convergence analysis for three parareal solvers)** covers a practically important case that had no theory. When the fine propagator is only A-stable rather than L-stable — most importantly the trapezoidal rule and the fourth-order two-stage Gauss Runge-Kutta method — $|R_f(z)|\to1$ as $z\to\infty$, so the classical argument producing the $\approx0.3$ constant breaks down and the contraction factor can approach $1$ when $\lambda_{\max}$ is large. Earlier analyses assumed either that $\mathcal F$ was the exact propagator or that it was L-stable, and so said nothing about the A-stable-only case, especially at small coarsening ratio $J$. The paper's remedy is a **threshold**: an only-A-stable fine propagator needs $J_{\min}=O(\log^2 z_{\max})$ to restore robust contraction, and the critical ratio for a third-order DIRK method is $J_{\rm cri}=4$.
- **20 (fast parareal iterations for fractional diffusion)** and **30 (parareal with local time integrators)** treat **two different** kinds of fractional nonlocality, and their technical difficulties do not overlap. Paper 20 is **space**-fractional: the nonlocality is in space, time remains local, the parareal structure survives intact, and the cost lands on the dense $A$ — which makes the coarse propagator the bottleneck. Paper 30 is **time**-fractional: the historical effect destroys the independence of the subintervals, the work grows monotonically along the time axis, and load balancing becomes part of the algorithm design itself.
- **77 (convergence analysis with a nonuniform fine time grid)** drops the standard assumption that the fine grid inside a coarse interval is uniform. That assumption fails routinely under adaptive time stepping, and dropping it breaks the structure the modal decoupling above relies on: the fine propagator is no longer the $J$-th power $R_f^J(z/J)$ of a single stability function but a product of distinct factors, and if the fine grid varies from one coarse interval to the next then $M_f(z)$ is not even Toeplitz, so the engine behind every convergence-factor argument stops working.

### The diagonalisation technique

- **31 (time-periodic fractional diffusion)** combines diagonalisation with multigrid. The all-at-once matrix of a time-periodic problem has a different structure from an initial-value problem, and that difference is exploitable: periodic coupling makes the time matrix a **genuine circulant**, diagonalised by the **unitary** DFT matrix, so $\mathrm{Cond}_2(V)=1$. **The entire price of this route drops to zero here** — no $\alpha$, no outer iteration, no roundoff amplification. The second half of the paper handles the new problem this creates: each decoupled mode leaves a complex-shifted dense fractional-Laplacian system.
- **39 (accelerating two-level MGRIT by diagonalisation)** fixes two inefficiencies at once. First, MGRIT almost always takes backward Euler as its coarse propagator, which is only first order and gives a contraction of about $0.1$ under FCF-relaxation; the paper switches to second-order Lobatto IIIC, the $(0,2)$ Padé approximant, which is A-stable and L-stable and damps high frequencies far more aggressively than backward Euler's $1/(1+z)$ — exactly what a coarse propagator needs. Second, and structurally, the coarse-grid correction is an inherently **sequential** sweep through $N_t$ coarse levels and becomes the serial bottleneck as the processor count grows; the paper introduces a head-tail coupling condition $\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0$, turning the coarse correction into a diagonalisable all-at-once system

  $$
  \bigl(C_\alpha\otimes I_x-I_t\otimes\Delta T A\bigr)\boldsymbol U^{k+1}=\boldsymbol g^k,
  $$

  where $C_\alpha$ is an $\alpha$-circulant matrix diagonalised by FFT. The two fixes are compatible: a better coarse propagator normally worsens the sequential bottleneck, and parallelising the coarse correction cancels that. The four directly comparable contraction constants (L-stable $\mathcal F$, $J=O(1)$) are $0.2984$ and $0.1115$ for a backward-Euler coarse propagator with parareal and MGRIT-FCF respectively, and $0.0817$ and $0.0197$ for Lobatto IIIC-2. **On an equal-cost accounting the robust gain comes from changing the coarse propagator, not from the extra relaxation**: FCF performs two fine solves per sweep, and $0.2984^2=0.0890<0.1115$ while $0.0817^2=0.0067<0.0197$.

- **46 (parabolic PDE-constrained optimisation)** applies diagonalisation to the forward-backward optimality system of an optimal control problem, where the temporal coupling runs in both directions. The decisive step is to use the time matrix of the **time-periodic** version as a preconditioner for the **initial-value** problem: the periodic matrix is circulant and hence freely diagonalisable, the initial-value matrix is a Jordan block and hence not diagonalisable at all, and the two differ in a single corner entry.
- **53 (parallel implementation of two-stage SDIRK methods)** carries $\alpha$-circulant preconditioning to a **multi-stage** integrator. The difficulty is that implicit Runge-Kutta methods have no ready-made all-at-once form to stack, so the paper rebuilds a block preconditioner "with completely different structures and different implementation details" and proves an $O(\alpha)$ contraction for two-stage SDIRK when $\gamma\ge1/4$ — a condition that coincides exactly with A-stability for that method, which is what prompted the conjecture that A-stability alone should suffice.

### All-at-once preconditioners and spectra

- **59 (a well-conditioned direct PinT algorithm)** returns the all-at-once solution with no outer iteration, focusing on keeping the conditioning of the matrices involved under control — which is exactly where the geometric-step diagonalisation route pays its price. The move is to **choose the time discretisation backwards**: adopt a **boundary value method** that cannot be run as a time-stepping scheme at all, which is diagonalisable on a uniform grid and, through a connection with Chebyshev polynomials, satisfies $\mathrm{Cond}_2(V)=O(n^2)$ — turning the geometric-step route's **exponential** roundoff degradation in $N_t$ into a **polynomial** one. The paper reports over 60 times speedup on 256 cores.
- **65 (a uniform spectral analysis for a preconditioned all-at-once system)** gives one spectral analysis covering both first- and second-order evolutionary problems. The spectral distribution governs Krylov convergence, so this is the step that turns "here is a preconditioner" into "here is why it works". Its contribution is an **exchange of hypotheses**: from "the time-stepping matrix is sparse, Toeplitz and diagonally dominant" to "the method is stable", which yields for every stable one-step integrator the mesh-independent bound $\frac{1}{1+\alpha}\le|\lambda(\mathcal P_\alpha^{-1}\mathcal K)|\le\frac{1}{1-\alpha}$ and the corollary $\rho\le\alpha/(1-\alpha)$.
- **71 (a PinT preconditioner for forward-backward evolutionary equations)** handles the forward-backward coupled systems arising in optimal control and source identification by circulantising the **several** Toeplitz blocks of the system; the authors themselves identify the result as a parallel version of the matching Schur complement preconditioner. The price is that $\alpha$ can no longer be fixed freely and must shrink with $N_t$.
- **84 (a PinT preconditioner for time spectral methods)** deals with time spectral discretisations, whose all-at-once temporal blocks are **dense** rather than bidiagonal, so the recipe "replace the Toeplitz matrix by its $\alpha$-circulant counterpart" has nothing to act on. The paper's answer is to **factorise first and circulantise second**: factor Toeplitz pieces out of the dense $M$ and apply the Strang-type $\alpha$-circulant substitution only to those.
- **85 (the _Acta Numerica_ survey)** organises this work and the wider literature into two classes: methods that remain effective for propagative problems, and methods designed primarily for dissipative problems.

> [!note] Coverage status
> The formulas of papers **12, 39, 59 and 65** have been checked equation by equation, so their close readings carry full derivations, theorems with hypotheses, and verifiable numerical results. For papers **30, 31, 46, 53 and 84** the bibliographic record and abstract are verified, so the close readings give the setting, the construction and the **form** of the theorems, but report no specific constants from inside them. The abstracts and full texts of papers **20 and 77** are unobtainable through public channels — SIAM does not deposit abstracts with Crossref and the other aggregators return nulls for those DOIs — and paper **71**'s abstract is rendered with all inline mathematics deleted in every public source, so its scaling law for $\alpha$ is recorded only as $\alpha=CN_t^{-\theta}$ **with no exponent $\theta$**. No theorems, constants or experiments are supplied on their behalf. The complete section-by-section reading of paper 85 is in the [[en/computational-mathematics/knowledge-notes/time-parallelization/index|time-parallelization topic]].

## What numerical evidence can be verified

The three close readings collect every experimental datum that could be verified into tables. Taken together, the numerical evidence across these papers is very unevenly distributed, and it helps to know that in advance:

| Paper                  | Strongest verifiable evidence                                                                                                                               | What is missing                                              |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| 12                     | convergence constants and the thresholds $J_{\min}=O(\log^2z_{\max})$, $J_{\rm cri}=4$; the examples only to the level of application class                 | meshes, iteration counts                                     |
| 39                     | four directly comparable constants and the equal-cost accounting; the survey's viscosity scan showing where the constants fail                              | the scale and wall-clock data of the paper's own experiments |
| 59                     | over $60\times$ speedup on 256 cores; the baseline's $n\approx20$–$25$ ceiling; observed $\mathrm{Cond}_2(V)=O(n^{1.75})$                                   | test PDEs, meshes, wall-clock tables                         |
| 65                     | three qualitative observations (nonlinear GMRES, the effect of $T$, the Numerov stability boundary)                                                         | iteration counts, mesh sizes                                 |
| 71                     | the outer solver is CG, so the preconditioned system is symmetric positive definite; quantitative agreement between predicted spectrum and measured CG rate | $\gamma$, $N_t$, $C$, iteration counts                       |
| 84                     | six figures, GMRES as the outer solver                                                                                                                      | test problems, spectral degrees, clustering radius           |
| 20, 30, 77, 31, 46, 53 | see the "numerical experiments" sections of each close reading; most can be confirmed only to the level of application class                                | all quantitative results                                     |

One practical way to read this: **the quantitative strength of these papers lies in their theoretical constants, not in performance tables.** The only one carrying parallel performance numbers is paper 59, the only one carrying a set of directly comparable convergence constants is paper 39, and for most of the rest all that can be confirmed is that experiments exist and which outer solver was used — which is itself informative, since choosing CG over GMRES says directly whether the object being iterated on is symmetric positive definite.

## One trade-off that runs through everything

The price of the diagonalisation route can be stated exactly. Geometrically increasing steps $\Delta t_j=\Delta t_1\tau^{\,j-1}$ with ratio $\tau>1$ make the eigenvalues of the temporal matrix distinct and hence diagonalisable; but the closer $\tau$ is to $1$, the worse the conditioning of the eigenvector matrix $V$ and the more roundoff is amplified, while the larger $\tau$ is, the coarser the late steps and the larger the truncation error. This is a genuine dilemma, and it explains why the second route exists: with uniform steps, an $\alpha$-circulant approximation of the temporal matrix is diagonalised by FFT with conditioning controlled by $\alpha$, at the price of introducing an approximation and therefore needing an outer Krylov iteration. **The first is a direct solver, the second a preconditioner**, and that is the division of labour between the two groups of papers here.

Ordering that price by how fast it grows with the length of the time window brings out the economics of the whole thread. The geometric-step roundoff bound carries a factor $\varrho^{-(N_t-1)}$, where $\varrho=\tau-1$ is the amount of stretching — an **exponential** degradation, which is why it is a short-window method; paper 59's boundary value method brings that down to $\mathcal O(n^2)$, a **polynomial** one; fixed-$\alpha$ $\alpha$-circulant preconditioning gives $\mathrm{Cond}_2(V)\le1/\alpha$, which **does not grow with $N_t$ at all**; and paper 71's forward-backward coupling pushes it back to polynomial, because $\alpha$ must shrink with $N_t$. **The fixed-$\alpha$ tier is the only one with no window-length cost, and it pays for that by being only an approximation and needing an outer iteration.**

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]; per-paper references appear at the end of each close-reading page. The complete reading of the survey itself is in [[en/computational-mathematics/knowledge-notes/time-parallelization/index|Time Parallelization for Hyperbolic and Parabolic Problems]].
