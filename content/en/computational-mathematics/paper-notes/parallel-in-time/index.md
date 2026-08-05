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

- **12 (convergence analysis for three parareal solvers)** covers a practically important case that had no theory. When the fine propagator is only A-stable rather than L-stable — most importantly the trapezoidal rule and the fourth-order two-stage Gauss Runge-Kutta method — $|R_f(z)|\to1$ as $z\to\infty$, so the classical argument producing the $\approx0.3$ constant breaks down and the contraction factor can approach $1$ when $\lambda_{\max}$ is large. Earlier analyses assumed either that $\mathcal F$ was the exact propagator or that it was L-stable, and so said nothing about the A-stable-only case, especially at small coarsening ratio $J$.
- **20 (fast parareal iterations for fractional diffusion)** and **30 (parareal with local time integrators)** carry this route to time-fractional problems. The nonlocality of a fractional derivative undermines the premise that fine propagation inside each coarse interval is independent, and handling that is the technical focus of both.
- **77 (convergence analysis with a nonuniform fine time grid)** drops the standard assumption that the fine grid inside a coarse interval is uniform. That assumption fails routinely in practice, and dropping it breaks the structure the modal decoupling above relies on.

### The diagonalisation technique

- **31 (time-periodic fractional diffusion)** combines diagonalisation with multigrid. The all-at-once matrix of a time-periodic problem has a different structure from an initial-value problem, and that difference is exploitable.
- **39 (accelerating two-level MGRIT by diagonalisation)** fixes two inefficiencies at once. First, MGRIT almost always takes backward Euler as its coarse propagator, which is only first order and gives a contraction of about $0.1$ under FCF-relaxation; the paper switches to second-order Lobatto IIIC, the $(0,2)$ Padé approximant, which is A-stable and L-stable and damps high frequencies far more aggressively than backward Euler's $1/(1+z)$ — exactly what a coarse propagator needs. Second, and structurally, the coarse-grid correction is an inherently **sequential** sweep through $N_t$ coarse levels and becomes the serial bottleneck as the processor count grows; the paper introduces a head-tail coupling condition $\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0$, turning the coarse correction into a diagonalisable all-at-once system

  $$
  \bigl(C_\alpha\otimes I_x-I_t\otimes\Delta T A\bigr)\boldsymbol U^{k+1}=\boldsymbol g^k,
  $$

  where $C_\alpha$ is an $\alpha$-circulant matrix diagonalised by FFT. The two fixes are compatible: a better coarse propagator normally worsens the sequential bottleneck, and parallelising the coarse correction cancels that.

- **46 (parabolic PDE-constrained optimisation)** applies diagonalisation to the forward-backward optimality system of an optimal control problem, where the temporal coupling runs in both directions.
- **53 (parallel implementation of two-stage SDIRK methods)** applies diagonalisation across the stages within a single time step rather than across time levels.

### All-at-once preconditioners and spectra

- **59 (a well-conditioned direct PinT algorithm)** returns the all-at-once solution with no outer iteration, focusing on keeping the conditioning of the matrices involved under control — which is exactly where the geometric-step diagonalisation route pays its price.
- **65 (a uniform spectral analysis for a preconditioned all-at-once system)** gives one spectral analysis covering both first- and second-order evolutionary problems. The spectral distribution governs Krylov convergence, so this is the step that turns "here is a preconditioner" into "here is why it works".
- **71 (a PinT preconditioner for forward-backward evolutionary equations)** handles the forward-backward coupled systems that arise in optimal control.
- **84 (a PinT preconditioner for time spectral methods)** deals with time spectral discretisations, whose all-at-once temporal blocks are **dense** rather than bidiagonal, which changes the premises of preconditioner design.
- **85 (the _Acta Numerica_ survey)** organises this work and the wider literature into two classes: methods that remain effective for propagative problems, and methods designed primarily for dissipative problems.

> [!note] Coverage status
> The formulas of papers 12, 39, 59 and 65 have been checked equation by equation. For papers 30, 31, 46, 53 and 84 the bibliographic record and abstract are verified, but the specific constants in their theorems could not be checked. The abstracts and full texts of papers 20 and 77 are unobtainable through public channels — SIAM does not deposit abstracts with Crossref and the other aggregators return nulls for those DOIs — and paper 71's abstract is rendered with all inline mathematics deleted in every public source, so this site reports no theorems or constants for those three. The complete section-by-section reading of paper 85 is in the [[en/computational-mathematics/knowledge-notes/time-parallelization/index|time-parallelization topic]].

## One trade-off that runs through everything

The price of the diagonalisation route can be stated exactly. Geometrically increasing steps $\tau_n=\tau_1\gamma^{n-1}$ with $\gamma>1$ make the eigenvalues of the temporal matrix distinct and hence diagonalisable; but the closer $\gamma$ is to $1$, the worse the conditioning of the eigenvector matrix $V$ and the more roundoff is amplified, while the larger $\gamma$ is, the coarser the late steps and the larger the truncation error. This is a genuine dilemma, and it explains why the second route exists: with uniform steps, an $\alpha$-circulant approximation of the temporal matrix is diagonalised by FFT with conditioning controlled by $\alpha$, at the price of introducing an approximation and therefore needing an outer Krylov iteration. **The first is a direct solver, the second a preconditioner**, and that is the division of labour between the two groups of papers here.

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]; per-paper references appear at the end of each close-reading page. The complete reading of the survey itself is in [[en/computational-mathematics/knowledge-notes/time-parallelization/index|Time Parallelization for Hyperbolic and Parabolic Problems]].
