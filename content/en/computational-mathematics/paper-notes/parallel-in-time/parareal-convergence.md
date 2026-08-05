---
title: Convergence Analysis for Parareal
description: Papers 12, 20, 30 and 77 - how the contraction factor depends on fine-propagator stability and on the fine grid
lang: en
translation: computational-mathematics/paper-notes/parallel-in-time/parareal-convergence
tags:
  - paper-notes
  - parallel-in-time
  - convergence-analysis
---

> [!note] Coverage of this page
> Papers **12** (_SIAM J. Sci. Comput._ 37(2), 2015), **20** (_J. Comput. Phys._ 329, 2017), **30** (_J. Comput. Phys._ 358, 2018) and **77** (_SIAM J. Numer. Anal._ 62(5), 2024). None has a preprint, and the text and abstract of papers 20 and 77 could not be obtained in full through public channels, so those sections state what could be verified.

## The standard analytical machinery

Consider the symmetric positive definite system $\boldsymbol u'(t)+A\boldsymbol u(t)=g(t)$ obtained by semi-discretising a parabolic PDE. The parareal iteration is

$$
\boldsymbol u_{n+1}^{k+1}=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^{k})
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k}),
$$

with $J=\Delta T/\Delta t\ge2$ fine steps per coarse interval. If $A=V_ADV_A^{-1}$ and both propagators are one-step methods with stability functions $R_g$ and $R_f$, the error obeys $\boldsymbol\xi^{k+1}(z)=M(z)\boldsymbol\xi^k(z)$ with

$$
M(z)=M_g^{-1}(z)\bigl[M_g(z)-M_f(z)\bigr]=I_t-M_g^{-1}(z)M_f(z),
$$

where $M_g(z)$ and $M_f(z)$ are lower-triangular Toeplitz matrices with unit diagonal and subdiagonal entries $-R_g(z)$ and $-R_f^{J}(z/J)$, and $z=\Delta T\lambda$. Factoring out the scalar gives

$$
M(z)=\bigl[R_f^J(z/J)-R_g(z)\bigr]\,\widetilde M\bigl(R_g(z)\bigr),
$$

with $\widetilde M(\beta)$ strictly lower-triangular Toeplitz carrying $1,\beta,\beta^2,\dots$ down its subdiagonals. The linear and superlinear factors are therefore

$$
\varrho_l(J,z)=\frac{\bigl|R_g(z)-R_f^{J}(z/J)\bigr|}{1-|R_g(z)|},
\qquad
\varrho_s(J,z,N_t,k)=\frac{\bigl|R_g(z)-R_f^J(z/J)\bigr|^k}{k!}\prod_{j=1}^k(N_t-j),
$$

the second vanishing at $k=N_t$, which is termination in exact arithmetic.

**One structural assumption is built into this machinery, and paper 77 is about removing it.** The action of $\mathcal F$ over a coarse interval is written as the $J$-th **power** of a single scalar stability function at a single argument, $R_f^{J}(z/J)$, and that step requires the $J$ fine steps inside the interval to be equal.

## 12: what happens when the fine propagator is only A-stable

### What the theory covered and what it did not

The best-understood case is parareal-Euler: backward Euler for both propagators, for which Mathew, Sarkis and Schaerer proved a mesh-independent, $T$-independent contraction factor of roughly $0.298$ valid for every $J\ge2$.

What was not understood is a higher-order $\mathcal F$ that is **only A-stable rather than L-stable**, most importantly the trapezoidal rule and the fourth-order two-stage Gauss Runge-Kutta method. With the convention $\boldsymbol u'+A\boldsymbol u=g$ so that $z=\Delta T\lambda$ with $\lambda>0$, the relevant stability functions are

$$
\text{backward Euler:}\ R(z)=\frac{1}{1+z},
\qquad
\text{trapezoidal:}\ R(z)=\frac{2-z}{2+z},
\qquad
\text{fourth-order Gauss:}\ R(z)=\frac{12-6z+z^2}{12+6z+z^2}.
$$

Backward Euler is L-stable with $R(\infty)=0$; the trapezoidal rule is A-stable but not L-stable with $R(\infty)=-1$; the fourth-order Gauss method is A-stable and symplectic with $R(\infty)=+1$.

For the latter two, $|R_f(z)|\to1$ as $z\to\infty$, so the argument producing the $\approx0.3$ constant breaks down: the numerator of $\varrho_l$ no longer decays with frequency, and the contraction factor can approach $1$ when $\lambda_{\max}$ is large. Earlier analyses assumed either the exact propagator $\exp(\Delta TA)$ or L-stability, and therefore said nothing about the practically important A-stable-only case, especially at small $J$.

This paper's contribution is a critical-ratio framework for that case: how large $J$ must be to push the high-frequency behaviour back under control. The resulting statements appear in Chapter 4 of the [[en/computational-mathematics/knowledge-notes/time-parallelization/index|time-parallelization survey reading]] as equation (4.8), with the finite-spectrum extension $J_{\min}=O(\log^2 z_{\max})$.

## 20 and 30: two kinds of "fractional", two different difficulties

These two papers are often cited together, but the difficulties they treat are **structurally different**.

### 20: space-fractional, where the coarse propagator becomes the bottleneck

The setting is a time-dependent PDE containing a fractional Laplacian $(-\Delta)^{\alpha}$, semi-discretised in space into a stiff system $\boldsymbol u'+A\boldsymbol u=g$. Here $A$ inherits the nonlocality of the operator and is therefore **dense** (Toeplitz-structured in one dimension) with a very large $\lambda_{\max}$. Two difficulties compound: each implicit solve with $A$ is far more expensive than in the classical Laplacian case, so the **coarse** propagator $\mathcal G$ — executed sequentially and therefore directly on the critical path — dominates the runtime; and the huge $\lambda_{\max}$ is exactly the regime in which the critical ratio of paper 12 grows like $O(\log^2(\Delta T\lambda_{\max}))$. The word "fast" in the title refers to accelerating the parareal **iteration itself** rather than merely bounding its contraction factor.

> [!warning] What could be verified
> The abstract could not be retrieved from Crossref, OpenAlex, Semantic Scholar, NASA ADS or the publisher page, and the full text is unavailable. What is confirmable: third-party literature records that the paper analyses parareal convergence for space-fractional diffusion **on constant time steps**. Its specific fractional discretisation, theorems, contraction factors and numerical examples are unverified here, so no numeric constant is reported.

### 30: time-fractional, where the load becomes unbalanced

The difficulty for time-fractional equations is structural: the Caputo or Riemann-Liouville derivative at time $t$ depends on the entire solution history on $[0,t]$. But the whole premise of parareal is that fine propagation on $[T_n,T_{n+1}]$ needs only the local initial value $\boldsymbol u_n^k$. Once the operator is nonlocal in time, advancing on the $n$-th subinterval also needs the history over $[0,T_n]$, so a naive extension gives processor $n$ a workload proportional to $n$: **the computational time is unbalanced across processes**, with the last processor doing $N_t$ times the work of the first, destroying parallel efficiency even when the iteration converges.

The remedy is to **localise** the fractional operator with two local time integrators: introduce auxiliary variables whose evolution is governed by local (ODE) dynamics, so the nonlocal history is encoded in a finite set of extra local unknowns rather than in an ever-growing convolution sum. Once the problem is written as an augmented **local** system in (solution variable, auxiliary variables), each fine-propagator call costs the same on every subinterval and the load is balanced.

The genuinely new algorithmic ingredient is a **mixed coarse-grid correction** in which the auxiliary variables and the solution variable are corrected **separately**, rather than applying the single parareal update to the whole augmented state vector. The paper reports that the algorithm admits a robust rate of convergence, and in this series "robust" consistently means independent of the eigenvalues of the spatial operator, of the coarsening ratio $J$ and of the number of coarse intervals $N_t$.

> [!note] What could be verified
> The problem setting, the localisation idea and the **existence** of the mixed correction are confirmable from the abstract. The names of the two local time integrators, the exact mixed-correction formula, the explicit convergence factor and its hypotheses, and whether the load-balancing claim is quantified, are all unverified here.

## 77: removing the equal-fine-step assumption

### Why this is not a routine extension

The paper studies parareal convergence with a **uniform coarse grid** and an **arbitrarily distributed, nonuniform fine grid**. The phrase "arbitrarily distributed" carries the substance: the fine grid inside a coarse interval is not assumed uniform, not assumed graded, and not even assumed the same from one coarse interval to the next.

Returning to the machinery above shows where the difficulty sits. Every classical linear convergence result — Gander and Vandewalle's bounds, the constant near $0.3$, the critical-ratio formulas of paper 12, the MGRIT factor of paper 39 — rests on the identity

$$
\mathcal F(T_n,T_{n+1},\cdot)\ \longleftrightarrow\ R_f^{J}\!\Bigl(\frac{z}{J}\Bigr),
$$

that the fine propagator over a coarse interval is the $J$-th power of a single scalar stability function at a single argument. With a nonuniform fine grid one instead gets a **product of distinct factors**,

$$
\mathcal F(T_n,T_{n+1},\cdot)\ \longleftrightarrow\ \prod_{i=1}^{J_n}R_f(\theta_{n,i}z),
\qquad
\theta_{n,i}=\frac{\Delta t_{n,i}}{\Delta T},\quad \sum_i\theta_{n,i}=1 .
$$

Two things then fail. First, $\varrho_l(J,z)$ is no longer even defined, because there is no single $J$. Second, and more seriously, if the fine grid **differs across coarse intervals** the subdiagonal entries of $M_f(z)$ vary from row to row, so $M_f(z)$ is **no longer Toeplitz**, the factorisation $M(z)=[R_f^J(z/J)-R_g(z)]\widetilde M(R_g(z))$ fails, and Gander and Vandewalle's bound on $\|\widetilde M^k\|_\infty$ — the engine of every convergence-factor result in this literature — cannot be applied. Recovering a robust, mesh-independent contraction factor in this setting is a genuinely new analysis.

The reason one wants nonuniform fine grids is clear: adaptive time stepping, graded meshes near $t=0$ for nonsmooth or incompatible initial data, and local refinement around fast transients all produce them. Restricting parareal theory to uniform fine grids excludes essentially all adaptive practice.

> [!warning] What could be verified
> The abstract could not be obtained: SIAM does not deposit abstracts in Crossref, and OpenAlex, Semantic Scholar and NASA ADS all return nothing for this DOI. Confirmable: the problem setting above, and the paper's position in the literature — its reference list is the classical parareal-convergence reading list (paper 12, Gander and Hairer, Bal, Gander and Vandewalle, Southworth and others), which identifies it as a pure analysis contribution rather than a new algorithm. This site does not report its theorems, convergence factors or thresholds.
>
> The calibration point is the constant near $0.3$ for a uniform fine grid with an L-stable $\mathcal F$ and backward Euler coarse propagator. Whether this paper attains a comparable constant on nonuniform fine grids is unverified here.

## One tension worth naming

Nonuniform time steps occupy opposite positions in the two routes. In the parareal branch they are an **obstacle to the theory**: they destroy the Toeplitz structure and with it every convergence-factor argument. In the direct diagonalisation branch (see [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|the diagonalisation technique]]) nonuniform steps are a **requirement**: Maday and Rønquist need the $\Delta t_n$ to be distinct for the temporal matrix to be diagonalisable. The same modelling freedom is a nuisance on one route and a prerequisite on the other.

## Coverage check

| Item                                                | Paper  | Status                                                           |
| --------------------------------------------------- | ------ | ---------------------------------------------------------------- |
| Parareal iteration and error matrix                 | 12     | iteration, $M(z)$, factorisation, both factors                   |
| Stability functions of four integrators             | 12     | backward Euler, trapezoidal, Gauss, and their $R(\infty)$        |
| Why the A-stable-only case escapes old theory       | 12     | non-decaying numerator, need for a critical ratio                |
| Space-fractional bottleneck in the coarse solve     | 20     | dense $A$, large $\lambda_{\max}$ (limited verification)         |
| Historical effect and unbalanced load               | 30     | broken premise, workload proportional to $n$                     |
| Localisation and mixed coarse correction            | 30     | auxiliary variables, augmented local system, separate correction |
| Nonuniform fine grid destroying Toeplitz structure  | 77     | power to product, both failures, why it is new                   |
| Opposite role of nonuniform steps in the two routes | 12, 77 | obstacle versus prerequisite                                     |

## Sources for this page

- S. Wu and T. Zhou, [_Convergence analysis for three parareal solvers_](https://doi.org/10.1137/140970756), SIAM J. Sci. Comput. 37(2) (2015), pp. A970-A992.
- S. Wu and T. Zhou, [_Fast parareal iterations for fractional diffusion equations_](https://doi.org/10.1016/j.jcp.2016.10.046), J. Comput. Phys. 329 (2017), pp. 210-226.
- S. Wu and T. Zhou, [_Parareal algorithms with local time-integrators for time fractional differential equations_](https://doi.org/10.1016/j.jcp.2017.12.029), J. Comput. Phys. 358 (2018), pp. 135-149.
- S.-L. Wu and T. Zhou, [_Convergence analysis of the parareal algorithm with nonuniform fine time grid_](https://doi.org/10.1137/23M1592481), SIAM J. Numer. Anal. 62(5) (2024), pp. 2308-2330.
