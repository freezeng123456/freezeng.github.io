---
title: Stability Theory for Discretisations
description: Papers 47 and 63 - a mean-square Lax equivalence theorem, and the first fully discrete error analysis
lang: en
translation: computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes
tags:
  - paper-notes
  - stochastic-differential-equations
  - stability-theory
---

> [!note] Coverage of this page
> Papers **47** (_SIAM J. Numer. Anal._ 58(4), 2020) and **63** (_SIAM J. Numer. Anal._ 60(4), 2022). Neither has a preprint and both texts require a subscription, so this page separates what the abstract confirms from what is reconstructed from context.

## 47: replacing scheme-by-scheme proofs with one equivalence theorem

### The problem

By 2020 the group had produced a proliferation of temporal discretisations for forward-backward systems: backward Euler, the theta-scheme family, the Zhao-Zhang-Ju multistep scheme, the multistep family of [[en/computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|paper 8]], and deferred correction. Each came with its own ad hoc convergence proof and its own stability condition, such as the empirically observed $k\le6$ root-condition window in paper 8. What was missing was a **single framework** in which stability is defined once and convergence follows abstractly.

### The main theorem

The paper supplies a general discretisation family and proves a **mean-square version of the Lax equivalence theorem**:

> For forward-backward stochastic differential equations, a consistent discretisation scheme is convergent if and only if it is stable.

This is the exact counterpart of the classical Lax equivalence theorem for linear evolution equations, with the norm now mean-square ($L^2(\Omega)$) and the problem nonlinear (Lipschitz generator). The abstract also states that applications of the analysis to existing schemes are discussed, that is, the abstract theorem is instantiated to recover or strengthen the known convergence results for backward Euler, theta-schemes and multistep schemes.

**The theorem gives the empirical window of paper 8 a theoretical home**: the root condition there is exactly the stability half, consistency comes from the construction, and together they give convergence.

### The two families the framework covers

One of the two unified families is the theta-scheme:

$$
Y^\pi_{t_k}=\mathbb E\bigl[Y^\pi_{t_{k+1}}\mid\mathcal F_{t_k}\bigr]
+\theta_1\Delta_n f\bigl(t_k,Y^\pi_{t_k},Z^\pi_{t_k}\bigr)
+(1-\theta_1)\Delta_n\,\mathbb E\bigl[f(t_{k+1},Y^\pi_{t_{k+1}},Z^\pi_{t_{k+1}})\mid\mathcal F_{t_k}\bigr],
$$

$$
Z^\pi_{t_k}=-\frac{1-\theta_2}{\theta_2}\mathbb E\bigl[Z^\pi_{t_{k+1}}\mid\mathcal F_{t_k}\bigr]
+\frac{1}{\theta_2\Delta_n}\mathbb E\bigl[Y^\pi_{t_{k+1}}\Delta W_k\mid\mathcal F_{t_k}\bigr]
+\frac{1-\theta_2}{\theta_2}\mathbb E\bigl[f(t_{k+1},Y^\pi_{t_{k+1}},Z^\pi_{t_{k+1}})\Delta W_k\mid\mathcal F_{t_k}\bigr],
$$

with $\theta_1\in[0,1]$ and $\theta_2\in(0,1]$. This scheme converges with order two when $\theta_1=\theta_2=1/2$, the Crank-Nicolson member, and order one otherwise. The other family is the multistep scheme of paper 8.

> [!note] What could be verified
> The **content** of the main theorem (consistency plus stability is equivalent to convergence) and the coverage claim (backward Euler, theta-schemes, various multistep schemes) are confirmable from the abstract. The exact form of the general discretisation family, the precise norms in the definitions of stability and consistency, and whether the framework covers coupled as well as decoupled systems, are unverified here. The theta-scheme above is transcribed in its standard form from an independent third-party survey.

## 63: the first fully discrete error analysis

### How the spatial error had been sidestepped

Every probabilistic BSDE scheme has two error sources: the **temporal** discretisation of the backward equation and the **spatial** approximation of the conditional expectations $\mathbb E_{t_n}[\cdot]$, which in the Markovian case are $d$-dimensional Gaussian integrals of functions known only at grid points.

The group's earlier pipeline (papers 8, 19, 25, 41) computed those integrals by **Gauss-Hermite quadrature**, whose nodes $x_n+\sqrt{\Delta t}\lambda_j$ do **not** fall on the spatial grid, so a local polynomial interpolation is unavoidable at every node at every time step. That costs three things: time; a ceiling on spatial accuracy at the interpolation order; and a **fully discrete** (time plus space) error analysis that becomes extremely awkward. The consequence is that before this paper essentially all rigorous BSDE error analyses were **semi-discrete**, bounding the temporal error while assuming exact conditional expectations.

The paper's own framing is that this seems to be the **first attempt at analysing fully discrete schemes for BSDEs**, with second-order convergence in time and exponential convergence in space.

### The two components

The temporal part uses the theta-scheme family above, with $\theta_1=\theta_2=1/2$ giving second order, matching the claimed second-order convergence in time.

The spatial part uses **Sinc approximation**. The key property of a Sinc quadrature rule is exponential convergence for a suitable class of functions, with nodes distributed by a $\mathrm{sinh}$-type transformation so the rule can be arranged compatibly with the grid — which is how the interpolation step is avoided. The paper's keywords list "backward stochastic differential equations; error estimates; conditional mathematical expectation; Sinc-theta schemes", consistent with this.

**Where this paper sits in the thread**: paper 47 unifies the stability theory in the temporal direction and paper 63 supplies the spatial direction, so the error analysis covers a fully discrete scheme for the first time. Both use the same "stability first" architecture.

> [!note] What could be verified
> The first-fully-discrete claim, the second-order-in-time and exponential-in-space rates, and the use of theta-schemes in time with Sinc approximation for the conditional expectations are all confirmable from the abstract. The theta-scheme above is transcribed from an independent source by the same group, while the concrete Sinc quadrature construction and its error constants are unverified here.

## How the two relate

| No. | Error sources covered | Main result                                                | Effect on the other papers                                                          |
| --- | --------------------- | ---------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| 47  | temporal              | mean-square Lax equivalence                                | gives paper 8's empirical window a theory; makes paper 68's reverse design possible |
| 63  | temporal and spatial  | first fully discrete analysis (order two plus exponential) | a second route of spatial improvement alongside paper 25                            |

One general judgement: **when a family of methods has scheme-by-scheme convergence proofs, it is worth stopping to define stability properly.** After paper 47, "design a convergent scheme" becomes "design a stable scheme", and the latter is an objective one can optimise — which is exactly what [[en/computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|paper 68]] does.

## Sources for this page

- J. Yang, W. Zhao, and T. Zhou, [_A unified probabilistic discretization scheme for FBSDEs: stability, consistency, and convergence analysis_](https://doi.org/10.1137/19M1260177), SIAM J. Numer. Anal. 58(4) (2020), pp. 2351-2375.
- X. Wang, W. Zhao, and T. Zhou, [_Sinc-theta schemes for backward stochastic differential equations_](https://doi.org/10.1137/21M1444679), SIAM J. Numer. Anal. 60(4) (2022), pp. 1799-1823.
