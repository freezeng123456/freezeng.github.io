---
title: Low-Rank Manifolds and Variable Separation
description: Papers 15, 17 and 83 - evolving basis and coefficients together, with the price written on the smallest singular value
lang: en
translation: computational-mathematics/paper-notes/spectral-and-reduced-order/low-rank-and-variable-separation
tags:
  - paper-notes
  - reduced-order-models
  - low-rank-approximation
---

> [!note] Coverage of this page
> Papers **15** (_SIAM J. Sci. Comput._ 37(2), 2015), **17** (_J. Comput. Phys._ 303, 2015) and **83** (_SIAM J. Sci. Comput._ 47(3), 2025). The journal text of paper 17 requires a subscription and has no preprint, so that section separates what the paper itself confirms from the general framework it adapts; the close reading of paper 83 is still to be written.

![Evolve a random solution directly on a low-rank manifold](assets/diagrams/tao-zhou-papers/en/low-rank-dynamics.svg)

## 15: an error theory for the dynamically orthogonal approximation

### Why the basis has to evolve

Time-dependent random PDEs are usually expanded in a **fixed** basis: generalised polynomial chaos in the random variables, proper orthogonal decomposition modes in space. The difficulty is that the solution's dependence on the random parameters changes substantially as time advances, so a fixed basis needs a steadily growing number of terms to hold accuracy; adaptive and greedy repairs (time-dependent polynomial chaos, generalised proper orthogonal decomposition) only partly help.

The dynamically orthogonal approximation of Sapsis and Lermusiaux evolves the spatial basis and the stochastic coefficients **together**, which sidesteps the growth but had essentially no error theory: existing results covered the closely related multiconfiguration time-dependent Hartree method and the dynamical low-rank matrix setting, not random PDEs. The paper's first goal is to establish a precise correspondence between the two so the existing theory can be imported.

### The ansatz and the gauge conditions

On a bounded $D\subset\mathbb R^d$ with $1\le d\le3$ and a complete probability space, for

$$
\frac{\partial u(x,t,\omega)}{\partial t}=\mathcal L\bigl(u(x,t,\omega),\omega\bigr),
$$

the ansatz is

$$
u_S(x,t,\omega)=\bar u_S(x,t)+\sum_{i=1}^{S}U_i(x,t)\,Y_i(t,\omega).
$$

Uniqueness of the representation is fixed by three conditions:

$$
\mathbb E[Y_i(t,\cdot)]=0,
\qquad
\langle U_i(\cdot,t),U_j(\cdot,t)\rangle=\delta_{ij},
\qquad
\Bigl\langle \frac{\partial U_i(\cdot,t)}{\partial t},\,U_j(\cdot,t)\Bigr\rangle=0 .
$$

The third is the **gauge condition**: the time derivative of the basis must be orthogonal to the span of the basis itself, so the basis is not allowed to rotate within its own span. That is exactly what removes the redundant freedom in the decomposition.

### The evolution equations and their three distinct roles

Combining Galerkin projection with the gauge conditions gives

$$
\frac{\partial\bar u_S}{\partial t}=\mathbb E\bigl[\mathcal L(u_S)\bigr],
$$

$$
\sum_{i=1}^{S}C_{ij}(t)\,\frac{\partial U_i(x,t)}{\partial t}
=\Pi^{\perp}_{U}\,\mathbb E\bigl[\mathcal L(u_S)\,Y_j(t,\cdot)\bigr],
\qquad j=1,\dots,S,
$$

$$
\frac{\partial Y_i(t,\omega)}{\partial t}
=\bigl\langle \mathcal L^{*}(u_S(\cdot,t,\omega),\omega),\,U_i(\cdot,t)\bigr\rangle,
\qquad i=1,\dots,S,
$$

with covariance matrix $C_{ij}(t)=\mathbb E[Y_iY_j]$, centred operator $\mathcal L^{*}(u,\omega)=\mathcal L(u,\omega)-\mathbb E[\mathcal L(u)]$ and the $L^2(D)$ projector

$$
\Pi^{\perp}_{U}[v]=v-\sum_{i=1}^{S}\langle v,U_i\rangle U_i
$$

onto the complement of $\mathcal U=\mathrm{span}\{U_1,\dots,U_S\}$.

The three equations play different roles. The mean equation is the plain averaged PDE. The basis equation is driven by the **part projected out of the span** of the correlation between the operator and the stochastic coefficients, and modes couple only through $C$. The coefficient equation is an ODE in $\omega$ obtained by testing the centred operator against each basis function.

### The tangent space and where $C^{-1}$ comes from

Writing the centred part as $u_S^{*}=\mathbf U^{T}\mathbf Y$ places it on the rank-$S$ manifold $\mathcal M_S$. The tangent space consists of $\delta u_S^{*}=\mathbf U^{T}\delta\mathbf Y+\delta\mathbf U^{T}\mathbf Y$, with the decomposition unique once the analogue of the gauge condition is imposed, $\langle\delta\mathbf U,\mathbf U^{T}\rangle=\mathbf 0$ and $\mathbb E[\delta\mathbf Y]=\mathbf 0$. The orthogonal projection onto the tangent space is

$$
P_{u_S^{*}}(v)=\mathbf U^{T}\langle v^{*},\mathbf U^{T}\rangle
+\Bigl(\Pi^{\perp}_{U}\bigl\{\mathbb E[v^{*}\mathbf Y^{T}]\bigr\}\,C^{-1}\Bigr)^{T}\mathbf Y,
$$

where $C=\mathbb E[\mathbf Y\mathbf Y^{T}]$ has full rank by the definition of a rank-$S$ function. **This is where the $\rho^{-1}$ that later appears in the error bound enters.** The complementary projector factorises as $P^{\perp}_{u_S^{*}}v=\Pi^{\perp}_{U}\otimes\Pi^{\perp}_{\tilde{\mathcal Y}}v$.

In Dirac-Frenkel form the method is the Galerkin condition

$$
\mathbb E\Bigl[\Bigl\langle \frac{\partial u_S}{\partial t}-\mathcal L(u_S),\,v\Bigr\rangle\Bigr]=0
\quad \forall v=\bar v+v^{*},\ v^{*}\in T_{u_S^{*}(t)}\mathcal M_S,
$$

equivalently

$$
\frac{\partial u_S}{\partial t}=\mathbb E[\mathcal L(u_S)]
+P_{u_S^{*}(t)}\bigl(\mathcal L^{*}(u_S)\bigr).
$$

That formulation compresses "evolve basis and coefficients together" into a single sentence: **project the right-hand side onto the tangent space of the manifold.**

### Theorem 4.1: quasi-optimality and its price

The reference object is the best rank-$S$ approximation at each time, the truncated Karhunen-Loève expansion

$$
z_S(x,\omega,t)=\bar u(x,t)+\sum_{i=1}^{S}\sqrt{\mu_i(t)}\,\gamma_i(t,\omega)\,Z_i(x,t),
$$

with $\{\mu_i,Z_i\}$ the eigenpairs of the covariance operator. The theorem assumes this best approximation exists and is continuously differentiable in $(H^2(D)\cap H^1_0(D))\otimes L^2(\Omega)$, with a uniform lower bound on its smallest singular value,

$$
\sigma\bigl(z_S(t)\bigr)\ \ge\ \rho\ >\ 0,\qquad \forall t\in[0,\bar t].
$$

Then there is $0<\hat t\le\bar t$ such that the dynamically orthogonal solution started at $u_S(0)=z_S(0)$ satisfies

$$
\|u_S(t)-z_S(t)\|_0^2+a_{\min}\int_0^t |u_S(\tau)-z_S(\tau)|_1^2\,\mathrm d\tau
\ \le\ 2\alpha\,e^{2\beta(t)}\int_0^t \|z_S(\tau)-u(\tau)\|_1^2\,\mathrm d\tau,
$$

$$
\beta(t)=4\rho^{-1}\int_0^t\Bigl(4\|\mathcal L^{*}(z_S)\|_0+\|\mathcal L^{*}(u)\|_0
+\|\mathcal L^{*}(u_S)\|_0+\|\dot z_S^{*}\|_0^2\Bigr)\mathrm d\tau,
\qquad
\alpha=\max\Bigl\{\frac{a_{\max}^2}{2a_{\min}},\ 4\rho^{-1}\Bigr\}.
$$

The conclusion reads on two levels. **The good news**: the dynamically orthogonal error is bounded by the best rank-$S$ error, so the method is quasi-optimal. **The price**: the constant grows like $e^{C/\rho}$, exponentially in the reciprocal of the smallest singular value. That is the quantitative form of the curvature obstruction — the manifold $\mathcal M_S$ has curvature of order $1/\sigma_S$, so as the $S$-th singular value collapses the tangent projection becomes ill-conditioned and the bound goes vacuous.

### One concrete failure mode and one implementation detail

The paper constructs an eigenvalue-crossing example: a deterministic Laplacian with random initial data whose first two Karhunen-Loève eigenvalues cross at a time $t^{*}$. Both the best rank-1 error and the dynamically orthogonal rank-1 error can then be computed in closed form, with different branches dominating before and after the crossing. The example shows the assumption is not technical: with insufficient rank the smallest singular value really does collapse.

There is also an easily missed implementation problem. The covariance matrix $C(t^n)$ can be singular or very ill-conditioned — any system with stochastic coefficients and **deterministic** initial data starts with $C\equiv0$. Applying the Moore-Penrose pseudoinverse $C^{\dagger}$ naively zeroes the "non-active" basis functions and thereby **prevents the rank from increasing**. The paper instead reformulates the basis equation as

$$
\frac{\partial\mathbf U}{\partial t}=C^{\dagger}\,\Pi^{\perp}_{U}\,
\mathbb E\bigl[\mathbf Y\,\mathcal L(u_S)\bigr],
$$

which agrees with the original when $C$ has full rank, and in practice diagonalises the covariance at each step to decouple the system, since the flow does not preserve uncorrelatedness of the $Y_i$ even when $C(0)$ is diagonal.

## 17: confining the eigensolve to the coarsest mesh

A Karhunen-Loève expansion needs the eigenpairs of the Fredholm integral operator built from the covariance kernel: given $\mathrm{Cov}(x,y)$, find

$$
\int_D \mathrm{Cov}(x,y)\,u(y)\,\mathrm dy=\lambda\,u(x),
\qquad \|u\|_{L^2(D)}=1 .
$$

Compared with an eigenvalue problem for a differential operator there is a structural difference: the integral operator $\mathcal T$ is **compact and smoothing**, its eigenvalues accumulate at $0$ rather than at $+\infty$, and discretisation produces a **dense** matrix. So there is no sparsity to exploit, and both assembly and solution grow sharply with the mesh resolution needed to resolve the kernel.

Earlier work attacked that cost from the discretisation side — wavelet-Galerkin schemes, generalised fast multipole methods, spectral elements, Legendre-Galerkin with tensor structure, multilevel augmentation for compact integral operators — but the eigensolve itself still had to happen on the finest mesh. This paper removes that: the eigenvalue problem becomes a sequence of **integral iterations** plus eigenvalue solves on the **coarsest** mesh only, so any efficient integration scheme can be plugged in and total work is comparable to one integration step on the finest mesh.

The multilevel correction framework it adapts (Lin and Xie) has a generic step: given a coarsest space $V_H$ and a nested hierarchy $V_H\subset V_{h_2}\subset\cdots$, first solve a **source problem** on the finer space,

$$
a(\tilde u_{h_{k+1}},v)=\lambda_{h_k}\,b(u_{h_k},v),\qquad \forall v\in V_{h_{k+1}},
$$

then solve the eigenvalue problem in the small augmented space $V_{H,h_{k+1}}=V_H+\mathrm{span}\{\tilde u_{h_{k+1}}\}$,

$$
a(u_{h_{k+1}},v)=\lambda_{h_{k+1}}\,b(u_{h_{k+1}},v),\qquad \forall v\in V_{H,h_{k+1}} .
$$

The point of the second step is that the eigenproblem is only ever solved in a space of dimension $\dim V_H+1$, the coarsest space plus one direction, never on the fine mesh. For an integral operator the "source problem" in the first step is not a linear solve but **one quadrature**, an application of $\mathcal T$, which is what the abstract means by integral iterations.

> [!warning] What could be verified
> The journal text requires a subscription and there is no preprint. Confirmable: the problem setting, the keywords (uncertainty quantification, Karhunen-Loève expansion, Fredholm eigenvalue problem, multigrid finite element) and the contribution as stated in the abstract. The generic correction step above comes from the Lin-Xie framework being adapted; **the specific adaptation to the integral operator, the constants and all numerical results are unverified here** and should be checked against the published article.

## 83: dynamical variable separation for parameter-dependent systems

Paper 83 belongs to the same family as paper 15: build a separated representation that evolves in time for a parameter-dependent dynamical system rather than expanding in a fixed basis. The published title is _A Dynamical Variable-Separation Method for Parameter-Dependent Dynamical Systems_, slightly different from the homepage listing; this site records the published version.

> [!note] Close reading pending
> This page has not yet checked the paper equation by equation, so it does not report the separated form, the update formulas or the theoretical results. The usable orientation is that it shares the stance of "let the representation evolve with the solution" with paper 15, in contrast to fixed-basis reduced-order models.

## How the three relate

| No. | What evolves                              | Optimal reference              | Principal risk                          |
| --- | ----------------------------------------- | ------------------------------ | --------------------------------------- |
| 15  | spatial basis and stochastic coefficients | truncated Karhunen-Loève       | collapse of the smallest singular value |
| 17  | nothing (a one-off eigensolve)            | the exact Karhunen-Loève basis | assembling and solving a dense matrix   |
| 83  | a separated representation (pending)      | pending                        | pending                                 |

The relation between papers 17 and 15 is worth naming: the error bound in paper 15 is stated against the truncated Karhunen-Loève approximation, and paper 17 solves the problem of **computing that reference object**. Together they give the full cost structure of this route: either pay for the dense eigenproblem that produces the Karhunen-Loève basis, or let the basis evolve and accept the risk carried by the smallest singular value.

## Coverage check

| Item                                               | Paper | Status                                                             |
| -------------------------------------------------- | ----- | ------------------------------------------------------------------ |
| Ansatz and the three gauge conditions              | 15    | form, conditions, meaning of the gauge                             |
| Three evolution equations and their roles          | 15    | mean, basis, coefficient equations, projector                      |
| Tangent projection and the origin of $C^{-1}$      | 15    | tangent space, projection formula, Dirac-Frenkel form              |
| Theorem 4.1 and the $e^{C/\rho}$ constant          | 15    | assumptions, conclusion, both readings, curvature                  |
| The eigenvalue-crossing failure mode               | 15    | construction and what it demonstrates                              |
| Handling a singular covariance                     | 15    | the pseudoinverse trap, reformulation, per-step diagonalisation    |
| Structural difference of the Fredholm problem      | 17    | compactness, spectral accumulation, dense matrix                   |
| Multilevel correction step and integral iterations | 17    | the two generic steps, dimension argument, quadrature versus solve |

## Sources for this page

- E. Musharbash, F. Nobile, and T. Zhou, [_Error analysis of the dynamically orthogonal approximation of time dependent random PDEs_](https://doi.org/10.1137/140967787), SIAM J. Sci. Comput. 37(2) (2015), pp. A776-A810.
- H. Xie and T. Zhou, [_A multilevel finite element method for Fredholm integral eigenvalue problems_](https://doi.org/10.1016/j.jcp.2015.09.043), J. Comput. Phys. 303 (2015), pp. 173-184.
- L. Chen, Y. Chen, Q. Li, and T. Zhou, [_A dynamical variable-separation method for parameter-dependent dynamical systems_](https://doi.org/10.1137/24M168427X), SIAM J. Sci. Comput. 47(3) (2025), pp. A1783-A1808.
