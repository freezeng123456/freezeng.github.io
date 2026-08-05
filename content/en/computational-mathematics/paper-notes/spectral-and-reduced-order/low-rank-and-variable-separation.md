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
> Papers **15** (_SIAM J. Sci. Comput._ 37(2), 2015), **17** (_J. Comput. Phys._ 303, 2015) and **83** (_SIAM J. Sci. Comput._ 47(3), 2025). The journal text of paper 17 requires a subscription and has no preprint, so that section separates what the paper itself confirms from the general framework it adapts; papers 15 and 83 are checked equation by equation against publicly available full texts.

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

## 83: enrich the basis one term at a time so every step decouples

### Three reductions, each freezing the wrong thing

The paper states its position sharply. Generalised polynomial chaos writes $u\approx\sum_i\zeta_i(\xi)g_i(x,t)$, **freezing the parametric coefficients in time** while the spatial modes move, and needs substantial prior information from the full model. Proper orthogonal decomposition writes $u\approx\sum_i\zeta_i(t;\xi)g_i(x)$, **freezing the spatial modes in both time and parameter**.

Whether proper orthogonal decomposition succeeds depends on the Kolmogorov widths of the solution manifold $\mathcal M=\{u(\cdot,t;\xi)\}$ decaying fast, and that **fails for first-order linear transport and hyperbolic problems** — the Kolmogorov barrier. The dynamically orthogonal method analysed in paper 15 breaks the barrier by letting both factors move, but it and its variants (DyBO, dual-DO, DBO) all solve a **coupled** system for all $N$ modes at once, and that is precisely where the covariance inverse and its conditioning problems come from.

This paper's angle is to keep the doubly time-dependent ansatz but build the modes **one at a time by greedy enrichment**, so each step decouples into two scalar-coefficient subproblems.

### The model and the affine assumption

$$
\frac{\partial u}{\partial t}(x,t;\xi)=F\bigl(u(x,t;\xi);\xi\bigr),
\qquad u(x,0;\xi)=\mu(x;\xi),
$$

with the structural assumption $F(u;\xi)=\mathcal C(\xi)+\mathcal A(u;\xi)+\mathcal H(u;\xi)$, where $\mathcal A$ is linear and $\mathcal H$ nonlinear, and **affine parameter dependence**:

$$
\mathcal C(\xi)=\sum_{i}\kappa^i_C(\xi)\mathcal C^i,
\qquad
\mathcal A(u;\xi)=\sum_{i}\kappa^i_A(\xi)\mathcal A^i(u),
\qquad
\mathcal H(u;\xi)=\sum_{i}\kappa^i_H(\xi)\mathcal H^i(u),
$$

together with a separable initial condition $\mu(x;\xi)=\sum_i p^i(\xi)q^i(x)$. **Affinity is exactly what makes the offline-online split possible**; the paper notes that when it fails, a variable-separation step can first produce an affine approximation with negligible loss.

### The ansatz: no mean field and no orthogonality gauge

$$
u(x,t;\xi)\approx u_N(x,t;\xi):=\sum_{i=1}^{N}\zeta_i(t;\xi)\,g_i(x,t),
$$

with $\{\zeta_i\}$ parameter-dependent, $\{g_i\}$ parameter-independent, and **both time-dependent**. Compare the dynamically orthogonal ansatz $u\approx\bar u(x,t)+\sum_i\zeta_i(t;\xi)g_i(x,t)$, which keeps a statistical mean field and requires $\{g_i(\cdot,t)\}$ orthonormal at every $t$. **This paper drops both the mean field and the orthonormality gauge**, and uniqueness comes instead from the sequential greedy construction.

### The greedy rule and the two decoupled subproblems

The first step picks $\xi_1$ arbitrarily and takes $g_1(x,t)$ to be the **full solution** at $\xi=\xi_1$, after which testing with $g_1$ gives a parameter-dependent ordinary differential equation for $\zeta_1(t;\xi)$. At step $k\ge2$, with the error $e:=u-u_{k-1}$, the selection rule is

$$
\xi_k\in\arg\max_{\xi\in\Xi}\triangle_k(\xi),
$$

where $\triangle_k$ is either $\|e\|_{L^2([0,T];V)}$ itself, if affordable, or the a posteriori bound below; the loop stops when $\triangle_k(\xi_k)<\varepsilon$ and otherwise removes $\xi_k$ from the candidate set. Rewriting the equation in terms of the error,

$$
\Bigl\langle \frac{\partial(e+u_{k-1})}{\partial t},v\Bigr\rangle
=\bigl\langle F\bigl((e+u_{k-1});\xi\bigr),v\bigr\rangle,
\qquad \forall v\in V,
$$

the two subproblems are: $g_k(x,t)$ solves this at the **single** parameter $\xi=\xi_k$, a parameter-independent partial differential equation; and setting $e=g_k\zeta_k$ with $v=g_k$ gives a scalar ordinary differential equation for $\zeta_k(t;\xi)$.

### Initial conditions: recovering the gauge implicitly at $t=0$

The paper calls the initial conditions "one of the most essential ingredients", and constructs them by $L^2$ matching against the current mode. Writing $g_{k,0}(x)=g_k(x,0)$, we have $g_1(x,0)=\mu(x;\xi_1)$ and

$$
\zeta_{1,0}(\xi)=\sum_{i}\frac{\langle q^i,g_{1,0}\rangle}{\langle g_{1,0},g_{1,0}\rangle}p^i(\xi);
$$

for $k\ge2$ the initial error $e_0(x;\xi)=\mu(x;\xi)-\sum_{j<k}g_{j,0}(x)\zeta_{j,0}(\xi)$ gives $g_{k,0}=e_0(\cdot;\xi_k)$ and

$$
\zeta_{k,0}(\xi)=\sum_{i}\frac{\langle q^i,g_{k,0}\rangle}{\langle g_{k,0},g_{k,0}\rangle}p^i(\xi)
-\sum_{j=1}^{k-1}\frac{\langle g_{j,0},g_{k,0}\rangle}{\langle g_{k,0},g_{k,0}\rangle}\zeta_{j,0}(\xi).
$$

**The second sum is a Gram-Schmidt-like correction.** It is where this method recovers, implicitly and only at $t=0$, part of what the dynamically orthogonal gauge condition supplies explicitly. One edge case: if $\mu\equiv0$ then $g_{k,0}=0$ and $\zeta_{k,0}=0$ at every step.

### The linear case: the new mode is driven by the previous residual

With $F(u;\xi)=\mathcal C(\xi)+\mathcal A(u;\xi)$ the error equation becomes

$$
\Bigl\langle\frac{\partial e}{\partial t},v\Bigr\rangle-\bigl\langle\mathcal A(e;\xi),v\bigr\rangle
=\bigl\langle r_k,v\bigr\rangle,
\qquad
r_k:=\mathcal A(u_{k-1};\xi)-\frac{\partial u_{k-1}}{\partial t}+\mathcal C(\xi),
$$

so **the new spatial mode is driven by the residual of the previous approximation**. The coefficient equation is the scalar ordinary differential equation

$$
\bigl\langle (g_k)_t,g_k\bigr\rangle\zeta_k
+\bigl\langle g_k,g_k\bigr\rangle(\zeta_k)_t
-\bigl\langle\mathcal A(g_k;\xi),g_k\bigr\rangle\zeta_k
=\bigl\langle r_k,g_k\bigr\rangle .
$$

Note the term $\langle(g_k)_t,g_k\rangle$: in the dynamically orthogonal method it is **forced to vanish** by the gauge condition, whereas here it is simply carried along. That reduces the difference between the two routes to one concrete term.

### The decisive point: the divisor is a scalar, not a matrix

Splitting $[0,T]$ into $N_t$ steps, using first-order differences and backward Euler, and substituting the affine expansions together with the separated form of $u_{k-1}$, gives a closed-form recursion

$$
\zeta_{k,n+1}(\xi)=\frac{c_{n+1}\,\zeta_{k,n}(\xi)+s_{n+1}(\xi)}{l_{n+1}(\xi)},
\qquad n=0,\dots,N_t-1,
$$

$$
c_{n+1}=\frac{\langle g_{k,n+1},g_{k,n+1}\rangle}{\tau},
\qquad
l_{n+1}(\xi)=2c_{n+1}-\frac{\langle g_{k,n},g_{k,n+1}\rangle}{\tau}
-\sum_{i}\kappa^i_A(\xi)\bigl\langle\mathcal A^i(g_{k,n+1}),g_{k,n+1}\bigr\rangle,
$$

with $s_{n+1}(\xi)$ assembled from the affine terms and the $k-1$ modes already computed. **This is the paper's answer to the covariance inversion of the dynamically orthogonal method: $l_{n+1}(\xi)$ is a scalar rather than a matrix, so neither $C^{-1}$ nor $C^{\dagger}$ appears anywhere.**

The offline stage stores only parameter-independent scalar inner products such as $\langle g_{k,n+1},g_{k,n+1}\rangle$ and $\langle\mathcal A^i(g_{j,n+1}),g_{k,n+1}\rangle$. Online, a new parameter $\bar\xi$ needs only a cheap scalar ordinary differential equation whose coefficients are affine in the parameter, so **the online cost is independent of the spatial discretisation** of the original problem.

In the nonlinear case (Burgers) the residual equation contains $\langle e_{n+1}\partial_xe_{n+1},v\rangle$, which the paper approximates semi-implicitly by $e_n\partial_xe_{n+1}$ to control cost, and states that the same procedure carries over to Allen-Cahn.

### The a posteriori bound, and the absence of a convergence theorem

The bound rests on the **local logarithmic Lipschitz constant** of $F$ at $u$,

$$
L_V[F](u):=\sup_{v\ne u}\frac{\langle v-u,\;F(v;\xi)-F(u;\xi)\rangle}{\|v-u\|_V^2}.
$$

Testing the error equation $\partial_te=F(u;\xi)-F(u_{k-1};\xi)+r_k$ with $e$ and applying a comparison lemma gives $\|e\|_V\le\delta_k(t;\xi)$ with

$$
\delta_k(t;\xi)=\int_0^{t}\alpha(s;\xi)e^{\int_s^{t}\beta(\tau;\xi)\mathrm d\tau}\mathrm ds
+e^{\int_0^{t}\beta(\tau;\xi)\mathrm d\tau}\|e(\cdot,0;\xi)\|_V,
$$

where $\alpha=\|r_k\|_V$ and $\beta=L_V[F](u_{k-1})$, and the selection criterion takes $\triangle_k(\xi)=\int_0^T\delta_k^2\mathrm dt$. **For dissipative $F$ the logarithmic Lipschitz constant can be negative**, in which case the exponential factors damp rather than amplify — which is what makes the estimator usable over long horizons.

> [!warning] No convergence theorem
> The conclusion explicitly lists "a rigorous convergence analysis of the proposed method under reasonable assumptions" as future work. The rigorous content is the a posteriori bound above together with the exact algebraic recursion; the claim of reduced complexity and improved efficiency over existing low-rank separation techniques is the abstract's, supported numerically rather than by theorem. The paper also states a limitation of its own: both storage and online computation still depend on the size of the space-time discretisation, and reducing that is left open.

### The two most informative numerical comparisons

**A head-to-head against a time-dependent reduced-basis method** (one-dimensional reaction-diffusion with a parameter-dependent boundary condition, $\xi\in[1,3]^4$, $M=10^3$ test samples) yields a frank trade-off. At $N=4$ the competitor reaches error $9.19\times10^{-6}$ against this method's $1.46\times10^{-4}$, so **the competitor's error falls faster**; the paper's own explanation is that the competitor updates all parametric coefficients at every iteration while this method never modifies the time-parameter basis functions computed at earlier steps. But **beyond $N=7$ the competitor's error rises quickly**, attributed to ill-conditioning of its linear systems for the parametric coefficients, while this method's online time is more than an order of magnitude smaller ($7.20\times10^{-3}$ s per sample at $N=7$ against $1.19\times10^{-1}$ s, with full-order finite elements plus backward Euler at $3.75$ s per sample). The paper's summary: the competitor gives higher accuracy under certain conditions, and this method does better on efficiency and error stability.

**A comparison against the static-basis variable-separation method** (two-dimensional heat equation) is the cleanest piece of evidence in the paper. As $N$ grows the static-basis method **diverges** (error $2.18$ at $N=8$, $9.64\times10^{-1}$ at $N=10$) while this method decreases monotonically ($4.52\times10^{-5}$ at $N=8$). The paper reads that as demonstrating the necessity of time-dependent basis functions for both the parametric and the spatial variable.

In the Burgers example the basis fields show a clear amplitude hierarchy, from $g_1\sim3.5\times10^{-1}$ down to $g_9\sim3\times10^{-5}$, so the first field carries the core information and the last few carry fine scales. In the Allen-Cahn example the error first decreases with the number of separated terms and then **reaches a floor**, unlike the other examples.

> [!note] Title difference
> The homepage lists this paper as _A dynamical variable-separation method for dynamical systems with random input_, whereas the published version, the preprint and third-party records all give _A Dynamical Variable-Separation Method for Parameter-Dependent Dynamical Systems_. This site records the published version; the preprint is [arXiv:2502.08464](https://arxiv.org/abs/2502.08464).

## How the three relate

| No. | What evolves                                 | How the modes relate                 | Principal risk                          |
| --- | -------------------------------------------- | ------------------------------------ | --------------------------------------- |
| 15  | spatial basis and stochastic coefficients    | evolve together, coupled through $C$ | collapse of the smallest singular value |
| 17  | nothing (a one-off eigensolve)               | not applicable                       | assembling and solving a dense matrix   |
| 83  | a greedily enriched separated representation | decouple into two equations per step | the greedy sequence need not be optimal |

The relation between papers 17 and 15 is worth naming: the error bound in paper 15 is stated against the truncated Karhunen-Loève approximation, and paper 17 solves the problem of **computing that reference object**. Together they give the full cost structure of this route: either pay for the dense eigenproblem that produces the Karhunen-Loève basis, or let the basis evolve and accept the risk carried by the smallest singular value.

Paper 83 supplies a third option: neither compute the optimal basis in advance nor evolve every mode at once, but enrich greedily term by term so each step decouples. Taken together the three map out the choice: **coupled evolution buys quasi-optimality against the best rank-$S$ approximation, decoupled greedy enrichment buys freedom from inverting a covariance, and a one-off eigensolve buys the optimal basis at the price of a dense matrix.**

## Coverage check

| Item                                                              | Paper | Status                                                                |
| ----------------------------------------------------------------- | ----- | --------------------------------------------------------------------- |
| Ansatz and the three gauge conditions                             | 15    | form, conditions, meaning of the gauge                                |
| Three evolution equations and their roles                         | 15    | mean, basis, coefficient equations, projector                         |
| Tangent projection and the origin of $C^{-1}$                     | 15    | tangent space, projection formula, Dirac-Frenkel form                 |
| Theorem 4.1 and the $e^{C/\rho}$ constant                         | 15    | assumptions, conclusion, both readings, curvature                     |
| The eigenvalue-crossing failure mode                              | 15    | construction and what it demonstrates                                 |
| Handling a singular covariance                                    | 15    | the pseudoinverse trap, reformulation, per-step diagonalisation       |
| Structural difference of the Fredholm problem                     | 17    | compactness, spectral accumulation, dense matrix                      |
| Multilevel correction step and integral iterations                | 17    | the two generic steps, dimension argument, quadrature versus solve    |
| What each of the three reductions freezes                         | 83    | polynomial chaos, proper orthogonal decomposition, Kolmogorov barrier |
| Model, affine assumption and ansatz                               | 83    | structural split, affinity, dropping mean field and gauge             |
| Greedy rule and the two decoupled subproblems                     | 83    | selection rule, stopping, origin of both equations                    |
| Initial conditions and the Gram-Schmidt-like correction           | 83    | both formulas and their relation to the gauge                         |
| Residual driving and the carried $\langle(g_k)_t,g_k\rangle$ term | 83    | both equations in the linear case, where the difference lands         |
| Closed-form recursion with a scalar divisor                       | 83    | recursion, $c_{n+1}$, $l_{n+1}$, absence of $C^{-1}$                  |
| Offline-online split and the online cost                          | 83    | stored scalars, scalar online equation, independence of the mesh      |
| A posteriori bound and the logarithmic Lipschitz constant         | 83    | definition, $\delta_k$, possible negativity and its effect            |
| Absence of a convergence theorem, stated limitation               | 83    | the conclusion's wording and the storage dependence                   |
| Both numerical comparisons                                        | 83    | the reduced-basis trade-off, the static-basis divergence              |

## Sources for this page

- E. Musharbash, F. Nobile, and T. Zhou, [_Error analysis of the dynamically orthogonal approximation of time dependent random PDEs_](https://doi.org/10.1137/140967787), SIAM J. Sci. Comput. 37(2) (2015), pp. A776-A810.
- H. Xie and T. Zhou, [_A multilevel finite element method for Fredholm integral eigenvalue problems_](https://doi.org/10.1016/j.jcp.2015.09.043), J. Comput. Phys. 303 (2015), pp. 173-184.
- L. Chen, Y. Chen, Q. Li, and T. Zhou, [_A dynamical variable-separation method for parameter-dependent dynamical systems_](https://doi.org/10.1137/24M168427X), SIAM J. Sci. Comput. 47(3) (2025), pp. A1783-A1808 (preprint [arXiv:2502.08464](https://arxiv.org/abs/2502.08464)).
