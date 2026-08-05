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
> Papers **15** (_SIAM J. Sci. Comput._ 37(2), 2015), **17** (_J. Comput. Phys._ 303, 2015) and **83** (_SIAM J. Sci. Comput._ 47(3), 2025). Papers 15 and 83 are checked equation by equation against publicly available full texts, and this page gives their complete derivations, the full hypotheses of their theorems, and their numerical experiments. The journal text of paper 17 requires a subscription and has no preprint, so that section reports only what the abstract and keywords confirm and keeps the general framework it adapts separate from the paper's own claims. That section is visibly shorter than the other two; this reflects the state of the evidence, not an omission.

![Evolve a random solution directly on a low-rank manifold](assets/diagrams/tao-zhou-papers/en/low-rank-dynamics.svg)

## A low-rank manifold is a reduction, not a truncation

The three papers share a starting point that is easy to misread, so it is worth stating first: **confining the solution to a rank-$S$ representation is not the same thing as cutting a fixed expansion off after $S$ terms.**

Truncation fixes the basis in advance — generalised polynomial chaos polynomials in the random variables, or proper orthogonal decomposition and Karhunen-Loève modes computed from snapshots or from the initial datum — keeps $S$ terms and discards the tail. The unknowns are then $S$ coefficients inside a **fixed linear subspace**, and accuracy is entirely governed by how fast the tail decays in that fixed basis. When the solution's dependence on the random parameters moves with time, the fixed subspace is no longer the right one and the only repair is to add terms: exactly the growth paper 15 opens with, and exactly the Kolmogorov barrier paper 83 opens with.

The low-rank manifold approach changes the state space. Let $\mathcal M_S$ be the set of all rank-$S$ fields. It is not a linear subspace but a smooth manifold, and the unknown is **a point moving on it**, with both factors of the representation unknown. Its dynamics come from a Galerkin (Dirac-Frenkel) condition: require the residual $\partial_t u_S-\mathcal L(u_S)$ to be orthogonal to the tangent space $T_{u_S}\mathcal M_S$ at the current point,

$$
\frac{\partial u_S}{\partial t}=P_{u_S}\bigl(\mathcal L(u_S)\bigr).
$$

Three consequences follow, and they are the whole subject matter of this page.

First, the subspace the solution occupies can rotate and deform **at fixed rank**, because the basis is itself part of the unknown. There is no need to enlarge $S$ merely to follow a moving subspace, and that is how the Kolmogorov barrier is dodged.

Second, the equations can be written entirely in the factors: one partial differential equation per basis function plus $S$ coefficient equations. Cost scales with $S$ and with the cost of one deterministic field solve, not with the dimension of the product space $L^2(D)\otimes L^2(\Omega)$; the statistics are read straight off the factors,

$$
\mathbb E[u(x,t,\cdot)]\approx\bar u_S(x,t),
\qquad
\mathrm{Var}_T[u](t)\approx\sum_{i=1}^{S}\mathbb E\bigl[Y_i^2(t)\bigr],
$$

with no sampling. **That is the substantive sense of "reduction": the ambient dimension never appears in the cost of the evolution.**

Third, the price. $\mathcal M_S$ is curved, its curvature at a point is of order $1/\sigma_S$ — the reciprocal of the smallest retained singular value — and the tangent projection contains the inverse of a covariance matrix explicitly. Everything else on this page grows out of that one inverse: paper 15's $e^{C/\rho}$ constant, its $C^{\dagger}$ workaround, and the entire design motivation of paper 83.

This is also why paper 17 belongs here. The pointwise-in-time best rank-$S$ object is the truncated Karhunen-Loève expansion, defined by the Fredholm integral eigenvalue problem $\mathcal T_u Z_i=\mu_i Z_i$ for the covariance operator. It is the reference object in paper 15's theorem, but computing it at every $t$ presupposes knowing $u(t)$, and $t\mapsto z_S(t)$ need not be continuously differentiable. So there are two routes: make that dense eigenproblem cheaper (paper 17), or replace the instantaneous eigendecomposition by a flow on the manifold (papers 15 and 83). Even the second route needs one eigensolve at $t=0$, because the dynamically orthogonal flow is launched from the Karhunen-Loève expansion of $u_0$.

## 15: an error theory for the dynamically orthogonal approximation

### The idea: why a fixed basis fails, and what dynamical low rank offers instead

Time-dependent random PDEs are usually expanded in a **fixed** basis: generalised polynomial chaos in the random variables, proper orthogonal decomposition modes in space. The difficulty is that the solution's dependence on the random parameters changes substantially as time advances, so a fixed basis needs a steadily growing number of terms to hold accuracy; adaptive and greedy repairs (time-dependent polynomial chaos, generalised proper orthogonal decomposition) only partly help.

The dynamically orthogonal approximation of Sapsis and Lermusiaux evolves the spatial basis and the stochastic coefficients **together**, which sidesteps the growth but had essentially no error theory: existing results covered the closely related multiconfiguration time-dependent Hartree method and the dynamical low-rank matrix setting, not random PDEs. The paper's first goal is therefore to **establish a precise correspondence between the dynamically orthogonal and dynamical low-rank approximations**, so the matrix theory can be imported wholesale. The correspondence reads as a dictionary:

| Dynamical low rank (Koch-Lubich, matrices)              | Dynamically orthogonal (random fields)                                        |
| ------------------------------------------------------- | ----------------------------------------------------------------------------- |
| $Y=USV^{T}\in\mathcal M_r$                              | $u_S=\bar u_S+\mathbf U^{T}\mathbf Y\in\mathcal M_S$                          |
| left factor $U$, orthonormal columns                    | spatial basis $\{U_i(\cdot,t)\}$, $\langle U_i,U_j\rangle=\delta_{ij}$        |
| right factor $V$ and core $S$                           | stochastic coefficients $\{Y_i(t,\cdot)\}$, $\mathbb E[Y_i]=0$                |
| gauge $U^{T}\dot U=0$                                   | gauge $\langle\partial_t U_i,U_j\rangle=0$                                    |
| $S^{-1}$ appears in the $\dot U$ and $\dot V$ equations | $C^{-1}$ appears in the tangent projection, $\sigma(u_S)=\sqrt{\mathrm{eig}(C)}$ |

The last row is the hinge of the whole paper: the $S^{-1}$ that stiffens the matrix flow as $\sigma_r\to0$ is, in the random-field setting, precisely $C^{-1}$.

### Setting: model, ansatz and the three gauge conditions

On an open bounded $D\subset\mathbb R^d$ with $1\le d\le3$ and a complete probability space $(\Omega,\mathcal A,\mathbb P)$,

$$
\frac{\partial u(x,t,\omega)}{\partial t}=\mathcal L\bigl(u(x,t,\omega),\omega\bigr),
\qquad x\in D,\ t\in[0,T],\ \omega\in\Omega,
$$

$$
u(x,0,\omega)=u_0(x,\omega),
\qquad
\mathcal B\bigl(u(\sigma,t,\omega)\bigr)=h(\sigma,t),\ \sigma\in\partial D .
$$

The error analysis specialises to the parabolic case, where the operator is

$$
\mathcal L(u):=\nabla\!\cdot\!(a\nabla u)+f,
\qquad
\mathcal L^{*}(\cdot):=\mathcal L(\cdot)-\mathbb E[\mathcal L(\cdot)] .
$$

The constants $a_{\min}$ and $a_{\max}$ in the theorem below belong to the diffusion coefficient $a$ here. The ansatz is

$$
u_S(x,t,\omega)=\bar u_S(x,t)+\sum_{i=1}^{S}U_i(x,t)\,Y_i(t,\omega),
$$

and uniqueness of the representation is fixed by three conditions,

$$
\mathbb E[Y_i(t,\cdot)]=0,
\qquad
\langle U_i(\cdot,t),U_j(\cdot,t)\rangle=\delta_{ij},
\qquad
\Bigl\langle \frac{\partial U_i(\cdot,t)}{\partial t},\,U_j(\cdot,t)\Bigr\rangle=0,
$$

for $1\le i,j\le S$ and all $t$, with $\langle u,v\rangle=\int_D uv\,\mathrm dx$. The third is the **gauge condition**: the time derivative of the spatial basis must be orthogonal to the span of the basis itself, so the basis is not allowed to rotate within its own span. That is what removes the redundant freedom in the decomposition — one point of $\mathcal M_S$ admits infinitely many factorisations, and without fixing that freedom the factor equations are underdetermined.

The asymmetry of these conditions is worth recording: only the spatial basis is constrained, the stochastic side is not. The paper says so itself and gives the **doubly dynamically orthogonal** form $u_S^{*}=\mathbf U^{T}A\tilde{\mathbf Y}$, which orthogonalises the stochastic side too and is what the analysis actually uses.

The boundary and initial conditions follow from the ansatz as well:

$$
\mathcal B\bigl(\bar u_S(\sigma,t)\bigr)=h(\sigma,t),
\qquad
\sum_{i=1}^{S}C_{ij}(t)\,\mathcal B\bigl(U_i(\sigma,t)\bigr)=0,
$$

$$
\bar u_S(x,0)=\mathbb E[u_0(x,\cdot)],
\qquad
U_i(x,0)=Z_{i0}(x),
\qquad
Y_i(0,\omega)=\bigl\langle u_0(\cdot,\omega)-\bar u_0,\,Z_{i0}\bigr\rangle,
$$

where $\{Z_{i0}\}$ are the Karhunen-Loève spatial modes of the initial datum: **the dynamically orthogonal flow is ignited from the Karhunen-Loève expansion of $u_0$.**

### Derivation: three evolution equations with three distinct roles

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

with covariance matrix $C_{ij}(t)=\mathbb E[Y_iY_j]$, centred operator $\mathcal L^{*}(u,\omega)=\mathcal L(u,\omega)-\mathbb E[\mathcal L(u)]$, and

$$
\Pi^{\perp}_{U}[v]=v-\Pi_U[v]=v-\sum_{i=1}^{S}\langle v,U_i\rangle U_i
$$

the $L^2(D)$ projector onto the complement of $\mathcal U=\mathrm{span}\{U_1,\dots,U_S\}$.

The three equations play different roles. The mean equation is the plain averaged PDE. The basis equation is driven by the **part projected out of the span** of the correlation between the operator and the stochastic coefficients, and modes couple only through $C$. The coefficient equation is an ODE in $\omega$ obtained by testing the centred operator against each basis function. Together they still amount to $S+1$ spatial fields and $S$ scalar equations parametrised by $\omega$ — the concrete form of "cost scales with $S$, not with the ambient dimension".

### Derivation: the tangent projection and where $C^{-1}$ comes from

Writing the centred part as $u_S^{*}=\mathbf U^{T}\mathbf Y$ places it on the rank-$S$ manifold $\mathcal M_S$. The tangent space consists of $\delta u_S^{*}=\mathbf U^{T}\delta\mathbf Y+\delta\mathbf U^{T}\mathbf Y$, with the decomposition unique once the analogue of the gauge is imposed, $\langle\delta\mathbf U,\mathbf U^{T}\rangle=\mathbf 0$ and $\mathbb E[\delta\mathbf Y]=\mathbf 0$. The orthogonal projection onto the tangent space is

$$
P_{u_S^{*}}(v)=P_{u_S^{*}}(v^{*})=\mathbf U^{T}\langle v^{*},\mathbf U^{T}\rangle
+\Bigl(\Pi^{\perp}_{U}\bigl\{\mathbb E[v^{*}\mathbf Y^{T}]\bigr\}\,C^{-1}\Bigr)^{T}\mathbf Y,
$$

where $C=\mathbb E[\mathbf Y\mathbf Y^{T}]$ has full rank by the definition of a rank-$S$ function. **This is where the $\rho^{-1}$ that later appears in the error bound enters.** The complementary projector factorises as $P^{\perp}_{u_S^{*}}v=\Pi^{\perp}_{U}\otimes\Pi^{\perp}_{\tilde{\mathcal Y}}v$, with $\tilde{\mathcal Y}=\mathrm{span}\{\tilde Y_1,\dots,\tilde Y_S\}$ coming from the doubly dynamically orthogonal form — the counterpart of $P^{\perp}(Y)B=P^{\perp}_UBP^{\perp}_V$ in the matrix setting.

$\rho$ is not an abstract quantity. The singular values of a rank-$S$ random field are defined by $\sigma(u_S):=\sigma(A)=\sqrt{\mathrm{eig}(C)}$, with $A$ the square root of $C$, so **$\rho$ is a lower bound on $\sqrt{\lambda_{\min}(C(t))}$** — something visible in the computation.

In Dirac-Frenkel form the method is the Galerkin condition

$$
\mathbb E\Bigl[\Bigl\langle \frac{\partial u_S}{\partial t}-\mathcal L(u_S),\,v\Bigr\rangle\Bigr]=0
\quad \forall v=\bar v+v^{*},\ (\bar v,v^{*})\in H\times T_{u_S^{*}(t)}\mathcal M_S,
$$

equivalently

$$
\frac{\partial u_S}{\partial t}=\mathbb E[\mathcal L(u_S)]
+P_{u_S^{*}(t)}\bigl(\mathcal L^{*}(u_S)\bigr).
$$

That formulation compresses "evolve basis and coefficients together" into a single sentence: **project the right-hand side onto the tangent space of the manifold.**

### Theorems: the curvature lemma

The first building block of the error analysis is a curvature lemma. Its matrix prototype (Koch-Lubich) reads: let $X\in\mathcal M_r$ have smallest nonzero singular value $\sigma_r(X)\ge\rho>0$ and let $Y\in\mathcal M_r$ satisfy $\|Y-X\|\le\tfrac18\rho$; then for every $B$,

$$
\bigl\|\bigl(P(Y)-P(X)\bigr)B\bigr\|\ \le\ 8\rho^{-1}\,\|Y-X\|\cdot\|B\|_2,
\qquad
\bigl\|P^{\perp}(Y)(Y-X)\bigr\|\ \le\ 4\rho^{-1}\,\|Y-X\|^{2},
$$

in the Frobenius norm, with $\|B\|_2$ the spectral norm. This paper's Lemma 3.1 is **the same statement with the same constants $8\rho^{-1}$ and $4\rho^{-1}$**, transported to $\mathcal M_S\subset L^2(D)\otimes L^2(\Omega)$; the paper attributes it to Conte and Lubich with small adjustments of notation. Both inequalities say one thing: how fast the tangent projection varies from point to point, and how far the manifold departs from its tangent plane, are both controlled by $\rho^{-1}$. **That is where "curvature is inversely proportional to the smallest singular value" becomes quantitative.**

### Theorem 4.1: quasi-optimality and its price

The reference object is the best rank-$S$ approximation at each time, the truncated Karhunen-Loève expansion

$$
z_S(x,\omega,t)=\bar u(x,t)+\sum_{i=1}^{S}\sqrt{\mu_i(t)}\,\gamma_i(t,\omega)\,Z_i(x,t),
$$

with $\{\mu_i,Z_i\}$ the eigenpairs of the covariance operator. The theorem has two hypotheses: this best approximation exists and is continuously differentiable in $(H^2(D)\cap H^1_0(D))\otimes L^2(\Omega)$ for $0\le t\le\bar t$, and its smallest singular value is uniformly bounded below,

$$
\sigma\bigl(z_S(t)\bigr)\ \ge\ \rho\ >\ 0,\qquad \forall t\in[0,\bar t].
$$

Then there is $0<\hat t\le\bar t$ such that the dynamically orthogonal solution $u_S=\bar u_S+u_S^{*}$ started at $u_S(0)=z_S(0)$ satisfies, for all $0<t\le\hat t$,

$$
\|u_S(t)-z_S(t)\|_0^2+a_{\min}\int_0^t |u_S(\tau)-z_S(\tau)|_1^2\,\mathrm d\tau
\ \le\ 2\alpha\,e^{2\beta(t)}\int_0^t \|z_S(\tau)-u(\tau)\|_1^2\,\mathrm d\tau,
$$

$$
\beta(t)=4\rho^{-1}\int_0^t\Bigl(4\|\mathcal L^{*}(z_S)\|_0+\|\mathcal L^{*}(u)\|_0
+\|\mathcal L^{*}(u_S)\|_0+\|\dot z_S^{*}\|_0^2\Bigr)\mathrm d\tau,
\qquad
\alpha=\max\Bigl\{\frac{a_{\max}^2}{2a_{\min}},\ 4\rho^{-1}\Bigr\},
$$

where $\|\cdot\|_1$ and $|\cdot|_1$ are the norm and seminorm of $H^1(D)\otimes L^2(\Omega)$, provided every term in the bound is well defined. The proof follows Lubich and coauthors' dynamical low-rank argument for time-dependent data matrices.

The conclusion reads on two levels. **The good news**: the dynamically orthogonal error is bounded by the best rank-$S$ error, so the method is quasi-optimal. **The price**: the constant grows like $e^{C/\rho}$, exponentially in the reciprocal of the smallest singular value. That is the quantitative form of the curvature obstruction — the $4\rho^{-1}$ of the curvature lemma passes untouched into $\beta(t)$ and $\alpha$, so as the $S$-th singular value collapses the tangent projection becomes ill-conditioned and the bound goes vacuous.

It is worth naming where each hypothesis enters. On **regularity**, the analysis needs $u\in L^2\bigl(\mathcal T,(H^2(D)\cap H^1_0(D))\otimes L^2(\Omega)\bigr)$, and the paper shows the corresponding energy bounds also hold for the dynamically orthogonal solution provided $\nabla a\in L^{\infty}(D\times\Omega)$ and $u(0),\dot u(0)\in H^1(D)\otimes L^2(\Omega)$. On **admissibility of the weak form** there is a key technical point (Prop. 3.5): $-\Delta u_S^{*}\in T_{u_S^{*}}\mathcal M_S$, so it is a legitimate test function in the dynamically orthogonal weak form — without it the parabolic energy argument could not be run on the manifold. **$\rho$** enters only through the curvature lemma, and **continuous differentiability** determines the time interval on which the conclusion holds; the next section explains why that hypothesis is not a technicality.

> [!note] Numbering and versions
> The theorem and proposition numbers in this section are those of the MATHICSE technical report, which circulated under the title _On the Dynamically Orthogonal approximation of time dependent random PDEs_. Whether the published SIAM version renumbers them has not been verified here.

### Why the hypothesis is not a technicality: eigenvalue crossing

The paper constructs an eigenvalue-crossing example: a deterministic Laplacian with random initial data whose first two Karhunen-Loève eigenvalues cross at a time $t^{*}$. Both the best rank-1 error and the dynamically orthogonal rank-1 error can then be computed in closed form:

$$
\epsilon_{KL}(t)=\min\bigl\{\mathbb E[\alpha_1^2]e^{-2t},\ \mathbb E[\alpha_2^2]e^{-8t}\bigr\},
\qquad
\epsilon_{KL}(t)=\mathbb E[\alpha_2^2]e^{-8t}\ \ (t>t^{*}),
$$

$$
\epsilon_{DO}(t)=\mathbb E[\alpha_1^2]e^{-2t}=
\begin{cases}
\epsilon_{KL}(t), & t\le t^{*},\\[4pt]
\dfrac{\mathbb E[\alpha_1^2]}{\mathbb E[\alpha_2^2]}\,e^{6t}\,\epsilon_{KL}(t), & t>t^{*}.
\end{cases}
$$

The reading is direct: the dynamically orthogonal error decays at the rate set by the **smallest** Laplacian eigenvalue, while after the crossing the Karhunen-Loève error decays at the rate set by the second, and their ratio grows like $e^{6t}$. **So the dynamically orthogonal error cannot be bounded uniformly by the Karhunen-Loève error.** The paper stresses that this does not contradict Theorem 4.1: at $t^{*}$ the rank-1 truncated Karhunen-Loève expansion fails to be continuously differentiable in time, so the theorem's hypothesis does not hold in the first place. Put differently, the differentiability hypothesis is not a technical artefact — **it is the no-crossing requirement itself.**

> [!note] The obstruction is an obstruction to the bound, not always to the method
> The original matrix literature adds two qualifications worth keeping in mind here. First, if the true $\varepsilon$-pseudorank is $q<r$ but one integrates at rank $r$, then even though the core matrix is ill-conditioned and $S^{-1}$ appears in the factor equations, the approximation is not severely affected. Second, no singularities arise where the singular values of $Y(t)$ coalesce, unlike in the ODEs for a smooth singular value decomposition. The genuine failure mode is a different one: when the rank is too small, a singular value **not** included in the approximation can grow past the ones being tracked without the method noticing — which is exactly the eigenvalue-crossing example this paper builds.

### Numerical experiments: three test problems

The paper reports three test problems, covering three different sources of randomness.

| Test                                     | Equation and domain                            | Randomness and initial datum                 | Purpose                                                                                          |
| ---------------------------------------- | ---------------------------------------------- | -------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| §6.2 linear parabolic equation           | one-dimensional interval                        | random initial condition                      | compare the dynamically orthogonal modes $U_i$ with the Karhunen-Loève modes $Z_i$, track coefficient variances |
| §6.3 linear parabolic equation           | $x\in[0,1]$, homogeneous Dirichlet              | random diffusion coefficient, deterministic initial datum | zero initial covariance, so it stresses the $C^{\dagger}$ treatment                              |
| §6.4 nonlinear reaction-diffusion equation | $D=[0,1]^2$, homogeneous Neumann              | random threshold potential, deterministic step initial datum | travelling wave whose speed is proportional to the excitation rate $\beta$                       |

The first test also includes a run with a regularisation threshold of order $10^{-16}$. The second is specified by

$$
\partial_t u-\mathrm{div}\bigl(a(x,\omega)\partial_x u\bigr)=0,
\qquad
a(x,\omega)=\bar a(x)+\sum_{i=1}^{2}\bigl(\eta_{2i-1}(\omega)\cos(i\pi x)+\eta_{2i}(\omega)\sin(i\pi x)\bigr),
$$

with $\bar a=1.45$, independent zero-mean uniform $\eta_1,\dots,\eta_4$ satisfying $\mathbb E[\eta_i^2]=\tfrac13\cdot10^{-i+1}$, and deterministic initial condition $u_0(x)=10\sin(\pi x)$. The third is

$$
\partial_t u-\mu\Delta u=f(u),
\qquad
f(u)=\beta\,u(u-1)\bigl(\alpha(\omega)-u\bigr),
$$

with a step initial condition ($1$ for $x_1\le0.5$, $0$ otherwise), uniformly distributed threshold potential $\alpha(\omega)$ and constant excitation rate $\beta$.

The qualitative outcome is that the tests confirm the theoretical bound and mark out where the method works and where it does not: as long as the singular values stay separated, the dynamically orthogonal solution tracks the truncated Karhunen-Loève expansion closely, and once the crossing and rank-deficiency regime is entered it degrades exactly as the analysis predicts.

> [!warning] Quantitative results not verified here
> Beyond the closed-form rates of §5, the specific error magnitudes and convergence slopes live in the figures and could not be recovered from the accessible text, so **this page gives no error tables for paper 15**. Concrete numbers should be read from the figures in the original.

### An implementation trap: singular covariance and the pseudoinverse

The second test above targets an easily missed implementation problem. The covariance matrix $C(t^n)$ can be singular or very ill-conditioned — any system with stochastic coefficients and **deterministic** initial data starts with $C\equiv0$. Applying the Moore-Penrose pseudoinverse $C^{\dagger}$ naively zeroes the "non-active" basis functions and thereby **prevents the rank from increasing**. The paper instead reformulates the basis equation as

$$
\frac{\partial\mathbf U}{\partial t}=C^{\dagger}\,\Pi^{\perp}_{U}\,
\mathbb E\bigl[\mathbf Y\,\mathcal L(u_S)\bigr],
$$

which agrees with the original when $C$ has full rank, and in practice diagonalises the covariance at each step to decouple the system, since the flow does not preserve uncorrelatedness of the $Y_i$ even when $C(0)$ is diagonal.

## 17: confining the eigensolve to the coarsest mesh

### The idea: why the eigenproblem for a compact operator is expensive

A Karhunen-Loève expansion needs the eigenpairs of the Fredholm integral operator built from the covariance kernel: given $\mathrm{Cov}(x,y)$, find

$$
\int_D \mathrm{Cov}(x,y)\,u(y)\,\mathrm dy=\lambda\,u(x),
\qquad \|u\|_{L^2(D)}=1 .
$$

Compared with an eigenvalue problem for a differential operator there is a structural difference: the integral operator $\mathcal T$ is **compact and smoothing**, its eigenvalues accumulate at $0$ rather than at $+\infty$ (so the interesting ones are the **largest**), and discretisation produces a **dense** matrix. There is no sparsity to exploit, and both assembly and solution grow sharply with the mesh resolution needed to resolve the kernel.

Earlier work attacked that cost from the discretisation side — wavelet-Galerkin schemes, generalised fast multipole methods, spectral elements, Legendre-Galerkin with tensor structure, multilevel augmentation for compact integral operators — but the eigensolve itself still had to happen on the finest mesh. This paper removes that: the eigenvalue problem becomes a sequence of **integral iterations** plus eigenvalue solves on the **coarsest** mesh only, so any efficient integration scheme can be plugged in, and the abstract states that total work is comparable to one integration step on the finest mesh.

### The framework it adapts

The multilevel correction framework being adapted (Lin and Xie) has a generic step: given a coarsest space $V_H$ and a nested hierarchy $V_H\subset V_{h_2}\subset\cdots\subset V_{h_n}$, first solve a **source problem** on the finer space,

$$
a(\tilde u_{h_{k+1}},v)=\lambda_{h_k}\,b(u_{h_k},v),\qquad \forall v\in V_{h_{k+1}},
$$

then solve the eigenvalue problem in the small augmented space $V_{H,h_{k+1}}=V_H+\mathrm{span}\{\tilde u_{h_{k+1}}\}$,

$$
a(u_{h_{k+1}},v)=\lambda_{h_{k+1}}\,b(u_{h_{k+1}},v),\qquad \forall v\in V_{H,h_{k+1}},
$$

normalised by $b(u_{h_{k+1}},u_{h_{k+1}})=1$. The point of the second step is that the eigenproblem is only ever solved in a space of dimension $\dim V_H+1$, the coarsest space plus one direction, never on the fine mesh. The full scheme solves the eigenproblem once on $V_H$, loops the correction step up the hierarchy, and finishes with one more source solve on $V_{h_n}$ and a Rayleigh quotient.

Two analytical results for that framework are publicly available as well. The first is the Rayleigh-quotient error identity: for the true eigenpair $(\lambda,u)$ and any $0\ne\psi\in V$, setting $\hat\lambda=a(\psi,\psi)/b(\psi,\psi)$,

$$
\hat\lambda-\lambda=\frac{a(u-\psi,u-\psi)}{b(\psi,\psi)}-\lambda\,\frac{b(u-\psi,u-\psi)}{b(\psi,\psi)} .
$$

**This identity is why an $O(\varepsilon)$ eigenfunction error becomes an $O(\varepsilon^2)$ eigenvalue error throughout this family of methods.** The second is the effect of one correction step and the overall outcome:

$$
\varepsilon_{h_{k+1}}(\lambda):=\eta_a(H)\,\varepsilon_{h_k}(\lambda)+\varepsilon_{h_k}^2(\lambda)+\delta_{h_{k+1}}(\lambda),
$$

$$
\|u_{h_n}-u\|_a\lesssim\varepsilon_{h_n}(\lambda),
\qquad
|\lambda_{h_n}-\lambda|\lesssim\varepsilon^2_{h_n}(\lambda),
\qquad
\varepsilon_{h_n}(\lambda)=\sum_{k=1}^{n}\eta_a(H)^{\,n-k}\delta_{h_k}(\lambda),
$$

where $\delta_h(\lambda)$ is the best-approximation error of the eigenfunction in $V_h$ and $\eta_a(H)$ is the coarse-space "smoothing" quantity that must be small enough. The geometric weights $\eta_a(H)^{n-k}$ explain why the coarse-space contributions do not pollute the final accuracy.

For an integral operator the "source problem" of the first step is not a linear solve but **one quadrature**, a single application of $\mathcal T$, which is what the abstract means by integral iterations; its natural form would be $\tilde u_{h_{k+1}}=\Pi_{h_{k+1}}\mathcal T u_{h_k}$, a power-method-like application followed by the same small augmented eigensolve. **That last display, however, is a reconstruction made here from the abstract's wording, not a formula from the paper.**

> [!warning] What could be verified
> The journal text requires a subscription and there is no preprint (Elsevier paywall, no arXiv preprint, no full text on either author's page). Confirmable: the problem setting, the keywords (uncertainty quantification, Karhunen-Loève expansion, Fredholm eigenvalue problem, multigrid finite element) and the contributions as stated in the abstract — error estimates are provided, the computational complexity is analysed, the total work is comparable to a single integration step on the finest mesh, and several numerical experiments validate the efficiency. The generic correction step, the Rayleigh-quotient identity and the error recursion above come from the Lin-Xie framework being adapted and **must not be cited as theorems of this paper**.
>
> Five things remain unreported here: the actual form of the correction step for the integral operator; the norms, constants and hypotheses of this paper's error estimates; the details of the complexity analysis; which covariance kernels and spatial dimensions the experiments used; and what speedups were measured. Given the reference list (Oliveira-Azevedo's spectral element method, Schwab-Todor's generalised fast multipole method, Phoon and coauthors' wavelet-Galerkin scheme), the natural comparison baselines would be direct finite-element or spectral discretisation of the same integral eigenproblem, but that is an inference rather than a verified statement.

## 83: enrich the basis one term at a time so every step decouples

### The idea: what decoupling buys, and what "dynamical" adds to a static separation

The paper states its position sharply. Generalised polynomial chaos writes $u\approx\sum_i\zeta_i(\xi)g_i(x,t)$, **freezing the parametric coefficients in time** while the spatial modes move, and needs substantial prior information from the full model. Proper orthogonal decomposition writes $u\approx\sum_i\zeta_i(t;\xi)g_i(x)$, **freezing the spatial modes in both time and parameter**.

Whether proper orthogonal decomposition succeeds depends on the Kolmogorov widths of the solution manifold $\mathcal M=\{u(\cdot,t;\xi):t\in[0,T],\xi\in\Omega\}$ decaying fast, and that **fails for first-order linear transport and hyperbolic problems** — the Kolmogorov barrier. The dynamically orthogonal method analysed in paper 15 breaks the barrier by letting both factors move, as do its variants (DyBO, dual-DO and DBO, the last equivalent to the first two up to an invertible matrix transformation); but all of them solve a **coupled** system for all $N$ modes at once, and that is precisely where the covariance inverse and its conditioning problems come from.

What "dynamical" adds over a static variable separation can be said concretely. Static variable separation takes $u\approx\sum_i\zeta_i(\xi)g_i(x,t)$: the parametric coefficients $\zeta_i(\xi)$ carry **no time dependence**, so how the modes are mixed at a given parameter value is decided once and holds for all time, and the space-time fields must carry the mode weights by themselves. This paper replaces them by $\zeta_i(t;\xi)$, so the mixture is re-decided at every instant and both factors of every term are time-dependent — the representation becomes a flow on a low-rank manifold rather than a fixed separated expansion. That is not a cosmetic change: the comparison in §5.1.2 below shows the static version **diverging** as $N$ grows while the dynamical one decreases monotonically.

The paper's angle is therefore to keep the doubly time-dependent ansatz but build the modes **one at a time by greedy enrichment**, so each step decouples into two scalar-coefficient subproblems and no covariance is ever inverted.

### Setting: the model and the affine assumption

$$
\frac{\partial u}{\partial t}(x,t;\xi)=F\bigl(u(x,t;\xi);\xi\bigr),
\qquad u(x,0;\xi)=\mu(x;\xi),
\qquad
\mathcal B\bigl(u(x,t;\xi)\bigr)=g(x,t;\xi),
$$

with variational form $\langle\partial_t u(\cdot,t;\xi),v\rangle=\langle F(u(\cdot,t;\xi);\xi),v\rangle$ for all $v\in V$. The structural assumption is $F(u;\xi)=\mathcal C(\xi)+\mathcal A(u;\xi)+\mathcal H(u;\xi)$, with $\mathcal A$ linear and $\mathcal H$ nonlinear, and **affine parameter dependence**:

$$
\mathcal C(\xi)=\sum_{i=1}^{N_C}\kappa^i_C(\xi)\mathcal C^i,
\qquad
\mathcal A(u;\xi)=\sum_{i=1}^{N_A}\kappa^i_A(\xi)\mathcal A^i(u),
\qquad
\mathcal H(u;\xi)=\sum_{i=1}^{N_H}\kappa^i_H(\xi)\mathcal H^i(u),
$$

together with a separable initial condition $\mu(x;\xi)=\sum_{i=1}^{N_{t_0}}p^i(\xi)q^i(x)$. **Affinity is exactly what makes the offline-online split possible**; the paper notes that when it fails, a variable-separation step can first produce an affine approximation with negligible loss of accuracy.

### Setting: an ansatz with no mean field and no orthogonality gauge

$$
u(x,t;\xi)\approx u_N(x,t;\xi):=\sum_{i=1}^{N}\zeta_i(t;\xi)\,g_i(x,t),
$$

with $\{\zeta_i\}$ parameter-dependent, $\{g_i\}$ parameter-independent and **both time-dependent**. Compare the dynamically orthogonal ansatz $u\approx\bar u(x,t)+\sum_i\zeta_i(t;\xi)g_i(x,t)$, which keeps a statistical mean field and requires $\{g_i(\cdot,t)\}$ orthonormal at every $t$. **This paper drops both the mean field and the orthonormality gauge**, and uniqueness comes instead from the sequential greedy construction. That also explains why paper 15's theorem cannot be carried over: that theory is built on exactly those two features.

### Derivation: the greedy rule and the two decoupled subproblems

The first step picks $\xi_1$ arbitrarily and takes $g_1(x,t)$ to be the **full solution** at $\xi=\xi_1$; testing with $g_1$ then gives a parameter-dependent ordinary differential equation for $\zeta_1(t;\xi)$,

$$
\Bigl\langle \frac{\partial\bigl(g_1(x,t)\zeta_1(t;\xi)\bigr)}{\partial t},\,g_1(x,t)\Bigr\rangle
=\bigl\langle F\bigl(g_1(x,t)\zeta_1(t;\xi);\xi\bigr),\,g_1(x,t)\bigr\rangle .
$$

At step $k\ge2$, with the error $e(x,t;\xi):=u-u_{k-1}$, the selection rule is

$$
\xi_k\in\arg\max_{\xi\in\Xi}\triangle_k(\xi),
$$

where $\triangle_k$ is either $\|e\|_{L^2([0,T];V)}$ itself, if affordable, or the a posteriori bound below; the loop stops when $\triangle_k(\xi_k)<\varepsilon$ and otherwise removes $\xi_k$ from the candidate set. Rewriting the equation in terms of the error,

$$
\Bigl\langle \frac{\partial(e+u_{k-1})}{\partial t},v\Bigr\rangle
=\bigl\langle F\bigl((e+u_{k-1});\xi\bigr),v\bigr\rangle,
\qquad \forall v\in V,
$$

the two subproblems are: $g_k(x,t)$ solves this at the **single** parameter $\xi=\xi_k$, a parameter-independent partial differential equation; and setting $\tilde e=g_k\zeta_k$ with $v=g_k$ gives a scalar ordinary differential equation for $\zeta_k(t;\xi)$,

$$
\Bigl\langle \frac{\partial\bigl(g_k\zeta_k+u_{k-1}\bigr)}{\partial t},\,g_k\Bigr\rangle
=\bigl\langle F\bigl(g_k\zeta_k+u_{k-1};\xi\bigr),\,g_k\bigr\rangle .
$$

**There is no coupling between the two: the spatial side solves at one parameter value, the parametric side solves one scalar equation.** That is all "each step decouples" means.

### Derivation: initial conditions recover the gauge at $t=0$

The paper calls the initial conditions "one of the most essential ingredients", and constructs them by $L^2$ matching against the current mode. Writing $g_{k,0}(x)=g_k(x,0)$, we have $g_1(x,0)=\mu(x;\xi_1)$, and the matching condition $\langle u_1(\cdot,0;\xi),g_{1,0}\rangle=\langle\mu(\cdot;\xi),g_{1,0}\rangle$ gives

$$
\zeta_{1,0}(\xi)=\sum_{i=1}^{N_{t_0}}\frac{\langle q^i,g_{1,0}\rangle}{\langle g_{1,0},g_{1,0}\rangle}p^i(\xi);
$$

for $k\ge2$ the initial error $e_0(x;\xi)=\mu(x;\xi)-\sum_{j=1}^{k-1}g_{j,0}(x)\zeta_{j,0}(\xi)$ gives

$$
g_{k,0}(x)=e_0(x;\xi_k)=\mu(x;\xi_k)-\sum_{j=1}^{k-1}g_{j,0}(x)\,\zeta_{j,0}(\xi_k),
$$

$$
\zeta_{k,0}(\xi)=\sum_{i=1}^{N_{t_0}}\frac{\langle q^i,g_{k,0}\rangle}{\langle g_{k,0},g_{k,0}\rangle}p^i(\xi)
-\sum_{j=1}^{k-1}\frac{\langle g_{j,0},g_{k,0}\rangle}{\langle g_{k,0},g_{k,0}\rangle}\zeta_{j,0}(\xi).
$$

**The second sum is a Gram-Schmidt-like correction.** It is where this method recovers, implicitly and only at $t=0$, part of what the dynamically orthogonal gauge condition supplies explicitly. One edge case (the paper's Remark 3.1): if $\mu\equiv0$ then $g_{k,0}=0$ and $\zeta_{k,0}=0$ at every step.

### Derivation: in the linear case the new mode is driven by the previous residual

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

With $u_0\equiv0$ and $r_1\equiv\mathcal C(\xi)$ it reduces to the first step's equation. Note the term $\langle(g_k)_t,g_k\rangle$: in the dynamically orthogonal method it is **forced to vanish** by the gauge condition, whereas here it is simply carried along. That reduces the difference between the two routes to one concrete term.

### Derivation: the discrete recursion divides by a scalar, not a matrix

Split $[0,T]$ into $N_t$ steps of size $\tau=T/N_t$, write $g_{k,n}=g_k(t_n)$ and $\zeta_{k,n}(\xi)=\zeta_k(t_n;\xi)$, use first-order differences and backward Euler, and substitute the affine expansions together with the separated form of $u_{k-1}$. The result is a closed-form recursion

$$
\zeta_{k,n+1}(\xi)=\frac{c_{n+1}\,\zeta_{k,n}(\xi)+s_{n+1}(\xi)}{l_{n+1}(\xi)},
\qquad n=0,\dots,N_t-1,
$$

$$
c_{n+1}=\frac{\langle g_{k,n+1},g_{k,n+1}\rangle}{\tau},
\qquad
l_{n+1}(\xi)=2c_{n+1}-\frac{\langle g_{k,n},g_{k,n+1}\rangle}{\tau}
-\sum_{i=1}^{N_A}\kappa^i_A(\xi)\bigl\langle\mathcal A^i(g_{k,n+1}),g_{k,n+1}\bigr\rangle,
$$

$$
\begin{aligned}
s_{n+1}(\xi)=&\sum_{i=1}^{N_A}\sum_{j=1}^{k-1}\kappa^i_A(\xi)\,\zeta_{j,n+1}(\xi)\bigl\langle\mathcal A^i(g_{j,n+1}),g_{k,n+1}\bigr\rangle
+\sum_{i=1}^{N_C}\kappa^i_C(\xi)\bigl\langle \mathcal C^i,g_{k,n+1}\bigr\rangle\\
&-\sum_{j=1}^{k-1}\frac{\zeta_{j,n+1}(\xi)-\zeta_{j,n}(\xi)}{\tau}\bigl\langle g_{j,n+1},g_{k,n+1}\bigr\rangle
-\sum_{j=1}^{k-1}\zeta_{j,n+1}(\xi)\Bigl\langle \frac{g_{j,n+1}-g_{j,n}}{\tau},\,g_{k,n+1}\Bigr\rangle .
\end{aligned}
$$

Every part of $s_{n+1}$ is assembled from the affine terms and the $k-1$ modes already computed. **This is the paper's answer to the covariance inversion of the dynamically orthogonal method: $l_{n+1}(\xi)$ is a scalar rather than a matrix, so neither $C^{-1}$ nor $C^{\dagger}$ appears anywhere.**

The offline stage does not store $\zeta_{k,n}(\xi)$ at all; it stores only the parameter-independent scalar inner products

$$
\langle g_{k,n+1},g_{k,n+1}\rangle,\quad
\langle g_{k,n},g_{k,n+1}\rangle,\quad
\langle\mathcal A^i(g_{k,n+1}),g_{k,n+1}\rangle,\quad
\langle\mathcal A^i(g_{j,n+1}),g_{k,n+1}\rangle,
$$

$$
\langle g_{j,n+1},g_{k,n+1}\rangle,\quad
\langle g_{j,n},g_{k,n+1}\rangle,\quad
\langle\mathcal C^i,g_{k,n+1}\rangle .
$$

Online, a new parameter $\bar\xi$ needs only one cheap scalar ordinary differential equation,

$$
\alpha(t)\,(\zeta_k)_t(t;\bar\xi)+\beta(t;\bar\xi)\,\zeta_k(t;\bar\xi)=\gamma(t;\bar\xi),
$$

$$
\alpha(t)=\langle g_k(t),g_k(t)\rangle,
\qquad
\beta(t;\bar\xi)=\langle (g_k)_t(t),g_k(t)\rangle-\sum_{i=1}^{N_A}\kappa^i_A(\bar\xi)\langle\mathcal A^i(g_k(t)),g_k(t)\rangle,
$$

with $\gamma$ assembled from the stored scalars. Both $\beta$ and $\gamma$ are affine in the parameter, so **the online cost is independent of the spatial discretisation** of the original problem.

### The nonlinear case: a semi-implicit treatment of Burgers

Take $\mathcal A(u;\xi)=\kappa(\xi)\partial_x^2u$ and $\mathcal H(u;\xi)=-u\partial_xu$. The discretised residual equation contains $\langle e_{n+1}\partial_xe_{n+1},v\rangle$, which the paper approximates **semi-implicitly** by $e_n\partial_xe_{n+1}$ to control cost, and states that the same procedure carries over to Allen-Cahn.

### Theorems: the a posteriori bound, and the absence of a convergence theorem

The bound rests on the **local logarithmic Lipschitz constant** of $F$ at $u$,

$$
L_V[F](u):=\sup_{v\in V,\ v\ne u}\frac{\langle v-u,\;F(v;\xi)-F(u;\xi)\rangle}{\|v-u\|_V^2}.
$$

Testing the error equation $\partial_te=F(u;\xi)-F(u_{k-1};\xi)+r_k$, with $r_k=F(u_{k-1};\xi)-\partial_tu_{k-1}$, against $e$ and applying a comparison lemma gives Proposition A.1: $\|e(x,t;\xi)\|_V\le\delta_k(t;\xi)$ with

$$
\delta_k(t;\xi)=\int_0^{t}\alpha(s;\xi)e^{\int_s^{t}\beta(\tau;\xi)\mathrm d\tau}\mathrm ds
+e^{\int_0^{t}\beta(\tau;\xi)\mathrm d\tau}\|e(x,0;\xi)\|_V,
$$

where $\alpha(t;\xi)=\|r_k(x,t;\xi)\|_V$ and $\beta(t;\xi)=L_V[F](u_{k-1}(x,t;\xi))$, and the selection criterion takes $\triangle_k(\xi)=\int_0^T\delta_k^2\,\mathrm dt$. **For dissipative $F$ the logarithmic Lipschitz constant can be negative**, in which case the exponential factors damp rather than amplify — which is what makes the estimator usable over long horizons.

> [!warning] No convergence theorem
> The conclusion explicitly lists "a rigorous convergence analysis of the proposed method under reasonable assumptions" as future work. The rigorous content is the a posteriori bound above together with the exact algebraic recursion; the claim of reduced complexity and improved efficiency over existing low-rank separation techniques is the abstract's, supported numerically rather than by theorem. The paper also states a limitation of its own: both storage and online computation still depend on the size of the space-time discretisation, and reducing that is left open.

### Numerical experiments: four examples

The accuracy metric is the average relative error $\epsilon$ over $M$ test samples, with reference solutions from finite elements in space plus backward Euler in time (FEM-BE below).

**Example 1: reaction-diffusion with a parameter-dependent boundary condition** ($D=[0,1]$)

$$
\partial_t u+\xi_1u=2\xi_2\,\partial_x^2u+\xi_3,
\qquad
u(x,0;\xi)=u|_{\partial D}=2(x+1)\xi_4,
$$

with $T=1$, $\xi\in[1,3]^4$, mesh $h=0.02$, $\tau=10^{-3}$, training set $|\Xi|=11$ and $M=10^3$ test samples. Two greedy strategies are compared, selecting $\xi_k$ by the true error $e$ and by the estimator $\delta_k$. The competitor is the time-dependent reduced-basis model order reduction of Billaud-Friess and Nouy (MTD below).

| $N$ | MTD $\epsilon$        | This method $\epsilon$ |
| --- | --------------------- | ---------------------- |
| 2   | $2.48\times10^{-4}$   | $3.43\times10^{-4}$    |
| 4   | $9.19\times10^{-6}$   | $1.46\times10^{-4}$    |
| 7   | $2.80\times10^{-5}$   | $4.66\times10^{-5}$    |

| $N$ | Method      | Offline (s) | Online (s) | Total (s)  | Online per sample (s) |
| --- | ----------- | ----------- | ---------- | ---------- | --------------------- |
| 2   | MTD         | $8.43$      | $19.93$    | $28.36$    | $1.99\times10^{-2}$   |
| 2   | This method | $8.01$      | $1.98$     | $9.99$     | $1.98\times10^{-3}$   |
| 4   | MTD         | $18.78$     | $22.48$    | $41.26$    | $2.25\times10^{-2}$   |
| 4   | This method | $16.51$     | $3.96$     | $20.47$    | $3.96\times10^{-3}$   |
| 7   | MTD         | $39.92$     | $118.52$   | $158.44$   | $1.19\times10^{-1}$   |
| 7   | This method | $29.86$     | $7.20$     | $37.06$    | $7.20\times10^{-3}$   |
| —   | FEM-BE      | —           | —          | $3752.14$  | $3.75$                |

These numbers give a frank trade-off. **MTD's error falls faster**: at $N=4$ it is $15.9$ times more accurate than this method. The paper's own explanation is that MTD updates all parametric coefficients at every iteration while this method never modifies the time-parameter basis functions computed at earlier steps. But MTD bottoms out at $N=4$ ($9.19\times10^{-6}$) and is already back up to $2.80\times10^{-5}$ at $N=7$, where its advantage is down to a factor of $1.7$; the paper reports that beyond $N=7$ its error rises quickly, attributed to ill-conditioning of its linear systems for the parametric coefficients. This method's error decreases monotonically across the three values of $N$, and its online time is more than an order of magnitude smaller ($7.20\times10^{-3}$ s per sample at $N=7$ against $1.19\times10^{-1}$ s). The paper's summary: the competitor gives higher accuracy under certain conditions, and this method does better on efficiency and error stability.

**Example 2: two-dimensional heat equation with a parameter-dependent source** ($D=[0,\pi]^2$, $\partial_tu=\kappa(\xi)\Delta u+f$, $\kappa(\xi)=\xi_1$, $T=1$, $u(x,0;\xi)=\sin x_1\sin x_2+1$, $u|_{\partial D}=1$), compared against the **static-basis** variable-separation method $u\approx\sum_i\zeta_i(\xi)g_i(x,t)$. The table is at $t=1$ with $M=10^3$; the time columns are the timings the paper reports.

| $N$ | This method $\epsilon$ | Time (s) | Static basis $\epsilon$ | Time (s) |
| --- | ---------------------- | -------- | ----------------------- | -------- |
| 2   | $9.81\times10^{-4}$    | $0.71$   | $1.24\times10^{-1}$     | $0.47$   |
| 4   | $3.66\times10^{-4}$    | $2.37$   | $4.91\times10^{-2}$     | $1.52$   |
| 6   | $1.77\times10^{-4}$    | $4.98$   | $6.37\times10^{-2}$     | $3.12$   |
| 8   | $4.52\times10^{-5}$    | $8.52$   | $2.18$                  | $5.28$   |
| 10  | $4.27\times10^{-5}$    | $13.01$  | $9.64\times10^{-1}$     | $8.00$   |
| —   | FEM-BE                 | $7.06\times10^{2}$ | —             | —        |

**This is the cleanest piece of evidence in the paper.** The static-basis method reaches its best at $N=4$ ($4.91\times10^{-2}$) and then turns around; by $N=8$ its average relative error of $2.18$ exceeds $1$, meaning the approximation carries no information at all. This method decreases monotonically to $4.27\times10^{-5}$. The paper reads that as demonstrating the necessity of time-dependent basis functions for both the parametric and the spatial variable — which is exactly the content of the word "dynamical" discussed above.

**Example 3: Burgers equation with a parameter-dependent initial condition** ($x\in[0,1]$, $T=2$)

$$
\partial_tu+u\,\partial_xu=\frac{\xi_1}{50}\,\partial_x^2u,
\qquad
u(x,0;\xi)=x(1-x)^2\xi_2,
\qquad
u(0,t)=u(1,t)=0,
$$

with $\xi\in[1,3]^2$, $h=0.01$, $\tau=10^{-4}$, $|\Xi|=12$ and $M=10^3$.

| $N$ | $\epsilon$ at $t=1$   | Online (s)            | $\epsilon$ at $t=2$   | Online (s)            |
| --- | --------------------- | --------------------- | --------------------- | --------------------- |
| 2   | $1.17\times10^{-2}$   | $1.63\times10^{-2}$   | $0.87\times10^{-2}$   | $2.92\times10^{-2}$   |
| 4   | $1.09\times10^{-3}$   | $5.48\times10^{-2}$   | $2.62\times10^{-3}$   | $9.87\times10^{-2}$   |
| 6   | $2.17\times10^{-4}$   | $1.23\times10^{-1}$   | $4.10\times10^{-4}$   | $2.17\times10^{-1}$   |
| 8   | $5.68\times10^{-5}$   | $2.18\times10^{-1}$   | $1.07\times10^{-4}$   | $3.82\times10^{-1}$   |
| 10  | $2.76\times10^{-5}$   | $3.44\times10^{-1}$   | $7.12\times10^{-5}$   | $5.92\times10^{-1}$   |
| —   | FEM-BE                | $10.72$               | FEM-BE                | $21.47$               |

The error decreases monotonically in $N$ at both times. Only at $N=2$ is the $t=2$ error smaller than the $t=1$ error; from $N=4$ on, the $t=2$ error is uniformly larger, so the accumulation over the longer horizon only becomes visible once enough terms are used. The basis fields show a clear amplitude hierarchy:

| Field     | $g_1$                | $g_2$              | $g_3$              | $g_5$              | $g_7$                | $g_9$              |
| --------- | -------------------- | ------------------ | ------------------ | ------------------ | -------------------- | ------------------ |
| Amplitude | $3.5\times10^{-1}$   | $2\times10^{-2}$   | $6\times10^{-3}$   | $1\times10^{-3}$   | $1.2\times10^{-4}$   | $3\times10^{-5}$   |

so the first field carries the core information of the solution and the last few carry the fine scales.

**Example 4: Allen-Cahn equation** ($D=[0,1]^2$, $T=1$)

$$
\partial_tu=\xi^2\Delta u-f(u),
\qquad
f=F',\quad F(u)=\tfrac14(u^2-1)^2,
$$

with $u(x,0;\xi)=\sqrt5\,(x_1^2-x_1)(x_2^2-x_2)$, homogeneous Dirichlet data, $\xi\in[0.1,0.2]$, $h_{x_1}=h_{x_2}=0.05$, $\tau=10^{-4}$, $|\Xi|=8$ and $M=10^3$. The outcome differs from the other three examples: the average relative error first decreases with the number of separated terms and then **stabilises**, that is, an accuracy floor appears. The magnitude of that floor could not be verified here, so no table is given.

> [!note] Title difference
> The homepage lists this paper as _A dynamical variable-separation method for dynamical systems with random input_, whereas the published version, the preprint and third-party records all give _A Dynamical Variable-Separation Method for Parameter-Dependent Dynamical Systems_. This site records the published version; the preprint is [arXiv:2502.08464](https://arxiv.org/abs/2502.08464). Equation and algorithm numbers in this section are those of the arXiv version, and whether the SIAM version renumbers them has not been verified here.

## How the three relate

| No. | What evolves                                 | How the modes relate                 | Principal risk                          |
| --- | -------------------------------------------- | ------------------------------------ | --------------------------------------- |
| 15  | spatial basis and stochastic coefficients    | evolve together, coupled through $C$ | collapse of the smallest singular value |
| 17  | nothing (a one-off eigensolve)               | not applicable                       | assembling and solving a dense matrix   |
| 83  | a greedily enriched separated representation | decouple into two equations per step | the greedy sequence need not be optimal |

The relation between papers 17 and 15 is worth naming: the error bound in paper 15 is stated against the truncated Karhunen-Loève approximation, and paper 17 solves the problem of **computing that reference object** — the singular values $\sqrt{\mu_i}$ that govern the constant in Theorem 4.1 are precisely the quantities paper 17 computes. Together they give the full cost structure of this route: either pay for the dense eigenproblem that produces the Karhunen-Loève basis, or let the basis evolve and accept the risk carried by the smallest singular value.

Between papers 83 and 15 there is a deliberate trade. Paper 15 shows the dynamically orthogonal method is quasi-optimal against the truncated Karhunen-Loève expansion but with a constant behaving like $e^{C/\rho}$, caused by the $C^{-1}=(\mathbb E[\mathbf Y\mathbf Y^T])^{-1}$ in the tangent projection and visible in paper 15's own $C^{\dagger}$ workaround. Paper 83 never forms $C$ at all: modes are built one at a time and the coefficient update divides by the scalar $l_{n+1}(\xi)$. The price is stated frankly by paper 83's own comparison — earlier modes are never revisited, so the error falls more slowly per mode than for methods that re-solve for all coefficients; it buys stability and cost rather than rate. Paper 83 also drops both features paper 15's theory is built around, the mean field $\bar u_S$ and the gauge $\langle\partial_tU_i,U_j\rangle=0$, which is why paper 15's theorem does not transfer and why paper 83 has no convergence theorem.

There is a secondary link between papers 83 and 17: both replace a large coupled problem with a sequence of one-direction-at-a-time enrichments of a small space (greedy enrichment in time here, mesh-hierarchy correction there), and both are motivated by uncertainty quantification. Taken together the three map out the choice: **coupled evolution buys quasi-optimality against the best rank-$S$ approximation, decoupled greedy enrichment buys freedom from inverting a covariance, and a one-off eigensolve buys the optimal basis at the price of a dense matrix.**

## Coverage check

| Item                                                              | Paper   | Status                                                                |
| ----------------------------------------------------------------- | ------- | --------------------------------------------------------------------- |
| Why a low-rank manifold is a reduction, not a truncation          | general | manifold state space, tangent condition, cost and curvature           |
| Dictionary between dynamical low rank and dynamical orthogonality | 15      | five rows, including $S^{-1}\leftrightarrow C^{-1}$                   |
| Model, parabolic specialisation and ansatz                        | 15      | operator form, three gauge conditions, boundary and initial data      |
| Three evolution equations and their roles                         | 15      | mean, basis, coefficient equations, projector                         |
| Tangent projection and the origin of $C^{-1}$                     | 15      | tangent space, projection formula, $\sigma=\sqrt{\mathrm{eig}(C)}$, Dirac-Frenkel form |
| Curvature lemma                                                   | 15      | both inequalities, constants $8\rho^{-1}$ and $4\rho^{-1}$            |
| Theorem 4.1 and the $e^{C/\rho}$ constant                         | 15      | assumptions, conclusion, both readings, where each hypothesis enters  |
| Closed-form eigenvalue-crossing failure mode                      | 15      | both closed-form rates, the $e^{6t}$ ratio, relation to the hypothesis |
| Three test problems                                               | 15      | setups complete; error values unverified and flagged as such          |
| Handling a singular covariance                                    | 15      | the pseudoinverse trap, reformulation, per-step diagonalisation       |
| Structural difference of the Fredholm problem                     | 17      | compactness, spectral accumulation, dense matrix                      |
| Multilevel correction step and integral iterations                | 17      | the two generic steps, dimension argument, quadrature versus solve    |
| Rayleigh-quotient identity and error recursion of the framework   | 17      | explicitly marked as framework results, not theorems of this paper    |
| The five unreported items                                         | 17      | correction step, constants, complexity details, kernels and dimensions, speedups |
| What each of the three reductions freezes                         | 83      | polynomial chaos, proper orthogonal decomposition, Kolmogorov barrier |
| What "dynamical" adds over static separation                      | 83      | time dependence of the parametric coefficients and its numerical consequence |
| Model, affine assumption and ansatz                               | 83      | structural split, affinity, dropping mean field and gauge             |
| Greedy rule and the two decoupled subproblems                     | 83      | selection rule, stopping, both equations in full                      |
| Initial conditions and the Gram-Schmidt-like correction           | 83      | all three formulas and their relation to the gauge                    |
| Residual driving and the carried $\langle(g_k)_t,g_k\rangle$ term | 83      | both equations in the linear case, where the difference lands         |
| Closed-form recursion with a scalar divisor                       | 83      | recursion, $c_{n+1}$, $l_{n+1}$, full $s_{n+1}$                       |
| Offline-online split and the online cost                          | 83      | the seven stored scalar families, scalar online equation, mesh independence |
| A posteriori bound and the logarithmic Lipschitz constant         | 83      | definition, $\delta_k$, possible negativity and its effect            |
| Absence of a convergence theorem, stated limitation               | 83      | the conclusion's wording and the storage dependence                   |
| Four numerical examples                                           | 83      | three result tables, amplitude hierarchy, the Allen-Cahn floor        |

## Sources for this page

- E. Musharbash, F. Nobile, and T. Zhou, [_Error analysis of the dynamically orthogonal approximation of time dependent random PDEs_](https://doi.org/10.1137/140967787), SIAM J. Sci. Comput. 37(2) (2015), pp. A776-A810.
- H. Xie and T. Zhou, [_A multilevel finite element method for Fredholm integral eigenvalue problems_](https://doi.org/10.1016/j.jcp.2015.09.043), J. Comput. Phys. 303 (2015), pp. 173-184.
- L. Chen, Y. Chen, Q. Li, and T. Zhou, [_A dynamical variable-separation method for parameter-dependent dynamical systems_](https://doi.org/10.1137/24M168427X), SIAM J. Sci. Comput. 47(3) (2025), pp. A1783-A1808 (preprint [arXiv:2502.08464](https://arxiv.org/abs/2502.08464)).
- The curvature lemma and the dynamical low-rank background come from O. Koch and C. Lubich, _Dynamical low-rank approximation_, SIAM J. Matrix Anal. Appl. 29(2) (2007), pp. 434-454.
- The generic multilevel correction framework cited in the section on paper 17 comes from Q. Lin and H. Xie, [_A multi-level correction scheme for eigenvalue problems_](https://doi.org/10.1090/s0025-5718-2014-02825-1), Math. Comp. 84(291) (2014), pp. 71-88 (preprint [arXiv:1107.0223](https://arxiv.org/abs/1107.0223)).
