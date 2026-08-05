---
title: Stochastic Galerkin and Collocation
description: Papers 1, 2, 3, 4, 5, 7 and 38 - algebra of the coupled system and regularity in the random variable
lang: en
translation: computational-mathematics/paper-notes/stochastic-approximation/stochastic-galerkin-and-collocation
tags:
  - paper-notes
  - uncertainty-quantification
  - stochastic-galerkin
---

> [!note] Coverage of this page
> Papers **1** (_J. Comput. Phys._ 229, 2010), **2** (_Commun. Comput. Phys._ 8, 2010), **3** (_J. Comput. Appl. Math._ 236, 2011), **4** (_J. Comput. Math._ 30, 2012), **5** (_J. Sci. Comput._ 51, 2012), **7** (_Adv. Appl. Math. Mech._ 6, 2014) and **38** (_J. Sci. Comput._ 80, 2019).

These seven belong to the **intrusive** route: expand the solution in generalised polynomial chaos and either apply Galerkin projection or solve decoupled problems at parameter points. Two technical questions recur: the algebraic properties of the coupled system, and the regularity of the solution **in the random variable**.

## 1: when is the coefficient matrix of the coupled system diagonally dominant

With a random field depending affinely on the parameters, $\kappa(x,y)=\kappa_0(x)+\sum_{i=1}^{N}\kappa_i(x)y_i$ with $\kappa_0>0$ and uniform ellipticity $\kappa\ge\kappa_{\min}>0$, the generalised polynomial chaos expansion $u\approx\sum_{m=0}^{M}v_m(x,t)\Phi_m(y)$ with $M=\binom{N+P}{N}$ and $\mathbb E[\Phi_m\Phi_n]=\delta_{mn}$ plus Galerkin projection gives $\partial_t v=\nabla\cdot(A\nabla_x v)+f$ with

$$
a_{jk}=\sum_{i=0}^{N}\kappa_i(x)\,e_{ijk},
\qquad
e_{ijk}=\int y_i\,\Phi_j(y)\Phi_k(y)\rho(y)\,\mathrm dy,
$$

and $A=A^{T}$ by construction.

Whether $A(x)$ is strictly diagonally dominant decides whether diagonal/off-diagonal splitting solvers — mixed explicit-implicit time stepping, Jacobi iteration, decoupled preconditioned conjugate gradients — are stable and converge at a mesh-independent rate. Xiu and Shen had proved $A$ positive definite pointwise with at most $2N+1$ nonzero entries per row (since $y_i\Phi_j$ couples only to nearest neighbours in each direction), and had proved a strengthened dominance $a_{jj}\ge\kappa_{\min}+\sum_{k\ne j}|a_{jk}|$ for **symmetric** Beta densities and for Gamma densities, leaving the general asymmetric Beta case as an open question.

### Reducing it to an inequality on three-term recurrence coefficients

The technical route is to reduce the matrix question to an inequality on Jacobi recurrence coefficients. Writing the normalised recurrence in direction $i$ as $xP_j=a^i_jP_{j+1}+b^i_jP_j+c^i_jP_{j-1}$,

$$
a_{jj}=\kappa_0(x)+\sum_{i=1}^{N}\kappa_i(x)\,b^i_j,
\qquad
\sum_{k\neq j}|a_{jk}|=\sum_{i=1}^{N}|\kappa_i(x)|\,|a^i_j+c^i_j| ,
$$

and since all $y_i$ share a distribution the superscript $i$ drops. **That step is the technical core of the paper**: it trades a question about an $M\times M$ matrix for an elementary inequality on two or three recurrence coefficients.

### Theorem 1

If the $y_i$ have an identical Beta density on $(-1,1)$, $\rho(y_i)=(1-y_i)^{\alpha}(1+y_i)^{\beta}$ with $\alpha,\beta>-1$ satisfying

$$
|\alpha|\ \ge\ \tfrac12,\qquad |\beta|\ \ge\ \tfrac12,
$$

or an identical Gamma density on $(0,+\infty)$, $\rho(y_i)=y_i^{\alpha}e^{-y_i}$ with $\alpha>-1$, then the corresponding $A(x)$ is strictly diagonally dominant for all $x\in\Omega$:

$$
a_{jj}\ \ge\ \kappa_{\min}+\sum_{k\ne j}|a_{jk}|,
\qquad 1\le j\le M .
$$

The open question is thus answered affirmatively but **conditionally**: asymmetric Beta parameters are allowed at the price of the explicit threshold $|\alpha|,|\beta|\ge1/2$. The paper does not claim dominance for all $\alpha\ne\beta$, and whether the threshold can be removed is not settled there.

One convention deserves note: with $\rho(y_i)=(1-y_i)^\alpha(1+y_i)^\beta$ and $\alpha,\beta>-1$, the condition $|\alpha|\ge1/2$ admits both $\alpha\ge1/2$ and $-1<\alpha\le-1/2$; **the uniform (Legendre) case $\alpha=\beta=0$ is not covered by Theorem 1** but is covered by the earlier symmetric result.

The paper reports no numerical experiments: its six pages are the problem restatement, the recurrence-coefficient algebra and the proof.

## 2 and 4: a random wave speed makes the solution barely regular in the random variable

### 2: convergence for collocation derived from hypotheses on the data

For the scalar hyperbolic equation $\partial_t u=c(y)\partial_x u$ on $x\in(-1,1)$, a random wave speed $c(y)$ that can **change sign** makes the inflow boundary switch ends of the domain: data are prescribed at $x=-1$ when $c(y)<0$ and at $x=1$ when $c(y)>0$. So even with smooth data the solution has low regularity in the random variable.

Gottlieb and Xiu's convergence proof in the Galerkin case assumed fast asymptotic decay of the expansion coefficients — an assumption about the **answer**, not about the **data**. This paper analyses collocation, where boundary and initial conditions are imposed pointwise in $y$ and therefore trivially, and derives convergence rates from hypotheses on the data.

The technical device is a set of weighted data hypotheses dividing the density by the wave speed:

$$
\int_\Gamma\!\!\int_D \rho(y)\bigl(\partial_y u_0(x;y)\bigr)^2\mathrm dx\,\mathrm dy<\infty,
$$

$$
\int_0^T\!\!\int_{\Gamma^+}\frac{\rho(y)}{c(y)}\bigl(\partial_y u_R(t;y)\bigr)^2\mathrm dy\,\mathrm dt<\infty,
\qquad
\int_0^T\!\!\int_{\Gamma^-}\frac{\rho(y)}{|c(y)|}\bigl(\partial_y u_L(t;y)\bigr)^2\mathrm dy\,\mathrm dt<\infty,
$$

where $\Gamma^{\pm}$ are the subsets on which $c$ is positive and negative. **The $1/|c|$ weight is exactly what couples the wave speed to the admissible boundary data**: parameter regions where the speed is near zero demand smoother boundary data there. Under $|c'(y)|\le C$ plus these finiteness conditions the solution has $H^1$ regularity in $y$ (energy estimate plus Gronwall), and second-derivative hypotheses give $H^2$.

### 4: why exponential convergence is observed

Papers 2 and 5 give finite-order regularity results ($H^k$, $BV$) with correspondingly **algebraic** rates, while numerically exponential convergence is observed and unexplained. Paper 4 closes that gap: it establishes **analytic** regularity in the random variable — an analytic extension to a complex neighbourhood of the parameter interval — and converts it into spectral (exponential) convergence. The paper states explicitly that it deliberately ignores the deterministic-solver error and studies only the random-space discretisation error.

The key hypothesis imposes **geometric growth** on all $y$-derivatives of the data:

$$
\max_{y\in\Gamma}|\partial_y^k c(y)|\le\gamma^k,
\qquad
\max|\partial_y^k u_R|\le\delta_R^k,
\qquad
\max_{y\in\Gamma}\|\partial_y^k u_0\|_V^2\le\eta^k,
$$

with the constants raised to the power $k$ — that geometric structure in $k$ is what makes the derivative bounds summable into an analytic extension. The resulting theorem gives

$$
\max_{\Gamma}\|\partial_y^k u(\cdot,t,\cdot)\|_V^2
\ \le\ C_k(T)\bigl(\delta_R^k+\delta_L^k+\eta^k\bigr)<+\infty,
$$

and an analytic extension in $\Sigma(\Gamma,\tau)=\{z\in\mathbb C:\mathrm{dist}(z,\Gamma)\le\tau\}$ for $0<\tau<1/\sqrt{\zeta}$. For stochastic collocation this yields

$$
\min_{v}\|u-v\|_{L^\infty[T,C^0(\Gamma,V)]}
\ \le\ \frac{2}{\varrho-1}\,e^{-p\log\varrho}
\max_{z\in\Sigma(\Gamma,\tau)}\|u(z)\|_{L^\infty(T,V)},
\qquad
\varrho=\tau+\sqrt{1+\tau^2}>1 ,
$$

exponential decay in the polynomial degree $p$ at a rate $\log\varrho$ set by the width $\tau$ of the analyticity strip.

> [!warning] An important qualification
> The paper states that unlike random elliptic or parabolic problems, solutions of random **hyperbolic** equations are **not analytic in general** with respect to the random parameters. The complex-analytic sharpening is therefore available only in special cases such as periodic boundary conditions with analytic data. The paper also calls its own stochastic Galerkin algebraic estimate
> $$\|u-u^N_{SG}\|_{L^2(\Gamma,V)}\le C_\Gamma\sqrt{C(T)}\bigl(\sqrt{\zeta}N\bigr)^{-m}$$
> "rather rough". Both qualifications belong with any citation of the exponential rate.

The numerical example takes $u_t=y\,u_x$ with $u(x,0;y)=\cos(y)$ and exact solution $u=\cos(x-yt)$, which lies in $H^{(m)}_y(-1,1)$ for every positive integer $m$, so exponential convergence is expected and is observed for both Galerkin and collocation. The paper also reports errors growing with $t$, the known long-time degradation of polynomial chaos and stochastic collocation.

## 5 and 3: decoupling the Galerkin system with bi-orthogonal polynomials

Paper 5 treats Galerkin methods for random hyperbolic problems and names a threefold difficulty: hyperbolic equations have the worst regularity of the three PDE types, so the parametric solution may be only $BV$ or in a low-order $H^k$; Galerkin projection couples the modes, so a decoupling technique is needed, and unlike collocation boundary conditions cannot be imposed pointwise; and if the random field is a truncated Karhunen-Loève expansion in $N$ variables, an isotropic tensor basis has size exponential in $N$.

The remedy is a **bi-orthogonal (double-orthogonal) polynomial basis**, defined by a pair of conditions:

$$
\int_{\Gamma_i}\rho(y_i)\phi_{j,i}\phi_{k,i}\,\mathrm dy_i=\delta_{jk},
\qquad
\int_{\Gamma_i} y_i\,\rho(y_i)\phi_{j,i}\phi_{k,i}\,\mathrm dy_i=C_{k,i}\,\delta_{jk} .
$$

The second, weighted orthogonality diagonalises the multiplication-by-$y_i$ operator and thereby decouples the Galerkin system into independent deterministic solves. That is how the structure of the matrix $A$ from paper 1 gets used: rather than asking when $A$ is diagonally dominant, change basis so that it is diagonal.

For the tensor-basis size, the paper uses a sensitivity estimate to choose anisotropic polynomial orders: a bound of the form

$$
\int\!\!\int\rho(y)u_{y_i}^2\,\mathrm dx\,\mathrm dy\le C(T)\bigl(\sqrt{\lambda_i}+\sqrt{\mu_i}+\nu_i\bigr)
$$

justifies larger order $r_i$ for small $i$, the directions with larger Karhunen-Loève spectrum.

Paper 3 applies the same machinery to elliptic **interface** problems with random input. Interface elliptic problems are already hard deterministically: standard finite elements lose accuracy unless the mesh fits the interface. Adding a random input makes the naive route doubly expensive, since projection produces a coupled system of deterministic interface problems each needing an interface-fitted mesh. The paper removes both couplings at once: bi-orthogonal polynomials decouple the random direction, and **immersed finite element** spaces let the mesh in the physical direction be independent of the interface geometry.

## 7 and 38: two application directions

- **7 (delay differential equations with random input)** handles delay structure. The delay term brings the history on $[t-\tau,t]$ into the equation, and that history itself depends on the parameters, giving the parameter dependence an extra piecewise character.
- **38 (Maxwell's equations with random inputs)** gives efficient stochastic Galerkin methods. Maxwell's system is a first-order hyperbolic system and therefore inherits all the regularity difficulties of papers 2, 4 and 5, while its divergence constraint adds a further structural requirement.

> [!note] Coverage status
> The theorems and key constructions of papers 1, 2, 4 and 5 have been checked against the sources. Papers 7 and 38 receive only their problem setting and positioning, at the level of abstract and metadata. Paper 3 is the one item that could not be verified: Crossref, OpenAlex, zbMATH and Semantic Scholar all lack its abstract and the publisher page returns 403, so this site reports none of its theorems or numerical results, and the discussion above gives only the constructional idea confirmable from zbMATH indexing keywords and independent sources by the same group.

## Coverage check

| Item                                            | Paper | Status                                                |
| ----------------------------------------------- | ----- | ----------------------------------------------------- |
| Galerkin coefficient matrix and its sparsity    | 1     | $a_{jk}$, $e_{ijk}$, symmetry, $2N+1$ nonzeros        |
| Reduction to three-term recurrence coefficients | 1     | both identities and their technical role              |
| Theorem 1 and the threshold $\|\alpha\|\ge1/2$  | 1     | both densities, conclusion, qualification, convention |
| Sign-changing speed and switching inflow        | 2     | problem setting and the regularity consequence        |
| Data hypotheses weighted by $1/\|c\|$           | 2     | all three integral conditions and the coupling        |
| Geometric growth and analytic extension         | 4     | hypothesis form, derivative bound, analyticity region |
| Exponential convergence for collocation         | 4     | the estimate, relation between $\varrho$ and $\tau$   |
| Hyperbolic solutions not analytic in general    | 4     | both qualifications the paper states                  |
| Both defining conditions of bi-orthogonality    | 5     | definition, decoupling mechanism, relation to paper 1 |
| Sensitivity basis for anisotropic orders        | 5     | the bound and its consequence                         |
| Double decoupling for interface problems        | 3     | remedy in the random and physical directions          |

## Sources for this page

- T. Zhou and T. Tang, [_Note on coefficient matrices from stochastic Galerkin methods for random diffusion equations_](https://doi.org/10.1016/j.jcp.2010.07.016), J. Comput. Phys. 229 (2010), pp. 8225-8230.
- T. Tang and T. Zhou, [_Convergence analysis for stochastic collocation methods to scalar hyperbolic equations with a random wave speed_](https://doi.org/10.4208/cicp.060109.130110a), Commun. Comput. Phys. 8 (2010), pp. 226-248.
- T. Zhou, [_Stochastic Galerkin methods for elliptic interface problems with random input_](https://doi.org/10.1016/j.cam.2011.05.033), J. Comput. Appl. Math. 236 (2011), pp. 782-792.
- T. Zhou and T. Tang, [_Convergence analysis for spectral approximation to a scalar transport equation with a random wave speed_](https://doi.org/10.4208/jcm.1206-m4012), J. Comput. Math. 30 (2012), pp. 643-656.
- T. Zhou and T. Tang, [_Galerkin methods for stochastic hyperbolic problems using bi-orthogonal polynomials_](https://doi.org/10.1007/s10915-011-9508-0), J. Sci. Comput. 51 (2012), pp. 274-292.
- T. Zhou, [_A stochastic collocation method for delay differential equations with random input_](https://doi.org/10.4208/aamm.2012.m38), Adv. Appl. Math. Mech. 6 (2014), pp. 403-418.
- Z. Feng, J. Li, T. Tang, and T. Zhou, [_Efficient stochastic Galerkin methods for Maxwell's equations with random inputs_](https://doi.org/10.1007/s10915-019-00936-z), J. Sci. Comput. 80 (2019), pp. 248-267.
