---
title: Stochastic Approximation and Collocation Design
description: Twenty-three papers that reduce construction problems in uncertainty quantification to one weighted sampling problem
lang: en
translation: computational-mathematics/paper-notes/stochastic-approximation
tags:
  - computational-mathematics
  - paper-notes
  - uncertainty-quantification
---

This is the earliest of the seven topics, starting in 2010, and the technical source of several later ones. Its 23 papers answer the most basic question in uncertainty quantification: **given an expensive parameter-to-solution map, build a stable polynomial approximation of it from as few forward solves as possible.**

![One pipeline for collocation design](assets/diagrams/tao-zhou-papers/en/sampling-design.svg)

## From intrusive to non-intrusive

The first five papers (2010-2012) treat **intrusive** methods: expand the solution in a generalised polynomial chaos basis and apply Galerkin projection. That route turns a random PDE into a coupled deterministic system, and the problem becomes the algebraic properties and convergence of that system.

Paper 1, for instance, studies exactly the coefficient matrix of that coupled system. With a random field depending affinely on the parameters, $\kappa(x,y)=\kappa_0(x)+\sum_{i}\kappa_i(x)y_i$, Galerkin projection gives $\partial_t v=\nabla\cdot(A\nabla_x v)+f$ with

$$
a_{jk}=\sum_{i=0}^{N}\kappa_i(x)\,e_{ijk},
\qquad
e_{ijk}=\int y_i\,\Phi_j(y)\Phi_k(y)\rho(y)\,\mathrm dy .
$$

Whether $A(x)$ is strictly diagonally dominant decides whether diagonal/off-diagonal splitting solvers — mixed explicit-implicit time stepping, Jacobi iteration, decoupled preconditioned conjugate gradients — are stable and converge at a mesh-independent rate. Xiu and Shen proved dominance only for **symmetric** Beta densities and for Gamma densities, and left the general asymmetric case as an open question. Paper 1 reduces the question to an inequality on Jacobi three-term recurrence coefficients and answers it affirmatively but conditionally: for $\rho(y_i)=(1-y_i)^{\alpha}(1+y_i)^{\beta}$, provided $|\alpha|\ge\tfrac12$ and $|\beta|\ge\tfrac12$,

$$
a_{jj}\ \ge\ \kappa_{\min}+\sum_{k\ne j}|a_{jk}|,
\qquad 1\le j\le M,\ \forall x\in\Omega .
$$

The paper does not claim dominance for all $\alpha\ne\beta$, and whether the threshold can be removed is not settled there.

From papers 6 and 9 onward the route turns **non-intrusive**: solve at a set of parameter points and fit by discrete least squares or sparse recovery. That turn replaces "algebraic properties" by "sampling design", and sampling design is the genuine technical core of this topic.

## Four close readings

| Close reading                                                                                                                                    | Papers                 | Technical core                                |
| ------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------- | --------------------------------------------- |
| [[en/computational-mathematics/paper-notes/stochastic-approximation/stochastic-galerkin-and-collocation\|Stochastic Galerkin and collocation]]   | 1, 2, 3, 4, 5, 7, 38   | algebra and convergence of the coupled system |
| [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares\|Discrete least-squares approximation]]               | 6, 9, 11, 13, 14       | collocation grids and randomised quadrature   |
| [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning\|Optimal sampling and preconditioning]] | 22, 24, 28, 45, 54     | Christoffel weights and greedy selection      |
| [[en/computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce\|Sparse recovery and data-driven chaos]] | 10, 21, 29, 32, 36, 44 | ℓ1 recovery and empirical measures            |

## The central object: the Christoffel function

These 23 papers belong to one topic because they all come down to one quantity. Set out that quantity and the internal order of the list becomes visible.

### What it is

Let $w$ be a weight (probability density) on $D$, let $\{\varphi_\alpha\}$ be orthonormal against $w$, let $\Lambda$ be a finite multi-index set with $N=|\Lambda|$, and let $\mathbb P_\Lambda=\mathrm{span}\{\varphi_\alpha:\alpha\in\Lambda\}$. The reproducing kernel diagonal is

$$
K(z)=K_\Lambda(z)=\sum_{\alpha\in\Lambda}\varphi_\alpha^{2}(z) ,
$$

the (normalised) Christoffel function is its reciprocal $N/K(z)$, and every weight in weighted least squares takes the form $1/K$. Both conventions circulate in the literature — some authors call $\sum_\alpha\varphi_\alpha^2$ itself the Christoffel function, others its reciprocal; the convention below is the former for $K$ and the latter for the Christoffel function, so check each paper's own convention before quoting it.

In vector form $K(z)=\varphi^{\mathsf T}\varphi$, so any orthogonal change of basis $\psi\leftarrow U\varphi$ leaves it unchanged: **$K$ is a property of the subspace $\mathbb P_\Lambda$, not of the basis.** Only a basis-independent quantity deserves to be the centre of a family of methods.

### Why it is the one

Draw $M$ points i.i.d. from a density $\rho=q^2w$ and form the weighted design and Gram matrices

$$
(A)_{m,n}=\frac{\varphi_n(z_m)}{\sqrt{M\,q^2(z_m)}} ,
\qquad
G=A^{\mathsf T}A ,
\qquad
\mathbb E\,G=I .
$$

The squared norm of row $m$ of $A$ is $\frac1M\sum_n(\varphi_n(z_m)/q(z_m))^2$, and $\mathbb E\,\mathrm{tr}\,G=N$, so these squared row norms **sum to $N$ on average**. The stability theorem asks precisely that none of them be an outlier: if

$$
\frac{M}{\log M}\ \ge\ C(r+1)\,\sup_{z\in D}\sum_{n=1}^{N}\Bigl(\frac{\varphi_n(z)}{q(z)}\Bigr)^{2},
\qquad C=\frac{2}{\log(27/8e)}\approx 9.24 ,
$$

then with probability at least $1-2M^{-r}$ we have $\|G-I\|_2\le\frac12$. The supremum on the right is the **stability factor**; its smallest possible value is $N$, and taking $q^2=K/N$ — that is, sampling from the **induced density**

$$
\rho(z)=\frac{1}{N}\sum_{n=1}^{N}\varphi_n^{2}(z)\,w(z)
$$

and weighting by $1/q^2=w/\rho$ — attains that lower bound **exactly**. In one sentence: **the Christoffel function measures how concentrated the polynomial space is at a point, and sampling weighted by it flattens the row norms of the design matrix; only with flat row norms is $G$ close to the identity and least squares stable.**

The same quantity arrives from a second direction: the weights of Gauss quadrature are exactly Christoffel function values. So "randomly subsample a tensor-product Gauss grid" is not an ad hoc device — it is Christoffel-weighted sampling in disguise, which is the construction of papers 13 and 21, treated on [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|Discrete least-squares approximation]].

### Three sampling routes

Once the criterion "$M/\log M$ beats the stability factor" is written down, the topic splits into three competing routes that share one pipeline and differ only in which density the points come from and what weight goes with it:

| Route                                   | Sampling density                                       | Weight            | Sample requirement                                                            |
| --------------------------------------- | ------------------------------------------------------ | ----------------- | ----------------------------------------------------------------------------- |
| plain Monte Carlo                       | the orthogonality measure $w$                          | none              | $M\sim N^2$ (tensor Legendre); $N^{\ln3/\ln2}\approx N^{1.585}$ for Chebyshev |
| Christoffel-weighted (induced) sampling | $\rho=(K/N)\,w$                                        | of the form $1/K$ | $M\gtrsim N\log N$, non-asymptotic, criterion depends on $N$ alone            |
| greedy deterministic selection          | chosen from a candidate set by QR with column pivoting | weighted          | pushes $M$ close to $N$                                                       |

The criterion on the second row depends only on $N=\dim\mathbb P_\Lambda$ and not on the dimension $d$, the domain $D$, the weight $w$, or which $N$-dimensional subspace is used; the price is that one must be able to sample from $\rho$, which does depend on all of those. The third row trades randomness for a selection procedure, and its details are on [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|Optimal sampling and preconditioning]].

Three qualifications, all of them made by the papers themselves:

- Using the **equilibrium measure** — the $k\to\infty$ limit of the induced measure — in place of the induced measure also gives log-linear counts, but only **asymptotically in the degree**, whereas the induced measure achieves it at every finite $N$. Paper 44 notes that at $K=20$ the two designs are still visibly different.
- **The advantage of induced sampling shrinks in high dimension at low degree**, because a low-degree space makes $\rho$ close to $w$; paper 45 finds only a modest advantage on its four-dimensional hyperbolic-cross PDE example. Induced sampling is consistently among the best and is the only route with a minimal-sample-count theorem, but it is not always the numerical winner.
- On a bounded domain the large-$N$ limit of the induced measure is the (tensorised) Chebyshev density, and that is a theorem. **The Gaussian counterpart is only a conjecture** in papers 22, 28, 36 and 45, never proved, and must not be cited as a theorem.

## Core idea of each paper

### Stochastic Galerkin and collocation

- **1** answers Xiu and Shen's open question on diagonal dominance of the Galerkin coefficient matrix, giving an explicit threshold condition for asymmetric Beta densities.
- **2** analyses stochastic collocation for a scalar hyperbolic equation with a random wave speed. The difficulty is that $c(y)$ can change sign, so the inflow boundary switches ends of the domain, and the solution's regularity **in the random variable** is low even for smooth data. The technical device is a set of weighted data hypotheses dividing the density by the wave speed, for example

  $$
  \int_0^T\!\!\int_{\Gamma^-}\frac{\rho(y)}{|c(y)|}\bigl(\partial_y u_L(t;y)\bigr)^2\,\mathrm dy\,\mathrm dt<\infty,
  $$

  and the $1/|c|$ weight is exactly what couples the wave speed to the admissible boundary data. Compared with the Galerkin analysis of Gottlieb and Xiu, the convergence rates here follow from hypotheses on the **data** rather than from an assumed asymptotic decay of the expansion coefficients.

- **3** treats stochastic Galerkin methods for elliptic interface problems with random input.
- **4** gives a convergence analysis for spectral approximation of a scalar transport equation with random wave speed.
- **5** uses bi-orthogonal polynomials for Galerkin methods on stochastic hyperbolic problems.
- **7** treats delay differential equations with random input, where the delay structure gives the parameter dependence an extra piecewise character.
- **38** gives efficient stochastic Galerkin methods for Maxwell's equations with random inputs.

### Discrete least-squares approximation

- **6** studies the choice of design points for least-squares polynomial approximation directly, and is where this topic turns from intrusive methods to sampling design.
- **9** introduces a new type of collocation grid for multivariate discrete least-squares approximation.
- **11** treats discrete least-squares projection with random evaluations on an unbounded domain and applies it to parametric uncertainty quantification. On an unbounded domain, "from which density do we sample" turns from a technical detail into the decisive question.
- **13** builds weighted discrete least-squares polynomial approximations using randomised quadratures.
- **14** performs stochastic collocation on unstructured multivariate meshes, for parameter domains without tensor-product structure.

### Optimal sampling and preconditioning

- **22** introduces the Christoffel function weighted least-squares algorithm and is the theoretical centre of this line.
- **24** gives a generalised sampling and preconditioning scheme for sparse approximation of polynomial chaos expansions, designing the sampling density and the preconditioner as a pair.
- **28** introduces weighted approximate Fekete points, using greedy selection to push the sample count close to the polynomial dimension.
- **45** is the _SIAM Review_ survey organising the construction of least-squares polynomial approximations, and the recommended entry point to this literature.
- **54** treats optimal design for kernel interpolation, carrying the same reasoning from polynomial spaces to reproducing kernel spaces.

### Sparse recovery and data-driven chaos

- **10** studies sparse interpolation and the design of deterministic interpolation points.
- **21** performs stochastic collocation by ℓ1 minimisation using randomised quadratures.
- **29** and **32** introduce gradient-enhanced ℓ1 recovery: one forward solve returns the value and every partial derivative, so the number of measurement rows grows by a factor $(1+d)$ with no new unknowns. That turns "information per expensive sample" into an object one can design.
- **36** and **44** treat data-driven polynomial chaos: when the input distribution is given only through samples with no analytic form, the orthogonal basis must be built from an empirical measure and the sampling design changes accordingly.

![Fold gradients and sampling density into the recovery problem](assets/diagrams/tao-zhou-papers/en/sparse-recovery.svg)

## Relation to the other topics

This topic supplies the machinery for several later ones. The weighted least squares with weight $w_i=M/\sum_m\Phi_m^2(z_i)$ used by the multi-fidelity polynomial chaos surrogates in the [[en/computational-mathematics/paper-notes/bayesian-inference/index|Bayesian topic]] is exactly the Christoffel weight from here, and "let a sampleable density decide where the collocation points go" in the [[en/computational-mathematics/paper-notes/scientific-machine-learning/index|scientific machine learning topic]] is the same strategy as "draw samples from a Christoffel-weighted density", realised on a different function class.

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]; per-paper references appear at the end of each close-reading page.
