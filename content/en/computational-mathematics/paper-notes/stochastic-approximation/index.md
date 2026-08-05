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

The object at the centre of this family is the Christoffel function. For an orthonormal basis $\{\phi_n\}_{n=1}^{N}$, its discrete form supplies the weight in weighted least squares,

$$
w(y)=\frac{N}{\sum_{n=1}^{N}\phi_n^{2}(y)} ,
$$

and its reciprocal supplies the sampling density. Its role can be stated in one sentence: **it measures how concentrated the polynomial space is at a point, and sampling weighted by it makes the row norms of the design matrix nearly uniform.** With uniform row norms the weighted Gram matrix $G=A^{\mathsf T}W^2A$ is close to the identity, and only then is least squares stable.

That observation improves the sample budget from roughly $M\gtrsim N^2$ for draws from the target measure to roughly $M\gtrsim N\log N$ for draws from a Christoffel-weighted density, while greedy selection (approximate Fekete points, computed by QR with column pivoting) pushes $M$ close to $N$. All three routes share one pipeline; only the sampling density differs.

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

> [!note] Coverage status
> Papers 1, 2, 4, 5, 9, 10, 11, 14, 21, 22, 24, 28, 29, 32, 36, 44 and 45 have been checked equation by equation against full texts, so their formulas, theorem hypotheses and constants are transcribed. Papers 6, 7 and 38 reach only abstract and metadata level. Papers 3 and 13 could not be verified — the publisher blocks the full text of both and no aggregator holds an abstract — so their sections give only what indexing keywords, sister papers and third-party literature confirm, and report none of their theorems, constants or numerical results.
>
> One terminological caution as well: this group of papers is not consistent about the Christoffel function, some writing $\sum_\alpha\varphi_\alpha^2$ and others its reciprocal. This site uses $K(z)=\sum_\alpha\varphi_\alpha^2(z)$ throughout and takes the Christoffel function to be $N/K(z)$, so every weight appearing on these pages is of the form $1/K$. Check each paper's own convention before quoting it.

## Relation to the other topics

This topic supplies the machinery for several later ones. The weighted least squares with weight $w_i=M/\sum_m\Phi_m^2(z_i)$ used by the multi-fidelity polynomial chaos surrogates in the [[en/computational-mathematics/paper-notes/bayesian-inference/index|Bayesian topic]] is exactly the Christoffel weight from here, and "let a sampleable density decide where the collocation points go" in the [[en/computational-mathematics/paper-notes/scientific-machine-learning/index|scientific machine learning topic]] is the same strategy as "draw samples from a Christoffel-weighted density", realised on a different function class.

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]; per-paper references appear at the end of each close-reading page.
