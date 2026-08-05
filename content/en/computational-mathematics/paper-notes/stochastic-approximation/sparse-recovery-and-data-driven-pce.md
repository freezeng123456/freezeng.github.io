---
title: Sparse Recovery and Data-Driven Chaos
description: Papers 10, 21, 29, 32, 36 and 44 - folding gradients and empirical measures into the recovery problem
lang: en
translation: computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce
tags:
  - paper-notes
  - uncertainty-quantification
  - sparse-recovery
---

> [!note] Coverage of this page
> Papers **10** (_SIAM J. Sci. Comput._ 36(4), 2014), **21** (_SIAM J. Sci. Comput._ 39(1), 2017), **29** (_Commun. Comput. Phys._ 24, 2018), **32** (_J. Comput. Phys._ 367, 2018), **36** (_J. Comput. Phys._ 381, 2019) and **44** (_Commun. Math. Res._ 36, 2020).

![Fold gradients and sampling density into the recovery problem](assets/diagrams/tao-zhou-papers/en/sparse-recovery.svg)

## Sparse recovery and least squares are governed by the same quantity

In the sparse-recovery framework the design matrix must satisfy the restricted isometry property in a bounded orthonormal system, and the required number of measurements is controlled by the uniformity of the row norms — **exactly the quantity that governs the [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|least-squares case]]**. That is why the work on this page shares the language of Christoffel weights and preconditioning with that page.

## 10: deterministic interpolation points and sparse interpolation

Paper 10 studies sparse interpolation and the design of deterministic interpolation points. The question is: if the target is known to be sparse in some basis, can it be recovered from far fewer points than the dimension of the space, with the point set given deterministically rather than drawn at random?

The technical tools come from the same source as [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 9]] (the Weil-sum machinery), so the conclusions are likewise deterministic rather than probabilistic. Paper **29 (gradient-enhanced ℓ1 recovery of sparse trigonometric polynomials)** applies the same machinery in the trigonometric basis.

## 21: randomly subsampling a Gauss grid, then minimising the ℓ1 norm

Paper 21 is the ℓ1 sibling of [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 13]]: the same randomised Gauss-quadrature subsampling design, with least squares there and sparse recovery here. It generalises the uniform-only result of Tang and Iaccarino to Beta, normal and exponential parameters.

Its central observation matches paper 13: **the Gauss weights are exactly Christoffel function values**, so subsampling a Gauss grid is implicitly Christoffel-weighted sampling. On a tensor Gauss grid $\Theta_{\mathbf n}$ the weights are

$$
w_{\mathbf k}=\lambda_{\mathbf n}(z_{\mathbf k})
=\prod_{i=1}^{d}\frac{1}{\sum_{k=0}^{n_i-1}\bigl[\phi^i_k(z^i_{k_i})\bigr]^2},
$$

and the uniform empirical probability measure on the grid,

$$
\nu_{\mathbf n}=\frac{1}{\prod_i n_i}\sum_{\mathbf k\le\mathbf n}\delta_{z_{\mathbf k}},
$$

makes i.i.d. sampling from $\nu_{\mathbf n}$ equivalent to uniform sampling from the tensor Gauss grid.

The paper further shows that after weighting the design matrix by the Christoffel function the resulting family $\{\psi_{k,\mathbf n}\}$ is **orthonormal under the discrete measure $\nu_{\mathbf n}$**, and the corresponding matrix $D=\bigotimes_iD^i$ with $D^i=(\Sigma^i)^{1/2}\Psi^i$ is orthogonal. That is the concrete realisation of "row normalisation" in a tensor structure.

The paper also quantifies the unbounded-domain difficulty (the order given is of the form $n^{2d/3}$), which is the sparse-recovery manifestation of the same difficulty that [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 11]] met in the least-squares setting.

## 29 and 32: gradient-enhanced ℓ1 recovery

### The idea

One expensive forward solve usually returns the value **and every partial derivative**, for instance through an adjoint method or automatic differentiation. Using only the values, each sample contributes one row of the design matrix; using the gradient as well, the same sample contributes $1+d$ rows while the number of unknowns is unchanged.

The measurement matrix therefore goes from

$$
A_{\text{value}}\in\mathbb R^{M\times N}
\qquad\text{to}\qquad
\begin{bmatrix}A_{\text{value}}\\ A_{\text{grad}}\end{bmatrix}\in\mathbb R^{M(1+d)\times N},
$$

and the recovery problem reads

$$
\min_c\ \|c\|_1
\qquad\text{s.t.}\qquad
\bigl\|W(Ac-b)\bigr\|_2\le\delta,
$$

with $W$ a diagonal weight matrix equalising the row norms. **This turns "information per expensive sample" into an object one can design**: rather than adding samples, extract more rows from the same sample.

### Two things that must be handled

First, **the row norms are no longer uniform**. Value rows and gradient rows have different magnitudes, since gradient rows carry the scale of a derivative, so preconditioning is no longer optional: without weighting, the restricted-isometry constant is dominated by the worst class of rows. That is why paper 32 shares the "design sampling and preconditioning as a pair" stance with [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|paper 24]].

Second, **gradients are not free**. If the gradient requires an extra solve rather than coming as a by-product of an adjoint method, the factor $(1+d)$ in rows must be weighed against that cost. The regime these papers target is the one where gradients are cheaply available.

Paper 29 treats recovery of sparse **trigonometric** polynomials and paper 32 treats sparse approximation of polynomial chaos expansions. Both share the same stacked measurement structure, differing in the basis and the corresponding bounded-orthonormal-system constant.

## 36 and 44: when the input distribution is given only through samples

### The problem

Standard polynomial chaos assumes the input distribution is known analytically, so the orthogonal basis is given by the orthogonal polynomials of that distribution. In the **data-driven** setting the input distribution is given only through finitely many samples: there is no analytic density and therefore no ready-made orthogonal basis.

Paper 36 (data-driven polynomial chaos expansions with a weighted least-squares approximation) and paper 44 (sparse approximation of data-driven expansions by an induced sampling approach) treat this setting. The point is that the orthogonal basis must be built from the **empirical measure**, for instance by Gram-Schmidt or a numerical three-term recurrence, and the sampling design changes accordingly: both induced sampling and Christoffel weighting need a basis, and the basis is now itself estimated.

That introduces an error source absent from the rest of this topic: **the estimation error of the basis**. It compounds with the approximation error and the sampling error, so analysis in the data-driven setting carries one extra layer compared with the known-distribution case.

> [!note] Coverage status
> The construction of paper 21 (Gauss weights as Christoffel values, discrete orthonormality, the orthogonal matrix structure) has been checked. The specific theorems, sample-complexity constants and experiments of papers 10, 29, 32, 36 and 44 have not been verified item by item here; the content above is limited to what abstracts, citation relations and cross-references in neighbouring papers confirm.

## How the six relate

| No. | Composition of the measurements         | Sampling design                 | Origin of the basis   |
| --- | --------------------------------------- | ------------------------------- | --------------------- |
| 10  | values, deterministic points            | Weil-sum point set              | known distribution    |
| 21  | values, subset of a Gauss grid          | implicitly Christoffel weighted | known distribution    |
| 29  | values plus gradients, trigonometric    | preconditioning required        | known distribution    |
| 32  | values plus gradients, polynomial chaos | preconditioning required        | known distribution    |
| 36  | values, weighted least squares          | set by the empirical measure    | **empirical measure** |
| 44  | values, sparse approximation            | induced sampling                | **empirical measure** |

One judgement runs through all six: **sample complexity in sparse recovery is decided by the uniformity of the row norms, so any operation that changes the row structure — adding gradient rows, changing basis, changing the sampling density — must come with a matching preconditioner.** That unifies this page with the [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|optimal sampling page]]: both handle the same quantity, one inside a least-squares framework and one inside an ℓ1 framework.

## Sources for this page

- Z. Xu and T. Zhou, [_On sparse interpolation and the design of deterministic interpolation points_](https://doi.org/10.1137/13094596X), SIAM J. Sci. Comput. 36(4) (2014), pp. A1752-A1769.
- L. Guo, A. Narayan, T. Zhou, and Y. Chen, [_Stochastic collocation methods via ℓ1 minimization using randomized quadratures_](https://doi.org/10.1137/16M1059680), SIAM J. Sci. Comput. 39(1) (2017), pp. A333-A359.
- Z. Xu and T. Zhou, [_A gradient-enhanced ℓ1 approach for the recovery of sparse trigonometric polynomials_](https://doi.org/10.4208/cicp.OA-2018-0006), Commun. Comput. Phys. 24 (2018), pp. 286-308.
- L. Guo, A. Narayan, and T. Zhou, [_A gradient enhanced ℓ1-minimization for sparse approximation of polynomial chaos expansions_](https://doi.org/10.1016/j.jcp.2018.04.026), J. Comput. Phys. 367 (2018), pp. 49-64.
- L. Guo, Y. Liu, and T. Zhou, [_Data-driven polynomial chaos expansions: a weighted least-square approximation_](https://doi.org/10.1016/j.jcp.2018.12.020), J. Comput. Phys. 381 (2019), pp. 129-145.
- L. Guo, A. Narayan, Y. Liu, and T. Zhou, [_Sparse approximation of data-driven polynomial chaos expansions: an induced sampling approach_](https://doi.org/10.4208/cmr.2020-0010), Commun. Math. Res. 36 (2020), pp. 128-153.
