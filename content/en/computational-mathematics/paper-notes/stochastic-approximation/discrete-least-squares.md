---
title: Discrete Least-Squares Approximation
description: Papers 6, 9, 11, 13 and 14 - from "which points" to "how many points before it is stable"
lang: en
translation: computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares
tags:
  - paper-notes
  - uncertainty-quantification
  - least-squares
---

> [!note] Coverage of this page
> Papers **6** (_Commun. Comput. Phys._ 16, 2014), **9** (_SIAM J. Sci. Comput._ 36(5), 2014), **11** (_SIAM J. Sci. Comput._ 36(5), 2014), **13** (_J. Comput. Phys._ 298, 2015) and **14** (_Commun. Comput. Phys._ 18, 2015).

![One pipeline for collocation design](assets/diagrams/tao-zhou-papers/en/sampling-design.svg)

## 6: answering "which points" by experiment first

Paper 6 is the empirical prologue to the whole sampling-design programme. It asks directly how design points for least-squares polynomial approximation should be chosen and answers by experiment. The specific worry it sets out to dispel is that random points spoil the convergence rate — a worry later explained away by the stability theory of papers 9, 22 and 45, which shows that the real issue is not the rate but **how many samples are needed before the design matrix is well conditioned**.

## 9: deterministic point sets and quadratic sample complexity

### The stability theorem

Paper 9 is the first in the sampling-design family to supply theory. It constructs a new type of collocation grid (built from Weil sums) and proves a **deterministic** stability result. With $N=\#\Lambda$ the dimension of the polynomial space and $A$ the corresponding design (Gram) matrix, if

$$
M\ \ge\ 4^{\,d+1}\cdot d^2\cdot N^2
$$

with $M$ prime, then in the spectral norm

$$
\Bigl|\!\Bigl|\!\Bigl|\frac{2^{d+1}}{M}A-I\Bigr|\!\Bigr|\!\Bigr|\le\frac12 .
$$

The proof is short: set $\delta=2^{d}((d-1)\sqrt{M}+1)/M$, apply Gershgorin to get $|\lambda_i-1|\le N\delta$, and check that $N\delta\le\frac12$ under the stated condition.

**This is the quadratic sample-complexity statement, $M\gtrsim N^2$, with an explicit dimension-dependent prefactor $4^{d+1}d^2$.** It is the benchmark that later papers attack.

Uniqueness follows immediately ($A$ is strictly diagonally dominant hence nonsingular), and so does optimal convergence: with $P^\Lambda f$ the best $L^2$ approximation in the Chebyshev measure and $P^\Lambda_m f$ the discrete least-squares solution on this point set,

$$
\|f-P^\Lambda_m f\|_{L^2_{\rho_c}}\ \le\ \Bigl(1+\frac{4}{d^2 N}\Bigr)\,\|f-P^\Lambda f\|_{L^\infty} .
$$

The factor $1+4/(d^2N)$ tends to $1$ as $N$ grows, so the discrete projection is asymptotically as good as the best $L^\infty$ approximation. **The estimate is deterministic, with no "with high probability".** That is what the quadratic sample count buys relative to random sampling.

### Extending to other measures

The paper introduces a domination condition: if there is a constant $C$ independent of $Y$ with $0<\rho(Y)\le C\rho_c(Y)$ for the Chebyshev density $\rho_c$, then

$$
\|f-P^\Lambda_m f\|_{L^2_\rho}\ \le\ \sqrt{C}\,\Bigl(1+\frac{4}{d^2 N}\Bigr)\|f-P^\Lambda f\|_{L^\infty} .
$$

This covers the uniform measure and every measure with $0<\rho_{\min}\le\rho\le\rho_{\max}$, and the paper notes it is useful for **epistemic** uncertainty: when the density of $Y$ is unknown, Chebyshev-based approximation remains efficient as long as the unknown density satisfies the domination condition.

**This paper establishes the template for the family**: an explicit point set or density, a stability theorem quantifying the required sample count, and a near-best-approximation corollary.

## 11: on unbounded domains the situation is qualitatively worse

### The problem

All existing least-squares stability theory concerned **bounded** parameter domains — uniform or Chebyshev measures on $[-1,1]^d$, where $M\sim N^2$ samples suffice and fewer for Chebyshev. Gaussian and Gamma parameters live on unbounded domains, and the paper shows the situation there is worse **in kind**, not merely in degree: with the natural Hermite or Laguerre polynomial chaos basis and samples drawn from the Gaussian or Gamma measure, the design matrix is well conditioned only when the sample count is **exponential** in the dimension of the approximation space.

A second, independent problem is the poor resolution of Hermite expansions: the paper quotes Gottlieb and Orszag's remark that resolving $M$ wavelengths of $\sin(x)$ needs nearly $M^2$ Hermite polynomials.

### The remedy and its place in the family

The paper changes both the **basis** (from polynomials to functions) and the **sampling density** (from Gaussian to a mapped uniform law), obtaining the first **log-linear** rather than quadratic sample requirement in this topic.

**This route and the later Christoffel weighting are two distinct choices**: paper 11 changes basis and density, while papers 22 and 45 keep the polynomial basis and instead change the density and add weights. Both reach log-linear complexity at different prices — the former gives up the orthogonal structure of the polynomial basis, the latter needs the Christoffel function computed.

## 13 and 14: randomised quadrature and unstructured meshes

### 13: Gauss weights are Christoffel function values

Paper 13 randomly subsamples a tensor Gauss grid. Its key observation is that the Gauss weights are exactly Christoffel function values:

$$
w_{\mathbf k}=\lambda_{\mathbf n}(z_{\mathbf k})
=\prod_{i=1}^{d}\lambda^i_{n_i}(z^i_{k_i})
=\prod_{i=1}^{d}\frac{1}{\sum_{k=0}^{n_i-1}\bigl[\phi^i_k(z^i_{k_i})\bigr]^2} .
$$

**Subsampling a Gauss grid is therefore implicitly Christoffel-weighted sampling.** That observation is the direct ancestor of papers 22 and 45, where the sampling density is taken proportional to the inverse Christoffel function rather than inherited from a quadrature rule.

As characterised in the third-party literature, the method is stable with the sample count growing **linearly** in the polynomial dimension — stronger than the $M\gtrsim N\log N$ of induced sampling and the $M\gtrsim N^2$ of plain Monte Carlo. But that is within the restricted design space of subsets of a tensor Gauss grid, and the precise hypotheses and constants are unverified here.

### 14: unstructured multivariate meshes

Paper 14 handles parameter domains without tensor-product structure, performing stochastic collocation on unstructured multivariate meshes. At 36 pages it is the most survey-like member of this family, and what it supplies is a framework in which the earlier constructions can be compared side by side.

## Sample complexity across the five papers

| No. | Point set or density                                | Sample requirement               | Type of result                       |
| --- | --------------------------------------------------- | -------------------------------- | ------------------------------------ |
| 6   | several candidate designs (compared experimentally) | not applicable                   | empirical                            |
| 9   | Weil-sum deterministic point set                    | $M\ge4^{d+1}d^2N^2$ (prime)      | deterministic                        |
| 11  | mapped uniform density plus a function basis        | log-linear                       | probabilistic (limited verification) |
| 13  | random subset of a tensor Gauss grid                | linear (restricted design space) | probabilistic (limited verification) |
| 14  | unstructured multivariate meshes                    | see the source                   | survey and construction              |

One judgement runs through all five: **the right form of the question "which points" is "how many points, drawn from which density".** Paper 6 asks the former; from paper 9 onward it is replaced by the latter, and the latter admits theorems. That reformulation is why the thread was able to advance.

> [!note] Coverage status
> The theorems and sample-complexity statements of papers 9 and 11 have been checked. The specific experiments, hypotheses and constants of papers 6, 13 and 14 have not been verified item by item here; the content above is limited to what abstracts, third-party literature and citation relations confirm.

## Sources for this page

- Z. Gao and T. Zhou, [_On the choice of design points for least square polynomial approximations with application to uncertainty quantification_](https://doi.org/10.4208/cicp.130813.060214a), Commun. Comput. Phys. 16 (2014), pp. 365-381.
- T. Zhou, A. Narayan, and Z. Xu, [_Multivariate discrete least-squares approximations with a new type of collocation grid_](https://doi.org/10.1137/130950434), SIAM J. Sci. Comput. 36(5) (2014), pp. A2401-A2422.
- T. Tang and T. Zhou, [_On discrete least-squares projection in unbounded domain with random evaluations and its application to parametric uncertainty quantification_](https://doi.org/10.1137/140961894), SIAM J. Sci. Comput. 36(5) (2014), pp. A2272-A2295.
- T. Zhou, A. Narayan, and D. Xiu, [_Weighted discrete least-squares polynomial approximation using randomized quadratures_](https://doi.org/10.1016/j.jcp.2015.06.042), J. Comput. Phys. 298 (2015), pp. 787-800.
- A. Narayan and T. Zhou, [_Stochastic collocation on unstructured multivariate meshes_](https://doi.org/10.4208/cicp.020215.070515a), Commun. Comput. Phys. 18 (2015), pp. 1-36.
