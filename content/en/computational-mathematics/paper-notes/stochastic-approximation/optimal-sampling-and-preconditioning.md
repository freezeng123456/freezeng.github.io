---
title: Optimal Sampling and Preconditioning
description: Papers 22, 24, 28, 45 and 54 - decoupling the sampling density from the orthogonality density
lang: en
translation: computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning
tags:
  - paper-notes
  - uncertainty-quantification
  - optimal-sampling
---

> [!note] Coverage of this page
> Papers **22** (_Math. Comput._ 86, 2017), **24** (_SIAM J. Sci. Comput._ 39(3), 2017), **28** (_SIAM J. Sci. Comput._ 40(1), 2018), **45** (_SIAM Rev._ 62(2), 2020) and **54** (_J. Comput. Phys._ 430, 2021).

![One pipeline for collocation design](assets/diagrams/tao-zhou-papers/en/sampling-design.svg)

## 22: the stability factor is something one can design away

### The exact form of the problem

Standard Monte Carlo least squares draws samples i.i.d. from the orthogonality density $w$. By the Cohen-Davenport-Leviatan analysis the required sample count is controlled by the **stability factor**

$$
\frac{\|K(z)\|_\infty}{N},
$$

where $K$ is the diagonal of the $L^2_w$ reproducing kernel of the polynomial space $\mathbb P_k$: with $\varphi$ the vector of $w$-orthonormal basis functions, $K(z)=\varphi^T\varphi$. The quantity is basis-independent, since any orthogonal change of basis $\psi\leftarrow U\varphi$ leaves $K$ unchanged.

For many weights of interest that factor blows up with the polynomial degree, forcing superlinear sample counts, often far worse than quadratic.

**The insight is that the factor is an artefact of a mismatch.** The density one samples from need not be the density that defines orthogonality. Decoupling the two — sampling from the pluripotential equilibrium measure and reweighting with the Christoffel function — drives the effective stability factor down to its optimal value $N$.

The quantity $N/K_k(z)$ is the (normalised) **Christoffel function**, which gives the method its name.

### The algorithm

On a compact domain $D$ the algorithm is:

1. draw $S$ i.i.d. samples $\{z_s\}$ from the equilibrium measure $\mu_D$;
2. assemble the data vector $u$ with $(u)_s=u(z_s)$;
3. compute least-squares weights $K$ with $(K)_{s,s}=N/K(z_s)$;
4. assemble the $S\times N$ Vandermonde-like matrix $V$ with $(V)_{s,n}=\varphi_n(z_s)$;
5. solve $c=\arg\min_{g}\bigl\|\sqrt{K}Vg-\sqrt{K}u\bigr\|$.

There is a corresponding variant for unbounded domains. The paper notes an equivalent reading of its own: **weighting by the Christoffel function is the same as normalising the rows of the design matrix.** That translates "optimal sampling" into a purely linear-algebraic operation and explains why uniform row norms are the decisive property.

The theoretical basis is the asymptotics of the Christoffel function (Berman, Bloom-Levenberg and others): under suitable pluripotential conditions, $\rho^{2k}(z)K^{(k)}_k(z)\,\mathrm dV(z)/N$ converges weakly to the equilibrium measure $\mathrm d\mu_{D,Q}(z)$. In other words, **in the high-degree limit the Christoffel function finds the right sampling density by itself**, with nothing imposed by hand.

## 45: writing the thread up as a survey

Paper 45, in _SIAM Review_, is the hub of this family. It restates and unifies papers 9, 10 and 11 together with external results (Rauhut, Rauhut-Ward, Yan-Guo-Xiu) and introduces the potential-theoretic vocabulary: Lebesgue constants, Fekete points, equilibrium measures, contraction factors.

Its sample-complexity comparison is the most informative part:

- **i.i.d. sampling from the orthogonality (uniform) measure** (Cohen-Davenport-Leviatan): for any $r>0$, if $M/\log M\ge C_rN^2$ then $\Pr[|\!|\!|\hat A-I|\!|\!|\ge\frac12]\le2M^{-r}$. This is the **quadratic** requirement $M\gtrsim N^2\log N$, and it extends to multidimensional spaces with arbitrary **lower** index sets under the same scaling.
- **Monte Carlo sampling from the Chebyshev measure**: the requirement falls to $M\sim N^{\log3/\log2}$, strictly better than $N^2$ but still superlinear.
- **Deterministic Weil point sets** (from paper 9): $M\ge C(d)N^2$ gives a unique solution and near-best approximation; that quadratic requirement is stronger than Chebyshev Monte Carlo, and the compensating advantage is determinism.

The paper is equally explicit about the gap between theory and practice: practitioners typically take $M\simeq cN$ with $c$ between 2 and 3, that is **linear**, and definitive theory for the linear regime "is not yet definitively available".

**That passage is the coordinate system for the whole thread.** The linear regime identified in papers 22 and 45 is exactly the target that paper 28 reaches by entirely different means.

## 28: greedy selection instead of random sampling

Paper 28 introduces **weighted approximate Fekete points**. Fekete points maximise the absolute value of the Vandermonde determinant and relate directly to the Lebesgue constant, but computing them exactly is a hard nonconvex problem. Approximate Fekete points are selected greedily from a candidate set by **QR with column pivoting**: the pivot order approximately maximises the determinant.

"Weighted" means the Christoffel weights are folded into the candidate matrix, so selection is optimal in the weighted sense. The gain is that the sample count can be pushed close to $N$; the price is that a candidate pool is needed and the independence structure of random sampling is lost, so the guarantees are no longer of with-high-probability type but rest on deterministic potential-theoretic arguments.

**Papers 22, 45 and 28 therefore form a complete spectrum:**

| Route                            | Sample requirement   | Type of guarantee                  | Extra structure needed      |
| -------------------------------- | -------------------- | ---------------------------------- | --------------------------- |
| sample the orthogonality measure | $M\gtrsim N^2\log N$ | probabilistic                      | none                        |
| Christoffel-weighted sampling    | log-linear           | probabilistic                      | equilibrium measure, $K(z)$ |
| greedy (approximate Fekete)      | close to $N$         | deterministic, potential-theoretic | candidate pool, pivoted QR  |

## 24: designing sampling and preconditioning as a pair

Paper 24 targets **sparse** approximation of polynomial chaos expansions with a generalised sampling and preconditioning scheme. The key recognition is that in the sparse-recovery framework (bounded orthonormal systems and the restricted isometry property) the quantity governing sample complexity is again the uniformity of the row norms, so the "sampling density" and the "preconditioner" should be designed as a pair rather than choosing the sampling first and patching with a preconditioner afterwards.

The paper explicitly advertises the framework as covering "bounded or unbounded" domains, so it also answers the unbounded-domain difficulty that [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 11]] met in the least-squares setting.

## 54: carrying the reasoning to reproducing kernel spaces

Paper 54 treats **optimal design for kernel interpolation**. Replacing the polynomial space by a reproducing kernel Hilbert space, the role of the Christoffel function passes to the diagonal of the kernel (the power function), while the principle of "sample weighted by how concentrated the space is at each point" is unchanged. The paper connects that design problem to applications in uncertainty quantification.

**This paper reveals the level of abstraction of the whole thread**: Christoffel-weighted sampling is not really about polynomials but about how concentrated a finite-dimensional function space is at each point. Wherever there is a reproducing kernel, the same design principle applies.

> [!note] Coverage status
> The core construction, the algorithm and the sample-complexity comparison of papers 22 and 45 have been checked. The specific theorems, constants and experiments of papers 24, 28 and 54 have not been verified item by item here; the content above is limited to what abstracts, algorithm descriptions and citation relations confirm.

## Coverage check

| Item                                            | Paper | Status                                                     |
| ----------------------------------------------- | ----- | ---------------------------------------------------------- |
| Stability factor and its basis independence     | 22    | $K(z)=\varphi^T\varphi$, $\|K\|_\infty/N$, invariance      |
| Decoupling sampling from orthogonality          | 22    | the insight and the definition of the Christoffel function |
| Five algorithm steps and row normalisation      | 22    | full steps and the equivalent reading                      |
| Christoffel asymptotics and equilibrium measure | 22    | weak convergence and what it means                         |
| Sample complexity for three sampling routes     | 45    | quadratic, $N^{\log3/\log2}$, deterministic quadratic      |
| The linear gap between theory and practice      | 45    | practice at $M\simeq cN$, state of theory                  |
| Weighted approximate Fekete points, pivoted QR  | 28    | selection mechanism, meaning of weighting, cost and gain   |
| Sampling and preconditioning designed as a pair | 24    | motivation in the sparse setting, scope                    |
| From polynomial spaces to kernel spaces         | 54    | the level of abstraction of the principle                  |

## Sources for this page

- A. Narayan, J. Jakeman, and T. Zhou, [_A Christoffel function weighted least squares algorithm for collocation approximations_](https://doi.org/10.1090/mcom/3192), Math. Comput. 86 (2017), pp. 1913-1947.
- J. Jakeman, A. Narayan, and T. Zhou, [_A generalized sampling and preconditioning scheme for sparse approximation of polynomial chaos expansions_](https://doi.org/10.1137/16M1063885), SIAM J. Sci. Comput. 39(3) (2017), pp. A1114-A1144.
- L. Guo, A. Narayan, L. Yan, and T. Zhou, [_Weighted approximate Fekete points: sampling for least-squares polynomial approximation_](https://doi.org/10.1137/17M1140960), SIAM J. Sci. Comput. 40(1) (2018), pp. A366-A387.
- L. Guo, A. Narayan, and T. Zhou, [_Constructing least-squares polynomial approximations_](https://doi.org/10.1137/18M1234151), SIAM Rev. 62(2) (2020), pp. 483-508.
- A. Narayan, L. Yan, and T. Zhou, [_Optimal design for kernel interpolation: applications to uncertainty quantification_](https://doi.org/10.1016/j.jcp.2020.110094), J. Comput. Phys. 430 (2021), 110094.
