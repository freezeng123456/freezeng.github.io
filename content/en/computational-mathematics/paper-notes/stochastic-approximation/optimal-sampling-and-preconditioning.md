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

### What the stability factor is replaced by

To put the qualitative claim above quantitatively: for compact $D$, sampling i.i.d. from the equilibrium measure with weights $N/K(z)$, the procedure is stable with high probability provided that, for any $r>0$,

$$
\frac{S}{N\log S}\ \ge\ C\,\frac{1+r}{\lambda_{\min}(R)} .
$$

**The decisive comparison is that the factor $\|K\|_\infty/N$ has been replaced by $1/\lambda_{\min}(R)$**, where $R$ is the relevant Gram matrix. The accompanying accuracy estimate is

$$
\mathbb E\Bigl[\bigl\|f-T_L(\tilde\Pi_Sf)\bigr\|_w^2\Bigr]
\le\|f-\Pi f\|_w^2
+\frac{\varepsilon(S)}{\lambda_{\min}(R)}\|f-\Pi f\|_{\tilde w}^2
+\frac{8L^2}{S^r}
+4\kappa^2(R)\,d^2(f),
$$

with $\varepsilon(S)=\frac{2-2\log2}{(1+r)\log S}\to0$ and $\kappa(R)=\lambda_{\max}(R)/\lambda_{\min}(R)$.

> [!warning] The paper's own reservations about its theory
> This deserves recording faithfully, because it is easily lost in restatement. The last term $4\kappa^2(R)d^2(f)$ **does not vanish as $S\to\infty$**, so taken at face value this theory is **weaker** than the established standard Monte Carlo theory. The paper says plainly that it cannot yet prove the two comparable, and offers two weaker forms of support instead: an argument that its bound is of the same magnitude as the Monte Carlo bounds, and numerical evidence that weighted least squares is **frequently**, not always, superior. It also notes that in one dimension $\lambda_{\min}(R)$ and $\kappa(R)$ are numerically well behaved for general weights.
>
> In other words the case for this route rests here more on experiment than on theorem; it is paper 45's induced-sampling criterion below that turns the linear regime into a clean theorem.

## 45: writing the thread up as a survey

Paper 45, in _SIAM Review_, is the hub of this family. It restates and unifies papers 9, 10 and 11 together with external results (Rauhut, Rauhut-Ward, Yan-Guo-Xiu) and introduces the potential-theoretic vocabulary: Lebesgue constants, Fekete points, equilibrium measures, contraction factors.

Its sample-complexity comparison is the most informative part:

- **i.i.d. sampling from the orthogonality (uniform) measure** (Cohen-Davenport-Leviatan): for any $r>0$, if $M/\log M\ge C_rN^2$ then $\Pr[|\!|\!|\hat A-I|\!|\!|\ge\frac12]\le2M^{-r}$. This is the **quadratic** requirement $M\gtrsim N^2\log N$, and it extends to multidimensional spaces with arbitrary **lower** index sets under the same scaling.
- **Monte Carlo sampling from the Chebyshev measure**: the requirement falls to $M\sim N^{\log3/\log2}$, strictly better than $N^2$ but still superlinear.
- **Deterministic Weil point sets** (from paper 9): $M\ge C(d)N^2$ gives a unique solution and near-best approximation; that quadratic requirement is stronger than Chebyshev Monte Carlo, and the compensating advantage is determinism.

The paper is equally explicit about the gap between theory and practice: practitioners typically take $M\simeq cN$ with $c$ between 2 and 3, that is **linear**, and definitive theory for the linear regime "is not yet definitively available".

### One criterion covering every sampling density

The survey's most important move is to gather the results above into a single criterion. With $q$ the sampling density and $\{v_n\}$ the basis, stability of the normal equations requires only

$$
\frac{M}{\log M}\ \ge\ C(r+1)\,
\sup_{x\in D}\sum_{n=1}^{N}\Bigl(\frac{v_n(x)}{q(x)}\Bigr)^{2},
\qquad
C=\frac{2}{\log(27/8e)}\approx9.24,
$$

and then $\|G-I\|_2\le\frac12$ holds with probability at least $1-2M^{-r}$, by a matrix Chernoff argument.

**The problem thus becomes a pure optimisation: choose $q$ to minimise that supremum.** Its lower bound is exactly $N$, attained when $q\propto\sum_n v_n^2$ — which is **induced sampling**. Substituting gives

$$
\frac{M}{\log M}\ \ge\ C(r+1)\,N,
$$

that is $M\sim N\log N$, optimal in $N$ up to the logarithm.

**This criterion depends on $N=\dim V$ alone.** It does not depend on the dimension $d$, on the domain $D$, on the weight $w$, or even on which $N$-dimensional subspace was chosen. The price is explicit: one must sample from the nonstandard density $\rho$, and $\rho$ **does** depend on $(V,w,D)$. The substance of this whole thread is therefore that **all the problem-dependence has been moved out of the sample complexity and into the sampling density** — the complexity becomes problem-independent, and every difficulty concentrates into the single question of how to sample from $\rho$.

**That passage is the coordinate system for the whole thread.** The linear regime identified in papers 22 and 45 is exactly the target that paper 28 reaches by entirely different means.

## 28: greedy selection instead of random sampling

Paper 28 introduces **weighted approximate Fekete points**. Fekete points maximise the absolute value of the Vandermonde determinant and relate directly to the Lebesgue constant, but computing them exactly is a hard nonconvex problem. Approximate Fekete points are selected greedily from a candidate set by **QR with column pivoting**: the pivot order approximately maximises the determinant.

"Weighted" means the Christoffel weights are folded into the candidate matrix, so selection is optimal in the weighted sense. The gain is that the sample count can be pushed close to $N$; the price is that a candidate pool is needed and the independence structure of random sampling is lost, so the guarantees are no longer of with-high-probability type but rest on deterministic potential-theoretic arguments.

### In one dimension the greedy choice is exactly optimal

This paper has a conclusion far stronger than "greedy is a good heuristic", and it deserves separating out. Take any probability density $\rho$ on $\Gamma=\mathbb R$ with $\Lambda=\{0,\dots,N-1\}$, write $\phi_N$ for the degree-$N$ orthonormal polynomial, and define the meromorphic function

$$
r_N(y)=\frac{\phi_N(y)}{\phi_{N-1}(y)} .
$$

For any $y\notin\phi_{N-1}^{-1}(0)$ the set $A_N(y)=r_N^{-1}\bigl(r_N(y)\bigr)$, taking the set-valued inverse, is **uniquely determined**, and it carries an $N$-point positive quadrature rule exact to degree $2N-2$ whose weights are precisely the **Christoffel weights**:

$$
\int_\Gamma p(z)\rho(z)\,\mathrm dz=\sum_{z\in A_N(y)}\frac{1}{K_\Lambda(z)}\,p(z),
\qquad \deg p\le2N-2 .
$$

In particular $A_N(y)$ is the $N$-point **Gauss** rule when $y\in\phi_N^{-1}(0)$. **The $1/K_\Lambda(z)$ appearing here is the same object as paper 22's weight $N/K(z)$**, so two apparently different routes meet at this point.

On that basis the paper proves that starting from such an $A_N(y)$, greedy selection and globally optimal selection **coincide exactly**, with

$$
\bigl|\det V\bigr|=1=\kappa(V),
$$

so the resulting Vandermonde matrix is **perfectly conditioned**. In one dimension the greedy choice is therefore not an approximation but attains the optimum. More importantly this holds for **any** univariate density, **including densities without compact support** — which is exactly the property approximate Fekete points lack, so it is also a positive answer to the unbounded-domain difficulty.

> [!note] A misprint
> The theorem's text prints the starting condition as $y\in\phi_{N-1}^{-1}(0)$, which contradicts the Lemma 3.1 and Corollary 3.1 it relies on, where $y$ must lie **outside** the zero set. It should read $y\notin\phi_{N-1}^{-1}(0)$. This page states the latter and flags the discrepancy here.

**Papers 22, 45 and 28 therefore form a complete spectrum:**

| Route                            | Sample requirement   | Type of guarantee                  | Extra structure needed      |
| -------------------------------- | -------------------- | ---------------------------------- | --------------------------- |
| sample the orthogonality measure | $M\gtrsim N^2\log N$ | probabilistic                      | none                        |
| Christoffel-weighted sampling    | log-linear           | probabilistic                      | equilibrium measure, $K(z)$ |
| greedy (approximate Fekete)      | close to $N$         | deterministic, potential-theoretic | candidate pool, pivoted QR  |

## 24: designing sampling and preconditioning as a pair

Paper 24 targets **sparse** approximation of polynomial chaos expansions with a generalised sampling and preconditioning scheme. The key recognition is that in the sparse-recovery framework (bounded orthonormal systems and the restricted isometry property) the quantity governing sample complexity is again the uniformity of the row norms, so the "sampling density" and the "preconditioner" should be designed as a pair rather than choosing the sampling first and patching with a preconditioner afterwards.

The paper explicitly advertises the framework as covering "bounded or unbounded" domains, so it also answers the unbounded-domain difficulty that [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 11]] met in the least-squares setting.

### The price of an unbounded domain can be written down

The main recovery theorem requires a sample count of the form

$$
M\ \ge\ L(n)\,\bigl\|R^{-1/2}\bigr\|_1^2\;s\log^3(s)\log(N),
$$

and **the entire difference between bounded and unbounded sits in the factor $L(n)$**:

- bounded, Chebyshev-like: $L(n)=C(\alpha,\beta)$ **uniformly** in $n$, so the requirement is just $M\gtrsim s\log^3(s)\log(N)$, near-linear in the sparsity $s$;
- unbounded: $L(n)=Cn^{\max\{1/\alpha,\,2/3\}}$, for instance $Cn^{2/3}$ when $\alpha\ge\frac32$, so the requirement degrades to $M\gtrsim s\,n^{2/3}$, with an extra factor growing in the polynomial degree.

The conclusion solves the weighted $\ell_1$ problem $\min_\alpha\|R^{1/2}\alpha\|_1$ subject to $\|\sqrt W\Phi\alpha-\sqrt Wf\|_2\le\varepsilon$, and the two terms of the error bound carry $\lambda_{\min}(R)$ and $\|R^{-1/2}\|_1$ as their coefficients — **as in paper 22, it is the spectrum of $R$ that ends up in the bound.**

> [!note] An internal inconsistency
> In the paper's third case, labelled CSA-c, the $\max$ expression reads $\max\{1/(2\alpha),2/3\}$ while the branch expansion immediately below prints $Cn^{1/\alpha}$ rather than $Cn^{1/(2\alpha)}$. This page reproduces the source as printed and flags it here rather than correcting it silently.

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
