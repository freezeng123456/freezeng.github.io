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

## Shared machinery: row norms, the Christoffel function, and Gauss weights

### The idea: stability depends only on how uniform the row norms of the design matrix are

On the surface these five papers do five different things — compare candidate point sets, construct a deterministic grid, change both basis and density, subsample a Gauss grid, collocate on unstructured meshes — but one quantity governs all of them. State it first and every later trade-off becomes obvious.

Let $\{\varphi_n\}_{n=1}^{N}$ be an orthonormal basis of $L^2_w$, let $q>0$ be a biasing function with $\|q\|_{L^2_w}=1$, and draw $x_1,\dots,x_M$ i.i.d. from the density $\rho:=q^2w$. The weighted design matrix and right-hand side are

$$
(A)_{m,n}=\frac{1}{\sqrt{M\,q^2(x_m)}}\,\varphi_n(x_m),
\qquad
(f)_m=\frac{1}{\sqrt{M\,q^2(x_m)}}\,f(x_m),
$$

and the Gram matrix $G=A^{\mathsf T}A$ satisfies $\mathbb E\,G=I$, with $G\to I$ almost surely for every admissible $q$. Taking $q\equiv1$ recovers the naive scheme: sample from the orthogonality measure, do not weight.

The point is that the squared norm of the $m$th **row** of $A$ is exactly

$$
\|A_{m,\cdot}\|_2^2=\frac1M\sum_{n=1}^{N}\Bigl(\frac{\varphi_n(x_m)}{q(x_m)}\Bigr)^2 ,
$$

and $\mathbb E\,\mathrm{tr}\,G=N$, so these $M$ squared row norms **sum to $N$ on average**. What the stability theorem asks is precisely that none of them be an outlier: if for some $r>0$

$$
\frac{M}{\log M}\ \ge\ C(r+1)\,\sup_{x\in D}\sum_{n=1}^{N}\Bigl(\frac{\varphi_n(x)}{q(x)}\Bigr)^2 ,
\qquad C=\frac{2}{\log(27/8e)}\approx 9.24 ,
$$

then with probability at least $1-2M^{-r}$ we have $\|G-I\|_2\le\frac12$. The supremum on the right is the **stability factor**, and its smallest possible value is $N$. The whole topic therefore reduces to one sentence: **the stability theorem asks how much worse the worst row is than the average row, and every sampling-design trick is an attempt to squeeze that ratio.** With row norms equal everywhere the stability factor attains its lower bound $N$ and the budget is $M\gtrsim N\log N$; a spike anywhere multiplies the budget by its height.

### The Christoffel function and this site's convention

When $q\equiv1$ the stability factor is the supremum of

$$
K(z)=K_\Lambda(z)=\sum_{\alpha\in\Lambda}\varphi_\alpha^2(z) .
$$

Writing it as $K(z)=\varphi^{\mathsf T}\varphi$ shows that any orthogonal change of basis $\psi\leftarrow U\varphi$ leaves it unchanged: **$K$ is a property of the subspace $\mathbb P_\Lambda$, not of the basis**. That is why it deserves to be the central object.

Conversely, take $q^2=K/N$, that is, sample from

$$
\rho(x)=\frac{1}{N}\sum_{n=1}^{N}\varphi_n^2(x)\,w(x)
$$

and weight by $1/q^2=w/\rho$. Then the stability factor is **identically** $N$ — the row norms have been flattened by construction. That $\rho$ really is a probability density follows from $\int_D\frac1N\sum_n\varphi_n^2\,w\,\mathrm dx=\frac1N\sum_n\|\varphi_n\|^2=1$. This is induced sampling, the content of the other branch of the topic (papers 22 and 45).

> [!warning] Convention
> This site uses $K(z)=\sum_\alpha\varphi_\alpha^2(z)$ throughout, takes the (normalised) Christoffel function to be $N/K(z)$, and therefore writes every weight on this page in the form $1/K$. The $\lambda_\Lambda(z)$ of the sources is $1/K_\Lambda(z)$. The papers themselves are not consistent: papers 22 and 28 write $K$ for $\sum_\alpha\varphi_\alpha^2$, the *reciprocal* of the Christoffel function; paper 24 writes $\lambda_\Lambda$ for $1/\sum\varphi_i^2$, the Christoffel function itself; paper 36 writes $K(\xi)$ for $N/\sum_j\Phi_j^2$; papers 44 and 45 write $\kappa$ and $q^2$ for the normalised reciprocal. **Check each paper's own convention before quoting it**, or the weight comes out inverted.

### Gauss weights are Christoffel function values

This is the fact to state before any formalism, because papers 13 and 14 are built entirely on it. In one dimension with $\Lambda=\{0,\dots,N-1\}$ there is an $N$-point positive quadrature rule, exact for polynomials of degree at most $2N-2$, whose weights are precisely Christoffel function values:

$$
\int_\Gamma p(z)\rho(z)\,\mathrm dz=\sum_{z\in A_N}\frac{1}{K_\Lambda(z)}\,p(z),
\qquad \deg p\le 2N-2 ,
$$

with the Gauss abscissae the special case $y\in\phi_N^{-1}(0)$. **So "randomly subsample a Gauss grid and weight by the Gauss weights" is not an ad hoc reweighting: it is Christoffel-weighted sampling.** Papers 13 and 21 exploit exactly this, and papers 22 and 45 later make the implicit choice explicit by taking the sampling density to be $K/N$ times the orthogonality measure.

### Three sampling routes and their budgets

Substituting particular bases into the criterion $M/\log M\gtrsim\lVert K\rVert_\infty$ produces the three budgets that recur throughout the topic:

| Sampling scheme                             | Growth of the stability factor            | Sample requirement            |
| ------------------------------------------- | ----------------------------------------- | ----------------------------- |
| uniform measure, tensor Legendre, unweighted | $\lVert K\rVert_\infty\sim N^2$           | $M\sim N^2$                   |
| Chebyshev measure, unweighted               | $\lVert K\rVert_\infty\sim N^{\ln3/\ln2}$ | $M\sim N^{1.585}$             |
| Christoffel-weighted (induced) sampling     | identically $N$                           | $M\gtrsim N\log N$, non-asymptotic |

The criterion on the third row depends only on $N=\dim\mathbb P_\Lambda$ and is independent of the dimension $d$, the domain $D$, the weight $w$, and which $N$-dimensional subspace is chosen; the price is that one must be able to sample from $\rho=q^2w$, which does depend on $(V,w,D)$. Greedy deterministic selection (approximate Fekete points) is a fourth route that pushes $M$ further towards $N$, and belongs to [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|Optimal sampling and preconditioning]].

> [!warning] The Gaussian asymptotic density is a conjecture
> On a bounded domain the large-$N$ limit of the induced measure is the (tensorised) Chebyshev density, and that is a theorem. **The Gaussian case is not.** For $D=\mathbb R^d$ with the Gaussian weight, $\lim_{k\to\infty}\rho(x/\sqrt k)=C(2-\|x\|_2^2)^{d/2}$ appears in papers 22, 28, 36 and 45 strictly as a **conjecture** and has never been proved; paper 45 says outright that in the multivariate case "much more is unknown, and in many cases we currently have only conjectures". Neither this page nor its neighbours cite it as an established result.

## 6: answering "which points" by experiment first

### The idea

Paper 6 is the empirical prologue to the whole sampling-design programme. The worry at the time was concrete: **do random points spoil the convergence rate?** The worry has intuitive backing — a structured grid (a sparse grid) is designed, whereas random points may clump and leave gaps, so they ought to be worse. The paper does not answer theoretically. It splits "which points are better" into three criteria that are usually studied separately and measures all three at once: convergence rate, stability (the condition number of the design matrix), and robustness when the function values carry numerical noise.

Its contribution is not a mathematical object but **the comparison protocol itself**. Later work showed the framing was half right: the real issue is not the rate but how many samples are needed before the design matrix is well conditioned — that is, the stability factor of the previous section.

### Setting

Discrete least squares on polynomial spaces, with three candidate families of design points: Sparse Grid (SG), Monte Carlo (MC) and Quasi Monte Carlo (QMC). The "design matrix" of the paper should in standard form be $A_{ij}=\phi_j(y^{(i)})$, the Vandermonde-like matrix of basis functions evaluated at the design points — **the abstract refers to "the design matrix" without reproducing it, and this form is supplied here from context**. It differs from the weighted form used by the later papers (previous section), which carries an extra $\sqrt{w}$ row scaling.

### Numerical evidence

The tests are several classical high-dimensional test functions together with a random ODE model. Three findings:

| Points | Convergence                                                                       | Stability                                     |
| ------ | --------------------------------------------------------------------------------- | --------------------------------------------- |
| MC     | introduces no low convergence rate; high order holds given regularity and enough points | —                                             |
| QMC    | the same, and a good choice in higher dimension                                   | deterministic, and better on conditioning too |
| SG     | better convergence only in very low dimension (the paper says $d\le2$)            | —                                             |

The qualifier "enough points" stays qualitative: no relation between $M$ and $N$ is given anywhere — **exactly the gap that papers 9, 11, 13, 22, 28 and 45 later close**. One more thing is worth noting: of the three criteria the paper sets out, the third (robustness to noise in the function values) has no matching entry among the three findings the abstract reports, so what it concluded there is not known here.

> [!note] This paper reaches abstract level only
> The full text was not obtained here. Specific dimensions beyond the abstract's $d\le2$, error magnitudes, condition-number values, and the configuration behind each finding are all unverified, so this page reports no numbers for it and invents neither theorems nor experiment tables.

### Relation to the others

It asks "which points" and answers by experiment. Everything afterwards answers a sharpened version of the same question, and does so with **constructed** rather than **selected** point sets: the new collocation grid of paper 9, random evaluations in unbounded domains in paper 11, randomised quadratures in papers 13 and 21, unstructured meshes in paper 14, and the Christoffel-weighted sampling and weighted approximate Fekete points of [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|the neighbouring page]].

## 9: deterministic point sets and quadratic sample complexity

### The idea

The least-squares stability theory that preceded this paper (Migliorati–Nobile–von Schwerin–Tempone; Cohen–Davenport–Leviatan) sampled i.i.d. from the orthogonality measure, so its conclusions were inescapably probabilistic: they hold "with high probability" and guarantee nothing for **any particular realisation**. This paper asks whether the qualifier can be removed.

The mechanism deserves a sentence of its own. A concentration inequality says that a random sum of $M$ terms is $O(\sqrt M)$ with high probability. **Square-root cancellation** in number theory says that a particular structured exponential sum **never** exceeds $(d-1)\sqrt M$. The second is stronger because it has no exceptional set. Applying Weil's exponential-sum bound to the off-diagonal entries of the Gram matrix yields not spectral concentration but a **deterministic entrywise** bound, hence strict diagonal dominance — a stronger statement than "the spectral norm is close to $1$ with high probability", and one that holds for every realisation.

### Setting

$Y=(Y_1,\dots,Y_d)^{\mathsf T}$ are mutually independent on $\Gamma\equiv[-1,1]^d$ with marginals $\rho_i$ and joint density $\rho(Y)=\prod_{i=1}^{d}\rho_i(Y_i)$; $f:\Gamma\to\mathbb R$ is approximated in $L^2_\rho$ with $\|f\|_{L^2_\rho}=\bigl(\int_\Gamma f^2\rho\,\mathrm dY\bigr)^{1/2}$; the polynomial space is $\mathbb P_\Lambda=\mathrm{span}\{\Phi_n\}_{n\in\Lambda}$ in tensor-product (TP) or total-degree (TD) form, with $N=\#\Lambda$.

**The new grid (Weil points).** For a prime $M>2q+1$,

$$
\Theta_M:=\Bigl\{\,y_j=\cos(x_j)\ :\ x_j=2\pi\bigl(j,\,j^2,\dots,j^d\bigr)/M,\quad j=0,\dots,\lfloor M/2\rfloor\,\Bigr\},
$$

so the $q$th coordinate of $x_j$ is $2\pi j^q/M$. The number of points is $m+1$ with $m=\lfloor M/2\rfloor$, and the paper notes that $\{y_j\}_{j=0}^{m}$ coincides with $\{y_j\}_{j=m+1}^{M}$, so half of them suffice. The construction is a cosine transplant of a Weyl-sum node set used earlier by Xu for deterministic sampling of sparse trigonometric polynomials and extended to sparse high-dimensional Chebyshev polynomials.

The design matrix uses the unweighted discrete inner product:

$$
A=\bigl(\langle\Phi_i,\Phi_j\rangle_m\bigr)_{1\le i,j\le N},
\qquad
\langle u,v\rangle_m=\sum_{k=0}^{m}u(y_k)v(y_k) .
$$

### Derivation

**Step one: Weil's formula (Theorem 2.1).** Let $M$ be prime and $f(x)=m_1x+m_2x^2+\cdots+m_dx^d$. If there is a $j$ with $1\le j\le d$ and $M\nmid m_j$, then

$$
\Bigl|\sum_{j=0}^{M-1}e^{\frac{2\pi i f(j)}{M}}\Bigr|\ \le\ (d-1)\sqrt{M} .
$$

This is where the cosine transplant earns its place: $y_j=\cos(x_j)$ turns Chebyshev basis values into cosines of integer multiples of the angle, so a product of two basis functions becomes an exponential sum whose exponent is a polynomial in $j$ determined by the multi-indices. For $n\ne k$ the difference of multi-indices has a component not divisible by $M$ (guaranteed by $M>2q+1$) and Weil applies; for $n=k$ the exponent vanishes identically and the sum degenerates into a count.

**Step two: the two bounds of Lemma 3.1.**

$$
\Bigl|\sum_{j=0}^{m}\Phi_n(y_j)\Phi_k(y_j)\Bigr|\ \le\ \frac{(d-1)\sqrt{M}+1}{2}\quad(n\ne k),
\qquad
\sum_{j=0}^{m}\Phi_n^2(y_j)\ \ge\ \frac{M}{2^{d+1}}-\frac{(d-1)\sqrt{M}}{2} .
$$

The leading term $M/2^{d+1}$ of the diagonal bound is the mean value $1/2$ of $\cos^2$ in each dimension times the $m+1\approx M/2$ points, multiplied across dimensions; the off-diagonal bound is entirely square-root cancellation.

**Step three: Gershgorin.** After scaling $A$ by $2^{d+1}/M$, the bound on the off-diagonal entries is exactly

$$
\delta=\frac{2^{d+1}}{M}\cdot\frac{(d-1)\sqrt{M}+1}{2}=\frac{2^{d}\bigl((d-1)\sqrt{M}+1\bigr)}{M} ,
$$

and the scaled diagonal entries lie within $\delta$ of $1$. Gershgorin gives $|\lambda_i-1|\le N\delta$, and all that remains is to check that $M\ge4^{d+1}d^2N^2$ forces $N\delta\le\frac12$: the condition reads $\sqrt M\ge2^{d+1}dN$, whence $N\delta\lesssim 2^d(d-1)N/\sqrt M\le\frac{d-1}{2d}<\frac12$. **The entire quantitative content of the paper rests on that one line.**

**Step four: weights for general measures (Section 4.2).** Weighted least squares reads $f_\Lambda=\arg\min_{v\in\mathbb P_\Lambda}\sum_{i=0}^{m}w_i(f(y_i)-v(y_i))^2$. Because $\Theta_M$ equidistributes to the Chebyshev measure (Theorem 4.3), the unweighted discrete norm emulates the Chebyshev norm; to emulate a $\rho$-weighted norm instead, take

$$
w_i=\frac{\rho(y_i)}{\rho_c(y_i)}=\pi^d\rho(y_i)\prod_{q=1}^{d}\bigl(1-(y_i^q)^2\bigr)^{1/2} ,
$$

which for the uniform density $\rho\equiv2^{-d}$ becomes $w_i=(\pi/2)^d\prod_{q=1}^{d}\bigl(1-(y_i^q)^2\bigr)^{1/2}$. Since $w_i$ multiplies a quadratic form, the effect is preconditioning by $\sqrt{w_i}$, that is $\Phi_i(y)\mapsto\prod_{q=1}^{d}(1-(y^q)^2)^{1/4}\Phi_i(y)$ — and **the paper points out that this is exactly the preconditioning known to make the design matrix well conditioned for $\ell_1$ minimisation with Legendre approximations**. When $\rho\propto\rho_c$ the weights are constant. Note that $w_i=\rho/\rho_c$ is the same change-of-measure weight $w(z_s)/v(z_s)$ as in the shared machinery above.

### Theorems

- **Lemma 3.1.** Hypotheses: $M$ prime, $M>2q+1$, $y_j=\cos(x_j)$ as above. Conclusion: the two bounds displayed above.
- **Theorem 3.2 (stability).** Hypotheses: $M\ge4^{d+1}d^2N^2$ and $M$ prime. Conclusion: $\bigl|\!\bigl|\!\bigl|\frac{2^{d+1}}{M}A-I\bigr|\!\bigr|\!\bigr|\le\frac12$ in the spectral norm. **This is the quadratic sample-complexity statement $M\gtrsim N^2$, with an explicit dimension-dependent prefactor $4^{d+1}d^2$.**
- **Corollary 3.3 (uniqueness).** Under the same hypotheses the minimiser of $\sum_{k=0}^{m}(p(y_k)-f(y_k))^2$ over $\mathbb P_\Lambda$ is unique, because $A$ is strictly diagonally dominant hence nonsingular.
- **Theorem 3.4 (optimal convergence in the Chebyshev measure).** With $P^\Lambda f=\arg\min_{p\in\mathbb P_\Lambda}\|f-p\|_{L^2_{\rho_c}}$ and $P^\Lambda_m f$ the discrete least-squares solution on $\Theta_M$, if $M\ge4^{d+1}d^2N^2$ is prime then

  $$
  \|f-P^\Lambda_m f\|_{L^2_{\rho_c}}\ \le\ \Bigl(1+\frac{4}{d^2 N}\Bigr)\,\|f-P^\Lambda f\|_{L^\infty} .
  $$

  The factor $1+4/(d^2N)$ tends to $1$ as $N$ grows, so the discrete projection is asymptotically as good as the best $L^\infty$ approximation. **The estimate is deterministic, with no "with high probability".** That is what the quadratic sample count buys.
- **Definition 3.5 and Corollary 3.6 (other measures).** Say $\rho$ is bounded by the Chebyshev density if there is a constant $C$ independent of $Y$ with $0<\rho(Y)\le C\rho_c(Y)$ for all $Y\in\Gamma$. Then for any $f\in L^2_{\rho_c}$,

  $$
  \|f-P^\Lambda_m f\|_{L^2_\rho}\ \le\ \sqrt{C}\,\Bigl(1+\frac{4}{d^2 N}\Bigr)\|f-P^\Lambda f\|_{L^\infty} ,
  $$

  again for prime $M\ge4^{d+1}d^2N^2$. Remark 3.8 notes this covers the uniform measure and every measure with $0<\rho_{\min}\le\rho\le\rho_{\max}$; Remark 3.7 notes the use for **epistemic** uncertainty: when the density of $Y$ is unknown, Chebyshev-based approximation stays efficient provided the unknown density satisfies the domination condition.
- **Theorem 4.3 (asymptotic equidistribution).** Let $M_K$ be the $K$th prime, $m_K=\lfloor M_K/2\rfloor+1$, and let $\nu_K:=\frac{1}{m_K}\sum_{j=1}^{m_K}\delta(y_{j,K})$ be the empirical measure, $\nu_c$ the normalised Chebyshev measure with density $\rho_c(y)=\pi^{-d}\prod_q(1-y_q^2)^{-1/2}$. Then $\nu_K\to\nu_c$ weakly. The proof runs Weil's formula through **Weyl's equidistribution criterion** (Theorem 4.1, Corollary 4.2). **This result is the licence for the weight formula**: without it the derivation of $w_i=\rho/\rho_c$ has no basis.

### Numerical experiments

Section 5 compares $\Theta_M$ against Monte Carlo grids in both TP and TD spaces, plotting deterministic results as dots and Monte Carlo results as squares.

| Figure    | Setup                                                       | Comparison                                                                       | Observation                                                            |
| --------- | ----------------------------------------------------------- | -------------------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| Fig. 1    | $d=2$ point distributions at $M=997$                        | one Monte Carlo realisation from the Chebyshev measure                           | both cluster near the boundary, as Theorem 4.3 predicts                |
| Figs. 2, 3 | two-dimensional TD and TP spaces, plus a three-dimensional TD case | linear rule $m=c(\#\Lambda)$ versus quadratic rule $m=c(\#\Lambda)^2$, several $c$ | design-matrix condition numbers; deterministic points behave like random ones |

The paper's own qualitative conclusion, stated in the abstract, is that the deterministic points perform similarly to randomly generated points — **the gain is not accuracy but the removal of probabilistic qualifiers**. It is also candid that its choices of $d$, $q$ and test functions are not special and that other choices behave similarly.

One gap between theory and experiment is worth naming. The prefactor the theorem demands is $4^3\cdot4=256$ at $d=2$ and $4^4\cdot9=2304$ at $d=3$, while the experiments run at $m=c(\#\Lambda)$ and $m=c(\#\Lambda)^2$ (the values of $c$ are unverified here). Unless $c$ reaches that magnitude, the experiments live in a regime the theorem does not cover — and the conditioning is good anyway. **That suggests Theorem 3.2 is sufficient rather than sharp**, a point the paper does not pursue.

### Relation to the others

**This paper establishes the template for the family**: an explicit point set or density, a stability theorem quantifying the required sample count, and a near-best-approximation corollary. Its $M\gtrsim N^2$ is the benchmark the later papers attack — papers 22 and 45 show that sampling from the Christoffel-weighted or induced density brings the requirement down to $M\sim N\log N$, and paper 28 replaces random sampling entirely with weighted approximate Fekete points. The weight $w_i=\rho(y_i)/\rho_c(y_i)$ introduced here is the direct ancestor of the Christoffel weights of paper 22 and the preconditioners of papers 24 and 32. Xu's Weil-sum machinery reappears in paper 10 (deterministic interpolation points) and paper 29 (sparse trigonometric recovery), on [[en/computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce|Sparse recovery and data-driven chaos]]. The empirical question of paper 6 is answered here with theory for the first time.

## 11: an unbounded domain pushes the sample count from polynomial to exponential

### The idea

All existing stability theory concerned **bounded** parameter domains: uniform or Chebyshev measures on $[-1,1]^d$, where $M\sim N^2$ suffices and Chebyshev needs fewer. Gaussian and Gamma parameters live on unbounded domains, and the paper shows the situation there is worse **in kind, not merely in constants**.

Why it is a difference in kind can be read straight off the shared criterion $M/\log M\gtrsim\lVert K\rVert_\infty$. On a bounded domain that supremum is taken over a **compact** set: the basis is continuous, the supremum is finite, and it grows only polynomially in $N$ ($N^2$ for Legendre, $N^{1.585}$ for Chebyshev). On an unbounded domain the set is no longer compact, and two things happen at once and reinforce each other: orthogonal polynomials have no uniform bound far out, while the Gaussian measure has almost no mass there. **Nearly every sample lands where the basis is small, and the occasional far-out sample produces an enormous row norm** — row norms are wildly non-uniform, which is exactly the shape of an exploding stability factor. No constant factor repairs that.

The diagnosis has internal corroboration: the paper's Lemma 3.1 gives a uniform decay bound for Hermite **functions** and not for Hermite **polynomials**, and the two differ only by a factor $e^{-y^2/2}$. **What that factor cures is precisely the unbounded growth of the polynomials far from the origin.**

Seen this way the remedy has a one-line summary: **move the weight out of the measure and into the basis**. Once the weight sits in the basis, the basis has a uniform bound and the stability factor is controlled; the sampling density is then free to be chosen for coverage rather than for orthogonality, and is taken to be a mapped uniform law. This is the same move as the Christoffel weighting of papers 22 and 45, performed differently — those keep the polynomial basis and instead change the density and add weights.

A second, independent problem is the poor resolution of Hermite expansions: the paper quotes Gottlieb and Orszag's remark that resolving $M$ wavelengths of $\sin(x)$ needs nearly $M^2$ Hermite polynomials. That is a problem of convergence speed rather than stability, and the paper meets it with a separate device, the scaling factor.

> [!warning] The notation here clashes with paper 9
> In this paper $m$ is the **number of random points**, $K$ is the **number of basis functions**, and $M$ is the effective support radius of the target function ($|f(y)|<\epsilon$ for $|y|>M$). In the section on paper 9, $M$ is the prime and $N$ is the dimension of the space. Keep them apart when reading the theorems.

### Setting

Normalised Hermite polynomials $\{H_k\}$ are orthonormal on $\mathbb R$ against $\rho_G(y)=e^{-y^2}$, that is $\int_{-\infty}^{+\infty}\rho_G(y)H_m(y)H_n(y)\,\mathrm dy=\delta_{mn}$; Laguerre polynomials are orthonormal against $\rho_E(y)=\prod_{i=1}^{d}e^{-y_i}$, the exponential special case of Gamma, with the general Gamma density $\rho_E(y)=\beta^\alpha y^{\alpha-1}e^{-\beta y}/\Gamma(\alpha)$ handled by generalised Laguerre chaos. Stability is measured by $\mathrm{cond}(A)=\sigma_{\max}(A)/\sigma_{\min}(A)$. The normal equations are $\bigl(\langle\Phi_i,\Phi_j\rangle_m\bigr)_{i,j=1,\dots,N}$ with right-hand side $f=D^{\mathsf T}b=\bigl(\langle f,\Phi_j\rangle_m\bigr)_j$, solved by QR on the design matrix or Cholesky on the normal equations.

### Derivation

**Diagnosis (Section 2.3).** Sampling from the Gaussian or Gamma measure with the corresponding polynomial chaos basis, the condition number grows **exponentially** with polynomial order under both the linear rule $m=c(\#\Lambda)$ and the quadratic rule $m=c(\#\Lambda)^2$. The paper writes the requirement as

$$
m=(\#\Lambda)^{\,c\,\#\Lambda} ,
$$

which it calls unacceptable for practical computation. **This is the paper's own characterisation of a numerically observed scaling, not a proven lower bound** — this page does not cite it as a theorem.

**Remedy one: Hermite/Laguerre functions rather than polynomials (eqs. (3.1)–(3.3)).**

$$
\tilde H_m(y)=e^{-\frac{y^2}{2}}H_m(y),\qquad
\tilde L_m(y)=e^{-\frac{y}{2}}L_m(y),\qquad m=0,1,\dots
$$

These are orthonormal against **Lebesgue** measure, $\int_{-\infty}^{+\infty}\tilde H_m\tilde H_n\,\mathrm dy=\delta_{mn}$ — the Gaussian weight has been absorbed into the basis. Multivariate versions are tensorised. They are no longer polynomials, but the paper keeps calling the index $q$ the "polynomial order". The justification that decaying bases are the right ones for UQ (eq. (3.4)): a quantity of interest has the form $\mathrm{QoI}=\int_\Gamma\rho(y)(g\circ f)(y)\,\mathrm dy$, and even when $g\circ f$ does not decay, $\rho\cdot(g\circ f)$ does provided $g\circ f$ grows more slowly than a Gaussian, so one approximates $\tilde f(y)=\rho(y)(g\circ f)(y)$ instead.

**Remedy two: mapped uniform random points (eq. (3.7)).** Points are drawn uniformly on a bounded interval and mapped to the unbounded domain with a mapping parameter $L$:

$$
y(\xi)=\begin{cases}\dfrac{L}{2}\log\dfrac{1+\xi}{1-\xi}, & r=0,\\ \dfrac{L\xi}{\sqrt{1-\xi^2}}, & r=1,\end{cases}
\qquad
\xi(y)=\begin{cases}\tanh\bigl(\tfrac{y}{L}\bigr), & r=0,\\ \dfrac{y/L}{\sqrt{y^2/L^2+1}}, & r=1.\end{cases}
$$

The $r=0$ **logarithmic** mapping makes the transformed points decay exponentially and is used for the Gaussian measure; the $r=1$ **algebraic** mapping is used for the Gamma measure. $L$ is the knob that controls conditioning.

**Remedy three: a scaling factor for convergence (eqs. (3.28), (3.32), (3.33)).** Let $f$ decay exponentially, with $|f(y)|<\epsilon$ for $|y|>M$. Write the expansion as

$$
f(y)=\sum_{n=0}^{K-1}c_n\tilde H_n(\alpha y)
\ \Longleftrightarrow\
f\Bigl(\frac{y}{\alpha}\Bigr)=\sum_{n=0}^{K-1}c_n\tilde H_n(y),
\qquad \alpha>0,
$$

so that the scaled points $y_i/\alpha$ sit inside the effective support of $f$. Coefficients are computed as $f_k=\langle f,H_k\rangle_m=\sum_{i=1}^{m}f(y_i/\alpha)H_k(y_i)$, so the requirement $\max_j|y_j|/\alpha\le M$ gives the naive rule $\alpha=\max_{1\le j\le m}\{|y_j|\}/M$. But the points are random, and a handful of extreme "bad points" over-scale, so the paper's **quasi-optimal scaling** discards them:

$$
\tilde\alpha=\max_{1\le j\le\tilde m}\{|y_j|\}/M,
\qquad \tilde m=\lfloor\mu m\rfloor,
$$

with $\mu$ close to 1; in practice $\mu\approx98\%$, so the most extreme 2% of samples are dropped. **This trimmed maximum is the paper's most useful practical contribution**, and in the experiments it is what decides success or failure.

**The route of the stability proof.** Matrix Chernoff: write the scaled design matrix as $\hat A=X_1+\cdots+X_m$ with i.i.d. copies of

$$
X=\frac{L}{m}\bigl(\tilde H_i(y)\tilde H_j(y)\bigr)_{i,j=0,\dots,K-1} ,
$$

$y$ a transformed uniform random variable, then apply a Chernoff bound for independent positive self-adjoint random matrices, which requires an almost-sure bound $\lambda_{\max}(X_i)=|\!|\!|X_i|\!|\!|\le R$. Since $X$ has rank one,

$$
\lambda_{\max}(X)=\frac{L}{m}\sum_{i=0}^{K-1}\tilde H_i^2(y) ,
$$

**which is precisely the squared row norm of the shared machinery**. The hypotheses $L>3\tau$ and $L>5\sqrt K$ are consumed exactly here: the first activates the decay bound of Lemma 3.1, the second makes the mapped interval wide enough to cover the range over which the first $K$ basis functions are active. The analysis is written for the one-dimensional Hermite case only; the Laguerre case is said to follow in a straightforward manner — **it is not written out in what was verified here**.

The two-parameter structure the paper stresses in its conclusion is worth remembering: $L$ controls **stability** and $\alpha$ controls the **rate of convergence**, and they are tuned for different purposes.

### Theorems

- **Lemma 3.1.** For any integer $K$ there is a constant $\tau>0$ such that $|\tilde H_k(y)|\le|y|^{-3/2}$ for all $0\le k\le K-1$ whenever $|y|>\tau$. (True because the factor $e^{-y^2/2}$ forces $|\tilde H_k(y)|\cdot|y|^t\to0$ for every $t>0$.) **The essential point is that the bound is uniform in $k$.**
- **Theorem 3.2 (stability).** With the Hermite functions (3.1) and the transformed uniform random points (3.7), the scaled design matrix $\hat A=LA$ satisfies, for every $r>0$,

  $$
  \Pr\Bigl\{\,\bigl|\!\bigl|\!\bigl|\hat A-I\bigr|\!\bigr|\!\bigr|\ \ge\ \tfrac58\,\Bigr\}\ \le\ 2m^{-r} ,
  $$

  **provided**

  $$
  K\ \le\ \kappa\,\frac{m}{\log m},
  \qquad
  \kappa:=\frac{4c_{1/2}}{3(1+r)},
  \quad c_{1/2}=\frac12+\frac12\log\frac12>0 ,
  $$

  and the mapping parameter satisfies $L>\max\{3\tau,\ 5\sqrt K\}$, where $m$ is the number of random points and $K$ the number of basis functions. **This is the log-linear sample complexity $K\lesssim m/\log m$, equivalently $m\gtrsim K\log K$ up to constants** — a qualitative improvement over the exponential requirement of the polynomial chaos formulation, and matching in form the best bounded-domain results.

  Converting the constant makes the size of the demand visible: $c_{1/2}=\frac12(1-\ln2)\approx0.153$, so $\kappa\approx0.205/(1+r)$, and at $r=1$ the condition reads $K\lesssim0.10\,m/\log m$. **The theorem is not cheap**: it asks for about an order of magnitude more samples than the bare rate suggests.

The $\sqrt K$ inside the condition $L>5\sqrt K$ is the same phenomenon as the contraction factor $k^{-1/r}$ of Theorem 14 in paper 14 (with $r=2$ in the Gaussian case). That is no coincidence: both say that the region in which $K$ basis functions are active expands like $\sqrt K$.

### Numerical experiments

Condition numbers are measured as $\mathrm{cond}(A)=\sigma_{\max}(A)/\sigma_{\min}(A)$ and averaged over **100 independent realisations**, since the matrix is random.

**Stability.**

| Figure | Basis and sampling                     | Dimension and sampling rule                                             | Observation                                                      |
| ------ | -------------------------------------- | ------------------------------------------------------------------------ | ----------------------------------------------------------------- |
| Fig. 1 | Hermite polynomial chaos, Gaussian samples | one dimension, $m=c(\#\Lambda)$ and $m=c(\#\Lambda)^2$                  | condition number grows exponentially with polynomial order — the negative result |
| Fig. 4 | Laguerre functions, mapped uniform points | 1D with $m=30(\#\Lambda)$; 2D with $m=6(\#\Lambda)^2$, in TD and TP spaces | mapping parameters $L=8$ and $L=64$ compared                     |

**Convergence.**

| Case          | Target function                                | Key parameters                                                                                    | Result                                                                                  |
| ------------- | ---------------------------------------------- | -------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- |
| Fig. 6        | $f(y)=2^{-6y^2}$ (fast decay)                  | $L=8$ fixed for stability                                                                          | $\alpha=1$ (no scaling) converges very slowly but stably; scaling accelerates markedly, best in that test around $\alpha=2.8$ |
| second 1D test | $\tilde f(y)=2^{-0.2y^2}$                      | $M=16$                                                                                             | —                                                                                        |
| 2D test       | $\tilde f(y)=e^{-4(y_1^2+y_2^2)}\sin(y_1+y_2)$ | $M=2.5$, TD space; $m=10(\#\Lambda)$ versus $m=2(\#\Lambda)^2$; no scaling / $\mu=0.985$ / $\mu=0.980$ | the trimmed scaling (3.33) converges fast, while the plain maximum (3.32) behaves almost like no scaling |

The consistent finding is that **what decides the outcome is discarding the bad points, not the scaling as such**. The parametric UQ applications use a random ODE model and an elliptic problem with lognormal random input. One further optimal scaling of about $\alpha=0.82$ is reported (which test it belongs to is unverified here).

What the experiments establish, and where they fall short, are worth separating. They do establish two things: the failure of polynomial chaos on unbounded domains is real and measurable, and after changing basis and adding the mapping the condition number is controlled with $L$ as an effective knob. What they do not establish is the convergence half — **no theorem about the choice of $\alpha$ was found in the source; the abstract's word is "discussed"** — and the two reported optima ($2.8$ and $0.82$) differ by more than a factor of three, so they are case-specific observations rather than general prescriptions. No error magnitudes were verified here, so no error table is given.

### Relation to the others

This is the unbounded-domain counterpart of the bounded-domain stability theory, and the first paper in the topic to reach a **log-linear** rather than quadratic sample requirement — but it buys that by changing both the **basis** (polynomials to functions) and the **sampling density** (Gaussian to mapped uniform), rather than by weighting. **That route and the later Christoffel weighting are two distinct choices**: papers 22 and 45 keep the polynomial basis and instead change the density and add weights, reaching log-linear complexity at a different price — the former gives up the orthogonal structure of the polynomial basis, the latter needs the Christoffel function computed. The unbounded-domain concern reappears in papers 24 and 32, whose preconditioner framework is explicitly advertised as covering bounded or unbounded domains. The bounded-domain quadratic benchmark $m\sim(\#\Lambda)^2$ that this paper cites is the one paper 9 established for deterministic points.

## 13: randomised quadrature, and a paper that could not be verified

### The idea

The two ends of the design space are clear. At one end, i.i.d. sampling from the orthogonality measure: simple to describe, expensive in samples. At the other, the full tensor-product Gauss grid: perfect discrete orthogonality, but cardinality $\prod_i n_i$, exponential in dimension. The idea in between is natural — **randomly subsample the Gauss grid and reweight**, keeping the excellent discrete-orthogonality structure of Gauss quadrature while paying for only $M$ expensive evaluations.

What makes it a principled reweighting rather than an arbitrary one is the identity from the shared machinery: **Gauss weights are Christoffel function values**, so subsampling a Gauss grid is implicitly Christoffel-weighted sampling.

### Setting

The construction below is taken from the sister paper 21, which uses an identical design and analyses it for $\ell_1$ recovery; **equation numbers and any variants specific to this paper are unverified here**.

Tensor-product Gauss grid: with $\Theta^i_{n}=\{z^i_1,\dots,z^i_n\}\subset\Gamma_i$ the $n$-point Gauss set in dimension $i$, the tensor set is $\Theta_{\mathbf n}=\Theta^1_{n_1}\otimes\cdots\otimes\Theta^d_{n_d}$ with $|\Theta_{\mathbf n}|=\prod_{i=1}^{d}n_i$. The Gauss weights are

$$
w_{\mathbf k}=\lambda_{\mathbf n}(z_{\mathbf k})
=\prod_{i=1}^{d}\lambda^i_{n_i}(z^i_{k_i})
=\prod_{i=1}^{d}\frac{1}{\sum_{k=0}^{n_i-1}\bigl[\phi^i_k(z^i_{k_i})\bigr]^2} .
$$

The uniform empirical measure on the grid is $\nu_{\mathbf n}=\bigotimes_{i=1}^{d}\nu^i_{n_i}=\frac{1}{\prod_i n_i}\sum_{\mathbf k\le\mathbf n}\delta_{z_{\mathbf k}}$, and drawing i.i.d. from $\nu_{\mathbf n}$ is exactly uniform sampling from the tensor-product Gauss grid. The weighted design matrix is $D=\sqrt W\Psi$ with $(\Psi)_{m,n}=\varphi_n(x_m)$ and $(W)_{m,m}=w_m>0$; the weighted problem replaces $\Psi c=f$ with $Dc=\sqrt W f$.

Indexing metadata (verifiable): keywords are least squares method, uncertainty quantification, orthogonal polynomials, generalized polynomial chaos; MSC 65C20, 65D30, 41A10, 65D15, 60G99. Its **reference list** — not its results — reveals the technical apparatus: two papers of Nevai on generalised Jacobi weights, Christoffel functions and Jacobi polynomials; a paper on the asymptotic expansion of Gaussian quadrature weights; Tang–Iaccarino on subsampled Gauss quadrature nodes for polynomial chaos; Cohen–Davenport–Leviatan; Migliorati et al. on discrete $L^2$ projection with random evaluations; and **four monographs on experimental design**. The Christoffel references confirm that the asymptotics of Gauss weights as Christoffel function values is the analytical backbone, and the design-of-experiments references point to an explicit link with optimal experimental design that none of the neighbouring papers has.

### What could not be verified

> [!warning] No theorem, constant or numerical result of this paper is reported here
> The publisher blocks the full text, researchr explicitly records that the abstract is missing, and Crossref, OpenAlex, MaRDI/zbMATH and Semantic Scholar hold no abstract either (Semantic Scholar's record notes that the abstract field has been elided by the publisher). ScienceDirect returns HTTP 403 for both the article and the accepted-manuscript URL. **No theorem statement, no constant, and no explicit sample-complexity relation could be verified from the paper itself, and neither could any numerical experiment.**
>
> The sole basis for the claim that the sample count grows **linearly** in the polynomial dimension is a third-party statement: Seshadri, Narayan and Sarkar, in _Quadrature Strategies for Constructing Polynomial Approximations_, write that Zhou et al. randomly subsample the rows and demonstrate stability of the least-squares problem with $m$ scaling linearly with $n$. One of its authors is also an author of this paper, so it is a credible corroboration, **but it is not the paper's own wording, and the precise hypotheses, probability statement and constants are all unverified**. No downstream citation should attribute a specific constant or probability bound to this paper.

What can be verified from the sister paper 21, the closest analogue, is this: for the $\ell_1$ / compressed-sensing version of the same sampling scheme, the sufficient sample count is $M\ge L(\mathbf n)C_1s\log^3(s)\log(N)$, and for a uniform random variable on a hypercube the earlier work of Tang and Iaccarino gives $L\le C^d$ with $C$ essentially $3$, that is $M\gtrsim3^d s$. **The $3^d$ deserves attention**: inside the restricted design space of tensor Gauss grids the dimension dependence has not disappeared, it has moved into the constant. "Linear" should therefore be read as linear in $N$ with a prefactor growing in $d$, not as unconditionally better than the $M\gtrsim N\log N$ of induced sampling.

### Relation to the others

This is the direct predecessor of paper 21 — the same group of authors, the same device of randomly subsampling Gauss quadrature, but with least squares instead of $\ell_1$ minimisation. Its Gauss weights are Christoffel function values, which links it forward to the explicit Christoffel-weighted framework of paper 22 and the induced-sampling theory of paper 45. It appears in the reference lists of papers 21, 28 and the review, paper 14.

## 14: putting all three reconstruction modes in one ledger

### The idea

By 2015 three families of collocation-based reconstruction were in wide use — least-squares regression, compressive sampling, and interpolation — each with its own theory and its own preferred point sets, and with no unified account of what makes a **geometrically unstructured** multivariate mesh good. This is the survey-and-synthesis paper of the programme: **it collects the stability and accuracy results of all three modes into one notation** and then uses them as design guidance for generating stochastic collocation grids in several dimensions. Being a review with new numerical comparisons, most of its theorems are attributed to prior work; the contribution is the organisation plus the head-to-head experiments.

The deeper answer it gives is that the requirement has the same shape in all three modes: **the empirical measure of the point set must approach a measure determined by the domain itself (the pluripotential equilibrium measure), and that is necessary but not sufficient.** This is why Chebyshev-like boundary clustering keeps reappearing in otherwise unrelated methods.

### Setting

A model $\mathcal L(u;t,x,\omega)=0$ whose randomness is parameterised by a $d$-dimensional random vector $Z(\omega)$, typically after a Karhunen–Loève truncation, with a gPC expansion of the response.

**Weil points** (eqs. (16)–(17)): the deterministic Weil-sum grid of papers 9 and 10, with $\bigl|\sum_j e^{2\pi if(j)/M}\bigr|\le(d-1)\sqrt M$ as the analytical engine.

**Structured random points** (Section 3.3): take a high-cardinality structured candidate set — a tensor-product Chebyshev grid, or a hyperbolic cross space $\mathcal H^d_k$ viewed as a subspace of the tensor space, $u(z)=\sum_{\alpha\in\Lambda^P_{d,k}}\hat c_\alpha\varphi_\alpha(z)$ with $\hat c_\alpha=0$ for $\alpha\in\Lambda^P_{d,k}\setminus\Lambda^H_{d,k}$ — and randomly extract a subset, yielding an essentially unstructured grid. **This is the design used by papers 13 and 21.**

**Least orthogonal interpolation** (Section 5.1): introduce the "least-$\rho$" operation $p_{\downarrow,\rho}=P_{\hat k}p$ with $\hat k=\min\{k\in\mathbb N:P_k\ne0\}$, the first nonvanishing "Taylor" contribution, which depends on $\rho$. For nodes $Z=\{z_1,\dots,z_M\}$ define

$$
\Pi_Z=\mathrm{span}\bigl\{g_{\downarrow,\rho}\ :\ g\in\mathrm{span}\{\delta_{z_1},\dots,\delta_{z_M}\}\bigr\} ,
$$

the least orthogonal polynomial space for interpolation. It is computed in practice by a combination of LU and QR on the rectangular design matrix, yielding the factorisation $PA=LUH$ (eq. (39)) with $L,U$ standard $M\times M$ triangular factors, $P$ a row permutation and $H$ rectangular; the operation count is asymptotically comparable to standard interpolatory matrix factorisations. When $\rho$ is the standard Gaussian density, $\Pi_Z$ coincides with the classical least interpolant space of de Boor and Ron.

### Derivation: converting every requirement into the same units

There is no single derivation running through the paper. Its move is to convert mutually incomparable results into one unit — $M$ against $N$, or against $s$ in the sparse case — and line them up:

| Reconstruction mode and sampling                | Sample requirement                     | Source                            |
| ----------------------------------------------- | -------------------------------------- | --------------------------------- |
| least squares, i.i.d. from the uniform (orthogonality) measure | $M/\log M\ge C_rN^2$                   | Theorem 1 (CDL)                   |
| least squares, Chebyshev Monte Carlo            | $M\sim N^{\log3/\log2}$                | the paper's contrast              |
| least squares, Hermite functions and mapped points | $M/\log M\gtrsim rN$, $L\gtrsim\sqrt N$ | Theorem 2 (from paper 11)         |
| least squares, deterministic Weil points        | $M\ge C(d)N^2$                         | Theorems 5, 6 (from paper 9)      |
| $\ell_1$, Chebyshev samples with preconditioner | $M>C\delta^{-2}L^2s\log^3(s)\log(N)$   | Theorems 8, 9                     |
| what practitioners actually do                  | $M\simeq cN$ with $c$ between 2 and 3  | the paper's account of practice   |

The last row is the point of the table: **practice runs in the linear regime, and by the paper's own statement the theory for that regime "is not yet definitively available".** The organisation of the whole paper is built around that gap.

### Theorems

- **Theorem 1 (Cohen–Davenport–Leviatan).** Hypotheses: $\rho$ uniform on $[-1,1]$ and $\Lambda=\Lambda^T_{1,N-1}$. If for some $r>0$ we have $\frac{M}{\log M}\ge C_rN^2$ with a universal constant $C$, then $\Pr\bigl[|\!|\!|\hat A-I|\!|\!|\ge\frac12\bigr]\le2M^{-r}$. **This is the quadratic requirement $M\gtrsim N^2\log N$ for i.i.d. sampling from the orthogonality measure.** The paper notes the extension to multidimensional spaces with arbitrary **lower** index sets under the same scaling.
- **Contrast.** For **Monte Carlo sampling from the Chebyshev measure** the requirement drops to $M\sim N^{\log3/\log2}\approx N^{1.585}$ — strictly better than $N^2$ but still superlinear.
- **Theorem 2 (from paper 11).** For any $r>0$, if $\frac{M}{\log M}\gtrsim rN$ and $L\gtrsim\sqrt N$, then the least-squares design matrix built from Hermite **functions** with mapped uniform points satisfies $\Pr\bigl[|\!|\!|\hat A-I|\!|\!|\ge\frac58\bigr]\le2M^{-r}$. The paper restates it as a **weighted** least-squares result, since Hermite functions are weighted Hermite polynomials — the restatement is itself a demonstration of the unified notation.
- **Theorem 4 (from paper 9).** The Weil points $W_{M_K}$ generated by the $K$th prime distribute asymptotically according to the Chebyshev measure, $\nu_K\to\nu_c$ in the weak-$*$ topology.
- **Theorems 5 and 6 (from paper 9).** If $M\ge C(d)N^2$ then $\bigl|\!\bigl|\!\bigl|\frac{2^{d+1}}{M}\hat A-I\bigr|\!\bigr|\!\bigr|\le\frac12$, hence a unique least-squares solution, and $\|f-P^N_Mf\|_{L^2_\rho}\le C\|f-p^*_N\|_{L^\infty}$ with $p^*_N$ the $L^\infty$-best polynomial. The paper is explicit that this quadratic requirement is **stronger** than the $M\sim N^{\log3/\log2}$ needed by Chebyshev Monte Carlo, and that the compensating advantage is determinism.
- **Theorem 8 (Rauhut, RIP for bounded orthonormal systems).** Let $\sup_n\|\varphi_n\|_\infty\le L$ with $L\ge1$, let $A$ be the interpolation matrix and $W$ diagonal with $w_{m,m}=(\pi/2)^{1/2}(1-z_m^2)^{1/4}$, the $z_m$ i.i.d. from the one-dimensional Chebyshev measure. If $M>C\delta^{-2}L^2s\log^3(s)\log(N)$ then with probability at least $1-N^{-\gamma\log^3(s)}$ the restricted isometry constant of $\frac{1}{\sqrt M}WA$ satisfies $\delta_s\le\delta$.
- **Theorem 9 (Rauhut–Ward).** Same sampling, with $A$ the Legendre design matrix and $W$ diagonal with $w_{m,m}=(\pi/2)^{-1/2}(1-z_m^2)^{1/4}$. If $M>Cs\log^3(s)\log(N)$, the solution of $\min\|c\|_1$ subject to $WAc=WA\tilde c$ satisfies $\Pr\bigl[\|c^{\#}-\tilde c\|_2\le C\sigma_{s,1}(\tilde c)/\sqrt s\bigr]\ge1-N^{-\gamma\log^3(s)}$. The paper explains what the $(1-z^2)^{1/4}$ weight is for: **it is exactly what makes the weighted Legendre polynomials a uniformly bounded system**. Theorem 10 (Yan–Guo–Xiu) extends this to high dimension, and Theorem 11 (from paper 10) gives the deterministic Weil-point counterpart via the incoherence parameter.
- **Theorem 12 (least orthogonal interpolation, from Narayan–Xiu).** $\Pi_Z$ has dimension exactly $M$; for any continuous $u$ there is a unique $p\in\Pi_Z$ with $u_m=p(z_m)$; and Lagrange functions $\ell_m\in\Pi_Z$ exist with $p(z)=\sum_{m=1}^{M}u_m\ell_m(z)$ and $\ell_m(z_n)=\delta_{m,n}$.
- **Theorem 13 (Bloom / Bos et al. / Berman–Boucksom, the unweighted trichotomy).** For an array $Z_{N_k}$ consider: (1) subexponential Lebesgue constant growth, $\lim_k(\Lambda(Z_{N_k}))^{1/k}=1$; (2) asymptotically Fekete, $\lim_k|\det A(Z_{N_k})|^{1/s^d_k}=\delta(D)$; (3) distribution according to the pluripotential equilibrium measure, $\lim_k\frac{1}{N_k}\sum_n\delta_{z_{n,N_k}}=\mu_D$. Then $1\Rightarrow2$ and $2\Rightarrow3$, and **the reverse implications are false**. The design lesson the paper draws is that a stable interpolation operator **requires** asymptotic sampling according to $\mu_D$ — necessary but not sufficient. On a tensor-product interval $\mu_D$ is the tensor-product arcsine measure, which is why Chebyshev-like clustering keeps reappearing. The paper also notes that the $d$-dimensional Weil points on a hypercube do distribute according to $\mu_D$ (Theorem 4), but whether they are asymptotically Fekete is **unknown**, and the restriction on their cardinality limits their usefulness for interpolation.
- **Theorem 14 (weighted version, from Berman–Boucksom and Narayan–Xiu).** On $D=\mathbb R^d$ with weights $\rho(z)=\exp(\|z\|^r)$, $r\ge1$, and log-weight $Q(z)\triangleq-\log\rho=\|z\|^r$, let $\mu_{D,Q}$ be the weighted equilibrium measure, whose support is compact even though $D$ is not. Consider: (1) subexponential weighted Lebesgue growth; (2) contracted asymptotically weighted Fekete, $\lim_k\bigl|\det A(k^{-1/r}Z_{N_k})\prod_n\rho^k(z_n)\bigr|^{1/s^d_k}=\delta_\rho(D)$; (2a) uniform contracted boundedness, i.e. a compact $S\supset\mathrm{supp}\,\mu_{\Omega,Q}$ with $k^{-1/r}Z_{N_k}\subset S$ for all $k$; (3) $\lim_k\frac{1}{N_k}\sum_n\delta_{k^{-1/r}z_{n,N_k}}=\mu_{D,Q}$. Then $1\Rightarrow2$ and $(2+2a)\Rightarrow3$, with the reverse implications false. **The essential ingredient is the contraction factor $k^{-1/r}$.** The paper draws out the counterintuitive part: one should **not** sample directly from $\mu_{D,Q}$, since its compact support makes it of limited use for polynomial approximation on an unbounded domain; what one wants are grids whose $k^{-1/r}$-contractions distribute according to $\mu_{D,Q}$. It also records an **open question**: whether multidimensional weighted approximate or discrete Fekete and Leja arrays distribute according to the weighted pluripotential equilibrium measure is unknown, though the affirmative one-dimensional answer gives hope.

> [!warning] The weight and log-weight of Theorem 14 are not sign-consistent
> The form recorded here is $\rho(z)=\exp(\|z\|^r)$ together with $Q(z)\triangleq-\log\rho=\|z\|^r$, and the two cannot both hold, since the first gives $-\log\rho=-\|z\|^r$. One of them carries a sign misprint; the standard convention in this literature is $\rho=\exp(-\|z\|^r)$ with $Q=\|z\|^r$. Whether the misprint is in the source or in transcription could not be determined here, so it is reproduced and flagged.

### Numerical experiments

| Section | Reconstruction mode | Point sets compared                                        | Sampling rules                   |
| ------- | ------------------- | ----------------------------------------------------------- | -------------------------------- |
| 3.4     | least squares       | subsampled Gauss grid, i.i.d. random points, Weil points    | $M=2.5N$ and $M=1.5N\log N$      |
| 4.3     | compressive sampling | Monte Carlo versus deterministic sampling strategies        | —                                |
| 5       | interpolation       | the least orthogonal interpolant                            | —                                |

Section 3.4 plots design-matrix condition numbers against polynomial order $k$; Figure 2 shows two-dimensional Weil grids from the prime seeds $M=359$ (179 points) and $M=751$ (375 points).

The choice of sampling rules is what deserves attention: $M=2.5N$ and $M=1.5N\log N$ are **linear** and **log-linear** respectively, so **the experiments are deliberately run below the quadratic regime that the paper's own Theorems 1, 5 and 6 require**. That is not an oversight; it is the empirical form of the last row of the ledger above. Practice works in the linear regime, and the theory covering it did not yet exist. No numbers are supplied here for these figures.

### Relation to the others

This is the hub of the whole list. It restates and unifies paper 9 (Theorems 4, 5, 6), paper 10 (Theorem 11), paper 11 (Theorem 2) and the external Rauhut, Rauhut–Ward and Yan–Guo–Xiu results, and it introduces the potential-theoretic vocabulary — Lebesgue constants, Fekete points, equilibrium measures, contraction factors — that the group then acts on. Paper 28 constructs **weighted approximate Fekete points** by pivoted QR; papers 22 and 45 pursue the Christoffel-weighted and induced-sampling routes towards the linear-in-$N$ regime this paper identifies as desirable but unavailable. Its discussion of structured random subsampling is the design used in papers 13 and 21.

## Sample complexity across the five papers

| No. | Point set or density                                | Sample requirement                                                              | Type of result                     |
| --- | --------------------------------------------------- | -------------------------------------------------------------------------------- | ---------------------------------- |
| 6   | SG, MC and QMC candidate designs                    | no relation between $M$ and $N$ given                                            | empirical                          |
| 9   | Weil-sum deterministic point set                    | $M\ge4^{d+1}d^2N^2$ with $M$ prime                                               | deterministic, no probabilistic qualifier |
| 11  | mapped uniform density with Hermite/Laguerre functions | $K\le\kappa m/\log m$, i.e. $m\gtrsim K\log K$; $L>\max\{3\tau,5\sqrt K\}$      | probabilistic, $\Pr\le2m^{-r}$     |
| 13  | random subset of a tensor Gauss grid                | linear (third-party statement only; unverified from the paper)                   | not verifiable                     |
| 14  | all of the above plus $\ell_1$ and interpolation    | $N^2\log N$ / $N^{1.585}$ / $N\log N$ / $s\log^3 s\log N$                        | survey and unification             |

One judgement runs through all five: **the right form of the question "which points" is "how many points, drawn from which density".** Paper 6 asks the former; from paper 9 onward it is replaced by the latter, and the latter admits theorems. Taken together the five go one step further, because the choice of density itself reduces to a single question — **how to flatten the row norms of the design matrix** — and the quantity that flattens them is the Christoffel function. The Gauss weights of paper 13 and the decaying basis of paper 11 are two faces of that same requirement in two different settings.

> [!note] Coverage status
> Papers 9, 11 and 14 have been checked against their full texts: the settings, derivation chains, theorem hypotheses and constants, and experimental configurations given on this page are transcribed, and for figures whose numbers the papers do not supply only the configuration is reported. Paper 6 reaches abstract and metadata level only, so no theorems or experimental numbers are given for it. Paper 13 could not be verified: the publisher blocks the full text and researchr explicitly records that the abstract is missing, so its construction here is inferred from its sister paper 21, which uses an identical design, and the claim of a linearly growing sample count rests on a third-party statement (Seshadri–Narayan–Sarkar) rather than on the paper itself. **Its constants and probability bounds are not reported.**
>
> Two further cautions are flagged in place: this group of papers uses conflicting conventions for the Christoffel function (this site takes $K=\sum_\alpha\varphi_\alpha^2$, the Christoffel function to be $N/K$ and every weight to be of the form $1/K$, and each paper's own convention must be checked before quoting it); and the asymptotic induced measure in the Gaussian case is only a **conjecture** in every paper that states it, never a theorem.

## Sources for this page

- Z. Gao and T. Zhou, [_On the choice of design points for least square polynomial approximations with application to uncertainty quantification_](https://doi.org/10.4208/cicp.130813.060214a), Commun. Comput. Phys. 16 (2014), pp. 365-381.
- T. Zhou, A. Narayan, and Z. Xu, [_Multivariate discrete least-squares approximations with a new type of collocation grid_](https://doi.org/10.1137/130950434), SIAM J. Sci. Comput. 36(5) (2014), pp. A2401-A2422.
- T. Tang and T. Zhou, [_On discrete least-squares projection in unbounded domain with random evaluations and its application to parametric uncertainty quantification_](https://doi.org/10.1137/140961894), SIAM J. Sci. Comput. 36(5) (2014), pp. A2272-A2295.
- T. Zhou, A. Narayan, and D. Xiu, [_Weighted discrete least-squares polynomial approximation using randomized quadratures_](https://doi.org/10.1016/j.jcp.2015.06.042), J. Comput. Phys. 298 (2015), pp. 787-800.
- A. Narayan and T. Zhou, [_Stochastic collocation on unstructured multivariate meshes_](https://doi.org/10.4208/cicp.020215.070515a), Commun. Comput. Phys. 18 (2015), pp. 1-36.
