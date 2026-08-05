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

The six papers can be read by asking what supplies the rows, what supplies the columns, and how the rows get normalised:

- Paper **10** replaces the rows by a **deterministic** point set, paying the $m\gtrsim s^2$ that the coherence route demands;
- Paper **21** takes the rows to be a random subset of a tensor Gauss grid, and observes that the Gauss weights already are Christoffel weights;
- Papers **29** and **32** extract **$d$ extra rows per sample** (the gradient), paying with row norms that are no longer uniform and so must be preconditioned;
- Papers **36** and **44** change the **columns**: the orthogonal basis is no longer handed down by a known density but has to be estimated from the empirical measure.

> [!note] Notation used on this page
> **Christoffel convention.** Throughout, this site writes the reproducing-kernel diagonal as $K(z)=\sum_{\alpha\in\Lambda}\varphi_\alpha^2(z)$, the Christoffel function as its reciprocal $\lambda_\Lambda(z)=1/K(z)$, and the normalised Christoffel function as $N/K(z)$ with $N=|\Lambda|$. "Christoffel weight" on this page therefore always means $1/K$; the $N/K$ that appears in least-squares objectives differs from $1/K$ by a constant independent of $z$ and does not change the minimiser. The papers themselves disagree: paper 36 writes $K$ for $N/\sum_j\Phi_j^2$ (the normalised Christoffel function itself), while paper 44 writes $\kappa$ for $\frac1N\sum_j\Phi_j^2$ (the reciprocal). Wherever an original formula is quoted it is converted to the convention above and the conversion is stated.
>
> **Samples versus basis functions.** Papers 10, 21, 36 and 44 put the sample count first: $m$ or $M$ samples and $N$ basis functions. Papers 29 and 32 do the opposite: $N$ samples and $M$ basis functions. Each section below follows the convention of the paper it discusses and says so at the start, but be careful when comparing across sections. Note also that $\Gamma$ denotes the frequency index set $\Gamma\subset\mathbb Z^d$ in paper 29 and the parameter domain in papers 21 and 32.

### Two routes to a recovery guarantee

Every theoretical guarantee on this page travels one of two roads. Identify which one first.

**The coherence route.** The mutual coherence of a design matrix $A$ is

$$
\mu(A):=\max_{k\ne j}\frac{|\langle A_k,A_j\rangle|}{\|A_k\|_2\|A_j\|_2},
$$

with $A_k$ the $k$th column. If $c_0$ is $s$-sparse and $\mu<\frac{1}{2s-1}$, then $\ell_1$ minimisation recovers $c_0$ exactly — the Donoho–Huo result for unions of two orthobases, extended by Fuchs and by Gribonval–Nielsen, and also sufficient for stable recovery under noise. Equivalently one can pass through the restricted isometry constant, $\delta_s\le(s-1)\mu(A)$, and invoke a recovery theorem requiring $\delta_s<\frac13$. The arithmetic consequence is unavoidable: $\mu$ is at best of order $1/\sqrt m$, so $\mu<\frac{1}{2s-1}$ forces $m\gtrsim s^2$ — **the sample count is quadratic in the sparsity**. Papers 10 and 32 take this route, as does Theorem 2.1 of paper 29.

**The restricted-isometry route for bounded orthonormal systems.** The restricted isometry constant $\delta_s<1$ is the smallest number for which

$$
(1-\delta_s)\|c\|_2^2\ \le\ \|Dc\|_2^2\ \le\ (1+\delta_s)\|c\|_2^2
$$

holds for all $\|c\|_0\le s$. If $\delta_s\le0.307$ then $c^\sharp=\arg\min\|c\|_1$ subject to $Dc=D\tilde c$ obeys $\|c^\sharp-\tilde c\|_2\le C\sigma_{s,1}(\tilde c)/\sqrt s$, where $\sigma_{s,p}(c)=\inf_{\|y\|_0\le s}\|y-c\|_p$ is the best $s$-term approximation error, with exact recovery when $\tilde c$ is $s$-sparse. A family $\{\psi_k\}$ orthonormal with respect to a density $\nu$ is a **bounded orthonormal system** when

$$
\max_{1\le k\le N}\|\psi_k\|_\infty^2=\max_{1\le k\le N}\ \sup_{x\in\operatorname{supp}\nu}|\psi_k(x)|^2\ \le\ L(N)<\infty,
$$

and uniformly bounded when the bound does not depend on $N$. The theorem of Rauhut and Rauhut–Ward then says: for $x_1,\dots,x_M$ drawn i.i.d. from $\nu$ and $d_{ij}=\psi_j(x_i)$,

$$
M\ \ge\ C\,\delta^{-2}\,L\,s\,\log^3(s)\,\log(N)
$$

implies that the restricted isometry constant of $\frac{1}{\sqrt M}D$ satisfies $\delta_s\le\delta$ with probability at least $1-N^{-\gamma\log^3(s)}$, with $C,\gamma$ universal. **This route brings the sample count down from $s^2$ to $s$ times logarithms, at the price of a probabilistic conclusion and of having to control $L$ — and $L$ is precisely "the uniformity of the row norms".** Paper 21 and Theorem 3.1 of paper 29 take this route.

Legendre polynomials are **not** uniformly bounded, but multiplying them by $(1-x^2)^{1/4}$ makes them so, and the biased measure that keeps the weighted system orthogonal is the Chebyshev measure. That is the Rauhut–Ward strategy, and it is the **reciprocal pairing** that recurs throughout this page: sample from a density $\propto w/(N\lambda_\Lambda)$ and precondition with $N\lambda_\Lambda$, so that the composite system is uniformly bounded. Paper 21 realises it with Gauss weights, paper 32 with $\sqrt{\rho^{(\alpha,\beta)}/\rho_c}$, and paper 44 by sampling from $\mu\propto\kappa$ and weighting by $1/\kappa$.

## 10: deterministic interpolation points and sparse interpolation

(Notation in this section: $m$ interpolation points, $N=\#\Lambda$ basis functions, sparsity $s$.)

### The idea

The paper answers two unrelated questions.

The first is pure existence: **how few points can possibly suffice for $s$-sparse interpolation?** The intuition is clean. If $f,g\in U_s$ agree at the points, then $f-g$ is at most $2s$-sparse and vanishes there, so "unisolvent" means "no nonzero $2s$-sparse function vanishes at all these points". Take any $2s$ basis functions and $2s$ points; the resulting square determinant, viewed as a function of the points, is not identically zero — that is exactly what the strong-linear-independence hypothesis asserts — so its zero set is null. There are only finitely many subsets of size $2s$, namely $\binom{N}{2s}$, a finite union of null sets is null, and any point of the complement works. **So $2s$ points always suffice, independently of $N$.**

The second question is constructive. The existence proof works by avoiding a null set and hands you no concrete point set; and the compressed-sensing theory underwriting practical $\ell_1$ recovery in high-dimensional Chebyshev bases uses **random** points, so its guarantees are probabilistic. The paper therefore exhibits an explicit **deterministic** point set and bounds its mutual coherence with Weil exponential sums. **The price is the $m\gtrsim s^2$ that the coherence route always charges; what is bought is the removal of the "with high probability" clause.**

### Setting

Let $\Omega\subset\mathbb R^d$ and let $\{B_j\}_{j\in\Lambda}$ be $N:=\#\Lambda$ complex-valued basis functions. The $s$-sparse class is

$$
U_s:=\Bigl\{f=\sum_{j\in T}c_jB_j\ :\ T\subset\Lambda,\ \#T\le s\Bigr\}.
$$

A set $\{x_1,\dots,x_m\}\subset\Omega$ is **unisolvent** for $U_s$ if $f(x_j)=g(x_j)$ for $j=1,\dots,m$ with $f,g\in U_s$ forces $f\equiv g$. The interpolation matrix is $A:=[B_j(x_t)]_{t=1,\dots,m;\,j\in\Lambda}$ with data $b:=[f(x_1),\dots,f(x_m)]^\top$.

### Derivation

**Step one: turn unisolvency into a measure-theoretic condition.** Call $f_1,\dots,f_k$ **strongly linearly independent** on $\Omega_0$ when the vanishing of $\sum_tc_tf_t$ on $\Omega_0$ forces $c_1=\cdots=c_k=0$. Lemma 3.1 proves this equivalent to the determinantal degeneracy set

$$
S:=\{x=(x_1,\dots,x_k)\in\Omega_0^{\,k}:\det(A_x)=0\},
\qquad A_x:=[f_t(x_j)]_{j,t=1,\dots,k},
$$

having zero outer Lebesgue measure, $\lambda^*_{d\cdot k}(S)=0$. **This is the mechanism of the whole existence argument**: it converts an algebraic condition into "$S$ is null", which can then be unioned and complemented.

**Step two: the union argument.** For each index set $T$ of size $2s$, the hypothesis makes $\{B_j\}_{j\in T}$ strongly linearly independent, so the corresponding $S_T$ is null. There are only $\binom{N}{2s}$ such $T$, a finite union stays null, hence $\Omega_0^{2s}\setminus\bigcup_TS_T$ is nonempty and any of its points gives $2s$ unisolvent interpolation points.

**Step three: the explicit point set.** The paper takes the Weil-sum grid also used in [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 9]],

$$
\Theta_M=\Bigl\{x_j=\cos(p_j)\ :\ p_j=\tfrac{2\pi}{M}\bigl(j,j^2,\dots,j^d\bigr),\ j=0,\dots,\lfloor M/2\rfloor\Bigr\},
\qquad m:=\#\Theta_M=\lfloor M/2\rfloor+1 .
$$

Stopping the index at $\lfloor M/2\rfloor$ is justified by Lemma 4.1: for any integer $M$ and $m=\lfloor M/2\rfloor+1$,

$$
\cos\bigl(2\pi j^k/M\bigr)=\cos\bigl(2\pi(M-j)^k/M\bigr)
\qquad\text{for all }k\in\mathbb N,\ 0\le j\le m-1,
$$

so the second half of the index range simply repeats the first and carries no new information.

**Step four: interpolation matrices and column normalisation.** With $\Phi_n(x_j)=C_n(p_j)$ and $C_n(p_j):=\prod_{t=1}^{d}\cos(2\pi n_tj^t/M)$, the matrix over the tensor-product index set $\Lambda^{q,d}_P$ is $A_P:=[C_n(p_j)]\in\mathbb R^{m\times(q+1)^d}$, and $A_D$ is its analogue over the total-degree set $\Lambda^{q,d}_D$. Recovery is performed on $A:=A_P\cdot C$ (or $A_D\cdot C$), where $C$ is the diagonal matrix standardising each column to unit $\ell_2$ norm. **That unremarkable column normalisation is the ancestor of the systematically designed preconditioners of [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|paper 24]].**

**Step five: the coherence chain.** **Weil's exponential sum theorem** for prime $p$ gives an $O(\sqrt p)$ bound on the relevant trigonometric sums; substituting into the definition of mutual coherence yields the explicit bound of Lemma 4.2; then $\delta_s\le(s-1)\mu(A)$ together with the recovery condition $\delta_s<\frac13$ is solved back for the requirement on $M$. The whole chain is explicit, with no probabilistic step.

### Theorems

**Theorem 3.2 (existence of $2s$ unisolvent points).** If $s\le N/2$ and there is a set $\Omega_0\subset\Omega\subset\mathbb R^d$ on which any $2s$ of the $\{B_j\}_{j\in\Lambda}$ are strongly linearly independent, then there exist $2s$ points $\{x_1,\dots,x_{2s}\}\subset\Omega_0$ that are unisolvent for $s$-sparse interpolation in the basis $\{B_j\}$. **The theoretical limit is $2s$, independently of $N$.**

**Theorem 3.5 (in a Chebyshev system, any $2s$ distinct points work).** The following are equivalent: (i) for $f,g\in U_s$ and **any** $2s$ distinct points $x_1,\dots,x_{2s}\in\Omega$, $f(x_j)=g(x_j)$ implies $f\equiv g$; (ii) for every index set $T$ with $\#T=2s$, the system $\{B_j\}_{j\in T}$ is a Chebyshev system. The paper also records that Chebyshev systems essentially do not exist on $[-1,1]^d$ for $d\ge2$, so **this clean univariate picture does not transfer to several variables** — which is why the second, deterministic-design half of the paper is needed at all.

**A negative result.** In a setting the paper exhibits, unisolvency forces $m\ge N$ in the other direction. **Sparsity alone does not always buy a reduction in the number of points.**

**Lemma 4.2 (coherence bound, tensor-product case).** For prime $M\ge\max\{2q+1,\ (2d(d-1))^2\}$,

$$
\mu(A_P)\ \le\ \frac{1}{\sqrt{M}}\cdot\frac{2^d\,d}{1-\dfrac{2d(d-1)}{\sqrt{M}}}.
$$

**Theorem 4.3 ($\ell_1$ recovery, tensor-product Chebyshev).** Let $M\ge\max\{2q+1,\ 9\cdot4^d\cdot d^2\cdot s^2\}$ be prime, let $f=\sum_{n\in\Lambda^{q,d}_P}c_n\Phi_n$, and let $c^{\#}$ solve the $\ell_1$ problem with matrix $A=A_P\cdot C$ and data $b=(f(x_1),\dots,f(x_m))^\top$, $x_j\in\Theta_M$, $m=\lfloor M/2\rfloor+1$. Then

$$
\|c^{\#}-c\|_2\ \lesssim\ \frac{\sigma_{s,1}(c)}{\sqrt{s}} .
$$

**Remark 4.4** gives the exact-recovery form: if $c$ is exactly $s$-sparse, recovery is exact once $m\ge\max\{q,\ \tfrac92\cdot4^d\cdot d^2\cdot s^2\}+1$ with $M$ prime.

**Lemma 4.5 and Theorem 4.6 (total-degree case).** Lemma 4.5 assumes $d\ge q$ and prime $M\ge\max\{2q+1,\ (2q(d-1))^2\}$. Theorem 4.6: for prime $M\ge9\cdot4^{\,q}\cdot d^2\cdot s^2$ and $f=\sum_{n\in\Lambda^{q,d}_D}c_n\Phi_n$, the $\ell_1$ solution with $A=A_D\cdot C$ satisfies the same bound $\|c^{\#}-c\|_2\lesssim\sigma_{s,1}(c)/\sqrt s$; Remark 4.7 gives exact recovery once $m\ge\tfrac92\cdot4^{\,q}\cdot d^2s^2+1$.

| Index set | Hypothesis of the coherence lemma | Points needed for exact recovery |
| ------------------------------- | ------------------------------------------------------ | -------------------------------------------- |
| Tensor product $\Lambda^{q,d}_P$ | $M$ prime, $M\ge\max\{2q+1,(2d(d-1))^2\}$ | $m\ge\max\{q,\ \tfrac92\cdot4^{d}d^2s^2\}+1$ |
| Total degree $\Lambda^{q,d}_D$ | $d\ge q$, $M$ prime, $M\ge\max\{2q+1,(2q(d-1))^2\}$ | $m\ge\tfrac92\cdot4^{\,q}d^2s^2+1$ |

**The only difference in that table that matters is the exponent moving from $4^d$ to $4^{\,q}$**: from exponential in the dimension to exponential in the polynomial degree. Uncertainty-quantification problems routinely have $d$ in the hundreds while $q$ stays small, so the total-degree version is the usable one.

**Remark 4.8** adds a geometric reading: when $4\mid M$, the first coordinates of $p_j$ over odd $j\in[1,M/2]$ form a set of Chebyshev nodes, so $\Theta_M$ can be seen as a high-dimensional extension of Chebyshev nodes.

### Numerical experiments

Section 5 compares $\Theta_M$ against interpolation points drawn i.i.d. uniformly on $[-1,1]^d$.

| Item | Setting |
| ---------- | -------------------------------------------------------------- |
| Comparison | the deterministic set $\Theta_M$ against i.i.d. uniform points |
| Spaces | two groups, tensor-product Chebyshev and total-degree Chebyshev |
| Test functions | support drawn uniformly over all size-$s$ subsets of $\Lambda$ |
| Coefficients | nonzeros i.i.d. standard Gaussian (mean 0, standard deviation 1) |
| Solver | SPGL1 in MATLAB |
| Statistics | 100 repetitions at each fixed sparsity $s$, giving an empirical success rate |

The reported outcome is that **the deterministic points perform similarly to random points**. What the experiment establishes is therefore not an accuracy advantage but the absence of an accuracy penalty for being deterministic — which is how the abstract itself puts it ("a similar performance").

> [!warning] There are no numbers here to transcribe
> The results are published as empirical success-rate curves; the material behind this page records the experimental configuration and the authors' qualitative conclusion but **no point-by-point success rates**. Do not infer any specific success probability from the table above. The gap to the theory is worth stating too: Theorem 4.3 demands $m\gtrsim4^dd^2s^2$, whereas the experiments use far fewer points, so what is being tested is how the method behaves outside the range the theorem covers, not the sharpness of the theorem.

### Relation to the others

Paper 10 is the sparse-recovery twin of [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 9]]: overlapping authors, the same Weil-sum point set $\Theta_M$, the same instinct to replace probabilistic guarantees with deterministic ones, but with $\ell_1$ minimisation and mutual coherence in place of least squares and Gershgorin. Its $m\gtrsim s^2$ is the quadratic barrier inherent to the coherence route, and the restricted-isometry route of papers 21, 24, 29 and 32 exists precisely to get around it. Xu and Zhou return to the same Weil machinery for trigonometric polynomials in paper **29**, there adding gradient measurements; the Gauss-sum computation $\mu(\Psi)=1/\sqrt p$ in Proposition 2.1 of paper 29 comes from the same source as the coherence estimates here. The column-normalisation matrix $C$ becomes the systematically designed preconditioner of [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|paper 24]].

## 21: randomly subsampling a Gauss grid, then minimising the ℓ1 norm

(Notation in this section: $M$ samples, $N=|\Lambda|$ basis functions.)

### The idea

Tang and Iaccarino had already shown that subsampling a tensor-product Gauss quadrature grid gives a good compressed-sensing design — but only for a **uniform** random variable on a hypercube, with sufficient sample count $M\gtrsim3^ds$. The open question was whether the distributions that actually matter in uncertainty quantification, especially the unbounded normal and exponential ones, work too.

The paper's key observation is an identity rather than a new construction: **the Gauss quadrature weights are exactly Christoffel function values**. Weighting the design matrix by $\sqrt w$ therefore does two things at once — it is the row normalisation that bounded-orthonormal-system theory demands, and it is the weight inherited naturally from the quadrature rule. After weighting, the family is orthonormal under the uniform empirical measure on the grid, and sampling i.i.d. from that measure just means "pick a grid point at random", which costs nothing to implement. **The whole question thus collapses to a purely analytic one: how large is the uniform bound $L$ for the Christoffel-weighted family?** For bounded (Beta) parameters $L$ is independent of the degree; for unbounded ones it grows like $n^{2/3}$. That is the paper's entire answer.

### Setting

Let $X=(X_1,\dots,X_d)^\top$ have mutually independent components with marginals $\rho_i$ on $\Gamma_i$ and joint density $\rho(x)=\prod_i\rho_i(x_i)$ on $\Gamma=\otimes_i\Gamma_i$. The univariate generalised polynomial chaos basis is defined by

$$
\mathbb E\bigl[\phi^i_n(X_i)\phi^i_\ell(X_i)\bigr]=\int_{\Gamma_i}\phi^i_n\phi^i_\ell\,\rho_i\,ds=\delta_{n,\ell},
$$

and the marginals may be Beta (Legendre / Chebyshev / Jacobi), normal (Hermite) or exponential (Laguerre). The index set $\Lambda$ is general; the tensor set is $\Lambda^P_{\mathbf n}=\{k\in\mathbb N_0^d:k\le\mathbf n\}$.

### Derivation

**The Gauss grid and its Christoffel weights.** Write $\Theta_{\mathbf n}=\Theta^1_{n_1}\otimes\cdots\otimes\Theta^d_{n_d}$ with $|\Theta_{\mathbf n}|=\prod_in_i$; it integrates every polynomial in the tensor space $\mathbb P_{2\mathbf n-1}$ exactly. Its weights are

$$
w_{\mathbf k}=\lambda_{\mathbf n}(z_{\mathbf k})
\triangleq\prod_{i=1}^{d}\lambda^i_{n_i}(z^i_{k_i})
=\prod_{i=1}^{d}\frac{1}{\sum_{k=0}^{n_i-1}\bigl[\phi^i_k(z^i_{k_i})\bigr]^2},
$$

**that is, the Gauss weights are exactly Christoffel function values** — in this page's convention, $1/K$ in tensorised form. The uniform empirical probability measure on the grid is

$$
\nu_{\mathbf n}=\bigotimes_{i=1}^{d}\nu^i_{n_i}
=\frac{1}{\prod_i n_i}\sum_{\mathbf k\le\mathbf n}\delta_{z_{\mathbf k}},
$$

so that i.i.d. sampling from $\nu_{\mathbf n}$ is equivalent to uniform sampling from the tensor Gauss grid.

**The weighted measurement matrix.** With $(\Psi)_{m,n}=\varphi_n(x_m)$ and $W$ diagonal, $(W)_{m,m}=w_m>0$, set $D=\sqrt W\Psi$. The unweighted problem is

$$
\arg\min\|c\|_1\quad\text{s.t.}\quad \Psi c=f,
$$

and the weighted one is

$$
\arg\min\|c\|_1\quad\text{s.t.}\quad Dc=\sqrt W f .
$$

**Discrete orthonormality.** The rows of $D$ are polynomials weighted by a Christoffel function. For any $\Lambda\subseteq\Lambda^P_{\mathbf n-1}$ with $N=|\Lambda|$, Lemma 4.1.D shows that the $N$ functions $\{\psi_{k,\mathbf n}(z)\}_{k\in\Lambda}$ are **orthonormal under the discrete measure $\nu_{\mathbf n}$**. Lemmas 3.1.A/B give the matrix-level version: $D^i=(\Sigma^i)^{1/2}\Psi^i$ is an orthogonal matrix and so is the tensor product $D=\bigotimes_iD^i$ — what Tang and Iaccarino called a "discrete orthogonal matrix". Lemma 3.1.C extends this to $\mathbf m\ge\mathbf n$, that is, more quadrature points than basis functions. **This is the concrete realisation of "row normalisation" in a tensor structure.**

**The algorithm.** (1) Given $\Lambda$, find $\mathbf n$ with $\Lambda\subseteq\Lambda^P_{\mathbf n-1}$; (2) generate the $n_i$-point Gauss rules — **the full tensor rule $\Theta_{\mathbf n}$ need never be constructed explicitly**; (3) choose $M$ points $\{(y_m,v_m)\}_{m=1}^M\subset\Theta_{\mathbf n}$ uniformly at random; (4) assemble $(D)_{m,n}=\sqrt{v_m}\,\varphi_{k(n)}(y_m)$ and $(W)_{m,m}=v_m$; (5) solve the weighted $\ell_1$ problem. Step (2) is what makes the design feasible in high dimension: of the $\prod_in_i$ grid points, only $M$ are ever generated.

### Theorems

**Three uniform-boundedness lemmas.** These are the technical core, since they determine the size of $L(\mathbf n)$.

| Marginal | Condition | $L_i(n)$ |
| ---------------------------------------------------------------------- | ------------------------------------ | ---------------------------------- |
| Beta / Jacobi, $B(\gamma+1,\delta+1)$, $\gamma,\delta\ge-\tfrac12$ | bounded interval $[-1,1]$ | $\le C(\gamma,\delta)$, **uniform in $n$** |
| two-sided exponential $\rho_i\propto e^{-\lvert x\rvert^\alpha}$, $\alpha>\tfrac32$ | $\mathbb R$ (normal is $\alpha=2$) | $\le C(\alpha)\,n^{2/3}$ |
| one-sided exponential $\rho_i\propto e^{-\lvert x\rvert^\alpha}$, $\alpha>\tfrac34$ | $[0,\infty)$ (exponential is $\alpha=1$) | $\le C(\alpha)\,n^{2/3}$ |

Remark 4.1 states that the conclusion is expected for the more general weight $\rho_i\propto x^\mu e^{-\lvert x\rvert^\alpha}$ with $\mu\ge-\frac12$, **but this is not proved here**.

**Theorem 4.1 (the main result).** Let $\Lambda$ be finite with $|\Lambda|=N$ and let $\mathbf n$ be the smallest multi-index with $\Lambda\subseteq\Lambda^P_{\mathbf n-1}$. Draw $M$ samples **without replacement** from $\nu_{\mathbf n}$. If

$$
M\ \ge\ L(\mathbf n)\,C_1\,s\log^3(s)\log(N),
\qquad
L(\mathbf n)=\prod_{i=1}^{d}L_i(n_i),
$$

with $C_1$ universal, then for any $c\in\mathbb R^N$ the solution $c^\sharp$ of the weighted $\ell_1$ problem satisfies

$$
\|c-c^\sharp\|_2\le C_2\,\frac{\sigma_{s,1}(c)}{\sqrt s}
$$

with probability at least $1-N^{-\gamma\log^3(s)}$, with $C_2,\gamma$ universal. **One easily overlooked detail of the setup: those $M$ samples are drawn from $\nu_{\mathbf n}$ without replacement**, not i.i.d. — the equivalence noted above concerns how the design is built, while the sampling step itself is without replacement.

> [!warning] The paper prints the probability inequality in the reversed direction
> The published statement of Theorem 4.1 reads
>
> $$
> \Pr\Bigl[\|c-c^\sharp\|_2\le C_2\tfrac{\sigma_{s,1}(c)}{\sqrt s}\Bigr]\ \le\ 1-N^{-\gamma\log^3(s)},
> $$
>
> that is, the probability is at **most** $1-N^{-\gamma\log^3 s}$. This contradicts the intended "with high probability" reading and also contradicts Theorem 2.2 of the same paper, which says "with probability at least $1-N^{-\gamma\log^3(s)}$". **This page uses the corrected direction $\Pr[\cdot]\ \ge\ 1-N^{-\gamma\log^3(s)}$** and records that the original prints it the other way. It is almost certainly a misprint, but cite the corrected direction.

**The three sample-complexity regimes**, in the paper's own summary form $M\gtrsim L(\mathbf n)s$:

| Parameter type | $L(\mathbf n)$ | Sample requirement |
| --------------------------------------------------- | ------------------------------ | ------------------------------------ |
| Beta (bounded) | $\le C^d$ | $M\gtrsim C^d\,s\log^3 s\log N$ |
| Normal (maximum degree $n-1$ the same in every dimension) | $\le(Cn)^{2d/3}$ | $M\gtrsim n^{2d/3}s\log^3 s\log N$ |
| One-sided exponential ($\rho\propto e^{-\lVert z\rVert_1}$ on $[0,\infty)^d$) | as above | as above |
| Normal, empirical observation at degree $n-1=9$ | $Cn^{2/3}\lesssim4$ | $M\ge4^d\,s\log^3 s\log N$ for $n\le10$ |

In the bounded case the constant of Tang and Iaccarino is essentially $C=3$; Figure 1 of the paper **suggests** that $C=2$ would be sharper, but that is **not proved**.

The paper is candid about the unbounded case: the $n^{2d/3}$ dependence "seems unpleasant", but it is **essentially sharp** given this analysis strategy and the insistence on subsampling a tensor-product Gauss grid. It then softens the statement empirically: at degree $n-1=9$ one observes $Cn^{2/3}\lesssim4$, so a degree-9 approximation with $n\le10$ requires $M\ge4^ds\log^3(s)\log(N)$, comparable to the bounded case — and high-dimensional problems use low degree anyway.

One closing remark deserves to be remembered on its own: one may subsample from $\nu_{\mathbf m}$ with $\mathbf m\ge\mathbf n$, and **for bounded $X$ the bound is unchanged**, but for exponential densities it degrades like $L\sim m^{2/3}>n^{2/3}$. The reason is concrete: larger Gauss rules place points where degree-$n$ polynomials weighted by $\lambda_{\mathbf m}$ decay rapidly to zero, and coefficients are hard to recover from there. **In other words, a larger grid is worse, not better — which is worth noting because it runs against intuition.**

The paper also restates Theorem 2.3 of Yan–Guo–Xiu: for multivariate Legendre on the total-degree space $T^d_n$, the preconditioning weights are $w_m=(\tfrac2\pi)^d\prod_{n=1}^{d}\bigl(1-(x^n_i)^2\bigr)^{1/2}$.

### Numerical experiments

The $\ell_1$ problems are solved with the Spectral Projected Gradient algorithm SPGL1 in MATLAB. Test coefficient vectors are built by fixing a sparsity $s$, drawing $s$ nonzero entries i.i.d. standard normal and zeroing the rest. Five recovery procedures are compared by empirical success probability:

| Label | Sampling | Problem solved |
| -------------- | ---------------------------------------------------- | -------------- |
| Random | i.i.d. from $\rho$ | unweighted |
| PreChebyshev | i.i.d. from the Chebyshev density $v(x)=\frac{1}{\pi^d\prod_i\sqrt{1-x_i^2}}$ | weighted (preconditioned) |
| **Gaussian** | subsampled Gauss grid with Gauss weights | weighted (**the method proposed here**) |
| Chebyshev | i.i.d. from the Chebyshev measure | unweighted |
| Uniform | i.i.d. from the uniform measure | unweighted |

Figure 2 plots recovery probability against sparsity $s$ for Legendre polynomials at fixed $M=85$, in two configurations:

| Configuration | $d$ | degree bound $n$ | $N$ | $M$ |
| ---------------------- | --- | ---------------- | --- | --- |
| low dimension, high degree | 2 | 21 | 231 | 85 |
| high dimension, low degree | 10 | 4 | 286 | 85 |

Figure 1 plots the theoretical bound $L$ itself: for Jacobi families with symmetric parameters $\gamma=\delta$, and for the two exponential-type densities (Hermite with $\rho=\exp(-x^2)$ on $\mathbb R$, Laguerre with $\rho=\exp(-x)$ on $[0,\infty)$).

> [!warning] What the experiments allow one to transcribe
> The configuration numbers above ($M=85$; $d=2,n=21,N=231$; $d=10,n=4,N=286$) are verifiable, but **the success-probability curves themselves are not transcribed here**. Note the gap to the theory as well: at $d=10$ Theorem 4.1 requires an $M$ carrying a factor like $C^{10}$, while the experiment uses $M=85$. **The experiments show that the method still works far below the theoretical sample count, not that the bound is sharp.** The paper likewise says only that Figure 1 "suggests" $C=2$, and draws no conclusion.

### Relation to the others

Paper 21 is the $\ell_1$ sibling of [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 13]]: the same randomised Gauss-quadrature subsampling design, with least squares there and sparse recovery here. It generalises the uniform-only result of Tang and Iaccarino to Beta, normal and exponential parameters. Its central observation — that Gauss weights **are** Christoffel function values, so subsampling a Gauss grid is implicitly Christoffel-weighted sampling — is developed systematically in [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|papers 22 and 45]], where the sampling density itself is taken proportional to the inverse Christoffel function rather than inherited from a quadrature rule.

All the dimension dependence sits inside $L(\mathbf n)$, the same structure as the $L(n)$ of [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|paper 24]], and the restricted-isometry / bounded-orthonormal-system machinery is shared with papers 10, 24, 29 and 32. The unbounded-domain difficulty it quantifies, of order $n^{2d/3}$, is the sparse-recovery manifestation of the difficulty [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 11]] met in the least-squares setting.

## Gradient enhancement: more rows, not more samples

Papers 29 and 32 share one starting point, worth stating before either set of technical details.

**One expensive forward solve usually returns the value and every partial derivative at once** — through a single adjoint solve, say, or by differentiating the whole solver automatically. Using only the values, each sample contributes one row of the design matrix; using the gradient as well, the same sample contributes $1+d$ rows while the number of unknowns is unchanged. The measurement matrix therefore goes from

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

with $W$ a diagonal weight matrix equalising the row norms. **It matters to be precise about where the gain lies: it is in rows per solve, not in samples.** The same $M$ sample points are used and the solver is still called $M$ times; what changes is how many linear equations are extracted from each call. **This turns "information per expensive sample" into an object one can design.**

Two things must then be handled.

First, **the row norms are no longer uniform**. Value rows and gradient rows have different magnitudes — gradient rows carry the scale of a derivative, and derivatives of orthogonal polynomials grow with the degree — so preconditioning stops being optional: without weighting, the restricted isometry constant or the mutual coherence is dominated by the worst class of rows. The numerics of paper 32 confirm this directly: **naively stacking gradient rows actually destroys the stability of the original matrix**. That is why paper 32 shares the "design sampling and preconditioning as a pair" stance with [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|paper 24]].

Second, **gradients are not free**. If the gradient requires an extra solve rather than coming as a by-product of an adjoint method, the factor $(1+d)$ in rows must be weighed against that cost. Paper 32 confronts this in its experiments with a "standard-double" baseline: give the standard method the same total budget of scalar data ($3N$ function values when $d=2$) and see whether the gradient method still wins.

The two papers differ in the **basis**, and the difference is essential:

- Paper 29 works in the Fourier basis. Since $D_{v}e^{ik\cdot z}=\langle v,k\rangle e^{ik\cdot z}$, **the gradient rows are scalar multiples of the value rows**, so the gradient block contributes no new randomness, only a known per-column rescaling. The coherence analysis is therefore purely algebraic — a refined Cauchy–Schwarz — and the conclusion is **deterministic**, valid for arbitrary sample points.
- Paper 32 works in an orthogonal polynomial basis. There $\frac{d}{dx}p^{(\alpha,\beta)}_n=c(n,\alpha,\beta)\,p^{(\alpha+1,\beta+1)}_{n-1}$ is a **different** family with a **different** orthogonality weight, and the gradient rows are not multiples of the value rows. Restoring isotropy requires a new preconditioner.

## 29: gradient-enhanced ℓ1 recovery of sparse trigonometric polynomials

(This section follows the paper: $N$ sample points, $M=\#\Gamma$ basis functions, $\Gamma\subset\mathbb Z^d$ the frequency index set, sparsity $s$. The paper writes $k$ both for the number of directions and for a frequency multi-index; this page renames the number of directions to $r$ to avoid the clash.)

### The idea

For sparse Fourier recovery, $\ell_1$ minimisation stands or falls on the mutual coherence $\mu(\Phi)$ or the restricted isometry constant. The paper asks two concrete questions: does appending **directional derivative** measurements at the same sample points strictly improve the recovery conditions? And — the more interesting question — can derivative measurements **replace** function values?

The mechanism behind the first is visible at a glance. After stacking, the column belonging to frequency $k$ becomes

$$
\bigl(a_k,\ \langle v_1,k\rangle a_k,\ \dots,\ \langle v_r,k\rangle a_k\bigr),
\qquad a_k=\bigl[e^{ik\cdot z_j}:z_j\in\Xi\bigr],
$$

so the inner product of two columns picks up a factor $1+\sum_t\langle v_t,k\rangle\langle v_t,k'\rangle$ and each norm picks up $1+\sum_t\langle v_t,k\rangle^2$. After normalisation the coherence is the largest cosine between the **augmented vectors** $(1,V^\top k)$ and $(1,V^\top k')$. As long as those two are not proportional, the angle is strictly larger than the original angle between $a_k$ and $a_{k'}$ — **and "not proportional" is exactly the admissibility condition the paper defines**. The improvement is thus purely geometric: adding the gradient lifts each frequency into a higher-dimensional space where distinct frequencies are pushed further apart.

**Prior work on gradient-enhanced $\ell_1$ was largely numerical, or (for Hermite polynomials with random samples) gave only an almost-sure strict inequality $\mu(\tilde\Phi)\le\mu(\Phi)$. The improvement claimed here is a deterministic bound $\mu(\tilde\Phi)\le\lambda\mu(\Phi)$ with an explicit $\lambda<1$ and no probabilistic qualifier.**

### Setting

Let $f(x)=\sum_{k\in\Gamma}c_ke^{ik\cdot x}$ on $x\in[-\pi,\pi)^d$ with $\Gamma\subset\mathbb Z^d$ finite, $M=\#\Gamma$; re-indexed, $f=\sum_{j=1}^Mc_j\varphi_j(x)$. Samples are $\Xi=\{z_j\}_{j=1}^N\subset[-\pi,\pi)^d$ and the interpolation matrix is $\Phi=(\varphi_t(z_j))\in\mathbb C^{N\times M}$. The support is $T=\{k:c_k\ne0\}$ with sparsity $s=\#T\ll M$, and $\Pi_s(\Gamma)=\bigcup_{\#T\le s}\Pi(T)$.

The paper poses two problems:

- **Problem 1:** recover $f\in\Pi_s(\Gamma)$ from $f(z_j)=f_j$ **and** $D_{v_t}f(z_j)=f'_{j,t}$ for $t=1,\dots,r$ with $r\le d$, where $D_{v}f(x)=\langle\nabla f(x),v\rangle$.
- **Problem 2:** recover $f$ from $D^{\tau_j}_{v_j}f(z_j)=y_j$, $j=1,\dots,N$, with $v_j\in\mathbb R^d$ and $\tau_j\in\mathbb Z_{\ge0}$ — that is, at each point one knows **either** a function value ($\tau_j=0$) **or** a directional derivative of some order, but not both.

### Derivation

**The gradient-enhanced matrix.** Set

$$
\Phi=\bigl(e^{ik\cdot z_j}\bigr)_{z_j\in\Xi,\,k\in\Gamma},
\qquad
\Phi_t=\bigl(D_{v_t}e^{ik\cdot z_j}\bigr)_{z_j\in\Xi,\,k\in\Gamma},
\qquad
\tilde\Phi:=\begin{pmatrix}\Phi\\ \Phi_1\\ \vdots\\ \Phi_r\end{pmatrix}.
$$

The key elementary identity is $D_{v_t}e^{ik\cdot z_j}=\langle v_t,k\rangle e^{ik\cdot z_j}$, so the columns of $\Phi_t$ are $\{\langle v_t,k\rangle a_k\}$. **The derivative rows are scalar multiples of the value rows — which is why the coherence analysis can avoid probability entirely.**

**The exact coherence ratio.** A direct computation gives

$$
\lambda=\max_{k\ne k'}\frac{\bigl|1+\sum_{t=1}^{r}\langle v_t,k\rangle\langle v_t,k'\rangle\bigr|}
{\sqrt{1+\sum_t\langle v_t,k\rangle^2}\ \sqrt{1+\sum_t\langle v_t,k'\rangle^2}},
$$

so the ratio in $\mu(\tilde\Phi)=\lambda\,\mu(\Phi)$ is exactly the supremum of the cosine between the augmented vectors.

**Admissibility — the paper's new hypothesis.** With $V:=(v_1,\dots,v_r)\in\mathbb R^{d\times r}$, define the distance

$$
\|k-k'\|_V:=\|V^\top k-V^\top k'\|_\infty .
$$

The directions $v_1,\dots,v_r$ are **admissible with respect to $\Gamma$** if $\|k-k'\|_V\ne0$ for all $k\ne k'$ in $\Gamma$, equivalently $V^\top k\ne V^\top k'$. If $\mathrm{span}\{v_1,\dots,v_r\}=\mathbb R^d$ then $V$ is admissible, but admissibility is **strictly weaker**: the paper's example takes $d=2$, $\Gamma=[-q,q]^2\cap\mathbb Z^2$ and $v_1=[1,\sqrt2]^\top$, so that $\{\eta\in\mathbb R^2:\langle\eta,v_1\rangle=0\}\cap\mathbb Z^2=\{0\}$ and **the single direction $v_1$ is already admissible even though it does not span $\mathbb R^2$**. Also set

$$
\Gamma_{\min}:=\min_{k\ne k',\ k,k'\in\Gamma}\|k-k'\|_V^2 .
$$

**A refined Cauchy–Schwarz (Lemma 2.1, quoted from the literature).** If $x,y\in\mathbb R^d$ are not proportional and $u$ satisfies $\langle u,x\rangle=0$ and $\langle u,y\rangle=1$, then

$$
\frac{\langle x,y\rangle^2}{\|x\|^2\|y\|^2}\le1-\frac{1}{\|y\|^2\|u\|^2}.
$$

Apply this to the augmented vectors $x=(1,V^\top k)$ and $y=(1,V^\top k')$: admissibility guarantees they are not proportional, so a $u$ meeting both constraints exists, and bounding $\|y\|^2$ and $\|u\|^2$ uniformly over $\Gamma$ produces the two maxima in the denominator of Theorem 2.1. (The material behind this page records the lemma's statement and the fact that it is the tool used, but **not the explicit construction of $u$**, so this step reproduces the shape of the argument rather than its every line.)

**The weighted variant.** Let $W$ be diagonal with value rows kept at weight 1 and **all derivative rows scaled by a single constant $\alpha$**. It is a preconditioner with exactly one free parameter, yet it improves the $q$-dependence of the coherence bound from $O(q^4)$ to $O(q)$, as recorded below.

**Column normalisation for Problem 2.** Set $\Psi=\tilde\Psi W$ with $\tilde\Psi:=\bigl(D^{\tau_j}_{v_j}\varphi_t(z_j)\bigr)\in\mathbb C^{N\times M}$ and

$$
W:=\mathrm{diag}\Bigl(\Bigl(\textstyle\sum_{j=1}^N|\langle v_j,k\rangle|^{2\tau_j}\Bigr)^{-1/2}:k\in\Gamma\Bigr)\in\mathbb C^{M\times M}.
$$

Solve $\arg\min\|c\|_1$ subject to $\Psi c=f$, then undo the scaling as $c^{\#}=Wc^{\#}$. The auxiliary quantities are

$$
Z_k:=\{j:\langle v_j,k\rangle^{\tau_j}=0,\ 1\le j\le N\},
\qquad
\kappa:=\max_{k\in\Gamma}\#Z_k,
\qquad
R_0:=\frac{\max_{j\in Z_k^c}|\langle v_j,k\rangle|^{\tau_j}}{\min_{j\in Z_k^c}|\langle v_j,k\rangle|^{\tau_j}},
$$

with the convention $\langle v_j,k\rangle^{\tau_j}=1$ when $\tau_j=0$. Here $Z_k$ collects the measurements that carry no information at all about frequency $k$, and $R_0$ measures how unevenly the measurements respond to a single frequency — **it is the constant that reappears, possibly large, in the sample complexity below**.

> [!warning] Two constants are typographically ambiguous in the full text used here
> First, the diagonal entry of $W$ above is printed as $1/\sum_j|\langle v_j,k\rangle|^{2\tau_j}$ with the extent of the radical impossible to determine from the extraction; since its stated purpose is to normalise the **column norms** of $\Psi$, it must be the inverse square root, which is how this page writes it. Second, the optimal $\alpha$ in Theorem 2.2 and its concrete value for $\Gamma=[-q,q]^d$ suffer from the same ambiguity about how far the radical extends. **Both are reconstructions; do not quote these constants without checking the typeset original.**

### Theorems

**Theorem 1.1 (the quoted recovery result).** If $\Phi\in\mathbb R^{N\times M}$ satisfies the $s$-order restricted isometry property with $\delta_s<\frac13$, then $\ell_1$ minimisation recovers. The paper notes that the older condition was $\delta_{3s}+3\delta_{4s}<2$ and that $\delta_s<\frac13$ is the recent improvement.

**Theorem 2.1 (the main coherence theorem).** If $V=(v_1,\dots,v_r)\in\mathbb R^{d\times r}$ is admissible with respect to $\Gamma$, then

$$
\mu(\tilde\Phi)\le\lambda\cdot\mu(\Phi),
\qquad
\lambda\le\left(1-\frac{\Gamma_{\min}}
{\max_{k\in\Gamma}\bigl(1+\|V^\top k\|^2_\infty\bigr)\ \max_{k\in\Gamma}\bigl(1+\|V^\top k\|^2_2\bigr)}\right)^{1/2}.
$$

**This is deterministic**: it holds for arbitrary sample points, with no probabilistic qualifier.

**Corollary 2.1.** For $\Gamma=[-q,q]^d\cap\mathbb Z^d$ and $V=(e_1,\dots,e_d)$ (the full gradient), $\Gamma_{\min}=1$, $\max_k\|V^\top k\|_\infty^2=q^2$ and $\max_k\|V^\top k\|_2^2=dq^2$, so

$$
\mu(\tilde\Phi)\le\left(1-\frac{1}{(1+q^2)(1+dq^2)}\right)^{1/2}\mu(\Phi).
$$

**Proposition 2.1 (admissibility is not decorative).** Take $d=2$, $\Gamma=[-q,q]^2\cap\mathbb Z^2$ and the deterministic quadratic-phase samples $\Xi=\{2\pi(j,j^2)/p:j=0,\dots,p-1\}$ with prime $p>2q+1$. The Gauss-sum formula gives $\mu(\Psi)=1/\sqrt p$, and appending the $e_1$-derivative block gives $\mu(\tilde\Psi)=\mu(\Psi)$ **exactly — no improvement whatsoever** — because $V=(e_1)$ is not admissible for the two-dimensional index set. **The proposition shows admissibility is a genuine hypothesis: without it, gradient information buys nothing in the coherence sense.**

**Theorem 2.2 and Corollary 2.2 (weighting strictly improves the bound).** Under the hypotheses of Theorem 2.1, the weighting above gives

$$
\mu(W\tilde\Phi)\le\lambda\mu(\Phi),
\qquad
\lambda\le\left(1-\frac{\alpha^2\Gamma_{\min}}
{\max_k\bigl(1+\alpha^2\|V^\top k\|_\infty^2\bigr)\ \max_k\bigl(1+\alpha^2\|V^\top k\|_2^2\bigr)}\right)^{1/2},
$$

minimised at $\alpha=\bigl(\max_{k\in\Gamma}\|V^\top k\|_\infty\cdot\max_{k\in\Gamma}\|V^\top k\|\bigr)^{-1/2}$, which is the second ambiguity flagged above. For $\Gamma=[-q,q]^d$, $V=(e_1,\dots,e_d)$ and $\alpha=1/(\sqrt d\,q)$,

$$
\mu(W\tilde\Phi)\le\left(1-\frac{1}{(1+\sqrt d)^2q}\right)^{1/2}\mu(\tilde\Phi),
$$

and Remark 2.1 verifies $\bigl(1-\frac{1}{(1+\sqrt d)^2q}\bigr)^{1/2}\le\bigl(1-\frac{1}{(1+q^2)(1+dq^2)}\bigr)^{1/2}$, so weighting strictly improves Corollary 2.1. **The substance of the improvement is the order in $q$: the denominator goes from $O(q^4)$ to $O(q)$.**

**Theorem 3.1 (restricted isometry for Problem 2).** Let $z_1,\dots,z_N$ be i.i.d. uniform on $[-\pi,\pi)^d$, and let $v_j\in\mathbb R^d$, $\tau_j\in\mathbb Z_{\ge0}$ satisfy $\sum_{j=1}^N\langle v_j,k\rangle^{2\tau_j}\ne0$ for every $k\in\Gamma$. For $0<\delta\le\frac12$ and $0<\epsilon<1$, if

$$
N\ \ge\ \frac{2(C_0R_0)^2}{\delta^2}\,s\,(\ln 100s)^2\,\ln(4M)\,\ln(10N)\,\ln\frac{\beta}{\epsilon}\ +\ \kappa,
$$

then $\mathbb P(\delta_s\le\delta)\ge1-\epsilon$, where $\delta_s$ is the $s$-order restricted isometry constant of $\Psi$, $M=\#\Gamma$, and $C_0,\beta$ are universal constants. Combined with Theorem 1.1 at $\delta=1/3$ this gives $\ell_1$ recovery with probability $\ge1-\epsilon$. **The authors' own reading: derivative evaluations play the same role as function evaluations, up to the constant $R_0$ — which the paper concedes "might be large".** The proof follows the Rudelson–Vershynin / Rauhut route: a symmetrisation lemma, a Rademacher chaos bound with the explicit literature constants $C_1'=94.81$, $C_2'\approx82.56$, $\beta=6.028$, and a moment-to-tail lemma.

**Corollary 3.1.** For $\Gamma=[1,q]^d\cap\mathbb Z^d$ with $q\ge2$, $\tau_j\in\{0,1\}$, $v_j\in\{e_1,\dots,e_d\}$ and uniform samples on $[-\pi,\pi)^d$, one has $\kappa=0$ and $R_0=q$, so the condition becomes

$$
N\ \ge\ \frac{2(C_0q)^2}{\delta^2}\,s\,(\ln100s)^2\,\ln(4M)\,\ln(10N)\,\ln\frac\beta\epsilon
\ \Longrightarrow\ \mathbb P(\delta_s\le\delta)\ge1-\epsilon .
$$

**Prior sample counts the paper records (citations, not new results).** Rauhut: uniform samples give $\delta_s\le\delta$ with probability $\ge1-\epsilon$ provided $N/\log N\ge C\delta^{-2}s\log^2(s)\log(M)\log(\epsilon^{-1})$. Kunis–Rauhut: $\mu<\frac{1}{2s-1}$ with probability $\ge1-\epsilon$ provided $N\ge C(2s-1)^2\log(4M^2/\epsilon)$. Xu (deterministic): for $\Gamma=[-q,q]^d$ with $d\ge2$, the points $z_j=2\pi(j,j^2,\dots,j^d)/N$ give $\mu(\Phi)<\frac{1}{2s-1}$ provided $N>\max\{(2s-1)^2(d-1)^2,\ 2q+1\}$ is prime.

> [!note] A conjecture the paper states explicitly
> Whether derivative information also **decreases the restricted isometry constant**, rather than only the coherence, is left as a **conjecture**: the numerics answer affirmatively but nothing is proved, for uniform samples on $[-\pi,\pi)^d$. Optimal choice of the directions $\{v_j\}$ and of $\alpha$ is also future work, as is extending from the Fourier basis to orthogonal polynomials — which is precisely what paper 32 does.

### Numerical experiments

The $\ell_1$ problems are solved with SPGL1 and every experiment uses $\Gamma=[-q,q]^d\cap\mathbb Z^d$; the labels "standard" and "gradient-enhanced" carry a percentage denoting the fraction of derivative components used.

| Example | $d$ | $q$ | Setting | Quantity reported |
| ------- | --- | --- | -------------------------------------------- | ---------------------------- |
| 4.1 | 2 | 10 | $v_1=e_1$, $v_2=e_2$; uniform samples; support drawn uniformly among size-$s$ subsets, nonzeros i.i.d. standard normal; 500 trials per configuration | recovery rate vs $N$ at $s=5$; vs $s$ at $N=20$ |
| 4.2 | 5 | 2 | 20% and 40% enhancement (partials in one and in two variables) | recovery rate vs $N$ at $s=6$; vs $s$ at $N=30$ |
| 4.3 | 2 | 5 | Problem 2 setting, $s=8$, $M=40$; $N/4$ function values and $3N/4$ values of $D_{e_1}f$, against standard $\ell_1$ with $N$ function values | recovery rate |
| 4.4 | 2 | 10 | function approximation: $g(x)=\sin\bigl(\cos\sum_jx_j\bigr)$ and $g(x)=\cos\bigl(\sum_jx_j\bigr)\exp\bigl(\sin(\sum_jx_j)\bigr)$, matching $f$ and $D_{e_t}f$ at the samples | discrete $L^2$ error vs $N$ |

The reported outcomes: in Examples 4.1 and 4.2 gradient information improves the recovery rate, **monotonically in how much of it is included**; in Example 4.3 **derivative values play a role similar to function values**, the empirical counterpart of Theorem 3.1; in Example 4.4 gradient information improves accuracy "dramatically".

> [!warning] Experimental data, and one configuration that does not add up
> No recovery rates or errors are transcribed here; the table records configurations and the authors' qualitative conclusions only. Separately, the configuration recorded for Example 4.3 — $q=5$, $d=2$, $M=40$ — is inconsistent with the paper's own convention $M=\#\Gamma$, since $\#([-5,5]^2\cap\mathbb Z^2)=121$. **Check that configuration against the typeset original.**
>
> The gap to the theory is also worth stating: Theorems 2.1 and 2.2 bound the coherence **ratio**, whereas Examples 4.1–4.2 measure recovery rates directly, and the step in between requires $\mu<\frac{1}{2s-1}$, which the experimental sizes $N=20$ and $N=30$ come nowhere near. **The experiments support the claim that gradients help; they do not establish sharpness of the theorems.**

### Relation to the others

Paper 29 is the trigonometric / Fourier twin of paper **32**, published the same year by an overlapping author set; paper 32 generalises the idea to preconditioned orthogonal polynomial chaos, which is exactly the extension paper 29 lists as future work. Both take the mutual-coherence route ($\mu<1/(2s-1)$), as does paper **10**. The Problem 2 formulation of paper 29 — mixed function and derivative data at **different** locations — reappears verbatim as an open problem in §4 of paper 32. The deterministic quadratic-phase samples $z_j=2\pi(j,j^2,\dots,j^d)/N$ it cites are precisely the **Weil points** of [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|papers 9 and 14]] and of paper 10 on this page; the Gauss-sum computation $\mu(\Psi)=1/\sqrt p$ in Proposition 2.1 is the same calculation that underlies the coherence estimates of paper 10.

## 32: gradient-enhanced ℓ1 minimisation for sparse polynomial chaos

(This section follows the paper: $N$ sample points, $M=\binom{d+n}{n}$ basis functions, $N\ll M$.)

### The idea

When only $N\ll M$ function evaluations are affordable, sparse polynomial chaos recovery via $\ell_1$ minimisation is the standard tool. Gradients $\partial_kf$ are often cheap — from adjoint solves, say — and appending them multiplies the row count by $d+1$ at little extra cost. **But naive stacking destroys the mean isotropy** $\mathbb E[\frac1N\Phi^\top\Phi]=I$ that underlies compressed-sensing guarantees, and the paper's own numerics show that it actually **degrades** the mutual coherence constant.

The root cause is one crucial difference from the Fourier basis: $\frac{d}{dx}p^{(\alpha,\beta)}_n=c(n,\alpha,\beta)\,p^{(\alpha+1,\beta+1)}_{n-1}$, so **the derivatives are a different Jacobi family with a different orthogonality weight, and the scaling constant $c(n,\alpha,\beta)$ grows with the degree**. The gradient rows are therefore neither multiples of the value rows nor orthonormal under the same measure. The paper's contribution is a general recipe: a **row preconditioner** $W$ that returns each block of rows to its own orthogonality weight, and a **column normalisation** $P$ that absorbs the growth of $c(n,\alpha,\beta)$, so that the composite matrix is isotropic again. The recipe applies on both bounded and unbounded domains.

### Setting

Let $x=(x_1,\dots,x_d)^\top$ have independent components with marginals $\rho_i$ on $\Gamma_i$ and $\rho(x)=\prod_i\rho_i(x_i)$ on $\Gamma=\otimes_i\Gamma_i$; univariate orthonormal $\phi^i_n$ satisfy $\int_{\Gamma_i}\phi^i_n\phi^i_\ell\rho_i=\delta_{n,\ell}$, and the multivariate basis is $\psi_n(x)=\prod_{i=1}^d\phi^i_{n_i}(x_i)$ with $\mathbb E[\psi_n\psi_j]=\delta_{n,j}$. The total-degree index set is $\Lambda^T_n=\{k\in\mathbb N_0^d:\sum_ik_i\le n\}$ of size $M=\binom{d+n}{n}$, and the expansion is $f_n=\sum_{j=1}^Mc_j\psi_j(x)$.

Standard $\ell_1$: with samples $\Xi=\{z^{(1)},\dots,z^{(N)}\}\subset\Gamma$ and $[\Phi]_{ij}=\psi_j(z^{(i)})$, $\Phi\in\mathbb R^{N\times M}$, solve $\arg\min\|c\|_1$ subject to $\Phi c=f$ (or $\|\Phi c-f\|_2\le\epsilon$). The mutual coherence $\mu(\Phi)$ and the criterion $\mu<\frac{1}{2s-1}$ are as described at the top of this page.

### Derivation

**The gradient-enhanced problem — the paper's central object.**

$$
\arg\min_{c\in\mathbb R^M}\|c\|_1
\quad\text{s.t.}\quad
W\tilde\Phi P\,c=W\tilde f,
$$

$$
\tilde f=\begin{pmatrix}f\\ f^\partial\end{pmatrix},\quad
\tilde\Phi=\begin{pmatrix}\Phi\\ \Phi^\partial\end{pmatrix},\quad
\Phi^\partial=\begin{pmatrix}\frac{\partial\Phi}{\partial x_1}\\ \vdots\\ \frac{\partial\Phi}{\partial x_d}\end{pmatrix},\quad
f^\partial=\begin{pmatrix}\frac{\partial f}{\partial x_1}\\ \vdots\\ \frac{\partial f}{\partial x_d}\end{pmatrix},
$$

with $\bigl[\frac{\partial\Phi}{\partial x_k}\bigr]_{ij}=\frac{\partial\psi_j}{\partial x_k}(z_i)$, so $\tilde\Phi\in\mathbb R^{N(d+1)\times M}$. **The two new matrices have distinct jobs: $W$ is a row preconditioner whose form depends on the chaos family and on how $\Xi$ was drawn, and $P$ is a column normalisation ensuring the composite $\hat\Phi:=W\tilde\Phi P$ is mean isotropic.**

**Legendre with Chebyshev sampling.** Take the product Chebyshev density $\rho_c(x)=\prod_{j=1}^d\frac{1}{\pi\sqrt{1-x_j^2}}$. Using the classical fact that derivatives of Legendre polynomials are orthogonal with respect to $\eta(x)=(1-x^2)$, the paper derives

$$
\mathbb E_c\Bigl[\frac{2^{-d}}{\rho_c(z)}\psi_i(z)\psi_j(z)
+\sum_{k=1}^d\frac{1-z_k^2}{\rho_c(z)}\frac{\partial\psi_i}{\partial x_k}(z)\frac{\partial\psi_j}{\partial x_k}(z)\Bigr]
=\delta_{ij}\Bigl(1+\sum_{k=1}^d c_k\,i_k(i_k+1)\Bigr).
$$

**This identity is the template for the whole design**: the two terms on the left carry the orthogonality weight of the value rows and of the gradient rows respectively ($2^{-d}$ for the uniform density, $1-z_k^2$ for the derivative family), the $\delta_{ij}$ says orthogonality survives the weighting, and the bracketed factor is exactly the degree-growing quantity that $P$ must divide out. Accordingly

$$
W^0_{n,n}=\Bigl(\tfrac{4}{\pi^2}\bigl(1-(z^{(n)}_j)^2\bigr)\Bigr)^{d/4},
\qquad
W^j_{n,n}=\frac{W^0_{n,n}}{\sqrt2}\bigl(1-(z^{(n)}_j)^2\bigr)^{1/2},
\qquad
P_{i,i}=\Bigl(1+\sum_{k=1}^d c_k\,i_k(i_k+1)\Bigr)^{-1/2},
$$

with $W$ block diagonal, $W=(W^0,W^1,\dots,W^d)$, giving $\mathbb E_c\bigl[\frac1N\hat\Phi^\top\hat\Phi\bigr]=I$.

> [!warning] Do not copy the two constants of the Legendre case
> The leading factor in the identity above is printed as "$2d$" in the full text used here, plainly having lost an exponent; the general Jacobi formula below shows it must be the density ratio $\rho/\rho_c$ with $\rho\equiv2^{-d}$ uniform, which is how this page writes it — **a reconstruction**. The expression for $W^0_{n,n}$ involves only the $j$th component yet carries the exponent $d/4$, which strongly suggests the extraction dropped a $\prod_{j=1}^d(\cdot)^{1/4}$. **Do not rely on either printed constant**; the general Jacobi version in the next paragraph is clean and internally consistent and should be used instead.

**General Jacobi with Chebyshev sampling — the clean general form.** The univariate Beta / Jacobi density is

$$
\rho^{(\alpha,\beta)}(x)=d^{(\alpha,\beta)}(1-x)^\alpha(1+x)^\beta,
\quad \alpha,\beta\ge-\tfrac12,
\qquad
d^{(\alpha,\beta)}=\frac{\Gamma(\alpha+\beta+2)}{\Gamma(\beta+1)\Gamma(\alpha+1)2^{\alpha+\beta+1}},
$$

with $\rho_c\equiv\rho^{(-1/2,-1/2)}$ and the multivariate density taken as a product. Then for $z\sim\rho_c$,

$$
\mathbb E\Bigl[\frac{\rho^{(\alpha,\beta)}(z)}{\rho_c(z)}\psi_i\psi_j
+\sum_{k=1}^d\frac{\rho^{(\alpha+e_k,\beta+e_k)}(z)}{\rho_c(z)}
\frac{\partial\psi_i}{\partial x_k}\frac{\partial\psi_j}{\partial x_k}\Bigr]
=\delta_{ij}\Bigl(1+\sum_{k=1}^d c^2(i_k,\alpha_k,\beta_k)\Bigr),
$$

with $e_j$ the $j$th cardinal vector, $e_0=0$, and the explicit normalisation constant

$$
c^2(i_k,\alpha_k,\beta_k)=\frac{i_k(i_k+\alpha_k+\beta_k+1)(\alpha_k+\beta_k+2)(\alpha_k+\beta_k+3)}{4(\alpha_k+1)(\beta_k+1)} .
$$

**Note that $c^2$ grows linearly in the degree $i_k$ — that is the quantitative reason the gradient rows would swamp the value rows without $P$.** Hence

$$
W^0_{n,n}=\sqrt{\frac{\rho^{(\alpha,\beta)}(z^{(n)})}{\rho_c(z^{(n)})}},
\qquad
W^j_{n,n}=\sqrt{\frac{\rho^{(\alpha+e_j,\beta+e_j)}(z^{(n)})}{\rho_c(z^{(n)})}},
\qquad
P_{i,i}=\Bigl(1+\sum_{k=1}^d c^2(i_k,\alpha_k,\beta_k)\Bigr)^{-1/2},
$$

which gives $\mathbb E\bigl[\frac1N\hat\Phi^\top\hat\Phi\bigr]=I$.

**The general recipe, as the paper states it.** (i) Draw $\Xi$ from a sampling measure that is **degree-asymptotically good** for the chaos basis; (ii) design $W$ so that the chaos basis is mean isotropic under that measure; (iii) choose $P$ so that the isotropy of the **gradient** rows is retained as well.

**Coherence parameters.** The standard approach uses $\mu_L(\Phi):=\sup_{i,\,z\in\Xi}|\Phi_i(z)|_2^2$; the gradient-enhanced approach uses

$$
\beta_L(\hat\Phi):=\sup_{i,\,z\in\Xi}\bigl\|\hat\Phi_i(z)\bigr\|_2,
\qquad
\hat\Phi_i(z)=\frac{1}{P_{i,i}}
\begin{pmatrix}
\sqrt{\rho^{(\alpha,\beta)}(z)/\rho_c(z)}\ \Phi_i(z)\\
\sqrt{\rho^{(\alpha+e_1,\beta+e_1)}(z)/\rho_c(z)}\ \partial_{x_1}\Phi_i(z)\\
\vdots\\
\sqrt{\rho^{(\alpha+e_d,\beta+e_d)}(z)/\rho_c(z)}\ \partial_{x_d}\Phi_i(z)
\end{pmatrix}.
$$

**Hermite with Gaussian sampling.** Derivatives of Hermite polynomials are orthogonal under the **same** Gaussian measure, so this case is much cleaner: for suitably normalised $\psi_j$ and standard multivariate normal $z$,

$$
\mathbb E\Bigl(\psi_i(z)\psi_j(z)+\sum_{k=1}^d\frac{\partial\psi_i}{\partial x_k}(z)\frac{\partial\psi_j}{\partial x_k}(z)\Bigr)
=\delta_{ij}\Bigl(1+\sum_{k=1}^d i_k\Bigr),
$$

so one takes $P_{i,i}=\bigl(1+\sum_ki_k\bigr)^{-1/2}$ and **$W=I$** — **no row preconditioner is needed at all on the unbounded domain**. (The paper attributes this case to its reference [35], where a result similar to its own Theorem 3.1 was shown.)

**A structural remark.** $P$ can be taken diagonal only because the derivatives of the classical families are again orthogonal families, and the paper notes that **Jacobi, Laguerre and Hermite are the only univariate families with that property**. For a non-classical chaos basis, $P$ must instead be any inverse square root of the Gramian of the polynomial derivatives, and is **no longer diagonal**.

### Theorems

**Theorem 3.1 (the main theorem).** Let $\Phi$ and $\hat\Phi$ be the design matrices for standard and gradient-enhanced $\ell_1$ via Jacobi expansions with Chebyshev sampling. Then

$$
\mu_L(\Phi)\le\prod_{j=1}^d 2e\Bigl(2+\sqrt{\alpha_j^2+\beta_j^2}\Bigr),
\qquad
\beta_L(\hat\Phi)\le C\prod_{j=1}^d 2e\Bigl(2+\sqrt{\alpha_j^2+\beta_j^2}\Bigr),
$$

where

$$
1\le C\le 1+\frac{\sqrt2}{2}\approx1.707 .
$$

The lower bound $C=1$ is attained when $\alpha_k=\beta_k=-\frac12$ for every $k$, and the upper bound occurs when some $k$ has $\alpha_k=\beta_k=0$. Furthermore $\mathcal N(\hat\Phi)\subset\mathcal N(\Phi)$, and this is **almost surely a strict subset when $\Phi$ is undersampled** — the precise sense in which adding gradients cannot hurt and generically helps.

**Supporting Lemma A.1 (a quoted Rauhut–Ward-type uniform Jacobi bound).** For all Jacobi weights with $\alpha\ge-\frac12$ and $\beta\ge-\frac12$,

$$
\sup_{x\in[-1,1]}\frac{\rho^{(\alpha,\beta)}(x)}{\rho_c(x)}\bigl[p^{(\alpha,\beta)}_n(x)\bigr]^2
\le 2e\Bigl(2+\sqrt{\alpha^2+\beta^2}\Bigr),
$$

uniformly in $n,\alpha,\beta$. The proof of Theorem 3.1 applies this bound to the two-row structure $\bigl(p_n^{(\alpha,\beta)}(z),\ \frac{d}{dx}p_n^{(\alpha,\beta)}(z)\bigr)=\bigl(p_n^{(\alpha,\beta)}(z),\ c(n,\alpha,\beta)p_{n-1}^{(\alpha+1,\beta+1)}(z)\bigr)$.

> [!warning] What the theorem does not deliver — stated candidly by the authors
> The paper writes: "Ideally we could show the gradient approach admits an improved (smaller) parameter $\beta_L$, i.e. $\beta_L(\hat\Phi)\le\mu_L(\Phi)$… Our analysis does not bear this fruit." What is actually shown is two things: (i) the coherence bound for both $\Phi$ and $\hat\Phi$ is a constant raised to the $d$th power and is **independent of the polynomial degree**; (ii) the extra factor $C\le1+\frac{\sqrt2}{2}$ is **dimension-independent and small**.
>
> One more point matters as much: **there is no sample-complexity theorem of the form $N\gtrsim s\log^k$ in this paper**, the guarantee running instead through the quoted criterion $\mu<\frac{1}{2s-1}$. Any claim about a sample-complexity exponent for this method has no source.

**Open directions the paper raises.** Partial or directional gradient data $D_{v_t}f(z^{(j)})=\langle\nabla f,v_t\rangle$ (with the honest caveat that it is unclear how such directional derivatives would be obtained in practice); and higher-order directional data $D^{\tau_j}_{v_j}f(z_j)=y_j$ where value locations and derivative locations **need not coincide** — which is exactly Problem 2 of paper 29. The paper also states that its preconditioners are **not claimed optimal** and suggests the Christoffel-weighted approach of [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|papers 22 and 24]] as a degree-asymptotically better alternative.

### Numerical experiments

The $\ell_1$ problems are solved with SPGL1. Three method labels run through all experiments: *standard* ($\ell_1$ on $\Phi$), *gradient-enhanced* (on $\hat\Phi$), and ***standard-double***, meaning the standard approach given **the same total budget of scalar data**: in $d=2$ the 100% gradient method uses $3N$ pieces of information ($N$ values plus $2N$ derivatives), so standard-double uses $3N$ function values. Percentages denote the fraction of derivative components used (50% in $d=2$ means $N$ values plus $N$ randomly chosen partials).

| Section | Configuration | Quantity reported |
| ------- | ---------------------------------------- | ------------------------------------- |
| §5.1 stability | $(d,n)=(2,30)$ and $(6,5)$ | mutual coherence of $\Phi$, $\tilde\Phi$, $\hat\Phi$ vs sample count; and vs $M$ at fixed $N=80$ |
| §5.2 fixed sparsity | $(d,n)=(2,20)$; uniform inputs, Chebyshev samples; exactly $s$-sparse $c$ with nonzeros i.i.d. standard normal; 100 trials; success when $\lVert c-\tilde c\rVert_\infty\le10^{-3}$ | recovery probability vs $N$ at $s=8$; vs $s$ at fixed $N$ |
| §5.2 fixed sparsity | $(d,n)=(10,3)$; 10% and 20% gradient enhancement (one and two partials) | recovery probability vs $N$ at $s=6$; vs $s$ at $N=70$ |
| §5.3 function approximation | $(d,n,M)=(2,20,231)$, $(10,3,286)$ and $(6,5)$; Legendre chaos with Chebyshev samples | discrete $L^2$ / root-mean-square error vs $N$ |

The three test functions of §5.3 are the sphere $f_1(x)=\sum_{i=1}^dx_i^2$, the Gaussian $f_2(x)=\exp\bigl(-\sum_{i=1}^d0.01(\tfrac12(x_i+1)-0.375)^2\bigr)$ and the sinusoids $f_3(x)=\sum_{i=1}^d0.3+\sin(\tfrac{16}{15}x_i-0.7)+\sin^2(\tfrac{16}{15}x_i-0.7)$.

**§5.1 is the most important experiment in the paper and the reason the whole preconditioning apparatus exists**: the preconditioned $\hat\Phi$ has a far better mutual coherence than $\Phi$, while **naively stacking derivative rows to form $\tilde\Phi$ actually destroys the stability $\Phi$ had**. Sections 5.2 and 5.3 conclude that including gradients improves the recovery rate, that more gradient information gives better recovery, and that approximation accuracy improves "dramatically".

> [!warning] The text and the figure captions disagree about sample counts
> In §5.2 the fixed sample count for the "recovery probability vs $s$" plot at $(d,n)=(2,20)$ is **$N=35$ in the text and $N=50$ in the figure caption**; the corresponding number at $(d,n)=(10,3)$ is **$N=70$ in most places and $N=50$ in one**. This page records both inconsistencies without picking a side; the table above uses the more frequently occurring value.
>
> On the gap to the theory: Theorem 3.1 bounds only the coherence parameters $\mu_L$ and $\beta_L$ and says explicitly that $\beta_L\le\mu_L$ was not proved, whereas the experiments measure recovery probabilities and approximation errors directly. **The experiments establish that the method works, not that the theorem predicted it would.** No specific success rates or errors are transcribed here.

### Relation to the others

Paper 32 is the polynomial-chaos, multivariate counterpart of paper **29**: the same authors' circle, the same gradient-enhanced $\ell_1$ idea, but with orthogonal polynomial chaos and preconditioning in place of trigonometric polynomials and random frequencies. It builds directly on the preconditioning framework of [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|paper 24]]: Chebyshev sampling with $\sqrt{\rho^{(\alpha,\beta)}/\rho_c}$ weighting is exactly CSA-a of paper 24 and the Rauhut–Ward asymptotic sampling restated in [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 14]], and the paper explicitly says the Christoffel-weighted schemes of papers 22 and 24 would be a degree-asymptotically better alternative to the preconditioners it uses. The mutual-incoherence route it takes is the one paper **10** takes for deterministic Chebyshev point sets, in contrast to the restricted-isometry route of papers **21** and 24. Its Lemma A.1, the uniform Jacobi bound $2e(2+\sqrt{\alpha^2+\beta^2})$, belongs to the same family of bounds cited in papers 21 and 24.

## Data-driven chaos: the basis has to be estimated too

Standard polynomial chaos assumes the input distribution is known analytically, so the orthogonal basis is given by the orthogonal polynomials of that distribution. In the **data-driven** setting the input distribution is given only through finitely many samples: there is no analytic density and therefore no ready-made orthogonal basis.

Papers 36 and 44 treat this setting, and the point is the same in both: **the orthogonal basis must be built from the empirical measure** — by Gram–Schmidt, or, more stably, by reading three-term recurrence coefficients off a Cholesky factorisation of the Hankel matrix of moments — and the sampling design changes accordingly, because induced sampling and Christoffel weighting both need a basis, **and the basis is now itself estimated**.

That introduces an error source absent from the rest of this topic: **the estimation error of the basis**. It compounds with the approximation error and the sampling error, so analysis in the data-driven setting carries one extra layer compared with the known-distribution case. Neither paper resolves that layer; both concede it. Paper 36 says plainly that a "density error" is introduced because the mean and variance are estimated, and leaves quantifying and controlling it to future work.

The division of labour is: paper 36 uses **weighted least squares** with **equilibrium-measure** sampling (asymptotically optimal), and paper 44 uses **preconditioned $\ell_1$** with **induced sampling** (exact at every finite degree, and supported only on the data).

## 36: data-driven polynomial chaos by weighted least squares

(Notation in this section: $M$ samples, $N$ data-driven basis functions, maximum degree $k$.)

### The idea

Classical generalised polynomial chaos presumes the input density $\rho$ is known — aleatory uncertainty. In many applications one has only **samples or moments** of the inputs, which is epistemic uncertainty, and assuming a wrong density biases the entire surrogate. Data-driven or arbitrary polynomial chaos (aPC) solves the first half by building the orthogonal basis directly from moments, but the existing post-processing was a sparse-grid collocation **whose points are themselves computed from the data-driven basis by matrix operations** — so the design must be recomputed every time the basis changes.

The paper's contribution is to replace that post-processing with one whose **sampling depends neither on the input density nor on the data-driven basis**, and can therefore be done offline, while improving stability from quadratic to quasi-linear in $N$. The mechanism is the one that recurs on this page: Christoffel weighting drives the stability factor $\kappa(N)=\max_\xi\sum_j\Phi_j^2(\xi)$ down to its theoretical minimum $N$. The price is that the weighted basis is no longer orthogonal with respect to $\rho$ but with respect to a transformed measure that depends on the polynomial space — and the paper removes that dependence by passing to a potential-theoretic limit (the equilibrium measure), **so stability holds only in the asymptotic sense $N\to\infty$**.

### Setting and derivation

**The Hankel matrix of moments and the moment-matching basis.** Given raw moments $\mu_0,\dots,\mu_{2k}$, set

$$
H=\begin{pmatrix}
\mu_0&\mu_1&\cdots&\mu_k\\
\mu_1&\mu_2&\cdots&\mu_{k+1}\\
\vdots&\vdots&\ddots&\vdots\\
\mu_k&\mu_{k+1}&\cdots&\mu_{2k}
\end{pmatrix}.
$$

When the moments come from $M$ samples one requires the set to be **determinate in the Hamburger sense**, that is $\det(H)>0$. Cholesky-factor $H=R^\top R$ with $R$ upper triangular; by the **Mysovskih theorem** the entries of $R$ furnish an orthogonal polynomial system, with the three-term recurrence

$$
\eta\,\phi_{j-1}(\eta)=b_{j-1}\phi_{j-2}(\eta)+a_j\phi_{j-1}(\eta)+b_j\phi_j(\eta),\quad j=1,\dots,k,
$$

$$
a_j=\frac{r_{j,j+1}}{r_{j,j}}-\frac{r_{j-1,j}}{r_{j-1,j-1}},
\qquad
b_j=\frac{r_{j+1,j+1}}{r_{j,j}},
\qquad r_{0,0}=1,\ r_{0,1}=0 .
$$

Multivariate bases are formed by tensorisation. The paper notes that the alternative route, inverting a Vandermonde system, is badly conditioned for large $k$.

**Weighted least squares.** With data-driven bases $\{\Phi_j\}_{j=1}^N$ spanning $\mathbb P_N$ (total-degree type, maximum order $k$),

$$
f_N:=\arg\min_{p\in\mathbb P_N}\frac1M\sum_{m=1}^M w_m\bigl(p(z_m)-f(z_m)\bigr)^2
\ \Longleftrightarrow\
c=\arg\min_{c\in\mathbb R^N}\bigl\|W^{1/2}Ac-W^{1/2}f\bigr\|_2^2,
$$

with $A=[\Phi_j(z_m)]\in\mathbb R^{M\times N}$ and $W=\mathrm{diag}(w_1,\dots,w_M)$. The weights are normalised Christoffel function values: in this page's convention $K(\xi)=\sum_{j=1}^N\Phi_j^2(\xi)$,

$$
w_m=\frac{N}{K(z_m)}=\frac{N}{\sum_{j=1}^N\Phi_j^2(z_m)} .
$$

**The five-step algorithm.** (1) Sample from the probability density $\hat\rho$ of an **equilibrium measure** (which depends on the input density $\rho$ and, on unbounded state spaces, on the maximum degree $k$, written $\hat\rho_k$); (2) evaluate $f$ at the samples; (3) form the $M\times N$ Vandermonde-like matrix $A$; (4) form the diagonal $W$ from the normalised Christoffel function; (5) solve the weighted least-squares problem. **Note the ordering of steps (1) and (3)–(4): sampling comes before the basis, which is exactly what "offline sampling" means.**

**On bounded domains the sampling measure is universal.** On $[-1,1]^d$, **whatever the input measure is, even if unknown**, sample from the tensor Chebyshev density

$$
\hat\rho(\xi)\sim\frac{1}{\pi^d\prod_{k=1}^d\sqrt{1-\xi_k^2}} .
$$

The only information required is that the random variable lives in a bounded domain — the paper calls the Chebyshev measure universal in the bounded setting. Sampling is trivial: draw $u_m$ uniform and set $z_m=\cos(u_m)$.

**Unbounded domains: conjectures with explicit samplers.** The paper says plainly that very few results are known and that "the results in what follows are our conjectures", validated only numerically.

> [!warning] The unbounded-domain equilibrium measures are conjectures, not theorems
> Both densities below are labelled **conjectures** in the original, and can only be cited that way. The Gaussian limiting induced / equilibrium measure $C(2-\|\xi\|^2)^{d/2}$ recurs across several papers in this topic, and **it is a conjecture in every one of them, never a theorem**; the Chebyshev limit on bounded domains, by contrast, is backed by a theorem. Do not conflate the two.

*(Gaussian, $\mathbb R^d$)*: estimate $(\hat\mu,\hat\sigma)$ from the data, standardise $\hat\xi=(\xi-\hat\mu)/\hat\sigma$, and use

$$
\hat\rho(\xi)=C\bigl(2-\|\xi\|^2\bigr)^{d/2}
$$

with $C$ a normalising constant, samples expanded by $\sqrt k$. The concrete sampler: compute $k$; draw $y=(y_1,\dots,y_d)$ i.i.d. standard normal; draw $\nu\sim\mathrm{Beta}(\alpha=d/2,\ \beta=d/2+1)$ on $[0,1]$; set

$$
z=\frac{y}{\|y\|_2}(2k\nu)^{1/2},
$$

producing samples on the Euclidean ball of radius $\sqrt{2k}$.

*(Exponential, $\mathbb R_+^d$)*:

$$
\hat\rho(\xi)=C\sqrt{\frac{\bigl(4-\sum_{i=1}^d\xi_i\bigr)^d}{\prod_{i=1}^d\xi_i}}
$$

(the placement of the radical over the whole fraction is **reconstructed** from the extraction and flagged as such). The sampler: compute $k$; draw a $(d+1)$-dimensional Dirichlet vector $y$ with parameters $\bigl(\frac12,\frac12,\dots,\frac12,\frac d2+1\bigr)$; discard the last entry; set $z=4ky$. Remark 4.1 concedes that for other unbounded densities **even conjectures are lacking**, and suggests domain truncation plus Chebyshev sampling as a fallback, with the truncation error left unaddressed.

### Theorems

**Theorem 3.1 (quoted: when the moment problem determines the basis).** Under Assumption 1 (all moments finite) and Assumption 2 (continuous distribution function $F_\eta$), if **any one** of the following holds, the moment problem is uniquely solvable and the moment-matched polynomials are dense in $L^2(\Omega,\sigma(\eta),P)$: (1) $F_\eta$ has compact support; (2) $\liminf_{k\to\infty}\frac{\sqrt[2k]{\mu_{2k}}}{2k}<\infty$; (3) exponential integrability, $\int_{\mathbb R}e^{a|x|}F_\eta(dx)<\infty$ for some $a>0$, equivalently a finite moment generating function near the origin; (4) **Carleman's condition** $\sum_{k=0}^\infty\frac{1}{\sqrt[2k]{\mu_{2k}}}=\infty$; (5) **Lin's condition**: the distribution has a symmetric, differentiable, strictly positive density $f_\eta$ and for some $x_0>0$, $\int_{-\infty}^\infty\frac{-\log f_\eta(x)}{1+x^2}dx=\infty$ and $-\frac{xf'_\eta(x)}{f_\eta(x)}\nearrow\infty$ as $x\to\infty$ for $x\ge x_0$. The paper notes the theory is open for more general settings such as discrete distributions, **even though the method still runs there** — and indeed the experiments use binomial and Poisson inputs.

**Theorem 4.1 (quoted least-squares stability).** For $f$ approximated in $\mathbb P_N=\mathrm{span}\{\Phi_j\}$ with orthogonality density $\rho$, samples $\{z_m\}_{m=1}^M$ drawn from $\rho$, and **unweighted** least squares $c=\arg\min\|Ac-f\|_2^2$,

$$
\Pr\Bigl\{\|\,\cdot\,-I\|\ge\tfrac12\Bigr\}\le 2M^{-r}
\qquad\text{provided}\qquad
\kappa(N):=\max_\xi\sum_{j=1}^N\Phi_j^2(\xi)\ \le\ \delta\,\frac{M}{\log M}.
$$

> [!warning] Theorem 4.1 needs two corrections
> First, the object inside the norm is printed as "$A$", but the content of the theorem requires it to be the **Gram matrix $A^\top A$** or a scaled version of it, and certainly not the rectangular design matrix; this page therefore writes a placeholder there. Second, the constant $\delta$ is printed as $1-\frac{\log2}{2-2r}$, whose fraction structure is ambiguous in the extraction used here; the corresponding constant in the Cohen–Davenport–Leviatan source has the form $\frac{1-\log2}{2+2r}$. **Do not quote this $\delta$ without checking the typeset original.**

**The motivating consequence — the paper's key argument.** Stability requires $M\gtrsim\kappa(N)$ up to a logarithmic factor. For Legendre polynomials $\kappa(N)\sim N^2$, so $M\ge CN^2$ — **unsatisfactory**. Introducing $W$ amounts to working with the rescaled basis

$$
\hat{\mathbb P}_N=\mathrm{span}\Bigl\{\hat\Phi_j=\sqrt{\tfrac{N}{K(\xi)}}\,\Phi_j\ \Big|\ 1\le j\le N\Bigr\},
\qquad
\hat\kappa(N):=\max_\xi\sum_{j=1}^N\hat\Phi_j^2(\xi)\equiv N,
$$

which is **optimal control of the stability factor** and gives the quasi-linear $M\gtrsim N\log N$. But the rescaled basis is orthogonal not under $\rho$ but under the transformed measure

$$
\tilde\rho(\xi)\ \propto\ K(\xi)\rho(\xi)=\Bigl(\sum_{j=1}^N\Phi_j^2(\xi)\Bigr)\rho(\xi),
$$

**which is exactly the induced measure**. It depends on the polynomial space and is not easy to sample. The paper's resolution comes from potential theory:

$$
\tilde\rho(\xi)\ \longrightarrow\ \hat\rho(\xi)\qquad (N\to\infty),
$$

so one samples from the equilibrium measure $\hat\rho$ instead and obtains stability **only in the asymptotic sense $N\to\infty$**. The practical payoff is that the design becomes independent of the polynomial space, which is valuable for adaptive schemes.

> [!warning] Equations (4.5) and (4.7) use $K$ in a way that conflicts with (4.3)
> The paper writes $K$ for $N/\sum_j\Phi_j^2$, the normalised Christoffel function, and the weights $w_m=N/\sum_j\Phi_j^2(z_m)$ of (4.3) are consistent with that. But (4.5) is printed as $\hat\Phi_j=\Phi_j/\sqrt{K}$ and (4.7) as $\tilde\rho\sim K\rho=N\rho/\sum_j\Phi_j^2$, and **read with that same $K$ those two give $\hat\kappa\ne N$ and do not give the induced measure** — they are self-consistent only if $K$ there means $\sum_j\Phi_j^2$, the reciprocal. This page writes the consistent forms above in the site convention $K=\sum_j\Phi_j^2$ and records the difference from the printed form. The test is simple: the correct rescaling must satisfy $\hat\Phi_j=\sqrt{w_m}\,\Phi_j$ with the weights of (4.3).

**What this paper actually proves.** Theorems 3.1 and 4.1 are **both quoted from prior literature**. The paper's own contributions are the **combination** (aPC basis plus Christoffel-weighted equilibrium sampling), the explicit samplers, and the numerics. **There is no new convergence theorem and no error analysis of the moment-estimation error**: the paper says explicitly that a "density error" is introduced because the mean and variance are estimated, and that "how to quantify and control such errors will be our future projects".

### Numerical experiments

All results are averaged over 100 independent trials, and two sampling rates are compared throughout: $M=CN$ (linear) and $M=CN\log N$ (log-linear). The input distributions used are a discrete binomial $\mathrm{Bino}(n,p)$ mapped to $[-1,1]$, a discrete Poisson $\mathrm{Pois}(\lambda)$ on $[-1,1]$, uniform $U[a,b]$, exponential $\mathrm{Exp}(\mu)$ on $(0,\infty)$, and normal $N(\mu,\sigma)$.

**§5.1 stability.** The condition number of $\hat A=W^{1/2}A$ (printed as "$W^2A$", **flagged as an extraction artifact**) is reported as a mean with 20% and 80% quantiles against polynomial degree, for four $d=2$ test types:

| Case | $\xi_1$ | $\xi_2$ |
| ---- | ----------------------- | ------------------- |
| 1 | $\mathrm{Bino}(20,1/2)$ | $U[-0.6,0.6]$ |
| 2 | $U[-0.8,0.8]$ | $U[-1,1]$ |
| 3 | $\mathrm{Bino}(20,1/2)$ | $\mathrm{Pois}(10)$ |
| 4 | $U[-0.6,0.6]$ | $N(0.1,1.2)$ |

**Case 4 mixes a bounded with an unbounded marginal, so a different equilibrium measure is used in each dimension** — a demonstration in itself of the per-coordinate sampler design. The rates tested are $M=1.5N$, $2N$, $N\log N$ and $1.5N\log N$; a five-dimensional set of test cases is also used.

**§5.2 accuracy**, measured by a discrete $\ell^2$ error over $L$ reference samples (the printed prefactor is "$L\sum_l$", which must be $\frac1L\sum_l$, **flagged**).

- §5.2.1 analytic test functions: $f_1(\xi)=\exp(\sum_k\xi_k)$; $f_2(\xi)=\sum_k0.3+\sin\bigl(\tfrac{16}{15}(\xi_k-0.7)\bigr)+\sin^2\bigl(\tfrac{16}{15}(\xi_k-0.7)\bigr)$; $f_3(\xi)=\exp\bigl(-\sum_kc_k^2(\xi_k-0.01)^2\bigr)$ with $c_k=\exp(-6k/d)$; $f_4(\xi)=\sin(\sum_k\cdot)$.
- §5.2.2 an **electrical resistor network** with $d=2p$ uncertain resistances driven by $V_0=1$, the quantity of interest being the voltage $V$. Tested with $d=2$ and $\xi_i\sim U[10,100]$, and with $d=4$ and $\xi_1\sim\mathrm{Exp}(0.9)$, $\xi_2\sim\mathrm{Exp}(1.1)$, $\xi_3\sim\mathrm{Exp}(0.8)$, $\xi_4\sim\mathrm{Exp}(1.0)$, with **moments estimated from 1000 samples**.
- §5.2.3 the **stochastic elliptic PDE** $-\nabla\cdot(a(y,\omega)\nabla u)=f$ on $D=[0,1]^2$ with $u|_{\partial D}=0$, deterministic load $f=\cos(y_1)\sin(y_2)$, and the Babuška–Nobile–Tempone log-normal-type random diffusion $\log(a_N(y,\omega)-0.5)=1+\xi_1(\omega)(\sqrt\pi L)^{1/2}+\sum_{i=2}^5\zeta_ig_i(y)\xi_i(\omega)$, with the spatial problem solved by standard finite elements.

The stated qualitative outcome throughout is that Christoffel-weighted least squares gives stable and accurate approximations at both the linear and the log-linear sampling rates.

> [!warning] What is transcribable, and the gap to the theory
> The configurations are transcribed here — distributions, dimensions, sampling rates, trial counts, the number of samples used to estimate moments — but **the condition numbers and errors themselves are not**. More importantly, §5.1 measures condition numbers at **finite degree**, whereas the guarantee holds only asymptotically in $N$ because $\hat\rho$ replaces $\tilde\rho$; and the discrete and unbounded inputs of cases 3 and 4 fall outside the reach of Theorem 3.1 or rest on the conjectured densities. **The experiments show the method works in these cases; they do not show the theorems cover them.**

### Relation to the others

Paper 36 is [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|paper 22]] transplanted into the data-driven, epistemic setting: the same Christoffel weight, the same equilibrium-measure sampling, the same asymptotic-in-$N$ caveat, but with the orthonormal family constructed from Hankel moment data rather than from a known density. Its conjectured Gaussian equilibrium density $C(2-\|\xi\|^2)^{d/2}$ is the same conjecture that appears in [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|papers 45 and 28]]. Its direct successor is paper **44**, which replaces asymptotic equilibrium sampling with exact **induced** sampling and least squares with $\ell_1$, thereby removing the $N\to\infty$ caveat. The stability criterion it quotes as Theorem 4.1 is the same Cohen–Davenport–Leviatan result restated as Theorem 6.1 of paper 45 and Theorem 1 of papers 9 and 14.

## 44: sparse approximation of data-driven chaos by induced sampling

(Notation in this section: $Q$ empirical data points forming $S$, from which $M$ training samples are drawn, with $N=|\Lambda^{TD}_K|$ basis functions and $M\ll N$.)

### The idea

Two difficulties are coupled here: the input distribution $\omega$ is known only through a finite sample set $S$ (epistemic uncertainty), and the training budget may be far smaller than the polynomial dimension ($M\ll N$), forcing a sparse recovery. Paper 36 handled the first with data-driven bases plus equilibrium-measure weighted least squares, but equilibrium sampling has two problems: it is optimal only in the degree limit, and — critically in the data-driven setting — **it can place samples where the empirical data has essentially no mass, producing sample locations at which the model was never queried**.

The substitution the paper makes is direct: **replace the equilibrium measure by an exact induced measure supported on $S$ itself, and replace least squares by preconditioned $\ell_1$.** The induced measure is the Christoffel weight used in reverse — sample **more** where the basis functions are large, then undo it by weighting with $1/\kappa$. That "sampling density $\propto\kappa$, weight $\propto1/\kappa$" pairing is one more instance of the uniform-boundedness principle stated at the top of this page, and in the discrete setting it buys one extra thing: **every sample lands on a genuine data point**.

### Setting and derivation

**The arbitrary polynomial chaos basis from moments.** The naive route is to solve a linear system in the moment matrix for the coefficients $\{\beta_j\}_{j=0}^K$, which by itself shows something interesting: **the true density $\rho$ need only be known through its first $2K$ moments to determine every polynomial up to degree $K$**. But that matrix is ill-conditioned for large $K$. The paper instead uses the Hankel matrix

$$
H=\begin{pmatrix}
\nu_0&\nu_1&\cdots&\nu_K\\
\nu_1&\nu_2&\cdots&\nu_{K+1}\\
\vdots&\vdots&\ddots&\vdots\\
\nu_K&\nu_{K+1}&\cdots&\nu_{2K}
\end{pmatrix},
$$

requires Hamburger determinacy $\det H>0$ (**a sufficient condition given is that the empirical set $\Xi$ contains at least $K+1$ distinct samples**), Cholesky-factors $H=R^\top R$, and reads off the three-term recurrence coefficients by the Mysovskih theorem,

$$
a_j=\frac{r_{j,j+1}}{r_{j,j}}-\frac{r_{j-1,j}}{r_{j-1,j-1}},
\qquad
b_j=\frac{r_{j+1,j+1}}{r_{j,j}},
\qquad r_{0,0}=1,\ r_{0,1}=0,
$$

exactly as in paper 36. This route is stated to be empirically far better conditioned than solving the linear system. **Algorithm 1** is this construction: input the data and weights $\{\xi_j,w_j\}_{j=1}^Q$ and the maximal degree $K$, output $\{a_j,b_j\}$, and optionally evaluate the polynomials through the recurrence.

**Tensorisation and its honest caveat.** With $S=\{z^{(1)},\dots,z^{(Q)}\}\subset\mathbb R^d$, marginalise to $\Xi_i=\{z^{(1)}_i,\dots,z^{(Q)}_i\}$ with uniform weights $w_j=1/Q$, run Algorithm 1 per coordinate, and tensorise. The resulting $\Phi_\lambda$ are orthonormal with respect to the **tensorised marginal** measure, not with respect to $S$:

$$
\int_\Gamma\Phi_\lambda(z)\Phi_\theta(z)\,d\nu(z)=\delta_{\lambda,\theta},
\qquad
d\nu(z)=\bigotimes_{i=1}^d d\nu_i(z_i),
\quad
d\nu_i(z_i):=\frac1Q\sum_{j=1}^Q\delta_{z^{(j)}_i}(z).
$$

The paper states plainly that "our multivariate polynomials $\Phi_\lambda$ constructed in this way do not respect the moments of $S$", incurring an extra penalty "essentially equivalent to tensorizing the set $S$" — **that is, any dependence between the input coordinates is discarded**. This is a real modelling loss and the paper does not hide it.

**The induced measure on the data — the paper's central object.** With $N=|\Lambda^{TD}_K|$, the normalised $\Lambda$-Christoffel quantity for the tensorial measure $\nu$ is

$$
\kappa\bigl(z;\Lambda^{TD}_K\bigr):=\frac1N\sum_{j=1}^N\Phi_j^2(z)
\ \ \Bigl(=\frac{K(z)}{N}\ \text{in this page's convention}\Bigr),
$$

and the **discrete induced measure supported on $S$** is

$$
d\mu(z):=\sum_{j=1}^Q\tilde\kappa_j\,\delta_{z^{(j)}},
\qquad
\tilde\kappa_j:=\frac{\kappa(z^{(j)})}{\sum_{q=1}^Q\kappa(z^{(q)})} .
$$

Note the direction: $\kappa$ is the **reciprocal** of the Christoffel function, large where the basis functions are large, so $\mu$ **up-weights** those regions.

**Preconditioned $\ell_1$.**

$$
\min\|c\|_1
\quad\text{s.t.}\quad
\sqrt W A c=\sqrt W b,
\qquad
A_{ij}=\Phi_j(z^{(i)}),\quad b_j=f(z^{(j)}),\quad W_j=1/\kappa\bigl(z^{(j)}\bigr),
$$

where the $z^{(j)}$ are $M$ i.i.d. samples from $\mu$. **Sampling density $\propto\kappa$, weight $\propto1/\kappa$ — an exact inverse pairing**, so the composite $\sqrt WA$ is unbiased with respect to the uniform measure on $S$ and **its entries are uniformly bounded regardless of $S$ and of the shape of the $\Phi_j$**. That uniform boundedness is the bounded-orthonormal-system hypothesis itself, and the penalty terms in compressed-sensing error bounds are proportional to the largest entry of the matrix — **so this choice of $W$ is not a tuning decision but one dictated by the theory**. Because $S$ is finite, drawing i.i.d. samples from $\mu$ is trivial by inverse transform sampling.

**Algorithm 2.** Input the discrete distribution or data set $\{z^{(j)},w_j\}_{j=1}^Q$ and $f$; output $c^*$ with $f\approx\hat f=\sum_jc^*_j\Phi_j$: (1) build the arbitrary polynomial space $\{\Phi_j\}$ with Algorithm 1 per dimension; (2) generate $M$ i.i.d. samples from the induced measure $\mu$; (3) assemble $b_i=f(z_i)$ and $A_{ij}=\Phi_j(z^{(i)})$; (4) form $W$ from $W_j=1/\kappa(z^{(j)})$; (5) solve the preconditioned $\ell_1$ problem.

**Relationship to the equilibrium measure, stated explicitly.** For fairly general $\omega$ there is a unique $\mu_\infty$ with $d\mu_\infty=\lim_{K\to\infty}\omega\kappa$ — the **weighted pluripotential equilibrium measure** — whose explicit form is known only in some cases; on the hypercube $[-1,1]^d$ it is the product Chebyshev measure $\frac{1}{\pi^d\prod_i\sqrt{1-z_i^2}}$. The paper's reasons for preferring $\mu$ over $\mu_\infty$ are:

1. $\mu$ is supported **only on $S$**, so every sample is one at which data exists;
2. the product Chebyshev measure "may be optimal in the limit, [but] it is not optimal for a finite degree $K$";
3. structurally, equilibrium sampling depends only mildly on $\omega$ and $K$, whereas the induced measure depends explicitly on both — **which is why the paper expects, and observes, better data-driven performance**.

### Theorems

**The paper states no new theorem.** This is an algorithms-and-experiments paper, with no numbered Theorem, Lemma or Proposition environments anywhere in it. Its justification is by citation: to established compressed-sensing theory for bounded orthonormal systems (error terms proportional to the largest entry of $\sqrt WA$, which is what dictates $W$), and to the induced-sampling least-squares theory surveyed in [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|paper 45]], together with the two structural advantages listed above.

> [!warning] Do not attach a sample-complexity exponent to this algorithm
> The paper states **no** sample complexity for its own method. Expressions such as $M\gtrsim s\log^3s\log N$ appear in papers 21, 24 and 29 with their own explicit hypotheses, but **none of them is a conclusion about Algorithm 2**. Any such exponent has no source in this paper.

### Numerical experiments

Three method labels run through all experiments: **induced distribution** (this paper), **CSA** ($\ell_1$ with equilibrium-measure sampling, i.e. the method of [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|paper 24]]), and **MC** (non-preconditioned $\ell_1$ with i.i.d. samples from $\omega$). All results are averaged over 100 independent trials, and accuracy errors are computed as a discrete norm over $E=10{,}000$ i.i.d. samples from $\omega$ (**the printed formula divides by $M$ rather than $E$, flagged as a typo**).

**Figure 1 (illustrative)** is the cleanest single picture in this topic of why the non-asymptotic and asymptotic designs differ. The data set is tensorial, $S=\{(z^{(j)}_1,z^{(k)}_2)\}_{1\le j,k\le24}$ with weights $w_{j,k}=u_jv_k$: the $z^{(j)}_1$ are 24 equispaced points on $[-1,1]$ with $u_j$ the probability mass function of $\mathrm{Binomial}(24,0.5)$, and the $z^{(k)}_2$ are 24 equispaced points on $[-1,1]$ with $v_k$ the first 24 values of the $\mathrm{Poisson}(10)$ mass function (truncated mass about $10^{-5}$, renormalised). The equilibrium (product Chebyshev) samples and the $K=20$ induced samples **visibly concentrate in different regions**, **confirming that $K=20$ is not yet in the asymptotic regime**.

**§4.1 exact sparse recovery.** Here $\omega$ is tensorial with isotropic marginals that are equal mixtures of uniform on $[-1,1]$, truncated $N(0.2,1.5)$ on $[-1,1]$ and truncated lognormal on $[0,1]$; $S$ consists of $Q=10^5$ i.i.d. samples from $\omega$; $c^*$ is randomly $s$-sparse with standard normal entries; success is declared when $\|c-c^*\|_\infty<10^{-3}$.

| Scenario | $d$ | $K$ | $N$ | $s$ | Outcome |
| ------------------------- | --- | --- | --- | --- | ------------------------------ |
| low dimension, high degree | 2 | 20 | 231 | 8 | induced sampling **significantly better** than MC |
| high dimension, low degree | 10 | 3 | 286 | 8 | the two are **similar** |

The explanation for the second row generalises: **a low-degree space makes $\mu$ close to $\omega$, so induced and naive sampling converge** — the same degree dependence reported in Example 8.1 of paper 45.

**§4.2 analytic function approximation.** The test functions are $f_1(z)=\exp\bigl(-\sum_iz_i\bigr)$; $f_2(z)=\sum_i(1-z_{i-1})^2+\sum_{i\ge2}100(z_i-z_{i-1}^2)^2$ (Rosenbrock); $f_3(z)=\sin(\sum_iz_i)$; and $f_4(z)=\bigl(1+\frac{1}{2d}\sum_ic_i(1+z_i)\bigr)^{-d-1}$ with $c_i=\frac{1+i}{4d}$ (two of them are printed with the same label $f_3$, **flagged as a numbering typo**). In the two-dimensional case with $S$ as in Figure 1 and $K=20$, **induced sampling is much better than both MC and CSA for larger $M$, while for relatively small $M$ MC moderately outperforms it**. Repeating with the mixture data of §4.1 and in $d=5$ with $K=7$, induced sampling again beats MC as $M$ grows.

**§4.3 parametric PDE.** A clamped **Kirchhoff plate bending** problem on a bounded polygonal $D\subset\mathbb R^2$: $-M_{IJ,IJ}(u)=f$ in $\Gamma\times D$ with $u=\partial_nu=0$ on $\Gamma\times\partial D$, $M_{IJ}(u)=D(z,x)\bigl((1-\nu)K_{IJ}(u)+\nu K_{LL}(u)\delta_{IJ}\bigr)$, $K_{IJ}(u)=-\partial_{IJ}u$, and rigid flexibility $D(z,x)=\frac{E(z,x)h^3}{12(1-\nu^2)}$. The **random Young's modulus** enters through a truncated Karhunen–Loève expansion of $Y(z,x)=\log(E(z,x)-100)$, namely $Y=1+Z_1(\sqrt\pi L)^{1/2}+\sum_{i=2}^d\zeta_ig_i(x)Z_i$, for a squared-exponential covariance $K(x_1,x_1')=\exp\bigl(-(x_1-x_1')^2/L_c^2\bigr)$ along $x_1$ only, with eigenfunctions $g_i(x)=\sin\bigl(\lfloor i/2\rfloor\pi x_1/L_p\bigr)$ for even $i$ and $\cos\bigl(\lfloor i/2\rfloor\pi x_1/L_p\bigr)$ for odd $i$ (**the signs as printed are ambiguous in the extraction, flagged**).

> [!warning] What the experiments establish, and what they do not
> No recovery rates or errors are transcribed here. More to the point: **this paper has no theorems**, so there is no "gap to the theory" to discuss, only the question of what the experiments support. They support the claim that in the data-driven setting, at finite degree and with $M$ not too small, induced sampling usually beats both equilibrium-measure sampling and naive Monte Carlo. They do **not** support the claim that induced sampling always wins — in §4.2 MC is moderately better at small $M$ — nor any quantitative statement about the required sample count.

### Relation to the others

Paper 44 is the direct successor of paper **36** — the same aPC/Hankel basis construction, the same authors plus Narayan — with two substitutions: exact induced sampling in place of asymptotic equilibrium sampling, and $\ell_1$ in place of least squares. Its comparison method "CSA" **is** [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|paper 24]], and the paper spells out the relationship itself: paper 24 "utilize precisely the same preconditioning matrix $W$… but choose to sample not from $\mu$, but instead from a $K$-asymptotic version of $\mu$, called the pluripotential equilibrium measure". The least-squares analogue of that asymptotic choice is paper 22, and paper 45 is cited as the review for induced-measure sampling. The bounded-orthonormal-system justification for $W$ uses the same restricted-isometry machinery as papers **21** and 24.

> [!note] Coverage status
> All six papers — 10, 21, 29, 32, 36 and 44 — were read in full text (10, 21, 32, 36 and 44 from arXiv, 29 from the CiCP record), so all six receive an intuition, a full derivation, theorem statements with their hypotheses, and concrete experimental configurations. Three limitations remain and should be stated:
>
> 1. **Experimental numbers.** All six publish their success-rate and error results as curves. This page transcribes the configurations (dimensions, degrees, sparsity levels, sample counts, trial counts, success thresholds, test functions and application problems) and the authors' qualitative conclusions, but **no point-by-point success rates or errors**.
> 2. **Typographic ambiguity.** Two constants in the Legendre case of paper 32 (the leading factor of eq. (14) and the exponent structure of $W^0_{n,n}$), two normalisation constants in paper 29 (the $W$ of Problem 2 and the optimal $\alpha$), and the $\delta$ of Theorem 4.1 together with the use of $K$ in (4.5)/(4.7) of paper 36, are ambiguous or self-contradictory in the full texts used here. Each is flagged in place with a self-consistent form supplied; **check the typeset originals before quoting them**.
> 3. **Limitations internal to the papers.** Paper 21 prints its main probability inequality in the reversed direction (corrected and flagged here); paper 32 states that $\beta_L\le\mu_L$ was not proved and contains **no sample-complexity theorem at all**; both theorems in paper 36 are quoted from others, its unbounded-domain equilibrium measures are **conjectures**, and it does not analyse the moment-estimation error; paper 44 contains **no new theorem**.

## How the six relate

| No. | Composition of the measurements | Sampling design | Origin of the basis |
| --- | --------------------------------------- | ------------------------------- | --------------------- |
| 10 | values, deterministic points | Weil-sum point set | known distribution |
| 21 | values, subset of a Gauss grid | implicitly Christoffel weighted | known distribution |
| 29 | values plus gradients, trigonometric | preconditioning required | known distribution |
| 32 | values plus gradients, polynomial chaos | preconditioning required | known distribution |
| 36 | values, weighted least squares | set by the empirical measure | **empirical measure** |
| 44 | values, sparse approximation | induced sampling | **empirical measure** |

The same six seen through their theoretical route and what it costs:

| No. | Route | Main guarantee | Price of that guarantee |
| --- | ---------------------------- | ---------------------------------------- | ------------------------------------ |
| 10 | coherence plus Weil sums | deterministic exact recovery | $m\gtrsim4^{\,q}d^2s^2$, quadratic in sparsity |
| 21 | restricted isometry for a bounded orthonormal system | $M\gtrsim L(\mathbf n)s\log^3s\log N$ | $L(\mathbf n)\sim n^{2d/3}$ for unbounded parameters |
| 29 | coherence (Thm 2.1) and restricted isometry (Thm 3.1) | deterministic $\mu(\tilde\Phi)\le\lambda\mu(\Phi)$ with $\lambda<1$ | admissible directions needed; $R_0$ may be large |
| 32 | coherence | mean isotropy and a degree-independent coherence bound | **no sample-complexity theorem**; $\beta_L\le\mu_L$ unproved |
| 36 | quoted least-squares stability | $M\gtrsim N\log N$ | **asymptotic in $N$ only**; unbounded densities are conjectures |
| 44 | quoted bounded-orthonormal-system theory | uniformly bounded entries of $\sqrt WA$ | **no new theorem**; tensorisation discards coordinate dependence |

One judgement runs through all six: **sample complexity in sparse recovery is decided by the uniformity of the row norms, so any operation that changes the row structure — adding gradient rows, changing basis, changing the sampling density — must come with a matching preconditioner.** That unifies this page with the [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|optimal sampling page]]: both handle the same quantity, one inside a least-squares framework and one inside an ℓ1 framework.

A second judgement concerns the trade between deterministic and probabilistic guarantees. Paper 10 and Theorem 2.1 of paper 29 are the only two **deterministic** results in this topic, and they pay for it respectively with a sample count quadratic in the sparsity and with an extra admissibility hypothesis on the directions. **In this field, removing the phrase "with high probability" is not free.**

## Sources for this page

- Z. Xu and T. Zhou, [_On sparse interpolation and the design of deterministic interpolation points_](https://doi.org/10.1137/13094596X), SIAM J. Sci. Comput. 36(4) (2014), pp. A1752-A1769.
- L. Guo, A. Narayan, T. Zhou, and Y. Chen, [_Stochastic collocation methods via ℓ1 minimization using randomized quadratures_](https://doi.org/10.1137/16M1059680), SIAM J. Sci. Comput. 39(1) (2017), pp. A333-A359.
- Z. Xu and T. Zhou, [_A gradient-enhanced ℓ1 approach for the recovery of sparse trigonometric polynomials_](https://doi.org/10.4208/cicp.OA-2018-0006), Commun. Comput. Phys. 24 (2018), pp. 286-308.
- L. Guo, A. Narayan, and T. Zhou, [_A gradient enhanced ℓ1-minimization for sparse approximation of polynomial chaos expansions_](https://doi.org/10.1016/j.jcp.2018.04.026), J. Comput. Phys. 367 (2018), pp. 49-64.
- L. Guo, Y. Liu, and T. Zhou, [_Data-driven polynomial chaos expansions: a weighted least-square approximation_](https://doi.org/10.1016/j.jcp.2018.12.020), J. Comput. Phys. 381 (2019), pp. 129-145. (Another source cites this article as pp. 110-128 of the same volume; the two disagree and 129-145 is used here.)
- L. Guo, A. Narayan, Y. Liu, and T. Zhou, [_Sparse approximation of data-driven polynomial chaos expansions: an induced sampling approach_](https://doi.org/10.4208/cmr.2020-0010), Commun. Math. Res. 36 (2020), pp. 128-153.
