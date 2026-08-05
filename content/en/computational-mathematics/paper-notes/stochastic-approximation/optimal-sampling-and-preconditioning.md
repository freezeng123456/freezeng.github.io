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

## The one idea running through the whole page

Every paper here rests on a single sentence: **the density you sample from need not be the density that defines orthogonality.**

Spelled out. Let $V$ be an $N$-dimensional subspace of $L^2_w$ with $w$-orthonormal basis $\{v_n\}$. The least-squares normal-equation matrix is an average of rank-one matrices,

$$
G=\frac1M\sum_{m=1}^{M}v_mv_m^{T},
\qquad
v_m=\frac{1}{\sqrt M}\bigl(v_1(x_m),\dots,v_N(x_m)\bigr)^{T} .
$$

However you sample, $\mathbb E G=I$. What actually decides how many samples are needed before $G$ is close to $I$ is whether those rank-one terms are **comparable in size**. If $\sum_n v_n^2(x)$ is orders of magnitude larger at some points than elsewhere, most samples contribute almost nothing to $G$ while a handful dominate it, and $M$ must be large before the average settles. That maximum,

$$
\bigl\|K\bigr\|_\infty=\sup_{x\in D}\sum_{n=1}^{N}v_n^2(x),
$$

is the **stability factor**. When one samples from the orthogonality density $w$ it frequently blows up with the polynomial degree, and that is exactly why standard Monte Carlo needs superlinear sample counts.

**There are two mutually inverse ways to level the rows**: put more samples where the row norms are large (change the sampling density), and shrink the large rows by their own size (weight them). Applied together in exactly reciprocal amounts, they leave every $M\to\infty$ limit untouched — it is just an importance-sampling change of variables — and alter only the **preasymptotic** behaviour, which is where the whole argument lives. The yardstick this needs is the quantity that measures **how concentrated the function space is at each point**: the (inverse) Christoffel function.

The five papers differ only in which side of the pair they apply it to. Papers 22, 24 and 28 sample from the equilibrium measure and weight by the Christoffel function, which is optimal in the limit of large degree. Paper 45 takes the inverse Christoffel function to be the sampling density itself, which is optimal at every finite $N$. Paper 28 goes further and replaces randomness altogether by greedy selection.

### Notation used on this page

Throughout,

$$
K(z)=K_\Lambda(z)=\sum_{\alpha\in\Lambda}\varphi_\alpha^2(z),
\qquad
\lambda_\Lambda(z)=\frac{1}{K_\Lambda(z)},
\qquad
N\lambda_\Lambda(z)=\frac{N}{K_\Lambda(z)},
$$

so $K$ is the **sum of squared orthonormal basis functions** and the (normalised) **Christoffel function** is $N/K(z)$; weights therefore always appear as $1/K$. The quantity $K$ is basis-independent: writing $\varphi$ for the vector of basis functions, $K(z)=\varphi^T\varphi$, so any orthogonal change of basis $\psi\leftarrow U\varphi$ leaves it unchanged. Equivalently, $\sum_n v_n^2$ is a property of the subspace $V$ and not of the chosen basis.

The sources do not agree on notation, and this is the commonest source of confusion when reading the family:

| Paper | Symbol            | Meaning                                                                          |
| ----- | ----------------- | -------------------------------------------------------------------------------- |
| 22    | $K$, $K_k$        | $\sum_\alpha\varphi_\alpha^2$, the reciprocal of the Christoffel function; $N/K_k$ is called "the (normalised) Christoffel function" |
| 24    | $\lambda_\Lambda$ | $1/\sum_{i\in\Lambda}\varphi_i^2$, the Christoffel function itself; the preconditioner uses $N\lambda_\Lambda$ |
| 28    | $K_\Lambda$       | $\sum_\alpha\psi_\alpha^2$, the reciprocal; the weight applied is $1/\sqrt{K_\Lambda}$ |
| 45    | $q^2$             | $\frac1N\sum_n v_n^2$, reciprocal and normalised; used directly as the sampling bias |

## 22: the stability factor is something one can design away

### The idea

Sampling from $w$ and doing least squares in $L^2_w$ looks self-evident, because "sampling" and "orthogonality" then use the same $w$. But the two play different roles: $w$ decides what is orthogonal to what, while the sampling density decides where the information comes from. Tying them together throws away the one knob that is free to turn.

**The insight is that the stability factor is an artefact of that mismatch.** Sampling from the pluripotential equilibrium measure and reweighting with the Christoffel function drives the effective stability factor down to its optimal value $N$. Put in a form closer to implementation: **weighting by the Christoffel function is the same as normalising the rows of the design matrix** — the paper's own equivalent reading, which translates "optimal sampling" into a purely linear-algebraic operation and explains why uniform row norms are the decisive property.

What goes wrong without it: for many weights of interest $\|K\|_\infty/N$ blows up with the polynomial degree, forcing superlinear sample counts, often far worse than quadratic.

### Setting

Let $K_k$ be the $L^2_w$ reproducing kernel diagonal of the polynomial space $\mathbb P_k$, with $N=\dim\mathbb P_k$. The quantity $N/K_k(z)$ is the (normalised) **Christoffel function**, which gives the method its name.

CLS in fact approximates in $L^2$ with the modified weight

$$
\tilde w(z)\triangleq\frac{N}{K(z)}\,v(z),
$$

where $v$ is the Lebesgue density of the equilibrium measure. The whole theory therefore measures how far $\tilde w$ is from $w$, using two **discrepancy objects**. The function-independent one is the $\tilde w$-Gramian of the $w$-orthonormal basis,

$$
(R)_{m,n}=\int_D\varphi_m(z)\varphi_n(z)\,\tilde w(z)\,\mathrm dz;
$$

the function-dependent one is $d(f)=\|\tilde\Pi f-\Pi f\|_w$, with $\Pi$ and $\tilde\Pi$ the $L^2_w$- and $L^2_{\tilde w}$-orthogonal projectors onto $P$. It vanishes for $f\in P$, and for every $f$ when $\tilde w=w$.

### Derivation

The first step is to write down the standard Monte Carlo criterion (Theorem 4.1): the sample requirement is that $\frac{S}{N\log S}$ exceed a constant times $\frac{\|K(z)\|_\infty}{N}$. **Everything rests on that factor.**

The second step is an asymptotic identity (the paper's eq. (2)): if $D$ is compact with nonvanishing interior and positive $d$-dimensional Lebesgue measure, and $w$ is continuous on the interior and admits an orthogonal polynomial family, then

$$
\lim_{k\to\infty}\frac{N}{K_k(z)}=\frac{w(z)}{v(z)}\quad\text{a.e. in }D,
$$

with $v$ the Lebesgue density of the pluripotential equilibrium measure of $D$. In $d=1$ on $D=[-1,1]$, $v$ is the arcsine (Chebyshev) density.

The third step turns that identity into a construction. Define the non-polynomial functions

$$
\psi_n(z)=\sqrt{\frac{N}{K_k(z)}}\,\varphi_n(z),
$$

a basis for $\frac{1}{\sqrt{K_k}}\mathbb P_k$, orthogonal in $L^2$ with the modified weight $\frac{wK_k}{N}$. By the limit above they are **approximately orthonormal with respect to $v$**, and their approximate reproducing kernel diagonal is

$$
\tilde K_k=\sum_n\psi_n^2=\frac{N K_k}{K_k}=N,
$$

**attaining exactly the smallest value the supremum can take**. So sampling from $v$ and doing least squares with these functions hits the best possible sample-count criterion. In the paper's own words, that is the CLS algorithm in a nutshell.

The fourth step writes the algorithm out. On a compact domain $D$ (the paper's Algorithm 2):

1. draw $S$ i.i.d. samples $\{z_s\}$ from the equilibrium measure $\mu_D$;
2. assemble the data vector $u$ with $(u)_s=u(z_s)$;
3. compute least-squares weights $K$ with $(K)_{s,s}=N/K(z_s)$;
4. assemble the $S\times N$ Vandermonde-like matrix $V$ with $(V)_{s,n}=\varphi_n(z_s)$;
5. solve $c=\arg\min_{g}\bigl\|\sqrt{K}Vg-\sqrt{K}u\bigr\|$.

Algorithm 3 is the unbounded-domain variant. One detail is worth keeping: for bounded $D$ the sampling density $v(z)=v_D(z)=\frac{\mathrm d\mu_D}{\mathrm dz}$ is **independent of the orthogonality density $w$** — on an interval CLS prescribes Chebyshev sampling whatever $w$ is, which the paper connects to the folklore that the Chebyshev measure on an interval is universal. For $D=[-1,1]^d$, $v_D$ is a product of univariate arcsine measures; the unit ball and other special domains have explicit formulas; on an unbounded conic $D$, $v$ is a scaled version of the $\sqrt w$-weighted pluripotential equilibrium measure.

The fifth step is the accounting. Since the approximation actually happens in $L^2_{\tilde w}$ rather than $L^2_w$, traces of $R$ and $d(f)$ must survive in the error bound. That is where the last term of the main theorem comes from.

### Theorems

**Theorem 4.1 (Cohen-Davenport-Leviatan, restated; the standard Monte Carlo benchmark).** Let $P$ be any $N$-dimensional subspace of $L^2_w$ and $\{z_s\}_{s=1}^S$ i.i.d. from $w$. If

$$
\frac{S}{N\log S}\ \ge\ \Bigl[\frac{1+r}{c_\delta}\Bigr]\frac{\|K(z)\|_\infty}{N},
\qquad c_\delta\triangleq\delta+(1-\delta)\log(1-\delta),\ \delta\in(0,1),\ r>0,
$$

then the discrete Gramian satisfies $\Pr\bigl[|\!|\!|G-I|\!|\!|>\delta\bigr]\le\frac{2}{S^r}$, together with an accuracy estimate for $f$ with $\|f\|_\infty\le L$.

**Theorem 4.2 (Berman, Bloom-Levenberg and others, restated; Christoffel asymptotics).** Let $D\subset\mathbb R^d$ be potential-theoretically admissible with smooth $\mathrm dV(z)=q(z)\mathrm dz$ and $\rho$ bounded continuous such that $\rho\,\mathrm dV$ defines an orthonormal family in $L^2_{q\rho}(D)$, and let $K^{(k)}_k$ be the $L^2_{q\rho^{2k}}$ reproducing kernel diagonal of $\mathbb P_k$. Then

$$
\lim_{k\to\infty}\frac{1}{N}\rho^{2k}(z)K^{(k)}_k(z)\,\mathrm dV(z)=\mathrm d\mu_{D,Q}(z)\quad\text{weakly}.
$$

**Corollary 4.1.** (1) Bounded: for continuous $w$ admitting an orthogonal basis on a compact connected $D$, with $\rho\equiv1$ and $\mathrm dV=w\,\mathrm dz$, $\lim_{k\to\infty}\frac1NK_k(z)=\frac{\mathrm d\mu_D}{\mathrm dV}=\frac{v_D(z)}{w(z)}$. (2) Unbounded convex cone $D$ with $w=\exp(-2Q)$, $\rho=\sqrt w$, $\mathrm dV=\mathrm dz$: $\lim_{k\to\infty}\frac1N\rho^{2k}(z)K^{(k)}_k(z)=v_{D,Q}(z)$.

**Theorem 4.3 (optimal measures, from Bos et al.).** If $\mu_k$ is an optimal measure for $\mathbb P_k$ on $D$ with weight $\rho=\exp(-Q)$, then (i) $\kappa_{k,\rho}=N$ $\mu_{D,Q}$-a.e., and (ii) $\mu_k\to\mu_{D,Q}$ weakly as $k\to\infty$. With Theorem 4.1 this gives **asymptotically simple log-linear scaling $S\log S\gtrsim N$, the best possible sample-count criterion**.

> [!warning] Three limitations of Theorem 4.3, listed by the paper itself
> First, the result gives **no optimality at fixed $k$**, so equilibrium sampling may be quite suboptimal when $k$ is small. Second, it is unclear how large $k$ must be relative to $d$ for the asymptotics to bite. Third, large $k$ is computationally infeasible in high dimension, since $\dim\mathbb P_k\sim k^d$. Those three gaps are exactly what paper 45's induced sampling removes.

**Theorem 5.1 (CLS stability).** For compact $D$ and admissible $w$ with index set $\Lambda$, if

$$
\frac{S}{N\log S}\ \ge\ \Bigl[\frac{1+r}{c_\delta}\Bigr]\frac{1}{\lambda_{\min}(R)},
$$

then the discrete Gramian $G$ of the CLS procedure is close to $R$ with high probability. **The comparison with Theorem 4.1 is the core of the whole thread: the factor $\|K\|_\infty/N$ has been replaced by $1/\lambda_{\min}(R)$.**

**Theorem 5.2 (CLS accuracy, the paper's headline theorem).** Let $D$ be compact, and let CLS mean i.i.d. sampling from the equilibrium measure $v$ with weights $N/K(z)$. For any $r>0$ the procedure is stable with high probability provided

$$
\frac{S}{N\log S}\ \ge\ C\,\frac{1+r}{\lambda_{\min}(R)},
$$

with $C$ an absolute constant. Moreover, for $|f|\le L$, with $T_L(x)=\mathrm{sgn}(x)\min\{|x|,L\}$ and $\tilde\Pi_Sf$ the $S$-sample CLS estimator,

$$
\mathbb E\Bigl[\bigl\|f-T_L(\tilde\Pi_Sf)\bigr\|_w^2\Bigr]
\le\|f-\Pi f\|_w^2
+\frac{\varepsilon(S)}{\lambda_{\min}(R)}\|f-\Pi f\|_{\tilde w}^2
+\frac{8L^2}{S^r}
+4\kappa^2(R)\,d^2(f),
$$

with $\varepsilon(S)=\frac{2-2\log2}{(1+r)\log S}\to0$ and $\kappa(R)=\lambda_{\max}(R)/\lambda_{\min}(R)$.

Theorems 5.3 and 5.4 are the unbounded-domain analogues. The paper says the analysis presents no great difficulty but that **implementation does**, because no explicit formula for the weighted pluripotential equilibrium measure is known for the weights of interest. It therefore **conjectures** forms for them and reports simulations supporting the conjectures (Section 6, Table 2). Those forms are conjectures, not results.

> [!warning] The paper's own reservations about its theory
> This deserves recording faithfully, because it is easily lost in restatement. The last term $4\kappa^2(R)d^2(f)$ **does not vanish as $S\to\infty$**, so taken at face value this theory is **weaker** than the established standard Monte Carlo theory. The paper says plainly that it cannot yet prove the two comparable, and offers two weaker forms of support instead: an argument that its bound is of the same magnitude as the Monte Carlo bounds (Section 5.1.1), and numerical evidence that weighted least squares is **frequently**, not always, superior (Section 6). It also notes that in one dimension $\lambda_{\min}(R)$ and $\kappa(R)$ are numerically well behaved for general weights (Figure 3).
>
> In other words the case for this route rests here more on experiment than on theorem; it is paper 45's induced-sampling criterion below that turns the linear regime into a clean theorem.

### Numerical experiments

What has been checked here is the **design of the experiments and the qualitative outcomes the paper states**. The numbers in the figures — sample counts, condition numbers, error magnitudes — have not been verified item by item, so none are quoted below.

| Experiment | Setup                                                                                          | Outcome as stated                                              |
| ---------- | ---------------------------------------------------------------------------------------------- | -------------------------------------------------------------- |
| Figure 1   | univariate symmetric Jacobi families ($\alpha=\beta$), $N=\dim\mathbb P_k=k+1$; plotted against $\beta$ and against $k$ | the standard Monte Carlo stability factor blows up with degree, which is the motivation |
| §6.1.1     | matrix stability on bounded domains: uniform/Legendre and other Jacobi cases                   | CLS is more stable than standard Monte Carlo                   |
| §6.1.2     | matrix stability on unbounded domains (Gaussian)                                                | supports the conjectured forms of the weighted equilibrium measure (Table 2) |
| §6.1.3     | non-total-degree $\ell_p$ polynomial spaces                                                     | CLS still performs well although the theory is framed for total-degree spaces |
| §6.2.1     | an algebraic function                                                                           | accuracy of CLS against Monte Carlo                            |
| §6.2.2     | a heterogeneous diffusion equation in one spatial dimension                                     | likewise                                                       |
| §6.2.3     | an electrical resistor network                                                                  | likewise                                                       |

The overall finding, in the abstract's wording, is that CLS is superior to standard Monte Carlo in many situations of interest; elsewhere the paper adds "many (but not all)".

**What the experiments establish.** First, Figure 1 confirms that the blow-up of the stability factor is real behaviour of the common weight families, not a worry on paper. Second, §6.1.3 shows the method reaches further than its theorem, since the theory is written for total-degree spaces and the $\ell_p$ spaces work too. Third, §6.1.2 is the **only** support the unbounded-domain conjectures have.

**What they do not establish.** They do not supply the missing piece of Theorem 5.2: the term $4\kappa^2(R)d^2(f)$ that does not vanish with $S$. Being "frequently superior" numerically is not the same as being comparable in theory, and the paper does not claim otherwise. The unbounded-domain experiments are moreover run on conjectured sampling densities, so what they verify is that the method works if the conjectures are right, not the conjectures themselves.

### Relation to the others

This paper is the theoretical centre of gravity of the sampling-design programme. It takes the potential-theoretic vocabulary catalogued in work such as [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 14]] — equilibrium measures, weighted equilibrium measures, Fekete arrays — and turns it into an algorithm, replacing the $\|K\|_\infty/N$ of the Cohen-Davenport-Leviatan bound with $1/\lambda_{\min}(R)$. Its weight $N/K(z)$ is the same object that appears implicitly as the Gauss weight in the randomised-quadrature papers. Using the **inverse** Christoffel function directly as a sampling density, rather than as a weight on equilibrium samples, is paper 45's induced-sampling refinement, which removes the asymptotic-in-$k$ caveat left open by Theorem 4.3. The deterministic counterpart, with the same Christoffel weight but pivoted QR in place of Monte Carlo, is paper 28.

## 45: writing the thread up as a survey

### The idea

The survey has one thesis: **the most "obvious" sampling density is frequently a bad choice.** Its Example 5.2 exists to make the point — a $d=2$ Gaussian problem where even $M=10N$ samples drawn from $w$ produce an estimator orders of magnitude worse than the best approximation.

The remedy is the one described above: bias the sampling to $\rho=q^2w$ and divide the bias back out with the weight $1/q^2$. **The decisive fact is that this pair changes no $M\to\infty$ limit** — the Gramian still tends to $I$, the right-hand side still tends to $\hat f_n$, the objective still tends to $\|g-f\|^2_{L^2_w}$ — so $q$ is a completely free design variable that may be chosen purely to make the preasymptotic Gramian concentrate quickly. And that choice problem has an exact solution.

### Setting

$D\subset\mathbb R^d$, $w:D\to[0,\infty)$ a **probability density**, and $L^2_w$ with $\langle u,v\rangle=\int_Duvw\,\mathrm dx$. $V\subset L^2_w$ with $\dim V=N$ and $L^2_w$-orthonormal basis $v_1,\dots,v_N$. The best approximation is

$$
f_N(x)=\sum_{n=1}^N\hat f_nv_n(x),\quad \hat f_n=\langle f,v_n\rangle,
\qquad f_N=\arg\min_{v\in V}\|f-v\|.
$$

Since $g_N-f_N\in V$ while $f-f_N\perp V$, Pythagoras gives a clean quality metric,

$$
\eta_N:=\frac{\|g_N-f_N\|}{\|f-f_N\|},
\qquad \|f-g_N\|^2=(1+\eta_N^2)\|f-f_N\|^2 .
$$

$\eta_N\approx1$ means the computed approximation is as good as the best one. **This single number is what every experiment in the survey reports.**

Multi-index spaces: $V(\Lambda):=\mathrm{span}\{x^\lambda:\lambda\in\Lambda\}$, the $\ell^p$ balls $\Lambda_p(k)$, the hyperbolic cross $\Lambda^{HC}(k):=\{\lambda:\|\log(\lambda+1)\|_1\le\log(k+1)\}$, and the specialisations $\Lambda^{TD}(k)=\Lambda_1(k)$ (total degree), $\Lambda^{ED}(k)=\Lambda_2(k)$ (Euclidean degree) and $\Lambda^{TP}(k)=\Lambda_\infty(k)$ (tensor product). Their dimensions are ordered:

| Index set        | Symbol            | Dimension asymptotics                                        |
| ---------------- | ----------------- | ------------------------------------------------------------ |
| hyperbolic cross | $\Lambda^{HC}(k)$ | $\sim(k+1)\log(k+1)^{d-1}$                                   |
| total degree     | $\Lambda^{TD}(k)$ | $(k+1)^d/\Gamma(d+1)$                                        |
| Euclidean degree | $\Lambda^{ED}(k)$ | $\bigl[\tfrac{\sqrt\pi}{2}(k+1)\bigr]^d/\Gamma(\tfrac d2+1)$ |
| tensor product   | $\Lambda^{TP}(k)$ | $(k+1)^d$                                                    |

Design matrix: unbiased, $(A)_{m,n}=\frac{1}{\sqrt M}v_n(x_m)$ and $(f)_m=\frac{1}{\sqrt M}f(x_m)$, with normal equations $Gc=g$, $G=A^TA$, $g=A^Tf$. The paper advises solving the least-squares system by QR rather than forming the normal equations, which are used only for analysis.

> [!note] A printing inconsistency
> Example 2.2 prints the Gaussian weight as $w(x)=(2\pi)^{-d}\exp(-\|x\|_2^2)$, which is inconsistent with the $w=\exp(-\|x\|_2^2)/\pi^{d/2}$ actually used later in Example 8.1. The numerics use the latter, and so does this page; the discrepancy is flagged here.

### Derivation

**Biased sampling.** Take any positive $q\in L^2_w$ with $\|q\|=1$, so that $\rho(x):=q^2(x)w(x)$ is again a probability density on $D$. Draw $x_1,\dots,x_M$ i.i.d. from $\rho$ and set

$$
(A)_{m,n}=\frac{1}{\sqrt{Mq^2(x_m)}}v_n(x_m),
\qquad
(f)_m=\frac{1}{\sqrt{Mq^2(x_m)}}f(x_m),
$$

equivalently $g_N=\arg\min_{g\in V}\frac1M\sum_m\frac{(g(x_m)-f(x_m))^2}{q(x_m)^2}$. Taking $q\equiv1$ recovers the unbiased case. The invariance stressed above lives here: sampling from $\rho=q^2w$ **and** weighting by $1/q^2$ leaves $(G)_{m,n}\to\delta_{m,n}$, $(g)_n\to\hat f_n$ and the objective $\to\|g-f\|^2_{L^2_w}$ all unchanged.

**Stability is proved by a matrix Chernoff argument**, worth reproducing in full because the constant comes out of it. Write $G=\sum_mV_m$ with $V_m=v_mv_m^T$ rank-one positive semidefinite and $v_m^T=\frac{1}{\sqrt Mq(x_m)}(v_1(x_m),\dots,v_N(x_m))^T$. Then $\mathbb EV_m=\frac1MI$, so $\mathbb EG=I$ and $\tau_{\min}=\tau_{\max}=1$. The summand norm bound is

$$
\lambda_{\max}(V_m)=\|v_m\|^2
=\frac1M\sum_{n=1}^N\Bigl(\frac{v_n(x_m)}{q(x_m)}\Bigr)^2
\le\frac1M\sup_{x\in D}\sum_{n=1}^N\Bigl(\frac{v_n(x)}{q(x)}\Bigr)^2=:Q .
$$

The two Chernoff bounds are $\Pr[\mathcal E_{\min}]\le N(2/e)^{\tau_{\min}/2Q}$ and $\Pr[\mathcal E_{\max}]\le N(8e/27)^{\tau_{\max}/2Q}$, so together

$$
\Pr[\mathcal E]\le 2N\exp\Bigl(-\frac{1}{2Q}\log\frac{27}{8e}\Bigr).
$$

Requiring this to be at most $2M^{-r}$ amounts to $\frac{\log(27/8e)}{2Q}\ge(r+1)\log M$, which rearranges into the criterion below and gives $\Pr[\mathcal E]\le2NM^{-r}M^{-1}\le2M^{-r}$. **That is where the constant $C=2/\log(27/8e)$ comes from**; it is not fitted.

**The optimality argument.** The only $w$-, $D$- and $d$-dependent quantity in the criterion is $\sup_x\sum_n(v_n/q)^2$, and it can never be smaller than $N$:

$$
\sup_{x\in D}\sum_{n=1}^N\Bigl(\frac{v_n(x)}{q(x)}\Bigr)^2
\ \ge\ \int_D\sum_{n=1}^N\Bigl(\frac{v_n(x)}{q(x)}\Bigr)^2\rho(x)\,\mathrm dx
\ =\ \sum_{n=1}^N\int_Dv_n^2(x)w(x)\,\mathrm dx\ =\ N .
$$

The first step is only "a supremum is at least an average"; in the second, substituting $\rho=q^2w$ makes $q^2$ cancel exactly — **that cancellation is the technical core of the entire thread**. And the choice

$$
q^2(x)=\frac1N\sum_{n=1}^Nv_n^2(x)
$$

**attains the lower bound exactly** (attributed to Cohen-Migliorati). It is the **induced distribution**

$$
\rho(x)=\frac1N\sum_{n=1}^Nv_n^2(x)\,w(x),
$$

that is $q^2=K_\Lambda/N=1/(N\lambda_\Lambda)$: the inverse Christoffel function itself. It is basis-independent, a property of the subspace $V$. The word "induced" is borrowed from Gautschi and Li's work on a set of orthogonal polynomials induced by a given orthogonal polynomial.

Substituting back gives $M/\log M\ge C(r+1)N$. **All the problem-dependence has thereby been moved out of the sample complexity and into the sampling density**: the complexity retains only $N$, and every difficulty concentrates into the single question of how to sample from $\rho$.

### Theorems

**Theorem 6.1 (a specialisation of Theorem 1 of Cohen-Davenport-Leviatan; stability of the normal equations).** With $A$ and $G$ as above, if for some $r>0$

$$
\frac{M}{\log M}\ \ge\ C(r+1)\,
\sup_{x\in D}\sum_{n=1}^{N}\Bigl(\frac{v_n(x)}{q(x)}\Bigr)^{2},
\qquad
C=\frac{2}{\log(27/8e)}\approx9.24,
$$

then $\|G-I\|_2\le\frac12$ with probability at least $1-2M^{-r}$.

**The sample-complexity criterion under induced sampling.** Substituting the optimal $q$,

$$
\frac{M}{\log M}\ \ge\ C(r+1)\,N,
$$

that is $M\sim N\log N$, optimal in $N$ up to the logarithm. **This criterion depends on $N=\dim V$ alone**: not on the dimension $d$, not on the domain $D$, not on the weight $w$, not even on which $N$-dimensional subspace was chosen. The price is explicit: one must sample from the nonstandard density $\rho$, and $\rho$ **does** depend on $(V,w,D)$.

**Lemma 5.1 (asymptotic consistency).** For any admissible $q$, $g_N\to f_N$ almost surely in $L^2_w$ as $M\to\infty$, by the strong law of large numbers together with $G^{-1}\to I$, $g\to\hat f$ and the continuous mapping theorem. **So the choice of $q$ matters only in the preasymptotic regime**, a point worth repeating.

**Theorem 7.2 (from Cohen-Migliorati; accuracy with truncation).** Assume $f\in L^2_w$ is bounded with $\sup_D|f|=L<\infty$ and set $T_L(y)=\mathrm{sign}(y)\min\{|y|,L\}$. Let the $x_m$ be i.i.d. from the induced density $\rho$ and let $M$ satisfy the criterion above. Then

$$
\mathbb E\,\|f-T_L\circ g_N\|^2\ \le\ \Bigl(1+\frac{4}{C(1+r)\log M}\Bigr)\|f-f_N\|^2\ +\ 8L^2M^{-r} .
$$

In expectation the truncated least-squares approximation therefore commits an error comparable to the best possible one. The paper notes that high-probability statements for the **untruncated** $g_N$ are also available.

**The sample-complexity comparison.** The most informative part of the survey places three routes side by side:

- **i.i.d. sampling from the orthogonality (uniform) measure** (Cohen-Davenport-Leviatan): for any $r>0$, if $M/\log M\ge C_rN^2$ then $\Pr[|\!|\!|\hat A-I|\!|\!|\ge\frac12]\le2M^{-r}$. This is the **quadratic** requirement $M\gtrsim N^2\log N$, and it extends to multidimensional spaces with arbitrary **lower** index sets under the same scaling.
- **Monte Carlo sampling from the Chebyshev measure**: the requirement falls to $M\sim N^{\log3/\log2}$, strictly better than $N^2$ but still superlinear.
- **Deterministic Weil point sets** (from paper 9): $M\ge C(d)N^2$ gives a unique solution and near-best approximation; that quadratic requirement is stronger than Chebyshev Monte Carlo, and the compensating advantage is determinism.

The paper is equally explicit about the gap between theory and practice: practitioners typically take $M\simeq cN$ with $c$ between 2 and 3, that is **linear**, and definitive theory for the linear regime "is not yet definitively available". **That passage is the coordinate system for the whole thread**: the linear regime identified in papers 22 and 45 is exactly the target that paper 28 reaches by entirely different means.

**Asymptotic induced measures (Section 8).** In one dimension with $w$ uniform on $[-1,1]$ and $V$ the polynomials of degree $\le N-1$, $\lim_{N\to\infty}\rho(x)=\rho_\infty(x)=\frac{1}{\pi\sqrt{1-x^2}}$ weakly, the arcsine (Chebyshev) density — which is why Chebyshev grids are the classical answer for polynomial approximation on an interval, and which is the large-$N$ optimal strategy identified in paper 22. For $D=[-1,1]^d$ with $w$ uniform and $V=V(\Lambda^{TD}(k))$, the limit is the tensorised Chebyshev density $\frac{1}{\pi^d\prod_j\sqrt{1-(x^{(j)})^2}}$.

> [!warning] The Gaussian case is only a conjecture
> For the Gaussian weight on $D=\mathbb R^d$, paper 22 **conjectures**
> $$\lim_{k\to\infty}\rho\bigl(x/\sqrt k\bigr)=C\bigl(2-\|x\|_2^2\bigr)^{d/2},$$
> with $C$ a normalising constant and the input scaled by $1/\sqrt k$. **This is a conjecture, not a theorem**, and it is treated as one everywhere it appears on this page. The survey itself is blunt about the multivariate case: "much more is unknown, and in many cases we currently have only conjectures", and the large-$k$ asymptotics of $\rho^{TP}$, $\rho^{ED}$, $\rho^{TD}$ and $\rho^{HC}$ are "largely understudied".

**Sampling from $\rho$ is cheap (Section 8.1).** General multivariate densities are onerous to sample, but the induced density is "an additive mixture of tensor-product densities" and can be sampled with complexity **linear in the dimension $d$**. Implementation is referenced to Narayan's work on computing induced orthogonal polynomial distributions, with software at `https://github.com/akilnarayan/induced-distributions`.

### Numerical experiments

All three experiments come with complete setups, and they are the firmest numerical evidence on this page.

**Example 5.2 — the model failure of standard sampling.**

| Item          | Setup                                                                    |
| ------------- | ------------------------------------------------------------------------ |
| domain, weight | $D=\mathbb R^2$, $w(x)=\exp(-\lVert x\rVert^2)/\pi$                     |
| space         | $\Lambda=\Lambda^{TD}(k)$, $N=\binom{k+2}{k}=(k+1)(k+2)/2$               |
| degrees       | $k=1,\dots,25$                                                           |
| test function | $f(x)=B\bigl(\lVert x/4-(0.2,-0.1)\rVert_2\bigr)$, $B$ a univariate bump function |
| samples       | $M=10N$                                                                  |
| trials        | 100                                                                      |
| result        | sampling from $w$ gives an "extremely inaccurate" $g_N$ with very large $\eta_N$; induced sampling at $k=20$ ($N=231$) gives $\eta_N\sim1$ at moderate $M/N$ |

**Example 8.1 — the effect of dimension.**

| Item          | Setup                                                                            |
| ------------- | -------------------------------------------------------------------------------- |
| domain, weight | $D=\mathbb R^d$, $w=\exp(-\lVert x\rVert_2^2)/\pi^{d/2}$                        |
| test function | $f(x)=\prod_{j=1}^d\exp\bigl([x^{(j)}]^2/j\bigr)$, in product form so that $f_N$ stays computable in high dimension |
| space         | $\Lambda=\Lambda^{HC}(k)$                                                        |
| $(d,k)$       | $(4,20)$, $(8,10)$, $(20,5)$                                                     |
| trials        | 100                                                                              |
| result        | induced sampling consistently outperforms standard sampling, **but the advantage diminishes as $d$ grows**; the qualitative behaviour of the induced-sampling $\eta_N$ is essentially unchanged with $d$, as Theorems 6.1 and 7.2 predict |

The paper is clear about why the advantage shrinks: large $d$ forces small $k$, and a low-degree space makes $\rho$ close to $w$ to begin with.

**Section 9 — a parametric thermal-diffusion equation.**

| Item             | Setup                                                                       |
| ---------------- | --------------------------------------------------------------------------- |
| equation         | $-\nabla\cdot(a(y,x)\nabla u)=S$ on the unit square, $u\rvert_{\partial\Omega}=0$ |
| source           | $S(y,x)=100\chi_F(y)$, $F$ a centred square subdomain of side $0.2$         |
| geometry         | four circular inclusions of radius $r=0.13$, symmetric about the centre     |
| coefficient      | $a(y,x)=1+\sum_{i=1}^4x^{(i)}\chi^i(y)$, so $d=4$                           |
| parameters       | $x^{(i)}\sim U(-0.99,0.2)$ i.i.d., giving a tensor Legendre basis           |
| quantity of interest | $f(x)=u\bigl((0.25,0.375),x\bigr)$                                      |
| discretisation   | $P_1$ finite elements in space; $V=V(\Lambda^{HC}(25))$                     |
| error estimate   | $Q=10{,}000$ reference samples, 100 trials                                  |
| result           | induced-distribution weighted least squares is "perhaps a bit better", but "the difference is not substantial" |

Two further figures: Figure 4 compares $w$, $\rho$ and $\rho_\infty$ at $N-1=19$ for three univariate weights (uniform on $[-1,1]$, exponential on $[0,\infty)$, Gaussian on $\mathbb R$), showing $\rho$ and $\rho_\infty$ visually close while $w$ and $\rho$ can differ substantially, with $\rho$ demanding samples where $w$ demands very few; Figure 5 contours the TD, ED, TP and HC induced densities for the $d=2$ Gaussian at $k=8$, showing that they differ substantially from each other and from $w$.

**What the experiments establish.** Example 5.2 turns "standard sampling can be disastrous" from a warning into a reproducible fact. Example 8.1 confirms the theory's most important property, that induced sampling **does not degrade with dimension**. Figures 4 and 5 show that changing the density is not a tweak but a change of magnitude.

**What they do not establish.** On the real PDE of Section 9 induced sampling is only slightly better. The paper does not dodge this and neither does this page: **the claim that survives is that induced sampling is consistently among the best and is the one with a minimal-sample-count theorem behind it, not that it always wins numerically.** Note also that Theorem 7.2 guarantees near-best accuracy for the **truncated** estimator in **expectation**, whereas the experiments report the distribution of $\eta_N$; neither implies the other directly.

### Relation to the others

Paper 45 is the pedagogical capstone of the thread. It presents paper 22 as the **asymptotic** optimal-sampling result — equilibrium-measure sampling with Christoffel weighting, optimal only as the degree tends to infinity — and induced sampling (Cohen-Migliorati) as the **non-asymptotic exact** version, which is why its criterion carries no asymptotic caveat where paper 22's Theorem 4.3 does. Its $q^2=\frac1N\sum v_n^2$ is precisely $1/(N\lambda_\Lambda)$ in paper 24's notation and $K_\Lambda/N$ in paper 28's: **all three manipulate the same function, and the difference is only whether it is used as a weight on equilibrium samples (22, 24, 28) or as the sampling density itself (45).** Its Chebyshev limit $\rho_\infty$ is exactly the Chebyshev preconditioning of paper 24's CSA-a, and the Gaussian conjecture it quotes has the same functional form as the candidate density of paper 28, §4.2.

## 28: greedy selection instead of random sampling

### The idea

Randomised designs give only with-high-probability guarantees, and classical **approximate Fekete points** — greedy maximisation of $|\det V|$ — are defined only for **compact** $\Gamma$, so they simply **do not exist** for Gaussian/Hermite problems. On top of that there are four objectives one might optimise — the global determinant, the global condition number, and the greedy version of each — with no a priori reason for them to agree.

**Inserting the Christoffel weight into the determinant objective resolves both at once.** The mechanism can be made concrete: in the weighted space $Q$, every row of the Vandermonde matrix has norm exactly 1, since the squared norm of row $j$ is $\sum_\alpha\psi_\alpha^2(y_j)/K_\Lambda(y_j)=1$. For a matrix with unit-norm rows Hadamard's inequality gives $|\det V|\le1$, with equality precisely when the rows are mutually orthogonal, that is when $\kappa(V)=1$. **So "maximal determinant" and "minimal condition number" must coincide at the optimum** — a property of the weighted space with no unweighted counterpart. (This reading of Theorem 3.1 is this page's own; the paper's proof has not been checked line by line here.) Note how exactly it matches paper 22's remark about itself: Christoffel weighting is row normalisation of the design matrix.

As for unbounded domains: $1/\sqrt{K_\Lambda}$ suppresses polynomial growth far out, so the weighted formulation is defined on non-compact $\Gamma$ where unweighted approximate Fekete points are not.

### Setting

$\Gamma=\prod_{i=1}^d\Gamma_i\subset\mathbb R^d$ is tensorial, $\rho(y)=\prod_i\rho_i(y_i)$ a tensor-product probability density with independent components, $\phi^i_n$ the degree-$n$ orthonormal polynomial for $\rho_i$, and the multivariate basis $\psi_\alpha(y)=\prod_j\phi^j_{\alpha_j}(y_j)$ satisfies $\langle\psi_\alpha,\psi_\beta\rangle=\delta_{\alpha,\beta}$. The index sets used are total degree $\Lambda^{TD}_k=\{\alpha:|\alpha|\le k\}$ and hyperbolic cross $\Lambda^{HC}_k=\{\alpha:\prod_j(\alpha_j+1)\le k+1\}$, with $P=\mathrm{span}\{\psi_\alpha:\alpha\in\Lambda\}$ and $N=|\Lambda|$.

The inverse Christoffel function and the weighted space are

$$
K_\Lambda(y)=\sum_{\alpha\in\Lambda}\psi_\alpha^2(y),
\qquad
Q=\mathrm{span}\Bigl\{\tfrac{\psi_\alpha}{\sqrt{K_\Lambda}}\ \Big|\ \alpha\in\Lambda\Bigr\}.
$$

The Vandermonde-like matrices are $(V(A_m,P))_{j,k}=\psi_{\alpha(k)}(y_j)$ and $(V(A_m,Q))_{j,k}=\psi_{\alpha(k)}(y_j)/\sqrt{K_\Lambda(y_j)}$, and the weighted least-squares problem solved is

$$
\hat f=\arg\min_{v\in\mathbb R^N}\bigl\|V(A_M,Q)\,v-Wf\bigr\|_2,
\qquad (W)_{m,m}=1/\sqrt{K_\Lambda(y_m)} .
$$

The determinant modulus of a rectangular matrix is defined by $|\det V|=\sqrt{|\det(VV^T)|}$ for $1\le m\le N$, coinciding with the usual determinant modulus when $m=N$. The unweighted Fekete and condition-number-optimal sets are

$$
A^F_N(P):=\arg\max_{A_N\in\Gamma^N}|\det V(A_N,P)|,
\qquad
A^C_N(P):=\arg\min_{A_N\in\Gamma^N}\kappa\bigl(V(A_N,P)\bigr),
$$

which are different sets in general. Fekete points give the classical Lebesgue-constant bound $\|I_N\|_{C(\Gamma)\to C(\Gamma)}\le N$, and the paper notes that logarithmic growth is what is observed in practice. Replacing $P$ by $Q$ throughout gives the weighted global objectives, with greedy iterations

$$
y^{F*}_{n+1}=\arg\max_{y\in\Gamma}\bigl|\det V\bigl(A^{F*}_n\cup y,\ Q\bigr)\bigr|,
\qquad
y^{C*}_{n+1}=\arg\min_{y\in\Gamma}\kappa\bigl(V(A^{C*}_n\cup y,\ Q)\bigr).
$$

The first is the **CFP method** (Christoffel-weighted approximate Fekete points), written compactly as $y_{n+1}=\arg\max_y\det|WVV^TW|$. The starting point $y^{F*}_1$ is arbitrary in $\mathbb R$ but **does affect the result**, and ties produce multiple "branches".

### Derivation

One dimension has a complete structure, and this is the prettiest part of the paper. Take any probability density $\rho$ on $\Gamma=\mathbb R$ with $\Lambda=\{0,\dots,N-1\}$, write $\phi_N$ for the degree-$N$ orthonormal polynomial, and define the meromorphic function

$$
r_N(y)=\frac{\phi_N(y)}{\phi_{N-1}(y)} .
$$

For any $y\notin\phi_{N-1}^{-1}(0)$ the set $A_N(y)=r_N^{-1}\bigl(r_N(y)\bigr)$, taking the set-valued inverse, is **uniquely determined**, and it carries an $N$-point positive quadrature rule exact to degree $2N-2$ whose weights are precisely the **Christoffel weights**:

$$
\int_\Gamma p(z)\rho(z)\,\mathrm dz=\sum_{z\in A_N(y)}\frac{1}{K_\Lambda(z)}\,p(z),
\qquad \deg p\le2N-2 .
$$

In particular $A_N(y)$ is the $N$-point **Gauss** rule when $y\in\phi_N^{-1}(0)$. **The $1/K_\Lambda(z)$ appearing here is the same object as paper 22's weight $N/K(z)$**, up to the factor $N$, so two apparently different routes meet at this point; it is also the source of the fact that Gauss quadrature weights are inverse-Christoffel weights.

The one-dimensional proof then runs as follows: first show that configurations with $|\det V(A_N,Q)|=1$ exist and are characterised by the level sets of $r_N$, then show that the greedy iteration started anywhere inside such a configuration rebuilds the whole of it.

Three implementation details cannot be skipped:

- **Oversampling by enriching the index set (§4.1).** CFP by itself produces exactly $N=\dim Q$ points, while least squares wants $M>N$. Given a downward-closed $\Lambda$, compute $n=\max\{|\alpha|:\alpha\in\Lambda\}$, set $S:=\Lambda^{TD}_n\setminus\Lambda$ (with $n\leftarrow n+1$ first if $\Lambda=\Lambda^{TD}_n$), order $S$ by total degree with ties broken in reverse lexicographic order, and append the first $\Delta N$ elements to get $\tilde\Lambda$ of size $M=N+\Delta N$, defining $\tilde Q:=Q_{\tilde\Lambda}$.
- **The continuum maximisation becomes a candidate pool (§4.2).** Maximisation over $\Gamma$ is replaced by maximisation over a finite candidate set $\tilde A$.
- **The greedy step is pivoted QR (§4.3).** Form $V(\tilde A,\tilde Q)$ and greedily choose the rows maximising the spanned volume; "this, in turn, is easily performed by a column-pivoted QR decomposition of $V^T$", and the ordered pivots give $A_M$. Note the transpose: pivoting on the **columns** of $V^T$ selects **rows** of the candidate Vandermonde matrix, that is, selects points.
- **The final solve (§4.4).** With $A_M$ fixed, the weighted problem is solved on the **original** index set $\Lambda$ of size $N$, not on the enriched one.

### Theorems

**Theorem 3.1 (determinant optimality $\Leftrightarrow$ condition-number optimality, and greedy attains the global optimum).** Let $\rho:\Gamma\to[0,\infty)$ be a probability density on $\mathbb R^d$ and $\Lambda$ an arbitrary multi-index set of size $N$ defining $Q$. A configuration $A_N$ satisfies $|\det V(A_N,Q)|=1$ **if and only if** $\kappa(V(A_N,Q))=1$, so the solutions of the two global problems attaining optimal objective values coincide. Moreover, if $A_N$ attains either (hence both), then the determinant-greedy iteration started at $y^{F*}_1\in A_N$ has a branch with $A^{F*}_N=A_N$, and the condition-number-greedy iteration started at $y^{C*}_1\in A_N$ has a branch with $A^{C*}_N=A_N$.

> [!warning] The premise of Theorem 3.1
> The theorem **presumes the existence** of a unit-condition-number configuration. The authors say plainly that in $d>1$ this existence is hard to verify, citing known nontrivial multidimensional examples. The multidimensional optimality conclusions are therefore **conditional**.

**Lemma 3.1 (the one-dimensional structure).** The five statements used in the derivation above: for $y\notin\phi_{N-1}^{-1}(0)$ a configuration $A_N=A_N(y)$ attaining the optimum exists; $A_N(y)$ is unique as a function of $y$; $A_N(y)=r_N^{-1}(r_N(y))$; $A_N(y)$ carries an $N$-point positive quadrature rule exact to degree $2N-2$ with Christoffel weights; and $A_N(y)$ is the Gauss rule when $y\in\phi_N^{-1}(0)$. The paper notes that $\rho$ must be regarded as a density on all of $\mathbb R$ even when its support is compact.

**Theorem 3.2 (in one dimension greedy = global = optimal, for almost every start).** With $\rho$ any density on $\Gamma=\mathbb R$ and $\Lambda=\{0,\dots,N-1\}$, fix $y\notin\phi_{N-1}^{-1}(0)$, take $A_N(y)$ from Lemma 3.1, set $A^F_N(Q)=A^C_N(Q)=A_N(y)$ and initialise $y^{F*}_1=y^{C*}_1=y$. Then all four point sets coincide and

$$
\bigl|\det V\bigr|=1=\kappa(V),
$$

so the resulting Vandermonde matrix is **perfectly conditioned**. In one dimension the greedy choice is therefore not an approximation but attains the optimum. More importantly this holds for **any** univariate density, **including densities without compact support** — which is exactly the property approximate Fekete points lack, so it is also a positive answer to the unbounded-domain difficulty.

> [!note] A misprint
> The theorem's text prints the starting condition as $y\in\phi_{N-1}^{-1}(0)$, which contradicts the Lemma 3.1 and Corollary 3.1 it relies on, where $y$ must lie **outside** the zero set. It should read $y\notin\phi_{N-1}^{-1}(0)$. This page states the latter and flags the discrepancy here.

**Corollary 3.1.** For $y\notin\phi_{N-1}^{-1}(0)$ the greedy iteration started at $y^{F*}_1=y$ produces the unique positive $L^2_\rho$ quadrature abscissae with optimal polynomial accuracy; if $y\in\phi_N^{-1}(0)$ it produces the $\rho$-weighted Gauss abscissae.

**What is not proved, as the paper states.** There is **no sample-complexity theorem** here, no bound of the form $M\gtrsim N\log N$. The oversampling procedure of Section 4 is explicitly called "largely ad hoc" and no optimality is claimed for it. In $d>1$ every optimality conclusion is conditional on the existence assumption of Theorem 3.1.

### Numerical experiments

All experiments share one protocol:

| Item          | Value                                                                                     |
| ------------- | ----------------------------------------------------------------------------------------- |
| candidate pool | $\tilde M=10^4$, half i.i.d. draws from $\rho$ and half from a degree-asymptotic density inspired by paper 22 |
| oversampling  | $\Delta N=\lfloor0.05N\rfloor$, that is $M=1.05N$ — **linear, 5% oversampling**            |
| trials        | 50 per configuration; mean condition number reported with 20% and 80% quantiles            |
| baselines     | MC (i.i.d. from $\rho$, unweighted $V(A_M,P)$), Fekete (AFP, unweighted), C-Fekete (CFP, weighted $V(A_M,Q)$) |

The second candidate ensemble is tensor-product Chebyshev for uniform $\rho$ on $[-1,1]^d$, and for Gaussian $\rho$ on $\mathbb R^d$ it is supported on the ball of radius $\sqrt{2n}$ with density

$$
C_d\Bigl(1-\tfrac{1}{2n}\sum_{k=1}^ds_k^2\Bigr)^{d/2},\qquad s\in\mathbb R^d,
$$

with $C_d$ a normalising constant. **This is the same functional form as the asymptotic Gaussian induced measure quoted in paper 45, and that form is only a conjecture**, so the candidate density here should likewise be read as a heuristic inspired by a conjecture. The authors note that weakly admissible meshes would be the natural candidates but that known constructions grow too fast in dimension.

| Experiment | Setup                                                                                           | Result                                                     |
| ---------- | ------------------------------------------------------------------------------------------------- | ----------------------------------------------------------- |
| §5.1.1     | matrix stability, bounded: Legendre on $[-1,1]^d$, $d=2,6,10$, both TD and HC index sets        | CFP gives notably more stable systems than AFP or MC       |
| §5.1.2     | matrix stability, unbounded: Gaussian $\rho\propto\exp(-\lVert y\rVert_2^2)$, tensor Hermite, $d=2,6,10,25$ | CFP is more stable in every case, but the authors state the improvement is **modest in high dimensions** |
| §5.2(a)    | accuracy: $f(y)=\exp\bigl(-\sum_{j=1}^dy_j^2\bigr)$, Legendre approximation, $d=2$              | CFP is comparable to but no worse than AFP, and much better than MC |
| §5.2(b)    | accuracy: the stochastic elliptic equation below, $d=2,6$ (TD) and $d=10,25$ (HC)               | likewise                                                    |

The equation in §5.2(b) is, on $(0,1)\times\mathbb R^d$,

$$
-\frac{\mathrm d}{\mathrm dx}\Bigl[\kappa(x,y)\frac{\mathrm du}{\mathrm dx}(x,y)\Bigr]=f,
\qquad u(0,y)=u(1,y)=0,\quad f=2,
$$

with a Karhunen-Loève-like diffusivity

$$
\kappa(x,y)=1+\sigma\sum_{k=1}^d\frac{1}{k^2\pi^2}\cos(2\pi kx)\,y_k,
$$

quantity of interest $u(0.5,y)$ and $\rho$ uniform on $[-1,1]^d$. (The value of $\sigma$ has not been verified here and is therefore not quoted.)

**What the experiments establish.** First, at only 5% oversampling — $M=1.05N$, extremely close to $N$ — CFP still yields usable condition numbers, which is precisely the regime paper 45 identifies as where practitioners live and theory cannot reach. Second, the unbounded row is itself a result, because unweighted AFP **is not even defined** there. Third, the two accuracy studies show that the weighting does not buy stability at the cost of accuracy.

**What they do not establish.** The reported quantity is the **condition number**, not the Lebesgue constant that motivated Fekete points in the first place. The experiments do not, and cannot, supply the missing sample-complexity theorem: $M=1.05N$ is a choice observed to work, not a choice proved to work. The authors themselves call the high-dimensional improvement over AFP modest. And the $d>1$ optimality theory still hangs on the hard-to-verify existence assumption of Theorem 3.1.

### Relation to the others

CFP is the deterministic, greedy counterpart of paper 22's randomised Christoffel-weighted least squares: the same weight $1/K_\Lambda$, and a candidate mesh literally drawn from paper 22's equilibrium-type densities, with the point selection done by pivoted QR rather than Monte Carlo. Where paper 45 and the induced-sampling line achieve non-asymptotic optimality by **exact induced sampling**, paper 28 achieves (conditional) exact optimality by **deterministic optimisation**; the paper's historical section explicitly places these as the two competing resolutions of the asymptotic-in-degree caveat left by paper 22.

**Papers 22, 45 and 28 therefore form a complete spectrum:**

| Route                            | Sample requirement   | Type of guarantee                  | Extra structure needed      |
| -------------------------------- | -------------------- | ---------------------------------- | --------------------------- |
| sample the orthogonality measure | $M\gtrsim N^2\log N$ | probabilistic                      | none                        |
| Christoffel-weighted sampling    | log-linear           | probabilistic                      | equilibrium measure, $K(z)$ |
| greedy (approximate Fekete)      | close to $N$         | deterministic, potential-theoretic | candidate pool, pivoted QR  |

## 24: designing sampling and preconditioning as a pair

### The idea

Paper 24 is the $\ell_1$ counterpart of paper 22. In sparse recovery the sample complexity is governed by the uniform bound $L$ of a **bounded orthonormal system**, and polynomial bases at high degree or on unbounded domains are **not** uniformly bounded, so Monte Carlo sampling from $w$ recovers badly once the degree is high.

The two existing fixes each fall short. Rauhut-Ward's "asymptotic sampling" — draw from Chebyshev and precondition the Vandermonde matrix — applies only to **bounded** random variables, and its accuracy degrades as the parameter dimension grows. Hampton-Doostan's coherence-optimal sampling is general but needs MCMC to draw from its measure. The paper wants one scheme that is explicit, applies to any admissible weight on a bounded domain and to a wide class of exponential weights on unbounded domains, and comes with a proved sample-count criterion.

**The key is the reciprocal structure**: the sampling density is $\propto w/(N\lambda)$, that is $w$ times the **reciprocal** of the Christoffel function, while the preconditioner weight is the Christoffel function **itself**. Sampling and weighting are inverse to each other, which is why the composite system becomes uniformly bounded — and why the "sampling density" and the "preconditioner" should be designed as a pair rather than choosing the sampling first and patching with a preconditioner afterwards.

### Setting

The Christoffel function and the preconditioner (the paper's eq. (9)) are

$$
(W)_{m,m}=N\lambda_\Lambda(Z^{(m)}),
\qquad
\lambda_\Lambda(Z)=\frac{1}{\sum_{i\in\Lambda}\varphi_i^2(Z)} .
$$

In the univariate case with $\Lambda=\{0,1,\dots,n\}$ and $N=n+1$ this is written $\lambda_{n+1}(z)=1/\sum_{k=0}^n\varphi_k^2(z)$.

> [!note] A printing inconsistency
> The paper uses $\lambda$ for the Christoffel function itself and $N\lambda$ for its scaled version, but line 4 of Algorithm 1 as printed reads $(W)_{m,m}=N/\lambda_\Lambda(Z^{(m)})$, which contradicts the $N\lambda_\Lambda$ of eq. (9); one of the two is a typo. This page reproduces the source as printed and flags it here rather than correcting it silently.

**Algorithm 1 (Christoffel Sparse Approximation, CSA)** has five steps: (1) generate $M$ i.i.d. samples $\{Z^{(m)}\}$ from the equilibrium density $v=\frac{\mathrm d\mu}{\mathrm dZ}$; (2) assemble $f$ with $f_m=f(Z^{(m)})$; (3) form the $M\times N(\Lambda)$ Vandermonde-like matrix $\Phi$ with $\Phi_{m,i}=\varphi_i(Z^{(m)})$; (4) compute the diagonal preconditioner from evaluations of the (scaled) Christoffel function; (5) solve the preconditioned basis-pursuit-denoising problem $\alpha^\star=\arg\min_\alpha\|\alpha\|_1$ subject to $\|\sqrt W\Phi\alpha-\sqrt Wf\|_2\le\varepsilon$.

**Three concrete regimes:**

- **CSA-a**: $Z$ Beta on $[-1,1]$ with shape parameters $\beta+1,\alpha+1\ge\frac12$ (Jacobi parameters $\alpha,\beta\ge-\frac12$) — sample from the **Chebyshev density**, independent of $n$, with $S_n\equiv[-1,1]$.
- **CSA-b**: $Z$ two-sided exponential on $\mathbb R$ with $w(z)=\exp(-|z|^\alpha)$, $\alpha>1$ — sample from an **expanded** equilibrium measure whose density $v_n$ and support $S_n$ both depend on $n$.
- **CSA-c**: $Z$ one-sided exponential on $[0,\infty)$ with $w(z)=\exp(-|z|^\alpha)$, $\alpha>\frac12$ — an analogous $n$-dependent expanded equilibrium measure.

In the multivariate unbounded construction the sampling ends with a rescaling $Z=4W_n$, producing samples on the set of points of $\mathbb R^d$ with $\ell_1$ norm at most $4n$.

The Gramian of the Christoffel-weighted basis (eq. (21)) is

$$
R_{k,\ell}=\int_{S_n}\varphi_{k-1}(z)\varphi_{\ell-1}(z)\,\bigl(N\lambda_n(z)\bigr)\,v_n(z)\,\mathrm dz,
\qquad 1\le k,\ell\le N=n+1 .
$$

$R$ is positive definite and each fixed entry converges to the corresponding entry of the identity as $n\to\infty$; $R^{1/2}$ is its unique symmetric positive definite square root, and $\|A\|_1$ is the induced $\ell_1$ matrix norm, the maximum column $\ell_1$ norm.

### Derivation

The design principle is one equation: the sampling density must satisfy $v\approx w/(N\lambda_N)$, with equality in the limit $N\to\infty$ given by the weighted pluripotential equilibrium measure — the asymptotic behaviour of the Christoffel function relative to $w$ is governed by that measure. This is the same asymptotic identity as in paper 22; only its use changes, from weighting to pairing.

Why the composite is uniformly bounded takes one line. The weighted basis $\sqrt{N\lambda_\Lambda}\,\varphi_i$ satisfies

$$
\sum_{i\in\Lambda}\bigl(\sqrt{N\lambda_\Lambda}\,\varphi_i\bigr)^2
=N\lambda_\Lambda\sum_{i\in\Lambda}\varphi_i^2=N
$$

**identically** — the same cancellation as $\tilde K_k\equiv N$ in paper 22. With uniform boundedness in hand, the restricted-isometry theory for bounded orthonormal systems takes over and the sample requirement has the shape $M\gtrsim L\,s\log^3(s)\log(N)$ with $L$ the uniform bound of the system. The entire difference between the paper's three regimes sits in that $L$.

The price sits in $R$: after weighting, the basis is orthonormal under $N\lambda_n v_n$ only once transformed by $R$. That is why the theorem below recovers $R^{1/2}\alpha$ rather than $\alpha$, and why the two terms of the error bound carry $\lambda_{\min}(R)$ and $\|R^{-1/2}\|_1$ as coefficients — **as in paper 22, it is the spectrum of $R$ that ends up in the bound.**

### Theorems

**Theorem 4.1 (the main recovery theorem, univariate).** Draw $M$ points i.i.d. from the equilibrium density $v_n$, form $\Phi_{ij}=\varphi_j(Z^{(i)})$ and the diagonal $W$ from eq. (9). Assume

$$
M\ \ge\ L(n)\,\bigl\|R^{-1/2}\bigr\|_1^2\;s\log^3(s)\log(N),
$$

where $L(n)$ is given by regime as follows:

| Regime | $L(n)$                                          | Branch expansion                                                             |
| ------ | ----------------------------------------------- | ----------------------------------------------------------------------------- |
| CSA-a  | $C(\alpha,\beta)$, uniformly in $n\ge1$         | independent of $n$                                                            |
| CSA-b  | $Cn^{\max\{1/\alpha,\,2/3\}}$, $C=C(\alpha)$    | $Cn^{2/3}$ for $\alpha\ge\tfrac32$; $Cn^{1/\alpha}$ for $1<\alpha<\tfrac32$   |
| CSA-c  | $Cn^{\max\{1/(2\alpha),\,2/3\}}$, $C=C(\alpha)$ | $Cn^{2/3}$ for $\alpha\ge\tfrac34$; printed as $Cn^{1/\alpha}$ for $\tfrac12<\alpha<\tfrac34$ |

Then with probability exceeding $1-N^{-\gamma\log^3(s)}$, for all polynomials $p(x)=\sum_j\alpha_j\varphi_j(x)$ and noisy samples $f=\Phi\alpha+\eta$ with $\|W\eta\|_\infty\le\varepsilon$, solving

$$
R^{1/2}\alpha^\star=\arg\min_\alpha\|R^{1/2}\alpha\|_1
\quad\text{s.t.}\quad
\|\sqrt W\Phi\alpha-\sqrt Wf\|_2\le\varepsilon
$$

gives

$$
\|\alpha-\alpha^\star\|_2\le C_1\frac{\sigma_s\bigl(R^{1/2}\alpha\bigr)_1}{\sqrt s\,\lambda_{\min}(R)}+C_2\frac{\varepsilon}{\sqrt{\lambda_{\min}(R)}},
$$

$$
\|\alpha-\alpha^\star\|_1\le D_1\sigma_s\bigl(R^{1/2}\alpha\bigr)_1\bigl\|R^{-1/2}\bigr\|_1+D_2\sqrt s\,\bigl\|R^{-1/2}\bigr\|_1\varepsilon .
$$

> [!note] An internal inconsistency
> In CSA-c the $\max$ expression reads $\max\{1/(2\alpha),2/3\}$ while the branch expansion immediately below prints $Cn^{1/\alpha}$ rather than $Cn^{1/(2\alpha)}$. This page reproduces the source as printed and flags it here rather than correcting it silently.

**The entire difference between bounded and unbounded sits in $L(n)$.** In the bounded, Chebyshev-like case $L(n)$ is uniform in $n$, so the requirement is just $M\gtrsim s\log^3(s)\log(N)$, near-linear in the sparsity $s$ and hence **optimal**; in the unbounded case it degrades to $M\gtrsim s\,n^{2/3}$, with an extra factor growing in the polynomial degree. Remark 4.2 states that within the bounded-orthonormal-system analysis used, the $n^{2/3}$ behaviour for $\alpha\ge\frac32$ and the $n^{1/\alpha}$ behaviour for $\alpha<\frac32$ are both **sharp**, with a similar statement for CSA-c.

**Corollary 4.1 (the clean Legendre case).** In CSA-a with $\alpha=\beta=0$ ($Z$ uniform, Legendre basis), $R=I$ for all $n$, so under $M\ge Cs\log^3(s)\log(N)$ the **untransformed** problem $\alpha^\star=\arg\min\|\alpha\|_1$ subject to $\|\sqrt W\Phi\alpha-\sqrt Wf\|_2\le\varepsilon$ already satisfies $\|\alpha-\alpha^\star\|_2\le C_1\sigma_s(\alpha)_1/\sqrt s+C_2\varepsilon$ and $\|\alpha-\alpha^\star\|_1\le D_1\sigma_s(\alpha)_1+D_2\sqrt s\varepsilon$.

> [!warning] Three gaps between the theorem and the implemented algorithm, stated by the paper
> First, **the theory is univariate only**. Remark 4.3 says the results generalise to tensor-product domains and weights if one samples tensor-product-wise from the respective univariate densities, in which case the criterion keeps its form but $L(n)$ becomes a **product of $d$ univariate factors** — the dimensional dependence is exponential. Remark 4.1 expects the Beta/Jacobi case to extend to almost any bounded weight on a compact interval.
> Second, the theorem recovers $R^{1/2}\alpha$ and **not** $\alpha$, because the PCE basis is orthonormal under $N\lambda_N(z)v_n(z)$ only after transformation by $R$ — whereas **the actual Algorithm 1 recovers $\alpha$ directly**, so the theorem does not literally analyse the implemented algorithm.
> Third, in the unbounded regimes the sampling density depends on $n$, so implementation needs the family of expanded equilibrium measures — the same class of difficulty paper 22 met on unbounded domains.

The paper explicitly advertises the framework as covering "bounded or unbounded" domains, so it also answers the unbounded-domain difficulty that [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|paper 11]] met in the least-squares setting.

### Numerical experiments

What has been checked here is the **subject of the four experiments and the qualitative outcomes the paper states**; the figures plot parameters against $N$, and the specific dimensions, sparsity levels, sample counts and error magnitudes have not been verified item by item, so no numbers are quoted.

| Experiment | Setup                                                          | Outcome as stated                                        |
| ---------- | ---------------------------------------------------------------- | ---------------------------------------------------------- |
| §6.1       | asymptotic (Chebyshev) sampling for Beta random variables      | validation of the CSA-a regime                            |
| §6.2       | Gaussian random variables, where the Chebyshev method **does not apply** | direct evidence that the unbounded case needs CSA-b       |
| §6.3       | manufactured sparse solutions                                   | a controlled check of recovery quality                    |
| §6.4       | an elliptic PDE with random inputs                              | an application-level check                                |

The abstract states two overall outcomes: CSA is superior to standard Monte Carlo in many situations of interest, and it gives **comparable or improved accuracy even against algorithms specialised to Legendre or Hermite bases** — a notable claim for a general-purpose scheme. The paper also notes that although its theory is univariate, the numerics show good performance in high-dimensional multivariate settings.

**What the experiments establish.** The most important point is that **all the multivariate evidence lives in the experiments**, since the theorem covers only the univariate case. And §6.2 shows the scheme works where Chebyshev preconditioning does not apply at all, which is precisely its selling point over Rauhut-Ward.

**What they do not establish.** They cannot repair the exponential dimensional dependence noted in Remark 4.3: the constant in the multivariate extension is a product of $d$ univariate factors, and experiments can only show good behaviour at the dimensions tried. Nor do they touch the gap between recovering $R^{1/2}\alpha$ and recovering $\alpha$.

### Relation to the others

Paper 24 is to $\ell_1$ what paper 22 is to least squares: the same authors, the same pairing of Christoffel weight and equilibrium measure, and the paper says explicitly that CSA "is based on a similar algorithm for discrete least-squares approximation that we introduced in [22]". Its bounds on Christoffel-weighted polynomials are used elsewhere to quantify sample counts for $\ell_1$ recovery when subsampling a tensor-product Gauss quadrature grid. Its preconditioning framework, general enough to cover bounded and unbounded domains, is the one later extended to include gradient measurements, and its CSA-a Chebyshev-sampling case coincides with Rauhut-Ward preconditioning. The reciprocal structure it makes explicit — sample from $\propto w/(N\lambda)$, weight by $N\lambda$ — is the same structure that paper 45's induced distribution realises **exactly** rather than asymptotically.

## 54: carrying the reasoning to reproducing kernel spaces

### The idea

The four papers above all work in polynomial spaces, but looking back at the idea stated at the top of this page, polynomials never really entered: what does the work is only "how concentrated a finite-dimensional function space is at each point". Wherever there is a reproducing kernel that quantity exists, and the same design principle ought to apply. **That is the level of abstraction paper 54 reveals**: replacing the polynomial space by a reproducing kernel Hilbert space, the role of the Christoffel function passes to the diagonal of the kernel (the power function), while the principle of "sample weighted by how concentrated the space is at each point" is unchanged. The paper connects that design problem to applications in uncertainty quantification.

### What has actually been verified here

> [!warning] This section is verified to a lower standard than the rest of the page
> The full text of this paper **has not been obtained** here, and the research notes behind this page contain no entry for it: no transcribed abstract, no theorems, no constants, no record of experiments. This section therefore keeps only what the title, the venue and the positioning above support, and **gives no setting, derivation, theorems or numerical experiments**. The statement that the kernel diagonal takes over the role of the Christoffel function is a reasonable placement at the level of the title, not something checked against the source. Bringing this section to the depth of papers 22 and 45 requires obtaining the paper.

### Relation to the others

It is the function-space generalisation of the thread: papers 22, 24 and 28 use the polynomial-space quantity $K_\Lambda=\sum_\alpha\psi_\alpha^2$, paper 45 uses the normalised reciprocal of that same quantity as a sampling density, and the kernel counterpart is the diagonal of the kernel. All three are instances of one design principle in different spaces.

> [!note] Coverage status
> Papers **22** and **45** are now developed in full from verified material: the construction, the chain of derivation (including paper 45's matrix Chernoff proof and the origin of the constant $C=2/\log(27/8e)$), the theorems with their hypotheses, and the experimental setups. Paper 45's three experiments (Example 5.2, Example 8.1 and the Section 9 thermal-diffusion PDE) have complete numerical setups and outcomes; paper **22**'s experiments have their design and qualitative outcomes only, since the numbers in its figures have not been verified here. Paper **28**'s Theorem 3.1, Lemma 3.1, Theorem 3.2 and Corollary 3.1, together with its full experimental protocol (a $10^4$ candidate pool, $M=1.05N$, 50 trials, three-way comparison), have been checked. Paper **24**'s main recovery theorem, the three cases of $L(n)$ and Corollary 4.1 have been checked, while its four experiments have their subject and qualitative outcomes only. For paper **54** the full text was never obtained, so only its title-level placement is kept, with no theorems, constants or experiments.
>
> Five source-level problems are preserved on this page: paper 22's three limitations on Theorem 4.3 and the non-vanishing $4\kappa^2(R)d^2(f)$ term in its main theorem; the inconsistent Gaussian normalising constant between paper 45's Example 2.2 and Example 8.1; the misprinted starting condition in paper 28's Theorem 3.2; the mismatch between line 4 of paper 24's Algorithm 1 and its eq. (9); and the internal inconsistency in paper 24's CSA-c branch expansion. The asymptotic Gaussian induced measure is treated as a **conjecture** everywhere it appears.

## Coverage check

| Item                                              | Paper | Status                                                              |
| ------------------------------------------------- | ----- | --------------------------------------------------------------------- |
| Stability factor and its basis independence       | 22    | $K(z)=\varphi^T\varphi$, $\lVert K\rVert_\infty/N$, invariance        |
| Decoupling sampling from orthogonality            | 22    | the insight, the $\tilde K_k\equiv N$ cancellation, the definition    |
| Five algorithm steps and row normalisation        | 22    | full steps and the equivalent reading                                 |
| Christoffel asymptotics and equilibrium measure   | 22    | Theorem 4.2, Corollary 4.1, Theorem 4.3 and its three limitations     |
| CLS stability and accuracy theorems               | 22    | Theorems 5.1 and 5.2, and the non-vanishing $4\kappa^2(R)d^2(f)$      |
| Experimental design and qualitative outcomes      | 22    | Figure 1, §6.1.1-6.1.3, §6.2.1-6.2.3; numbers not verified            |
| Sample complexity for three sampling routes       | 45    | quadratic, $N^{\log3/\log2}$, deterministic quadratic                 |
| Matrix Chernoff proof and origin of the constant  | 45    | the full argument for Theorem 6.1 and $C=2/\log(27/8e)$               |
| Induced distribution and the optimality argument  | 45    | Definition 7.1, attainment of the lower bound $N$, $M/\log M\ge C(r+1)N$ |
| Accuracy with truncation                          | 45    | Theorem 7.2 and Lemma 5.1                                             |
| Three numerical experiments                       | 45    | full setups and outcomes for Example 5.2, Example 8.1, the Section 9 PDE |
| Asymptotic induced measures                       | 45    | univariate and multivariate Chebyshev limits; Gaussian as conjecture  |
| Weighted approximate Fekete points, pivoted QR    | 28    | selection mechanism, meaning of weighting, enrichment and candidate pool |
| Row normalisation and equivalence of objectives   | 28    | Theorem 3.1, its existence premise, and this page's reading of it     |
| Greedy is optimal in one dimension                | 28    | Lemma 3.1, Theorem 3.2, Corollary 3.1 and the misprint                |
| Three-way comparison experiments                  | 28    | $\tilde M=10^4$, $M=1.05N$, 50 trials, four setups and outcomes       |
| Sampling and preconditioning designed as a pair   | 24    | the reciprocal structure, Algorithm 1, three regimes, the Gramian $R$ |
| Main recovery theorem, bounded versus unbounded   | 24    | Theorem 4.1, the three $L(n)$ cases, Corollary 4.1, two inconsistencies |
| Univariate theory against multivariate experiments | 24   | Remark 4.3's exponential dependence and four qualitative outcomes     |
| From polynomial spaces to kernel spaces           | 54    | the level of abstraction; full text never obtained, no theorems or experiments |

## Sources for this page

- A. Narayan, J. Jakeman, and T. Zhou, [_A Christoffel function weighted least squares algorithm for collocation approximations_](https://doi.org/10.1090/mcom/3192), Math. Comput. 86 (2017), pp. 1913-1947.
- J. Jakeman, A. Narayan, and T. Zhou, [_A generalized sampling and preconditioning scheme for sparse approximation of polynomial chaos expansions_](https://doi.org/10.1137/16M1063885), SIAM J. Sci. Comput. 39(3) (2017), pp. A1114-A1144.
- L. Guo, A. Narayan, L. Yan, and T. Zhou, [_Weighted approximate Fekete points: sampling for least-squares polynomial approximation_](https://doi.org/10.1137/17M1140960), SIAM J. Sci. Comput. 40(1) (2018), pp. A366-A387.
- L. Guo, A. Narayan, and T. Zhou, [_Constructing least-squares polynomial approximations_](https://doi.org/10.1137/18M1234151), SIAM Rev. 62(2) (2020), pp. 483-508.
- A. Narayan, L. Yan, and T. Zhou, [_Optimal design for kernel interpolation: applications to uncertainty quantification_](https://doi.org/10.1016/j.jcp.2020.110094), J. Comput. Phys. 430 (2021), 110094.
