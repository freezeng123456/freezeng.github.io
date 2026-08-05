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

| Paper | Symbol            | Meaning                                                                                                                                                                                                                  |
| ----- | ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 22    | $K$, $K_k$        | $\sum_\alpha\varphi_\alpha^2$, the reciprocal of the Christoffel function; $N/K_k$ is called "the (normalised) Christoffel function"                                                                                     |
| 24    | $\lambda_\Lambda$ | $1/\sum_{i\in\Lambda}\varphi_i^2$, the Christoffel function itself; the preconditioner uses $N\lambda_\Lambda$                                                                                                           |
| 28    | $K_\Lambda$       | $\sum_\alpha\psi_\alpha^2$, the reciprocal; the weight applied is $1/\sqrt{K_\Lambda}$                                                                                                                                   |
| 45    | $q^2$             | $\frac1N\sum_n v_n^2$, reciprocal and normalised; used directly as the sampling bias                                                                                                                                     |
| 54    | $K$, $\Lambda_N$  | here $K$ is the **kernel function itself** (the matrix is $\mathbf A=K(\Xi,\Xi)$) and $\Lambda_N$ is the **Lebesgue constant** — both collide with the rows above, so switch conventions wholesale when reading this one |

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

| Experiment | Setup                                                                                                                   | Outcome as stated                                                                       |
| ---------- | ----------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| Figure 1   | univariate symmetric Jacobi families ($\alpha=\beta$), $N=\dim\mathbb P_k=k+1$; plotted against $\beta$ and against $k$ | the standard Monte Carlo stability factor blows up with degree, which is the motivation |
| §6.1.1     | matrix stability on bounded domains: uniform/Legendre and other Jacobi cases                                            | CLS is more stable than standard Monte Carlo                                            |
| §6.1.2     | matrix stability on unbounded domains (Gaussian)                                                                        | supports the conjectured forms of the weighted equilibrium measure (Table 2)            |
| §6.1.3     | non-total-degree $\ell_p$ polynomial spaces                                                                             | CLS still performs well although the theory is framed for total-degree spaces           |
| §6.2.1     | an algebraic function                                                                                                   | accuracy of CLS against Monte Carlo                                                     |
| §6.2.2     | a heterogeneous diffusion equation in one spatial dimension                                                             | likewise                                                                                |
| §6.2.3     | an electrical resistor network                                                                                          | likewise                                                                                |

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

| Item           | Setup                                                                                                                                                        |
| -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| domain, weight | $D=\mathbb R^2$, $w(x)=\exp(-\lVert x\rVert^2)/\pi$                                                                                                          |
| space          | $\Lambda=\Lambda^{TD}(k)$, $N=\binom{k+2}{k}=(k+1)(k+2)/2$                                                                                                   |
| degrees        | $k=1,\dots,25$                                                                                                                                               |
| test function  | $f(x)=B\bigl(\lVert x/4-(0.2,-0.1)\rVert_2\bigr)$, $B$ a univariate bump function                                                                            |
| samples        | $M=10N$                                                                                                                                                      |
| trials         | 100                                                                                                                                                          |
| result         | sampling from $w$ gives an "extremely inaccurate" $g_N$ with very large $\eta_N$; induced sampling at $k=20$ ($N=231$) gives $\eta_N\sim1$ at moderate $M/N$ |

**Example 8.1 — the effect of dimension.**

| Item           | Setup                                                                                                                                                                                                                                     |
| -------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| domain, weight | $D=\mathbb R^d$, $w=\exp(-\lVert x\rVert_2^2)/\pi^{d/2}$                                                                                                                                                                                  |
| test function  | $f(x)=\prod_{j=1}^d\exp\bigl([x^{(j)}]^2/j\bigr)$, in product form so that $f_N$ stays computable in high dimension                                                                                                                       |
| space          | $\Lambda=\Lambda^{HC}(k)$                                                                                                                                                                                                                 |
| $(d,k)$        | $(4,20)$, $(8,10)$, $(20,5)$                                                                                                                                                                                                              |
| trials         | 100                                                                                                                                                                                                                                       |
| result         | induced sampling consistently outperforms standard sampling, **but the advantage diminishes as $d$ grows**; the qualitative behaviour of the induced-sampling $\eta_N$ is essentially unchanged with $d$, as Theorems 6.1 and 7.2 predict |

The paper is clear about why the advantage shrinks: large $d$ forces small $k$, and a low-degree space makes $\rho$ close to $w$ to begin with.

**Section 9 — a parametric thermal-diffusion equation.**

| Item                 | Setup                                                                                                          |
| -------------------- | -------------------------------------------------------------------------------------------------------------- |
| equation             | $-\nabla\cdot(a(y,x)\nabla u)=S$ on the unit square, $u\rvert_{\partial\Omega}=0$                              |
| source               | $S(y,x)=100\chi_F(y)$, $F$ a centred square subdomain of side $0.2$                                            |
| geometry             | four circular inclusions of radius $r=0.13$, symmetric about the centre                                        |
| coefficient          | $a(y,x)=1+\sum_{i=1}^4x^{(i)}\chi^i(y)$, so $d=4$                                                              |
| parameters           | $x^{(i)}\sim U(-0.99,0.2)$ i.i.d., giving a tensor Legendre basis                                              |
| quantity of interest | $f(x)=u\bigl((0.25,0.375),x\bigr)$                                                                             |
| discretisation       | $P_1$ finite elements in space; $V=V(\Lambda^{HC}(25))$                                                        |
| error estimate       | $Q=10{,}000$ reference samples, 100 trials                                                                     |
| result               | induced-distribution weighted least squares is "perhaps a bit better", but "the difference is not substantial" |

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

| Item           | Value                                                                                                         |
| -------------- | ------------------------------------------------------------------------------------------------------------- |
| candidate pool | $\tilde M=10^4$, half i.i.d. draws from $\rho$ and half from a degree-asymptotic density inspired by paper 22 |
| oversampling   | $\Delta N=\lfloor0.05N\rfloor$, that is $M=1.05N$ — **linear, 5% oversampling**                               |
| trials         | 50 per configuration; mean condition number reported with 20% and 80% quantiles                               |
| baselines      | MC (i.i.d. from $\rho$, unweighted $V(A_M,P)$), Fekete (AFP, unweighted), C-Fekete (CFP, weighted $V(A_M,Q)$) |

The second candidate ensemble is tensor-product Chebyshev for uniform $\rho$ on $[-1,1]^d$, and for Gaussian $\rho$ on $\mathbb R^d$ it is supported on the ball of radius $\sqrt{2n}$ with density

$$
C_d\Bigl(1-\tfrac{1}{2n}\sum_{k=1}^ds_k^2\Bigr)^{d/2},\qquad s\in\mathbb R^d,
$$

with $C_d$ a normalising constant. **This is the same functional form as the asymptotic Gaussian induced measure quoted in paper 45, and that form is only a conjecture**, so the candidate density here should likewise be read as a heuristic inspired by a conjecture. The authors note that weakly admissible meshes would be the natural candidates but that known constructions grow too fast in dimension.

| Experiment | Setup                                                                                                       | Result                                                                                                   |
| ---------- | ----------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| §5.1.1     | matrix stability, bounded: Legendre on $[-1,1]^d$, $d=2,6,10$, both TD and HC index sets                    | CFP gives notably more stable systems than AFP or MC                                                     |
| §5.1.2     | matrix stability, unbounded: Gaussian $\rho\propto\exp(-\lVert y\rVert_2^2)$, tensor Hermite, $d=2,6,10,25$ | CFP is more stable in every case, but the authors state the improvement is **modest in high dimensions** |
| §5.2(a)    | accuracy: $f(y)=\exp\bigl(-\sum_{j=1}^dy_j^2\bigr)$, Legendre approximation, $d=2$                          | CFP is comparable to but no worse than AFP, and much better than MC                                      |
| §5.2(b)    | accuracy: the stochastic elliptic equation below, $d=2,6$ (TD) and $d=10,25$ (HC)                           | likewise                                                                                                 |

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

| Regime | $L(n)$                                          | Branch expansion                                                                              |
| ------ | ----------------------------------------------- | --------------------------------------------------------------------------------------------- |
| CSA-a  | $C(\alpha,\beta)$, uniformly in $n\ge1$         | independent of $n$                                                                            |
| CSA-b  | $Cn^{\max\{1/\alpha,\,2/3\}}$, $C=C(\alpha)$    | $Cn^{2/3}$ for $\alpha\ge\tfrac32$; $Cn^{1/\alpha}$ for $1<\alpha<\tfrac32$                   |
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

| Experiment | Setup                                                                    | Outcome as stated                                   |
| ---------- | ------------------------------------------------------------------------ | --------------------------------------------------- |
| §6.1       | asymptotic (Chebyshev) sampling for Beta random variables                | validation of the CSA-a regime                      |
| §6.2       | Gaussian random variables, where the Chebyshev method **does not apply** | direct evidence that the unbounded case needs CSA-b |
| §6.3       | manufactured sparse solutions                                            | a controlled check of recovery quality              |
| §6.4       | an elliptic PDE with random inputs                                       | an application-level check                          |

The abstract states two overall outcomes: CSA is superior to standard Monte Carlo in many situations of interest, and it gives **comparable or improved accuracy even against algorithms specialised to Legendre or Hermite bases** — a notable claim for a general-purpose scheme. The paper also notes that although its theory is univariate, the numerics show good performance in high-dimensional multivariate settings.

**What the experiments establish.** The most important point is that **all the multivariate evidence lives in the experiments**, since the theorem covers only the univariate case. And §6.2 shows the scheme works where Chebyshev preconditioning does not apply at all, which is precisely its selling point over Rauhut-Ward.

**What they do not establish.** They cannot repair the exponential dimensional dependence noted in Remark 4.3: the constant in the multivariate extension is a product of $d$ univariate factors, and experiments can only show good behaviour at the dimensions tried. Nor do they touch the gap between recovering $R^{1/2}\alpha$ and recovering $\alpha$.

### Relation to the others

Paper 24 is to $\ell_1$ what paper 22 is to least squares: the same authors, the same pairing of Christoffel weight and equilibrium measure, and the paper says explicitly that CSA "is based on a similar algorithm for discrete least-squares approximation that we introduced in [22]". Its bounds on Christoffel-weighted polynomials are used elsewhere to quantify sample counts for $\ell_1$ recovery when subsampling a tensor-product Gauss quadrature grid. Its preconditioning framework, general enough to cover bounded and unbounded domains, is the one later extended to include gradient measurements, and its CSA-a Chebyshev-sampling case coincides with Rauhut-Ward preconditioning. The reciprocal structure it makes explicit — sample from $\propto w/(N\lambda)$, weight by $N\lambda$ — is the same structure that paper 45's induced distribution realises **exactly** rather than asymptotically.

## 54: carrying the reasoning to reproducing kernel spaces

### The idea

The four papers above all work in polynomial spaces, but looking back at the idea stated at the top of this page, polynomials never really entered: what does the work is only "how concentrated a finite-dimensional function space is at each point". Wherever there is a reproducing kernel that quantity exists, and the same design principle ought to apply. **That is the level of abstraction paper 54 reveals**: replacing the polynomial space by a reproducing kernel Hilbert space, the role of the Christoffel function passes to the diagonal of the kernel, while the design principle — put the points where the space concentrates — is unchanged.

To make that concrete one has to say **which** diagonal, or the analogy comes out distorted. For a translation-invariant radial kernel $\mathbf A_{ii}=\Phi(0)$ for every $i$, so **the diagonal of the kernel itself is constant** and there is nothing to level out. What actually corresponds to the polynomial-space quantity $K_\Lambda=\sum_\alpha\psi_\alpha^2$ is the reproducing-kernel diagonal of **the subspace spanned by the points already chosen**, namely $k(z)^{T}\mathbf A_N^{-1}k(z)$, where $k(z)$ is the vector of kernel values between $z$ and the chosen points. Its complement

$$
P_N^2(z)=\Phi(0)-k(z)^{T}\mathbf A_N^{-1}k(z)
$$

is the squared **power function**, equivalently the Gaussian-process posterior variance given the data already placed. The paper picks the next point by maximising exactly that. **So "put the points where the space concentrates" reads, in kernel language, as "put the next point where the current span cannot reach" — two faces of one rule.** The paper positions itself the same way: the greedy iteration is called approximate Fekete points in the polynomial literature, and "for kernel interpolation it is a power function maximization approach".

**Why instability grows with the number of points.** The paper gives a concrete mechanism (Remark 2): when two centres sit too close together, identically shaped basis functions centred at them take nearly equal values at all the nodes, so the interpolation matrix acquires two nearly identical columns. On a bounded parameter domain, growing $N$ forces the fill distance down, so clustering is unavoidable — the ill-conditioning is not bad luck but the necessary consequence of adding points. The paper observes it at $N\sim\mathcal O(10^2)$ for all the usual choices: Sobol' points, low-discrepancy points, randomly distributed points. The shape parameter offers no way out, being squeezed from both sides: too large and the interpolant is inaccurate, too small and the matrix is so ill-conditioned that linear solvers perform poorly, and how to choose the best $\epsilon$ is, the paper says plainly, still an open problem. The gradient-enhanced case is worse still: once derivative observations are folded in, the condition number of the system matrix $\mathbf B$ grows **much faster** than that of the ordinary interpolation matrix $\mathbf A$.

**Why greedy selection from a large candidate pool addresses it.** The paper trades an intractable question for a solvable one: rather than optimising $\epsilon$ or searching a continuum for optimal centres, it asks how to choose an "optimal" subset given a large set of candidate points. Two things follow. First, the selection involves only linear algebra, amounting essentially to Cholesky-type decompositions of the design matrix. Second — the practically decisive point — the selection is **data-independent**: $\mathbf A$ and $\mathbf B$ depend only on where the points are, not on the values of $u$, so the whole design can be computed offline before any expensive simulation is run.

### Setting

The parametric problem is stated in general form. With $Z=(Z^{(1)},\dots,Z^{(d)})\in I_Z\subseteq\mathbb R^d$ the uncertain inputs,

$$
\begin{cases}
u_t(x,t,Z)=\mathcal L(u), & \mathcal D\times(0,T]\times I_Z,\\
\mathcal B(u)=0, & \partial\mathcal D\times(0,T]\times I_Z,\\
u=u_0, & \mathcal D\times\{t=0\}\times I_Z,
\end{cases}
$$

with $\mathcal D\subset\mathbb R^l$, $l=1,2,3$, the physical domain. Fixing $x$ and $t$ and writing $u_j=u(x,t,z_j)$, the goal is to build $u_N(Z)\approx u(Z)$ from the pairs $(z_j,u_j)$, $j=1,\dots,N$.

**Kernel interpolation.** Pick a translation-invariant radial kernel $K=\Phi(\epsilon\lVert\cdot-\cdot\rVert)$, with $\Phi$ the radial basis function and $\epsilon$ a shape parameter. With centres $\Xi=\{z_j\}_{j=1}^N$ the trial space is

$$
\mathcal U_\Xi=\mathrm{span}\bigl\{\Phi(\epsilon\lVert\cdot-z_j\rVert)\ \big|\ z_j\in\Xi\bigr\},
\qquad
u_N(Z)=\sum_{j=1}^Nc_j\,\Phi(\epsilon\lVert Z-z_j\rVert).
$$

The paper restricts itself to **positive definite** kernels, positive definiteness guaranteeing existence and uniqueness of the approximation; the three used are the Gaussian $\Phi(r)=\exp(-r^2)$, the inverse multiquadric (IMQ) $\Phi(r)=1/\sqrt{1+r^2}$, and a compactly supported (CS) family. The interpolation conditions $u_N(z_i)=u_i$ give

$$
\mathbf{Ac}=\mathbf u,\qquad \mathbf A=K(\Xi,\Xi),\quad \mathbf A_{ij}=\Phi(\epsilon\lVert z_i-z_j\rVert),
$$

with $\mathbf A$ symmetric and $\mathbf c$ unique when $\mathbf A$ is invertible. Remark 1 records the Gaussian-process reading: finding the best linear unbiased estimator for a Gaussian process with covariance $K$ from the realisation $\mathbf u$ on $\Xi$ produces the same function (Appendix A derives simple kriging and concludes $u_N(Z)=K(Z,\Xi)^TK(\Xi,\Xi)^{-1}\mathbf u=\hat f(Z)$, the two estimators being **identical**). Remark 2 is the ill-conditioning mechanism quoted above.

**LOOCV for the shape parameter.** The paper declines to analyse the choice of $\epsilon$ systematically and uses leave-one-out cross validation instead. The error vector has entries

$$
e_i(\epsilon)=\frac{c_i}{\mathbf A^{-1}_{ii}},
\qquad
\lVert e(\epsilon^*)\rVert=\min_\epsilon\lVert e(\epsilon)\rVert,
$$

where $\mathbf A^{-1}_{ii}$ is the $i$th diagonal entry of $\mathbf A^{-1}$. The point of this form is that the leave-one-out error can be read off the **already computed** $\mathbf c$ and the diagonal of $\mathbf A^{-1}$, with no need to actually solve $N$ reduced systems.

**The gradient-enhanced formulation.** Given $\{z_i,u'_m(z_i)\}_{i=1}^N$ for $m=1,\dots,d$ alongside $\{z_i,u(z_i)\}_{i=1}^N$, where $u'_m=\partial u/\partial z^{(m)}$, the interpolant becomes

$$
u_N(Z)=\sum_{j=1}^Nc_j\Phi(\epsilon\lVert Z-z_j\rVert)-\sum_{m=1}^d\sum_{j=1}^N\beta_{m,j}\Phi'_m(\epsilon\lVert Z-z_j\rVert),
$$

required to satisfy the $N(d+1)$ generalised interpolation conditions $\lambda_iu=\lambda_iu_N$, where each $\lambda_i$ is either point evaluation at $z_i$ or evaluation of some derivative at $z_i$. The system matrix is block structured,

$$
\mathbf B=\begin{pmatrix}
\mathbf A_{0,0}&\mathbf A_{0,1}&\cdots&\mathbf A_{0,d}\\
\mathbf A_{1,0}&\mathbf A_{1,1}&\cdots&\mathbf A_{1,d}\\
\vdots&\vdots&\ddots&\vdots\\
\mathbf A_{d,0}&\mathbf A_{d,1}&\cdots&\mathbf A_{d,d}
\end{pmatrix},
\qquad
(\mathbf A_{m,n})_{i,j}=\begin{cases}
\Phi(\epsilon\lVert z_i-z_j\rVert), & m=n=0,\\
-\Phi'_n(\epsilon\lVert z_i-z_j\rVert), & m=0,\ n\neq0,\\
\Phi'_m(\epsilon\lVert z_i-z_j\rVert), & m\neq0,\ n=0,\\
-\Phi''_{m,n}(\epsilon\lVert z_i-z_j\rVert), & m\neq0,\ n\neq0,
\end{cases}
$$

and $\mathbf B$ is a symmetric $(d+1)N\times(d+1)N$ matrix. Appendix B notes that when the RBF kernel and the covariance kernel coincide, $\mathbf B$ and the joint covariance matrix of the gradient-enhanced GP emulator are **identical**. That identity is what gives the greedy criterion of the next section its meaning: the quantity being maximised acquires a probabilistic reading.

### Derivation

**Starting from Fekete points.** Write the interpolant in cardinal form,

$$
u_N(Z)=\sum_{j=1}^Nc_j\Phi(\epsilon\lVert Z-z_j\rVert)=\sum_{j=1}^Nd_j\,\ell_j(Z),
\qquad
\ell_j(Z)=\frac{\det\mathbf A(z_1,\dots,z_{j-1},Z,z_{j+1},\dots,z_N)}{\det\mathbf A(z_1,\dots,z_N)},
$$

with Lebesgue constant $\Lambda_N:=\max_{Z\in I_Z}\sum_{j=1}^N\lvert\ell_j(Z)\rvert$. **Fekete points are the configuration maximising the Vandermonde-like determinant,**

$$
\Xi^*=\arg\max_{\Xi=\{z_1,\dots,z_N\}\subset I_Z}\bigl\lvert\det\mathbf A(z_1,\dots,z_N)\bigr\rvert .
$$

By the cardinal formula, at $\Xi^*$ the numerator of each $\ell_j$ is a determinant with one $z_j$ replaced by $Z$, so its modulus cannot exceed the denominator; hence $\lvert\ell_j(Z)\rvert\le1$ and $\Lambda_N\le N$. The paper adds that in practice the observed growth of the Lebesgue constant for Fekete points is frequently sublinear. The difficulty is that outside special polynomial cases there is no known way to characterise or compute these points explicitly, and solving the optimisation directly is "a daunting task". The standard relaxation is to go greedy, adding one point at a time so as to maximise the determinant.

> [!note] A misprint in the greedy formula of §3.1
> The paper prints the greedy iteration as $z_{N+1}=\arg\max_{z\in I_Z}\lvert\det\mathbf A(z_1,\dots,z_N)\rvert$, in which **$z$ does not occur in the maximand at all**, so read literally it is a constant in the optimisation variable. It should be $\det\mathbf A(z_1,\dots,z_N,z)$, the determinant after the new point is folded in — which is exactly the form of the corresponding expression written for $\mathbf B$ in §3.2. This page states it in the corrected form and flags it here. (A lighter one: the caption of Figure 4 reads "condition numbers with respect to shape parameters with respect to the number of sample points $N$", while the abscissa of that figure is $N$ and $\epsilon=3,5$ are the fixed parameters.)

The paper says explicitly that the polynomial literature calls this the **approximate Fekete point** approach, and that **for kernel interpolation it is a power function maximisation approach** — the source of the equivalence stated above.

**What one greedy step actually maximises.** §3.2 writes the greedy iteration on the gradient-enhanced matrix:

$$
z_{N+1}:=\arg\max_{z\in I_Z}\det\boldsymbol B(z_1,\dots,z_N,z).
$$

Evaluating that determinant from the definition would mean rebuilding the whole of $\boldsymbol B$, which is too expensive, so the paper takes it apart with a Schur complement. Let $\boldsymbol P_{N+1}$ be the permutation moving the $d+1$ entries contributed by the new point to the end while keeping the other $N(d+1)$ in order. Then

$$
\boldsymbol P_{N+1}\boldsymbol B\boldsymbol P_{N+1}^{-T}
=\begin{pmatrix}\boldsymbol B_N&\boldsymbol W(z)\\\boldsymbol W^T(z)&\boldsymbol B(z)\end{pmatrix},
\qquad
\boldsymbol B_N:=\boldsymbol B(z_1,\dots,z_N),
$$

so that

$$
\det\boldsymbol B(z_1,\dots,z_N,z)
=\det\boldsymbol B_N\cdot\det\bigl(\boldsymbol B(z)-\boldsymbol W^T\boldsymbol B_N^{-1}\boldsymbol W\bigr).
$$

**The factor $\det\boldsymbol B_N$ does not depend on $z$ and drops out entirely**, so the greedy step is equivalent to

$$
z_{N+1}=\arg\max_{z\in I_Z}F(z),
\qquad
F(z):=\det\bigl(\boldsymbol B(z)-\boldsymbol W^T\boldsymbol B_N^{-1}\boldsymbol W\bigr).
$$

That answers what each pivot maximises: **the determinant of the Schur complement, that is, of the $(d+1)\times(d+1)$ conditional covariance block at the candidate $z$ given every value and gradient observation already placed.** By the identity of Appendix B this block is precisely the posterior covariance of the gradient-enhanced Gaussian process at $z$, so the greedy rule reads "put the next point where the current model is least certain". Specialising to the scalar, gradient-free case — $\boldsymbol B(z)$ becomes $\Phi(0)$ and $\boldsymbol W$ becomes $k(z)$ — the Schur complement is $\Phi(0)-k(z)^T\mathbf A_N^{-1}k(z)=P_N^2(z)$, the squared power function from the top of this section.

**Why this is a pivoted Cholesky factorisation.** The dominant cost of evaluating $F$ is solving $\boldsymbol B_N\boldsymbol X=\boldsymbol W$ with $\boldsymbol X\in\mathbb R^{N(d+1)\times(d+1)}$; the remaining $(d+1)\times(d+1)$ determinant is comparatively negligible. The paper keeps the Cholesky factor of the permuted $\boldsymbol B_N$,

$$
\boldsymbol P_N\boldsymbol B_N\boldsymbol P_N^T=\boldsymbol L_N\boldsymbol L_N^T,
$$

and stores $\boldsymbol L_N^{-1}$, so each evaluation of $F$ needs one application of $\boldsymbol L_N^{-1}$ at cost $\mathcal O(d^3N^2)$ and the optimisation becomes

$$
z_{N+1}=\arg\max_{z\in I_Z}\det\bigl(\boldsymbol B(z)-\boldsymbol V^T\boldsymbol V\bigr),
\qquad
\boldsymbol V:=\boldsymbol L_N^{-1}\boldsymbol W(z).
$$

Once $z_{N+1}$ is fixed the factor is updated blockwise rather than refactorised:

$$
\boldsymbol L_{N+1}=\begin{pmatrix}\boldsymbol L_N&\boldsymbol 0\\\boldsymbol V^T&\widetilde{\boldsymbol L}\end{pmatrix},
\qquad
\widetilde{\boldsymbol L}\widetilde{\boldsymbol L}^T=\boldsymbol B(z_{N+1})-\boldsymbol V^T\boldsymbol V,
$$

$$
\boldsymbol L_{N+1}^{-1}=\begin{pmatrix}\boldsymbol L_N^{-1}&\boldsymbol 0\\-\widetilde{\boldsymbol L}^{-1}\boldsymbol V^T\boldsymbol L_N^{-1}&\widetilde{\boldsymbol L}^{-1}\end{pmatrix}.
$$

A step therefore consists of three things: (i) compute $\boldsymbol V(z_{N+1})$; (ii) invert the $(d+1)\times(d+1)$ Cholesky factor $\widetilde{\boldsymbol L}$ of the Schur complement; (iii) assemble $\boldsymbol L_{N+1}^{-1}$, which costs one further application of $\boldsymbol L_N^{-1}$. **This is a pivoted Cholesky factorisation and nothing else**: $\boldsymbol L_N$ is the Cholesky factor of the kernel matrix on the selected points, and each step promotes the "largest" remaining block of the Schur complement to be the next pivot. The paper's own summary describes the method in exactly those terms — an efficient Cholesky decomposition with pivoting on the Vandermonde-like interpolation matrix, used to choose the sample points. The last implementation concession is that maximisation over the continuum $I_Z$ is replaced by maximisation over a discrete candidate set; the paper emphasises that in practice this is what is always computed.

### Theorems

**This paper contains no theorems.** There is not a single theorem, lemma, proposition, corollary or proof environment anywhere in the text. The only general statements it makes are of two kinds: the classical Fekete-point consequence $\lvert\ell_j\rvert\le1\Rightarrow\Lambda_N\le N$, and the identities of Appendices A and B between kernel interpolation and (gradient-enhanced) Gaussian-process estimation. Both are asserted and attributed to the literature rather than proved here.

This has to be said clearly, because "optimal design" in the title and "quasi-optimal" in the body both invite a quantitative reading. **They carry none.** The paper gives no comparison between the greedy set and the true maximiser, no bound on the condition number, no Lebesgue-constant bound for the selected points, and no convergence rate. "Quasi-optimal" is simply the name for a greedy relaxation of a combinatorial optimisation that cannot be solved directly. All the quantitative evidence lives in the figures.

> [!note] This is where the contrast with paper 28 is sharpest
> In the polynomial setting paper 28 proves two substantive results about a structurally identical greedy iteration: that determinant optimality and condition-number optimality coincide in the weighted space (Theorem 3.1, under an existence premise that is hard to verify for $d>1$), and that in one dimension the greedy iteration attains the optimum exactly (Theorem 3.2). Paper 54 carries the same strategy into kernel spaces **without bringing any theorem along**. Whether the kernel analogues hold is neither proved nor discussed.

### Numerical experiments

**Shared protocol.** Every experiment uses the same setup.

| Item           | Value                                                                                                            |
| -------------- | ---------------------------------------------------------------------------------------------------------------- |
| Candidate pool | $M=10^4$ centres taken at random in $I_Z$; "uniform" (evenly spaced in $I_Z$) also compared                      |
| Selection      | $N$ RBF centres chosen from the pool by the §3.2 greedy iteration (called "Cholesky" in the paper)               |
| Baselines      | random, Sobol' and Halton points; the PDE example adds sparse grids                                              |
| Error          | $E_{\ell_2}=\bigl(\tfrac1Q\sum_{i=1}^Q\lvert s(z_i)-u(z_i)\rvert^2\bigr)^{1/2}$ over $Q=1000$ random test points |
| Trials         | 50 per configuration; mean reported with 20% and 80% quantiles                                                   |
| Kernels        | Gaussian $\Phi(r)=e^{-r^2}$, IMQ $\Phi(r)=1/\sqrt{1+r^2}$, Wendland CS with $l=\lfloor d/2\rfloor+4$             |

> [!warning] No condition numbers or errors can be quoted in this section
> **Every numerical result in this paper is presented as a figure; there is not one data table.** Figures 1-14 are curves, and no condition number or error value appears as a number anywhere in the body. (The HTML rendering of the paper contains forty-odd table elements, but every one of them is a typesetting container for a displayed equation rather than a table of data.) What follows therefore reports the experimental setups and the paper's stated conclusions, both of which can be checked item by item, and quotes no condition-number magnitudes or error values.

**Condition-number experiments.** The first two figures establish the problem before the method is proposed; the last three test the method.

| Figure | Object                                                          | Abscissa                   | Fixed parameters          |
| ------ | --------------------------------------------------------------- | -------------------------- | ------------------------- |
| Fig. 1 | $\mathbf A$ (Lagrange, solid) and $\mathbf B$ (Hermite, dashed) | shape parameter $\epsilon$ | $N=50,100,300$; $d=2$     |
| Fig. 2 | the same                                                        | number of points $N$       | $\epsilon=0.1,1,3$; $d=2$ |
| Fig. 3 | $\mathbf A$, four selection methods compared                    | $\epsilon$                 | $N=100,300$; $d=2$        |
| Fig. 4 | $\mathbf A$, four selection methods compared                    | $N$                        | $\epsilon=3,5$; $d=2$     |
| Fig. 5 | $\mathbf A$, different candidate pools                          | $N$                        | $\epsilon=3,5$            |

The paper's conclusion from Figures 1 and 2 is that the design matrix becomes more singular both as the shape parameter approaches zero and as the number of points grows, with $\mathbf B$ degrading far faster than $\mathbf A$. From Figures 3 and 4: **for both the Gaussian and the IMQ kernel the proposed algorithm is much more stable than the other sampling methods** — that sentence in the source names only those two kernels, and the compactly supported kernel, although plotted in the same figures, is not covered by the claim. Figure 5 concludes that the choice of candidate points, random or evenly spaced, does not dramatically affect the method's performance.

**Accuracy experiments.** The first test function is Franke's benchmark on $[0,1]^2$,

$$
\begin{aligned}
u(z)=&\tfrac34e^{-((9z^{(1)}-2)^2+(9z^{(2)}-2)^2)/4}+\tfrac34e^{-(9z^{(1)}+1)^2/49-(9z^{(2)}+1)^2/10}\\
&+\tfrac12e^{-((9z^{(1)}-7)^2+(9z^{(2)}-3)^2)/4}-\tfrac15e^{-(9z^{(1)}-4)^2-(9z^{(2)}-7)^2}.
\end{aligned}
$$

> [!note] The second exponent is misprinted
> The paper prints that term as $-(9z^{(2)}+1)^2)/10$, with unmatched parentheses; moreover the standard Franke function is **linear** in this term, $-(9z^{(2)}+1)/10$, not squared. The display above transcribes the paper (with the parentheses balanced) and flags the discrepancy here.

The second is a stochastic elliptic equation in one spatial dimension — **the same benchmark used in §5.2(b) of paper 28**:

$$
-\frac{\mathrm d}{\mathrm dx}\Bigl[\kappa(x,z)\frac{\mathrm du}{\mathrm dx}(x,z)\Bigr]=f,
\quad(x,z)\in(0,1)\times\mathbb R^d,
\qquad u(0,z)=u(1,z)=0,\quad f=2,
$$

$$
\kappa(x,z)=1+\sigma\sum_{k=1}^d\frac{1}{k^2\pi^2}\cos(2\pi kx)\,z^{(k)},
$$

with quantity of interest $u(z)=u(0.5,z)$ and $z^{(i)}\sim U[-1,1]$ i.i.d. (The value of $\sigma$ is not given in the paper, so none is quoted.)

| Figure  | Test problem                 | Dimension | Baselines                                           | Notes                               |
| ------- | ---------------------------- | --------- | --------------------------------------------------- | ----------------------------------- |
| Fig. 6  | Franke's function, $[0,1]^2$ | $d=2$     | random, Sobol', Halton                              | $\epsilon=3,5$                      |
| Fig. 7  | the same, other candidates   | $d=2$     | uniform and other candidate pools                   | $\epsilon=3,5$                      |
| Fig. 8  | the same                     | $d=2$     | several $\epsilon$, plus $\epsilon$ from LOOCV      | left condition number, right error  |
| Fig. 9  | stochastic elliptic equation | $d=3$     | **plus sparse grids** (Legendre, total order $k=8$) | left Gaussian, middle IMQ, right CS |
| Fig. 10 | stochastic elliptic equation | $d=6$     | **plus sparse grids** (Legendre, total order $k=4$) | left Gaussian, middle IMQ, right CS |

The paper draws two conclusions from Figure 6, the second more interesting than the first: the proposed algorithm is more accurate than the other sampling methods, and **it shows a clear convergence pattern as $N$ grows, whereas the error profiles of the other methods make it clear that supplying more points does not always improve accuracy**. Figure 8 concludes that smaller $\epsilon$ makes the matrix ill-conditioned, and that using LOOCV to pick $\epsilon$ gives very good results.

**Gradient-enhanced experiments.** Three test functions:

| Figure  | Test function                                                                            | Dimension | Domain     |
| ------- | ---------------------------------------------------------------------------------------- | --------- | ---------- |
| Fig. 12 | corner peak, $u(z)=\bigl(1+\sum_{i=1}^d\omega_iz^{(i)}\bigr)^{-(d+1)}$, $\omega_i=1/i^2$ | $d=2$     | $[0,1]^2$  |
| Fig. 13 | Rastrigin, $u(z)=20+\sum_{i=1}^2\bigl((z^{(i)})^2-10\cos(2\pi z^{(i)})\bigr)$            | $d=2$     | $[-4,4]^2$ |
| Fig. 14 | Friedman, $u(z)=10\sin(\pi z^{(1)}z^{(2)})+20(z^{(3)}-0.5)^2+10z^{(4)}+5z^{(5)}$         | $d=5$     | not stated |

(Figure 11 plots the Rastrigin function itself and carries no results.) The conclusion is that under the proposed algorithm the gradient-enhanced design matrix $\mathbf B$ can be well conditioned, and that accuracy again exceeds the other sampling methods.

**What the sparse-grid comparison actually shows.** The abstract says the method "can outperform sparse grid methods in many interesting cases". Checked section by section, the comparison against sparse grids **occurs exactly once**: the stochastic elliptic equation of §4.1.2, at $d=3$ (Legendre, total order $k=8$) and $d=6$ ($k=4$), reported in Figures 9 and 10 with one panel per kernel. The paper's own sentence is that the RBF approximation methods are "notably superior to the sparse grids method". **Nowhere in the paper is a case reported in which the method loses to sparse grids, and nowhere is one reported in which it loses to random, Sobol' or Halton points either.**

So the hedge means something different from what it sounds like. It is not "we tried many cases and won in many of them", with the implication that some were lost; it is "we ran one comparison, won it, and are careful not to promote that into a general claim". **The evidence behind it is one benchmark problem, two parameter dimensions and three kernels.** Nor does the paper put the two methods' **costs** side by side: the abscissa in Figures 9 and 10 is the number of points $N$, while the greedy selection itself must scan $M=10^4$ candidates at $\mathcal O(d^3N^2)$ per candidate per step, which is not counted; and sparse grids do not offer an arbitrary node count either, so the two curves are not really aligned on that axis.

**What these experiments establish.** First, Figures 1 and 2 turn "ill-conditioning grows with the number of points" from a worry into a reproducible fact, and show that gradient enhancement makes it substantially worse — that is the paper's starting point and the reason the greedy iteration is run on $\mathbf B$ rather than $\mathbf A$. Second, Figures 3 and 4 show that at a fixed number of points the choice of points alone opens a visible gap in the condition number, and does so without any information about $u$. Third, Figures 5 and 7 show insensitivity to how the candidate pool is drawn, which is a necessary robustness check for a two-stage "pick candidates, then pick points" method. Fourth, the observation in Figure 6 that the other methods do not reliably improve as points are added is itself the empirical form of the stability argument.

**What they do not establish.** The central one is stated above: there is no theorem, so "quasi-optimal" has no quantitative content and every conclusion is a curve in a figure. Next, the shape parameter is fixed at $\epsilon\in\{3,5\}$ everywhere except Figure 8, even though the paper itself calls the optimal $\epsilon$ an open problem, so LOOCV is tested in exactly one figure. Third, every test has $d\le6$ while the stated motivation is high-dimensional UQ; the candidate pool is fixed at $M=10^4$ without adjustment for $d$, which is thin coverage at $d=5$ or $6$, and beyond "random versus evenly spaced" the paper does not study the effect of $M$. Fourth, the largest point count appearing explicitly in the text is $N=300$ (Figures 1 and 3), and the figures with $N$ on the abscissa give no range in the body or the captions, so how far the claim of "postponing the instability" extends cannot be determined from the text.

### Relation to the others

**The thing most worth seeing is the structural identity with paper 28.** The two papers do the same thing: instead of optimising over a continuum, assemble a large candidate pool and then peel the points off one at a time with a single pivoted matrix factorisation.

| Item               | Paper 28                                   | Paper 54                                     |
| ------------------ | ------------------------------------------ | -------------------------------------------- |
| Function space     | polynomial space (weighted $Q$)            | reproducing kernel Hilbert space             |
| Quantity maximised | Vandermonde-like determinant               | kernel-matrix determinant (Schur complement) |
| Factorisation      | **column-pivoted QR** of $V^T$             | **pivoted Cholesky** of $\boldsymbol B$      |
| Candidate pool     | $\tilde M=10^4$                            | $M=10^4$                                     |
| Statistics         | 50 trials, 20%/80% quantiles               | 50 trials, 20%/80% quantiles                 |
| Theorems           | Theorems 3.1, 3.2 (conditional optimality) | none                                         |

Both factorisations can implement the greedy iteration for the same reason: **each writes the determinant as a product over pivots** — QR gives $\prod\lvert r_{ii}\rvert$, Cholesky gives the product of the diagonal blocks — so "add the point that maximises the determinant" automatically becomes "take the next pivot". The only difference is the function space, and with it the kind of factorisation the matrix admits: a non-symmetric design matrix calls for QR, a symmetric positive definite kernel matrix for Cholesky.

The difference in how each handles weighting follows the same line. Paper 28 must insert the Christoffel weight $1/\sqrt{K_\Lambda}$, because polynomial bases concentrate very unevenly and only after weighting do the rows of the Vandermonde matrix have norm exactly 1, which is what lets Hadamard's inequality bite. On the kernel side, translation invariance makes $\mathbf A_{ii}=\Phi(0)$ constant — **that comes free** — so there is no weight to insert. All the remaining spatial variation moves into the residual diagonal after conditioning, the power function, which is exactly what paper 54's greedy iteration maximises. One paper is greedy on a **weighted** determinant, the other greedy on a **residual** diagonal; they are one design principle in two geometries.

The two even share a PDE benchmark: the one-dimensional stochastic elliptic equation with KL-type diffusivity $\kappa=1+\sigma\sum_k\frac{1}{k^2\pi^2}\cos(2\pi kx)z^{(k)}$, quantity of interest $u(0.5,\cdot)$, and $z$ uniform on $[-1,1]^d$. The author lists overlap and so do the experimental protocols.

**The gradient-enhanced side points to another page.** One adjoint solve delivers the function value together with all $d$ partial derivatives, multiplying the number of rows of the measurement matrix by $d+1$ at very little extra cost — the same fact exploited by papers 29 and 32 on [[en/computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce|Sparse recovery and data-driven chaos]], except that there the extra information feeds $\ell_1$ sparse recovery while here it feeds a Gaussian process emulator. **The price paid is of the same kind on both routes.** Paper 32 finds that naively stacking gradient rows destroys mean isotropy and has to be repaired with row preconditioning and column normalisation; paper 54 meets the conditioning version of the same damage, its Figures 1 and 2 recording that the condition number of $\mathbf B$ grows much faster than that of $\mathbf A$. The two repairs sit at opposite ends: paper 32 fixes the **matrix** by preconditioning it, paper 54 fixes the **points** by choosing where they go.

**Within this page it is the function-space generalisation of the thread.** Papers 22, 24 and 28 use the polynomial-space quantity $K_\Lambda=\sum_\alpha\psi_\alpha^2$, paper 45 uses the normalised reciprocal of that same quantity as a sampling density, and the kernel counterpart is the reproducing-kernel diagonal of the already-selected subspace together with its complement, the power function. All four are instances of one design principle in different spaces — **with the single difference that only paper 54 arrives without a theorem attached.**

> [!note] Coverage status
> Papers **22** and **45** are now developed in full from verified material: the construction, the chain of derivation (including paper 45's matrix Chernoff proof and the origin of the constant $C=2/\log(27/8e)$), the theorems with their hypotheses, and the experimental setups. Paper 45's three experiments (Example 5.2, Example 8.1 and the Section 9 thermal-diffusion PDE) have complete numerical setups and outcomes; paper **22**'s experiments have their design and qualitative outcomes only, since the numbers in its figures have not been verified here. Paper **28**'s Theorem 3.1, Lemma 3.1, Theorem 3.2 and Corollary 3.1, together with its full experimental protocol (a $10^4$ candidate pool, $M=1.05N$, 50 trials, three-way comparison), have been checked. Paper **24**'s main recovery theorem, the three cases of $L(n)$ and Corollary 4.1 have been checked, while its four experiments have their subject and qualitative outcomes only. Paper **54** has now been verified in full from its preprint (arXiv:2104.06291): the setting, the complete Schur-complement and pivoted-Cholesky derivation, and every experimental setup (candidate pool $M=10^4$, $Q=1000$ test points, 50 trials, three kernels, four classes of test problem) are taken from the source. But that paper **contains no theorems**, and all of its numerical results appear as figures with no data tables, so no condition numbers or error values are quoted here.
>
> Six source-level problems are preserved on this page: paper 22's three limitations on Theorem 4.3 and the non-vanishing $4\kappa^2(R)d^2(f)$ term in its main theorem; the inconsistent Gaussian normalising constant between paper 45's Example 2.2 and Example 8.1; the misprinted starting condition in paper 28's Theorem 3.2; the mismatch between line 4 of paper 24's Algorithm 1 and its eq. (9); the internal inconsistency in paper 24's CSA-c branch expansion; and the missing new point $z$ in the greedy formula of paper 54's §3.1. The asymptotic Gaussian induced measure is treated as a **conjecture** everywhere it appears.

## Coverage check

| Item                                               | Paper | Status                                                                         |
| -------------------------------------------------- | ----- | ------------------------------------------------------------------------------ |
| Stability factor and its basis independence        | 22    | $K(z)=\varphi^T\varphi$, $\lVert K\rVert_\infty/N$, invariance                 |
| Decoupling sampling from orthogonality             | 22    | the insight, the $\tilde K_k\equiv N$ cancellation, the definition             |
| Five algorithm steps and row normalisation         | 22    | full steps and the equivalent reading                                          |
| Christoffel asymptotics and equilibrium measure    | 22    | Theorem 4.2, Corollary 4.1, Theorem 4.3 and its three limitations              |
| CLS stability and accuracy theorems                | 22    | Theorems 5.1 and 5.2, and the non-vanishing $4\kappa^2(R)d^2(f)$               |
| Experimental design and qualitative outcomes       | 22    | Figure 1, §6.1.1-6.1.3, §6.2.1-6.2.3; numbers not verified                     |
| Sample complexity for three sampling routes        | 45    | quadratic, $N^{\log3/\log2}$, deterministic quadratic                          |
| Matrix Chernoff proof and origin of the constant   | 45    | the full argument for Theorem 6.1 and $C=2/\log(27/8e)$                        |
| Induced distribution and the optimality argument   | 45    | Definition 7.1, attainment of the lower bound $N$, $M/\log M\ge C(r+1)N$       |
| Accuracy with truncation                           | 45    | Theorem 7.2 and Lemma 5.1                                                      |
| Three numerical experiments                        | 45    | full setups and outcomes for Example 5.2, Example 8.1, the Section 9 PDE       |
| Asymptotic induced measures                        | 45    | univariate and multivariate Chebyshev limits; Gaussian as conjecture           |
| Weighted approximate Fekete points, pivoted QR     | 28    | selection mechanism, meaning of weighting, enrichment and candidate pool       |
| Row normalisation and equivalence of objectives    | 28    | Theorem 3.1, its existence premise, and this page's reading of it              |
| Greedy is optimal in one dimension                 | 28    | Lemma 3.1, Theorem 3.2, Corollary 3.1 and the misprint                         |
| Three-way comparison experiments                   | 28    | $\tilde M=10^4$, $M=1.05N$, 50 trials, four setups and outcomes                |
| Sampling and preconditioning designed as a pair    | 24    | the reciprocal structure, Algorithm 1, three regimes, the Gramian $R$          |
| Main recovery theorem, bounded versus unbounded    | 24    | Theorem 4.1, the three $L(n)$ cases, Corollary 4.1, two inconsistencies        |
| Univariate theory against multivariate experiments | 24    | Remark 4.3's exponential dependence and four qualitative outcomes              |
| From polynomial spaces to kernel spaces            | 54    | power function replaces Christoffel; Schur-complement greedy, pivoted Cholesky |
| Structural identity with paper 28                  | 54    | both greedy over a candidate pool; QR pivots versus Cholesky pivots            |
| Four experiment classes, the sparse-grid claim     | 54    | $M=10^4$, $Q=1000$, 50 trials; sparse grids compared exactly once              |
| The source proves nothing                          | 54    | no theorem/lemma/proof environment; "quasi-optimal" carries no bound           |

## Sources for this page

- A. Narayan, J. Jakeman, and T. Zhou, [_A Christoffel function weighted least squares algorithm for collocation approximations_](https://doi.org/10.1090/mcom/3192), Math. Comput. 86 (2017), pp. 1913-1947.
- J. Jakeman, A. Narayan, and T. Zhou, [_A generalized sampling and preconditioning scheme for sparse approximation of polynomial chaos expansions_](https://doi.org/10.1137/16M1063885), SIAM J. Sci. Comput. 39(3) (2017), pp. A1114-A1144.
- L. Guo, A. Narayan, L. Yan, and T. Zhou, [_Weighted approximate Fekete points: sampling for least-squares polynomial approximation_](https://doi.org/10.1137/17M1140960), SIAM J. Sci. Comput. 40(1) (2018), pp. A366-A387.
- L. Guo, A. Narayan, and T. Zhou, [_Constructing least-squares polynomial approximations_](https://doi.org/10.1137/18M1234151), SIAM Rev. 62(2) (2020), pp. 483-508.
- A. Narayan, L. Yan, and T. Zhou, [_Optimal design for kernel interpolation: applications to uncertainty quantification_](https://doi.org/10.1016/j.jcp.2020.110094), J. Comput. Phys. 430 (2021), 110094.
