---
title: Stability Theory for Discretisations
description: Papers 47 and 63 - a mean-square Lax equivalence theorem, and the first fully discrete error analysis
lang: en
translation: computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes
tags:
  - paper-notes
  - stochastic-differential-equations
  - stability-theory
---

> [!note] Coverage of this page
> Papers **47** (_SIAM J. Numer. Anal._ 58(4), 2020) and **63** (_SIAM J. Numer. Anal._ 60(4), 2022).
>
> Neither paper has a preprint, and the SIAM full texts are unreachable from this environment (`epubs.siam.org` is behind a Cloudflare challenge), so **neither article body was read here**. The method of this page is therefore: state verbatim what the abstracts confirm; verify the **independent external sources** that supply the technical components separately and attribute them (the Chessari-Kawai-Shinozaki-Yamada survey of BSDE numerics, the Sun-Zhao CSIAM-AM article on Sinc quadrature, the Wang-Zhao Sinc multistep sequel, and the stability papers of Tang-Xiong and Chassagneux); and mark everything else as reconstructed from context. **No numerical data from either paper was obtained here, so this page carries no experiment tables for them.**

## 47: replacing scheme-by-scheme proofs with one equivalence theorem

### The idea

By 2020 this group's collection of temporal discretisations for FBSDEs had grown into a small zoo: backward Euler, the $\theta$-scheme family, the Zhao-Zhang-Ju multistep scheme, the multistep family of [[en/computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|paper 8]], deferred correction. Each carried its own convergence proof and its own stability condition — paper 8's $1\le k\le6$ is a numerically computed root-condition window, paper 23's "no barrier" is only a numerical assertion, and paper 33 proves stability first and then derives error estimates.

That situation is itself a signal: **when the convergence proofs of a family of methods are ad hoc one at a time, it means nobody has yet isolated and defined "stability" on its own.** Deterministic numerical analysis has had the template for a long time — the Lax equivalence theorem for linear evolution equations says that a consistent scheme converges if and only if it is stable. Its value is not that proofs get shorter but that it swaps the vague target "convergence" for the property "stability", which can be checked scheme by scheme and even optimised.

This paper carries that template over to FBSDEs. Two things make the transfer non-trivial: the norm must be mean-square ($L^2(\Omega)$) because the solution is a stochastic process, and the problem is nonlinear (the generator is only assumed Lipschitz) while the classical Lax theorem is about linear equations.

### Setting

Consider the decoupled Markovian backward equation

$$
Y_t=\varphi(X_T)+\int_t^Tf(s,X_s,Y_s,Z_s)\,\mathrm ds-\int_t^TZ_s\,\mathrm dW_s .
$$

The abstract is explicit: what is proposed is "a **general discretisation framework** for the numerical solution of FBSDEs", and that framework "covers several existing probabilistic temporal discretisation schemes in the literature". The description in Chessari, Kawai, Shinozaki and Yamada's survey of BSDE numerics agrees independently and is more specific about the scope: "a general framework is constructed in [311] to study stability, consistency and convergence of the discretisation schemes for FBSDEs in a unified way, **including the backward Euler method, theta-schemes and various multistep methods**."

### Derivation

**Step one: what the two unified families look like.** One is the $\theta$-scheme. The two recursions below are transcribed verbatim from the survey, which attributes them to Zhao, Wang and Peng, and to Zhao, Li and Zhang: for $k=n-1,\dots,0$,

$$
Y^\pi_{t_k}=\mathbb E\bigl[Y^\pi_{t_{k+1}}\mid\mathcal F_{t_k}\bigr]
+\theta_1\Delta_n f\bigl(t_k,Y^\pi_{t_k},Z^\pi_{t_k}\bigr)
+(1-\theta_1)\Delta_n\,\mathbb E\bigl[f(t_{k+1},Y^\pi_{t_{k+1}},Z^\pi_{t_{k+1}})\mid\mathcal F_{t_k}\bigr],
$$

$$
Z^\pi_{t_k}=-\frac{1-\theta_2}{\theta_2}\mathbb E\bigl[Z^\pi_{t_{k+1}}\mid\mathcal F_{t_k}\bigr]
+\frac{1}{\theta_2\Delta_n}\mathbb E\bigl[Y^\pi_{t_{k+1}}\Delta W_k\mid\mathcal F_{t_k}\bigr]
+\frac{1-\theta_2}{\theta_2}\mathbb E\bigl[f(t_{k+1},Y^\pi_{t_{k+1}},Z^\pi_{t_{k+1}})\Delta W_k\mid\mathcal F_{t_k}\bigr],
$$

with parameters $\theta_1\in[0,1]$, $\theta_2\in(0,1]$ and terminal values $Y^\pi_{t_n}=\Phi(X^\pi_{t_n})$, $Z^\pi_{t_n}=(\nabla\Phi(X^\pi_{t_n}))\sigma(t_n,X^\pi_{t_n})$. It comes from the exact identity $Y_{t_k}=\mathbb E[Y_{t_{k+1}}|\mathcal F_{t_k}]+\int_{t_k}^{t_{k+1}}\mathbb E[f(s,Y_s,Z_s)|\mathcal F_{t_k}]\mathrm ds$, with the integral approximated by a convex combination of the implicit ($t_k$) and explicit ($t_{k+1}$) endpoint values.

**One notational point is easy to read backwards: $\theta_1$ multiplies the generator at the current level $t_k$, and in a backward recursion the current level is the unknown one.** Hence:

| $(\theta_1,\theta_2)$ | Name and character                                                                                       |
| --------------------- | -------------------------------------------------------------------------------------------------------- |
| $\theta_1=1$          | left-endpoint rectangle rule, generator implicit (backward-Euler-like)                                   |
| $(1,1)$               | as above, with the $Z$ recursion collapsing to the single term $\mathbb E_n[Y_{n+1}\Delta W_n]/\Delta_n$ |
| $\theta_1=0$          | right-endpoint rule, fully explicit                                                                      |
| $(1/2,1/2)$           | the Crank-Nicolson member, the only second-order one                                                     |

The $(1,1)$ member is exactly the scheme paper 26 applies to its adjoint BSDE (see the [[en/computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control|second-order FBSDEs and control page]]); it is nominally implicit yet needs no iteration, because the adjoint generator there is linear in $p$. The other family is the multistep scheme of paper 8, whose form is on the [[en/computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|multistep schemes page]].

**Step two: the shape of the unified scheme.** A $K$-step probabilistic discretisation able to hold both families would be a pair of recursions roughly of the form

$$
\sum_{j=0}^{K}a_j\,\mathbb E_{t_n}\bigl[Y^{n+j}\bigr]
=\Delta t\sum_{j=0}^{K}b_j\,\mathbb E_{t_n}\bigl[f(t_{n+j},X^{n+j},Y^{n+j},Z^{n+j})\bigr],
$$

$$
\sum_{j=0}^{K}c_j\,\mathbb E_{t_n}\bigl[Z^{n+j}\bigr]
=\frac{1}{\Delta t}\sum_{j=0}^{K}d_j\,\mathbb E_{t_n}\bigl[Y^{n+j}\Delta W_{n,j}\bigr]
+\sum_{j=0}^{K}e_j\,\mathbb E_{t_n}\bigl[f(t_{n+j},\cdots)\Delta W_{n,j}\bigr],
$$

with coefficient vectors $(a_j),(b_j),(c_j),(d_j),(e_j)$. **This pair is a shape inferred here from the two families it must cover, not transcribed from the paper**; the index ranges, normalisations and notation the authors actually use were not verified here. What is verifiable is the **coverage claim** itself: the framework covers backward Euler, $\theta$-schemes and various multistep methods.

**Step three: what the three concepts look like in the mean-square sense.** Write $\|\cdot\|=\|\cdot\|_{L^2(\Omega)}$. The three definitions below are reconstructed here from the phrase "mean-square version of the Lax equivalence theorem", so **do not cite constants or exact norms from them**:

- **Consistency of order $p$.** Substituting the **exact solution** $(Y_{t_n},Z_{t_n})$ into the discrete operator leaves local truncation residuals $(R^Y_n,R^Z_n)$ with
  $$
  \Bigl(\sum_n\|R^Y_n\|^2+\Delta t\sum_n\|R^Z_n\|^2\Bigr)^{1/2}=O\bigl((\Delta t)^p\bigr)\xrightarrow[\Delta t\to0]{}0 .
  $$
  **The extra factor of $\Delta t$ on the $Z$ term is standard in BSDE analysis**; it reflects that $Z$ has "half a derivative" less regularity than $Y$.
- **Stability.** Perturbing the discrete data by $(\delta^Y_n,\delta^Z_n)$ (including, for a $k$-step scheme, the starting values it requires), the perturbation of the discrete solution is controlled by a constant **independent of $\Delta t$** times the accumulated perturbation, uniformly in $N$:
  $$
  \max_n\|\delta Y^n\|^2+\Delta t\sum_n\|\delta Z^n\|^2\le C\bigl(\|\delta Y^N\|^2+\text{accumulated data perturbation}\bigr).
  $$
- **Convergence.** $\max_n\|Y_{t_n}-Y^n\|\to0$ together with the corresponding $Z$ quantity as $\Delta t\to0$; if the rate is $O((\Delta t)^p)$ the scheme converges with order $p$.

### Theorems

**Main theorem (content verifiable, exact statement not verified): a consistent FBSDE discretisation scheme converges if and only if it is stable.** The abstract's own wording is that "in particular, we prove a **stochastic mean-square version of the Lax equivalence theorem** — showing that a consistent discretisation scheme for FBSDEs is convergent if and only if it is stable." This corresponds exactly to the classical Lax equivalence theorem for well-posed linear evolution equations, the differences being that the norm is mean-square and the problem is nonlinear (Lipschitz generator).

The abstract adds that "applications of the analysis to existing numerical schemes are also discussed", that is, the abstract theorem is instantiated to recover or strengthen the known convergence results for backward Euler, $\theta$-schemes and multistep schemes. **This is precisely where paper 8's empirical window acquires its theoretical home**: the root condition there on $P(\lambda)=\alpha_{k,0}\lambda^k+\sum_{j=1}^k\lambda^{k-j}$ is the stability half of the equivalence, consistency comes from the construction, and together they give convergence.

> [!note] What "stable" means concretely in this family — three independently verified corroborations
> The article body was not read here, but three external works from the same period pin the concept down from the outside, and all three are interlocutors of this paper.
>
> **First, the root condition.** Applying paper 8's weights to a deterministic ODE gives $\alpha_{k,0}Y^n+\sum_{j=1}^k\alpha_{k,j}Y^{n+j}=f(t_n,Y^n)$, whose characteristic polynomial $P(\lambda)=\alpha_{k,0}\lambda^k+\sum_{j=1}^k\lambda^{k-j}=0$ must satisfy the classical root condition ($|\lambda|\le1$, with roots of modulus one simple). The largest root moduli it reports (excluding the common root $1.0$) are $k=2$: 0.3333, 3: 0.4264, 4: 0.5608, 5: 0.7087, 6: 0.8633, 7: **1.0222**, 8: **1.1839**, so $k\ge7$ is unstable.
>
> **Second, the root condition really does imply mean-square stability, and this is proved by external work.** The abstract of Tang and Xiong (_IMA J. Numer. Anal._ 42(2) (2022) 1789-1805) says verbatim: "Under the classical root condition, we prove that the general linear multistep methods for decoupled FBSDEs, whose generators depend on both $y$ and $z$, are **mean-square (zero) stable**. Based on this stability result, we further establish a fundamental convergence theorem." **"Mean-square (zero) stable" is almost certainly the concept paper 47 formalises.**
>
> **Third, the same architecture appears elsewhere.** Chassagneux (_SIAM J. Numer. Anal._ 52(6) (2014) 2815-2836) states verbatim: "We prove that these schemes enjoy a **fundamental stability property** under a sufficient condition on the coefficients. Combining this result with a truncation error analysis allows us to design approximations with arbitrary order of convergence." That is the same "stability plus truncation error yields order" architecture, which paper 47 axiomatises.
>
> Two already-recorded concrete stability windows can serve as calibration: Zhao-Zhang-Ju's $Y$ equation is stable only for $K_y\in\{1,\dots,7,9\}$ and its $Z$ equation only for $K_z\in\{1,2,3\}$; papers 8 and 19 both give $1\le k\le6$.

> [!note] What could be verified
> The **content** of the main theorem (consistency plus stability is equivalent to convergence) and the coverage claim (backward Euler, $\theta$-schemes, various multistep schemes) are confirmable from the abstract and independently corroborated by the Chessari et al. survey. The exact form of the general discretisation family, the precise norms in the definitions of stability and consistency, the explicit constants, the precise hypotheses (very likely uniform Lipschitz continuity of $f$ in $(y,z)$, square-integrable terminal data, and regularity of the forward diffusion), and whether the framework covers **coupled** rather than only decoupled systems, are all unverified here. The $\theta$-scheme above is transcribed in its standard form from an independent third-party survey, and the shape of the unified scheme is inferred here.

### Numerical experiments

**The abstract mentions no numerical experiments at all**, only "applications of the analysis to existing numerical schemes". This page therefore reports no experimental data for this paper. The abstract suggests a purely analytical paper, but that is unverified here — the body was not read, and the presence of examples inside it cannot be ruled out.

### Relation to the others

**This is the theoretical keystone of the whole series.** Papers 8 and 19 propose multistep schemes and give truncation error estimates but no convergence theorem; papers 23 and 35 propose deferred correction schemes, again with no FBSDE convergence theorem; paper 33 proves stability first for one specific family and then derives error estimates. Paper 47 abstracts this recurring pattern into a framework and proves that stability plus consistency is equivalent to convergence. It directly explains the $k\le6$ barrier observed empirically in papers 8 and 19 and asserted numerically in paper 25.

**It also makes paper 68's reverse design possible**: once stability is established as the pivotal property, "design a convergent scheme" becomes "design a stable scheme", and the latter is an optimisation problem with a clear objective — which is exactly how paper 68 obtains new sufficient conditions on the coefficients and optimises strong-stability-preserving schemes of orders one to five. Paper 63 pushes the analysis from **semi-discrete** (time only, conditional expectations assumed exact) to **fully discrete** (time plus space), and the two share the "stability first" architecture. The first author, Jie Yang, is the same as in papers 35 and 61; the earlier Yang-Zhao paper _Convergence of recent multistep schemes for a forward-backward stochastic differential equation_ (EAJAM 5 (2015) 387-404, found in paper 41's reference list) is its predecessor.

## 63: the first fully discrete error analysis

### The idea

Every probabilistic BSDE scheme has two error sources: the **temporal** discretisation of the backward equation, and the **spatial** approximation of the conditional expectations $\mathbb E_{t_n}[\cdot]$, which in the Markovian case are $d$-dimensional Gaussian integrals whose integrand is known only at grid points.

This group's earlier pipeline (papers 8, 19, 25, 41) computed those integrals with Gauss-Hermite quadrature. Its nodes are $x_n+\sqrt{\Delta t}\,\lambda_j$ with $\lambda_j$ the roots of a Hermite polynomial — **generally irrational, and so never on the spatial grid**, which forces a local polynomial interpolation at every node at every time step. That step costs three things: it takes time; it caps spatial accuracy at the interpolation order; and it makes a **fully discrete** error analysis extremely difficult. The consequence is that before this paper, essentially every rigorous BSDE error analysis was **semi-discrete**: conditional expectations assumed exact, only the temporal error bounded.

**The mechanism of this paper fits in one sentence: change the quadrature rule so that its nodes land on the grid by themselves.** Sinc quadrature has equispaced nodes $kh$, and $h$ is a **free parameter**. So tuning $h$ until the node spacing equals the spatial step makes the interpolation step disappear entirely. As a bonus, Sinc quadrature is **exponentially convergent** for a suitable class of functions, far above any interpolation order.

### Setting

A standard (first-order) BSDE, with **no** $\Gamma$ process — this paper is not about 2FBSDEs. The paper's own framing is that this seems to be the **first attempt at analysing fully discrete schemes for BSDEs**, achieving second-order convergence in time and exponential convergence in space. The zbMATH keywords are "backward stochastic differential equations; error estimates; conditional mathematical expectation; Sinc-$\theta$ schemes", with MSC codes 60H10, 60H35, 65C30.

### Derivation

**Step one: the $\theta$-scheme family in time.** The abstract confirms verbatim that "we consider the $\theta$-scheme for the temporal discretisation". The two recursions of that family are the same as those quoted for paper 47 above (again transcribed from the Chessari et al. survey). The survey also records that the scheme converges with order two when $\theta_1=\theta_2=1/2$ and order one otherwise — **which matches this paper's claimed second-order rate in time exactly, meaning the Crank-Nicolson member is the one used.**

There is a detail about the order worth recording. Zhao, Wang and Peng (_DCDS-B_ 12(4) (2009) 905-924) studied a **single-parameter** $\theta$-scheme with a generator **not depending on $z$**, and their abstract says: order one for $y$ for general $\theta$; order two for $y$ but **order one for $z$** when $\theta=\frac12$. Reaching second order for $Y$ **and** $Z$ simultaneously when the generator depends on $(y,z)$ requires the generalisation of Zhao, Li and Zhang (_DCDS-B_ 17(5) (2012) 1585-1603), which "introduces more parameters". **Both papers are in paper 63's citation list**, which is why its "second order in time" is consistent with $\theta_1=\theta_2=1/2$.

**Step two: Sinc approximation in space.** The abstract confirms verbatim that "the Sinc approximation is then adopted to approximate the associated conditional mathematical expectations". The construction below is transcribed from another article by **the same group** (Y. Sun and W. Zhao, _CSIAM Trans. Appl. Math._ 6(1) (2025) 176-206, freely downloadable, its Section 4.1):

$$
\mathrm{sinc}(x)=
\begin{cases}
\dfrac{\sin(\pi x)}{\pi x},&x\neq0,\\[6pt]
1,&x=0 .
\end{cases}
$$

Let $B(h)$ be the class of entire functions $g$ with $g\in L^2(\mathbb R)$ on the real axis and $|g(z)|\le K\exp(\pi|z|/h)$ for all $z\in\mathbb C$. **Whittaker cardinal expansion (Stenger)**: if $g\in B(h)$ then $g(z)=\sum_{k=-\infty}^{\infty}g(kh)\,\mathrm{sinc}\bigl(\frac{z-kh}{h}\bigr)$; and if $\sum_kg(kh)$ converges then $\int_{\mathbb R}g(x)\mathrm dx=h\sum_{k=-\infty}^{\infty}g(kh)$ for sufficiently small $h$. This defines the **Sinc quadrature rule**

$$
T_M(g,h)=h\sum_{k=-M}^{M}g(kh),
\qquad
\eta_M(g,h)=\int_{\mathbb R}g(x)\,\mathrm dx-T_M(g,h),
$$

whose error theorem reads: if $g$ is bounded then, for sufficiently small $h$ and any $\gamma_0>0$ with $\gamma_0\le Mh^2$,

$$
|\eta_M(g,h)|\ \le\ C_{\gamma_0,g}\,h\,\exp\Bigl(-\frac{M^2h^2}{2}\Bigr),
$$

with $C_{\gamma_0,g}$ depending only on $\gamma_0$ and $\|g\|_{L^\infty}$. **This is where "exponential convergence in space" comes from.**

**Step three: applying it to Gaussian conditional expectations.** For $X^{t,x}_r=x+\sigma_0(W_r-W_t)$ and smooth $v:\mathbb R^d\to\mathbb R$,

$$
\mathbb E\bigl[v(X^{t,x}_r)\bigr]
=\int_{\mathbb R^d}v\bigl(x+\sigma_0\sqrt{r-t}\,p\bigr)\Bigl(\tfrac{1}{\sqrt{2\pi}}\Bigr)^{d}e^{-p^\top p/2}\,\mathrm dp
\ \approx\ \sum_{\mathbf k=-M}^{M}v\bigl(x+\sigma_0\sqrt{r-t}\,h\mathbf k\bigr)\,\alpha_{\mathbf k},
$$

$$
\alpha_{\mathbf k}=\prod_{i=1}^{d}\alpha_{k_i},
\qquad
\alpha_{k_i}=\frac{h}{\sqrt{2\pi}}\exp\Bigl(-\frac{k_i^2h^2}{2}\Bigr),
\qquad
\beta^M_{\mathbf k}=\frac{\alpha_{\mathbf k}}{\sum_{\mathbf k=-M}^{M}\alpha_{\mathbf k}} .
$$

The $\beta^M_{\mathbf k}$ are the **renormalised** weights, which sum to one — making the row sums of the resulting matrix equal to one, and so preserving boundedness. **(The display above is transcribed from Sun and Zhao; the corresponding step inside paper 63 was not verified here and is inferred by analogy.)**

**Step four: why Sinc removes the interpolation.** The abstract confirms verbatim that "by choosing appropriate parameters in the Sinc quadrature rule, our scheme is shown to be very efficient, as **no spatial interpolation is needed**". The mechanism is stated most plainly in Sun and Zhao's Remark 4.2, transcribed verbatim:

> "Except for the spectral accuracy, the main reason we choose the Sinc quadrature rule rather than other quadrature rules is that its quadrature nodes are **uniform** and there is a **free parameter $h$** in it. Thus by setting different values of $h$ for different time increments between different intermediate time levels, we can obtain different uniform quadrature nodes, and consequently **avoid using interpolation** when approximating conditional expectations."

Concretely: the Sinc nodes for the increment $r-t$ sit at $x+\sigma_0\sqrt{r-t}\,h\mathbf k$, an **arithmetic progression** of spacing $\sigma_0\sqrt{r-t}\,h$. Choosing $h$ so that $\sigma_0\sqrt{r-t}\,h=\Delta x$ (the spatial step) puts every quadrature node **exactly on a grid point**. **Gauss-Hermite cannot do this, because its nodes are the roots of Hermite polynomials.** The same authors' Sinc multistep sequel (X. Wang and W. Zhao, _Adv. Appl. Math. Mech._ 15(3) (2023)) states the same mechanism more explicitly: "by using the **change of integral variables** and **choosing the space step parameter properly** in the Sinc quadrature rule, no spatial interpolation is needed."

**Step five: what it amounts to together.** Substituting the Sinc rule into the four conditional expectations of the $\theta$-scheme turns each into a finite, grid-aligned weighted sum, so the scheme becomes a **matrix recursion** on the grid values $\{Y^n_{\mathbf i}\},\{Z^n_{\mathbf i}\}$, with banded/Toeplitz transition matrices whose row sums are one. **This paragraph is a reconstruction here**: whether the paper writes it in matrix form, how the truncation of the spatial domain is handled (truncating the Sinc sum at $\pm M$ is itself a domain truncation), and the exact relations imposed among $h,M,\Delta t,\Delta x$ were all not verified here.

The treatment of $Z$ deserves its own note (also reconstructed): as the $\theta$-scheme above shows, $Z^n$ comes from $\mathbb E_n[Y^{n+1}\Delta W_n]/(\theta_2\Delta t)$ plus correction terms in $\mathbb E_n[Z^{n+1}]$ and $\mathbb E_n[f^{n+1}\Delta W_n]$, so **no derivative of an interpolant is ever taken**. That matters, because differentiating an interpolant is the usual bottleneck for the accuracy of $Z$ elsewhere — it is exactly the price deferred correction (paper 23) pays.

### Theorems

- **Main result (verbatim from the abstract):** "**Rigorous stability analysis and error estimates** are carried out, which seem to be the first attempt at analysing fully discrete schemes for BSDEs, achieving a **second-order convergence rate in time** and **exponential convergence in space**."
- The headline estimate therefore has roughly the shape
  $$
  \max_n\bigl\|Y_{t_n}-Y^{n,M,h}\bigr\|\ \le\ C\Bigl(\underbrace{(\Delta t)^2}_{\theta_1=\theta_2=1/2}+\underbrace{C'\,h\,e^{-M^2h^2/2}}_{\text{Sinc}}\Bigr),
  $$
  possibly with a domain-truncation term as well. **This shape is assembled here from two verified rates rather than transcribed**; the exact statement, norms, constants, required regularity, whether the same rate holds for $Z$, and the precise coupling conditions among $M,h,\Delta t$ are all unverified here.
- **The stability analysis certainly was carried out** (confirmable from the abstract), but the precise meaning of "stable" here was not verified. Given that the paper cites paper 47, it is almost certainly the mean-square stability notion of that framework, now extended to include the **spatial perturbation introduced by the quadrature** — and that extension is exactly where the technical novelty lies.

> [!note] An independently verified sibling result, usable as calibration
> The abstract of the Sinc **multistep** sequel (X. Wang and W. Zhao, AAMM 2023) says verbatim: "The stability and the $K$-order error estimates in time of the $K$-step Sinc multistep scheme are theoretically proved ($1\le K\le6$). This seems to be the first result of this kind for fully discrete multistep schemes for FBSDEs in space and time." Its contribution list adds that the optimal $K$-order holds "for FBSDEs (1.1) with the drift term $b=0$ and the diffusion term $\sigma$ a non-singular matrix", obtained "by using **numerical algebra theory** and properties of the approximation operators". **The $1\le K\le6$ here is exactly paper 8's root-condition window** — the same barrier reappearing in a fully discrete analysis, which shows it comes from the temporal template rather than from the spatial approximation.

> [!note] What could be verified
> "First fully discrete analysis", "second order in time and exponential in space", "the $\theta$-scheme in time with Sinc approximation for the conditional expectations", and "no spatial interpolation needed" are all confirmable verbatim from the abstract. The concrete form of the $\theta$-scheme is transcribed from a third-party survey; the definition of Sinc quadrature, its error theorem and the form used for Gaussian conditional expectations are transcribed from the same group's CSIAM-AM article; the "uniform nodes plus a free parameter, hence no interpolation" mechanism is transcribed from that article's Remark 4.2 and from the Sinc multistep sequel's abstract. **The body of this paper itself was not read here**, so the exact statements, constants and hypotheses of its theorems, and the concrete form of the fully discrete scheme, are not reported.

### Numerical experiments

The abstract says "several numerical experiments are provided to verify the theoretical results and to show the efficiency and accuracy of our scheme". **The concrete test problems and observed orders were not verified here, so no table is given.** What the theoretical results would lead one to expect is: a temporal refinement table showing a rate near 2 for $\theta_1=\theta_2=1/2$ and near 1 otherwise; a semi-log plot in $M$ or $h$ showing the spatial error decaying along a straight line (exponentially); and a CPU-time comparison against a "Gauss-Hermite plus interpolation" scheme. **None of these three was verified here, and none should be taken as a statement of fact.**

### Relation to the others

**This is the spatial-discretisation counterpart of the whole research programme.** Papers 8, 19, 25 and 41 all approximate the conditional expectations with Gauss-Hermite quadrature plus local Lagrange interpolation; paper 63 replaces both with a single Sinc rule and makes the nodes coincide with the grid. It is the first paper on this list to prove a **fully discrete** (time plus space) error estimate.

It cites paper 47 and inherits its "stability first" architecture, and it also cites papers 8, 18, 25 and 33, placing itself opposite the entire multistep/$\theta$ thread. Paper 68 cites it (as its reference 39), and the two 2022-2023 articles are companions on the stability theme. Paper 93 (the deep random difference method) also cites it, showing that the classical Sinc work fed back into the deep-learning branch. Its direct sequel is Wang and Zhao's _Sinc-Multistep Schemes for FBSDEs_ (AAMM 2023), which combines the Sinc spatial rule with paper 8-style multistep temporal rules and proves order $K$ for $1\le K\le6$.

**A spillover worth recording:** the same Sinc machinery was later used by Sun and Zhao for **maximum-bound-principle-preserving** stochastic Runge-Kutta schemes for the Allen-Cahn equation (CSIAM-AM 2024/2025, parts I and II), a crossover from BSDE numerics into structure-preserving parabolic PDE solvers. Tao Zhou is also a co-author of the BDF2/Allen-Cahn article cited in that thread (H. Liao, T. Tang and T. Zhou, _SIAM J. Numer. Anal._ 58 (2020) 2294-2314).

## How the two relate

| No. | Error sources covered | Main result                                                | Effect on the other papers                                                          | Verification here              |
| --- | --------------------- | ---------------------------------------------------------- | ----------------------------------------------------------------------------------- | ------------------------------ |
| 47  | temporal              | mean-square Lax equivalence                                | gives paper 8's empirical window a theory; makes paper 68's reverse design possible | abstract + independent survey  |
| 63  | temporal and spatial  | first fully discrete analysis (order two plus exponential) | a second route of spatial improvement alongside paper 25                            | abstract + components verified |

One general judgement: **when a family of methods has scheme-by-scheme convergence proofs, it is worth stopping to define stability properly.** After paper 47, "design a convergent scheme" becomes "design a stable scheme", and the latter is an objective one can optimise — which is exactly what [[en/computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|paper 68]] does.

A second judgement, about paper 63: **what it changes is not the accuracy but the location of the nodes.** Gauss-Hermite is not inaccurate; the problem is that its nodes are determined by the roots of a polynomial and cannot be aligned with a grid, which forces an interpolation step that consumes time, accuracy and analysability at once. Swapping in a quadrature rule with equispaced nodes and adjustable spacing removes all three costs together. **This is an improvement by replacing a tool rather than by adding one, which is rarer in numerical analysis than it looks.**

## Coverage checklist

| Item                                                                                | Paper | Coverage status                                             |
| ----------------------------------------------------------------------------------- | ----- | ----------------------------------------------------------- |
| The two $\theta$-scheme recursions and four named members                           | 47    | transcribed from a third-party survey                       |
| The shape of the unified scheme                                                     | 47    | inferred here, flagged as such                              |
| Mean-square definitions of consistency, stability, convergence                      | 47    | reconstructed here, flagged as such                         |
| The mean-square Lax equivalence theorem and the coverage claim                      | 47    | verbatim from the abstract, corroborated by survey          |
| Three external corroborations of "stable" (root condition, Tang-Xiong, Chassagneux) | 47    | each independently verified                                 |
| No numerical experiments                                                            | 47    | not mentioned in the abstract; noted                        |
| The order of the $\theta$-scheme and the roles of the two DCDS-B papers             | 63    | independently verified                                      |
| Definition of Sinc quadrature, its error theorem, Gaussian form                     | 63    | transcribed from the group's CSIAM-AM article               |
| The "uniform nodes plus free parameter, hence no interpolation" mechanism           | 63    | verbatim from Remark 4.2 and the sequel's abstract          |
| The shape of the fully discrete estimate and the stability analysis                 | 63    | two rates verifiable; shape assembled here, flagged         |
| The sequel's $1\le K\le6$ and its $b=0$ condition                                   | 63    | independently verified                                      |
| Numerical experiments                                                               | 63    | existence stated in the abstract; data unverified, no table |

## Sources for this page

- J. Yang, W. Zhao, and T. Zhou, [_A unified probabilistic discretization scheme for FBSDEs: stability, consistency, and convergence analysis_](https://doi.org/10.1137/19M1260177), SIAM J. Numer. Anal. 58(4) (2020), pp. 2351-2375.
- X. Wang, W. Zhao, and T. Zhou, [_Sinc-theta schemes for backward stochastic differential equations_](https://doi.org/10.1137/21M1444679), SIAM J. Numer. Anal. 60(4) (2022), pp. 1799-1823.
- External sources used for cross-checking: C. Chessari, R. Kawai, Y. Shinozaki, and T. Yamada, [_Numerical methods for backward stochastic differential equations: a survey_](https://doi.org/10.1214/23-PS18), Probab. Surveys 20 (2023), pp. 486-567 ([arXiv:2101.08936](https://arxiv.org/abs/2101.08936); source of the $\theta$-scheme transcription and of the independent description of paper 47); Y. Sun and W. Zhao, _Stochastic Runge-Kutta methods for preserving maximum bound principle of semilinear parabolic equations, part II: Sinc quadrature rule_, [CSIAM Trans. Appl. Math. 6(1) (2025), pp. 176-206](https://doi.org/10.4208/csiam-am.SO-2024-0012) (source of the Sinc quadrature definition, error theorem and interpolation-free mechanism); X. Wang and W. Zhao, [_Sinc-multistep schemes for forward backward stochastic differential equations_](https://doi.org/10.4208/aamm.OA-2022-0073), Adv. Appl. Math. Mech. 15(3) (2023) (the multistep sequel to paper 63); X. Tang and J. Xiong, [_Stability analysis of general multistep methods for Markovian backward stochastic differential equations_](https://doi.org/10.1093/imanum/drab023), IMA J. Numer. Anal. 42(2) (2022), pp. 1789-1805; J.-F. Chassagneux, [_Linear multistep schemes for BSDEs_](https://doi.org/10.1137/120902951), SIAM J. Numer. Anal. 52(6) (2014), pp. 2815-2836; W. Zhao, J. Wang, and S. Peng, [_Error estimates of the theta-scheme for backward stochastic differential equations_](https://doi.org/10.3934/dcdsb.2009.12.905), Discrete Contin. Dyn. Syst. Ser. B 12(4) (2009), pp. 905-924; W. Zhao, Y. Li, and G. Zhang, _A generalized theta-scheme for solving backward stochastic differential equations_, Discrete Contin. Dyn. Syst. Ser. B 17(5) (2012), pp. 1585-1603.
