---
title: Unbounded Domains and Fractional Operators
description: Papers 27, 42 and 92 - making the far-field behaviour of the basis compatible with the action of a fractional operator
lang: en
translation: computational-mathematics/paper-notes/spectral-and-reduced-order/unbounded-domain-spectral
tags:
  - paper-notes
  - spectral-methods
  - fractional-operators
---

> [!note] Coverage of this page
> Papers **27** (_Commun. Comput. Phys._ 24(4), 2018), **42** (_SIAM J. Sci. Comput._ 42(2), 2020) and **92** (_SIAM J. Sci. Comput._ 48(1), 2026). Papers 27 and 42 were checked equation by equation against the preprints and the published texts, so this page gives their full derivations, the complete hypotheses of their theorems, and their numerical setups; the error magnitudes of both papers live mostly in figures and rate tables, and wherever a number could not be recovered from the available text this page says so instead of filling the gap. The full text of paper 92 cannot be reached through public channels — the publisher's site returns a Cloudflare challenge to every client available here, and there is no preprint and no open-access copy — so that section reports only what the abstract, keywords and reference list confirm and **gives no derivation, theorem, constant or numerical experiment for it**. That section is therefore markedly shorter than the other two, which reflects the evidence rather than an omission.

The model problem is

$$
(-\Delta)^{\alpha/2}u(x)+\rho\,u(x)=f(x),\quad x\in\mathbb R^d,
\qquad u(x)\to0,\ |x|\to\infty,
$$

with $\rho>0$ and $\alpha\in(0,2)$. The fractional Laplacian has two equivalent definitions, the hypersingular integral

$$
(-\Delta)^{\alpha/2}u(x)=C_{d,\alpha}\,\mathrm{p.v.}\!\int_{\mathbb R^d}
\frac{u(x)-u(y)}{|x-y|^{\,d+\alpha}}\,\mathrm dy,
\qquad
C_{d,\alpha}=\alpha\,2^{\alpha-1}\,
\frac{\Gamma\bigl(\tfrac{d+\alpha}{2}\bigr)}{\pi^{d/2}\,\Gamma\bigl(\tfrac{2-\alpha}{2}\bigr)},
$$

and the pseudodifferential form $(-\Delta)^{\alpha/2}u=\mathcal F^{-1}[|\xi|^{\alpha}\mathcal F[u]]$.

## On an unbounded domain, what has to be matched is the decay rate

Choosing a basis on a bounded domain means choosing smoothness. Choosing one on an unbounded domain means choosing a **far-field decay rate**. How a basis function goes to zero as $|x|\to\infty$ decides how fast a finite expansion can consume the solution's tail; once the two decay types differ, accuracy fails in the tail first and then drags the global convergence rate down with it.

The two candidate families have tails of different kinds. Hermite functions carry a Gaussian factor $e^{-x^2/2}$ and therefore suit solutions with **exponential-type decay**; mapped rational (Gegenbauer) bases decay algebraically like $|x|^{-(\lambda+1)}$ with the exponent set by the parameter $\lambda$, and therefore suit **algebraically decaying** solutions. That single sentence is the whole difference between paper 27 and paper 42.

What forces the choice is the fractional operator itself. The hypersingular kernel $|x-y|^{-(d+\alpha)}$ decays slowly and the operator is nonlocal: values far away keep feeding back through the kernel, so **a solution of a fractional problem generally does not decay like a Gaussian even when its data does**. Paper 42 turns this into measurement: with an exponentially decaying source $e^{-x^2/2}(1+x)$ and an algebraically decaying one $(1+x^2)^{-2}$, both solutions decay along the same power law,

$$
\lvert u(x)\rvert\sim\lvert x\rvert^{-\alpha-1} .
$$

**The tail exponent is set by the operator, not by the data** — the ordinary Laplacian shows nothing of the kind, and this is the root of the triple failure described in the next section.

The achievements and the limits of both papers can be read off in advance from this. Paper 27's Hermite basis converges spectrally for exponentially decaying solutions, yet degrades to algebraic convergence on its own eigenvalue example, where the eigenfunctions decay only algebraically; paper 42 replaces the basis with rational functions whose tail exponent is tunable, precisely so that the basis tail can be aligned with $|x|^{-\alpha-1}$.

## All three standard unbounded-domain strategies fail here

Paper 42 states the difficulty precisely: solutions of fractional operators decay only **algebraically** at infinity, following a power law. That single fact defeats three standard strategies at once.

Domain truncation suits rapidly decaying solutions but not these; worse, naive truncation introduces nonphysical singularities at the interface where the unbounded domain is terminated. Transparent boundary conditions and sponge layers are highly nontrivial because the operator is nonlocal. Orthogonal functions on unbounded domains — Hermite and Laguerre — are tuned to **exponential** decay and are mismatched to power-law tails.

Mapped Jacobi ("rational") bases were already known to beat Hermite and Laguerre for **integer-order** problems with algebraically decaying solutions. Extending the mapping technique to the fractional case is far from trivial, and for a structural reason: the ordinary Laplacian maps a rational basis function to another function of the same type, while the fractional Laplacian maps it to **a class of functions of a completely different nature**.

## 27: Hermite bases with explicit differentiation matrices

### The idea: do the hypersingular integral once, keep an offline matrix

On bounded domains the successful approach is to pick a basis adapted to the fractional operator — Jacobi poly-fractonomials (eigenfunctions of a fractional Sturm-Liouville problem) and generalised Jacobi functions — for which a fractional derivative of a basis function is again a basis function with shifted parameters, making the operator effectively local. No such device exists on an unbounded domain.

Without that device the cost is immediate: the fractional Laplacian is nonlocal, so any discretisation is dense to begin with, and every entry is itself a hypersingular integral that would have to be evaluated repeatedly during the solve. The contemporaneous alternative (the Hermite spectral method of Mao and Shen) collocates in frequency space, so every evaluation costs a forward and a backward Hermite transform, which is awkward for nonlinear problems. Paper 27 collocates directly in **physical space** using explicit closed-form differentiation matrices.

What makes that possible is an unglamorous property: the Fourier image of a Hermite function is still a monomial in $\xi$ times a Gaussian. Multiplying by $|\xi|^{\alpha}$ and transforming back therefore lands inside a family of integrals with known closed forms, and the answer is a confluent hypergeometric ${}_1F_1$. **The hypersingular integral is thus done once, and done analytically**: the matrix is independent of the data, so it can be precomputed entirely offline, while the nonlinear term is evaluated pointwise at the collocation nodes and handled by Newton iteration. That is the substantive difference from the frequency-space route — there, every residual evaluation costs two Hermite transforms; here it costs one matrix-vector product.

### Setting

The paper's model problem carries a reaction term, linear or nonlinear:

$$
(-\Delta)^{\alpha/2}u(x)+\gamma f(u)=g(x),\quad x\in\mathbb R^d,
\qquad u(x)=0,\ |x|\to\infty .
$$

The fractional Laplacian is taken in the hypersingular form given at the top of this page and equivalently in the pseudodifferential form $\mathcal F[(-\Delta)^{\alpha/2}u](\xi)=|\xi|^{\alpha}\mathcal F[u](\xi)$, with $0<\alpha<2$ and the ordinary Laplacian recovered as $\alpha\to2$. **All derivations use the pseudodifferential definition**; the hypersingular form only serves to exhibit the nonlocality.

Hermite polynomials come from the three-term recurrence $H_0=1$, $H_1(x)=2x$, $H_{n+1}(x)=2xH_n(x)-2nH_{n-1}(x)$ with $\int_{\mathbb R}H_mH_n e^{-x^2}\mathrm dx=\gamma_n\delta_{mn}$, $\gamma_n=\sqrt{\pi}\,2^n n!$. Two bases are built from them. The normalised Hermite functions

$$
\widehat H_n(x)=\frac{1}{\sqrt{2^n n!}}\,e^{-x^2/2}H_n(x),
\qquad
\int_{\mathbb R}\widehat H_m\widehat H_n\,\mathrm dx=\sqrt{\pi}\,\delta_{mn}
$$

are orthogonal with respect to the weight $\omega\equiv1$; the over-scaled (Brinkman) basis

$$
\widetilde H_n(x)=e^{-x^2/2}\widehat H_n(x)=\frac{1}{\sqrt{2^n n!}}e^{-x^2}H_n(x),
\qquad
\int_{\mathbb R}\widetilde H_m\widetilde H_n\,e^{x^2}\mathrm dx=\sqrt{\pi}\,\delta_{mn}
$$

was introduced by Brinkman for Fokker-Planck equations, where the velocity part of the probability density is expanded in these functions, and has since become one of the standard bases there. The two differ in whether $e^{-x^2/2}$ sits in the basis or in the weight, and that choice determines the form of the discrete inner product and the differentiation matrices — and, as will be seen, the way the condition number grows, which is the paper's second theme.

### Derivation: the fractional image of the over-scaled basis

The collocation scheme expands the solution in the over-scaled basis, $u(x)\approx u_N(x)=\sum_{n=0}^{N-1}c_n\widetilde H_n(x)$, and imposes the equation at the roots $\{x_i\}_{i=0}^{N-1}$ of the $N$-th Hermite polynomial:

$$
\sum_{n=0}^{N-1}c_n\,(-\Delta)^{\alpha/2}\widetilde H_n(x_i)+\gamma f\bigl(u_N(x_i)\bigr)=g(x_i),
$$

that is, $\widetilde{\mathcal D}^{\alpha}c+\gamma F(c)=g$ with $\widetilde{\mathcal D}^{\alpha}_{i,j}=(-\Delta)^{\alpha/2}\widetilde H_j(x_i)$. **Everything therefore reduces to one task: evaluating $(-\Delta)^{\alpha/2}\widetilde H_n$.**

The chain has three steps. First, write down the Fourier image of the over-scaled basis, a monomial times a Gaussian:

$$
\mathcal F\bigl[\widetilde H_{2n}\bigr](\xi)
=\frac{(-1)^n}{\sqrt2\,\sqrt{2^{2n}(2n)!}}\,\xi^{2n}e^{-\xi^2/4}.
$$

Second, multiply by $|\xi|^{\alpha}$ per the pseudodifferential definition and invert; by evenness the inverse transform collapses to a cosine integral,

$$
\int_{\mathbb R^+}\xi^{2n+\alpha}e^{-\xi^2/4}\cos(x\xi)\,\mathrm d\xi .
$$

Third, finish with the known evaluation of that integral, which produces a ${}_1F_1$. The result is Theorems 3.1 and 3.2: for $0<\alpha<2$,

$$
(-\Delta)^{\alpha/2}\widetilde H_{2n}(x)=\frac{2^{\alpha}(-1)^n}{\sqrt{(2n)!}}\cdot\frac{1}{2^n n!}\cdot
\frac{\Gamma\bigl(n+\tfrac{\alpha}{2}+\tfrac12\bigr)}{\Gamma\bigl(n+\tfrac12\bigr)}\;
{}_1F_1\!\Bigl(n+\tfrac{\alpha}{2}+\tfrac12,\ \tfrac12;\ -x^2\Bigr),
$$

$$
(-\Delta)^{\alpha/2}\widetilde H_{2n+1}(x)=\frac{2^{\alpha+1}(-1)^n}{\sqrt{(2n+1)!}}\cdot\frac{1}{2^{\,n+\frac12}n!}\cdot
\frac{\Gamma\bigl(n+\tfrac{\alpha}{2}+\tfrac32\bigr)}{\Gamma\bigl(n+\tfrac32\bigr)}\;
x\,{}_1F_1\!\Bigl(n+\tfrac{\alpha}{2}+\tfrac32,\ \tfrac32;\ -x^2\Bigr).
$$

Note that the even and odd cases differ only in the second ${}_1F_1$ parameter ($\tfrac12$ against $\tfrac32$) and the extra factor $x$ in front, which is the trace of the fractional operator preserving parity.

> [!warning] One typographic uncertainty
> In the two formulas above, the grouping of $\sqrt{(2n)!}$ relative to $2^n n!$ cannot be pinned down uniquely from the text available here, which came from a PDF extraction; it should be checked against the typeset paper when the constants matter. The $\Gamma$-ratios, the three ${}_1F_1$ parameters and the prefactors $2^{\alpha}$ and $2^{\alpha+1}$ are unambiguous in the source, and the proof is fully legible.

One more saving matters in practice. Calling ${}_1F_1$ for every entry is too expensive, so the paper uses the closed forms only for $j=0,1,2,3$ — those do require a fast and accurate ${}_1F_1$ routine — and fills the remaining $4<j\le N-1$ from the contiguous relation

$$
(2a-b+x)\,{}_1F_1(a,b;x)=a\,{}_1F_1(a+1,b;x)-(b-a)\,{}_1F_1(a-1,b;x).
$$

The matrix does not depend on the data, so the whole table is computed once, offline.

### Derivation: the normalised basis goes through Fourier eigenfunctions

The normalised basis needs a different entry point, because it consists of eigenfunctions of the Fourier transform:

$$
\mathcal F\bigl[\widehat H_n\bigr](\xi)
=\frac{1}{\sqrt{2\pi}}\int_{\mathbb R}\widehat H_n(x)e^{-i\xi x}\mathrm dx
=(-i)^n\widehat H_n(\xi).
$$

Expanding $\widehat H_n$ itself as $\widehat H_n(\xi)=\sum_{k=0}^{n}\hat a_{n,k}e^{-\xi^2/2}\xi^k$ with $\hat a_{n,k}=a_{n,k}/\sqrt{2^n n!}$, where the polynomial coefficients are generated by

$$
a_{0,0}=1,\quad a_{1,0}=0,\quad a_{1,1}=2,
\qquad
a_{n+1,0}=-a_{n,1},
\qquad
a_{n+1,k}=2a_{n,k-1}-(k+1)a_{n,k+1}\ (k>0),
$$

confines the fractional operator to the single family $e^{-\xi^2/2}\xi^k$. Setting $F_k(x)=\mathcal F^{-1}\bigl[e^{-\xi^2/2}\xi^k|\xi|^{\alpha}\bigr](x)$, both parities are available in closed form:

$$
F_{2m}(x)=\frac{2^{\,m+\frac{\alpha}{2}}}{\sqrt{\pi}}\,
\Gamma\!\Bigl(\frac{2m+1+\alpha}{2}\Bigr)\,
{}_1F_1\!\Bigl(\frac{2m+1+\alpha}{2},\ \frac12;\ -\frac{x^2}{2}\Bigr),
$$

$$
F_{2m+1}(x)=\frac{2^{\,m+1+\frac{\alpha}{2}}\,i}{\sqrt{\pi}}\,
\Gamma\!\Bigl(\frac{2m+3+\alpha}{2}\Bigr)\,
x\,{}_1F_1\!\Bigl(\frac{2m+3+\alpha}{2},\ \frac32;\ -\frac{x^2}{2}\Bigr),
$$

and the differentiation matrix becomes a single sum,

$$
\widehat{\mathcal D}^{\alpha}_{mn}=(-\Delta)^{\alpha/2}\widehat H_n(x_m)
=(-i)^n\sum_{k=0}^{N-1}\hat a_{n,k}F_k(x_m),
\qquad 0\le n,m\le N-1,
$$

with the convention $a_{n,k}=0$ for $k>n$. Note that the argument of $F_k$ is $-x^2/2$ where the over-scaled formulas have $-x^2$; the two sets of formulas are not interchangeable.

The same approximation space also admits a Lagrange-type nodal basis $h_j(x_k)=\delta_{jk}$, which yields **a different** differentiation matrix assembled from the same $(-\Delta)^{\alpha/2}\widehat H_k(x_i)$ values — the two matrices have different condition numbers, a fact the next subsection uses.

### Derivation: scaling factor, multi-term models and two dimensions

**Scaling factor.** Hermite-type spectral convergence degrades when the decay rate of the solution and that of the basis are mismatched. If $|u(x)|\approx0$ for all $|x|>M$, expand instead

$$
u(x)=\sum_{n=0}^{N-1}c_n\widetilde H_n(rx)
\iff
u(x/r)=\sum_{n=0}^{N-1}c_n\widetilde H_n(x),\qquad r>0,
$$

so that the scaled collocation nodes $\{x_k/r\}$ sit inside the effective support of $u$. The suggested choice is

$$
\max_{0\le k\le N-1}|x_k|/r\le M
\ \Longrightarrow\
r=\max_{0\le k\le N-1}|x_k|/M ,
$$

and the paper notes explicitly that $M$ itself is nontrivial to find in practice, so the optimal $r$ is generally hard to obtain. With scaling, the matrix entries are the same formulas evaluated at $z_i=rx_i$; no new derivation is needed.

**Multi-term and distributed-order models.** For $\sum_{j=1}^{J}(-\Delta)^{\alpha_j/2}u+\gamma f(u)=g$ — motivated by quadrature approximation of distributed-order models — the collocation system is simply

$$
\mathcal Jc+\gamma F(c)=g,
\qquad
\mathcal J=\sum_{j=1}^{J}\widetilde{\mathcal D}^{\alpha_j},
$$

so the differentiation matrices just add. **That is a real dividend of the closed-form route: one more fractional term costs nothing structurally.**

**Two dimensions.** With the tensor-product basis $\{\widetilde H_n(x)\widetilde H_m(y)\}$ and $\mathcal F$ factorising,

$$
(-\Delta)^{\alpha/2}\bigl[\widetilde H_n(x)\widetilde H_m(y)\bigr]
=\frac{1}{2\pi}\int_{\mathbb R}\int_{\mathbb R}(\xi^2+\eta^2)^{\alpha/2}
\mathcal F[\widetilde H_n](\xi)\,\mathcal F[\widetilde H_m](\eta)\,
e^{ix\xi}e^{iy\eta}\,\mathrm d\xi\,\mathrm d\eta .
$$

The angular integral is removed with Bessel identities — the key ones being $J_{-1/2}(x)=\sqrt{2/(\pi x)}\cos x$ and a Sonine-Gegenbauer-type formula — leaving each entry a single ${}_1F_1$ value. The double-even block reads

$$
\widetilde{\mathcal D}^{\alpha}_{(i-1)N+j,\,(p-1)N+q}
=\frac{2^{\alpha}\,\Gamma\bigl(n+m+\tfrac{\alpha}{2}+1\bigr)}
{\sqrt{2^{2n}(2n)!}\,\sqrt{2^{2m}(2m)!}\;\Gamma(2n+2m+1)}\;
{}_1F_1\!\Bigl(n+m+\tfrac{\alpha}{2}+1;\ 2n+2m+1;\ -(x_i^2+y_j^2)\Bigr)
$$

for $p=2n$, $q=2m$, with three analogous cases. There is a pleasant structure here: **the two-dimensional entries depend on the nodes only through $x_i^2+y_j^2$**, so the second dimension does not inflate the complexity of the formulas.

### Theorems: closed forms, but no convergence theorem

This is an algorithmic paper and **it contains no convergence theorem**. Its rigorous content is exactly the closed-form evaluations above (Theorems 3.1-3.2 for the over-scaled basis and the $F_k$ formulas for the normalised one). Spectral convergence is an empirical statement: the introduction claims that "both methods admit spectral convergence for solutions with exponential decay at infinity", and that is an assertion rather than a proved theorem, which is how this page records it.

Where the paper does deliver quantitative conclusions is conditioning. For the over-scaled basis $\{\widetilde H_n\}$ the condition number of $\widetilde{\mathcal D}^{\alpha}$ "grows very fast with respect to $N$" (Figure 1 spans $10^0$ to $10^{10}$ for $\alpha=0.4,1,1.6$), which the paper attributes to the poor conditioning of that basis and which motivates the normalised alternative. For the normalised basis the growth is **algebraic**, with fitted slopes reported in the figures:

| Differentiation matrix                     | $\alpha=0.4$ | $\alpha=1$ | $\alpha=1.6$ |
| ------------------------------------------ | ------------ | ---------- | ------------ |
| $\widehat{\mathcal D}^{\alpha}$ (Figure 2) | $0.446$      | $1.065$    | $1.684$      |
| Lagrange nodal-basis matrix (Figure 3)     | $0.428$      | $1.054$    | $1.659$      |

The three slopes nearly equal the three values of $\alpha$, which reads as growth like $N^{\alpha}$; **but identifying the exponent with $\alpha$ is this site's reading of those slopes, not a statement of the paper.** The paper does state that efficient preconditioners should be designed for practical use.

### Numerical experiments: six examples

Section 6 gives six groups of experiments. Error magnitudes and convergence slopes appear there only in figures, and could not be recovered from the text available here, so **the table below reports the setups and the verifiable qualitative conclusions and gives no error values**.

| Example                          | Equation and exact solution                                                                                                   | Scaling factor                                         | What is verifiable                                                                                                                                  |
| -------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| §6.1 1-D benchmark               | exact solution $u(x)=e^{-x^2}\sin x$                                                                                          | over-scaled $r=1$ (none needed); normalised $r=\sqrt2$ | spectral convergence for $\alpha=0.4,1,1.6$ in a weighted and in the maximum norm; the over-scaled run is polluted at large $N$ by condition growth |
| §6.2 scaling (over-scaled basis) | exact solution $u(x)=e^{-x^2/2}x^2\cos x$                                                                                     | $r=1/\sqrt2$ against $r=1$                             | the right scaling factor gives a visibly faster rate                                                                                                |
| §6.2 scaling (normalised basis)  | exact solution $u(x)=e^{-2x^2}x^2\cos x$                                                                                      | optimal $r=2$ against $r=1$                            | same conclusion                                                                                                                                     |
| §6.3 two dimensions              | $(-\Delta)^{\alpha/2}u+2u=g$, exact solution $u=e^{-(x^2+y^2)}\sin(x+y)$                                                      | tensor-product basis                                   | only the setup could be verified                                                                                                                    |
| §6.4 multi-term model            | $\sum_{j=1}^{4}(-\Delta)^{\alpha_j/2}u=g$, exact solution $u=e^{-3x^2/2}(\sin x+x^6+x^2\cos x)$                               | $r=\sqrt{1.5},\sqrt{1.3},1$ compared                   | $r=\sqrt{1.5}$ clearly best                                                                                                                         |
| §6 nonlinear                     | $(-\Delta)^{\alpha/2}u+u^2=g$, exact solutions $e^{-x^2}(\sin x+x^2)$ (over-scaled) and $e^{-x^2/2}(\sin x+x^2)$ (normalised) | not reported                                           | Newton iteration with tolerance $10^{-16}$; spectral convergence observed                                                                           |
| §6 eigenvalue problem            | $\bigl[(-\Delta)^{\alpha/2}+x^2\bigr]u=\lambda u$                                                                             | not reported                                           | only **algebraic** decay of the error                                                                                                               |

The four orders in the multi-term model are transformed Legendre-Gauss points:

| $j$        | 1       | 2       | 3       | 4       |
| ---------- | ------- | ------- | ------- | ------- |
| $\alpha_j$ | $0.139$ | $0.660$ | $1.340$ | $1.861$ |

The eigenvalue example is the most diagnostic of the six because it has an independent reference: for $\alpha=1$ the spectrum is known, $\lambda_{2k-1}=-a'_k$ and $\lambda_{2k}=-a_k$, where $a_k$ and $a'_k$ are the decreasingly ordered roots of the Airy function

$$
A(x)=\frac1\pi\int_0^{\infty}\cos\Bigl(\frac{t^3}{3}+xt\Bigr)\mathrm dt
$$

and of its derivative. The reference values used are

| Eigenvalue  | Reference value ($\alpha=1$) |
| ----------- | ---------------------------- |
| $\lambda_1$ | $1.01879297164747$           |
| $\lambda_2$ | $2.33810741045976$           |
| $\lambda_3$ | $3.24819758217983$           |

and here only algebraic convergence is observed, which the paper attributes to the merely **algebraic** decay of the eigenfunctions.

What this set of experiments establishes: physical-space collocation with a precomputed matrix does deliver spectral convergence for exponentially decaying solutions, uniformly across the three values of $\alpha$; the scaling factor is not cosmetic but sets the rate; the multi-term extension really does amount to adding matrices; and nonlinear problems can be attacked with Newton directly, which is the case the frequency-space method handles awkwardly.

It falls short of the theory in three places. First, nothing underwrites the observed spectral convergence, and the conditioning results are fitted slopes; since the over-scaled basis has its convergence "polluted" at large $N$, **the accuracy actually reachable is set by conditioning rather than by approximation properties**, and the paper can only list preconditioning as future work. Second, the optimal scaling factor depends on the effective support $M$, which the paper admits is hard to find, so the well-chosen $r$ values of §6.2 may be unavailable on a real problem. Third, the eigenvalue example draws the method's boundary: as soon as the solution decays only algebraically, spectral convergence collapses to algebraic. **That is precisely the problem paper 42 sets out to fix, and paper 27 has already demonstrated it.**

## 42: decomposing the rational basis into two elementary shapes

### The idea: align the basis tail with the tail the operator produces

Every design decision in paper 42 follows from the empirical fact stated at the top of this page: the solution's tail is $|x|^{-\alpha-1}$. For a spectral approximation to consume such a tail, the basis functions must themselves decay by a power law, and the power must be tunable — anything frozen to Gaussian or exponential type will not do, because $\alpha$ is a parameter of the problem and the tail exponent moves with it. Modified mapped Gegenbauer functions are exactly such a family: they decay like $|x|^{-(\lambda+1)}$ with $\lambda$ as the knob on the tail exponent, and a separate scaling parameter $\mu$ controls the width of the profile. The two decay rates can then genuinely be aligned instead of being brute-forced by raising $N$.

The price shows up at the other end, and it is structural. The ordinary Laplacian sends a rational basis function to another rational function of the same type, so the whole mapped-spectral machinery is self-contained at integer order; the fractional Laplacian is not self-contained, it sends the rational basis into hypergeometric functions and, as shown below, even the decay rate need not follow the basis parameter. So the paper's task is concrete: **decompose the basis into two elementary rational shapes, work out $(-\Delta)^{s}$ for those two once and for all, and recombine linearly.** Whether those two closed forms exist is what the whole paper turns on — with them one gets an offline differentiation matrix isomorphic to paper 27's; without them the mapping technique simply does not reach the fractional setting.

### Setting: the map and the basis

The one-to-one map and its identities are

$$
x=\frac{t}{\sqrt{1-t^2}}
\ \Longleftrightarrow\
t=\frac{x}{\sqrt{1+x^2}},
\qquad
1-t^2=\frac{1}{1+x^2},
\qquad
\frac{\mathrm dx}{\mathrm dt}=\frac{1}{(1-t^2)^{3/2}} .
$$

This is the **singular** (unbounded-interval) algebraic map, and it is exactly what produces power-law tails in the basis. The modified mapped Gegenbauer functions are

$$
R^{\lambda}_n(x):=\bigl(1+x^2\bigr)^{-\frac{\lambda+1}{2}}
C^{\lambda}_n\!\Bigl(\frac{x}{\sqrt{1+x^2}}\Bigr),
\qquad x\in\mathbb R,\ \lambda>-\tfrac12 .
$$

"Modified" means precisely that the basis **absorbs the square root of the weight function**: writing $R^{\lambda}_n(x)=S(t)C^{\lambda}_n(t)$ with

$$
S(t)=\sqrt{\omega_{\lambda}(t)\frac{\mathrm dt}{\mathrm dx}}=(1-t^2)^{\frac{\lambda+1}{2}} ,
$$

these functions are orthogonal with respect to the **uniform** weight on $\mathbb R$:

$$
\int_{-\infty}^{\infty}R^{\lambda}_n(x)R^{\lambda}_m(x)\,\mathrm dx=\gamma^{\lambda}_n\delta_{nm},
\qquad
\gamma^{\lambda}_n=\frac{\pi\,2^{1-2\lambda}\Gamma(n+2\lambda)}{n!\,(n+\lambda)\,\Gamma^2(\lambda)} .
$$

The basis satisfies a three-term recurrence

$$
n\,R^{\lambda}_n(x)=\frac{2x}{\sqrt{1+x^2}}(n+\lambda-1)R^{\lambda}_{n-1}(x)
-(n+2\lambda-2)R^{\lambda}_{n-2}(x),\qquad n\ge2,
$$

$$
R^{\lambda}_0(x)=\frac{1}{(1+x^2)^{\frac{\lambda+1}{2}}},
\qquad
R^{\lambda}_1(x)=\frac{2\lambda x}{(1+x^2)^{1+\frac{\lambda}{2}}} .
$$

The far-field behaviour is the point:

$$
\lim_{x\to+\infty}(1+x^2)^{\frac{\lambda+1}{2}}R^{\lambda}_n(x)=\frac{(2\lambda)_n}{n!},
$$

so every basis function decays like $|x|^{-(\lambda+1)}$ and $\lambda$ is the knob that tunes the basis to the solution's power law.

### Derivation: decomposition into two elementary shapes

Rewriting the Gegenbauer polynomial by the Pfaff transformation as a ${}_2F_1$ in the variable $(1+x^2)^{-1}$ expands the basis functions into finite sums:

$$
R^{\lambda}_{2n}(x)=a^{\lambda}_n\sum_{k=0}^{n}
\frac{(-n)_k(n+\lambda)_k}{(\lambda+\tfrac12)_k\,k!}\cdot
\frac{1}{(1+x^2)^{\,k+\frac{\lambda+1}{2}}},
$$

$$
R^{\lambda}_{2n+1}(x)=b^{\lambda}_n\sum_{k=0}^{n}
\frac{(-n)_k(n+\lambda+1)_k}{(\lambda+\tfrac12)_k\,k!}\cdot
\frac{x}{(1+x^2)^{\,k+\frac{\lambda}{2}+1}} .
$$

The two constants $a^{\lambda}_n$ and $b^{\lambda}_n$ are built from $(-1)^n$, $(\lambda+n)$ and the Beta function $B(\lambda,n+1)$, but **their exact expressions are corrupted in the text available here and could not be verified**. That does not affect the structural conclusion: whatever the constants are, the expansion contains only two shapes, so

**every basis function is a finite linear combination of just two elementary shapes**:

$$
\frac{1}{(1+x^2)^{\gamma}}\ \bigl(\gamma=k+\tfrac{\lambda+1}{2}\bigr),
\qquad
\frac{x}{(1+x^2)^{\gamma}}\ \bigl(\gamma=k+\tfrac{\lambda}{2}+1\bigr).
$$

It therefore suffices to differentiate those two fractionally. For real $s>0$:

$$
(-\Delta)^{s}\Bigl\{\frac{1}{(1+x^2)^{\gamma}}\Bigr\}
=A^{\gamma}_s\ {}_2F_1\!\Bigl(s+\gamma,\ s+\tfrac12;\ \tfrac12;\ -x^2\Bigr),
\qquad \gamma>0,
$$

$$
(-\Delta)^{s}\Bigl\{\frac{x}{(1+x^2)^{\gamma}}\Bigr\}
=(2s+1)A^{\gamma}_s\,x\ {}_2F_1\!\Bigl(s+\gamma,\ s+\tfrac32;\ \tfrac32;\ -x^2\Bigr),
\qquad \gamma>\tfrac12,
$$

$$
A^{\gamma}_s=\frac{2^{2s}\,\Gamma(s+\gamma)\,\Gamma\bigl(s+\tfrac12\bigr)}{\sqrt{\pi}\,\Gamma(\gamma)} .
$$

The proof route is worth recording because it shows where modified Bessel functions enter. Write $v(x)=(1+x^2)^{-\gamma}$ as ${}_2F_1(\gamma,\tfrac12;\tfrac12;-x^2)$; a cosine-transform identity gives a modified Bessel function of the second kind,

$$
\mathcal F[v](\xi)=\frac{2^{1-\gamma}}{\Gamma(\gamma)}\,\xi^{\gamma-\frac12}
K_{\gamma-\frac12}(\xi),\qquad \xi>0,
$$

and then $(-\Delta)^s v=\mathcal F^{-1}[|\xi|^{2s}\mathcal F[v]]$ becomes a Bessel integral,

$$
(-\Delta)^{s}v=\frac{2^{1-\gamma}}{\sqrt{2\pi}\,\Gamma(\gamma)}
\int_0^{\infty}\cos(x\xi)\,\xi^{2s+\gamma-\frac12}K_{\gamma-\frac12}(\xi)\,\mathrm d\xi,
$$

which the corresponding integral identity finishes with $\lambda=2s+\gamma-\tfrac12$, $\mu=\gamma-\tfrac12$ and $b=x$. The odd shape follows the same route through $\mathcal F[xv](\xi)=-i\frac{2^{1-\gamma}}{\Gamma(\gamma)}\xi^{\gamma-\frac12}K_{\gamma-\frac32}(\xi)$ together with $\Gamma(z+1)=z\Gamma(z)$. **Modified Bessel functions are unavoidable on this route**, which is why the abstract names integral identities related to them as the source of the representations.

The same result has an equivalent series form obtained by the Pfaff/Euler transformation (Corollary 3.3); the paper notes that the first form is more convenient for computation and the second more suitable for analysis:

$$
(-\Delta)^{s}\Bigl\{\frac{1}{(1+x^2)^{\gamma}}\Bigr\}
=\frac{A^{\gamma}_s}{(1+x^2)^{s+\gamma}}\
{}_2F_1\!\Bigl(-s,\ s+\gamma;\ \tfrac12;\ \tfrac{x^2}{1+x^2}\Bigr),
$$

$$
(-\Delta)^{s}\Bigl\{\frac{x}{(1+x^2)^{\gamma}}\Bigr\}
=\frac{(2s+1)A^{\gamma}_s\,x}{(1+x^2)^{s+\gamma}}\
{}_2F_1\!\Bigl(-s,\ s+\gamma;\ \tfrac32;\ \tfrac{x^2}{1+x^2}\Bigr).
$$

The second form pulls the decay factor $(1+x^2)^{-(s+\gamma)}$ out in front and confines the ${}_2F_1$ argument to $[0,1)$ — **which is the entry point for the decay analysis in the next subsection**, since the whole difference between integer and fractional order lies in whether that ${}_2F_1$ terminates.

### Derivation: where integer and fractional order diverge

This is the most informative structural observation in the paper. For positive integer $s$ both hypergeometric functions terminate and the operator **gains $s$ powers of decay**:

$$
(-\Delta)^{s}\Bigl\{\tfrac{1}{(1+x^2)^{\gamma}}\Bigr\}\sim\frac{1}{(1+x^2)^{s+\gamma}} .
$$

For noninteger $s$ the hypergeometric functions can diverge as $|x|\to\infty$ and the picture changes. For $\gamma>1/2$,

$$
(-\Delta)^{s}\bigl\{(1+x^2)^{-\gamma}\bigr\}\sim(1+x^2)^{-(s+1/2)},
$$

so **the decay rate is independent of $\gamma$**; at $\gamma=1/2$ a logarithm appears, $\ln(1+x^2)/(1+x^2)^{s+1/2}$; and only $0<\gamma<1/2$ retains the integer-order behaviour. At the level of the basis,

$$
(-\Delta)^{s}R^{\lambda}_{2n}(x)\sim
\begin{cases}
(1+x^2)^{-\left(s+\frac{\lambda+1}{2}\right)}, & -\tfrac12<\lambda<0,\\[3pt]
\ln(1+x^2)\,(1+x^2)^{-\left(s+\frac12\right)}, & \lambda=0,\\[3pt]
(1+x^2)^{-\left(s+\frac12\right)}, & \lambda>0 .
\end{cases}
$$

In other words, the fractional Laplacian applied to this basis **does not always gain the factor $1/(1+x^2)^s$**. That is the quantitative version of the "completely different nature" statement in the introduction.

### The fractional image of the basis, and two implementation details

Substituting the two closed forms back into the Pfaff expansion gives what the abstract calls the "main building block", namely the fractional image of the basis functions themselves (Theorem 3.4, for real $s>0$ and $\lambda>-1/2$):

$$
(-\Delta)^{s}R^{\lambda}_{2n}(x)=a^{\lambda}_n\sum_{k=0}^{n}
\frac{(-n)_k(n+\lambda)_k}{(\lambda+\tfrac12)_k\,k!}\;
A^{\,k+\frac{\lambda+1}{2}}_{s}\
{}_2F_1\!\Bigl(s+k+\tfrac{\lambda+1}{2},\ s+\tfrac12;\ \tfrac12;\ -x^2\Bigr),
$$

$$
(-\Delta)^{s}R^{\lambda}_{2n+1}(x)=(2s+1)\,b^{\lambda}_n\,x\sum_{k=0}^{n}
\frac{(-n)_k(n+\lambda+1)_k}{(\lambda+\tfrac12)_k\,k!}\;
A^{\,k+\frac{\lambda}{2}+1}_{s}\
{}_2F_1\!\Bigl(s+k+\tfrac{\lambda}{2}+1,\ s+\tfrac32;\ \tfrac32;\ -x^2\Bigr).
$$

The formulas for the fractional Laplacian of the basis are less compact than in the Hermite case, as the paper itself says, but they evaluate efficiently by recurrence, and the paper reports that the fractional Laplacian of the rational basis can be computed accurately up to degrees of order $10^3$ with Maple or Mathematica. With $F_k(x)={}_2F_1(a,b;c;-x^2)$, $a=s+k+\tfrac{\lambda+1}{2}$, $b=s+\tfrac12$ and $c=\tfrac12$,

$$
F_{k+1}(x)=\frac{c-a}{a(1+x^2)}F_{k-1}(x)
+\frac{(2a-c)+(a-b)x^2}{a(1+x^2)}F_k(x),\qquad k\ge1 .
$$

A scaling parameter $\mu>0$ improves resolution: the map becomes $x=\mu t/\sqrt{1-t^2}$ and the basis becomes

$$
R^{\lambda}_{n,\mu}(x)=\frac{\mu^{\lambda+\frac12}}{(\mu^2+x^2)^{\frac{\lambda+1}{2}}}
C^{\lambda}_n\!\Bigl(\frac{x}{\sqrt{\mu^2+x^2}}\Bigr)
=\mu^{-\frac12}R^{\lambda}_n\!\Bigl(\frac{x}{\mu}\Bigr).
$$

This is the rational analogue of the Hermite scaling factor in paper 27, with an important difference: here $\mu$ controls the **width** of the algebraic profile while $\lambda$ independently controls the **exponent** of the tail. Separating the two is the practical advantage of the rational basis over the Hermite one.

### Approximation space and two schemes

With $\mathcal V^{\lambda}_N=\mathrm{span}\{R^{\lambda}_n\}=\{S(t)P(t):P\in\mathcal P_N\}$, the change of variables $\breve u(x)=u(x)/s(x)=\breve U(t)$ with $s(x)=(1+x^2)^{-(\lambda+1)/2}=S(t)$ gives $\hat u^{\lambda}_n=\widehat{\breve U}_n$ exactly, so modified mapped Gegenbauer approximation on $\mathbb R$ **is** Gegenbauer approximation on $(-1,1)$. The mapped Gauss nodes and weights are

$$
x^{\lambda}_j=\frac{t^{\lambda}_j}{\sqrt{1-(t^{\lambda}_j)^2}},
\qquad
\omega^{\lambda}_j=\bigl(1+(t^{\lambda}_j)^2\bigr)^{-\lambda}\rho^{\lambda}_j .
$$

The Galerkin scheme uses $\tilde a_s(u,v)=((-\Delta)^{s/2}u,(-\Delta)^{s/2}v)+\rho(u,v)$ over $\mathcal V^{\lambda}_N$ with the source term interpolated, which is why the error bound carries a separate $f$ term. The collocation scheme uses the Lagrange basis on the mapped Gauss grid to build fractional differentiation matrices

$$
\mathcal D_{i,j}=(-\Delta)^{\alpha/2}l_j(x^{\lambda}_i)
=\sum_{k=0}^{N-1}b^j_k\,(-\Delta)^{\alpha/2}R^{\lambda}_k(x^{\lambda}_i),
\qquad
b^j_k=\frac{R^{\lambda}_k(x^{\lambda}_j)\,\omega^{\lambda}_j}{\gamma^{\lambda}_k} .
$$

The paper states explicitly that convergence analysis of the collocation scheme "seems nontrivial and largely open", and that self-assessment is worth preserving.

The multidimensional case moves to the Fourier domain, where the problem is diagonal:

$$
(|\xi|^{\alpha}+\rho)\hat u(\xi)=\hat f(\xi),
\qquad
\hat u(\xi)=a(\xi)\hat f(\xi),\quad a(\xi)=\frac{1}{|\xi|^{\alpha}+\rho} .
$$

The algorithm interpolates $f$ on the tensor grid, transforms, multiplies by $a(\xi)$ and projects back. A self-dual structure helps: the inverse transform of $R^{\lambda}_n$ is computed by the **same** expansion formulas. There is an extra reason to choose a rational basis here — a function's Fourier transform typically decays more slowly than the function itself, so an algebraically decaying basis is exactly the right tool in the transformed variable.

### Theorems: optimal estimates for projection, interpolation and Galerkin

The error theory uses two spaces. One is the usual fractional Sobolev space

$$
H^r(\mathbb R)=\Bigl\{u\in L^2:\int_{\mathbb R}(1+|\xi|^2)^r\bigl|\mathcal F[u](\xi)\bigr|^2\mathrm d\xi<\infty\Bigr\};
$$

the other is tailored to the rational basis and measures regularity through **mapped derivatives** rather than ordinary ones:

$$
D_xu:=a(x)\frac{\mathrm d\breve u}{\mathrm dx}=\frac{\mathrm d\breve U}{\mathrm dt},
\qquad
D^2_xu=a(x)\frac{\mathrm d}{\mathrm dx}\Bigl\{a(x)\frac{\mathrm d\breve u}{\mathrm dx}\Bigr\}=\frac{\mathrm d^2\breve U}{\mathrm dt^2},\ \dots
$$

where $a(x)$ is the factor supplied by the map and the second equality identifies $D_x$ as differentiation in the reference variable $t$. The seminorm is

$$
\lvert u\rvert_{\mathbb B^m_{\lambda}(\mathbb R)}
=\bigl\|(1+x^2)^{-\frac{\lambda+m+1}{2}}D^m_xu\bigr\|_{L^2(\mathbb R)} .
$$

**This choice is why the estimates come out at optimal order**: regularity in $\mathbb B^m_{\lambda}$ is not ordinary smoothness of the solution on $\mathbb R$ but smoothness after pulling it back to $(-1,1)$, and a power-law tail is smooth once pulled back.

The estimates below all follow one route. The basic tool is a space-interpolation lemma (Lemma 4.1): for $r=(1-\theta)r_0+\theta r_1$ with $\theta\in[0,1]$, $\|u\|_{H^r}\le\|u\|^{1-\theta}_{H^{r_0}}\|u\|^{\theta}_{H^{r_1}}$, proved by Hölder's inequality with exponents $p=1/(1-\theta)$ and $q=1/\theta$ on the Fourier side.

**Theorem 4.2 (optimal $L^2$-projection error).** Let $u\in H^s(\mathbb R)\cap\mathbb B^m_{\lambda}(\mathbb R)$ with integer $1\le m\le N+1$, $s\in(0,1)$ and $\lambda>-1/2$. Then

$$
\|\pi^{\lambda}_Nu-u\|_{H^s(\mathbb R)}\le c\,N^{\,s-m}\,|u|_{\mathbb B^m_{\lambda}(\mathbb R)},
$$

with $c$ independent of $N$ and $u$. The proof route is the equivalence established above: transfer the problem to Gegenbauer approximation on $(-1,1)$, apply the known estimate $\|(\Pi^{\lambda}_N\Phi-\Phi)^{(l)}\|_{L^2_{\omega_{\lambda+l}}}\le cN^{\,l-m}\|\Phi^{(m)}\|_{L^2_{\omega_{\lambda+m}}}$ for $l=0,1$, then interpolate between $L^2$ and $H^1$ to reach $H^s$.

**Theorem 4.3 ($H^s$-orthogonal projection).** With $a_s(u,v)=\bigl((-\Delta)^{s/2}u,(-\Delta)^{s/2}v\bigr)+(u,v)$ and the induced projection $\pi^{s}_{N,\lambda}$, that projection is the best approximation in $H^s$ and has the same order:

$$
\|\pi^{s}_{N,\lambda}u-u\|_{H^s}=\inf_{\phi\in\mathcal V^{\lambda}_N}\|\phi-u\|_{H^s}
\le cN^{\,s-m}|u|_{\mathbb B^m_{\lambda}} .
$$

**Theorem 4.4 (interpolation).** The mapped Gegenbauer-Gauss interpolant has the same order: $\|I^{\lambda}_Nu-u\|_{H^s(\mathbb R)}\le cN^{\,s-m}|u|_{\mathbb B^m_{\lambda}(\mathbb R)}$.

**Theorem 5.1 (optimal convergence of the Galerkin scheme).** Let $u\in H^s(\mathbb R)\cap\mathbb B^m_{\lambda}(\mathbb R)$ and $f\in\mathbb B^{k}_{\lambda}(\mathbb R)$ with integers $1\le m,k\le N+1$, $s=\alpha/2\in(0,1)$ and $\lambda>-1/2$. Then

$$
\|u-u_N\|_{H^s(\mathbb R)}\le c\,N^{\,s-m}\,|u|_{\mathbb B^m_{\lambda}(\mathbb R)}
+c\,N^{-k}\,|f|_{\mathbb B^{k}_{\lambda}(\mathbb R)},
$$

with $c$ independent of $N$, $u$ and $f$. The two terms have different origins: the first is the approximation error of the solution, the second comes from interpolating the source term. **Section 5.2 points out that the second term usually dominates in practice**, and the numerical results in the next subsection are the direct consequence of that.

In several dimensions the bound has the same shape (eq. (7.21)): if $a\hat f\in\mathbb B^m_{\lambda}(\mathbb R^d)$ and $f\in\mathbb B^{m'}_{\lambda}(\mathbb R^d)$ with $m\ge0$ and $m'\ge d$, then

$$
\|u_{NM}-u\|_{L^2(\mathbb R^d)}\le c\,M^{-m}\,|a\hat f|_{\mathbb B^m_{\lambda}(\mathbb R^d)}
+c\,N^{-m'}\,|f|_{\mathbb B^{m'}_{\lambda}(\mathbb R^d)} ,
$$

where $M$ counts degrees of freedom on the frequency side and $N$ those used to interpolate $f$; the two terms rest on $\|\pi^{\lambda}_Mu-u\|_{L^2}\le cM^{-m}|u|_{\mathbb B^m_{\lambda}}$ and $\|I^{\lambda}_Nu-u\|_{L^2}\le cN^{-m}|u|_{\mathbb B^m_{\lambda}}$ for $m\ge d$. Remark 7.2 records something the paper does not pursue: replacing $\pi^{\lambda}_M$ by $I^{\lambda}_M$ would require stability of the interpolant in norms involving partial derivatives of order $d$.

> [!warning] Two qualifications
> First, the exponent of the weight in the $\mathbb B^m_{\lambda}$ seminorm is ambiguous in the text available here, so the $-\frac{\lambda+m+1}{2}$ written above is a **reconstruction**; only the structure — some power of $(1+x^2)$ times $D^m_xu$, measured in $L^2$ — is certain. Second, all of these theorems cover the Galerkin scheme only. **The collocation scheme has no convergence theory**; the paper says plainly that it "seems nontrivial and largely open", and the reason is visible: the fractional Laplacian pushes the rational basis out of the natural class of the approximation space, which is exactly what the divergence result above quantifies.

### Numerical experiments: two kinds of source, two values of $\lambda$, one 2-D example

Common setup: errors are reported in the discrete $H^{\alpha/2}$-norm with $\rho=1$; only two values of $\lambda$ are tested, $\lambda=0$ (modified mapped Chebyshev rational functions) and $\lambda=0.5$ (modified mapped Legendre functions); $\alpha$ is $0.4$, $1$ and $1.6$ throughout.

| Example                                     | Data                                                                  | Rational-basis parameters         | Baseline                                      |
| ------------------------------------------- | --------------------------------------------------------------------- | --------------------------------- | --------------------------------------------- |
| Example 1 (exponentially decaying solution) | exact $u_e(x)=e^{-x^2}$, closed-form source                           | $\mu=5$                           | none (manufactured solution)                  |
| Example 1 (algebraically decaying solution) | exact $u_a(x)=(1+x^2)^{-r}$, $r=2.3$, closed-form source              | $\mu=3$                           | none (manufactured solution)                  |
| Example 2 (exponentially decaying source)   | $f(x)=e^{-x^2/2}(1+x)$, no closed-form solution                       | $\mu=5$, $\lambda\in\{0,0.5\}$    | Hermite functions (Mao-Shen), scaling $1/0.4$ |
| Example 3 (algebraically decaying source)   | $f(x)=(1+x^2)^{-2}$                                                   | $\mu=3$, both values of $\lambda$ | Hermite functions, scaling $1/0.7$            |
| Collocation: multi-term model               | sources $e^{-x^2/2}(1+x)$ and $(1+x^2)^{-1.8}$                        | not reported                      | Hermite functions                             |
| Two dimensions (§7)                         | $f(x,y)=e^{-\sqrt{x^2+y^2}}$, $\mathcal F[f]=(1+\xi^2+\eta^2)^{-3/2}$ | not reported                      | Hermite collocation                           |

Both sources in Example 1 are available in closed form, which is what makes it a benchmark:

$$
f_e(x)=\rho e^{-x^2}+\frac{2^{\alpha}\Gamma\bigl(\tfrac{\alpha+1}{2}\bigr)}{\Gamma\bigl(\tfrac12\bigr)}\
{}_1F_1\!\Bigl(\tfrac{\alpha+1}{2};\tfrac12;-x^2\Bigr),
$$

$$
f_a(x)=\rho(1+x^2)^{-r}
+\frac{2^{\alpha}\Gamma\bigl(\tfrac{\alpha}{2}+r\bigr)\Gamma\bigl(\tfrac{\alpha+1}{2}\bigr)}
{\Gamma(r)\,\Gamma\bigl(\tfrac12\bigr)}\
{}_2F_1\!\Bigl(\tfrac{\alpha}{2}+r,\ \tfrac{\alpha+1}{2};\ \tfrac12;\ -x^2\Bigr).
$$

The outcome is the most instructive part of the section: **even when the exact solution is as benign as $e^{-x^2}$, the error decays only algebraically**. The second term of Theorem 5.1 explains why — the source $f_e$ decays only like $(1+x^2)^{-\frac{\alpha+1}{2}}$ (the text available here drops the minus sign on that exponent, corrected here to match the decay it describes), so the error is dominated by the **interpolation error of the source**. Theorem 5.1 plus a direct calculation predicts $O\bigl(N^{-(\alpha+\frac12)}\bigr)$, which the paper reports agrees well with the numerics. Substituting the three values of $\alpha$ tested:

| $\alpha$ | Predicted rate |
| -------- | -------------- |
| $0.4$    | $N^{-0.9}$     |
| $1$      | $N^{-1.5}$     |
| $1.6$    | $N^{-2.1}$     |

Those three rows are only $N^{-(\alpha+1/2)}$ evaluated at the tested $\alpha$, not separately measured values. For $u_a$ with $r=2.3$ the rate is the same, because $f_a\sim(1+x^2)^{-\min(r,\frac{\alpha+1}{2})}$ and the minimum still sits on the $\frac{\alpha+1}{2}$ side.

Examples 2 and 3 are the head-to-head comparison against the Hermite approach. Example 2 has no closed-form solution, so the $N=600$ numerical solution serves as reference. Both give the same verdict: **the rational basis outperforms the Hermite approximation in every case, with markedly higher convergence rates**, tabulated in the paper's Table 1 (for $\alpha=1$) and Table 2. Those tables' numbers could not be recovered from the text available here, so this page does not reproduce the rate values.

Figure 4 of the paper reports the asymptotics of the computed solutions, and that is the explanatory group of results as well as the source of the fact stated at the top of this page: two sources with completely different decay produce solutions with the same tail, $|u(x)|\sim|x|^{-\alpha-1}$. The paper uses that single observation to explain "why MMGFs have a better performance than the Hermite functions".

The collocation scheme is tested on a multi-term (distributed-order type) model $\sum_{j=1}^{4}\rho_j(-\Delta)^{\alpha_j/2}u=f$ with

| $j$        | 1       | 2       | 3       | 4       |
| ---------- | ------- | ------- | ------- | ------- |
| $\alpha_j$ | $0$     | $0.5$   | $1.5$   | $2$     |
| $\rho_j$   | $\pi/6$ | $\pi/3$ | $\pi/3$ | $\pi/6$ |

The four weights are Simpson-type weights on $[0,2]$ scaled by $\pi/2$, consistent with a quadrature discretisation of a distributed-order model — **that last observation is this site's reading of the parameters, not a statement of the paper**. The conclusion matches the Galerkin case: the rational basis is much better than the Hermite approach in all cases. The paper also notes that collocation is the more practical choice for variable-coefficient and nonlinear problems.

The two-dimensional example $f(x,y)=e^{-\sqrt{x^2+y^2}}$ is chosen with care: it decays exponentially in $x$ while its Fourier transform $(1+\xi^2+\eta^2)^{-3/2}$ decays only algebraically, a clean illustration of the claim that a rational basis is the right tool in the transformed variable. The rational collocation method is more accurate and converges faster than the Hermite one.

The experiments establish three things: the rational basis beats the Hermite approach for both kinds of source, in both the Galerkin and the collocation scheme, and in one as well as two dimensions; the predicted $N^{-(\alpha+1/2)}$ rate is realised, so the term of Theorem 5.1 that actually binds is the source term; and the power-law tail $|x|^{-\alpha-1}$ is confirmed numerically, with the exponent set by the operator rather than the data.

They fall short of the theory in four places. First, **the convergence theorem covers Galerkin only, while the multi-term and two-dimensional experiments both use collocation**, so the most heavily tested scheme is the unanalysed one. Second, every example is bottlenecked by the regularity of the source, so the first term of Theorem 5.1 — the solution's own approximation error, carrying $N^{s-m}$ — is never really probed; the experiments confirm the weaker half of the theory. Third, only two values of $\lambda$ and three of $\alpha$ are tested, and the scaling $\mu$ is set by hand ($5$ or $3$) with no selection rule given — the same gap that appears in paper 27 as "$M$ is hard to find, so the optimal $r$ is hard to get". Fourth, the reference for Example 2 is the authors' own $N=600$ run rather than an independent exact solution.

## 92: uniformity rather than speed

Paper 92 closes a **uniformity** gap. The abstract states that existing methods do not work well as $\alpha$ approaches $0$ or $1$, and the structural reason is clear: the standard integral representations carry a prefactor $\sin(\alpha\pi)$ that vanishes at integer $\alpha$, and the strength of the integrand's endpoint behaviour is governed by $\alpha$. The Balakrishnan form and its negative-power counterpart are

$$
A^{\alpha}=\frac{\sin(\alpha\pi)}{\alpha\pi}\,A\int_0^{\infty}\bigl(t^{1/\alpha}I+A\bigr)^{-1}\mathrm dt,
\qquad
\mathcal L^{-\alpha}=\frac{2\sin(\alpha\pi)}{\pi}\int_0^{\infty}t^{2\alpha-1}
\bigl(\mathcal I+t^{2}\mathcal L\bigr)^{-1}\mathrm dt .
$$

In the best-uniform-rational-approximation family the error asymptotics carry the same factor: the best $(k,k)$ rational approximation to $t^{\beta-\alpha}$ on $[0,1]$ satisfies

$$
\lim_{k\to\infty}e^{2\pi\sqrt{(\beta-\alpha)k}}E_{\alpha}(k,k;\beta)
=4^{\,1+\beta-\alpha}\,\bigl|\sin\pi(\beta-\alpha)\bigr| ,
$$

root-exponential convergence in which $|\sin\pi(\beta-\alpha)|$ is exactly the degeneracy at $\alpha\to0,1$.

The paper targets cost and accuracy uniform in $\alpha\in(0,1)$ for two objects: $A^{-\alpha}$ with $A$ positive definite, and $(q\mathcal I+A^{\alpha})^{-1}$ with $q>0$. The second is the resolvent needed for an implicit time step of a space-fractional evolution equation, which explains the applications listed.

> [!warning] What could be verified
> The full text of this paper could not be reached through public channels and no preprint exists. What is confirmable: the journal keywords are trapezoidal rule, Laguerre-Gauss quadrature, fractional Laplacian, fractional Poisson equation and space-fractional Allen-Cahn equation, and the reference list includes the AAA rational approximation algorithm, double-exponential quadrature for fractional diffusion, the double-exponential formula for matrix fractional powers, and the exponentially convergent trapezoidal rule. That identifies the method as a **quadrature-based scheme** compared against double-exponential quadrature and rational approximation. This site does not report its theorems, convergence rates or the exact sense of "uniformly fast"; those require the published article.
>
> On the numerical side only existence can be reported: the abstract claims "sufficient numerical simulations" and "comparisons with the state-of-the-art methods", the keywords name three application problems (the fractional Laplacian, the fractional Poisson equation and the space-fractional Allen-Cahn equation), and the acknowledgments record that the right-hand part of one figure was produced by externally supplied code, so at least one comparison rests on someone else's implementation. **Matrix sizes and conditioning, the range of $\alpha$ actually tested, error levels and timings are all unverifiable, and this page lists no numerical-experiment table for this paper.** The two integral representations and the rational-approximation error asymptotics above come from other publicly checkable sources and are used only to exhibit the structural origin of the degeneracy; they must not be cited as this paper's formulas.

## Side-by-side comparison

| No. | Basis or representation           | Far-field decay                 | How the operator is handled                                    |
| --- | --------------------------------- | ------------------------------- | -------------------------------------------------------------- |
| 27  | Hermite and over-scaled functions | Gaussian                        | explicit physical-space differentiation matrices               |
| 42  | modified mapped Gegenbauer        | $\lvert x\rvert^{-(\lambda+1)}$ | decomposition into two elementary rational shapes              |
| 92  | no basis (a matrix function)      | not applicable                  | integral representation plus quadrature (limited verification) |

Paper 42 says explicitly that it works in the spirit of paper 27, and the two share one design pattern: derive closed forms for the fractional image of the basis, precompute the fractional differentiation matrices from them, and collocate directly in physical space. The multi-term and distributed-order application, the trick of filling the matrices cheaply through hypergeometric contiguous relations, and the idea of a scaling parameter are all inherited too. **The upgrade is only the basis**: $e^{-x^2/2}H_n$ becomes $(1+x^2)^{-(\lambda+1)/2}C^{\lambda}_n(x/\sqrt{1+x^2})$, an exponential tail becomes the power-law tail $\lvert x\rvert^{-(\lambda+1)}$, so that it can match the $\lvert x\rvert^{-\alpha-1}$ the fractional Laplacian actually produces.

Their theoretical status is complementary in the same way. Paper 27 has no convergence theorem; paper 42 proves optimal estimates for its Galerkin scheme — but **the collocation scheme is unanalysed in both**, and it is the scheme both papers actually promote. Their limitations mirror each other: paper 27 degrades to algebraic convergence on an algebraically decaying eigenproblem, and paper 42 converges algebraically at $N^{-(\alpha+1/2)}$ whenever the regularity of the source is the bottleneck.

Paper 92 is the discrete-algebra side of the same operator. Papers 27 and 42 make $(-\Delta)^{\alpha/2}$ computable by choosing basis functions with closed-form fractional images; paper 92 accepts an already discretised $A$ and computes $A^{-\alpha}$ and $(q\mathcal I+A^{\alpha})^{-1}$ directly. The latter is the matrix version of the Fourier multiplier $1/(|\xi|^{\alpha}+\rho)$ in paper 42's multidimensional scheme, and paper 92's listed applications (the fractional Laplacian, the fractional Poisson equation) are the discretised form of paper 42's model problem. The two routes are complementary rather than competing: **paper 42's route wins when you get to choose the discretisation, paper 92's when $A$ is handed to you.**

## Coverage check

| Item                                                | Paper  | Status                                                                                                  |
| --------------------------------------------------- | ------ | ------------------------------------------------------------------------------------------------------- |
| Decay-rate matching as the page's through-line      | all    | tail type of both families, operator-set tail exponent                                                  |
| Model problem and both operator definitions         | 27, 42 | hypersingular integral, constant, symbol form                                                           |
| Why all three standard routes fail                  | 42     | truncation, transparent conditions, exponential bases                                                   |
| The idea: what a closed-form matrix buys            | 27     | cost of dense and hypersingular, contrast with frequency space                                          |
| Model, both Hermite bases and orthogonality         | 27     | reaction term, normalised and over-scaled forms, Brinkman origin                                        |
| Fractional image of the over-scaled basis           | 27     | Fourier image, cosine integral, ${}_1F_1$ closed forms (typography flagged)                             |
| $F_k$ formulas and matrix for the normalised basis  | 27     | Fourier eigenfunctions, coefficient recursion, both closed forms, eq. (4.3)                             |
| Recurrence filling, scaling, multi-term, 2-D        | 27     | contiguous relation, $r=\max\lvert x_k\rvert/M$, matrices add, $x_i^2+y_j^2$ structure                  |
| No convergence theorem, condition-number results    | 27     | two sets of fitted slopes; slope $\approx\alpha$ is this site's reading                                 |
| Six numerical examples                              | 27     | full setups, Airy reference values; error magnitudes only in figures                                    |
| Map, modified basis and orthogonality               | 42     | map identities, absorbed weight, recurrence, far-field limit                                            |
| Pfaff expansion and the unverified constants        | 42     | both finite sums; $a^{\lambda}_n,b^{\lambda}_n$ flagged unverified                                      |
| Fractional Laplacian of the two shapes              | 42     | closed forms, the constant $A^{\gamma}_s$, Bessel proof route, series form                              |
| Divergence between integer and fractional order     | 42     | the three-case conclusion and its meaning                                                               |
| Fractional image of the basis (Theorem 3.4)         | 42     | both parities and the reachable degree                                                                  |
| Recurrence evaluation and scaling parameter         | 42     | contiguous relation, roles of $\mu$ and $\lambda$                                                       |
| Approximation space, nodes, both schemes            | 42     | equivalence, mapped Gauss data, matrices, open problem                                                  |
| $\mathbb B^m_{\lambda}$, the lemma, three estimates | 42     | full hypotheses and proof route; seminorm exponent reconstructed                                        |
| Theorem 5.1 and the multidimensional bound          | 42     | origin of both terms, which dominates, Remark 7.2's caveat                                              |
| Multidimensional Fourier algorithm, self-duality    | 42     | diagonalisation, algorithm steps, reason for the basis                                                  |
| Six groups of numerical experiments                 | 42     | setups, closed-form sources, predicted rate, power-law tail, baselines; rate-table values not recovered |
| Uniformity gap and $\sin(\alpha\pi)$                | 92     | both representations, rational error asymptotics, both targets                                          |
| What paper 92 leaves unreported                     | 92     | theorems, rates, matrix sizes, $\alpha$ range, errors, timings                                          |

## Sources for this page

- T. Tang, H. Yuan, and T. Zhou, [_Hermite spectral collocation methods for fractional PDEs in unbounded domains_](https://doi.org/10.4208/cicp.2018.hh80.12), Commun. Comput. Phys. 24(4) (2018), pp. 1143-1168 (preprint [arXiv:1801.09073](https://arxiv.org/abs/1801.09073)).
- T. Tang, L.-L. Wang, H. Yuan, and T. Zhou, [_Rational spectral methods for PDEs involving fractional Laplacian in unbounded domains_](https://doi.org/10.1137/19M1244299), SIAM J. Sci. Comput. 42(2) (2020), pp. A585-A611 (preprint [arXiv:1905.02476](https://arxiv.org/abs/1905.02476)).
- Y. Duan, F. Zeng, H. Zhang, and T. Zhou, [_Fast computation of the fractional power of a matrix_](https://doi.org/10.1137/25M1757411), SIAM J. Sci. Comput. 48(1) (2026), pp. A309-A334.
