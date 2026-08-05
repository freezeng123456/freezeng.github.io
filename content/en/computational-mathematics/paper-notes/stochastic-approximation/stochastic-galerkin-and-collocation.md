---
title: Stochastic Galerkin and Collocation
description: Papers 1, 2, 3, 4, 5, 7 and 38 - algebra of the coupled system and regularity in the random variable
lang: en
translation: computational-mathematics/paper-notes/stochastic-approximation/stochastic-galerkin-and-collocation
tags:
  - paper-notes
  - uncertainty-quantification
  - stochastic-galerkin
---

> [!note] Coverage of this page
> Papers **1** (_J. Comput. Phys._ 229, 2010), **2** (_Commun. Comput. Phys._ 8, 2010), **3** (_J. Comput. Appl. Math._ 236, 2011), **4** (_J. Comput. Math._ 30, 2012), **5** (_J. Sci. Comput._ 51, 2012), **7** (_Adv. Appl. Math. Mech._ 6, 2014) and **38** (_J. Sci. Comput._ 80, 2019).

These seven belong to the **intrusive** route: expand the solution in generalised polynomial chaos and either apply Galerkin projection or solve decoupled problems at parameter points. Two technical questions recur: the algebraic properties of the coupled system, and the regularity of the solution **in the random variable**.

## The mechanism first: why a random wave speed leaks accuracy as time grows

Most of this page revolves around one model problem: the scalar transport equation $\partial_t u=c(y)\partial_x u$ with a wave speed $c$ depending on a random parameter $y$. In the compatible case the solution translates along characteristics,

$$
u(x,t;y)=u_0\bigl(x+c(y)t\bigr),
$$

which is exactly the explicit solution paper 5 uses in its eq. (3.14) to estimate parametric sensitivity. Differentiating it in $y$, the chain rule emits a factor $t$ every time:

$$
\partial_y u=t\,c'(y)\,u_0'\bigl(x+c(y)t\bigr),
\qquad
\partial_y^k u=\bigl(t\,c'(y)\bigr)^k u_0^{(k)}\bigl(x+c(y)t\bigr)+\cdots,
$$

the omitted terms carrying lower powers of $t$ and higher derivatives of $c$. This formula is the physical intuition for the whole page: **two parameter values $y_1\ne y_2$ carry characteristics that separate at the speed difference $c(y_1)-c(y_2)$, so by time $t$ they sit a distance $|c(y_1)-c(y_2)|t$ apart**. Seen at a fixed $(x,t)$, therefore, the map $y\mapsto u(x,t;y)$ oscillates on a parameter scale of about $1/t$: the longer you integrate, the finer the structure a polynomial in the random variable has to resolve. A degree-$p$ polynomial cannot resolve a feature of scale $1/t$ until $p\gtrsim t$.

That mechanism explains three things that recur below:

1. every constant in every theorem on this page carries a $T$ ($C(T)$, $C_k(T)$, $C_\Sigma C(T)$) and **none is uniform in time** — the convergence analyses have to be time-dependent;
2. the experiments of papers 4 and 5 both report the error growing with $t$ at fixed discretisation, paper 5 describing the dependence as typically linear;
3. the same machinery does worse on hyperbolic problems than on elliptic or parabolic ones, because the latter have no time variable to amplify the parametric derivatives step after step. Paper 4 gives the sharp version: solutions of random **hyperbolic** equations are **not analytic in general** with respect to the random parameters.

> [!note] What this passage is
> The differentiation above is elementary calculus on the model problem these papers share, not a theorem from any one of them. The papers report the time growth as an observation and encode it in $T$-dependent constants; the reading "the parameter scale is about $1/t$" is offered here as the thread tying the results together. It can be made quantitative on paper 4's own example, as done below.

## 1: when the stochastic Galerkin coefficient matrix is strictly diagonally dominant

### The idea

Generalised polynomial chaos Galerkin projection replaces one random PDE by a **coupled** system of $M+1$ deterministic PDEs, and all of the coupling sits in a single matrix $A(x)$. To advance that system cheaply one normally splits $A$ into its diagonal and off-diagonal parts: mixed explicit-implicit time stepping with the diagonal implicit, Jacobi iteration, or conjugate gradients preconditioned by the diagonal blocks. The rate of any such splitting is governed by the size of $D^{-1}(A-D)$, and **strict diagonal dominance is precisely the condition that keeps that contraction factor below one**.

The point is that dominance holds pointwise: the conclusion is asserted for every $x\in\Omega$, so it survives spatial discretisation and the splitting solvers converge at a mesh-independent rate. Moreover the paper aims at the strengthened form

$$
a_{jj}\ \ge\ \kappa_{\min}+\sum_{k\ne j}|a_{jk}|,
$$

whose margin is the ellipticity constant $\kappa_{\min}$ rather than merely "positive" — a gap that does not degrade with $j$ or $x$, which is where mesh-independence comes from. Without it the splitting iteration can diverge, or its rate can deteriorate as the mesh is refined or the polynomial order $P$ is raised.

Xiu and Shen had proved this only for **symmetric** Beta densities ($\alpha=\beta$) and for Gamma densities, and explicitly left the general asymmetric Beta case ($\alpha\ne\beta$) open. The entire content of this six-page note is the answer to that question: no new discretisation, no new solver.

### Setting

The model problem is

$$
\partial_t u(x,y,t)=\nabla\cdot\bigl(\kappa(x,y)\nabla_x u(x,y,t)\bigr)+f(x,y,t),
\qquad x\in\Omega\subset\mathbb R^d\ (d=1,2,3),\ t\in(0,T],
$$

with $u(x,y,0)=u_0(x,y)$ and $u(\cdot,y,t)|_{\partial\Omega}=0$; the steady counterpart is $\nabla\cdot(\kappa\nabla_x u)=f$. The random field depends **affinely** on the parameters,

$$
\kappa(x,y)=\kappa_0(x)+\sum_{i=1}^{N}\kappa_i(x)\,y_i,
\qquad \kappa_0(x)>0,
\qquad \kappa(x,y)\ge\kappa_{\min}>0,
$$

with $y=(y_1,\dots,y_N)$ having i.i.d. components. The generalised polynomial chaos expansion is

$$
u\approx\sum_{m=0}^{M}v_m(x,t)\Phi_m(y),
\qquad
f\approx\sum_{m=0}^{M}f_m(x,t)\Phi_m(y),
\qquad
M=\binom{N+P}{N},
$$

with $\Phi_m(y)=\phi_{m_1}(y_1)\cdots\phi_{m_N}(y_N)$, $m_1+\cdots+m_N\le P$, and orthonormal univariate factors $\int\phi_j(y_i)\phi_k(y_i)\rho_i(y_i)\,\mathrm dy_i=\delta_{jk}$, so that $\mathbb E[\Phi_m\Phi_n]=\delta_{mn}$ with $\rho(y)=\prod_{i=1}^{N}\rho_i(y_i)$.

### Derivation

Substituting the expansion, multiplying by $\Phi_k$ and taking expectations gives the coupled system in component form,

$$
\partial_t v_k=\sum_{j=1}^{M}\nabla\cdot\bigl(a_{jk}(x)\nabla v_j\bigr)+f_k(x,t),
\qquad k=1,\dots,M,
$$

that is $\partial_t v=\nabla\cdot(A\nabla_x v)+f$ with $A(x)=(a_{jk})_{1\le j,k\le M}$ and

$$
a_{jk}=\sum_{i=0}^{N}\kappa_i(x)\,e_{ijk},
\qquad
e_{ijk}=\int y_i\,\Phi_j(y)\Phi_k(y)\rho(y)\,\mathrm dy .
$$

By construction $A=A^{T}$. With the convention $y_0\equiv1$ one has $e_{0jk}=\int\Phi_j\Phi_k\rho\,\mathrm dy=\delta_{jk}$, so the $i=0$ term only adds $\kappa_0(x)$ to the diagonal; all genuine coupling comes from the terms with $i\ge1$, that is from the operator **multiplication by $y_i$**.

The sparsity is read off the three-term recurrence. Writing the normalised recurrence in direction $i$ as

$$
y_i\,\phi_{m}(y_i)=a^i_{m}\,\phi_{m+1}(y_i)+b^i_{m}\,\phi_{m}(y_i)+c^i_{m}\,\phi_{m-1}(y_i),
$$

one gets, for a multi-index $j=(j_1,\dots,j_N)$,

$$
y_i\Phi_j=\bigl(a^i_{j_i}\phi_{j_i+1}+b^i_{j_i}\phi_{j_i}+c^i_{j_i}\phi_{j_i-1}\bigr)\prod_{l\ne i}\phi_{j_l},
$$

so $e_{ijk}$ is nonzero only when $k$ agrees with $j$ in every direction except possibly $i$, where $k_i\in\{j_i-1,j_i,j_i+1\}$. Taking the union over $i=1,\dots,N$ leaves at most $2N+1$ nonzeros per row (the diagonal counted once) — that is where the Xiu-Shen sparsity statement comes from.

The same expansion yields the two identities that carry the whole note:

$$
a_{jj}=\kappa_0(x)+\sum_{i=1}^{N}\kappa_i(x)\,b^i_j,
\qquad
\sum_{k\neq j}|a_{jk}|=\sum_{i=1}^{N}|\kappa_i(x)|\,|a^i_j+c^i_j| .
$$

The first follows from $e_{ijj}=b^i_{j_i}$. The second collects the two neighbour contributions $\kappa_i a^i_j$ and $\kappa_i c^i_j$ in direction $i$, and combining them into the single absolute value $|a^i_j+c^i_j|$ uses the sign structure of the normalised recurrence coefficients. Since all $y_i$ share a distribution, the superscript $i$ drops. **That step is the technical core of the note**: it trades a question about an $M\times M$ matrix for an elementary inequality on two or three recurrence coefficients.

What remains is Jacobi-polynomial bookkeeping. With $\varpi:=2n+\alpha+\beta$, orthogonality reads

$$
\int_{-1}^{1}(1-x)^{\alpha}(1+x)^{\beta}P^{\alpha,\beta}_mP^{\alpha,\beta}_n\,\mathrm dx
=\frac{2^{\alpha+\beta+1}}{2n+\alpha+\beta+1}\frac{\Gamma(n+\alpha+1)\Gamma(n+\beta+1)}{\Gamma(n+\alpha+\beta+1)\,n!}\,\delta_{nm}
=:h^{\alpha,\beta}_n\delta_{nm},
$$

and the unnormalised recurrence coefficients are

$$
a_n=\frac{2(n+1)(n+\alpha+\beta+1)}{(\varpi+1)(\varpi+2)},
\qquad
b_n=\frac{\beta^2-\alpha^2}{\varpi(\varpi+2)},
\qquad
c_n=\frac{2(n+\alpha)(n+\beta)}{\varpi(\varpi+1)} .
$$

Normalising by $\widetilde P^{\alpha,\beta}_n=P^{\alpha,\beta}_n/\sqrt{h^{\alpha,\beta}_n}$ gives $x\widetilde P_n=\tilde a_n\widetilde P_{n+1}+\tilde b_n\widetilde P_n+\tilde c_n\widetilde P_{n-1}$, where the numerator of $\tilde a_n$ is $2\sqrt{(n+1)(n+\alpha+1)(n+\beta+1)(n+\alpha+\beta+1)}$.

> [!warning] One formula that could not be verified
> The **denominator** of $\tilde a_n$ — a product of $\varpi$-shifts — was corrupted in the text available here, so its exact form cannot be confirmed and is not transcribed. The numerator is verified.

The expression for $b_n$ explains why the symmetric case fell first: $b_n=0$ exactly when $\alpha^2=\beta^2$. For $\alpha=\beta$ the diagonal entry is then exactly $\kappa_0(x)$ and dominance reduces to bounding the off-diagonal row sum by $\kappa_0-\kappa_{\min}$. For $\alpha\ne\beta$ the coefficient $b_n$ is nonzero and **lowers** the diagonal wherever $\kappa_i b_j<0$, which is precisely the new difficulty. The threshold $|\alpha|,|\beta|\ge1/2$ is the condition under which the note's estimate closes; the note does not claim it is sharp.

### Theorems

**Lemma 1 (quoted from Xiu-Shen)**: $A(x)$ is positive definite for every $x$, and each row has at most $2N+1$ nonzero entries.

**Lemma 2 (Xiu-Shen)**: if the $y_i$ have an identical **symmetric** Beta density on $(-1,1)$, $\rho(y_i)=(1-y_i)^{\alpha}(1+y_i)^{\alpha}$, or an identical Gamma density on $(0,+\infty)$, $\rho(y_i)=y_i^{\alpha}e^{-y_i}$ with $\alpha>-1$ (scaling constants omitted), then $A(x)$ is strictly diagonally dominant for all $x\in\Omega$ in the strengthened sense $a_{jj}\ge\kappa_{\min}+\sum_{k\ne j}|a_{jk}|$.

**Theorem 1 (the new result)**: suppose (I) the $y_i$ have an identical Beta density on $(-1,1)$, $\rho(y_i)=(1-y_i)^{\alpha}(1+y_i)^{\beta}$ with $\alpha,\beta>-1$ satisfying

$$
|\alpha|\ \ge\ \tfrac12,\qquad |\beta|\ \ge\ \tfrac12,
$$

or (II) the $y_i$ have an identical Gamma density on $(0,+\infty)$, $\rho(y_i)=y_i^{\alpha}e^{-y_i}$ with $\alpha>-1$ (scaling constants omitted). Then the matrices $A(x)$ built from the corresponding generalised polynomial chaos basis are strictly diagonally dominant for all $x\in\Omega$:

$$
a_{jj}\ \ge\ \kappa_{\min}+\sum_{k\ne j}|a_{jk}|,
\qquad 1\le j\le M,\ \forall x\in\Omega .
$$

The open question is thus answered affirmatively but **conditionally**: asymmetric Beta parameters are allowed at the price of the explicit threshold $|\alpha|,|\beta|\ge1/2$. The note does not claim dominance for all $\alpha\ne\beta$ with $\alpha,\beta>-1$, and whether the threshold can be removed is not settled there.

One convention deserves note: with $\rho(y_i)=(1-y_i)^\alpha(1+y_i)^\beta$ and $\alpha,\beta>-1$, the condition $|\alpha|\ge1/2$ admits both $\alpha\ge1/2$ and $-1<\alpha\le-1/2$; **the uniform (Legendre) case $\alpha=\beta=0$ is not covered by Theorem 1** but is covered by the symmetric case of Lemma 2.

### Numerical experiments

None. This is a purely analytical short note: its six pages are the problem restatement, the recurrence-coefficient algebra and the proof, with no numerical example anywhere.

### Relation to the others

This is the earliest and most narrowly algebraic entry in the group — it studies the **linear algebra produced by** stochastic Galerkin projection rather than sampling or approximation theory. It shares the intrusive-Galerkin theme with papers 3, 5 and 38 and pairs with papers 2 and 4, which analyse the non-intrusive collocation alternative. It is essentially disjoint from the later sampling-design line (the [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|least-squares page]] and the [[en/computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce|sparse-recovery page]]), which abandons Galerkin projection in favour of least squares and $\ell_1$ recovery from point evaluations.

## 2: convergence rates for collocation derived from hypotheses on the data

### The idea

For the scalar hyperbolic equation $\partial_t u=c(y)\partial_x u$, a random wave speed $c(y)$ that can **change sign** makes the inflow boundary switch ends of the domain: data are prescribed at $x=-1$ when $c(y)<0$ and at $x=1$ when $c(y)>0$. Near the parameter value where $c$ crosses zero the type of boundary condition changes, so **even with data smooth in both $x$ and $y$, the map $y\mapsto u(x,t;y)$ can have very low regularity**.

This hits the two method families differently. At a collocation node $y_j$ the speed $c(y_j)$ has a definite sign, so one simply imposes the well-posed boundary condition of that scalar problem; imposing boundary and initial conditions pointwise in $y$ is therefore trivial. Galerkin projection instead produces a symmetric hyperbolic system with eigenvalues of both signs, and boundary conditions must be imposed consistently on the characteristic variables — the difficulty Gottlieb and Xiu had to handle.

The second motivation concerns where the hypotheses sit. Gottlieb and Xiu's convergence proof in the Galerkin case assumed fast asymptotic decay of the expansion coefficients — an assumption about the **answer**, not about the **data**. This paper wants rates derived from $c$, $u_0$, $u_L$ and $u_R$.

### Setting

$$
\partial_t u(x,t;y)=c(y)\,\partial_x u(x,t;y),
\qquad x\in D\equiv(-1,1),\ t>0,
$$

with $u(x,0;y)=u_0(x;y)$ and the well-posed boundary conditions $u(-1,t;y)=u_L(t;y)$ when $c(y)<0$ and $u(1,t;y)=u_R(t;y)$ when $c(y)>0$. The random variable $y$ has density $\rho(y)$ on $\Gamma\equiv[-1,1]$.

The collocation scheme takes Gauss nodes $\{y_i\}_{i=0}^{N}$ (roots of $\Phi_{N+1}$; Hermite for Gaussian, Legendre for uniform, Laguerre for Gamma), solves $N+1$ **decoupled** deterministic problems, and interpolates:

$$
u^N(x,t;y)=\mathcal I^y_N u:=\sum_{k=0}^{N}u(x,t;y_k)\,F_k(y),
\qquad F_k\in\mathbb P_N,\ F_i(y_k)=\delta_{ik}.
$$

Two error functionals are used: the mean-square error $e_{ms}(u-u^N):=M[u-u^N]$ and the mean error $e_{mean}(u-u^N):=\mathbb E[|u-u^N|]$. Regularity in the random direction is measured in $H^1$, $H^2$ and $BV$, using the tensor-space isomorphism $L^2\otimes H^k(D)\simeq L^2(\Gamma;H^k(D))\simeq H^k(D;L^2(\Gamma))$.

The technical device is a set of weighted data hypotheses (eq. (2.8)) that divide the density by the wave speed:

$$
\int_\Gamma\!\!\int_D \rho(y)\bigl(\partial_y u_0(x;y)\bigr)^2\mathrm dx\,\mathrm dy<\infty,
$$

$$
\int_0^T\!\!\int_{\Gamma^+}\frac{\rho(y)}{c(y)}\bigl(\partial_y u_R(t;y)\bigr)^2\mathrm dy\,\mathrm dt<\infty,
\qquad
\int_0^T\!\!\int_{\Gamma^-}\frac{\rho(y)}{|c(y)|}\bigl(\partial_y u_L(t;y)\bigr)^2\mathrm dy\,\mathrm dt<\infty,
$$

where $\Gamma^{\pm}$ are the subsets on which $c$ is positive and negative. **The $1/|c|$ weight is exactly what couples the wave speed to the admissible boundary data**: parameter regions where the speed is near zero demand smoother boundary data there.

### Derivation

Regularity comes from an energy estimate after differentiating in $y$. With $w=\partial_y u$, differentiating the equation gives

$$
\partial_t w=c(y)\,\partial_x w+c'(y)\,\partial_x u .
$$

Multiplying by $w$ and integrating over $D$ turns the transport term into a boundary flux:

$$
\frac{\mathrm d}{\mathrm dt}\frac12\int_D w^2\,\mathrm dx
=\frac{c(y)}{2}\Bigl[w^2\Bigr]_{x=-1}^{x=1}
+\int_D c'(y)\,\partial_x u\,w\,\mathrm dx .
$$

When $c(y)>0$ the inflow end is $x=1$, where $w(1,t;y)=\partial_y u_R(t;y)$ is given data while the outflow end enters with the favourable sign; when $c(y)<0$ the ends swap. Multiplying by $\rho(y)$ and integrating over $\Gamma$ and $[0,T]$, the boundary contributions land on the quantities that (2.8b-c) are there to control; the cross term $\int c'\,\partial_x u\,w$ is handled by Cauchy-Schwarz, which needs $|c'(y)|\le C$ and an a priori bound on $\partial_x u$ (besides (2.8a-c) the theorem's hypotheses also invoke the earlier finiteness conditions (2.1)), and Gronwall closes the estimate, which is where $C(T)$ comes from. The second derivative repeats the same argument.

> [!note] On where the $1/|c|$ weight comes from
> The weight in (2.8b-c) is $\rho/|c|$, whereas the boundary flux in the energy identity above carries $c$ itself — a discrepancy of $c^2$. A natural reconciliation: at the inflow end the equation itself gives $\partial_x u=\partial_t u/c$, so trading the time derivative on the boundary for the space derivative costs exactly one factor $1/c$. Paper 4's hypotheses contain $d_R:=\partial_t u_R/c(y)$ and $d_L:=\partial_t u_L/c(y)$ outright, consistent with that exchange. The full derivation at this point was not available here; the observation is recorded only as a consistency check, and (2.8) itself is transcribed as printed.

The rates follow by combining that regularity with a standard interpolation estimate. **Lemma 3.1** (cited to Canuto et al., p. 289): for $w^{(m)}\in L^2(-1,1)$ and its interpolant $\mathcal I_N w$ at $N+1$ Gauss, Gauss-Radau or Gauss-Lobatto points,

$$
\|w-\mathcal I_N w\|_{L^2(D)}\le C\,N^{-m}\,\|w^{(m)}\|_{L^2(-1,1)},
\qquad m\le N .
$$

Taking $m$ to be the regularity index in the random direction gives Theorems 3.1 and 3.2: **the order of convergence equals the order of regularity**.

The $BV$ case travels a different road. No interpolation estimate is available, so the proof uses the Gauss-quadrature remainder for functions of bounded variation (**Lemma 3.2**, quoted from the literature): if $f$ has total variation $V(f)$ on $[-1,1]$, the $N$-point Gauss quadrature remainder satisfies

$$
|R_N(f)|\le\frac{\pi}{2N+1}V(f).
$$

The mean error $\mathbb E[|u-u^N|]$ is itself an integral against $\rho$ and can be estimated by that remainder directly, whereas the mean-square error needs pointwise control that $BV$ does not supply. **That is why Theorem 3.3 controls the mean and not the mean square.**

### Theorems

**Theorem 2.1 ($H^1$ regularity)**: under $|c'(y)|\le C$ a.e. in $\Gamma$ (boundedness in the distribution sense) plus the finiteness hypotheses (2.1) and (2.8a-c), the solution has bounded $H^1$ regularity in $y$, uniformly for $0<t\le T$.

**Theorem 2.2 ($H^2$ regularity)**: if in addition $c''(y)$ is bounded in the distribution sense and the analogous second-derivative data hypotheses hold, the solution has $H^2$ regularity in $y$.

**Theorem 2.3 ($BV$ regularity)**: under $|c'(y)|\le C$ and bounded-variation hypotheses on the data (eq. (2.30a) onward, e.g. $\int_\Gamma\int_D\rho(y)|\partial_x u_0(x;y)|\,\mathrm dx\,\mathrm dy<+\infty$), the solution is of bounded variation in $y$.

**Theorem 3.1**: under the hypotheses of Theorem 2.1 ($H^1$ in $y$),

$$
e_{ms}(u-u^N)\le C(T)\,N^{-1},
\qquad
e_{mean}(u-u^N)\le C(T)\,N^{-1},
\qquad 0<t\le T,
$$

with $C(T)$ depending on $T$ but independent of $N$.

**Theorem 3.2**: under the hypotheses of Theorem 2.2 ($H^2$ in $y$), $e_{ms}\le C(T)N^{-2}$ and $e_{mean}\le C(T)N^{-2}$ for $0<t\le T$.

**Theorem 3.3**: under the hypotheses of Theorem 2.3 ($BV$ in $y$), **only the mean error** is controlled: $e_{mean}=\mathbb E[|u-u^N|]\le C(T)N^{-1}$ for $0<t\le T$.

The overall message is that **the rate is governed by the regularity of the solution in $y$, and that regularity is jointly determined by the random wave speed and the initial and boundary data — smooth data alone do not buy it.**

### Numerical experiments

The experiments take $y$ uniform on $\Gamma$ with Legendre-Gauss collocation points. Section 4.1 solves $u_t=y\,u_x$ with three initial conditions $u(x,0;y)=\sin(x)+4\,\mathrm{sgn}(y)\,y^{k}$, increasing $k$ raising the regularity in $y$:

| The term $4\,\mathrm{sgn}(y)y^{k}$ | Its shape           | Regularity in $y$ | Basis                     | Predicted order |
| ---------------------------------- | ------------------- | ----------------- | ------------------------- | --------------- |
| $k=1$                              | $\lvert y\rvert$    | $H^1$             | Theorem 2.1 + Theorem 3.1 | $N^{-1}$        |
| $k=2$                              | $y\lvert y\rvert$   | $H^2$             | Theorem 2.2 + Theorem 3.2 | $N^{-2}$        |
| $k=3$                              | $y^2\lvert y\rvert$ | $H^3$             | Lemma 3.1 with $m=3$      | $N^{-3}$        |

The paper reports that the observed order tracks the regularity index. (The theorems it states stop at $H^2$; the predicted order in the $H^3$ row comes from Lemma 3.1 for general $m$, not from a separately stated theorem. Pointwise error values were not available here, so only the orders the paper states are reported.)

A further example takes data satisfying the $BV$ hypothesis (2.32a) but **not** the $H^1$ hypothesis (2.8a):

| Error      | Theory                                | Observed in the paper           |
| ---------- | ------------------------------------- | ------------------------------- |
| $e_{mean}$ | $O(N^{-1})$ by Theorem 3.3            | first order, matching the theory |
| $e_{ms}$   | Theorem 3.1 does not apply, no prediction | **not** first order, roughly half order |

This is the point of the experiment: **the hypotheses of Theorem 3.1 are sharp**. Drop the $H^1$ hypothesis and the mean error that Theorem 3.3 protects stays first order while the mean-square error immediately falls away. The paper also investigates an accuracy-enhancement technique following the multi-element collocation method of Foo, Wan and Karniadakis (_J. Comput. Phys._ 227 (2008) 9572-9595).

### Relation to the others

This is the collocation half of a pair with paper 4, which redoes the analysis for the same transport model in a spectral/Galerkin setting. Paper 5 builds a Galerkin scheme with bi-orthogonal polynomials to sidestep the coupling that makes Galerkin awkward here. It shares with paper 1 the Gottlieb-Xiu / Xiu-Shen intrusive-methods context. Methodologically it belongs to the "structured nodes, one-dimensional parameter" stage that the [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|least-squares page]] and the [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|optimal-sampling page]] leave behind in favour of random or designed multivariate samples.

## 4: analytic regularity, spectral convergence, and why it degrades with time

### The idea

Papers 2 and 5 give finite-order regularity results ($H^k$, $BV$) with correspondingly **algebraic** rates, while the experiments show exponential convergence that the theory cannot explain. The gap is clear: exponential convergence does not need "many derivatives" but "derivatives growing at most geometrically", which is exactly analyticity. So the route here is to assume geometric growth of the data derivatives, propagate it along the equation to the solution, sum the Taylor series into an analytic extension to a complex neighbourhood, and then invoke a classical best-approximation estimate. The paper states explicitly that it deliberately ignores the deterministic-solver error and studies only the random-space discretisation error.

### Setting

The model is the same as in paper 2, with $V=L^2(D)$. **Assumption 2.1** (eqs. (2.7a)-(2.7d)) imposes geometric-growth bounds on all $y$-derivatives of the data. With $d_R:=\partial_t u_R/c(y)$ and $d_L:=\partial_t u_L/c(y)$, for every $k\in\mathbb N$:

$$
\max_{y\in\Gamma}|\partial_y^k c(y)|\le\gamma^k;
\qquad
\max_{\Gamma\otimes T}|\partial_y^k d_R|\le\delta_R^k,\quad
\max_{\Gamma\otimes T}|\partial_y^k u_R|\le\delta_R^k;
$$

$$
\max_{\Gamma\otimes T}|\partial_y^k d_L|\le\delta_L^k,\quad
\max_{\Gamma\otimes T}|\partial_y^k u_L|\le\delta_L^k;
\qquad
\max_{y\in\Gamma}\|\partial_y^k u_0'\|_V^2\le\eta^k,
\quad
\max_{y\in\Gamma}\|\partial_y^k u_0\|_V^2\le\eta^k,
$$

where $u'=\partial_x u$ and, without loss of generality, $\gamma\ge\max\{\delta_R,\delta_L,\eta\}$. **The constants appear raised to the power $k$** — that geometric structure in $k$ is what makes the derivative bounds summable into an analytic extension.

The quantity $d_R$ is worth a remark: at the inflow end $x=1$ the equation itself gives $\partial_x u=\partial_t u/c$, so $d_R=\partial_t u_R/c$ is the spatial derivative of the solution on that boundary. That also explains why $c$ sits in a denominator, and matches the position of the $1/|c|$ weight in paper 2.

### Derivation

Differentiating the equation $k$ times in $y$, the Leibniz rule gives

$$
\partial_t\bigl(\partial_y^k u\bigr)
=c(y)\,\partial_x\partial_y^k u
+\sum_{l=1}^{k}\binom{k}{l}\bigl(\partial_y^l c\bigr)\,\partial_x\partial_y^{k-l}u .
$$

The first term is the original transport operator acting on $\partial_y^k u$; every term in the sum involves lower-order $y$-derivatives, controlled by the induction hypothesis, with coefficients bounded through $|\partial_y^l c|\le\gamma^l$. Running the energy estimate as in paper 2, the boundary fluxes are given by $\partial_y^k u_R$ and $\partial_y^k u_L$ and are bounded by $\delta_R^k$ and $\delta_L^k$, the initial term by $\eta^k$, and Gronwall closes in time and produces the constant $C_k(T)$ depending on $k$, $\gamma$ and $T$. The paper's proof is exactly this induction on $k$.

Given that family of derivative bounds, the analytic extension is built from the power series

$$
u(z,x,t)=\sum_{k=0}^{\infty}(z-y)^k\,\partial_y^k u(y,x,t)
$$

with region of convergence

$$
\Sigma(\Gamma,\tau)\equiv\{z\in\mathbb C:\ \mathrm{dist}(z,\Gamma)\le\tau\},
\qquad 0<\tau<1/\sqrt{\zeta}.
$$

Two things match up here. Theorem 2.1 controls the **square** of the derivative, $\|\partial_y^k u\|_V^2\lesssim\zeta^{k}$, hence $\|\partial_y^k u\|_V\lesssim\zeta^{k/2}$, so the series has ratio $\sqrt{\zeta}\,|z-y|$ and radius of convergence exactly $1/\sqrt{\zeta}$ — **that is where the condition $\tau<1/\sqrt\zeta$ in Theorem 2.2 comes from**. It also shows that writing the series without a $1/k!$ is not a slip: with $1/k!$ the geometric derivative bounds would give convergence everywhere and the finite radius $1/\sqrt\zeta$ would have no meaning.

> [!note] Two qualifications
> The series is transcribed as printed, without the $1/k!$ that the usual Taylor form carries. The constant $\zeta$ is introduced earlier in the paper; its exact definition could not be verified here, only that it is a data-dependent constant coming from Assumption 2.1.

### Theorems

**Theorem 2.1**: under Assumption 2.1, for every $k\in\mathbb N$,

$$
\max_{\Gamma}\|\partial_y^k u(\cdot,t,\cdot)\|_V^2
\ \le\ C_k(T)\bigl(\delta_R^k+\delta_L^k+\eta^k\bigr)\ <\ +\infty,
$$

with $C_k(T)$ depending on $k$, $\gamma$ and $T$ (proof by induction on $k$). The paper notes the consequence $u\in L^\infty[T,C^0(\Gamma,V)]$.

**Theorem 2.2**: the solution as a function of $y$, $u:\Gamma\to L^\infty(T,V)$, admits an analytic extension $u(x,t;z)$, $z\in\mathbb C$, in the region $\Sigma(\Gamma,\tau)$ with $0<\tau<1/\sqrt\zeta$.

**Theorem 3.1 (stochastic collocation; the paper attributes the proof to the literature and omits it)**: for $u\in L^\infty[T,C^0(\Gamma,V)]$ admitting an analytic extension in $\Sigma(\Gamma,\tau)$ for some $\tau>0$,

$$
\min_{v\in L^1[T,\mathbb P_p(\Gamma)\otimes V]}\|u-v\|_{L^\infty[T,C^0(\Gamma,V)]}
\ \le\ \frac{2}{\varrho-1}\,e^{-p\log\varrho}
\max_{z\in\Sigma(\Gamma,\tau)}\|u(z)\|_{L^\infty(T,V)},
\qquad
\varrho=\tau+\sqrt{1+\tau^2}>1 ,
$$

exponential decay in the polynomial degree $p$ at a rate $\log\varrho$ set by the width $\tau$ of the analyticity strip.

**Lemma 3.3 (projection error, quoted from the spectral-methods literature)**: for $u\in H^m(I)$, $m\ge1$, and $N$ sufficiently large, $\|u-\Pi_N u\|_{L^2(I)}\le C_I N^{-m}|u|_{H^m(I)}$ with $|u|_{H^m(I)}=\|\partial^m u\|_{L^2(I)}$ and $C_I$ independent of $N$.

**Theorem 3.2 (stochastic Galerkin, algebraic form)**:

$$
\|u-u^N_{SG}(\cdot,t,\cdot)\|_{L^2(\Gamma,V)}
\ \le\ C_\Gamma\sqrt{C(T)}\,\bigl(\sqrt{\zeta}\,N\bigr)^{-m},
$$

where $m$ is an integer index tied to the solution's regularity in the random space.

**Lemma 3.4 (quoted from paper 2)**: for any finite time $t$, the generalised polynomial chaos error satisfies $\mathbb E\|u-u^p\|_2^2\le C(T)\sum_{k=p+1}^{\infty}\|\tilde u_k\|_1^2$ with $u^p=\sum_{k=0}^{p}\tilde u_kP_k(y)$ and $\|u\|_1^2=\int_D(u^2+u_x^2)\,\mathrm dx$.

**Lemma 3.5 (periodic special case, $c(z)$ and $u_0$ analytic)**:

$$
\|\tilde u_n\|_1\ \le\ \frac{C_\Sigma}{2^n}\sqrt{\frac{2n+1}{2}}\int_{-1}^{1}\left(\frac{1-y^2}{1-|y|+\tau}\right)^{n}\mathrm dy .
$$

**Theorem 3.3 (spectral convergence of generalised polynomial chaos, as printed)**:

$$
\bigl(\mathbb E\|u-u^p\|_2^2\bigr)^{1/2}
\ \le\ \sqrt{C_\Sigma C(T)}\left(\frac{\sqrt{\pi}}{\sqrt{1-r^2}}+O\!\left(\frac{1}{p^{1/3}}\right)\right)\frac{r^{p+1}}{\sqrt{1-r^2}},
\qquad \xi=-1-\tau<-1 ,
$$

the geometric factor $r^{p+1}$ with $0<r<1$ being the spectral rate. The closed form of $r$ was corrupted in the text available here beyond $r\equiv1/(|\xi|+\cdots)$ with $\xi=-1-\tau$; the remainder is **unverified**, though its role is clear — a Bernstein-ellipse-type parameter determined by the analyticity width $\tau$.

> [!warning] An important qualification
> The paper states that unlike random elliptic or parabolic problems, solutions of random **hyperbolic** equations are **not analytic in general** with respect to the random parameters. The complex-analytic sharpening is therefore available only in special cases such as periodic boundary conditions with analytic data. The paper also calls its own stochastic Galerkin algebraic estimate
> $$\|u-u^N_{SG}\|_{L^2(\Gamma,V)}\le C_\Gamma\sqrt{C(T)}\bigl(\sqrt{\zeta}N\bigr)^{-m}$$
> "rather rough". Both qualifications belong with any citation of the exponential rate.

### Numerical experiments

Section 4 has a single example, reused from paper 2:

| Item        | Setting                                                                 |
| ----------- | ----------------------------------------------------------------------- |
| Equation    | $u_t(x,t;y)=y\,u_x(x,t;y)$, $0<x<2\pi$, $t>0$                           |
| Initial data | $u(x,0;y)=\cos(y)$ (transcribed as printed, see the note below)        |
| Boundary    | chosen so that the exact solution is $u(x,t;y)=\cos(x-yt)$              |
| Swept       | number of projection terms (Galerkin) / number of collocation nodes     |
| Recorded    | mean-square error against the swept quantity (Figs. 4.1, 4.2), several time levels |
| Observed    | exponential convergence for both families; errors grow with $t$         |

This solution lies in $H^{(m)}_y(-1,1)$ for every positive integer $m$, so exponential convergence is expected and is observed for both Galerkin and collocation. The growth of the error with $t$ is the known long-time degradation of polynomial chaos and stochastic collocation.

> [!warning] An inconsistency in the example's data
> The three pieces of information transcribed above are mutually incompatible: $u(x,0;y)=\cos(y)$ does not match $u=\cos(x-yt)$ at $t=0$ (the latter gives $\cos(x)$), and $\cos(x-yt)$ solves $u_t=-y\,u_x$, not $u_t=y\,u_x$. The plausible repair is $u_0=\cos(x)$ with the sign convention made consistent ($u_t=y\,u_x$ pairs with $u=\cos(x+yt)$). It is transcribed here as available and flagged rather than silently corrected; under either repair the point of the example is unchanged — the solution lies in $H^{(m)}_y$ for every order, so exponential convergence is expected and observed.

This example also makes the time mechanism from the top of the page quantitative. For $u=\cos(x\mp yt)$ one has $|\partial_y^k u|\le t^k$, so the power series $\sum(z-y)^k\partial_y^k u$ converges for $|z-y|<1/t$: **the analyticity width $\tau$ available shrinks like $1/t$**. For small $\tau$, $\varrho=\tau+\sqrt{1+\tau^2}\approx1+\tau$, so $\log\varrho\approx\tau\approx1/t$ and the factor $e^{-p\log\varrho}$ in Theorem 3.1 becomes roughly $e^{-p/t}$: to hold the error at a fixed level the polynomial degree must grow in proportion to $t$. The equivalent reading puts the effect in the prefactor — continued to a strip of width $\tau$, $|\cos(x\mp zt)|$ reaches order $\cosh(\tau t)$, so the rate $\log\varrho$ is time-independent while $\max_{z\in\Sigma}\|u(z)\|$ grows with $t$. Both readings agree with the error growth the paper reports. This quantification is a computation carried out here on the paper's own example; the paper reports only the qualitative observation.

### Relation to the others

This is the theoretical capstone of the transport-equation trilogy 2, 4, 5: paper 2 supplies $H^1$/$H^2$/$BV$ regularity and the matching $N^{-1}$/$N^{-2}$ algebraic rates for collocation, paper 5 supplies the decoupled Galerkin scheme, and this paper supplies analytic regularity and spectral rates covering **both** families, the paper itself describing it as a theoretical complement to those numerical papers. Its message that the achievable rate is dictated by parametric regularity becomes a different question in the later [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|least-squares]] and [[en/computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce|sparse-recovery]] papers: how many samples are needed, rather than what rate is attainable.

## 5: bi-orthogonal polynomials decouple the Galerkin system exactly

### The idea

Paper 1 asks when the coefficient matrix $A$ is diagonally dominant. Paper 5 asks a better question: **can one change basis so that it is simply diagonal?**

This works because the coupling has a single source. The wave field depends affinely on the parameters, so the only coupling operator in the Galerkin system is multiplication by $y_i$. An ordinary orthogonal basis diagonalises only the $L^2_\rho$ inner product — it makes the mass matrix the identity and does nothing at all about multiplication by $y_i$. A **bi-orthogonal basis diagonalises both quadratic forms at once**, so the coupling term becomes diagonal too and the whole Galerkin system falls apart into unrelated scalar transport problems.

It is worth naming the matrix being diagonalised: its entries are $\int y_i\phi_j\phi_k\rho\,\mathrm dy_i$, the very same object as $e_{ijk}$ in paper 1, that is the Jacobi matrix of the three-term recurrence whose entries are paper 1's $a^i_j$, $b^i_j$ and $c^i_j$. **The two papers handle literally the same matrix: paper 1 bounds its off-diagonal part, paper 5 diagonalises it away.**

The cost of not doing so is real: Galerkin projection couples the modes; unlike collocation, boundary conditions cannot be imposed pointwise (the coupled system is a symmetric hyperbolic system with eigenvalues of both signs); and if the random field is a truncated Karhunen-Loève expansion in $N$ variables, an isotropic tensor basis has size exponential in $N$. Add that hyperbolic equations have the worst regularity of the three PDE types, so the parametric solution may be only $BV$ or in a low-order $H^k$, and one has the threefold difficulty the paper attacks head on.

### Setting

$$
\partial_t u(x,t,y)=c(x,y)\,\partial_x u(x,t,y),\qquad x\in D\equiv[-1,1],
$$

with $u(1,t;y)=u_R(t,y)$ where $c(1,y)>0$, $u(-1,t;y)=u_L(t,y)$ where $c(-1,y)<0$, and $u(x,0;y)=u_0(x,y)$. The crucial generalisation is that the wave field $c$ (equivalently written $\kappa$) depends on both $x$ and a **vector** $y=(y_1,\dots,y_N)$ — a genuine random field, not just a random variable as in paper 2.

(On notation: the number of random variables is written $N$ in the model problem and $M$ in the basis construction; they denote the same count.)

The **bi-orthogonal (double-orthogonal) basis** is the central object. For each direction $i$, the univariate basis $\{\phi_{j,i}\}_{j=0}^{r_i}$ of $\mathbb P_{r_i}$ must satisfy **two** orthogonality relations simultaneously:

$$
\int_{\Gamma_i}\rho(y_i)\phi_{j,i}\phi_{k,i}\,\mathrm dy_i=\delta_{jk},
\qquad
\int_{\Gamma_i} y_i\,\rho(y_i)\phi_{j,i}\phi_{k,i}\,\mathrm dy_i=C_{k,i}\,\delta_{jk},
$$

where the $\{C_{k,i}\}_{k=0}^{r_i}$ are nonzero constants. The paper attributes the space to the literature and notes that constructing such a basis amounts to solving **eigenproblems**, at a cost negligible against solving the coupled system when the $r_j$ are not large. The $C_{k,i}$ are the eigenvalues of multiplication by $y_i$, equivalently Gauss quadrature nodes, which is why Remark 3.2 states that the resulting Galerkin method is **equivalent to a certain type of collocation method**, while retaining a Galerkin framework better suited to analysis.

The multivariate basis is a tensor product: with $\mathbb P_{\mathbf r}=\mathbb P_{r_1}\otimes\cdots\otimes\mathbb P_{r_M}\subset L^2(\Gamma,\rho)$ and multi-index $\mathbf i=(i_1,\dots,i_M)\le\mathbf r$,

$$
\Phi_{\mathbf i}(y)=\prod_{k=1}^{M}\phi_{i_k,k}(y_k),
\qquad
\int_\Gamma y_k\,\rho(y)\,\Phi_{\mathbf i}(y)\Phi_{\mathbf j}(y)\,\mathrm dy=C_{i_k,k}\,\delta_{\mathbf i\mathbf j},
$$

giving $N_y=\prod_{i=1}^{M}(r_i+1)$ basis functions.

### Derivation

Write $u=\sum_{\mathbf i\le\mathbf r}u_{\mathbf i}\Phi_{\mathbf i}(y)$, substitute the affine wave field and project against $\Phi_{\mathbf j}$. The $\kappa_0$ term yields $\delta_{\mathbf i\mathbf j}$ by the first orthogonality and the $\kappa_k y_k$ term yields $C_{i_k,k}\delta_{\mathbf i\mathbf j}$ by the second — **both diagonal** — so one is left with $N_y$ mutually independent deterministic equations:

$$
\frac{\partial u_{\mathbf i}}{\partial t}=\kappa_{\mathbf i}(x)\frac{\partial u_{\mathbf i}}{\partial x},
\qquad
u_{\mathbf i}(x,0)=u_{0\mathbf i}(x),
\qquad
\kappa_{\mathbf i}(x)=\kappa_0(x)+\sum_{k=1}^{M}\kappa_k(x)\,C_{i_k}.
$$

The coupled Galerkin system collapses into a family of scalar transport problems, each with its own **effective wave speed** $\kappa_{\mathbf i}(x)$. Boundary conditions can then be imposed mode by mode: $u_{\mathbf i}(1,t)=u_{R\mathbf i}(t)$ when $\kappa_{\mathbf i}(1)>0$, and $u_{\mathbf i}(-1,t)=u_{L\mathbf i}(t)$ when $\kappa_{\mathbf i}(-1)<0$.

This step requires the sign tests on the effective speeds to be meaningful, and the paper's reason is clean. Setting $j=k$ in the second defining relation gives

$$
C_{k,i}=\int_{\Gamma_i}y_i\,\rho(y_i)\,\phi_{k,i}^2\,\mathrm dy_i,
$$

while the first defining relation says $\rho\phi_{k,i}^2$ is a probability density, integrating to $1$. **So $C_{k,i}$ is the mean of $y_i$ under some probability density and must lie between the lower and upper bounds of $y_i$, that is inside the convex hull of its support.** The effective speeds therefore stay within the range of the original speed and the sign tests match the original problem.

The statistics come out for free:

$$
\mathbb E(u^M)=\sum_{\mathbf i\le\mathbf r}u_{\mathbf i}(x,t)\int_\Gamma\rho(y)\Phi_{\mathbf i}\,\mathrm dy,
\qquad
\int_\Gamma\rho(y)\Bigl(\sum_{\mathbf i\le\mathbf r}u_{\mathbf i}\Phi_{\mathbf i}\Bigr)^2\mathrm dy=\sum_{\mathbf i\le\mathbf r}\bigl(u_{\mathbf i}(x,t)\bigr)^2 .
$$

What remains is the size of the tensor basis, and the paper trades it for a sensitivity-driven anisotropic choice of orders. Karhunen-Loève eigenvalues typically decay as $\lambda_i\sim i^{-2m}$, $m\ge1$, with the decay index set by the covariance function; and for a compatible problem with constant $\kappa_i$ the exact solution $u(x,y,t)=u_0(x+\kappa(y)t)$ gives $|\partial u/\partial y_i|\sim i^{-m}$ and $|\partial^2u/\partial y_i\partial y_j|\sim i^{-m}j^{-m}$ (eq. (3.14)). The leading error contributions therefore come from the first derivatives in the low-index directions. **The rule is: relatively large $r_i$ for small $i$, small $r_i$ for large $i$, subject to a total-degree constraint $\sum_{i=1}^{M}r_i\le P$.** This pushes $N_y$ well below the isotropic $r^M$.

Remark 3.1 adds that the construction extends directly to several space dimensions, for instance $\partial_t u=\kappa_1(x_1,x_2,y)\partial_{x_1}u+\kappa_2(x_1,x_2,y)\partial_{x_2}u$ with independent random wave fields, by building the basis in each direction and multiplying; extension to hyperbolic **systems** is also stated to be straightforward.

### Theorems

**Lemma 3.1**: let $\kappa(x,y)=\kappa_0(x)+\sum_{i=1}^{N}\sqrt{\lambda_i}\,\kappa_i(x)y_i$ satisfy $0<\kappa_{\min}\le\kappa(x,y)<\kappa_{\max}$ and $|\partial_x\kappa(x,y)|<\bar\kappa_{\max}$. If

$$
\int_\Gamma\!\!\int_D\rho(y)\bigl(\partial_x u_0(x;y)\bigr)^2\mathrm dx\,\mathrm dy<\infty,
\qquad
\int_0^T\!\!\int_\Gamma\rho(y)\bigl(\partial_t u_R(t;y)\bigr)^2\mathrm dy\,\mathrm dt<\infty
$$

and the analogous condition for $u_L$ hold, then $\int_\Gamma\int_D\rho(y)u_x^2\,\mathrm dx\,\mathrm dy\le C(T)<\infty$ for $0<t\le T$.

**Theorem 3.1 (the sensitivity estimate, eq. (3.17))**: with deterministic $u_0$ and boundary data in Karhunen-Loève form

$$
u_R(y,t)=u_{R0}(t)+\sum_{i=1}^{N}\sqrt{\mu_i}\,u_{Ri}(t)y_i,
\qquad
u_L(y,t)=u_{L0}(t)+\sum_{i=1}^{N}\sqrt{\nu_i}\,u_{Li}(t)y_i,
$$

where $\int_0^T u_{Ri}^2\,\mathrm dt<\infty$, and under assumption (3.15),

$$
\int_\Gamma\!\!\int_D \rho(y)\,u_{y_i}^2\,\mathrm dx\,\mathrm dy
\ \le\ C(T)\bigl(\sqrt{\lambda_i}+\sqrt{\mu_i}+\nu_i\bigr),
\qquad 0<t\le T .
$$

A companion mixed-derivative bound (eq. (3.18)) gives $\int_\Gamma\int_D\rho(y)(\partial^2u/\partial y_i\partial y_j)^2\mathrm dx\,\mathrm dy\lesssim C(T)(\lambda_i+\mu_i+\lambda_i\mu_i+\lambda_i\nu_i)$ under further assumptions.

> [!warning] An asymmetry in the exponents
> Equation (3.17) is transcribed as printed: the first two eigenvalue families appear under square roots and the third, $\nu_i$, does not. This is likely a typographical inconsistency in the paper, so the exponent on $\nu_i$ is **verified as printed but unverified in intent**. The conclusion drawn from the estimate — low-index directions matter more — does not depend on that exponent, but any quantitative use of the bound should be checked against the original.

The paper lists one further contribution: sufficient conditions on the data guaranteeing that the solution lies in appropriate random spaces ($BV$ and $H^k$), generalising paper 2 from a random *variable* wave speed to a random *field* wave speed.

### Numerical experiments

Example 4.1 has the following setup and measurements:

| Item        | Setting                                                                     |
| ----------- | --------------------------------------------------------------------------- |
| Equation    | $\partial_t u=\kappa(y)\partial_x u$, $x\in[-1,1]$                          |
| Wave speed  | $\kappa(y)=0.5\sum_{n=1}^{4}n^{-2}y_n$, each $y_n$ uniform on $(-1,1)$      |
| Initial data | $u(x,0,y)=\sin(x)$                                                         |
| Boundary    | $u(\pm1,t,y)=\sin(\pm1+\kappa(y)t)$, imposed by the sign of $\kappa$        |
| Errors      | $e_{mean}$, $e_{std}$, $e_2$ (defined below)                                |
| Table 1     | nonuniform (anisotropic) convergence results at $t=2$                       |

The three error measures are

$$
e_{mean}=\max_x\bigl|\mathbb E(u_{num})-\mathbb E(u_{exa})\bigr|,
\quad
e_{std}=\max_x\bigl|\sigma_{u_{num}}-\sigma_{u_{exa}}\bigr|,
\quad
e_2=\max_x\bigl(\mathbb E(u_{num}-u_{exa})^2\bigr)^{1/2}.
$$

Uniform distributions are used throughout, on the grounds that the sensitivity results (3.17)-(3.18) are distribution-independent.

The range of this wave speed can be computed directly: $0.5\sum_{n=1}^{4}n^{-2}\approx0.5\times1.4236\approx0.712$, so $\kappa$ sweeps $[-0.712,\,0.712]$ over the parameter domain and **does change sign**. The mode-by-mode sign tests above are therefore genuinely exercised, not a formality.

Treating the $\kappa_n$ as constants gives $\sqrt{\lambda_n}\propto n^{-2}$, so the direction weights in the sensitivity estimate (3.17) decay like $n^{-2}$, consistent with the anisotropic orders the paper uses:

| $i$ | Coefficient $0.5\,i^{-2}$ | Relative weight | Order $r_i$ | Basis count $r_i+1$ |
| --- | ------------------------- | --------------- | ----------- | ------------------- |
| 1   | $0.5$                     | $1$             | 3           | 4                   |
| 2   | $0.125$                   | $1/4$           | 2           | 3                   |
| 3   | $\approx0.0556$           | $1/9$           | 1           | 2                   |
| 4   | $0.03125$                 | $1/16$          | 1           | 2                   |

The findings are that numerical and exact solutions agree well even at the very low anisotropic orders $(r_1,r_2,r_3,r_4)=(3,2,1,1)$; that isotropic refinement converges fast; and that the error grows with $t$, the paper describing the dependence on time as typically linear and noting that this long-time-integration degradation is a known problem for generalised polynomial chaos and collocation alike.

By the paper's own formula $N_y=\prod_i(r_i+1)$, that order tuple gives $4\times3\times2\times2=48$ basis functions, whereas an isotropic tensor basis at the same top order $r=3$ would need $4^4=256$ — **about a fifth of the degrees of freedom**. (These two numbers are computed from the paper's formula and its stated order tuple; they are not figures the paper reports.) The specific error values in Table 1 were not available here and are therefore not reported.

### Relation to the others

Paper 5 is the constructive Galerkin counterpart to paper 2 and shares its model problem; paper 4 then supplies analytic regularity and spectral rates covering both. Paper 3 applies the identical bi-orthogonal decoupling device to elliptic interface problems. The "double orthogonality equals hidden collocation" observation in Remark 3.2 is a conceptual bridge to the sampling-based papers that follow, and the anisotropic degree-selection rule prefigures the anisotropic and sparse index sets used in the [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|least-squares]], [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|optimal-sampling]] and [[en/computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce|sparse-recovery]] papers.

## 3: the same device for elliptic interface problems (unverifiable)

> [!warning] This paper could not be verified
> Crossref, OpenAlex, MaRDI/zbMATH and Semantic Scholar all carry no abstract for this DOI (the Semantic Scholar record explicitly states the abstract field has been elided by the publisher), and ScienceDirect returns HTTP 403 to automated retrieval. This site therefore **reports no theorem, hypothesis, constant, convergence rate or numerical result from this paper**, and quotes no text as its abstract. What follows is limited to what zbMATH indexing keywords, the paper's reference list, and independent papers by the same group in the same period confirm.

The verifiable record: Tao Zhou, sole author, _J. Comput. Appl. Math._ 236(5) (1 October 2011), pp. 782-792, DOI `10.1016/j.cam.2011.05.033`. The zbMATH/MaRDI record (Q651914) lists the indexing keywords "immersed finite element", "random elliptic interface problems", "bi-orthogonal polynomials" and "Galerkin method". Its reference list includes Li-Lin-Lin-Rogers (_Numer. Methods Partial Differential Equations_ **20** (2004) 338-367) and He-Lin-Lin (same journal, **24** (2008) 1265-1300) on immersed finite element spaces and their approximation capability. Semantic Scholar additionally carries a machine-generated summary — **not the author's words and no substitute for the abstract** — which is not quoted here; it is recorded only because it corroborates, consistently with every other source checked, two structural claims: bi-orthogonal polynomials are the decoupling device, and the output is an **uncoupled** system of deterministic interface problems.

The background itself is standard. Elliptic problems with an interface are hard enough deterministically, since standard finite elements lose accuracy unless the mesh is fitted to the interface; adding a random input makes the naive stochastic Galerkin route doubly expensive, because projection produces a coupled system of deterministic interface problems and each would need an interface-fitted mesh. The two verifiable claims together say that the paper removes both couplings at once: bi-orthogonal polynomials decouple the random direction (the mechanism is set out under paper 5), and immersed finite element spaces let the mesh in the physical direction be independent of the interface geometry.

What is not reported needs stating plainly. Whether the paper gives an error estimate, in what form, with what constants and at what rate, is unknown here. Secondary records indicate numerical experiments accompany the method, but the test problems and observed outcomes could not be confirmed. The exact form of the sensitivity estimate used to pick anisotropic polynomial orders in this paper is likewise unverified — the verified Theorem 3.1 and eqs. (3.17)-(3.18) of paper 5 show what such a device looks like, but cannot stand in for this paper's own statement.

## 7: delay differential equations with random input (abstract level)

> [!note] Verification level
> This section is limited to the abstract, keywords and article front matter. Equation-level detail, hypotheses, constants and experimental data are **unverified** here and are therefore not reported.

### The idea

Delay differential equations with random coefficients stack two difficulties. The delay term makes the solution's smoothness in **time** break at successive multiples of the lag — the classical propagation-of-derivative-discontinuities phenomenon — while the random coefficients raise the usual question of whether the solution is smooth enough **in the parameter** for a spectral method in the random space to pay off. There is real tension between the two: a Legendre spectral collocation method in time is fast only if the solution is smooth on the interval, and the delay destroys that smoothness at specific instants. Establishing parametric regularity is therefore a prerequisite, not a formality, and the paper indeed does the regularity first and proposes the scheme afterwards.

### Construction

The random space is handled by **stochastic collocation**: solve the deterministic delay differential equation independently at each collocation node in parameter space, then interpolate. The MaRDI keyword list includes "sparse grid", so with several random inputs the nodes are Smolyak-type sparse grids. Each resulting deterministic delay differential equation is solved in time by a **Legendre spectral collocation method**. Regularity in the random space is established first, "provided that the given data satisfy some reasonable assumptions" — the same pattern as papers 2, 4 and 5, where data hypotheses are converted into parametric regularity and then into a rate. The precise assumptions and the form of the regularity estimate in this paper are unverified here.

### Main results (abstract level)

Under reasonable assumptions on the data the exact solution admits good regularity in the random space; the paper gives a convergence analysis of the proposed method; and the method attains **the familiar exponential order of convergence in both the random space and the time space**. No explicit constant, explicit exponent or sample-complexity relation is recoverable from the sources verified here, so no specific rate constant is given. How the time discretisation handles the derivative discontinuities at multiples of the lag likewise cannot be reported.

### Numerical experiments

The abstract states that numerical examples are given to illustrate the theoretical results. Which delay problems are solved, what lag values are used and what error magnitudes are observed are unverified and not reported.

### Relation to the others

Methodologically this is the last of Zhou's "regularity first, then spectral rate" papers, in the same mould as 2, 4 and 5 but for a delay ODE rather than a transport PDE, and the only one in the list where sparse grids are the sampling device. It sits at the hinge of the collection: papers up to 7 ask what accuracy a **given** structured node set achieves, and papers from 9 onward ask how to **design** the node set — see the [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|least-squares page]], where paper 6, published the same year and also using sparse grids as one candidate, is the transitional study finding structured sparse grids uncompetitive beyond very low dimension.

## 38: Maxwell's equations with random inputs (abstract level)

> [!note] Verification level
> Only the abstract, metadata and reference list were accessible; the full text is behind a paywall and no preprint was found. The construction and results below are reported as the abstract states them; hypotheses, constants, CFL conditions and orders of convergence are **unverified** and are not reported.

### The idea

Time-dependent Maxwell's equations form a first-order hyperbolic system and therefore inherit all the regularity difficulties of papers 2, 4 and 5. They also carry something worth preserving: an **energy conservation law**. Whether a Galerkin truncation inherits such a quadratic invariant is decided by the structure of the coupled system the projection produces — the same class of coefficient matrices paper 1 studies. Those are exactly the paper's two concerns: whether the projection is energy preserving, and how to advance the coupled Galerkin system efficiently.

The contrast with paper 5 is clean: both take a hyperbolic system, project it with generalised polynomial chaos, and then engineer the resulting coupled deterministic system to decouple, but **paper 5 achieves exact decoupling by changing basis (bi-orthogonal polynomials) while paper 38 achieves it by modifying the time-stepping**.

### Construction (abstract level)

Generalised polynomial chaos first converts the random Maxwell system into a deterministic **Galerkin system** for the expansion coefficients; a **finite element method** discretises the physical space; and three time-domain schemes are constructed: a **Crank-Nicolson** scheme, a **classical leap-frog** scheme, and a **modified leap-frog type** scheme designed so that the coupled Galerkin system can be advanced in a decoupled way.

### Main results (abstract level)

- The stochastic Galerkin approach **preserves the energy conservation law** at the generalised-polynomial-chaos-projected level.
- Error estimates are presented for the finite element approach used to solve the Galerkin system. **Order of convergence unverified.**
- For the **Crank-Nicolson** scheme, the **fully discrete** scheme is shown to be energy preserving.
- For the **classical leap-frog** scheme, a **conditional** energy stability property is shown. **The explicit CFL-type condition is unverified.**
- The **modified leap-frog type** scheme decouples the Galerkin system, the abstract calling it "a very efficient numerical approach". **Any quantitative speed-up is unverified.**

### Numerical experiments

The abstract states that numerical examples are presented to support the theoretical findings, and the article contains three figures. The specific test problems and quantitative outcomes are unverified and not reported.

### Relation to the others

This is the last of the intrusive-Galerkin strand and cites four papers from this collection directly: paper 3 (stochastic Galerkin for elliptic interface problems), paper 5 (bi-orthogonal-polynomial Galerkin for stochastic hyperbolic problems), paper 2 (convergence analysis for stochastic collocation for scalar hyperbolic equations) and paper 14 (stochastic collocation on unstructured meshes, see the [[en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|least-squares page]]). Its reference list also includes Xiu-Shen's _Efficient stochastic Galerkin methods for random diffusion equations_, the source of the coefficient matrices paper 1 analyses.

> [!note] Coverage status
> Papers **1, 2, 4 and 5** have been checked against their full texts, and this page gives their intuition, derivation chains, theorems with hypotheses and experimental setups (paper 1 has no numerical experiments at all — it is analysis throughout). Papers **7 and 38** reach only the level of abstract and metadata: their constructions and results are reported as the abstracts state them, and their hypotheses, constants, CFL conditions, convergence orders and experimental data are not reported. Paper **3** is the one item that could not be verified at all: Crossref, OpenAlex, MaRDI/zbMATH and Semantic Scholar all lack its abstract (Semantic Scholar records explicitly that the abstract field was elided by the publisher) and ScienceDirect returns 403 to automated retrieval, so this page gives none of its theorems, constants, rates or numerical results, and reports only the constructional idea confirmable from zbMATH indexing keywords, the reference list, and contemporaneous papers by the same group.
> Several formulas transcribed as printed but problematic are flagged in place: the denominator of the normalised recurrence coefficient $\tilde a_n$ in paper 1, the analytic-extension series and the closed form of $r$ in paper 4, the mutual inconsistency of the three data items in paper 4's numerical example, and the exponent on $\nu_i$ in paper 5's eq. (3.17).

## Coverage check

| Item                                                | Paper | Status                                                     |
| --------------------------------------------------- | ----- | ---------------------------------------------------------- |
| Separating characteristics and accuracy lost in time | 2/4/5 | explicit solution, $t^k$ growth of $\partial_y^k u$, consequences |
| Galerkin coefficient matrix and its sparsity        | 1     | $a_{jk}$, $e_{ijk}$, symmetry, derivation of $2N+1$ nonzeros |
| Reduction to three-term recurrence coefficients     | 1     | both identities, their derivation and technical role        |
| Jacobi coefficients and the $b_n=0$ criterion       | 1     | $h_n$, $a_n$, $b_n$, $c_n$, unverified denominator flagged  |
| Theorem 1 and the threshold $\lvert\alpha\rvert\ge1/2$ | 1  | both densities, conclusion, qualification, convention       |
| Sign-changing speed and switching inflow            | 2     | problem setting and the regularity consequence              |
| Data hypotheses weighted by $1/\lvert c\rvert$      | 2     | all three conditions, energy derivation, weight reconciled  |
| Three regularity and three convergence theorems     | 2     | $H^1$/$H^2$/$BV$ and $N^{-1}$/$N^{-2}$/mean only            |
| Numerical evidence that the hypotheses are sharp    | 2     | orders for three initial data, half order in the $BV$ case  |
| Geometric growth and analytic extension             | 4     | hypothesis form, Leibniz induction, origin of $1/\sqrt\zeta$ |
| Exponential convergence for collocation             | 4     | the estimate, relation between $\varrho$ and $\tau$         |
| Hyperbolic solutions not analytic in general        | 4     | both qualifications and the "rather rough" self-assessment  |
| The example and a quantitative reading of time decay | 4    | setup table, inconsistency flagged, $\tau\sim1/t$           |
| Both defining conditions of bi-orthogonality        | 5     | definition, eigenproblem, same matrix as paper 1            |
| Exact decoupling and sign tests on effective speeds | 5     | derivation of (3.8)-(3.10), $C_{k,i}$ inside the support    |
| Sensitivity basis for anisotropic orders            | 5     | (3.17), the $\nu_i$ exponent flagged, (3.18)                |
| Setup of Example 4.1 and the saving in unknowns     | 5     | sign-changing speed range, weight table, $48$ against $256$ |
| Double decoupling for interface problems            | 3     | only the idea confirmable from keywords and references      |
| Collocation and spectral time stepping for delay ODEs | 7   | abstract-level construction and "exponential in both spaces" |
| Energy preservation and decoupled stepping for Maxwell | 38 | abstract-level schemes and their structural properties      |

## Sources for this page

- T. Zhou and T. Tang, [_Note on coefficient matrices from stochastic Galerkin methods for random diffusion equations_](https://doi.org/10.1016/j.jcp.2010.07.016), J. Comput. Phys. 229 (2010), pp. 8225-8230.
- T. Tang and T. Zhou, [_Convergence analysis for stochastic collocation methods to scalar hyperbolic equations with a random wave speed_](https://doi.org/10.4208/cicp.060109.130110a), Commun. Comput. Phys. 8 (2010), pp. 226-248.
- T. Zhou, [_Stochastic Galerkin methods for elliptic interface problems with random input_](https://doi.org/10.1016/j.cam.2011.05.033), J. Comput. Appl. Math. 236 (2011), pp. 782-792.
- T. Zhou and T. Tang, [_Convergence analysis for spectral approximation to a scalar transport equation with a random wave speed_](https://doi.org/10.4208/jcm.1206-m4012), J. Comput. Math. 30 (2012), pp. 643-656.
- T. Zhou and T. Tang, [_Galerkin methods for stochastic hyperbolic problems using bi-orthogonal polynomials_](https://doi.org/10.1007/s10915-011-9508-0), J. Sci. Comput. 51 (2012), pp. 274-292.
- T. Zhou, [_A stochastic collocation method for delay differential equations with random input_](https://doi.org/10.4208/aamm.2012.m38), Adv. Appl. Math. Mech. 6 (2014), pp. 403-418.
- Z. Feng, J. Li, T. Tang, and T. Zhou, [_Efficient stochastic Galerkin methods for Maxwell's equations with random inputs_](https://doi.org/10.1007/s10915-019-00936-z), J. Sci. Comput. 80 (2019), pp. 248-267.
