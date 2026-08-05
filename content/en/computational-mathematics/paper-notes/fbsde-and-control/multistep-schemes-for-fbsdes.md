---
title: Multistep Schemes
description: Papers 8, 18, 23, 33, 35, 61 and 68 - raising the temporal order of the backward equation with several future levels
lang: en
translation: computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes
tags:
  - paper-notes
  - stochastic-differential-equations
  - multistep-schemes
---

> [!note] Coverage of this page
> Papers **8** (_SIAM J. Sci. Comput._ 36(4), 2014), **18** (_J. Sci. Comput._ 69(2), 2016), **23** (_Numer. Math. Theor. Meth. Appl._ 10(2), 2017), **33** (_SIAM J. Numer. Anal._ 56(4), 2018), **35** (_J. Sci. Comput._ 79(3), 2019), **61** (_J. Comput. Math._ 40(4), 2022) and **68** (_J. Sci. Comput._ 94:53, 2023).
>
> Of these, **8** (from the full arXiv preprint) and **23** (from the published text) were checked equation by equation, so this page gives their complete derivations, theorems and experiment tables. The bodies of **18, 33, 35, 61 and 68** could not be obtained — none has a preprint, Springer and SIAM require subscriptions, and the Global Science Press PDF endpoint returns 403 or is unreachable — so for those five only what the abstracts and reference lists support is written here. Their schemes, theorems, convergence orders and experimental data are not reported on this site.

![How multistep schemes raise the backward accuracy](assets/diagrams/tao-zhou-papers/en/fbsde-multistep.svg)

## 8: moving the entire high-order requirement onto the backward direction

### The idea

Solving a forward-backward system numerically, one's instinct is that both directions must be accurate: if the forward $X$ is wrong, the backward $Y,Z$ cannot be right. But the backward equation only ever sees the forward process through **conditional expectations** — what appears in a scheme is never a particular trajectory of $X$ but a quantity like $\mathbb E^x_{t_n}[\,\cdot\,]$. A conditional expectation needs only the **transition law**, not the paths. And on the scale of a single time step, the first-order change of the transition law of a diffusion started at $(t_n,x)$ is determined by **the coefficients at that point alone**: that is what a generator is.

This observation is the entire mechanism of paper 8. It means the true forward diffusion may be replaced by one whose coefficients are frozen at $(t_n,x)$ — an Euler step — without losing any order in the backward direction. The paper is explicit about what happens without the observation: to get order $k$ in $Y,Z$ one would have to build an order-$k$ strongly convergent scheme for the forward SDE as well, and such schemes (Itô-Taylor type) require repeated multiple stochastic integrals, which are expensive and nearly unimplementable in the coupled case, where the forward coefficients depend on the not-yet-computed $(Y,Z)$.

What is left is a purely backward question: how to raise the temporal order using information from several future time levels. The paper's answer is not to interpolate the integrand but to **differentiate the integral identities in $t$ into two ordinary differential equations, then approximate the derivative with a forward-looking BDF-type stencil**. That choice pays an unexpected dividend: at $t=t_n$ every term of the form $f\cdot\Delta W_{t_n,t}$ vanishes because $\Delta W_{t_n,t_n}=0$, so $Z^n$ comes out **explicit** and only $Y^n$ needs iteration.

### Setting

The coupled Markovian forward-backward system is

$$
X_t=X_0+\int_0^t b(s,X_s,Y_s,Z_s)\,\mathrm ds+\int_0^t\sigma(s,X_s,Y_s,Z_s)\,\mathrm dW_s,
$$

$$
Y_t=\xi+\int_t^T f(s,X_s,Y_s,Z_s)\,\mathrm ds-\int_t^T Z_s\,\mathrm dW_s,
\qquad \xi=\varphi(X_T),
$$

with $b:\Omega\times[0,T]\times\mathbb R^q\times\mathbb R^p\times\mathbb R^{p\times d}\to\mathbb R^q$, $\sigma:\cdots\to\mathbb R^{q\times d}$, generator $f:\cdots\to\mathbb R^p$, and the unknown triple $(X_t,Y_t,Z_t)$ valued in $\mathbb R^q\times\mathbb R^p\times\mathbb R^{p\times d}$. The system is decoupled when $b,\sigma$ do not depend on $(Y_t,Z_t)$.

The paper describes the state of the art as follows: most existing schemes are Euler-type with rate $1/2$; the high-order methods for decoupled systems rely "on the high order approaches for both the forward SDE and the high order scheme for the backward SDE", and high-order approaches for the forward equation "require large amounts of computations and are often difficult to be applied". In the coupled case the forward coefficients depend on the backward unknowns, so "it seems not easy to design high-order (yet efficient) numerical schemes". The paper poses and answers affirmatively the explicit question: **can one still expect high-order accurate numerical solutions of the backward equation if the Euler method is used to solve the forward equation?**

### Derivation

**Step one: the generator only looks at the left endpoint.** For $X_s$ solving $\mathrm dX_s=b\,\mathrm ds+\sigma\,\mathrm dW_s$ the generator is

$$
\mathcal A_t^x g(t,x)=\lim_{s\downarrow t}\frac{\mathbb E_t^x[g(s,X_s)]-g(t,x)}{s-t},
\qquad
\mathcal L^0_{t,x}=\frac{\partial}{\partial t}+\sum_i b_i\frac{\partial}{\partial x_i}
+\frac12\sum_{i,j}(\sigma\sigma^{\top})_{i,j}\frac{\partial^2}{\partial x_i\partial x_j},
$$

with $\mathcal A_t^x f=\mathcal L^0_{t,x}f$ on $C^{1,2}$ functions. Theorem 2 of the paper states that if $f\in C^{1,2}$ and $\mathbb E^{x_0}_{t_0}[\mathcal L^0_{t,X_t}f(t,X_t)]<\infty$, then for $t\ge t_0$

$$
\frac{\mathrm d\,\mathbb E^{x_0}_{t_0}[f(t,X_t)]}{\mathrm dt}
=\mathbb E^{x_0}_{t_0}\bigl[\mathcal A_t^{X_t}f(t,X_t)\bigr],
$$

and that **at $t=t_0$** this derivative is unchanged if $X_t$ is replaced by **any** diffusion $\bar X_t=x+\int_{t_0}^t\bar b_s\,\mathrm ds+\int_{t_0}^t\bar\sigma_s\,\mathrm dW_s$ whose coefficients merely match at the left endpoint:

$$
\bar b(t_0,\bar X_{t_0};t_0,x_0)=b(t_0,x_0),
\qquad
\bar\sigma(t_0,\bar X_{t_0};t_0,x_0)=\sigma(t_0,x_0)
\quad\Longrightarrow\quad
\left.\frac{\mathrm d\,\mathbb E[f(t,X_t)]}{\mathrm dt}\right|_{t_0}
=\left.\frac{\mathrm d\,\mathbb E[f(t,\bar X_t)]}{\mathrm dt}\right|_{t_0}.
$$

The paper notes explicitly that one may simply take $\bar b(s,\bar X_s;t_0,x_0)=b(t_0,x_0)$ and $\bar\sigma(s,\bar X_s;t_0,x_0)=\sigma(t_0,x_0)$ for all $s\in[t_0,t]$ — that is the Euler step. **The general lesson deserves separate emphasis:** when the two directions of a coupled system carry different accuracy requirements, first check whether the error from the cheaper direction actually reaches the quantity of interest. Here it does not, because what reaches the quantity of interest is only a one-step derivative approximation, and that sees only the left endpoint.

**Step two: the weights for the derivative approximation.** For $u\in C_b^{k+1}$ and nodes $t_0<t_1<\cdots<t_k$ with $\Delta t_i=t_i-t_0$, the weights are fixed by the moment conditions

$$
\sum_{i=0}^{k}\alpha_{k,i}\frac{(\Delta t_i)^j}{j!}=\delta_{j1}
=\begin{cases}1,&j=1\\0,&j\ne1\end{cases},
\qquad j=0,1,\dots,k
$$

giving $\frac{\mathrm du}{\mathrm dt}(t_0)=\sum_{i=0}^k\alpha_{k,i}u(t_i)+\mathrm{Err}$ with $\mathrm{Err}=O\bigl(\sum_i\alpha_{k,i}(\Delta t_i)^{k+1}\bigr)$. On a uniform grid $\Delta t_i=i\Delta t$ this is a Vandermonde system

$$
\begin{pmatrix}
1&1&\cdots&1\\ 0&1&\cdots&k\\ 0&1^2&\cdots&k^2\\ \vdots&&&\vdots\\ 0&1^k&\cdots&k^k
\end{pmatrix}
\begin{pmatrix}\alpha_{k,0}\Delta t\\ \alpha_{k,1}\Delta t\\ \vdots\\ \alpha_{k,k}\Delta t\end{pmatrix}
=\begin{pmatrix}0\\1\\0\\\vdots\\0\end{pmatrix}.
$$

The resulting $\alpha_{k,i}\Delta t$ are exactly the **BDF coefficients read in reverse** (a backward differentiation formula pointing forward in time):

| $k$ | $i=0$     | $1$ | $2$     | $3$    | $4$     | $5$   | $6$    |
| --- | --------- | --- | ------- | ------ | ------- | ----- | ------ |
| 1   | $-1$      | $1$ |         |        |         |       |        |
| 2   | $-3/2$    | $2$ | $-1/2$  |        |         |       |        |
| 3   | $-11/6$   | $3$ | $-3/2$  | $1/3$  |         |       |        |
| 4   | $-25/12$  | $4$ | $-3$    | $4/3$  | $-1/4$  |       |        |
| 5   | $-137/60$ | $5$ | $-5$    | $10/3$ | $-5/4$  | $1/5$ |        |
| 6   | $-49/20$  | $6$ | $-15/2$ | $20/3$ | $-15/4$ | $6/5$ | $-1/6$ |

**Step three: the two reference ordinary differential equations.** This is the structural core of the paper. Taking $\mathbb E^x_{t_n}[\cdot]$ of the backward equation gives the integral identity $\mathbb E^x_{t_n}[Y_t]=\mathbb E^x_{t_n}[\xi]+\int_t^T\mathbb E^x_{t_n}[f(s,X_s,Y_s,Z_s)]\,\mathrm ds$, and differentiating in $t$:

$$
\frac{\mathrm d\,\mathbb E^x_{t_n}[Y_t]}{\mathrm dt}
=-\,\mathbb E^x_{t_n}\bigl[f(t,X_t,Y_t,Z_t)\bigr],\qquad t\in[t_n,T].
$$

Multiplying $Y_{t_n}=Y_t+\int_{t_n}^tf\,\mathrm ds-\int_{t_n}^tZ_s\,\mathrm dW_s$ by $(\Delta W_{t_n,t})^{\top}$ and taking $\mathbb E^x_{t_n}[\cdot]$ — the Itô isometry makes $\mathbb E^x_{t_n}[Y_{t_n}(\Delta W_{t_n,t})^\top]=0$ and turns the stochastic integral term into $\int_{t_n}^t\mathbb E^x_{t_n}[Z_s]\mathrm ds$ — gives

$$
0=\mathbb E^x_{t_n}\bigl[Y_t(\Delta W_{t_n,t})^{\top}\bigr]
+\int_{t_n}^t\mathbb E^x_{t_n}\bigl[f(s,X_s,Y_s,Z_s)(\Delta W_{t_n,s})^{\top}\bigr]\mathrm ds
-\int_{t_n}^t\mathbb E^x_{t_n}[Z_s]\,\mathrm ds,
$$

and differentiating again in $t$ gives the second reference equation:

$$
\frac{\mathrm d\,\mathbb E^x_{t_n}\bigl[Y_t(\Delta W_{t_n,t})^{\top}\bigr]}{\mathrm dt}
=-\,\mathbb E^x_{t_n}\bigl[f(t,X_t,Y_t,Z_t)(\Delta W_{t_n,t})^{\top}\bigr]+\mathbb E^x_{t_n}[Z_t].
$$

**Evaluating at $t=t_n$ is the decisive step**: $\Delta W_{t_n,t_n}=0$ annihilates the terms containing $f\cdot\Delta W$, leaving only $-f(t_n,x,Y_{t_n},Z_{t_n})$ and $Z_{t_n}$ on the right-hand sides. That is why the schemes below carry **no $f$ term at the future levels**, and why the $Z$ equation is fully explicit.

**Step four: the semi-discrete scheme.** Applying the derivative weights to the left-hand sides of the two reference equations and invoking Theorem 2 to swap in the frozen diffusion:

$$
\sum_{i=0}^k\alpha_{k,i}\,\mathbb E^x_{t_n}\bigl[\bar Y_{t_{n+i}}\bigr]=-f(t_n,x,Y_{t_n},Z_{t_n})+R^k_{y,n},
\qquad
\sum_{i=1}^k\alpha_{k,i}\,\mathbb E^x_{t_n}\bigl[\bar Y_{t_{n+i}}(\Delta W_{n,i})^{\top}\bigr]=Z_{t_n}+R^k_{z,n}.
$$

Dropping the truncation terms and taking the Euler choice $\bar b(s,\cdot)=b(t_n,x)$, $\bar\sigma(s,\cdot)=\sigma(t_n,x)$ yields the paper's Scheme 2:

$$
X^{n+j}=X^n+b(t_n,X^n)\Delta t_{n,j}+\sigma(t_n,X^n)\Delta W_{n,j},\qquad j=1,\dots,k,
$$

$$
Z^n=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}(\Delta W_{n,j})^{\top}\bigr],
\qquad
\alpha_{k,0}Y^n=-\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}\bigr]-f(t_n,X^n,Y^n,Z^n).
$$

One detail is easy to miss: $X^{n+j}$ is produced for **all** $j=1,\dots,k$ by a **single** Euler step from $t_n$, with step $\Delta t_{n,j}=t_{n+j}-t_n$ and increment $\Delta W_{n,j}=W_{t_{n+j}}-W_{t_n}$, not by $j$ successive Euler steps. That is precisely what Theorem 2 permits, and it is why the scheme is cheap.

**Step five: the fully discrete scheme.** Introduce a spatial grid $\mathcal D^n_h\subset\mathbb R^q$ with density $h_n=\max_x\mathrm{dist}(x,\mathcal D^n_h)$, local neighbour sets $\mathcal D^n_{h,x}$ with $\#\mathcal D^n_{h,x}\le N_e$, a local interpolation operator $\mathcal I^n_{h,x}$, and a quadrature approximation $\mathbb E^{x,h}_{t_n}[\cdot]$ of the conditional expectation. Given $Y^{N-i},Z^{N-i}$ for $i=0,\dots,k-1$, for $n=N-k,\dots,0$ and each $x\in\mathcal D^n_h$:

$$
Z^n=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{x,h}_{t_n}\bigl[\mathcal I^{n+j}_{h,X^{n+j}}Y^{n+j}(\Delta W_{n,j})^{\top}\bigr],
$$

$$
\alpha_{k,0}Y^{n}=-\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{x,h}_{t_n}\bigl[\mathcal I^{n+j}_{h,X^{n+j}}Y^{n+j}\bigr]-f(t_n,x,Y^n,Z^n).
$$

The order of operations at each level is: Euler-step to get $X^{n+j}$, compute $Z^n$ explicitly, then solve for $Y^n$ implicitly by the Picard iteration $\alpha_{k,0}Y^{n,l+1}=-\sum_j\alpha_{k,j}\mathbb E^{x,h}_{t_n}[\mathcal I Y^{n+j}]-f(t_n,x,Y^{n,l},Z^n)$ until $|Y^{n,l+1}-Y^{n,l}|\le\epsilon_0$. **The cost should be recorded here too: a $k$-step scheme needs $k$ starting levels and visits $k$ future levels at every step.** That is the structural dividing line between it and deferred correction (paper 23).

**Step six: the coupled case.** Scheme 4 replaces the Euler step by $X^{n+j}=X^n+b(t_n,X^n,Y^n,Z^n)\Delta t_{n,j}+\sigma(t_n,X^n,Y^n,Z^n)\Delta W_{n,j}$, so that $X^{n+j}$ and $(Y^n,Z^n)$ now depend on each other. Scheme 5 unties this with an outer Picard iteration: set $Y^{n,0}=Y^{n+1}$, $Z^{n,0}=Z^{n+1}$, and for $l=0,1,\dots$ take the Euler step using $(Y^{n,l},Z^{n,l})$, then compute $Z^{n,l+1}$ explicitly and $Y^{n,l+1}$ implicitly, until $\max(|Y^{n,l+1}-Y^{n,l}|,|Z^{n,l+1}-Z^{n,l}|)<\epsilon_0$. If $b,\sigma$ do not depend on $(Y,Z)$, Scheme 5 degenerates to Scheme 3.

**Step seven: how the conditional expectations are computed.** The Euler step makes $\Delta W_{n,j}\sim\sqrt{\Delta t_{n,j}}\,N(0,I_d)$, so every conditional expectation is an integral against a known Gaussian density. The paper uses Gauss-Hermite quadrature

$$
\int_{-\infty}^{+\infty}e^{-x^2}g(x)\,\mathrm dx\approx\sum_{j=1}^{L}\omega_jg(a_j),
\qquad
\omega_j=\frac{2^{L+1}L!\sqrt\pi}{\bigl(H_L'(a_j)\bigr)^2},
$$

with $\{a_j\}$ the roots of the degree-$L$ Hermite polynomial $H_L$ and truncation error $R(g,L)=\frac{L!\sqrt\pi}{2^L(2L)!}g^{(2L)}(\eta)$, hence exactness for polynomials of degree $\le2L-1$. In $d$ dimensions a tensor rule $\omega_{\mathbf j}=\prod_i\omega_{j_i}$ is used, converted through $\mathbb E[g(N)]=\pi^{-d/2}\int g(\sqrt2x)e^{-x^\top x}\mathrm dx$. Spatially, $\mathcal I^n_{h,x}$ is local Lagrange interpolation. **The quadrature nodes are roots of a Hermite polynomial and in general do not sit on the spatial grid, so the interpolation step cannot be avoided** — which is exactly what paper 63 later removes by switching to Sinc quadrature.

### Theorems

- **Generator theorem (Theorem 2).** Hypotheses: $f\in C^{1,2}$ and $\mathbb E^{x_0}_{t_0}[\mathcal L^0_{t,X_t}f(t,X_t)]<\infty$. Conclusion as above: the time derivative of the conditional expectation equals the conditional expectation of the generator, and at the left endpoint it depends only on the left-endpoint coefficients. This is the licence for the whole construction.
- **Local truncation error.** At the semi-discrete level $\bar R^k_{y,n}=\bar R^k_{z,n}=O((\Delta t)^k)$. At the fully discrete level the truncation error splits into six pieces: $R^k_{y,n},R^k_{z,n}$ (derivative approximation), $R^{k,I_h}_{y,n},R^{k,I_h}_{z,n}$ (interpolation) and $R^{k,E}_{y,n},R^{k,E}_{z,n}$ (quadrature for the conditional expectations). Under smooth data and interpolating polynomials of degree $r$,
  $$
  R^k_{y,n}=O\bigl((\Delta t_n)^k\bigr),\quad R^k_{z,n}=O\bigl((\Delta t_n)^k\bigr),\quad
  R^{k,I_h}_{y,n}=O\bigl(h^{r+1}\bigr),\quad R^{k,I_h}_{z,n}=O\bigl(h^{r+1}\bigr).
  $$
- **Error balancing.** To keep the temporal and spatial contributions of the same size the paper takes $h=(\Delta t)^{(k+1)/(r+1)}$.
- **Zero stability sets the order ceiling.** The window $1\le k\le6$ is not a vague empirical observation but a check anyone can redo. Applying the same weights to the deterministic ODE $\mathrm dY/\mathrm dt=f(t,Y)$ gives $\alpha_{k,0}Y^n+\sum_{j=1}^k\alpha_{k,j}Y^{n+j}=f(t_n,Y^n)$, with characteristic polynomial
  $$
  P(\lambda)=\alpha_{k,0}\lambda^{k}+\sum_{j=1}^{k}\lambda^{k-j}=0,
  $$
  subject to the root condition $|\lambda_{k,j}|\le1$, simple where equality holds. The maximum root moduli the paper reports, excluding the common root $1.0$, are

  | $k$              | 2      | 3      | 4      | 5      | 6      | 7          | 8          |
  | ---------------- | ------ | ------ | ------ | ------ | ------ | ---------- | ---------- |
  | max root modulus | 0.3333 | 0.4264 | 0.5608 | 0.7087 | 0.8633 | **1.0222** | **1.1839** |

  **The modulus crosses $1$ between $k=6$ and $k=7$, so the scheme is unstable for $k\ge7$**, which is why the tabulation stops at $k=6$.

> [!warning] This paper has no convergence theorem
> The arXiv version gives **no** rigorous convergence theorem with explicit constants. The high-order claim rests on the truncation-error estimates above plus numerical evidence; the paper's own words about Scheme 4 are only that "we can expect that Scheme 4 is a high order numerical scheme … we shall numerically verify this". Whether the published SIAM version adds a theorem is unverified here. Rigorous stability and convergence analysis for this family arrives with paper 47, on the [[en/computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes|stability theory page]].

### Numerical experiments

Common setup: $T=1.0$, uniform partitions, 8 Gauss-Hermite nodes (so quadrature error is negligible), FORTRAN 95 implementation. The tables below report observed convergence rates, not absolute errors.

**Example 1 (decoupled, analytic solution).** Exact solution $Y_t=\dfrac{e^{t+X_t}}{1+e^{t+X_t}}$, $Z_t=\dfrac{(e^{t+X_t})^2}{(1+e^{t+X_t})^3}$, with $x=1.0$ and $N=16,\dots,256$:

| $k$        | 1     | 2     | 3     | 4     | 5     | 6             | 7     | 8        |
| ---------- | ----- | ----- | ----- | ----- | ----- | ------------- | ----- | -------- |
| $Y$ rate   | 1.000 | 1.973 | 3.002 | 3.922 | 5.196 | 5.116 (6.273) | 4.382 | $-5.487$ |
| $Z$ rate   | 1.000 | 2.021 | 2.893 | 3.919 | 5.017 | 5.687 (6.256) | 4.759 | $-7.170$ |

The bracketed figures are the rates over $N=16,\dots,128$, that is, before double-precision round-off takes over — **the nominal order at $k=6$ is only visible on the coarser meshes; on the finer ones round-off has already eaten it**. At $k=7$ the rate degrades to between $4$ and $5$ rather than approaching $7$, and at $k=8$ it goes negative (divergence). This matches the root-condition computation above.

**Example 2 (European call under Black-Scholes).** Parameters $b=0.05$, $\sigma=0.2$, $r=0.03$, $d=0.04$, $T=1.0$, $K=S_0=100$:

| $k$      | 1     | 2     | 3     | 4     |
| -------- | ----- | ----- | ----- | ----- |
| $Y$ rate | 1.002 | 1.964 | 2.935 | 3.957 |
| $Z$ rate | 1.002 | 1.998 | 2.819 | 3.991 |

**Example 3 (coupled, two systems).** System (5.10) has exact solution $Y_t=\sin(t+X_t)$, $Z_t=\sqrt2\cos(t+X_t)\sin^2(t+X_t)$, and its diffusion coefficient $\sigma=\sqrt2\,Y_s\sin(s+X_s)$ **deliberately violates** the paper's uniform ellipticity assumption (2.11); system (5.11) has $Z_t=\sqrt2\cos(t+X_t)(\sin^2(t+X_t)+1)$ and satisfies it.

| $k$             | 1     | 2     | 3     | 4     |
| --------------- | ----- | ----- | ----- | ----- |
| (5.10) $Y$ rate | 0.984 | 1.984 | 2.935 | 4.006 |
| (5.10) $Z$ rate | 0.977 | 2.064 | 2.976 | 4.055 |
| (5.11) $Y$ rate | 1.081 | 2.273 | 2.916 | 3.822 |
| (5.11) $Z$ rate | 1.023 | 1.948 | 2.970 | 4.213 |

Example 4 extends $\sigma$ to depend on $(X_s,Y_s,Z_s)$; this site confirms the example exists but did not transcribe its data table.

**What these experiments establish, and where they fall short.** They establish that on smooth one-dimensional problems, including a coupled one that **violates the ellipticity assumption**, the $k$-step scheme really does exhibit order $k$, and that instability sets in exactly where the deterministic root condition predicts. Three shortfalls. First, every example has low spatial dimension and the conditional expectations are computed with tensor Gauss-Hermite quadrature, so these tables say nothing about what happens as dimension grows — that has to wait for the sparse grids of paper 25. Second, the $k=6$ rate is polluted by round-off, which shows that **in double precision the practical ceiling of this family sits below the theoretical window's ceiling**. Third, observing order $k$ in Example 3 despite the violated ellipticity is a numerical observation, not a theorem; the paper does not, and cannot, relax the hypothesis on that basis.

### Relation to the others

This is the **root paper** of the whole thread: papers 16, 18, 19, 23, 25, 41, 61 and 68 all describe themselves as extensions of it. Papers 16 and 19 push into the second-order (fully nonlinear) setting, 18 adds jumps, 25 swaps the spatial discretisation for sparse grids, 23 and 35 swap the order mechanism for deferred correction, 47 supplies the rigorous stability and convergence framework, 63 swaps in Sinc interpolation and quadrature, 68 designs schemes backwards from stability, and 61 and 33 carry the construction into the mean-field setting.

### Distinguishing it from the other multistep route

Two constructions both called "multistep" need separating here, or their conclusions get attached to the wrong scheme.

- **Interpolate then integrate** (Zhao, Zhang and Ju, _SIAM J. Numer. Anal._ 48(4) 2010): Lagrange-interpolate the integrand of the reference integral identities over several future time levels, then integrate, producing Newton-Cotes-type weights of the form $h\sum_jb_j\mathbb E_i[f(t_{i+j},Y_{i+j},Z_{i+j})]$. Its stability windows **differ between the two directions**: the reference equation for $Y$ is stable only for $K_y\in\{1,\dots,7,9\}$ — note that $K_y=8$ is excluded — while the one for $Z$ is stable only for $K_z\in\{1,2,3\}$.
- **Differentiate into reference ODEs** (paper 8): it does **not** interpolate the integrand and integrate. It differentiates the reference integral identities in $t$, turning them into two reference ordinary differential equations, then discretises with the derivative-approximation weights above. Hence a single window $1\le k\le6$.

In the language of ODE solvers: Zhao-Zhang-Ju is Adams-type (interpolate the integrand and integrate), paper 8 is BDF-type (difference the derivative). Both are covered by the unified framework of paper 47.

One common misattribution is worth correcting in passing: the multistep, interpolate-then-integrate construction is due to Zhao-Zhang-Ju (2010), whereas Zhao, Chen and Peng (_SIAM J. Sci. Comput._ 28(4) 2006) is the origin of the **$\theta$-scheme**, which the multistep scheme extends. The two should not be conflated.

## 18: forward-backward systems with jumps

### The idea

Carrying the mechanism of paper 8 to jump diffusions runs into a combinatorial obstacle: a $k$-step scheme spans the interval $[t_n,t_{n+k}]$, and over that span the number of jumps has no a priori bound, so enumerating jump combinations one by one makes the cost explode. The paper defuses this with two observations: first, the generator of a jump diffusion is equally **local** — now an integro-differential operator $\mathcal L^0+\int_E[\cdot]\,\nu(\mathrm de)$ — so the forward direction can still use Euler alone; second, each time step accounts for **one jump only**, treating multi-jump events as higher-order contributions in $\Delta t$ and discarding them.

### Setting

The general form of a decoupled forward-backward system with jumps is

$$
X_t=X_0+\int_0^tb(s,X_{s^-})\,\mathrm ds+\int_0^t\sigma(s,X_{s^-})\,\mathrm dW_s
+\int_0^t\!\!\int_E\gamma(s,X_{s^-},e)\,\tilde\mu(\mathrm ds,\mathrm de),
$$

$$
Y_t=\varphi(X_T)+\int_t^Tf(s,X_s,Y_s,Z_s,\Gamma_s)\,\mathrm ds
-\int_t^TZ_s\,\mathrm dW_s-\int_t^T\!\!\int_EU_s(e)\,\tilde\mu(\mathrm ds,\mathrm de),
$$

with $\tilde\mu$ the compensated jump measure and $U$ the extra unknown attached to the jump part. This is the standard Barles-Buckdahn-Pardoux setting, and the generalised Feynman-Kac formula links it to a **partial integro-differential equation** whose nonlocal term comes from the jump measure.

### What could be verified

The abstract confirms three points verbatim: (i) "inspired by the local property of the jump diffusion processes, the Euler method is used to solve the corresponding forward stochastic differential equation with jumps, and this admits a dramatic reduction of the entire computational complexity, however, the high order rate of convergence is still maintained for the quantities of interest in the backward equation with jumps"; (ii) in each time step, the computational procedure involves **only one jump**, again dramatically reducing complexity; (iii) via the generalised Feynman-Kac formula the method applies readily to partial integro-differential equations (and certain nonlocal PDE models).

> [!note] What could be verified
> The body of the paper is paywalled with no preprint, so this site read only the abstract and reference list. The system above is written in the standard Barles-Buckdahn-Pardoux setting; **the notation the paper actually uses is unverified**. The reference equation matching the jump component $U$ in the multistep template, the formula-level implementation of "only one jump", and the specific order, hypotheses and constants behind the claimed high order are all unreported here. The abstract says "high order rate of convergence" without stating the order. For numerics, the abstract mentions "several numerical experiments", but the test problems and observed orders are unverified here.

### Relation to the others

A direct sibling of paper 8: the same "Euler for the forward direction" trick, the same $\alpha_{k,j}$ multistep template, extended to jumps and nonlocal equations. It is one of three orthogonal directions off the root paper, alongside papers 16 and 19 (towards second-order FBSDEs) and paper 25 (towards higher dimension). The first author, Yu Fu, is also an author of papers 8, 25 and 41.

## 23: trading one expensive scheme for a cheap one repeated

### The idea

The multistep way to gain order is to **build a wider stencil**: order $k$ costs $k$ future levels, $k$ starting values, and a root condition that only holds for $k\le6$. Deferred correction takes the other route: **use only the cheapest Euler scheme, but use it to solve a sequence of residual equations**. The key observation is that if the interpolant $\mathcal I u$ of an already-computed low-order solution is substituted back into the original equation, the error $\delta=y-\mathcal Iu$ satisfies an equation of **the same shape** as the original, so the same Euler scheme can solve it; adding the result back raises the order by one. Repeat $J$ times and the order rises by $J$, up to what the subgrid can support.

The appeal of this route is structural: no $\alpha_{k,i}$ stencil, no starting-value problem, and therefore no $k\le6$ stability barrier. The paper's own diagnosis is that "due to the involvement of randomness and the coupling of the forward and backward equations, it is difficult to design high order and relatively 'clean' numerical schemes for FBSDEs." The selling point of deferred correction is precisely "simplicity and robustness".

### Setting and the derivation in the ODE case

Start from the ODE $y'(t)=f(t,y(t))$, $y(0)=y_0$. Take a partition $0=t_0<\cdots<t_N=T$ and refine each $I_n=[t_n,t_{n+1}]$ into $\mathcal G^n_K=\{t_{n,k}\}_{k=0}^K$ with substep $\delta t=(t_{n+1}-t_n)/K$. Given low-order values $\{u_{n,k}\}$ and their continuous interpolant $\mathcal Iu(t)$, the error $\delta(t)=y(t)-\mathcal Iu(t)$ satisfies the **residual equation**

$$
\delta'(t)=f\bigl(t,\delta(t)+\mathcal Iu(t)\bigr)-\frac{\mathrm d}{\mathrm dt}\mathcal Iu(t),
\qquad\delta(0)=0 .
$$

The first term on the right has the same shape as the original equation and the second is a known function, so **the same low-order scheme** solves it directly, after which one updates $u_{n,k}^{\rm new}=u_{n,k}+\delta_k$. Repeating $J$ times gives the rate (Hairer 1978)

$$
O\bigl((\delta t)^{\min(J,K)+1}\bigr).
$$

The form $\min(J,K)$ says two things: the number of corrections $J$ cannot exceed the order the subgrid resolution $K$ can support, and conversely.

### Derivation: carrying it to FBSDEs

On $I_n$, starting from $Y_t=Y_{t_{n+1}}+\int_t^{t_{n+1}}f(s,X_s,Y_s,Z_s)\mathrm ds-\int_t^{t_{n+1}}Z_s\mathrm dW_s$, define the error processes $\delta Y_t=Y_t-\mathcal I_hY_t$ and $\delta Z_t=Z_t-\mathcal I_hZ_t$. They satisfy the **residual BSDE**

$$
\delta Y_t=\delta Y_{t_{n+1}}+\int_t^{t_{n+1}}F(s,X_s,\delta Y_s,\delta Z_s)\,\mathrm ds
-\int_t^{t_{n+1}}\delta Z_s\,\mathrm dW_s+E(t),
$$

$$
F(s,X_s,\delta Y_s,\delta Z_s)=f\bigl(s,X_s,\delta Y_s+\mathcal I_hY_s,\ \delta Z_s+\mathcal I_hZ_s\bigr),
\qquad
E(t)=\mathcal I_hY_{t_{n+1}}-\int_t^{t_{n+1}}\mathcal I_hZ_s\,\mathrm dW_s-\mathcal I_hY_t .
$$

The low-order scheme is the frozen-coefficient one-step method — exactly Euler, and exactly the $k=1$ member of the family of paper 8, as the paper says explicitly. For $k=K-1,\dots,0$,

$$
X^{k+1}=X^k+b(\tau_k,X^k)\delta t+\sigma(\tau_k,X^k)\Delta W_{\tau_k,\tau_{k+1}},
$$

$$
Z^k=\mathbb E^{X^k}_{\tau_k}\bigl[\bar Y^{k+1}(\Delta W_{\tau_k,\tau_{k+1}})^{\top}\bigr]\big/\delta t,
\qquad
Y^k=\mathbb E^{X^k}_{\tau_k}\bigl[\bar Y^{k+1}\bigr]+\delta t\cdot f(\tau_k,X^k,Y^k,Z^k).
$$

The corrector is: set $\delta Y^K=\delta Z^K=0$ and for $k=K-1,\dots,0$,

$$
\delta Z^k=\mathbb E^{X^k}_{\tau_k}\bigl[\delta\bar Y^{k+1}(\Delta W_{\tau_k,\tau_{k+1}})^{\top}\bigr]\big/\delta t
\;-\;Z^k\;+\;\nabla(\mathcal I_hY_{\tau_k})\,\sigma(\tau_k,X^k),
$$

$$
\delta Y^k=\mathbb E^{X^k}_{\tau_k}\bigl[\delta\bar Y^{k+1}\bigr]
+\delta t\Bigl(f\bigl(\tau_k,X^k,\delta Y^k+Y^k,\ \delta Z^k+Z^k\bigr)
+\mathcal L^0_{\tau_k,X^k}(\mathcal I_hY_{\tau_k})\Bigr).
$$

The two extra terms at the end come from the residual $E(t)$ and are supplied by the paper's Lemma 2.1, which is the generator theorem of paper 8:

$$
\left.\frac{\mathrm d\,\mathbb E^{X^k}_{\tau_k}[\mathcal I_hY_t]}{\mathrm dt}\right|_{t=\tau_k}
=\mathcal L^0_{\tau_k,X^k}(\mathcal I_hY_{\tau_k}),
\qquad
\left.\frac{\mathrm d\,\mathbb E^{X^k}_{\tau_k}[\mathcal I_hY_t(\Delta W_{\tau_k,\tau_{k+1}})^{\top}]}{\mathrm dt}\right|_{t=\tau_k}
=\nabla(\mathcal I_hY_{\tau_k})\,\sigma(\tau_k,X^k).
$$

> [!warning] The real cost of deferred correction sits in the interpolation operator
> The paper says so itself: "the high order accuracy of the deferred correction scheme depends heavily on the approximation quality of $\partial(\mathcal I_hY_t)/\partial t$, $\partial(\mathcal I_hY_t)/\partial x$ and $\partial^2(\mathcal I_hY_t)/\partial x^2$", because $\mathcal L^0$ contains a second spatial derivative. In other words, the multistep route puts the burden on the temporal stencil while deferred correction puts it on **the interpolation operator having to be twice differentiable and still accurate after differentiation**. This is not a free simplification but a transfer of cost.

The full algorithm: given $Y^N_i,Z^N_i$, for $n=N-1,\dots,0$ set $Y^{n,K}_i=Y^{n+1}_i$ and $Z^{n,K}_i=Z^{n+1}_i$; then for $j=1,\dots,J$ in turn (i) run the low-order scheme backwards from $k=K-1$ to obtain $Y^{n,k,[j]}_i,Z^{n,k,[j]}_i$, (ii) set $\delta Y^{K,[j]}_i=\delta Z^{K,[j]}_i=0$ and solve for $\delta Y^{k,[j]}_i,\delta Z^{k,[j]}_i$ with the same low-order scheme, (iii) update $Y^{n,k,[j+1]}_i=Y^{n,k,[j]}_i+\delta Y^{k,[j]}_i$ and $Z^{n,k,[j+1]}_i=Z^{n,k,[j]}_i+\delta Z^{k,[j]}_i$. Finally take $Y^n_i=Y^{n,0,[J]}_i$ and $Z^n_i=Z^{n,0,[J]}_i$. The spatial framework is the same as in paper 8: grid $\mathcal D_h=\{x_i\}$, density $h=\max_x\mathrm{dist}(x,\mathcal D_h)$, neighbour sets with $\#\mathcal D_{h,x}\le N_e$.

### Theorems

- The ODE rate $O((\delta t)^{\min(J,K)+1})$ is a restatement of Hairer's known result.
- For FBSDEs the paper **claims** order $K$ and verifies it numerically; this site found **no** FBSDE convergence theorem with explicit constants in the text.
- Stability is a numerical assertion: "stable and efficient, and is a $K$th order method at least for $K=1,\dots,4$". The point is that deferred correction has **no $k\le6$ barrier**: Test 2 goes all the way to $K=12$.

### Numerical experiments

**Test 1, decoupled system (5.1).** $\mathrm dX_t=\frac{1}{1+2e^{t+X_t}}\mathrm dt+\frac{e^{t+X_t}}{1+e^{t+X_t}}\mathrm dW_t$, $-\mathrm dY_t=\bigl(-\frac{2Y_t}{1+2e^{t+X_t}}-\frac12(\frac{Y_tZ_t}{1+e^{t+X_t}}-Y_t^2Z_t)\bigr)\mathrm dt-Z_t\mathrm dW_t$, $Y_T=\frac{e^{T+X_T}}{1+e^{T+X_T}}$, $x=1$; exact solution $Y_t=\frac{e^{t+X_t}}{1+e^{t+X_t}}$, $Z_t=\frac{(e^{t+X_t})^2}{(1+e^{t+X_t})^3}$. Rates over $N=4,6,8,10,12$:

| $K$      | 1     | 2     | 3     | 4     |
| -------- | ----- | ----- | ----- | ----- |
| $Y$ rate | 0.995 | 1.993 | 3.109 | 4.097 |
| $Z$ rate | 0.994 | 1.980 | 2.982 | 4.024 |

**Note how small $N$ is.** These are 4 to 12 large time steps, each subdivided into $K$ substeps — deferred correction moves the accuracy from "number of time steps" to "number of corrections per step".

**Test 1, coupled system (5.3).** $\mathrm dX_t=\frac{1}{1+e^{t+X_t}}\cdot\frac{1}{1+Y_t}\mathrm dt+Y_t\,\mathrm dW_t$, with the same backward equation and terminal condition, $x=0$. The rates at $K=1$ are $0.937/1.008$; higher $K$ follows the same order-$K$ pattern (the $K=2$ errors start at $3.666\times10^{-4}$), but this site transcribed the table only partially, so it is not tabulated here.

**Test 2, system (5.4).** $\mathrm dX_t=\sin(t+X_t)\mathrm dt+\frac{3}{10}\cos(t+X_t)\mathrm dW_t$, $-\mathrm dY_t=\bigl(\frac{3}{20}Y_tZ_t-\cos(t+X_t)(1+Y_t)\bigr)\mathrm dt-Z_t\mathrm dW_t$, $Y_T=\sin(T+X_T)$, $x=0.5$; exact solution $Y_t=\sin(t+X_t)$, $Z_t=\frac{3}{10}\cos^2(t+X_t)$. The paper concludes that "our deferred correction method is a $K$th order method ($K=1,2,\dots,12$), stable, efficient, with very high rate of convergence (up to 12)". Sampled rates are $0.906/1.001$ at $K=1$ and $2.193/2.031$ at $K=2$. The efficiency claim is that the error at $N=4$ with $K=2$ is **far smaller** than the error at $N=12$ with $K=1$, that is, the Euler scheme.

**What these experiments establish, and where they fall short.** They establish that deferred correction really does converge at order $K$ on smooth one-dimensional problems, including a coupled one, and that it **reaches $K=12$, far past the $k\le6$ barrier of the multistep family** — the most concrete evidence for this route over that one. Two shortfalls. First, all of this is numerical observation; the paper proves no FBSDE convergence theorem, so "why there is no barrier" remains blank in this paper (the framework of paper 47 provides the language for discussing it, but its specific conclusions for the deferred-correction family are unverified here). Second, this site transcribed only part of the tables, omitting the coupled test and the full error tables at high $K$, so the table above cannot support any cross-$K$ efficiency comparison; the only efficiency claim recorded is the paper's own comparison of $N=4,K=2$ against $N=12,K=1$.

### Relation to the others

This is the **methodological alternative** to paper 8: the same "Euler for the forward direction" trick and the same spatial framework, but with order coming from deferred-correction iterations rather than a $k$-step stencil. The $\min(J,K)+1$ rate structure and the absence of a $k\le6$ barrier are what set deferred correction apart from the multistep family, and they are part of what later motivates the strong-stability-preserving work of paper 68. Paper 35 is its direct sequel: deferred correction applied to the second-order equation class of paper 19, made explicit.

## 33: explicit $\theta$-schemes for mean-field backward equations

### The idea

The generator of a mean-field (McKean-Vlasov type) backward equation depends on the **law** of the solution itself, so every evaluation of the generator requires another layer of expectation. That layer is a new source of cost and the main obstacle to carrying existing schemes over: any implicit treatment turns "solve a nonlinear equation" into "solve a nonlinear equation containing an unknown distribution". This is the point of an "explicit" $\theta$-scheme — evaluate the mean-field generator at **already-known future levels**, so no nonlinear solve is needed per step.

### Setting

The standard form of a mean-field backward stochastic differential equation (Buckdahn-Djehiche-Li-Peng and Buckdahn-Li-Peng) is

$$
Y_t=\xi+\int_t^T \mathbb E'\bigl[f(s,X'_s,Y'_s,Z'_s,X_s,Y_s,Z_s)\bigr]\mathrm ds
-\int_t^T Z_s\,\mathrm dW_s,
$$

where $\mathbb E'$ is expectation over an independent copy $(X',Y',Z')$ of $(X,Y,Z)$.

### What could be verified

The abstract gives three things, confirmed verbatim here. First, what is proposed is a family of **explicit $\theta$-schemes**. Second, "we first prove a rigorous **stability result**", and on that basis "sharp error estimates" are given "showing that the proposed $\theta$-scheme admits a **second order convergence rate**". Third, "this seems to be the **first attempt to design high order numerical schemes for mean-field backward stochastic differential equations**" — earlier mean-field schemes were first order.

The **architecture** of the analysis deserves separate mention: stability first, then error estimates. That is the same order paper 47 later abstracts into a general framework, and the same order paper 63 follows. All three are papers by Weidong Zhao and Tao Zhou in SIAM J. Numer. Anal.

> [!note] What could be verified
> The paper has no preprint and SIAM blocks automated access from this environment, so nothing beyond the abstract is verified here. Specifically: the exact placement of $\theta$, how the law or the expectation $\mathbb E'$ is discretised (particle system? quadrature? nested expectation over Gaussians?), the treatment of $Z$, the exact form of the stability inequality, the hypotheses and constants in the error estimate, which values of $\theta$ give second order, and whether that second order holds for $Y$ alone or for $(Y,Z)$ jointly — none of these is reported here. For numerics, the abstract says "several numerical experiments are carried out to verify our theoretical results", but the test problems and observed orders are unverified.

### Relation to the others

This opens the mean-field branch of the research programme. Paper 61 is its explicit sequel — it describes itself as "one of our series works on numerical methods for mean-field forward-backward stochastic differential equations" — moving from mean-field **B**SDEs to mean-field **FB**SDEs and from $\theta$-schemes to multistep schemes, and it cites paper 33 directly. The $\theta$-scheme family is later shown by paper 47 to be a special case of its unified scheme. Papers 97 (DeepSPoC) and 108 attack the same problem class from the entirely different direction of deep learning and propagation of chaos.

## 35: deferred correction applied to second-order equations

### The idea

The corrector of paper 23 is **implicit** ($\delta Y^k$ appears on the right-hand side of its own equation), so each correction requires solving a nonlinear equation. If the correction can be made explicit, the cost per iteration becomes comparable to a single Euler sweep while the order still rises with each round. That is what this paper claims, and it claims it for the second-order (fully nonlinear) equation class of papers 16 and 19.

### Setting

The target equation class is the same as in paper 19 (see the [[en/computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control|second-order page]]): $\mathrm dX=b\,\mathrm dt+\sigma\,\mathrm dW$, $-\mathrm dY=f(t,X,Y,Z,\Gamma)\mathrm dt-Z\,\mathrm dW$, $\mathrm dZ=A\,\mathrm dt+\Gamma\,\mathrm dW$, whose solution represents a fully nonlinear parabolic equation.

### What could be verified

The abstract gives four points verbatim. First, "this is the second part of our series papers on deferred correction methods for forward backward stochastic differential equations", continuing paper 23 directly. Second, "we propose a class of **explicit** deferred correction schemes for 2FBSDEs" — "explicit" being the key difference from paper 23. Third, the key feature is using "the simple Euler scheme as the initialisation", after which "a simple deferred correction iterative scheme" yields "very high accuracy" approximations. Fourth, "in each iteration, the computational complexity is always comparable to the Euler solver".

The paper also gives a self-assessment of scope worth keeping verbatim: "we believe the schemes proposed in this work are promising for solving 2FBSDEs in **moderate dimensions**" — an explicit acknowledgement that the method does not reach high dimension.

> [!note] What could be verified
> The paper is paywalled with no preprint, so this site read only the abstract and reference list. The exact form of the residual 2FBSDE, the relation between the number of correction sweeps $J$ and the order attained, and how $\Gamma$ and $A$ are handled in the corrector are all unverified here. The abstract claims "very high accuracy" and a constant cost per round but **gives no order, hypotheses or constants**, and this site read no theorem. For numerics, the abstract says "numerical examples are given to show the effectiveness of the proposed schemes", but the test problems and observed orders are unverified.

### Relation to the others

It sits at the intersection of paper 23 (the deferred-correction methodology) and papers 19 and 16 (the second-order equation class). The first author, Jie Yang, is also an author of papers 47 and 61, and the general stability and consistency framework of paper 47 is exactly the language that accommodates both the multistep family and the deferred-correction and $\theta$ families. The claim that each iteration costs about as much as an Euler solve is the structural selling point of the deferred-correction family against the multistep one, whose $k$-step schemes cost roughly $k$ Euler steps per level and additionally need $k$ starting values.

## 61: an explicit multistep scheme for mean-field forward-backward systems

### The idea

Paper 33 used $\theta$-schemes to push mean-field BSDEs to second order. The natural way to go higher is to swap in a multistep stencil — exactly what paper 8 did in the non-mean-field case. The difficulty is still that layer of expectation over the law, $\mathbb E'$: it has to be discretised into something computable, and that discretisation must not destroy the temporal order.

### Setting

The McKean-Vlasov type forward-backward system is

$$
X_t=X_0+\int_0^t\mathbb E'\bigl[b(s,X'_s,X_s)\bigr]\mathrm ds+\int_0^t\mathbb E'\bigl[\sigma(s,X'_s,X_s)\bigr]\mathrm dW_s,
$$

$$
Y_t=\mathbb E'\bigl[\varphi(X'_T,X_T)\bigr]+\int_t^T\mathbb E'\bigl[f(s,X'_s,Y'_s,Z'_s,X_s,Y_s,Z_s)\bigr]\mathrm ds-\int_t^TZ_s\,\mathrm dW_s .
$$

### What could be verified

The abstract and keywords confirm verbatim that the scheme is an **explicit multistep** one, that it is "easy to implement", that it has "high order rate of convergence", and that "rigorous error estimates of the proposed multistep scheme are given". The keywords are "mean-field forward backward stochastic differential equations; explicit multistep scheme; error estimates", with MSC codes 60H35, 65C20 and 60H10.

The composition of the reference list is itself verifiable evidence of what this paper stands on: Zhao-Chen-Peng (SISC 28 (2006) 1563-1581, the origin of the Lagrange-interpolation multistep idea), Zhao-Zhang-Ju (SINUM 48 (2010) 1369-1394, the stable multistep scheme; and NMTMA 9 (2016) 262-288, multistep schemes for decoupled FBSDEs), papers 8, 33, 35, 47 and 25, Sun-Yang-Zhao (NMTMA 10 (2017) 798-828, Itô-Taylor schemes for mean-field SDEs), Kloeden-Shardlow (SISC 39 (2017) A2784-A2807, a Gauss-quadrature method for one-dimensional mean-field SDEs), and Buckdahn-Djehiche-Li-Peng and Buckdahn-Li-Peng. **The presence of Kloeden-Shardlow makes a quadrature treatment of $\mathbb E'$ a reasonable guess, but it is only a guess and this site asserts nothing.**

> [!note] What could be verified
> The body is paywalled with no preprint, so this site read only the landing page (which carries the full 38-item reference list), OpenAlex and Semantic Scholar records. The multistep coefficients, the number of steps $k$ used, and the single most important ingredient of any mean-field scheme — how $\mathbb E'$ is actually discretised — are all unverified here. The specific order, norms and hypotheses behind the claimed high order, and whether a root-condition window like the $1\le k\le6$ of paper 8 appears, are likewise unverified. For numerics, the abstract says "numerical experiments are carried out to illustrate the efficiency and accuracy of the proposed scheme", but the test problems and observed orders are unverified.
>
> One bibliographic correction: entries commonly give the pages as 40 (2022) 519-543, whereas the Global Science Press landing page and OpenAlex both give **517-540**, which is what this page uses.

### Relation to the others

A direct sequel to paper 33, which it cites: 33 reaches second order for mean-field BSDEs with $\theta$-schemes, 61 claims high order for mean-field FBSDEs with multistep schemes. Its methodological parent is paper 8. The second author, Jie Yang, is shared with papers 35 and 47, and the unified framework of paper 47 is cited here, which suggests its error analysis most likely follows the same stability-first architecture — again only a suggestion.

## 68: once stability is the pivot, design schemes to maximise it

### The idea

By 2022 the group had two things in hand: a family of high-order multistep schemes whose stability was characterised only by an empirically observed root-condition window (for instance the $1\le k\le6$ of paper 8), and a general stability, consistency and convergence framework (paper 47) that proved a mean-square Lax equivalence theorem but said nothing about how to **build** schemes with good stability. The natural step connecting them: since stability is the necessary and sufficient other half of convergence, make stability the **design objective** and, among the coefficient sets meeting a required consistency order, pick the ones with the best stability. That is exactly the strong-stability-preserving idea from the numerics of hyperbolic conservation laws, where one maximises the CFL-like coefficient for which a convex-combination (contractivity) property survives.

### What could be verified

The abstract states that the authors first perform a comprehensive analysis of a general type of multistep scheme for FBSDEs, on that basis present **new sufficient conditions on the coefficients** so the associated schemes are stable and enjoy a certain order of consistency, and then propose a **practical way to design high-order strong-stability-preserving multistep schemes**. The appendix titles on the Springer page are also verifiable: Appendix A, "Additional Optimal SSPM Schemes", contains Table 7, "SSPM schemes with uniform time partition (part 2)", giving "the coefficients of the optimal SSPM schemes with **order up to 5**". So the paper's practical output is a table of optimal coefficients for orders one through five.

The scheme template analysed has the same shape as the family unified in paper 47: for a $k$-step method with $\mathbb E_n[\cdot]=\mathbb E[\cdot\mid\mathcal F_{t_n}]$,

$$
\sum_{i=0}^{k}\alpha_i\,\mathbb E_n\bigl[Y^{n+i}\bigr]
=\Delta t\sum_{i=0}^{k}\beta_i\,\mathbb E_n\bigl[f(t_{n+i},X^{n+i},Y^{n+i},Z^{n+i})\bigr],
$$

with a paired recursion for $Z^n$ built from $\mathbb E_n[Y^{n+i}\Delta W]/\Delta t$ terms. The strong-stability question is then: for which coefficient vectors $(\alpha_i),(\beta_i)$ can the scheme be rewritten as a **convex combination** of backward-Euler-like steps, so any monotonicity or contractivity property of the base step is inherited, with the largest possible step-size coefficient? That is the exact analogue of the strong-stability-preserving linear multistep theory of Lenferink (_Numer. Math._ 55 (1989) 213-223; _Math. Comp._ 56 (1991) 177-199) and of Spiteri and Ruuth (_SIAM J. Numer. Anal._ 40 (2002) 469-491), all three of which are cited in the paper.

> [!note] What could be verified
> The positioning, the three-step statement in the abstract and the **existence** of the appendix table with orders one through five are all confirmable. The scheme template above is reconstructed from the same-family schemes of papers 47 and 8; **it is not transcribed from this paper**. The exact form of the sufficient conditions, the precise definition of "strongly stable" in the FBSDE sense (the most important technical definition in the paper, and one the abstract does not recover), the stability functional used, the optimisation problem behind "optimal", and the numerical coefficient values are all unverified here; nor can the spatial discretisation and the treatment of conditional expectations be inferred from the abstract.
>
> For numerics, the abstract says "numerical experiments are carried out to demonstrate the **strong stability** of our SSPM schemes", and the Springer page shows at least one figure. The test problems and observed orders are unverified here; from the abstract's emphasis these experiments are more likely stability demonstrations than order-verification tables, but that too is only conjecture.

### Relation to the others

A direct successor to papers 8 and 47, both of which it cites. It also cites paper 63 (another 2022 SINUM article from the same group), as well as Tang and Xiong (_IMA J. Numer. Anal._ 42 (2022) 1789-1805, stability analysis of general multistep methods for Markovian BSDEs) and Chassagneux (_SIAM J. Numer. Anal._ 52 (2014) 2815-2836, linear multistep schemes for BSDEs) — two external works that paper 47 is also in dialogue with. The first author, Shuixin Fang, appears here for the first time and later becomes the principal collaborator on the deep-learning papers 86, 93, 96, 100 and 108; **this paper is his classical-numerical-analysis entry point, and the martingale deep-learning papers inherit its stability instincts.**

## How the seven relate

| No. | Route to higher order             | Problem class           | Basis of stability                          | Verification here      |
| --- | --------------------------------- | ----------------------- | ------------------------------------------- | ---------------------- |
| 8   | differentiate into reference ODEs | coupled FBSDE           | root-condition window $k\le6$               | full text, equations   |
| 18  | multistep plus jump handling      | FBSDE with jumps        | unverified                                  | abstract and reference list |
| 23  | deferred correction               | FBSDE                   | numerical assertion, no barrier             | full text, equations   |
| 33  | explicit $\theta$-schemes         | mean-field BSDE         | proved in the paper (unverified here)       | abstract               |
| 35  | explicit deferred correction      | second-order FBSDE      | unverified                                  | abstract and reference list |
| 61  | explicit multistep                | mean-field FBSDE        | unverified                                  | abstract and reference list |
| 68  | design backwards from stability   | general multistep FBSDE | new sufficient conditions plus optimisation | abstract and appendix titles |

The shape of the thread is worth summarising: **construct schemes first (papers 8 through 35), then unify the analysis (paper 47), then design backwards from the analysis (paper 68).** Paper 68 is possible precisely because paper 47 had established stability as the necessary and sufficient other half of convergence, which turns "design a stable scheme" into an optimisation problem with a clear target rather than trial and error.

The other shape is the division of labour between the two order-raising mechanisms. The multistep route puts the burden on the **temporal stencil**: $k$ starting values, $k$ future levels, one root-condition window. Deferred correction puts it on the **interpolation operator**, which must be smooth enough to differentiate twice because $\mathcal L^0(\mathcal I_hY)$ appears in the residual equation. Their measured ceilings differ too — $k=6$ for paper 8 (in practice lower in double precision) against $K=12$ for paper 23.

## Coverage check

| Item                                                              | Paper      | Status                                                       |
| ----------------------------------------------------------------- | ---------- | ------------------------------------------------------------ |
| Coupled FBSDE, prior state of the art and the question posed      | 8          | complete                                                     |
| Generator theorem, derivative weights, two reference ODEs, five schemes | 8    | complete derivation                                          |
| Truncation-error split, error balancing, root-condition window    | 8          | complete, with the missing convergence theorem flagged       |
| Three experiment groups and their limitations                     | 8          | rate tables for Examples 1 through 3; Example 4 not transcribed |
| System with jumps and the three abstract claims                   | 18         | only what the abstract supports                              |
| Deferred correction: residual equation, low-order scheme, corrector | 23       | complete derivation with both generator identities           |
| Two experiment groups and $K=12$                                  | 23         | Test 1 decoupled complete; coupled and Test 2 partly transcribed |
| Mean-field BSDE form, "first high order" and the second-order claim | 33       | only what the abstract supports                              |
| Four abstract claims for explicit deferred correction and the dimension caveat | 35 | only what the abstract supports                          |
| Mean-field FBSDE form, error-estimate claim and reference composition | 61     | only what the abstract and reference list support; includes a page-number correction |
| Strong-stability design idea, template and existence of the coefficient table | 68 | only what the abstract and appendix titles support        |

## Sources for this page

- W. Zhao, Y. Fu, and T. Zhou, [_New kinds of high-order multistep schemes for coupled forward backward stochastic differential equations_](https://doi.org/10.1137/130941274), SIAM J. Sci. Comput. 36(4) (2014), pp. A1731-A1751 (preprint [arXiv:1310.5307](https://arxiv.org/abs/1310.5307)).
- Y. Fu, W. Zhao, and T. Zhou, [_Multistep schemes for forward backward stochastic differential equations with jumps_](https://doi.org/10.1007/s10915-016-0212-y), J. Sci. Comput. 69(2) (2016), pp. 651-672.
- T. Tang, W. Zhao, and T. Zhou, [_Deferred correction methods for forward backward stochastic differential equations_](https://doi.org/10.4208/nmtma.2017.s02), Numer. Math. Theor. Meth. Appl. 10(2) (2017), pp. 222-242.
- Y. Sun, W. Zhao, and T. Zhou, [_Explicit theta-schemes for mean-field backward stochastic differential equations_](https://doi.org/10.1137/17M1161944), SIAM J. Numer. Anal. 56(4) (2018), pp. 2672-2697.
- J. Yang, W. Zhao, and T. Zhou, [_Explicit deferred correction methods for second-order forward backward stochastic differential equations_](https://doi.org/10.1007/s10915-018-00896-w), J. Sci. Comput. 79(3) (2019), pp. 1409-1432.
- Y. Sun, J. Yang, W. Zhao, and T. Zhou, [_An explicit multistep scheme for mean-field forward-backward stochastic differential equations_](https://doi.org/10.4208/jcm.2011-m2019-0205), J. Comput. Math. 40(4) (2022), pp. 517-540.
- S. Fang, W. Zhao, and T. Zhou, [_Strong stability preserving multistep schemes for forward backward stochastic differential equations_](https://doi.org/10.1007/s10915-023-02111-x), J. Sci. Comput. 94(3) (2023), 53.
- External sources used for cross-checking: W. Zhao, G. Zhang, and L. Ju, [_A stable multistep scheme for solving backward stochastic differential equations_](https://doi.org/10.1137/09076979X), SIAM J. Numer. Anal. 48(4) (2010), pp. 1369-1394; W. Zhao, L. Chen, and S. Peng, [_A new kind of accurate numerical method for backward stochastic differential equations_](https://doi.org/10.1137/05063341X), SIAM J. Sci. Comput. 28(4) (2006), pp. 1563-1581; L. Teng, A. Lapitckii, and M. Günther, _A multi-step scheme based on cubic spline for solving backward stochastic differential equations_, [arXiv:1809.00324](https://arxiv.org/abs/1809.00324) (the source from which this site transcribed the Zhao-Zhang-Ju reference equations and stability windows).
