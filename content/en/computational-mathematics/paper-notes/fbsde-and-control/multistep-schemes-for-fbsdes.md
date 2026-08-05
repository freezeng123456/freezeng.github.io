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
> Papers **8** (_SIAM J. Sci. Comput._ 36(4), 2014), **18** (_J. Sci. Comput._ 69(2), 2016), **23** (_Numer. Math. Theor. Meth. Appl._ 10(2), 2017), **33** (_SIAM J. Numer. Anal._ 56(4), 2018), **35** (_J. Sci. Comput._ 79, 2019), **61** (_J. Comput. Math._ 40, 2022) and **68** (_J. Sci. Comput._ 94:53, 2023). Most appear in SIAM or Springer journals with no preprint, so several technical details could only be confirmed from abstracts and reference lists; those places are marked.

![How multistep schemes raise the backward accuracy](assets/diagrams/tao-zhou-papers/en/fbsde-multistep.svg)

## 8: moving the entire high-order requirement onto the backward direction

### The question the paper poses

The coupled Markovian forward-backward system is

$$
X_t=X_0+\int_0^t b(s,X_s,Y_s,Z_s)\,\mathrm ds+\int_0^t\sigma(s,X_s,Y_s,Z_s)\,\mathrm dW_s,
$$

$$
Y_t=\xi+\int_t^T f(s,X_s,Y_s,Z_s)\,\mathrm ds-\int_t^T Z_s\,\mathrm dW_s,
\qquad \xi=\varphi(X_T).
$$

The paper describes the state of the art as follows: most existing schemes are Euler-type with rate $1/2$; the high-order methods for decoupled systems rely "on the high order approaches for both the forward SDE and the high order scheme for the backward SDE", and high-order approaches for the forward equation "require large amounts of computations and are often difficult to be applied". In the coupled case the forward coefficients depend on the backward unknowns, so "it seems not easy to design high-order (yet efficient) numerical schemes".

The paper poses and answers affirmatively the explicit question: **can one still expect high-order accurate numerical solutions of the backward equation if the Euler method is used to solve the forward equation?**

### The theorem licensing frozen forward coefficients

For $X_s$ solving $\mathrm dX_s=b\,\mathrm ds+\sigma\,\mathrm dW_s$ the generator is

$$
\mathcal A_t^x g(t,x)=\lim_{s\downarrow t}\frac{\mathbb E_t^x[g(s,X_s)]-g(t,x)}{s-t},
\qquad
\mathcal L^0_{t,x}=\frac{\partial}{\partial t}+\sum_i b_i\frac{\partial}{\partial x_i}
+\frac12\sum_{i,j}(\sigma\sigma^{\top})_{i,j}\frac{\partial^2}{\partial x_i\partial x_j}.
$$

The key theorem states that if $f\in C^{1,2}$ and $\mathbb E^{x_0}_{t_0}[\mathcal L^0_{t,X_t}f(t,X_t)]<\infty$, then

$$
\frac{\mathrm d\,\mathbb E^{x_0}_{t_0}[f(t,X_t)]}{\mathrm dt}
=\mathbb E^{x_0}_{t_0}\bigl[\mathcal A_t^{X_t}f(t,X_t)\bigr],
$$

and that **at $t=t_0$** this derivative is unchanged if $X_t$ is replaced by **any** diffusion $\bar X_t$ whose coefficients merely match at the left endpoint:

$$
\bar b(t_0,\bar X_{t_0};t_0,x_0)=b(t_0,x_0),
\qquad
\bar\sigma(t_0,\bar X_{t_0};t_0,x_0)=\sigma(t_0,x_0)
\quad\Longrightarrow\quad
\left.\frac{\mathrm d\,\mathbb E[f(t,X_t)]}{\mathrm dt}\right|_{t_0}
=\left.\frac{\mathrm d\,\mathbb E[f(t,\bar X_t)]}{\mathrm dt}\right|_{t_0}.
$$

This is the licence for replacing the true forward diffusion by a **frozen-coefficient** (Euler) one; the paper notes explicitly that one can simply take $\bar b(s,\bar X_s;t_0,x_0)=b(t_0,x_0)$ and $\bar\sigma(s,\bar X_s;t_0,x_0)=\sigma(t_0,x_0)$ for all $s\in[t_0,t]$.

**The general lesson deserves separate emphasis:** when the two directions of a coupled system carry different accuracy requirements, first check whether the error from the cheaper direction actually reaches the quantity of interest. Here it does not — at the level of the one-step derivative approximation the forward diffusion need only match at the left endpoint.

### The derivative weights come from moment conditions

For $u\in C_b^{k+1}$ and nodes $t_0<t_1<\cdots<t_k$ with $\Delta t_i=t_i-t_0$, the weights satisfy

$$
\sum_{i=0}^{k}\alpha_{k,i}\frac{(\Delta t_i)^j}{j!}=\delta_{j1},
\qquad j=0,1,\dots,k .
$$

These conditions turn "approximate a first derivative from $k+1$ levels" into a linear system of the same shape as the coefficient conditions of a multistep ODE method.

### The order limit comes from a root condition

The window $1\le k\le6$ is not a vague empirical observation but a zero-stability check anyone can redo. Applying the same weights to the deterministic ODE $\mathrm dY/\mathrm dt=f(t,Y)$ gives

$$
\alpha_{k,0}Y^n+\sum_{j=1}^{k}\alpha_{k,j}Y^{n+j}=f(t_n,Y^n),
$$

with characteristic polynomial

$$
P(\lambda)=\alpha_{k,0}\lambda^{k}+\sum_{j=1}^{k}\lambda^{k-j}=0,
$$

subject to the root condition $|\lambda_{k,j}|\le1$, simple where equality holds. The maximum root moduli the paper reports, excluding the common root $1.0$, are

| $k$              | 2      | 3      | 4      | 5      | 6      | 7          | 8          |
| ---------------- | ------ | ------ | ------ | ------ | ------ | ---------- | ---------- |
| max root modulus | 0.3333 | 0.4264 | 0.5608 | 0.7087 | 0.8633 | **1.0222** | **1.1839** |

**The modulus crosses $1$ between $k=6$ and $k=7$, so the scheme is unstable for $k\ge7$**, which is why the tabulation stops at $k=6$. What paper 8 lacks is not this computation but a **stochastic** stability theory; that arrives with paper 47 (see the [[en/computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes|stability theory page]]), where the root condition is precisely the stability half of an equivalence.

### Distinguishing it from the other multistep route

Two constructions both called "multistep" need separating here, or their conclusions get attached to the wrong scheme.

- **Interpolate then integrate** (Zhao, Zhang and Ju, _SIAM J. Numer. Anal._ 48(4) 2010): Lagrange-interpolate the integrand of the reference integral identities over several future time levels, then integrate, producing Newton-Cotes-type weights of the form $h\sum_jb_j\mathbb E_i[f(t_{i+j},Y_{i+j},Z_{i+j})]$. Its stability windows **differ between the two directions**: the reference equation for $Y$ is stable only for $K_y\in\{1,\dots,7,9\}$ — note that $K_y=8$ is excluded — while the one for $Z$ is stable only for $K_z\in\{1,2,3\}$.
- **Differentiate into reference ODEs** (paper 8): it does **not** interpolate the integrand and integrate. It differentiates the reference integral identities in $t$, turning them into two reference ODEs, then discretises with the derivative-approximation weights above. Hence a single window $1\le k\le6$.

One common misattribution is worth correcting in passing: the multistep, interpolate-then-integrate construction is due to Zhao-Zhang-Ju (2010), whereas Zhao, Chen and Peng (_SIAM J. Sci. Comput._ 28(4) 2006) is the origin of the **$\theta$-scheme**, which the multistep scheme extends. The two should not be conflated.

## 18, 23 and 35: three different routes to higher order

- **18 (multistep schemes with jumps)** carries the same idea to jump processes. Jumps add a compensated Poisson-measure integral to the backward equation, so the conditional expectations associated with jumps need separate treatment.
- **23 (deferred correction)** raises order by deferred correction rather than multistep interpolation: compute a low-order approximation, then solve a sequence of correction equations that lift the order step by step. Unlike multistep methods it needs no extra starting values, at the price of one extra sweep per order gained.
- **35 (explicit deferred correction for second-order FBSDEs)** applies deferred correction to the second-order (fully nonlinear) setting, where a second-order process $\Gamma$ joins $Y$ and $Z$ (see the [[en/computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control|second-order page]]).

## 33 and 61: the mean-field case

Mean-field backward equations are of McKean-Vlasov type, with a generator depending on the **law** of the solution:

$$
Y_t=\xi+\int_t^T \mathbb E'\bigl[f(s,X'_s,Y'_s,Z'_s,X_s,Y_s,Z_s)\bigr]\mathrm ds
-\int_t^T Z_s\,\mathrm dW_s,
$$

where $\mathbb E'$ is expectation over an independent copy $(X',Y',Z')$ of $(X,Y,Z)$.

The abstract of paper 33 states that this "seems to be the **first attempt to design high order numerical schemes for mean-field backward stochastic differential equations**" — earlier mean-field schemes were first order. The construction is a class of **explicit** theta-schemes, where "explicit" means the mean-field generator is evaluated at already-known future levels so no nonlinear solve in $Y^n$ is required. Paper 61 gives an explicit multistep scheme for the mean-field forward-backward system, carrying the route of paper 8 into the mean-field setting.

> [!note] What could be verified
> The texts of papers 33 and 61 could not be retrieved (no preprint, publisher blocking automated access). The problem setting above and the claims of "first high order" and "explicit theta-schemes" are confirmable from the abstracts; the exact placement of the theta parameters, how the law or the expectation $\mathbb E'$ is discretised (particle system, quadrature or nested expectation), and the treatment of $Z$ are unverified here.

## 68: once stability is the pivot, design schemes to maximise it

The position of paper 68 is clear. By 2022 the group had, first, a family of high-order multistep schemes whose stability was characterised only through an empirically observed root-condition window such as $1\le k\le6$ in paper 8, and second, a general stability, consistency and convergence framework (paper 47) proving a mean-square Lax equivalence theorem but saying nothing about how to **build** schemes with good stability. This paper closes the loop by importing **strong stability preserving** design — the Gottlieb-Shu-Tadmor and Gottlieb-Ketcheson-Shu technology from hyperbolic conservation laws, where one maximises the CFL-like coefficient for which a convex-combination (contractivity) property survives — into the FBSDE setting.

The abstract states that the authors first perform a comprehensive analysis of a general type of multistep scheme for FBSDEs, and on that basis present **new sufficient conditions on the coefficients** so the associated schemes are stable and enjoy a certain order of consistency, then propose a **practical way to design high-order strong-stability-preserving multistep schemes**. An appendix supplies a table of optimal coefficient sets with **orders up to 5** under uniform time partitions.

The scheme template analysed has the same shape as the family unified in paper 47: for a $k$-step method with $\mathbb E_n[\cdot]=\mathbb E[\cdot\mid\mathcal F_{t_n}]$,

$$
\sum_{i=0}^{k}\alpha_i\,\mathbb E_n\bigl[Y^{n+i}\bigr]
=\Delta t\sum_{i=0}^{k}\beta_i\,\mathbb E_n\bigl[f(t_{n+i},X^{n+i},Y^{n+i},Z^{n+i})\bigr],
$$

with a paired recursion for $Z^n$ built from $\mathbb E_n[Y^{n+i}\Delta W]/\Delta t$ terms. The strong-stability question is then: for which coefficient vectors $(\alpha_i),(\beta_i)$ can the scheme be rewritten as a **convex combination** of backward-Euler-like steps, so any monotonicity or contractivity property of the base step is inherited, with the largest possible step-size coefficient? That is the exact analogue of the strong-stability-preserving linear multistep theory of Lenferink and of Spiteri and Ruuth, both cited in the paper.

> [!note] What could be verified
> The positioning, the three-step statement in the abstract and the **existence** of the appendix table with orders one through five are confirmable. The exact form of the sufficient conditions, the stability functional used, the optimisation problem behind "optimal", and the numerical coefficient values are unverified here.

## How the seven relate

| No. | Route to higher order             | Problem class           | Basis of stability                          |
| --- | --------------------------------- | ----------------------- | ------------------------------------------- |
| 8   | differentiate into reference ODEs | coupled FBSDE           | root-condition window $k\le6$               |
| 18  | multistep plus jump handling      | FBSDE with jumps        | same family                                 |
| 23  | deferred correction               | FBSDE                   | same family                                 |
| 33  | explicit theta-schemes            | mean-field BSDE         | proved in the paper                         |
| 35  | deferred correction               | second-order FBSDE      | same family                                 |
| 61  | explicit multistep                | mean-field FBSDE        | same family                                 |
| 68  | design backwards from stability   | general multistep FBSDE | new sufficient conditions plus optimisation |

The shape of the thread is worth summarising: **construct schemes first (papers 8 through 35), then unify the analysis (paper 47), then design backwards from the analysis (paper 68).** Paper 68 is possible precisely because paper 47 had established stability as the necessary and sufficient other half of convergence, which turns "design a stable scheme" into an optimisation problem with a clear target rather than trial and error.

## Coverage check

| Item                                            | Paper      | Status                                                                                     |
| ----------------------------------------------- | ---------- | ------------------------------------------------------------------------------------------ |
| Coupled FBSDE and the question posed            | 8          | system, prior state of the art, explicit question                                          |
| Generator theorem and left-endpoint matching    | 8          | generator, theorem, licence, general lesson                                                |
| Moment conditions for the derivative weights    | 8          | the linear conditions                                                                      |
| The root-condition window                       | 8          | the computed window and where its stochastic theory arrives                                |
| Jumps, deferred correction, second order        | 18, 23, 35 | each route and its cost                                                                    |
| Mean-field form and the first-high-order claim  | 33         | McKean-Vlasov form, meaning of explicit (limited verification)                             |
| Strong-stability design and the scheme template | 68         | origin, three-step statement, template, convex-combination question (limited verification) |

## Sources for this page

- W. Zhao, Y. Fu, and T. Zhou, [_New kinds of high-order multistep schemes for coupled forward backward stochastic differential equations_](https://doi.org/10.1137/130941274), SIAM J. Sci. Comput. 36(4) (2014), pp. A1731-A1751 (preprint [arXiv:1310.5307](https://arxiv.org/abs/1310.5307)).
- Y. Fu, W. Zhao, and T. Zhou, [_Multistep schemes for forward backward stochastic differential equations with jumps_](https://doi.org/10.1007/s10915-016-0212-y), J. Sci. Comput. 69(2) (2016), pp. 651-672.
- T. Tang, W. Zhao, and T. Zhou, [_Deferred correction methods for forward backward stochastic differential equations_](https://doi.org/10.4208/nmtma.2017.s02), Numer. Math. Theor. Meth. Appl. 10(2) (2017), pp. 222-242.
- Y. Sun, W. Zhao, and T. Zhou, [_Explicit theta-schemes for mean-field backward stochastic differential equations_](https://doi.org/10.1137/17M1161944), SIAM J. Numer. Anal. 56(4) (2018), pp. 2672-2697.
- J. Yang, W. Zhao, and T. Zhou, [_Explicit deferred correction methods for second-order forward backward stochastic differential equations_](https://doi.org/10.1007/s10915-018-00896-w), J. Sci. Comput. 79 (2019), pp. 1409-1432.
- Y. Sun, J. Yang, W. Zhao, and T. Zhou, [_An explicit multistep scheme for mean-field forward-backward stochastic differential equations_](https://doi.org/10.4208/jcm.2011-m2019-0205), J. Comput. Math. 40 (2022), pp. 519-543.
- S. Fang, W. Zhao, and T. Zhou, [_Strong stability preserving multistep schemes for forward backward stochastic differential equations_](https://doi.org/10.1007/s10915-023-02111-x), J. Sci. Comput. 94 (2023), 53.
