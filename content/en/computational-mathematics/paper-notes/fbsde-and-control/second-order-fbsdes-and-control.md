---
title: Second-Order FBSDEs and Stochastic Control
description: Papers 16, 19, 25, 26, 41, 50 and 51 - fully nonlinear equations, sparse grids and control iterations
lang: en
translation: computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control
tags:
  - paper-notes
  - stochastic-optimal-control
  - fully-nonlinear-equations
---

> [!note] Coverage of this page
> Papers **16** (_Commun. Comput. Phys._ 18(5), 2015), **19** (_Commun. Comput. Phys._ 21(3), 2017), **25** (_Discrete Contin. Dyn. Syst. Ser. B_ 22(9), 2017), **26** (_SIAM J. Numer. Anal._ 55(6), 2017), **41** (_Numer. Math. Theor. Meth. Appl._ 13(2), 2020), **50** (_J. Sci. Comput._ 85(2), 2020) and **51** (_SIAM J. Control Optim._ 58(6), 2020).
>
> Of these, **19** (full arXiv preprint), **25** (full arXiv preprint), **26** (full author-accepted manuscript), **50** (full arXiv preprint) and **51** (published text) were checked equation by equation, so this page gives their complete derivations, theorems and numerical experiments. The bodies of **16** and **41** could not be obtained — neither has a preprint, and the Global Science Press PDF endpoints are unreachable or return 403 to direct download — so for those two only what the abstracts and reference lists support is written here. Where scheme details are inferred from their sibling papers, that is flagged in place.

## Why a "second-order" FBSDE is needed

A first-order forward-backward system corresponds through the nonlinear Feynman-Kac relation to a **semilinear** parabolic equation: $Y_t$ gives the solution and $Z_t$ gives the ($\sigma$-weighted) gradient. Second derivatives can only enter **linearly**, through the quadratic-variation term of Itô's formula; there is no other doorway. So covering **fully nonlinear** equations, where the Hessian enters nonlinearly, requires one more process in the probabilistic representation.

The second-order backward stochastic differential equation of Cheridito, Soner, Touzi and Victoir is designed for exactly that. Its definition: for $\mathrm dX_t=\mu(X_t)\mathrm dt+\sigma(X_t)\mathrm dW_t$ and a quadruple $(Y,Z,\Gamma,A)$ valued in $\mathbb R\times\mathbb R^d\times\mathbb S^d\times\mathbb R^d$ ($\mathbb S^d$ the real symmetric $d\times d$ matrices),

$$
\mathrm dY_t=f(t,X_t,Y_t,Z_t,\Gamma_t)\,\mathrm dt+Z_t^{\top}\circ\mathrm dX_t,
\qquad
\mathrm dZ_t=A_t\,\mathrm dt+\Gamma_t\,\mathrm dX_t,
\qquad
Y_T=g(X_T),
$$

where $Z_t^\top\circ\mathrm dX_t$ is Fisk-Stratonovich integration, related to Itô integration by $Z_t^\top\circ\mathrm dX_t=Z_t^\top\mathrm dX_t+\frac12\operatorname{Tr}[\Gamma_t\sigma\sigma^\top]\mathrm dt$. The associated partial differential equation is

$$
-v_t(t,x)+f\bigl(t,x,v(t,x),Dv(t,x),D^2v(t,x)\bigr)=0,\qquad v(T,x)=g(x),
$$

with the representation $Y_t=v(t,X_t)$, $Z_t=Dv(t,X_t)$, $\Gamma_t=D^2v(t,X_t)$, $A_t=\mathcal L\,Dv(t,X_t)$, where $\mathcal L\varphi=\varphi_t+\frac12\operatorname{Tr}[D^2\varphi\,\sigma\sigma^\top]$ is the Dynkin operator with the drift term removed. The theory requires $f$ to be **decreasing in $\Gamma$** — that is the ellipticity or parabolicity condition — and a viscosity comparison principle for the PDE before the 2BSDE can conversely pin down $v$ uniquely.

**There is a design freedom here that matters enormously for numerics.** Right after stating the PDE, CSTV note that its form does **not** depend on $\mu$ and $\sigma$, so in principle one could restrict to $\mu\equiv0$, $\sigma\equiv I_d$; but keeping the freedom to choose $\mu$ and $\sigma$ "provides additional flexibility in the design of the Monte Carlo schemes". **That sentence is the theoretical basis of paper 16's headline claim**, and it dovetails exactly with Theorem 2 of paper 8 (the one-sided derivative sees only the left-endpoint values of $\bar b,\bar\sigma$): in the semilinear case the forward SDE is dictated by the PDE; in the fully nonlinear case it is a free design parameter.

Papers 16, 19 and 35 write the system slightly differently from CSTV, and the two differences must be kept in mind when matching formulas across the two literatures:

1. **Driving noise.** CSTV drives $Z$ by $\mathrm dX$; paper 19 drives $Z$ by $\mathrm dW$. Since $\mathrm dX=\mu\,\mathrm dt+\sigma\,\mathrm dW$, the two $\Gamma$'s differ by a factor of $\sigma$ and the two $A$'s absorb the drift term.
2. **Weighting by $\sigma$.** CSTV sets $Z=Dv$, $\Gamma=D^2v$; the representation theorem of paper 19 sets $Z_t=(\nabla_xu\,\sigma)(t,X_t)$, $\Gamma_t=(\nabla_x(\nabla_xu\,\sigma)\sigma)(t,X_t)$, $A_t=(\mathcal L(\nabla_xu\,\sigma))(t,X_t)$. **This preserves the same $\sigma$-weighting convention as the first-order case, which is exactly why the multistep stencils of paper 8 carry over verbatim.**

## 16: probabilistic high-order schemes for fully nonlinear parabolic PDEs

### The idea

Probabilistic solvers for fully nonlinear parabolic equations existed only at low order (Fahim-Touzi-Warin, Guo-Zhang-Zhuo, Tan), while the multistep machinery of paper 8 had been built only for first-order FBSDEs. Joining the two is the natural next step, but the join brings in something the semilinear case does not have: **the forward SDE is no longer dictated by the equation**. In the Feynman-Kac representation of a semilinear equation, $b$ and $\sigma$ are the equation's coefficients; for a fully nonlinear equation they are free. This paper treats that freedom as a design parameter to be optimised.

### Setting

Cauchy problems for fully nonlinear parabolic equations, converted by the CSTV result into a 2FBSDE:

$$
\begin{cases}
\mathrm dX_t=b(t,X_t)\,\mathrm dt+\sigma(t,X_t)\,\mathrm dW_t,\\
-\mathrm dY_t=f(t,X_t,Y_t,Z_t,\Gamma_t)\,\mathrm dt-Z_t\,\mathrm dW_t,\\
\mathrm dZ_t=A_t\,\mathrm dt+\Gamma_t\,\mathrm dW_t,
\end{cases}
\qquad Y_T=g(X_T).
$$

**This system is written here by inference from equation (2.9) of the sibling paper 19**; the numbering and notation actually printed in this paper are unverified.

### What could be verified

The abstract confirms two things verbatim. First, the method is obtained by extending the results of paper 8 ("by extending our previous results [W. Zhao, Y. Fu and T. Zhou, SIAM J. Sci. Comput., 36 (2014), pp. A1731-A1751]"). Second, and flagged as the specific contribution: "in our numerical schemes, one has the **flexibility to choose the associated forward SDE**, and a suitable choice can significantly reduce the computational complexity."

Where that freedom comes from was explained above: the CSTV PDE does not depend on $\mu,\sigma$, and Theorem 2 of paper 8 says the one-step derivative approximation sees only the left-endpoint coefficients. Stacked together, they permit picking a convenient auxiliary diffusion. **One concrete payoff (this site's reasoning, not the paper's words) is that with constant coefficients the same Gauss-Hermite nodes can be reused at every grid point instead of being recomputed pointwise.**

> [!note] What could be verified
> The publisher returns 403 to direct PDF download, `doc.global-sci.org` is unreachable, and the paper has no preprint, so this site read only the landing page (title, authors, volume/issue/pages, DOI, dates, abstract and the full 28-item reference list) and the OpenAlex record. The four reference ordinary differential equations, the $\alpha_{k,i}$ weights, the Euler step for the forward SDE, local Lagrange interpolation and Gauss-Hermite quadrature can all be read in the sibling paper 19, but **the scheme numbering and notation actually printed in this paper, and whether its spatial discretisation differs from that of paper 19, are unverified here**. No theorem statement was read; by the conventions of the family a $k$-step scheme should be order $k$ in $(Y,Z,\Gamma,A)$, but that is unverified for this paper. For numerics, the abstract says "various numerical examples **including the HJB equations** are presented" — **the "including the HJB equations" part is verifiable** — but the specific test problems and observed orders are not.

### Relation to the others

Papers 16 and 19 are near-twins: **16 is framed PDE-first (fully nonlinear parabolic PDE, then 2FBSDE, then scheme), 19 is framed 2FBSDE-first with the stochastic-control application attached.** The reference lists show this batch was written simultaneously: 16 cites 19 as "arXiv:1502.03206, 2015", 25 as "submitted, 2015" and 23 as "submitted, 2015". Paper 35 later replaces the multistep stencil here with explicit deferred correction for the same 2FBSDE class.

## 19: carrying the multistep machinery to second-order equations, and on to control

### The idea

In the second-order setting the number of objects to discretise rises from two to four: $Y,Z,\Gamma,A$. The two extra ones are not decoration — $\Gamma$ corresponds to the Hessian and is the only doorway for the fully nonlinear term. On the face of it this means four times as many conditional expectations, a prohibitive cost.

But the structure of paper 8 recurs here almost unchanged, and more neatly than one might expect. The key is that **the same "multiply by $\Delta W$ and condition" operation drops one level each time it is applied**: applied to $Y$ it recovers $Z$, applied to $Z$ it recovers $\Gamma$; and differentiating the plain $Z$ reference equation recovers $A$. So all four quantities come from four isomorphic reference ordinary differential equations using the same $\alpha_{k,j}$ weights. Better still, at $t=t_n$ the identity $\Delta W_{t_n,t_n}=0$ kills every term containing $\Delta W$, so $Z^n,\Gamma^n,A^n$ are **all three explicit** and only $Y^n$ is implicit. **The cost does not quadruple, because each extra quantity needs only one more sum.**

### Setting

The coupled second-order forward-backward system (the paper's equation (1.1)) is

$$
\begin{cases}
X_t=x+\displaystyle\int_0^tb(s,\Theta_s)\,\mathrm ds+\int_0^t\sigma(s,\Theta_s)\,\mathrm dW_s,\\[4pt]
Y_t=g(X_T)+\displaystyle\int_t^Tf(s,\Theta_s)\,\mathrm ds-\int_t^TZ_s\,\mathrm dW_s,\\[4pt]
Z_t=Z_0+\displaystyle\int_0^tA_s\,\mathrm ds+\int_0^t\Gamma_s\,\mathrm dW_s,
\end{cases}
\qquad
\Theta_t=(X_t,Y_t,Z_t,A_t,\Gamma_t)\in\mathbb R^m\times\mathbb R\times\mathbb R^d\times\mathbb S^d,
$$

with $b:[0,T]\times\mathbb R^m\times\mathbb R\times\mathbb R^d\times\mathbb S^d\to\mathbb R^m$, $\sigma:\cdots\to\mathbb R^{m\times d}$, $f:\cdots\to\mathbb R$ and $g:\mathbb R^m\to\mathbb R$. **The third equation is everything that is new in the second-order setting**: it expands $Z$ itself as an Itô process, and its diffusion coefficient $\Gamma$ is the process corresponding to the Hessian. Decoupled means $b,\sigma$ independent of $(Y_t,Z_t,A_t,\Gamma_t)$.

The paper's reading of the state of the art: existing high-order FBSDE schemes handle only low dimension, and existing high-dimensional schemes are low order — it cites a 12-dimensional coupled FBSDE example converging at order 1, against the sparse-grid schemes of paper 25 reaching dimension 6 with rates up to 3. Its stated gap: "To the best of our knowledge, there is no related studies for high order numerical methods for 2FBSDEs."

### Derivation

**Step one: the representation theorem (Theorem 2.2).** Let $u(t,x)$ solve the fully nonlinear equation

$$
\mathcal Lu+f\bigl(t,x,u,\nabla_xu\,\sigma,\ \nabla_x(\nabla_xu\,\sigma)\sigma\bigr)=0,\qquad u(T,x)=g(x),
$$

and let $(X_t,Y_t,Z_t,\Gamma_t,A_t)$ solve the decoupled 2FBSDE. Then

$$
Y_t=u(t,X_t),\quad
Z_t=(\nabla_xu\,\sigma)(t,X_t),\quad
\Gamma_t=\bigl(\nabla_x(\nabla_xu\,\sigma)\sigma\bigr)(t,X_t),\quad
A_t=\bigl(\mathcal L(\nabla_xu\,\sigma)\bigr)(t,X_t),
$$

with $\mathcal L$ the generator. This is the fully nonlinear generalisation of the nonlinear Feynman-Kac formula.

**Step two: four reference ordinary differential equations.** Conditioning gives two integral identities,

$$
\mathbb E^x_{t_n}[Y_t]=\mathbb E^x_{t_n}[g(X_T)]+\int_t^T\mathbb E^x_{t_n}[f(s,\Theta_s)]\,\mathrm ds,
\qquad
\mathbb E^x_{t_n}[Z_t]=\mathbb E^x_{t_n}[Z_{t_n}]+\int_{t_n}^t\mathbb E^x_{t_n}[A_s]\,\mathrm ds,
$$

and multiplying each by $\Delta W^{\top}_{t_n,t}$ before conditioning gives two more. Differentiating in $t$:

$$
\frac{\mathrm d\,\mathbb E^x_{t_n}[Y_t]}{\mathrm dt}=-\mathbb E^x_{t_n}[f(t,\Theta_t)],
\qquad
\frac{\mathrm d\,\mathbb E^x_{t_n}[Z_t]}{\mathrm dt}=\mathbb E^x_{t_n}[A_t],
$$

$$
\frac{\mathrm d\,\mathbb E^x_{t_n}\bigl[Y_t\Delta W^{\top}_{t_n,t}\bigr]}{\mathrm dt}
=-\mathbb E^x_{t_n}\bigl[f(t,\Theta_t)\Delta W^{\top}_{t_n,t}\bigr]+\mathbb E^x_{t_n}[Z_t],
\qquad
\frac{\mathrm d\,\mathbb E^x_{t_n}\bigl[Z^{\top}_t\Delta W^{\top}_{t_n,t}\bigr]}{\mathrm dt}
=\mathbb E^x_{t_n}\bigl[A^{\top}_t\Delta W^{\top}_{t_n,t}\bigr]+\mathbb E^x_{t_n}[\Gamma_t].
$$

**The pattern is worth seeing clearly**: $(Y,\Delta W)$ recovers $Z$, $(Z,\Delta W)$ recovers $\Gamma$, and the plain $Z$ equation recovers $A$. At $t=t_n$ every $(\cdot)\Delta W^{\top}_{t_n,t_n}$ term vanishes, so $Z^n,\Gamma^n,A^n$ are explicit and only $Y^n$ is implicit.

**Step three: semi-discrete Scheme 2.** Taking the Euler choice $\bar b(s,\cdot)=b(t_n,x)$, $\bar\sigma(s,\cdot)=\sigma(t_n,x)$ for all $s\in[t_n,T]$:

$$
X^{n,j}=X^n+b(t_n,X^n)\Delta t_{n,j}+\sigma(t_n,X^n)\Delta W_{n,j},\qquad j=1,\dots,k,
$$

$$
Z^n=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}\Delta W^{\top}_{n,j}\bigr],
\qquad
A^n=\sum_{j=0}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Z^{n+j}\bigr],
$$

$$
\Gamma^n=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[(\bar Z^{n+j})^{\top}\Delta W^{\top}_{n,j}\bigr],
\qquad
-\alpha_{k,0}Y^n=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}\bigr]+f(t_n,X^n,Y^n,Z^n,\Gamma^n).
$$

Note that **the $A^n$ sum starts at $j=0$** — it uses the plain $Z$-derivative stencil, which needs the current level — while the other three start at $j=1$.

**Step four: fully discrete Scheme 3.** The same with $x\in\mathcal D^n_h$, interpolation $\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}$ and quadrature $\hat{\mathbb E}^{n,x}$:

$$
X^{n,j}=x+b(t_n,x)\Delta t_{n,j}+\sigma(t_n,x)\Delta W_{n,j},
\qquad
Z^n=\sum_{j=1}^{k}\alpha_{k,j}\hat{\mathbb E}^{n,x}\bigl[\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}Y^{n+j}\Delta W^{\top}_{n,j}\bigr],
$$

$$
\Gamma^n=\sum_{j=1}^{k}\alpha_{k,j}\hat{\mathbb E}^{n,x}\bigl[\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}(Z^{n+j})^{\top}\Delta W^{\top}_{n,j}\bigr],
\qquad
A^n=\sum_{j=0}^{k}\alpha_{k,j}\hat{\mathbb E}^{n,x}\bigl[\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}Z^{n+j}\bigr],
$$

$$
-\alpha_{k,0}Y^n=\sum_{j=1}^{k}\alpha_{k,j}\hat{\mathbb E}^{n,x}\bigl[\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}Y^{n+j}\bigr]+f(t_n,x,Y^n,Z^n,\Gamma^n).
$$

The order of operations at each $(t_n,x)$ is: Euler step, then $Z^n,\Gamma^n,A^n$ explicitly, then $Y^n$ implicitly. **Because grid points are independent of one another, the paper states explicitly that "the process can be completely parallel".** The implicit $Y^n$ can be solved by Picard iteration or, when $f$ is differentiable in $y$, by the Newton iteration

$$
Y^{n,l+1}=Y^{n,l}-\frac{\alpha_{k,0}Y^{n,l}+\sum_{j=1}^k\alpha_{k,j}\hat{\mathbb E}^{n,x}\bigl[\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}Y^{n+j}\bigr]+f(t_n,x,Y^{n,l},Z^n,\Gamma^n)}{\alpha_{k,0}+f_y(t_n,x,Y^{n,l},Z^n,\Gamma^n)} .
$$

**Step five: the coupled case and the spatial discretisation.** Scheme 4 replaces the Euler step by $X^{n,j}=x+b(t_n,x,Y^n,Z^n,\Gamma^n)\Delta t_{n,j}+\sigma(\cdots)\Delta W_{n,j}$, and Scheme 5 wraps it in a Picard iteration initialised by $(Y^{n,0},Z^{n,0},\Gamma^{n,0})=(Y^{n+1},Z^{n+1},\Gamma^{n+1})$; if $b,\sigma$ do not depend on $(Y,Z,\Gamma)$, Scheme 5 reduces to Scheme 3. Spatially: a uniform tensor grid $\mathcal D_h=\{x_0+j\cdot h\}$, local Lagrange interpolation, 10-point Gauss-Hermite quadrature per dimension, and error balancing $h=(\Delta t)^{(k+1)/(r+1)}$ with $r$ the interpolation degree.

### Theorems

- **Local truncation errors (equation (3.22)).** If $\mathbb E^x_{t_n}[Y_t]$, $\mathbb E^x_{t_n}[Y_t\Delta W^{\top}_{t_n,t}]$, $\mathbb E^x_{t_n}[Z_t]$, $\mathbb E^x_{t_n}[Z^{\top}_t\Delta W^{\top}_{t_n,t}]$ and their $t$-derivatives up to order $k+1$ are bounded, then
  $$
  \bar R^k_{y,n}=O((\Delta t)^k),\quad
  \bar R^k_{z,n}=O((\Delta t)^k),\quad
  \bar R^k_{A,n}=O((\Delta t)^k),\quad
  \bar R^k_{\Gamma,n}=O((\Delta t)^k).
  $$
  The fully discrete truncation error splits into **twelve** terms: for each $\bullet\in\{y,z,A,\Gamma\}$ a derivative-approximation term $R^k_{\bullet,n}$, a quadrature term $R^{k,E}_{\bullet,n}$ and an interpolation term $R^{k,I}_{\bullet,n}$.
- **Stability range.** Section 5.1 states that the method "admits a $k$-order convergence rate, and it remains stable for $1\le k\le6$, which is coincide with the classic numerical ODEs theory and our previous results" — the same $k\le6$ barrier as paper 8.

> [!warning] This paper likewise has no convergence theorem
> The arXiv version proves **no** convergence theorem with explicit constants; high order rests on the truncation estimates above plus numerics. Whether the published CiCP version adds one is unverified here.

### Numerical experiments

Common setup: $T=1$, uniform grids, 10 Gauss-Hermite points per dimension, more than 6 Lagrange points, **quadruple precision** (`real(16)`, 34 significant digits), FORTRAN 95 with OpenMP on an Intel Xeon E5-2620 v2. The quadruple precision is worth noticing — paper 8 lost its $k=6$ order to double-precision round-off, and switching to quadruple precision sidesteps exactly that trap.

**Example (5.1), decoupled 2FBSDE.** $\mathrm dX_t=\sin(t+X_t)\mathrm dt+c\cos(t+X_t)\mathrm dW_t$, $-\mathrm dY_t=\bigl(-\cos(t+X_t)\frac1cZ_t-\cos(t+X_t)(Y_t^2+Y_t)-\frac14\Gamma_t\bigr)\mathrm dt-Z_t\mathrm dW_t$, $\mathrm dZ_t=A_t\mathrm dt+\Gamma_t\mathrm dW_t$, $Y_T=\sin(T+X_T)$, $x=0.5$. The exact solution is

$$
Y_t=\sin(t+X_t),\quad
Z_t=c\cos^2(t+X_t),\quad
\Gamma_t=-2c^2\sin(t+X_t)\cos^2(t+X_t),
$$

$$
A_t=-c\sin(2t+2X_t)\bigl(1+\sin(t+X_t)\bigr)-c^3\cos(2t+2X_t)\cos^2(t+X_t).
$$

Order $k$ is reported, stable for $1\le k\le6$. The efficiency comparison is the most persuasive part of this example:

| Scheme | $N$ needed to reach $\lvert Y^0-Y_0\rvert=5.172\times10^{-5}$ | Time    |
| ------ | ------------------------------------------------------------- | ------- |
| $k=1$  | 8192                                                          | 135.0 s |
| $k=2$  | 2048                                                          | 10.82 s |

**Going from first to second order saves more than an order of magnitude in time**, with the same code and the same Euler forward step.

**Example 4 (Table 5.4), convergence rates of all four components.**

| $K$ | $Y$  | $Z$  | $\Gamma$ | $A$  | Time     |
| --- | ---- | ---- | -------- | ---- | -------- |
| 1   | 1.00 | 0.99 | 0.97     | 0.97 | 21.54 s  |
| 2   | 1.98 | 1.95 | 1.99     | 1.97 | 70.10 s  |
| 3   | 2.98 | 2.98 | 3.11     | 2.99 | 259.11 s |

**This table is the paper's most substantial evidence**: not only $Y$ and $Z$ but also the second-order processes $\Gamma$ and $A$ converge at rate $k$. It also shows the cost slope — going from $k=1$ to $k=3$ multiplies the time on the same mesh by roughly 12, consistent with a $k$-step scheme visiting $k$ future levels per step.

**The stochastic-control example of Section 5.2, "tracking a particle under the microscope".** State $\mathrm dX_t=\beta\alpha_t\mathrm dt+\sigma\mathrm dW_t$, cost $J(\alpha)=\mathbb E\bigl[p\int_0^TX_t^2\mathrm dt+q\int_0^T\alpha^2\mathrm dt\bigr]$. The HJB equation is

$$
0=\partial_tV+\inf_{\alpha\in\mathbb R}\Bigl\{\tfrac{\sigma^2}{2}\partial^2_{xx}V+\beta\alpha\partial_xV+px^2+q\alpha^2\Bigr\}
=\partial_tV+\tfrac{\sigma^2}{2}\partial^2_{xx}V-\tfrac{\beta}{4q}(\partial_xV)^2+px^2,
\qquad V(T,x)=0,
$$

with optimal control $\alpha^\ast_t=-\frac{\beta}{2q}\partial_xV(t,x)$. The 2FBSDE actually solved is

$$
\mathrm dX_t=\beta c\,\mathrm dt+\sigma\,\mathrm dW_t,
\qquad
-\mathrm dY_t=\Bigl(-\tfrac{\beta^2}{4q\sigma^2}Z_t^2-\tfrac{\beta c}{\sigma}Z_t+pX_t^2\Bigr)\mathrm dt-Z_t\,\mathrm dW_t,
\qquad
\mathrm dZ_t=A_t\,\mathrm dt+\Gamma_t\,\mathrm dW_t,
$$

with the control recovered by $\alpha^n=-\frac{\beta}{2q\sigma}Z^n$ and parameters $\mu=0.1$, $r=0.03$, $\sigma=0.5$, $c=0.1$. The paper concludes that "the approach is of high order accuracy, both for the 2FBSDE solution and the optimal control $\alpha$"; **this site confirms the claim exists but did not transcribe its data table.**

The paper also states a general recipe: for $\mathrm dX_t=b(t,X_t,\alpha_t)\mathrm dt+\sigma(t,X_t,\alpha_t)\mathrm dW_t$ with $J(\alpha)=\mathbb E[\int_0^Tf(t,X_t,\alpha_t)\mathrm dt+g(X_t)]$, the HJB is $\partial_tV+\inf_{\alpha\in U}\{\frac{\sigma^2}{2}\partial^2_{xx}V+b\,\partial_xV-f\}=0$ with $\alpha^\ast(t,x)=\arg\inf_\alpha\{\cdots\}$; substituting gives $\partial_tV+G(t,x,\partial_xV,\partial^2_{xx}V)=0$ with $G(t,x,p,P)=\frac{\sigma(t,x,\alpha^\ast)^2}{2}P+b(t,x,\alpha^\ast)p-f(t,x,\alpha^\ast)$; the corresponding 2FBSDE is solved and the control recovered as $\alpha^\ast_t=g(X_t,Y_t,Z_t,\Gamma_t)$. **The essential premise of the recipe is that $\inf_\alpha$ can be computed analytically** — precisely what the deep-learning route of papers 86 and 96 later sets out to avoid.

**What these experiments establish, and where they fall short.** They establish that all four components attain order $k$ in the second-order setting, which is the paper's central claim, and Table 5.4 supports it directly; the efficiency comparison also shows that high order pays off in measured time. Three shortfalls. First, the examples remain low-dimensional and smooth with tensor Gauss-Hermite quadrature; the paper's own introduction frames the gap as "high order only in low dimension, high dimension only at low order", and this paper closes the first half only. Second, quadruple-precision arithmetic is a substantive implementation choice, which leaves open whether this family can reach $k=6$ stably in double precision (paper 8's $k=6$ in double precision was polluted by round-off). Third, the high order of the control example is only asserted in prose; this site did not transcribe its data table, so no order is reported here for it.

### Relation to the others

A direct extension of paper 8 (its abstract: "we extend (with non-trivial updates) our multistep schemes"). A near-twin of paper 16 (same authors, same year, same 2FBSDE machinery, differing only in PDE-first versus 2FBSDE-first framing). Paper 35 later replaces the multistep stencil here with explicit deferred correction. **Its route from stochastic control through the HJB equation to a 2FBSDE is the dynamic-programming alternative to the Pontryagin/adjoint route of papers 26 and 41**; papers 86 and 96 later attack the same HJB problem in high dimension with neural networks, and notably avoid the explicit $\inf_\alpha$ computed analytically here. Its introduction cites paper 25 for reaching dimension 6 with rates up to 3, the benchmark used to frame the low-order/high-dimension trade-off.

## 25: pushing the cost of dimension from exponential to polynomial

### The idea

The time discretisation of paper 8 has no dimension problem: the $\alpha_{k,j}$ weights do not depend on dimension. What does have a dimension problem is the **spatial** half — storing functions on a uniform tensor grid, computing conditional expectations with tensor Gauss-Hermite quadrature, evaluating with local Lagrange interpolation. All three costs grow exponentially with dimension. The paper's diagnosis is blunt: this combination "is less efficient for high dimensional FBSDEs, as the required computational work increases exponentially as the dimension increases."

So this paper **leaves the time discretisation entirely alone** and replaces only the spatial half, with Smolyak sparse grids. The sparse-grid idea is to take not all nodes of a tensor product but only those sub-tensor-products whose multi-index sum is bounded; for functions with mixed smoothness very little accuracy is lost, while the node count drops from exponential to polynomial-times-logarithmic. The paper adds one more layer: **hierarchical bases** turn the sparse interpolant into a single spectral expansion whose coefficients can be recovered by a dimension-by-dimension fast transform.

### Setting

The time discretisation is verbatim that of paper 8: the same two reference ordinary differential equations, the same $\alpha_{k,i}$ Vandermonde system, the same local-generator theorem, and the same iterative semi-discrete algorithm for coupled systems,

$$
X^{n,j}=X^n+b(t_n,X^n,Y^{n,l},Z^{n,l})\Delta t_{n,j}+\sigma(t_n,X^n,Y^{n,l},Z^{n,l})\Delta W_{n,j},\quad j=1,\dots,k,
$$

$$
Z^{n,l+1}(X^n)=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}(\Delta W_{n,j})^{\top}\bigr],
\qquad
\alpha_{k,0}Y^{n,l+1}(X^n)=-\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}\bigr]-f(t_n,X^n,Y^{n,l+1},Z^{n,l+1}),
$$

iterated until $\max\{|Y^{n,l+1}-Y^{n,l}|,|Z^{n,l+1}-Z^{n,l}|\}<\epsilon_0$. The paper restates the known property: "the $k$-th step scheme admits a $k$ order convergence rate, provided that $1\le k\le6$." It describes itself as "the second part in a series of papers on multi-step schemes for solving coupled FBSDEs".

### Derivation

**Step one: sparse grids.** From a one-dimensional sequence $\chi_i=\{x^i_0,\dots,x^i_{N_i-1}\}$ (typically $N_i=2^i+1$), the $q$-dimensional Smolyak sparse grid is

$$
\chi^p_q=\bigcup_{q\le|\mathbf i|_1\le p}\chi_{i_1}\otimes\chi_{i_2}\otimes\cdots\otimes\chi_{i_q},\qquad p\ge q .
$$

For a **nested** sequence $\chi_1\subset\chi_2\subset\cdots$ this collapses to $\chi^p_q=\bigcup_{q\le|\mathbf i|_1\le p}\tilde\chi_{i_1}\otimes\cdots\otimes\tilde\chi_{i_q}$ with $\tilde\chi_1=\chi_1$ and $\tilde\chi_i=\chi_i\setminus\chi_{i-1}$ for $i>1$. Two one-dimensional families are used, with a clear division of labour:

- **Chebyshev-Gauss-Lobatto (CGL)**, $C_i=\{x^i_j=\cos(j\pi/2^i),\ j=0,\dots,2^i\}$ — **nested**, used for function approximation and interpolation;
- **Gauss-Hermite (GH)**, $G_i=\{x^i_j,\ j=1,\dots,2^i-1\}$ with $x^i_j$ the roots of the Hermite polynomial of order $2^i-1$ — explicitly noted as **not nested**, used for the quadrature that evaluates conditional expectations.

**This split is not arbitrary: interpolation needs nesting to be written as hierarchical increments, quadrature does not, so quadrature can keep the non-nested but more accurate Gauss-Hermite nodes.**

**Step two: sparse interpolation.** With $\Delta_1=\mathcal I_1$ and $\Delta_i=\mathcal I_i-\mathcal I_{i-1}$ for $i>1$, the Smolyak interpolation operator is

$$
\mathcal I^p_q[f]=\sum_{q\le|\mathbf i|_1\le p}\Delta_{i_1}\otimes\cdots\otimes\Delta_{i_q}[f]
=\sum_{p-q<|\mathbf i|_1\le p}(-1)^{p-|\mathbf i|_1}\binom{q-1}{p-|\mathbf i|_1}\ \mathcal I_{i_1}\otimes\cdots\otimes\mathcal I_{i_q}[f].
$$

**Step three: where "spectral" enters — hierarchical bases plus a fast transform.** A basis $\{\tilde\phi_k\}$ is called **hierarchical** (Definition 3.1) if $\tilde\phi_k(x_j)=0$ for all $j\in I_i$ and $k\notin I_i$. With nested grids and hierarchical bases the expansion coefficients become **independent of the level index**, so the interpolant collapses to a single spectral expansion

$$
\mathcal I^p_q[f](x)=\sum_{\mathbf k\in I^p_q}b_{\mathbf k}\,\tilde\phi_{\mathbf k}(x),
\qquad
\tilde\phi_{\mathbf k}=\prod_{i=1}^{q}\tilde\phi_{k_i},
$$

whose coefficients solve $f(x_{\mathbf j})=\sum_{\mathbf k\in I^p_q}b_{\mathbf k}\tilde\phi_{\mathbf k}(x_{\mathbf j})$ for all $\mathbf j\in I^p_q$. That system is solved by a **dimension-by-dimension fast transform** (Algorithm 2, `FastTran`), applying the one-dimensional inverse matrix $T=(\tilde\phi_k(x_j))^{-1}_{k,j}$ along one dimension at a time over the hierarchical blocks — the FFT-type acceleration advertised in the abstract. The basis used is a transformed Chebyshev basis $T_k(x)$.

**Step four: sparse-grid Gauss-Hermite quadrature.** For the conditional expectations, the tensor GH rule is replaced by the Smolyak combination $\sum\Delta_{k_1}\otimes\cdots\otimes\Delta_{k_q}$ of one-dimensional GH quadrature operators $Q_k$.

**Step five: the full algorithm.** Inputs are $Y^{N-i}(x),Z^{N-i}(x)$ on $C^{p_i}_q[a_i,b_i]$ for $i=0,\dots,k-1$; each time level runs `FastTran` to obtain spectral coefficients $\{\beta^{N-i}_{\mathbf j}\},\{\gamma^{N-i}_{\mathbf j}\}$, then evaluates the multistep formulas with sparse-grid GH quadrature. The paper contrasts the two approaches in a table:

| Method           | Meshes                 | Conditional expectations     | Approximation and interpolation |
| ---------------- | ---------------------- | ---------------------------- | ------------------------------- |
| SSG (this paper) | sparse grid            | sparse-grid GH quadrature    | sparse-grid interpolation       |
| LTG (= paper 8)  | tensor-product uniform | tensor-product GH quadrature | Lagrangian                      |

### Theorems

**This site found no self-contained convergence theorem in the paper.** It transports the order-$k$ claim from paper 8 and quotes standard sparse-grid quadrature error results from the literature (its words: it "obtain[s] the following result of the sparse grid quadrature for functions defined on a high-dimensional cube", but **the exact constants were not transcribed here**). Its concrete contribution is at the complexity level: running time grows "in certain polynomial level (non-exponential)" with dimension — **a numerical observation (Figure 2), not a theorem.**

### Numerical experiments

FORTRAN 95 on 16 Intel Xeon E5620 CPUs (2.40 GHz) with 3.0 GB RAM.

**Example 1 (two-dimensional coupled, periodic).** $b_i=\cos(4(x_i+t))/4-1$, $\sigma$ diagonal with $\sigma_{ii}=\cos(4(x_i+t))\sin(4(x_i+t))/4$, and $f$ chosen so that the exact solution is

$$
Y_t=\sin\bigl(4(X_{t,1}+t)\bigr)\sin\bigl(4(X_{t,2}+t)\bigr),
\qquad
Z_{t,i}=\Bigl(\prod_{k=1}^{2}\sin\bigl(4(X_{t,k}+t)\bigr)\Bigr)\cos^2\bigl(4(t+X_{t,i})\bigr).
$$

Periodicity lets the problem be posed on $[-\pi,\pi]^2$. SSG uses the CGL sparse grid $C^7_2$ for space and the GH sparse grid $G^3_2$ for quadrature; LTG uses $\Delta x=(\Delta t)^{(k+1)/(n+1)}$ to balance errors. Rates over $N=8,\dots,128$:

| Scheme       | $E_Y$ rate | $E_Z$ rate |
| ------------ | ---------- | ---------- |
| SSG, 1 step  | 0.982      | 0.986      |
| SSG, 2 steps | 1.987      | 1.955      |
| SSG, 3 steps | 2.632      | 2.955      |
| LTG, 1 step  | 1.174      | —          |

Running times show SSG substantially cheaper at comparable accuracy.

**Example 2 ($q$-dimensional decoupled).** $b_i=\frac1qx_ie^{-x_i^2}$, $\sigma_{ii}=\frac1qe^{-x_i^2}$, with exact solution $Y_t=\frac1q\sum_{j=1}^q\bigl(X_{t,j}^2\prod_{k\ne j}(X_{t,k}+t)\bigr)$ and a matching $Z_{t,i}$. Solved for $q=3,4,5,6$: "the proposed multi-step schemes admit high order convergence rates even for the 6-dimensional problem", and running time "grows in certain polynomial level (non-exponential)" in $q$.

**Example 3 ($q$-dimensional coupled).** $b_i=\frac t2\cos^2(y+x_i)$, $\sigma_{ii}=\frac t2\sin^2(y+x_i)$.

**What these experiments establish, and where they fall short.** They establish that sparse grids carry the same time discretisation to six dimensions while keeping high order, and that the time curve is no longer exponential — exactly what the paper set out to show. Three shortfalls. First, **six dimensions is the measured ceiling of this route, not its floor**; the later deep-learning papers (86, 93, 96, 97, 100) reach thousands of dimensions by abandoning grids entirely. Second, the three-step scheme in Example 1 shows only 2.632 for $E_Y$, clearly below the nominal 3, against 2.955 for $E_Z$ — **the two components are out of step in the error table, the paper does not explain it, and this page will not explain it on the paper's behalf.** Third, Example 1 relies on periodicity to compress the domain into a bounded box, which is what makes the sparse grid usable; how to truncate a general unbounded problem is not covered by the experiments.

### Relation to the others

Self-described as "the second part in a series" after paper 8, keeping its time discretisation verbatim and replacing only the spatial machinery — **the clearest single demonstration that the temporal and spatial lines of improvement in this thread are entirely independent.** Paper 19 cites it for reaching dimension 6 with rates up to 3, using it as the benchmark for the "high order and low dimension versus low order and high dimension" trade-off. Paper 16's reference list carries it as "submitted, 2015". Paper 63 takes a different route to the same problem for BSDEs, using Sinc quadrature and interpolation, which has the advantage of **removing the need for spatial interpolation altogether**. Papers 86, 93, 96, 97 and 100 abandon grids completely and reach thousands of dimensions; **paper 25 marks the practical ceiling of the deterministic-grid approach.**

## 26: organising the control iteration as gradient projection

### The idea

A constrained optimal control problem is hard in two places: the constraint and the gradient. The constraint has a clean treatment — write "minimise over the convex set $K$" as the fixed-point equation $u^\ast=P_K(u^\ast-\rho J'(u^\ast))$, so the algorithm becomes "compute a gradient, take a step, project back". That step concentrates every remaining difficulty in the gradient.

And the gradient is hard because computing $J'$ directly requires the Gâteaux derivative $Dx^u_t(v)$ of the state process, that is, one variational SDE per direction $v$ — unaffordable. The adjoint method exists to kill that direction dependence **once and for all**: introduce an adjoint process $(p,q)$, apply Itô's formula to the product $p^u_tDx^u_t(v)$, and every term containing $Dx^u_t(v)$ cancels, leaving an expression in which $v$ appears only explicitly and linearly. One adjoint solve then gives the gradient in all directions.

**There is a technical dividend here that is easy to skip past: the generator of the adjoint BSDE is linear in $(p,q)$.** So although its left-point rectangle scheme is nominally implicit, one rearrangement solves it, with no iteration at all — in contrast to general BSDE schemes, which must run Picard iterations.

### Setting

The constrained stochastic optimal control problem is

$$
\min_{u\in K}J(u)=\mathbb E\Bigl[\int_0^T\bigl(h(x_t^u)+j(u(t))\bigr)\mathrm dt+k(x_T^u)\Bigr],
\qquad
\mathrm dx_t^u=b(x_t^u,u(t))\,\mathrm dt+\sigma(x_t^u,u(t))\,\mathrm dW_t,\quad x|_{t=0}=x_0 .
$$

The paper sorts existing numerical routes into four classes: reduction to finite-dimensional stochastic programming; dynamic programming, that is, solving the HJB equation ("one of the most widely used numerical methods"); martingale-based methods; and methods based on the **stochastic maximum principle**. Its stated gap is specific: while the stochastic maximum principle is a popular tool for theoretical studies, "it has not been widely used in the numerical setting". That is the gap this paper fills.

An easily missed part of the setting: the control space $U=L^2([0,T];\mathbb R)$ consists of **deterministic** square-integrable controls, which the paper defends as what future planning in engineering and financial applications needs; the adapted (feedback) case $U_{\mathbb F}=L^2_{\mathbb F}([0,T]\times\Omega;\mathbb R)$ is handled separately in Section 5. The exposition is one-dimensional, with the statement that "the whole framework applies easily to multi-dimensional cases".

### Derivation

**Step one: projection gives a fixed-point characterisation.** The first-order optimality condition is the variational inequality $(J'(u^\ast),v-u^\ast)\ge0$ for all $v\in K$, with $(J'(u),v)=\lim_{\rho\downarrow0}\frac{J(u+\rho v)-J(u)}{\rho}$. Define the projection $P_K\omega=\arg\min_{u\in K}\|u-\omega\|$, equivalently characterised by $(P_K\omega-\omega,v-P_K\omega)\ge0$ for all $v\in K$. Substituting $\omega=u^\ast-\rho J'(u^\ast)$ and comparing with the variational inequality gives, for any $\rho>0$, the **fixed-point characterisation**

$$
u^\ast=P_K\bigl(u^\ast-\rho J'(u^\ast)\bigr).
$$

**Step two: discretised control space and the iteration.** With $0=t_0<\cdots<t_N=T$, $\Delta t=T/N$ and $I^N_n=[t_{n-1},t_n)$, set $U_N=\{u\in U\mid u=\sum_{n=1}^N\alpha_n\chi_{I^N_n}\ \text{a.e.},\ \alpha_n\in\mathbb R\}$ and $K_N=K\cap U_N$ (convex and closed). Then $u^{\ast,N}=P_{K_N}(u^{\ast,N}-\rho J'(u^{\ast,N}))$, and the algorithm is the **gradient projection iteration**

$$
u^{i+1,N}=P_{K_N}\bigl(u^{i,N}-\rho_iJ'_N(u^{i,N})\bigr),\qquad i=1,2,\dots,
$$

with $J'_N$ a numerical approximation of $J'$.

**Step three: the gradient via an adjoint BSDE.** Written out directly,

$$
(J'(u),v)=\mathbb E\Bigl[\int_0^Th'(x^u_t)Dx^u_t(v)\,\mathrm dt+\int_0^Tj'(u(t))v(t)\,\mathrm dt+k'(x^u_T)Dx^u_T(v)\Bigr],
\qquad
Dx^u_t(v)=\lim_{\rho\downarrow0}\frac{x^{u+\rho v}_t-x^u_t}{\rho},
$$

where $Dx^u_t(v)$ solves the variational SDE $\mathrm dDx^u_t(v)=(b'_xDx^u_t(v)+b'_uv(t))\mathrm dt+(\sigma'_xDx^u_t(v)+\sigma'_uv(t))\mathrm dW_t$. To eliminate $Dx^u_t(v)$, introduce the **adjoint BSDE**

$$
-\mathrm dp_t^u=f(x_t^u,p_t^u,q_t^u,u(t))\,\mathrm dt-q_t^u\,\mathrm dW_t,
\qquad p_T^u=g(x_T^u)=k'(x_T^u),
$$

$$
f(x,p,q,u)=h'(x)+p\,b_x'(x,u)+q\,\sigma_x'(x,u).
$$

Applying Itô's formula to $p_t^uDx_t^u(v)$ cancels the $Dx^u_t(v)$ terms and leaves the **gradient representation**

$$
J'(u)\big|_t=\mathbb E\bigl[p_t^u\,b_u'(x_t^u,u(t))+q_t^u\,\sigma_u'(x_t^u,u(t))\bigr]+j'(u(t)).
$$

**Note that this generator is linear in $p,q$** — which is why the scheme below, though nominally implicit, needs no iteration.

> [!warning] A substantive difference from an earlier route
> The paper's Remark 1 draws the contrast deliberately: the adjoint equation of its reference [12] "is an anticipating integrand stochastic differential equation, where the solution is required to be backward-adapted instead of the classic forward adapted. However, such a requirement is not true in general" — that is, its well-posedness is unclear, whereas the BSDE above is well posed by standard theory. **This is a methodological difference of substance, not a matter of technical preference.**

**Step four: the FBSDE and its Feynman-Kac form.** Together the system and its nonlinear Feynman-Kac representation are

$$
\begin{cases}
\mathrm dx^u_t=b(x^u_t,u(t))\,\mathrm dt+\sigma(x^u_t,u(t))\,\mathrm dW_t,\quad x|_{t=0}=x_0,\\
-\mathrm dp^u_t=f(x^u_t,p^u_t,q^u_t,u(t))\,\mathrm dt-q^u_t\,\mathrm dW_t,\quad p^u_T=g(x^u_T),
\end{cases}
$$

$$
p_t=\eta(t,x_t),\qquad q_t=\sigma(x_t,u(t))\,\partial_x\eta(t,x_t),
$$

with $\eta$ solving $\mathcal L^0\eta(t,x)=-f(x,\eta,\sigma\partial_x\eta,u(t))$, $\eta(T,x)=g(x)$, and $\mathcal L^0\eta=\partial_t\eta+b\,\partial_x\eta+\frac12\sigma^2\partial_{xx}\eta$.

**Step five: the Euler scheme.** Integrating the backward equation over $[t_n,t_{n+1}]$, taking $\mathbb E^x_{t_n}[\cdot]$ and applying the **left-point rectangular rule** to the generator:

$$
p^x_{t_n}=\mathbb E^x_{t_n}[p_{t_{n+1}}]+\Delta t\,f(x,p^x_{t_n},q^x_{t_n},u(t_n))+\bar R^x_{p,n},
\qquad
q^x_{t_n}=\frac{1}{\Delta t}\Bigl(\mathbb E^x_{t_n}[p_{t_{n+1}}\Delta W_{n+1}]+\bar R^x_{q,n}\Bigr),
$$

and dropping the truncation terms gives the semi-discrete scheme

$$
p^x_n=\mathbb E^x_{t_n}[p_{n+1}]+\Delta t\,f(x,p^x_n,q^x_n,u(t_n)),
\qquad
q^x_n=\frac{1}{\Delta t}\mathbb E^x_{t_n}[p_{n+1}\Delta W_{n+1}].
$$

In the language of the $\theta$-scheme family this is the member $\theta_1=\theta_2=1$: the generator is evaluated at the **current** level $t_n$ with the unknowns $(p^x_n,q^x_n)$, and $q^x_n$ comes from the single term $\mathbb E_n[p_{n+1}\Delta W_{n+1}]/\Delta t$. **It is formally implicit, but because the adjoint generator is linear in $p$, one rearrangement gives $p^x_n$ in closed form with no iteration.**

**Step six: the conditional expectations.** The Euler state $\tilde x^{t_n,x}_{t_{n+1}}=x+b(x,u(t_n))\Delta t+\sigma(x,u(t_n))\Delta W_{n+1}$ turns the conditional expectations into Gaussian integrals against $\rho(\xi)=\frac{1}{\sqrt{2\pi}}e^{-\xi^2/2}$, approximated by the $L$-point Gauss-Hermite rule:

$$
\hat{\mathbb E}^x_{t_n}[p_{t_{n+1}}]=\sum_{\ell=1}^{L}\mathcal I_hp_{t_{n+1}}\bigl(x+b(x,u(t_n))\Delta t+\sigma(x,u(t_n))\sqrt{\Delta t}\,\xi_\ell\bigr)\,\omega_\ell,
$$

$$
\hat{\mathbb E}^x_{t_n}[p_{t_{n+1}}\Delta W_{n+1}]=\sum_{\ell=1}^{L}\mathcal I_hp_{t_{n+1}}\bigl(x+b(x,u(t_n))\Delta t+\sigma(x,u(t_n))\sqrt{\Delta t}\,\xi_\ell\bigr)\sqrt{\Delta t}\,\xi_\ell\,\omega_\ell,
$$

with $\mathcal I_h$ **linear interpolation** onto the uniform spatial grid $\mathcal R_h=\{x_k\}$, finite in practice with $|k|\le P$. The fully discrete scheme substitutes these, solving for $(p^k_n,q^k_n)$ for $n=N-1,\dots,0$ at every $x_k$, and then $J'_N(u)|_{t_n}=\hat{\mathbb E}[p_nb'_u+q_n\sigma'_u]+j'(u(t_n))$. Algorithm 1 loops: pick $u^0\in U_N$ and a tolerance $\epsilon_0$; (i) set the terminal condition $p^k_N=g(x_k)$; (ii) march backwards for $(p_n,q_n)$; (iii) compute $J'_N$; (iv) update $u$ by the projection step; repeat.

**Step seven: feedback controls (Algorithm 2).** When the control may be a feedback $\bar u(t,x)$ and the constraint set $K$ is **pointwise** in both time and space, the expectation disappears from the gradient:

$$
J'(u)^k_n=p^k_n\,b'_u(x_k,\bar u(t_n,x_k))+q^k_n\,\sigma'_u(x_k,\bar u(t_n,x_k))+j'(\bar u(t_n,x_k)),
$$

and the projection degenerates to a pointwise scalar projection $\bar u^\ast(t_n,x)=P_C(\cdots)$. **The key consequence**, as the paper puts it: "we no longer need the history information before time $t$ to compute $J'(u)_t$, but only the information at time instance $t$", so the whole algorithm runs in a **single backward sweep** with the control update nested inside the time loop, cutting storage sharply.

### Theorems

- **Theorem 1 (convergence of the iteration).** Assume $J'$ is Lipschitz and uniformly monotone near $u^\ast$ and $u^{\ast,N}$: there exist $c,C>0$ such that for all $v\in K$
  $$
  \|J'(u^\ast)-J'(v)\|\le C\|u^\ast-v\|,
  \qquad
  \bigl(J'(u^\ast)-J'(v),\,u^\ast-v\bigr)\ge c\|u^\ast-v\|^2,
  $$
  and the same two inequalities hold with $(u^{\ast,N},K_N)$ in place of $(u^\ast,K)$. Assume also $\epsilon_N=\sup_i\|J'(u^{i,N})-J'_N(u^{i,N})\|\to0$ as $N\to\infty$. If the step sizes $\rho_i$ satisfy
  $$
  0<1-2c\rho_i+(1+2C)\rho_i^2\le\delta^2\quad\text{for some }0<\delta<1,
  $$
  then $\|u^\ast-u^{i,N}\|\to0$ as $i,N\to\infty$.
- **Corollary 1.** Under Theorem 1's conditions, if additionally $u^\ast$ and $J'(u^\ast)$ are Lipschitz in $U$, then $\epsilon_N\sim O(\Delta t)$ implies $\|u^\ast-u^{i,N}\|\sim O(\Delta t)$ as $i\to\infty$.
- **Theorem 2 (the main error estimate).** Under Assumption 1 and the assumptions of Lemmas 1 through 3, writing $(\mu_n,\nu_n)$ for the errors in $(p^k_n,q^k_n)$,
  $$
  \hat{\mathbb E}\bigl[(\mu_n)^2\bigr]+\Delta t\sum_{n=0}^{N-1}\hat{\mathbb E}\bigl[(\nu_n)^2\bigr]
  =O\bigl((\Delta t)^2\bigr)+O\bigl((\Delta x)^4/(\Delta t)^2\bigr),
  $$
  whence
  $$
  \epsilon_N=\sup_i\|J'(u^{N,i})-J'_N(u^{N,i})\|=O(\Delta t)+O\bigl((\Delta x)^2/\Delta t\bigr),
  $$
  and **in particular, taking $\Delta x=\Delta t$ gives $\epsilon_N=O(\Delta t)$**, so by Corollary 1, $\|u^\ast-u^{N,i}\|=O(\Delta t)$ as $i\to\infty$. That is, **the entire algorithm — control iteration, BSDE solver, quadrature and interpolation together — is first order.** Supporting results: Proposition 2 requires $m\ge2$, $L\ge2$ and $\Delta x=O(\sqrt{\Delta t})$ for the discrete conditional-expectation operator; Lemma 1 assumes $b,\sigma\in C^{0,4}_b$; the local truncation errors satisfy $\tilde R^k_{q,n}=O((\Delta t)^2)$.

> [!warning] $\Delta x$ and $\Delta t$ cannot be refined independently
> The form of the second term, $O((\Delta x)^4/(\Delta t)^2)$, is worth remembering: **the spatial interpolation error is divided by $(\Delta t)^2$, so refining the spatial mesh at fixed $\Delta t$ does not keep helping, and refining the time step at fixed $\Delta x$ makes this term blow up.** Balancing the two forces the spatial mesh to be refined together with the time step, the characteristic cost of interpolation-based conditional-expectation approximations. The $\Delta x=\Delta t$ balancing condition is the low-order special case of the $h=(\Delta t)^{(k+1)/(r+1)}$ rule in papers 8, 19 and 25.

### Numerical experiments

The step sizes are $\rho_i=1/\sqrt i$, with the paper's explanation that a small $\rho$ helps the estimate converge but too small a $\rho$ slows the iteration.

**Example 1.** $J(u)=\frac12\int_0^T\mathbb E[(x_t-x^\ast(t))^2]\mathrm dt+\frac12\int_0^Tu^2(t)\mathrm dt$, $K=U$, state $\mathrm dx_t=u(t)x_t\mathrm dt+\sigma x_t\mathrm dW_t$. Setup $x_0=1$, $T=1$, $\sigma=0.1$, $M=10^5$ samples, tolerance $\epsilon_0=10^{-5}$; first-order convergence is observed. **This site transcribed only part of the exact control expression for this example, so it is not reproduced here.** A second test with $\sigma=0.1$, $M=10^5$, $\epsilon_0=10^{-5}$ and $N=40,50,\dots,100$ concludes that "the numerical solution matches the exact solution very well and first order convergence rate is observed".

**Example 3 — feedback control.** The same state and cost, but with $K$ a set of **stochastic** controls. Comparing Algorithm 1 (deterministic control) with Algorithm 2 (feedback control), the values of $J(u)$ make the paper's clearest single table:

| $N$ | Algorithm 1 (deterministic) | Algorithm 2 (feedback) |
| --- | --------------------------- | ---------------------- |
| 100 | 0.84833                     | 0.62535                |
| 200 | 0.84797                     | 0.64507                |
| 400 | 0.84777                     | 0.65509                |
| 800 | 0.84770                     | 0.66013                |

The paper concludes: "the use of feedback control can indeed improve the results (produces a smaller value of objective functional), and this is reasonable as we are minimizing the objective functional within a larger control set." **Note in passing that the two columns move in opposite directions**: Algorithm 1 decreases monotonically towards about $0.8477$ while Algorithm 2 increases monotonically towards about $0.660$ — both converging, but approaching their respective limits from opposite sides.

**Example 4 — portfolio problem.** The reference optimal value is $J(u)=15023$; with $N=1000,2000,4000,8000$ and $M=N^2/10$, "the method admits a first order rate of convergence".

**What these experiments establish, and where they fall short.** They establish that the first order predicted by Theorem 2 is observed in all four examples, and the feedback comparison adds a conclusion independent of order, about the size of the control set. Three shortfalls. First, first order here is the design target rather than evidence of a limitation — the paper deliberately uses only the Euler scheme with matching first-order analysis; whether higher order is attainable is the subject of paper 41. Second, this site transcribed Example 1's exact solution and Example 2's content only partially, so no error tables for those two appear on this page. Third, every example is one-dimensional, while the paper only asserts that the framework extends easily to several dimensions, without multi-dimensional experiments.

### Relation to the others

**This is the Pontryagin, or stochastic-maximum-principle, counterpart of the dynamic-programming route taken by paper 19.** Here the FBSDE is the adjoint/Hamiltonian system, not a Feynman-Kac representation of the value function. It deliberately uses only a first-order Euler scheme with matching first-order analysis; paper 41 is the direct sequel that upgrades exactly this component to the high-order multistep schemes of paper 8 while keeping the gradient-projection outer loop. Its statement that it "follow[s] closely" earlier work for the Euler scheme points to the Zhao-school BSDE-scheme literature that papers 8 and 47 belong to. Paper 50 also takes the adjoint/SMP route, but in a partially observed, data-driven feedback setting.

## 41: replacing the inner solver with a high-order one

### The idea

The overall accuracy of paper 26 is first order, and the bottleneck is clear: the adjoint BSDE uses an Euler scheme and the outer loop uses fixed-step gradient projection. Both can be replaced — the inner one by a high-order member of the family of paper 8, the outer one by a quasi-Newton update that uses curvature information.

**The headline claim inherits the central phenomenon of paper 8**: even with the state equation still discretised by Euler, the whole thing reaches second order. The reason is the same as in paper 8 — what is actually approximated in a control problem is the adjoint quantities and the gradient, and those see the state process only through conditional expectations.

### Setting

The same stochastic optimal control problem as paper 26. This paper converts it into "an equivalent stochastic optimality system of FBSDEs": the forward state SDE, the backward adjoint BSDE, and the variational inequality or stationarity condition in the control, combined into one coupled system.

### What could be verified

The abstract gives three points verbatim. First, the stochastic optimal control problem is converted into "an equivalent stochastic optimality system of FBSDEs". Second, the authors design "an efficient **second order FBSDE solver**" and "an quasi-Newton type optimization solver" for the resulting system. Third, "it is noticed that our approach admits the **second order rate of convergence even when the state equation is approximated by the Euler scheme**". The keywords are "FBSDEs; stochastic optimal control; stochastic maximum principle; projected quasi-Newton methods", with AMS classes 60H35, 93E20, 93E25, 49M29, 65C20 and 65K15.

The composition of the reference list is verifiable evidence too: it cites papers 8, 19, 25, 26 and 23, as well as Zhao-Zhang-Ju (SINUM 2010), Zhao-Chen-Peng (SISC 2006), Peng's stochastic maximum principle, Pontryagin et al., Yong and Zhou, Powell-Yuan trust regions, Dai-Yuan conjugate gradients and He's projection-and-contraction methods. **That list says directly where each of its two components comes from: the high-order FBSDE solver from the multistep and deferred-correction toolbox, the optimiser from the nonlinear-programming literature.** It also cites Yang and Zhao's "Convergence of recent multistep schemes for a forward-backward stochastic differential equation" (EAJAM 2015), the convergence theory that paper 47 later generalises.

> [!note] What could be verified
> `doc.global-sci.org` is unreachable, `global-sci.org` returns 403 to direct download, and the paper has no preprint, so this site read only the landing page (with the full 46-item reference list) and the OpenAlex record. Specifically, **whether the second-order FBSDE solver is the $k=2$ member of the paper 8 / paper 47 family or an equivalent $\theta$ or Crank-Nicolson-type scheme is unverified here**; so are the exact quasi-Newton update (BFGS? limited memory?) and the globalisation strategy (line search? trust region?). Whether the second-order claim is proved as a theorem or demonstrated numerically, and under what hypotheses and constants, is likewise unverified — the venue and reference list make a proof plausible, but this site read none. For numerics, the abstract says "several numerical examples are presented to illustrate the effectiveness and the accuracy of the proposed numerical schemes", but the test problems and observed orders are unverified.

### Relation to the others

**The direct sequel to paper 26**, which it cites: the same problem and the same SMP/adjoint-BSDE route, but second order instead of first and quasi-Newton instead of fixed-step gradient projection. It also cites papers 8, 19, 25 and 23 — the whole multistep and deferred-correction toolbox, which is where the high-order FBSDE solver comes from. **The division of labour between papers 26 and 41 is a recurring pattern in this thread: one paper improves the outer iteration, the other the inner solve.** It contrasts with the HJB/dynamic-programming route of paper 19 and with the neural-network route of papers 86 and 96; notably, its reference list already includes Han and E's "Deep learning approximation for stochastic control problems" and Pereira et al.'s "Learning deep stochastic optimal control policies using forward-backward SDEs", showing that the authors were aware of the deep-learning alternative in 2019.

## 50: once observations enter, filtering and control cannot be separated

### The idea

Every paper above assumes the state is fully observed: knowing $X_t$, the control follows from the optimality condition. Once only a noisy indirect observation $M_t$ is available, the control must be measurable with respect to the observation filtration, and the structure of the problem changes. The classical remedy is the **separation principle** — filter the state, then treat the estimate as if it were the truth — but it holds rigorously only for time-invariant linear systems with linear observations; the general nonlinear problem remains open.

Even settling for gluing filtering and control together is prohibitively expensive: the maximum-principle route requires solving the adjoint FBSDE over the whole state space at every gradient step, plus a Zakai-equation filter, which does not run beyond one dimension.

**This paper's move is to abandon the goal of solving the FBSDE accurately at all.** It observes that the FBSDE is here only a **vehicle for the gradient**, and stochastic gradient descent only ever needed unbiased noisy gradients. So each gradient step draws **one** trajectory: one particle from the cloud as the starting point, one sequence of Gaussian increments, one backward sweep. Two layers of expectation — over particles and over Monte-Carlo paths — collapse into one draw. **This eliminates the spatial representation of $(Y,Z)$ as functions of $X$ entirely**, which is why the method can run in real-time feedback, and also why it no longer claims to solve the FBSDE accurately.

### Setting

State, cost and observation are

$$
\mathrm dX_t=b(t,X_t,u_t)\,\mathrm dt+\sigma(t,X_t,u_t)\,\mathrm dW_t,\quad X_0=\xi,
\qquad
J(u)=\mathbb E\Bigl[\int_0^Tf(t,X_t,u_t)\,\mathrm dt+h(X_T)\Bigr],
$$

$$
\mathrm dM_t=g(X_t)\,\mathrm dt+\mathrm dB_t,\qquad M_0=0,
$$

with $X$ valued in $\mathbb R^d$, $u$ valued in $U\subseteq\mathbb R^m$, $g:\mathbb R^d\to\mathbb R^\ell$, and $B$ an $\ell$-dimensional Brownian motion independent of $W$. **The control enters both the drift and the diffusion.** Writing $\mathbb F^M$ for the augmented filtration generated by $M$, the admissible set is $\mathcal U_{\rm ad}[0,T]=\{u:\ \mathbb F^M\text{-progressively measurable, }U\text{-valued}\}$ and Problem (C\*) is $J^\ast(u^\ast)=\inf_{u^M\in\mathcal U_{\rm ad}[0,T]}J^\ast(u^M)$. The authors call this a **data driven feedback control** problem. They also note that under the Girsanov change of measure $\mathrm d\mathbb P^M=\Theta^T_t\,\mathrm d\mathbb P$ with $\Theta^T_t=\exp(-\int_t^Tg(X_s)\mathrm dB_s-\int_t^T\frac12|g(X_s)|^2\mathrm ds)$, the observation $M$ becomes a standard Brownian motion.

### Derivation

**Step one: the gradient from the stochastic maximum principle.** When $u^\ast$ is interior to $\mathcal U_{\rm ad}$,

$$
(J^\ast)'_u(u^\ast_t)=\mathbb E\Bigl[b_u(t,X^\ast_t,u^\ast_t)^\top Y_t+\sigma_u(t,X^\ast_t,u^\ast_t)^\top Z_t+f_u(t,X^\ast_t,u^\ast_t)^\top\ \Big|\ \mathcal F^M_t\Bigr],
$$

where $(Y,Z,\zeta)$ solves the **adjoint (Pontryagin) FBSDE system**

$$
\begin{cases}
\mathrm dX^\ast_t=b(t,X^\ast_t,u^\ast_t)\,\mathrm dt+\sigma(t,X^\ast_t,u^\ast_t)\,\mathrm dW_t, & X_0=\xi,\\[2pt]
\mathrm dM^\ast_t=g(X^\ast_t)\,\mathrm dt+\mathrm dB_t, & M_0=0,\\[2pt]
\mathrm dY_t=\bigl(-b_x^\top Y_t-\sigma_x^\top Z_t-f_x^\top\bigr)\mathrm dt+Z_t\,\mathrm dW_t+\zeta_t\,\mathrm dB_t, & Y_T=h_x(X^\ast_T)^\top .
\end{cases}
$$

Here $Z$ is the martingale-representation integrand of $Y$ against $W$ and $\zeta$ that against $B$. **Note that the backward equation carries two integrands, because there are two sources of noise.**

**Step two: conditional gradient descent.** The plain iteration is $u^{l+1,M}_t=u^{l,M}_t-\rho\,(J^\ast)'_u(u^{l,M}_t)$. But at the current time $t$ the future observations $\{\mathcal F^M_s\}_{s>t}$ are unavailable, so the control is projected onto $\mathcal F^M_t$: writing $u^{l,M}_s|_t:=\mathbb E[u^{l+1,M}_s|\mathcal F^M_t]$,

$$
u^{l+1,M}_s\big|_t=u^{l,M}_s\big|_t-\rho\,\mathbb E\bigl[(J^\ast)'_u\bigl(u^{l,M}_s|_t\bigr)\ \big|\ \mathcal F^M_t\bigr],
\qquad l=0,1,\dots,\ s\in[t,T].
$$

The observation equation is **deliberately excluded** from the driving FBSDE, since $M$ does not appear in the gradient expression and future data are unavailable.

**Step three: time discretisation.** With $\mathbb E_i[\cdot]=\mathbb E[\cdot\mid\mathcal F^{X,B}_{t_i}]$, $\Delta t_i=t_{i+1}-t_i$ and $\Delta W_{t_i}=W_{t_{i+1}}-W_{t_i}$, for $i=N_T-1,\dots,n$:

$$
X_{i+1}=X_i+b(t_i,X_i,u^{l,M}_{t_i}|_{t_n})\Delta t_i+\sigma(t_i,X_i)\Delta W_{t_i},
$$

$$
Y_i=\mathbb E_i[Y_{i+1}]+\mathbb E_i\Bigl[b_x^\top Y_{i+1}+\sigma_x^\top Z_{i+1}+f_x^\top\Bigr]\Delta t_i,
\qquad
Z_i=\mathbb E_i\bigl[Y_{i+1}\Delta W_{t_i}\bigr]\cdot(\Delta t_i)^{-1},
$$

with $b_x,\sigma_x,f_x$ evaluated at level $t_{i+1}$. That is: the forward SDE uses the **left-point (Euler-Maruyama)** rule; $Y$ uses the **right-point** rule on the drift integral after conditioning, so the backward step is **explicit** ($Y_{i+1},Z_{i+1}$ are already known); and $Z$ uses the left-point rule with the classical $\mathbb E_i[Y_{i+1}\Delta W_{t_i}]/\Delta t_i$ representation. **Remark 3.1 notes that $\zeta$ appears neither in the scheme nor in $(J^\ast)'_u$, so no scheme for it is needed at all.**

**Step four: conditional expectations by Monte Carlo rather than quadrature.** Unlike papers 8, 19, 25 and 26, this one uses plain Monte Carlo with $K$ standard-Gaussian samples $\{\omega^k_i\}$ and $\Delta W_{t_i}\approx\sqrt{\Delta t_i}\,\omega^k_i$:

$$
Y_i=\frac1K\sum_{k=1}^KY^k_{i+1}+\frac{\Delta t_i}{K}\sum_{k=1}^K\Bigl[b_x(t_{i+1},X^k_{i+1},\cdot)^\top Y^k_{i+1}+\sigma_x(t_{i+1},X^k_{i+1})^\top Z^k_{i+1}+f_x(t_{i+1},X^k_{i+1},\cdot)^\top\Bigr],
$$

$$
Z_{t_i}=\frac{1}{\Delta t_i}\sum_{k=1}^K\frac{Y^k_{i+1}\sqrt{\Delta t_i}\,\omega^k_i}{K}.
$$

The authors' stated reason is effectiveness in high dimension.

**Step five: the particle filter and the single-realisation trick.** The conditional distribution $p(X_t|\mathcal F^M_t)$ is represented by a particle cloud $\{x^{(s)}_n\}_{s=1}^S$ propagated by a particle filter. The paper's central computational idea is to collapse **two** layers of expectation (the $S$ particles and the $\Lambda$ Monte-Carlo paths) into **one** randomly drawn realisation per gradient step, that is, into stochastic gradient descent:

$$
Y^{(\hat l,\hat s)}_i=Y^{(\hat l,\hat s)}_{i+1}+\Bigl[b_x^\top Y^{(\hat l,\hat s)}_{i+1}+\sigma_x^\top Z^{(\hat l,\hat s)}_{i+1}+f_x^\top\Bigr]\Delta t_i,
\qquad
Z^{(\hat l,\hat s)}_i=Y^{(\hat l,\hat s)}_{i+1}\,\omega^{(\hat l,\hat s)}_i\,(\Delta t_i)^{-1/2},
$$

with the starting point $X^{(\hat l,\hat s)}_{t_n}=x^{(\hat s)}_n$ drawn from the particle cloud and $\omega^{(\hat l,\hat s)}_i\sim N(0,1)$.

> [!warning] The paper itself concedes that this does not solve the FBSDE accurately
> The authors state explicitly that the scheme **does not solve the FBSDE accurately**; the FBSDE is only a vehicle for the gradient process, and what they appeal to is the usual stochastic-gradient-descent justification of unbiased noisy gradients. **This is exactly where the spatial approximation of $(Y,Z)$ as functions of $X$ is eliminated, and it is where the paper parts ways with paper 26 in objective — 26 buys accuracy and theory, 50 buys speed and dimension.**

The overall algorithm (PF-SGD) is structured as follows: initialise the particle cloud and the iteration count $L$; for $n=0,1,\dots,N_T$, initialise the estimated control process and the step size, and for $l=0,\dots,L$ simulate **one** trajectory, sweep backwards for $\{Y^{(l,s)}_i\}$, and update the control; take $\hat u^\ast(t_n)=u^{L,M}_{t_n}|_{t_n}$; propagate the particles with the particle filter using it.

### Theorems

**The preprint states no convergence theorem for PF-SGD.** It cites external analyses instead (numerical analysis of the backward scheme and its extensions in prior Bao-Cao-Webster-type and Zhao-school work) along with the standard stochastic-gradient-descent convergence literature. **This is therefore an algorithmic and computational paper, not an error-analysis one.** The one clearly stated theoretical ingredient is the gradient representation derived in the appendix from the stochastic maximum principle with Gâteaux derivatives.

### Numerical experiments

**Example 1 — an LQ benchmark with nonlinear observations.** $\mathrm dX_t=A(t)X_t\mathrm dt+BU_t\mathrm dt+C\,\mathrm dW_t$, cost $J(U)=\mathbb E\bigl[\frac12\int_0^T(\langle QX_t,X_t\rangle+\langle RU_t,U_t\rangle)\mathrm dt+\frac12\langle FX_T,X_T\rangle\bigr]$, whose **fully observed** optimal control is the Riccati feedback $\bar U_t=-R^{-1}B^\top P(t)X_t$ with $\dot P=-PA-A^\top P+PBR^{-1}B^\top P-Q$, $P(T)=F$. The observations are **nonlinear**: $M_t=\sin(X_t)+\eta_t$ with Gaussian noise of covariance $\Gamma$. The purpose is to check that the data-driven control tracks the analytical one. **No convergence order is reported.**

**Example 2 — a one-dimensional nonlinear, non-quadratic control problem, used for the efficiency study.** $\mathrm dX_t=\arctan(X_t+u_t)\mathrm dt+\sigma X_t\mathrm dW_t$ with $\sigma=0.05$; cost $J^\ast(u^M)=\mathbb E\bigl[\frac12\int_0^T\sin^2(X_t+u^M_t)\mathrm dt\bigr]$; observation $M_t=X_t+\eta_t$ with noise standard deviation $0.05$. The benchmark is a "full solution method" (Zakai-equation filter plus a fully solved FBSDE with grid interpolation) on three meshes. $T=1$; PF-SGD uses $\Delta t=0.02$ ($N_T=50$), 500 particles and 1000 SGD iterations; the full-solution meshes are $\Delta t=0.1/\Delta x=0.1$, $\Delta t=0.05/\Delta x=\frac{\sqrt2}{2}\cdot0.1$ and $\Delta t=0.025/\Delta x=0.05$ over $[3,6]$.

| Method                | Coarse mesh | Finer mesh | Finest mesh | PF-SGD      |
| --------------------- | ----------- | ---------- | ----------- | ----------- |
| Overall cost $J^\ast$ | 0.0481      | 0.0318     | 0.0076      | **0.00095** |
| CPU time (s)          | 29.78       | 220.47     | 1560.15     | **0.93**    |

**This table is the strongest evidence in the paper**: PF-SGD reaches a cost an order of magnitude below the finest-mesh full solution in about $1/1700$ of the time. A repeated-experiment average cost $\hat J^\ast_t(\hat u)=\frac{1}{M_{\rm rept}}\sum_m\frac12\int_0^t\sin^2(\hat X^{(m)}_s+\hat u^{(m)}_s)\mathrm ds$ over $M_{\rm rept}=50$ runs confirms the same ordering over time.

**Example 3 — a Dubins vehicle with bearings-only observations.** $\mathrm dX_t=\sin\theta_t\mathrm dt+\sigma\mathrm dW_t$, $\mathrm dY_t=\cos\theta_t\mathrm dt+\sigma\mathrm dW_t$, $\mathrm d\theta_t=u_t\mathrm dt+\sigma_2\mathrm dW_t$ with $\sigma=0.2$; cost $J^\ast(u^M)=\mathbb E\bigl[\int_0^T\frac12(u^M_t)^2\mathrm dt+\delta((X_T-X_P)^2+(Y_T-Y_P)^2)\bigr]$ with $\delta=10$ and target $(X_P,Y_P)=(5,3)$. Two bearing detectors at $(6,1)$ and $(-1,4)$ give $M_t=[\arctan\frac{X_t-6}{Y_t-1},\ \arctan\frac{X_t+1}{Y_t-4}]^\top+\eta_t$ with noise standard deviation $0.01I_2$. $T=1$, $\Delta t=0.02$ ($N_T=50$), 1000 particles, 1000 SGD steps. Qualitative trajectory plots only.

**What these experiments establish, and where they fall short.** They establish that on partially observed problems with nonlinear observations PF-SGD produces a usable feedback control, and that in Example 2, where a comparison against the full solution method is available, it wins on both cost and time. Three shortfalls. First, **no convergence order is reported and no convergence theorem is proved**, so the lower cost is an observation on this set of problems, not a general guarantee. Second, Example 3 offers only trajectory plots with no quantitative metric, so its three-dimensional behaviour cannot be connected to Example 2's conclusion. Third, Example 2 compares against one particular implementation of the full solution method, which uses grid interpolation — **the table says that on this problem abandoning the spatial representation beats refining it, not that PF-SGD is more accurate in general.**

### Relation to the others

**The partially observed counterpart of paper 26.** Both start from the same adjoint-FBSDE gradient representation; 26 does full-observation gradient projection with Gauss-Hermite quadrature and interpolation plus a proved convergence rate, while 50 replaces every conditional expectation by a one-sample stochastic gradient and adds a particle filter for the observation filtration. In its treatment of conditional expectations it contrasts sharply with papers 8, 19, 25 and 41: those use quadrature plus interpolation on a (sparse) grid, while 50 deliberately avoids any spatial representation of $(Y,Z)$. Its backward step is the $k=1$ member of the multistep family of paper 8, and the paper **does not** pursue high-order time accuracy — its objective is real-time feedback. Feng Bao and Richard Archibald are outside the Shandong/CAS multistep line; **this paper is where Tao Zhou's FBSDE-based control machinery meets the nonlinear-filtering and data-assimilation community.** It also points forward to papers 86 and 96, which likewise abandon grid-based spatial representations but replace Monte-Carlo and particle representations with neural networks.

## 51: block iteration under nonsmooth constraints

### The idea

This paper is not about FBSDEs but about dynamic nonlinear complementarity problems, which couple an ordinary differential equation with a complementarity condition. After backward-Euler discretisation, every time point requires solving a coupled nonlinear system whose two blocks are entirely different in character: one smooth differential system and one nonsmooth complementarity system. **Existing approaches either merge the two blocks (direct elimination, semismooth Newton) or come with no convergence guarantee; the idea here is not to merge them but to alternate.**

The Gauss-Seidel-style alternation is the most naive one imaginable: solve the complementarity system with the old $x$ to get a new $y$, then solve the differential system with the new $y$ to get a new $x$. Each subsystem can then be handed to an off-the-shelf solver. The paper's real contribution is the half that had always been missing: **why it converges, and how fast.** And the answer has an attractive structure — for a fixed time interval the convergence is **superlinear** with a rate **independent of the step size $h$**; for a fixed number of time points it is linear with rate $O(h)$.

### Setting

The dynamic nonlinear complementarity problem is

$$
\dot x(t)=F(t,x(t),y(t)),
\qquad
0\le y(t)\ \perp\ G(t,x(t),y(t))\ge0,
\qquad t\in(0,T),\ x(0)=x_0,
$$

with $x(t)\in\mathbb R^m$ and $y(t)\in\mathbb R^n_+$. Two subclasses are the **differential semiaffine system** $0\le y\perp Nx+My+g(t)\ge0$ and the **dynamic linear complementarity problem** $\dot x=Ax+By+f(t)$, $0\le y\perp Nx+My+g\ge0$. Backward Euler gives

$$
0\le y_j\ \perp\ G(t_j,x_j,y_j)\ge0,
\qquad
x_j=x_{j-1}+hF(t_j,x_j,y_j),
\qquad j=1,\dots,N_t,\ h=T/N_t .
$$

The paper identifies definite defects in both mainstream approaches. **Direct elimination** works only in the linear case: after substituting out $x_j$, the reduced linear complementarity system has matrix $M_h:=hN(I-hA)^{-1}B+M$, which **may fail to be a P-matrix even when $M$ is**, and forming it requires $n$ large linear solves. **Semismooth Newton** is only locally convergent and requires the Clarke generalised Jacobian, which is expensive and awkward for large problems.

### Derivation

**Step one: the single-point Gauss-Seidel iteration.**

$$
0\le y^{k+1}_j\ \perp\ G\bigl(t_j,\,x^{k}_j,\,y^{k+1}_j\bigr)\ge0,
\qquad
x^{k+1}_j=x_{j-1}+hF\bigl(t_j,\,x^{k+1}_j,\,y^{k+1}_j\bigr).
$$

The complementarity system uses the **old** $x^k_j$ and the differential system uses the **new** $y^{k+1}_j$ — that is what makes it Gauss-Seidel-style. The two subsystems are decoupled, so any existing NCP solver (linearisation-based methods, the PATH solver) and any Newton solver for the smooth $F$-system plug in unchanged.

**Step two: the multipoint (window) version, the parallel one.** Split $\{t_1,\dots,t_{N_t}\}$ into $P$ groups of $J=N_t/P$ points and, within one group,

$$
0\le y^{k+1}_j\ \perp\ G\bigl(t_j,x^k_j,y^{k+1}_j\bigr)\ge0,
\qquad
x^{k+1}_j=x^{k+1}_{j-1}+hF\bigl(t_j,x^{k+1}_j,y^{k+1}_j\bigr),
\qquad j=1,\dots,J .
$$

**The point is that the $\{y^{k+1}_j\}_{j=1}^J$ are mutually independent**, so all $J$ complementarity systems can be solved **in parallel**; the $x$-recursion is then sequential in $j$. The paper spells out several special cases: the semiaffine version; the dynamic-linear-complementarity version, where each iteration costs exactly one linear solve $(I-hA)x^{k+1}_j=x_{j-1}+hBy^{k+1}_j+hf_j$; and, for $M$ a Z-matrix, $y^{k+1}_j$ obtainable from the linear program $\min\|y\|_1$ subject to $y\ge0$ and $My+Nx^k_j+g_j\ge0$.

**Step three: the hypotheses.** $G$ is a **uniform P-function of $y$**,

$$
\max_{1\le l\le n}(\bar y_l-\tilde y_l)\bigl(G_l(t,x,\bar y)-G_l(t,x,\tilde y)\bigr)\ \ge\ L_0\|\bar y-\tilde y\|_2^2,\qquad L_0>0,
$$

together with $\|G(t,x,\bar y)-G(t,x,\tilde y)\|_2\le L_G\|\bar y-\tilde y\|_2$; and $F$ satisfies a **one-sided Lipschitz condition in $x$** and an ordinary Lipschitz condition in $y$:

$$
\langle F(t,\bar x,y)-F(t,\tilde x,y),\bar x-\tilde x\rangle\le L_1\|\bar x-\tilde x\|_2^2,
\qquad
\|F(t,x,\bar y)-F(t,x,\tilde y)\|_2\le L_2\|\bar y-\tilde y\|_2,
$$

with $L_1\in(-\infty,\infty)$, which **may be negative**, and $L_2>0$. Lemma 2.1 gives a unique $\mathcal Y(x)$ solving the static NCP, Lipschitz in $x$ with constant $\eta_t$.

### Theorems

- **Lemma 2.2 (a combinatorial identity).** With $\psi(r,J,k)=\sum_{j_1=1}^{J}\sum_{j_2=1}^{j_1}\cdots\sum_{j_k=1}^{j_{k-1}}r^{J-j_k}$,
  $$
  \psi(r,J,k)=
  \begin{cases}
  \dfrac{1}{(1-r)^k}-\displaystyle\sum_{l=1}^{k}\binom{J+k-l-1}{k-l}\dfrac{r^{J}}{(1-r)^l}, & r\neq1,\\[10pt]
  \dbinom{J+k-1}{k}, & r=1 .
  \end{cases}
  $$
  **This identity is the technical pivot of every convergence result**: the error recursion of the multipoint iteration is a $k$-fold nested sum, and $\psi$ is what evaluates it in closed form.
- **Theorem 2.3 (basic error bound).** Under the hypotheses above and $hL_1<1$, the error $e^k_j=x_j-x^k_j$ of the multipoint method satisfies
  $$
  \max_{0\le j\le J}\|e^k_j\|_2\le
  \begin{cases}
  (h\tilde\eta)^k\,\psi(1,J,k)\,\max_{0\le j\le J}\|e^0_j\|_2, & L_1=0,\\[6pt]
  \psi\bigl((1-hL_1)^{-1},J,k\bigr)\Bigl(\dfrac{h\tilde\eta}{1-hL_1}\Bigr)^{k}\max_{1\le j\le J}\|e^0_j\|_2, & L_1\neq0,
  \end{cases}
  $$
  with $\tilde\eta=L_2\eta$ and $x_j$ the converged (exact backward-Euler) solution.
- **Theorem 2.4 (fixed number of time points gives linear rate $O(h)$).** Under the assumptions of Theorem 2.3, for $k\gg1$ with $J$ fixed and $h(L_1+\tilde\eta)<1$,
  $$
  \max_{0\le j\le J}\|e^k_j\|_2\le\rho^k\max_{0\le j\le J}\|e^0_j\|_2,
  \qquad
  \rho=\frac{h\tilde\eta}{1-hL_1}=O(h),\quad\tilde\eta=\eta L_2 .
  $$
- **Theorem 2.5 (fixed time interval gives superlinear, $h$-independent convergence).**
  $$
  \|x^k(t)-x(t)\|_2\le\max\{1,e^{L_1t}\}\,\frac{(t\tilde\eta)^k}{k!}\,\sup_{t\in[0,T]}\|x^0(t)-x(t)\|_2,
  \qquad t\in(0,T).
  $$
  **The $k!$ in the denominator makes the convergence superlinear at a rate independent of $h$** (the paper's Remark 2.2).
- **Remark 2.3.** A **negative** $L_1$ accelerates convergence; for dynamic linear complementarity problems this happens when $A$ is orthogonally similar to its Jordan form with eigenvalues of negative real part — **the typical situation when the ODE comes from semi-discretising a parabolic PDE**, for instance a parabolic Signorini problem.
- **The $P_0$ case.** If $G$ is only a $P_0$-function ($L_0=0$), the results still apply after the Tikhonov regularisation $\widehat G(t,x,y)=G(t,x,y)+\varepsilon y$, with $\varepsilon$ chosen comparable to the temporal error $O(h)$.
- **The linear complementarity case (Section 3, $G(t,x,y)=My(t)+\widetilde G(t,x(t))$).** Theorem 3.2 gives **unique existence of a least-norm solution** $(x,y)\in C^1(0,t)\times C(0,t)$ on $(0,T^\ast)$ for positive semidefinite $M$, with the explicit horizon
  $$
  T^\ast=
  \begin{cases}
  T, & L_1\le-\dfrac{\tilde\eta\beta+C_0}{\beta},\\[8pt]
  \min\Bigl\{T,\ \dfrac1{L_1}\log\Bigl(1+\dfrac{L_1\beta}{\tilde\eta\beta+C_0}\Bigr)\Bigr\}, & L_1>-\dfrac{\tilde\eta\beta+C_0}{\beta},
  \end{cases}
  $$
  where $\tilde\eta=L_2\eta_0\eta_1$, $C_0=\max_{t\in[0,T]}\|F(x_0,\mathcal Y(\widetilde G(t,x_0)))\|_2$ and $\mathcal B(x_0,\beta)=\{v:\|v-x_0\|_2\le\beta\}$, assuming $\mathrm{FEA}(M,\widetilde G(t,v))\neq\emptyset$ on that ball. Theorem 3.3 shows the least-norm iterative method is **well defined** (feasibility is preserved along the iteration) under $hL_1<1$ and an explicit restriction on the number of usable time points $J^\ast$. The case of $M$ a Z-matrix is treated separately in Remark 3.1.

### Numerical experiments

The stopping rule is $\max_{0\le j\le J}\|x^k_j-x_j\|_2\le10^{-8}$ against the converged solution $\{x_j\}$. The two examples are a **four-diode bridge wave rectifier** (a nonsmooth circuit with a nonlinear resistor and a capacitor of **random** value — this is where Tao Zhou's uncertainty-quantification interest enters) and a **projected dynamic system arising from the spatial price equilibrium problem**.

The reported conclusion is that the method is superior to existing methods in robustness, complexity and computation time, and that both theoretical regimes — $h$-robust superlinear convergence at fixed $T$ and $O(h)$ linear convergence at fixed $J$ — are observed. **This site did not transcribe the per-experiment iteration counts and timings, so no quantitative comparison appears on this page.**

### Relation to the others

**It is essentially disjoint from the FBSDE line**: no BSDE, no conditional expectation, no probabilistic scheme. Its only structural kinship is the **splitting and decoupling** philosophy — solving two coupled subsystems alternately, exactly as the coupled FBSDE solvers of papers 8, 23 and 26 iterate between the forward SDE and the backward equation. The gradient-projection outer loop of paper 26 is likewise a fixed-point iteration on a coupled system with a projection onto a convex admissible set, making it **the closest mathematical relative in this list**. Shu-Lin Wu is a parallel-in-time (waveform relaxation and parareal) specialist, and the multipoint window version is a waveform-relaxation-type idea, which is why parallelism over the $J$ time points matters. The randomised circuit parameters of Section 4.1 connect to Tao Zhou's uncertainty-quantification work rather than to the FBSDE work.

## How the seven fit together

| No. | What it treats                                                | Position relative to the others                      | Verification here           |
| --- | ------------------------------------------------------------- | ---------------------------------------------------- | --------------------------- |
| 16  | fully nonlinear parabolic PDEs                                | start of the second-order setting, PDE-first framing | abstract and reference list |
| 19  | second-order FBSDEs and stochastic control                    | attaches the second-order setting to control         | full text, equations        |
| 25  | cost of evaluating multi-dimensional conditional expectations | improvement on the spatial side                      | full text, equations        |
| 26  | constrained stochastic optimal control                        | outer iteration (gradient projection), first order   | full text, equations        |
| 41  | high-accuracy schemes for stochastic control                  | inner solve (high-order FBSDE), second order         | abstract and reference list |
| 50  | data-driven feedback control                                  | partial observation, filtering coupled to control    | full text, equations        |
| 51  | dynamic nonlinear complementarity problems                    | block iteration under nonsmooth constraints          | full text, equations        |

The page has two shapes. The first is **two attacks on the same control problem**: paper 19 goes through dynamic programming (HJB, then 2FBSDE, with $\inf_\alpha$ solved analytically), while papers 26 and 41 go through Pontryagin (adjoint BSDE, then gradient, then projection or quasi-Newton). The former needs $\inf_\alpha$ to be available in closed form, the latter needs a computable gradient; the difficulties sit in completely different places. The second shape is **where the cost is saved**: paper 25 saves it inside the grid (sparse grids compress exponential into polynomial), while paper 50 discards the grid entirely (one-sample stochastic gradients) and gives up accuracy guarantees in exchange.

## Coverage check

| Item                                                                            | Paper      | Status                                                                          |
| ------------------------------------------------------------------------------- | ---------- | ------------------------------------------------------------------------------- |
| CSTV second-order setting, the role of $\Gamma$, two notation differences       | background | complete                                                                        |
| The design freedom to choose the forward SDE                                    | 16         | only what the abstract supports; equations inferred from paper 19               |
| Representation theorem, four reference ODEs, five schemes, Newton iteration     | 19         | complete derivation                                                             |
| Truncation errors, the twelve-term split, the $k\le6$ range                     | 19         | complete, with the missing convergence theorem flagged                          |
| Four-component rate table and the efficiency comparison                         | 19         | Table 5.4 and Example 5.1 complete; the control example's table not transcribed |
| Sparse grids, hierarchical bases, fast transform, sparse quadrature             | 25         | complete derivation                                                             |
| Two-dimensional and $q$-dimensional experiments, polynomial growth              | 25         | Example 1 complete; Examples 2 and 3 setup and conclusions only                 |
| Fixed-point characterisation, adjoint BSDE, left-point rectangle scheme         | 26         | complete derivation, including the Remark 1 criticism                           |
| Theorem 1, Corollary 1, Theorem 2 and error balancing                           | 26         | complete, with hypotheses                                                       |
| Feedback comparison table and first-order convergence                           | 26         | Example 3 complete; Examples 1, 2 and 4 partly transcribed                      |
| Second-order solver and quasi-Newton optimiser                                  | 41         | only what the abstract and reference list support                               |
| Partially observed setting, adjoint system, the PF-SGD single-realisation trick | 50         | complete derivation, including the paper's own caveat                           |
| Efficiency comparison table and three examples                                  | 50         | Example 2 complete; Examples 1 and 3 setup and qualitative conclusions only     |
| Gauss-Seidel iteration, parallel multipoint version, four theorems              | 51         | complete                                                                        |
| Two numerical examples                                                          | 51         | setup and conclusions only; data tables not transcribed                         |

## Sources for this page

- T. Kong, W. Zhao, and T. Zhou, [_Probabilistic high order numerical schemes for fully nonlinear parabolic PDEs_](https://doi.org/10.4208/cicp.240515.280815a), Commun. Comput. Phys. 18(5) (2015), pp. 1482-1503.
- T. Kong, W. Zhao, and T. Zhou, [_High order numerical schemes for second-order FBSDEs with applications to stochastic optimal control_](https://doi.org/10.4208/cicp.OA-2016-0056), Commun. Comput. Phys. 21(3) (2017), pp. 808-834 (preprint [arXiv:1502.03206](https://arxiv.org/abs/1502.03206)).
- Y. Fu, W. Zhao, and T. Zhou, [_Efficient spectral sparse grid approximations for solving multi-dimensional forward backward SDEs_](https://doi.org/10.3934/dcdsb.2017174), Discrete Contin. Dyn. Syst. Ser. B 22(9) (2017), pp. 3439-3458 (preprint [arXiv:1607.06897](https://arxiv.org/abs/1607.06897)).
- B. Gong, W. Liu, T. Tang, W. Zhao, and T. Zhou, [_An efficient gradient projection method for stochastic optimal control problems_](https://doi.org/10.1137/17M1123559), SIAM J. Numer. Anal. 55(6) (2017), pp. 2982-3005.
- Y. Fu, W. Zhao, and T. Zhou, [_Highly accurate numerical schemes for stochastic optimal control via FBSDEs_](https://doi.org/10.4208/nmtma.OA-2019-0137), Numer. Math. Theor. Meth. Appl. 13(2) (2020), pp. 296-319.
- R. Archibald, F. Bao, J. Yong, and T. Zhou, [_An efficient numerical algorithm for solving data driven feedback control problems_](https://doi.org/10.1007/s10915-020-01358-y), J. Sci. Comput. 85(2) (2020), 58 (preprint [arXiv:2006.03047](https://arxiv.org/abs/2006.03047)).
- S. Wu, T. Zhou, and X. Chen, [_A Gauss-Seidel type method for dynamic nonlinear complementarity problems_](https://doi.org/10.1137/19M1268884), SIAM J. Control Optim. 58(6) (2020), pp. 3389-3412.
- External source used for cross-checking: P. Cheridito, H. M. Soner, N. Touzi, and N. Victoir, [_Second-order backward stochastic differential equations and fully nonlinear parabolic PDEs_](https://doi.org/10.1002/cpa.20168), Comm. Pure Appl. Math. 60(7) (2007), pp. 1081-1110 (the source from which this site transcribed the definition of the second-order setting, the PDE correspondence, and the original statement that the forward diffusion may be chosen freely).
