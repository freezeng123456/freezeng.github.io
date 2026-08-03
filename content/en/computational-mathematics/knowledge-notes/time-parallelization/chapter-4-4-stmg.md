---
title: "4.6: Space–Time Multigrid (STMG)"
description: Complete derivation from the all-at-once system and temporal block-Jacobi smoothing to local Fourier analysis and nonlinear FAS
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-4-stmg
tags:
  - parallel-in-time
  - space-time multigrid
  - local Fourier analysis
---

> [!note] Reading scope
> This page follows Section 4.6 (pp. 472–481). It covers equations (4.30)–(4.44), Theorem 4.9, Figures 4.18–4.22, and Table 4.1. The linear two-level cycle, smoothing symbol, damping choice, time-integrator dependence, and nonlinear FAS are all retained.

## 4.6 Space–time multigrid (STMG)

Space–time multigrid goes back to the parabolic multigrid of Hackbusch (1984) and its improvement by Horton and Vandewalle (1995). The STMG of this section is due to Gander and Neumüller (2016), whose central finding is that replacing the temporal smoother by block Jacobi makes the method as effective as multigrid applied to Poisson problems, using only standard multigrid components. Related work includes Janssen and Vandewalle (1996), Van Lent and Vandewalle (2002), and Chaudet-Dumas, Gander and Pogozelskyte (2024).

### All-at-once system (4.30)–(4.31)

As in the ParaDiag methods of Section 3.5, all time unknowns are first collected into one system. Spatial discretization of heat or advection–diffusion gives

$$
\boldsymbol u'=A\boldsymbol u+\boldsymbol f.
$$

A generic one-step method is

$$
r_1\boldsymbol u_{n+1}
=r_2\boldsymbol u_n+\widetilde{\boldsymbol f}_n,
\qquad n=0,\ldots,N_t-1, \tag{4.30}
$$

where $r_1,r_2$ are matrix polynomials in $\Delta tA$. Backward Euler has $r_1=I_x-\Delta tA,r_2=I_x$; the trapezoidal rule has $r_1=I_x-\frac12\Delta tA,r_2=I_x+\frac12\Delta tA$.

Stacking the time points gives

$$
\underbrace{
\begin{bmatrix}
r_1\\-r_2&r_1\\&\ddots&\ddots\\&&-r_2&r_1
\end{bmatrix}}_{K}
\underbrace{
\begin{bmatrix}
\boldsymbol u_1\\\boldsymbol u_2\\\vdots\\\boldsymbol u_{N_t}
\end{bmatrix}}_{\boldsymbol U}
=\boldsymbol b. \tag{4.31}
$$

STMG builds coarse grids in both space and time to solve this coupled system.

### Time-parallel block-Jacobi smoother

The damped block-Jacobi smoother is

$$
\boldsymbol U^{\mathrm{new}}
=S_\eta(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s):
\left\{
\begin{aligned}
\boldsymbol U^0&=\boldsymbol U^{\mathrm{ini}},\\
(I_t\otimes r_1)\Delta\boldsymbol U^j
&=\eta(\boldsymbol b-K\boldsymbol U^j),
&&j=0,\ldots,s-1,\\
\boldsymbol U^{j+1}&=\boldsymbol U^j+\Delta\boldsymbol U^j,\\
\boldsymbol U^{\mathrm{new}}&=\boldsymbol U^s.
\end{aligned}
\right. \tag{4.32}
$$

$I_t\otimes r_1$ is block diagonal, so all $N_t$ spatial corrections in one sweep are independent. Eta is the damping parameter and $s$ the sweep count.

For $N_x=7$, linear interpolation and restriction are

$$
P_x=
\begin{bmatrix}
1/2&0&0\\1&0&0\\1/2&1/2&0\\0&1&0\\0&1/2&1/2\\0&0&1\\0&0&1/2
\end{bmatrix},
\qquad
R_x=\frac12P_x^\top. \tag{4.33}
$$

The time operators $P_t,R_t$ are analogous.

### Two-level cycle (4.34)

Let `Mat` reshape an all-at-once vector into a space-by-time matrix and `Vec` reverse it. One two-level iteration is

$$
\left\{
\begin{aligned}
\boldsymbol U^{k+1/3}&=S_\eta(\boldsymbol b,\boldsymbol U^k,s_1),\\
\boldsymbol r&=\boldsymbol b-K\boldsymbol U^{k+1/3},\\
\boldsymbol r_c&=[R_x\operatorname{Mat}(\boldsymbol r)]R_t^\top,\\
\boldsymbol e_c&=K_c^{-1}\operatorname{Vec}(\boldsymbol r_c),\\
\boldsymbol e&=[P_x\operatorname{Mat}(\boldsymbol e_c)]P_t^\top,\\
\boldsymbol U^{k+2/3}&=\boldsymbol U^{k+1/3}+\operatorname{Vec}(\boldsymbol e),\\
\boldsymbol U^{k+1}&=S_\eta(\boldsymbol b,\boldsymbol U^{k+2/3},s_2).
\end{aligned}
\right. \tag{4.34}
$$

$K_c$ is rediscretized with $\Delta T=2\Delta t$ and $\Delta X=2\Delta x$ and has the same block structure as $K$. For $N_x=2^{l_x}-1,N_t=2^{l_t}-1$, the next sizes are $2^{l_x-1}-1$ and $2^{l_t-1}-1$. Recursion gives the multilevel method.

### Difference from early parabolic multigrid

Hackbusch-style parabolic multigrid uses the following pointwise,
time-sequential Gauss–Seidel smoother:

$$
\boldsymbol U^{\mathrm{new}}
=S_{GS}(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s):
\left\{
\begin{aligned}
&\text{for }n=0,\ldots,N_t-1,\\
&\qquad \boldsymbol u_{n+1}^0=\boldsymbol u_{n+1}^{\mathrm{ini}},\\
&\qquad \text{for }j=0,\ldots,s-1,\\
&\qquad\quad
(D+L)\Delta\boldsymbol u_{n+1}^j
=\widetilde{\boldsymbol f}_n+r_2\boldsymbol u_n^s
-r_1\boldsymbol u_{n+1}^j,\\
&\qquad\quad
\boldsymbol u_{n+1}^{j+1}
=\boldsymbol u_{n+1}^j+\Delta\boldsymbol u_{n+1}^j,\\
&\qquad \boldsymbol u_{n+1}^{\mathrm{new}}
=\boldsymbol u_{n+1}^s,
\end{aligned}
\right. \tag{4.35}
$$

Here $\boldsymbol u_0^s=\boldsymbol u_0$, and both
$\boldsymbol U^{\mathrm{ini}}$ and $\boldsymbol U^{\mathrm{new}}$
include the initial value. The matrices $D$ and $L$ are the diagonal
and **upper** triangular parts of $r_1$. The next time point depends on
the completed smoothed value at the current point. This works rapidly
for heat when only space is coarsened but becomes slow under
simultaneous space–time coarsening. Horton and Vandewalle (1995)
improved that behavior by reading time as a strongly advective
direction and designing special components. STMG's decisive change
(Gander and Neumüller 2016) is temporal block Jacobi: all time points
are smoothed concurrently using otherwise standard multigrid
components.

### Starting point of local Fourier analysis

Ignoring initial and boundary conditions, consider the error mode

$$
u_{n,m}^j=C_{\omega,\xi}^j
e^{i\omega n\Delta t}e^{i\xi m\Delta x}. \tag{4.36}
$$

For the one-dimensional heat equation,

$$
A=\frac1{\Delta x^2}\operatorname{Tri}(1,-2,1),
\qquad r_1=I_x-\Delta tA,\quad r_2=I_x.
$$

Block Jacobi on the homogeneous error is

$$
r_1(\boldsymbol u_{n+1}^{j+1}-\boldsymbol u_{n+1}^{j})
=-\eta(r_1\boldsymbol u_{n+1}^j-r_2\boldsymbol u_n^j). \tag{4.37}
$$

The spatial stencil has symbol

$$
Au_{n+1,m}^l
=\frac{2(\cos(\xi\Delta x)-1)}{\Delta x^2}
C_{\omega,\xi}^l
e^{i\omega(n+1)\Delta t}e^{i\xi m\Delta x},
\qquad l=j,j+1. \tag{4.38}
$$

Substitution gives $C^{j+1}=\rho C^j$ with

$$
\rho(\omega,\xi,\eta)
=1-\eta\left(
1-\frac{e^{-i\omega\Delta t}}
{1+\frac{2\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))}
\right). \tag{4.39}
$$

The amplification factor is $|\rho|$. To decide whether time
coarsening is admissible, define

$$
\Theta_t^{\mathrm{high}}
=(-\pi,-\pi/2)\cup(\pi/2,\pi)
$$

and the smoothing factor

$$
\mu_t(\eta)=
\sup_{\substack{\xi\Delta x\in(-\pi,\pi)\\
\omega\Delta t\in\Theta_t^{\mathrm{high}}}}
|\rho(\omega,\xi,\eta)|.
$$

Including the zero mode in the maximum would give
$\rho(0,0,\eta)=1$ for every $\eta$, so it cannot determine a unique
damping parameter.

### Theorem 4.9: optimal damping for backward Euler

Theorem 4.9 is proved in Gander and Lunet (2024). “Optimal” means
minimizing the high-temporal-frequency factor $\mu_t(\eta)$ above. For
centered differences in space and backward Euler in time, the damping
that always permits time coarsening is

$$
\eta_{\mathrm{opt}}=\frac12.
$$

Every temporal high frequency $\omega\in\pm(\pi/(2\Delta t),\pi/\Delta t)$ has amplification at most $1/\sqrt2$. If

$$
\frac{\Delta t}{\Delta x^2}\ge\frac1{\sqrt2},
$$

the same bound holds for spatial high frequencies $\xi\in\pm(\pi/(2\Delta x),\pi/\Delta x)$, so space can also be coarsened.

For centered advection–diffusion,

$$
A=\frac{\nu}{\Delta x^2}\operatorname{Tri}(1,-2,1)
+\frac1{2\Delta x}\operatorname{Tri}(-1,0,1),
$$

and

$$
\rho(\omega,\xi,\eta)
=1-\eta\left(
1-\frac{e^{-i\omega\Delta t}}
{1+\frac{2\nu\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))
+i\frac{\Delta t}{\Delta x}\sin(\xi\Delta x)}
\right). \tag{4.40}
$$

> [!warning] Source check: the ADE matrix versus (4.40)
> The $A$ above and (4.40) cannot both be right. Under the paper's own convention, $\operatorname{Tri}(-1,0,1)/(2\Delta x)$ acting on a Fourier mode gives $+i\sin(\xi\Delta x)/\Delta x$, so the symbol of $r_1=I_x-\Delta tA$ carries the opposite imaginary sign from the $+i\frac{\Delta t}{\Delta x}\sin(\xi\Delta x)$ in the denominator of (4.40). Model problem (2.5) is $\partial_tu+\partial_xu-\nu\partial_{xx}u=g$, so $A$ discretizes $-\partial_x+\nu\partial_{xx}$ and should read $A=\frac{\nu}{\Delta x^2}\operatorname{Tri}(1,-2,1)-\frac1{2\Delta x}\operatorname{Tri}(-1,0,1)$, with which (4.40) holds as printed. The display above keeps the paper's typesetting. Because $\rho_{\max}$ is a maximum over $\xi\Delta x\in(-\pi,\pi)$ and $\sin$ is odd, the sign does not affect any number quoted on this page.

![Original Figure 4.18: maximum high-frequency smoothing factors for ADE at three viscosities](assets/papers/time-parallelization/source-figures/figure-4-18.svg)

The panels from left to right use $\nu=0.1,0.01,0.001$; each compares
$\Delta x=\Delta t=1/64,1/128,1/256$. The vertical quantity should be
read as

$$
\rho_{\max}=
\max_{\substack{\Delta x\xi\in(-\pi,\pi)\\
\Delta t\omega\in\Theta_t^{\mathrm{high}}}}
|\rho(\omega,\xi,\eta)|.
$$

The journal caption shows only the positive half of the temporal set
and omits the modulus; conjugate symmetry gives the expression above.
Within each panel the three grid minima nearly coincide. As $\nu$
falls, the minimizer drifts slowly right of $0.5$ and the minimum rises
from roughly $0.71$ to $0.79$. Thus $\eta=\tfrac12$ is a robust
starting point, not the exact optimizer for every viscosity.

### Damping, sweep count, and time-integrator dependence

![Original Figure 4.19: error after five, ten, and fifteen cycles versus damping](assets/papers/time-parallelization/source-figures/figure-4-19.svg)

With one block-Jacobi sweep per cycle, panel (a) is heat and panel (b) is ADE at $\nu=0.01$; each reports errors after 5, 10, and 15 cycles. More cycles sharpen the low-error valley. The heat valley is broad, while the 15-cycle ADE curve has a flat basin running from roughly $0.3$ to $0.5$. Figure 4.19 therefore supports $\eta=1/2$ as a robust heuristic, not as the exact finite-grid optimizer for every fixed cycle count.

![Original Figure 4.20: one versus three block-Jacobi smoothing sweeps](assets/papers/time-parallelization/source-figures/figure-4-20.svg)

Panels (a) and (b) use one and three block-Jacobi sweeps, respectively; each compares heat with ADE at $\nu=0.1,0.01,0.001$. Three sweeps cost more per cycle and sharply reduce the cycle count. Several ADE curves in the right panel steepen at later cycles, revealing a superlinear stage. ADE remains slower than heat, but the additional smoothing substantially reduces its viscosity sensitivity.

![Original Figure 4.21: STMG with the trapezoidal rule across damping values and sweep counts](assets/papers/time-parallelization/source-figures/figure-4-21.svg)

With trapezoidal time integration, the heat row, group (a), uses 3, 5, and 10 sweeps. Every damping scan either remains at large error or becomes unstable; more smoothing does not recover the backward-Euler behavior. The ADE row, group (b), uses 2, 3, and 4 sweeps at $\nu=0.01$. Its error valley deepens with additional smoothing. The paper reports $\eta\approx0.8$ as the useful value for that row; read off the plots, the per-panel optimum drifts from about $0.79$ to about $0.92$ as the sweep count rises. The backward-Euler value $1/2$ depends on L-stable high-frequency damping and does not transfer directly.

![Original Table 4.1: weak and strong scaling of STMG on a three-dimensional heat equation](assets/papers/time-parallelization/source-figures/table-4-1.svg)

The table is taken from Gander and Neumüller (2016). Weak scaling grows from 1 core, 2 time steps, and 59,768 degrees of freedom to 262,144 cores, 524,288 steps, and 15,667,822,592 degrees of freedom. The iteration count remains seven and wall time stays near 28.8–30.0 seconds. The reference column, classical time stepping with the best possible parallelization in space only, grows from 19.0 to 4,988,060 seconds. There are two strong-scaling blocks: 512 time steps with 15,300,608 degrees of freedom drops from 7,635.2 to 30.0 seconds, and 524,288 time steps with 15,667,822,592 degrees of freedom drops from 15,205.9 to 30.0 seconds. The result combines space–time concurrency with a grid-independent cycle count.

### Nonlinear system and FAS

Consider

$$
\boldsymbol u'=f(\boldsymbol u),
\qquad \boldsymbol u(0)=\boldsymbol u_0,
\quad t\in(0,T). \tag{4.41}
$$

For the stacked state, define explicitly

$$
\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top,
\qquad
f(\boldsymbol U)=
\left(
f(\boldsymbol u_1)^\top,\ldots,f(\boldsymbol u_{N_t})^\top
\right)^\top.
$$

The journal display mistakenly writes the state blocks themselves in
place of the $f(\boldsymbol u_n)$ blocks.

The linear-theta method gives

$$
\underbrace{(B\otimes I_x)\boldsymbol U
-\Delta t(\widetilde B\otimes I_x)f(\boldsymbol U)}_{K(\boldsymbol U)}
=\boldsymbol b, \tag{4.42}
$$

with

$$
B=\begin{bmatrix}1\\-1&1\\&\ddots&\ddots\\&&-1&1\end{bmatrix},
\qquad
\widetilde B=
\begin{bmatrix}\theta\\1-\theta&\theta\\&\ddots&\ddots\\&&1-\theta&\theta\end{bmatrix},
$$

$$
\boldsymbol b=
(\boldsymbol u_0^\top+\Delta t(1-\theta)f(\boldsymbol u_0)^\top,0,\ldots,0)^\top.
$$

The nonlinear block-Jacobi smoother is

$$
\boldsymbol U^{\mathrm{new}}
=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s):
\left\{
\begin{aligned}
\widetilde{\boldsymbol U}^0&=\boldsymbol U^{\mathrm{ini}},\\
\Delta\widetilde{\boldsymbol U}^j
-\Delta t\theta f(\Delta\widetilde{\boldsymbol U}^j)
&=\eta[\boldsymbol b-K(\widetilde{\boldsymbol U}^j)],\\
\widetilde{\boldsymbol U}^{j+1}
&=\widetilde{\boldsymbol U}^j+\Delta\widetilde{\boldsymbol U}^j,\\
\boldsymbol U^{\mathrm{new}}&=\widetilde{\boldsymbol U}^s.
\end{aligned}
\right. \tag{4.43, as printed}
$$

> [!warning] Source check: what “block Jacobi” means in (4.43)
> If $\Delta\widetilde{\boldsymbol U}^j$ is the correction added to the
> current state, a consistent nonlinear correction for
> $D(\boldsymbol U)=\boldsymbol U-\Delta t\theta f(\boldsymbol U)$
> would satisfy
>
> $$
> \Delta\widetilde{\boldsymbol U}^j
> -\Delta t\theta\left[
> f(\widetilde{\boldsymbol U}^j+\Delta\widetilde{\boldsymbol U}^j)
> -f(\widetilde{\boldsymbol U}^j)
> \right]
> =\eta[\boldsymbol b-K(\widetilde{\boldsymbol U}^j)].
> $$
>
> The printed $f(\Delta\widetilde{\boldsymbol U}^j)$ agrees only in
> special cases such as a linear $f$. Equation (4.43) is therefore
> safer to read as a separately defined nonlinear preconditioned
> iteration, not the exact inverse of the nonlinear block diagonal.

Each temporal nonlinear block can use an inner Newton solve and run
concurrently. LFA no longer supplies an optimized eta.

The two-level full approximation scheme is

$$
\left\{
\begin{aligned}
\boldsymbol U^{k+1/3}&=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^k,s_1),\\
\boldsymbol r&=\boldsymbol b-K(\boldsymbol U^{k+1/3}),\\
\boldsymbol r_c&=\operatorname{Vec}\!\left(
[R_x\operatorname{Mat}(\boldsymbol r)]R_t^\top
\right),\\
\boldsymbol U_c^{k+1/3}&=\operatorname{Vec}\!\left(
[R_x\operatorname{Mat}(\boldsymbol U^{k+1/3})]R_t^\top
\right),\\
K_c(\boldsymbol U_c^{k+2/3})
&=\boldsymbol r_c+K_c(\boldsymbol U_c^{k+1/3}),\\
\boldsymbol e_c&=\boldsymbol U_c^{k+2/3}-\boldsymbol U_c^{k+1/3},\\
\boldsymbol e&=\operatorname{Vec}\!\left(
[P_x\operatorname{Mat}(\boldsymbol e_c)]P_t^\top
\right),\\
\boldsymbol U^{k+2/3}&=\boldsymbol U^{k+1/3}+\boldsymbol e,\\
\boldsymbol U^{k+1}&=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^{k+2/3},s_2).
\end{aligned}
\right. \tag{4.44}
$$

FAS (Brandt 1977) solves for a full coarse approximation and uses
$\boldsymbol r_c+K_c(\boldsymbol U_c)$ to maintain nonlinear
consistency, unlike a linear coarse error equation. The explicit
`Vec` operations above keep the coarse residual and every argument of
$K_c$ as vectors; the journal omits these conversions between matrix
and vector views.

![Original Figure 4.22: Burgers STMG with two block-Jacobi sweeps](assets/papers/time-parallelization/source-figures/figure-4-22.svg)

The experiment uses two sweeps and the empirically best $\eta=1/4$. The $\nu=1$ curve crosses the plotted discretization-error line after roughly four cycles and continues toward $10^{-4}$; the $\nu=0.1$ curve remains above that line after 16 cycles. Nonlinear STMG therefore retains the viscosity sensitivity seen in the linear case (the $\nu$ sweep of Figure 4.20).

Among the surveyed parabolic methods, STMG is judged the most
effective, although it is more intrusive than Parareal. Figures 4.20
and 4.22 actually test low-diffusion ADE and viscous Burgers rather
than strictly hyperbolic equations; they directly establish
deterioration near a hyperbolic limit. The authors extrapolate from
this evidence that truly hyperbolic problems still need further work.
Figure 4.21 separately leaves dependence on the time integrator as an
open question even for parabolic problems.

## Equation, theorem, and figure audit

| Source item                             | Paper section | Coverage                                                             |
| --------------------------------------- | ------------- | -------------------------------------------------------------------- |
| (4.30)–(4.31)                           | 4.6           | generic one-step formula and all-at-once matrix                      |
| (4.32)–(4.34)                           | 4.6           | parallel block Jacobi, space–time transfer, full two-level cycle     |
| (4.35)                                  | 4.6           | early sequential Gauss–Seidel smoother and contrast                  |
| (4.36)–(4.40), Theorem 4.9, Figure 4.18 | 4.6           | full LFA, heat/ADE symbols, optimal damping, coarsening condition    |
| Figures 4.19–4.21, Table 4.1            | 4.6           | damping, sweep count, integrator dependence, strong/weak scaling     |
| (4.41)–(4.44), Figure 4.22              | 4.6           | nonlinear all-at-once system, concurrent smoother, FAS, Burgers test |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4.6, pp. 472–481.
