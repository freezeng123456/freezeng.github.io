---
title: "4.3–4.4: PFASST and MGRiT"
description: Complete structures from collocation, SDC smoothing, and transfer operators to FCF as overlapping Parareal and a work-normalized comparison
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-2-pfasst-mgrit
tags:
  - parallel-in-time
  - PFASST
  - MGRiT
---

> [!note] Reading scope
> This page follows Sections 4.3–4.4 (pp. 452–461). It covers equations (4.10)–(4.13), Theorems 4.5–4.6, and Figures 4.6–4.11. PFASST is taken down to its collocation matrices, transfer operators, and SDC approximation; MGRiT is derived as overlapping Parareal.

## 4.3 PFASST

### Origin and two node sets

PFASST was introduced by Emmett and Minion (2012). Its precursor is due to Minion (2010), who replaced a complete Parareal fine solve with one spectral deferred-correction sweep (SDC; Dutt et al. 2000). Minion et al. (2015) gave the algebraic description, and Bolten, Moser and Speck (2017) read PFASST as temporal multigrid, with a convergence analysis in Bolten, Moser and Speck (2018): collocation on the fine level, low-order SDC as a smoother. Applications appear in Speck et al. (2012, 2014). The paper also cautions that a clear description and theoretical analysis of PFASST are rather challenging.

Split $(0,T)$ into $N_t$ intervals. This section reuses
$\Delta t=T_{n+1}-T_n$ for one PFASST interval; it is not the fine
Parareal step from the preceding page. On every $[T_n,T_{n+1}]$,
define $M_f$ fine and $M_c$ coarse collocation unknowns

$$
t_{n,m}^{f}=T_n+\tau_m^f\Delta t,
\qquad
t_{n,m}^{c}=T_n+\tau_m^c\Delta t,
$$

where $\tau_0^{f,c}=0$ is the separately listed initial location,
$\tau_{M_{f,c}}^{f,c}=1$, and $M_f>M_c$. Thus $M_f=3$ means three
collocation unknowns plus $\tau_0^f$, for four local time locations.

### Collocation equation (4.10)

For $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$, either level satisfies

$$
\boldsymbol u_{n,m}
=\boldsymbol u_{n,0}
+\Delta t\sum_{j=1}^{M}q_{m,j}
[A\boldsymbol u_{n,j}+\boldsymbol g(t_{n,j})],
\qquad m=1,\ldots,M. \tag{4.10}
$$

Let

$$
Q=(q_{m,j}),\qquad
\boldsymbol u_n=(\boldsymbol u_{n,1}^\top,\ldots,\boldsymbol u_{n,M}^\top)^\top,
$$

$$
\boldsymbol g_n=
\left(
\boldsymbol g(t_{n,1})^\top,\ldots,
\boldsymbol g(t_{n,M})^\top
\right)^\top,
\qquad
\boldsymbol b_n=(Q\otimes I_x)\boldsymbol g_n,
$$

$$
\chi=
\begin{bmatrix}
0&\cdots&0&1\\
\vdots&&\vdots&\vdots\\
0&\cdots&0&1
\end{bmatrix},
\qquad \boldsymbol\chi=\chi\otimes I_x.
$$

The copying matrix passes the terminal node of one interval to every collocation node of the next:

$$
\boldsymbol u_n
=\Delta t(Q\otimes A)\boldsymbol u_n
+\boldsymbol\chi\boldsymbol u_{n-1}
+\Delta t\boldsymbol b_n.
$$

Define

$$
\Phi_f=I_f-\Delta tQ_f\otimes A,
\qquad
\Phi_c=I_c-\Delta tQ_c\otimes A.
$$

Here

$$
I_f=I_{M_f}\otimes I_x,
\qquad I_c=I_{M_c}\otimes I_x.
$$

Then

$$
\boldsymbol u_n^f
=\Phi_f^{-1}(\boldsymbol\chi_f\boldsymbol u_{n-1}^f+\Delta t\boldsymbol b_n^f),
\qquad
\boldsymbol u_n^c
=\Phi_c^{-1}(\boldsymbol\chi_c\boldsymbol u_{n-1}^c+\Delta t\boldsymbol b_n^c).
$$

### Lagrange transfer between levels

Coarse-node values define

$$
p^c(\tau;\boldsymbol u^c)
=\sum_{m=1}^{M_c}u_m^cL_m^c(\tau),
\qquad
L_m^c(\tau)=
\prod_{\substack{j=1\\j\ne m}}^{M_c}
\frac{\tau-\tau_j^c}{\tau_m^c-\tau_j^c}.
$$

Evaluation at the fine nodes yields $T_{c\to f}$; evaluating fine Lagrange basis functions at coarse nodes yields $T_{f\to c}$. Both matrices are tensored with $I_x$, so transfer acts only in the collocation-node direction.

> [!warning] Source check: the Lagrange basis
> The journal and arXiv versions misprint both the node index in the
> numerator and the denominator orientation. The numerical transfer
> matrices that follow correspond to the standard basis above, not to
> the printed expression.

The block-iteration form of PFASST is

$$
\boldsymbol u_{n+1}^{k+1}
=\mathbf B_1^0\boldsymbol u_{n+1}^k
+\mathbf B_0^1(\boldsymbol\chi\boldsymbol u_n^{k+1}+\Delta t\boldsymbol b_n^f)
+\mathbf B_0^0(\boldsymbol\chi\boldsymbol u_n^k+\Delta t\boldsymbol b_n^f),
$$

with

$$
\begin{aligned}
\mathbf B_1^0&=[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]
(I_f-\widetilde\Phi_f^{-1}\Phi_f),\\
\mathbf B_0^1&=T_{c\to f}\Phi_c^{-1}T_{f\to c},\\
\mathbf B_0^0&=[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]\widetilde\Phi_f^{-1}.
\end{aligned}
$$

The subscript records the time offset and the superscript the iteration offset, following the paper's notation. These blocks represent fine SDC smoothing, new-iterate coarse propagation, and old-iterate correction. The matrix $\widetilde\Phi_f$ is an inexpensive approximation of the fine collocation matrix.

### SDC approximation (4.11) and Figure 4.6

The paper constructs $\widetilde\Phi_f$ from implicit Euler between adjacent fine nodes:

$$
\frac{\boldsymbol u_{n,m+1}-\boldsymbol u_{n,m}}
{\Delta t(\tau_{m+1}^f-\tau_m^f)}
=A\boldsymbol u_{n,m+1}+\boldsymbol g(t_{n,m+1}^f),
\quad m=0,\ldots,M_f-1. \tag{4.11}
$$

Its matrix is the Kronecker product of a lower bidiagonal difference matrix with $I_x$, minus the Kronecker product of the node-gap diagonal matrix with $\Delta tA$:

$$
\widetilde\Phi_f=
\begin{bmatrix}
1\\-1&1\\&\ddots&\ddots\\&&-1&1
\end{bmatrix}\otimes I_x
-\Delta t
\begin{bmatrix}
\tau_1^f-\tau_0^f\\
&\ddots\\
&&\tau_{M_f}^f-\tau_{M_f-1}^f
\end{bmatrix}\otimes A.
$$

Solving it is one low-order SDC sweep rather than a full collocation solve.

The experiment uses $T=3$, periodic boundaries, zero initial data, source (2.4) with $\sigma=1000$, $\Delta x=1/128$, and $\Delta t=1/64$. Fine Radau IIA nodes are

$$
\left\{0,\frac{4-\sqrt6}{10},\frac{4+\sqrt6}{10},1\right\},
$$

that is $M_f=3$; the coarse nodes are $\{0,1/3,1\}$, that is
$M_c=2$, with Radau IIA on both levels. The weight matrices are

$$
Q_f=
\begin{bmatrix}
\frac{88-7\sqrt6}{360}
&\frac{296-169\sqrt6}{1800}
&\frac{-2+3\sqrt6}{225}\\
\frac{296+169\sqrt6}{1800}
&\frac{88+7\sqrt6}{360}
&\frac{-2-3\sqrt6}{225}\\
\frac{16-\sqrt6}{36}
&\frac{16+\sqrt6}{36}
&\frac19
\end{bmatrix},
\qquad
Q_c=
\begin{bmatrix}
\frac5{12}&-\frac1{12}\\
\frac34&\frac14
\end{bmatrix}.
$$

The same nodes give the numerical transfer matrices

$$
T_{c\to f}=
\begin{bmatrix}1.2674&-0.2674\\0.5325&0.4674\\0&1\end{bmatrix}
\otimes I_x,
\qquad
T_{f\to c}=
\begin{bmatrix}0.5018&0.6833&-0.1851\\0&0&1\end{bmatrix}
\otimes I_x,
$$

all determined by the node sets and Lagrange interpolation.

![Original Figure 4.6: PFASST error for the heat equation and advection–diffusion at three viscosities](assets/papers/time-parallelization/source-figures/figure-4-6.svg)

The three advection–diffusion curves use
$\nu=0.1,10^{-3},10^{-4}$, with the iteration index running to $300$.
Heat and $\nu=0.1$ reach $10^{-12}$ near $k\approx90$–95,
$\nu=10^{-3}$ near $120$, and $\nu=10^{-4}$ near $290$. As viscosity
weakens, the coarse propagator no longer represents the dominant
advection well enough, so the iteration count rises.

> [!note] Site supplement: relation to Parareal
> The source already compares the weak-diffusion deterioration
> qualitatively with Parareal and MGRiT. Interpreting it as persistent
> high frequencies becoming harder for the coarse collocation level to
> represent is this site's more specific mechanism-level explanation.

## 4.4 MGRiT

### As overlapping Parareal

MGRiT was introduced by Falgout et al. (2014) and also admits algebraic-multigrid and block-iteration interpretations. Two-level FCF relaxation on a nonlinear system is

$$
\boldsymbol u_0^{k+1}=\boldsymbol u_0,
\qquad
\boldsymbol u_1^{k+1}=\mathcal F(T_0,T_1,\boldsymbol u_0),
$$

$$
\begin{aligned}
\boldsymbol u_{n+1}^{k+1}
={}&\mathcal F\!\left(T_n,T_{n+1},
\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)\right)\\
&+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})\\
&-\mathcal G\!\left(T_n,T_{n+1},
\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)\right),
\quad n=1,\ldots,N_t-1. \tag{4.12}
\end{aligned}
$$

One iteration uses two fine solves, compared with one for Parareal.

![Original Figure 4.7: FCF-MGRiT as Parareal with one coarse-interval overlap](assets/papers/time-parallelization/source-figures/figure-4-7.svg)

The extra CF relaxation advances old-iterate information across one additional coarse interval. The global error therefore vanishes in at most $\lceil N_t/2\rceil$ iterations. More generally, $F(CF)^\nu$ corresponds to overlap $\nu\Delta T$. These results are Theorem 5, Theorem 6 and Corollary 1 of Gander et al. (2018b); Theorem 6 is a superlinear convergence result for nonlinear problems.

### Theorem 4.5: long-time modal factor

Theorem 4.5 is taken from Dobrev et al. (2017). Under the notation and
assumptions of Theorem 4.2, $\mathcal F$ and $\mathcal G$ are one-step
integrators, $\boldsymbol u_n$ is the sequential fine solution,
$\boldsymbol e_n^k=V_A(\boldsymbol u_n^k-\boldsymbol u_n)$, and $A$ is
diagonalizable with $\sigma(A)\subset\mathbb C_-$. If $|R_g(z)|<1$,

$$
\max_n\|\boldsymbol e_n^k\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}\varrho_l^k(J,z)
\max_n\|\boldsymbol e_n^0\|_\infty,
$$

$$
\varrho_l(J,z)=
\frac{|R_f^J(z/J)|\,|R_g(z)-R_f^J(z/J)|}
{1-|R_g(z)|}. \tag{4.13}
$$

Hence

$$
\varrho_{l,\mathrm{MGRiT}}
=|R_f^J(z/J)|\,\varrho_{l,\mathrm{Parareal}}.
$$

Provided $|R_f(z/J)|\le1$ (equivalently,
$|R_f^J(z/J)|\le1$), the extra CF relaxation contributes one
fine-propagation contraction and costs one additional parallel fine
solve. A-stability of the fine method on the entire left half-plane is
sufficient.

![Original Figure 4.8: complex-plane convergence regions for one MGRiT iteration and two work-matched Parareal iterations](assets/papers/time-parallelization/source-figures/figure-4-8.svg)

Figure 4.8 uses backward Euler $R_g(z)=1/(1-z)$ and exact fine propagation $R_f(z)=e^z$. The top row, panel group (a), is MGRiT; the bottom row, group (b), plots $\varrho_{l,\mathrm{Parareal}}^2$ at the same cost of two fine solves. The columns from left to right use $\widehat\varrho=0.2,0.4,0.6$. The shaded contours are close within every column, showing why fine-propagation count is the appropriate unit for this comparison. SDIRK fine methods give the same qualitative result.

### Theorem 4.6: work-normalized constants

Theorem 4.6 is taken from Wu and Zhou (2019). For an L-stable fine
method and $J=\Delta T/\Delta t=O(1)$, let $s=-z\ge0$. With backward
Euler coarse propagation, $\max_{s\ge0}\varrho_l(J,-s)$ is

$$
\max\varrho_l\approx
\begin{cases}
0.2984,&\text{Parareal},\\
0.1115,&\text{FCF-MGRiT}.
\end{cases}
$$

With second-order Lobatto IIIC coarse propagation, the same
$\max_{s\ge0}\varrho_l(J,-s)$ is

$$
\max\varrho_l\approx
\begin{cases}
0.0817,&\text{Parareal},\\
0.0197,&\text{FCF-MGRiT}.
\end{cases}
$$

At equal cost in fine solves, one MGRiT iteration is slightly slower than two Parareal iterations: $0.2984^2=0.0890<0.1115$ and $0.0817^2=0.0067<0.0197$. A steeper per-iteration curve does not by itself imply better efficiency per fine solve.

### Figures 4.9–4.11: linear and nonlinear tests

The linear tests use homogeneous Dirichlet data for heat, periodic data for ADE, and

$$
u_0(x)=\sin^2(8\pi(1-x)^2),\quad
T=5,\quad J=20,\quad \Delta T=1/8,\quad \Delta x=1/160,
$$

with backward Euler coarse and SDIRK22 fine propagation.

![Original Figure 4.9: modal-factor distributions for heat and ADE at two viscosities](assets/papers/time-parallelization/source-figures/figure-4-9.svg)

Figure 4.9 places MGRiT in the top row, group (a), and Parareal in the bottom row, group (b). The columns are heat, ADE with $\nu=0.1$, and ADE with $\nu=0.01$. The annotated maxima are $0.08375,0.2718,0.9021$ for MGRiT and $0.2822,0.4453,0.9986$ for Parareal. In the first two columns, the MGRiT factor is roughly comparable to the square of the Parareal factor. In the third column, both are close to one; the reason the paper gives is that the coarse propagator is no longer good enough for small $\nu$ once advection dominates.

![Original Figure 4.10: measured Parareal and MGRiT errors after matching the number of fine solves](assets/papers/time-parallelization/source-figures/figure-4-10.svg)

The four panels in Figure 4.10 are heat, ADE with $\nu=0.1$, ADE with $\nu=0.01$, and ADE with $\nu=0.002$. Each plotted Parareal unit contains two iterations. The first two panels show closely matched curves. At $\nu=0.01$, both are slow and Parareal deteriorates more because it inserts the inaccurate coarse solve after each fine solve. A dash-dotted line marks the discretization truncation error $\max\{\Delta t^2,\Delta x^2\}$, beyond which one would not iterate in practice. At $\nu=0.002$, both diverge; the maximum modal factor is $1.4211$ for Parareal and $1.2812$ for MGRiT (these two numbers appear only in the paper's running text, since Figure 4.9 plots just the first three diffusion regimes).

The nonlinear case requires nonlinear solvers inside both
$\mathcal F$ and $\mathcal G$. The comparison also assumes Lipschitz
conditions on $\mathcal F$, $\mathcal G$, and their difference, so the
“one FCF iteration behaves like two Parareal iterations” statement is
not unconditional. The Burgers test uses homogeneous Dirichlet data,
the same initial condition, $T=5$, $\Delta T=1/16$,
$\Delta x=1/160$, $J=10$, centered differences, backward Euler coarse
propagation, and SDIRK22 fine propagation.

![Original Figure 4.11: work-matched comparison on Burgers' equation at three viscosities](assets/papers/time-parallelization/source-figures/figure-4-11.svg)

The three panels from left to right use $\nu=0.5,0.01,0.002$.
Provided the coarse solver remains reasonably accurate, one FCF-MGRiT
iteration behaves like two Parareal iterations at all three
viscosities, as in Figure 4.10. At $\nu=0.002$ both curves pass through
a long slow phase before their later rapid decrease, so the two methods
deteriorate together as viscosity falls.

## Equation, theorem, and figure audit

| Source item                                  | Paper section | Coverage                                                                  |
| -------------------------------------------- | ------------- | ------------------------------------------------------------------------- |
| (4.10)                                       | 4.3           | fine/coarse nodes, copying matrix, collocation system                     |
| Lagrange transfer and PFASST block iteration | 4.3           | $T_{c\to f},T_{f\to c}$ and $\mathbf B_1^0,\mathbf B_0^1,\mathbf B_0^0$   |
| (4.11), Figure 4.6                           | 4.3           | implicit-Euler SDC, Radau nodes, PFASST test                              |
| (4.12), Figure 4.7                           | 4.4           | FCF update, two fine solves, overlap, finite-step property                |
| (4.13), Theorem 4.5, Figure 4.8              | 4.4           | MGRiT factor, Parareal relation, work-matched regions                     |
| Theorem 4.6                                  | 4.4           | four worst factors and squared-factor comparison                          |
| Figures 4.9–4.11                             | 4.4           | linear spectra, measured error, nonlinear Burgers, weak-diffusion failure |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 4.3–4.4, pp. 452–461.
