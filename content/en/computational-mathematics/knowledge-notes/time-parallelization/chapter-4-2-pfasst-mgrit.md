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
> This page follows Sections 4.3–4.4 (pp. 452–460). It covers equations (4.10)–(4.13), Theorems 4.5–4.6, and Figures 4.6–4.11. PFASST is taken down to its collocation matrices, transfer operators, and SDC approximation; MGRiT is derived as overlapping Parareal.

## 4.3 PFASST

### Origin and two node sets

PFASST was introduced by Emmett and Minion (2012). Its precursor replaced a complete Parareal fine solve with one spectral deferred-correction sweep. The later algebraic description views PFASST as temporal multigrid with collocation on the fine level and low-order SDC as a smoother.

Split $(0,T)$ into $N_t$ large intervals. On every $[T_n,T_{n+1}]$, define $M_f$ fine and $M_c$ coarse nodes

$$
t_{n,m}^{f}=T_n+\tau_m^f\Delta t,
\qquad
t_{n,m}^{c}=T_n+\tau_m^c\Delta t,
$$

with $\tau_0^{f,c}=0$, $\tau_{M_{f,c}}^{f,c}=1$, and $M_f>M_c$.

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

The block-iteration form of PFASST is

$$
\boldsymbol u_{n+1}^{k+1}
=B_{01}\boldsymbol u_{n+1}^k
+B_{10}(\boldsymbol\chi\boldsymbol u_n^{k+1}+\Delta t\boldsymbol b_n)
+B_{00}(\boldsymbol\chi\boldsymbol u_n^k+\Delta t\boldsymbol b_n),
$$

with

$$
\begin{aligned}
B_{01}&=[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]
(I_f-\widetilde\Phi_f^{-1}\Phi_f),\\
B_{10}&=T_{c\to f}\Phi_c^{-1}T_{f\to c},\\
B_{00}&=[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]\widetilde\Phi_f^{-1}.
\end{aligned}
$$

These blocks represent fine SDC smoothing, new-iterate coarse propagation, and old-iterate correction. The matrix $\widetilde\Phi_f$ is an inexpensive approximation of the fine collocation matrix.

### SDC approximation (4.11) and Figure 4.6

The paper constructs $\widetilde\Phi_f$ from implicit Euler between adjacent fine nodes:

$$
\frac{\boldsymbol u_{n,m+1}-\boldsymbol u_{n,m}}
{\Delta t(\tau_{m+1}^f-\tau_m^f)}
=A\boldsymbol u_{n,m+1}+\boldsymbol g(t_{n,m+1}^f),
\quad m=0,\ldots,M_f-1. \tag{4.11}
$$

Its matrix combines node-spacing diagonal blocks with a lower bidiagonal difference matrix and $\Delta tA$. Solving it is one low-order SDC sweep rather than a full collocation solve.

The experiment uses $T=3$, periodic boundaries, zero initial data, source (2.4) with $\sigma=1000$, $\Delta x=1/128$, and $\Delta t=1/64$. Fine Radau IIA nodes are

$$
\left\{0,\frac{4-\sqrt6}{10},\frac{4+\sqrt6}{10},1\right\},
$$

and coarse nodes are $\{0,1/3,1\}$. The paper gives the resulting $Q_f,Q_c$ and numerical transfer matrices, all determined by the node sets and Lagrange interpolation.

![Original Figure 4.6: PFASST error for the heat equation and advection–diffusion at three viscosities](assets/papers/time-parallelization/source-figures/figure-4-6.svg)

Heat converges fastest. As advection–diffusion viscosity decreases, persistent high frequencies become harder for the coarse collocation level to represent, reproducing the coarse–fine mismatch mechanism of Parareal.

## 4.4 MGRiT

### As overlapping Parareal

MGRiT also admits algebraic-multigrid and block-iteration interpretations. Two-level FCF relaxation on a nonlinear system is

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

The extra F relaxation advances old-iterate information across one additional coarse interval. The global error therefore vanishes in at most $\lceil N_t/2\rceil$ iterations. More generally, $F(CF)^\nu$ corresponds to overlap $\nu\Delta T$.

### Theorem 4.5: long-time modal factor

Under the notation of Theorem 4.2 and $|R_g(z)|<1$,

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

The extra F relaxation contributes one fine-propagation contraction and costs one additional parallel fine solve.

![Original Figure 4.8: complex-plane convergence regions for one MGRiT iteration and two work-matched Parareal iterations](assets/papers/time-parallelization/source-figures/figure-4-8.svg)

Figure 4.8 uses backward Euler $R_g(z)=1/(1-z)$ and exact fine propagation $R_f(z)=e^z$. The top row, panel group (a), is MGRiT; the bottom row, group (b), plots $\varrho_{l,\mathrm{Parareal}}^2$ at the same cost of two fine solves. The columns from left to right use $\widehat\varrho=0.2,0.4,0.6$. The shaded contours are close within every column, showing why fine-propagation count is the appropriate unit for this comparison. SDIRK fine methods give the same qualitative result.

### Theorem 4.6: work-normalized constants

For an L-stable fine method and $J=O(1)$, the worst negative-real-axis factors with backward Euler coarse propagation are

$$
\max\varrho_l\approx
\begin{cases}
0.2984,&\text{Parareal},\\
0.1115,&\text{FCF-MGRiT}.
\end{cases}
$$

With second-order Lobatto IIIC coarse propagation they are

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

Figure 4.9 places MGRiT in the top row, group (a), and Parareal in the bottom row, group (b). The columns are heat, ADE with $\nu=0.1$, and ADE with $\nu=0.01$. The annotated maxima are $0.08375,0.2718,0.9021$ for MGRiT and $0.2822,0.4453,0.9986$ for Parareal. In the first two columns, the MGRiT factor is roughly comparable to the square of the Parareal factor. In the third column, both are close to one because the coarse propagator no longer represents the long-lived weak-diffusion modes accurately.

![Original Figure 4.10: measured Parareal and MGRiT errors after matching the number of fine solves](assets/papers/time-parallelization/source-figures/figure-4-10.svg)

The four panels in Figure 4.10 are heat, ADE with $\nu=0.1$, ADE with $\nu=0.01$, and ADE with $\nu=0.002$. Each plotted Parareal unit contains two iterations. The first two panels show closely matched curves. At $\nu=0.01$, both are slow and Parareal deteriorates more because it inserts the inaccurate coarse solve after each fine solve. At $\nu=0.002$, both diverge; the maximum modal factors are $1.4211$ and $1.2812$. The modal quantities in Figure 4.9 and measured curves in Figure 4.10 therefore correspond case by case across the four diffusion regimes.

The nonlinear Burgers test uses homogeneous Dirichlet data, the same initial condition, $T=5$, $\Delta T=1/16$, $\Delta x=1/160$, $J=10$, centered differences, backward Euler coarse propagation, and SDIRK22 fine propagation.

![Original Figure 4.11: work-matched comparison on Burgers' equation at three viscosities](assets/papers/time-parallelization/source-figures/figure-4-11.svg)

The three panels from left to right use $\nu=0.5,0.01,0.002$. In the first two, one FCF-MGRiT iteration behaves like two Parareal iterations. At $\nu=0.002$, both curves pass through a long slow phase before their later rapid decrease. The two methods deteriorate together as viscosity falls, matching the nonlinear Lipschitz analysis.

## Equation, theorem, and figure audit

| Source item                                  | Paper section | Coverage                                                                  |
| -------------------------------------------- | ------------- | ------------------------------------------------------------------------- |
| (4.10)                                       | 4.3           | fine/coarse nodes, copying matrix, collocation system                     |
| Lagrange transfer and PFASST block iteration | 4.3           | $T_{c\to f},T_{f\to c}$ and $B_{01},B_{10},B_{00}$                        |
| (4.11), Figure 4.6                           | 4.3           | implicit-Euler SDC, Radau nodes, PFASST test                              |
| (4.12), Figure 4.7                           | 4.4           | FCF update, two fine solves, overlap, finite-step property                |
| (4.13), Theorem 4.5, Figure 4.8              | 4.4           | MGRiT factor, Parareal relation, work-matched regions                     |
| Theorem 4.6                                  | 4.4           | four worst factors and squared-factor comparison                          |
| Figures 4.9–4.11                             | 4.4           | linear spectra, measured error, nonlinear Burgers, weak-diffusion failure |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 4.3–4.4, pp. 452–460.
