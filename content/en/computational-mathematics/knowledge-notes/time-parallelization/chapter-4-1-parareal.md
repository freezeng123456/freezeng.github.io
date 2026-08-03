---
title: "4.1–4.2: Historical Context and Parareal"
description: A close reading from multiple shooting and coarse/fine propagators to the linear and nonlinear convergence bounds
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-1-parareal
tags:
  - parallel-in-time
  - Parareal
  - parabolic equations
---

> [!note] Reading scope
> This page follows Sections 4, 4.1, and 4.2 (pp. 443–452). It covers equations (4.1)–(4.9), Theorems 4.1–4.4, Remark 4.1, and Figures 4.1–4.5. The reasoning after each formula follows the proof order in the paper.

## 4.1 Historical development

### Why the paper separates methods designed for parabolic problems

Chapter 2 showed that parabolic solutions are local in time except at very low frequencies, whereas every hyperbolic frequency may propagate over a long horizon. The Chapter 3 methods address long-range temporal coupling and often work on parabolic problems as well. Their nonlinear forms retain specific difficulties: an effective OSWR Robin parameter can be hard to find, and the outer Newton iterations in ParaExp and ParaDiag may slow down or fail on long windows.

Chapter 4 exploits dissipation-induced temporal locality. Parareal, PFASST, MGRiT, and STMG handle linear and nonlinear problems effectively when diffusion is sufficient. Their performance degrades continuously toward weakly diffusive and hyperbolic regimes.

Parareal has roots in multiple shooting, waveform relaxation, and the noniterative precursor of Nievergelt (1964). Lions et al. (2001) introduced the modern method independently. Later descendants include PITA, PFASST, MGRiT, and Parareal–ParaDiag combinations. Space–time multigrid followed another line: early schemes could not coarsen time effectively, while the temporal block-Jacobi smoother of Gander and Neumüller (2016) enabled scalable STMG.

## 4.2 Parareal

### The update and the two temporal grids

Partition $[0,T]$ by $0=T_0<T_1<\cdots<T_{N_t}=T$. Let $\mathcal F$ be an accurate, expensive fine propagator and $\mathcal G$ a cheap coarse propagator. Starting from interface guesses $\boldsymbol u_n^0$, Parareal computes

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k). \tag{4.1}
$$

All fine solves at iteration $k$ are concurrent. The coarse propagation carrying the new index $k+1$ advances sequentially in $n$. The last two terms form a coarse prediction–correction and can also be viewed as a finite-difference Jacobian in a multiple-shooting Newton method.

![Original Figure 4.1: each coarse time step contains J fine time steps](assets/papers/time-parallelization/source-figures/figure-4-1.svg)

The analysis uses uniform grids with $\Delta T/\Delta t=J\ge2$, though nonuniform grids are possible. The target is the sequential discrete solution produced by $\mathcal F$; fixing $\mathcal F$ fixes the result to which Parareal converges.

### Theorem 4.1: a modal Toeplitz error iteration

For $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$, let $A=V_ADV_A^{-1}$ and let $R_g,R_f$ be the coarse and fine one-step stability functions. Fine propagation over one coarse interval is $R_f^J(z/J)$, with $z=\Delta T\lambda(A)$.

If $|R_g(z)|\le1$, then

$$
\max_{1\le n\le N_t}
\|V_A(\boldsymbol u_n^k-\boldsymbol u_n)\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}\|M^k(z)\|_\infty
\max_{1\le n\le N_t}
\|V_A(\boldsymbol u_n^0-\boldsymbol u_n)\|_\infty, \tag{4.2}
$$

where $\boldsymbol u_n$ is the sequential fine solution and

$$
M(z)=M_g^{-1}(z)[M_g(z)-M_f(z)],
$$

$$
M_g(z)=
\begin{bmatrix}
1\\-R_g(z)&1\\&\ddots&\ddots\\&&-R_g(z)&1
\end{bmatrix},
\qquad
M_f(z)=
\begin{bmatrix}
1\\-R_f^J(z/J)&1\\&\ddots&\ddots\\&&-R_f^J(z/J)&1
\end{bmatrix}. \tag{4.3}
$$

#### Proof chain

Equation (4.1) becomes

$$
\boldsymbol u_{n+1}^{k+1}
=R_f^J(\Delta TA/J)\boldsymbol u_n^k
+R_g(\Delta TA)\boldsymbol u_n^{k+1}
-R_g(\Delta TA)\boldsymbol u_n^k.
$$

The sequential fine solution satisfies the same identity after adding and subtracting the coarse term. Thus $\boldsymbol e_n^k=\boldsymbol u_n-\boldsymbol u_n^k$ obeys

$$
\boldsymbol e_{n+1}^{k+1}
=R_g(\Delta TA)\boldsymbol e_n^{k+1}
+[R_f^J(\Delta TA/J)-R_g(\Delta TA)]\boldsymbol e_n^k.
$$

Diagonalization separates the modes:

$$
\xi_{n+1}^{k+1}(z)
=R_g(z)\xi_n^{k+1}(z)
+[R_f^J(z/J)-R_g(z)]\xi_n^k(z),
$$

with

$$
\|V_A\boldsymbol e_n^k\|_\infty
=\max_{z\in\sigma(\Delta TA)}|\xi_n^k(z)|. \tag{4.4}
$$

Stacking the temporal components gives

$$
M_g(z)\boldsymbol\xi^{k+1}
=[M_g(z)-M_f(z)]\boldsymbol\xi^k,
$$

so $\boldsymbol\xi^k=M^k(z)\boldsymbol\xi^0$, which proves (4.2) after taking both maxima.

> [!note] Remark 4.1: preconditioner interpretation
> $M(z)=I_t-M_g^{-1}(z)M_f(z)$. The equation $M_fU=b$ is the fine all-at-once system, and $M_g$ is a coarse temporal preconditioner. Fine residual blocks are computed concurrently; application of the coarse preconditioner is sequential. Section 4.5 replaces this serial correction by a diagonalizable one.

The strictly lower triangular structure also proves exactness at the first $k$ coarse points after iteration $k$ and termination in at most $N_t$ iterations in exact arithmetic.

### Theorem 4.2: short-time superlinear and long-time linear convergence

The lower triangular entries of $M_g^{-1}$ are powers of $R_g(z)$, so

$$
M(z)=[R_f^J(z/J)-R_g(z)]\widetilde M(R_g(z)),
$$

where the first subdiagonal of $\widetilde M(\beta)$ is $1$ and deeper diagonals contain $\beta,\beta^2,\ldots$. Bounds on its powers give two regimes. For short horizons,

$$
\max_n\|\boldsymbol e_n^k\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}
\varrho_s(J,z,N_t,k)
\max_n\|\boldsymbol e_n^0\|_\infty,
$$

$$
\varrho_s(J,z,N_t,k)
=\frac{|R_g(z)-R_f^J(z/J)|^k}{k!}
\prod_{j=1}^{k}(N_t-j). \tag{4.5a}
$$

The product vanishes at $k=N_t$. If $|R_g(z)|<1$, a horizon-independent estimate is

$$
\max_n\|\boldsymbol e_n^k\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}\varrho_l^k(J,z)
\max_n\|\boldsymbol e_n^0\|_\infty,
$$

$$
\varrho_l(J,z)=
\frac{|R_g(z)-R_f^J(z/J)|}{1-|R_g(z)|}. \tag{4.5b}
$$

The numerator measures coarse–fine mismatch; the denominator measures the dissipation margin of the coarse method. A nearly nondissipative coarse mode combined with a large mismatch makes the factor approach or exceed one.

![Original Figure 4.2: short-time superlinear and long-time linear convergence](assets/papers/time-parallelization/source-figures/figure-4-2.svg)

Figure 4.2 uses the periodic heat equation, zero source, $u_0(x)=\sin^2(2\pi x)$, $\Delta x=1/5$, backward Euler on both levels, and $J=10$. For $T=0.02,N_t=6$, $\varrho_s$ predicts the superlinear decrease. A longer horizon gives an approximately fixed slope described by $\varrho_l$. Refining to $\Delta x=1/8$ also brings the linear regime forward.

### Theorem 4.3: nonlinear superlinear estimate

Let $\mathcal F$ be exact and $\mathcal G$ an order-$p$ method with local error at most $C_3\Delta T^{p+1}$. Assume

$$
\|\mathcal G(T_n,T_n+\Delta T,\boldsymbol v)
-\mathcal G(T_n,T_n+\Delta T,\boldsymbol w)\|
\le(1+C_2\Delta T)\|\boldsymbol v-\boldsymbol w\|,
$$

and

$$
\mathcal F(T_n,T_{n+1},\boldsymbol v)
-\mathcal G(T_n,T_{n+1},\boldsymbol v)
=c_{p+1}(\boldsymbol v)\Delta T^{p+1}
+c_{p+2}(\boldsymbol v)\Delta T^{p+2}+\cdots,
$$

with continuously differentiable coefficients. Then

$$
\|\boldsymbol u(T_n)-\boldsymbol u_n^k\|
\le
\frac{C_3\Delta T^{p+1}(C_1\Delta T^{p+1})^{k+1}}{(k+1)!}
(1+C_2\Delta T)^{n-k-1}
\prod_{j=0}^{k}(n-j). \tag{4.6}
$$

The product vanishes once $k\ge n$, and every iteration introduces another factor of $\Delta T^{p+1}$ on short intervals.

### Theorem 4.4: the parabolic long-time factor near 0.3

If $\mathcal G$ is backward Euler and $\mathcal F$ an L-stable Runge–Kutta method, then for some $J_{\min}=O(1)$,

$$
\max_{z\in\mathbb R_-}\varrho_l(J,z)\approx0.3,
\qquad J\ge J_{\min}. \tag{4.7}
$$

This factor is independent of $T$ and $N_t$. Suitable Radau IIA combinations can reduce the worst factor to about $0.068$.

For an A-stable but non-L-stable fine method such as the trapezoidal rule, one instead has over a bounded spectrum

$$
\max_{z\in[0,z_{\max}]}\varrho_l(J,z)\approx0.3,
\qquad
J\ge J_{\min}=O(\log_2z_{\max}). \tag{4.8}
$$

More fine steps are needed to resolve the dissipative high-frequency physics. Equation (4.9) gives the two schemes explicitly:

$$
\begin{array}{c|cc}
\gamma&\gamma&0\\
1&1-\gamma&\gamma\\ \hline
&1-\gamma&\gamma
\end{array}
\quad \gamma=\frac{2-\sqrt2}{2}
\qquad\text{(SDIRK22)},
$$

$$
\begin{array}{c|cc}
\gamma&\gamma&0\\
1-\gamma&-1/\sqrt3&\gamma\\ \hline
&1/2&1/2
\end{array}
\quad \gamma=\frac{3+\sqrt3}{6}
\qquad\text{(SDIRK23)}. \tag{4.9}
$$

The bound holds from $J_{\min}=2$ for SDIRK22 and from $J_{\min}=4$ for SDIRK23.

![Original Figure 4.3: Parareal convergence under different fine propagators and coarsening factors](assets/papers/time-parallelization/source-figures/figure-4-3.svg)

Figure 4.3 uses the periodic heat equation with $\Delta x=1/256$, $\Delta T=0.1$, $T=4$, and diffusion coefficient $0.1$. The panels from left to right use $J=2,10,50$. At $J=2$, the trapezoidal curve stalls near $10^{-4}$ and SDIRK23 is visibly slower than SDIRK22. At $J=10$, the two SDIRK curves nearly coincide while the trapezoidal rule retains a slow tail. At $J=50$, all three follow the $0.3^k$ reference slope. The transition verifies the qualification in Theorem 4.4: the fine integrator's stability class and the coarse-to-fine ratio must be assessed together.

### Deterioration as diffusion weakens

The remaining tests fix $T=4$, $\Delta T=0.1$, $\Delta x=1/128$, and $J=32$, with backward Euler coarse propagation and SDIRK22 fine propagation.

![Original Figure 4.4: long-time factors for every advection–diffusion eigenvalue at three viscosities](assets/papers/time-parallelization/source-figures/figure-4-4.svg)

As viscosity falls, the advection–diffusion spectrum spreads from the negative real axis toward the imaginary axis and $\max\varrho_l$ approaches one.

![Original Figure 4.5: deterioration of Parareal on advection–diffusion and Burgers' equation as viscosity falls](assets/papers/time-parallelization/source-figures/figure-4-5.svg)

Figure 4.5(a) confirms the modal prediction. Burgers' equation lacks an equally sharp modal theory, but panel (b) shows the same trend. Standard iteration typically diverges near $\nu\le10^{-3}$. The finite-step property remains algebraically true but loses practical value. Wave equations generally fail as well, marking the intended boundary of the methods in this chapter.

## Equation, theorem, and figure audit

| Source item                          | Paper section | Coverage                                                                 |
| ------------------------------------ | ------------- | ------------------------------------------------------------------------ |
| Section 4 introduction and 4.1       | 4, 4.1        | temporal locality, four method families, historical lines                |
| (4.1), Figure 4.1                    | 4.2           | update, parallel fine solves, sequential coarse correction, two grids    |
| (4.2)–(4.4), Theorem 4.1, Remark 4.1 | 4.2           | modal reduction, Toeplitz matrix, complete proof, preconditioner view    |
| (4.5), Theorem 4.2, Figure 4.2       | 4.2           | superlinear and long-time linear bounds, both regimes                    |
| (4.6), Theorem 4.3                   | 4.2           | nonlinear assumptions, bound, finite-step implication                    |
| (4.7)–(4.9), Theorem 4.4, Figure 4.3 | 4.2           | L/A stability, 0.3 factor, SDIRK comparison                              |
| Figures 4.4–4.5                      | 4.2           | advection–diffusion spectrum, Burgers experiment, weak-diffusion failure |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 4–4.2, pp. 443–452.
