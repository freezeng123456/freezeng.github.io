---
title: "4.5: Diagonalization-Based Parareal"
description: The two complete routes of parallel coarse-grid correction and an interval-local diagonalized coarse propagator
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-3-diagonalized-parareal
tags:
  - parallel-in-time
  - Parareal
  - ParaDiag
---

> [!note] Reading scope
> This page follows Section 4.5 (pp. 461–472). It covers equations (4.14)–(4.29), Theorems 4.7–4.8, Remark 4.2, and Figures 4.12–4.17. Both variants use diagonalization, but at different locations: the first parallelizes the global coarse-grid correction; the second defines a special coarse propagator inside each coarse interval.

## 4.5.1 Distinguishing the two constructions

- **Diagonalized CGC (Section 4.5.1):** modifies the serial correction across the $N_t$ coarse points. Its concurrency lies across coarse points, and its convergence mechanism remains close to standard Parareal, so its main range is parabolic.
- **Diagonalized coarse propagator (Section 4.5.2):** retains the outer Parareal form and uses ParaDiag on the $J$ fine points inside every $[T_n,T_{n+1}]$. Coarse and fine propagation use the same integrator and step size. This construction transports long-lived modes and can also handle hyperbolic problems.

## 4.5.1.1 From serial CGC to head–tail coupling

Standard coarse-grid correction is

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\boldsymbol b_{n+1}^k,
\qquad n=0,\ldots,N_t-1, \tag{4.14}
$$

$$
\boldsymbol b_{n+1}^k
=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k),
$$

starting from $\boldsymbol u_0^{k+1}=\boldsymbol u_0$. Wu (2018) replaced this initial condition by

$$
\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0.
$$

To preserve the original fixed point, define

$$
\widetilde{\boldsymbol u}_n^k=
\begin{cases}
\boldsymbol u_0,&n=0,\\
\boldsymbol u_n^k,&n\ge1,
\end{cases}
$$

and iterate

$$
\left\{
\begin{aligned}
\boldsymbol u_{n+1}^{k+1}
={}&\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol u}_n^k)
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k),\\
\boldsymbol u_0^{k+1}&=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0.
\end{aligned}
\right. \tag{4.15}
$$

This predates the more natural ParaDiag-II condition (3.55); the modified value $\widetilde{\boldsymbol u}_0^k$ maintains consistency at convergence.

## 4.5.1.2 Linear all-at-once system and three-stage solve

For $\boldsymbol u'=A\boldsymbol u$ and backward Euler coarse propagation, substitution of the head–tail condition yields

$$
(C_\alpha\otimes I_x-I_t\otimes\Delta TA)\boldsymbol U^{k+1}
=\boldsymbol g^k, \tag{4.16}
$$

where

$$
C_\alpha=
\begin{bmatrix}
1&&&-\alpha\\-1&1\\&\ddots&\ddots\\&&-1&1
\end{bmatrix},
$$

$$
\boldsymbol g^k=
\begin{bmatrix}
\boldsymbol u_0+(I_x-\Delta TA)\boldsymbol b_1^k\\
(I_x-\Delta TA)\boldsymbol b_2^k\\\vdots\\
(I_x-\Delta TA)\boldsymbol b_{N_t}^k
\end{bmatrix}.
$$

The paper writes the diagonalized solve as

$$
\left\{
\begin{aligned}
\boldsymbol U^{a,k+1}&=(F\otimes I_x)\boldsymbol g^k,\\
(\lambda_nI_x-\Delta TA)\boldsymbol u_n^{b,k+1}
&=\boldsymbol u_n^{a,k+1},
&&n=1,\ldots,N_t,\\
\boldsymbol U^{k+1}&=(F^*\otimes I_x)\boldsymbol U^{b,k+1}.
\end{aligned}
\right. \tag{4.17}
$$

An alpha-circulant implementation also applies the corresponding diagonal scalings. As in Section 3.5.2, the essential stages are an FFT-like transform, independent shifted spatial solves, and the inverse transform.

## 4.5.1.3 Theorem 4.7: threshold for matching standard Parareal

As $\alpha\to0$, equation (4.15) approaches standard CGC, while alpha-circulant roundoff increases. Let $\rho$ be the standard Parareal factor and $\rho_{\mathrm{new}}$ the new factor. For stable coarse propagation and a linear system with negative real eigenvalues,

$$
\rho_{\mathrm{new}}=\rho,
\qquad
\alpha\le\frac{\rho}{1+\rho}. \tag{Theorem 4.7}
$$

The practical choice is the threshold itself: reducing alpha further does not improve the asymptotic rate and increases roundoff exposure. Typically both $\rho$ and alpha are $O(10^{-1})$.

![Original Figure 4.12: standard and diagonalized CGC for heat and ADE](assets/papers/time-parallelization/source-figures/figure-4-12.svg)

The test uses periodic data, $u_0(x)=\sin(2\pi x)$, backward Euler coarse and SDIRK22 fine propagation, $T=4$, $J=10$, $\Delta T=0.1$, and $\Delta x=1/128$. Heat gives $\rho\approx0.22$ and threshold $0.18$; ADE at $\nu=0.1$ gives $\rho\approx0.39$ and threshold $0.28$. Alpha above the threshold slows diagonalized CGC.

## 4.5.1.4 Nonlinear all-at-once quasi-Newton solve

With backward Euler coarse propagation, define

$$
\boldsymbol b_{n+1}^k
=\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol u}_n^k)
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k).
$$

The correction is

$$
(C_\alpha\otimes I_x)\boldsymbol U^{k+1}
-\Delta TF(\boldsymbol U^{k+1})=\boldsymbol g^k. \tag{4.18}
$$

The $n$th block of $F$ is $f(\boldsymbol u_n^{k+1}-\boldsymbol b_n^k)$. The inner quasi-Newton iteration is

$$
P_\alpha^{k+1,l}\Delta\boldsymbol U^{k+1,l}
=\boldsymbol g^k-(C_\alpha\otimes I_x)\boldsymbol U^{k+1,l}
+\Delta TF(\boldsymbol U^{k+1,l}),
$$

$$
\boldsymbol U^{k+1,l+1}
=\boldsymbol U^{k+1,l}+\Delta\boldsymbol U^{k+1,l}, \tag{4.19a}
$$

$$
P_\alpha^{k+1,l}
=C_\alpha\otimes I_x-I_t\otimes\Delta TA^{k+1,l}, \tag{4.19b}
$$

where $A^{k+1,l}$ averages the temporal Jacobian blocks. The matrix has the structure of (4.16), so every increment uses (4.17). A nearest Kronecker approximation could replace the average Jacobian but is not pursued in the paper.

![Original Figure 4.13: the two CGCs for Burgers' equation at two viscosities](assets/papers/time-parallelization/source-figures/figure-4-13.svg)

The nonlinear experiment displays the same alpha threshold: a suitably small alpha makes diagonalized CGC track standard CGC.

### Remark 4.2: MGRIT needs a consistent head–tail condition

Directly transplanting (4.15) into MGRIT diverges for every alpha. The consistent condition is

$$
\boldsymbol u_1^{k+1}
=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_1,
$$

which leads to

$$
\left\{
\begin{aligned}
\boldsymbol u_0^{k+1}&=\boldsymbol u_0,\\
\boldsymbol u_1^{k+1}
&=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_1,\\
\boldsymbol u_{n+1}^{k+1}
&=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\widetilde{\boldsymbol b}_{n+1}^k.
\end{aligned}
\right. \tag{4.20}
$$

Here $\widetilde{\boldsymbol b}_{n+1}^k=\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol s}_n^k)-\mathcal G(T_n,T_{n+1},\widetilde{\boldsymbol s}_n^k)$ and $\widetilde{\boldsymbol s}_n^k=\mathcal F(T_{n-1},T_n,\widetilde{\boldsymbol u}_{n-1}^k)$. For small alpha this variant matches original MGRIT, with the same threshold mechanism as Theorem 4.7.

## 4.5.2.1 Fine propagation and a head–tail coarse propagator

Inside each coarse interval, both propagators use the same linear-theta method and $\Delta t=\Delta T/J$. Fine propagation advances sequentially:

$$
\boldsymbol v_{j+1}-\boldsymbol v_j
=\Delta t[\theta f(\boldsymbol v_{j+1})
+(1-\theta)f(\boldsymbol v_j)],
\quad j=0,\ldots,J-1,
\quad \boldsymbol v_0=\boldsymbol u_n. \tag{4.21}
$$

$\theta=1$ is backward Euler and $\theta=1/2$ is trapezoidal. The special coarse propagator $\mathcal F_\alpha^*$ changes only the condition to

$$
\boldsymbol v_0=\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n. \tag{4.22}
$$

The $J$ fine steps are now a head–tail system solvable in parallel by ParaDiag.

## 4.5.2.2 Nonlinear all-at-once system and quasi-Newton iteration

With $\boldsymbol V=(\boldsymbol v_1^\top,\ldots,\boldsymbol v_J^\top)^\top$,

$$
\underbrace{(C_\alpha\otimes I_x)\boldsymbol V
-\Delta tF(\boldsymbol V)}_{K(\boldsymbol V)}
=\boldsymbol b(\boldsymbol u_n), \tag{4.23}
$$

$$
\boldsymbol b(\boldsymbol u_n)
=((1-\alpha)\boldsymbol u_n^\top,0,\ldots,0)^\top. \tag{4.24}
$$

The first block of $F$ contains $\theta f(\boldsymbol v_1)$ and $(1-\theta)f(\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n)$; later blocks are $\theta f(\boldsymbol v_j)+(1-\theta)f(\boldsymbol v_{j-1})$. The quasi-Newton update is

$$
P_\alpha(\boldsymbol V^l)\Delta\boldsymbol V^l
=\boldsymbol b(\boldsymbol u_n)-K(\boldsymbol V^l),
\qquad
\boldsymbol V^{l+1}=\boldsymbol V^l+\Delta\boldsymbol V^l, \tag{4.25a}
$$

$$
P_\alpha(\boldsymbol V^l)
=C_\alpha\otimes I_x
-\Delta t\widetilde C_{\alpha,\theta}\otimes\overline{\nabla f}(\boldsymbol V^l), \tag{4.25b}
$$

where

$$
\widetilde C_{\alpha,\theta}=
\begin{bmatrix}
\theta&&&(1-\theta)\alpha\\
1-\theta&\theta\\&\ddots&\ddots\\&&1-\theta&\theta
\end{bmatrix},
$$

and $\overline{\nabla f}$ averages the $J$ Jacobian blocks. The two temporal matrices are simultaneously diagonalizable. The outer iteration is

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal F_\alpha^*(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
-\mathcal F_\alpha^*(T_n,T_{n+1},\boldsymbol u_n^k). \tag{4.26}
$$

## 4.5.2.3 Linear system, concurrency, and limiting cases

For $f(\boldsymbol u)=A\boldsymbol u$,

$$
(C_\alpha\otimes I_x
-\widetilde C_{\theta,\alpha}\otimes\Delta tA)\boldsymbol V
=\boldsymbol b(\boldsymbol u_n), \tag{4.27}
$$

$$
\boldsymbol b(\boldsymbol u_n)
=([(I_x+\Delta t(1-\theta)A)(1-\alpha)\boldsymbol u_n]^\top,0,\ldots,0)^\top.
$$

The coarse output is $\mathcal F_\alpha^*=(H_J\otimes I_x)\boldsymbol V=\boldsymbol v_J$, equivalently

$$
\left\{
\begin{aligned}
\boldsymbol v_{j+1}-\boldsymbol v_j
&=\Delta tA[\theta\boldsymbol v_{j+1}+(1-\theta)\boldsymbol v_j],\\
\boldsymbol v_0&=\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n^k.
\end{aligned}
\right. \tag{4.28}
$$

At $\alpha=0$, coarse propagation equals sequential fine propagation and outer Parareal converges in one iteration with no speedup. For $0<\alpha<1$, all $J$ points are solved at once; with enough spatial-solve resources, wall time is approximately $1/J$ of sequential fine propagation.

## 4.5.2.4 Theorem 4.8: parabolic and hyperbolic spectra

For a stable one-step Runge–Kutta method, let

$$
e^k=\max_{1\le n\le N_t}\|\boldsymbol u_n-\boldsymbol u_n^k\|_\infty.
$$

Then

$$
e^k\le\rho^ke^0,
\qquad
\rho=
\begin{cases}
\alpha,&\sigma(A)\subset\mathbb R_-,\\[4pt]
\dfrac{2\alpha N_t}{1+\alpha},&\sigma(A)\subset i\mathbb R.
\end{cases} \tag{4.29}
$$

The heat-equation factor is independent of the number of coarse intervals. The imaginary-spectrum bound grows linearly with $N_t$ and may be loose when $\alpha N_t$ is large.

![Original Figure 4.14: sharp rho=alpha prediction for the heat equation](assets/papers/time-parallelization/source-figures/figure-4-14.svg)

The heat test uses homogeneous Dirichlet data, $u_0=\sin^2(2\pi x)$, trapezoidal integration, $\Delta T=1/2$, $J=10$, and $\Delta x=1/100$. Both choices of $N_t$ follow the alpha slope.

![Original Figure 4.15: joint influence of alpha and the coarse-interval count on the wave equation](assets/papers/time-parallelization/source-figures/figure-4-15.svg)

The wave test uses periodic data, $u_0=\sin^2(2\pi x)$, and $u_t(0)=0$. At $\alpha=0.01$, increasing $N_t$ slows convergence. At smaller alpha, increasing $N_t$ from 24 to 960 costs only about two extra iterations to reach $\max\{\Delta t^2,\Delta x^2\}$.

![Original Figure 4.16: the linear bound is sharp for small alpha Nt and conservative when the product is large](assets/papers/time-parallelization/source-figures/figure-4-16.svg)

At $\alpha=10^{-4},N_t=24$, (4.29) is accurate. Larger-product cases display superlinear decay.

![Original Figure 4.17: iterations needed to reach 1e-8 on Burgers' equation](assets/papers/time-parallelization/source-figures/figure-4-17.svg)

The Burgers test uses periodic data, $u_0=\sin^2(2\pi x)$, $\Delta T=0.1$, $J=10$, and $\Delta x=1/100$. At $N_t=40$, small alpha accelerates convergence and reduces viscosity sensitivity. At $\alpha=10^{-3}$, iteration counts remain robust from $N_t=10$ to $160$. Nonlinear theory gives $\rho=O(\alpha)$ under exact solution of (4.23) and suitable Lipschitz assumptions.

## Final comparison

| Property                  | Diagonalized CGC (4.15)    | Diagonalized coarse propagator (4.26)                                |
| ------------------------- | -------------------------- | -------------------------------------------------------------------- |
| diagonalization direction | global $N_t$ coarse points | $J$ fine points within each coarse interval                          |
| modified component        | CGC                        | coarse propagator                                                    |
| integrators               | coarse and fine may differ | identical integrator and step size                                   |
| main range                | parabolic problems         | parabolic and hyperbolic problems                                    |
| alpha rule                | $\alpha\le\rho/(1+\rho)$   | parabolic factor $\alpha$; hyperbolic bound $2\alpha N_t/(1+\alpha)$ |

## Equation, theorem, and figure audit

| Source item                            | Location here | Coverage                                                                 |
| -------------------------------------- | ------------- | ------------------------------------------------------------------------ |
| (4.14)–(4.17)                          | §§4.5.1.1–2   | standard/head–tail CGC, linear all-at-once matrix, three stages          |
| Theorem 4.7, Figure 4.12               | §4.5.1.3      | alpha threshold, roundoff tradeoff, heat and ADE                         |
| (4.18)–(4.19), Figure 4.13             | §4.5.1.4      | nonlinear system, average-Jacobian quasi-Newton, Burgers                 |
| Remark 4.2, (4.20)                     | §4.5.1.4      | consistent MGRIT head–tail condition and convergent variant              |
| (4.21)–(4.26)                          | §§4.5.2.1–2   | equal-integrator fine/coarse propagation, nonlinear system, outer update |
| (4.27)–(4.28)                          | §4.5.2.3      | linear form, terminal extraction, alpha-zero limit, J-way concurrency    |
| Theorem 4.8, (4.29), Figures 4.14–4.17 | §4.5.2.4      | negative-real/imaginary bounds and all heat, wave, Burgers figures       |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4.5, pp. 461–472.
