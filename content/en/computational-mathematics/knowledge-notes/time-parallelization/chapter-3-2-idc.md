---
title: "3.3: Parallel Integral Deferred Correction"
description: A step-by-step derivation of IDC and a complete account of PIDC, RIDC, and their regularity requirement
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-3-2-idc
tags:
  - parallel-in-time
  - IDC
  - high-order-time-integration
---

> [!note] Reading scope
> This page covers Section 3.3 of the paper (pp. 405–411), including equations (3.5)–(3.12), Theorem 3.3, and Figures 3.4–3.6. Every algebraic step in the correction formula is retained, and the PIDC and RIDC schedules are explained separately.

## 3.3.1 Starting point of IDC

Write the nonlinear initial-value problem in integral form:

$$
\boldsymbol u(t)=\boldsymbol u_0+
\int_0^t\boldsymbol f(\boldsymbol u(\tau),\tau)\,d\tau,
\qquad t\in(0,T]. \tag{3.5}
$$

Let $\widetilde{\boldsymbol u}(t)$ be a rough approximation, either the constant trajectory $\boldsymbol u_0$ or the output of a low-accuracy integrator. Define the true error and the integral residual by

$$
\boldsymbol e(t)=\boldsymbol u(t)-\widetilde{\boldsymbol u}(t)
$$

and

$$
\boldsymbol r(t)=\boldsymbol u_0+
\int_0^t\boldsymbol f(\widetilde{\boldsymbol u}(\tau),\tau)\,d\tau
-\widetilde{\boldsymbol u}(t). \tag{3.6}
$$

The residual measures the defect of the current trajectory in the integral equation. Substitute $\boldsymbol u=\widetilde{\boldsymbol u}+\boldsymbol e$ into (3.5), then add and subtract $\int_0^t\boldsymbol f(\widetilde{\boldsymbol u},\tau)d\tau$:

$$
\begin{aligned}
\boldsymbol e(t)
&=\boldsymbol u_0+
\int_0^t\boldsymbol f(\widetilde{\boldsymbol u}(\tau)+\boldsymbol e(\tau),\tau)\,d\tau
-\widetilde{\boldsymbol u}(t)\\
&=\boldsymbol r(t)+
\int_0^t\left[
\boldsymbol f(\widetilde{\boldsymbol u}(\tau)+\boldsymbol e(\tau),\tau)
-\boldsymbol f(\widetilde{\boldsymbol u}(\tau),\tau)
\right]d\tau. \tag{3.7}
\end{aligned}
$$

Differentiating gives the error equation

$$
\boldsymbol e'(t)-\boldsymbol r'(t)
=\boldsymbol f(\widetilde{\boldsymbol u}(t)+\boldsymbol e(t),t)
-\boldsymbol f(\widetilde{\boldsymbol u}(t),t). \tag{3.8}
$$

The task of improving the solution has become an error solve. The residual integral carries high-order information, while the differential error equation may still use a simple stepper.

## 3.3.2 Discrete correction formula

Choose nodes

$$
0=t_0<t_1<\cdots<t_M=T,
\qquad \boldsymbol u_m^k\approx\boldsymbol u(t_m),
$$

and let $\Delta t_m=t_{m+1}-t_m$. Applying the linear $\theta$ method to (3.8) gives

$$
\begin{aligned}
\boldsymbol e_{m+1}-\boldsymbol e_m
={}&\boldsymbol r_{m+1}-\boldsymbol r_m\\
&+\Delta t_m(1-\theta)
\left[\boldsymbol f(\boldsymbol u_m^{k+1},t_m)
-\boldsymbol f(\boldsymbol u_m^k,t_m)\right]\\
&+\Delta t_m\theta
\left[\boldsymbol f(\boldsymbol u_{m+1}^{k+1},t_{m+1})
-\boldsymbol f(\boldsymbol u_{m+1}^k,t_{m+1})\right]. \tag{3.9}
\end{aligned}
$$

The residual definition gives

$$
\boldsymbol r_{m+1}-\boldsymbol r_m
=\int_{t_m}^{t_{m+1}}\boldsymbol f(\boldsymbol u^k(\tau),\tau)\,d\tau
-\left(\boldsymbol u_{m+1}^k-\boldsymbol u_m^k\right).
$$

Lagrange interpolation over the complete node set approximates the local integral:

$$
\int_{t_m}^{t_{m+1}}\boldsymbol f(\boldsymbol u^k(\tau),\tau)\,d\tau
\approx\sum_{j=1}^{M}\omega_{m,j}
\boldsymbol f(\boldsymbol u_j^k,t_j), \tag{3.10a}
$$

$$
\omega_{m,j}=
\int_{t_m}^{t_{m+1}}
\prod_{\substack{i=1\\i\ne j}}^{M}
\frac{\tau-t_i}{t_j-t_i}\,d\tau. \tag{3.10b}
$$

Finally use $\boldsymbol u_m^{k+1}=\boldsymbol u_m^k+\boldsymbol e_m$ to eliminate the error variable. The node update is

$$
\begin{aligned}
\boldsymbol u_{m+1}^{k+1}
={}&\boldsymbol u_m^{k+1}
+\Delta t_m(1-\theta)
\left[\boldsymbol f(\boldsymbol u_m^{k+1},t_m)
-\boldsymbol f(\boldsymbol u_m^k,t_m)\right]\\
&+\Delta t_m\theta
\left[\boldsymbol f(\boldsymbol u_{m+1}^{k+1},t_{m+1})
-\boldsymbol f(\boldsymbol u_{m+1}^k,t_{m+1})\right]\\
&+\sum_{j=1}^{M}\omega_{m,j}
\boldsymbol f(\boldsymbol u_j^k,t_j),
\qquad m=0,\ldots,M-1. \tag{3.11}
\end{aligned}
$$

Every correction level $k$ sweeps from left to right through $m=0,1,\ldots,M-1$. The choice $\theta=1$ gives a backward-Euler correction, and $\theta=1/2$ gives a trapezoidal correction.

> [!tip] Roles of the terms in (3.11)
> The first line advances the new trajectory with a simple integrator. The middle terms correct the difference between the new and old dynamics at the local endpoints. The final quadrature term injects high-order integral information from the complete old node set. The last term is what lets the correction exceed the order of the base stepper.

## 3.3.3 Theorem 3.3: order gained per correction

If the base integrator has order $p$ and $M$ equally spaced nodes are used, the approximation after correction $k$ has order

$$
O\!\left(\Delta t^{\min\{M,(k+1)p\}}\right).
$$

Each correction therefore adds at most $p$ orders until the quadrature ceiling $M$ is reached. Backward Euler has $p=1$, and the trapezoidal rule has $p=2$. Gauss–Lobatto nodes can reach order $2J-1$. This form of IDC is spectral deferred correction (SDC), the central integration component of PFASST.

A single high-order polynomial becomes ineffective on a long interval. The paper partitions $[0,T]$ into windows

$$
I_n=[T_{n-1},T_n],
\qquad n=1,\ldots,N_t.
$$

Standard IDC completes every correction on $I_n$ before passing its endpoint to $I_{n+1}$. Windows are sequential, and the node updates within a window are also sequential in $m$.

## 3.3.4 PIDC: a pipeline across windows

PIDC applies the pipeline idea of Womble (1990). After the first sweep on $I_n$, the rough endpoint $\boldsymbol u_{n,M}^{1}$ is available. Window $I_{n+1}$ immediately starts its first sweep from this value while $I_n$ performs its second. In general, when $I_n$ executes sweep $k$, $I_{n+1}$ can execute sweep $k-1$, $I_{n+2}$ sweep $k-2$, and so on through the first sweep on $I_{n+k-1}$.

![Source Figure 3.4: PIDC pipeline startup and steady state on four time windows](assets/papers/time-parallelization/source-figures/figure-3-4.svg)

Figure 3.4 uses $M=6$ and $k_{\max}=4$. The first four stages fill the pipeline. Four sweeps on four windows then run concurrently. Black dashed lines record completed sweep histories, red lines with circles mark the sweeps currently running, and the solid black lines are the exact solution.

Every PIDC window receives a rough initial value that changes as the upstream window is corrected. Additional corrections therefore need not reduce the error monotonically. The regularity experiments below quantify the issue.

## 3.3.5 Periodic advection–diffusion experiment and matrix (3.12)

The paper uses periodic boundaries and $\Delta x=1/64$. The semidiscrete matrix is

$$
A=\frac{\nu}{\Delta x^2}A_{xx}
-\frac{1}{2\Delta x}A_x,
$$

where

$$
A_{xx}=
\begin{bmatrix}
-2&1&&&1\\
1&-2&1&&\\
&\ddots&\ddots&\ddots&\\
&&1&-2&1\\
1&&&1&-2
\end{bmatrix},
\qquad
A_x=
\begin{bmatrix}
0&1&&&-1\\
-1&0&1&&\\
&\ddots&\ddots&\ddots&\\
&&-1&0&1\\
1&&&-1&0
\end{bmatrix}. \tag{3.12}
$$

The remaining parameters are $T=3$, window width $\Delta T=0.1$, $M=5$, and backward Euler as the base stepper. Error on window $n$ after sweep $k$ is

$$
\operatorname{err}_n^k=
\frac{\max_m\lVert\boldsymbol u_{\mathrm{ref}}^{n,m}
-\boldsymbol u_k^{n,m}\rVert_\infty}
{\max_{n,m}\lVert\boldsymbol u_{\mathrm{ref}}^{n,m}\rVert_\infty}.
$$

MATLAB ODE45 supplies the reference trajectory with relative and absolute tolerances both set to $10^{-13}$. The initial trajectory on the next window is held constant at the first-sweep endpoint of the previous window: $\boldsymbol u_{n+1,m}^0\equiv\boldsymbol u_{n,M}^1$.

![Source Figure 3.5: windowwise errors of IDC and PIDC under two regularity and viscosity settings](assets/papers/time-parallelization/source-figures/figure-3-5.svg)

The four panels in Figure 3.5 form two controlled comparisons.

- Panels (a) and (b) use the sharp source with $\sigma=1000$; (c) and (d) use the smooth source with $\sigma=5$.
- Panels (a) and (c) use $\nu=1$; (b) and (d) use $\nu=10^{-3}$.
- Every panel reports the initial error and the IDC/PIDC errors after one and two sweeps.

With low regularity, the second sweep does not reduce the error further. With a smooth source and large viscosity, both IDC and PIDC improve again on the second sweep and remain comparable. Small viscosity weakens the order increase; in panel (d), the second PIDC sweep is clearly worse than sequential IDC. The paper concludes that PIDC is poorly suited to low-regularity hyperbolic solutions because it relies on high-order correction.

## 3.3.6 RIDC: a sliding quadrature window

RIDC moves the concurrency from time windows to correction levels. The first processor advances continuously with a low-order integrator. Once it has produced the first $M$ values, a second processor starts the first IDC correction. After reaching step $M$, the second processor keeps advancing and slides its quadrature nodes from $1,\ldots,M$ to $2,\ldots,M+1$, then to $3,\ldots,M+2$. A third processor starts when sufficient first-correction data is available and follows the same sliding rule.

In steady state, each processor owns one correction level and all levels advance adjacent time steps concurrently. The pipeline must retain several historical values and handle both startup and drain phases.

![Source Figure 3.6: windowwise errors of IDC and RIDC on the same advection–diffusion tests](assets/papers/time-parallelization/source-figures/figure-3-6.svg)

Figure 3.6 reuses the PDE, grid, source, and viscosity settings of Figure 3.5, replacing PIDC by RIDC. The conclusion is unchanged. A smooth solution permits order growth across correction levels. A sharp source and small viscosity remove much of the expected benefit of high-order interpolation and quadrature. RIDC changes scheduling and memory organization; it retains IDC's regularity requirement.

## Equation and figure coverage

| Source item        | Location here | Coverage                                                                              |
| ------------------ | ------------- | ------------------------------------------------------------------------------------- |
| (3.5)–(3.8)        | 3.3.1         | integral equation, residual, integral error equation, and differential error equation |
| (3.9)–(3.11)       | 3.3.2         | $\theta$ discretization, quadrature weights, and final IDC update                     |
| Theorem 3.3        | 3.3.3         | base order, correction count, and quadrature ceiling                                  |
| Figure 3.4         | 3.3.4         | complete PIDC startup and steady-state schedule                                       |
| (3.12), Figure 3.5 | 3.3.5         | periodic matrices, all parameters, error definition, and four-panel interpretation    |
| Figure 3.6         | 3.3.6         | RIDC sliding window and regularity experiment                                         |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3.3, pp. 405–411.
