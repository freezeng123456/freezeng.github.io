---
title: "3.4: ParaExp"
description: A complete derivation of exact linear superposition, nonlinear iteration, and the connection to Parareal
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-3-3-paraexp
tags:
  - parallel-in-time
  - ParaExp
  - matrix-exponential
---

> [!note] Reading scope
> This page covers Section 3.4 of the paper (pp. 412–415), including equations (3.13)–(3.21), Theorem 3.4, and Figures 3.7–3.8. The induction proof of the linear reconstruction and the finite-step property of the nonlinear method are both expanded.

## 3.4 ParaExp

### Two families of parallel subproblems

ParaExp solves

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol g(t),
\qquad \boldsymbol u(0)=\boldsymbol u_0.
$$

Partition time into $[T_{n-1},T_n]$, $n=1,\ldots,N_t$. The first, “red,” family keeps the source and uses a zero initial value on every local interval:

$$
\boldsymbol v_n'(t)=A\boldsymbol v_n(t)+\boldsymbol g(t),
\quad t\in(T_{n-1},T_n],
\qquad
\boldsymbol v_n(T_{n-1})=0. \tag{3.13}
$$

All red problems are independent. The second, “blue,” family removes the source and propagates the red endpoint from the previous interval to the global final time:

$$
\boldsymbol w_n'(t)=A\boldsymbol w_n(t),
\quad t\in(T_{n-1},T],
\qquad
\boldsymbol w_n(T_{n-1})=\boldsymbol v_{n-1}(T_{n-1}), \tag{3.14}
$$

where $\boldsymbol v_0(T_0)=\boldsymbol u_0$. The blue problems are also mutually independent.

![Source Figure 3.7: local red forced problems and homogeneous blue tail propagation in ParaExp](assets/papers/time-parallelization/source-figures/figure-3-7.svg)

The red curves in Figure 3.7 remain local, while every blue dashed curve extends from one interface to $T$. On interval $n$, the exact solution is the local red response plus all blue responses that have already started:

$$
\boldsymbol u(t)=\boldsymbol v_n(t)
+\sum_{j=1}^{n}\boldsymbol w_j(t),
\qquad t\in[T_{n-1},T_n]. \tag{3.15}
$$

![ParaExp separates local forced responses from global homogeneous propagation](assets/diagrams/pint/en/paraexp-decomposition.svg)

### Induction proof of (3.15)

On the first interval, add (3.13) and (3.14):

$$
(\boldsymbol v_1+\boldsymbol w_1)'
=A(\boldsymbol v_1+\boldsymbol w_1)+\boldsymbol g,
$$

with

$$
\boldsymbol v_1(0)+\boldsymbol w_1(0)
=0+\boldsymbol u_0.
$$

Uniqueness gives $\boldsymbol u=\boldsymbol v_1+\boldsymbol w_1$, proving (3.15) for $n=1$.

Assume the formula holds on interval $n$. At the right endpoint,

$$
\boldsymbol u(T_n)=\boldsymbol v_n(T_n)
+\sum_{j=1}^{n}\boldsymbol w_j(T_n).
$$

Equation (3.14) gives $\boldsymbol w_{n+1}(T_n)=\boldsymbol v_n(T_n)$, hence

$$
\boldsymbol u(T_n)=\sum_{j=1}^{n+1}\boldsymbol w_j(T_n).
$$

Let $\boldsymbol w=\sum_{j=1}^{n+1}\boldsymbol w_j$. On $(T_n,T_{n+1}]$, it satisfies $\boldsymbol w'=A\boldsymbol w$ with initial value $\boldsymbol u(T_n)$. Adding the zero-initial forced solution $\boldsymbol v_{n+1}$ gives the unique solution of the original initial-value problem on the next interval. The induction is complete.

### Why the long blue tails can remain inexpensive

The homogeneous solution is

$$
\boldsymbol w_n(t)=
\exp\!\left((t-T_{n-1})A\right)
\boldsymbol v_{n-1}(T_{n-1}),
\qquad t\in[T_{n-1},T]. \tag{3.16}
$$

An algorithm for $e^{\tau A}\boldsymbol b$ can jump directly to the requested time; its cost need not scale with the number of intermediate steps. The paper lists rational Krylov and Chebyshev expansions for large sparse matrices, and scaling-and-squaring with Padé approximation for smaller matrices. MATLAB R2023b and later provide `expmv` for a matrix-exponential action. REXI and early Laplace-transform PinT methods use related exponential-approximation ideas.

The wave-equation experiment of Gander and Güttel (2013) reported time-parallel efficiency up to about 80%. This number belongs to that implementation and depends on the exponential algorithm, matrix structure, partition, and hardware.

### Nonlinear splitting

Assume the nonlinear dynamics can be written as

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B(\boldsymbol u(t))+\boldsymbol g(t). \tag{3.17}
$$

Direct superposition no longer applies because $B$ recouples the intervals. The initial construction writes $\boldsymbol u=\boldsymbol w+\boldsymbol v$, evolves $\boldsymbol w'=A\boldsymbol w$, and evaluates $B(\boldsymbol v+\boldsymbol w)$ in the $\boldsymbol v$ equation. Time parallelism then requires an iteration whose blue problems start from the previous interface iterate.

Explicitly evaluating the full blue sum $\sum_j\boldsymbol w_j^k(t)$ inside every nonlinear local solve would repeat expensive work when $A$ is large. The paper substitutes $\boldsymbol v_n^k=\boldsymbol u_n^k-\sum_{j=1}^n\boldsymbol w_j^k$ and obtains the following two-stage iteration.

First construct the homogeneous propagations for $n=1,\ldots,N_t$:

$$
\begin{aligned}
(\boldsymbol w_n^k)'(t)&=A\boldsymbol w_n^k(t),
&&t\in[T_{n-1},T],\\
\boldsymbol w_n^k(T_{n-1})
&=\boldsymbol u_{n-1}^{k-1}(T_{n-1})
-\sum_{j=1}^{n-1}\boldsymbol w_j^{k-1}(T_{n-1}),
&&\boldsymbol w_1^k(T_0)=\boldsymbol u_0.
\end{aligned} \tag{3.18}
$$

Then solve the complete nonlinear problems concurrently on all local intervals:

$$
\begin{aligned}
(\boldsymbol u_n^k)'(t)
&=A\boldsymbol u_n^k(t)+B(\boldsymbol u_n^k(t))+\boldsymbol g(t),
&&t\in[T_{n-1},T_n],\\
\boldsymbol u_n^k(T_{n-1})
&=\sum_{j=1}^{n}\boldsymbol w_j^k(T_{n-1}).
\end{aligned} \tag{3.19}
$$

The global iterate equals $\boldsymbol u_n^k(t)$ on interval $n$.

### Theorem 3.4: finite-step convergence and Parareal equivalence

**Finite-step result.** After iteration $k$, $\boldsymbol u^k(t)$ agrees with the exact solution on $[0,T_k]$. A time-interval induction explains the result. The first interval always starts from the true initial value. If the first $k-1$ intervals are exact, (3.18) constructs the exact initial value at the next interface, and (3.19) advances exactness by one more interval.

At the coarse points, the method is equivalent to

$$
\boldsymbol U_n^k
=\mathcal G(T_{n-1},T_n,\boldsymbol U_{n-1}^k)
+\mathcal F(T_{n-1},T_n,\boldsymbol U_{n-1}^{k-1})
-\mathcal G(T_{n-1},T_n,\boldsymbol U_{n-1}^{k-1}), \tag{3.20a}
$$

where the coarse propagator solves only the homogeneous linear problem

$$
\boldsymbol u'=A\boldsymbol u,
\qquad
\boldsymbol u(T_{n-1})=\boldsymbol U,
\qquad t\in[T_{n-1},T_n], \tag{3.20b}
$$

and the fine propagator solves the complete nonlinear problem

$$
\boldsymbol u'=A\boldsymbol u+B(\boldsymbol u)+\boldsymbol g,
\qquad
\boldsymbol u(T_{n-1})=\boldsymbol U,
\qquad t\in[T_{n-1},T_n]. \tag{3.20c}
$$

Standard Parareal normally lets its coarse propagator approximate the complete nonlinear problem as well. Here $\mathcal G$ retains only $A$, so this is a simplified Parareal. Nonlinear ParaExp succeeds when that split captures the dominant dynamics.

### Figure 3.8: splitting failure on Burgers' equation

The experiment uses

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B\boldsymbol u^2(t),
\qquad t\in(0,2), \tag{3.21}
$$

obtained from centered differences for periodic Burgers' equation. The mesh is $\Delta x=1/100$, with $A=\nu A_{xx}/\Delta x^2$ and $B=-A_x/(2\Delta x)$; the matrices $A_{xx}$ and $A_x$ are those in (3.12). The fine propagator in both methods is backward Euler with step $0.01/20$. Standard Parareal also uses backward Euler as its coarse propagator, with step $0.01$. ParaExp uses MATLAB `expmv` for the linear coarse propagation.

![Source Figure 3.8: errors of nonlinear ParaExp and standard Parareal at three viscosities](assets/papers/time-parallelization/source-figures/figure-3-8.svg)

The panels from left to right use $\nu=1,0.1,0.02$. Their horizontal lines mark the truncation level $\max\{\Delta t,\Delta x^2\}$, beyond which further iteration is unnecessary in practice. At $\nu=1$, $A$ captures the dominant dynamics and ParaExp is substantially faster than standard Parareal. At $\nu=0.1$, standard Parareal becomes faster while ParaExp still decays slowly. At $\nu=0.02$, the ParaExp error grows at every displayed step, whereas standard Parareal still crosses the truncation line. Standard Parareal eventually fails as viscosity is reduced further. The three panels show a transfer of dominant dynamics between the two parts of the split, not merely a uniform slowdown under weaker diffusion.

> [!important] Applicability boundary
> Linear ParaExp equation (3.15) is an exact algebraic decomposition, subject only to errors in the exponential action and local forced solves. Nonlinear ParaExp is iterative and depends on the $A+B$ split. An exact exponential cannot repair a linear part that misses the dominant propagation mechanism.

## Equation and figure coverage

| Source item                  | Paper section | Coverage                                                                                                        |
| ---------------------------- | ------------- | --------------------------------------------------------------------------------------------------------------- |
| (3.13)–(3.16), Figure 3.7    | 3.4           | both subproblem families, exact reconstruction, induction proof, exponential action, and complete source figure |
| (3.17)–(3.19)                | 3.4           | nonlinear split, source of redundant work, and reformulated iteration                                           |
| (3.20a)–(3.20c), Theorem 3.4 | 3.4           | finite-step result, induction argument, and Parareal propagators                                                |
| (3.21), Figure 3.8           | 3.4           | Burgers discretization, all temporal parameters, truncation line, and three-regime interpretation               |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3.4, pp. 412–415.
