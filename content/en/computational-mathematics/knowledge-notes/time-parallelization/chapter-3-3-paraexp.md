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

An algorithm for $e^{\tau A}\boldsymbol b$ can jump directly to the
requested time; its cost need not scale with the number of intermediate
steps. Rational Krylov and Chebyshev expansions suit large sparse
matrices, while scaling-and-squaring with Padé approximation suits
smaller ones (Higham 2008; Moler and Van Loan 2003). MATLAB R2023b and
later provide `expmv` for a matrix-exponential action. The REXI method of
Schreiber, Peixoto, Haut and Wingate (2018) and early Laplace-transform
PinT methods use related exponential-approximation ideas.

The wave-equation experiment of Gander and Güttel (2013) on problem (2.7) reported time-parallel efficiency up to 80%, and the paper draws an explicit conclusion from it: **ParaExp is an excellent time parallelization method for linear hyperbolic problems**. The particular number still depends on the exponential algorithm, matrix structure, partition, and hardware.

### Nonlinear splitting

Assume the nonlinear dynamics can be written as

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B(\boldsymbol u(t))+\boldsymbol g(t). \tag{3.17}
$$

The nonlinear extension is due to Gander, Güttel and Petcu (2018a).
Writing $\boldsymbol u=\boldsymbol w+\boldsymbol v$, the exact split
required to preserve the original equation is

$$
\boldsymbol w'(t)=A\boldsymbol w(t),
\qquad
\boldsymbol v'(t)=A\boldsymbol v(t)
+B\!\left(\boldsymbol v(t)+\boldsymbol w(t)\right)+\boldsymbol g(t),
\qquad \boldsymbol v(0)=\boldsymbol 0,
$$

because $(\boldsymbol w+\boldsymbol v)'=
A(\boldsymbol w+\boldsymbol v)+B(\boldsymbol w+\boldsymbol v)+\boldsymbol g$.
Unlike (3.14), the time intervals are now coupled: the initial value of
$\boldsymbol w$ at $T_{n-1}$ depends on $\boldsymbol v(T_{n-1})$.

> [!warning] Source check: the nonlinear split
> The journal and arXiv versions first omit $A\boldsymbol v$ from the
> $\boldsymbol v$ equation. In the following unnumbered interval
> iteration they then print $A\boldsymbol u_n^k$ and
> $\boldsymbol u_n^k(T_{n-1})=0$ where consistency requires
> $A\boldsymbol v_n^k$ and $\boldsymbol v_n^k(T_{n-1})=0$. Those
> expressions do not satisfy
> $\boldsymbol u=\boldsymbol v+\sum_j\boldsymbol w_j$ as printed; this
> page follows the variable definitions.

Starting the homogeneous problems from the previous interface iterate
gives the first parallel formulation:

$$
\begin{aligned}
(\boldsymbol w_n^k)'&=A\boldsymbol w_n^k,
&&t\in[T_{n-1},T],\\
\boldsymbol w_1^k(T_0)&=\boldsymbol u_0,
&
\boldsymbol w_n^k(T_{n-1})
&=\boldsymbol v_{n-1}^{k-1}(T_{n-1}),
&&n=2,\ldots,N_t,\\
(\boldsymbol v_n^k)'&=A\boldsymbol v_n^k+
B\!\left(\boldsymbol v_n^k+\sum_{j=1}^n\boldsymbol w_j^k\right)
+\boldsymbol g,
&&t\in[T_{n-1},T_n],\\
\boldsymbol v_n^k(T_{n-1})&=\boldsymbol 0,\\
\boldsymbol u_n^k&=\boldsymbol v_n^k+\sum_{j=1}^n\boldsymbol w_j^k,
&&n=1,\ldots,N_t.
\end{aligned}
$$

Explicitly evaluating the full blue sum
$\sum_j\boldsymbol w_j^k(t)$ inside every nonlinear local solve repeats
expensive work when $A$ is large. Substituting
$\boldsymbol v_n^k=\boldsymbol u_n^k-\sum_{j=1}^n\boldsymbol w_j^k$
eliminates this explicit dependence and yields (3.18)–(3.19).

First construct the homogeneous propagations for all
$n=1,\ldots,N_t$:

$$

Every initial value contains only iteration-$k-1$ data, so the
different values of $n$ are independent during iteration $k$ and can
be solved concurrently.

> [!warning] Source check: `sequentially`
> The journal text says “sequentially” before (3.18), but the right-hand
> side of (3.18) contains only previous-iteration data. That word also
> conflicts with the “in parallel” descriptions immediately before and
> after it. This page follows the algorithm's data dependencies.
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

**Finite-step result.** Theorem 3.4, taken from Gander et al. (2018a), states that after iteration $k$ the iterate $\boldsymbol u^k(t)$ agrees with the exact solution on $[0,T_k]$, so the iteration converges in a finite number of steps.

> [!note] Site supplement: the induction behind it
> The paper states only the result. Inducting on the iteration counter $k$ makes the mechanism visible: at $k=1$ the first interval starts from the true initial value, so $[0,T_1]$ is exact; if $[0,T_{k-1}]$ is exact after iteration $k-1$, then (3.18) constructs the exact initial value at $T_{k-1}$ and (3.19) advances the exact interval by one more subinterval to $T_k$.

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

Standard Parareal normally lets its coarse propagator approximate the
complete nonlinear problem (3.20c) as well. Here $\mathcal G$ retains
only $A$, so this is a simplified Parareal. The equivalence holds only
at the coarse nodes:

$$
\boldsymbol u^k(T_n)=\boldsymbol u_n^k(T_n)=\boldsymbol U_n^k,
\qquad n=0,1,\ldots,N_t,\qquad \boldsymbol U_0^k=\boldsymbol u_0.
$$

The local trajectories do not coincide pointwise with Parareal
trajectories. Chapter 4 discusses standard Parareal in detail.

### Figure 3.8: splitting failure on Burgers' equation

Periodic Burgers' equation tests the nonlinear split:

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B\boldsymbol u^2(t),
\qquad t\in(0,2), \tag{3.21}
$$

Here $\boldsymbol u^2$ is componentwise and $A_{xx},A_x$ are the
dimensionless stencils in (3.12). With $\Delta x=1/100$, a consistent
semidiscretization is

$$
A=\frac{\nu}{\Delta x^2}A_{xx},
\qquad
B=-\frac{1}{4\Delta x}A_x.
$$

The second coefficient follows from
$-\frac12\partial_x(u^2)\approx
-\frac1{4\Delta x}A_x\boldsymbol u^2$. Both methods use backward Euler
as the fine propagator with step $0.01/20$. Standard Parareal also uses
backward Euler as its coarse propagator with step $0.01$; ParaExp uses
MATLAB `expmv` for the linear coarse propagation.

> [!warning] Source check: spatial scaling in Figure 3.8
> The journal and arXiv versions write $A=A_{xx}$ and $B=-A_x/2$,
> omitting the factors involving $\nu$ and $\Delta x$. Read literally,
> $A$ would not depend on $\nu$, making the three-viscosity experiment
> impossible. The display above follows model equation (2.6) and the
> stencil definitions in (3.12).

![Source Figure 3.8: errors of nonlinear ParaExp and standard Parareal at three viscosities](assets/papers/time-parallelization/source-figures/figure-3-8.svg)

The panels from left to right use $\nu=1,0.1,0.02$. Their horizontal
lines mark the truncation level $\max\{\Delta t,\Delta x^2\}$, beyond
which further iteration is unnecessary. At $\nu=1$, the linear part
captures the dominant dynamics and ParaExp is substantially faster than
standard Parareal. At $\nu=0.1$ the ordering reverses. At $\nu=0.02$
the ParaExp error grows at every displayed step, whereas standard
Parareal still crosses the truncation line; it too eventually fails as
viscosity falls further. The important change is not a uniform slowdown
under weaker diffusion but a transfer of dominant dynamics between the
two parts of the split: the simplified coarse propagation loses
representative power before standard Parareal does.

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
