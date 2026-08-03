---
title: "3.4: ParaExp"
description: A complete derivation of exact linear superposition in ParaExp, a matrix-exponential-action toolbox, nonlinear iteration, and the equivalence with Parareal, with supplementary references
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-3-3-paraexp
tags:
  - parallel-in-time
  - ParaExp
  - matrix-exponential
---

> [!note] Reading scope
> This page corresponds to Section 3.4 of the paper (pp. 412–415) and covers equations (3.13)–(3.21), Theorem 3.4, and Figures 3.7–3.8. The induction proof of the linear reconstruction, the algorithmic toolbox for the matrix-exponential action, the finite-step property of the nonlinear version, and its equivalence with Parareal are all developed point by point. The original method is due to Gander and Güttel (2013), and the nonlinear generalization to Gander et al. (2018a).

## 3.4 ParaExp

### Two families of parallel subproblems for the linear problem

ParaExp targets

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol g(t),
\qquad \boldsymbol u(0)=\boldsymbol u_0.
$$

ParaExp is a **direct** time-parallel method: rather than approaching the solution iteratively step by step, it reconstructs the solution algebraically in one shot (by contrast, Parareal, waveform relaxation, and the like are iterative). This idea of “splitting a long time interval, solving in parallel, and stitching the pieces back together” goes back to the earliest visions of time parallelism (Nievergelt 1964; see the supplementary references); the contribution of ParaExp is to find a split whose stitching is **exact** and whose required propagation is cheap enough.

Partition the time interval into $[T_{n-1},T_n]$, $n=1,\ldots,N_t$. The first, “red,” family keeps the source and uses a zero initial value on every subinterval:

$$
\boldsymbol v_n'(t)=A\boldsymbol v_n(t)+\boldsymbol g(t),
\quad t\in(T_{n-1},T_n],
\qquad
\boldsymbol v_n(T_{n-1})=0. \tag{3.13}
$$

These problems are mutually independent and can all run in parallel. The second, “blue,” family removes the source and propagates the red endpoint from the previous interval forward to the global final time:

$$
\boldsymbol w_n'(t)=A\boldsymbol w_n(t),
\quad t\in(T_{n-1},T],
\qquad
\boldsymbol w_n(T_{n-1})=\boldsymbol v_{n-1}(T_{n-1}), \tag{3.14}
$$

where $\boldsymbol v_0(T_0)=\boldsymbol u_0$. All blue problems are likewise mutually independent.

![Source Figure 3.7: local red forced problems and homogeneous blue tail propagation in ParaExp](assets/papers/time-parallelization/source-figures/figure-3-7.svg)

The red curves in Figure 3.7 cover only their local subinterval, while every blue dashed curve extends from each interface all the way to $T$. On interval $n$, the exact solution is the superposition of the local red response and all blue responses that have already started:

$$
\boldsymbol u(t)=\boldsymbol v_n(t)
+\sum_{j=1}^{n}\boldsymbol w_j(t),
\qquad t\in[T_{n-1},T_n]. \tag{3.15}
$$

![ParaExp separates local forced responses from global homogeneous propagation](assets/diagrams/pint/en/paraexp-decomposition.svg)

> [!tip] Insight
> These two families of subproblems are exactly the parallelization of the “variation-of-constants” structure of the linear ODE solution: the red $\boldsymbol v_n$ carries the **forced particular part** (right-hand side $\boldsymbol g$, zero initial value), and the blue $\boldsymbol w_n$ carries the **homogeneous initial-value propagation part** (no right-hand side, carrying the interface value). The red zero initial value makes the intervals independent of one another, so they parallelize immediately; the price is that each interval discards its “history from the past,” and that history is exactly what the blue tails restore. Understanding this division of labor makes it clear that (3.15) is not an approximation but an identity (proved in the next section), and it also explains why the real cost concentrates in the long-range blue propagation (the following section explains why it is still cheap).

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

Uniqueness of the solution gives $\boldsymbol u=\boldsymbol v_1+\boldsymbol w_1$, so (3.15) holds for $n=1$.

Assume it holds on interval $n$. Then the interface value satisfies

$$
\boldsymbol u(T_n)=\boldsymbol v_n(T_n)
+\sum_{j=1}^{n}\boldsymbol w_j(T_n).
$$

By (3.14), $\boldsymbol w_{n+1}(T_n)=\boldsymbol v_n(T_n)$, hence

$$
\boldsymbol u(T_n)=\sum_{j=1}^{n+1}\boldsymbol w_j(T_n).
$$

Let $\boldsymbol w=\sum_{j=1}^{n+1}\boldsymbol w_j$. On $(T_n,T_{n+1}]$, $\boldsymbol w'=A\boldsymbol w$ with initial value $\boldsymbol u(T_n)$. Adding the zero-initial forced solution $\boldsymbol v_{n+1}$ yields the unique solution of the original initial-value problem on the next interval. The induction is complete.

> [!tip] Insight
> Each induction step uses only two facts: **linearity** (adding two subproblems still satisfies the same equation) and **uniqueness of the solution** (whatever satisfies the equation with the correct initial value is the true solution). This shows that (3.15) is an **algebraic identity**, not a limit in the sense of convergence—as long as each subproblem is solved exactly, the superposition is exact in one shot, and there is no “number of iterations” dimension. The key interface hand-off used repeatedly in the induction is supplied by $\boldsymbol w_{n+1}(T_n)=\boldsymbol v_n(T_n)$: the history accumulated by the red solution on interval $n$ becomes precisely the “starting point” of the $(n+1)$-th blue tail, so the entire history is relayed interval by interval to the end of the time domain. For this reason, the error of linear ParaExp can only come from two sources: the approximation error of the matrix-exponential action, and the local discretization error of the red forced subproblems—and **not from any splitting assumption**.

### Why the long blue tails are still inexpensive

The homogeneous problem has the closed-form expression

$$
\boldsymbol w_n(t)=
\exp\!\left((t-T_{n-1})A\right)
\boldsymbol v_{n-1}(T_{n-1}),
\qquad t\in[T_{n-1},T]. \tag{3.16}
$$

At first glance the earlier blue problems must integrate all the way to the later time $T$, which seems as expensive as the original problem. The key is that (3.16) is homogeneous: all that is needed is the **matrix-exponential action** $e^{\tau A}\boldsymbol b$ ($\tau=t-T_{n-1}$, $\boldsymbol b=\boldsymbol v_{n-1}(T_{n-1})$), which can **jump directly** to the target time, with a cost determined by the spectral properties of $A$ and the required accuracy, and **not proportional to the number of intermediate time steps**. By contrast, the cost of ordinary time-stepping accumulates linearly with the interval length, which is exactly why the long tail looks expensive. This “jumping” ability decouples the cost of long-range propagation from the interval length, and it is the fundamental reason for ParaExp’s high parallel efficiency.

The paper lists mature tools for computing $e^{\tau A}\boldsymbol b$ (for surveys see Higham 2008; Moler and Van Loan 2003), to be chosen according to the size and structure of the matrix:

- **Rational Krylov methods**: aimed at large sparse $A$. They approximate $e^{\tau A}\boldsymbol b$ in a rational Krylov subspace, with each step requiring the solution of a shifted linear system $(A-\sigma I)\boldsymbol x=\boldsymbol y$; convergence is especially fast for stiff/parabolic operators whose spectrum lies near the negative real axis.
- **Chebyshev expansion**: when $A$ is symmetric (or its spectrum falls in a known real interval, as for wave/diffusion operators discretized by centered differences), the scalar function $e^{\tau\lambda}$ is approximated uniformly by Chebyshev polynomials, requiring only repeated matrix–vector products and no linear solves.
- **Scaling-and-squaring with Padé approximation**: a dense algorithm that forms the **complete** matrix exponential $e^{\tau A}$, suitable for smaller dense matrices.
- **`expmv`-type action algorithms**: they compute the action $e^{\tau A}\boldsymbol b$ directly without explicitly forming the whole exponential matrix (the underlying idea is in Al-Mohy and Higham 2011; see the supplementary references), and MATLAB R2023b and later include `expmv`.
- **REXI / early Laplace-transform-based PinT**: they write $e^{\tau A}\boldsymbol b$ as a weighted sum of several shifted resolvents $(A-\sigma_\ell I)^{-1}\boldsymbol b$, whose terms are mutually independent and admit an additional layer of parallelism; this is especially well suited to oscillatory (hyperbolic) problems (for REXI see Schreiber et al. 2018).

> [!tip] Insight
> Choosing among these tools is essentially a trade-off between “whether to form the whole $e^{\tau A}$” and “the spectral geometry of $A$”: small dense matrices can afford scaling-and-squaring to compute everything at once; large sparse matrices should avoid generating the dense exponential and instead use action-type methods that need only matrix–vector products or shifted solves (Krylov / Chebyshev / expmv). The location of the spectrum refines this further: when the spectrum leans toward the negative real axis (strong diffusion), the polynomial approximations of Krylov/Chebyshev converge quickly; when the spectrum hugs the imaginary axis (strong oscillation, waves), resolvent sums such as REXI are more robust. Notably, REXI’s “resolvent sum” is itself another layer of parallel structure, so ParaExp can nest a further level of parallelism internally—this is precisely the technical foundation that lets it approach high parallel efficiency on hyperbolic problems.

The wave-equation experiment of Gander and Güttel (2013) once reported time-parallel efficiency of about 80%. This number depends on the exponential algorithm, matrix structure, partition, and hardware; it expresses the result of that implementation and is not a universal upper bound for the method. With a different exponential-action implementation or a different spectral structure, the efficiency will vary accordingly.

### Nonlinear splitting

Suppose the nonlinear term can be written as

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B(\boldsymbol u(t))+\boldsymbol g(t). \tag{3.17}
$$

The direct superposition of the linear case fails here, because $B$ recouples the intervals. The initial construction following the linear idea writes $\boldsymbol u=\boldsymbol w+\boldsymbol v$, lets the homogeneous part $\boldsymbol w'=A\boldsymbol w$ carry the initial value $\boldsymbol w(0)=\boldsymbol u_0$, and lets the nonlinear part satisfy $\boldsymbol v'=B(\boldsymbol v+\boldsymbol w)+\boldsymbol g$ with zero initial value $\boldsymbol v(0)=0$; their sum still solves (3.17). The trouble is that $\boldsymbol v+\boldsymbol w$ appears inside the $B$ term, so the subintervals become entangled through the interface value $\boldsymbol v(T_{n-1})$ and cannot be parallelized in one shot as in the linear case—they must be turned into an iteration.

If one explicitly computes the full blue tail $\sum_j\boldsymbol w_j^k(t)$ in every nonlinear local problem, then each interval must integrate the linear problem over the whole of $[T_{n-1},T]$, and a large $A$ causes redundancy and waste. The paper instead uses $\boldsymbol v_n^k=\boldsymbol u_n^k-\sum_{j=1}^n\boldsymbol w_j^k$ to eliminate this redundancy, obtaining a two-stage iteration that requires only local integration.

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

Then solve the complete nonlinear problems in parallel over all time subintervals:

$$
\begin{aligned}
(\boldsymbol u_n^k)'(t)
&=A\boldsymbol u_n^k(t)+B(\boldsymbol u_n^k(t))+\boldsymbol g(t),
&&t\in[T_{n-1},T_n],\\
\boldsymbol u_n^k(T_{n-1})
&=\sum_{j=1}^{n}\boldsymbol w_j^k(T_{n-1}).
\end{aligned} \tag{3.19}
$$

The global approximation at iteration $k$ takes $\boldsymbol u^k(t)=\boldsymbol u_n^k(t)$ on interval $n$.

> [!tip] Insight
> The cleverness of the substitution $\boldsymbol v_n^k=\boldsymbol u_n^k-\sum_{j=1}^n\boldsymbol w_j^k$ lies in replacing an “expensive and redundant” quantity: in the original form $B$ depends on $\sum_j\boldsymbol w_j^k$, forcing the linear tail to be re-integrated to the end of the interval in **every** local nonlinear problem; after the substitution, the nonlinear problem (3.19) solves the complete equation $\boldsymbol u'=A\boldsymbol u+B(\boldsymbol u)+\boldsymbol g$ only locally on $[T_{n-1},T_n]$, while the linear tail (3.18) does a one-shot propagation independently. Note that the blue initial value in (3.18) takes the previous iteration’s $\boldsymbol u_{n-1}^{k-1}-\sum_{j}\boldsymbol w_j^{k-1}$, which is exactly the eliminated $\boldsymbol v_{n-1}^{k-1}(T_{n-1})$; the nonlinear initial value in (3.19) takes the current iteration’s $\sum_j\boldsymbol w_j^k(T_{n-1})$. The linear propagation and the nonlinear solve are thus completely decoupled into two separately parallelizable steps, which is also the structural basis for writing it in Parareal form in the next section.

### Theorem 3.4: finite-step convergence and Parareal equivalence

**Finite-step result.** After iteration $k$, $\boldsymbol u^k(t)$ coincides with the exact solution on $[0,T_k]$; that is, iterative ParaExp converges in a finite number of steps. The reason follows by induction over time intervals: the first interval always starts from the true initial value $\boldsymbol u_0$, so it is exact on $[0,T_1]$ from iteration 1 onward; if the first $k-1$ iterations are already exact on $[0,T_{k-1}]$, then (3.18) constructs the exact initial value at the next interface, and (3.19) advances correctness by one more interval, so iteration $k$ reaches $T_k$. Information advances one coarse interval along the time direction each iteration, so at most $N_t$ iterations yield the exact solution—this is fundamentally different from the “asymptotic convergence” of iterative methods; it is **finite-step termination**.

> [!tip] Insight
> The finite-step property is really the “afterimage” of the linear exact superposition (3.15) in the nonlinear case: the nonlinear $B$ destroys the one-shot superposition, but each iteration can still “nail down” the solution on one coarse interval to its exact value, so convergence advances interval by interval rather than decaying globally. This also explains the shape of the convergence curve—the error typically drops in a staircase pattern until it bottoms out at iteration $N_t$. In practice one need not run the full $N_t$ iterations: one can stop as soon as the error drops to the level of the discretization truncation error (see the horizontal lines in Figure 3.8 below).

At the coarse points $T_n$, the iteration is equivalent to

$$
\boldsymbol U_n^k
=\mathcal G(T_{n-1},T_n,\boldsymbol U_{n-1}^k)
+\mathcal F(T_{n-1},T_n,\boldsymbol U_{n-1}^{k-1})
-\mathcal G(T_{n-1},T_n,\boldsymbol U_{n-1}^{k-1}), \tag{3.20a}
$$

where the coarse propagator solves only the linear homogeneous problem

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

This is the first appearance of Parareal in the paper (detailed in Section 4). The coarse propagator $\mathcal G$ of standard Parareal usually **also approximately solves the complete nonlinear problem (3.20c)** (just with a coarser step size), so it is cheap yet still captures the main dynamics. Here $\mathcal G$ retains only the linear part $A$ and discards $B$ (i.e., it solves (3.20b) rather than (3.20c)), which is a **simplified version**. In other words, the “coarse model” of ParaExp is a linearization of the original problem, not a coarse discretization.

> [!tip] Insight
> Recognizing nonlinear ParaExp as a Parareal in which “$\mathcal G$ solves only the linear problem while $\mathcal F$ solves the complete nonlinear one” makes its key to success or failure immediately clear: the convergence of Parareal depends on how well the coarse propagator **approximates** the fine propagator. Standard Parareal discretizes the **same** nonlinear operator with a coarse step, so its error comes mainly from temporal resolution; ParaExp’s $\mathcal G$ simply throws away $B$, so its error comes from **model mismatch**—as long as $B$ carries the main propagation mechanism (such as strong advection), the linear $\mathcal G$ is far from $\mathcal F$, the correction term $\mathcal F-\mathcal G$ is no longer small, and the iteration slows down or even diverges. Hence “whether the $A+B$ split places the dominant dynamics into $A$” directly determines the performance of nonlinear ParaExp, a point that shows up very intuitively in the Burgers experiment below. The paper also hints in advance that standard Parareal is already poor at hyperbolic problems, so the simplified version can hardly be expected to work in that regime.

### Figure 3.8: splitting failure on Burgers' equation

The paper uses

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B\boldsymbol u^2(t),
\qquad t\in(0,2), \tag{3.21}
$$

which comes from centered differences of the periodic Burgers equation, with $\Delta x=1/100$, $A=\nu A_{xx}/\Delta x^2$, $B=-A_x/(2\Delta x)$, and $A_{xx},A_x$ as in (3.12). Here the linear part $A$ corresponds to the viscous diffusion term, scaled by the viscosity $\nu$, and the nonlinear part $B\boldsymbol u^2$ corresponds to the advection term. The fine propagators of both ParaExp and standard Parareal use backward Euler with fine step $0.01/20$. Standard Parareal still uses backward Euler as its coarse propagator, with coarse step $0.01$; ParaExp’s linear coarse propagation calls MATLAB `expmv`.

![Source Figure 3.8: errors of nonlinear ParaExp and standard Parareal at three viscosities](assets/papers/time-parallelization/source-figures/figure-3-8.svg)

The three panels take, from left to right, $\nu=1,0.1,0.02$; the horizontal line marks the discretization truncation error $\max\{\Delta t,\Delta x^2\}$, and in practice the iteration can stop once it reaches this line. At $\nu=1$ (strong diffusion), $A$ captures the main dynamics and ParaExp is clearly faster than standard Parareal; at $\nu=0.1$, standard Parareal is instead faster, while ParaExp still decays slowly; at $\nu=0.02$ (advection-dominated), **ParaExp fails first**—its error keeps growing, i.e., diverges, whereas standard Parareal can still cross the truncation-error line. On further reducing the viscosity, standard Parareal will eventually fail too (see Section 4 for details). These three panels display a shift of the dominant term of the split as the viscosity changes, and cannot be summarized merely as “convergence slows as diffusion weakens.”

> [!tip] Insight
> The three values of $\nu$ sweep exactly across the phase transition “$A$ dominant → evenly matched → $B$ dominant”: at $\nu=1$ diffusion (in $A$) dominates, ParaExp’s linear coarse model captures almost all the dynamics and is therefore fastest; at $\nu=0.02$ advection (in $B$) dominates and is entirely ignored by ParaExp’s $\mathcal G$, so it **diverges earlier than standard Parareal**—because standard Parareal’s coarse propagator, after all, discretizes the complete equation including advection, just at a coarse resolution. This confirms the judgment of the previous section: whoever puts the dominant term into the coarse model wins on the corresponding regime. The engineering implication is that, to use ParaExp on advection–diffusion-type problems, one should try to include advection in $A$ too (for example by re-splitting $A+B$ around a linearized operator), otherwise no matrix exponential, however accurate, can make up for the model mismatch.

> [!important] Applicability boundary
> The linear ParaExp equation (3.15) is an exact algebraic decomposition, affected only by the exponential-action approximation and the local forced-solve error. Nonlinear ParaExp has become an iterative method, and its convergence depends on the $A+B$ split. If the linear part does not cover the dominant propagation mechanism, no amount of exponential accuracy can compensate for the model mismatch.

## Equation and figure coverage

| Source item                  | Paper section | Coverage                                                                                                                 |
| ---------------------------- | ------------- | ------------------------------------------------------------------------------------------------------------------------ |
| (3.13)–(3.16), Figure 3.7    | 3.4           | both subproblem families, exact superposition, induction proof, exponential-action toolbox, and complete source figure   |
| (3.17)–(3.19)                | 3.4           | nonlinear split, initial decoupling, source of redundancy, and the two-stage iteration after substitution                |
| (3.20a)–(3.20c), Theorem 3.4 | 3.4           | finite-step result, induction argument, and the difference from the standard Parareal coarse propagator                  |
| (3.21), Figure 3.8           | 3.4           | Burgers discretization, all temporal parameters, truncation-error line, and the three-regime phase-transition conclusion |

## Further reading (supplementary references)

- The earliest vision of direct time parallelism: J. Nievergelt, _Parallel methods for integrating ordinary differential equations_, _Comm. ACM_ 7 (1964).
- The original linear ParaExp method and wave-equation efficiency: M. J. Gander and S. Güttel, _ParaExp: A parallel integrator for linear initial-value problems_, _SIAM J. Sci. Comput._ 35 (2013).
- Surveys and algorithms for the matrix-exponential action: N. J. Higham, _Functions of Matrices_ (SIAM 2008); C. Moler and C. Van Loan, _Nineteen dubious ways to compute the exponential of a matrix, twenty-five years later_, _SIAM Review_ 45 (2003); A. Al-Mohy and N. J. Higham, _Computing the action of the matrix exponential_, _SIAM J. Sci. Comput._ 33 (2011).
- Exponential-type PinT and resolvent sums: M. Schreiber, P. Peixoto, T. Haut and B. Wingate (2018) (REXI).

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3.4, pp. 412–415.
