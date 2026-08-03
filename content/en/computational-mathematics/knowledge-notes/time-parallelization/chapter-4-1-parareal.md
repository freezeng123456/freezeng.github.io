---
title: "4.1–4.2: Historical Context and Parareal"
description: A step-by-step close reading from multiple shooting and virtual control to coarse/fine propagators, covering Parareal's linear/nonlinear convergence bounds, the 0.3 parabolic factor, and the weak-diffusion/hyperbolic failure mechanism
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-1-parareal
tags:
  - parallel-in-time
  - Parareal
  - parabolic equations
---

> [!note] Reading scope
> This page corresponds to Sections 4, 4.1, and 4.2 of the paper (pp. 443–452), covering equations (4.1)–(4.9), Theorems 4.1–4.4, Remark 4.1, and Figures 4.1–4.5. The reasoning after each formula follows the proof order in the paper, with added "why" explanations and references at key points.

## 4.1 Historical development

### Why the paper separates methods designed for parabolic problems

Chapter 2 showed that, except at very low frequencies, parabolic equations exhibit strong temporal locality, whereas every frequency of a hyperbolic equation can propagate over long times. The Chapter 3 methods treat long-range coupling in time simultaneously, so they can often handle parabolic problems as well; yet they retain clear limitations in nonlinear settings: the optimized Robin parameter of OSWR is hard to determine, and the outer Newton iterations of ParaExp and the two kinds of ParaDiag slow down—or even fail—as the time window grows.

The Chapter 4 methods actively exploit the temporal locality induced by dissipation. Parareal, PFASST, MGRiT, and STMG are effective on both linear and nonlinear problems when diffusion is sufficient; transplanting them directly to weakly diffusive or hyperbolic problems makes convergence slow down continuously and possibly diverge. This page traces that central theme—"the stronger the diffusion the better, the stronger the propagation the worse"—from the update formula all the way through every step of the convergence factor.

The idea of Parareal can be traced back to multiple shooting, waveform relaxation, and the noniterative precursor of Nievergelt (1964). Lions et al. (2001) independently proposed the modern algorithm. Later came PITA, PFASST, MGRiT, and Parareal–ParaDiag combinations. Another line is space–time multigrid: early methods struggled to coarsen time effectively; Gander and Neumüller (2016) re-established scalable STMG through a temporal block-Jacobi smoother.

> [!tip] Insight
> The paper points out specifically that Lions–Maday–Turinici (2001) did not propose Parareal within a "multiple shooting" framework, but arrived at the same scheme in the context of **virtual control**. That two distinct origins converge tells us one thing: the $\mathcal F+\mathcal G-\mathcal G$ structure of Parareal is not an accident of some particular derivation, but the inevitable form of the idea "predict with a cheap model, correct with an expensive model." The most effortless way to understand it is to read it as a finite-difference approximation of the Jacobian in the multiple-shooting Newton iteration (Gander and Vandewalle 2007).

## 4.2 Parareal

### The update formula and the two-level temporal grid

Partition $[0,T]$ by $0=T_0<T_1<\cdots<T_{N_t}=T$. Here $\mathcal F$ is the expensive, accurate fine propagator, and $\mathcal G$ is the cheap coarse propagator. Starting from interface initial guesses $\boldsymbol u_n^0$, Parareal iterates

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k). \tag{4.1}
$$

Within the same iteration, all $\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)$ can be computed in parallel; the coarse propagation carrying the new iteration index $k+1$ must advance sequentially in $n$. The last two terms form a prediction–correction on the coarse grid, and can also be interpreted as the finite-difference Jacobian of a multiple-shooting Newton method.

It is worth making the role of each term fully explicit. $\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)$ uses a small step $\Delta t$ and integrates the interval $[T_n,T_{n+1}]$ to the end with $\boldsymbol u_n^k$ as initial value—costly but accurate; $\mathcal G$ completes the same interval with a large step $\Delta T$ (or a simpler model)—cheap but crude. Rearranging (4.1) into

$$
\boldsymbol u_{n+1}^{k+1}
=\underbrace{\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})}_{\text{prediction (serial, cheap)}}
+\underbrace{\big[\mathcal F(T_n,T_{n+1},\boldsymbol u_n^{k})-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k})\big]}_{\text{correction (parallel, using last iteration)}},
$$

makes it clear: the coarse propagation "sweeps" a serial predicted trajectory along time, while the expensive fine–coarse difference depends only on the already-known $\boldsymbol u_n^k$ from the previous iteration, and can therefore be computed simultaneously for all $n$. At convergence $\boldsymbol u_n^{k+1}\to\boldsymbol u_n^k$, the two $\mathcal G$ terms in the bracket cancel, and what remains at the interfaces is precisely the sequential solution of $\mathcal F$—this explains why the fixed point is exactly the sequential fine solution, independent of the accuracy of the coarse propagator $\mathcal G$.

> [!tip] Insight
> The reason the "difference" $\mathcal F-\mathcal G$ can serve as a Jacobian: multiple shooting writes interface continuity as the nonlinear system $\boldsymbol u_{n+1}=\mathcal F(T_n,T_{n+1},\boldsymbol u_n)$, and the Newton step needs $\partial\mathcal F/\partial\boldsymbol u_n$. Computing this Jacobian directly is as expensive as computing $\mathcal F$ itself, and it is serial. Parareal approximates the product of the Jacobian with the increment by "the first-order change of the cheap propagator across adjacent iterations," $\mathcal G(\cdot,\boldsymbol u_n^{k+1})-\mathcal G(\cdot,\boldsymbol u_n^{k})$ (Gander and Vandewalle 2007). This both avoids the expensive tangent propagation and carries Newton's **local superlinear** property into Parareal—the factorial/product factors appearing later in (4.5a) and (4.6) are precisely the quantification of this superlinearity.

![Source Figure 4.1: each coarse time step contains J fine time steps](assets/papers/time-parallelization/source-figures/figure-4-1.svg)

The paper mainly discusses uniform grids, setting $\Delta T/\Delta t=J\ge2$; nonuniform grids can also be used (Gander 2017; Maday and Mula 2020; Wu and Zhou 2024). The target solution is the discrete solution obtained by running $\mathcal F$ sequentially, and Parareal does not additionally change the discrete target beyond a fixed $\mathcal F$. In other words, $J\ge2$ guarantees that one coarse step encloses at least two fine steps, so that the coarse level can be "cheaper" than the fine level; it also means that the fine propagation mode appearing later must be written as the **$J$-th power of the single-step fine stability function** $R_f^J(z/J)$, because within one coarse interval the fine propagator must take $J$ consecutive steps. This is the bookkeeping basis of the entire linear analysis, and it also shows that Parareal is **nonintrusive**: it only repeatedly calls the existing $\mathcal F$ and $\mathcal G$, without altering their internal discrete schemes.

### Theorem 4.1: reducing the error to a mode-wise Toeplitz iteration

Consider $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$, and let $A=V_ADV_A^{-1}$, with coarse and fine single-step stability functions $R_g,R_f$. Within one coarse interval the fine propagator takes $J$ steps, so the fine propagation mode is $R_f^J(z/J)$, where $z=\Delta T\lambda(A)$.

If $|R_g(z)|\le1$, then

$$
\max_{1\le n\le N_t}
\|V_A(\boldsymbol u_n^k-\boldsymbol u_n)\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}\|M^k(z)\|_\infty
\max_{1\le n\le N_t}
\|V_A(\boldsymbol u_n^0-\boldsymbol u_n)\|_\infty. \tag{4.2}
$$

Here $\boldsymbol u_n$ is the sequential solution of the fine propagator,

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

Here $|R_g(z)|\le1$ is the only structural assumption, namely **coarse-propagator stability**. It guarantees that $M_g^{-1}(z)$ (lower triangular, with entries $R_g^j$) does not amplify the error through $|R_g|>1$; the spectrum of parabolic problems lies near the negative real axis, and backward Euler exactly satisfies $|R_g|\le1$, so this assumption holds almost automatically.

#### Proof chain

At the matrix level, (4.1) gives

$$
\boldsymbol u_{n+1}^{k+1}
=R_f^J(\Delta TA/J)\boldsymbol u_n^k
+R_g(\Delta TA)\boldsymbol u_n^{k+1}
-R_g(\Delta TA)\boldsymbol u_n^k.
$$

The sequential fine solution satisfies the same identity. Let $\boldsymbol e_n^k=\boldsymbol u_n-\boldsymbol u_n^k$; subtracting gives

$$
\boldsymbol e_{n+1}^{k+1}
=R_g(\Delta TA)\boldsymbol e_n^{k+1}
+[R_f^J(\Delta TA/J)-R_g(\Delta TA)]\boldsymbol e_n^k.
$$

The key point is that **the sequential fine solution also satisfies the same identity**: adding and then subtracting $R_g(\Delta TA)\boldsymbol u_n$ gives $\boldsymbol u_{n+1}=R_f^J(\Delta TA/J)\boldsymbol u_n+R_g(\Delta TA)\boldsymbol u_n-R_g(\Delta TA)\boldsymbol u_n$. When the two equations are subtracted, the common $R_f^J\boldsymbol u_n$ and $R_g\boldsymbol u_n$ cancel, so the error recursion contains **no source term $\boldsymbol g$** and no initial value—because the initial value at $T_0$ is known, $\boldsymbol e_0^k=0$ holds for all $k$. This step is why the whole analysis can be reduced to a homogeneous linear iteration.

After diagonalizing $A$, each mode $\xi_n^k(z)$ independently satisfies

$$
\xi_{n+1}^{k+1}(z)
=R_g(z)\xi_n^{k+1}(z)
+[R_f^J(z/J)-R_g(z)]\xi_n^k(z),
$$

and

$$
\|V_A\boldsymbol e_n^k\|_\infty
=\max_{z\in\sigma(\Delta TA)}|\xi_n^k(z)|. \tag{4.4}
$$

Because $R_g(\Delta TA)=V_AR_g(\Delta TD)V_A^{-1}$ and $R_f^J(\Delta TA/J)=V_AR_f^J(\Delta TD/J)V_A^{-1}$, in the eigenbasis the matrix recursion **decouples eigenvalue by eigenvalue**, leaving on each $z=\Delta T\lambda$ only a scalar doubly-indexed (in $n$ and in $k$) Toeplitz iteration. This is precisely the mathematical portrait of parabolic problems "converging mode by mode": slow modes ($z$ near $0$) and fast modes (very negative $z$) use the same scheme but converge at entirely different rates.

Stacking over all $n$ into $\boldsymbol\xi^k$, we have

$$
M_g(z)\boldsymbol\xi^{k+1}
=[M_g(z)-M_f(z)]\boldsymbol\xi^k,
$$

so $\boldsymbol\xi^k=M^k(z)\boldsymbol\xi^0$, and taking the maximum over modes and temporal nodes yields (4.2). Thus $\|M^k(z)\|_\infty$ is exactly the convergence factor of applying Parareal to the Dahlquist test equation $u'(t)=\lambda u(t)+g(t)$; the entire analysis then reduces to "studying the powers of the Toeplitz matrix $M(z)$ at a single complex number $z$."

> [!note] Remark 4.1: preconditioner interpretation
> $M(z)=I_t-M_g^{-1}(z)M_f(z)$. The equation $M_fU=b$ is the fine propagator's all-at-once system, and $M_g$ is a coarse-propagation preconditioner. Each fine propagation in the residual can be computed in parallel, while the coarse preconditioning solve proceeds sequentially in time. This viewpoint leads directly to the diagonalized coarse correction of Section 4.5.

Writing Remark 4.1 in more detail: one step of Parareal is exactly one **preconditioned Richardson/stationary iteration** for the all-at-once system $M_f(z)U=b$,

$$
M_g(z)\,\Delta U^k=r^k:=b-M_f(z)U^k,\qquad U^{k+1}=U^k+\Delta U^k,
$$

where the residual component

$$
r_n^k=b_n-\big(u_n^k-\mathcal F(T_{n-1},T_n,u_{n-1}^k)\big)
=b_n-\big(u_n^k-R_f^J(z/J)\,u_{n-1}^k\big)
$$

depends only on $u_{n-1}^k,u_n^k$ from the previous iteration, so it can be computed **simultaneously** for all $n$ (the source of parallelism); while the $M_g$ (lower bidiagonal) solve is a single forward-substitution sweep in time (serial, cheap). This reading places Parareal in the general framework of "the expensive operator $M_f$ forms the residual, the cheap operator $M_g$ does the preconditioning," and it also explains why Section 4.5 need only replace the serial $M_g$ with a diagonalizable coarse correction to parallelize further.

The strictly lower triangular structure further shows: after iteration $k$, the first $k$ coarse nodes already agree with the sequential fine solution, so in exact arithmetic the method terminates in at most $N_t$ iterations. The reason is that $M(z)$ is strictly lower triangular, and after the $k$-th power its first $k$ diagonals vanish, so the error is "locked" to exact values node by node starting from the initial time. This is Parareal's **finite-step termination** property, and the "product containing a zero factor" in (4.5a) and (4.6) is its quantitative version.

### Theorem 4.2: short-time superlinear and long-time linear convergence

The lower triangular entries of $M_g^{-1}$ are $R_g^j(z)$. Therefore

$$
M(z)=[R_f^J(z/J)-R_g(z)]\widetilde M(R_g(z)),
$$

where the first subdiagonal of $\widetilde M(\beta)$ is $1$, and the deeper diagonals are successively $\beta,\beta^2,\ldots$. Taking the infinity norm of its powers gives two kinds of bounds.

Splitting $M(z)$ into "scalar magnitude $\times$ structural matrix" is the technical core of this theorem: all information about $z$ is concentrated in the scalar $R_f^J(z/J)-R_g(z)$, while all information about the window length $N_t$ and iteration count $k$ is concentrated in the powers of the purely structural matrix $\widetilde M(\beta)$. Thus

$$
\|M^k(z)\|_\infty=|R_f^J(z/J)-R_g(z)|^k\,\|\widetilde M^{\,k}(R_g(z))\|_\infty .
$$

Gander and Vandewalle (2007, Lemma 4.4) gives the key estimate for the norm of the powers of $\widetilde M$ (supplied here from the references, to make clear the origin of (4.5a) and (4.5b)):

$$
\|\widetilde M^{\,k}(R_g(z))\|_\infty\le
\begin{cases}
\min\!\left\{\left(\dfrac{1-|R_g(z)|^{N_t-1}}{1-|R_g(z)|}\right)^{\!k},\ \dbinom{N_t-1}{k}\right\}, & |R_g(z)|<1,\\[2ex]
\dbinom{N_t-1}{k}, & |R_g(z)|=1.
\end{cases}
$$

Taking the binomial branch of the $\min$ on the right gives the short-time bound; taking the geometric-series branch (letting $N_t\to\infty$) gives the long-time bound. The boundary between the two branches is precisely the criterion for "how long is a short window and how long is a long window."

For short times or small $N_t$,

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

When $k=N_t$ the product contains a zero, explicitly exhibiting finite-step convergence. Here $\binom{N_t-1}{k}=\frac1{k!}\prod_{j=1}^{k}(N_t-j)$, and when $k=N_t$ (or $k\ge N_t$) the product contains the factor $N_t-N_t=0$, directly giving a zero bound, mutually corroborating the lower-triangular argument after Remark 4.1. The $k$-th power in the numerator combined with the factorial decay $1/k!$ is the source of **superlinearity**: as long as the coarse–fine difference $|R_g-R_f^J|$ is bounded, as $k$ increases $|R_g-R_f^J|^k/k!$ tends to zero faster than any geometric series, so on short windows the error "steepens as it iterates."

If $|R_g(z)|<1$, one also has the $N_t$-independent long-time bound

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

The numerator is the coarse–fine propagation difference, and the denominator is the dissipation margin of the coarse propagator. When the coarse propagator is close to unit modulus and yet cannot approximate the fine propagator, this factor approaches or exceeds $1$.

> [!tip] Insight
> (4.5a) and (4.5b) are not two independent results but two branches of the same $\min$, corresponding to two "bookkeeping methods." On short windows $N_t$ is small and the binomial $\binom{N_t-1}{k}$ is "driven to zero" as $k\to N_t$, so **finite-step, superlinear** behavior dominates; on long windows $N_t$ is large and the binomial branch fails, the geometric series $\big(\tfrac{1}{1-|R_g|}\big)^k$ takes over, so each iteration has a fixed contraction factor $\varrho_l$ and convergence is **linear**. The denominator $1-|R_g(z)|$ is the key to the whole chapter: it is the "dissipation margin" of the coarse propagator on that mode. When diffusion is strong, high modes have $|R_g|\to0$, the margin is close to $1$, and $\varrho_l$ is small; but once the spectrum approaches the imaginary axis (propagation-dominated) so that $|R_g|\to1$, the denominator tends to zero and $\varrho_l\to\infty$—this is the unified explanation for the weak-diffusion/hyperbolic failure discussed later.

![Source Figure 4.2: short-time superlinear and long-time linear convergence stages](assets/papers/time-parallelization/source-figures/figure-4-2.svg)

Figure 4.2 uses the periodic heat equation, zero source, $u_0(x)=\sin^2(2\pi x)$, $\Delta x=1/5$, backward Euler on both levels, and $J=10$. For $T=0.02,N_t=6$, $\varrho_s$ accurately describes the superlinear decrease; over a longer interval the error decreases at an approximately fixed slope, and $\varrho_l$ is more appropriate. Refining the spatial grid to $\Delta x=1/8$ also enters the linear stage earlier. Intuitively, the larger $N_t$ is, the later the available "drive-to-zero" step, so the error already goes through a long linear stage before reaching finite-step termination; and spatial refinement raises the maximum $|z|$, letting the spectrum reach the region where the linear factor dominates earlier, so it likewise exhibits linear convergence sooner.

### Theorem 4.3: nonlinear superlinear bound

Let $\mathcal F$ be exact propagation and $\mathcal G$ an order-$p$ method with local truncation error at most $C_3\Delta T^{p+1}$. Assume

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

with each coefficient continuously differentiable in $\boldsymbol v$. Then

$$
\|\boldsymbol u(T_n)-\boldsymbol u_n^k\|
\le
\frac{C_3\Delta T^{p+1}(C_1\Delta T^{p+1})^{k+1}}{(k+1)!}
(1+C_2\Delta T)^{n-k-1}
\prod_{j=0}^{k}(n-j). \tag{4.6}
$$

When $k\ge n$ the product contains a zero, and the first $n$ interfaces are already exact. For small $\Delta T$, each additional iteration brings another factor of $\Delta T^{p+1}$, which is the source of the nonlinear superlinear stage.

This nonlinear result (Gander and Hairer 2008, Theorem 1; see also Gander and Lunet 2024, Theorem 2.6) is proved using **generating functions**, matching the three features of (4.5a) one-to-one on a general nonlinear ODE: the product $\prod_{j=0}^k(n-j)$ together with the factorial $1/(k+1)!$ gives the same type of **finite-step + superlinear** behavior as in the linear case; the newly added $C_1\Delta T^{p+1}$ per iteration comes from the order-$p$ consistency of $\mathcal F-\mathcal G$ (the more accurate the coarse propagator, the higher $p$, and the smaller $\Delta T$, the stronger the gain per iteration); the power $(1+C_2\Delta T)^{n-k-1}$ of the Lipschitz constant $1+C_2\Delta T$ is the **amplification of error propagating along time**—over the finite time $n\Delta T\le T$ it is uniformly controlled by $e^{C_2T}$, so it does not break convergence, which is precisely the nonlinear generalization of the linear requirement $|R_g|\le1$ (the coarse propagator does not amplify error). A more refined backward-error analysis on Hamiltonian systems appears in Gander and Hairer (2014).

> [!tip] Insight
> In (4.6), the fact that the **power of $\Delta T^{p+1}$ grows linearly with $k$** (the $k$-th iteration has order roughly $\Delta T^{(p+1)(k+1)}$) is another viewpoint for understanding "why Parareal is so efficient on diffusive problems": as long as $\Delta T$ lies in the range that makes $C_1\Delta T^{p+1}<1$, a few iterations can push the error below fine-solution accuracy—and this is exactly the condition possessed by problems where dissipation quickly smooths high modes, making the nonlinear local behavior "close to linear and mild." Conversely, if the solution contains long-lived, undamped oscillations, both $C_1$ (depending on the coefficients of $\mathcal F-\mathcal G$) and $C_2$ (the Lipschitz constant) grow, and the superlinear window narrows accordingly.

### Theorem 4.4: the parabolic long-time factor near 0.3

If the coarse propagator $\mathcal G$ uses backward Euler and the fine propagator $\mathcal F$ uses an L-stable Runge–Kutta method, then there exists $J_{\min}=O(1)$ such that

$$
\max_{z\in\mathbb R_-}\varrho_l(J,z)\approx0.3,
\qquad J\ge J_{\min}. \tag{4.7}
$$

This constant does not grow with $T$ or $N_t$. It comes from the negative real spectrum, the high-frequency dissipation of the coarse level, and the bounded coarse–fine difference. Using combinations such as Radau IIA, the worst factor can drop further to about $0.068$.

It is worth making clear "why it stabilizes precisely at $\approx0.3$, independent of $T,N_t$." The parabolic discrete matrix $A$ is negative semidefinite, with spectrum in $z\in\mathbb R_-$. Looking at the two ends of $\varrho_l(J,z)=|R_g(z)-R_f^J(z/J)|/(1-|R_g(z)|)$:

- **Low-frequency end $z\to0^-$:** both coarse and fine stability functions satisfy consistency $R\approx 1+z+O(z^2)$, so $R_g(z)\approx R_f^J(z/J)$, and the numerator $\to0$ at a higher order than the denominator, so $\varrho_l\to0$. The slow modes are essentially corrected in one iteration.
- **High-frequency end $z\to-\infty$:** L-stability means $R_f(z)\to0$, hence $R_f^J(z/J)\to0$; backward Euler has $R_g(z)=1/(1-z)\to0$. So the numerator $|R_g-R_f^J|\to0$, the denominator $1-|R_g|\to1$, and $\varrho_l\to0$. The high modes are "zeroed out in place" by the strong dissipation of the coarse level, with no need for an exact match.
- **Mid-frequency range:** the worst value occurs at some finite $z$ between the two ends, numerically about $0.3$, independent of whether the spectrum extends to $-\infty$, so it holds for arbitrarily long $T,N_t$.

The role of $J\ge J_{\min}=O(1)$ is to ensure that $R_f^J(z/J)$ has already entered its L-stable asymptotic regime (each fine step's $z/J$ is small enough, and the $J$-th power "flattens" the high frequencies enough), so that the supremum in (4.7) stabilizes near $0.3$. When $\mathcal F$ is exact propagation $\mathcal F=\exp(\Delta TA)$, $J\ge2$ already achieves $\approx0.3$. The proof of Theorem 4.4 is scattered across the literature by case: $\mathcal F$ as backward Euler in Mathew, Sarkis and Schaerer (2010); the trapezoidal rule/BDF2/two SDIRK methods in Wu (2015), Wu and Zhou (2015); a general L-stable $\mathcal F$ in Yang, Yuan and Zhou (2023); and the earlier continuous-level roots in Gander and Vandewalle (2007, Table 5.1).

If the fine propagator is only A-stable, such as the trapezoidal rule, high frequencies do not vanish as $|z|\to\infty$. In this case, over a bounded spectral interval $[0,z_{\max}]$ one still has

$$
\max_{z\in[0,z_{\max}]}\varrho_l(J,z)\approx0.3,
\qquad
J\ge J_{\min}=O(\log_2 z_{\max}). \tag{4.8}
$$

The fine level needs more small steps for its high-frequency behavior to match the physical dissipation. The WHY here is: the trapezoidal rule $R_f(z)=\dfrac{1+z/2}{1-z/2}$ has $|R_f|\to1$ (approaching $-1$) as $z\to-\infty$, so high frequencies are neither dissipated nor free of sign oscillation, and $R_f^J(z/J)$ no longer tends to zero at large $|z|$; thus the numerator $|R_g-R_f^J|$ stays $O(1)$ and $\varrho_l$ cannot be pushed down to $0.3$. The remedy is to increase $J$: shrinking each fine step's argument to $z/J$ so that, once $z/J$ falls back into the region where the trapezoidal rule behaves well, $R_f^J$ can again approximate dissipative behavior; the required $J$ grows logarithmically with the spectral width $z_{\max}$, as $O(\log_2 z_{\max})$ (Wu and Zhou 2015, proved for the trapezoidal rule and fourth-order Gauss–RK). This also explains the essential difference between A-stability and L-stability: L-stability naturally "kills high frequencies," while A-stability must "borrow" high-frequency dissipation by refining the fine steps.

The paper also gives two SDIRK methods:

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

For SDIRK22, (4.7) holds at $J_{\min}=2$; SDIRK23 requires $J_{\min}=4$. The notation "SDIRK$sp$" denotes a singly diagonally implicit RK method with $s$ stages and order $p$: SDIRK22 is L-stable, two-stage and second-order, and immediately satisfies $\approx0.3$; SDIRK23 (Wu and Zhou 2015 give $J_{\min}=4$) is higher order but needs a larger coarse-to-fine ratio to pull the worst factor back to $0.3$. This agrees with the intuition of (4.8)—the fine propagator's **stability type** and the **coarse-to-fine ratio $J$** must be judged together, and order alone is insufficient to predict Parareal's convergence rate.

![Source Figure 4.3: differences in Parareal convergence from different fine propagators and coarsening ratios](assets/papers/time-parallelization/source-figures/figure-4-3.svg)

Figure 4.3 uses the periodic heat equation, $\Delta x=1/256$, $\Delta T=0.1$, $T=4$, diffusion coefficient $0.1$. The three panels from left to right correspond to $J=2,10,50$. At $J=2$, the trapezoidal rule stalls near $10^{-4}$, and SDIRK23 is also clearly slower than SDIRK22; at $J=10$, the two SDIRK curves already nearly coincide, while the trapezoidal rule still retains a slower tail; at $J=50$, all three fine propagators hug the $0.3^k$ reference slope. The changes in the figure verify the qualifying conditions of Theorem 4.4: the fine propagator's stability type and the coarse-to-fine ratio must be judged together. The trapezoidal rule "stalling" at small $J$ is precisely a portrait of (4.8)—it is only A-stable, the high frequencies are not dissipated, and Parareal tries to approximate a high-frequency-inaccurate fine solution with backward Euler (which gives strong, "physically correct" dissipation for high frequencies); the contradiction between the two produces an error plateau. Increasing $J$ brings the fine solution's high-frequency behavior back to the dissipative side, and the plateau disappears.

### Degradation as diffusion weakens

The paper then fixes $T=4$, $\Delta T=0.1$, $\Delta x=1/128$, $J=32$, with backward Euler on the coarse level and SDIRK22 on the fine level, and examines the periodic advection–diffusion and Burgers' equations.

![Source Figure 4.4: long-time factors for each advection–diffusion eigenvalue under three viscosities](assets/papers/time-parallelization/source-figures/figure-4-4.svg)

The spectrum of advection–diffusion gradually spreads from near the negative real axis toward the imaginary axis. As $\nu$ decreases, $\max\varrho_l$ approaches $1$, showing that the coarse propagator finds it increasingly difficult to correct the long-lived propagation modes. Reading through the denominator of (4.5b) is the most direct: the eigenvalues of the advection–diffusion operator are $\lambda\approx-\nu\kappa^2+\mathrm i c\kappa$ ($\kappa$ the wavenumber), where viscosity $\nu$ provides the real part (dissipation) and advection $c$ provides the imaginary part (propagation). When $\nu$ is large the real part dominates, $|R_g|$ is far from $1$, the dissipation margin $1-|R_g|$ is ample, and $\varrho_l$ is small; when $\nu$ decreases the spectrum moves toward the imaginary axis, $|R_g|\to1$, the denominator $\to0$, and even with a bounded numerator $\varrho_l$ is pushed toward $1$. This is the experimental confirmation of the "the denominator is the key" insight of Theorem 4.2.

![Source Figure 4.5: degradation of Parareal on advection–diffusion and Burgers' equations as viscosity is lowered](assets/papers/time-parallelization/source-figures/figure-4-5.svg)

Figure 4.5(a) agrees with the spectral-factor prediction. Burgers' equation lacks an equally precise modal analysis, but Figure 4.5(b) exhibits the same trend; around $\nu\le10^{-3}$ the standard iteration diverges. The finite-step property guaranteed by the strictly lower triangular structure still holds, only the number of iterations needed loses practical value. The wave equation also typically does not converge (already noted by Gander and Vandewalle 2007; for a more detailed analysis see Gander and Lunet 2020a,b, Gander, Lunet and Pogoželskytė 2023a, Gander, Lunet, Ruprecht and Speck 2023b), which is precisely the boundary of the applicable range of the Chapter 4 methods.

For the failure mechanism on hyperbolic/propagation problems, the paper explains using the picture of Figure 2.4 in Chapter 2: **an arbitrarily small high-frequency component can propagate arbitrarily far in space and time**. Therefore it is extremely difficult to make the cheap coarse propagator $\mathcal G$ approximate the accuracy of the fine propagator $\mathcal F$ in both space and time; and once $\mathcal G$ is made very accurate for the sake of precision, the coarse correction itself becomes expensive, and Parareal loses its acceleration purpose.

> [!tip] Insight
> To state "why Parareal is excellent for diffusion and fails for propagation" in one breath: **in parabolic problems high frequencies are quickly dissipated and only a few slow modes survive long-term**, so the coarse propagator only needs to be accurate enough on these retained slow modes (low-frequency end $R_g\approx R_f^J$), while for already-decayed high modes the strong dissipation of the coarse level "zeroes them out in place"—both numerator and denominator favor convergence, giving the window-length-independent $\approx0.3$. **In propagation/hyperbolic problems the modes neither decay nor stop carrying phase across space and time**, so the coarse propagator's **phase (dispersion) error** does not vanish over time but accumulates iteration by iteration; here $|R_g|\to1$ makes the denominator of (4.5b) tend to zero, and $\varrho_l\to1$ or even diverges. "Phase-optimized" coarse propagators such as semi-Lagrangian (see below) are the mainstream remedy precisely because they target phase rather than amplitude. This agrees with Ruprecht's (2018) analysis of Parareal's wave-propagation properties: the mismatch of coarse and fine phase speeds is the root cause of convergence failure.

In the MGRiT community (MGRiT is the multilevel generalization of Parareal, see Section 4.4), much work is devoted to making MGRiT handle advection equations (Howse et al. 2019; De Sterck et al. 2021; De Sterck, Falgout, Krzysik and Schroder 2023b; De Sterck, Falgout and Krzysik 2023a and references therein). The core idea is to design a phase-optimized coarse propagator using a **semi-Lagrangian discretization**: tracking information along characteristics and directly aligning phase speeds, which therefore performs well on the **linear** advection equation; but semi-Lagrangian is a characteristic-based method, and in nonlinear problems characteristics interact and are hard to implement directly, so the nonlinear case remains an **open problem**. Another line, proposed by Gander and Wu (2020), is also aimed at hyperbolic problems and has the advantage of handling nonlinearity relatively easily; it is developed in Section 4.5.

## Equation, theorem, and figure coverage audit

| Source item                          | Paper section | Coverage                                                                                                                         |
| ------------------------------------ | ------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| Section 4 introduction and 4.1       | 4, 4.1        | parabolic locality, four method families, historical lines, virtual-control origin                                               |
| (4.1), Figure 4.1                    | 4.2           | update formula, prediction–correction, parallel fine solves, sequential coarse correction, two-level grid, nonintrusiveness      |
| (4.2)–(4.4), Theorem 4.1, Remark 4.1 | 4.2           | modal reduction, Toeplitz error matrix, complete proof, preconditioner/residual-parallel interpretation, finite-step termination |
| (4.5), Theorem 4.2, Figure 4.2       | 4.2           | Lemma 4.4 norm bound, superlinear and long-time linear, dissipation-margin denominator, two experimental regimes                 |
| (4.6), Theorem 4.3                   | 4.2           | nonlinear assumptions, generating functions, three-factor reading of the error bound and finite-step meaning                     |
| (4.7)–(4.9), Theorem 4.4, Figure 4.3 | 4.2           | L/A-stable fine propagators, 0.3 and 0.068 factors, $J_{\min}$, SDIRK comparison                                                 |
| Figures 4.4–4.5                      | 4.2           | advection–diffusion spectrum, Burgers experiment, weak-diffusion/hyperbolic failure mechanism, semi-Lagrangian remedy            |

## Source of this page

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 4–4.2, pp. 443–452.
