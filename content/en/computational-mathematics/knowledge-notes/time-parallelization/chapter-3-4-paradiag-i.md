---
title: "3.5.1: Direct ParaDiag (ParaDiag-I)"
description: A complete derivation spanning the all-at-once system, the truncation/roundoff balance of geometric time grids, and BVM elimination through nonlinear quasi-Newton iteration (including NKA)
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i
tags:
  - parallel-in-time
  - ParaDiag
  - all-at-once systems
---

> [!note] Reading scope
> This page follows Sections 3.5–3.5.1 (pp. 415–431) of the paper. It covers the overall ParaDiag classification, equations (3.22)–(3.48), Theorems 3.5–3.7, Figures 3.9–3.14, and Tables 3.1–3.2. The direct method's truncation error, roundoff error, geometric grid, BVM, and nonlinear quasi-Newton branches are all developed step by step, with an emphasis on reconstructing the "discretization error vs. roundoff error" design trade-off that runs throughout the section.

## 3.5.1 Direct ParaDiag methods (ParaDiag-I)

ParaDiag is the last class of methods introduced in this chapter. Its core idea is to **diagonalize the time-stepping matrix** (or an approximation of it). Parallelism follows because, once the time direction is diagonalized, the coupling between time levels is removed and each level degenerates into an independent system that involves only the spatial operator. Depending on how the time matrix is treated, ParaDiag splits into two branches.

### The boundary between ParaDiag-I and ParaDiag-II

ParaDiag-I (Maday and Rønquist 2008) diagonalizes the time-stepping matrix **exactly**, so like ParaExp it is a direct (non-iterative) parallel-in-time solver. Diagonalization requires one hard prerequisite: the time matrix must be diagonalizable. For the most naive uniform-step backward Euler, the time matrix is a defective Jordan block and cannot be diagonalized at all, so ParaDiag-I must **either use non-uniform time steps** **or switch to a different time-integration formula at the last step** (such as a boundary-value method). The price is this: the error analysis for non-uniform steps (Gander et al. 2019) shows that the number of steps that can be processed in parallel within a single time window is limited by roundoff, usually only about twenty steps in double precision; moreover, systematic results currently exist only for low-order integrators such as backward Euler and the trapezoidal rule, and are hard to extend to higher-order methods such as Runge–Kutta. Switching to a boundary-value-method-type discretization can significantly relax the step-count ceiling, but the usable integrators are still restricted to backward Euler and the trapezoidal rule (Liu et al. 2022).

ParaDiag-II (Gander et al. 2021c) instead **approximates** the time matrix so that it becomes well-conditioned and diagonalizable, and then places this approximation inside a stationary iteration or uses it as a preconditioner for a Krylov method, thus "trading directness for iteration." Its preconditioner design follows two principles: first, the diagonalization must be well-conditioned (the eigenvector matrix must have a small condition number) so as to suppress the roundoff error introduced by the diagonalization step—precisely what ParaDiag-I lacks; second, the iteration must converge fast (small spectral radius, or a preconditioned spectrum tightly clustered near 1 to benefit Krylov acceleration). This iterative ParaDiag was first proposed at the discrete level by McDonald, Pestana and Wathen (2018) and independently given at the continuous level by Gander and Wu (2019). A full discussion of ParaDiag-II is left to the next page.

> [!tip] Insight
> The split between the two branches is essentially an **accuracy–robustness** trade-off: ParaDiag-I insists on an exact factorization, buying a clean directness in which "there is no outer iteration and the error is just the discretization error," but it concentrates all the risk on the condition number of the eigenvector matrix $V$; once $\operatorname{Cond}(V)$ explodes as the window grows, directness turns into a shackle. ParaDiag-II deliberately gives up exactness, trading a well-conditioned approximation for arbitrarily long windows and arbitrary integrators, at the cost of only a controllable handful of iterations. Every technical detail that follows—geometric grids, BVM, and the choice of $\varrho$—is fighting to secure as large a usable step count as possible for ParaDiag-I's "exact but fragile" route.

### Backward Euler as an all-at-once system

Apply variable-step backward Euler to $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$:

$$
\frac{\boldsymbol u_n-\boldsymbol u_{n-1}}{\Delta t_n}
=A\boldsymbol u_n+\boldsymbol g_n,
\qquad n=1,\ldots,N_t. \tag{3.22}
$$

The traditional approach marches level by level in $n$ and is inherently serial. ParaDiag's first shift is to **abandon level-by-level solving** and stack the $N_t$ difference equations into a single "all-at-once" system, solving all time levels at once. With $\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top$, this gives

$$
K\boldsymbol U=\boldsymbol b,
\qquad
K=B\otimes I_x-I_t\otimes A, \tag{3.23a}
$$

where

$$
B=
\begin{bmatrix}
\Delta t_1^{-1}\\
-\Delta t_2^{-1}&\Delta t_2^{-1}\\
&\ddots&\ddots\\
&&-\Delta t_{N_t}^{-1}&\Delta t_{N_t}^{-1}
\end{bmatrix},
\qquad
\boldsymbol b=
\begin{bmatrix}
\Delta t_1^{-1}\boldsymbol u_0+\boldsymbol g_1\\
\boldsymbol g_2\\
\vdots\\
\boldsymbol g_{N_t}
\end{bmatrix}. \tag{3.23b}
$$

Here $I_x\in\mathbb R^{N_x\times N_x}$ and $I_t\in\mathbb R^{N_t\times N_t}$ are identity matrices. The Kronecker structure captures the "time $\times$ space" tensor decomposition exactly: $B\otimes I_x$ acts only in the time direction ($B$ is the backward-Euler time-difference operator, lower bidiagonal, with diagonal $1/\Delta t_n$ and subdiagonal $-1/\Delta t_n$), while $I_t\otimes A$ acts only in the space direction. The initial value $\boldsymbol u_0$ appears only in the right-hand side of the first block, showing that a pure initial-value problem has been flattened into a system coupled across the spatial dimension.

The reason diagonalization is possible is that **only the time-direction matrix $B$ needs to be treated**. If all $\Delta t_n$ are distinct, the diagonal entries $1/\Delta t_n$ of $B$ are pairwise distinct, so as a triangular matrix $B$ has $N_t$ distinct eigenvalues and is therefore diagonalizable:

$$
B=VDV^{-1},
\qquad
D=\operatorname{diag}\!\left(
\frac1{\Delta t_1},\ldots,\frac1{\Delta t_{N_t}}
\right). \tag{3.24}
$$

Substituting $B=VDV^{-1}$ into $K$ and using the mixed-product property of the Kronecker product $(V\otimes I_x)(D\otimes I_x)(V^{-1}\otimes I_x)=VDV^{-1}\otimes I_x=B\otimes I_x$, together with $I_t\otimes A=(V\otimes I_x)(I_t\otimes A)(V^{-1}\otimes I_x)$ (because $V I_t V^{-1}=I_t$), gives the block factorization

$$
K=(V\otimes I_x)
(D\otimes I_x-I_t\otimes A)
(V^{-1}\otimes I_x).
$$

The middle factor $D\otimes I_x-I_t\otimes A$ is **block diagonal**: its $n$-th diagonal block is exactly $\tfrac1{\Delta t_n}I_x-A$, with no remaining coupling between time levels. Solving $K\boldsymbol U=\boldsymbol b$ therefore splits into three stages:

$$
\left\{
\begin{aligned}
\boldsymbol U^a&=(V^{-1}\otimes I_x)\boldsymbol b,
&&\text{inverse time transform},\\
\left(\frac1{\Delta t_n}I_x-A\right)\boldsymbol u_n^b
&=\boldsymbol u_n^a,
&&n=1,\ldots,N_t,\\
\boldsymbol U&=(V\otimes I_x)\boldsymbol U^b,
&&\text{forward time transform}.
\end{aligned}
\right. \tag{3.25}
$$

Here $\boldsymbol U^a=((\boldsymbol u_1^a)^\top,\ldots,(\boldsymbol u_{N_t}^a)^\top)^\top$, and similarly for $\boldsymbol U^b$. The first stage (step a) and the third stage (step c) merely apply $V^{-1}$ and $V$ to the data along the time direction, i.e. one dense $N_t\times N_t$ matrix multiplication per spatial degree of freedom, which is relatively cheap. The real workload is the second stage (step b): the $N_t$ shifted spatial systems $\bigl(\tfrac1{\Delta t_n}I_x-A\bigr)\boldsymbol u_n^b=\boldsymbol u_n^a$. **They are completely independent of one another** and can be distributed to $N_t$ processors and solved in parallel—this is the entire mechanism behind "diagonalization = parallel in time."

![The time transform, independent spatial solves, and inverse transform in ParaDiag](assets/diagrams/pint/en/paradiag-three-stage.svg)

> [!tip] Insight
> The three-stage method can be viewed as a "spectral transform—diagonal solve—inverse transform" along the time direction, exactly like solving a constant-coefficient circular convolution with the FFT: step a/step c are the transform and its inverse, and step b is a pointwise (eigenvalue-by-eigenvalue) solve in the "frequency" domain. The difference is that a circulant matrix uses the unitary Fourier basis (whose condition number is always 1), whereas here $V$ is a general triangular Toeplitz eigenvector matrix whose **condition number can be enormous**. All of ParaDiag-I's risk is therefore concentrated in step a/step c: an ill-conditioned $V$ re-pollutes the clean step-b solution back into the time domain. This is why every effort that follows revolves around "how to keep $V$ as well-conditioned as possible."

### Geometric grids and the two competing errors

For arbitrary steps $\{\Delta t_n\}$, $V$ generally has no closed form and can only be obtained by a numerical `eig`. This is not expensive in practice ($N_t$ need not be large), but it leaves the full roundoff analysis and parameter selection with no foothold. To guarantee distinct eigenvalues while permitting analytic treatment, Maday and Rønquist (2008) adopt a **geometric grid** $\Delta t_n=\mu^{n-1}\Delta t_1$ with $\mu>1$. On a one-dimensional heat equation they take $\mu=1.2$ and obtain an almost perfect speedup. From $\sum_n\Delta t_n=\sum_n\mu^{n-1}\Delta t_1=T$ the first step is fixed as $\Delta t_1=T/\sum_{j=1}^{N_t}\mu^{j-1}$, so that

$$
\Delta t_n=\frac{\mu^{n-1}}{\sum_{j=1}^{N_t}\mu^{j-1}}T. \tag{3.26}
$$

The geometric grid turns "diagonalizability" into a continuously tunable question. Writing $\mu=1+\varrho$, the quantity $\varrho$ directly measures how far the grid deviates from uniform, and each extreme has its own failure mechanism:

- $\varrho$ **too large**: the final step $\Delta t_{N_t}=\mu^{N_t-1}\Delta t_1$ swells exponentially relative to the first step, the grid becomes severely non-uniform, and the **truncation (discretization) error** of the later time levels rises sharply.
- $\varrho$ **too small** ($\mu\to1$): all $\Delta t_n$ approach equality, the eigenvalues $1/\Delta t_n$ of $B$ crowd together, and $B$ approaches a single Jordan block. A Jordan block is a defective matrix whose eigenvector matrix $V$ has a diverging condition number, so during diagonalization the **roundoff error** is amplified by $\operatorname{Cond}(V)$.

Parameter selection in direct ParaDiag is essentially a balance between these two opposing errors.

> [!tip] Insight
> There is an intriguing counterintuitive point hidden here: numerical analysis usually encourages **uniform refinement** to reduce truncation error, but in ParaDiag-I the uniform grid is exactly the worst choice—it degenerates $B$ into a non-diagonalizable Jordan block and the method collapses outright. The geometric grid is a "grid distortion deliberately introduced for the sake of diagonalizability": the larger $\varrho$ is, the farther $B$ is from a Jordan block (the better conditioned $V$ is), but the more it departs from equal spacing (the larger the truncation error). Thus $\varrho$ is not an ordinary discretization parameter in the "smaller is more accurate" sense, but a robustness parameter that must be **propped up from below**—which is exactly what makes it intuitively fundamentally different from a conventional step-size parameter.

### Theorem 3.5: balance formula for a first-order problem

The trade-off above was made precise for first-order parabolic problems in Gander et al. (2016a). Assume $\sigma(A)\subset\mathbb R_-$ and $|\lambda(A)|\leq\lambda_{\max}$. Let $\boldsymbol u_{N_t}(\varrho)$ and $\boldsymbol u_{N_t}(0)$ denote the geometric-grid and uniform-grid backward-Euler solutions at $t=T$, and let $\widetilde{\boldsymbol u}_n(\varrho)$ denote the value actually computed by the diagonalization method (3.25). Then

$$
\left\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u_{N_t}(0)\right\|
\lesssim C(\lambda_*T,N_t)\varrho^2,
$$

$$
\left\|\widetilde{\boldsymbol u}_n(\varrho)-\boldsymbol u_n(\varrho)\right\|
\lesssim
\epsilon\,
\frac{N_t^2(2N_t+1)(N_t+\lambda_{\max}T)}{\phi(N_t)}
\varrho^{-(N_t-1)}. \tag{3.27}
$$

Here

$$
C(x,N_t)=\frac{N_t(N_t^2-1)}{24}r(x/N_t,N_t),
\qquad
r(\widetilde x,N_t)=
\left(\frac{\widetilde x}{1+\widetilde x}\right)^2
(1+\widetilde x)^{-N_t},
$$

$$
\phi(N_t)=
\begin{cases}
\left(\dfrac{N_t}{2}\right)!\left(\dfrac{N_t}{2}-1\right)!,&N_t\ \text{even},\\[6pt]
\left[\left(\dfrac{N_t-1}{2}\right)!\right]^2,&N_t\ \text{odd},
\end{cases}
$$

$\widetilde x_*$ maximizes $r(\widetilde x,N_t)$ on $[0,\infty)$, and $\lambda_*=N_t\widetilde x_*/T$. The two estimates are exactly the two error classes described above: **the first is proportional to $\varrho^2$**, the extra truncation error the geometric grid pays relative to the uniform grid (increasing with $\varrho$); **the second is proportional to $\varrho^{-(N_t-1)}$**, the roundoff error introduced by diagonalization (exploding as $\varrho$ decreases), and it carries the machine epsilon $\epsilon$ explicitly. Setting the two equal, i.e. solving $C\varrho^2=\epsilon(\cdots)\varrho^{-(N_t-1)}$, gives the optimal parameter

$$
\varrho_{\mathrm{opt}}=
\left(
\epsilon
\frac{N_t^2(2N_t+1)(N_t+\lambda_{\max}T)}
{\phi(N_t)C(\lambda_*T,N_t)}
\right)^{1/(N_t+1)}. \tag{3.28}
$$

The idea of the proof is to **decouple the coupled system eigenvalue by eigenvalue**: for $\lambda\in\sigma(A)$, the all-at-once system degenerates along that eigen-direction into the Dahlquist test equation $y'=\lambda y$. The first estimate comes from Gander et al. (2016a, Theorem 2), comparing the discretization error of the geometric and uniform grids for that scalar equation; the second comes from Theorem 6 of the same paper, estimating the roundoff introduced by the eigenvector transforms, with the worst case attained at $|\lambda|=\lambda_{\max}$ (the mode with the largest spectral radius is the most sensitive). The first estimate also yields the truncation error to the exact solution: by the triangle inequality
$\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u(T)\|\leq\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u_{N_t}(0)\|+\|\boldsymbol u_{N_t}(0)-\boldsymbol u(T)\|$,
where the last term is the error of standard uniform backward Euler, which is well understood and not dominant, so the extra cost of the geometric grid is captured entirely by the $\varrho^2$ term.

On a geometric grid, $V$ and $V^{-1}$ have closed forms, both lower-triangular Toeplitz matrices (Gander et al. 2016a):

$$
V=\mathbb T(p_1,\ldots,p_{N_t-1}),
\qquad
p_n=\frac1{\prod_{j=1}^{n}(1-\varrho^j)},
$$

$$
V^{-1}=\mathbb T(q_1,\ldots,q_{N_t-1}),
\qquad
q_n=(-1)^n\varrho^{n(n-1)/2}p_n, \tag{3.29a}
$$

where the lower-triangular Toeplitz operator is

$$
\mathbb T(a_1,\ldots,a_{N_t-1})=
\begin{bmatrix}
1\\
a_1&1\\
\vdots&\ddots&\ddots\\
a_{N_t-1}&\cdots&a_1&1
\end{bmatrix}. \tag{3.29b}
$$

The value of these closed forms lies not in the actual solve but in **analytically estimating $\operatorname{Cond}(V)$**, so that the roundoff constant in the second estimate of (3.27) can be written as an explicit function of $N_t$ (the factorials in $\phi(N_t)$ come precisely from the products defining $p_n$ and $q_n$). When (3.25) is actually implemented, these formulas are not used; instead one calls `eig`, which automatically rescales the eigenvectors to improve the condition number, obtaining a better practical $\operatorname{Cond}(V)$ than the theoretical closed form. The single- and double-precision machine epsilons are approximately $1.19\times10^{-7}$ and $2.22\times10^{-16}$, respectively.

> [!tip] Insight
> In (3.28), the exponent $1/(N_t+1)$ in $\varrho_{\mathrm{opt}}\sim\epsilon^{1/(N_t+1)}$ reveals the mathematical root of ParaDiag-I's "about-twenty-step ceiling." On one hand, as $N_t$ grows, $\epsilon^{1/(N_t+1)}\to1$, forcing the optimal $\varrho$ to increase and the grid to become ever more distorted; on the other hand, the roundoff constant contains $\phi(N_t)$ with **factorial** growth, and $\varrho^{-(N_t-1)}$ is an exponential amplification. The superposition of these two forces makes the window width that can be pushed below machine precision only of order $O(\log(1/\epsilon))$—which in double precision lands right around twenty steps. In other words, the ceiling is not an implementation flaw but the intrinsic information-theoretic limit of the strategy "trading time parallelism for an ill-conditioned similarity transform."

### Figures 3.9–3.10: even the optimal parameter has a parallel-width ceiling

![Original Figure 3.9: error versus the geometric-grid parameter for the heat and advection–diffusion equations](assets/papers/time-parallelization/source-figures/figure-3-9.svg)

Figure 3.9 uses homogeneous Dirichlet boundaries, $u_0(x)=\sin(2\pi x)$, $\Delta x=1/50$, and $T=0.2$. The advection–diffusion viscosity in the right panel is $10^{-2}$. Five values of $N_t$ each scan $\varrho\in[10^{-2},1]$, the error is the maximum $L^\infty$ error over all time nodes, and the reference solution is given by the exponential integrator $\boldsymbol u(t_n)=e^{-At_n}\boldsymbol u_0$. Every curve shows a clear V shape—the left branch is roundoff-dominated, the right branch is truncation-dominated, the trough is the optimal $\varrho$, and the star is the theoretical prediction of (3.28). The prediction is accurate for advection–diffusion and slightly off for the heat equation at small $N_t$, consistent with (3.28) being only an **asymptotic** balance of the two estimates.

![Original Figure 3.10: error first falls and then rises with the time-step count at the optimal geometric parameter](assets/papers/time-parallelization/source-figures/figure-3-10.svg)

Figure 3.10 instead takes $T=0.5$ and $N_t=2^4,2^5,\ldots,2^{10}$, using the numerically optimal $\varrho_{\mathrm{num}}$ for each $N_t$. The error of uniform backward Euler continues to decrease monotonically with $N_t$ (finer time steps are more accurate); the error of ParaDiag-I, by contrast, first falls and then rises, crossing a threshold at fewer than 100 steps and deteriorating rapidly. This non-monotone curve is the direct experimental manifestation of the ceiling described above: once the roundoff error amplified by $\operatorname{Cond}(V)$ overtakes the decrease in time-discretization error, adding steps becomes harmful.

### Wave equations and the trapezoidal rule

The same framework extends to second-order systems. Consider the spatially discretized wave equation

$$
\boldsymbol u''(t)=A\boldsymbol u(t),
\quad
\boldsymbol u(0)=\boldsymbol u_0,
\quad
\boldsymbol u'(0)=\widetilde{\boldsymbol u}_0 \tag{3.30}
$$

(in the wave equation $A$ is the discrete Laplace operator). To apply ParaDiag-I, first reduce to a first-order system

$$
\boldsymbol w'=\mathbb A\boldsymbol w,
\qquad
\boldsymbol w=(\boldsymbol u^\top,(\boldsymbol u')^\top)^\top,
\qquad
\mathbb A=
\begin{bmatrix}0&I_x\\A&0\end{bmatrix}. \tag{3.31}
$$

Backward Euler has excessive numerical dissipation on wave problems, destroying energy conservation and introducing dispersive artifacts. We therefore switch to the **trapezoidal rule**

$$
\frac{\boldsymbol w_n-\boldsymbol w_{n-1}}{\Delta t_n}
=\frac{\mathbb A}{2}(\boldsymbol w_n+\boldsymbol w_{n-1}), \tag{3.32}
$$

which is energy-conserving ($\|\boldsymbol w_n\|^2=\|\boldsymbol w_0\|^2$) and hence suitable for long-time wave propagation. The steps are still taken from the geometric grid (3.26). The two-point average of the trapezoidal rule makes the right-hand side contain both $\boldsymbol w_n$ and $\boldsymbol w_{n-1}$, so besides the time-difference matrix $B$ the all-at-once system also acquires a **time-averaging matrix** $\widetilde B$:

$$
K\boldsymbol W=\boldsymbol b,
\qquad
K=B\otimes I_x-\widetilde B\otimes\mathbb A, \tag{3.33a}
$$

$$
\widetilde B=\frac12
\begin{bmatrix}
1\\
1&1\\
&\ddots&\ddots\\
&&1&1
\end{bmatrix}. \tag{3.33b}
$$

Now $K$ contains two different time matrices $B$ and $\widetilde B$, so single-matrix diagonalization cannot be applied directly. The trick is to **left-multiply by $\widetilde B^{-1}\otimes I_x$**, normalizing $\widetilde B$ to the identity and restoring the standard form of (3.23a):

$$
\mathcal K\boldsymbol W=\widetilde{\boldsymbol b},
\qquad
\mathcal K=\widetilde B^{-1}B\otimes I_x-I_t\otimes\mathbb A,
\qquad
\widetilde{\boldsymbol b}=(\widetilde B^{-1}\otimes I_x)\boldsymbol b. \tag{3.34}
$$

Now only the single time matrix $\widetilde B^{-1}B$ needs to be diagonalized. Gander et al. (2019) give its closed-form diagonalization:

$$
\widetilde B^{-1}B
=V\operatorname{diag}\!\left(
\frac2{\Delta t_1},\ldots,\frac2{\Delta t_{N_t}}
\right)V^{-1}, \tag{3.35a}
$$

$$
\begin{aligned}
V&=\mathbb T(p_1,\ldots,p_{N_t-1}),
&p_n&=\prod_{j=1}^{n}\frac{1+\mu^j}{1-\mu^j},\\
V^{-1}&=\mathbb T(q_1,\ldots,q_{N_t-1}),
&q_n&=\mu^{-n}\prod_{j=1}^{n}
\frac{1+\mu^{-j+2}}{1-\mu^{-j}}.
\end{aligned} \tag{3.35b}
$$

The eigenvalues $2/\Delta t_n$ are still distinct (guaranteed by the geometric grid), $V$ is still lower-triangular Toeplitz, and so (3.34) can still be solved by the three stages of (3.25).

### Theorem 3.6: balance formula for a wave problem

For $\lambda(A)\leq0$, the geometric-grid and uniform-grid trapezoidal solutions, together with the diagonalization-computed value, satisfy

$$
\left\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u_{N_t}(0)\right\|
\lesssim\frac{N_t(N_t^2-1)}{15}\varrho^2,
$$

$$
\left\|\widetilde{\boldsymbol u}_n(\varrho)-\boldsymbol u_n(\varrho)\right\|
\lesssim
\epsilon\frac{2^{2N_t-1/2}N_t}{(N_t-1)!}
\varrho^{-(N_t-1)}. \tag{3.36}
$$

The structure is identical to the first-order case—a truncation term $\propto\varrho^2$ and a roundoff term $\propto\varrho^{-(N_t-1)}$—and setting them equal gives

$$
\varrho_{\mathrm{opt}}=
\left(
\epsilon\frac{15\times2^{2N_t-1/2}}
{(N_t^2-1)(N_t-1)!}
\right)^{1/(N_t+1)}. \tag{3.37}
$$

The proof again proceeds eigenvalue by eigenvalue: for each eigenvalue $\lambda>0$ of $-A$, consider the scalar oscillator $u''+\lambda u=0$. By Gander et al. (2019, Theorem 2.1), the difference between the geometric- and uniform-grid truncation errors is $O\!\bigl(\tfrac{N_t(N_t^2-1)}{6}r_1(\tfrac{\lambda T}{2N_t})\varrho^2\bigr)$, where

$$
r_1(s)=\frac{s^3}{(1+s^2)^2},
\qquad r_1(s)\leq\frac25\ (s\geq0),
$$

and substituting gives the first estimate in (3.36) (since $\tfrac{N_t(N_t^2-1)}{6}\cdot\tfrac25=\tfrac{N_t(N_t^2-1)}{15}$). The roundoff error is given by Theorem 2.11 of the same paper and contains the factor $r_2(s)=1/(1+s^2)\leq1$; a uniform bound gives the second estimate. The two **uniform bounds** $r_1\leq2/5$ and $r_2\leq1$ eliminate the dependence on the specific eigenvalue $\lambda$, leaving the bound depending only on $N_t,\varrho,\epsilon$.

![Original Figure 3.11: optimal parameters and step-count threshold for geometric-trapezoidal ParaDiag-I on the wave equation](assets/papers/time-parallelization/source-figures/figure-3-11.svg)

Figure 3.11 uses homogeneous Dirichlet boundaries, $\Delta x=1/20$, and $T=0.2$. In panel (a), each $N_t$ has a minimum-error point, and the star from (3.37) is close to the measured optimum; in panel (b), with the numerically optimal parameter, the error deteriorates rapidly for $N_t>32$, consistent with the ceiling phenomenon of the first-order case.

![Original Table 3.1: eigenvector condition numbers for the backward-Euler and trapezoidal time matrices](assets/papers/time-parallelization/source-figures/table-3-1.svg)

Table 3.1 measures the culprit behind the failure directly—$\operatorname{Cond}(V)$. For $N_t=5,10,20,30,60,100$, the $B$ of backward Euler grows from $1.7\times10^3$ to $4.8\times10^6$, and the $\widetilde B^{-1}B$ of the trapezoidal rule grows from $4.7\times10^3$ to $4.1\times10^9$. The condition number rises rapidly with $N_t$, in agreement with the roundoff analysis of (3.27) and (3.36), and it explains the error deterioration of Figures 3.10–3.11. An interesting observation is that when the **numerically optimal** $\varrho_{\mathrm{num}}$ is used (rather than the theoretical $\varrho_{\mathrm{opt}}$), $\operatorname{Cond}(V)$ exhibits a **plateau** in the later range rather than continual explosion—showing that the automatic rescaling of `eig` suppresses the condition number better in practice than the closed-form analysis, a phenomenon the authors note deserves further study.

### BVM: fix the step size and change the terminal formula

The dilemma of the geometric grid is that diagonalizability requires distorting the grid, and the distortion itself limits the step count. Liu et al. (2022) take a different route—**keep the step size uniform** and instead let the **last step use a different formula**, thereby sidestepping the Jordan block. The first $N_t-1$ steps use a second-order centered difference, and the last step uses backward Euler:

$$
\left\{
\begin{aligned}
\frac{\boldsymbol u_{n+1}-\boldsymbol u_{n-1}}{2\Delta t}
&=A\boldsymbol u_n+\boldsymbol g_n,
&&n=1,\ldots,N_t-1,\\
\frac{\boldsymbol u_{N_t}-\boldsymbol u_{N_t-1}}{\Delta t}
&=A\boldsymbol u_{N_t}+\boldsymbol g_{N_t}.
\end{aligned}
\right. \tag{3.38}
$$

The centered difference is unstable in time-stepping mode (it uses both $\boldsymbol u_{n-1}$ and $\boldsymbol u_{n+1}$ and cannot march level by level); it works here because the whole system is **solved simultaneously**, making it a **boundary-value method (BVM)**. Its stability cannot be judged like an ordinary step-by-step scheme—the first-order backward-Euler "cap" at the last step supplies the missing boundary condition for the entire centered-difference skeleton (see Gander 2015, Section 5.2 for the stability discussion). Axelsson and Verwer (1985) used exactly this boundary-value technique to circumvent the famous **Dahlquist convergence–stability barrier**: they proved that even with only a first-order last step, the numerical solution obtained by solving simultaneously is still **uniformly second-order accurate** in the general nonlinear case (Axelsson and Verwer 1985, Theorem 4). Earlier, Fox (1954) and Fox and Mitchell (1957) had already used such discretizations, only with BDF2 at the last step:

$$
\frac{3\boldsymbol u_{N_t}-4\boldsymbol u_{N_t-1}+\boldsymbol u_{N_t-2}}{2\Delta t}=A\boldsymbol u_{N_t}+\boldsymbol g_{N_t}.
$$

Method (3.38) is the prototypical example of what later became known as BVM, and the well-posedness of its all-at-once system was rigorously established by Brugnano, Mazzia and Trigiante (1993) (see also Brugnano and Trigiante 2003); in the BVM tradition such all-at-once systems are usually solved by constructing effective preconditioned iterations, whereas ParaDiag-I solves them directly by diagonalization.

The all-at-once form is still $K\boldsymbol U=\boldsymbol b$, $K=B\otimes I_x-I_t\otimes A$, only with $B$ replaced by the centered-difference skeleton plus a backward-Euler last row:

$$
K\boldsymbol U=\boldsymbol b,
\qquad
K=B\otimes I_x-I_t\otimes A, \tag{3.39a}
$$

$$
B=\frac1{\Delta t}
\begin{bmatrix}
0&\tfrac12\\
-\tfrac12&0&\tfrac12\\
&\ddots&\ddots&\ddots\\
&&-\tfrac12&0&\tfrac12\\
&&&-1&1
\end{bmatrix},
\qquad
\boldsymbol b=
\begin{bmatrix}
\boldsymbol u_0/(2\Delta t)+\boldsymbol g_1\\
\boldsymbol g_2\\
\vdots\\
\boldsymbol g_{N_t}
\end{bmatrix}. \tag{3.39b}
$$

Since the step size is uniform, only the initial value $\boldsymbol u_0$ is needed, and all time levels are solved at once. The key is the well-conditioning of its diagonalization:

**Theorem 3.7.** The matrix $B$ in (3.39b) admits a factorization $B=VDV^{-1}$ with $\operatorname{Cond}(V)=O(N_t^2)$ (the closed forms of $V$, $V^{-1}$, and $D$ are given in Liu et al. 2022, Section 3).

This is BVM's decisive advantage over the geometric grid: the geometric grid's $\operatorname{Cond}(V)$ deteriorates **exponentially** with $N_t$ (Table 3.1), whereas BVM grows only **polynomially** ($O(N_t^2)$), so the number of steps that can be parallelized in a single window is greatly increased.

> [!tip] Insight
> The geometric grid and BVM are two complementary routes to "diagonalizability": the geometric grid keeps **the same integrator** and separates the eigenvalues by **perturbing the step sizes**, at the cost of exponentially ill-conditioned eigenvectors; BVM keeps **uniform steps** and separates the eigenvalues by **changing the integrator in one row**, at the cost of only a lower-order last step (compensated by the Axelsson–Verwer uniform-second-order result). The former puts the distortion on the grid, the latter on the scheme—and putting the distortion in the terminal formula is far gentler than putting it across the whole grid, which is the essential source of the $O(N_t^2)$-versus-exponential improvement.

A second-order system can likewise be reduced to first order and then solved by BVM. Set $\boldsymbol v=\boldsymbol u'$ and $\boldsymbol w=(\boldsymbol u^\top,\boldsymbol v^\top)^\top$, and apply the same BVM to $\boldsymbol w$:

$$
\left\{
\begin{aligned}
\frac{\boldsymbol w_{n+1}-\boldsymbol w_{n-1}}{2\Delta t}
&=\mathbb A\boldsymbol w_n,
&&n=1,\ldots,N_t-1,\\
\frac{\boldsymbol w_{N_t}-\boldsymbol w_{N_t-1}}{\Delta t}
&=\mathbb A\boldsymbol w_{N_t}.
\end{aligned}
\right. \tag{3.40}
$$

But introducing the velocity variable $\boldsymbol v$ doubles the spatial storage at each time point, a considerable cost for high dimensions or fine meshes. We can therefore **eliminate $\boldsymbol v$** and retain only $\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top$. To see this, write (3.40) at the discrete level separately in terms of $\{\boldsymbol u_n\}$ and $\{\boldsymbol v_n\}$ (the upper half-row of $\mathbb A$ gives $\boldsymbol u'=\boldsymbol v$, the lower half-row gives $\boldsymbol v'=A\boldsymbol u$):

$$
\left\{
\begin{aligned}
\frac{\boldsymbol u_{n+1}-\boldsymbol u_{n-1}}{2\Delta t}&=\boldsymbol v_n,\quad n<N_t,\\
\frac{\boldsymbol u_{N_t}-\boldsymbol u_{N_t-1}}{\Delta t}&=\boldsymbol v_{N_t},
\end{aligned}
\right.
\qquad
\left\{
\begin{aligned}
\frac{\boldsymbol v_{n+1}-\boldsymbol v_{n-1}}{2\Delta t}&=A\boldsymbol u_n,\quad n<N_t,\\
\frac{\boldsymbol v_{N_t}-\boldsymbol v_{N_t-1}}{\Delta t}&=A\boldsymbol u_{N_t}.
\end{aligned}
\right.
$$

Using the $B$ of (3.39b), these two groups are exactly

$$
(B\otimes I_x)\boldsymbol U-\boldsymbol V=\boldsymbol b_1,
\qquad
(B\otimes I_x)\boldsymbol V-A\boldsymbol U=\boldsymbol b_2,
$$

where $\boldsymbol b_1=(\tfrac{\boldsymbol u_0^\top}{2\Delta t},0,\ldots,0)^\top$ and $\boldsymbol b_2=(\tfrac{\widetilde{\boldsymbol u}_0^\top}{2\Delta t},0,\ldots,0)^\top$ (the initial values enter only the first row). Solving the first equation for $\boldsymbol V=(B\otimes I_x)\boldsymbol U-\boldsymbol b_1$ and substituting into the second gives $(B\otimes I_x)^2\boldsymbol U-A\boldsymbol U=\boldsymbol b_2+(B\otimes I_x)\boldsymbol b_1$. Using $(B\otimes I_x)^2=B^2\otimes I_x$ and computing the right-hand side directly as $\boldsymbol b_2+(B\otimes I_x)\boldsymbol b_1=\boldsymbol b$, we obtain the all-at-once system in $\boldsymbol U$ only,

$$
(B^2\otimes I_x-I_t\otimes A)\boldsymbol U=\boldsymbol b, \tag{3.41}
$$

$$
\boldsymbol b=
\left(
\frac{\widetilde{\boldsymbol u}_0^\top}{2\Delta t},
-\frac{\boldsymbol u_0^\top}{4\Delta t^2},
0,\ldots,0
\right)^\top.
$$

The time matrix becomes $B^2$, whose eigenvector matrix is the same as that of $B$ (since $B^2=VD^2V^{-1}$), so it still inherits the $O(N_t^2)$ good conditioning of Theorem 3.7 and can be solved in three stages as usual, with no doubling of storage.

![Original Figure 3.12: errors and condition numbers for the geometric time grid and the BVM on the wave equation](assets/papers/time-parallelization/source-figures/figure-3-12.svg)

Figure 3.12 uses $T=0.5$, $\Delta x=1/40$, and homogeneous Dirichlet boundaries to compare geometric-trapezoidal ParaDiag-I with BVM ParaDiag-I ($N_t=2^2,\ldots,2^8$). The geometric-trapezoidal method begins to show the typical roundoff deterioration near $N_t\approx32$; the BVM retains $O(\Delta t^2)$ accuracy with no deterioration, matching serial trapezoidal integration. The right panel gives the corresponding eigenvector condition numbers, with the BVM markedly lower, directly confirming the source of the "no deterioration."

### Nonlinear all-at-once equations and quasi-Newton iteration

Everything above concerns linear problems. For the nonlinear $\boldsymbol u'=\boldsymbol f(\boldsymbol u,t)$ (second-order problems are handled similarly), define

$$
F(\boldsymbol U)=
\left(
\boldsymbol f(\boldsymbol u_1,t_1)^\top,
\ldots,
\boldsymbol f(\boldsymbol u_{N_t},t_{N_t})^\top
\right)^\top,
$$

and the all-at-once equation is

$$
(B\otimes I_x)\boldsymbol U-F(\boldsymbol U)=\boldsymbol b, \tag{3.42}
$$

where $B$ can be taken either from the geometric grid (3.23b) or from the BVM (3.39b). Applying Newton's method to the nonlinear system, the raw update is
$(B\otimes I_x-\nabla F(\boldsymbol U^k))(\boldsymbol U^{k+1}-\boldsymbol U^k)=\boldsymbol b-((B\otimes I_x)\boldsymbol U^k-F(\boldsymbol U^k))$,
rearranged into a form convenient for diagonalization,

$$
\left(B\otimes I_x-\nabla F(\boldsymbol U^k)\right)\boldsymbol U^{k+1}
=\boldsymbol b-\left(
\nabla F(\boldsymbol U^k)\boldsymbol U^k-F(\boldsymbol U^k)
\right), \tag{3.43a}
$$

$$
\nabla F(\boldsymbol U^k)=
\operatorname{blkdiag}\left(
\nabla f(\boldsymbol u_1^k,t_1),\ldots,
\nabla f(\boldsymbol u_{N_t}^k,t_{N_t})
\right). \tag{3.43b}
$$

The problem is that $\nabla F$ is a block-diagonal matrix that **differs from one time level to another**, so $K=B\otimes I_x-\nabla F$ no longer has the clean Kronecker structure $I_t\otimes(\cdot)$, cannot be diagonalized by a single $V$, and time parallelism is thereby lost. The remedy (inspired by Gander and Halpern 2017) is to replace all blocks by a **single averaged Jacobian** $A_k$:

$$
A_k=\frac1{N_t}\sum_{n=1}^{N_t}\nabla f(\boldsymbol u_n^k,t_n),
\qquad\text{or}\qquad
A_k=\nabla f\!\left(
\frac1{N_t}\sum_{n=1}^{N_t}\boldsymbol u_n^k,
\frac{T}{N_t}
\right). \tag{3.44}
$$

Then $\nabla F(\boldsymbol U^k)\approx I_t\otimes A_k$, and substituting back into (3.43a) gives the quasi-Newton iteration

$$
(B\otimes I_x-I_t\otimes A_k)\boldsymbol U^{k+1}
=\boldsymbol b-\left((I_t\otimes A_k)\boldsymbol U^k-F(\boldsymbol U^k)\right). \tag{3.45}
$$

The Kronecker structure is restored, and after diagonalizing $B=VDV^{-1}$ each round is still solved in three stages:

$$
\left\{
\begin{aligned}
\boldsymbol U^a&=(V^{-1}\otimes I_x)\boldsymbol r^k,\\
(\lambda_n I_x-A_k)\boldsymbol u_n^b&=\boldsymbol u_n^a,
&&n=1,\ldots,N_t,\\
\boldsymbol U^{k+1}&=(V\otimes I_x)\boldsymbol U^b,
\end{aligned}
\right. \tag{3.46}
$$

where $\boldsymbol r^k=\boldsymbol b-((I_t\otimes A_k)\boldsymbol U^k-F(\boldsymbol U^k))$ is the right-hand side of (3.45). In the linear case $A_k=A$ and $\boldsymbol r^k=\boldsymbol b$, and (3.46) reduces to (3.25). The convergence of such quasi-Newton methods based on an approximate Jacobian has been thoroughly studied (Deuflhard 2004, Theorem 2.5; Ortega and Rheinboldt 2000). Their convergence rate depends on how well the single $A_k$ approximates all $N_t$ Jacobian blocks: if $\nabla f(\boldsymbol u_n^k,t_n)$ varies strongly in time, no single matrix can approximate them all at once, in which case one should shorten the time window and process windows serially.

![Original Figure 3.13: convergence of BVM ParaDiag-I for Burgers' equation at two viscosities and several window lengths](assets/papers/time-parallelization/source-figures/figure-3-13.svg)

Figure 3.13 uses the periodic Burgers equation, $\Delta x=0.01$, and keeps $N_t=T/\Delta t=200$. The horizontal line is the space–time discretization error $\max\{\Delta t^2,\Delta x^2\}=10^{-4}$. At $\nu=0.1$, the convergence rate changes little from $T=0.1$ to $1.6$, showing that the Jacobian varies mildly in time and a single $A_k$ is good enough; at $\nu=0.002$, longer windows clearly deteriorate, and $T=0.8,1.6$ do not converge to the target line—an example of how a violently varying Jacobian defeats the averaged approximation.

![Original Table 3.2: Jacobian-solve counts for serial trapezoidal integration and parallel BVM ParaDiag-I](assets/papers/time-parallelization/source-figures/table-3-2.svg)

Table 3.2 quantifies the parallel gain. With $N_t$ processors, the $N_t$ Jacobian systems in (3.46) are solved simultaneously, so the **parallel Jacobian-solve count equals the number of outer quasi-Newton rounds $k$**; serial time-stepping, by contrast, solves $\sum_n It_n$ systems ($It_n$ being the number of Newton steps at step $n$). In the table, at $\nu=0.1$ serial trapezoidal needs 401–443 solves and ParaDiag-I only 5–7 rounds; at $\nu=0.002$ ParaDiag-I grows from 7 to 22 rounds and fails on longer windows—consistent with the conclusion of Figure 3.13: when convergence is fast, the parallel advantage is extremely significant.

### Nearest Kronecker approximation

A single $I_t\otimes A_k$ approximates every time level with the **same** $A_k$, ignoring how the Jacobian's **magnitude** varies in time. The improvement proposed by Liu and Wu (2022b, Section 3.3) (originally for accelerating nonlinear ParaDiag-II, but equally applicable here) retains the tensor structure while allowing level-by-level scaling: approximate $\nabla F$ by $\Phi_k\otimes A_k$, where $\Phi_k=\operatorname{diag}(\phi_1,\ldots,\phi_{N_t})$ is determined by minimization,

$$
\min_{\Phi_k\ \mathrm{diagonal}}
\left\|\nabla F(\boldsymbol U^k)-\Phi_k\otimes A_k\right\|. \tag{3.47}
$$

Under the Frobenius norm, this is exactly the **nearest Kronecker approximation (NKA)**, with solution (Van Loan and Pitsianis 1993, Theorem 3, assuming $\operatorname{trace}(A_k^\top A_k)>0$)

$$
\phi_n=
\frac{\operatorname{trace}\!\left(
\nabla f(\boldsymbol u_n^k,t_n)A_k^\top
\right)}
{\operatorname{trace}(A_k^\top A_k)},
\qquad n=1,\ldots,N_t. \tag{3.48}
$$

That is, each Jacobian block is orthogonally projected onto the direction of $A_k$, and $\phi_n$ is the optimal scalar coefficient. The new quasi-Newton equation is

$$
(B\otimes I_x-\Phi_k\otimes A_k)\boldsymbol U^{k+1}
=\boldsymbol b-\left((\Phi_k\otimes A_k)\boldsymbol U^k-F(\boldsymbol U^k)\right).
$$

Since the time matrix is no longer $B$ but a combination of $B$ and $\Phi_k$, left-multiplying by $B^{-1}\otimes I_x$ turns it into
$(I_t\otimes I_x-B^{-1}\Phi_k\otimes A_k)\boldsymbol U^{k+1}=(B^{-1}\otimes I_x)(\boldsymbol b+F(\boldsymbol U^k))-(B^{-1}\Phi_k\otimes A_k)\boldsymbol U^k$,
so one only needs to diagonalize $B^{-1}\Phi_k=V\operatorname{diag}(\lambda_1,\ldots,\lambda_{N_t})V^{-1}$, and the second stage becomes $(I_x-\lambda_nA_k)\boldsymbol u_n^b=\boldsymbol u_n^a$. There is as yet no theory guaranteeing that $B^{-1}\Phi_k$ is diagonalizable in general, but in experiments this matrix is usually diagonalizable with a well-conditioned $V$.

Computing $\phi_n$ involves matrix–matrix multiplications and is fairly expensive, so it should not be recomputed at every Newton round. The practical recommendation is to compute the scale factors just once **offline** before the iteration, using a **coarse spatial model** (such as a coarse-grid semi-discrete ODE). The fine grid in Figure 3.14 is $\Delta x=1/200$, and $\{\phi_n\}$ is obtained offline from a coarse trapezoidal model with $\Delta X=1/20$.

![Original Figure 3.14: convergence comparison of averaged-Jacobian and NKA quasi-Newton iterations for Burgers' equation](assets/papers/time-parallelization/source-figures/figure-3-14.svg)

Figure 3.14 uses the periodic Burgers equation to compare the two quasi-Newton variants (averaged Jacobian vs. NKA), with the two panels corresponding to two viscosities, each containing $T=0.7$ and $T=1.3$. NKA is faster in all settings, with the clearest gain on the long window $T=1.3$—because the longer the window, the larger the variation of the Jacobian's magnitude in time, and the greater the information gain of the level-by-level scaling $\Phi_k$ over a single $A_k$.

> [!tip] Insight
> Going from $I_t\otimes A_k$ to $\Phi_k\otimes A_k$ is a "rank-one upgrade": the former assumes all time levels share the same Jacobian, while the latter admits that they share a direction but differ in magnitude. The cleverness of NKA is to unify the hard constraint of "preserving diagonalizability" (must be Kronecker-structured) with the soft goal of "approximating the true Jacobian as closely as possible" (Frobenius-optimal) in a single explicit closed form—$\phi_n$ is precisely the projection coefficient of $\nabla f_n$ onto $A_k$. This suggests a general recipe: when the exact operator destroys the parallel structure, rather than retreating to a single approximation, one can make an optimal projection within the **structure-preserving subspace**, trading a minimal extra offline cost for a substantial convergence speedup.

## Equation-and-figure audit

| Source item                                        | Paper section | Coverage                                                                                                                    |
| -------------------------------------------------- | ------------- | --------------------------------------------------------------------------------------------------------------------------- |
| (3.22)–(3.25)                                      | 3.5.1         | variable-step backward Euler, all-at-once matrix, Kronecker factorization, three-stage solve and parallel mechanism         |
| (3.26)–(3.29), Theorem 3.5                         | 3.5.1         | geometric grid, Jordan-block cause, two error classes, balance parameter, Toeplitz eigenvectors and proof route             |
| Figures 3.9–3.10                                   | 3.5.1         | all original figures, V-shaped parameter scan and parallel-width threshold                                                  |
| (3.30)–(3.37), Theorem 3.6, Figure 3.11, Table 3.1 | 3.5.1         | second-order reformulation, trapezoidal normalization, error balance, $r_1/r_2$ uniform bounds and condition-number plateau |
| (3.38)–(3.41), Theorem 3.7, Figure 3.12            | 3.5.1         | BVM, Dahlquist barrier, uniform second order, $O(N_t^2)$ conditioning and full elimination derivation                       |
| (3.42)–(3.46), Figure 3.13, Table 3.2              | 3.5.1         | nonlinear Newton, structure breakdown, averaged Jacobian, parallel solves and cost comparison                               |
| (3.47)–(3.48), Figure 3.14                         | 3.5.1         | NKA, Frobenius-optimal projection, offline coarse model and convergence improvement                                         |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 3.5–3.5.1, pp. 415–431.
