---
title: "4.6: Space–Time Multigrid (STMG)"
description: A complete derivation and set of insights spanning the all-at-once system, temporal block-Jacobi smoothing, space–time transfer operators and local Fourier analysis, through optimal damping, integrator dependence, and nonlinear FAS
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-4-stmg
tags:
  - parallel-in-time
  - space-time multigrid
  - local Fourier analysis
---

> [!note] Reading scope
> This page corresponds to Section 4.6 of the paper (pp. 472–481) and covers equations (4.30)–(4.44), Theorem 4.9, Figures 4.18–4.22, and Table 4.1. The linear two-level cycle, smoothing symbol, damping choice, integrator dependence, and nonlinear FAS are all worked out in full, with step-by-step derivations at the "why" level.

## 4.6 Space–time multigrid (STMG)

Space–time multigrid (STMG) is the last parallel method introduced in this chapter, and its core idea is to apply multigrid (MG) simultaneously in both the spatial and temporal directions. After the early foundational work of Hackbusch (1984) and Horton & Vandewalle (1995), people gradually realized that a **block-Jacobi smoother in the time direction** is the key component that makes the whole method both parallel and efficient (Gander & Neumüller 2016). With it, STMG uses only standard multigrid building blocks yet achieves, on parabolic problems, an efficiency comparable to applying multigrid to the Poisson equation—that is, an iteration count that is essentially grid-independent. The following develops the method layer by layer, in the order of the paper, starting from the all-at-once system.

### All-at-once system (4.30)–(4.31)

Discretizing the heat equation or the advection–diffusion equation in space gives

$$
\boldsymbol u'=A\boldsymbol u+\boldsymbol f.
$$

Here $A\in\mathbb R^{N_x\times N_x}$ is the discrete matrix of the spatial derivatives: for the heat equation it is the Laplace operator $\partial_{xx}$, and for the advection–diffusion equation (ADE) it is $-\partial_x+\nu\partial_{xx}$, both discretized with mesh spacing $\Delta x$. Applying a generic one-step integrator to this system of ODEs gives

$$
r_1\boldsymbol u_{n+1}
=r_2\boldsymbol u_n+\widetilde{\boldsymbol f}_n,
\qquad n=0,\ldots,N_t-1, \tag{4.30}
$$

where $\boldsymbol u_0$ is the initial value and $r_1,r_2$ are matrix polynomials in $\Delta tA$. Backward Euler corresponds to $r_1=I_x-\Delta tA,\ r_2=I_x$; the trapezoidal rule corresponds to $r_1=I_x-\frac12\Delta tA,\ r_2=I_x+\frac12\Delta tA$.

As with [[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-5-paradiag-ii|ParaDiag]] (Section 3.5), we stack all $N_t$ difference equations into a single "all-at-once" system:

$$
\underbrace{
\begin{bmatrix}
r_1\\-r_2&r_1\\&\ddots&\ddots\\&&-r_2&r_1
\end{bmatrix}}_{K}
\underbrace{
\begin{bmatrix}
\boldsymbol u_1\\\boldsymbol u_2\\\vdots\\\boldsymbol u_{N_t}
\end{bmatrix}}_{\boldsymbol U}
=\boldsymbol b. \tag{4.31}
$$

Here $\boldsymbol b$ is the appropriate vector obtained by folding the initial value and the source term into the right-hand side. The matrix $K$ has a block lower bidiagonal structure: the main diagonal blocks are all $r_1$ and the subdiagonal blocks are all $-r_2$, which is exactly the algebraic manifestation of time stepping's "current step depends on the previous step." Rather than solving this system sequentially step by step, STMG builds a coarse-grid hierarchy in **space and time simultaneously** and solves for $\boldsymbol U$ as a whole with a single multigrid cycle.

> [!tip] Insight
> The shift from time stepping to the all-at-once system (4.31) is the watershed for understanding all parallel-in-time methods. Sequential time marching amounts to solving the block bidiagonal $K$ block by block via block forward substitution, whose data dependency is inherently serial; once we view $K\boldsymbol U=\boldsymbol b$ as a single monolithic linear system, we gain the freedom to redistribute computation across the entire space–time grid. The differences between STMG, ParaDiag, and Parareal are essentially different choices of "how to approximately solve the same $K$."

### Time-parallel block-Jacobi smoother (4.32)

STMG solves for $\boldsymbol U$ within the multigrid framework, and its smoother uses a **damped block-Jacobi iteration**. Starting from an initial approximation $\boldsymbol U^{\mathrm{ini}}$:

$$
\boldsymbol U^{\mathrm{new}}
=S_\eta(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s):
\left\{
\begin{aligned}
\boldsymbol U^0&=\boldsymbol U^{\mathrm{ini}},\\
(I_t\otimes r_1)\Delta\boldsymbol U^j
&=\eta(\boldsymbol b-K\boldsymbol U^j),
&&j=0,\ldots,s-1,\\
\boldsymbol U^{j+1}&=\boldsymbol U^j+\Delta\boldsymbol U^j,\\
\boldsymbol U^{\mathrm{new}}&=\boldsymbol U^s.
\end{aligned}
\right. \tag{4.32}
$$

Here $s$ is the number of smoothing sweeps and $\eta$ is the damping parameter. It is called "block Jacobi" because the iteration matrix takes the block-diagonal part $I_t\otimes r_1$ of $K$ (discarding the subdiagonal $-r_2$ coupling) and then damps by $\eta$. This step is the entire source of STMG's parallelism:

> [!tip] Insight: why "block Jacobi in the time direction" is the key to parallelism
> The subdiagonal block $-r_2$ of $K$ encodes the coupling "time step $n+1$ depends on time step $n$." The block-Jacobi smoother **throws these coupling terms wholesale into the residual** $\boldsymbol b-K\boldsymbol U^j$, while the left-hand side of the correction equation retains only the block diagonal $I_t\otimes r_1$. Since $I_t\otimes r_1$ is block diagonal, solving for $\Delta\boldsymbol U^j$ splits into $N_t$ mutually independent spatial systems $r_1\Delta\boldsymbol u_{n+1}^j=\eta(\cdots)_{n+1}$, which can be solved **fully in parallel** across all time points. By contrast, if one used pointwise Gauss–Seidel (see (4.35)) as in classical parabolic multigrid, the left-hand side would contain the temporal coupling and would have to be solved in time order, so parallelism would vanish. In other words, STMG demotes the "temporal coupling" from "an operator that must be solved serially" to "an explicit term that appears only in the residual," and it is precisely this demotion that buys time parallelism. The multigrid framework then compensates for the approximation: block Jacobi efficiently suppresses high-frequency error, and the remaining low-frequency error is left to coarse-grid correction.

For the smoother to work within multigrid, we also need restriction and prolongation operators in both space and time. Taking spatial linear interpolation with $N_x=7$ as an example:

$$
P_x=
\begin{bmatrix}
1/2&0&0\\
1&0&0\\
1/2&1/2&0\\
0&1&0\\
0&1/2&1/2\\
0&0&1\\
0&0&1/2
\end{bmatrix},
\qquad
R_x=\frac12P_x^\top. \tag{4.33}
$$

$P_x\in\mathbb R^{7\times3}$ is standard one-dimensional linear interpolation (prolongation), which interpolates the 3 coarse-grid interior points onto the 7 fine-grid interior points; the restriction is $R_x=\frac12P_x^\top$, i.e. full-weighting restriction. The time direction defines $P_t,R_t$ analogously. These are all **standard multigrid building blocks** with no special modification for the space–time problem—it is precisely the block-Jacobi smoother that makes standard building blocks sufficient.

### Two-level cycle (4.34)

Let `Mat` reshape an all-at-once vector into a "space $\times$ time" matrix, and let `Vec` do the inverse operation (in practice this is Matlab's `reshape`). This lets the spatial operator act by left multiplication and the temporal operator by right multiplication separately. One two-level cycle from iteration $k$ to $k+1$ is

$$
\left\{
\begin{aligned}
\boldsymbol U^{k+1/3}&=S_\eta(\boldsymbol b,\boldsymbol U^k,s_1),\\
\boldsymbol r&=\boldsymbol b-K\boldsymbol U^{k+1/3},\\
\boldsymbol r_c&=[R_x\operatorname{Mat}(\boldsymbol r)]R_t^\top,\\
\boldsymbol e_c&=K_c^{-1}\operatorname{Vec}(\boldsymbol r_c),\\
\boldsymbol e&=[P_x\operatorname{Mat}(\boldsymbol e_c)]P_t^\top,\\
\boldsymbol U^{k+2/3}&=\boldsymbol U^{k+1/3}+\operatorname{Vec}(\boldsymbol e),\\
\boldsymbol U^{k+1}&=S_\eta(\boldsymbol b,\boldsymbol U^{k+2/3},s_2).
\end{aligned}
\right. \tag{4.34}
$$

Its structure is the classical two-level multigrid V-cycle: first do $s_1$ pre-smoothing sweeps to suppress high-frequency error; then compute the residual and restrict it to the coarse grid simultaneously in space (left multiplication by $R_x$) and time (right multiplication by $R_t^\top$); solve the error equation $K_c\boldsymbol e_c=\boldsymbol r_c$ on the coarse level; prolongate the coarse error back to the fine level in space (left multiplication by $P_x$) and time (right multiplication by $P_t^\top$) and apply the correction; and finally do $s_2$ post-smoothing sweeps.

The coarse-level operator $K_c$ has exactly the same block structure as the fine-level $K$, only **rediscretized** on the doubled space–time steps $\Delta T=2\Delta t$, $\Delta X=2\Delta x$:

$$
K_c=
\begin{bmatrix}
r_1^c\\-r_2^c&r_1^c\\&\ddots&\ddots\\&&-r_2^c&r_1^c
\end{bmatrix},\qquad(N_t^c\ \text{blocks})
$$

where $r_1^c,r_2^c$ are matrix polynomials in $\Delta TA_c$, and $A_c\in\mathbb R^{N_x^c\times N_x^c}$ is the coarse discrete matrix of the spatial derivatives on $\Delta X$, for example

$$
\begin{cases}
r_1^c=I_x^c-\Delta TA_c,\quad r_2^c=I_x^c, & \text{backward Euler},\\[2pt]
r_1^c=I_x^c-\tfrac12\Delta TA_c,\quad r_2^c=I_x^c+\tfrac12\Delta TA_c, & \text{trapezoidal rule}.
\end{cases}
$$

In practice one takes $N_x=2^{l_x}-1,\ N_t=2^{l_t}-1$ (with $l_x,l_t\ge2$), so the coarse-level sizes are $N_x^c=2^{l_x-1}-1$ and $N_t^c=2^{l_t-1}-1$. Applying the two-level construction recursively to the coarse level naturally yields the full multilevel STMG.

> [!tip] Insight
> $K_c$ is obtained by "rediscretization" rather than by the "Galerkin triple product" $R K P$. This is because the time direction of STMG can be viewed as a strong advection term (see the historical remarks below), and the Galerkin coarse operator does not necessarily give a stable coarse-level integrator in this nonsymmetric, hyperbolic-like setting; rediscretizing directly on $2\Delta t,2\Delta x$ guarantees that the coarse level remains a reasonable time-stepping operator with the same block structure as the fine level.

### Difference from early parabolic multigrid (4.35)

STMG differs from the parabolic multigrid proposed by Hackbusch (1984) four decades ago in one key respect: the latter uses a **pointwise Gauss–Seidel** smoother,

$$
\boldsymbol U^{\mathrm{new}}
=S_{GS}(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s):
\left\{
\begin{aligned}
&\text{for }n=0,\ldots,N_t-1:\\
&\quad \boldsymbol u_{n+1}^0=\boldsymbol u_{n+1}^{\mathrm{ini}},\\
&\quad \text{for }j=0,\ldots,s-1:\\
&\qquad (D+L)\Delta\boldsymbol u_{n+1}^j=\widetilde{\boldsymbol f}_n+r_2\boldsymbol u_n^s-r_1\boldsymbol u_{n+1}^j,\\
&\qquad \boldsymbol u_{n+1}^{j+1}=\boldsymbol u_{n+1}^j+\Delta\boldsymbol u_{n+1}^j,\\
&\quad \boldsymbol u_{n+1}^{\mathrm{new}}=\boldsymbol u_{n+1}^s,
\end{aligned}
\right. \tag{4.35}
$$

where $\boldsymbol u_0^s=\boldsymbol u_0$, and $D,L$ are the diagonal and upper-triangular parts of $r_1$, respectively. Here $\boldsymbol U^{\mathrm{ini}}=(\boldsymbol u_0^\top,(\boldsymbol u_1^{\mathrm{ini}})^\top,\ldots)^\top$, and $\boldsymbol U^{\mathrm{new}}$ is correspondingly composed of $\boldsymbol u_0$ and $\boldsymbol u_1^{\mathrm{new}},\ldots,\boldsymbol u_{N_t}^{\mathrm{new}}$. This smoother is **strictly sequential in time**: one must finish smoothing at time step $n$ to obtain $\boldsymbol u_n^s$ before it can be used as the right-hand side to smooth time step $n+1$. After smoothing, the residual $\boldsymbol b-K\boldsymbol U^{\mathrm{new}}$ is likewise restricted to the coarse grid in space–time, the coarse problem is solved, and the process recurses.

Hackbusch (1984) coarsened only in **space** at the time and found that parabolic multigrid converged very fast for the heat equation. Gander & Lunet (2024) examined the two-level version that coarsens in space and time simultaneously and found that convergence **only slows down** in that case. This slow-convergence phenomenon had in fact long had a remedy: Horton & Vandewalle (1995) significantly improved convergence under space–time coarsening by **interpreting the time direction as a strong advection term** and designing multigrid building blocks accordingly (anisotropic/semi-coarsening-style smoothing and transfer operators); for related multigrid waveform relaxation variants, see Janssen & Vandewalle (1996) and Van Lent & Vandewalle (2002). STMG takes a different, complementary route—keeping the standard multigrid building blocks while concentrating the improvement in the temporal block-Jacobi smoother of (4.32): it makes all time points parallel and consistently suppresses high-frequency error before coarsening.

> [!tip] Insight
> The opposition between Gauss–Seidel and block Jacobi corresponds exactly to the trade-off between "sequential accuracy" and "parallel throughput." $S_{GS}$ sweeps forward in time and in a single pass propagates information from early times to late times, so a single smoothing pass is stronger and closer to the sequential solution, but it is completely non-parallel in time; $S_\eta$ gives up this in-time information propagation (the temporal coupling is thrown into the residual), so a single smoothing pass is weaker but buys $N_t$-way parallelism. Multigrid precisely makes up for the latter's weakness: the long-range temporal coupling that block Jacobi ignores is essentially low frequency, and it can be efficiently recovered by coarse-grid correction, so the combination of "weak smoothing + coarsening" preserves parallelism without sacrificing overall efficiency.

### Starting point of local Fourier analysis (4.36)–(4.39)

The basic principle for designing an effective smoother is to eliminate as much of the **high-frequency** error component as possible with as few smoothing sweeps as possible, so that the remaining low-frequency error can be accurately represented and eliminated on the coarse grid. The sharp tool for analyzing this is **local Fourier analysis (LFA)**—ignoring initial and boundary conditions and looking only at how the finite-difference stencil acts on a single Fourier mode of the error

$$
u_{n,m}^j=C_{\omega,\xi}^j
e^{i\omega n\Delta t}e^{i\xi m\Delta x}, \tag{4.36}
$$

where $\boldsymbol u_n^j:=(u_{n,1}^j,\ldots,u_{n,N_x}^j)^\top$, $i=\sqrt{-1}$, $\omega$ is the temporal frequency, and $\xi$ is the spatial frequency. For the one-dimensional heat equation with centered differences (space) + backward Euler (time):

$$
A=\frac1{\Delta x^2}\operatorname{Tri}(1,-2,1),
\qquad r_1=I_x-\Delta tA,\quad r_2=I_x.
$$

Setting the right-hand side $\boldsymbol b$ of (4.32) to zero and viewing $\boldsymbol U^j$ as the error at iteration $j$, each temporal-block equation of block Jacobi is

$$
r_1(\boldsymbol u_{n+1}^{j+1}-\boldsymbol u_{n+1}^{j})
=-\eta(r_1\boldsymbol u_{n+1}^j-r_2\boldsymbol u_n^j). \tag{4.37}
$$

First compute the result of the spatial discrete operator $A$ acting on the Fourier mode $u_{n+1,m}^l$ ($l=j,j+1$). Using the phase factors $e^{\pm i\xi\Delta x}$ from neighboring grid points and simplifying $e^{-i\xi\Delta x}-2+e^{i\xi\Delta x}=2(\cos(\xi\Delta x)-1)$:

$$
Au_{n+1,m}^l
=\frac{2(\cos(\xi\Delta x)-1)}{\Delta x^2}
C_{\omega,\xi}^l
e^{i\omega(n+1)\Delta t}e^{i\xi m\Delta x},
\qquad l=j,j+1. \tag{4.38}
$$

Thus the left-hand side is

$$
r_1(u_{n+1}^{j+1}-u_{n+1}^j)
=\Big(1-\tfrac{2\Delta t(\cos(\xi\Delta x)-1)}{\Delta x^2}\Big)
(C_{\omega,\xi}^{j+1}-C_{\omega,\xi}^j)\,e^{i\omega(n+1)\Delta t}e^{i\xi x_h},
$$

while the residual part on the right-hand side (note that $u_n^j$ differs from $u_{n+1}^j$ by a temporal phase $e^{-i\omega\Delta t}$) is

$$
r_1u_{n+1}^j-r_2u_n^j
=\Big(1-e^{-i\omega\Delta t}-\tfrac{2\Delta t(\cos(\xi\Delta x)-1)}{\Delta x^2}\Big)
C_{\omega,\xi}^j\,e^{i\omega(n+1)\Delta t}e^{i\xi x_h}.
$$

Substituting back into (4.37) and canceling the common phase factor on both sides yields the amplitude recurrence $C^{j+1}=\rho C^j$, where the **convergence factor** is

$$
\rho(\omega,\xi,\eta)
=1-\eta\left(
1-\frac{e^{-i\omega\Delta t}}
{1+\frac{2\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))}
\right), \tag{4.39}
$$

with $\omega\Delta t\in(-\pi,\pi)$ and $\xi\Delta x\in(-\pi,\pi)$. This symbol directly answers: given $\eta$, by how much is the error at any temporal/spatial frequency amplified or reduced after one smoothing sweep. The denominator $1+\frac{2\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))$ is exactly the dissipation factor brought by the implicitness of backward Euler—the higher the spatial frequency ($\cos(\xi\Delta x)\to-1$), the larger the denominator, the smaller the fraction, and the closer $\rho$ is to $1-\eta$.

### Theorem 4.9: optimal damping for backward Euler

Taking the maximum of (4.39) over $\xi,\omega$ and then minimizing that maximum (i.e. a min–max optimization) proves the following result (Gander & Lunet 2024, Ch. 4; a full analysis for more general discretizations is in Gander & Neumüller 2016):

> [!note] Theorem 4.9
> For the centered-difference–backward-Euler discretization of the one-dimensional heat equation, the optimal damping that makes the damped Jacobi smoother (4.32) **always permit time coarsening** is
> $$\eta_{\mathrm{opt}}=\tfrac12.$$
> In this case all temporal high frequencies $\omega\in\pm\big(\tfrac{\pi}{2\Delta t},\tfrac{\pi}{\Delta t}\big)$ are suppressed by at least a factor of $\tfrac1{\sqrt2}$. If the grid parameters further satisfy
> $$\frac{\Delta t}{\Delta x^2}\ge\frac1{\sqrt2},$$
> then the spatial high frequencies $\xi\in\pm\big(\tfrac{\pi}{2\Delta x},\tfrac{\pi}{\Delta x}\big)$ are also suppressed by at least a factor of $\tfrac1{\sqrt2}$, so one can coarsen in space **simultaneously**.

A more refined analysis of optimality is given by Chaudet-Dumas, Gander & Pogozelskyte (2024).

> [!tip] Insight: where does the coarsening condition $\Delta t/\Delta x^2\ge1/\sqrt2$ come from
> The prerequisite for multigrid to converge is that **every high frequency the coarse grid cannot represent must be suppressed by the smoother itself**. Temporal high frequencies can always be pushed below $1/\sqrt2$ by $\eta=1/2$, so time coarsening is "unconditionally" feasible. Spatial high frequencies are different—by (4.39), the spatial frequency acts only through the dissipation factor $\frac{2\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))$ in the denominator, whose strength is determined by the dimensionless number $\Delta t/\Delta x^2$. Only when $\Delta t/\Delta x^2\ge1/\sqrt2$, i.e. when implicit backward Euler dissipates spatial high frequencies sufficiently within one step, can block Jacobi push the spatial high frequencies below $1/\sqrt2$ as well, and only then does one dare to coarsen space simultaneously. If $\Delta t/\Delta x^2$ is too small (the time step is too fine relative to the space step), the spatial high frequencies are not smoothed adequately, and one should coarsen only in time (semi-coarsening). This condition is a direct, interpretable product of LFA, not an empirical fit.

For centered-difference ADE,

$$
A=\frac{\nu}{\Delta x^2}\operatorname{Tri}(1,-2,1)
+\frac1{2\Delta x}\operatorname{Tri}(-1,0,1),
$$

computing the result of $A$ acting on the Fourier mode in the same way, the diffusion term gives the real part $\frac{2\nu(\cos(\xi\Delta x)-1)}{\Delta x^2}$ and the centered difference of the advection term gives the purely imaginary part $\frac{i\sin(\xi\Delta x)}{\Delta x}$, so the symbol becomes

$$
\rho(\omega,\xi,\eta)
=1-\eta\left(
1-\frac{e^{-i\omega\Delta t}}
{1+\frac{2\nu\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))
+i\frac{\Delta t}{\Delta x}\sin(\xi\Delta x)}
\right). \tag{4.40}
$$

The extra imaginary part in the denominator is precisely the contribution of advection. Heuristically, one can still take $\eta=\tfrac12$ as the damping for time coarsening, as shown in the figure below.

![Source Figure 4.18: maximum high-frequency smoothing factors for the ADE at three viscosities](assets/papers/time-parallelization/source-figures/figure-4-18.svg)

Figure 4.18 plots the maximum of the convergence factor over the high-frequency range, $\rho_{\max}=\max_{(\xi\Delta x,\omega\Delta t)}\rho(\omega,\xi,\eta)$, as a function of $\eta$. From left to right it takes $\nu=0.1,0.01,0.001$, and each panel compares $\Delta x=\Delta t=1/64,1/128,1/256$. The minima on all three grids lie near $\eta\approx1/2$, indicating that this choice is fairly robust to the grid; as viscosity decreases, the minimum rises from about $0.71$ to about $0.79$, so the achievable high-frequency contraction weakens accordingly—the more advection dominates, the harder standard Jacobi smoothing works.

### Damping, sweep count, and integrator dependence

![Source Figure 4.19: error after five, ten, and fifteen cycles versus the damping parameter](assets/papers/time-parallelization/source-figures/figure-4-19.svg)

Figure 4.19 further validates $\eta=1/2$ with numerical experiments: two-level STMG does just one block-Jacobi sweep per cycle. Panel (a) is the heat equation and panel (b) is the ADE at $\nu=0.01$; each reports the error curve versus $\eta$ after 5, 10, and 15 cycles. The more iterations, the clearer the low-error valley near $\eta=1/2$. The heat-equation valley is broad, while the lowest point of the 15-cycle ADE curve lies somewhat toward $0.4$, so Figure 4.19 supports "$\eta=1/2$ as a robust empirical value" rather than claiming it is the exact optimum on every finite grid and fixed cycle count.

![Source Figure 4.20: error for one and three block-Jacobi smoothing sweeps](assets/papers/time-parallelization/source-figures/figure-4-20.svg)

Figure 4.20 fixes $\eta=1/2$ and examines the influence of the sweep count. Panels (a) and (b) use one and three block-Jacobi sweeps, respectively, and each compares the heat equation with the ADE at $\nu=0.1,0.01,0.001$. Two-level STMG converges faster with more smoothing sweeps for both equations; increasing to three sweeps raises the per-cycle cost but markedly reduces the number of cycles. The ADE converges worse overall than the heat equation, but an interesting phenomenon appears when the sweep count is large: several ADE curves in the right panel steepen in their later parts, entering a **superlinear** convergence stage, and their sensitivity to viscosity $\nu$ drops noticeably.

![Source Figure 4.21: STMG with the trapezoidal rule for different sweep counts and damping values](assets/papers/time-parallelization/source-figures/figure-4-21.svg)

Figure 4.21 reveals that the STMG convergence rate **depends on the time integrator**. Switching to the trapezoidal rule, the top row (a) for the heat equation uses 3, 5, and 10 sweeps in turn, and all damping scans stall at a large error or even diverge—even at 10 sweeps, the backward-Euler behavior cannot be recovered. The bottom row (b) is the ADE at $\nu=0.01$, using 2, 3, and 4 sweeps in turn; the error valley deepens as the sweep count increases, and the better damping falls near $\eta\approx0.8$, clearly different from the $\eta=1/2$ of backward Euler.

> [!tip] Insight: why the $\eta=1/2$ conclusion for backward Euler cannot be transplanted to the trapezoidal rule
> Theorem 4.9 holds precisely because backward Euler is **strongly dissipative** for high frequencies: in the denominator of (4.39), a temporal high frequency is significantly attenuated after just one backward-Euler step, so block Jacobi needs only a light touch. The trapezoidal rule is A-stable but **not L-stable**: its amplification factor tends to $-1$ rather than $0$ as $\Delta t\lambda\to-\infty$, so the highest-frequency error is barely dissipated by the time discretization itself and instead oscillates with sign flips at each step. These "persistent temporal high frequencies" cannot be represented on the coarse grid nor helped by the integrator's dissipation, so block-Jacobi smoothing fails—this is exactly the root cause of the top-row heat equation in Figure 4.21 "converging poorly no matter how $\eta$ is tuned." The ADE fares slightly better because of the extra phase mixing from advection, but its optimal damping is pushed up to about $0.8$. This shows that the STMG smoother and integrator must be designed to match, and it is the quantitative source of the method's "sensitivity to the integrator."

![Source Table 4.1: weak and strong scaling of STMG for the three-dimensional heat equation](assets/papers/time-parallelization/source-figures/table-4-1.svg)

Table 4.1 (from Gander & Neumüller 2016) shows full STMG solving the three-dimensional heat equation on a modern supercomputer with weak and strong scaling. Weak scaling grows from 1 core, 2 time steps, and 59,768 degrees of freedom all the way to 262,144 cores, 524,288 time steps, and 15,667,822,592 degrees of freedom; the **iteration count stays at 7 throughout**, and the wall-clock time increases only slightly from 28.8 seconds to about 30.0 seconds. By contrast, sequential time marching with space-only parallelism (column `fwd. sub.`) is estimated to explode from 19.0 seconds to 4,988,060 seconds. On the right, strong scaling at fixed problem size drives the wall-clock time down step by step from about 7,635.2 seconds to 30.0 seconds. This table shows that STMG's value comes from achieving **space–time parallelism** and a **grid-independent iteration count** at the same time—without either, constant wall-clock time cannot be maintained at this scale.

### Nonlinear systems and FAS (4.41)–(4.44)

STMG can be generalized to nonlinear problems

$$
\boldsymbol u'=f(\boldsymbol u),
\qquad \boldsymbol u(0)=\boldsymbol u_0,
\quad t\in(0,T), \tag{4.41}
$$

where $f:\mathbb R^{N_x}\to\mathbb R^{N_x}$ comes from the spatial discretization of some PDE. Applying the linear-$\theta$ method to (4.41) gives the nonlinear all-at-once system

$$
\underbrace{(B\otimes I_x)\boldsymbol U
-\Delta t(\widetilde B\otimes I_x)f(\boldsymbol U)}_{K(\boldsymbol U)}
=\boldsymbol b, \tag{4.42}
$$

where

$$
B=
\begin{bmatrix}
1\\-1&1\\&\ddots&\ddots\\&&-1&1
\end{bmatrix},
\qquad
\widetilde B=
\begin{bmatrix}
\theta\\1-\theta&\theta\\&\ddots&\ddots\\&&1-\theta&\theta
\end{bmatrix},
$$

$$
\boldsymbol b=
(\boldsymbol u_0^\top+\Delta t(1-\theta)f(\boldsymbol u_0)^\top,0,\ldots,0)^\top,
\qquad
\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top.
$$

$B$ encodes the bidiagonal coupling of the time difference, and $\widetilde B$ weights the nonlinear term $f$ at two adjacent times according to $\theta$. By analogy with (4.32), we first define the **nonlinear block-Jacobi smoother** $\boldsymbol U^{\mathrm{new}}=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s)$:

$$
\left\{
\begin{aligned}
\widetilde{\boldsymbol U}^0&=\boldsymbol U^{\mathrm{ini}},\\
\text{solve }\ \Delta\widetilde{\boldsymbol U}^j
-\Delta t\theta f(\Delta\widetilde{\boldsymbol U}^j)
&=\eta[\boldsymbol b-K(\widetilde{\boldsymbol U}^j)],
&&j=0,\ldots,s-1,\\
\widetilde{\boldsymbol U}^{j+1}
&=\widetilde{\boldsymbol U}^j+\Delta\widetilde{\boldsymbol U}^j,\\
\boldsymbol U^{\mathrm{new}}&=\widetilde{\boldsymbol U}^s.
\end{aligned}
\right. \tag{4.43}
$$

The correction $\Delta\widetilde{\boldsymbol U}^j$ is obtained by an inner solver (such as Newton iteration); since the left-hand side is still block diagonal (nonlinear coupling only within each temporal block), the $N_t$ local nonlinear corrections can still run in parallel across time. But nonlinearity makes LFA break down, so one cannot give a theoretically optimal value of $\eta$ as in Theorem 4.9 and can only choose it from experiments or other analyses.

Following Brandt (1977), we define nonlinear two-level STMG using the **full approximation scheme (FAS)**:

$$
\left\{
\begin{aligned}
\boldsymbol U^{k+1/3}&=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^k,s_1),\\
\boldsymbol r&=\boldsymbol b-K(\boldsymbol U^{k+1/3}),\\
\boldsymbol r_c&=[R_x\operatorname{Mat}(\boldsymbol r)]R_t^\top,\\
\boldsymbol U_c^{k+1/3}&=[R_x\operatorname{Mat}(\boldsymbol U^{k+1/3})]R_t^\top,\\
K_c(\boldsymbol U_c^{k+2/3})
&=\boldsymbol r_c+K_c(\boldsymbol U_c^{k+1/3}),\\
\boldsymbol e_c&=\boldsymbol U_c^{k+2/3}-\boldsymbol U_c^{k+1/3},\\
\boldsymbol e&=[P_x\operatorname{Mat}(\boldsymbol e_c)]P_t^\top,\\
\boldsymbol U^{k+2/3}&=\boldsymbol U^{k+1/3}+\operatorname{Vec}(\boldsymbol e),\\
\boldsymbol U^{k+1}&=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^{k+2/3},s_2).
\end{aligned}
\right. \tag{4.44}
$$

> [!tip] Insight: why the FAS correction $\boldsymbol r_c+K_c(\boldsymbol U_c^{k+1/3})$ is necessary
> In the linear case $K(\boldsymbol U+\boldsymbol e)-K(\boldsymbol U)=K\boldsymbol e$, so the coarse level can directly solve the error equation $K_c\boldsymbol e_c=\boldsymbol r_c$. A nonlinear operator has no such superposition property: the coarse level cannot represent the "error" on its own. FAS's approach is to also carry the restricted current solution $\boldsymbol U_c^{k+1/3}$ to the coarse level, solve for the **full solution** $K_c(\boldsymbol U_c^{k+2/3})=\boldsymbol r_c+K_c(\boldsymbol U_c^{k+1/3})$, and then take the **increment** of the coarse solution $\boldsymbol e_c=\boldsymbol U_c^{k+2/3}-\boldsymbol U_c^{k+1/3}$ as the correction. The added $K_c(\boldsymbol U_c^{k+1/3})$ on the right-hand side is exactly the "$\tau$-correction" that keeps the coarse and fine operators consistent at the current solution—without it, the coarse level would linearize at the wrong operating point and destroy convergence.

![Source Figure 4.22: Burgers STMG with two block-Jacobi smoothing sweeps](assets/papers/time-parallelization/source-figures/figure-4-22.svg)

Figure 4.22 is two-level nonlinear STMG on the Burgers equation (2.6), using two smoothing sweeps and the empirically optimal $\eta=1/4$. The $\nu=1$ curve crosses the discretization-error line in the figure after about 4 cycles and continues down toward $10^{-4}$; the $\nu=0.1$ curve is still above that line at cycle 16. As one can see, the viscosity dependence of nonlinear STMG matches the linear case (Figures 4.20–4.21): when viscosity is large enough, STMG also works well in the nonlinear setting, and as viscosity decreases and the equation becomes closer to hyperbolic, convergence deteriorates markedly.

Overall, for parabolic problems STMG is currently the **most efficient** parallel-in-time solver; but unlike Parareal, it is **highly intrusive** (it requires rewriting the linear/nonlinear solver, the transfer operators, and the all-at-once assembly). Moreover, as the top row of Figure 4.21 shows, even for parabolic problems the STMG convergence rate **depends on the time integrator used**, and this dependence still needs further study; and on hyperbolic problems (Figures 4.20, 4.22) STMG's efficiency drops noticeably, indicating that this area still needs more work.

## Equation, theorem, and figure coverage audit

| Source item                             | Paper section | Coverage status                                                                                |
| --------------------------------------- | ------------- | ---------------------------------------------------------------------------------------------- |
| (4.30)–(4.31)                           | 4.6           | generic one-step formula and all-at-once matrix $K$ (block lower bidiagonal)                   |
| (4.32)–(4.34)                           | 4.6           | parallel block Jacobi, space–time transfer $P_x,R_x$, full two-level cycle and $K_c$           |
| (4.35)                                  | 4.6           | early sequential Gauss–Seidel smoother and parallelism contrast, Horton–Vandewalle improvement |
| (4.36)–(4.40), Theorem 4.9, Figure 4.18 | 4.6           | full LFA derivation, heat/ADE symbols, optimal damping $\eta=1/2$, coarsening condition        |
| Figures 4.19–4.21, Table 4.1            | 4.6           | damping scan, sweep count, integrator dependence, weak/strong scaling                          |
| (4.41)–(4.44), Figure 4.22              | 4.6           | nonlinear all-at-once system, parallel nonlinear smoother, FAS consistency, Burgers test       |

## Further references

- Foundations and key components of space–time multigrid: W. Hackbusch (1984) parabolic multigrid; G. Horton & S. Vandewalle (1995) interpreting the time direction as a strong advection term to improve space–time coarsening convergence; the multigrid waveform relaxation variants of C. Janssen & S. Vandewalle (1996) and S. Van Lent & S. Vandewalle (2002); M. J. Gander & M. Neumüller (2016) establishing the temporal block-Jacobi smoother and giving large-scale scaling results.
- Convergence and optimality analysis: the two-level LFA and proof of Theorem 4.9 in M. J. Gander & T. Lunet (2024); for a refined analysis of optimality, see B. Chaudet-Dumas, M. J. Gander & A. Pogozelskyte (2024).
- Nonlinear framework: the full approximation scheme (FAS) of A. Brandt (1977).

## Source of this page

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4.6, pp. 472–481.
