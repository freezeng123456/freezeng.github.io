---
title: "4.5: Diagonalization-Based Parareal"
description: Two complete routes—parallel coarse-grid correction and an interval-local diagonalized coarse propagator—with step-by-step derivations of the α-circulant all-at-once solve, the Theorem 4.7 threshold, and the Theorem 4.8 parabolic/hyperbolic spectral bounds
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-3-diagonalized-parareal
tags:
  - parallel-in-time
  - Parareal
  - ParaDiag
---

> [!note] Reading scope
> This page corresponds to Section 4.5 of the paper (pp. 460–472) and covers equations (4.14)–(4.29), Theorems 4.7–4.8, Remark 4.2, and Figures 4.12–4.17. Both methods use diagonalization, but at different locations and with different ranges of applicability: the first parallelizes the CGC across coarse points, and the second builds a special, parallelizable coarse propagator inside each coarse interval.

## 4.5 Diagonalization-based Parareal

This section presents a third class of Parareal variant, one that injects ParaDiag (the diagonalization technique in the time direction) into Parareal's coarse-grid correction (CGC). The motivation is that the bottleneck of standard Parareal lies precisely at the coarse level: the coarse correction (4.14) must advance sequentially through the $N_t$ coarse points, and this serial dependence limits parallel scaling. ParaDiag offers a way out—it uses an $\alpha$-circulant matrix to rewrite the sequential coupling as an all-at-once system that can be solved by a single diagonalization. Two essentially different routes appear in the literature:

- The first (Wu 2018; Wu and Zhou 2019) uses a head–tail coupling condition to recast the CGC itself into a form that ParaDiag can solve all at once;
- The second (Gander and Wu 2020) designs a special coarse propagator that closely mirrors the fine propagator and is solved cheaply by ParaDiag inside each large interval $[T_n,T_{n+1}]$.

The two routes differ in their mechanism, convergence properties, and range of applicability, and this contrast is the recurring thread of the section.

### Distinguishing the two routes first

- **Diagonalized CGC (Section 4.5.1):** modifies Parareal's sequential coarse correction across the $N_t$ coarse points. Its parallel width comes from the coarse time points; coarse and fine propagation may still use different integrators, the convergence mechanism remains close to standard Parareal, and it is mainly suited to parabolic problems.
- **Diagonalized coarse propagator (Section 4.5.2):** keeps the outer form of the standard Parareal coarse correction and, inside every $[T_n,T_{n+1}]$, uses ParaDiag to process $J$ fine steps simultaneously. Coarse and fine propagation use the same integrator and step size; this construction can transport long-lived frequencies and can therefore also handle hyperbolic problems.

> [!tip] Insight
> The dividing line between the two routes is whether the coarse propagator faithfully reproduces the fine propagator. Route one preserves Parareal's "cheap coarse, expensive fine" layering and merely parallelizes the sequential solve of the CGC, so it inherits the convergence factor and range of applicability of standard Parareal (see $\rho_{\mathrm{new}}=\rho$ in Theorem 4.7). Route two does the opposite: it makes the coarse propagator use exactly the same integrator and step size as the fine propagator, buying parallelism only through the head–tail coupling. It can therefore carry all frequency components as accurately as the fine propagator, breaking through the difficulty that standard Parareal faces on hyperbolic problems. All later derivations can be understood by returning to this contrast.

## 4.5.1 Diagonalization-based CGC

### From sequential CGC to head–tail coupling

The standard Parareal coarse-grid correction is

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

and it advances sequentially from $\boldsymbol u_0^{k+1}=\boldsymbol u_0$. Here $\boldsymbol b_{n+1}^k$ is the "fine minus coarse" residual, known from the previous iteration, and the initial value $\boldsymbol u_0^{k+1}=\boldsymbol u_0$ at $n=0$ stays fixed. It is exactly this chain of a fixed initial value plus point-by-point recursion that forces the CGC to execute serially. Wu's (2018) key idea is to break this chain: replace the fixed initial value with a head–tail coupling condition

$$
\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0,
$$

so that the "head" depends explicitly on the "tail." The recursion then becomes head-to-tail connected, forming a cyclic structure that can be diagonalized as a whole, but at the cost that the initial value at $n=0$ is no longer exactly $\boldsymbol u_0$. To make the converged limit still satisfy the original initial-value problem, the argument at node $0$ in the residual must be redefined to be the true initial value, i.e.

$$
\widetilde{\boldsymbol u}_n^k=
\begin{cases}
\boldsymbol u_0,&n=0,\\
\boldsymbol u_n^k,&n\ge1,
\end{cases}
$$

and one uses

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

The correction $\widetilde{\boldsymbol u}_0^k=\boldsymbol u_0$ is essential: without this substitution the converged fixed point would deviate from the solution of the original ODE; with it, one can verify that the converged limit exactly satisfies the discrete initial-value problem. It is worth noting that the head–tail condition used in (4.15) appeared a year before the "natural" head–tail condition (3.55) of ParaDiag-II, and its form differs slightly—(3.55) uses the difference-type condition $\boldsymbol u_0^{k+1}=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_0$, proposed a year later by Gander and Wu (2019)—but in the context of this section the two are equivalent in effect. This detail of "appearing early and later being endorsed by a more natural version" becomes important again in Remark 4.2.

### Linear all-at-once system and three-stage solve

Consider first the linear ODE $\boldsymbol u'=A\boldsymbol u$, $\boldsymbol u(0)=\boldsymbol u_0$, $t\in(0,T)$. Taking backward Euler at the coarse level gives $\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})=(I_x-\Delta TA)^{-1}\boldsymbol u_n^{k+1}$. Substituting this back into (4.15) (multiplying both sides by $I_x-\Delta TA$) yields the $N_t$ linear equations

$$
\begin{aligned}
(I_x-\Delta TA)\boldsymbol u_1^{k+1}&=\boldsymbol u_0^{k+1}+(I_x-\Delta TA)\boldsymbol b_1^k,\\
(I_x-\Delta TA)\boldsymbol u_2^{k+1}&=\boldsymbol u_1^{k+1}+(I_x-\Delta TA)\boldsymbol b_2^k,\\
&\ \,\vdots\\
(I_x-\Delta TA)\boldsymbol u_{N_t}^{k+1}&=\boldsymbol u_{N_t-1}^{k+1}+(I_x-\Delta TA)\boldsymbol b_{N_t}^k,
\end{aligned}
$$

together with $\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0$. These equations cannot be solved one by one, because the head–tail condition ties the right-hand side of the first equation to the solution of the last. Substituting the head–tail condition into the first coarse equation gives the all-at-once system

$$
(C_\alpha\otimes I_x-I_t\otimes\Delta TA)\boldsymbol U^{k+1}
=\boldsymbol g^k, \tag{4.16}
$$

where $\boldsymbol U^{k+1}=((\boldsymbol u_1^{k+1})^\top,\ldots,(\boldsymbol u_{N_t}^{k+1})^\top)^\top$,

$$
C_\alpha=
\begin{bmatrix}
1&&&-\alpha\\
-1&1\\
&\ddots&\ddots\\
&&-1&1
\end{bmatrix},
$$

$$
\boldsymbol g^k=
\begin{bmatrix}
\boldsymbol u_0+(I_x-\Delta TA)\boldsymbol b_1^k\\
(I_x-\Delta TA)\boldsymbol b_2^k\\
\vdots\\
(I_x-\Delta TA)\boldsymbol b_{N_t}^k
\end{bmatrix}.
$$

The $-\alpha$ in the top-right corner of $C_\alpha$ is precisely the cyclic entry left by the head–tail coupling; it turns the bidiagonal sequential structure into an $\alpha$-circulant matrix, which can be diagonalized by the discrete Fourier transform. Setting $C_\alpha=F\,\mathrm{diag}(\lambda_1,\ldots,\lambda_{N_t})\,F^*$ ($F$ the discrete Fourier matrix, $\lambda_n$ the eigenvalues, see (3.50)–(3.51)), solving (4.16) proceeds in three stages:

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

The three stages are: an FFT along time (stage a); solving in the frequency domain $N_t$ mutually independent, fully parallel shifted spatial systems $(\lambda_nI_x-\Delta TA)\,\cdot=\cdot$ (stage b); and an inverse FFT back to the time domain (stage c). A practical $\alpha$-circulant implementation also includes the diagonal scaling $\mathrm{diag}(\alpha^{(n-1)/N_t})$ needed to turn the $\alpha$-circulant into a standard circulant, exactly as in Section 3.5.2. Through (4.17), the new CGC (4.15) can be completed simultaneously at all $N_t$ coarse points, and the serial bottleneck is entirely removed.

> [!tip] Insight
> The middle stage of (4.17), "$N_t$ independent shifted spatial solves," is the computational kernel of the whole route, and the nonlinear quasi-Newton step (4.19) reuses it verbatim. In other words, $\alpha$-circulant diagonalization decouples a large "$N_t\times N_t$ time coupling plus space" system into "$N_t$ purely spatial systems" plus two FFTs—the cost in the time dimension drops from $O(N_t)$ sequential steps to the $O(\log N_t)$ depth of an FFT, which is exactly the fundamental advantage of ParaDiag over sequential time stepping.

### Theorem 4.7: threshold for keeping the standard Parareal speed

The parameter $\alpha$ faces a pair of competing constraints. On one hand, (4.15) shows that as $\alpha\to0$ the head–tail condition degenerates to $\boldsymbol u_0^{k+1}=\boldsymbol u_0$ and the whole method returns to standard CGC (4.14), so a "sufficiently small $\alpha$" should converge as fast as the original Parareal. On the other hand, the smaller $\alpha$ is, the larger the roundoff error introduced by diagonalizing the $\alpha$-circulant matrix $C_\alpha$ (see the analysis in Section 3.5.2), which is especially dangerous at low working precisions such as single or half precision. Fortunately, $\alpha$ need not be taken extremely small to match standard CGC: there is an explicit threshold.

**Theorem 4.7 (Wu 2018).** Let $\rho$ be the convergence factor of standard Parareal (4.14) and $\rho_{\mathrm{new}}$ the convergence factor of the new variant (4.15), with the coarse propagator $\mathcal G$ a stable time integrator. Then

$$
\rho_{\mathrm{new}}=\rho,
\qquad
\alpha\le\frac{\rho}{1+\rho}. \tag{Theorem 4.7}
$$

That is, as long as $\alpha$ does not exceed the threshold $\rho/(1+\rho)$, the asymptotic convergence factor of the new method is exactly the same as that of standard Parareal; reducing $\alpha$ further does not improve the asymptotic speed and only amplifies the roundoff risk. The optimal choice is therefore the threshold itself, $\alpha=\rho/(1+\rho)$. In practice $\rho=O(10^{-1})$, so $\alpha=\rho/(1+\rho)=O(10^{-1})$, and at this magnitude the roundoff error of $\alpha$-circulant diagonalization is negligible.

> [!tip] Insight
> One must distinguish the theorem's proven range from its numerically observed range. Theorem 4.7 is proved rigorously for the linear problem $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$ with $A$ having **negative real eigenvalues** (typical parabolic case). For other situations—for example $A$ with **complex eigenvalues**—only numerical experiments currently support the conclusion; there is no proof. It is therefore inaccurate to write "complex eigenvalues also apply" into the theorem's assumptions: this is a numerically observed, conjectural extension rather than a proven result. In engineering practice one may accordingly try the same threshold choice on complex-spectrum problems such as advection-dominated ones, but should be aware that the theoretical guarantee covers only the negative-real spectrum.

The intuitive meaning of the threshold is a "free lunch": the CGC is turned from strictly sequential into fully parallel, without sacrificing convergence speed. The only price is the negligible roundoff introduced by an $O(10^{-1})$ cyclic-coupling term. This is why route one is so attractive for parabolic problems.

![Original Figure 4.12: errors of standard and diagonalized CGC on the heat equation and the ADE](assets/papers/time-parallelization/source-figures/figure-4-12.svg)

The experiment uses periodic boundaries, $u_0(x)=\sin(2\pi x)$, backward Euler at the coarse level, SDIRK22 at the fine level, $T=4$, $J=10$, $\Delta T=0.1$, and $\Delta x=1/128$. For the diagonalized CGC three values of $\alpha$ are taken to observe how the convergence factor changes. Panel (a) is the heat equation, where standard CGC measures $\rho\approx0.22$, so Theorem 4.7 gives the threshold $\rho/(1+\rho)\approx0.18$; hence $\alpha=0.25,0.4$ exceed the threshold and are slower than standard CGC, while $\alpha=0.1$ is within the threshold and coincides with the standard curve. Panel (b) is the ADE at $\nu=0.1$, with $\rho\approx0.39$ and a threshold of about $0.28$; here $\alpha=0.1,0.25$ are both within the threshold and keep up with standard CGC, while $\alpha=0.4$ crosses the threshold and is clearly slower. The two panels each verify the same threshold formula under two spectral structures, and the theoretical prediction agrees closely with the measurements.

### Nonlinear all-at-once quasi-Newton

For the nonlinear problem $\boldsymbol u'=f(\boldsymbol u)$, $\boldsymbol u(0)=\boldsymbol u_0$, the coarse level again uses backward Euler. Still define

$$
\boldsymbol b_{n+1}^k
=\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol u}_n^k)
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k).
$$

Expanding the coarse correction $\boldsymbol u_{n+1}^{k+1}=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})+\boldsymbol b_{n+1}^k$ of (4.15) with backward Euler can be arranged into $(\boldsymbol u_{n+1}^{k+1}-\boldsymbol b_{n+1}^k-\boldsymbol u_n^{k+1})/\Delta T=f(\boldsymbol u_{n+1}^{k+1}-\boldsymbol b_{n+1}^k)$; combining this with the head–tail condition $\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0$ gives the nonlinear all-at-once system

$$
(C_\alpha\otimes I_x)\boldsymbol U^{k+1}
-\Delta T F(\boldsymbol U^{k+1})=\boldsymbol g^k, \tag{4.18}
$$

where $\boldsymbol U^{k+1}$ and $C_\alpha$ are as in (4.16), the $n$th block of $F$ is $f(\boldsymbol u_n^{k+1}-\boldsymbol b_n^k)$, and the first block of $\boldsymbol g^k$ is $\boldsymbol b_1^k+\boldsymbol u_0$ with the rest being $\boldsymbol b_2^k,\ldots,\boldsymbol b_{N_t}^k$. System (4.18) is solved by the same quasi-Newton iteration as nonlinear ParaDiag (Section 3.5.1):

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

where $A^{k+1,l}=\frac1J\sum_{j=1}^{J}\nabla f(\boldsymbol u_n^{k+1,l}-\boldsymbol b_n^k)$ is the average of the Jacobians at the time nodes, and $I_t\otimes A^{k+1,l}$ approximates the full block-diagonal Jacobian $\nabla F(\boldsymbol U^{k+1,l})$ with it.

> [!tip] Insight
> Using an "average Jacobian" rather than the "pointwise true Jacobian" is deliberate: after averaging, $P_\alpha^{k+1,l}$ retains exactly the $C_\alpha\otimes I_x-I_t\otimes(\cdot)$ tensor structure of (4.16), so the increment $\Delta\boldsymbol U^{k+1,l}$ can still be solved in parallel by the three-stage diagonalization of (4.17). With a pointwise Jacobian, the time blocks would no longer share the same spatial matrix, and the FFT diagonalization would fail. This is a typical tradeoff of "trading quasi-Newton for a diagonalizable structure." The paper notes that the nearest Kronecker product approximation (NKA) of Section 3.5.1 can give a more accurate approximation of $\nabla F$ than the average, but does not develop it for brevity. The nonlinear convergence analysis (Wu 2018, Section 4) shows that when $\alpha$ is suitably small, the convergence speed matches that of standard-CGC Parareal, with the same threshold mechanism.

![Original Figure 4.13: the two types of CGC on the viscous Burgers equation at two viscosities](assets/papers/time-parallelization/source-figures/figure-4-13.svg)

Figure 4.13 uses the Burgers equation, with the same problem setup and discretization parameters as the heat and ADE experiments above; the left and right panels correspond to $\nu=1$ and $\nu=0.01$, and each compares $\alpha=0.4,0.25,0.1$ with standard CGC. At both viscosities, $\alpha=0.4$ is slowest and $\alpha=0.1$ is closest to the standard curve (and even slightly faster in the weak-diffusion panel). The nonlinear case therefore shows the same effect of $\alpha$ on the convergence rate as the linear case, retaining the threshold structure of Figure 4.12.

### Remark 4.2: MGRiT needs a consistent head–tail condition

The underlying mechanism of MGRiT (4.12) is the same as Parareal, and its CGC can also be written in the form (4.14). But if the head–tail condition $\boldsymbol u_1^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_1$ of (4.15) is transplanted directly into MGRiT, it **diverges for every $\alpha$**. The reason is that this condition is not self-consistent at convergence—it is cruder than the "natural" condition given a year later by Gander and Wu (2019). The consistent condition for convergence should instead use the difference type

$$
\boldsymbol u_1^{k+1}
=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_1,
$$

which, at convergence ($\boldsymbol u_{N_t}^{k+1}=\boldsymbol u_{N_t}^k$), automatically degenerates to the exact relation $\boldsymbol u_1^{k+1}=\boldsymbol u_1$ and is therefore "consistent at convergence." Wu and Zhou (2019) used it to obtain the convergent MGRiT variant

$$
\left\{
\begin{aligned}
\boldsymbol u_0^{k+1}&=\boldsymbol u_0,\\
\boldsymbol u_1^{k+1}
&=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_1,\\
\boldsymbol u_{n+1}^{k+1}
&=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\widetilde{\boldsymbol b}_{n+1}^k,\quad n=1,\ldots,N_t-1.
\end{aligned}
\right. \tag{4.20}
$$

Here $\widetilde{\boldsymbol b}_{n+1}^k=\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol s}_n^k)-\mathcal G(T_n,T_{n+1},\widetilde{\boldsymbol s}_n^k)$, $\widetilde{\boldsymbol s}_n^k=\mathcal F(T_{n-1},T_n,\widetilde{\boldsymbol u}_{n-1}^k)$, and $\widetilde{\boldsymbol u}_n^k=\boldsymbol u_n$ ($n=0,1$), $\widetilde{\boldsymbol u}_n^k=\boldsymbol u_n^k$ ($n\ge2$). For suitably small $\alpha$, this variant is as fast as the original MGRiT (4.12), and the threshold mechanism of Theorem 4.7 applies in a similar way. Parareal itself can also replace (4.15) with this consistent difference-type head–tail condition $\boldsymbol u_0^{k+1}=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_0$, with the same convergence factor as in Theorem 4.7.

> [!tip] Insight
> This remark reveals an easily overlooked principle: when migrating a ParaDiag head–tail condition to a different iterative scheme, "consistency at convergence" is a prerequisite for convergence. Parareal tolerates the inconsistent (4.15) (because its error-propagation structure is different), yet MGRiT diverges completely because of it. The difference-type condition $\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)$ is more "natural" precisely because it vanishes automatically at the fixed point, introducing no bias into the converged limit.

## 4.5.2 Diagonalization-based coarse propagator

Gander and Wu (2020) proposed an idea essentially different from Section 4.5.1. Its key innovation is that the coarse and fine propagators use **the same time integrator and the same step size**, with the coarse propagator merely realized through diagonalization. Because the coarse propagator is no longer a cheap "large-step, strongly dissipative" approximation but is instead frequency-by-frequency consistent with the fine propagator, it can faithfully transport all frequency components over very long times—this is the root of its ability to handle hyperbolic problems, and the sharpest contrast with route one.

### Fine propagation and head–tail-coupled coarse propagation

Inside each large interval $[T_n,T_{n+1}]$ a linear-$\theta$ method is used with step size $\Delta t=\Delta T/J$ (the generalization to $s$-stage Runge–Kutta is given in the appendix of Gander and Wu 2020). The fine propagator $\mathcal F(T_n,T_{n+1},\boldsymbol u_n)=\boldsymbol v_J$ executes $J$ fine steps sequentially:

$$
\boldsymbol v_{j+1}-\boldsymbol v_j
=\Delta t[\theta f(\boldsymbol v_{j+1})
+(1-\theta)f(\boldsymbol v_j)],
\quad j=0,\ldots,J-1,
\quad \boldsymbol v_0=\boldsymbol u_n. \tag{4.21}
$$

$\theta=1$ is backward Euler, and $\theta=1/2$ is the trapezoidal rule. The special coarse propagator $\mathcal F_\alpha^*$ uses the **exact same** difference formula and only replaces the sequential initial value $\boldsymbol v_0=\boldsymbol u_n$ with the head–tail-coupled initial value

$$
\boldsymbol v_0=\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n. \tag{4.22}
$$

It is precisely this step that turns the $J$ fine steps inside the interval from a "sequential chain" into a "head-to-tail cyclic system," which can then be solved simultaneously by ParaDiag. Note the contrast with route one: route one performs head–tail coupling on the $N_t$ **coarse points**, while route two performs it on the $J$ **fine points** within a single interval.

### Nonlinear all-at-once system and quasi-Newton

Let $\boldsymbol V=(\boldsymbol v_1^\top,\ldots,\boldsymbol v_J^\top)^\top$; then (4.21)–(4.22) are written as the nonlinear all-at-once system

$$
\underbrace{(C_\alpha\otimes I_x)\boldsymbol V
-\Delta tF(\boldsymbol V)}_{K(\boldsymbol V)}
=\boldsymbol b(\boldsymbol u_n), \tag{4.23}
$$

$$
\boldsymbol b(\boldsymbol u_n)
=((1-\alpha)\boldsymbol u_n^\top,0,\ldots,0)^\top. \tag{4.24}
$$

Here $C_\alpha$ is as before, the first block of $F$ contains, because of the head–tail coupling, both $\theta f(\boldsymbol v_1)$ and $(1-\theta)f(\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n)$, and the remaining blocks are $\theta f(\boldsymbol v_j)+(1-\theta)f(\boldsymbol v_{j-1})$. The quasi-Newton update is

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
1-\theta&\theta\\
&\ddots&\ddots\\
&&1-\theta&\theta
\end{bmatrix},
$$

and $\overline{\nabla f}=\frac1J\big[\sum_{j=1}^{J-1}\nabla f(\boldsymbol v_j^l)+\nabla f(\alpha\boldsymbol v_J^l+(1-\alpha)\boldsymbol u_n)\big]$ is the average of the $J$ Jacobian blocks. $P_\alpha(\boldsymbol V^l)$ is a block $\alpha$-circulant matrix, used as an approximation of the true Jacobian $\nabla K(\boldsymbol V^l)=C_\alpha\otimes I_x-\Delta t(\widetilde C_{\theta,\alpha}\otimes I_x)\nabla F(\boldsymbol V^l)$.

> [!tip] Insight
> Compared with (4.19b) of route one, here there is an additional $\theta$-weighted circulant matrix $\widetilde C_{\alpha,\theta}$. The key is that $C_\alpha$ and $\widetilde C_{\alpha,\theta}$ belong to **the same family of circulant matrices** and can be **simultaneously diagonalized** by the same Fourier matrix, so the inner spatial systems of (4.25) can still be fully parallelized—the implicit/explicit weighting of the linear-$\theta$ method does not destroy the diagonalizable structure. This explains why one must choose a "circulant-friendly" integrator such as the linear-$\theta$ method (or its RK generalization).

The outer Parareal is still written in standard form

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal F_\alpha^*(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
-\mathcal F_\alpha^*(T_n,T_{n+1},\boldsymbol u_n^k). \tag{4.26}
$$

Note that (4.26) has exactly the same outer form as the standard Parareal coarse correction, only with the coarse propagator replaced by $\mathcal F_\alpha^*$—this is the embodiment of "keeping the CGC form and changing only the coarse propagator," complementing route one's "changing the CGC itself."

### Linear system, concurrency, and limiting cases

When $f(\boldsymbol u)=A\boldsymbol u$, (4.23) reduces to

$$
(C_\alpha\otimes I_x
-\widetilde C_{\theta,\alpha}\otimes\Delta tA)\boldsymbol V
=\boldsymbol b(\boldsymbol u_n), \tag{4.27}
$$

$$
\boldsymbol b(\boldsymbol u_n)
=([(I_x+\Delta t(1-\theta)A)(1-\alpha)\boldsymbol u_n]^\top,0,\ldots,0)^\top.
$$

The coarse propagation only takes the terminal point, $\mathcal F_\alpha^*=(H_J\otimes I_x)\boldsymbol V=\boldsymbol v_J$ ($H_J=(0,\ldots,0,1)\in\mathbb R^{1\times J}$). It is equivalent to solving

$$
\left\{
\begin{aligned}
\boldsymbol v_{j+1}-\boldsymbol v_j
&=\Delta tA[\theta\boldsymbol v_{j+1}+(1-\theta)\boldsymbol v_j],\\
\boldsymbol v_0&=\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n^k.
\end{aligned}
\right. \tag{4.28}
$$

Two extremes reveal the role of $\alpha$:

- At $\alpha=0$, the initial value of (4.28) degenerates to $\boldsymbol v_0=\boldsymbol u_n^k$, the coarse propagator equals the sequential fine propagator $\mathcal F$, and Parareal (4.26) **converges in a single iteration**—but because the $J$ steps must be solved sequentially, there is **no speedup at all**.
- At $0<\alpha<1$, (4.28) becomes a head–tail-coupled system that can be solved all at once by diagonalization, with the $J$ fine points in parallel; if the spatial systems have enough computing resources, the wall-clock cost of coarse propagation is about $1/J$ of sequential fine propagation.

> [!tip] Insight
> Here $\alpha$ is an "accuracy–parallelism" knob, with a different meaning than in route one. Route one's $\alpha$ mainly trades off "cyclic-coupling roundoff vs. keeping the standard convergence rate"; route two's $\alpha$ directly trades off "the fidelity of the coarse propagator to the fine propagator ($\alpha\to0$ is fully faithful but has no parallelism) vs. parallel speedup ($\alpha>0$ buys a $1/J$ cost but introduces a convergence factor)." Understanding this makes it clear why $\rho$ in Theorem 4.8 varies directly with $\alpha$.

### Theorem 4.8: parabolic and hyperbolic spectra

**Theorem 4.8 (Gander and Wu 2020).** For the linear initial-value problem $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$, $\boldsymbol u(0)=\boldsymbol u_0$, $A\in\mathbb C^{N_x\times N_x}$, let $\{\boldsymbol u_n^k\}$ be the $k$th iterate of variant (4.26) and $\{\boldsymbol u_n\}$ the converged solution. If both $\mathcal F$ and $\mathcal F_\alpha^*$ use a stable one-step Runge–Kutta method, and letting

$$
e^k=\max_{1\le n\le N_t}\|\boldsymbol u_n-\boldsymbol u_n^k\|_\infty,
$$

then

$$
e^k\le\rho^ke^0,
\qquad
\rho=
\begin{cases}
\alpha,&\sigma(A)\subset\mathbb R_-,\\[4pt]
\dfrac{2\alpha N_t}{1+\alpha},&\sigma(A)\subset i\mathbb R.
\end{cases} \tag{4.29}
$$

The two spectra give strikingly different behavior, which is the core conclusion of this section:

- **Negative real spectrum (parabolic):** $A$ comes from the semidiscrete heat equation $A\approx\Delta$ with $\sigma(A)\subset\mathbb R_-$, so here $\rho=\alpha$, **independent of the number of coarse intervals $N_t$**. Choosing $\alpha=O(10^{-2})$ already gives extremely fast and robust convergence.
- **Purely imaginary spectrum (hyperbolic):** the second-order wave equation (2.7), the Schrödinger equation, and similar problems have $\sigma(A)\subset i\mathbb R$ after semidiscretization, so here $\rho=2\alpha N_t/(1+\alpha)$, **growing linearly with $N_t$**. But this is an upper bound and can be very loose when $\alpha N_t$ is large, so it does not necessarily mean the actual convergence degrades—as long as $\alpha$ is small enough and $N_t$ is not too large, fast convergence is retained.

> [!tip] Insight
> Here appears the most crucial contrast of the whole section: why can route two handle hyperbolic problems while route one (Theorem 4.7, $\rho_{\mathrm{new}}=\rho$) stays confined to parabolic ones? The root lies in the construction of the coarse propagator. Route one's coarse propagator $\mathcal G$ is a **large-step, cheap integrator** with strong dissipation or phase error on high-frequency/oscillatory components; it inherits the old defect of standard Parareal on hyperbolic problems—"the coarse propagator cannot faithfully reproduce the fine propagator"—so its convergence factor can only equal that of standard Parareal (parabolic-friendly, hyperbolic-hard). Route two's coarse propagator $\mathcal F_\alpha^*$ uses **the same step size as the fine propagator** and transports every frequency accurately, with the only approximation coming from the head–tail-coupling perturbation introduced by $\alpha$; hence the hyperbolic convergence factor degrades only to the mild $2\alpha N_t/(1+\alpha)$—and the degree of degradation is controlled directly by $\alpha$ rather than being uncontrollable as with a large-step coarse propagator. In a word: route two downgrades "hyperbolic failure" from a structural defect to a tunable-parameter problem.

![Original Figure 4.14: the sharp rho=alpha prediction on the heat equation](assets/papers/time-parallelization/source-figures/figure-4-14.svg)

The heat equation uses homogeneous Dirichlet boundaries, $u_0=\sin^2(2\pi x)$, the trapezoidal rule, $\Delta T=1/12$, $J=10$, and $\Delta x=1/100$. The left and right panels take $N_t=36$ and $72$, and each compares $\alpha=10^{-1},10^{-2},10^{-3}$. The measured dashed curves are almost parallel to the theoretical dotted curves ($\rho=\alpha$), and doubling $N_t$ does not change the slope determined by $\rho=\alpha$—directly confirming that the negative-real-spectrum factor is independent of $N_t$ and that the bound is sharp.

![Original Figure 4.15: the joint influence of alpha and the number of coarse intervals on the wave equation](assets/papers/time-parallelization/source-figures/figure-4-15.svg)

The wave equation uses periodic boundaries, $u_0=\sin^2(2\pi x)$, $u_t(x,0)=0$; after semidiscretization the eigenvalues of $A\approx\Delta$ are all purely imaginary, so by (4.29) the convergence rate worsens as $N_t$ grows. Panel (a) fixes $\alpha=0.01$ (relatively large) and compares $N_t=24,48,96$: increasing the number of intervals slows convergence noticeably, verifying the linear trend of $2\alpha N_t/(1+\alpha)$. Panel (b) fixes $\alpha=10^{-4}$ (very small) and compares $N_t=24,48,96,960$: increasing $N_t$ from 24 to 960 needs only about two more iterations to reach the discretization-error magnitude $\max\{\Delta t^2,\Delta x^2\}$. The two panels separate the joint effect of $\alpha N_t$: the negative impact of $N_t$ can be effectively suppressed by a small $\alpha$.

![Original Figure 4.16: the theoretical factor is sharp for small alpha Nt and superlinear decay appears for a large product](assets/papers/time-parallelization/source-figures/figure-4-16.svg)

Unlike the heat equation, whether the factor in (4.29) is sharp for the wave equation depends on the product $\alpha N_t$. Figure 4.16 examines three pairs $(\alpha,N_t)$: only the small-product combination $\alpha=10^{-4},N_t=24$ hugs the dotted upper bound of (4.29); the other two measured curves exhibit **superlinear** decay, and the linear upper bound is clearly conservative. This shows that $2\alpha N_t/(1+\alpha)$ is a safe pessimistic estimate, and the actual performance is often better.

![Original Figure 4.17: iterations needed to reach 1e-8 on the Burgers equation](assets/papers/time-parallelization/source-figures/figure-4-17.svg)

The nonlinear convergence analysis (Gander and Wu 2020, Section 4), under the two assumptions of "exactly solving the all-at-once system (4.23)" and "$f$ satisfying a Lipschitz condition," gives $\rho=O(\alpha)$ for small $\alpha$, consistent with the linear case. The Burgers experiment uses periodic boundaries, $u_0=\sin^2(2\pi x)$, $\Delta T=0.1$, $J=10$, and $\Delta x=1/100$, with three curves for $\nu=1,0.01,10^{-4}$. Panel (a) fixes $N_t=40$ and shows how the number of iterations needed to reach a global error of $10^{-8}$ varies with $\alpha$: small $\alpha$ accelerates convergence; the dependence on $\nu$ shows that for small $\alpha$ the effect of viscosity is very weak, while for large $\alpha$ decreasing viscosity slows convergence. Panel (b) fixes $\alpha=10^{-3}$ and, as $N_t$ ranges from 10 to 160, the iteration count varies only between 2 and 5, showing that the convergence rate is robust to $N_t$.

## Final comparison of the two routes

The two variants of this section, (4.15) and (4.26), inject ParaDiag into standard Parareal in different ways: the former diagonalizes over the $N_t$ coarse time points and modifies the CGC; the latter diagonalizes over the $J$ fine time points within each large interval $[T_n,T_{n+1}]$, defining a special coarse propagator while keeping the CGC unchanged. Their ranges of applicability differ—the first, like standard Parareal, is mainly suited to parabolic problems, while the second is effective for both parabolic and hyperbolic problems.

| Problem                   | Diagonalized CGC (4.15)    | Diagonalized coarse propagator (4.26)                                |
| ------------------------- | -------------------------- | -------------------------------------------------------------------- |
| diagonalization direction | global $N_t$ coarse points | $J$ fine points per coarse interval                                  |
| modified component        | CGC                        | coarse propagator                                                    |
| coarse/fine integrators   | may differ                 | same integrator and step size                                        |
| main range                | parabolic problems         | parabolic and hyperbolic problems                                    |
| key parameter             | $\alpha\le\rho/(1+\rho)$   | parabolic factor $\alpha$; hyperbolic bound $2\alpha N_t/(1+\alpha)$ |

## Equation, theorem, and figure coverage audit

| Source item                            | Paper section | Coverage                                                                                                                                 |
| -------------------------------------- | ------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| (4.14)–(4.17)                          | 4.5.1         | standard/head–tail CGC, linear all-at-once matrix, three-stage parallel solve                                                            |
| Theorem 4.7, Figure 4.12               | 4.5.1         | $\alpha$ threshold, roundoff tradeoff, negative-real spectrum proven vs. complex-spectrum numerically observed, heat and ADE experiments |
| (4.18)–(4.19), Figure 4.13             | 4.5.1         | nonlinear system, average-Jacobian quasi-Newton, Burgers experiment                                                                      |
| Remark 4.2, (4.20)                     | 4.5.1         | MGRiT's consistent head–tail condition and convergent variant                                                                            |
| (4.21)–(4.26)                          | 4.5.2         | same-integrator fine/coarse propagation, nonlinear all-at-once system, quasi-Newton, outer update                                        |
| (4.27)–(4.28)                          | 4.5.2         | linearization, terminal extraction, $\alpha=0$ limit and $J$-way concurrency                                                             |
| Theorem 4.8, (4.29), Figures 4.14–4.17 | 4.5.2         | negative-real/purely-imaginary spectral bounds; all original heat, wave, and Burgers figures                                             |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4.5, pp. 460–472.
