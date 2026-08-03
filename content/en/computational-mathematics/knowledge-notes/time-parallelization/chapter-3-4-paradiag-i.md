---
title: "3.5.1: Direct ParaDiag (ParaDiag-I)"
description: A complete derivation from the all-at-once system to geometric time grids, boundary-value methods, and nonlinear quasi-Newton solvers
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i
tags:
  - parallel-in-time
  - ParaDiag
  - all-at-once systems
---

> [!note] Reading scope
> This page follows Sections 3.5–3.5.1 (pp. 415–430) of the paper. It covers the ParaDiag classification, equations (3.22)–(3.48), Theorems 3.5–3.7, Figures 3.9–3.14, and Tables 3.1–3.2. The truncation-error, roundoff-error, boundary-value-method, and nonlinear quasi-Newton branches are all retained.

## Section 3.5 introduction: the boundary between ParaDiag-I and ParaDiag-II

ParaDiag-I diagonalizes the time-stepping matrix directly. The paper notes that this is possible only when unequal time steps are used or when the formula is changed at the last time step. No outer iteration is required, but there are two limitations. Roundoff limits the number of time steps that can be processed in one window — with a geometric grid in double precision, the practical limit is usually about twenty steps. And the construction has **only been explored for a few low-order integrators**, namely backward Euler and the trapezoidal rule; it does not generalize easily to higher-order methods such as Runge–Kutta. A boundary-value-method discretization greatly improves the window length (Liu, Wang, Wu and Zhou 2022), but again only backward Euler and the trapezoidal rule apply.

ParaDiag-II replaces the time matrix by a diagonalizable approximation and uses that approximation in a stationary iteration or as a Krylov preconditioner. It gives up a direct solve in exchange for higher-order integrators, longer windows, and better-conditioned transforms. The next close-reading page treats ParaDiag-II.

## 3.5.1 Direct ParaDiag methods (ParaDiag-I)

### Backward Euler as an all-at-once system

Apply variable-step backward Euler to $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$:

$$
\frac{\boldsymbol u_n-\boldsymbol u_{n-1}}{\Delta t_n}
=A\boldsymbol u_n+\boldsymbol g_n,
\qquad n=1,\ldots,N_t. \tag{3.22}
$$

With $\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top$,

$$
K\boldsymbol U=\boldsymbol b,
\qquad K=B\otimes I_x-I_t\otimes A, \tag{3.23a}
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

If the $\Delta t_n$ are distinct, then

$$
B=VDV^{-1},
\qquad
D=\operatorname{diag}\!\left(
\frac1{\Delta t_1},\ldots,\frac1{\Delta t_{N_t}}
\right). \tag{3.24}
$$

The Kronecker structure yields

$$
K=(V\otimes I_x)
(D\otimes I_x-I_t\otimes A)
(V^{-1}\otimes I_x),
$$

and hence the three-stage solve

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

The $N_t$ shifted spatial systems in the middle stage are independent. The first and last stages are dense matrix multiplications in time.

![The time transform, independent spatial solves, and inverse transform in ParaDiag](assets/diagrams/pint/en/paradiag-three-stage.svg)

### Geometric grids and the two competing errors

Maday and Rønquist (2008) set $\Delta t_n=\mu^{n-1}\Delta t_1$ with $\mu>1$. Since the steps sum to $T$,

$$
\Delta t_n=\frac{\mu^{n-1}}{\sum_{j=1}^{N_t}\mu^{j-1}}T. \tag{3.26}
$$

Write $\mu=1+\varrho$. A large $\varrho$ enlarges the last steps and raises the truncation error. A small $\varrho$ makes $B$ close to a Jordan block and makes the eigenvector matrix ill-conditioned. Parameter selection in direct ParaDiag balances these effects.

### Theorem 3.5: balance formula for a first-order problem

Theorem 3.5 is taken from Gander, Halpern, Ryan and Tran (2016a), Theorems 2 and 6. Assume $\sigma(A)\subset\mathbb R_-$ and $|\lambda(A)|\leq\lambda_{\max}$. Let $\boldsymbol u_{N_t}(\varrho)$ and $\boldsymbol u_{N_t}(0)$ be the geometric-grid and uniform-grid backward Euler solutions at $T$, and let $\widetilde{\boldsymbol u}_n(\varrho)$ denote the value computed through diagonalization. Then

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

$\widetilde x_*$ maximizes $r(\widetilde x,N_t)$ on $[0,\infty)$ and $\lambda_*=N_t\widetilde x_*/T$. Balancing the two terms in (3.27) gives

$$
\varrho_{\mathrm{opt}}=
\left(
\epsilon
\frac{N_t^2(2N_t+1)(N_t+\lambda_{\max}T)}
{\phi(N_t)C(\lambda_*T,N_t)}
\right)^{1/(N_t+1)}. \tag{3.28}
$$

The proof reduces every spatial eigenmode to the Dahlquist equation $y'=\lambda y$. Its first estimate compares the two grids; its second bounds the roundoff introduced by the eigenvector transforms, with the worst case attained at $|\lambda|=\lambda_{\max}$. On the geometric grid,

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

where

$$
\mathbb T(a_1,\ldots,a_{N_t-1})=
\begin{bmatrix}
1\\
a_1&1\\
\vdots&\ddots&\ddots\\
a_{N_t-1}&\cdots&a_1&1
\end{bmatrix}. \tag{3.29b}
$$

These closed forms support the condition-number analysis. In computation, a numerical eigensolver is preferable because it rescales eigenvectors. Single- and double-precision machine epsilons are approximately $1.19\times10^{-7}$ and $2.22\times10^{-16}$.

### Figures 3.9–3.10: a finite parallel width remains

![Original Figure 3.9: error versus the geometric-grid parameter for the heat and advection–diffusion equations](assets/papers/time-parallelization/source-figures/figure-3-9.svg)

Figure 3.9 uses homogeneous Dirichlet data, $u_0(x)=\sin(2\pi x)$, $\Delta x=1/50$, and $T=0.2$; the advection–diffusion experiment has viscosity $10^{-2}$. Five values of $N_t$ are tested over $\varrho\in[10^{-2},1]$. The error is the maximum $L^\infty$ error over all time nodes. Each curve has a minimum, and the star shows (3.28). The prediction is close for advection–diffusion and less accurate for heat at small $N_t$.

![Original Figure 3.10: error first falls and then rises with the time-step count at the optimal geometric parameter](assets/papers/time-parallelization/source-figures/figure-3-10.svg)

Figure 3.10 sets $T=0.5$ and $N_t=2^4,2^5,\ldots,2^{10}$. Uniform-step backward Euler continues to improve. ParaDiag-I with the numerically optimal $\varrho_{\mathrm{num}}$ improves only up to fewer than one hundred steps, after which roundoff overwhelms discretization error.

### Wave equations and the trapezoidal rule

For

$$
\boldsymbol u''(t)=A\boldsymbol u(t),
\quad \boldsymbol u(0)=\boldsymbol u_0,
\quad \boldsymbol u'(0)=\widetilde{\boldsymbol u}_0, \tag{3.30}
$$

introduce

$$
\boldsymbol w'=\mathbb A\boldsymbol w,
\qquad
\boldsymbol w=(\boldsymbol u^\top,(\boldsymbol u')^\top)^\top,
\qquad
\mathbb A=\begin{bmatrix}0&I_x\\A&0\end{bmatrix}. \tag{3.31}
$$

The trapezoidal rule reduces wave-dispersion error:

$$
\frac{\boldsymbol w_n-\boldsymbol w_{n-1}}{\Delta t_n}
=\frac{\mathbb A}{2}(\boldsymbol w_n+\boldsymbol w_{n-1}). \tag{3.32}
$$

Its all-at-once system is

$$
K\boldsymbol W=\boldsymbol b,
\qquad K=B\otimes I_{2N_x}-\widetilde B\otimes\mathbb A, \tag{3.33a}
$$

$$
\widetilde B=\frac12
\begin{bmatrix}
1\\1&1\\&\ddots&\ddots\\&&1&1
\end{bmatrix}. \tag{3.33b}
$$

Multiplication by $\widetilde B^{-1}\otimes I_x$ gives

$$
\mathcal K\boldsymbol W=\widetilde{\boldsymbol b},
\quad
\mathcal K=\widetilde B^{-1}B\otimes I_{2N_x}-I_t\otimes\mathbb A,
\quad
\widetilde{\boldsymbol b}=(\widetilde B^{-1}\otimes I_x)\boldsymbol b. \tag{3.34}
$$

The time matrix has the factorization

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

### Theorem 3.6: balance formula for a wave problem

Theorem 3.6 is taken from Gander, Halpern, Rannou and Ryan (2019), Theorems 2.1 and 2.11. For $\lambda(A)\leq0$,

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

Balancing them yields

$$
\varrho_{\mathrm{opt}}=
\left(
\epsilon\frac{15\times2^{2N_t-1/2}}
{(N_t^2-1)(N_t-1)!}
\right)^{1/(N_t+1)}. \tag{3.37}
$$

The proof analyzes $u''+\lambda u=0$ mode by mode. It uses

$$
r_1(s)=\frac{s^3}{(1+s^2)^2}\leq\frac25,
\qquad r_2(s)=\frac1{1+s^2}\leq1,
$$

to obtain bounds uniform in the spatial eigenvalue.

![Original Figure 3.11: optimal geometric parameters and the step-count threshold for wave-equation ParaDiag-I](assets/papers/time-parallelization/source-figures/figure-3-11.svg)

Figure 3.11 uses homogeneous Dirichlet data, $\Delta x=1/20$, and $T=0.2$. Panel (a) confirms the minima predicted by (3.37). With a numerically optimal parameter, panel (b) shows rapid deterioration beyond $N_t=32$.

![Original Table 3.1: eigenvector condition numbers for backward Euler and the trapezoidal rule](assets/papers/time-parallelization/source-figures/table-3-1.svg)

For $N_t=5,10,20,30,60,100$, Table 3.1 shows that $\operatorname{Cond}(V)$ grows from $1.7\times10^3$ to $4.8\times10^6$ for backward Euler and from $4.7\times10^3$ to $4.1\times10^9$ for the trapezoidal rule. This explains the roundoff plateau and subsequent degradation in Figures 3.10–3.11.

### BVM: fixed steps and a different terminal formula

Liu et al. (2022) use a uniform $\Delta t$, the centered formula for the first $N_t-1$ equations, and backward Euler at the last step:

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

This is a boundary-value method (BVM), solved globally in time. Its stability is therefore not the stability of a step-by-step centered scheme. Axelsson and Verwer (1985) proved uniform second-order accuracy for the nonlinear case even though the last equation alone is first order.

$$
K\boldsymbol U=\boldsymbol b,
\qquad K=B\otimes I_x-I_t\otimes A, \tag{3.39a}
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
\quad
\boldsymbol b=
\begin{bmatrix}
\boldsymbol u_0/(2\Delta t)+\boldsymbol g_1\\
\boldsymbol g_2\\\vdots\\\boldsymbol g_{N_t}
\end{bmatrix}. \tag{3.39b}
$$

**Theorem 3.7.** The matrix $B$ admits $B=VDV^{-1}$ with $\operatorname{Cond}(V)=O(N_t^2)$. Closed forms for $V$, $V^{-1}$ and $D$ are given in Liu et al. (2022, Section 3).

> [!note] Site supplement: comparison with the geometric grid
> The paper does not place Theorem 3.7 next to the bounds of Theorems 3.5 and 3.6. Comparing them, the geometric grid's roundoff amplification carries $\varrho^{-(N_t-1)}$ and $2^{2N_t-1/2}/(N_t-1)!$, which grow very fast in $N_t$, whereas the BVM construction grows only polynomially as $O(N_t^2)$. Note, however, that the condition numbers measured in Table 3.1 under $\varrho_{\mathrm{num}}$ plateau as $N_t$ increases, a phenomenon the paper explicitly says merits further study, so the theoretical bounds should not be read as measured behaviour.

For a second-order problem, apply the same BVM to $\boldsymbol w$:

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

The velocity can be eliminated:

$$
(B\otimes I_x)\boldsymbol U-\boldsymbol V=\boldsymbol b_1,
\qquad
(B\otimes I_x)\boldsymbol V-A\boldsymbol U=\boldsymbol b_2,
$$

$$
(B^2\otimes I_x-I_t\otimes A)\boldsymbol U=\boldsymbol b, \tag{3.41}
$$

$$
\boldsymbol b=
\left(
\frac{\widetilde{\boldsymbol u}_0^\top}{2\Delta t},
-\frac{\boldsymbol u_0^\top}{4\Delta t^2},0,\ldots,0
\right)^\top.
$$

![Original Figure 3.12: errors and condition numbers for a geometric grid and the BVM on the wave equation](assets/papers/time-parallelization/source-figures/figure-3-12.svg)

Figure 3.12 uses $T=0.5$, $\Delta x=1/40$, and homogeneous Dirichlet data. Geometric-grid trapezoidal ParaDiag-I becomes roundoff-limited near $N_t=32$. The BVM retains $O(\Delta t^2)$ accuracy and tracks serial trapezoidal integration; its eigenvector condition number is far smaller.

### Nonlinear all-at-once equations and quasi-Newton iteration

For $\boldsymbol u'=\boldsymbol f(\boldsymbol u,t)$, define

$$
F(\boldsymbol U)=
\left(
\boldsymbol f(\boldsymbol u_1,t_1)^\top,\ldots,
\boldsymbol f(\boldsymbol u_{N_t},t_{N_t})^\top
\right)^\top.
$$

The nonlinear all-at-once equation is

$$
(B\otimes I_x)\boldsymbol U-F(\boldsymbol U)=\boldsymbol b. \tag{3.42}
$$

Newton's method can be written

$$
\left(B\otimes I_x-\nabla F(\boldsymbol U^k)\right)\boldsymbol U^{k+1}
=\boldsymbol b-\left(
\nabla F(\boldsymbol U^k)\boldsymbol U^k-F(\boldsymbol U^k)
\right), \tag{3.43a}
$$

$$
\nabla F(\boldsymbol U^k)=
\operatorname{blkdiag}\!\left(
\nabla f(\boldsymbol u_1^k,t_1),\ldots,
\nabla f(\boldsymbol u_{N_t}^k,t_{N_t})
\right). \tag{3.43b}
$$

The varying diagonal blocks destroy Kronecker separability. Following an idea of Gander and Halpern (2017), the paper approximates them by one average matrix:

$$
A_k=\frac1{N_t}\sum_{n=1}^{N_t}\nabla f(\boldsymbol u_n^k,t_n),
\qquad\text{or}\qquad
A_k=\nabla f\!\left(
\frac1{N_t}\sum_{n=1}^{N_t}\boldsymbol u_n^k,\frac{T}{N_t}
\right). \tag{3.44}
$$

Then

$$
(B\otimes I_x-I_t\otimes A_k)\boldsymbol U^{k+1}
=\boldsymbol b-\left((I_t\otimes A_k)\boldsymbol U^k-F(\boldsymbol U^k)\right). \tag{3.45}
$$

After $B=VDV^{-1}$, each iteration has the same three stages:

$$
\left\{
\begin{aligned}
\boldsymbol U^a&=(V^{-1}\otimes I_x)\boldsymbol r^k,\\
(\lambda_n I_x-A_k)\boldsymbol u_n^b&=\boldsymbol u_n^a,
&&n=1,\ldots,N_t,\\
\boldsymbol U^{k+1}&=(V\otimes I_x)\boldsymbol U^b.
\end{aligned}
\right. \tag{3.46}
$$

Here $\boldsymbol r^k$ is the right-hand side of (3.45). If the Jacobian changes strongly in time, one $A_k$ cannot approximate every block; the remedy is to shorten the window and process windows serially.

![Original Figure 3.13: BVM ParaDiag-I convergence for Burgers' equation at two viscosities and several window lengths](assets/papers/time-parallelization/source-figures/figure-3-13.svg)

Figure 3.13 uses periodic Burgers, $\Delta x=0.01$, and $N_t=T/\Delta t=200$. The horizontal line is $\max\{\Delta t^2,\Delta x^2\}=10^{-4}$. At $\nu=0.1$, convergence changes little between $T=0.1$ and $1.6$. At $\nu=0.002$, longer windows slow or defeat the iteration; $T=0.8$ and $1.6$ do not reach the target.

![Original Table 3.2: Jacobian-solve counts for serial trapezoidal integration and parallel BVM ParaDiag-I](assets/papers/time-parallelization/source-figures/table-3-2.svg)

With $N_t$ processors, all Jacobian systems in one ParaDiag-I iteration run concurrently, so the parallel count equals the outer iteration count; the corresponding serial trapezoidal count is $\sum_{n=1}^{N_t}\mathrm{It}_n$, the sum of the step-by-step Newton iterations. For $\nu=0.1$, Table 3.2 reports 401–443 serial solves and only 5–7 ParaDiag-I iterations. For $\nu=0.002$, the serial counts are 400/446/476/460/526 while the parallel count grows from 7 to 12 to 22, and the windows $T=0.8$ and $T=1.6$ fail. The paper explains the near-constant parallel count at $\nu=0.1$ by noting that $\nabla f$ varies little along the trajectory there, so a single $A_k$ represents all the blocks well.

### Nearest Kronecker approximation

This acceleration is due to Liu and Wu (2022, Section 3.3). A single $I_t\otimes A_k$ misses time-dependent amplitude changes. Use $\Phi_k\otimes A_k$, where $\Phi_k=\operatorname{diag}(\phi_1,\ldots,\phi_{N_t})$, and solve

$$
\min_{\Phi_k\ \mathrm{diagonal}}
\left\|\nabla F(\boldsymbol U^k)-\Phi_k\otimes A_k\right\|. \tag{3.47}
$$

Provided $\operatorname{trace}(A_k^\top A_k)>0$, the nearest Kronecker approximation under the Frobenius norm has the closed form (Van Loan and Pitsianis 1993, Theorem 3)

$$
\phi_n=
\frac{\operatorname{trace}\!\left(
\nabla f(\boldsymbol u_n^k,t_n)A_k^\top
\right)}
{\operatorname{trace}(A_k^\top A_k)},
\qquad n=1,\ldots,N_t. \tag{3.48}
$$

The quasi-Newton equation becomes

$$
(B\otimes I_x-\Phi_k\otimes A_k)\boldsymbol U^{k+1}
=\boldsymbol b-\left((\Phi_k\otimes A_k)\boldsymbol U^k-F(\boldsymbol U^k)\right).
$$

After multiplication by $B^{-1}\otimes I_x$, one diagonalizes $B^{-1}\Phi_k$ and solves $(I_x-\lambda_nA_k)\boldsymbol u_n^b=\boldsymbol u_n^a$. General diagonalizability has not been proved, but the matrices in the reported experiments are diagonalizable with well-conditioned eigenvectors.

Computing every $\phi_n$ on the fine grid is expensive. The paper therefore proposes a one-time offline computation on a coarse spatial model. Figure 3.14 uses $\Delta x=1/200$ for the fine problem and obtains the scale factors on $\Delta X=1/20$.

![Original Figure 3.14: convergence of average-Jacobian and nearest-Kronecker quasi-Newton iterations for Burgers' equation](assets/papers/time-parallelization/source-figures/figure-3-14.svg)

The two panels use different viscosities and compare $T=0.7$ with $T=1.3$. The nearest Kronecker approximation is faster in every case, with its clearest gain on the long window.

## Equation-and-figure audit

| Source item                                        | Paper section | Coverage                                                                             |
| -------------------------------------------------- | ------------- | ------------------------------------------------------------------------------------ |
| (3.22)–(3.25)                                      | 3.5.1         | variable-step backward Euler, all-at-once matrix, diagonalization, three-stage solve |
| (3.26)–(3.29), Theorem 3.5                         | 3.5.1         | geometric grid, error balance, optimal parameter, Toeplitz eigenvectors, proof route |
| Figures 3.9–3.10                                   | 3.5.1         | original figures, full settings, finite parallel width                               |
| (3.30)–(3.37), Theorem 3.6, Figure 3.11, Table 3.1 | 3.5.1         | first-order reformulation, trapezoidal all-at-once system, bounds, conditioning      |
| (3.38)–(3.41), Theorem 3.7, Figure 3.12            | 3.5.1         | BVM, second-order accuracy, polynomial conditioning, velocity elimination            |
| (3.42)–(3.46), Figure 3.13, Table 3.2              | 3.5.1         | nonlinear Newton system, average Jacobian, concurrent solves, cost comparison        |
| (3.47)–(3.48), Figure 3.14                         | 3.5.1         | nearest Kronecker approximation, coarse-model scaling, convergence gain              |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 3.5–3.5.1, pp. 415–430.
