---
title: "4.1–4.2: Historical Context and Parareal"
description: A close reading from multiple shooting and coarse/fine propagators to the linear and nonlinear convergence bounds
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-1-parareal
tags:
  - parallel-in-time
  - Parareal
  - parabolic equations
---

> [!note] Reading scope
> This page follows Sections 4, 4.1, and 4.2 (pp. 443–452). It covers equations (4.1)–(4.9), Theorems 4.1–4.4, Remark 4.1, and Figures 4.1–4.5. The reasoning after each formula follows the proof order in the paper.

## 4.1 Historical development

### Why the paper separates methods designed for parabolic problems

Chapter 2 showed that parabolic solutions are local in time except at very low frequencies, whereas every hyperbolic frequency may propagate over a long horizon. The Chapter 3 methods address long-range temporal coupling and often work on parabolic problems as well. Their nonlinear forms retain specific difficulties: an effective OSWR Robin parameter can be hard to find, and the outer Newton iterations in ParaExp and ParaDiag may slow down or fail on long windows.

Chapter 4 exploits dissipation-induced temporal locality. Parareal, PFASST, MGRiT, and STMG handle linear and nonlinear problems effectively when diffusion is sufficient. Their performance degrades continuously toward weakly diffusive and hyperbolic regimes.

Parareal has roots in the noniterative precursor of Nievergelt (1964) and in multiple shooting (Bellen and Zennaro 1989; Chartier and Philippe 1993). Saha, Stadel and Tremaine (1997) had already presented the algorithm with a coarse model and noted its relation to waveform relaxation. Lions, Maday and Turinici (2001) introduced the modern method independently in the context of virtual control, emphasizing that it is non-intrusive. Convergence theory appears in Gander and Vandewalle (2007), Gander and Hairer (2008, 2014), and Gander and Lunet (2024). Later descendants include PITA (Farhat and Chandesris 2003; Farhat et al. 2006; Cortial and Farhat 2009), PFASST (Minion 2011; Emmett and Minion 2012; Minion et al. 2015), MGRiT (Falgout et al. 2014; Dobrev et al. 2017; Hessenthaler et al. 2020), and Parareal–ParaDiag combinations (Wu 2018; Gander and Wu 2020). Space–time multigrid followed another line, beginning with the parabolic multigrid of Hackbusch (1984) and the multigrid waveform relaxation of Lubich and Ostermann (1987); early schemes could not coarsen time effectively, while the temporal block-Jacobi smoother of Gander and Neumüller (2016) enabled scalable STMG.

The paper adds two framing statements: Parareal can be regarded as a template for developing more efficient PinT methods, and Parareal-based methods use two grids (or more) in time while using just one grid in space.

## 4.2 Parareal

### The update and the two temporal grids

Partition $[0,T]$ by $0=T_0<T_1<\cdots<T_{N_t}=T$. Let $\mathcal F$ be an accurate, expensive fine propagator and $\mathcal G$ a cheap coarse propagator. Starting from interface guesses $\boldsymbol u_n^0$, Parareal computes

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k). \tag{4.1}
$$

All fine solves at iteration $k$ are concurrent. The coarse propagation carrying the new index $k+1$ advances sequentially in $n$. The last two terms form a coarse prediction–correction and can also be viewed as a finite-difference Jacobian in a multiple-shooting Newton method.

![Original Figure 4.1: each coarse time step contains J fine time steps](assets/papers/time-parallelization/source-figures/figure-4-1.svg)

The analysis uses uniform grids with $\Delta T/\Delta t=J\ge2$, though nonuniform grids are possible. The target is the sequential discrete solution produced by $\mathcal F$; fixing $\mathcal F$ fixes the result to which Parareal converges.

### Theorem 4.1: a modal Toeplitz error iteration

For $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$, let $A=V_ADV_A^{-1}$ and let $R_g,R_f$ be the coarse and fine one-step stability functions. Fine propagation over one coarse interval is $R_f^J(z/J)$, with $z=\Delta T\lambda(A)$.

If $|R_g(z)|\le1$, then

$$
\max_{1\le n\le N_t}
\|V_A(\boldsymbol u_n^k-\boldsymbol u_n)\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}\|M^k(z)\|_\infty
\max_{1\le n\le N_t}
\|V_A(\boldsymbol u_n^0-\boldsymbol u_n)\|_\infty, \tag{4.2}
$$

where $\boldsymbol u_n$ is the sequential fine solution and

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

#### Proof chain

Equation (4.1) becomes

$$
\boldsymbol u_{n+1}^{k+1}
=R_f^J(\Delta TA/J)\boldsymbol u_n^k
+R_g(\Delta TA)\boldsymbol u_n^{k+1}
-R_g(\Delta TA)\boldsymbol u_n^k.
$$

The sequential fine solution satisfies the same identity after adding and subtracting the coarse term. Thus $\boldsymbol e_n^k=\boldsymbol u_n-\boldsymbol u_n^k$ obeys

$$
\boldsymbol e_{n+1}^{k+1}
=R_g(\Delta TA)\boldsymbol e_n^{k+1}
+[R_f^J(\Delta TA/J)-R_g(\Delta TA)]\boldsymbol e_n^k.
$$

Diagonalization separates the modes:

$$
\xi_{n+1}^{k+1}(z)
=R_g(z)\xi_n^{k+1}(z)
+[R_f^J(z/J)-R_g(z)]\xi_n^k(z),
$$

with

$$
\|V_A\boldsymbol e_n^k\|_\infty
=\max_{z\in\sigma(\Delta TA)}|\xi_n^k(z)|. \tag{4.4}
$$

Stacking the temporal components gives

$$
M_g(z)\boldsymbol\xi^{k+1}
=[M_g(z)-M_f(z)]\boldsymbol\xi^k,
$$

so $\boldsymbol\xi^k=M^k(z)\boldsymbol\xi^0$, which proves (4.2) after taking both maxima.

> [!note] Remark 4.1: preconditioner interpretation
> $M(z)=I_t-M_g^{-1}(z)M_f(z)$. If $M_fU=b$ is the fine
> all-at-once system, one preconditioned correction is
>
> $$
> M_g(z)\Delta U^k=r^k:=b-M_f(z)U^k,
> \qquad U^{k+1}=U^k+\Delta U^k.
> $$
>
> Its $n$th residual block is
>
> $$
> r_n^k=b_n-
> \left[u_n^k-\mathcal F(T_{n-1},T_n,u_{n-1}^k)\right]
> =b_n-u_n^k+R_f^J(z/J)u_{n-1}^k.
> $$
>
> Fine propagation in the residual blocks is concurrent; application
> of the coarse preconditioner remains sequential. Section 4.5 replaces
> this serial correction by a diagonalizable one.

The strictly lower triangular structure also proves exactness at the first $k$ coarse points after iteration $k$ and termination in at most $N_t$ iterations in exact arithmetic.

### Theorem 4.2: short-time superlinear and long-time linear convergence

The error symbol here differs slightly from Theorem 4.1: $\boldsymbol e_n^k=V_A(\boldsymbol u_n^k-\boldsymbol u_n)$ is the modal error. The lower triangular entries of $M_g^{-1}$ are powers of $R_g(z)$, so

$$
M(z)=[R_f^J(z/J)-R_g(z)]\widetilde M(R_g(z)),
$$

where the first subdiagonal of $\widetilde M(\beta)$ is $1$ and deeper diagonals contain $\beta,\beta^2,\ldots$. Gander and Vandewalle (2007, Lemma 4.4) bound its powers by

$$
\|\widetilde M^k(R_g)\|_\infty\le
\begin{cases}
\min\left\{
\left(\dfrac{1-|R_g|^{N_t-1}}{1-|R_g|}\right)^{\!k},\
\dbinom{N_t-1}{k}
\right\},
&|R_g|<1,\\[10pt]
\dbinom{N_t-1}{k},
&|R_g|=1,
\end{cases}
$$

and the two regimes below come from the two branches: the binomial coefficient produces the $\prod_{j=1}^k(N_t-j)/k!$ of (4.5a), and the geometric sum produces the $N_t$-independent (4.5b). For short horizons,

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

The product vanishes at $k=N_t$. If $|R_g(z)|<1$, a horizon-independent estimate is

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

The numerator measures coarse–fine mismatch; the denominator measures the dissipation margin of the coarse method. A nearly nondissipative coarse mode combined with a large mismatch makes the factor approach or exceed one.

![Original Figure 4.2: short-time superlinear and long-time linear convergence](assets/papers/time-parallelization/source-figures/figure-4-2.svg)

Figure 4.2 uses the periodic heat equation, zero source,
$u_0(x)=\sin^2(2\pi x)$, $\Delta x=1/5$, backward Euler on both
levels, and $J=10$. The left panel has $(T,N_t)=(0.02,6)$ and is
predicted by $\varrho_s$; the right has $(T,N_t)=(0.5,64)$ and shows
the fixed slope described by $\varrho_l$. A mesh refined to
$\Delta x=1/8$ also converges linearly, but the source does not claim
that it enters the same two-stage process earlier.

### Theorem 4.3: nonlinear superlinear estimate

Theorem 4.3 is due to Gander and Hairer (2008, Theorem 1; see also
Gander and Lunet 2024, Theorem 2.6). Let $\mathcal F$ be exact and
$\mathcal G$ an order-$p$ method with local error at most
$C_3\Delta T^{p+1}$. On the relevant bounded set of states, assume

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

with continuously differentiable coefficients. Their derivative bounds
give a defect Lipschitz constant $C_1$:

$$
\|(\mathcal F-\mathcal G)(\boldsymbol v)
-(\mathcal F-\mathcal G)(\boldsymbol w)\|
\le C_1\Delta T^{p+1}\|\boldsymbol v-\boldsymbol w\|.
$$

For sufficiently small $\Delta T$,

$$
\|\boldsymbol u(T_n)-\boldsymbol u_n^k\|
\le
\frac{C_3}{C_1}
\frac{(C_1\Delta T^{p+1})^{k+1}}{(k+1)!}
(1+C_2\Delta T)^{n-k-1}
\prod_{j=0}^{k}(n-j). \tag{4.6}
$$

The product vanishes once $k\ge n$, and every iteration introduces another factor of $\Delta T^{p+1}$ on short intervals.

> [!warning] Source check: the extra $\Delta T^{p+1}$ in (4.6)
> The journal and arXiv versions print
> $C_3\Delta T^{p+1}(C_1\Delta T^{p+1})^{k+1}$.
> The cited Gander–Hairer theorem instead has
> $(C_3/C_1)(C_1\Delta T^{p+1})^{k+1}$, as displayed above. The
> one-step case $n=1,k=0$ makes the issue immediate: a coarse local
> error is $O(\Delta T^{p+1})$, not $O(\Delta T^{2p+2})$.

### Theorem 4.4: the parabolic long-time factor near 0.3

If $\mathcal G$ is backward Euler and $\mathcal F$ an L-stable Runge–Kutta method, then for some $J_{\min}=O(1)$,

$$
\max_{z\in\mathbb R_-}\varrho_l(J,z)\approx0.3,
\qquad J\ge J_{\min}. \tag{4.7}
$$

The theorem also assumes $A$ is negative semi-definite. This factor is
independent of $T$ and $N_t$. Backward Euler for $\mathcal F$ is
treated by Mathew, Sarkis and Schaerer (2010), and general L-stable
Runge–Kutta methods by Yang, Yuan and Zhou (2023). Wu (2015) and Wu and
Zhou (2015) also analyze BDF2 and the SDIRK methods below. The
trapezoidal rule is not an L-stable case of this theorem; it belongs to
the bounded-spectrum extension (4.8).

The result originates in Gander and Vandewalle (2007, Table 5.1) **at the continuous level** (that is, with $\mathcal F$ the exact propagator $\exp(\Delta TA)$) and **for other coarse propagators**, where the contraction can be even better — about $0.068$ for Radau IIA. Note that what changes there is the coarse propagator, not the fine one.

For an A-stable but non-L-stable fine method such as the trapezoidal
rule, note that $z=\Delta T\lambda(A)\le0$. The bounded-spectrum
statement should therefore be written

$$
\max_{z\in[-z_{\max},0]}\varrho_l(J,z)\approx0.3,
\qquad
J\ge J_{\min}=O(\log^2z_{\max}). \tag{4.8}
$$

> [!warning] Source check: the spectral interval in (4.8)
> The published formula uses $z\in[0,z_{\max}]$, although this section
> defines $z=\Delta T\lambda(A)$ with negative semidefinite $A$. An
> equivalent positive variable is $s=-z$, giving
> $\max_{s\in[0,z_{\max}]}\varrho_l(J,-s)$.

More fine steps are needed to resolve the dissipative high-frequency physics. This differs sharply from the case where $\mathcal F$ is the exact propagator $\exp(\Delta TA)$, for which a rate around $0.3$ holds already for $J\ge2$. Equation (4.8) was proved for the trapezoidal rule and a fourth-order Gauss Runge–Kutta method in Wu and Zhou (2015). Equation (4.9) gives the two schemes explicitly:

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

The bound holds from $J_{\min}=2$ for SDIRK22 and from $J_{\min}=4$ for SDIRK23 (Wu 2015; Wu and Zhou 2015).

> [!warning] Source check: the SDIRK22 Butcher tableau
> Both the journal version and the arXiv preprint print the second SDIRK22 node as $c_2=1-\gamma$ and the weights as $b=(1-\gamma,1-\gamma)$ in (4.9). With $\gamma=(2-\sqrt2)/2$ those weights sum to $\sqrt2\ne1$, and $c_2$ contradicts the row sum $a_{21}+a_{22}=1$. The tableau above uses the standard L-stable, stiffly accurate SDIRK22, with $c_2=1$ and $b$ equal to the last row of $A$. The adjacent SDIRK23 matches the paper and needs no correction.

![Original Figure 4.3: Parareal convergence under different fine propagators and coarsening factors](assets/papers/time-parallelization/source-figures/figure-4-3.svg)

Figure 4.3 uses the periodic heat equation with $\Delta x=1/256$, $\Delta T=0.1$, $T=4$, and diffusion coefficient $0.1$. The panels from left to right use $J=2,10,50$. At $J=2$, the trapezoidal curve stalls near $10^{-4}$ and SDIRK23 is visibly slower than SDIRK22. At $J=10$, the two SDIRK curves nearly coincide while the trapezoidal rule retains a slow tail. At $J=50$, all three follow the $0.3^k$ reference slope. The transition verifies the qualification in Theorem 4.4: the fine integrator's stability class and the coarse-to-fine ratio must be assessed together.

### Deterioration as diffusion weakens

The remaining tests fix $T=4$, $\Delta T=0.1$, $\Delta x=1/128$, and $J=32$, with backward Euler coarse propagation and SDIRK22 fine propagation.

![Original Figure 4.4: long-time factors for every advection–diffusion eigenvalue at three viscosities](assets/papers/time-parallelization/source-figures/figure-4-4.svg)

The experiment uses a zero source term and $u(x,0)=\sin(2\pi x)$. As viscosity falls, $\max\varrho_l$ approaches one: the three panels report $\varrho_{l,\max}=0.23$ at $\nu=1$, $0.39$ at $\nu=0.1$, and $0.79$ at $\nu=0.02$, so the coarse propagator becomes progressively less able to correct long-lived propagating modes.

![Original Figure 4.5: deterioration of Parareal on advection–diffusion and Burgers' equation as viscosity falls](assets/papers/time-parallelization/source-figures/figure-4-5.svg)

Figure 4.5 uses the same three viscosities as Figure 4.4; panel (a) confirms the modal prediction. Burgers' equation lacks an equally sharp modal theory, but panel (b) shows the same trend. Standard iteration typically diverges near $\nu\le10^{-3}$. The finite-step property remains algebraically true but loses practical value.

The mechanism the paper gives for the hyperbolic case is that arbitrarily small high-frequency components propagate arbitrarily far in space and time, so it is very hard to make $\mathcal G$ comparable to $\mathcal F$; and once $\mathcal G$ is made accurate enough, the coarse-grid correction itself becomes so expensive that the speed-up disappears. Analyses appear in Gander and Vandewalle (2007), Gander and Lunet (2020a,b), Gander, Lunet and Pogoželskytė (2023a), and Gander, Lunet, Ruprecht and Speck (2023b). This marks the intended boundary of the methods in this chapter.

Section 4.2 closes by noting that MGRiT is the multilevel generalization of Parareal and that considerable effort has gone into making it work for advection (Howse et al. 2019; De Sterck et al. 2021, 2023a, 2023b). One route is a semi-Lagrangian optimized coarse solver, which remains difficult in the nonlinear case; another is the diagonalization-based coarse solver of Gander and Wu (2020), covered in Section 4.5 of this chapter.

## Equation, theorem, and figure audit

| Source item                          | Paper section | Coverage                                                                 |
| ------------------------------------ | ------------- | ------------------------------------------------------------------------ |
| Section 4 introduction and 4.1       | 4, 4.1        | temporal locality, four method families, historical lines                |
| (4.1), Figure 4.1                    | 4.2           | update, parallel fine solves, sequential coarse correction, two grids    |
| (4.2)–(4.4), Theorem 4.1, Remark 4.1 | 4.2           | modal reduction, Toeplitz matrix, complete proof, preconditioner view    |
| (4.5), Theorem 4.2, Figure 4.2       | 4.2           | superlinear and long-time linear bounds, both regimes                    |
| (4.6), Theorem 4.3                   | 4.2           | nonlinear assumptions, bound, finite-step implication                    |
| (4.7)–(4.9), Theorem 4.4, Figure 4.3 | 4.2           | L/A stability, 0.3 factor, SDIRK comparison                              |
| Figures 4.4–4.5                      | 4.2           | advection–diffusion spectrum, Burgers experiment, weak-diffusion failure |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 4–4.2, pp. 443–452.
