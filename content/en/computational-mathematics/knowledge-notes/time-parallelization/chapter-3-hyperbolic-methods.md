---
title: "Chapter 3: Effective PinT Methods for Hyperbolic Problems"
description: A Chapter 3 overview—the challenge of long-range propagation in hyperbolic problems, the shared principles and historical lineages of the four hyperbolic-effective PinT methods (SWR, IDC, ParaExp, ParaDiag), together with complete derivations and numerical interpretation
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-3-hyperbolic-methods
tags:
  - parallel-in-time
  - hyperbolic-PDE
---

> [!note] Reading scope
> This chapter corresponds to Section 3 of the paper (pp. 396–443). The main text strictly retains the original paper's hierarchy: the Section 3 introduction, 3.1 historical development, 3.2 SWR, 3.3 IDC, 3.4 ParaExp, and the two ParaDiag families in 3.5.1/3.5.2. Equations, theorems, and paper experiments are explained in the order of the argument; Python results, parameter comparisons, and coverage audits are marked separately as site supplements and do not take up paper section numbers.

## Source-to-page map

This page retains the chapter-level overview and the site's reproduction experiments. The complete equation-by-equation, theorem-by-theorem, and figure-by-figure derivations are split into the following pages, for reading in the paper's order:

| Paper section      | Source pages | Close-reading page                                                                                                                                        | Coverage                                                          |
| ------------------ | ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| Section 3, 3.1–3.2 | pp. 396–405  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-1-history-and-swr\|Historical development and Schwarz waveform relaxation]] | (3.1)–(3.4), Theorems 3.1–3.2, Figures 3.1–3.3                    |
| 3.3                | pp. 405–412  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-2-idc\|Parallel integral deferred correction]]                              | (3.5)–(3.12), Theorem 3.3, Figures 3.4–3.6                        |
| 3.4                | pp. 412–415  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-3-paraexp\|ParaExp]]                                                        | (3.13)–(3.21), Theorem 3.4, Figures 3.7–3.8                       |
| 3.5, 3.5.1         | pp. 415–431  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i\|Direct ParaDiag]]                                             | (3.22)–(3.48), Theorems 3.5–3.7, Figures 3.9–3.14, Tables 3.1–3.2 |
| 3.5.2              | pp. 431–443  | [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-5-paradiag-ii\|Iterative ParaDiag]]                                         | (3.49)–(3.68), Theorems 3.8–3.9, Figures 3.15–3.18                |

## Section 3 introduction and 3.1 Historical development

### From parabolic to hyperbolic: why temporal locality disappears

Chapter 2 showed that the solution of a parabolic problem is highly local in time (Figure 2.1): under Dirichlet boundaries almost all information is quickly "forgotten," while under Neumann and periodic boundaries only the lowest-frequency component (the constant) survives over long times. Once a transport term appears and dominates (Figures 2.2–2.3), and all the more so in the hyperbolic limit, the situation changes qualitatively—for the second-order wave equation (Figure 2.4) every frequency component can propagate arbitrarily far in space and time, accompanied by multidirectional propagation and reflection. It is precisely this "propagation of fine information over long times" that makes time parallelization far harder for hyperbolic problems than for parabolic ones, and that dictates the use of PinT techniques different from the parabolic methods.

> [!tip] Insight
> "Temporal locality" versus "temporal globality" can serve as the first criterion for choosing a PinT method. The dissipation of a parabolic problem rapidly damps distant coupling, so a coarse propagator only needs to capture the low frequencies to drive a Parareal-type iteration (see Chapter 4); a hyperbolic problem has no such damping, and any coarse model that loses high frequencies or phase becomes distorted over a long time window. The shared goal of this chapter's four families is therefore not to "construct a good, strongly dissipative coarse propagator," but to "make the inversion of the all-at-once coupling system $K\boldsymbol u=\boldsymbol b$ spanning the entire time interval parallelizable in itself."

### Why this group of methods can handle long-range propagation

The paper groups four families as "hyperbolic-effective": Schwarz waveform relaxation (SWR), parallel integral deferred correction (IDC), ParaExp, and ParaDiag. They have different origins but share one principle—rather than approximating a strongly dissipative coarse propagator, they recast the inversion of the all-at-once coupling system into a parallelizable form, relying respectively on iteration, decomposition, and transformation:

- **Iteration plus decomposition**: SWR decomposes the space–time domain into subdomains, lets each subdomain solve a full time window at once, and couples them only through iterated interface waveforms;
- **Error decomposition and pipelining**: IDC splits high-order accuracy into "a low-order prediction plus several integral-residual corrections," and pipelines different time windows and correction levels in parallel;
- **Exact decomposition**: ParaExp splits the linear problem into "local zero-initial-value forced responses" and "global homogeneous propagation," using the matrix exponential to jump directly to any later time;
- **Transform diagonalization**: ParaDiag performs an (approximate) diagonalization along the time direction, turning the inverse of the all-at-once matrix into several mutually independent complex-shifted spatial solves.

The paper stresses a noteworthy phenomenon: many methods designed for hyperbolic problems are equally effective, or even better, on parabolic problems; the sole exception is mapped tent pitching (MTP), introduced at the end of 3.2—its construction directly exploits the finite propagation speed of hyperbolic problems and is therefore genuinely "hyperbolic-only." Conversely, PinT methods designed for parabolic problems (Chapter 4) usually perform poorly on hyperbolic problems. This asymmetry is the fundamental reason for devoting a separate chapter to these four families.

Finite propagation speed is the physical fact this group of methods exploits repeatedly. SWR and tent pitching explain it through characteristic cones: one iteration can advance correct interface information only across a finite propagation distance, and after enough iterations the characteristic cones cover the entire space–time domain (see 3.2.2). ParaExp expresses homogeneous propagation exactly with the matrix exponential, at a cost that does not grow linearly with the number of intermediate time steps. The $\alpha$-circulant of ParaDiag weakens the periodic head–tail closure and so avoids carrying the head–tail mismatch produced by propagation throughout the entire time domain (see 3.5.2). From different angles, all four preserve the phase and high-frequency information on which hyperbolic solutions depend.

The equation-by-equation and theorem-by-theorem derivations of each lineage are split into the close-reading pages: [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-1-history-and-swr\|Historical development and Schwarz waveform relaxation]], [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-2-idc\|Parallel integral deferred correction]], [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-3-paraexp\|ParaExp]], [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i\|Direct ParaDiag]], and [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-5-paradiag-ii\|Iterative ParaDiag]].

### 3.1 Historical development of the four method lineages

The paper organizes this group of methods along four mutually independent lineages and validates each of their main theoretical properties on the four PDEs introduced in Chapter 2.

**Space–time continuous subproblems: from the confluence of DD and WR to SWR.** The first lineage originates in solving on overlapping or nonoverlapping "space–time continuous" subproblems, first proposed for parabolic problems by Gander (1999) and independently introduced by Giladi and Keller (2002). It merges two traditions: domain decomposition (DD), a classical technique for solving PDEs in parallel that traces back to Schwarz (1870); and waveform relaxation (WR), which first appeared in circuit simulation (Lelarasmee et al. 1982). Gander (1997) developed and analyzed both parabolic and hyperbolic problems simultaneously, and the name "Schwarz waveform relaxation (SWR)" was coined by Gander et al. (1999). Further results for nonlinear parabolic problems appear in Gander (1999) and Gander and Rohde (2005). Optimized SWR (OSWR) with more effective transmission conditions: for parabolic problems see Gander and Halpern (2007), Bennequin, Gander and Halpern (2009), and Bennequin, Gander, Gouarin and Halpern (2016); for hyperbolic problems see Gander, Halpern and Nataf (2003) and Gander and Halpern (2004); for nonlinear advection–diffusion see Gander, Lunowa and Rohde (2023c). The recent unmapped tent pitching (UTP) technique (Ciaramella, Gander and Mazzieri 2023) is built directly on SWR. There are also Dirichlet–Neumann and Neumann–Neumann waveform-relaxation variants (Gander, Kwok and Mandal 2016b, 2021b).

**Time parallelization of integral deferred correction.** Another lineage is based on the time parallelization of IDC. IDC for evolution problems was first introduced by Böhmer and Stetter (1984) and later identified by Dutt, Greengard and Rokhlin (2000) as a dedicated time integrator—by treating the corresponding integrals accurately it can in principle generate numerical solutions of arbitrarily high order. Revisionist IDC (RIDC, Christlieb et al. 2010) is a time-parallelizable technique; there is also the pipelined parallel version (PIDC) proposed by Guibert and Tromeur-Dervout (2007), detailed in 3.3.

**ParaExp.** The ParaExp method introduced in 3.4 was proposed by Gander and Güttel (2013) about a decade ago, and its core strategy is to treat the initial value and the source term separately; see also Merkel, Niyonzima and Schöps (2017) and Kooij, Botchev and Geurts (2017), with a nonlinear variant in Gander, Güttel and Petcu (2018a).

**The ParaDiag family.** ParaDiag is introduced in 3.5. Diagonalization-based time-parallel methods were first proposed by Maday and Rønquist (2008) as a noniterative direct time-parallel solver; Gander, Halpern, Ryan and Tran (2016a) carried out a more detailed study of parabolic problems, with a nonlinear variant in Gander and Halpern (2017) and hyperbolic problems in Gander, Halpern, Rannou and Ryan (2019). Iterative variants soon followed and can be embedded either in WR (Gander and Wu 2019) or in Parareal (Gander and Wu 2020). Approximate ParaDiag has also been used directly as a Krylov preconditioner for the all-at-once system by McDonald et al. (2018) and Liu and Wu (2020). A systematic survey of ParaDiag appears in Gander et al. (2019); the family's influence in the PinT field then expanded rapidly, with new developments such as the interpolatory direct ParaDiag of Kressner, Massei and Zhu (2023) and the variant of Gander and Palitta (2024) combining the Sherman–Morrison–Woodbury formula with Krylov techniques.

## 3.2 Schwarz waveform relaxation

### From classical waveform relaxation to space–time decomposition

For the linear system $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$, classical waveform relaxation first splits $A=M+N$ and then iterates

$$
\boldsymbol u^{k+1\prime}-M\boldsymbol u^{k+1}
=N\boldsymbol u^k+\boldsymbol g.
$$

If $M$ has a block-diagonal or colored structure, several components can be computed in parallel. This approach depends directly on the algebraic splitting; an unsuitable $M$ can converge very slowly or even diverge. Classical spatial domain decomposition usually solves an elliptic problem at each time step and also ties every subdomain to the same temporal discretization.

SWR decomposes the continuous spatial domain first and then lets each subdomain solve an entire time window at once. Neighboring subdomains exchange whole time-dependent waveforms over the overlap or interface. Each subdomain may use a discretization suited to its local equation and mesh, and the transmission condition can be designed from the PDE's propagation mechanism.

![Schwarz waveform relaxation exchanges interface waveforms over a complete time window](assets/diagrams/pint/en/schwarz-waveform-relaxation.svg)

The iterated unknown of this method is an interface function. Dirichlet conditions exchange only solution values; Robin, Ventcel, or convolution conditions additionally approximate normal fluxes and a fuller Dirichlet-to-Neumann map. Versions with optimized transmission conditions are usually called OSWR.

### 3.2.1 First-order parabolic problems: overlap width and the Robin parameter

The paper illustrates the mechanism on a one-dimensional advection–diffusion equation. The domain is split into two overlapping subdomains with overlap width $l$. Each iteration solves both complete time-window problems in parallel and imposes Robin conditions on the artificial boundaries,

$$
(\partial_x+p)u_1^{k+1}=(\partial_x+p)u_2^k,
\qquad
(\partial_x-p)u_2^{k+1}=(\partial_x-p)u_1^k. \tag{3.1}
$$

Here $p>0$ controls the transmission condition, and $p\to\infty$ corresponds to Dirichlet exchange. The error analysis is carried out in Laplace/Fourier frequency space. Theorem 3.1 supplies the optimized Robin parameter and a worst-frequency convergence estimate for the continuous problem. The Dirichlet case satisfies a bound of the form

$$
\rho_{\mathrm D}\leq
\exp\!\left(-\frac{l\pi}{\nu T}\right)
$$

The bound clearly displays three trends: more overlap accelerates convergence; a longer time window makes the cross-window coupling harder; and as $\nu$ decreases, directed transport carries information through the overlap faster, so SWR can in fact speed up in this particular model.

Figure 3.1 uses $L=8.2$, $T=5$, $\Delta t=0.01$, $\Delta x=0.02$, $l=2\Delta x$, and a Gaussian initial condition. The numerical curves confirm that convergence improves as viscosity decreases, and also show that Robin conditions clearly outperform Dirichlet ones. In the four-subdomain experiment at $\nu=0.1$, Dirichlet and optimized Robin transmission require about 92 and 28 iterations respectively, while the continuous two-subdomain theory predicts about 32 and 4. The gap comes from the differences between continuous/discrete operators, two/four subdomains, and boundary configurations, so the theoretical bound should not be taken as the exact iteration count for this discrete experiment.

Higher-order Ventcel conditions can improve asymptotic convergence further. If a temporal convolution is used to approximate an exact transparent boundary, a mesh-independent convergence factor is theoretically attainable, at the cost of a more complicated, temporally nonlocal interface operator.

### 3.2.2 Second-order hyperbolic problems: finite-step propagation

For the wave equation, the properties of SWR are more direct. Each iteration advances correct interface information into the neighboring subdomain by one finite propagation distance. Theorem 3.2 states that, using Dirichlet transmission on two overlapping subdomains, once

$$
k>\frac{Tc}{\beta-\alpha},
$$

the interface error over the whole time window vanishes after iteration $k$. Here $\beta-\alpha$ is the overlap expressed as a fraction of the domain length (a dimensionless fraction), the true physical overlap width being $(\beta-\alpha)L$, where $L$ is the domain length; $c$ is the wave speed. The result follows from finite propagation speed: one iteration can extend correctness only over the region covered by a characteristic cone, and after enough iterations the cones cover the complete space–time domain.

Figure 3.2 visualizes this process with characteristic cones. A portion of each subdomain solution already agrees with the exact solution, and the next iteration enlarges that correct region through the interface data. This geometric interpretation also motivates red–black SWR: neighboring space–time blocks are computed concurrently by color, allowing some redundant work in exchange for greater concurrency.

### Tent pitching, MTP, and UTP

Tent pitching builds inclined space–time elements according to finite propagation speed. Once data on a tent's lower boundary is available, its interior can be computed independently. Mapped tent pitching (MTP) maps each inclined tent to a regular cylinder so that standard solvers can be used; the mapping increases implementation cost and can also cause order reduction.

Unmapped tent pitching (UTP) retains the original space–time geometry and can be interpreted as red–black SWR or restricted additive Schwarz on the all-at-once system. A residual determines how far a given tent can still advance, dispensing with explicit mapping and the associated order loss. Parabolic equations have infinite propagation speed and cannot form strictly independent tents; when diffusion is small, SWR/UTP can still correct the cross-tent influence through a few iterations.

## 3.3 Time-parallel IDC

### The residual error equation of IDC

Consider the initial-value problem

$$
u'(t)=f(u(t),t),\qquad u(0)=u_0.
$$

For a current approximation $u^k(t)$, the integral defect or residual can be written as

$$
r^k(t)=u_0+\int_0^t f(u^k(s),s)\,ds-u^k(t). \tag{3.6}
$$

Let the true error be $e^k=u-u^k$. Substituting $u=u^k+e^k$ back into the integral equation gives the integral form of the error; differentiating in time yields

$$
e^{k\prime}(t)
=f(u^k(t)+e^k(t),t)-f(u^k(t),t)+r^{k\prime}(t). \tag{3.8}
$$

IDC first obtains a predicted solution with a low-order method, then discretely solves this error equation and updates $u^{k+1}=u^k+e^k$. The paper gives the discrete correction formula for a $\theta$ method and the weights produced by interpolatory integration, ultimately forming the nodal recurrence (3.11). The correction uses the integral residual, avoiding direct numerical differentiation of a noisy residual.

![Correction and pipeline structure of IDC, PIDC, and RIDC](assets/diagrams/pint/en/idc-pipeline.svg)

Theorem 3.3 states that if the base integrator has order $p$, then after $k$ corrections on $M$ uniform nodes the overall order reaches

$$
\mathcal O\!\left(\Delta t^{\min\{M,(k+1)p\}}\right).
$$

The correction order is eventually limited by the node count. Spectral deferred correction (SDC) with $J$ Gauss–Lobatto nodes can reach order $2J-1$; PFASST in Chapter 4 places SDC inside a multilevel time-parallel structure.

### Why standard IDC remains sequential

A long time interval is usually cut into several windows. After the prediction and several corrections are completed on one window, its endpoint value becomes the initial value of the next window. Standard IDC still processes windows sequentially, and each correction level also contains a nodal recurrence, so high-order accuracy by itself brings no automatic time parallelism.

### 3.3.1 Pipeline IDC (PIDC)

PIDC lets different windows execute correction sweeps of different numbers at the same time. Once the pipeline is full, the prediction level works on a later window, the first correction level works on the preceding window, and higher correction levels lag further behind. The degree of concurrency is roughly determined by the number of correction levels; the fill and drain phases reduce efficiency for short jobs.

This scheduling carries an accuracy risk. The initial value of a later window keeps changing as the previous window receives higher-order corrections. A higher correction level may start from an initial value that is not yet smooth or settled. IDC's order-raising theory requires sufficient temporal regularity; when the initial data is irregular, the diffusion is too weak, or the solution contains high-frequency structure, the expected high order is lost.

Figure 3.5 uses the periodic advection–diffusion equation with $\Delta x=1/64$, $T=3$, window length $\Delta T=0.1$, $M=5$ nodes per window, and backward Euler as the base integrator. A narrow source with $\sigma=1000$ produces low regularity, while $\sigma=5$ gives smoother input. The results separate into three regimes: at low regularity, repeated IDC/PIDC corrections cannot reliably achieve high order; with smooth data and stronger diffusion, the corrections are clearly effective; with smooth data but very small diffusion, long-lived high frequencies still spoil the ideal order-raising.

### 3.3.2 Revisionist IDC (RIDC)

RIDC maintains a sliding $M$-node window for each correction level. The levels advance across successive time steps like an assembly line, without waiting for a complete window to finish. This reduces global synchronization and is well suited to keeping different correction levels running continuously on several cores.

RIDC improves scheduling but does not remove the regularity requirement. Figure 3.6 again shows that low regularity and weak diffusion limit the accuracy gain from the correction levels. For hyperbolic problems, preserving high frequencies is both a physical advantage and a numerical hazard for IDC order-raising.

## 3.4 ParaExp

### Exact decomposition for linear problems

ParaExp constructs a direct time-parallel solution for

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol g(t)
$$

Partition $[0,T]$ into $N$ subintervals $[T_{n-1},T_n]$. The first step solves the zero-initial-value inhomogeneous problems in parallel on all subintervals,

$$
\boldsymbol v_n'(t)=A\boldsymbol v_n(t)+\boldsymbol g(t),
\qquad \boldsymbol v_n(T_{n-1})=0. \tag{3.13}
$$

The second step propagates each segment's endpoint contribution to later times through the homogeneous equation:

$$
\boldsymbol w_n'(t)=A\boldsymbol w_n(t),
\qquad \boldsymbol w_n(T_n)=\boldsymbol v_n(T_n). \tag{3.14}
$$

Linear superposition gives the exact reconstruction

$$
\boldsymbol u(t)=\boldsymbol v_j(t)+
\sum_{n=0}^{j-1}\boldsymbol w_n(t),
\qquad t\in[T_{j-1},T_j]. \tag{3.15}
$$

![ParaExp separates local forced responses from global homogeneous propagation](assets/diagrams/pint/en/paraexp-decomposition.svg)

The key point lies in the homogeneous tails:

$$
\boldsymbol w_n(t)=e^{(t-T_n)A}\boldsymbol v_n(T_n).
$$

The matrix-exponential action can jump directly to any later time; its cost depends on the matrix and the target accuracy and usually does not grow linearly with the number of intermediate time steps. Large sparse systems commonly use rational Krylov or polynomial/Chebyshev approximations, while small dense matrices can use scaling-and-squaring with Padé. A wave-equation experiment cited in the paper reached about 80% parallel efficiency; that figure depends on the specific exponential algorithm, partition, and hardware.

### Nonlinear extension and its limits

For

$$
\boldsymbol u'=A\boldsymbol u+B(\boldsymbol u)+\boldsymbol g, \tag{3.17}
$$

the paper continues to hand the linear homogeneous part to exponential propagation and places $B(\boldsymbol u)$ into iterative, inhomogeneous local problems. Theorem 3.4 gives two important conclusions: after iteration $k$, the first $k$ time subintervals already agree with the sequential fine solution; and at the coarse time points, the iteration is equivalent to a simplified Parareal, whose coarse propagator is the linear homogeneous evolution and whose fine propagator resolves the full nonlinearity.

Figure 3.8 compares ParaExp with standard Parareal on Burgers' equation, with spatial step $0.01$, $T=2$, fine step $0.01/20$, and standard-Parareal coarse step $0.01$. At large viscosity, linear diffusion is the dominant dynamics and ParaExp's coarse model is very effective; as viscosity decreases, nonlinear transport grows in importance and standard Parareal is faster for a while; at $\nu=0.02$ ParaExp diverges. Moving still closer to the hyperbolic limit, standard Parareal also fails. The linear ParaExp construction is strong, and the quality of the nonlinear split determines the applicable range of the extended version.

## 3.5 ParaDiag: diagonalization along time

### Two ParaDiag routes

The goal of ParaDiag is to turn the all-at-once coupling system into several independent spatial systems. The paper distinguishes two classes:

1. **ParaDiag-I** exactly diagonalizes a specially designed time discretization to form a direct solver;
2. **ParaDiag-II** approximates the original matrix with a circulant or $\alpha$-circulant time matrix and then performs a stationary iteration or Krylov preconditioning.

Both require a well-conditioned temporal eigenvector matrix and that every complex-shifted spatial system can be solved efficiently. Their three-stage structure is the same: transform in time, solve the shifted problems in parallel, and apply the inverse transform.

![ParaDiag time transform, independent spatial solves, and inverse transform](assets/diagrams/pint/en/paradiag-three-stage.svg)

### 3.5.1 Direct ParaDiag methods (ParaDiag-I)

#### Backward Euler on a geometric variable-step mesh

For the linear system (2.1), the all-at-once matrix of variable-step backward Euler has the form

$$
K=B\otimes I_x-I_t\otimes A. \tag{3.23}
$$

If the time steps grow geometrically as $\Delta t_n=\mu^{n-1}\Delta t_1$, then $B$ can be written as $B=VDV^{-1}$. Hence

$$
K^{-1}
=(V\otimes I_x)
(D\otimes I_x-I_t\otimes A)^{-1}
(V^{-1}\otimes I_x), \tag{3.25}
$$

and the middle block-diagonal system contains $N_t$ mutually independent shifted spatial problems.

The parameter $\mu=1+\rho$ exposes the core conflict of the direct method. A larger $\rho$ makes the step variation pronounced, with truncation error of order $\mathcal O(\rho^2)$; when $\rho$ is too small, $B$ approaches a nondiagonalizable Jordan structure and roundoff error is amplified like $\epsilon\rho^{-(N_t-1)}$. Theorem 3.5 gives the balance between the two and the scale of the optimal $\rho$. Figures 3.9–3.10 confirm the U-shaped error curve and also show that the number of usable time steps is very limited in double precision. Directly increasing $N_t$ soon lets conditioning and roundoff dominate.

#### The wave equation and the trapezoidal rule

The wave equation is first converted to a first-order system and then integrated with a variable-step trapezoidal rule to preserve its energy property. The corresponding all-at-once system is again diagonalizable. Theorem 3.6 still obtains the competition between truncation error and $\epsilon\rho^{-(N_t-1)}$ roundoff amplification. Figure 3.11 and Table 3.1 show that as $N_t$ grows the eigenvector condition number increases rapidly, and the error starts to rise beyond roughly $N_t>32$.

#### Boundary-value methods mitigate the ill-conditioning

A boundary-value method (BVM) uses centered differences at the first $N_t-1$ nodes and backward Euler at the last node to close the all-at-once system:

$$
\frac{\boldsymbol u_{n+1}-\boldsymbol u_{n-1}}{2\Delta t}
=A\boldsymbol u_n+\boldsymbol g_n,
$$

which, together with the terminal discretization, forms a diagonalizable matrix. Even though the terminal formula is first order, the whole scheme still attains second order. Theorem 3.7 gives $\operatorname{Cond}(V)=\mathcal O(N_t^2)$, far more stable than the geometric variable-step scheme. Figure 3.12 shows second-order error decay as the uniform time step shrinks, without the rapid deterioration seen earlier. For the second-order wave equation, one can also directly construct a system containing $B^2\otimes I_x-I_t\otimes A$, avoiding the doubled storage caused by first-order conversion.

#### Nonlinear ParaDiag-I

The nonlinear all-at-once system is first linearized by Newton. The true Jacobian differs at each time point, losing a single Kronecker structure. The paper approximates all time blocks with an average Jacobian to obtain a quasi-Newton system; after the time transform, the shifted Jacobian problems remain parallel. If the solution varies too much over a long time interval, the average Jacobian degrades, and one can switch to several sequential windows.

The Burgers experiments in Figure 3.13 and Table 3.2 show that at $\nu=0.1$ the number of parallel Jacobian solves is far below that of sequential Newton; as viscosity decreases, convergence becomes increasingly sensitive to the total time $T$; and at $\nu=0.002$ with a longer time window the method fails.

The approximate Jacobian can also be extended into a low-rank Kronecker sum

$$
J\approx\sum_{q=1}^r \Phi_q\otimes A_q. \tag{3.47}
$$

NKA chooses the temporal scaling matrices $\Phi_q$ offline on a coarse model, retaining more of the Jacobian's variation in time. Figure 3.14 shows that it clearly improves quasi-Newton convergence on longer windows such as $T=1.3$.

### 3.5.2 Iterative ParaDiag methods (ParaDiag-II)

#### Strang circulant preconditioner

The all-at-once system of a linear multistep method is written as

$$
K=B_1\otimes I_x-B_2\otimes\Delta t A.
$$

ParaDiag-II replaces the Toeplitz-like matrices $B_1,B_2$ with Strang circulant matrices $C_1,C_2$ to construct

$$
P=C_1\otimes I_x-C_2\otimes\Delta t A.
$$

Circulant matrices are simultaneously diagonalized by the discrete Fourier matrix, so applying $P^{-1}$ is again realized by "FFT, independent shifted spatial solves, inverse FFT." If $P$ is very close to $K$, a stationary iteration can be used directly; when the convergence factor is poor, $P$ is better suited as a GMRES preconditioner. Even when the spectral radius of the stationary iteration exceeds one, Krylov methods can still exploit clustered spectra for fast convergence.

Theorem 3.8 gives a structural result for symmetric negative definite $A$: the preconditioned matrix has only a finite set of nonunit eigenvalues, so GMRES has a finite-step upper bound in exact arithmetic. This bound grows with the spatial dimension and alone cannot guarantee fast convergence on large problems.

Figure 3.15 uses $T=2$, $\Delta t=1/50$, $\Delta x=1/100$. The circulant preconditioner is very effective for the heat equation and for advection–diffusion at $\nu=10^{-3}$; as $\nu$ falls further the spectral clustering worsens; for the wave equation the nonunit eigenvalues spread across the complex plane and the outer iterations increase markedly. This set of experiments connects directly to the conclusions of Chapter 2: dissipation localizes the head–tail mismatch caused by cyclic closure, while persistent propagation carries that mismatch throughout the entire time domain.

#### $\alpha$-circulants and a waveform-relaxation interpretation

The paper also derives the $\alpha$-circulant matrix from continuous-time head–tail waveform relaxation. The initial and terminal values of a new time window are weakly coupled through the parameter $\alpha$; after discretization, the eigenvector matrix has the form $\Lambda_\alpha F^*$ and can still use FFTs.

$\alpha=1$ corresponds to the standard circulant approximation. Periodic propagation problems can become singular at this value; $0<\alpha<1$ breaks the exact periodic head–tail closure and often significantly improves advection and wave problems. Figure 3.16 shows strong convergence without Krylov acceleration.

For stable one-step methods and symmetric two-step methods, Theorem 3.9 gives spectral bounds:

$$
\frac{1}{1+\alpha}
\le |\lambda(P_\alpha^{-1}K)|
\le \frac{1}{1-\alpha},
\qquad
\rho(I-P_\alpha^{-1}K)
\le \frac{\alpha}{1-\alpha}.
$$

These bounds are independent of the spatial grid and the number of time steps, but require the time integrator to satisfy the corresponding stability. A critical Numerov-parameter experiment shows that $\gamma=1/120$ satisfies the condition, whereas slightly crossing the stability bound to $1/120.01$ destroys the prediction.

Decreasing $\alpha$ improves the iteration factor but at the same time amplifies the roundoff error of the time transform to roughly $\epsilon/\alpha$. Figure 3.18 displays this tradeoff. Updating through an error equation avoids repeatedly bringing a large solution vector into an ill-conditioned transform, thereby reducing roundoff contamination. A more general multistep Volterra analysis likewise concludes that the eigenvalues deviate from one by $\mathcal O(\alpha)$.

#### Nonlinear ParaDiag-II

Nonlinear problems use Newton–Krylov. Each Newton linearization produces an all-at-once Jacobian, and GMRES is then preconditioned by $P_\alpha$ built from an average Jacobian. Stationary iteration often fails on long nonlinear windows, while Krylov methods can still exploit eigenvalue clustering. Shortening the time window makes the Jacobian more uniform, and NKA can also retain more temporal variation; both improve preconditioner quality.

## Site numerical supplement: Python validation of Figure 3.15

### Baseline experiment

The recomputed baseline uses $N_x=N_t=100$, $T=2$, $\nu=10^{-6}$, $\alpha=1$, and GMRES tolerance $10^{-12}$. The algorithm converges after 13 iterations with a true relative residual of $1.152\times10^{-14}$.

![GMRES convergence baseline for ParaDiag-II on advection–diffusion](assets/pint/paradiag-baseline.svg)

### Paper-grid validation

| Problem                            |                                         Python result | Interpretation                                                            |
| ---------------------------------- | ----------------------------------------------------: | ------------------------------------------------------------------------- |
| heat                               |                                      2 Krylov updates | eigenvalues tightly clustered near one                                    |
| advection–diffusion, $\nu=10^{-3}$ |                                                     3 | circulant preconditioner is close to the original all-at-once system      |
| advection–diffusion, $\nu=10^{-6}$ |                                                    13 | weak diffusion lets the head–tail closure error propagate for a long time |
| wave                               | preconditioned residual below $10^{-11}$ at update 89 | nonunit eigenvalues spread along $\operatorname{Re}\lambda=0.5$           |

![ParaDiag-II spectra and GMRES convergence for the heat, advection–diffusion, and wave equations](assets/pint/paradiag-figure-3-15.svg)

The ADE counts of 3 and 13 agree with Figure 3.15(c,d). The wave experiment uses the paper's $\gamma=1/100$ and $\alpha=1$; single-threaded SciPy reaches a preconditioned residual of $10^{-11}$ at update 89, close to the paper curve's endpoint of about 88. If stopped at SciPy's true relative residual of $10^{-12}$, it needs 103 updates. MATLAB and SciPy differ in residual normalization, restart, and stopping rules, so one should compare the spectral geometry and the convergence phase rather than a single iteration count.

## Site method comparison: parameters, applicability, and implementation cost

| Method      | Key parameter or choice                                  | Property determined                                   | Main risk                                                               |
| ----------- | -------------------------------------------------------- | ----------------------------------------------------- | ----------------------------------------------------------------------- |
| SWR         | overlap width, time window, transmission condition       | speed of characteristic information across subdomains | poor tuning or overly expensive interface operator                      |
| PIDC/RIDC   | node count, correction levels, window and pipeline depth | formal order and concurrency                          | loss of order under low regularity                                      |
| ParaExp     | exponential-action algorithm, linear/nonlinear split     | cost of homogeneous tail propagation                  | expensive matrix exponential or inaccurate nonlinear split              |
| ParaDiag-I  | time discretization, $N_t$, eigenvector conditioning     | direct concurrency scale                              | conflict between truncation and roundoff error                          |
| ParaDiag-II | $\alpha$, outer Krylov method, shifted-solve tolerance   | preconditioner clustering and stability               | small $\alpha$ amplifies roundoff, large $\alpha$ weakens approximation |

The public [ActaPinT-Python](https://github.com/freezeng123456/ActaPinT-Python) project already covers the Heat, ADE, and Wave experiments of Figure 3.15(a–f) and saves SVG, PNG, and JSON outputs. The upstream MATLAB repository also contains SWR, PIDC/RIDC, ParaExp, direct ParaDiag, and wave-domain-decomposition scripts; the current formal Python results do not yet cover these scripts, so this page does not fabricate numerical curves for them.

## Source coverage audit

| Source location                 | This page | Material covered                                                                                                                                                                                                                              |
| ------------------------------- | --------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Sections 3 and 3.1, pp. 396–398 | 3.1       | the challenge of long-range hyperbolic propagation and the temporal-locality criterion, the shared principle of the four families, the historical development and literature of the four lineages, MTP/UTP and finite-propagation-speed clues |
| Section 3.2, pp. 398–405        | 3.2       | WR versus SWR, Robin/Dirichlet/Ventcel/convolution conditions, Theorems 3.1–3.2, Figures 3.1–3.3, MTP/UTP                                                                                                                                     |
| Section 3.3, pp. 405–412        | 3.3       | IDC integral residual and recurrence, Theorem 3.3, SDC, PIDC/RIDC scheduling, Figures 3.4–3.6, regularity limitation                                                                                                                          |
| Section 3.4, pp. 412–415        | 3.4       | ParaExp linear decomposition and exponential action, nonlinear iteration, Theorem 3.4, Figures 3.7–3.8                                                                                                                                        |
| Sections 3.5–3.5.1, pp. 415–431 | 3.5.1     | ParaDiag-I three-stage method, geometric steps, Theorems 3.5–3.7, wave/BVM, nonlinear quasi-Newton and NKA, Figures 3.9–3.14, Tables 3.1–3.2                                                                                                  |
| Section 3.5.2, pp. 431–443      | 3.5.2     | Strang and $\alpha$-circulant preconditioners, Theorems 3.8–3.9, Figures 3.15–3.18, stability/roundoff tradeoff, nonlinear Newton–Krylov                                                                                                      |

## Summary

The four families of this chapter preserve long-range information in different ways. SWR propagates waveforms along interfaces and characteristics; IDC arranges high-order error corrections into a pipeline; ParaExp uses the matrix exponential to propagate the linear homogeneous response directly; and ParaDiag transforms the all-at-once coupling into parallel spatial solves. All of them avoid exclusive reliance on a strongly dissipative coarse propagator, and each carries its own constraint: SWR needs suitable transmission conditions, IDC needs temporal regularity, ParaExp depends on the exponential action and the linear split, and ParaDiag needs stable temporal diagonalization and scalable shifted solvers.

## Source of this chapter

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3, pp. 396–443.
