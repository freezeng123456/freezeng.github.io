---
title: "3.1–3.2: Historical Development and Schwarz Waveform Relaxation"
description: Starting from waveform relaxation and domain decomposition, a complete derivation of OSWR, the Robin minimax parameter, finite-step convergence, and tent pitching, with the multi-subdomain Dirichlet theory filled in
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-3-1-history-and-swr
tags:
  - parallel-in-time
  - SWR
  - hyperbolic-PDE
---

> [!note] Reading scope
> This page corresponds to Sections 3.1–3.2 of the paper (pp. 396–405) and covers the historical development, equations (3.1)–(3.4), Theorems 3.1–3.2, and Figures 3.1–3.3. The source graphics are kept as-is; the exposition, derivation notes, and intuition are all newly written.

> [!info] Figure license
> The paper is distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). This page reproduces the complete graphics of Figures 3.1–3.3 and provides an English figure explanation and attribution in the adjacent text.

## 3.1 Historical development

### How the four lines of work took shape

The paper discusses together four classes of PinT methods that are effective for hyperbolic problems, and uses the four PDEs from Chapter 2 to check their theoretical properties. It is worth tracing the origins first, because the concurrency of each class comes from a different mathematical structure, and understanding the source is the key to understanding where each method excels and where it fails.

1. **Schwarz waveform relaxation (SWR)** arises from combining domain decomposition with waveform relaxation. The continuous space–time subproblem construction was first proposed for parabolic problems by Gander (1999), with Giladi and Keller (2002) independently giving a similar construction. Its two roots are: domain decomposition (DD), a classical technique for parallel PDE solving that can be traced back to Schwarz (1870); and waveform relaxation (WR), which originated in the circuit simulation of Lelarasmee et al. (1982). Gander (1997) developed and analyzed both the parabolic and hyperbolic cases, and the name SWR was formally introduced by Gander et al. (1999). Later results for nonlinear parabolic problems appear in Gander (1999) and Gander and Rohde (2005). The optimized version using more effective transmission conditions (OSWR) was developed in the parabolic case by Gander and Halpern (2007), Bennequin, Gander and Halpern (2009), and Bennequin, Gander, Gouarin and Halpern (2016), and in the hyperbolic case by Gander, Halpern and Nataf (2003) and Gander and Halpern (2004); nonlinear advection–diffusion is treated in Gander, Lunowa and Rohde (2023). There are also Dirichlet–Neumann and Neumann–Neumann waveform relaxation variants (Gander, Kwok and Mandal 2016b, 2021b). The unmapped tent pitching (UTP) of Ciaramella, Gander and Mazzieri (2023) is also built on SWR.
2. **Parallel integral deferred correction (IDC)** descends from the correction idea for evolution problems in Böhmer and Stetter (1984). Dutt, Greengard and Rokhlin (2000) identified it as a dedicated time integrator that raises the order round by round by treating the integral term exactly. Guibert and Tromeur-Dervout (2007) proposed the pipelined version PIDC, and Christlieb et al. (2010) proposed the revisionist version RIDC; both parallelize in the time direction.
3. **ParaExp** was proposed by Gander and Güttel (2013). It handles initial-value propagation and the source-term response separately, which subsequently led to linear implementations (Merkel et al. 2017; Kooij et al. 2017) and a nonlinear iterative extension (Gander, Güttel and Petcu 2018a).
4. **ParaDiag** developed from the direct diagonalization method of Maday and Rønquist (2008) as a direct, iteration-free parallel-in-time solver. Later work extended it to parabolic (Gander et al. 2016a), hyperbolic (Gander et al. 2019), and nonlinear problems (Gander and Halpern 2017), and produced implementations based on waveform relaxation (Gander and Wu 2019), Parareal (Gander and Wu 2020), Krylov preconditioning (McDonald et al. 2018; Liu and Wu 2020), interpolation (Kressner, Massei and Zhu 2023), and Sherman–Morrison–Woodbury plus Krylov (Gander and Palitta 2024).

The role of this history is to mark the source of concurrency. SWR exploits complete waveforms on the spatial interfaces, IDC exploits a pipeline between correction levels, ParaExp exploits linear superposition, and ParaDiag exploits the spectral structure of the all-at-once time matrix. Because the four rely on different structures, their behavior on hyperbolic and parabolic problems also differs systematically—precisely what the following sections verify one by one.

> [!tip] Insight: why SWR comes first
> Among the four classes, only SWR builds "parallelism" directly on the physical propagation mechanism of the PDE: it first cuts the domain in space, then lets each subdomain solve the entire time window independently, and exchanges time functions on the artificial interfaces. For hyperbolic problems, the finite speed of propagation makes this interface exchange converge exactly in a finite number of steps (see 3.2.2); for parabolic problems, the infinite speed of propagation makes optimized transmission conditions necessary for acceleration. This "fit" between the method and the PDE structure is the fundamental reason it spans both parabolic and hyperbolic problems.

## 3.2 Schwarz waveform relaxation (SWR)

### From WR to SWR

SWR fuses the strengths of classical Schwarz DD and WR while overcoming the limitations of each. To see this clearly, one must first understand separately why each predecessor is insufficient.

### The limitation of classical spatial domain decomposition

In the context of evolution PDEs, classical Schwarz DD typically first performs a uniform implicit time discretization and then, at each time step, uses DD techniques to solve the resulting elliptic problem sequentially (Cai 1991; Meurant 1991; Cai 1994). This carries two structural costs: first, all subdomains must wait until the DD iteration at the current time step has converged before advancing together to the next step, so the time direction is fully serialized; second, all subdomains must share the same time discretization, which erases exactly the most valuable advantage of DD—tailoring the numerical treatment separately for each subdomain (different spatial grids, time steps, or integrators).

### The splitting difficulty of waveform relaxation

Classical waveform relaxation takes a different starting point: from the first-order ODE system obtained by spatial discretization of the evolution PDE,

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol f(t),
\qquad A=M+N
$$

it uses a dynamic iteration that resembles Picard iteration but with a block splitting of the system:

$$
\frac{d\boldsymbol u^k}{dt}-M\boldsymbol u^k
=N\boldsymbol u^{k-1}+\boldsymbol f,
\qquad \boldsymbol u^k(0)=\boldsymbol u_0.
$$

Here $k\geq1$ is the iteration index, $\boldsymbol u^k(0)=\boldsymbol u_0$ holds for all $k$, and $(M,N)$ is a consistent splitting of $A=M+N$. The way of splitting directly determines the parallel structure:

- **Jacobi (diagonal) type $M$**: solving for $\boldsymbol u^k$ reduces to a set of mutually independent scalar ODEs and can be fully parallelized. It becomes a PinT method because the "future" of every unknown is approximated before the future of its neighboring unknowns is known—the time axis is advanced ahead of schedule.
- **Gauss–Seidel (triangular) type splitting**: this likewise reduces to scalar ODEs and can obtain parallelism through red–black or other colorings.
- **Cyclic reduction**: this can add further levels of parallelism (Worley 1991; Horton, Vandewalle and Worley 1995; Simoens and Vandewalle 2000).

The real difficulty of WR lies not in the parallelism but in **finding an effective splitting that guarantees fast convergence**. The convergence rate depends entirely on whether $A=M+N$ keeps the strong coupling in $M$ and leaves only the weak coupling to $N$. As Nevanlinna (1989) emphasized, what one really cares about is "what kind of subsystem partition makes the iteration converge fast … how to split so that the coupling stays weak is an important question." This sentence pinpoints the crux: the splitting is usually treated as "given," but a bad splitting throws the strongly coupled terms into $N$, making convergence arbitrarily slow or even divergent, which renders WR unusable in practice.

> [!tip] Insight: how SWR "sidesteps" the splitting difficulty
> The key step in SWR is to change the object being decoupled. Rather than first discretizing in space and then searching for a good algebraic splitting, it **first cuts the spatial domain at the continuous level** and then, like WR, solves the continuous space–time PDE independently on each subdomain. The advantage is this: since we know exactly through what mechanism neighboring subdomains couple under a specific PDE, we can directly design **transmission conditions** to approximate that coupling and make the iteration converge fast. So Nevanlinna's difficulty of "how to split so that the coupling is weak" is replaced by "how to design transmission conditions so that the subproblems are approximately decoupled"—and the latter has the physical structure of the PDE to rely on, so it is solvable. This is exactly the starting point of OSWR.

SWR partitions the continuous domain before spatial discretization. Each subdomain solves the PDE independently over the complete time window and exchanges a time function at the artificial boundary. This both permits choosing different space–time discretizations for each subdomain and allows the transmission conditions to be designed according to the PDE's propagation mechanism. Optimized Robin, Ventcel, or convolution conditions try to approximate the Dirichlet-to-Neumann (DtN) map of the continuous problem: the DtN map is precisely the operator that condenses the entire influence of a subdomain's exterior on the interface, and if the transmission condition matches it exactly, the iteration converges exactly in one step, so the closer the transmission condition is to the DtN map, the faster the convergence. OSWR has been used for many PDEs, for example the shallow-water equations (Martin 2009), time-domain Maxwell (Courvoisier and Gander 2013), Schrödinger (Halpern and Szeftel 2010; Besse and Xing 2017; Antoine and Lorin 2017), the primitive equations of the ocean (Audusse, Dreyfuss and Merlet 2010), quantum wave problems (Antoine and Lorin 2016), fractional diffusion (Wu 2017), and coupled Ekman boundary layers (Thery et al. 2022).

![Schwarz waveform relaxation exchanges interface waveforms over a complete time window](assets/diagrams/pint/en/schwarz-waveform-relaxation.svg)

OSWR has completely different convergence properties for first-order parabolic problems (such as the advection–diffusion equation (2.5) and the nonlinear Burgers equation (2.6)) and for second-order hyperbolic problems (such as the wave equation (2.7)), which are discussed separately below.

## 3.2.1 First-order parabolic problems

### The complete OSWR iteration

Consider the advection–diffusion equation on the interval $(0,L)$ with homogeneous Dirichlet boundaries $u(0,t)=u(L,t)=0$ and initial value $u(x,0)=u_0(x)$, and write

$$
\mathcal L=\partial_x-\nu\partial_{xx}.
$$

Take the overlapping subdomains $\Omega_1=(0,\beta L)$ and $\Omega_2=(\alpha L,L)$, where $0<\alpha\leq\beta<1$. Here $\alpha,\beta$ are dimensionless coordinates relative to the interval length $L$, so $\beta-\alpha$ is the **overlap fraction in units of the interval length**, while the true **physical overlap width** is $l=(\beta-\alpha)L$. The paper's Robin OSWR is written as

$$
\left\{
\begin{aligned}
\partial_tu_1^k+\mathcal Lu_1^k&=0,
&& (x,t)\in\Omega_1\times(0,T],\\
u_1^k(0,t)&=0,\\
\frac1p\partial_xu_1^k(\beta L,t)+u_1^k(\beta L,t)
&=\frac1p\partial_xu_2^{k-1}(\beta L,t)+u_2^{k-1}(\beta L,t),
\end{aligned}
\right.
$$

$$
\left\{
\begin{aligned}
\partial_tu_2^k+\mathcal Lu_2^k&=0,
&& (x,t)\in\Omega_2\times(0,T],\\
\frac1p\partial_xu_2^k(\alpha L,t)-u_2^k(\alpha L,t)
&=\frac1p\partial_xu_1^{k-1}(\alpha L,t)-u_1^{k-1}(\alpha L,t),\\
u_2^k(L,t)&=0.
\end{aligned}
\right. \tag{3.1}
$$

Every round restores the physical initial value: $u_i^k(x,0)=u_0(x)$. The two interface waveforms at $k=0$, namely $\{u_1^0(\alpha L,t),u_2^0(\beta L,t)\}$, may be chosen arbitrarily. The parameter $p>0$ determines the Robin transmission strength: the left subdomain on its right interface $x=\beta L$ and the right subdomain on its left interface $x=\alpha L$ each use one Robin condition to bring in the neighbor's previous-round solution (a combination of value and derivative). As $p\to\infty$, the $\tfrac1p\partial_x$ terms vanish and the Robin condition reduces to the classical Dirichlet exchange $u_1^k(\beta L,t)=u_2^{k-1}(\beta L,t)$ and $u_2^k(\alpha L,t)=u_1^{k-1}(\alpha L,t)$. The multi-subdomain generalization repeats the same construction along every artificial boundary, directly and without extra theory; the nonlinear version keeps the iteration structure and merely replaces $\mathcal L$ by the corresponding nonlinear operator.

### Theorem 3.1: the minimax choice of the Robin parameter

Theorem 3.1 (Gander and Halpern 2007) analyzes (3.1) under the simplifying assumptions of an **unbounded spatial domain, two subdomains, and continuous space and time**. Let $l>0$ be the physical overlap width, and set

$$
y_0=\frac{l}{\nu},
\qquad
y=\frac{l\omega}{\nu},
\qquad
\omega\in\left[\frac\pi T,\frac\pi{\Delta t}\right].
$$

The optimized parameter has the scaling $p^*=\widetilde p^*/\nu$ (an equivalent form of the paper's original expression $p^*=\widetilde p^*\nu/s$). Define the single-frequency convergence factor

$$
R_0(y,\widetilde p,y_0)
=\frac{(y-\widetilde p)^2+y^2-y_0^2}
{(y+\widetilde p)^2+y^2-y_0^2}e^{-y},
$$

and the interior extremum location

$$
\bar y(y_0,\widetilde p)
=\sqrt{\frac{
y_0^2+2\widetilde p
+\sqrt{\widetilde p\left(-\widetilde p^3-4\widetilde p^2
+(4+2y_0^2)\widetilde p+8y_0^2\right)}}{2}}.
$$

$R_0$ is the OSWR convergence factor obtained in Fourier space, where $y$ corresponds to a single Fourier mode $\omega\in[\pi/T,\pi/\Delta t]$, that is $y=l\omega/\nu$. The goal of optimizing $p$ is to make the maximum of $R_0$ over the entire relevant frequency interval as small as possible, so it splits into two cases:

If $y_0<y_c$, where $y_c=1.618386576\ldots$, then $\widetilde p^*$ is the unique solution of

$$
R_0(y_0,\widetilde p^*,y_0)
=R_0\!\left(\bar y(y_0,\widetilde p^*),\widetilde p^*,y_0\right). \tag{3.2a}
$$

If $y_0\geq y_c$, then $\widetilde p^*$ is uniquely determined by

$$
y_0=\widetilde p^*\sqrt{\frac{\widetilde p^*}{4+\widetilde p^*}} \tag{3.2b}
$$

Setting

$$
y_{\min}=\frac{l\pi}{\nu T},
\qquad
y_{\max}=\frac{l\pi}{\nu\Delta t},
$$

then with the optimized Robin parameter $p^*$, the worst convergence factor over all relevant frequencies satisfies

$$
\rho:=\max_{y\in[y_{\min},y_{\max}]}
R_0(y,\widetilde p^*,y_0)
\leq
R_0\!\left(\bar y(y_0,\widetilde p^*),\widetilde p^*,y_0\right). \tag{3.2c}
$$

Classical Dirichlet exchange corresponds to $p=\infty$. Substituting $p=\infty$ into $R_0$ (the rational coefficient tends to $1$) leaves only the exponential decay across the overlap, so

$$
\rho\leq e^{-y_{\min}}
=\exp\!\left(-\frac{l\pi}{\nu T}\right). \tag{3.3}
$$

> [!tip] Derivation intuition: why the equal-peak condition and two branches appear
> After a Fourier transform in time, every frequency $\omega$ becomes an independent spatial interface-error propagation problem. Crossing the overlap gives the decay $e^{-y}$, and the Robin condition gives the rational reflection coefficient in front, $\big[(y-\widetilde p)^2+y^2-y_0^2\big]/\big[(y+\widetilde p)^2+y^2-y_0^2\big]$. Optimizing $p$ is equivalent to pushing down the highest peak over the entire frequency interval (a minimax problem). When the two candidate peaks are equally high, lowering one necessarily raises the other; equation (3.2a) is exactly this "equal-peak" (equioscillation) balance condition: the peak at the interval endpoint $y_0$ equals the peak at the interior extremum $\bar y$. The value $y_0=y_c\approx1.618$ is a structural dividing point: when the overlap relative to diffusion is large enough ($y_0\geq y_c$), the interior extremum is no longer an active constraint, and the optimal $p$ is instead determined by the simpler algebraic relation (3.2b). Equation (3.3) then shows that Dirichlet relies only on overlap decay, with a factor varying with $l/(\nu T)$, clearly weaker than the Robin case.

### Figure 3.1: the distance between continuous theory and discrete experiment

![Source Figure 3.1: theoretical OSWR convergence factors and four-subdomain iteration counts](assets/papers/time-parallelization/source-figures/figure-3-1.svg)

The paper uses $L=8.2$, $T=5$, $\Delta t=0.01$, $\Delta x=0.02$, $l=2\Delta x$, centered differences in space, backward Euler in time, and the initial value

$$
u_0(x)=e^{-10(x-L/2)^2}.
$$

Figure 3.1(a) plots the theoretical convergence factors of the Dirichlet and optimized Robin conditions as a function of the diffusion parameter $\nu$. The smaller $\nu$ is, the more the advection term dominates, the smaller the factor, and the faster the method converges. Figure 3.1(b) divides $(0,L)$ into four subdomains, starts from a random interface guess, and stops when the difference between the iterate and the converged solution is smaller than $10^{-8}$. The measured iteration count likewise decreases as $\nu$ decreases, in good agreement with the theoretical prediction in (a).

At $\nu=0.1$, the measured Dirichlet and optimized Robin counts are 92 and 28, respectively; whereas the two-subdomain continuous theory predicts $\rho$-based counts of 32 and 4, clearly smaller. The discrepancy comes from three changes of setting: the theory uses an **unbounded domain**, the experiment uses a bounded one; the theory has only **two subdomains**, the experiment has four; and the experiment also discretizes space and time. On this basis the paper cautions that Figure 3.1(a) describes trends and ideal factors and cannot serve directly as an iteration-count formula for the multi-subdomain discrete implementation.

The current state of the multi-subdomain theory must be stated on two levels: **a multi-subdomain convergence analysis for Dirichlet transmission conditions does exist**, see Gander and Stuart (1998) and Wu, Huang and Huang (2012); **whereas a complete convergence analysis for Robin transmission conditions in the multi-subdomain case is still missing**. For Robin conditions, a two-subdomain semi-discrete convergence analysis appears in Wu and Al-Khaleel (2014); a detailed comparison between continuous and discrete, and between bounded and unbounded domains, for the steady-state case appears in Gander, Halpern, Hubert and Krell (2020, 2021a).

> [!tip] Insight: reading Figure 3.1 correctly
> Figure 3.1(a) is the ideal "two-subdomain, unbounded, continuous" factor, while Figure 3.1(b) is the real "four-subdomain, bounded, discrete" iteration count, and the two are separated by three simplifications. So the correct reading is: the theoretical curves are used to judge **trends and relative merit** (Robin is systematically better than Dirichlet, and small $\nu$ is faster), while the **absolute iteration count must be taken from the discrete experiment**. This also explains why, beyond filling in the multi-subdomain Dirichlet theory (Gander and Stuart 1998; Wu, Huang and Huang 2012), it remains necessary to point out the gap in the multi-subdomain Robin theory—the latter is precisely the largest remaining unclosed distance between theory and Figure 3.1(b).

### More accurate transmission conditions

Beyond Dirichlet and Robin, one can further accelerate OSWR with Ventcel conditions (Bennequin et al. 2016). Ventcel conditions are essentially a local high-order approximation of the optimal transmission condition—the latter being given in Fourier (or Laplace) space by Gander and Halpern (2007), of the form

$$
\partial_x-\frac{1}{2\nu}
\mathcal F^{-1}\!\left(1+\sqrt{1+4i\nu\omega}\right),
$$

where $i^2=-1$ and $\mathcal F^{-1}$ is the inverse transform with respect to the Fourier mode $\omega$. This operator is nonlocal (it contains a square root and couples all frequencies) and cannot be implemented directly, so Ventcel uses a local differential operator to approximate it. In the asymptotic sense, if $l=C_1\Delta x$, $\Delta t=C_1\Delta x^\beta$, and $\Delta x$ is small, the convergence factor satisfies $\rho=1-O(\Delta x^\gamma)$, with the exponent $\gamma>0$ determined by $\beta$ (Gander and Halpern 2007; Bennequin et al. 2009, 2016). This means the finer the mesh, the closer the factor is to $1$, and the slower the convergence becomes. The convolution transmission condition of Wu and Xu (2017) instead gives a **mesh-independent** constant factor $\rho=1-C$, $C\in(0,1)$, thereby avoiding this degeneration with mesh refinement; it is especially suited to evolution equations with nonlocal terms, such as Volterra-type partial integro-differential equations. The cost is that the convolution interface operator must retain the time history.

> [!tip] Insight: the trade-offs among the three tiers of transmission conditions
> From Dirichlet to Robin to Ventcel/convolution is a spectrum that becomes "ever closer to the DtN map, and ever more expensive." Dirichlet is simplest but relies only on overlap decay; Robin uses a single optimizable scalar $p$ to buy significant acceleration; Ventcel uses a local high-order operator to approximate the optimal operator more closely, but its factor still tends to $1$ under mesh refinement; the convolution condition buys a mesh-independent constant factor at the cost of retaining the time history. Which tier to choose depends on how much interface complexity and memory one can afford, and on whether the mesh will be significantly refined.

## 3.2.2 Second-order hyperbolic problems

### Dirichlet conditions can also converge in a finite number of steps

Unlike first-order parabolic problems, for second-order hyperbolic problems (such as the wave equation (2.7)), SWR converges in a **finite number of steps** even with simple Dirichlet transmission conditions. On two overlapping subdomains the SWR reads

$$
\left\{
\begin{aligned}
\partial_{tt}u_1^k&=c^2\partial_{xx}u_1^k+g,
&& (x,t)\in\Omega_1\times(0,T],\\
u_1^k(x,0)&=u_0(x),
&\partial_tu_1^k(x,0)&=\widetilde u_0(x),\\
u_1^k(0,t)&=0,
&u_1^k(\beta L,t)&=u_2^{k-1}(\beta L,t),
\end{aligned}
\right.
$$

$$
\left\{
\begin{aligned}
\partial_{tt}u_2^k&=c^2\partial_{xx}u_2^k+g,
&& (x,t)\in\Omega_2\times(0,T],\\
u_2^k(x,0)&=u_0(x),
&\partial_tu_2^k(x,0)&=\widetilde u_0(x),\\
u_2^k(\alpha L,t)&=u_1^{k-1}(\alpha L,t),
&u_2^k(L,t)&=0.
\end{aligned}
\right. \tag{3.4}
$$

Here $c>0$ is the wave speed and $0<\alpha<\beta<1$. Both subdomains restore the physical initial values $u_0,\widetilde u_0$, and only on the artificial interfaces $x=\beta L$ and $x=\alpha L$ do they exchange the neighbor's previous-round interface values in Dirichlet fashion.

**Theorem 3.2 (Gander 1997, Theorem 6.3.3).** For $0<\alpha<\beta<1$, as long as

$$
k>\frac{Tc}{\beta-\alpha},
$$

the error at the two artificial interfaces at round $k$ is zero throughout $(0,T)$:

$$
u_1^k(\alpha L,t)-u(\alpha L,t)=0,
\qquad
u_2^k(\beta L,t)-u(\beta L,t)=0.
$$

This conclusion comes from the finite speed of propagation inherent to hyperbolic problems. The error itself satisfies a wave equation with zero source and zero initial data, so it can only propagate along characteristics at speed $c$. Each round of interface exchange can only push the "known-to-be-correct" region inward along the characteristic cone by one physical width $l=(\beta-\alpha)L$; the time height corresponding to this width is about $(\beta-\alpha)L/c$. After $k$ rounds, the cumulative time height by which the correct region has advanced is about $k(\beta-\alpha)L/c$. Note that here $\beta-\alpha$ is the overlap as a fraction of the interval length, and $T$ has been made dimensionless in the same way, so when the cumulative height exceeds the entire time window, the condition is exactly $k>Tc/(\beta-\alpha)$, after which the whole interface waveform is already exact. Corresponding results hold for many subdomains, higher-dimensional decompositions (Gander and Halpern 2004), and one-dimensional nonlinear conservation laws (Gander and Rohde 2005).

> [!tip] Insight: the watershed between parabolic and hyperbolic lies in the propagation speed
> The contrast between Theorem 3.1 and Theorem 3.2 is the pivot of the whole chapter. Parabolic problems have infinite propagation speed, so any interface error instantly affects the entire subdomain, and one can only get linear convergence of the form "multiply by a factor $<1$ each round," with a factor depending on $\nu,l,T$. Hyperbolic problems have finite propagation speed, so the error is locked inside the characteristic cone, and each round "exactly" enlarges a correct region, giving finite-step convergence whose step count depends only on the geometric quantity $Tc/(\beta-\alpha)$. The larger the overlap fraction $\beta-\alpha$, the more each round advances and the fewer steps are needed—which directly inspires the red–black construction below that trades "large overlap" for parallelism.

### Figure 3.2: how red–black SWR advances the correct region

![Source Figure 3.2: the four geometric stages of red–black SWR with large overlap](assets/papers/time-parallelization/source-figures/figure-3-2.svg)

Theorem 3.2 reveals an exploitable property (already pointed out in Gander et al. 2003, Figure 3.1): within the characteristic cone that is "influenced only by the initial data and not yet contaminated by possibly wrong interface data," each subdomain computes exactly the exact solution. On this basis one can deliberately choose space–time subdomains and use **directly parallel subdomain solves** instead of repeated iteration. Figure 3.2 uses five interleaved subdomains, paired with generous overlap, and computes in the following order.

1. Solve in parallel on the three red subdomains $\Omega_1,\Omega_3,\Omega_5$ up to $T_1$. The interior interfaces located at $x_2,x_4$ are unknown at this point and may first be given arbitrary values. By the finite speed of propagation, the triangular tent at the bottom of each subdomain bounded by characteristics is already exact; the two end subdomains, whose outer boundaries are known, each gain an additional small exact region. This step also computes a "not-yet-correct" approximation above the correct tents—precisely the redundant computation that Nievergelt advocated in exchange for parallelism.
2. Solve in parallel on the two black subdomains $\Omega_2,\Omega_4$ up to $T_2$. Since the blue regions from the previous round are already correct, they obtain correct interface data from them, thereby producing two correct rhombic tents, again accompanied by some redundant computation.
3. Return to the red subdomains and advance the time to the interval $(T_1,T_3)$.
4. The red–black alternation continues, the correct region is pushed up along time, and eventually the coverage reaches $T_4$.

Each subdomain also computes a temporarily unreliable region above the correct tent, and this redundant work buys parallelism. Nievergelt's idea reappears here in the geometric form of characteristic cones: it is better to compute some regions destined to be discarded than to prevent the subdomains from advancing independently of one another within the same step.

### MTP, UTP, and Figure 3.3

The red–black SWR above is in fact one simple way to realize one of the most powerful current hyperbolic space–time solvers—mapped tent pitching (MTP, Gopalakrishnan, Schöberl and Wintersteiger 2017; for the time-domain Maxwell application see Gopalakrishnan et al. 2020). MTP maps the inclined tents that appear in red–black SWR into space–time cylinders, solves inside the cylinders with classical time stepping, and maps back to the original geometry, thereby **avoiding redundant computation**. There are two costs: first, one must compute the mapping itself; second, after mapping the size of the computational domain is comparable to the space–time subdomains in red–black SWR, so the computational cost is comparable. More troublesome still, the mapping introduces **order reduction**, which requires a specially designed time integrator to compensate.

By contrast, red–black SWR (now also called unmapped tent pitching, UTP) keeps the original coordinates and therefore **suffers no order reduction**, and it can be implemented very simply: it is equivalent to applying the restricted additive Schwarz (RAS) technique from DD directly to the all-at-once space–time system (see Gander 2008 for the explanation), an equivalence that lets higher-dimensional implementations reuse mature RAS infrastructure.

![Source Figure 3.3: the round-by-round elimination of wave-equation error regions by UTP](assets/papers/time-parallelization/source-figures/figure-3-3.svg)

Figure 3.3 shows red–black SWR / UTP solving the wave-equation model problem (2.7). Panel (a) is the initial error produced by a random interface guess; panels (b), (c), and (d) show the fourth red update, the eighth black update, and the twelfth red update, respectively. One sees that UTP constructs the exact solution inside the red and black tents without needing to know the tent structure in advance, advancing in exactly the same way as MTP. The blue zero-error region rises tent by tent. UTP can likewise be applied easily to nonlinear hyperbolic problems; if the tent height is not known in advance, one need only observe the residual of the computed solution—the height at which the residual drops to zero in the time direction naturally indicates the current tent height, from which the time-window length $T_i-T_{i-1}$ can be chosen adaptively.

The original MTP relies on the finite speed of propagation and therefore **cannot be applied directly to parabolic problems**: parabolic problems have infinite propagation speed, and there is no tent "within which the solution is exact." Nevertheless SWR and UTP can still be used for parabolic equations, where the optimized SWR variants (Gander and Halpern 2007; Bennequin et al. 2009) are especially effective. For weakly diffusive, advection-dominated problems (such as the advection–diffusion model here), one may consider using UTP just the same, only adding one or two extra iterations within each time slab to correct the influence across tents.

> [!tip] Insight: UTP = RAS on the all-at-once system
> The most practical point about UTP is that it reduces "tent pitching"—a seemingly special kind of space–time advancement—to an object people already know very well: restricted additive Schwarz on the all-at-once space–time system. This brings three benefits: no need to explicitly construct a mapping, no order reduction, and the ability to directly borrow existing RAS code to scale to higher dimensions. Combined with the trick of "adapting the tent height from the height at which the residual drops to zero," UTP retains the exact advancement structure of MTP while avoiding its mapping cost and order issues; for weakly diffusive parabolic problems, it can be carried over with only one or two extra iterations per slab.

## Equation and figure coverage check

| Source item                         | Paper section    | Coverage status                                                                                                                   |
| ----------------------------------- | ---------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| historical development, pp. 396–398 | 3.1              | origins, principal branches, and key references of all four families listed, including DN/NN variants                             |
| WR splitting and SWR motivation     | 3.2 introduction | continuous decomposition limits, $A=M+N$ splitting parallelism, bad-splitting dilemma, transmission bypass                        |
| (3.1)                               | 3.2.1            | two-subdomain PDEs, physical boundaries, Robin exchange, the $p\to\infty$ limit, and initial conditions                           |
| (3.2a)–(3.2c), (3.3), Theorem 3.1   | 3.2.1            | parameter scaling, single-mode factor, the two $y_c$ branches, minimax equal-peak, worst-frequency bound, Dirichlet limit         |
| Figure 3.1                          | 3.2.1            | complete source figure, parameters, stopping test, threefold discrepancy, and the state of Dirichlet/Robin multi-subdomain theory |
| Ventcel and convolution conditions  | 3.2.1            | optimal DtN operator form, asymptotic factor, mesh-independent convolution factor, and nonlocal cost                              |
| (3.4), Theorem 3.2                  | 3.2.2            | two-subdomain iteration, finite-step condition $k>Tc/(\beta-\alpha)$, finite-propagation-speed argument                           |
| Figures 3.2–3.3                     | 3.2.2            | complete source figures, red–black progression, Nievergelt redundancy, MTP/UTP, RAS equivalence, residual adaptivity              |

## Source of this page

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 3.1–3.2, pp. 396–405.
