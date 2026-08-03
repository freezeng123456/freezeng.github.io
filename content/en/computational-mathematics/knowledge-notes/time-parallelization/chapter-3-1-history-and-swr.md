---
title: "3.1–3.2: Historical Development and Schwarz Waveform Relaxation"
description: A complete account of OSWR, finite-step convergence, and tent pitching from waveform relaxation and domain decomposition
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-3-1-history-and-swr
tags:
  - parallel-in-time
  - SWR
  - hyperbolic-PDE
---

> [!note] Reading scope
> This page covers Sections 3.1–3.2 of the paper (pp. 396–405), including the historical development, equations (3.1)–(3.4), Theorems 3.1–3.2, and Figures 3.1–3.3. The source graphics are unchanged. The exposition, derivation notes, and intuition are newly written.

> [!info] Figure license
> The paper is distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). This page reproduces the complete graphics in Figures 3.1–3.3 and supplies an independent explanation next to each one.

## Section 3 introduction: how this group is positioned

Before turning to individual algorithms, the paper states three judgements. First, PinT methods designed for hyperbolic problems generally work just as well for parabolic problems, often better. Second, mapped tent pitching is the exception, because it is built on finite propagation speed, which parabolic problems do not have. Third, the converse fails: the methods designed for parabolic problems in Chapter 4 generally break down on hyperbolic problems.

## 3.1 Historical development

### How the four lines of work developed

The paper groups four PinT families that remain effective for hyperbolic problems and tests their main theoretical properties on the four PDEs from Chapter 2.

1. **Schwarz waveform relaxation (SWR)** combines domain decomposition and waveform relaxation. Gander (1999) proposed the continuous space–time subproblem construction, and Giladi and Keller (2002) developed a related approach independently. The domain-decomposition lineage reaches back to Schwarz (1870), while waveform relaxation originated in the circuit simulations of Lelarasmee et al. (1982). Gander (1997) treated parabolic and hyperbolic problems, and Gander et al. (1999) introduced the name SWR. Later work covered nonlinear parabolic equations, optimized transmission conditions, hyperbolic equations, and Dirichlet–Neumann and Neumann–Neumann variants. The unmapped tent pitching method of Ciaramella et al. (2023) is also based on SWR.
2. **Parallel integral deferred correction (IDC)** descends from the correction process for evolution problems in Böhmer and Stetter (1984). Dutt et al. (2000) formulated IDC as an integrator whose order can increase with successive corrections. Guibert and Tromeur-Dervout (2007) proposed PIDC, and Christlieb et al. (2010) proposed RIDC.
3. **ParaExp**, introduced by Gander and Güttel (2013), separates initial-value propagation from the source response. Subsequent work developed practical linear implementations and a nonlinear iterative extension.
4. **ParaDiag** began with the direct diagonalization method of Maday and Rønquist (2008). It was extended to parabolic, hyperbolic, and nonlinear problems and later appeared as waveform relaxation, a Parareal component, a Krylov preconditioner, an interpolation method, and a Sherman–Morrison–Woodbury plus Krylov construction.

The history identifies the source of concurrency in each family. SWR exchanges complete interface waveforms. IDC pipelines correction levels. ParaExp uses linear superposition. ParaDiag exposes the spectral structure of the all-at-once time matrix.

## 3.2 Schwarz waveform relaxation (SWR)

### Limitation of classical spatial domain decomposition

A traditional algorithm first fixes a common time discretization and then applies domain decomposition to the elliptic problem at every time step. All subdomains wait for the DD iteration at the current step to converge before advancing. Sharing one time grid also restricts local adaptivity and the use of different integrators in different subdomains.

Classical waveform relaxation starts from

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol f(t),
\qquad A=M+N,
$$

and iterates

$$
\frac{d\boldsymbol u^k}{dt}-M\boldsymbol u^k
=N\boldsymbol u^{k-1}+\boldsymbol f,
\qquad \boldsymbol u^k(0)=\boldsymbol u_0.
$$

A Jacobi-type $M$ permits simultaneous scalar or block ODE solves. A Gauss–Seidel splitting can use red–black coloring, and cyclic reduction exposes further concurrency. Convergence depends completely on whether the splitting $A=M+N$ retains the strong couplings in $M$. Nevanlinna (1989) identified the construction of weakly coupled subsystems as the central practical problem in WR. A poor splitting may converge arbitrarily slowly or diverge.

SWR partitions the continuous spatial domain before discretization. Each subdomain solves its PDE over a complete time window and exchanges a time-dependent function on the artificial boundary. This construction supports different space–time discretizations on different subdomains and allows the transmission operator to reflect the PDE's propagation mechanism. Optimized Robin, Ventcel, and convolution conditions approximate increasingly complete Dirichlet-to-Neumann information.

![Schwarz waveform relaxation exchanges interface waveforms over a complete time window](assets/diagrams/pint/en/schwarz-waveform-relaxation.svg)

## 3.2.1 First-order parabolic problems

### The full OSWR iteration

For advection–diffusion on $(0,L)$, define

$$
\mathcal L=\partial_x-\nu\partial_{xx}.
$$

Let $\Omega_1=(0,\beta L)$ and $\Omega_2=(\alpha L,L)$ with $0<\alpha\leq\beta<1$. Their overlap has width $l=(\beta-\alpha)L$. The Robin OSWR iteration in the paper is

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

Every iteration restores the physical initial condition, $u_i^k(x,0)=u_0(x)$. The two interface waveforms at $k=0$ may be chosen arbitrarily. The parameter $p>0$ controls Robin transmission. As $p\to\infty$, the derivative terms vanish and the method reduces to classical Dirichlet exchange. A multi-subdomain version repeats the same construction at every artificial boundary. A nonlinear version keeps the iteration pattern and replaces $\mathcal L$ by the corresponding nonlinear operator.

### Theorem 3.1: the minimax Robin parameter

Theorem 3.1, taken from Gander and Halpern (2007), analyzes (3.1) on an unbounded spatial domain, with two subdomains and continuous space and time. Let

$$
y_0=\frac{l}{\nu},
\qquad
y=\frac{l\omega}{\nu},
\qquad
\omega\in\left[\frac\pi T,\frac\pi{\Delta t}\right].
$$

The optimized parameter scales as

$$
p^*=\frac{\widetilde p^*s}{\nu}.
$$

Define the single-frequency factor

$$
R_0(y,\widetilde p,y_0)
=\frac{(y-\widetilde p)^2+y^2-y_0^2}
{(y+\widetilde p)^2+y^2-y_0^2}e^{-y}
$$

and the interior extremum

$$
\bar y(y_0,\widetilde p)
=\sqrt{\frac{
y_0^2+2\widetilde p
+\sqrt{\widetilde p\left(-\widetilde p^3-4\widetilde p^2
+(4+2y_0^2)\widetilde p+8y_0^2\right)}}{2}}.
$$

If $y_0<y_c$, where $y_c=1.618386576\ldots$, then $\widetilde p^*$ is the unique solution of

$$
R_0(y_0,\widetilde p^*,y_0)
=R_0\!\left(\bar y(y_0,\widetilde p^*),\widetilde p^*,y_0\right). \tag{3.2a}
$$

If $y_0\geq y_c$, it is instead determined uniquely by

$$
y_0=\widetilde p^*\sqrt{\frac{\widetilde p^*}{4+\widetilde p^*}}. \tag{3.2b}
$$

With

$$
y_{\min}=\frac{l\pi}{\nu T},
\qquad
y_{\max}=\frac{l\pi}{\nu\Delta t},
$$

the worst factor over the relevant frequencies satisfies

$$
\rho:=\max_{y\in[y_{\min},y_{\max}]}
R_0(y,\widetilde p^*,y_0)
\leq
R_0\!\left(\bar y(y_0,\widetilde p^*),\widetilde p^*,y_0\right). \tag{3.2c}
$$

Dirichlet exchange is the limit $p=\infty$, which yields

$$
\rho\leq e^{-y_{\min}}
=\exp\!\left(-\frac{l\pi}{\nu T}\right). \tag{3.3}
$$

> [!warning] Source check: the symbol $s$ in $p^*$
> The journal version prints the scaling as $p^*=\widetilde p^*s/\nu$, but $s$ is never defined anywhere else in the paper; reading it as the unit advection speed of model problem (2.5) reduces it to $p^*=\widetilde p^*/\nu$. The arXiv preprint prints the reciprocal factor, $p^*=\widetilde p^*\nu/s$. This page follows the journal version and records the discrepancy here so that readers comparing the two are not confused.

> [!tip] Derivation intuition: why the theorem uses an equal-peak condition
> A temporal Fourier transform turns every frequency $\omega$ into an independent spatial interface-error problem. Crossing the overlap contributes $e^{-y}$, while Robin transmission contributes the rational reflection factor in front of it. Optimizing $p$ minimizes the largest peak over the frequency interval. Once the two candidate peaks have equal height, lowering either one raises the other; equation (3.2a) states this minimax balance. The value $y_0=y_c$ marks a change in the interior-peak structure. Note that the equal-peak condition holds at $y=y_0$ and $y=\bar y$, whereas (3.2c) maximizes over $[y_{\min},y_{\max}]$; the paper does not comment on the fact that these intervals do not coincide exactly.

### Figure 3.1: continuous theory versus a discrete experiment

![Source Figure 3.1: theoretical OSWR factors and four-subdomain iteration counts](assets/papers/time-parallelization/source-figures/figure-3-1.svg)

The experiment uses $L=8.2$, $T=5$, $\Delta t=0.01$, $\Delta x=0.02$, and $l=2\Delta x$. Space is discretized by centered differences and time by backward Euler. The initial value is

$$
u_0(x)=e^{-10(x-L/2)^2}.
$$

Figure 3.1(a) plots the theoretical factors for Dirichlet and optimized Robin transmission. As $\nu$ decreases and advection becomes dominant, both factors fall. Figure 3.1(b) divides the interval into four subdomains, starts from a random interface guess, and stops when the error relative to the converged solution falls below $10^{-8}$. The measured iteration count follows the same trend.

At $\nu=0.1$, the measured Dirichlet and optimized Robin counts are 92 and 28. The continuous two-subdomain theory predicts 32 and 4. Three differences explain the gap: the theory uses an unbounded domain, the experiment uses a bounded one; the theory has two subdomains, the experiment has four; and the experiment discretizes both space and time. Figure 3.1(a) therefore predicts a trend and an idealized factor, not an exact count for the discrete multi-subdomain implementation.

The theoretical coverage is uneven. Multi-subdomain convergence analyses for Dirichlet transmission conditions appear in Gander and Stuart (1998) and Wu, Huang and Huang (2012); the two-subdomain semi-discrete Robin case is treated in Wu and Al-Khaleel (2014), and steady-state variants in Gander, Halpern, Hubert and Krell (2020, 2021a). A comprehensive convergence analysis for Robin transmission conditions in the multi-subdomain case is still lacking, a gap the paper states explicitly.

### More accurate transmission conditions

Ventcel conditions (Bennequin et al. 2016) use a higher-order local operator to approximate optimal nonlocal transmission. The Fourier-space optimal operator quoted in the paper (Gander and Halpern 2007, Section 3) has the form

$$
\partial_x-\frac{1}{2\nu}
\mathcal F^{-1}\!\left(1+\sqrt{1+4i\nu\omega}\right),
$$

where $i^2=-1$. If $l=C_1\Delta x$ and $\Delta t=C_1\Delta x^\beta$, the asymptotic factor has the form $\rho=1-O(\Delta x^\gamma)$, with $\gamma$ determined by $\beta$. The convolution condition of Wu and Xu (2017) gives a mesh-independent factor $\rho=1-C$, $C\in(0,1)$, and applies naturally to evolution equations with Volterra-type nonlocal terms. Its interface operator must retain temporal history.

## 3.2.2 Second-order hyperbolic problems

### Finite-step convergence with Dirichlet data

For the wave equation, two-subdomain SWR is

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

**Theorem 3.2** (Gander 1997, Theorem 6.3.3)**.** For $0<\alpha<\beta<1$, if

$$
k>\frac{Tc}{\beta-\alpha},
$$

then both artificial-interface errors vanish over the complete interval $(0,T)$:

$$
u_1^k(\alpha L,t)-u(\alpha L,t)=0,
\qquad
u_2^k(\beta L,t)-u(\beta L,t)=0.
$$

Finite propagation speed drives the result. One exchange moves correct interface data only a finite distance along characteristics. Every iteration adds another exact space–time cone. Once their cumulative height exceeds $T$, the complete interface waveform is exact. Higher-dimensional decompositions are covered in Gander and Halpern (2004) and one-dimensional nonlinear conservation laws in Gander and Rohde (2005).

The paper then draws the key consequence. Since each subdomain already computes the exact solution inside the cone influenced only by the initial condition, one can choose the space–time subdomains so as to **avoid iterations entirely** and advance directly with parallel space–time subdomain solves. This observation goes back to Gander et al. (2003). Tent pitching, described next, is how that route is realized.

### Figure 3.2: red–black SWR advances exact regions

![Source Figure 3.2: four geometric stages of red–black SWR with generous overlap](assets/papers/time-parallelization/source-figures/figure-3-2.svg)

Figure 3.2 uses five interleaved subdomains.

1. Solve concurrently on the red subdomains $\Omega_1,\Omega_3,\Omega_5$ up to $T_1$. The unknown interior interfaces may initially receive arbitrary data. Finite propagation makes a triangular region at the bottom of each subdomain exact; the known physical boundaries contribute two additional exact side regions.
2. Solve concurrently on the black subdomains $\Omega_2,\Omega_4$ up to $T_2$. Their boundary data comes from the exact blue regions, producing two exact rhombi.
3. Return to the red subdomains and advance to $T_3$, lifting the apex of the exact region to $T_3$.
4. Panel (d) solves on the black subdomains over $(T_2,T_4)$; alternating colors in the same way lifts the exact region layer by layer.

Each solve also computes a region above the exact tent whose data is not yet reliable. This redundant work buys concurrency and realizes Nievergelt's idea through characteristic geometry.

### MTP, UTP, and Figure 3.3

Mapped tent pitching (MTP; Gopalakrishnan, Schöberl and Wintersteiger 2017, with a time-domain Maxwell application in Gopalakrishnan, Hochsteger, Schöberl and Wintersteiger 2020) maps an inclined tent to a space–time cylinder, applies a classical time stepper, and maps the solution back. It avoids the redundant regions of red–black SWR but adds the cost of computing the mapping; the paper states explicitly that after the mapping the computational domains have the same size as the space–time subdomains of red–black SWR, so the **computational cost is comparable** — MTP is not cheaper for having removed the redundancy. The mapping can also reduce the observed order, which motivates specialized integrators.

Unmapped tent pitching (UTP; Ciaramella, Gander and Mazzieri 2023) executes red–black SWR directly and can also be viewed as restricted additive Schwarz on the all-at-once system (Gander 2008). It retains the original coordinates, advances exactly like MTP, avoids mapping-induced order reduction, and is just as easily applied to nonlinear hyperbolic problems. Higher-dimensional implementations can reuse established RAS infrastructure.

![Source Figure 3.3: successive elimination of wave-equation error regions by UTP](assets/papers/time-parallelization/source-figures/figure-3-3.svg)

Figure 3.3 solves model problem (2.7), whose solution is shown in Figure 2.4(d). It starts from the error produced by a random interface guess. Panels (b), (c), and (d) show the fourth red update, the eighth black update, and the twelfth red update. The blue zero-error region rises tent by tent. An implementation need not identify the tents in advance. The time at which the residual first ceases to be zero gives the effective tent height and can be used to adapt $T_i-T_{i-1}$.

MTP requires finite propagation speed and does not extend directly to a parabolic equation. SWR and UTP still apply to parabolic problems, particularly with optimized transmission. A weakly diffusive advection problem may need one or two additional iterations in every time slab to correct influence across the nominal tents.

## Equation and figure coverage

| Source item                         | Paper section    | Coverage                                                                                                  |
| ----------------------------------- | ---------------- | --------------------------------------------------------------------------------------------------------- |
| historical development, pp. 396–398 | 3.1              | origins and principal branches of all four families                                                       |
| WR splitting and SWR motivation     | 3.2 introduction | continuous decomposition, splitting difficulty, and OSWR objective                                        |
| (3.1)                               | 3.2.1            | PDEs, physical boundaries, Robin exchange, and initial conditions for both subdomains                     |
| (3.2a)–(3.2c), (3.3), Theorem 3.1   | 3.2.1            | parameter scaling, single-mode factor, both parameter regimes, worst-frequency bound, and Dirichlet limit |
| Figure 3.1                          | 3.2.1            | complete source figure, parameters, stopping test, and theory–experiment discrepancy                      |
| Ventcel and convolution conditions  | 3.2.1            | optimal operator, asymptotic factor, and nonlocal cost                                                    |
| (3.4), Theorem 3.2                  | 3.2.2            | two-subdomain iteration, finite-step condition, and interface-error conclusion                            |
| Figures 3.2–3.3                     | 3.2.2            | complete source figures, red–black progression, MTP/UTP, and residual adaptivity                          |

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 3.1–3.2, pp. 396–405.
