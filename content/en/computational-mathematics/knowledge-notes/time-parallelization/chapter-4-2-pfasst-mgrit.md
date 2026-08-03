---
title: "4.3–4.4: PFASST and MGRiT"
description: From collocation equations, SDC smoothing, and inter-level interpolation to FCF overlapping Parareal, a complete derivation of the matrix structure, convergence factors, and equal-cost comparison of both algorithm classes
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-2-pfasst-mgrit
tags:
  - parallel-in-time
  - PFASST
  - MGRiT
---

> [!note] Reading scope
> This page follows Sections 4.3–4.4 (pp. 452–460), covering equations (4.10)–(4.13), Theorems 4.5–4.6, and Figures 4.6–4.11. The collocation equations of PFASST, its inter-level transfer, its SDC approximation, and the overlapping interpretation of MGRiT are all expanded to an implementable matrix form, while the reasoning chain of "why it is built this way" is filled in step by step.

## 4.3 PFASST

### Origin and two node sets

PFASST (parallel full approximation scheme in space–time) was formally introduced by Emmett and Minion (2012). Its core idea had in fact appeared two years earlier: to reduce the cost of a single Parareal iteration, Minion (2010) replaced the expensive complete fine propagation $\mathcal F$ with **one SDC (spectral deferred correction) iteration**. SDC itself originates from Dutt, Greengard and Rokhlin (2000), a framework that successively approximates a high-order collocation solution through low-order sweeps. From the outset, then, PFASST carried the genetic trait of "trading a cheap approximate propagator for parallelism."

PFASST long lacked a clear description and theoretical analysis, because it superimposes three structures at once: the time-block Parareal iteration, the SDC sweeps inside each block, and the fine/coarse two-level collocation. Only when Bolten, Moser and Speck (2017), using the **algebraic representation of SDC** given by Minion et al. (2015), reinterpreted PFASST as a **temporal multigrid** method—and provided a convergence analysis in Bolten, Moser and Speck (2018)—did the algorithm acquire an analyzable algebraic skeleton. Gander et al. (2023b) further gave an exact characterization for the model problem in block-iteration form, and it is this formulation that this page adopts.

> [!tip] Insight
> Reading PFASST as "multilevel SDC" is the key to understanding it. A single SDC sweep is essentially a **low-order preconditioner solve** $\widetilde\Phi^{-1}$ of the dense collocation system $\Phi$ (see (4.11)), playing the role of the "smoother" in multigrid; while the fine/coarse two-level collocation plus the Lagrange transfer operators play the "coarse-grid correction." Thus PFASST = a time-block Parareal outer layer + a two-level SDC V-cycle inside each block. The triple loop can be parallelized because the most expensive fine-level sweep is decoupled across time blocks, and only the cheap coarse level propagates serially in time.

Split $(0,T)$ into $N_t$ large intervals $[T_n,T_{n+1}]$. On each interval, set up $M_f$ fine nodes and $M_c$ coarse nodes:

$$
t_{n,m}^{f}=T_n+\tau_m^f\Delta t,
\qquad
t_{n,m}^{c}=T_n+\tau_m^c\Delta t,
$$

where $\tau_0^{f,c}=0$, $\tau_{M_{f,c}}^{f,c}=1$, and $M_f>M_c$. Requiring $M_f>M_c$ lets the two levels play different roles: the fine level has many nodes and high collocation order, and is responsible for accuracy; the coarse level has few nodes and a small system, and is responsible for cheaply propagating the correction in time. The superscripts $f$ and $c$ mark the fine and coarse time grids, respectively.

### Collocation equation (4.10)

For $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$, at the $M$ nodes of either level, numerical quadrature (collocation) gives

$$
\boldsymbol u_{n,m}
=\boldsymbol u_{n,0}
+\Delta t\sum_{j=1}^{M}q_{m,j}
[A\boldsymbol u_{n,j}+\boldsymbol g(t_{n,j})],
\qquad m=1,\ldots,M. \tag{4.10}
$$

Here $\boldsymbol u_{n,j}$ is the approximation of $\boldsymbol u$ at $t=t_{n,j}$, with $M=M_f$ or $M_c$ and $t_{n,j}=t_{n,j}^f$ or $t_{n,j}^c$. Equation (4.10) is the result of approximating the interval integral $\int_{T_n}^{t_{n,m}}$ by the collocation weights $q_{m,j}$—each node value depends on **all** nodes, so this is a dense implicit system coupling the $M$ nodes together, and solving it directly is very expensive. This is precisely why SDC sweeps are used later to approximate the solution.

Let

$$
Q=(q_{m,j}),\qquad
\boldsymbol u_n=(\boldsymbol u_{n,1}^\top,\ldots,\boldsymbol u_{n,M}^\top)^\top,
$$

$$
\chi=
\begin{bmatrix}
0&\cdots&0&1\\
\vdots&&\vdots&\vdots\\
0&\cdots&0&1
\end{bmatrix},
\qquad \boldsymbol\chi=\chi\otimes I_x,
$$

so that the terminal node is copied as the initial value for all collocation nodes of the next interval:

$$
\boldsymbol u_n
=\Delta t(Q\otimes A)\boldsymbol u_n
+\boldsymbol\chi\boldsymbol u_{n-1}
+\Delta t\boldsymbol b_n.
$$

Here $\boldsymbol b_n=(Q\otimes I_x)\boldsymbol g_n$, and every row of the copying matrix $\chi$ is $(0,\ldots,0,1)$: it extracts the **terminal node** $\boldsymbol u_{n-1,M}$ of the previous interval (i.e. the solution at $t=T_n$, since $\tau_M=1$) and uses it as the common initial value $\boldsymbol u_{n,0}$ for all collocation nodes of the current interval. This "copying" action is the only serial coupling between blocks, and is exactly the Parareal-style sequential dependence expressed at the collocation level.

The fine and coarse levels are denoted, respectively,

$$
\Phi_f=I_f-\Delta tQ_f\otimes A,
\qquad
\Phi_c=I_c-\Delta tQ_c\otimes A,
$$

$$
\boldsymbol u_n^f
=\Phi_f^{-1}(\boldsymbol\chi_f\boldsymbol u_{n-1}^f+\Delta t\boldsymbol b_n^f),
\qquad
\boldsymbol u_n^c
=\Phi_c^{-1}(\boldsymbol\chi_c\boldsymbol u_{n-1}^c+\Delta t\boldsymbol b_n^c).
$$

where $I_c=I_{M_c}\otimes I_x$ and $I_f=I_{M_f}\otimes I_x$. The matrices $\Phi_f,\Phi_c$ are the system matrices of the collocation system: given the terminal value of the previous interval, solving $\Phi^{-1}$ solves the whole-block collocation solution on that interval. $\Phi_f$ is dense and large ($M_f$ nodes), the most expensive step in the whole pipeline, and precisely the object that needs to be approximately replaced by SDC sweeps.

### Lagrange inter-level transfer

To move data between the two collocation levels, one needs a prolongation $T_{c\to f}$ (coarse→fine) and a restriction $T_{f\to c}$ (fine→coarse). Both are defined by Lagrange interpolation: data on the coarse nodes spans an interpolation polynomial

$$
p^c(\tau;\boldsymbol u^c)
=\sum_{m=1}^{M_c}u_m^cL_m^c(\tau),
\qquad
L_m^c(\tau)=
\prod_{\substack{j=1\\j\ne m}}^{M_c}
\frac{\tau-\tau_j^c}{\tau_m^c-\tau_j^c}.
$$

Evaluating at the fine nodes $\{\tau_m^f\}$ gives the prolongation matrix $T_{c\to f}$; conversely, evaluating the fine Lagrange basis $L_m^f$ at the coarse nodes gives the restriction matrix $T_{f\to c}$. Both matrices are Kronecker-multiplied with $I_x$, so they perform polynomial interpolation only along the collocation-node direction and **do not touch the spatial variable at all**—which contrasts with spatial multigrid, where restriction/prolongation operators resample spatial degrees of freedom, and is a hallmark of "temporal multigrid."

> [!tip] Insight
> The transfer operators are chosen as Lagrange interpolation rather than simple injection or averaging because collocation nodes are usually taken to be non-equidistant Gauss-type points such as Radau/Gauss. Converting between the two levels using the interpolation polynomial defined by the nodes keeps the polynomial representation implied by the collocation solution consistent, avoiding the introduction of extra low-order truncation error that would contaminate the high-order fine solution during inter-level transfer.

In the block-iteration formulation of Gander et al. (2023b), PFASST is written as

$$
\boldsymbol u_{n+1}^{k+1}
=B_{01}\boldsymbol u_{n+1}^k
+B_{10}(\boldsymbol\chi\boldsymbol u_n^{k+1}+\Delta t\boldsymbol b_n)
+B_{00}(\boldsymbol\chi\boldsymbol u_n^k+\Delta t\boldsymbol b_n),
$$

$$
\begin{aligned}
B_{01}&=[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]
(I_f-\widetilde\Phi_f^{-1}\Phi_f),\\
B_{10}&=T_{c\to f}\Phi_c^{-1}T_{f\to c},\\
B_{00}&=[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]\widetilde\Phi_f^{-1}.
\end{aligned}
$$

The three subscripts follow the block-iteration convention: the first subscript marks whether the block acts on the current iterate of the present interval (subscript 1) or on a neighboring interval (subscript 0), and the second subscript distinguishes the new-iterate $k{+}1$ (subscript 0) from the old-iterate $k$ (subscript 1) input. Understanding its structure block by block reveals the multigrid skeleton of PFASST:

- The common factor $[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]$ is the **coarse-grid correction operator**. It restricts the fine-level residual to the coarse level via $T_{f\to c}$, cheaply solves the coarse collocation system with $\Phi_c^{-1}$, then prolongs back to the fine level via $T_{c\to f}$ and subtracts from the identity—exactly the standard "restrict–coarse-solve–prolong" structure of two-level multigrid.
- The factor $(I_f-\widetilde\Phi_f^{-1}\Phi_f)$ is the **iteration matrix of the fine-level SDC smoother**: $\widetilde\Phi_f^{-1}$ is one cheap sweep (preconditioner), and $(I_f-\widetilde\Phi_f^{-1}\Phi_f)$ is its error-propagation operator after one application. Thus $B_{01}$ = one fine-level SDC smoothing on the old iterate followed by coarse-grid correction, which is exactly a two-level V-cycle.
- $B_{10}=T_{c\to f}\Phi_c^{-1}T_{f\to c}$ uses the neighboring-interval data of the **new iterate** (including the copied initial value $\boldsymbol\chi\boldsymbol u_n^{k+1}$) to propagate one step forward on the coarse level—this serial chain is precisely the coarse-level "information highway" that makes the overall scheme converge in finitely many steps.
- $B_{00}$ uses the neighboring data of the **old iterate** for correction, and $\widetilde\Phi_f^{-1}$ again reflects the cheap fine sweep; its difference from $B_{10}$ plays the role of the difference between new and old coarse propagations in the Parareal correction.

Here $\widetilde\Phi_f$ is an easily solved approximation (preconditioner) of the fine collocation matrix.

> [!tip] Insight
> Comparing these three blocks with the classical Parareal correction $\mathcal G_{new}-\mathcal G_{old}+\mathcal F_{old}$, one sees that PFASST merely replaces the "complete fine propagation $\mathcal F$" with "fine-level SDC smoothing $\widetilde\Phi_f^{-1}$ + coarse-grid correction." This explains why it is cheaper than Parareal: a single iteration no longer requires solving the dense $\Phi_f^{-1}$, only one triangular sweep $\widetilde\Phi_f^{-1}$ plus a small coarse system $\Phi_c^{-1}$.

### SDC approximation (4.11) and Figure 4.6

In practice the preconditioner $\widetilde\Phi_f$ is constructed from implicit Euler between fine nodes:

$$
\frac{\boldsymbol u_{n,m+1}-\boldsymbol u_{n,m}}
{\Delta t(\tau_{m+1}^f-\tau_m^f)}
=A\boldsymbol u_{n,m+1}+\boldsymbol g(t_{n,m+1}^f),
\quad m=0,\ldots,M_f-1. \tag{4.11}
$$

that is,

$$
\widetilde\Phi_f=
\begin{bmatrix}
1&&&\\
-1&1&&\\
&\ddots&\ddots&\\
&&-1&1
\end{bmatrix}\otimes I_x
-\Delta t
\begin{bmatrix}
\tau_1^f-\tau_0^f&&&\\
&\tau_2^f-\tau_1^f&&\\
&&\ddots&\\
&&&\tau_{M_f}^f-\tau_{M_f-1}^f
\end{bmatrix}\otimes A.
$$

It is composed of a "diagonal block of node differences" times $I_x$, minus the Kronecker term of a "lower-bidiagonal difference matrix" with $\Delta tA$. The key point is that $\widetilde\Phi_f$ is **lower bidiagonal (block lower triangular)**: adjacent nodes advance segment by segment via implicit Euler, so $\widetilde\Phi_f^{-1}$ can be solved in one forward sweep, node by node, each node requiring only the solution of a spatial subsystem of the form $(I-\Delta t(\tau_{m+1}^f-\tau_m^f)A)$. This contrasts sharply with the dense $\Phi_f$: one $\widetilde\Phi_f^{-1}$ corresponds to one round of low-order SDC sweep, far cheaper than a full collocation solve, yet in the multigrid sense sufficient to serve as a high-frequency smoother.

> [!tip] Insight
> Why does the implicit-Euler difference make a good smoother? Because $\widetilde\Phi_f$ shares the same set of collocation nodes with $\Phi_f$, merely replacing the dense quadrature weights $Q_f$ with a lower-triangular "approximate integration" built from node differences. The difference between the two concentrates on the part describing the high-order coupling between nodes, exactly the error components that are high-frequency in space and high-order in time; implicit Euler strongly damps these components, so $(I_f-\widetilde\Phi_f^{-1}\Phi_f)$ can efficiently attenuate them, leaving the residual low-frequency error for the coarse-grid correction to handle. This is the origin of the "smoother + coarse correction" division of labor.

The numerical experiment takes $T=3$, periodic boundaries, zero initial data, the source (2.4) with $\sigma=1000$, $\Delta x=1/128$, and $\Delta t=1/64$. The fine level uses three-node Radau IIA ($M_f=3$):

$$
\left\{0,\frac{4-\sqrt6}{10},\frac{4+\sqrt6}{10},1\right\},
$$

and the coarse level uses two-node Radau IIA ($M_c=2$) $\{0,1/3,1\}$. Radau IIA is chosen because its right end contains the node $\tau_M=1$ (matching the copying matrix $\chi$ that extracts the interval's terminal value) and it is L-stable, suitable for stiff problems. The corresponding quadrature matrices are

$$
Q_f=
\begin{bmatrix}
\frac{88-7\sqrt6}{360}&\frac{296-169\sqrt6}{1800}&\frac{-2+3\sqrt6}{225}\\[2pt]
\frac{296+169\sqrt6}{1800}&\frac{88+7\sqrt6}{360}&\frac{-2-3\sqrt6}{225}\\[2pt]
\frac{16-\sqrt6}{36}&\frac{16+\sqrt6}{36}&\frac19
\end{bmatrix},
\qquad
Q_c=
\begin{bmatrix}
\frac{5}{12}&-\frac{1}{12}\\[2pt]
\frac34&\frac14
\end{bmatrix},
$$

from which the two transfer matrices are

$$
T_{c\to f}=
\begin{bmatrix}
1.2674&-0.2674\\
0.5325&0.4674\\
0&1
\end{bmatrix}\otimes I_x,
\qquad
T_{f\to c}=
\begin{bmatrix}
0.5018&0.6833&-0.1851\\
0&0&1
\end{bmatrix}\otimes I_x.
$$

These matrices are entirely determined by the node positions and Lagrange interpolation, containing no tunable parameters. Note that the last row of both transfer matrices is $(0,\ldots,0,1)$: the terminal node $\tau=1$ exists on both levels, so interpolation degenerates to the identity, guaranteeing that the interval's terminal value is transferred between levels without distortion.

![Source Figure 4.6: PFASST error on the heat equation and three viscous advection–diffusion equations](assets/papers/time-parallelization/source-figures/figure-4-6.svg)

The heat equation converges fastest. As the viscosity of the advection–diffusion equation decreases, the persistently propagating high frequencies become increasingly hard to represent on the coarse collocation level, and convergence worsens in the same direction. This is consistent with the coarse–fine propagation-mismatch mechanism observed earlier in Parareal and MGRiT: the coarse collocation level (with few nodes) cannot accurately represent the long-lived, nearly undamped high-frequency modes under weak diffusion, and the coarse-grid correction correspondingly fails.

## 4.4 MGRiT

### Three interpretations: as Parareal with overlap

Multigrid reduction in time (MGRiT) was introduced by Falgout et al. (2014) as yet another variant of Parareal. It admits at least three equivalent viewpoints:

1. **Algebraic multigrid + FCF relaxation**: view the time levels as the coarse/fine degrees of freedom of AMG, with F-relaxation updating fine points and C-relaxation updating coarse points; FCF is one composite F–C–F relaxation round;
2. **Block iteration** (Gander et al. 2023b; Gander and Lunet 2024, Chapter 4.6);
3. **Parareal with overlap** (Gander et al. 2018b, Theorem 4 and Corollary 1).

This page adopts the third viewpoint, directly giving the two-level FCF-relaxation update for the nonlinear ODE system (2.2) $\boldsymbol u'=\boldsymbol f(\boldsymbol u,t)$, $\boldsymbol u(0)=\boldsymbol u_0$:

$$
\boldsymbol u_0^{k+1}=\boldsymbol u_0,
\qquad
\boldsymbol u_1^{k+1}=\mathcal F(T_0,T_1,\boldsymbol u_0),
$$

$$
\begin{aligned}
\boldsymbol u_{n+1}^{k+1}
={}&\mathcal F\!\left(T_n,T_{n+1},
\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)\right)\\
&+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})\\
&-\mathcal G\!\left(T_n,T_{n+1},
\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)\right),
\quad n=1,\ldots,N_t-1. \tag{4.12}
\end{aligned}
$$

Here $\mathcal G,\mathcal F$ are the coarse and fine propagators of Parareal (see (4.1)). Comparing with the Parareal correction $\mathcal F(\boldsymbol u_{n-1}^k)+\mathcal G(\boldsymbol u_n^{k+1})-\mathcal G(\boldsymbol u_n^k)$, one sees the only difference: MGRiT replaces the corrected "previous-interval value" from $\boldsymbol u_{n-1}^k$ to **one fine propagation first**, $\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)$. Therefore (4.12) contains two fine propagations per round (the inner $\mathcal F$ plus the outer $\mathcal F$), whereas Parareal has only one per round.

![Source Figure 4.7: FCF-MGRiT is Parareal with an overlap width of one coarse time step](assets/papers/time-parallelization/source-figures/figure-4-7.svg)

In Figure 4.7 the dark circles are the coarse points where the coarse propagation $\mathcal G$ runs. The extra inner F-relaxation first pushes the old-iterate information forward by one coarse interval $\Delta T$, then enters the usual Parareal-style correction—geometrically equivalent to overlapping two adjacent Parareal "windows" by one coarse step. Hence two-level FCF-MGRiT is **Parareal with overlap width $\Delta T$**, and likewise converges exactly in finitely many steps: the global error vanishes after at most $k=\lceil N_t/2\rceil$ rounds (Gander et al. 2018b, Theorem 5). More generally, $F(CF)^\nu$ relaxation corresponds to Parareal with overlap width $\nu\Delta T$ (Gander et al. 2018b, Theorem 6 and Corollary 1 give the corresponding superlinear convergence results).

> [!tip] Insight
> Why does overlap reduce the exact-convergence step count from Parareal's $N_t$ to $\lceil N_t/2\rceil$? Parareal advances the "already exactly solved" temporal front by only one interval per round; whereas the extra F-relaxation of FCF makes the front cross one additional coarse interval per round, doubling the advance speed, so the number of rounds needed to reach global exactness is halved. The cost is one extra fine solve per round—this exactly trades "halving the number of iterations" for "doubling the per-round cost," and the two break even in total number of fine solves. This is the intuitive root of the equal-cost comparison below.

### Theorem 4.5: long-time modal factor

The convergence estimate in the linear case comes from Dobrev et al. (2017). Consider $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$, with $A$ diagonalizable and spectrum $\sigma(A)\subset\mathbb C^-$. Following the notation of Theorem 4.2, if the coarse propagator is stable, $|R_g(z)|<1$, then two-level FCF-MGRiT satisfies

$$
\max_n\|\boldsymbol e_n^k\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}\varrho_l^k(J,z)
\max_n\|\boldsymbol e_n^0\|_\infty,
$$

$$
\varrho_l(J,z)=
\frac{|R_f^J(z/J)|\,|R_g(z)-R_f^J(z/J)|}
{1-|R_g(z)|}. \tag{4.13}
$$

Here $R_g,R_f$ are the stability functions of the coarse propagator $\mathcal G$ and the fine propagator $\mathcal F$, $\boldsymbol e_n^k:=V_A(\boldsymbol u_n^k-\boldsymbol u_n)$, and $J=\Delta T/\Delta t$ is the number of fine steps within each coarse interval. $R_f^J(z/J)$ denotes continuously taking $J$ fine steps with step size $\Delta t=\Delta T/J$ (corresponding to parameter $z/J$), i.e. the entire fine propagation over one coarse interval. Comparing with the Parareal convergence factor (4.5b)

$$
\varrho_{l,\mathrm{Parareal}}
=\frac{|R_g(z)-R_f^J(z/J)|}{1-|R_g(z)|}
$$

immediately gives

$$
\varrho_{l,\mathrm{MGRiT}}
=|R_f^J(z/J)|\,\varrho_{l,\mathrm{Parareal}}.
$$

> [!tip] Insight
> Why does FCF multiply the Parareal factor by exactly one $|R_f^J(z/J)|$? Look at (4.12): what is fed into the "coarse correction difference $\mathcal G_{new}-\mathcal G_{old}$" is no longer the old error itself, but the error after a full fine propagation $\mathcal F$. After linearization, this fine propagation multiplies the error by $R_f^J(z/J)$, followed by the Parareal set of correction factors $|R_g-R_f^J|/(1-|R_g|)$. Thus the effect of the extra F-relaxation appears precisely as "one additional fine-propagation amplitude." Because the fine propagator is stable, $|R_f^J(z/J)|\le 1$, so this step only contracts and never amplifies—but it does genuinely spend one extra expensive parallel fine solve. In other words, FCF buys "a smaller convergence factor" with "more fine solves," and whether this pays off depends on the unit of comparison.

Moreover, the denominator $1-|R_g(z)|$ of (4.13) reveals a soft spot common to both method classes: when $z\in\sigma(\Delta TA)$ approaches the imaginary axis (weak diffusion, advection-dominated), $|R_g(z)|\to 1$, the denominator tends to zero, and the factor blows up. This is precisely the analytic reason for the sharp deterioration of convergence as viscosity weakens in the subsequent experiments.

![Source Figure 4.8: complex-plane convergence regions of MGRiT and two Parareal iterations at equal fine-solve cost](assets/papers/time-parallelization/source-figures/figure-4-8.svg)

Figure 4.8 uses backward Euler coarse propagation $R_g(z)=1/(1-z)$ and exact fine propagation $R_f(z)=e^z$ to compare one FCF-MGRiT iteration and two Parareal iterations. The top row (a) is the region where $\varrho_{l,\mathrm{MGRiT}}\le\widehat\varrho$, and the bottom row (b) is the region where $\varrho_{l,\mathrm{Parareal}}^2\le\widehat\varrho$ drawn at the same cost of two fine solves; the three columns from left to right correspond to $\widehat\varrho=0.2,0.4,0.6$. The shaded-region contours in each column are quite close, showing that **the fair unit of comparison should be the number of fine propagations**: one FCF-MGRiT iteration (two fine solves) and two Parareal iterations (also two fine solves) have roughly comparable convergence regions in the complex plane. Switching to the two SDIRK fine propagators in (4.9) gives the same conclusion.

### Theorem 4.6: equal-cost constant comparison

Wu and Zhou (2019) give a more quantitative comparison for $z\in\mathbb R^-$. If the fine propagator is L-stable and $J=\Delta T/\Delta t=O(1)$, the worst factor on the negative real axis under backward Euler coarse propagation is approximately

$$
\max_{z\ge0}\varrho_l\approx
\begin{cases}
0.2984,&\text{Parareal},\\
0.1115,&\text{FCF-MGRiT}.
\end{cases}
$$

Switching to second-order Lobatto IIIC (LIIIC-2) coarse propagation gives

$$
\max_{z\ge0}\varrho_l\approx
\begin{cases}
0.0817,&\text{Parareal},\\
0.0197,&\text{FCF-MGRiT}.
\end{cases}
$$

Comparing after balancing by two fine solves: one MGRiT iteration (two fine solves) should be aligned with two Parareal iterations (two fine solves), i.e. compare $\varrho_{l,\mathrm{Parareal}}^2$ against $\varrho_{l,\mathrm{MGRiT}}$. The result is

$$
0.2984^2=0.0890<0.1115,\qquad
0.0817^2=0.0067<0.0197.
$$

Therefore, under these two L-stable coarse propagators, one FCF-MGRiT iteration is actually **slightly slower** than two Parareal iterations.

> [!tip] Insight
> This conclusion is often misread as "MGRiT is worse than Parareal." The real meaning is: **considering only the contraction per fine solve, the extra F-relaxation of FCF does not automatically dominate**. The reason is that the single-round FCF curve is steeper—it reduces the error in one round to the level Parareal reaches only in two rounds, so convergence comes more "in bulk"; but when cost is strictly amortized over each fine solve, the composition of two successive independent Parareal corrections slightly beats one overlapping MGRiT correction. Which to choose depends on whether the practical bottleneck is the number of iterations (large synchronization/communication overhead → favors MGRiT's fewer rounds) or the total number of fine solves (→ the two are comparable, with Parareal's constant slightly better).

### Figures 4.9–4.11: linear and nonlinear experiments

The linear experiments take homogeneous Dirichlet boundaries for the heat equation and periodic boundaries for ADE, sharing

$$
u_0(x)=\sin^2(8\pi(1-x)^2),\quad
T=5,\quad J=20,\quad \Delta T=1/8,\quad \Delta x=1/160,
$$

with backward Euler on the coarse level and SDIRK22 from (4.9) on the fine level.

![Source Figure 4.9: modal-factor distributions on the heat equation and ADE at two viscosities](assets/papers/time-parallelization/source-figures/figure-4-9.svg)

The top row (a) of Figure 4.9 is MGRiT and the bottom row (b) is Parareal; the three columns correspond in turn to the heat equation, ADE with $\nu=0.1$, and ADE with $\nu=0.01$, plotting the distribution of $\varrho_l(J,z)$ over $z\in\sigma(\Delta TA)$. The maximum factors $\varrho_{l,\max}=\max_{z}\varrho_l$ annotated in the figure are $0.08375,0.2718,0.9021$ for MGRiT and $0.2822,0.4453,0.9986$ for Parareal. In the first two columns, the MGRiT factor is roughly equal to the square of the Parareal factor ($0.2822^2\approx0.0796\approx0.08375$, $0.4453^2\approx0.198\lesssim0.2718$), confirming Theorem 4.6's statement that "one MGRiT round ≈ two Parareal rounds with doubled cost"; in the third column both approach $1$, and the coarse propagator can no longer represent the long-lived advection modes under weak diffusion.

![Source Figure 4.10: measured Parareal and MGRiT errors after balancing by two fine solves](assets/papers/time-parallelization/source-figures/figure-4-10.svg)

The four panels of Figure 4.10 are arranged as "heat equation, $\nu=0.1$, $\nu=0.01$, $\nu=0.002$," with the Parareal abscissa counting every two rounds as one (aligning two fine solves), and the dash-dotted line marking the truncation-error level $\max\{\Delta t^2,\Delta x^2\}$ (below which one would not actually iterate). In the first two panels the two curves are close; at $\nu=0.01$ both are very slow, and Parareal degrades more markedly—because it follows **each** fine propagation immediately with an already-distorted coarse propagation (see (4.1)), whereas FCF-MGRiT performs two fine solves in a row with no coarse propagation in between (see (4.12)), thus suffering less from coarse-level distortion. At $\nu=0.002$ both diverge, with maximum modal factors of $1.4211$ for Parareal and $1.2812$ for MGRiT. Thus the per-mode spectral indicators of Figure 4.9 and the measured curves of Figure 4.10 correspond one-to-one across the four diffusion strengths.

The nonlinear Burgers experiment uses homogeneous Dirichlet boundaries, the same initial value, $T=5$, $\Delta T=1/16$, $\Delta x=1/160$, $J=10$, centered differences in space, backward Euler on the coarse level, and SDIRK22 on the fine level. The convergence analysis for the nonlinear case is in Gander et al. (2018b), which proves, under certain Lipschitz conditions on $\mathcal G,\mathcal F$ and their difference: as long as the coarse propagator is sufficiently accurate, the contraction of one FCF-MGRiT round (two fine solves) is comparable to that of two Parareal rounds (two fine solves).

![Source Figure 4.11: comparison at equal fine-solve cost on the Burgers equation at three viscosities](assets/papers/time-parallelization/source-figures/figure-4-11.svg)

The three panels from left to right take $\nu=0.5,0.01,0.002$, again drawing two Parareal rounds as one to align the number of fine solves. In the first two panels, the decrease of one FCF-MGRiT round roughly corresponds to two Parareal rounds; at $\nu=0.002$ both curves first go through a relatively long slow phase before entering a rapid decrease. Lowering the viscosity makes both methods deteriorate in step, which is consistent with the conclusion of the nonlinear Lipschitz analysis and echoes the mechanism in the linear spectral analysis of "weak diffusion, advection-dominated → $|R_g|\to1$ → factor approaching 1."

## Equation, theorem, and figure coverage check

| Source item                                  | Paper section | Coverage status                                                                           |
| -------------------------------------------- | ------------- | ----------------------------------------------------------------------------------------- |
| PFASST origin (Minion 2010, etc.)            | 4.3           | tracing SDC replacing fine propagation and the multilevel-SDC viewpoint                   |
| (4.10)                                       | 4.3           | fine/coarse collocation nodes, copying matrix $\chi$, system $\Phi$                       |
| Lagrange transfer and PFASST block iteration | 4.3           | multigrid reading of $T_{c\to f},T_{f\to c}$ and $B_{01},B_{10},B_{00}$                   |
| (4.11), Figure 4.6                           | 4.3           | implicit-Euler SDC preconditioner, Radau nodes, $Q_f/Q_c$, transfer matrices, PFASST test |
| (4.12), Figure 4.7                           | 4.4           | FCF update, two fine solves, overlap $\Delta T$, finite-step property                     |
| (4.13), Theorem 4.5, Figure 4.8              | 4.4           | MGRiT factor, extra $R_f^J$ contraction factor, equal-cost region                         |
| Theorem 4.6                                  | 4.4           | four worst factors and squared (equal-cost) comparison for two coarse propagators         |
| Figures 4.9–4.11                             | 4.4           | linear spectra, measured error, nonlinear Burgers, weak-diffusion failure                 |

## Source of this page

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 4.3–4.4, pp. 452–460.

### Supplementary references

- P. F. Emmett and M. L. Minion, _Toward an efficient parallel in time method for partial differential equations_, Commun. Appl. Math. Comput. Sci. **7** (2012), 105–132. (introduction of PFASST)
- M. L. Minion, _A hybrid parareal spectral deferred corrections method_, Commun. Appl. Math. Comput. Sci. **5** (2010), 265–301. (precursor of replacing fine propagation with one SDC iteration)
- A. Dutt, L. Greengard and V. Rokhlin, _Spectral deferred correction methods for ordinary differential equations_, BIT **40** (2000), 241–266. (origin of the SDC framework)
- M. Bolten, D. Moser and R. Speck, _A multigrid perspective on the parallel full approximation scheme in space and time_, Numer. Linear Algebra Appl. **24** (2017), e2110; and its convergence analysis (2018). (interpreting PFASST as temporal multigrid)
- R. D. Falgout, S. Friedhoff, Tz. V. Kolev, S. P. MacLachlan and J. B. Schroder, _Parallel time integration with multigrid_, SIAM J. Sci. Comput. **36** (2014), C635–C661. (introduction of MGRiT)
- V. A. Dobrev, Tz. Kolev, N. A. Petersson and J. B. Schroder, _Two-level convergence theory for multigrid reduction in time (MGRiT)_, SIAM J. Sci. Comput. **39** (2017), S501–S527. (the linear convergence estimate of Theorem 4.5)
