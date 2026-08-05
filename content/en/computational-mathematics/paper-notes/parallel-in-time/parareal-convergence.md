---
title: Convergence Analysis for Parareal
description: Papers 12, 20, 30 and 77 - how the contraction factor depends on fine-propagator stability and on the fine grid
lang: en
translation: computational-mathematics/paper-notes/parallel-in-time/parareal-convergence
tags:
  - paper-notes
  - parallel-in-time
  - convergence-analysis
---

## The standard analytical machinery

### The iteration and how to read it

Consider the symmetric positive definite system $\boldsymbol u'(t)+A\boldsymbol u(t)=g(t)$ obtained by semi-discretising a parabolic PDE. With $0=T_0<T_1<\dots<T_{N_t}=T$ and $T_n=n\Delta T$, the parareal iteration is

$$
\boldsymbol u_{n+1}^{k+1}=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^{k})
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k}),
$$

with $J=\Delta T/\Delta t\ge2$ fine steps per coarse interval. The expensive $\mathcal F$ evaluations use only iterate-$k$ data, so they are independent across $n$ and run in parallel; the two $\mathcal G$ terms form a cheap sequential sweep. The iteration can equally be read as a **preconditioned stationary iteration**: writing $M_f$ and $M_g$ for the all-at-once matrices of the fine and coarse solvers, it is exactly

$$
M_g(z)\,\Delta U^k=b-M_f(z)U^k,
\qquad U^{k+1}=U^k+\Delta U^k,
$$

that is, a block lower-bidiagonal $M_g$ preconditioning a block lower-bidiagonal $M_f$. All the parallelism lives in forming the residual $b-M_fU^k$. This reading is the bridge to the other two pages of this topic: replace $M_g$ by a matrix that FFT can diagonalise and you get the [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|diagonalisation]] and [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|all-at-once preconditioning]] routes.

### Error propagation and the two convergence factors

If $A=V_ADV_A^{-1}$ and both propagators are one-step methods with stability functions $R_g$ and $R_f$, the error decouples by mode and obeys $\boldsymbol\xi^{k+1}(z)=M(z)\boldsymbol\xi^k(z)$ with

$$
M(z)=M_g^{-1}(z)\bigl[M_g(z)-M_f(z)\bigr]=I_t-M_g^{-1}(z)M_f(z),
$$

where $M_g(z)$ and $M_f(z)$ are lower-triangular Toeplitz matrices with unit diagonal and subdiagonal entries $-R_g(z)$ and $-R_f^{J}(z/J)$, and $z=\Delta T\lambda$. Note the scaling: the fine solver takes $J$ steps of size $\Delta T/J$ inside a coarse interval, hence $R_f^{J}(z/J)$ rather than $R_f(z)$. Factoring out the scalar gives

$$
M(z)=\bigl[R_f^J(z/J)-R_g(z)\bigr]\,\widetilde M\bigl(R_g(z)\bigr),
$$

with $\widetilde M(\beta)$ strictly lower-triangular Toeplitz carrying $1,\beta,\beta^2,\dots$ down its subdiagonals. **This step is the engine of the whole literature**: it reduces an estimate for a matrix power to a scalar factor times a fixed Toeplitz matrix depending only on $R_g$, whose $\|\widetilde M^k\|_\infty$ is controlled by Gander and Vandewalle's lemma. The linear and superlinear factors follow,

$$
\varrho_l(J,z)=\frac{\bigl|R_g(z)-R_f^{J}(z/J)\bigr|}{1-|R_g(z)|},
\qquad
\varrho_s(J,z,N_t,k)=\frac{\bigl|R_g(z)-R_f^J(z/J)\bigr|^k}{k!}\prod_{j=1}^k(N_t-j),
$$

together with $\max_n\|\boldsymbol e_n^k\|_\infty\le\bigl(\max_{z\in\sigma(\Delta TA)}\varrho_l\bigr)^k\max_n\|\boldsymbol e_n^0\|_\infty$. The factor $\varrho_s$ vanishes at $k=N_t$, which is termination in exact arithmetic; $\varrho_l$ requires $|R_g(z)|<1$ and is the one that governs long time intervals.

The two parts of $\varrho_l$ each have a clear meaning: **the numerator is how badly $\mathcal G$ mismatches $\mathcal F$ over one coarse step, and the denominator is how strongly $\mathcal G$ damps.** The first says how large the correction is, the second says how much that correction is compressed as it accumulates along the sequential coarse sweep.

### The nonlinear case

With $\mathcal F$ exact, $\mathcal G$ of order $p$ with local truncation error at most $C_3\Delta T^{p+1}$, a Lipschitz condition $\|\mathcal G(\cdot,\boldsymbol v)-\mathcal G(\cdot,\boldsymbol w)\|\le(1+C_2\Delta T)\|\boldsymbol v-\boldsymbol w\|$ and an expansion $\mathcal F-\mathcal G=c_{p+1}(\boldsymbol v)\Delta T^{p+1}+\cdots$ with continuously differentiable $c_j$, Gander and Hairer (2008) give

$$
\|\boldsymbol u(T_n)-\boldsymbol u_n^k\|\le
\frac{C_3\Delta T^{p+1}\bigl(C_1\Delta T^{p+1}\bigr)^{k+1}}{(k+1)!}
\,(1+C_2\Delta T)^{n-k-1}\prod_{j=0}^{k}(n-j).
$$

This is a superlinear-type bound structurally identical to $\varrho_s$: the factorial beats the polynomial, so convergence is very fast on short windows, and long windows have to rely on a linear bound of $\varrho_l$ type.

### Where this machinery stops working in practice

The constants the linear theory produces are accurate on dissipative problems and fail on propagative ones, and the failure is continuous rather than abrupt. One directly reported set of runs serves as a yardstick: for the advection-diffusion equation with $T=4$, $\Delta T=0.1$, $\Delta x=1/128$, $J=32$, $\mathcal G$ backward Euler and $\mathcal F$ a two-stage second-order SDIRK method, parareal degrades monotonically as the viscosity $\nu$ decreases and **diverges at roughly $\nu\le10^{-3}$**; for the second-order wave equation it does not converge at all. At the other end of the scale, Gander and Vandewalle's continuous-level analysis (with $\mathcal F$ exact) gives contractions far better than $0.3$, for instance about $0.068$ for Radau IIA.

Those two figures come from the Gander-Wu-Zhou _Acta Numerica_ survey (paper 85) and from Gander-Vandewalle (2007, Table 5.1) rather than from the four papers on this page; they serve only to give the constants below a comparable scale.

**One further structural assumption is built into this machinery, and paper 77 is about removing it.** The action of $\mathcal F$ over a coarse interval is written as the $J$-th **power** of a single scalar stability function at a single argument, $R_f^{J}(z/J)$, and that step requires the $J$ fine steps inside the interval to be equal.

## 12: what happens when the fine propagator is only A-stable

### The idea

Why does parareal converge at all? Because on a **stiff** problem the coarse and fine propagators agree on precisely the dangerous modes. High-frequency modes are killed by both, so their difference is small; low-frequency modes are integrated accurately by both, so their difference is small there too; only a band of intermediate frequencies sees the coarse propagator perform visibly worse, and that band has bounded width. That is where the mesh-independent, $T$-independent constant $\approx0.3$ comes from: it is a peak of $\varrho_l(J,\cdot)$ in the middle of the frequency range, with both ends pushed down.

Once the fine propagator is only A-stable rather than L-stable, **the high-frequency end is no longer pushed down.** The trapezoidal rule has $R_f\to-1$ as $z\to\infty$ and the fourth-order Gauss method has $R_f\to+1$: they do not damp high frequencies, they preserve their amplitude (with an alternating sign in the trapezoidal case). The coarse propagator, backward Euler, still damps them. The difference between the two therefore tends to $1$ at high frequency, so $\varrho_l$ tends to $1$, and parareal **stops contracting on exactly the modes that make the problem stiff**.

The paper's rescue is to accept that a discrete problem has a bounded spectrum: one needs a uniformly small contraction on $z\in[0,z_{\max}]$, not on the whole axis. Given $z_{\max}$, a large enough coarsening ratio $J$ pulls the argument $z/J$ that the fine propagator sees at each fine step back into the region where it still damps, and the $J$-th power restores the decay. The remaining question is only **how large "large enough" is** — which is exactly the paper's critical ratio $J_{\rm cri}$.

### Setting

$A$ is symmetric positive definite, obtained by semi-discretising a parabolic PDE in space, and may contain a fractional Laplacian (in which case $A$ is dense and $\lambda_{\max}$ is very large). The coarse propagator is fixed to backward Euler. Three fine propagators are considered: the trapezoidal rule, a third-order two-stage DIRK method, and the fourth-order two-stage Gauss RK method. With the convention $\boldsymbol u'+A\boldsymbol u=g$ so that $z=\Delta T\lambda$ with $\lambda>0$, the relevant stability functions are

$$
\text{backward Euler:}\ R(z)=\frac{1}{1+z},
\qquad
\text{trapezoidal:}\ R(z)=\frac{2-z}{2+z},
\qquad
\text{fourth-order Gauss:}\ R(z)=\frac{12-6z+z^2}{12+6z+z^2}.
$$

Backward Euler is L-stable with $R(\infty)=0$; the trapezoidal rule is A-stable but not L-stable with $R(\infty)=-1$; the fourth-order Gauss method is A-stable and symplectic with $R(\infty)=+1$. The two-stage third-order DIRK method used in this line of work carries the Butcher tableau

$$
\begin{array}{c|cc}
\gamma & \gamma & 0\\
1-\gamma & -\tfrac{1}{\sqrt3} & \gamma\\ \hline
 & \tfrac12 & \tfrac12
\end{array},
\qquad \gamma=\frac{3+\sqrt3}{6},
$$

whose stability function is a $(2,2)$ rational function, A-stable but not L-stable.

### Derivation: the high-frequency limit shows why the old argument fails

Substitute $z\to\infty$ into $\varrho_l$. With $\mathcal G$ backward Euler, $R_g(z)=1/(1+z)\to0$ so the denominator tends to $1$, while $R_f^J(z/J)\to R_f(\infty)^J$. Hence

$$
\lim_{z\to\infty}\varrho_l(J,z)=\bigl|R_f(\infty)\bigr|^{J}.
$$

Four fine propagators make the situation immediately visible:

| Fine propagator $\mathcal F$          | $R_f(\infty)$ | $\lim_{z\to\infty}\varrho_l(J,z)$ |
| ------------------------------------- | ------------- | --------------------------------- |
| backward Euler (L-stable)             | $0$           | $0$                               |
| exact propagator $e^{-z}$             | $0$           | $0$                               |
| trapezoidal rule (A-stable only)      | $-1$          | $1$                               |
| fourth-order Gauss RK (A-stable only) | $+1$          | $1$                               |

(Computed directly from the expression for $\varrho_l$ and the three stability functions, with $\mathcal G$ backward Euler.)

The top half of the table is why the classical results hold: an L-stable or exact $\mathcal F$ makes the supremum of $\varrho_l$ over the whole real axis a mid-frequency phenomenon, so $\max_{z\in\mathbb R^-}\varrho_l\approx0.3$ for all $J\ge J_{\min}=O(1)$ — and already for all $J\ge2$ when $\mathcal F$ is exact. The bottom half is why the same argument fails outright for the trapezoidal rule and the fourth-order Gauss method: the supremum of $\varrho_l$ is $1$ independently of $J$, so no choice of $J$ improves the whole-axis supremum.

The only way out is therefore to restrict $z$ to the finite range that actually occurs, $[0,z_{\max}]$ with $z_{\max}=\Delta T\lambda_{\max}$. Then $J$ matters again: the larger $J$ is, the closer the evaluation point $z/J$ sits to the region where the fine propagator still damps, and the $J$-th power amplifies that damping. What the paper supplies is the quantitative form of this threshold.

### Theorems

**(A-stable but not L-stable fine propagator.)** With $\mathcal G$ backward Euler and $\mathcal F$ either the trapezoidal rule or the fourth-order Gauss RK method,

$$
\max_{z\in[0,z_{\max}]}\varrho_l(J,z)\approx0.3
\qquad\text{for all }J\ge J_{\min}=O\bigl(\log^2 (z_{\max})\bigr).
$$

The **logarithmic-squared** growth of the threshold in $z_{\max}=\Delta T\lambda_{\max}$ is the precise form of the paper's statement that $J_{\rm cri}$ depends on $\Delta T$, $\Delta t$ and $\lambda_{\max}$, and it is the paper's headline technical content.

**(L-stable contrast.)** If $\mathcal F$ is L-stable, then $\max_{z\in\mathbb R^-}\varrho_l\approx0.3$ for all $J\ge J_{\min}=O(1)$; and if $\mathcal F=\exp(\Delta TA)$ is exact, that rate already holds for all $J\ge2$.

**(The third-order DIRK method is special.)** For the third-order two-stage DIRK fine propagator, $J_{\rm cri}=4$, **independently of $\Delta T$, $\Delta t$ and $\lambda_{\max}$**. It is the only one of the three conclusions with no parameter dependence, and the most usable result in this line: a single fixed coarsening ratio that can be written straight into code. For the trapezoidal rule and the fourth-order Gauss RK method the paper supplies concise closed-form formulas for computing $J_{\rm cri}$.

### Numerical experiments

The numerical examples fall into two classes — problems with fractional PDEs, and uncertainty quantification applications — chosen to support the theoretical predictions. Neither class is arbitrary:

- **Time-dependent PDEs with a fractional Laplacian $(-\Delta)^{\alpha}$.** Semi-discretising in space gives a dense symmetric positive definite $A$ (Toeplitz-structured in one dimension) with a very large $\lambda_{\max}$, hence a very large $z_{\max}=\Delta T\lambda_{\max}$ — exactly the regime where the $J_{\rm cri}$ formulas make a practical difference. If the threshold were only $O(1)$, this class of example would carry no information.
- **Uncertainty quantification.** Many parametric realisations of the same ODE system must each be integrated over a long time, which is where parallel-in-time actually pays. It simultaneously tests the claimed **robustness** of the contraction factor with respect to the spectrum of the coefficient matrix, since different realisations give different $\sigma(A)$ while the bound is supposed to be independent of it.

The contraction constants that circulate around this paper in the literature all sit near $0.3$: Mathew-Sarkis-Schaerer's parareal-Euler gives $\approx0.298$ for $J\ge2$, and the authors' own survey rounds this paper's conclusion to $\approx0.3$ and records two concrete thresholds — $J_{\min}=2$ for the two-stage second-order SDIRK method (Wu, IMA J. Numer. Anal. 2015) and $J_{\min}=4$ for the two-stage third-order SDIRK method (this paper), the latter being exactly the $J_{\rm cri}=4$ above.

### Relation to the others

This is the **foundational** parareal-analysis paper of the series and the base against which the later work is measured. Papers 20 and 30 below carry the same "which $\mathcal F/\mathcal G$ pair contracts robustly" question into fractional problems; paper 77 removes the uniform-fine-grid assumption underlying $R_f^J(z/J)$. Changing the object of analysis instead gives [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|paper 39]], which analyses two-level MGRIT and produces the directly comparable constants $0.2984$ (parareal) and $0.1115$ (MGRIT with FCF-relaxation). The corresponding statements appear as equation (4.8) in Chapter 4 of the [[en/computational-mathematics/knowledge-notes/time-parallelization/index|time-parallelization survey reading]].

## 20: space-fractional problems make the coarse propagator the bottleneck

This paper analyses parareal **on constant time steps** for a time-dependent PDE containing a fractional Laplacian $(-\Delta)^{\alpha}$. Semi-discretising in space gives a dense, symmetric positive definite $A$ with a very large $\lambda_{\max}$ (Toeplitz-structured in one dimension), and that breaks parareal's cost model: a single implicit solve becomes expensive in itself, and the coarse propagator has to perform it sequentially on all $N_t$ coarse intervals of every sweep, directly on the critical path, so even when three or four sweeps suffice the runtime can be consumed by that sequential chain. At the same time the enormous $\lambda_{\max}$ places the problem exactly where paper 12's threshold grows like $O(\log^2(\Delta T\lambda_{\max}))$, and the two difficulties compound: the more expensive coarse propagation is, the more one wants a large $\Delta T$, and a large $\Delta T$ raises $z_{\max}$ and the required $J$. The word "fast" in the title therefore refers to accelerating the parareal **iteration itself** rather than merely bounding its contraction factor.

Two contrasts fix its position. Paper 12 is the direct predecessor: it established the $J_{\rm cri}$ framework for symmetric positive definite $A$ and already used fractional PDEs as a motivating application. Paper 30 is the sibling for the other kind of fractionality — paper 20 treats **space**-fractional problems (nonlocal in space, local in time, so the parareal structure survives) and paper 30 treats **time**-fractional problems (nonlocal in time, destroying independence across subintervals). [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|Paper 31]] faces the same fractional Laplacian but abandons iteration entirely, turning to direct diagonalisation of the time-periodic problem.

## 30: time-fractional problems destroy the load balance

### The idea

The whole premise of parareal is that fine propagation on $[T_n,T_{n+1}]$ needs only the local initial value $\boldsymbol u_n^k$. A time-fractional derivative violates it: the Caputo or Riemann-Liouville derivative at time $t$ depends on the entire solution history on $[0,t]$.

What is worth noticing is that **the failure mode here is not non-convergence but the disappearance of parallelism.** The naive extension can still be written down: when advancing on the $n$-th subinterval, carry the history over $[0,T_n]$ as well. The iteration still converges, but processor $n$ does work proportional to $n$, and the last processor does $N_t$ times the work of the first. Total work grows from $O(N_t)$ to

$$
\sum_{n=1}^{N_t}n=\frac{N_t(N_t+1)}{2},
$$

and since parallel time is set by the slowest processor, the critical path is also $O(N_t)$ fine steps — **the same order as sequential time stepping**. That is what the paper means by unbalanced computational time across processes: the algorithm is parallel on paper and sequential on the machine.

The fix is not to change the iteration but to change **how the equation is written**: compress the nonlocal history into a finite number of auxiliary variables governed by local ODEs. Once that is done every fine-propagator call costs the same on every subinterval, the load balances automatically, and parareal's original cost model is restored.

### Setting

Time-fractional differential equations, whose fractional operator has a historical effect. The goal is a parareal algorithm whose cost is **balanced across processes** and whose convergence rate is robust. The route the paper states is to adopt two recently developed local time-integrators for time-fractional operators, both of which introduce auxiliary variables to localise the fractional operator.

### Derivation

**Step one: localisation.** Let the augmented state be $(\boldsymbol u_n^k,\boldsymbol z_n^k)$ with $\boldsymbol z$ the auxiliary variables, whose evolution is governed by local dynamics. The fractional convolution is replaced by a finite set of extra local unknowns, so fine propagation on the $n$-th subinterval no longer needs the history over $[0,T_n]$, only $(\boldsymbol u_n^k,\boldsymbol z_n^k)$. This turns "work proportional to $n$" back into "work independent of $n$".

**Step two: the parareal update cannot simply be carried over.** Applying the single parareal formula to the whole augmented state, that is

$$
\begin{pmatrix}\boldsymbol u_{n+1}^{k+1}\\ \boldsymbol z_{n+1}^{k+1}\end{pmatrix}
=\mathcal F\!\begin{pmatrix}\boldsymbol u_{n}^{k}\\ \boldsymbol z_{n}^{k}\end{pmatrix}
+\mathcal G\!\begin{pmatrix}\boldsymbol u_{n}^{k+1}\\ \boldsymbol z_{n}^{k+1}\end{pmatrix}
-\mathcal G\!\begin{pmatrix}\boldsymbol u_{n}^{k}\\ \boldsymbol z_{n}^{k}\end{pmatrix},
$$

is precisely what the paper says it does **not** do. The genuinely new algorithmic ingredient is a **mixed coarse-grid correction** in which the auxiliary variables and the solution variable are corrected **separately**, with different rules for the two blocks.

**Why they have to be separated.** The two kinds of variable play different dynamical roles. The auxiliary variables exist in order to approximate a memory kernel, so their modes span a wide range of time scales — that spread is the point of kernel compression. Correcting them with the same coarse propagator as the solution variable means using one $\Delta T$ across time scales separated by orders of magnitude. Separating the corrections lets each class of variable be treated appropriately. (That is a mechanistic reading of the separation, not an argument the paper itself gives.)

### Theorems

The paper's conclusion is that the proposed parareal algorithm **admits a robust rate of convergence**. In this series of papers "robust" consistently means independent of the eigenvalues of the spatial operator, of the coarsening ratio $J$ and of the number of coarse intervals $N_t$.

### Relation to the others

Paper 20 is the space-fractional sibling, and the two technical difficulties do not overlap: the load imbalance caused by the historical effect has no analogue in paper 20. Paper 12 supplies the analytical style — robust, parameter-independent contraction factors — that paper 30 claims to reproduce in the fractional setting. [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|Paper 31]] attacks the same class of problem from the opposite direction: rather than repairing the iteration it uses direct diagonalisation, sidestepping convergence-factor questions entirely. The manoeuvre of localising a nonlocal operator with auxiliary variables also reappears structurally in [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|paper 71]], where the nonlocality comes from coupling a forward and a backward evolution rather than from a fractional kernel, but where the remedy is again to enlarge the state and precondition the enlarged system.

## 77: removing the equal-fine-step assumption

This is a pure analysis paper. It studies parareal convergence with a **uniform coarse grid** and an **arbitrarily distributed, nonuniform fine grid**, leaving the iteration itself unchanged. The phrase "arbitrarily distributed" carries the substance: the fine grid inside a coarse interval is not assumed uniform, not assumed graded, and not even assumed the same from one coarse interval to the next. Adaptive time stepping, graded meshes near $t=0$ for nonsmooth or incompatible initial data, and local refinement around fast transients all produce such grids, so restricting the theory to uniform fine grids excludes essentially all adaptive practice.

What makes this more than a change of notation is that every classical linear convergence result — Gander and Vandewalle's bounds, the constant near $0.3$, the critical-ratio formulas of paper 12, the MGRIT factor of paper 39 — rests on fine propagation over a coarse interval being the $J$-th power $R_f^{J}(z/J)$ of a single scalar stability function at a single argument, which requires the $J$ fine steps inside the interval to be equal. With a nonuniform fine grid one gets instead a **product of distinct factors**,

$$
\mathcal F(T_n,T_{n+1},\cdot)\ \longleftrightarrow\ \prod_{i=1}^{J_n}R_f(\theta_{n,i}z),
\qquad
\theta_{n,i}=\frac{\Delta t_{n,i}}{\Delta T},\quad \sum_i\theta_{n,i}=1 ,
$$

and two things fail. First, $\varrho_l(J,z)$ is no longer even defined, because there is no single $J$. Second, and more seriously, if the fine grid differs across coarse intervals then the subdiagonal entries of $M_f(z)$ vary from row to row, $M_f(z)$ is **no longer Toeplitz**, and the factorisation of the error matrix into a scalar times a fixed Toeplitz matrix — together with Gander and Vandewalle's bound on $\|\widetilde M^k\|_\infty$, which acts on $\widetilde M$ — has nothing left to act on. Recovering a robust, mesh-independent contraction factor in this setting is therefore a genuinely new analysis.

That is exactly the relation to paper 12: paper 12 fixed $\mathcal G$ to backward Euler, allowed three choices of $\mathcal F$, and found the critical ratio $J_{\rm cri}$ **on a uniform fine grid**, with everything depending on the scalar $J$; paper 77 keeps the uniform coarse grid but frees the fine grid entirely, and there is no scalar $J$ left to depend on. [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|Paper 39]] inherits the same uniform-fine-grid restriction through the factor $|R_f^J(z/J)|$ in the MGRIT convergence factor, so paper 77's techniques are the natural route to a nonuniform-grid MGRIT theory; papers 20 and 30 both analyse parareal in constant time steps, and paper 77 is the one that removes that standing assumption from this group of work.

## One tension worth naming

Nonuniform time steps occupy opposite positions in the two routes. In the parareal branch they are an **obstacle to the theory**: they destroy the Toeplitz structure and with it every convergence-factor argument. In the direct diagonalisation branch (see [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|the diagonalisation technique]]) nonuniform steps are a **requirement**: Maday and Rønquist need the $\Delta t_n$ to be distinct for the temporal matrix to be diagonalisable. The same modelling freedom is a nuisance on one route and a prerequisite on the other.

Following the contrast further shows that the two routes have different cost structures. The parareal branch pays in a **contraction factor** determined by the stability functions of $\mathcal F$ and $\mathcal G$, which therefore varies violently with the type of problem — dissipative or propagative. The diagonalisation branch pays in a **condition number** determined by the structure of the temporal matrix, largely indifferent to the physics. That is why parareal collapses around $\nu\approx10^{-3}$ as viscosity decreases while $\alpha$-circulant preconditioning remains usable on hyperbolic problems.

## Sources for this page

- S. Wu and T. Zhou, [_Convergence analysis for three parareal solvers_](https://doi.org/10.1137/140970756), SIAM J. Sci. Comput. 37(2) (2015), pp. A970-A992.
- S. Wu and T. Zhou, [_Fast parareal iterations for fractional diffusion equations_](https://doi.org/10.1016/j.jcp.2016.10.046), J. Comput. Phys. 329 (2017), pp. 210-226.
- S. Wu and T. Zhou, [_Parareal algorithms with local time-integrators for time fractional differential equations_](https://doi.org/10.1016/j.jcp.2017.12.029), J. Comput. Phys. 358 (2018), pp. 135-149.
- S.-L. Wu and T. Zhou, [_Convergence analysis of the parareal algorithm with nonuniform fine time grid_](https://doi.org/10.1137/23M1592481), SIAM J. Numer. Anal. 62(5) (2024), pp. 2308-2330.
