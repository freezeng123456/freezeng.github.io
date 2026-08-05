---
title: The Diagonalisation Technique
description: Papers 31, 39, 46 and 53 - making the temporal matrix diagonalisable and paying for it in conditioning
lang: en
translation: computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique
tags:
  - paper-notes
  - parallel-in-time
  - diagonalisation
---

![Replace the sequential recurrence by a diagonalisable time matrix](assets/diagrams/tao-zhou-papers/en/pint-diagonalization.svg)

## Shared machinery: the three-step pattern and its price

### The idea: why a diagonalisable time matrix turns a sequential sweep into concurrent solves

Time stepping is sequential because the solution at level $n$ depends on level $n-1$. Stacking all time levels into one long vector, that dependence appears as a **block lower-triangular** matrix, and sequential time stepping is precisely block forward substitution on it. Forward substitution cannot be parallelised, because it is a chain by construction.

Diagonalisation changes **the coordinates in which the chain is viewed**. If the temporal matrix factors as $B=VDV^{-1}$, then in the coordinates given by $V^{-1}$ the couplings between temporal degrees of freedom come apart completely: "level $n$ waits for level $n-1$" becomes "**mode** $n$ is solved independently". Each mode carries a complex scalar $\lambda_n$ and requires one complex-shifted spatial problem $(\lambda_nI_x+\cdots)x=y$, and those $N_t$ problems are unrelated, so they can run on $N_t$ processes at once. Multiplying by $V^{-1}$ and by $V$ acts only in the time direction, and if $V$ is a (scaled) DFT matrix those two steps are FFTs costing $O(N_xN_t\log N_t)$, negligible against the spatial solves.

**The whole method therefore succeeds or fails on one sentence: is that temporal matrix diagonalisable, and what is the condition number of $V$?** The second half is the substance. $V^{-1}$ and $V$ are applied in floating point, so the computed solution carries roundoff pollution of order $O(\epsilon\,\mathrm{Cond}_2(V))$ with $\epsilon$ the machine precision. A temporal matrix that is diagonalisable in theory but has $\mathrm{Cond}_2(V)=10^{14}$ is of no use in double precision. The four papers here can be ordered by how they handle that price.

### The all-at-once system and $\alpha$-circulant diagonalisation

Stacking the unknowns at all time levels gives an all-at-once system in Kronecker form,

$$
\bigl(B_1\otimes I_x+B_2\otimes A\bigr)\boldsymbol u=\boldsymbol b,
\qquad
\boldsymbol u=(U_1^\top,\dots,U_{N_t}^\top)^\top .
$$

The organising principle of the whole ParaDiag family fits in one sentence: **every algorithm acts on $B_1$ and $B_2$ and leaves the spatial matrices untouched.** There are two ways to act, and they are the two halves of this literature: diagonalise $B:=B_2^{-1}B_1$ **exactly** (direct solvers), or replace $B_1,B_2$ by $\alpha$-**circulant** matrices to obtain a diagonalisable **preconditioner** (iterative solvers). All four papers here live on the second branch or at its boundary.

The key property of $\alpha$-circulant matrices is that they are **simultaneously** diagonalisable. Let $\mathbb F$ be the unitary DFT matrix with $\omega=e^{2\pi\mathrm i/N_t}$, and for $\alpha\in(0,1]$ set

$$
\Gamma_\alpha:=\mathrm{diag}\bigl(1,\ \alpha^{1/N_t},\ \alpha^{2/N_t},\dots,\ \alpha^{(N_t-1)/N_t}\bigr).
$$

Then **any** two $\alpha$-circulant matrices are diagonalised by the same eigenvector matrix:

$$
C_j^{(\alpha)}=VD_jV^{-1},
\qquad
V=\Gamma_\alpha^{-1}\mathbb F^{*},
\qquad
D_j=\mathrm{diag}\bigl(\sqrt{N_t}\,\mathbb F\,\Gamma_\alpha\,C_j^{(\alpha)}(:,1)\bigr).
$$

Note that $D_j$ needs only the **first column** of $C_j^{(\alpha)}$ — the eigenvalues come from a single length-$N_t$ FFT, with no eigenvalue solver anywhere. The solve becomes three steps,

$$
\text{(a) } S_1=(\mathbb F\Gamma_\alpha\otimes I_x)\boldsymbol b,
\qquad
\text{(b) } S_{2,n}=(\lambda_{1,n}I_x+\lambda_{2,n}A)^{-1}S_{1,n},
\qquad
\text{(c) } \boldsymbol u=(\Gamma_\alpha^{-1}\mathbb F^*\otimes I_x)S_2 ,
$$

where (a) and (c) are scaled FFTs and their inverses at $O(N_xN_t\log N_t)$ and (b) is fully parallel across the $N_t$ time levels. At $\alpha=1$ we have $\Gamma_\alpha=I$ and $V=\mathbb F^*$, recovering the ordinary circulant case.

### The price, in one inequality

$$
\mathrm{Cond}_2(V)=\mathrm{Cond}_2(\Gamma_\alpha^{-1}\mathbb F^{*})
\le\mathrm{Cond}_2(\Gamma_\alpha^{-1})\,\mathrm{Cond}_2(\mathbb F^{*})
=\mathrm{Cond}_2(\Gamma_\alpha^{-1})\le\frac1\alpha .
$$

The unitary factor $\mathbb F^*$ is free; all the conditioning comes from the scaling $\Gamma_\alpha^{-1}$, whose diagonal ranges from $1$ to $\alpha^{-(N_t-1)/N_t}$, giving the bound $1/\alpha$. **This single inequality is the master trade-off of the route**: smaller $\alpha$ brings the $\alpha$-circulant matrix closer to the original lower-triangular Toeplitz matrix (a better approximation and a faster outer contraction) and worsens the conditioning of $V$ (more roundoff amplification).

Writing both sides together shows how $\alpha$ should be chosen. On the contraction side, use the standard bound of this route, $\rho\le\alpha/(1-\alpha)$ (established for the trapezoidal rule by Gander and Wu 2019, and for all stable one-step methods by [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|paper 65]]); on the roundoff side use $\epsilon\,\mathrm{Cond}_2(V)\le\epsilon/\alpha$ with $\epsilon=2.22\times10^{-16}$ in double precision:

| $\alpha$   | $\rho\le\alpha/(1-\alpha)$ | $\mathrm{Cond}_2(V)\le1/\alpha$ | accuracy floor $\approx\epsilon/\alpha$ | iterations to reach $10^{-10}$ |
| ---------- | -------------------------- | ------------------------------- | --------------------------------------- | ------------------------------ |
| $10^{-1}$  | $0.111$                    | $10$                            | $2.2\times10^{-15}$                     | $11$                           |
| $10^{-2}$  | $0.0101$                   | $10^{2}$                        | $2.2\times10^{-14}$                     | $6$                            |
| $10^{-3}$  | $0.00100$                  | $10^{3}$                        | $2.2\times10^{-13}$                     | $4$                            |
| $10^{-13}$ | $10^{-13}$                 | $10^{13}$                       | $2.2\times10^{-3}$                      | $1$                            |

(The first two columns are the two bounds above. The accuracy floor and the iteration counts are computed directly from them, the latter as $\lceil 10/\log_{10}(1/\rho)\rceil$.)

The last row shows why $\alpha$ cannot simply be made small: $\alpha=10^{-13}$ really does converge in one sweep, but only to an accuracy of about $10^{-3}$. The ParaDiag survey's practical recommendation is exactly $\alpha=10^{-2}$ and $10^{-3}$, with an explicit warning against $\alpha=10^{-13}$; the table puts the reason for that advice in one place.

### Three ways of writing the same method

This route appears in three languages across the literature, and recognising them as the same object saves a great deal of confusion.

1. **A preconditioned stationary iteration.** $\mathcal P_\alpha\Delta\boldsymbol U^{k}=\boldsymbol b-\mathcal K\boldsymbol U^{k}$ with $\boldsymbol U^{k+1}=\boldsymbol U^{k}+\Delta\boldsymbol U^{k}$. The case $\alpha=1$ is the discrete construction of McDonald, Pestana and Wathen (2018); $\alpha\in(0,1)$ is the parallel method of Banjai and Peterseim (2012).
2. **Krylov acceleration.** Apply GMRES, MINRES, BiCGStab or CG directly to $\mathcal P_\alpha^{-1}\mathcal K\boldsymbol U=\mathcal P_\alpha^{-1}\boldsymbol b$. This is the stationary iteration written at its fixed point, and it **can work even when $\rho(\mathcal P_\alpha^{-1}\mathcal K)\ge1$** — which matters for nonlinear problems.
3. **Head-tail coupled waveform relaxation** (Gander and Wu, _Numer. Math._ 143 (2019) 489-527), at the continuous level:

   $$
   \boldsymbol u^k_t(t)=A\boldsymbol u^k(t)+\boldsymbol g(t),
   \qquad
   \boldsymbol u^k(0)=\alpha\bigl[\boldsymbol u^k(T)-\boldsymbol u^{k-1}(T)\bigr]+\boldsymbol u_0 .
   $$

   At convergence the tail term cancels and the original initial-value problem is recovered. Discretising this **periodic-like** problem gives exactly $\mathcal P_\alpha\boldsymbol U^k=\boldsymbol b^k$.

The third form carries the conceptual content of the whole route: **imposing an $\alpha$-periodic boundary condition in time is what makes the temporal matrix circulant and hence FFT-diagonalisable.** That one idea runs through papers 31, 39, 46 and 53 here, and 65, 71 and 84 on the next page.

## 31: the time-periodic case makes the price zero

### The idea

The table above describes a transaction that has to be paid for. Paper 31's observation is that **for one class of problem the payment is exactly zero.**

Time-periodic problems are bad news for time stepping: there is no initial condition to march from, and the periodicity condition $u(0)=u(T)$ couples the last time level back to the first, so all discrete solutions must be considered at once. Whether or not parallelism is the goal, one is forced into an all-at-once formulation. From the diagonalisation point of view this is good news: the periodic coupling makes the time-discretisation matrix a genuine **circulant**, and circulants are diagonalised by the **unitary** DFT matrix. Then $\mathrm{Cond}_2(V)=1$, the parameter $\alpha$ never has to appear, and there is no approximation for an outer iteration to repair. The result is a **direct** (non-iterative) parallel-in-time algorithm.

The second half of the paper deals with the problem this creates. After decoupling, every temporal mode leaves behind a **complex-coefficient** spatial problem whose operator is a dense fractional Laplacian. Direct factorisation is out of reach, so an iterative solver is needed, and the complex shift means standard multigrid theory does not carry over unchanged.

### Setting

Time-periodic diffusion equations with a fractional Laplacian. Semi-discretisation in space gives a dense $A$; time discretisation uses the $\theta$-method with a periodicity condition. The all-at-once system is the Kronecker form above, with $B_1$ and $B_2$ genuine circulants — in the language of the $\alpha$-circulant family, the special value $\alpha=1$:

$$
C_1^{(\alpha)}=\frac{1}{\Delta t}
\begin{bmatrix}1&&&-\alpha\\-1&1&&\\&\ddots&\ddots&\\&&-1&1\end{bmatrix},
\qquad
C_2^{(\alpha)}=
\begin{bmatrix}\theta&&&(1-\theta)\alpha\\1-\theta&\theta&&\\&\ddots&\ddots&\\&&1-\theta&\theta\end{bmatrix}.
$$

### Derivation

**Layer one: direct diagonalisation.** At $\alpha=1$ we have $\Gamma_\alpha=I$ and $V=\mathbb F^*$ unitary, so $\mathrm{Cond}_2(V)=1$. The shared three-step pattern applies as usual but with three differences: no $\alpha$ has to be chosen, no outer Krylov iteration is needed, and there is no roundoff amplification. **The entire trade-off that constrains the other diagonalisation routes disappears here.** That is why the time-periodic case admits a direct parallel-in-time algorithm.

**Layer two: multigrid for the complex-shifted systems.** Step (b) leaves a series of independent complex-coefficient systems $(\lambda_{1,n}I_x+\lambda_{2,n}A)x=y$. Two features determine the choice of solver. First, $A$ is a dense fractional Laplacian, so any smoother needing triangular solves or sparsity (damped Jacobi, Gauss-Seidel) is unsuitable, while **damped Richardson iteration**, which needs only matrix-vector products, fits exactly. Second, the shift is complex, so the multigrid convergence theory for the real symmetric positive definite case cannot simply be invoked and the convergence factor has to be proved afresh. The paper proves that the multigrid solver has a **mesh-independent** convergence factor.

**Layer three: optimising the damping parameter.** The convergence factor of damped Richardson is a function of the damping parameter, and the paper optimises that parameter to **minimise this constant convergence factor**.

### Theorems

The paper establishes three things: the diagonalisation yields a **direct parallel-in-time computation of all discrete solutions**, non-iterative in time; the multigrid solver for the complex-shifted spatial systems has a **mesh-independent convergence factor**; and there is a choice of damping parameter minimising that constant factor.

The three sit in a hierarchy worth naming. The first layer needs no numerical support, being exact and non-iterative; the second and third are the paper's analysis, and "mesh-independent" by itself only says the convergence factor does not degrade under refinement, not that it is small.

### Relation to the others

This is the **pivot** of the series: the first paper to use diagonalisation rather than the parareal iteration, and everything afterwards builds on it. Paper 46 is the decisive generalisation — it diagonalises $B_{\rm per}$ for the time-periodic optimal control problem by the same mechanism and then takes the crucial extra step of using a periodic-like approximation $\widehat B_{\rm per}$ as a preconditioner for the **initial-value** problem. Papers 39 and 53 here, and 65, 71 and 84 on the next page, all use $\alpha\in(0,1)$, the branch where the price is $\mathrm{Cond}_2(V)\le1/\alpha$ and the reward is an $O(\alpha)$ contraction. [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|Paper 59]] attacks the **other** diagonalisation obstruction — the ill-conditioning of $V$ in the Maday-Rønquist geometric-step route — by choosing a time integrator whose $V$ satisfies $\mathrm{Cond}_2(V)=O(n^2)$.

## 39: fixing the coarse propagator and the coarse correction together

### The idea

MGRIT is the multilevel generalisation of parareal and the algorithm behind the XBraid package. It performs one extra F-relaxation compared with parareal, and its convergence factor carries one extra factor as a result. But there is an accounting point that is easy to miss: **that extra F-relaxation also costs a fine solve**, so "MGRIT contracts faster than parareal" is not by itself a gain.

The paper starts from two concrete inefficiencies that pull against each other. First, the coarse propagator is almost always backward Euler, which is only first order; a coarse propagator with better damping ought to improve the contraction substantially, but nobody had asked. Second, the coarse-grid correction is a **sequential** sweep through $N_t$ coarse levels, so by Amdahl's law it becomes the serial bottleneck as processor counts grow. **The two are in conflict**: a better coarse propagator is usually a multi-stage implicit Runge-Kutta method costing several backward-Euler steps, so fixing the first worsens the second.

The paper fixes both and shows they are compatible: switch the coarse propagator to second-order Lobatto IIIC, and use a head-tail coupling to turn the sequential coarse correction into an all-at-once system that FFT can diagonalise. The extra cost of the first change is absorbed by the parallelism that the second one creates. That is why each change alone delivers less than the two together.

### Setting

A linear ODE system with a symmetric positive definite coefficient matrix. Two-level MGRIT with FCF-relaxation reads

$$
\begin{aligned}
\boldsymbol u_0^{k+1}&=\boldsymbol u_0,\qquad \boldsymbol u_1^{k+1}=\mathcal F(T_0,T_1,\boldsymbol u_0),\\
\boldsymbol u_{n+1}^{k+1}&=\mathcal F\bigl(T_n,T_{n+1},\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^{k})\bigr)
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G\bigl(T_n,T_{n+1},\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^{k})\bigr),
\end{aligned}
$$

for $n=1,\dots,N_t-1$. On terminology: **F-relaxation** applies the fine propagator on the F-points inside each coarse interval, starting from the C-point value at its left end, and those solves are independent and parallel; **C-relaxation** updates the C-point values from the last F-point; **FCF-relaxation** is F, then C, then F again, which is why the iteration costs two fine solves. Comparing with the parareal iteration shows the exact relationship: **MGRIT-FCF is parareal with every $\boldsymbol u_n^k$ replaced by $\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)$**, that is, parareal with one coarse interval of overlap. Two-level MGRIT with F-relaxation only **is** parareal.

The corresponding linear convergence factor is

$$
\varrho_l(J,z)=\frac{\bigl|R_f^{J}(z/J)\bigr|\,\bigl|R_g(z)-R_f^{J}(z/J)\bigr|}{1-|R_g(z)|}
\qquad\Longleftrightarrow\qquad
\varrho_{l,\mathrm{MGRIT}}=\bigl|R_f^J(z/J)\bigr|\times\varrho_{l,\mathrm{parareal}} .
$$

The extra $|R_f^J(z/J)|\le1$ is the entire gain from the additional F-relaxation: **the extra fine solve buys exactly one more application of the fine stability function.**

### Derivation

**Strategy one: change the coarse propagator.** Take the second-order Lobatto IIIC method as $\mathcal G$, with Butcher tableau and stability function

$$
\begin{array}{c|cc}
0 & \tfrac12 & -\tfrac12\\
1 & \tfrac12 & \tfrac12\\ \hline
 & \tfrac12 & \tfrac12
\end{array},
\qquad
R_g(z)=\frac{1}{1+z+z^2/2}
$$

(the latter written in the convention $\boldsymbol u'+A\boldsymbol u=g$). This is the $(0,2)$ Padé approximant: A-stable, **L-stable** ($R_g(\infty)=0$) and second order. What matters is not the order but the **rate of decay**: compared with backward Euler's $1/(1+z)$, the function $1/(1+z+z^2/2)$ decays like $z^{-2}$ rather than $z^{-1}$ at high frequency. Reading $\varrho_l$ as before — numerator is the coarse-fine mismatch, denominator is the dissipation margin — stronger high-frequency damping simultaneously shrinks the numerator ($R_g$ is closer to the equally decaying $R_f^J$) and raises the denominator ($1-|R_g|$ is closer to $1$). The two effects push the same way, and that is the mechanism behind the drop from $0.2984$ to $0.0817$.

**Strategy two: parallelise the coarse correction.** Replace the sequential coarse correction by one with a **head-tail coupling condition**:

$$
\boldsymbol u_{n+1}^{k+1}=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\underbrace{\mathcal F(T_n,T_{n+1},\tilde{\boldsymbol u}_n^{k})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k})}_{=:\,\boldsymbol b_{n+1}^k},
\qquad
\boldsymbol u_0^{k+1}=\alpha\,\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0 ,
$$

with $\tilde{\boldsymbol u}_n^k=\boldsymbol u_n^k$ for $n\ge1$ and $\tilde{\boldsymbol u}_0^k=\boldsymbol u_0$; that redefinition is necessary or the fixed point is no longer the true solution. For backward Euler as $\mathcal G$ and $\boldsymbol u'=A\boldsymbol u$, the coupled coarse correction is the all-at-once system

$$
\bigl(C_\alpha\otimes I_x-I_t\otimes\Delta T A\bigr)\boldsymbol U^{k+1}=\boldsymbol g^k ,
\qquad
C_\alpha=\begin{bmatrix}1&&&-\alpha\\-1&1&&\\&\ddots&\ddots&\\&&-1&1\end{bmatrix},
$$

that is, an $\alpha$-circulant matrix replacing the original lower-bidiagonal Toeplitz one, so the shared three-step pattern applies: FFT, $N_t$ parallel complex-shifted spatial solves, inverse FFT. Note that $\alpha\to0$ recovers the standard sequential coarse correction, and $\mathrm{Cond}_2(V)\le1/\alpha$ is the roundoff price. This is an instance of the third formulation above: $\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0$ is literally "make the problem $\alpha$-periodic".

**Strategy two, addendum: make a two-stage method cost as much as a one-stage one.** Within the parallel coarse-correction framework the cost of the implicit two-stage Runge-Kutta method LIIIC-2 can be reduced to **the cost of backward Euler**, again by a suitable application of the diagonalisation technique. This is what makes the two strategies genuinely compatible: without it, the cost of the better coarse propagator would cancel the gain from parallelising the coarse correction.

### Theorems

**(Robustness.)** The paper proves a **robust** convergence factor for MGRIT, **independent of the eigenvalues of the coefficient matrix and of the ratio $J=\Delta T/\Delta t$**.

**(Four directly comparable constants.)** With $\mathcal F$ an L-stable integrator and $J=O(1)$, the values of $\max_{z\ge0}\varrho_l$ are

| Coarse propagator $\mathcal G$ | Parareal | MGRIT (FCF-relaxation) |
| ------------------------------ | -------- | ---------------------- |
| backward Euler                 | 0.2984   | 0.1115                 |
| Lobatto IIIC (second order)    | 0.0817   | 0.0197                 |

This is the precise form of the paper's "from 0.1 to 0.02".

**(How to choose $\alpha$.)** For a suitable choice of $\alpha$, the new algorithm with parallel coarse correction has the **same convergence rate** as the original. For the parareal case there is an explicit threshold, due to Wu (SISC 2018): provided

$$
\alpha\le\frac{\rho}{1+\rho},
$$

one has $\rho_{\text{new}}=\rho$. Since $\rho=O(10^{-1})$ in practice, $\alpha=O(10^{-1})$ suffices, and the diagonalisation roundoff amplification $\mathrm{Cond}_2(V)\le1/\alpha$ is then negligible. That explicit threshold belongs to the parareal case; for MGRIT paper 39 states only that a suitable $\alpha$ preserves the convergence rate.

> [!warning] Equal-cost accounting: one MGRIT-FCF iteration is slightly worse than two parareal iterations
> The two columns above **cannot be compared across the row**, because their per-iteration costs differ. FCF-relaxation performs two fine solves per sweep, so the right comparison for MGRIT is against two parareal sweeps, and
>
> $$
> 0.2984^2=0.0890<0.1115,
> \qquad
> 0.0817^2=0.0067<0.0197 .
> $$
>
> Two parareal sweeps beat one MGRIT-FCF sweep in both cases. **The robust gain comes from changing the coarse propagator, not from changing the relaxation.**

Converting the four constants into an equivalent contraction per fine solve (the last column is $\varrho_l$ taken to the power one over the number of fine solves per sweep) makes the point plainer:

| Coarse propagator | Method    | $\varrho_l$ | Fine solves per sweep | Equivalent contraction per fine solve |
| ----------------- | --------- | ----------- | --------------------- | ------------------------------------- |
| backward Euler    | parareal  | 0.2984      | 1                     | 0.2984                                |
| backward Euler    | MGRIT-FCF | 0.1115      | 2                     | 0.3339                                |
| LIIIC-2           | parareal  | 0.0817      | 1                     | 0.0817                                |
| LIIIC-2           | MGRIT-FCF | 0.0197      | 2                     | 0.1404                                |

Both parareal rows beat the corresponding MGRIT rows. Within either column, by contrast, moving from backward Euler to Lobatto IIIC lowers the contraction factor by roughly a factor of $3.7$ for parareal and $5.7$ for MGRIT-FCF — and within the parallel coarse-correction framework the cost of that better coarse propagator is pushed back down to backward-Euler level, which is where the paper's net gain comes from.

### Numerical experiments

The two test problems each have a clear purpose:

- **Advection-diffusion equations with uncertain coefficients.** This is a parametric family of ODE systems (a UQ setting) in which every realisation must be integrated over a long time, exactly where parallel-in-time pays. It also tests the robustness claim above, since different realisations give different $\sigma(A)$ while the bound is supposed not to depend on it.
- **The Gray-Scott model.** A stiff nonlinear reaction-diffusion system with pattern formation. It pushes the algorithm outside the reach of the linear theory and tests whether the four constants retain predictive value on a nonlinear problem.

The measured factors match $\approx0.1$ (backward Euler coarse propagator) and $\approx0.02$ (LIIIC-2 coarse propagator), and the parallel coarse correction does not degrade the rate.

Where these constants stop applying is fixed by a different set of numbers. The viscosity sweep on the heat and advection-diffusion equations in the Gander-Wu-Zhou _Acta Numerica_ survey (paper 85) gives:

| Viscosity $\nu$ | $\varrho_{l,\max}$ (parareal) | $\varrho_{l,\max}$ (MGRIT-FCF)      | State                                      |
| --------------- | ----------------------------- | ----------------------------------- | ------------------------------------------ |
| $0.1$           | —                             | $\approx\varrho_{l,\rm parareal}^2$ | the squaring relation holds                |
| $0.01$          | close to $1$                  | close to $1$                        | both on the verge of failing               |
| $0.002$         | $1.4211$                      | $1.2812$                            | both diverge, MGRIT degrades less severely |

The survey's explanation for the last row is that MGRIT performs two consecutive fine solves without an intervening coarse solve that would by then be unhelpful. These numbers mark the boundary of the relation $\varrho_{l,\rm MGRIT}\approx\varrho_{l,\rm parareal}^2$: once the viscosity is low enough that both are near divergence, the squaring relation goes the way of the linear theory behind it.

### Relation to the others

[[en/computational-mathematics/paper-notes/parallel-in-time/parareal-convergence|Paper 12]] supplies the analytical apparatus ($\varrho_l$, stability functions, robustness in $z$ and $J$) and the reference constant $0.2984$ for parareal; paper 39 is the MGRIT counterpart and the one that pins down all four constants side by side. Paper 31 supplied the diagonalisation idea, and paper 39 is the first in the list to use it **inside an iterative method**, as a device for parallelising the coarse correction rather than as a standalone direct solver — the $C_\alpha$ here is the same object that becomes the preconditioner $\mathcal P_\alpha$ in paper 53 and in [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|papers 65, 71 and 84]]. Paper 46 makes the same "replace the sequential structure by a circulant one and diagonalise" move in the optimal-control setting, where the sequential structure is a forward-backward pair rather than a single forward sweep. Paper 65 is the eventual generalisation, replacing case-by-case analyses with a uniform spectral analysis of $\mathcal P_\alpha^{-1}\mathcal K$; [[en/computational-mathematics/paper-notes/parallel-in-time/parareal-convergence|paper 77]] relaxes the uniform-fine-grid hypothesis that makes $R_f^J(z/J)$ meaningful.

## 46: a forward-backward system has no single direction of propagation

### The idea

Everything so far treats evolution in **one direction**: information flows from $t=0$ to $t=T$, the all-at-once matrix is block lower triangular, and every parallel-in-time method is a way of avoiding forward substitution on it. Optimal control removes that premise outright: the state equation runs forward from an initial condition, the adjoint runs backward from a terminal condition, and the two are coupled pointwise in time.

**The key difference is that parallel-in-time here is no longer an accelerator but a necessity.** Sequential time stepping does not apply to the coupled system at all — you cannot march forward and backward simultaneously — and parareal-style remedies do not apply either, because there is no single direction of propagation to iterate along. The question is not how to make it faster but how to solve it.

Paper 46's observation builds on paper 31 and then takes a decisive extra step. Paper 31 says that a time-periodic problem has a circulant time matrix and is therefore diagonalisable for free. Paper 46 first applies that to the time-periodic optimal control problem, obtaining a direct algorithm, and then asks: **if periodic problems are this convenient, why not approximate a non-periodic problem by a periodic one?** The time matrix of the initial-value problem is a non-diagonalisable Jordan block, but as long as it appears inside the **preconditioner** rather than the system itself, it can be replaced by a nearby diagonalisable periodic-type matrix and the error handed to an outer Krylov iteration. That step is the conceptual heart of the entire ParaDiag-II family.

### Setting

Distributed control of a parabolic PDE: minimise

$$
\mathcal J(y,q)=\tfrac12\|y-y_d\|^2_{L^2(\Omega\times(0,T))}
+\tfrac{\gamma}{2}\|q\|^2_{L^2(\Omega\times(0,T))}
\quad\text{s.t.}\quad
\partial_t y-\Delta y=f+q,
$$

with $y=0$ on $\partial\Omega\times(0,T)$ and one of two end conditions: $y(\cdot,0)=y(\cdot,T)$ (the **time-periodic** model problem) or $y(\cdot,0)=y_0$ (the **initial-value** model problem). Eliminating the control via $\gamma q=p$, the first-order optimality (KKT) system is the forward-backward pair

$$
\begin{cases}
\partial_t y-\Delta y-\tfrac1\gamma p=f, & y(\cdot,0)=y_0\ \text{or}\ y(\cdot,0)=y(\cdot,T),\\[2pt]
-\partial_t p-\Delta p+y=y_d, & p(\cdot,T)=0\ \text{or}\ p(\cdot,0)=p(\cdot,T).
\end{cases}
$$

After space-time discretisation, with $A$ the discrete negative Laplacian and $B\in\{B_{\rm per},B_{\rm ini}\}$ the time-discretisation matrix, the all-at-once KKT system takes the form

$$
\begin{bmatrix}
B\otimes I_x+I_t\otimes A & -\tfrac{1}{\gamma}I_t\otimes I_x\\[3pt]
I_t\otimes I_x & B^\top\otimes I_x+I_t\otimes A
\end{bmatrix}
\begin{bmatrix}\boldsymbol y\\ \boldsymbol p\end{bmatrix}
=\begin{bmatrix}\boldsymbol f\\ \boldsymbol y_d\end{bmatrix}.
$$

The **transpose** $B^\top$ in the $(2,2)$ block is the algebraic signature of the backward adjoint: $B$ is lower triangular (causal, forward), $B^\top$ upper triangular (anti-causal, backward).

The paper's main idea sits in one place: **careful handling of the two time discretisation matrices** $B_{\rm per}$ and $B_{\rm ini}$.

### Derivation

**The decisive structural fact.** For a uniform step size and backward Euler, the two time matrices are

$$
B_{\rm ini}=\frac{1}{\Delta t}\begin{bmatrix}1&&&\\-1&1&&\\&\ddots&\ddots&\\&&-1&1\end{bmatrix},
\qquad
B_{\rm per}=\frac{1}{\Delta t}\begin{bmatrix}1&&&-1\\-1&1&&\\&\ddots&\ddots&\\&&-1&1\end{bmatrix}.
$$

$B_{\rm ini}$ is lower-bidiagonal Toeplitz, its only eigenvalue is $1/\Delta t$ with algebraic multiplicity $N_t$, and $B_{\rm ini}-\tfrac1{\Delta t}I$ is strictly lower bidiagonal of rank $N_t-1$, so the geometric multiplicity is $1$. **It is a single Jordan block: defective and not diagonalisable.** By contrast $B_{\rm per}$ is circulant, diagonalised by the unitary DFT matrix with $\mathrm{Cond}_2(V)=1$. This dichotomy is exactly what the paper exploits.

**Algorithm one: a direct parallel-in-time method for the time-periodic problem.** Diagonalise $B_{\rm per}=VDV^{-1}$ with $V$ the DFT matrix. The KKT system then decouples in time into $N_t$ independent complex systems of size $2N_x\times2N_x$, one per temporal Fourier mode, each coupling the state and the adjoint at that mode, all solvable in parallel. **No iteration, no convergence factor, no $\alpha$ parameter and no roundoff penalty.**

**Algorithm two: precondition the initial-value problem by a periodic surrogate.** Since $B_{\rm ini}$ is defective and cannot be diagonalised, replace it **inside the preconditioner** by a nearby periodic-type matrix $\widehat B_{\rm per}$ that FFT can diagonalise, so that algorithm one's machinery applies the preconditioner. The error introduced by the substitution is absorbed by an outer Krylov method. The two matrices differ in a single corner entry, which is why the substitution can be expected to make a good preconditioner.

### Theorems

**(Directness.)** For the time-periodic problem the algorithm is **direct** (non-iterative) and parallel across all time steps.

**(Clustering.)** For the initial-value problem, and **for both the backward-Euler method and the trapezoidal rule**, the clustering of the **eigenvalues and the singular values** of the preconditioned matrix is justified. Two things about this pair of statements deserve emphasis:

- Proving **singular-value** clustering as well as eigenvalue clustering is both stronger and more useful. The preconditioned matrix is generally non-normal, and for non-normal matrices eigenvalue clustering **does not** control GMRES convergence; it is the singular-value information that yields genuine residual bounds.
- Covering the **trapezoidal rule** (A-stable but not L-stable) rather than backward Euler alone mirrors the L-stable versus A-stable-only distinction that drove [[en/computational-mathematics/paper-notes/parallel-in-time/parareal-convergence|paper 12]]. The same dividing line appears on the iterative route as a threshold on the contraction factor and on the preconditioning route as whether clustering holds at all.

**(Comparison with Schur-complement preconditioners.)** Compared with existing preconditioners designed by **approximating the Schur complement** of the discrete KKT system, the new preconditioner gives **much faster convergence** for certain Krylov subspace solvers, namely **GMRES and BiCGStab**.

### Numerical experiments

Two model problems (time-periodic and initial-value parabolic control), a baseline of Schur-complement-based preconditioners, and GMRES and BiCGStab as outer solvers.

The choice of outer solver is itself informative: using GMRES and BiCGStab rather than CG means the preconditioned system being iterated on is **not symmetric positive definite**. That contrasts sharply with [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|paper 71]], which uses CG and therefore has a symmetric positive definite preconditioned system — consistent with that paper being identified as a parallel version of the matching Schur complement preconditioner. **Two papers on the same class of physical problem choosing different outer solvers tells us their preconditioned objects differ algebraically.**

### Relation to the others

Paper 31 is the direct ancestor: it observed that time-periodicity makes the time matrix circulant and hence gives a free, perfectly conditioned diagonalisation. Paper 46 carries that into the optimal-control setting and takes the extra step of using the periodic structure as an **approximation** for the non-periodic problem. Paper 39 makes the same move one level down — its head-tail coupling $\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0$ is literally "make the problem $\alpha$-periodic so the time matrix becomes $\alpha$-circulant" — while paper 46 applies it to a coupled forward-backward pair rather than a single forward sweep. [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|Paper 71]] is the mature successor: it treats general forward-backward evolutionary equations, gives the preconditioner explicitly as "replace the Toeplitz matrices by $\alpha$-circulant matrices", and proves a mesh-independent spectral interval for **any** one-step stable integrator. Where paper 46 proves clustering for two specific integrators, paper 71 proves an interval for all of them — exactly as paper 65 does for paper 53.

## 53: multi-stage methods have no ready-made all-at-once form

### The idea

The $\alpha$-circulant recipe looks like a one-liner: replace the Toeplitz blocks of the all-at-once matrix by $\alpha$-circulant blocks. But the recipe has an implicit prerequisite — **there has to be an all-at-once matrix in the first place.**

For linear multistep methods that prerequisite is automatic: the difference equation at each step involves the solution at a few time levels, so stacking them gives a block Toeplitz matrix. **Multi-stage Runge-Kutta methods do not satisfy it.** A Runge-Kutta step couples **stage values**, and stage values are not time-level unknowns; they are intermediate quantities within a step, normally eliminated inside it. The paper puts it directly: for multi-stage integrators "we can not directly formulate the difference equation into an all-at-once system". There is therefore no Toeplitz matrix to circulantise and the recipe has nothing to act on.

What the paper does is **rebuild** the all-at-once form — keeping the stage values as unknowns, reorganising them into a new block structure, and constructing the $\alpha$-circulant preconditioner on that new structure. It is the first paper in this list to do so for a multi-stage method. Its conclusion has an elegant shape: the algebraic hypothesis needed for a robust convergence rate turns out to be exactly the method's classical A-stability condition.

### Setting

Two-stage singly diagonally implicit Runge-Kutta (SDIRK) methods. The two Butcher tableaux used in this line of work are

$$
\underbrace{\begin{array}{c|cc}\gamma&\gamma&0\\ 1-\gamma&1-\gamma&\gamma\\ \hline &1-\gamma&1-\gamma\end{array}}_{\text{SDIRK22},\ \gamma=\frac{2-\sqrt2}{2}}
\qquad
\underbrace{\begin{array}{c|cc}\gamma&\gamma&0\\ 1-\gamma&\frac{-1}{\sqrt3}&\gamma\\ \hline &\frac12&\frac12\end{array}}_{\text{SDIRK23},\ \gamma=\frac{3+\sqrt3}{6}}
$$

"Singly diagonally implicit" means both diagonal entries equal the same $\gamma$, so both stage solves use the **same** shifted matrix $I+\gamma\Delta tA$ and a single factorisation serves both. That $\gamma$ is exactly what the paper calls the principle element. Both tableaux come from the authors' own survey.

### Derivation

The construction fits in one sentence: the preconditioner $\mathcal P_\alpha$ built here is **also a block $\alpha$-circulant matrix, but with completely different structures and different implementation details from the linear-multistep case**. It is $\alpha$-circulant at the block level, $\alpha\in(0,1)$, and applied via block-Fourier diagonalisation so the time steps decouple. For contrast, the linear-multistep machinery being extended is set out in full in the shared section at the top of this page — simultaneous diagonalisation, $D_j$ from the first column alone, $\mathrm{Cond}_2(V)\le1/\alpha$, and the reference contraction $\rho\le\alpha/(1-\alpha)$.

### Theorems

**(Main theorem.)** For $\alpha\in(0,1)$ the spectral radius of the iteration matrix satisfies

$$
\rho(\text{iteration matrix})\le\frac{\alpha}{1-\alpha},
\qquad\text{provided}\qquad \gamma\ge\tfrac14 ,
$$

where $\gamma$ is the principal diagonal element of the two-stage SDIRK method. This is exactly the bound of the trapezoidal-rule case in the linear-multistep theory, so the multi-stage extension loses nothing in rate.

**(The punchline.)** The condition $\gamma\ge1/4$ **coincides exactly with the A-stable condition** of the SDIRK method. The algebraic hypothesis needed for a robust parallel-in-time convergence rate is therefore not an extra artificial restriction but a classical stability condition. Both tableaux above satisfy it: $(2-\sqrt2)/2\approx0.2929\ge0.25$ and $(3+\sqrt3)/6\approx0.7887\ge0.25$.

**(The conjecture as a research programme.)** The paper frames the result as "a preliminary validation of our conjecture: **the A-stable condition of an implicit Runge-Kutta method can guarantee a robust convergence rate $O(\alpha)$ for the preconditioned PinT iterative algorithm**". That conjecture is the explicit programme advanced by [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|papers 65 and 71]]: paper 65 proves it for all stable one-step methods, and paper 71 confirms the same pattern in the coupled forward-backward setting at the price of letting $\alpha$ vary with $N_t$.

### Numerical experiments

Both test problems are pointedly chosen:

- **Advection-dominated diffusion equations.** A hard case for parallel-in-time. Parareal-type methods work because dissipation dominance makes the coarse and fine propagators agree at high frequency, and advection dominance destroys exactly that (see the data on the [[en/computational-mathematics/paper-notes/parallel-in-time/parareal-convergence|parareal convergence]] page, where parareal diverges once the viscosity drops to around $10^{-3}$). The $\alpha$-circulant bound $\rho\le\alpha/(1-\alpha)$ **does not depend on $\sigma(A)$**, so it ought not to degrade on such problems — which is what this example tests.
- **The viscous Burgers equation.** Nonlinear, testing applicability beyond the linear theory.

Both support the theoretical finding. Reading iteration counts of this kind requires the value of $\alpha$ alongside them: the trade-off table at the top of this page shows that the convergence rate and the attainable accuracy are determined entirely by $\alpha$. At $\alpha=10^{-2}$ the bound gives $\rho\le0.0101$, about $6$ sweeps to $10^{-10}$, with an accuracy floor near $2\times10^{-14}$.

### Relation to the others

Papers 39 and 31 supply the $\alpha$-circulant and diagonalisation toolkit; paper 53 is the first in the list to apply it to a **multi-stage** integrator, which required rebuilding the all-at-once formulation from scratch. [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|Paper 65]] is the natural next step on the single-step side: rather than one integrator at a time, it gives a uniform spectral analysis of $\mathcal P_\alpha^{-1}\mathcal K$ valid for **all** stable single-step integrators (first-order problems) and a large class of symmetric two-step methods (second-order problems). Papers 53 and 65 are therefore the two halves of the "stop doing case-by-case studies" programme, with 53 the last and hardest case study and 65 the general theorem. [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|Paper 71]] extends the same $\mathcal P_\alpha$-for-$\mathcal K$ idea from a single forward evolution to a coupled forward-backward pair, and [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|paper 84]] goes to the opposite extreme of discretisation type — time-spectral methods, whose all-at-once temporal blocks are **dense** rather than sparse blocks scrambled by stage coupling.

## The four papers ordered by price

| No. | What is diagonalised                     | Role of $\alpha$                                              | Conditioning cost                     | What it buys                                        |
| --- | ---------------------------------------- | ------------------------------------------------------------- | ------------------------------------- | --------------------------------------------------- |
| 31  | the circulant of a time-periodic problem | $\alpha=1$, exact periodicity                                 | $\mathrm{Cond}_2(V)=1$, free          | a direct (non-iterative) algorithm                  |
| 39  | the head-tail coupled coarse correction  | $\alpha\in(0,1)$, introduced                                  | $\mathrm{Cond}_2(V)\le1/\alpha$       | the rate is preserved when $\alpha\le\rho/(1+\rho)$ |
| 46  | the forward-backward optimality system   | a periodic surrogate $\widehat B_{\rm per}$ as preconditioner | absorbed by an outer Krylov iteration | eigenvalue and singular-value clustering            |
| 53  | a rebuilt multi-stage all-at-once system | $\alpha\in(0,1)$                                              | $\mathrm{Cond}_2(V)\le1/\alpha$       | $\rho\le\alpha/(1-\alpha)$ when $\gamma\ge1/4$      |

The table makes the logic of this route visible: **diagonalisation is not the goal in itself; making some matrix diagonalisable with a well-conditioned eigenvector matrix is.** The time-periodic problem supplies that for free. Everything else must trade a controllable approximation for it, and the approximation error is then handed to an outer Krylov iteration — which is the subject of the [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|all-at-once preconditioners]] page.

## Sources for this page

- S. Wu, H. Zhang, and T. Zhou, [_Solving time-periodic fractional diffusion equations via diagonalization technique and multigrid_](https://doi.org/10.1002/nla.2178), Numer. Linear Algebra Appl. 25(5) (2018), e2178.
- S. Wu and T. Zhou, [_Acceleration of the two-level MGRIT algorithm via the diagonalization technique_](https://doi.org/10.1137/18M1207697), SIAM J. Sci. Comput. 41(5) (2019), pp. A3421-A3448.
- S. Wu and T. Zhou, [_Diagonalization-based parallel-in-time algorithms for parabolic PDE-constrained optimization problems_](https://doi.org/10.1051/cocv/2020012), ESAIM Control Optim. Calc. Var. 26 (2020), 88.
- S. Wu and T. Zhou, [_Parallel implementation for the two-stage SDIRK methods via diagonalization_](https://doi.org/10.1016/j.jcp.2020.110076), J. Comput. Phys. 428 (2021), 110076.
