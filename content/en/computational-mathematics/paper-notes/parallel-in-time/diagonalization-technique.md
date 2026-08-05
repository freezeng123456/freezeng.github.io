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

> [!note] Coverage of this page
> Papers **31** (_Numer. Linear Algebra Appl._ 25(5), 2018), **39** (_SIAM J. Sci. Comput._ 41(5), 2019), **46** (_ESAIM Control Optim. Calc. Var._ 26, 2020) and **53** (_J. Comput. Phys._ 428, 2021). None has a preprint; the texts of papers 46 and 53 could not be retrieved, so those sections separate what the abstract confirms from the general machinery they rely on.

![Replace the sequential recurrence by a diagonalisable time matrix](assets/diagrams/tao-zhou-papers/en/pint-diagonalization.svg)

## The three-step pattern and its price

Stacking the unknowns at all time levels gives an all-at-once system in Kronecker form,

$$
\bigl(B_1\otimes I_x+B_2\otimes A\bigr)\boldsymbol u=\boldsymbol b,
\qquad
\boldsymbol u=(U_1^\top,\dots,U_{N_t}^\top)^\top .
$$

If the temporal matrix diagonalises as $C^{(\alpha)}_j=VD_jV^{-1}$, the solve becomes three steps:

$$
\text{(a) } S_1=(\mathbb F\Gamma_\alpha\otimes I_x)\boldsymbol b,
\qquad
\text{(b) } S_{2,n}=(\lambda_{1,n}I_x+\lambda_{2,n}A)^{-1}S_{1,n},
\qquad
\text{(c) } \boldsymbol u=(\Gamma_\alpha^{-1}\mathbb F^*\otimes I_x)S_2 ,
$$

with $\mathbb F$ the unitary DFT matrix and $\Gamma_\alpha=\mathrm{diag}(1,\alpha^{1/N_t},\dots,\alpha^{(N_t-1)/N_t})$. Steps (a) and (c) are FFTs and step (b) is fully parallel across the $N_t$ time levels.

The price sits in one place: $\mathrm{Cond}_2(V)\le1/\alpha$. Smaller $\alpha$ brings the $\alpha$-circulant matrix closer to the original lower-triangular Toeplitz matrix (a better approximation) and worsens the conditioning of $V$ (more roundoff amplification). The four papers here can be ordered by how they handle that price.

## 31: the time-periodic case makes the price zero

### An obstacle becomes an opportunity

Time-periodic diffusion equations are structurally different for time stepping: there is no initial condition to march from, and the periodicity condition $u(0)=u(T)$ couples the last time level back to the first. As the abstract puts it, one has to consider **all the discrete solutions at once** rather than one by one. Whether or not parallelism is the goal, one is forced into an all-at-once formulation. Combined with a fractional Laplacian, which makes the spatial operator dense, the space-time system is both large and expensive.

The observation is that this apparent obstacle is an opportunity. The periodic coupling makes the time-discretisation matrix a genuine **circulant** rather than a lower-triangular Toeplitz matrix, so it is unitarily diagonalisable by FFT. In the language of the $\alpha$-circulant family this is the special value $\alpha=1$:

$$
C_1^{(\alpha)}=\frac{1}{\Delta t}
\begin{bmatrix}1&&&-\alpha\\-1&1&&\\&\ddots&\ddots&\\&&-1&1\end{bmatrix},
\qquad
C_2^{(\alpha)}=
\begin{bmatrix}\theta&&&(1-\theta)\alpha\\1-\theta&\theta&&\\&\ddots&\ddots&\\&&1-\theta&\theta\end{bmatrix}.
$$

At $\alpha=1$ we have $\Gamma_\alpha=I$ and $V=\mathbb F^*$ is **unitary**, so $\mathrm{Cond}_2(V)=1$ and the roundoff obstruction that limits the other diagonalisation routes disappears entirely. That is why the time-periodic case admits a **direct** parallel-in-time algorithm with no parameter trade-off.

### Multigrid for the complex-shifted systems

Step (b) leaves a series of independent linear systems with **complex** coefficients, $(\lambda_{1,n}I_x+\lambda_{2,n}A)x=y$ with $A$ the discrete fractional Laplacian. The paper solves each with multigrid using **damped Richardson iteration** as the smoother. That choice has a concrete reason: the fractional-Laplacian matrix is dense, so a smoother needing only matrix-vector products is preferable to damped Jacobi or Gauss-Seidel, which need triangular solves.

The paper proves the linear solver has a mesh-independent convergence factor and optimises the Richardson damping parameter to minimise that constant.

> [!note] What could be verified
> The problem setting, the direct diagonalisation, the multigrid treatment of the complex-shifted systems, the mesh-independent convergence factor and the **existence** of the damping optimisation are all confirmable from the abstract. The optimal damping value, the minimised factor, and whether the analysis is uniform in the fractional order or in $N_t$ are unverified here.

## 39: fixing the coarse propagator and the coarse correction together

Paper 39 treats two-level MGRIT, whose convergence factor carries one extra factor relative to parareal:

$$
\varrho_l(J,z)=\frac{\bigl|R_f^{J}(z/J)\bigr|\,\bigl|R_g(z)-R_f^{J}(z/J)\bigr|}{1-|R_g(z)|}
\qquad\Longleftrightarrow\qquad
\varrho_{l,\mathrm{MGRIT}}=\bigl|R_f^J(z/J)\bigr|\times\varrho_{l,\mathrm{parareal}} .
$$

The extra $|R_f^J(z/J)|\le1$ is precisely the gain from an additional F-relaxation. On terminology: F-relaxation applies the fine propagator on the F-points inside each coarse interval, starting from the C-point value at its left end, and those solves are independent and parallel; C-relaxation updates the C-point values from the last F-point; FCF-relaxation is F, then C, then F again, which is why the iteration costs two fine solves. Two-level MGRIT with F-relaxation only **is** parareal.

### Two inefficiencies

First, the coarse propagator $\mathcal G$ is almost always backward Euler, which is only first order and gives a contraction of about $0.1$ under FCF-relaxation (the exact values appear below); nobody had asked whether a better $\mathcal G$ could do materially better **without making the coarse solve dominant**. The paper switches to the second-order Lobatto IIIC method, whose stability function is

$$
R_g(z)=\frac{1}{1+z+z^2/2}
$$

in the convention $\boldsymbol u'+A\boldsymbol u=g$. This is the $(0,2)$ Padé approximant: A-stable, **L-stable** ($R_g(\infty)=0$), second order, and stiffly accurate in the damping sense — it damps high-frequency modes far more aggressively than backward Euler's $1/(1+z)$, which is exactly what a coarse propagator needs.

Second, and structurally, the coarse-grid correction is an inherently **sequential** sweep through $N_t$ coarse levels, so by Amdahl's law it becomes the serial bottleneck as the processor count grows. A better $\mathcal G$ normally worsens this, since a higher-order implicit Runge-Kutta coarse solver has several stages and costs several backward-Euler steps.

### A head-tail coupling makes the coarse correction diagonalisable

The paper replaces the sequential coarse correction by one with a **head-tail coupling condition**:

$$
\boldsymbol u_{n+1}^{k+1}=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\underbrace{\mathcal F(T_n,T_{n+1},\tilde{\boldsymbol u}_n^{k})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k})}_{=:\,\boldsymbol b_{n+1}^k},
\qquad
\boldsymbol u_0^{k+1}=\alpha\,\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0 ,
$$

with $\tilde{\boldsymbol u}_n^k=\boldsymbol u_n^k$ for $n\ge1$ and $\tilde{\boldsymbol u}_0^k=\boldsymbol u_0$; that redefinition is necessary or the fixed point is no longer the true solution. For backward Euler as $\mathcal G$, the coupled coarse correction is the all-at-once system

$$
\bigl(C_\alpha\otimes I_x-I_t\otimes\Delta T A\bigr)\boldsymbol U^{k+1}=\boldsymbol g^k ,
$$

diagonalised by FFT.

That the two fixes are compatible is worth emphasising: improving the coarse propagator alone worsens the sequential bottleneck, and parallelising the coarse correction cancels that worsening. Each change on its own delivers less than the two together.

### Four directly comparable constants

With $\mathcal F$ an L-stable integrator and $J=\Delta T/\Delta t=O(1)$, the values of $\max_{z\ge0}\varrho_l$ are

| Coarse propagator $\mathcal G$ | Parareal | MGRIT (FCF-relaxation) |
| ------------------------------ | -------- | ---------------------- |
| backward Euler                 | 0.2984   | 0.1115                 |
| Lobatto IIIC (second order)    | 0.0817   | 0.0197                 |

This is the precise form of the abstract's "from 0.1 to 0.02". The bound is also **robust**: it depends neither on the eigenvalues of the coefficient matrix nor on the coarsening ratio $J$.

One comparison has to be stated plainly or the gain will be overestimated: **at equal cost, one MGRIT-FCF iteration is slightly worse than two parareal iterations.** FCF-relaxation performs two fine solves per sweep, so the right comparison is against two parareal sweeps, and $0.2984^2=0.0890<0.1115$ while $0.0817^2=0.0067<0.0197$. The two columns therefore cannot be compared directly to decide which is better, because their per-iteration costs differ. What does hold firmly is the gain from **changing $\mathcal G$**: within either column, moving from backward Euler to Lobatto IIIC lowers the contraction factor by roughly a factor of 3.7 for parareal and 5.7 for MGRIT-FCF.

As for choosing the coupling parameter $\alpha$, the paper's conclusion is that a suitable choice leaves the new algorithm with the **same** convergence rate as the original. For the parareal case there is an explicit threshold, due to Wu (SISC 2018): provided

$$
\alpha\le\frac{\rho}{1+\rho},
$$

one has $\rho_{\text{new}}=\rho$. Since $\rho=O(10^{-1})$ in practice, $\alpha=O(10^{-1})$ suffices, and the diagonalisation roundoff amplification $\mathrm{Cond}_2(V)\le1/\alpha$ is then negligible. **Whether paper 39 establishes the same threshold for MGRIT has not been verified here**; only the statement that a suitable $\alpha$ preserves the convergence rate is confirmed.

## 46: a forward-backward system has no single direction of propagation

The difficulty in parabolic PDE-constrained optimisation is stated plainly in the abstract: the computation has to take the discrete time points **all-at-once**, and the coupled system has **opposite evolution directions** — the state equation runs forward from an initial condition, the adjoint runs backward from a terminal condition. That is the defining obstruction. Sequential time stepping does not apply to the coupled system at all, and parareal-style remedies do not either, because there is no single direction of propagation to iterate along.

The standard distributed-control problem minimises

$$
\mathcal J(y,q)=\tfrac12\|y-y_d\|^2_{L^2(\Omega\times(0,T))}
+\tfrac{\gamma}{2}\|q\|^2_{L^2(\Omega\times(0,T))}
\quad\text{s.t.}\quad
\partial_t y-\Delta y=f+q,
$$

and eliminating the control gives the forward-backward first-order optimality pair

$$
\begin{cases}
\partial_t y-\Delta y-\tfrac1\gamma p=f,\\
-\partial_t p-\Delta p+y=y_d,
\end{cases}
$$

with either initial conditions (the initial-value model) or periodicity (the time-periodic model).

The paper's organising objects are two **time discretisation matrices**, $B_{\mathrm{per}}$ for the time-periodic problem and $B_{\mathrm{ini}}$ for the initial-value problem, and the abstract states that the main idea lies in carefully handling them. Earlier work built preconditioners by **approximating the Schur complement** of the discrete KKT system, and the paper's point of departure is that those converge much more slowly than a diagonalisation-based construction can.

The relation to paper 31 can be stated in one sentence: **paper 31 is the observation that $\alpha=1$ is free, and paper 46 is the decision to approximate the non-periodic problem by a periodic one.** The former diagonalises $B_{\mathrm{per}}$ itself; the latter uses a periodic-like approximation $\widehat B_{\mathrm{per}}$ as a **preconditioner** for the initial-value matrix $B_{\mathrm{ini}}$.

> [!note] What could be verified
> The obstruction as stated, the organising role of the two time-discretisation matrices, and the comparison against Schur-complement preconditioners are confirmable from the abstract. The paper's exact cost-functional notation, the precise form of the preconditioner, its theorems and numerical results are unverified here; the publisher's full text returned HTTP 403 to automated retrieval.

## 53: diagonalisation across stages instead of time levels

Paper 53 applies the same mechanism in a different direction: not across time levels but across the **stages within a single time step**. A two-stage singly diagonally implicit Runge-Kutta method solves two implicit stages per step, and those stages are normally solved in sequence. Diagonalisation decouples them so both can be solved concurrently.

> [!note] What could be verified
> The text could not be retrieved. Confirmable: the object named in the title, a parallel implementation of two-stage SDIRK methods by diagonalisation. Its specific diagonalisation form, conditioning analysis and parallel-efficiency results are unverified here.

## The four papers ordered by price

| No. | What is diagonalised                     | Role of $\alpha$ or the step sizes         | Conditioning cost               |
| --- | ---------------------------------------- | ------------------------------------------ | ------------------------------- |
| 31  | the circulant of a time-periodic problem | $\alpha=1$, exact periodicity              | $\mathrm{Cond}_2(V)=1$, free    |
| 39  | the head-tail coupled coarse correction  | $\alpha\in(0,1)$, introduced               | $\mathrm{Cond}_2(V)\le1/\alpha$ |
| 46  | the forward-backward optimality system   | a periodic approximation as preconditioner | absorbed by an outer iteration  |
| 53  | the stages within one step               | no $\alpha$ involved                       | set by the stage structure      |

The table makes the logic of this route visible: **diagonalisation is not the goal in itself; making some matrix diagonalisable with a well-conditioned eigenvector matrix is.** The time-periodic problem supplies that for free. Everything else must trade a controllable approximation for it, and the approximation error is then handed to an outer Krylov iteration — which is the subject of the [[en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|all-at-once preconditioners]] page.

## Coverage check

| Item                                             | Paper | Status                                                                         |
| ------------------------------------------------ | ----- | ------------------------------------------------------------------------------ |
| All-at-once Kronecker system, three steps        | 31    | system, three steps, FFTs and parallel step                                    |
| $\alpha$-circulant matrices and $\alpha=1$       | 31    | both matrix forms, unitarity, unit condition number                            |
| Multigrid for complex shifts, smoother choice    | 31    | reason for damped Richardson, mesh-independent factor                          |
| MGRIT convergence factor and FCF relaxation      | 39    | the extra factor, all three relaxations, relation to parareal                  |
| Lobatto IIIC as coarse propagator                | 39    | stability function, Padé order, L-stability and damping                        |
| Head-tail coupling and diagonalisable correction | 39    | coupling condition, why the redefinition is needed, all-at-once system         |
| Forward-backward obstruction, two matrices       | 46    | opposite directions, $B_{\mathrm{per}}$ and $B_{\mathrm{ini}}$                 |
| Logical relation to paper 31                     | 46    | from "free $\alpha=1$" to "approximate the non-periodic by a periodic problem" |

## Sources for this page

- S. Wu, H. Zhang, and T. Zhou, [_Solving time-periodic fractional diffusion equations via diagonalization technique and multigrid_](https://doi.org/10.1002/nla.2178), Numer. Linear Algebra Appl. 25(5) (2018), e2178.
- S. Wu and T. Zhou, [_Acceleration of the two-level MGRIT algorithm via the diagonalization technique_](https://doi.org/10.1137/18M1207697), SIAM J. Sci. Comput. 41(5) (2019), pp. A3421-A3448.
- S. Wu and T. Zhou, [_Diagonalization-based parallel-in-time algorithms for parabolic PDE-constrained optimization problems_](https://doi.org/10.1051/cocv/2020012), ESAIM Control Optim. Calc. Var. 26 (2020), 88.
- S. Wu and T. Zhou, [_Parallel implementation for the two-stage SDIRK methods via diagonalization_](https://doi.org/10.1016/j.jcp.2020.110076), J. Comput. Phys. 428 (2021), 110076.
