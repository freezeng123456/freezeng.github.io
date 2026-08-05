---
title: IMEX and Relaxation Schemes
description: Papers 78, 91 and 104 - treating the nonlinearity explicitly while keeping the energy argument
lang: en
translation: computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation
tags:
  - paper-notes
  - phase-field
  - implicit-explicit-schemes
---

> [!note] Coverage of this page
> Papers **78** (_J. Comput. Phys._ 515, 2024), **91** (_Math. Comput._ 95(359), 2026) and **104** (submitted to _SIAM J. Numer. Anal._, [arXiv:2605.05619](https://arxiv.org/abs/2605.05619)).

Implicit-explicit schemes treat the stiff linear part implicitly and the nonlinearity explicitly, so no inner iteration is needed per stage. The price is that the energy argument becomes harder: an explicitly treated nonlinearity no longer supplies a sign-definite contribution automatically. These three papers give three different answers.

## 91: making the average dissipation rate independent of the discretisation parameters

### Two defects

The Cahn-Hilliard model is the $H^{-1}$ gradient flow of the Ginzburg-Landau free energy,

$$
E[\Phi]=\int_\Omega\Bigl[\tfrac{\epsilon^2}{2}|\nabla\Phi|^2+F(\Phi)\Bigr]\mathrm d\mathbf x,
\qquad F(\Phi)=\tfrac14(\Phi^2-1)^2,
$$

$$
\partial_t\Phi=\Delta\bigl[F'(\Phi)-\epsilon^2\Delta\Phi\bigr],
\qquad
\frac{\mathrm dE}{\mathrm dt}
=-\bigl((-\Delta)^{-1}\partial_t\Phi,\partial_t\Phi\bigr)_{L^2}\le0 .
$$

Proving that implicit-explicit Runge-Kutta methods inherit the **original** energy dissipation law requires knowing first that the **stage solutions are uniformly bounded in the maximum norm** before the stabilised-splitting energy argument can close. The earlier unified framework simply **assumed** the nonlinear bulk $F'$ is globally Lipschitz, which the Cahn-Hilliard quartic $F(\Phi)=\frac14(\Phi^2-1)^2$ is not.

The second defect is practical: the average dissipation rate of a generic implicit-explicit Runge-Kutta method depends on $\tau_n\overline\lambda_{\mathrm{ML}}$, the time step times the average eigenvalue of the discrete operator, so as soon as an adaptive algorithm takes a large step the effective dissipation changes and the computed energy curve drifts.

### Differential form and the refinement condition

Rewriting with the stabilisation $L_\kappa\Phi=-\epsilon^2\Delta\Phi+\kappa\Phi$ and $f_\kappa(\Phi)=\kappa\Phi-F'(\Phi)$, the $s$-stage scheme is

$$
u_h^{n,i}=u_h^{n,1}
+\tau_n\sum_{j=1}^{i}a_{ij}\Delta_hL_{\kappa,h}u_h^{n,j}
-\tau_n\sum_{j=1}^{i-1}\hat a_{ij}\Delta_h f_\kappa(u_h^{n,j}),
$$

with the implicit part $A$ a stiffly accurate diagonally implicit Runge-Kutta method with explicit first stage (first same as last), $\widehat A$ strictly lower triangular, and the **node condition** $\hat{\mathbf c}=\mathbf c$ (equivalently $A\mathbf 1=\widehat A\mathbf 1$), which makes the method consistent at all stages and preserves equilibria.

The paper rewrites the scheme in differential form and introduces a matrix

$$
D(z)=D_{\mathrm E}-zD_{\mathrm{EI}},
\qquad
D_{\mathrm E}=A_{\mathrm E}^{-1}E_{s_{\mathrm I}},
\qquad
D_{\mathrm{EI}}=A_{\mathrm E}^{-1}A_{\mathrm I}E_{s_{\mathrm I}}-E_{s_{\mathrm I}}+\tfrac12 I_{s_{\mathrm I}},
$$

with $E_{s_{\mathrm I}}$ the lower-triangular all-ones matrix, and calls $D$ positive (semi-)definite when its symmetric part is. **This $D$ plays exactly the role that the discrete orthogonal convolution kernels play in the [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|BDF family]]**, with the index changed from time levels to stages; the corresponding orthogonality identity is

$$
\sum_{i=j}^{k}d^{(R)}_{k,i}\,\underline{\hat a}_{i+1,j}\equiv\delta_{kj},
\qquad
(\underline{\hat a}_{i+1,j})=E_{s_{\mathrm I}}^{-1}A_{\mathrm E} .
$$

Given positive semi-definiteness of $D_{\mathrm E}$ and $D_{\mathrm{EI}}$, the stage energy law carries the average dissipation rate

$$
\mathcal R=\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}\frac1{\hat a_{k+1,k}}
+\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}
\Bigl(\frac{a_{k+1,k+1}}{\hat a_{k+1,k}}-\frac12\Bigr)\tau_n\overline\lambda_{\mathrm{ML}}
\ \ge 0 .
$$

A method is "good" when $\mathcal R$ stays close to $1$ over a wide range of $\tau_n\overline\lambda_{\mathrm{ML}}$.

**The central observation is that $\mathcal R$ is independent of $\tau_n\overline\lambda_{\mathrm{ML}}$ if and only if $D_{\mathrm{EI}}=\mathbf 0$**, that is

$$
A_{\mathrm I}=A_{\mathrm E}P_{s_{\mathrm I}},
\qquad
P_{s_{\mathrm I}}=I_{s_{\mathrm I}}-\tfrac12E_{s_{\mathrm I}}^{-1} .
$$

Then $D_{\mathrm R}=D_{\mathrm E}=A_{\mathrm E}^{-1}E_{s_{\mathrm I}}$ is independent of $z$ and therefore of both mesh parameters, with

$$
\mathcal R_{\mathrm R}=\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}\frac1{\hat a_{k+1,k}},
$$

and the scheme collapses to the compact refined form

$$
u_h^{n,i+1}=u_h^{n,1}+\tau_n\sum_{j=1}^{i}\hat a_{i+1,j}\Delta_h
\Bigl[L_{\kappa,h}u_h^{n,j+\frac12}-f_\kappa(u_h^{n,j})\Bigr].
$$

### A structural exclusion

The node condition forces the first implicit column $\mathbf a_1=(\tfrac12\hat a_{21},\dots,\tfrac12\hat a_{s1})^T\ne\mathbf 0$, so such methods are necessarily **Lobatto-type and never Radau- or ARS-type**. The paper states this as a proposition: no Radau-type or ARS-type implicit-explicit Runge-Kutta method has an average dissipation rate independent of $\tau_n\overline\lambda_{\mathrm{ML}}$.

Results of this kind are valuable in method design: the conclusion is not "we chose Lobatto-type" but "this property forces Lobatto-type".

Two concrete methods: the two-stage case with $\theta=1/2$ gives $\mathcal R_{\mathrm R}=1$ exactly and a Crank-Nicolson-type scheme

$$
\delta_\tau\phi_h^{n}=\tau_n\Delta_h\Bigl[\tfrac12L_{\kappa,h}(\phi_h^{n}+\phi_h^{n-1})
-f_\kappa(\phi_h^{n-1})\Bigr];
$$

at second order no three-stage method exists, and the four-stage one-parameter family sets $\hat a_{32}=\hat a_{43}=c_2$ and $\hat a_{42}=\frac1{2c_2}-c_3$ with $c_3$ a root of

$$
c_3^2-\Bigl(\frac1{2c_2}+c_2\Bigr)c_3+(c_2-1)^2=0,
$$

valid for $0<c_2\le\frac{2+\sqrt6}{2}$, with positive definiteness decided by two explicit determinant conditions.

## 104: replacing a family-by-family search for multipliers

### The problem

For nonlinear parabolic equations, implicit-explicit multistep methods treat the stiff linear operator implicitly and the nonlinearity explicitly. Many families exist — weighted BDF, modified BDF, generalised BDF, NIMEX — all designed to enlarge absolute stability regions. But third- and higher-order variants are **not A-stable**, so rigorous **discrete energy** stability was open.

Both available routes require bespoke objects found family by family. The first is Nevanlinna-Odeh-type multipliers via Dahlquist G-stability, which needs a different multiplier family for each of weighted, modified and generalised BDF and NIMEX. The second is an implicit-part decomposition of Huang-Shen type,

$$
\sum_{j=0}^{k-1}b_{\mathrm G,j}^{(k,\beta)}v^{n-j}
=\eta_k(\beta)\sum_{j=0}^{k-1}c_{\mathrm G,j}^{(k,\beta)}v^{n-j}
+\sum_{j=0}^{k-1}d_j^{(k,\beta)}v^{n-j},
$$

with $\eta_2=\frac{\beta-1}{\beta}$, $\eta_3=\frac{\beta-1}{\beta+1}$, $\eta_4=\frac{\beta-1}{\beta+3}$, giving stability under $\eta_k(\beta)>\mu_0/\varpi$. The limitations are concrete: it degenerates at $\beta=1$ and therefore **never covers plain BDF-$k$**; it exists only for $2\le k\le4$; its refined three-term version works only at fixed $\beta_k=3,6,9$; and there is no analogue for generalised BDF5 or for weighted BDF, modified BDF and NIMEX.

**This paper replaces all of that with one computation: three extremal values of three explicit rational functions on the unit circle.**

### The abstract setting and the global discrete energy method

On a Hilbert triple $V\subset H=H'\subset V'$ consider $u_t+\varpi\mathcal Lu=\mathcal F(u)$ with $\mathcal L$ positive definite, self-adjoint, bounded and linear, $\mathcal F$ possibly nonlinear, and the **local** Lipschitz hypothesis

$$
\|\mathcal F(v)-\mathcal F(w)\|_{\star}\le\mu_0\|v-w\|_V+\mu_1\|v-w\|_H,
\qquad \mu_0\in(0,\varpi).
$$

The framework also covers a non-self-adjoint $\mathcal L=\mathcal L_s+\mathcal L_a$ by moving the low-order anti-self-adjoint part into the explicit term.

A general $k$-step implicit-explicit multistep method on a uniform mesh reads

$$
\sum_{j=0}^{k-1}a_j^{(k)}\partial_\tau u^{n-j}
+\varpi\sum_{j=0}^{k}b_j^{(k)}\mathcal Lu^{n-j}
=\sum_{j=0}^{k-1}c_j^{(k)}\mathcal F(u^{n-j-1})+\mathfrak C^{(k)}_n(u^0),
$$

so a method is a triad $(\vec a^{(k)},\vec b^{(k)},\vec c^{(k)})$. Using the discrete orthogonal convolution kernels of $\vec a^{(k)}$,

$$
a_0^{(-1,k)}=\frac{1}{a_0^{(k)}},
\qquad
a_j^{(-1,k)}=-\frac{1}{a_0^{(k)}}\sum_{i=1}^{j}a^{(-1,k)}_{j-i}a_i^{(k)},
$$

which satisfy two-sided orthogonality $\sum_{\ell=j}^{n}a^{(-1,k)}_{n-\ell}a^{(k)}_{\ell-j}\equiv\delta_{nj}$, the scheme becomes a "differential" form

$$
\partial_\tau u^{n}+\varpi\sum_{\ell=1}^{n}\hat b^{(k)}_{n-\ell}\mathcal Lu^{\ell}
=\sum_{\ell=1}^{n}\hat c^{(k)}_{n-\ell}\mathcal F(u^{\ell-1})
+\sum_{\ell=1}^{n}a^{(-1,k)}_{n-\ell}\mathfrak C^{(k)}_\ell(u^0),
$$

with composite kernels $\hat b_j^{(k)}=\sum_{i=0}^{j}a^{(-1,k)}_{j-i}b_i^{(k)}$ and $\hat c_j^{(k)}=\sum_{i=0}^{j}a^{(-1,k)}_{j-i}c_i^{(k)}$. In matrix form $\widehat B_{L,k}=A_{L,k}^{-1}B_{L,k}$ and $\widehat C_{L,k}=A_{L,k}^{-1}C_{L,k}$, both lower-triangular Toeplitz and commuting. Because nothing is discarded, the authors call this the **global discrete energy method**.

The route is the same as in [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|paper 58]]: use discrete orthogonal convolution kernels to invert the multistep formula back to a single first difference and recover a textbook energy argument. The difference is that here three kernel sequences are in play at once — implicit, explicit and differential — so the criterion must be composite, which is exactly what the semi-generating function approach handles.

## 78: linear relaxation with a regularised energy reformulation

Paper 78 takes a different route: **do not improve the time discretisation, rewrite the energy.** The paper names its method RRER, for relaxation with regularized energy reformulation.

### The problem: IEQ and SAV both differentiate the auxiliary variable in time

Invariant energy quadratisation (IEQ) and the scalar auxiliary variable (SAV) method both rewrite a nonlinear energy as a quadratic form in an auxiliary variable, so that each step needs only a linear solve. IEQ typically produces a **coupled** system with **time-dependent** coefficients; SAV produces a **decoupled** system with **constant** coefficients. Both, however, share one move: the evolution equation for the auxiliary variable is obtained by **differentiating that variable in time**, which introduces its own truncation error and makes the discrete system less faithful to the original equation. The paper's own phrasing is that it does not need to take time derivatives of the auxiliary variables.

Its answer is to let the auxiliary variable be defined by an **algebraic** relation only, never differentiated in time.

### The regularised auxiliary variable

With free energy and gradient flow

$$
E(\phi)=\tfrac12(\mathcal L\phi,\phi)+\bigl(F(\phi),1\bigr),
\qquad
\frac{\partial\phi}{\partial t}=-\mathcal G\bigl(\mathcal L\phi+F'(\phi)\bigr),
$$

where $\mathcal L$ is linear and $\mathcal G\ge0$ is the mobility operator, one has $\frac{\mathrm d}{\mathrm dt}E=-\bigl(\mathcal L\phi+F'(\phi),\mathcal G(\mathcal L\phi+F'(\phi))\bigr)\le0$. Choose $C_0$ so that $F(\phi)+C_0\ge0$, introduce a **stabilisation parameter** $\gamma$, and set

$$
q=\sqrt{4\bigl(F(\phi)+C_0\bigr)}-\gamma
\qquad\Longrightarrow\qquad
F(\phi)=\tfrac14q^2+\tfrac12\gamma q+\tfrac{\gamma^2}{4}-C_0 .
$$

The equivalent system and its energy are

$$
\frac{\partial\phi}{\partial t}=-\mathcal G\mu,
\quad
\mu=\mathcal L\phi+\tfrac12qq'+\tfrac{\gamma}{2}q',
\quad
q=\sqrt{4\bigl(F(\phi)+C_0\bigr)}-\gamma,
$$

$$
\widehat E(\phi,q)=\tfrac12(\mathcal L\phi,\phi)
+\Bigl(\tfrac14q^2+\tfrac12\gamma q+\tfrac{\gamma^2}{4}-C_0,\,1\Bigr),
\qquad
\frac{\mathrm d}{\mathrm dt}\widehat E=-(\mu,\mathcal G\mu)\le0 .
$$

**One point must be emphasised: $\widehat E(\phi,q)$ is exactly equal to $E(\phi)$ at the continuous level, not merely an approximation to it.** The paper derives that identity line by line. This is unlike the modified energies of SAV and IEQ, which relate to the original energy only approximately, and it matters for the placement discussion below.

### A staggered time grid makes the scheme linear

Take the molecular beam epitaxy model **with** slope selection — **note that this is a different model from the one in paper 52**, see below:

$$
E(\phi)=\int_\Omega\Bigl(\frac{\epsilon^2}{2}(\Delta\phi)^2+\frac14\bigl(|\nabla\phi|^2-1\bigr)^2\Bigr)\mathrm d\mathbf x,
\qquad
\phi_t=-\epsilon^2\Delta^2\phi+\nabla\cdot\bigl((|\nabla\phi|^2-1)\nabla\phi\bigr),
$$

with $q=|\nabla\phi|^2-1-\gamma$ for $\gamma>0$. The Crank-Nicolson scheme lives on a **staggered** grid, with $\phi$ at integer levels and $q$ at half-integer levels:

$$
\begin{aligned}
&\text{(a)}\ \ \frac{\phi^{n+1}-\phi^{n}}{\delta t}
=-\epsilon^2\Delta g^{n+\frac12}
+\nabla\cdot\Bigl(q^{n+\frac12}\nabla\frac{\phi^{n+1}+\phi^{n}}{2}\Bigr)
+\gamma g^{n+\frac12},\\
&\text{(b)}\ \ g^{n+\frac12}=\Delta\frac{\phi^{n+1}+\phi^{n}}{2},
\qquad
\text{(c)}\ \ \frac{q^{n+\frac12}+q^{n-\frac12}}{2}=|\nabla\phi^{n}|^2-1-\gamma .
\end{aligned}
$$

**What makes the scheme linear is (c): it is algebraic, and its right-hand side involves only the known $\phi^n$.** Each step therefore solves one linear algebraic system. The scheme is second-order in time and carries no step-ratio restriction at all — the steps are uniform and the stability is unconditional.

The paper treats the phase-field crystal model the same way, with $q=\phi^2-b_0-\gamma$, and extends RRER to a ternary phase-field model and a grain-growth model.

### Two theorems

For the molecular beam epitaxy model the scheme **conserves mass**, with $\int_\Omega\phi^{n+1}=\int_\Omega\phi^n$, and is **unconditionally energy stable**. The numerical section uses $P_1$ finite elements in FreeFEM: a manufactured-solution test gives temporal rates of $2.00,2.00,2.00$ against $1.94,1.97,1.99$ for IEQ; for the phase-field crystal model RRER is more accurate than both IEQ and exponential SAV while also the cheapest in CPU time; and the coarsening, stripe and triangle patterns, and the tests on spheres and tori all agree with the phase diagrams in the literature.

### Why this paper is the exception in this topic

This is worth stating plainly, because it bears on a tension running through the whole topic. Paper 78 has no Liao co-authorship, no DOC or DCC kernels, no variable steps and no step-ratio analysis; it belongs to the IEQ/SAV/relaxation branch rather than the convolution-kernel branch.

More importantly, **the type of its energy statement differs**: paper 78 proves dissipation of a **modified** energy $\widehat E(\phi,q)$, whereas papers 40, 48, 52, 57 and 91 work hard precisely to **avoid** that, establishing dissipation of the original or variational energy instead. The tension is explicit in the literature: paper 91's introduction criticises SAV-based high-order schemes for establishing stability only with respect to a modified energy involving the auxiliary variable. What mitigates it here is that $\widehat E\equiv E$ holds exactly at the continuous level, so this modified energy is not a new object but another way of writing the same one.

## The three approaches side by side

| No. | What is treated explicitly | How the energy argument is kept                               |
| --- | -------------------------- | ------------------------------------------------------------- |
| 78  | the nonlinear potential    | reformulate the energy (regularised quadratisation)           |
| 91  | the nonlinear bulk         | design the scheme so the dissipation rate is mesh-independent |
| 104 | the nonlinear term         | supply one criterion (extrema of three rational functions)    |

What the three share is that none tries to prove directly that the explicit term contributes with a definite sign; each substitutes a controllable object instead. Paper 78 substitutes the energy, paper 91 substitutes a structural condition on the Butcher tableau, and paper 104 substitutes the form of the criterion. **That matches the topic's overall stance: locate the algebraic object carrying the whole difficulty, then build independent criteria for it.**

## Coverage check

| Item                                              | Paper | Status                                                                    |
| ------------------------------------------------- | ----- | ------------------------------------------------------------------------- |
| Cahn-Hilliard energy and $H^{-1}$ gradient flow   | 91    | energy, equation, dissipation law                                         |
| The two defects                                   | 91    | failed global Lipschitz assumption; drifting rate                         |
| Stabilisation, node condition, differential form  | 91    | $L_\kappa$, $\hat{\mathbf c}=\mathbf c$, $D(z)$                           |
| Correspondence between $D$ and DOC kernels        | 91    | orthogonality identity with the index changed                             |
| Average dissipation rate and refinement condition | 91    | $\mathcal R$, $D_{\mathrm{EI}}=\mathbf 0$, $A_{\mathrm I}=A_{\mathrm E}P$ |
| Why Lobatto-type is forced                        | 91    | nonzero first column and the exclusion proposition                        |
| Both concrete methods                             | 91    | two-stage Crank-Nicolson type; four-stage family and its range            |
| Concrete limits of the two existing routes        | 104   | family-by-family multipliers; decomposition degenerating at $\beta=1$     |
| Abstract setting and local Lipschitz condition    | 104   | Hilbert triple, $\mu_0\in(0,\varpi)$, non-self-adjoint case               |
| Three kernel sequences and the global method      | 104   | the triad, DOC kernels, composite kernels, matrix form                    |

## Sources for this page

- J. Zhang, X. Guo, M. Jiang, T. Zhou, and J. Zhao, [_Linear relaxation method with regularized energy reformulation for phase field models_](https://doi.org/10.1016/j.jcp.2024.113225), J. Comput. Phys. 515 (2024), 113225.
- H.-l. Liao, T. Tang, X. Wang, and T. Zhou, [_A class of refined implicit-explicit Runge-Kutta methods with robust time adaptability and unconditional convergence for the Cahn-Hilliard model_](https://doi.org/10.1090/mcom/4090), Math. Comput. 95(359) (2026), pp. 1293-1325 (preprint [arXiv:2412.07321](https://arxiv.org/abs/2412.07321)).
- H.-l. Liao, C. Quan, T. Tang, and T. Zhou, _A semi-generating function approach to the stability of implicit-explicit multistep methods for nonlinear parabolic equations_, [arXiv:2605.05619](https://arxiv.org/abs/2605.05619), submitted to SIAM J. Numer. Anal.
