---
title: All-at-Once Preconditioners and Spectra
description: Papers 59, 65, 71, 84 and 85 - turning "here is a preconditioner" into "here is why its spectrum is controlled"
lang: en
translation: computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners
tags:
  - paper-notes
  - parallel-in-time
  - preconditioning
---

> [!note] Coverage of this page
> Papers **59** (_Adv. Comput. Math._ 48:16, 2022), **65** (_SIAM J. Matrix Anal. Appl._ 43(3), 2022), **71** (_SIAM J. Matrix Anal. Appl._ 44(4), 2023), **84** (_J. Sci. Comput._ 103:82, 2025) and **85** (_Acta Numer._ 34, 2025). The section-by-section reading of paper 85 is a separate topic, [[en/computational-mathematics/knowledge-notes/time-parallelization/index|time parallelization]]; this page only places it in this thread. The abstract of paper 71 has all mathematical symbols stripped in every public source, which that section states.

## 59: change the time discretisation so that $V$ is well conditioned by construction

### The ceiling of the Maday-Rønquist route

For $u'(t)+Au(t)=g(t)$ one forms the all-at-once system $\mathcal M\boldsymbol u:=(B\otimes I_x+I_t\otimes A)\boldsymbol u=\boldsymbol b$, and if $B=VDV^{-1}$ then

$$
\mathcal M=(V\otimes I_x)\bigl(D\otimes I_x+I_t\otimes A\bigr)(V^{-1}\otimes I_x),
$$

whose middle step is $n$ fully decoupled spatial solves. **The obstruction is that a standard uniform-step discretisation gives a $B$ that is not diagonalisable at all**: for backward Euler with uniform $\Delta t$, $B$ is lower bidiagonal Toeplitz, a single Jordan block, and for multistep methods $B$ is lower-triangular Toeplitz and likewise defective.

Maday and Rønquist's fix is **distinct** step sizes $\{\Delta t_j\}$, which makes the $n$ diagonal entries distinct and hence $B$ diagonalisable. But then the roundoff error obeys

$$
\texttt{roundoff error}=\mathcal O\bigl(\epsilon\,\mathrm{Cond}_2(V)\bigr),
\qquad \epsilon=\text{machine precision},
$$

and with geometrically increasing steps $\Delta t_j=\Delta t_1\tau^{\,n-j}$ the parameter $\tau>1$ is caught in a vice: $\tau\to1$ makes $B$ nearly defective and $\mathrm{Cond}_2(V)$ explode, while $\tau\gg1$ makes the step sizes grow exponentially and destroys the discretisation accuracy. Balancing the two errors limits the method to roughly $n\approx20$ to $25$ time points, which caps the achievable parallelism. Removing that cap is the paper's objective.

### A boundary value method: use an unstable scheme only all-at-once

The paper's time discretisation is hybrid: a centred (leapfrog) difference for the first $n-1$ steps and implicit Euler for the last step only,

$$
\begin{cases}
\dfrac{u_{j+1}-u_{j-1}}{2\Delta t}+Au_j=g_j, & j=1,2,\dots,n-1,\\[6pt]
\dfrac{u_n-u_{n-1}}{\Delta t}+Au_n=g_n. &
\end{cases}
$$

This is a **boundary value method**: it must **not** be run as a time-stepping scheme, since the centred scheme is unstable that way; it is only meaningful solved all-at-once. The construction is due to Axelsson and Verwer (1985), who proved that the simultaneously obtained solutions have **uniform second-order accuracy** even though the last step is only first order. Brugnano, Mazzia and Trigiante (1993) solved the resulting all-at-once system iteratively; this paper solves it **directly** by diagonalisation.

The key gain is that the steps can be **uniform**, so there is no step-ratio parameter at all:

$$
B=\frac{1}{\Delta t}
\begin{bmatrix}
0&\tfrac12&&&\\
-\tfrac12&0&\tfrac12&&\\
&\ddots&\ddots&\ddots&\\
&&-\tfrac12&0&\tfrac12\\
&&&-1&1
\end{bmatrix}.
$$

### Why $V$ is conditioned like $\mathcal O(n^2)$: the Chebyshev connection

The result is easy to remember, but the mechanism is the substance of this paper: **the eigenproblem reduces to Chebyshev polynomials, which yields explicit formulas for both $V$ and $V^{-1}$.** Work with the rescaled $\mathcal B:=\Delta t\,B$ and write the Chebyshev polynomials of the first and second kind,

$$
T_n(x)=\cos(n\arccos x),
\qquad
U_n(x)=\frac{\sin\bigl[(n+1)\arccos x\bigr]}{\sin(\arccos x)} .
$$

The $n$ eigenvalues of $\mathcal B$ are then $\lambda_j=\mathrm ix_j$, where $\{x_j\}_{j=1}^n$ are the $n$ roots of the **characteristic equation**

$$
U_{n-1}(x)-\mathrm i\,T_n(x)=0,
$$

and the components of the corresponding eigenvector are

$$
p_{j,k}=\mathrm i^{\,k}\,U_k(x_j),\qquad k=0,1,\dots,n-1,\qquad p_{j,0}=1 .
$$

The structure of these roots supplies exactly the three things needed: all of them are **simple**, so $\mathcal B$ is diagonalisable; all are complex with negative imaginary part, so the eigenvalues sit in the correct half-plane; and all have modulus below $1+1/\sqrt{2n}$, which is what controls $\mathrm{Cond}_2(V)$. If $x$ is a root then so is $-\bar x$. The proof substitutes $y=x+\sqrt{x^2-1}$, so that $T_n=\tfrac12(y^n+y^{-n})$ and $U_{n-1}=(y^n-y^{-n})/(y-y^{-1})$, reducing the characteristic equation to

$$
y^{2n}=-\frac{(y-\mathrm i)^2}{(y+\mathrm i)^2},
$$

and then uses $|y|>1$ together with the Chebyshev Pythagorean identity $T_n^2(x)+(1-x^2)U_{n-1}^2(x)=1$.

The eigenvector matrix accordingly factors into a **unitary** diagonal term times a Chebyshev-Vandermonde-like matrix,

$$
V=\underbrace{\mathrm{diag}\bigl(\mathrm i^0,\mathrm i^1,\dots,\mathrm i^{\,n-1}\bigr)}_{\text{unitary}}\;
\underbrace{\begin{bmatrix}
U_0(x_1)&\cdots&U_0(x_n)\\
\vdots&&\vdots\\
U_{n-1}(x_1)&\cdots&U_{n-1}(x_n)
\end{bmatrix}}_{=:\,\Phi},
$$

so $\mathrm{Cond}_2(V)=\mathrm{Cond}_2(\Phi)$ — the unitary factor is free, and the conditioning question reduces entirely to $\Phi$. From this one proves

$$
\mathrm{Cond}_2(V)=\mathcal O(n^2),
$$

**polynomial** rather than exponential growth in the number of time points, which replaces the $n\approx20$ ceiling with a mild polynomial cost. The paper additionally designs a fast algorithm exploiting this structure to compute the spectral decomposition of $\mathcal B$, and reports over 60 times speedup on 256 cores.

One structural convenience is worth recording as well: the all-at-once matrix for second-order (wave-type) equations is $B^2$, which shares **the same** $V$ as $B$, so the same conditioning analysis applies unchanged. That is why a single paper covers both the first- and second-order cases.

**The methodology of this paper deserves separate emphasis.** Earlier work took the time discretisation as given and looked for ways to cope with the conditioning of $V$; this paper inverts that and **chooses the time discretisation so that $V$ is well conditioned**. The price is giving up the time-stepping interpretation and accepting a scheme that only makes sense all-at-once.

## 65: replace every structural assumption with the classical one, stability

### The shape of the earlier theory

By 2022 block $\alpha$-circulant preconditioning had become a leading parallel-in-time method, especially for hyperbolic problems, because applying $\mathcal P_\alpha^{-1}$ reduces via FFT to independent solves at all time levels. But the theory had grown as a **pile of case-by-case studies**: each paper picked one time integrator and, exploiting that integrator's particular structure, worked out the spectrum of $\mathcal P_\alpha^{-1}\mathcal K$. The survey is blunt about it: the analyses are intricate and rely heavily on special properties of the time-stepping matrix such as **sparsity, Toeplitz structure and diagonal dominance**.

What was missing was one theorem covering all integrators, whose hypothesis is a **classical** property (stability) rather than an ad hoc structural one. A second gap: second-order (wave-type) problems were usually handled by rewriting $u''=Au+g$ as a first-order system, doubling the memory per time step — bad for fine spatial meshes and high dimensions — so a direct analysis of two-step discretisations of the second-order form was wanted.

### One notation for both settings

The first-order case writes a general one-step integrator in two-matrix form,

$$
r_1(\Delta t A)\,\boldsymbol u_n=r_2(\Delta t A)\,\boldsymbol u_{n-1}+\tilde{\boldsymbol g}_n,
\qquad n=1,\dots,N_t,
$$

so the method's stability function is $r_1^{-1}(z)r_2(z)$ (backward Euler: $r_1(z)=I-z$, $r_2(z)=I$). Stacking gives the all-at-once matrix and its block $\alpha$-circulant preconditioner

$$
\mathcal K=I_t\otimes r_1(\Delta tA)-B\otimes r_2(\Delta tA),
\qquad
\mathcal P_\alpha=I_t\otimes r_1(\Delta tA)-C_\alpha\otimes r_2(\Delta tA),
$$

with $B$ the strict down-shift and $C_\alpha=B+\alpha e_1e_{N_t}^\top$. In other words, **the preconditioner is obtained from the all-at-once matrix by putting $\alpha$ in the top-right corner**, replacing the (block) Toeplitz matrix by an $\alpha$-circulant one. Setting $\alpha=1$ gives an ordinary circulant and $\alpha\to0$ recovers $\mathcal K$ itself.

The second-order case avoids reduction to first order and uses a symmetric two-step method,

$$
r_1(\Delta t^2A)\,\boldsymbol u_{n+1}-r_2(\Delta t^2A)\,\boldsymbol u_n
+r_1(\Delta t^2A)\,\boldsymbol u_{n-1}=\tilde{\boldsymbol g}_n,
\qquad n=1,\dots,N_t-1,
$$

with all-at-once matrix and preconditioner

$$
\mathcal K=\tilde B\otimes r_1(\Delta t^2A)-B\otimes r_2(\Delta t^2A),
\qquad
\mathcal P_\alpha=\tilde C_\alpha\otimes r_1(\Delta t^2A)-C_\alpha\otimes r_2(\Delta t^2A),
$$

where the two-step stencil puts $\alpha$ in **two** corner positions.

### The main theorem, assuming only stability

Let $\mathcal K$ come from the one-step integrator above and let $\mathcal P_\alpha$ be the block $\alpha$-circulant preconditioner with $\alpha\in(0,1)$. If the integrator is **stable**, that is

$$
\bigl|r_1^{-1}(z)\,r_2(z)\bigr|\le1
\qquad \forall z\in\sigma(\Delta t A)\subset\mathbb C^-,
$$

then every eigenvalue of the preconditioned matrix satisfies the **mesh-independent** bound

$$
\frac{1}{1+\alpha}\ \le\ \bigl|\lambda(\mathcal P_\alpha^{-1}\mathcal K)\bigr|\ \le\ \frac{1}{1-\alpha} .
$$

The value of this is the exchange of hypotheses: from "the time-stepping matrix is sparse, Toeplitz and diagonally dominant" to "the method is stable". The latter is a property every textbook checks, so the theorem applies directly to any stable one-step method without redoing a spectral analysis for each.

Applying $\mathcal P_\alpha^{-1}$ uses $C_\alpha=VDV^{-1}$ with $V=\Gamma_\alpha^{-1}\mathbb F^*$ and $\Gamma_\alpha=\mathrm{diag}(1,\alpha^{1/N_t},\dots,\alpha^{(N_t-1)/N_t})$, so steps (a) and (c) are scaled FFTs and step (b) is $N_t$ independent complex spatial solves, with $\mathrm{Cond}_2(V)\le1/\alpha$ the roundoff price of small $\alpha$.

> [!warning] The direction of the bound
> The survey restates this theorem with the endpoints transposed, printing $\frac{1}{1-\alpha}\le|\lambda|\le\frac{1}{1+\alpha}$, which is impossible for $\alpha\in(0,1)$ since $1/(1-\alpha)>1/(1+\alpha)$. The display above has the correct orientation, for the reason given by the direct computation below.

### The direct computation in a single scalar channel

The orientation can be settled by calculation rather than taken on trust. Take one $z\in\sigma(\Delta tA)$ at a time and write $r:=r_1^{-1}(z)r_2(z)$, so that stability means $|r|\le1$. In that scalar channel $\mathcal K=I_t-rB$ and $\mathcal P_\alpha=I_t-rC_\alpha$ with $C_\alpha=B+\alpha e_1e_{N_t}^{\top}$, hence

$$
\mathcal K=\mathcal P_\alpha+\alpha r\,e_1e_{N_t}^{\top}
\qquad\Longrightarrow\qquad
\mathcal P_\alpha^{-1}\mathcal K=I_t+\alpha r\,\bigl(\mathcal P_\alpha^{-1}e_1\bigr)e_{N_t}^{\top},
$$

a **rank-one update of the identity**. Its eigenvalues are therefore $1$ with multiplicity $N_t-1$ together with the single value $1+\alpha r\,e_{N_t}^{\top}\mathcal P_\alpha^{-1}e_1$. Expanding $\mathcal P_\alpha^{-1}=\sum_{j\ge0}r^jC_\alpha^{\,j}$ and using $C_\alpha^{\,j}e_1=e_{1+j}$ for $j<N_t$ together with $C_\alpha^{N_t}=\alpha I_t$,

$$
e_{N_t}^{\top}\mathcal P_\alpha^{-1}e_1=\sum_{m\ge0}r^{\,N_t-1+mN_t}\alpha^m=\frac{r^{\,N_t-1}}{1-\alpha r^{\,N_t}},
\qquad
\lambda=1+\frac{\alpha r^{\,N_t}}{1-\alpha r^{\,N_t}}=\frac{1}{1-\alpha r^{\,N_t}} .
$$

Since $|r|\le1$ gives $|\alpha r^{N_t}|\le\alpha$ and hence $|1-\alpha r^{N_t}|\in[1-\alpha,1+\alpha]$, this is exactly $\frac1{1+\alpha}\le|\lambda|\le\frac1{1-\alpha}$; the same computation also gives the survey's own corollary $\rho(\mathcal M)=\bigl|\alpha r^{N_t}/(1-\alpha r^{N_t})\bigr|\le\alpha/(1-\alpha)$. Both follow from one step, which fixes the orientation.

The computation shows two further things in passing. First, in each scalar channel every eigenvalue but **one** equals $1$ exactly, so across the whole system at most $N_x$ eigenvalues differ from $1$ — the counterpart of McDonald, Pestana and Wathen's clustering theorem for $\alpha=1$. Second, the deviation carries the factor $r^{N_t}$, so for a **strictly** contractive channel ($|r|<1$, that is, with dissipation) the clustering is far tighter than the worst-case bound. **This explains why the method performs markedly better on parabolic than on hyperbolic problems**: in the hyperbolic case $|r|\approx1$, so $r^{N_t}$ does not decay and the worst-case bound is nearly attained.

> [!note] What this computation is
> The above is a derivation carried out for this page in order to settle the direction of the bound. It covers only the scalar channel and is not the paper's own proof. The paper appeared in _SIAM J. Matrix Anal. Appl._ and its theorem may be stated with sharper, $z$-dependent endpoints, so the journal text should still be consulted before quoting the paper's original formulation.

## 71: the forward-backward case has several Toeplitz blocks

In problems governed by a forward-backward evolutionary equation, the dominant computational cost is solving one large linear system whose central object is the all-at-once matrix of the **forward subproblem** after space-time discretisation. An efficient solver therefore needs a good preconditioner. The difficulty specific to the forward-backward case: the forward and backward evolutions are coupled and run in opposite directions, so neither can be time-stepped independently, and the resulting system is a saddle-point or Schur-complement system rather than a single block lower-triangular Toeplitz system. The $\alpha$-circulant theory of paper 65 covers a single forward evolution only.

The preconditioner construction can be stated in one line: **replace the Toeplitz matrices in the all-at-once matrix by $\alpha$-circulant matrices.** Note the plural — "the Toeplitz matrices" — consistent with a two-by-two block system in which several Toeplitz blocks (at least the forward matrix and its transpose) are each circulantised. The survey describes the result as ParaDiag-II applied to the forward-backward system producing a **parallel version** of the matching Schur complement preconditioner of Pearson, Stoll and Wathen.

The paper's two applications are PDE-constrained optimal control (the parabolic KKT system) and parabolic source identification (recovering an unknown source from observations, again giving a forward-backward optimality system).

> [!warning] What could be verified
> In every public source for this abstract — the OpenAlex inverted index, Crossref-derived aggregators, and the publisher page as relayed by search — **all inline mathematics is deleted**, producing sentences such as "Solving the linear system ___ is often the major computational burden … where ___ is the so-called all-at-once matrix". The prose is reliable; the symbols cannot be recovered from open sources. In particular the scaling law for $\alpha$ is unverified here, so no relation between $\alpha$, $N_t$ and the error is reported.

## 84: what to do when there is no Toeplitz structure to circulantise

Time spectral methods approximate the solution by a combination of basis functions such as polynomials, a natural companion of spectral discretisation in space and very accurate in time. But all the combination coefficients must be computed **in one shot** by solving an all-at-once system, with no time-stepping option, exactly as for the time-periodic problems of paper 31. The system has Kronecker tensor structure between a time matrix and a space matrix, and the crux is:

**the matrix arising from time spectral methods, denoted $M$, is a non-structured matrix.**

That is the precise obstruction. Every earlier paper in the ParaDiag-II family relies on the temporal matrix being (block) **Toeplitz** — usually block bidiagonal Toeplitz — so the preconditioner follows from the one-line recipe "replace the Toeplitz matrix by its $\alpha$-circulant counterpart" (papers 65 and 71). A time-spectral $M$ is dense and carries no Toeplitz structure, so **there is nothing to circulantise** and the recipe simply does not apply. The paper states that this brings significant challenges for practical computation.

This paper therefore marks the boundary of the whole route: $\alpha$-circulant preconditioning is not a universal trick, its premise is Toeplitz structure in time, and once the time discretisation becomes spectral one must find a different diagonalisable surrogate.

## 85: placing all of this in one framework

Paper 85 is the 2025 _Acta Numerica_ survey, organising this work and the wider literature into two classes: methods that remain effective for propagative problems (Schwarz waveform relaxation, integral deferred correction, ParaExp, ParaDiag) and methods designed primarily for dissipative problems (parareal, PFASST, MGRIT, diagonalised parareal, space-time multigrid). The section-by-section reading on this site, with all 48 original graphical assets and reproducible Python experiments, is a separate topic: [[en/computational-mathematics/knowledge-notes/time-parallelization/index|Time Parallelization for Hyperbolic and Parabolic Problems]].

## How the five relate

| No. | Obstruction addressed                                                           | Means                                                                                      |
| --- | ------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| 59  | $B$ not diagonalisable at uniform steps, $V$ ill conditioned at geometric steps | a boundary value method with $\mathrm{Cond}_2(V)=\mathcal O(n^2)$                          |
| 65  | case-by-case spectra relying on ad hoc structure                                | assume only stability, obtain a modulus bound in $[\frac{1}{1+\alpha},\frac{1}{1-\alpha}]$ |
| 71  | a forward-backward system is not one Toeplitz system                            | circulantise several Toeplitz blocks                                                       |
| 84  | a time-spectral matrix has no Toeplitz structure                                | identify and confront the premise of the recipe                                            |
| 85  | a scattered literature                                                          | organise it by the type of dynamics                                                        |

One judgement runs through the thread: **every advance here takes the form of replacing one assumption with a weaker or more classical one.** Paper 59 replaces "distinct step sizes" with "uniform steps and a different scheme"; paper 65 replaces "Toeplitz plus diagonal dominance" with "stable"; paper 84 identifies when the one remaining assumption, Toeplitz structure, stops holding.

## Coverage check

| Item                                                 | Paper | Status                                                                                   |
| ---------------------------------------------------- | ----- | ---------------------------------------------------------------------------------------- |
| Three-step factorisation and why $B$ is defective    | 59    | factorisation, Jordan block, defective Toeplitz                                          |
| Geometric-step dilemma and the $n\approx20$ ceiling  | 59    | roundoff formula, both limits on $\tau$, parallelism cap                                 |
| Boundary value method and uniform second order       | 59    | hybrid scheme, why not time stepping, $B$ and its conditioning                           |
| One two-matrix notation for both settings            | 65    | $r_1,r_2$, $\mathcal K$ and $\mathcal P_\alpha$, $\alpha$ in two corners                 |
| Main theorem and the exchange of hypotheses          | 65    | stability condition, modulus bound, orientation check, applying $\mathcal P_\alpha^{-1}$ |
| Multiple-block circulantisation for forward-backward | 71    | obstruction, plural recipe, relation to matching Schur complement (limited verification) |
| Non-structured time-spectral matrix                  | 84    | one-shot solve, Kronecker structure, exact reason the recipe fails                       |

## Sources for this page

- J. Liu, X.-S. Wang, S.-L. Wu, and T. Zhou, [_A well-conditioned direct PinT algorithm for first- and second-order evolutionary equations_](https://doi.org/10.1007/s10444-022-09928-4), Adv. Comput. Math. 48 (2022), 16 (preprint [arXiv:2108.01716](https://arxiv.org/abs/2108.01716)).
- S.-L. Wu, T. Zhou, and Z. Zhou, [_A uniform spectral analysis for a preconditioned all-at-once system from first-order and second-order evolutionary problems_](https://doi.org/10.1137/21M145358X), SIAM J. Matrix Anal. Appl. 43(3) (2022), pp. 1331-1353.
- S.-L. Wu, Z. Wang, and T. Zhou, [_PinT preconditioner for forward-backward evolutionary equations_](https://doi.org/10.1137/22M1516476), SIAM J. Matrix Anal. Appl. 44(4) (2023), pp. 1771-1798.
- C. Tang, S.-L. Wu, T. Zhou, and Y. Zhou, [_Parallel-in-time preconditioner for the time spectral methods_](https://doi.org/10.1007/s10915-025-02899-w), J. Sci. Comput. 103 (2025), 82.
- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time parallelization for hyperbolic and parabolic problems_](https://doi.org/10.1017/S0962492924000072), Acta Numer. 34 (2025), pp. 385-489.
