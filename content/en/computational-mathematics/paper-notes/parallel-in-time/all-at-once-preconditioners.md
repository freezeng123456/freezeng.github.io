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

## 59: change the time discretisation so that $V$ is well conditioned by construction

### The idea

Time stepping is sequential because level $n$ depends on level $n-1$. Stack all the time levels into one long vector and that dependence becomes a block lower-triangular matrix, on which time stepping is exactly block forward substitution — which cannot be parallelised, being a chain by construction.

Diagonalisation changes **the coordinates in which the chain is viewed**. If the temporal matrix factors as $B=VDV^{-1}$, then in the coordinates given by $V^{-1}$ the couplings between temporal degrees of freedom come apart completely: "level $n$ waits for level $n-1$" becomes "**mode** $n$ is solved independently". Each mode leaves one complex scalar $\lambda_n$ and one complex-shifted spatial problem $(\lambda_nI_x+\cdots)x=y$, and those $N_t$ problems are unrelated, so they run at once.

The whole price sits in one sentence: $V^{-1}$ and $V$ are applied in floating point, so the computed result carries roundoff pollution of order $O(\epsilon\,\mathrm{Cond}_2(V))$ with $\epsilon$ the machine precision. **A temporal matrix that is diagonalisable in theory but has $\mathrm{Cond}_2(V)=10^{14}$ is useless in double precision.**

That is the bind of the Maday-Rønquist route. The $B$ produced by a standard uniform-step discretisation is not diagonalisable at all (it is a single Jordan block), and their fix is to take pairwise distinct step sizes so that the eigenvalues separate. But the closer the steps are to each other, the closer $B$ is to defective and the larger $\mathrm{Cond}_2(V)$; the more they differ, the worse the discretisation accuracy. Squeezed from both sides, the method is limited to roughly $20$ to $25$ time points.

**This paper turns the design question around.** Earlier work asks "given the time discretisation, how should the step sizes be arranged so that $V$ is usable?"; this paper asks "what time discretisation produces a well-conditioned $V$ **by construction**?" The answer is a **boundary value method**: a scheme that cannot be run as a time-stepping method at all and is meaningful only all-at-once. It is diagonalisable on a uniform grid, and its $\mathrm{Cond}_2(V)$ grows only like $n^2$. The price is giving up the time-stepping interpretation — but since the goal was a one-shot all-at-once solve anyway, that price is free.

### Setting

For $u'(t)+Au(t)=g(t)$ one forms the all-at-once system $\mathcal K\boldsymbol u:=(B\otimes I_x+I_t\otimes A)\boldsymbol u=\boldsymbol b$, and if $B=VDV^{-1}$ then

$$
\mathcal K=(V\otimes I_x)\bigl(D\otimes I_x+I_t\otimes A\bigr)(V^{-1}\otimes I_x),
$$

giving a three-step solve: apply $(V^{-1}\otimes I_x)$ to the right-hand side, perform $n$ fully decoupled spatial solves, apply $(V\otimes I_x)$ to come back. **The obstruction is that the $B$ of a standard uniform-step discretisation is not diagonalisable at all**: for uniform-step backward Euler $B$ is lower bidiagonal Toeplitz, a single Jordan block; for multistep methods $B$ is lower-triangular Toeplitz, likewise defective.

Maday and Rønquist (2008) take **pairwise distinct** step sizes $\{\Delta t_j\}$. The diagonal entries $1/\Delta t_j$ are then distinct, $B$ has $n$ distinct eigenvalues and is diagonalisable, with $D=\mathrm{diag}(1/\Delta t_1,\dots,1/\Delta t_n)$. But the roundoff error obeys

$$
\texttt{roundoff error}=\mathcal O\bigl(\epsilon\,\mathrm{Cond}_2(V)\bigr),
\qquad \epsilon=\text{machine precision}\ (2.22\times10^{-16}\text{ in double}),
$$

and with geometrically increasing steps $\Delta t_j=\Delta t_1\tau^{\,j-1}$ the ratio $\tau>1$ is caught in a vice: $\tau\to1$ makes $B$ nearly defective and $\mathrm{Cond}_2(V)$ explode; $\tau\gg1$ makes the steps grow exponentially and destroys the discretisation accuracy.

This is not merely qualitative — both sides have closed forms. With geometric steps and the linear $\theta$-method, $V$ and $V^{-1}$ are **both unit lower-triangular Toeplitz**,

$$
V=\mathbb T(p_1,\dots,p_{n-1}),
\qquad
V^{-1}=\mathbb T(q_1,\dots,q_{n-1}),
\qquad
\mathbb T(a_1,\dots,a_{n-1}):=
\begin{bmatrix}1\\a_1&1\\\vdots&\ddots&\ddots\\a_{n-1}&\cdots&a_1&1\end{bmatrix},
$$

and for backward Euler ($\theta=1$),

$$
p_j=\frac{1}{\prod_{i=1}^{j}\bigl(1-\tau^{\,i}\bigr)},
\qquad
q_j=(-1)^j\,\tau^{\,j(j-1)/2}\,p_j .
$$

> [!warning] The parameter here is the geometric **ratio**, not the increment
> The survey transcribes this closed form as $p_j=1/\prod_i(1-\varrho^{\,i})$ with $\varrho=\tau-1$ the **increment**. That cannot be right: written with the increment, $\varrho\to0$ would send $p_j\to1$ and make $V$ well conditioned, which is precisely the limit everyone agrees is catastrophic. The $(2,1)$ entry of $BV=VD$ gives $p_1=1/(1-\tau)$ directly, confirming the ratio.

Quantifying both errors lets one solve for the optimal stretching. Write $\varrho=\tau-1$. For $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$ with $\sigma(A)\subset\mathbb R^-$ and $|\lambda(A)|\le\lambda_{\max}$, writing $\boldsymbol u_{N_t}(\varrho)$ for backward Euler on the geometric grid, $\boldsymbol u_{N_t}(0)$ for the uniform grid, and $\tilde{\boldsymbol u}_n(\varrho)$ for what the three-step diagonalisation actually computes in floating point,

$$
\underbrace{\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u_{N_t}(0)\|\lesssim C(\lambda_*T,N_t)\,\varrho^{2}}_{\text{truncation penalty of stretching, grows with }\varrho},
\qquad
\underbrace{\|\tilde{\boldsymbol u}_n(\varrho)-\boldsymbol u_n(\varrho)\|\lesssim\epsilon\,\frac{N_t^2(2N_t+1)(N_t+\lambda_{\max}T)}{\phi(N_t)}\,\varrho^{-(N_t-1)}}_{\text{roundoff of the diagonalisation, blows up as }\varrho\to0},
$$

with $C(x,N_t):=\tfrac{N_t(N_t^2-1)}{24}r(x/N_t,N_t)$, $r(\tilde x,N_t):=\bigl(\tfrac{\tilde x}{1+\tilde x}\bigr)^2(1+\tilde x)^{-N_t}$, $\lambda_*:=N_t\tilde x_*/T$ where $\tilde x_*$ maximises $r(\cdot,N_t)$ on $[0,\infty)$, and

$$
\phi(N_t):=
\begin{cases}
\bigl(\tfrac{N_t}{2}\bigr)!\bigl(\tfrac{N_t}{2}-1\bigr)!, & N_t\ \text{even},\\[4pt]
\Bigl(\bigl(\tfrac{N_t-1}{2}\bigr)!\Bigr)^2, & N_t\ \text{odd}.
\end{cases}
$$

Balancing the two gives

$$
\varrho_{\rm opt}=\left(\epsilon\,\frac{N_t^{2}(2N_t+1)(N_t+\lambda_{\max}T)}{\phi(N_t)\,C(\lambda_*T,N_t)}\right)^{\frac{1}{N_t+1}} .
$$

**The crux is that $\varrho^{-(N_t-1)}$: the roundoff bound degrades exponentially in the number of time steps.** The Maday-Rønquist route is therefore intrinsically a **short-window** method, practical for about $N_t\approx20$ to $30$, with long horizons handled window by window. Removing that cap is paper 59's objective.

### Derivation

**Step one: adopt a scheme that only makes sense all-at-once.** The paper's time discretisation is hybrid: a centred (leapfrog) difference for the first $n-1$ steps and implicit Euler for the last step only,

$$
\begin{cases}
\dfrac{u_{j+1}-u_{j-1}}{2\Delta t}+Au_j=g_j, & j=1,2,\dots,n-1,\\[6pt]
\dfrac{u_n-u_{n-1}}{\Delta t}+Au_n=g_n. &
\end{cases}
$$

This is a **boundary value method**: it must **not** be run as a time-stepping scheme, since the centred scheme used that way is unstable; it is meaningful only when solved all-at-once. The construction is due to Axelsson and Verwer (1985), who proved **uniform second-order accuracy** of the simultaneously obtained solutions even though the last step is only first order; earlier, Fox (1954) and Fox-Mitchell (1957) used a BDF2 last step. Brugnano, Mazzia and Trigiante (1993) solved the resulting all-at-once system iteratively; this paper solves it **directly** by diagonalisation.

The decisive gain is that the steps can be **uniform**, so there is no step-ratio parameter at all:

$$
B=\frac{1}{\Delta t}
\begin{bmatrix}
0&\tfrac12&&&\\
-\tfrac12&0&\tfrac12&&\\
&\ddots&\ddots&\ddots&\\
&&-\tfrac12&0&\tfrac12\\
&&&-1&1
\end{bmatrix},
\qquad
\boldsymbol b=\begin{bmatrix}\tfrac{u_0}{2\Delta t}+g_1\\ g_2\\ \vdots\\ g_n\end{bmatrix}.
$$

Only $u_0$ enters the right-hand side and all time steps are obtained in one shot. Note that $B$ is skew-symmetric Toeplitz apart from its last row — **it is precisely that last row which breaks the symmetry, and precisely that last row which makes the spectrum computable in closed form**. For the second-order (wave-type) equation $u''+Au=g$ the same discretisation gives $\bigl(B^2\otimes I_x+I_t\otimes A\bigr)\boldsymbol u=\boldsymbol b$, and $B^2$ shares the **same** $V$ with $B$, so the same conditioning analysis applies verbatim. That is why one paper covers first- and second-order problems together.

**Step two: reduce the eigenproblem to Chebyshev polynomials.** The conclusion is easy to remember; the mechanism is the substance: **the eigenproblem reduces to Chebyshev polynomials, so both $V$ and $V^{-1}$ have explicit formulas.** Work with the rescaled $\mathcal B:=\Delta t\,B$ and write the Chebyshev polynomials of the first and second kind as

$$
T_n(x)=\cos(n\arccos x),
\qquad
U_n(x)=\frac{\sin\bigl[(n+1)\arccos x\bigr]}{\sin(\arccos x)} .
$$

Then the $n$ eigenvalues of $\mathcal B$ are $\lambda_j=\mathrm ix_j$, where $\{x_j\}_{j=1}^n$ are the $n$ roots of the **characteristic equation**

$$
U_{n-1}(x)-\mathrm i\,T_n(x)=0,
$$

and the components of the corresponding eigenvector are

$$
p_{j,k}=\mathrm i^{\,k}\,U_k(x_j),\qquad k=0,1,\dots,n-1,\qquad p_{j,0}=1 .
$$

**Step three: establish the root structure.** The roots supply exactly the three properties needed: all of them are **simple** (so $\mathcal B$ is diagonalisable), all are complex with negative imaginary part (so the eigenvalues lie in the correct half-plane), and all have modulus below $1+1/\sqrt{2n}$ (this is what controls $\mathrm{Cond}_2(V)$); moreover if $x$ is a root then so is $-\bar x$. The proof substitutes

$$
y=x+\sqrt{x^2-1}
\qquad\Longrightarrow\qquad
T_n=\tfrac12\bigl(y^n+y^{-n}\bigr),
\qquad
U_{n-1}=\frac{y^n-y^{-n}}{y-y^{-1}},
$$

so that the characteristic equation $U_{n-1}-\mathrm iT_n=0$ becomes

$$
y^{2n}=-\frac{(y-\mathrm i)^2}{(y+\mathrm i)^2},
$$

and then uses $|y|>1$ together with the Chebyshev Pythagorean identity $T_n^2(x)+(1-x^2)U_{n-1}^2(x)=1$ to bound the modulus. **This is the technical core of the paper**: the $y$-substitution turns a transcendental equation in two kinds of Chebyshev polynomial into a purely algebraic equation of degree $2n$ whose right-hand side is the square of a Möbius-type quantity of modulus one, so $|y|^{2n}$ is pinned down and with it the location of the roots.

**Step four: reduce the conditioning to a Vandermonde-like matrix.** The eigenvector matrix factors into a **unitary** diagonal term times a Chebyshev-Vandermonde matrix:

$$
V=\underbrace{\mathrm{diag}\bigl(\mathrm i^0,\mathrm i^1,\dots,\mathrm i^{\,n-1}\bigr)}_{=:\,\Theta,\ \text{unitary}}\;
\underbrace{\begin{bmatrix}
U_0(x_1)&\cdots&U_0(x_n)\\
\vdots&&\vdots\\
U_{n-1}(x_1)&\cdots&U_{n-1}(x_n)
\end{bmatrix}}_{=:\,\Phi},
$$

so $\mathrm{Cond}_2(V)=\mathrm{Cond}_2(\Theta\Phi)=\mathrm{Cond}_2(\Phi)$ — the unitary factor is free and the conditioning question is entirely reduced to $\Phi$. The final step uses the **Christoffel-Darboux formula** and properties of the relevant orthogonal polynomials to estimate $\|\Phi\|_2$ and $\|\Phi^{-1}\|_2$, yielding $\mathrm{Cond}_2(\Phi)=\mathcal O(n^2)$.

### Theorems

**(Main theorem: polynomial growth of the condition number.)** Let $B$ come from the uniform-step boundary value method above and let $V$ be the eigenvector matrix of $\mathcal B=\Delta tB$, normalised as above with $p_{j,0}=1$. Then

$$
\mathrm{Cond}_2(V)=\mathcal O(n^2) .
$$

The bound holds for $n\ge8$. Combined with $\texttt{roundoff error}=\mathcal O(\epsilon\,\mathrm{Cond}_2(V))$, this means the roundoff error grows only **polynomially** in the number of time points, in contrast with the exponential blow-up of the geometric-step route. The immediate corollary is that, compared with other direct parallel-in-time algorithms, a **much larger $n$** may be used to obtain satisfactory parallelism.

**(The bound is conservative.)** $\mathcal O(n^2)$ is a **conservative** upper bound: numerically one observes $\mathrm{Cond}_2(V)=\mathcal O(n^{1.75})$, which the authors could not prove.

**(A fast spectral-decomposition algorithm.)** The paper also designs an $\mathcal O(n^2)$ structure-exploiting algorithm for computing the spectral decomposition of $\mathcal B$, in particular $V^{-1}=\Phi^{-1}\Theta^{*}$, based on the three-term recurrence $2yU_j(y)=U_{j+1}(y)+U_{j-1}(y)$ of the second-kind Chebyshev polynomials. It is reported to be much faster than MATLAB's `eig`.

**(Second-order problems come for free.)** The second-order all-at-once matrix is $B^2\otimes I_x+I_t\otimes A$, and $B^2$ shares $V$ with $B$, so no new analysis is required.

Placing the two ParaDiag-I routes side by side makes the change visible:

| Route                            | Grid        | Source of diagonalisability          | Growth of $\mathrm{Cond}_2(V)$          | Practical $n$    |
| -------------------------------- | ----------- | ------------------------------------ | --------------------------------------- | ---------------- |
| Maday-Rønquist / geometric steps | variable    | distinct steps separate the spectrum | exponential blow-up as $\tau\to1$       | about $20$–$25$  |
| Paper 59 / boundary value method | **uniform** | the scheme itself gives simple roots | $\mathcal O(n^2)$ (observed $n^{1.75}$) | no intrinsic cap |

### Numerical experiments

The paper reports numerical results on a parallel machine with **over 60 times speedup achieved on 256 cores**; §4.1 compares against the geometric-step direct algorithm of Gander et al. and documents that method's limitation to about $n\approx20$ to $25$ time points.

A parallel efficiency of $60/256\approx23\%$ is a reasonable order of magnitude for a direct method, and the loss has two structural sources: steps (a) and (c) are global transforms in the time direction and hence communication-heavy, while the $N_t$ complex-shifted spatial systems in step (b) all carry different shifts, so an iterative spatial solver gives different iteration counts on each and a natural load imbalance.

The cap that paper 59 sets out to remove has a complete experimental characterisation in the survey (paper 85; see the [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i|ParaDiag-I chapter]]), which serves as the quantitative background to this paper's motivation:

| Item                      | Setting                                                                         |
| ------------------------- | ------------------------------------------------------------------------------- |
| Equations                 | 1D heat equation and advection-diffusion (viscosity $10^{-2}$)                  |
| Boundary and initial data | homogeneous Dirichlet, $u_0(x)=\sin(2\pi x)$                                    |
| Spatial mesh              | $\Delta x=1/50$                                                                 |
| Time window               | $T=0.2$ (sweeping $\varrho$); $T=0.5$ (sweeping $N_t$)                          |
| Sweep range               | $\varrho\in[10^{-2},1]$ at five values of $N_t$; and $N_t=2^4,2^5,\dots,2^{10}$ |
| Error measure             | maximum $L^\infty$ error over all time nodes                                    |

| Quantity observed                                       | Result                                                                                   |
| ------------------------------------------------------- | ---------------------------------------------------------------------------------------- |
| Error as a function of $\varrho$                        | every curve has a unique minimum                                                         |
| Predictive power of $\varrho_{\rm opt}$                 | very accurate for advection-diffusion; slightly off for the heat equation at small $N_t$ |
| Error versus $N_t$ at the numerically optimal $\varrho$ | falls, then crosses a threshold below $100$ steps and rises sharply                      |
| Control: uniform-step backward Euler                    | error decreases monotonically in $N_t$                                                   |

The last two rows are the point of the experiment: **the cap is a roundoff cap, not a convergence cap.** The same time discretisation on a uniform grid keeps improving as $N_t$ grows; only after geometric stretching plus diagonalisation does the error turn around near $100$ steps, so the loss is entirely attributable to the diagonalisation itself. This is exactly the mechanism paper 59 removes by changing the time discretisation.

### Relation to the others

This paper is the decisive repair of the **ParaDiag-I** (direct, Maday-Rønquist) branch. [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|Papers 31 and 46]] exploit the diagonalisation that time-**periodic** problems supply for free ($V$ unitary, $\mathrm{Cond}_2(V)=1$); paper 59 obtains $\mathcal O(n^2)$ for genuine **initial-value** problems, by changing the **time integrator** rather than the step sizes.

The contrast with the **ParaDiag-II** branch (papers 39, 53, 65, 71, 84) is the main axis of this thread: there one keeps uniform steps and a standard integrator and instead uses an $\alpha$-circulant matrix as a **preconditioner**, paying $\mathrm{Cond}_2(V)\le1/\alpha$ and iterating. Paper 59's selling point is that it is **direct** — no iteration, no $\alpha$ — with only polynomial growth of the conditioning.

Paper 65 shares the "first- and second-order evolutionary problems in one framework" ambition and an overlapping author set: paper 65 does it on the preconditioned iterative branch, paper 59 on the direct branch. Finally, analysing a Vandermonde-like eigenvector matrix through Chebyshev polynomials and Christoffel-Darboux has no second instance in this group; it is the most classically "special functions" paper of the set.

## 65: replace every structural hypothesis by one classical property, stability

### The idea

By 2022 block $\alpha$-circulant preconditioning had become a leading parallel-in-time method, especially for hyperbolic problems, because applying $\mathcal P_\alpha^{-1}$ reduces via FFT to independent solves at all time levels. But its theory had accumulated **case by case**: each paper picks a time integrator and, exploiting that integrator's particular structure, works out the spectrum of $\mathcal P_\alpha^{-1}\mathcal K$. The survey is blunt about this: the analyses are intricate and rely heavily on special properties of the time-stepping matrix such as **sparsity, Toeplitz structure and diagonal dominance**.

What was missing was one theorem covering all integrators, with a hypothesis that is a **classical** property rather than an ad hoc structural one. The property this paper settles on is **stability**.

Why is stability exactly the right hypothesis? Look first at what the preconditioner actually changes. $\mathcal P_\alpha$ differs from the all-at-once matrix $\mathcal K$ in a single top-right entry $\alpha$: it wraps the last time level back onto the first with weight $\alpha$. So the error committed by preconditioning is precisely **one trip of the propagator around the whole time window**. If the method is stable, one complete trip cannot amplify — in scalar language $|r|^{N_t}\le1$ — so the size of that wrap-around perturbation is capped by $\alpha$ itself, **independently of the number of steps, the mesh sizes and the stiffness of $A$**. The derivation below writes this as $\lambda=1/(1-\alpha r^{N_t})$, and the mesh-independence is the direct consequence of $r^{N_t}$ being controlled by $|r|\le1$.

The second gap is memory. Second-order (wave-type) problems are usually handled by rewriting $u''=Au+g$ as a first-order system, doubling the memory per time level — bad for fine spatial meshes and high dimensions — so a direct analysis of **two-step** discretisations of the second-order form was wanted.

### Setting

In the first-order case a general one-step method is written in two-matrix form,

$$
r_1(\Delta t A)\,\boldsymbol u_n=r_2(\Delta t A)\,\boldsymbol u_{n-1}+\tilde{\boldsymbol g}_n,
\qquad n=1,\dots,N_t,
$$

so the method's stability function is $r_1^{-1}(z)r_2(z)$ (backward Euler, for instance, takes $r_1(z)=I-z$ and $r_2(z)=I$). Stacking gives the all-at-once matrix and its block $\alpha$-circulant preconditioner,

$$
\mathcal K=I_t\otimes r_1(\Delta tA)-B\otimes r_2(\Delta tA),
\qquad
\mathcal P_\alpha=I_t\otimes r_1(\Delta tA)-C_\alpha\otimes r_2(\Delta tA),
$$

written out as

$$
\mathcal P_\alpha=
\begin{bmatrix}
r_1(\Delta tA)&&&-\alpha\,r_2(\Delta tA)\\
-r_2(\Delta tA)&r_1(\Delta tA)&&\\
&\ddots&\ddots&\\
&&-r_2(\Delta tA)&r_1(\Delta tA)
\end{bmatrix},
\qquad
C_\alpha=
\begin{bmatrix}0&&&\alpha\\1&0&&\\&\ddots&\ddots&\\&&1&0\end{bmatrix},
$$

where $B$ is the strict down-shift and $C_\alpha=B+\alpha e_1e_{N_t}^\top$. In other words, **the preconditioner is obtained from the all-at-once matrix by putting $\alpha$ in the top-right corner**, i.e. by replacing the (block) Toeplitz matrix with an $\alpha$-circulant one. Setting $\alpha=1$ gives an ordinary circulant (Strang type); $\alpha\to0$ recovers $\mathcal K$ itself.

The second-order case avoids first-order reduction and uses a symmetric two-step method,

$$
r_1(\Delta t^2A)\,\boldsymbol u_{n+1}-r_2(\Delta t^2A)\,\boldsymbol u_n
+r_1(\Delta t^2A)\,\boldsymbol u_{n-1}=\tilde{\boldsymbol g}_n,
\qquad n=1,\dots,N_t-1,
$$

for example the Numerov-type method with $r_1=I_x-\tfrac{\Delta t^2A}{12}+\tfrac{10\gamma(\Delta t^2A)^2}{12}$ and $r_2=2I_x+\tfrac{10\Delta t^2A}{12}+\tfrac{20\gamma(\Delta t^2A)^2}{12}$. The corresponding all-at-once matrix and preconditioner are

$$
\mathcal K=\tilde B\otimes r_1(\Delta t^2A)-B\otimes r_2(\Delta t^2A),
\qquad
\mathcal P_\alpha=\tilde C_\alpha\otimes r_1(\Delta t^2A)-C_\alpha\otimes r_2(\Delta t^2A),
$$

$$
\tilde B=\begin{bmatrix}1&&&&\\0&1&&&\\1&0&1&&\\&\ddots&\ddots&\ddots&\\&&1&0&1\end{bmatrix},
\qquad
\tilde C_\alpha=\begin{bmatrix}1&&&\alpha&\\0&1&&&\alpha\\1&0&1&&\\&\ddots&\ddots&\ddots&\\&&1&0&1\end{bmatrix},
$$

where the two-step stencil puts $\alpha$ in **two** corner positions.

$\mathcal P_\alpha^{-1}$ is applied exactly as elsewhere in the $\alpha$-circulant family: $C_\alpha=VDV^{-1}$ with $V=\Gamma_\alpha^{-1}\mathbb F^*$ and $\Gamma_\alpha=\mathrm{diag}(1,\alpha^{1/N_t},\dots,\alpha^{(N_t-1)/N_t})$, so steps (a) and (c) are (scaled) FFTs, step (b) is $N_t$ independent complex spatial solves, and $\mathrm{Cond}_2(V)\le1/\alpha$ is the roundoff price of small $\alpha$. That machinery is derived in full on the [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|diagonalisation technique]] page.

### Derivation

The orientation can be settled by direct computation rather than taken on trust. Take one $z\in\sigma(\Delta tA)$ at a time and write $r:=r_1^{-1}(z)r_2(z)$, so stability means $|r|\le1$. In that scalar channel $\mathcal K=I_t-rB$ and $\mathcal P_\alpha=I_t-rC_\alpha$ with $C_\alpha=B+\alpha e_1e_{N_t}^{\top}$, hence

$$
\mathcal K=\mathcal P_\alpha+\alpha r\,e_1e_{N_t}^{\top}
\qquad\Longrightarrow\qquad
\mathcal P_\alpha^{-1}\mathcal K=I_t+\alpha r\,\bigl(\mathcal P_\alpha^{-1}e_1\bigr)e_{N_t}^{\top},
$$

a **rank-one update of the identity**. Its eigenvalues are therefore $1$ with multiplicity $N_t-1$, together with the single value $1+\alpha r\,e_{N_t}^{\top}\mathcal P_\alpha^{-1}e_1$. Expanding $\mathcal P_\alpha^{-1}=\sum_{j\ge0}r^jC_\alpha^{\,j}$ and using $C_\alpha^{\,j}e_1=e_{1+j}$ for $j<N_t$ together with $C_\alpha^{N_t}=\alpha I_t$,

$$
e_{N_t}^{\top}\mathcal P_\alpha^{-1}e_1=\sum_{m\ge0}r^{\,N_t-1+mN_t}\alpha^m=\frac{r^{\,N_t-1}}{1-\alpha r^{\,N_t}},
\qquad
\lambda=1+\frac{\alpha r^{\,N_t}}{1-\alpha r^{\,N_t}}=\frac{1}{1-\alpha r^{\,N_t}} .
$$

Since $|r|\le1$ gives $|\alpha r^{N_t}|\le\alpha$ and hence $|1-\alpha r^{N_t}|\in[1-\alpha,1+\alpha]$, this yields exactly $\frac1{1+\alpha}\le|\lambda|\le\frac1{1-\alpha}$; and the same computation gives the survey's own corollary $\rho(\mathcal M)=\bigl|\alpha r^{N_t}/(1-\alpha r^{N_t})\bigr|\le\alpha/(1-\alpha)$. Both conclusions come from one step, which fixes the orientation.

The computation also settles two side facts. First, in each scalar channel all but **one** eigenvalue equal $1$ exactly, so the full system has at most $N_x$ eigenvalues away from $1$ — the exact analogue of McDonald, Pestana and Wathen's clustering theorem at $\alpha=1$. Second, the deviation carries the factor $r^{N_t}$, so channels that are **strictly** contractive ($|r|<1$, i.e. dissipative) cluster far more tightly than the worst-case bound. **This is why the method performs markedly better on parabolic than on hyperbolic problems**: in the hyperbolic case $|r|\approx1$, $r^{N_t}$ does not decay, and the worst case is nearly attained.

The argument runs **channel by channel**: it assumes $A$ is diagonalisable so that each $z\in\sigma(\Delta tA)$ can be treated separately. The paper's own theorem may take a sharper, $z$-dependent form at the endpoints.

### Theorems

**(Main theorem, first-order case.)** Let $\mathcal K$ come from the one-step integrator above and let $\mathcal P_\alpha$ be the block $\alpha$-circulant preconditioner with $\alpha\in(0,1)$. If the integrator is **stable**, that is

$$
\bigl|r_1^{-1}(z)\,r_2(z)\bigr|\le1
\qquad \forall z\in\sigma(\Delta t A)\subset\mathbb C^-,
$$

then every eigenvalue of the preconditioned matrix satisfies the **mesh-independent** bound

$$
\frac{1}{1+\alpha}\ \le\ \bigl|\lambda(\mathcal P_\alpha^{-1}\mathcal K)\bigr|\ \le\ \frac{1}{1-\alpha} .
$$

> [!warning] Orientation of the bound
> The survey transcribes this theorem with the two endpoints transposed (as $\frac{1}{1-\alpha}\le|\lambda|\le\frac{1}{1+\alpha}$), which is impossible for $\alpha\in(0,1)$ since $1/(1-\alpha)>1/(1+\alpha)$. The display above has the correct orientation, for the reason given by the direct computation above.

**(Main theorem, second-order case.)** For the symmetric two-step method above, if

$$
\bigl|r_1^{-1}(z)r_2(z)\bigr|\le2
\qquad\forall z\in\sigma(\Delta t^2A)\subset\mathbb R^-,\ \text{with equality only at }z=0,
$$

then the same two bounds hold. The "$\le2$" is the natural stability condition for a symmetric two-step scheme, matching the coefficient of the middle term $-r_2\boldsymbol u_n$.

**(Corollary: contraction of the stationary iteration.)** For the iteration matrix $\mathcal M=\mathcal I-\mathcal P_\alpha^{-1}\mathcal K$,

$$
\rho(\mathcal M)\le\frac{\alpha}{1-\alpha},
$$

an $\mathcal O(\alpha)$ contraction that is **independent of $\Delta t$, $\Delta x$, $N_t$, $T$ and $\sigma(A)$**. This is exactly the rate conjectured in [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|paper 53]], now proved for all stable one-step methods.

Writing the three statements as numbers makes the choice of $\alpha$ obvious:

| $\alpha$  | $1/(1+\alpha)$ | $1/(1-\alpha)$ | $\rho(\mathcal M)\le\alpha/(1-\alpha)$ | $\mathrm{Cond}_2(V)\le1/\alpha$ |
| --------- | -------------- | -------------- | -------------------------------------- | ------------------------------- |
| $10^{-1}$ | $0.9091$       | $1.1111$       | $0.1111$                               | $10$                            |
| $10^{-2}$ | $0.9901$       | $1.0101$       | $0.0101$                               | $10^{2}$                        |
| $10^{-3}$ | $0.9990$       | $1.0010$       | $0.0010$                               | $10^{3}$                        |

(The first four columns are computed directly from the theorem and its corollary; the last is the generic bound for $\alpha$-circulant diagonalisation.) Adding the second fact from the derivation — that only **one** eigenvalue per scalar channel deviates from $1$ — shows how friendly this theorem is to Krylov methods: across the full $N_xN_t$-dimensional system at most $N_x$ eigenvalues differ from $1$, and all the rest sit exactly at $1$.

**(Sharpness of the hypothesis.)** Stability is **sufficient**; the survey reports that numerically it is also **necessary**, illustrating with a Numerov-type method at $\gamma=1/120$ (unconditionally stable, fourth order) versus $\gamma=1/120.01$ (unstable), where $\sigma(\mathcal M)$ escapes the circle of radius $\alpha/(1-\alpha)$.

**(Scope and its boundary.)** For first-order problems the analysis works for **all stable single-step time-integrators**, and for second-order problems for **a large class of symmetric two-step methods which could be arbitrarily high-order**. But for **general multistep** methods, beyond one-step and symmetric two-step, the bound need not hold: for instance the $B$ arising from Volterra partial integro-differential equations is a **dense** lower-triangular Toeplitz matrix, for which only $|\lambda(\mathcal P_\alpha^{-1}\mathcal K)|=1+\mathcal O(\alpha)$ is obtained, and only under positivity or monotonicity conditions on the quadrature weights. **That boundary is exactly where paper 84 begins.**

### Numerical experiments

The survey cites this paper's experiments in two places, both concerning situations the theory does **not** cover:

| Experiment cited                         | Behaviour observed                                                                                                                   |
| ---------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| Nonlinear problems, Newton + GMRES       | clustering of $\mathcal P_\alpha^{-1}\mathcal J$; $\rho>1$ can occur, so the stationary iteration fails while GMRES still works well |
| Effect of the time window $T$            | a shorter $T$ gives more clustered eigenvalues and faster GMRES convergence                                                          |
| Numerov parameter (the survey's own run) | $\gamma=1/120$ stays inside the circle of radius $\alpha/(1-\alpha)$; $\gamma=1/120.01$ escapes                                      |

The first row separates the theorem from practice: what the theorem guarantees is $\rho(\mathcal M)\le\alpha/(1-\alpha)$, a statement about the **stationary iteration**. Replace $\mathcal K$ by the Jacobian $\mathcal J$ of a nonlinear problem and the hypothesis no longer holds, so $\rho$ may exceed $1$; but GMRES needs the spectrum to be **clustered**, not $\rho<1$, so it still converges quickly. This is the concrete reason for choosing a Krylov method over the stationary iteration.

The second row needs care, because it appears to run against the linear theory. The linear bound $\frac1{1+\alpha}\le|\lambda|\le\frac1{1-\alpha}$ contains no $T$ at all, and Gander and Wu's refinement for backward Euler — written in the sign convention $\boldsymbol u'+A\boldsymbol u=\boldsymbol g$, opposite to this section's, with $\Re\lambda(A)\ge\kappa\ge0$ —

$$
\rho\le\frac{\alpha e^{-T\kappa}}{1-\alpha e^{-T\kappa}},
$$

actually **improves as $T$ grows**. Neither bound supports "shorter $T$ clusters better": one contains no $T$ at all and the other runs the other way. So the observation cannot be a phenomenon of the linear picture. The reading compatible with both bounds is the nonlinear one: $\mathcal P_\alpha$ is built from a frozen Jacobian, and a shorter window keeps that Jacobian closer to constant.

All three rows are qualitative. The paper's contribution is a theorem rather than a new algorithm, so its experiments corroborate rather than measure.

### Relation to the others

Paper 53 stated the conjecture — that the A-stability of an implicit Runge-Kutta method suffices for a robust $\mathcal O(\alpha)$ convergence rate — and verified it for one family (two-stage SDIRK under $\gamma\ge1/4$); paper 65 proves the corresponding statement for the whole class of stable one-step methods. Together they are the **case study and the theorem**.

[[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|Papers 39 and 31]] built the $\alpha$-circulant machinery, and paper 65 is its spectral theory. Paper 46 proved eigenvalue **and** singular-value clustering for two specific integrators in the **forward-backward** (KKT) setting; paper 71 on this page is the forward-backward analogue of paper 65 — the same "one theorem for all stable one-step integrators" ambition applied to a coupled system, at the price of a scaling law tying $\alpha$ to $N_t$. Paper 59 shares the "first- and second-order in one framework" scope but on the direct branch: paper 59 bounds $\mathrm{Cond}_2(V)=\mathcal O(n^2)$ for an exact diagonalisation, paper 65 bounds the spectrum of an approximate-but-cheap preconditioned system. Paper 84 pushes the analysis to time-spectral discretisations, where the temporal blocks are dense and $\mathcal K$ is not block-bidiagonal at all — one of the cases explicitly outside this paper's hypotheses.

## 71: the forward-backward case has several Toeplitz blocks

Everything so far treats **one** direction of evolution. Optimal control and inverse problems instead produce a pair of evolutions running in opposite directions and coupled to each other: the state equation runs forward from $t=0$, the adjoint equation backward from $t=T$, and the two are tied together through the control. Here not even a sequential algorithm exists — neither direction can march on its own, the coupled system must be solved in one shot, and what results is a saddle-point-type or Schur-complement system rather than a single block lower-triangular Toeplitz system. **Parallel-in-time is not an accelerator here; it is the only way to organise the computation.** The paper's two applications are parabolic PDE-constrained optimal control and parabolic source identification, whose optimality systems are of this same type.

The mechanism is that of paper 65: add a wrap-around of weight $\alpha$ in the time direction, turning Toeplitz into $\alpha$-circulant and hence FFT-diagonalisable. What differs is **how many** wrap-arounds are needed: the coupled system contains not one Toeplitz matrix but several (at minimum the forward matrix and its transpose), each of which must be circulantised. The resulting preconditioner is identified by the authors themselves as a **parallel version of the matching Schur complement preconditioner** of Pearson, Stoll and Wathen, so the contribution is not inventing a preconditioner but converting a mature sequential one into a diagonalisable, hence time-parallel, form. The matching Schur complement is symmetric positive definite, which is consistent with this paper using **conjugate gradients** as its outer solver; by contrast paper 46 uses GMRES and BiCGStab and paper 65 uses GMRES, both of which imply non-symmetric preconditioned systems.

The theorem has the same shape as paper 65's — for any one-step stable time-integrator the eigenvalues of the preconditioned matrix lie in a mesh-independent interval — with one added restriction: **$\alpha$ can no longer be fixed freely and must scale weakly downwards with the number of time steps.** That has a direct practical consequence. Paper 65 permits a fixed $\alpha$, so its roundoff floor $\epsilon\,\mathrm{Cond}_2(V)\le\epsilon/\alpha$ does not grow with $N_t$ at all; once $\alpha$ has to shrink with $N_t$, that floor grows again with the length of the time window. The coupling costs you the freedom to fix $\alpha$.

Three contrasts place the paper. [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|Paper 46]] is the direct predecessor: the same class of physical problem and the same strategic move (replace the non-diagonalisable time matrix by a circulant-type surrogate and use it as a preconditioner), but paper 46 proves clustering for only two specific integrators, backward Euler and the trapezoidal rule, where paper 71 proves a mesh-independent interval for all stable one-step integrators. **Paper 71 is to paper 46 what paper 65 is to paper 53: one general theorem replacing a set of case studies.** Paper 65 is the sibling result for the uncoupled forward problem, and the contrast is the sharpest single sentence of the pair: for a single forward evolution, $|\lambda(\mathcal P_\alpha^{-1}\mathcal K)|\in[\tfrac{1}{1+\alpha},\tfrac{1}{1-\alpha}]$ holds for **any fixed** $\alpha\in(0,1)$; for the forward-backward system, $\alpha$ must shrink with $N_t$. Papers 31 and 39 supply the underlying $\alpha$-circulant and FFT machinery, and paper 84 is the latest member of the same programme, moving on to dense time-spectral blocks.

## 84: what to do when there is no Toeplitz structure to circulantise

### The idea

Time spectral methods approximate the solution by a combination of basis functions (polynomials, for instance) and are a natural companion of spectral discretisation in space, giving very high accuracy in time. But all the combination coefficients must be computed **in one shot** by solving an all-at-once system — there is no time-stepping option at all, structurally the same situation as the time-periodic problems of paper 31: the discretisation itself offers no sequential algorithm, so an all-at-once solve is mandatory rather than optional.

The trouble lies in the one-line recipe this whole route runs on. Every earlier paper in the ParaDiag-II family relies on the temporal matrix being (block) **Toeplitz** — usually block bidiagonal Toeplitz — so that the preconditioner follows from a single sentence: "replace the Toeplitz matrix by its $\alpha$-circulant counterpart". A time-spectral $M$ is dense and carries no Toeplitz structure, so **there is nothing to circulantise and the recipe simply fails**.

The paper's answer fits in four words: **factorise first, then circulantise.** The dense, unstructured $M$ is first factored into pieces some of which **are** Toeplitz, and only then is the $\alpha$-circulant substitution applied — to those factors, not to $M$ itself. $M$ need not have structure; it only needs to have structure extractable from it.

### Setting

The all-at-once system is still of Kronecker tensor form,

$$
\bigl(M\otimes I_x+I_t\otimes A\bigr)\boldsymbol u=\boldsymbol b,
\qquad
M\in\mathbb R^{N_t\times N_t}\ \text{dense and unstructured}
$$

Compared with the earlier papers, the only thing that changes is $M$: one-step methods give a block bidiagonal Toeplitz $B$, symmetric two-step methods give the two-diagonal $\tilde B$, and a time-spectral method gives a full matrix.

Why the time matrix becomes dense is clearest in collocation-type methods. Integral deferred correction (IDC) on nodes $t_0<t_1<\dots<t_N$ has an update containing the term

$$
\int_{t_m}^{t_{m+1}}f\bigl(\boldsymbol u^k(\tau),\tau\bigr)\,d\tau
\approx\sum_{j=1}^{N}\omega_{m,j}\,f(\boldsymbol u_j^k,t_j),
\qquad
\omega_{m,j}=\int_{t_m}^{t_{m+1}}\Bigl(\prod_{i\ne j}\frac{\tau-t_i}{t_j-t_i}\Bigr)d\tau,
$$

that is, the integral of the Lagrange interpolant through **all** the nodes. The matrix $\Omega=[\omega_{m,j}]$ is therefore full: the update on sub-interval $m$ depends on **every** node, including nodes to its right, and the lower-triangular structure is destroyed. Time spectral methods are the extreme form of the same phenomenon — the basis functions are supported on the whole interval, so all coefficients are fully coupled from the outset.

### Derivation

**The key construction.** The preconditioner "is obtained by a **novel factorization of $M$** and then a **replacement of the Toeplitz matrices in such a factorization** by the corresponding $\alpha$-circulant matrix of **Strang-type**."

Two terms need unpacking. **Strang-type** refers to the classical Strang circulant preconditioner for Toeplitz systems, in which the circulant is built from the central diagonals of the Toeplitz matrix; the $\alpha$-circulant generalisation replaces the wrap-around $1$s by $\alpha$. That lineage goes back to Strang (1986), and the spectral analysis of $\sigma(C^{-1}B)$ has three decades of literature behind it (Chan-Ng 1996; Ng 2004; Bini-Latouche-Meini 2005). What is **new** in ParaDiag-II is that the blocks $r_1(\Delta tA)$ and $r_2(\Delta tA)$ are themselves **not** Toeplitz, so the classical block-Toeplitz theory does not transfer — which is precisely the gap papers 65, 71 and 84 fill, with paper 84 facing the most extreme version of it: not even the time direction is Toeplitz any more.

**Scope.** The paper exhibits the factorisation concretely for the **Legendre dual-Petrov-Galerkin method** and claims it holds for other widely used time spectral methods as well. The background for that time discretisation is Shen's dual-Petrov-Galerkin method (SINUM 2003), Shen-Wang's Legendre and Chebyshev dual-Petrov-Galerkin methods for hyperbolic equations, Kong-Shen-Wang-Xiang's eigenvalue analysis of Legendre dual-Petrov-Galerkin methods for initial value problems (Adv. Comput. Math. 2024), Tang-Ma's single- and multi-interval Legendre $\tau$-methods, and Yang-Wang's Chebyshev-Gauss spectral collocation method.

**The design criterion.** The preconditioner "permits **well-conditioned diagonalization** and thus each preconditioning step can be solved in an efficient time parallel manner". Note that this criterion is identical to paper 59's, only imposed on the **preconditioner** rather than on a direct solver: paper 59 requires the exact diagonalisation to be well conditioned, paper 84 requires the approximate one to be, and hands the approximation error to GMRES.

### Theorems

The paper establishes three things: (i) the spectral analysis of the preconditioned matrix reveals **highly clustered** eigenvalues, promoting rapid convergence of **GMRES**; (ii) the factorisation is exhibited for the Legendre dual-Petrov-Galerkin method and claimed to hold for other widely used time spectral methods; (iii) the preconditioner admits a well-conditioned diagonalisation, so applying $\mathcal P^{-1}$ is parallel across all time levels.

> [!warning] Paper 65's interval must not be transplanted here
> Paper 65's bound $[\tfrac{1}{1+\alpha},\tfrac{1}{1-\alpha}]$ rests on the hypothesis that $\mathcal K$ is block bidiagonal Toeplitz, and that hypothesis is exactly what fails in the time-spectral setting — which is why this paper exists. "Highly clustered" and an interval with endpoints are not the same kind of statement either: the latter converts directly into a bound on the GMRES iteration count, the former does not.

### Relation to the others

Paper 65 is the theoretical backbone being extended, and is cited. Paper 65 covers all stable **one-step** methods (block bidiagonal Toeplitz $\mathcal K$) and symmetric **two-step** methods; the survey notes explicitly that for general multistep and non-Toeplitz temporal structures its bounds need not hold. A time-spectral $M$ is the extreme case of that gap: fully dense, with no Toeplitz structure at all.

The closest methodological precedent is [[en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|paper 53]], which met the same kind of difficulty one step earlier — a **multi-stage** integrator whose difference equations cannot be stacked directly into a Toeplitz all-at-once system — and solved it by constructing an $\alpha$-circulant preconditioner "with completely different structures and different implementation details". Paper 84 repeats the manoeuvre for a **dense** temporal matrix, but by factorising rather than by re-deriving the all-at-once form. Paper 31's pattern recurs here too: a discretisation that **forces** an all-at-once solve (time-periodicity there, spectral coefficients here) is an opportunity for parallel-in-time rather than an obstacle. Paper 59 shares the "well-conditioned diagonalisation" criterion and is cited, the difference being direct solver (paper 59, $\mathrm{Cond}_2(V)=\mathcal O(n^2)$) versus preconditioned GMRES (paper 84). Paper 71 is the group's immediately preceding paper in the same line and is also cited.

Chronologically this is the most recent paper in the list and reads as the current frontier of the programme: having covered the one-step, two-step, multi-stage and forward-backward cases, the remaining frontier was unstructured, dense temporal operators.

## 85: placing the above in a unified framework

Paper 85 is the 2025 _Acta Numerica_ survey, which organises this work together with the wider literature into two classes: methods that remain effective for **propagation**-type problems (Schwarz waveform relaxation, integral deferred correction, ParaExp, ParaDiag) and methods designed mainly for **dissipative** problems (parareal, PFASST, MGRIT, diagonalisation-based parareal, space-time multigrid). The dichotomy is itself a judgement: the $\alpha$-circulant route of this page falls in the first class because its contraction $\alpha/(1-\alpha)$ does not depend on $\sigma(A)$ and so does not fail on hyperbolic problems the way parareal does — although, as paper 65's derivation shows, in the hyperbolic case it merely **does not fail**, with clustering still markedly weaker than in the parabolic case.

A section-by-section reading of it is a separate topic: [[en/computational-mathematics/knowledge-notes/time-parallelization/index|Time Parallelization for Hyperbolic and Parabolic Problems]].

## How the five relate

| Paper | Obstruction removed                                                              | Means                                                                                        |
| ----- | -------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| 59    | $B$ not diagonalisable on a uniform grid, $V$ ill conditioned on a geometric one | switch to a boundary value method, giving $\mathrm{Cond}_2(V)=\mathcal O(n^2)$               |
| 65    | case-by-case spectral analysis resting on ad hoc structural hypotheses           | assume only stability, obtaining the modulus bound $[\frac{1}{1+\alpha},\frac{1}{1-\alpha}]$ |
| 71    | a forward-backward system is not a single Toeplitz system                        | circulantise the several Toeplitz blocks of the system                                       |
| 84    | a time-spectral matrix has no Toeplitz structure                                 | factor out Toeplitz factors first, then circulantise                                         |
| 85    | a scattered literature                                                           | organise it by the type of dynamics                                                          |

One judgement runs through all of them: **every advance on this route takes the form of replacing one hypothesis by a weaker or more classical one.** Paper 59 replaces "distinct step sizes" by "uniform steps but a different scheme"; paper 65 replaces "Toeplitz plus diagonal dominance" by "stable"; paper 71 replaces "a single forward evolution" by "a forward-backward coupling", at the cost of losing the freedom to fix $\alpha$; paper 84 identifies when the last remaining hypothesis (Toeplitz structure) stops holding, and how to get around it.

Ordering the five by how much roundoff the diagonalisation costs makes the economics of the whole thread visible:

| Route                             | Time matrix                                       | Roundoff amplification of the diagonalisation                    | Dependence on $N_t$ | Representative |
| --------------------------------- | ------------------------------------------------- | ---------------------------------------------------------------- | ------------------- | -------------- |
| ParaDiag-I, geometric steps       | variable-step backward Euler                      | bound carries a factor $\varrho^{-(N_t-1)}$                      | **exponential**     | Maday-Rønquist |
| ParaDiag-I, boundary value method | uniform-step BVM                                  | $\mathcal O(n^2)$ (observed $\mathcal O(n^{1.75})$)              | polynomial          | Paper 59       |
| ParaDiag-II, fixed $\alpha$       | $\alpha$-circulant                                | $\le1/\alpha$                                                    | **none**            | Paper 65       |
| ParaDiag-II, forward-backward     | $\alpha$-circulant, $\alpha$ shrinking with $N_t$ | $\le1/\alpha$, growing with $N_t$                                | grows               | Paper 71       |
| ParaDiag-II, time spectral        | circulantised after factorisation                 | the preconditioner must admit a well-conditioned diagonalisation | —                   | Paper 84       |

The third row is the centre of this table: **fixed-$\alpha$ $\alpha$-circulant preconditioning is the only member of the family whose roundoff cost does not grow with the time window at all**, the price being that it is only an approximation and needs an outer Krylov iteration. The other four rows all explain when that ideal cannot be reached, and how fast it degrades.

## Sources for this page

- J. Liu, X.-S. Wang, S.-L. Wu, and T. Zhou, [_A well-conditioned direct PinT algorithm for first- and second-order evolutionary equations_](https://doi.org/10.1007/s10444-022-09928-4), Adv. Comput. Math. 48 (2022), 16 (preprint [arXiv:2108.01716](https://arxiv.org/abs/2108.01716)).
- S.-L. Wu, T. Zhou, and Z. Zhou, [_A uniform spectral analysis for a preconditioned all-at-once system from first-order and second-order evolutionary problems_](https://doi.org/10.1137/21M145358X), SIAM J. Matrix Anal. Appl. 43(3) (2022), pp. 1331-1353.
- S.-L. Wu, Z. Wang, and T. Zhou, [_PinT preconditioner for forward-backward evolutionary equations_](https://doi.org/10.1137/22M1516476), SIAM J. Matrix Anal. Appl. 44(4) (2023), pp. 1771-1798.
- C. Tang, S.-L. Wu, T. Zhou, and Y. Zhou, [_Parallel-in-time preconditioner for the time spectral methods_](https://doi.org/10.1007/s10915-025-02899-w), J. Sci. Comput. 103 (2025), 82.
- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time parallelization for hyperbolic and parabolic problems_](https://doi.org/10.1017/S0962492924000072), Acta Numer. 34 (2025), pp. 385-489.
