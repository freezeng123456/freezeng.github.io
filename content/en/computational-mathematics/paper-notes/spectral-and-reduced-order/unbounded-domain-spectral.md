---
title: Unbounded Domains and Fractional Operators
description: Papers 27, 42 and 92 - making the far-field behaviour of the basis compatible with the action of a fractional operator
lang: en
translation: computational-mathematics/paper-notes/spectral-and-reduced-order/unbounded-domain-spectral
tags:
  - paper-notes
  - spectral-methods
  - fractional-operators
---

> [!note] Coverage of this page
> Papers **27** (_Commun. Comput. Phys._ 24(4), 2018), **42** (_SIAM J. Sci. Comput._ 42(2), 2020) and **92** (_SIAM J. Sci. Comput._ 48(1), 2026). The full text of paper 92 could not be reached through public channels, so that section reports only what the abstract, keywords and reference list confirm.

The model problem is

$$
(-\Delta)^{\alpha/2}u(x)+\rho\,u(x)=f(x),\quad x\in\mathbb R^d,
\qquad u(x)\to0,\ |x|\to\infty,
$$

with $\rho>0$ and $\alpha\in(0,2)$. The fractional Laplacian has two equivalent definitions, the hypersingular integral

$$
(-\Delta)^{\alpha/2}u(x)=C_{d,\alpha}\,\mathrm{p.v.}\!\int_{\mathbb R^d}
\frac{u(x)-u(y)}{|x-y|^{\,d+\alpha}}\,\mathrm dy,
\qquad
C_{d,\alpha}=\alpha\,2^{\alpha-1}\,
\frac{\Gamma\bigl(\tfrac{d+\alpha}{2}\bigr)}{\pi^{d/2}\,\Gamma\bigl(\tfrac{2-\alpha}{2}\bigr)},
$$

and the pseudodifferential form $(-\Delta)^{\alpha/2}u=\mathcal F^{-1}[|\xi|^{\alpha}\mathcal F[u]]$.

## All three standard unbounded-domain strategies fail here

Paper 42 states the difficulty precisely: solutions of fractional operators decay only **algebraically** at infinity, following a power law. That single fact defeats three standard strategies at once.

Domain truncation suits rapidly decaying solutions but not these; worse, naive truncation introduces nonphysical singularities at the interface where the unbounded domain is terminated. Transparent boundary conditions and sponge layers are highly nontrivial because the operator is nonlocal. Orthogonal functions on unbounded domains — Hermite and Laguerre — are tuned to **exponential** decay and are mismatched to power-law tails.

Mapped Jacobi ("rational") bases were already known to beat Hermite and Laguerre for **integer-order** problems with algebraically decaying solutions. Extending the mapping technique to the fractional case is far from trivial, and for a structural reason: the ordinary Laplacian maps a rational basis function to another function of the same type, while the fractional Laplacian maps it to **a class of functions of a completely different nature**.

## 27: Hermite bases with explicit differentiation matrices

On bounded domains the successful approach is to pick a basis adapted to the fractional operator — Jacobi poly-fractonomials (eigenfunctions of a fractional Sturm-Liouville problem) and generalised Jacobi functions — for which a fractional derivative of a basis function is again a basis function with shifted parameters, making the operator effectively local. No such device exists on an unbounded domain.

The contemporaneous alternative (the Hermite spectral method of Mao and Shen) collocates in frequency space, so every evaluation costs a forward and a backward Hermite transform, which is awkward for nonlinear problems. Paper 27 collocates directly in **physical space** using explicit closed-form differentiation matrices.

The paper supplies two Hermite bases. The normalised Hermite functions

$$
\widehat H_n(x)=\frac{1}{\sqrt{2^n n!}}\,e^{-x^2/2}H_n(x),
\qquad
\int_{\mathbb R}\widehat H_m\widehat H_n\,\mathrm dx=\sqrt{\pi}\,\delta_{mn}
$$

are orthogonal with respect to the weight $\omega\equiv1$; the over-scaled (Brinkman) basis

$$
\widetilde H_n(x)=e^{-x^2/2}\widehat H_n(x)=\frac{1}{\sqrt{2^n n!}}e^{-x^2}H_n(x),
\qquad
\int_{\mathbb R}\widetilde H_m\widetilde H_n\,e^{x^2}\mathrm dx=\sqrt{\pi}\,\delta_{mn}
$$

was introduced by Brinkman for Fokker-Planck equations, where the velocity part of the probability density is expanded in these functions, and has since become one of the standard bases there. The two differ in whether $e^{-x^2/2}$ sits in the basis or in the weight, and that choice determines the form of the discrete inner product and the differentiation matrices.

## 42: decomposing the rational basis into two elementary shapes

### The map and the basis

The one-to-one map and its identities are

$$
x=\frac{t}{\sqrt{1-t^2}}
\ \Longleftrightarrow\
t=\frac{x}{\sqrt{1+x^2}},
\qquad
1-t^2=\frac{1}{1+x^2},
\qquad
\frac{\mathrm dx}{\mathrm dt}=\frac{1}{(1-t^2)^{3/2}} .
$$

This is the **singular** (unbounded-interval) algebraic map, and it is exactly what produces power-law tails in the basis. The modified mapped Gegenbauer functions are

$$
R^{\lambda}_n(x):=\bigl(1+x^2\bigr)^{-\frac{\lambda+1}{2}}
C^{\lambda}_n\!\Bigl(\frac{x}{\sqrt{1+x^2}}\Bigr),
\qquad x\in\mathbb R,\ \lambda>-\tfrac12 .
$$

"Modified" means precisely that the basis **absorbs the square root of the weight function**: writing $R^{\lambda}_n(x)=S(t)C^{\lambda}_n(t)$ with

$$
S(t)=\sqrt{\omega_{\lambda}(t)\frac{\mathrm dt}{\mathrm dx}}=(1-t^2)^{\frac{\lambda+1}{2}} ,
$$

these functions are orthogonal with respect to the **uniform** weight on $\mathbb R$:

$$
\int_{-\infty}^{\infty}R^{\lambda}_n(x)R^{\lambda}_m(x)\,\mathrm dx=\gamma^{\lambda}_n\delta_{nm},
\qquad
\gamma^{\lambda}_n=\frac{\pi\,2^{1-2\lambda}\Gamma(n+2\lambda)}{n!\,(n+\lambda)\,\Gamma^2(\lambda)} .
$$

The basis satisfies a three-term recurrence

$$
n\,R^{\lambda}_n(x)=\frac{2x}{\sqrt{1+x^2}}(n+\lambda-1)R^{\lambda}_{n-1}(x)
-(n+2\lambda-2)R^{\lambda}_{n-2}(x),\qquad n\ge2,
$$

$$
R^{\lambda}_0(x)=\frac{1}{(1+x^2)^{\frac{\lambda+1}{2}}},
\qquad
R^{\lambda}_1(x)=\frac{2\lambda x}{(1+x^2)^{1+\frac{\lambda}{2}}} .
$$

The far-field behaviour is the point:

$$
\lim_{x\to+\infty}(1+x^2)^{\frac{\lambda+1}{2}}R^{\lambda}_n(x)=\frac{(2\lambda)_n}{n!},
$$

so every basis function decays like $|x|^{-(\lambda+1)}$ and $\lambda$ is the knob that tunes the basis to the solution's power law.

### Two elementary shapes

Expanding the basis by the Pfaff transformation shows that **every basis function is a finite linear combination of just two elementary shapes**:

$$
\frac{1}{(1+x^2)^{\gamma}}\ \bigl(\gamma=k+\tfrac{\lambda+1}{2}\bigr),
\qquad
\frac{x}{(1+x^2)^{\gamma}}\ \bigl(\gamma=k+\tfrac{\lambda}{2}+1\bigr).
$$

It therefore suffices to differentiate those two fractionally. For real $s>0$:

$$
(-\Delta)^{s}\Bigl\{\frac{1}{(1+x^2)^{\gamma}}\Bigr\}
=A^{\gamma}_s\ {}_2F_1\!\Bigl(s+\gamma,\ s+\tfrac12;\ \tfrac12;\ -x^2\Bigr),
\qquad \gamma>0,
$$

$$
(-\Delta)^{s}\Bigl\{\frac{x}{(1+x^2)^{\gamma}}\Bigr\}
=(2s+1)A^{\gamma}_s\,x\ {}_2F_1\!\Bigl(s+\gamma,\ s+\tfrac32;\ \tfrac32;\ -x^2\Bigr),
\qquad \gamma>\tfrac12,
$$

$$
A^{\gamma}_s=\frac{2^{2s}\,\Gamma(s+\gamma)\,\Gamma\bigl(s+\tfrac12\bigr)}{\sqrt{\pi}\,\Gamma(\gamma)} .
$$

The proof route is worth recording because it shows where modified Bessel functions enter. Write $v(x)=(1+x^2)^{-\gamma}$ as ${}_2F_1(\gamma,\tfrac12;\tfrac12;-x^2)$; a cosine-transform identity gives a modified Bessel function of the second kind,

$$
\mathcal F[v](\xi)=\frac{2^{1-\gamma}}{\Gamma(\gamma)}\,\xi^{\gamma-\frac12}
K_{\gamma-\frac12}(\xi),\qquad \xi>0,
$$

and then $(-\Delta)^s v=\mathcal F^{-1}[|\xi|^{2s}\mathcal F[v]]$ becomes a Bessel integral finished by the corresponding integral identity. The odd shape follows the same route through $\mathcal F[xv](\xi)=-i\frac{2^{1-\gamma}}{\Gamma(\gamma)}\xi^{\gamma-\frac12}K_{\gamma-\frac32}(\xi)$.

### Where integer and fractional order diverge

This is the most informative structural observation in the paper. For positive integer $s$ both hypergeometric functions terminate and the operator **gains $s$ powers of decay**:

$$
(-\Delta)^{s}\Bigl\{\tfrac{1}{(1+x^2)^{\gamma}}\Bigr\}\sim\frac{1}{(1+x^2)^{s+\gamma}} .
$$

For noninteger $s$ the hypergeometric functions can diverge as $|x|\to\infty$ and the picture changes. For $\gamma>1/2$,

$$
(-\Delta)^{s}\bigl\{(1+x^2)^{-\gamma}\bigr\}\sim(1+x^2)^{-(s+1/2)},
$$

so **the decay rate is independent of $\gamma$**; at $\gamma=1/2$ a logarithm appears, $\ln(1+x^2)/(1+x^2)^{s+1/2}$; and only $0<\gamma<1/2$ retains the integer-order behaviour. At the level of the basis,

$$
(-\Delta)^{s}R^{\lambda}_{2n}(x)\sim
\begin{cases}
(1+x^2)^{-\left(s+\frac{\lambda+1}{2}\right)}, & -\tfrac12<\lambda<0,\\[3pt]
\ln(1+x^2)\,(1+x^2)^{-\left(s+\frac12\right)}, & \lambda=0,\\[3pt]
(1+x^2)^{-\left(s+\frac12\right)}, & \lambda>0 .
\end{cases}
$$

In other words, the fractional Laplacian applied to this basis **does not always gain the factor $1/(1+x^2)^s$**. That is the quantitative version of the "completely different nature" statement in the introduction.

### Two implementation details that matter

The formulas for the fractional Laplacian of the basis are less compact than in the Hermite case, but they evaluate efficiently by recurrence. With $F_k(x)={}_2F_1(a,b;c;-x^2)$, $a=s+k+\tfrac{\lambda+1}{2}$, $b=s+\tfrac12$ and $c=\tfrac12$,

$$
F_{k+1}(x)=\frac{c-a}{a(1+x^2)}F_{k-1}(x)
+\frac{(2a-c)+(a-b)x^2}{a(1+x^2)}F_k(x),\qquad k\ge1 .
$$

A scaling parameter $\mu>0$ improves resolution: the map becomes $x=\mu t/\sqrt{1-t^2}$ and the basis becomes

$$
R^{\lambda}_{n,\mu}(x)=\frac{\mu^{\lambda+\frac12}}{(\mu^2+x^2)^{\frac{\lambda+1}{2}}}
C^{\lambda}_n\!\Bigl(\frac{x}{\sqrt{\mu^2+x^2}}\Bigr)
=\mu^{-\frac12}R^{\lambda}_n\!\Bigl(\frac{x}{\mu}\Bigr).
$$

This is the rational analogue of the Hermite scaling factor in paper 27, with an important difference: here $\mu$ controls the **width** of the algebraic profile while $\lambda$ independently controls the **exponent** of the tail. Separating the two is the practical advantage of the rational basis over the Hermite one.

### Approximation space and two schemes

With $\mathcal V^{\lambda}_N=\mathrm{span}\{R^{\lambda}_n\}=\{S(t)P(t):P\in\mathcal P_N\}$, the change of variables $\breve u(x)=u(x)/s(x)=\breve U(t)$ with $s(x)=(1+x^2)^{-(\lambda+1)/2}=S(t)$ gives $\hat u^{\lambda}_n=\widehat{\breve U}_n$ exactly, so modified mapped Gegenbauer approximation on $\mathbb R$ **is** Gegenbauer approximation on $(-1,1)$. The mapped Gauss nodes and weights are

$$
x^{\lambda}_j=\frac{t^{\lambda}_j}{\sqrt{1-(t^{\lambda}_j)^2}},
\qquad
\omega^{\lambda}_j=\bigl(1+(t^{\lambda}_j)^2\bigr)^{-\lambda}\rho^{\lambda}_j .
$$

The Galerkin scheme uses $\tilde a_s(u,v)=((-\Delta)^{s/2}u,(-\Delta)^{s/2}v)+\rho(u,v)$ over $\mathcal V^{\lambda}_N$ with the source term interpolated, which is why the error bound carries a separate $f$ term. The collocation scheme uses the Lagrange basis on the mapped Gauss grid to build fractional differentiation matrices

$$
\mathcal D_{i,j}=(-\Delta)^{\alpha/2}l_j(x^{\lambda}_i)
=\sum_{k=0}^{N-1}b^j_k\,(-\Delta)^{\alpha/2}R^{\lambda}_k(x^{\lambda}_i),
\qquad
b^j_k=\frac{R^{\lambda}_k(x^{\lambda}_j)\,\omega^{\lambda}_j}{\gamma^{\lambda}_k} .
$$

The paper states explicitly that convergence analysis of the collocation scheme "seems nontrivial and largely open", and that self-assessment is worth preserving.

The multidimensional case moves to the Fourier domain, where the problem is diagonal:

$$
(|\xi|^{\alpha}+\rho)\hat u(\xi)=\hat f(\xi),
\qquad
\hat u(\xi)=a(\xi)\hat f(\xi),\quad a(\xi)=\frac{1}{|\xi|^{\alpha}+\rho} .
$$

The algorithm interpolates $f$ on the tensor grid, transforms, multiplies by $a(\xi)$ and projects back. A self-dual structure helps: the inverse transform of $R^{\lambda}_n$ is computed by the **same** expansion formulas. There is an extra reason to choose a rational basis here — a function's Fourier transform typically decays more slowly than the function itself, so an algebraically decaying basis is exactly the right tool in the transformed variable.

## 92: uniformity rather than speed

Paper 92 closes a **uniformity** gap. The abstract states that existing methods do not work well as $\alpha$ approaches $0$ or $1$, and the structural reason is clear: the standard integral representations carry a prefactor $\sin(\alpha\pi)$ that vanishes at integer $\alpha$, and the strength of the integrand's endpoint behaviour is governed by $\alpha$. The Balakrishnan form and its negative-power counterpart are

$$
A^{\alpha}=\frac{\sin(\alpha\pi)}{\alpha\pi}\,A\int_0^{\infty}\bigl(t^{1/\alpha}I+A\bigr)^{-1}\mathrm dt,
\qquad
\mathcal L^{-\alpha}=\frac{2\sin(\alpha\pi)}{\pi}\int_0^{\infty}t^{2\alpha-1}
\bigl(\mathcal I+t^{2}\mathcal L\bigr)^{-1}\mathrm dt .
$$

In the best-uniform-rational-approximation family the error asymptotics carry the same factor: the best $(k,k)$ rational approximation to $t^{\beta-\alpha}$ on $[0,1]$ satisfies

$$
\lim_{k\to\infty}e^{2\pi\sqrt{(\beta-\alpha)k}}E_{\alpha}(k,k;\beta)
=4^{\,1+\beta-\alpha}\,\bigl|\sin\pi(\beta-\alpha)\bigr| ,
$$

root-exponential convergence in which $|\sin\pi(\beta-\alpha)|$ is exactly the degeneracy at $\alpha\to0,1$.

The paper targets cost and accuracy uniform in $\alpha\in(0,1)$ for two objects: $A^{-\alpha}$ with $A$ positive definite, and $(q\mathcal I+A^{\alpha})^{-1}$ with $q>0$. The second is the resolvent needed for an implicit time step of a space-fractional evolution equation, which explains the applications listed.

> [!warning] What could be verified
> The full text of this paper could not be reached through public channels and no preprint exists. What is confirmable: the journal keywords are trapezoidal rule, Laguerre-Gauss quadrature, fractional Laplacian, fractional Poisson equation and space-fractional Allen-Cahn equation, and the reference list includes the AAA rational approximation algorithm, double-exponential quadrature for fractional diffusion, the double-exponential formula for matrix fractional powers, and the exponentially convergent trapezoidal rule. That identifies the method as a **quadrature-based scheme** compared against double-exponential quadrature and rational approximation. This site does not report its theorems, convergence rates or the exact sense of "uniformly fast"; those require the published article.

## Side-by-side comparison

| No. | Basis or representation           | Far-field decay        | How the operator is handled                                    |
| --- | --------------------------------- | ---------------------- | -------------------------------------------------------------- |
| 27  | Hermite and over-scaled functions | Gaussian               | explicit physical-space differentiation matrices               |
| 42  | modified mapped Gegenbauer        | $\|x\|^{-(\lambda+1)}$ | decomposition into two elementary rational shapes              |
| 92  | no basis (a matrix function)      | not applicable         | integral representation plus quadrature (limited verification) |

## Coverage check

| Item                                             | Paper  | Status                                                         |
| ------------------------------------------------ | ------ | -------------------------------------------------------------- |
| Model problem and both operator definitions      | 27, 42 | hypersingular integral, constant, symbol form                  |
| Why all three standard routes fail               | 42     | truncation, transparent conditions, exponential bases          |
| Both Hermite bases                               | 27     | normalised and over-scaled forms with orthogonality            |
| Map, modified basis and orthogonality            | 42     | map identities, absorbed weight, recurrence, far-field limit   |
| Fractional Laplacian of the two shapes           | 42     | closed forms, the constant $A^{\gamma}_s$, Bessel proof route  |
| Divergence between integer and fractional order  | 42     | the three-case conclusion and its meaning                      |
| Recurrence evaluation and scaling parameter      | 42     | contiguous relation, roles of $\mu$ and $\lambda$              |
| Approximation space, nodes, both schemes         | 42     | equivalence, mapped Gauss data, matrices, open problem         |
| Multidimensional Fourier algorithm, self-duality | 42     | diagonalisation, algorithm steps, reason for the basis         |
| Uniformity gap and $\sin(\alpha\pi)$             | 92     | both representations, rational error asymptotics, both targets |

## Sources for this page

- T. Tang, H. Yuan, and T. Zhou, [_Hermite spectral collocation methods for fractional PDEs in unbounded domains_](https://doi.org/10.4208/cicp.2018.hh80.12), Commun. Comput. Phys. 24(4) (2018), pp. 1143-1168 (preprint [arXiv:1801.09073](https://arxiv.org/abs/1801.09073)).
- T. Tang, L.-L. Wang, H. Yuan, and T. Zhou, [_Rational spectral methods for PDEs involving fractional Laplacian in unbounded domains_](https://doi.org/10.1137/19M1244299), SIAM J. Sci. Comput. 42(2) (2020), pp. A585-A611 (preprint [arXiv:1905.02476](https://arxiv.org/abs/1905.02476)).
- Y. Duan, F. Zeng, H. Zhang, and T. Zhou, [_Fast computation of the fractional power of a matrix_](https://doi.org/10.1137/25M1757411), SIAM J. Sci. Comput. 48(1) (2026), pp. A309-A334.
