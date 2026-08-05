---
title: Variable-Step BDF and Convolution Kernels
description: Papers 48, 52, 58, 67, 69 and 74 - collecting the history terms into one object whose positivity can be decided
lang: en
translation: computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf
tags:
  - paper-notes
  - phase-field
  - discrete-energy
---

> [!note] Coverage of this page
> Papers **48** (_SIAM J. Numer. Anal._ 58(4), 2020), **52** (_Sci. China Math._ 64, 2021), **58** (_CSIAM Trans. Appl. Math._ 3, 2022), **67** (_J. Comput. Math._ 41, 2023), **69** (_BIT Numer. Math._ 63:39, 2023) and **74** (_Sci. China Math._ 67, 2024).
>
> All but paper 69 have been checked equation by equation against preprint full texts. Paper 69 has no preprint — confirmed by sweeping all 34 preprints under Liao's name — and neither the publisher nor any aggregator supplies the full text, so this page reports only what the complete abstract and reference list confirm: the step-ratio condition, the ingredients (a one-leg reformulation, a discrete gradient structure, and two new classes of DOC kernel) and the three results. Its explicit modified energy, the definitions of those kernels, the convergence order and the numerical tests are not verified here.

## 48: two step-ratio constants for two different properties

### The obstruction

The Allen-Cahn equation satisfies both energy dissipation and the maximum bound principle at the continuous level, and its solutions span two widely separated time scales, so variable steps are a practical necessity. At the time of writing, however, the numerical analysis of variable-step schemes was, in the authors' words, far from complete even for linear and semilinear parabolic equations. The specific gap: nobody had proved a **discrete maximum principle** for a second-order variable-step scheme for Allen-Cahn, and the sharpness of the step-ratio conditions in existing nonuniform BDF2 energy results was unclear.

The classical difficulty is that variable-step BDF2 is not self-adjoint and its kernels are not of convolution (Toeplitz) type, so neither the uniform-step energy test-function trick nor the standard positivity argument applies.

### The scheme and two named conditions

The variable-step BDF2 formula is

$$
D_{2}v^{n}=\frac{1+2r_{n}}{\tau_{n}(1+r_{n})}\nabla_{\tau}v^{n}
-\frac{r_{n}^{2}}{\tau_{n}(1+r_{n})}\nabla_{\tau}v^{n-1},
\qquad
r_{k}=\frac{\tau_{k}}{\tau_{k-1}},
$$

and as a discrete convolution $D_{2}v^{n}=\sum_{k=1}^{n}b^{(n)}_{n-k}\nabla_{\tau}v^{k}$ the kernels are

$$
b^{(n)}_0=\frac{1+2r_n}{\tau_n(1+r_n)},
\qquad
b^{(n)}_1=-\frac{r_n^2}{\tau_n(1+r_n)},
\qquad
b^{(n)}_j=0\ (j\ge2),
$$

and $b^{(n)}_1<0$ is exactly why the kernels are not non-negative and need recombining. The fully discrete scheme is **fully implicit** with no stabilisation and no convex splitting:

$$
D_{2}u^{n}=\varepsilon^{2}\Lambda_{h}u^{n}-f(u^{n}),
\qquad f(u^n)=(u^n)^{.3}-u^n .
$$

The paper uses two named conditions:

- **S1**: $0<r_{k}<\dfrac{3+\sqrt{17}}{2}\approx3.561$, for **energy stability**. The constant does not originate here: it comes from Lemma 2.1 of Liao and Zhang (_Math. Comp._ 90 (2021) 1207-1226) as the condition for positive definiteness of the BDF2 kernels, and both papers 48 and 52 import it as S1. Paper 52's own assessment of it is blunt — it calls it an artificial constant that is due to the condition S1.
- **S0**: $0<r_{k}<1+\sqrt{2}\approx2.414$, for the **discrete maximum principle** and max-norm convergence; the paper notes explicitly that this coincides with Grigorieff's (1983) zero-stability condition for ODE problems.

### The modified energy and the exact origin of 3.561

What decays is not the plain discrete energy but a modified one,

$$
\widehat{E}[u^{k}]=E[u^{k}]+\frac{r_{k+1}\tau_{k}}{2(1+r_{k+1})}
\sum_{i}\bigl(\partial_{\tau}u_{i}^{k}\bigr)^{2},
\qquad
\widehat{E}[u^{0}]=E[u^{0}],
$$

$$
E[u^{k}]=-\frac{\varepsilon^{2}}{2}(u^{k})^{T}\Lambda_{h}u^{k}
+\frac{1}{4}\sum_{i}\bigl(1-(u_{i}^{k})^{2}\bigr)^{2}.
$$

The correction is $O(\tau)$, so $\widehat E\to E$ as $\tau\to0$. Note it uses $r_{k+1}$, the **next** step ratio; that look-ahead is what makes the telescoping work.

Using $2a(a-b)=a^2-b^2+(a-b)^2$ and the kernel definitions,

$$
D_{2}u_{i}^{n}\,(\nabla_\tau u_{i}^{n})
\ \ge\
\frac{r_{n+1}\tau_{n}}{2(1+r_{n+1})}(\partial_{\tau}u_{i}^{n})^{2}
-\frac{r_{n}\tau_{n-1}}{2(1+r_{n})}(\partial_{\tau}u_{i}^{n-1})^{2}
+\Bigl(\frac{2+4r_{n}-r_{n}^{2}}{1+r_{n}}-\frac{r_{n+1}}{1+r_{n+1}}\Bigr)
\frac{\tau_{n}}{2}(\partial_{\tau}u_{i}^{n})^{2}.
$$

The first two terms telescope into $\widehat E$ and the last must be non-negative after absorbing the nonlinear contribution, which gives

$$
\frac{r_{k+1}}{1+r_{k+1}}<\frac{r_{s}}{1+r_{s}}=\frac{\sqrt{17}-1}{4}\approx0.78,
\qquad
r_{s}=\frac{3+\sqrt{17}}{2}\ \text{the positive root of}\ 2+3r-r^{2}=0 .
$$

The paper also derives the corresponding step-size bounds. With $h(x)=\dfrac{2+4x-x^{2}}{1+x}$ one has $h'(x)=\dfrac{x+1+\sqrt3}{(1+x)^{2}}(\sqrt3-1-x)$, so $h$ increases on $(0,\sqrt3-1]$ and decreases afterwards, with $h(0)=2$ and $h(\sqrt2+1)=1+\frac{\sqrt2}{2}$. Three regimes follow, and in the third ($\sqrt2+1<r_k<r_s$) one must additionally control the **next** step ratio.

### Where $1+\sqrt2$ comes from: kernel recombination

The maximum-principle threshold comes from a different technique, **kernel recombination and complementarity**. Introducing a parameter $\eta$ and

$$
\bar{v}^{0}=v^{0},\qquad \bar{v}^{k}=v^{k}-\eta\,v^{k-1}\ (k\ge1),
\qquad
v^{k}=\sum_{\ell=0}^{k}\eta^{k-\ell}\bar{v}^{\ell},
$$

substituting and exchanging the order of summation gives a recombined BDF2 formula whose kernels become non-negative for a suitable $\eta$. The zero-stability condition $r_k<1+\sqrt2$ is exactly the range in which that recombination is available.

**The two constants belong to two different properties, and citations must keep them apart.** Energy stability comes from positivity of the quadratic form (S1); the maximum principle comes from non-negativity after recombination (S0). They are not interchangeable.

## 52: the same analysis for the molecular beam epitaxy model

Paper 52 applies the variable-step BDF2 analysis to the molecular beam epitaxy model **without slope selection**, whose free energy is

$$
E_m[\phi]=\frac{\varepsilon}{2}\|\Delta\phi\|^{2}
-\frac{1}{2\varepsilon}\bigl\langle\ln|1+|\nabla\phi|^{2}|,1\bigr\rangle .
$$

The nonlinearity $-\nabla\phi/(1+|\nabla\phi|^{2})$ is bounded, so there is no maximum bound principle to exploit and the energy argument becomes the only tool. The paper establishes energy stability together with $L^2$ stability and convergence under the same threshold $r_k<(3+\sqrt{17})/2$.

## 58: bringing non-A-stable BDF-$k$ back to a textbook argument

### The problem

A-stable BDF1 and BDF2 admit a direct textbook discrete energy proof, testing with $u^n$ in $L^2$. For $3\le k\le5$ the BDF-$k$ formulas are **not** A-stable and that direct argument fails. Since Lubich, Mansour and Venkataraman the standard remedy has been the **Nevanlinna-Odeh multiplier technique**, which rests on Dahlquist's equivalence of A-stability and G-stability and tests with $u^n-\sum_i\eta_iu^{n-i}$ instead of $u^n$; another is Liu's telescope formulas. Both introduce artificial multipliers and, crucially, force a **stronger norm on the starting data** — $H^1$-type quantities appear in the multiplier-based estimates. The paper asks and answers whether a straightforward discrete energy analysis exists for BDF-$k$ with $3\le k\le5$.

### DOC kernels and orthogonality

Writing BDF-$k$ as a discrete convolution

$$
D_{k}v^n=\frac1{\tau}\sum_{j=1}^n b_{n-j}^{(k)}\,\nabla_\tau v^{j},
\qquad n\ge k,
$$

the kernels come from the generating function

$$
\sum_{\ell=1}^{k}\frac{1}{\ell}(1-\zeta)^{\ell-1}
=\sum_{\ell=0}^{k-1}b_{\ell}^{(k)}\zeta^{\ell},
$$

for example $(b_0,b_1,b_2)=(11/6,-7/6,1/3)$ for $k=3$ and $(137/60,-163/60,137/60,-21/20,1/5)$ for $k=5$.

The **discrete orthogonal convolution** kernels $\theta_j^{(k)}$ are defined recursively by

$$
\theta_{0}^{(k)}=\frac{1}{b_{0}^{(k)}},
\qquad
\theta_{n-j}^{(k)}=-\frac{1}{b_{0}^{(k)}}\sum_{\ell=j+1}^{n}\theta_{n-\ell}^{(k)}b_{\ell-j}^{(k)},
$$

and their key property is the orthogonality identity

$$
\sum_{\ell=j}^{n}\theta_{n-\ell}^{(k)}\,b_{\ell-j}^{(k)}\equiv\delta_{nj},
\qquad k\le j\le n,
$$

together with mutual orthogonality $\sum_{\ell=j}^{n}b_{n-\ell}^{(k)}\theta_{\ell-j}^{(k)}\equiv\delta_{nj}$. Consequently

$$
\sum_{j=k}^{n}\theta_{n-j}^{(k)}D_{k}u^{j}
=\frac1{\tau}u_{\mathrm I}^{(k,n)}+\partial_{\tau}u^{n},
\qquad
u_{\mathrm I}^{(k,n)}=\sum_{\ell=1}^{k-1}\nabla_\tau u^{\ell}
\sum_{j=k}^{n}\theta_{n-j}^{(k)}b_{j-\ell}^{(k)} .
$$

In other words, **acting with the DOC kernels inverts the BDF-$k$ convolution back to a single first difference**, up to a starting-value remainder. Because mutual orthogonality holds too, this is a **reversible** discrete transform with no loss of information. Taking $w=2\tau u^j$ in the transformed scheme and summing over $j$ then gives the classical energy inequality directly.

### Positivity decided by a generating function

The paper establishes that $\{b_j^{(k)}\}$ is positive (semi-)definite if and only if $\{\theta_j^{(k)}\}$ is, then decides positivity through the generating function of the symmetrised Toeplitz matrix,

$$
\mathrm g^{(k)}(\varphi)=2\sum_{j=0}^{k-1}b_j^{(k)}\cos(j\varphi),
$$

with Grenander-Szegő giving $\mathrm g^{(k)}_{\min}\le\lambda_{\min}(B_k)\le\lambda_{\max}(B_k)\le\mathrm g^{(k)}_{\max}$. For instance

$$
\mathrm g^{(3)}(\varphi)=\tfrac13\bigl(11-7\cos\varphi+2\cos2\varphi\bigr)
=\tfrac43\bigl(\cos\varphi-\tfrac78\bigr)^2+\tfrac{95}{48},
$$

so completing the square exhibits a positive lower bound. This yields, for $3\le k\le5$ and any real sequence, the quadratic-form bound

$$
2\sum_{m=k}^{n}w_m\sum_{j=k}^{m}b_{m-j}^{(k)}w_j
\ \ge\ \sigma_k\sum_{j=k}^{n}w_j^2
$$

with explicit constants $\sigma_k$.

> [!note] The mesh setting
> The time grid in this paper is **uniform**. Variable steps for $k\ge3$ are listed explicitly as an open problem in the concluding remarks — which is what paper 67 takes up.

## 67: variable-step BDF3

### Why a new analysis is needed

Variable-step BDF3 previously had essentially one classical result: Calvo and Grigorieff (2002) proved $L^2$ stability under the step-ratio condition $r_k<1.199$, with an estimate of the form

$$
\|u^n\|\le C\exp(C\Gamma_n)\Bigl(\|u_0\|+\sum_{j=1}^{n}\tau_j\|f^j\|\Bigr),
\qquad \Gamma_n=\sum_{k=2}^{n}|r_k-r_{k-1}| .
$$

**The prefactor $\exp(C\Gamma_n)$ is not mesh robust.** For the alternating sequence $\{\tau_1,\mu\tau_1,\tau_1,\mu\tau_1,\dots\}$ with $\mu\ne1$ and $T=\frac M2(1+\mu)\tau_1$ fixed, $\Gamma_M=(M-1)|\mu-\mu^{-1}|\to\infty$ as $\tau_1\to0$, so the bound degenerates precisely in the adaptive regime it is meant to cover. This paper replaces it with an analysis whose constants are independent of the ratios entirely, at the price of a ratio threshold.

### Variable-step BDF3 and its DOC kernels

The variable-step BDF3 formula is

$$
D_3v^n=d_0(r_n,r_{n-1})\partial_\tau v^{n}
+d_1(r_n,r_{n-1})\partial_\tau v^{n-1}
+d_2(r_n,r_{n-1})\partial_\tau v^{n-2},
$$

$$
d_0(x,y)=\frac{1+2x}{1+x}+\frac{xy}{1+y+xy},
\qquad
d_2(x,y)=\frac{xy^2}{1+y+xy}\cdot\frac{1+x}{1+y},
$$

with $d_1$ the combination that makes the formula consistent. The variable-step DOC kernels are

$$
\vartheta_0^{(n)}=\frac{1}{d_0^{(n)}},
\qquad
\vartheta_{n-j}^{(n)}=-\frac{1}{d_0^{(j)}}\sum_{i=j+1}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-j},
$$

satisfying the orthogonality identity $\sum_{i=j}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-j}\equiv\delta_{nj}$, in matrix form $\Theta_3D_3=I$, and since $D_3\Theta_3=I$ also mutual orthogonality. The paper establishes a discrete gradient structure, an energy dissipation law, and $L^2$ stability and convergence under the threshold $r_k<1.4877$.

### Where $1.4877$ comes from

The threshold is not a natural constant but the outcome of a **parameter trade-off**, which is worth spelling out. With $\gamma=7/10$, $R_e$ is the **unique positive root** of

$$
d_1(R_e,0)+\tfrac{7}{10}\sqrt{R_e}\,d_2(R_e,R_e)=0
\qquad\Longleftrightarrow\qquad
\frac{10}{7(R_e+1)}-\frac{R_e^2\sqrt{R_e}}{R_e^2+R_e+1}=0,
$$

numerically $R_e\approx1.4877$.

The value $\gamma=7/10$ arises as follows. The discrete gradient decomposition requires two conditions, $q_{n+1}\ge0$ and $p_{n+1}>0$, in the five variables $r_{n+1},r_n,r_{n-1},\gamma,R_e$, which cannot be solved exactly. The paper works on a constant-ratio grid instead: $q_{n+1}\ge0$ with $r_{n-1}=0$ and $r_{n+1}=r_n=r$ forces $\gamma\le-d_1(r,0)/(\sqrt r\,d_2(r,r))$, and imposing the second condition on a constant-ratio grid gives $\bar R_e\approx1.4965$ with $\bar\gamma\approx0.6924$. The authors then **fix $\gamma=7/10$**, close to that $0.6924$, for tractability, which yields $R_e\approx1.4877$. The justification is that $q_{n+1}\ge0$ is necessary and sharp whereas $p_{n+1}>0$ can be relaxed.

The sharpness of this threshold can therefore be **quantified**: numerically $R_e<1.69$ is necessary, while the theory delivers $R_e<1.4877$ as sufficient, so the gap is small. That contrasts instructively with paper 52's verdict on $3.561$ as an artificial constant — both are chosen for convenience, but here the distance to necessity has been measured.

## 69 and 74: making the tool itself the object of study

**Paper 69** analyses the stability and convergence of the variable-step time-filtered backward Euler scheme. Time filtering is a post-processing technique that raises the accuracy of a low-order scheme, and under variable steps its stability again reduces to this family of convolution-kernel analyses. Its energy stability and $L^2$ error estimate hold under the step-ratio condition

$$
\tfrac12\le r_k\le2,
$$

which is **not the same kind of object** as the other thresholds on this page, for two separate reasons.

First, it is **two-sided**. Every other result here restricts only the **upper** step ratio ($1+\sqrt2$, $3.561$, $1.4877$), whereas the ratio is bounded below by $1/2$ as well, so the steps may not be **cut** too abruptly either. This reflects the loss of A-stability for the variable-step filtered backward Euler scheme.

Second, $[1/2,2]$ is **not a sharp analytic threshold** but the standard heuristic safeguard used in adaptive codes; the paper itself calls it a practical constraint. It should therefore not be quoted alongside $1+\sqrt2$, $3.561$ or $1.4877$ as though they were quantities of the same type.

**Paper 74** studies **the tool itself**: the positive definiteness of the real quadratic forms produced by variable-step L1-type approximations of convolution operators. It supplies the algebraic criteria that the other papers invoke repeatedly — [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|paper 57]], for example, cites three inequalities on the kernel sequences from it when proving monotonicity of its DOC kernels.

**Singling out paper 74 makes the point of the whole series visible**: the technical core is not any one scheme but the algebraic question of when the real quadratic form collected from the history terms is positive definite. Once that question has independent criteria, the same argument transfers to third-order BDF, time-filtered Euler, fractional L1 rules and implicit-explicit Runge-Kutta methods.

## How the six relate

| No. | Object                             | Mesh     | Threshold                        | Core tool                                   |
| --- | ---------------------------------- | -------- | -------------------------------- | ------------------------------------------- |
| 48  | Allen-Cahn, BDF2                   | variable | S1 $3.561$; S0 $1+\sqrt2$        | modified energy; kernel recombination       |
| 52  | molecular beam epitaxy, BDF2       | variable | $3.561$                          | modified energy                             |
| 58  | linear reaction-diffusion, BDF-$k$ | uniform  | none for $3\le k\le5$            | DOC kernels; generating-function positivity |
| 67  | diffusion, BDF3                    | variable | $1.4877$                         | variable-step DOC kernels                   |
| 69  | parabolic, filtered Euler          | variable | $[1/2,2]$ (two-sided, heuristic) | same kernel analysis                        |
| 74  | the quadratic form itself          | variable | supplies the criteria            | algebraic criteria                          |

## Coverage check

| Item                                             | Paper | Status                                                      |
| ------------------------------------------------ | ----- | ----------------------------------------------------------- |
| Variable-step BDF2 formula and kernel signs      | 48    | formula, convolution kernels, meaning of $b^{(n)}_1<0$      |
| Conditions S0 and S1 with their properties       | 48    | both constants, both properties, agreement with Grigorieff  |
| Modified energy and the look-ahead term          | 48    | form, $O(\tau)$, why $r_{k+1}$ appears                      |
| Exact derivation of 3.561 and step bounds        | 48    | lower-bound inequality, positive root, monotonicity of $h$  |
| Kernel recombination and $1+\sqrt2$              | 48    | the $\eta$ transform, substitution, origin of the threshold |
| MBE free energy and the absent maximum principle | 52    | energy form and the resulting limitation                    |
| DOC kernels for BDF-$k$ and reversibility        | 58    | generating function, kernel table, orthogonality, inversion |
| Generating function and Grenander-Szegő          | 58    | $\mathrm g^{(k)}$, completed square, quadratic-form bound   |
| Non-robustness of the Calvo-Grigorieff prefactor | 67    | alternating-mesh counterexample, divergence of $\Gamma_M$   |
| Variable-step BDF3 and its DOC kernels           | 67    | $d_0,d_2$, DOC recursion, both orthogonalities, threshold   |
| Quadratic-form positivity as its own object      | 74    | its role across the whole thread                            |

## Sources for this page

- H.-l. Liao, T. Tang, and T. Zhou, [_On energy stable, maximum-principle preserving, second-order BDF scheme with variable steps for the Allen-Cahn equation_](https://doi.org/10.1137/19M1289157), SIAM J. Numer. Anal. 58(4) (2020), pp. 2294-2314 (preprint [arXiv:2003.00421](https://arxiv.org/abs/2003.00421)).
- H.-l. Liao, X. Song, T. Tang, and T. Zhou, [_Analysis of the second-order BDF scheme with variable steps for the molecular beam epitaxial model without slope selection_](https://doi.org/10.1007/s11425-020-1817-4), Sci. China Math. 64 (2021), pp. 887-902.
- H.-l. Liao, T. Tang, and T. Zhou, [_A new discrete energy technique for multi-step backward difference formulas_](https://doi.org/10.4208/csiam-am.SO-2021-0032), CSIAM Trans. Appl. Math. 3 (2022), pp. 318-334 (preprint [arXiv:2102.04644](https://arxiv.org/abs/2102.04644)).
- H.-l. Liao, T. Tang, and T. Zhou, [_Discrete energy analysis of the third-order variable-step BDF time-stepping for diffusion equations_](https://doi.org/10.4208/jcm.2207-m2022-0020), J. Comput. Math. 41 (2023), pp. 325-344 (preprint [arXiv:2204.12742](https://arxiv.org/abs/2204.12742)).
- H.-l. Liao, T. Tang, and T. Zhou, [_Stability and convergence of the variable-step time filtered backward Euler scheme for parabolic equations_](https://doi.org/10.1007/s10543-023-00982-y), BIT Numer. Math. 63 (2023), 39.
- H.-l. Liao, T. Tang, and T. Zhou, [_Positive definiteness of real quadratic forms resulting from the variable-step L1-type approximations of convolution operators_](https://doi.org/10.1007/s11425-022-2229-5), Sci. China Math. 67 (2024), pp. 237-252.
