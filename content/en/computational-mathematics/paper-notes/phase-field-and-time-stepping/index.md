---
title: Phase-Field Models and Variable-Step Time Discretisation
description: Twelve papers that reduce energy arguments under adaptive steps to the positive definiteness of one real quadratic form
lang: en
translation: computational-mathematics/paper-notes/phase-field-and-time-stepping
tags:
  - computational-mathematics
  - paper-notes
  - phase-field
---

This topic holds 12 papers written between 2019 and 2026, mostly with Hong-lin Liao and Tao Tang. They address one concrete and difficult problem: **phase-field equations need adaptive time steps, and adaptive steps destroy the telescoping identity that energy arguments rely on.**

![Energy arguments for variable steps rest on a quadratic form](assets/diagrams/tao-zhou-papers/en/variable-step-energy.svg)

## Where the difficulty comes from

The Allen-Cahn equation

$$
\partial_{t}u=\varepsilon^{2}\Delta u-f(u),\qquad f(u)=u^{3}-u
$$

is the $L^2$ gradient flow of the Ginzburg-Landau free energy

$$
E[u](t)=\int_{\Omega}\Bigl(\tfrac12\varepsilon^{2}|\nabla u|^{2}+F[u]\Bigr)\mathrm d\mathbf x,
\qquad F[u]=\tfrac14(1-u^{2})^{2},
$$

so at the continuous level it satisfies both energy dissipation $\mathrm dE/\mathrm dt\le0$ and the maximum bound principle ($|u(\cdot,0)|\le1$ implies $|u(\cdot,t)|\le1$). Its solutions evolve on two widely separated time scales: fast initial dynamics followed by extremely slow coarsening. Variable time steps are therefore a practical necessity.

Written as a discrete convolution $D_{2}v^{n}=\sum_{k=1}^{n}b^{(n)}_{n-k}\nabla_{\tau}v^{k}$, the variable-step BDF2 formula

$$
D_{2}v^{n}=\frac{1+2r_{n}}{\tau_{n}(1+r_{n})}\nabla_{\tau}v^{n}
-\frac{r_{n}^{2}}{\tau_{n}(1+r_{n})}\nabla_{\tau}v^{n-1},
\qquad r_{k}:=\frac{\tau_{k}}{\tau_{k-1}}
$$

has kernels

$$
b^{(n)}_0=\frac{1+2r_n}{\tau_n(1+r_n)},
\qquad
b^{(n)}_1=-\frac{r_n^2}{\tau_n(1+r_n)},
\qquad
b^{(n)}_j=0\ (j\ge2).
$$

That $b^{(n)}_1<0$ is the root of the whole difficulty: the kernels are not non-negative, and with variable steps BDF2 is not self-adjoint and its kernels are not of convolution (Toeplitz) type, so neither the standard energy test-function trick nor the standard positivity argument applies.

## Three close readings

| Close reading                                                                                                                           | Papers                 | What they share                                    |
| --------------------------------------------------------------------------------------------------------------------------------------- | ---------------------- | -------------------------------------------------- |
| [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field\|Time-fractional phase field]]     | 40, 43, 57             | energy laws for fractional gradient flows          |
| [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf\|Variable-step BDF and convolution kernels]] | 48, 52, 58, 67, 69, 74 | kernel recombination and quadratic-form positivity |
| [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation\|IMEX and relaxation schemes]]             | 78, 91, 104            | linearisation and time adaptability                |

## Three step-ratio thresholds and where each belongs

Several specific constants recur in this literature, each attached to a **different** property. Mixing them up produces incorrect citations.

| Threshold                                 | Source and role                                                                                                                                        |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| $r_k<1+\sqrt2\approx2.4142$               | condition S0 of paper 48, for the **discrete maximum principle** and max-norm convergence; coincides with Grigorieff's (1983) zero-stability condition |
| $r_k<\tfrac{3+\sqrt{17}}{2}\approx3.5616$ | condition S1 of paper 48, for **energy stability**; paper 52 uses the same threshold for energy and $L^2$ stability                                    |
| $r_k<1.4877$                              | paper 67, for the discrete gradient structure and energy dissipation law of variable-step BDF3                                                         |

The origin of $3.5616$ is an exact positivity computation. After bounding $D_{2}u^{n}\nabla_\tau u^{n}$ from below, the first two terms telescope into a modified energy

$$
\widehat{E}[u^{k}]=E[u^{k}]+\frac{r_{k+1}\tau_{k}}{2(1+r_{k+1})}
\sum_{i}\bigl(\partial_{\tau}u_{i}^{k}\bigr)^{2},
$$

and the remaining term requires

$$
\frac{r_{k+1}}{1+r_{k+1}}<\frac{r_{s}}{1+r_{s}}=\frac{\sqrt{17}-1}{4}\approx0.78,
\qquad
r_{s}=\frac{3+\sqrt{17}}{2}\ \text{the positive root of}\ 2+3r-r^{2}=0 .
$$

Note that the modified energy uses $r_{k+1}$, the **next** step ratio. That look-ahead is exactly what makes the telescoping work.

> [!warning] About the constant $4.8645$
> Another frequently cited threshold in this literature is $r^\ast$, the positive root of $1+2r-r^{3/2}=0$, exactly $4.864536512317583$. On checking, it appears in **none** of the papers in this topic — a search for `4.86` across the ten available full texts returns no hit. It comes from work by Liao, Ji, Wang and Zhang on convex-splitting BDF2 for the Cahn-Hilliard equation (_J. Sci. Comput._ 92 (2022) 52, preprint [arXiv:2102.03731](https://arxiv.org/abs/2102.03731)), where it delivers **A-stability** of variable-step BDF2, a mesh-robust $L^2$ estimate and modified-energy dissipation, improving Liao and Zhang's $3.561$.
>
> Two details are where the confusion originates. First, **that paper prints $4.864$**; the more widely circulated form $4.8645$ comes from Zhang and Zhao (_J. of Math._ (PRC) 41(6) (2021) 471-488) for linear reaction-diffusion. Second, the technique is often credited to Lemma A.1 of Liao, Ji and Zhang (_IMA J. Numer. Anal._, phase field crystal), but that paper states only $3.561$ throughout. Papers 69 and 91 both cite Liao-Ji-Wang-Zhang, which is most likely how the constant entered this reading list. It should not be attributed to the papers here.
>
> Note also that paper 57 imposes **no step-ratio restriction at all** for its energy law, and only a step-size bound for the maximum bound principle.

## Core idea of each paper

### Time-fractional phase field

- **40** establishes an energy dissipation theory and numerical stability for time-fractional phase-field equations, proving $E[u](t)\le E[u](0)$ and the maximum bound principle. That energy statement is weaker than the classical dissipation law $\mathrm dE/\mathrm dt+\|\delta E/\delta u\|^2=0$ and is not a **differential** law, so it cannot be used the way the classical law is — which is exactly what paper 57 takes up.
- **43** gives a second-order nonuniform-step maximum-principle-preserving scheme for the time-fractional Allen-Cahn equation.
- **57** changes the question. Instead of asking whether the original energy dissipates, it defines a **variational energy** that does. The key move is to rewrite the equation as $\partial_t u=-{}^{R}\!\partial_t^{1-\alpha}(\delta E/\delta u)$, making the time derivative on the left local and moving the nonlocality onto the variational derivative, where a positivity property is available. The variational energy and its law are

  $$
  \mathcal{E}_{\alpha}[u]=E[u]+\tfrac12\,\mathcal{I}_t^{\alpha}\bigl\|\delta E/\delta u\bigr\|^{2},
  \qquad
  \frac{\mathrm d\mathcal{E}_{\alpha}}{\mathrm dt}
  +\tfrac12\,\omega_{\alpha}(t)\bigl\|\delta E/\delta u\bigr\|^{2}\le0 ,
  $$

  and as $\alpha\to1$ this recovers the classical dissipation law. The discrete version uses a variable-step L1$_R$ formula whose kernel positivity needs **no step-ratio restriction whatsoever**.

### Variable-step BDF and convolution kernels

- **48** is the paper that establishes both famous step-ratio constants, and establishes them for two different properties: S1 for energy stability and S0 for the discrete maximum principle and max-norm convergence. The scheme is fully implicit, with no stabilisation and no convex splitting.
- **52** carries the same analysis to the molecular beam epitaxy model without slope selection.
- **58** introduces a new discrete energy technique for multistep backward difference formulas.
- **67** pushes the discrete energy analysis to variable-step third-order BDF, giving a gradient structure and energy dissipation law under $r_k<1.4877$.
- **69** analyses the stability and convergence of the variable-step time-filtered backward Euler scheme.
- **74** studies the positive definiteness of the real quadratic forms produced by variable-step L1-type approximations of convolution operators, supplying the algebraic criteria that the earlier papers depend on.

### IMEX and relaxation schemes

- **78** treats phase-field models by a linear relaxation method with regularised energy reformulation, in the family that improves on scalar auxiliary variable and invariant energy quadratisation approaches.
- **91** gives a class of refined implicit-explicit Runge-Kutta methods with robust time adaptability and unconditional convergence for the Cahn-Hilliard model.
- **104** analyses the stability of implicit-explicit multistep methods for nonlinear parabolic equations by a semi-generating function approach.

> [!note] Coverage status
> Papers 40, 43, 48, 52, 57, 58, 67, 74, 91 and 104 have been checked equation by equation against the authors' own arXiv sources, and paper 78 against the full text of the published PDF, so for those eleven the equations, theorem statements and constants are transcribed rather than paraphrased. The one exception is paper 69: it has no preprint and neither the publisher nor any aggregator supplies the full text, so only what the complete abstract and reference list confirm is reported, leaving its explicit modified energy, kernel definitions, convergence order and numerical tests unverified.

## One transferable judgement

This group compresses "variable-step stability", which looks scheme-specific, into a single algebraic object: **the history terms collect into a real quadratic form whose positive definiteness depends only on the step-ratio sequence.** Once that positivity has algebraic criteria — which is what paper 74 supplies — the same toolkit transfers to third-order BDF, time-filtered Euler, fractional L1 rules and implicit-explicit Runge-Kutta methods. It is the same strategy as "write the sequential structure as one all-at-once operator and study its spectrum" in the [[en/computational-mathematics/paper-notes/parallel-in-time/index|parallel-in-time topic]]: first locate the one algebraic object that carries the whole difficulty.

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]; per-paper references appear at the end of each close-reading page.
