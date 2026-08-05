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

All three pages treat each paper in the same order: the intuition before any formalism, then the setting, the full derivation, the theorems with their hypotheses, and finally the concrete configuration and measured numbers of the experiments. Eleven of the twelve papers can be covered at that depth; the exception is paper 69, for the reason given in the coverage note below.

## Three step-ratio thresholds and where each belongs

Several specific constants recur in this literature, each attached to a **different** property. Mixing them up produces incorrect citations.

| Threshold                                 | Source and role                                                                                                                                        |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| $r_k<1+\sqrt2\approx2.4142$               | condition S0 of paper 48, for the **discrete maximum principle** and max-norm convergence; coincides with Grigorieff's (1983) zero-stability condition |
| $r_k<\tfrac{3+\sqrt{17}}{2}\approx3.5616$ | condition S1 of paper 48, for **energy stability**; paper 52 uses the same threshold for energy and $L^2$ stability                                    |
| $r_k<R_e\approx1.4877$                    | paper 67, for the discrete gradient structure and energy dissipation law of variable-step BDF3                                                         |

The three govern different properties, and this is worth repeating: $1+\sqrt2$ governs the maximum principle and the maximum norm, $3.5616$ governs the energy, and although both come from the same paper (48) they are not interchangeable; $1.4877$ belongs to a formula of a different order altogether. Paper 52 calls $3.5616$ outright "an artificial constant that is due to the condition S1" — that is, a product of the analysis rather than an intrinsic limit of the scheme.

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

Paper 67's $R_e$ has an exact definition: it is the unique positive root of $d_1(R_e,0)+\tfrac{7}{10}\sqrt{R_e}\,d_2(R_e,R_e)=0$, equivalently of

$$
\frac{10}{7(R_e+1)}-\frac{R_e^{2}\sqrt{R_e}}{R_e^{2}+R_e+1}=0 .
$$

The $\gamma=7/10$ here is a tractable stand-in for the $0.6924$ that the constant-ratio analysis produces (which would give $R_e\approx1.4965$). The substitution is legitimate because the condition $q_{n+1}\ge0$ is necessary and sharp whereas $p_{n+1}>0$ can be relaxed. Numerically $R_e<1.69$ is **necessary** while the theory delivers $1.4877$ as **sufficient**, so the gap is small and quantified. For bearings, $1.4877$ sits below the BDF3 threshold $R_3\approx2.553$ that Li and Liao obtain for nonlinear ODEs, and above Calvo and Grigorieff's (2002) $1.199$ for BDF3 in $L^2$ — whose error prefactor $\exp(C\Gamma_n)$ is not mesh-robust, and removing it is exactly what paper 67 achieves.

> [!warning] Paper 69's $[1/2,2]$ is not the same kind of object
> The step-ratio condition in the variable-step time-filtered backward Euler paper is a **two-sided** interval $r_k\in[1/2,2]$, and it is the standard practical safeguard of adaptive codes rather than a sharp algebraic threshold determined by the root of some polynomial. It reflects the loss of A-stability for variable-step FiBE. Listing it alongside $1+\sqrt2$, $3.5616$ and $1.4877$ in one "table of thresholds" is misleading: each of those three is the exact boundary of a positivity or zero-stability condition, and this one is not.

## Core idea of each paper

### Time-fractional phase field

- **40** establishes an energy dissipation theory and numerical stability for time-fractional phase-field equations, proving $E[u](t)\le E[u](0)$ and the maximum bound principle. That energy statement is weaker than the classical dissipation law $\mathrm dE/\mathrm dt+\|\delta E/\delta u\|^2=0$ and is not a **differential** law, so it cannot be used the way the classical law is — which is exactly what paper 57 takes up.
- **43** gives a second-order nonuniform-step maximum-principle-preserving scheme for the time-fractional Allen-Cahn equation, answering head-on the "higher-order schemes are needed" problem that paper 40's concluding remarks pose. It brings the discrete complementary convolution (DCC) kernels and the fractional Grönwall machinery into this line of work, and proves the sharp convergence rate $\min\{\gamma\sigma,2\}$, in which the grading parameter $\gamma$ and the solution's regularity exponent $\sigma$ jointly set the order, so the initial layer is no longer an obstruction. Note that its step-ratio convention is $\rho_k=1/r_{k+1}$, the reciprocal of the one used in the other papers; quoting its numbers directly will produce errors.
- **57** changes the question. Instead of asking whether the original energy dissipates, it defines a **variational energy** that does. The key move is to rewrite the equation as $\partial_t u=-{}^{R}\!\partial_t^{1-\alpha}(\delta E/\delta u)$, making the time derivative on the left local and moving the nonlocality onto the variational derivative, where a positivity property is available. The variational energy and its law are

  $$
  \mathcal{E}_{\alpha}[u]=E[u]+\tfrac12\,\mathcal{I}_t^{\alpha}\bigl\|\delta E/\delta u\bigr\|^{2},
  \qquad
  \frac{\mathrm d\mathcal{E}_{\alpha}}{\mathrm dt}
  +\tfrac12\,\omega_{\alpha}(t)\bigl\|\delta E/\delta u\bigr\|^{2}\le0 ,
  $$

  and as $\alpha\to1$ this recovers the classical dissipation law. The discrete version uses a variable-step L1$_R$ formula whose kernel positivity needs **no step-ratio restriction whatsoever**.

### Variable-step BDF and convolution kernels

- **48** is the paper that establishes both famous step-ratio constants, and establishes them for two different properties: S1 for energy stability and S0 for the discrete maximum principle and max-norm convergence. The scheme is fully implicit, with no stabilisation and no convex splitting. Its technical core is **kernel recombination and complementarity** (KRC): the sign-indefinite BDF2 kernels $\{b^{(n)}_0>0,\,b^{(n)}_1<0\}$ are traded for a family of non-negative decreasing kernels $\{d^{(n)}_j\}$, after which the classical energy argument becomes available again.
- **52** carries the same analysis to the molecular beam epitaxy model **without slope selection**, obtaining energy stability together with $L^2$ stability and convergence at the same $3.5616$ threshold, and upgrading the tool from KRC to **discrete orthogonal convolution** (DOC) kernels. The difficulty here is that the nonlinearity of the logarithmic energy $-\tfrac12\ln(1+|\nabla\phi|^2)$ is not polynomially bounded; the conclusions state that the technique does not apply to the model with slope selection, which is precisely what paper 78 treats.
- **58** introduces a new discrete energy technique for BDF-$k$, $3\le k\le5$: DOC kernels plus the Grenander-Szegő theorem settle the positivity of a symmetrised Toeplitz form directly, replacing the traditional family-by-family search for Nevanlinna-Odeh multipliers. This is the uniform-grid branch, so no step ratio enters.
- **67** pushes the discrete energy analysis to variable-step third-order BDF, giving a discrete gradient structure and energy dissipation law under $r_k<R_e\approx1.4877$, and removing the non-mesh-robust prefactor $\exp(C\Gamma_n)$ from the Calvo-Grigorieff error estimate.
- **69** analyses the stability and convergence of the variable-step time-filtered backward Euler scheme. It belongs to the same programme as the papers above — a one-leg reformulation makes the left-hand side variable-step BDF2, so the gradient-structure machinery transfers — but it is the one paper in this topic whose full text could not be obtained; see the coverage note below.
- **74** studies the positive definiteness of the real quadratic forms produced by variable-step L1-type approximations of convolution operators, setting out in one place the definitions of the DOC and DCC kernels, their two identities ($\equiv\delta_{nk}$ for DOC, $\equiv1$ for DCC) and the algebraic criteria C1-C4 for positivity, which is what the earlier papers depend on. It is purely analytical, with no numerical experiments.

### IMEX and relaxation schemes

- **78** treats phase-field models by a linear relaxation method with regularised energy reformulation (RRER), in the family that improves on scalar auxiliary variable (SAV) and invariant energy quadratisation (IEQ) approaches: each step solves only a constant-coefficient linear problem, and introducing the auxiliary variable requires neither shifting the energy nor taking a square root of it. The models it treats are molecular beam epitaxy **with slope selection** and phase field crystal, which are not the models of paper 52.
- **91** gives a class of refined implicit-explicit Runge-Kutta methods with robust time adaptability and unconditional convergence for the Cahn-Hilliard model. Its central quantity is an average dissipation rate $\mathcal R$: for ordinary IERK methods $\mathcal R$ depends on $\tau$ and on the stabilisation parameter, so energy curves drift under adaptivity; the "refinement" condition makes $\mathcal R$ a constant independent of the discretisation parameters, and the energy curves become insensitive to $\tau_{\max}$.
- **104** analyses the stability of implicit-explicit multistep methods for nonlinear parabolic equations by a semi-generating function approach. It carries paper 58's DOC-plus-Grenander-Szegő route into the implicit-explicit setting: the one-sided ("semi") generating function lets the ratios $b/a$ and $c/a$ be handled together, and stability ends up controlled by three computable quantities, $\lambda_{\mathrm I}$, $\sigma_{\mathrm E}$ and $\sigma_{\mathrm F}$.

> [!warning] Paper 78 is the exception in this topic
> It has no Liao co-authorship, uses no DOC/DCC kernels, involves no variable steps and contains no step-ratio analysis; the connection to the rest is Tao Zhou. More importantly, its energy statement is of **modified-energy** type — an $\widehat E$ built from $(\phi,q)$ — which is exactly what papers 40, 48, 52, 57 and 91 work hard to **avoid**, since those prove decay of the original or variational energy. Paper 91's introduction criticises SAV-based high-order schemes explicitly for establishing stability "with respect to a modified energy involving the auxiliary variable."
>
> What saves paper 78 from that criticism is that its modified energy coincides with the original **exactly at the continuous level** ($\widehat E(\phi,q)\equiv E(\phi)$), so it is not an arbitrary substitute. The tension is worth keeping: it is the dividing line that makes this topic cohere — whether to run the stability argument on the original energy, or to swap in another object and then argue about its relation to the original.

## What the experiments keep showing

Eight of the twelve papers have experiments that can be reported item by item. Papers 58, 74 and 104 are purely analytical, their sources containing only plots of kernels or generating functions; paper 69's abstract says numerical tests are presented, but the specific examples cannot be verified. Taken together, the eight that can be reported deliver one and the same judgement: **the theoretical thresholds are sufficient conditions, and fairly conservative ones.**

| Paper | What the experiment tests directly | Outcome                                                                                            |
| ----- | ---------------------------------- | -------------------------------------------------------------------------------------------------- |
| 48    | random time meshes                 | second order in time observed; random steps stress step-ratio robustness head-on                   |
| 52    | random time meshes                 | second order survives $\max r_k=850.80$ and $N_1=49$ levels violating $3.5616$                     |
| 67    | periodic and random meshes         | third order ($2.98$–$3.00$) with about half the levels at $r_k\ge R_e$                             |
| 91    | adaptive stepping                  | at $T=450$, fixed $\tau=10^{-4}$ costs $6724.63$ s against $2.62$ s adaptively ($\tau_{\max}=0.5$) |

Paper 52's table is the most explicit of these: it deliberately reports $\max r_k$ and $N_1$, the number of levels violating the threshold, which amounts to the authors themselves showing that $3.5616$ is far from necessary. Paper 67 is the finer case — numerically $R_e<1.69$ is necessary, so there the gap is small and known. That contrast is itself one of the topic's conclusions: the BDF2 threshold barely constrains practice, while the BDF3 one is close to the truth.

> [!note] Coverage status
> Papers 40, 43, 48, 52, 57, 58, 67, 74, 91 and 104 have been checked equation by equation against the authors' own arXiv sources, and paper 78 against the full text of the published PDF, so for those eleven the equations, theorem statements and constants are transcribed rather than paraphrased. The three close readings now give all eleven an intuition, a full derivation and theorems with their hypotheses, and for the eight that ran experiments, the concrete configuration and the measured numbers. Papers 58, 74 and 104 ran none, and the pages say so rather than implying otherwise.
>
> The one exception is paper 69: it has no preprint (confirmed by sweeping all 34 preprints under Liao's name), and neither Springer, Semantic Scholar nor zbMATH supplies the full text, so only what the complete abstract and the 22-item reference list confirm is reported. Its step-ratio condition, its three ingredients (a one-leg reformulation, a discrete gradient structure, two new classes of DOC kernel) and its three results are verifiable; its explicit modified energy, its DOC kernel definitions, its convergence order and its numerical tests are not, and the page does not fill those in. The variable-step FiBE scheme shown there is reconstructed from the method source (DeCaria-Guzel-Layton-Li, [arXiv:1810.06670](https://arxiv.org/abs/1810.06670)) and is flagged as such on the page.

## One transferable judgement

This group compresses "variable-step stability", which looks scheme-specific, into a single algebraic object: **the history terms collect into a real quadratic form whose positive definiteness depends only on the step-ratio sequence.** Once that positivity has algebraic criteria — which is what paper 74 supplies — the same toolkit transfers to third-order BDF, time-filtered Euler, fractional L1 rules and implicit-explicit Runge-Kutta methods. It is the same strategy as "write the sequential structure as one all-at-once operator and study its spectrum" in the [[en/computational-mathematics/paper-notes/parallel-in-time/index|parallel-in-time topic]]: first locate the one algebraic object that carries the whole difficulty.

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]; per-paper references appear at the end of each close-reading page.
