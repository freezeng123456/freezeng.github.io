---
title: Densities as Invertible Maps - From Random Fields to Fractional Equations
description: Papers 62, 64, 72 and 87 - writing an unknown density as the pushforward of a reference distribution and training it on the equation residual
lang: en
translation: computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities
tags:
  - paper-notes
  - scientific-machine-learning
  - normalizing-flows
---

> [!note] Coverage of this page
> Papers **62** (_J. Comput. Phys._ 461, 2022), **64** (_Commun. Comput. Phys._ 32(2), 2022), **72** (_J. Sci. Comput._ 97:68, 2023) and **87** (_SIAM J. Sci. Comput._ 47(6), 2025).

![Turning a density equation into a parameterised problem with an invertible map](assets/diagrams/tao-zhou-papers/en/density-flow-solvers.svg)

The shared judgement of these four papers is that when the solution of an equation **is itself a probability density**, approximating it by values on a grid throws away structure. A density must be non-negative and integrate to one; in grid methods those are constraints, in an invertible map they are identities.

The pushforward formula

$$
p_{X}(x)=p_{Z}\bigl(f(x)\bigr)\,\bigl|\det\nabla_x f(x)\bigr|
$$

delivers non-negativity for free (the right-hand side is a density times a positive determinant) and normalisation for free (change of variables). Only two questions remain: how to parameterise $f$ so that it is expressive while $\det\nabla_x f$ stays computable, and what loss to train it on when no labelled data exists. The four papers each push on a different **target object**: random fields, time-dependent densities, fractional equations, bounded support.

## 62: a random field as a flow with Karhunen-Loève structure

### The idea

A random field is a family of correlated random variables indexed by spatial position. Modelling one means controlling two things at once: the **spatial correlation structure** and the **non-Gaussianity of the pointwise marginal law**. Each has a mature tool, and neither tool suffices alone.

The Karhunen-Loève expansion handles correlation: it writes the field as orthogonal modes times uncorrelated random coefficients. That is the standard representation for a Gaussian field, but it requires knowing the covariance kernel, and in a data-driven setting the covariance is precisely what is unknown; moreover, once the truncation order rises, building a surrogate in the parametric space runs into the curse of dimensionality — the arbitrary-polynomial-chaos deep-network methods the paper cites have term counts that grow exponentially with the effective dimension.

A normalising flow handles non-Gaussianity: an invertible network pushes a simple distribution into a complicated one. But a flow carries no spatial structure of its own; it acts on a vector of fixed dimension and knows nothing about the notion of a field.

Paper 62 chains the two together: **first build a Gaussian reference field carrying a Karhunen-Loève expansion structure, but let networks learn the expansion coefficients; then push that Gaussian field through a flow to obtain the target field.** The first half supplies spatial correlation without needing a covariance prior; the second half supplies non-Gaussianity.

Chaining them buys one more thing, and this is the substantive advantage over physics-informed GANs: **the flow is a bijection, so the likelihood can be written down explicitly.** The training criterion can be plain maximum likelihood rather than adversarial. A second benefit is easy to overlook: the model gives a distribution over field values at arbitrary locations rather than over a fixed-length vector, so **sensor positions need not be held fixed across observations** — which in real measurement campaigns is the norm, not the exception.

### Setting

The problem class is

$$
\mathcal N_x\bigl[u(x;\omega);k(x;\omega)\bigr]=f(x;\omega),
\quad x\in\mathcal D,\ \omega\in\Omega,
\qquad
\mathcal B_x\bigl[u(x;\omega)\bigr]=0,\quad x\in\Gamma,
$$

where $k$ is a random input field (a diffusivity, say), $f$ a random forcing, and $u$ the random solution field. Only scattered sensor readings are observed, and the covariance structure of the randomness is unknown.

### Derivation

The paper states the construction as three explicit steps:

1. Construct a reference Gaussian random field $z(x,\omega)$ with a **truncated Karhunen-Loève expansion structure**, whose expansion coefficients are parameterised by deep neural networks;
2. Construct a bijective transformation between the reference field and the target stochastic field — a normalising flow;
3. Train all parameters by maximising the sum of the log-likelihood over the scattered measurements.

For SDEs and stochastic PDEs the known physics enters as a residual loss added to the total, producing the **physics-informed** version. In other words the same representation works in pure-data mode (step 3 alone) and in physics-constrained mode (likelihood plus residual), while forward, inverse and mixed problems differ only in which quantities are observed and which are constrained by the residual. That is the concrete content of the paper's claim to a "unified framework".

The paper lists three advantages of its own: one framework for forward, inverse and mixed problems; no requirement that sensor positions stay fixed across observations; and relief from the curse of dimensionality that afflicts polynomial chaos.

> [!warning] How far this could be verified
> The algebraic form of the Karhunen-Loève-structured reference field, the form of the flow's coupling layers, and the combined likelihood-plus-physics loss could not be confirmed equation by equation from the verification material behind this page. The three steps above are the paper's own statement and can be quoted safely; **the specific formulas are not given on this page** and should be read directly from Section 3 of the original.

### Theorems

The abstract claims no theorem, and no convergence or approximation result was found in the body. This is worth saying plainly: the contribution is a representation and a training framework, not an analysis.

### Numerical experiments

Section 4 falls into three groups:

| Group | Content                                                          |
| ----- | ---------------------------------------------------------------- |
| One   | learning stochastic processes, including non-Gaussian and **mixed** non-Gaussian fields |
| Two   | forward stochastic elliptic equations                            |
| Three | inverse stochastic elliptic equations                            |

The mixed non-Gaussian group is the most informative of the three: it tests exactly whether the division of labour holds — Karhunen-Loève structure for correlation, flow for non-Gaussianity — because an insufficiently expressive flow would smear a mixture into a single mode.

The abstract's qualitative conclusion is that the model can learn non-Gaussian processes and solve several types of stochastic PDE.

> [!warning] How far the numerical results were verified
> The composition of the three groups is verified; **the error magnitudes are not**, so this page gives no numbers.

### Relation to the others

This is the earliest normalising-flow paper in the topic and it sets the template for papers 64, 72, 87 and for [[en/computational-mathematics/paper-notes/scientific-machine-learning/spectral-bias-and-generative-solvers|105]] and [[en/computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning|107]]: **represent an unknown distribution or field by an invertible network and train it against a physics residual rather than against labelled data.** It is the direct ancestor of paper 105 (FLUID), which replaces the Karhunen-Loève-structured reference field with a recurrent summary statistic and conditional flows. Paper 64 cites it in the text as prior work on flow-based surrogates for uncertainty quantification.

## 64: time is not an extra dimension

### The idea

The solution of a time-dependent Fokker-Planck equation is not a density but **a family of densities indexed by time**. There are three apparently reasonable ways to make a flow represent such a family.

The first is to train one flow per time slice. Computationally that is a disaster, and it guarantees nothing about continuity in time — two adjacent instants could be learned as entirely unrelated objects.

The second is to treat time as a $(d+1)$-st dimension and let the flow learn a joint space-time density. This route is **mathematically wrong**, for a blunt reason: the density integrates to one over space at each fixed time, but integrating over space and time together gives the length of the time interval, $\int p(x,t)\,\mathrm dx\,\mathrm dt\neq1$. The pushforward formula guarantees *joint* normalisation while the equation demands *per-instant* normalisation, and the two are incompatible.

The third is the paper's choice: **condition, do not augment.** One flow, but its coupling layers take $t$ as an extra input; the latent time coordinate is pinned to real time, $t^\ast=t$. The Jacobian then degenerates to $1$ in the time direction, only the spatial block remains, per-instant normalisation holds automatically, and one set of weights covers the whole time interval.

A further benefit has nothing to do with the flow itself: the training signal is the equation residual, so the method needs **neither labelled data nor sample paths of the stochastic differential equation**. Path-based methods rely on many Monte Carlo trajectories to estimate the density, with the attendant cost and noise; the residual route removes that entire block.

### Setting

The underlying SDE and the time-dependent Fokker-Planck equation are

$$
\mathrm d\bm X_t=\bm\mu(\bm X_t,t)\,\mathrm dt+\bm\sigma(\bm X_t,t)\,\mathrm d\bm W_t,
$$

$$
\frac{\partial p(\bm x,t)}{\partial t}
=-\sum_{i=1}^d\frac{\partial}{\partial x_i}\bigl[\mu_i(\bm x,t)p(\bm x,t)\bigr]
+\sum_{i=1}^d\sum_{j=1}^d\frac{\partial^2}{\partial x_i\partial x_j}
\bigl[D_{ij}(\bm x,t)p(\bm x,t)\bigr],
\qquad D=\tfrac12\sigma\sigma^{\mathsf T}.
$$

The stationary version drops the time derivative:

$$
\sum_i\frac{\partial}{\partial x_i}\bigl[\mu_i p\bigr]
+\sum_{i,j}\frac{\partial^2}{\partial x_i\partial x_j}\bigl[D_{ij}p\bigr]=0 .
$$

### Derivation

**From augmentation to conditioning.** Write the naive version first, with $\widehat{\bm x}=(\bm x,t)$ and $\widehat{\bm z}=(\bm z,t^\ast)$:

$$
p_{\widehat{\bm X}}(\widehat{\bm x})=p_{\widehat{\bm Z}}(\widehat{\bm z})\,|\det J|,
\qquad \widehat{\bm z}=f(\widehat{\bm x}).
$$

Now pin the latent time, $t^\ast=t$. That pinning makes the last row of the Jacobian $(0,\dots,0,1)$, and the determinant expands blockwise:

$$
\det J=\begin{vmatrix}\partial\bm z/\partial\bm x & \partial\bm z/\partial t\\ 0 & 1\end{vmatrix}
=\Bigl|\frac{\partial \bm z(\bm x,t)}{\partial\bm x}\Bigr| .
$$

Note that the upper-right block $\partial\bm z/\partial t$ is generally non-zero — the flow really does depend on time — but it falls outside the expansion and never enters the result. Hence

$$
p_{\widehat{\bm X}}(\bm x,t)=p_{\widehat{\bm Z}}(\bm z,t)\,
\Bigl|\frac{\partial\bm z}{\partial\bm x}\Bigr|,
\qquad \bm z=f(\bm x,t).
$$

This deserves its own line: **the flow is time-conditioned, not time-augmented.** One word decides whether normalisation holds.

**Composition and Jacobian.** The flow composes layerwise and the Jacobian multiplies layerwise:

$$
\bm z=f(\bm x,t)=f_{[L]}\circ\cdots\circ f_{[1]}(\bm x,t),
\qquad
\bm x=f^{-1}(\bm z,t),
$$

$$
\bigl|\det\nabla_{\bm x}f(\cdot,t)\bigr|
=\prod_{i=1}^L\bigl|\det\nabla_{\bm x_{[i-1]}}f_{[i]}(\cdot,t)\bigr| .
$$

**Architecture.** The paper describes its model as "a simplified extension of KRnet from the spatial domain to the temporal-spatial domain". Each $f_{[i]}$ is an Actnorm layer followed by a modified **time-dependent** affine coupling layer; the last layer is a polynomial spline transformation that increases modelling power.

The Actnorm (scale-and-bias) layer and its inverse are

$$
\bm y_{[i]}=\bm a_i\odot \bm x_{[i]}+\bm b_i,
\qquad
\bm x_{[i]}=\frac{\bm y_{[i]}-\bm b_i}{\bm a_i},
$$

with $\bm a_i,\bm b_i$ initialised data-dependently from mini-batch statistics.

The **time-dependent affine coupling layer** is the object to quote when contrasting with real NVP and with the plain KRnet:

$$
\bm x_{[i],1}=\bm x_{[i-1],1},
$$

$$
\bm x_{[i],2}=\bm x_{[i-1],2}\odot\Bigl(\bm 1_{d-m}
+\beta\tanh\bigl(\bm s_i(\bm x_{[i-1],1},t)\bigr)\Bigr)
+e^{\bm\zeta_i}\odot\tanh\bigl(\bm q_i(\bm x_{[i-1],1},t)\bigr),
$$

where $|\beta|<1$ is user-specified (for instance $0.6$), $\bm s_i,\bm q_i:\mathbb R^{m+1}\to\mathbb R^{d-m}$, and $\bm\zeta_i\in\mathbb R^{d-m}$ is trainable.

There are two differences from standard real NVP, each with a clear purpose. First, $\tanh$ together with $\beta$ clamps the scale factor inside $(1-\beta,1+\beta)$, whereas real NVP's scale is $\exp(\cdot)$ and can be arbitrarily large or arbitrarily close to zero; clamping keeps the log-determinant of the Jacobian from blowing up and stabilises training. Second, the input dimension of $\bm s_i,\bm q_i$ is $m+1$ rather than $m$, and the extra coordinate is $t$ — that is how one set of weights covers the whole time interval.

**Loss.** The training signal is the time-dependent Fokker-Planck residual, with no labelled data.

> [!note] The discrete form of the loss
> The explicit discrete residual expression was **reconstructed from prose** in the verification material behind this page rather than transcribed equation by equation: it is the mean-square Fokker-Planck residual evaluated at space-time collocation points $(\bm x^{(i)},t^{(i)})$, matching the practice of the ADDA line. Check the original before quoting the precise form.

### Theorems

The paper claims no theorem. Its abstract argues at the level of methodology: the scheme is mesh-free, needs no labelled data, and extends readily to high dimensions.

### Numerical experiments

The paper's own feature list amounts to four steps: model $p(\bm x,t)$ by a temporal normalising flow with latent time pinned to real time; build a physics-informed residual from the equation, using neither sample paths nor labels; minimise the mean-square residual over space-time collocation points; then adaptively refresh the collocation set from the current flow and repeat.

The experiments cover three categories:

| Category         | Content                                                       |
| ---------------- | ------------------------------------------------------------- |
| Linear drift     | time-dependent Fokker-Planck equations with a linear drift    |
| Nonlinear drift  | the same, with a nonlinear drift term                          |
| High dimension   | high-dimensional time-dependent problems                       |

The linear-drift group serves as calibration: such problems often have analytic or high-accuracy reference solutions, so one can tell whether the error comes from the flow's expressiveness or from the residual discretisation. The nonlinear-drift and high-dimensional groups are the real targets.

> [!warning] How far the numerical results were verified
> The composition of the three categories is verified; **the specific dimensions and error magnitudes are not**, nor is the exact schedule for adaptively refreshing the collocation set, so this page gives no numbers.

### Relation to the others

Paper 64 is the time-dependent bridge between the stationary ADDA-KRnet framework and the later flow papers. Paper 72 then treats **fractional** time-dependent Fokker-Planck equations with a very similar time-dependent coupling layer (compare paper 72's MCTNF layer, which additionally multiplies the $\tanh$ arguments by $t$). Paper 87 replaces the unbounded KRnet with a bounded one, and [[en/computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning|paper 107]] replaces "one equation, one flow" with an **operator over initial conditions**. The temporal normalising flow idea itself is credited by the paper to Both and Kusters.

## 72: two ways to handle a fractional operator

### The idea

Paper 72 pushes this line to fractional Fokker-Planck equations, that is, to densities of particles driven by Lévy noise as well as Gaussian noise. Here **three difficulties appear at once**: the solution lives on an unbounded domain, the dimension may be high, and the fractional Laplacian is non-local. Mesh methods fail on all three.

The flow disposes of the first two at a stroke: a pushforward density is automatically defined and normalised on all of $\mathbb R^d$, and dimension affects only network width rather than producing a grid explosion. The third needs separate treatment, and the paper offers two routes whose trade-off is the most memorable part of the work.

**Route one accepts the randomness.** Write the fractional Laplacian as an expectation and estimate it from a few samples. The cost is dimension-independent, but every evaluation carries noise.

**Route two avoids the randomness.** The fractional Laplacian of a Gaussian has a closed form (a confluent hypergeometric function), so switch to an auxiliary model built from a mixture of Gaussian radial basis functions, evaluate the fractional Laplacian exactly on it, and tie the auxiliary model to the flow with a consistency penalty. Evaluation is noise-free, but the number of basis centres grows with dimension.

One route is dimension-friendly but noisy, the other noise-free but dimension-unfriendly, and they complement each other exactly. This manoeuvre of **moving a hard operator onto a model where the operator is analytic** has no counterpart elsewhere in this topic and is the paper's most distinctive idea.

### Setting

The SDE with Lévy noise and the corresponding fractional Fokker-Planck equation are

$$
\mathrm d X_t=\mu(X_t,t)\,\mathrm dt+\sigma(X_t,t)\,\mathrm dW_t+\mathrm dL^\alpha_t,
\qquad
\frac{\partial p}{\partial t}=\mathcal L p-(-\Delta)^{\alpha/2}p,
$$

$$
\mathcal L p=-\nabla\cdot(p\mu)+\tfrac12\nabla\cdot\nabla\cdot(\sigma\sigma^{\mathsf T}p),
\qquad
(-\Delta)^{\alpha/2}p=C_{d,\alpha}\,\mathrm{P.V.}\!
\int_{\mathbb R^d\setminus\{0\}}\frac{p(x)-p(y)}{|x-y|_2^{d+\alpha}}\,\mathrm dy .
$$

The density model is a flow, $p_X(x)=p_Z(f(x))\,|\det\nabla_x f(x)|$, with $f$ a "simplified KRnet" composed of Actnorm layers

$$
y_{[i]}=a_i\odot x_{[i]}+b_i
$$

and affine coupling layers

$$
x_{[i],1}=x_{[i-1],1},
\qquad
x_{[i],2}=x_{[i-1],2}\odot\Bigl(\bm 1_{d-m}+\beta\tanh\bigl(s_i(x_{[i-1],1})\bigr)\Bigr)
+e^{\zeta_i}\odot\tanh\bigl(q_i(x_{[i-1],1})\bigr).
$$

> [!warning] A notation clash
> In the preprint rendering this coupling-layer equation writes the $\tanh$ scale bound as $\alpha$, but $\alpha$ is the fractional exponent everywhere else in the paper. Read it as $\beta$, matching the corresponding equation in paper 64 and the original KRnet layer.

### Derivation

**Route one: Monte Carlo (MCNF).** The loss is the mean-square fractional Fokker-Planck residual:

$$
L(p_{\mathrm{KRnet},\theta})=\frac{1}{N_S}\sum_{i=1}^{N_S}\bigl|R_\theta(x^i)\bigr|^2,
\qquad
R_\theta(x)=\bigl(\mathcal L-(-\Delta)^{\alpha/2}\bigr)p_{\mathrm{KRnet},\theta}(x).
$$

The fractional Laplacian inside it comes from the inner-outer split estimator of [[en/computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|paper 66]], appearing here as Lemma 3.1:

$$
(-\Delta)^{\alpha/2}u(\bm x)=C_{d,\alpha}\frac{|S^{d-1}|r_0^{2-\alpha}}{2(2-\alpha)}\,
\mathbb E_{\bm\xi\sim U(S^{d-1}),\,r_1\sim f_{\mathrm I}}\!
\left[\frac{2u(\bm x)-u(\bm x-r_1\bm\xi)-u(\bm x+r_1\bm\xi)}{r_1^2}\right]
$$

$$
\qquad+\;C_{d,\alpha}\frac{|S^{d-1}|r_0^{-\alpha}}{2\alpha}\,
\mathbb E_{\bm\eta\sim U(S^{d-1}),\,r_2\sim f_{\mathrm O}}
\bigl[2u(\bm x)-u(\bm x-r_2\bm\eta)-u(\bm x+r_2\bm\eta)\bigr],
$$

$$
f_{\mathrm I}(r)=\frac{2-\alpha}{r_0^{2-\alpha}}r^{1-\alpha}\mathbf 1_{r\in[0,r_0]},
\qquad
f_{\mathrm O}(r)=\alpha r_0^\alpha r^{-1-\alpha}\mathbf 1_{r\in[r_0,\infty)},
$$

$$
r_1/r_0\sim\mathrm{Beta}(2-\alpha,1),
\qquad
r_0/r_2\sim\mathrm{Beta}(\alpha,1),
\qquad
r_\epsilon=\max\{\epsilon,r_1\}.
$$

What it acts on differs from paper 66: there the target was a generic PINN surrogate, here it is a **density**, so non-negativity and normalisation are guaranteed by the flow and need no extra constraint.

One observation is easy to overlook and matters in practice: $\mathrm{Beta}(a,1)$ concentrates at the origin as $a\to0$, and the inner radius obeys $r_1/r_0\sim\mathrm{Beta}(2-\alpha,1)$, so as $\alpha$ approaches $2$ the first parameter $2-\alpha$ approaches zero, the samples $r_1$ pile up near zero, and $1/r_1^2$ amplifies cancellation error more readily — a **larger** floor $r_\epsilon$ is then needed for stability. This turns "parameter values affect numerical stability" into an actionable rule rather than a generic warning.

> [!note] A small difference from paper 66
> Paper 72 uses two **independent** direction variables $\bm\xi$ (inner) and $\bm\eta$ (outer), whereas the corresponding formula in paper 66 reuses a single $\xi$. Small, but real.

**Route two: an analytic auxiliary model (GRBFNF).** The starting point is Lemma 3.2: for a Gaussian $u(\bm x)=\exp(-\sigma^{-2}|\bm x-\bm x_0|_2^2)$,

$$
(-\Delta)^{\alpha/2}u(\bm x)=c_{\alpha,d}\,|\sigma|^{-\alpha}\;
{}_1F_1\!\Bigl(\frac{d+\alpha}{2};\frac{d}{2};-\sigma^{-2}|\bm x-\bm x_0|^2_2\Bigr),
\qquad
c_{\alpha,d}=\frac{2^\alpha\Gamma\bigl(\frac{d+\alpha}{2}\bigr)}{\Gamma\bigl(\frac d2\bigr)},
$$

with ${}_1F_1$ the confluent hypergeometric function. A Gaussian radial-basis mixture is then introduced as an auxiliary model, with centres $S_{\text{center}}=\{\tilde{\bm x}_i\}_{i=1}^M$:

$$
p_{\mathrm{GRBF},\tilde\theta}(\bm x)=\sum_{i=1}^M w_i\,\mathcal N(\tilde{\bm x}_i,\sigma_i^2\mathbf I)(\bm x),
\qquad 0\le w_i\le1,\ \sum_i w_i=1,
$$

$$
\mathcal N(\tilde{\bm x}_i,\sigma_i^2\mathbf I)(\bm x)
=(2\pi)^{-d/2}\sigma_i^{-d}\exp\!\Bigl(-\frac{|\bm x-\tilde{\bm x}_i|_2^2}{2\sigma_i^2}\Bigr),
$$

with both $w_i$ and $\sigma_i$ trainable. Applying Lemma 3.2 term by term and collecting constants gives a closed-form fractional Laplacian:

$$
(-\Delta)^{\alpha/2}p_{\mathrm{GRBF},\tilde\theta}(\bm x)
=c_{\alpha,d}\,\pi^{-d/2}2^{-\frac{d+\alpha}{2}}\sum_{i=1}^M w_i\,|\sigma_i|^{-(d+\alpha)}\,
{}_1F_1\!\Bigl(\frac{d+\alpha}{2};\frac d2;-\frac{|\bm x-\tilde{\bm x}_i|_2^2}{2\sigma_i^2}\Bigr).
$$

The coupled loss ties the flow to the auxiliary model. Its first term puts the **local part of the flow** and the **non-local part of the auxiliary model** into one residual; the second is a consistency penalty:

$$
\tilde L\bigl(p_{\mathrm{KRnet},\theta},p_{\mathrm{GRBF},\tilde\theta}\bigr)
=\frac{1}{N_S}\sum_{i=1}^{N_S}\Bigl(\mathcal L p_{\mathrm{KRnet},\theta}(x^i)
-(-\Delta)^{\alpha/2}p_{\mathrm{GRBF},\tilde\theta}(x^i)\Bigr)^2
+\frac{\beta_m}{N_S}\sum_{i=1}^{N_S}\Bigl(p_{\mathrm{KRnet},\theta}(x^i)
-p_{\mathrm{GRBF},\tilde\theta}(x^i)\Bigr)^2 .
$$

The construction of the first term deserves a pause: without the second term, the two models could drift onto unrelated functions while the first term stayed small. The consistency penalty $\beta_m$ is what turns "the auxiliary model must represent the same density" into an optimisable quantity. This is the paper's most distinctive idea: **move a hard operator onto a model where the operator is known, then pull the two representations together with a consistency penalty.**

**A detail in the time-dependent version (MCTNF).** Its coupling layer is

$$
x_{[i],2}=x_{[i-1],2}\odot\Bigl(\bm 1_{d-m}
+\beta\tanh\bigl(t\,s_{i,t}(x_{[i-1],1},t)\bigr)\Bigr)
+e^{\zeta_i}\odot\tanh\bigl(t\,q_{i,t}(x_{[i-1],1},t)\bigr).
$$

Note the explicit factor $t$ in front of the network outputs. Since $\tanh(0)=0$, at $t=0$ the scale factor becomes $\bm 1$ and the shift becomes $\bm 0$, so the whole layer degenerates to the identity and the flow at the initial instant is the reference distribution itself. This means **the initial condition is satisfied by construction**, not approximated by a penalty. Writing constraints into the architecture is the same instinct that makes the pushforward formula guarantee normalisation.

### Theorems

The paper's formal results are two lemmas, neither a convergence statement about networks but both identities about operator evaluation:

**Lemma 3.1 (Monte Carlo representation of the fractional Laplacian).** As above, under a splitting radius $r_0>0$, $0<\alpha<2$, and the two Beta sampling laws; $r_\epsilon$ is only a numerical safeguard and does not enter the identity itself.

**Lemma 3.2 (analytic fractional Laplacian of a Gaussian).** As above, for functions of the form $\exp(-\sigma^{-2}|\bm x-\bm x_0|^2)$. The paper attributes it to a reference rather than proving it.

Two sentences frequently repeated beyond the abstract — that "MCNF and GRBFNF improve accuracy by at least an order of magnitude over non-adaptive methods", and that "GRBFNF is better in low dimension while MCNF is better in high dimension" — **could not be confirmed in wording or magnitude in the accessible text** by the verification material behind this page. The second is a natural reading of the two constructions (basis-centre counts degrade with $d$, Monte Carlo does not), but it should be treated as a hypothesis to check rather than as the paper's conclusion.

### Numerical experiments

The adaptive strategy itself is short (stated for MCNF; GRBFNF is analogous):

1. Draw an initial training set from a simple distribution, uniform or Gaussian;
2. Train the flow on the current training set by minimising the residual loss (MCNF with the Monte Carlo fractional Laplacian; GRBFNF with the analytic fractional Laplacian of the auxiliary model plus the consistency term);
3. Draw a **new** training set by pushing reference samples through the inverse flow, that is, sample from the current approximate density;
4. Continue training on the refreshed set, alternating steps 2 and 3 until the stopping criterion.

This is exactly the "refine the training set and the approximate solution alternately" of the abstract. What distinguishes it from the [[en/computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|paper 70]] line is the quantity it follows: there, the residual failure probability; here, **the solution density itself**.

The examples cover four classes:

| Example                                          | What it tests                                       |
| ------------------------------------------------ | ---------------------------------------------------- |
| 2D fractional FPE driven by the fractional Laplacian alone | the non-local term in isolation from drift and diffusion |
| a bimodal target distribution                    | flow expressiveness, and whether adaptivity collapses onto one mode |
| higher-dimensional stationary FPEs               | where the two routes diverge with dimension          |
| time-dependent fractional FPE with a Cauchy solution | time conditioning and heavy tails; Cauchy is the natural stationary law for $\alpha$-stable noise |

The last choice is deliberate: the Cauchy distribution is heavy-tailed with no second moment, exactly the object a finite-domain grid method handles worst, whereas a flow is defined on all of $\mathbb R^d$ and heavy tails are just a shape the reference distribution takes after mapping.

> [!warning] How far the numerical results were verified
> The composition of the four example classes is verified; **the specific dimensions (4D/6D/8D) and the error magnitudes are not**, so this page gives no numbers.

### Relation to the others

Paper 72 equals the Monte Carlo fractional operator of [[en/computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|paper 66]] times the ADDA/KRnet density framework of paper 64 and Tang-Wan-Liao. Its coupling layer is in the same family as paper 64's and its adaptive refresh is the same idea as ADDA's. Paper 87 then replaces the unbounded KRnet with a B-KRnet suited to bounded and mixed domains, while [[en/computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning|paper 107]] replaces the goal of "solve one Fokker-Planck equation" with "learn a Fokker-Planck solution operator over initial conditions".

> [!note] Author order
> The homepage lists the authors as Xiaoliang Wan, Li Zeng, Tao Zhou; both the preprint and the journal version give **Li Zeng, Xiaoliang Wan, Tao Zhou**. This site follows the published version.

## 87: bounded support needs a different coupling layer

### The idea

KRnet and normalising flows generally map to a Gaussian reference on $\mathbb R^d$, so the induced density has support on the whole space. For a density genuinely supported on a bounded set that is wrong — the paper's everyday example is human age, where a model should not assign positive probability to $-3$ or $700$ years — and it is equally wrong for a PDE solution posed on a hyperrectangle.

The obvious remedy, restricting the flow to a box by force, destroys exact invertibility: a point mapped outside the box has no preimage, the inverse map is no longer everywhere defined, and the pushforward formula fails with it.

The paper instead **reasons backwards from the requirement**. What is needed is a monotone bijection $[-1,1]\to[-1,1]$ whose inverse is explicitly computable. Which family satisfies that naturally? **Cumulative distribution functions.** The CDF of any strictly positive density on $[-1,1]$ maps $[-1,1]$ monotonically onto $[0,1]$, and an affine rescaling brings it back to $[-1,1]$. All that remains is to choose a density family that is expressive and explicitly invertible.

The paper chooses **piecewise linear** densities. Their CDF is piecewise quadratic, so inverting means solving a quadratic with a radical formula. This choice turns invertibility from a property requiring careful maintenance into one line of algebra.

The price is equally clear, and the paper puts it in its own conclusion: the model density is piecewise linear, hence **only first-order differentiable**. A second-order PDE cannot act on it directly and must first be recast as a first-order system.

### Setting

The flow density and its composite structure are

$$
p_{\bm X}(\bm x)=p_{\bm Z}\bigl(f(\bm x)\bigr)\,\bigl|\det\nabla_{\bm x}f(\bm x)\bigr|,
\qquad
f=f_{[L]}\circ\cdots\circ f_{[1]},
$$

$$
f^{-1}=f^{-1}_{[1]}\circ\cdots\circ f^{-1}_{[L]},
\qquad
\bigl|\det\nabla_{\bm x}f\bigr|=\prod_i\bigl|\det\nabla_{\bm x_{[i-1]}}f_{[i]}\bigr| .
$$

The key difference is the reference: $\bm Z$ is **uniform on $[-1,1]^d$** rather than Gaussian on $\mathbb R^d$.

The structural justification is the Knothe-Rosenblatt rearrangement (the paper's Prop. 2.2): there is a triangular map $T$ with $T_\#\rho=\pi$,

$$
\bm z=T(\bm x)=\bigl(T_1(x_1),\,T_2(x_1,x_2),\dots,T_d(x_1,\dots,x_d)\bigr)^{\mathsf T},
$$

which factors as a composition of $d$ maps each updating a single coordinate.

### Derivation

**Pseudo-triangular (block) structure.** Partition $\bm x=\bigl((\bm x^{(1)})^{\mathsf T},\dots,(\bm x^{(K)})^{\mathsf T}\bigr)^{\mathsf T}$ with $\sum_i d_i=d$ and relax the per-coordinate triangular structure to a per-block one:

$$
\bm z=f_{\text{KR}}(\bm x)=
\begin{pmatrix}\tilde f_1\\ \bm x^{(2:K)}\end{pmatrix}\circ
\begin{pmatrix}\tilde f_2\\ \bm x^{(3:K)}\end{pmatrix}\circ\cdots\circ
\begin{pmatrix}\tilde f_{K-1}\\ \bm x^{(K)}\end{pmatrix}\circ \tilde f_K(\bm x),
$$

with $\tilde f_k:[-1,1]^{\sum_{i\le k}d_i}\to[-1,1]^{\sum_{i\le k}d_i}$. The outer loop has $K$ stages, $f^{outer}_{[k]}=\tilde f_{K-k+1}$, and one block is frozen after each stage (the paper calls it a squeezing operation that deactivates some dimensions); the inner loop composes several CDF coupling layers, $f^{outer}_{[k]}=L_{\text{CDF},[k,l_k]}\circ\cdots\circ L_{\text{CDF},[k,1]}$.

**The CDF coupling layer.** On a mesh $-1=s_0<s_1<\cdots<s_n=1$, define a piecewise linear density

$$
p(s)=\frac{w_{i+1}-w_i}{h_i}(s-s_i)+w_i,
\quad s\in[s_i,s_{i+1}],
\quad p(s_i)=w_i\ge0,\ h_i=s_{i+1}-s_i .
$$

Integrating piece by piece — the linear term gives a quadratic, the constant term gives a linear one, and earlier pieces contribute accumulated trapezoid areas — gives a **piecewise quadratic** cumulative distribution function:

$$
F(s)=\frac{w_{i+1}-w_i}{2h_i}(s-s_i)^2+w_i(s-s_i)
+\sum_{k=0}^{i-1}\frac{w_k+w_{k+1}}{2}(s_{k+1}-s_k),
\quad s\in[s_i,s_{i+1}].
$$

Inverting means solving one quadratic per piece. With $q_0=0$ and $q_i=\sum_{k=0}^{i-1}\frac{w_k+w_{k+1}}{2}(s_{k+1}-s_k)$,

$$
F^{-1}(q)=s_i+\frac{-w_i+\sqrt{w_i^2+2(w_{i+1}-w_i)(q-q_i)/h_i}}{(w_{i+1}-w_i)/h_i}
=s_i+\frac{2(q-q_i)}{w_i+\sqrt{w_i^2+2(w_{i+1}-w_i)(q-q_i)/h_i}} .
$$

The two forms are algebraically equivalent, the second obtained by rationalising the numerator. It is the numerically stable one when $w_{i+1}\approx w_i$ (a nearly constant piece, where the first form's denominator tends to zero) — an implementation necessity, not a stylistic preference.

Setting $\tilde F(s)=2F(s)-1$ gives a bijection $[-1,1]\to[-1,1]$, applied componentwise as $F_i(y_i)=\tilde F(y_i;\bm\theta_i)$. For $\bm y=(\bm y_1^{\mathsf T},\bm y_2^{\mathsf T})^{\mathsf T}\in[-1,1]^l$ with $\bm y_1\in\mathbb R^m$, the coupling layer is

$$
\begin{pmatrix}\hat{\bm y}_1\\ \hat{\bm y}_2\end{pmatrix}
=L_{\text{CDF},\bm\theta}\begin{pmatrix}\bm y_1\\ \bm y_2\end{pmatrix}:
\qquad
\hat{\bm y}_1=\bm y_1,
\qquad
\hat{\bm y}_2=\bm F\bigl(\bm y_2;\bm\theta(\bm y_1)\bigr),
$$

where $\bm\theta=(\bm s_1^{\mathsf T},\dots,\bm s_{n-1}^{\mathsf T},\bm w_0^{\mathsf T},\dots,\bm w_n^{\mathsf T})^{\mathsf T}\in\mathbb R^{2n(l-m)}$ is produced by a network $\mathrm{NN}(\bm y_1)$ and reparameterised to keep the mesh strictly increasing and the weights non-negative:

$$
\bm s_1=-1+\frac2n\bigl(1+\beta_1\tanh(\hat{\bm s}_1)\bigr),
\qquad
\bm s_{i+1}=\bm s_i+\Bigl(\frac{1-\bm s_i}{n-i}\Bigr)\bigl(1+\beta_{i+1}\tanh(\hat{\bm s}_{i+1})\bigr),
\qquad
\bm w_i=\frac{1+\gamma_i\tanh(\hat{\bm w}_i)}{C},
$$

with $C$ the normalisation constant and, in practice, $\beta_1=\tfrac{65}{66}$, $\beta_i=0.97$ for $i\ge2$, and $\gamma_i=0.99$. All of these are strictly below $1$, which clamps $1+\beta\tanh(\cdot)$ inside $(1-\beta,1+\beta)\subset(0,2)$: mesh increments stay positive so the $s_i$ increase strictly, and the weights $w_i$ stay positive so the density does not degenerate. A second layer with the roles of $\bm y_1,\bm y_2$ swapped follows, and layers alternate by parity of the index so that every coordinate gets updated.

Degrees of freedom accumulate layerwise, $\mathrm{DOFs}(p_{\text{B-KRnet},\bm\theta})=\sum_{k=1}^{K}\sum_{i=1}^{l_k}\mathrm{DOFs}(L_{\text{CDF},[k,i]})$, with $\mathrm{DOFs}(f^{outer}_{[K]})=0$ when $d_1=1$.

**An edge case (Remark 2.3).** A CDF coupling layer needs input dimension at least $2$ (otherwise there is no $\bm y_1$ to produce the parameters), so when $\bm x^{(1)}=x_1\in\mathbb R$ the map $f_{\text{KR}}$ consists of only $K-1$ bijections and the outermost $f^{outer}_{[K]}$ is taken to be the identity. The paper's justification is that $\bm x^{(1)}_{[K-1]}$ is by then already close to uniform, so $2F(\cdot)-1$ is close to the identity anyway and little is lost.

### Four differences from KRnet

| Item             | KRnet                     | B-KRnet                            |
| ---------------- | ------------------------- | ---------------------------------- |
| Reference        | Gaussian on $\mathbb R^d$ | uniform on $[-1,1]^d$              |
| Coupling layer   | affine (scale and shift)  | monotone piecewise-quadratic CDF   |
| Scale-bias layer | present                   | **removed** (the domain is bounded) |
| Smoothness       | set by the activation     | piecewise linear density, first-order differentiable only |

The second row is the essential one: an affine map sends $[-1,1]^m$ elsewhere, a CDF map keeps it invariant. The third follows from the second — the scale-bias layer exists to bring data to a sensible range, and on a bounded domain the domain itself does that. The fourth has a direct consequence: **second-order PDEs must be recast as first-order systems.**

A general product domain $\prod_i[a_i,b_i]$ is normalised to $[-1,1]^d$ by the linear map $\bm y=\hat{\bm a}\odot\bm x+\hat{\bm b}$ with $\hat a_i=2/(b_i-a_i)$ and $\hat b_i=-(b_i+a_i)/(b_i-a_i)$.

### Theorems

**Prop. 2.4 (weak convergence of pushforwards).** If $f_n\to T$ in $L^p([-1,1]^d)$ with $T_\#\rho=\pi$, then $(f_n)_\#\rho\rightharpoonup T_\#\rho$ and $(f_n^{-1})_\#\pi\rightharpoonup(T^{-1})_\#\pi$. The paper gives a proof sketch.

**Prop. 2.6 (universality on bounded domains).** For any absolutely continuous $\rho$ on a bounded $\Omega$ there is a sequence of B-KRnet transformations $f_k$ with $(f_k)_\#\rho$ converging weakly to the uniform measure on $\Omega$. The paper states this and defers the proof to another article.

The paper also notes that in the limit $l_k=1$, $n\to\infty$ the CDF coupling layer is a direct approximation of the Knothe-Rosenblatt rearrangement — which pins down the relation between architecture and classical construction: not a heuristic analogy but an approximation with a limit behind it.

### Numerical experiments

**Solving PDEs with B-KRnet.** For a second-order operator $\mathcal N[\bm x;\{\partial^{\bm\alpha}_{\bm x}p,|\bm\alpha|\le2\}]=0$ with $\int_\Omega p=1$ and $p\ge0$ (or with mass $M$ absorbed into a trainable scalar $\zeta$ so that $\zeta\hat p\approx p$), introduce $\bm g(\bm x)=\nabla p$ approximated by a separate network $\bm g_{\text{NN}}$, reducing to a first-order system, and minimise

$$
\mathcal L\bigl(p_{\text{B-KRnet},\bm\theta},\bm g_{\text{NN}}\bigr)
=\lambda_{pde}\mathcal L_{pde}+\lambda_b\mathcal L_b+\lambda_{\bm g}\mathcal L_{\bm g},
$$

$$
\mathcal L_{pde}=\mathbb E_{\bm x\sim\rho}\bigl[\bigl|\mathcal N[\bm x;p_{\text{B-KRnet},\bm\theta},\bm g_{\text{NN}}]\bigr|^2\bigr],
\qquad
\mathcal L_b=\mathbb E_{\bm x\sim\rho_b}\bigl[|\mathcal B[\cdots]|^2\bigr],
$$

$$
\mathcal L_{\bm g}=\mathbb E_{\bm x\sim\rho}
\bigl[\|\bm g_{\text{NN}}(\bm x)-\nabla p_{\text{B-KRnet},\bm\theta}(\bm x)\|_2^2\bigr].
$$

The third term is a consistency penalty pinning $\bm g_{\text{NN}}$ to the true gradient, the same kind of construction as the term that pins the GRBF model to the flow in paper 72. Non-negativity and mass conservation hold **by construction** and never enter the penalties — which is the entire reason for representing the density by a flow.

**Mixed domains.** The solution of a kinetic Fokker-Planck equation is a joint density of position $\bm x\in\Omega$ (bounded) and velocity $\bm v\in\mathbb R^d$ (unbounded). The paper factors $p(\bm x,\bm v)=h(\bm v|\bm x)q(\bm x)$ and models

$$
p_{\bm\theta}(\bm x,\bm v)=h_{\text{KRnet},\bm\theta_1}(\bm v|\bm x)\cdot q_{\text{B-KRnet},\bm\theta_2}(\bm x),
$$

that is, the unbounded KRnet as a conditional flow for velocity times B-KRnet as the density of position. The conditioning mechanism here is exactly the time conditioning of paper 64 with position substituted for time — the same technical part playing different roles in different papers.

**The adaptive strategy** (Section 4.2's four steps plus Algorithm 1):

1. Set $\rho(\bm x)=1/|\Omega|$ and draw $C^0_{pde}=\{\bm x^{i,0}\}_{i=1}^{N_{pde}}\sim\mathrm{Unif}(\Omega)$ and $C_b=\{\bm x^i_b\}_{i=1}^{N_b}\sim\mathrm{Unif}(\partial\Omega)$;
2. Train B-KRnet on $C^0_{pde}$, $C_b$ to obtain $\bm\theta^{*,0}$;
3. Update $\rho(\bm x)\leftarrow(1-\gamma)\rho(\bm x)+\gamma\,p_{\text{B-KRnet},\bm\theta^{*,0}}(\bm x)$ with $N_{new}=\gamma N_{pde}$, generating the new points by the inverse flow, $\bm z^{i,1}\sim\mathrm{Unif}(-1,1)^d$ and $\bm x^{i,1}=(f_{\bm\theta^{*,0}})^{-1}(\bm z^{i,1})$, keeping the remaining points;
4. Repeat steps 2 and 3 a total of $N_{\mathrm{adaptive}}$ times.

Algorithm 1 wraps this with mini-batching, Adam, and a step-decay learning-rate schedule ($l_r\leftarrow\eta\,l_r$ every $n_s$ steps). The mixture form in step 3 is worth noticing: the new sampling density is not simply replaced by the current flow but blended with the old density at rate $\gamma$ — the same class of safety valve as the density floor in [[en/computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|paper 80]] and the cosine-annealed uniform-to-adaptive transition in paper 73.

**The example list.** Three density-estimation examples (Section 3) and three PDE examples (Section 6):

| Section | Problem                                                  | What it tests                        |
| ------- | -------------------------------------------------------- | ------------------------------------- |
| 3.1     | density estimation on an annulus                          | non-convex support                    |
| 3.2     | a mixture of Gaussians                                    | multimodal expressiveness             |
| 3.3     | a logistic distribution with holes                        | disconnected support                  |
| 6.1     | a four-dimensional $-\Delta p+p=f$                        | whether the first-order recast works at moderate dimension |
| 6.2     | the stationary Keller-Segel system (two coupled densities) | several densities at once             |
| 6.3     | a 2D stationary kinetic Fokker-Planck equation, KRnet ⊗ B-KRnet | mixing bounded and unbounded dimensions |

The shared settings are: `Tanh` activation; **three subintervals per CDF coupling layer**; $\mathrm{NN}(\bm y_1)$ a fully connected network with two hidden layers; Adam in PyTorch; **update rate $\gamma=0.8$**.

The qualitative claim of the abstract and conclusion is that samples generated by B-KRnet show excellent agreement with the ground truth and that the adaptive density-approximation scheme is effective on the four-dimensional problem, the Keller-Segel equation and the kinetic Fokker-Planck equation. The paper concedes two limitations: B-KRnet's smoothness stops at first order and improving it would need higher-order polynomials; and adaptive point updating inside a deep Ritz framework is left open — the latter picked up by [[en/computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|paper 80]].

> [!warning] How far the numerical results were verified
> The six examples, the shared settings and both conceded limitations are verified; **the specific error values are not**, so this page gives no numbers.

> [!note] Version differences
> The v1 abstract of the preprint says B-KRnet "consists of a series of coupling layers with progressively fewer active transformation dimensions, inspired by the triangular structure of the Knothe-Rosenblatt rearrangement" and calls the reference measure a base distribution; v3 says it "adapts the pseudo-triangular structure into a normalizing flow model" and calls it a prior distribution. Note which version you cite.

### Relation to the others

B-KRnet is the bounded member of the KRnet family that runs through this whole topic. It is the density model used by [[en/computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|paper 80]] (adaptive importance sampling for deep Ritz). Its conditional-flow trick for mixed domains comes from paper 64, with position substituted for time as the conditioning variable. It is the bounded counterpart of the KRnet used in paper 72 and in the original ADDA-FPE work of Tang-Wan-Liao, and [[en/computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning|paper 107]] generalises the conditional-flow-for-Fokker-Planck idea into a whole operator over initial conditions.

## How the four build on one another

| No. | Target object                        | Reference                    | Training signal                            | New technical part                      |
| --- | ------------------------------------ | ---------------------------- | ------------------------------------------ | --------------------------------------- |
| 62  | random field from scattered sensors  | Gaussian field with KL structure | log-likelihood (plus optional physics residual) | KL expansion with network coefficients |
| 64  | time-dependent Fokker-Planck         | unbounded Gaussian           | equation residual                          | latent time pinned to real time         |
| 72  | fractional Fokker-Planck             | unbounded Gaussian           | equation residual                          | Monte Carlo or analytic auxiliary model |
| 87  | densities with bounded support       | uniform on $[-1,1]^d$        | equation residual (recast as a first-order system) | piecewise-quadratic CDF coupling layer |

One judgement runs through all four: **structural constraints should be guaranteed by the architecture, not approximated by penalties.** Normalisation is guaranteed by change of variables, non-negativity by the pushforward formula, the initial condition by multiplying the coupling layer by $t$, and bounded support by the CDF coupling layer. Each such guarantee removes one penalty term whose weight would otherwise need tuning.

A second judgement concerns the price: every constructive guarantee is bought with expressiveness or smoothness. The factor $t$ restricts how much the layer can vary near $t=0$; the CDF coupling layer buys exact invertibility but flattens the model density into a piecewise linear one, so second-order equations have to be reduced in order. **Structure is not free.**

## Coverage check

| Item                                             | Paper | Status                                                    |
| ------------------------------------------------ | ----- | ---------------------------------------------------------- |
| Division of labour between KL structure and flow | 62    | intuition, three-step construction, three conceded advantages |
| Numerical experiments of 62                      | 62    | three groups; **formulas and error magnitudes unverified**  |
| Time conditioning and Jacobian collapse          | 64    | three candidate routes, the wrong one, the right one, consequences |
| Time-dependent affine coupling layer             | 64    | formula, roles of $\beta$ and the $t$ input, two differences from real NVP |
| Numerical experiments of 64                      | 64    | four-step algorithm and three categories; **dimensions and errors unverified** |
| Fractional FPE under Lévy noise                  | 72    | equation, density model, coupling layer and the notation clash |
| Lemma 3.1 and the Beta-concentration rule        | 72    | full representation, two independent directions, $r_\epsilon$ dependence |
| Lemma 3.2, GRBF and the coupled loss             | 72    | closed form, mixture model, necessity of the consistency penalty |
| Multiplying by $t$ for the initial condition     | 72    | coupling-layer form and its meaning                        |
| Numerical experiments of 72                      | 72    | four adaptive steps, four example classes; **dimensions and errors unverified**; two often-repeated conclusions marked as unconfirmed |
| Piecewise linear density and quadratic CDF       | 87    | piecewise integration, the stable inverse, reparameterisation constants and their role |
| Pseudo-triangular structure, freezing, $d_1=1$   | 87    | block decomposition, squeezing, the reason for the edge case |
| Prop. 2.4 and Prop. 2.6                          | 87    | both statements and their proof status                     |
| First-order recast and the three-term loss       | 87    | introducing $\bm g$, the three terms, the two properties held by construction |
| Numerical experiments of 87                      | 87    | six examples, shared settings, four adaptive steps, two conceded limitations; **errors unverified** |

## Sources for this page

- L. Guo, H. Wu, and T. Zhou, _Normalizing field flows: solving forward and inverse stochastic differential equations using physics-informed flow models_, J. Comput. Phys. 461 (2022), 111202 (preprint [arXiv:2108.12956](https://arxiv.org/abs/2108.12956)).
- X. Feng, L. Zeng, and T. Zhou, [_Solving time dependent Fokker-Planck equations via temporal normalizing flow_](https://doi.org/10.4208/cicp.OA-2022-0090), Commun. Comput. Phys. 32(2) (2022), pp. 401-423 (preprint [arXiv:2112.14012](https://arxiv.org/abs/2112.14012)).
- L. Zeng, X. Wan, and T. Zhou, [_Adaptive deep density approximation for fractional Fokker-Planck equations_](https://doi.org/10.1007/s10915-023-02379-z), J. Sci. Comput. 97 (2023), 68 (preprint [arXiv:2210.14402](https://arxiv.org/abs/2210.14402)).
- L. Zeng, X. Wan, and T. Zhou, _Bounded KRnet and its applications to density estimation and approximation_, SIAM J. Sci. Comput. 47(6) (2025), pp. C1294-C1318 (preprint [arXiv:2305.09063](https://arxiv.org/abs/2305.09063)).
