---
title: Where to Put the Samples - Failure Probability and Importance Sampling
description: Papers 66, 70, 73, 76 and 80 - collocation design as reliability analysis and variance reduction
lang: en
translation: computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns
tags:
  - paper-notes
  - scientific-machine-learning
  - adaptive-sampling
---

![Let reliability analysis choose the next collocation batch](assets/diagrams/tao-zhou-papers/en/failure-informed-sampling.svg)

The accuracy of a physics-informed neural network depends strongly on where the collocation points sit. The shared judgement across these papers is that **point placement is not an implementation detail but an object that can be driven by an explicit mathematical target**. The four targets differ — the probability that the residual is too large, the variance of a variational loss, the unbiasedness of a conditional expectation — but each reduces to the same recipe: define an estimable scalar, then generate points from it.

## 66: write the operator itself as an expectation

### The idea

The fractional Laplacian is non-local: computing its value at a single point requires, in principle, knowledge of the solution over **all of space**. The earlier fPINN approach discretised that integral on an auxiliary grid, so every collocation point had to consult many neighbouring grid nodes, and both storage and work grew exponentially with dimension.

Paper 66 starts from a plain observation about the singular integral

$$
\mathrm{P.V.}\!\int_{\mathbb R^d}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\,\mathrm dy .
$$

It **already is a weighted average** — the weight $\|x-y\|^{-d-\alpha}$ simply is not normalised. Normalise it into a probability density and the integral becomes an expectation, and an expectation can be estimated without bias from a handful of random samples. The fractional operator then needs neither assembly nor storage: each collocation point only draws a few directions and a few radii and evaluates the network at a few shifted points. Dimension enters only through "draw a $d$-dimensional direction", and that step costs linearly.

There is a trap in this plan, and it is precisely the ingredient the abstract advertises. A physics-informed loss wants the **square** of the residual, and the square of a noisy estimate is biased: $\mathbb E[\widehat L^2]=(\mathbb E\widehat L)^2+\mathrm{Var}(\widehat L)$. Squaring directly means the optimiser is no longer minimising the true residual but the residual plus the estimator variance — a quantity unrelated to the PDE which can nevertheless be reduced by flattening the network's second differences. The paper's fix is to draw **two groups** of independent random numbers and multiply two independent estimates instead of squaring one.

### Setting

The model equation is a fractional advection-diffusion equation,

$$
L[u(x,t)]:=\frac{\partial^\gamma u(x,t)}{\partial t^\gamma}
+c\,(-\Delta)^{\alpha/2}u(x,t)+v\cdot\nabla u(x,t)=f(x,t),
$$

whose Caputo time derivative and integral fractional Laplacian are

$$
\frac{\partial^\gamma u(x,t)}{\partial t^\gamma}
\triangleq\frac{1}{\Gamma(1-\gamma)}\int_0^t (t-\tau)^{-\gamma}
\frac{\partial u(x,\tau)}{\partial\tau}\,\mathrm d\tau,
\qquad 0<\gamma<1,
$$

$$
(-\Delta)^{\alpha/2}u(x)\triangleq C_{d,\alpha}\,\mathrm{P.V.}\!
\int_{\mathbb R^d}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\,\mathrm dy .
$$

The training set splits into three parts $\mathcal D=(\mathcal D_f,\mathcal D_g,\mathcal D_u)$ — equation points, boundary or initial points, and observation points — with total loss

$$
\mathcal{LOSS}(\theta)=w_{equ}L_{equ}(\theta)+w_g L_g(\theta)+w_u L_u(\theta).
$$

### Derivation

**Inner-outer splitting of the fractional Laplacian.** Cut $\mathbb R^d$ into a ball of radius $r_0$ around $x$ and its complement, and treat the two pieces separately, because only the inner piece is singular.

Inside the ball, pass to polar coordinates $y=x+r\xi$ with $\mathrm dy=r^{d-1}\mathrm dr\,\mathrm dS(\xi)$, so that $\|x-y\|^{-d-\alpha}$ becomes $r^{-d-\alpha}$ and the two powers of $r$ cancel down to $r^{-1-\alpha}$. Then symmetrise under $\xi\mapsto-\xi$ (the sphere $S^{d-1}$ is invariant under this involution), which in the principal-value sense replaces the first difference by a second difference:

$$
u(x)-u(x+r\xi)\ \longrightarrow\ \tfrac12\bigl[2u(x)-u(x+r\xi)-u(x-r\xi)\bigr].
$$

This step is the technical core of the paper: **symmetrisation kills the $O(r)$ first-order term, and the remaining second difference is $O(r^2)$, exactly cancelling the residual singularity of $r^{-1-\alpha}$ when $\alpha<2$.** Factoring out $r^2$ leaves a radial weight $r^{1-\alpha}$, integrable on $[0,r_0]$ with normalising constant $\int_0^{r_0}r^{1-\alpha}\mathrm dr=r_0^{2-\alpha}/(2-\alpha)$. Hence

$$
\int_{y\in B_{r_0}(x)}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\mathrm dy
=\frac{|S^{d-1}|\,r_0^{2-\alpha}}{2(2-\alpha)}\;
\mathbb E_{\xi,\,r\sim f_I}\!\left[\frac{2u(x)-u(x-r\xi)-u(x+r\xi)}{r^2}\right],
$$

$$
f_I(r)=\frac{2-\alpha}{r_0^{2-\alpha}}\,r^{1-\alpha}\mathbf 1_{r\in[0,r_0]},
\qquad r/r_0\sim\mathrm{Beta}(2-\alpha,1).
$$

The Beta law can be read straight off the cumulative distribution function: the antiderivative of $f_I$ is $(r/r_0)^{2-\alpha}$, and the cumulative distribution function of $\mathrm{Beta}(a,1)$ is exactly $x^a$, so $r/r_0$ is sampled by a single power transform without rejection.

Outside the ball there is no singularity, so no division by $r^2$ is needed; the radial weight is $r^{-1-\alpha}$ with normalising constant $\int_{r_0}^\infty r^{-1-\alpha}\mathrm dr=r_0^{-\alpha}/\alpha$:

$$
\int_{y\notin B_{r_0}(x)}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\mathrm dy
=\frac{|S^{d-1}|\,r_0^{-\alpha}}{2\alpha}\;
\mathbb E_{\xi,\,r\sim f_O}\bigl[2u(x)-u(x-r\xi)-u(x+r\xi)\bigr],
$$

$$
f_O(r)=\alpha\,r_0^{\alpha}\,r^{-1-\alpha}\mathbf 1_{r\in[r_0,\infty)},
\qquad r_0/r\sim\mathrm{Beta}(\alpha,1).
$$

**Numerical safeguard.** The second difference quotient is well defined as $r\to0$,

$$
\lim_{r\to0}\frac{2u(x)-u(x-r\xi)-u(x+r\xi)}{r^2}
=\partial^2_r u(x+r\xi)\big|_{r=0},
$$

but in floating point the numerator subtracts two nearly equal quantities and $1/r^2$ amplifies the cancellation error. The implementation therefore replaces the radius by $r_\epsilon=\max\{\epsilon,r_I\}$. Assembled, the estimator on the network reads

$$
(-\Delta)^{\alpha/2}u_{NN}(x)=C_{d,\alpha}\frac{|S^{d-1}|r_0^{2-\alpha}}{2(2-\alpha)}
\mathbb E_{\xi,r_I\sim f_I}\!\left[
\frac{2u_{NN}(x)-u_{NN}(x-r_\epsilon\xi)-u_{NN}(x+r_\epsilon\xi)}{r_\epsilon^2}\right]
$$

$$
\qquad+\;C_{d,\alpha}\frac{|S^{d-1}|r_0^{-\alpha}}{2\alpha}\,
\mathbb E_{\xi,r_o\sim f_O}\bigl[2u_{NN}(x)-u_{NN}(x-r_o\xi)-u_{NN}(x+r_o\xi)\bigr].
$$

**The Caputo time derivative** is handled by the same manoeuvre, producing one expectation plus one explicit term:

$$
\frac{\partial^\gamma u_{NN}(x,t)}{\partial t^\gamma}
=\frac{\gamma}{1-\gamma}\,t^{1-\gamma}\,
\mathbb E_{\tau\sim f_{I,t}}\!\left[\frac{u_{NN}(x,t)-u_{NN}(x,t-\tau)}{\tau}\right]
+\frac{u_{NN}(x,t)-u_{NN}(x,0)}{t^{\gamma}},
$$

$$
f_{I,t}(\tau)=(1-\gamma)\tau^{-\gamma}\mathbf 1_{\tau\in[0,1]},
\qquad \tau\sim\mathrm{Beta}(1-\gamma,1),
\qquad \tau_\epsilon=\max\{\tau,\epsilon_t t^{-1}\}.
$$

The second term carries the initial value $u_{NN}(x,0)$ explicitly, consistent with the Caputo definition, which subtracts the initial value before differentiating fractionally. The three pieces — fractional Laplacian, Caputo derivative, and the advection term $v\cdot\nabla u_{NN}$ supplied by automatic differentiation — assemble into one random operator $\widehat L$.

**The unbiased construction.** Replace the square of $\widehat L$ by the product of two independent copies:

$$
\hat L_{equ}(\theta)=\frac{1}{mN_u}\sum_{i,j}
\widehat L\bigl[u_{NN}(x_i,t_i;\theta);\epsilon,\epsilon_t,\tau_j,\xi_j,r_{Ij},r_{oj}\bigr]\cdot
\widehat L\bigl[u_{NN}(x_i,t_i;\theta);\epsilon,\epsilon_t,\tau'_j,\xi'_j,r'_{Ij},r'_{oj}\bigr].
$$

The two groups of random parameters are independent, so the expectation of the product equals the product of the expectations, which is the true squared residual; the variance term disappears.

### Theorems

The paper proves no convergence theorem. Its theoretical content is the **exactness** of the two expectation identities above, together with a statement of unbiasedness: provided $\epsilon=\epsilon_t=0$ (that is, with the numerical safeguard switched off) and round-off is ignored,

$$
\mathbb E\bigl[\hat L_{equ}(\theta)\bigr]=L_{equ}(\theta).
$$

The condition is worth remembering: in practice $\epsilon,\epsilon_t>0$, so unbiasedness in the strict sense holds only in the limit, and $\epsilon$ doubles as a knob trading stability against bias.

### Numerical experiments

The algorithm itself is short (Algorithm 1 in the paper): given $\mathcal D_f,\mathcal D_g,\mathcal D_u$, draw $N$ snapshots from the training data at each step; for each residual point draw $x_i,t_i$ uniformly from $\Omega\subset\mathbb R^d$ and the time interval, then draw two independent groups $\{\tau_i,\xi_i,r_{Ii},r_{oi}\}_{i=1}^m$ and $\{\tau'_i,\xi'_i,r'_{Ii},r'_{oi}\}_{i=1}^m$; assemble the loss as above and update with Adam until convergence. Note the two layers of randomness: the randomness of the collocation points and the randomness of the operator estimate, the latter redrawn at every step.

The experiments fall into three groups:

| Group | Target                                                                        |
| ----- | ----------------------------------------------------------------------------- |
| One   | high-dimensional integral fractional Laplacian equations (forward problems)   |
| Two   | parameter identification for time-space fractional PDEs, recovering the order |
| Three | fractional diffusion equations with random inputs                             |

The paper's quantitative claim is that the **overall computational cost is lower than fPINN**, which is what opens high-dimensional fractional PDEs at all.

### Relation to the others

This inner-outer estimator is imported verbatim by paper 72 as its Lemma 3.1 (same $f_I$, same $f_O$, same Beta sampling, same $r_\epsilon$ safeguard); what differs is the object it acts on. Paper 66 applies it to a generic PINN surrogate $u_{NN}$, whereas paper 72 applies it to a **normalising-flow density**, which brings non-negativity and normalisation into play. Paper 72 also offers an alternative route that avoids randomness entirely, a Gaussian radial basis auxiliary model; see the [[en/computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|density-flow page]].

## 70: define "the residual is too large" as a failure event

### The idea

Standard residual-based adaptive refinement takes the $m$ points of largest residual from a uniform candidate pool. The scenario in which it fails is very concrete: when the high-residual region occupies a tiny fraction of the domain — a peak, an interior layer, a solution with concentrated support on an unbounded domain — **almost no uniform candidate lands inside it**, so "take the largest $m$" still returns a pile of irrelevant points.

Paper 70 borrows its tools from another discipline. Structural reliability analysis poses a standard problem: given a load distribution and a criterion for "the structure fails", estimate the failure probability, which is typically so small that naive Monte Carlo is useless and importance sampling must push the proposal density toward the failure region. Replace "the structure fails" by "the residual exceeds a tolerance" and the whole apparatus transfers word for word: limit-state function, failure set, self-adaptive importance sampling. **The by-product is a posterior error indicator** — the failure probability itself — whose status matches that of a posterior error estimate in adaptive finite elements, and which (see the theorem below) genuinely appears in the error bound.

### Setting

The physics-informed loss is

$$
\mathcal L(\theta)=\mathcal L_c(\theta)+\lambda\mathcal L_b(\theta),
\qquad
\mathcal L_c(\theta)=\frac{1}{N_c}\sum_{i=1}^{N_c}\bigl|r(x^c_i;\theta)\bigr|^2,
\qquad
\mathcal L_b(\theta)=\frac{1}{N_b}\sum_{i=1}^{N_b}\bigl|b(x^b_i;\theta)\bigr|^2 .
$$

A general limit-state function is written $g(x)=\mathcal Q(x)-\epsilon_r$, where $\mathcal Q$ maps the domain to some quantity of interest; for a PINN one takes $\mathcal Q(x)=|r(x;\theta)|$, so that

$$
g(x)=\bigl|r(x;\theta)\bigr|-\epsilon_r,
$$

whose zero level set splits the domain into a safe set $\Omega_{\mathcal S}=\{g<0\}$ and a failure set $\Omega_{\mathcal F}=\{g>0\}$. Under a prior $\omega(x)$, the failure probability

$$
P_{\mathcal F}=\int_\Omega \omega(x)\,\mathbb I_{\Omega_{\mathcal F}}(x)\,\mathrm dx
$$

serves as the posterior error indicator: the network is declared reliable once $P_{\mathcal F}<\epsilon_p$.

One concrete difference from residual-based refinement is worth naming: Monte Carlo enrichment adds **every** candidate point that falls in $\Omega_{\mathcal F}$, so the number $m$ added varies from round to round, whereas residual refinement adds a fixed $m$ each time.

### Derivation

The naive estimate $\hat P^{MC}_{\mathcal F}=\frac{1}{|\mathcal S|}\sum_{x\in\mathcal S}\mathbb I_{\Omega_{\mathcal F}}(x)$ with $\mathcal S\sim\omega$ is almost all zeros when $P_{\mathcal F}$ is small. The paper uses self-adaptive importance sampling to push the proposal toward the failure region step by step. Start from $h_1=\omega$; at step $k$, draw $N_1$ samples, sort them by limit-state value in **descending** order to get $\widetilde x^k_1,\dots,\widetilde x^k_{N_1}$, and set

$$
N_\eta=\max_{1\le i\le N_1}\{i:\ g(\widetilde x_i)>0\},
\qquad
N_p=\lfloor p_0N_1\rfloor .
$$

Here $N_\eta$ counts the samples that land in the failure set and $N_p$ is a fixed quantile. If $N_\eta<N_p$ (too few failure samples, meaning the proposal is not yet concentrated enough), update the truncated Gaussian proposal from the moments of the top $N_p$ samples:

$$
\mu_{k+1}=\frac{1}{N_p}\sum_{i=1}^{N_p}\widetilde x^k_i,
\qquad
\Sigma_{k+1}=\frac{1}{N_p-1}\sum_{i=1}^{N_p}
\bigl(\widetilde x^k_i-\mu_{k+1}\bigr)\otimes\bigl(\widetilde x^k_i-\mu_{k+1}\bigr),
$$

and set $h_{k+1}=\mathcal N_T(\mu_{k+1},\Sigma_{k+1})$, a Gaussian truncated to $\Omega$. This is a textbook cross-entropy-method step: use the most extreme handful of samples under the current proposal to define the next one. Conversely, once $N_\eta\ge N_p$ there are enough failure samples and the loop terminates.

On termination the final proposal switches to a **prior-weighted** mean:

$$
\mu_{opt}=\frac{\sum_{i=1}^{N_p}\widetilde x_i\,\omega(\widetilde x_i)}{\sum_{i=1}^{N_p}\omega(\widetilde x_i)},
\qquad
\Sigma_{opt}=\frac{1}{N_p-1}\sum_{i=1}^{N_p}
\bigl(\widetilde x_i-\mu_{opt}\bigr)\otimes\bigl(\widetilde x_i-\mu_{opt}\bigr).
$$

The importance-sampling estimate of the failure probability is

$$
\hat P^{SAIS}_{\mathcal F}=\frac{1}{N_2}\sum_{i=1}^{N_2}
\frac{\omega(x_i)}{\hat h_{opt}(x_i)}\,\mathbb I_{\Omega_{\mathcal F}}(x_i),
\qquad x_i\sim\hat h_{opt}.
$$

Those samples in the same batch that fall in $\Omega_{\mathcal F}$ serve directly as the new collocation points, so the indicator and the new points are **two outputs of one sampling pass**. The paper reports that $p_0=0.1$ makes this loop self-terminate quickly with good accuracy. Its Remark 1 adds that if $\Omega$ is unbounded the intermediate proposals are taken as plain (untruncated) Gaussians consistent with the form of the prior $\omega$; the paper itself treats only the truncated-Gaussian case.

> [!note] The intermediate and final updates are asymmetric
> The intermediate $\mu_{k+1}$ is an **unweighted** average while the final $\mu_{opt}$ is $\omega$-**weighted**. The intermediate steps only need to push the proposal toward the failure region; the final step needs it to approximate the zero-variance optimal density $\mathbb I_{\Omega_{\mathcal F}}\omega/P_{\mathcal F}$, which carries the prior factor, so weighting is the natural choice there.

### Theorems

**Assumption 4.1 (stability).** There exist constants $C_1,C_2>0$ independent of $v$ giving a two-sided bound between $\|v\|$ and $\|\mathcal Av\|_{2,\Omega}+\|\mathcal Bv\|_{2,\partial\Omega}$.

**Assumption 4.2.** The boundary residual satisfies $\|\mathcal B(u-u(\cdot;\theta^\ast))\|_{2,\partial\Omega}\le\epsilon_b$.

**Assumption 4.3.** The residual is bounded on the closed domain, $M:=\max_{x\in\Omega}|r(x;\theta^\ast)|<\infty$.

**Theorem 4.4.** Let $\Omega$ be bounded and let the FI-PINNs solution $u(x;\theta^\ast)$ satisfy Assumptions 4.1-4.3. Then

$$
\bigl\|u(x)-u(x;\theta^\ast)\bigr\|_{2,\Omega}
\le\sqrt2\,C_1^{-1}\Bigl(S_\Omega\bigl(M^2\epsilon_p+\epsilon_r^2\bigr)+\epsilon_b^2\Bigr)^{1/2},
$$

where $S_\Omega$ is the area of $\Omega$ and $\epsilon_r,\epsilon_p$ are the two **prescribed** tolerances.

The proof is short and worth remembering. Split the $L^2$ norm of the residual across the safe and failure sets,

$$
\|r\|^2_{2,\Omega}=\int_{\Omega_{\mathcal F}}r^2+\int_{\Omega_{\mathcal S}}r^2 .
$$

On the failure set the residual may be large but the area is small, and $S_{\Omega_{\mathcal F}}<S_\Omega\epsilon_p$ gives $\int_{\Omega_{\mathcal F}}r^2\le M^2S_\Omega\epsilon_p$; on the safe set the definition of the limit-state function gives $|r|<\epsilon_r$ and hence $\int_{\Omega_{\mathcal S}}r^2\le S_\Omega\epsilon_r^2$. So $\|r\|^2_{2,\Omega}\le S_\Omega(M^2\epsilon_p+\epsilon_r^2)$, and the lower bound in Assumption 4.1 converts the residual norm into the solution error norm while Assumption 4.2 handles the boundary term. Both prescribed tolerances therefore enter the error bound directly, which is the quantitative content of the claim that the failure probability is a posterior error indicator rather than an analogy for one.

> [!warning] A prior factor inside the proof
> The proof writes $P_{\mathcal F}=\int_\Omega\mathbb I_{\Omega_{\mathcal F}}(x)\,\mathrm dx=S_{\Omega_{\mathcal F}}/S_\Omega$, dropping the $\omega(x)$ present in the definition. The two agree only for the uniform prior $\omega\equiv1/S_\Omega$. For a non-uniform prior the step $S_{\Omega_{\mathcal F}}<S_\Omega\epsilon_p$ needs an extra condition.

### Numerical experiments

The outer loop (Algorithm 1) takes as input the network solution $u(x;\theta)$, the boundary set $\mathcal D_b$, the collocation set $\mathcal D_c$, a maximum round count $M$, the residual tolerance $\epsilon_r$ and the failure-probability tolerance $\epsilon_p$. Each round trains on the current sets, then runs self-adaptive importance sampling (Algorithm 2) to obtain $\hat P_{\mathcal F}$ and a new point set $\mathcal D_{adaptive}$ at once; if $\hat P_{\mathcal F}<\epsilon_p$ it stops, otherwise $\mathcal D_c\leftarrow\mathcal D_c\cup\mathcal D_{adaptive}$ and it continues. **The training set grows monotonically**, which is exactly what Part II changes.

The examples of Section 5 cover five classes of problem, all chosen so that uniform sampling is at a disadvantage:

| Example                                       | Why it is hard                                             |
| --------------------------------------------- | ---------------------------------------------------------- |
| 2D Poisson with a peaked or singular solution | the high-residual region is a tiny fraction of the domain  |
| Burgers' equation                             | the solution develops a steepening interior layer          |
| high-dimensional Poisson problem              | uniform candidates barely cover the active region          |
| Poisson on an unbounded 2D domain             | uniform sampling is not even definable on an unbounded set |
| time-dependent problem on an unbounded domain | both difficulties at once                                  |

The baselines are uniform sampling and residual-based adaptive refinement. The qualitative conclusion is that self-adaptive importance sampling concentrates points in the high-residual region and the error decays faster than under either baseline. The authors' group implementation is at [SEU-YL-UQ/FI-PINNs](https://github.com/SEU-YL-UQ/FI-PINNs).

### Relation to the others

This is the first part of a trilogy, continued by paper 73 (Part II: resampling plus subset simulation) and paper 76 (Part III: inverse problems, truncated Gaussian **mixtures**). It contrasts with density-driven adaptivity (papers 72, 80, 87): FI-PINNs adds points where the **residual** is large, whereas the ADDA/B-KRnet line samples from **the solution density itself**. The paper explicitly lists KRnet-based generative sampling as a third class of adaptive strategy alongside its own.

## 73: fixed-size training sets and subset simulation

### The idea

Part I leaves two visible gaps. First, the training set grows monotonically, so each round costs more than the last; by the later rounds most of the points added early no longer carry information yet still participate in every gradient computation. Second, the failure probability is estimated with a **single** truncated Gaussian — a unimodal, ellipsoidal proposal. When the failure region is multimodal, or simply oddly shaped (distributed along a curve, say), a unimodal proposal either misses some modes or becomes too broad to cover them all, and in either case the variance of the importance weights explodes.

Part II addresses these separately: replace "accumulate" by "fixed-size resampling", and replace "a single truncated Gaussian" by "subset simulation".

### Derivation

**Extension one: fixed-size resampling under cosine annealing.** The training set keeps a constant size while the composition of the collocation points shifts gradually from uniform to adaptive under a cosine annealing schedule. The motivation is that early in training the network is still poor, its residual field is not trustworthy, and adding points according to it merely amplifies noise; only as training proceeds does the residual field become a meaningful indicator. The annealing therefore encodes "when may the adaptive signal be trusted" as a rule that varies with training progress.

**Extension two: subset simulation as the posterior model.** Subset simulation writes a small probability as a product of larger conditional probabilities over a sequence of **nested intermediate failure levels**: take thresholds $\epsilon_r^{(1)}>\epsilon_r^{(2)}>\cdots$, so that the conditional probability at each level is not too small and can be estimated with a moderate sample size, with the samples inside each level generated by MCMC started from the samples of the previous level. This yields both an estimate of $P_{\mathcal F}$ and samples in the failure region, and it presupposes nothing about the shape of that region — which is precisely its advantage over a single truncated Gaussian.

### Numerical experiments

The paper reports a significant improvement over Part I on several challenging problems: the subset-simulation posterior model estimates the failure probability more efficiently and generates effective new training points in the failure region.

### Relation to the others

The middle part of the trilogy. Relative to Part I what changes is the **posterior model** (truncated Gaussian to subset simulation) and the **set-management strategy** (accumulation to resampling); Part III changes the posterior model once more (to a truncated Gaussian mixture) and turns to inverse problems. The cosine-annealed transition from uniform to adaptive is conceptually the same device as the "start uniform, then sample from the current model" schedule in the ADDA line (papers 64, 72, 80, 87).

## 76: a truncated Gaussian mixture, aimed at inverse problems

### The idea

Inverse problems amplify the difficulty of Part I by a further step. Two networks are now trained jointly — the state $u(x;\theta_u)$ and the coefficient $\gamma(x;\theta_\gamma)$ — and beyond the equation and boundary residuals the loss carries a data-misfit term. The two networks constrain each other, so the residual field splits more readily into **several disconnected high-residual regions**: a coefficient misestimated in one place drives the state off in another. A unimodal proposal is not enough for that geometry.

Part III's choice is to replace the proposal by a truncated Gaussian **mixture** fitted by EM. A mixture handles multimodality naturally, and EM is the standard tool for fitting one at controllable cost.

### Setting

The inverse problem reads

$$
\mathcal A[x;u(x),\gamma(x)]=0\ \text{in}\ \Omega,
\qquad
\mathcal B[x;u(x),\gamma(x)]=0\ \text{on}\ \partial\Omega,
$$

with indirect observations $y(x)=\mathcal G[x;u(x),\gamma(x)]$ for $x\in D_{\mathrm{indirect}}$. The two networks are trained on

$$
\mathcal L(\theta)=\mathcal L_c(\theta)+\lambda\mathcal L_b(\theta)+\mu\mathcal L_d(\theta),
\qquad
\mathcal L_c=\|r(x;\theta)\|_R^2,\quad
\mathcal L_b=\|b(x;\theta)\|_B^2,\quad
\mathcal L_d=|d(\theta)|_D^2
$$

with $d(x_i;\theta)=y(x_i)-\mathcal G[x_i;u(x_i,\theta_u),\gamma(x_i,\theta_\gamma)]$. The limit-state function becomes $g(x)=|r(x;\theta)|_{\tilde R}-\varepsilon_r$, where $\tilde R$ sums the absolute values of the terms entering the $R$ norm; for the impedance-tomography example below it reduces to $|r(x;\theta)|$.

### Derivation

The paper writes out the zero-variance optimal proposal explicitly,

$$
h_{\mathrm{opt}}(x)=\frac{\mathbb I_{\Omega_{\mathcal F}}(x)\,\omega(x)}{P_{\mathcal F}}
=\frac{\mathbb I_{\{g(x)>0\}}(x)\,\omega(x)}{\int_\Omega\mathbb I_{\{g(x)>0\}}(x)\omega(x)\,\mathrm dx},
$$

and approximates it by a truncated Gaussian mixture,

$$
h_{\mathrm{opt}}\approx h(x;\eta)=\sum_{k=1}^{K}\pi_k\,\mathcal N(x;\mu_k,\Sigma_k),
\qquad \pi_k\ge0,\ \sum_k\pi_k=1 .
$$

The fitting criterion is $\min_\eta D_{\mathrm{KL}}(h_{\mathrm{opt}}\|h(\cdot;\eta))$. Expanding the KL divergence, only the cross-entropy term depends on $\eta$, so after discretisation it becomes a maximum-likelihood problem over $x_j\sim h_{\mathrm{opt}}$,

$$
\max_\eta\ \log\prod_{j=1}^{N_c}h(x_j,\eta),
$$

and drawing $x_j\sim h_{\mathrm{opt}}$ only requires keeping those samples of the current proposal that land in the failure set. The parameters come from EM iterations whose E step computes the responsibilities

$$
q^{(t)}_{k,j}=\frac{\pi^{(t)}_k\,\mathcal N\bigl(x_j;\mu^{(t)}_k,\Sigma^{(t)}_k\bigr)}
{\sum_{k=1}^{K}\pi^{(t)}_k\,\mathcal N\bigl(x_j;\mu^{(t)}_k,\Sigma^{(t)}_k\bigr)},
$$

and whose M step updates the three parameter groups

$$
\pi^{(t+1)}_k=\frac1{N_c}\sum_{j}q^{(t)}_{k,j},
\qquad
\mu^{(t+1)}_k=\frac{\sum_j q^{(t)}_{k,j}x_j}{\sum_j q^{(t)}_{k,j}},
\qquad
\Sigma^{(t+1)}_k=\frac{\sum_j q^{(t)}_{k,j}\bigl(x_j-\mu^{(t+1)}_k\bigr)\bigl(x_j-\mu^{(t+1)}_k\bigr)^{\!\top}}{\sum_j q^{(t)}_{k,j}} .
$$

The implementation calls `sklearn.mixture` directly (Remark 3.1 of the paper).

Truncation is realised by projection rather than by a truncated density formula: $\mathrm{Proj}_\Omega(x)=\arg\min_{y\in\bar\Omega}\|x-y\|_2$. The EM update at round $l$ is written $\eta^{\star,l+1}=\mathrm{EM}(\mathrm{Proj}_\Omega(x^l),K^l)$, and the final fit uses the $N_p$ failure points with the largest limit-state values. The importance-sampling estimate becomes

$$
\hat P^{\mathrm{SAIS}}_{\mathcal F}=\frac1{N_2}\sum_{i=1}^{N_2}
\frac{\omega\bigl(\mathrm{Proj}_\Omega(x_i)\bigr)}{\hat h_{\mathrm{opt}}\bigl(\mathrm{Proj}_\Omega(x_i)\bigr)}
\,\mathbb I_{\Omega_{\mathcal F}}\bigl(\mathrm{Proj}_\Omega(x_i)\bigr).
$$

### Theorems

The paper proves no new convergence theorem. Its contribution is the proposal model plus the fitting route: compared with the alternative of fitting a mixture by maximising a risk function, which is forced to diagonalise $\Sigma_k$ in high dimension, EM permits full covariances. Remark 3.3 supplies the cost judgement: the time EM adds is negligible next to network training.

### Numerical experiments

Two classes of inverse problem. The first is the **inverse conductivity problem of electrical impedance tomography**:

$$
-\nabla\!\cdot\!(\gamma\nabla u)-f=0\ \text{in}\ \Omega,
\qquad
u=u_b,\ \ \gamma=\gamma_b,\ \ \partial_{\vec n}u=u_n\ \text{on}\ \partial\Omega,
$$

with loss $\mathcal L_c+\lambda_{u_b}\mathcal L_{u_b}+\lambda_{\gamma_b}\mathcal L_{\gamma_b}+\lambda_{\partial u}\mathcal L_{\partial u}$, where $\mathcal L_c=\|-\nabla\cdot(\gamma\nabla u)-f\|^2_{L^2(\Omega)}$ and the three boundary terms are $L^2(\partial\Omega)$ norms. The weighted norm is defined by $\|f\|^2_{L^2(\Omega)}=\int_\Omega f^2(x)\omega(x)\,\mathrm dx$, and this example takes $\omega\equiv1$ — worth noting, because the caveat about the prior factor above is not an issue when $\omega\equiv1$. The second class is an **inverse source problem in a parabolic system**.

### Relation to the others

The third part of the trilogy, and the only one aimed at inverse problems with two coupled networks. The differences between the three parts sit almost entirely in the **posterior model** that estimates $P_{\mathcal F}$ and generates new points: a single truncated Gaussian (70), then subset simulation (73), then an EM-fitted truncated Gaussian mixture (76). The variational principle "approximate the optimal density by minimising a KL divergence" is the same one papers 80 and 87 use when they replace the mixture by a normalising flow.

## 80: a variational loss has no residual, so use the variance instead

### The idea

The deep Ritz method minimises a variational energy rather than a pointwise residual, so none of the machinery above transfers: with no $r(x;\theta)$ there is no limit-state function and no failure set.

The paper changes the angle of attack. A deep Ritz loss is an integral, and that integral is approximated by Monte Carlo quadrature, so the discretisation error naturally splits into two pieces: whether the network can represent the true solution (the approximation error), and whether finitely many samples can evaluate the integral accurately (the **statistical error**). The second dominates when the integrand has low regularity. The paper turns this into a number with an example one can do by hand: take $G(u(x))=\frac{1}{\sqrt{2\pi}}e^{-x^2/2\sigma^2}$ on $[-1,1]$, for which the relative Monte Carlo error is $C(\sigma N)^{-1/2}$, so $O(1/\sigma)$ uniform samples are needed for $O(1)$ relative accuracy. When $\sigma$ is small — that is, when the integrand concentrates in a narrow peak — uniform sampling is catastrophic.

The target of adaptivity should therefore be **variance reduction in the variational loss**, and importance sampling is precisely the tool designed for that: treat the integrand as an unnormalised density, learn it, then sample from it.

### Setting

The penalised deep Ritz problem is

$$
\min_{u\in V} J(u),
\qquad
J(u)=\int_\Omega G\bigl(u(x)\bigr)\,\mathrm dx
+\beta\,\bigl\|B\bigl(x,u(x)\bigr)\bigr\|^2_{\partial\Omega,2},
$$

discretised on uniform collocation sets $S_\Omega=\{x_i\}_{i=1}^{N_v}$ and $S_{\partial\Omega}=\{x_j\}_{j=1}^{N_b}$ as

$$
J_N\bigl(u(x,\theta)\bigr)=\frac1{N_v}\sum_{i=1}^{N_v}G\bigl(u(x_i,\theta)\bigr)
+\frac{\beta}{N_b}\sum_{j=1}^{N_b}B\bigl(x_j,u(x_j,\theta)\bigr)^2 .
$$

The importance-sampling estimator is

$$
I(u)=\int_\Omega G\bigl(u(x,\theta)\bigr)\mathrm dx
=\mathbb E_p\Bigl[\frac{G(u(X,\theta))}{p(X)}\Bigr]
\approx\frac1{N_v}\sum_{i=1}^{N_v}\frac{G(u(x_i,\theta))}{p(x_i)},
\qquad x_i\sim p .
$$

### Derivation

**A sign-changing integrand: adaptivity cannot reach zero variance.** If $G\ge0$, the optimal proposal $p^\star=G/\mu$ with $\mu=\int_\Omega G$ makes the integrand $G/p$ identically the constant $\mu$ and the variance zero. But a variational integrand need not be non-negative — the Dirichlet energy minus the source term already changes sign. In that case the optimal choice is

$$
p^\star(x)=\frac{|G(u(x,\theta))|}{\mu},
\qquad \mu=\int_\Omega|G(u(x,\theta))|\,\mathrm dx,
$$

which by Cauchy-Schwarz minimises the variance among all densities but leaves a strictly positive residual:

$$
\sigma_{p^\star}=\Bigl(\int_\Omega|G(u(x,\theta))|\mathrm dx\Bigr)^2
-\Bigl(\int_\Omega G(u(x,\theta))\mathrm dx\Bigr)^2>0 .
$$

The right-hand side is the squared integral of $|G|$ minus the squared integral of $G$, which vanishes only when $G$ does not change sign. This is a structural conclusion: **when the integrand changes sign, adaptivity alone cannot remove the statistical error and the sample count** $N_v$ **must also grow**, because the Monte Carlo error is $\sigma_p/\sqrt{N_v}\ge\sigma_{p^\star}/\sqrt{N_v}>0$. It draws a clear boundary around what adaptive sampling can deliver.

**Two networks updated alternately.** The density model is the bounded KRnet (see the [[en/computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|density-flow page]]):

$$
p_{\text{bKRnet}}(x,\theta_f)=p_Z\bigl(f_{\text{bKRnet}}(x,\theta_f)\bigr)\,
\bigl|\det\nabla_x f_{\text{bKRnet}}\bigr|,
\qquad Z\sim\mathrm{Unif}([-1,1]^d),
$$

with sampling by the inverse flow $X=f^{-1}_{\text{bKRnet}}(Z)$. In the training criterion $\min_{\theta_f}D_{\mathrm{KL}}(p^\star\|p_{\text{bKRnet}})$ only the cross-entropy term depends on $\theta_f$, so it reduces to minimising $H(p^\star,p_{\text{bKRnet}})$; and since $p^\star$ can only be evaluated in the unnormalised form $|G|$, the cross entropy itself must be estimated by importance sampling:

$$
H(p^\star,p_{\text{bKRnet}})\approx
-\frac1{N_v}\sum_{i=1}^{N_v}
\frac{|G(u(x_i,\theta))|\,\log p_{\text{bKRnet}}(x_i,\theta_f)}{\mu\,\tilde p(x_i)},
\qquad x_i\sim\tilde p .
$$

The two networks alternate at adaptivity round $k$, the solution network performing an importance-weighted Ritz minimisation under the current density and the density network a cross-entropy minimisation under the new solution:

$$
\theta^{k+1}=\arg\min_\theta\ \frac1{N_v}\sum_{i}
\frac{G\bigl(u(x^k_{\Omega,i},\theta)\bigr)}{p_{\text{bKRnet}}(x^k_{\Omega,i},\theta_f^k)}
+\frac{\beta}{N_b}\sum_{j}B^2\bigl(x^k_{\partial\Omega,j},u(x^k_{\partial\Omega,j},\theta)\bigr),
$$

$$
\theta_f^{k+1}=\arg\min_{\theta_f}\
-\frac1{N_B}\sum_{l}
\frac{\bigl|G\bigl(u(x^k_{B,l},\theta^{k+1})\bigr)\bigr|\,
\log p_{\text{bKRnet}}(x^k_{B,l},\theta_f)}
{\mu^{k+1}\,p_{\text{bKRnet}}(x^k_{B,l},\theta_f^k)} .
$$

The initial density $p_{\text{bKRnet}}(\cdot,\theta_f^0)$ is uniform, and $\mu^{k+1}$ can be dropped as a constant.

**Two mixture models give the density a floor.** The ratio $G/p$ blows up wherever the learned density is small, so the paper puts a floor under it. Model one mixes with the uniform density,

$$
p_{\mathrm{mixture}}(x,\theta_f^{k+1})
=\epsilon\,p_{\text{bKRnet}}(x,\theta_f^{k+1})+(1-\epsilon)\,p_{\mathrm{uniform}}(x)
\ \ge\ \frac{1-\epsilon}{|\Omega|} .
$$

Model two mixes recursively with **every** previous flow,

$$
p_{\mathrm{mixture}}\bigl(x,\{\theta_f^t\}_{t\le k+1}\bigr)
=\epsilon\,p_{\text{bKRnet}}(x,\theta_f^{k+1})
+(1-\epsilon)\,p_{\mathrm{mixture}}\bigl(x,\{\theta_f^t\}_{t\le k}\bigr)
=\sum_{t=1}^{k+1}\epsilon(1-\epsilon)^{k+1-t}p_{\text{bKRnet}}(x,\theta_f^t)
+\frac{(1-\epsilon)^{k+1}}{|\Omega|},
$$

with floor $(1-\epsilon)^{k+1}/|\Omega|$, which shrinks exponentially with the round count — the price of model two, paid in exchange for a proposal that covers every region ever judged important.

### Theorems

The paper proves no convergence theorem for the coupled adaptive iteration. Its theoretical content is the variance ordering: by Cauchy-Schwarz, $\sigma_{p^\star}\le\sigma_p$ for every density $p$; combined with the explicit expression showing $\sigma_{p^\star}>0$ in the sign-changing case, this gives a lower bound on the Monte Carlo error,

$$
\mathrm{Var}_p^{1/2}[\cdot]=\frac{\sigma_p}{\sqrt{N_v}}\ \ge\ \frac{\sigma_{p^\star}}{\sqrt{N_v}}>0 .
$$

### Numerical experiments

Algorithm 1 (the basic version) is structured as follows: input an initial solution network $u(x,\theta^0)$, an initial density $p_{\text{bKRnet}}(x,\theta_f^0)$, a maximum round count $N_e$, a batch size $m$, and three initial training sets $S^0_\Omega$, $S^0_{\partial\Omega}$, $S^0_B$; the outer loop runs over $k=0,\dots,N_{\mathrm{adaptive}}-1$, and the inner loop first solves the PDE by stochastic gradient descent, then trains the bounded KRnet by stochastic gradient descent, then draws new sets $S^{k+1}_\Omega$ and $S^{k+1}_B$ from $p_{\text{bKRnet}}(x,\theta_f^{k+1})$. Algorithm 2 makes the new sets the union of uniform points and flow samples and evaluates densities with mixture model one; Algorithm 3 stores $\theta_f^t$ from **every** round, resamples from all the flows, and evaluates densities with mixture model two.

Four examples are compared side by side under four sampling strategies:

| Example     | Problem                            |
| ----------- | ---------------------------------- |
| Section 4.1 | a two-dimensional peak problem     |
| Section 4.2 | a two-peak problem                 |
| Section 4.3 | a problem with a singularity       |
| Section 4.4 | a high-dimensional Poisson problem |

The four strategies are conventional deep Ritz (fresh uniform points each round), Algorithm 1, Algorithm 2 and Algorithm 3. The shared setup is: error measured as a relative discrete $L_2$ error on a tensor grid, replaced by uniform samples in high dimension; a ResNet-like architecture with $\sin^3(x)$ activations; Adam for both networks; and a bounded KRnet with the same architecture as in paper 87. The paper's conclusion is that the adaptive methods improve accuracy over plain deep Ritz, **particularly on low-regularity and high-dimensional problems**.

### Relation to the others

Paper 80 is the deep Ritz counterpart of the FI-PINNs line (papers 70, 73, 76) — the same "learn a density, resample, retrain" loop, but with target density $|G(u)|/\mu$, the variational integrand, rather than a residual-based failure indicator. Its density model is exactly the bounded KRnet of paper 87, and the conclusion of paper 87 explicitly lists this coupling as future work. Mixing with the uniform density (model one) is the same kind of safety valve as the update rate $\gamma$ in paper 87 and the cosine-annealed uniform-to-adaptive transition in paper 73. It also complements paper 60, which repairs **how boundary conditions are imposed** in the same deep Ritz objective.

## The shared structure

| No. | Target being driven                    | Proposal or density model         | Training-set update                          |
| --- | -------------------------------------- | --------------------------------- | -------------------------------------------- |
| 66  | unbiasedness of the operator estimate  | two independent random draws      | not applicable (sampling is in the operator) |
| 70  | residual failure probability           | a single truncated Gaussian       | monotone accumulation                        |
| 73  | residual failure probability           | subset simulation                 | fixed-size resampling with annealing         |
| 76  | residual failure probability plus data | truncated Gaussian mixture via EM | monotone accumulation                        |
| 80  | variance of the variational loss       | bounded KRnet mixed with uniform  | resample from the current density            |

Three transferable judgements:

- **The indicator should enter the error bound.** Theorem 4.4 of paper 70 puts both prescribed tolerances into the final estimate, which is the quantitative version of "the failure probability is a posterior error indicator" rather than an analogy.
- **The expressiveness of the proposal sets the ceiling.** From a single truncated Gaussian to subset simulation to a mixture model, what improves is always the estimation of the same quantity.
- **Adaptivity has a ceiling of its own.** The residual variance in paper 80 shows that with a sign-changing integrand the sample count remains a necessary independent resource.

## Sources for this page

- L. Guo, H. Wu, X. Yu, and T. Zhou, _Monte Carlo fPINNs: deep learning method for forward and inverse problems involving high dimensional fractional partial differential equations_, Comput. Methods Appl. Mech. Engrg. 400 (2022), 115523 (preprint [arXiv:2203.08501](https://arxiv.org/abs/2203.08501)).
- Z. Gao, L. Yan, and T. Zhou, [_Failure-informed adaptive sampling for PINNs_](https://doi.org/10.1137/22M1527763), SIAM J. Sci. Comput. 45(4) (2023), pp. A1971-A1994 (preprint [arXiv:2210.00279](https://arxiv.org/abs/2210.00279)).
- Z. Gao, T. Tang, L. Yan, and T. Zhou, [_Failure-informed adaptive sampling for PINNs, Part II: combining with re-sampling and subset simulation_](https://doi.org/10.1007/s42967-023-00312-7), Commun. Appl. Math. Comput. 6 (2024), pp. 1720-1741 (preprint [arXiv:2302.01529](https://arxiv.org/abs/2302.01529)).
- W. Liu, L. Yan, T. Zhou, and Y. Zhou, [_Failure-informed adaptive sampling for PINNs, Part III: applications to inverse problems_](https://doi.org/10.4208/csiam-am.SO-2023-0059), CSIAM Trans. Appl. Math. 5(3) (2024), pp. 636-670.
- X. Wan, T. Zhou, and Y. Zhou, _Adaptive importance sampling for deep Ritz_, Commun. Appl. Math. Comput. 7(3) (2025), pp. 929-953 (preprint [arXiv:2310.17185](https://arxiv.org/abs/2310.17185)).
