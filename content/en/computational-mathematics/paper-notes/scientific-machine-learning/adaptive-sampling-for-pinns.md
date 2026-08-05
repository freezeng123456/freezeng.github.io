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

> [!note] Coverage of this page
> Papers **66** (_Comput. Methods Appl. Mech. Engrg._ 400, 2022), **70** (_SIAM J. Sci. Comput._ 45(4), 2023), **73** (_Commun. Appl. Math. Comput._ 6, 2024), **76** (_CSIAM Trans. Appl. Math._ 5(3), 2024) and **80** (_Commun. Appl. Math. Comput._ 7(3), 2025).

![Let reliability analysis choose the next collocation batch](assets/diagrams/tao-zhou-papers/en/failure-informed-sampling.svg)

The accuracy of a physics-informed neural network depends strongly on where the collocation points sit. The shared judgement across these papers is that **point placement is not an implementation detail but an object that can be driven by an explicit mathematical target**. The four targets differ — the probability that the residual is too large, the variance of a variational loss, the unbiasedness of a conditional expectation — but each reduces to the same recipe: define an estimable scalar, then generate points from it.

## 66: write the operator itself as an expectation

Before collocation points are even discussed, paper 66 removes a different obstruction: assembling a fractional operator is too expensive. The earlier fPINN approach discretises the fractional derivative on an auxiliary grid, which does not scale beyond very low dimension. This paper's choice is to **never discretise the fractional operator at all** and instead write it as an expectation estimated by Monte Carlo inside the loss.

The model equation is a fractional advection-diffusion equation,

$$
L[u(x,t)]:=\frac{\partial^\gamma u(x,t)}{\partial t^\gamma}
+c\,(-\Delta)^{\alpha/2}u(x,t)+v\cdot\nabla u(x,t)=f(x,t),
$$

with the integral fractional Laplacian

$$
(-\Delta)^{\alpha/2}u(x)\triangleq C_{d,\alpha}\,\mathrm{P.V.}\!
\int_{\mathbb R^d}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\,\mathrm dy .
$$

### Inner-outer splitting and two Beta distributions

Split the integral at radius $r_0$ around $x$. Inside the ball, with direction $\xi\sim\mathrm{Unif}(S^{d-1})$ and radius $r\sim f_I$,

$$
\int_{y\in B_{r_0}(x)}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\mathrm dy
=\frac{|S^{d-1}|\,r_0^{2-\alpha}}{2(2-\alpha)}\;
\mathbb E_{\xi,\,r\sim f_I}\!\left[\frac{2u(x)-u(x-r\xi)-u(x+r\xi)}{r^2}\right],
$$

$$
f_I(r)=\frac{2-\alpha}{r_0^{2-\alpha}}\,r^{1-\alpha}\mathbf 1_{r\in[0,r_0]},
\qquad r/r_0\sim\mathrm{Beta}(2-\alpha,1).
$$

Outside the ball, with $r\sim f_O$,

$$
\int_{y\notin B_{r_0}(x)}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\mathrm dy
=\frac{|S^{d-1}|\,r_0^{-\alpha}}{2\alpha}\;
\mathbb E_{\xi,\,r\sim f_O}\bigl[2u(x)-u(x-r\xi)-u(x+r\xi)\bigr],
$$

$$
f_O(r)=\alpha\,r_0^{\alpha}\,r^{-1-\alpha}\mathbf 1_{r\in[r_0,\infty)},
\qquad r_0/r\sim\mathrm{Beta}(\alpha,1).
$$

What these identities accomplish is that the singularity of the integrand is absorbed into the sampling density, leaving only a second difference inside the expectation. Because

$$
\lim_{r\to0}\frac{2u(x)-u(x-r\xi)-u(x+r\xi)}{r^2}
=\partial^2_r u(x+r\xi)\big|_{r=0},
$$

a floor $r_\epsilon=\max\{\epsilon,r_I\}$ is still needed numerically so that $1/r^2$ does not amplify round-off.

The Caputo time derivative is handled the same way:

$$
\frac{\partial^\gamma u(x,t)}{\partial t^\gamma}
=\frac{\gamma}{1-\gamma}\,t^{1-\gamma}\,
\mathbb E_{\tau\sim f_{I,t}}\!\left[\frac{u(x,t)-u(x,t-\tau)}{\tau}\right]
+\frac{u(x,t)-u(x,0)}{t^{\gamma}},
\qquad
f_{I,t}(\tau)=(1-\gamma)\tau^{-\gamma}\mathbf 1_{\tau\in[0,1]} .
$$

### Unbiasedness through two independent draws

This is the ingredient the abstract advertises. Squaring a noisy operator estimate inflates the residual loss **systematically** by the estimator variance, since $\mathbb E[\widehat L^2]=(\mathbb E\widehat L)^2+\mathrm{Var}(\widehat L)$. The fix is to draw **two independent groups** of random parameters per collocation point and multiply the two independent estimates instead of squaring one:

$$
\hat L_{equ}(\theta)=\frac{1}{mN_u}\sum_{i,j}
\widehat L\bigl[u_{NN}(x_i,t_i;\theta);\tau_j,\xi_j,r_{Ij},r_{oj}\bigr]\cdot
\widehat L\bigl[u_{NN}(x_i,t_i;\theta);\tau'_j,\xi'_j,r'_{Ij},r'_{oj}\bigr].
$$

Provided $\epsilon=\epsilon_t=0$ and round-off is ignored, $\mathbb E[\hat L_{equ}(\theta)]=L_{equ}(\theta)$.

The experiments fall into three groups: high-dimensional integral fractional Laplacian equations, parameter identification (inverse problems) for time-space fractional PDEs, and fractional diffusion equations with random inputs.

> [!note] Title and source reconciliation
> The preprint is titled "Monte Carlo PINNs: deep learning approach for forward and inverse problems involving high dimensional fractional partial differential equations" while the journal version reads "Monte Carlo fPINNs: Deep learning method for ...". Separately, the prefactor in the Caputo definition renders as $\Gamma(1-\alpha)$ in the preprint while the exponent uses $\gamma$; with $0<\gamma<1$ it should be $\Gamma(1-\gamma)$. Check the journal version before quoting a prefactor.

## 70: define "the residual is too large" as a failure event

### Why the largest-residual points are not enough

Standard residual-based adaptive refinement takes the $m$ points of largest residual from a uniform candidate pool. When the high-residual region occupies a small fraction of the domain, almost no uniform candidate lands in it and the strategy fails. Paper 70 changes the frame and treats point placement as a **reliability analysis** problem.

The PINN loss is

$$
\mathcal L(\theta)=\mathcal L_c(\theta)+\lambda\mathcal L_b(\theta),
\qquad
\mathcal L_c(\theta)=\frac{1}{N_c}\sum_{i=1}^{N_c}\bigl|r(x^c_i;\theta)\bigr|^2,
\qquad
\mathcal L_b(\theta)=\frac{1}{N_b}\sum_{i=1}^{N_b}\bigl|b(x^b_i;\theta)\bigr|^2 .
$$

Define the limit-state function

$$
g(x)=\bigl|r(x;\theta)\bigr|-\epsilon_r,
$$

whose zero level set splits the domain into a safe set $\Omega_{\mathcal S}=\{g<0\}$ and a failure set $\Omega_{\mathcal F}=\{g>0\}$. Under a prior $\omega(x)$, the failure probability

$$
P_{\mathcal F}=\int_\Omega \omega(x)\,\mathbb I_{\Omega_{\mathcal F}}(x)\,\mathrm dx
$$

plays the role of a **posterior error indicator**, exactly as in adaptive finite elements: the network is declared reliable once $P_{\mathcal F}<\epsilon_p$.

One concrete difference from residual-based refinement is worth naming: Monte Carlo enrichment adds **every** candidate point that falls in $\Omega_{\mathcal F}$, so the number added varies per round, whereas residual refinement adds a fixed $m$.

### Self-adaptive importance sampling

When the failure probability is small, a plain Monte Carlo estimate is useless. The paper pushes the proposal density toward the failure region step by step. Start from $h_1=\omega$; at step $k$, draw $N_1$ samples, sort them by limit-state value in **descending** order to get $\widetilde x^k_1,\dots,\widetilde x^k_{N_1}$, and set

$$
N_\eta=\max_{1\le i\le N_1}\{i:\ g(\widetilde x_i)>0\},
\qquad
N_p=\lfloor p_0N_1\rfloor .
$$

If $N_\eta<N_p$ (too few failure samples), refine the truncated Gaussian proposal from the top $N_p$ samples:

$$
\mu_{k+1}=\frac{1}{N_p}\sum_{i=1}^{N_p}\widetilde x^k_i,
\qquad
\Sigma_{k+1}=\frac{1}{N_p-1}\sum_{i=1}^{N_p}
\bigl(\widetilde x^k_i-\mu_{k+1}\bigr)\otimes\bigl(\widetilde x^k_i-\mu_{k+1}\bigr),
$$

and set $h_{k+1}=\mathcal N_T(\mu_{k+1},\Sigma_{k+1})$. On termination the final proposal uses a **prior-weighted** mean:

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

The paper reports that $p_0=0.1$ makes this loop self-terminate quickly with good accuracy.

> [!note] The intermediate and final updates are asymmetric
> The intermediate $\mu_{k+1}$ is an **unweighted** average while the final $\mu_{opt}$ is $\omega$-**weighted**. The intermediate steps only need to push the proposal toward the failure region; the final step needs it to approximate the zero-variance optimal density $\mathbb I_{\Omega_{\mathcal F}}\omega/P_{\mathcal F}$, which carries the prior factor, so weighting is the natural choice there.

### Theorem 4.4: both tolerances enter the error bound

Under three assumptions — (4.1) constants $C_1,C_2>0$ independent of $v$ giving a two-sided bound between $\|v\|$ and $\|\mathcal Av\|_{2,\Omega}+\|\mathcal Bv\|_{2,\partial\Omega}$; (4.2) boundary residual $\|\mathcal B(u-u(\cdot;\theta^\ast))\|_{2,\partial\Omega}\le\epsilon_b$; (4.3) $M:=\max_{x\in\Omega}|r(x;\theta^\ast)|<\infty$ — a bounded domain satisfies

$$
\bigl\|u(x)-u(x;\theta^\ast)\bigr\|_{2,\Omega}
\le\sqrt2\,C_1^{-1}\Bigl(S_\Omega\bigl(M^2\epsilon_p+\epsilon_r^2\bigr)+\epsilon_b^2\Bigr)^{1/2},
$$

with $S_\Omega$ the area of $\Omega$.

The proof is short and worth remembering. Split the residual norm across the two sets,

$$
\|r\|^2_{2,\Omega}=\int_{\Omega_{\mathcal F}}r^2+\int_{\Omega_{\mathcal S}}r^2 .
$$

On the failure set the residual may be large but the area is small, and $S_{\Omega_{\mathcal F}}<S_\Omega\epsilon_p$ gives $\int_{\Omega_{\mathcal F}}r^2\le M^2S_\Omega\epsilon_p$; on the safe set $|r|<\epsilon_r$ gives $\int_{\Omega_{\mathcal S}}r^2\le S_\Omega\epsilon_r^2$. Both prescribed tolerances $\epsilon_p$ and $\epsilon_r$ therefore appear in the final estimate, which is the quantitative content of the claim that the failure probability is a posterior error indicator rather than an analogy for one.

> [!warning] A prior factor inside the proof
> The proof writes $P_{\mathcal F}=\int_\Omega\mathbb I_{\Omega_{\mathcal F}}(x)\,\mathrm dx=S_{\Omega_{\mathcal F}}/S_\Omega$, dropping the $\omega(x)$ present in the definition. The two agree only for the uniform prior $\omega\equiv1/S_\Omega$. For a non-uniform prior the step $S_{\Omega_{\mathcal F}}<S_\Omega\epsilon_p$ needs an extra condition.

The experiments cover a two-dimensional Poisson equation with a peaked or singular solution, Burgers' equation, a high-dimensional Poisson problem, a Poisson problem on an unbounded two-dimensional domain, and a time-dependent unbounded problem, against uniform sampling and residual-based refinement. The authors' group repository is [SEU-YL-UQ/FI-PINNs](https://github.com/SEU-YL-UQ/FI-PINNs).

## 73 and 76: replacing the posterior model

Across the trilogy, what changes is almost entirely **which model estimates the failure probability and generates new points**.

### 73: fixed-size resampling and subset simulation

Paper 70 grows the training set monotonically ($\mathcal D_c\leftarrow\mathcal D_c\cup\mathcal D_{adaptive}$), so training cost rises with each round, and a single truncated Gaussian is a crude proposal when the failure region is multimodal or oddly shaped. Part II addresses both: the training set stays a constant size while its composition shifts from uniform to adaptive under a cosine annealing schedule, and the failure probability is estimated by subset simulation, which writes a small probability as a product of larger conditional probabilities over a nested sequence of intermediate failure levels with MCMC inside each level.

> [!note] What could be verified
> The exact annealing formula, the construction of the intermediate levels and the MCMC kernel could not be confirmed equation by equation from the public material used here, so this page reports only their role and position.

### 76: a truncated Gaussian mixture, aimed at inverse problems

Part III targets **inverse problems**, where a state network and a coefficient network are trained jointly and the residual can concentrate in several disconnected regions, making a unimodal proposal insufficient. The loss is

$$
\mathcal L(\theta)=\mathcal L_c(\theta)+\lambda\mathcal L_b(\theta)+\mu\mathcal L_d(\theta),
$$

with an equation residual, a boundary residual and a data misfit $d(x_i;\theta)=y(x_i)-\mathcal G[x_i;u(x_i,\theta_u),\gamma(x_i,\theta_\gamma)]$.

The paper writes out the zero-variance optimal proposal explicitly,

$$
h_{\mathrm{opt}}(x)=\frac{\mathbb I_{\Omega_{\mathcal F}}(x)\,\omega(x)}{P_{\mathcal F}}
=\frac{\mathbb I_{\{g(x)>0\}}(x)\,\omega(x)}{\int_\Omega\mathbb I_{\{g(x)>0\}}(x)\omega(x)\,\mathrm dx},
$$

and approximates it by a truncated Gaussian mixture,

$$
h_{\mathrm{opt}}\approx h(x;\eta)=\sum_{k=1}^{K}\pi_k\,\mathcal N(x;\mu_k,\Sigma_k),
\qquad \pi_k\ge0,\ \sum_k\pi_k=1,
$$

fitted by $\min_\eta D_{\mathrm{KL}}(h_{\mathrm{opt}}\|h(\cdot;\eta))$, which discretises to the maximum-likelihood problem $\max_\eta\log\prod_{j=1}^{N_c}h(x_j,\eta)$ over $x_j\sim h_{\mathrm{opt}}$. The parameters come from EM iterations:

$$
q^{(t)}_{k,j}=\frac{\pi^{(t)}_k\,\mathcal N\bigl(x_j;\mu^{(t)}_k,\Sigma^{(t)}_k\bigr)}
{\sum_{k=1}^{K}\pi^{(t)}_k\,\mathcal N\bigl(x_j;\mu^{(t)}_k,\Sigma^{(t)}_k\bigr)},
$$

$$
\pi^{(t+1)}_k=\frac1{N_c}\sum_{j}q^{(t)}_{k,j},
\qquad
\mu^{(t+1)}_k=\frac{\sum_j q^{(t)}_{k,j}x_j}{\sum_j q^{(t)}_{k,j}},
\qquad
\Sigma^{(t+1)}_k=\frac{\sum_j q^{(t)}_{k,j}\bigl(x_j-\mu^{(t+1)}_k\bigr)\bigl(x_j-\mu^{(t+1)}_k\bigr)^{\!\top}}{\sum_j q^{(t)}_{k,j}} .
$$

Truncation is realised by projection, $\mathrm{Proj}_\Omega(x)=\arg\min_{y\in\bar\Omega}\|x-y\|_2$, so the importance-sampling estimate becomes

$$
\hat P^{\mathrm{SAIS}}_{\mathcal F}=\frac1{N_2}\sum_{i=1}^{N_2}
\frac{\omega\bigl(\mathrm{Proj}_\Omega(x_i)\bigr)}{\hat h_{\mathrm{opt}}\bigl(\mathrm{Proj}_\Omega(x_i)\bigr)}
\,\mathbb I_{\Omega_{\mathcal F}}\bigl(\mathrm{Proj}_\Omega(x_i)\bigr).
$$

The paper contrasts this route with fitting a mixture by maximising a risk function, which is forced to diagonalise $\Sigma_k$ in high dimension, and notes that the time EM adds is negligible next to network training. The experiments are the inverse conductivity problem of electrical impedance tomography and an inverse source problem in a parabolic system.

## 80: a variational loss has no residual, so use the variance instead

### Why the previous machinery does not transfer

The deep Ritz method minimises a variational energy rather than a pointwise residual, so residual-based adaptivity has nothing to work with. The paper observes that the discretisation error of a variational loss splits into an approximation error and a **statistical (Monte Carlo quadrature) error**, and that the latter dominates when the integrand has low regularity. Its worked example puts $G(u(x))=\frac{1}{\sqrt{2\pi}}e^{-x^2/2\sigma^2}$ on $[-1,1]$, where the relative Monte Carlo error is $C(\sigma N)^{-1/2}$, so $O(1/\sigma)$ uniform samples are needed for $O(1)$ relative accuracy. The target of adaptivity should therefore be **variance reduction in the variational loss**.

The penalised deep Ritz problem is

$$
J(u)=\int_\Omega G\bigl(u(x)\bigr)\,\mathrm dx
+\beta\,\bigl\|B\bigl(x,u(x)\bigr)\bigr\|^2_{\partial\Omega,2},
$$

with importance-sampling estimator

$$
I(u)=\int_\Omega G\bigl(u(x,\theta)\bigr)\mathrm dx
=\mathbb E_p\Bigl[\frac{G(u(X,\theta))}{p(X)}\Bigr]
\approx\frac1{N_v}\sum_{i=1}^{N_v}\frac{G(u(x_i,\theta))}{p(x_i)},
\qquad x_i\sim p .
$$

### A sign-changing integrand: adaptivity cannot reach zero variance

If $G\ge0$, the optimal proposal $p^\star=G/\mu$ with $\mu=\int_\Omega G$ gives zero variance. But a variational integrand need not be non-negative. In that case the optimal choice is

$$
p^\star(x)=\frac{|G(u(x,\theta))|}{\mu},
\qquad \mu=\int_\Omega|G(u(x,\theta))|\,\mathrm dx,
$$

which minimises the variance among all densities by Cauchy-Schwarz but leaves a strictly positive residual:

$$
\sigma_{p^\star}=\Bigl(\int_\Omega|G(u(x,\theta))|\mathrm dx\Bigr)^2
-\Bigl(\int_\Omega G(u(x,\theta))\mathrm dx\Bigr)^2>0 .
$$

This is a structural conclusion: **when the integrand changes sign, adaptivity alone cannot remove the statistical error and the sample count** $N_v$ **remains an independent resource**. It draws a clear boundary around what adaptive sampling can deliver.

### Two networks updated alternately

The density model is the bounded KRnet described on the [[en/computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|density-flow page]]:

$$
p_{\text{bKRnet}}(x,\theta_f)=p_Z\bigl(f_{\text{bKRnet}}(x,\theta_f)\bigr)\,
\bigl|\det\nabla_x f_{\text{bKRnet}}\bigr|,
\qquad Z\sim\mathrm{Unif}([-1,1]^d),
$$

trained by $\min_{\theta_f}D_{\mathrm{KL}}(p^\star\|p_{\text{bKRnet}})$, which reduces to minimising a cross entropy written under importance sampling as

$$
H(p^\star,p_{\text{bKRnet}})\approx
-\frac1{N_v}\sum_{i=1}^{N_v}
\frac{|G(u(x_i,\theta))|\,\log p_{\text{bKRnet}}(x_i,\theta_f)}{\mu\,\tilde p(x_i)},
\qquad x_i\sim\tilde p .
$$

The two networks alternate at adaptivity round $k$:

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

### Two mixture models give the density a floor

The ratio $G/p$ blows up wherever the learned density is small, so the paper puts a floor under it. Model one mixes with the uniform density,

$$
p_{\mathrm{mixture}}(x,\theta_f^{k+1})
=\epsilon\,p_{\text{bKRnet}}(x,\theta_f^{k+1})+(1-\epsilon)\,p_{\mathrm{uniform}}(x)
\ \ge\ \frac{1-\epsilon}{|\Omega|} .
$$

Model two mixes recursively with **every** previous flow,

$$
p_{\mathrm{mixture}}\bigl(x,\{\theta_f^t\}_{t\le k+1}\bigr)
=\sum_{t=1}^{k+1}\epsilon(1-\epsilon)^{k+1-t}p_{\text{bKRnet}}(x,\theta_f^t)
+\frac{(1-\epsilon)^{k+1}}{|\Omega|},
$$

with floor $(1-\epsilon)^{k+1}/|\Omega|$. The three algorithm variants correspond to no mixing, mixture one and mixture two, and four sampling strategies including conventional deep Ritz are compared throughout: a two-dimensional peak problem, a two-peak problem, a problem with a singularity, and a high-dimensional Poisson problem. The error metric is a relative discrete $L_2$ error on a tensor grid, or on uniform samples in high dimension, using a ResNet-like architecture with $\sin^3(x)$ activations.

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

## Coverage check

| Item                                                     | Paper | Status                                                            |
| -------------------------------------------------------- | ----- | ----------------------------------------------------------------- |
| Inner-outer split estimator for the fractional Laplacian | 66    | both identities, both Beta laws, the $r_\epsilon$ floor           |
| Expectation form of the Caputo derivative                | 66    | formula and sampling law (with a prefactor caveat)                |
| Squaring bias and two independent draws                  | 66    | the unbiased construction and its conditions                      |
| Limit-state function and failure probability             | 70    | definitions, safe and failure sets, contrast with refinement      |
| Self-adaptive importance sampling                        | 70    | iterative update, termination, weighted final proposal, estimator |
| Theorem 4.4 with its proof skeleton                      | 70    | three assumptions, bound, split estimate, prior-factor check      |
| Fixed-size resampling and subset simulation              | 73    | role and position, with the verification limits stated            |
| Zero-variance proposal and the mixture model             | 76    | optimal density, KL criterion, EM updates, projection             |
| Statistical error of a variational loss                  | 80    | error split, worked example, the $O(1/\sigma)$ requirement        |
| Residual variance for a sign-changing integrand          | 80    | optimal density and the strictly positive floor                   |
| Two-network alternation and both mixtures                | 80    | both objectives, both density floors                              |

## Sources for this page

- L. Guo, H. Wu, X. Yu, and T. Zhou, _Monte Carlo fPINNs: deep learning method for forward and inverse problems involving high dimensional fractional partial differential equations_, Comput. Methods Appl. Mech. Engrg. 400 (2022), 115523 (preprint [arXiv:2203.08501](https://arxiv.org/abs/2203.08501)).
- Z. Gao, L. Yan, and T. Zhou, [_Failure-informed adaptive sampling for PINNs_](https://doi.org/10.1137/22M1527763), SIAM J. Sci. Comput. 45(4) (2023), pp. A1971-A1994 (preprint [arXiv:2210.00279](https://arxiv.org/abs/2210.00279)).
- Z. Gao, T. Tang, L. Yan, and T. Zhou, [_Failure-informed adaptive sampling for PINNs, Part II: combining with re-sampling and subset simulation_](https://doi.org/10.1007/s42967-023-00312-7), Commun. Appl. Math. Comput. 6 (2024), pp. 1720-1741 (preprint [arXiv:2302.01529](https://arxiv.org/abs/2302.01529)).
- W. Liu, L. Yan, T. Zhou, and Y. Zhou, [_Failure-informed adaptive sampling for PINNs, Part III: applications to inverse problems_](https://doi.org/10.4208/csiam-am.SO-2023-0059), CSIAM Trans. Appl. Math. 5(3) (2024), pp. 636-670.
- X. Wan, T. Zhou, and Y. Zhou, _Adaptive importance sampling for deep Ritz_, Commun. Appl. Math. Comput. 7(3) (2025), pp. 929-953 (preprint [arXiv:2310.17185](https://arxiv.org/abs/2310.17185)).
