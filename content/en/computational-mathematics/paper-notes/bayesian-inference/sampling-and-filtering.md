---
title: Samplers, Filters and Function-Space Flows
description: Papers 55, 56, 82, 88, 99 and 106 - replacing the sampler instead of the surrogate
lang: en
translation: computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering
tags:
  - paper-notes
  - bayesian-inverse-problems
  - variational-inference
---

> [!note] Coverage of this page
> Papers **55** (_J. Comput. Math._ 39(6), 2021), **56** (_Comput. Methods Appl. Mech. Engrg._ 386, 2021), **82** (_Comput. Phys. Commun._ 311, 2025), **88** (_Int. J. Mech. Sci._ 313, 2026), **99** (submitted to _Math. Comput._, [arXiv:2411.13277](https://arxiv.org/abs/2411.13277)) and **106** (submitted to _SIAM J. Sci. Comput._, [arXiv:2605.29373](https://arxiv.org/abs/2605.29373)).
>
> Sections 2 to 4 of paper 82 require a subscription. For that paper this page reports only what the abstract, introduction and section openings confirm, and marks anything unverified.

The four papers on the [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|previous page]] fix the sampler and replace the surrogate. The six here do the reverse: the surrogate's role is comparatively stable and what changes is **how one advances on the posterior**.

## 55: putting a surrogate inside an optimisation proposal

### The structure of randomize-then-optimize

RTO turns posterior sampling into repeated solution of a randomly perturbed nonlinear least-squares problem. In whitened coordinates (prior $u\sim\mathcal N(u_{\mathrm{pr}},\Gamma_{\mathrm{pr}})$ with $u=S_{\mathrm{pr}}v+u_{\mathrm{pr}}$), the linearisation point is

$$
v_{\mathrm{ref}}=\arg\min_{v}\tfrac12\|H(v)\|^{2},
$$

and the proposal is

$$
v^{(i)}_{\mathrm{prop}}=\arg\min_{v}\ \tfrac12\bigl\|\bar Q^{T}H(v)-\xi^{(i)}\bigr\|^{2},
\qquad \xi^{(i)}\sim\mathcal N(0,I_n).
$$

Its density is explicit,

$$
\pi_{\mathrm{RTO}}(v)=(2\pi)^{-n/2}\bigl|{\det}\bigl(Q^{T}\nabla H(v)\bigr)\bigr|
\exp\Bigl(-\tfrac12\bigl\|Q^{T}H(v)\bigr\|^{2}\Bigr),
$$

so an independence-proposal Metropolis correction is available. The acceptance ratio collapses to a ratio of weights:

$$
\frac{\pi_{\mathrm{tar}}(v^{(i)}_{\mathrm{prop}})\,\pi_{\mathrm{RTO}}(v^{(i-1)})}
{\pi_{\mathrm{tar}}(v^{(i-1)})\,\pi_{\mathrm{RTO}}(v^{(i)}_{\mathrm{prop}})}
=\frac{w\bigl(v^{(i)}_{\mathrm{prop}}\bigr)}{w\bigl(v^{(i-1)}\bigr)},
$$

$$
w(v)=\bigl|{\det}\bigl(Q^{T}\nabla H(v)\bigr)\bigr|^{-1}
\exp\Bigl(-\tfrac12\|H(v)\|^{2}+\tfrac12\|Q^{T}H(v)\|^{2}\Bigr).
$$

The advantage of RTO is an independent proposal with far lower autocorrelation than a random walk. The price is that every sample needs many evaluations of the forward model **and its Jacobian**, both inside the optimiser and inside the weight.

The paper also restates the scalable RTO of Bardsley and coauthors: take a rank-$r$ reduced singular value decomposition $\nabla f(v_{\mathrm{ref}})=\Psi\Lambda\Phi^{T}$, split $v=\Phi v_r+v^{\perp}$, and then

$$
v^{\perp}=(I_n-\Phi\Phi^{T})\xi,
\qquad
v_r=\arg\min_{z}\bigl\|(\Lambda^{2}+I_r)^{-1/2}z+\Lambda\Psi^{T}f(v^{\perp}+\Phi z)-\Phi^{T}\xi\bigr\|^{2},
$$

which replaces an $n$-dimensional optimisation by an $r$-dimensional one, while the determinant simplifies to

$$
\bigl|{\det}(\tilde Q^{T}\nabla H(v))\bigr|
=\bigl|{\det}(\Lambda^{2}+I_r)^{-1/2}\bigr|\cdot
\bigl|{\det}\bigl(I_r+\Lambda\Psi^{T}\nabla f(v)\Phi\bigr)\bigr| .
$$

### The substantive contribution: training points from an approximate posterior

Once the surrogate is a network, $\nabla_v\mathcal{NN}(v;\theta)$ comes from back-propagation and both the optimisation and the weight become cheap. What actually decides the outcome is where the training points sit: the paper draws them from an **approximate posterior** $\widetilde\pi_{\mathrm{pos}}$ rather than from the prior. The ablation baseline "NN-RTO-pr" is the same algorithm with prior-drawn points.

The difference from the [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|four papers on the previous page]] is that the surrogate here is **trained once offline with no online refinement**. That makes it the static member of the family and the natural contrast case when arguing for the value of online refinement.

Efficiency is measured by CPU-time-adjusted effective sample size:

$$
\mathrm{ESS}=\frac{n_{\mathrm{samps}}}{1+2\sum_{k=1}^{K}\rho(k)} .
$$

### Numerical evidence

A benchmark Darcy-type elliptic inverse problem: recover the permeability $\kappa(x)$ from noisy pointwise pressure measurements with source $f(x)=100\sin(\pi x_1)\sin(\pi x_2)$, solved by bilinear Galerkin finite elements on a uniform $40\times40$ grid. Three methods are compared: direct RTO on the true model, prior-trained NN-RTO-pr, and posterior-trained NN-RTO. The network has 3 hidden layers of 40 neurons, trained with Adam at learning rate $5\times10^{-4}$, and the inner optimisation uses a trust-region-reflective Newton solver. NN-RTO reproduces the one- and two-dimensional posterior marginals of direct RTO already at $N=50$ training points while NN-RTO-pr does not, and that ablation isolates the value of placing training points near the posterior.

The paper proves no theorem; the abstract's "converges to the direct RTO approach" is an argued and empirical statement rather than a proved proposition.

## 56: particle flow combined with online refinement

### Why SVGD needs a gradient

Stein variational gradient descent transports a particle ensemble to the target along a kernelised functional gradient flow. The particle update is $x_i\leftarrow x_i+\varepsilon\varphi(x_i)$ and the optimal direction solves a functional optimisation

$$
\varphi^{*}=\arg\max_{\varphi\in\mathcal S}
\Bigl\{-\tfrac{d}{d\varepsilon}\mathrm{KL}\bigl(q_{[\varepsilon\varphi]}\,\|\,\pi\bigr)\Bigr\}.
$$

The Stein operator identity turns that derivative into an expectation:

$$
-\frac{d}{d\varepsilon}\mathrm{KL}\bigl(q_{[\varepsilon\varphi]}\|\pi\bigr)\Big|_{\varepsilon=0}
=\mathbb E_{x\sim q}\bigl[\mathcal A_{\pi}^{\top}\varphi(x)\bigr],
\qquad
\mathcal A_{\pi}^{\top}\varphi(x)=\nabla_x\log\pi(x)^{\top}\varphi(x)+\nabla_x^{\top}\varphi(x).
$$

Maximising over the unit ball of a reproducing kernel Hilbert space, $\mathcal S=\{\varphi\in\mathcal H^{d}:\|\varphi\|_{\mathcal H^{d}}\le1\}$, gives a closed form,

$$
\varphi^{*}(\cdot)\ \propto\ \mathbb E_{x\sim q}
\bigl[\nabla_x\log\pi(x)\,\kappa(x,\cdot)+\nabla_x\kappa(x,\cdot)\bigr],
$$

with empirical update

$$
x^{(l+1)}_{i}\leftarrow x^{(l)}_{i}+\varepsilon_l Q_l\bigl(x^{(l)}_i\bigr),
\qquad
Q_l(x)=\frac1N\sum_{j=1}^{N}
\Bigl[\nabla_{x^{(l)}_j}\log\pi\bigl(x^{(l)}_j\bigr)\kappa\bigl(x^{(l)}_j,x\bigr)
+\nabla_{x^{(l)}_j}\kappa\bigl(x^{(l)}_j,x\bigr)\Bigr].
$$

The two terms divide the labour cleanly: the first is a kernel-weighted ascent of $\log\pi$ and attracts, the second is a kernel gradient and repels, which is what prevents the particles from collapsing onto a point. The kernel is the radial basis $\kappa(x,x')=\exp(-h\|x-x'\|^{2})$. Step sizes use AdaGrad, and the paper notes explicitly that a constant step can diverge in high dimension.

The obstruction is $\nabla_x\log\pi(x)$: for a PDE-constrained posterior that means one adjoint solve per particle per step. If the gradient is unavailable (legacy code, non-differentiable solver) or simply too expensive, vanilla SVGD is unusable.

### The test point and the selection constraint

The paper uses the particle mean as the design point,

$$
x^{*}=\frac1M\sum_{i=1}^{M}x^{(t+1)}_{i},
$$

and the indicator is the relative $\ell^2$ error there:

$$
\mathrm{err}(x^{*})=\frac{\bigl\|f(x^{*})-\tilde f(x^{*})\bigr\|_{2}}{\bigl\|f(x^{*})\bigr\|_{2}} .
$$

How new training points are chosen deserves separate emphasis: they are taken **from the existing particles rather than from fresh random draws**, subject to a separation constraint,

$$
x_i=\arg\min_{x'\in X^{(t+1)}}\|x'-x^{*}\|_{2}
\quad\text{s.t.}\quad \|x'-x\|_{2}\ge R\ \ \forall x\in\mathcal D .
$$

The constraint exists to prevent degeneracy: without separation from the existing design $\mathcal D$, the points cluster and the training problem becomes ill-conditioned. If no feasible minimiser exists the loop breaks and the radius shrinks, $R\leftarrow\rho R$ with $0<\rho<1$. Online retraining uses transfer learning from the previous network.

The outer loop runs $T$ inner SVGD steps, tests for refinement, and repeats at most $I_{\max}$ times. Total online true-solve cost is $N_{\mathrm{eval}}=\sum_t q_t$, where $q_t$ is the number of points actually accepted at outer step $t$.

### Numerical evidence

Three problems. The first is a two-dimensional "double banana" posterior: a fixed prior-trained network (10 training points, 3 hidden layers of 20 neurons) drives the particles to the **wrong** high-probability region, while the locally refined version spreads them over the true support after 100 iterations. The second is two-dimensional heat source inversion and the third estimates a diffusion coefficient in a time-fractional PDE. Default parameters are $Q=5$, $R=0.2$, tolerance $10^{-2}$, $\rho=0.8$, $I_{\max}=30$ outer and $T=10$ inner iterations, momentum $0.9$; sample quality is scored by maximum mean discrepancy with a shared median-heuristic bandwidth. In the third example the locally refined version reaches error $0.1732$ with 50 online true solves at $n_t=100$, against $0.3968$ for the fixed network under the same conditions; at $n_t=500$ the figures are 80 solves with error $0.1155$ against $0.2941$.

The abstract's "refine the local approximation online without destroying the convergence" states the design goal; the paper proves no corresponding theorem.

## 82: two-tier learning in a time-dependent setting

> [!warning] What could be verified
> Sections 2 to 4 of this paper require a subscription. What follows comes from the abstract, the full introduction and the section openings; it contains **no verifiable formulas**, and the algorithm outline is reconstructed from the abstract and introduction rather than read from an algorithm box.

### The problem

Joint state and parameter estimation for time-dependent parametric PDEs by ensemble Kalman filtering needs one PDE solve per ensemble member per assimilation step. Shrinking the ensemble introduces sampling error in the covariance estimate, which localisation only partly mitigates. Classical reduced-order models rely on linear mode superposition, but the solution manifold of a parametric dynamical system often has slowly decaying Kolmogorov width, so the reduced dimension grows and the benefit evaporates. Data-driven deep reduced-order models escape the linear-subspace restriction, but their offline snapshot generation is itself expensive. That last cost is the paper's target.

### The two tiers

The first tier runs operator inference on a **coarse grid**: the reduced model is given a polynomial form mirroring the structure of the governing equations, and the reduced operators are learned from simulation data by least-squares regression rather than by intrusive projection. Snapshots are generated on a coarse grid $H$ with $H\gg h$, where $h$ is the fine grid of the full-order model. That is simultaneously why the offline stage is cheap and why the reduced model is biased.

The second tier learns the **model error** with networks that map time and parameters to the reduced-order discrepancy, justified by the observation that these errors stay relatively stable across the temporal domain. The corrected model is embedded in the ensemble Kalman filter, with the correction acting on the observation output side so that the Kalman update sees a debiased predicted observation.

The paper's novelty claim is that this is the first operator-inference-based ensemble Kalman approach for simultaneous state and parameter estimation in nonlinear time-dependent PDEs. Error metrics are $\mathrm{rel}_{L^2}(\mu)=\|\mu-\mu^{*}\|_2/\|\mu^{*}\|_2$ and an analogous relative $L^2$ error for the state.

### Where it sits relative to the rest of the topic

Structurally this is paper 34 with the ingredients swapped: a cheap model inside a Kalman method plus an explicit correction for the surrogate's error. The difference is **when** the correction is learned. Papers 34, 37, 49, 56 and 79 all refine online, driven by an indicator; this paper learns the discrepancy **offline** as a function of time and parameters and applies it online with no further refinement, which puts it alongside paper 55 in that respect. The benchmarks are Burgers' equation, the FitzHugh-Nagumo model and advection-diffusion-reaction systems, in a three-way comparison against the uncorrected reduced-order filter and the full-order filter. The abstract reports considerable speedup without loss of accuracy; a specific speedup factor does not appear in the verifiable material.

## 88: hierarchy and staging against non-identifiability

This is the one paper in the topic that is not a PDE inverse problem, and that is exactly what makes it useful for comparison.

### The structure of the problem

Aligning a mesoscopic red blood cell model with experiment is an inverse problem with two structural difficulties. First, the data come from different experimental platforms — interference microscopy for geometry, optical tweezers for stretching, fluctuation and relaxation assays for viscoelasticity — with mutually inconsistent cross-platform discrepancies. Second, the parameters are not simultaneously identifiable from any single dataset. On top of that the forward model is a dissipative particle dynamics simulation, far too slow for direct MCMC.

The inferred parameters are $\vartheta_{\mathrm{in}}=(A_0,v,\mu_{sh},k_b,\eta_m)$: desired surface area, reduced volume, shear stiffness, bending stiffness and membrane viscosity. The simulation outputs are $y=(D_{eq},h_{\max},h_{\min},D_{ax},D_{tr},t_c,W_{fl})$: equilibrium diameter, maximum and minimum thickness, axial and transverse deformation under stretching, relaxation characteristic time and a membrane-fluctuation summary. The observation model lumps all model error and uncertainty into one zero-mean Gaussian term, $y_j=\text{model}_j(x;\vartheta)+\sigma\epsilon_j$.

### Three structural decisions

**Dynamic annealing fixes the stress-free baseline.** Cells at different reduced volumes $v$ are deformed from the biconcave shape at $v=0.64$; at larger $v$ the fixed triangular mesh retains large residual bond forces, so whenever $|v-v_{\mathrm{desired}}|$ exceeds a tolerance the equilibrium bond length is moved toward the current bond length. This step defines the stress-free baseline that all four experiments are referenced to.

**Eight surrogate networks.** Four experiments times two cell populations (healthy cells and _Plasmodium falciparum_-infected cells), each with three hidden layers, `tanh` activations, mean-squared-error loss and Adam with a step learning-rate schedule, each trained on 10,000 simulation results; widths were chosen per experiment guided by sensitivity analysis. The abstract reports sub-$10^{-2}$ prediction errors.

**Two-stage hierarchical inference.** Stage I performs hierarchical inference over the single-level models for the equilibrium and stretching tests, giving stable distributions for the geometric parameters and the shear modulus. Stage II carries that information forward, adds the fluctuation and relaxation tests, and identifies the full parameter set including the viscoelastic parameters. Hyperparameters encode information sharing across datasets. Experimental anchors come from the published literature: for healthy cells $D_{eq}=7.82\pm0.62\,\mu\mathrm m$, $h_{\max}=2.58\pm0.27\,\mu\mathrm m$, $h_{\min}=0.81\pm0.35\,\mu\mathrm m$; for infected cells, inferred from surface-area and volume measurements, $D_{eq}=6.9\,\mu\mathrm m$ and $h_{\max}=h_{\min}=3.2\,\mu\mathrm m$.

The reported outcome is statistically robust posteriors, pathological cells inferred to be stiffer and more viscous, and cross-platform data fusion mitigating the multi-source uncertainty that defeats single-level inference. The paper proves no mathematical theorem.

> [!note] Title difference
> The preprint is titled _An RBC-MsUQ Framework for Red Blood Cell Morpho-Mechanics_, with the same authors and abstract; the journal title is the one listed on the homepage, and this site records the journal version.

## 99: variational inference directly in function space

### The central difficulty is measure-theoretic

Discretising first and then doing variational inference makes the algorithm degrade as the mesh is refined. Working directly in function space fixes that but raises a new obstruction: probability measures on infinite-dimensional spaces are generically **mutually singular**, so an arbitrary flow transformation of the prior may have no Radon-Nikodym derivative with respect to it, and then the KL loss is not even defined.

The infinite-dimensional Bayes formula reads $\dfrac{d\mu}{d\mu_0}(u)=\dfrac{1}{Z_\mu}\exp(-\Phi(u))$, and every measure in the approximating family $\mathcal M(\mathcal H_u)$ is required to be equivalent to the prior $\mu_0$ so that

$$
D_{\mathrm{KL}}(\nu\|\mu)=\int_{\mathcal H_u}
\Bigl[\ln\frac{d\nu}{d\mu_0}(u)-\ln\frac{d\mu}{d\mu_0}(u)\Bigr]\nu(du)
$$

makes sense.

### Theorem 2.3: equivalence with an explicit Radon-Nikodym derivative

For a composition $f_\theta=f^{(N)}_{\theta_N}\circ\cdots\circ f^{(1)}_{\theta_1}$ whose layers are bijective and satisfy the paper's conditions, $\mu_{f_\theta}\sim\mu_0$ and

$$
\frac{d\mu_{f_\theta}}{d\mu_0}\bigl(f_\theta(u)\bigr)
=\prod_{n=1}^{N}\bigl|{\det}_{1}\bigl(Df^{(n)}_{\theta_n}(u_{n-1})\bigr)\bigr|^{-1}
\exp\Bigl(\tfrac12\bigl\langle f_\theta(u)-u,\ f_\theta(u)-u\bigr\rangle_{\mathcal H}
+\bigl\langle u,\ u-f_\theta(u)\bigr\rangle_{\mathcal H}\Bigr),
$$

where ${\det}_1$ is the Carleman-Fredholm determinant, the correct object in infinite dimensions. The loss decomposition actually used pushes every expectation back onto the **prior**:

$$
D_{\mathrm{KL}}(\mu_{f_\theta}\|\mu)
=\mathbb E_{\mu_0}\ln\frac{d\mu_{f_\theta}}{d\mu_0}\bigl(f_\theta(u)\bigr)
-\mathbb E_{\mu_{f_\theta}}\ln\frac{d\mu}{d\mu_0}(u),
$$

and the prior can be sampled directly, so the Monte Carlo gradient is computable: draw $N$ samples from $\mu_0$, push them through the flow and update with

$$
\nabla_{\theta_k}L(\theta_k)\approx
\frac1N\sum_{i=1}^{N}\nabla_{\theta_k}\ln\frac{d\mu_{f_{\theta_k}}}{d\mu_0}\bigl(f_{\theta_k}(u_i)\bigr)
-\frac1N\sum_{i=1}^{N}\nabla_{\theta_k}\ln\frac{d\mu}{d\mu_0}\bigl(f_{\theta_k}(u_i)\bigr),
$$

using Adam in the experiments.

### Four flow layers and why low rank is mandatory

The paper gives four concrete transformations satisfying its conditions: a functional Householder flow (linear, one-dimensional image), a functional projected transformation flow (linear, $M$-dimensional image and therefore more expressive than Householder), a functional planar flow (nonlinear) and a functional Sylvester flow (nonlinear). Invertibility of the linear and nonlinear families is established by two lemmas.

The functional planar flow has the form $f_n(u)=u+u_n h\bigl(\langle u,w_n\rangle_{\mathcal H_u}+b\bigr)$ with the low-rank parameterisation

$$
u_n=\sum_{i=1}^{r}\lambda_i\alpha_i\phi_i,
\qquad
w_n=\sum_{i=1}^{r}\lambda_i\beta_i\phi_i
$$

built from the eigenpairs $\{\lambda_i,\phi_i\}$ of the prior covariance $\mathcal C_0$, with $\{\alpha_i,\beta_i\}$ trainable. The paper stresses that low rank here is not a computational convenience but a **theoretical necessity**: without it, measure equivalence fails and the KL loss is no longer finite. That is the substantive difference from low-rank tricks in neural operators.

The paper also proves a discretisation-invariance proposition (if $\mathcal H_u$ embeds continuously in $C(D)$, all layers of all four flows are discretisation invariant) and introduces a conditional variant: for a fixed forward model, different data induce different posteriors, so the unconditional version must be retrained per dataset, while the conditional version conditions the flow on the data and handles measurement vectors of varying dimension.

The experiments cover three inverse problems — a one-dimensional smooth equation, two-dimensional steady-state Darcy flow, and electrical impedance tomography — with pCN as an MCMC baseline, targeting agreement with the theory, efficiency relative to pCN, and empirical discretisation invariance (the same flow trained at different discretisation levels behaves the same).

> [!note] Version drift
> Preprint v1 reports two inverse problems (the smooth equation and steady-state Darcy flow); v2 and v3 add electrical impedance tomography. This page follows v3.

## 106: a latent-variable flow with an adaptive prior

### Three coupled failure modes

In high-dimensional PDE-governed inversion, non-Gaussian and often multimodal posteriors defeat Gaussian-approximation samplers (unscented and ensemble Kalman inversion) and strain even pCN; a neural surrogate pre-trained on prior samples is out of distribution once the posterior concentrates elsewhere; and the prior mean is usually misspecified, while hand-tuning it is exactly the manual intervention one wants to remove. Standard normalizing flows cannot help with dimensionality because they are bijections and therefore dimension preserving.

### Variational Flow

The architecture is VAE-style nonlinear dimension reduction (latent $z\in\mathbb R^{k}$, data $x\in\mathbb R^{d}$, $k<d$) augmented by **dual** flows: one on the latent prior and one turning the encoder from a diagonal Gaussian into a conditional normalizing flow. The sampling map is

$$
z=f^{-1}_{\mathrm{pr},\beta}(v),\ v\sim\mathcal N(0,I)
\quad\Longrightarrow\quad
\xi=\mu_{\mathrm{de},\theta}(z)+\sigma_{\mathrm{de},\theta}(z)\odot\epsilon,\ \epsilon\sim\mathcal N(0,I),
$$

and the target is the unnormalised posterior $\hat p(x)=\exp(-\Phi(\xi,y))\pi_0(\xi)$.

The paper proves a strict evidence-lower-bound improvement: under either of two conditions the model attains a strictly higher bound than a standard VAE, and when both hold, the flow prior and the conditional flow encoder each contribute a strictly positive improvement. The proof uses

$$
\mathrm{ELBO}=\text{const}-\mathbb E_{p_x}D_{\mathrm{KL}}\bigl(q_{z|x,\alpha}\,\|\,p_{z|x,\theta^{*},\tilde\beta}\bigr)
$$

in two steps: first fix the optimal VAE encoder and decoder and replace the latent prior by a flow prior, then fix the decoder and flow prior and extend the encoder from a diagonal Gaussian to a conditional flow. The paper states explicitly that the adaptive loop itself has no convergence guarantee.

### Momentum updating of the prior mean

$$
\mu^{(k,i)}_{\mathrm{prior}}=\alpha\,\mu^{(k,i)}_{\mathrm{post}}+(1-\alpha)\,\mu^{(k-1)}_{\mathrm{prior}},
\qquad
\mu^{(k,i)}_{\mathrm{post}}=\frac1M\sum_{j=1}^{M}\xi^{(j)},
\quad \xi^{(j)}\sim p^{(k,i-1)}_{\mathrm{VF}}(\xi).
$$

Two time scales operate: the posterior mean is re-estimated every epoch, while the anchor $\mu^{(k-1)}_{\mathrm{prior}}$ stays fixed within a stage and is updated at the end of it. The covariance is **deliberately held** at $\Sigma_0$ to prevent mode collapse. Smaller $\alpha$ means stronger regularisation, $\alpha=1$ removes it, prior work cited in the paper finds $\alpha\approx0.5$ close to optimal, and initialisation is $\mu_0=0$, $\Sigma_0=I$.

### Surrogate fine-tuning and aggressive data replacement

Posterior samples are perturbed before being used as training data,

$$
\hat\xi^{(j)}=\xi^{(j)}_{\mathrm{post}}+\gamma\,\nu^{(j)},
\qquad \nu^{(j)}\sim\mathcal N(0,I),\ \gamma>1,
$$

with the perturbation applied in Karhunen-Loève coefficient space where the prior is standard normal, so $\gamma$ is measured in prior standard deviations. The paper contrasts itself with paper 79 explicitly: instead of accumulating a growing dataset or running expensive greedy filtering, it **discards** the old data and fine-tunes the Fourier neural operator only on the newly generated local set.

The stopping criterion is the relative change in the data misfit at the prior mean:

$$
\frac{\bigl|\Phi\bigl(\mu^{(k-1)}_{\mathrm{prior}},y\bigr)-\Phi\bigl(\mu^{(k)}_{\mathrm{prior}},y\bigr)\bigr|}
{\Phi\bigl(\mu^{(k-1)}_{\mathrm{prior}},y\bigr)}<\epsilon .
$$

This is the same class of goal-oriented quantity as $e_D$ in paper 79; only the test point differs, an anchor selected by the true model there against the evolving prior mean here.

### Numerical evidence

Four problems. The first is a 100-dimensional Rosenbrock inverse problem used to test posterior approximation quality alone, so the prior-updating module is **deliberately switched off**, and two-dimensional marginals are compared against a VAE, MCMC, SVGD and unscented Kalman inversion; the model is reported most accurate and significantly better than a vanilla VAE, which the authors read as empirical support for the bound theorem. The other three are one-dimensional Darcy flow, two-dimensional Darcy flow and two-dimensional Navier-Stokes (recovering the initial vorticity from the vorticity field at $T=1$), each at Karhunen-Loève truncation $d\in\{32,64\}$ and noise $\delta\in\{1\%,5\%,10\%\}$, averaged over three runs, against pCN, SVGD and unscented Kalman inversion with both full-order and operator surrogates. The reported outcomes: in one-dimensional Darcy the method beats all baselines, most clearly at high noise; in two-dimensional Darcy it is best in most medium- and high-noise cases, and when the reference field is drawn out of distribution from a uniform law the prior updating is what lets it find the true posterior, while pCN and SVGD stay near the initial $\mathcal N(0,I)$ prior and get the scale of the field visibly wrong; in Navier-Stokes it attains the lowest inversion error at every tested truncation dimension and noise level.

> [!note] Quantitative tables
> The specific error values in this paper could not be transcribed reliably from the available material, so this page reports only the confirmable qualitative conclusions and experimental configurations.

## Side-by-side comparison

| No. | Posterior approximation                  | Surrogate                        | Refinement                   | Theory                                   |
| --- | ---------------------------------------- | -------------------------------- | ---------------------------- | ---------------------------------------- |
| 55  | optimisation-based independence proposal | feed-forward network             | none (offline once)          | none                                     |
| 56  | particle flow                            | feed-forward network             | online, particle mean        | none                                     |
| 82  | ensemble Kalman filter                   | operator-inference reduced model | none (error learned offline) | none                                     |
| 88  | staged hierarchical MCMC                 | eight feed-forward networks      | none                         | none                                     |
| 99  | function-space variational flow          | no surrogate                     | not applicable               | measure equivalence, explicit derivative |
| 106 | latent-variable variational flow         | Fourier neural operator          | online, replacement          | strict evidence-bound improvement        |

## Coverage check

| Item                                                | Paper | Status                                                                                   |
| --------------------------------------------------- | ----- | ---------------------------------------------------------------------------------------- |
| RTO proposal, density and acceptance weight         | 55    | linearisation point, proposal, density, weight, scalable variant                         |
| Ablation on posterior-drawn training points         | 55    | contrast with the prior-trained baseline and its conclusion                              |
| Stein operator and the closed-form direction        | 56    | functional optimisation, identity, closed form, attraction and repulsion                 |
| Point selection with a separation constraint        | 56    | constraint form, degeneracy motivation, radius shrinking                                 |
| Two-tier reduction plus an error network            | 82    | coarse-grid operator inference, error network, output-side correction (with caveats)     |
| Staged hierarchical architecture and annealing      | 88    | parameters and outputs, annealing, eight surrogates, stage split                         |
| Measure equivalence, the derivative and four flows  | 99    | theorem conclusion, loss decomposition, necessity of low rank, discretisation invariance |
| Latent flow, prior updating, aggressive replacement | 106   | sampling map, momentum update, perturbation, stopping rule                               |

## Sources for this page

- L. Yan and T. Zhou, [_An acceleration strategy for randomize-then-optimize sampling via deep neural networks_](https://doi.org/10.4208/jcm.2102-m2020-0339), J. Comput. Math. 39(6) (2021), pp. 848-864 (preprint [arXiv:2104.06285](https://arxiv.org/abs/2104.06285)).
- L. Yan and T. Zhou, [_Stein variational gradient descent with local approximations_](https://doi.org/10.1016/j.cma.2021.114087), Comput. Methods Appl. Mech. Engrg. 386 (2021), 114087 (preprint [arXiv:2104.06276](https://arxiv.org/abs/2104.06276)).
- Y. Wang, L. Yan, and T. Zhou, [_Deep learning-enhanced reduced-order ensemble Kalman filter for efficient Bayesian data assimilation of parametric PDEs_](https://doi.org/10.1016/j.cpc.2025.109544), Comput. Phys. Commun. 311 (2025), 109544.
- S. Wang, L. Ma, L. Guo, X. Li, and T. Zhou, [_Multi-stage uncertainty quantification framework for red blood cell morpho-mechanics_](https://doi.org/10.1016/j.ijmecsci.2026.111352), Int. J. Mech. Sci. 313 (2026), 111352 (preprint [arXiv:2508.06852](https://arxiv.org/abs/2508.06852)).
- Y. Zhao, H. Lu, J. Jia, and T. Zhou, _Functional normalizing flow for statistical inverse problems of partial differential equations_, [arXiv:2411.13277](https://arxiv.org/abs/2411.13277), submitted to Math. Comput.; reference implementation [jjx323/FunctionalNormalizingFlow](https://github.com/jjx323/FunctionalNormalizingFlow).
- Y. Wang, X. Wang, K. Tang, X. Wan, T. Zhou, and C. Yang, _Deep adaptive dimension reduction for Bayesian inference in inverse problems_, [arXiv:2605.29373](https://arxiv.org/abs/2605.29373), submitted to SIAM J. Sci. Comput.
