---
title: Scientific Machine Learning
description: Twenty-two papers on the same three design freedoms - sample placement, spectral content, and where constraints live
lang: en
translation: computational-mathematics/paper-notes/scientific-machine-learning
tags:
  - computational-mathematics
  - paper-notes
  - scientific-machine-learning
---

This is the largest of the seven topics, 22 papers concentrated after 2022. They do not evaluate neural networks as a monolith; they separate three design freedoms that can be controlled independently: **where the collocation points go, which frequencies the network can represent, and in what form constraints enter the loss.**

![Three places uncertainty enters operator learning](assets/diagrams/tao-zhou-papers/en/operator-learning-uq.svg)

## Why those three freedoms deserve separating

A physics-informed loss is usually written

$$
\mathcal L(\theta)=\underbrace{\frac{1}{N_c}\sum_i\bigl|r(x_i;\theta)\bigr|^2}_{\text{equation residual}}
+\lambda\underbrace{\frac{1}{N_b}\sum_j\bigl|b(x_j;\theta)\bigr|^2}_{\text{boundary residual}} .
$$

Three choices hide inside it, usually left at their defaults. Where the $x_i$ come from (most implementations draw them uniformly at random); which frequencies the function class parameterised by $\theta$ can represent (most implementations take a standard fully connected network, which is biased toward low frequencies); and how $\lambda$ is set (most implementations tune it by hand). Every paper in this topic challenges one of those defaults and supplies an explicit target to drive it.

## Five close readings

| Close reading                                                                                                                                       | Papers                    | Default being challenged              |
| --------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------- | ------------------------------------- |
| [[en/computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns\|Where to put the samples]]                      | 66, 70, 73, 76, 80        | uniform random collocation            |
| [[en/computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities\|Densities as invertible maps]]              | 62, 64, 72, 87            | unknowns as grid values               |
| [[en/computational-mathematics/paper-notes/scientific-machine-learning/variational-and-basis-networks\|How constraints enter the loss]]             | 60, 90, 102               | constraints as weighted penalties     |
| [[en/computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning\|Operator learning with UQ]]             | 75, 95, 98, 107           | a prediction is a single number       |
| [[en/computational-mathematics/paper-notes/scientific-machine-learning/spectral-bias-and-generative-solvers\|Spectral bias and generative solvers]] | 81, 89, 94, 101, 103, 105 | all frequencies equally representable |

## Core idea of each paper

### Where to put the samples

- **66 (Monte Carlo fPINNs)** writes the fractional Laplacian and the Caputo derivative as expectations, estimated by an inner-outer split with two Beta laws, so the fractional operator is never discretised. The key supplement is to multiply two independent estimates instead of squaring one, removing the variance bias in the residual loss.
- **70 (failure-informed sampling, Part I)** declares "the residual exceeds a tolerance" a failure event, uses the failure probability as a posterior error indicator, and generates new points in the failure region by self-adaptive importance sampling. Its theorem writes both prescribed tolerances directly into an $L^2$ error bound.
- **73 (Part II)** keeps the training set at constant size and shifts its composition from uniform to adaptive under cosine annealing, replacing the proposal with subset simulation to handle multimodal failure regions.
- **76 (Part III)** targets inverse problems, where a state network and a coefficient network train jointly and the residual can split into disconnected regions, so the proposal becomes a truncated Gaussian mixture fitted by EM with truncation realised by projection.
- **80 (adaptive importance sampling for deep Ritz)** has no pointwise residual to work with and therefore targets the **statistical error**: treat the variational integrand as an unnormalised density, learn it with a bounded KRnet, and use importance sampling. A structural conclusion is that a sign-changing integrand leaves a strictly positive residual variance, so adaptivity cannot substitute for sample count.

### Densities as invertible maps

- **62 (normalizing field flows)** builds a bijection between a Gaussian reference field with truncated Karhunen-Loève structure and the target stochastic field, with the expansion coefficients parameterised by networks and training by maximum likelihood on scattered measurements; adding an equation residual gives the physics-informed version.
- **64 (temporal normalizing flow)** solves time-dependent Fokker-Planck equations, whose solution is a density. Treating time as an extra dimension breaks normalisation because mass is not conserved along the time axis; pinning the latent time to the real time collapses the Jacobian to its spatial part, so the flow is time-**conditioned**.
- **72 (adaptive deep density approximation for fractional Fokker-Planck)** offers two routes for the nonlocal operator: a Monte Carlo estimate, or a Gaussian radial-basis auxiliary model whose fractional Laplacian is a closed-form confluent hypergeometric function, tied to the flow by a consistency penalty.
- **87 (bounded KRnet)** builds an exactly invertible coupling layer on $[-1,1]$ from a piecewise-quadratic cumulative distribution function, giving a density model with bounded support. The price is that the density is piecewise linear and only $C^1$, so second-order equations must be recast as first-order systems.

### How constraints enter the loss

- **60 (augmented Lagrangian deep learning)** moves the essential boundary condition from a penalty to a multiplier. The technical obstruction is that the infinite-dimensional multiplier update cannot act on network parameters, which the paper resolves by a least-squares projection in parameter space.
- **90 (adaptive neural network basis)** freezes randomly preset hidden layers so training becomes linear least squares. For localised loss of regularity it isolates a peak in a subdomain by non-overlapping decomposition and recentres and rescales the local basis; the scale coefficient comes from brute-force search over the integers 1 to 10, affordable because each trial is a small least-squares solve.
- **102 (decoupled divergence-free network basis)** makes the divergence condition hold by construction through $u=\mathbf{curl}\,\phi$ and supplies advection identities in two and three dimensions, decoupling velocity and pressure into two sequential least-squares problems; the single-variable fourth-order formulation of three-dimensional Navier-Stokes is the most substantive step.

### Operator learning with uncertainty quantification

- **75 (IB-UQ)** quantifies uncertainty with an information bottleneck. A confidence-aware encoder $z=\mathrm{diag}(m(x))\bar z(x)+\mathrm{diag}(\mathbf 1-m(x))z_0$ lets the latent variable collapse to pure noise when the input is far from the training distribution; replacing $I(Z;X)$ by $I(\tilde Z;\tilde X)$ over a flattened input distribution is the decisive change for out-of-distribution behaviour.
- **95 (deep set based operator learning)** keeps the conditional variational-autoencoder bound but replaces the confidence gate with a permutation-invariant set encoder, so sensor count and placement may vary.
- **98 (LVM-GP)** follows the confidence-aware encoder idea while replacing the noise term with a Gaussian process prior and the decoder mean with a neural operator, giving uncertainty with correlation structure.
- **107 (operator learning for Fokker-Planck equations)** generalises "one equation, one flow" into learning the **solution operator** over initial conditions, in contrast to the single-equation setting of papers 64, 72 and 87.

### Spectral bias and generative solvers

- **81 (frequency-adaptive multi-scale networks)** notes that the down-scaling factors of a multi-scale network are normally fixed a priori at $2^{i-1}$ while the ideal scales depend on the target's frequency content. The paper proves two bounds (down-scaling shrinks the $kh$ factor in the numerator; a band-limited target embedded with its true dominant modes has a bound with no explicit frequency dependence), proposes a hybrid linear-plus-sinusoidal feature embedding, and captures the dominant modes a posteriori from the discrete Fourier transform of the current solution.
- **94 (deep learning for the semi-classical limit)** applies this directly to the Schrödinger equation with a smooth potential near the semi-classical limit, whose localised solution oscillates in space and time on wavelength $O(\varepsilon)$: precisely the multi-scale target.
- **101 (overcoming spectral bias via cross-attention)** attacks the same problem from another angle, keeping a fixed bank of multiscale random Fourier features and learning input-dependent attention weights over it instead of rebuilding the network around captured modes.
- **89 (energy-based diffusion generator)** samples a Boltzmann distribution with known unnormalised energy. It is a generative-model paper rather than a PDE solver, but it shares the stance of representing a distribution by an invertible or sampleable model.
- **103 (PI-DOSnet)** writes the structure of operator splitting into the network for evolution equations.
- **105 (FLUID)** carries the random-field flow route of paper 62 into a unified inference framework for dynamics, replacing the Karhunen-Loève-structured reference field with conditional flows.

> [!note] Coverage status
> Papers 66, 70, 73, 76, 80, 62, 64, 72, 87, 60, 90, 102, 75, 81 and 94 have close-reading content checked equation by equation. The close readings for 89, 95, 98, 101, 103, 105 and 107 are still being filled in: public material for these is limited, and this page gives only the core ideas confirmable from titles, abstracts and cross-references in neighbouring papers, without expanding the formulas.

## Three transferable judgements

### A sampling distribution is an object that a target can drive

Papers 70, 73, 76 and 80 share one structure: define an estimable scalar target (failure probability, variance of the variational loss), then approximate its optimal density with a sampleable model (truncated Gaussian, mixture, bounded flow). This turns "add more points" into "add points from which density".

### Structural constraints belong in the architecture

Normalisation comes from the change of variables (papers 64, 72, 87), positivity from the push-forward, the initial condition from the factor $t$ in a coupling layer (paper 72), the divergence condition from the $\mathbf{curl}$ ansatz (paper 102), and the essential boundary condition from a multiplier (paper 60). Each guarantee removes one weight that would otherwise need tuning.

### Frequency content should be measured, not assumed

The approach in paper 81 deserves separate emphasis: it does not assume the target's frequency distribution but takes a discrete Fourier transform of the **current network solution**, keeps the coefficients of largest modulus as next round's features, and decides how many modes to keep by an energy-fraction condition

$$
\sum_{j=1}^{N_0}\bigl|\hat f_{\text{net},k_j}\bigr|^2\ \ge\ (1-\delta)\,\|f_{\text{net}}\|^2_{L^2}.
$$

That is the same design as "let an error indicator decide when to solve exactly" in the [[en/computational-mathematics/paper-notes/bayesian-inference/index|Bayesian topic]].

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]; per-paper references appear at the end of each close-reading page.
