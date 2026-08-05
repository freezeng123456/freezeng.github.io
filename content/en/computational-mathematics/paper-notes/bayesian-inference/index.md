---
title: Bayesian Inverse Problems and Data Assimilation
description: Ten papers on one tension - a surrogate makes sampling cheap and makes the posterior wrong
lang: en
translation: computational-mathematics/paper-notes/bayesian-inference
tags:
  - computational-mathematics
  - paper-notes
  - bayesian-inverse-problems
---

This topic holds 10 papers spanning 2019 to 2026. They all attack one concrete tension: **posterior sampling under a PDE constraint needs thousands of forward solves, and replacing the forward model by a surrogate makes sampling cheap while writing the surrogate's error straight into the posterior.**

![Surrogate and sampler refine together](assets/diagrams/tao-zhou-papers/en/bayesian-surrogate-loop.svg)

## The precise shape of the problem

With unknown parameter $m$, forward operator $\mathcal G$ and observations

$$
y=\mathcal G(m)+\eta,\qquad \eta\sim\mathcal N(0,\Sigma_\eta),
$$

Stuart's infinite-dimensional formulation writes the posterior through its density against the prior:

$$
\frac{d\nu}{d\mu_0}(m)=\frac{1}{Z}\exp\bigl(-\Phi(m;y)\bigr),
\qquad
\Phi(m;y)=\tfrac12\bigl\|y-\mathcal G(m)\bigr\|_{\Sigma_\eta}^2 .
$$

Replacing $\mathcal G$ by a surrogate $\widehat{\mathcal G}$ gives an approximate posterior $\widehat\nu$. The difficulty is not the surrogate's global accuracy but its **location**: a surrogate trained on prior samples is accurate on the prior support, while the posterior usually concentrates on a thin subset of it, and when the data carry information beyond the prior that subset can sit in a low-probability region of the prior. The surrogate is therefore least accurate exactly where accuracy matters most.

This observation recurs throughout the topic as a design justification, and its quantitative form is a bound that asks only for **local** accuracy. Define the $\epsilon$-feasible set

$$
\mathcal M(\epsilon)=\bigl\{m:\ \|\mathcal G(m)-\widehat{\mathcal G}(m)\|_\infty\le\epsilon\bigr\},
\qquad
\mathcal M^{\perp}(\epsilon)=\mathcal M\setminus\mathcal M(\epsilon).
$$

Then there are constants $K_1,K_2>0$ with

$$
D_{\mathrm{KL}}\bigl(\widehat\nu\,\|\,\nu\bigr)
\le\Bigl(K_1\epsilon+K_2\,\nu\bigl(\mathcal M^{\perp}(\epsilon)\bigr)\Bigr)^{2},
$$

and the corresponding Hellinger bound drops the square. The bound says: as long as almost all **posterior** mass lies where the surrogate is accurate, the surrogate's behaviour elsewhere does not change the order of the KL distance. Every algorithm in this topic is an implementation of that sentence.

## How the ten papers group

| Close reading                                                                                                                      | Papers                  | What they share                                                               |
| ---------------------------------------------------------------------------------------------------------------------------------- | ----------------------- | ----------------------------------------------------------------------------- |
| [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates\|Multi-fidelity surrogates refined online]] | 34, 37, 49, 79          | an error indicator decides when to call the true model and retrain            |
| [[en/computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering\|Samplers, filters and function-space flows]] | 55, 56, 82, 88, 99, 106 | the sampler itself changes: optimisation proposals, particles, filters, flows |

## Lineage

| Year | No. | What changes                                                                          | What stays                                      |
| ---- | --- | ------------------------------------------------------------------------------------- | ----------------------------------------------- |
| 2019 | 37  | the origin: multi-fidelity polynomial chaos with Metropolis-Hastings                  | an $\ell^\infty$ trigger at an accepted state   |
| 2019 | 34  | the sampler becomes a regularising ensemble Kalman smoother                           | the same low-order multi-fidelity correction    |
| 2020 | 49  | the surrogate becomes a composite multi-fidelity network                              | the same retraining inside a local ball         |
| 2021 | 55  | the sampler becomes randomize-then-optimize, trained once offline                     | training points from an approximate posterior   |
| 2021 | 56  | the sampler becomes Stein variational gradient descent                                | trigger at the particle mean, points reused     |
| 2024 | 79  | the surrogate becomes an operator network with a goal-oriented indicator              | a closed refinement loop with transfer learning |
| 2025 | 82  | the setting becomes sequential data assimilation, error learned offline               | a cheap model plus an explicit error correction |
| 2025 | 99  | variational inference directly in function space                                      | variational rather than sampling                |
| 2026 | 88  | the setting becomes particle-simulation calibration, staged and hierarchical          | network surrogate plus Bayesian inference       |
| 2026 | 106 | the posterior approximation becomes a latent-variable flow, and the prior mean adapts | the same relative-misfit stopping rule as 79    |

## Three recurring design decisions

### The indicator must be cheap enough to evaluate every step

The binding constraint in this family is not accuracy but the **cost of the indicator itself**. The honest local model error

$$
e_M(t)=\Bigl(\int \bigl|\mathcal G(m)-\widehat{\mathcal G}_t(m)\bigr|^2\,\nu_t(dm)\Bigr)^{1/2}
$$

requires a high-dimensional integral and is not implementable. Every paper therefore retreats to a single-point quantity: a relative $\ell^\infty$ error at an accepted chain state (papers 37 and 49), a relative $\ell^2$ error at the particle mean (paper 56), or the data misfit itself at an anchor point (papers 79 and 106). The last of these is goal-oriented: it measures not whether the surrogate is accurate but whether its inaccuracy has damaged the fit.

### New training points must not degenerate

If refinement draws points randomly inside a small ball, the points cluster and the least-squares or training problem degenerates. Three remedies appear across the papers: shrink the ball radius each round as $R\leftarrow\rho R$ (papers 37, 49, 56); impose a separation constraint $\|x'-x\|_2\ge R$ on candidates (paper 56); and run a greedy diversity selection in the surrogate's output space with a coefficient pulling selections back toward the anchor (paper 79). Paper 106 goes the other way and discards the historical data entirely, fine-tuning only on the newly generated local set.

### Surrogate error must not enter the posterior directly

If the whole chain uses the surrogate's acceptance probability, the result is the surrogate's posterior. These algorithms therefore insert the true model at one decisive place: papers 37 and 49 compute the trigger point using the **true-model acceptance probability**, so the point being tested sits closer to the posterior bulk; paper 79 uses the true model to select an anchor from a candidate pool. The cost is a few true solves per round, and it is exactly this step that makes refinement move toward the true posterior instead of the surrogate's.

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]. Per-paper references appear at the end of each close-reading page.
