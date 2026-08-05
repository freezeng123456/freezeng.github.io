---
title: Bayesian Inverse Problems and Data Assimilation
description: Ten papers on one tension - a surrogate makes sampling cheap and makes the posterior wrong - and the theorems that make that quantitative
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

## Verification tier

How much detail the two close-reading pages carry is decided entirely by how much of each source could be checked. The summary is below; the per-paper gaps are recorded on the pages themselves.

| No. | Source checked                                                     | Tier                       |
| --- | ------------------------------------------------------------------ | -------------------------- |
| 34  | arXiv:1809.08931 full text plus journal abstract                   | Full text                  |
| 37  | arXiv:1807.00618 full text                                         | Full text                  |
| 49  | arXiv:1911.08926 full text                                         | Full text                  |
| 55  | arXiv:2104.06285 full text                                         | Full text                  |
| 56  | arXiv:2104.06276 full text                                         | Full text                  |
| 79  | arXiv:2310.17844 v3 full text                                      | Full text                  |
| 82  | ScienceDirect abstract, complete introduction, section openings    | **Abstract and introduction only** |
| 88  | arXiv:2508.06852 full text                                         | Full text                  |
| 99  | arXiv:2411.13277 v3 full text                                      | Full text                  |
| 106 | arXiv:2605.29373 v1 full text                                      | Full text                  |

Nine are full text. One, paper 82, is paywalled with no preprint, so **not a single equation from its Sections 2 to 4 is verifiable**; its close reading covers only what the abstract and introduction support and says where the line falls.

Full text does not mean every number could be transcribed, and this is the topic's largest remaining gap. **Only paper 34 left a complete, legible set of error and cost figures.** Paper 56 left a quantitative comparison for one example. Papers 37, 49, 55, 79, 88 and 99 left experimental configurations and qualitative conclusions. All of paper 106's tables were destroyed in extraction. The experiments sections below therefore give the reproducible configuration first, then say plainly which conclusions have no numbers behind them.

## Turning surrogate error into posterior error

"A surrogate only has to be accurate near the posterior" is not a slogan; it has three layers of quantitative support, in decreasing order of strength: a classical bound that demands global accuracy, a bound that demands only local accuracy, and the single bound that carries the error all the way to what the sampler actually returns. The differences between them decide what this family of algorithms is entitled to claim.

### The global bound: uniformly close potentials suffice

Stuart's approximation theorem is where every surrogate paper starts. Let $\mu$ and $\mu^{N}$ both be absolutely continuous with respect to the prior $\mu_0$, with potentials $\Phi$ and $\Phi^{N}$ satisfying, **uniformly in $N$**, a lower bound (for every $\varepsilon,r>0$ there is $M$ with $\Phi(u;y)\ge M-\varepsilon\|u\|_X^{2}$) and an upper bound on bounded sets. If there is $\psi(N)\to0$ with

$$
\bigl|\Phi(u)-\Phi^{N}(u)\bigr|\ \le\ K\exp\bigl(\varepsilon\|u\|_X^{2}\bigr)\,\psi(N),
$$

then there is $C$ independent of $N$ with $d_{\mathrm{Hell}}(\mu,\mu^{N})\le C\psi(N)$. The same conclusion has a version stated directly on the forward operator: replace the condition on $\Phi$ by one on $\mathcal G-\mathcal G^{N}$, provided both satisfy the corresponding exponential-growth assumption uniformly in $N$.

**Why Hellinger and not total variation or KL.** Because two measures absolutely continuous with respect to a common Gaussian reference and $\varepsilon$ apart in Hellinger have expectations of polynomially bounded functionals $O(\varepsilon)$ apart — in particular their means and covariance operators are $O(\varepsilon)$ apart. A Hellinger bound therefore translates directly into error bars on the quantities people actually report, which the other two distances do not.

Stuart also separates the error into two kinds: representing $u$ in a finite-dimensional basis, and approximating the forward operator. **All ten papers here address only the second.**

The price of this bound is the word "uniformly": $\psi(N)$ must hold for every $u$. Making $\psi(N)$ small requires a globally accurate surrogate, and global accuracy is precisely the cost this family refuses to pay. Hence the second bound.

### The local bound: the $\epsilon$-feasible set

Split the parameter space into the part where the surrogate meets tolerance and its complement, and measure both with the **posterior**:

$$
\mathcal M(\epsilon)=\bigl\{m:\ \|\mathcal G(m)-\widehat{\mathcal G}(m)\|_\infty\le\epsilon\bigr\},
\qquad
\mathcal M^{\perp}(\epsilon)=\mathcal M\setminus\mathcal M(\epsilon).
$$

**Theorem (paper 37's Theorem 2; the same statement is paper 79's Theorem 2.1, both citing the same source).** Assume the forward operator is uniformly bounded, $\sup_{m}\|\mathcal G(m)\|<\infty$, with the surrogate satisfying the same bound uniformly in its accuracy parameter, and that the observational noise is i.i.d. Gaussian. Then for a given $\epsilon>0$ there exist constants $K_1,K_2>0$ with

$$
D_{\mathrm{KL}}\bigl(\widehat\nu\,\|\,\nu\bigr)
\le\Bigl(K_1\epsilon+K_2\,\nu\bigl(\mathcal M^{\perp}(\epsilon)\bigr)\Bigr)^{2}.
$$

**Corollary (Hellinger version, which paper 37 attributes to its reference [8]).** Under the same hypotheses the square disappears:

$$
d_{\mathrm{Hell}}\bigl(\widehat\nu,\nu\bigr)\le K_1\epsilon+K_2\,\nu\bigl(\mathcal M^{\perp}(\epsilon)\bigr).
$$

The two are consistent in order: Hellinger distance is controlled by the square root of the KL divergence, so a KL bound of order $\epsilon^{2}$ gives a Hellinger bound of order $\epsilon$, which is exactly the shape above. The KL is taken with the approximation first and the truth second; Hellinger distance is symmetric, so no ordering question arises there.

**What the two terms are.** The first term $K_1\epsilon$ is the price paid on the feasible set, where the difference of potentials is controlled by $\epsilon$. The second term $K_2\nu(\mathcal M^{\perp}(\epsilon))$ is the price paid on the complement, where the surrogate may be arbitrarily bad, so the integrand can only be held down by the uniform-boundedness assumption and then multiplied by the **posterior** mass of that region. That is why an apparently gratuitous uniform bound is needed. (This decomposition is this page's reading of the bound's structure; both papers state only the conclusion.)

**What it licenses.** The second term uses posterior measure, not prior measure and not Lebesgue measure. However badly the surrogate behaves over most of the prior support, that behaviour never enters the bound as long as those regions carry little posterior mass. "Locally accurate near the posterior" thus turns from an intuition into a target one can pursue: drive $\epsilon$ down and drive $\nu(\mathcal M^{\perp}(\epsilon))$ down. Paper 37 writes the algorithmic goal out explicitly — **if sampling is good enough that $\nu(\mathcal M^{\perp}(\epsilon))\le\epsilon$, the KL distance is characterised entirely by $\epsilon^{2}$.**

> [!warning] The bound is a theorem; the convergence is not
> The bound above is rigorously proved. But "adaptive refinement drives $\nu(\mathcal M^{\perp}(\epsilon))$ to zero" is **an argument, not a theorem, in every paper in this topic**. The argument runs: once a candidate point lands in $\mathcal M^{\perp}(\epsilon)$ the algorithm refines in its neighbourhood, so the bad set decays asymptotically. No paper says when or how fast it goes to zero, and no paper measures the bad set's measure in its experiments. Do not lend the bound's rigour to the convergence claim when reporting these results.

### End to end: the one bound that reaches the sampler's output

Both bounds above compare **two posterior measures**. Neither says how far what the algorithm hands back — samples, an ensemble, a Gaussian approximation — sits from the true posterior. Of the ten papers only 79 closes that step, and only in the linear case.

There are three hypotheses: (3.2) for any $\epsilon$ the *linear* neural operator $\widehat G:\mathbb R^{N_m}\to\mathbb R^{N_y}$ can be trained so that $\|\widehat G-G\|_2<\epsilon$; (3.3) $\|G\|_2<H$; (3.4) $G^{T}\Sigma_\eta^{-1}G\succ0$ with $\|G^{T}\Sigma_\eta^{-1}G\|_2>C_1$. Lemma 3.5 then gives the surrogate-side counterpart a positive lower bound, $\|\widehat G^{T}\Sigma_\eta^{-1}\widehat G\|_2>C_2$, through the perturbation estimate $\|G^{T}\Sigma_\eta^{-1}G-\widehat G^{T}\Sigma_\eta^{-1}\widehat G\|_2\le2\epsilon H\|\Sigma_\eta^{-1}\|_2$.

**Theorem 3.6.** Add the range conditions $\mathrm{Range}(G^{T})=\mathrm{Range}(\widehat G^{T})=\mathbb R^{N_m}$ together with $\Sigma_\omega\succ0$ and $\Sigma_\eta\succ0$. Then the surrogate-driven unscented Kalman inversion fixed point $(\widehat r_\infty,\widehat C^{-1}_\infty)$ converges to the full-model fixed point $(r_\infty,C^{-1}_\infty)$, with

$$
\bigl\|\widehat C^{-1}_\infty-C^{-1}_\infty\bigr\|_2\le\frac{2\epsilon HH_\eta}{1-\beta},
\qquad
\bigl\|\widehat r_\infty-r_\infty\bigr\|_2\le\frac{K_1H_\eta H_y}{C_1}
\Bigl(1+\frac{2(1+\alpha\beta)K_2H_\eta H^{2}}{(1-\beta)C_2}\Bigr)\epsilon,
$$

all constants positive and bounded. Both are $O(\epsilon)$, that is, linear in the surrogate's operator-norm error.

> [!warning] The theorem does not cover the network used in the experiments
> Paper 79's Remark 1 says so directly: to satisfy Theorem 3.6 one may drop the nonlinear activations in the branch net while keeping them in the trunk net, making $\widehat G_\theta$ a linear operator. The theorem therefore governs a restricted DeepONet, while every numerical experiment uses the nonlinear version. The paper points this gap out itself.

### What is proved and what is only argued

| Statement                                                                     | Source                       | Scope and qualification                                        |
| ----------------------------------------------------------------------------- | ---------------------------- | -------------------------------------------------------------- |
| Uniformly close potentials $\Rightarrow$ small Hellinger distance             | Stuart, Acta Numer. 2010     | demands global accuracy, the very cost this family avoids      |
| Local accuracy + small bad-set posterior mass $\Rightarrow$ small KL          | 37 Thm 2 = 79 Thm 2.1        | compares two posteriors only, says nothing about the output    |
| Hellinger version, without the square                                         | 37 Cor. 3, attributed to [8] | as above                                                        |
| Surrogate operator error $\epsilon$ $\Rightarrow$ inversion fixed point $O(\epsilon)$ | 79 Thm 3.6            | linear operator network only, unscented Kalman inversion only  |
| Local polynomial regression error $O(R^{2})/O(R^{3})$ in a ball of radius $R$ | Conrad et al. 2016, background | needs a controlled poisedness constant, which a network lacks   |
| Asymptotic exactness under unbounded refinement                               | Conrad et al. 2016, background | none of the ten papers takes this route                        |
| Convergence of the adaptive loop, or bad-set measure going to zero            | **none**                     | argued as a mechanism in every paper                            |

The last two rows are the key to reading this family. The Conrad–Marzouk–Pillai–Smith route buys asymptotic exactness with a Metropolis correction: the samples eventually follow the true posterior. This topic's route gives up that guarantee. **The subchains are run against the surrogate, so the samples follow the surrogate-induced posterior**, and correctness rests entirely on the KL bound above showing the two posteriors are close. The two guarantees differ in strength and in price — the first needs several leave-one-out refits per step, the second one or two true solves per outer iteration.

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

requires a high-dimensional integral and is not implementable. Every paper therefore retreats to a single-point quantity. Here is the whole family, with the baseline it departs from on the first row:

| No.               | Indicator                                                | Evaluated at                                    | Trigger                              |
| ----------------- | -------------------------------------------------------- | ----------------------------------------------- | ------------------------------------ |
| Conrad et al.     | leave-one-out sensitivity $\epsilon^{\pm}$ of the acceptance probability | both $\theta^{-}$ and $\theta^{+}$ | $>\gamma_t$, refine at the larger    |
| 37                | **absolute** $\ell^\infty$ error                         | a state accepted by a true-model Metropolis step | $>\epsilon$                          |
| 34                | relative $\ell^\infty$ error                             | the ensemble mean                               | $>\mathrm{tol}$                      |
| 49                | relative $\ell^\infty$ error                             | an accepted state                               | $>\mathrm{tol}$                      |
| 56                | relative $\ell^2$ error                                  | the particle mean                               | $>\mathrm{tol}$                      |
| 79                | $e_D(t)=\Phi(r_t;y)$                                     | the anchor: best true-model fit among posterior samples | relative change $>\epsilon$  |
| 106               | $\Phi$ at the current prior mean                         | the evolving prior mean                         | relative change $<\epsilon$ means **stop** |

Two design principles run through the table. The first is **goal orientation**: the indicator should measure error in the quantity that actually drives the sampler, not in the forward model as such. Conrad and coauthors put it most sharply — the scale of an error in forward-model outputs or in the log-likelihood cannot be known in advance, whereas an error in the acceptance probability reads directly as an additive error in a probability, which lets a user set the threshold with problem-independent intuition. Papers 79 and 106 arrive at the same place from the other side: they use the data misfit $\Phi$ itself, the objective the inversion is trying to reduce. Papers 34, 37, 49 and 56 compromise, measuring raw model error but **only at one point believed to sit in the posterior bulk**. The second principle is **locality through the design point**: the expensive high-fidelity evaluation is always spent at a single point — an accepted chain state, the ensemble mean, the particle mean, the best-fitting posterior sample, the current prior mean. That is the cost-control device itself: the Conrad cross-validation indicator needs $N$ leave-one-out refits per step, while these need one forward solve per outer iteration.

### New training points must not degenerate

If refinement draws points randomly inside a small ball, the points cluster and the least-squares or training problem degenerates. Three remedies appear across the papers: shrink the ball radius each round as $R\leftarrow\rho R$ (papers 37, 49, 56); impose a separation constraint $\|x'-x\|_2\ge R$ on candidates (paper 56); and run a greedy diversity selection in the surrogate's output space with a coefficient pulling selections back toward the anchor (paper 79). Paper 106 goes the other way and discards the historical data entirely, fine-tuning only on the newly generated local set.

### Surrogate error must not enter the posterior directly

If the whole chain uses the surrogate's acceptance probability, the result is the surrogate's posterior. These algorithms therefore insert the true model at one decisive place: papers 37 and 49 compute the trigger point using the **true-model acceptance probability**, so the point being tested sits closer to the posterior bulk; paper 79 uses the true model to select an anchor from a candidate pool. The cost is a few true solves per round, and it is exactly this step that makes refinement move toward the true posterior instead of the surrogate's.

Read that carefully against the last section: what the true model corrects here is **where refinement happens**, not the samples themselves. The samples still come from the surrogate-induced posterior; it is the KL bound that licenses using them.

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]. Per-paper references appear at the end of each close-reading page.
