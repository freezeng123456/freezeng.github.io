---
title: Multi-Fidelity Surrogates Refined Online
description: Papers 34, 37, 49 and 79 - one loop rebuilt four times, with full derivations, theorem hypotheses and experiments
lang: en
translation: computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates
tags:
  - paper-notes
  - bayesian-inverse-problems
  - surrogate-modelling
---

> [!note] Coverage of this page
> Papers **37** (_J. Comput. Phys._ 381, 2019), **34** (_Int. J. Uncertain. Quantif._ 9(3), 2019), **49** (_Commun. Comput. Phys._ 28, 2020) and **79** (_SIAM/ASA J. Uncertain. Quantif._ 12(4), 2024). All four share one closed-loop skeleton and replace one component at a time.

## Verification tier

Before writing, the level at which each paper could be checked was fixed, and the depth of what follows is allocated strictly by that table: **only papers verified in full text get complete derivations and experimental configurations**, and everything else is limited to what a source actually supports.

| No. | Source checked                          | Tier          | Remaining gap                                                                                       |
| --- | --------------------------------------- | ------------- | --------------------------------------------------------------------------------------------------- |
| 37  | arXiv:1807.00618 (ar5iv full text)      | **Full text** | no transcribable error values are recorded, so the experiments section reports setup and abstract-level conclusions only |
| 34  | arXiv:1809.08931 (full text) + abstract | **Full text** | in the cost accounting $J_1+Q_1+J_2Q_2$ the definitions of $J_1,J_2$ were lost and are unverified   |
| 49  | arXiv:1911.08926 (full text)            | **Full text** | equation (10) appears to have the acceptance-ratio numerator and denominator swapped, see the warning below; no CPU timings for the adaptive version |
| 79  | arXiv:2310.17844 v3 (full text)         | **Full text** | no substantive gap; Theorem 3.6 applies only to a linearised operator network, see the warning below |

All four are at the full-text tier, so this page gives all four a derivation chain, full theorem hypotheses and a reproducible experimental configuration. What is missing is **the measured values themselves**: only paper 34 left a complete set of error and cost figures in the verifiable material, while papers 37, 49 and 79 left configurations and qualitative conclusions. Each experiments section says so in place.

## The skeleton the four share

Read together, the four answer one set of three questions and simply give different answers each time:

1. **What surrogate** — a low-order polynomial chaos expansion with a high-fidelity correction (37, 34), a composite network (49), or an operator network (79).
2. **When to admit the surrogate is not good enough** — by comparing surrogate against true model at one carefully chosen point (37, 34, 49), or by asking whether the fit itself has stopped improving (79).
3. **Where to add points once triggered** — randomly inside a local ball (37, 49), inside a ball (34), or greedily in the surrogate's output space (79).

The loop itself is fixed: run some sampling against the current surrogate, evaluate an indicator at one point, and if the indicator exceeds its threshold spend a few true solves retraining before returning to sampling. **Cost control lives entirely in the second step**: the indicator is allowed one (or a few) true solves per round, so it has to be a single-point quantity, and for a single-point quantity to mean anything that point has to sit near the posterior bulk. Everything that distinguishes the four papers reduces to how they find that point.

## 37: the prototype of the loop

### The idea

A prior-based polynomial chaos surrogate is built like this: draw parameters from the prior, solve the true problem at each, and fit the input-to-output map by a polynomial. The surrogate is accurate where those prior samples fell — and the posterior is not there. Once the data are informative the posterior contracts onto a thin subset of the prior support, and the more informative the data, the thinner that subset, until it can sit in the prior's tail. The surrogate's accuracy and the posterior's mass are then **very nearly orthogonal**: the computational effort went into a region the posterior will never visit.

Two standard escapes exist and the paper rejects both. Raise the polynomial order until the surrogate is globally accurate — the total-degree space has dimension $\binom{N+n_z}{n_z}$, which grows combinatorially in the parameter dimension and is out of reach already at moderate dimension. Or pick a prior closer to the posterior — but the posterior's location is not known in advance, so that is circular.

The paper takes a third route: **keep the order low and let the surrogate follow the chain**. The chain walks into the posterior bulk on its own, so spend a true solve where the chain went and repair the surrogate there. The surrogate is then never globally accurate; it is accurate only where the chain cares — and Theorem 2 below says precisely that this is enough.

One problem remains. The chain is driven by the surrogate, so where it walks is the bulk of the **surrogate's** posterior, not necessarily of the true one. Trusting the surrogate's acceptance probability completely would leave the chain circling inside the surrogate's own high-probability region, where the indicator will of course pass, refinement will never trigger, and the error will have certified itself as absent. The paper's patch is to insert, once per outer iteration, a Metropolis decision that uses the **true-model acceptance probability**, and to use the point it accepts as the test point. This costs two true solves per round (one at $z^-$ and one at $z^+$) and buys a test point chosen by the true posterior rather than by the surrogate's. **That is what "a surrogate inside MCMC needs the acceptance step to correct for it" concretely means in this family** — and note that what it corrects is where refinement happens, not the samples themselves, on which see the honest note below.

### Setting

Parameter $z\in\Gamma\subset\mathbb R^{n_z}$ with prior density $\pi(z)$, data $d\in\mathbb R^{n_d}$, forward operator $u^{H}:\Gamma\to\mathbb R^{n_d}$ (high fidelity, i.e. the true model) and surrogate $u^{L}$ (low fidelity). Observational noise is i.i.d. Gaussian. The true posterior is $\pi^{d}$ and the surrogate-induced posterior is $\widetilde\pi^{d}_{N}$.

The low-fidelity model is a prior-based polynomial chaos expansion in a basis $\{\Phi_m\}_{m=1}^{M}$ orthogonal with respect to the prior measure, with $M=\dim\Lambda_N$ the size of the total-degree index set $\Lambda_N$.

### Derivation

**Step one: coefficients by weighted discrete least squares.** At $Q$ sample points $\{z_i\}$,

$$
\{c_m\}_{m=1}^{M}=\arg\min_{c}\sum_{i=1}^{Q}
\Bigl[G(z_i)-\sum_{m=1}^{M}c_m\Phi_m(z_i)\Bigr]^{2},
\qquad
w_i=\frac{M}{\sum_{m=1}^{M}\Phi_m^{2}(z_i)} .
$$

In matrix form, with $\Psi\in\mathbb R^{Q\times M}$, $\Psi_{im}=\Phi_m(z_i)$ and $W=\mathrm{diag}(w_i)$, the problem solved is $c^{\#}=\arg\min_c\|\sqrt W\Psi c-\sqrt W b\|_2^2$.

The weight is not arbitrary. Writing the reproducing-kernel diagonal $K(z)=\sum_{m=1}^{M}\Phi_m^{2}(z)$ gives $w_i=M/K(z_i)$, and then row $i$ of $\sqrt W\Psi$ has squared norm

$$
\sum_{m=1}^{M}w_i\Phi_m^{2}(z_i)=\frac{M}{K(z_i)}\cdot K(z_i)=M,
$$

**independently of where $z_i$ landed**. The weighting normalises every row of the design matrix to the same length $\sqrt M$, so the trace of $\frac1Q\Psi^{T}W\Psi$ is identically $M$; the only remaining job of the sampling density is to decide how that fixed trace is distributed across the eigenvalues.

The weight $w_i$ is the discrete Christoffel function, the same object studied on the [[en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|optimal sampling and preconditioning]] page: the stability factor there is $\|K\|_\infty/N$, which is the supremum of the reciprocal of this weight, and induced sampling $q\propto\sum_n v_n^{2}$ is exactly the choice that drives that supremum to its lower bound. The sampling density follows degree-asymptotic designs: tensor-product Chebyshev for a uniform measure, and for a Gaussian measure a density supported on the ball of radius $\sqrt{2N}$ (Narcowich–Ward-style).

**Step two: the correction is expanded only on low-order indices.** Define the discrepancy between fidelities,

$$
C(z)=u^{H}(z)-u^{L}(z)\ \approx\ \sum_{\mathbf m\in\Lambda_{N_C}}u^{C}_{\mathbf m}\Phi_{\mathbf m}(z),
\qquad N_C\ll N .
$$

**Step three: merging is only a regrouping.** Add $C$ to $u^{L}=\sum_{\alpha\in\Lambda^{d}_{N}}u^{L}_{\alpha}\Psi_{\alpha}$; since $\Lambda^{d}_{N_C}\subset\Lambda^{d}_{N}$, collecting index by index gives

$$
f_{M}(\theta)=\sum_{\alpha\in\Lambda^{d}_{N_C}}\bigl(u^{L}_{\alpha}+u^{C}_{\alpha}\bigr)\Psi_{\alpha}
+\sum_{\alpha\in\Lambda^{d}_{N}\setminus\Lambda^{d}_{N_C}}u^{L}_{\alpha}\Psi_{\alpha}.
$$

Nothing is approximated here; the sum is rearranged. Its value is that it makes the structure visible: **only the low-order coefficients are corrected by high-fidelity data, and the high-order coefficients keep their prior-fitted values**. That is deliberate rather than an omission — each round affords only $Q=2\binom{N_C+n_z}{n_z}$ new true solves, twice the dimension of the correction space, which is the usual oversampling factor for weighted least squares, and spreading that many samples over the high-order indices would overfit.

**Step four: the trigger.** Refinement is decided at an accepted point using an **absolute** $\ell^\infty$ error:

$$
\mathrm{err}(y)=\bigl\|u^{H}(y)-u^{L}(y)\bigr\|_{\infty}.
$$

Here $y$ is not an arbitrary chain state but a point accepted with the high-fidelity acceptance probability, for the reason given above.

**Step five: the local refinement design.** If $\mathrm{err}(y)>\epsilon$, draw $Q$ random points in the $\ell^\infty$ ball $B(y,R)=\{x:\|x-y\|_{\infty}\le R\}$ and refit. The radius shrinks each round by an input constant $\rho$, because early accepted points may still be far from the posterior bulk: a large ball is needed at first to cover an uncertain location, and a small ball later to buy local accuracy.

**The rhythm of the algorithm.** Each outer round first runs $m-1$ steps of standard Metropolis-Hastings against the current low-fidelity model (the subchain exists to decorrelate its last state from the starting state), then proposes $z^{*}\sim q(\cdot\mid z_{m-1})$, uses the true-model acceptance probability between $z^{-}=z_{m-1}$ and $z^{+}=z^{*}$ to produce $y$, evaluates the indicator, refines if required, and repeats at most $I_{\max}$ times. **The final posterior sample is the union of all subchain samples.**

> [!warning] Whose posterior this algorithm samples
> That last sentence deserves a pause. The subchains run against the current low-fidelity model, so the samples in the union are distributed according to the **surrogate-induced posterior** $\widetilde\pi^{d}_{N}$, not the true posterior $\pi^{d}$. The true model appears in exactly two places: choosing the test point, and the $Q$ evaluations used to refine. The correctness argument is therefore not that a Metropolis correction pulls the samples back onto the true posterior — that is the asymptotic-exactness route of Conrad–Marzouk–Pillai–Smith — but Theorem 2 below: if $\widetilde\pi^{d}_{N}$ is close enough to $\pi^{d}$ in KL, using the former in place of the latter is good enough. The two routes give guarantees of different strength and should not be conflated.

### Theorems

**Assumption 1.** The forward operator is uniformly bounded: $\sup_{z\in\Gamma}\|u^{H}(z)\|=:C_H<\infty$.

Define the $\epsilon$-feasible set and its posterior measure,

$$
\Gamma_{N}(\epsilon)=\bigl\{y\in\Gamma:\ \|u^{H}(y)-u^{L}(y)\|_{\infty}\le\epsilon\bigr\},
\qquad
\mu\bigl(\Gamma_{N}(\epsilon)\bigr)=\int_{\Gamma_{N}(\epsilon)}\pi^{d}(z)\,dz,
$$

with complement $\Gamma^{\perp}_{N}(\epsilon)=\Gamma\setminus\Gamma_{N}(\epsilon)$.

**Theorem 2 (KL bound).** Assume $u^{H}$ and $u^{L}$ satisfy Assumption 1 uniformly in $N$ and that the observational error is i.i.d. Gaussian. Then for a given $\epsilon>0$ there exist constants $K_1,K_2>0$ with

$$
D_{\mathrm{KL}}\bigl(\widetilde\pi^{d}_{N}\,\|\,\pi^{d}\bigr)
\le\Bigl(K_1\epsilon+K_2\,\mu\bigl(\Gamma^{\perp}_{N}(\epsilon)\bigr)\Bigr)^{2}.
$$

**Corollary 3 (Hellinger bound, attributed by the paper to its reference [8]).** Under the same hypotheses,

$$
D_{\mathrm{Hell}}\bigl(\widetilde\pi^{d}_{N}\,\|\,\pi^{d}\bigr)
\le K_1\epsilon+K_2\,\mu\bigl(\Gamma^{\perp}_{N}(\epsilon)\bigr),
$$

that is, without the square.

**Corollary (the paper's equation 21, which states the algorithmic goal directly).** If sampling is good enough that $\mu(\Gamma^{\perp}_{N}(\epsilon))\le\epsilon$, the KL distance is characterised entirely by $\epsilon^{2}$.

The shape of the bound records where it comes from: the two terms correspond to splitting the KL integral over $\Gamma_{N}(\epsilon)$ and its complement. On the feasible set the potentials differ by at most $\epsilon$, which gives $K_1\epsilon$; on the complement the surrogate may be arbitrarily bad, so the integrand can only be controlled by the uniform bound of Assumption 1, weighted by the posterior mass of that region, which gives $K_2\mu(\Gamma^{\perp}_{N}(\epsilon))$. That is why an apparently redundant uniform-boundedness assumption is needed at all. (This reading of the bound's structure is this page's; the paper states only the conclusion.)

> [!warning] What is proved and what is only argued
> The paper deliberately presents the link between algorithm and theorem as a **mechanism** rather than a theorem: whenever a candidate falls in $\Gamma^{\perp}_{N}(\epsilon)$, the algorithm refines near it, so $\mu(\Gamma^{\perp}_{N}(\epsilon))$ decays asymptotically with refinement. **That is an argument, not a rate, and it says nothing about when or how fast $\mu(\Gamma^{\perp}_{N}(\epsilon))$ reaches zero.** Every paper in this topic is in the same position on this point.

### Numerical experiments

Two nonlinear PDE inverse problems:

| Ex. | Problem                                    | Parameter dimension | Design intent                                                                     |
| --- | ------------------------------------------ | ------------------- | --------------------------------------------------------------------------------- |
| 1   | two-dimensional heat source inversion      | $n_z=2$             | small enough to afford a genuinely accurate high-order prior surrogate, so what is bought here is accuracy |
| 2   | diffusion coefficient of an elliptic PDE   | $n_z=9$             | a globally accurate prior surrogate is expensive, so both accuracy and cost are bought |

The division of labour is deliberate. Example 1 rules out the reading that the adaptive method wins only because its opponent is weak, since at $n_z=2$ the opponent can be made globally accurate; example 2 is the method's actual target. The abstract's order-of-magnitude claim is several orders of efficiency gain over MCMC on the true model alone.

> [!note] The numerical evidence for this paper stops here
> The verifiable material contains the configuration above and the abstract's order-of-magnitude statement, and **no transcribable error values, chain lengths, CPU timings or true-solve counts**. This page therefore prints no results table for paper 37. For this method with numbers attached, see paper 34 in the next section: the same correction inside a different sampler, with errors and solve counts reported in full.

**What the experiments establish and where they fall short.** They establish that a low-order surrogate plus local correction recovers accuracy on two nonlinear PDE inverse problems. Three things are missing. First, the two examples have parameter dimensions 2 and 9, while the paper's own stated motivation is the combinatorial blow-up in high dimension, which is therefore never tested directly (paper 49 exists precisely for that). Second, Theorem 2 carries undetermined constants $K_1,K_2$, so the experiments neither do nor can calibrate $\epsilon$ from it, and the tolerance is supplied as a tuning parameter. Third, the mechanism claim that $\mu(\Gamma^{\perp}_{N}(\epsilon))$ decays with refinement is never measured — what is measured is whether the posterior marginals match, not whether the measure of the bad set is shrinking.

### Relation to the others

This is where the family starts. Paper 34 transplants it into ensemble Kalman inversion, paper 49 into neural networks (while criticising this paper's two weaknesses explicitly: polynomials handle limited regularity poorly and suffer the curse of dimensionality), and papers [[en/computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|55 and 56]] carry the "train where the posterior lives" principle into randomize-then-optimize and Stein variational gradient descent. Theorem 2 here is Theorem 2.1 of paper 79 (both cite the same source), and it is the $\epsilon$-feasible-set bound on the [[en/computational-mathematics/paper-notes/bayesian-inference/index|topic index]].

## 34: the same correction inside ensemble Kalman inversion

### The idea

Ensemble Kalman inversion needs no derivative of the forward operator, which is decisive under a PDE constraint: the adjoint may not be implemented, or the solver may simply not be differentiable. The price is that it substitutes the empirical covariance of a particle ensemble for derivative information, and an empirical covariance is only trustworthy if the ensemble is not small — $N_e$ members over $J$ iterations is roughly $N_eJ$ forward solves.

Replacing the forward model by a prior-based polynomial chaos surrogate drives the online cost to nearly nothing. But ensemble Kalman inversion fails differently from MCMC, and this paper's experimental design goes straight at that difference: **the ensemble actively moves in the direction the data point, and if the data point outside the prior, the ensemble walks out of the surrogate's reliable region while the surrogate raises no alarm — it just keeps returning smooth, wrong outputs**. Drawing the truth from $\log\theta_i\sim U(-4,4)$ against a prior of $\log\theta_i\sim\mathcal N(0,1)$ is that failure mode turned up to maximum.

Choosing the test point is more natural here than in MCMC. MCMC needs an extra true-model acceptance decision to produce a trustworthy test point; ensemble Kalman inversion has a candidate ready-made — the **ensemble mean**, which is the ensemble's current best guess at where the truth is. That removes the extra acceptance decision of paper 37.

### Setting

Forward model and observations (the paper's equation 2.1):

$$
y=f(\theta)+\xi,\qquad \xi\sim\mathcal N(0,\Gamma),
\qquad \theta\in\mathbb R^{d},\ y\in\mathbb R^{m}.
$$

The ensemble has $N_e$ members; member $j$ at iteration $n$ is $\theta^{(j)}_n$ with output $\omega^{(j)}_n=f(\theta^{(j)}_n)$, and the data are perturbed as $y^{(j)}=y+\xi^{(j)}$, $\xi^{(j)}\sim\mathcal N(0,\Gamma)$.

### Derivation

**Step one: where the update comes from.** The ensemble Kalman inversion of Iglesias–Law–Stuart works in the **extended state space** $Z=X\times Y$ with $z=(u,p)$, $p=\mathcal G(u)$ and projection $H=[0,I]$. The prediction step advances each particle under the artificial dynamics, $\hat z^{(j)}_{n+1}=(u^{(j)}_n,\mathcal G(u^{(j)}_n))$, and the analysis step applies the Kalman gain

$$
K_n=C_nH^{T}\bigl(HC_nH^{T}+\Gamma\bigr)^{-1},
\qquad
z^{(j)}_{n+1}=\hat z^{(j)}_{n+1}+K_{n+1}\bigl(y^{(j)}_{n+1}-H\hat z^{(j)}_{n+1}\bigr).
$$

Splitting $C_n$ into blocks according to $z=(u,p)$, substituting $H=[0,I]$ and projecting out the parameter component with $H^{\perp}=[I,0]$ turns $HC_nH^{T}$ into $C^{pp}_n$ and the parameter block of $C_nH^{T}$ into $C^{up}_n$, so that

$$
u^{(j)}_{n+1}=u^{(j)}_{n}+C^{up}_{n+1}\bigl(C^{pp}_{n+1}+\Gamma\bigr)^{-1}
\bigl(y^{(j)}_{n+1}-\mathcal G(u^{(j)}_{n})\bigr).
$$

(This reduction is the standard one and is given here as such; Iglesias and coauthors write the pre-block form.) Paper 34 uses exactly this, with one regularisation parameter $\alpha_n$ multiplying $\Gamma$:

$$
\theta^{(j)}_{n+1}=\theta^{(j)}_{n}
+C^{\theta\omega}_{n}\bigl(C^{\omega\omega}_{n}+\alpha_{n}\Gamma\bigr)^{-1}\bigl(y^{(j)}-\omega^{(j)}_{n}\bigr),
\qquad j=1,\dots,N_e,
$$

with empirical covariances

$$
C^{\theta\omega}_{n}=\frac{1}{N_e-1}\sum_{j=1}^{N_e}\bigl(\theta^{(j)}_{n}-\bar\theta_{n}\bigr)\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)^{\!T},
\qquad
C^{\omega\omega}_{n}=\frac{1}{N_e-1}\sum_{j=1}^{N_e}\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)^{\!T}.
$$

**Step two: why this is called iterative regularisation.** Three things act at once. The first is structural: the ensemble never leaves the linear span of the initial ensemble, $\mathcal A=\mathrm{span}\{\psi^{(j)}\}_{j=1}^{J}$, so a misfit minimisation that is ill-posed on $X$ is in fact carried out on the compact set $\mathcal A$. The second is a comparison: in the linear case $\mathcal G(u)=Gu$, one ensemble update converges as $J\to\infty$ to the Tikhonov–Phillips solution

$$
u_{TP}=\bar u+CG^{*}\bigl(GCG^{*}+\Gamma\bigr)^{-1}(y-G\bar u),
$$

which is **exactly the linear-Gaussian posterior mean** (Stuart's (3.4) with finite-dimensional data). Ensemble Kalman inversion is therefore not an approximation in the linear case but exact, at the sole cost of not needing the derivative of $\mathcal G$. The third is algorithmic: the iteration must be stopped by the discrepancy principle.

**Step three: regularisation parameter and stopping.** $\alpha_n$ is the value attached to the first integer $N$ satisfying

$$
\alpha^{N}_{n}\bigl\|\Gamma^{1/2}\bigl(C^{\omega\omega}_{n}+\alpha^{N}_{n}\Gamma\bigr)^{-1}\bigl(y^{(j)}-\bar\omega_{n}\bigr)\bigr\|
\ \ge\ \rho\,\bigl\|\Gamma^{-1/2}\bigl(y^{(j)}-\bar\omega_{n}\bigr)\bigr\| ,
$$

and the iteration stops by the discrepancy principle

$$
\bigl\|\Gamma^{-1/2}(y-\bar\omega_{n})\bigr\|\le\tau\eta,
\qquad \rho<1,\ \tau\ge1/\rho .
$$

Both come from Iglesias' regularising ensemble Kalman smoother. **Iglesias and coauthors state explicitly that a complete convergence and regularisation analysis is beyond their scope and that the discrepancy principle is supported numerically, not proved.** The stopping rule therefore remains an empirical device throughout this family. The same source records a counterintuitive empirical fact: switching off the data perturbation ($\eta^{(j)}_{n+1}=0$) made results **worse**, and the authors conjecture the noise helps the algorithm explore $\mathcal A$.

**Step four: a relative indicator, tested at the ensemble mean.**

$$
\mathrm{err}=\frac{\bigl\|f(\bar\theta_{n+1})-f_{M}(\bar\theta_{n+1})\bigr\|_{\infty}}{\bigl\|f(\bar\theta_{n+1})\bigr\|_{\infty}},
\qquad
\bar\theta_{n+1}=\frac{1}{N_e}\sum_{j}\theta^{(j)}_{n+1}.
$$

If $\mathrm{err}\le\mathrm{tol}$ the current surrogate is kept; otherwise $f_M$ is rebuilt by the multi-fidelity least-squares routine. The sample budget is explicit: initialisation uses $Q_1=2\binom{N+d}{d}$ prior samples and each refinement uses $Q_2=2\binom{N_C+d}{d}$ new points.

### Theorems

**This paper proves no convergence theorem for the adaptive scheme, and all its claims are computational.** A total-cost accounting of the form $J_1+Q_1+J_2Q_2$ appears in the text, but the definitions of $J_1$ and $J_2$ were lost in the verifiable material and are not repeated here. The linear-case Tikhonov–Phillips limit and the discrepancy principle quoted above are results of Iglesias and coauthors, not of this paper.

### Numerical experiments

The forward problem is a two-dimensional **time-fractional** inverse diffusion problem (Caputo derivative of order $0<\alpha<1$).

| Item                    | Setting                                              |
| ----------------------- | ---------------------------------------------------- |
| time discretisation     | finite differences, $\Delta t=0.01$                  |
| space discretisation    | spectral, polynomial degree $P=6$                    |
| data generation         | $P=10$, to avoid an inverse crime                    |
| example 1 unknown       | nine radial-basis weights in the permeability ($d=9$) |
| example 2 unknown       | a Karhunen-Loève parameterised log-permeability field |

**Example 1: the truth is deliberately placed outside the prior.** The prior is $\log\theta_i\sim\mathcal N(0,1)$ while the truth is drawn from $\log\theta_i\sim U(-4,4)$; $N_e=100$ and the noise level is $\sigma=10^{-3}$.

| Method                                  | PC order | Relative error | True-model evaluations | CPU time     |
| --------------------------------------- | -------- | -------------- | ---------------------- | ------------ |
| conventional ensemble Kalman inversion  | —        | 0.0461         | 2000 (online)          | about 56.71 s |
| fixed prior polynomial chaos surrogate  | $N=4$    | 0.7921         | 1430 (offline)         | about 0.82 s  |
| fixed prior polynomial chaos surrogate  | $N=6$    | 0.2892         | 10010 (offline)        | not reported  |
| adaptive multi-fidelity (tol $10^{-2}$) | $N=2$    | not reported   | 250                    | not reported  |
| adaptive multi-fidelity (tol $10^{-3}$) | $N=2$    | not reported   | 575                    | not reported  |

Several points are worth pulling out.

- The offline evaluation counts **match the budget formula**: with $d=9$, $Q_1=2\binom{N+d}{d}$ gives $2\binom{13}{9}=1430$ at $N=4$ and $2\binom{15}{9}=10010$ at $N=6$, exactly the printed figures. That confirms both that the parameter dimension really is 9 and that the budget formula is implemented literally.
- The 2000 evaluations of the conventional method correspond to 20 iterations at $N_e=100$, consistent with the $N_eJ$ cost model.
- The relative error of the adaptive version has **no numerical value** in the verifiable material, only the qualitative statement that the local correction essentially recovers full-model accuracy. Those two cells are therefore marked as not reported rather than filled with $0.0461$.
- The comparison on solve counts is unambiguous: 250 against 2000 is an eighth, 575 against 2000 is under a third, and both are far below the 10010 needed to push the fixed surrogate to $N=6$.

**Example 2: a high-dimensional random field.** Retaining 95% of the prior energy gives $d=22$ Karhunen-Loève modes, with $N_e=300$ and noise $\mathcal N(0,0.01^2)$.

| Method                                 | PC order | Relative error |
| -------------------------------------- | -------- | -------------- |
| fixed prior polynomial chaos surrogate | $N=2$    | 0.3430         |
| fixed prior polynomial chaos surrogate | $N=3$    | 0.2146         |
| adaptive multi-fidelity                | $N=1$    | 0.0889         |

**What the experiments establish and where they fall short.** They establish a rather sharp conclusion: **when the truth is out of prior, raising the order of a global surrogate is an inefficient repair.** Going from $N=4$ to $N=6$ raises the offline evaluation count from 1430 to 10010, a factor of seven, and only moves the error from 0.7921 to 0.2892 — still six times worse than the conventional method — while the adaptive method needs 250 solves on an expansion of order two. Example 2 is more extreme still: the adaptive method reaches 0.0889 on a **first-order** expansion, an order of magnitude better than the fixed second- and third-order surrogates.

Four things are missing. No CPU time is reported for the adaptive version, so "cheap" rests on true-solve counts rather than measured wall-clock time, and the cost of the multi-fidelity refits themselves never enters the comparison. No numerical relative error is reported for the adaptive version either, so "essentially recovers full-model accuracy" cannot be checked quantitatively. Both examples live in the same PDE family (time-fractional diffusion), so transfer across families is untested. And the paper has no theorem, while the out-of-prior setting is precisely where a Theorem 2-type bound is hardest to apply: $\mu(\Gamma^{\perp}_{N}(\epsilon))$ starts large when the truth is out of prior, and how it decays remains an argument.

### Relation to the others

This is the ensemble Kalman version of paper 37: the same multi-fidelity correction, the same family of $\ell^\infty$ triggers, with the driver changed from an MCMC chain to an ensemble iteration, the indicator from absolute to relative, and the test point from an accepted state to the ensemble mean. Paper 49 swaps in a composite network while keeping the trigger; paper 79 swaps the surrogate again (DeepONet), the sampler again (unscented Kalman inversion), and the indicator for a data misfit. Paper [[en/computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|82]] moves "a cheap model plus a learned correction inside a Kalman method" into time-dependent data assimilation.

## 49: replacing polynomials with a composite network

### The idea

There is nothing wrong with the loop of paper 37; the problem is its surrogate. Polynomial chaos has two weaknesses that bite here: poor approximation of maps with limited regularity, and a basis count that explodes combinatorially with parameter dimension. Both block the move to high dimension directly.

Switching to a network raises a new and very concrete problem: **only about ten true solves are affordable per round, and training a network on ten samples will overfit.** The paper's answer is not to shrink the network until it is useless but to change what is being learned — not $z\mapsto f^{H}(z)$, but where the existing surrogate is wrong. And not, as in paper 37, as an additive correction term: the already-trained low-fidelity network is fed into a new shallow network **as an input variable**.

The distinction matters. An additive correction $f^{H}\approx f^{L}+f_{\mathrm{CORR}}$ assumes the discrepancy between fidelities is itself an easy function to fit. The composite form $f^{H}\approx\mathcal{NN}(z,\mathcal{NN}^{L}(z))$ makes no such assumption; it lets the high-fidelity output depend on the low-fidelity output in any way at all, including multiplicative rescaling or piecewise switching. Because the two models are highly correlated, that map is itself simple — **which is why a shallow network suffices, and a shallow network sufficing is why ten samples do not overfit**. The capacity constraint is a link in that chain, not a regularisation applied afterwards.

### Setting

Parameter $z$, data $d$, likelihood $\mathcal L(d,\cdot)$, prior $\pi$. The low-fidelity model $f^{L}=\mathcal{NN}^{L}$ is a deep network trained on prior samples (4 hidden layers of 40 neurons on 50 training points in the experiments); the high-fidelity model $f^{H}$ is the true solve.

### Derivation

**Step one: the composite multi-fidelity network.** Posit an unknown nonlinear relation between the two fidelities,

$$
f^{H}(z)=\mathcal F\bigl(z,f^{L}(z)\bigr)=\mathcal F\bigl(z,\mathcal{NN}^{L}(z)\bigr),
$$

and approximate that relation with a network:

$$
f^{H}(z)\approx\mathcal{NN}^{H}(z;\theta):=\mathcal{NN}\bigl(z,\mathcal{NN}^{L}(z);\theta\bigr).
$$

The paper draws the contrast with paper 37 itself: the polynomial chaos correction $f^{H}\approx f^{L}_{\mathrm{PCE}}+f_{\mathrm{CORR}}$ is a **linear** superposition, whereas the display above learns a **nonlinear** correlation between fidelities.

**Step two: building the training set.** Choose $Q$ points $\{z_k\}_{k=1}^{Q}$, evaluate $\mathcal{NN}^{L}$ (cheap) and $f^{H}$ (expensive) at each, and form

$$
\mathcal D=\Bigl\{\bigl((z_k,\mathcal{NN}^{L}(z_k)),\ f^{H}(z_k)\bigr)\Bigr\}_{k=1}^{Q},
$$

then train $\mathcal{NN}^{H}(z;\theta)$ on $\mathcal D$. The input is a **pair**: the low-fidelity prediction enters as an extra input channel, and that is how the composition is realised.

**Step three: the capacity chain.** $Q$ must be small (true solves are expensive; $Q=10$ in the experiments) $\Rightarrow$ the correction network must be shallow (one hidden layer of 50 neurons in the experiments, described in the paper as at most two hidden layers) $\Rightarrow$ this is only viable when the two fidelities are highly correlated and the map to be learned is itself simple. **The paper describes the shallow network explicitly as an anti-overfitting constraint rather than a performance optimisation.**

**Step four: indicator and refinement.**

$$
\mathrm{err}(\tilde z)=\frac{\bigl\|f^{H}(\tilde z)-f^{L}(\tilde z)\bigr\|_{\infty}}{\bigl\|f^{H}(\tilde z)\bigr\|_{\infty}} .
$$

The absolute error of paper 37 becomes relative. If $\mathrm{err}(\tilde z)>\mathrm{tol}$, draw $Q$ uniform points in $B(\tilde z,R)=\{z:\|z-\tilde z\|_\infty\le R\}$, retrain $\mathcal{NN}^{H}$, and set $f^{L}\leftarrow\mathcal{NN}^{H}$.

**Step five: the composition nests round by round.** That last assignment has a consequence that is easy to miss: the next round's "low-fidelity model" is the previous composite network, so the round-$k$ network takes the round-$(k-1)$ network as input, which in turn takes round $k-2$. Composition depth grows with the outer loop even though every correction layer is shallow. The paper does not discuss the stability of that nesting.

**The outer loop.** Inputs are the initial prior-trained surrogate $f^{L}=\mathcal{NN}^{L}$, a proposal density $q$, a subchain length $m$ ($m=1000$ in the experiments) and a maximum number of corrections $I_{\max}$. Each round runs $m-1$ subchain steps, proposes $z^{*}$, refines if required, computes the acceptance probability and accepts or rejects, and finally returns the pooled posterior samples.

> [!warning] Checking the source: numerator and denominator of the acceptance probability
> Equation (10) of the paper prints the high-fidelity acceptance probability as
> $\beta=\min\{1,\ \mathcal L(d,f^{H}(z^{-}))\pi(z^{-})/\mathcal L(d,f^{H}(z^{+}))\pi(z^{+})\}$.
> In Metropolis-Hastings with a symmetric proposal the proposed state $z^{+}$ belongs in the numerator. As printed, $\beta$ would be a **decreasing** function of the proposed state's posterior density: a proposal that fits the data better would be accepted less often and the chain would drift toward low-probability regions — which directly contradicts the paper's own stated purpose for this point, namely to place the test point in the posterior bulk. The printed form therefore cannot be what was intended. In context it should read
> $\beta=\min\{1,\ \mathcal L(d,f^{H}(z^{+}))\pi(z^{+})/\mathcal L(d,f^{H}(z^{-}))\pi(z^{-})\}$.
> This is a reader's reconciliation and not the source text; check the journal version before quoting.

### Theorems

**This paper proves nothing**, and states explicitly that it inherits the analytical setting of paper 37. The one verifiable design claim is a premise rather than a result: $Q$ must stay small, therefore the correction network's capacity must be restricted, as an anti-overfitting constraint.

### Numerical experiments

A benchmark elliptic inverse problem in two configurations.

| Item                       | Example 1                                                              |
| -------------------------- | ---------------------------------------------------------------------- |
| unknown                    | nine-parameter permeability                                            |
| prior surrogate            | 4 hidden layers of 40 neurons, 50 training points                      |
| correction network         | one hidden layer of 50 neurons                                         |
| refinement points per round | $Q=10$                                                                 |
| tolerance                  | $\mathrm{tol}\in\{0.1,\ 0.05\}$                                        |
| regularisation             | $\lambda=0$                                                            |
| subchain length            | $m=1000$                                                               |
| baselines                  | Direct (true-model MCMC), DNN (fixed prior-trained network), ADNN (this paper) |

Example 2 replaces the unknown by a high-dimensional random-field permeability parameterised by a Karhunen-Loève expansion; the verifiable material does not separately record its network sizes or tolerances.

**Results.** The fixed prior-trained network gives **visibly wrong** posterior marginals; the adaptive version recovers the true-model MCMC marginals; tightening the tolerance from $0.1$ to $0.05$ tightens the agreement further.

> [!note] This paper likewise has no transcribable numbers
> The verifiable material gives the configuration above and the qualitative comparison just stated, and **no error values, effective sample sizes or CPU timings**. Even the claim that $0.1\to0.05$ improves agreement is qualitative, with no accompanying numerical metric. This page therefore prints no results table for paper 49.

**What the experiments establish and where they fall short.** They establish an ablation-grade conclusion: same chain, same initial surrogate, the only difference being whether online refinement is on, and the version with it on matches true-model MCMC marginals while the version with it off does not. That attributes the improvement to the refinement mechanism rather than to network architecture or training tricks. What is missing: the motivation is high dimension (the polynomial basis blow-up), yet example 1 has only nine parameters and the dimension of example 2 is not recorded in the verifiable material; there is no quantitative metric at all, so trade-offs between tolerances or values of $Q$ cannot be compared; and the composite nesting deepens with the outer loop with its effect entirely unmeasured.

### Relation to the others

A direct successor to paper 37: the same outer MCMC skeleton and the same "retrain inside an $\ell^\infty$ ball around an accepted point" design, with polynomial chaos replaced by a composite network and the absolute indicator replaced by a relative one. The composite idea — feeding the previous surrogate into the next as an input — is this paper's signature contribution, and neither paper 79 nor [[en/computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|106]] adopts it; both fine-tune one network by transfer learning instead. The paper itself notes that extensions to paper 34 (ensemble Kalman) and [[en/computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|55]] (randomize-then-optimize) "are also possible".

## 79: a goal-oriented indicator and a greedy design

### The idea

The first three papers all ask the same question of their indicators: is the surrogate accurate at this point? That question has a hidden defect — **an inaccurate surrogate need not matter.** What inversion actually cares about is whether the data can be fitted. If the surrogate errs in a direction that does not move the data misfit, correcting it is wasted effort; conversely, a small error in a direction that dominates the misfit does matter. The first three papers align these two things indirectly, by putting the test point near the posterior bulk. This paper changes what is measured instead: **it measures the data misfit itself.**

Two consequences follow. On the good side, the indicator now shares units and scale with the inversion objective, so a threshold can be set on the relative change of the misfit without guessing the magnitude of the surrogate's error. On the cost side, a misfit is a scalar and does not say *in which direction* the surrogate is wrong, so the point-selection rule has to supply that information itself. That is what the greedy design does: it picks points that are far apart in the **surrogate's output space**, that is, points the surrogate itself treats as most different and is therefore most likely to be wrong about.

The second piece of intuition concerns the sampler. Unscented Kalman inversion needs only $2N_m+1$ forward evaluations per iteration and typically converges in a dozen or so, so the whole inversion costs a few hundred forward solves — two or three orders of magnitude below MCMC. That makes "one true solve per round for the indicator" a non-trivial fraction of the budget and puts a premium on keeping the indicator cheap. The greedy score uses only surrogate predictions on the candidate pool, so selection is essentially free, and true solves are spent only on the $Q$ points finally chosen.

### Setting

The infinite-dimensional Bayesian formulation: unknown $m\in\mathcal M$, forward operator $\mathcal G$, observations $y=\mathcal G(m)+\eta$ with $\eta\sim\mathcal N(0,\Sigma_\eta)$, and potential $\Phi(m;y)=\tfrac12\|y-\mathcal G(m)\|^{2}_{\Sigma_\eta}$. The surrogate $\widehat{\mathcal G}_t=\mathcal O\circ\mathcal F_\theta$ composes a DeepONet $\mathcal F_\theta$ with the observation operator $\mathcal O$, and $\nu_t$ is the surrogate-induced posterior approximation at round $t$.

### Derivation

**Step one: write the honest indicator, then reject it.** The local model error should be

$$
e_{M}(t):=\mathbb E_{\nu_t}\bigl\|\mathcal G-\widehat{\mathcal G}_t\bigr\|
=\Bigl(\int_{\mathcal M}\bigl|\mathcal G(m)-\widehat{\mathcal G}_t(m)\bigr|^{2}\,\nu_t(dm)\Bigr)^{1/2}.
$$

The paper rejects this as an **implementable** indicator because it needs a high-dimensional integral. The step is worth keeping because it fixes the status of every single-point indicator that follows: they are cheap stand-ins for this quantity, not estimators of it.

**Step two: the anchor.** Draw $T$ samples $\mathcal M^{(t)}=\{m^{(t)}_k\}_{k=1}^{T}$ from the current surrogate posterior $\nu_t$ and use the **true model** to pick the one with the best data fit:

$$
r_t=\arg\min_{m\in\mathcal M^{(t)}}\ \tfrac12\bigl\|y-\mathcal G(m)\bigr\|^{2}_{\Sigma_\eta}.
$$

**Step three: indicator and trigger.** The indicator is the data misfit at the anchor and the trigger is its **relative change**:

$$
e_{D}(t):=\Phi(r_t;y)=\tfrac12\bigl\|y-\mathcal G(r_t)\bigr\|^{2}_{\Sigma_\eta},
\qquad
\frac{\bigl|e_D(t)-e_D(t-1)\bigr|}{e_D(t)}\ >\ \epsilon .
$$

If it fails, refinement stops. Using a relative change rather than an absolute value carries a point: $e_D$ is bounded below by the noise level and can never reach zero, so the criterion has to be "is this still improving?" rather than "is this small enough?". The experiments use $\epsilon=0.01$ and $I_{\max}=10$.

**Step four: the greedy design.** Draw a large candidate pool $\Gamma=\{m_1,\dots,m_K\}$ from $\nu_t$ and grow a subset $\gamma_Q=\{\hat m_1,\dots,\hat m_Q\}\subset\Gamma$ one point at a time:

$$
\hat m_{j+1}=\arg\max_{m\in\Gamma\setminus\gamma_j}
\Bigl\{d\bigl(\widehat{\mathcal G}_t(m),\widehat{\mathcal G}^{j}_t\bigr)-\lambda\|m-r_t\|_{2}\Bigr\},
\qquad
d\bigl(\widehat{\mathcal G}_t(\cdot),\widehat{\mathcal G}^{j}_t\bigr)
=\max_{\hat m\in\gamma_j}\bigl\|\widehat{\mathcal G}_t(\cdot)-\widehat{\mathcal G}_t(\hat m)\bigr\|_{2},
$$

where $\widehat{\mathcal G}^{j}_t:=\{\widehat{\mathcal G}_t(\hat m_i)\}_{i=1}^{j}$ and $\lambda=1$ in every experiment. The two terms are deliberately adversarial: the first demands separation in the **surrogate's output space**, serving generalisation, and the second pulls selections back toward the anchor, serving local accuracy. The score uses only surrogate predictions, so selection is essentially free next to the full-order solves that follow.

**Step five: retraining is transfer learning.** The DeepONet parameters start from the previous round's weights rather than from scratch.

**Step six: unscented Kalman inversion.** The posterior approximation comes from a stochastic dynamical system

$$
m_{n+1}=r_0+\alpha(m_n-r_0)+\omega_{n+1},\quad \omega_{n+1}\sim\mathcal N(0,\Sigma_\omega),
\qquad
y_{n+1}=\mathcal G(m_{n+1})+\eta_{n+1},\quad\eta_{n+1}\sim\mathcal N(0,\Sigma_\eta),
$$

with regularisation parameter $\alpha\in(0,1]$. The prediction step is

$$
\hat r_{n+1}=\alpha r_n+(1-\alpha)r_0,
\qquad
\hat C_{n+1}=\alpha^{2}C_n+\Sigma_\omega ,
$$

and the analysis step follows from the Gaussian joint:

$$
r_{n+1}=\hat r_{n+1}+\hat C^{my}_{n+1}\bigl(\hat C^{yy}_{n+1}\bigr)^{-1}\bigl(y_{n+1}-\hat y_{n+1}\bigr),
\qquad
C_{n+1}=\hat C_{n+1}-\hat C^{my}_{n+1}\bigl(\hat C^{yy}_{n+1}\bigr)^{-1}\bigl(\hat C^{my}_{n+1}\bigr)^{T},
$$

with $\hat C^{yy}_{n+1}=\mathrm{Cov}[\mathcal G(m_{n+1})\mid Y_n]+\Sigma_\eta$. Expectations are evaluated by a modified unscented transform on $2N_m+1$ symmetric $\sigma$-points,

$$
m^{0}=r,\qquad m^{j}=r+c_j[\sqrt C]_j,\qquad m^{j+N_m}=r-c_j[\sqrt C]_j,\qquad 1\le j\le N_m .
$$

Each iteration costs $2N_m+1$ forward evaluations, and the paper reports convergence typically in $O(10)$ iterations.

**The closed loop.** Train $\mathcal F_\theta$ offline on a small number of prior samples; run unscented Kalman inversion against the current surrogate to obtain $\nu_t$ and $T$ samples from it (this step touches only the surrogate and the data); compute the anchor and $e_D(t)$ and test the relative change; if triggered, draw the pool, select $\gamma_Q$ greedily, run the full-order model on it, and fine-tune from the current weights; repeat at most $I_{\max}$ times.

### Theorems

**Theorem 2.1 (quoted from the paper's reference [43]).** Identical in shape to Theorem 2 of paper 37. With $\nu$ the true posterior, $\widehat\nu$ the surrogate-induced posterior, $\mathcal M(\epsilon)$ the $\epsilon$-feasible set and $\mathcal M^{\perp}(\epsilon)$ its complement, there exist $K_1,K_2>0$ with

$$
D_{\mathrm{KL}}(\widehat\nu\,\|\,\nu)\ \le\ \bigl(K_1\epsilon+K_2\,\nu(\mathcal M^{\perp}(\epsilon))\bigr)^{2}.
$$

This is the justification for demanding only local accuracy, and the paper cites it as motivation.

**Assumptions 3.2–3.4.**

- (3.2) for any $\epsilon$, a **linear** neural operator $\widehat G:\mathbb R^{N_m}\to\mathbb R^{N_y}$ can be trained with $\|\widehat G-G\|_2<\epsilon$;
- (3.3) $\|G\|_2<H$;
- (3.4) $G^{T}\Sigma_\eta^{-1}G\succ0$ with $\|G^{T}\Sigma_\eta^{-1}G\|_2>C_1$.

**Lemma 3.5.** Under 3.2–3.4, $\|\widehat G^{T}\Sigma_\eta^{-1}\widehat G\|_2>C_2$, by way of

$$
\bigl\|G^{T}\Sigma_\eta^{-1}G-\widehat G^{T}\Sigma_\eta^{-1}\widehat G\bigr\|_2\le 2\epsilon H\bigl\|\Sigma_\eta^{-1}\bigr\|_2 .
$$

The step rests on an elementary identity,

$$
A^{T}\Sigma^{-1}A-B^{T}\Sigma^{-1}B=A^{T}\Sigma^{-1}(A-B)+(A-B)^{T}\Sigma^{-1}B,
$$

with $A=G$ and $B=\widehat G$; bounding each term by $H\|\Sigma_\eta^{-1}\|_2\epsilon$ gives the result, and the reverse triangle inequality then yields the existence of $C_2$. (The identity and this route are routine steps supplied by this page; the paper states only the conclusion and asserts abstractly that some $C_2>0$ exists.)

**Theorem 3.6 (convergence in the linear case).** Hypotheses: Assumptions 3.2–3.4, $\mathrm{Range}(G^{T})=\mathrm{Range}(\widehat G^{T})=\mathbb R^{N_m}$, $\Sigma_\omega\succ0$, $\Sigma_\eta\succ0$. Conclusion: the surrogate-driven unscented Kalman inversion fixed point $(\widehat r_\infty,\widehat C_\infty^{-1})$ converges to the full-model fixed point $(r_\infty,C_\infty^{-1})$, with

$$
\bigl\|\widehat C^{-1}_\infty-C^{-1}_\infty\bigr\|_2\le\frac{2\epsilon H H_\eta}{1-\beta},
\qquad
\bigl\|\widehat r_\infty-r_\infty\bigr\|_2
\le\frac{K_1H_\eta H_y}{C_1}\Bigl(1+\frac{2(1+\alpha\beta)K_2H_\eta H^{2}}{(1-\beta)C_2}\Bigr)\epsilon ,
$$

where $\beta,C_1,C_2,K_1,K_2,H_\eta,H_y,H$ are positive bounded constants. Both bounds are $O(\epsilon)$, linear in the surrogate's operator-norm error.

**This is the only convergence result in the family that concerns the sampler's output.** Theorem 2 of paper 37 and Theorem 2.1 here bound the distance between posteriors; this one bounds the mean and covariance the algorithm actually computes.

> [!warning] The scope of the theorem
> Remark 1 of the paper explains that $\widehat G_\theta$ can be made linear, as Theorem 3.6 requires, by dropping the nonlinear activations in the branch net while keeping them in the trunk net. The theorem therefore covers a restricted DeepONet, not the network used in the nonlinear experiments. Quote it with that qualification attached.

### Numerical experiments

| Item                    | Setting                                                                          |
| ----------------------- | -------------------------------------------------------------------------------- |
| benchmarks              | Darcy flow, heat source inversion (Cases I and II), reaction-diffusion           |
| surrogate               | DeepONet with branch and trunk nets of 5 hidden layers × 100 neurons, `tanh`     |
| offline training        | $1\times10^{5}$ iterations on $N_{\mathrm{prior}}=1000$ Gaussian-random-field prior samples |
| data                    | $y_{\mathrm{obs}}=y_{\mathrm{ref}}+\max\{\lvert y_{\mathrm{ref}}\rvert\}\delta\xi$ |
| indicator threshold     | $\epsilon=0.01$                                                                   |
| retraining cap          | $I_{\max}=10$                                                                     |
| greedy weight           | $\lambda=1$                                                                       |
| sampler cost per round  | $2N_m+1$ forward evaluations, typically $O(10)$ rounds                           |
| baselines               | FEM-UKI (full model), DeepOnet-UKI-Direct (fixed surrogate), DeepOnet-UKI-Adaptive (this paper) |

Case I of the heat source inversion is verifiable in detail: a two-parameter source-location problem, DeepONet trained on $[0.5,1]\times[0.5,1]$ from 500 uniform samples, unscented Kalman inversion started at $[0.6,0.6]$; Case II is a higher-dimensional variant.

**The key experimental design is the pair of in-distribution and out-of-distribution truths.** The out-of-distribution case manufactures the failure mode "the surrogate is unreliable where the posterior sits" explicitly, and it is exactly where the fixed surrogate fails.

**Results (qualitative; the verifiable material contains no numerical table).** The fixed surrogate gives only a rough estimate and a visibly wrong inversion trajectory; the adaptive version drives the model error down monotonically with refinements and approaches full-order accuracy. In the third example refinement terminated after six iterations by the stopping criterion — evidence that the criterion actually fires rather than merely being cut off by $I_{\max}$.

**What the experiments establish and where they fall short.** They establish three things: the goal-oriented indicator can stop itself; greedy selection works across three structurally different PDEs; and under an out-of-distribution truth a fixed operator surrogate really does fail while the adaptive version recovers accuracy. What is missing: all experiments use the nonlinear DeepONet while Theorem 3.6 covers only a linearised branch net, so **there is a gap between theorem and experiment that the paper itself points out**; the monotone decrease of the model error is described qualitatively with no accompanying numerical sequence; and unscented Kalman inversion delivers a Gaussian approximation, so applicability to multimodal posteriors is untested — which is exactly where paper 106 picks up.

### Relation to the others

The most theoretically complete member of the family: it cites paper 37's $\epsilon$-feasible-set KL bound as motivation *and* adds a genuine convergence theorem for the surrogate-driven sampler in the linear case, which papers 34, 37, 49 and [[en/computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|55 and 56]] all lack. It replaces the pointwise $\ell^\infty$ indicators of papers 34, 37 and 49 with a goal-oriented indicator built from the data misfit itself, and random-in-a-ball sampling with an explicit diversity-versus-proximity trade-off. It is also the direct methodological ancestor of [[en/computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|106]], which keeps the adaptive fine-tuning loop while swapping unscented Kalman inversion for a variational flow, DeepONet for a Fourier neural operator, and greedy filtering for aggressive data replacement.

## Side-by-side comparison

| Component        | 37                              | 34                              | 49                               | 79                                    |
| ---------------- | ------------------------------- | ------------------------------- | -------------------------------- | ------------------------------------- |
| Surrogate        | multi-fidelity polynomial chaos | multi-fidelity polynomial chaos | composite multi-fidelity network | DeepONet                              |
| Sampler          | Metropolis-Hastings             | regularising ensemble Kalman    | Metropolis-Hastings              | unscented Kalman inversion            |
| Indicator        | absolute $\ell^\infty$          | relative $\ell^\infty$          | relative $\ell^\infty$           | relative change in the data misfit    |
| Test point       | accepted via the true model     | ensemble mean                   | accepted via the true model      | anchor chosen by the true model       |
| New points       | random in a shrinking ball      | random in a ball                | uniform in a ball                | greedy in output space plus proximity |
| Retraining       | refit least squares             | refit least squares             | train a shallow correction net   | transfer learning from prior weights  |
| Theory           | KL and Hellinger bounds         | none                            | none                             | fixed-point convergence, linear case  |
| Numerical results | setup and qualitative only     | full errors and solve counts    | setup and qualitative only       | setup and qualitative only            |

## Coverage check

| Item                                                        | Paper  | Status                                                                                     |
| ----------------------------------------------------------- | ------ | ------------------------------------------------------------------------------------------ |
| Verification tier and remaining gaps                        | all    | per-paper source, tier and unverified items                                                |
| Weighted least squares, Christoffel weight, row normalisation | 37   | definition of the weight, derivation that every row has norm $\sqrt M$, identity with the optimal-sampling topic |
| Index structure of the multi-fidelity correction            | 37     | correction term, merging as a pure regrouping, sample budget and the overfitting constraint |
| That the algorithm samples the surrogate's posterior        | 37     | role of the subchain, where the true model appears, contrast with asymptotic exactness     |
| $\epsilon$-feasible set with KL and Hellinger bounds        | 37, 79 | assumptions, conclusions, corollary, how the two terms match the two regions               |
| Full regularising ensemble Kalman update                    | 34     | extended-state reduction, update, covariances, regularisation choice, stopping rule        |
| The linear-case Tikhonov–Phillips limit                     | 34     | coincidence with the linear-Gaussian posterior mean, and the cost structure of being derivative-free |
| Quantitative comparison with out-of-prior truth             | 34     | errors, solve counts, CPU times, and the consistency check against the budget formula      |
| Composite multi-fidelity network                            | 49     | coupling form, training-set construction, contrast with a linear correction, capacity chain |
| Checking the acceptance-probability ratio                   | 49     | why the printed form cannot be what was intended                                           |
| Goal-oriented indicator and greedy design                   | 79     | the rejected honest indicator, anchor, relative change, greedy score and its cost structure |
| Unscented Kalman inversion and linear convergence           | 79     | prediction and analysis steps, $\sigma$-points, the identity behind Lemma 3.5, both $O(\epsilon)$ bounds and their scope |

## Sources for this page

- L. Yan and T. Zhou, [_Adaptive multi-fidelity polynomial chaos approach to Bayesian inference in inverse problems_](https://doi.org/10.1016/j.jcp.2018.12.025), J. Comput. Phys. 381 (2019), pp. 110-128 (preprint [arXiv:1807.00618](https://arxiv.org/abs/1807.00618)).
- L. Yan and T. Zhou, [_An adaptive multifidelity PC-based ensemble Kalman inversion for inverse problems_](https://doi.org/10.1615/Int.J.UncertaintyQuantification.2019029059), Int. J. Uncertain. Quantif. 9(3) (2019), pp. 205-220 (preprint [arXiv:1809.08931](https://arxiv.org/abs/1809.08931)).
- L. Yan and T. Zhou, [_An adaptive surrogate modeling based on deep neural networks for large-scale Bayesian inverse problems_](https://doi.org/10.4208/cicp.OA-2020-0186), Commun. Comput. Phys. 28 (2020), pp. 2180-2205 (preprint [arXiv:1911.08926](https://arxiv.org/abs/1911.08926)).
- Z. Gao, L. Yan, and T. Zhou, [_Adaptive operator learning for infinite-dimensional Bayesian inverse problems_](https://doi.org/10.1137/24M1643815), SIAM/ASA J. Uncertain. Quantif. 12(4) (2024), pp. 1389-1423 (preprint [arXiv:2310.17844](https://arxiv.org/abs/2310.17844)).
- Background sources: M. A. Iglesias, K. J. H. Law, and A. M. Stuart, [_Ensemble Kalman methods for inverse problems_](https://doi.org/10.1088/0266-5611/29/4/045001), Inverse Problems 29 (2013), 045001; A. M. Stuart, [_Inverse problems: a Bayesian perspective_](https://doi.org/10.1017/S0962492910000061), Acta Numerica 19 (2010), pp. 451-559.
