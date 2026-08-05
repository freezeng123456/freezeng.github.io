---
title: Samplers, Filters and Function-Space Flows
description: Papers 55, 56, 82, 88, 99 and 106 - replacing the sampler, with derivations, theorem hypotheses and experiments
lang: en
translation: computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering
tags:
  - paper-notes
  - bayesian-inverse-problems
  - variational-inference
---

> [!note] Coverage of this page
> Papers **55** (_J. Comput. Math._ 39(6), 2021), **56** (_Comput. Methods Appl. Mech. Engrg._ 386, 2021), **82** (_Comput. Phys. Commun._ 311, 2025), **88** (_Int. J. Mech. Sci._ 313, 2026), **99** (submitted to _Math. Comput._, [arXiv:2411.13277](https://arxiv.org/abs/2411.13277)) and **106** (submitted to _SIAM J. Sci. Comput._, [arXiv:2605.29373](https://arxiv.org/abs/2605.29373)).

The four papers on the [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|previous page]] fix the sampler and replace the surrogate. The six here do the reverse: the surrogate's role is comparatively stable and what changes is **how one advances on the posterior**.

## Verification tier

| No. | Source checked                                       | Tier                       | Remaining gap                                                                                        |
| --- | ---------------------------------------------------- | -------------------------- | ------------------------------------------------------------------------------------------------------ |
| 55  | arXiv:2104.06285 (full text)                         | **Full text**              | the explicit form of the approximate posterior $\widetilde\pi_{\mathrm{pos}}$ was lost and is unverified; the proposal density's normalising constant carries a typo |
| 56  | arXiv:2104.06276 (full text)                         | **Full text**              | whether the kernel exponent carries an extra $1/2$ is uncertain; the cost table's headers are partly damaged and the CPU-second columns are unverified |
| 82  | ScienceDirect abstract, full introduction, section openings | **Abstract and introduction only** | Sections 2 to 4 require a subscription and there are **no verifiable formulas**; no preprint exists |
| 88  | arXiv:2508.06852 (full text)                         | **Full text**              | the indexing and normalisation of the likelihood equation are reconstructed; the journal and preprint titles differ |
| 99  | arXiv:2411.13277 v3 (full text)                      | **Full text**              | the hypothesis list of Theorem 2.3 and the scalar coefficient of the Householder layer were not transcribed verbatim; experimental error values are unverified |
| 106 | arXiv:2605.29373 v1 (full text)                      | **Full text**              | conditions (i) and (ii) of the evidence-bound theorem were not transcribed verbatim; the displayed forms of equations (38) and (68) are unverified; **all numerical tables are illegible** |

**Paper 82 is the only one on this page below the full-text tier**, so it receives only the space its verifiable range supports: no formulas, no algorithm box, no derivation. The other five get full derivations and experimental configurations, though paper 106's experiments are limited to configuration and qualitative conclusions because its numerical tables could not be transcribed reliably.

## 55: putting a surrogate inside an optimisation proposal

### The idea

Randomize-then-optimize does something quite concrete geometrically, and it is worth saying first, because the formulas on their own are opaque.

In whitened coordinates the posterior density has the form $p(\theta\mid y)\propto\exp\bigl(-\tfrac12\|f(\theta)-y\|^{2}\bigr)$ with $f:\mathbb R^{n}\to\mathbb R^{m}$, $m\ge n$ (the prior term already folded into $f$). Write the residual map as $r(\theta)=f(\theta)-y$. As $\theta$ sweeps $\mathbb R^{n}$, $r(\theta)$ traces out an $n$-dimensional surface in $\mathbb R^{m}$, and **the posterior density at $\theta$ is just the standard Gaussian density on $\mathbb R^{m}$ evaluated at the corresponding point $r(\theta)$ of that surface.**

If the surface were a flat $n$-dimensional affine subspace, sampling would have a completely elementary solution: draw a standard Gaussian vector in $\mathbb R^{m}$, project it orthogonally onto the subspace, and read off the corresponding $\theta$. Because the orthogonal projection of a standard Gaussian onto a subspace is a standard Gaussian on that subspace, the resulting $\theta$ is an exact posterior sample.

RTO carries that recipe over to a curved surface. Take the thin QR factorisation of the Jacobian at a linearisation point and let $\bar Q$ be an orthonormal basis for the surface's tangent space there. Given a random vector $\xi$, RTO solves

$$
\bar Q^{T}\bigl(f(\theta)-(y+\xi)\bigr)=0,
$$

that is, it finds the point on the surface whose residual has the same $\bar Q$-coordinates as $\xi$. **When the surface is flat this is exact sampling** — which is precisely Remark 1 of Bardsley and coauthors: in the linear case the correction factor $c(\theta)$ is constant, nonlinear RTO reduces to linear RTO, and the proposal is the posterior itself.

Curvature spoils this in exactly two places, and those two places are exactly the two factors of $c(\theta)$:

$$
c(\theta)=\bigl|\bar Q^{T}J(\theta)\bigr|
\exp\Bigl(\tfrac12\bigl\|\tilde Q^{T}\bigl(f(\theta)-y\bigr)\bigr\|^{2}\Bigr).
$$

The first factor is **volume distortion**: where the surface tilts away from $\bar Q$, equal volume elements in $\xi$ map to unequal volume elements in $\theta$. The second is the **residual component nobody controlled**: the solve matches only the $\bar Q$-coordinates of the residual to $\xi$, the part orthogonal to $\bar Q$ is left free, and the posterior does penalise it. A Metropolis independence step then removes the discrepancy exactly, because the proposal density is explicit and the acceptance ratio collapses to a ratio of $c$.

Where the cost is. Every proposal requires solving a nonlinear least-squares problem, with $f$ and its Jacobian at each optimiser step, and the weight $c(\theta)$ needs another Jacobian and a determinant. Once the surrogate is a network the Jacobian comes from back-propagation and both become cheap. But the paper's real judgement is about something else: **where the surrogate should be trained.** RTO proposals already concentrate near the posterior — that is its entire advantage over a random walk — so a prior-trained surrogate is exactly mismatched. The paper therefore draws training points from an approximate posterior and uses "the same algorithm with prior-drawn points" as its ablation baseline.

### Setting

Gaussian prior $u\sim\mathcal N(u_{\mathrm{pr}},\Gamma_{\mathrm{pr}})$, whitened variable $u=S_{\mathrm{pr}}v+u_{\mathrm{pr}}$, with sampling carried out in $v$ and mapped back. $H(v)$ is the whitened stacked residual, $f$ the forward model, $n_{\mathrm{samps}}$ the number of samples.

### Derivation

**Step one: the linearisation point.**

$$
v_{\mathrm{ref}}=\arg\min_{v}\tfrac12\|H(v)\|^{2}.
$$

**Step two: the proposal.**

$$
v^{(i)}_{\mathrm{prop}}=\arg\min_{v}\ \tfrac12\bigl\|\bar Q^{T}H(v)-\xi^{(i)}\bigr\|^{2},
\qquad \xi^{(i)}\sim\mathcal N(0,I_n).
$$

**Step three: the proposal density is explicit.**

$$
\pi_{\mathrm{RTO}}(v)=(2\pi)^{-n/2}\bigl|{\det}\bigl(Q^{T}\nabla H(v)\bigr)\bigr|
\exp\Bigl(-\tfrac12\bigl\|Q^{T}H(v)\bigr\|^{2}\Bigr).
$$

**Step four: the acceptance ratio collapses to a ratio of weights.** Because the proposal does not depend on the current state, target and proposal densities pair up in the ratio:

$$
\frac{\pi_{\mathrm{tar}}(v^{(i)}_{\mathrm{prop}})\,\pi_{\mathrm{RTO}}(v^{(i-1)})}
{\pi_{\mathrm{tar}}(v^{(i-1)})\,\pi_{\mathrm{RTO}}(v^{(i)}_{\mathrm{prop}})}
=\frac{w\bigl(v^{(i)}_{\mathrm{prop}}\bigr)}{w\bigl(v^{(i-1)}\bigr)},
\qquad
w(v)=\bigl|{\det}\bigl(Q^{T}\nabla H(v)\bigr)\bigr|^{-1}
\exp\Bigl(-\tfrac12\|H(v)\|^{2}+\tfrac12\|Q^{T}H(v)\|^{2}\Bigr).
$$

> [!warning] Two notational hazards
> **First, a convention clash.** Up to the whitening, the $w$ used here is the **reciprocal** of the $c$ of Bardsley and coauthors, so the proposed sample sits in the **numerator** in the $w$ formulation while the previous sample sits in the numerator in Bardsley's ratio $r=\min\{1,\ c(\theta^{k-1})/c(\theta^{*})\}$. Both are correct in their own notation, and mixing them yields a reversed acceptance rule.
>
> **Second, a typo.** The preprint prints the normalising constant as $(2\pi)^{-\pi/2}$, plainly a slip for $(2\pi)^{-n/2}$. The display above uses the latter; that single symbol is this page's correction and everything else is as printed.

**Step five: scalable RTO reduces an $n$-dimensional optimisation to an $r$-dimensional one.** The paper restates the construction of Bardsley–Cui–Marzouk–Wang: take the rank-$r$ reduced singular value decomposition $\nabla f(v_{\mathrm{ref}})=\Psi\Lambda\Phi^{T}$, split $v_r=\Phi^{T}v$ and $v=\Phi v_r+v^{\perp}$, and then for each $\xi\sim\mathcal N(0,I_n)$,

$$
v^{\perp}=(I_n-\Phi\Phi^{T})\xi,
\qquad
v_r=\arg\min_{z}\bigl\|(\Lambda^{2}+I_r)^{-1/2}z+\Lambda\Psi^{T}f(v^{\perp}+\Phi z)-\Phi^{T}\xi\bigr\|^{2},
$$

with the determinant simplifying to

$$
\bigl|{\det}(\tilde Q^{T}\nabla H(v))\bigr|
=\bigl|{\det}(\Lambda^{2}+I_r)^{-1/2}\bigr|\cdot
\bigl|{\det}\bigl(I_r+\Lambda\Psi^{T}\nabla f(v)\Phi\bigr)\bigr| .
$$

The meaning is worth spelling out: **the orthogonal-complement direction needs no optimisation at all, being a plain Gaussian projection**, because the linearised model is insensitive to the data there; only the $r$-dimensional data-informed subspace, whose dimension is the numerical rank of the Jacobian, requires a solve. Cost is governed by $r$, not by $n$.

**Step six: goal-oriented training design.** The training points $\{v_i\}_{i=1}^{N}$ are drawn from an approximate posterior $\widetilde\pi_{\mathrm{pos}}$ rather than the prior. The ablation baseline "NN-RTO-pr" is the same algorithm with prior-drawn points. The explicit form of $\widetilde\pi_{\mathrm{pos}}$ was lost in the verifiable material and is not repeated here.

**Step seven: the parallel structure of the online stage.** Proposals and weights are independent across samples and can be computed fully in parallel; only the subsequent independence-Metropolis accept-reject has to run in series. That is RTO's second advantage over a random walk, distinct from its low autocorrelation.

Efficiency is measured by CPU-time-adjusted effective sample size:

$$
\mathrm{ESS}=\frac{n_{\mathrm{samps}}}{1+2\sum_{k=1}^{K}\rho(k)} .
$$

### Theorems

**This paper contains no theorem**, and the abstract's "converges to the direct RTO approach" is an argued and empirical statement rather than a proved proposition. The available theorems come from the work it restates, and are worth recording with their hypotheses.

**Theorem 3.1 (Bardsley–Solonen–Haario–Laine, the RTO density).** Hypotheses: $p(\theta\mid y)\propto\exp(-\tfrac12\|f(\theta)-y\|^{2})$; $f$ continuously differentiable with Jacobian $J(\theta)\in\mathbb R^{m\times n}$ of rank $n$ for every $\theta$; the MAP estimator $\bar\theta$ unique; $J(\bar\theta)=[\bar Q,\tilde Q]\bigl(\begin{smallmatrix}\bar R\\0\end{smallmatrix}\bigr)$; and $\bar Q^{T}J(\theta)$ invertible for all $\theta$. Conclusion: the induced density is $p_{\bar\theta}(\theta)\propto c(\theta)\,p(\theta\mid y)$ with $c$ as given in the intuition section above.

**The condition for uniform ergodicity.** By the general theory of independence Metropolis-Hastings, RTO-MH is uniformly ergodic if and only if $p(u\mid y)\le M\,p_{\mathrm{RTO}}(u\mid y)$ for some $M>0$, which holds exactly when $c(u)$ is bounded away from zero. **This is a rare checkable convergence statement in this family**, and the cost is that it does not carry over to the surrogate version: replacing $c$ by the surrogate's $c$ also replaces the target being sampled.

Two implementation notes from the same source: the determinant makes $c(\theta)$ numerically delicate and it can be very large, which is why RTO-MH is more stable in practice than RTO importance sampling; and samples whose residual exceeds a tolerance are rejected, with a tolerance of $10^{-8}$ in all the original examples.

### Numerical experiments

| Item              | Setting                                                                     |
| ----------------- | --------------------------------------------------------------------------- |
| forward problem   | Darcy-type elliptic inversion, recovering $\kappa(x)$ from noisy pointwise pressure |
| source            | $f(x)=100\sin(\pi x_1)\sin(\pi x_2)$                                        |
| discretisation    | bilinear Galerkin finite elements on a uniform $40\times40$ grid            |
| surrogate         | 3 hidden layers × 40 neurons                                                |
| training          | Adam, learning rate $5\times10^{-4}$                                        |
| training points   | $N\in\{50,100\}$                                                            |
| inner optimiser   | MATLAB `lsqnonlin` (trust-region-reflective Newton)                         |
| methods compared  | direct RTO (true model), NN-RTO-pr (prior-trained), NN-RTO (posterior-trained, this paper) |
| efficiency metric | CPU-time-adjusted effective sample size                                     |

**Results.** NN-RTO reproduces the one- and two-dimensional posterior marginals of direct RTO already at $N=50$ training points; NN-RTO-pr does not.

> [!note] The verifiable material contains no numbers
> The table above is the configuration and the result is qualitative: there are **no effective-sample-size values, no CPU timings and no error values** in the verifiable material, even though the paper defines a CPU-time-adjusted effective sample size as its efficiency metric. This page therefore prints no results table for paper 55.

**What the experiments establish and where they fall short.** They establish a clean ablation: three methods on the same forward problem, the same sampler and the same network architecture, with the sole difference being whether training points come from the prior or from an approximate posterior — and that one change decides whether the posterior marginals match. The credit goes to the location of the training points rather than to the network or the optimiser.

Three things are missing. First, the paper's headline claim is efficiency, yet its own well-defined efficiency metric has no reported values, so "significantly outperforms traditional RTO" can only be recorded here as an abstract-level claim. Second, the surrogate is trained once offline with no online refinement, so if the approximate posterior is itself badly off there is no self-correcting mechanism — untested here, and exactly the dividing line against the online-refinement methods of the [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|previous page]]. Third, Theorem 3.1 and the ergodicity condition above are stated for the true model $f$; whether hypotheses like "$\bar Q^{T}J(\theta)$ invertible everywhere" survive the replacement of $f$ by a network is not discussed.

### Relation to the others

It shares the "train the surrogate where the posterior lives" principle with papers [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|37, 49 and 79]] and with 56 and 106 on this page, but attaches it to an optimisation-based independence sampler rather than a random-walk chain or a particle method. The substantive difference is that **the surrogate here is trained once offline with no online refinement**. That makes it the static member of the family and the natural control case when arguing for the value of online refinement.

## 56: particle flow combined with online refinement

### The idea

Stein variational gradient descent differs from MCMC in more than efficiency. MCMC is one chain traversing the posterior in time; SVGD is a set of particles covering the posterior in space simultaneously, pushed there by a deterministic (kernelised) gradient flow. For many posteriors it needs far fewer samples than MCMC.

The price is $\nabla_x\log\pi(x)$ at every particle at every step. Under a PDE constraint that means one adjoint solve per particle per step, and if the gradient is unavailable at all (legacy code, non-differentiable solver) vanilla SVGD simply cannot run. Replacing the true model by a differentiable surrogate gets the gradient from back-propagation and removes that cost.

What is subtle is the failure mode. **A gradient flow does not hesitate.** MCMC with a bad surrogate is at least still random-walking and may wander elsewhere; SVGD with a bad surrogate will confidently push every particle into the **surrogate's** high-probability region and converge there beautifully. The paper's first experiment, a two-dimensional "double banana" posterior, puts exactly this on display: the fixed prior-trained network drives the particles into the wrong high-probability region, looking entirely converged. Online refinement is therefore not a refinement of the method but a precondition for it working at all.

Choosing the particle mean as the design point has a different justification here than on the previous page. The MCMC methods there need an extra true-model acceptance decision to obtain a trustworthy test point; here the particle set is itself a sample from the current approximate posterior, so its mean naturally sits in the high-probability region. By the same token new training points need not be drawn afresh — **they can be taken from the existing particles**, which are already in the right region. The only hazard is that they cluster and make the training problem ill-conditioned, which is what the separation constraint against the existing design prevents.

### Setting

Target density $\pi$ (unnormalised is enough), particles $\{x_i\}_{i=1}^{N}$ whose empirical measure $q(x)=\frac1N\sum_i\delta(x-x_i)$ should converge weakly to $\pi$. Forward model $f$, surrogate $\tilde f$, design set $\mathcal D$, inner-loop length $T$, outer-loop cap $I_{\max}$.

### Derivation

**Step one: write sampling as transport.** The particle update is

$$
x_i\ \leftarrow\ x_i+\varepsilon\varphi(x_i),\qquad i=1,\dots,N,
$$

with step size $\varepsilon$ and velocity field $\varphi:\mathbb R^{d}\to\mathbb R^{d}$. The question becomes how to choose $\varphi$.

**Step two: choose the steepest descent direction.**

$$
\varphi^{*}=\arg\max_{\varphi\in\mathcal S}
\Bigl\{-\tfrac{d}{d\varepsilon}\mathrm{KL}\bigl(q_{[\varepsilon\varphi]}\,\|\,\pi\bigr)\Bigr\},
$$

where $q_{[\varepsilon\varphi]}$ is the law of $x'=x+\varepsilon\varphi(x)$.

**Step three: the Stein operator identity turns that derivative into an expectation.** This is the hinge of the whole derivation:

$$
-\frac{d}{d\varepsilon}\mathrm{KL}\bigl(q_{[\varepsilon\varphi]}\|\pi\bigr)\Big|_{\varepsilon=0}
=\mathbb E_{x\sim q}\bigl[\mathcal A_{\pi}^{\top}\varphi(x)\bigr],
\qquad
\mathcal A_{\pi}^{\top}\varphi(x)=\nabla_x\log\pi(x)^{\top}\varphi(x)+\nabla_x^{\top}\varphi(x).
$$

Two consequences. The right-hand side is a **linear functional** of $\varphi$, so maximising it over a convex set has a chance of a closed form. And $\mathcal A_\pi$ depends on $\pi$ only through $\nabla_x\log\pi$, so **the normalising constant never appears** — which is why SVGD applies directly to unnormalised posteriors.

**Step four: maximising over the unit ball of a reproducing kernel Hilbert space gives a closed form.** With $\mathcal S=\{\varphi\in\mathcal H^{d}:\|\varphi\|_{\mathcal H^{d}}\le1\}$,

$$
\varphi^{*}(\cdot)\ \propto\ \mathbb E_{x\sim q}\bigl[\mathcal A_{\pi}\kappa(x,\cdot)\bigr]
=\mathbb E_{x\sim q}
\bigl[\nabla_x\log\pi(x)\,\kappa(x,\cdot)+\nabla_x\kappa(x,\cdot)\bigr],
$$

and the maximal value of the objective is the kernelised Stein discrepancy between $q$ and $\pi$.

**Step five: substitute the empirical measure to get an executable iteration.**

$$
x^{(l+1)}_{i}\leftarrow x^{(l)}_{i}+\varepsilon_l Q_l\bigl(x^{(l)}_i\bigr),
\qquad
Q_l(x)=\frac1N\sum_{j=1}^{N}
\Bigl[\nabla_{x^{(l)}_j}\log\pi\bigl(x^{(l)}_j\bigr)\kappa\bigl(x^{(l)}_j,x\bigr)
+\nabla_{x^{(l)}_j}\kappa\bigl(x^{(l)}_j,x\bigr)\Bigr].
$$

The two terms divide the labour cleanly: the first is a kernel-weighted ascent of $\log\pi$ and attracts, the second is a kernel gradient and repels, preventing collapse. A limiting case confirms the reading: with $N=1$ and a kernel satisfying $\nabla_x\kappa(x,x)=0$ the repulsion vanishes and SVGD degenerates to gradient ascent on the MAP point.

The kernel is the radial basis $\kappa(x,x')=\exp(-h\|x-x'\|^{2})$. Step sizes use AdaGrad, and the paper notes explicitly that a constant step can diverge in high dimension.

> [!warning] The kernel's constant
> The extracted source writes the kernel with a $/2$ attached, so whether the exponent carries an extra $1/2$ is uncertain in the verifiable material and the constant above is this page's reconstruction. It affects no conclusion ($h$ is a tunable bandwidth) but should be checked against the source before quoting verbatim. In the experiments the radial-basis bandwidth used for maximum-mean-discrepancy scoring is set by the median heuristic on the reference samples.

**Step six: which term is actually expensive in an inverse problem.** The log-posterior gradient splits as

$$
\nabla_x\log\pi(x)=\nabla_x\log p_0(x)-\nabla_x\Phi(x;y),
$$

with the prior term usually available in closed form. The expensive part is $\nabla_x\Phi$, which needs the Jacobian of the forward map. **That is exactly what the surrogate replaces.**

**Step seven: design point and indicator.** The design point is the particle mean and the indicator is the relative $\ell^2$ error there:

$$
x^{*}=\frac1M\sum_{i=1}^{M}x^{(t+1)}_{i},
\qquad
\mathrm{err}(x^{*})=\frac{\bigl\|f(x^{*})-\tilde f(x^{*})\bigr\|_{2}}{\bigl\|f(x^{*})\bigr\|_{2}} .
$$

(The printed formula mixes two index bounds $N$ and $M$; the display above follows the evident intent, the mean over the current particle set.)

**Step eight: selection with a separation constraint.** New training points are taken from the existing particles:

$$
x_i=\arg\min_{x'\in X^{(t+1)}}\|x'-x^{*}\|_{2}
\quad\text{s.t.}\quad \|x'-x\|_{2}\ge R\ \ \forall x\in\mathcal D .
$$

The objective asks for closeness to $x^{*}$ and the constraint asks for separation of at least $R$ from the whole existing design; together they push selections as close to the high-probability region as non-degeneracy allows. If no feasible minimiser exists the loop breaks and the radius shrinks, $R\leftarrow\rho R$ with $0<\rho<1$ — note this shrinking means the opposite of the previous page's: here it fires because no new point can be found at the current radius, which says the design is already dense enough.

**Step nine: the outer loop and its cost.** Run $T$ inner SVGD steps, test for refinement, and repeat at most $I_{\max}$ times; online retraining is transfer learning from the previous network. Total online true-solve cost is

$$
N_{\mathrm{eval}}=\sum_t q_t,
$$

where $q_t$ is the number of points actually accepted from the ball $\mathcal B(x^{*},R_t)$ at outer step $t$.

### Theorems

**This paper proves no convergence theorem.** The abstract's "refine the local approximation online without destroying the convergence of the resulting SVGD" states the design goal; the substantive claim is that refinement is infinite in principle, so the particles migrate to the true target rather than to the surrogate's.

That deserves comparison with the framework it imitates. The local-approximation MCMC of Conrad–Marzouk–Pillai–Smith **does prove asymptotic exactness**, and what carries the proof is a set of trust-region error bounds for local polynomial regression (from Conn–Scheinberg–Vicente). For local linear and quadratic fits on a ball $\mathcal B(\theta,R)$,

$$
\bigl|f_i(\theta')-\bigl(\mathcal L^{\sim j}_{\mathcal B(\theta,R)}(\theta')\bigr)_i\bigr|\le\kappa_l(\nu_1,\lambda,d)\,R^{2},
\qquad
\bigl|f_i(\theta')-\bigl(\mathcal Q^{\sim j}_{\mathcal B(\theta,R)}(\theta')\bigr)_i\bigr|\le\kappa_q(\nu_2,\lambda,d)\,R^{3},
$$

where $\nu_1,\nu_2$ are Lipschitz constants and $\lambda$ is a **poisedness** constant measuring the geometry of the design. As long as $\lambda$ stays bounded, shrinking $R$ provably reduces the error — that is the trust-region mechanism.

**Paper 56 keeps that skeleton — its Algorithm 2 is explicitly a sketch in the same style — while replacing the components that make the bounds work**: local polynomial regression becomes a globally parameterised network with local retraining, and the cross-validation indicator becomes a far cheaper single-point relative error. A network has no analogue of a poisedness constant, so the $R^{2}$ and $R^{3}$ bounds do not transfer, and neither does the asymptotic exactness. That judgement of the relationship is this page's, and it explains why the paper can only state "without destroying convergence" as a goal rather than as a theorem.

### Numerical experiments

Three problems:

| Ex. | Problem                                          | Role                                                                    |
| --- | ------------------------------------------------ | ------------------------------------------------------------------------- |
| 1   | two-dimensional "double banana" posterior        | exhibits the failure mode of a gradient flow converging confidently to the wrong mode |
| 2   | two-dimensional heat source inversion            | a standard PDE inverse problem                                          |
| 3   | diffusion coefficient in a time-fractional PDE   | the only example that left transcribable numbers                        |

Default parameters:

| Parameter              | Value                                                                        |
| ---------------------- | ---------------------------------------------------------------------------- |
| refinement points per round | $Q=5$                                                                   |
| initial radius         | $R=0.2$                                                                      |
| tolerance              | $10^{-2}$                                                                    |
| radius shrink factor   | $\rho=0.8$                                                                   |
| outer iteration cap    | $I_{\max}=30$                                                                |
| inner iterations       | $T=10$                                                                       |
| momentum               | $0.9$                                                                        |
| sample quality score   | maximum mean discrepancy, radial-basis bandwidth by median heuristic, shared across methods |

**Configuration and outcome of example 1.** The fixed prior-trained network uses 10 training points and 3 hidden layers of 20 neurons; it drives the particles into the **wrong** high-probability region, while the locally refined version spreads them over the true support after 100 iterations and matches at the last iteration.

**Quantitative comparison in example 3.**

| Configuration | Method                       | Online true solves | Error  |
| ------------- | ---------------------------- | ------------------ | ------ |
| $n_t=100$     | fixed prior-trained network  | 0                  | 0.3968 |
| $n_t=100$     | local refinement (this paper) | 50                | 0.1732 |
| $n_t=500$     | fixed prior-trained network  | 0                  | 0.2941 |
| $n_t=500$     | local refinement (this paper) | 80                | 0.1155 |

Each configuration is an internally matched comparison: $n_t$ is equal within a row pair and the only difference is whether online refinement is on. The result is more than a halving of the error ($0.3968\to0.1732$ and $0.2941\to0.1155$) for 50 and 80 online true solves.

> [!warning] How far this table can be read
> Three qualifications have to travel with it. First, the definition of $n_t$ could not be recovered from the verifiable material; what is confirmed is that it takes the same value across the two methods being compared, so the row-wise comparison is meaningful. Second, by the paper's stated scoring the error column should be a maximum mean discrepancy, but the extracted headers are partly damaged, so that identification is an **inference** and not a transcription. Third, the same table's cost cell for vanilla SVGD extracts as "$500\times300$" alongside a total of about $3710$, two figures that do not agree (the product is far larger than the stated total), so this page **does not use it as a baseline** and does not repeat the number. The CPU-second columns are likewise unverifiable.

**What the experiments establish and where they fall short.** They establish two things: running SVGD against a fixed prior-trained surrogate converges to the wrong answer with no visible sign of trouble (example 1 demonstrates this cleanly), and a few dozen online true solves halve the error (example 3 gives the numbers). What is missing: only one example left numbers at all; how the benefit of refinement scales with $n_t$ cannot be judged from two data points; the paper's own key phrase, "without destroying convergence", has no corresponding measurement, since a lower maximum mean discrepancy is not the same as convergence to the true target; and, as argued above, the trust-region bounds that make convergence provable in the CMPS framework do not apply to a network surrogate.

### Relation to the others

The particle-method member of the family, and the closest of all of them in spirit to Conrad–Marzouk–Pillai–Smith: it copies the "sketch of an approximate algorithm, refine, repeat" structure while replacing local polynomial regression by a globally parameterised network with local retraining, and the cross-validation indicator by a much cheaper single-point relative error. The design-point choice, the particle mean, plays the role that the accepted chain state plays in papers [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|37 and 49]] and the anchor plays in [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|79]]. Paper 106 turns the relationship around and uses SVGD as a **baseline** rather than a tool.

## 82: two-tier learning in a time-dependent setting

> [!warning] What could be verified
> Sections 2 to 4 of this paper require a subscription and no preprint was found. What follows comes from the abstract, the full introduction and ScienceDirect's section openings; it contains **no verifiable formulas**, and the algorithm outline is reconstructed from the abstract and introduction rather than read from an algorithm box. This section therefore has no derivation and no theorems subsection: those are not omitted for brevity, there is simply nothing in the verifiable material to put in them.

### The idea

Cost accounting in data assimilation differs from the static case: ensemble Kalman filtering needs one PDE solve per member per assimilation step, and assimilation steps keep arriving, so the total grows linearly in time. Shrinking the ensemble introduces sampling error in the covariance estimate, which localisation only partly mitigates.

Reduced-order models are the standard escape, and for parametric dynamical systems they often fail: classical reduction rests on linear mode superposition, while the solution manifolds of such systems tend to have slowly decaying Kolmogorov width, so the reduced dimension needed for a given accuracy keeps growing and the benefit evaporates. Data-driven deep reduced-order models are not restricted to a linear subspace, but **they move the cost rather than remove it** — generating the offline snapshots still means many runs of the expensive full-order solver.

That last cost is this paper's target. The idea is to generate snapshots on a **coarse grid**, accept the resulting bias, and learn that bias separately from a small number of **fine-grid, short-horizon** runs. What makes this viable is an empirical observation the paper reports: the reduced-order error stays relatively stable across the temporal domain, so a short fine-grid window suffices to characterise it. The correction is applied on the **observation output** side, which is where the Kalman update actually consumes the model, so the analysis step sees a debiased predicted observation.

### Setting

Joint state and parameter estimation for time-dependent parametric PDEs. The full-order model lives on a fine grid $h$ while the reduced model's training snapshots are generated on a coarse grid $H$ with $H\gg h$. The error metric is

$$
\mathrm{rel}_{L^2}(\mu)=\frac{\|\mu-\mu^{*}\|_2}{\|\mu^{*}\|_2},
$$

together with an analogous relative $L^2$ error for the state $u$ (the source snippet is truncated mid-formula).

### The two tiers

**Tier one: operator inference on the coarse grid.** The reduced model is given a **polynomial** form mirroring the structure of the governing equations, and the reduced operators are learned from simulation data by least-squares regression rather than by intrusive projection. The only departure from standard operator inference is that the snapshots come from the coarse grid. That is simultaneously why the offline stage is cheap and why the reduced model is biased.

**Tier two: a model error network.** Networks map **time and parameters** to the reduced-order discrepancy, justified by the stability of those errors across the temporal domain. The corrected object is called D-ROM, and embedding it in the ensemble Kalman filter gives DR-EnKF; the correction acts on the output side, consistent with the prior work the paper cites (Pagani and coauthors, who used kriging interpolation to calibrate a reduced-order input-to-output map).

The paper's novelty claim is that this is the first operator-inference-based ensemble Kalman approach for simultaneous state and parameter estimation in nonlinear time-dependent PDEs.

**Outline (reconstructed from the abstract and introduction, not read from an algorithm box):** long-horizon coarse-grid simulations feed operator inference to learn the reduced operators; short-horizon fine-grid simulations give the reduced-order discrepancy and train the error network; online, the reduced model propagates the ensemble while the error network corrects each member's predicted observation before the analysis step.

### Numerical experiments

| Item             | Content                                                                                       |
| ---------------- | ----------------------------------------------------------------------------------------------- |
| benchmarks       | Burgers' equation, the FitzHugh-Nagumo model, advection-diffusion-reaction systems             |
| comparison       | DR-EnKF (this paper), R-EnKF (uncorrected operator-inference reduced filter), full-order ensemble Kalman filter |
| metrics          | relative $L^2$ error for the parameter and for the state                                       |
| reported outcome | DR-EnKF matches full-order accuracy at a fraction of the cost, while R-EnKF is biased          |

**This section stops here, and has to.** The abstract reports "considerable computational speedup without compromising accuracy"; **the specific speedup factor, error values, grid sizes, ensemble sizes and network architectures are all outside the verifiable range.** Any more detailed account would be invention. Obtaining the paper's operator-inference least-squares loss or its error-network architecture requires the published full text.

### Relation to the others

Structurally this is paper [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|34]] with the ingredients swapped: a cheap model inside a Kalman method plus an explicit correction for the surrogate's error. The difference is **when** the correction is learned. Papers 34, 37, 49, 56 and 79 all refine online, driven by an indicator; this one learns the discrepancy **offline** as a function of time and parameters and applies it online with no further refinement, which puts it alongside paper 55 on this page. The paper cites paper 34, and companion work by the same group (Wang–Li–Yan on adaptive ensemble Kalman inversion with statistical linearisation, and a paper on continuous-data-assimilation-enhanced reduced-order filtering) makes the lineage explicit.

## 88: hierarchy and staging against non-identifiability

This is the one paper in the topic that is not a PDE inverse problem, and that is exactly what makes it useful for comparison.

### The idea

Two structural difficulties block single-level inference when a mesoscopic red blood cell model is aligned with experiment, and they need two different remedies.

**The first is cross-platform.** The data come from different experimental platforms — interference microscopy for geometry, optical tweezers for stretching, fluctuation and relaxation assays for viscoelasticity — whose systematic discrepancies are mutually inconsistent. Treating them as independent observations under one likelihood forces the model into a compromise between contradictory biases, and that compromise corresponds to no real cell. **Hierarchical structure addresses this one**: hyperparameters let each dataset keep its own discrepancy while information about the underlying physical parameters is shared.

**The second is non-identifiability.** No single dataset pins down all five parameters at once: the equilibrium shape fixes the geometry but is nearly insensitive to viscosity, while the relaxation assay is sensitive to viscosity but depends on the geometry being known. **Staging addresses this one**: fix the parameters that the cleanest experiments do determine (equilibrium shape and stretching pin down geometry and shear modulus), then carry them forward as information into a second stage that determines the viscoelastic parameters from fluctuation and relaxation. The order is not arbitrary; it follows identifiability.

On top of this the forward model is a dissipative particle dynamics simulation, far too slow for direct MCMC, so a surrogate is unavoidable — but here the surrogate is the plainest kind: trained offline, never refined, with no error indicator. **The value of this paper is not in its surrogate machinery but in showing that the surrogate-accelerated Bayesian skeleton moves outside PDE inversion.**

### Setting

The inferred parameters are

$$
\vartheta_{\mathrm{in}}=(A_0,\ v,\ \mu_{sh},\ k_b,\ \eta_m),
$$

that is desired surface area, reduced volume, shear stiffness, bending stiffness and membrane viscosity. The simulation outputs are

$$
y=(D_{eq},\ h_{\max},\ h_{\min},\ D_{ax},\ D_{tr},\ t_c,\ W_{fl}),
$$

that is equilibrium diameter, maximum and minimum thickness, axial and transverse deformation under stretching, relaxation characteristic time and a membrane-fluctuation summary. The observation model lumps all model error and uncertainty into one zero-mean Gaussian term, $y_j=\text{model}_j(x;\vartheta)+\sigma\epsilon_j$ with $\epsilon_j$ of standard deviation $\sigma$. (The indexing and normalisation of that equation are reconstructed.)

The forward model is the viscoelastic red blood cell model of Pivkin and coauthors, built on dissipative particle dynamics: the cell is a two-dimensional triangular network of $N_v$ vertices and $N_b$ bonds (500 vertices in the simulations) governed by viscoelastic WLC-POW bonds, harmonic dihedrals, and area and volume constraints. Inference uses the Korali package.

### Three structural decisions

**Dynamic annealing fixes the stress-free baseline.** Cells at different reduced volumes $v$ are deformed from the biconcave shape at $v=0.64$; at larger $v$ the fixed triangular mesh retains large residual bond forces, so whenever $|v-v_{\mathrm{desired}}|$ exceeds a tolerance $\epsilon_v$ the equilibrium bond length $l^{\mathrm{bond}}_{eq}$ is moved toward the current bond length $l^{\mathrm{bond}}_i$. This defines the stress-free baseline all four experiments are referenced to — without it the four experiments would effectively be compared against different reference configurations.

**Eight surrogate networks.** Four experiments times two cell populations (healthy hRBCs and _Plasmodium falciparum_-infected Pf-RBCs).

**Two-stage hierarchical inference.** Stage I performs hierarchical inference over the single-level models for the equilibrium and stretching tests, giving stable distributions for the geometric parameters and the shear modulus. Stage II carries that information forward, adds the fluctuation and relaxation tests, and identifies the full parameter set including the viscoelastic parameters. Hyperparameters $\psi$ encode information sharing across datasets and $\sigma_{\mathrm{out}}$ are the output standard deviations.

The full pipeline: establish priors from microscopic simulations and literature; apply dynamic annealing to reach a genuine stress-free baseline; run dissipative particle dynamics simulations replicating the four representative experiments and collect $(\vartheta_{\mathrm{in}},y)$ pairs; train the eight surrogates and validate against held-out simulations; run Stage I; run Stage II; validate the posterior predictive $y_{\mathrm{new}}$ against the experimental datasets.

### Numerical experiments

Four experiments: the stress-free equilibrium shape, membrane fluctuation on an adhered substrate, deformation under optical-tweezer stretching, and relaxation after force release. Two cell populations: healthy cells and _Plasmodium falciparum_-infected cells.

**Surrogate configuration:**

| Item              | Setting                                        |
| ----------------- | ---------------------------------------------- |
| number of networks | 8 (4 experiments × 2 cell types)              |
| architecture      | three hidden layers, `tanh`, mean-squared-error loss |
| optimiser         | Adam with a step learning-rate schedule        |
| training data     | 10,000 simulation results per network          |
| width selection   | per experiment, guided by sensitivity analysis |
| reported accuracy | prediction errors below $10^{-2}$              |

**Experimental anchors (from the published literature):**

| Quantity   | Healthy cells (Evans & Fung)  | Infected cells (from surface-area and volume measurements) |
| ---------- | ----------------------------- | ---------------------------------------------------------- |
| $D_{eq}$   | $7.82\pm0.62\ \mu\mathrm m$   | $6.9\ \mu\mathrm m$                                        |
| $h_{\max}$ | $2.58\pm0.27\ \mu\mathrm m$   | $3.2\ \mu\mathrm m$                                        |
| $h_{\min}$ | $0.81\pm0.35\ \mu\mathrm m$   | $3.2\ \mu\mathrm m$                                        |

The stretching data ($D_{ax}$ and $D_{tr}$ against applied force) come from Mills and coauthors. That $h_{\max}=h_{\min}$ for the infected cells is itself a morphological statement: the cell moves from a biconcave disc toward a sphere.

**Outcome.** Posteriors are statistically robust; pathological cells are inferred to be **stiffer and more viscous**; cross-platform data fusion mitigates the multi-source uncertainty that defeats single-level inference; posterior predictions are consistent with the experimental observations. The paper proves no mathematical theorem.

**What the experiments establish and where they fall short.** They establish that the framework runs on a real, non-PDE calibration problem and yields physically meaningful conclusions. What is missing: the verifiable material gives only the order of magnitude of surrogate accuracy (below $10^{-2}$) and the qualitative conclusions above, with **no posterior interval values, no per-parameter identifiability comparison and no quantitative contrast against single-level inference**. The central claim that hierarchy and staging are necessary can therefore only be recorded here as the paper's reported conclusion, not checked quantitatively.

> [!note] Title difference
> The preprint is titled _An RBC-MsUQ Framework for Red Blood Cell Morpho-Mechanics_, with the same authors and abstract; the journal title is the one listed on the homepage. This site records the journal version, and both should be given when citing.

### Relation to the others

The outlier of the set, and useful for exactly that reason. It shares only the surrogate-accelerated Bayesian skeleton with papers [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|34, 37, 49 and 79]] and with 55 and 56 on this page — the surrogate is a plain feed-forward network trained offline, with no online refinement, no error indicator and no multi-fidelity correction. What it adds that none of the others has is **hierarchical** structure across heterogeneous datasets and a **staged** decomposition of an otherwise non-identifiable inverse problem. The closest cited relative is Economides and coauthors on hierarchical Bayesian uncertainty quantification for a red blood cell model.

## 99: variational inference directly in function space

### The idea

Discretising first and then doing variational inference makes the algorithm's behaviour depend on the mesh: the same hyperparameters behave differently on a finer grid, and iteration counts and accuracy drift as the mesh is refined. Defining the algorithm directly on function space removes that dependence but runs into an obstruction with no finite-dimensional analogue.

In finite dimensions any smooth bijection pushes one density to another, Lebesgue measure supplies a common reference, and the KL divergence is always defined. On an infinite-dimensional space there is no Lebesgue measure, so "density" can only mean a Radon-Nikodym derivative against some reference measure — and **two probability measures on an infinite-dimensional space are generically mutually singular**: they put their mass on disjoint sets and the derivative does not exist. The KL loss is then not large but **undefined**; the training signal is not poor, it is absent.

The paper's diagnosis is specific: the Radon-Nikodym derivative contains $\mathcal H$-norm terms such as $\langle f_\theta(u)-u,\ f_\theta(u)-u\rangle_{\mathcal H}$, and if $f_\theta(u)-u$ is not regular enough these are infinite. So requiring the transformed measure to be equivalent to the prior really functions as a device that **forces every flow layer to be regular enough**. That explains a design choice that otherwise looks odd: the flow layers are restricted to be low rank. Low rank is a computational convenience in neural operators; here it is a theoretical necessity, because it confines $f_\theta(u)-u$ to a finite-dimensional smooth subspace spanned by prior covariance eigenfunctions, which is what keeps those inner products finite.

The second piece of intuition is about how the loss is computed. It appears to require expectations under the transformed measure $\mu_{f_\theta}$, which is only reachable through the flow. But $\mu_{f_\theta}$ is the pushforward of $\mu_0$ under $f_\theta$, so a change of variables moves every expectation back onto the **prior**, which is a Gaussian measure one can sample directly. Each training step therefore only needs prior samples pushed through the flow and a Radon-Nikodym derivative evaluated — the posterior is never touched.

### Setting

An unknown $u$ in a Hilbert space $\mathcal H_u$ with Gaussian prior $\mu_0$, potential $\Phi$, and posterior $\mu$ given by the infinite-dimensional Bayes formula

$$
\frac{d\mu}{d\mu_0}(u)=\frac{1}{Z_\mu}\exp\bigl(-\Phi(u)\bigr).
$$

Every measure in the approximating family $\mathcal M(\mathcal H_u)$ is required to be **equivalent** to $\mu_0$ (following Pinski and coauthors), so that the Radon-Nikodym derivative is defined.

### Derivation

**Step one: write the KL objective as a difference of two derivatives against the prior.**

$$
D_{\mathrm{KL}}(\nu\|\mu)=\int_{\mathcal H_u}\ln\frac{d\nu}{d\mu}(u)\,\nu(du)
=\int_{\mathcal H_u}
\Bigl[\ln\frac{d\nu}{d\mu_0}(u)-\ln\frac{d\mu}{d\mu_0}(u)\Bigr]\nu(du),
$$

minimised over $\nu\in\mathcal M(\mathcal H_u)$. This only inserts the common reference $\mu_0$, but it turns the problem into "both derivatives against the prior must exist and be computable".

**Step two: Theorem 2.3 supplies both.** For a composition $f_\theta=f^{(N)}_{\theta_N}\circ\cdots\circ f^{(1)}_{\theta_1}$ whose layers are bijective and satisfy the paper's conditions, $\mu_{f_\theta}\sim\mu_0$ and

$$
\frac{d\mu_{f_\theta}}{d\mu_0}\bigl(f_\theta(u)\bigr)
=\prod_{n=1}^{N}\bigl|{\det}_{1}\bigl(Df^{(n)}_{\theta_n}(u_{n-1})\bigr)\bigr|^{-1}
\exp\Bigl(\tfrac12\bigl\langle f_\theta(u)-u,\ f_\theta(u)-u\bigr\rangle_{\mathcal H}
+\bigl\langle u,\ u-f_\theta(u)\bigr\rangle_{\mathcal H}\Bigr),
$$

where ${\det}_1$ is the Carleman-Fredholm determinant, the object that replaces the usual Jacobian determinant in infinite dimensions. The product form says the layers contribute multiplicatively as in finite-dimensional normalizing flows; what is new are the two $\mathcal H$ inner products, which are specific to infinite dimensions and are the source of the regularity requirement.

**Step three: change of variables pushes the expectations onto the prior.** Since $\mu_{f_\theta}$ is the pushforward of $\mu_0$ under $f_\theta$, for any integrable $g$ one has $\mathbb E_{\mu_{f_\theta}}[g(u)]=\mathbb E_{\mu_0}[g(f_\theta(u))]$. Applied to the first term:

$$
D_{\mathrm{KL}}(\mu_{f_\theta}\|\mu)
=\mathbb E_{\mu_{f_\theta}}\ln\frac{d\mu_{f_\theta}}{d\mu_0}(u)
-\mathbb E_{\mu_{f_\theta}}\ln\frac{d\mu}{d\mu_0}(u),
\qquad
\mathbb E_{\mu_{f_\theta}}\ln\frac{d\mu_{f_\theta}}{d\mu_0}(u)
=\mathbb E_{\mu_0}\ln\frac{d\mu_{f_\theta}}{d\mu_0}\bigl(f_\theta(u)\bigr).
$$

**This step is why the algorithm is executable at all**: the prior can be sampled directly.

**Step four: Monte Carlo gradient and iteration.** Initialise $\theta_0$; at each step draw $N$ samples $u_i\sim\mu_0$, push them through the flow, and update $\theta_{k+1}=\theta_k-\alpha_k\nabla_{\theta_k}L(\theta_k)$ with

$$
\nabla_{\theta_k}L(\theta_k)\approx
\frac1N\sum_{i=1}^{N}\nabla_{\theta_k}\ln\frac{d\mu_{f_{\theta_k}}}{d\mu_0}\bigl(f_{\theta_k}(u_i)\bigr)
-\frac1N\sum_{i=1}^{N}\nabla_{\theta_k}\ln\frac{d\mu}{d\mu_0}\bigl(f_{\theta_k}(u_i)\bigr).
$$

The experiments use Adam, with RMSprop and AdaGrad named as alternatives.

**Step five: four flow layers that satisfy the conditions.**

| Flow layer                       | Type      | Image dimension | Note                                                       |
| -------------------------------- | --------- | --------------- | ------------------------------------------------------------ |
| functional Householder flow      | linear    | one             | a reflection with parameters $v_n\in\mathcal H_u$ and $b_n\in\mathbb R$ |
| functional projected transformation flow | linear | $M$        | hence more expressive than Householder                     |
| functional planar flow           | nonlinear | —               | low-rank parameterisation, below                           |
| functional Sylvester flow        | nonlinear | —               | —                                                            |

Invertibility of the linear and nonlinear families is guaranteed by Lemma 2.1 and Lemma 2.2 respectively.

**Step six: why low rank is mandatory.** The functional planar flow has the form $f_n(u)=u+u_n h\bigl(\langle u,w_n\rangle_{\mathcal H_u}+b\bigr)$ with the low-rank parameterisation

$$
u_n=\sum_{i=1}^{r}\lambda_i\alpha_i\phi_i,
\qquad
w_n=\sum_{i=1}^{r}\lambda_i\beta_i\phi_i
$$

built from the eigenpairs $\{\lambda_i,\phi_i\}$ of the prior covariance $\mathcal C_0$, with $\{\alpha_i,\beta_i\}_{i=1}^{r}$ trainable. The paper stresses that low rank here is not a computational convenience but a **theoretical necessity**: without it, measure equivalence fails and the KL loss is no longer finite. That is the substantive difference from low-rank tricks in neural operators, which exist purely to save computation.

**Step seven: the conditional variant.** For a fixed forward model, different data induce different posteriors, so the unconditional version must be retrained per dataset. CNF-iVI conditions the flow on the data and handles measurement vectors of varying dimension at minimal extra cost.

### Theorems

**Theorem 2.3 (measure equivalence with an explicit Radon-Nikodym derivative).** The conclusion is as above. **Its hypothesis list was not transcribed verbatim from the verifiable material**: the conclusion and the existence of four verification propositions are confirmed, but not the conditions the theorem body imposes on each layer $f^{(n)}_{\theta_n}$. Quoting the theorem verbatim requires going back to the source.

**Lemmas 2.1 and 2.2.** Invertibility conditions for linear and nonlinear layers in infinite dimensions. Their existence and what they assert are confirmed; **their precise hypotheses are unverified**.

**Four propositions.** Each verifies that one of the concrete flows satisfies the conditions of Theorem 2.3.

**Discretisation-invariance proposition.** If $\mathcal H_u$ embeds continuously in $C(D)$, then for any $n$ the layers $\mathcal F^{(n)}_{\theta_n}:\mathcal H_u\to\mathcal H_u$ of all four flows are discretisation invariant.

The paper's own framing of the regularity issue is worth keeping as stated: without sufficient regularity, the $\mathcal H$-norm terms in the Radon-Nikodym derivative "are likely to be infinite", and the equivalence requirement of Theorem 2.3 is exactly the condition that forces enough regularity.

### Numerical experiments

| Item                 | Content                                                                             |
| -------------------- | ------------------------------------------------------------------------------------- |
| inverse problems     | a one-dimensional smooth equation, two-dimensional steady-state Darcy flow, electrical impedance tomography |
| baseline             | pCN (MCMC)                                                                          |
| verification targets | agreement with the theory, efficiency relative to pCN, empirical discretisation invariance |
| reference implementation | the author's repository contains a dedicated `commen_flows_dis_inv.py` for the discretisation-invariance test |

Discretisation invariance is tested by training the same flow at different discretisation levels and comparing behaviour, which is the direct experimental contrast with the discretise-then-infer route and the most distinctive piece of evidence in the paper.

> [!note] What is verifiable about the results
> **The specific error values could not be confirmed in the verifiable material**, so this page reports only the composition of the experiments and their verification targets, without a results table. Of the three inverse problems, electrical impedance tomography is the most strongly nonlinear and is also the one added in later versions.

**What the experiments establish and where they fall short.** They establish that the three verification targets are well chosen: agreement with theory, an efficiency comparison and discretisation invariance each correspond to a clear falsifiable claim, and discretisation invariance has independent implementation support. What is missing is that no numbers could be verified here, so "more efficient than pCN" can only be recorded as the paper's reported conclusion; and the paper has no formal publication record yet.

> [!note] Version drift
> Preprint v1 reports two inverse problems (the smooth equation and steady-state Darcy flow); v2 and v3 add electrical impedance tomography. This page follows v3. The paper is listed as submitted to _Math. Comput._ and **no published DOI was found**.

### Relation to the others

The one paper in the topic that works in function space by construction rather than by discretise-then-Bayesianise, and the only one whose central difficulty is measure-theoretic rather than computational. It shares the infinite-dimensional Bayes formulation with paper [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|79]] (both use $d\mu/d\mu_0\propto\exp(-\Phi)$ with a Gaussian prior) and the variational-rather-than-sampling stance with paper 106; but 106 does its flow modelling in a finite-dimensional latent space obtained by VAE-style dimension reduction, whereas this paper refuses to reduce and instead constrains the flow so that it stays measure-equivalent. The two make a natural pair for a discussion of discretise-then-Bayesianise against Bayesianise-then-discretise.

## 106: a latent-variable flow with an adaptive prior

### The idea

Three failure modes appear together in high-dimensional PDE-governed inversion, and this paper's judgement is that they have to be solved together.

**First, the posterior is non-Gaussian and often multimodal.** Gaussian-approximation samplers (unscented and ensemble Kalman inversion) answer such a posterior by averaging its modes away, and pCN strains in high dimension too. The posterior approximation therefore needs a more expressive family.

**Second, the surrogate falls out of distribution.** A neural operator pre-trained on prior samples is queried in the wrong region once the posterior concentrates elsewhere. This is the same diagnosis as on the [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|previous page]].

**Third, the prior mean is misspecified.** Every earlier paper handles this by hand, and hand-tuning the prior is exactly the intervention one wants to remove.

Standard normalizing flows help with the first but not with dimensionality, because a flow is a **bijection** and necessarily preserves dimension. The paper's answer combines VAE-style nonlinear dimension reduction with flows: the VAE supplies a non-bijective map from a $k$-dimensional latent space to the $d$-dimensional parameter space, but the VAE's own two components are too rigid — a standard normal latent prior and a diagonal-Gaussian encoder — so a flow is attached at each. **That is the whole content of the name "Variational Flow"**, and the strict evidence-lower-bound theorem is the statement that neither replacement can hurt and, under its conditions, each strictly helps.

One design choice in the prior-mean update deserves separate mention: **the covariance is deliberately not moved.** If mean and covariance both contracted toward the current posterior estimate, then one bad step would lock the effective prior onto the wrong region with no way back — mode collapse. Holding $\Sigma_0$ fixed keeps the search width constant while only the centre moves, so a wrong move remains recoverable.

Perturbing samples before fine-tuning the surrogate is motivated the same way: posterior samples are **concentrated**, and fine-tuning on them alone yields a surrogate accurate on a nearly null set and worse elsewhere, while the posterior is still moving. The perturbation $\gamma\nu^{(j)}$ with $\gamma>1$ spreads the training set back over a neighbourhood with width. It is applied in Karhunen-Loève coefficient space, where the prior is standard normal, so $\gamma$ is measured in prior standard deviations — which gives that hyperparameter a scale-free meaning.

### Setting

Parameter $\xi\in\mathbb R^{d}$ (Karhunen-Loève coefficients), latent $z\in\mathbb R^{k}$ with $k<d$, data $y$, potential $\Phi(\xi,y)$, prior $\pi_0$. The surrogate is a Fourier neural operator $\mathcal F_\vartheta$ pre-trained on prior samples. The unnormalised target is the posterior

$$
\hat p(x)=\exp\bigl(-\Phi(\xi,y)\bigr)\pi_0(\xi),
$$

where the flow's data variable $x$ is identified with the unknown parameter vector $\xi$.

### Derivation

**Step one: the sampling map.**

$$
z=f^{-1}_{\mathrm{pr},\beta}(v),\ v\sim\mathcal N(0,I)
\quad\Longrightarrow\quad
\xi=\mu_{\mathrm{de},\theta}(z)+\sigma_{\mathrm{de},\theta}(z)\odot\epsilon,\ \epsilon\sim\mathcal N(0,I).
$$

The latent prior is generated from a standard normal by a flow $f_{\mathrm{pr},\beta}$, and the encoder $q_{z|x,\alpha}$ is a conditional normalizing flow rather than a diagonal Gaussian. The loss $\mathcal L_{\mathrm{VF}}$ is a KL divergence between the model's joint distribution and a target joint constructed with the encoder (its displayed form is unverified).

**Step two: the evidence-lower-bound decomposition, which is the tool the theorem uses.** Maximising the bound over the encoder is equivalent to minimising

$$
\mathbb E_{p_x}\bigl[D_{\mathrm{KL}}\bigl(q_{z|x,\alpha}\,\|\,p_{z|x,\theta^{*},\tilde\beta}\bigr)\bigr],
\qquad\text{that is}\qquad
\mathrm{ELBO}=\text{const}-\mathbb E_{p_x}D_{\mathrm{KL}}\bigl(q_{z|x,\alpha}\,\|\,p_{z|x,\theta^{*},\tilde\beta}\bigr).
$$

**Step three: momentum updating of the prior mean.**

$$
\mu^{(k,i)}_{\mathrm{prior}}=\alpha\,\mu^{(k,i)}_{\mathrm{post}}+(1-\alpha)\,\mu^{(k-1)}_{\mathrm{prior}},
\qquad
\mu^{(k,i)}_{\mathrm{post}}=\frac1M\sum_{j=1}^{M}\xi^{(j)},
\quad \xi^{(j)}\sim p^{(k,i-1)}_{\mathrm{VF}}(\xi).
$$

Two time scales operate: $\mu_{\mathrm{post}}$ is re-estimated at every epoch $i$ of stage $k$, while the anchor $\mu^{(k-1)}_{\mathrm{prior}}$ is fixed within a stage and updated only at its end as $\mu^{(k)}_{\mathrm{prior}}\leftarrow\mu^{(k,N_e)}_{\mathrm{prior}}$. The effective prior is $\pi^{(k,i)}_0(\xi)=\mathcal N(\mu^{(k,i)}_{\mathrm{prior}},\Sigma_0)$, with the covariance held at $\Sigma_0$ to prevent mode collapse. Here $\alpha\in(0,1]$ controls update strength: smaller means stronger regularisation, $\alpha=1$ removes it, and prior work cited in the paper finds $\alpha\approx0.5$ close to optimal. Initialisation is $\mu_0=0$, $\Sigma_0=I$.

**Step four: perturbation before fine-tuning.**

$$
\hat\xi^{(j)}=\xi^{(j)}_{\mathrm{post}}+\gamma\,\nu^{(j)},
\qquad \nu^{(j)}\sim\mathcal N(0,I),\ \gamma>1 .
$$

**Step five: aggressive data replacement.** The paper contrasts itself with paper [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|79]] explicitly: instead of accumulating a growing historical dataset or running expensive greedy filtering, it **discards** the old data and fine-tunes the Fourier neural operator only on the newly generated local set $\mathcal D_k$.

**Step six: the stopping criterion is the relative change in the data misfit at the prior mean.**

$$
\frac{\bigl|\Phi\bigl(\mu^{(k-1)}_{\mathrm{prior}},y\bigr)-\Phi\bigl(\mu^{(k)}_{\mathrm{prior}},y\bigr)\bigr|}
{\Phi\bigl(\mu^{(k-1)}_{\mathrm{prior}},y\bigr)}\ <\ \epsilon .
$$

This is the same class of goal-oriented quantity as $e_D$ in paper 79; only the test point differs, an anchor selected by the true model from a candidate pool there against the evolving prior mean here. **The cost difference is not small**: the anchor of paper 79 requires running the true model on $T$ samples, while this needs $\Phi$ at a single point.

**Step seven: the order of the closed loop.** For $k=1,\dots,K$: run $N_e$ inner epochs, each re-estimating the posterior mean, updating the prior mean, setting the effective prior and updating the variational flow by minimising $\mathcal L_{\mathrm{VF}}$ against that prior; after the inner loop, update the stage anchor and test the stopping criterion; if not stopped, draw $M$ posterior samples from the flow, perturb them, solve the PDE for each perturbed sample, **replace** the dataset with the new one and fine-tune the Fourier neural operator. The output is the trained flow and operator, posterior samples, and the posterior mean.

**Evaluation metrics.** The relative inversion error $e_{\mathcal I}$ is computed from the flow's posterior mean; the relative **surrogate fitting error** $e_{\mathcal S}$ is computed on a local test set of $N=100$ samples $\xi^{(i)}=\xi_{\mathrm{ref}}+\eta^{(i)}$ with $\eta^{(i)}\sim\mathcal N(0,I)$, as the average relative $L^{2}$ discrepancy against the finite-difference state field. The second metric is worth noting: it measures surrogate quality near the truth directly, so "the surrogate got better" can be separated from "the inversion got more accurate".

### Theorems

**Strict evidence-lower-bound improvement.** Under either of two conditions (i) and (ii), the Variational Flow model attains a **strictly higher** evidence lower bound than a standard VAE; when both hold, the flow prior and the conditional flow encoder each contribute a strictly positive improvement. The proof runs in two steps: fix the optimal VAE encoder and decoder and replace the latent prior by a flow prior to obtain improvement (i); then fix the decoder and flow prior and extend the encoder from a diagonal Gaussian to a conditional flow to obtain improvement (ii), using the decomposition of step two above; if only (ii) holds, take $\tilde\beta$ to be the identity map.

> [!warning] The theorem's conditions were not transcribed
> **The precise statements of conditions (i) and (ii) could not be confirmed in the verifiable material.** What is confirmed is the form of the conclusion and the two-step proof structure, not the conditions themselves. Quoting the theorem verbatim requires going back to the source.

**No convergence guarantee.** The paper states explicitly that future work will address the current lack of theoretical convergence guarantees. The adaptive loop therefore has no convergence theorem, consistent with every paper in this topic except 79.

### Numerical experiments

Four problems:

| Problem                        | Configuration                                                              |
| ------------------------------ | ---------------------------------------------------------------------------- |
| 100-dimensional Rosenbrock     | tests posterior approximation quality alone, so prior updating is **deliberately switched off** |
| one-dimensional Darcy flow     | $\Omega=[0,1]$                                                              |
| two-dimensional Darcy flow     | $[0,1]^2$                                                                   |
| two-dimensional Navier-Stokes  | recover the initial vorticity $\omega_0$ from the vorticity field $\omega(\cdot,T)$ at $T=1$ |

Shared settings for the three PDE problems:

| Item             | Setting                                                            |
| ---------------- | -------------------------------------------------------------------- |
| truncation       | $d\in\{32,64\}$                                                     |
| noise levels     | $\delta\in\{1\%,5\%,10\%\}$                                         |
| repetitions      | averaged over 3 runs                                                |
| baselines        | pCN, SVGD, UKI-FDM (full order), UKI-FNO (operator surrogate)      |
| metrics          | relative inversion error $e_{\mathcal I}$, relative surrogate fitting error $e_{\mathcal S}$ |

The Rosenbrock problem additionally compares two-dimensional marginals on the $\xi_1$–$\xi_2$ plane against a known ground-truth density, with a VAE, MCMC, SVGD and unscented Kalman inversion.

**Reported outcomes.**

- Rosenbrock: the Variational Flow is the most accurate and significantly better than a vanilla VAE, which the authors read as empirical support for the bound theorem. Switching off prior updating is deliberate: what wins can then only be the generative model itself.
- One-dimensional Darcy: better than every baseline, most clearly at high noise.
- Two-dimensional Darcy: best in most medium- and high-noise cases and competitive at low noise. When the reference field is drawn out of distribution from a uniform law, **prior updating is what lets it find the true posterior**, while pCN and SVGD, lacking that module, stay near the initial $\mathcal N(0,I)$ prior and get the scale of the field visibly wrong.
- Navier-Stokes: the lowest inversion error at **every** tested truncation dimension and noise level.
- Across all three PDE problems, the surrogate fitting error $e_{\mathcal S}$ converges to lower values in fewer adaptive stages than the baselines.

> [!note] Quantitative tables
> **Every numerical table in this paper is illegible in the available material**, so what appears above is configuration and qualitative conclusions only, and this page prints no error values from it. All the "best" and "lowest" statements are the paper's reported conclusions and cannot be checked quantitatively here.

**What the experiments establish and where they fall short.** They establish modular evidence: Rosenbrock tests the generative model in isolation (prior updating off), the out-of-distribution two-dimensional Darcy case tests prior updating in isolation (against pCN and SVGD, which lack the module), and $e_{\mathcal S}$ tests the surrogate fine-tuning in isolation. Splitting the three modules apart is more convincing than a single aggregate accuracy comparison. What is missing: no numbers could be verified here, so the margin behind "competitive or superior" is unknown; the paper concedes it has no convergence guarantee; the update strength $\alpha$ is taken as the near-optimal value reported in cited prior work with no sensitivity study of its own; and the superiority of aggressive replacement over greedy filtering is asserted as a position rather than shown by a direct ablation.

### Relation to the others

The most direct descendant of paper [[en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|79]]: the same closed loop (approximate posterior, select points, fine-tune the surrogate, repeat), the same goal-oriented relative-change-in-$\Phi$ stopping rule, and the same out-of-distribution diagnosis of why fixed operator surrogates fail. What changes: unscented Kalman inversion becomes a variational flow, so the posterior approximation is no longer Gaussian — exactly paper 79's acknowledged limitation; DeepONet becomes a Fourier neural operator; greedy diversity selection becomes perturb-and-replace; and prior-mean adaptation is added, which none of the earlier papers has. It also uses paper 56's SVGD as a **baseline** rather than a tool. The dimension-reduction thread connects to the earlier deep adaptive sampling work of Kejun Tang and Xiaoliang Wan.

## Side-by-side comparison

| No. | Posterior approximation                  | Surrogate                        | Refinement                   | Theory                                   | Numerical results                  |
| --- | ---------------------------------------- | -------------------------------- | ---------------------------- | ---------------------------------------- | ---------------------------------- |
| 55  | optimisation-based independence proposal | feed-forward network             | none (offline once)          | none in the paper itself                 | setup and qualitative only         |
| 56  | particle flow                            | feed-forward network             | online, particle mean        | none                                     | one example with complete numbers  |
| 82  | ensemble Kalman filter                   | operator-inference reduced model | none (error learned offline) | none                                     | problem list and qualitative only  |
| 88  | staged hierarchical Bayesian inference   | eight feed-forward networks      | none                         | none                                     | configuration and anchors          |
| 99  | function-space variational flow          | no surrogate                     | not applicable               | measure equivalence, explicit derivative | setup and verification targets only |
| 106 | latent-variable variational flow         | Fourier neural operator          | online, replacement          | strict evidence-bound improvement        | setup and qualitative only         |

## Coverage check

| Item                                              | Paper    | Status                                                                                   |
| ------------------------------------------------- | -------- | ------------------------------------------------------------------------------------------ |
| Verification tier and remaining gaps              | all      | per-paper source, tier and unverified items                                              |
| Geometric reading of RTO                          | 55       | residual surface, tangent projection, exactness in the linear case, the two effects of curvature matching the two factors of $c$ |
| RTO proposal, density, weight and scalable variant | 55      | linearisation point, proposal, density, weight, rank-$r$ split and determinant simplification |
| Convention clash and normalisation typo           | 55       | $w$ and $c$ as reciprocals with reversed numerators; the correction of $(2\pi)^{-\pi/2}$ |
| Bardsley's theorem and the ergodicity condition   | 55       | full hypotheses, conclusion, the positive lower bound on $c$ and why it fails for the surrogate version |
| Stein operator and the closed-form direction      | 56       | functional optimisation, identity, absence of the normalising constant, closed form, attraction and repulsion |
| Point selection with a separation constraint      | 56       | division of labour between objective and constraint, degeneracy motivation, meaning of radius shrinking |
| Relation to the CMPS trust-region bounds          | 56       | the $R^{2}$ and $R^{3}$ bounds with poisedness, and why they do not transfer to a network |
| Quantitative table for example 3 and its limits   | 56       | errors and online solve counts, undefined $n_t$, damaged headers, self-inconsistent vanilla-SVGD cell |
| Two-tier reduction plus an error network          | 82       | coarse-grid operator inference, error network, output-side correction, with the verifiability limit stated |
| Staged hierarchical architecture and annealing    | 88       | the two structural difficulties matched to two remedies, parameters and outputs, annealing, eight surrogates, experimental anchors |
| Measure equivalence, the derivative and four flows | 99      | the singularity obstruction, theorem conclusion, change of variables onto the prior, necessity of low rank, discretisation invariance |
| Latent flow, prior updating, aggressive replacement | 106    | sampling map, evidence-bound decomposition, momentum update with fixed covariance, perturbation, stopping rule |
| Record of two untranscribed theorem hypothesis lists | 99, 106 | the condition list of Theorem 2.3 and conditions (i), (ii) of the evidence-bound theorem |

## Sources for this page

- L. Yan and T. Zhou, [_An acceleration strategy for randomize-then-optimize sampling via deep neural networks_](https://doi.org/10.4208/jcm.2102-m2020-0339), J. Comput. Math. 39(6) (2021), pp. 848-864 (preprint [arXiv:2104.06285](https://arxiv.org/abs/2104.06285)).
- L. Yan and T. Zhou, [_Stein variational gradient descent with local approximations_](https://doi.org/10.1016/j.cma.2021.114087), Comput. Methods Appl. Mech. Engrg. 386 (2021), 114087 (preprint [arXiv:2104.06276](https://arxiv.org/abs/2104.06276)).
- Y. Wang, L. Yan, and T. Zhou, [_Deep learning-enhanced reduced-order ensemble Kalman filter for efficient Bayesian data assimilation of parametric PDEs_](https://doi.org/10.1016/j.cpc.2025.109544), Comput. Phys. Commun. 311 (2025), 109544.
- S. Wang, L. Ma, L. Guo, X. Li, and T. Zhou, [_Multi-stage uncertainty quantification framework for red blood cell morpho-mechanics_](https://doi.org/10.1016/j.ijmecsci.2026.111352), Int. J. Mech. Sci. 313 (2026), 111352 (preprint titled _An RBC-MsUQ Framework for Red Blood Cell Morpho-Mechanics_, [arXiv:2508.06852](https://arxiv.org/abs/2508.06852)).
- Y. Zhao, H. Lu, J. Jia, and T. Zhou, _Functional normalizing flow for statistical inverse problems of partial differential equations_, [arXiv:2411.13277](https://arxiv.org/abs/2411.13277), submitted to Math. Comput.; reference implementation [jjx323/FunctionalNormalizingFlow](https://github.com/jjx323/FunctionalNormalizingFlow).
- Y. Wang, X. Wang, K. Tang, X. Wan, T. Zhou, and C. Yang, _Deep adaptive dimension reduction for Bayesian inference in inverse problems_, [arXiv:2605.29373](https://arxiv.org/abs/2605.29373), submitted to SIAM J. Sci. Comput.
- Background sources: J. M. Bardsley, A. Solonen, H. Haario, and M. Laine, [_Randomize-then-optimize: a method for sampling from posterior distributions in nonlinear inverse problems_](https://doi.org/10.1137/140964023), SIAM J. Sci. Comput. 36(4) (2014), pp. A1895-A1910; Q. Liu and D. Wang, _Stein variational gradient descent: a general purpose Bayesian inference algorithm_, NeurIPS 2016 ([arXiv:1608.04471](https://arxiv.org/abs/1608.04471)); P. R. Conrad, Y. M. Marzouk, N. S. Pillai, and A. Smith, [_Accelerating asymptotically exact MCMC for computationally intensive models via local approximations_](https://doi.org/10.1080/01621459.2015.1096787), J. Amer. Statist. Assoc. 111(516) (2016), pp. 1591-1607.
