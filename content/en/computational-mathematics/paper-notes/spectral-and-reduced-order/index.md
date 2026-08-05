---
title: Spectral Methods and Reduced-Order Representations
description: Six papers on one question - align the structure of the representation with the structure of the solution
lang: en
translation: computational-mathematics/paper-notes/spectral-and-reduced-order
tags:
  - computational-mathematics
  - paper-notes
  - spectral-methods
---

This topic holds only 6 papers but spans two apparently unrelated technical lines: spectral bases on unbounded domains, and low-rank representations for random problems. They belong together because both answer one question: **the structure of the representation should be dictated by the structure of the solution, not by algorithmic convenience.**

On an unbounded domain the decisive structure is the **far-field decay rate**. For random low-rank problems it is the **rank of the solution manifold and its smallest singular value**. Both settings show the same consequence of a mismatch: spectral accuracy or the error constant degrades badly, and adding degrees of freedom does not recover it.

![Match the decay of the basis to the decay of the solution](assets/diagrams/tao-zhou-papers/en/unbounded-spectral.svg)

## Two close readings

| Close reading                                                                                                                                        | Papers     | Decisive structure                     |
| ---------------------------------------------------------------------------------------------------------------------------------------------------- | ---------- | -------------------------------------- |
| [[en/computational-mathematics/paper-notes/spectral-and-reduced-order/unbounded-domain-spectral\|Unbounded domains and fractional operators]]        | 27, 42, 92 | far-field decay and the operator image |
| [[en/computational-mathematics/paper-notes/spectral-and-reduced-order/low-rank-and-variable-separation\|Low-rank manifolds and variable separation]] | 15, 17, 83 | rank and the smallest singular value   |

## Core idea of each paper

- **27 (Hermite spectral collocation)** collocates directly in physical space on an unbounded domain, using explicit differentiation matrices for the fractional Laplacian applied to Hermite functions. The contemporaneous alternative collocates in frequency space, so every evaluation costs a forward and a backward Hermite transform, which is awkward for nonlinear problems. The paper supplies both the normalised Hermite functions and the over-scaled (Brinkman) basis, the latter already standard for Fokker-Planck equations.
- **42 (rational spectral methods)** names the core issue: solutions of fractional operators decay only **algebraically** at infinity, while Hermite and Laguerre bases are tuned to exponential decay. The paper uses the singular algebraic map $x=t/\sqrt{1-t^2}$ to pull the real line onto $(-1,1)$ and builds modified mapped Gegenbauer functions decaying like $|x|^{-(\lambda+1)}$, so $\lambda$ becomes a knob tuned to the solution's power law. The technical obstruction is structural: the ordinary Laplacian maps a rational basis function to another of the same type, while the fractional Laplacian maps it to a class of a completely different nature. The paper resolves this by writing each basis function as a finite combination of two elementary rational shapes and giving the fractional Laplacian of those shapes in closed hypergeometric form.
- **92 (fast computation of a fractional matrix power)** closes a **uniformity** gap rather than a raw speed gap: existing methods do not work well as $\alpha$ approaches $0$ or $1$. The structural reason is clear — the standard integral representations carry a prefactor $\sin(\alpha\pi)$ that vanishes at integer $\alpha$, and the integrand's endpoint behaviour degenerates with $\alpha$; the error asymptotics of best uniform rational approximation carry the same factor. The paper targets cost and accuracy uniform in $\alpha\in(0,1)$ for both $A^{-\alpha}$ and the resolvent $(q\mathcal I+A^{\alpha})^{-1}$ needed for implicit time steps.
- **15 (error analysis of the dynamically orthogonal approximation)** addresses the fact that the dependence of a time-dependent random PDE solution on its random parameters changes substantially as time advances, so a fixed basis (generalised polynomial chaos or proper orthogonal decomposition modes) needs a steadily growing number of terms. The dynamically orthogonal approximation evolves the spatial basis and the stochastic coefficients together, avoiding that growth, but previously had essentially no error theory. The paper's first goal is to establish a precise correspondence with the dynamical low-rank approximation so the existing theory can be imported.
- **17 (multilevel finite elements for Fredholm integral eigenvalue problems)** targets the Karhunen-Loève expansion, which needs the eigenpairs of the Fredholm integral operator built from a covariance kernel. Discretisation gives a **dense** matrix eigenvalue problem, with no sparsity to exploit unlike a differential operator. The paper replaces the fine-mesh eigensolve by a sequence of integral iterations plus eigenvalue solves on the coarsest mesh only, bringing total work down to roughly one integration step on the finest mesh.
- **83 (dynamical variable separation for parameter-dependent dynamical systems)** belongs to the same family as paper 15 but inverts its structure: instead of evolving all modes of a rank-$S$ representation together, it enriches the reduced basis one term at a time by a greedy algorithm, so each step splits into two decoupled evolution equations, a parameter-independent PDE and a parameter-dependent ODE. No covariance has to be inverted, and the computation divides into an offline basis construction and an online low-rank evaluation.

> [!note] Coverage status
> Papers **15, 27, 42 and 83** have close-reading content checked equation by equation against publicly available full texts, including the hypotheses of their theorems and the setups of their numerical experiments. The text of paper **17** sits behind the Elsevier paywall, with no preprint on arXiv or on either author's page, so its section reports only what the abstract and keywords confirm and keeps the Lin-Xie framework it adapts explicitly separate from the paper's own claims; the actual form of its correction step for the integral operator, the constants and hypotheses of its error estimates, and all of its numerical results are unverified here. The full text of paper **92** could not be reached at all (the SIAM site is blocked by a Cloudflare challenge, and there is no arXiv preprint or open-access copy), so this site reports only what its abstract, keywords and reference list confirm.

## Decisive structure one: the far-field decay rate

### Match the tail of the basis to the tail of the solution

A spectral method on a bounded domain has to cope with boundaries; on an unbounded domain there are no boundaries to cope with — **there are tails**. That single fact drives everything: spectral accuracy comes from the basis reproducing the solution's behaviour at infinity with a few terms, and if the basis decays in a different manner from the solution, then however many terms are added, the expansion spends them repairing the tail and convergence drops from spectral to algebraic.

Papers 27 and 42 stand on opposite sides of that criterion, having chosen two different kinds of tail.

**Paper 27 takes the exponential side.** Hermite functions carry a Gaussian factor $e^{-x^2/2}$, so their tails are exponential; on that basis the paper claims spectral convergence for solutions that decay exponentially at infinity, and supplies both the normalised Hermite functions and the over-scaled (Brinkman) basis so that the scaling factor becomes a knob — a scaling factor is essential when decay rates are mismatched. The limit of that route is visible in the paper's own experiments: one eigenvalue example shows only **algebraic** convergence, which the paper attributes to the eigenfunctions of that problem decaying algebraically rather than exponentially.

**Paper 42 takes the algebraic side, making that limit its starting point.** It observes that solutions of fractional operators follow a power law at infinity and decay only algebraically, so Hermite and Laguerre bases, designed for exponential decay, are mismatched by construction; the same reasoning rules out domain truncation as a shortcut. Its replacement is a family of modified mapped Gegenbauer functions, obtained by pulling the real line onto $(-1,1)$, decaying like $|x|^{-(\lambda+1)}$ — so **the decay rate itself becomes a parameter $\lambda$** that can be tuned directly to the solution's power law.

That is why two very different bases appear inside one topic: not a matter of taste, but two papers aiming at two kinds of tail. Paper 92 does not sit on this axis — it works on the algebraic side of the same fractional operators, computing $A^{-\alpha}$ and $(q\mathcal I+A^{\alpha})^{-1}$ uniformly in $\alpha\in(0,1)$.

### The image of the basis decides whether the method is feasible

In papers 27 and 42 the real obstacle is not the basis but **the image of the basis under the operator**. If the image stays in the same family — Jacobi poly-fractonomials or generalised Jacobi functions on bounded domains — the fractional operator becomes effectively local. If it leaves the family, the basis must be decomposed into a few elementary shapes whose images are known in closed form. Paper 42 takes exactly that route, and that is how it makes concrete its statement that extending the mapping technique to the fractional setting is far from trivial.

## Decisive structure two: the low-rank manifold

### A reduction, not a truncation

The reason papers 15, 17 and 83 belong to one topic is sharper than "they all involve low rank": all three deal with a single object — **the manifold $\mathcal M_S$ of rank-$S$ fields** — approached from three different directions.

The distinction that matters is between reduction and truncation. Truncation fixes the basis in advance (polynomial chaos in the random variables, or proper orthogonal decomposition modes computed from snapshots), keeps $S$ terms, and leaves as unknowns just $S$ coefficients in a **fixed linear subspace**; accuracy is governed by how fast the tail decays in that fixed basis, and once the solution rotates out of the subspace the only repair is more terms. The low-rank manifold approach instead treats $\mathcal M_S$ itself as the state space. It is not a linear subspace but a curved manifold, the unknown is a point moving on it, and both factors of the representation are unknown. The dynamics come from a Galerkin (Dirac-Frenkel) condition requiring the residual to be orthogonal to the tangent space at the current point,

$$
\frac{\partial u_S}{\partial t}=P_{u_S}\bigl(\mathcal L(u_S)\bigr).
$$

Two things follow that truncation cannot deliver. First, the subspace occupied by the solution can rotate and deform **at fixed rank**, because the basis is itself evolving, so no terms need to be added merely to chase a moving subspace — that is how the Kolmogorov barrier which proper orthogonal decomposition meets on transport and hyperbolic problems is dodged. Second, the equations are written entirely in the factors, so cost scales with $S$ and with one deterministic field solve rather than with the dimension of the ambient product space, and the mean and total variance are read straight off the factors without sampling. **That is the substantive sense of "reduction": the ambient dimension never appears in the cost of the evolution.**

The three papers occupy three corners of that picture. Paper **15** supplies the error theory of the flow: the dynamically orthogonal solution is quasi-optimal against the pointwise-in-time best rank-$S$ approximation, but with a constant of the form $e^{C/\rho}$. Paper **17** deals with that "best rank-$S$ approximation" itself — it is the truncated Karhunen-Loève expansion, defined by the Fredholm integral eigenvalue problem for the covariance kernel, which discretises to a **dense** matrix eigenproblem; even on the manifold route it has to be solved once at $t=0$, because the dynamically orthogonal flow is ignited from the Karhunen-Loève expansion of the initial datum. Paper **83** builds a trajectory on the manifold a different way: greedy enrichment one term at a time instead of evolving all modes together, so the coefficient update divides by a scalar rather than a matrix.

### Rank is not a free parameter, and the smallest singular value is a risk

Low-rank evolution confines the solution to a rank-$R$ manifold. The error has two parts: the distance from the true solution to the best rank-$R$ field, and the drift introduced by the tangent-space projection. The second is tied to the curvature of the manifold, which grows as the smallest singular value approaches zero. Choosing the rank is therefore not only an accuracy-cost trade-off but a stability decision.

Paper 15 makes that judgement quantitative: the constants in its curvature lemma are multiples of $\rho^{-1}$, with $\rho$ a lower bound on the smallest singular value, and they pass untouched into the quasi-optimality theorem, where the error constant grows like $e^{C/\rho}$. In a computation that same $\rho$ is a lower bound on $\sqrt{\lambda_{\min}(C(t))}$, with $C$ the covariance matrix of the stochastic coefficients. Paper 83 is designed squarely against that inverse — building modes one at a time reduces the divisor to a scalar, so $C^{-1}$ never appears; the price is that earlier modes are never revisited, so each added term reduces the error more slowly. **Both routes trade the same pair of goods: a guarantee of quasi-optimality against the stability of never inverting anything.**

![Evolve a random solution directly on a low-rank manifold](assets/diagrams/tao-zhou-papers/en/low-rank-dynamics.svg)

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]; per-paper references appear at the end of each close-reading page.
