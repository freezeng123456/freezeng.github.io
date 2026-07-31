---
title: Rotationally Symmetric Stein Estimators
description: Bandwidth-stable high-order stochastic differentiation for smoothed PINNs
tags:
  - computational-mathematics
  - PINN
  - Monte-Carlo
---

Repository: [freezeng123456/RSSE](https://github.com/freezeng123456/RSSE)

## Problem

Randomized smoothing can replace a high-order derivative in a physics-informed neural network by an expectation involving function values and a Stein weight. This avoids the explicit construction of a high-order automatic-differentiation tensor. The Stein weight, however, contains a factor of order $t^{-p}$ for a derivative of order $p$. If lower-order Taylor terms are not cancelled exactly, their contribution is amplified as the smoothing bandwidth $t$ approaches zero.

The standard two-point antithetic construction is effective for several second-order operators, but it does not provide enough degrees of freedom to remove every Taylor mode below a general derivative of order three or higher.

## Main insight

Let $\omega=\exp(2\pi i/N)$. For one Gaussian direction $\xi$, evaluate an analytic network at the complex phases

$$
x+t\omega^j\xi,\qquad j=0,\ldots,N-1,
$$

and combine the values with weights $\omega^{-jp}$. The discrete rotational symmetry acts as a Taylor-mode filter: it preserves the target mode of order $p$ and cancels all lower-order modes.

The principal result has the form

$$
\mathbb E[\mathcal E_{\alpha,t,N}u]
=\partial^\alpha u+O(t^N),\qquad
\operatorname{Var}(\mathcal E_{\alpha,t,N}u)
=V_\alpha(u)+O(t^N),
$$

where the leading variance $V_\alpha(u)$ is independent of $t$. The distinction from an ordinary antithetic estimator is therefore not only a higher-order bias expansion. It is the removal of the bandwidth singularity from the leading variance term.

## Algorithmic structure

1. Fix a target multi-index $\alpha$, set $p=|\alpha|$, and choose $N\ge p$.
2. Draw $\xi\sim\mathcal N(0,I)$ for each Monte Carlo sample.
3. Evaluate the network at the $N$ rotationally symmetric complex phases.
4. Apply the roots-of-unity weights that extract the target Taylor mode.
5. Multiply by the multivariate Hermite weight and average over samples.
6. Insert the estimator into the PINN residual, together with an exact boundary parameterization or a Leibniz expansion when required.

## Selected results

RSSE was evaluated on triharmonic and biharmonic problems in dimensions $d=20,50,100$, and on a sixth-order nonlinear Cahn-Hilliard problem in dimensions $d=5,10,20$. Under the reported protocols, it attained lower relative $L^2$ errors than the tuned antithetic baseline in each of these settings.

For the biharmonic experiments, the mean relative errors of RSSE were approximately
$1.10\times10^{-2}$, $7.17\times10^{-3}$, and $4.29\times10^{-3}$ for the three dimensions. The corresponding baseline errors were approximately
$1.55\times10^{-2}$, $9.83\times10^{-3}$, and $7.35\times10^{-3}$.

![[assets/research/rsse-k4.png]]

## Scope

- This note does not disclose the complete weight construction, the proofs, or implementation-specific optimizations.
- The analysis assumes an analytic continuation of the network to a complex neighborhood. Implementations therefore use an entire activation function.
- Complex evaluation increases the cost of one forward pass. The observed benefit arises from variance control and a wider usable bandwidth range, not from a universal reduction in wall-clock time.
- The empirical comparison concerns the stated PDEs, architectures, and sampling budgets. It does not establish superiority for arbitrary differential operators or neural networks.

## Relation to the apolarity project

[[computational-mathematics/research/apolarity|Apolarity-Guided Taylor Jets]] seek the shortest deterministic directional schedule for one exact mixed derivative. RSSE instead removes low-order modes from a stochastic smoothing estimator. Both use roots of unity, complex directions, and Taylor coefficient extraction, but they solve different optimization problems.
