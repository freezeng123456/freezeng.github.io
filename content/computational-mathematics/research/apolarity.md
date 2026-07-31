---
title: Apolarity-Guided Taylor Jets
description: Waring-optimal directional schedules for prescribed mixed derivatives
tags:
  - computational-mathematics
  - automatic-differentiation
  - algebraic-geometry
---

Repository: [freezeng123456/apolarity](https://github.com/freezeng123456/apolarity)

## Problem

Nested automatic differentiation builds an increasingly deep derivative graph for a prescribed mixed derivative of order $p$. Taylor-mode automatic differentiation can propagate a truncated jet along one direction, but this observation leaves a combinatorial question: which directions recover the target mixed derivative exactly, and what is the minimum number of such directions?

## Main insight

The directional Taylor coefficient

$$
T_p(x;v)=\frac1{p!}\frac{d^p}{d\tau^p}u(x+\tau v)\bigg|_{\tau=0}
$$

is a homogeneous polynomial of degree $p$ in the direction $v$. Expressing a prescribed mixed derivative as a linear combination of $T_p(x;v_r)$ is equivalent to expressing the associated monomial $z^\alpha$ as a sum of $p$th powers of linear forms. This is a Waring decomposition.

If the active exponents are ordered as $a_0\le\cdots\le a_n$, the complex Waring rank gives the minimum number of directions:

$$
R_{\mathbb C}(z^\alpha)=\prod_{j=1}^{n}(a_j+1).
$$

A roots-of-unity construction attains this rank and produces an executable schedule of complex directions.

## Computational structure

```text
multi-index
  -> active exponents
  -> roots-of-unity Waring schedule
  -> directional Taylor jets
  -> weighted exact derivative
  -> custom reverse rule
  -> parameter-gradient training
```

The forward pass propagates Taylor coefficients through linear layers and the recurrence for $\sinh$ up to order $p$. A custom reverse rule differentiates the resulting derivative evaluation with respect to the network parameters.

## Selected results

- The computed derivative values and parameter gradients agree with nested automatic differentiation in the reported validation tests.
- For mixed derivatives with repeated indices, the complex Waring schedule can use fewer directions than real polarization. For square-free derivatives, both constructions require $2^{p-1}$ directions, so the complex construction has no directional-count advantage in that case.
- Under a common network width, a common 1200-second budget, and five random seeds, the complex-$\sinh$ model attained the lowest mean relative $L^2$ error in the 12 reported PDE configurations. The test set covers polyharmonic, radial-chirp, and lossy-Maxwell problems.

![[assets/research/apolarity-chirp.png]]

## Scope

The construction is an exact backend for one prescribed monomial derivative. It is not a universal acceleration for sums of operators such as the Laplacian or $\Delta^m$. A shorter algebraic schedule also does not imply a proportional wall-clock improvement, because complex arithmetic, batching, numerical precision, and memory layout affect the measured cost.

The PDE results provide evidence for the complex-$\sinh$ representation under the reported settings. They do not isolate the Waring schedule as the sole cause of the accuracy gain. Matched-parameter, runtime, and peak-memory comparisons remain necessary for that attribution.
