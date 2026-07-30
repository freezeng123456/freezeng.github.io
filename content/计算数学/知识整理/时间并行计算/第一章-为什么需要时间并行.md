---
title: "Chapter 1: Why Parallelize in Time?"
description: Causality, all-at-once systems, and criteria for useful parallelism
tags:
  - parallel-in-time
  - computational-mathematics
---

## The limit of spatial parallelism

Domain decomposition and mesh parallelism are usually the first sources of concurrency in a PDE solver. Once each spatial subdomain becomes too small, communication and synchronization dominate additional local computation. The temporal grid remains large, but conventional time stepping is sequential:

$$
u_{n+1}=\Phi_{\Delta t}(u_n).
$$

Parallel-in-time methods do not eliminate causal dependence. They introduce redundant work, global iteration, decomposition, or a transform so that multiple time intervals can be processed concurrently.

## The all-at-once formulation

For a linear one-step discretization,

$$
u_{n+1}-\Phi u_n=g_{n+1}.
$$

Stacking all temporal unknowns gives

$$
\begin{bmatrix}
I \\
-\Phi&I\\
&-\Phi&I\\
&&\ddots&\ddots
\end{bmatrix}
\begin{bmatrix}u_0\\u_1\\u_2\\\vdots\end{bmatrix}
=
\begin{bmatrix}g_0\\g_1\\g_2\\\vdots\end{bmatrix}.
$$

Sequential time stepping is exact forward substitution. The principal PinT families replace it by different approximations to the inverse:

- Parareal and MGRiT use inexpensive coarse propagation.
- STMG uses a multilevel approximation on the complete space-time grid.
- ParaDiag replaces the triangular temporal coupling by a matrix that can be diagonalized by FFT.
- SWR solves space-time subdomain problems and exchanges interface waveforms.

## Why hyperbolic problems are more difficult

For a dissipative propagator, many error modes satisfy $|\lambda(\Phi)|<1$. The rapidly damped components do not need to be represented accurately on the coarsest temporal level. For wave and transport problems, relevant modes remain close to the unit circle. Their error is primarily a phase error rather than an amplitude error. A small phase mismatch between fine and coarse propagation can accumulate over many time intervals.

This mechanism explains why:

1. decreasing the viscosity $\nu$ makes advection-diffusion increasingly difficult for standard temporal coarsening;
2. Parareal may exhibit transient error growth in transport-dominated regimes;
3. characteristic transmission, phase correction, waveform relaxation, or global diagonalization can be more effective than a merely more dissipative coarse solver.

## Evaluation criteria

Convergence histories are necessary but insufficient. A computational study should report:

| Quantity             | Question                                                                     |
| -------------------- | ---------------------------------------------------------------------------- |
| convergence          | How many global iterations reach the discretization-error scale?             |
| coarse-to-fine cost  | Is the coarse propagator substantially cheaper than the fine propagator?     |
| temporal concurrency | How many fine or local solves can execute simultaneously?                    |
| communication        | Does each iteration require neighbor exchange or a global collective?        |
| memory               | Must all temporal states be stored?                                          |
| parameter robustness | Does convergence degrade toward the hyperbolic limit?                        |
| scalability          | Does the iteration count remain controlled as the number of intervals grows? |

For Parareal with $N$ time intervals, a simplified idealized speedup bound is

$$
S\lesssim \frac{NC_F}{K(C_F+NC_G)},
$$

when communication is neglected. Useful acceleration requires both $K\ll N$ and $C_G\ll C_F$. Algorithmic convergence alone does not imply either condition.

## Reproduction policy used in these notes

The companion project has two execution levels:

```bash
python3 run_experiments.py all --quick --output-dir results/quick
python3 run_experiments.py all --output-dir results/formal
```

The quick configuration preserves the algorithmic structure with smaller grids and serves only as a smoke test. Every figure and numerical value in Chapters 2-4 comes from the formal configuration. Each experiment writes a PNG figure and a JSON record containing its parameters and metrics.

## Summary

PinT is best viewed as a structured approximation to the inverse of a global temporal system. The first diagnostic question is whether the dynamics are dominated by dissipation or propagation. The second is which physical information the parallel approximation preserves: amplitude, phase, characteristics, or slowly decaying modes.
