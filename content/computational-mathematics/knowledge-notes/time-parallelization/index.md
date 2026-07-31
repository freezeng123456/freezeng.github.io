---
title: Time Parallelization for Hyperbolic and Parabolic Problems
description: Chapter-by-chapter notes and reproducible experiments for parallel-in-time methods
tags:
  - computational-mathematics
  - parallel-in-time
---

These notes follow M. J. Gander, S.-L. Wu, and T. Zhou, _Time Parallelization for Hyperbolic and Parabolic Problems_, Acta Numerica 34 (2025), pp. 385-489. The survey distinguishes methods that remain effective for propagative problems from methods designed primarily for dissipative problems.

## Reading sequence

1. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-1-why-parallelize-in-time|Chapter 1: Why Parallelize in Time?]] introduces causality, the all-at-once formulation, and performance criteria.
2. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-2-model-problems|Chapter 2: Model Problems]] compares the heat, advection-diffusion, Burgers, and wave equations. The chapter includes all three recomputed solution experiments.
3. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-hyperbolic-methods|Chapter 3: Methods Effective for Hyperbolic Problems]] discusses SWR, PIDC/RIDC, ParaExp, and ParaDiag. The chapter includes both ParaDiag-II experiments.
4. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-4-parabolic-methods|Chapter 4: Methods Designed for Parabolic Problems]] covers Parareal, PFASST, MGRiT, diagonalization-based Parareal, and STMG. It contains the recomputed convergence studies.
5. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view|Chapter 5: A Unified View and Method Selection]] compares the algorithms and records the complete experiment ledger and reproduction protocol.

## Method map

| Method    | Parallel unit                             | Mechanism                              | Natural regime                 |
| --------- | ----------------------------------------- | -------------------------------------- | ------------------------------ |
| SWR       | overlapping space-time subdomains         | waveform transmission                  | transport and wave problems    |
| PIDC/RIDC | correction levels and time nodes          | deferred correction and pipelining     | initial-value problems         |
| ParaExp   | inhomogeneous and homogeneous subproblems | exponential propagation                | linear systems                 |
| ParaDiag  | all-at-once temporal matrix               | circulant approximation and FFT        | linear or linearized systems   |
| Parareal  | coarse time intervals                     | coarse prediction plus fine correction | moderate to strong dissipation |
| PFASST    | collocation nodes across time steps       | SDC with multilevel correction         | high-order integration         |
| MGRiT     | a hierarchy of temporal grids             | relaxation and coarse-grid correction  | long time intervals            |
| STMG      | the full space-time grid                  | space-time smoothing and coarsening    | parabolic systems              |

## Three organizing principles

### Causality is reformulated, not removed

For a one-step method, $u_{n+1}=\Phi_{\Delta t}(u_n)$ is sequential. A parallel-in-time method instead constructs a parallel approximation to the inverse of the coupled all-at-once system. Parallel work occurs inside an iteration, a decomposition, or a transform.

### Dissipation determines whether a coarse representation is informative

Parabolic dynamics damp high-frequency error. A coarse propagator can therefore reproduce the slowly varying components that remain relevant over long intervals. Hyperbolic dynamics preserve phase information. Small phase errors in the coarse problem may then accumulate across time intervals.

### Iteration count is not parallel efficiency

A simplified cost model is

$$
T_{\mathrm{parallel}}
\approx K(C_G+C_F/P)+C_{\mathrm{comm}}+C_{\mathrm{setup}},
$$

where $K$ is the iteration count, $C_G$ and $C_F$ are the coarse and fine propagation costs, and $P$ is the temporal concurrency. The numerical experiments in these notes measure convergence, not end-to-end parallel speedup.

> [!note] Numerical provenance
> All displayed experiment figures were regenerated on 31 July 2026 from the Python/SciPy reproduction project on the original T4 host. The implementation uses CPU sparse linear algebra and FFTs. The presence of a T4 accelerator does not imply that the reported runs use GPU kernels.

## Primary sources

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), Acta Numerica 34 (2025), pp. 385-489.
- Original MATLAB examples: [wushulin/ActaPinT](https://github.com/wushulin/ActaPinT).
