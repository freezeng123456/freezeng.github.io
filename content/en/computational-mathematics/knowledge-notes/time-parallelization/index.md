---
title: Time Parallelization for Hyperbolic and Parabolic Problems
description: Chapter-by-chapter notes and reproducible experiments for parallel-in-time methods
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization
tags:
  - computational-mathematics
  - parallel-in-time
---

These notes follow M. J. Gander, S.-L. Wu, and T. Zhou, _Time Parallelization for Hyperbolic and Parabolic Problems_, Acta Numerica 34 (2025), pp. 385-489. The survey distinguishes methods that remain effective for propagative problems from methods designed primarily for dissipative problems.

## Reading sequence

1. [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-1-why-parallelize-in-time|Chapter 1: Why Parallelize in Time?]] follows the abstract and introduction paragraph by paragraph, covering the hardware context, the causal chain, four historical lineages, the two-way classification, and the all-at-once system.
2. [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-2-model-problems|Chapter 2: Model Problems]] compares boundary conditions and temporal memory across heat, advection–diffusion, Burgers, and wave equations, then presents three recomputed solutions.
3. [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-hyperbolic-methods|Chapter 3: Methods Effective for Hyperbolic Problems]] derives SWR, PIDC/RIDC, ParaExp, and ParaDiag and includes ParaDiag-II experiments for heat, ADE, and wave problems.
4. [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-parabolic-methods|Chapter 4: Methods Designed Primarily for Parabolic Problems]] analyzes Parareal, PFASST, MGRiT, two diagonalized variants, and STMG, with recomputed convergence studies.
5. [[en/computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view|Chapter 5: Conclusions, a Unified View, and Reproduction Boundaries]] separates the paper conclusion from the site's method selection, experiment inventory, GPU performance, and reporting rules.

## Paper-coverage progress

“Paragraph-level complete” means that the claims, equations, figures, historical links, qualifications, and section relationships have all been checked against the source. Existing experiments remain in place after the corresponding source discussion.

| Source range           | Website chapter | Current status                                                                                                                                               |
| ---------------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Abstract and Section 1 | Chapter 1       | **Paragraph-level complete**: every argument on pp. 385–388 is covered, and Figure 1.1 has been redrawn                                                      |
| Sections 2.1–2.4       | Chapter 2       | **Paragraph-level complete**: four models, every boundary setting, all Figure 2.1–2.4 observation groups, and three supplemental experiments                 |
| Sections 3.1–3.5.2     | Chapter 3       | **Paragraph-level complete**: history, SWR, PIDC/RIDC, ParaExp, ParaDiag-I/II derivations, theorems, and all numerical discussions                           |
| Sections 4.1–4.6       | Chapter 4       | **Paragraph-level complete**: Parareal, PFASST, MGRiT, diagonalized variants, STMG, Theorems 4.1–4.9, and Figures 4.1–4.22                                   |
| Section 5              | Chapter 5       | **Paragraph-level complete**: paper conclusions are covered; unified analysis, GPU work, and the experiment ledger are explicitly marked as site supplements |

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

### How causal constraints enter a parallel computation

For a one-step method, $u_{n+1}=\Phi_{\Delta t}(u_n)$ forms a sequential recurrence. A parallel-in-time method exposes the coupling among all temporal unknowns and constructs a parallel approximation to the inverse. Concurrency then occurs inside an iteration, decomposition, or transform.

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
> All displayed experiment figures were regenerated on 31 July 2026 from the Python reproduction project on the new experiment server. The initial formal results used the SciPy CPU path. A subsequent CuPy/T4 hybrid backend batches the independent Burgers fine propagators on the GPU, preserves the Figure 4.5 stopping iterations, and reduces the complete paper suite from 263.57 to 67.92 seconds.

## Primary sources

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), Acta Numerica 34 (2025), pp. 385-489.
- Original MATLAB examples: [wushulin/ActaPinT](https://github.com/wushulin/ActaPinT).
- Python conversion, extensions, and formal results: [freezeng123456/ActaPinT-Python](https://github.com/freezeng123456/ActaPinT-Python).
