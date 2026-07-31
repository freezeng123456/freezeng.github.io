---
title: "Chapter 5: A Unified View and Method Selection"
description: A common algebraic interpretation, an experiment ledger, and reproducibility limits
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view
tags:
  - parallel-in-time
  - methodology
---

## A common algebraic form

Many PinT iterations can be written as

$$
U^{k+1}=U^k+M^{-1}(b-AU^k),
$$

where $A$ is the all-at-once space-time operator and $M^{-1}$ is a parallel approximation to its inverse.

| Method   | Source of $M^{-1}$                                            |
| -------- | ------------------------------------------------------------- |
| Parareal | block lower-triangular inverse defined by a coarse propagator |
| MGRiT    | a temporal multilevel cycle                                   |
| STMG     | a multilevel cycle on the complete space-time grid            |
| ParaDiag | FFT inversion of a circulant temporal approximation           |
| SWR      | local space-time inverses coupled by waveform transmission    |

This formulation makes the central questions explicit: which error modes does $M^{-1}$ reduce, which physical information does it preserve, and which operations are actually concurrent?

## Method selection

![A decision map for parallel-in-time methods](assets/diagrams/pint/en/method-selection.svg)

- **Strongly dissipative, low- or moderate-order integration:** MGRiT and STMG are natural candidates; Parareal provides a simple initial prototype.
- **High-order temporal accuracy:** PFASST incorporates collocation and SDC but requires a more complex schedule.
- **Large linear systems:** ParaDiag is attractive when shifted spatial solves are scalable; ParaExp is attractive when matrix-exponential actions are efficient.
- **Transport- or wave-dominated problems:** characteristic transmission, phase correction, or SWR can preserve propagation more faithfully than a strongly dissipative coarse solver.

## Complete experiment ledger

The Python reproduction exposes eight baseline experiments and one composite paper-validation entry. The composite entry generates six paper-specific figures. All 14 SVG/PNG figure sets and their numerical conclusions are now assigned to Chapters 2-4.

| Python output               | Website location               | Machine-readable result                              |
| --------------------------- | ------------------------------ | ---------------------------------------------------- |
| `solution_heat_ade`         | Chapter 2, advection-diffusion | [[assets/pint/data/solution_heat_ade.json            | JSON]] |
| `solution_burgers`          | Chapter 2, Burgers             | [[assets/pint/data/solution_burgers.json             | JSON]] |
| `solution_wave`             | Chapter 2, wave equation       | [[assets/pint/data/solution_wave.json                | JSON]] |
| `parareal_heat_ade`         | Chapter 4, Parareal            | [[assets/pint/data/parareal_heat_ade.json            | JSON]] |
| `parareal_burgers`          | Chapter 4, Parareal            | [[assets/pint/data/parareal_burgers.json             | JSON]] |
| `mgrit_heat_ade`            | Chapter 4, MGRiT               | [[assets/pint/data/mgrit_heat_ade.json               | JSON]] |
| `iterative_paradiag_ade`    | Chapter 3, ParaDiag            | [[assets/pint/data/iterative_paradiag_ade.json       | JSON]] |
| `stmg_heat_ade`             | Chapter 4, STMG                | [[assets/pint/data/stmg_heat_ade.json                | JSON]] |
| Figure 3.15 validation      | Chapter 3, ParaDiag            | [[assets/pint/data/figure_3_15_validation.json       | JSON]] |
| Figure 4.5 validation       | Chapter 4, Parareal            | [[assets/pint/data/figure_4_5_validation.json        | JSON]] |
| Figures 4.9-4.10 validation | Chapter 4, MGRiT               | [[assets/pint/data/figure_4_10_validation.json       | JSON]] |
| Figure 4.19 validation      | Chapter 4, STMG                | [[assets/pint/data/figures_4_19_4_20_validation.json | JSON]] |
| Figure 4.20 validation      | Chapter 4, STMG                | [[assets/pint/data/figures_4_19_4_20_validation.json | JSON]] |

The compact cross-experiment record is available as [[assets/pint/data/paper_validation_summary.json|paper_validation_summary.json]].

The upstream MATLAB repository contains additional scripts for direct ParaDiag, diagonalized Parareal, ParaExp, SWR, IDC/PIDC, and wave-domain decomposition. They have been inventoried but not all ported to the current Python experiment interface; the wave ParaDiag case in Figure 3.15 is now implemented. The ledger therefore claims complete coverage of generated Python result artifacts, not complete reproduction of every MATLAB script in the upstream repository.

## Formal reproduction

```bash
python3.11 -m pip install -e ".[test]"
actapint all --quick --output-dir results/quick
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  actapint paper_validation --output-dir results/paper-full
```

The formal run on 31 July 2026 used Python 3.11, NumPy 2.4.6, SciPy 1.17.1, and Matplotlib 3.11.1 on the new experiment server. The paper suite took about 4 minutes 24 seconds of wall time and less than 280 MiB peak resident memory. The code path uses CPU sparse factorizations, Krylov methods, and FFTs.

Each experiment writes:

1. an editable SVG and a high-resolution PNG generated from the same arrays used for analysis;
2. a JSON record of the grid, physical parameters, tolerance, and metric;
3. deterministic random initialization where an all-at-once initial iterate is required.

## Interpretation limits

- The experiments measure numerical convergence and reproduce selected trends and values from the paper. They do not measure temporal strong or weak scaling.
- No figure reports MPI process count, communication volume, setup cost, or wall-clock speedup.
- The original MATLAB random generator and NumPy do not produce identical initial arrays. Convergence factors, iteration counts, and final states are the relevant comparisons.
- Sparse factorization, FFT ordering, and GMRES reduction can produce ordinary differences near $10^{-14}$ to $10^{-16}$.
- In `MGRiT_Heat_ADE.m`, the invalid expression `nu=0.002max;` is interpreted as $\nu=0.002$, which is consistent with the surrounding branches and the paper.
- For the STMG paper validation, backward Euler and the original MATLAB residual convention are retained. A consistent post-smoothing residual is also stored in JSON for diagnostic comparison.

## Minimum reporting standard for future experiments

A performance-oriented PinT study should report the spatial and temporal discretizations, fine reference, number of time intervals, hardware allocation, coarse-to-fine cost ratio, stopping criterion, communication and setup time, parameter sweeps, and strong or weak scaling. Error versus iteration alone supports an algorithmic convergence statement, not a parallel-efficiency claim.
