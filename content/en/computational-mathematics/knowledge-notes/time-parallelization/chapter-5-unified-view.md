---
title: "Chapter 5: Conclusions"
description: Paper conclusion (the hyperbolic/parabolic temporal-memory criterion), methods recommended by class, full experiment inventory, GPU optimization, and minimum reporting standards
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view
tags:
  - parallel-in-time
  - methodology
---

> [!note] Content boundary
> Section 5 of the paper (p. 481) contains no numbered subsections or numbered equations. This page first explains the source conclusion paragraph by paragraph. The subsequent unified algebraic view, method selection, Python experiments, T4 GPU performance, and reporting standards are marked as site supplements and do not use 5.x numbering.

## Section 5: Conclusions

The paper's conclusion condenses the whole analysis into a single criterion: what decides whether a given parallel-in-time method applies is not the algorithm's name, nor a traditional classification such as "iterative/direct" or "Section 3/Section 4", but the problem's own **temporal memory**—how much fine information the solution retains as it evolves.

**Parabolic problems forget information quickly, and their solutions are temporally local.** The heat equation and sufficiently diffusive reaction–diffusion equations dissipate high-frequency components rapidly, so the solution at any instant depends mainly on recent history. Even a cheap coarse temporal model that discards high frequencies can still capture the dominant slow modes, and the error between coarse prediction and fine correction decays across iterations. The paper notes that this class admits many highly effective PinT methods, including Parareal, space–time multigrid (STMG), ParaExp, ParaDiag, and waveform relaxation (WR) built on domain decomposition (DD).

**Hyperbolic problems retain fine structure, phase, and propagation paths over very long horizons.** Waves, transport, and low-viscosity conservation laws do not dissipate high frequencies, so errors lack a natural decay mechanism; any phase or amplitude mismatch between coarse and fine propagation accumulates along characteristics and amplifies as the number of time intervals grows. Only a subset of methods therefore remain effective: ParaExp, ParaDiag, and Schwarz waveform relaxation (SWR)—the paper particularly emphasizes the connection between SWR and tent pitching. What these methods share is that they do not rely on a dissipative coarse temporal model; instead they move long-range information through exact matrix-exponential propagation, all-at-once algebraic (frequency) structure, or subdomain solves along the characteristic cone.

> [!tip] Insight
> Temporal memory can be read off from the spectrum of the discrete evolution operator. For parabolic discretizations the eigenvalues cluster near the negative real axis, $e^{\lambda\Delta t}$ strongly damps high-frequency modes, and coarse propagation only needs to track the few slowly decaying low-frequency modes; cross-iteration error is suppressed by the same decay mechanism. For hyperbolic discretizations the eigenvalues sit close to the imaginary axis, modes rotate approximately as $e^{i\omega t}$ with almost no decay, and the phase difference between coarse and fine propagation is not absorbed by damping but instead accumulates along characteristics. This is the algebraic expression of "parabolic forgets, hyperbolic remembers", and it explains why a coarse model that only seeks dissipative smoothing fails on hyperbolic problems.

**Further reading and code.** The paper recommends the research monograph by Gander and Lunet (2024), _Time Parallel Time Integration_ (SIAM), as systematic reading: it provides, for each PinT method, historical context, concise yet fully self-contained convergence analyses, and short directly runnable MATLAB programs. The paper also states that **the code used to generate all of the paper's results** is public at [wushulin/ActaPinT](https://github.com/wushulin/ActaPinT); that is, every experiment and figure in the entire paper can be reproduced from that repository, not just a few figures.

This conclusion offers a practical first-pass principle: first judge whether the problem quickly forgets high-frequency information, then choose the parallel-in-time structure. An algorithm's name or a traditional classification cannot replace this dynamical judgment by itself.

> [!tip] Insight
> This criterion is orthogonal to traditional classifications such as "iterative/direct" and "Section 3/Section 4". Methods within the same classification may fall on opposite sides of the criterion: the iterative SWR is hyperbolic-friendly, whereas the equally iterative standard Parareal/MGRiT, which rely on a dissipative coarse level, tend to fail on hyperbolic problems; conversely, ParaExp and ParaDiag, because they do not rely on a dissipative coarse model (moving long-range information through exact matrix-exponential propagation and all-at-once frequency structure, respectively), are usable on both the parabolic and hyperbolic sides. On the parabolic side, the Section 4 multilevel methods in the paper (PFASST, MGRiT) and the Section 3 methods (SWR, IDC/PIDC/RIDC, ParaExp, ParaDiag) apply equally well; on the hyperbolic side one can additionally add methods such as parallel IDC that pipeline along a high-order error equation without forcing dissipative coarsening. The grouping by class above is a site synthesis; the paper's conclusion itself explicitly lists only the few representative methods named earlier.

> [!tip] Insight
> The natural questions after the conclusion center on how to construct a "memory-friendly" coarse level for hyperbolic problems and how to turn iterative convergence into real wall-clock gains:
>
> - Transport-oriented coarse operators: replace naive coarse propagation with semi-Lagrangian or corrected coarse-grid operators to ease phase mismatch (e.g., the MGRiT work of De Sterck et al. on linear advection);
> - Coarse-propagation-free Parareal: explore variants that do not rely on sequential coarse prediction, to shorten the serial tail (Gander, Ohlberger, Rave 2024);
> - Mixed precision and roundoff control: balance $\alpha$, outer Krylov, and floating-point error in ParaDiag (Wu, Yang, Zhou 2025);
> - Phase awareness and characteristic preservation: combine tent pitching, SWR transmission conditions, and all-at-once algebraic structure for waves and conservation laws;
> - Multilevel and multi-device scaling: design temporal concurrency, spatial concurrency, and communication overlap together, and build strong/weak scaling evidence.
>   These directions are a site synthesis based on the paper's main text and references; the paper's conclusion itself gives only directional guidance and leaves the systematic development to the Gander and Lunet (2024) monograph.

## Site synthesis: representative literature for the methods in the conclusion

The table below maps the methods mentioned in the conclusion to representative works in the paper's references, for convenient lookup by class. This is a site-compiled search aid, not part of the original conclusion text; class membership follows the paper's "temporal memory" criterion.

| Method            | Applicable class                                      | Representative literature (all in the paper's references)                             |
| ----------------- | ----------------------------------------------------- | ------------------------------------------------------------------------------------- |
| Parareal          | parabolic (usable as a hyperbolic baseline)           | Lions, Maday & Turinici (2001); convergence analysis Gander & Vandewalle (2007)       |
| STMG              | parabolic                                             | Gander & Neumüller (2016); Horton & Vandewalle (1995)                                 |
| PFASST            | parabolic (high-order collocation)                    | Emmett & Minion (2012); Bolten, Moser & Speck (2017)                                  |
| MGRiT             | parabolic; advection needs a modified coarse operator | Falgout et al. (2014); De Sterck et al. (2021, 2023a)                                 |
| ParaExp           | parabolic and hyperbolic                              | Gander & Güttel (2013); nonlinear Gander, Güttel & Petcu (2018)                       |
| ParaDiag          | parabolic and hyperbolic                              | Gander et al. (2021c); Gander & Wu (2020); McDonald, Pestana & Wathen (2018)          |
| SWR / OSWR        | parabolic and hyperbolic                              | Gander & Stuart (1998); Gander & Halpern (2007); waves Gander, Halpern & Nataf (2003) |
| tent pitching     | hyperbolic                                            | Gopalakrishnan, Schöberl & Wintersteiger (2017); Ciaramella, Gander & Mazzieri (2023) |
| IDC / PIDC / RIDC | parabolic; parallel IDC also for hyperbolic           | Dutt, Greengard & Rokhlin (2000); Christlieb, Macdonald & Ong (2010)                  |

## Site synthesis: solving one all-at-once system

For a linear all-at-once discretization

$$
A\boldsymbol U=\boldsymbol b,
$$

many PinT iterations can be written as

$$
\boldsymbol U^{k+1}
=\boldsymbol U^k+M^{-1}
(\boldsymbol b-A\boldsymbol U^k).
$$

$M^{-1}$ is a parallel approximation to $A^{-1}$. Each method chooses a different locality, hierarchy, or transform:

| Method    | Main source of $M^{-1}$                           | Concurrent work                              | Carrier of long-range information            |
| --------- | ------------------------------------------------- | -------------------------------------------- | -------------------------------------------- |
| SWR       | space–time subdomain inverses and transmission    | complete subdomain waveforms                 | interface waveforms and characteristic cones |
| PIDC/RIDC | integral residual correction                      | window or correction-level pipeline          | high-order error equation                    |
| ParaExp   | local inhomogeneous solves and exponential action | local forced responses and homogeneous tails | exact $e^{tA}$ propagation                   |
| ParaDiag  | circulant/diagonalizable time operator            | shifted spatial systems after FFT            | global temporal frequencies                  |
| Parareal  | lower-triangular coarse propagator inverse        | fine solves on large intervals               | sequential coarse prediction                 |
| PFASST    | multilevel preconditioner for collocation         | SDC sweeps on time steps                     | FAS coarse collocation correction            |
| MGRiT     | temporal multilevel cycle                         | F relaxation and coarse levels               | C points and overlapping relaxation          |
| STMG      | full space–time multilevel cycle                  | time-block Jacobi                            | coarse space–time grids                      |

This common form poses three questions. Which error modes does $M^{-1}$ reduce? Does it preserve phase, mean value, and shock position? During one application of $M^{-1}$, which operations are genuinely concurrent and which remain sequential?

> [!tip] Insight
> Viewed through this unified form, the conclusion's criterion can be restated as a spectral requirement on $M^{-1}$: parabolic problems allow $M^{-1}$ to be accurate only at low frequencies and to strongly damp high frequencies, because the true solution's high frequencies are also decaying; hyperbolic problems require $M^{-1}$ to preserve phase and amplitude across the frequency band of interest, and any approximation that "only smooths without transporting phase" drives the spectral radius of $I-M^{-1}A$ toward or even beyond 1 near the hyperbolic band. ParaExp and ParaDiag are usable across both classes precisely because their $M^{-1}$ comes from exact exponentials and all-at-once frequency diagonalization, respectively, rather than from a dissipative coarse level.

## Site synthesis: method-selection map

![Map for selecting parallel-in-time methods](assets/diagrams/pint/en/method-selection.svg)

### Strong dissipation and low-to-moderate temporal order

Heat and sufficiently diffusive reaction–diffusion systems are natural candidates for MGRiT and STMG. Parareal provides a quick test of coarse-propagator quality and a nonintrusive baseline. STMG has greater scalability potential but requires access to the all-at-once operator, smoother, and grid transfers.

### High-order collocation

PFASST is attractive when high temporal order is required and a collocation solve is expensive. Nodes, SDC sweeps, the coarse collocation level, and spatial parallel resources need to be designed together. High formal order does not guarantee an efficient temporal pipeline.

### Large linear or linearized systems

ParaExp accurately transports long-range linear information when exponential action scales. ParaDiag removes temporal forward substitution through FFTs when complex shifted spatial systems have capable solvers. ParaDiag-I is limited by conditioning of the temporal eigenvectors. ParaDiag-II must also balance $\alpha$, outer Krylov convergence, and roundoff.

### Transport, waves, and low-viscosity nonlinearity

Prioritize structures that represent characteristics and phase, including SWR/OSWR, tent pitching, ParaExp, $\alpha$-ParaDiag, and phase-aware coarse propagation. Standard Parareal, MGRiT, and STMG remain useful diagnostic baselines. If fine/coarse phase mismatch grows with frequency, increasing the interval count is likely to amplify the difficulty.

### Six questions before choosing a method

1. Do important modes lie near the negative real axis or close to the unit circle/imaginary axis?
2. Do boundaries permit outflow, recirculate a periodic signal, or reflect waves?
3. Does nonlinearity continuously create shocks and high frequencies?
4. Which reusable component is available: a time stepper, shifted spatial solver, exponential action, or all-at-once operator?
5. Is the goal lower iteration count, higher single-node throughput, or multi-node strong/weak scaling?
6. How much intrusive modification, global transformation, and all-time storage is acceptable?

## Site synthesis: parameter reference

| Parameter                            | Location                       | Direct role                              | Quantities to monitor together                            |
| ------------------------------------ | ------------------------------ | ---------------------------------------- | --------------------------------------------------------- |
| SWR overlap and Robin $p$            | subdomain interface            | controls waveform transfer               | window length, viscosity, interface cost                  |
| IDC node count $M$ and corrections   | error equation                 | limit formal order and pipeline depth    | regularity, fill and drain time                           |
| Parareal intervals $N$ and ratio $J$ | fine/coarse propagation        | set concurrency and mismatch             | iteration count, sequential coarse cost, phase error      |
| MGRiT coarsening and CF count        | time hierarchy                 | set overlap contraction and fine work    | total factor at equal fine-solve cost                     |
| ParaDiag $\alpha$                    | cyclic head–tail approximation | smaller values improve the approximation | $\epsilon/\alpha$ roundoff and shifted-solve stability    |
| STMG damping $\eta$                  | time-block Jacobi              | controls high-frequency smoothing        | time integrator, cycle cost, spatial coarsening condition |

Each parameter interacts with the physical spectrum, discrete stability function, and machine cost. Minimizing iteration count in isolation can move the work into a much more expensive iteration.

## Site reproduction: experiment inventory

The Python reproduction project provides eight baseline experiments and one combined paper-validation entry point. The combined entry point generates six paper-matched plots. Chapters 2–4 reference fourteen SVG/PNG result groups with corresponding JSON records.

| Python output               | Page location                  | Machine-readable record                              |
| --------------------------- | ------------------------------ | ---------------------------------------------------- |
| `solution_heat_ade`         | Chapter 2, advection–diffusion | [[assets/pint/data/solution_heat_ade.json            | JSON]] |
| `solution_burgers`          | Chapter 2, Burgers             | [[assets/pint/data/solution_burgers.json             | JSON]] |
| `solution_wave`             | Chapter 2, wave                | [[assets/pint/data/solution_wave.json                | JSON]] |
| `parareal_heat_ade`         | Chapter 4, Parareal            | [[assets/pint/data/parareal_heat_ade.json            | JSON]] |
| `parareal_burgers`          | Chapter 4, Parareal            | [[assets/pint/data/parareal_burgers.json             | JSON]] |
| `mgrit_heat_ade`            | Chapter 4, MGRiT               | [[assets/pint/data/mgrit_heat_ade.json               | JSON]] |
| `iterative_paradiag_ade`    | Chapter 3, ParaDiag            | [[assets/pint/data/iterative_paradiag_ade.json       | JSON]] |
| `stmg_heat_ade`             | Chapter 4, STMG                | [[assets/pint/data/stmg_heat_ade.json                | JSON]] |
| Figure 3.15 validation      | Chapter 3, ParaDiag            | [[assets/pint/data/figure_3_15_validation.json       | JSON]] |
| Figure 4.5 validation       | Chapter 4, Parareal            | [[assets/pint/data/figure_4_5_validation.json        | JSON]] |
| Figures 4.9–4.10 validation | Chapter 4, MGRiT               | [[assets/pint/data/figure_4_10_validation.json       | JSON]] |
| Figure 4.19 validation      | Chapter 4, STMG                | [[assets/pint/data/figures_4_19_4_20_validation.json | JSON]] |
| Figure 4.20 validation      | Chapter 4, STMG                | [[assets/pint/data/figures_4_19_4_20_validation.json | JSON]] |
| T4 GPU validation           | this chapter, GPU acceleration | [[assets/pint/data/gpu_benchmark_t4.json             | JSON]] |

The cross-experiment summary is [[assets/pint/data/paper_validation_summary.json|paper_validation_summary.json]].

The upstream MATLAB repository also contains direct ParaDiag, diagonalized Parareal, ParaExp, SWR, IDC/PIDC, and wave-domain-decomposition scripts. They are registered in the Python migration inventory but do not yet all have formal Python results. "Complete" here means that every Python artifact cited by the site has a matched parameter record, plot, and JSON file. It does not claim that every upstream MATLAB script has been ported.

## Site reproduction: formal run workflow

```bash
python3.11 -m pip install -e ".[test]"
actapint all --quick --output-dir results/quick
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  actapint paper_validation --output-dir results/paper-full
```

The formal run on July 31, 2026 used Python 3.11, NumPy 2.4.6, SciPy 1.17.1, and Matplotlib 3.11.1. The CPU paper suite took about 4 minutes 24 seconds wall time and stayed below 280 MiB peak resident memory. The code path includes CPU sparse factorization, Krylov methods, and FFTs.

Each experiment writes:

1. editable SVG and high-resolution PNG produced from the same analysis arrays;
2. JSON containing grid, physical parameters, tolerance, stopping convention, and metrics;
3. deterministic seeds where an all-at-once random initial vector is needed;
4. a separate validation summary for paper-matched experiments.

## Site reproduction: GPU acceleration and profiling

Function-level profiling attributes 43.06 of 62.95 seconds in the quick paper suite to Figure 4.5, including 38.11 seconds in Burgers fine propagation. All FFTs take only 0.007 seconds and all GMRES calls total 0.251 seconds. The first CUDA backend therefore batches the 40 independent Burgers fine propagations in each Parareal iteration.

The CuPy backend keeps spatial operators resident on the GPU and assembles and solves 40 independent Newton systems in batches. Causal coarse propagation and the remaining experiments stay on the CPU, yielding a hybrid CPU/GPU implementation.

| T4 double-precision test                 |      CPU |     GPU | Speedup |
| ---------------------------------------- | -------: | ------: | ------: |
| 40 Burgers fine propagators, 32 substeps |  2.893 s | 0.246 s |  11.76× |
| complete paper-validation suite          | 263.57 s | 67.92 s |   3.88× |

The maximum absolute CPU/GPU difference for one batch is $2.33\times10^{-15}$. Figure 4.5 retains stopping iterations ADE 14/24/35 and Burgers 14/21/25. Continuing beyond the $10^{-10}$ target toward machine precision produces normal rounding differences between CPU SuperLU and GPU batched LU because their floating-point operation orders differ.

```bash
python3.11 -m pip install -e ".[gpu,test]"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  actapint paper_validation --backend gpu \
  --output-dir results/paper-gpu
```

### Next optimization opportunities

1. Replace dense batched LU with a cyclic tridiagonal CUDA kernel so storage and work scale linearly with $N_x$.
2. Keep the complete Parareal state on the GPU and investigate parallel-prefix/scan coarse propagation to reduce host transfers and the sequential tail.
3. Implement batched complex shifted banded solves on larger ParaDiag grids. FFT and GMRES work on the current $100\times100$ grid is too small for stable GPU benefit.
4. Separate spatial concurrency, temporal concurrency, and communication overlap in multi-GPU strong- and weak-scaling tests.

The machine-readable record is [[assets/pint/data/gpu_benchmark_t4.json|gpu_benchmark_t4.json]]. Current data demonstrate single-GPU kernel and end-to-end acceleration and do not establish multi-GPU scaling.

## Site reproduction: interpretation boundaries

- Current experiments measure numerical convergence and single-node CPU/GPU performance, without time-dimensional MPI strong or weak scaling.
- Site figures do not report MPI rank count, network volume, initialization cost, or cross-node wall-clock speedup.
- MATLAB and NumPy random generators do not produce identical initial arrays; convergence factors, phases, and final states are better comparison targets.
- Sparse factorization, FFT ordering, and GMRES reductions can differ normally near $10^{-14}$ to $10^{-16}$.
- The invalid expression `nu=0.002max;` in `MGRiT_Heat_ADE.m` is interpreted as $\nu=0.002$ from its branch context and the paper.
- Source Figure 4.9 annotates the maximum Parareal factor at $\nu=0.01$ as $0.9986$, while equation (4.5b) and the upstream-script stability function give $1.0501$ in Python. The site retains both and records this as a reproduction discrepancy.
- STMG validation retains the upstream backward-Euler residual convention and stores a consistent postsmoothing residual separately in JSON.
- Convergence to the sequential fine solution establishes algorithmic consistency and gives no wall-clock speedup guarantee.
- The large-scale STMG data in Table 4.1 belong to the cited three-dimensional parallel implementation and are not measurements from the present Python project.

## Site standard: minimum reporting requirements for future experiments

An algorithmic-convergence experiment should state:

- PDE, boundaries, initial data, and source;
- spatial and temporal discretizations and the fine-grid reference;
- fine/coarse propagators, interval count, and coarsening;
- error norm, residual definition, stopping threshold, and iteration cap;
- parameter sweeps and observed failure points.

A parallel-performance experiment should additionally state:

- CPU/GPU models, precision, and rank/thread/device allocation;
- fine/coarse cost ratio, concurrent tasks per iteration, and load balance;
- initialization, transfer, communication, synchronization, and I/O time;
- total wall time, speedup, efficiency, and baseline implementation;
- strong or weak scaling as the interval count, spatial size, and device count vary.

Error decay by iteration supports a numerical-convergence claim but cannot support a parallel-efficiency claim on its own. Method comparisons must also normalize fine-propagation count or total work so curves with different per-iteration costs are not ranked directly.

## Site-wide coverage table

| Source range                  | Site chapter    | Completeness statement                                                                                                                    |
| ----------------------------- | --------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| Section 2, pp. 388–396        | Chapter 2       | four models, all boundary settings, every Figure 2.1–2.4 observation group, and PinT implications                                         |
| Sections 3.1–3.2, pp. 396–405 | Chapter 3.1–3.2 | history, WR/SWR, Theorems 3.1–3.2, OSWR, MTP/UTP                                                                                          |
| Sections 3.3–3.4, pp. 405–415 | Chapter 3.3–3.4 | IDC/PIDC/RIDC derivation and regularity tests, linear and nonlinear ParaExp                                                               |
| Section 3.5, pp. 415–443      | Chapter 3.5     | ParaDiag-I/II, Theorems 3.5–3.9, BVM, NKA, circulant and $\alpha$-circulant experiments                                                   |
| Sections 4.1–4.4, pp. 443–460 | Chapter 4.1–4.4 | Parareal, PFASST, MGRiT, Theorems 4.1–4.6, Figures 4.1–4.11                                                                               |
| Sections 4.5–4.6, pp. 460–481 | Chapter 4.5–4.6 | both diagonalized Parareal variants, STMG, Theorems 4.7–4.9, Figures 4.12–4.22, Table 4.1                                                 |
| Section 5, p. 481             | Chapter 5       | hyperbolic/parabolic temporal-memory summary, methods recommended by class, the monograph, and the public code that generates all results |

Each chapter ends with a more granular source-page audit. Site supplements occupy explicitly labeled sections and are not blended into claims attributed to the paper.

## Summary

The starting point for PinT selection is the temporal memory of the dynamics. Strong diffusion permits coarse temporal levels to represent the remaining slow modes. Transport, waves, and low-viscosity nonlinear systems require preservation of phase, characteristics, and shock position. The unified all-at-once view helps compare algorithms, while a successful implementation must satisfy three conditions together: iterations to tolerance remain controlled, concurrent work dominates runtime, and communication and memory scale. A complete reproduction reports evidence for each layer separately.

## Source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 5, p. 481.
