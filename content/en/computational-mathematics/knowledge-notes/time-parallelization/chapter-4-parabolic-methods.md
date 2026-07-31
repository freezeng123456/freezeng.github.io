---
title: "Chapter 4: Methods Designed for Parabolic Problems"
description: Parareal, PFASST, MGRiT, diagonalized coarse correction, and space-time multigrid
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-4-parabolic-methods
tags:
  - parallel-in-time
  - parabolic-PDE
---

Parabolic dynamics damp high-frequency components and make coarse propagation informative. This property supports the most widely studied coarse-grid and multilevel PinT methods.

## 4.1 Parareal

Let $F$ be an accurate propagator and $G$ an inexpensive coarse propagator over one time interval. Parareal applies

$$
U_{n+1}^{k+1}
=G(U_n^{k+1})+F(U_n^k)-G(U_n^k).
$$

The fine solves for a fixed iteration are independent across intervals. The new coarse propagation remains sequential. In exact arithmetic the method has a finite-step property, but convergence in nearly $N$ iterations gives little opportunity for speedup over $N$ time intervals.

### Recomputed baseline experiments

| Problem             | Formal parameters                       | Iterations |   Final maximum error |
| ------------------- | --------------------------------------- | ---------: | --------------------: |
| advection-diffusion | $N_x=128$, $N_t=40$, $J=32$, $\nu=0.02$ |         39 | $5.940\times10^{-15}$ |
| Burgers             | $N_x=128$, $N_t=40$, $J=32$, $\nu=1$    |         16 | $3.577\times10^{-12}$ |

![[assets/pint/parareal-ade-baseline.png]]

![[assets/pint/parareal-burgers-baseline.png]]

The iteration counts refer to convergence toward the serial fine solution. They are not wall-clock speedups.

### Diffusion sweep corresponding to Figure 4.5

The formal sweep fixes $T=4$, $\Delta T=0.1$, $\Delta x=1/128$, and $J=32$. The number of iterations required to reduce the maximum error below $10^{-10}$ is:

| Equation            | $\nu=1$ | $\nu=0.1$ | $\nu=0.02$ |
| ------------------- | ------: | --------: | ---------: |
| advection-diffusion |      14 |        24 |         35 |
| Burgers             |      14 |        21 |         25 |

![[assets/pint/parareal-figure-4-5.png]]

The data show slower convergence as diffusion weakens for both equations under this fixed protocol.

## 4.2 PFASST

PFASST combines spectral deferred correction with a multilevel full approximation scheme:

1. each time step contains several collocation nodes;
2. SDC sweeps approach the high-order collocation solution;
3. a coarse level transfers slowly converging components;
4. sweeps on different time steps execute as a pipeline.

The method is appropriate when high-order temporal accuracy is required, but its performance depends on the collocation nodes, number of sweeps, coarse operator, and coordination with spatial parallelism.

## 4.3 MGRiT

MGRiT separates the temporal grid into C-points and intervening F-points. F-relaxation fills the F-points in parallel. C-relaxation and coarse-grid correction transmit long-range information. FCF relaxation performs an F, C, and F sequence and is usually more robust than F relaxation alone, at a higher cost per iteration.

### Baseline near-hyperbolic experiment

The baseline comparison uses $N_x=160$, $N_t=40$, $J=20$, $T=5$, and $\nu=0.002$. After the reported iteration window, the Parareal maximum error is $2.895\times10^2$, whereas two-level MGRiT reaches $4.441\times10^{-16}$ through its finite-step behavior.

![[assets/pint/mgrit-baseline.png]]

The late finite-step drop does not make the preceding iteration useful as a scalable solver. The transient growth remains the relevant observation for a practical stopping tolerance.

### Equal-cost comparison corresponding to Figures 4.9-4.10

One FCF-MGRiT iteration uses approximately two fine propagations and is therefore compared with two Parareal iterations:

| Problem                          | Parareal factor | FCF-MGRiT factor |
| -------------------------------- | --------------: | ---------------: |
| heat equation                    |          0.2824 |           0.0835 |
| advection-diffusion, $\nu=0.1$   |          0.4453 |           0.2719 |
| advection-diffusion, $\nu=0.01$  |          1.0501 |           0.9021 |
| advection-diffusion, $\nu=0.002$ |          1.4211 |           1.2812 |

![[assets/pint/mgrit-figure-4-10.png]]

Both methods deteriorate as the spectrum approaches the hyperbolic regime. At $\nu=0.002$, the factors exceed one during the practically relevant linear phase.

## 4.4 Diagonalization-based Parareal

Diagonalization-based variants apply a temporal FFT to the coarse-grid correction or replace the sequential coarse solve by a circulant approximation. They reduce the serial fraction but introduce complex shifted spatial solves, global transforms, and approximation error from the cyclic coupling.

## 4.5 Space-time multigrid

STMG treats the complete discretization as a multilevel space-time system. Its effectiveness depends on the smoother, damping $\eta$, coarsening strategy, and time integrator.

### Baseline trapezoidal-rule experiment

The baseline experiment uses $N_x=N_t=255$, $\nu=10^{-3}$, the trapezoidal rule, and three pre- and post-smoothing steps. In the sampled grid, the minimum error after 15 cycles occurs near $\eta=0.98$.

![[assets/pint/stmg-baseline.png]]

This value is not comparable to the backward-Euler value below because the time integrator and residual construction differ.

### Paper-grid damping validation

Figures 4.19-4.20 use backward Euler. Reproducing that setting gives:

| Problem                         | Best sampled $\eta$ after 15 iterations |
| ------------------------------- | --------------------------------------: |
| heat equation                   |                                   0.500 |
| advection-diffusion, $\nu=0.01$ |                                   0.372 |

![[assets/pint/stmg-figure-4-19.png]]

At fixed $\eta=0.5$, three pre- and post-smoothing steps reduce the error in fewer cycles than one step for the tested heat and advection-diffusion problems.

![[assets/pint/stmg-figure-4-20.png]]

The comparison is in cycle count. A performance decision should instead minimize

$$
\text{work to tolerance}
=
(\text{cost per cycle})\times(\text{number of cycles}),
$$

with communication and memory traffic included.

## Parameter summary

| Method   | Parameter                              | Effect                                                           |
| -------- | -------------------------------------- | ---------------------------------------------------------------- |
| Parareal | number of intervals $N$                | controls concurrency and the maximum finite-step iteration count |
| Parareal | coarse-to-fine step ratio              | trades coarse cost against phase and damping mismatch            |
| PFASST   | collocation nodes and SDC sweeps       | control formal accuracy and work per pipeline stage              |
| MGRiT    | coarsening factor and F/FCF relaxation | control coarse-model quality and smoothing cost                  |
| STMG     | damping $\eta$ and smoothing steps     | control high-frequency reduction and work per cycle              |
