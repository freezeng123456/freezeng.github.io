---
title: "Chapter 1: Why Parallelize in Time?"
description: Motivation, causal structure, method lineages, and the hyperbolic/parabolic reclassification
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-1-why-parallelize-in-time
tags:
  - parallel-in-time
  - computational-mathematics
---

> [!note] Reading scope
> Sections 1.1-1.6 follow the argument of the paper's abstract and Section 1 in sequence. Sections 1.7 onward add an algebraic interpretation, performance criteria, and the reproduction policy used in these notes. This separation preserves the source argument while placing the later algorithms in a common framework.

## 1.1 The central question

Parallel-in-time (PinT) methods seek additional concurrency for large evolution problems. Conventional solvers first partition the spatial domain. Once each spatial subdomain becomes small, communication, synchronization, and load imbalance consume a growing fraction of the run time. The temporal grid still contains many unknowns and therefore offers another direction for parallel work.

Evolution has a one-way flow of information: later states are determined by earlier states. Direct time stepping consequently forms a causal chain. The paper asks how several time intervals can perform useful work concurrently while the final solution continues to satisfy the evolution equation.

The literature is traditionally organized into four algorithmic lineages:

1. multiple-shooting methods;
2. domain-decomposition and waveform-relaxation methods;
3. space-time multigrid methods;
4. direct time-parallel methods.

This classification records how the algorithms were constructed. It is less informative about their behavior on different dynamics. Recent numerical experience has shown that several methods that work well for parabolic problems deteriorate on hyperbolic or advection-dominated problems. The paper therefore reorganizes the subject by problem class: one group remains effective for hyperbolic dynamics, while the other was designed around parabolic dissipation.

Four PDEs provide a common set of examples: the heat equation, the advection-diffusion equation, Burgers' equation, and the second-order wave equation. Together they form a progression from strong dissipation to persistent propagation.

> [!info] Subject classification
> The paper lists the following 2020 Mathematics Subject Classification: Primary 65M55, 65M12, 65M15, 65Y05; Secondary 65M06, 65L10.

## 1.2 Why hardware made the time direction important

The paper places the rapid growth of PinT research in the context of computer architecture. Processor clock rates approached physical and power limits, so further performance increasingly came from additional cores. The abstract identifies 2004 as the point at which research activity accelerated. Multicore designs became standard even in small devices such as smartphones, while high-performance systems grew to very large core counts. A solver that exploits only spatial parallelism can use those cores effectively only while the spatial problem provides enough local work.

The idea of parallel time integration predates modern multicore hardware. In 1964, Nievergelt proposed a prototype that exchanged redundant computation for concurrency. His conclusion emphasized the principle behind the construction and anticipated its value once machines capable of many simultaneous computations became available.

Several methods followed over the next decades. The Parareal algorithm of Lions, Maday, and Turinici (2001) brought time parallelism to the center of numerical research. The paper points readers to three complementary accounts: Gander (2015) for the historical development, Ong and Schröder (2020) for applications, and Gander and Lunet (2024) for a systematic monograph.

## 1.3 How causality creates a sequential chain

Consider an ordinary differential equation and its forward Euler discretization:

$$
\partial_t u=f(u),\qquad u(0)=u_0,\qquad
u_{n+1}=u_n+\Delta t\,f(u_n). \tag{1.1}
$$

The recurrence exposes the data dependency directly. The value $u_n$ must be available before $u_{n+1}$ can be computed. If $u_9$ is still unknown, an idle processor cannot use the same forward Euler recurrence to compute accurate values of $u_{10}$, $u_{11}$, and $u_{12}$. This is the sequential structure illustrated by Figure 1.1 of the paper.

![The one-way causal chain in forward Euler time integration](assets/diagrams/pint/en/sequential-time-stepping.svg)

PinT algorithms preserve the causal constraint of the original problem while changing how that constraint is enforced computationally. They may assign provisional interface values and correct them globally, exchange complete interface waveforms, introduce a hierarchy of temporal grids, or transform the coupled time problem into concurrent subproblems. Redundant work, iteration, and communication are the costs paid for this additional concurrency.

## 1.4 Four historical lineages based on algorithmic mechanism

The paper first reviews the conventional four-way classification. Each lineage uses a different construction to break one long recurrence into coupled pieces.

| Lineage                                      | Origins and representative methods                                                                               | Treatment of coupling across time                                                                              |
| -------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| Multiple shooting                            | Nievergelt (1964), Bellen-Zennaro (1989), Chartier-Philippe (1993), Saha-Stadel-Tremaine (1997), Parareal (2001) | Split the interval into initial-value subproblems and correct their initial data through continuity conditions |
| Domain decomposition and waveform relaxation | Schwarz decomposition (1870), waveform relaxation (1982), Bjørhus (1995), SWR (1999)                             | Solve on space-time subdomains and exchange time-dependent interface data                                      |
| Multigrid                                    | Parabolic multigrid (Hackbusch, 1984), STMG (Gander-Neumüller, 2016)                                             | Smooth error and apply coarse-grid correction across a space-time hierarchy                                    |
| Direct time parallelism                      | Parallel time stepping (Miranker-Liniger, 1967), RIDC (2010), ParaExp (2013), time diagonalization and ParaDiag  | Exploit deferred correction, matrix exponentials, or diagonalizable temporal structure                         |

Multiple shooting replaces one long initial-value problem with several short ones. The subproblems can run concurrently, while their starting values are corrected until the interfaces become continuous. Parareal is the most influential modern representative of this lineage.

SWR combines spatial domain decomposition with waveform relaxation. Each subdomain solves over a complete time window, and the iteration exchanges time functions along interfaces. This organization has a different communication granularity from step-by-step exchange and can exploit finite propagation speed in hyperbolic equations.

Space-time multigrid places all spatial and temporal unknowns in one hierarchy. Smoothers reduce local error, and coarse levels transmit long-range information. STMG follows this approach, and MGRiT, discussed in Chapter 4, is closely related to the same multilevel viewpoint.

Direct time-parallel methods use special algebraic structure in the discretization. RIDC forms a pipeline of deferred corrections, ParaExp separates inhomogeneous response from homogeneous propagation, and ParaDiag diagonalizes a temporal matrix so that shifted spatial systems can be solved concurrently.

## 1.5 Why the four classes overlap

The boundaries of the conventional classification are deliberately loose. The paper identifies several important connections:

- the first three lineages usually produce iterative methods, while the fourth initially emphasized direct methods;
- ParaDiag began as a direct solver based on diagonalizing the time-stepping matrix, then acquired waveform-relaxation and Parareal-based iterative variants;
- an approximate ParaDiag operator can define a stationary iteration or precondition the all-at-once system for a Krylov solver;
- Parareal can be interpreted as a temporal multigrid method with aggressive coarsening;
- MGRiT can also be viewed as an overlapping generalization of Parareal.

One algorithm may therefore have features of shooting, waveform relaxation, multigrid, and all-at-once linear algebra. Classification by technical ancestry is useful for tracing ideas. Its performance under propagation, dissipation, and nonlinearity requires a separate analysis.

## 1.6 The paper's two-way classification by problem class

Sections 3 and 4 are organized by the type of dynamics that a method can handle effectively.

### Methods that remain effective for hyperbolic problems

This group contains Schwarz waveform relaxation and its relation to tent pitching, parallel integral deferred correction, ParaExp, and ParaDiag. These methods generally work for parabolic problems as well. Their constructions transmit propagative information or reduce dependence on a dissipative coarse propagator through global algebraic structure.

### Methods designed primarily for parabolic problems

This group contains Parareal, PFASST, MGRiT, and STMG. These algorithms exploit error damping and scale separation and can be highly effective on problems such as the heat equation. Near a hyperbolic limit, coarse propagation and coarse-grid correction have difficulty preserving phase and propagation paths, and convergence may slow substantially.

Section 2 first develops the distinction through four model problems. Section 3 studies the first group, Section 4 studies the second, and Section 5 summarizes the conclusions. MATLAB code for the experiments in Sections 2-4 is available in [wushulin/ActaPinT](https://github.com/wushulin/ActaPinT). The Python conversion and extensions used here are maintained in [freezeng123456/ActaPinT-Python](https://github.com/freezeng123456/ActaPinT-Python).

## 1.7 From a recurrence to an all-at-once system

The following algebraic view supplements the introduction. It preserves the causal relation in equation (1.1) and writes all temporal unknowns in a single system.

For a linear one-step discretization, let

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

Sequential time stepping is exact forward substitution on this block lower-triangular system. The main PinT families construct different parallel approximations:

- Parareal and MGRiT use inexpensive coarse propagation and coarse temporal levels;
- STMG applies multilevel correction on the complete space-time grid;
- ParaDiag approximates the temporal matrix by a circulant or nearly circulant matrix that is diagonalizable by FFT;
- SWR solves on space-time subdomains and exchanges interface waveforms;
- deferred-correction methods organize corrections into an overlapping pipeline.

Parallel work takes place inside an iteration, decomposition, transform, or correction pipeline. The converged result still satisfies the temporal relations imposed by the discretized evolution equation.

## 1.8 Why hyperbolic problems are harder: a first diagnostic

For a dissipative propagator, many error modes satisfy $|\lambda(\Phi)|<1$. Rapidly damped components have little long-time influence, so a coarse temporal model can remain useful even when it represents those components inaccurately.

Important modes of wave and transport problems stay closer to the unit circle. Their error often appears as a phase shift or a displacement of the propagated profile. A small phase mismatch between fine and coarse propagation can accumulate across many intervals. This observation suggests three expectations:

1. standard temporal coarsening becomes harder as the viscosity $\nu$ decreases in advection-diffusion;
2. Parareal may experience transient error growth in transport-dominated regimes;
3. characteristic transmission, phase correction, waveform relaxation, and global diagonalization deserve priority in such regimes.

Chapter 2 develops this diagnostic through the heat, advection-diffusion, Burgers, and wave equations.

## 1.9 Criteria for useful parallelism

A convergence curve describes the reduction of algorithmic error. Parallel efficiency also depends on coarse and fine costs, communication, memory, and scalability.

| Quantity             | Question                                                                                  |
| -------------------- | ----------------------------------------------------------------------------------------- |
| Convergence          | How many global iterations reach the discretization-error scale?                          |
| Coarse-to-fine cost  | Is the coarse propagator substantially cheaper than the fine propagator?                  |
| Temporal concurrency | How many fine or local solves can execute simultaneously?                                 |
| Communication        | Does each iteration use neighbor exchange or a global collective?                         |
| Memory               | Must the solver store all temporal states?                                                |
| Parameter robustness | How does convergence change near the hyperbolic limit?                                    |
| Scalability          | Do iteration counts and communication remain controlled as the number of intervals grows? |

For Parareal with $N$ time intervals, an idealized speedup bound that neglects communication is

$$
S\lesssim \frac{NC_F}{K(C_F+NC_G)}.
$$

Here $C_F$ and $C_G$ are the fine and coarse propagation costs, and $K$ is the number of global iterations. Useful speedup requires $K\ll N$ and $C_G\ll C_F$. Convergence alone does not ensure either condition.

## 1.10 Reproduction policy used in these notes

The companion project provides two execution levels:

```bash
python3 run_experiments.py all --quick --output-dir results/quick
python3 run_experiments.py all --output-dir results/formal
```

The quick configuration uses smaller grids to test code paths and data interfaces. Figures and numerical values in Chapters 2-4 come from the formal configuration. Each experiment stores SVG and PNG figures together with a JSON record of parameters and metrics.

## 1.11 Source-coverage audit

| Source location                           | Corresponding sections | Material covered                                                                                                                                    |
| ----------------------------------------- | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| Abstract, pp. 385-386                     | 1.1 and 1.6            | Definition of PinT, the 1964/2004 timeline, four conventional classes, the hyperbolic difficulty, the two-way reclassification, and four model PDEs |
| Introduction, paragraphs 1-2, pp. 386-387 | 1.2                    | Multicore hardware, Nievergelt's redundancy principle, the post-Parareal expansion, and three complementary surveys                                 |
| Causality paragraph, p. 387               | 1.3                    | Forward Euler equation (1.1), the dependency from $u_9$ to $u_{10}$, and the meaning of Figure 1.1                                                  |
| Four method families, pp. 387-388         | 1.4                    | Historical origins and coupling mechanisms of multiple shooting, SWR, STMG, RIDC, ParaExp, and ParaDiag                                             |
| Overlapping classifications, p. 388       | 1.5                    | Direct and iterative ParaDiag, all-at-once preconditioning, Parareal as multigrid, and the MGRiT-Parareal connection                                |
| Organization of the paper, p. 388         | 1.6                    | Hyperbolic-effective methods, parabolic-designed methods, the section plan, and the original code repository                                        |

## Summary

Section 1 of the paper establishes three foundations. It connects the motivation for PinT to hardware concurrency, uses forward Euler to expose the causal chain, and replaces the conventional four-way taxonomy with a two-way organization by problem class. The later chapters repeatedly ask how an algorithm communicates information across time and whether that mechanism preserves dissipation, phase, and finite propagation speed.

## Primary source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), pp. 385-489, abstract and Section 1.
