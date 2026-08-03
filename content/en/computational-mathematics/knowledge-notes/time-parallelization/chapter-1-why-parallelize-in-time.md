---
title: "Chapter 1: Why Parallelize in Time?"
description: Motivation, the causal structure of time, four algorithmic lineages and their crossovers, the hyperbolic/parabolic reclassification, and four model PDEs
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-1-why-parallelize-in-time
tags:
  - parallel-in-time
  - computational-mathematics
---

> [!note] Reading scope
> Sections 1.1-1.6 unfold, in order, each cluster of arguments in the paper's abstract and its Section 1 (Introduction), restoring the source's compressed statements into complete causal reasoning. From Section 1.7 onward, these notes add an algebraic interpretation, a leading diagnostic, performance metrics, and a reproduction policy. This keeps the source's line of argument while placing the algorithms of later chapters into a single unified framework.

## 1.1 The central question of the paper

Parallel-in-time (PinT) research asks how to keep increasing concurrency in large-scale evolution problems. Traditional parallelism almost always develops along the spatial direction: the computational mesh is cut into subdomains, and each processor handles one spatial region. This approach has an inherent ceiling—once the subdomains are already small enough, adding more processors makes the genuine floating-point work increasingly uneconomical relative to subdomain-boundary communication, global synchronization, and load imbalance. In other words, purely spatial parallelism has a **strong-scaling saturation point**: beyond some core count, the speedup all but stalls.

Evolution problems still stack up a large number of unknowns along the time direction. For a simulation that runs to a long time window, the product "number of time points × spatial degrees of freedom per time point" is often far larger than the degrees of freedom of a single time level. Since concurrency in the spatial dimension is exhausted, the time dimension becomes the only parallel resource that has not yet been exploited. This is the direct motivation behind the phrase "parallelize in time" in the paper's title.

> [!tip] Insight
> The value of time parallelism is not that it is "faster than spatial parallelism," but that it provides **one additional scalable dimension**. In scenarios where a supercomputer already has millions of cores and spatial parallelism has long saturated, the payoff of PinT is to raise the total concurrency from "number of spatial cores" to "number of spatial cores × number of time segments." Therefore, when evaluating a PinT method, the real question is not "can it converge," but "after spatial parallelism has saturated, how much extra concurrency can it squeeze out."

The difficulty comes from the one-way flow of information in evolution problems: a later state is determined by an earlier one, and not the reverse. Direct time stepping thus forms a causal chain that is inherently sequential. The paper develops precisely this tension: how can several time intervals produce useful computation simultaneously while preserving the correctness of the evolution equation?

The authors first review the four research lineages that formed around algorithmic mechanisms:

1. multiple-shooting methods;
2. domain-decomposition and waveform-relaxation methods;
3. space-time multigrid methods;
4. direct time-parallel methods.

This classification is well suited to tracing the origins of the algorithms, but it is hard to use directly to judge how a method behaves on different dynamics. Recent computational experience shows that many algorithms that perform well on parabolic problems degrade noticeably once they enter hyperbolic or transport-dominated regimes. The paper therefore proposes a second organization, by "which class of problems can be solved effectively": one group of methods remains effective on hyperbolic problems, while the other group mainly exploits the dissipative nature of parabolic dynamics.

Four partial differential equations run through the whole discussion: the heat equation, the advection-diffusion equation, Burgers' equation, and the second-order wave equation. They pass gradually from strong dissipation to persistent propagation—the heat equation is purely parabolic dissipation, the advection-diffusion equation contains both dissipation and transport (and can be tuned continuously between them via the viscosity coefficient), Burgers' equation introduces nonlinear transport and shocks, and the wave equation is undamped, purely hyperbolic propagation. This continuous thread from dissipation to propagation provides a unified physical axis for comparing time-parallel algorithms.

> [!info] Subject classification
> The paper gives the following 2020 Mathematics Subject Classification: Primary 65M55, 65M12, 65M15, 65Y05; Secondary 65M06, 65L10. Here 65M55 corresponds precisely to "multigrid and domain-decomposition methods (evolution equations)," and 65Y05 to "parallel computation," consistent with the theme of the paper.

## 1.2 Why hardware change pushed the time direction to the fore

The paper understands the PinT revival of the past two decades in the context of hardware evolution. Early computers gained much of their performance growth from higher clock frequency (clock rate): the same serial code would automatically run faster on a machine with a higher clock. But the increase in clock frequency hit physical and power constraints in the early 2000s and could no longer be sustained. From then on, single-core performance nearly stalled, and performance growth shifted to relying on **more cores**. This turning point changed "parallelizing an algorithm" from an optional optimization into the only route to greater speed. The paper regards the period around 2004 as the point at which PinT research truly accelerated, precisely because "add cores, not clock" had become the pervasive reality of hardware.

> [!tip] Insight
> The claim in the abstract that "clock frequency peaked around 2004" corresponds to the widely known **failure of Dennard scaling** and the **power wall** in computer architecture: transistors could still keep shrinking (Moore's law did not immediately fail), but leakage and heat made it infeasible to keep trading frequency for performance, so the industry turned to multicore. Understanding this helps grasp the historical inevitability of PinT—it is not a fashionable topic internal to numerical analysis, but the parallelization question that evolution problems must answer once hardware stopped offering "free speedup."

This trend appeared simultaneously at every scale: even small devices such as smartphones have become multicore, while the core counts of high-performance computers reach the order of millions. If an evolution problem is only parallelized spatially, the number of cores it can use effectively is limited by mesh size and spatial communication (see the strong-scaling saturation point in 1.1). To truly keep these cores fed, one must also create concurrency along the time direction.

It is worth emphasizing that the idea of time parallelism predates modern multicore hardware by far. As early as 1964, Nievergelt proposed a prototype method that **trades redundant computation for parallelism**: guess initial values at several interior time points, let the integrations of the individual time segments proceed in parallel, and then stitch and correct these independent trajectories using continuity conditions. He stated explicitly in the original text that the value of such methods lay not in their practicality at the time, but in their potential as a "prototype"—when computers could truly execute many tasks simultaneously, algorithms improved along the same idea would become practically meaningful (Nievergelt 1964). This quotation is therefore often regarded as the historical starting point of the idea of time parallelism: the correct algorithmic idea appeared first, waiting for the hardware to catch up forty years later.

In the decades that followed, a variety of time-parallel methods appeared. After Lions, Maday, and Turinici proposed **Parareal** in 2001, the field entered a phase of rapid development—Parareal uses a cheap coarse propagator for global prediction and an expensive fine propagator for parallel correction, is concise in form and easy to analyze, and became the trigger for modern PinT research. The paper also points to three survey resources for readers with different needs: Gander (2015) traces the historical development, Ong and Schröder (2020) emphasize engineering applications, and Gander and Lunet (2024) is a systematic research monograph.

> [!info] Background: the cross-field origins of waveform relaxation and Schwarz
> Among the four lineages discussed below, waveform relaxation was not originally proposed for evolution PDEs but was invented in the 1980s for circuit simulation in very-large-scale integration (VLSI) (Lelarasmee, Ruehli, and Sangiovanni-Vincentelli 1982); the idea of domain decomposition can be traced further back to Schwarz's (1870) alternating method used to prove existence. Many of the tools of time parallelism were "borrowed" from other fields and reinterpreted, which also explains why its classification boundaries are inherently blurred (see 1.5).

## 1.3 How causality forms a sequential computation chain

At first glance, parallelizing evolution problems seems impossible because of the **causality principle**: the solution at a later time is determined by earlier times, not the reverse. To make this obstacle clear, the paper uses the simplest ordinary differential equation and its forward Euler discretization as an illustration:

$$
\partial_t u=f(u),\qquad u(0)=u_0,\qquad
u_{n+1}=u_n+\Delta t\,f(u_n). \tag{1.1}
$$

The recurrence directly exposes the **data dependency**. To compute $u_{n+1}$, one must first know $u_n$—because the right-hand side $f(u_n)$ depends on the current state. Stringing this dependency along the time axis yields an unbreakable chain: $u_0\to u_1\to\cdots\to u_{12}$. The paper points out the consequence with a concrete instance: if $u_9$ has not yet been computed, then even with idle processors available, one cannot accurately obtain $u_{10}$, $u_{11}$, $u_{12}$ from the same forward Euler recurrence. This is exactly the sequentiality that Figure 1.1 emphasizes.

![Original Figure 1.1: sequential dependence in forward Euler time integration](assets/papers/time-parallelization/source-figures/figure-1-1.svg)

**Reading Figure 1.1.** The horizontal axis lists $t_0$ through $t_{12}$, and each node on the polyline represents the numerical solution $u_n$ at the corresponding discrete time. The segments connecting neighboring nodes draw the recurrence dependency as one continuous path: computation can only advance point by point in the order $u_0,u_1,\ldots,u_{12}$, and each step must wait for the previous one to finish. The shape of the curve in the figure does not express the exact solution of any particular equation; it uses a single schematic trajectory only to highlight the structural fact that "time stepping is ordered."

It must be emphasized that this sequential chain lies at the **level of data dependency**, not at the level of numerical stability. Even if forward Euler is replaced by any other one-step or multistep scheme, as long as it is defined in the manner of "deriving the unknown future from the known past," it carries the same causal chain. The challenge of time parallelism is therefore not a defect of any particular discretization scheme, but the structure of the evolution problem itself.

> [!tip] Insight
> Figure 1.1 can also be read in reverse: if we **allow temporarily using wrong initial values** to launch the integration of segments 10, 11, and 12, then the computation of those segments no longer has to wait for $u_9$ and can begin in parallel immediately—at the cost that those segments start out "wrong" and need subsequent iterations to correct the interface values round by round until they are self-consistent. Almost all PinT methods do the same thing: they replace "compute correctly once, sequentially" with "first compute in parallel but wrong, then correct iteratively." This is precisely the modern form of Nievergelt's principle of "trading redundancy for parallelism."

PinT algorithms preserve the causal constraint of the original problem while changing the computational way in which that constraint is satisfied. Common approaches include: providing provisional values at the time interfaces and correcting them iteratively (shooting/Parareal), exchanging waveforms over the entire time window (waveform relaxation), building coarse time levels to propagate long-range information (multigrid), or transforming the fully time-coupled problem into a set of subproblems that can be solved concurrently (diagonalization). All of these introduce extra computation, iteration, and communication—precisely the price paid to buy time concurrency, in the same spirit as Nievergelt's judgment sixty years ago.

## 1.4 Four historical lineages based on algorithmic mechanism

The paper first follows the traditional four-way division. The common task of these four classes is to **break the long recurrence chain of (1.1)**; the difference lies in the mechanism used to break it. Understanding the historical starting point of each lineage helps to see which body of convergence theory each one inherits.

| Lineage                                      | Historical origin and representative methods                                                                     | How coupling across time is organized                                                                                     |
| -------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| Multiple shooting                            | Nievergelt (1964), Bellen-Zennaro (1989), Chartier-Philippe (1993), Saha-Stadel-Tremaine (1997), Parareal (2001) | Split the long interval into initial-value subproblems and correct each segment's initial data via continuity conditions  |
| Domain decomposition and waveform relaxation | Schwarz decomposition (1870), waveform relaxation (1982), Bjørhus (1995), SWR (1999)                             | Solve on space-time subdomains and exchange whole time waveforms in overlaps or on interfaces                             |
| Multigrid                                    | Parabolic multigrid (Hackbusch, 1984), STMG (Gander-Neumüller, 2016)                                             | Smooth error and apply coarse-grid correction between levels of a space-time mesh                                         |
| Direct time parallelism                      | Parallel time stepping (Miranker-Liniger, 1967), RIDC (2010), ParaExp (2013), time diagonalization and ParaDiag  | Solve concurrently using deferred correction, matrix exponentials, or the diagonalizable structure of the temporal matrix |

**(1) Multiple shooting.** This line splits one long initial-value problem into several short ones, so that the segments can be integrated simultaneously, while the initial value at the start of each segment must be corrected round by round through a global continuity condition. Its conceptual source is Nievergelt's (1964) parallel guessing of initial values along the time axis; it was then developed by Bellen-Zennaro (1989) and Chartier-Philippe (1993) into more systematic time-parallel/multiple-shooting schemes, Saha-Stadel-Tremaine (1997) applied a similar idea to long-time integration in celestial mechanics, and it finally culminated in the Parareal algorithm of Lions et al. (2001), which spawned a large number of variants.

> [!info] Background: the predictor-corrector structure of Parareal
> The Parareal method proposed by Lions, Maday, and Turinici (2001) can be written as a one-step predictor-corrector iteration
> $$U_{n+1}^{k+1}=\mathcal{G}(U_n^{k+1})+\mathcal{F}(U_n^{k})-\mathcal{G}(U_n^{k}),$$
> where $\mathcal{F}$ is the expensive but accurate fine propagator (evaluable in parallel over the time segments) and $\mathcal{G}$ is the cheap coarse propagator (which propagates global information sequentially). When the iteration count $k$ reaches the number of time segments, it reduces to the exact serial solution, so effective speedup requires early stopping of the iteration (see 1.9). This structure is also the bridge for viewing Parareal as "multigrid coarsening" or "multiple shooting" (see 1.5).

**(2) Domain decomposition and waveform relaxation (SWR).** This line combines spatial domain decomposition (Schwarz 1870) with waveform relaxation (Lelarasmee et al. 1982)—the two were merged in Bjørhus (1995), forming Schwarz waveform relaxation (Gander, Halpern, and Nataf 1999). Its distinctive feature is that each subdomain **solves the complete time window at once**, and what is exchanged during the iteration is a whole time function (waveform) on the interface, rather than a scalar at each individual time step. This brings two benefits: first, the communication granularity is coarser, well suited to packing whole-segment time communication; second, it naturally matches the **finite propagation speed** of hyperbolic equations—information can only propagate a finite spatial distance within one iteration step, so the convergence of interface exchange can be analyzed in terms of propagation speed.

**(3) Space-time multigrid.** This line places all time points and space points into a single unified space-time mesh hierarchy, using a **smoother** to damp local high-frequency error and **coarse-grid correction** to propagate long-range low-frequency information. It goes back to Hackbusch's (1984) parabolic multigrid and developed into the fully parallel space-time multigrid STMG (Gander and Neumüller 2016). MGRiT, to be discussed in Chapter 4, is closely related to this idea.

**(4) Direct (non-iterative) time parallelism.** This line performs no outer iteration but directly exploits the special algebraic structure of the discrete system to solve it in parallel all at once. It began with the parallel time stepping of Miranker and Liniger (1967); its modern representatives are revisionist integral deferred correction RIDC (Christlieb, Macdonald, and Ong 2010, arranging deferred corrections into a pipeline), ParaExp (Gander and Güttel 2013, separating the inhomogeneous response from homogeneous propagation and handling the latter in parallel with matrix exponentials), and diagonalization-based parallelism (Maday and Rønquist 2008), the last of which further developed into ParaDiag (Gander et al. 2021c, using the diagonalizable structure of the time-stepping matrix and an FFT to decouple the all-time system into spatial systems that can be solved concurrently).

> [!tip] Insight
> The four lineages appear to go their separate ways, but they all answer the same algebraic question: sequential time stepping is equivalent to forward substitution on a **block lower-triangular** all-time system (see 1.7). Multiple shooting blocks it by time segment and iteratively fills in the interfaces; waveform relaxation blocks it by space-time subdomain and exchanges interface waveforms; multigrid approximates its inverse with a hierarchy; and direct methods look for approximations that can be quickly diagonalized/factorized. Once "breaking the causal chain" is translated into "approximately solving the inverse of a block lower-triangular system," all four classes fall onto the same map.

## 1.5 Why the four classes cross over into one another

The paper points out in particular that the boundaries of the traditional four classes are **not strict**. The first three classes usually appear in iterative form, and the fourth class was initially dominated by direct (non-iterative) methods—but this "iterative vs. direct" line was quickly broken. The paper uses ParaDiag as the most striking example and lists a series of interwoven connections:

- ParaDiag was originally just a **direct** solver constructed by diagonalizing the time-stepping matrix; **iterative** forms soon appeared—both a version embedded in second-class waveform relaxation (Gander and Wu 2019) and a version embedded in first-class Parareal (Gander and Wu 2020);
- the **approximate ParaDiag** operator can also serve as a stationary iteration, or as a preconditioner acting directly on the all-at-once system arising from the third-class space-time discretization, which is then handed to a Krylov method (McDonald, Pestana, and Wathen 2018; Liu and Wu 2020);
- Parareal can be interpreted as a multigrid method using **aggressive coarsening** along the time direction (Gander and Vandewalle 2007), connecting the first class with the third;
- conversely, MGRiT can be understood as an **overlapping** generalization of Parareal (Gander, Kwok, and Zhang 2018b).

These connections show that the same algorithm can simultaneously carry the multiple identities of shooting, waveform relaxation, multigrid, and all-time algebraic solving. Classifying by technical origin certainly helps trace the development of ideas, but a method's **actual behavior** in the face of propagation, dissipation, and nonlinearity cannot be read off directly from its "ancestry" and must be analyzed separately.

> [!tip] Insight
> The practical value of the crossover classification is that **convergence theories can be borrowed from one another**: since Parareal is a multigrid with aggressive coarsening, one can use multigrid's two-grid local Fourier analysis to predict its convergence factor; since ParaDiag has a waveform-relaxation form, one can use waveform relaxation's propagation-speed analysis to understand its behavior on hyperbolic problems. This also foreshadows a hidden thread running through the whole book—it is precisely these "identity swaps" that make it clear that what really decides success or failure is not which class an algorithm belongs to, but whether the underlying dynamics is dissipative or propagative (see 1.6, 1.8).

## 1.6 The two-way classification adopted by the paper

Based on the observations above, the paper abandons the old "by technical origin" map and instead organizes Sections 3 and 4 "by which class of problems can be solved effectively." The core judgment behind this reorganization is that **parabolic and hyperbolic problems place essentially different demands on time parallelism**, and the old classification precisely obscured this.

### Methods that remain effective for hyperbolic problems (Section 3)

This group contains four methods: Schwarz waveform relaxation (and its connection to tent pitching), parallel integral deferred correction, ParaExp, and ParaDiag. Their common feature is that they **can handle persistently propagating information**—either by explicitly respecting the finite propagation speed of the characteristics (SWR/tent pitching) or by bypassing a "coarse dissipative propagator" and relying on global algebraic structure (the matrix exponential of ParaExp, the diagonalization of ParaDiag). The paper emphasizes that this group of methods **also works very well on parabolic problems**, so their range of applicability is broader.

### Methods designed primarily for parabolic problems (Section 4)

This group contains Parareal, PFASST, MGRiT, and STMG. They fully exploit the **dissipation** of parabolic problems—error decays rapidly over time and scale separation is pronounced—so that a cheap coarse propagator/coarse-level correction is enough to drive convergence, and they are often highly efficient on problems such as the heat equation. But once one enters the hyperbolic limit, dissipation disappears, coarse propagation and coarse-level correction struggle to preserve **phase and propagation paths**, and convergence slows substantially or even exhibits transient error growth (see 1.8).

> [!tip] Insight
> The two classifications are not a replacement relation but **two orthogonal axes**: the horizontal axis is "which mechanism is used to break the causal chain" (the four classes of 1.4), and the vertical axis is "does it rely on dissipation or respect propagation" (the two groups of 1.6). The paper chooses the vertical axis to organize its main text because the vertical axis is the line that decides whether a method can cross the hyperbolic/parabolic boundary; the horizontal axis is more about differences at the historical and implementation level. Reading the later chapters with these two axes in mind lets one simultaneously answer "how does it parallelize" and "why is it (in)effective on this class of problems."

Section 2 first uses four model problems to explain physically the difference between parabolic and hyperbolic dynamics; Section 3 introduces the first group (hyperbolic-effective) methods, which are generally also better on parabolic problems; Section 4 discusses the second group (parabolic-designed) methods, which generally run into difficulty on hyperbolic problems; and Section 5 gives the conclusions. The MATLAB code for Sections 2-4 of the paper is publicly available at [wushulin/ActaPinT](https://github.com/wushulin/ActaPinT). The Python conversion and extended version used on this site is at [freezeng123456/ActaPinT-Python](https://github.com/freezeng123456/ActaPinT-Python).

## 1.7 From a recurrence to an all-at-once system

The following adds one algebraic layer. It does not change the causal relation of equation (1.1); it merely writes the unknowns at all time points simultaneously into a single linear system—this is exactly the concrete form of the "unified map" at the end of 1.4.

For a linear one-step discretization, let one step of advancement be given by a propagation matrix $\Phi$:

$$
u_{n+1}-\Phi u_n=g_{n+1}.
$$

Stacking all temporal unknowns in time order yields the block lower-triangular linear system

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

Sequential time stepping is exactly equivalent to performing **exact forward substitution** on this block lower-triangular system—solving row by row downward is precisely the causal chain of (1.1). The mainstream PinT method families can therefore be understood uniformly as: constructing different **parallel approximate solvers** for this all-at-once system.

- Parareal and MGRiT approximate its inverse with low-cost coarse propagators and coarse temporal levels;
- STMG constructs multilevel corrections on the complete space-time mesh;
- ParaDiag approximates the original system with a circulant or nearly circulant temporal matrix that can be diagonalized by FFT, thereby replacing the block lower-triangular solve with a batch of spatial solves that can run concurrently;
- SWR solves on space-time subdomains and exchanges interface waveforms, which amounts to iterating this system in a block Jacobi/Gauss-Seidel manner;
- deferred-correction methods arrange the correction process into a pipeline that can execute with overlap.

This viewpoint clarifies a frequently misunderstood fact: **the concurrent work of time parallelism takes place inside the iteration, decomposition, transform, or correction pipeline**, while the final solution still exactly satisfies the time relations prescribed by the discrete evolution equation—the block lower-triangular structure itself is not destroyed; only its solution is approximated in parallel.

## 1.8 Why hyperbolic problems are harder: a leading diagnostic

The reorganization of 1.6 raises "dissipation vs. propagation" to a central position. Here, starting from the propagation matrix $\Phi$ of 1.7, we give a quantifiable leading diagnostic that sets the stage for the physical discussion in Chapter 2.

For **dissipative** propagators (such as the discretization of the heat equation), many error modes satisfy $|\lambda(\Phi)|<1$, i.e. the spectral radius rapidly compresses the corresponding components. Rapidly decaying components have negligible influence on long time scales, so even if the coarse temporal model cannot represent them accurately, it will not significantly pollute the long-time solution—**the coarse propagator only needs to capture the slow, decaying, large-scale components**, which is precisely the source of the efficiency of Parareal/MGRiT/STMG on parabolic problems.

For **hyperbolic/transport** problems (such as wave propagation or pure advection), the eigenvalues of the important modes are closer to the **unit circle** $|\lambda(\Phi)|\approx 1$: energy neither decays nor grows, and error mainly manifests as a **phase shift** or a **deviation in the propagation position**. In this case, even a very small phase mismatch between the fine and coarse propagators will not be smoothed out by dissipation but instead **accumulates** across many time intervals, destroying convergence.

> [!tip] Insight
> In one sentence: parabolic problems forgive a coarse propagator (error is automatically absorbed by $|\lambda|<1$), while hyperbolic problems do not (error accumulates losslessly on $|\lambda|\approx 1$). Therefore hyperbolic-effective methods must explicitly "align phase/characteristics" (SWR, tent pitching) or "bypass the coarse propagator" (the matrix exponential of ParaExp, the diagonalization of ParaDiag), rather than relying, as parabolic methods do, on a cheap but inaccurate coarse propagator.

From this, one can anticipate three conclusions that will be verified later:

1. as the viscosity $\nu$ of the advection-diffusion equation decreases and it gradually approaches pure advection, standard temporal coarsening becomes increasingly difficult;
2. Parareal in transport-dominated regimes may experience transient error growth (the convergence curve rises before it falls);
3. methods that respect characteristic transport, explicit phase correction, waveform relaxation, and global diagonalization deserve priority.

Chapter 2 will, through the heat equation, the advection-diffusion equation, Burgers' equation, and the wave equation, ground this judgment—moving from spectral-radius intuition to concrete physical and numerical phenomena.

## 1.9 How to judge the actual parallel payoff

A convergence curve can only answer "how does algorithmic error decrease with iteration," but cannot alone answer "how much faster is it, really." The true parallel payoff also depends on the combined effect of coarse and fine solve costs, communication, memory, and scalability. The table below lists the questions that should each be answered when evaluating a PinT method:

| Metric               | Question to answer                                                                              |
| -------------------- | ----------------------------------------------------------------------------------------------- |
| Convergence          | How many global iterations are needed to reach the discretization-error scale?                  |
| Coarse-to-fine cost  | Is the coarse propagator substantially cheaper than the fine propagator?                        |
| Temporal concurrency | How many fine or local solves can execute simultaneously?                                       |
| Communication        | Does each iteration use neighbor exchange or a global collective?                               |
| Memory               | Must all temporal states be stored?                                                             |
| Parameter robustness | How does convergence change as the dynamics approaches the hyperbolic limit?                    |
| Scalability          | As the number of time intervals grows, do iteration counts and communication remain controlled? |

For Parareal with $N$ time intervals, an idealized speedup bound that neglects communication overhead is

$$
S\lesssim \frac{NC_F}{K(C_F+NC_G)}.
$$

Here $C_F$ and $C_G$ are the costs of the fine and coarse propagators respectively, and $K$ is the number of global iterations needed to reach the target accuracy. To obtain meaningful speedup, one must satisfy both $K\ll N$ (early stopping of the iteration) and $C_G\ll C_F$ (the coarse propagator is cheap enough).

> [!tip] Insight
> Reading 1.8 together with this formula makes clear where hyperbolic problems are "slow": phase-error accumulation makes the $K$ needed for convergence larger, directly depressing the speedup in the denominator; and to suppress the phase mismatch one often has to use a more accurate (more expensive) coarse propagator, raising $C_G$. Both factors push $S$ in the wrong direction at once. This explains why "parabolic-designed" methods on hyperbolic problems are not merely a little slower but may lose speedup entirely—convergence itself does not guarantee that $K\ll N$ and $C_G\ll C_F$ hold.

## 1.10 The reproduction policy of these notes

The companion project provides two execution levels:

```bash
python3 run_experiments.py all --quick --output-dir results/quick
python3 run_experiments.py all --output-dir results/formal
```

The quick configuration shrinks the mesh size and is used to quickly check that the code paths and data interfaces are working. The figures and numerical values shown in Chapters 2-4 come from the formal configuration. Each experiment saves SVG/PNG figures together with a JSON record of parameters and metrics, making it easy to trace the exact numerical settings behind each figure and ensuring the results are reproducible.

## 1.11 Source-coverage audit

| Source location                                | Corresponding sections | Material covered                                                                                                                                             |
| ---------------------------------------------- | ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Abstract, pp. 385-386                          | 1.1 and 1.6            | Definition of PinT, the 1964/2004 timeline, the traditional four classes, the hyperbolic difficulty, the two new groups, four model PDEs, MSC classification |
| Introduction, paragraphs 1-2, pp. 386-387      | 1.2                    | Multicore hardware background and the power wall, Nievergelt's "redundancy for parallelism," the Parareal trigger point, and three types of surveys          |
| Causality paragraph, p. 387                    | 1.3                    | Forward Euler equation (1.1), the $u_9\to u_{10}$ data dependency, the original Figure 1.1 and its figure-level reading                                      |
| Four method families, pp. 387-388              | 1.4                    | The complete lineages of multiple shooting, SWR, STMG, and direct methods (RIDC/ParaExp/ParaDiag), representative references, and coupling mechanisms        |
| Classification crossover paragraph, p. 388     | 1.5                    | The loosening of iterative vs. direct, ParaDiag's direct/iterative/preconditioner three forms, Parareal↔multigrid, and the MGRiT↔Parareal connection         |
| Organization of the paper, p. 388              | 1.6                    | The four hyperbolic-effective methods, the four parabolic-designed methods, the section plan, and the original MATLAB code repository                        |
| Model-problem sentence in the abstract, p. 386 | 1.1 and 1.8            | The four model PDEs and their continuous thread from dissipation to propagation                                                                              |

## Summary

Section 1 of the paper accomplishes three things: it explains the research motivation from hardware parallelism (the power wall, multicore, saturation of spatial parallelism); it uses the data dependency of forward Euler to show why the temporal causal chain is inherently sequential; and it shifts the material from the traditional four classes "by technical origin" to a two-way classification "by problem nature," pointing out that the crossover classification lets convergence theories be borrowed from one another. The core judgment running through the whole book is: how an algorithm transmits information across time, and whether this mechanism can preserve the key dynamical features of dissipation, phase, and finite propagation speed—this is also the true watershed distinguishing the "parabolic-designed" group of methods from the "hyperbolic-effective" one.

## Primary source

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), pp. 385-489, abstract and Section 1.
