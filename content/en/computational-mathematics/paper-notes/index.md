---
title: Close Readings of Tao Zhou's Work
description: Bilingual paper notes organising the publication list on Tao Zhou's homepage into seven topics
lang: en
translation: computational-mathematics/paper-notes
tags:
  - computational-mathematics
  - paper-notes
  - uncertainty-quantification
---

These notes organise all 108 items listed on the homepage of Tao Zhou (Institute of Computational Mathematics and Scientific/Engineering Computing, Chinese Academy of Sciences): 95 published papers and 13 submissions or preprints. The numbering and topic assignment follow the [homepage publication list](https://lsec.cc.ac.cn/~tzhou/).

> [!info] What these pages are and are not
> They are **third-party reading notes**, not a restatement of the originals and not the authors' own words. Each page gives the problem setting, the construction, the main results and their qualifications, and closes with the source reference; proofs belong to the papers. Technical content is assembled from public abstracts, public preprints and public full texts. Anything that could not be checked against public material is marked as unverified rather than filled in by guesswork.

![Seven topics on one technical spine](assets/diagrams/tao-zhou-papers/en/research-map.svg)

## Why read by topic instead of by year

The homepage list is in reverse chronological order, which is convenient for lookup and inconvenient for understanding. Reading by date forces repeated jumps between unrelated problems: discrete least squares in 2014, parareal convergence analysis in 2015, multistep schemes for forward-backward stochastic differential equations also in 2015. Each belongs to a different technical tradition.

Reading by topic exposes something else: the same technical move recurs across problems. Write a sequential time recurrence as one all-at-once operator, then build a parallel approximate inverse for it. Write a high-dimensional integral over random inputs as a weighted least-squares problem, then design the sampling density for it. Write a PDE residual as a martingale property along paths, then train two networks against each other. These belong nominally to parallel-in-time methods, uncertainty quantification and scientific machine learning, but they share one strategy: **express the intractable structure explicitly as one operator equation, design a weighted and adaptive approximate inversion for it, and let an error indicator decide where true information gets added.**

## The seven topics

| Topic                                                                                                                 | Papers | Central question                                                           | Governing quantity                                   |
| --------------------------------------------------------------------------------------------------------------------- | ------ | -------------------------------------------------------------------------- | ---------------------------------------------------- |
| [[en/computational-mathematics/paper-notes/stochastic-approximation/index\|Stochastic approximation and collocation]] | 23     | build a stable polynomial surrogate from as few forward solves as possible | Christoffel function, induced sampling density       |
| [[en/computational-mathematics/paper-notes/spectral-and-reduced-order/index\|Spectral and reduced-order methods]]     | 6      | keep spectral accuracy on unbounded domains and low-rank manifolds         | far-field decay rate, smallest singular value        |
| [[en/computational-mathematics/paper-notes/bayesian-inference/index\|Bayesian inverse problems]]                      | 10     | sample a posterior under an expensive forward model without surrogate bias | posterior error indicator, acceptance correction     |
| [[en/computational-mathematics/paper-notes/scientific-machine-learning/index\|Scientific machine learning]]           | 22     | when networks solve PDEs, sample placement and spectrum decide the outcome | failure probability, frequency-adaptive scales       |
| [[en/computational-mathematics/paper-notes/fbsde-and-control/index\|FBSDEs and stochastic optimal control]]           | 22     | avoid high-dimensional grids while keeping high temporal order             | conditional expectation, martingale property         |
| [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/index\|Phase field and variable steps]]      | 12     | preserve the discrete energy law and maximum bound under adaptive steps    | convolution kernels, step-ratio thresholds           |
| [[en/computational-mathematics/paper-notes/parallel-in-time/index\|Parallel-in-time algorithms]]                      | 13     | turn sequential time stepping into a concurrent algebraic problem          | eigenvector conditioning, all-at-once preconditioner |

The full bibliography is in the [[en/computational-mathematics/paper-notes/catalog|catalogue]], which lists the number, title, venue and corresponding close-reading page for all 108 items by topic.

## Distribution over time

![Topic distribution from 2010 to 2026](assets/diagrams/tao-zhou-papers/en/research-timeline.svg)

The timeline shows three shifts of emphasis. From 2010 to 2015 the work concentrates on stochastic Galerkin methods, stochastic collocation and discrete least squares, all answering one question: how to characterise parametric dependence from a limited number of forward solves. From 2016 to 2021 two lines run in parallel. One pushes sampling design towards optimal sampling and sparse recovery and then into Bayesian inverse problems. The other turns to the time direction, covering energy analysis for variable-step phase-field schemes and diagonalisation-based parallel-in-time algorithms. After 2022 scientific machine learning dominates, but its technical core comes from the earlier two phases: adaptive sampling inherits collocation design, density flows inherit the change of variables, and martingale methods inherit the probabilistic representation of forward-backward stochastic differential equations.

## Relation to the other topics here

[[en/computational-mathematics/knowledge-notes/time-parallelization/index|Time Parallelization for Hyperbolic and Parabolic Problems]] is a section-by-section close reading of item 85, the 2025 _Acta Numerica_ survey by Gander, Wu and Zhou, including the original figures and reproducible experiments. The [[en/computational-mathematics/paper-notes/parallel-in-time/index|parallel-in-time]] page here does not repeat that reading; it explains what each of the twelve original papers underneath the survey actually solved.

## How to read this

- For the overall strategy: start with the three recurring moves above, then open any topic index.
- To follow one technical line: use the "lineage" table on a topic index, which lists the internal dependencies in chronological order.
- To find one specific paper: search the [[en/computational-mathematics/paper-notes/catalog|catalogue]] by its homepage number.

## Three things only visible once all 108 are read together

None of the following is obvious paper by paper. Read across the whole corpus they recur, and each bears directly on whether a stated result can be used as it stands. This site flags them on every page they apply to rather than leaving them to be rediscovered.

### One: "proved" and "demonstrated" are frequently not the same claim

Several papers' own experiments fall **outside** the regime their theorems cover. Paper 9's theorem carries a prefactor of $256$ at $d=2$ and $2304$ at $d=3$, so its examples almost certainly lie outside what it proves; paper 14 tests at $M=2.5N$, below the quadratic requirement of its own theorems; and paper 22 concedes outright that one term of its error bound does not vanish as the sample count grows, so its theory is on its face weaker than standard Monte Carlo theory and it argues from magnitude and experiment instead.

**This is not sloppiness — it is evidence that these sufficient conditions are not sharp.** Anyone sizing a sample count or a step ratio from these bounds should know they are usually conservative. Paper 67 is the one place the gap is quantified: numerically $R_e<1.69$ is necessary while the theory delivers $1.4877$ as sufficient.

### Two: some papers contain no theorem, and some contain no experiment

Both exist here, and both are worth knowing before citing:

- **No theorem.** Paper 54 contains no theorem, lemma, proposition or proof anywhere in its text, so "optimal design" in its title and "quasi-optimal" in its body carry no quantitative content — there is no comparison between the greedy point set and the true maximiser, and no bound on the condition number or the Lebesgue constant. Paper 27 likewise proves no convergence theorem.
- **No experiment.** Papers 58, 74 and 104 are purely analytical, their sources containing only plots of kernels or generating functions; papers 1 and 47 also report none. **The consequence is worth carrying: paper 58's claim of $k$-th order convergence has no measured data behind it.**

### Three: a paper's stated conclusion sometimes disagrees with its own tables

Paper 8 gives two mutually inconsistent rates for $k=6$ ($5.116$ and $6.273$, depending on the mesh range); paper 25's first example reports $2.632$ in the $Y$ component against its third-order claim; paper 93 concedes it trails PINN and RS-PINN yet beats RS-PINN on Allen-Cahn at $d=10^5$ in its own table; and paper 105 claims smoothing uniformly improves on filtering while its Table 9 shows the reverse at $K=50$.
