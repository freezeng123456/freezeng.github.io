---
title: "MCoTTA: Momentum-Guided Continual Test-Time Adaptation"
description: Three-seed paired results, causal gradient-retention diagnostics, and long-horizon limits
lang: en
translation: 计算机科学/个人科研/mcotta
tags:
  - Computer-Science
  - Test-Time-Adaptation
  - AAAI
---

Project: [freezeng123456/MCoTTA](https://github.com/freezeng123456/MCoTTA)

> [!info] Result snapshot
> This page reflects the repository state on July 31, 2026. The previously omitted ViT-B/16 seed-3 and CIFAR-100-C 10-cycle experiment branches have been merged into the main branch. The public account retains the methodological insight, algorithmic structure, and verifiable results without exposing all implementation details of the paper.

## Research Question

In continual test-time adaptation, LCoTTA collects entropy-gradient snapshots, obtains a low-rank subspace through sliding-window PCA, and projects the current gradient into that subspace before updating the model. The difficulty is that a single-batch raw gradient contains substantial transient noise. The PCA queue may therefore keep following short-lived directions instead of tracking slowly varying structure in the corruption stream.

The central observation is that PCA selects directions according to historical variance, whereas continual adaptation needs directions that remain informative for subsequent batches. High variance does not imply temporal persistence. A brief but large gradient change can consume the limited subspace capacity while smoother directions are excluded.

## Algorithm

The original method writes the current filtered gradient $g_t$ directly to the PCA queue:

$$
\mathcal Q_t\leftarrow g_t.
$$

MCoTTA maintains an exponential momentum:

$$
m_t=\beta m_{t-1}+(1-\beta)g_t,
\qquad
\mathcal Q_t\leftarrow m_t.
$$

Centered rank-$r$ PCA on the momentum history gives a projection matrix $P_t^M$. The current gradient is still used for the actual update:

$$
\tilde g_t=P_t^M g_t,
\qquad
\theta_{t+1}=\theta_t-\eta\tilde g_t.
$$

The algorithm deliberately separates two paths:

- **Estimation path**: $g_t\rightarrow m_t\rightarrow\mathcal Q_t\rightarrow P_t^M$;
- **Update path**: $g_t\rightarrow P_t^M g_t\rightarrow\theta_{t+1}$.

The momentum state determines which directions are available for adaptation, while the current batch still determines the step along those directions. The method adds no teacher–student network, memory bank, auxiliary loss, or trainable parameters.

## Theoretical Insight

Write the gradient as a temporally persistent low-rank component and a faster-varying residual:

$$
g_t=s_t+\epsilon_t.
$$

Exponential averaging reduces the relative covariance contribution of uncorrelated residuals while retaining components that are correlated across batches. The core claim is therefore not that “the momentum gradient is more accurate,” but that:

> A momentum history is better suited to estimating an update subspace that will remain relevant, while the current gradient should still drive the parameter update.

The theory establishes a persistent-energy ordering at the population-covariance level. It explains why a smoothed history can be more suitable than a raw-gradient history for subspace estimation, but it does not equate “retaining more future energy” with improved classification accuracy in every direction.

## Main Experimental Protocol

All three tasks use corruption severity 5, all 15 corruption types, a single continuous pass without resets at corruption boundaries, batch size 64, and paired LCoTTA/MCoTTA runs with the same random stream for each seed.

| Task                     | Samples per corruption | $W$ | $r$ | $K$ | Main $\beta$ |
| ------------------------ | ---------------------: | --: | --: | --: | -----------: |
| CIFAR-100-C / ResNeXt-29 |                 10,000 | 100 |  25 |  50 |         0.99 |
| ImageNet-C / ResNet-50   |                 50,000 | 100 |  25 |  50 |         0.99 |
| ImageNet-C / ViT-B/16    |                 50,000 | 100 |  50 | 100 |         0.99 |

Here $W$ is the PCA queue length, $r$ is the subspace rank, and $K$ is the interval between gradient snapshots.

## Three-Seed Paired Results

After the omitted branches were merged, all three tasks have complete three-seed paired results. The table reports mean accuracy ± sample standard deviation:

| Task                     |        LCoTTA |            MCoTTA | Mean paired gain |
| ------------------------ | ------------: | ----------------: | ---------------: |
| CIFAR-100-C / ResNeXt-29 | 65.37 ± 0.21% | **66.00 ± 0.09%** |        **+0.63** |
| ImageNet-C / ResNet-50   | 36.79 ± 0.17% | **43.97 ± 0.38%** |        **+7.18** |
| ImageNet-C / ViT-B/16    | 60.18 ± 0.40% | **62.98 ± 0.54%** |        **+2.80** |

The seed-wise paired gains are:

| Task                     | seed 1 | seed 2 | seed 3 |
| ------------------------ | -----: | -----: | -----: |
| CIFAR-100-C / ResNeXt-29 |  +0.54 |  +0.94 |  +0.42 |
| ImageNet-C / ResNet-50   |  +7.75 |  +6.90 |  +6.89 |
| ImageNet-C / ViT-B/16    |  +2.94 |  +2.80 |  +2.66 |

The strongest evidence is concentrated in the two ImageNet-C settings. The CIFAR-100-C improvement has the same direction for all three seeds but is markedly smaller; it is therefore not presented as a benefit of comparable magnitude across tasks.

## Causal Mechanism Diagnostic

Classification outcomes reflect both subspace quality and diverging update trajectories. To isolate the former, shared-gradient replay fixes the ImageNet-C / ResNet-50 gradient trajectory, freezes the model, and feeds exactly the same gradient sequence to each tracker. Evaluation uses only the next 1–49 batches and excludes the current snapshot $k=0$.

With $W=100$, $r=25$, $K=50$, and 6,595 scored records:

| Subspace history        | Mean future-gradient energy retained |
| ----------------------- | -----------------------------------: |
| Raw gradient            |                              1.3275% |
| Momentum, $\beta=0.5$   |                              1.7891% |
| Momentum, $\beta=0.7$   |                              2.5831% |
| Momentum, $\beta=0.9$   |                              3.3133% |
| Momentum, $\beta=0.99$  |                          **3.6882%** |
| Momentum, $\beta=0.999$ |                              3.6491% |

The retention rate at $\beta=0.99$ is approximately 2.78 times that of the raw-gradient history. It no longer increases at $\beta=0.999$, indicating that the benefit begins to saturate under excessive smoothing.

![[assets/research/mcotta-retention.png]]

## Parameters and Ablations

### Momentum coefficient

The complete single-seed $\beta$ sweep shows that stronger momentum matters more in the two ImageNet-C settings:

| $\beta$ | CIFAR-100-C / ResNeXt-29 | ImageNet-C / ResNet-50 | ImageNet-C / ViT-B/16 |
| ------: | -----------------------: | ---------------------: | --------------------: |
|    0.10 |                   65.69% |                 36.94% |                44.18% |
|    0.30 |                   64.95% |                 38.02% |                60.40% |
|    0.50 |                   65.00% |                 38.52% |                60.80% |
|    0.70 |                   65.66% |                 40.33% |                60.31% |
|    0.90 |                   66.06% |                 43.03% |                62.63% |
|    0.99 |                   66.04% |             **44.39%** |                63.56% |
|   0.999 |               **66.16%** |                 43.77% |            **64.07%** |

$\beta=0.99$ is not the isolated optimum for every task, but it is stable across all three and is therefore used in the main three-seed experiments.

### Subspace rank

The rank ablation shows that a larger subspace is not automatically better:

| Task                     |     $r=10$ | $r=50$ |     $r=75$ |
| ------------------------ | ---------: | -----: | ---------: |
| CIFAR-100-C / ResNeXt-29 | **66.52%** | 65.99% |     66.17% |
| ImageNet-C / ResNet-50   | **43.01%** | 42.26% |     40.31% |
| ImageNet-C / ViT-B/16    |     62.92% | 63.46% | **63.51%** |

ResNet-50 degrades clearly as rank increases, suggesting that additional directions can reintroduce transient variation. ViT is more tolerant of higher rank, so the main configuration uses $r=50$.

### Queue length

Under the common ablation setting $r=25$, $K=50$:

| Task                     |     $W=35$ |     $W=50$ |     $W=75$ |
| ------------------------ | ---------: | ---------: | ---------: |
| CIFAR-100-C / ResNeXt-29 |     66.04% | **66.21%** |     65.73% |
| ImageNet-C / ResNet-50   |     42.26% |     43.07% | **43.50%** |
| ImageNet-C / ViT-B/16    | **64.12%** |     64.09% |     63.39% |

The ViT window ablation uses $r=25,K=50$, unlike the main experiment with $r=50,K=100$. Its entries therefore cannot be compared with the main $W=100$ result as a single-variable experiment.

## Long-Horizon Behavior

The ImageNet-C / ResNet-50 experiment at $\beta=0.9$ completed 10 cycles. LCoTTA declined from 36.83% to 35.79%, while MCoTTA declined from 42.43% to 39.76%. MCoTTA remained higher in every cycle, but the gap narrowed over time, and neither method prevented long-term drift.

The paired CIFAR-100-C / ResNeXt-29 experiment with seed 2 and $\beta=0.9$ also completed 10 cycles:

| Method |    cycle 1 |       peak |   cycle 10 |
| ------ | ---------: | ---------: | ---------: |
| LCoTTA |     65.47% |     65.78% |     65.68% |
| MCoTTA | **66.16%** | **66.83%** | **66.10%** |

The classification gain is smaller on this task, but the 10-cycle result does not show the pronounced continuing degradation observed on ImageNet-C / ResNet-50.

## Scientific Boundaries

- The three-seed results are complete, but the sample size remains insufficient to replace larger-scale statistical testing.
- The mechanistic conclusion is that the momentum subspace retains more future-gradient energy, not that every retained direction benefits classification.
- Improvements are substantial in the two ImageNet-C settings and smaller on CIFAR-100-C; an equal-magnitude benefit across all datasets is not claimed.
- The long-horizon experiments show that MCoTTA improves update directions but does not solve drift under indefinitely continuing adaptation.
- Raw-PCA has an in-sample self-inclusion effect when the current gradient is inserted before that same gradient is projected. A strict causal comparison must project first and insert the snapshot afterward; conclusions from that comparison must be reported separately from the existing experiments.
- The page does not disclose the paper's complete proofs, implementation techniques, or all per-corruption results.
