---
title: GFlowGR Reproduction Guide
description: Public datasets, baselines, hyperparameters, training steps, and acceptance gates
lang: en
translation: 推荐系统理论/GFlowGR/实验复现
tags:
  - GFlowGR
  - Paper-Reproduction
  - Experiment-Protocol
---

## Reproduction Goal

On public data, verify two claims:

1. attaching GFlowGR to TIGER / LETTER improves over SFT and common RLFT baselines;
2. the key ablations hold: CM curriculum sampling, multi-signal rewards, and reverse-V curves in $N$ / $\lambda$.

The official repository currently ships only a README; code is still being prepared. Reproduce from the paper on top of existing GR scaffolds rather than waiting for private scripts.

![A six-step path from data to evaluation](assets/diagrams/gflowgr/en/reproduction-pipeline.svg)

## 1. Data and Splits

| Dataset     |  Users |  Items | Sparsity | Avg.length |
| ----------- | -----: | -----: | -------: | ---------: |
| Yelp        | 30,431 | 20,032 |   99.85% |       9.39 |
| Beauty      | 22,362 | 12,101 |   99.92% |       7.87 |
| Instruments | 24,772 |  9,922 |   99.99% |       7.32 |

Preprocess as in LETTER:

- Beauty / Instruments: tokenize from title + description;
- Yelp: tokenize from title + description + brand + category;
- split: last interaction for test, second-to-last for validation, the rest for training.

Implementation checks:

- user sequences must be time-ordered, with no test-item leakage into prompts;
- tokenizer training may use only training-period item text features;
- train / val / test targets for the same user must not overlap.

## 2. Backbones and Baselines

Install GFlowGR on two GR backbones:

- **TIGER**: classic RQ-VAE identifiers plus generative retrieval;
- **LETTER**: tokenization with collaborative embeddings and diversity constraints.

Compare three fine-tuning families:

| Family          | Methods                   | Reproduction notes                                                                 |
| --------------- | ------------------------- | ---------------------------------------------------------------------------------- |
| SFT             | next-token prediction     | Must land first; every gain is measured against it                                 |
| On-policy RLFT  | GRPO                      | Relative advantages inside a group; weak coverage outside the set                  |
| Off-policy RLFT | DPO / S-DPO / SPRec / IPA | S-DPO consumes multiple negatives; IPA should reuse GFlowGR’s collaborative reward |

Metrics: R@5, R@10, N@5, N@10. Repeat each run with seeds `{42,43,44}` and report the mean.

## 3. Default Hyperparameters

Default ensemble from the paper:

- trajectory sampling: CM-based curriculum;
- reward: $R=r_a+r_{\hat{y}}+r_{\mathrm{sim}}$;
- loss: $\mathcal{L}_{\mathrm{GR}}+\lambda\sum\mathcal{L}_{\mathrm{GFN}}$;
- stopping: Transformer Trainer early stopping.

Sensitivity grids:

| Hyperparameter             | Grid                    | Empirical sweet spot |
| -------------------------- | ----------------------- | -------------------- |
| Augmented trajectories $N$ | $\{0,1,3,5\}$           | around $3$           |
| Loss weight $\lambda$      | $\{0.01,0.1,1,10,100\}$ | around $1$           |

$N=0$ almost collapses to a single positive; $N=5$ injects noise; $\lambda\ge 10$ overpowers SFT and hurts basic recommendation.

## 4. Minimal Implementation Checklist

Land modules one by one instead of rewriting a full GR stack:

1. **Tokenizer / Backbone**
   - Reuse TIGER or LETTER training scripts;
   - fix identifier length $L=3$ and export the item → token map.
2. **Collaborative Model**
   - train a scorer $\mathrm{CM}(U,v)$;
   - cache scores for both sampling and rewards.
3. **Trajectory Sampler**
   - input a positive, output $N$ token trajectories;
   - implement at least random and CM curriculum modes for ablation.
4. **Reward Model**
   - emit $r_a$, $r_{\hat{y}}$, and $r_{\mathrm{sim}}$;
   - start with Sum, then optionally Weighted Sum.
5. **GFN Loss**
   - read $P_{\mathrm{GR}}$ along each trajectory from LLM logits;
   - implement both DB and TB; use TB as the primary result.
6. **Trainer**
   - each batch: positive SFT + set-wise GFN;
   - log $\mathcal{L}_{\mathrm{GR}}$, $\mathcal{L}_{\mathrm{GFN}}$, mean reward, and sampling difficulty.

Pseudocode:

```python
for batch in loader:
    prompt, positive, aug_items = sample_trajectories(batch, N, strategy="cm_curriculum")
    trajectories = [positive, *aug_items]
    sft_loss = next_token_loss(model, prompt, positive)
    gfn_loss = 0.0
    for traj in trajectories:
        logp = token_logprobs(model, prompt, traj)
        reward = r_a(traj) + r_cm(traj) + r_sim(traj, positive)
        gfn_loss += trajectory_balance(logp, reward, Z)  # or detailed_balance(...)
    loss = sft_loss + lambda_gfn * gfn_loss
    loss.backward(); optimizer.step()
```

## 5. Suggested Experiment Order

1. **SFT baseline** on Beauty + TIGER;
2. **GFlowGR-TB default** with CM sampling, Sum reward, $\lambda=1$, $N=3$;
3. **GFlowGR-DB** under the same setting;
4. **hyperparameter curves** over $N$ and $\lambda$;
5. **ablations**: w/o Ada, w/o Traj, w/o CM, single-signal rewards;
6. **extension** to LETTER and to Instruments / Yelp;
7. **optional enhancements**: LLM confidence sampling and LLM probability reward.

## 6. Acceptance Gates

Treat the public mainline as reproduced when:

1. on the same backbone, GFlowGR-TB beats SFT on all four metrics;
2. TB is at least as strong as DB, especially on Yelp / TIGER where intermediate-state ambiguity is larger;
3. $N$ and $\lambda$ show the reverse-V shape, with optima in the middle of the grid;
4. removing augmented trajectories, switching to random sampling, or dropping any reward signal hurts;
5. the inference path matches the SFT checkpoint exactly; only the fine-tuning objective changes.

## 7. Common Failure Modes

| Symptom                 | Likely cause                                 | Fix                                               |
| ----------------------- | -------------------------------------------- | ------------------------------------------------- |
| No gain over SFT        | $\lambda$ too large or unnormalized rewards  | lock $\lambda=1$ and standardize $r_{\hat{y}}$    |
| Late collapse           | negatives too easy or too hard               | enable CM curriculum and bound the candidate pool |
| TB worse than SFT       | $Z$ not learned, or positives ignored        | keep positive SFT and verify $R>0$                |
| High variance           | inconsistent seeds / early stop / beam width | fix three seeds and one generation config         |
| Gains look like CM only | need the w/o CM variant                      | confirm robustness without collaborative signals  |

## 8. Non-Reproducible Pieces

These can be cited from the paper but not rebuilt in public settings:

- full Taobao search-ads training data and features;
- production collaborative models and auction coupling;
- online A/B traffic splits and business metric definitions.

The public reproduction endpoint is: **recover the relative conclusions of Table 2 and the key ablation trends on three public datasets and two backbones**.
