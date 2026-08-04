---
title: GFlowGR Results and Analysis
description: Main results, ablations, diversity observations, and industrial deployment takeaways
lang: en
translation: 推荐系统理论/GFlowGR/结果与分析
tags:
  - GFlowGR
  - Experimental-Results
  - Industrial-Deployment
---

## 1. How to Read the Main Results

On Beauty, Instruments, and Yelp, the paper attaches GFlowGR-DB / TB to TIGER and LETTER and compares against SFT, GRPO, DPO, S-DPO, SPRec, and IPA. The stable public conclusions are:

1. **GFlowGR beats SFT across the board**: augmented trajectories, flow matching, and value alignment all contribute;
2. **TB usually beats DB**: intermediate identifier states do not uniquely identify items, so whole-trajectory balance fits better;
3. **On-policy GRPO is unstable**: it mainly boosts the positive inside a set and struggles to cover high-value regions outside that set;
4. **Among off-policy methods, S-DPO is closest**: multiple negatives matter, but token-level flow supervision is still missing.

Example numbers from LETTER + Beauty:

| Method     |    R@5 |   R@10 |    N@5 |   N@10 |
| ---------- | -----: | -----: | -----: | -----: |
| SFT        | 0.0344 | 0.0540 | 0.0227 | 0.0290 |
| S-DPO      | 0.0390 | 0.0607 | 0.0259 | 0.0336 |
| GFlowGR-DB | 0.0407 | 0.0630 | 0.0266 | 0.0337 |
| GFlowGR-TB | 0.0433 | 0.0672 | 0.0286 | 0.0363 |

When reproducing, prioritize **relative ordering** and **gain direction** over absolute decimals. Framework differences, beam settings, and early stopping can shift the baseline.

## 2. Hyperparameters: Why the Reverse V

Both $N$ and $\lambda$ follow a reverse-V pattern:

- too-small $N$: sparse data lacks diversity and set-wise learning never forms;
- too-large $N$: noisy augmentations start to mislead training;
- too-small $\lambda$: the GFN term barely matters;
- too-large $\lambda$: $\mathcal{L}_{\mathrm{GR}}$ is washed out and basic next-item competence drops.

A practical start is $N=3$, $\lambda=1$, then local refinement. That is stabler than a blind grid search.

## 3. Ablations: Which Modules Are Necessary

Ablations on LETTER + Beauty give four signals:

1. **dropping any reward signal hurts**, and the interaction signal $r_a$ matters most;
2. **random sampling is weaker than CM curriculum**, so augmentation quality drives stability;
3. **$N=1$ (no augmented trajectories) generalizes worse**, confirming that a single positive is not enough;
4. **removing CM entirely still retains strong results**, so the method is not merely “inject collaborative scores into an LLM”.

Reward fusion comparison:

| Strategy         | Takeaway                                 |
| ---------------- | ---------------------------------------- |
| Sum              | strong default and simplest to ship      |
| Weighted Sum     | further top-5 gains; good second default |
| Deep Integration | over-parameterized with no clear upside  |

## 4. LLM Enhancements and Diversity

Two lightweight LLM enhancements help:

- **LLM sampling**: select from a CM-filtered pool by generation confidence; low-confidence samples aid exploration and help R@5 / N@5;
- **LLM reward**: add normalized generation probability as an auxiliary reward with no extra inference.

Diversity analysis on TIGER + Beauty shows that, relative to SFT, GFlowGR’s generation distribution is closer to normal with a smaller standard deviation, while collaborative-score range is wider. In the case study, SFT tends to repeat recent same-scent items; GFlowGR expands to same-brand but different-function / scent combinations with higher utility.

## 5. Industrial Deployment Takeaways

The paper reports GFlowGR usage in Taobao search-ads GR training since May 2025. Publicly citable points include:

- offline gains over SFT / GRPO on H@20 and N@20;
- about 2.1× training cost versus SFT at the same 80× sample budget;
- a 15-day, 10% traffic A/B with roughly +0.43% total revenue;
- larger gains on new items and long-tail queries, plus a small coverage increase.

These numbers show **large-scale validation**. They are not acceptance criteria for public reproduction. Public reproduction still ends on the three academic datasets.

## 6. Engineering Implications

If attaching GFlowGR to an existing generative-retrieval path:

1. **stabilize SFT and tokenizer versions first**, then add GFN;
2. **keep reward scales interpretable**: log interaction level, collaborative score, and token similarity separately;
3. **make the sampling curriculum replayable**: store difficulty quantiles so data drift does not silently break the schedule;
4. **keep inference zero-invasive**: compare SFT and GFlowGR checkpoints under the same beam / constrained decoding;
5. **when comparing with OneRec-style GRPO**, watch whether probabilities are forced to normalize inside the candidate set. GFlowGR emphasizes $P(\tau)\propto R$ and does not require the set probabilities to sum to one.

## 7. Next Reading

- For implementation details, return to [[en/recommender-systems/gflowgr/method\|Method]];
- To run the public protocol, follow the gates in [[en/recommender-systems/gflowgr/reproduction\|Reproduction]];
- For another reward-aware industrial route, compare [[en/recommender-systems/onerec/training-and-sample-generation\|OneRec training and sample generation]].
