---
title: Reproducing GFlowGR
description: Using GFlowNets to move generative recommendation from point-wise SFT to value-aware set-wise fine-tuning
lang: en
translation: 推荐系统理论/GFlowGR
tags:
  - Recommender-Systems
  - Generative-Recommendation
  - GFlowNet
  - Paper-Reproduction
---

> [!note] Source and public scope
> This series is based on the public paper [GFlowGR: Fine-tuning Generative Recommendation Frameworks with Generative Flow Networks](https://arxiv.org/abs/2506.16114) (SIGIR ’26), with reproduction boundaries checked against the official repository [Applied-Machine-Learning-Lab/SIGIR26_GFlowGR](https://github.com/Applied-Machine-Learning-Lab/SIGIR26_GFlowGR). It keeps only publicly shareable method structure, formulas, experimental protocol, and result interpretation. It excludes internal business details, credentials, and unauthorized data.

## Executive Summary

GFlowGR is not another tokenizer redesign. It targets a mismatch between **GR fine-tuning** and **online serving**:

1. SFT fits one ground-truth item, while production needs a set of high-value candidates;
2. Interactions such as impression, click, and purchase have different utilities, yet SFT treats them as equal positives;
3. Reward-based methods such as DPO / GRPO can use sets or preferences, but usually provide only item-level rewards and miss token-level supervision.

GFlowGR treats autoregressive item-identifier generation as a GFlowNet trajectory and enforces

$$
P(\tau)\propto R(s_L).
$$

Higher-value candidates therefore receive higher generation probability, while every token transition gets a flow-consistent gradient.

![The three modules of GFlowGR](assets/diagrams/gflowgr/en/framework.svg)

## Page Map

1. [[en/recommender-systems/gflowgr/method\|Method]] covers problem formulation, trajectory definition, and how the sampler, reward model, and loss fit together.
2. [[en/recommender-systems/gflowgr/reproduction\|Reproduction]] covers datasets, baselines, hyperparameters, training steps, and acceptance gates.
3. [[en/recommender-systems/gflowgr/results\|Results]] covers main tables, ablations, deployment takeaways, and differences from SFT / RLFT.

## Reproduction Boundary

| Item | Status | Implication |
| ---- | ------ | ----------- |
| Method and formulas | Public | The training objective and algorithm can be rebuilt |
| Public experiment setup | Public | Beauty / Instruments / Yelp with TIGER / LETTER |
| Official code | README says still being prepared | Implement from the paper on top of existing GR scaffolds |
| Taobao production data and configs | Not public | Cite relative online metrics only; do not claim production replay |

## Relation to Nearby Notes

- Unlike [[en/recommender-systems/long-sequence-modeling\|long-sequence modeling]], GFlowGR focuses on the **fine-tuning objective of generative recommendation**, not long-history encoding.
- It sits next to [[en/recommender-systems/onerec/index\|OneRec]]: both discuss set-wise / reward-aware training. OneRec emphasizes an industrial closed loop; GFlowGR emphasizes how a GFlowNet objective attaches to a GR backbone.
