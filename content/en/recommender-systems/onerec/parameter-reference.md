---
title: OneRec Parameter Reference
description: Parameter sources, current values, constraints, effects, and troubleshooting
lang: en
translation: 推荐系统理论/onerec/参数手册
tags:
  - OneRec
  - Parameters
  - Internal-Review
---

Parameters must be grouped by call path. The same name appearing in sampling, training, and online serving does not imply the same semantics.

## Online Generative Retrieval

| Parameter          |         Current value/default | Hard constraint                           | Primary effect                             |
| ------------------ | ----------------------------: | ----------------------------------------- | ------------------------------------------ |
| `pam_beam_width`   |         New-path fallback 180 | $\le180$                                  | Decoder prefix-search width                |
| legacy beam width  |                           200 | Constrained by legacy model/engine        | Legacy path only                           |
| `beam_search_topk` |             Applied when $>0$ | $\le$ output candidates                   | Final SID truncation                       |
| beam shard ratio   |        Strategy configuration | Shards must cover request budget in total | Shard load and retrieval distribution      |
| tree level         |               From trie/model | Must agree with model and SID protocol    | Number of token-generation levels          |
| multi-level        |                          true | Requires a multilevel model               | Whether generation proceeds level by level |
| unique filter      |            Strategy-dependent | Consistent SID format                     | Removes duplicate SIDs                     |
| KV/output quota    | Strategy value; fallback 1000 | Nonnegative and bounded                   | TID truncation after SID expansion         |

See [[en/recommender-systems/onerec/beam-search\|Beam Search]] for detailed beam semantics.

## Online RPCs and Timeouts

| Parameter             |         Current default/review value | Explanation                                   |
| --------------------- | -----------------------------------: | --------------------------------------------- |
| Datahub timeout       |                              2000 ms | Retrieves the sequence and generative context |
| infer timeout         |                              3000 ms | Overall or segmented Encoder/Decoder/KV limit |
| Ranking Fetch timeout | Must be no smaller than pipeline P99 | Waits for Retrieval Cache `fill_done`         |

One review snapshot had a Ranking Fetch timeout of only 140 ms while the integration pipeline took approximately 400–500 ms. These figures represent a configuration snapshot that still required correction; the value must be recalibrated against current P99 before release. A timeout cannot simply be increased because it shares the total main-path latency budget. The priority is to reduce pipeline latency and allow fast degradation of a failing side path.

## Sampling and Policy Optimization

| Parameter           |  Example | Role                                   | Tuning risk                                               |
| ------------------- | -------: | -------------------------------------- | --------------------------------------------------------- |
| `return_ad_count`   |       40 | Maximum ads returned by generator      | Too small truncates; too large raises sample/ranking cost |
| `grpo_m_items_kept` |       20 | Candidates retained within a group     | Must not exceed valid candidates                          |
| `creative_top_k`    |     1024 | Upstream candidate pool                | Excessive values admit long-tail noise                    |
| `creative_sample_n` |       20 | Number of creatives sampled            | Affects diversity and variance                            |
| `sampling_percent`  |      100 | Overall sample rate                    | Changes volume and distribution                           |
| `llm_sampling_rate` |     0.05 | Coverage of expensive semantic rewards | Cost–variance tradeoff                                    |
| `gen_type`          | `MID_TD` | Output format                          | Trainer must use the same schema                          |

## RQ-VAE

| Parameter      |                            Current baseline | Effect                                            |
| -------------- | ------------------------------------------: | ------------------------------------------------- |
| input dim      |                                        1024 | Input-representation contract                     |
| latent dim     |                                         128 | Quantization-space capacity                       |
| codebook sizes |                                 256/256/256 | Capacity of the three SID levels                  |
| batch size     |                                        1024 | Stability of code-use statistics                  |
| learning rate  |                                   $10^{-3}$ | Convergence                                       |
| weight decay   |                                   $10^{-4}$ | Regularization                                    |
| loss           |                          MSE + quantization | Balance between reconstruction and discretization |
| K-means init   | Disabled by default, at most 100 iterations | Codebook initialization                           |
| EMA codebook   |                         Disabled by default | Code-update method                                |
| Sinkhorn eps   |                                       0/0/0 | Balanced-assignment switch/strength               |

The key is not memorizing the values but synchronizing codebook/version metadata with online SID parsing.

## Sequence Model and Online Schema

Parameters that require versioning but should not be tuned in isolation include:

- USER/CONTEXT slot order;
- maximum lengths of primary and auxiliary sequences;
- time cutoff;
- token levels and bit width;
- model `data_name`;
- padding, BOS/EOS, and mask;
- Encoder output dimension;
- graph buckets supported by the Decoder.

These are interface contracts. Changing them changes the protocol and requires an end-to-end release.

## Parameter-Source Precedence

A unified rule is recommended:

```text
request/experiment override
  > strategy/quota plan
  > model metadata
  > service config
  > code fallback
```

Each online request should report the final resolved values rather than only the raw experiment parameters. Otherwise, engineers cannot determine which value actually took effect after multiple defaults and overrides.

## Four Questions for Every Parameter

1. Who writes it, and which module parses it?
2. What is the fallback when it is absent?
3. Is the hard boundary imposed by code, the model, or hardware?
4. Which metric changes first when the parameter changes?

For beam width: Entry parses it, the fallback is 180, the model maximum is 180, and Decoder latency/SID count change first, followed by KV, coarse ranking, and the final ad count.
