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

| Parameter                        | Effective source/fallback                        | Hard constraint                                   | Primary effect                        |
| -------------------------------- | ------------------------------------------------ | ------------------------------------------------- | ------------------------------------- |
| `gpr_pam_beam_width`             | Position experiment; Entry falls back to 180     | Current Decoder maximum 180                       | Decoder search width, SID count, P99  |
| Shared-structure initializer     | Source contains 200, but Entry overwrites it     | Must not be treated as the online effective value | Reveals default-override risk only    |
| `gpr_beam_search_topk`           | Position experiment; effective when positive     | No greater than output SIDs                       | Final SID truncation                  |
| `gpr_pam_beam_width_shard_ratio` | Position experiment; safe missing-value behavior | Shards must cover the request budget in aggregate | Shard load and candidate distribution |
| Trie level/version               | Trie/model metadata                              | Must agree with model, SID packing, and KV        | Token levels and KV hit rate          |
| Uniqueness filter                | Quota/policy plan                                | Consistent SID format                             | Duplicate-SID removal                 |
| KV/output quota                  | Policy value; fallback 1000                      | Nonnegative and bounded                           | TID truncation after SID expansion    |

See [[en/recommender-systems/onerec/beam-search\|Beam Search]] for detailed beam semantics.

## Online RPCs and Timeouts

| Parameter                  | Source and precedence                                             | Explanation                                                                                |
| -------------------------- | ----------------------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| Datahub timeout            | Position/service configuration                                    | Retrieves sequence and generative context; RPC failure and business emptiness are separate |
| Encoder/Decoder/KV timeout | Per-service stage configuration                                   | Each stage needs its own budget rather than one opaque infer timeout                       |
| Ranking Fetch timeout      | Position-level `gpr_server_hub_timeout_ms`, then service fallback | Waits for Retrieval Cache `fill_done` and fails open at expiry                             |
| Cache normal/hard expiry   | Cache-service configuration                                       | Cleans filled-but-unconsumed nodes and producers that never complete                       |

Fixed values in older review notes represent integration snapshots and must not be presented as current production settings. Fetch timeout has to cover an acceptable cache-ready percentile while remaining inside the total Ranking deadline. Increasing it may recover side-path coverage but can directly become final P99.

## Ranking, Merge, and Reranking

| Parameter                     | Effective rule                                             | Primary effect                                       |
| ----------------------------- | ---------------------------------------------------------- | ---------------------------------------------------- |
| `ranking_enable_fetch_gpr`    | Service capability and position experiment must both match | Whether the side path enters Ranking                 |
| `gpr_bind_adtable_batch_size` | Code fallback 10                                           | Ad-table binding batches, throughput, tail latency   |
| GPR prediction / pCTR gates   | Controlled by position experiment and Ranking capability   | Side-path fine-ranking score and cost                |
| OneRec relevance gate         | Global gate + position experiment + designated scene       | Whether generative relevance enters common reranking |
| Relevance threshold           | Score is clamped to $[0,1]$ after thresholding             | Factor coverage and reranking strength               |

These parameters need main-only, side-only, AID-overlap, and creative-overlap buckets. Final ad count alone hides Fetch gating, missing binding fields, DocWash removal, and merge deduplication.

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

For beam width: Entry parses it, the fallback is 180, the current Decoder maximum is 180, and Decoder latency/SID count change first, followed by KV, coarse ranking, Ranking Fetch coverage, and final ad count.

See [[en/recommender-systems/onerec/source-level-recall-to-reranking|Recall, Coarse Ranking, Fine Ranking, and Reranking]] for the complete effect chain.
