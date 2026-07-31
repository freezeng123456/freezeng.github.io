---
title: OneRec Training and Sample Generation
description: Sample Service, candidate selection, reward features, and sequence-model training
lang: en
translation: 推荐系统理论/onerec/训练与样本生成
tags:
  - OneRec
  - Training
  - Reinforcement-Learning
  - Internal-Review
---

## Main Sample Service Flow

The sample service can read impression and feedback data from a message queue or local input. Its core generator executes:

```text
RebuildMixerRequest
  → Mixer simulation
  → ParseFeatures
  → PrepareSequence
  → PickItems
  → GenerateFeatures/Rewards
  → SequenceExample or MID_TD
  → downstream storage
```

Mixer simulation is not intended to reproduce online latency exactly. It reconstructs candidates, positions, and experimental context so that the training samples match the online decision space.

## Candidate Selection

The current OneRec generator first takes the top $M$ items and requires at least $M$ candidates. Important parameters are:

| Parameter           | Example value | Meaning                                                       |
| ------------------- | ------------: | ------------------------------------------------------------- |
| `return_ad_count`   |            40 | Maximum number of ads considered or returned by the generator |
| `grpo_m_items_kept` |            20 | Number $M$ of items retained for GRPO                         |
| `creative_top_k`    |          1024 | Top-k entering the creative candidate pool                    |
| `creative_sample_n` |            20 | Number sampled from the candidate pool                        |
| `sampling_percent`  |           100 | Sample-traffic percentage                                     |
| `gen_type`          |      `MID_TD` | Downstream training-data format                               |

These quantities form a hierarchy: `creative_top_k` defines the upstream pool, `creative_sample_n` the sample size, `return_ad_count` the ad-output boundary, and $M$ the count retained for optimization or within-group comparison. Any upstream value smaller than downstream demand causes implicit truncation.

## Reward and Training Features

The generator writes several classes of signals:

- value: basic, quality, and raised eCPM, plus economic value;
- ranking: rank score and recommendation score;
- prediction: pCTR- and CVR-related terms;
- experience: relevance, negative feedback, repetition, and watch time;
- semantics: LLM relevance and long-term-interest relevance;
- business constraints and audience filtering.

The principle is to separate the model's generation probability from business value and constraints. The generative model learns the sequential or semantic distribution, while policy optimization uses multidimensional rewards. Reward scales must be normalized or explicitly weighted; otherwise, fields with large numeric magnitudes dominate within-group advantages.

## External Reward RPCs

LLM evaluation and other reward services can run concurrently to reduce sample-generation latency. Multiple callbacks, however, must not write to the same protobuf feature map at the same time. The current implementation separates concurrent RPC calls from serialized final feature writes to avoid map data races.

In the example configuration, LLM evaluation is enabled for only about 5% of samples. `llm_sampling_rate=0.05` controls expensive reward-label coverage, not the overall sampling rate of training data.

## GRPO/MCTS-Style Parameters

The configuration supports:

- top-k sampling;
- probability thresholds and min-p;
- MCTS candidate and branching controls;
- secondary sampling;
- forced retention of impressed items;
- retention of $M$ items within a group.

These switches must be versioned with the data format. If $M$, candidate counts, or reward fields change while the trainer still reads an old shape, the best outcome is an explicit error; the worst is silent misalignment.

## Sequence-Generation Model

Training data feed USER, CONTEXT, and behavior/ad sequences into the sequence model, which predicts SID tokens level by level:

$$
p(c_1,\ldots,c_L\mid x)
=\prod_{l=1}^{L}p(c_l\mid x,c_{<l}).
$$

Online, the Encoder computes a shared user representation and the Decoder expands tokens level by level. Training must therefore cover the online conventions for:

- token levels and masking;
- prefix conditions;
- padding and truncation;
- special BOS/EOS or business tokens;
- model `data_name`;
- sequence-schema version.

## Data-Quality Checks

| Check                                         | Meaning of failure                                |
| --------------------------------------------- | ------------------------------------------------- |
| Missing rate and length distribution per slot | Feature-routing or USD-completion drift           |
| SID frequency and code usage at each level    | RQ-VAE dead codes or distribution shift           |
| Full-SID collision                            | Insufficient quantization granularity             |
| Reward quantiles and correlations             | Reward-service or unit drift                      |
| Candidate-count funnel                        | Inconsistent `top_k`, sampling, filtering, or $M$ |
| Online/offline reader checksum                | Schema or model-version mismatch                  |

## Training Release Gate

After offline model metrics pass, the release should still include:

1. Encoder/Decoder golden tests on fixed samples;
2. reproducibility and maximum-width stress tests for beam search;
3. replay of SID→TID coverage;
4. shadow comparison of the old and new models on the same requests;
5. low-traffic online rollout that monitors retrieval volume, empty results, latency, and business guardrails.
