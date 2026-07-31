---
title: The OneRec Production Pipeline
description: The complete loop from behavior sequences, semantic IDs, and training to online generative retrieval
lang: en
translation: 推荐系统理论/onerec
tags:
  - Recommender-Systems
  - OneRec
  - Internal-Review
---

> [!warning] Temporary internal-review draft
> This section supports internal review and contains module boundaries, parameters, and fault-localization methods. The public version removes environment addresses, keys, real business samples, identities, and uninterpretable experiment IDs. It should still be moved to a private knowledge base after review.

## One-Sentence Definition

OneRec maps ads or content to multilevel discrete semantic IDs, organizes user history into a unified sequence, and uses a generative model to predict candidate SIDs level by level. Online services map those SIDs back to eligible creatives and connect them to existing filtering, coarse-ranking, fine-ranking, and bidding stages.

It does not “replace every recommendation stage with a large model.” It adds a **generative retrieval and side-path ranking route** to an existing production system and ultimately remains subject to inventory, creative, budget, filtering, ranking, and auction constraints.

## Complete Loop

![The complete OneRec loop from semantic IDs to online feedback](assets/diagrams/onerec/en/overview.svg)

## Four Flows That Must Be Understood Together

| Flow           | Origin                                  | Destination                  | Typical risk                                                     |
| -------------- | --------------------------------------- | ---------------------------- | ---------------------------------------------------------------- |
| Data           | Behavior logs, user features, inventory | Training `SequenceExample`   | Temporal misalignment, slot misalignment, missing feature routes |
| Model          | RQ-VAE, sequence model                  | Encoder/Decoder services     | Mismatched versions, vocabulary, or tree depth                   |
| Online control | Mixer request                           | Cache, Ranking Fetch, merge  | Reused parameters across ad positions, timeout, empty result     |
| Feedback       | Impressions, clicks, value signals      | Rewards and the next samples | Attribution delay, selection bias, reward drift                  |

## Page Map

- [[en/recommender-systems/onerec/end-to-end-pipeline\|End-to-End Pipeline]]: service topology and the boundary between old and new retrieval paths.
- [[en/recommender-systems/onerec/data-and-features\|Data and Features]]: USER/CONTEXT/sequence slots and online–offline consistency.
- [[en/recommender-systems/onerec/semantic-id-and-rq-vae\|Semantic IDs and RQ-VAE]]: continuous vectors to hierarchical tokens.
- [[en/recommender-systems/onerec/training-and-sample-generation\|Training and Sample Generation]]: Sample Service, reward features, and GRPO/MCTS-related sampling.
- [[en/recommender-systems/onerec/online-retrieval\|Online Retrieval]]: six core operators and the common downstream path.
- [[en/recommender-systems/onerec/beam-search\|Beam Search]]: beam width, top-k, graph buckets, and cost.
- [[en/recommender-systems/onerec/parameter-reference\|Parameter Reference]]: sources, constraints, and tuning of important parameters.
- [[en/recommender-systems/onerec/stability-and-observability\|Stability and Observability]]: rollout, rollback, funnels, and alerts.

## Conclusions That Matter Most in Review

1. The new path does not bypass fine ranking. It uses an independent Ranking Fetch and side-path ranker before merging with the main path.
2. `beam width` is not the only candidate-volume parameter. Request width, per-level width, final top-k, shard quota, and KV truncation must be considered together.
3. USER/CONTEXT entries in online samples are not globally ordered by time. They are assembled in a fixed slot order, and training must preserve the same semantics.
4. After the old and new retrieval paths enter common filtering and coarse ranking, they can be compared under the same `recall_path` convention. Inside the new path, a segmented SID→TID→AID funnel is required.
5. Rollback requires both request-level and service-level switches, allowing traffic to return to the old path without a release.
