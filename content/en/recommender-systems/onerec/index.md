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

It adds a **generative retrieval and side-path ranking route** to the existing production system. Main retrieval continues concurrently in production traffic; the paths merge late in Ranking and then share inventory, creative, budget, filtering, reranking, and auction constraints.

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
- [[en/recommender-systems/onerec/source-level-recall-to-reranking\|Recall, Coarse Ranking, Fine Ranking, and Reranking]]: a source-driven account of concurrent routing, asynchronous cache, dual ranking paths, and ad/creative merge.
- [[en/recommender-systems/onerec/data-and-features\|Data and Features]]: USER/CONTEXT/sequence slots and online–offline consistency.
- [[en/recommender-systems/onerec/semantic-id-and-rq-vae\|Semantic IDs and RQ-VAE]]: continuous vectors to hierarchical tokens.
- [[en/recommender-systems/onerec/training-and-sample-generation\|Training and Sample Generation]]: Sample Service, reward features, and GRPO/MCTS-related sampling.
- [[en/recommender-systems/onerec/online-retrieval\|Online Retrieval]]: six core operators and the common downstream path.
- [[en/recommender-systems/onerec/beam-search\|Beam Search]]: beam width, top-k, graph buckets, and cost.
- [[en/recommender-systems/onerec/parameter-reference\|Parameter Reference]]: sources, constraints, and tuning of important parameters.
- [[en/recommender-systems/onerec/stability-and-observability\|Stability and Observability]]: rollout, rollback, funnels, and alerts.

## Conclusions That Matter Most in Review

1. In production and gray environments, the main Retrieval Proxy and the GprHub side path normally run concurrently. The generative-versus-legacy choice occurs inside GprHub.
2. The new path passes retrieval-side Scoring, Ranking Fetch, side-path DocWash, and independent prediction/pCTR before merging with the main path.
3. `beam width` is one layer of the candidate budget. Per-level width, final top-k, shard quota, KV truncation, coarse top-n, and Ranking merge must be considered together.
4. USER/CONTEXT entries in online samples follow a fixed slot order, which training must preserve.
5. New and legacy GPR implementations can be compared by `recall_path`; the main-versus-side analysis also needs main-only, side-only, and overlap attribution.
6. Rollback needs a side-path request gate, a generative-implementation gate, and a service-capability gate, while every side-path failure remains fail-open for the main path.
