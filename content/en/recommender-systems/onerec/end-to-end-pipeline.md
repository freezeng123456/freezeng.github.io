---
title: OneRec End-to-End Pipeline
description: Training loop, main/side online topology, and merge boundaries
lang: en
translation: 推荐系统理论/onerec/端到端链路
tags:
  - OneRec
  - Architecture
  - Internal-Review
---

> [!warning] Temporary internal-review draft
> The names below explain module boundaries. Environments, addresses, experiment IDs, and real business fields have been redacted.

## 1. Offline Training Loop

```text
Behavior/impression logs
  → LogToFeature / USD sequence completion
  → OneRec SequenceExample
  → Sample filtering and reward features
  → RQ-VAE semantic IDs
  → HSTU/sequence-generation model
  → checkpoint + model metadata
  → Encoder / Decoder release
```

RQ-VAE and the sequence model are two version-coupled components. Token encodings in sequence-training data must correspond to the same codebook, number of levels, and SID-packing rules. A model release must synchronize not only the checkpoint but also:

- `model_version` and `data_name`;
- token/codebook version;
- trie or KV-creative version;
- bit width per token level and tree depth;
- maximum beam supported by the Encoder/Decoder.

## 2. Online Routing

After Mixer builds the retrieval request, the main Retrieval Proxy and the GprHub side path become sibling tasks. Production and gray environments normally keep both active; a GprHub experiment hit does not stop main retrieval.

Inside GprHub, a position-level generative switch selects:

- generative implementation: `GenerativeRecall` enters the six-operator DAG;
- control implementation: the side path calls a Retrieval Proxy flow compatible with legacy GPR.

Only an isolated simulation mode skips main retrieval under dedicated conditions. GprHub splits multi-position requests at the service layer and gives every position its own `pos_id`, experiment parameters, cache node, and monitoring tag.

## 3. New Generative Retrieval DAG

![The OneRec production path from entry to retrieval cache](assets/diagrams/onerec/en/retrieval-pipeline.svg)

The six generative nodes have the following responsibilities:

| Node                | Input                                          | Output                                            | Failure symptom                                   |
| ------------------- | ---------------------------------------------- | ------------------------------------------------- | ------------------------------------------------- |
| Entry               | Single-position request, experiment parameters | `GenerativeData`, quota, beam parameters          | Missing parameters or route not selected          |
| Datahub             | User and context identifiers                   | Generative behavior sequence, model-input context | Empty sequence or RPC timeout                     |
| Encoder             | Structured sequence                            | Shared user embedding, model version              | Input schema or model-version error               |
| Decoder             | Embedding, tree/model parameters               | SID sequences and scores                          | Empty SID or beam out of range                    |
| KV Creative         | SID                                            | TID list                                          | Version mismatch, SID miss, empty shard response  |
| Creative List Build | TIDs and scores                                | Creative list deduplicated by AID                 | Missing TID→AID mapping or complete deduplication |

The result then enters GprHub's common quota, property-filtering, creative-service, and coarse-ranking graph. The generative and legacy-control graphs are broadly isomorphic, although the current static graphs retain a small number of node differences that must be accounted for in attribution.

## 4. Ranking and Main-Path Merge

The main path continues through conventional retrieval, coarse ranking, and primary Ranking. In GprHub, OneRec completes common filters and retrieval-side Scoring before writing Retrieval Cache. Ranking Fetch waits for `fill_done`, then performs candidate binding, lightweight DocWash, independent prediction/pCTR, and bid-context enrichment.

![Where the OneRec side path merges into the main path](assets/diagrams/onerec/en/merge-paths.svg)

This topology provides isolation:

- before merge, OneRec does not consume the main DocWash and main-prediction candidate set;
- side-path results can be abandoned on failure without blocking the main path;
- the final merge combines and deduplicates at both AID and creative granularity before common reranking, constraints, and auction.

See [[en/recommender-systems/onerec/source-level-recall-to-reranking|Recall, Coarse Ranking, Fine Ranking, and Reranking]] for the complete execution semantics.

## 5. Cache Synchronization Semantics

Mixer's retrieval RPC and Ranking Fetch cooperate asynchronously:

1. GprHub creates a cache node for the request and position.
2. When the generative pipeline completes, it writes the ad list, `recall_path`, and experiment tags, then sets `fill_done`.
3. If Ranking Fetch reaches the node before completion, it waits until its timeout.
4. On success, it returns ads and context; failure distinguishes key miss, timeout, and empty result.

“The pipeline itself is slow” and “Fetch waits too long” must be monitored separately. The former is diagnosed through segmented GprHub RPC/Op metrics, while the latter uses `fetch_wait_latency` and cache-key hit rate.

## 6. Critical Invariants

- The main and GprHub side paths may run concurrently in production; inside GprHub, one position selects either the generative or legacy-GPR control implementation.
- Encoder, Decoder, RQ-VAE codebook, and KV-creative versions agree.
- SID packing and unpacking use the same number of levels and bit width at each level.
- Decoder scores remain aligned with their SIDs.
- Quota truncation occurs after Z-order/shard merging so that one shard cannot dominate.
- Deduplication occurs by AID after TID→AID mapping; multiple creatives for one ad must not be counted as distinct ads.
- Cache write volume and Ranking Fetch return volume can be cross-checked.
- AID overlap, creative overlap, and final global-creative deduplication can be reconciled independently.
- OneRec failure must not escalate into failure of the main path.
