---
title: OneRec End-to-End Pipeline
description: Training loop, online topology, and the boundary between new and legacy paths
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

For each ad position, Mixer reads experimental parameters and selects one of two paths:

- New path: `GenerativeRecall` enters GprHub's `GENERATIVE_FLOW`;
- Legacy path: `RetrievalProxy` enters Retrieval Proxy, which calls the legacy GPR service through `RemoteGprOp`.

One position must use only one path to prevent duplicated retrieval and duplicated quota consumption. Requests containing multiple ad positions must be split at the service layer, with the correct `pos_id` and `exp_param` copied for each position instead of reusing those of the first position.

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

The path then meets the legacy path near `AdListMerge` and enters structurally equivalent quota, property-filtering, creative-service, and coarse-ranking subgraphs.

## 4. Ranking and Main-Path Merge

The main path continues through conventional retrieval, coarse ranking, and the primary Ranking stage. Within GprHub, OneRec writes results to Retrieval Cache. On the Ranking side, `FetchGenerativeX` waits for `fill_done`, reads the result, and executes independent side-path fine ranking.

![Where the OneRec side path merges into the main path](assets/diagrams/onerec/en/merge-paths.svg)

This topology provides isolation:

- OneRec does not consume the primary fine-ranking candidate quota;
- side-path results can be abandoned on failure without blocking the main path;
- the final merge still applies cross-path deduplication, common constraints, and the auction.

## 5. Cache Synchronization Semantics

Mixer's retrieval RPC and Ranking Fetch cooperate asynchronously:

1. GprHub creates a cache node for the request and position.
2. When the generative pipeline completes, it writes the ad list, `recall_path`, and experiment tags, then sets `fill_done`.
3. If Ranking Fetch reaches the node before completion, it waits until its timeout.
4. On success, it returns ads and context; failure distinguishes key miss, timeout, and empty result.

“The pipeline itself is slow” and “Fetch waits too long” must be monitored separately. The former is diagnosed through segmented GprHub RPC/Op metrics, while the latter uses `fetch_wait_latency` and cache-key hit rate.

## 6. Critical Invariants

- Each position selects exactly one of the new and legacy retrieval paths.
- Encoder, Decoder, RQ-VAE codebook, and KV-creative versions agree.
- SID packing and unpacking use the same number of levels and bit width at each level.
- Decoder scores remain aligned with their SIDs.
- Quota truncation occurs after Z-order/shard merging so that one shard cannot dominate.
- Deduplication occurs by AID after TID→AID mapping; multiple creatives for one ad must not be counted as distinct ads.
- Cache write volume and Ranking Fetch return volume can be cross-checked.
- OneRec failure must not escalate into failure of the main path.
