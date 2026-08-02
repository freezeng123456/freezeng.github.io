---
title: OneRec Online Retrieval
description: The six-node GprHub DAG, SID/TID funnel, and common downstream path
lang: en
translation: 推荐系统理论/onerec/在线召回实现
tags:
  - OneRec
  - Online-Serving
  - Internal-Review
---

## Service Entry and Routing

Mixer creates concurrent tasks for the main Retrieval Proxy and the GprHub side path. After the side path is selected, GprHub reads a position-level generative switch:

- enabled: call GprHub's generative-retrieval RPC;
- disabled: call the legacy-GPR control RPC within the side path;
- side path not selected or call failure: the main Retrieval Proxy still proceeds.

Production and gray environments do not stop the main path on a GprHub hit; skipping main retrieval is reserved for isolated simulation under dedicated conditions. GprHub splits the request by position and executes asynchronously. `BuildSinglePosRequest` must rewrite the corresponding `pos_id` and `exp_param`; otherwise, multiple positions will reuse one another's parameters.

## Shared `GenerativeData` Context

The six operators in the new DAG pass state through `GenerativeData`, avoiding repeated request parsing. Its core fields fall into:

- request: position, experiment parameters, policy, and quota;
- user: model-input context and behavior sequence;
- model: encoder embedding, model version, and data name;
- search: beam, tree level, uniqueness filter, SIDs, and scores;
- inventory: TIDs, AIDs, and creative list;
- monitoring: experiment tag, retrieval path, and segmented status.

This is a strongly typed pipeline. A field change must be checked at its producer, consumer, cache serialization, and monitoring sites.

## Six Generative Operators

### 1. Entry

The position-level experiment parameters and quota plan provide:

- tree/trie/KV version;
- `pam_beam_width`;
- `beam_search_topk`;
- beam-sharding ratio;
- retrieval strategy and output quota;
- multilevel/any-recall mode.

The new path defaults `pam_beam_width` to 180 and requires it not to exceed the model/Decoder maximum of 180. Although the shared structure has an initializer of 200, Entry overwrites it at runtime with the position value or the safe value of 180.

### 2. Datahub

`QueryGenerativeToken` retrieves the generative behavior sequence and context. It does not redundantly fetch all Mixer user features; it prepares a model-specific sequence. Monitoring must distinguish an empty source sequence, emptiness after business filtering, and RPC failure.

### 3. Encoder

The Encoder assembles a multistage model request, runs the first-stage encoding, and returns:

- a shared user embedding;
- the model version;
- model context needed by the Decoder.

Encoding is reusable within one request. Repeating it for multiple searches wastes a principal source of compute.

### 4. Decoder

The Decoder performs level-wise beam search and returns SIDs with cumulative scores. Its request carries:

- shared user embedding;
- `data_name` / model version;
- multilevel tree parameters;
- beam width/top-k;
- uniqueness filter.

Tokens at each level are packed into `sid_key` according to the protocol. See [[en/recommender-systems/onerec/beam-search\|Beam Search]] for the search semantics.

### 5. KV Creative

The system queries multiple KV shards by SID and obtains TID lists. Shard outputs are interleaved in Z-order to prevent early results from one shard from consuming the entire quota, then truncated.

If an explicit quota is absent, the code has a fallback limit of 1000. This fallback prevents unbounded returns and must not replace the policy-level quota.

### 6. Creative List Build

TIDs are mapped to AIDs and creative objects. The system deduplicates by AID, retains generation scores and substrategies, and emits a common `MultiAdCreativeList`. From this point, generated candidates can enter the existing common ad-retrieval pipeline.

## Common Downstream Path

```text
AdListMerge
  → OutputQuotaControl (+ diversity)
  → RecallExit
  → CreativePropertyFilter
  → CreativeServer
  → Sharding
  → Scoring
  → ScoringExit
  → RetrievalCache
```

The generative implementation and the legacy-GPR control inside GprHub have broadly isomorphic common downstream graphs, forming the basis for A/B attribution. A small number of nodes still differ between the current static DAGs and require a structured diff during review:

1. compare the generative exit volume before the paths meet;
2. compare quota truncation, property filtering, creative serving, and coarse-ranking input/output segment by segment;
3. finally compare cache writes with Ranking Fetch returns.

“New versus legacy” here refers only to implementations inside GprHub. The main Retrieval Proxy is a concurrent production baseline and requires separate main-only, side-only, and overlap accounting.

## SID→TID→AID Funnel

| Stage            | Count                            | Meaning of a decline                |
| ---------------- | -------------------------------- | ----------------------------------- |
| Decoder          | Unique SID count                 | Beam/model output                   |
| KV pre-Z-order   | Total TIDs from all shards       | SID coverage and KV hits            |
| KV post-Z-order  | TIDs after truncation            | Quota/interleaving impact           |
| Creative Build   | TIDs successfully mapped to AIDs | Missing forward-index mapping       |
| AID dedupe       | Unique ad count                  | Multiple creatives merged           |
| Common filtering | Ads/TIDs after filtering         | Inventory and policy constraints    |
| Scoring top-n    | Coarse-ranking output            | Ranking quota                       |
| Cache/Fetch      | Final ad count                   | Asynchronous write/read consistency |

When empty results rise, the first collapse in the funnel must be located instead of examining only the final `no_ad` metric.

## Multiple Positions and Concurrency

- Positions can execute asynchronously.
- Each position requires independent parameters, cache node, and monitoring tag.
- A request-level shared Encoder can be reused only when schemas and models are identical.
- The merge callback must support partial success rather than discarding other positions when one fails.
- The cache key must contain sufficient request and position dimensions to prevent Ranking Fetch from reading the wrong result.

See [[en/recommender-systems/onerec/source-level-recall-to-reranking|Recall, Coarse Ranking, Fine Ranking, and Reranking]] for Scoring, Ranking Fetch, dual fine-ranking branches, and merge behavior.
