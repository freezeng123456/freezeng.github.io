---
title: OneRec Recall, Coarse Ranking, Fine Ranking, and Reranking
description: A source-driven walkthrough from concurrent Mixer routing to Ranking merge and common reranking
lang: en
translation: 推荐系统理论/onerec/召回粗排精排重排源码链路
tags:
  - OneRec
  - Online-Pipeline
  - Source-Walkthrough
  - Internal-Review
---

> [!warning] Redacted review draft verified against source
> This page follows the control flow, data flow, and failure semantics of the currently accessible serving branch. Internal repository paths, endpoints, experiment IDs, real fields, and code excerpts are omitted for a public-site safety boundary. Module names are retained only where they are needed to identify an execution stage. The implementation-level version should still move to a private knowledge base after review.

## 1. Executive Conclusions

The current online topology can be summarized as follows:

1. Mixer builds one retrieval request and schedules both the main Retrieval Proxy and the GprHub side path. In production and gray environments, the paths normally coexist.
2. GprHub has a second A/B decision inside the side path. With the generative switch on, it enters the six-operator OneRec DAG; with the switch off, it uses the legacy GPR control implementation.
3. The generative RPC first returns a side-path handle to Mixer. Retrieval, filtering, and coarse ranking continue asynchronously, and their result is written to a cache isolated by request and position.
4. Ranking starts the side-path Fetch early. The main and side paths independently bind, wash, predict, compute pCTR, and enrich bid context.
5. `MergeMultiAds` combines the paths at both ad and creative granularity. Prediction post-processing, style reranking, template reranking, the dynamic strategy graph, auction logic, and response construction are shared after that point.
6. A OneRec-specific relevance score is an optional factor in common reranking when its global gate, position experiment, scene, and threshold all match. It does not replace the reranker.

The design is a concurrent, late-merge, fail-open side path. A OneRec timeout, cache miss, or empty response does not have to block the main result.

## 2. Three Comparison Baselines

“Main path,” “legacy GPR,” and “isolated simulation” refer to different branches in the code and must be analyzed separately.

| Comparison                                 | Subjects                                           | Execution relationship                                                 | Question it answers                                                                                 |
| ------------------------------------------ | -------------------------------------------------- | ---------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| Production main path vs GprHub side path   | Retrieval Proxy and GprHub                         | Usually concurrent in production/gray; merged late in Ranking          | How much incremental final value does OneRec contribute?                                            |
| New vs legacy implementation inside GprHub | Six-operator generative DAG and legacy GPR control | Mutually exclusive within one position's side-path action              | How do the new generative implementation's quality, latency, and stability compare with legacy GPR? |
| Isolated gprsim validation                 | Selected side path only                            | Main retrieval may be skipped; a simulation-specific route may be used | Can a single path be replayed and traced deterministically?                                         |

The statement “new and legacy paths are mutually exclusive for a position” applies only to the **implementation choice inside GprHub**. It does not describe the production main and side paths, which run together.

![The complete online topology of the main and OneRec side paths](assets/diagrams/onerec/en/full-online-path.svg)

## 3. End-to-End Request Timeline

### 3.1 Mixer: one request, two sibling tasks

The Mixer DAG executes `BuildRetrievalRequest` and then starts two sibling tasks:

- the main task sends the request to Retrieval Proxy;
- the GprHub task evaluates global capability, position-level experiment, and routing conditions before creating a side-path action.

The GprHub task then reads the generative switch:

- enabled: call the generative entry point;
- disabled: call the legacy GPR control entry point;
- side-path experiment not matched: create no GprHub action while the main path proceeds normally.

`UpdateLocalCache` waits for the main retrieval task and the GprHub entry call, after which retrieval-response processing and `RankingAction` run. A completed GprHub entry call only means that the side-path handle has been returned; the generative pipeline may still be running in the background.

### 3.2 GprHub: acknowledge first, produce candidates per position asynchronously

For a request containing multiple positions, the generative service:

1. splits it into single-position subrequests;
2. retains the matching `pos_id`, experiment parameters, and required context in each subrequest;
3. writes the side-path handle needed by Ranking Fetch into the entry response;
4. completes the entry RPC immediately;
5. launches one independent asynchronous task per position for retrieval, common filtering, and coarse ranking;
6. writes the completed result to Retrieval Cache and wakes any waiting Ranking Fetch.

Entry-RPC latency and candidate-production latency are therefore separate metrics. Mixer-to-GprHub RPC P99 alone omits Decoder, KV, coarse ranking, and Fetch wait time.

### 3.3 Ranking: continue the main path while Fetch waits in parallel

After parsing the request and merging base contexts, Ranking starts the side-path Fetch early. Main-path ad-table binding, DocWash, and prediction continue concurrently. Fetch returns immediately if the cache is ready, waits on a condition variable within the position-level budget if production is still running, and yields an empty side path on key miss, timeout, or an empty result. None of those side-path outcomes must stop the main path.

## 4. Recall: Responsibilities of the Main and Generative Paths

### 4.1 Main Retrieval Proxy

The main path executes multiple recall sources and then performs quota control, inventory and targeting filters, creative processing, sharding, and retrieval-side Scoring in the same DAG. Its candidates enter Ranking through the regular Mixer response.

A GprHub experiment hit does not automatically disable the main path in production or gray. Only an isolated simulation mode skips the main task under dedicated conditions. Any review of online traces must first establish the environment so that a simulation topology is not generalized to production.

### 4.2 The six generative operators in GprHub

| Operator                   | Principal source-level action                                                   | Main product                            | First failure symptom                                              |
| -------------------------- | ------------------------------------------------------------------------------- | --------------------------------------- | ------------------------------------------------------------------ |
| Entry                      | Resolve position experiment, quota plan, trie/KV version, and search parameters | Shared `GenerativeData` context         | Route not matched, missing version, invalid width                  |
| Datahub                    | Retrieve the behavioral sequence and context needed by the generation model     | Model input sequence                    | RPC failure, empty source sequence, empty sequence after filtering |
| Encoder                    | Build the encoding request and compute a reusable user representation           | Shared user embedding and model context | Schema or model-version mismatch                                   |
| Decoder / AnyRecall Search | Generate token paths level by level and compute SID scores                      | SID set and SID→score mapping           | Beam out of range, empty output, score misalignment                |
| KV Creative                | Query sharded KV by SID, collect TIDs, interleave, and truncate                 | TIDs, SID/TID mapping, model scores     | KV-version mismatch, shard miss, excessive quota loss              |
| Creative List Build        | Map TID→AID, deduplicate by AID, and build the common ad list                   | `MultiAdCreativeList`                   | Missing forward-index entry, complete deduplication                |

The Decoder packs level-wise tokens into SID keys under a fixed protocol. The current implementation combines fixed-width levels into a 64-bit key, and a complete path uses its final-level score as the SID score. The model, trie, KV, and packing protocol require a version-bound release. A drift in level count or bit width appears immediately as KV misses.

KV Creative first orders results by shard cost and then interleaves TID groups by SID before applying the quota. This prevents one SID or one shard from consuming the full budget. Creative List Build then performs TID→AID mapping and AID deduplication while preserving generation scores, unified strategy, and substrategy on the common ad object.

### 4.3 New and legacy DAGs within GprHub

The generative DAG reuses most of the Retrieval Proxy's common graph and replaces the legacy remote-GPR node with the six generative operators. A static diff of the current graphs still shows a few differences: the main graph retains additional AnyRecall/LTR, GNN-creative, and early-trace nodes that are not fully mirrored in the generative graph.

This drift must remain an explicit review item:

- intentional differences need a business rationale and a defined metric boundary;
- accidental graph drift needs a machine-readable DAG diff and tests;
- new-versus-legacy attribution should compare the generative exit, common-filter exit, and Scoring exit rather than only the final ad count.

## 5. Coarse Ranking: the Actual Boundary of Retrieval-Side Scoring

After generated candidates are materialized as a common ad list, they continue through:

```text
AdListMerge
  → targeting and inventory filters
  → offline / live / allowlist filters
  → OutputQuota and diversity quota
  → CreativeFilter / PropertyFilter
  → CreativeServer
  → Sharding
  → Scoring subflow
  → ScoringExit
```

The Scoring subflow then:

1. initializes the coarse-ranking context;
2. retrieves or assembles cross features;
3. builds the scoring request;
4. calls the coarse-ranking model;
5. lets `ScoringExit` package top-n candidates and strategy context into the retrieval response.

This is the OneRec side path's coarse-ranking stage. It precedes the GprHub cache write and is independent of fine-ranking prediction inside Ranking. Review metrics must cover:

- Scoring input and output counts;
- top-n truncation;
- cross-feature missing rate;
- model failure, timeout, and empty response;
- serialization and cache-write latency between `ScoringExit` and cache completion.

## 6. Asynchronous Cache: Producer–Consumer Handshake

### 6.1 Cache-node state

Each position creates an independent cache node when it enters generative retrieval. The node contains:

- a redacted request-and-position composite key;
- processing start and end times;
- experiment tag and recall path;
- the post-coarse-ranking retrieval response;
- `fill_done` state;
- a fiber mutex and condition variable.

The cache uses a sharded map to reduce lock contention. It has configurable normal expiry for filled nodes and hard expiry for a producer that never finishes. A successful Fetch removes its node proactively; background cleanup covers absent consumers and abnormal exits.

### 6.2 Producer write order

1. Create and insert the node.
2. Run the single-position retrieval DAG.
3. Under the node lock, write the response, end time, experiment tag, and recall path.
4. Set `fill_done=true` last.
5. Release the lock and notify waiters.

Writing `fill_done` last prevents Ranking from observing completion while the response is only partially populated.

### 6.3 Ranking Fetch read order

1. Reconstruct the same composite key.
2. Return an empty side path immediately on cache miss.
3. On a hit that is not ready, wait within the position-level timeout.
4. Recheck the completion predicate after wake-up.
5. Copy the response and delete the node on success.
6. Report key miss, wait timeout, empty response, successful ad count, and end-to-end latency separately.

### 6.4 Race windows to observe independently

| Window                                   | Symptom                                    | Evidence                                            |
| ---------------------------------------- | ------------------------------------------ | --------------------------------------------------- |
| Ranking arrives before the producer      | Fetch wait rises but may succeed           | Cache hit followed by `fill_done`                   |
| Producer finishes before Ranking         | Fetch returns quickly                      | Cache age with near-zero wait                       |
| Node hard-expires before Ranking         | Key miss                                   | Production time, cleanup count, Fetch arrival time  |
| Producer exits abnormally                | Wait until timeout or hard expiry          | Node remains incomplete; last successful operator   |
| Composite-key collision across positions | Candidates or tags from the wrong position | Position validation and response-tag reconciliation |

## 7. Fine Ranking: Two Branches inside Ranking

### 7.1 Fetch overlaps main-path preprocessing

`FetchGprAction` starts after the base contexts are merged; it does not wait for full main-path DocWash. It fetches each position concurrently and reconstructs the compressed coarse-ranking response as native Ranking candidates.

Side-path binding enriches those candidates with ad-table data, bid information, original recall sources, creative model data, and strategy data. It also marks GPR/OneRec origin and side-only state, which allows the merge to distinguish main-only, side-only, and overlapping objects.

### 7.2 DocWash is asymmetric

| Capability                                         | Main DocWash                            | GPR/OneRec DocWash                                                    |
| -------------------------------------------------- | --------------------------------------- | --------------------------------------------------------------------- |
| Candidate split and timeout truncation             | Full                                    | Basic split                                                           |
| Contract, RTA, same-ad, and material filters       | Full multistage chain                   | Lightweight chain that relies on part of the retrieval-side filtering |
| Shuffle, aggregation-page, and special-ad policies | Applied according to main configuration | Not fully replicated                                                  |
| Creative flatten and merge-back                    | Full                                    | Dedicated flatten/get/merge chain                                     |
| Top-n and bucket construction                      | Full main implementation                | Preserves the side-path candidate set                                 |

This asymmetry affects experimental interpretation. If removal rates differ, the investigation must distinguish retrieval-side filtering, Ranking DocWash, and missing binding fields. A main-to-side field-contract test and comparable removal-reason counters are required.

### 7.3 Separate prediction and pCTR paths

The main branch runs its regular `PredictionRequest → PctrAction`. The side branch has a GPR-specific prediction request and pCTR action. Each branch also enriches RTA, bid, ecosystem, eGMV, pacing, and style context before meeting around common bid/ecosystem merge operations.

OneRec's fine-ranking increment should therefore be decomposed into:

1. coverage from incremental candidates;
2. ordering changes from side-specific prediction scores;
3. additional model and strategy information retained on overlapping candidates after merge.

## 8. Merge and Reranking: Concrete Late-Merge Semantics

![Candidate binding, merge, and common reranking for OneRec inside Ranking](assets/diagrams/onerec/en/ranking-merge-detail.svg)

### 8.1 Ad-level merge in `MergeMultiAds`

The merge operator first indexes main candidates by AID and then processes every side-path candidate:

- unseen AID: append the side-path ad as an incremental candidate;
- existing AID: keep one ad object while merging main/GPR recall sources, policies, and quota metadata;
- overlapping ad: clear the side-only flag and retain both normal and GPR attribution.

The final ad count is consequently not the sum of the two input counts. Monitoring needs `main_count`, `gpr_count`, `overlap_aid_count`, and `merged_unique_aid_count`.

### 8.2 Creative-level merge and global deduplication

Within an ad, a creative is matched by creative size and creative ID:

- merge GPR model, policy, and score data for the same creative;
- append a creative that exists only on the side path;
- after ad-level merge, deduplicate again by the global creative identifier;
- on conflict, prefer the greater bid-adjusted value and then the more recent creative.

This stage can change the final creative owner. Missing bid, ecosystem, or rerank-factor fields on side candidates must be monitored because they can bias the deduplication winner.

### 8.3 Common reranking

After merge, execution proceeds through:

```text
MergeMultiPredictionRecommendContext
  → PredictionPostProcess
  → UpdateAdxRecommendContext
  → RerankingPreProcess
  → StyleReranking
  → TemplateReranking
  → Dynamic Reranking Graph
  → RerankingPostProcess
  → AdLoad / strategy execution
  → Response
```

When a global gate, position experiment, and designated scene all match, OneRec relevance optimization reads a generative experience score, applies a threshold, clamps it to $[0,1]$, and uses it as a relevance factor. It is an optional input to the shared reranking graph. Experiments must bucket this factor separately so that generative-recall gain is not conflated with reranking-factor gain.

## 9. Key Parameters and Their Source-Level Effect

| Layer                | Parameter                                      | Resolution and effective behavior                                                 | First-order effect                            |
| -------------------- | ---------------------------------------------- | --------------------------------------------------------------------------------- | --------------------------------------------- |
| Mixer routing        | GprHub enable / call gate                      | Global capability, position experiment, and routing must all match                | Side-path QPS                                 |
| A/B inside side path | `mixer_enable_generative_recall`               | Selects generative RPC or legacy-GPR control RPC                                  | Retrieval implementation and `recall_path`    |
| Generative search    | `gpr_pam_beam_width`                           | Entry resolves it; effective fallback is 180 and the current Decoder limit is 180 | Decoder work, SID count, P99                  |
| Final SID truncation | `gpr_beam_search_topk`                         | Takes effect when positive; no production value should be inferred when absent    | Decoder-exit SID count                        |
| Beam sharding        | `gpr_pam_beam_width_shard_ratio`               | Resolved per position and must agree with the total beam budget                   | Shard load and candidate distribution         |
| Tree and KV          | trie level / KV-creative version               | From model/trie configuration and must match the SID protocol                     | KV hit rate                                   |
| Candidate uniqueness | uniqueness filter                              | Supplied by quota/policy plan                                                     | Duplicate SID rate                            |
| TID quota            | policy quota                                   | Quota plan is authoritative; 1000 is only a missing-value safety cap              | Post-KV TID count and coarse-rank cost        |
| Ranking Fetch        | `ranking_enable_fetch_gpr`                     | Controlled jointly by service capability and position experiment                  | Whether the side path enters Ranking          |
| Fetch wait           | `gpr_server_hub_timeout_ms` / service fallback | Position value takes precedence and must fit the total main-path budget           | Timeout, side coverage, Ranking P99           |
| Candidate binding    | `gpr_bind_adtable_batch_size`                  | Current code fallback is 10                                                       | Binding batch count, throughput, tail latency |
| Rerank enhancement   | OneRec relevance gate / threshold              | Requires global gate, position experiment, and designated scene                   | Relevance factor in common reranking          |

The source also contains a structure initializer of 200 for beam width, but Entry overwrites the effective value with the position parameter or the safe fallback of 180. Online analysis must report the resolved value rather than inferring runtime behavior from an inactive initializer. See [[en/recommender-systems/onerec/beam-search|Beam Search]] and [[en/recommender-systems/onerec/parameter-reference|Parameter Reference]] for related details.

## 10. Stage-by-Stage Comparison

| Stage             | Main path                    | OneRec / GprHub side path                             | Fair comparison                                                   |
| ----------------- | ---------------------------- | ----------------------------------------------------- | ----------------------------------------------------------------- |
| Trigger           | Default main retrieval task  | Global + position experiment + route                  | Same request, position, and time window                           |
| Recall            | Multi-source Retrieval Proxy | Six-operator generative DAG or legacy-GPR control     | Report “main vs side” separately from “new vs legacy inside side” |
| Filters and quota | Main retrieval DAG           | Cloned common GprHub DAG with small graph differences | Per-node input/output and removal reasons                         |
| Coarse rank       | Main Scoring subflow         | Side Scoring subflow                                  | Feature completeness, top-n, model version, P99                   |
| Ranking handoff   | Regular retrieval response   | Async cache + Fetch                                   | Request-to-cache, wait, miss, timeout                             |
| Candidate binding | Main ad-table chain          | Dedicated binding after Fetch                         | AID/TID/creative field coverage                                   |
| DocWash           | Full main flow               | Lightweight dedicated flow                            | Comparable reasons and ownership of prior filters                 |
| Fine rank         | Main prediction / pCTR       | Side prediction / pCTR                                | Calibration, missing features, inference latency                  |
| Merge             | Baseline candidates          | Merge AID overlaps; append side-only candidates       | Overlap AID/creative and field conflict                           |
| Rerank            | Shared with side path        | Shared with main; optional OneRec factor              | Factor-gate buckets and common-policy version                     |
| Failure           | Affects primary response     | Miss/timeout/empty fails open                         | Stable main-path success rate                                     |

## 11. End-to-End Funnel and Latency Ledger

### 11.1 Candidate funnel

```text
request hits GprHub
  → valid Datahub sequence
  → Decoder unique SIDs
  → KV pre-interleave TIDs
  → KV post-quota TIDs
  → successful TID→AID mapping
  → AID deduplication
  → candidates after common filtering
  → Scoring top-n
  → cache fill
  → Ranking Fetch return
  → GPR DocWash output
  → GPR fine-rank output
  → Merge side-only / overlap
  → common-rerank winner
  → final exposure
```

Every edge requires absolute volume, conversion rate, failure reason, and position bucket. Final exposure alone conflates model, inventory, cache, and Ranking failures.

### 11.2 Latency ledger

Side-path readiness can be written as:

$$
T_{\text{side-ready}} = T_{\text{datahub}} + T_{\text{encoder}} + T_{\text{decoder}} + T_{\text{kv}} + T_{\text{filter}} + T_{\text{coarse}} + T_{\text{cache-write}}.
$$

The Fetch wait is approximately:

$$
T_{\text{wait}} = \max\left(0,\;T_{\text{side-ready}}-T_{\text{ranking-arrival}}\right),
$$

and is capped by the position-level Fetch timeout. Main-path Ranking overlaps asynchronous GprHub production, so summing component P99 values does not produce user-visible P99. A request timeline, critical path, and wait-time distribution are all required.

## 12. Source-Verification Index

The table retains reviewable logical anchors without exposing internal paths.

| Logical anchor           | Verified behavior                                                                      | Documentation conclusion                                       |
| ------------------------ | -------------------------------------------------------------------------------------- | -------------------------------------------------------------- |
| Mixer DAG                | Main retrieval and GprHub depend on the same request-build node                        | Both paths coexist in production                               |
| GprHub task              | The generative switch changes only the side-path RPC type                              | New/legacy exclusion is internal to the side path              |
| Retrieval Proxy guard    | Only isolated simulation skips the main task under dedicated conditions                | Simulation topology cannot be generalized to production        |
| Generative service entry | Split positions, return entry response first, continue in background fibers            | Entry-RPC and candidate-ready latencies are distinct           |
| Retrieval cache          | Sharded map, completion predicate, notification, success deletion, normal/hard cleanup | Fetch is a synchronous consumer of an asynchronous producer    |
| Generative DAG           | Six operators replace legacy GPR and then enter common filters and Scoring             | The side path includes complete recall and coarse ranking      |
| Ranking Fetch            | Starts before full main DocWash completes and runs per position                        | Fetch wait overlaps main preprocessing                         |
| GPR Ranking branch       | Dedicated binding, lightweight DocWash, separate prediction/pCTR                       | The side path has independent fine rank and asymmetric washing |
| Multi-path merge         | AID merge, creative merge, global creative dedup                                       | Input candidate counts cannot simply be added                  |
| Common reranking         | Style, template, and dynamic policies are shared; OneRec factor is conditional         | Retrieval gain and rerank-factor gain need separate buckets    |

## 13. Highest-Priority Review Closures

1. **DAG drift:** build a machine-readable node/dependency diff and enumerate allowed differences.
2. **DocWash asymmetry:** add field-contract and removal-reason reconciliation, confirming that mandatory safety and business filters remain covered.
3. **Timeout budget:** recompute Fetch timeout from current per-position P95/P99 while protecting the main-path deadline.
4. **Effective-parameter visibility:** report resolved beam, top-k, quota, trie/KV version, Fetch timeout, and model version for every request.
5. **Multi-position isolation:** add end-to-end tests for parameters, cache keys, tags, timeouts, and partial success.
6. **Merge-field completeness:** monitor bid, ecosystem, model score, policy quota, and rerank-factor missing rates on overlapping ads.
7. **Gain decomposition:** retain main-only, side-only, overlap, and OneRec-rerank-factor on/off dimensions.
8. **Environment labeling:** mark production, gray, and gprsim at the trace root so isolated behavior cannot leak into production postmortems.

See [[en/recommender-systems/onerec/stability-and-observability|Stability and Observability]] for dashboards and alert design.
