---
title: OneRec Stability and Observability
description: Rollout, rollback, segmented funnels, alerts, and comparison of old and new paths
lang: en
translation: 推荐系统理论/onerec/稳定性与可观测性
tags:
  - OneRec
  - Stability
  - Observability
  - Internal-Review
---

## Stability Objectives

| Objective    | Requirement                                                                          |
| ------------ | ------------------------------------------------------------------------------------ |
| Rollback     | Disable the GprHub side path or switch its internal implementation without a release |
| Comparison   | Main/side and new/legacy-within-side comparisons use separately aligned conventions  |
| Localization | Attribute empty/slow requests to an operator, RPC, or funnel stage                   |
| Degradation  | OneRec failure does not block the main path                                          |
| Auditability | Final parameters and model/KV versions are traceable                                 |

## Three Layers of Switches

1. **Side-path enable gate**: decides whether to create a GprHub action in addition to main Retrieval Proxy.
2. **Side-path implementation gate**: selects the generative RPC or legacy-GPR control RPC inside the GprHub action.
3. **GprHub service-capability gate**: determines whether the generative flow runner is created or enabled.

Rollback sequence:

- abnormal side-path metrics: disable the GprHub action and verify the continuing main result;
- generative-only failure: route the side path to the legacy-GPR control RPC;
- new-DAG code failure: disable the set-level switch and restart the affected set;
- Fetch timeout: stop expansion and compare pipeline P99 with the Ranking budget;
- model/KV version error: roll back both the model package and KV version; changing the beam alone is insufficient.

## Rollout Cadence

```text
Deploy (switch off by default)
  → 1%
  → 5%
  → 50%
  → 100%
  → legacy GPR control QPS inside the side path reaches zero
  → evaluate removal of the side-path legacy implementation
```

This cadence covers the generative replacement inside GprHub. Retiring the main Retrieval Proxy is a separate architecture decision and must not follow automatically from full side-path rollout.

At every step, compare at least:

- request/success/failure QPS;
- GprHub and Fetch P95/P99;
- cache-write and Fetch-return ad counts;
- SID→TID→AID funnel;
- quota/filter/coarse-ranking removal volume;
- business guardrails.

## Three Dashboards

### 1. Mixer → GprHub

Core metrics:

- retrieval request/success/failure;
- GprHub retrieval latency P99;
- cache ad count;
- request-to-cache latency;
- no-ad count.

This dashboard answers whether the GprHub action is selected and writes cache results on time. Main Retrieval Proxy success and candidate volume remain a concurrent baseline.

### 2. Generative pipeline

Each operator reports:

- entry count;
- latency/P99;
- failure count;
- empty output;
- input/output ad count.

Datahub, Encoder, Decoder, and KV separately report RPC request/failure/latency. Add the SID→TID funnel:

```text
sid_key
  → kv_tid_pre_zorder
  → decoder_tids
  → creative_build_tid
  → unique_aid
```

This dashboard identifies where the new path loses candidates or becomes slow.

### 3. Ranking Fetch

Core metrics:

- Fetch QPS, failures, timeouts, and key misses;
- Fetch latency and wait latency;
- response ad count / no-ad;
- end-to-end delay from completed cache write to successful Fetch.

This dashboard answers whether Ranking receives the side-path result in time.

## Two Comparison Layers

The first layer compares the generative and legacy-GPR control implementations inside GprHub. After the retrieval DAG's `AdListMerge`, compare by `recall_path` + `exp_tag`:

| Stage              | Focus                                                     |
| ------------------ | --------------------------------------------------------- |
| OutputQuotaControl | Truncated ad count and activation count                   |
| PropertyFilter     | TID input/output and removal rate                         |
| CreativeServer     | Empty shards, request/response ad count, preparation time |
| Scoring            | Input TIDs/AIDs, top-n, response size                     |
| Flow node          | Per-operator latency, input/output, failure/empty result  |

The second layer compares the production main path with the GprHub side path at Ranking merge. It adds:

- main-only, side-only, and AID overlap;
- creative overlap, side-added creatives, and global creative deduplication;
- main/GPR DocWash removal reasons;
- main/GPR prediction, pCTR, and bid-context field completeness;
- OneRec rerank-factor on/off buckets.

Attribution must proceed from upstream to downstream. If generative-exit volume already differs, the final ad difference cannot be assigned directly to coarse ranking. If retrieval-side volumes match but GPR DocWash loses candidates, Ranking fields and filter contracts are the next target.

## Review Thresholds

Thresholds must be derived from current online baselines and position-level budgets. Absolute values from old integration snapshots must not be carried forward:

| Guardrail                  | Condition to pause expansion                                            |
| -------------------------- | ----------------------------------------------------------------------- |
| GprHub cache-ready P99     | Exceeds the position's side-path budget throughout the alert window     |
| Segmented RPC failure rate | Significantly above the same-machine, same-position experiment baseline |
| Fetch timeout              | Side coverage drops while wait consumes the main-path deadline          |
| Empty result               | Significantly above the legacy-GPR control or same-experiment baseline  |
| Main-path health           | Regression in main success, P99, or final response volume               |

Thresholds must be bucketed by traffic, machine type, and position. Global averages obscure long-tail ad positions.

## Fault Tree

![The OneRec diagnostic path when no final ads remain](assets/diagrams/onerec/en/debugging-tree.svg)

## Gate for Removing the Side-Path Legacy Implementation

The legacy-GPR control inside GprHub can be retired only after all of the following hold:

1. the new path is at full traffic and remains stable throughout the observation window;
2. side-path rollback no longer depends on the legacy implementation, or an equivalent disaster-recovery path exists;
3. legacy-GPR control QPS remains zero;
4. new/old dashboards, alerts, and the on-call handbook are complete;
5. the legacy Datahub workflow, side-path `RemoteGprOp`, service routing, and graph configuration have an explicit removal list;
6. a reverse traffic-switching drill is completed before deletion.

Retiring the main Retrieval Proxy requires a separate proof that multi-source retrieval, full DocWash, field contracts, and disaster recovery have all been replaced. It is outside the scope of side-path legacy cleanup.
