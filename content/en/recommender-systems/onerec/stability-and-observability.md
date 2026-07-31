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

| Objective    | Requirement                                                               |
| ------------ | ------------------------------------------------------------------------- |
| Rollback     | Instantly return a request to the legacy path without a release           |
| Comparison   | Old and new paths use the same experiment tags and downstream conventions |
| Localization | Attribute empty/slow requests to an operator, RPC, or funnel stage        |
| Degradation  | OneRec failure does not block the main path                               |
| Auditability | Final parameters and model/KV versions are traceable                      |

## Two Levels of Switches

1. **Request-level Mixer switch**: selects either the new GprHub generative RPC or legacy Retrieval Proxy.
2. **GprHub set-level switch**: determines whether the new `generative_flow_runner` is created or enabled.

Rollback sequence:

- abnormal online metrics: first route experimental traffic back to the legacy RPC;
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
  → legacy service QPS reaches zero
  → remove legacy path
```

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

This dashboard answers whether the old and new paths write healthy cache results.

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

## Comparing Old and New Paths Downstream

After `AdListMerge`, compare by `recall_path` + `exp_tag`:

| Stage              | Focus                                                     |
| ------------------ | --------------------------------------------------------- |
| OutputQuotaControl | Truncated ad count and activation count                   |
| PropertyFilter     | TID input/output and removal rate                         |
| CreativeServer     | Empty shards, request/response ad count, preparation time |
| Scoring            | Input TIDs/AIDs, top-n, response size                     |
| Flow node          | Per-operator latency, input/output, failure/empty result  |

Attribution must proceed from upstream to downstream. If exit volumes already differ before the paths meet, the final ad difference should not be attributed to coarse ranking. If they match before the merge but diverge after PropertyFilter, field or inventory alignment is the next target.

## Review Thresholds

These are current review suggestions and must not be hardened without validation against online baselines:

| Guardrail                  | Condition to pause expansion                                 |
| -------------------------- | ------------------------------------------------------------ |
| GprHub cache-write P99     | $>150$ ms for 5 minutes                                      |
| Segmented RPC failure rate | $>1\%$ for 3 minutes                                         |
| Fetch timeout              | Clearly higher than baseline; review draft uses a 5% change  |
| Empty result               | Significantly higher than the legacy path or same experiment |

Thresholds must be bucketed by traffic, machine type, and position. Global averages obscure long-tail ad positions.

## Fault Tree

![The OneRec diagnostic path when no final ads remain](assets/diagrams/onerec/en/debugging-tree.svg)

## Gate for Removing the Legacy Path

The legacy path can be retired only after all of the following hold:

1. the new path is at full traffic and remains stable throughout the observation window;
2. request-level rollback no longer depends on the legacy service, or an equivalent disaster-recovery path exists;
3. legacy-service QPS remains zero;
4. new/old dashboards, alerts, and the on-call handbook are complete;
5. the legacy Datahub workflow, Proxy `RemoteGprOp`, service routing, and graph configuration have an explicit removal list;
6. a reverse traffic-switching drill is completed before deletion.
