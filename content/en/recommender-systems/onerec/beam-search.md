---
title: OneRec Beam Search
description: Beam width, level-wise retention, top-k, graph buckets, and tuning tradeoffs
lang: en
translation: 推荐系统理论/onerec/beam-search
tags:
  - OneRec
  - Beam-Search
  - Parameters
  - Internal-Review
---

> [!important] The most easily misunderstood parameter
> `beam width` is neither “the final number of returned ads” nor “the GPU graph bucket.” It is the maximum number of prefix states retained at each semantic-ID level. Final SID count is further changed by top-k, deduplication, KV quota, SID→TID expansion, and downstream filtering.

## Search over a Multilevel Semantic Tree

The model predicts $L$ token levels. At level $l$, it expands the prefixes retained from the previous level:

$$
\mathcal C_l
=\{(s_{1:l-1},c):s_{1:l-1}\in\mathcal B_{l-1},c\in V_l\},
$$

accumulates a log-score for every candidate, and selects the best $B_l$:

$$
\mathcal B_l=\operatorname{TopK}_{B_l}(\mathcal C_l).
$$

`UpdateBeamState` must inherit the parent prefix, append a token, accumulate the score, and maintain length. Only the final complete token sequence is packed into an SID.

## Four Easily Confused Widths

### 1. Request beam width

`pam_beam_width` is the primary width passed from a strategy or experiment to the Decoder. In the new AnyRecall refactoring path:

- default fallback: 180;
- hard constraint: at most the model's exported maximum width of 180.

The legacy retrieval engine has a default of 200. The two values belong to different call paths and model constraints. A legacy configuration of 200 does not justify changing the new path to 200.

### 2. Per-level beam width

The Decoder supports a different width $B_1,\ldots,B_L$ at each level. A single `pam_beam_width` may be expanded to a constant sequence or combined with level-specific tree configuration.

One design retains more semantic clusters in shallow levels and narrows them in deeper levels. Another uses a narrow first level and allows more fine-grained branches later. The array produced at the concrete callsite is authoritative.

### 3. `beam_search_topk`

This parameter truncates output at the final or another specified stage. When it is $>0$, the final retained count can override the default beam output; otherwise, the beam width is typically reused.

```text
beam width = number of prefixes retained during search
top-k      = number of complete SIDs retained after search
```

Reducing top-k lowers KV-query and downstream costs but cannot recover Decoder computation already performed at earlier levels.

### 4. GPU graph/bucket width

Inference engines often precapture execution graphs for several supported widths. Values such as `90, 200`, or larger-width lists can denote **reusable graph buckets**, not the search width requested by the current strategy.

At runtime, the system generally selects the smallest bucket that can hold the request and pads to it. A request width of 180 that lands in bucket 200 may execute kernels with a shape of 200. Latency and memory should therefore be estimated from the bucket rather than the logical width of 180.

## Compute and Memory Cost

If the vocabulary at level $l$ has size $|V_l|$, naive expansion has approximate cost:

$$
O\!\left(\sum_{l=1}^{L} B_{l-1}|V_l|\right).
$$

Memory grows mainly with the largest active beam, hidden/cache state per beam, and candidate logits:

$$
M\approx O(B_{\max}d_{\text{state}})+O(B_{\max}|V_l|).
$$

Actual kernels reduce constants through fused top-k, shared prefixes, and graph capture, but increasing width generally raises Decoder compute, KV query volume, and downstream candidate volume approximately linearly.

## What Increasing Beam Width Can Do

Potential benefits:

- reduce search-truncation error;
- improve long-tail SID coverage;
- provide more candidates for KV and downstream ranking.

Potential costs:

- increased Decoder latency and memory;
- stepwise latency growth when the request moves to a larger graph bucket;
- more SIDs, KV fan-out, and network-response volume;
- wasted quota/filter work when low-scoring SIDs hit many TIDs;
- higher ranking load without necessarily increasing the final ad count.

Beam width should therefore be evaluated by whether marginal SIDs become marginal valid ads, not only by Decoder top-k metrics.

## Relationship to Shards and Quota

Beam-sharding ratio parameters allocate the total search or return budget across parallel shards. Outputs are interleaved in the KV stage and truncated by quota. The allocation should satisfy:

$$
\sum_s B_s \gtrsim B_{\text{request}},
$$

but $B$ returns from each of $S$ shards do not automatically produce $S\times B$ final SIDs because deduplication, empty shards, and top-k reduce the count.

## Interpreting Current Configurations

| Scenario                        |          Observed value | Correct interpretation                                                     |
| ------------------------------- | ----------------------: | -------------------------------------------------------------------------- |
| New generative refactoring path |        fallback/max 180 | Default and model maximum for the request beam                             |
| Legacy retrieval engine         |             default 200 | Legacy-path default; not transferable across paths                         |
| One sample/experiment version   |           90, 200, etc. | Use the parameter-ID callsite to distinguish beam, top-k, and graph bucket |
| Execution-graph list            | Several discrete widths | Runtime buckets, not strategy values                                       |

## Tuning Experiment

Scan $B\in\{30,60,90,120,150,180\}$ together with final top-k and record:

1. Decoder P50/P95/P99 and memory;
2. SID count, unique SIDs, and mean score;
3. SID→TID hit rate and TIDs per SID;
4. TIDs/AIDs before and after quota;
5. coarse- and fine-ranking pass rates;
6. final impressions, value, and diversity;
7. actual hit ratio of each GPU bucket.

The selected operating point should be where valid-ad gain begins to saturate and latency has not yet jumped to the next bucket, rather than mechanically choosing the maximum of 180.

## Troubleshooting

| Symptom                                   | Possible cause                                                 |
| ----------------------------------------- | -------------------------------------------------------------- |
| Request fails immediately                 | Beam exceeds model maximum, or level/width array is invalid    |
| Latency jumps suddenly                    | A larger graph bucket is selected                              |
| SID count rises but ad count does not     | KV miss, duplicate SID, AID deduplication, or quota truncation |
| A smaller beam performs better            | Noisy low-score tail consumes downstream quota                 |
| Shard results are skewed                  | Uneven sharding ratios, Z-order, or shard inventory            |
| Same-named parameter differs across paths | Different default, unit, or callsite semantics                 |
