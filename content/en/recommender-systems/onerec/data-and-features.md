---
title: OneRec Data and Feature Pipeline
description: Online/offline samples, USER/CONTEXT routing, and sequence slots
lang: en
translation: 推荐系统理论/onerec/数据与特征
tags:
  - OneRec
  - Feature-Engineering
  - Internal-Review
---

## Online Path

Online logs pass through prefiltering and feature decoding, then enter OneRec sample generation after an attribution window:

```text
DatahubLog
  → PreFilter
  → LogToFeature
  → attribution delay
  → OneRecGenerateSequence
  → USD enrichment
  → SequenceExample / training data
```

The input is not a flat feature map. It combines several semantic sources:

- compressed user-feature data;
- feature lists;
- PageView/position;
- BizTimeline context;
- WUID long sequences and session completion;
- content/ad ID embeddings.

## USER and CONTEXT Routing

Features first enter two base entries according to their source and allowlist.

### USER

USER receives stable user profiles and the direct user stream. Only features allowed by `userIntFeatureIdSet` enter USER, while a small number of special features are explicitly rerouted to CONTEXT.

### CONTEXT

CONTEXT aggregates request, page, position, label, mix, and list features, filtered through the `contextFeature` set. It describes the circumstances of this decision rather than long-term user attributes.

A hard gate requires the first USER entry in initial training data to be nonempty. If allowlisting or decoding errors empty USER, the sample is rejected even when CONTEXT is complete.

## USD Sequence Completion

After base USER/CONTEXT construction, USD adds:

- the primary behavior sequence;
- auxiliary item/content sequences;
- in-session content sequences;
- session completion for specific business forms;
- ID-embedding completion;
- time truncation and maximum-length limits.

An important implementation fact is that final entries are **not globally sorted by event time**. They are concatenated in a fixed slot order, for example:

```text
USER
CONTEXT
[optional] CONTENT_INTRA_SESSION
ITEM_AUX
CONTENT_AUX_*
CONTENT_AUX
ITEM / CONTENT main sequence
```

Within-slot order can remain temporal, but cross-slot array positions carry semantics. A training feature reader that treats the entries as one timeline and reorders them breaks the schema.

## Offline Path

Offline Burn reads base features from persisted WUID and user columns, then reuses USD sequence-construction logic. The goal is not identical code but identical output semantics:

| Dimension           | Online                     | Offline alignment                            |
| ------------------- | -------------------------- | -------------------------------------------- |
| USER allowlist      | Real-time stream filtering | Same-version set                             |
| CONTEXT source      | Live request               | Reconstruct the same fields from logs        |
| Delayed attribution | Fixed window               | Same cutoff                                  |
| Sequence truncation | Online USD                 | Same length/time rules                       |
| Slot order          | Fixed entry layout         | Exactly identical                            |
| ID embedding        | Online completion          | Same version or explicit missing-data policy |

Some external embedding completion may be disabled in the batch offline path. Missing markers or version fields must make that distinction visible; silently filling zeros and pretending equivalence is unsafe.

## Risk from Dual Implementations

A parallel sample-generation implementation currently targets the same semantics through different configuration and filtering entry points. The principal migration risk is not compilation, but the slow divergence of allowlists and slot definitions.

A golden sample should:

1. fix a redacted raw input;
2. have both implementations emit a `SequenceExample`;
3. compare entry type, count, order, feature IDs, sequence lengths, and cutoffs field by field;
4. maintain an explicit allowlist for accepted differences;
5. require every schema change to update both the reader and the golden sample.

## Common Fault Localization

| Symptom                                        | First checks                                             |
| ---------------------------------------------- | -------------------------------------------------------- |
| Sudden sample-volume decline                   | Empty-USER gate, prefiltering, attribution delay         |
| Online is good, offline is poor                | Allowlist version, sequence cutoff, embedding completion |
| One context-feature class becomes zero         | USER/CONTEXT routing sets                                |
| Sequence length looks normal but quality drops | Slot order, time units, cross-session concatenation      |
| Training runs but online Decoder returns empty | Schema/model-version mismatch                            |
