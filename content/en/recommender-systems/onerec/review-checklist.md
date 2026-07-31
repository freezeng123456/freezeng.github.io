---
title: OneRec Architecture Review Checklist
description: Entry points for reviewing code, parameters, data, stability, and security
lang: en
translation: 推荐系统理论/onerec/评审清单
tags:
  - OneRec
  - Review
  - Internal-Review
---

## High-Level Conclusion

OneRec closes the loop from semantic IDs, training samples, and a sequence model through online generative retrieval, the existing advertising pipeline, and feedback. The new path isolates retrieval and fine-ranking resources from the main path, remains fail-open as a side path before the merge, and applies common constraints downstream.

The principal review question is not whether the algorithm can generate SIDs, but whether the following engineering loops are reliable:

- training/online schema and model/KV versions;
- hierarchy of beam and candidate budgets;
- parameter isolation across ad positions;
- asynchronous cache and Ranking timeout budget;
- comparable monitoring and instantaneous rollback for old and new paths.

## Itemized Confirmation

### Data

- [ ] USER/CONTEXT allowlists have one version source.
- [ ] Online/offline slot order and cutoff have a golden diff.
- [ ] The empty-USER gate, delayed attribution, and sequence completion have scale metrics.
- [ ] Accepted differences between parallel sample-generation implementations are documented.
- [ ] Schema versions are included in sample and model metadata.

### Semantic IDs

- [ ] Codebook hash, levels, bit width, and pack version are bound to the release.
- [ ] Collision, code usage, dead codes, and KV coverage meet requirements.
- [ ] RQ-VAE checkpoint selection considers more than reconstruction.
- [ ] SID→TID version mismatch is directly identifiable on dashboards.

### Training

- [ ] Candidate pool, sample count, return count, and $M$ have a valid size hierarchy.
- [ ] Reward units, normalization, missing values, and versions are explicit.
- [ ] Concurrent reward RPCs do not concurrently write to a shared protobuf map.
- [ ] The model reader and online Encoder schema checksums agree.
- [ ] Frozen replay / causal control is reported separately from the main experiment.

### Online Retrieval

- [ ] New and old paths are mutually exclusive for one position.
- [ ] Parameters, cache nodes, and tags are independent across positions.
- [ ] Input/output/latency/failure are observable for all six operators.
- [ ] The SID/TID/AID/ad funnel is complete.
- [ ] Common downstream stages are compared by `recall_path`.
- [ ] Partial OneRec failure does not expand into main-path failure.

### Beam

- [ ] New-path `beam_width ≤ 180`.
- [ ] The legacy default of 200 does not enter the new path.
- [ ] Per-level beam, final top-k, and graph bucket are reported separately.
- [ ] Tuning considers valid-ad gain, P99, and memory together.
- [ ] Beam growth does not overwhelm KV quota/Z-order.

### Stability

- [ ] Request-level and set-level rollback have both been drilled.
- [ ] Fetch timeout exceeds an acceptable pipeline P99 without violating the main-path budget.
- [ ] Rollout guardrails are bucketed by position and machine type.
- [ ] Cache writes and Fetch returns can be cross-validated.
- [ ] Reverse traffic switching is drilled before retiring the legacy path.

## Code-Reading Index

To avoid exposing internal repository paths on the site, only module-level entry points are listed:

| Topic               | Module/class                                                                  |
| ------------------- | ----------------------------------------------------------------------------- |
| Generative DAG      | `GenerativeEntry/Datahub/Encoder/Decoder/KvCreative/CreativeListBuildOp`      |
| Service entry/cache | retrieval service implementation, `DoGenerativeRetrieval`, `FetchGenerativeX` |
| Legacy path         | `RemoteGprOp` and Retrieval Proxy flow                                        |
| Common downstream   | OutputQuota, PropertyFilter, CreativeServer, Scoring                          |
| Beam state          | multilevel handler, SID decoder, `UpdateBeamState`                            |
| RQ-VAE              | encoder, residual vector quantizer, trainer                                   |
| Sample generation   | `ExampleGenerator`, OneRec generator, RL-agent configuration                  |
| Data features       | LogToFeature, OneRec sequence generator, USD enrichment                       |

## Pre-Publication Security Check

- [ ] Search for all IPv4 addresses, domains, accounts, tokens, passwords, cookies, and keys.
- [ ] Search for experiment IDs, names, group identifiers, and real business fields.
- [ ] Inspect image metadata and SVG text layers.
- [ ] Inspect Git history, not only the current files.
- [ ] Confirm that OneRec content is authorized for public GitHub Pages.
- [ ] After internal review, migrate it to a private repository and rewrite public history.

> [!danger] Publication gate
> Even after technical redaction, OneRec remains implementation-level internal material. Business and security authorization is required before public release. Without authorization, it should be previewed only locally or deployed from a private repository.
