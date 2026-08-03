---
title: Recommender Systems
description: From long-sequence modeling to the OneRec production pipeline
lang: en
translation: 推荐系统理论
tags:
  - Recommender-Systems
---

This section has two layers. Knowledge notes establish design frameworks for recommender methods, while production studies record module boundaries, important parameters, and reliability constraints in an end-to-end system.

## Knowledge Notes

- [[en/recommender-systems/long-sequence-modeling\|Long-Sequence Modeling for Recommender Systems: Methods and Design Choices]] explains how retrieval, compression, hierarchy, sparse interaction, and linear models trade information fidelity against latency and memory. It also provides a parameter reference, an experimental plan, and an engineering review checklist.

## Production Study

[[en/recommender-systems/onerec/index\|The OneRec production pipeline]] starts with user behavior sequences and advertising semantic IDs, then covers sample construction, model training, generative retrieval, creative mapping, ranking and merging, and online reliability.

> [!warning] Temporary internal-review draft
> The OneRec pages contain implementation-level information and are intended only for internal architecture review. The site excludes credentials, internal addresses, identities, and real samples. After review, this material should be moved to a private repository and thoroughly removed from the public Git history.

Suggested OneRec reading order:

1. [[en/recommender-systems/onerec/end-to-end-pipeline\|End-to-End Pipeline]]
2. [[en/recommender-systems/onerec/data-and-features\|Data and Features]]
3. [[en/recommender-systems/onerec/semantic-id-and-rq-vae\|Semantic IDs and RQ-VAE]]
4. [[en/recommender-systems/onerec/training-and-sample-generation\|Training and Sample Generation]]
5. [[en/recommender-systems/onerec/online-retrieval\|Online Retrieval]]
6. [[en/recommender-systems/onerec/beam-search\|Beam Search]]
7. [[en/recommender-systems/onerec/parameter-reference\|Parameter Reference]]
8. [[en/recommender-systems/onerec/stability-and-observability\|Stability and Observability]]
9. [[en/recommender-systems/onerec/review-checklist\|Review Checklist]]
