---
title: "Long-Sequence Modeling for Recommender Systems: Methods and Design Choices"
description: An industrial design space spanning retrieval, compression, hierarchy, sparse interaction, and linear sequence models
lang: en
translation: computational-mathematics/paper-notes/推荐系统/长序列建模
tags:
  - Recommender-Systems
  - Long-Sequence-Modeling
  - Sparse-Attention
  - Multimodal-Learning
---

> [!note] Scope
> Method details follow the public papers cited for each route. The route-level assessments establish a design framework; they are not universal judgments about any one method.

## Executive Summary

Long-sequence recommendation addresses a concrete problem. Highly active users can accumulate thousands or tens of thousands of interactions, yet an online model must still identify the long-term interests relevant to the current candidate under strict latency and memory budgets. Keeping every event reduces information loss, while simultaneously increasing attention compute, feature reads, cache demand, and operational risk.

Industrial approaches can be organized into five main routes:

1. **Retrieval** selects a small set of target-relevant events before detailed interaction;
2. **Compression** encodes a long history into a fixed number of interest representations, decoupling online cost from the original sequence length;
3. **Hierarchical modeling** processes local segments first and propagates long-range information through global representations;
4. **Sparse interaction** retains only high-value attention edges and directly controls the size of the interaction graph;
5. **Linear sequence models** replace full quadratic attention with state-space recurrences or mixers.

No route dominates across all settings. Retrieval emphasizes target relevance, compression captures stable long-term interests, sparse interaction preserves more event-level detail, hierarchy exploits temporal locality, and linear models scale more gently with length. A practical architecture today combines two complementary paths: **compressed representations maintain long-term memory, while target-aware semantic sparse interactions retain candidate-specific detail**.

![The design space of industrial long-sequence recommendation](assets/diagrams/long-sequence-recommendation/en/design-space.svg)

## 1. How the Problem Arises

Let the user history be

$$
H=(h_1,h_2,\ldots,h_L),
$$

and let $q$ represent the current candidate. Standard self-attention constructs an $L\times L$ interaction matrix, so its runtime and intermediate activations scale approximately as $O(L^2)$. If each candidate also interacts with the complete history, the additional cost normally grows at least linearly with $L$. An online system faces constraints that extend beyond this complexity expression:

- **Feature I/O:** IDs, context, and multimodal features for thousands of events must be read, assembled, and transferred;
- **Memory and caches:** activations, KV caches, candidate batches, and model parameters share device memory;
- **Tail latency:** history length typically has a long-tailed distribution, so exceptionally long samples can dominate P95 and P99 latency;
- **Training stability:** longer dependency paths, sparse topologies, and altered negative-sample structures affect optimization;
- **Online consistency:** offline indexes, interest vectors, and semantic adjacency must remain synchronized with online features and refresh promptly.

A long-sequence design must therefore answer three questions:

1. **Which information is retained?** Temporal order, semantic similarity, interest transitions, and infrequent long-term preferences require different levels of fidelity.
2. **When does the target enter?** The candidate can participate during retrieval, compression, sparse-graph construction, or final ranking.
3. **Where is cost paid?** Work may move to offline indexing, training preprocessing, online GPU kernels, or a lightweight service, but it still has to be accounted for.

## 2. Route-Level Overview

| Route              | Core operation                                            | Online-cost intuition                                    | Main advantage                                                 | Main risk                                                               |
| ------------------ | --------------------------------------------------------- | -------------------------------------------------------- | -------------------------------------------------------------- | ----------------------------------------------------------------------- |
| Retrieval          | Select $K_r$ events from $L$, then refine                 | Depends on $K_r$, plus retrieval                         | Strong target awareness, controllable cost, mature engineering | Prefiltering may lose transitions; coarse and exact stages can disagree |
| Compression        | Encode $L$ events into $M$ interest tokens                | Depends mainly on $M$ after compression                  | Reusable representation; online cost can be decoupled from $L$ | Lossy bottleneck may average away local target detail                   |
| Hierarchy          | Interact within local blocks, then aggregate globally     | Determined by block size, block count, and global tokens | Preserves temporal structure and local continuity              | Cross-time semantic jumps may fall into different blocks                |
| Sparse interaction | Compute only window, semantic-neighbor, or directed edges | Often targets $O(LK_s)$                                  | Retains more raw events and has a high capacity ceiling        | Sparse rules, indexes, and kernels require joint design                 |
| Linear model       | Aggregate with recurrent state or token/channel mixing    | Usually linear in $L$                                    | Length-friendly and memory-efficient                           | Target awareness and production-scale gains need further validation     |
| Dense end to end   | Preserve the sequence and all interactions                | Approximately $O(L^2)$                                   | Most complete information paths and a direct objective         | Highest compute, memory, and deployment cost                            |

The evolution reflects a shift in emphasis. Early systems focused on locating candidate-relevant events in a very long history. Fixed-slot interest compression and hierarchical structures then reduced representation cost. Recent work combines semantic adjacency, directed attention, and system-level kernel optimization so that a model can retain more detail at longer lengths.

## 3. Retrieval: Filter First, Interact Precisely

Retrieval decomposes long-sequence modeling into two stages. A coarse stage scores the history against candidate $q$ and selects

$$
S(q,H)=\operatorname{TopK}_{i\in[1,L]} s(q,h_i).
$$

The exact model processes only the $K_r$ events in $S$. Retrieval speed matters, but alignment between the coarse score and the final ranking objective matters just as much.

### 3.1 SIM: The Two-Stage Pattern

[SIM](https://arxiv.org/abs/2006.05639) separates a General Search Unit (GSU) from an Exact Search Unit (ESU). The GSU applies hard or soft search to produce a shorter Search-based Interest Sequence, and the ESU performs target-aware attention over that sequence. This pattern has had lasting engineering value: most events are removed cheaply, and expensive interactions are concentrated on a small set.

Its boundary is equally clear. Hard rules are easy to index, yet can miss semantically related events from different categories. Soft retrieval is more flexible, while a coarse representation or objective trained separately from the final task can still create stage mismatch.

### 3.2 ETA: Training Hash Retrieval with the Task

[ETA](https://arxiv.org/abs/2108.04468) maps targets and behaviors to binary codes with SimHash/LSH, filters them by Hamming distance, and trains the retrieval mechanism with the final task. This reduces the gap between handcrafted rules and ranking objectives while making high-dimensional similarity search less expensive online.

Important parameters include hash width, the number of hash tables, bucket capacity, and final $K_r$. Short codes increase collisions; very long codes can reduce recall. These parameters should be evaluated jointly through retrieval recall, the end objective, and tail latency.

### 3.3 TWIN: Preserving Coarse-to-Exact Consistency

[TWIN](https://arxiv.org/abs/2302.02352) shares a target-aware relevance structure between the GSU and ESU. Feature splitting and cached static projections lower coarse-stage cost: stable features are projected offline, and the online path computes only a lightweight candidate-dependent interaction bias. This preserves target awareness and reduces the loss caused by distinct scoring spaces in the two stages.

### 3.4 MUSE: Improving Retrieval with Multimodal Semantics

[MUSE](https://arxiv.org/abs/2512.07216) introduces multimodal item representations into long-sequence retrieval. Semantic and behavioral contrastive learning produce embeddings suited to retrieval, followed by target-aware mechanisms in both coarse filtering and exact modeling. The method illustrates a broader direction: once visual, textual, and behavioral meaning share a searchable space, history selection can cross discrete ID and category boundaries.

This approach depends on a stable multimodal representation service, explicit version management, and a reliable index-refresh path. A semantic-model update can also create version skew between stored history embeddings and online candidate embeddings, which should be monitored directly.

### 3.5 Applicability and Checks

Retrieval fits settings with strong candidate dependence, high value variation, and a clear online budget, including advertising and tightly related recommendations. Before launch, evaluate at least:

- recall of target-relevant history at Top-$K_r$;
- correlation between coarse and exact scores;
- index versions, incremental-update delay, and fallback behavior;
- whether retrieval time is smaller than the interaction time it removes;
- whether discarded events contain transitions required to explain interest drift.

## 4. Compression: Turning History into a Few Interest Tokens

Compression uses $M$ learned queries or interest anchors to encode a sequence of length $L$:

$$
Z=\operatorname{CrossAttention}(Q,H),\qquad Q\in\mathbb{R}^{M\times d},
$$

where $M\ll L$. The final model interacts only with $Z$. If interests are precomputed or updated less frequently, online cost depends mainly on $M$.

### 4.1 Structural Inspiration from BLIP-2

[BLIP-2](https://arxiv.org/abs/2301.12597) is a vision-language model rather than a recommender. Its reported Q-Former setup uses 32 learned query tokens that cross-attend to features from a frozen image encoder. The relevant structural insight is that a small learned query set can form a representation bottleneck between a long input and a downstream model.

The value 32 does not establish a suitable number of recommender interests. $M$ should be selected from user-interest diversity, history length, target quality, and online cost.

### 4.2 CFormer: Returning to Original Behavior Representations

[CFormer](https://dl.acm.org/doi/10.1145/3746252.3761294) first applies learned queries and a Transformer encoder to derive a small set of interest centers. A decoder or reconstruction objective constrains compression quality during training. It then uses final-layer encoder attention to assign each historical event to one center and pools the original event embeddings in each group. The resulting interests retain information from the original feature space before entering the ranking model.

The design has four steps:

1. learn query embeddings and form interest tokens;
2. use reconstruction to prevent a collapsed compression representation;
3. hard-assign behaviors to interest centers using attention;
4. pool original behavior representations within each group and interact them with the candidate.

An implementation should define interest count $M$, assignment temperature or hardening policy, reconstruction-loss weight, empty-cluster handling, and refresh frequency. If interests are generated offline, the evaluation should measure how stale interests affect short-term preferences.

### 4.3 UxSID: Interest Anchors Coordinated with Semantic IDs

[UxSID](https://arxiv.org/abs/2605.09040) applies Interest-Aware Information Compression to map a long history to $M$ interest anchors, then uses the candidate semantic ID (SID) for hierarchical probing:

- **Global probing** lets the target SID query raw behaviors and capture broad relevance;
- **Local probing** uses the global interest to modulate the SID before querying compressed anchors;
- **Fusion** combines global and local results into a target-aware user representation.

This structure puts reusable long-term interests and candidate semantics on the same path. It also requires a dependable SID generation, storage, versioning, and online-query pipeline.

### 4.4 The Right Role for Compression

Fixed-length interests fit commerce, utility products, and content domains with relatively stable preferences. They also provide a useful memory layer in a larger long-sequence system. A compressed representation is reusable across candidates and can substantially reduce per-candidate cost. The lossy bottleneck can weaken rare behaviors and abrupt short-term changes, so a candidate-specific path should remain available to supply local detail.

## 5. Hierarchy: Dividing Work Between Local Segments and Global State

A hierarchical method partitions history into chunks. Each chunk is modeled locally before summary tokens enter a global layer:

$$
z_j=f_{\text{local}}(h_{(j-1)B+1:jB}),\qquad
u=f_{\text{global}}(z_1,\ldots,z_{\lceil L/B\rceil}).
$$

Here $B$ is the block length. The local stage preserves adjacent detail, while the global stage communicates across blocks. HAT-like structures can be understood through this abstraction. They are natural for continuous viewing sessions and periodic behavior because nearby interactions often share context.

The risk lies in the partition assumption. Semantically related events may be far apart, and a single summary per block compresses a returning interest twice. Block boundaries, overlap width, global-token count, and positional encoding therefore need separate sensitivity studies.

[LONGER](https://arxiv.org/abs/2505.04421) combines global tokens, token merging, lightweight InnerTransformers, and hybrid attention. Its design also addresses mixed precision, activation recomputation, KV caching, and integrated training and serving systems. The method shows that practical hierarchical capacity depends on algorithms and runtime design together. In a hybrid model, hierarchy can perform temporal compression so that a semantic sparse layer focuses on cross-block relations.

## 6. Sparse Interaction: Designing High-Value Attention Edges

Sparse methods preserve the original behavior nodes while evaluating only selected positions of the attention matrix:

$$
A_{ij}=
\begin{cases}
\operatorname{softmax}(q_i^\top k_j), & (i,j)\in\mathcal E,\\
0, & \text{otherwise},
\end{cases}
$$

where $\mathcal E$ can be defined by direction, time windows, semantic neighbors, or domain structure. If each node connects to only $K_s$ neighbors, the central computation can fall from $O(L^2)$ to $O(LK_s)$.

![A hybrid architecture combining semantic sparse interaction with compressed memory](assets/diagrams/long-sequence-recommendation/en/hybrid-memory.svg)

### 6.1 Sparse Edges from Structural Priors

A [public WeChat Channels technical article](https://mp.weixin.qq.com/s/xPidLQfNEF-fCCVksT9u_w) presents one structured approach to long-sequence recommendation. The candidate applies one-way cross-attention to the history, while internal sequence interactions are limited through perceived similarity, perceived relevance, and perceived neighborhood.

- **Perceived similarity** converts high-dimensional behaviors into candidate-similarity and attribute signals;
- **Perceived relevance** clusters history into $K$ interests before the candidate reads the aggregated sequence;
- **Perceived neighborhood** uses a sliding window to retain local context and adjust target attention;
- **Multi-field interaction** uses two-dimensional or cross-mask structures to constrain directions between fields and behaviors.

Such rules encode domain knowledge directly in the interaction graph and remain relatively interpretable. Masks can become complicated as fields and tasks grow, so each edge type should have separate coverage, latency, and incremental-gain measurements.

### 6.2 EST: Top-$K$ Adjacency from Multimodal Semantics

[EST](https://arxiv.org/abs/2602.10811) contributes two central components:

1. **LCA (Lightweight Cross-Attention)** preserves the high-value Non-behavior-to-Behavior interaction direction and removes lower-yield directions;
2. **CSA (Content-aware Sparse Attention)** uses frozen multimodal embeddings to construct semantic Top-$K$ neighbors for every behavior.

The paper uses $K=5$ semantic neighbors, turning the internal sequence from a complete graph into a row-sparse adjacency graph. This is a reported configuration, not a universal online default. Engineering studies should jointly scan $K_s$, semantic-model version, adjacency refresh interval, and the fallback for events without neighbors.

### 6.3 ULTRA-HSTU: Semi-Local Attention and System Optimization

[ULTRA-HSTU](https://arxiv.org/abs/2602.16986) uses semi-local attention. Non-behavior tokens can read the full behavior sequence, while behavior-to-behavior interactions mainly use local windows. The work also includes low-precision execution, specialized kernels, input organization, and distributed execution in its optimization scope.

The practical lesson is straightforward: sparse matrices produce end-to-end gains only when the runtime skips invalid edges. Constructing a dense matrix and multiplying by a mask does not realize the theoretical complexity reduction. An implementation should inspect:

- whether the sparse topology can be blocked or compiled into a regular kernel;
- whether adjacency-index reads become a new bandwidth bottleneck;
- numerical stability and calibration under FP8 or BF16;
- whether training and inference use identical sparse semantics;
- load balancing and batching efficiency under variable sequence lengths.

## 7. Linear Sequence Models: Folding History into Recurrent State

### 7.1 Mamba4Rec

[Mamba4Rec](https://arxiv.org/abs/2403.03900) introduces a selective state-space model to sequential recommendation. A continuous form is

$$
h'(t)=A(t)h(t)+B(t)x(t),\qquad
y(t)=C(t)h(t)+D(t)x(t).
$$

The discrete implementation updates a finite-dimensional state along the sequence, with complexity approximately linear in $L$. Selectivity makes the step size and some state parameters input-dependent, allowing the model to choose what to preserve or forget according to content.

This route is attractive under tight memory budgets, very long inputs, and sequential scans. Ranking often still benefits from explicit target-history interaction, so a lightweight target-aware layer can be added beyond the state-space backbone. Evaluation should also confirm that the framework has an efficient selective-scan kernel; a naive loop can erase much of the theoretical advantage.

### 7.2 TokenMixer-Large

[TokenMixer-Large](https://arxiv.org/abs/2602.06563) alternates token mixing and channel mixing over a fixed set of heterogeneous feature tokens, with mixing-and-reverting, residual paths, and sparse Mixture-of-Experts. It is best understood as an efficient feature-interaction model, rather than a complete replacement for arbitrary-length behavior modeling.

Its engineering challenges include MoE routing, expert load balancing, low-precision training, and token parallelism. Within a long-sequence stack, a suitable position may be the feature-interaction layer after retrieval or compression.

## 8. Choosing a Route

![A decision framework for long-sequence recommendation](assets/diagrams/long-sequence-recommendation/en/decision-framework.svg)

### 8.1 Begin with the Business Signal

| Setting                                                                         | Preferred starting point                  | Reason                                                                                    |
| ------------------------------------------------------------------------------- | ----------------------------------------- | ----------------------------------------------------------------------------------------- |
| Advertising, strongly related candidates, explicit real-time target             | Retrieval                                 | Finds directly relevant history first, prioritizing precision under a controlled budget   |
| Short video or content feeds, extreme length, rapid drift, rich multimodal data | Semantic sparse interaction + compression | Sparse edges retain candidate detail; compressed state supplies stable long-term memory   |
| Commerce or utility products, stable categories and preferences                 | Compression                               | Interest tokens are reusable across candidates and can decouple cost from raw length      |
| Strong session continuity and local temporal structure                          | Hierarchy + local windows                 | Chunks preserve continuous patterns while a global layer carries slow state               |
| Tight memory, offline research, or rapid context-length scaling                 | Linear model                              | Complexity and activations grow gently with length, making a useful long-context baseline |

### 8.2 Account for Team Resources

| Compute | Engineering capacity | Suggested starting point                                                                           |
| ------- | -------------------- | -------------------------------------------------------------------------------------------------- |
| Limited | Strong               | Two-stage retrieval; begin with rules or hashes, then add semantic representations                 |
| Ample   | Strong               | Semantic sparse interaction with indexes, kernels, low precision, and monitoring designed together |
| Limited | Limited              | Fixed-slot compression, preferably decoupled from the main ranker, with reuse measured first       |
| Ample   | Limited              | Offline and small-traffic validation first; postpone custom kernels or a full pipeline redesign    |

### 8.3 A Deployable Hybrid Architecture

A practical hybrid system contains two parallel paths and a fusion layer:

1. **Long-term memory path:** compress all history into $M$ interest tokens, update them at a lower frequency, and reuse them across candidates;
2. **Candidate-detail path:** use multimodal semantics and temporal windows to retrieve $K_r$ raw events or create $K_s$ edges per event;
3. **Fusion:** let the candidate read both interest tokens and the sparse behavior subgraph before the ranking head.

This separation gives long-term preference and immediate relevance distinct capacities and refresh rates. If the semantic index is unavailable, the system can fall back to interest memory plus a recent window. If an interest cache is stale, the candidate-detail path can still read the latest events.

## 9. Important Parameters and Experimental Design

Parameters should be managed in three groups: information capacity, compute budget, and data freshness.

| Parameter              | Meaning                                   | Main effect                                            | Recommended measurements                                              |
| ---------------------- | ----------------------------------------- | ------------------------------------------------------ | --------------------------------------------------------------------- |
| $L_{\max}$             | Maximum history length                    | Coverage, feature I/O, and tail latency                | Quality-versus-length curve; P50/P95/P99 lengths and latency          |
| $K_r$                  | Events retained by retrieval              | Target-event recall and exact-stage cost               | Top-$K$ recall, final quality, retrieval and interaction time         |
| $M$                    | Compressed interest-token count           | Interest capacity and per-candidate cost               | Reconstruction, slot utilization, empty-slot rate, quality saturation |
| $B$                    | Hierarchical block length                 | Local detail and number of global tokens               | Cross-block loss, boundary sensitivity, local-stage time              |
| $W$                    | Local attention window                    | Temporal continuity and sparsity                       | Contribution by time span and window-boundary effects                 |
| $K_s$                  | Semantic neighbors per node               | Semantic coverage and sparse-kernel cost               | Neighbor recall, active-edge ratio, sparse-operator throughput        |
| $d$                    | Sequence representation width             | Capacity, memory, and bandwidth                        | Parameters, activations, throughput, and degradation at small widths  |
| $T_{\text{refresh}}$   | Index or interest refresh period          | Freshness and update cost                              | Staleness, version skew, and incremental-job latency                  |
| $\lambda_{\text{aux}}$ | Reconstruction or contrastive-loss weight | Balance between representation structure and main task | Main objective, auxiliary loss, and gradient ratio                    |

A useful ablation order is:

1. hold model capacity fixed and scan $L_{\max}$ to establish that additional history contains useful signal;
2. build separate retrieval, compression, local-window, and linear-model baselines;
3. compare routes under a fixed latency or FLOPs budget, rather than parameter count alone;
4. scan $K_s$ for semantic sparsity and compare with random, temporal, and category neighbors;
5. scan $M$ for compression and inspect whether interest slots receive balanced use;
6. test interest drift, cold start, missing history, stale semantic indexes, and exceptional sequence lengths;
7. report training throughput, online batch size, mean and tail latency, peak memory, and feature-service cost.

> [!important] Theoretical sparsity still requires end-to-end measurement
> FLOPs, attention-edge count, and parameter count explain only part of the cost. Index lookup, feature retrieval, data movement, kernel launch, padding, batching, and fallback traffic all affect online performance. Every route should be measured on the same hardware, batch regime, and input distribution with complete latency accounting.

## 10. Engineering Review Checklist

### Information and Objectives

- Which truncation, filtering, clustering, and pooling operations affect each behavior type before it reaches the model?
- At which layer does the target participate, and are training and inference aligned?
- Does the representation retain order, repeated consumption, negative feedback, and returning interests?
- Do multimodal embeddings improve retrieval of task-relevant behaviors, or only offline similarity scores?

### Systems and Cost

- Does the sparse operator actually skip invalid edges, or does it still construct a dense matrix?
- Where are interest tokens, hash buckets, and semantic adjacency generated, stored, and refreshed?
- How are extreme histories rate-limited, bucketed, or processed asynchronously?
- Can an online failure fall back to a recent window, a short-sequence model, or cached interests?

### Consistency and Observability

- Do the semantic model, index, interest cache, and ranker carry traceable versions?
- Are neighbor coverage, interest-slot utilization, truncation rate, and index staleness monitored?
- Do offline samples reproduce online truncation, refresh intervals, and fallback behavior?
- After release, are results segmented by history length, activity, interest drift, and cold-start state?

## 11. Research Assessment

Long-sequence recommendation is moving from competition among isolated architectures toward joint algorithm-system design. Three trends deserve continued attention:

1. **Hybrid architectures are likely to become more common.** Compressed memory supplies stable, reusable long-term state; sparse interaction performs candidate-specific detailed reads.
2. **Multimodal semantics are becoming infrastructure.** They can support retrieval, graph construction, clustering, cold start, and semantic IDs; their ceiling depends on representation quality and version consistency.
3. **Runtime design determines the production ceiling.** Sparse attention, low precision, caching, indexing, and dynamic batching must work together. Reducing symbolic complexity in isolation is rarely sufficient.

Current public evidence supports the feasibility of these directions, while leaving open which route wins in a particular domain. Reliable selection still requires domain-specific signals, equal-cost ablations, and end-to-end online measurements.

## Primary References

- [SIM: Search-based User Interest Modeling with Lifelong Sequential Behavior Data for Click-Through Rate Prediction](https://arxiv.org/abs/2006.05639)
- [ETA: End-to-End User Behavior Retrieval in Click-Through Rate Prediction Model](https://arxiv.org/abs/2108.04468)
- [TWIN: TWo-stage Interest Network for Lifelong User Behavior Modeling in CTR Prediction at Kuaishou](https://arxiv.org/abs/2302.02352)
- [MUSE: A Simple Yet Effective Multimodal Search-Based Framework for Lifelong User Interest Modeling](https://arxiv.org/abs/2512.07216)
- [BLIP-2: Bootstrapping Language-Image Pre-training with Frozen Image Encoders and Large Language Models](https://arxiv.org/abs/2301.12597)
- [Transformers are Good Clusterers for Lifelong User Behavior Sequence Modeling (CFormer)](https://dl.acm.org/doi/10.1145/3746252.3761294)
- [UxSID: Semantic-Aware User Interests Modeling for Ultra-Long Sequence](https://arxiv.org/abs/2605.09040)
- [LONGER: Scaling Up Long Sequence Modeling in Industrial Recommenders](https://arxiv.org/abs/2505.04421)
- [EST: Towards Efficient Scaling Laws in Click-Through Rate Prediction via Unified Modeling](https://arxiv.org/abs/2602.10811)
- [Bending the Scaling Law Curve in Large-Scale Recommendation Systems (ULTRA-HSTU)](https://arxiv.org/abs/2602.16986)
- [Mamba4Rec: Towards Efficient Sequential Recommendation with Selective State Space Models](https://arxiv.org/abs/2403.03900)
- [TokenMixer-Large: Scaling Up Large Ranking Models in Industrial Recommenders](https://arxiv.org/abs/2602.06563)
