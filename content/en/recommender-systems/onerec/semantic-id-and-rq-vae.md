---
title: Semantic IDs and RQ-VAE
description: Turning continuous ad representations into generatable and retrievable multilevel tokens
lang: en
translation: 推荐系统理论/onerec/语义id与rq-vae
tags:
  - OneRec
  - RQ-VAE
  - Semantic-ID
  - Internal-Review
---

## Why Semantic IDs?

A direct softmax over a vast TID/AID vocabulary is expensive, handles cold starts poorly, and makes it difficult for similar ads to share statistical strength. OneRec first quantizes a continuous representation into multilevel discrete tokens:

$$
\text{item embedding}
\rightarrow(c_1,c_2,\ldots,c_L)
\rightarrow \text{SID}.
$$

Similar content may share shallow tokens, while deeper tokens progressively distinguish finer details. The generative model consequently searches an implicit semantic tree rather than enumerating every item at once.

## RQ-VAE Structure

The current code has the following abstract structure:

```text
input embedding
  → MLP encoder
  → residual vector quantizer
  → code indices c₁...cL
  → MLP decoder
  → reconstructed embedding
```

Residual quantization proceeds level by level. Let $r_0=z_e$. Level $l$ selects the nearest code $e_{c_l}^{(l)}$ and updates:

$$
r_l=r_{l-1}-e_{c_l}^{(l)}.
$$

The sum of the multilevel codes approximates the encoder output. The training objective is:

$$
\mathcal L
=\mathcal L_{\text{recon}}
+\lambda_q\mathcal L_{\text{quant}},
$$

where the quantization term includes codebook loss and a $\beta$-weighted commitment loss.

## Current Baseline Parameters

| Parameter            | Current default/baseline | Role                                             |
| -------------------- | -----------------------: | ------------------------------------------------ |
| input dimension      |                     1024 | Dimension of the input ad/content vector         |
| codebooks            |                        3 | Number of SID levels                             |
| entries per codebook |              256/256/256 | Vocabulary size at each level                    |
| latent dimension     |                      128 | Quantization-space dimension                     |
| encoder hidden       |             1024→512→128 | Continuous-representation compression            |
| batch size           |                     1024 | Training throughput and code-use statistics      |
| learning rate        |                $10^{-3}$ | Initial AdamW learning rate                      |
| weight decay         |                $10^{-4}$ | Regularization                                   |
| training epochs      |                    20000 | Upper limit; validation must determine selection |
| reconstruction loss  |                      MSE | Continuous-representation reconstruction         |

The repository also contains smaller YAML structures for experiments; these are not production defaults. Documentation and release processes must record the configuration actually applied rather than infer a model from a filename.

## Checkpoint-Selection Metrics

Reconstruction loss alone is insufficient. The trainer also monitors:

- reconstruction loss;
- collision rate: whether too many different items map to the same full SID;
- codebook usage at each level;
- empty/dead-code ratio;
- long-tailed token frequency.

The current selection logic emphasizes collision rate because online retrieval ultimately performs SID→TID mapping. If many dissimilar items collide, a low vector-reconstruction error can still produce inflated KV candidate lists and mixed semantics.

## SID Packing

The online response packs the token at each level into a `sid_key` with a fixed bit width. The current implementation uses 12 bits per level, so it must validate:

$$
0\le c_l < 2^{12}
$$

and the Encoder/Decoder, KV construction, and Creative Lookup must use the same level count and order. A codebook with only 256 entries does not permit the protocol to be changed arbitrarily to 8 bits: 12 bits are part of both the service contract and the future capacity boundary.

## Version Invariants

A releasable model package must bind at least:

```text
rqvae_version
codebook_hash
tree_depth
bits_per_level
sequence_model_version
kv_creative_version
sid_pack_version
```

The most dangerous failure mode is that every RPC succeeds while version misalignment causes widespread SID misses. Service error rates do not rise; instead, the SID→TID funnel collapses. Version dimensions must therefore be included in monitoring labels.

## Quality–System Cost Tradeoff

- Larger codebooks can reduce collisions but raise Decoder branching and graph-bucket cost at each level.
- More levels increase representational capacity but lengthen beam search and allow errors to accumulate.
- A larger latent dimension improves reconstruction capacity but increases RQ-VAE and sequence-model input cost.
- Lower collision does not necessarily imply online gain: overly fine tokens can weaken semantic sharing and generalization.

Evaluation must combine quantization metrics, token metrics from the generative model, SID→TID coverage, retrieval quality, and online latency.
