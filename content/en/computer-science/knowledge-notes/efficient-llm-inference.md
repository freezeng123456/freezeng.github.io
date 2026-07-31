---
title: "Efficient LLM Inference: From Token Scheduling to GPU Kernels"
description: Continuous batching, PagedAttention, Transformer tensor flow, sampling, and preemption through the lens of vLLM
lang: en
translation: 计算机科学/知识整理/efficient-llm-inference
tags:
  - LLM-Inference
  - vLLM
  - GPU
  - Systems
---

> [!note] Source and scope
> This page is an independently structured note based on the KM article “[AI Infra Introduction: How Large Models Perform Efficient Inference](https://km.woa.com/articles/show/659449?kmref=profile_feeds)” (updated May 12, 2026), not a reproduction of the original. Example parameters use Llama 3 8B. vLLM behavior follows the v1 implementation analyzed by the article and should be checked against source code after version upgrades. The public version omits identities, organizations, comments, and internal operational information.

## One Picture: The Serving Loop

Efficient inference is not a single Transformer forward pass. It is a continuously running resource-scheduling loop. At each iteration, the system selects the tokens that need computation, executes the model, samples new tokens, and returns unfinished requests to the next iteration.

![The token-level LLM serving loop](assets/diagrams/ml-inference/en/serving-loop.svg)

Here:

- $T$ is the total number of tokens scheduled across all requests in the current iteration;
- $R$ is the number of requests requiring sampling in the current iteration;
- $H$ is the hidden dimension and $V$ is the vocabulary size;
- the Transformer backbone computes over flattened tokens, while Attention reconstructs request isolation from sequence boundaries.

## Continuous Batching: Schedule Work, Not Whole Requests

Static batching requires requests in a batch to begin and end together. Because autoregressive output lengths are unpredictable, completed short requests leave idle slots. Continuous batching moves scheduling granularity down to tokens: every step can admit new requests, advance running requests, and remove completed requests.

The vLLM scheduler can be understood through two counters:

$$
d_r
=
\texttt{num\_tokens}_r
-
\texttt{num\_computed\_tokens}_r ,
$$

where $d_r$ is the number of tokens in request $r$ that have not yet been computed. Each scheduling iteration attempts to make `num_computed_tokens` catch up with `num_tokens` while satisfying:

| Constraint             | What it controls                                     |
| ---------------------- | ---------------------------------------------------- |
| `max_num_seqs`         | Maximum number of concurrently running requests      |
| scheduled-token budget | Total computation per step and kernel shape          |
| `max_model_len`        | Maximum context accepted for one request             |
| free KV-cache blocks   | Whether the current tokens have writable cache space |

All selected valid tokens are concatenated into `[T]` rather than padded to `[batch, max_seq_len]`. This removes padding and allows requests to share the cost of reading model weights.

## PagedAttention: Make KV Memory Allocatable

Preallocating an entire KV-cache region per request creates substantial external fragmentation because output lengths are unknown. PagedAttention divides the KV cache into fixed-size physical blocks and stores the mapping from logical request blocks to GPU physical blocks in a `block_table`:

$$
\text{KV cache shape}
=
[L,2,N_b,B,H_{kv},D_h].
$$

Here $L$ is the number of layers, 2 denotes K/V, $N_b$ is the number of physical blocks, $B$ is block size, $H_{kv}$ is the number of KV heads, and $D_h$ is head dimension.

![PagedAttention mapping from a logical token to a physical KV slot](assets/diagrams/ml-inference/en/paged-kv-cache.svg)

`slot_mapping` tells the kernel where to write new K/V values, while `block_table` tells the Attention kernel which physical blocks contain the historical K/V values. Neither needs a layer dimension: the same token uses the same physical block number and offset at every layer, while accessing a different cache slice for each layer.

This design removes external fragmentation and makes prefix caching, chunked prefill, and dynamic batching more natural. The tradeoff is that block traversal introduces indirect addressing and reduces memory-access continuity. A common `block_size = 16` balances both sides: access remains coalesced within a block, while allocation granularity remains sufficiently fine.

## Tensor-Shape Ledger Through Llama 3 8B

Let $T=\texttt{num\_sched\_tokens}$. The following parameters are used only to develop scale intuition:

| Symbol         | Llama 3 8B example | Meaning                       |
| -------------- | -----------------: | ----------------------------- |
| $H$            |               4096 | hidden size                   |
| $H_q / H_{kv}$ |             32 / 8 | query heads / KV heads in GQA |
| $D_h$          |                128 | head dimension                |
| $I$            |              14336 | FFN intermediate size         |
| $V$            |             128256 | vocabulary size               |
| $L$            |                 32 | Transformer layers            |

The tensor flow through one Transformer block is:

| Stage                      | Input         | Output                          | Systems interpretation                                                     |
| -------------------------- | ------------- | ------------------------------- | -------------------------------------------------------------------------- |
| RMSNorm                    | `[T, 4096]`   | `[T, 4096]`                     | Per-channel scaling, fusible with residual addition                        |
| Fused QKV                  | `[T, 4096]`   | `[T, 6144]`                     | Combines three Q/K/V projections into one wide GEMM                        |
| Split + RoPE               | `[T, 6144]`   | Q `[T,32,128]`; K/V `[T,8,128]` | Q/K use precomputed sin/cos tables; K/V are written through `slot_mapping` |
| Attention                  | Q + paged K/V | `[T,32,128]`                    | Uses `cu_seqlens` to isolate variable-length requests                      |
| Concatenate + O projection | `[T,32,128]`  | `[T,4096]`                      | Output projection after a view transformation                              |
| Fused gate/up              | `[T,4096]`    | `[T,28672]`                     | Combines two FFN projections into one GEMM                                 |
| SiLU-and-multiply          | `[T,28672]`   | `[T,14336]`                     | Fused element-wise kernel                                                  |
| Down projection            | `[T,14336]`   | `[T,4096]`                      | Returns to the residual stream                                             |

The backbone always retains a flattened token view, but Attention must not mix KV values from different requests. `cu_seqlens` supplies request boundaries through cumulative sequence lengths, enabling a kernel to process variable-length sequences without constructing a padded three-dimensional batch.

## Why Prefill and Decode Have Different Bottlenecks

| Phase   |  Query length | Dominant operation                                          | Typical bottleneck |
| ------- | ------------: | ----------------------------------------------------------- | ------------------ |
| Prefill | prompt length | Large GEMMs over multiple tokens                            | compute-bound      |
| Decode  |     usually 1 | GEMV-like projections; reads weights and a growing KV cache | memory-bound       |

Prefill reuses weights across tokens within a request, while continuous batching also lets multiple requests share weights in one step. Decode usually adds only one token per request, but must still read extensive weights and historical KV values; memory bandwidth therefore tends to become limiting before FLOPs.

FlashAttention addresses a separate I/O problem. Through tiling, kernel fusion, and online softmax, it keeps intermediate QK, softmax, and value-multiplication results in on-chip storage instead of repeatedly writing the full attention matrix to HBM.

> [!important] Do not conflate the three optimizations
> Continuous batching primarily amortizes model-weight reads and raises concurrency. PagedAttention manages dynamic KV memory. FlashAttention reduces HBM traffic for Attention intermediates. They address different bottlenecks, but jointly determine end-to-end throughput.

## LM Head and Sampling

In ordinary generation, the LM head need not produce vocabulary logits for every prefill position. The system first gathers the hidden state used by each request to predict its next token:

$$
[T,H]\longrightarrow [R,H]\longrightarrow [R,V].
$$

The sampling pipeline then proceeds as follows:

![The sampling pipeline from logits to the next token](assets/diagrams/ml-inference/en/sampling-pipeline.svg)

Each parameter should be understood by the layer it changes:

| Parameter           | Effect                                                                                |
| ------------------- | ------------------------------------------------------------------------------------- |
| `temperature`       | Controls distribution sharpness; approaches greedy behavior near 0                    |
| `top_k`             | Retains only the $k$ highest-scoring candidates                                       |
| `top_p`             | Retains the smallest candidate set whose cumulative probability reaches the threshold |
| `min_p`             | Adaptively prunes relative to the highest probability                                 |
| grammar / allowlist | Restricts the legal token set at each step                                            |
| repetition penalty  | Reduces the attraction of tokens already present in the output                        |

In the implementation described by the source article, raw logprobs requested by the caller should be snapshotted before penalties, temperature scaling, and truncation; otherwise the returned values incorporate the sampling policy.

## KV Pressure and Preemption

When a RUNNING request can no longer acquire a KV block, vLLM v1 preempts a low-priority running request: it releases that request's KV cache, resets `num_computed_tokens` to zero, returns the request to WAITING, and later reconstructs its state through prefill.

![Request states and preemption under KV pressure](assets/diagrams/ml-inference/en/preemption-states.svg)

Preemption indicates that cache capacity is already tight, so that step generally admits no new WAITING requests; admitting them could trigger repeated preemption. The v1 preference for recomputation over swapping KV to the host reflects a systems tradeoff: PCIe transfers, host-memory use, and a more complex state machine can cost more than repeating prefill.

## A Compact Performance Model

Three categories of traffic provide a quick way to identify the likely optimization direction:

$$
\text{step time}
\approx
\max\left(
\frac{\text{FLOPs}}{\text{GPU compute}},
\frac{\text{model bytes}+\text{KV bytes}+\text{intermediate bytes}}{\text{HBM bandwidth}}
\right).
$$

- Slow prefill: inspect GEMM shapes, Tensor Core utilization, chunked prefill, and parallel partitioning.
- Slow decode: inspect batch token count, weight/KV reads, context length, and KV layout.
- OOM or a batch-size ceiling: inspect the number of KV blocks, block size, maximum context, and scheduling budget.
- High throughput but poor per-request latency: inspect the continuous-batching token budget and fairness.
- Frequent preemption: KV capacity, concurrency, and the context-length distribution are mismatched; merely raising the batch limit will not solve it.

## What to Remember

1. The unifying vLLM view is not “prefill followed by decode,” but that each request's `num_computed_tokens` catches up with its `num_tokens`.
2. Flattened tokens remove padding; request boundaries are reconstructed inside the Attention kernel from variable-length sequence metadata.
3. PagedAttention exchanges address mapping for fine-grained KV allocation. Its principal benefit is capacity utilization rather than faster individual memory accesses.
4. Fused QKV, fused gate/up, fused residual-norm, and FlashAttention all reduce kernel launches or HBM round trips.
5. Prefill is more often compute-bound and decode more often bandwidth-bound. Optimization should begin only after identifying the phase and its actual bottleneck.

## Primary References

- [Efficient Memory Management for Large Language Model Serving with PagedAttention](https://arxiv.org/abs/2309.06180)
- [Orca: A Distributed Serving System for Transformer-Based Generative Models](https://www.usenix.org/conference/osdi22/presentation/yu)
- [FlashAttention: Fast and Memory-Efficient Exact Attention with IO-Awareness](https://arxiv.org/abs/2205.14135)
- [SARATHI: Efficient LLM Inference by Piggybacking Decodes with Chunked Prefills](https://arxiv.org/abs/2308.16369)
