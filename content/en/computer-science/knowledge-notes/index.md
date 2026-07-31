---
title: Machine Learning Knowledge Notes
description: Machine-learning notes from algorithmic abstractions to training, inference, and hardware execution
lang: en
translation: 计算机科学/知识整理
tags:
  - Machine-Learning
draft: false
---

This section organizes machine-learning knowledge around problems, mechanisms, data flow, and performance boundaries instead of accumulating disconnected links.

## AI Systems

### [[en/computer-science/knowledge-notes/efficient-llm-inference|Efficient LLM Inference: From Token Scheduling to GPU Kernels]]

Using vLLM as the main thread, four flowcharts and a tensor-shape ledger explain:

- how continuous batching moves scheduling granularity from requests to tokens;
- how PagedAttention manages a dynamic KV cache with block tables;
- how data for Llama 3 8B flows through Attention, the FFN, the LM head, and the sampler;
- which hardware and resource constraints dominate prefill, decode, and preemption.

### [[en/computer-science/knowledge-notes/how-gpus-work|How GPUs Work: From the CUDA Lifecycle to SIMT]]

This note expands the GPU from an abstract “accelerator” into a traceable execution path:

- how the CPU submits work through the runtime, driver, command queues, and a doorbell;
- how CUDA grids, blocks, and threads map to SMs, warps, and execution units;
- how registers, shared memory, and thread limits jointly constrain occupancy;
- why high occupancy does not guarantee high performance, and how SIMD, SIMT, and warp divergence differ.

## Organization

- **Foundations**: optimization, generalization, representation learning, and online learning;
- **Systems**: training systems, inference serving, data pipelines, and observability;
- **Experimental methods**: controls, ablations, statistical uncertainty, and reproducible protocols.
