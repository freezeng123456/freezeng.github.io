---
title: Computer Science Knowledge Notes
description: Computer science notes from algorithmic abstractions to systems implementation
lang: en
translation: 计算机科学/知识整理
tags:
  - Computer-Science
draft: false
---

This section organizes computer science knowledge around problems, mechanisms, data flow, and performance boundaries instead of accumulating disconnected links.

## AI Systems

### [[en/computer-science/knowledge-notes/efficient-llm-inference|Efficient LLM Inference: From Token Scheduling to GPU Kernels]]

Using vLLM as the main thread, four flowcharts and a tensor-shape ledger explain:

- how continuous batching moves scheduling granularity from requests to tokens;
- how PagedAttention manages a dynamic KV cache with block tables;
- how data for Llama 3 8B flows through Attention, the FFN, the LM head, and the sampler;
- which hardware and resource constraints dominate prefill, decode, and preemption.

## Organization

- **Foundations**: optimization, generalization, representation learning, and online learning;
- **Systems**: training systems, inference serving, data pipelines, and observability;
- **Experimental methods**: controls, ablations, statistical uncertainty, and reproducible protocols.
