---
title: 计算机科学知识整理
description: 从算法抽象到系统实现的计算机科学学习笔记
lang: zh
translation: en/computer-science/knowledge-notes
tags:
  - 计算机科学
draft: false
---

这一部分按“问题—机制—数据流—性能边界”组织计算机科学知识，不堆叠零散链接。

## AI 系统

### [[计算机科学/知识整理/efficient-llm-inference|高效大模型推理：从 Token 调度到 GPU Kernel]]

以 vLLM 为主线，用四张流程图和张量维度表解释：

- continuous batching 如何把调度粒度从请求下沉到 token；
- PagedAttention 如何用 block table 管理动态 KV Cache；
- Llama 3 8B 的数据如何流经 Attention、FFN、LM head 与 sampler；
- prefill、decode、抢占分别受什么硬件与资源约束。

## 内容组织

- **基础理论**：优化、泛化、表示学习、在线学习；
- **系统实现**：训练系统、推理服务、数据管线、可观测性；
- **实验方法**：对照、消融、统计不确定性、可复现协议。
