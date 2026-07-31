---
title: 机器学习
description: 持续测试时适应、机器学习系统与 AI Infra 知识整理
lang: zh
translation: en/computer-science
tags:
  - 机器学习
---

## 个人科研

[[计算机科学/个人科研/MCoTTA\|MCoTTA]]研究 continual test-time adaptation 中如何用动量梯度改善 LCoTTA 的低秩更新子空间。当前页面只介绍最小算法改动、机制证据、主要实验与科学边界。

## 知识整理

[[计算机科学/知识整理/index\|知识整理]]从机器学习系统与 AI Infra 开始：

- [[计算机科学/知识整理/efficient-llm-inference\|高效大模型推理]]以 vLLM 为主线，可视化 continuous batching、PagedAttention、Transformer 张量流、采样与抢占；
- [[计算机科学/知识整理/how-gpus-work\|GPU 如何工作]]解释 CPU–GPU 协作、CUDA 程序生命周期、Grid/Block/Warp/SM 映射、占用率与分支发散。
