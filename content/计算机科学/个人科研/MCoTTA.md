---
title: MCoTTA：动量引导的持续测试时适应
description: 用梯度动量估计更稳定、更前瞻的 PCA 更新子空间
tags:
  - 计算机科学
  - Test-Time-Adaptation
  - AAAI
---

项目：[freezeng123456/MCoTTA](https://github.com/freezeng123456/MCoTTA)

## 问题

LCoTTA 在 continual test-time adaptation 中收集熵梯度快照，使用滑动窗口 PCA 得到低秩子空间，再把当前梯度投影到这个子空间中更新模型。问题在于：单批次原始梯度含有较强瞬时噪声，PCA 队列可能持续追逐短期方向，而不是追踪污染流中的慢变化结构。

## 最小改动

原方法把当前过滤梯度 $g_t$ 直接写入 PCA 队列：

$$
\mathcal Q_t\leftarrow g_t.
$$

MCoTTA 先维护指数动量

$$
m_t=\beta m_{t-1}+(1-\beta)g_t,
\qquad
\mathcal Q_t\leftarrow m_t.
$$

PCA 得到投影矩阵 $P_t$ 后，真正用于更新的仍是当前原始梯度：

$$
\tilde g_t=P_t g_t.
$$

因此方法没有引入 teacher–student、memory bank、routing network 或额外适应目标；它只改变**用什么观测来估计子空间**。

## 为什么可能有效

EMA 相当于一阶低通滤波。若梯度可分为缓慢变化信号 $s_t$ 与零均值噪声 $\epsilon_t$，则

$$
g_t=s_t+\epsilon_t,
$$

动量可以降低队列中噪声方向的方差，让主成分更接近跨批次持久方向。更重要的机制问题是：这个子空间是否能解释**未来**到达的梯度，而不仅是刚被放入队列的当前梯度？

冻结模型的 shared-replay 诊断给出了直接证据。窗口 $W=100$、秩 $r=25$、采样间隔 $K=50$ 时，原始梯度 PCA 对未来 1–49 个 batch 的平均能量保留约 1.33%；$\beta=0.9$ 与 $0.99$ 的动量 PCA 分别约 3.31% 与 3.69%。

![[assets/research/mcotta-retention.png]]

## 重要参数

| 参数       | 含义           | 当前理解                                        |
| ---------- | -------------- | ----------------------------------------------- |
| $\beta$    | EMA 动量系数   | 越大越平滑，但过大可能滞后；主实验使用 0.99     |
| $W$        | PCA 快照窗口   | 大窗口更稳、内存与分解成本更高；主实验为 100    |
| $r$        | 子空间秩       | 控制适应自由度；ResNet 主设置 25，ViT 主设置 50 |
| $K$        | 快照间隔       | 控制子空间更新频率；ResNet 为 50，ViT 为 100    |
| batch size | 单次梯度样本量 | 主实验 64                                       |

ImageNet-C/ResNet-50 的 rank 消融中，$r=10,50,75$ 的准确率约为 43.01%、42.26%、40.31%；更高秩没有自动带来更好结果，说明“保留更多方向”也会重新引入噪声和不稳定更新。

## 核心结果

主协议采用 severity 5、15 种 corruption 连续输入、边界不 reset、batch size 64。

- ImageNet-C / ResNet-50：MCoTTA 三个已完成种子的平均准确率约 43.97%。
- CIFAR-100-C / ResNeXt-29：MCoTTA 三个种子的平均准确率约 66.00%。
- ImageNet-C / ViT-B/16：两个已完成种子的暂定平均约 63.23%，第三个种子尚未纳入。
- 已完成的 seed 2、3 配对结果中，ResNet-50 上 MCoTTA 均高于对应 LCoTTA；seed 1 缺少配对的 LCoTTA 日志，因此不混入配对结论。

## 科学边界

- 机制结论是“动量子空间保留更多未来梯度能量”，不是证明所有未来梯度都对适应有益。
- ViT 三种子结果尚未完整，不能把暂定均值当最终结论。
- 长时间 10-cycle 实验中两种方法后期都有退化；MCoTTA 在每个已完成 cycle 仍较高，但并未解决无限期持续适应的漂移。
- 原始 Raw-PCA 在“先写入当前梯度，再投影当前梯度”时具有 in-sample 自包含效应。严格因果对照应先投影、后写入快照；相关结论必须与现有实验分开陈述。
- 页面不公开论文的全部技术细节、完整消融与实现技巧。
