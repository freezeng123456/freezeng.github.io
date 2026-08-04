---
title: GFlowGR 论文复现
description: 用 GFlowNet 把生成式推荐从单点 SFT 推进到价值感知的 set-wise 微调
lang: zh
translation: en/recommender-systems/gflowgr
tags:
  - 推荐系统
  - 生成式推荐
  - GFlowNet
  - 论文复现
---

> [!note] 来源与公开范围
> 本专题基于公开论文 [GFlowGR: Fine-tuning Generative Recommendation Frameworks with Generative Flow Networks](https://arxiv.org/abs/2506.16114)（SIGIR ’26）整理，并结合官方仓库 [Applied-Machine-Learning-Lab/SIGIR26_GFlowGR](https://github.com/Applied-Machine-Learning-Lab/SIGIR26_GFlowGR) 的发布状态说明复现边界。文中只保留可公开的方法结构、关键公式、实验设定与结果解读；不包含内部业务细节、凭据或未授权数据。

## 先看结论

GFlowGR 要解决的不是“再换一个 tokenizer”，而是 **GR 微调目标与线上服务目标不对齐**：

1. SFT 只拟合单个 ground-truth item，线上却要产出一组高价值候选；
2. 不同交互（曝光、点击、购买）效用不同，SFT 却把它们当成等权正样本；
3. DPO / GRPO 一类 reward-based 微调能引入集合或偏好，但通常只给 item 级奖励，缺少 token 级监督。

GFlowGR 的做法是把 item identifier 的自回归生成看成一条 GFlowNet 轨迹，并强制

$$
P(\tau)\propto R(s_L).
$$

这样，高价值候选自然获得更高生成概率，同时每个 token 转移都获得流平衡带来的细粒度梯度。

![GFlowGR 三组件框架](assets/diagrams/gflowgr/zh/framework.svg)

## 页面地图

1. [[推荐系统理论/GFlowGR/方法框架\|方法框架]]：问题建模、轨迹定义，以及采样器 / 奖励 / 损失如何拼成一条训练链路。
2. [[推荐系统理论/GFlowGR/实验复现\|实验复现]]：数据集、基线、超参、训练步骤和验收门槛。
3. [[推荐系统理论/GFlowGR/结果与分析\|结果与分析]]：主结果、消融、线上部署结论，以及和 SFT / RLFT 的差异。

## 复现边界

| 项目                   | 状态                    | 对复现的含义                                 |
| ---------------------- | ----------------------- | -------------------------------------------- |
| 论文方法与公式         | 已公开                  | 可完整重建训练目标与算法流程                 |
| 公开实验设定           | 已公开                  | Beauty / Instruments / Yelp + TIGER / LETTER |
| 官方代码               | README 标明仍在准备发布 | 当前应按论文与已有 GR 开源骨架自建实现       |
| 淘宝生产数据与线上配置 | 未公开                  | 只能引用论文给出的相对指标，不可复刻生产流量 |

## 和现有栏目的关系

- 与 [[推荐系统理论/长序列建模\|长序列建模]] 不同：GFlowGR 关注的是 **生成式推荐的微调目标**，不是长历史编码。
- 与 [[推荐系统理论/OneRec/index\|OneRec]] 相邻：两者都讨论 set-wise / reward-aware 训练；OneRec 更偏工业闭环，GFlowGR 更偏 GFlowNet 目标如何接到 GR backbone。
