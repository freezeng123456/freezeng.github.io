---
title: 语义 ID 与 RQ-VAE
description: 连续广告表征如何变成可生成、可检索的多级 token
lang: zh
translation: en/recommender-systems/onerec/semantic-id-and-rq-vae
tags:
  - OneRec
  - RQ-VAE
  - 语义ID
  - 内部评审
---

## 为什么需要语义 ID

直接对海量 TID/AID 做 softmax，词表大、冷启动差，也难把相似广告共享统计强度。OneRec 先把连续表征量化为多级离散 token：

$$
\text{item embedding}
\rightarrow(c_1,c_2,\ldots,c_L)
\rightarrow \text{SID}.
$$

相似内容可能共享较浅层 token，深层 token 再逐步区分细节。生成模型由此在一棵隐式语义树上搜索，而不是一次性枚举全部 item。

## RQ-VAE 结构

当前代码的抽象结构是：

```text
input embedding
  → MLP encoder
  → residual vector quantizer
  → code indices c₁...cL
  → MLP decoder
  → reconstructed embedding
```

残差量化逐层工作。令 $r_0=z_e$，第 $l$ 层选择最近 code $e_{c_l}^{(l)}$，并更新

$$
r_l=r_{l-1}-e_{c_l}^{(l)}.
$$

多层 code 之和近似 encoder 输出。训练损失为

$$
\mathcal L
=\mathcal L_{\text{recon}}
+\lambda_q\mathcal L_{\text{quant}},
$$

量化项包含 codebook loss 与 $\beta$ 加权 commitment loss。

## 当前基线参数

| 参数                 | 当前默认/基线 | 作用                     |
| -------------------- | ------------: | ------------------------ |
| input dimension      |          1024 | 输入广告/内容向量维度    |
| codebooks            |             3 | SID 层数                 |
| entries per codebook |   256/256/256 | 每层词表规模             |
| latent dimension     |           128 | 量化空间维度             |
| encoder hidden       |  1024→512→128 | 连续表征压缩             |
| batch size           |          1024 | 训练吞吐与 code 使用统计 |
| learning rate        |     $10^{-3}$ | AdamW 初始学习率         |
| weight decay         |     $10^{-4}$ | 正则                     |
| training epochs      |         20000 | 上限，需结合验证选择     |
| reconstruction loss  |           MSE | 连续表征重建             |

仓库中还存在面向实验的较小 YAML 结构，它不是生产默认。文档和发布流程必须记录实际生效配置，不能只看文件名推断模型。

## 选择 checkpoint 的指标

仅看 reconstruction loss 不够。训练器同时监控：

- reconstruction loss；
- collision rate：不同 item 是否过多落到同一完整 SID；
- 每层 codebook usage；
- 空码/死码比例；
- token 频率长尾。

当前选择逻辑偏重 collision rate。原因是在线召回最终需要 SID→TID；若大量不相似 item 冲突，即使向量重建误差小，也会导致 KV 候选膨胀和语义混杂。

## SID 打包

在线响应把每级 token 按固定位宽打包进 `sid_key`。当前实现以每级 12 bit 的格式处理，因此必须验证：

$$
0\le c_l < 2^{12}
$$

且 Encoder/Decoder、KV 构建和 Creative Lookup 使用相同层数与顺序。codebook 当前只有 256 个 entry，不代表协议可以随意改成 8 bit；12 bit 是服务接口与未来容量边界的一部分。

## 版本不变量

一个可发布模型包至少要绑定：

```text
rqvae_version
codebook_hash
tree_depth
bits_per_level
sequence_model_version
kv_creative_version
sid_pack_version
```

最危险的故障是“所有 RPC 都成功，但版本错位导致 SID 大量 miss”。它不会表现为服务错误率上升，只会表现为 SID→TID 漏斗断崖。因此版本维度必须进入监控标签。

## 质量与系统成本的权衡

- codebook 更大：碰撞可能下降，但 Decoder 每层分支和图档位成本上升；
- 层数更深：表达容量增加，但 beam 搜索更长、误差逐层累积；
- latent 更大：重建能力提高，但 RQ-VAE 与序列模型输入成本增加；
- collision 更低不必然带来线上收益：过细 token 可能削弱语义共享与泛化。

最终应联合评估量化指标、生成模型 token 指标、SID→TID 覆盖率、召回质量和线上时延。
