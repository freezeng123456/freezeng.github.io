---
title: OneRec Beam Search
description: beam width、逐层保留、top-k、图档位与调优权衡
lang: zh
translation: en/recommender-systems/onerec/beam-search
tags:
  - OneRec
  - Beam-Search
  - 参数
  - 内部评审
---

> [!important] 最容易误解的参数
> `beam width` 不是“最终返回多少广告”，也不是“GPU 图的档位”。它是生成语义 ID 时，每一层最多保留多少个前缀状态。最终 SID 数还会被 top-k、去重、KV quota、SID→TID 扩展和下游过滤继续改变。

## 多级语义树上的搜索

模型分层预测 $L$ 个 token。第 $l$ 层把上一层保留的前缀扩展：

$$
\mathcal C_l
=\{(s_{1:l-1},c):s_{1:l-1}\in\mathcal B_{l-1},c\in V_l\},
$$

为每个候选累加 log-score，再选择最高的 $B_l$ 个：

$$
\mathcal B_l=\operatorname{TopK}_{B_l}(\mathcal C_l).
$$

`UpdateBeamState` 要同时继承父前缀、追加 token、累积分数并维护长度。最终完整 token 序列才被打包为 SID。

## 四个容易混淆的“宽度”

### 1. 请求 beam width

`gpr_pam_beam_width` 是策略/实验传给 Decoder 的主宽度。新 AnyRecall 重构路径中：

- 缺省回退值：180；
- 硬约束：不得超过模型导出的最大宽度 180。

共享 `GenerativeData` 结构中存在 200 的初始值，Entry 随后会用位置实验值或 180 的安全回退覆盖它。这个 200 不能当作旧链路的线上有效配置，也不能据此把新路径改为 200。

### 2. 每层 beam width

Decoder 支持每一层不同宽度 $B_1,\ldots,B_L$。单个 `pam_beam_width` 可能被展开成统一宽度，也可能结合树层配置得到层级数组。

常见设计是浅层保留更多语义簇、深层逐渐收窄，或反过来在首层较窄、深层允许更多细粒度分支。应以具体 callsite 生成的数组为准。

### 3. `beam_search_topk`

这是最终或指定阶段的输出截断。若值 $>0$，最终保留数可覆盖默认 beam 输出；否则通常沿用 beam width。

因此：

```text
beam width = 搜索过程中保留的前缀数
top-k      = 搜索结束对完整 SID 的输出数
```

把 top-k 调小可以减少 KV 查询与下游成本，却不能完全回收前面各层已经发生的 Decoder 计算。

### 4. GPU graph/bucket width

推理引擎常预捕获若干支持宽度的执行图。配置中看到的 `90, 200` 或更大宽度列表，可能是**可复用图档位**，不是当前请求要搜索的宽度。

运行时通常选择能容纳请求的最小档位并 padding。若请求宽度 180 落到 bucket 200，实际 kernel 可能按 200 的 shape 执行；时延与显存应按 bucket 而非逻辑 180 估计。

## 计算与内存成本

若每层词表大小为 $|V_l|$，朴素展开成本近似

$$
O\!\left(\sum_{l=1}^{L} B_{l-1}|V_l|\right).
$$

显存主要随最大活动 beam、每个 beam 的 hidden/cache 和候选 logits 增长：

$$
M\approx O(B_{\max}d_{\text{state}})+O(B_{\max}|V_l|).
$$

实际 kernel 会通过 fused top-k、共享前缀和图捕获降低常数，但宽度增加通常仍近似线性增加 Decoder 算量、KV 查询量与后续候选量。

## beam width 调大带来什么

潜在收益：

- 降低搜索截断误差；
- 增加长尾 SID 覆盖；
- 为 KV 和下游排序提供更多可选候选。

潜在代价：

- Decoder 时延与显存上升；
- 可能切换到更大的执行图档位，产生阶跃式时延；
- SID 数增加，KV fan-out 与网络回包增大；
- 低分 SID 命中大量 TID 后，quota/过滤做了无效工作；
- 排序负担上升，最终广告数未必增加。

因此 beam width 的有效性要看“边际 SID 是否转化成边际有效广告”，而不是只看 Decoder top-k 指标。

## 与分片和 quota 的关系

存在 beam 分片比例参数，用于把总搜索/返回预算分配到并行分片。分片输出在 KV 阶段交织并受 quota 截断。需要保持：

$$
\sum_s B_s \gtrsim B_{\text{request}},
$$

但不能简单认为每个分片返回 $B$ 就能得到 $S\times B$ 个最终 SID，因为去重、空分片和 top-k 会收缩。

## 当前配置解读

| 场景                |       观察到的值 | 正确解释                                      |
| ------------------- | ---------------: | --------------------------------------------- |
| 新生成式重构路径    | fallback/max 180 | Entry 最终解析的请求 beam 与 Decoder 上限     |
| 共享结构初值        |              200 | 会被 Entry 覆盖，不能当作运行值               |
| 某版本样本/实验配置 |       90、200 等 | 必须按参数 callsite 区分 beam、top-k 或图档位 |
| 执行图列表          |     多个离散宽度 | runtime bucket，不代表策略值                  |

## 调优实验

建议同时扫描 $B\in\{30,60,90,120,150,180\}$ 与最终 top-k，记录：

1. Decoder P50/P95/P99 与显存；
2. SID 数、unique SID 与平均分数；
3. SID→TID 命中率和每 SID TID 数；
4. quota 前后 TID/AID；
5. 粗排与精排通过率；
6. 最终曝光、价值与多样性；
7. 各 GPU bucket 的实际命中比例。

选择点应位于“有效广告增益开始饱和、时延尚未跳档”的拐点，而不是机械取最大 180。

## 排障

| 症状                     | 可能原因                                 |
| ------------------------ | ---------------------------------------- |
| 请求直接失败             | beam 超过模型 max、层数/宽度数组不合法   |
| 时延突然阶跃             | 命中更大 graph bucket                    |
| SID 增加但广告不增       | KV miss、重复 SID、AID 去重或 quota 截断 |
| 小 beam 反而效果更好     | 低分长尾引入噪声，挤占下游 quota         |
| 分片结果偏斜             | 分片比例、Z-order 或分片库存分布不均     |
| 新旧路径同名参数效果不同 | 默认值、单位或 callsite 语义不同         |
