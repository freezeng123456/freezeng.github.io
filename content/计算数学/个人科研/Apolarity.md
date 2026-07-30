---
title: Apolarity-guided Taylor Jets
description: 用 Waring 分解计算固定高阶混合偏导
tags:
  - 计算数学
  - 自动微分
  - 代数几何
---

项目：[freezeng123456/apolarity](https://github.com/freezeng123456/apolarity)

## 问题

对一个固定的 $p$ 阶混合偏导，嵌套自动微分会形成越来越深的导数图。Taylor-mode AD 可以沿一个方向传播截断 jet，却没有回答：**究竟需要评估哪些方向，才能用最少次数精确恢复目标混合偏导？**

## 核心洞见

方向 Taylor 系数

$$
T_p(x;v)=\frac1{p!}\frac{d^p}{d\tau^p}u(x+\tau v)\bigg|_{\tau=0}
$$

是方向 $v$ 的 $p$ 次齐次多项式。用若干方向的 $T_p(x;v_r)$ 线性组合提取指定混合偏导，等价于把对应单项式 $z^\alpha$ 写成线性形式的 $p$ 次幂之和，也就是 Waring 分解。

若活跃指数排序为 $a_0\le\cdots\le a_n$，复 Waring 秩给出最短方向数

$$
R_{\mathbb C}(z^\alpha)=\prod_{j=1}^{n}(a_j+1).
$$

根单位构造不仅证明这个上界，而且直接给出可执行的复方向调度。

## 实现骨架

```text
multi-index
  → active exponents
  → roots-of-unity Waring schedule
  → directional Taylor jets
  → weighted exact derivative
  → custom reverse rule
  → parameter-gradient training
```

前向通过 Linear 与 $\sinh$ 的 Taylor 递推传播到 $p$ 阶，反向用自定义规则保持对网络参数可微。

## 核心结果

- 导数值和参数梯度均与嵌套自动微分对齐。
- 对含重复索引的混合偏导，复 Waring 方向数可少于实 polarization；对 square-free 模式，两者均需 $2^{p-1}$ 个方向，没有方向数优势。
- 在共同宽度、共同 1200 秒预算、五个随机种子的 12 个 PDE 设置中，complex-$\sinh$ 模型的平均相对 $L^2$ 误差均最低，覆盖 polyharmonic、radial chirp 与 lossy Maxwell。

![[assets/research/apolarity-chirp.png]]

## 结论边界

这是针对**一个固定单项式偏导**的精确后端，不是对 Laplacian、$\Delta^m$ 等算子和的通用加速。代数上更短的方向表也不自动等于 wall-clock 更快：复数算术、批处理、精度和显存布局都会改变实际收益。

PDE 结果支持 complex-$\sinh$ 表示的潜力，但不能把全部精度提升归因于 Waring 方向缩减；参数量匹配与运行时/峰值显存对照仍是必要的后续实验。
