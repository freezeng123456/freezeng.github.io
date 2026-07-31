---
title: 旋转对称 Stein 估计器
description: 面向平滑 PINN 的带宽稳定高阶随机微分
lang: zh
translation: en/computational-mathematics/research/rsse
tags:
  - 计算数学
  - PINN
  - Monte-Carlo
---

代码仓库：[freezeng123456/RSSE](https://github.com/freezeng123456/RSSE)

## 问题

在物理信息神经网络中，随机平滑可以把高阶导数替换成由函数值与 Stein 权重构成的期望，从而避免显式构造高阶自动微分张量。然而，对于 $p$ 阶导数，Stein 权重包含量级为 $t^{-p}$ 的因子。如果低阶 Taylor 项没有被精确消去，当平滑带宽 $t$ 趋于零时，这些项会被放大。

标准的两点反对称构造对若干二阶算子有效，但它没有足够的自由度消去一般三阶及更高阶导数之前的所有 Taylor 模态。

## 核心 insight

令 $\omega=\exp(2\pi i/N)$。对一个高斯方向 $\xi$，在复相位

$$
x+t\omega^j\xi,\qquad j=0,\ldots,N-1
$$

处计算解析网络，并以权重 $\omega^{-jp}$ 组合函数值。离散旋转对称性相当于一个 Taylor 模态滤波器：目标 $p$ 阶模态被保留，所有低阶模态被消去。

主要结果具有如下形式：

$$
\mathbb E[\mathcal E_{\alpha,t,N}u]
=\partial^\alpha u+O(t^N),\qquad
\operatorname{Var}(\mathcal E_{\alpha,t,N}u)
=V_\alpha(u)+O(t^N),
$$

其中主导方差 $V_\alpha(u)$ 与 $t$ 无关。因此，它与普通反对称估计器的差别不只体现在更高阶的偏差展开，更在于消除了主导方差对带宽的奇异依赖。

## 算法结构

1. 给定目标多重指标 $\alpha$，令 $p=|\alpha|$，并选择 $N\ge p$。
2. 对每个 Monte Carlo 样本生成 $\xi\sim\mathcal N(0,I)$。
3. 在 $N$ 个旋转对称复相位处计算网络。
4. 施加单位根权重，提取目标 Taylor 模态。
5. 乘以多元 Hermite 权重，并对样本求平均。
6. 将估计器放入 PINN 残差；必要时结合精确边界参数化或 Leibniz 展开。

## 代表性结果

RSSE 在 $d=20,50,100$ 的三调和与双调和问题，以及 $d=5,10,20$ 的六阶非线性 Cahn-Hilliard 问题上进行了评估。在论文给出的实验协议下，它在每组设置中都取得了低于调优反对称基线的相对 $L^2$ 误差。

在双调和实验中，RSSE 在三个维度上的平均相对误差约为
$1.10\times10^{-2}$、$7.17\times10^{-3}$ 和 $4.29\times10^{-3}$；对应基线约为
$1.55\times10^{-2}$、$9.83\times10^{-3}$ 和 $7.35\times10^{-3}$。

![[assets/research/rsse-k4.png]]

## 适用边界

- 本页不公开完整权重构造、证明和实现层面的专门优化。
- 分析要求网络能解析延拓到一个复邻域，因此实现中使用整函数激活。
- 复数计算会提高单次前向传播成本。收益来自方差控制和更宽的可用带宽区间，而不是对墙钟时间的普遍降低。
- 实验比较只对应给定 PDE、网络结构和采样预算，不能据此推断对任意微分算子或神经网络都占优。

## 与 Apolarity 项目的关系

[[computational-mathematics/research/apolarity|Apolarity 引导的 Taylor Jet]]寻找一个精确混合偏导数的最短确定性方向计划；RSSE 则消去随机平滑估计器中的低阶模态。两者都使用单位根、复方向和 Taylor 系数提取，但解决的是不同优化问题。
