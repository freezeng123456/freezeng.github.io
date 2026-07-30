---
title: RSSE：旋转对称 Stein 估计器
description: 高阶随机平滑 PINN 的带宽稳定估计
tags:
  - 计算数学
  - PINN
  - Monte-Carlo
---

项目：[freezeng123456/RSSE](https://github.com/freezeng123456/RSSE)

## 问题

随机平滑 PINN 用 Stein 恒等式把高阶导数改写为函数值期望，避免显式构造高阶自动微分张量。但当平滑带宽 $t\to0$ 时，Stein 权重含 $t^{-p}$；只要低阶 Taylor 模式没有完全抵消，方差就会被这个奇异因子放大。

经典二点 antithetic 对二阶算子很有效，却没有足够自由度消除三阶及以上目标之前的全部 Taylor 模式。

## 核心洞见

令 $\omega=\exp(2\pi i/N)$。沿同一个高斯方向 $\xi$，在复相位

$$
x+t\omega^j\xi,\qquad j=0,\ldots,N-1
$$

上评估解析网络，并用 $\omega^{-jp}$ 加权。离散旋转对称性相当于一个 Taylor 模式滤波器：目标 $p$ 阶模式被保留，而所有更低阶模式被系统消去。

主结果可概括为

$$
\mathbb E[\mathcal E_{\alpha,t,N}u]
=\partial^\alpha u+O(t^N),\qquad
\operatorname{Var}(\mathcal E_{\alpha,t,N}u)
=V_\alpha(u)+O(t^N),
$$

其中主导方差 $V_\alpha(u)$ 与带宽 $t$ 无关。这是方法区别于普通 antithetic 估计器的关键。

## 算法骨架

1. 固定目标多重指标 $\alpha$，令 $p=|\alpha|$，选择 $N\ge p$；
2. 对每个 Monte Carlo 样本抽取 $\xi\sim\mathcal N(0,I)$；
3. 计算 $N$ 个旋转对称复相位网络值；
4. 用单位根权重提取目标 Taylor 模式；
5. 乘以多元 Hermite 权重并对样本平均；
6. 在 PINN 残差中结合精确边界参数化或 Leibniz 展开。

## 核心结果

在三调和与四调和球问题 $d=20,50,100$ 上，RSSE 相对调优的 antithetic 基线均得到更低的相对 $L^2$ 误差；在六阶非线性 Cahn–Hilliard 算例 $d=5,10,20$ 上也保持优势。

四调和实验的三种维度中，RSSE 的平均相对误差分别约为
$1.10\times10^{-2}$、$7.17\times10^{-3}$、$4.29\times10^{-3}$，对应基线约为
$1.55\times10^{-2}$、$9.83\times10^{-3}$、$7.35\times10^{-3}$。

![[assets/research/rsse-k4.png]]

## 有意保留的边界

- 这里不公开完整权重构造、证明细节与工程优化。
- 理论依赖网络在复邻域中的解析延拓；实现通常选择 entire activation。
- 复数评估会增加单次前向成本，优势来自稳定性与更宽的可用带宽区间，不等于任何硬件上都更快。
- 结果应理解为特定协议下的估计器优势，而不是对所有 PDE、网络和采样预算的普遍结论。

## 与另一项工作的关系

[[计算数学/个人科研/Apolarity\|Apolarity]]追求“精确混合偏导的最短确定性方向调度”，RSSE 追求“随机平滑估计的低阶模式消除”。两者共享单位根、复方向和 Taylor 系数提取的代数结构，但优化目标不同。
