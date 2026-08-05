---
title: 稀疏恢复与数据驱动混沌
description: 编号 10、21、29、32、36、44：把梯度信息与经验测度写进恢复问题
lang: zh
translation: en/computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce
tags:
  - 论文笔记
  - 不确定性量化
  - 稀疏恢复
---

> [!note] 本页覆盖
> 编号 **10**（_SIAM J. Sci. Comput._ 36(4), 2014）、**21**（_SIAM J. Sci. Comput._ 39(1), 2017）、**29**（_Commun. Comput. Phys._ 24, 2018）、**32**（_J. Comput. Phys._ 367, 2018）、**36**（_J. Comput. Phys._ 381, 2019）、**44**（_Commun. Math. Res._ 36, 2020）。

![把梯度信息与采样密度一起写进恢复问题](assets/diagrams/tao-zhou-papers/zh/sparse-recovery.svg)

## 稀疏恢复与最小二乘共享同一个决定量

在稀疏恢复框架下，设计矩阵被要求满足有界正交系下的限制等距性质，而所需测量数由行范数的均匀性控制——**与[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|最小二乘情形]]完全相同的量**。这一点解释了为什么本页的工作与那一页共享 Christoffel 权与预条件的语言。

## 10：确定性插值点与稀疏插值

编号 10 研究稀疏插值与确定性插值点的设计。稀疏插值的问题是：若已知目标在某个基下稀疏，能否用远少于空间维数的点把它恢复出来，且点集是确定性给出的（而非随机抽取）？

论文用的技术工具与[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 9]] 同源（Weil 和机制），因此其结论同样是确定性而非概率性的。编号 **29（稀疏三角多项式的梯度增强 ℓ1 恢复）** 用同一套机制处理三角多项式基。

## 21：对 Gauss 网格做随机子采样后做 ℓ1 极小化

编号 21 是[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 13]] 的 ℓ1 兄弟篇：同样的随机子采样 Gauss 求积设计，那里是最小二乘，这里是稀疏恢复。它把 Tang 与 Iaccarino 只覆盖均匀分布的结论推广到 Beta、正态与指数参数。

它的中心观察与编号 13 相同：**Gauss 权恰好是 Christoffel 函数值**，因此对 Gauss 网格做子采样隐含地就是 Christoffel 加权采样。相应地，张量 Gauss 网格 $\Theta_{\mathbf n}$ 上的权为

$$
w_{\mathbf k}=\lambda_{\mathbf n}(z_{\mathbf k})
=\prod_{i=1}^{d}\frac{1}{\sum_{k=0}^{n_i-1}\bigl[\phi^i_k(z^i_{k_i})\bigr]^2},
$$

而网格上的均匀经验概率测度

$$
\nu_{\mathbf n}=\frac{1}{\prod_i n_i}\sum_{\mathbf k\le\mathbf n}\delta_{z_{\mathbf k}}
$$

使「从 $\nu_{\mathbf n}$ 独立同分布采样」等价于「从张量 Gauss 网格均匀采样」。

论文进一步指出，设计矩阵按 Christoffel 函数加权后，得到的函数族 $\{\psi_{k,\mathbf n}\}$ 在离散测度 $\nu_{\mathbf n}$ 下**正交归一**，相应的矩阵 $D=\bigotimes_iD^i$（$D^i=(\Sigma^i)^{1/2}\Psi^i$）是正交矩阵。这正是「行归一化」在张量结构下的具体实现。

论文同时量化了无界域的困难（论文给出的量级形如 $n^{2d/3}$），这是[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 11]] 在最小二乘情形遇到的同一困难在稀疏恢复框架下的表现。

## 29 与 32：梯度增强的 ℓ1 恢复

### 想法

一次昂贵的正问题求解通常同时给出函数值与**全部偏导数**（例如通过伴随方法或自动微分）。若只用函数值，每个样本贡献设计矩阵的 1 行；若同时用梯度，同一个样本贡献 $1+d$ 行，而未知量个数不变。

因此测量矩阵由

$$
A_{\text{value}}\in\mathbb R^{M\times N}
\qquad\text{变为}\qquad
\begin{bmatrix}A_{\text{value}}\\ A_{\text{grad}}\end{bmatrix}\in\mathbb R^{M(1+d)\times N},
$$

对应的恢复问题为

$$
\min_c\ \|c\|_1
\qquad\text{s.t.}\qquad
\bigl\|W(Ac-b)\bigr\|_2\le\delta,
$$

其中 $W$ 是把各行范数拉平的对角权矩阵。**这就把「每个昂贵样本的信息量」变成一个可以设计的对象**：不是增加样本数，而是从同一个样本里取出更多行。

### 两个必须处理的问题

第一，**行范数不再一致**。值行与梯度行的量级不同（梯度行带一个导数的尺度），因此预条件不再是可选项：若不加权，限制等距性质的常数会被最坏的一类行支配。这正是编号 32 与[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 24]] 共享「采样与预条件成对设计」这一取向的原因。

第二，**梯度并非免费**。若梯度需要额外的求解（而不是伴随方法的副产品），则 $(1+d)$ 倍的行数要与相应的代价对比。论文的适用场景正是梯度可廉价获得的场合。

编号 29 处理稀疏**三角**多项式的恢复，编号 32 处理多项式混沌展开的稀疏逼近。二者共享同一个测量堆叠结构，差别在于基与相应的有界正交系常数。

## 36 与 44：输入分布只通过样本给出时

### 问题

标准多项式混沌假设输入分布已知且有解析形式，因此正交基可以由该分布的正交多项式给出。**数据驱动**情形下，输入分布只通过有限样本给出：没有解析密度，因此没有现成的正交基。

编号 36（数据驱动多项式混沌展开：加权最小二乘逼近）与编号 44（数据驱动多项式混沌展开的稀疏逼近：诱导采样路线）处理这一设定。要点是：正交基必须从**经验测度**构造（例如用 Gram-Schmidt 或三项递推的数值版本），而采样设计随之改变——诱导采样与 Christoffel 加权都需要一个基，而基本身现在是估计出来的。

这引出一个本专题内其他工作不涉及的误差来源：**基的估计误差**。它与逼近误差、采样误差叠加，因此数据驱动情形的分析比已知分布情形多一层。

> [!note] 覆盖进度
> 编号 21 的构造（Gauss 权即 Christoffel 值、离散正交归一性、正交矩阵结构）已核对。编号 10、29、32、36、44 的具体定理、样本复杂度常数与实验本站未逐项核实；上述内容限于可从摘要、参考关系与相邻论文交叉引用确认的部分。

## 六篇的关系

| 编号 | 测量的构成             | 采样设计              | 基的来源     |
| ---- | ---------------------- | --------------------- | ------------ |
| 10   | 函数值，确定性点       | Weil 和点集           | 已知分布     |
| 21   | 函数值，Gauss 网格子集 | 隐含 Christoffel 加权 | 已知分布     |
| 29   | 值 + 梯度，三角多项式  | 需预条件              | 已知分布     |
| 32   | 值 + 梯度，多项式混沌  | 需预条件              | 已知分布     |
| 36   | 函数值，加权最小二乘   | 由经验测度决定        | **经验测度** |
| 44   | 函数值，稀疏逼近       | 诱导采样              | **经验测度** |

一条贯穿的判断：**稀疏恢复的样本复杂度由行范数的均匀性决定，因此凡是改变行结构的操作（加梯度行、换基、换采样密度）都必须配一个相应的预条件。** 这条判断把本页与[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|最优采样一页]]统一起来：两者处理的是同一个量，只是一个在最小二乘框架下、一个在 ℓ1 框架下。

## 本页原文

- Z. Xu and T. Zhou, [_On sparse interpolation and the design of deterministic interpolation points_](https://doi.org/10.1137/13094596X), SIAM J. Sci. Comput. 36(4) (2014), pp. A1752-A1769。
- L. Guo, A. Narayan, T. Zhou, and Y. Chen, [_Stochastic collocation methods via ℓ1 minimization using randomized quadratures_](https://doi.org/10.1137/16M1059680), SIAM J. Sci. Comput. 39(1) (2017), pp. A333-A359。
- Z. Xu and T. Zhou, [_A gradient-enhanced ℓ1 approach for the recovery of sparse trigonometric polynomials_](https://doi.org/10.4208/cicp.OA-2018-0006), Commun. Comput. Phys. 24 (2018), pp. 286-308。
- L. Guo, A. Narayan, and T. Zhou, [_A gradient enhanced ℓ1-minimization for sparse approximation of polynomial chaos expansions_](https://doi.org/10.1016/j.jcp.2018.04.026), J. Comput. Phys. 367 (2018), pp. 49-64。
- L. Guo, Y. Liu, and T. Zhou, [_Data-driven polynomial chaos expansions: a weighted least-square approximation_](https://doi.org/10.1016/j.jcp.2018.12.020), J. Comput. Phys. 381 (2019), pp. 129-145。
- L. Guo, A. Narayan, Y. Liu, and T. Zhou, [_Sparse approximation of data-driven polynomial chaos expansions: an induced sampling approach_](https://doi.org/10.4208/cmr.2020-0010), Commun. Math. Res. 36 (2020), pp. 128-153。
