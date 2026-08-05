---
title: 时间并行算法
description: 十三项工作把串行时间推进改写成一个可对角化或可预条件的代数问题
lang: zh
translation: en/computational-mathematics/paper-notes/parallel-in-time
tags:
  - 计算数学
  - 论文笔记
  - 时间并行
---

这个专题包含 13 项工作，2015 至 2025 年间与吴淑林等合作完成，并在 2025 年汇总为一篇 _Acta Numerica_ 综述（编号 85）。它们回答同一个问题：**时间推进天然是串行的，如何把这个串行结构改写成一个可以并发求解的代数问题。**

![把串行时间递推换成可对角化的时间矩阵](assets/diagrams/tao-zhou-papers/zh/pint-diagonalization.svg)

> [!info] 与站内综述精读的关系
> 编号 85 的逐节精读是本站另一个专题：[[computational-mathematics/knowledge-notes/time-parallelization/index|双曲与抛物问题的时间并行方法]]，含原论文全部图表与可复现的 Python 实验。本页不重复那份精读，而是说明构成该综述基础的十二项原始工作各自解决了什么问题。

## 两条技术路线

### 迭代路线：粗细传播子

parareal 迭代

$$
\boldsymbol u_{n+1}^{k+1}=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^{k})
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k})
$$

把昂贵的细传播 $\mathcal F$ 并行化，代价是一个仍需顺序推进的粗校正。对线性问题 $\boldsymbol u'+A\boldsymbol u=g$，误差按模态解耦，收缩由

$$
\varrho_l(J,z)=\frac{\bigl|R_g(z)-R_f^{J}(z/J)\bigr|}{1-|R_g(z)|},
\qquad z=\Delta T\lambda,\ \lambda\in\sigma(A)
$$

控制，其中 $J=\Delta T/\Delta t$ 是粗细比，$R_g,R_f$ 是粗细传播子的稳定函数。分子是粗细差，分母是粗传播的耗散裕量。

### 直接路线：全时间系统的对角化

把所有时间层的未知量堆成一个向量，一步法离散给出一个块下三角 Toeplitz 系统。若时间方向的矩阵可对角化 $B=VDV^{-1}$，则乘 $V^{-1}$ 后系统解耦为若干**独立的复移位空间问题**，可完全并发求解，再乘回 $V$。这条路线没有外层迭代，代价转移到 $V$ 的条件数上。

## 三组精读

| 精读页                                                                                                       | 论文               | 技术核心                   |
| ------------------------------------------------------------------------------------------------------------ | ------------------ | -------------------------- |
| [[computational-mathematics/paper-notes/parallel-in-time/parareal-convergence\|parareal 的收敛分析]]         | 12, 20, 30, 77     | 收缩因子与细传播子的稳定性 |
| [[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique\|对角化技术]]             | 31, 39, 46, 53     | 可对角化的时间矩阵         |
| [[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners\|全时间预条件与谱分析]] | 59, 65, 71, 84, 85 | 预条件系统的谱分布         |

## 逐篇核心思想

### parareal 的收敛分析

- **12（三种 parareal 求解器的收敛分析）** 处理一个此前没有理论覆盖的实际情形。当细传播子只是 A-稳定而非 L-稳定时——最重要的是梯形规则与四阶两级 Gauss Runge-Kutta 方法——$|R_f(z)|\to1$（$z\to\infty$），因此给出经典 $\approx0.3$ 常数的论证失效，$\lambda_{\max}$ 很大时收缩因子会趋近 $1$。此前的分析要么假设 $\mathcal F$ 是精确传播子，要么假设 L-稳定，对实践中重要的仅 A-稳定情形（尤其小粗细比 $J$）无话可说。
- **20（分数阶扩散方程的快速 parareal 迭代）** 与 **30（带局部时间积分器的 parareal）** 把这条路线用于时间分数阶问题。分数阶导数的非局部性使「每个粗区间内独立细传播」这个前提本身需要重新处理，这是这两篇的技术焦点。
- **77（非均匀细时间网格下的收敛分析）** 放弃「粗区间内细网格均匀」这一标准假设。这在实际计算中普遍存在，但会破坏上面的模态解耦所依赖的结构。

### 对角化技术

- **31（时间周期分数阶扩散方程）** 把对角化技术与多重网格结合。时间周期问题的全时间矩阵结构与初值问题不同，这一差别是可利用的。
- **39（两级 MGRIT 的对角化加速）** 同时修两处低效。第一，MGRIT 的粗传播子几乎总取后向 Euler，只有一阶，FCF 松弛下收缩约 $0.1$；论文改用二阶 Lobatto IIIC，它是 $(0,2)$ Padé 逼近，A-稳定且 L-稳定，对高频的衰减远强于后向 Euler 的 $1/(1+z)$，这正是粗传播子需要的性质。第二，也是结构性瓶颈：粗网格校正本质上是沿 $N_t$ 个粗时间层的**顺序**扫描，处理器数增加后它成为串行瓶颈；论文引入首尾耦合条件 $\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0$，把粗校正变成一个可对角化的全时间系统

  $$
  \bigl(C_\alpha\otimes I_x-I_t\otimes\Delta T A\bigr)\boldsymbol U^{k+1}=\boldsymbol g^k,
  $$

  其中 $C_\alpha$ 是 $\alpha$-循环矩阵，可由 FFT 对角化。两处修改是相容的：更好的粗传播子通常会加重顺序瓶颈，而并行化的粗校正正好抵消这一点。

- **46（抛物 PDE 约束最优化）** 把对角化用于最优控制的前向-后向最优性系统，那里时间方向的耦合是双向的。
- **53（两级 SDIRK 方法的并行实现）** 把对角化用在单个时间步的级之间，而不是时间层之间。

### 全时间预条件与谱分析

- **59（一、二阶演化方程的良条件直接 PinT 算法）** 直接给出全时间解而不做外层迭代，重点是让所涉矩阵的条件数受控——这正是几何步长对角化路线的代价所在。
- **65（预条件全时间系统的统一谱分析）** 对一阶与二阶演化问题给出统一的谱分析。谱分布决定 Krylov 迭代的收敛，因此这是把「构造一个预条件子」变成「证明它有效」的关键一步。
- **71（前向-后向演化方程的 PinT 预条件子）** 处理最优控制中出现的前向-后向耦合系统。
- **84（时间谱方法的 PinT 预条件子）** 时间谱离散给出的全时间时间块是**稠密**的，而不是前面几篇的双对角结构，因此预条件子的设计前提改变。
- **85（_Acta Numerica_ 综述）** 把上述工作与更广的文献统一组织为两类：对传播型问题仍然有效的方法，与主要为耗散问题设计的方法。

> [!note] 覆盖进度
> 编号 12、39、59、65 的公式已逐式核对。编号 30、31、46、53、84 的书目记录与摘要已核实，但定理中的具体常数未能核对。编号 20、77 的摘要与正文无法通过公开渠道获取（SIAM 不向 Crossref 提交摘要，其余聚合站对相应 DOI 返回空值），编号 71 的摘要在所有公开来源中都把行内公式整段丢失，因此这三篇的定理与常数本站不报告。编号 85 的完整逐节精读见[[computational-mathematics/knowledge-notes/time-parallelization/index|时间并行专题]]。

## 一条贯穿的权衡

对角化路线的代价可以精确说出来。用几何递增步长 $\tau_n=\tau_1\gamma^{n-1}$（$\gamma>1$）使时间矩阵的特征值互不相同，从而可对角化；但 $\gamma$ 越接近 $1$，特征向量矩阵 $V$ 的条件数越大，舍入误差被放大越严重；$\gamma$ 越大，末端时间步越粗，截断误差越大。这是一个真实的两难，并且解释了另一条路线的存在理由：等步长下用 $\alpha$-循环矩阵近似时间矩阵，对角化由 FFT 完成，条件数受 $\alpha$ 控制，代价是引入一个近似，因此需要外层 Krylov 迭代。**前者是直接解法，后者是预条件子**，这正是本专题两组工作的分工。

## 本专题的原文

编号与题录见[[computational-mathematics/paper-notes/catalog|论文总目录]]，各篇具体出处列在对应精读页末尾。综述本身的完整精读见[[computational-mathematics/knowledge-notes/time-parallelization/index|双曲与抛物问题的时间并行方法]]。
