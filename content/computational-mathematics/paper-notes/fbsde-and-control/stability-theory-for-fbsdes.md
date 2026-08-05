---
title: 离散化的稳定性理论
description: 编号 47、63：均方意义下的 Lax 等价定理，以及第一个全离散误差分析
lang: zh
translation: en/computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes
tags:
  - 论文笔记
  - 随机微分方程
  - 稳定性理论
---

> [!note] 本页覆盖
> 编号 **47**（_SIAM J. Numer. Anal._ 58(4), 2020）、**63**（_SIAM J. Numer. Anal._ 60(4), 2022）。两篇均无预印本且正文需订阅访问，因此本页把「摘要可确认的结论」与「按上下文重构的形式」分开标注。

## 47：把逐格式证明换成一条等价定理

### 问题

到 2020 年，这一组已经产出了大量 FBSDE 时间离散：后向 Euler、$\theta$ 格式族、Zhao-Zhang-Ju 多步格式、[[computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|编号 8]] 的多步族、延迟校正。每一个都有自己特设的收敛证明，每一个都有自己的稳定性条件（例如编号 8 中经验观察到的 $k\le6$ 根条件窗口）。缺少的是一个**单一框架**，在其中稳定性只定义一次，收敛性可以抽象地推出。

### 主定理

论文给出一般离散化族，并证明**均方意义下的 Lax 等价定理**：

> 对 FBSDE 而言，一个相容的离散化格式收敛，当且仅当它稳定。

这与线性演化方程的经典 Lax 等价定理完全对应，差别在于这里的范数是均方（$L^2(\Omega)$）的，而问题是非线性的（生成元 Lipschitz）。摘要同时说明「分析结果对已有数值格式的应用也做了讨论」，即抽象定理被实例化，用来恢复或加强后向 Euler、$\theta$ 格式与多步格式的已知收敛结果。

**这条定理给编号 8 的经验窗口提供了理论归宿**：那里的根条件正是「稳定性」这一半，而相容性由构造保证，两者合起来给出收敛。

### 框架覆盖的两族

被统一的两族之一是 $\theta$ 格式：

$$
Y^\pi_{t_k}=\mathbb E\bigl[Y^\pi_{t_{k+1}}\mid\mathcal F_{t_k}\bigr]
+\theta_1\Delta_n f\bigl(t_k,Y^\pi_{t_k},Z^\pi_{t_k}\bigr)
+(1-\theta_1)\Delta_n\,\mathbb E\bigl[f(t_{k+1},Y^\pi_{t_{k+1}},Z^\pi_{t_{k+1}})\mid\mathcal F_{t_k}\bigr],
$$

$$
Z^\pi_{t_k}=-\frac{1-\theta_2}{\theta_2}\mathbb E\bigl[Z^\pi_{t_{k+1}}\mid\mathcal F_{t_k}\bigr]
+\frac{1}{\theta_2\Delta_n}\mathbb E\bigl[Y^\pi_{t_{k+1}}\Delta W_k\mid\mathcal F_{t_k}\bigr]
+\frac{1-\theta_2}{\theta_2}\mathbb E\bigl[f(t_{k+1},Y^\pi_{t_{k+1}},Z^\pi_{t_{k+1}})\Delta W_k\mid\mathcal F_{t_k}\bigr],
$$

$\theta_1\in[0,1]$、$\theta_2\in(0,1]$。当 $\theta_1=\theta_2=1/2$（Crank-Nicolson 成员）时该格式二阶收敛，否则一阶。另一族是编号 8 的多步格式。

> [!note] 可核实范围
> 主定理的**内容**（相容加稳定等价于收敛）与覆盖声明（后向 Euler、$\theta$ 格式、各类多步格式）可从摘要确认。一般离散化族的确切写法、稳定性与相容性的精确范数定义、以及是否覆盖耦合而非仅解耦系统，本站未核实。上面给出的 $\theta$ 格式是从独立的第三方综述转录的标准形式。

## 63：第一个全离散误差分析

### 空间方向的误差此前被绕过

每个概率型 BSDE 格式有两个误差来源：倒向方程的**时间**离散，以及条件期望 $\mathbb E_{t_n}[\cdot]$ 的**空间**逼近——在 Markov 情形下它们是 $d$ 维高斯积分，而被积函数只在网格点上已知。

这一组此前的流水线（编号 8、19、25、41）用 **Gauss-Hermite 求积**计算这些积分，其节点 $x_n+\sqrt{\Delta t}\lambda_j$ **不落在**空间网格上，因此每个节点、每个时间步都必须做一次局部多项式插值。这一步有三重代价：它本身花时间；它把空间精度限制在插值阶；它使**全离散**（时间加空间）误差分析极其困难。后果是，在本文之前，几乎所有严格的 BSDE 误差分析都是**半离散**的：假设条件期望精确，只界定时间误差。

论文的自述是：这似乎是**首次尝试分析 BSDE 的全离散格式**，达到时间二阶收敛与空间指数收敛。

### 两个部件

时间方向用 $\theta$ 格式族（如上），$\theta_1=\theta_2=1/2$ 给出二阶，与论文声明的「时间二阶」一致。

空间方向用 **Sinc 逼近**。Sinc 求积的关键性质是它对适当类的函数具有指数收敛，且节点按 $\mathrm sinh$ 型变换分布，因此可以设计成与网格相容——这正是绕开插值步的方式。论文的关键词列出「倒向随机微分方程、误差估计、条件数学期望、Sinc-$\theta$ 格式」，与此一致。

**这一篇在整条线索中的位置**：编号 47 把时间方向的稳定性理论统一起来，而编号 63 补上空间方向，使误差分析第一次覆盖全离散格式。二者共用「先确立稳定性」的架构。

> [!note] 可核实范围
> 「首次全离散分析」「时间二阶、空间指数收敛」「用 $\theta$ 格式做时间离散、用 Sinc 逼近处理条件期望」均可从摘要确认。$\theta$ 格式的具体形式按同组独立来源转录，Sinc 求积的具体构造与误差常数本站未核实。

## 两篇的关系

| 编号 | 覆盖的误差来源 | 主要结论                      | 对其他篇的作用                                            |
| ---- | -------------- | ----------------------------- | --------------------------------------------------------- |
| 47   | 时间           | 均方 Lax 等价定理             | 给编号 8 的经验窗口理论归宿；使编号 68 的反向设计成为可能 |
| 63   | 时间 + 空间    | 首个全离散分析（二阶 + 指数） | 补上编号 25 之外的另一条空间方向改进                      |

一条一般判断：**当一族方法的收敛证明逐个特设时，值得停下来把「稳定性」定义清楚。** 编号 47 之后，「设计一个收敛的格式」变成「设计一个稳定的格式」，而后者是一个可以优化的目标——[[computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|编号 68]] 正是这样做的。

## 本页原文

- J. Yang, W. Zhao, and T. Zhou, [_A unified probabilistic discretization scheme for FBSDEs: stability, consistency, and convergence analysis_](https://doi.org/10.1137/19M1260177), SIAM J. Numer. Anal. 58(4) (2020), pp. 2351-2375。
- X. Wang, W. Zhao, and T. Zhou, [_Sinc-theta schemes for backward stochastic differential equations_](https://doi.org/10.1137/21M1444679), SIAM J. Numer. Anal. 60(4) (2022), pp. 1799-1823。
