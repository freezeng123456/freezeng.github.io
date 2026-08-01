---
title: 3.3：并行积分延迟校正
description: 从积分残差逐步推导 IDC 更新，并完整解释 PIDC、RIDC 与正则性限制
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-2-idc
tags:
  - 时间并行
  - IDC
  - 高阶时间积分
---

> [!note] 阅读范围
> 本页对应论文 Section 3.3（pp. 405–411），覆盖公式 (3.5)–(3.12)、Theorem 3.3 和 Figures 3.4–3.6。推导保留每个代数步骤，PIDC 与 RIDC 的并行调度分别说明。

## 3.3.1 IDC 从哪里开始

考虑非线性初值问题。积分形式为

$$
\boldsymbol u(t)=\boldsymbol u_0+
\int_0^t\boldsymbol f(\boldsymbol u(\tau),\tau)\,d\tau,
\qquad t\in(0,T]. \tag{3.5}
$$

先给出一个粗近似 $\widetilde{\boldsymbol u}(t)$。它可以是常函数 $\boldsymbol u_0$，也可以来自低精度积分。定义真实误差和积分残差

$$
\boldsymbol e(t)=\boldsymbol u(t)-\widetilde{\boldsymbol u}(t),
$$

$$
\boldsymbol r(t)=\boldsymbol u_0+
\int_0^t\boldsymbol f(\widetilde{\boldsymbol u}(\tau),\tau)\,d\tau
-\widetilde{\boldsymbol u}(t). \tag{3.6}
$$

$\boldsymbol r$ 衡量当前轨道违反积分方程的程度。将 $\boldsymbol u=\widetilde{\boldsymbol u}+\boldsymbol e$ 代入 (3.5)，再加减 $\int_0^t\boldsymbol f(\widetilde{\boldsymbol u},\tau)d\tau$，得到

$$
\begin{aligned}
\boldsymbol e(t)
&=\boldsymbol u_0+
\int_0^t\boldsymbol f(\widetilde{\boldsymbol u}(\tau)+\boldsymbol e(\tau),\tau)\,d\tau
-\widetilde{\boldsymbol u}(t)\\
&=\boldsymbol r(t)+
\int_0^t\left[
\boldsymbol f(\widetilde{\boldsymbol u}(\tau)+\boldsymbol e(\tau),\tau)
-\boldsymbol f(\widetilde{\boldsymbol u}(\tau),\tau)
\right]d\tau. \tag{3.7}
\end{aligned}
$$

对时间求导，得到误差方程

$$
\boldsymbol e'(t)-\boldsymbol r'(t)
=\boldsymbol f(\widetilde{\boldsymbol u}(t)+\boldsymbol e(t),t)
-\boldsymbol f(\widetilde{\boldsymbol u}(t),t). \tag{3.8}
$$

这一步把“改进解”转化为“求解误差方程”。残差积分使用高阶求积，误差方程本身可以继续采用简单时间步进器。

## 3.3.2 离散校正公式

取节点

$$
0=t_0<t_1<\cdots<t_M=T,
\qquad \boldsymbol u_m^k\approx\boldsymbol u(t_m),
$$

并令 $\Delta t_m=t_{m+1}-t_m$。对 (3.8) 使用线性 $\theta$ 方法，得到

$$
\begin{aligned}
\boldsymbol e_{m+1}-\boldsymbol e_m
={}&\boldsymbol r_{m+1}-\boldsymbol r_m\\
&+\Delta t_m(1-\theta)
\left[\boldsymbol f(\boldsymbol u_m^{k+1},t_m)
-\boldsymbol f(\boldsymbol u_m^k,t_m)\right]\\
&+\Delta t_m\theta
\left[\boldsymbol f(\boldsymbol u_{m+1}^{k+1},t_{m+1})
-\boldsymbol f(\boldsymbol u_{m+1}^k,t_{m+1})\right]. \tag{3.9}
\end{aligned}
$$

由残差定义，

$$
\boldsymbol r_{m+1}-\boldsymbol r_m
=\int_{t_m}^{t_{m+1}}\boldsymbol f(\boldsymbol u^k(\tau),\tau)\,d\tau
-\left(\boldsymbol u_{m+1}^k-\boldsymbol u_m^k\right).
$$

局部积分由全组节点上的 Lagrange 插值求积：

$$
\int_{t_m}^{t_{m+1}}\boldsymbol f(\boldsymbol u^k(\tau),\tau)\,d\tau
\approx\sum_{j=1}^{M}\omega_{m,j}
\boldsymbol f(\boldsymbol u_j^k,t_j), \tag{3.10a}
$$

$$
\omega_{m,j}=
\int_{t_m}^{t_{m+1}}
\prod_{\substack{i=1\\i\ne j}}^{M}
\frac{\tau-t_i}{t_j-t_i}\,d\tau. \tag{3.10b}
$$

再用 $\boldsymbol u_m^{k+1}=\boldsymbol u_m^k+\boldsymbol e_m$ 消去误差变量。逐节点校正公式为

$$
\begin{aligned}
\boldsymbol u_{m+1}^{k+1}
={}&\boldsymbol u_m^{k+1}
+\Delta t_m(1-\theta)
\left[\boldsymbol f(\boldsymbol u_m^{k+1},t_m)
-\boldsymbol f(\boldsymbol u_m^k,t_m)\right]\\
&+\Delta t_m\theta
\left[\boldsymbol f(\boldsymbol u_{m+1}^{k+1},t_{m+1})
-\boldsymbol f(\boldsymbol u_{m+1}^k,t_{m+1})\right]\\
&+\sum_{j=1}^{M}\omega_{m,j}
\boldsymbol f(\boldsymbol u_j^k,t_j),
\qquad m=0,\ldots,M-1. \tag{3.11}
\end{aligned}
$$

每个校正层 $k$ 都从左向右扫过 $m=0,1,\ldots,M-1$。$\theta=1$ 给出后向 Euler 型校正，$\theta=1/2$ 给出梯形型校正。

> [!tip] 公式 (3.11) 的三项分工
> 第一行沿用简单积分器推进新解；中间两项修正新旧轨道在局部端点上的动力学差；最后的求积项注入旧轨道在整个节点集上的高阶积分信息。校正能够升阶，原因就在于最后一项比基础步进器看到了更完整的时间区间。

## 3.3.3 Theorem 3.3：每轮能提高多少阶

若基础积分器阶数为 $p$，使用 $M$ 个等距节点，则第 $k$ 次校正后的误差阶为

$$
O\!\left(\Delta t^{\min\{M,(k+1)p\}}\right).
$$

因此每轮最多增加 $p$ 阶，并受求积的最高阶 $M$ 限制。后向 Euler 对应 $p=1$，梯形规则对应 $p=2$。Gauss–Lobatto 节点可达到 $2J-1$ 阶；这类 IDC 通常称为 spectral deferred correction（SDC），也是 PFASST 的核心部件。

长时间区间不适合用一个高阶多项式整体近似。论文将 $[0,T]$ 划分为窗口

$$
I_n=[T_{n-1},T_n],
\qquad n=1,\ldots,N_t.
$$

普通 IDC 先把 $I_n$ 的全部校正完成，再把末值传给 $I_{n+1}$。窗口之间串行，窗口内部的节点更新也按 $m$ 串行。

## 3.3.4 PIDC：按窗口铺开的流水线

PIDC 采用 Womble (1990) 的流水思想。$I_n$ 完成第一遍 sweep 后，已经得到一个粗末值 $\boldsymbol u_{n,M}^{1}$。$I_{n+1}$ 可以立即用它启动第一遍 sweep，同时 $I_n$ 继续第二遍。一般地，在 $I_n$ 上做第 $k$ 遍时，可以同时在 $I_{n+1}$ 做第 $k-1$ 遍，在 $I_{n+2}$ 做第 $k-2$ 遍，直到 $I_{n+k-1}$ 的第一遍。

![原论文 Figure 3.4：四个时间窗口上的 PIDC 流水线启动和稳态阶段](assets/papers/time-parallelization/source-figures/figure-3-4.svg)

Figure 3.4 取 $M=6$、$k_{\max}=4$。前四个阶段依次填满流水线；随后四个窗口上的四个 sweep 可以同时运行。黑色虚线记录已经走过的 sweep，带红色圆点的红线表示当前并行 sweep，黑色实线表示精确解。

PIDC 的每个窗口会收到不断变化的粗初值。校正次数增加时，误差不一定单调下降。这个现象需要通过正则性实验判断。

## 3.3.5 周期对流扩散实验与离散矩阵 (3.12)

论文使用周期边界和 $\Delta x=1/64$。半离散矩阵写成

$$
A=\frac{\nu}{\Delta x^2}A_{xx}
-\frac{1}{2\Delta x}A_x,
$$

其中

$$
A_{xx}=
\begin{bmatrix}
-2&1&&&1\\
1&-2&1&&\\
&\ddots&\ddots&\ddots&\\
&&1&-2&1\\
1&&&1&-2
\end{bmatrix},
\qquad
A_x=
\begin{bmatrix}
0&1&&&-1\\
-1&0&1&&\\
&\ddots&\ddots&\ddots&\\
&&-1&0&1\\
1&&&-1&0
\end{bmatrix}. \tag{3.12}
$$

取 $T=3$、窗口宽度 $\Delta T=0.1$、$M=5$，基础积分器为后向 Euler。第 $n$ 个窗口、第 $k$ 遍后的误差定义为

$$
\operatorname{err}_n^k=
\frac{\max_m\lVert\boldsymbol u_{\mathrm{ref}}^{n,m}
-\boldsymbol u_k^{n,m}\rVert_\infty}
{\max_{n,m}\lVert\boldsymbol u_{\mathrm{ref}}^{n,m}\rVert_\infty}.
$$

参考解由 MATLAB ODE45 计算，相对和绝对容差均为 $10^{-13}$。下一个窗口的初始轨道固定为上一窗口第一遍的末值：$\boldsymbol u_{n+1,m}^0\equiv\boldsymbol u_{n,M}^1$。

![原论文 Figure 3.5：IDC 与 PIDC 在不同正则性和黏性下的逐窗口误差](assets/papers/time-parallelization/source-figures/figure-3-5.svg)

Figure 3.5 的四个面板形成两组对照。

- (a)、(b) 使用 $\sigma=1000$ 的尖锐源项，解的时间正则性不足；(c)、(d) 使用 $\sigma=5$ 的平滑源项。
- (a)、(c) 取 $\nu=1$；(b)、(d) 取 $\nu=10^{-3}$。
- 每个面板同时给出初始误差、第一遍和第二遍之后的 IDC/PIDC 误差。

低正则性时，第二遍无法继续降低误差。平滑源项加大黏性后，IDC 与 PIDC 的第二遍都继续改善，二者接近。小黏性会削弱升阶效果；在 (d) 中，PIDC 第二遍明显落后于串行 IDC。论文据此判断，低正则双曲解并不适合依赖高阶校正的 PIDC。

## 3.3.6 RIDC：滑动求积窗口

RIDC 把并行粒度从“时间窗口”细化到“校正层”。第一个处理器用低阶积分器连续向前推进。它产生前 $M$ 个值后，第二个处理器开始第一层 IDC 校正。第二个处理器到达第 $M$ 步后继续前进，并将求积节点从 $1,\ldots,M$ 滑到 $2,\ldots,M+1$，再滑到 $3,\ldots,M+2$。第三个处理器在第二层数据足够后启动，并采用相同滑动规则。

稳定阶段中，每个处理器负责一个校正层，所有层同时推进相邻的时间步。流水线需要保存若干历史值，并处理启动和排空阶段。

![原论文 Figure 3.6：IDC 与 RIDC 在相同对流扩散数据上的逐窗口误差](assets/papers/time-parallelization/source-figures/figure-3-6.svg)

Figure 3.6 沿用 Figure 3.5 的 PDE、网格、源项和黏性设置，将 PIDC 换成 RIDC。结论保持一致：平滑解允许校正层升阶；尖锐源项和小黏性使高阶插值与求积失去预期收益。RIDC 改变调度和内存结构，没有消除 IDC 对正则性的依赖。

## 公式与图表覆盖核对

| 原文项目           | 本页位置 | 覆盖状态                                     |
| ------------------ | -------- | -------------------------------------------- |
| (3.5)–(3.8)        | 3.3.1    | 积分方程、残差、积分误差方程和微分误差方程   |
| (3.9)–(3.11)       | 3.3.2    | $\theta$ 离散、求积权重和最终 IDC 更新       |
| Theorem 3.3        | 3.3.3    | 基础阶数、校正次数和求积阶上限               |
| Figure 3.4         | 3.3.4    | 完整 PIDC 启动与稳态流水图                   |
| (3.12), Figure 3.5 | 3.3.5    | 周期差分矩阵、全部参数、误差定义与四面板结论 |
| Figure 3.6         | 3.3.6    | RIDC 滑动窗口和正则性实验                    |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3.3, pp. 405–411.
