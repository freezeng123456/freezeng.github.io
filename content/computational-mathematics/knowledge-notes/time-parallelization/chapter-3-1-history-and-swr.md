---
title: 3.1–3.2：历史脉络与 Schwarz 波形松弛
description: 从波形松弛与区域分解出发，完整推导 OSWR、有限步收敛和 tent pitching
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-1-history-and-swr
tags:
  - 时间并行
  - SWR
  - 双曲-PDE
---

> [!note] 阅读范围
> 本页对应论文 Sections 3.1–3.2（pp. 396–405），覆盖历史发展、公式 (3.1)–(3.4)、Theorems 3.1–3.2 和 Figures 3.1–3.3。原论文图像保持原样；正文、推导说明与 intuition 均重新撰写。

> [!info] 图表许可
> 原论文以 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) 发布。本页复制 Figure 3.1–3.3 的完整图形，并在相邻文字中给出中文图解和出处。

## 3.1 历史发展

### 四条方法线索怎样形成

论文把四类对双曲问题有效的 PinT 方法放在一起讨论，并用第二章的四个 PDE 检查它们的理论性质。

1. **Schwarz 波形松弛（SWR）**来自区域分解与波形松弛的结合。连续时空子问题的做法由 Gander (1999) 提出，Giladi and Keller (2002) 独立给出相近构造。其区域分解根源可追溯到 Schwarz (1870)，波形松弛则来自 Lelarasmee et al. (1982) 的电路仿真。Gander (1997) 同时研究了抛物与双曲情形，Gander et al. (1999) 引入 SWR 这一名称。后续工作覆盖非线性抛物问题、优化传输条件、双曲问题以及 Dirichlet–Neumann 和 Neumann–Neumann 变体。Ciaramella et al. (2023) 的 unmapped tent pitching（UTP）也建立在 SWR 上。
2. **并行积分延迟校正（IDC）**源于 Böhmer and Stetter (1984) 的演化问题校正思想。Dutt et al. (2000) 将其整理成能够逐轮提高阶数的积分器。Guibert and Tromeur-Dervout (2007) 提出 PIDC，Christlieb et al. (2010) 提出 RIDC。
3. **ParaExp**由 Gander and Güttel (2013) 提出。它把初值传播和源项响应分开计算，随后产生了线性实现和非线性迭代扩展。
4. **ParaDiag**从 Maday and Rønquist (2008) 的直接对角化方法发展而来。后续研究扩展到抛物、双曲和非线性问题，也形成波形松弛、Parareal、Krylov 预条件、插值以及 Sherman–Morrison–Woodbury 加 Krylov 等实现路径。

这段历史的作用是标明并行性的来源。SWR 利用空间接口上的整段波形，IDC 利用校正层之间的流水线，ParaExp 利用线性叠加，ParaDiag 利用全时间矩阵的谱结构。

## 3.2 Schwarz 波形松弛（SWR）

### 从 WR 到 SWR

### 经典空间区域分解的限制

传统做法先统一离散时间，再在每个时间步上用区域分解求解椭圆问题。所有子域必须等当前时间步的 DD 迭代收敛，才能进入下一步。各子域还要共享同一时间网格，局部自适应和不同积分器的空间受到限制。

经典波形松弛从半离散系统

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol f(t),
\qquad A=M+N
$$

出发，并执行

$$
\frac{d\boldsymbol u^k}{dt}-M\boldsymbol u^k
=N\boldsymbol u^{k-1}+\boldsymbol f,
\qquad \boldsymbol u^k(0)=\boldsymbol u_0.
$$

Jacobi 型 $M$ 让标量或块 ODE 同时求解；Gauss–Seidel 型分裂可配合红黑着色。循环约化还能增加并行层次。收敛速度完全取决于 $A=M+N$ 是否把强耦合保留在 $M$ 中。Nevanlinna (1989) 已指出，如何找到弱耦合分块是 WR 的核心难题；不合适的分裂可能极慢，甚至发散。

SWR 在空间离散前切分连续区域。每个子域独立求解完整时间窗上的 PDE，并交换人工边界处的时间函数。这样既能为各子域选择不同的时空离散，也能根据 PDE 的传播机制设计传输条件。优化后的 Robin、Ventcel 或卷积条件试图逼近连续问题的 Dirichlet-to-Neumann 映射。

![Schwarz 波形松弛在完整时间窗上交换界面波形](assets/diagrams/pint/zh/schwarz-waveform-relaxation.svg)

## 3.2.1 一阶抛物问题

### OSWR 的完整迭代

考虑区间 $(0,L)$ 上的对流扩散方程，记

$$
\mathcal L=\partial_x-\nu\partial_{xx}.
$$

取重叠子域 $\Omega_1=(0,\beta L)$、$\Omega_2=(\alpha L,L)$，其中 $0<\alpha\leq\beta<1$。重叠宽度为 $l=(\beta-\alpha)L$。论文的 Robin OSWR 写成

$$
\left\{
\begin{aligned}
\partial_tu_1^k+\mathcal Lu_1^k&=0,
&& (x,t)\in\Omega_1\times(0,T],\\
u_1^k(0,t)&=0,\\
\frac1p\partial_xu_1^k(\beta L,t)+u_1^k(\beta L,t)
&=\frac1p\partial_xu_2^{k-1}(\beta L,t)+u_2^{k-1}(\beta L,t),
\end{aligned}
\right.
$$

$$
\left\{
\begin{aligned}
\partial_tu_2^k+\mathcal Lu_2^k&=0,
&& (x,t)\in\Omega_2\times(0,T],\\
\frac1p\partial_xu_2^k(\alpha L,t)-u_2^k(\alpha L,t)
&=\frac1p\partial_xu_1^{k-1}(\alpha L,t)-u_1^{k-1}(\alpha L,t),\\
u_2^k(L,t)&=0.
\end{aligned}
\right. \tag{3.1}
$$

每一轮都恢复物理初值：$u_i^k(x,0)=u_0(x)$。$k=0$ 的两条界面波形可以任意给定。参数 $p>0$ 决定 Robin 传输；$p\to\infty$ 时导数项消失，得到经典 Dirichlet 交换。多子域推广沿每个人工边界重复相同构造。非线性版本保留迭代结构，并把 $\mathcal L$ 换成相应的非线性算子。

### Theorem 3.1：Robin 参数的 minimax 选择

Theorem 3.1 在无界空间、两个子域和连续时空层面分析 (3.1)。令 $l>0$ 为重叠宽度，

$$
y_0=\frac{l}{\nu},
\qquad
y=\frac{l\omega}{\nu},
\qquad
\omega\in\left[\frac\pi T,\frac\pi{\Delta t}\right].
$$

优化参数具有 $p^*=\widetilde p^*/\nu$ 的缩放。定义单频收敛因子

$$
R_0(y,\widetilde p,y_0)
=\frac{(y-\widetilde p)^2+y^2-y_0^2}
{(y+\widetilde p)^2+y^2-y_0^2}e^{-y},
$$

以及内部极值位置

$$
\bar y(y_0,\widetilde p)
=\sqrt{\frac{
y_0^2+2\widetilde p
+\sqrt{\widetilde p\left(-\widetilde p^3-4\widetilde p^2
+(4+2y_0^2)\widetilde p+8y_0^2\right)}}{2}}.
$$

若 $y_0<y_c$，其中 $y_c=1.618386576\ldots$，则 $\widetilde p^*$ 是下式的唯一解：

$$
R_0(y_0,\widetilde p^*,y_0)
=R_0\!\left(\bar y(y_0,\widetilde p^*),\widetilde p^*,y_0\right). \tag{3.2a}
$$

若 $y_0\geq y_c$，则 $\widetilde p^*$ 由

$$
y_0=\widetilde p^*\sqrt{\frac{\widetilde p^*}{4+\widetilde p^*}} \tag{3.2b}
$$

唯一确定。令

$$
y_{\min}=\frac{l\pi}{\nu T},
\qquad
y_{\max}=\frac{l\pi}{\nu\Delta t},
$$

则所有相关频率上的最坏收敛因子满足

$$
\rho:=\max_{y\in[y_{\min},y_{\max}]}
R_0(y,\widetilde p^*,y_0)
\leq
R_0\!\left(\bar y(y_0,\widetilde p^*),\widetilde p^*,y_0\right). \tag{3.2c}
$$

Dirichlet 交换对应 $p=\infty$，于是

$$
\rho\leq e^{-y_{\min}}
=\exp\!\left(-\frac{l\pi}{\nu T}\right). \tag{3.3}
$$

> [!tip] 推导 intuition：为什么出现等峰条件
> 对时间作 Fourier 变换后，每个频率 $\omega$ 都变成一个独立的空间界面误差传播问题。穿越重叠区给出 $e^{-y}$，Robin 条件给出前面的有理反射系数。优化 $p$ 等价于压低整个频率区间上最高的峰。当两个候选峰一样高时，继续降低其中一个会抬高另一个，公式 (3.2a) 正是这个 minimax 平衡条件。$y_0=y_c$ 是内部峰结构发生变化的分界点。

### Figure 3.1：连续理论与离散实验的距离

![原论文 Figure 3.1：OSWR 理论收敛因子和四子域迭代数](assets/papers/time-parallelization/source-figures/figure-3-1.svg)

论文使用 $L=8.2$、$T=5$、$\Delta t=0.01$、$\Delta x=0.02$、$l=2\Delta x$，空间用中心差分，时间用后向 Euler，初值为

$$
u_0(x)=e^{-10(x-L/2)^2}.
$$

Figure 3.1(a) 画出 Dirichlet 和优化 Robin 条件的理论因子。$\nu$ 越小，对流越占主导，因子越小。Figure 3.1(b) 将区间分成四个子域，从随机界面猜测开始，误差低于 $10^{-8}$ 时停止。实测迭代数同样随 $\nu$ 减小而下降。

在 $\nu=0.1$ 时，Dirichlet 与优化 Robin 实测分别需要 92 和 28 轮；二子域连续理论预测 32 和 4 轮。差异来自三项设置变化：理论使用无界域，实验使用有界域；理论只有两个子域，实验有四个；实验还离散了空间和时间。论文据此提醒，Figure 3.1(a) 描述趋势和理想因子，不能直接充当多子域离散实现的迭代数公式。Robin 多子域的完整收敛理论仍有缺口。

### 更精确的传输条件

Ventcel 条件用局部高阶算子逼近最优的非局部传输。论文给出的 Fourier 空间最优算子具有形式

$$
\partial_x-\frac{1}{2\nu}
\mathcal F^{-1}\!\left(1+\sqrt{1+4i\nu\omega}\right),
$$

其中 $i^2=-1$。若 $l=C_1\Delta x$、$\Delta t=C_1\Delta x^\beta$，渐近因子可写为 $\rho=1-O(\Delta x^\gamma)$，$\gamma$ 由 $\beta$ 决定。Wu and Xu (2017) 的卷积条件给出与网格无关的常数因子 $\rho=1-C$，$C\in(0,1)$，并适合含 Volterra 型非局部项的演化方程。代价是接口算子需要时间历史。

## 3.2.2 二阶双曲问题

### Dirichlet 条件也能有限步收敛

对波动方程，两个重叠子域上的 SWR 为

$$
\left\{
\begin{aligned}
\partial_{tt}u_1^k&=c^2\partial_{xx}u_1^k+g,
&& (x,t)\in\Omega_1\times(0,T],\\
u_1^k(x,0)&=u_0(x),
&\partial_tu_1^k(x,0)&=\widetilde u_0(x),\\
u_1^k(0,t)&=0,
&u_1^k(\beta L,t)&=u_2^{k-1}(\beta L,t),
\end{aligned}
\right.
$$

$$
\left\{
\begin{aligned}
\partial_{tt}u_2^k&=c^2\partial_{xx}u_2^k+g,
&& (x,t)\in\Omega_2\times(0,T],\\
u_2^k(x,0)&=u_0(x),
&\partial_tu_2^k(x,0)&=\widetilde u_0(x),\\
u_2^k(\alpha L,t)&=u_1^{k-1}(\alpha L,t),
&u_2^k(L,t)&=0.
\end{aligned}
\right. \tag{3.4}
$$

**Theorem 3.2.** 对 $0<\alpha<\beta<1$，只要

$$
k>\frac{Tc}{\beta-\alpha},
$$

第 $k$ 轮在两个人工界面上的误差便在整个 $(0,T)$ 内为零：

$$
u_1^k(\alpha L,t)-u(\alpha L,t)=0,
\qquad
u_2^k(\beta L,t)-u(\beta L,t)=0.
$$

该结论来自有限传播速度。一轮交换只能把正确界面数据沿特征线推进有限距离。每轮增加一个已知正确的时空锥；当锥的累计高度超过 $T$ 时，整段界面波形都已精确。多子域、高维分解和一维非线性守恒律也有相应结果。

### Figure 3.2：红黑 SWR 如何推进正确区域

![原论文 Figure 3.2：带大重叠的红黑 SWR 四阶段几何过程](assets/papers/time-parallelization/source-figures/figure-3-2.svg)

Figure 3.2 采用五个交错子域。计算顺序如下。

1. 在 $\Omega_1,\Omega_3,\Omega_5$ 三个红色子域上并行求解到 $T_1$。未知的内部界面可先给任意值。有限传播速度保证每个子域底部的三角形区域已经精确，外边界附近还各多出一小块精确区域。
2. 在 $\Omega_2,\Omega_4$ 两个黑色子域上并行求解到 $T_2$。它们从上一轮已正确的蓝色区域取得边界数据，因而生成两个正确的菱形区域。
3. 再回到红色子域，把时间推进到 $T_3$。
4. 红黑交替继续进行，正确区域最终覆盖到 $T_4$。

每个子域在正确 tent 上方还计算了暂时不可靠的区域，这些冗余工作换来了并行性。Nievergelt 的思想在这里以特征锥的几何方式出现。

### MTP、UTP 与 Figure 3.3

Mapped tent pitching（MTP）把倾斜 tent 映射成时空柱体，再用经典时间步进求解，最后映回原几何。它省去红黑 SWR 中的冗余区域，同时增加映射成本。映射可能导致阶数下降，因此需要专门设计积分器。

Unmapped tent pitching（UTP）直接执行红黑 SWR，也可以理解为全时间系统上的限制加性 Schwarz。它保留原坐标，不会引入映射导致的阶数下降。高维实现也可以沿用成熟的 RAS 基础设施。

![原论文 Figure 3.3：UTP 对二阶波动方程误差区域的逐轮消除](assets/papers/time-parallelization/source-figures/figure-3-3.svg)

Figure 3.3 从随机界面猜测产生的初始误差开始。(b)、(c)、(d) 分别显示第 4 次红更新、第 8 次黑更新和第 12 次红更新。蓝色零误差区沿 tent 逐层上升。实现不需要预先显式标出 tent；残差在哪个时间高度降为零，就能判断当前的有效 tent 高度，并据此自适应选择 $T_i-T_{i-1}$。

MTP 依赖有限传播速度，无法直接用于传播速度无限的抛物问题。SWR 与 UTP 仍可用于抛物方程，优化传输条件尤其有效。对弱扩散的对流问题，可以在每个时间 slab 增加一两轮迭代来修正跨 tent 影响。

## 公式与图表覆盖核对

| 原文项目                          | 论文小节 | 覆盖状态                                                      |
| --------------------------------- | -------- | ------------------------------------------------------------- |
| 历史发展，pp. 396–398             | 3.1      | 四类方法的来源与主要分支均已列出                              |
| WR 分裂与 SWR 动机                | 3.2 导论 | 连续区域分解、WR 分裂困难、OSWR 设计目标                      |
| (3.1)                             | 3.2.1    | 两个子域的 PDE、物理边界、Robin 交换和初值                    |
| (3.2a)–(3.2c), (3.3), Theorem 3.1 | 3.2.1    | 参数缩放、单频因子、两种参数区间、最坏频率界和 Dirichlet 极限 |
| Figure 3.1                        | 3.2.1    | 完整原图、参数、停止准则、理论与实验差异                      |
| Ventcel 与卷积条件                | 3.2.1    | 最优算子形式、渐近因子、非局部代价                            |
| (3.4), Theorem 3.2                | 3.2.2    | 双子域迭代、有限步条件与误差结论                              |
| Figures 3.2–3.3                   | 3.2.2    | 完整原图、红黑推进、MTP/UTP 和残差自适应                      |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 3.1–3.2, pp. 396–405.
