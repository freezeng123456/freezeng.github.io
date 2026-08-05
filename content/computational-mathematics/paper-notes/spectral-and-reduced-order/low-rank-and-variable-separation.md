---
title: 低秩流形与变量分离
description: 编号 15、17、83：让基与系数同时随时间演化，代价写在最小奇异值上
lang: zh
translation: en/computational-mathematics/paper-notes/spectral-and-reduced-order/low-rank-and-variable-separation
tags:
  - 论文笔记
  - 降阶模型
  - 低秩逼近
---

> [!note] 本页覆盖
> 编号 **15**（_SIAM J. Sci. Comput._ 37(2), 2015）、**17**（_J. Comput. Phys._ 303, 2015）、**83**（_SIAM J. Sci. Comput._ 47(3), 2025）。编号 17 的期刊正文需订阅访问且无预印本，该节区分「本文可确认的内容」与「所改造的通用框架」；编号 83 的精读待补。

![在低秩流形上直接演化随机解](assets/diagrams/tao-zhou-papers/zh/low-rank-dynamics.svg)

## 15：动态正交逼近的误差理论

### 为什么需要演化的基

时间依赖随机 PDE 通常用**固定**基展开：随机变量方向用广义多项式混沌，空间方向用本征正交分解模态。困难在于解对随机参数的依赖会随时间显著改变，因此固定基需要不断增加项数才能维持精度；自适应或贪心修补（时间依赖多项式混沌、推广的本征正交分解）只能部分缓解。

Sapsis 与 Lermusiaux 的动态正交逼近让空间基与随机系数**同时**演化，从而规避这一增长，但此前几乎没有误差理论：已有结果针对相近的多构型时间依赖 Hartree 方法与动力低秩矩阵情形，而非随机 PDE。本文的第一个目标就是在动态正交逼近与动力低秩逼近之间建立精确对应，从而把已有理论搬过来。

### 拟设与规范条件

在有界 $D\subset\mathbb R^d$（$1\le d\le3$）与完备概率空间上考虑

$$
\frac{\partial u(x,t,\omega)}{\partial t}=\mathcal L\bigl(u(x,t,\omega),\omega\bigr),
$$

拟设为

$$
u_S(x,t,\omega)=\bar u_S(x,t)+\sum_{i=1}^{S}U_i(x,t)\,Y_i(t,\omega).
$$

表示的唯一性由三条条件固定：

$$
\mathbb E[Y_i(t,\cdot)]=0,
\qquad
\langle U_i(\cdot,t),U_j(\cdot,t)\rangle=\delta_{ij},
\qquad
\Bigl\langle \frac{\partial U_i(\cdot,t)}{\partial t},\,U_j(\cdot,t)\Bigr\rangle=0 .
$$

第三条是**规范条件**：基的时间导数被要求与基自身的张成空间正交，也就是不允许基在自己的张成空间内旋转。这一条消去的正是分解的多余自由度。

### 演化方程与三处不同的角色

把 Galerkin 投影与规范条件结合，得到

$$
\frac{\partial\bar u_S}{\partial t}=\mathbb E\bigl[\mathcal L(u_S)\bigr],
$$

$$
\sum_{i=1}^{S}C_{ij}(t)\,\frac{\partial U_i(x,t)}{\partial t}
=\Pi^{\perp}_{U}\,\mathbb E\bigl[\mathcal L(u_S)\,Y_j(t,\cdot)\bigr],
\qquad j=1,\dots,S,
$$

$$
\frac{\partial Y_i(t,\omega)}{\partial t}
=\bigl\langle \mathcal L^{*}(u_S(\cdot,t,\omega),\omega),\,U_i(\cdot,t)\bigr\rangle,
\qquad i=1,\dots,S,
$$

其中 $C_{ij}(t)=\mathbb E[Y_iY_j]$ 是协方差矩阵，$\mathcal L^{*}(u,\omega)=\mathcal L(u,\omega)-\mathbb E[\mathcal L(u)]$，而

$$
\Pi^{\perp}_{U}[v]=v-\sum_{i=1}^{S}\langle v,U_i\rangle U_i
$$

是到 $\mathcal U=\mathrm{span}\{U_1,\dots,U_S\}$ 的正交补的 $L^2(D)$ 投影。

三个方程的角色不同：均值方程就是平均后的 PDE；基方程由算子与随机系数的相关性中**被投影出张成空间的那部分**驱动，各模态之间仅通过 $C$ 耦合；系数方程是把中心化算子对每个基函数作检验得到的、以 $\omega$ 为参数的常微分方程。

### 切空间与 $C^{-1}$ 的来源

把中心化部分写成 $u_S^{*}=\mathbf U^{T}\mathbf Y$，它落在秩 $S$ 流形 $\mathcal M_S$ 上。切空间由 $\delta u_S^{*}=\mathbf U^{T}\delta\mathbf Y+\delta\mathbf U^{T}\mathbf Y$ 构成，在规范条件的类比 $\langle\delta\mathbf U,\mathbf U^{T}\rangle=\mathbf 0$ 与 $\mathbb E[\delta\mathbf Y]=\mathbf 0$ 下分解唯一。切空间上的正交投影为

$$
P_{u_S^{*}}(v)=\mathbf U^{T}\langle v^{*},\mathbf U^{T}\rangle
+\Bigl(\Pi^{\perp}_{U}\bigl\{\mathbb E[v^{*}\mathbf Y^{T}]\bigr\}\,C^{-1}\Bigr)^{T}\mathbf Y,
$$

其中 $C=\mathbb E[\mathbf Y\mathbf Y^{T}]$ 按秩 $S$ 的定义满秩。**这就是后来出现在误差界里的 $\rho^{-1}$ 的来源。** 互补投影可以因子化为 $P^{\perp}_{u_S^{*}}v=\Pi^{\perp}_{U}\otimes\Pi^{\perp}_{\tilde{\mathcal Y}}v$。

用 Dirac-Frenkel 形式写出来，动态正交逼近就是 Galerkin 条件

$$
\mathbb E\Bigl[\Bigl\langle \frac{\partial u_S}{\partial t}-\mathcal L(u_S),\,v\Bigr\rangle\Bigr]=0
\quad \forall v=\bar v+v^{*},\ v^{*}\in T_{u_S^{*}(t)}\mathcal M_S,
$$

等价地

$$
\frac{\partial u_S}{\partial t}=\mathbb E[\mathcal L(u_S)]
+P_{u_S^{*}(t)}\bigl(\mathcal L^{*}(u_S)\bigr).
$$

这个写法把「同时演化基与系数」这件事翻译成一句话：**把右端投影到流形的切空间上。**

### Theorem 4.1：拟最优性与它的代价

参照对象是每个时刻的最佳 $S$ 秩逼近，即截断的 Karhunen-Loève 展开

$$
z_S(x,\omega,t)=\bar u(x,t)+\sum_{i=1}^{S}\sqrt{\mu_i(t)}\,\gamma_i(t,\omega)\,Z_i(x,t),
$$

其中 $\{\mu_i,Z_i\}$ 是协方差算子的特征对。定理假设该最佳逼近连续可微存在于 $(H^2(D)\cap H^1_0(D))\otimes L^2(\Omega)$ 内，且其最小奇异值一致有下界

$$
\sigma\bigl(z_S(t)\bigr)\ \ge\ \rho\ >\ 0,\qquad \forall t\in[0,\bar t].
$$

则存在 $0<\hat t\le\bar t$，使动态正交解在 $u_S(0)=z_S(0)$ 下满足

$$
\|u_S(t)-z_S(t)\|_0^2+a_{\min}\int_0^t |u_S(\tau)-z_S(\tau)|_1^2\,\mathrm d\tau
\ \le\ 2\alpha\,e^{2\beta(t)}\int_0^t \|z_S(\tau)-u(\tau)\|_1^2\,\mathrm d\tau,
$$

$$
\beta(t)=4\rho^{-1}\int_0^t\Bigl(4\|\mathcal L^{*}(z_S)\|_0+\|\mathcal L^{*}(u)\|_0
+\|\mathcal L^{*}(u_S)\|_0+\|\dot z_S^{*}\|_0^2\Bigr)\mathrm d\tau,
\qquad
\alpha=\max\Bigl\{\frac{a_{\max}^2}{2a_{\min}},\ 4\rho^{-1}\Bigr\}.
$$

结论的读法有两层。**上层是好消息**：动态正交误差被最佳秩 $S$ 误差控制，方法是拟最优的。**下层是代价**：常数按 $e^{C/\rho}$ 增长，即随最小奇异值倒数指数增长。这就是曲率障碍的定量形式——流形 $\mathcal M_S$ 的曲率约为 $1/\sigma_S$，当第 $S$ 个奇异值塌缩时切空间投影变得病态，界失去意义。

### 一个具体的失效场景与一个实现细节

论文构造了一个特征值交叉的例子：确定性 Laplacian 加随机初值，前两个 Karhunen-Loève 特征值在某时刻 $t^{*}$ 交叉。此时最佳秩 1 误差与动态正交秩 1 误差都可以闭式算出，交叉前后由不同分支主导。这个例子说明上面的假设不是技术性的：秩不足时最小奇异值确实会塌缩。

实现层面还有一个容易忽略的问题。协方差矩阵 $C(t^n)$ 可能奇异或极度病态——例如任何带随机系数但**确定性初值**的系统在起始时刻都有 $C\equiv0$。直接用 Moore-Penrose 伪逆 $C^{\dagger}$ 会把「非活跃」的基函数置零，从而**阻止秩增加**。论文改为把基方程直接写成

$$
\frac{\partial\mathbf U}{\partial t}=C^{\dagger}\,\Pi^{\perp}_{U}\,
\mathbb E\bigl[\mathbf Y\,\mathcal L(u_S)\bigr],
$$

在 $C$ 满秩时它与原方程一致；实际计算中每步对协方差矩阵做对角化以解耦系统，因为即便 $C(0)$ 对角，动态正交流也不保持 $Y_i$ 的不相关性。

## 17：让特征求解只在最粗网格上进行

Karhunen-Loève 展开需要协方差核对应的 Fredholm 积分算子的特征对：给定 $\mathrm{Cov}(x,y)$，求

$$
\int_D \mathrm{Cov}(x,y)\,u(y)\,\mathrm dy=\lambda\,u(x),
\qquad \|u\|_{L^2(D)}=1 .
$$

与微分算子的特征问题相比，这里有一处结构性差别：积分算子 $\mathcal T$ 是**紧的、光滑化的**，特征值在 $0$ 处聚集（而不是趋于 $+\infty$），离散化给出**稠密**矩阵。因此既没有稀疏性可用，装配与求解的代价又随分辨协方差核所需的网格分辨率急剧上升。

先前工作从离散化一侧攻这个代价——小波 Galerkin、推广的快速多极、谱元、带张量结构的 Legendre-Galerkin、紧积分算子的多层增广方法——但特征求解本身仍必须在最细网格上完成。本文的贡献正是取消这一点：把特征值问题转化为一列**积分迭代**加上只在**最粗网格**上的特征求解，因此任何高效积分方案都可以插进来，总工作量与最细网格上的一次积分相当。

其所改造的多层修正框架（Lin 与 Xie）的通用一步是：给定最粗空间 $V_H$ 与嵌套层次 $V_H\subset V_{h_2}\subset\cdots$，先在细空间上解一个**源问题**

$$
a(\tilde u_{h_{k+1}},v)=\lambda_{h_k}\,b(u_{h_k},v),\qquad \forall v\in V_{h_{k+1}},
$$

再在小的增广空间 $V_{H,h_{k+1}}=V_H+\mathrm{span}\{\tilde u_{h_{k+1}}\}$ 内解特征值问题

$$
a(u_{h_{k+1}},v)=\lambda_{h_{k+1}}\,b(u_{h_{k+1}},v),\qquad \forall v\in V_{H,h_{k+1}} .
$$

第二步的意义在于：特征问题只在维数为 $\dim V_H+1$ 的空间里解，即最粗空间加一个方向，从不在细网格上解。对积分算子而言，第一步的「源问题」不是线性求解而是**一次求积**（即 $\mathcal T$ 的一次作用），这正是摘要所说的积分迭代。

> [!warning] 可核实范围
> 该文正文需订阅访问且无预印本。可确认的是问题设定、关键词（不确定性量化、Karhunen-Loève 展开、Fredholm 特征值问题、多重网格有限元）与摘要中的贡献表述。上面给出的通用修正步来自被改造的 Lin-Xie 框架；**该框架对积分算子的具体改造形式、常数与数值结果本站均未核实**，引用时应回查期刊版。

## 83：贪心地逐项加基，使每一步都解耦

编号 83 与编号 15 属于同一族问题——为参数依赖的动力系统构造随时间演化的分离表示，而不是在固定基上展开——但走的是一条结构上不同的路线。

编号 15 一次性演化一个秩 $S$ 的表示：全部模态同时推进，彼此通过协方差矩阵 $C$ 耦合，唯一性由规范条件固定。编号 83 反过来，用**贪心算法逐项加基**：每一步只往约化基里加一项，从而把问题重写成该步上的**两个解耦演化方程**，其中一个是**与参数无关的偏微分方程**，另一个是**依赖参数的常微分方程**。两个方程都直接由原动力系统与之前已分离出的表示项导出。

这条路线的收益可以对着编号 15 的误差界读出来。那里的常数按 $e^{C/\rho}$ 增长，$\rho$ 是最佳秩 $S$ 场的最小奇异值下界，而 $\rho^{-1}$ 的来源正是切空间投影里的 $C^{-1}$。**逐项贪心加基使每步的两个方程解耦，因此不出现需要求逆的协方差矩阵**，最小奇异值塌缩这一风险源随之不出现在同一位置。代价是贪心：逐项选出的序列不保证是最优的秩 $S$ 表示，而编号 15 的界恰恰是以最佳秩 $S$ 逼近为参照的。

第二处结构收益是计算被切成两段：**离线**阶段构造约化基函数，**在线**阶段使用得到的低秩表示。编号 15 没有这个划分——它的基与系数在同一次时间推进中一起演化，因此无法把代价前移。论文据此声称相对许多已有低秩分离技术降低了计算复杂度并提高了效率，并给出线性与非线性参数依赖动力系统的数值结果。

> [!note] 覆盖范围
> 上述问题设定、贪心逐项加基、每步两个解耦方程（一个与参数无关的偏微分方程加一个依赖参数的常微分方程）、以及离线在线划分，均可从预印本摘要确认。分离表示的具体形式、两个方程的显式写法、贪心准则与收敛性结果本站尚未逐式核对。
>
> 出版题名为 _A Dynamical Variable-Separation Method for Parameter-Dependent Dynamical Systems_，与主页所列略有差别，本站按出版版本记录；预印本为 [arXiv:2502.08464](https://arxiv.org/abs/2502.08464)（2025 年 2 月 12 日）。

## 三篇的关系

| 编号 | 演化对象                 | 模态之间的关系        | 主要风险                 |
| ---- | ------------------------ | --------------------- | ------------------------ |
| 15   | 空间基与随机系数         | 同时演化，经 $C$ 耦合 | 最小奇异值塌缩           |
| 17   | 不演化（一次性特征求解） | 不适用                | 稠密矩阵的装配与求解代价 |
| 83   | 逐项贪心加入的分离表示   | 每步解耦成两个方程    | 贪心序列不保证最优       |

编号 17 与编号 15 的关系值得指出：编号 15 的误差界以截断 Karhunen-Loève 逼近为参照，而编号 17 解决的正是**如何算出那个参照对象**。两篇一起给出了这条路线的完整代价结构：要么承担计算 Karhunen-Loève 基的稠密特征问题，要么让基自己演化并承担最小奇异值的风险。

编号 83 提供了第三种取法：既不预先算出最优基，也不让全部模态同时演化，而是逐项贪心地加入并让每步解耦。三篇合起来把这条线索的选择空间划清了——**耦合演化换来对最佳秩 $S$ 逼近的拟最优性，解耦贪心换来不出现协方差求逆，一次性特征求解换来最优基但要付稠密矩阵的代价。**

## 覆盖核对

| 内容                             | 论文 | 覆盖状态                                 |
| -------------------------------- | ---- | ---------------------------------------- |
| 拟设与三条规范条件               | 15   | 形式、条件与规范条件的含义               |
| 三个演化方程与各自的角色         | 15   | 均值、基、系数方程及投影算子             |
| 切空间投影与 $C^{-1}$ 的来源     | 15   | 切空间、投影公式、Dirac-Frenkel 形式     |
| Theorem 4.1 与 $e^{C/\rho}$ 常数 | 15   | 假设、结论、两层读法与曲率解释           |
| 特征值交叉的失效场景             | 15   | 例子构造与其说明的问题                   |
| 奇异协方差的处理                 | 15   | 伪逆的陷阱、改写形式、逐步对角化         |
| Fredholm 特征问题的结构差别      | 17   | 紧性、谱聚集、稠密矩阵                   |
| 多层修正步与积分迭代             | 17   | 通用两步、维数论证、积分与线性求解的差别 |
| 贪心逐项加基与每步的两个方程     | 83   | 加基方式、解耦结构、与编号 15 的对照     |
| 离线在线划分                     | 83   | 两阶段分工及编号 15 为何没有这一划分     |

## 本页原文

- E. Musharbash, F. Nobile, and T. Zhou, [_Error analysis of the dynamically orthogonal approximation of time dependent random PDEs_](https://doi.org/10.1137/140967787), SIAM J. Sci. Comput. 37(2) (2015), pp. A776-A810。
- H. Xie and T. Zhou, [_A multilevel finite element method for Fredholm integral eigenvalue problems_](https://doi.org/10.1016/j.jcp.2015.09.043), J. Comput. Phys. 303 (2015), pp. 173-184。
- L. Chen, Y. Chen, Q. Li, and T. Zhou, [_A dynamical variable-separation method for parameter-dependent dynamical systems_](https://doi.org/10.1137/24M168427X), SIAM J. Sci. Comput. 47(3) (2025), pp. A1783-A1808（预印本 [arXiv:2502.08464](https://arxiv.org/abs/2502.08464)）。
