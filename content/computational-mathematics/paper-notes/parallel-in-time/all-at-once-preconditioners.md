---
title: 全时间预条件与谱分析
description: 编号 59、65、71、84、85：把「构造一个预条件子」变成「证明它的谱受控」
lang: zh
translation: en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners
tags:
  - 论文笔记
  - 时间并行
  - 预条件
---

> [!note] 本页覆盖
> 编号 **59**（_Adv. Comput. Math._ 48:16, 2022）、**65**（_SIAM J. Matrix Anal. Appl._ 43(3), 2022）、**71**（_SIAM J. Matrix Anal. Appl._ 44(4), 2023）、**84**（_J. Sci. Comput._ 103:82, 2025）、**85**（_Acta Numer._ 34, 2025）。编号 85 的逐节精读见[[computational-mathematics/knowledge-notes/time-parallelization/index|时间并行专题]]，本页只给出它在这条线索中的位置。编号 71 的摘要在所有公开来源中数学符号均被剥离，相应小节明确标出。

## 59：换一个时间离散，使 $V$ 天生条件良好

### Maday-Rønquist 路线的天花板

对 $u'(t)+Au(t)=g(t)$ 形成全时间系统 $\mathcal M\boldsymbol u:=(B\otimes I_x+I_t\otimes A)\boldsymbol u=\boldsymbol b$，若 $B=VDV^{-1}$ 则

$$
\mathcal M=(V\otimes I_x)\bigl(D\otimes I_x+I_t\otimes A\bigr)(V^{-1}\otimes I_x),
$$

中间一步是 $n$ 个完全解耦的空间求解。**障碍在于标准等步长离散给出的 $B$ 根本不可对角化**：等步长后向 Euler 的 $B$ 是下二对角 Toeplitz 矩阵，即单个 Jordan 块；多步法的 $B$ 是下三角 Toeplitz 矩阵，同样亏损。

Maday 与 Rønquist 的对策是取**互不相同**的步长 $\{\Delta t_j\}$，使 $n$ 个对角元互异从而 $B$ 可对角化。但此时舍入误差满足

$$
\texttt{舍入误差}=\mathcal O\bigl(\epsilon\,\mathrm{Cond}_2(V)\bigr),
\qquad \epsilon=\text{机器精度},
$$

而对几何递增步长 $\Delta t_j=\Delta t_1\tau^{\,n-j}$，参数 $\tau>1$ 被夹在两边：$\tau\to1$ 使 $B$ 接近亏损、$\mathrm{Cond}_2(V)$ 爆炸；$\tau\gg1$ 使步长指数增长、离散精度被破坏。两个误差平衡下来，方法被限制在约 $n\approx20$ 到 $25$ 个时间点，可达的并行度因此有一个上限。取消这个上限就是本文的目标。

### 边值方法：把不稳定的格式只在全时间层面使用

论文的时间离散是混合的：前 $n-1$ 步用中心（跳蛙）差分，只有最后一步用隐式 Euler：

$$
\begin{cases}
\dfrac{u_{j+1}-u_{j-1}}{2\Delta t}+Au_j=g_j, & j=1,2,\dots,n-1,\\[6pt]
\dfrac{u_n-u_{n-1}}{\Delta t}+Au_n=g_n. &
\end{cases}
$$

这是一个**边值方法**：它**不能**作为时间推进格式运行，因为中心格式那样用是不稳定的；只有一次性全时间求解才有意义。这一构造出自 Axelsson 与 Verwer（1985），他们证明即使最后一步只有一阶，同时得到的全部解仍具有**一致二阶精度**。Brugnano、Mazzia 与 Trigiante（1993）用迭代方法求解相应的全时间系统，本文改为用对角化**直接**求解。

关键收益是步长可以**等长**，因此完全没有步长比参数：

$$
B=\frac{1}{\Delta t}
\begin{bmatrix}
0&\tfrac12&&&\\
-\tfrac12&0&\tfrac12&&\\
&\ddots&\ddots&\ddots&\\
&&-\tfrac12&0&\tfrac12\\
&&&-1&1
\end{bmatrix}.
$$

这个 $B$ 的特征向量矩阵条件数可证为 $\mathrm{Cond}_2(V)=\mathcal O(n^2)$，即随时间点数多项式增长而不是指数增长。这就把 $n\approx20$ 的天花板换成了一个温和的多项式代价。

**这一篇的方法论值得单独指出：** 前面的工作都在给定时间离散后设法处理 $V$ 的条件数，而这一篇反过来——**为了让 $V$ 条件良好而重新选择时间离散**。代价是必须放弃时间推进的解释，接受一个只在全时间意义下成立的格式。

## 65：用「稳定」这一条经典性质换掉全部结构性假设

### 此前理论的形态

到 2022 年，块 $\alpha$-循环预条件已成为主流的时间并行方法，尤其对双曲问题，因为 $\mathcal P_\alpha^{-1}$ 的作用经 FFT 归约为各时间层上的独立求解。但它的理论是**逐情形堆积**起来的：每篇论文选定一个时间积分器，利用该积分器的特殊结构算出 $\mathcal P_\alpha^{-1}\mathcal K$ 的谱。综述对此的评价很直接：这些分析繁复，并且严重依赖时间推进矩阵的特殊性质，例如**稀疏性、Toeplitz 结构与对角占优**。

缺少的是一条覆盖全部积分器的定理，而且其假设应当是一条**经典**性质（稳定性）而不是某个特设的结构性质。第二个缺口是：二阶（波型）问题通常通过把 $u''=Au+g$ 改写成一阶系统来处理，这会把每个时间层的内存翻倍——对精细空间网格与高维问题不利——因此需要对二阶形式的两步离散做直接分析。

### 两个设定的统一写法

一阶情形把一般一步法写成两矩阵形式

$$
r_1(\Delta t A)\,\boldsymbol u_n=r_2(\Delta t A)\,\boldsymbol u_{n-1}+\tilde{\boldsymbol g}_n,
\qquad n=1,\dots,N_t,
$$

因此方法的稳定函数是 $r_1^{-1}(z)r_2(z)$（例如后向 Euler 取 $r_1(z)=I-z$、$r_2(z)=I$）。堆叠后得到全时间矩阵与其块 $\alpha$-循环预条件子

$$
\mathcal K=I_t\otimes r_1(\Delta tA)-B\otimes r_2(\Delta tA),
\qquad
\mathcal P_\alpha=I_t\otimes r_1(\Delta tA)-C_\alpha\otimes r_2(\Delta tA),
$$

其中 $B$ 是严格下移矩阵、$C_\alpha=B+\alpha e_1e_{N_t}^\top$。也就是说，**预条件子是把全时间矩阵的右上角放上 $\alpha$ 得到的**，即用 $\alpha$-循环矩阵替换（块）Toeplitz 矩阵。$\alpha=1$ 给出普通循环矩阵，$\alpha\to0$ 回到 $\mathcal K$ 本身。

二阶情形不做一阶化，而用对称两步法

$$
r_1(\Delta t^2A)\,\boldsymbol u_{n+1}-r_2(\Delta t^2A)\,\boldsymbol u_n
+r_1(\Delta t^2A)\,\boldsymbol u_{n-1}=\tilde{\boldsymbol g}_n,
\qquad n=1,\dots,N_t-1,
$$

相应的全时间矩阵与预条件子为

$$
\mathcal K=\tilde B\otimes r_1(\Delta t^2A)-B\otimes r_2(\Delta t^2A),
\qquad
\mathcal P_\alpha=\tilde C_\alpha\otimes r_1(\Delta t^2A)-C_\alpha\otimes r_2(\Delta t^2A),
$$

其中两步模板使 $\alpha$ 出现在**两个**角落位置。

### 主定理：假设只有稳定性

设 $\mathcal K$ 来自上述一步积分器，$\mathcal P_\alpha$ 是 $\alpha\in(0,1)$ 的块 $\alpha$-循环预条件子。若该积分器**稳定**，即

$$
\bigl|r_1^{-1}(z)\,r_2(z)\bigr|\le1
\qquad \forall z\in\sigma(\Delta t A)\subset\mathbb C^-,
$$

则预条件矩阵的每个特征值满足**与网格无关**的界

$$
\frac{1}{1+\alpha}\ \le\ \bigl|\lambda(\mathcal P_\alpha^{-1}\mathcal K)\bigr|\ \le\ \frac{1}{1-\alpha} .
$$

这条结论的价值在于把假设换掉了：从「时间推进矩阵稀疏、Toeplitz、对角占优」换成「方法是稳定的」。后者是每本教科书都会检验的性质，因此定理可以对任意稳定的一步法直接使用，而不必逐个重做谱分析。

$\mathcal P_\alpha^{-1}$ 的作用方式是 $C_\alpha=VDV^{-1}$，$V=\Gamma_\alpha^{-1}\mathbb F^*$，$\Gamma_\alpha=\mathrm{diag}(1,\alpha^{1/N_t},\dots,\alpha^{(N_t-1)/N_t})$，因此步骤 (a) 与 (c) 是（带缩放的）FFT，步骤 (b) 是 $N_t$ 个独立的复空间求解，而 $\mathrm{Cond}_2(V)\le1/\alpha$ 是小 $\alpha$ 的舍入代价。

> [!warning] 界的方向
> 综述中转述这条定理时把两个端点印反了（写成 $\frac{1}{1-\alpha}\le|\lambda|\le\frac{1}{1+\alpha}$），这对 $\alpha\in(0,1)$ 不可能成立，因为 $1/(1-\alpha)>1/(1+\alpha)$。上式给出的是正确方向；单通道的直接计算可以确认这一点。

## 71：前向-后向情形有多个 Toeplitz 块

在由前向-后向演化方程支配的问题中，主要计算代价是求解一个大型线性系统，其中心对象是**前向子问题**经时空离散后的全时间矩阵。因此高效求解器需要一个好的预条件子。前向-后向情形特有的困难是：前向与后向演化耦合且方向相反，二者都不能独立时间推进，得到的是鞍点型或 Schur 补型系统而不是单个块下三角 Toeplitz 系统。编号 65 的 $\alpha$-循环理论只覆盖单一前向演化。

论文的预条件子构造可以一句话概述：**把全时间矩阵中的 Toeplitz 矩阵替换为 $\alpha$-循环矩阵。** 注意这里是复数——「那些 Toeplitz 矩阵」——与一个二乘二块系统相符，其中若干个 Toeplitz 块（至少包括前向矩阵及其转置）各自被循环化。综述把这一构造的结果表述为：ParaDiag-II 作用在前向-后向系统上，产生了 Pearson、Stoll 与 Wathen 的匹配 Schur 补预条件子的**并行版本**。

论文的两类应用是 PDE 约束最优控制问题（抛物 KKT 系统）与抛物源识别问题（从观测恢复未知源，同样给出前向-后向最优性系统）。

> [!warning] 可核实范围
> 该文摘要在所有公开来源中（OpenAlex 倒排索引、Crossref 衍生聚合站、经检索中继获取的出版社页面）**内联数学符号均被删除**，产生形如「求解线性系统 ___ 通常是主要计算负担……其中 ___ 是所谓的全时间矩阵」这样的句子。散文部分可靠，符号无法从公开来源恢复。特别是 $\alpha$ 的取值定标律本站未核实，因此不报告任何关于 $\alpha$ 与 $N_t$ 或误差之间关系的公式。

## 84：没有 Toeplitz 结构可循环化时怎么办

时间谱方法用基函数（例如多项式）的组合逼近解，是空间谱离散的自然搭配，时间方向精度很高。但全部组合系数必须**一次性**通过求解一个全时间系统得到，完全没有时间推进的选项，这一点与编号 31 的时间周期问题相同。系统在时间矩阵与空间矩阵之间具有 Kronecker 张量结构，而症结在于：

**时间谱方法给出的矩阵 $M$ 是一个非结构化矩阵。**

这是精确的障碍。ParaDiag-II 一族此前的每篇论文都依赖时间矩阵是（块）**Toeplitz** 的——通常还是块二对角 Toeplitz——从而预条件子可以由一句话的配方得到：「把 Toeplitz 矩阵换成它的 $\alpha$-循环对应物」（见编号 65 与 71）。时间谱的 $M$ 稠密且不带 Toeplitz 结构，因此**没有可循环化的对象**，配方直接失效。论文的表述是这为实际计算带来显著挑战。

这一篇因此界定了整条路线的适用边界：$\alpha$-循环预条件不是一个普适技巧，它的前提是时间方向的 Toeplitz 结构；一旦时间离散换成谱型，就需要另找一个可对角化的近似对象。

## 85：把上述工作放进一个统一框架

编号 85 是 2025 年的 _Acta Numerica_ 综述，把上述工作与更广的文献统一组织为两类：对传播型问题仍然有效的方法（Schwarz 波形松弛、积分延迟校正、ParaExp、ParaDiag），与主要为耗散问题设计的方法（parareal、PFASST、MGRIT、对角化 parareal、时空多重网格）。本站对它的逐节精读——含原论文全部 48 个图表资产与可复现的 Python 实验——是一个独立专题：[[computational-mathematics/knowledge-notes/time-parallelization/index|双曲与抛物问题的时间并行方法]]。

## 五篇的位置关系

| 编号 | 被解决的障碍                                 | 手段                                                                |
| ---- | -------------------------------------------- | ------------------------------------------------------------------- |
| 59   | 等步长下 $B$ 不可对角化，几何步长下 $V$ 病态 | 换用边值方法，使 $\mathrm{Cond}_2(V)=\mathcal O(n^2)$               |
| 65   | 逐情形的谱分析依赖特设结构假设               | 只假设稳定性，得到 $[\frac{1}{1+\alpha},\frac{1}{1-\alpha}]$ 的模界 |
| 71   | 前向-后向系统不是单个 Toeplitz 系统          | 循环化系统中的多个 Toeplitz 块                                      |
| 84   | 时间谱矩阵没有 Toeplitz 结构                 | 界定并突破 $\alpha$-循环配方的前提                                  |
| 85   | 文献分散                                     | 按动力学类型统一组织                                                |

一条贯穿的判断：**这条路线的每一步进展都表现为「把一个假设换成一个更弱或更经典的假设」。** 编号 59 把「步长互异」换成「等步长但换格式」；编号 65 把「Toeplitz 加对角占优」换成「稳定」；编号 84 指出剩下的那个假设（Toeplitz 结构）在什么时候不再成立。

## 覆盖核对

| 内容                                 | 论文 | 覆盖状态                                                            |
| ------------------------------------ | ---- | ------------------------------------------------------------------- |
| 三步分解与 $B$ 不可对角化的原因      | 59   | 分解式、Jordan 块、下三角 Toeplitz 的亏损性                         |
| 几何步长的两难与 $n\approx20$ 天花板 | 59   | 舍入误差公式、$\tau$ 的两侧限制、并行度上限                         |
| 边值方法及其一致二阶精度             | 59   | 混合格式、不可作为时间推进、$B$ 与条件数增长                        |
| 一阶与二阶设定的统一两矩阵写法       | 65   | $r_1,r_2$、$\mathcal K$ 与 $\mathcal P_\alpha$、两个角落的 $\alpha$ |
| 主定理与其假设的替换                 | 65   | 稳定性条件、模界、方向核对、$\mathcal P_\alpha^{-1}$ 的作用         |
| 前向-后向系统的多块循环化            | 71   | 障碍、配方的复数形式、与匹配 Schur 补的关系（限定核实）             |
| 时间谱矩阵的非结构性                 | 84   | 一次性求解、Kronecker 结构、配方失效的精确原因                      |

## 本页原文

- J. Liu, X.-S. Wang, S.-L. Wu, and T. Zhou, [_A well-conditioned direct PinT algorithm for first- and second-order evolutionary equations_](https://doi.org/10.1007/s10444-022-09928-4), Adv. Comput. Math. 48 (2022), 16（预印本 [arXiv:2108.01716](https://arxiv.org/abs/2108.01716)）。
- S.-L. Wu, T. Zhou, and Z. Zhou, [_A uniform spectral analysis for a preconditioned all-at-once system from first-order and second-order evolutionary problems_](https://doi.org/10.1137/21M145358X), SIAM J. Matrix Anal. Appl. 43(3) (2022), pp. 1331-1353。
- S.-L. Wu, Z. Wang, and T. Zhou, [_PinT preconditioner for forward-backward evolutionary equations_](https://doi.org/10.1137/22M1516476), SIAM J. Matrix Anal. Appl. 44(4) (2023), pp. 1771-1798。
- C. Tang, S.-L. Wu, T. Zhou, and Y. Zhou, [_Parallel-in-time preconditioner for the time spectral methods_](https://doi.org/10.1007/s10915-025-02899-w), J. Sci. Comput. 103 (2025), 82。
- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time parallelization for hyperbolic and parabolic problems_](https://doi.org/10.1017/S0962492924000072), Acta Numer. 34 (2025), pp. 385-489。
