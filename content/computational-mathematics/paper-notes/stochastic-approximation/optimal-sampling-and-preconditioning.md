---
title: 最优采样与预条件
description: 编号 22、24、28、45、54：把「采样密度」与「正交性密度」解耦
lang: zh
translation: en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning
tags:
  - 论文笔记
  - 不确定性量化
  - 最优采样
---

> [!note] 本页覆盖
> 编号 **22**（_Math. Comput._ 86, 2017）、**24**（_SIAM J. Sci. Comput._ 39(3), 2017）、**28**（_SIAM J. Sci. Comput._ 40(1), 2018）、**45**（_SIAM Rev._ 62(2), 2020）、**54**（_J. Comput. Phys._ 430, 2021）。

![配点设计的统一流程](assets/diagrams/tao-zhou-papers/zh/sampling-design.svg)

## 22：稳定性因子是一个可以设计掉的东西

### 问题的精确形式

标准 Monte Carlo 最小二乘从正交性密度 $w$ 独立同分布采样。按 Cohen、Davenport 与 Leviatan 的分析，所需样本数由**稳定性因子**

$$
\frac{\|K(z)\|_\infty}{N}
$$

控制，其中 $K$ 是多项式空间 $\mathbb P_k$ 在 $L^2_w$ 中的再生核对角线：以 $\varphi$ 记 $w$-正交归一基向量，则 $K(z)=\varphi^T\varphi$。这个量与基的选取无关：任何正交变换 $\psi\leftarrow U\varphi$ 都不改变 $K$。

对许多常用权而言，这个因子随多项式次数增长而爆炸，从而迫使样本数超线性增长（常常远差于二次）。

**论文的洞察是：这个因子是一个错配造成的假象。** 采样所用的密度不必等于定义正交性的密度。把两者解耦——从多势论平衡测度采样，再用 Christoffel 函数重新加权——就能把有效稳定性因子降到它的最优值 $N$。

量 $N/K_k(z)$ 即（归一化的）**Christoffel 函数**，方法由此得名。

### 算法

论文给出的算法在紧域 $D$ 上是：

1. 从平衡测度 $\mu_D$ 独立同分布抽 $S$ 个样本 $\{z_s\}$；
2. 组装数据向量 $u$，$(u)_s=u(z_s)$；
3. 计算最小二乘权 $K$，$(K)_{s,s}=N/K(z_s)$；
4. 组装 $S\times N$ 的类 Vandermonde 矩阵 $V$，$(V)_{s,n}=\varphi_n(z_s)$；
5. 解 $c=\arg\min_{g}\bigl\|\sqrt{K}Vg-\sqrt{K}u\bigr\|$。

无界域有相应的变体。论文自己指出一个等价说法：**按 Christoffel 函数加权等价于把设计矩阵的各行归一化。** 这句话把「最优采样」翻译成一个纯线性代数的操作，也解释了为什么行范数一致是关键。

理论依据是 Christoffel 函数的渐近性质（Berman、Bloom-Levenberg 等）：在适当的多势论条件下，$\rho^{2k}(z)K^{(k)}_k(z)\,\mathrm dV(z)/N$ 弱收敛到平衡测度 $\mathrm d\mu_{D,Q}(z)$。也就是说，**Christoffel 函数在高次极限下自动「找到」正确的采样密度**，无需人为设定。

### 稳定性因子被换成了什么

把上面的定性说法写成定量的：对紧域 $D$，从平衡测度独立同分布采样并取权 $N/K(z)$，则对任意 $r>0$，只要

$$
\frac{S}{N\log S}\ \ge\ C\,\frac{1+r}{\lambda_{\min}(R)},
$$

该过程即以高概率稳定。**关键对比是：$\|K\|_\infty/N$ 这个因子被换成了 $1/\lambda_{\min}(R)$**，其中 $R$ 是相应的 Gram 矩阵。精度估计相应为

$$
\mathbb E\Bigl[\bigl\|f-T_L(\tilde\Pi_Sf)\bigr\|_w^2\Bigr]
\le\|f-\Pi f\|_w^2
+\frac{\varepsilon(S)}{\lambda_{\min}(R)}\|f-\Pi f\|_{\tilde w}^2
+\frac{8L^2}{S^r}
+4\kappa^2(R)\,d^2(f),
$$

其中 $\varepsilon(S)=\frac{2-2\log2}{(1+r)\log S}\to0$，$\kappa(R)=\lambda_{\max}(R)/\lambda_{\min}(R)$。

> [!warning] 论文对自己理论的保留意见
> 这一点值得如实记下，因为它容易在转述中被抹掉。上式最后一项 $4\kappa^2(R)d^2(f)$ **不随 $S\to\infty$ 消失**，因此就表面形式而言，这套理论比已有的标准 Monte Carlo 理论要**弱**。论文明确表示目前无法严格证明两者可比，只能给出两条较弱的支持：一是论证该界与 Monte Carlo 界在量级上相同，二是数值上表明加权最小二乘**经常**（而非总是）更优。论文另外指出，一维情形下 $\lambda_{\min}(R)$ 与 $\kappa(R)$ 在一般权下数值表现都很好。
>
> 换言之，这条路线的说服力此处更多来自实验而非定理，而下文编号 45 的诱导采样判据才把「线性区域」变成一条干净的定理。

## 45：把这条线索写成综述

编号 45 发表在 _SIAM Review_，是这一族的枢纽。它重述并统一编号 9、10、11 与外部（Rauhut、Rauhut-Ward、Yan-Guo-Xiu 等）结果，并引入势论词汇——Lebesgue 常数、Fekete 点、平衡测度、收缩因子。

论文给出的样本复杂度对照最有信息量：

- **从正交性（均匀）测度独立同分布采样**（Cohen-Davenport-Leviatan）：对任意 $r>0$，若 $M/\log M\ge C_rN^2$，则 $\Pr[|\!|\!|\hat A-I|\!|\!|\ge\frac12]\le2M^{-r}$。这是**二次**要求 $M\gtrsim N^2\log N$，且在任意**下**指标集与多维情形下同样成立。
- **从 Chebyshev 测度 Monte Carlo 采样**：要求降到 $M\sim N^{\log3/\log2}$，严格优于 $N^2$ 但仍超线性。
- **确定性 Weil 点集**（来自编号 9）：$M\ge C(d)N^2$ 给出唯一解与近最佳逼近；这个二次要求比 Chebyshev Monte Carlo 更强，其补偿是确定性。

论文同时明确指出实践与理论之间的落差：实践者通常取 $M\simeq cN$，$c$ 在 $2$ 到 $3$ 之间，即**线性**，而线性区域的确定性理论「尚未确定地可得」。

### 一条判据统摄全部采样密度

综述最要紧的一步是把上面几条结果收进同一条判据。设 $q$ 为采样密度、$\{v_n\}$ 为基，则法方程的稳定性只要求

$$
\frac{M}{\log M}\ \ge\ C(r+1)\,
\sup_{x\in D}\sum_{n=1}^{N}\Bigl(\frac{v_n(x)}{q(x)}\Bigr)^{2},
\qquad
C=\frac{2}{\log(27/8e)}\approx9.24,
$$

此时 $\|G-I\|_2\le\frac12$ 以至少 $1-2M^{-r}$ 的概率成立（证明是矩阵 Chernoff 论证）。

**于是问题变成一个纯粹的优化问题：选 $q$ 使那个上确界最小。** 而上确界的下界恰是 $N$，在 $q\propto\sum_n v_n^2$ 时达到——这正是**诱导采样**。代入即得

$$
\frac{M}{\log M}\ \ge\ C(r+1)\,N,
$$

即 $M\sim N\log N$，在 $N$ 上除对数因子外已是最优。

**这条判据只依赖 $N=\dim V$。** 它不依赖维数 $d$、不依赖域 $D$、不依赖权 $w$，甚至不依赖具体取了哪个 $N$ 维子空间。代价则很明确：必须从非标准的密度 $\rho$ 抽样，而 $\rho$ **确实**依赖 $(V,w,D)$。也就是说，这条线索的实质是**把问题的全部依赖从样本复杂度搬到了采样密度里**——复杂度变得与问题无关，而所有困难集中到「怎么从 $\rho$ 抽样」这一件事上。

**这一段是整条线索的坐标系。** 编号 22 与 45 指出的线性区域正是编号 28 用完全不同的手段去够的目标。

## 28：用贪心选点代替随机采样

编号 28 提出**加权近似 Fekete 点**。Fekete 点是使 Vandermonde 行列式绝对值最大的点集，与 Lebesgue 常数直接相关，但精确计算是一个困难的非凸问题。近似 Fekete 点用**列主元 QR 分解**在一个候选点集上贪心地选点：QR 的主元顺序恰好近似地最大化行列式。

「加权」指把 Christoffel 权并入候选矩阵，因此选点在加权的意义下最优。这条路线的收益是样本数可以压到接近 $N$；代价是它需要一个候选池，并且失去了随机采样的独立性结构（因此不再有 with-high-probability 型的界，而要用势论的确定性论证）。

### 一维情形下贪心恰好就是最优

这一篇有一个比「贪心是个好启发式」强得多的结论，值得单列。取 $\Gamma=\mathbb R$ 上任意概率密度 $\rho$、$\Lambda=\{0,\dots,N-1\}$，记 $\phi_N$ 为 $N$ 次正交归一多项式，并定义亚纯函数

$$
r_N(y)=\frac{\phi_N(y)}{\phi_{N-1}(y)} .
$$

则对任意 $y\notin\phi_{N-1}^{-1}(0)$，集合 $A_N(y)=r_N^{-1}\bigl(r_N(y)\bigr)$（取集值逆）是**唯一确定**的，且它承载一个 $N$ 点正权求积公式，对 $2N-2$ 次以下精确，其权恰是 **Christoffel 权**：

$$
\int_\Gamma p(z)\rho(z)\,\mathrm dz=\sum_{z\in A_N(y)}\frac{1}{K_\Lambda(z)}\,p(z),
\qquad \deg p\le2N-2 .
$$

特别地，$y\in\phi_N^{-1}(0)$ 时 $A_N(y)$ 就是 $N$ 点 Gauss 公式。**这里出现的 $1/K_\Lambda(z)$ 与编号 22 的权 $N/K(z)$ 是同一个对象**，两条看起来不同的路线在这一点上合流。

在此基础上，论文证明：以这样的 $A_N(y)$ 为起点，贪心选点与全局最优选点**完全重合**，且

$$
\bigl|\det V\bigr|=1=\kappa(V),
$$

即所得 Vandermonde 矩阵是**完美条件**的。也就是说一维情形下贪心不是近似，而恰好取到最优。更要紧的是这条结论对**任意**一维密度成立，**包括非紧支撑的密度**——而这正是近似 Fekete 点所不具备的性质，因此它也是对无界域困难的一个正面回答。

> [!note] 一处印刷错误
> 该定理正文把起点条件印作 $y\in\phi_{N-1}^{-1}(0)$，这与它所依赖的引理 3.1 及推论 3.1 相矛盾（那里要求 $y$ **不在**零集中），应为 $y\notin\phi_{N-1}^{-1}(0)$。本页按后者叙述并在此标出。

**编号 22、45 与 28 因此构成一个完整的谱**：

| 路线                 | 样本要求             | 结论类型    | 需要的额外结构     |
| -------------------- | -------------------- | ----------- | ------------------ |
| 从正交性测度采样     | $M\gtrsim N^2\log N$ | 概率        | 无                 |
| Christoffel 加权采样 | 对数线性             | 概率        | 平衡测度 + $K(z)$  |
| 贪心（近似 Fekete）  | 接近 $N$             | 确定性/势论 | 候选池 + 列主元 QR |

## 24：把采样与预条件作为一对来设计

编号 24 面向多项式混沌的**稀疏**逼近，给出推广的采样与预条件方案。这里的关键认识是：在稀疏恢复的框架下（有界正交系与限制等距性质），决定样本复杂度的量同样是行范数的均匀性，因此「采样密度」与「预条件矩阵」应当作为一对来设计，而不是先定采样再补预条件。

论文明确把框架的适用范围写成「有界或无界」域，因此它同时回应了[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 11]] 在最小二乘情形遇到的无界域困难。

### 无界域的代价是可以写出来的

主恢复定理的样本要求形如

$$
M\ \ge\ L(n)\,\bigl\|R^{-1/2}\bigr\|_1^2\;s\log^3(s)\log(N),
$$

而**有界与无界的差别全部落在因子 $L(n)$ 上**：

- 有界（Chebyshev 型）情形：$L(n)=C(\alpha,\beta)$ 关于 $n$ **一致**有界，因此要求就是 $M\gtrsim s\log^3(s)\log(N)$——在稀疏度 $s$ 上近乎线性；
- 无界情形：$L(n)=Cn^{\max\{1/\alpha,\,2/3\}}$，例如 $\alpha\ge\frac32$ 时为 $Cn^{2/3}$，于是要求退化为 $M\gtrsim s\,n^{2/3}$，多出一个随多项式次数增长的因子。

结论求解的是加权的 $\ell_1$ 问题 $\min_\alpha\|R^{1/2}\alpha\|_1$ 受约束于 $\|\sqrt W\Phi\alpha-\sqrt Wf\|_2\le\varepsilon$，误差界的两项分别以 $\lambda_{\min}(R)$ 与 $\|R^{-1/2}\|_1$ 为系数——**与编号 22 一样，$R$ 的谱是最终出现在界里的量。**

> [!note] 一处内部不一致
> 论文的第三类情形（记作 CSA-c）中，$\max$ 表达式写作 $\max\{1/(2\alpha),2/3\}$，而紧随其后分支展开时印的是 $Cn^{1/\alpha}$（而非 $Cn^{1/(2\alpha)}$）。本页按原文照录并在此标出，不作静默更正。

## 54：把同一思路搬到再生核空间

编号 54 处理**核插值的最优设计**。多项式空间被换成再生核 Hilbert 空间后，Christoffel 函数的角色由核的对角线（幂函数）承担，而「按空间的集中程度加权采样」这一原则不变。论文把这一设计问题与不确定性量化的应用连起来。

**这一篇说明了这条线索的抽象层次**：Christoffel 加权采样的实质不是关于多项式的，而是关于「一个有限维函数空间在各点的集中程度」的。凡是有再生核的地方，同一设计原则就适用。

> [!note] 覆盖进度
> 编号 22 与 45 的核心构造、算法与样本复杂度对照已核对。编号 24、28、54 的具体定理、常数与实验本站未逐项核实；上述内容限于可从摘要、算法描述与参考关系确认的部分。

## 覆盖核对

| 内容                          | 论文 | 覆盖状态                                          |
| ----------------------------- | ---- | ------------------------------------------------- |
| 稳定性因子与其基无关性        | 22   | $K(z)=\varphi^T\varphi$、$\|K\|_\infty/N$、不变性 |
| 解耦采样密度与正交性密度      | 22   | 洞察表述与 Christoffel 函数的定义                 |
| 算法五步与行归一化的等价说法  | 22   | 完整步骤与等价解读                                |
| Christoffel 渐近与平衡测度    | 22   | 弱收敛结论及其含义                                |
| 三种采样的样本复杂度对照      | 45   | 二次、$N^{\log3/\log2}$、确定性二次               |
| 理论与实践之间的线性落差      | 45   | 实践取 $M\simeq cN$ 与理论现状                    |
| 加权近似 Fekete 点与列主元 QR | 28   | 选点机制、加权含义、代价与收益                    |
| 采样与预条件成对设计          | 24   | 稀疏恢复框架下的动机与适用范围                    |
| 从多项式空间到再生核空间      | 54   | 设计原则的抽象层次                                |

## 本页原文

- A. Narayan, J. Jakeman, and T. Zhou, [_A Christoffel function weighted least squares algorithm for collocation approximations_](https://doi.org/10.1090/mcom/3192), Math. Comput. 86 (2017), pp. 1913-1947。
- J. Jakeman, A. Narayan, and T. Zhou, [_A generalized sampling and preconditioning scheme for sparse approximation of polynomial chaos expansions_](https://doi.org/10.1137/16M1063885), SIAM J. Sci. Comput. 39(3) (2017), pp. A1114-A1144。
- L. Guo, A. Narayan, L. Yan, and T. Zhou, [_Weighted approximate Fekete points: sampling for least-squares polynomial approximation_](https://doi.org/10.1137/17M1140960), SIAM J. Sci. Comput. 40(1) (2018), pp. A366-A387。
- L. Guo, A. Narayan, and T. Zhou, [_Constructing least-squares polynomial approximations_](https://doi.org/10.1137/18M1234151), SIAM Rev. 62(2) (2020), pp. 483-508。
- A. Narayan, L. Yan, and T. Zhou, [_Optimal design for kernel interpolation: applications to uncertainty quantification_](https://doi.org/10.1016/j.jcp.2020.110094), J. Comput. Phys. 430 (2021), 110094。
