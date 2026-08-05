---
title: 对角化技术
description: 编号 31、39、46、53：让时间方向的矩阵可对角化，并为此付出条件数
lang: zh
translation: en/computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique
tags:
  - 论文笔记
  - 时间并行
  - 对角化
---

> [!note] 本页覆盖
> 编号 **31**（_Numer. Linear Algebra Appl._ 25(5), 2018）、**39**（_SIAM J. Sci. Comput._ 41(5), 2019）、**46**（_ESAIM Control Optim. Calc. Var._ 26, 2020）、**53**（_J. Comput. Phys._ 428, 2021）。四篇均无预印本；编号 46、53 的正文未能获取，相应小节区分摘要可确认的内容与所依赖的通用机制。

![把串行时间递推换成可对角化的时间矩阵](assets/diagrams/tao-zhou-papers/zh/pint-diagonalization.svg)

## 三步模式与它的价格

把全部时间层的未知量堆成一个向量，得到 Kronecker 形式的全时间系统

$$
\bigl(B_1\otimes I_x+B_2\otimes A\bigr)\boldsymbol u=\boldsymbol b,
\qquad
\boldsymbol u=(U_1^\top,\dots,U_{N_t}^\top)^\top .
$$

若时间矩阵可对角化 $C^{(\alpha)}_j=VD_jV^{-1}$，求解就变成三步：

$$
\text{(a) } S_1=(\mathbb F\Gamma_\alpha\otimes I_x)\boldsymbol b,
\qquad
\text{(b) } S_{2,n}=(\lambda_{1,n}I_x+\lambda_{2,n}A)^{-1}S_{1,n},
\qquad
\text{(c) } \boldsymbol u=(\Gamma_\alpha^{-1}\mathbb F^*\otimes I_x)S_2 ,
$$

其中 $\mathbb F$ 是单位 DFT 矩阵、$\Gamma_\alpha=\mathrm{diag}(1,\alpha^{1/N_t},\dots,\alpha^{(N_t-1)/N_t})$。步骤 (a) 与 (c) 是 FFT 与其逆，步骤 (b) 在 $N_t$ 个时间层上完全并行。

价格写在一处：$\mathrm{Cond}_2(V)\le1/\alpha$。$\alpha$ 越小，$\alpha$-循环矩阵越接近原来的下三角 Toeplitz 矩阵（近似越好），而 $V$ 的条件数越差（舍入放大越严重）。本页四篇论文可以按它们如何处理这个价格来排列。

## 31：时间周期问题让价格归零

### 障碍变成机会

时间周期扩散方程对时间推进而言结构上不同：没有初值可以出发，而周期条件 $u(0)=u(T)$ 把最后一个时间层耦合回第一个。摘要的说法是，必须**一次性考虑全部离散解**，而不是逐个求解。因此无论是否以并行为目标，都被迫进入全时间表述。再加上分数阶 Laplacian 使空间算子稠密，时空系统既大又贵。

论文的观察是：这个表面上的障碍实际上是机会。周期耦合使时间离散矩阵成为**真正的循环矩阵**，而不是下三角 Toeplitz 矩阵，因此可由 FFT 酉对角化。用 $\alpha$-循环族的语言说，这对应特殊值 $\alpha=1$：

$$
C_1^{(\alpha)}=\frac{1}{\Delta t}
\begin{bmatrix}1&&&-\alpha\\-1&1&&\\&\ddots&\ddots&\\&&-1&1\end{bmatrix},
\qquad
C_2^{(\alpha)}=
\begin{bmatrix}\theta&&&(1-\theta)\alpha\\1-\theta&\theta&&\\&\ddots&\ddots&\\&&1-\theta&\theta\end{bmatrix}.
$$

在 $\alpha=1$ 处 $\Gamma_\alpha=I$、$V=\mathbb F^*$ 是**酉矩阵**，故 $\mathrm{Cond}_2(V)=1$，限制其他对角化路线的舍入障碍完全消失。这就是时间周期情形允许一个**直接**的时间并行算法、且没有参数权衡的原因。

### 复移位系统的多重网格

步骤 (b) 留下一列独立的**复系数**线性系统 $(\lambda_{1,n}I_x+\lambda_{2,n}A)x=y$，其中 $A$ 是离散分数阶 Laplacian。论文用多重网格求解，光滑子取**带阻尼的 Richardson 迭代**。这个选择有具体理由：分数阶 Laplacian 矩阵稠密，因此只需要矩阵向量乘的光滑子优于需要三角求解的阻尼 Jacobi 或 Gauss-Seidel。

论文证明该线性求解器具有与网格无关的收敛因子，并对 Richardson 阻尼参数做优化以最小化这个常数因子。

> [!note] 可核实范围
> 上述问题设定、直接对角化、复移位系统的多重网格、与网格无关的收敛因子、以及阻尼参数优化的**存在**均可从摘要确认。最优阻尼值、最小化后的因子、分析是否关于分数阶阶数或 $N_t$ 一致，本站均未核实。

## 39：同时修好粗传播子与粗校正

编号 39 处理两级 MGRIT。它的收敛因子比 parareal 多一个因子：

$$
\varrho_l(J,z)=\frac{\bigl|R_f^{J}(z/J)\bigr|\,\bigl|R_g(z)-R_f^{J}(z/J)\bigr|}{1-|R_g(z)|}
\qquad\Longleftrightarrow\qquad
\varrho_{l,\mathrm{MGRIT}}=\bigl|R_f^J(z/J)\bigr|\times\varrho_{l,\mathrm{parareal}} .
$$

多出的 $|R_f^J(z/J)|\le1$ 正是额外一次 F-松弛带来的收益。术语上：F-松弛是在每个粗区间内部的 F 点上应用细传播子，从该区间左端的 C 点值出发，这些求解相互独立、可并行；C-松弛是用最后一个 F 点更新 C 点值；FCF-松弛是 F、C、F 依次进行，因此每轮代价是两次细求解。仅用 F-松弛的两级 MGRIT **就是** parareal。

### 两处低效

第一，粗传播子 $\mathcal G$ 几乎总取后向 Euler，只有一阶，FCF-松弛下收缩约 $0.1$（精确值见下）；此前没有人问过更好的 $\mathcal G$ 能否显著改进而**不让粗求解主导**。论文改用二阶 Lobatto IIIC 方法，其稳定函数为

$$
R_g(z)=\frac{1}{1+z+z^2/2}
$$

（对 $\boldsymbol u'+A\boldsymbol u=g$ 的约定）。这是 $(0,2)$ 阶 Padé 逼近：A-稳定、**L-稳定**（$R_g(\infty)=0$）、二阶，并且在阻尼意义上刚性精确——它对高频模态的衰减远强于后向 Euler 的 $1/(1+z)$，这正是粗传播子需要的性质。

第二，也是与 parareal 共有的结构瓶颈：粗网格校正本质上是沿 $N_t$ 个粗时间层的**顺序**扫描，处理器数增加后它按 Amdahl 定律成为串行瓶颈。而更好的 $\mathcal G$ 通常会加重这一点，因为高阶隐式 Runge-Kutta 粗求解器有多个级，代价是后向 Euler 步的若干倍。

### 首尾耦合把粗校正变成可对角化的系统

论文用一个**首尾耦合条件**替换顺序粗校正：

$$
\boldsymbol u_{n+1}^{k+1}=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\underbrace{\mathcal F(T_n,T_{n+1},\tilde{\boldsymbol u}_n^{k})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k})}_{=:\,\boldsymbol b_{n+1}^k},
\qquad
\boldsymbol u_0^{k+1}=\alpha\,\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0 ,
$$

其中 $\tilde{\boldsymbol u}_n^k=\boldsymbol u_n^k$（$n\ge1$）而 $\tilde{\boldsymbol u}_0^k=\boldsymbol u_0$——这一重定义是必要的，否则不动点不再是真解。取 $\mathcal G$ 为后向 Euler 时，耦合后的粗校正就是全时间系统

$$
\bigl(C_\alpha\otimes I_x-I_t\otimes\Delta T A\bigr)\boldsymbol U^{k+1}=\boldsymbol g^k ,
$$

由 FFT 对角化。

两处修改是相容的，这一点值得强调：单独改进粗传播子会加重顺序瓶颈，而并行化粗校正正好抵消这一加重。也就是说，两个改动各自的收益都不如合在一起时明显。

### 四个可直接比较的常数

设 $\mathcal F$ 为 L-稳定积分器、$J=\Delta T/\Delta t=O(1)$，则 $\max_{z\ge0}\varrho_l$ 为

| 粗传播子 $\mathcal G$ | Parareal | MGRIT（FCF-松弛） |
| --------------------- | -------- | ----------------- |
| 后向 Euler            | 0.2984   | 0.1115            |
| Lobatto IIIC（二阶）  | 0.0817   | 0.0197            |

这就是摘要中「从 0.1 降到 0.02」的精确形式。而且这条界是**稳健**的：它既不依赖系数矩阵的特征值，也不依赖粗细比 $J$。

有一处比较必须写清楚，否则容易高估收益：**按等代价折算，一次 MGRIT-FCF 迭代略差于两次 parareal 迭代。** FCF-松弛每轮做两次细求解，所以对应的比较对象是两轮 parareal，而 $0.2984^2=0.0890<0.1115$、$0.0817^2=0.0067<0.0197$。因此这两列常数不能直接对比来判断优劣——它们的迭代代价不同。真正稳固的结论是那个**换 $\mathcal G$** 的收益：同一列内从后向 Euler 换到 Lobatto IIIC，收缩因子降低约 3.7 倍（parareal）与 5.7 倍（MGRIT-FCF）。

至于耦合参数 $\alpha$ 的取法，论文的结论是适当选取时新算法与原算法**收敛率相同**。parareal 情形有一个明确阈值（归于 Wu, SISC 2018）：只要

$$
\alpha\le\frac{\rho}{1+\rho},
$$

就有 $\rho_{\text{new}}=\rho$。由于实践中 $\rho=O(10^{-1})$，取 $\alpha=O(10^{-1})$ 即足够，此时对角化的舍入放大 $\mathrm{Cond}_2(V)\le1/\alpha$ 可以忽略。**编号 39 对 MGRIT 是否给出同一个阈值本页未核实**，只确认了「适当选取 $\alpha$ 可保持收敛率」这一表述。

## 46：前向-后向系统没有单一传播方向

抛物 PDE 约束最优化的困难在摘要中说得很清楚：计算必须**一次性**考虑全部离散时间点，而待解的耦合系统具有**相反的演化方向**——状态方程从初值向前，伴随方程从终值向后。这是定义性的障碍。顺序时间推进对耦合系统根本不适用（无法同时向前与向后推进），parareal 式的补救也不适用，因为没有单一的传播方向可供迭代。

标准的分布控制问题是极小化

$$
\mathcal J(y,q)=\tfrac12\|y-y_d\|^2_{L^2(\Omega\times(0,T))}
+\tfrac{\gamma}{2}\|q\|^2_{L^2(\Omega\times(0,T))}
\quad\text{s.t.}\quad
\partial_t y-\Delta y=f+q,
$$

消去控制后一阶最优性系统是前向-后向对

$$
\begin{cases}
\partial_t y-\Delta y-\tfrac1\gamma p=f,\\
-\partial_t p-\Delta p+y=y_d,
\end{cases}
$$

配以初值条件（初值问题）或周期条件（时间周期问题）。

论文的组织对象是两个**时间离散矩阵**：$B_{\mathrm{per}}$（时间周期）与 $B_{\mathrm{ini}}$（初值）。摘要的表述是主要思想在于「小心处理相关的时间离散矩阵」。此前工作通过**近似离散 KKT 系统的 Schur 补**构造预条件子，而本文的出发点是这类预条件子的收敛远慢于对角化构造所能达到的水平。

编号 31 与编号 46 的关系可以一句话说清：**编号 31 是「$\alpha=1$ 是免费的」这一观察，编号 46 是「那就用一个周期问题去近似非周期问题」。** 前者对角化 $B_{\mathrm{per}}$ 本身，后者把周期型近似 $\widehat B_{\mathrm{per}}$ 用作初值问题矩阵 $B_{\mathrm{ini}}$ 的**预条件子**。

> [!note] 可核实范围
> 上述障碍表述、两个时间离散矩阵的组织地位、与 Schur 补预条件子的对比均可从摘要确认。该文的具体代价函数记号、预条件子的确切形式、定理与数值结果本站未核实；出版社全文对自动化访问返回 403。

## 53：把对角化用在级之间

编号 53 把同一机制用在一个不同的方向上：不是在时间层之间，而是在**单个时间步的级之间**。两级奇异对角隐式 Runge-Kutta 方法每步要解两个隐式级，而这两个级通常顺序求解。对角化把它们解耦，使两个级可以并发求解。

> [!note] 可核实范围
> 该文正文未能获取。可确认的是题名所示的对象（两级 SDIRK 方法的并行实现，途径是对角化）。其具体的对角化形式、条件数分析与并行效率结果本站未核实。

## 四篇按「价格」排列

| 编号 | 被对角化的对象         | $\alpha$ 或步长的角色      | 条件数代价                      |
| ---- | ---------------------- | -------------------------- | ------------------------------- |
| 31   | 时间周期的循环矩阵     | $\alpha=1$，精确周期       | $\mathrm{Cond}_2(V)=1$，无代价  |
| 39   | 首尾耦合后的粗校正矩阵 | $\alpha\in(0,1)$，人为引入 | $\mathrm{Cond}_2(V)\le1/\alpha$ |
| 46   | 前向-后向最优性系统    | 周期近似用作预条件子       | 通过外层迭代吸收                |
| 53   | 单步内的级             | 不涉及 $\alpha$            | 由方法的级结构决定              |

这张表把这条路线的逻辑显示得比较清楚：**对角化本身不是目的，让某个矩阵可对角化并且其特征向量矩阵条件良好才是。** 时间周期问题免费提供了这一点；其余情形必须用一个可控的近似去换，而近似的误差再交给外层 Krylov 迭代。后一条正是[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|全时间预条件]]一页的主题。

## 覆盖核对

| 内容                             | 论文 | 覆盖状态                                               |
| -------------------------------- | ---- | ------------------------------------------------------ |
| 全时间 Kronecker 系统与三步模式  | 31   | 系统、三步、FFT 与并行步                               |
| $\alpha$-循环矩阵与 $\alpha=1$   | 31   | 两个矩阵形式、酉性、条件数为 1                         |
| 复移位系统的多重网格与光滑子选择 | 31   | 阻尼 Richardson 的理由、与网格无关的因子               |
| MGRIT 收敛因子与 FCF 松弛        | 39   | 多出的因子、三种松弛的定义、与 parareal 的关系         |
| Lobatto IIIC 作为粗传播子        | 39   | 稳定函数、Padé 阶、L-稳定与阻尼性质                    |
| 首尾耦合与可对角化的粗校正       | 39   | 耦合条件、重定义的必要性、全时间系统                   |
| 前向-后向障碍与两个时间矩阵      | 46   | 相反演化方向、$B_{\mathrm{per}}$ 与 $B_{\mathrm{ini}}$ |
| 与编号 31 的逻辑关系             | 46   | 「免费的 $\alpha=1$」到「用周期近似非周期」            |

## 本页原文

- S. Wu, H. Zhang, and T. Zhou, [_Solving time-periodic fractional diffusion equations via diagonalization technique and multigrid_](https://doi.org/10.1002/nla.2178), Numer. Linear Algebra Appl. 25(5) (2018), e2178。
- S. Wu and T. Zhou, [_Acceleration of the two-level MGRIT algorithm via the diagonalization technique_](https://doi.org/10.1137/18M1207697), SIAM J. Sci. Comput. 41(5) (2019), pp. A3421-A3448。
- S. Wu and T. Zhou, [_Diagonalization-based parallel-in-time algorithms for parabolic PDE-constrained optimization problems_](https://doi.org/10.1051/cocv/2020012), ESAIM Control Optim. Calc. Var. 26 (2020), 88。
- S. Wu and T. Zhou, [_Parallel implementation for the two-stage SDIRK methods via diagonalization_](https://doi.org/10.1016/j.jcp.2020.110076), J. Comput. Phys. 428 (2021), 110076。
