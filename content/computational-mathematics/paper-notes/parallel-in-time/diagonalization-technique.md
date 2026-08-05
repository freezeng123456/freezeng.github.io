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
> 编号 **31**（_Numer. Linear Algebra Appl._ 25(5), 2018）、**39**（_SIAM J. Sci. Comput._ 41(5), 2019）、**46**（_ESAIM Control Optim. Calc. Var._ 26, 2020）、**53**（_J. Comput. Phys._ 428, 2021）。四篇均无预印本。编号 39 的四个收缩常数已逐式核对，可以给出完整推导与等代价折算；编号 31、46、53 的摘要可核实，因此问题设定、构造思路与定理的**形式**可以给出，但其中的具体常数（最优阻尼值、聚集半径、实测迭代数）均未核实。

![把串行时间递推换成可对角化的时间矩阵](assets/diagrams/tao-zhou-papers/zh/pint-diagonalization.svg)

## 共用装置：三步模式与它的价格

### 直觉：为什么可对角化的时间矩阵能把顺序扫描变成并发求解

时间推进之所以串行，是因为第 $n$ 层的解依赖第 $n-1$ 层。把所有时间层堆成一个大向量后，这种依赖表现为一个**块下三角**矩阵——顺序时间推进就是对它做块前代。前代无法并行，因为它天生是一条链。

对角化改变的是**看这条链的坐标系**。若时间方向的矩阵 $B$ 可写成 $VDV^{-1}$，那么在 $V^{-1}$ 给出的新坐标下，各时间自由度之间的耦合被完全解开：原来「第 $n$ 层等着第 $n-1$ 层」变成「第 $n$ 个**模态**独立求解」。每个模态对应一个复标量 $\lambda_n$，要解的是一个复移位的空间问题 $(\lambda_nI_x+\cdots)x=y$，而这 $N_t$ 个问题彼此无关，可以放在 $N_t$ 个进程上同时做。乘 $V^{-1}$ 与乘 $V$ 只在时间方向上作用，若 $V$ 是（缩放的）DFT 矩阵，这两步就是 FFT，代价 $O(N_xN_t\log N_t)$，相对于空间求解可以忽略。

**因此整个方法的成败落在一句话上：那个时间矩阵是否可对角化，以及 $V$ 的条件数是多少。** 后半句是实质。$V^{-1}$ 与 $V$ 是在浮点数下乘上去的，因此计算出的解带有 $O(\epsilon\,\mathrm{Cond}_2(V))$ 量级的舍入污染（$\epsilon$ 为机器精度）。一个理论上可对角化但 $V$ 条件数为 $10^{14}$ 的时间矩阵在双精度下没有任何用处。本页四篇论文可以按它们如何处理这个价格来排列。

### 全时间系统与 $\alpha$-循环对角化

把全部时间层的未知量堆成一个向量，得到 Kronecker 形式的全时间系统

$$
\bigl(B_1\otimes I_x+B_2\otimes A\bigr)\boldsymbol u=\boldsymbol b,
\qquad
\boldsymbol u=(U_1^\top,\dots,U_{N_t}^\top)^\top .
$$

ParaDiag 一族的组织原则可以一句话说清：**所有算法都只处理 $B_1$ 与 $B_2$，而让空间矩阵保持不动。** 处理的方式有两种，正好是这条文献的两半：把 $B:=B_2^{-1}B_1$ **精确**对角化（直接解法），或把 $B_1,B_2$ 换成 $\alpha$-**循环矩阵**从而得到一个可对角化的**预条件子**（迭代解法）。本页四篇都在后一支或它的边界上。

$\alpha$-循环矩阵的关键性质是**同时**可对角化。记 $\mathbb F$ 为单位 DFT 矩阵、$\omega=e^{2\pi\mathrm i/N_t}$，并对 $\alpha\in(0,1]$ 取

$$
\Gamma_\alpha:=\mathrm{diag}\bigl(1,\ \alpha^{1/N_t},\ \alpha^{2/N_t},\dots,\ \alpha^{(N_t-1)/N_t}\bigr),
$$

则**任意**两个 $\alpha$-循环矩阵被同一个特征向量矩阵对角化：

$$
C_j^{(\alpha)}=VD_jV^{-1},
\qquad
V=\Gamma_\alpha^{-1}\mathbb F^{*},
\qquad
D_j=\mathrm{diag}\bigl(\sqrt{N_t}\,\mathbb F\,\Gamma_\alpha\,C_j^{(\alpha)}(:,1)\bigr).
$$

注意 $D_j$ 只需要 $C_j^{(\alpha)}$ 的**第一列**——特征值可以由一次长度 $N_t$ 的 FFT 直接算出，不需要任何特征值求解器。求解因此变成三步：

$$
\text{(a) } S_1=(\mathbb F\Gamma_\alpha\otimes I_x)\boldsymbol b,
\qquad
\text{(b) } S_{2,n}=(\lambda_{1,n}I_x+\lambda_{2,n}A)^{-1}S_{1,n},
\qquad
\text{(c) } \boldsymbol u=(\Gamma_\alpha^{-1}\mathbb F^*\otimes I_x)S_2 ,
$$

步骤 (a) 与 (c) 是带缩放的 FFT 与其逆（$O(N_xN_t\log N_t)$），步骤 (b) 在 $N_t$ 个时间层上完全并行。$\alpha=1$ 时 $\Gamma_\alpha=I$、$V=\mathbb F^*$ 退化为普通循环矩阵的情形。

### 价格：一条不等式

$$
\mathrm{Cond}_2(V)=\mathrm{Cond}_2(\Gamma_\alpha^{-1}\mathbb F^{*})
\le\mathrm{Cond}_2(\Gamma_\alpha^{-1})\,\mathrm{Cond}_2(\mathbb F^{*})
=\mathrm{Cond}_2(\Gamma_\alpha^{-1})\le\frac1\alpha .
$$

酉矩阵 $\mathbb F^*$ 是免费的，全部条件数来自缩放因子 $\Gamma_\alpha^{-1}$，其对角元从 $1$ 张到 $\alpha^{-(N_t-1)/N_t}$，因此上界 $1/\alpha$。**这一条不等式是这条路线的主权衡**：$\alpha$ 越小，$\alpha$-循环矩阵越接近原来的下三角 Toeplitz 矩阵（近似越好、外层迭代收缩越快），而 $V$ 的条件数越差（舍入放大越严重）。

把两边写在一起就能看清楚该怎么选 $\alpha$。收缩侧用这条路线的标准界 $\rho\le\alpha/(1-\alpha)$（梯形规则情形由 Gander 与 Wu 2019 给出，对全部稳定一步法由[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|编号 65]] 证明），舍入侧用 $\epsilon\,\mathrm{Cond}_2(V)\le\epsilon/\alpha$（双精度 $\epsilon=2.22\times10^{-16}$）：

| $\alpha$   | $\rho\le\alpha/(1-\alpha)$ | $\mathrm{Cond}_2(V)\le1/\alpha$ | 精度地板 $\approx\epsilon/\alpha$ | 达到 $10^{-10}$ 所需迭代 |
| ---------- | -------------------------- | ------------------------------- | --------------------------------- | ------------------------ |
| $10^{-1}$  | $0.111$                    | $10$                            | $2.2\times10^{-15}$               | $11$                     |
| $10^{-2}$  | $0.0101$                   | $10^{2}$                        | $2.2\times10^{-14}$               | $6$                      |
| $10^{-3}$  | $0.00100$                  | $10^{3}$                        | $2.2\times10^{-13}$               | $4$                      |
| $10^{-13}$ | $10^{-13}$                 | $10^{13}$                       | $2.2\times10^{-3}$                | $1$                      |

（表中的收缩与条件数两列是上面两条**已核实**的界；精度地板与迭代次数是由它们直接算出的量级估计，迭代次数取 $\lceil 10/\log_{10}(1/\rho)\rceil$。）

最后一行说明为什么不能把 $\alpha$ 一味取小：$\alpha=10^{-13}$ 确实一轮就收敛，但收敛到的只是 $10^{-3}$ 量级的精度。ParaDiag 综述给出的实践建议正是 $\alpha=10^{-2}$ 与 $10^{-3}$，并明确警告不要取 $\alpha=10^{-13}$——上表把这条建议的理由摆在了一起。

### 同一个方法的三种写法

这条路线在不同论文中以三种语言出现，识别它们是同一个东西可以省去很多混淆。

1. **预条件的定常迭代。** $\mathcal P_\alpha\Delta\boldsymbol U^{k}=\boldsymbol b-\mathcal K\boldsymbol U^{k}$，$\boldsymbol U^{k+1}=\boldsymbol U^{k}+\Delta\boldsymbol U^{k}$。$\alpha=1$ 是 McDonald、Pestana 与 Wathen（2018）的离散构造；$\alpha\in(0,1)$ 是 Banjai 与 Peterseim（2012）的并行方法。
2. **Krylov 加速。** 直接对 $\mathcal P_\alpha^{-1}\mathcal K\boldsymbol U=\mathcal P_\alpha^{-1}\boldsymbol b$ 用 GMRES/MINRES/BiCGStab/CG。它就是上面那个定常迭代写在不动点处，而且**即使 $\rho(\mathcal P_\alpha^{-1}\mathcal K)\ge1$ 也可能有效**——这一点在非线性问题上会用到。
3. **首尾耦合的波形松弛**（Gander 与 Wu, _Numer. Math._ 143 (2019) 489-527），在连续层面：

   $$
   \boldsymbol u^k_t(t)=A\boldsymbol u^k(t)+\boldsymbol g(t),
   \qquad
   \boldsymbol u^k(0)=\alpha\bigl[\boldsymbol u^k(T)-\boldsymbol u^{k-1}(T)\bigr]+\boldsymbol u_0 .
   $$

   收敛时尾项相消，恢复原初值问题。把这个**周期型**问题离散，得到的恰好是 $\mathcal P_\alpha\boldsymbol U^k=\boldsymbol b^k$。

第三种写法给出了整条路线的概念内核：**在时间方向加一个 $\alpha$-周期边界条件，就是让时间矩阵成为循环矩阵、从而可 FFT 对角化的那一步。** 「把初值问题装扮成周期问题」这一个念头，贯穿本页的编号 31、39、46、53 与下一页的编号 65、71、84。

## 31：时间周期问题让价格归零

### 直觉

上面那张表描述的是一场需要付款的交易。编号 31 的观察是：**有一类问题，这笔款项恰好为零。**

时间周期问题对时间推进而言结构上是坏消息：没有初值可以出发，周期条件 $u(0)=u(T)$ 把最后一个时间层耦合回第一个，因此必须一次性考虑全部离散解。无论是否以并行为目标，都被迫进入全时间表述。但从对角化的角度看这是好消息：周期耦合使时间离散矩阵成为**真正的循环矩阵**，而循环矩阵由**酉**的 DFT 矩阵对角化。$\mathrm{Cond}_2(V)=1$，$\alpha$ 这个参数根本不需要出现，也没有近似需要外层迭代去修。于是得到一个**直接**（非迭代）的时间并行算法。

论文的第二半处理由此产生的新问题：解耦后每个时间模态留下一个**复系数**的空间问题，而空间算子是稠密的分数阶 Laplacian。直接分解不可行，必须用迭代解法，而复移位使标准多重网格理论不能照搬。

### 问题设定

含分数阶 Laplacian 的时间周期扩散方程。空间半离散给出稠密的 $A$；时间离散用 $\theta$-方法配周期条件。全时间系统即前一节的 Kronecker 形式，其中 $B_1,B_2$ 是真正的循环矩阵——用 $\alpha$-循环族的语言说，这对应特殊值 $\alpha=1$：

$$
C_1^{(\alpha)}=\frac{1}{\Delta t}
\begin{bmatrix}1&&&-\alpha\\-1&1&&\\&\ddots&\ddots&\\&&-1&1\end{bmatrix},
\qquad
C_2^{(\alpha)}=
\begin{bmatrix}\theta&&&(1-\theta)\alpha\\1-\theta&\theta&&\\&\ddots&\ddots&\\&&1-\theta&\theta\end{bmatrix}.
$$

### 推导

**第一层：直接对角化。** 在 $\alpha=1$ 处 $\Gamma_\alpha=I$、$V=\mathbb F^*$ 是酉矩阵，故 $\mathrm{Cond}_2(V)=1$。共用装置的三步模式照常适用，但三点不同：不需要选 $\alpha$；不需要外层 Krylov 迭代；不存在舍入放大。**限制其他对角化路线的整个权衡在这里消失。** 这就是时间周期情形允许一个直接时间并行算法的原因。

**第二层：复移位系统的多重网格。** 步骤 (b) 留下一列独立的复系数线性系统 $(\lambda_{1,n}I_x+\lambda_{2,n}A)x=y$。两个特征决定了求解器的选择：其一，$A$ 是稠密的分数阶 Laplacian，任何需要三角求解或稀疏结构的光滑子（阻尼 Jacobi、Gauss-Seidel）都不合适，而只需要矩阵向量乘的**阻尼 Richardson 迭代**正好合适；其二，移位是复的，因此不能直接引用实对称正定情形的多重网格收敛理论，收敛因子必须重新证明。论文的做法是证明该多重网格具有**与网格无关**的收敛因子。

**第三层：优化阻尼参数。** 阻尼 Richardson 的收敛因子是阻尼参数的函数，论文对该参数做优化，以**最小化这个常数收敛因子**。

### 定理

摘要可确认三条：（i）对角化给出**全部离散解的直接并行计算**，算法在时间方向非迭代；（ii）复移位空间系统的多重网格求解器具有**与网格无关的收敛因子**；（iii）存在一个使该常数因子最小的阻尼参数选择，论文把它算了出来。

> [!note] 可核实范围
> 上述三条的**存在**均可从摘要确认。**最优阻尼值、最小化后的收敛因子数值、分析是否关于分数阶阶数或 $N_t$ 一致，本站均未核实**，因此本页不给出任何具体常数。

### 数值实验

摘要只确认「给出了数值结果以支持理论分析」。**具体的测试方程、分数阶阶数、空间网格规模、$N_t$、实测多重网格收敛因子与最优阻尼值，本站均未核实。**

这里缺失的信息里有一项特别值得指出：本文的核心量化主张是那个**被最小化的常数收敛因子**，而它是一个纯数字。没有这个数字，读者无法判断「与网格无关」是与 $0.1$ 无关还是与 $0.9$ 无关——两者在实用上差别巨大。第一层的对角化不需要实验支持（它是精确的、非迭代的），第二、三层则完全依赖它。

### 与其他论文的关系

这是本系列的**转折点**：第一篇用对角化技术而非 parareal 迭代的工作，之后的每一篇都建立在它之上。编号 46 是它的决定性推广——用同一机制对角化时间周期最优控制问题的 $B_{\rm per}$，再迈出关键一步，把周期型近似 $\widehat B_{\rm per}$ 用作**初值**问题的预条件子。编号 39、53 与下一页的编号 65、71、84 都使用 $\alpha\in(0,1)$ 的 $\alpha$-循环机制，即代价为 $\mathrm{Cond}_2(V)\le1/\alpha$、回报为 $O(\alpha)$ 收缩的那一支。[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|编号 59]] 则攻击对角化的**另一个**障碍——Maday-Rønquist 几何步长路线中 $V$ 的病态——办法是选一个使 $\mathrm{Cond}_2(V)=O(n^2)$ 的时间积分器。

## 39：同时修好粗传播子与粗校正

### 直觉

MGRIT 是 parareal 的多层推广，也是 XBraid 软件包背后的算法。它比 parareal 多做一次 F-松弛，收缩因子相应多一个因子。但这里有一个容易被忽略的记账问题：**多出来的那次 F-松弛也要花一次细求解**，所以「MGRIT 比 parareal 收缩得快」这件事本身不构成收益。

论文的出发点是两个具体的低效，而且它们互相牵制。第一，粗传播子几乎总取后向 Euler，它只有一阶；换一个衰减性质更好的粗传播子应当能显著改进收缩，但没有人问过。第二，粗网格校正是沿 $N_t$ 个粗时间层的**顺序**扫描，处理器数增加后它按 Amdahl 定律成为串行瓶颈。**这两者是冲突的**：更好的粗传播子通常是多级隐式 Runge-Kutta 方法，代价是后向 Euler 步的若干倍，因此改善第一处会加重第二处。

论文的做法是同时修两处，并说明它们相容：把粗传播子换成二阶 Lobatto IIIC，同时用首尾耦合把顺序粗校正变成一个可 FFT 对角化的全时间系统。第一处改动的额外代价，正好由第二处改动带来的并行性吸收掉。这也是为什么两个改动分开看都不如合起来明显。

### 问题设定

线性 ODE 系统，系数矩阵对称正定。两级 MGRIT 配 FCF-松弛的迭代式为

$$
\begin{aligned}
\boldsymbol u_0^{k+1}&=\boldsymbol u_0,\qquad \boldsymbol u_1^{k+1}=\mathcal F(T_0,T_1,\boldsymbol u_0),\\
\boldsymbol u_{n+1}^{k+1}&=\mathcal F\bigl(T_n,T_{n+1},\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^{k})\bigr)
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G\bigl(T_n,T_{n+1},\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^{k})\bigr),
\end{aligned}
$$

$n=1,\dots,N_t-1$。术语上：**F-松弛**是在每个粗区间内部的 F 点上应用细传播子，从该区间左端的 C 点值出发，这些求解相互独立、可并行；**C-松弛**是用最后一个 F 点更新 C 点值；**FCF-松弛**是 F、C、F 依次进行，因此每轮代价是两次细求解。把上式与 parareal 迭代对照可以看出：**MGRIT-FCF 就是把 parareal 中每处 $\boldsymbol u_n^k$ 换成 $\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)$**，即带一个粗区间重叠的 parareal。仅用 F-松弛的两级 MGRIT **就是** parareal。

相应的线性收敛因子为

$$
\varrho_l(J,z)=\frac{\bigl|R_f^{J}(z/J)\bigr|\,\bigl|R_g(z)-R_f^{J}(z/J)\bigr|}{1-|R_g(z)|}
\qquad\Longleftrightarrow\qquad
\varrho_{l,\mathrm{MGRIT}}=\bigl|R_f^J(z/J)\bigr|\times\varrho_{l,\mathrm{parareal}} .
$$

多出的 $|R_f^J(z/J)|\le1$ 正是额外一次 F-松弛带来的全部收益：**多花的那次细求解买到的恰好是多作用一次细传播子的稳定函数。**

### 推导

**策略一：换粗传播子。** 取二阶 Lobatto IIIC 方法作 $\mathcal G$，其 Butcher 表与稳定函数为

$$
\begin{array}{c|cc}
0 & \tfrac12 & -\tfrac12\\
1 & \tfrac12 & \tfrac12\\ \hline
 & \tfrac12 & \tfrac12
\end{array},
\qquad
R_g(z)=\frac{1}{1+z+z^2/2}
$$

（后者按 $\boldsymbol u'+A\boldsymbol u=g$ 的约定写出）。这是 $(0,2)$ 阶 Padé 逼近：A-稳定、**L-稳定**（$R_g(\infty)=0$）、二阶。关键不在阶数而在**衰减速度**：与后向 Euler 的 $1/(1+z)$ 相比，$1/(1+z+z^2/2)$ 在高频上按 $z^{-2}$ 而不是 $z^{-1}$ 衰减。回到 $\varrho_l$ 的读法——分子是粗细失配、分母是粗传播的耗散裕量——更强的高频阻尼同时压小分子（$R_g$ 更接近同样衰减的 $R_f^J$）并抬高分母（$1-|R_g|$ 更接近 $1$）。两个效应同向，这就是常数从 $0.2984$ 降到 $0.0817$ 的机制。

**策略二：把粗校正并行化。** 用一个**首尾耦合条件**替换顺序粗校正：

$$
\boldsymbol u_{n+1}^{k+1}=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\underbrace{\mathcal F(T_n,T_{n+1},\tilde{\boldsymbol u}_n^{k})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k})}_{=:\,\boldsymbol b_{n+1}^k},
\qquad
\boldsymbol u_0^{k+1}=\alpha\,\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0 ,
$$

其中 $\tilde{\boldsymbol u}_n^k=\boldsymbol u_n^k$（$n\ge1$）而 $\tilde{\boldsymbol u}_0^k=\boldsymbol u_0$——这一重定义是必要的，否则不动点不再是真解。取 $\mathcal G$ 为后向 Euler 并写 $\boldsymbol u'=A\boldsymbol u$ 时，耦合后的粗校正就是全时间系统

$$
\bigl(C_\alpha\otimes I_x-I_t\otimes\Delta T A\bigr)\boldsymbol U^{k+1}=\boldsymbol g^k ,
\qquad
C_\alpha=\begin{bmatrix}1&&&-\alpha\\-1&1&&\\&\ddots&\ddots&\\&&-1&1\end{bmatrix},
$$

即用一个 $\alpha$-循环矩阵替换原来的下二对角 Toeplitz 矩阵，于是共用装置的三步模式适用：FFT、$N_t$ 个并行的复移位空间求解、逆 FFT。注意 $\alpha\to0$ 恢复标准的顺序粗校正，$\mathrm{Cond}_2(V)\le1/\alpha$ 是舍入代价。这正是上一节第三种写法的实例：$\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0$ 就是「把问题装扮成 $\alpha$-周期问题」。

**策略二之补：让两级方法的代价降到一级。** 摘要明确宣称，在并行粗校正的框架内，隐式两级 Runge-Kutta 方法 LIIIC-2 的代价可以**降到后向 Euler 的水平**，途径仍然是对角化的一个恰当应用。

> [!note] 这一条的机制未核实
> 一个与该宣称相符的自然机制是：LIIIC-2 的 Butcher 矩阵 $A_{\rm RK}=\bigl[\begin{smallmatrix}1/2&-1/2\\1/2&1/2\end{smallmatrix}\bigr]$ 的特征值是共轭对 $\tfrac12(1\pm\mathrm i)$，对角化它使两个隐式级解耦成两个独立的复移位空间求解，而由于两个移位互为共轭，实际只需解**一个**复系统，另一个取共轭即可，因此每步代价与一次实的后向 Euler 求解相当。**这是与宣称相符的重建，不是论文的机制**；论文如何实现这一点本站未核实。

### 定理

**（稳健性。）** 论文证明 MGRIT 的收敛因子是**稳健**的：**与系数矩阵的特征值无关，也与粗细比 $J=\Delta T/\Delta t$ 无关**。

**（四个可直接比较的常数。）** 设 $\mathcal F$ 为 L-稳定积分器、$J=O(1)$，则 $\max_{z\ge0}\varrho_l$ 为

| 粗传播子 $\mathcal G$ | Parareal | MGRIT（FCF-松弛） |
| --------------------- | -------- | ----------------- |
| 后向 Euler            | 0.2984   | 0.1115            |
| Lobatto IIIC（二阶）  | 0.0817   | 0.0197            |

这就是摘要中「从 0.1 降到 0.02」的精确形式。

**（参数 $\alpha$ 的取法。）** 适当选取 $\alpha$ 时，带并行粗校正的新算法与原算法**收敛率相同**。parareal 情形有一个明确阈值（归于 Wu, SISC 2018）：只要

$$
\alpha\le\frac{\rho}{1+\rho},
$$

就有 $\rho_{\text{new}}=\rho$。由于实践中 $\rho=O(10^{-1})$，取 $\alpha=O(10^{-1})$ 即足够，此时对角化的舍入放大 $\mathrm{Cond}_2(V)\le1/\alpha$ 可以忽略。**编号 39 对 MGRIT 是否给出同一个阈值本站未核实**，只确认了「适当选取 $\alpha$ 可保持收敛率」这一表述。

> [!warning] 等代价折算：一次 MGRIT-FCF 迭代略差于两次 parareal 迭代
> 上表的两列**不能直接横向比较**，因为它们的每轮代价不同。FCF-松弛每轮做两次细求解，所以 MGRIT 的比较对象是两轮 parareal，而
>
> $$
> 0.2984^2=0.0890<0.1115,
> \qquad
> 0.0817^2=0.0067<0.0197 .
> $$
>
> 两次 parareal 都比一次 MGRIT-FCF 更小。**真正稳固的收益来自换粗传播子，而不是来自换松弛方式。**

把四个常数折算成「每次细求解的等效收缩」（最后一列由 $\varrho_l$ 按每轮细求解次数开方得到），这一点看得更清楚：

| 粗传播子   | 方法      | $\varrho_l$ | 每轮细求解 | 等效单次细求解收缩 |
| ---------- | --------- | ----------- | ---------- | ------------------ |
| 后向 Euler | parareal  | 0.2984      | 1          | 0.2984             |
| 后向 Euler | MGRIT-FCF | 0.1115      | 2          | 0.3339             |
| LIIIC-2    | parareal  | 0.0817      | 1          | 0.0817             |
| LIIIC-2    | MGRIT-FCF | 0.0197      | 2          | 0.1404             |

两个 parareal 行都优于对应的 MGRIT 行。相反，同一列内从后向 Euler 换到 Lobatto IIIC，收缩因子降低约 $3.7$ 倍（parareal）与 $5.7$ 倍（MGRIT-FCF），而粗传播子的代价在并行粗校正框架下被压回后向 Euler 的量级——这才是本文净收益的来源。

### 数值实验

摘要确认的两个测试问题各有明确的检验目的：

- **带不确定系数的对流扩散方程。** 这是一族参数化的 ODE 系统（UQ 场景），每个实现都要长时间积分，正是时间并行能兑现收益的地方；同时它检验上面那条稳健性主张——不同实现给出不同的 $\sigma(A)$，而理论宣称的界不依赖它。
- **Gray-Scott 模型。** 刚性非线性反应扩散系统，具有斑图生成行为。它把算法推到线性理论的适用范围之外，检验四个常数在非线性问题上是否仍有预测力。

定性结论是实验「支持本文的结论」，即实测因子与 $\approx0.1$（后向 Euler 粗传播子）和 $\approx0.02$（LIIIC-2 粗传播子）一致，且并行粗校正不使收敛率退化。**实测的加速比、进程数与迭代次数本站未核实。**

作为对照，可以引用一组**独立于本文**的、可核实的数据，说明这套常数在什么条件下开始失效。综述对热方程与对流扩散方程做的黏性扫描给出：

| 黏性 $\nu$ | parareal 的 $\varrho_{l,\max}$ | MGRIT-FCF 的 $\varrho_{l,\max}$     | 状态                       |
| ---------- | ------------------------------ | ----------------------------------- | -------------------------- |
| $0.1$      | —                              | $\approx\varrho_{l,\rm parareal}^2$ | 平方关系成立               |
| $0.01$     | 接近 $1$                       | 接近 $1$                            | 两者都濒临失效             |
| $0.002$    | $1.4211$                       | $1.2812$                            | 两者都发散，MGRIT 退化较轻 |

综述对最后一行的解释是：MGRIT 连做两次细求解，中间没有插入一次（此时已无帮助的）粗求解，因此退化较轻。

> [!note] 这组数据的出处
> 上表来自 Gander、Wu 与 Zhou 的 _Acta Numerica_ 综述（编号 85），**不是**编号 39 自己的实验。本页引用它是为了给「$\varrho_{l,\rm MGRIT}\approx\varrho_{l,\rm parareal}^2$」这条关系一个可核实的适用边界，而不是把这些数字归给编号 39。

**这些实验建立了什么、又缺什么。** 能建立的是：四个常数在耗散问题上有预测力，并行粗校正不损害收敛率，且方法在非线性问题上仍然可用。缺的是本文最核心主张的直接证据——**并行粗校正究竟消除了多少 Amdahl 瓶颈**。这需要在固定问题上做处理器数扫描，报告强可扩展曲线；没有可核实的这类数据，本站只能确认收敛率不退化，不能确认并行效率的改善幅度。

### 与其他论文的关系

[[computational-mathematics/paper-notes/parallel-in-time/parareal-convergence|编号 12]] 提供了分析装置（$\varrho_l$、稳定函数、在 $z$ 与 $J$ 上的稳健性）与 parareal 的参考常数 $0.2984$；编号 39 是它的 MGRIT 对应物，也是把四个常数并排钉死的那一篇。编号 31 提供了对角化思想，而编号 39 是本列表中第一篇把它用在**迭代方法内部**（作为并行化粗校正的装置）而不是作为独立直接解法的工作——这里的 $C_\alpha$ 就是后来在编号 53、[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|65、71、84]] 中成为预条件子 $\mathcal P_\alpha$ 的同一个对象。编号 46 把同样的「用循环结构替换顺序结构再对角化」的动作带进最优控制设定，那里的顺序结构是一对前向-后向演化而不是单个前向扫描。编号 65 是最终的推广：不再逐情形分析，而对全部稳定一步法统一给出 $\mathcal P_\alpha^{-1}\mathcal K$ 的谱分析。[[computational-mathematics/paper-notes/parallel-in-time/parareal-convergence|编号 77]] 则放松使 $R_f^J(z/J)$ 有意义的均匀细网格假设。

## 46：前向-后向系统没有单一传播方向

### 直觉

前面所有方法都在处理**一个方向**的演化：信息从 $t=0$ 流向 $t=T$，全时间矩阵因此是块下三角的，而各种时间并行方法都是在设法避免对它做前代。最优控制把这个前提彻底取消：状态方程从初值向前，伴随方程从终值向后，两者逐点耦合。

**这里的关键区别是：时间并行不再是加速手段，而是必需品。** 顺序时间推进对耦合系统根本不适用——你无法同时向前和向后推进——parareal 式的补救也不适用，因为没有单一的传播方向可供迭代。所以问题不是「怎样让它更快」，而是「怎样解它」。

编号 46 的观察建立在编号 31 之上，但多走了决定性的一步。编号 31 说：时间周期问题的时间矩阵是循环的，因此免费可对角化。编号 46 先把这个观察用到时间周期最优控制问题上（得到一个直接算法），然后问：**既然周期问题这么好，为什么不用一个周期问题去近似非周期问题？** 初值问题的时间矩阵是不可对角化的 Jordan 块，但只要它出现在**预条件子**里而不是原系统里，就可以换成一个邻近的、可对角化的周期型矩阵，误差交给外层 Krylov 迭代。这一步是整个 ParaDiag-II 家族的概念核心。

### 问题设定

抛物 PDE 约束的分布控制问题：极小化

$$
\mathcal J(y,q)=\tfrac12\|y-y_d\|^2_{L^2(\Omega\times(0,T))}
+\tfrac{\gamma}{2}\|q\|^2_{L^2(\Omega\times(0,T))}
\quad\text{s.t.}\quad
\partial_t y-\Delta y=f+q,
$$

配 $y=0$ 于 $\partial\Omega\times(0,T)$，以及两种端点条件之一：$y(\cdot,0)=y(\cdot,T)$（**时间周期**模型问题）或 $y(\cdot,0)=y_0$（**初值**模型问题）。经由 $\gamma q=p$ 消去控制，一阶最优性（KKT）系统是前向-后向对

$$
\begin{cases}
\partial_t y-\Delta y-\tfrac1\gamma p=f, & y(\cdot,0)=y_0\ \text{或}\ y(\cdot,0)=y(\cdot,T),\\[2pt]
-\partial_t p-\Delta p+y=y_d, & p(\cdot,T)=0\ \text{或}\ p(\cdot,0)=p(\cdot,T).
\end{cases}
$$

时空离散后，记 $A$ 为离散负 Laplacian、$B\in\{B_{\rm per},B_{\rm ini}\}$ 为时间离散矩阵，全时间 KKT 系统具有形式

$$
\begin{bmatrix}
B\otimes I_x+I_t\otimes A & -\tfrac{1}{\gamma}I_t\otimes I_x\\[3pt]
I_t\otimes I_x & B^\top\otimes I_x+I_t\otimes A
\end{bmatrix}
\begin{bmatrix}\boldsymbol y\\ \boldsymbol p\end{bmatrix}
=\begin{bmatrix}\boldsymbol f\\ \boldsymbol y_d\end{bmatrix}.
$$

$(2,2)$ 块中的**转置** $B^\top$ 是后向伴随方程在代数上的签名：$B$ 下三角（因果、向前），$B^\top$ 上三角（反因果、向后）。（这个 Kronecker 形式是按摘要描述重建的标准写法，不是从论文读到的排版。）

摘要把论文的组织对象说得很清楚：主要思想在于**小心处理相关的时间离散矩阵** $B_{\rm per}$ 与 $B_{\rm ini}$。

### 推导

**决定性的结构事实。** 取等步长与后向 Euler，两个时间矩阵分别是

$$
B_{\rm ini}=\frac{1}{\Delta t}\begin{bmatrix}1&&&\\-1&1&&\\&\ddots&\ddots&\\&&-1&1\end{bmatrix},
\qquad
B_{\rm per}=\frac{1}{\Delta t}\begin{bmatrix}1&&&-1\\-1&1&&\\&\ddots&\ddots&\\&&-1&1\end{bmatrix}.
$$

$B_{\rm ini}$ 是下二对角 Toeplitz 矩阵，其唯一特征值是 $1/\Delta t$、代数重数 $N_t$；而 $B_{\rm ini}-\tfrac1{\Delta t}I$ 是严格下二对角矩阵，秩为 $N_t-1$，故几何重数只有 $1$。**它是一个单个 Jordan 块，亏损，不可对角化。** 相反，$B_{\rm per}$ 是循环矩阵，由酉的 DFT 矩阵对角化，$\mathrm{Cond}_2(V)=1$。论文利用的正是这一对立。

**算法一：时间周期问题的直接时间并行算法。** 直接对角化 $B_{\rm per}=VDV^{-1}$（$V$ 为 DFT 矩阵）。KKT 系统随之在时间方向解耦为 $N_t$ 个独立的 $2N_x\times2N_x$ 复系统，每个对应一个时间 Fourier 模态、把该模态上的状态与伴随耦合在一起，全部可并行求解。**没有迭代，没有收敛因子，没有 $\alpha$ 参数，也没有舍入代价。**

**算法二：用周期型代理预条件初值问题。** 由于 $B_{\rm ini}$ 亏损、无法对角化，改为在**预条件子**中把它替换为一个邻近的周期型（循环 / $\alpha$-循环）矩阵 $\widehat B_{\rm per}$，后者可由 FFT 对角化，于是算法一的机器可以直接用来施加预条件；替换带来的误差由外层 Krylov 方法吸收。

> [!note] 未核实的部分
> $\widehat B_{\rm per}$ 的确切定义——特别是它是普通循环矩阵（$\alpha=1$）还是 $\alpha\in(0,1)$ 的 $\alpha$-循环矩阵，以及角落元素如何选取——**本站未核实**。上面的 Jordan 块论证与 $2N_x$ 复系统的计数是从已核实的结构直接推出的，但论文的记号与实现细节未读到。

### 定理

**（直接性。）** 对时间周期问题，算法是**直接的**（非迭代）且在全部时间步上并行。

**（聚集性。）** 对初值问题，**对后向 Euler 与梯形规则两者**，预条件矩阵的**特征值与奇异值**的聚集性都被证明。这一对陈述里有两处值得强调：

- 同时证明**奇异值**聚集比只证特征值聚集更强，也更有用。预条件后的矩阵一般是非正规的，而对非正规矩阵，特征值聚集**不能**控制 GMRES 的收敛；奇异值信息才是真正给出残量界的东西。
- 覆盖**梯形规则**（A-稳定但不 L-稳定）而不只是后向 Euler，正好对应驱动[[computational-mathematics/paper-notes/parallel-in-time/parareal-convergence|编号 12]] 的那条 L-稳定 / 仅 A-稳定的分界。同一条分界在迭代路线上表现为收缩因子的门槛，在预条件路线上表现为聚集性是否成立。

**（与 Schur 补预条件子的比较。）** 与通过**近似离散 KKT 系统的 Schur 补**构造的现有预条件子相比，新预条件子对某些 Krylov 子空间解法——摘要点名 **GMRES 与 BiCGStab**——给出**快得多的收敛**。

> [!note] 未核实的部分
> **显式的聚集区间、聚集半径对 $\gamma$、$\Delta t$、$N_t$ 或 $\alpha$ 的依赖，以及定理编号，本站均未核实，因此不引用任何界。** 出版社全文对自动化访问返回 403。

### 数值实验

摘要确认给出了数值结果以说明所提时间并行算法的优势。从摘要的表述可以确定实验的**骨架**：两个模型问题（时间周期与初值抛物控制），基线是基于 Schur 补近似的预条件子，外层解法是 GMRES 与 BiCGStab。**具体的 PDE、正则化参数 $\gamma$ 的取值、网格规模、迭代次数与计时，本站均未核实。**

外层解法的选择本身是一条可核实且有信息量的线索：用 GMRES 与 BiCGStab 而不是 CG，意味着被迭代的预条件系统**不是对称正定的**。这与[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|编号 71]] 形成明确对照，后者用 CG，因此其预条件系统是对称正定的——这也与编号 71 被认定为「匹配 Schur 补预条件子的并行版本」相符。**同一类物理问题上两篇论文选用不同的外层解法，说明它们预条件的对象在代数结构上不同**，这一点比任何一个未核实的迭代次数都更有判别力。

### 与其他论文的关系

编号 31 是直接祖先：它发现时间周期性使时间矩阵成为循环矩阵，从而免费提供一个完美条件的对角化；编号 46 把这一观察带进最优控制设定，并额外迈出「用周期结构去**近似**非周期问题」的一步。编号 39 在低一层做同样的动作——它的首尾耦合条件 $\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0$ 字面上就是「让问题变成 $\alpha$-周期的，于是时间矩阵成为 $\alpha$-循环矩阵」；编号 46 把它用在一对耦合的前向-后向演化而不是单个前向扫描上。[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|编号 71]] 是成熟的后继：它处理一般的前向-后向演化方程，明确给出预条件子「把 Toeplitz 矩阵换成 $\alpha$-循环矩阵」，并对**任意**稳定一步法证明一个与网格无关的谱区间。编号 46 对两个具体积分器证明聚集性，编号 71 对全体积分器证明区间——正如编号 65 之于编号 53。

## 53：多级方法没有现成的全时间形式

### 直觉

$\alpha$-循环预条件的配方看上去只有一行：把全时间矩阵中的 Toeplitz 块换成 $\alpha$-循环块。但这条配方有一个隐含前提——**必须先有一个全时间矩阵**。

对线性多步法，这个前提自动成立：每个时间步的差分方程只涉及若干个时间层的解，把它们堆起来就是块 Toeplitz 矩阵。**多级 Runge-Kutta 方法不满足这个前提。** 一个 Runge-Kutta 步耦合的是若干个**级值**，而级值本身不是时间层未知量；它们是每步内部的中间量，通常在步内被消去。摘要把这一点说得很直白：对多级积分器「我们不能直接把差分方程写成全时间系统」。因此这里没有 Toeplitz 矩阵可以循环化，配方无从下手。

论文的做法是把全时间形式**重建**出来——保留级值作为未知量，重新组织成一个块结构，再在这个新结构上构造 $\alpha$-循环预条件子。它是本列表中第一篇对多级方法这么做的工作。它的结论有一个漂亮的形式：使收敛速率稳健所需的代数假设，恰好就是该方法的经典 A-稳定条件。

### 问题设定

两级奇异对角隐式 Runge-Kutta（SDIRK）方法。这条线索中常用的两张 Butcher 表是

$$
\underbrace{\begin{array}{c|cc}\gamma&\gamma&0\\ 1-\gamma&1-\gamma&\gamma\\ \hline &1-\gamma&1-\gamma\end{array}}_{\text{SDIRK22},\ \gamma=\frac{2-\sqrt2}{2}}
\qquad
\underbrace{\begin{array}{c|cc}\gamma&\gamma&0\\ 1-\gamma&\frac{-1}{\sqrt3}&\gamma\\ \hline &\frac12&\frac12\end{array}}_{\text{SDIRK23},\ \gamma=\frac{3+\sqrt3}{6}}
$$

「奇异对角隐式」指两个对角元相等，因此两个级求解使用**同一个**移位矩阵 $I+\gamma\Delta tA$，一次分解可服务两级。这个 $\gamma$ 正是摘要中的「principle element」。（这两张表出自作者自己的综述；编号 53 是否逐字采用它们本站未核实。）

### 推导

可核实的构造只有一句话，但它是关键的一句：本文构造的预条件子 $\mathcal P_\alpha$ **也是块 $\alpha$-循环矩阵，但结构与实现细节完全不同**。可确认的是：它在块层面是 $\alpha$-循环的；$\alpha\in(0,1)$；它通过块 Fourier 对角化施加，从而时间步解耦。作为对照，被推广的线性多步法机制在本页开头的共用装置一节已完整给出（同时对角化、$D_j$ 只需第一列、$\mathrm{Cond}_2(V)\le1/\alpha$、参考收缩 $\rho\le\alpha/(1-\alpha)$）。

> [!note] 未核实的部分
> **确切的块布局——特别是两个 SDIRK 级变量相对于时间下标如何排序、生成块是什么——本站未核实。** 摘要只保证它与线性多步法情形「结构完全不同」，不告诉我们那个结构是什么。

### 定理

**（主定理。）** 对 $\alpha\in(0,1)$，迭代矩阵的谱半径满足

$$
\rho(\text{迭代矩阵})\le\frac{\alpha}{1-\alpha},
\qquad\text{只要}\qquad \gamma\ge\tfrac14 ,
$$

其中 $\gamma$ 是两级 SDIRK 方法的主对角元。这与线性多步法理论中梯形规则情形的界**完全相同**，因此多级推广在速率上没有任何损失。

**（punchline。）** 条件 $\gamma\ge1/4$ **恰好就是该 SDIRK 方法的 A-稳定条件**。也就是说，得到稳健时间并行收敛率所需的代数假设，不是一条额外的人为限制，而是一条经典的稳定性条件。上面两张表都满足它：$(2-\sqrt2)/2\approx0.2929\ge0.25$，$(3+\sqrt3)/6\approx0.7887\ge0.25$。

**（作为研究纲领的猜想。）** 论文把这条结果表述为「对我们的猜想的一个初步验证：**隐式 Runge-Kutta 方法的 A-稳定条件足以保证预条件时间并行迭代算法具有稳健的 $O(\alpha)$ 收敛率**」。这条猜想是后续[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|编号 65 与 71]] 明确推进的研究纲领：编号 65 对全部稳定一步法证明了它，编号 71 在前向-后向耦合情形下确认了同一模式（代价是 $\alpha$ 必须随 $N_t$ 变化）。

**定理编号、界是否紧、以及 $\gamma<1/4$ 时的行为，本站均未核实。**

### 数值实验

摘要确认的两个测试问题都选得有针对性：

- **对流占优的扩散方程。** 这对时间并行是困难情形。parareal 一类方法之所以有效，靠的是耗散主导使粗传播子与细传播子在高频上一致；对流占优正好破坏这一点（参见[[computational-mathematics/paper-notes/parallel-in-time/parareal-convergence|parareal 收敛分析]]页中黏性减小到 $10^{-3}$ 量级时 parareal 发散的数据）。而 $\alpha$-循环预条件的界 $\rho\le\alpha/(1-\alpha)$ **不依赖 $\sigma(A)$**，因此在这类问题上应当不退化——这正是该算例要检验的。
- **黏性 Burgers 方程。** 非线性，检验线性理论之外的适用性。

两者都「支持理论结论」。**网格规模、$\alpha$ 的取值、Péclet 数、迭代次数与加速比，本站均未核实。**

这里未核实的量中最重要的是 $\alpha$ 的取值：本页开头那张权衡表说明，收敛速率与可达精度完全由 $\alpha$ 决定，因此一组不报告 $\alpha$ 的迭代次数是没有意义的。可以确认的是理论侧的形状——若取 $\alpha=10^{-2}$，则界给出 $\rho\le0.0101$，约 $6$ 轮到 $10^{-10}$，精度地板约 $2\times10^{-14}$。**这是把已核实的界代入算术的结果，不是论文报告的测量值。**

### 与其他论文的关系

编号 39 与编号 31 提供了 $\alpha$-循环与对角化的工具箱；编号 53 是本列表中第一篇把它用在**多级**积分器上的工作，为此不得不从头重建全时间形式。[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|编号 65]] 是单步方法方向上的自然下一步：不再一次一个积分器，而是对**全部**稳定单步法（一阶问题）与一大类对称两步法（二阶问题）统一给出 $\mathcal P_\alpha^{-1}\mathcal K$ 的谱分析。编号 53 与编号 65 因此是「停止逐情形研究」这一纲领的两半：前者是最后一个、也是最难的一个逐情形研究，后者是那条一般定理。[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|编号 71]] 把同一个 $\mathcal P_\alpha$-预条件-$\mathcal K$ 的思路从单个前向演化扩展到耦合的前向-后向对；[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|编号 84]] 则走向离散类型的另一个极端——时间谱方法，其全时间时间块是**稠密**的，而不是像这里一样被级耦合搅乱的稀疏块。

## 四篇按「价格」排列

| 编号 | 被对角化的对象         | $\alpha$ 的角色                              | 条件数代价                      | 换来的结论                                   |
| ---- | ---------------------- | -------------------------------------------- | ------------------------------- | -------------------------------------------- |
| 31   | 时间周期的循环矩阵     | $\alpha=1$，精确周期                         | $\mathrm{Cond}_2(V)=1$，无代价  | 直接（非迭代）算法                           |
| 39   | 首尾耦合后的粗校正矩阵 | $\alpha\in(0,1)$，人为引入                   | $\mathrm{Cond}_2(V)\le1/\alpha$ | $\alpha\le\rho/(1+\rho)$ 时保持收敛率        |
| 46   | 前向-后向最优性系统    | 周期近似 $\widehat B_{\rm per}$ 用作预条件子 | 通过外层 Krylov 吸收            | 特征值与奇异值聚集                           |
| 53   | 重建后的多级全时间系统 | $\alpha\in(0,1)$                             | $\mathrm{Cond}_2(V)\le1/\alpha$ | $\gamma\ge1/4$ 时 $\rho\le\alpha/(1-\alpha)$ |

这张表把这条路线的逻辑显示得比较清楚：**对角化本身不是目的，让某个矩阵可对角化并且其特征向量矩阵条件良好才是。** 时间周期问题免费提供了这一点；其余情形必须用一个可控的近似去换，而近似的误差再交给外层 Krylov 迭代。后一条正是[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|全时间预条件]]一页的主题。

## 覆盖核对

| 内容                               | 论文 | 覆盖状态                                                           |
| ---------------------------------- | ---- | ------------------------------------------------------------------ |
| 对角化为何能把顺序扫描变成并发求解 | —    | 坐标变换读法、复移位空间问题、FFT 代价                             |
| $\alpha$-循环矩阵的同时对角化      | —    | $V=\Gamma_\alpha^{-1}\mathbb F^*$、$D_j$ 只需第一列                |
| 主权衡与它的量化                   | —    | $\mathrm{Cond}_2(V)\le1/\alpha$ 的证明、四行权衡表                 |
| 同一方法的三种写法                 | —    | 定常迭代、Krylov、首尾耦合波形松弛                                 |
| 全时间 Kronecker 系统与三步模式    | 31   | 系统、三步、FFT 与并行步                                           |
| $\alpha=1$ 使价格归零              | 31   | 循环性、酉性、$\mathrm{Cond}_2(V)=1$、无参数权衡                   |
| 复移位系统的多重网格与光滑子选择   | 31   | 稠密 $A$ 与阻尼 Richardson、与网格无关的因子、阻尼优化（限定核实） |
| MGRIT-FCF 与 parareal 的精确关系   | 39   | 迭代式、三种松弛、重叠读法、多出的因子                             |
| Lobatto IIIC 为何改善收缩          | 39   | Butcher 表、$(0,2)$ Padé、$z^{-2}$ 衰减对分子分母的双向作用        |
| 首尾耦合与可对角化的粗校正         | 39   | 耦合条件、重定义的必要性、$C_\alpha$、$\alpha\to0$ 的极限          |
| 四个常数与等代价折算               | 39   | 常数表、平方比较、等效单次细求解收缩表                             |
| $\alpha$ 阈值与其归属              | 39   | $\alpha\le\rho/(1+\rho)$（parareal，Wu 2018）、MGRIT 情形未核实    |
| 两个测试问题各检验什么             | 39   | UQ 对流扩散、Gray-Scott；黏性扫描（出自综述）                      |
| 前向-后向障碍与两个时间矩阵        | 46   | KKT 系统、$B^\top$ 的含义、Jordan 块与循环矩阵的对立               |
| 两个算法与聚集性定理               | 46   | 直接算法、周期型代理、特征值与奇异值聚集、为何奇异值更强           |
| 外层解法的判别力                   | 46   | GMRES/BiCGStab 对 CG，与编号 71 的对照                             |
| 多级方法为何没有现成全时间形式     | 53   | 级值不是时间层未知量、配方失去作用对象                             |
| A-稳定条件与收敛率的重合           | 53   | $\gamma\ge1/4$、两张表的验算、$\rho\le\alpha/(1-\alpha)$、猜想     |

## 本页原文

- S. Wu, H. Zhang, and T. Zhou, [_Solving time-periodic fractional diffusion equations via diagonalization technique and multigrid_](https://doi.org/10.1002/nla.2178), Numer. Linear Algebra Appl. 25(5) (2018), e2178。
- S. Wu and T. Zhou, [_Acceleration of the two-level MGRIT algorithm via the diagonalization technique_](https://doi.org/10.1137/18M1207697), SIAM J. Sci. Comput. 41(5) (2019), pp. A3421-A3448。
- S. Wu and T. Zhou, [_Diagonalization-based parallel-in-time algorithms for parabolic PDE-constrained optimization problems_](https://doi.org/10.1051/cocv/2020012), ESAIM Control Optim. Calc. Var. 26 (2020), 88。
- S. Wu and T. Zhou, [_Parallel implementation for the two-stage SDIRK methods via diagonalization_](https://doi.org/10.1016/j.jcp.2020.110076), J. Comput. Phys. 428 (2021), 110076。
