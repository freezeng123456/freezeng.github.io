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

本页六篇论文可以按「用什么当行、用什么当列、行怎么归一化」来读：

- 编号 **10** 把行换成一组**确定性**的点，代价是相干性路线带来的 $m\gtrsim s^2$；
- 编号 **21** 把行取成张量 Gauss 网格的随机子集，发现 Gauss 权本身就是 Christoffel 权；
- 编号 **29**、**32** 在每个样本处**多取 $d$ 行**（梯度），代价是行范数不再齐性，必须补预条件；
- 编号 **36**、**44** 换的是**列**：正交基不再由已知密度给出，而要从经验测度估计出来。

> [!note] 本页的符号约定
> **Christoffel 约定。** 全页统一采用站内约定：再生核对角线记为 $K(z)=\sum_{\alpha\in\Lambda}\varphi_\alpha^2(z)$，Christoffel 函数是它的倒数 $\lambda_\Lambda(z)=1/K(z)$，归一化 Christoffel 函数是 $N/K(z)$，$N=|\Lambda|$。因此本页所说的「Christoffel 权」一律指 $1/K$；最小二乘目标里出现的 $N/K$ 与 $1/K$ 只差一个与 $z$ 无关的常数，不改变极小元。原文各家的符号并不一致：编号 36 用符号 $K$ 表示 $N/\sum_j\Phi_j^2$（即归一化 Christoffel 函数本身），编号 44 用 $\kappa$ 表示 $\frac1N\sum_j\Phi_j^2$（即倒数）。凡引用原式处本页都换算到上面的约定并注明。
>
> **样本数与基函数个数。** 编号 10、21、36、44 采用「样本数在前」：$m$ 或 $M$ 个样本、$N$ 个基函数；编号 29、32 恰好相反：$N$ 个样本、$M$ 个基函数。本页在每一节内沿用该论文自己的记号并在节首标明，跨节比较时请留意。另外编号 29 的 $\Gamma$ 指频率指标集 $\Gamma\subset\mathbb Z^d$，而编号 21、32 的 $\Gamma$ 指参数域。

### 两条恢复路线

六篇论文的理论保证只走两条路，读的时候先分清是哪一条。

**相干性路线。** 设计矩阵 $A$ 的互相干性为

$$
\mu(A):=\max_{k\ne j}\frac{|\langle A_k,A_j\rangle|}{\|A_k\|_2\|A_j\|_2},
$$

$A_k$ 为第 $k$ 列。若 $c_0$ 是 $s$-稀疏的且 $\mu<\frac{1}{2s-1}$，则 $\ell_1$ 极小化精确恢复 $c_0$（Donoho–Huo 关于两组正交基之并的结果，经 Fuchs 与 Gribonval–Nielsen 推广，含噪时给出稳定恢复）。等价地可以经限制等距常数走：$\delta_s\le(s-1)\mu(A)$，再配合 $\delta_s<\frac13$ 的恢复定理。这条路线的算术后果是硬性的：$\mu$ 至少是 $O(1/\sqrt m)$ 量级，所以 $\mu<\frac{1}{2s-1}$ 迫使 $m\gtrsim s^2$，**样本数在稀疏度上是二次的**。编号 10 与 32 走的是这条路，编号 29 的定理 2.1 也是。

**有界正交系的限制等距路线。** 限制等距常数 $\delta_s<1$ 是使

$$
(1-\delta_s)\|c\|_2^2\ \le\ \|Dc\|_2^2\ \le\ (1+\delta_s)\|c\|_2^2
$$

对一切 $\|c\|_0\le s$ 成立的最小数。若 $\delta_s\le0.307$，则 $c^\sharp=\arg\min\|c\|_1$（约束 $Dc=D\tilde c$）满足 $\|c^\sharp-\tilde c\|_2\le C\sigma_{s,1}(\tilde c)/\sqrt s$，其中 $\sigma_{s,p}(c)=\inf_{\|y\|_0\le s}\|y-c\|_p$ 是最佳 $s$ 项逼近误差；$\tilde c$ 恰为 $s$-稀疏时精确恢复。一族关于密度 $\nu$ 正交归一的函数 $\{\psi_k\}$ 若满足

$$
\max_{1\le k\le N}\|\psi_k\|_\infty^2=\max_{1\le k\le N}\ \sup_{x\in\operatorname{supp}\nu}|\psi_k(x)|^2\ \le\ L(N)<\infty,
$$

就称为**有界正交系**；若该界与 $N$ 无关则称一致有界。Rauhut 与 Rauhut–Ward 的定理给出：$x_1,\dots,x_M$ 独立同分布取自 $\nu$，$d_{ij}=\psi_j(x_i)$，则

$$
M\ \ge\ C\,\delta^{-2}\,L\,s\,\log^3(s)\,\log(N)
$$

蕴含 $\frac{1}{\sqrt M}D$ 的限制等距常数 $\delta_s\le\delta$，概率至少 $1-N^{-\gamma\log^3(s)}$，$C,\gamma$ 为普适常数。**这条路线把样本数从 $s^2$ 降到 $s$ 的对数倍数，代价是概率性结论，以及必须把 $L$ 控制住——而 $L$ 正是「行范数的均匀性」。** 编号 21 与编号 29 的定理 3.1 走的是这条路。

Legendre 多项式**不是**一致有界的，但乘上 $(1-x^2)^{1/4}$ 之后是；保持加权系正交性的偏置测度恰是 Chebyshev 测度。这就是 Rauhut–Ward 策略，也是本页反复出现的**倒数配对**原则：按 $\propto w/(N\lambda_\Lambda)$ 的密度采样，再用 $N\lambda_\Lambda$ 做预条件，使复合系一致有界。编号 21 用 Gauss 权实现它，编号 32 用 $\sqrt{\rho^{(\alpha,\beta)}/\rho_c}$ 实现它，编号 44 用 $\mu\propto\kappa$ 采样、$1/\kappa$ 加权实现它。

## 10：确定性插值点与稀疏插值

（本节记号：$m$ 个插值点，$N=\#\Lambda$ 个基函数，$s$ 为稀疏度。）

### 直觉

论文同时回答两个互不相干的问题。

第一个是纯存在性的：**$s$-稀疏插值最少需要多少个点？** 直觉很干净。若 $f,g\in U_s$ 在这些点上取值相同，则 $f-g$ 至多 $2s$-稀疏且在这些点上为零；于是「唯一可解」等价于「没有非零的 $2s$-稀疏函数在这些点上全零」。而任取 $2s$ 个基函数与 $2s$ 个点组成的方阵，其行列式作为点的函数不恒为零（这正是「强线性无关」这一假设的内容），因此它的零点集是一个零测集。$2s$ 元子集只有有限多个（$\binom{N}{2s}$ 个），有限个零测集之并仍是零测集，取补集中的点即可。**结论是 $2s$ 个点总够用，且与 $N$ 无关。**

第二个问题是构造性的：存在性证明是「避开一个零测集」，它不给出任何具体点集；而支撑高维 Chebyshev 基下 $\ell_1$ 恢复的压缩感知理论用的是**随机**点，保证是概率性的。论文因此给出一个显式的**确定性**点集，并用 Weil 指数和估计它的互相干性。**代价是走相干性路线必然带来的 $m\gtrsim s^2$；换来的是没有「以高概率」这一前提。**

### 问题设定

设 $\Omega\subset\mathbb R^d$，$\{B_j\}_{j\in\Lambda}$ 为 $N:=\#\Lambda$ 个复值基函数，$s$-稀疏类为

$$
U_s:=\Bigl\{f=\sum_{j\in T}c_jB_j\ :\ T\subset\Lambda,\ \#T\le s\Bigr\}.
$$

称 $\{x_1,\dots,x_m\}\subset\Omega$ 对 $U_s$ **唯一可解**（unisolvent），若由 $f(x_j)=g(x_j)$（$j=1,\dots,m$）与 $f,g\in U_s$ 可推出 $f\equiv g$。插值矩阵为 $A:=[B_j(x_t)]_{t=1,\dots,m;\,j\in\Lambda}$，数据为 $b:=[f(x_1),\dots,f(x_m)]^\top$。

### 推导

**第一步：把唯一可解性换成一个测度论条件。** 称 $f_1,\dots,f_k$ 在 $\Omega_0$ 上**强线性无关**，若 $\sum_tc_tf_t$ 在 $\Omega_0$ 上恒为零蕴含 $c_1=\cdots=c_k=0$。引理 3.1 证明这等价于行列式退化集

$$
S:=\{x=(x_1,\dots,x_k)\in\Omega_0^{\,k}:\det(A_x)=0\},
\qquad A_x:=[f_t(x_j)]_{j,t=1,\dots,k},
$$

的外 Lebesgue 测度为零，$\lambda^*_{d\cdot k}(S)=0$。**这一步是整套存在性论证的机关**：它把一个代数条件变成「$S$ 是零测集」，于是可以做并集与取补。

**第二步：并集论证。** 对每个大小为 $2s$ 的指标集 $T$，由假设 $\{B_j\}_{j\in T}$ 强线性无关，故对应的 $S_T$ 零测。$T$ 只有 $\binom{N}{2s}$ 个，有限并仍零测，故 $\Omega_0^{2s}\setminus\bigcup_TS_T$ 非空，其中任一点给出 $2s$ 个唯一可解的插值点。

**第三步：显式点集。** 论文取与[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 9]] 同源的 Weil 和网格

$$
\Theta_M=\Bigl\{x_j=\cos(p_j)\ :\ p_j=\tfrac{2\pi}{M}\bigl(j,j^2,\dots,j^d\bigr),\ j=0,\dots,\lfloor M/2\rfloor\Bigr\},
\qquad m:=\#\Theta_M=\lfloor M/2\rfloor+1 .
$$

指标只跑到 $\lfloor M/2\rfloor$ 是有理由的：引理 4.1 指出，对任意整数 $M$ 与 $m=\lfloor M/2\rfloor+1$，

$$
\cos\bigl(2\pi j^k/M\bigr)=\cos\bigl(2\pi(M-j)^k/M\bigr)
\qquad\text{对一切 }k\in\mathbb N,\ 0\le j\le m-1,
$$

即后一半指标完全重复前一半，取了也没有新信息。

**第四步：插值矩阵与列归一化。** 记 $\Phi_n(x_j)=C_n(p_j)$，$C_n(p_j):=\prod_{t=1}^{d}\cos(2\pi n_tj^t/M)$。张量积指标集 $\Lambda^{q,d}_P$ 上的矩阵是 $A_P:=[C_n(p_j)]\in\mathbb R^{m\times(q+1)^d}$，全次数指标集 $\Lambda^{q,d}_D$ 上的是 $A_D$。恢复问题解在 $A:=A_P\cdot C$（或 $A_D\cdot C$）上，其中 $C$ 是把各列 $\ell_2$ 范数标准化为 1 的对角矩阵。**这个不起眼的列归一化，就是[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 24]] 里系统化设计的预条件矩阵的祖先。**

**第五步：相干性链。** 对素数 $p$ 的 **Weil 指数和定理**给出三角和的 $O(\sqrt p)$ 上界，代入互相干性定义得到引理 4.2 的显式界；再用 $\delta_s\le(s-1)\mu(A)$ 与 $\delta_s<\frac13$ 的恢复条件，反解出对 $M$ 的要求。整条链是显式的，没有任何概率步骤。

### 定理

**定理 3.2（$2s$ 个唯一可解点的存在性）。** 若 $s\le N/2$，且存在 $\Omega_0\subset\Omega\subset\mathbb R^d$ 使 $\{B_j\}_{j\in\Lambda}$ 中任意 $2s$ 个函数在 $\Omega_0$ 上强线性无关，则存在 $2s$ 个点 $\{x_1,\dots,x_{2s}\}\subset\Omega_0$ 对基 $\{B_j\}$ 下的 $s$-稀疏插值唯一可解。**$2s$ 是理论下限，与 $N$ 无关。**

**定理 3.5（Chebyshev 系下任意 $2s$ 个相异点都行）。** 以下两条等价：(i) 对 $f,g\in U_s$ 与**任意** $2s$ 个相异点 $x_1,\dots,x_{2s}\in\Omega$，由 $f(x_j)=g(x_j)$ 可推出 $f\equiv g$；(ii) 对每个 $\#T=2s$ 的指标集 $T$，$\{B_j\}_{j\in T}$ 构成 Chebyshev 系。论文同时指出：$d\ge2$ 时 $[-1,1]^d$ 上本质上不存在 Chebyshev 系，**因此这幅干净的一维图景不能移植到多维**——这正是论文另一半（确定性设计）存在的理由。

**一条负面结果。** 论文举出的一个设定中，唯一可解性反过来迫使 $m\ge N$。也就是说，**仅有稀疏性并不总能换来点数的减少**。

**引理 4.2（张量积情形的相干性界）。** 设 $M\ge\max\{2q+1,\ (2d(d-1))^2\}$ 为素数，则

$$
\mu(A_P)\ \le\ \frac{1}{\sqrt{M}}\cdot\frac{2^d\,d}{1-\dfrac{2d(d-1)}{\sqrt{M}}}.
$$

**定理 4.3（张量积 Chebyshev 的 $\ell_1$ 恢复）。** 设 $M\ge\max\{2q+1,\ 9\cdot4^d\cdot d^2\cdot s^2\}$ 为素数，$f=\sum_{n\in\Lambda^{q,d}_P}c_n\Phi_n$，$c^{\#}$ 是以 $A=A_P\cdot C$ 与 $b=(f(x_1),\dots,f(x_m))^\top$（$x_j\in\Theta_M$，$m=\lfloor M/2\rfloor+1$）为数据的 $\ell_1$ 问题的解，则

$$
\|c^{\#}-c\|_2\ \lesssim\ \frac{\sigma_{s,1}(c)}{\sqrt{s}} .
$$

**注 4.4** 给出精确恢复的形式：$c$ 恰为 $s$-稀疏时，只要 $m\ge\max\{q,\ \tfrac92\cdot4^d\cdot d^2\cdot s^2\}+1$ 且 $M$ 为素数，恢复精确。

**引理 4.5 与定理 4.6（全次数情形）。** 引理 4.5 假设 $d\ge q$ 且 $M\ge\max\{2q+1,\ (2q(d-1))^2\}$ 为素数。定理 4.6：对素数 $M\ge9\cdot4^{\,q}\cdot d^2\cdot s^2$ 与 $f=\sum_{n\in\Lambda^{q,d}_D}c_n\Phi_n$，以 $A=A_D\cdot C$ 求解的 $\ell_1$ 解满足同样的 $\|c^{\#}-c\|_2\lesssim\sigma_{s,1}(c)/\sqrt s$；注 4.7 给出精确恢复条件 $m\ge\tfrac92\cdot4^{\,q}\cdot d^2s^2+1$。

| 指标集                   | 相干性引理的前提                                     | 精确恢复所需点数                             |
| ------------------------ | ---------------------------------------------------- | -------------------------------------------- |
| 张量积 $\Lambda^{q,d}_P$ | $M$ 为素数，$M\ge\max\{2q+1,(2d(d-1))^2\}$           | $m\ge\max\{q,\ \tfrac92\cdot4^{d}d^2s^2\}+1$ |
| 全次数 $\Lambda^{q,d}_D$ | $d\ge q$，$M$ 为素数，$M\ge\max\{2q+1,(2q(d-1))^2\}$ | $m\ge\tfrac92\cdot4^{\,q}d^2s^2+1$           |

**这张表里唯一重要的差别是指数从 $4^d$ 变成 $4^{\,q}$**：从维数上的指数变成多项式次数上的指数。不确定性量化里 $d$ 动辄上百而 $q$ 通常很小，因此全次数版本才是可用的那个。

**注 4.8** 补了一个几何解释：当 $4\mid M$ 时，$p_j$ 的第一个坐标在奇数 $j\in[1,M/2]$ 上构成一组 Chebyshev 节点，所以 $\Theta_M$ 可以读作 Chebyshev 节点的一种高维推广。

### 数值实验

第 5 节把 $\Theta_M$ 与 $[-1,1]^d$ 上独立同分布均匀随机点作对照。

| 项目     | 设置                                             |
| -------- | ------------------------------------------------ |
| 对照     | 确定性点集 $\Theta_M$ 对独立同分布均匀随机点     |
| 空间     | 张量积 Chebyshev 空间与全次数 Chebyshev 空间两组 |
| 测试函数 | 支撑在 $\Lambda$ 的全部 $s$ 元子集上均匀抽取     |
| 系数     | 非零系数独立同分布标准正态（均值 0，标准差 1）   |
| 求解器   | MATLAB 中的 SPGL1                                |
| 统计     | 每个固定稀疏度 $s$ 重复 100 次，统计经验成功率   |

作者报告的结论是：**确定性点集与随机点表现相当**。也就是说，实验建立的不是精度上的优势，而是「确定性并不需要付出精度代价」这一点——论文自己在摘要里也是这么说的（"a similar performance"）。

> [!warning] 这里没有可转录的数字
> 论文以经验成功率曲线的形式发表结果，本站掌握的材料只记录了实验配置与作者陈述的定性结论，**没有逐点的成功率数值**。不要从上表推断任何具体的成功概率。与理论的差距也应当说明：定理 4.3 要求 $m\gtrsim4^dd^2s^2$，实验中实际使用的点数远小于这个量级，因此实验检验的是方法在理论保证之外的实际表现，而不是定理本身的紧性。

### 与其他论文的关系

编号 10 是[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 9]] 的稀疏恢复孪生篇：同一批作者、同一个 Weil 和点集 $\Theta_M$、同一种「用确定性替换概率性」的取向，只是把最小二乘与 Gershgorin 换成了 $\ell_1$ 与互相干性。它的 $m\gtrsim s^2$ 正是相干性路线的二次壁垒，而编号 21、24、29、32 用限制等距路线绕开的就是这道壁垒。Xu 与 Zhou 在编号 **29** 里回到同一套 Weil 机制处理三角多项式，并在那里加上梯度测量；编号 29 的命题 2.1 用的 Gauss 和计算 $\mu(\Psi)=1/\sqrt p$，与本文的相干性估计出自同一处。列归一化矩阵 $C$ 则成为[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 24]] 中系统设计的预条件。

## 21：对 Gauss 网格做随机子采样后做 ℓ1 极小化

（本节记号：$M$ 个样本，$N=|\Lambda|$ 个基函数。）

### 直觉

Tang 与 Iaccarino 已经证明：对张量积 Gauss 求积网格做子采样，可以得到一个好的压缩感知设计——但只对超立方体上的**均匀**随机变量成立，充分样本数为 $M\gtrsim3^ds$。悬而未决的问题是：不确定性量化里真正要用的那些分布（尤其是无界的正态与指数）是否同样可行？

论文的关键观察是一个恒等式而不是一个新构造：**Gauss 求积权恰好就是 Christoffel 函数值**。于是「按 $\sqrt{w}$ 给设计矩阵加权」这件事，同时是两件事——它是有界正交系理论所要的行归一化，也是从 Gauss 求积继承下来的自然权。加权之后的函数族在网格上的均匀经验测度下正交归一，而从这个测度独立采样就是「从张量 Gauss 网格里等概率挑一个点」，实现起来毫不费力。**整个问题于是被压缩成一个纯粹的分析问题：Christoffel 加权后的族的一致上界 $L$ 有多大？** 有界（Beta）参数下 $L$ 与次数无关，无界参数下 $L$ 以 $n^{2/3}$ 增长——这就是论文全部的答案。

### 问题设定

$X=(X_1,\dots,X_d)^\top$ 分量相互独立，边缘密度 $\rho_i$ 定义在 $\Gamma_i$ 上，联合密度 $\rho(x)=\prod_i\rho_i(x_i)$ 定义在 $\Gamma=\otimes_i\Gamma_i$ 上。一维广义多项式混沌基由

$$
\mathbb E\bigl[\phi^i_n(X_i)\phi^i_\ell(X_i)\bigr]=\int_{\Gamma_i}\phi^i_n\phi^i_\ell\,\rho_i\,ds=\delta_{n,\ell}
$$

定义；边缘分布可以是 Beta（Legendre / Chebyshev / Jacobi）、正态（Hermite）或指数（Laguerre）。指标集 $\Lambda$ 一般，张量指标集记 $\Lambda^P_{\mathbf n}=\{k\in\mathbb N_0^d:k\le\mathbf n\}$。

### 推导

**Gauss 网格与它的 Christoffel 权。** 记 $\Theta_{\mathbf n}=\Theta^1_{n_1}\otimes\cdots\otimes\Theta^d_{n_d}$，$|\Theta_{\mathbf n}|=\prod_in_i$，它精确积分张量空间 $\mathbb P_{2\mathbf n-1}$ 中的一切多项式。它的权为

$$
w_{\mathbf k}=\lambda_{\mathbf n}(z_{\mathbf k})
\triangleq\prod_{i=1}^{d}\lambda^i_{n_i}(z^i_{k_i})
=\prod_{i=1}^{d}\frac{1}{\sum_{k=0}^{n_i-1}\bigl[\phi^i_k(z^i_{k_i})\bigr]^2},
$$

**即 Gauss 权恰好是 Christoffel 函数值**（在本页约定下就是 $1/K$ 的张量形式）。网格上的均匀经验概率测度是

$$
\nu_{\mathbf n}=\bigotimes_{i=1}^{d}\nu^i_{n_i}
=\frac{1}{\prod_i n_i}\sum_{\mathbf k\le\mathbf n}\delta_{z_{\mathbf k}},
$$

于是「从 $\nu_{\mathbf n}$ 独立同分布采样」等价于「从张量 Gauss 网格均匀采样」。

**加权测量矩阵。** 令 $(\Psi)_{m,n}=\varphi_n(x_m)$，$W$ 为对角阵、$(W)_{m,m}=w_m>0$，置 $D=\sqrt W\Psi$。未加权问题是

$$
\arg\min\|c\|_1\quad\text{s.t.}\quad \Psi c=f,
$$

加权问题是

$$
\arg\min\|c\|_1\quad\text{s.t.}\quad Dc=\sqrt W f .
$$

**离散正交归一性。** $D$ 的行是被 Christoffel 函数加权过的多项式。对任意 $\Lambda\subseteq\Lambda^P_{\mathbf n-1}$、$N=|\Lambda|$，引理 4.1.D 给出：这 $N$ 个函数 $\{\psi_{k,\mathbf n}(z)\}_{k\in\Lambda}$ 在离散测度 $\nu_{\mathbf n}$ 下**正交归一**。引理 3.1.A/B 进一步给出矩阵层面的版本：$D^i=(\Sigma^i)^{1/2}\Psi^i$ 是正交矩阵，张量积 $D=\bigotimes_iD^i$ 也是——Tang 与 Iaccarino 称之为「离散正交矩阵」；引理 3.1.C 把结论推广到 $\mathbf m\ge\mathbf n$（求积点多于基函数）的情形。**这正是「行归一化」在张量结构下的具体实现。**

**算法。** (1) 给定 $\Lambda$，找出使 $\Lambda\subseteq\Lambda^P_{\mathbf n-1}$ 的 $\mathbf n$；(2) 生成各维的 $n_i$ 点 Gauss 规则——**完整的张量规则 $\Theta_{\mathbf n}$ 无需显式构造**；(3) 从 $\Theta_{\mathbf n}$ 中均匀随机取 $M$ 个点 $\{(y_m,v_m)\}_{m=1}^M$；(4) 组装 $(D)_{m,n}=\sqrt{v_m}\,\varphi_{k(n)}(y_m)$ 与 $(W)_{m,m}=v_m$；(5) 解加权 $\ell_1$ 问题。第 (2) 步是这个设计在高维下可行的原因：$\prod_in_i$ 个点里只有 $M$ 个会被真正生成。

### 定理

**一致有界性的三个引理。** 这是全文的技术核心，它决定了 $L(\mathbf n)$ 的大小。

| 边缘分布                                                                 | 条件                                  | $L_i(n)$                                |
| ------------------------------------------------------------------------ | ------------------------------------- | --------------------------------------- |
| Beta / Jacobi，$B(\gamma+1,\delta+1)$，$\gamma,\delta\ge-\tfrac12$       | 有界区间 $[-1,1]$                     | $\le C(\gamma,\delta)$，**与 $n$ 无关** |
| 双边指数型 $\rho_i\propto e^{-\lvert x\rvert^\alpha}$，$\alpha>\tfrac32$ | $\mathbb R$（含正态 $\alpha=2$）      | $\le C(\alpha)\,n^{2/3}$                |
| 单边指数型 $\rho_i\propto e^{-\lvert x\rvert^\alpha}$，$\alpha>\tfrac34$ | $[0,\infty)$（含指数分布 $\alpha=1$） | $\le C(\alpha)\,n^{2/3}$                |

注 4.1 说明：更一般的权 $\rho_i\propto x^\mu e^{-|x|^\alpha}$（$\mu\ge-\frac12$）预期也成立，**但本文没有证明**。

**定理 4.1（主定理）。** 设 $\Lambda$ 有限、$|\Lambda|=N$，$\mathbf n$ 是使 $\Lambda\subseteq\Lambda^P_{\mathbf n-1}$ 的最小多重指标。从 $\nu_{\mathbf n}$ 中**不放回**随机抽取 $M$ 个样本。若

$$
M\ \ge\ L(\mathbf n)\,C_1\,s\log^3(s)\log(N),
\qquad
L(\mathbf n)=\prod_{i=1}^{d}L_i(n_i),
$$

$C_1$ 普适，则对任意 $c\in\mathbb R^N$，加权 $\ell_1$ 问题的解 $c^\sharp$ 以至少 $1-N^{-\gamma\log^3(s)}$ 的概率满足

$$
\|c-c^\sharp\|_2\le C_2\,\frac{\sigma_{s,1}(c)}{\sqrt s},
$$

$C_2,\gamma$ 普适。**一处容易被忽略的设定：这 $M$ 个样本是从 $\nu_{\mathbf n}$ 中不放回抽取的**，而不是独立同分布——上面那条等价说法讲的是设计的构造，抽样环节本身用的是不放回。

> [!warning] 原文印出的概率不等式方向是反的
> 论文正文把定理 4.1 的结论印成
>
> $$
> \Pr\Bigl[\|c-c^\sharp\|_2\le C_2\tfrac{\sigma_{s,1}(c)}{\sqrt s}\Bigr]\ \le\ 1-N^{-\gamma\log^3(s)},
> $$
>
> 即「概率**不超过** $1-N^{-\gamma\log^3 s}$」。这与「以高概率成立」的本意相反，而且与同一篇论文的定理 2.2（"with probability at least $1-N^{-\gamma\log^3(s)}$"）自相矛盾。**本页采用修正后的方向 $\Pr[\cdot]\ \ge\ 1-N^{-\gamma\log^3(s)}$**，并在此标明原文是这样印的。这几乎可以肯定是排印错误，但引用时请以修正后的方向为准。

**三类参数的样本复杂度。** 论文自己把结论总结成 $M\gtrsim L(\mathbf n)s$：

| 参数类型                                                          | $L(\mathbf n)$      | 样本要求                               |
| ----------------------------------------------------------------- | ------------------- | -------------------------------------- |
| Beta（有界）                                                      | $\le C^d$           | $M\gtrsim C^d\,s\log^3 s\log N$        |
| 正态（各维最高次数同为 $n-1$）                                    | $\le(Cn)^{2d/3}$    | $M\gtrsim n^{2d/3}s\log^3 s\log N$     |
| 单边指数（$\rho\propto e^{-\lVert z\rVert_1}$ 于 $[0,\infty)^d$） | 同上                | 同上                                   |
| 正态，次数 $n-1=9$ 的经验观察                                     | $Cn^{2/3}\lesssim4$ | $M\ge4^d\,s\log^3 s\log N$（$n\le10$） |

有界情形下 Tang 与 Iaccarino 给出的常数本质上是 $C=3$；论文的图 1 **提示** $C=2$ 可能更紧，但这一点**没有被证明**。

论文对无界情形非常坦率：$n^{2d/3}$ 这个依赖「看上去不太愉快」（"seems unpleasant"），但在这一套分析策略、并且坚持子采样张量 Gauss 网格的前提下，它**本质上是紧的**。随后论文用经验数据把它软化：次数 $n-1=9$ 时观察到 $Cn^{2/3}\lesssim4$，于是 $n\le10$ 的九次逼近所需样本数读作 $M\ge4^ds\log^3(s)\log(N)$，与有界情形可比；而高维问题本来也只用低次逼近。

最后一条注记值得单独记住：可以从 $\nu_{\mathbf m}$（$\mathbf m\ge\mathbf n$）中子采样，**有界情形下界不变**，但指数型密度下界退化成 $L\sim m^{2/3}>n^{2/3}$。原因很具体：更大的 Gauss 规则会把点放到那些「$n$ 次多项式被 $\lambda_{\mathbf m}$ 加权后迅速衰减到零」的区域，从那里恢复系数是困难的。**换句话说，网格取得越大越不划算——这与直觉相反，值得注意。**

论文还转述了 Yan–Guo–Xiu 的定理 2.3：全次数空间 $T^d_n$ 上的多维 Legendre 预条件权为 $w_m=(\tfrac2\pi)^d\prod_{n=1}^{d}\bigl(1-(x^n_i)^2\bigr)^{1/2}$。

### 数值实验

$\ell_1$ 问题用 MATLAB 的谱投影梯度算法 SPGL1 求解。测试系数向量的构造是：固定稀疏度 $s$，取 $s$ 个独立同分布标准正态的非零分量，其余置零。论文比较五种恢复流程的经验成功概率：

| 标签         | 采样                                                                     | 求解的问题           |
| ------------ | ------------------------------------------------------------------------ | -------------------- |
| Random       | 从 $\rho$ 独立同分布                                                     | 未加权               |
| PreChebyshev | 从 Chebyshev 密度 $v(x)=\frac{1}{\pi^d\prod_i\sqrt{1-x_i^2}}$ 独立同分布 | 加权（预条件）       |
| **Gaussian** | 子采样 Gauss 网格，用 Gauss 权                                           | 加权（**本文方法**） |
| Chebyshev    | 从 Chebyshev 测度独立同分布                                              | 未加权               |
| Uniform      | 从均匀测度独立同分布                                                     | 未加权               |

图 2 在固定 $M=85$ 下画 Legendre 基的恢复概率对稀疏度 $s$ 的曲线，两种配置：

| 配置     | $d$ | 各维次数上限 $n$ | $N$ | $M$ |
| -------- | --- | ---------------- | --- | --- |
| 低维高次 | 2   | 21               | 231 | 85  |
| 高维低次 | 10  | 4                | 286 | 85  |

图 1 则画理论界 $L$ 本身：对称参数 $\gamma=\delta$ 的 Jacobi 族，以及两个指数型密度（$\mathbb R$ 上的 $\rho=\exp(-x^2)$ 对应 Hermite、$[0,\infty)$ 上的 $\rho=\exp(-x)$ 对应 Laguerre）。

> [!warning] 实验的可转录范围
> 上表的配置数字（$M=85$，$d=2,n=21,N=231$；$d=10,n=4,N=286$）是可核实的，**但成功概率曲线的具体数值本站没有转录**。此外要注意实验与理论之间的落差：定理 4.1 在 $d=10$ 时要求的 $M$ 含有 $C^{10}$ 这样的因子，而实验用的是 $M=85$。**实验说明的是方法在远小于理论样本数时仍然可用，而不是理论界是紧的。** 论文对 $C=2$ 是否成立也只说图 1「提示」，没有下结论。

### 与其他论文的关系

编号 21 是[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 13]] 的 $\ell_1$ 兄弟篇：同样的随机子采样 Gauss 求积设计，那里是最小二乘，这里是稀疏恢复。它把 Tang 与 Iaccarino 只覆盖均匀分布的结论推广到 Beta、正态与指数参数。它的中心观察——Gauss 权**就是** Christoffel 函数值，因此子采样 Gauss 网格隐含地就是 Christoffel 加权采样——在[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 22 与 45]] 中被系统展开：那里采样密度本身取成 Christoffel 函数的倒数，而不是从求积规则里继承。

维数依赖全部集中在 $L(\mathbf n)$ 里，这与[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 24]] 的 $L(n)$ 是同一种结构；限制等距 / 有界正交系机制则与编号 10、24、29、32 共享。论文量化的无界域困难（$n^{2d/3}$），是[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 11]] 在最小二乘情形遇到的同一困难在稀疏恢复框架下的表现。

## 梯度增强的公共直觉：多的是行，不是样本

编号 29 与 32 共享同一个出发点，值得在进入各自的技术细节之前先说清楚。

**一次昂贵的正问题求解通常同时给出函数值与全部偏导数**——例如通过伴随方法求解一次伴随方程，或者对整个求解器做自动微分。若只用函数值，每个样本贡献设计矩阵的 1 行；若同时用梯度，同一个样本贡献 $1+d$ 行，而未知量个数不变。因此测量矩阵由

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

其中 $W$ 是把各行范数拉平的对角权矩阵。**说清楚收益在哪里很重要：收益是「每次求解得到的行数」，不是「样本数」。** 样本点还是那 $M$ 个，求解器还是被调用 $M$ 次；变的是从每次调用里榨出的线性方程个数。**这就把「每个昂贵样本的信息量」变成一个可以设计的对象。**

有两个问题必须处理。

第一，**行范数不再一致**。值行与梯度行的量级不同——梯度行带一个导数的尺度，而正交多项式的导数随次数增长——因此预条件不再是可选项：若不加权，限制等距常数或互相干性会被最坏的一类行支配。编号 32 的数值实验直接证实了这一点：**朴素堆叠梯度行反而破坏原矩阵的稳定性**。这正是编号 32 与[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 24]] 共享「采样与预条件成对设计」这一取向的原因。

第二，**梯度并非免费**。若梯度需要额外的求解（而不是伴随方法的副产品），则 $(1+d)$ 倍的行数要与相应的代价对比。编号 32 在实验里用 "standard-double" 这一对照直面这个问题：给标准方法同样多的标量数据（$d=2$ 时是 $3N$ 个函数值），看梯度方法是否仍占优。

两篇的差别在**基**，而这个差别是本质性的：

- 编号 29 用 Fourier 基。由于 $D_{v}e^{ik\cdot z}=\langle v,k\rangle e^{ik\cdot z}$，**梯度行是值行的标量倍**，所以梯度块不带来任何新的随机性，只带来一个逐列的已知缩放。相干性分析因此是纯代数的（一个精细化的 Cauchy–Schwarz），结论**确定性**成立，对任意采样点都成立。
- 编号 32 用正交多项式基。$\frac{d}{dx}p^{(\alpha,\beta)}_n=c(n,\alpha,\beta)\,p^{(\alpha+1,\beta+1)}_{n-1}$ 是**另一族**多项式、带**另一个**正交权，梯度行不是值行的倍数。因此必须重新设计预条件才能恢复各向同性。

## 29：稀疏三角多项式的梯度增强 ℓ1 恢复

（本节记号沿用论文：$N$ 个采样点，$M=\#\Gamma$ 个基函数，$\Gamma\subset\mathbb Z^d$ 是频率指标集，$s$ 为稀疏度。论文用 $k$ 同时表示方向个数与频率多重指标，本页把方向个数改记为 $r$ 以免混淆。）

### 直觉

对稀疏 Fourier 恢复而言，$\ell_1$ 极小化的成败取决于互相干性 $\mu(\Phi)$ 或限制等距常数。论文问的是两个具体问题：在同样的采样点上追加**方向导数**测量，能否严格改善恢复条件？以及——这是更有意思的一问——导数测量能否**替代**函数值？

第一问的机制可以一眼看穿。频率 $k$ 对应的列在堆叠后变成

$$
\bigl(a_k,\ \langle v_1,k\rangle a_k,\ \dots,\ \langle v_r,k\rangle a_k\bigr),
\qquad a_k=\bigl[e^{ik\cdot z_j}:z_j\in\Xi\bigr],
$$

于是两列的内积多出因子 $1+\sum_t\langle v_t,k\rangle\langle v_t,k'\rangle$，两列的范数多出因子 $1+\sum_t\langle v_t,k\rangle^2$。归一化之后，相干性等于两个**增广向量** $(1,V^\top k)$ 与 $(1,V^\top k')$ 之间夹角余弦的最大值。只要这两个向量不成比例，夹角就严格大于原来的、由 $a_k$ 与 $a_{k'}$ 决定的夹角——**而「不成比例」正好就是论文所定义的可容许性**。所以改进是纯几何的：加梯度等于把每个频率提升到一个更高维的空间，在那里不同频率被推得更开。

**先前的工作要么只有数值结果，要么（Hermite 多项式加随机样本）只给出几乎必然的严格不等式 $\mu(\tilde\Phi)\le\mu(\Phi)$。本文的改进是给出确定性的 $\mu(\tilde\Phi)\le\lambda\mu(\Phi)$，$\lambda<1$ 显式，且不带任何概率前提。**

### 问题设定

$f(x)=\sum_{k\in\Gamma}c_ke^{ik\cdot x}$，$x\in[-\pi,\pi)^d$，$\Gamma\subset\mathbb Z^d$ 有限、$M=\#\Gamma$；重新编号后 $f=\sum_{j=1}^Mc_j\varphi_j(x)$。采样点 $\Xi=\{z_j\}_{j=1}^N\subset[-\pi,\pi)^d$，插值矩阵 $\Phi=(\varphi_t(z_j))\in\mathbb C^{N\times M}$。支撑 $T=\{k:c_k\ne0\}$，稀疏度 $s=\#T\ll M$，$\Pi_s(\Gamma)=\bigcup_{\#T\le s}\Pi(T)$。

论文分两个问题：

- **问题 1：** 从 $f(z_j)=f_j$ **以及** $D_{v_t}f(z_j)=f'_{j,t}$（$t=1,\dots,r$，$r\le d$，$D_{v}f(x)=\langle\nabla f(x),v\rangle$）恢复 $f\in\Pi_s(\Gamma)$。
- **问题 2：** 从 $D^{\tau_j}_{v_j}f(z_j)=y_j$（$j=1,\dots,N$，$v_j\in\mathbb R^d$，$\tau_j\in\mathbb Z_{\ge0}$）恢复 $f$——即每个点上**要么**知道函数值（$\tau_j=0$），**要么**知道某阶方向导数，但不同时知道。

### 推导

**梯度增强矩阵。** 记

$$
\Phi=\bigl(e^{ik\cdot z_j}\bigr)_{z_j\in\Xi,\,k\in\Gamma},
\qquad
\Phi_t=\bigl(D_{v_t}e^{ik\cdot z_j}\bigr)_{z_j\in\Xi,\,k\in\Gamma},
\qquad
\tilde\Phi:=\begin{pmatrix}\Phi\\ \Phi_1\\ \vdots\\ \Phi_r\end{pmatrix}.
$$

关键的初等恒等式是 $D_{v_t}e^{ik\cdot z_j}=\langle v_t,k\rangle e^{ik\cdot z_j}$，所以 $\Phi_t$ 的列就是 $\{\langle v_t,k\rangle a_k\}$。**导数行是函数值行的标量倍——这是相干性分析可以完全避开概率论的原因。**

**精确的相干性比值。** 直接计算给出

$$
\lambda=\max_{k\ne k'}\frac{\bigl|1+\sum_{t=1}^{r}\langle v_t,k\rangle\langle v_t,k'\rangle\bigr|}
{\sqrt{1+\sum_t\langle v_t,k\rangle^2}\ \sqrt{1+\sum_t\langle v_t,k'\rangle^2}},
$$

即 $\mu(\tilde\Phi)=\lambda\,\mu(\Phi)$ 的比值恰是增广向量之间夹角的余弦上确界。

**可容许性——论文引入的新假设。** 记 $V:=(v_1,\dots,v_r)\in\mathbb R^{d\times r}$，定义距离

$$
\|k-k'\|_V:=\|V^\top k-V^\top k'\|_\infty .
$$

称方向 $v_1,\dots,v_r$ 关于 $\Gamma$ **可容许**，若对 $\Gamma$ 中一切 $k\ne k'$ 有 $\|k-k'\|_V\ne0$（等价地 $V^\top k\ne V^\top k'$）。若 $\mathrm{span}\{v_1,\dots,v_r\}=\mathbb R^d$ 则 $V$ 可容许；但可容许**严格更弱**：论文举 $d=2$、$\Gamma=[-q,q]^2\cap\mathbb Z^2$、$v_1=[1,\sqrt2]^\top$，此时 $\{\eta\in\mathbb R^2:\langle\eta,v_1\rangle=0\}\cap\mathbb Z^2=\{0\}$，**单个方向 $v_1$ 就已可容许，尽管它张不成 $\mathbb R^2$**。再记

$$
\Gamma_{\min}:=\min_{k\ne k',\ k,k'\in\Gamma}\|k-k'\|_V^2 .
$$

**精细化的 Cauchy–Schwarz（引理 2.1，引自文献）。** 若 $x,y\in\mathbb R^d$ 不成比例，且 $u$ 满足 $\langle u,x\rangle=0$、$\langle u,y\rangle=1$，则

$$
\frac{\langle x,y\rangle^2}{\|x\|^2\|y\|^2}\le1-\frac{1}{\|y\|^2\|u\|^2}.
$$

把它用到增广向量 $x=(1,V^\top k)$ 与 $y=(1,V^\top k')$ 上：可容许性保证二者不成比例，因而存在满足两个约束的 $u$；再把 $\|y\|^2$ 与 $\|u\|^2$ 在 $\Gamma$ 上一致地放大，就得到定理 2.1 中分母上那两个极大值。（本站掌握的材料记录了引理的陈述与「用它去做」这一技术，但**没有记录 $u$ 的显式构造**，所以上面这一步是对论证路径的复述而非逐字重现。）

**加权变体。** 令 $W$ 为对角阵，值行的权取 1，**所有导数行统一乘以同一个常数 $\alpha$**。这是一个只有一个自由参数的极简预条件，却足以把相干性界在 $q$ 上的阶数从 $O(q^4)$ 改善到 $O(q)$（见下）。

**问题 2 的列归一化。** 置 $\Psi=\tilde\Psi W$，其中 $\tilde\Psi:=\bigl(D^{\tau_j}_{v_j}\varphi_t(z_j)\bigr)\in\mathbb C^{N\times M}$，

$$
W:=\mathrm{diag}\Bigl(\Bigl(\textstyle\sum_{j=1}^N|\langle v_j,k\rangle|^{2\tau_j}\Bigr)^{-1/2}:k\in\Gamma\Bigr)\in\mathbb C^{M\times M}.
$$

解 $\arg\min\|c\|_1$（约束 $\Psi c=f$）后再还原 $c^{\#}=Wc^{\#}$。辅助量为

$$
Z_k:=\{j:\langle v_j,k\rangle^{\tau_j}=0,\ 1\le j\le N\},
\qquad
\kappa:=\max_{k\in\Gamma}\#Z_k,
\qquad
R_0:=\frac{\max_{j\in Z_k^c}|\langle v_j,k\rangle|^{\tau_j}}{\min_{j\in Z_k^c}|\langle v_j,k\rangle|^{\tau_j}},
$$

约定 $\tau_j=0$ 时 $\langle v_j,k\rangle^{\tau_j}=1$。$Z_k$ 收集那些「对频率 $k$ 完全不提供信息」的测量，$R_0$ 度量各测量对同一频率的灵敏度差异——**它就是后面样本复杂度里那个可能很大的常数**。

> [!warning] 两处常数在本站所用的全文提取中排版有歧义
> 其一，上式 $W$ 的对角元在原文中印作 $1/\sum_j|\langle v_j,k\rangle|^{2\tau_j}$，根号覆盖到哪里无法从提取结果判断；由于它的既定用途是把 $\Psi$ 的**列范数**归一化，必须是**负一半次幂**，本页据此写出。其二，下面定理 2.2 的最优 $\alpha$ 与 $\Gamma=[-q,q]^d$ 的具体取值同样存在根号范围的歧义。**这两处都是重构，不要直接引用其常数形式，请核对排版原文。**

### 定理

**定理 1.1（引用的恢复结果）。** 若 $\Phi\in\mathbb R^{N\times M}$ 满足 $s$ 阶限制等距性质且 $\delta_s<\frac13$，则 $\ell_1$ 极小化恢复成功。论文注明较早的条件是 $\delta_{3s}+3\delta_{4s}<2$，$\delta_s<\frac13$ 是近期的改进。

**定理 2.1（主相干性定理）。** 若 $V=(v_1,\dots,v_r)\in\mathbb R^{d\times r}$ 关于 $\Gamma$ 可容许，则

$$
\mu(\tilde\Phi)\le\lambda\cdot\mu(\Phi),
\qquad
\lambda\le\left(1-\frac{\Gamma_{\min}}
{\max_{k\in\Gamma}\bigl(1+\|V^\top k\|^2_\infty\bigr)\ \max_{k\in\Gamma}\bigl(1+\|V^\top k\|^2_2\bigr)}\right)^{1/2}.
$$

**这是确定性的**：对任意采样点成立，不带任何概率前提。

**推论 2.1。** 取 $\Gamma=[-q,q]^d\cap\mathbb Z^d$ 与 $V=(e_1,\dots,e_d)$（完整梯度），则 $\Gamma_{\min}=1$、$\max_k\|V^\top k\|_\infty^2=q^2$、$\max_k\|V^\top k\|_2^2=dq^2$，故

$$
\mu(\tilde\Phi)\le\left(1-\frac{1}{(1+q^2)(1+dq^2)}\right)^{1/2}\mu(\Phi).
$$

**命题 2.1（可容许性不是装饰）。** 取 $d=2$、$\Gamma=[-q,q]^2\cap\mathbb Z^2$，以及确定性的二次相位样本 $\Xi=\{2\pi(j,j^2)/p:j=0,\dots,p-1\}$（$p>2q+1$ 为素数）。Gauss 和公式给出 $\mu(\Psi)=1/\sqrt p$；追加 $e_1$ 方向的导数块后 $\mu(\tilde\Psi)=\mu(\Psi)$ **精确相等，毫无改进**——因为 $V=(e_1)$ 对二维指标集不可容许。**这条命题说明可容许性是真正的假设，不满足它时梯度信息在相干性意义上一无所获。**

**定理 2.2 与推论 2.2（加权严格改进）。** 在定理 2.1 的假设下，用上面的 $W$ 得到

$$
\mu(W\tilde\Phi)\le\lambda\mu(\Phi),
\qquad
\lambda\le\left(1-\frac{\alpha^2\Gamma_{\min}}
{\max_k\bigl(1+\alpha^2\|V^\top k\|_\infty^2\bigr)\ \max_k\bigl(1+\alpha^2\|V^\top k\|_2^2\bigr)}\right)^{1/2},
$$

在 $\alpha=\bigl(\max_{k\in\Gamma}\|V^\top k\|_\infty\cdot\max_{k\in\Gamma}\|V^\top k\|\bigr)^{-1/2}$ 处最优（此式即上文警告中的第二处歧义）。对 $\Gamma=[-q,q]^d$、$V=(e_1,\dots,e_d)$、$\alpha=1/(\sqrt d\,q)$，

$$
\mu(W\tilde\Phi)\le\left(1-\frac{1}{(1+\sqrt d)^2q}\right)^{1/2}\mu(\tilde\Phi),
$$

且注 2.1 验证了 $\bigl(1-\frac{1}{(1+\sqrt d)^2q}\bigr)^{1/2}\le\bigl(1-\frac{1}{(1+q^2)(1+dq^2)}\bigr)^{1/2}$，即加权严格改进推论 2.1。**改进的实质在 $q$ 的阶：分母里的 $O(q^4)$ 变成 $O(q)$。**

**定理 3.1（问题 2 的限制等距性质）。** 设 $z_1,\dots,z_N$ 独立同分布均匀取自 $[-\pi,\pi)^d$，且 $v_j\in\mathbb R^d$、$\tau_j\in\mathbb Z_{\ge0}$ 满足对每个 $k\in\Gamma$ 都有 $\sum_{j=1}^N\langle v_j,k\rangle^{2\tau_j}\ne0$。对 $0<\delta\le\frac12$ 与 $0<\epsilon<1$，若

$$
N\ \ge\ \frac{2(C_0R_0)^2}{\delta^2}\,s\,(\ln 100s)^2\,\ln(4M)\,\ln(10N)\,\ln\frac{\beta}{\epsilon}\ +\ \kappa,
$$

则 $\mathbb P(\delta_s\le\delta)\ge1-\epsilon$，其中 $\delta_s$ 是 $\Psi$ 的 $s$ 阶限制等距常数，$M=\#\Gamma$，$C_0,\beta$ 为普适常数。结合定理 1.1（取 $\delta=1/3$）即得 $\ell_1$ 恢复成功的概率 $\ge1-\epsilon$。**作者给出的解读是：导数测量与函数值测量起同样的作用，差别只在常数 $R_0$——而论文自己承认 $R_0$「可能很大」。** 证明沿 Rudelson–Vershynin / Rauhut 路线：对称化引理、带显式常数（文献值 $C_1'=94.81$、$C_2'\approx82.56$、$\beta=6.028$）的 Rademacher 混沌界、以及矩到尾概率的引理。

**推论 3.1。** 取 $\Gamma=[1,q]^d\cap\mathbb Z^d$（$q\ge2$）、$\tau_j\in\{0,1\}$、$v_j\in\{e_1,\dots,e_d\}$，样本均匀取自 $[-\pi,\pi)^d$，则 $\kappa=0$、$R_0=q$，样本条件化为

$$
N\ \ge\ \frac{2(C_0q)^2}{\delta^2}\,s\,(\ln100s)^2\,\ln(4M)\,\ln(10N)\,\ln\frac\beta\epsilon
\ \Longrightarrow\ \mathbb P(\delta_s\le\delta)\ge1-\epsilon .
$$

**论文记录的对照样本数（均为引用，不是新结果）。** Rauhut：均匀样本下 $N/\log N\ge C\delta^{-2}s\log^2(s)\log(M)\log(\epsilon^{-1})$ 给出 $\delta_s\le\delta$，概率 $\ge1-\epsilon$。Kunis–Rauhut：$N\ge C(2s-1)^2\log(4M^2/\epsilon)$ 给出 $\mu<\frac{1}{2s-1}$，概率 $\ge1-\epsilon$。Xu（确定性）：对 $\Gamma=[-q,q]^d$、$d\ge2$，点 $z_j=2\pi(j,j^2,\dots,j^d)/N$（$N>\max\{(2s-1)^2(d-1)^2,\ 2q+1\}$ 为素数）给出 $\mu(\Phi)<\frac{1}{2s-1}$。

> [!note] 论文明确留下的猜想
> 导数信息是否也**降低限制等距常数**（而不只是相干性），论文列为**猜想**：数值实验给出肯定的迹象，但没有证明，对象是 $[-\pi,\pi)^d$ 上的均匀样本。方向 $\{v_j\}$ 与权 $\alpha$ 的最优选取同样留作未来工作；从 Fourier 基推广到正交多项式基这一条，正是编号 32 所做的事。

### 数值实验

$\ell_1$ 用 SPGL1 求解，全部实验取 $\Gamma=[-q,q]^d\cap\mathbb Z^d$；图例中的 "standard" 与 "gradient-enhanced" 分别指标准方法与梯度增强方法，百分比表示所用导数分量的比例。

| 例  | $d$ | $q$ | 具体设置                                                                                                                                                | 报告的量                                      |
| --- | --- | --- | ------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------- |
| 4.1 | 2   | 10  | $v_1=e_1$，$v_2=e_2$；均匀样本；支撑在 $s$ 元子集上均匀抽取，非零系数独立同分布标准正态；每种配置 500 次试验                                            | 恢复率对 $N$（$s=5$）；恢复率对 $s$（$N=20$） |
| 4.2 | 5   | 2   | 测试 20% 与 40% 增强（一个与两个变量的偏导）                                                                                                            | 恢复率对 $N$（$s=6$）；恢复率对 $s$（$N=30$） |
| 4.3 | 2   | 5   | 问题 2 设定，$s=8$，$M=40$；用 $N/4$ 个函数值与 $3N/4$ 个 $D_{e_1}f$ 值，对照用 $N$ 个函数值的标准 $\ell_1$                                             | 恢复率                                        |
| 4.4 | 2   | 10  | 函数逼近：$g(x)=\sin\bigl(\cos\sum_jx_j\bigr)$ 与 $g(x)=\cos\bigl(\sum_jx_j\bigr)\exp\bigl(\sin(\sum_jx_j)\bigr)$；约束为在样本点匹配 $f$ 与 $D_{e_t}f$ | 离散 $L^2$ 误差对 $N$                         |

作者报告的结论：例 4.1 与 4.2 中梯度信息提高恢复率，且**随所含导数比例单调改善**；例 4.3 中**导数值起到与函数值类似的作用**——这是定理 3.1 的经验对应；例 4.4 中梯度信息「显著地」（"dramatically"）改善精度。

> [!warning] 实验数据与一处对不上的配置
> 本站没有转录任何恢复率或误差的具体数值，上表只记录配置与作者陈述的定性结论。另外，例 4.3 记录的 $q=5$、$d=2$ 与 $M=40$ 在论文自己的约定 $M=\#\Gamma$ 下对不上（$\#([-5,5]^2\cap\mathbb Z^2)=121$），**这一处配置请以排版原文为准**。
>
> 与理论的落差同样值得说明：定理 2.1、2.2 讲的是相干性**比值**的改善，而例 4.1–4.2 直接测的是恢复率；两者之间还隔着 $\mu<\frac{1}{2s-1}$ 这一步，实验用的 $N=20$、$N=30$ 远达不到该条件所需的规模。**实验支持的是「梯度有用」这一论断，不是定理的紧性。**

### 与其他论文的关系

编号 29 是编号 **32** 的三角 / Fourier 孪生篇，同年由重叠的作者群发表；编号 32 把同一想法推广到带预条件的正交多项式混沌，而这恰是编号 29 自己列出的未来工作。两篇都走互相干性路线（$\mu<1/(2s-1)$），与编号 **10** 相同。编号 29 的问题 2 表述——函数值与导数值在**不同位置**混合给出——在编号 32 的第 4 节中被原样列为公开问题。它引用的确定性二次相位样本 $z_j=2\pi(j,j^2,\dots,j^d)/N$ 正是[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 9、14]] 与本页编号 10 使用的 **Weil 点**；命题 2.1 中 $\mu(\Psi)=1/\sqrt p$ 的 Gauss 和计算与编号 10 的相干性估计同出一源。

## 32：多项式混沌稀疏逼近的梯度增强 ℓ1 极小化

（本节记号沿用论文：$N$ 个采样点，$M=\binom{d+n}{n}$ 个基函数，$N\ll M$。）

### 直觉

当只有 $N\ll M$ 次函数求值可用时，经 $\ell_1$ 极小化做稀疏多项式混沌恢复是标准工具。梯度 $\partial_kf$ 常常可以廉价获得（例如伴随求解），把它们接到测量矩阵上使行数乘以 $d+1$ 而额外代价很小——**但朴素堆叠会破坏均值各向同性** $\mathbb E[\frac1N\Phi^\top\Phi]=I$，而这正是压缩感知保证的基石；论文自己的数值实验显示，朴素堆叠**实际上恶化了**互相干性常数。

问题的根源在于多项式基与 Fourier 基的一处关键差异：$\frac{d}{dx}p^{(\alpha,\beta)}_n=c(n,\alpha,\beta)\,p^{(\alpha+1,\beta+1)}_{n-1}$，**导数是另一族 Jacobi 多项式、带另一个正交权，而且缩放常数 $c(n,\alpha,\beta)$ 随次数增长**。所以梯度行既不是值行的倍数，也不在同一个测度下正交归一。论文的贡献就是一套通用配方：用**行预条件** $W$ 把每一块行放回它自己的正交权下，用**列归一化** $P$ 吸收 $c(n,\alpha,\beta)$ 的增长，使复合矩阵重新各向同性。配方对有界与无界域都适用。

### 问题设定

$x=(x_1,\dots,x_d)^\top$ 分量独立，边缘密度 $\rho_i$ 于 $\Gamma_i$，$\rho(x)=\prod_i\rho_i(x_i)$ 于 $\Gamma=\otimes_i\Gamma_i$；一维正交归一基 $\phi^i_n$ 满足 $\int_{\Gamma_i}\phi^i_n\phi^i_\ell\rho_i=\delta_{n,\ell}$，多维基 $\psi_n(x)=\prod_{i=1}^d\phi^i_{n_i}(x_i)$，$\mathbb E[\psi_n\psi_j]=\delta_{n,j}$。全次数指标集 $\Lambda^T_n=\{k\in\mathbb N_0^d:\sum_ik_i\le n\}$，$M=\binom{d+n}{n}$，展开 $f_n=\sum_{j=1}^Mc_j\psi_j(x)$。

标准 $\ell_1$：样本 $\Xi=\{z^{(1)},\dots,z^{(N)}\}\subset\Gamma$，$[\Phi]_{ij}=\psi_j(z^{(i)})$，$\Phi\in\mathbb R^{N\times M}$，解 $\arg\min\|c\|_1$ 使 $\Phi c=f$（或 $\|\Phi c-f\|_2\le\epsilon$）。互相干性 $\mu(\Phi)$ 与 $\mu<\frac{1}{2s-1}$ 的恢复判据如本页开头所述。

### 推导

**梯度增强问题（论文的中心对象）。**

$$
\arg\min_{c\in\mathbb R^M}\|c\|_1
\quad\text{s.t.}\quad
W\tilde\Phi P\,c=W\tilde f,
$$

$$
\tilde f=\begin{pmatrix}f\\ f^\partial\end{pmatrix},\quad
\tilde\Phi=\begin{pmatrix}\Phi\\ \Phi^\partial\end{pmatrix},\quad
\Phi^\partial=\begin{pmatrix}\frac{\partial\Phi}{\partial x_1}\\ \vdots\\ \frac{\partial\Phi}{\partial x_d}\end{pmatrix},\quad
f^\partial=\begin{pmatrix}\frac{\partial f}{\partial x_1}\\ \vdots\\ \frac{\partial f}{\partial x_d}\end{pmatrix},
$$

其中 $\bigl[\frac{\partial\Phi}{\partial x_k}\bigr]_{ij}=\frac{\partial\psi_j}{\partial x_k}(z_i)$，故 $\tilde\Phi\in\mathbb R^{N(d+1)\times M}$。**两个新矩阵各司其职：$W$ 是行预条件，形式取决于混沌族与 $\Xi$ 的抽取方式；$P$ 是列归一化，保证复合矩阵 $\hat\Phi:=W\tilde\Phi P$ 均值各向同性。**

**Legendre 基配 Chebyshev 采样。** 取乘积 Chebyshev 密度 $\rho_c(x)=\prod_{j=1}^d\frac{1}{\pi\sqrt{1-x_j^2}}$。利用 Legendre 多项式的导数关于 $\eta(x)=(1-x^2)$ 正交这一经典事实，论文推出

$$
\mathbb E_c\Bigl[\frac{2^{-d}}{\rho_c(z)}\psi_i(z)\psi_j(z)
+\sum_{k=1}^d\frac{1-z_k^2}{\rho_c(z)}\frac{\partial\psi_i}{\partial x_k}(z)\frac{\partial\psi_j}{\partial x_k}(z)\Bigr]
=\delta_{ij}\Bigl(1+\sum_{k=1}^d c_k\,i_k(i_k+1)\Bigr).
$$

**这一式是全套设计的模板**：左边两项分别是值行与梯度行各自的正交权（$2^{-d}$ 对应均匀密度，$1-z_k^2$ 对应导数族的权），右边的 $\delta_{ij}$ 说明加权后仍然正交，而括号里的因子就是要用 $P$ 除掉的那个随次数增长的量。相应地

$$
W^0_{n,n}=\Bigl(\tfrac{4}{\pi^2}\bigl(1-(z^{(n)}_j)^2\bigr)\Bigr)^{d/4},
\qquad
W^j_{n,n}=\frac{W^0_{n,n}}{\sqrt2}\bigl(1-(z^{(n)}_j)^2\bigr)^{1/2},
\qquad
P_{i,i}=\Bigl(1+\sum_{k=1}^d c_k\,i_k(i_k+1)\Bigr)^{-1/2},
$$

$W$ 为分块对角 $(W^0,W^1,\dots,W^d)$，结果是 $\mathbb E_c\bigl[\frac1N\hat\Phi^\top\hat\Phi\bigr]=I$。

> [!warning] Legendre 情形的两个常数不要照抄
> 上面第一式的首项系数在本站所用的全文提取中印作 "$2d$"，显然丢了指数；由下面的一般 Jacobi 公式可知它必须是密度比 $\rho/\rho_c$ 中的 $\rho\equiv2^{-d}$（均匀），故本页写作 $2^{-d}$——**这是重构**。$W^0_{n,n}$ 的表达式里只出现第 $j$ 个分量却带 $d/4$ 次幂，强烈提示提取过程丢掉了一个 $\prod_{j=1}^d(\cdot)^{1/4}$。**不要依赖这两处印出的常数**；一般 Jacobi 版本（下一段）是干净且自洽的，应以它为准。

**一般 Jacobi 配 Chebyshev 采样——干净的通用形式。** 一维 Beta / Jacobi 密度

$$
\rho^{(\alpha,\beta)}(x)=d^{(\alpha,\beta)}(1-x)^\alpha(1+x)^\beta,
\quad \alpha,\beta\ge-\tfrac12,
\qquad
d^{(\alpha,\beta)}=\frac{\Gamma(\alpha+\beta+2)}{\Gamma(\beta+1)\Gamma(\alpha+1)2^{\alpha+\beta+1}},
$$

$\rho_c\equiv\rho^{(-1/2,-1/2)}$，多维取乘积。则对 $z\sim\rho_c$，

$$
\mathbb E\Bigl[\frac{\rho^{(\alpha,\beta)}(z)}{\rho_c(z)}\psi_i\psi_j
+\sum_{k=1}^d\frac{\rho^{(\alpha+e_k,\beta+e_k)}(z)}{\rho_c(z)}
\frac{\partial\psi_i}{\partial x_k}\frac{\partial\psi_j}{\partial x_k}\Bigr]
=\delta_{ij}\Bigl(1+\sum_{k=1}^d c^2(i_k,\alpha_k,\beta_k)\Bigr),
$$

$e_j$ 为第 $j$ 个基本向量、$e_0=0$，归一化常数显式为

$$
c^2(i_k,\alpha_k,\beta_k)=\frac{i_k(i_k+\alpha_k+\beta_k+1)(\alpha_k+\beta_k+2)(\alpha_k+\beta_k+3)}{4(\alpha_k+1)(\beta_k+1)} .
$$

**注意 $c^2$ 关于次数 $i_k$ 是线性增长的——这就是不加 $P$ 时梯度行会压过值行的定量原因。** 于是

$$
W^0_{n,n}=\sqrt{\frac{\rho^{(\alpha,\beta)}(z^{(n)})}{\rho_c(z^{(n)})}},
\qquad
W^j_{n,n}=\sqrt{\frac{\rho^{(\alpha+e_j,\beta+e_j)}(z^{(n)})}{\rho_c(z^{(n)})}},
\qquad
P_{i,i}=\Bigl(1+\sum_{k=1}^d c^2(i_k,\alpha_k,\beta_k)\Bigr)^{-1/2},
$$

给出 $\mathbb E\bigl[\frac1N\hat\Phi^\top\hat\Phi\bigr]=I$。

**论文陈述的通用配方（三步）。** (i) 从对该混沌基**次数渐近良好**的采样测度抽取 $\Xi$；(ii) 设计 $W$ 使多项式混沌基在该测度下均值各向同性；(iii) 选取 $P$ 使**梯度行**的各向同性也被保持。

**相干性参数。** 标准方法用 $\mu_L(\Phi):=\sup_{i,\,z\in\Xi}|\Phi_i(z)|_2^2$，梯度增强方法用

$$
\beta_L(\hat\Phi):=\sup_{i,\,z\in\Xi}\bigl\|\hat\Phi_i(z)\bigr\|_2,
\qquad
\hat\Phi_i(z)=\frac{1}{P_{i,i}}
\begin{pmatrix}
\sqrt{\rho^{(\alpha,\beta)}(z)/\rho_c(z)}\ \Phi_i(z)\\
\sqrt{\rho^{(\alpha+e_1,\beta+e_1)}(z)/\rho_c(z)}\ \partial_{x_1}\Phi_i(z)\\
\vdots\\
\sqrt{\rho^{(\alpha+e_d,\beta+e_d)}(z)/\rho_c(z)}\ \partial_{x_d}\Phi_i(z)
\end{pmatrix}.
$$

**Hermite 配 Gauss 采样。** Hermite 多项式的导数在**同一个** Gauss 测度下正交，因此这一情形干净得多：对适当归一化的 $\psi_j$ 与标准多维正态 $z$，

$$
\mathbb E\Bigl(\psi_i(z)\psi_j(z)+\sum_{k=1}^d\frac{\partial\psi_i}{\partial x_k}(z)\frac{\partial\psi_j}{\partial x_k}(z)\Bigr)
=\delta_{ij}\Bigl(1+\sum_{k=1}^d i_k\Bigr),
$$

于是取 $P_{i,i}=\bigl(1+\sum_ki_k\bigr)^{-1/2}$ 而 **$W=I$**——**无界域上根本不需要行预条件**。（论文注明这一情形与其参考文献 [35] 中的一个类似结果相当。）

**一条结构性注记。** $P$ 之所以能取成对角阵，唯一的原因是经典族的导数仍然是正交族；论文指出 **Jacobi、Laguerre、Hermite 是仅有的具备此性质的一维族**。因此对非经典的混沌基，$P$ 必须取成多项式导数的 Gram 矩阵的任一逆平方根，**不再是对角阵**。

### 定理

**定理 3.1（主定理）。** 设 $\Phi$ 与 $\hat\Phi$ 分别是 Jacobi 展开配 Chebyshev 采样下标准与梯度增强 $\ell_1$ 的设计矩阵，则

$$
\mu_L(\Phi)\le\prod_{j=1}^d 2e\Bigl(2+\sqrt{\alpha_j^2+\beta_j^2}\Bigr),
\qquad
\beta_L(\hat\Phi)\le C\prod_{j=1}^d 2e\Bigl(2+\sqrt{\alpha_j^2+\beta_j^2}\Bigr),
$$

其中

$$
1\le C\le 1+\frac{\sqrt2}{2}\approx1.707 .
$$

下界 $C=1$ 在 $\alpha_k=\beta_k=-\frac12$（全部 $k$）时取到，上界在存在 $k$ 使 $\alpha_k=\beta_k=0$ 时出现。此外 $\mathcal N(\hat\Phi)\subset\mathcal N(\Phi)$，且**当 $\Phi$ 欠采样时几乎必然是真子集**——这就是「加梯度不会更差、且一般会更好」的精确含义。

**支撑引理 A.1（引自 Rauhut–Ward 型的 Jacobi 一致界）。** 对一切 $\alpha\ge-\frac12$、$\beta\ge-\frac12$ 的 Jacobi 权，

$$
\sup_{x\in[-1,1]}\frac{\rho^{(\alpha,\beta)}(x)}{\rho_c(x)}\bigl[p^{(\alpha,\beta)}_n(x)\bigr]^2
\le 2e\Bigl(2+\sqrt{\alpha^2+\beta^2}\Bigr),
$$

关于 $n,\alpha,\beta$ 一致。定理 3.1 的证明就是把这条界用到两行结构 $\bigl(p_n^{(\alpha,\beta)}(z),\ \frac{d}{dx}p_n^{(\alpha,\beta)}(z)\bigr)=\bigl(p_n^{(\alpha,\beta)}(z),\ c(n,\alpha,\beta)p_{n-1}^{(\alpha+1,\beta+1)}(z)\bigr)$ 上。

> [!warning] 定理没有做到什么——作者自己讲得很清楚
> 论文写道：「理想情况下我们希望证明梯度方法给出更小的参数 $\beta_L$，即 $\beta_L(\hat\Phi)\le\mu_L(\Phi)$……我们的分析没有结出这个果实。」实际证到的只有两条：(i) $\Phi$ 与 $\hat\Phi$ 的相干性界都是「一个常数的 $d$ 次幂」，且**与多项式次数无关**；(ii) 多出的因子 $C\le1+\frac{\sqrt2}{2}$ **与维数无关且很小**。
>
> 还有一条更重要：**本文没有 $N\gtrsim s\log^k$ 形式的样本复杂度定理**，保证走的是从文献引用的 $\mu<\frac{1}{2s-1}$ 判据。任何关于本文方法样本复杂度指数的说法都没有出处。

**论文提出的公开方向。** 部分 / 方向梯度数据 $D_{v_t}f(z^{(j)})=\langle\nabla f,v_t\rangle$（并坦承「实际中如何获得这类方向导数并不清楚」）；以及更高阶的方向数据 $D^{\tau_j}_{v_j}f(z_j)=y_j$，其中函数值位置与导数位置**不必重合**（这恰是编号 29 的问题 2）。论文还明确表示它的预条件**不宣称最优**，并建议[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 22 / 24]] 的 Christoffel 加权方案是次数渐近意义下更好的替代。

### 数值实验

$\ell_1$ 用 SPGL1 求解。三个方法标签贯穿全部实验：_standard_（在 $\Phi$ 上做 $\ell_1$）、_gradient-enhanced_（在 $\hat\Phi$ 上）、以及 _**standard-double**_——给标准方法**同样多的标量数据总预算**：$d=2$ 时 100% 梯度方法用掉 $3N$ 条信息（$N$ 个函数值加 $2N$ 个导数），故 standard-double 用 $3N$ 个函数值。百分比表示所用导数分量的比例（$d=2$ 的 50% 指 $N$ 个函数值加 $N$ 个随机挑选的偏导）。

| 小节          | 配置                                                                                                                                           | 报告的量                                                                      |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- |
| §5.1 稳定性   | $(d,n)=(2,30)$ 与 $(6,5)$                                                                                                                      | $\Phi$、$\tilde\Phi$、$\hat\Phi$ 的互相干性对样本数；以及固定 $N=80$ 时对 $M$ |
| §5.2 定稀疏度 | $(d,n)=(2,20)$；均匀输入、Chebyshev 样本；恰 $s$-稀疏、非零系数独立同分布标准正态；100 次试验；判定 $\lVert c-\tilde c\rVert_\infty\le10^{-3}$ | 恢复概率对 $N$（$s=8$）；恢复概率对 $s$（固定 $N$）                           |
| §5.2 定稀疏度 | $(d,n)=(10,3)$；测试 10% 与 20% 梯度增强（一个与两个偏导）                                                                                     | 恢复概率对 $N$（$s=6$）；恢复概率对 $s$（$N=70$）                             |
| §5.3 函数逼近 | $(d,n,M)=(2,20,231)$、$(10,3,286)$ 与 $(6,5)$；Legendre 混沌配 Chebyshev 样本                                                                  | 离散 $L^2$ / 均方根误差对 $N$                                                 |

§5.3 的三个测试函数是球面型 $f_1(x)=\sum_{i=1}^dx_i^2$、Gauss 型 $f_2(x)=\exp\bigl(-\sum_{i=1}^d0.01(\tfrac12(x_i+1)-0.375)^2\bigr)$、以及正弦型 $f_3(x)=\sum_{i=1}^d0.3+\sin(\tfrac{16}{15}x_i-0.7)+\sin^2(\tfrac{16}{15}x_i-0.7)$。

**§5.1 是全篇最重要的实验，它才是整套预条件装置的存在理由**：加了预条件的 $\hat\Phi$ 互相干性远好于 $\Phi$，而**朴素堆叠导数行得到的 $\tilde\Phi$ 反而破坏了 $\Phi$ 原有的稳定性**。§5.2 与 §5.3 的结论是：含梯度信息提高恢复率，导数信息越多恢复越好；函数逼近精度「显著」改善。

> [!warning] 论文正文与图注在样本数上不一致
> §5.2 中 $(d,n)=(2,20)$ 一组「恢复概率对 $s$」的固定样本数，**正文写 $N=35$，图注写 $N=50$**；$(d,n)=(10,3)$ 一组的对应数字**大部分处写 $N=70$，另有一处写 $N=50$**。本站保留这两处不一致而不选边。上表按较常出现的值填写。
>
> 与理论的落差：定理 3.1 只界定了相干性参数 $\mu_L$ 与 $\beta_L$ 的上界，且明说没有证出 $\beta_L\le\mu_L$；而实验直接测的是恢复概率与逼近误差。**实验建立的是方法有效，不是定理预测了这种有效性。** 本站同样没有转录任何具体的成功率或误差数值。

### 与其他论文的关系

编号 32 是编号 **29** 在多项式混沌 / 多变量情形的对应物：同一作者圈、同一个梯度增强 $\ell_1$ 想法，但编号 29 处理三角多项式与随机频率，编号 32 处理正交多项式混沌与预条件。它直接建立在[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 24]] 的预条件框架之上：Chebyshev 采样配 $\sqrt{\rho^{(\alpha,\beta)}/\rho_c}$ 加权就是编号 24 的 CSA-a，也就是[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 14]] 中转述的 Rauhut–Ward 渐近采样；论文并且明说编号 22 / 24 的 Christoffel 加权方案（次数渐近最优）会是比它所用预条件更好的选择。它走的互相干性路线与编号 **10** 在确定性 Chebyshev 点集上走的是同一条，与编号 **21**、24 的限制等距路线相对。引理 A.1 的一致 Jacobi 界 $2e(2+\sqrt{\alpha^2+\beta^2})$ 与编号 21、24 引用的是同一族界。

## 数据驱动：要估计的不只是系数，还有基

标准多项式混沌假设输入分布已知且有解析形式，因此正交基可以由该分布的正交多项式给出。**数据驱动**情形下，输入分布只通过有限样本给出：没有解析密度，因此没有现成的正交基。

编号 36 与 44 处理这一设定，要点相同：**正交基必须从经验测度构造**——用 Gram–Schmidt，或者更稳定地，用矩的 Hankel 矩阵经 Cholesky 分解读出三项递推系数——而采样设计随之改变，因为诱导采样与 Christoffel 加权都需要一个基，**而基本身现在是估计出来的**。

这引出一个本专题内其他工作不涉及的误差来源：**基的估计误差**。它与逼近误差、采样误差叠加，因此数据驱动情形的分析比已知分布情形多一层。两篇论文对这一层的处理都是坦白承认而非解决：编号 36 明说因为均值与方差是估计出来的而引入了「密度误差」，如何量化与控制留作未来工作。

两篇的分工是：编号 36 用**加权最小二乘**配**平衡测度**采样（渐近最优），编号 44 用**预条件 $\ell_1$** 配**诱导采样**（在每个有限次数下都精确，且只落在数据上）。

## 36：数据驱动多项式混沌展开的加权最小二乘

（本节记号：$M$ 个样本，$N$ 个数据驱动基函数，$k$ 为最高次数。）

### 直觉

经典广义多项式混沌预设输入密度 $\rho$ 已知（偶然不确定性）。很多应用里手上只有输入的**样本或矩**（认知不确定性），而假设一个错的密度会让整个代理模型有偏。数据驱动 / 任意多项式混沌（aPC）直接从矩构造正交基解决了前半个问题，但既有的后处理是一种稀疏网格配点，**其配点本身要由数据驱动基经矩阵运算算出**——于是基一变，设计就得重算。

论文的贡献正是替换这一后处理：让**采样既不依赖输入密度、也不依赖数据驱动基**，因而可以离线完成；同时把稳定性从 $N$ 的二次改进到拟线性。达成这一点的机制是本页反复出现的那个：Christoffel 加权把稳定性因子 $\kappa(N)=\max_\xi\sum_j\Phi_j^2(\xi)$ 压到它的理论下限 $N$。代价则是加权后的基不再关于 $\rho$ 正交，而是关于一个依赖多项式空间的变换测度正交——论文用位势论的极限（平衡测度）把这个依赖抹掉，**因而稳定性只在 $N\to\infty$ 的渐近意义下成立**。

### 问题设定与推导

**矩的 Hankel 矩阵与矩匹配基。** 给定原点矩 $\mu_0,\dots,\mu_{2k}$，

$$
H=\begin{pmatrix}
\mu_0&\mu_1&\cdots&\mu_k\\
\mu_1&\mu_2&\cdots&\mu_{k+1}\\
\vdots&\vdots&\ddots&\vdots\\
\mu_k&\mu_{k+1}&\cdots&\mu_{2k}
\end{pmatrix}.
$$

当矩来自 $M$ 个样本时，要求该样本集在 **Hamburger 意义下是定的**，即 $\det(H)>0$。对 $H$ 作 Cholesky 分解 $H=R^\top R$（$R$ 上三角），由 **Mysovskih 定理**，$R$ 的元素直接给出一个正交多项式系，三项递推为

$$
\eta\,\phi_{j-1}(\eta)=b_{j-1}\phi_{j-2}(\eta)+a_j\phi_{j-1}(\eta)+b_j\phi_j(\eta),\quad j=1,\dots,k,
$$

$$
a_j=\frac{r_{j,j+1}}{r_{j,j}}-\frac{r_{j-1,j}}{r_{j-1,j-1}},
\qquad
b_j=\frac{r_{j+1,j+1}}{r_{j,j}},
\qquad r_{0,0}=1,\ r_{0,1}=0 .
$$

多维基由张量化得到。论文指出另一条路——反解 Vandermonde 系统——在 $k$ 较大时条件数很差。

**加权最小二乘。** 以数据驱动基 $\{\Phi_j\}_{j=1}^N$ 张成 $\mathbb P_N$（全次数型，最高阶 $k$），

$$
f_N:=\arg\min_{p\in\mathbb P_N}\frac1M\sum_{m=1}^M w_m\bigl(p(z_m)-f(z_m)\bigr)^2
\ \Longleftrightarrow\
c=\arg\min_{c\in\mathbb R^N}\bigl\|W^{1/2}Ac-W^{1/2}f\bigr\|_2^2,
$$

$A=[\Phi_j(z_m)]\in\mathbb R^{M\times N}$，$W=\mathrm{diag}(w_1,\dots,w_M)$。权取归一化 Christoffel 函数值：在本页约定 $K(\xi)=\sum_{j=1}^N\Phi_j^2(\xi)$ 下，

$$
w_m=\frac{N}{K(z_m)}=\frac{N}{\sum_{j=1}^N\Phi_j^2(z_m)} .
$$

**五步算法。** (1) 从**平衡测度**的概率密度 $\hat\rho$ 采样（它依赖输入密度 $\rho$，无界情形还依赖最高次数 $k$，记 $\hat\rho_k$）；(2) 在样本处求 $f$；(3) 组装 $M\times N$ 的 Vandermonde 型矩阵 $A$；(4) 由归一化 Christoffel 函数组装对角阵 $W$；(5) 解加权最小二乘。**注意第 (1) 步与第 (3)(4) 步的顺序：采样在前，基在后，这正是「离线采样」的含义。**

**有界域上采样测度是通用的。** 在 $[-1,1]^d$ 上，**无论输入测度是什么（哪怕未知）**，都从张量 Chebyshev 密度

$$
\hat\rho(\xi)\sim\frac{1}{\pi^d\prod_{k=1}^d\sqrt{1-\xi_k^2}}
$$

采样。唯一需要的信息是「随机变量落在有界域内」——论文称之为「有界情形下 Chebyshev 测度是通用的」。采样也很简单：抽均匀的 $u_m$，取 $z_m=\cos(u_m)$。

**无界域：带显式采样器的猜想。** 论文明说这方面已知结果极少，「以下结果是我们的猜想」，仅由数值验证。

> [!warning] 无界域的平衡测度是猜想，不是定理
> 下面两条密度在原文中就被标为**猜想**，本页也只能这样引用。Gauss 情形的极限诱导 / 平衡测度 $C(2-\|\xi\|^2)^{d/2}$ 在本专题的多篇论文中反复出现，**它在任何一处都是猜想，从来不是定理**；有界情形的 Chebyshev 极限则是有定理支撑的。引用时请勿混为一谈。

_（正态，$\mathbb R^d$）_：先从数据估计 $(\hat\mu,\hat\sigma)$，标准化 $\hat\xi=(\xi-\hat\mu)/\hat\sigma$，用

$$
\hat\rho(\xi)=C\bigl(2-\|\xi\|^2\bigr)^{d/2}
$$

（$C$ 为归一化常数），样本按 $\sqrt k$ 放大。具体采样器：算出 $k$；抽独立标准正态的 $y=(y_1,\dots,y_d)$；抽 $[0,1]$ 上的 $\nu\sim\mathrm{Beta}(\alpha=d/2,\ \beta=d/2+1)$；令

$$
z=\frac{y}{\|y\|_2}(2k\nu)^{1/2},
$$

得到半径 $\sqrt{2k}$ 的欧氏球上的样本。

_（指数，$\mathbb R_+^d$）_：

$$
\hat\rho(\xi)=C\sqrt{\frac{\bigl(4-\sum_{i=1}^d\xi_i\bigr)^d}{\prod_{i=1}^d\xi_i}}
$$

（根号覆盖整个分式这一点是从提取结果**重构**的，标此存疑）。采样器：算出 $k$；抽参数为 $\bigl(\frac12,\frac12,\dots,\frac12,\frac d2+1\bigr)$ 的 $(d+1)$ 维 Dirichlet 向量 $y$；丢掉最后一个分量；令 $z=4ky$。注 4.1 承认：对其他无界密度**连猜想都没有**，退路是截断区域再用 Chebyshev 采样，而截断误差论文没有处理。

### 定理

**定理 3.1（引用：矩问题何时确定基）。** 在假设 1（一切矩有限）与假设 2（分布函数 $F_\eta$ 连续）下，若以下**任一条**成立，则矩问题唯一可解，且矩匹配多项式在 $L^2(\Omega,\sigma(\eta),P)$ 中稠密：(1) $F_\eta$ 有紧支撑；(2) $\liminf_{k\to\infty}\frac{\sqrt[2k]{\mu_{2k}}}{2k}<\infty$；(3) 指数可积性，存在 $a>0$ 使 $\int_{\mathbb R}e^{a|x|}F_\eta(dx)<\infty$（等价于矩母函数在原点附近有限）；(4) **Carleman 条件** $\sum_{k=0}^\infty\frac{1}{\sqrt[2k]{\mu_{2k}}}=\infty$；(5) **Lin 条件**：分布有对称、可微、严格正的密度 $f_\eta$，且存在 $x_0>0$ 使 $\int_{-\infty}^\infty\frac{-\log f_\eta(x)}{1+x^2}dx=\infty$ 且 $-\frac{xf'_\eta(x)}{f_\eta(x)}\nearrow\infty$（$x\to\infty$，$x\ge x_0$）。论文指出：对更一般的设定（例如离散分布）理论仍是开放的，**尽管方法照样跑得动**——实验里确实用了二项与 Poisson 输入。

**定理 4.1（引用的最小二乘稳定性）。** 设 $f$ 在 $\mathbb P_N=\mathrm{span}\{\Phi_j\}$ 中逼近、正交密度为 $\rho$、样本 $\{z_m\}_{m=1}^M$ 取自 $\rho$，作**未加权**最小二乘 $c=\arg\min\|Ac-f\|_2^2$，则

$$
\Pr\Bigl\{\|\,\cdot\,-I\|\ge\tfrac12\Bigr\}\le 2M^{-r}
\qquad\text{只要}\qquad
\kappa(N):=\max_\xi\sum_{j=1}^N\Phi_j^2(\xi)\ \le\ \delta\,\frac{M}{\log M}.
$$

> [!warning] 定理 4.1 有两处需要更正
> 其一，范数里的对象在原文印作 "$A$"，但按定理的内容那必须是 **Gram 矩阵 $A^\top A$**（或其缩放版本），不可能是矩形设计矩阵。本页因此把该位置写成占位符。其二，常数 $\delta$ 原文印作 $1-\frac{\log2}{2-2r}$，其分式结构在本站所用的提取中有歧义；Cohen–Davenport–Leviatan 原始文献中对应的常数形如 $\frac{1-\log2}{2+2r}$。**不要在未核对排版原文的情况下引用这个 $\delta$。**

**推论式的动机（本文的关键论证）。** 稳定性要求 $M\gtrsim\kappa(N)$（差一个对数因子）。对 Legendre 多项式 $\kappa(N)\sim N^2$，于是 $M\ge CN^2$——**不能接受**。引入 $W$ 等价于改用重标基

$$
\hat{\mathbb P}_N=\mathrm{span}\Bigl\{\hat\Phi_j=\sqrt{\tfrac{N}{K(\xi)}}\,\Phi_j\ \Big|\ 1\le j\le N\Bigr\},
\qquad
\hat\kappa(N):=\max_\xi\sum_{j=1}^N\hat\Phi_j^2(\xi)\equiv N,
$$

即**把稳定性因子压到它的理论下限**，从而得到拟线性的 $M\gtrsim N\log N$。但重标基不再关于 $\rho$ 正交，而是关于变换测度

$$
\tilde\rho(\xi)\ \propto\ K(\xi)\rho(\xi)=\Bigl(\sum_{j=1}^N\Phi_j^2(\xi)\Bigr)\rho(\xi)
$$

正交——**这恰恰就是诱导测度**。它依赖多项式空间，抽样并不平凡。论文的化解办法来自位势论：

$$
\tilde\rho(\xi)\ \longrightarrow\ \hat\rho(\xi)\qquad (N\to\infty),
$$

于是改从平衡测度 $\hat\rho$ 采样，**稳定性只在 $N\to\infty$ 的渐近意义下成立**；实际收益是设计变得与多项式空间无关，这对自适应方案很有价值。

> [!warning] 原文 (4.5) 与 (4.7) 中 $K$ 的用法与 (4.3) 冲突
> 论文用符号 $K$ 表示 $N/\sum_j\Phi_j^2$（归一化 Christoffel 函数），(4.3) 的权 $w_m=N/\sum_j\Phi_j^2(z_m)$ 与之一致。但 (4.5) 印作 $\hat\Phi_j=\Phi_j/\sqrt{K}$、(4.7) 印作 $\tilde\rho\sim K\rho=N\rho/\sum_j\Phi_j^2$，**这两式若按同一个 $K$ 去读，会给出 $\hat\kappa\ne N$，也不会给出诱导测度**——它们只有在把 $K$ 读作 $\sum_j\Phi_j^2$（即倒数）时才自洽。本页统一采用站内约定 $K=\sum_j\Phi_j^2$ 写出自洽的形式（上面两式），并在此标明与原文印刷形式的差异。判据很简单：正确的重标必须满足 $\hat\Phi_j=\sqrt{w_m}\,\Phi_j$，与 (4.3) 的权一致。

**本文自己证了什么。** 定理 3.1 与定理 4.1 **都是引用他人结果**。论文自身的贡献是**组合**（aPC 基 + Christoffel 加权的平衡测度采样）、显式采样器、以及数值实验。**没有新的收敛定理，也没有对矩估计误差的分析**：论文明说因为均值与方差是估计出来的而引入了「密度误差」，「如何量化与控制这类误差是我们未来的工作」。

### 数值实验

全部结果对 100 次独立试验取平均；全程比较两种采样率：$M=CN$（线性）与 $M=CN\log N$（对数线性）。用到的输入分布包括映射到 $[-1,1]$ 的离散二项 $\mathrm{Bino}(n,p)$、$[-1,1]$ 上的离散 Poisson $\mathrm{Pois}(\lambda)$、均匀 $U[a,b]$、$(0,\infty)$ 上的指数 $\mathrm{Exp}(\mu)$，以及正态 $N(\mu,\sigma)$。

**§5.1 稳定性。** 报告 $\hat A=W^{1/2}A$ 的条件数（原文印作 "$W^2A$"，**标为提取造成的排版讹误**），以均值配 20% 与 80% 分位数的形式对多项式次数作图。四个 $d=2$ 的测试类型为：

| 情形 | $\xi_1$                 | $\xi_2$             |
| ---- | ----------------------- | ------------------- |
| 1    | $\mathrm{Bino}(20,1/2)$ | $U[-0.6,0.6]$       |
| 2    | $U[-0.8,0.8]$           | $U[-1,1]$           |
| 3    | $\mathrm{Bino}(20,1/2)$ | $\mathrm{Pois}(10)$ |
| 4    | $U[-0.6,0.6]$           | $N(0.1,1.2)$        |

**情形 4 混合了有界与无界的边缘分布，因此两个维度用的是不同的平衡测度**——这一点本身就是「采样器按维分别选取」这一设计的演示。测试的采样率为 $M=1.5N$、$2N$、$N\log N$、$1.5N\log N$；另有一组五维测试用例。

**§5.2 精度**，用 $L$ 个参考样本上的离散 $\ell^2$ 误差度量（原文印出的前置因子是 "$L\sum_l$"，按定义必须是 $\frac1L\sum_l$，**标为讹误**）。

- §5.2.1 解析测试函数：$f_1(\xi)=\exp(\sum_k\xi_k)$；$f_2(\xi)=\sum_k0.3+\sin\bigl(\tfrac{16}{15}(\xi_k-0.7)\bigr)+\sin^2\bigl(\tfrac{16}{15}(\xi_k-0.7)\bigr)$；$f_3(\xi)=\exp\bigl(-\sum_kc_k^2(\xi_k-0.01)^2\bigr)$，$c_k=\exp(-6k/d)$；$f_4(\xi)=\sin(\sum_k\cdot)$。
- §5.2.2 **电阻网络**：$d=2p$ 个不确定电阻，驱动电压 $V_0=1$，关心量是电压 $V$。$d=2$ 时 $\xi_i\sim U[10,100]$；$d=4$ 时 $\xi_1\sim\mathrm{Exp}(0.9)$、$\xi_2\sim\mathrm{Exp}(1.1)$、$\xi_3\sim\mathrm{Exp}(0.8)$、$\xi_4\sim\mathrm{Exp}(1.0)$，**矩由 1000 个样本估计**。
- §5.2.3 **随机椭圆方程**：$-\nabla\cdot(a(y,\omega)\nabla u)=f$ 于 $D=[0,1]^2$，$u|_{\partial D}=0$，确定性载荷 $f=\cos(y_1)\sin(y_2)$，扩散系数取 Babuška–Nobile–Tempone 的对数正态型 $\log(a_N(y,\omega)-0.5)=1+\xi_1(\omega)(\sqrt\pi L)^{1/2}+\sum_{i=2}^5\zeta_ig_i(y)\xi_i(\omega)$；空间方向用标准有限元求解。

作者给出的定性结论是：Christoffel 加权最小二乘在线性与对数线性采样率下都给出稳定而精确的逼近。

> [!warning] 实验的可转录范围与理论的落差
> 本站转录的是配置（分布、维数、采样率、试验次数、矩的样本数），**条件数与误差的具体数值没有转录**。更要紧的是，§5.1 测的是**有限次数**下的条件数，而定理保证只在 $N\to\infty$ 的渐近意义下成立（因为用 $\hat\rho$ 替代了 $\tilde\rho$）；情形 3、4 里的离散与无界输入更是落在定理 3.1 覆盖范围之外或只落在猜想密度上。**实验能说明方法在这些情形下可用，不能说明定理覆盖了这些情形。**

### 与其他论文的关系

编号 36 是[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 22]] 移植到数据驱动 / 认知不确定性设定的版本：同样的 Christoffel 权、同样的平衡测度采样、同样的「只在 $N$ 渐近意义下」的告诫，区别在于正交族由 Hankel 矩数据构造而不是由已知密度给出。它猜想的 Gauss 平衡密度 $C(2-\|\xi\|^2)^{d/2}$ 与[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 45、28]] 中出现的候选密度是同一条。它的直接后继是编号 **44**：用精确的**诱导**采样替换渐近的平衡测度采样，并把最小二乘换成 $\ell_1$，从而去掉 $N\to\infty$ 这一条告诫。它引用的稳定性判据（定理 4.1）与编号 45 的定理 6.1、编号 9 / 14 的定理 1 是同一条 Cohen–Davenport–Leviatan 结果。

## 44：数据驱动多项式混沌的稀疏逼近——诱导采样路线

（本节记号：$Q$ 个经验数据点组成 $S$，从中抽 $M$ 个训练样本，$N=|\Lambda^{TD}_K|$ 个基函数，$M\ll N$。）

### 直觉

两个困难交织在一起：输入分布 $\omega$ 只通过有限样本集 $S$ 已知（认知不确定性），且训练预算可能远小于多项式维数（$M\ll N$），必须做稀疏恢复。编号 36 用数据驱动基加平衡测度加权最小二乘处理了前一个困难，但平衡采样有两个问题——它只在次数极限下最优，而且**在数据驱动情形下它可能把样本放到经验数据几乎没有质量的地方，于是要求在从未查询过模型的位置上求值**。

论文的替换很直接：**用支撑在 $S$ 本身上的精确诱导测度替换平衡测度，用预条件 $\ell_1$ 替换最小二乘**。诱导测度就是把 Christoffel 权反过来用——在基函数大的地方**多**采样，再用 $1/\kappa$ 加权把它抵消掉。这个「采样密度 $\propto\kappa$、权 $\propto1/\kappa$」的倒数配对，正是本页开头那条一致有界性原则的又一次实现，而在离散情形下它还额外买到一件事：**每个样本都落在真实数据点上**。

### 问题设定与推导

**从矩构造任意多项式混沌基。** 朴素的路线是在矩矩阵中解一个关于系数 $\{\beta_j\}_{j=0}^K$ 的线性方程组——这条路线本身说明了一件有意思的事：**真密度 $\rho$ 只需通过它的前 $2K$ 阶矩已知，就足以确定到 $K$ 次为止的全部多项式**。但该矩阵在 $K$ 大时条件数很差。论文改用 Hankel 矩阵

$$
H=\begin{pmatrix}
\nu_0&\nu_1&\cdots&\nu_K\\
\nu_1&\nu_2&\cdots&\nu_{K+1}\\
\vdots&\vdots&\ddots&\vdots\\
\nu_K&\nu_{K+1}&\cdots&\nu_{2K}
\end{pmatrix},
$$

要求 Hamburger 定性 $\det H>0$（**给出的一个充分条件是经验集 $\Xi$ 至少含 $K+1$ 个互异样本**），作 Cholesky 分解 $H=R^\top R$，再由 Mysovskih 定理读出三项递推系数

$$
a_j=\frac{r_{j,j+1}}{r_{j,j}}-\frac{r_{j-1,j}}{r_{j-1,j-1}},
\qquad
b_j=\frac{r_{j+1,j+1}}{r_{j,j}},
\qquad r_{0,0}=1,\ r_{0,1}=0
$$

（与编号 36 完全相同）。论文称这条路线在经验上的条件数远好于直接解线性方程组。**算法 1** 即为这一构造：输入数据与权 $\{\xi_j,w_j\}_{j=1}^Q$ 与最高次数 $K$，输出 $\{a_j,b_j\}$，可选地经递推求值。

**张量化与它坦白的代价。** 设 $S=\{z^{(1)},\dots,z^{(Q)}\}\subset\mathbb R^d$，按坐标边缘化得 $\Xi_i=\{z^{(1)}_i,\dots,z^{(Q)}_i\}$，取均匀权 $w_j=1/Q$，逐坐标跑算法 1，再张量化。所得 $\Phi_\lambda$ 关于**张量化的边缘测度**正交归一，而**不是**关于 $S$：

$$
\int_\Gamma\Phi_\lambda(z)\Phi_\theta(z)\,d\nu(z)=\delta_{\lambda,\theta},
\qquad
d\nu(z)=\bigotimes_{i=1}^d d\nu_i(z_i),
\quad
d\nu_i(z_i):=\frac1Q\sum_{j=1}^Q\delta_{z^{(j)}_i}(z).
$$

论文直言：「按此法构造的多维多项式 $\Phi_\lambda$ 并不尊重 $S$ 的矩」，由此产生的额外代价「本质上等价于把集合 $S$ 张量化」——**也就是说，输入坐标之间的任何相关性都被丢掉了**。这是本方法一处实打实的建模损失，论文没有回避。

**数据上的诱导测度（论文的中心对象）。** 记 $N=|\Lambda^{TD}_K|$，关于张量测度 $\nu$ 的归一化 $\Lambda$-Christoffel 量为

$$
\kappa\bigl(z;\Lambda^{TD}_K\bigr):=\frac1N\sum_{j=1}^N\Phi_j^2(z)
\ \ \Bigl(=\frac{K(z)}{N}\ \text{在本页约定下}\Bigr),
$$

定义支撑在 $S$ 上的**离散诱导测度**

$$
d\mu(z):=\sum_{j=1}^Q\tilde\kappa_j\,\delta_{z^{(j)}},
\qquad
\tilde\kappa_j:=\frac{\kappa(z^{(j)})}{\sum_{q=1}^Q\kappa(z^{(q)})} .
$$

注意方向：$\kappa$ 是 Christoffel 函数的**倒数**，在基函数大的地方大，因此 $\mu$ **上调**那些区域的采样概率。

**预条件 $\ell_1$。**

$$
\min\|c\|_1
\quad\text{s.t.}\quad
\sqrt W A c=\sqrt W b,
\qquad
A_{ij}=\Phi_j(z^{(i)}),\quad b_j=f(z^{(j)}),\quad W_j=1/\kappa\bigl(z^{(j)}\bigr),
$$

其中 $z^{(j)}$ 是从 $\mu$ 独立同分布抽取的 $M$ 个样本。**采样密度 $\propto\kappa$、权 $\propto1/\kappa$——正好互为倒数**，于是复合矩阵 $\sqrt WA$ 关于 $S$ 上的均匀测度无偏，且**其元素一致有界，与 $S$ 和 $\Phi_j$ 的具体形状无关**。这一条一致有界性正是有界正交系假设本身，而压缩感知误差界中的惩罚项与矩阵最大元成正比——**所以 $W$ 的这个取法不是调参，而是被理论指定的**。由于 $S$ 有限，从 $\mu$ 独立采样只需逆变换抽样，代价可以忽略。

**算法 2。** 输入离散分布 / 数据集 $\{z^{(j)},w_j\}_{j=1}^Q$ 与 $f$，输出 $c^*$ 使 $f\approx\hat f=\sum_jc^*_j\Phi_j$：(1) 用算法 1 逐维建立任意多项式空间 $\{\Phi_j\}$；(2) 从诱导测度 $\mu$ 生成 $M$ 个独立同分布样本；(3) 组装 $b_i=f(z_i)$ 与 $A_{ij}=\Phi_j(z^{(i)})$；(4) 由 $W_j=1/\kappa(z^{(j)})$ 组装 $W$；(5) 解预条件 $\ell_1$ 问题。

**与平衡测度的关系（论文明确写出的三点）。** 对相当一般的 $\omega$ 存在唯一的 $\mu_\infty$ 使 $d\mu_\infty=\lim_{K\to\infty}\omega\kappa$——即**加权多重位势平衡测度**——但它的显式形式只在部分情形已知；在超立方体 $[-1,1]^d$ 上它是乘积 Chebyshev 测度 $\frac{1}{\pi^d\prod_i\sqrt{1-z_i^2}}$。论文主张 $\mu$ 优于 $\mu_\infty$ 的理由是：

1. $\mu$ **只支撑在 $S$ 上**，因此每个样本都是真实存在数据的位置；
2. 乘积 Chebyshev 测度「也许在极限意义下最优，但对有限次数 $K$ 并不最优」；
3. 结构上，平衡采样只轻微依赖 $\omega$ 与 $K$，而诱导测度对二者都显式依赖——**这正是论文预期（并观察到）它在数据驱动情形下表现更好的原因**。

### 定理

**本文不含任何新定理。** 这是一篇算法与实验的论文（全文没有编号的 Theorem / Lemma / Proposition 环境）。它的论证方式是引用：一是有界正交系的压缩感知理论（误差项与 $\sqrt WA$ 的最大元成正比，因而决定了 $W$ 的取法），二是[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 45]] 综述中的诱导采样最小二乘理论，再加上上面列出的两条结构性优势。

> [!warning] 不要给这个算法安一个样本复杂度指数
> 论文**没有**陈述本算法的样本复杂度。诸如「$M\gtrsim s\log^3s\log N$」之类的式子在编号 21、24、29 中各有其明确的假设，但它们**都不是**本文对算法 2 的结论。任何此类指数在本文中都无出处。

### 数值实验

三个方法标签贯穿全部实验：**induced distribution**（本文方法）、**CSA**（平衡测度采样的 $\ell_1$，即[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 24]] 的方法）、**MC**（从 $\omega$ 独立同分布采样、不加预条件的 $\ell_1$）。全部结果对 100 次独立试验取平均；精度误差在 $E=10{,}000$ 个从 $\omega$ 独立抽取的样本上以离散范数计算（**原文印出的公式除以 $M$ 而非 $E$，标为讹误**）。

**图 1（示意）** 是本专题里关于「非渐近设计与渐近设计有何不同」最干净的一张图。数据集取张量形式 $S=\{(z^{(j)}_1,z^{(k)}_2)\}_{1\le j,k\le24}$，权 $w_{j,k}=u_jv_k$：$z^{(j)}_1$ 是 $[-1,1]$ 上 24 个等距点、$u_j$ 取 $\mathrm{Binomial}(24,0.5)$ 的概率质量函数；$z^{(k)}_2$ 也是 24 个等距点、$v_k$ 取 $\mathrm{Poisson}(10)$ 概率质量函数的前 24 个值（截断掉的质量约 $10^{-5}$，再归一化）。结果是平衡（乘积 Chebyshev）样本与 $K=20$ 的诱导样本**明显聚集在不同区域**，**证实 $K=20$ 还远未进入渐近区**。

**§4.1 精确稀疏恢复。** $\omega$ 取张量形式、各向同性的边缘分布是三者的等权混合：$[-1,1]$ 上的均匀、$[-1,1]$ 上截断的 $N(0.2,1.5)$、$[0,1]$ 上截断的对数正态。$S$ 取 $Q=10^5$ 个从 $\omega$ 独立同分布的样本；$c^*$ 随机 $s$-稀疏、非零项标准正态；判定成功的准则是 $\|c-c^*\|_\infty<10^{-3}$。

| 场景     | $d$ | $K$ | $N$ | $s$ | 结论                    |
| -------- | --- | --- | --- | --- | ----------------------- |
| 低维高次 | 2   | 20  | 231 | 8   | 诱导采样**显著优于** MC |
| 高维低次 | 10  | 3   | 286 | 8   | 两者**接近**            |

高维低次一栏的解释是可推广的：**低次空间使 $\mu$ 接近 $\omega$，诱导采样与朴素采样自然趋同**——编号 45 在它的例 8.1 中报告了同样的次数依赖性。

**§4.2 解析函数逼近。** 测试函数为 $f_1(z)=\exp\bigl(-\sum_iz_i\bigr)$；$f_2(z)=\sum_i(1-z_{i-1})^2+\sum_{i\ge2}100(z_i-z_{i-1}^2)^2$（Rosenbrock）；$f_3(z)=\sin(\sum_iz_i)$；$f_4(z)=\bigl(1+\frac{1}{2d}\sum_ic_i(1+z_i)\bigr)^{-d-1}$，$c_i=\frac{1+i}{4d}$（原文把两个函数都编号为 $f_3$，**标为编号讹误**）。二维情形用图 1 的 $S$ 与 $K=20$：**$M$ 较大时诱导采样明显优于 MC 与 CSA，而 $M$ 相对较小时 MC 反而略胜**。改用 §4.1 的混合数据，以及 $d=5$、$K=7$ 的情形重复：$M$ 增大时诱导采样同样胜过 MC。

**§4.3 参数化偏微分方程。** 有界多边形区域 $D\subset\mathbb R^2$ 上的固支 **Kirchhoff 板弯曲**问题：$-M_{IJ,IJ}(u)=f$ 于 $\Gamma\times D$，$u=\partial_nu=0$ 于 $\Gamma\times\partial D$，$M_{IJ}(u)=D(z,x)\bigl((1-\nu)K_{IJ}(u)+\nu K_{LL}(u)\delta_{IJ}\bigr)$，$K_{IJ}(u)=-\partial_{IJ}u$，抗弯刚度 $D(z,x)=\frac{E(z,x)h^3}{12(1-\nu^2)}$。**随机 Young 模量**经截断 Karhunen–Loève 展开进入：$Y(z,x)=\log(E(z,x)-100)$，$Y=1+Z_1(\sqrt\pi L)^{1/2}+\sum_{i=2}^d\zeta_ig_i(x)Z_i$，协方差只沿 $x_1$ 方向取平方指数型 $K(x_1,x_1')=\exp\bigl(-(x_1-x_1')^2/L_c^2\bigr)$，特征函数为 $g_i(x)=\sin\bigl(\lfloor i/2\rfloor\pi x_1/L_p\bigr)$（$i$ 偶）与 $\cos\bigl(\lfloor i/2\rfloor\pi x_1/L_p\bigr)$（$i$ 奇）（**符号在本站所用提取中有歧义，标此存疑**）。

> [!warning] 实验说明了什么，没说明什么
> 本站没有转录任何恢复率或误差数值。更重要的是：**本文没有定理**，因此这里不存在「实验与理论的差距」问题，只存在「实验支持的论断范围」问题。实验支持的是：在数据驱动、有限次数、$M$ 不太小的情形下，诱导采样通常优于平衡测度采样与朴素蒙特卡罗。它**不**支持：诱导采样总是更好（§4.2 中 $M$ 小时 MC 略胜）；也**不**支持任何关于所需样本数的定量结论。

### 与其他论文的关系

编号 44 是编号 **36** 的直接后继——同样的 aPC / Hankel 基构造、同一批作者再加上 Narayan——只做了两处替换：用精确的诱导采样替换渐近的平衡测度采样，用 $\ell_1$ 替换最小二乘。它实验中的对照方法 "CSA" **就是**[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|编号 24]]，论文自己把两者的关系写得很清楚：编号 24「使用的正是同一个预条件矩阵 $W$……但选择不从 $\mu$ 采样，而是从 $\mu$ 的 $K$-渐近版本、即多重位势平衡测度采样」；该渐近选择在最小二乘框架下的对应物是编号 22。它引用编号 45 作为诱导测度采样的综述。$W$ 的有界正交系论证与编号 **21**、24 用的是同一套限制等距机制。

> [!note] 覆盖进度
> 编号 10、21、29、32、36、44 六篇均已通读全文（10、21、32、36、44 为 arXiv 全文，29 为 CiCP 全文），因此本页对六篇都给出了直觉、完整推导、带假设的定理陈述，以及数值实验的具体配置。仍有三类限制需要说明：
>
> 1. **实验数值。** 六篇的成功率曲线与误差曲线都以图的形式发表。本页转录了实验配置（维数、次数、稀疏度、样本数、试验次数、判定阈值、测试函数与算例）与作者陈述的定性结论，**没有逐点的成功率或误差数值**。
> 2. **排版歧义。** 编号 32 的 Legendre 情形两个常数（式 (14) 的首项系数与 $W^0_{n,n}$ 的指数结构）、编号 29 的两处归一化常数（问题 2 的 $W$ 与最优 $\alpha$）、编号 36 的定理 4.1 中的 $\delta$ 与 (4.5)/(4.7) 里 $K$ 的用法，在本站所用的全文提取中排版有歧义或自相矛盾。这些位置都已就地以警告标出并给出自洽形式，**引用前请核对排版原文**。
> 3. **各篇自身的限制。** 编号 21 主定理印出的概率不等式方向是反的（已修正并标明）；编号 32 明说没有证出 $\beta_L\le\mu_L$，且**不含任何样本复杂度定理**；编号 36 的两条定理都是引用他人结果，其无界域平衡测度是**猜想**，且未分析矩估计误差；编号 44 **不含任何新定理**。

## 六篇的关系

| 编号 | 测量的构成             | 采样设计              | 基的来源     |
| ---- | ---------------------- | --------------------- | ------------ |
| 10   | 函数值，确定性点       | Weil 和点集           | 已知分布     |
| 21   | 函数值，Gauss 网格子集 | 隐含 Christoffel 加权 | 已知分布     |
| 29   | 值 + 梯度，三角多项式  | 需预条件              | 已知分布     |
| 32   | 值 + 梯度，多项式混沌  | 需预条件              | 已知分布     |
| 36   | 函数值，加权最小二乘   | 由经验测度决定        | **经验测度** |
| 44   | 函数值，稀疏逼近       | 诱导采样              | **经验测度** |

按理论路线与所付代价再看一遍：

| 编号 | 理论路线                                   | 主要保证                                                 | 该保证的代价                                     |
| ---- | ------------------------------------------ | -------------------------------------------------------- | ------------------------------------------------ |
| 10   | 互相干性 + Weil 和                         | 确定性精确恢复                                           | $m\gtrsim4^{\,q}d^2s^2$，稀疏度上二次            |
| 21   | 限制等距 / 有界正交系                      | $M\gtrsim L(\mathbf n)s\log^3s\log N$                    | 无界参数下 $L(\mathbf n)\sim n^{2d/3}$           |
| 29   | 互相干性（定理 2.1）+ 限制等距（定理 3.1） | $\mu(\tilde\Phi)\le\lambda\mu(\Phi)$，$\lambda<1$ 确定性 | 需可容许方向；$R_0$ 可能很大                     |
| 32   | 互相干性                                   | 均值各向同性 + 与次数无关的相干性界                      | **没有样本复杂度定理**；未证 $\beta_L\le\mu_L$   |
| 36   | 引用的最小二乘稳定性                       | $M\gtrsim N\log N$                                       | **只在 $N\to\infty$ 渐近成立**；无界域密度是猜想 |
| 44   | 引用的有界正交系理论                       | $\sqrt WA$ 元素一致有界                                  | **无新定理**；张量化丢掉坐标间相关性             |

一条贯穿的判断：**稀疏恢复的样本复杂度由行范数的均匀性决定，因此凡是改变行结构的操作（加梯度行、换基、换采样密度）都必须配一个相应的预条件。** 这条判断把本页与[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|最优采样一页]]统一起来：两者处理的是同一个量，只是一个在最小二乘框架下、一个在 ℓ1 框架下。

第二条判断关于确定性与概率性的取舍：编号 10 与编号 29 的定理 2.1 是本专题里仅有的两处**确定性**结论，它们的代价分别是稀疏度上的二次样本数与对方向可容许性的额外假设。**在这个领域里，去掉「以高概率」这个前提是要付费的。**

## 本页原文

- Z. Xu and T. Zhou, [_On sparse interpolation and the design of deterministic interpolation points_](https://doi.org/10.1137/13094596X), SIAM J. Sci. Comput. 36(4) (2014), pp. A1752-A1769。
- L. Guo, A. Narayan, T. Zhou, and Y. Chen, [_Stochastic collocation methods via ℓ1 minimization using randomized quadratures_](https://doi.org/10.1137/16M1059680), SIAM J. Sci. Comput. 39(1) (2017), pp. A333-A359。
- Z. Xu and T. Zhou, [_A gradient-enhanced ℓ1 approach for the recovery of sparse trigonometric polynomials_](https://doi.org/10.4208/cicp.OA-2018-0006), Commun. Comput. Phys. 24 (2018), pp. 286-308。
- L. Guo, A. Narayan, and T. Zhou, [_A gradient enhanced ℓ1-minimization for sparse approximation of polynomial chaos expansions_](https://doi.org/10.1016/j.jcp.2018.04.026), J. Comput. Phys. 367 (2018), pp. 49-64。
- L. Guo, Y. Liu, and T. Zhou, [_Data-driven polynomial chaos expansions: a weighted least-square approximation_](https://doi.org/10.1016/j.jcp.2018.12.020), J. Comput. Phys. 381 (2019), pp. 129-145。（另有文献把此文的页码引作同卷 pp. 110-128，两者不一致，此处采用 129-145。）
- L. Guo, A. Narayan, Y. Liu, and T. Zhou, [_Sparse approximation of data-driven polynomial chaos expansions: an induced sampling approach_](https://doi.org/10.4208/cmr.2020-0010), Commun. Math. Res. 36 (2020), pp. 128-153。
