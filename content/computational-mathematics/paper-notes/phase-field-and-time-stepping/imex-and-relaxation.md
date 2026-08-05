---
title: IMEX 与松弛型格式
description: 编号 78、91、104：把非线性项显式处理，同时保住能量论证
lang: zh
translation: en/computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation
tags:
  - 论文笔记
  - 相场模型
  - 隐显格式
---

> [!note] 本页覆盖
> 编号 **78**（_J. Comput. Phys._ 515, 2024）、**91**（_Math. Comput._ 95(359), 2026）、**104**（投稿 _SIAM J. Numer. Anal._，[arXiv:2605.05619](https://arxiv.org/abs/2605.05619)）。

隐显格式把线性刚性部分隐式处理、非线性部分显式处理，因此每级不需要内层迭代。代价是能量论证变难：显式处理的非线性项不再自动提供符号确定的贡献。这三篇给出三种不同的对策。

## 91：让平均耗散率与时空参数无关

### 两个缺陷

Cahn-Hilliard 模型是 Ginzburg-Landau 自由能的 $H^{-1}$ 梯度流：

$$
E[\Phi]=\int_\Omega\Bigl[\tfrac{\epsilon^2}{2}|\nabla\Phi|^2+F(\Phi)\Bigr]\mathrm d\mathbf x,
\qquad F(\Phi)=\tfrac14(\Phi^2-1)^2,
$$

$$
\partial_t\Phi=\Delta\bigl[F'(\Phi)-\epsilon^2\Delta\Phi\bigr],
\qquad
\frac{\mathrm dE}{\mathrm dt}
=-\bigl((-\Delta)^{-1}\partial_t\Phi,\partial_t\Phi\bigr)_{L^2}\le0 .
$$

要证明隐显 Runge-Kutta 方法继承**原始**能量耗散律，必须先知道**各级解在最大范数下一致有界**，稳定化分裂的能量论证才能闭合。而此前的统一框架直接**假设**非线性主体 $F'$ 全局 Lipschitz，Cahn-Hilliard 的四次势 $F(\Phi)=\frac14(\Phi^2-1)^2$ 并不满足。

第二个缺陷是实践性的：一般隐显 Runge-Kutta 方法的平均耗散率依赖 $\tau_n\overline\lambda_{\mathrm{ML}}$（时间步乘离散算子的平均特征值），因此自适应算法一旦取大步长，方法的有效耗散就改变，计算出的能量曲线随之漂移。

### 差分形式与「精化」条件

以稳定化 $L_\kappa\Phi=-\epsilon^2\Delta\Phi+\kappa\Phi$、$f_\kappa(\Phi)=\kappa\Phi-F'(\Phi)$ 改写方程后，$s$ 级隐显 Runge-Kutta 格式为

$$
u_h^{n,i}=u_h^{n,1}
+\tau_n\sum_{j=1}^{i}a_{ij}\Delta_hL_{\kappa,h}u_h^{n,j}
-\tau_n\sum_{j=1}^{i-1}\hat a_{ij}\Delta_h f_\kappa(u_h^{n,j}),
$$

其中隐式部分 $A$ 是刚性精确的对角隐式 Runge-Kutta（第一级显式，首尾相同），$\widehat A$ 严格下三角，并要求**节点条件** $\hat{\mathbf c}=\mathbf c$（等价于 $A\mathbf 1=\widehat A\mathbf 1$），这使方法在所有级上相容并保持平衡态。

论文把格式改写成差分形式，引入一个矩阵

$$
D(z)=D_{\mathrm E}-zD_{\mathrm{EI}},
\qquad
D_{\mathrm E}=A_{\mathrm E}^{-1}E_{s_{\mathrm I}},
\qquad
D_{\mathrm{EI}}=A_{\mathrm E}^{-1}A_{\mathrm I}E_{s_{\mathrm I}}-E_{s_{\mathrm I}}+\tfrac12 I_{s_{\mathrm I}},
$$

其中 $E_{s_{\mathrm I}}$ 是下三角全一矩阵。$D$ 正（半）定指其对称部分正（半）定。**这个 $D$ 承担的正是[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|BDF 一族]]中离散正交卷积核的角色**，只是索引从时间层换成级：相应的正交性恒等式为

$$
\sum_{i=j}^{k}d^{(R)}_{k,i}\,\underline{\hat a}_{i+1,j}\equiv\delta_{kj},
\qquad
(\underline{\hat a}_{i+1,j})=E_{s_{\mathrm I}}^{-1}A_{\mathrm E} .
$$

在 $D_{\mathrm E}$ 与 $D_{\mathrm{EI}}$ 半正定的前提下，各级能量律带有平均耗散率

$$
\mathcal R=\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}\frac1{\hat a_{k+1,k}}
+\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}
\Bigl(\frac{a_{k+1,k+1}}{\hat a_{k+1,k}}-\frac12\Bigr)\tau_n\overline\lambda_{\mathrm{ML}}
\ \ge 0 .
$$

一个方法「好」的标准是 $\mathcal R$ 在很大的 $\tau_n\overline\lambda_{\mathrm{ML}}$ 范围内都接近 $1$。

**论文的中心观察是：$\mathcal R$ 与 $\tau_n\overline\lambda_{\mathrm{ML}}$ 无关当且仅当 $D_{\mathrm{EI}}=\mathbf 0$**，即

$$
A_{\mathrm I}=A_{\mathrm E}P_{s_{\mathrm I}},
\qquad
P_{s_{\mathrm I}}=I_{s_{\mathrm I}}-\tfrac12E_{s_{\mathrm I}}^{-1} .
$$

此时 $D_{\mathrm R}=D_{\mathrm E}=A_{\mathrm E}^{-1}E_{s_{\mathrm I}}$ 与 $z$（因而与两个网格参数）无关，且

$$
\mathcal R_{\mathrm R}=\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}\frac1{\hat a_{k+1,k}} .
$$

格式塌缩为紧凑的精化形式

$$
u_h^{n,i+1}=u_h^{n,1}+\tau_n\sum_{j=1}^{i}\hat a_{i+1,j}\Delta_h
\Bigl[L_{\kappa,h}u_h^{n,j+\frac12}-f_\kappa(u_h^{n,j})\Bigr].
$$

### 一条结构性的排除

节点条件迫使隐式部分的第一列 $\mathbf a_1=(\tfrac12\hat a_{21},\dots,\tfrac12\hat a_{s1})^T\ne\mathbf 0$，因此这类方法必然是 **Lobatto 型，绝不可能是 Radau 型或 ARS 型**。论文把这一点写成命题：不存在平均耗散率与 $\tau_n\overline\lambda_{\mathrm{ML}}$ 无关的 Radau 型或 ARS 型隐显 Runge-Kutta 方法。

这类结论在方法设计中很有价值：它不是「我们选了 Lobatto 型」，而是「要这条性质就只能是 Lobatto 型」。

两个具体方法为：两级情形取 $\theta=1/2$，得到 $\mathcal R_{\mathrm R}=1$ 且格式为 Crank-Nicolson 型

$$
\delta_\tau\phi_h^{n}=\tau_n\Delta_h\Bigl[\tfrac12L_{\kappa,h}(\phi_h^{n}+\phi_h^{n-1})
-f_\kappa(\phi_h^{n-1})\Bigr];
$$

二阶情形不存在三级方法，四级单参数族由 $\hat a_{32}=\hat a_{43}=c_2$、$\hat a_{42}=\frac1{2c_2}-c_3$ 给出，其中 $c_3$ 是

$$
c_3^2-\Bigl(\frac1{2c_2}+c_2\Bigr)c_3+(c_2-1)^2=0
$$

的根，参数范围 $0<c_2\le\frac{2+\sqrt6}{2}$，正定性由两个显式行列式条件判定。

## 104：用半生成函数取代逐族寻找乘子

### 问题

对非线性抛物方程，隐显多步方法把线性刚性算子隐式、非线性显式处理。已有的族很多——加权 BDF、修正 BDF、推广 BDF、NIMEX——都以扩大绝对稳定域为设计目标。但三阶及更高阶的变体**不 $A$-稳定**，因此严格的**离散能量**稳定性此前是开放的。

两条现有路线都要求逐族寻找特设对象。第一条是基于 Dahlquist $G$-稳定的 Nevanlinna-Odeh 型乘子：加权 BDF、修正 BDF、推广 BDF、NIMEX 各需要一族不同的乘子。第二条是 Huang 与 Shen 型的隐式部分分解

$$
\sum_{j=0}^{k-1}b_{\mathrm G,j}^{(k,\beta)}v^{n-j}
=\eta_k(\beta)\sum_{j=0}^{k-1}c_{\mathrm G,j}^{(k,\beta)}v^{n-j}
+\sum_{j=0}^{k-1}d_j^{(k,\beta)}v^{n-j},
$$

其中 $\eta_2=\frac{\beta-1}{\beta}$、$\eta_3=\frac{\beta-1}{\beta+1}$、$\eta_4=\frac{\beta-1}{\beta+3}$，稳定性条件为 $\eta_k(\beta)>\mu_0/\varpi$。这条路线的局限很具体：它在 $\beta=1$ 处退化（因此**从不覆盖普通 BDF-$k$**），只对 $2\le k\le4$ 存在，而其三项精化版本只在固定的 $\beta_k=3,6,9$ 处有效；对推广 BDF5 以及加权 BDF、修正 BDF、NIMEX 都没有对应物。

**本文把上述全部替换为一次计算：三个显式有理函数在单位圆上的三个极值。**

### 抽象设定与全局离散能量方法

在 Hilbert 三元组 $V\subset H=H'\subset V'$ 上考虑 $u_t+\varpi\mathcal Lu=\mathcal F(u)$，$\mathcal L$ 正定自伴有界线性，$\mathcal F$ 可非线性，并设**局部** Lipschitz 条件

$$
\|\mathcal F(v)-\mathcal F(w)\|_{\star}\le\mu_0\|v-w\|_V+\mu_1\|v-w\|_H,
\qquad \mu_0\in(0,\varpi).
$$

该框架也覆盖非自伴 $\mathcal L=\mathcal L_s+\mathcal L_a$，办法是把低阶的反自伴部分移入显式项。

均匀网格上的一般 $k$ 步隐显多步方法写成

$$
\sum_{j=0}^{k-1}a_j^{(k)}\partial_\tau u^{n-j}
+\varpi\sum_{j=0}^{k}b_j^{(k)}\mathcal Lu^{n-j}
=\sum_{j=0}^{k-1}c_j^{(k)}\mathcal F(u^{n-j-1})+\mathfrak C^{(k)}_n(u^0),
$$

因此一个方法就是一个三元组 $(\vec a^{(k)},\vec b^{(k)},\vec c^{(k)})$。用 $\vec a^{(k)}$ 的离散正交卷积核

$$
a_0^{(-1,k)}=\frac{1}{a_0^{(k)}},
\qquad
a_j^{(-1,k)}=-\frac{1}{a_0^{(k)}}\sum_{i=1}^{j}a^{(-1,k)}_{j-i}a_i^{(k)},
$$

满足双向正交性 $\sum_{\ell=j}^{n}a^{(-1,k)}_{n-\ell}a^{(k)}_{\ell-j}\equiv\delta_{nj}$，作用后格式化为「差分」形式

$$
\partial_\tau u^{n}+\varpi\sum_{\ell=1}^{n}\hat b^{(k)}_{n-\ell}\mathcal Lu^{\ell}
=\sum_{\ell=1}^{n}\hat c^{(k)}_{n-\ell}\mathcal F(u^{\ell-1})
+\sum_{\ell=1}^{n}a^{(-1,k)}_{n-\ell}\mathfrak C^{(k)}_\ell(u^0),
$$

复合核为 $\hat b_j^{(k)}=\sum_{i=0}^{j}a^{(-1,k)}_{j-i}b_i^{(k)}$、$\hat c_j^{(k)}=\sum_{i=0}^{j}a^{(-1,k)}_{j-i}c_i^{(k)}$。矩阵形式下 $\widehat B_{L,k}=A_{L,k}^{-1}B_{L,k}$、$\widehat C_{L,k}=A_{L,k}^{-1}C_{L,k}$，均为下三角 Toeplitz 且可交换。由于没有丢弃信息，作者称之为**全局离散能量方法**。

这条思路与[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 58]] 完全同源：用离散正交卷积核把多步格式逆回单个一阶差分，从而恢复教科书式的能量论证。差别在于这里同时有三组核（隐式、显式、差分），因此需要的判据也是复合的——这正是「半生成函数」要处理的对象。

## 78：线性松弛与正则化能量重构

编号 78 处理的是另一类思路：**不改进时间离散，而改写能量本身。** 论文把自己的方法命名为 RRER（relaxation with regularized energy reformulation）。

### 问题：IEQ 与 SAV 都要对辅助变量求时间导数

不变能量二次化（IEQ）与标量辅助变量（SAV）都把非线性能量写成带辅助变量的二次形式，从而每步只解线性系统。IEQ 通常给出**耦合**且系数**依赖时间**的系统，SAV 则给出**解耦**且系数为**常数**的系统。但两者有一处共同的做法：辅助变量的演化方程是靠**对该变量求时间导数**得到的，这本身引入截断误差，也使离散系统对原方程的忠实度下降。论文的原话是它「不需要对辅助变量求时间导数」。

本文的做法是让辅助变量只由一个**代数**关系定义，从不参与时间求导。

### 正则化辅助变量

设自由能与梯度流为

$$
E(\phi)=\tfrac12(\mathcal L\phi,\phi)+\bigl(F(\phi),1\bigr),
\qquad
\frac{\partial\phi}{\partial t}=-\mathcal G\bigl(\mathcal L\phi+F'(\phi)\bigr),
$$

$\mathcal L$ 线性、$\mathcal G\ge0$ 为迁移算子，于是 $\frac{\mathrm d}{\mathrm dt}E=-\bigl(\mathcal L\phi+F'(\phi),\mathcal G(\mathcal L\phi+F'(\phi))\bigr)\le0$。取 $C_0$ 使 $F(\phi)+C_0\ge0$，并引入**稳定化参数** $\gamma$，定义

$$
q=\sqrt{4\bigl(F(\phi)+C_0\bigr)}-\gamma
\qquad\Longrightarrow\qquad
F(\phi)=\tfrac14q^2+\tfrac12\gamma q+\tfrac{\gamma^2}{4}-C_0 .
$$

等价系统与相应的能量为

$$
\frac{\partial\phi}{\partial t}=-\mathcal G\mu,
\quad
\mu=\mathcal L\phi+\tfrac12qq'+\tfrac{\gamma}{2}q',
\quad
q=\sqrt{4\bigl(F(\phi)+C_0\bigr)}-\gamma,
$$

$$
\widehat E(\phi,q)=\tfrac12(\mathcal L\phi,\phi)
+\Bigl(\tfrac14q^2+\tfrac12\gamma q+\tfrac{\gamma^2}{4}-C_0,\,1\Bigr),
\qquad
\frac{\mathrm d}{\mathrm dt}\widehat E=-(\mu,\mathcal G\mu)\le0 .
$$

**有一处必须强调：$\widehat E(\phi,q)$ 在连续层面与 $E(\phi)$ 精确相等，而不只是近似。** 论文把这个恒等式逐行推了出来。这与 SAV、IEQ 的修正能量不同，后者与原能量之间只有近似关系。这一点在下面的定位讨论中是关键。

### 交错时间网格使格式线性

以带斜率选择的分子束外延模型为例（**注意与编号 52 的模型不同**，见下），

$$
E(\phi)=\int_\Omega\Bigl(\frac{\epsilon^2}{2}(\Delta\phi)^2+\frac14\bigl(|\nabla\phi|^2-1\bigr)^2\Bigr)\mathrm d\mathbf x,
\qquad
\phi_t=-\epsilon^2\Delta^2\phi+\nabla\cdot\bigl((|\nabla\phi|^2-1)\nabla\phi\bigr),
$$

取 $q=|\nabla\phi|^2-1-\gamma$（$\gamma>0$）。Crank-Nicolson 型格式放在**交错**网格上：$\phi$ 在整数层、$q$ 在半整数层，

$$
\begin{aligned}
&\text{(a)}\ \ \frac{\phi^{n+1}-\phi^{n}}{\delta t}
=-\epsilon^2\Delta g^{n+\frac12}
+\nabla\cdot\Bigl(q^{n+\frac12}\nabla\frac{\phi^{n+1}+\phi^{n}}{2}\Bigr)
+\gamma g^{n+\frac12},\\
&\text{(b)}\ \ g^{n+\frac12}=\Delta\frac{\phi^{n+1}+\phi^{n}}{2},
\qquad
\text{(c)}\ \ \frac{q^{n+\frac12}+q^{n-\frac12}}{2}=|\nabla\phi^{n}|^2-1-\gamma .
\end{aligned}
$$

**格式线性的原因就在 (c)：它是代数的，且右端只含已知的 $\phi^n$。** 因此每步只解一个线性代数系统。时间二阶精度，并且没有任何步长比限制——步长均匀，稳定性是无条件的。

论文对相场晶体模型作同样处理（取 $q=\phi^2-b_0-\gamma$），并把 RRER 推广到三元相场模型与晶粒生长模型。

### 两条定理

对分子束外延模型，格式**守质量**（$\int_\Omega\phi^{n+1}=\int_\Omega\phi^n$）且**无条件能量稳定**。数值部分用 $P_1$ 有限元在 FreeFEM 中实现：制造解检验给出的时间收敛率为 $2.00,2.00,2.00$（IEQ 为 $1.94,1.97,1.99$），相场晶体情形 RRER 比 IEQ 与指数型 SAV 都更准，且 CPU 时间最省；粗化、条纹与三角图样、球面与环面上的算例都与文献相图一致。

### 这一篇为什么是本专题的例外

值得把这一点写清楚，因为它关系到整个专题的内在张力。编号 78 没有 Liao 的合作、没有 DOC/DCC 核、没有变步长、也没有步长比分析——它属于 IEQ/SAV/松弛这一支，而不是卷积核那一支。

更要紧的是**能量陈述的类型不同**：编号 78 证明的是**修正**能量 $\widehat E(\phi,q)$ 的耗散，而编号 40、48、52、57、91 恰恰在设法**避免**这一点，它们要的是原能量或变分能量的耗散。这个张力在文献里是明写着的：编号 91 的引言正是批评基于 SAV 的高阶格式只能建立「关于含辅助变量的修正能量」的稳定性。编号 78 的缓解之处在于 $\widehat E\equiv E$ 在连续层面精确成立，因此它的修正能量不是一个新对象，而是同一个能量的另一种写法。

## 三篇的对策对照

| 编号 | 被显式处理的对象 | 保住能量论证的手段                 |
| ---- | ---------------- | ---------------------------------- |
| 78   | 非线性势         | 重构能量（正则化的二次化）         |
| 91   | 非线性主体       | 设计格式使平均耗散率与网格参数无关 |
| 104  | 非线性项         | 给出统一判据（三个有理函数的极值） |

三条路线的共同点是它们都不试图直接证明「显式项贡献符号确定」，而是各自换一个可控对象：编号 78 换能量，编号 91 换 Butcher 表的结构条件，编号 104 换判据的形式。**这与本专题的整体取向一致：先找到承载全部困难的那个代数对象，再为它建立独立的判据。**

## 覆盖核对

| 内容                                 | 论文 | 覆盖状态                                                                  |
| ------------------------------------ | ---- | ------------------------------------------------------------------------- |
| Cahn-Hilliard 能量与 $H^{-1}$ 梯度流 | 91   | 能量、方程、耗散律                                                        |
| 两个缺陷                             | 91   | 全局 Lipschitz 假设失效；耗散率随步长漂移                                 |
| 稳定化、节点条件与差分形式           | 91   | $L_\kappa$、$\hat{\mathbf c}=\mathbf c$、$D(z)$                           |
| $D$ 与 DOC 核的对应                  | 91   | 正交性恒等式与索引的替换                                                  |
| 平均耗散率与精化条件                 | 91   | $\mathcal R$、$D_{\mathrm{EI}}=\mathbf 0$、$A_{\mathrm I}=A_{\mathrm E}P$ |
| Lobatto 型的必然性                   | 91   | 第一列非零与排除性命题                                                    |
| 两个具体方法                         | 91   | 两级 Crank-Nicolson 型；四级单参数族与其范围                              |
| 现有两条路线的具体局限               | 104  | 乘子逐族寻找；分解在 $\beta=1$ 退化、$k$ 受限                             |
| 抽象设定与局部 Lipschitz 条件        | 104  | Hilbert 三元组、$\mu_0\in(0,\varpi)$、非自伴处理                          |
| 三组核与全局离散能量方法             | 104  | 隐显多步的三元组、DOC 核、复合核、矩阵形式                                |

## 本页原文

- J. Zhang, X. Guo, M. Jiang, T. Zhou, and J. Zhao, [_Linear relaxation method with regularized energy reformulation for phase field models_](https://doi.org/10.1016/j.jcp.2024.113225), J. Comput. Phys. 515 (2024), 113225。
- H.-l. Liao, T. Tang, X. Wang, and T. Zhou, [_A class of refined implicit-explicit Runge-Kutta methods with robust time adaptability and unconditional convergence for the Cahn-Hilliard model_](https://doi.org/10.1090/mcom/4090), Math. Comput. 95(359) (2026), pp. 1293-1325（预印本 [arXiv:2412.07321](https://arxiv.org/abs/2412.07321)）。
- H.-l. Liao, C. Quan, T. Tang, and T. Zhou, _A semi-generating function approach to the stability of implicit-explicit multistep methods for nonlinear parabolic equations_, [arXiv:2605.05619](https://arxiv.org/abs/2605.05619)，投稿 SIAM J. Numer. Anal.
