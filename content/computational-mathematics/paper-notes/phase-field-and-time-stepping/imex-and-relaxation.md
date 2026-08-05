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

编号 78 处理的是另一类思路：**不改进时间离散，而改写能量本身。** 标量辅助变量与不变能量二次化方法把非线性能量重写为一个带辅助变量的二次形式，从而每步只需解线性系统。本文在这一族中给出带正则化能量重构的线性松弛方法，用于相场模型。

> [!note] 覆盖进度
> 本页尚未对编号 78 做逐式核对，因此不报告其松弛形式、正则化项与稳定性结论。它在本专题中的定位是：与编号 91、104 一样以「显式处理非线性、保住能量论证」为目标，但手段是重构能量而不是设计格式或寻找判据。

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
