---
title: 离散最小二乘逼近
description: 编号 6、9、11、13、14：从「取哪些点」到「取多少个点才稳定」
lang: zh
translation: en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares
tags:
  - 论文笔记
  - 不确定性量化
  - 最小二乘
---

> [!note] 本页覆盖
> 编号 **6**（_Commun. Comput. Phys._ 16, 2014）、**9**（_SIAM J. Sci. Comput._ 36(5), 2014）、**11**（_SIAM J. Sci. Comput._ 36(5), 2014）、**13**（_J. Comput. Phys._ 298, 2015）、**14**（_Commun. Comput. Phys._ 18, 2015）。

![配点设计的统一流程](assets/diagrams/tao-zhou-papers/zh/sampling-design.svg)

## 共用装置：行范数、Christoffel 函数与 Gauss 权

### 直觉：稳定性只看设计矩阵各行的范数是否均匀

这五篇论文表面上在做五件事——比较候选点、构造确定性网格、换基与换密度、子采样 Gauss 网格、在非结构网格上配点——但它们受同一个量支配，先把这个量说清楚，后面每一篇的取舍都会变得显然。

设 $\{\varphi_n\}_{n=1}^{N}$ 是 $L^2_w$ 上的一组标准正交基，取一个偏置函数 $q>0$ 满足 $\|q\|_{L^2_w}=1$，从密度 $\rho:=q^2w$ 独立抽 $x_1,\dots,x_M$，加权设计矩阵与右端为

$$
(A)_{m,n}=\frac{1}{\sqrt{M\,q^2(x_m)}}\,\varphi_n(x_m),
\qquad
(f)_m=\frac{1}{\sqrt{M\,q^2(x_m)}}\,f(x_m),
$$

Gram 矩阵 $G=A^{\mathsf T}A$ 满足 $\mathbb E\,G=I$，且对任意可容许的 $q$ 都有 $G\to I$ 几乎必然成立。取 $q\equiv1$ 就退化为「从正交测度直接采样、不加权」的朴素做法。

关键在于 $A$ 的**第 $m$ 行的范数平方**恰好是

$$
\|A_{m,\cdot}\|_2^2=\frac1M\sum_{n=1}^{N}\Bigl(\frac{\varphi_n(x_m)}{q(x_m)}\Bigr)^2 ,
$$

而 $\mathbb E\,\mathrm{tr}\,G=N$，也就是说这 $M$ 个行范数平方**加起来平均是 $N$**。稳定性定理要求的正是这些行范数不要有异常值：若对某个 $r>0$

$$
\frac{M}{\log M}\ \ge\ C(r+1)\,\sup_{x\in D}\sum_{n=1}^{N}\Bigl(\frac{\varphi_n(x)}{q(x)}\Bigr)^2 ,
\qquad C=\frac{2}{\log(27/8e)}\approx 9.24 ,
$$

则以不小于 $1-2M^{-r}$ 的概率有 $\|G-I\|_2\le\frac12$。右端那个上确界称为**稳定因子**，它的最小可能值就是 $N$。于是整个专题可以一句话概括：**稳定性定理问的是「最坏的一行比平均的一行差多少」，而所有采样设计的花招都是在压这个比值。** 行范数处处相等时，稳定因子取到下界 $N$，样本预算就是 $M\gtrsim N\log N$；某处有尖峰时，尖峰有多高，预算就贵多少倍。

### Christoffel 函数与本站的记法约定

$q\equiv1$ 时稳定因子就是

$$
K(z)=K_\Lambda(z)=\sum_{\alpha\in\Lambda}\varphi_\alpha^2(z)
$$

的上确界 $\|K\|_\infty$。写成向量形式 $K(z)=\varphi^{\mathsf T}\varphi$ 就能看出：任何正交变换 $\psi\leftarrow U\varphi$ 都不改变 $K$，因此**它是子空间 $\mathbb P_\Lambda$ 的性质，与基的选取无关**。这也是它值得被当作中心量的原因。

反过来，取 $q^2=K/N$，即从密度

$$
\rho(x)=\frac{1}{N}\sum_{n=1}^{N}\varphi_n^2(x)\,w(x)
$$

采样并以 $1/q^2=w/\rho$ 加权，稳定因子就**恒等于** $N$——行范数被强行拉平。$\rho$ 确实是概率密度，因为 $\int_D\frac1N\sum_n\varphi_n^2\,w\,\mathrm dx=\frac1N\sum_n\|\varphi_n\|^2=1$。这就是诱导采样，本专题另一支（编号 22、45）的全部内容。

> [!warning] 记法约定
> 本站统一取 $K(z)=\sum_\alpha\varphi_\alpha^2(z)$，把（归一化的）Christoffel 函数取为 $N/K(z)$，因此本页出现的所有权都写成 $1/K$ 的形式。原文中的 $\lambda_\Lambda(z)$ 即 $1/K_\Lambda(z)$。这一族论文的记法并不统一：编号 22 与 28 用 $K$ 记 $\sum_\alpha\varphi_\alpha^2$（即 Christoffel 函数的倒数），编号 24 用 $\lambda_\Lambda$ 记 $1/\sum\varphi_i^2$（即 Christoffel 函数本身），编号 36 用 $K(\xi)$ 记 $N/\sum_j\Phi_j^2$，编号 44 与 45 用 $\kappa$、$q^2$ 记归一化的倒数。**引用任何一篇原文之前都要核对它自己的约定**，否则权会取成倒数。

### Gauss 权就是 Christoffel 函数值

这是本页最该先说的一件事，因为编号 13 与 14 的构造完全建立在它之上。在一维、$\Lambda=\{0,\dots,N-1\}$ 的情形，存在一族 $N$ 点正权求积公式，对次数不超过 $2N-2$ 的多项式精确，而它的权恰恰是 Christoffel 函数值：

$$
\int_\Gamma p(z)\rho(z)\,\mathrm dz=\sum_{z\in A_N}\frac{1}{K_\Lambda(z)}\,p(z),
\qquad \deg p\le 2N-2 ,
$$

Gauss 点是其中 $y\in\phi_N^{-1}(0)$ 的特例。**于是「在 Gauss 网格上做随机子采样并按 Gauss 权加权」不是一个凑出来的方案，它本身就是 Christoffel 加权采样。** 编号 13 与 21 用的正是这一点，而编号 22 与 45 后来把这个隐含的做法显式化，直接把采样密度取成 $K/N$ 倍的正交测度。

### 三条采样路线的预算

把上面的判据 $M/\log M\gtrsim\|K\|_\infty$ 代入具体的基，就得到本专题反复出现的三档预算：

| 采样方式                         | 稳定因子的增长                            | 样本要求                   |
| -------------------------------- | ----------------------------------------- | -------------------------- |
| 均匀测度 + 张量 Legendre，不加权 | $\lVert K\rVert_\infty\sim N^2$           | $M\sim N^2$                |
| Chebyshev 测度，不加权           | $\lVert K\rVert_\infty\sim N^{\ln3/\ln2}$ | $M\sim N^{1.585}$          |
| Christoffel 加权（诱导采样）     | 恒为 $N$                                  | $M\gtrsim N\log N$，非渐近 |

第三行的判据只依赖 $N=\dim\mathbb P_\Lambda$，与维数 $d$、区域 $D$、权 $w$ 以及具体取哪个 $N$ 维子空间都无关；代价是必须能从 $\rho=q^2w$ 采样，而这个密度依赖 $(V,w,D)$。贪心确定性选点（近似 Fekete 点）走的是第四条路，把 $M$ 进一步压向 $N$，那属于[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|最优采样与预条件]]那一页。

> [!warning] 高斯情形的渐近密度是猜想
> 有界区域上诱导测度的大 $N$ 极限是（张量化的）Chebyshev 密度，这是定理。**高斯情形不是。** $D=\mathbb R^d$ 配高斯权时，$\lim_{k\to\infty}\rho(x/\sqrt k)=C(2-\|x\|_2^2)^{d/2}$ 在编号 22、28、36、45 中一律以**猜想**的身份出现，从未被证明；编号 45 直接说多元情形「更多的是未知，很多情况下我们只有猜想」。本页与相邻各页都不把它当作已证结论引用。

## 6：先用实验回答「取哪些点」

### 直觉

编号 6 是整个采样设计纲领的经验序幕。当时的实际担忧很具体：**随机点会不会破坏收敛率**。这个担忧有直观基础——结构化网格（稀疏网格）是被设计出来的，随机点则可能挤成一团、留下空洞，看上去理应更差。论文没有从理论上回答，而是把「哪种点更好」拆成三条彼此独立、通常被分开研究的判据，一起测：收敛率、稳定性（设计矩阵的条件数）、以及函数值带数值噪声时的稳健性。

它的贡献不是一个数学对象，而是**这套比较协议本身**。后来的进展说明这个提法只对了一半：真正的问题不是收敛率，而是要让设计矩阵条件良好需要多少个样本——正是上一节那个稳定因子。

### 问题设定

多项式空间上的离散最小二乘，三族候选设计点：稀疏网格（SG）、Monte Carlo（MC）、拟 Monte Carlo（QMC）。论文中的「设计矩阵」按标准形式应为 $A_{ij}=\phi_j(y^{(i)})$，即基函数在设计点上取值构成的 Vandermonde 型矩阵——**摘要只说「设计矩阵」而没有写出它，这个形式是本站按上下文补的**，与后面几篇使用的加权形式（见上一节）不同，那里多了 $\sqrt{w}$ 的行缩放。

### 数值证据

测例是若干经典高维测试函数加上一个随机 ODE 模型。三条结论：

| 点集 | 收敛                                               | 稳定性                   |
| ---- | -------------------------------------------------- | ------------------------ |
| MC   | 不引入低收敛率；函数有足够正则性且点够多时仍是高阶 | —                        |
| QMC  | 同上；在较高维上是好选择                           | 确定性，条件数方面也更好 |
| SG   | 只在很低维（论文原话是 $d\le2$）时收敛更好         | —                        |

「点够多」这个判断在论文里是定性的，没有给出 $M$ 与 $N$ 的任何关系式——**这恰恰是编号 9、11、13、22、28、45 后来要补的洞**。另外值得一提：论文列出的三条判据里，第三条（函数值带噪声时的稳健性）在摘要给出的三条结论里没有对应的一条，本站因此不知道它的结论是什么。

> [!note] 这一篇只到摘要层面
> 本站未取得该文全文。具体维数（除了摘要里的 $d\le2$）、误差量级、条件数数值、以及每条结论背后的测例配置都未核实，因此本页不报告任何数字，也不给它编造定理或实验表格。

### 与其他论文的关系

它问的是「取哪些点」，答案靠实验。此后每一篇都在回答同一个问题的加强版，而且是用**构造出来的**点集而非**挑选出来的**点集：编号 9 的新配点网格、编号 11 无界区域上的随机求值、编号 13 与 21 的随机求积、编号 14 的非结构网格，以及[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|另一页]]上的 Christoffel 加权采样与加权近似 Fekete 点。

## 9：确定性点集与二次样本复杂度

### 直觉

在此之前的最小二乘稳定性理论（Migliorati–Nobile–von Schwerin–Tempone；Cohen–Davenport–Leviatan）都从正交测度独立采样，结论无可避免地是概率性的：「以高概率成立」，对**任何一次具体的实现**都不给保证。这篇问的是能不能去掉这个限定词。

它的机制值得单独说一句。概率方法说的是「$M$ 项随机和以高概率是 $O(\sqrt M)$」；数论里的**平方根消去**说的是「某个结构化的指数和**永远**不超过 $(d-1)\sqrt M$」。后者比前者强，因为它没有例外集。把 Weil 的指数和界用在 Gram 矩阵的非对角元上，得到的不是谱范数的集中，而是**逐元的、确定的**上界，进而是严格对角占优——这是一个比「谱范数以高概率接近 1」更强、且对每一次实现都成立的陈述。

### 问题设定

$Y=(Y_1,\dots,Y_d)^{\mathsf T}$ 在 $\Gamma\equiv[-1,1]^d$ 上相互独立，边缘密度 $\rho_i$，联合密度 $\rho(Y)=\prod_{i=1}^{d}\rho_i(Y_i)$；在 $L^2_\rho$ 中逼近 $f:\Gamma\to\mathbb R$，范数 $\|f\|_{L^2_\rho}=\bigl(\int_\Gamma f^2\rho\,\mathrm dY\bigr)^{1/2}$；多项式空间 $\mathbb P_\Lambda=\mathrm{span}\{\Phi_n\}_{n\in\Lambda}$ 取张量积（TP）或全次数（TD）形式，$N=\#\Lambda$。

**新网格（Weil 点）。** 对素数 $M>2q+1$，

$$
\Theta_M:=\Bigl\{\,y_j=\cos(x_j)\ :\ x_j=2\pi\bigl(j,\,j^2,\dots,j^d\bigr)/M,\quad j=0,\dots,\lfloor M/2\rfloor\,\Bigr\},
$$

即 $x_j$ 的第 $q$ 个分量是 $2\pi j^q/M$。点数为 $m+1$，$m=\lfloor M/2\rfloor$；论文指出 $\{y_j\}_{j=0}^{m}$ 与 $\{y_j\}_{j=m+1}^{M}$ 重合，所以取一半就够。这个构造是把 Xu 用于稀疏三角多项式确定性采样、并推广到高维稀疏 Chebyshev 多项式的 Weyl 和点集，做了一次余弦移植。

设计矩阵取不加权的离散内积：

$$
A=\bigl(\langle\Phi_i,\Phi_j\rangle_m\bigr)_{1\le i,j\le N},
\qquad
\langle u,v\rangle_m=\sum_{k=0}^{m}u(y_k)v(y_k) .
$$

### 推导

**第一步：Weil 公式（定理 2.1）。** 设 $M$ 为素数，$f(x)=m_1x+m_2x^2+\cdots+m_dx^d$。若存在 $1\le j\le d$ 使 $M\nmid m_j$，则

$$
\Bigl|\sum_{j=0}^{M-1}e^{\frac{2\pi i f(j)}{M}}\Bigr|\ \le\ (d-1)\sqrt{M} .
$$

余弦移植的作用就在这里：$y_j=\cos(x_j)$ 使 Chebyshev 基在这些点上的取值化为整数倍角的余弦，两个基函数的乘积化成指数和，而指数中的 $j$ 的次数由多重指标决定。$n\ne k$ 时多重指标之差有分量不被 $M$ 整除（$M>2q+1$ 保证了这一点），Weil 公式适用；$n=k$ 时指数恒为零，求和退化为计数。

**第二步：引理 3.1 的两个界。** 由此得

$$
\Bigl|\sum_{j=0}^{m}\Phi_n(y_j)\Phi_k(y_j)\Bigr|\ \le\ \frac{(d-1)\sqrt{M}+1}{2}\quad(n\ne k),
\qquad
\sum_{j=0}^{m}\Phi_n^2(y_j)\ \ge\ \frac{M}{2^{d+1}}-\frac{(d-1)\sqrt{M}}{2} .
$$

对角项的主项 $M/2^{d+1}$ 正是 $\cos^2$ 在每一维上的平均值 $1/2$ 乘以 $m+1\approx M/2$ 个点，逐维相乘的结果；非对角项则整个由平方根消去控制。

**第三步：Gershgorin。** 把 $A$ 按 $2^{d+1}/M$ 缩放后，非对角元的上界正好是

$$
\delta=\frac{2^{d+1}}{M}\cdot\frac{(d-1)\sqrt{M}+1}{2}=\frac{2^{d}\bigl((d-1)\sqrt{M}+1\bigr)}{M} ,
$$

而缩放后的对角元落在 $1$ 的 $\delta$ 邻域内。Gershgorin 定理给出 $|\lambda_i-1|\le N\delta$，剩下的只是验证 $M\ge4^{d+1}d^2N^2$ 蕴含 $N\delta\le\frac12$：此条件即 $\sqrt M\ge2^{d+1}dN$，代入得 $N\delta\lesssim 2^d(d-1)N/\sqrt M\le\frac{d-1}{2d}<\frac12$。**整篇的定量结论就压在这一行不等式上。**

**第四步：推广到其他测度的权（4.2 节）。** 加权最小二乘为 $f_\Lambda=\arg\min_{v\in\mathbb P_\Lambda}\sum_{i=0}^{m}w_i(f(y_i)-v(y_i))^2$。因为 $\Theta_M$ 渐近等分布到 Chebyshev 测度（定理 4.3），不加权的离散范数模拟的是 Chebyshev 范数；要模拟 $\rho$ 加权的范数，就取

$$
w_i=\frac{\rho(y_i)}{\rho_c(y_i)}=\pi^d\rho(y_i)\prod_{q=1}^{d}\bigl(1-(y_i^q)^2\bigr)^{1/2} ,
$$

均匀密度 $\rho\equiv2^{-d}$ 时化为 $w_i=(\pi/2)^d\prod_{q=1}^{d}\bigl(1-(y_i^q)^2\bigr)^{1/2}$。由于 $w_i$ 乘的是一个二次型，实际效果是按 $\sqrt{w_i}$ 做预条件，即 $\Phi_i(y)\mapsto\prod_{q=1}^{d}(1-(y^q)^2)^{1/4}\Phi_i(y)$——**论文自己指出这正是 Legendre 逼近下使 $\ell_1$ 极小化的设计矩阵条件良好的那个预条件**。$\rho\propto\rho_c$ 时权是常数。注意这个 $w_i=\rho/\rho_c$ 与共用装置一节里的换测度权 $w(z_s)/v(z_s)$ 是同一个东西。

### 定理

- **引理 3.1。** 假设：$M$ 为素数且 $M>2q+1$，$y_j=\cos(x_j)$ 如上。结论：上面两个界。
- **定理 3.2（稳定性）。** 假设：$M\ge4^{d+1}d^2N^2$ 且 $M$ 为素数。结论：谱范数下 $\bigl|\!\bigl|\!\bigl|\frac{2^{d+1}}{M}A-I\bigr|\!\bigr|\!\bigr|\le\frac12$。**这就是二次样本复杂度 $M\gtrsim N^2$，并带一个显式的、依赖维数的前因子 $4^{d+1}d^2$。**
- **推论 3.3（唯一性）。** 同样假设下，$\sum_{k=0}^{m}(p(y_k)-f(y_k))^2$ 在 $\mathbb P_\Lambda$ 上的极小元唯一，因为 $A$ 严格对角占优故非奇异。
- **定理 3.4（Chebyshev 测度下的最优收敛）。** 设 $P^\Lambda f=\arg\min_{p\in\mathbb P_\Lambda}\|f-p\|_{L^2_{\rho_c}}$，$P^\Lambda_m f$ 为 $\Theta_M$ 上的离散最小二乘解。假设 $M\ge4^{d+1}d^2N^2$ 为素数，则

  $$
  \|f-P^\Lambda_m f\|_{L^2_{\rho_c}}\ \le\ \Bigl(1+\frac{4}{d^2 N}\Bigr)\,\|f-P^\Lambda f\|_{L^\infty} .
  $$

  因子 $1+4/(d^2N)$ 随 $N$ 增大趋于 $1$，故离散投影渐近上与最佳 $L^\infty$ 逼近同样好。**这个估计是确定性的，不带「以高概率」。** 这是二次样本数换来的东西。

- **定义 3.5 与推论 3.6（推广到其他测度）。** 若存在与 $Y$ 无关的常数 $C$ 使 $0<\rho(Y)\le C\rho_c(Y)$ 对一切 $Y\in\Gamma$ 成立（称 $\rho$ 被 Chebyshev 密度支配），则对任意 $f\in L^2_{\rho_c}$，

  $$
  \|f-P^\Lambda_m f\|_{L^2_\rho}\ \le\ \sqrt{C}\,\Bigl(1+\frac{4}{d^2 N}\Bigr)\|f-P^\Lambda f\|_{L^\infty} ,
  $$

  仍要求 $M\ge4^{d+1}d^2N^2$ 为素数。注 3.8 指出这覆盖均匀测度以及一切满足 $0<\rho_{\min}\le\rho\le\rho_{\max}$ 的测度；注 3.7 指出它对**认知不确定性**有用：$Y$ 的密度未知时，只要该未知密度满足支配条件，基于 Chebyshev 的逼近仍然有效。

- **定理 4.3（渐近等分布）。** 设 $M_K$ 为第 $K$ 个素数，$m_K=\lfloor M_K/2\rfloor+1$，经验测度 $\nu_K:=\frac{1}{m_K}\sum_{j=1}^{m_K}\delta(y_{j,K})$，$\nu_c$ 为归一化 Chebyshev 测度（密度 $\rho_c(y)=\pi^{-d}\prod_q(1-y_q^2)^{-1/2}$）。则 $\nu_K\to\nu_c$ 弱收敛。证明把 Weil 公式送进 **Weyl 等分布判别法**（定理 4.1、推论 4.2）。**这条结论才是权公式的许可证**：没有它，$w_i=\rho/\rho_c$ 的推导无从谈起。

### 数值实验

第 5 节把 $\Theta_M$ 与 Monte Carlo 网格在 TP 与 TD 两种空间中对比，确定性结果画点、Monte Carlo 结果画方块。

| 图      | 设定                                            | 对照                                                                  | 观察                                     |
| ------- | ----------------------------------------------- | --------------------------------------------------------------------- | ---------------------------------------- |
| 图 1    | $d=2$，$M=997$ 的点分布                         | Chebyshev 测度下的一次 Monte Carlo 实现                               | 两者都在边界附近聚集，与定理 4.3 一致    |
| 图 2、3 | $d=2$ 的 TD 与 TP 空间，另有一个 $d=3$ 的 TD 例 | 线性规则 $m=c(\#\Lambda)$ 与二次规则 $m=c(\#\Lambda)^2$，取若干个 $c$ | 设计矩阵条件数；确定性点与随机点表现相当 |

论文自己在摘要里给的定性结论就是「确定性点的表现与随机点相似」——**收益不在精度，而在去掉概率限定词**。它也坦承所取的 $d$、$q$ 与测试函数并不特殊，换成别的也差不多。

实验与定理之间有一个应当点破的距离：定理要求的前因子在 $d=2$ 时是 $4^3\cdot4=256$，在 $d=3$ 时是 $4^4\cdot9=2304$，而实验跑的是 $m=c(\#\Lambda)$ 与 $m=c(\#\Lambda)^2$（$c$ 的具体取值本站未核实）。除非 $c$ 取到这个量级，实验所在的区间就根本不在定理的覆盖范围内，而条件数照样是好的。**这说明定理 3.2 是充分条件而非紧的条件**，这一点论文没有展开。

### 与其他论文的关系

**这篇确立了这一族的模板**：一个显式的点集或密度，一条量化所需样本数的稳定性定理，一条近最佳逼近的推论。它的 $M\gtrsim N^2$ 是后面几篇要攻击的基准——编号 22 与 45 证明从 Christoffel 加权／诱导密度采样可把要求降到 $M\sim N\log N$，编号 28 干脆用加权近似 Fekete 点取代随机采样。这里引入的权 $w_i=\rho(y_i)/\rho_c(y_i)$ 是编号 22 的 Christoffel 权与编号 24、32 的预条件的直系祖先。Xu 的 Weil 和装置在编号 10（确定性插值点）与编号 29（稀疏三角恢复）中再次出现，见[[computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce|稀疏恢复与数据驱动混沌]]。编号 6 提出的经验问题，在这里第一次得到理论回答。

## 11：无界区域把样本复杂度从多项式推到指数

### 直觉

已有的稳定性理论都针对**有界**参数域：$[-1,1]^d$ 上均匀或 Chebyshev 测度，$M\sim N^2$ 够用（Chebyshev 更少）。高斯与 Gamma 参数生活在无界域上，而论文指出那里的情况**在质上更差，不只是常数变坏**。

为什么是质的差别，可以从共用装置那条判据 $M/\log M\gtrsim\|K\|_\infty$ 直接读出来。有界域上这个上确界是在**紧集**上取的：基函数连续，上确界有限，而且随 $N$ 只是多项式地长（Legendre 是 $N^2$，Chebyshev 是 $N^{1.585}$）。无界域上取上确界的集合不再紧，于是有两件事同时发生且互相加剧：正交多项式在远处没有一致界，而高斯测度在远处几乎没有质量。**样本几乎全都落在基函数很小的地方，偶尔一个落到远处的样本就带来一个巨大的行范数**——行范数极不均匀，正是稳定因子爆炸的形态。这不是把常数放大若干倍能补救的。

这个诊断有一个内部证据：论文的引理 3.1 给出的一致衰减界只对 Hermite **函数**成立，对 Hermite **多项式**不成立，而两者只差一个 $e^{-y^2/2}$ 因子。**被这个因子治好的病，正是多项式在远处的无界增长。**

由此看，论文的修法可以概括成一句话：**把权从测度里搬进基里**。权一旦在基里，基就有一致界，稳定因子就受控；采样密度于是可以为「覆盖」而非为「正交」来选，选成映射后的均匀分布。这与编号 22、45 的 Christoffel 加权是同一件事的两种做法——后者保留多项式基，改密度并加权。

第二个独立的问题是 Hermite 展开的分辨率很差：论文直接引用 Gottlieb 与 Orszag 的说法——分辨 $\sin(x)$ 的 $M$ 个波长需要接近 $M^2$ 个 Hermite 多项式。这个问题不是稳定性，而是收敛速度，论文用另一件装置（缩放因子）来对付。

> [!warning] 本节的记号与编号 9 冲突
> 这一篇里 $m$ 是**随机点个数**，$K$ 是**基函数个数**，而 $M$ 是被逼近函数的有效支集半径（$|f(y)|<\epsilon$ 当 $|y|>M$）。编号 9 一节里 $M$ 是素数、$N$ 是空间维数。读定理时务必分清。

### 问题设定

标准化 Hermite 多项式 $\{H_k\}$ 在 $\mathbb R$ 上对 $\rho_G(y)=e^{-y^2}$ 标准正交，$\int_{-\infty}^{+\infty}\rho_G(y)H_m(y)H_n(y)\,\mathrm dy=\delta_{mn}$；Laguerre 多项式对 $\rho_E(y)=\prod_{i=1}^{d}e^{-y_i}$（Gamma 的指数特例，一般 Gamma 密度 $\rho_E(y)=\beta^\alpha y^{\alpha-1}e^{-\beta y}/\Gamma(\alpha)$ 用广义 Laguerre 混沌处理）。稳定性用 $\mathrm{cond}(A)=\sigma_{\max}(A)/\sigma_{\min}(A)$ 度量。正规方程为 $\bigl(\langle\Phi_i,\Phi_j\rangle_m\bigr)_{i,j=1,\dots,N}$，右端 $f=D^{\mathsf T}b=\bigl(\langle f,\Phi_j\rangle_m\bigr)_j$，可用设计矩阵的 QR 或正规方程的 Cholesky 求解。

### 推导

**诊断（2.3 节）。** 从高斯／Gamma 测度采样、用对应的多项式混沌基，条件数在线性规则 $m=c(\#\Lambda)$ 与二次规则 $m=c(\#\Lambda)^2$ 下都随多项式阶数**指数增长**。论文把所需样本数写成

$$
m=(\#\Lambda)^{\,c\,\#\Lambda} ,
$$

并称其在实际计算中不可接受。**这是论文对数值观察到的标度的自己的刻画，不是一条被证明的下界**——本页不把它当作定理引用。

**修法之一：用 Hermite／Laguerre 函数而不是多项式（式 (3.1)–(3.3)）。**

$$
\tilde H_m(y)=e^{-\frac{y^2}{2}}H_m(y),\qquad
\tilde L_m(y)=e^{-\frac{y}{2}}L_m(y),\qquad m=0,1,\dots
$$

它们对 **Lebesgue 测度**标准正交，$\int_{-\infty}^{+\infty}\tilde H_m\tilde H_n\,\mathrm dy=\delta_{mn}$——高斯权被吸收进了基里。多元情形做张量化。它们不再是多项式，但论文仍把指标 $q$ 叫作「多项式阶」。为什么衰减基对 UQ 是对的（式 (3.4)）：感兴趣量的形式是 $\mathrm{QoI}=\int_\Gamma\rho(y)(g\circ f)(y)\,\mathrm dy$，即使 $g\circ f$ 本身不衰减，只要它增长慢于高斯，$\rho\cdot(g\circ f)$ 就衰减，因此改为逼近 $\tilde f(y)=\rho(y)(g\circ f)(y)$。

**修法之二：映射后的均匀随机点（式 (3.7)）。** 在有界区间上均匀取点，再用带映射参数 $L$ 的映射送到无界域：

$$
y(\xi)=\begin{cases}\dfrac{L}{2}\log\dfrac{1+\xi}{1-\xi}, & r=0,\\ \dfrac{L\xi}{\sqrt{1-\xi^2}}, & r=1,\end{cases}
\qquad
\xi(y)=\begin{cases}\tanh\bigl(\tfrac{y}{L}\bigr), & r=0,\\ \dfrac{y/L}{\sqrt{y^2/L^2+1}}, & r=1.\end{cases}
$$

$r=0$ 的**对数映射**使变换后的点指数衰减，用于高斯测度；$r=1$ 的**代数映射**用于 Gamma 测度。$L$ 是控制条件数的旋钮。

**修法之三：为收敛速度引入缩放因子（式 (3.28)、(3.32)、(3.33)）。** 设 $f$ 指数衰减，$|f(y)|<\epsilon$ 当 $|y|>M$。把展开写成

$$
f(y)=\sum_{n=0}^{K-1}c_n\tilde H_n(\alpha y)
\ \Longleftrightarrow\
f\Bigl(\frac{y}{\alpha}\Bigr)=\sum_{n=0}^{K-1}c_n\tilde H_n(y),
\qquad \alpha>0,
$$

使缩放后的点 $y_i/\alpha$ 落在 $f$ 的有效支集内。系数按 $f_k=\langle f,H_k\rangle_m=\sum_{i=1}^{m}f(y_i/\alpha)H_k(y_i)$ 计算，于是要求 $\max_j|y_j|/\alpha\le M$，给出朴素规则 $\alpha=\max_{1\le j\le m}\{|y_j|\}/M$。但点是随机的，少数几个极端「坏点」会把 $\alpha$ 撑得过大，所以论文的**准最优缩放**把它们丢掉：

$$
\tilde\alpha=\max_{1\le j\le\tilde m}\{|y_j|\}/M,
\qquad \tilde m=\lfloor\mu m\rfloor,
$$

$\mu$ 接近 1，实践中取 $\mu\approx98\%$，即丢掉最极端的 2% 样本。**这个截尾最大值是本文最有实用价值的一处贡献**，数值实验里它也确实是成败所在。

**稳定性证明的路线。** 走矩阵 Chernoff：把缩放后的设计矩阵写成 $\hat A=X_1+\cdots+X_m$，其中 $X$ 是

$$
X=\frac{L}{m}\bigl(\tilde H_i(y)\tilde H_j(y)\bigr)_{i,j=0,\dots,K-1}
$$

的独立同分布拷贝，$y$ 是变换后的均匀随机变量，再用独立正自伴随机矩阵的 Chernoff 界，该界需要几乎必然的上界 $\lambda_{\max}(X_i)=|\!|\!|X_i|\!|\!|\le R$。注意 $X$ 是秩一的，所以

$$
\lambda_{\max}(X)=\frac{L}{m}\sum_{i=0}^{K-1}\tilde H_i^2(y) ,
$$

**这正是共用装置一节里的行范数平方**。假设 $L>3\tau$ 与 $L>5\sqrt K$ 就是在这里用掉的：前者让引理 3.1 的衰减界生效，后者保证映射区间宽到能盖住前 $K$ 个基函数的活动范围。分析只对一维 Hermite 情形写出，Laguerre 情形论文说「可直接照搬」——**本站未见其写出**。

论文在结论里强调的两参数结构值得记住：$L$ 控制**稳定性**，$\alpha$ 控制**收敛速度**，两者为不同目的调节。

### 定理

- **引理 3.1。** 对任意整数 $K$，存在常数 $\tau>0$，使得当 $|y|>\tau$ 时对一切 $0\le k\le K-1$ 有 $|\tilde H_k(y)|\le|y|^{-3/2}$。（成立的原因是 $e^{-y^2/2}$ 因子迫使 $|\tilde H_k(y)|\cdot|y|^t\to0$ 对每个 $t>0$ 成立。）**关键是这个界对 $k$ 一致。**
- **定理 3.2（稳定性）。** 用 Hermite 函数 (3.1) 与变换后的均匀随机点 (3.7)，缩放设计矩阵 $\hat A=LA$ 满足：对每个 $r>0$，

  $$
  \Pr\Bigl\{\,\bigl|\!\bigl|\!\bigl|\hat A-I\bigr|\!\bigr|\!\bigr|\ \ge\ \tfrac58\,\Bigr\}\ \le\ 2m^{-r} ,
  $$

  **前提是**

  $$
  K\ \le\ \kappa\,\frac{m}{\log m},
  \qquad
  \kappa:=\frac{4c_{1/2}}{3(1+r)},
  \quad c_{1/2}=\frac12+\frac12\log\frac12>0 ,
  $$

  且映射参数满足 $L>\max\{3\tau,\ 5\sqrt K\}$，其中 $m$ 是随机点个数、$K$ 是基函数个数。**这就是对数线性样本复杂度 $K\lesssim m/\log m$，等价于 $m\gtrsim K\log K$（相差常数）**，相对多项式混沌形式的指数要求是质的改进，形式上也与有界域的最好结果一致。

  换算一下常数便于估量：$c_{1/2}=\frac12(1-\ln2)\approx0.153$，故 $\kappa\approx0.205/(1+r)$，$r=1$ 时 $K\lesssim0.10\,m/\log m$。**定理并不廉价**，它要求的样本数比朴素读法给人的印象大一个量级左右。

条件 $L>5\sqrt K$ 里的 $\sqrt K$ 与编号 14 定理 14 中无界域收缩因子 $k^{-1/r}$（高斯情形 $r=2$）在量级上是同一件事，这不是巧合：两处说的都是「$K$ 个基函数的活动范围随 $\sqrt K$ 膨胀」。

### 数值实验

条件数一律用 $\mathrm{cond}(A)=\sigma_{\max}(A)/\sigma_{\min}(A)$，并对 **100 次独立实现取平均**，因为矩阵是随机的。

**稳定性部分。**

| 图   | 基与采样                     | 维数与采样规则                                                 | 观察                                 |
| ---- | ---------------------------- | -------------------------------------------------------------- | ------------------------------------ |
| 图 1 | Hermite 多项式混沌，高斯采样 | 一维，$m=c(\#\Lambda)$ 与 $m=c(\#\Lambda)^2$                   | 条件数随多项式阶数指数增长——负面结果 |
| 图 4 | Laguerre 函数，映射均匀点    | 一维 $m=30(\#\Lambda)$；二维 $m=6(\#\Lambda)^2$，TD 与 TP 空间 | 对比映射参数 $L=8$ 与 $L=64$         |

**收敛部分。**

| 算例       | 目标函数                                       | 关键参数                                                                                      | 结果                                                                        |
| ---------- | ---------------------------------------------- | --------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| 图 6       | $f(y)=2^{-6y^2}$（快速衰减）                   | $L=8$ 固定以保证稳定                                                                          | $\alpha=1$（不缩放）收敛很慢但稳定；加缩放显著加速，该例最好约 $\alpha=2.8$ |
| 一维第二例 | $\tilde f(y)=2^{-0.2y^2}$                      | $M=16$                                                                                        | —                                                                           |
| 二维例     | $\tilde f(y)=e^{-4(y_1^2+y_2^2)}\sin(y_1+y_2)$ | $M=2.5$，TD 空间；$m=10(\#\Lambda)$ 对 $m=2(\#\Lambda)^2$；不缩放 / $\mu=0.985$ / $\mu=0.980$ | 截尾缩放 (3.33) 收敛快，朴素最大值缩放 (3.32) 几乎与不缩放无异              |

一致的发现是：**决定成败的是丢掉坏点，而不是缩放本身**。参数化 UQ 的应用用了一个随机 ODE 模型和一个带对数正态随机输入的椭圆问题。另有一处报告的最优缩放约为 $\alpha=0.82$（本站未核实它对应哪一个算例）。

实验支持了什么、又差在哪里，值得分开说。它确实支持了两件事：多项式混沌在无界域上的失败是真实的、可测量的；换基加映射后条件数受控，且 $L$ 是有效的旋钮。它没有支持的是收敛那一半——**关于 $\alpha$ 的选取，本站在原文中未见任何定理，摘要用的词是「讨论」**，而两个报告出来的最优值（$2.8$ 与 $0.82$）相差三倍以上，说明它们是算例相关的观察，不是通用的取值建议。此外本站未核实任何误差量级，因此不给误差表。

### 与其他论文的关系

这是有界域稳定性理论的无界域对应物，也是本专题第一篇达到**对数线性**而非二次样本要求的论文——但它是靠同时更换**基**（多项式→函数）与**采样密度**（高斯→映射均匀）买到的，而不是靠加权。**这条路线与后来的 Christoffel 加权是两个不同的选择**：编号 22 与 45 保留多项式基，改密度并加权，同样达到对数线性；代价不同——前者放弃了多项式基的正交结构，后者需要计算 Christoffel 函数。无界域的关切在编号 24 与 32 中再次出现，那里的预条件框架明确宣称覆盖「有界或无界」区域。本文引用的有界域二次基准 $m\sim(\#\Lambda)^2$ 正是编号 9 为确定性点集建立的那一条。

## 13：随机求积——一篇无法核实的论文

### 直觉

设计空间的两端是清楚的：一端是从正交测度独立采样，描述简单但样本要得多；另一端是完整的张量积 Gauss 网格，离散正交性完美但基数 $\prod_i n_i$ 关于维数呈指数。中间的想法很自然——**从 Gauss 网格里随机抽一小部分并重新加权**，保住 Gauss 求积那套优良的离散正交结构，而只付 $M$ 次昂贵求值的代价。

它之所以不是一个任意的重加权方案，全靠共用装置一节里那个恒等式：**Gauss 权就是 Christoffel 函数值**，所以子采样 Gauss 网格隐含地就是 Christoffel 加权采样。

### 问题设定

以下构造取自设计完全相同的姊妹篇编号 21（它把同一套采样方案用于 $\ell_1$ 恢复并给出分析）；**本文特有的式号与变体本站未核实**。

张量积 Gauss 网格：设 $\Theta^i_{n}=\{z^i_1,\dots,z^i_n\}\subset\Gamma_i$ 是第 $i$ 维的 $n$ 点 Gauss 点集，张量集为 $\Theta_{\mathbf n}=\Theta^1_{n_1}\otimes\cdots\otimes\Theta^d_{n_d}$，$|\Theta_{\mathbf n}|=\prod_{i=1}^{d}n_i$。Gauss 权：

$$
w_{\mathbf k}=\lambda_{\mathbf n}(z_{\mathbf k})
=\prod_{i=1}^{d}\lambda^i_{n_i}(z^i_{k_i})
=\prod_{i=1}^{d}\frac{1}{\sum_{k=0}^{n_i-1}\bigl[\phi^i_k(z^i_{k_i})\bigr]^2} .
$$

网格上的均匀经验测度 $\nu_{\mathbf n}=\bigotimes_{i=1}^{d}\nu^i_{n_i}=\frac{1}{\prod_i n_i}\sum_{\mathbf k\le\mathbf n}\delta_{z_{\mathbf k}}$，从 $\nu_{\mathbf n}$ 独立抽样就是在张量 Gauss 网格上均匀取点。加权设计矩阵为 $D=\sqrt W\Psi$，$(\Psi)_{m,n}=\varphi_n(x_m)$，$(W)_{m,m}=w_m>0$，加权问题把 $\Psi c=f$ 换成 $Dc=\sqrt W f$。

索引元数据（可核实）：关键词为最小二乘、不确定性量化、正交多项式、广义多项式混沌；MSC 为 65C20、65D30、41A10、65D15、60G99。其**引用文献**（不是它的结论）透露了技术装置：Nevai 关于广义 Jacobi 权、Christoffel 函数与 Jacobi 多项式的两篇，一篇关于 Gauss 求积权渐近展开的文章，Tang–Iaccarino 关于子采样 Gauss 求积点估计多项式混沌，Cohen–Davenport–Leviatan，Migliorati 等人关于随机求值下的离散 $L^2$ 投影，以及**四本实验设计专著**。Christoffel 函数那几条印证了「Gauss 权的渐近就是 Christoffel 函数值」是它的分析主干；实验设计那四本则指向一条相邻各篇都没有的、与最优实验设计的显式联系。

### 无法核实的部分

> [!warning] 这一篇的定理、常数与数值结果本站一概不报告
> 出版社封锁正文，researchr 明确记录「摘要缺失」，Crossref、OpenAlex、MaRDI/zbMATH、Semantic Scholar 均无摘要（Semantic Scholar 的记录直接注明摘要字段被出版方删除），ScienceDirect 对正文与录用稿两个地址都返回 403。**本站未能从论文本身核实任何一条定理陈述、任何一个常数、任何一条显式的样本复杂度关系，也未核实任何数值实验。**
>
> 「样本数关于多项式维数**线性**增长」这一说法的唯一依据是第三方表述：Seshadri、Narayan 与 Sarkar 在 _Quadrature Strategies for Constructing Polynomial Approximations_ 中写道，Zhou 等人随机子采样各行并证明了最小二乘问题在 $m$ 关于 $n$ 线性增长时的稳定性。该文作者之一也是本文作者之一，所以这是可信的旁证，**但它不是本文自己的措辞，其精确假设、概率陈述与常数一律未核实**。任何下游引用都不应把某个具体常数或概率界归到这一篇名下。

从姊妹篇编号 21 可以核实的、最接近的一条是：同一采样方案的 $\ell_1$／压缩感知版本，充分样本数为 $M\ge L(\mathbf n)C_1s\log^3(s)\log(N)$；而对超立方体上的均匀随机变量，Tang 与 Iaccarino 早前的工作给出 $L\le C^d$，其中 $C$ 实质上是 $3$，即 $M\gtrsim3^d s$。**这个 $3^d$ 值得留意**：在张量 Gauss 网格这个受限设计空间里，维数依赖并没有消失，它搬进了常数。因此「线性」这个说法应当理解为「关于 $N$ 线性、但前因子随 $d$ 增长」，而不是无条件优于诱导采样的 $M\gtrsim N\log N$。

### 与其他论文的关系

它是编号 21 的直接前身——同一批作者、同一个「随机子采样 Gauss 求积」的装置，只是用最小二乘而非 $\ell_1$ 极小化。它的 Gauss 权就是 Christoffel 函数值，因此向前接上编号 22 的显式 Christoffel 加权框架与编号 45 的诱导采样理论。它被编号 21、28 以及综述编号 14 的参考文献列表引用。

## 14：把三种重构方式放进同一个账本

### 直觉

到 2015 年，基于配点的重构已有三大家族在广泛使用——最小二乘回归、压缩采样、插值——每一家有自己的理论、自己偏好的点集，而对「什么样的**几何上非结构**的多元网格算好网格」没有统一说法。这篇是这个研究纲领的综述与综合：**把三种重构模式下的稳定性与精度结果收拢到一套记号里**，再把它们当作生成多维随机配点网格的设计指南。因为是带新数值对比的综述，多数定理都归属于前人工作，贡献在于组织方式加上正面交锋的实验。

它给出的深层答案是：三种模式下的要求形状相同——**点集的经验测度必须渐近趋于由区域本身决定的某个测度（多位势平衡测度），而这只是必要条件，不是充分条件**。这句话解释了为什么 Chebyshev 式的边界聚集在完全不同的方法里反复出现。

### 问题设定

模型 $\mathcal L(u;t,x,\omega)=0$，随机性由 $d$ 维随机向量 $Z(\omega)$ 参数化（通常经 Karhunen–Loève 截断得到），响应用 gPC 展开。

**Weil 点**（式 (16)–(17)）：编号 9／10 的确定性 Weil 和网格，分析引擎是 $\bigl|\sum_j e^{2\pi if(j)/M}\bigr|\le(d-1)\sqrt M$。

**结构化随机点**（3.3 节）：先取一个高基数的结构化候选集（例如张量积 Chebyshev 网格，或把双曲叉空间 $\mathcal H^d_k$ 视为张量空间的子空间，$u(z)=\sum_{\alpha\in\Lambda^P_{d,k}}\hat c_\alpha\varphi_\alpha(z)$ 且对 $\alpha\in\Lambda^P_{d,k}\setminus\Lambda^H_{d,k}$ 有 $\hat c_\alpha=0$），再从中随机抽子集，得到实质上非结构的网格。**这正是编号 13 与 21 使用的设计。**

**最小正交插值**（5.1 节）：引入「least-$\rho$」运算 $p_{\downarrow,\rho}=P_{\hat k}p$，$\hat k=\min\{k\in\mathbb N:P_k\ne0\}$，即第一个不消失的「Taylor」贡献，它依赖 $\rho$。对节点 $Z=\{z_1,\dots,z_M\}$ 定义

$$
\Pi_Z=\mathrm{span}\bigl\{g_{\downarrow,\rho}\ :\ g\in\mathrm{span}\{\delta_{z_1},\dots,\delta_{z_M}\}\bigr\} ,
$$

这就是用于插值的最小正交多项式空间。实际计算用 LU 与 QR 的组合，得到分解 $PA=LUH$（式 (39)）：$L$、$U$ 是标准的 $M\times M$ 三角因子，$P$ 是行置换，$H$ 是长方形的；运算量渐近上与标准插值矩阵分解相当。当 $\rho$ 取标准高斯密度时，$\Pi_Z$ 与 de Boor–Ron 的经典「least interpolant」空间重合。

### 推导：把所有要求换成同一种单位

这篇没有一条贯穿的推导，它的动作是把彼此不可比的结论换算成同一个单位——$M$ 对 $N$（稀疏情形对 $s$）——然后排在一起：

| 重构方式与采样                    | 样本要求                                | 来源                  |
| --------------------------------- | --------------------------------------- | --------------------- |
| 最小二乘，均匀（正交）测度 i.i.d. | $M/\log M\ge C_rN^2$                    | 定理 1（CDL）         |
| 最小二乘，Chebyshev Monte Carlo   | $M\sim N^{\log3/\log2}$                 | 该文对照说明          |
| 最小二乘，Hermite 函数 + 映射点   | $M/\log M\gtrsim rN$，$L\gtrsim\sqrt N$ | 定理 2（自编号 11）   |
| 最小二乘，Weil 确定性点           | $M\ge C(d)N^2$                          | 定理 5、6（自编号 9） |
| $\ell_1$，Chebyshev 采样 + 预条件 | $M>C\delta^{-2}L^2s\log^3(s)\log(N)$    | 定理 8、9             |
| 实践中的做法                      | $M\simeq cN$，$c$ 取 2 到 3             | 该文对实践的陈述      |

最后一行是这张表真正的看点：**实践跑在线性区间，而这一栏的理论按该文自己的说法「尚无定论」。** 整篇的组织正是围绕这个缺口展开的。

### 定理

- **定理 1（Cohen–Davenport–Leviatan）。** 假设 $\rho$ 是 $[-1,1]$ 上的均匀分布，$\Lambda=\Lambda^T_{1,N-1}$。若对某 $r>0$ 有 $\frac{M}{\log M}\ge C_rN^2$（$C$ 为普适常数），则 $\Pr\bigl[|\!|\!|\hat A-I|\!|\!|\ge\frac12\bigr]\le2M^{-r}$。**这就是从均匀（正交）测度 i.i.d. 采样的二次要求 $M\gtrsim N^2\log N$。** 该文指出结论可推广到多维、任意**下集**（lower set）指标集，标度不变。
- **对照。** 从 **Chebyshev 测度**做 Monte Carlo 时要求降为 $M\sim N^{\log3/\log2}\approx N^{1.585}$，严格好于 $N^2$ 但仍是超线性。
- **定理 2（取自编号 11）。** 若 $\frac{M}{\log M}\gtrsim rN$ 且 $L\gtrsim\sqrt N$，则用 Hermite **函数**与映射均匀点建立的最小二乘设计矩阵满足 $\Pr\bigl[|\!|\!|\hat A-I|\!|\!|\ge\frac58\bigr]\le2M^{-r}$。该文把它重述为一条**加权**最小二乘结论，因为 Hermite 函数就是加权的 Hermite 多项式——这个重述本身就是全篇统一记号的示范。
- **定理 4（取自编号 9）。** 第 $K$ 个素数生成的 Weil 点 $W_{M_K}$ 渐近按 Chebyshev 测度分布，$\nu_K\to\nu_c$ 弱 * 收敛。
- **定理 5、6（取自编号 9）。** 若 $M\ge C(d)N^2$，则 $\bigl|\!\bigl|\!\bigl|\frac{2^{d+1}}{M}\hat A-I\bigr|\!\bigr|\!\bigr|\le\frac12$，从而最小二乘解唯一；且 $\|f-P^N_Mf\|_{L^2_\rho}\le C\|f-p^*_N\|_{L^\infty}$，$p^*_N$ 为 $L^\infty$ 最佳多项式。该文明说这个二次要求**强于** Chebyshev Monte Carlo 所需的 $M\sim N^{\log3/\log2}$，补偿是确定性。
- **定理 8（Rauhut，有界正交系的 RIP）。** 设 $\sup_n\|\varphi_n\|_\infty\le L$（$L\ge1$），$A$ 为插值矩阵，$W$ 为对角矩阵且 $w_{m,m}=(\pi/2)^{1/2}(1-z_m^2)^{1/4}$，$\{z_m\}$ 从一维 Chebyshev 测度 i.i.d. 抽取。若 $M>C\delta^{-2}L^2s\log^3(s)\log(N)$，则以不小于 $1-N^{-\gamma\log^3(s)}$ 的概率，$\frac{1}{\sqrt M}WA$ 的限制等距常数满足 $\delta_s\le\delta$。
- **定理 9（Rauhut–Ward）。** 同上采样，$A$ 为 Legendre 设计矩阵，$W$ 对角且 $w_{m,m}=(\pi/2)^{-1/2}(1-z_m^2)^{1/4}$。若 $M>Cs\log^3(s)\log(N)$，则 $\min\|c\|_1$ s.t. $WAc=WA\tilde c$ 的解满足 $\Pr\bigl[\|c^{\#}-\tilde c\|_2\le C\sigma_{s,1}(\tilde c)/\sqrt s\bigr]\ge1-N^{-\gamma\log^3(s)}$。该文解释了 $(1-z^2)^{1/4}$ 这个权的作用：**它正好让加权 Legendre 多项式成为一致有界系**。定理 10（Yan–Guo–Xiu）把它推到高维，定理 11（取自编号 10）给出 Weil 点的确定性采样对应物，走的是非相干参数。
- **定理 12（最小正交插值，取自 Narayan–Xiu）。** $\Pi_Z$ 的维数恰为 $M$；对任意连续 $u$，存在唯一的 $p\in\Pi_Z$ 使 $u_m=p(z_m)$；且存在 Lagrange 函数 $\ell_m\in\Pi_Z$ 使 $p(z)=\sum_{m=1}^{M}u_m\ell_m(z)$，$\ell_m(z_n)=\delta_{m,n}$。
- **定理 13（Bloom／Bos 等／Berman–Boucksom，非加权三分律）。** 对点阵列 $Z_{N_k}$ 考虑三个性质：(1) Lebesgue 常数次指数增长，$\lim_k(\Lambda(Z_{N_k}))^{1/k}=1$；(2) 渐近 Fekete，$\lim_k|\det A(Z_{N_k})|^{1/s^d_k}=\delta(D)$；(3) 按多位势平衡测度分布，$\lim_k\frac{1}{N_k}\sum_n\delta_{z_{n,N_k}}=\mu_D$。则 $1\Rightarrow2$、$2\Rightarrow3$，而**逆命题都不成立**。该文抽出的设计教训是：要得到稳定的插值算子，就**必须**渐近地按 $\mu_D$ 采样——必要而不充分。在张量积区间上 $\mu_D$ 就是张量积反正弦测度，这解释了 Chebyshev 式聚集为何反复出现。该文还指出：超立方体上 $d$ 维 Weil 点确实按 $\mu_D$ 分布（定理 4），但它们是否渐近 Fekete**未知**，而且其基数受限，用于插值的价值有限。
- **定理 14（加权版本，取自 Berman–Boucksom 与 Narayan–Xiu）。** 在 $D=\mathbb R^d$ 上取权 $\rho(z)=\exp(\|z\|^r)$、$r\ge1$，对数权 $Q(z)\triangleq-\log\rho=\|z\|^r$，加权平衡测度为 $\mu_{D,Q}$（尽管 $D$ 无界，它的支集是紧的）。考虑：(1) 加权 Lebesgue 常数次指数增长；(2) 收缩后渐近加权 Fekete，$\lim_k\bigl|\det A(k^{-1/r}Z_{N_k})\prod_n\rho^k(z_n)\bigr|^{1/s^d_k}=\delta_\rho(D)$；(2a) 一致收缩有界，即存在紧集 $S\supset\mathrm{supp}\,\mu_{\Omega,Q}$ 使对一切 $k$ 有 $k^{-1/r}Z_{N_k}\subset S$；(3) $\lim_k\frac{1}{N_k}\sum_n\delta_{k^{-1/r}z_{n,N_k}}=\mu_{D,Q}$。则 $1\Rightarrow2$、$(2+2a)\Rightarrow3$，逆命题不成立。**要害是收缩因子 $k^{-1/r}$。** 该文点出了其中的反直觉之处：**不应该**直接从 $\mu_{D,Q}$ 采样，因为它支集紧，对无界域上的多项式逼近用处有限；要的是那些**经 $k^{-1/r}$ 收缩后**按 $\mu_{D,Q}$ 分布的网格。它同时记下一个**公开问题**：多维加权近似／离散 Fekete、Leja 阵列是否按加权多位势平衡测度分布，尚不清楚，尽管一维情形的肯定答案给了希望。

> [!warning] 定理 14 的权与对数权符号不自洽
> 本站记录到的形式是 $\rho(z)=\exp(\|z\|^r)$ 与 $Q(z)\triangleq-\log\rho=\|z\|^r$，这两式不能同时成立（由前者应得 $-\log\rho=-\|z\|^r$）。二者之一有一个符号误植，按该领域的标准写法应为 $\rho=\exp(-\|z\|^r)$、$Q=\|z\|^r$，但本站无法判定误植出在原文还是在转录，故照录并标出。

### 数值实验

| 节  | 重构方式 | 点集对照                                 | 采样规则                   |
| --- | -------- | ---------------------------------------- | -------------------------- |
| 3.4 | 最小二乘 | Gauss 网格子采样、i.i.d. 随机点、Weil 点 | $M=2.5N$ 与 $M=1.5N\log N$ |
| 4.3 | 压缩采样 | Monte Carlo 与确定性采样策略             | —                          |
| 5   | 插值     | 最小正交插值                             | —                          |

3.4 节把设计矩阵条件数对多项式阶数 $k$ 作图；图 2 展示由素数种子 $M=359$（179 个点）与 $M=751$（375 个点）生成的二维 Weil 网格。

这里最该指出的一点是实验规则的选取：$M=2.5N$ 与 $M=1.5N\log N$ 分别是**线性**与**对数线性**，也就是说，**该文的实验刻意跑在它自己的定理 1、5、6 所要求的二次区间之下**。这不是疏忽，正是那张表最后一行的实证形态：实践在线性区间里就能工作，而覆盖该区间的理论当时还没有。本站不为这些图补任何数字。

### 与其他论文的关系

这是整份清单的枢纽。它重述并统一了编号 9（定理 4、5、6）、编号 10（定理 11）、编号 11（定理 2）以及 Rauhut、Rauhut–Ward、Yan–Guo–Xiu 等外部结果，并引入了这一群体随后据以行动的位势论词汇——Lebesgue 常数、Fekete 点、平衡测度、收缩因子。编号 28 用列主元 QR 构造**加权近似 Fekete 点**；编号 22 与 45 走 Christoffel 加权与诱导采样两条路，去够本文点名「可欲而不可得」的线性区间。它关于结构化随机子采样的讨论，正是编号 13 与 21 采用的设计。

## 五篇的样本复杂度对照

| 编号 | 点集或密度                             | 样本要求                                                                 | 结论类型              |
| ---- | -------------------------------------- | ------------------------------------------------------------------------ | --------------------- |
| 6    | SG／MC／QMC 三族候选点                 | 未给出 $M$ 与 $N$ 的关系                                                 | 经验                  |
| 9    | Weil 和确定性点集                      | $M\ge4^{d+1}d^2N^2$，$M$ 为素数                                          | 确定性，无概率限定词  |
| 11   | 映射均匀密度 + Hermite/Laguerre 函数基 | $K\le\kappa m/\log m$，即 $m\gtrsim K\log K$；$L>\max\{3\tau,5\sqrt K\}$ | 概率，$\Pr\le2m^{-r}$ |
| 13   | 张量 Gauss 网格的随机子集              | 线性（仅第三方表述；本文未核实）                                         | 无法核实              |
| 14   | 汇总以上并加入 $\ell_1$ 与插值         | $N^2\log N$／$N^{1.585}$／$N\log N$／$s\log^3 s\log N$                   | 综述与统一            |

一条贯穿的判断：**「取哪些点」这个问题的正确形式是「取多少个点，从什么密度取」。** 编号 6 问的是前者，编号 9 起把它换成后者，而后者有定理可证。更进一步，这五篇合起来指出后者还能再化简一次——从什么密度取，取决于**如何把设计矩阵的行范数拉平**，而拉平它的那个量就是 Christoffel 函数。编号 13 的 Gauss 权与编号 11 的衰减基，是同一个要求在两种设定下的两副面孔。

> [!note] 覆盖进度
> 编号 9、11、14 已按原文全文核对：本页给出的问题设定、推导链条、定理假设与常数、以及实验配置均为转录，其中论文未提供数值的图表只报告配置而不报告数字。编号 6 只到摘要与元数据层面，因此本页不给它定理或实验数值。编号 13 无法核实：出版社封锁正文，researchr 明确记录「摘要缺失」，因此其构造是按设计完全相同的姊妹篇（编号 21）反推的，而「样本数线性增长」这一说法只有第三方文献（Seshadri–Narayan–Sarkar）的表述支持，并非出自该文本身；**其常数与概率界本站不报告**。
>
> 另有两处提醒已就地标出：这一族论文对 Christoffel 函数的记法互相冲突（本站统一取 $K=\sum_\alpha\varphi_\alpha^2$、Christoffel 函数为 $N/K$、权为 $1/K$，引用原文前须核对该文自己的约定）；以及高斯情形的渐近诱导测度在所有陈述它的论文中都只是**猜想**，不得当作定理引用。

## 本页原文

- Z. Gao and T. Zhou, [_On the choice of design points for least square polynomial approximations with application to uncertainty quantification_](https://doi.org/10.4208/cicp.130813.060214a), Commun. Comput. Phys. 16 (2014), pp. 365-381。
- T. Zhou, A. Narayan, and Z. Xu, [_Multivariate discrete least-squares approximations with a new type of collocation grid_](https://doi.org/10.1137/130950434), SIAM J. Sci. Comput. 36(5) (2014), pp. A2401-A2422。
- T. Tang and T. Zhou, [_On discrete least-squares projection in unbounded domain with random evaluations and its application to parametric uncertainty quantification_](https://doi.org/10.1137/140961894), SIAM J. Sci. Comput. 36(5) (2014), pp. A2272-A2295。
- T. Zhou, A. Narayan, and D. Xiu, [_Weighted discrete least-squares polynomial approximation using randomized quadratures_](https://doi.org/10.1016/j.jcp.2015.06.042), J. Comput. Phys. 298 (2015), pp. 787-800。
- A. Narayan and T. Zhou, [_Stochastic collocation on unstructured multivariate meshes_](https://doi.org/10.4208/cicp.020215.070515a), Commun. Comput. Phys. 18 (2015), pp. 1-36。
