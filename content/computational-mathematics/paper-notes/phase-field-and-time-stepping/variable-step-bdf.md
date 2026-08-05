---
title: 变步长 BDF 与卷积核
description: 编号 48、52、58、67、69、74：把历史项汇成一个可判定正定性的对象
lang: zh
translation: en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf
tags:
  - 论文笔记
  - 相场模型
  - 离散能量
---

> [!note] 本页覆盖
> 编号 **48**（_SIAM J. Numer. Anal._ 58(4), 2020）、**52**（_Sci. China Math._ 64, 2021）、**58**（_CSIAM Trans. Appl. Math._ 3, 2022）、**67**（_J. Comput. Math._ 41, 2023）、**69**（_BIT Numer. Math._ 63:39, 2023）、**74**（_Sci. China Math._ 67, 2024）。
>
> 除编号 69 外均已按预印本全文逐式核对。编号 69 无预印本（已对 Liao 的 34 篇预印本逐一排查确认），出版社与各聚合站均不提供正文，因此本页只报告可从完整摘要与参考文献确认的内容：步长比条件、所用成分（单腿改写、离散梯度结构、两类新的 DOC 核）与三条结果。其显式修正能量、DOC 核的定义、收敛阶与数值实验本站未核实。

## 48：两个步长比常数，两条不同的性质

### 障碍

Allen-Cahn 方程在连续层面同时具有能量耗散与最大值原理，而解跨越两个差别极大的时间尺度，因此实践上必须用变步长。但在本文写作时，即便对线性与半线性抛物方程，变步长格式的数值分析用作者的话说仍「远未完成」。具体的缺口是：没有人对 Allen-Cahn 方程的二阶变步长格式证明过**离散最大值原理**；而已有的非均匀 BDF2 能量稳定性结果所依赖的步长比条件，其尖锐性也不清楚。

经典困难是：变步长 BDF2 不自伴，其核不是卷积（Toeplitz）型，因此常步长下的能量测试函数技巧与正定性论证都不适用。

### 格式与两条条件

变步长 BDF2 公式为

$$
D_{2}v^{n}=\frac{1+2r_{n}}{\tau_{n}(1+r_{n})}\nabla_{\tau}v^{n}
-\frac{r_{n}^{2}}{\tau_{n}(1+r_{n})}\nabla_{\tau}v^{n-1},
\qquad
r_{k}=\frac{\tau_{k}}{\tau_{k-1}},
$$

写成离散卷积 $D_{2}v^{n}=\sum_{k=1}^{n}b^{(n)}_{n-k}\nabla_{\tau}v^{k}$ 时核为

$$
b^{(n)}_0=\frac{1+2r_n}{\tau_n(1+r_n)},
\qquad
b^{(n)}_1=-\frac{r_n^2}{\tau_n(1+r_n)},
\qquad
b^{(n)}_j=0\ (j\ge2),
$$

$b^{(n)}_1<0$ 正是核不非负、需要重排的原因。全离散格式是**完全隐式**的（无稳定化、无凸分裂）：

$$
D_{2}u^{n}=\varepsilon^{2}\Lambda_{h}u^{n}-f(u^{n}),
\qquad f(u^n)=(u^n)^{.3}-u^n .
$$

论文用两个命名条件：

- **S1**：$0<r_{k}<\dfrac{3+\sqrt{17}}{2}\approx3.561$，用于**能量稳定性**。这个常数并非本文首创：它出自 Liao 与 Zhang（_Math. Comp._ 90 (2021) 1207-1226）的 Lemma 2.1，作为 BDF2 核正定性的条件，编号 48 与 52 都是把它作为 S1 引入的。编号 52 对它的评价很直接——称其为「由条件 S1 带来的人为常数」。
- **S0**：$0<r_{k}<1+\sqrt{2}\approx2.414$，用于**离散最大值原理**与最大范数收敛；论文明确指出它与 Grigorieff（1983）对常微分方程的零稳定性条件一致。

### 修正能量与 3.561 的精确来源

耗散的对象不是普通离散能量，而是修正能量

$$
\widehat{E}[u^{k}]=E[u^{k}]+\frac{r_{k+1}\tau_{k}}{2(1+r_{k+1})}
\sum_{i}\bigl(\partial_{\tau}u_{i}^{k}\bigr)^{2},
\qquad
\widehat{E}[u^{0}]=E[u^{0}],
$$

$$
E[u^{k}]=-\frac{\varepsilon^{2}}{2}(u^{k})^{T}\Lambda_{h}u^{k}
+\frac{1}{4}\sum_{i}\bigl(1-(u_{i}^{k})^{2}\bigr)^{2}.
$$

修正项是 $O(\tau)$，故 $\tau\to0$ 时 $\widehat E\to E$。注意它用的是 $r_{k+1}$，即**下一个**步长比；这个前瞻正是让求和电报式成立的原因。

用 $2a(a-b)=a^2-b^2+(a-b)^2$ 与核的定义可得

$$
D_{2}u_{i}^{n}\,(\nabla_\tau u_{i}^{n})
\ \ge\
\frac{r_{n+1}\tau_{n}}{2(1+r_{n+1})}(\partial_{\tau}u_{i}^{n})^{2}
-\frac{r_{n}\tau_{n-1}}{2(1+r_{n})}(\partial_{\tau}u_{i}^{n-1})^{2}
+\Bigl(\frac{2+4r_{n}-r_{n}^{2}}{1+r_{n}}-\frac{r_{n+1}}{1+r_{n+1}}\Bigr)
\frac{\tau_{n}}{2}(\partial_{\tau}u_{i}^{n})^{2}.
$$

前两项电报式抵消进 $\widehat E$，最后一项必须非负（在吸收非线性项贡献之后），这给出条件

$$
\frac{r_{k+1}}{1+r_{k+1}}<\frac{r_{s}}{1+r_{s}}=\frac{\sqrt{17}-1}{4}\approx0.78,
\qquad
r_{s}=\frac{3+\sqrt{17}}{2}\ \text{是}\ 2+3r-r^{2}=0\ \text{的正根}.
$$

论文还给出对应的步长上界讨论：令 $h(x)=\dfrac{2+4x-x^{2}}{1+x}$，其导数为 $h'(x)=\dfrac{x+1+\sqrt3}{(1+x)^{2}}(\sqrt3-1-x)$，因此 $h$ 在 $(0,\sqrt3-1]$ 上递增、之后递减，$h(0)=2$、$h(\sqrt2+1)=1+\frac{\sqrt2}{2}$。据此分三段给出足够的步长上界，第三段（$\sqrt2+1<r_k<r_s$）还需要控制**下一个**步长比。

### $1+\sqrt2$ 的来源：核重排

最大值原理的门槛来自一个不同的技术：**核重排与互补**。引入参数 $\eta$ 与

$$
\bar{v}^{0}=v^{0},\qquad \bar{v}^{k}=v^{k}-\eta\,v^{k-1}\ (k\ge1),
\qquad
v^{k}=\sum_{\ell=0}^{k}\eta^{k-\ell}\bar{v}^{\ell},
$$

代入并交换求和顺序，得到重排后的 BDF2 公式，其核在适当选择 $\eta$ 后变为非负。零稳定性条件 $r_k<1+\sqrt2$ 正是这一重排可行的范围。

**两个常数对应两条不同性质，这一点在引用时必须区分。** 能量稳定性来自二次型的正定性（S1），最大值原理来自核重排后的非负性（S0），二者不能互换。

## 52：同一套分析用于分子束外延模型

编号 52 把变步长 BDF2 的分析用于**无斜率选择**的分子束外延模型，其自由能为

$$
E_m[\phi]=\frac{\varepsilon}{2}\|\Delta\phi\|^{2}
-\frac{1}{2\varepsilon}\bigl\langle\ln|1+|\nabla\phi|^{2}|,1\bigr\rangle .
$$

这个模型的非线性项 $-\nabla\phi/(1+|\nabla\phi|^{2})$ 有界，因此没有最大值原理可用，能量论证成为唯一工具。论文在同一门槛 $r_k<(3+\sqrt{17})/2$ 下给出能量稳定性与 $L^2$ 稳定性及收敛性。

## 58：把不 A-稳定的 BDF-$k$ 拉回教科书式论证

### 问题

$A$-稳定的 BDF1 与 BDF2 允许直接的教科书式离散能量证明（用 $u^n$ 在 $L^2$ 中作检验）。对 $3\le k\le5$，BDF-$k$ **不** $A$-稳定，这条直接论证失效。自 Lubich、Mansour 与 Venkataraman 以来的标准对策是 **Nevanlinna-Odeh 乘子技术**，它基于 Dahlquist 关于 $A$-稳定与 $G$-稳定等价的结论，用 $u^n-\sum_i\eta_iu^{n-i}$ 而不是 $u^n$ 作检验；另一条是 Liu 的伸缩公式。两者都引入人为乘子，并且关键地**对启动数据要求更强的范数**（乘子型稳定性估计中会出现 $H^1$ 型量）。论文提出并回答的问题是：对 $3\le k\le5$ 的 BDF-$k$，是否存在一个直接的离散能量分析？

### DOC 核与正交性

把 BDF-$k$ 写成离散卷积

$$
D_{k}v^n=\frac1{\tau}\sum_{j=1}^n b_{n-j}^{(k)}\,\nabla_\tau v^{j},
\qquad n\ge k,
$$

其核由生成函数

$$
\sum_{\ell=1}^{k}\frac{1}{\ell}(1-\zeta)^{\ell-1}
=\sum_{\ell=0}^{k-1}b_{\ell}^{(k)}\zeta^{\ell}
$$

给出，例如 $k=3$ 时 $(b_0,b_1,b_2)=(11/6,-7/6,1/3)$，$k=5$ 时 $(137/60,-163/60,137/60,-21/20,1/5)$。

**离散正交卷积核** $\theta_j^{(k)}$ 由

$$
\theta_{0}^{(k)}=\frac{1}{b_{0}^{(k)}},
\qquad
\theta_{n-j}^{(k)}=-\frac{1}{b_{0}^{(k)}}\sum_{\ell=j+1}^{n}\theta_{n-\ell}^{(k)}b_{\ell-j}^{(k)}
$$

递归定义，其关键性质是正交性恒等式

$$
\sum_{\ell=j}^{n}\theta_{n-\ell}^{(k)}\,b_{\ell-j}^{(k)}\equiv\delta_{nj},
\qquad k\le j\le n,
$$

以及相互正交性 $\sum_{\ell=j}^{n}b_{n-\ell}^{(k)}\theta_{\ell-j}^{(k)}\equiv\delta_{nj}$。因此

$$
\sum_{j=k}^{n}\theta_{n-j}^{(k)}D_{k}u^{j}
=\frac1{\tau}u_{\mathrm I}^{(k,n)}+\partial_{\tau}u^{n},
\qquad
u_{\mathrm I}^{(k,n)}=\sum_{\ell=1}^{k-1}\nabla_\tau u^{\ell}
\sum_{j=k}^{n}\theta_{n-j}^{(k)}b_{j-\ell}^{(k)} .
$$

也就是说，**用 DOC 核作用把 BDF-$k$ 卷积逆回单个一阶差分**，只差一个启动值余项。由于相互正交性同时成立，这是一个**可逆**的离散变换：没有信息损失。变换后的格式取 $w=2\tau u^j$ 并对 $j$ 求和，直接给出经典能量不等式。

### 正定性由生成函数判定

论文把 $\{b_j^{(k)}\}$ 的正（半）定性与 $\{\theta_j^{(k)}\}$ 的正定性建立等价，再用对称化 Toeplitz 矩阵的生成函数

$$
\mathrm g^{(k)}(\varphi)=2\sum_{j=0}^{k-1}b_j^{(k)}\cos(j\varphi)
$$

与 Grenander-Szegő 定理判定：$\mathrm g^{(k)}_{\min}\le\lambda_{\min}(B_k)\le\lambda_{\max}(B_k)\le\mathrm g^{(k)}_{\max}$。例如

$$
\mathrm g^{(3)}(\varphi)=\tfrac13\bigl(11-7\cos\varphi+2\cos2\varphi\bigr)
=\tfrac43\bigl(\cos\varphi-\tfrac78\bigr)^2+\tfrac{95}{48},
$$

配方后下界显然为正。由此得到对 $3\le k\le5$ 与任意实序列的二次型下界

$$
2\sum_{m=k}^{n}w_m\sum_{j=k}^{m}b_{m-j}^{(k)}w_j
\ \ge\ \sigma_k\sum_{j=k}^{n}w_j^2 ,
$$

带显式常数 $\sigma_k$。

> [!note] 步长设定
> 本文的时间网格是**均匀**的。$k\ge3$ 的变步长情形在结论部分被明确列为公开问题——这正是编号 67 接手的对象。

## 67：变步长 BDF3

### 为什么需要新分析

变步长 BDF3 此前基本只有一条经典结果：Calvo 与 Grigorieff（2002）在步长比条件 $r_k<1.199$ 下证明 $L^2$ 稳定性，估计形如

$$
\|u^n\|\le C\exp(C\Gamma_n)\Bigl(\|u_0\|+\sum_{j=1}^{n}\tau_j\|f^j\|\Bigr),
\qquad \Gamma_n=\sum_{k=2}^{n}|r_k-r_{k-1}| .
$$

**前因子 $\exp(C\Gamma_n)$ 不是网格稳健的。** 取交替步长 $\{\tau_1,\mu\tau_1,\tau_1,\mu\tau_1,\dots\}$（$\mu\ne1$）并固定 $T=\frac M2(1+\mu)\tau_1$，则 $\Gamma_M=(M-1)|\mu-\mu^{-1}|\to\infty$（$\tau_1\to0$），界恰好在它本应覆盖的自适应区域内退化。本文用一个常数与步长比完全无关的分析取代它，代价是一个步长比门槛。

### 变步长 BDF3 与其 DOC 核

变步长 BDF3 写成

$$
D_3v^n=d_0(r_n,r_{n-1})\partial_\tau v^{n}
+d_1(r_n,r_{n-1})\partial_\tau v^{n-1}
+d_2(r_n,r_{n-1})\partial_\tau v^{n-2},
$$

$$
d_0(x,y)=\frac{1+2x}{1+x}+\frac{xy}{1+y+xy},
\qquad
d_2(x,y)=\frac{xy^2}{1+y+xy}\cdot\frac{1+x}{1+y},
$$

$d_1$ 为使公式相容的相应组合。变步长 DOC 核由

$$
\vartheta_0^{(n)}=\frac{1}{d_0^{(n)}},
\qquad
\vartheta_{n-j}^{(n)}=-\frac{1}{d_0^{(j)}}\sum_{i=j+1}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-j}
$$

定义，满足正交性 $\sum_{i=j}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-j}\equiv\delta_{nj}$，等价于矩阵形式 $\Theta_3D_3=I$；由 $D_3\Theta_3=I$ 又得相互正交性。论文在门槛 $r_k<1.4877$ 下给出离散梯度结构、能量耗散律与 $L^2$ 稳定性及收敛性。

### $1.4877$ 是怎么来的

这个门槛不是一个天然常数，而是一次**参数取舍**的结果，值得写清楚。取 $\gamma=7/10$，则 $R_e$ 是

$$
d_1(R_e,0)+\tfrac{7}{10}\sqrt{R_e}\,d_2(R_e,R_e)=0
\qquad\Longleftrightarrow\qquad
\frac{10}{7(R_e+1)}-\frac{R_e^2\sqrt{R_e}}{R_e^2+R_e+1}=0
$$

的**唯一正根**，数值上 $R_e\approx1.4877$。

$\gamma=7/10$ 的来历是这样：离散梯度分解需要两个条件 $q_{n+1}\ge0$ 与 $p_{n+1}>0$ 同时成立，而它们含 $r_{n+1},r_n,r_{n-1},\gamma,R_e$ 五个变量，无法精确求解。论文的做法是在等比网格上处理：由 $q_{n+1}\ge0$ 在 $r_{n-1}=0$、$r_{n+1}=r_n=r$ 时得到 $\gamma\le-d_1(r,0)/(\sqrt r\,d_2(r,r))$；把第二个条件也放到等比网格上，得 $\bar R_e\approx1.4965$ 与 $\bar\gamma\approx0.6924$。作者随后**固定 $\gamma=7/10$**（与 $0.6924$ 很接近）以换取可处理性，由此得到 $R_e\approx1.4877$。这个取法的依据是：$q_{n+1}\ge0$ 是「必要且尖锐」的，而 $p_{n+1}>0$ 可以放松。

因此这个门槛的**尖锐程度是可以估量的**：数值实验给出 $R_e<1.69$ 是必要的，而理论给出 $R_e<1.4877$ 是充分的，两者之间的间隙很小。这与编号 52 对 $3.561$ 的自我评价（「由条件 S1 带来的人为常数」）恰好构成对照——同样是人为选定的常数，这里的距离是被量化了的。

## 69 与 74：把工具本身作为研究对象

**编号 69** 分析变步长时间滤波后向 Euler 格式的稳定性与收敛性。时间滤波是在低阶格式上加一个后处理以提高精度的技术，而变步长下它的稳定性同样落回这一套卷积核分析。论文的能量稳定性与 $L^2$ 误差估计在步长比条件

$$
\tfrac12\le r_k\le2
$$

下成立，而这个条件与本页其余门槛**不是同一类对象**，有两点必须分开。

第一，它是**双边**的。本页其他所有结果都只限制步长比的**上界**（$1+\sqrt2$、$3.561$、$1.4877$），而这里步长比还被 $1/2$ 从下方限制——也就是说步长既不能放大太快，也不能**收缩**太快。这反映的是变步长滤波后向 Euler 失去了 A-稳定性。

第二，$[1/2,2]$ **不是一个尖锐的解析门槛**，而是自适应程序里常用的经验保护区间，论文本身称之为「实用」约束。因此不应把它与 $1+\sqrt2$、$3.561$、$1.4877$ 并列引用，好像它们是同一种量。

**编号 74** 直接研究**工具本身**：变步长 L1 型卷积算子逼近所得实二次型的正定性。这一篇提供的是代数判据，被前面几篇反复调用——例如[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 57]] 在证明其 DOC 核单调性时就引用了本文给出的关于核序列的三条不等式。

**把编号 74 单列出来的意义在于**：这一系列工作的技术核心不是某个格式，而是「历史项汇成的实二次型何时正定」这一代数问题。一旦这个问题有了独立的判据，同一套论证就可以搬到三阶 BDF、时间滤波 Euler、分数阶 L1 逼近与隐显 Runge-Kutta 上。

## 六篇的关系

| 编号 | 对象                  | 网格   | 门槛                      | 核心工具               |
| ---- | --------------------- | ------ | ------------------------- | ---------------------- |
| 48   | Allen-Cahn，BDF2      | 变步长 | S1 $3.561$；S0 $1+\sqrt2$ | 修正能量；核重排       |
| 52   | 分子束外延，BDF2      | 变步长 | $3.561$                   | 修正能量               |
| 58   | 线性反应扩散，BDF-$k$ | 均匀   | 无（$3\le k\le5$）        | DOC 核；生成函数正定性 |
| 67   | 扩散方程，BDF3        | 变步长 | $1.4877$                  | 变步长 DOC 核          |
| 69   | 抛物方程，滤波 Euler  | 变步长 | $[1/2,2]$（双边，经验）   | 同族卷积核分析         |
| 74   | 二次型本身            | 变步长 | 给出判据                  | 代数判据               |

## 覆盖核对

| 内容                            | 论文 | 覆盖状态                                   |
| ------------------------------- | ---- | ------------------------------------------ |
| 变步长 BDF2 公式与核的符号      | 48   | 公式、卷积核、$b^{(n)}_1<0$ 的含义         |
| 条件 S0 与 S1 及各自的性质      | 48   | 两个常数、两条性质、与 Grigorieff 的一致性 |
| 修正能量与前瞻项                | 48   | 形式、$O(\tau)$、为何用 $r_{k+1}$          |
| 3.561 的精确推导与步长上界讨论  | 48   | 下界不等式、正根、$h(x)$ 的单调性分段      |
| 核重排与 $1+\sqrt2$             | 48   | $\eta$ 变换、替换公式、门槛来源            |
| 分子束外延自由能与无最大值原理  | 52   | 能量形式与工具的限制                       |
| BDF-$k$ 的 DOC 核与可逆变换     | 58   | 生成函数、核表、正交性、逆回一阶差分       |
| 生成函数与 Grenander-Szegő 判定 | 58   | $\mathrm g^{(k)}$、配方、二次型下界        |
| Calvo-Grigorieff 前因子的不稳健 | 67   | 交替网格反例与 $\Gamma_M$ 的发散           |
| 变步长 BDF3 与其 DOC 核         | 67   | $d_0,d_2$、DOC 递归、双向正交性、门槛      |
| 二次型正定性作为独立对象        | 74   | 其在整条线索中的地位                       |

## 本页原文

- H.-l. Liao, T. Tang, and T. Zhou, [_On energy stable, maximum-principle preserving, second-order BDF scheme with variable steps for the Allen-Cahn equation_](https://doi.org/10.1137/19M1289157), SIAM J. Numer. Anal. 58(4) (2020), pp. 2294-2314（预印本 [arXiv:2003.00421](https://arxiv.org/abs/2003.00421)）。
- H.-l. Liao, X. Song, T. Tang, and T. Zhou, [_Analysis of the second-order BDF scheme with variable steps for the molecular beam epitaxial model without slope selection_](https://doi.org/10.1007/s11425-020-1817-4), Sci. China Math. 64 (2021), pp. 887-902。
- H.-l. Liao, T. Tang, and T. Zhou, [_A new discrete energy technique for multi-step backward difference formulas_](https://doi.org/10.4208/csiam-am.SO-2021-0032), CSIAM Trans. Appl. Math. 3 (2022), pp. 318-334（预印本 [arXiv:2102.04644](https://arxiv.org/abs/2102.04644)）。
- H.-l. Liao, T. Tang, and T. Zhou, [_Discrete energy analysis of the third-order variable-step BDF time-stepping for diffusion equations_](https://doi.org/10.4208/jcm.2207-m2022-0020), J. Comput. Math. 41 (2023), pp. 325-344（预印本 [arXiv:2204.12742](https://arxiv.org/abs/2204.12742)）。
- H.-l. Liao, T. Tang, and T. Zhou, [_Stability and convergence of the variable-step time filtered backward Euler scheme for parabolic equations_](https://doi.org/10.1007/s10543-023-00982-y), BIT Numer. Math. 63 (2023), 39。
- H.-l. Liao, T. Tang, and T. Zhou, [_Positive definiteness of real quadratic forms resulting from the variable-step L1-type approximations of convolution operators_](https://doi.org/10.1007/s11425-022-2229-5), Sci. China Math. 67 (2024), pp. 237-252。
