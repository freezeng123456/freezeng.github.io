---
title: 4.3–4.4：PFASST 与 MGRiT
description: 从配置方程、SDC 平滑和层间插值到 FCF 重叠 Parareal 的完整结构与成本比较
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-2-pfasst-mgrit
tags:
  - 时间并行
  - PFASST
  - MGRiT
---

> [!note] 阅读范围
> 本页对应论文 Sections 4.3–4.4（pp. 452–461），覆盖公式 (4.10)–(4.13)、Theorems 4.5–4.6 和 Figures 4.6–4.11。PFASST 的配置方程、层间传递、SDC 近似以及 MGRiT 的重叠解释均展开到可实现的矩阵形式。

## 4.3 PFASST

### 来源与两层节点

PFASST 由 Emmett 与 Minion（2012）提出。早期思路来自 Minion
（2010）：用一次 SDC 迭代（Dutt et al. 2000）替代 Parareal 中
昂贵的完整细传播。Minion 等（2015）给出代数表述，Bolten、Moser 与
Speck（2017, 2018）进一步把它解释为时间多重网格：配置方程是细层，
低阶 SDC 是平滑器。由于配置节点、SDC 平滑和跨时间步层间传递同时
出现，下面先固定记号，再拆解一次完整迭代。

把 $(0,T)$ 分成 $N_t$ 个区间 $[T_n,T_{n+1}]$。本节重新用
$\Delta t=T_{n+1}-T_n$ 表示一个 PFASST 区间的长度；它不是前页
Parareal 里的细步长。每个区间内设置 $M_f$ 个细配置未知点和
$M_c$ 个粗配置未知点：

$$
t_{n,m}^{f}=T_n+\tau_m^f\Delta t,
\qquad
t_{n,m}^{c}=T_n+\tau_m^c\Delta t,
$$

其中 $\tau_0^{f,c}=0$ 是另列的初值位置，
$\tau_{M_{f,c}}^{f,c}=1$，且 $M_f>M_c$。因此 $M_f=3$ 表示
三个配置未知点，再加一个 $\tau_0^f$，本地时间位置共四个。

### 配置方程 (4.10)

对 $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$，在任一层的 $M$ 个节点上写成

$$
\boldsymbol u_{n,m}
=\boldsymbol u_{n,0}
+\Delta t\sum_{j=1}^{M}q_{m,j}
[A\boldsymbol u_{n,j}+\boldsymbol g(t_{n,j})],
\qquad m=1,\ldots,M. \tag{4.10}
$$

令

$$
Q=(q_{m,j}),\qquad
\boldsymbol u_n=(\boldsymbol u_{n,1}^\top,\ldots,\boldsymbol u_{n,M}^\top)^\top,
$$

$$
\boldsymbol g_n=
\left(
\boldsymbol g(t_{n,1})^\top,\ldots,
\boldsymbol g(t_{n,M})^\top
\right)^\top,
\qquad
\boldsymbol b_n=(Q\otimes I_x)\boldsymbol g_n,
$$

$$
\chi=
\begin{bmatrix}
0&\cdots&0&1\\
\vdots&&\vdots&\vdots\\
0&\cdots&0&1
\end{bmatrix},
\qquad \boldsymbol\chi=\chi\otimes I_x,
$$

则末节点被复制为下一区间所有配置节点的初值：

$$
\boldsymbol u_n
=\Delta t(Q\otimes A)\boldsymbol u_n
+\boldsymbol\chi\boldsymbol u_{n-1}
+\Delta t\boldsymbol b_n.
$$

细、粗层分别记为

$$
\Phi_f=I_f-\Delta tQ_f\otimes A,
\qquad
\Phi_c=I_c-\Delta tQ_c\otimes A,
$$

其中

$$
I_f=I_{M_f}\otimes I_x,
\qquad I_c=I_{M_c}\otimes I_x.
$$

$$
\boldsymbol u_n^f
=\Phi_f^{-1}(\boldsymbol\chi_f\boldsymbol u_{n-1}^f+\Delta t\boldsymbol b_n^f),
\qquad
\boldsymbol u_n^c
=\Phi_c^{-1}(\boldsymbol\chi_c\boldsymbol u_{n-1}^c+\Delta t\boldsymbol b_n^c).
$$

### Lagrange 层间传递

粗层节点上的数据定义插值多项式

$$
p^c(\tau;\boldsymbol u^c)
=\sum_{m=1}^{M_c}u_m^cL_m^c(\tau),
\qquad
L_m^c(\tau)=
\prod_{\substack{j=1\\j\ne m}}^{M_c}
\frac{\tau-\tau_j^c}{\tau_m^c-\tau_j^c}.
$$

在细节点 $\{\tau_m^f\}$ 取值便得到 $T_{c\to f}$。反向以细层 Lagrange 基函数在粗节点取值，得到 $T_{f\to c}$。两个矩阵都与 $I_x$ 做 Kronecker 积，因此只在配置节点方向插值，不改变空间变量。

> [!warning] 原文公式核对：Lagrange 基函数
> 期刊版与 arXiv 版把分子中的节点下标和分母方向都排错了；其后给出的
> 数值传递矩阵对应的正是上面的标准基函数，而不是印刷式。本页保留与
> 数值矩阵一致的修正式。

在 Gander et al. (2023b) 的块迭代表述中，PFASST 为

$$
\boldsymbol u_{n+1}^{k+1}
=\mathbf B_1^0\boldsymbol u_{n+1}^k
+\mathbf B_0^1(\boldsymbol\chi\boldsymbol u_n^{k+1}+\Delta t\boldsymbol b_n^f)
+\mathbf B_0^0(\boldsymbol\chi\boldsymbol u_n^k+\Delta t\boldsymbol b_n^f),
$$

$$
\begin{aligned}
\mathbf B_1^0&=[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]
(I_f-\widetilde\Phi_f^{-1}\Phi_f),\\
\mathbf B_0^1&=T_{c\to f}\Phi_c^{-1}T_{f\to c},\\
\mathbf B_0^0&=[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]\widetilde\Phi_f^{-1}.
\end{aligned}
$$

下标记时间偏移，上标记迭代偏移，与原文记号一致。

三个块分别描述细层 SDC 平滑、粗层新迭代传播和旧迭代校正。$\widetilde\Phi_f$ 是易解的细层配置矩阵近似。

### SDC 近似 (4.11) 与 Figure 4.6

用相邻细节点之间的隐式 Euler 构造 $\widetilde\Phi_f$：

$$
\frac{\boldsymbol u_{n,m+1}-\boldsymbol u_{n,m}}
{\Delta t(\tau_{m+1}^f-\tau_m^f)}
=A\boldsymbol u_{n,m+1}+\boldsymbol g(t_{n,m+1}^f),
\quad m=0,\ldots,M_f-1. \tag{4.11}
$$

它的矩阵是下双对角差分矩阵与 $I_x$ 的 Kronecker 积，再减去由节点差 $\tau^f_{m+1}-\tau^f_m$ 构成的对角矩阵与 $\Delta tA$ 的 Kronecker 积：

$$
\widetilde\Phi_f=
\begin{bmatrix}
1\\-1&1\\&\ddots&\ddots\\&&-1&1
\end{bmatrix}\otimes I_x
-\Delta t
\begin{bmatrix}
\tau_1^f-\tau_0^f\\
&\ddots\\
&&\tau_{M_f}^f-\tau_{M_f-1}^f
\end{bmatrix}\otimes A.
$$

一次求解对应一轮低阶 SDC sweep，成本低于完整配置求解。

数值实验取 $T=3$、周期边界、零初值，源项用 (2.4) 且 $\sigma=1000$，$\Delta x=1/128$、$\Delta t=1/64$。细层采用三节点 Radau IIA：

$$
\left\{0,\frac{4-\sqrt6}{10},\frac{4+\sqrt6}{10},1\right\},
$$

即 $M_f=3$；粗层采用 $\{0,1/3,1\}$，即 $M_c=2$，两层都使用
Radau IIA。对应的权重矩阵是

$$
Q_f=
\begin{bmatrix}
\frac{88-7\sqrt6}{360}
&\frac{296-169\sqrt6}{1800}
&\frac{-2+3\sqrt6}{225}\\
\frac{296+169\sqrt6}{1800}
&\frac{88+7\sqrt6}{360}
&\frac{-2-3\sqrt6}{225}\\
\frac{16-\sqrt6}{36}
&\frac{16+\sqrt6}{36}
&\frac19
\end{bmatrix},
\qquad
Q_c=
\begin{bmatrix}
\frac5{12}&-\frac1{12}\\
\frac34&\frac14
\end{bmatrix}.
$$

由同一组节点得到数值传递矩阵

$$
T_{c\to f}=
\begin{bmatrix}1.2674&-0.2674\\0.5325&0.4674\\0&1\end{bmatrix}
\otimes I_x,
\qquad
T_{f\to c}=
\begin{bmatrix}0.5018&0.6833&-0.1851\\0&0&1\end{bmatrix}
\otimes I_x;
$$

这些矩阵完全由节点和 Lagrange 插值确定。

![原论文 Figure 4.6：热方程与三种黏性对流扩散方程上的 PFASST 误差](assets/papers/time-parallelization/source-figures/figure-4-6.svg)

三条对流扩散曲线取 $\nu=0.1,10^{-3},10^{-4}$，横轴迭代指标到
$300$。热方程与 $\nu=0.1$ 在 $k\approx90$–$95$ 时降到
$10^{-12}$，$\nu=10^{-3}$ 约需 $120$ 轮，$\nu=10^{-4}$ 约需
$290$ 轮。黏性减弱后粗传播子不再足以表示主导对流，迭代数随之增加。

> [!note] 本站补充：与 Parareal 的对照
> 原文已把弱扩散下的退化与 Parareal、MGRiT 作现象层面的对照。
> 进一步用第 4.2 节的语言解释，则是持续传播的高频越来越难由粗配置层
> 表示；这一机制化解读属于本站补充。

## 4.4 MGRiT

### 作为带重叠的 Parareal

MGRiT 由 Falgout 等（2014）提出，也可解释为代数多重网格、块迭代或 Parareal 的重叠版本。两层 FCF 松弛对非线性系统的更新为

$$
\boldsymbol u_0^{k+1}=\boldsymbol u_0,
\qquad
\boldsymbol u_1^{k+1}=\mathcal F(T_0,T_1,\boldsymbol u_0),
$$

$$
\begin{aligned}
\boldsymbol u_{n+1}^{k+1}
={}&\mathcal F\!\left(T_n,T_{n+1},
\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)\right)\\
&+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})\\
&-\mathcal G\!\left(T_n,T_{n+1},
\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)\right),
\quad n=1,\ldots,N_t-1. \tag{4.12}
\end{aligned}
$$

每轮有两个细传播，Parareal 每轮只有一个。

![原论文 Figure 4.7：FCF-MGRiT 是重叠宽度为一个粗时间步的 Parareal](assets/papers/time-parallelization/source-figures/figure-4-7.svg)

Figure 4.7 中深色圆是粗点。额外的 CF 松弛把旧迭代信息先向前推进一个粗区间，相当于重叠宽度 $\Delta T$。因此全局误差至多 $\lceil N_t/2\rceil$ 轮归零。一般 $F(CF)^\nu$ 对应 $\nu\Delta T$ 的重叠。这些结果见 Gander 等（2018b）的 Theorem 5、Theorem 6 和 Corollary 1，其中 Theorem 6 给出的是非线性问题的超线性收敛结论。

### Theorem 4.5：长时间模态因子

Theorem 4.5 引自 Dobrev 等（2017）。沿用 Theorem 4.2 的记号：
$\mathcal F,\mathcal G$ 都是单步积分器，
$\boldsymbol u_n$ 是顺序细传播得到的目标解，
$\boldsymbol e_n^k=V_A(\boldsymbol u_n^k-\boldsymbol u_n)$，
$A$ 可对角化且 $\sigma(A)\subset\mathbb C_-$。若
$|R_g(z)|<1$，两层 FCF-MGRiT 满足

$$
\max_n\|\boldsymbol e_n^k\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}\varrho_l^k(J,z)
\max_n\|\boldsymbol e_n^0\|_\infty,
$$

$$
\varrho_l(J,z)=
\frac{|R_f^J(z/J)|\,|R_g(z)-R_f^J(z/J)|}
{1-|R_g(z)|}. \tag{4.13}
$$

因此

$$
\varrho_{l,\mathrm{MGRiT}}
=|R_f^J(z/J)|\,\varrho_{l,\mathrm{Parareal}}.
$$

若 $|R_f(z/J)|\le1$（等价地，
$|R_f^J(z/J)|\le1$），额外的 CF 松弛提供一次细传播的收缩，同时也
增加一次昂贵的并行细求解。细积分器在整个左半平面 A-稳定时，这个
条件自动成立。

![原论文 Figure 4.8：等细求解成本下 MGRiT 与两轮 Parareal 的复平面收敛域](assets/papers/time-parallelization/source-figures/figure-4-8.svg)

Figure 4.8 用后向 Euler 粗传播 $R_g(z)=1/(1-z)$ 与精确细传播 $R_f(z)=e^z$，比较一轮 FCF-MGRiT 和两轮 Parareal。上排 (a) 是 MGRiT，下排 (b) 是按相同两次细求解成本绘制的 $\varrho_{l,\mathrm{Parareal}}^2$；三列从左到右对应 $\widehat\varrho=0.2,0.4,0.6$。每一列的着色区域轮廓都很接近，说明公平比较的单位应是细传播次数。SDIRK 细传播也给出相同结论。

### Theorem 4.6：等成本常数比较

Theorem 4.6 引自 Wu 与 Zhou（2019）。若细传播 L-稳定且
$J=\Delta T/\Delta t=O(1)$，令 $s=-z\ge0$。后向 Euler 粗传播下
$\max_{s\ge0}\varrho_l(J,-s)$ 约为

$$
\max\varrho_l\approx
\begin{cases}
0.2984,&\text{Parareal},\\
0.1115,&\text{FCF-MGRiT}.
\end{cases}
$$

二阶 Lobatto IIIC 粗传播下，同一个
$\max_{s\ge0}\varrho_l(J,-s)$ 为

$$
\max\varrho_l\approx
\begin{cases}
0.0817,&\text{Parareal},\\
0.0197,&\text{FCF-MGRiT}.
\end{cases}
$$

按两次细求解配平，一轮 MGRiT 稍慢于两轮 Parareal：$0.2984^2=0.0890<0.1115$，$0.0817^2=0.0067<0.0197$。FCF 的单轮曲线更陡，单位细求解成本的优势并不自动成立。

### Figures 4.9–4.11：线性与非线性实验

线性实验取热方程齐次 Dirichlet 边界、ADE 周期边界，共用

$$
u_0(x)=\sin^2(8\pi(1-x)^2),\quad
T=5,\quad J=20,\quad \Delta T=1/8,\quad \Delta x=1/160,
$$

粗层后向 Euler，细层 SDIRK22。

![原论文 Figure 4.9：热方程与两种黏性 ADE 上的模态因子分布](assets/papers/time-parallelization/source-figures/figure-4-9.svg)

Figure 4.9 先比较逐模态因子：上排是 MGRiT，下排是 Parareal，
三列依次对应热方程、$\nu=0.1$ 和 $0.01$ 的 ADE。最大因子分别为
MGRiT 的 $0.08375,0.2718,0.9021$，以及 Parareal 的
$0.2822,0.4453,0.9986$。前两列大致满足“一轮 FCF 对应两轮
Parareal”；第三列都接近 $1$，表明小黏性下粗传播已不足以代表主导
对流。

![原论文 Figure 4.10：按两次细求解配平后的 Parareal 与 MGRiT 实测误差](assets/papers/time-parallelization/source-figures/figure-4-10.svg)

Figure 4.10 的四个面板按“热方程、$\nu=0.1$、$\nu=0.01$、$\nu=0.002$”排列，每个 Parareal 横坐标单位包含两轮，点划线标出离散截断误差 $\max\{\Delta t^2,\Delta x^2\}$，实践中不会迭代到它以下。前两幅中两条曲线接近；$\nu=0.01$ 时都很慢，Parareal 退化稍强，因为它在每次细传播后都插入一次失真的粗传播；$\nu=0.002$ 时两者发散，最大模态因子 Parareal 为 $1.4211$、MGRiT 为 $1.2812$（这两个数字只出现在原文正文中，Figure 4.9 只画了前三种扩散强度）。

非线性情形需要分别为 $\mathcal F$、$\mathcal G$ 配置非线性求解器。
理论比较还要求 $\mathcal F$、$\mathcal G$ 及二者差满足相应
Lipschitz 条件，所以“一轮 FCF 类似两轮 Parareal”并非无条件结论。
下面的 Burgers 实验使用齐次 Dirichlet 边界、同一初值、$T=5$、
$\Delta T=1/16$、$\Delta x=1/160$、$J=10$，中心空间差分、
后向 Euler 粗层和 SDIRK22 细层。

![原论文 Figure 4.11：三种黏性 Burgers 方程上等细求解成本的比较](assets/papers/time-parallelization/source-figures/figure-4-11.svg)

三个面板从左到右取 $\nu=0.5,0.01,0.002$。只要粗求解器仍足够
准确，三种黏性下的一轮 FCF-MGRiT 都大致对应两轮 Parareal，
与 Figure 4.10 的线性结果一致。$\nu=0.002$ 时两条曲线都先经历
较长的缓慢阶段，随后才快速下降，说明黏性降低会让两种方法同步恶化。

## 公式、定理与图表覆盖核对

| 原文项目                        | 论文小节 | 覆盖状态                                                               |
| ------------------------------- | -------- | ---------------------------------------------------------------------- |
| (4.10)                          | 4.3      | 细/粗配置节点、复制矩阵、配置系统                                      |
| Lagrange 传递与 PFASST 块迭代   | 4.3      | $T_{c\to f},T_{f\to c}$ 和 $\mathbf B_1^0,\mathbf B_0^1,\mathbf B_0^0$ |
| (4.11), Figure 4.6              | 4.3      | 隐式 Euler SDC、Radau 节点、PFASST 实验                                |
| (4.12), Figure 4.7              | 4.4      | FCF 更新、两次细求解、重叠和有限步性质                                 |
| (4.13), Theorem 4.5, Figure 4.8 | 4.4      | MGRiT 因子、与 Parareal 的关系、等成本收敛域                           |
| Theorem 4.6                     | 4.4      | 两类粗传播的四个最坏因子及平方比较                                     |
| Figures 4.9–4.11                | 4.4      | 线性谱、实测误差、非线性 Burgers 与弱扩散失效                          |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 4.3–4.4, pp. 452–461.
