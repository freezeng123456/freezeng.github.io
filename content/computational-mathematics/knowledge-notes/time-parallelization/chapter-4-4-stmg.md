---
title: 4.6：时空多重网格（STMG）
description: 从全时间系统、时间块 Jacobi 平滑与局部 Fourier 分析到非线性 FAS 的完整推导
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-4-stmg
tags:
  - 时间并行
  - 时空多重网格
  - 局部-Fourier-分析
---

> [!note] 阅读范围
> 本页对应论文 Section 4.6（pp. 472–481），覆盖公式 (4.30)–(4.44)、Theorem 4.9、Figures 4.18–4.22 和 Table 4.1。线性两层循环、平滑符号、阻尼选择、积分器依赖和非线性 FAS 全部展开。

## 4.6 时空多重网格（STMG）

时空多重网格的历史可上溯到 Hackbusch（1984）的抛物多重网格和 Horton 与 Vandewalle（1995）的改进。本节的 STMG 来自 Gander 与 Neumüller（2016），其核心发现是：只要把时间方向的平滑换成块 Jacobi，就能只用标准多重网格分量达到 Poisson 问题上的效率。相关工作还有 Janssen 与 Vandewalle（1996）、Van Lent 与 Vandewalle（2002）以及 Chaudet-Dumas、Gander 与 Pogozelskyte（2024）。

### 全时间系统 (4.30)–(4.31)

与 Section 3.5 的 ParaDiag 一样，先把所有时间点的未知量收进一个系统。空间离散热方程或对流扩散方程后得到

$$
\boldsymbol u'=A\boldsymbol u+\boldsymbol f.
$$

一般单步积分器写为

$$
r_1\boldsymbol u_{n+1}
=r_2\boldsymbol u_n+\widetilde{\boldsymbol f}_n,
\qquad n=0,\ldots,N_t-1, \tag{4.30}
$$

$r_1,r_2$ 是 $\Delta tA$ 的矩阵多项式。后向 Euler 对应 $r_1=I_x-\Delta tA,r_2=I_x$；梯形规则对应 $r_1=I_x-\frac12\Delta tA,r_2=I_x+\frac12\Delta tA$。

叠起所有时间点：

$$
\underbrace{
\begin{bmatrix}
r_1\\-r_2&r_1\\&\ddots&\ddots\\&&-r_2&r_1
\end{bmatrix}}_{K}
\underbrace{
\begin{bmatrix}
\boldsymbol u_1\\\boldsymbol u_2\\\vdots\\\boldsymbol u_{N_t}
\end{bmatrix}}_{\boldsymbol U}
=\boldsymbol b. \tag{4.31}
$$

STMG 在空间和时间同时建立粗网格，并以多重网格循环解该全时间系统。

### 时间并行块 Jacobi 平滑器

阻尼块 Jacobi 定义为

$$
\boldsymbol U^{\mathrm{new}}
=S_\eta(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s):
\left\{
\begin{aligned}
\boldsymbol U^0&=\boldsymbol U^{\mathrm{ini}},\\
(I_t\otimes r_1)\Delta\boldsymbol U^j
&=\eta(\boldsymbol b-K\boldsymbol U^j),
&&j=0,\ldots,s-1,\\
\boldsymbol U^{j+1}&=\boldsymbol U^j+\Delta\boldsymbol U^j,\\
\boldsymbol U^{\mathrm{new}}&=\boldsymbol U^s.
\end{aligned}
\right. \tag{4.32}
$$

$I_t\otimes r_1$ 块对角，所以每次平滑的 $N_t$ 个空间系统完全独立。$\eta$ 控制阻尼，$s$ 是平滑次数。

空间线性插值在 $N_x=7$ 的示例中为

$$
P_x=
\begin{bmatrix}
1/2&0&0\\
1&0&0\\
1/2&1/2&0\\
0&1&0\\
0&1/2&1/2\\
0&0&1\\
0&0&1/2
\end{bmatrix},
\qquad
R_x=\frac12P_x^\top. \tag{4.33}
$$

时间方向类似定义 $P_t,R_t$。

### 两层循环 (4.34)

记 `Mat` 把全时间向量重排成空间×时间矩阵，`Vec` 做逆操作。一次两层循环为

$$
\left\{
\begin{aligned}
\boldsymbol U^{k+1/3}&=S_\eta(\boldsymbol b,\boldsymbol U^k,s_1),\\
\boldsymbol r&=\boldsymbol b-K\boldsymbol U^{k+1/3},\\
\boldsymbol r_c&=[R_x\operatorname{Mat}(\boldsymbol r)]R_t^\top,\\
\boldsymbol e_c&=K_c^{-1}\operatorname{Vec}(\boldsymbol r_c),\\
\boldsymbol e&=[P_x\operatorname{Mat}(\boldsymbol e_c)]P_t^\top,\\
\boldsymbol U^{k+2/3}&=\boldsymbol U^{k+1/3}+\operatorname{Vec}(\boldsymbol e),\\
\boldsymbol U^{k+1}&=S_\eta(\boldsymbol b,\boldsymbol U^{k+2/3},s_2).
\end{aligned}
\right. \tag{4.34}
$$

$K_c$ 在 $\Delta T=2\Delta t$、$\Delta X=2\Delta x$ 上重新离散，块结构与 $K$ 相同。若 $N_x=2^{l_x}-1,N_t=2^{l_t}-1$，粗层大小为 $2^{l_x-1}-1$ 与 $2^{l_t-1}-1$。递归应用两层构造便得到完整 STMG。

### 与早期抛物多重网格的差别

Hackbusch 型抛物多重网格采用逐时间点、逐空间点的 Gauss–Seidel
平滑。完整循环是

$$
\boldsymbol U^{\mathrm{new}}
=S_{GS}(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s):
\left\{
\begin{aligned}
&\text{for }n=0,\ldots,N_t-1,\\
&\qquad \boldsymbol u_{n+1}^0=\boldsymbol u_{n+1}^{\mathrm{ini}},\\
&\qquad \text{for }j=0,\ldots,s-1,\\
&\qquad\quad
(D+L)\Delta\boldsymbol u_{n+1}^j
=\widetilde{\boldsymbol f}_n+r_2\boldsymbol u_n^s
-r_1\boldsymbol u_{n+1}^j,\\
&\qquad\quad
\boldsymbol u_{n+1}^{j+1}
=\boldsymbol u_{n+1}^j+\Delta\boldsymbol u_{n+1}^j,\\
&\qquad \boldsymbol u_{n+1}^{\mathrm{new}}
=\boldsymbol u_{n+1}^s,
\end{aligned}
\right. \tag{4.35}
$$

这里 $\boldsymbol u_0^s=\boldsymbol u_0$，且
$\boldsymbol U^{\mathrm{ini}}$、$\boldsymbol U^{\mathrm{new}}$ 都包含初值
$\boldsymbol u_0$；$D$ 与 $L$ 是 $r_1$ 的对角部分和**上**三角
部分。下一个时间点依赖当前点平滑完成后的值，故时间方向严格顺序。
它在只粗化空间的热方程上很快，同时粗化时空时会变慢。Horton 与
Vandewalle（1995）曾把时间方向解释成强对流项，并用专门的多重网格
分量改善慢收敛。STMG 的关键变化则是 (4.32) 的时间块 Jacobi
（Gander 与 Neumüller 2016）：所有时间点并行平滑，只用标准多重
网格分量也能达到类似 Poisson 问题的效果。

### 局部 Fourier 分析的起点

忽略初边值，在无限规则网格上考察误差模态

$$
u_{n,m}^j=C_{\omega,\xi}^j
e^{i\omega n\Delta t}e^{i\xi m\Delta x}. \tag{4.36}
$$

$\omega$ 是时间频率，$\xi$ 是空间频率。对一维热方程，

$$
A=\frac1{\Delta x^2}\operatorname{Tri}(1,-2,1),
\qquad r_1=I_x-\Delta tA,\quad r_2=I_x.
$$

齐次误差上的块 Jacobi 为

$$
r_1(\boldsymbol u_{n+1}^{j+1}-\boldsymbol u_{n+1}^{j})
=-\eta(r_1\boldsymbol u_{n+1}^j-r_2\boldsymbol u_n^j). \tag{4.37}
$$

空间离散算子作用在 Fourier 模态上给出

$$
Au_{n+1,m}^l
=\frac{2(\cos(\xi\Delta x)-1)}{\Delta x^2}
C_{\omega,\xi}^l
e^{i\omega(n+1)\Delta t}e^{i\xi m\Delta x},
\qquad l=j,j+1. \tag{4.38}
$$

代回 (4.37)，振幅满足 $C^{j+1}=\rho C^j$，其中

$$
\rho(\omega,\xi,\eta)
=1-\eta\left(
1-\frac{e^{-i\omega\Delta t}}
{1+\frac{2\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))}
\right). \tag{4.39}
$$

放大因子是 $|\rho|$。为判断能否做时间粗化，定义时间高频集合

$$
\Theta_t^{\mathrm{high}}
=(-\pi,-\pi/2)\cup(\pi/2,\pi)
$$

以及平滑因子

$$
\mu_t(\eta)=
\sup_{\substack{\xi\Delta x\in(-\pi,\pi)\\
\omega\Delta t\in\Theta_t^{\mathrm{high}}}}
|\rho(\omega,\xi,\eta)|.
$$

若把零频也纳入最大值，则 (4.39) 在 $(\omega,\xi)=(0,0)$ 恒有
$\rho=1$，无法用来选择唯一阻尼。

### Theorem 4.9：后向 Euler 下的最优阻尼

Theorem 4.9 的证明见 Gander 与 Lunet（2024）。这里的“最优”是
最小化高时间频率平滑因子 $\mu_t(\eta)$。对一维热方程的
中心差分–后向 Euler 离散，始终允许时间粗化的最优阻尼为

$$
\eta_{\mathrm{opt}}=\frac12.
$$

时间高频 $\omega\in\pm(\pi/(2\Delta t),\pi/\Delta t)$ 的放大因子不超过 $1/\sqrt2$。若

$$
\frac{\Delta t}{\Delta x^2}\ge\frac1{\sqrt2},
$$

空间高频 $\xi\in\pm(\pi/(2\Delta x),\pi/\Delta x)$ 也不超过 $1/\sqrt2$，因此可同时粗化空间。

对中心差分 ADE，

$$
A=\frac{\nu}{\Delta x^2}\operatorname{Tri}(1,-2,1)
+\frac1{2\Delta x}\operatorname{Tri}(-1,0,1),
$$

Fourier 符号变为

$$
\rho(\omega,\xi,\eta)
=1-\eta\left(
1-\frac{e^{-i\omega\Delta t}}
{1+\frac{2\nu\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))
+i\frac{\Delta t}{\Delta x}\sin(\xi\Delta x)}
\right). \tag{4.40}
$$

> [!warning] 原文公式核对：ADE 矩阵与 (4.40) 的符号
> 上面的 $A$ 与 (4.40) 不能同时正确。按论文自己的约定，$\operatorname{Tri}(-1,0,1)/(2\Delta x)$ 作用在 Fourier 模态上给出 $+i\sin(\xi\Delta x)/\Delta x$，于是 $r_1=I_x-\Delta tA$ 的符号与 (4.40) 分母中的 $+i\frac{\Delta t}{\Delta x}\sin(\xi\Delta x)$ 相反。模型问题 (2.5) 是 $\partial_tu+\partial_xu-\nu\partial_{xx}u=g$，所以 $A$ 离散的是 $-\partial_x+\nu\partial_{xx}$，正确写法应为 $A=\frac{\nu}{\Delta x^2}\operatorname{Tri}(1,-2,1)-\frac1{2\Delta x}\operatorname{Tri}(-1,0,1)$，此时 (4.40) 如印刷所示成立。上式保留原文排印。由于 $\rho_{\max}$ 是对 $\xi\Delta x\in(-\pi,\pi)$ 取最大而 $\sin$ 是奇函数，这个符号不影响本页引用的任何数值。

![原论文 Figure 4.18：三种黏性 ADE 的高频最大平滑因子](assets/papers/time-parallelization/source-figures/figure-4-18.svg)

Figure 4.18 从左到右取 $\nu=0.1,0.01,0.001$，每幅都比较
$\Delta x=\Delta t=1/64,1/128,1/256$。纵轴应理解为

$$
\rho_{\max}=
\max_{\substack{\Delta x\xi\in(-\pi,\pi)\\
\Delta t\omega\in\Theta_t^{\mathrm{high}}}}
|\rho(\omega,\xi,\eta)|,
$$

即对全部空间频率、但只对时间高频取模后最大化。期刊图注只写了
正半区且省略模长；由共轭对称性可补成上式。同一面板内三种网格的谷底
基本重合；随 $\nu$ 减小，最优点从约 $0.5$ 缓慢右移，谷底值从约
$0.71$ 抬升到约 $0.79$。因此 $\eta=\tfrac12$ 是稳健起点，
不保证对每个 $\nu$ 都精确最优。

### 阻尼、平滑次数和积分器依赖

![原论文 Figure 4.19：五、十、十五轮后误差随阻尼参数变化](assets/papers/time-parallelization/source-figures/figure-4-19.svg)

两层 STMG 每轮只做一次块 Jacobi。(a) 是热方程，(b) 是 $\nu=0.01$ 的 ADE；每幅分别报告 5、10、15 轮后的误差。迭代越多，$\eta$ 附近的低误差谷越清楚。热方程的谷底较宽，ADE 的 15 轮曲线谷底在约 $0.3$ 到 $0.5$ 之间相当平坦，因此 Figure 4.19 支持的是“$\eta=1/2$ 为稳健经验值”，并没有声称它在每个有限网格和固定轮数下都是精确最优值。

![原论文 Figure 4.20：一次与三次块 Jacobi 平滑的误差](assets/papers/time-parallelization/source-figures/figure-4-20.svg)

(a)、(b) 分别使用一次和三次块 Jacobi 平滑，每幅都比较热方程与 $\nu=0.1,0.01,0.001$ 的 ADE。增加到三次平滑会提高单轮成本，同时显著减少循环数；右图中多条 ADE 曲线后段变陡，出现超线性阶段。ADE 仍比热方程慢，但对黏性的敏感性明显下降。

![原论文 Figure 4.21：梯形规则下不同平滑次数和阻尼的 STMG](assets/papers/time-parallelization/source-figures/figure-4-21.svg)

换成梯形规则后，上排热方程依次使用 3、5、10 次平滑，所有阻尼扫描
都停在较大误差或发生发散；增加平滑次数没有恢复后向 Euler 的效果。
下排是 $\nu=0.01$ 的 ADE，依次使用 2、3、4 次平滑，误差谷随平滑
次数增加而加深。原文把下排概括为经验值 $\eta\approx0.8$；逐面板
读取则显示最优点从约 $0.79$ 移到 $0.92$。前者适合作为初始选择，
后者是各固定面板的离散最小值。后向 Euler 的 $\eta=1/2$ 依赖其
高频耗散，不能直接移植。

![原论文 Table 4.1：三维热方程 STMG 的弱扩展和强扩展](assets/papers/time-parallelization/source-figures/table-4-1.svg)

该表取自 Gander 与 Neumüller（2016）。弱扩展从 1 核、2 个时间步、59,768 个自由度增长到 262,144 核、524,288 个时间步、15,667,822,592 个自由度；迭代数始终为 7，墙钟时间从 28.8 秒维持到约 30.0 秒。作为参照的“仅空间并行的经典时间推进”从 19.0 秒增至 4,988,060 秒。强扩展有两组：512 个时间步、15,300,608 个自由度的一组从 7,635.2 秒降到 30.0 秒；524,288 个时间步、15,667,822,592 个自由度的一组从 15,205.9 秒降到 30.0 秒。这张表说明 STMG 的价值来自时空并行与网格无关迭代数的同时实现。

### 非线性系统与 FAS

考虑

$$
\boldsymbol u'=f(\boldsymbol u),
\qquad \boldsymbol u(0)=\boldsymbol u_0,
\quad t\in(0,T). \tag{4.41}
$$

对堆叠向量明确约定

$$
\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top,
\qquad
f(\boldsymbol U)=
\left(
f(\boldsymbol u_1)^\top,\ldots,f(\boldsymbol u_{N_t})^\top
\right)^\top.
$$

期刊版在这里把 $f(\boldsymbol U)$ 误排成了
$(\boldsymbol u_1,\ldots,\boldsymbol u_{N_t})$ 本身。

线性-$\theta$ 离散形成

$$
\underbrace{(B\otimes I_x)\boldsymbol U
-\Delta t(\widetilde B\otimes I_x)f(\boldsymbol U)}_{K(\boldsymbol U)}
=\boldsymbol b, \tag{4.42}
$$

其中

$$
B=
\begin{bmatrix}
1\\-1&1\\&\ddots&\ddots\\&&-1&1
\end{bmatrix},
\qquad
\widetilde B=
\begin{bmatrix}
\theta\\1-\theta&\theta\\&\ddots&\ddots\\&&1-\theta&\theta
\end{bmatrix},
$$

$$
\boldsymbol b=
(\boldsymbol u_0^\top+\Delta t(1-\theta)f(\boldsymbol u_0)^\top,0,\ldots,0)^\top.
$$

非线性块 Jacobi 为

$$
\boldsymbol U^{\mathrm{new}}
=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s):
\left\{
\begin{aligned}
\widetilde{\boldsymbol U}^0&=\boldsymbol U^{\mathrm{ini}},\\
\Delta\widetilde{\boldsymbol U}^j
-\Delta t\theta f(\Delta\widetilde{\boldsymbol U}^j)
&=\eta[\boldsymbol b-K(\widetilde{\boldsymbol U}^j)],\\
\widetilde{\boldsymbol U}^{j+1}
&=\widetilde{\boldsymbol U}^j+\Delta\widetilde{\boldsymbol U}^j,\\
\boldsymbol U^{\mathrm{new}}&=\widetilde{\boldsymbol U}^s.
\end{aligned}
\right. \tag{4.43, as printed}
$$

> [!warning] 原文公式核对：(4.43) 的“块 Jacobi”含义
> 若 $\Delta\widetilde{\boldsymbol U}^j$ 是随后加到当前解上的修正，
> 对角块
> $D(\boldsymbol U)=\boldsymbol U-\Delta t\theta f(\boldsymbol U)$
> 的一致非线性修正应满足
> $$
> \Delta\widetilde{\boldsymbol U}^j
> -\Delta t\theta\left[
> f(\widetilde{\boldsymbol U}^j+\Delta\widetilde{\boldsymbol U}^j)
> -f(\widetilde{\boldsymbol U}^j)
> \right]
> =\eta[\boldsymbol b-K(\widetilde{\boldsymbol U}^j)].
> $$
> 原文使用 $f(\Delta\widetilde{\boldsymbol U}^j)$，只在线性等特殊
> 情形下与上式相同。因此 (4.43) 更安全的理解是一种另行定义的
> 非线性预条件迭代，而不是精确块 Jacobi 逆。

每个时间块的非线性修正都可用内层 Newton，并在时间上并行。非线性使
LFA 不再适用，$\eta$ 需要实验或其他分析选择。

FAS 两层循环为

$$
\left\{
\begin{aligned}
\boldsymbol U^{k+1/3}&=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^k,s_1),\\
\boldsymbol r&=\boldsymbol b-K(\boldsymbol U^{k+1/3}),\\
\boldsymbol r_c&=\operatorname{Vec}\!\left(
[R_x\operatorname{Mat}(\boldsymbol r)]R_t^\top
\right),\\
\boldsymbol U_c^{k+1/3}&=\operatorname{Vec}\!\left(
[R_x\operatorname{Mat}(\boldsymbol U^{k+1/3})]R_t^\top
\right),\\
K_c(\boldsymbol U_c^{k+2/3})
&=\boldsymbol r_c+K_c(\boldsymbol U_c^{k+1/3}),\\
\boldsymbol e_c&=\boldsymbol U_c^{k+2/3}-\boldsymbol U_c^{k+1/3},\\
\boldsymbol e&=\operatorname{Vec}\!\left(
[P_x\operatorname{Mat}(\boldsymbol e_c)]P_t^\top
\right),\\
\boldsymbol U^{k+2/3}&=\boldsymbol U^{k+1/3}+\boldsymbol e,\\
\boldsymbol U^{k+1}&=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^{k+2/3},s_2).
\end{aligned}
\right. \tag{4.44}
$$

FAS（Brandt 1977）在粗网格上求完整近似，并通过
$\boldsymbol r_c+K_c(\boldsymbol U_c)$ 保持非线性一致性；这与线性
情形只解粗误差方程不同。上式显式补出 `Vec`，使 $K_c$ 的输入和
粗残差都保持为向量；期刊版在矩阵与向量之间省略了这一步。

![原论文 Figure 4.22：两次块 Jacobi 平滑的 Burgers STMG](assets/papers/time-parallelization/source-figures/figure-4-22.svg)

实验用两次平滑和经验最优 $\eta=1/4$。$\nu=1$ 的曲线在约 4 轮后越过图中的离散误差线，并继续降到 $10^{-4}$ 附近；$\nu=0.1$ 到第 16 轮仍高于该线。非线性 STMG 对黏性的依赖因此与线性情形（Figure 4.20 的 $\nu$ 扫描）一致。

在综述覆盖的抛物方法中，STMG 被评价为最有效，但侵入性也高于
Parareal。Figures 4.20 和 4.22 实际测试的是低扩散 ADE 与黏性
Burgers，而非严格双曲方程；它们直接支持的是**近双曲极限下的退化**。
综述作者据此外推：真正双曲问题上的效率仍需进一步研究。Figure 4.21
则表明，即使对抛物问题，时间积分器依赖也仍是开放问题。

## 公式、定理与图表覆盖核对

| 原文项目                                | 论文小节 | 覆盖状态                                        |
| --------------------------------------- | -------- | ----------------------------------------------- |
| (4.30)–(4.31)                           | 4.6      | 一般单步公式与全时间矩阵                        |
| (4.32)–(4.34)                           | 4.6      | 并行块 Jacobi、时空传递和完整两层循环           |
| (4.35)                                  | 4.6      | 早期顺序 Gauss–Seidel 平滑及差别                |
| (4.36)–(4.40), Theorem 4.9, Figure 4.18 | 4.6      | 完整 LFA、热/ADE 符号、最优阻尼和粗化条件       |
| Figures 4.19–4.21, Table 4.1            | 4.6      | 阻尼扫描、平滑次数、积分器依赖、强弱扩展        |
| (4.41)–(4.44), Figure 4.22              | 4.6      | 非线性全时间系统、并行平滑、FAS 和 Burgers 实验 |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4.6, pp. 472–481.
