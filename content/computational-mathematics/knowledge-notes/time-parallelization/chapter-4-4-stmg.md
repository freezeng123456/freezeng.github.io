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

## 4.6.1 全时间系统 (4.30)–(4.31)

空间离散热方程或对流扩散方程后得到

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

## 4.6.2 时间并行块 Jacobi 平滑器

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

## 4.6.3 两层循环 (4.34)

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

## 4.6.4 与早期抛物多重网格的差别

Hackbusch 型抛物多重网格采用逐时间点、逐空间点的 Gauss–Seidel 平滑：

$$
\boldsymbol U^{\mathrm{new}}
=S_{GS}(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s), \tag{4.35}
$$

每个时间步内求 $(D+L)\Delta\boldsymbol u_{n+1}^j$，并立即更新 $\boldsymbol u_{n+1}^{j+1}$。下一个时间点依赖当前点平滑完成后的值，故时间方向严格顺序。它在只粗化空间的热方程上很快，同时粗化时空时会变慢。STMG 的关键变化是 (4.32) 的时间块 Jacobi，使全体时间点并行，并让高频误差在粗化前得到一致压制。

## 4.6.5 局部 Fourier 分析的起点

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

这一符号直接回答：给定 $\eta$ 后，时间/空间高频能被一轮平滑压低多少。

## 4.6.6 Theorem 4.9：后向 Euler 下的最优阻尼

对一维热方程的中心差分–后向 Euler 离散，始终允许时间粗化的最优阻尼为

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

![原论文 Figure 4.18：三种黏性 ADE 的高频最大平滑因子](assets/papers/time-parallelization/source-figures/figure-4-18.svg)

Figure 4.18 扫描空间与时间高频。即使含对流，$\eta=1/2$ 仍是时间粗化的合理经验选择；黏性降低会抬高最坏因子。

## 4.6.7 阻尼、平滑次数和积分器依赖

![原论文 Figure 4.19：五、十、十五轮后误差随阻尼参数变化](assets/papers/time-parallelization/source-figures/figure-4-19.svg)

两层 STMG 每轮只做一次块 Jacobi。热方程和 $\nu=0.01$ 的 ADE 都在 $\eta=1/2$ 附近达到较小误差，验证 LFA 选择。

![原论文 Figure 4.20：一次与三次块 Jacobi 平滑的误差](assets/papers/time-parallelization/source-figures/figure-4-20.svg)

增加到三次平滑会提高单轮成本，同时显著减少循环数。ADE 比热方程慢；平滑次数增加后对黏性的敏感性下降，并出现超线性阶段。

![原论文 Figure 4.21：梯形规则下不同平滑次数和阻尼的 STMG](assets/papers/time-parallelization/source-figures/figure-4-21.svg)

换成梯形规则后，热方程即使做十次平滑也有明显收敛问题。ADE 可以收敛，较多平滑有帮助，采样中的较优阻尼约为 $0.8$。后向 Euler 的 $\eta=1/2$ 结论依赖其高频耗散，不能直接移植。

![原论文 Table 4.1：三维热方程 STMG 的弱扩展和强扩展](assets/papers/time-parallelization/source-figures/table-4-1.svg)

弱扩展从 1 核、2 个时间步、59,768 个自由度增长到 262,144 核、524,288 个时间步、15,667,822,592 个自由度；迭代数始终为 7，墙钟时间从 28.8 秒维持到约 30.0 秒。只做空间并行的顺序时间推进估计从 19.0 秒增至 4,988,060 秒。强扩展在固定问题规模下也从约 7,635.2 秒降到 30.0 秒。这张表说明 STMG 的价值来自时空并行与网格无关迭代数的同时实现。

## 4.6.8 非线性系统与 FAS

考虑

$$
\boldsymbol u'=f(\boldsymbol u),
\qquad \boldsymbol u(0)=\boldsymbol u_0,
\quad t\in(0,T). \tag{4.41}
$$

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
\right. \tag{4.43}
$$

每个时间块的非线性修正可用内层 Newton，并在时间上并行。非线性使 LFA 不再适用，$\eta$ 需要实验或其他分析选择。

FAS 两层循环为

$$
\left\{
\begin{aligned}
\boldsymbol U^{k+1/3}&=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^k,s_1),\\
\boldsymbol r&=\boldsymbol b-K(\boldsymbol U^{k+1/3}),\\
\boldsymbol r_c&=[R_x\operatorname{Mat}(\boldsymbol r)]R_t^\top,\\
\boldsymbol U_c^{k+1/3}&=[R_x\operatorname{Mat}(\boldsymbol U^{k+1/3})]R_t^\top,\\
K_c(\boldsymbol U_c^{k+2/3})
&=\boldsymbol r_c+K_c(\boldsymbol U_c^{k+1/3}),\\
\boldsymbol e_c&=\boldsymbol U_c^{k+2/3}-\boldsymbol U_c^{k+1/3},\\
\boldsymbol e&=[P_x\operatorname{Mat}(\boldsymbol e_c)]P_t^\top,\\
\boldsymbol U^{k+2/3}&=\boldsymbol U^{k+1/3}+\operatorname{Vec}(\boldsymbol e),\\
\boldsymbol U^{k+1}&=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^{k+2/3},s_2).
\end{aligned}
\right. \tag{4.44}
$$

FAS 在粗网格上求完整近似，并通过 $\boldsymbol r_c+K_c(\boldsymbol U_c)$ 保持非线性一致性；这与线性情形只解粗误差方程不同。

![原论文 Figure 4.22：两次块 Jacobi 平滑的 Burgers STMG](assets/papers/time-parallelization/source-figures/figure-4-22.svg)

实验用两次平滑和经验最优 $\eta=1/4$。扩散充分时收敛很快；黏性降低后明显恶化。STMG 是论文所述抛物问题中最强的时间并行求解器，同时侵入性高于 Parareal，对积分器与方程类型也更敏感。

## 公式、定理与图表覆盖核对

| 原文项目                                | 本页位置  | 覆盖状态                                        |
| --------------------------------------- | --------- | ----------------------------------------------- |
| (4.30)–(4.31)                           | §4.6.1    | 一般单步公式与全时间矩阵                        |
| (4.32)–(4.34)                           | §§4.6.2–3 | 并行块 Jacobi、时空传递和完整两层循环           |
| (4.35)                                  | §4.6.4    | 早期顺序 Gauss–Seidel 平滑及差别                |
| (4.36)–(4.40), Theorem 4.9, Figure 4.18 | §§4.6.5–6 | 完整 LFA、热/ADE 符号、最优阻尼和粗化条件       |
| Figures 4.19–4.21, Table 4.1            | §4.6.7    | 阻尼扫描、平滑次数、积分器依赖、强弱扩展        |
| (4.41)–(4.44), Figure 4.22              | §4.6.8    | 非线性全时间系统、并行平滑、FAS 和 Burgers 实验 |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4.6, pp. 472–481.
