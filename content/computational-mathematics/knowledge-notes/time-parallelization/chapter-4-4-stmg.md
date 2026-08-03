---
title: 4.6：时空多重网格（STMG）
description: 从全时间系统、时间块 Jacobi 平滑、时空传递算子与局部 Fourier 分析，到最优阻尼、积分器依赖与非线性 FAS 的完整推导与洞见
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-4-stmg
tags:
  - 时间并行
  - 时空多重网格
  - 局部-Fourier-分析
---

> [!note] 阅读范围
> 本页对应论文 Section 4.6（pp. 472–481），覆盖公式 (4.30)–(4.44)、Theorem 4.9、Figures 4.18–4.22 和 Table 4.1。线性两层循环、平滑符号、阻尼选择、积分器依赖和非线性 FAS 全部展开，并逐步给出「为什么」层面的推导。

## 4.6 时空多重网格（STMG）

时空多重网格（space-time multigrid, STMG）是本章最后引入的并行方法，其核心是把多重网格（MG）同时用在空间与时间两个方向上。经过 Hackbusch (1984) 和 Horton & Vandewalle (1995) 的早期奠基工作后，人们逐步认识到：**时间方向上的块 Jacobi 平滑器**是让整个方法既并行又高效的关键部件（Gander & Neumüller 2016）。有了它，STMG 只用标准多重网格构件，就能在抛物问题上达到与把多重网格用于 Poisson 方程时相当的效率——即迭代数几乎与网格无关。下面按论文顺序，从全时间系统出发逐层展开。

### 全时间系统 (4.30)–(4.31)

空间离散热方程或对流扩散方程后得到

$$
\boldsymbol u'=A\boldsymbol u+\boldsymbol f.
$$

这里 $A\in\mathbb R^{N_x\times N_x}$ 是空间导数的离散矩阵：热方程取 Laplace 算子 $\partial_{xx}$，对流扩散方程（ADE）取 $-\partial_x+\nu\partial_{xx}$，二者都以网格步长 $\Delta x$ 离散。对这组常微分方程施加一般单步积分器，写为

$$
r_1\boldsymbol u_{n+1}
=r_2\boldsymbol u_n+\widetilde{\boldsymbol f}_n,
\qquad n=0,\ldots,N_t-1, \tag{4.30}
$$

其中 $\boldsymbol u_0$ 为初值，$r_1,r_2$ 是 $\Delta tA$ 的矩阵多项式。后向 Euler 对应 $r_1=I_x-\Delta tA,\ r_2=I_x$；梯形规则对应 $r_1=I_x-\frac12\Delta tA,\ r_2=I_x+\frac12\Delta tA$。

与 [[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-5-paradiag-ii|ParaDiag]]（Section 3.5）一样，把全部 $N_t$ 个差分方程叠成一个「一次性（all-at-once）」系统：

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

$\boldsymbol b$ 是把初值与源项归并到右端得到的合适向量。矩阵 $K$ 是块下双对角（block lower bidiagonal）结构：主对角块全是 $r_1$，次对角块全是 $-r_2$，这正是时间步进「当前步依赖前一步」的代数体现。STMG 不再逐步顺序求解这个系统，而是在**空间和时间同时**建立粗网格层级，用一次多重网格循环整体求解 $\boldsymbol U$。

> [!tip] 本站洞见
> 从时间步进转到全时间系统 (4.31) 是理解所有时间并行方法的分水岭。顺序时间推进相当于用块前代（forward substitution）逐块求解块双对角 $K$，其数据依赖天然串行；一旦把 $K\boldsymbol U=\boldsymbol b$ 视为一个整体线性系统，我们就获得了在整个时空网格上重新分配计算的自由。STMG、ParaDiag、Parareal 的差异，本质上都是「如何近似求解同一个 $K$」的不同选择。

### 时间并行块 Jacobi 平滑器 (4.32)

STMG 在多重网格框架内求解 $\boldsymbol U$，其平滑器采用**阻尼块 Jacobi 迭代**。从初始近似 $\boldsymbol U^{\mathrm{ini}}$ 出发：

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

$s$ 是平滑次数，$\eta$ 是阻尼参数。之所以称为「块 Jacobi」，是因为迭代矩阵取 $K$ 的块对角部分 $I_t\otimes r_1$（丢弃次对角的 $-r_2$ 耦合），再用 $\eta$ 阻尼。这一步是 STMG 并行性的全部来源：

> [!tip] 本站洞见：为什么「时间方向的块 Jacobi」是并行的钥匙
> $K$ 的次对角块 $-r_2$ 编码了「时间步 $n+1$ 依赖时间步 $n$」的耦合。块 Jacobi 平滑器把这些耦合项**整体丢进残差** $\boldsymbol b-K\boldsymbol U^j$，而校正方程的左端只保留块对角 $I_t\otimes r_1$。由于 $I_t\otimes r_1$ 是块对角矩阵，求解 $\Delta\boldsymbol U^j$ 就分裂成 $N_t$ 个彼此独立的空间系统 $r_1\Delta\boldsymbol u_{n+1}^j=\eta(\cdots)_{n+1}$，可以在所有时间点上**完全并行**。相比之下，若像经典抛物多重网格那样用逐点 Gauss–Seidel（见 (4.35)），左端会包含时间耦合，必须按时间顺序求解，并行性荡然无存。换言之，STMG 把「时间耦合」从「必须串行求解的算子」降级为「只出现在残差里的显式项」，正是这一降级换来了时间并行。多重网格框架则负责补偿这一近似：块 Jacobi 高效压制高频误差，剩下的低频误差交给粗网格校正。

平滑器要在多重网格中工作，还需要空间和时间上的限制与延拓算子。以 $N_x=7$ 的空间线性插值为例：

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

$P_x\in\mathbb R^{7\times3}$ 是标准的一维线性插值（延拓），把粗网格 3 个内点插到细网格 7 个内点上；限制取 $R_x=\frac12P_x^\top$，即满权（full-weighting）限制。时间方向类似定义 $P_t,R_t$。这些都是**标准多重网格构件**，没有为时空问题做特殊改造——正是块 Jacobi 平滑器让标准构件足够用。

### 两层循环 (4.34)

记 `Mat` 把全时间向量重排成「空间 $\times$ 时间」矩阵，`Vec` 做逆操作（实现上即 Matlab 的 `reshape`）。这样空间算子从左乘、时间算子从右乘就能分别作用。从第 $k$ 到第 $k+1$ 迭代的一次两层循环为

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

其结构就是经典的两层多重网格 V-循环：先做 $s_1$ 次前平滑压制高频误差；再算残差并同时在空间（左乘 $R_x$）和时间（右乘 $R_t^\top$）限制到粗网格；在粗层求解误差方程 $K_c\boldsymbol e_c=\boldsymbol r_c$；把粗误差在空间（左乘 $P_x$）和时间（右乘 $P_t^\top$）延拓回细层并作校正；最后做 $s_2$ 次后平滑。

粗层算子 $K_c$ 与细层 $K$ 块结构完全相同，只是在放大一倍的时空步长 $\Delta T=2\Delta t$、$\Delta X=2\Delta x$ 上**重新离散**：

$$
K_c=
\begin{bmatrix}
r_1^c\\-r_2^c&r_1^c\\&\ddots&\ddots\\&&-r_2^c&r_1^c
\end{bmatrix},\qquad(N_t^c\ \text{块})
$$

其中 $r_1^c,r_2^c$ 是 $\Delta TA_c$ 的矩阵多项式，$A_c\in\mathbb R^{N_x^c\times N_x^c}$ 是空间导数在 $\Delta X$ 上的粗离散矩阵，例如

$$
\begin{cases}
r_1^c=I_x^c-\Delta TA_c,\quad r_2^c=I_x^c, & \text{后向 Euler},\\[2pt]
r_1^c=I_x^c-\tfrac12\Delta TA_c,\quad r_2^c=I_x^c+\tfrac12\Delta TA_c, & \text{梯形规则}.
\end{cases}
$$

实践中取 $N_x=2^{l_x}-1,\ N_t=2^{l_t}-1$（$l_x,l_t\ge2$），于是粗层大小为 $N_x^c=2^{l_x-1}-1$、$N_t^c=2^{l_t-1}-1$。把两层构造递归应用到粗层，就自然得到完整的多层 STMG。

> [!tip] 本站洞见
> $K_c$ 是「重新离散」（rediscretization）而非「Galerkin 三重积」$R K P$。这是因为 STMG 的时间方向可以看成一个强对流项（见下文历史评述），Galerkin 粗算子在这种非对称、类双曲的情形下不一定给出稳定的粗层积分器；直接在 $2\Delta t,2\Delta x$ 上重新离散能保证粗层仍是一个合理的时间步进算子，块结构与细层一致。

### 与早期抛物多重网格的差别 (4.35)

STMG 与四十年前 Hackbusch (1984) 提出的抛物多重网格有一处关键差别：后者用**逐点 Gauss–Seidel** 作平滑器，

$$
\boldsymbol U^{\mathrm{new}}
=S_{GS}(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s):
\left\{
\begin{aligned}
&\text{for }n=0,\ldots,N_t-1:\\
&\quad \boldsymbol u_{n+1}^0=\boldsymbol u_{n+1}^{\mathrm{ini}},\\
&\quad \text{for }j=0,\ldots,s-1:\\
&\qquad (D+L)\Delta\boldsymbol u_{n+1}^j=\widetilde{\boldsymbol f}_n+r_2\boldsymbol u_n^s-r_1\boldsymbol u_{n+1}^j,\\
&\qquad \boldsymbol u_{n+1}^{j+1}=\boldsymbol u_{n+1}^j+\Delta\boldsymbol u_{n+1}^j,\\
&\quad \boldsymbol u_{n+1}^{\mathrm{new}}=\boldsymbol u_{n+1}^s,
\end{aligned}
\right. \tag{4.35}
$$

其中 $\boldsymbol u_0^s=\boldsymbol u_0$，$D,L$ 分别是 $r_1$ 的对角与上三角部分。这里 $\boldsymbol U^{\mathrm{ini}}=(\boldsymbol u_0^\top,(\boldsymbol u_1^{\mathrm{ini}})^\top,\ldots)^\top$，$\boldsymbol U^{\mathrm{new}}$ 相应由 $\boldsymbol u_0$ 和 $\boldsymbol u_1^{\mathrm{new}},\ldots,\boldsymbol u_{N_t}^{\mathrm{new}}$ 组成。这个平滑器**在时间上严格顺序**：必须先在时间步 $n$ 完成平滑得到 $\boldsymbol u_n^s$，才能用它作为右端去平滑时间步 $n+1$。平滑之后，同样把残差 $\boldsymbol b-K\boldsymbol U^{\mathrm{new}}$ 在时空上限制到粗网格、求粗问题并递归。

Hackbusch (1984) 当时只在**空间**粗化，发现抛物多重网格对热方程收敛得非常快。Gander & Lunet (2024) 考察了同时在空间与时间粗化的两层版本，发现此时**收敛只会变慢**。这一慢收敛现象其实早已有解决之道：Horton & Vandewalle (1995) 通过采用把**时间方向解读为强对流项**、并据此专门设计的多重网格构件（各向异性/半粗化式的平滑与传递算子），显著改善了时空粗化下的收敛；相关的多重网格波形松弛（multigrid waveform relaxation）变体见 Janssen & Vandewalle (1996) 与 Van Lent & Vandewalle (2002)。STMG 走的是另一条互补路线——保留标准多重网格构件，而把改进集中在 (4.32) 的时间块 Jacobi 平滑器上：它让全体时间点并行，并在粗化前一致地压制高频误差。

> [!tip] 本站洞见
> Gauss–Seidel 与块 Jacobi 的对立恰好对应「串行精度」与「并行吞吐」的权衡。$S_{GS}$ 沿时间前扫，一遍就把信息从早期时刻传播到晚期时刻，单遍平滑更强、更接近顺序解，但完全不可时间并行；$S_\eta$ 放弃了这种时间内的信息传播（时间耦合被丢进残差），单遍平滑更弱，却换来 $N_t$ 路并行。多重网格恰恰弥补了后者的短板：被块 Jacobi 忽略的长程时间耦合本质是低频的，交给粗网格校正即可高效恢复，因此「弱平滑 + 粗化」的组合在保持并行的同时不损失整体效率。

### 局部 Fourier 分析的起点 (4.36)–(4.39)

设计有效平滑器的基本原则是：用尽量少的平滑次数消去尽量多的**高频**误差分量，使残余的低频误差能在粗网格上被准确表示并消去。分析这一点的利器是**局部 Fourier 分析（LFA）**——忽略初边值，只看有限差分模板如何作用于误差中的一个 Fourier 模态

$$
u_{n,m}^j=C_{\omega,\xi}^j
e^{i\omega n\Delta t}e^{i\xi m\Delta x}, \tag{4.36}
$$

其中 $\boldsymbol u_n^j:=(u_{n,1}^j,\ldots,u_{n,N_x}^j)^\top$，$i=\sqrt{-1}$，$\omega$ 是时间频率，$\xi$ 是空间频率。对一维热方程用中心差分（空间）+ 后向 Euler（时间）：

$$
A=\frac1{\Delta x^2}\operatorname{Tri}(1,-2,1),
\qquad r_1=I_x-\Delta tA,\quad r_2=I_x.
$$

把 (4.32) 右端 $\boldsymbol b$ 置零、视 $\boldsymbol U^j$ 为第 $j$ 次迭代的误差，块 Jacobi 的每个时间块方程为

$$
r_1(\boldsymbol u_{n+1}^{j+1}-\boldsymbol u_{n+1}^{j})
=-\eta(r_1\boldsymbol u_{n+1}^j-r_2\boldsymbol u_n^j). \tag{4.37}
$$

先算空间离散算子 $A$ 作用在 Fourier 模态 $u_{n+1,m}^l$（$l=j,j+1$）上的结果。利用相邻网格点带来的相位因子 $e^{\pm i\xi\Delta x}$ 并化简 $e^{-i\xi\Delta x}-2+e^{i\xi\Delta x}=2(\cos(\xi\Delta x)-1)$：

$$
Au_{n+1,m}^l
=\frac{2(\cos(\xi\Delta x)-1)}{\Delta x^2}
C_{\omega,\xi}^l
e^{i\omega(n+1)\Delta t}e^{i\xi m\Delta x},
\qquad l=j,j+1. \tag{4.38}
$$

于是左端

$$
r_1(u_{n+1}^{j+1}-u_{n+1}^j)
=\Big(1-\tfrac{2\Delta t(\cos(\xi\Delta x)-1)}{\Delta x^2}\Big)
(C_{\omega,\xi}^{j+1}-C_{\omega,\xi}^j)\,e^{i\omega(n+1)\Delta t}e^{i\xi x_h},
$$

而右端残差部分（注意 $u_n^j$ 相对 $u_{n+1}^j$ 差一个时间相位 $e^{-i\omega\Delta t}$）

$$
r_1u_{n+1}^j-r_2u_n^j
=\Big(1-e^{-i\omega\Delta t}-\tfrac{2\Delta t(\cos(\xi\Delta x)-1)}{\Delta x^2}\Big)
C_{\omega,\xi}^j\,e^{i\omega(n+1)\Delta t}e^{i\xi x_h}.
$$

代回 (4.37)，两端约去公共相位因子，得到振幅递推 $C^{j+1}=\rho C^j$，其中**收敛因子**

$$
\rho(\omega,\xi,\eta)
=1-\eta\left(
1-\frac{e^{-i\omega\Delta t}}
{1+\frac{2\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))}
\right), \tag{4.39}
$$

$\omega\Delta t\in(-\pi,\pi)$，$\xi\Delta x\in(-\pi,\pi)$。这一符号直接回答：给定 $\eta$ 后，任意时间/空间频率的误差经过一轮平滑被放大或压低多少。分母 $1+\frac{2\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))$ 恰是后向 Euler 隐式性带来的耗散因子——空间频率越高（$\cos(\xi\Delta x)\to-1$），分母越大，该分式越小，$\rho$ 越接近 $1-\eta$。

### Theorem 4.9：后向 Euler 下的最优阻尼

对 (4.39) 关于 $\xi,\omega$ 取最大再对该最大值取极小（即 min–max 优化），可证明如下结果（Gander & Lunet 2024, Ch. 4；更一般离散的完整分析见 Gander & Neumüller 2016）：

> [!note] Theorem 4.9
> 对一维热方程的中心差分–后向 Euler 离散，使阻尼 Jacobi 平滑器 (4.32) **始终允许时间粗化**的最优阻尼为
> $$\eta_{\mathrm{opt}}=\tfrac12.$$
> 此时所有时间高频 $\omega\in\pm\big(\tfrac{\pi}{2\Delta t},\tfrac{\pi}{\Delta t}\big)$ 都被压制至少 $\tfrac1{\sqrt2}$ 倍。若网格参数还满足
> $$\frac{\Delta t}{\Delta x^2}\ge\frac1{\sqrt2},$$
> 则空间高频 $\xi\in\pm\big(\tfrac{\pi}{2\Delta x},\tfrac{\pi}{\Delta x}\big)$ 也被压制至少 $\tfrac1{\sqrt2}$ 倍，因此可以**同时**在空间粗化。

关于最优性的一个更精细的分析见 Chaudet-Dumas, Gander & Pogozelskyte (2024)。

> [!tip] 本站洞见：粗化条件 $\Delta t/\Delta x^2\ge1/\sqrt2$ 从何而来
> 多重网格能收敛的前提是：**凡是粗网格无法表示的高频，平滑器都必须自己压下去**。时间高频总能被 $\eta=1/2$ 压到 $1/\sqrt2$ 以下，所以时间粗化「无条件」可行。空间高频却不然——由 (4.39)，空间频率只通过分母里的耗散因子 $\frac{2\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))$ 起作用，其强弱由无量纲数 $\Delta t/\Delta x^2$ 决定。只有当 $\Delta t/\Delta x^2\ge1/\sqrt2$、即隐式后向 Euler 在一步内对空间高频有足够耗散时，块 Jacobi 才能把空间高频也压到 $1/\sqrt2$ 以下，这时才敢同时粗化空间。若 $\Delta t/\Delta x^2$ 太小（时间步相对空间步过细），空间高频得不到充分平滑，就只应做时间粗化（半粗化）。这条件是 LFA 直接的、可解释的产物，而非经验拟合。

对中心差分 ADE，

$$
A=\frac{\nu}{\Delta x^2}\operatorname{Tri}(1,-2,1)
+\frac1{2\Delta x}\operatorname{Tri}(-1,0,1),
$$

同样计算 $A$ 作用在 Fourier 模态上的结果，扩散项给出实部 $\frac{2\nu(\cos(\xi\Delta x)-1)}{\Delta x^2}$，对流项的中心差分给出纯虚部 $\frac{i\sin(\xi\Delta x)}{\Delta x}$，于是符号变为

$$
\rho(\omega,\xi,\eta)
=1-\eta\left(
1-\frac{e^{-i\omega\Delta t}}
{1+\frac{2\nu\Delta t}{\Delta x^2}(1-\cos(\xi\Delta x))
+i\frac{\Delta t}{\Delta x}\sin(\xi\Delta x)}
\right). \tag{4.40}
$$

分母多出的虚部正是对流的贡献。启发式地，仍可取 $\eta=\tfrac12$ 作时间粗化的阻尼，见下面的图示。

![原论文 Figure 4.18：三种黏性 ADE 的高频最大平滑因子](assets/papers/time-parallelization/source-figures/figure-4-18.svg)

Figure 4.18 画的是收敛因子在高频区间上的最大值 $\rho_{\max}=\max_{(\xi\Delta x,\omega\Delta t)}\rho(\omega,\xi,\eta)$ 随 $\eta$ 的变化。从左到右取 $\nu=0.1,0.01,0.001$，每幅都比较 $\Delta x=\Delta t=1/64,1/128,1/256$。三种网格上的谷底都位于 $\eta\approx1/2$ 附近，说明这项选择对网格较稳健；黏性降低时，谷底从约 $0.71$ 抬升到约 $0.79$，可获得的高频收缩随之减弱——对流越占主导，标准 Jacobi 平滑越吃力。

### 阻尼、平滑次数和积分器依赖

![原论文 Figure 4.19：五、十、十五轮后误差随阻尼参数变化](assets/papers/time-parallelization/source-figures/figure-4-19.svg)

Figure 4.19 用数值实验进一步验证 $\eta=1/2$：两层 STMG 每轮只做一次块 Jacobi。(a) 是热方程，(b) 是 $\nu=0.01$ 的 ADE；每幅分别报告 5、10、15 轮后的误差随 $\eta$ 的曲线。迭代越多，$\eta=1/2$ 附近的低误差谷越清楚。热方程的谷底较宽，ADE 的 15 轮曲线最低点略偏向 $0.4$，因此 Figure 4.19 支持的是「$\eta=1/2$ 为稳健经验值」，而非声称它在每个有限网格和固定轮数下都是精确最优。

![原论文 Figure 4.20：一次与三次块 Jacobi 平滑的误差](assets/papers/time-parallelization/source-figures/figure-4-20.svg)

Figure 4.20 固定 $\eta=1/2$，考察平滑次数的影响。(a)、(b) 分别使用一次和三次块 Jacobi 平滑，每幅都比较热方程与 $\nu=0.1,0.01,0.001$ 的 ADE。两方程的两层 STMG 都随平滑次数增加而收敛更快；增加到三次平滑虽提高单轮成本，但显著减少循环数。ADE 的收敛率整体差于热方程，但当平滑次数较大时出现有趣现象：右图中多条 ADE 曲线后段变陡，进入**超线性**收敛阶段，对黏性 $\nu$ 的敏感性明显下降。

![原论文 Figure 4.21：梯形规则下不同平滑次数和阻尼的 STMG](assets/papers/time-parallelization/source-figures/figure-4-21.svg)

Figure 4.21 揭示 STMG 收敛率**依赖时间积分器**。换成梯形规则后，上排 (a) 的热方程依次使用 3、5、10 次平滑，所有阻尼扫描都停在较大误差甚至发散——即便加到 10 次平滑，也恢复不了后向 Euler 的效果。下排 (b) 是 $\nu=0.01$ 的 ADE，依次使用 2、3、4 次平滑，误差谷随平滑次数增加而加深，较优阻尼落在约 $\eta\approx0.8$，与后向 Euler 的 $\eta=1/2$ 明显不同。

> [!tip] 本站洞见：为什么后向 Euler 的 $\eta=1/2$ 结论不能移植到梯形规则
> Theorem 4.9 之所以成立，关键在后向 Euler 对高频**强耗散**：在 (4.39) 的分母里，时间高频经过一步后向 Euler 就被显著衰减，块 Jacobi 只需轻轻一压即可。梯形规则是 A-稳定但**非 L-稳定**的：其放大因子在 $\Delta t\lambda\to-\infty$ 时趋于 $-1$ 而非 $0$，最高频误差几乎不被时间离散本身耗散，反而每步反号振荡。这类「持久的时间高频」无法靠粗网格表示，又得不到积分器帮忙耗散，块 Jacobi 平滑于是失效——正是 Figure 4.21 上排热方程「无论怎样调 $\eta$ 都收敛困难」的根源。ADE 因为有对流带来的额外相位混合，情况稍好，但最优阻尼被推高到约 $0.8$。这说明 STMG 的平滑器与积分器必须匹配设计，是该方法「对积分器敏感」的定量来源。

![原论文 Table 4.1：三维热方程 STMG 的弱扩展和强扩展](assets/papers/time-parallelization/source-figures/table-4-1.svg)

Table 4.1（取自 Gander & Neumüller 2016）展示完整 STMG 在现代超算上求解三维热方程的弱扩展与强扩展。弱扩展从 1 核、2 个时间步、59,768 个自由度，一路增长到 262,144 核、524,288 个时间步、15,667,822,592 个自由度；**迭代数始终为 7**，墙钟时间从 28.8 秒仅微增到约 30.0 秒。作为对照，只做空间并行的顺序时间推进（列 `fwd. sub.`）估计从 19.0 秒暴涨到 4,988,060 秒。右侧强扩展在固定问题规模下把墙钟时间从约 7,635.2 秒逐级降到 30.0 秒。这张表说明 STMG 的价值来自**时空并行**与**网格无关迭代数**的同时实现——两者缺一，都无法在这种规模上保持常数墙钟时间。

### 非线性系统与 FAS (4.41)–(4.44)

STMG 可推广到非线性问题

$$
\boldsymbol u'=f(\boldsymbol u),
\qquad \boldsymbol u(0)=\boldsymbol u_0,
\quad t\in(0,T), \tag{4.41}
$$

其中 $f:\mathbb R^{N_x}\to\mathbb R^{N_x}$ 来自某个 PDE 的空间离散。对 (4.41) 施加线性-$\theta$ 方法，得到非线性全时间系统

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
(\boldsymbol u_0^\top+\Delta t(1-\theta)f(\boldsymbol u_0)^\top,0,\ldots,0)^\top,
\qquad
\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top.
$$

$B$ 编码时间差分的双对角耦合，$\widetilde B$ 按 $\theta$ 把非线性项 $f$ 在相邻两时刻加权。与 (4.32) 类比，先定义**非线性块 Jacobi 平滑器** $\boldsymbol U^{\mathrm{new}}=S_{\mathrm{non},\eta}(\boldsymbol b,\boldsymbol U^{\mathrm{ini}},s)$：

$$
\left\{
\begin{aligned}
\widetilde{\boldsymbol U}^0&=\boldsymbol U^{\mathrm{ini}},\\
\text{解 }\ \Delta\widetilde{\boldsymbol U}^j
-\Delta t\theta f(\Delta\widetilde{\boldsymbol U}^j)
&=\eta[\boldsymbol b-K(\widetilde{\boldsymbol U}^j)],
&&j=0,\ldots,s-1,\\
\widetilde{\boldsymbol U}^{j+1}
&=\widetilde{\boldsymbol U}^j+\Delta\widetilde{\boldsymbol U}^j,\\
\boldsymbol U^{\mathrm{new}}&=\widetilde{\boldsymbol U}^s.
\end{aligned}
\right. \tag{4.43}
$$

校正量 $\Delta\widetilde{\boldsymbol U}^j$ 由内层求解器（如 Newton 迭代）得到；由于左端仍是块对角（只在每个时间块内非线性耦合），$N_t$ 个局部非线性修正仍可在时间上并行。但非线性使 LFA 失效，故无法像 Theorem 4.9 那样给出 $\eta$ 的理论最优值，只能由实验或其他分析选取。

依照 Brandt (1977)，用**全近似格式（FAS）**定义非线性两层 STMG：

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

> [!tip] 本站洞见：FAS 的 $\boldsymbol r_c+K_c(\boldsymbol U_c^{k+1/3})$ 修正为何必要
> 线性情形下 $K(\boldsymbol U+\boldsymbol e)-K(\boldsymbol U)=K\boldsymbol e$，粗层可以直接解误差方程 $K_c\boldsymbol e_c=\boldsymbol r_c$。非线性算子没有这种叠加性：粗层无法单独表示「误差」。FAS 的办法是把限制来的当前解 $\boldsymbol U_c^{k+1/3}$ 也带到粗层，求解**完整解** $K_c(\boldsymbol U_c^{k+2/3})=\boldsymbol r_c+K_c(\boldsymbol U_c^{k+1/3})$，再取粗层解的**增量** $\boldsymbol e_c=\boldsymbol U_c^{k+2/3}-\boldsymbol U_c^{k+1/3}$ 作校正。右端加上的 $K_c(\boldsymbol U_c^{k+1/3})$ 正是使粗、细算子在当前解处保持一致的「$\tau$-修正」——没有它，粗层会在错误的工作点上线性化，破坏收敛。

![原论文 Figure 4.22：两次块 Jacobi 平滑的 Burgers STMG](assets/papers/time-parallelization/source-figures/figure-4-22.svg)

Figure 4.22 是 Burgers 方程 (2.6) 上的两层非线性 STMG，用两次平滑和经验最优 $\eta=1/4$。$\nu=1$ 的曲线在约 4 轮后越过图中离散误差线并继续降到 $10^{-4}$ 附近；$\nu=0.1$ 到第 16 轮仍高于该线。可见非线性 STMG 对黏性的依赖与线性情形（Figures 4.20–4.21）一致：黏性足够大时 STMG 在非线性设定下也很好用，黏性变小、方程更接近双曲时收敛显著恶化。

综合来看，对抛物问题，STMG 是当前**最有效**的时间并行求解器；但与 Parareal 不同，它**侵入性强**（需重写线性/非线性求解器、传递算子与全时间装配）。此外，如 Figure 4.21 上排所示，即便对抛物问题，STMG 的收敛率也**依赖所用的时间积分器**，这一依赖仍需进一步研究；而在双曲问题上（Figures 4.20、4.22），STMG 效率明显下降，表明该领域还需更多努力。

## 公式、定理与图表覆盖核对

| 原文项目                                | 论文小节 | 覆盖状态                                                       |
| --------------------------------------- | -------- | -------------------------------------------------------------- |
| (4.30)–(4.31)                           | 4.6      | 一般单步公式与全时间矩阵 $K$（块下双对角结构）                 |
| (4.32)–(4.34)                           | 4.6      | 并行块 Jacobi、时空传递 $P_x,R_x$、完整两层循环与粗算子 $K_c$  |
| (4.35)                                  | 4.6      | 早期顺序 Gauss–Seidel 平滑及并行性差别、Horton–Vandewalle 改进 |
| (4.36)–(4.40), Theorem 4.9, Figure 4.18 | 4.6      | 完整 LFA 推导、热/ADE 符号、最优阻尼 $\eta=1/2$、粗化条件      |
| Figures 4.19–4.21, Table 4.1            | 4.6      | 阻尼扫描、平滑次数、积分器依赖、弱/强扩展                      |
| (4.41)–(4.44), Figure 4.22              | 4.6      | 非线性全时间系统、并行非线性平滑、FAS 一致性、Burgers 实验     |

## 参考补充

- 时空多重网格奠基与关键部件：W. Hackbusch (1984) 抛物多重网格；G. Horton & S. Vandewalle (1995) 将时间方向解读为强对流项以改善时空粗化收敛；C. Janssen & S. Vandewalle (1996)、S. Van Lent & S. Vandewalle (2002) 的多重网格波形松弛变体；M. J. Gander & M. Neumüller (2016) 确立时间块 Jacobi 平滑并给出大规模扩展结果。
- 收敛与最优性分析：M. J. Gander & T. Lunet (2024) 的两层 LFA 与 Theorem 4.9 证明；关于最优性的精细分析见 B. Chaudet-Dumas, M. J. Gander & A. Pogozelskyte (2024)。
- 非线性框架：A. Brandt (1977) 的全近似格式（FAS）。

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4.6, pp. 472–481.
