---
title: 3.5.1：直接 ParaDiag（ParaDiag-I）
description: 从全时间系统到几何时间网格、BVM 和非线性准 Newton 的完整推导
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i
tags:
  - 时间并行
  - ParaDiag
  - 全时间系统
---

> [!note] 阅读范围
> 本页对应论文 Sections 3.5–3.5.1（pp. 415–430），覆盖 ParaDiag 总体分类、公式 (3.22)–(3.48)、Theorems 3.5–3.7、Figures 3.9–3.14 和 Tables 3.1–3.2。直接法的截断误差、舍入误差、BVM 以及非线性准 Newton 分支均保留。

## Section 3.5 导论：ParaDiag-I 与 ParaDiag-II 的分界

ParaDiag-I 以特殊、可对角化的时间离散换取无外层迭代的直接并行解。
为此必须使用变时间步长，或在末步换用另一种积分公式。代价有两项：
舍入误差限制单窗内的并行步数，双精度几何网格通常只能稳定处理约
二十步；现有构造也只覆盖后向 Euler 和梯形规则等少数低阶方法，
不易推广到多级 Runge–Kutta。BVM 离散可显著扩大时间窗
（Liu、Wang、Wu 与 Zhou 2022），但积分器范围仍然受限。

ParaDiag-II 近似时间矩阵，并把对角化放入定常迭代或 Krylov 预条件。它牺牲直接性，换取更一般的积分器、更大的时间窗和条件良好的变换。下一页再完整讨论 ParaDiag-II。

## 3.5.1 直接 ParaDiag 方法（ParaDiag-I）

### 后向 Euler 的全时间系统

对 $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$ 使用变步长后向 Euler：

$$
\frac{\boldsymbol u_n-\boldsymbol u_{n-1}}{\Delta t_n}
=A\boldsymbol u_n+\boldsymbol g_n,
\qquad n=1,\ldots,N_t. \tag{3.22}
$$

令 $\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top$，得到

$$
K\boldsymbol U=\boldsymbol b,
\qquad
K=B\otimes I_x-I_t\otimes A, \tag{3.23a}
$$

其中

$$
B=
\begin{bmatrix}
\Delta t_1^{-1}\\
-\Delta t_2^{-1}&\Delta t_2^{-1}\\
&\ddots&\ddots\\
&&-\Delta t_{N_t}^{-1}&\Delta t_{N_t}^{-1}
\end{bmatrix},
\qquad
\boldsymbol b=
\begin{bmatrix}
\Delta t_1^{-1}\boldsymbol u_0+\boldsymbol g_1\\
\boldsymbol g_2\\
\vdots\\
\boldsymbol g_{N_t}
\end{bmatrix}. \tag{3.23b}
$$

若所有 $\Delta t_n$ 互不相同，$B$ 可对角化：

$$
B=VDV^{-1},
\qquad
D=\operatorname{diag}\!\left(
\frac1{\Delta t_1},\ldots,\frac1{\Delta t_{N_t}}
\right). \tag{3.24}
$$

Kronecker 结构给出

$$
K=(V\otimes I_x)
(D\otimes I_x-I_t\otimes A)
(V^{-1}\otimes I_x).
$$

于是求解分成三步：

$$
\left\{
\begin{aligned}
\boldsymbol U^a&=(V^{-1}\otimes I_x)\boldsymbol b,
&&\text{时间逆变换},\\
\left(\frac1{\Delta t_n}I_x-A\right)\boldsymbol u_n^b
&=\boldsymbol u_n^a,
&&n=1,\ldots,N_t,\\
\boldsymbol U&=(V\otimes I_x)\boldsymbol U^b,
&&\text{时间正变换}.
\end{aligned}
\right. \tag{3.25}
$$

第二步的 $N_t$ 个移位空间系统完全独立，是主要计算。首尾两步只做时间方向的稠密矩阵乘法。

![ParaDiag 的时间变换、独立空间求解与逆变换](assets/diagrams/pint/zh/paradiag-three-stage.svg)

任意互异步长都可以调用数值 `eig` 得到 $V$，额外成本通常不大；
但没有闭式特征向量，就无法完整分析舍入误差和参数选择。所以下面的
理论专门采用几何网格。Maday 与 Rønquist（2008）曾以 $\mu=1.2$
测试一维热方程，并获得接近理想的加速。

### 几何网格与两类误差

Maday and Rønquist (2008) 采用 $\Delta t_n=\mu^{n-1}\Delta t_1$，$\mu>1$。由 $\sum_n\Delta t_n=T$，

$$
\Delta t_n=\frac{\mu^{n-1}}{\sum_{j=1}^{N_t}\mu^{j-1}}T. \tag{3.26}
$$

令 $\mu=1+\varrho$。$\varrho$ 过大时，末端步长变大，非均匀网格的截断误差上升；$\varrho$ 过小时，$B$ 接近 Jordan 块，$V$ 的条件数变坏，舍入误差被放大。直接 ParaDiag 的参数选择就是平衡这两项。

### Theorem 3.5：一阶问题的平衡公式

Theorem 3.5 引自 Gander、Halpern、Ryan 与 Tran（2016a）的 Theorems 2 和 6。假设 $\sigma(A)\subset\mathbb R_-$ 且 $|\lambda(A)|\leq\lambda_{\max}$。记 $\boldsymbol u_{N_t}(\varrho)$ 和 $\boldsymbol u_{N_t}(0)$ 为几何网格与均匀网格在 $T$ 的后向 Euler 解，$\widetilde{\boldsymbol u}_n(\varrho)$ 为对角化计算值。则

$$
\left\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u_{N_t}(0)\right\|
\lesssim C(\lambda_*T,N_t)\varrho^2,
$$

$$
\left\|\widetilde{\boldsymbol u}_n(\varrho)-\boldsymbol u_n(\varrho)\right\|
\lesssim
\epsilon\,
\frac{N_t^2(2N_t+1)(N_t+\lambda_{\max}T)}{\phi(N_t)}
\varrho^{-(N_t-1)}. \tag{3.27}
$$

这里

$$
C(x,N_t)=\frac{N_t(N_t^2-1)}{24}r(x/N_t,N_t),
\qquad
r(\widetilde x,N_t)=
\left(\frac{\widetilde x}{1+\widetilde x}\right)^2
(1+\widetilde x)^{-N_t},
$$

$$
\phi(N_t)=
\begin{cases}
\left(\dfrac{N_t}{2}\right)!\left(\dfrac{N_t}{2}-1\right)!,&N_t\ \text{为偶数},\\[6pt]
\left[\left(\dfrac{N_t-1}{2}\right)!\right]^2,&N_t\ \text{为奇数},
\end{cases}
$$

$\widetilde x_*$ 是 $r(\widetilde x,N_t)$ 在 $[0,\infty)$ 上的最大点，$\lambda_*=N_t\widetilde x_*/T$。平衡 (3.27) 的两项得到

$$
\varrho_{\mathrm{opt}}=
\left(
\epsilon
\frac{N_t^2(2N_t+1)(N_t+\lambda_{\max}T)}
{\phi(N_t)C(\lambda_*T,N_t)}
\right)^{1/(N_t+1)}. \tag{3.28}
$$

第一条界只比较几何与均匀网格。相对于精确解的总截断误差还要分成

$$
\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u(T)\|
\le
\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u_{N_t}(0)\|
+\|\boldsymbol u_{N_t}(0)-\boldsymbol u(T)\|.
$$

第二项是熟知的均匀网格后向 Euler 误差，通常不占主导。证明把每个
空间特征值化为 Dahlquist 方程 $y'=\lambda y$；(3.27) 的第二条界
估计特征向量变换引入的舍入误差，并在
$|\lambda|=\lambda_{\max}$ 达到最坏值。几何网格下

$$
V=\mathbb T(p_1,\ldots,p_{N_t-1}),
\qquad
p_n=\frac1{\prod_{j=1}^{n}(1-\mu^j)},
$$

$$
V^{-1}=\mathbb T(q_1,\ldots,q_{N_t-1}),
\qquad
q_n=(-1)^n\mu^{n(n-1)/2}p_n, \tag{3.29a}
$$

> [!warning] 原文公式核对：(3.29a) 的网格比
> 期刊版与 arXiv 版把这里的 $\mu$ 误排成 $\varrho$。由
> $BV=VD$ 的 $(2,1)$ 元素直接得到
> $p_1=1/(1-\mu)$；若写成 $1/(1-\varrho)$，在
> $\varrho\to0$ 时反而保持有界，也与 (3.27) 的
> $\varrho^{-(N_t-1)}$ 放大矛盾。

其中下三角 Toeplitz 算子为

$$
\mathbb T(a_1,\ldots,a_{N_t-1})=
\begin{bmatrix}
1\\
a_1&1\\
\vdots&\ddots&\ddots\\
a_{N_t-1}&\cdots&a_1&1
\end{bmatrix}. \tag{3.29b}
$$

封闭形式用于分析 $\operatorname{Cond}(V)$；实际实现用数值 `eig`，因为它会缩放特征向量来改善条件数。单精度与双精度机器精度分别约为 $1.19\times10^{-7}$ 和 $2.22\times10^{-16}$。

### Figures 3.9–3.10：最优参数仍有并行宽度上限

![原论文 Figure 3.9：热方程与对流扩散方程上误差随几何参数变化](assets/papers/time-parallelization/source-figures/figure-3-9.svg)

Figure 3.9 使用齐次 Dirichlet 边界、$u_0(x)=\sin(2\pi x)$、$\Delta x=1/50$、$T=0.2$。右图的对流扩散黏性为 $10^{-2}$。五组 $N_t$ 分别扫描 $\varrho\in[10^{-2},1]$，误差取所有时间节点上的最大 $L^\infty$ 误差。每条曲线都有一个最小点，星号是 (3.28) 的理论预测。预测对对流扩散很准，对小 $N_t$ 的热方程略有偏差。

![原论文 Figure 3.10：最优几何参数下误差随时间步数先降后升](assets/papers/time-parallelization/source-figures/figure-3-10.svg)

Figure 3.10 改取 $T=0.5$ 和 $N_t=2^4,2^5,\ldots,2^{10}$。均匀后向 Euler 的误差继续随 $N_t$ 下降；ParaDiag-I 使用数值最优 $\varrho_{\mathrm{num}}$ 后，误差先降，在不足 100 步处越过阈值并迅速上升。舍入误差最终压过时间离散误差。

### 波动方程与梯形规则

二阶系统

$$
\boldsymbol u''(t)=A\boldsymbol u(t),
\quad
\boldsymbol u(0)=\boldsymbol u_0,
\quad
\boldsymbol u'(0)=\widetilde{\boldsymbol u}_0 \tag{3.30}
$$

改写为

$$
\boldsymbol w'=\mathbb A\boldsymbol w,
\qquad
\boldsymbol w=(\boldsymbol u^\top,(\boldsymbol u')^\top)^\top,
\qquad
\mathbb A=
\begin{bmatrix}0&I_x\\A&0\end{bmatrix}. \tag{3.31}
$$

为减小波动色散，时间离散采用梯形规则

$$
\frac{\boldsymbol w_n-\boldsymbol w_{n-1}}{\Delta t_n}
=\frac{\mathbb A}{2}(\boldsymbol w_n+\boldsymbol w_{n-1}). \tag{3.32}
$$

当 $A=A^\top\preceq0$ 时，梯形规则保持线性波动系统的二次能量。对
$\boldsymbol w=(\boldsymbol u,\boldsymbol u')$，未缩放状态的物理不变量是

$$
\|\boldsymbol u_n'\|_2^2-\boldsymbol u_n^\top A\boldsymbol u_n,
$$

而不是一般意义下的 $\|\boldsymbol w_n\|_2$；把位移分量按
$(-A)^{1/2}$ 缩放后，才可写成普通 Euclidean 范数守恒。

全时间系统为

$$
K\boldsymbol W=\boldsymbol b,
\qquad
K=B\otimes I_{2N_x}-\widetilde B\otimes\mathbb A, \tag{3.33a}
$$

$$
\widetilde B=\frac12
\begin{bmatrix}
1\\
1&1\\
&\ddots&\ddots\\
&&1&1
\end{bmatrix}. \tag{3.33b}
$$

左乘 $\widetilde B^{-1}\otimes I_{2N_x}$：

$$
\mathcal K\boldsymbol W=\widetilde{\boldsymbol b},
\qquad
\mathcal K=\widetilde B^{-1}B\otimes I_{2N_x}-I_t\otimes\mathbb A,
\qquad
\widetilde{\boldsymbol b}=(\widetilde B^{-1}\otimes I_{2N_x})\boldsymbol b. \tag{3.34}
$$

时间矩阵满足

$$
\widetilde B^{-1}B
=V\operatorname{diag}\!\left(
\frac2{\Delta t_1},\ldots,\frac2{\Delta t_{N_t}}
\right)V^{-1}, \tag{3.35a}
$$

$$
\begin{aligned}
V&=\mathbb T(p_1,\ldots,p_{N_t-1}),
&p_n&=\prod_{j=1}^{n}\frac{1+\mu^j}{1-\mu^j},\\
V^{-1}&=\mathbb T(q_1,\ldots,q_{N_t-1}),
&q_n&=\mu^{-n}\prod_{j=1}^{n}
\frac{1+\mu^{-j+2}}{1-\mu^{-j}}.
\end{aligned} \tag{3.35b}
$$

### Theorem 3.6：波动问题的平衡公式

Theorem 3.6 引自 Gander、Halpern、Rannou 与 Ryan（2019）的 Theorems 2.1 和 2.11。对 $\lambda(A)\leq0$，几何网格与均匀网格的梯形解、以及对角化计算值满足

$$
\left\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u_{N_t}(0)\right\|
\lesssim\frac{N_t(N_t^2-1)}{15}\varrho^2,
$$

$$
\left\|\widetilde{\boldsymbol u}_n(\varrho)-\boldsymbol u_n(\varrho)\right\|
\lesssim
\epsilon\frac{2^{2N_t-1/2}N_t}{(N_t-1)!}
\varrho^{-(N_t-1)}. \tag{3.36}
$$

平衡两项给出

$$
\varrho_{\mathrm{opt}}=
\left(
\epsilon\frac{15\times2^{2N_t-1/2}}
{(N_t^2-1)(N_t-1)!}
\right)^{1/(N_t+1)}. \tag{3.37}
$$

证明对任意 $\lambda\in\sigma(-A)$、$\lambda>0$ 考察
$u''+\lambda u=0$，并令

$$
s=\frac{\lambda T}{2N_t}.
$$

几何和均匀网格的差先得到更细的中间界

$$
O\!\left(
\frac{N_t(N_t^2-1)}6r_1(s)\varrho^2
\right),
\qquad
r_1(s)=\frac{s^3}{(1+s^2)^2}\le\frac25.
$$

舍入项则为

$$
O\!\left(
\epsilon
\frac{2^{2N_t-1/2}N_t}{(N_t-1)!}
r_2(s)\varrho^{-(N_t-1)}
\right),
\qquad
r_2(s)=\frac1{1+s^2}\le1.
$$

对 $r_1,r_2$ 取一致上界便得到 (3.36)。

![原论文 Figure 3.11：波动方程上几何梯形 ParaDiag-I 的最优参数与步数阈值](assets/papers/time-parallelization/source-figures/figure-3-11.svg)

Figure 3.11 使用齐次 Dirichlet 边界、$\Delta x=1/20$、$T=0.2$。(a) 中每组 $N_t$ 都有误差最小点，(3.37) 的星号接近实测最优值；(b) 使用数值最优参数后，$N_t>32$ 时误差迅速恶化。

![原论文 Table 3.1：后向 Euler 和梯形规则时间矩阵的特征向量条件数](assets/papers/time-parallelization/source-figures/table-3-1.svg)

Table 3.1 表明，在 $N_t=5,10,20,30,60,100$ 时，后向 Euler 的 $\operatorname{Cond}(V)$ 从 $1.7\times10^3$ 增至 $4.8\times10^6$；梯形规则从 $4.7\times10^3$ 增至 $4.1\times10^9$。数值最优参数使后段出现平台，但整体增长仍解释了 Figures 3.10–3.11 的舍入恶化。

### BVM：固定步长并改变末步公式

为了扩大并行宽度，Liu et al. (2022) 使用统一 $\Delta t$，前 $N_t-1$ 步采用中心公式，最后一步采用后向 Euler：

$$
\left\{
\begin{aligned}
\frac{\boldsymbol u_{n+1}-\boldsymbol u_{n-1}}{2\Delta t}
&=A\boldsymbol u_n+\boldsymbol g_n,
&&n=1,\ldots,N_t-1,\\
\frac{\boldsymbol u_{N_t}-\boldsymbol u_{N_t-1}}{\Delta t}
&=A\boldsymbol u_{N_t}+\boldsymbol g_{N_t}.
\end{aligned}
\right. \tag{3.38}
$$

它是边值方法（BVM）：只需给定 $\boldsymbol u_0$，所有时间未知量
一次同时求出，稳定性不能按普通逐步中心格式判断。Axelsson 与
Verwer（1985）以这类技术绕开收敛性与稳定性之间的 Dahlquist
屏障，并证明一般非线性情形可取得一致二阶精度，尽管末步本身只有
一阶。Fox（1954）与 Fox 和 Mitchell（1957）更早采用过 BDF2 末步：

$$
\frac{
3\boldsymbol u_{N_t}-4\boldsymbol u_{N_t-1}+\boldsymbol u_{N_t-2}
}{2\Delta t}
=A\boldsymbol u_{N_t}+\boldsymbol g_{N_t}.
$$

全时间形式仍为

$$
K\boldsymbol U=\boldsymbol b,
\qquad
K=B\otimes I_x-I_t\otimes A, \tag{3.39a}
$$

$$
B=\frac1{\Delta t}
\begin{bmatrix}
0&\tfrac12\\
-\tfrac12&0&\tfrac12\\
&\ddots&\ddots&\ddots\\
&&-\tfrac12&0&\tfrac12\\
&&&-1&1
\end{bmatrix},
\qquad
\boldsymbol b=
\begin{bmatrix}
\boldsymbol u_0/(2\Delta t)+\boldsymbol g_1\\
\boldsymbol g_2\\
\vdots\\
\boldsymbol g_{N_t}
\end{bmatrix}. \tag{3.39b}
$$

**Theorem 3.7.** $B=VDV^{-1}$，且 $\operatorname{Cond}(V)=O(N_t^2)$。$V$、$V^{-1}$ 和 $D$ 的闭式表达见 Liu 等（2022, Section 3）。

> [!note] 本站补充：与几何网格的对比
> 原文没有把 Theorem 3.7 与 Theorem 3.5/3.6 的界并排比较。
> 从两组结果看，几何网格的舍入放大因子含
> $\varrho^{-(N_t-1)}$ 和 $2^{2N_t-1/2}/(N_t-1)!$，随
> $N_t$ 增长很快，而 BVM 构造给出 $O(N_t^2)$ 的多项式增长。
> 但 Table 3.1 中按 $\varrho_{\mathrm{num}}$ 实测的条件数最终趋于
> 平台，这一现象仍待研究，所以不能把理论界直接当作实测行为。

二阶系统可先写成一阶系统，并对 $\boldsymbol w$ 使用同一 BVM：

$$
\left\{
\begin{aligned}
\frac{\boldsymbol w_{n+1}-\boldsymbol w_{n-1}}{2\Delta t}
&=\mathbb A\boldsymbol w_n,
&&n=1,\ldots,N_t-1,\\
\frac{\boldsymbol w_{N_t}-\boldsymbol w_{N_t-1}}{\Delta t}
&=\mathbb A\boldsymbol w_{N_t}.
\end{aligned}
\right. \tag{3.40}
$$

为避免速度变量使存储翻倍，可以消去 $\boldsymbol V=(\boldsymbol v_1^\top,\ldots)^\top$。离散关系为

$$
(B\otimes I_x)\boldsymbol U-\boldsymbol V=\boldsymbol b_1,
\qquad
(B\otimes I_x)\boldsymbol V-A\boldsymbol U=\boldsymbol b_2.
$$

代入 $\boldsymbol V=(B\otimes I_x)\boldsymbol U-\boldsymbol b_1$，得到

$$
(B^2\otimes I_x-I_t\otimes A)\boldsymbol U=\boldsymbol b, \tag{3.41}
$$

$$
\boldsymbol b=
\left(
\frac{\widetilde{\boldsymbol u}_0^\top}{2\Delta t},
-\frac{\boldsymbol u_0^\top}{4\Delta t^2},
0,\ldots,0
\right)^\top.
$$

![原论文 Figure 3.12：波动方程上几何时间网格与 BVM 的误差和条件数](assets/papers/time-parallelization/source-figures/figure-3-12.svg)

Figure 3.12 使用 $T=0.5$、$\Delta x=1/40$ 和齐次 Dirichlet 边界。几何梯形 ParaDiag-I 在 $N_t\approx32$ 开始受舍入影响；BVM 保持 $O(\Delta t^2)$，与串行梯形规则一致。右图显示 BVM 的特征向量条件数低得多。

### 非线性全时间方程与准 Newton

对 $\boldsymbol u'=\boldsymbol f(\boldsymbol u,t)$，定义

$$
F(\boldsymbol U)=
\left(
\boldsymbol f(\boldsymbol u_1,t_1)^\top,
\ldots,
\boldsymbol f(\boldsymbol u_{N_t},t_{N_t})^\top
\right)^\top.
$$

非线性全时间方程为

$$
(B\otimes I_x)\boldsymbol U-F(\boldsymbol U)=\boldsymbol b. \tag{3.42}
$$

这里的 $B$ 可以取变步长矩阵 (3.23b)，也可以取 BVM 矩阵
(3.39b)；非线性二阶系统可按同样思路处理。

精确 Newton 更新整理成

$$
\left(B\otimes I_x-\nabla F(\boldsymbol U^k)\right)\boldsymbol U^{k+1}
=\boldsymbol b-\left(
\nabla F(\boldsymbol U^k)\boldsymbol U^k-F(\boldsymbol U^k)
\right), \tag{3.43a}
$$

$$
\nabla F(\boldsymbol U^k)=
\operatorname{blkdiag}\left(
\nabla f(\boldsymbol u_1^k,t_1),\ldots,
\nabla f(\boldsymbol u_{N_t}^k,t_{N_t})
\right). \tag{3.43b}
$$

时间变化的 Jacobian 块会破坏 Kronecker 可分离性。借鉴 Gander 与
Halpern（2017）的思路，可用单一平均矩阵恢复近似结构：

$$
A_k=\frac1{N_t}\sum_{n=1}^{N_t}\nabla f(\boldsymbol u_n^k,t_n),
\qquad\text{或}\qquad
A_k=\nabla f\!\left(
\frac1{N_t}\sum_{n=1}^{N_t}\boldsymbol u_n^k,
\frac{T}{N_t}
\right). \tag{3.44}
$$

于是

$$
(B\otimes I_x-I_t\otimes A_k)\boldsymbol U^{k+1}
=\boldsymbol b-\left((I_t\otimes A_k)\boldsymbol U^k-F(\boldsymbol U^k)\right). \tag{3.45}
$$

对角化 $B=VDV^{-1}$ 后，每轮仍按三步求解：

$$
\left\{
\begin{aligned}
\boldsymbol U^a&=(V^{-1}\otimes I_x)\boldsymbol r^k,\\
(\lambda_n I_x-A_k)\boldsymbol u_n^b&=\boldsymbol u_n^a,
&&n=1,\ldots,N_t,\\
\boldsymbol U^{k+1}&=(V\otimes I_x)\boldsymbol U^b,
\end{aligned}
\right. \tag{3.46}
$$

其中 $\boldsymbol r^k$ 是 (3.45) 的右端。线性情形 $A_k=A$，该式退化为 (3.25)。若 Jacobian 沿时间变化很大，一个 $A_k$ 无法同时逼近全部块，应缩短时间窗并按窗串行。

![原论文 Figure 3.13：BVM ParaDiag-I 在两组黏性和多个时间窗长度下的 Burgers 收敛](assets/papers/time-parallelization/source-figures/figure-3-13.svg)

Figure 3.13 使用周期 Burgers 方程、$\Delta x=0.01$，并保持 $N_t=T/\Delta t=200$。横线为 $\max\{\Delta t^2,\Delta x^2\}=10^{-4}$。$\nu=0.1$ 时，$T=0.1$ 到 $1.6$ 的收敛变化较小；$\nu=0.002$ 时，时间窗变长会明显恶化，$T=0.8,1.6$ 未收敛到目标线。

![原论文 Table 3.2：串行梯形规则与并行 BVM ParaDiag-I 的 Jacobian 求解次数](assets/papers/time-parallelization/source-figures/table-3-2.svg)

若有 $N_t$ 个处理器，每轮的 $N_t$ 个 Jacobian 系统同时求解，
并行 Jacobian 求解次数就等于外层轮数；串行梯形规则的计数则是
$\sum_{n=1}^{N_t}\mathrm{It}_n$。Table 3.2 中，$\nu=0.1$ 时
串行梯形需要 401–443 次，ParaDiag-I 只需 5–7 轮，说明 Jacobian
沿轨道变化较小，单个 $A_k$ 足以代表全部块。$\nu=0.002$ 时，
串行计数为 400/446/476/460/526，ParaDiag-I 从 7 增至 12、22 轮，
并在 $T=0.8,1.6$ 失效；长窗下的平均 Jacobian 已失去代表性。

### 最近 Kronecker 近似

这一加速来自 Liu 与 Wu（2022, Section 3.3）。单一 $I_t\otimes A_k$ 忽略了 Jacobian 随时间的幅值变化。改用 $\Phi_k\otimes A_k$，其中 $\Phi_k=\operatorname{diag}(\phi_1,\ldots,\phi_{N_t})$，解

$$
\min_{\Phi_k\ \mathrm{diagonal}}
\left\|\nabla F(\boldsymbol U^k)-\Phi_k\otimes A_k\right\|. \tag{3.47}
$$

在 $\operatorname{trace}(A_k^\top A_k)>0$ 的前提下，Frobenius 范数下的最近 Kronecker 近似有闭式（Van Loan 与 Pitsianis 1993, Theorem 3）：

$$
\phi_n=
\frac{\operatorname{trace}\!\left(
\nabla f(\boldsymbol u_n^k,t_n)A_k^\top
\right)}
{\operatorname{trace}(A_k^\top A_k)},
\qquad n=1,\ldots,N_t. \tag{3.48}
$$

新准 Newton 方程为

$$
(B\otimes I_x-\Phi_k\otimes A_k)\boldsymbol U^{k+1}
=\boldsymbol b-\left((\Phi_k\otimes A_k)\boldsymbol U^k-F(\boldsymbol U^k)\right).
$$

左乘 $B^{-1}\otimes I_x$ 后，只需对角化 $B^{-1}\Phi_k$；第二步变成 $(I_x-\lambda_nA_k)\boldsymbol u_n^b=\boldsymbol u_n^a$。尚无一般对角化理论，但实验中该矩阵通常可对角化且特征向量条件良好。

$\phi_n$ 的矩阵乘法代价较高，实用上可在粗空间模型上离线计算一次。
Figure 3.14 的精细网格为 $\Delta x=1/200$，缩放因子由
$\Delta X=1/20$ 的粗模型得到。

![原论文 Figure 3.14：平均 Jacobian 与 NKA 准 Newton 的 Burgers 收敛比较](assets/papers/time-parallelization/source-figures/figure-3-14.svg)

两幅图分别对应两种黏性，并比较 $T=0.7$ 和 $T=1.3$。NKA 在所有设置下更快，长窗口 $T=1.3$ 的收益最明显。

## 公式与图表覆盖核对

| 原文项目                                           | 论文小节 | 覆盖状态                                                  |
| -------------------------------------------------- | -------- | --------------------------------------------------------- |
| (3.22)–(3.25)                                      | 3.5.1    | 变步长后向 Euler、全时间矩阵、对角化和三步求解            |
| (3.26)–(3.29), Theorem 3.5                         | 3.5.1    | 几何网格、两类误差、平衡参数、Toeplitz 特征向量与证明路线 |
| Figures 3.9–3.10                                   | 3.5.1    | 全部原图、参数扫描和并行宽度阈值                          |
| (3.30)–(3.37), Theorem 3.6, Figure 3.11, Table 3.1 | 3.5.1    | 二阶转写、梯形全时间系统、误差平衡、证明界和条件数        |
| (3.38)–(3.41), Theorem 3.7, Figure 3.12            | 3.5.1    | BVM、统一二阶精度、条件数和消元推导                       |
| (3.42)–(3.46), Figure 3.13, Table 3.2              | 3.5.1    | 非线性 Newton、平均 Jacobian、并行求解与成本对比          |
| (3.47)–(3.48), Figure 3.14                         | 3.5.1    | NKA、缩放矩阵、离线粗模型和收敛改善                       |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 3.5–3.5.1, pp. 415–430.
