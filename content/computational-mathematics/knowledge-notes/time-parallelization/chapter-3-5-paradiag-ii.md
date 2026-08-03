---
title: 3.5.2：迭代 ParaDiag（ParaDiag-II）
description: 从循环预条件、FFT 同时对角化、首尾耦合波形松弛，到 Krylov 与定常迭代的区别、谱界、α 的收敛–舍入权衡和非线性 Newton–Krylov 的完整推导
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-5-paradiag-ii
tags:
  - 时间并行
  - ParaDiag
  - 循环预条件
---

> [!note] 阅读范围
> 本页对应论文 Section 3.5.2（pp. 431–443），覆盖公式 (3.49)–(3.68)、Theorems 3.8–3.9 和 Figures 3.15–3.18。文中从离散循环预条件和连续首尾耦合两条路线出发，随后证明两条路线落到同一个 ParaDiag-II 迭代。

## 3.5.2 迭代 ParaDiag 方法（ParaDiag-II）

### 从多步公式到循环预条件

ParaDiag-I 很难直接推广到多级 Runge–Kutta 等高阶积分器。ParaDiag-II 改为近似时间步进矩阵，并在定常迭代或 Krylov 方法中求解全时间系统。McDonald et al. (2018) 从离散系统出发；Gander and Wu (2019) 同期从连续波形松弛出发。两种表述在离散后等价。

先用 $m$ 步线性多步法离散 (2.1)：

$$
\sum_{l=0}^{m}a_l\boldsymbol u_{n-l}
=\Delta t\sum_{l=0}^{m}b_lA\boldsymbol u_{n-l}
+\bar{\boldsymbol g}_n,
\qquad n=1,\ldots,N_t.
$$

已知初始历史 $\{\boldsymbol u_{-(m-1)},\ldots,\boldsymbol u_0\}$。令

$$
\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top,
\qquad
K=B_1\otimes I_x-B_2\otimes(\Delta t A),
$$

其中 $B_1,B_2$ 是由 $\{a_l\}$ 和 $\{b_l\}$ 生成的下三角带状 Toeplitz 矩阵。把两者换成 Strang 型循环矩阵 $C_1,C_2$，得到

$$
P=C_1\otimes I_x-C_2\otimes(\Delta t A).
$$

若 $P$ 已足够接近 $K$，可直接使用定常迭代

$$
P\Delta\boldsymbol U^k=\boldsymbol r^k
:=\boldsymbol b-K\boldsymbol U^k,
\qquad
\boldsymbol U^{k+1}=\boldsymbol U^k+\Delta\boldsymbol U^k,
\qquad k=0,1,\ldots. \tag{3.49}
$$

其渐近速度由 $\rho(P^{-1}K)$ 控制。即使谱半径不小，$P^{-1}K$ 的谱若高度聚集，预条件 GMRES 仍可能很快。

> [!tip] 本站洞见：定常迭代看谱半径，Krylov 看谱的聚集
> 定常迭代 (3.49) 的误差按 $\mathcal M=I-P^{-1}K$ 的幂衰减，第 $k$ 步误差范数约为 $\rho(\mathcal M)^k$（在 $\mathcal M$ 可对角化且特征向量条件数适中时）。因此只要 $\rho(\mathcal M)=\rho(I-P^{-1}K)\ge1$——即 $P^{-1}K$ 有特征值远离 $1$——定常迭代就不再收敛。GMRES 的行为完全不同：其第 $k$ 步相对残差满足多项式极小化界
> $$\frac{\|\boldsymbol r^k\|}{\|\boldsymbol r^0\|}\le\kappa(V)\min_{\substack{p\in\Pi_k,\;p(0)=1}}\ \max_{\lambda\in\sigma(P^{-1}K)}|p(\lambda)|,$$
> 其中 $V$ 是 $P^{-1}K$ 的特征向量矩阵、$\kappa(V)$ 是其条件数。这个上界只问一件事：能否用一个低次首项归一多项式在整个谱上取小值。若谱聚成少数几团（尤其绝大多数特征值精确等于 $1$），只需次数约等于"离群团数"的多项式就能让 $\max|p(\lambda)|$ 骤降，与 $\rho$ 是否 $\ge1$ 无关。这正是下面 Theorem 3.8 的意义所在：即便 $\rho(\mathcal M)\ge1$ 使定常迭代失效，只要非 $1$ 特征值个数少，预条件 GMRES 依旧在很少步内收敛。

### Fourier 对角化与三步并行求解

循环矩阵彼此交换，可由同一个离散 Fourier 矩阵同时对角化：

$$
C_l=F^*D_lF,\qquad l=1,2,
$$

$$
F=\frac1{\sqrt{N_t}}
\begin{bmatrix}
1&1&\cdots&1\\
1&\omega&\cdots&\omega^{N_t-1}\\
\vdots&\vdots&&\vdots\\
1&\omega^{N_t-1}&\cdots&\omega^{(N_t-1)^2}
\end{bmatrix},
\qquad
\omega=\exp\!\left(\frac{2\pi i}{N_t}\right). \tag{3.50}
$$

若 $C_l(:,1)$ 表示第一列，则

$$
D_l=\operatorname{diag}\!\left(\sqrt{N_t},F C_l(:,1)\right),
\qquad l=1,2. \tag{3.51}
$$

于是

$$
P=(F^*\otimes I_x)
\left(D_1\otimes I_x-D_2\otimes(\Delta t A)\right)
(F\otimes I_x),
$$

一次预条件作用分为

$$
\left\{
\begin{aligned}
\boldsymbol U^a&=(F\otimes I_x)\boldsymbol r^k,\\
(\lambda_{1,n}I_x-\lambda_{2,n}\Delta t A)\boldsymbol u_n^b
&=\boldsymbol u_n^a,
&&n=1,\ldots,N_t,\\
\boldsymbol U&=(F^*\otimes I_x)\boldsymbol U^b.
\end{aligned}
\right. \tag{3.52}
$$

首尾两步由 FFT 完成，复杂度为 $O(N_xN_t\log N_t)$；中间的 $N_t$ 个空间系统相互独立。

> [!tip] 本站洞见：为什么循环矩阵能"同时"对角化，成本又落到 $O(N_xN_t\log N_t)$
> 关键在于所有循环矩阵共享同一组特征向量——离散 Fourier 矩阵 $F$ 的列。任意 $N_t\times N_t$ 循环矩阵都是循环位移算子 $S$ 的多项式，而 $S=F^*\Omega F$，$\Omega=\operatorname{diag}(1,\omega,\ldots,\omega^{N_t-1})$。于是 $C_1=p_1(S)=F^*p_1(\Omega)F$、$C_2=p_2(S)=F^*p_2(\Omega)F$ 被同一个 $F$ 对角化，两者自动交换 $C_1C_2=C_2C_1$。在 Kronecker 结构里，$P=C_1\otimes I_x-C_2\otimes(\Delta tA)$ 的时间因子被同一个 $F\otimes I_x$ 块对角化，而空间矩阵 $A$ 根本不参与这一步——它只出现在被解耦出来的 $N_t$ 个空间系统 $\lambda_{1,n}I_x-\lambda_{2,n}\Delta tA$ 里。成本因此分三段：(step a)、(step c) 是沿时间方向对每个空间自由度做一次长度 $N_t$ 的 FFT，共 $N_x$ 条，代价 $O(N_xN_t\log N_t)$；(step b) 是 $N_t$ 个规模 $N_x$ 的独立空间解，可完全并行。这与直接求解在时间上稠密耦合的全时间系统 $K$ 相比，把耦合"对角化掉"正是 ParaDiag（parallel-in-time by diagonalization）名字的由来。反过来，如果两个时间矩阵不交换（例如换成一般非循环的 Toeplitz），就没有公共特征向量，这套同时对角化随即失效——这也是必须用循环（或 $\alpha$-循环）矩阵去近似 Toeplitz 矩阵的根本原因。

### Theorem 3.8：精确聚集的上界

若 $A\in\mathbb R^{N_x\times N_x}$ 对称负定，则 $P^{-1}K$ 至多有 $mN_x$ 个特征值不等于 $1$。因此精确算术下，预条件 GMRES 至多 $mN_x+1$ 步终止。该结论对大 $N_x$ 并未给出很小的迭代上界；$A$ 失去对称性后，谱聚集也会变弱。

> [!tip] 本站洞见：$mN_x$ 从哪来，为什么大 $N_x$ 会稀释它、非对称 $A$ 会破坏它
> Strang 循环矩阵 $C_l$ 与 Toeplitz 矩阵 $B_l$ 只在少数几条"环绕"对角线（右上角的循环补项）上不同：$C_l-B_l$ 的非零元集中在角上，秩至多为 $m$（$m$ 步公式的带宽）。于是
> $$P-K=(C_1-B_1)\otimes I_x-(C_2-B_2)\otimes(\Delta tA)$$
> 的秩至多 $mN_x$，从而 $P^{-1}K=I-P^{-1}(P-K)$ 是单位阵加一个秩 $\le mN_x$ 的扰动，最多 $mN_x$ 个特征值偏离 $1$。把它代入上一节的 GMRES 多项式界：取一个在这 $\le mN_x$ 个离群值处为零、在 $\lambda=1$ 处为 $1$ 的多项式（次数 $\le mN_x$），即得精确算术下至多 $mN_x+1$ 步终止。要点是这个上界随空间自由度 $N_x$ **线性增长**：细化空间会成比例抬高最坏步数，因此它是"有限步终止"的定性保证，而非"网格无关"的快速收敛保证。$A$ 对称负定时空间特征值落在负实轴，$P^{-1}K$ 的离群特征值排布规整、易被低次多项式统一压制；一旦 $A$ 非对称（强对流、波动），离群特征值散入复平面、聚集变松，实际 GMRES 步数明显上升——这正是接下来 Figure 3.15 从热方程到波动方程的连续退化。

### Figure 3.15：三类 PDE 上的谱与 GMRES

论文比较热方程、两种黏性的对流扩散方程，以及波动方程。全部使用齐次 Dirichlet 边界和 $u(x,0)=\sin(2\pi x)$；波动方程另取 $u_t(x,0)=0$。空间中心差分后，波动方程为

$$
\boldsymbol u''(t)=A\boldsymbol u(t),
\qquad
\boldsymbol u(0)=\boldsymbol u_0,
\quad \boldsymbol u'(0)=0,
\quad t\in(0,T], \tag{3.53}
$$

其中 $A=\operatorname{Tri}[1,-2,1]/\Delta x^2$。一阶问题使用梯形规则，二阶问题使用带参数的 Numerov 型方法：

$$
\left\{
\begin{aligned}
\widetilde{\boldsymbol u}_n-\boldsymbol u_n
+\gamma\Delta t^2A
(\boldsymbol u_{n+1}-2\boldsymbol u_n+\boldsymbol u_{n-1})&=0,\\
\boldsymbol u_{n+1}-2\boldsymbol u_n+\boldsymbol u_{n-1}
-\frac{\Delta t^2A}{12}
(\boldsymbol u_{n+1}+10\widetilde{\boldsymbol u}_n+\boldsymbol u_n)&=0.
\end{aligned}
\right. \tag{3.54}
$$

$\gamma=0$ 给出经典 Numerov 四阶方法，但仅条件稳定；$\gamma\ge1/120$ 时保持四阶并达到无条件稳定。

> [!tip] 本站洞见：参数 $\gamma$ 如何把条件稳定"抬"成无条件稳定，及稳定性为何是谱界的必要条件
> 经典 Numerov（$\gamma=0$）是紧致四阶格式，但绝对稳定区间有限：对波动这类特征值贴近虚轴（把 $A$ 取为 Laplacian 时贴近负实轴）的问题，$\Delta t$ 必须足够小才能让放大因子 $|r_1^{-1}(z)r_2(z)|$ 不越界。(3.54) 中的 $\gamma\Delta t^2A(\boldsymbol u_{n+1}-2\boldsymbol u_n+\boldsymbol u_{n-1})$ 项相当于加一个可调的隐式"人工质量"，把 $r_1,r_2$ 的高阶项配平：当 $\gamma\ge1/120$ 时，对全部 $z=\Delta t^2\lambda\in\mathbb R_-$ 都有 $|r_1^{-1}(z)r_2(z)|\le2$ 且仅在 $z=0$ 取等号（恰是 Theorem 3.9 二阶情形所需的稳定条件），从而无条件稳定；同时因为只添加了 $O(\Delta t^2)$ 量级的高阶修正，四阶精度不受影响。这个 $1/120$ 阈值取自 Chawla (1983) 对 Numerov 型方法的稳定性分析（补充参考）。稳定性之所以是谱界的**必要**（而非仅充分）条件：Theorem 3.9 的界靠 $|r_1^{-1}r_2|\le1$（一阶）或 $\le2$（二阶）把每个解耦子块的放大因子压住；一旦格式失稳，某些 $z$ 上 $|r_1^{-1}r_2|$ 越界，对应子块的特征值就逃出半径 $\alpha/(1-\alpha)$ 的圆——Figure 3.17 下排的 $\gamma=1/120.01$ 就是这种越界的直接可视化。

![原论文 Figure 3.15：三类 PDE 上循环预条件后的谱和 GMRES 残差](assets/papers/time-parallelization/source-figures/figure-3-15.svg)

Figure 3.15 使用 $T=2$、$\Delta t=1/50$、$\Delta x=1/100$、$\gamma=1/100$。左列 (a)、(c)、(e) 画 $P^{-1}K$ 的谱，右列 (b)、(d)、(f) 画对应的 GMRES 残差；三行依次是热方程、对流扩散方程和波动方程。(a)、(b) 中热方程的特征值几乎全聚在 $1$，残差三轮左右便降至机器精度。(c)、(d) 同时比较 $\nu=10^{-3}$ 与 $10^{-6}$：前者仍高度聚集，后者的非单位特征值明显散开并需要更多迭代。(e)、(f) 的波动谱沿虚方向大幅展开，GMRES 进入长尾阶段。六个面板共同说明退化过程是连续的：耗散减弱先破坏谱聚集，最终使循环预条件在双曲极限下失去效率。

### 连续首尾耦合与 $\alpha$-循环矩阵

Gander and Wu (2019) 在连续层面引入

$$
\boldsymbol u_t^k(t)=A\boldsymbol u^k(t)+\boldsymbol g(t),
\qquad
\boldsymbol u^k(0)=
\alpha[\boldsymbol u^k(T)-\boldsymbol u^{k-1}(T)]+\boldsymbol u_0, \tag{3.55}
$$

其中 $\alpha\in\mathbb C$。收敛后括号内的差消失，初值回到 $\boldsymbol u_0$。

用一类单步公式离散：

$$
\left\{
\begin{aligned}
r_1(\Delta t A)\boldsymbol u_n^k
&=r_2(\Delta t A)\boldsymbol u_{n-1}^k+\widetilde{\boldsymbol g}_n,
&&n=1,\ldots,N_t,\\
\boldsymbol u_0^k
&=\alpha(\boldsymbol u_{N_t}^k-\boldsymbol u_{N_t}^{k-1})+\boldsymbol u_0.
\end{aligned}
\right. \tag{3.56}
$$

后向 Euler 与梯形规则分别对应

$$
\left\{
\begin{aligned}
r_1&=I_x-\Delta tA,&r_2&=I_x,
&&\text{后向 Euler},\\
r_1&=I_x-\tfrac12\Delta tA,&r_2&=I_x+\tfrac12\Delta tA,
&&\text{梯形规则}.
\end{aligned}
\right. \tag{3.57}
$$

把首尾条件代入第一步，可写成

$$
P_\alpha\boldsymbol U^k=\boldsymbol b^k, \tag{3.58a}
$$

$$
\boldsymbol b^k=\boldsymbol b-
\alpha
\begin{bmatrix}
r_2(\Delta tA)\boldsymbol u_{N_t}^{k-1}\\0\\\vdots\\0
\end{bmatrix},
\qquad
\boldsymbol b=
\begin{bmatrix}
r_2(\Delta tA)\boldsymbol u_0+\widetilde{\boldsymbol g}_1\\
\widetilde{\boldsymbol g}_2\\\vdots\\
\widetilde{\boldsymbol g}_{N_t}
\end{bmatrix}. \tag{3.58b}
$$

$$
P_\alpha=
\begin{bmatrix}
r_1& & &-\alpha r_2\\
-r_2&r_1\\
&\ddots&\ddots\\
&&-r_2&r_1
\end{bmatrix}
=I_t\otimes r_1(\Delta tA)-C_\alpha\otimes r_2(\Delta tA), \tag{3.58c}
$$

$$
C_\alpha=
\begin{bmatrix}
0&&&\alpha\\
1&0\\
&\ddots&\ddots\\
&&1&0
\end{bmatrix}.
$$

对 Strang 型 $\alpha$-循环矩阵，

$$
C_\alpha=V_\alpha D_\alpha V_\alpha^{-1}, \tag{3.59a}
$$

$$
D_\alpha=\operatorname{diag}\!\left(\sqrt{N_t},F\Lambda_\alpha C_\alpha(:,1)\right),
\quad
V_\alpha=\Lambda_\alpha F^*,
\quad
\Lambda_\alpha=\operatorname{diag}
(1,\alpha^{-1/N_t},\ldots,\alpha^{-(N_t-1)/N_t}). \tag{3.59b}
$$

因此求解 (3.58a) 仍是三步：

$$
\left\{
\begin{aligned}
\boldsymbol U^a&=(V_\alpha^{-1}\otimes I_x)\boldsymbol b^k,\\
\left(r_1(\Delta tA)-\lambda_n r_2(\Delta tA)\right)\boldsymbol u_n^b
&=\boldsymbol u_n^a,
&&n=1,\ldots,N_t,\\
\boldsymbol U^k&=(V_\alpha\otimes I_x)\boldsymbol U^b.
\end{aligned}
\right. \tag{3.60}
$$

$\alpha=1$ 时 $V_\alpha=F^*$，该算法与 (3.52) 完全一致。$\alpha\ne1$ 时，$V_\alpha$ 只是 Fourier 变换前后多乘一个对角缩放，FFT 仍可使用。

> [!tip] 本站洞见：$\alpha$ 只是给 FFT 前后套一层对角缩放
> 由 $V_\alpha=\Lambda_\alpha F^*$：(step a) 用 $V_\alpha^{-1}=(\Lambda_\alpha F^*)^{-1}=F\Lambda_\alpha^{-1}$，(step c) 用 $V_\alpha=\Lambda_\alpha F^*$。写进 Kronecker 结构就是
> $$V_\alpha^{-1}\otimes I_x=(F\otimes I_x)(\Lambda_\alpha^{-1}\otimes I_x),\qquad V_\alpha\otimes I_x=(\Lambda_\alpha\otimes I_x)(F^*\otimes I_x).$$
> 也就是：整套变换 = 逐时间点乘 $\alpha^{\pm n/N_t}$ 的对角缩放 + 一次标准 FFT/IFFT。$\alpha=1$ 时 $\Lambda_\alpha=I$，缩放消失，退化为 McDonald et al. (2018) 的纯 FFT 对角化 (3.52)；$0<\alpha<1$ 时只是多两次 $O(N_xN_t)$ 的对角乘法，FFT 的 $O(N_xN_t\log N_t)$ 成本与并行结构分毫不变。这层缩放正是后文舍入分析里 $\operatorname{Cond}_2(V_\alpha)=1/\alpha$ 的来源——收敛速度与数值稳定性都被同一个 $\alpha$ 牵着走。

![原论文 Figure 3.16：两种 alpha 下首尾耦合波形松弛的误差衰减](assets/papers/time-parallelization/source-figures/figure-3-16.svg)

Figure 3.16 沿用 Figure 3.15 的 PDE 参数并改为周期边界。左、右面板分别取 $\alpha=0.1$ 与 $10^{-3}$，每个面板同时比较 $\nu=10^{-2}$、$10^{-6}$ 的对流扩散方程和波动方程。此时 $\alpha=1$ 的 $P$ 奇异，需要取 $\alpha<1$。$\alpha=0.1$ 时三条曲线均线性下降，波动方程最慢；减到 $10^{-3}$ 后，三类问题都在更少迭代内达到 $10^{-12}$ 附近。该图直接展示 $\alpha$ 对收敛斜率的控制，也表明强对流和波动问题无需 Krylov 加速即可收敛。

### 两条路线的离散等价

首尾耦合离散可以写成预条件定常迭代

$$
P_\alpha\Delta\boldsymbol U^{k-1}
=\boldsymbol r^{k-1}:=\boldsymbol b-K\boldsymbol U^{k-1},
\qquad
\boldsymbol U^k=\boldsymbol U^{k-1}+\Delta\boldsymbol U^{k-1}. \tag{3.61}
$$

这里

$$
K=I_t\otimes r_1(\Delta tA)-B\otimes r_2(\Delta tA), \tag{3.62a}
$$

$$
B=
\begin{bmatrix}
0\\1&0\\&\ddots&\ddots\\&&1&0
\end{bmatrix}. \tag{3.62b}
$$

由 $\boldsymbol b^k=(P_\alpha-K)\boldsymbol U^{k-1}+\boldsymbol b$，(3.58a) 立即化为 (3.61)。$\alpha=1$ 时，它就是 McDonald et al. 的 (3.49)；$0<\alpha<1$ 时，对应 Banjai and Peterseim (2012) 的并行方法。

> [!tip] 本站洞见：一张"谁等于谁"的路线图
> 离散层面只有一个对象——全时间方程 $K\boldsymbol U=\boldsymbol b$ 及其块 $\alpha$-循环预条件器 $P_\alpha$——两条历史路线只是它的不同"入口"：
>
> - **离散入口**（McDonald–Pestana–Wathen 2018）：直接把 Toeplitz 时间矩阵换成 Strang 循环矩阵，对应 $\alpha=1$，得到 (3.49)。
> - **连续入口**（Gander & Wu 2019）：在连续层面写首尾耦合波形松弛 (3.55)，离散后自然得到带自由参数 $\alpha$ 的 $\alpha$-循环 $P_\alpha$（(3.58)）。
>
> 代数上，$\boldsymbol b^k=(P_\alpha-K)\boldsymbol U^{k-1}+\boldsymbol b$ 把定值迭代 (3.58a) 改写成误差–预条件形式 (3.61)，两者恒等。于是：$\alpha=1$ 落到 McDonald 等人的循环预条件定常迭代；$0<\alpha<1$ 落到 Banjai and Peterseim (2012) 的并行方法（补充参考）。周期边界下 $\alpha=1$ 的 $P$ 奇异（循环矩阵的零频行退化），必须取 $\alpha<1$——这正是连续入口引入自由参数 $\alpha$ 反而更通用的原因，也解释了 Figure 3.16 为何只画 $\alpha<1$。

一般单步公式

$$
r_1(\Delta tA)\boldsymbol u_n
=r_2(\Delta tA)\boldsymbol u_{n-1}+\widetilde{\boldsymbol g}_n,
\qquad n=1,\ldots,N_t, \tag{3.63}
$$

生成全时间方程 $K\boldsymbol U=\boldsymbol b$，$P_\alpha$ 正是其广义块 $\alpha$-循环预条件器。

### 直接离散二阶系统

把二阶系统改写为一阶会让每个时间点的存储翻倍。论文因此考虑对称两步格式

$$
r_1(\Delta t^2A)\boldsymbol u_{n+1}
-r_2(\Delta t^2A)\boldsymbol u_n
+r_1(\Delta t^2A)\boldsymbol u_{n-1}
=\widetilde{\boldsymbol g}_n,
\quad n=1,\ldots,N_t-1, \tag{3.64}
$$

并假设第二个初值 $\boldsymbol u_1$ 已知。Numerov 型方法对应

$$
r_1(z)=I_x-\frac{z}{12}+\frac{10\gamma z^2}{12},
\qquad
r_2(z)=2I_x+\frac{10z}{12}+\frac{20\gamma z^2}{12},
\quad z=\Delta t^2A.
$$

全时间矩阵及预条件器为

$$
K=\widetilde B\otimes r_1(\Delta t^2A)-B\otimes r_2(\Delta t^2A),
\qquad
P_\alpha=\widetilde C_\alpha\otimes r_1(\Delta t^2A)
-C_\alpha\otimes r_2(\Delta t^2A), \tag{3.65a}
$$

其中 $B,C_\alpha$ 如上，

$$
\widetilde B=
\begin{bmatrix}
1\\0&1\\1&0&1\\&\ddots&\ddots&\ddots\\&&1&0&1
\end{bmatrix},
\qquad
\widetilde C_\alpha=
\begin{bmatrix}
1&&\alpha\\0&1&&\alpha\\1&0&1\\&\ddots&\ddots&\ddots\\&&1&0&1
\end{bmatrix}. \tag{3.65b}
$$

$C_\alpha$ 与 $\widetilde C_\alpha$ 可同时对角化，所以 (3.60) 的三步求解也适用于二阶预条件器。

### Theorem 3.9：稳定性、谱界与收敛因子

对于一阶问题，如果 (3.63) 稳定，即

$$
|r_1^{-1}(z)r_2(z)|\le1,
\qquad z\in\sigma(\Delta tA)\subset\mathbb C_-,
$$

论文印刷版给出

$$
\frac1{1-\alpha}
\le |\lambda(P_\alpha^{-1}K)|
\le\frac1{1+\alpha},
\qquad 0<\alpha<1. \tag{3.66, as printed}
$$

二阶两步格式若满足 $|r_1^{-1}(z)r_2(z)|\le2$，且等号只在 $z=0$ 取得，也声称满足同一界。

> [!warning] 原文公式核对
> 对 $0<\alpha<1$，印刷式左端大于右端，区间为空。这里完整保留论文排版，并把疑似校正记为 $1/(1+\alpha)\le|\lambda|\le1/(1-\alpha)$。后文给出的迭代矩阵界也与“分母交换后的次序”相容。此处是读者核对后的注释，不改写原论文陈述。

这一"任意稳定单步/对称两步格式都成立"的一般谱分析出自 Wu, Zhou and Zhou (2022)（补充参考）；此前针对具体积分器与具体方程类型的分析散见于抛物型（如 Gu and Wu 2020、Lin and Ng 2021 等）和双曲型（如 Liu and Wu 2020 等）的一系列工作。

定常迭代矩阵为

$$
\mathcal M=I-P_\alpha^{-1}K, \tag{3.67}
$$

论文据此得到 $\rho(\mathcal M)\le\alpha/(1-\alpha)$。小 $\alpha$ 加快收敛，Figure 3.16 的实验与此吻合。底层积分器稳定是谱界成立的充分条件；论文的数值结果也表明它接近必要条件。

![原论文 Figure 3.17：稳定阈值两侧的 Numerov 型方法迭代谱](assets/papers/time-parallelization/source-figures/figure-3-17.svg)

Figure 3.17 取 $\Delta t=1/16$、$\Delta x=1/128$、$\alpha=0.02$，并比较 $T=0.5,10,20$。上排 $\gamma=1/120$ 位于无条件稳定阈值，所有特征值都在半径 $\alpha/(1-\alpha)$ 的虚线圆内；下排 $\gamma=1/120.01$ 略低于阈值，长时间窗下谱界失效。

> [!tip] 本站洞见：为什么阈值只差 $0.01$，却要长时间窗才暴露
> $\gamma=1/120$ 恰在 Chawla (1983) 的无条件稳定阈值上，全部 $z=\Delta t^2\lambda\in\mathbb R_-$ 都满足 $|r_1^{-1}(z)r_2(z)|\le2$（仅 $z=0$ 取等号），Theorem 3.9 的二阶条件成立，$\mathcal M$ 的特征值悉数落进半径 $\alpha/(1-\alpha)=0.02/0.98\approx0.0204$ 的圆内，且 $T=0.5,10,20$ 无论多长都保持——谱界与时间窗长度无关。$\gamma=1/120.01$ 只低一点点，却让最"硬"的空间模态（$|z|$ 最大、对应最高频）越出稳定区，其 $|r_1^{-1}r_2|>2$。为什么要大 $T$ 才看得出：$\Delta t$ 固定时 $N_t=T/\Delta t$ 随 $T$ 增大，$\alpha$-循环特征值 $\lambda_n$ 越密集地扫过复平面、越容易"命中"那个失稳子块，被它撑破的特征值也就逐个冒出圆外；$T=0.5$ 时 $N_t$ 太少，可能"漏采"失稳频率，谱界看似仍成立。这既解释了 Figure 3.17 上下两排随 $T$ 的对照，也把"稳定性是必要条件"落到了图上。

### 小 $\alpha$ 的舍入限制与实现方式

$\alpha$ 不能无限减小。由 $V_\alpha=\Lambda_\alpha F^*$、$\operatorname{Cond}_2(F^*)=1$，

$$
\operatorname{Cond}_2(V_\alpha)=\frac1\alpha,
\qquad
\mathrm{err}_{\mathrm{ro}}
=O\!\left(\epsilon\operatorname{Cond}_2(V_\alpha)\right)
=O\!\left(\frac\epsilon\alpha\right).
$$

这一矩阵分解误差并不必然等量传入最终迭代。直接由 (3.58) 求 $\boldsymbol U^k$ 会显露 $O(\epsilon/\alpha)$ 的放大；先由 (3.61) 求误差增量 $\Delta\boldsymbol U^{k-1}$，再更新解，能显著减轻舍入误差。

> [!tip] 本站洞见：$\alpha$ 的双向权衡，与增量形式为何能"吃掉"放大
> 谱侧与舍入侧对 $\alpha$ 的诉求正好相反，合成一个 U 形权衡：
>
> - **收敛侧**：$\rho(\mathcal M)\le\alpha/(1-\alpha)$，$\alpha$ 越小，定常迭代每步压缩越猛（Figure 3.16）。
> - **舍入侧**：$\operatorname{Cond}_2(V_\alpha)=1/\alpha$，对角缩放 $\Lambda_\alpha$ 把末端分量放大约 $\alpha^{-(N_t-1)/N_t}\approx1/\alpha$ 倍，浮点分解误差按 $\mathrm{err}_{\mathrm{ro}}=O(\epsilon/\alpha)$ 进入结果。
>
> 因此存在一个令总误差最小的中间 $\alpha$（经验上 $\alpha\in[10^{-3},10^{-1}]$）。增量形式为什么能缓解：直接解 (3.58) 时 $O(\epsilon/\alpha)$ 的**绝对**误差直接叠加在整解 $\boldsymbol U^k$ 上，一旦逼近该量级就无法再降。改用 (3.61)，每步求的是**误差增量** $\Delta\boldsymbol U^{k-1}$，其范数随迭代收敛而变小；舍入放大作用在越来越小的增量上，产生的绝对扰动 $O(\epsilon\|\Delta\boldsymbol U^{k-1}\|/\alpha)$ 同步收缩，于是最终精度不再被 $\epsilon/\alpha$ 锁死。这与经典迭代精化（iterative refinement）同源：把大量级的解与小量级的修正分开算。Figure 3.18 左右两图正是这一机制的对照。

![原论文 Figure 3.18：三种 alpha 下两种实现方式的舍入误差](assets/papers/time-parallelization/source-figures/figure-3-18.svg)

Figure 3.18 同时比较 $\alpha=10^{-3},10^{-6},10^{-11}$。左图直接求 $\boldsymbol U^k$：三条曲线分别停在约 $10^{-12}$、$10^{-9}$ 和 $10^{-4}$，$\alpha$ 越小，$\epsilon/\alpha$ 放大越早显现。右图先求 $\Delta\boldsymbol U^{k-1}$ 再更新，三个参数最终都降到约 $10^{-13}$。这组面板把“算法的谱收敛”和“实现的舍入稳定性”分开了；减小 $\alpha$ 只有配合增量形式才会转化为实际精度。论文把更系统的舍入分析指向 Wu, Yang and Zhou (2025)。

对于一般多步法，(3.66) 未必成立。Volterra 偏积分微分方程会产生稠密下三角 Toeplitz 时间矩阵；若积分权重正且单调，已有结果给出 $|\lambda(P_\alpha^{-1}K)|=1+O(\alpha)$。

### 非线性 Newton–Krylov

对 $\boldsymbol u'=f(t,\boldsymbol u)$ 使用后向 Euler，先在非线性全时间方程上做 Newton：

$$
J\Delta\boldsymbol U^l=\boldsymbol b-F(\boldsymbol U^l),
\qquad
\boldsymbol U^{l+1}=\boldsymbol U^l+\Delta\boldsymbol U^l. \tag{3.68}
$$

其中

$$
J=B\otimes I_x-\nabla F_l,
\qquad
\nabla F_l=\operatorname{blkdiag}\!\left(
\nabla f(\boldsymbol u_1^l,t_1),\ldots,
\nabla f(\boldsymbol u_{N_t}^l,t_{N_t})
\right),
$$

$$
B=\frac1{\Delta t}
\begin{bmatrix}
1\\-1&1\\&\ddots&\ddots\\&&-1&1
\end{bmatrix}.
$$

内层 GMRES 使用

$$
P_\alpha=C_\alpha\otimes I_x-I_t\otimes A_l,
$$

$A_l$ 为所有 $\nabla f(\boldsymbol u_n^l,t_n)$ 的平均。一般有 $\rho(P_\alpha^{-1}J)>1$，所以定常迭代 (3.61) 不适合该 Jacobian 系统；预条件谱仍会聚集，GMRES 因而有效。缩短时间窗会增强聚集并降低内层迭代数，Section 3.5.1 的最近 Kronecker 近似还可进一步改善预条件。

> [!tip] 本站洞见：非线性情形正是"$\rho\ge1$ 但 GMRES 仍快"的典型
> 这里 $P_\alpha$ 用一个**平均** Jacobian $A_l$ 代替随时间变化的 $\nabla f(\boldsymbol u_n^l,t_n)$，近似误差通常大到让 $\rho(P_\alpha^{-1}J)>1$——定常迭代 (3.61) 直接发散。但只要各时间点的 Jacobian 变化不剧烈，$P_\alpha^{-1}J$ 的谱仍聚在 $1$ 附近的一小团，按本页开头的多项式界，GMRES 依旧几步收敛。这正呼应了 Krylov 与定常迭代的分野：定常迭代被 $\rho$ 判死刑，Krylov 只看谱能否被低次多项式压住。时间窗越短，$A_l$ 越接近真实 Jacobian、聚集越紧、内层步数越少；Section 3.5.1 的最近 Kronecker 近似则通过更贴近 $\nabla F_l$ 的低秩结构进一步收紧谱。

## 公式与图表覆盖核对

| 原文项目                                             | 论文小节 | 覆盖状态                                                              |
| ---------------------------------------------------- | -------- | --------------------------------------------------------------------- |
| 多步全时间系统、循环替换、(3.49)–(3.52), Theorem 3.8 | 3.5.2    | 离散起点、FFT 同时对角化、三步求解、谱聚集与 Krylov/定常之别          |
| (3.53)–(3.54), Figure 3.15                           | 3.5.2    | 三类 PDE、Numerov 稳定阈值与必要性、谱与 GMRES                        |
| (3.55)–(3.60), Figure 3.16                           | 3.5.2    | 连续首尾条件、$\alpha$-循环矩阵、对角缩放式 FFT 实现与实验            |
| (3.61)–(3.65)                                        | 3.5.2    | 两条路线等价（McDonald / Banjai–Peterseim）、一阶与直接二阶全时间系统 |
| Theorem 3.9, (3.66)–(3.67), Figure 3.17              | 3.5.2    | 稳定性条件、原文矛盾标注、$\rho(\mathcal M)$ 界与阈值/时间窗实验      |
| Figure 3.18                                          | 3.5.2    | $\alpha$ 的收敛–舍入权衡、$O(\epsilon/\alpha)$、直接解与增量更新      |
| (3.68)                                               | 3.5.2    | 非线性 Newton、$\alpha$-循环预条件 GMRES、$\rho>1$ 仍快与时间窗影响   |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3.5.2, pp. 431–443.
