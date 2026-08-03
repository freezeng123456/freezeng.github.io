---
title: 3.5.1：直接 ParaDiag（ParaDiag-I）
description: 从全时间系统、几何时间网格的截断/舍入平衡，到 BVM 消元与非线性准 Newton（含 NKA）的完整推导
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i
tags:
  - 时间并行
  - ParaDiag
  - 全时间系统
---

> [!note] 阅读范围
> 本页对应论文 Sections 3.5–3.5.1（pp. 415–431），覆盖 ParaDiag 总体分类、公式 (3.22)–(3.48)、Theorems 3.5–3.7、Figures 3.9–3.14 和 Tables 3.1–3.2。直接法的截断误差、舍入误差、几何网格、BVM 以及非线性准 Newton 分支均逐步展开，重点还原「离散误差 vs 舍入误差」这一贯穿全节的设计取舍。

## 3.5.1 直接 ParaDiag 方法（ParaDiag-I）

ParaDiag 是本章介绍的最后一类方法，核心思想是把**时间步进矩阵**（或其近似）对角化。之所以能带来并行性，是因为一旦时间方向被对角化，各个时间层之间的耦合就被解除，每一层退化成一个仅含空间算子的独立系统。按照对时间矩阵的处理方式，ParaDiag 分成两支。

### ParaDiag-I 与 ParaDiag-II 的分界

ParaDiag-I（Maday and Rønquist 2008）**精确**对角化时间步进矩阵，因此和 ParaExp 一样是一种直接（非迭代）时间并行求解器。对角化能成立有一个硬性前提：时间矩阵必须可对角化。对最朴素的均匀步长后向 Euler，时间矩阵是一个退化的 Jordan 块，根本不可对角化，所以 ParaDiag-I 必须**要么使用非均匀时间步**，**要么在最后一步换用另一种时间积分公式**（如边值方法）。它的代价是：非均匀步长下的误差分析（Gander et al. 2019）表明，单个时间窗内可并行的步数受舍入误差限制，双精度下通常只有约二十步；而且目前只对后向 Euler 和梯形规则这类低阶积分器有系统结果，很难推广到 Runge–Kutta 等高阶方法。改用边值方法型离散能显著放宽步数上限，但可用的积分器仍局限于后向 Euler 与梯形规则（Liu et al. 2022）。

ParaDiag-II（Gander et al. 2021c）则**近似**时间矩阵，使其变得良态可对角化，再把这个近似放进定常迭代或用作 Krylov 方法的预条件，于是「以迭代换直接性」。其预条件的设计有两条原则：一是对角化必须良态（特征向量矩阵条件数小），以压住对角化步骤引入的舍入误差——这正是 ParaDiag-I 所欠缺的；二是迭代要收敛快（谱半径小，或预条件后谱紧密聚集在 1 附近以利 Krylov 加速）。这一迭代型 ParaDiag 最早由 McDonald, Pestana and Wathen (2018) 在离散层面提出，并由 Gander and Wu (2019) 在连续层面独立给出。ParaDiag-II 的完整讨论留到下一页。

> [!tip] 本站洞见
> 两支的分野本质上是一次**精度—稳健性**的权衡：ParaDiag-I 坚持精确分解，换来「无外层迭代、误差就是离散误差」的干净直接性，但把全部风险压在特征向量矩阵 $V$ 的条件数上，一旦 $\operatorname{Cond}(V)$ 随窗口增大而爆炸，直接性反而成了枷锁；ParaDiag-II 主动放弃精确性，用一个良态近似换取任意长窗口和任意积分器，代价只是可控的几步迭代。后文关于几何网格、BVM 和 $\varrho$ 选取的每一个技术细节，都是在为 ParaDiag-I 这条「精确但脆弱」的路线争取尽可能大的可用步数。

### 后向 Euler 的全时间系统

对 $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$ 使用变步长后向 Euler：

$$
\frac{\boldsymbol u_n-\boldsymbol u_{n-1}}{\Delta t_n}
=A\boldsymbol u_n+\boldsymbol g_n,
\qquad n=1,\ldots,N_t. \tag{3.22}
$$

传统做法是按 $n$ 逐层前进，天然串行。ParaDiag 的第一步转变是**放弃逐层求解**，把 $N_t$ 个差分方程堆叠成一个「全时间（all-at-once）」系统，一次性求解所有时间层。令 $\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top$，得到

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

这里 $I_x\in\mathbb R^{N_x\times N_x}$、$I_t\in\mathbb R^{N_t\times N_t}$ 是单位阵。Kronecker 结构精确刻画了「时间 $\times$ 空间」的张量分解：$B\otimes I_x$ 只作用在时间方向（$B$ 是后向 Euler 的时间差分算子，下双对角，对角为 $1/\Delta t_n$、次对角为 $-1/\Delta t_n$），$I_t\otimes A$ 只作用在空间方向。初值 $\boldsymbol u_0$ 只出现在第一块右端，说明这是一个纯粹的初值问题被摊平成了空间维度上的耦合系统。

之所以能对角化，关键在于**只有时间方向的 $B$ 需要处理**。若所有 $\Delta t_n$ 互不相同，$B$ 的对角元 $1/\Delta t_n$ 两两互异，作为一个三角矩阵它就有 $N_t$ 个互异特征值，从而可对角化：

$$
B=VDV^{-1},
\qquad
D=\operatorname{diag}\!\left(
\frac1{\Delta t_1},\ldots,\frac1{\Delta t_{N_t}}
\right). \tag{3.24}
$$

把 $B=VDV^{-1}$ 代入 $K$，并利用 Kronecker 积的混合积性质 $(V\otimes I_x)(D\otimes I_x)(V^{-1}\otimes I_x)=VDV^{-1}\otimes I_x=B\otimes I_x$，同时 $I_t\otimes A=(V\otimes I_x)(I_t\otimes A)(V^{-1}\otimes I_x)$（因为 $V I_t V^{-1}=I_t$），得到块状分解

$$
K=(V\otimes I_x)
(D\otimes I_x-I_t\otimes A)
(V^{-1}\otimes I_x).
$$

中间因子 $D\otimes I_x-I_t\otimes A$ 是**块对角**的：第 $n$ 个对角块正是 $\tfrac1{\Delta t_n}I_x-A$，时间层之间再无耦合。于是求解 $K\boldsymbol U=\boldsymbol b$ 分成三步：

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

其中 $\boldsymbol U^a=((\boldsymbol u_1^a)^\top,\ldots,(\boldsymbol u_{N_t}^a)^\top)^\top$，$\boldsymbol U^b$ 同理。第一步（step a）与第三步（step c）只是把 $V^{-1}$、$V$ 沿时间方向作用在数据上，即对每个空间自由度做一次 $N_t\times N_t$ 的稠密矩乘，成本相对低廉。真正的主计算是第二步（step b）：$N_t$ 个移位空间系统 $\bigl(\tfrac1{\Delta t_n}I_x-A\bigr)\boldsymbol u_n^b=\boldsymbol u_n^a$。**它们彼此完全独立**，可以分配到 $N_t$ 个处理器上并行求解——这正是「对角化 = 时间并行」的全部机制。

![ParaDiag 的时间变换、独立空间求解与逆变换](assets/diagrams/pint/zh/paradiag-three-stage.svg)

> [!tip] 本站洞见
> 三步法可以类比成时间方向上的一次「谱变换—对角求解—逆变换」，与用 FFT 求解常系数循环卷积如出一辙：step a/step c 是变换与反变换，step b 是频域里的逐点（逐特征值）求解。差别在于循环矩阵用的是酉的 Fourier 基（条件数恒为 1），而这里的 $V$ 是一般三角 Toeplitz 特征向量矩阵，**条件数可能极大**。ParaDiag-I 的全部风险因此集中在 step a/step c：一个病态的 $V$ 会把 step b 里干净的解重新污染回时间域。这解释了为什么后文所有努力都围绕「如何让 $V$ 尽量良态」展开。

### 几何网格与两类误差

对任意步长 $\{\Delta t_n\}$，$V$ 一般没有闭式，只能靠数值 `eig` 得到。这在实践中并不昂贵（$N_t$ 不必很大），但也让完整的舍入误差分析和参数选取无从下手。为了既保证互异特征值、又能解析分析，Maday and Rønquist (2008) 采用**几何网格** $\Delta t_n=\mu^{n-1}\Delta t_1$，$\mu>1$。他们在一维热方程上取 $\mu=1.2$，得到近乎完美的加速比。由 $\sum_n\Delta t_n=\sum_n\mu^{n-1}\Delta t_1=T$ 定出首步 $\Delta t_1=T/\sum_{j=1}^{N_t}\mu^{j-1}$，于是

$$
\Delta t_n=\frac{\mu^{n-1}}{\sum_{j=1}^{N_t}\mu^{j-1}}T. \tag{3.26}
$$

几何网格把「是否可对角化」变成一个连续可调的问题。令 $\mu=1+\varrho$，则 $\varrho$ 直接度量网格偏离均匀的程度，两端各有一个失效机制：

- $\varrho$ **过大**：末端步长 $\Delta t_{N_t}=\mu^{N_t-1}\Delta t_1$ 相对首步指数级膨胀，网格严重非均匀，后段时间层的**截断（离散）误差**急剧上升。
- $\varrho$ **过小**（$\mu\to1$）：所有 $\Delta t_n$ 趋于相等，$B$ 的特征值 $1/\Delta t_n$ 相互靠拢，$B$ 逼近一个单一 Jordan 块。Jordan 块是缺陷（defective）矩阵，其特征向量矩阵 $V$ 的条件数发散，对角化时**舍入误差**被 $\operatorname{Cond}(V)$ 放大。

直接 ParaDiag 的参数选择，本质就是在这两项相互对立的误差之间取平衡。

> [!tip] 本站洞见
> 这里藏着一个耐人寻味的反直觉：数值分析通常鼓励**均匀细分**以减小截断误差，但在 ParaDiag-I 里均匀网格恰恰是最坏选择——它让 $B$ 退化成不可对角化的 Jordan 块，方法直接崩溃。几何网格是一种「为了可对角化而故意引入的网格畸变」，$\varrho$ 越大越远离 Jordan 块（$V$ 越良态），却也越偏离等距（截断误差越大）。因此 $\varrho$ 不是普通意义上「越小越准」的离散参数，而是一个必须**从下方托住**的稳健性参数，这正是它与常规步长参数在直觉上的根本区别。

### Theorem 3.5：一阶问题的平衡公式

上述取舍在 Gander et al. (2016a) 中对一阶抛物问题被精确量化。假设 $\sigma(A)\subset\mathbb R_-$ 且 $|\lambda(A)|\leq\lambda_{\max}$。记 $\boldsymbol u_{N_t}(\varrho)$ 和 $\boldsymbol u_{N_t}(0)$ 为几何网格与均匀网格在 $t=T$ 的后向 Euler 解，$\widetilde{\boldsymbol u}_n(\varrho)$ 为对角化方法 (3.25) 的实际计算值。则

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

$\widetilde x_*$ 是 $r(\widetilde x,N_t)$ 在 $[0,\infty)$ 上的最大点，$\lambda_*=N_t\widetilde x_*/T$。两条估计正是前述两类误差：**第一条正比于 $\varrho^2$**，是几何网格相对均匀网格多付出的截断误差（随 $\varrho$ 增大而增大）；**第二条正比于 $\varrho^{-(N_t-1)}$**，是对角化带来的舍入误差（随 $\varrho$ 减小而爆炸），其中显式带有机器精度 $\epsilon$。让两项相等，即解 $C\varrho^2=\epsilon(\cdots)\varrho^{-(N_t-1)}$，得到最优参数

$$
\varrho_{\mathrm{opt}}=
\left(
\epsilon
\frac{N_t^2(2N_t+1)(N_t+\lambda_{\max}T)}
{\phi(N_t)C(\lambda_*T,N_t)}
\right)^{1/(N_t+1)}. \tag{3.28}
$$

证明的思路是把耦合系统**逐特征值解耦**：取 $\lambda\in\sigma(A)$，全时间系统在该特征方向上退化为 Dahlquist 检验方程 $y'=\lambda y$。第一条估计来自 Gander et al. (2016a, Theorem 2)，比较几何网格与均匀网格对该标量方程的离散误差；第二条来自同文 Theorem 6，估计特征向量变换引入的舍入误差，并且当 $|\lambda|=\lambda_{\max}$ 时取到最坏值（谱半径最大的模态最敏感）。第一条同时给出了到精确解的截断误差：由三角不等式
$\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u(T)\|\leq\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u_{N_t}(0)\|+\|\boldsymbol u_{N_t}(0)-\boldsymbol u(T)\|$，
其中后一项是标准均匀后向 Euler 的误差、已被充分理解且不占主导，因此几何网格的额外代价完全由 $\varrho^2$ 项刻画。

几何网格下 $V$、$V^{-1}$ 有闭式，均为下三角 Toeplitz 矩阵（Gander et al. 2016a）：

$$
V=\mathbb T(p_1,\ldots,p_{N_t-1}),
\qquad
p_n=\frac1{\prod_{j=1}^{n}(1-\varrho^j)},
$$

$$
V^{-1}=\mathbb T(q_1,\ldots,q_{N_t-1}),
\qquad
q_n=(-1)^n\varrho^{n(n-1)/2}p_n, \tag{3.29a}
$$

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

这些闭式的价值不在于实际求解，而在于**解析估计 $\operatorname{Cond}(V)$**，从而把 (3.27) 第二条的舍入常数写成 $N_t$ 的显式函数（$\phi(N_t)$ 里的阶乘正来自 $p_n$、$q_n$ 的连乘）。真正实现 (3.25) 时并不套用这些公式，而是调用 `eig`：它会自动缩放特征向量来改善条件数，得到比理论闭式更好的实际 $\operatorname{Cond}(V)$。单精度与双精度的机器精度分别约为 $1.19\times10^{-7}$ 和 $2.22\times10^{-16}$。

> [!tip] 本站洞见
> (3.28) 里 $\varrho_{\mathrm{opt}}\sim\epsilon^{1/(N_t+1)}$ 的指数 $1/(N_t+1)$ 揭示了 ParaDiag-I「约二十步天花板」的数学根源。一方面，随 $N_t$ 增大，$\epsilon^{1/(N_t+1)}\to1$，最优 $\varrho$ 被迫增大、网格越来越畸变；另一方面，舍入常数里的 $\phi(N_t)$ 含**阶乘级**增长，$\varrho^{-(N_t-1)}$ 又是指数放大。两股力量叠加，使得能压到机器精度以下的窗口宽度只能是 $O(\log(1/\epsilon))$ 量级——双精度下恰好落在二十步附近。换句话说，天花板不是实现缺陷，而是「用一个病态相似变换换取时间并行」这一策略的内在信息论极限。

### Figures 3.9–3.10：最优参数仍有并行宽度上限

![原论文 Figure 3.9：热方程与对流扩散方程上误差随几何参数变化](assets/papers/time-parallelization/source-figures/figure-3-9.svg)

Figure 3.9 使用齐次 Dirichlet 边界、$u_0(x)=\sin(2\pi x)$、$\Delta x=1/50$、$T=0.2$。右图的对流扩散黏性为 $10^{-2}$。五组 $N_t$ 分别扫描 $\varrho\in[10^{-2},1]$，误差取所有时间节点上的最大 $L^\infty$ 误差，参照解由指数积分子 $\boldsymbol u(t_n)=e^{-At_n}\boldsymbol u_0$ 给出。每条曲线都呈明显的 V 形——左支是舍入主导、右支是截断主导，谷底即最优 $\varrho$，星号是 (3.28) 的理论预测。预测对对流扩散很准，对小 $N_t$ 的热方程略有偏差，这与 (3.28) 只是**渐近**平衡两条估计相符。

![原论文 Figure 3.10：最优几何参数下误差随时间步数先降后升](assets/papers/time-parallelization/source-figures/figure-3-10.svg)

Figure 3.10 改取 $T=0.5$ 和 $N_t=2^4,2^5,\ldots,2^{10}$，并对每个 $N_t$ 用数值最优 $\varrho_{\mathrm{num}}$。均匀后向 Euler 的误差继续随 $N_t$ 单调下降（更细的时间步更准）；而 ParaDiag-I 的误差先降后升，在不足 100 步处越过阈值并迅速恶化。这条非单调曲线是前述天花板的直接实验体现：一旦 $\operatorname{Cond}(V)$ 放大的舍入误差压过时间离散误差的下降，增加步数反而有害。

### 波动方程与梯形规则

同一框架可推广到二阶系统。考虑波动方程空间离散后的

$$
\boldsymbol u''(t)=A\boldsymbol u(t),
\quad
\boldsymbol u(0)=\boldsymbol u_0,
\quad
\boldsymbol u'(0)=\widetilde{\boldsymbol u}_0 \tag{3.30}
$$

（波动方程中 $A$ 是离散 Laplace 算子）。为套用 ParaDiag-I，先降为一阶系统

$$
\boldsymbol w'=\mathbb A\boldsymbol w,
\qquad
\boldsymbol w=(\boldsymbol u^\top,(\boldsymbol u')^\top)^\top,
\qquad
\mathbb A=
\begin{bmatrix}0&I_x\\A&0\end{bmatrix}. \tag{3.31}
$$

后向 Euler 在波动问题上有过强的数值耗散、破坏能量守恒，会引入色散伪影。为此改用**梯形规则**

$$
\frac{\boldsymbol w_n-\boldsymbol w_{n-1}}{\Delta t_n}
=\frac{\mathbb A}{2}(\boldsymbol w_n+\boldsymbol w_{n-1}), \tag{3.32}
$$

它是能量守恒的（$\|\boldsymbol w_n\|^2=\|\boldsymbol w_0\|^2$），因而适合长时间波传播。步长仍取几何网格 (3.26)。梯形规则的两点平均使右端同时含 $\boldsymbol w_n$ 与 $\boldsymbol w_{n-1}$，因此全时间系统里除了时间差分矩阵 $B$，还多出一个**时间平均矩阵** $\widetilde B$：

$$
K\boldsymbol W=\boldsymbol b,
\qquad
K=B\otimes I_x-\widetilde B\otimes\mathbb A, \tag{3.33a}
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

此时 $K$ 含两个不同的时间矩阵 $B$、$\widetilde B$，不能直接套用单矩阵对角化。技巧是**左乘 $\widetilde B^{-1}\otimes I_x$**，把 $\widetilde B$ 归一化成单位阵，从而恢复 (3.23a) 的标准形：

$$
\mathcal K\boldsymbol W=\widetilde{\boldsymbol b},
\qquad
\mathcal K=\widetilde B^{-1}B\otimes I_x-I_t\otimes\mathbb A,
\qquad
\widetilde{\boldsymbol b}=(\widetilde B^{-1}\otimes I_x)\boldsymbol b. \tag{3.34}
$$

现在只需对角化单一时间矩阵 $\widetilde B^{-1}B$。Gander et al. (2019) 给出其闭式对角化：

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

特征值 $2/\Delta t_n$ 依旧互异（几何网格所保证），$V$ 依旧是下三角 Toeplitz，于是 (3.34) 仍可按三步 (3.25) 求解。

### Theorem 3.6：波动问题的平衡公式

对 $\lambda(A)\leq0$，几何网格与均匀网格的梯形解、以及对角化计算值满足

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

结构与一阶情形完全一致——截断项 $\propto\varrho^2$、舍入项 $\propto\varrho^{-(N_t-1)}$，令二者相等给出

$$
\varrho_{\mathrm{opt}}=
\left(
\epsilon\frac{15\times2^{2N_t-1/2}}
{(N_t^2-1)(N_t-1)!}
\right)^{1/(N_t+1)}. \tag{3.37}
$$

证明同样逐特征值进行：对 $-A$ 的每个特征值 $\lambda>0$ 考察标量振子 $u''+\lambda u=0$。由 Gander et al. (2019, Theorem 2.1)，几何与均匀网格的截断误差之差为 $O\!\bigl(\tfrac{N_t(N_t^2-1)}{6}r_1(\tfrac{\lambda T}{2N_t})\varrho^2\bigr)$，其中

$$
r_1(s)=\frac{s^3}{(1+s^2)^2},
\qquad r_1(s)\leq\frac25\ (s\geq0),
$$

代入即得 (3.36) 第一条（$\tfrac{N_t(N_t^2-1)}{6}\cdot\tfrac25=\tfrac{N_t(N_t^2-1)}{15}$）。舍入误差由同文 Theorem 2.11 给出，含因子 $r_2(s)=1/(1+s^2)\leq1$，一致放缩即得第二条。两个**一致上界** $r_1\leq2/5$、$r_2\leq1$ 消去了对具体特征值 $\lambda$ 的依赖，使界只与 $N_t,\varrho,\epsilon$ 有关。

![原论文 Figure 3.11：波动方程上几何梯形 ParaDiag-I 的最优参数与步数阈值](assets/papers/time-parallelization/source-figures/figure-3-11.svg)

Figure 3.11 使用齐次 Dirichlet 边界、$\Delta x=1/20$、$T=0.2$。(a) 中每组 $N_t$ 都有误差最小点，(3.37) 的星号接近实测最优值；(b) 使用数值最优参数后，$N_t>32$ 时误差迅速恶化，与一阶情形的天花板现象一致。

![原论文 Table 3.1：后向 Euler 和梯形规则时间矩阵的特征向量条件数](assets/papers/time-parallelization/source-figures/table-3-1.svg)

Table 3.1 直接测量了失效的元凶——$\operatorname{Cond}(V)$。在 $N_t=5,10,20,30,60,100$ 时，后向 Euler 的 $B$ 从 $1.7\times10^3$ 增至 $4.8\times10^6$；梯形规则的 $\widetilde B^{-1}B$ 从 $4.7\times10^3$ 增至 $4.1\times10^9$。条件数随 $N_t$ 快速上升，与 (3.27)、(3.36) 的舍入分析吻合，也解释了 Figures 3.10–3.11 的误差恶化。一个有趣的观察是：使用**数值最优** $\varrho_{\mathrm{num}}$（而非理论 $\varrho_{\mathrm{opt}}$）时，$\operatorname{Cond}(V)$ 在后段出现**平台**而非持续爆炸——这说明 `eig` 的自动缩放在实践中比闭式分析更能压住条件数，作者指出这一现象值得进一步研究。

### BVM：固定步长并改变末步公式

几何网格的困境在于：为了可对角化必须畸变网格，而畸变本身又限制了步数。Liu et al. (2022) 换了一条思路——**保持均匀步长**，转而让**最后一步换用不同公式**，从而绕开 Jordan 块。前 $N_t-1$ 步用二阶中心差分，最后一步用后向 Euler：

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

中心差分在时间步进模式下是不稳定的（它同时用到 $\boldsymbol u_{n-1}$ 和 $\boldsymbol u_{n+1}$，无法逐层前进）；这里之所以可行，是因为整个系统**同时求解**，属于**边值方法（boundary value method, BVM）**。它的稳定性不能按普通逐步格式判断——用末步的一阶后向 Euler「封口」，为整个中心差分骨架提供了缺失的边界条件（见 Gander 2015, Section 5.2 对稳定性的讨论）。Axelsson and Verwer (1985) 正是用边值技术绕过著名的 **Dahlquist 收敛—稳定性壁垒**：他们证明即便末步只有一阶，同时求解得到的数值解在一般非线性情形下仍是**一致二阶精度**（Axelsson and Verwer 1985, Theorem 4）。更早的 Fox (1954) 与 Fox and Mitchell (1957) 已用过这类离散，只是末步换成 BDF2：

$$
\frac{3\boldsymbol u_{N_t}-4\boldsymbol u_{N_t-1}+\boldsymbol u_{N_t-2}}{2\Delta t}=A\boldsymbol u_{N_t}+\boldsymbol g_{N_t}.
$$

方法 (3.38) 是后来所谓 BVM 的典型代表，其全时间系统的适定性由 Brugnano, Mazzia and Trigiante (1993)（另见 Brugnano and Trigiante 2003）严格论证；在 BVM 传统里，这类全时间系统通常靠构造有效预条件迭代求解，而 ParaDiag-I 则用对角化直接求解。

全时间形式仍为 $K\boldsymbol U=\boldsymbol b$，$K=B\otimes I_x-I_t\otimes A$，只是 $B$ 换成中心差分骨架加后向 Euler 末行：

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

由于步长均匀，只需初值 $\boldsymbol u_0$，全部时间层一次求解。关键在于其对角化的良态性：

**Theorem 3.7.** (3.39b) 的 $B$ 可分解为 $B=VDV^{-1}$，且 $\operatorname{Cond}(V)=O(N_t^2)$（$V$、$V^{-1}$、$D$ 的闭式见 Liu et al. 2022, Section 3）。

这是 BVM 相对几何网格的决定性优势：几何网格的 $\operatorname{Cond}(V)$ 随 $N_t$ **指数级**恶化（Table 3.1），而 BVM 只是**多项式**（$O(N_t^2)$）增长，因此单窗口可并行步数大幅提高。

> [!tip] 本站洞见
> 几何网格与 BVM 是达成「可对角化」的两条互补路线：几何网格保持**同一积分器**、靠**扰动步长**拉开特征值，代价是特征向量指数病态；BVM 保持**均匀步长**、靠**改一格积分器**拉开特征值，代价仅是末步降阶（且被 Axelsson–Verwer 的一致二阶结论补偿）。前者把畸变加在网格上，后者把畸变加在格式上——而把畸变放进末步公式远比放进整条网格温和，这正是 $O(N_t^2)$ 对指数的本质改进来源。

二阶系统同样可先降为一阶再用 BVM。设 $\boldsymbol v=\boldsymbol u'$、$\boldsymbol w=(\boldsymbol u^\top,\boldsymbol v^\top)^\top$，对 $\boldsymbol w$ 用同一 BVM：

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

但引入速度变量 $\boldsymbol v$ 会让每个时间点的空间存储翻倍，在高维或细网格下代价可观。为此可**消去 $\boldsymbol v$**，只保留 $\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top$。为看清这一点，在离散层面把 (3.40) 按 $\{\boldsymbol u_n\}$ 与 $\{\boldsymbol v_n\}$ 分别写出（$\mathbb A$ 的上半行给出 $\boldsymbol u'=\boldsymbol v$，下半行给出 $\boldsymbol v'=A\boldsymbol u$）：

$$
\left\{
\begin{aligned}
\frac{\boldsymbol u_{n+1}-\boldsymbol u_{n-1}}{2\Delta t}&=\boldsymbol v_n,\quad n<N_t,\\
\frac{\boldsymbol u_{N_t}-\boldsymbol u_{N_t-1}}{\Delta t}&=\boldsymbol v_{N_t},
\end{aligned}
\right.
\qquad
\left\{
\begin{aligned}
\frac{\boldsymbol v_{n+1}-\boldsymbol v_{n-1}}{2\Delta t}&=A\boldsymbol u_n,\quad n<N_t,\\
\frac{\boldsymbol v_{N_t}-\boldsymbol v_{N_t-1}}{\Delta t}&=A\boldsymbol u_{N_t}.
\end{aligned}
\right.
$$

用 (3.39b) 的 $B$，这两组恰好写成

$$
(B\otimes I_x)\boldsymbol U-\boldsymbol V=\boldsymbol b_1,
\qquad
(B\otimes I_x)\boldsymbol V-A\boldsymbol U=\boldsymbol b_2,
$$

其中 $\boldsymbol b_1=(\tfrac{\boldsymbol u_0^\top}{2\Delta t},0,\ldots,0)^\top$、$\boldsymbol b_2=(\tfrac{\widetilde{\boldsymbol u}_0^\top}{2\Delta t},0,\ldots,0)^\top$（初值只进入第一行）。由第一式解出 $\boldsymbol V=(B\otimes I_x)\boldsymbol U-\boldsymbol b_1$，代入第二式得 $(B\otimes I_x)^2\boldsymbol U-A\boldsymbol U=\boldsymbol b_2+(B\otimes I_x)\boldsymbol b_1$。利用 $(B\otimes I_x)^2=B^2\otimes I_x$，并直接计算右端 $\boldsymbol b_2+(B\otimes I_x)\boldsymbol b_1=\boldsymbol b$，即得只含 $\boldsymbol U$ 的全时间系统

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

时间矩阵变成 $B^2$，其特征向量矩阵与 $B$ 同（$B^2=VD^2V^{-1}$），因此仍继承 Theorem 3.7 的 $O(N_t^2)$ 良态，可照常三步求解，且存储不再翻倍。

![原论文 Figure 3.12：波动方程上几何时间网格与 BVM 的误差和条件数](assets/papers/time-parallelization/source-figures/figure-3-12.svg)

Figure 3.12 使用 $T=0.5$、$\Delta x=1/40$ 和齐次 Dirichlet 边界，比较几何梯形 ParaDiag-I 与 BVM ParaDiag-I（$N_t=2^2,\ldots,2^8$）。几何梯形在 $N_t\approx32$ 开始出现典型的舍入恶化；BVM 则保持 $O(\Delta t^2)$、毫无恶化，与串行梯形规则一致。右图给出对应特征向量条件数，BVM 显著更低，直接印证了「无恶化」的来源。

### 非线性全时间方程与准 Newton

以上都是线性问题。对非线性 $\boldsymbol u'=\boldsymbol f(\boldsymbol u,t)$（二阶可类似处理），定义

$$
F(\boldsymbol U)=
\left(
\boldsymbol f(\boldsymbol u_1,t_1)^\top,
\ldots,
\boldsymbol f(\boldsymbol u_{N_t},t_{N_t})^\top
\right)^\top,
$$

全时间方程为

$$
(B\otimes I_x)\boldsymbol U-F(\boldsymbol U)=\boldsymbol b, \tag{3.42}
$$

其中 $B$ 既可取几何网格的 (3.23b)，也可取 BVM 的 (3.39b)。对非线性系统施加 Newton 法，原始更新为
$(B\otimes I_x-\nabla F(\boldsymbol U^k))(\boldsymbol U^{k+1}-\boldsymbol U^k)=\boldsymbol b-((B\otimes I_x)\boldsymbol U^k-F(\boldsymbol U^k))$，
整理成便于对角化的形式

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

问题在于：$\nabla F$ 是**逐时间层不同**的块对角阵，$K=B\otimes I_x-\nabla F$ 不再具有 $I_t\otimes(\cdot)$ 的干净 Kronecker 结构，无法用单一 $V$ 对角化，时间并行随之失效。解决办法（受 Gander and Halpern 2017 启发）是用**单一平均 Jacobian** $A_k$ 替换所有块：

$$
A_k=\frac1{N_t}\sum_{n=1}^{N_t}\nabla f(\boldsymbol u_n^k,t_n),
\qquad\text{或}\qquad
A_k=\nabla f\!\left(
\frac1{N_t}\sum_{n=1}^{N_t}\boldsymbol u_n^k,
\frac{T}{N_t}
\right). \tag{3.44}
$$

于是 $\nabla F(\boldsymbol U^k)\approx I_t\otimes A_k$，代回 (3.43a) 得准 Newton 迭代

$$
(B\otimes I_x-I_t\otimes A_k)\boldsymbol U^{k+1}
=\boldsymbol b-\left((I_t\otimes A_k)\boldsymbol U^k-F(\boldsymbol U^k)\right). \tag{3.45}
$$

Kronecker 结构恢复，对角化 $B=VDV^{-1}$ 后每轮仍按三步求解：

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

其中 $\boldsymbol r^k=\boldsymbol b-((I_t\otimes A_k)\boldsymbol U^k-F(\boldsymbol U^k))$ 是 (3.45) 的右端。线性情形 $A_k=A$、$\boldsymbol r^k=\boldsymbol b$，(3.46) 退化为 (3.25)。这类以近似 Jacobian 求解的准 Newton 法收敛性已被充分研究（Deuflhard 2004, Theorem 2.5；Ortega and Rheinboldt 2000）。其收敛速度取决于单一 $A_k$ 对全部 $N_t$ 个 Jacobian 块的逼近好坏：若 $\nabla f(\boldsymbol u_n^k,t_n)$ 随时间剧烈变化，任何单矩阵都无法同时逼近，此时应缩短时间窗、按窗串行。

![原论文 Figure 3.13：BVM ParaDiag-I 在两组黏性和多个时间窗长度下的 Burgers 收敛](assets/papers/time-parallelization/source-figures/figure-3-13.svg)

Figure 3.13 使用周期 Burgers 方程、$\Delta x=0.01$，并保持 $N_t=T/\Delta t=200$。横线为空时离散误差 $\max\{\Delta t^2,\Delta x^2\}=10^{-4}$。$\nu=0.1$ 时，$T=0.1$ 到 $1.6$ 的收敛速度变化很小，说明 Jacobian 沿时间变化温和、单一 $A_k$ 足够好；$\nu=0.002$ 时，时间窗变长明显恶化，$T=0.8,1.6$ 未收敛到目标线——正是 Jacobian 剧烈变化导致平均近似失效的例证。

![原论文 Table 3.2：串行梯形规则与并行 BVM ParaDiag-I 的 Jacobian 求解次数](assets/papers/time-parallelization/source-figures/table-3-2.svg)

Table 3.2 量化了并行收益。若有 $N_t$ 个处理器，(3.46) 中 $N_t$ 个 Jacobian 系统同时求解，因此**并行 Jacobian 求解次数就等于外层准 Newton 轮数 $k$**；而串行时间步进要解 $\sum_n It_n$ 个系统（$It_n$ 为第 $n$ 步的 Newton 步数）。表中 $\nu=0.1$ 时串行梯形需 401–443 次，ParaDiag-I 只需 5–7 轮；$\nu=0.002$ 时 ParaDiag-I 从 7 增至 22 轮，并在更长窗口失效——与 Figure 3.13 的结论一致：收敛快时并行优势极其显著。

### 最近 Kronecker 近似

单一 $I_t\otimes A_k$ 用**同一** $A_k$ 逼近每个时间层，忽略了 Jacobian 的**幅值**随时间的变化。Liu and Wu (2022b, Section 3.3) 提出的改进（原为加速非线性 ParaDiag-II，同样适用于此）保留张量结构但允许逐层缩放：用 $\Phi_k\otimes A_k$ 逼近 $\nabla F$，其中 $\Phi_k=\operatorname{diag}(\phi_1,\ldots,\phi_{N_t})$ 由极小化确定

$$
\min_{\Phi_k\ \mathrm{diagonal}}
\left\|\nabla F(\boldsymbol U^k)-\Phi_k\otimes A_k\right\|. \tag{3.47}
$$

在 Frobenius 范数下，这正是**最近 Kronecker 积近似（NKA）**，解为（Van Loan and Pitsianis 1993, Theorem 3，设 $\operatorname{trace}(A_k^\top A_k)>0$）

$$
\phi_n=
\frac{\operatorname{trace}\!\left(
\nabla f(\boldsymbol u_n^k,t_n)A_k^\top
\right)}
{\operatorname{trace}(A_k^\top A_k)},
\qquad n=1,\ldots,N_t. \tag{3.48}
$$

即把每个 Jacobian 块正交投影到 $A_k$ 方向上，$\phi_n$ 是最优标量系数。新准 Newton 方程为

$$
(B\otimes I_x-\Phi_k\otimes A_k)\boldsymbol U^{k+1}
=\boldsymbol b-\left((\Phi_k\otimes A_k)\boldsymbol U^k-F(\boldsymbol U^k)\right).
$$

由于时间矩阵不再是 $B$ 而是 $B$ 与 $\Phi_k$ 的组合，左乘 $B^{-1}\otimes I_x$ 化为
$(I_t\otimes I_x-B^{-1}\Phi_k\otimes A_k)\boldsymbol U^{k+1}=(B^{-1}\otimes I_x)(\boldsymbol b+F(\boldsymbol U^k))-(B^{-1}\Phi_k\otimes A_k)\boldsymbol U^k$，
于是只需对角化 $B^{-1}\Phi_k=V\operatorname{diag}(\lambda_1,\ldots,\lambda_{N_t})V^{-1}$，第二步相应变成 $(I_x-\lambda_nA_k)\boldsymbol u_n^b=\boldsymbol u_n^a$。目前尚无 $B^{-1}\Phi_k$ 一般可对角化的理论，但实验中该矩阵通常可对角化且 $V$ 良态。

$\phi_n$ 的计算含矩阵—矩阵乘法，代价较高，因此不宜每轮 Newton 都重算。实践建议只在迭代前**离线**用**粗空间模型**（如粗网格半离散 ODE）算一次缩放因子。Figure 3.14 的精细网格为 $\Delta x=1/200$，$\{\phi_n\}$ 由 $\Delta X=1/20$ 的粗梯形模型离线得到。

![原论文 Figure 3.14：平均 Jacobian 与 NKA 准 Newton 的 Burgers 收敛比较](assets/papers/time-parallelization/source-figures/figure-3-14.svg)

Figure 3.14 用周期 Burgers 方程比较两种准 Newton（平均 Jacobian vs NKA），两幅图分别对应两种黏性，各含 $T=0.7$ 与 $T=1.3$。NKA 在所有设置下都更快，长窗口 $T=1.3$ 的收益最明显——因为窗口越长、Jacobian 幅值沿时间的变化越大，逐层缩放 $\Phi_k$ 相对单一 $A_k$ 的信息增益也越大。

> [!tip] 本站洞见
> 从 $I_t\otimes A_k$ 到 $\Phi_k\otimes A_k$ 是一次「秩一升级」：前者假设所有时间层共享同一 Jacobian，后者承认它们方向相同但幅值不同。NKA 的巧妙之处在于把「保持可对角化」这一硬约束（必须是 Kronecker 结构）与「尽量逼近真实 Jacobian」这一软目标（Frobenius 最优）统一在一个显式闭式里——$\phi_n$ 就是 $\nabla f_n$ 在 $A_k$ 上的投影系数。这提示了一条通用配方：当精确算子破坏并行结构时，不必退回单一近似，而可在**保结构子空间**内做最优投影，用极小的额外离线成本换取显著的收敛加速。

## 公式与图表覆盖核对

| 原文项目                                           | 论文小节 | 覆盖状态                                                                 |
| -------------------------------------------------- | -------- | ------------------------------------------------------------------------ |
| (3.22)–(3.25)                                      | 3.5.1    | 变步长后向 Euler、全时间矩阵、Kronecker 分解、三步求解与并行机制         |
| (3.26)–(3.29), Theorem 3.5                         | 3.5.1    | 几何网格、Jordan 块诱因、两类误差、平衡参数、Toeplitz 特征向量与证明路线 |
| Figures 3.9–3.10                                   | 3.5.1    | 全部原图、V 形参数扫描和并行宽度阈值                                     |
| (3.30)–(3.37), Theorem 3.6, Figure 3.11, Table 3.1 | 3.5.1    | 二阶转写、梯形归一化、误差平衡、$r_1/r_2$ 一致界和条件数平台             |
| (3.38)–(3.41), Theorem 3.7, Figure 3.12            | 3.5.1    | BVM、Dahlquist 壁垒、一致二阶、$O(N_t^2)$ 条件数与消元全推导             |
| (3.42)–(3.46), Figure 3.13, Table 3.2              | 3.5.1    | 非线性 Newton、结构破坏、平均 Jacobian、并行求解与成本对比               |
| (3.47)–(3.48), Figure 3.14                         | 3.5.1    | NKA、Frobenius 最优投影、离线粗模型和收敛改善                            |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 3.5–3.5.1, pp. 415–431.
