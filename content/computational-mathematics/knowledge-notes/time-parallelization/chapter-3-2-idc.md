---
title: 3.3：并行积分延迟校正
description: 从积分残差逐步推导 IDC 更新，解释延迟校正为何升阶，并完整说明 PIDC、RIDC 的并行调度与正则性限制
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-2-idc
tags:
  - 时间并行
  - IDC
  - 高阶时间积分
---

> [!note] 阅读范围
> 本页对应论文 Section 3.3（pp. 405–412），覆盖公式 (3.5)–(3.12)、Theorem 3.3 和 Figures 3.4–3.6。推导保留每个代数步骤，PIDC 与 RIDC 的并行调度分别说明。IDC 由 Dutt、Greengard 与 Rokhlin (2000) 提出；PIDC 由 Guibert 与 Tromeur-Dervout (2007) 引入；RIDC 由 Christlieb、Macdonald 与 Ong (2010) 提出。

## 3.3 时间并行 IDC

积分延迟校正（integral deferred correction，IDC）由 Dutt、Greengard 与 Rokhlin (2000) 提出，用一轮又一轮的迭代校正把低阶时间积分器逐步提升为高阶方法。原始 IDC 在时间上是串行的；本节介绍它的两种时间并行改造：Guibert 与 Tromeur-Dervout (2007) 的流水线 IDC（PIDC）和 Christlieb 等人 (2010) 的修正型 IDC（RIDC）。两者的并行结构与原始 IDC 差别很大，因此要先讲清楚 IDC 本身如何工作。

> [!tip] 本站洞见
> 延迟校正的思想更早可追溯到 Böhmer 与 Stetter (1984) 对离散化误差的迭代修正，Dutt 等人 (2000) 把它重新表述为在**积分方程**上做残差校正，从而得到数值稳定、易于升阶的谱延迟校正（SDC）。理解 IDC 的关键，是把它看作 Picard 型不动点迭代：每一轮都用上一轮的轨道估计一次残差，再解一个结构相同、但右端更精确的误差方程。

### IDC 从哪里开始

考虑非线性初值问题。与其直接对微分方程做时间步进，IDC 先把 ODE 写成等价的积分形式：

$$
\boldsymbol u(t)=\boldsymbol u_0+
\int_0^t\boldsymbol f(\boldsymbol u(\tau),\tau)\,d\tau,
\qquad t\in(0,T]. \tag{3.5}
$$

选积分形式而非微分形式有两个原因。其一，积分算子是**光滑化**的：对 $\boldsymbol f$ 求积分比求导更稳定，误差不会被放大。其二，一旦有了轨道近似，就能用高阶求积公式精确地评估右端积分，这正是后面升阶的来源。

先给出一个粗近似 $\widetilde{\boldsymbol u}(t)$。它可以是常函数 $\boldsymbol u_0$（即 $\widetilde{\boldsymbol u}(t)\equiv\boldsymbol u_0$），也可以来自一次低精度积分。定义真实误差和积分残差

$$
\boldsymbol e(t)=\boldsymbol u(t)-\widetilde{\boldsymbol u}(t),
$$

$$
\boldsymbol r(t)=\boldsymbol u_0+
\int_0^t\boldsymbol f(\widetilde{\boldsymbol u}(\tau),\tau)\,d\tau
-\widetilde{\boldsymbol u}(t). \tag{3.6}
$$

残差 $\boldsymbol r$ 衡量当前轨道 $\widetilde{\boldsymbol u}$ 违反积分方程 (3.5) 的程度：若把 $\widetilde{\boldsymbol u}$ 代回积分方程右端得到的结果与 $\widetilde{\boldsymbol u}$ 本身相差多少。当 $\widetilde{\boldsymbol u}$ 恰好是精确解时 $\boldsymbol r\equiv 0$；因此 $\boldsymbol r$ 越小，轨道越接近真解。注意 $\boldsymbol r(0)=\boldsymbol u_0-\widetilde{\boldsymbol u}(0)=\boldsymbol 0$，且误差同样满足 $\boldsymbol e(0)=\boldsymbol 0$，两者的初值都被固定。

将 $\boldsymbol u=\widetilde{\boldsymbol u}+\boldsymbol e$ 代入 (3.5)，再加减 $\int_0^t\boldsymbol f(\widetilde{\boldsymbol u},\tau)\,d\tau$，得到

$$
\begin{aligned}
\boldsymbol e(t)
&=\boldsymbol u_0+
\int_0^t\boldsymbol f(\widetilde{\boldsymbol u}(\tau)+\boldsymbol e(\tau),\tau)\,d\tau
-\widetilde{\boldsymbol u}(t)\\
&=\boldsymbol r(t)+
\int_0^t\left[
\boldsymbol f(\widetilde{\boldsymbol u}(\tau)+\boldsymbol e(\tau),\tau)
-\boldsymbol f(\widetilde{\boldsymbol u}(\tau),\tau)
\right]d\tau. \tag{3.7}
\end{aligned}
$$

第一行只是把 $\boldsymbol u=\widetilde{\boldsymbol u}+\boldsymbol e$ 代进积分方程；第二行把 $\int_0^t\boldsymbol f(\widetilde{\boldsymbol u},\tau)\,d\tau$ 从积分里析出，正好凑成残差 $\boldsymbol r(t)$ 的定义 (3.6)，剩下的是被积函数关于 $\boldsymbol e$ 的**增量**。这一步把误差 $\boldsymbol e$ 显式地表示为“已知残差 $\boldsymbol r$”加上“一个关于 $\boldsymbol e$ 的积分方程”。

对时间求导，得到误差方程

$$
\boldsymbol e'(t)-\boldsymbol r'(t)
=\boldsymbol f(\widetilde{\boldsymbol u}(t)+\boldsymbol e(t),t)
-\boldsymbol f(\widetilde{\boldsymbol u}(t),t). \tag{3.8}
$$

(3.7) 是积分形式的误差方程，(3.8) 是它的微分形式。之所以要同时保留两种形式，是因为下一步的离散化会**区别对待**两项：残差项 $\boldsymbol r'$（等价于 $\boldsymbol r$ 的积分增量）用高阶求积公式处理，而误差方程的 $\boldsymbol f$ 增量项则可以继续用最简单的低阶时间步进器。

> [!tip] 本站洞见
> (3.8) 的右端是 $\boldsymbol f(\widetilde{\boldsymbol u}+\boldsymbol e,\cdot)-\boldsymbol f(\widetilde{\boldsymbol u},\cdot)$，当 $\boldsymbol e$ 已经很小时它约等于 $\partial_{\boldsymbol u}\boldsymbol f\cdot\boldsymbol e$，是一个关于误差的“线性化”方程；其非齐次项来自残差 $\boldsymbol r$。这解释了为什么误差方程可以用与原方程相同的廉价步进器来解——它的解本身就是小量，低阶步进器的**绝对**误差因此也小。升阶的全部功劳落在残差积分的高阶求积上，而不在步进器上。

### 离散校正公式

取节点

$$
0=t_0<t_1<\cdots<t_M=T,
\qquad \boldsymbol u_m^k\approx\boldsymbol u(t_m),
$$

其中上标 $k$ 是校正轮次（sweep），下标 $m$ 是节点编号；$\boldsymbol u_m^k$ 记第 $k$ 轮在节点 $t_m$ 上的近似。令 $\Delta t_m=t_{m+1}-t_m$。对 (3.8) 使用线性 $\theta$ 方法（即用参数 $\theta\in[0,1]$ 在区间两端对 $\boldsymbol f$ 加权），得到

$$
\begin{aligned}
\boldsymbol e_{m+1}-\boldsymbol e_m
={}&\boldsymbol r_{m+1}-\boldsymbol r_m\\
&+\Delta t_m(1-\theta)
\left[\boldsymbol f(\boldsymbol u_m^{k+1},t_m)
-\boldsymbol f(\boldsymbol u_m^k,t_m)\right]\\
&+\Delta t_m\theta
\left[\boldsymbol f(\boldsymbol u_{m+1}^{k+1},t_{m+1})
-\boldsymbol f(\boldsymbol u_{m+1}^k,t_{m+1})\right]. \tag{3.9}
\end{aligned}
$$

这里把新轨道记作 $\boldsymbol u^{k+1}=\boldsymbol u^k+\boldsymbol e$，把 (3.8) 右端的 $\widetilde{\boldsymbol u}+\boldsymbol e$ 与 $\widetilde{\boldsymbol u}$ 分别对应到 $\boldsymbol u^{k+1}$ 与 $\boldsymbol u^k$。线性 $\theta$ 方法对左端 $\boldsymbol e'$ 用差商 $(\boldsymbol e_{m+1}-\boldsymbol e_m)/\Delta t_m$，对右端在 $t_m$ 处取权重 $1-\theta$、在 $t_{m+1}$ 处取权重 $\theta$。取 $\theta=1$ 得后向 Euler，$\theta=0$ 得前向 Euler，$\theta=1/2$ 得梯形规则。

由残差定义 (3.6)，相邻两节点的残差增量为

$$
\boldsymbol r_{m+1}-\boldsymbol r_m
=\int_{t_m}^{t_{m+1}}\boldsymbol f(\boldsymbol u^k(\tau),\tau)\,d\tau
-\left(\boldsymbol u_{m+1}^k-\boldsymbol u_m^k\right).
$$

也就是说，残差增量等于旧轨道在小区间 $[t_m,t_{m+1}]$ 上的**精确积分**减去旧轨道自身的端点跳变。前者用高阶求积近似，后者是已知量，因此残差增量完全可算。局部积分由全组节点上的 Lagrange 插值求积：

$$
\int_{t_m}^{t_{m+1}}\boldsymbol f(\boldsymbol u^k(\tau),\tau)\,d\tau
\approx\sum_{j=1}^{M}\omega_{m,j}
\boldsymbol f(\boldsymbol u_j^k,t_j), \tag{3.10a}
$$

$$
\omega_{m,j}=
\int_{t_m}^{t_{m+1}}
\prod_{\substack{i=1\\i\ne j}}^{M}
\frac{\tau-t_i}{t_j-t_i}\,d\tau. \tag{3.10b}
$$

关键之处在于：求和 (3.10a) 的下标 $j$ 跑遍**全部** $M$ 个节点，而不是只用小区间 $[t_m,t_{m+1}]$ 的两个端点。权重 $\omega_{m,j}$ 由整段 $[0,T]$ 上的 Lagrange 基函数在小区间上积分得到，因而 (3.10a) 是一个可达 $M$ 阶的高阶求积公式。正是这一步让廉价的低阶步进器“看到”了整段区间的曲率信息。

再用 $\boldsymbol u_m^{k+1}=\boldsymbol u_m^k+\boldsymbol e_m$ 消去误差变量 $\boldsymbol e$：把 (3.9) 中的 $\boldsymbol e_{m+1}-\boldsymbol e_m$ 换成 $(\boldsymbol u_{m+1}^{k+1}-\boldsymbol u_{m+1}^k)-(\boldsymbol u_m^{k+1}-\boldsymbol u_m^k)$，并把残差增量代入，$\boldsymbol u_{m+1}^k-\boldsymbol u_m^k$ 恰好与残差里的同名项抵消。逐节点校正公式为

$$
\begin{aligned}
\boldsymbol u_{m+1}^{k+1}
={}&\boldsymbol u_m^{k+1}
+\Delta t_m(1-\theta)
\left[\boldsymbol f(\boldsymbol u_m^{k+1},t_m)
-\boldsymbol f(\boldsymbol u_m^k,t_m)\right]\\
&+\Delta t_m\theta
\left[\boldsymbol f(\boldsymbol u_{m+1}^{k+1},t_{m+1})
-\boldsymbol f(\boldsymbol u_{m+1}^k,t_{m+1})\right]\\
&+\sum_{j=1}^{M}\omega_{m,j}
\boldsymbol f(\boldsymbol u_j^k,t_j),
\qquad m=0,\ldots,M-1. \tag{3.11}
\end{aligned}
$$

每个校正层 $k=0,1,\ldots,k_{\max}-1$ 都从左向右扫过 $m=0,1,\ldots,M-1$。$\theta=1$ 给出后向 Euler 型校正，$\theta=1/2$ 给出梯形型校正。注意 (3.11) 中只有 $\boldsymbol u_{m+1}^{k+1}$（当 $\theta\ne 0$ 时出现在 $\boldsymbol f$ 里）是未知量，其余全部来自已算好的旧轨道或已推进的新轨道端点，因此每一步只需解一个与基础步进器同样规模的（可能非线性的）代数方程。

> [!tip] 本站洞见（公式 (3.11) 的三项分工）
> 第一行 $\boldsymbol u_m^{k+1}$ 沿用简单积分器向前推进新解；中间两项以 $\theta$ 权重修正新旧轨道在局部端点上的动力学差 $\boldsymbol f(\boldsymbol u^{k+1},\cdot)-\boldsymbol f(\boldsymbol u^k,\cdot)$；最后的求积项 $\sum_j\omega_{m,j}\boldsymbol f(\boldsymbol u_j^k,t_j)$ 注入旧轨道在整个节点集上的高阶积分信息。校正能够升阶，原因就在于最后一项比基础步进器看到了更完整的时间区间——它把一个 $M$ 阶精确的积分作为“目标”，而前两项负责把新解推向这个目标。若新旧轨道相同，中间两项为零，公式退化为“低阶步进 + 高阶求积残差”的纯延迟校正。

### Theorem 3.3：每轮能提高多少阶

若基础积分器阶数为 $p$，使用 $M$ 个等距节点，则第 $k$ 次校正后的误差阶为

$$
O\!\left(\Delta t^{\min\{M,(k+1)p\}}\right).
$$

其含义可分两部分解读。$\;(k+1)p$ 项表明：初始猜测（$k=0$ 之前）达到阶 $p$，此后每完成一轮校正就再增加 $p$ 阶——这来自误差方程本身也用 $p$ 阶步进器求解，每轮把当前的**主导误差项**再消去一次。$\min\{M,\cdot\}$ 则是求积精度的天花板：(3.10) 用 $M$ 个节点的 Lagrange 插值，其求积误差为 $O(\Delta t^{M})$，无论再迭代多少轮，整体精度都不可能超过求积公式本身的阶。因此后向 Euler（$p=1$）每轮升 $1$ 阶、梯形规则（$p=2$）每轮升 $2$ 阶，直到饱和于 $M$。

Dutt 等人 (2000) 的原始 IDC 采用 Gauss 型节点，使求积上限更高。例如 Gauss–Lobatto 节点可达到 $2J-1$ 阶（$J$ 为节点数）。这类以谱精度求积的 IDC 通常称为谱延迟校正（spectral deferred correction，SDC），也是 PFASST 算法的核心部件（见 4.3 节；PFASST 由 Emmett 与 Minion (2012) 提出，其思想源于 Minion (2010) 用一步 SDC 替代 Parareal 的细解算子）。

长时间区间不适合用一个高阶多项式整体近似——单个高次多项式在大区间上会出现剧烈振荡（Runge 现象），求积精度反而下降。论文将 $[0,T]$ 划分为窗口

$$
I_n=[T_{n-1},T_n],
\qquad n=1,\ldots,N_t,
$$

其中 $T_0=0$、$T_{N_t}=T$。窗口足够小时，低次多项式就能给出精确求积。普通 IDC 先把 $I_n$ 的全部校正完成，再把末值传给 $I_{n+1}$。这一过程**完全串行**：$I_{n+1}$ 的初值在 $I_n$ 算完之前是未知的，所以第 $n+1$ 个窗口必须等待第 $n$ 个窗口结束；而且窗口内部沿 $m$ 的节点更新也是逐步进行的。正是这种双重串行依赖，促使人们寻找 PIDC 与 RIDC 两种并行改造。

## 3.3.1 流水线 IDC（PIDC）

### 按窗口铺开的流水线

第一种并行 IDC 即 PIDC，由 Guibert 与 Tromeur-Dervout (2007) 提出。它把 IDC 组织成流水线，其核心是一个适用于任意时间演化计算的简单观察，早在 Womble (1990) 就已提出：只要 $I_n$ 上算出一个**初步的**末端初值，下游窗口 $I_{n+1}=[T_n,T_{n+1}]$ 就能立即启动，不必等 $I_n$ 完全收敛。

具体地，$I_n$ 完成第一遍 sweep 后，就在 $t=T_n$ 得到一个粗末值 $\boldsymbol u_{n,M}^{1}$（即窗口 $I_n$ 一遍 sweep 后最右端的解）。$I_{n+1}$ 可以立即用 $\boldsymbol u_{n,M}^{1}$ 作初值启动它自己的第一遍 sweep，与此同时 $I_n$ 继续做第二遍。这一步之后又可以在 $I_{n+2}$ 启动，同时 $I_{n+1}$、$I_n$ 各自推进，三遍 sweep 并行。一般地，在 $I_n$ 上做第 $k$ 遍时，可以同时在 $I_{n+1}$ 做第 $k-1$ 遍，在 $I_{n+2}$ 做第 $k-2$ 遍，依此类推，直到 $I_{n+k-1}$ 上的第一遍。这样各窗口沿对角线错位，构成一条不断加宽的流水线。

![原论文 Figure 3.4：四个时间窗口上的 PIDC 流水线启动和稳态阶段](assets/papers/time-parallelization/source-figures/figure-3-4.svg)

Figure 3.4 取 $M=6$、$k_{\max}=4$。(a)–(d) 依次启动前四个窗口，流水线宽度从一条 sweep 增至四条，这是**填充（bootstrap）阶段**：只有先算出上游窗口的初步末值，下游才有初值可用。(e)、(f) 展示填满后的连续两个时刻，此时 $I_n$ 到 $I_{n+3}$ 上的四个不同校正层可以同时运行，这是**稳态阶段**，并发宽度达到 $k_{\max}=4$。黑色虚线记录已经完成的 sweep 历史，带红色圆点的红线表示当前并行运行的 sweep，黑色实线表示精确解。因此，这六幅图同时包含了启动成本和稳态并发宽度，不能只把 (e)、(f) 当成一般的时间推进示意——流水线的加速比受填充/排空两端的开销折损。

> [!tip] 本站洞见
> PIDC 的并行度不是免费的：稳态并发宽度最多为 $k_{\max}$，即校正轮数，而不是窗口数 $N_t$。当 $N_t\gg k_{\max}$ 时，绝大部分时间处于稳态，加速比趋近 $k_{\max}$；但每个窗口拿到的初值是不断被上游改写的“粗初值”，这带来下面的收敛性隐患。

PIDC 的每个窗口在 $n\ge 1$ 时都从一个粗糙且**不断变化**的初值开始 sweep。因为下游窗口的初值来自上游尚未收敛的解，随着校正轮次推进，所得解的精度**不保证单调下降**。这与串行 IDC（每个窗口拿到的是上游已收敛的精确初值）本质不同，是 PIDC 特有的现象，需要通过正则性实验来判断其影响。

### 周期对流扩散实验与离散矩阵 (3.12)

为检验上述现象，论文对对流扩散方程 (2.5) 应用 IDC 与 PIDC，取两个扩散系数 $\nu=1$ 与 $\nu=10^{-3}$。使用周期边界，中心差分离散，网格 $\Delta x=1/64$，得到线性 ODE 系统 $\boldsymbol u'(t)=A\boldsymbol u(t)$，其中半离散矩阵写成

$$
A=\frac{\nu}{\Delta x^2}A_{xx}
-\frac{1}{2\Delta x}A_x,
$$

第一项 $\tfrac{\nu}{\Delta x^2}A_{xx}\approx\nu\partial_{xx}$ 逼近扩散，第二项 $\tfrac{1}{2\Delta x}A_x\approx\partial_x$ 逼近对流。两个差分矩阵为

$$
A_{xx}=
\begin{bmatrix}
-2&1&&&1\\
1&-2&1&&\\
&\ddots&\ddots&\ddots&\\
&&1&-2&1\\
1&&&1&-2
\end{bmatrix},
\qquad
A_x=
\begin{bmatrix}
0&1&&&-1\\
-1&0&1&&\\
&\ddots&\ddots&\ddots&\\
&&-1&0&1\\
1&&&-1&0
\end{bmatrix}. \tag{3.12}
$$

两个矩阵右上角与左下角的额外元素来自**周期边界**：第一格与最后一格互为邻居，故循环补上环绕项。$A_{xx}$ 对称（扩散是耗散的），$A_x$ 反对称（对流是保守的）。

取 $T=3$、窗口宽度 $\Delta T=1/10$（故共 $30$ 个窗口）、$M=5$，基础积分器为后向 Euler。第 $n$ 个窗口、第 $k$ 遍后的误差定义为

$$
\operatorname{err}_n^k=
\frac{\max_m\lVert\boldsymbol u_{\mathrm{ref}}^{n,m}
-\boldsymbol u_k^{n,m}\rVert_\infty}
{\max_{n,m}\lVert\boldsymbol u_{\mathrm{ref}}^{n,m}\rVert_\infty}.
$$

分子取窗口 $n$ 内所有节点上数值解与参考解之差的最大范数，分母用全局参考解的最大范数归一化，因此这是一个**逐窗口的相对误差**。参考解 $\boldsymbol u_{\mathrm{ref}}$ 由 MATLAB 内置 ODE45 计算，相对和绝对容差均设为 $10^{-13}$，可视为“精确解”。下一个窗口的初始轨道固定为上一窗口第一遍 sweep 的末值：$\boldsymbol u_{n+1,m}^0\equiv\boldsymbol u_{n,M}^1$（对所有 $m=0,1,\ldots,M$）——这正是 PIDC 流水线中下游拿到的粗初值。

![原论文 Figure 3.5：IDC 与 PIDC 在不同正则性和黏性下的逐窗口误差](assets/papers/time-parallelization/source-figures/figure-3-5.svg)

Figure 3.5 的四个面板形成两组对照。横轴是奇数编号的时间窗口（$1,3,\ldots,29$），纵轴是相对误差（对数刻度），每个面板同时给出初始误差、第一遍和第二遍 sweep 之后的 IDC/PIDC 误差；(a) 的图例适用于其余面板。

- (a)、(b) 使用 $\sigma=1000$ 的尖锐源项 $g(x,t)$（见 (2.4)），近似 $\delta$ 函数型源，解的时间正则性不足；(c)、(d) 使用 $\sigma=5$ 的平滑源项，解足够正则。
- (a)、(c) 取大黏性 $\nu=1$；(b)、(d) 取小黏性 $\nu=10^{-3}$。

逐面板结论：

- **(a) 低正则、强扩散**：IDC 与 PIDC 表现相近，但第一遍校正之后误差**不再下降**——解不够正则，高阶近似无从发挥。
- **(b) 低正则、弱扩散**：黏性减小后，连第一遍 IDC 的改善都远不如 (a)，再迭代也帮助不大，PIDC 同样如此。
- **(c) 高正则、强扩散**：解足够光滑，IDC 与 PIDC 的**第二遍都继续降低误差**，且 PIDC 与串行 IDC 表现相当——这是升阶机制正常工作的理想情形。
- **(d) 高正则、弱扩散**：小黏性削弱了升阶效果；PIDC 第二遍明显落后于串行 IDC，因为弱扩散下算子接近双曲，误差沿特征线传播而非被耗散抹平，粗初值的污染更难被后续校正修复。

论文据此判断：对双曲型问题，若解的正则性不足，PIDC 并不适合用于 PinT（parallel-in-time）计算。升阶方法的收益本质上依赖解的光滑性。

## 3.3.2 修正型 IDC（RIDC）

### 滑动求积窗口

RIDC 由 Christlieb、Macdonald 与 Ong (2010) 提出，把并行粒度从“时间窗口”细化到“校正层”，用一个**滑动的 IDC 求积区间**实现更细粒度的并行。为此取等距节点的求积公式。

其调度按处理器展开：

1. **第一个处理器**用低阶时间步进器（如后向 Euler）连续向前推进，和普通 IDC 一样；但它算完前 $M$ 步后**不停下**，而是继续推进第 $M+1$、$M+2$ 步……源源不断地产生低阶轨道。
2. **第二个处理器**在第一处理器凑够前 $M$ 个值后启动第一层 IDC 校正。它算完第一个 IDC 区间（对第 $M$ 步做校正）后也不停下，而是把它的**求积区间和求积节点整体向右滑动一个细时间步**：原来用第一处理器的第 $1,2,\ldots,M$ 步，现在改用第 $2,3,\ldots,M+1$ 步作为 IDC 区间与求积节点，据此算出第 $M+1$ 步的校正；再滑到第 $3,4,\ldots,M+2$ 步算第 $M+2$ 步，依此类推。
3. **第三个处理器**在第二层数据足够后启动，采用同样的滑动规则，对第二处理器的输出再做一层校正。更高层处理器依此类推。

稳定阶段中，每个处理器负责一个固定的校正层，各层同时推进相邻的时间步，形成沿时间轴平移的流水线。与按窗口分块的 PIDC 不同，RIDC 的每一步都用一个“居中”的滑动求积窗口，因此需要保存若干历史值（窗口宽 $M$），并像所有流水线一样处理启动（fill）和排空（drain）阶段。

![原论文 Figure 3.6：IDC 与 RIDC 在相同对流扩散数据上的逐窗口误差](assets/papers/time-parallelization/source-figures/figure-3-6.svg)

Figure 3.6 沿用 Figure 3.5 的 PDE、网格、源项和黏性设置，只把 PIDC 换成 RIDC，初始条件取 $u(x,0)=0$。四个面板的对应关系保持不变：上排是 $\sigma=1000$（低正则），下排是 $\sigma=5$（高正则）；左列是 $\nu=1$（强扩散），右列是 $\nu=10^{-3}$（弱扩散）。因此 (a)–(d) 仍分别表示“低正则/强扩散”“低正则/弱扩散”“高正则/强扩散”和“高正则/弱扩散”。下排允许校正层继续降低误差，上排和右列则暴露出正则性不足与弱扩散带来的困难。

> [!tip] 本站洞见
> RIDC 只改变了**调度与内存结构**（滑动窗口、按层分配处理器），并没有改变 IDC 的数学内核，因此它继承了 IDC 对正则性的依赖：RIDC 仍是高阶近似技术，对低正则的双曲解，误差同样无法靠增加校正层持续下降。换言之，PIDC 与 RIDC 解决的是“如何并行”，而不是“如何在低正则下升阶”——后者是 IDC 家族共同的物理限制。

## 公式与图表覆盖核对

| 原文项目           | 论文小节 | 覆盖状态                                                          |
| ------------------ | -------- | ----------------------------------------------------------------- |
| (3.5)–(3.8)        | 3.3 导论 | 积分方程、残差、积分误差方程和微分误差方程，含初值与光滑化动机    |
| (3.9)–(3.11)       | 3.3 导论 | $\theta$ 离散、残差增量、Lagrange 求积权重和最终 IDC 更新         |
| Theorem 3.3        | 3.3 导论 | 基础阶数、每轮升 $p$ 阶、求积阶上限 $M$、SDC/Gauss–Lobatto/PFASST |
| 窗口划分           | 3.3 导论 | 长时间分窗、串行依赖来源                                          |
| Figure 3.4         | 3.3.1    | 完整 PIDC 启动与稳态流水图、bootstrap、$M=6/k_{\max}=4$           |
| (3.12), Figure 3.5 | 3.3.1    | 周期差分矩阵、全部参数、误差定义与四面板结论                      |
| Figure 3.6         | 3.3.2    | RIDC 滑动窗口逐处理器机制和正则性实验                             |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3.3, pp. 405–412.

参考文献补充（论文引用，作者与年份）：

- A. Dutt, L. Greengard and V. Rokhlin (2000), _Spectral deferred correction methods for ordinary differential equations_ — IDC/SDC 原始文献。
- D. Guibert and D. Tromeur-Dervout (2007), _Parallel deferred correction method for CFD problems_ — PIDC 提出。
- D. E. Womble (1990), _A time-stepping algorithm for parallel computers_, SIAM J. Sci. — PIDC 流水线思想来源。
- A. J. Christlieb, C. B. Macdonald and B. W. Ong (2010), _Parallel high-order integrators_ — RIDC 提出。
- M. L. Minion (2010), _A hybrid parareal spectral deferred corrections method_；M. Emmett and M. L. Minion (2012)，PFASST — SDC 与 Parareal 的结合。
- K. Böhmer and H. J. Stetter (1984) — 延迟/缺陷校正思想的更早来源。
