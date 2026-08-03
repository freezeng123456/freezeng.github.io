---
title: 4.1–4.2：历史脉络与 Parareal
description: 从多重打靶、虚拟控制到粗细传播子，逐式精读 Parareal 的线性/非线性收敛界、0.3 抛物因子与弱扩散/双曲失效机理
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-1-parareal
tags:
  - 时间并行
  - Parareal
  - 抛物方程
---

> [!note] 阅读范围
> 本页对应论文 Sections 4、4.1 和 4.2（pp. 443–452），覆盖公式 (4.1)–(4.9)、Theorems 4.1–4.4、Remark 4.1 和 Figures 4.1–4.5。公式后的推理按论文证明顺序展开，并在关键处补充“为什么”与参考文献。

## 4.1 历史发展

### 为何另设“面向抛物问题的方法”

第二章说明，抛物方程除极低频外具有较强的时间局部性；双曲方程的全部频率都能跨越长时间传播。第三章的方法在时间上同时处理长程耦合，因此也常能处理抛物问题；它们在非线性场景中仍有明显限制：OSWR 的优化 Robin 参数难以确定，ParaExp 与两类 ParaDiag 的外层 Newton 迭代会随时间窗增长而变慢，甚至失效。

第四章的方法主动利用耗散带来的时间局部性。Parareal、PFASST、MGRiT 和 STMG 在扩散充分时对线性、非线性问题都有效；把它们直接搬到弱扩散或双曲问题上，收敛会持续变慢，最后可能发散。本页把这条“扩散越强越好、传播越强越坏”的主线，从更新公式一路推到收敛因子的每一步。

Parareal 的思想可追溯到多重打靶、波形松弛和 Nievergelt (1964) 的非迭代先驱工作。Lions et al. (2001) 独立提出现代算法。随后出现 PITA、PFASST、MGRiT、Parareal–ParaDiag 组合。另一条路线是时空多重网格：早期方法难以有效粗化时间；Gander and Neumüller (2016) 通过时间块 Jacobi 平滑重新建立了可扩展的 STMG。

> [!tip] 本站洞见
> 论文特别指出，Lions–Maday–Turinici (2001) 并不是在“多重打靶”框架里提出 Parareal 的，而是在**虚拟控制（virtual control）**的语境下得到同一格式。两条来源殊途同归说明了一件事：Parareal 的 $\mathcal F+\mathcal G-\mathcal G$ 结构不是某一种特定推导的偶然产物，而是“用便宜模型做预测、用昂贵模型做校正”这一思想的必然形态。理解它最省力的路径，是把它读成多重打靶 Newton 迭代中 Jacobian 的有限差分近似（Gander and Vandewalle 2007）。

## 4.2 Parareal

### 更新公式与两层时间网格

将 $[0,T]$ 分成 $0=T_0<T_1<\cdots<T_{N_t}=T$。$\mathcal F$ 是昂贵而准确的细传播子，$\mathcal G$ 是便宜的粗传播子。Parareal 从界面初值猜测 $\boldsymbol u_n^0$ 出发：

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k). \tag{4.1}
$$

同一轮中，所有 $\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)$ 可并行；带新迭代指标 $k+1$ 的粗传播必须沿 $n$ 顺序推进。最后两项形成粗网格上的预测–校正，也可由多重打靶 Newton 的有限差分 Jacobian 解释。

这里有必要把每一项的角色讲透。$\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)$ 用小步长 $\Delta t$、以 $\boldsymbol u_n^k$ 为初值把区间 $[T_n,T_{n+1}]$ 积到底，代价高但精度高；$\mathcal G$ 用大步长 $\Delta T$（或更简单的模型）完成同样区间，便宜但粗糙。把 (4.1) 重排成

$$
\boldsymbol u_{n+1}^{k+1}
=\underbrace{\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})}_{\text{预测(串行、便宜)}}
+\underbrace{\big[\mathcal F(T_n,T_{n+1},\boldsymbol u_n^{k})-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k})\big]}_{\text{校正(并行、用上一轮的值)}},
$$

就能看清：粗传播沿时间“扫”出一条串行的预测轨迹，而昂贵的细–粗之差只依赖上一轮已知的 $\boldsymbol u_n^k$，因而所有 $n$ 上可以同时算。收敛时 $\boldsymbol u_n^{k+1}\to\boldsymbol u_n^k$，括号里两个 $\mathcal G$ 相消，界面上剩下的正是 $\mathcal F$ 的顺序解——这解释了为什么不动点恰好是串行细解，而与粗传播 $\mathcal G$ 的精度无关。

> [!tip] 本站洞见
> 之所以 $\mathcal F-\mathcal G$ 这一“差分”能当 Jacobian 用：多重打靶把界面连续性写成 $\boldsymbol u_{n+1}=\mathcal F(T_n,T_{n+1},\boldsymbol u_n)$ 的非线性方程组，Newton 步需要 $\partial\mathcal F/\partial\boldsymbol u_n$。直接算这个 Jacobian 和算 $\mathcal F$ 本身一样贵、且串行。Parareal 用“便宜传播子在相邻迭代上的一阶变化” $\mathcal G(\cdot,\boldsymbol u_n^{k+1})-\mathcal G(\cdot,\boldsymbol u_n^{k})$ 近似 Jacobian 与增量之积（Gander and Vandewalle 2007）。这既省去了昂贵的切线传播，又把 Newton 的**局部超线性**性质带进 Parareal——后面 (4.5a)、(4.6) 的阶乘/连乘因子正是这种超线性的量化。

![原论文 Figure 4.1：每个粗时间步包含 J 个细时间步](assets/papers/time-parallelization/source-figures/figure-4-1.svg)

论文主要讨论均匀网格，令 $\Delta T/\Delta t=J\ge2$；非均匀网格也可使用（Gander 2017；Maday and Mula 2020；Wu and Zhou 2024）。目标解是 $\mathcal F$ 顺序运行得到的离散解，Parareal 并不在固定 $\mathcal F$ 之外额外改变离散目标。换言之，$J\ge2$ 保证一个粗步至少套着两个细步，粗层才可能比细层“便宜”；同时也意味着后面出现的细传播模态必须写成**单步细稳定函数的 $J$ 次幂** $R_f^J(z/J)$，因为一个粗区间里细传播要连走 $J$ 步。这是全部线性分析的记账基础，也说明 Parareal 是**非侵入式**的：它只反复调用现成的 $\mathcal F$、$\mathcal G$，不改动它们内部的离散格式。

### Theorem 4.1：把误差化成逐模态 Toeplitz 迭代

考虑 $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$，设 $A=V_ADV_A^{-1}$，粗、细单步稳定函数为 $R_g$、$R_f$。一个粗区间内细走 $J$ 步，因此细传播模态为 $R_f^J(z/J)$，其中 $z=\Delta T\lambda(A)$。

若 $|R_g(z)|\le1$，则

$$
\max_{1\le n\le N_t}
\|V_A(\boldsymbol u_n^k-\boldsymbol u_n)\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}\|M^k(z)\|_\infty
\max_{1\le n\le N_t}
\|V_A(\boldsymbol u_n^0-\boldsymbol u_n)\|_\infty. \tag{4.2}
$$

这里 $\boldsymbol u_n$ 是细传播子的顺序解，

$$
M(z)=M_g^{-1}(z)[M_g(z)-M_f(z)],
$$

$$
M_g(z)=
\begin{bmatrix}
1\\-R_g(z)&1\\&\ddots&\ddots\\&&-R_g(z)&1
\end{bmatrix},
\qquad
M_f(z)=
\begin{bmatrix}
1\\-R_f^J(z/J)&1\\&\ddots&\ddots\\&&-R_f^J(z/J)&1
\end{bmatrix}. \tag{4.3}
$$

这里 $|R_g(z)|\le1$ 是唯一的结构性假设，也就是**粗传播稳定**。它保证 $M_g^{-1}(z)$（下三角、元素 $R_g^j$）不会因 $|R_g|>1$ 而把误差放大；抛物问题谱落在负实轴附近，后向 Euler 恰好满足 $|R_g|\le1$，因此这个假设几乎自动成立。

#### 证明链条

在矩阵层面，(4.1) 给出

$$
\boldsymbol u_{n+1}^{k+1}
=R_f^J(\Delta TA/J)\boldsymbol u_n^k
+R_g(\Delta TA)\boldsymbol u_n^{k+1}
-R_g(\Delta TA)\boldsymbol u_n^k.
$$

顺序细解满足同一恒等式。令 $\boldsymbol e_n^k=\boldsymbol u_n-\boldsymbol u_n^k$，相减得到

$$
\boldsymbol e_{n+1}^{k+1}
=R_g(\Delta TA)\boldsymbol e_n^{k+1}
+[R_f^J(\Delta TA/J)-R_g(\Delta TA)]\boldsymbol e_n^k.
$$

关键是**顺序细解也满足同一恒等式**：给它加上又减去 $R_g(\Delta TA)\boldsymbol u_n$，即得 $\boldsymbol u_{n+1}=R_f^J(\Delta TA/J)\boldsymbol u_n+R_g(\Delta TA)\boldsymbol u_n-R_g(\Delta TA)\boldsymbol u_n$。两式相减时，公共的 $R_f^J\boldsymbol u_n$ 与 $R_g\boldsymbol u_n$ 抵消，误差递推里**不含源项 $\boldsymbol g$**，也不含初值——因为 $T_0$ 处初值已知，故 $\boldsymbol e_0^k=0$ 对所有 $k$ 成立。这一步是整个分析能约化成齐次线性迭代的原因。

对角化 $A$ 后，每个模态 $\xi_n^k(z)$ 独立满足

$$
\xi_{n+1}^{k+1}(z)
=R_g(z)\xi_n^{k+1}(z)
+[R_f^J(z/J)-R_g(z)]\xi_n^k(z),
$$

并且

$$
\|V_A\boldsymbol e_n^k\|_\infty
=\max_{z\in\sigma(\Delta TA)}|\xi_n^k(z)|. \tag{4.4}
$$

因为 $R_g(\Delta TA)=V_AR_g(\Delta TD)V_A^{-1}$、$R_f^J(\Delta TA/J)=V_AR_f^J(\Delta TD/J)V_A^{-1}$，在特征基下矩阵递推**逐特征值解耦**，每个 $z=\Delta T\lambda$ 上只剩一个标量二重（对 $n$、对 $k$）Toeplitz 迭代。这正是抛物问题“按模态各自收敛”的数学写照：慢模态（$z$ 靠近 $0$）和快模态（$z$ 很负）用同一格式却有完全不同的收敛速度。

把所有 $n$ 叠成 $\boldsymbol\xi^k$，有

$$
M_g(z)\boldsymbol\xi^{k+1}
=[M_g(z)-M_f(z)]\boldsymbol\xi^k,
$$

于是 $\boldsymbol\xi^k=M^k(z)\boldsymbol\xi^0$，再取模态和时间节点的最大值便得 (4.2)。由此 $\|M^k(z)\|_\infty$ 就是把 Parareal 用到 Dahlquist 测试方程 $u'(t)=\lambda u(t)+g(t)$ 上的收敛因子；整份分析于是归结为“研究单个复数 $z$ 上的 Toeplitz 矩阵 $M(z)$ 的幂”。

> [!note] Remark 4.1：预条件视角
> $M(z)=I_t-M_g^{-1}(z)M_f(z)$。$M_fU=b$ 是细传播子的全时间系统，$M_g$ 是粗传播预条件器。残差中的每个细传播可以并行计算，粗预条件求解则沿时间顺序进行。这一视角直接导向 Section 4.5 的对角化粗校正。

把 Remark 4.1 写细一点：Parareal 一步就是对全时间系统 $M_f(z)U=b$ 做一次**预条件 Richardson/静态迭代**

$$
M_g(z)\,\Delta U^k=r^k:=b-M_f(z)U^k,\qquad U^{k+1}=U^k+\Delta U^k,
$$

其中残差分量

$$
r_n^k=b_n-\big(u_n^k-\mathcal F(T_{n-1},T_n,u_{n-1}^k)\big)
=b_n-\big(u_n^k-R_f^J(z/J)\,u_{n-1}^k\big)
$$

只依赖上一轮的 $u_{n-1}^k,u_n^k$，因此对所有 $n$ 可**同时**计算（并行来源）；而 $M_g$（下双对角）求解是沿时间的一次前代扫掠（串行、便宜）。这一读法把 Parareal 归入“昂贵算子 $M_f$ 做残差、便宜算子 $M_g$ 做预条件”的一般框架，也解释了为何 Section 4.5 只要把串行的 $M_g$ 换成可对角化的粗校正，就能进一步并行化。

严格下三角结构还说明：第 $k$ 轮后前 $k$ 个粗节点已与顺序细解一致，精确算术下至多 $N_t$ 轮终止。原因是 $M(z)$ 严格下三角、$k$ 次幂后前 $k$ 条对角线全为零，误差从时间起点开始一节一节被“锁定”为精确值。这是 Parareal 的**有限步终止**性质，后面 (4.5a) 与 (4.6) 的“连乘含零因子”正是它的定量版本。

### Theorem 4.2：短时间超线性与长时间线性

$M_g^{-1}$ 的下三角元素为 $R_g^j(z)$。因此

$$
M(z)=[R_f^J(z/J)-R_g(z)]\widetilde M(R_g(z)),
$$

其中 $\widetilde M(\beta)$ 的第一条次对角线为 $1$，再往下依次为 $\beta,\beta^2,\ldots$。对其幂取无穷范数得到两类界。

把 $M(z)$ 拆成“标量幅度 $\times$ 结构矩阵”是本定理的技术核心：所有关于 $z$ 的信息浓缩进标量 $R_f^J(z/J)-R_g(z)$，而所有关于时间窗长度 $N_t$、迭代次数 $k$ 的信息浓缩进纯结构矩阵 $\widetilde M(\beta)$ 的幂。于是

$$
\|M^k(z)\|_\infty=|R_f^J(z/J)-R_g(z)|^k\,\|\widetilde M^{\,k}(R_g(z))\|_\infty .
$$

Gander and Vandewalle (2007, Lemma 4.4) 给出 $\widetilde M$ 的幂范数的关键估计（此处按参考文献补入，便于看清 (4.5a)、(4.5b) 的来历）：

$$
\|\widetilde M^{\,k}(R_g(z))\|_\infty\le
\begin{cases}
\min\!\left\{\left(\dfrac{1-|R_g(z)|^{N_t-1}}{1-|R_g(z)|}\right)^{\!k},\ \dbinom{N_t-1}{k}\right\}, & |R_g(z)|<1,\\[2ex]
\dbinom{N_t-1}{k}, & |R_g(z)|=1.
\end{cases}
$$

取右端 $\min$ 里的二项式支路即得短时间界，取几何级数支路（令 $N_t\to\infty$）即得长时间界。两支路的分界正是“时间窗多长算短、多长算长”的判据。

短时间或较小 $N_t$ 下，

$$
\max_n\|\boldsymbol e_n^k\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}
\varrho_s(J,z,N_t,k)
\max_n\|\boldsymbol e_n^0\|_\infty,
$$

$$
\varrho_s(J,z,N_t,k)
=\frac{|R_g(z)-R_f^J(z/J)|^k}{k!}
\prod_{j=1}^{k}(N_t-j). \tag{4.5a}
$$

$k=N_t$ 时乘积含零，显式体现有限步收敛。这里 $\binom{N_t-1}{k}=\frac1{k!}\prod_{j=1}^{k}(N_t-j)$，$k=N_t$（或 $k\ge N_t$）时连乘出现因子 $N_t-N_t=0$，直接给出零界，与 Remark 4.1 后的下三角论证互为印证。分子的 $k$ 次幂配上 $1/k!$ 的阶乘衰减，是**超线性**的来源：只要粗细差 $|R_g-R_f^J|$ 有界，$k$ 增大时 $|R_g-R_f^J|^k/k!$ 会比任何几何级数更快趋零，因此短窗上误差“越迭代越陡”。

若 $|R_g(z)|<1$，还可用与 $N_t$ 无关的长时间界

$$
\max_n\|\boldsymbol e_n^k\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}\varrho_l^k(J,z)
\max_n\|\boldsymbol e_n^0\|_\infty,
$$

$$
\varrho_l(J,z)=
\frac{|R_g(z)-R_f^J(z/J)|}{1-|R_g(z)|}. \tag{4.5b}
$$

分子是粗细传播差，分母是粗传播的耗散裕量。粗传播接近单位模且又不能逼近细传播时，该因子会接近或超过 $1$。

> [!tip] 本站洞见
> (4.5a) 与 (4.5b) 不是两个独立结论，而是同一个 $\min$ 的两条支路，对应两种“记账方式”。短窗上 $N_t$ 小、二项式 $\binom{N_t-1}{k}$ 在 $k\to N_t$ 时被“撞到零”，于是**有限步、超线性**主导；长窗上 $N_t$ 大、二项式支路失效，几何级数 $\big(\tfrac{1}{1-|R_g|}\big)^k$ 接管，于是每轮固定收缩因子 $\varrho_l$、**线性**收敛。分母 $1-|R_g(z)|$ 是理解全章的钥匙：它是粗传播在该模态上的“耗散裕量”。扩散强时高模 $|R_g|\to0$、裕量接近 $1$，$\varrho_l$ 小；而一旦谱靠近虚轴（传播主导）使 $|R_g|\to1$，分母趋零、$\varrho_l\to\infty$，这就是后面弱扩散/双曲失效的统一解释。

![原论文 Figure 4.2：短时间超线性与长时间线性两种收敛阶段](assets/papers/time-parallelization/source-figures/figure-4-2.svg)

Figure 4.2 使用周期热方程、零源、$u_0(x)=\sin^2(2\pi x)$、$\Delta x=1/5$，粗细层都用后向 Euler，$J=10$。$T=0.02,N_t=6$ 时 $\varrho_s$ 准确描述超线性下降；更长区间中，误差按近似固定斜率下降，$\varrho_l$ 更合适。把空间网格细化到 $\Delta x=1/8$ 后，也会更早进入线性阶段。直觉上，$N_t$ 越大、可用的“撞零”步数越靠后，误差在到达有限步终止前已经先经历了漫长的线性阶段；而空间加密提高了最大 $|z|$，让谱更早覆盖到线性因子占优的区域，因此同样更早呈现线性收敛。

### Theorem 4.3：非线性超线性界

令 $\mathcal F$ 为精确传播，$\mathcal G$ 为 $p$ 阶方法，局部截断误差不超过 $C_3\Delta T^{p+1}$。假设

$$
\|\mathcal G(T_n,T_n+\Delta T,\boldsymbol v)
-\mathcal G(T_n,T_n+\Delta T,\boldsymbol w)\|
\le(1+C_2\Delta T)\|\boldsymbol v-\boldsymbol w\|,
$$

并且

$$
\mathcal F(T_n,T_{n+1},\boldsymbol v)
-\mathcal G(T_n,T_{n+1},\boldsymbol v)
=c_{p+1}(\boldsymbol v)\Delta T^{p+1}
+c_{p+2}(\boldsymbol v)\Delta T^{p+2}+\cdots,
$$

各系数对 $\boldsymbol v$ 连续可微，则

$$
\|\boldsymbol u(T_n)-\boldsymbol u_n^k\|
\le
\frac{C_3\Delta T^{p+1}(C_1\Delta T^{p+1})^{k+1}}{(k+1)!}
(1+C_2\Delta T)^{n-k-1}
\prod_{j=0}^{k}(n-j). \tag{4.6}
$$

当 $k\ge n$ 时乘积出现零，前 $n$ 个界面已经精确。小 $\Delta T$ 下，每多一轮又带来一份 $\Delta T^{p+1}$，这就是非线性超线性阶段的来源。

这一非线性结果（Gander and Hairer 2008, Theorem 1；亦见 Gander and Lunet 2024, Theorem 2.6）用**生成函数**证明，把 (4.5a) 的三个特征在一般非线性 ODE 上一一对应：连乘 $\prod_{j=0}^k(n-j)$ 与阶乘 $1/(k+1)!$ 联手给出与线性情形同型的**有限步 + 超线性**；每轮新增的 $C_1\Delta T^{p+1}$ 来自 $\mathcal F-\mathcal G$ 的 $p$ 阶一致性（粗传播越准、$p$ 越高、$\Delta T$ 越小，则每轮增益越强）；Lipschitz 系数 $1+C_2\Delta T$ 的幂 $ (1+C_2\Delta T)^{n-k-1}$ 则是**误差沿时间传播的放大**——在有限时间 $n\Delta T\le T$ 内它被 $e^{C_2T}$ 一致控制，故不破坏收敛，这正对应线性情形要求 $|R_g|\le1$（粗传播不放大误差）的非线性推广。Hamilton 系统上更精细的后向误差分析见 Gander and Hairer (2014)。

> [!tip] 本站洞见
> (4.6) 里 $\Delta T^{p+1}$ 的**幂次随 $k$ 线性增长**（第 $k$ 轮量级约 $\Delta T^{(p+1)(k+1)}$）是理解“为什么 Parareal 在扩散问题上如此高效”的另一视角：只要 $\Delta T$ 落在使 $C_1\Delta T^{p+1}<1$ 的范围内，几轮迭代就能把误差压到细解精度以下，而这正是耗散把高模迅速抹平、使非线性局部行为“接近线性且温和”的问题所具备的条件。反过来，若解里存在长期不衰减的振荡，$C_1$（依赖 $\mathcal F-\mathcal G$ 的系数）与 $C_2$（Lipschitz 常数）都会变大，超线性的窗口随之收窄。

### Theorem 4.4：约 0.3 的抛物长时间因子

若粗传播 $\mathcal G$ 用后向 Euler，细传播 $\mathcal F$ 用 L-稳定 Runge–Kutta 方法，则存在 $J_{\min}=O(1)$，使

$$
\max_{z\in\mathbb R_-}\varrho_l(J,z)\approx0.3,
\qquad J\ge J_{\min}. \tag{4.7}
$$

这个常数不随 $T,N_t$ 增长。它来自负实谱、粗层的高频耗散和有限的粗细传播差。用 Radau IIA 等组合时最坏因子还可降到约 $0.068$。

这里值得把“为什么恰好稳定在 $\approx0.3$、且与 $T,N_t$ 无关”讲清楚。抛物离散矩阵 $A$ 负半定，谱落在 $z\in\mathbb R_-$。看 $\varrho_l(J,z)=|R_g(z)-R_f^J(z/J)|/(1-|R_g(z)|)$ 的两端：

- **低频端 $z\to0^-$：** 粗细稳定函数都满足一致性 $R\approx 1+z+O(z^2)$，故 $R_g(z)\approx R_f^J(z/J)$，分子 $\to0$ 的阶高于分母，$\varrho_l\to0$。慢模态几乎一轮就修正好。
- **高频端 $z\to-\infty$：** L-稳定意味着 $R_f(z)\to0$，从而 $R_f^J(z/J)\to0$；后向 Euler $R_g(z)=1/(1-z)\to0$。于是分子 $|R_g-R_f^J|\to0$、分母 $1-|R_g|\to1$，$\varrho_l\to0$。高模被粗层强耗散“就地清零”，无需精确匹配。
- **中频段：** 最坏值出现在两端之间的某个有限 $z$ 上，数值上约为 $0.3$，与谱是否延伸到 $-\infty$ 无关，故对任意长的 $T,N_t$ 都成立。

$J\ge J_{\min}=O(1)$ 的作用是保证 $R_f^J(z/J)$ 已经进入其 L-稳定的渐近区（每个细步的 $z/J$ 足够小、$J$ 次幂足够“压平”高频），从而 (4.7) 的上确界稳定在 $0.3$ 附近。当 $\mathcal F$ 取精确传播 $\mathcal F=\exp(\Delta TA)$ 时，$J\ge2$ 即可达到 $\approx0.3$。Theorem 4.4 的证明分情形散见于文献：$\mathcal F$ 为后向 Euler 见 Mathew, Sarkis and Schaerer (2010)；梯形规则/BDF2/两种 SDIRK 见 Wu (2015)、Wu and Zhou (2015)；一般 L-稳定 $\mathcal F$ 见 Yang, Yuan and Zhou (2023)；连续层面的更早根源见 Gander and Vandewalle (2007, Table 5.1)。

若细传播仅 A-稳定，例如梯形规则，高频不会随 $|z|\to\infty$ 消失。此时在有限谱区间 $[0,z_{\max}]$ 上仍有

$$
\max_{z\in[0,z_{\max}]}\varrho_l(J,z)\approx0.3,
\qquad
J\ge J_{\min}=O(\log_2 z_{\max}). \tag{4.8}
$$

细层需要更多小步，才能让其高频行为与物理耗散相符。这里的 WHY 在于：梯形规则 $R_f(z)=\dfrac{1+z/2}{1-z/2}$ 在 $z\to-\infty$ 时 $|R_f|\to1$（趋于 $-1$），高频既不被耗散又带来符号振荡，$R_f^J(z/J)$ 在大 $|z|$ 处不再趋零，于是分子 $|R_g-R_f^J|$ 保持 $O(1)$，$\varrho_l$ 无法被压到 $0.3$。补救办法是增大 $J$：把每个细步的自变量缩到 $z/J$，只要 $z/J$ 落回梯形规则表现良好的区域，$R_f^J$ 就能重新近似耗散行为，所需的 $J$ 随谱宽 $z_{\max}$ 对数增长 $O(\log_2 z_{\max})$（Wu and Zhou 2015，对梯形规则与四阶 Gauss–RK 证明）。这也解释了 A-稳定与 L-稳定的本质差别：L-稳定天然“杀高频”，A-稳定则要靠加密细步“借来”高频耗散。

论文还给出两个 SDIRK 方法：

$$
\begin{array}{c|cc}
\gamma&\gamma&0\\
1&1-\gamma&\gamma\\ \hline
&1-\gamma&\gamma
\end{array}
\quad \gamma=\frac{2-\sqrt2}{2}
\qquad\text{(SDIRK22)},
$$

$$
\begin{array}{c|cc}
\gamma&\gamma&0\\
1-\gamma&-1/\sqrt3&\gamma\\ \hline
&1/2&1/2
\end{array}
\quad \gamma=\frac{3+\sqrt3}{6}
\qquad\text{(SDIRK23)}. \tag{4.9}
$$

SDIRK22 的 (4.7) 在 $J_{\min}=2$ 成立，SDIRK23 需要 $J_{\min}=4$。记号 “SDIRK$sp$” 指 $s$ 级、$p$ 阶的单对角隐式 RK：SDIRK22 是 L-稳定，二级、二阶，能立刻满足 $\approx0.3$；SDIRK23（Wu and Zhou 2015 给出 $J_{\min}=4$）阶更高但需要更大的粗细比才把最坏因子拉回 $0.3$。这与 (4.8) 的直觉一致——细传播的**稳定性类型**与**粗细比 $J$** 必须一起判断，单看阶数不足以预测 Parareal 的收敛速度。

![原论文 Figure 4.3：不同细传播子与粗细比带来的 Parareal 收敛差异](assets/papers/time-parallelization/source-figures/figure-4-3.svg)

Figure 4.3 取周期热方程、$\Delta x=1/256$、$\Delta T=0.1$、$T=4$、扩散系数 $0.1$。三个面板从左到右对应 $J=2,10,50$。$J=2$ 时，梯形规则在约 $10^{-4}$ 处停滞，SDIRK23 也明显慢于 SDIRK22；$J=10$ 时两种 SDIRK 曲线已经接近，而梯形规则仍保留较慢的尾段；$J=50$ 时三种细传播都贴近 $0.3^k$ 参考斜率。图中变化验证了 Theorem 4.4 的限定条件：细传播子的稳定性类型和粗细比需要一起判断。梯形规则在小 $J$ 处“停滞”正是 (4.8) 的写照——它只是 A-稳定，高频未被耗散，Parareal 试图用后向 Euler（对高频给出“正确物理”的强耗散）去逼近一个高频不准的细解，二者矛盾导致误差平台；增大 $J$ 后细解高频行为回到耗散侧，平台消失。

### 扩散减弱后的退化

论文随后固定 $T=4$、$\Delta T=0.1$、$\Delta x=1/128$、$J=32$，粗层用后向 Euler，细层用 SDIRK22，考察周期对流扩散和 Burgers 方程。

![原论文 Figure 4.4：三种黏性下每个对流扩散特征值对应的长时间因子](assets/papers/time-parallelization/source-figures/figure-4-4.svg)

对流扩散的谱由负实轴附近逐渐向虚轴展开。$\nu$ 减小时，$\max\varrho_l$ 逼近 $1$，说明粗传播越来越难修正长期存在的传播模态。用 (4.5b) 的分母来读最直接：对流扩散算子的特征值 $\lambda\approx-\nu\kappa^2+\mathrm i c\kappa$（$\kappa$ 为波数），黏性 $\nu$ 提供实部（耗散）、对流 $c$ 提供虚部（传播）。$\nu$ 大时实部主导，$|R_g|$ 远离 $1$、耗散裕量 $1-|R_g|$ 充足，$\varrho_l$ 小；$\nu$ 减小时谱向虚轴靠拢，$|R_g|\to1$、分母 $\to0$，即使分子有界，$\varrho_l$ 也被推向 $1$。这正是 Theorem 4.2 洞见里“分母是钥匙”的实验确证。

![原论文 Figure 4.5：对流扩散与 Burgers 方程上黏性降低导致的 Parareal 退化](assets/papers/time-parallelization/source-figures/figure-4-5.svg)

Figure 4.5(a) 与谱因子预测一致。Burgers 方程缺少同样精确的模态分析，但 Figure 4.5(b) 呈现相同趋势；约在 $\nu\le10^{-3}$ 时常规迭代会发散。严格下三角结构保证的有限步性质仍在，只是需要的轮数失去实用价值。波动方程也通常不收敛（Gander and Vandewalle 2007 已指出；更细致的分析见 Gander and Lunet 2020a,b、Gander, Lunet and Pogoželskytė 2023a、Gander, Lunet, Ruprecht and Speck 2023b），这正是第四章方法适用范围的边界。

对双曲/传播问题失效的机理，论文用第二章 Figure 2.4 的图景解释：**任意小的高频分量都能在时空中传播任意远**。因此想让便宜的粗传播 $\mathcal G$ 在时空两方向都逼近细传播 $\mathcal F$ 的精度极其困难；一旦为了精度把 $\mathcal G$ 做得很准，粗校正本身就变得昂贵，Parareal 也就失去加速意义。

> [!tip] 本站洞见
> 把“为何 Parareal 对扩散极好、对传播失效”一句话讲透：**抛物问题里高频被迅速耗散、只有少数慢模态长期存活**，粗传播只需在这些被保留的慢模态上够准（低频端 $R_g\approx R_f^J$），而对已经衰亡的高模，粗层的强耗散把它们“就地清零”即可——分子与分母都利于收敛，于是得到与窗长无关的 $\approx0.3$。**传播/双曲问题里模态既不衰减又不断跨时空搬运相位**，粗传播的**相位（色散）误差**不随时间消失、反而逐轮累积；此时 $|R_g|\to1$ 使 (4.5b) 分母趋零，$\varrho_l\to1$ 乃至发散。半 Lagrange 等“相位优化”粗传播（见下）之所以是主流补救方向，正是因为它直接瞄准相位而非幅度。这与 Ruprecht (2018) 对 Parareal 波传播特性的分析结论一致：粗细传播的相速度不匹配是收敛失败的根源。

在 MGRiT 社区（MGRiT 是 Parareal 的多层推广，见 Section 4.4），大量工作致力于让 MGRiT 处理对流方程（Howse et al. 2019；De Sterck et al. 2021；De Sterck, Falgout, Krzysik and Schroder 2023b；De Sterck, Falgout and Krzysik 2023a 及其中文献）。核心想法是用**半 Lagrange 离散**设计一个相位优化的粗传播子：沿特征线追踪信息、直接对齐相速度，因此对**线性**对流方程表现良好；但半 Lagrange 是基于特征线的方法，非线性问题里特征线相互作用、难以直接实现，故非线性情形仍是**开放问题**。另一条思路由 Gander and Wu (2020) 提出，同样面向双曲问题，其优点是相对容易处理非线性，将在 Section 4.5 展开。

## 公式、定理与图表覆盖核对

| 原文项目                             | 论文小节 | 覆盖状态                                                             |
| ------------------------------------ | -------- | -------------------------------------------------------------------- |
| Section 4 引言与 4.1                 | 4、4.1   | 抛物局部性、四类算法、历史路线、虚拟控制来源                         |
| (4.1), Figure 4.1                    | 4.2      | 更新式、预测–校正、并行细传播、顺序粗校正、双层时间网格、非侵入性    |
| (4.2)–(4.4), Theorem 4.1, Remark 4.1 | 4.2      | 模态化、Toeplitz 误差矩阵、完整证明、预条件/残差并行解释、有限步终止 |
| (4.5), Theorem 4.2, Figure 4.2       | 4.2      | Lemma 4.4 范数界、超线性与长时间线性、分母耗散裕量、两种实验阶段     |
| (4.6), Theorem 4.3                   | 4.2      | 非线性假设、生成函数、误差界三因子解读与有限步含义                   |
| (4.7)–(4.9), Theorem 4.4, Figure 4.3 | 4.2      | L/A 稳定细传播、0.3 与 0.068 因子、$J_{\min}$、SDIRK 对比            |
| Figures 4.4–4.5                      | 4.2      | 对流扩散谱、Burgers 实验、弱扩散/双曲失效机理、半 Lagrange 补救      |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 4–4.2, pp. 443–452.
