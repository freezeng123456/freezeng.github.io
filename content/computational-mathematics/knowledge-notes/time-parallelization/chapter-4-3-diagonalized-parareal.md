---
title: 4.5：基于对角化的 Parareal
description: 并行粗网格校正与区间内对角化粗传播子的两条完整路线，含 α-循环全时间求解、Theorem 4.7 阈值与 Theorem 4.8 抛物/双曲谱界的逐步推导
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-3-diagonalized-parareal
tags:
  - 时间并行
  - Parareal
  - ParaDiag
---

> [!note] 阅读范围
> 本页对应论文 Section 4.5（pp. 460–472），覆盖公式 (4.14)–(4.29)、Theorems 4.7–4.8、Remark 4.2 和 Figures 4.12–4.17。两种方法都使用对角化，作用位置和适用范围不同：第一种并行化跨粗点的 CGC，第二种在每个粗区间内构造可并行的特殊粗传播子。

## 4.5 基于对角化的 Parareal

本节介绍第三类 Parareal 变体，它把 ParaDiag（时间方向的对角化技术）注入到 Parareal 的粗网格校正（coarse grid correction, CGC）里。之所以要这样做，是因为标准 Parareal 的瓶颈恰恰在粗层：粗校正 (4.14) 必须沿 $N_t$ 个粗点顺序推进，这段顺序依赖限制了并行标度。ParaDiag 提供了一条出路——用 $\alpha$-循环矩阵把顺序耦合改写成可一次性对角化求解的全时间系统。文献里出现了两条本质不同的路线：

- 第一条（Wu 2018；Wu and Zhou 2019）用一个首尾耦合条件，把 CGC 自身改造成能被 ParaDiag 一次性求解的形式；
- 第二条（Gander and Wu 2020）设计一个紧贴细传播子的特殊粗传播子，它在每个大区间 $[T_n,T_{n+1}]$ 内借助 ParaDiag 低成本求解。

两条路线的机制、收敛性质和适用范围都不同，这也是本节反复强调的对照主线。

### 两条路线先分清

- **对角化 CGC（Section 4.5.1）**：修改 Parareal 跨 $N_t$ 个粗点的顺序粗校正。并行宽度来自粗时间点；粗、细传播仍可用不同积分器，收敛机制仍接近标准 Parareal，主要适合抛物问题。
- **对角化粗传播子（Section 4.5.2）**：保留标准 Parareal 粗校正的外形，在每个 $[T_n,T_{n+1}]$ 内，用 ParaDiag 同时处理 $J$ 个细步。粗、细传播使用同一个积分器和步长；该构造能传递长寿命频率，因此也能处理双曲问题。

> [!tip] 本站洞见
> 两条路线的分水岭在于「粗传播子是否忠实还原细传播子」。路线一保留了 Parareal 「便宜粗、昂贵细」的分层，只是把 CGC 的顺序求解并行化，因此它继承标准 Parareal 的收敛因子与适用范围（见 Theorem 4.7 的 $\rho_{\mathrm{new}}=\rho$）。路线二反其道而行，让粗传播子与细传播子使用完全相同的积分器和步长，仅靠首尾耦合换取并行度，于是它能像细传播子一样精确搬运所有频率分量，突破了标准 Parareal 在双曲问题上的困境。后文所有推导都可以放回这一对照来理解。

## 4.5.1 基于对角化的 CGC

### 从顺序 CGC 到首尾耦合

标准 Parareal 的粗网格校正为

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\boldsymbol b_{n+1}^k,
\qquad n=0,\ldots,N_t-1, \tag{4.14}
$$

$$
\boldsymbol b_{n+1}^k
=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k),
$$

并从 $\boldsymbol u_0^{k+1}=\boldsymbol u_0$ 顺序推进。这里的 $\boldsymbol b_{n+1}^k$ 是上一轮迭代已知的「细减粗」残差，$n=0$ 时的初值 $\boldsymbol u_0^{k+1}=\boldsymbol u_0$ 固定不变，正是这条初值 + 逐点递推链条迫使 CGC 只能串行执行。Wu (2018) 的核心想法是打破这条链：把固定初值换成一个首尾耦合条件

$$
\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0,
$$

让「头」显式依赖「尾」。这样一来递推系统首尾相连，成为一个循环结构，可以整体对角化，但代价是 $n=0$ 处的初值不再精确等于 $\boldsymbol u_0$。为了让收敛极限仍满足原初值问题，需要把残差里第 $0$ 个节点的自变量重新定义为真实初值，即

$$
\widetilde{\boldsymbol u}_n^k=
\begin{cases}
\boldsymbol u_0,&n=0,\\
\boldsymbol u_n^k,&n\ge1,
\end{cases}
$$

并使用

$$
\left\{
\begin{aligned}
\boldsymbol u_{n+1}^{k+1}
={}&\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol u}_n^k)
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k),\\
\boldsymbol u_0^{k+1}&=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0.
\end{aligned}
\right. \tag{4.15}
$$

$\widetilde{\boldsymbol u}_0^k=\boldsymbol u_0$ 这个修正是关键：若不作此替换，收敛不动点会偏离原 ODE 的解；替换后可验证收敛极限恰好满足离散初值问题。值得一提的是，(4.15) 使用的首尾条件比 ParaDiag-II 中「自然」的首尾条件 (3.55) 早出现一年，形式也略有差别——(3.55) 用的是差分型条件 $\boldsymbol u_0^{k+1}=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_0$，一年后由 Gander and Wu (2019) 提出，但在本节这个语境里两者效果相当。这一「早出现、后被更自然版本追认」的细节在 Remark 4.2 里会再次变得重要。

### 线性全时间系统与三步解法

先看线性 ODE $\boldsymbol u'=A\boldsymbol u$，$\boldsymbol u(0)=\boldsymbol u_0$，$t\in(0,T)$。粗层取后向 Euler，则 $\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})=(I_x-\Delta TA)^{-1}\boldsymbol u_n^{k+1}$，把它代回 (4.15)（两边同乘 $I_x-\Delta TA$）得到 $N_t$ 个线性方程

$$
\begin{aligned}
(I_x-\Delta TA)\boldsymbol u_1^{k+1}&=\boldsymbol u_0^{k+1}+(I_x-\Delta TA)\boldsymbol b_1^k,\\
(I_x-\Delta TA)\boldsymbol u_2^{k+1}&=\boldsymbol u_1^{k+1}+(I_x-\Delta TA)\boldsymbol b_2^k,\\
&\ \,\vdots\\
(I_x-\Delta TA)\boldsymbol u_{N_t}^{k+1}&=\boldsymbol u_{N_t-1}^{k+1}+(I_x-\Delta TA)\boldsymbol b_{N_t}^k,
\end{aligned}
$$

再加上 $\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0$。这些方程不能逐个求解，因为首尾条件把第一式的右端和最后一式的解绑在了一起。把首尾条件代入第一条粗方程，就得到全时间（all-at-once）系统

$$
(C_\alpha\otimes I_x-I_t\otimes\Delta TA)\boldsymbol U^{k+1}
=\boldsymbol g^k, \tag{4.16}
$$

其中 $\boldsymbol U^{k+1}=((\boldsymbol u_1^{k+1})^\top,\ldots,(\boldsymbol u_{N_t}^{k+1})^\top)^\top$，

$$
C_\alpha=
\begin{bmatrix}
1&&&-\alpha\\
-1&1\\
&\ddots&\ddots\\
&&-1&1
\end{bmatrix},
$$

$$
\boldsymbol g^k=
\begin{bmatrix}
\boldsymbol u_0+(I_x-\Delta TA)\boldsymbol b_1^k\\
(I_x-\Delta TA)\boldsymbol b_2^k\\
\vdots\\
(I_x-\Delta TA)\boldsymbol b_{N_t}^k
\end{bmatrix}.
$$

$C_\alpha$ 右上角的 $-\alpha$ 正是首尾耦合留下的循环元；它把双对角的顺序结构变成 $\alpha$-循环矩阵，从而可用离散 Fourier 变换对角化。设 $C_\alpha=F\,\mathrm{diag}(\lambda_1,\ldots,\lambda_{N_t})\,F^*$（$F$ 为离散 Fourier 矩阵，$\lambda_n$ 为特征值，见 (3.50)–(3.51)），则求解 (4.16) 分三步：

$$
\left\{
\begin{aligned}
\boldsymbol U^{a,k+1}&=(F\otimes I_x)\boldsymbol g^k,\\
(\lambda_nI_x-\Delta TA)\boldsymbol u_n^{b,k+1}
&=\boldsymbol u_n^{a,k+1},
&&n=1,\ldots,N_t,\\
\boldsymbol U^{k+1}&=(F^*\otimes I_x)\boldsymbol U^{b,k+1}.
\end{aligned}
\right. \tag{4.17}
$$

三步分别是：沿时间做 FFT（step a），在频域解 $N_t$ 个彼此独立、可完全并行的移位空间系统 $(\lambda_nI_x-\Delta TA)\,\cdot=\cdot$（step b），再逆 FFT 回到时间域（step c）。实际 $\alpha$-循环实现还包含把 $\alpha$-循环化为标准循环所需的对角缩放 $\mathrm{diag}(\alpha^{(n-1)/N_t})$，与 Section 3.5.2 完全一致。经过 (4.17)，新 CGC (4.15) 就能在全部 $N_t$ 个粗点上同时完成，顺序瓶颈被彻底消除。

> [!tip] 本站洞见
> (4.17) 的中间那步「$N_t$ 个独立移位空间求解」是整条路线的计算内核，后面非线性准 Newton 的 (4.19) 会一字不差地复用它。换句话说，$\alpha$-循环对角化把一个「$N_t\times N_t$ 时间耦合 + 空间」的大系统，解耦成「$N_t$ 个纯空间系统」外加两次 FFT——时间维的成本从 $O(N_t)$ 顺序步降到 $O(\log N_t)$ 的 FFT 深度，这正是 ParaDiag 相对顺序时间步的根本优势。

### Theorem 4.7：保持标准 Parareal 速度的阈值

参数 $\alpha$ 面临一对相互拉扯的约束。一方面，由 (4.15) 可见 $\alpha\to0$ 时首尾条件退化为 $\boldsymbol u_0^{k+1}=\boldsymbol u_0$，整个方法回到标准 CGC (4.14)，所以「足够小的 $\alpha$」应当收敛得和原始 Parareal 一样快。另一方面，$\alpha$ 越小，$\alpha$-循环矩阵 $C_\alpha$ 对角化引入的舍入误差越大（见 Section 3.5.2 的分析），在单精度、半精度等低工作精度下尤其危险。好在并不需要把 $\alpha$ 取得极小就能追平标准 CGC：存在一个明确的阈值。

**Theorem 4.7（Wu 2018）.** 设 $\rho$ 为标准 Parareal (4.14) 的收敛因子，$\rho_{\mathrm{new}}$ 为新变体 (4.15) 的收敛因子，粗传播子 $\mathcal G$ 为稳定的时间积分器，则

$$
\rho_{\mathrm{new}}=\rho,
\qquad
\alpha\le\frac{\rho}{1+\rho}. \tag{Theorem 4.7}
$$

也就是说，只要 $\alpha$ 不超过阈值 $\rho/(1+\rho)$，新方法的渐近收敛因子就与标准 Parareal 完全相同；再往下减小 $\alpha$ 并不会改善渐近速度，只会放大舍入风险。因此最优选择正是阈值本身 $\alpha=\rho/(1+\rho)$。实际中 $\rho=O(10^{-1})$，于是 $\alpha=\rho/(1+\rho)=O(10^{-1})$，此量级下 $\alpha$-循环对角化的舍入误差可忽略。

> [!tip] 本站洞见
> 需要区分该定理的「已证范围」与「数值观察范围」。Theorem 4.7 是针对线性问题 $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$、且 $A$ 具有**负实特征值**（典型抛物型）严格证明的。对于其他情形——例如 $A$ 具有**复特征值**——目前只有数值实验支持其结论成立，尚无证明。因此把「复特征值也适用」写进定理假设并不准确：它是一个数值观察到的猜测性推广，而非已证结论。工程上可以据此在对流占优等复谱问题上尝试同一阈值选取，但应意识到理论保证仅覆盖负实谱。

阈值的直观含义是一场「免费午餐」：把 CGC 从严格顺序变成完全并行，却不牺牲收敛速度。付出的仅是一个 $O(10^{-1})$ 的循环耦合项带来的可忽略舍入。这也是路线一在抛物问题上极具吸引力的原因。

![原论文 Figure 4.12：标准与对角化 CGC 在热方程和 ADE 上的误差](assets/papers/time-parallelization/source-figures/figure-4-12.svg)

实验使用周期边界、$u_0(x)=\sin(2\pi x)$、后向 Euler 粗层、SDIRK22 细层、$T=4$、$J=10$、$\Delta T=0.1$、$\Delta x=1/128$。对角化 CGC 各取三个 $\alpha$ 值以观察收敛因子如何随之变化。(a) 是热方程，标准 CGC 测得 $\rho\approx0.22$，由 Theorem 4.7 得阈值 $\rho/(1+\rho)\approx0.18$；因此 $\alpha=0.25,0.4$ 超过阈值、慢于标准 CGC，$\alpha=0.1$ 在阈值内、与标准曲线重合。(b) 是 $\nu=0.1$ 的 ADE，$\rho\approx0.39$、阈值约 $0.28$；这里 $\alpha=0.1,0.25$ 都在阈值内、能跟上标准 CGC，$\alpha=0.4$ 越过阈值、明显变慢。两个面板分别验证同一个阈值公式在两种谱结构下的位置，理论预测与实测高度吻合。

### 非线性全时间准 Newton

对非线性问题 $\boldsymbol u'=f(\boldsymbol u)$，$\boldsymbol u(0)=\boldsymbol u_0$，粗层继续用后向 Euler。仍定义

$$
\boldsymbol b_{n+1}^k
=\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol u}_n^k)
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k).
$$

把 (4.15) 的粗校正 $\boldsymbol u_{n+1}^{k+1}=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})+\boldsymbol b_{n+1}^k$ 用后向 Euler 展开，可整理成 $(\boldsymbol u_{n+1}^{k+1}-\boldsymbol b_{n+1}^k-\boldsymbol u_n^{k+1})/\Delta T=f(\boldsymbol u_{n+1}^{k+1}-\boldsymbol b_{n+1}^k)$；再与首尾条件 $\boldsymbol u_0^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_0$ 合并，得到非线性全时间系统

$$
(C_\alpha\otimes I_x)\boldsymbol U^{k+1}
-\Delta T F(\boldsymbol U^{k+1})=\boldsymbol g^k, \tag{4.18}
$$

其中 $\boldsymbol U^{k+1}$、$C_\alpha$ 与 (4.16) 一致，$F$ 的第 $n$ 个块为 $f(\boldsymbol u_n^{k+1}-\boldsymbol b_n^k)$，$\boldsymbol g^k$ 首块为 $\boldsymbol b_1^k+\boldsymbol u_0$、其余为 $\boldsymbol b_2^k,\ldots,\boldsymbol b_{N_t}^k$。对 (4.18) 采用与非线性 ParaDiag（Section 3.5.1）相同的准 Newton 迭代：

$$
P_\alpha^{k+1,l}\Delta\boldsymbol U^{k+1,l}
=\boldsymbol g^k-(C_\alpha\otimes I_x)\boldsymbol U^{k+1,l}
+\Delta TF(\boldsymbol U^{k+1,l}),
$$

$$
\boldsymbol U^{k+1,l+1}
=\boldsymbol U^{k+1,l}+\Delta\boldsymbol U^{k+1,l}, \tag{4.19a}
$$

$$
P_\alpha^{k+1,l}
=C_\alpha\otimes I_x-I_t\otimes\Delta TA^{k+1,l}, \tag{4.19b}
$$

其中 $A^{k+1,l}=\frac1J\sum_{j=1}^{J}\nabla f(\boldsymbol u_n^{k+1,l}-\boldsymbol b_n^k)$ 是各时间节点 Jacobian 的平均，$I_t\otimes A^{k+1,l}$ 以之近似完整块对角 Jacobian $\nabla F(\boldsymbol U^{k+1,l})$。

> [!tip] 本站洞见
> 这里用「平均 Jacobian」而非「逐点真实 Jacobian」是有意为之：平均后 $P_\alpha^{k+1,l}$ 恰好保持 (4.16) 的 $C_\alpha\otimes I_x-I_t\otimes(\cdot)$ 张量结构，从而增量 $\Delta\boldsymbol U^{k+1,l}$ 仍能用 (4.17) 的三步对角化并行求解。若换成逐点 Jacobian，各时间块不再共享同一空间矩阵，FFT 对角化就失效了。这是一种「以准 Newton 换可对角化结构」的典型折中。论文指出 Section 3.5.1 的最近 Kronecker 积近似（NKA）能给出比平均更精确的 $\nabla F$ 近似，但为简洁未展开。非线性收敛分析（Wu 2018, Section 4）表明：$\alpha$ 取得适当小时，收敛速度与标准 CGC 的 Parareal 一致，阈值机制照旧。

![原论文 Figure 4.13：两种黏性 Burgers 方程上的两类 CGC](assets/papers/time-parallelization/source-figures/figure-4-13.svg)

Figure 4.13 取 Burgers 方程，问题设置与离散参数同前面热与 ADE 实验，左、右面板分别对应 $\nu=1$ 与 $\nu=0.01$，每幅都比较 $\alpha=0.4,0.25,0.1$ 与标准 CGC。两种黏性下，$\alpha=0.4$ 都最慢，$\alpha=0.1$ 最接近标准曲线（弱扩散面板中甚至略快）。可见非线性情形里 $\alpha$ 对收敛率的影响与线性情形一致，仍保留 Figure 4.12 的阈值结构。

### Remark 4.2：MGRiT 需要一致的首尾条件

MGRiT (4.12) 的底层机制与 Parareal 一致，其 CGC 也可写成 (4.14) 的形式。但若把 (4.15) 的首尾条件 $\boldsymbol u_1^{k+1}=\alpha\boldsymbol u_{N_t}^{k+1}+\boldsymbol u_1$ 直接照搬到 MGRiT，则**对任意 $\alpha$ 都发散**。原因在于这个条件在收敛处并不自洽——它比 Gander and Wu (2019) 一年后给出的「自然」条件粗糙。收敛的一致条件应改用差分型

$$
\boldsymbol u_1^{k+1}
=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_1,
$$

其在收敛（$\boldsymbol u_{N_t}^{k+1}=\boldsymbol u_{N_t}^k$）时自动退化为精确关系 $\boldsymbol u_1^{k+1}=\boldsymbol u_1$，因而「收敛处一致」。Wu and Zhou (2019) 用它得到收敛的 MGRiT 变体

$$
\left\{
\begin{aligned}
\boldsymbol u_0^{k+1}&=\boldsymbol u_0,\\
\boldsymbol u_1^{k+1}
&=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_1,\\
\boldsymbol u_{n+1}^{k+1}
&=\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\widetilde{\boldsymbol b}_{n+1}^k,\quad n=1,\ldots,N_t-1.
\end{aligned}
\right. \tag{4.20}
$$

这里 $\widetilde{\boldsymbol b}_{n+1}^k=\mathcal F(T_n,T_{n+1},\widetilde{\boldsymbol s}_n^k)-\mathcal G(T_n,T_{n+1},\widetilde{\boldsymbol s}_n^k)$，$\widetilde{\boldsymbol s}_n^k=\mathcal F(T_{n-1},T_n,\widetilde{\boldsymbol u}_{n-1}^k)$，且 $\widetilde{\boldsymbol u}_n^k=\boldsymbol u_n$（$n=0,1$）、$\widetilde{\boldsymbol u}_n^k=\boldsymbol u_n^k$（$n\ge2$）。$\alpha$ 适当小时，该变体与原 MGRiT (4.12) 同速，Theorem 4.7 的阈值机制以类似方式适用。Parareal 本身也可改用这条一致的差分首尾条件 $\boldsymbol u_0^{k+1}=\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)+\boldsymbol u_0$ 替换 (4.15)，收敛因子与 Theorem 4.7 相同。

> [!tip] 本站洞见
> 这条 Remark 揭示了一个易被忽略的原则：把 ParaDiag 首尾条件迁移到不同迭代格式时，「收敛处一致性」是收敛的前提。Parareal 对不一致的 (4.15) 尚且容忍（因其误差传播结构不同），MGRiT 却会因此彻底发散——差分型条件 $\alpha(\boldsymbol u_{N_t}^{k+1}-\boldsymbol u_{N_t}^k)$ 之所以更「自然」，正是因为它在不动点处自动消失，不给收敛极限引入偏差。

## 4.5.2 基于对角化的粗传播子

Gander and Wu (2020) 提出了一个与 Section 4.5.1 本质不同的想法。它的关键创新是：粗、细传播子使用**同一个时间积分器和同一个步长**，只是粗传播子通过对角化来实现。因为粗传播子不再是「大步长、强耗散」的廉价近似，而是与细传播子逐频一致，它能在很长的时间上忠实搬运所有频率分量——这正是它能处理双曲问题的根源，也是与路线一最尖锐的对照。

### 细传播与首尾耦合粗传播

在每个大区间 $[T_n,T_{n+1}]$ 内使用线性-$\theta$ 方法，步长 $\Delta t=\Delta T/J$（$s$ 级 Runge–Kutta 的推广见 Gander and Wu 2020 附录）。细传播子 $\mathcal F(T_n,T_{n+1},\boldsymbol u_n)=\boldsymbol v_J$ 顺序执行 $J$ 个细步

$$
\boldsymbol v_{j+1}-\boldsymbol v_j
=\Delta t[\theta f(\boldsymbol v_{j+1})
+(1-\theta)f(\boldsymbol v_j)],
\quad j=0,\ldots,J-1,
\quad \boldsymbol v_0=\boldsymbol u_n. \tag{4.21}
$$

$\theta=1$ 为后向 Euler，$\theta=1/2$ 为梯形规则。特殊粗传播 $\mathcal F_\alpha^*$ 使用**完全相同**的差分式，只把顺序初值 $\boldsymbol v_0=\boldsymbol u_n$ 换成首尾耦合初值

$$
\boldsymbol v_0=\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n. \tag{4.22}
$$

正是这一步把区间内的 $J$ 个细步从「顺序链」变成「首尾相连的循环系统」，从而可由 ParaDiag 同时求解。注意与路线一的对照：路线一在 $N_t$ 个**粗点**上做首尾耦合，路线二在单个区间的 $J$ 个**细点**上做首尾耦合。

### 非线性全时间系统与准 Newton

令 $\boldsymbol V=(\boldsymbol v_1^\top,\ldots,\boldsymbol v_J^\top)^\top$，(4.21)–(4.22) 写成非线性全时间系统

$$
\underbrace{(C_\alpha\otimes I_x)\boldsymbol V
-\Delta tF(\boldsymbol V)}_{K(\boldsymbol V)}
=\boldsymbol b(\boldsymbol u_n), \tag{4.23}
$$

$$
\boldsymbol b(\boldsymbol u_n)
=((1-\alpha)\boldsymbol u_n^\top,0,\ldots,0)^\top. \tag{4.24}
$$

其中 $C_\alpha$ 与前同，$F$ 的首块因首尾耦合而同时含 $\theta f(\boldsymbol v_1)$ 与 $(1-\theta)f(\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n)$，其余块为 $\theta f(\boldsymbol v_j)+(1-\theta)f(\boldsymbol v_{j-1})$。准 Newton 更新为

$$
P_\alpha(\boldsymbol V^l)\Delta\boldsymbol V^l
=\boldsymbol b(\boldsymbol u_n)-K(\boldsymbol V^l),
\qquad
\boldsymbol V^{l+1}=\boldsymbol V^l+\Delta\boldsymbol V^l, \tag{4.25a}
$$

$$
P_\alpha(\boldsymbol V^l)
=C_\alpha\otimes I_x
-\Delta t\widetilde C_{\alpha,\theta}\otimes\overline{\nabla f}(\boldsymbol V^l), \tag{4.25b}
$$

其中

$$
\widetilde C_{\alpha,\theta}=
\begin{bmatrix}
\theta&&&(1-\theta)\alpha\\
1-\theta&\theta\\
&\ddots&\ddots\\
&&1-\theta&\theta
\end{bmatrix},
$$

$\overline{\nabla f}=\frac1J\big[\sum_{j=1}^{J-1}\nabla f(\boldsymbol v_j^l)+\nabla f(\alpha\boldsymbol v_J^l+(1-\alpha)\boldsymbol u_n)\big]$ 是 $J$ 个 Jacobian 块的平均。$P_\alpha(\boldsymbol V^l)$ 是块 $\alpha$-循环矩阵，作为真实 Jacobian $\nabla K(\boldsymbol V^l)=C_\alpha\otimes I_x-\Delta t(\widetilde C_{\theta,\alpha}\otimes I_x)\nabla F(\boldsymbol V^l)$ 的近似。

> [!tip] 本站洞见
> 与路线一的 (4.19b) 相比，这里多了一个 $\theta$-加权的循环矩阵 $\widetilde C_{\alpha,\theta}$。关键在于 $C_\alpha$ 与 $\widetilde C_{\alpha,\theta}$ 是**同一族循环矩阵**、可被同一个 Fourier 矩阵**同时对角化**，因此 (4.25) 的内层空间系统仍可完全并行——线性-$\theta$ 方法的隐式/显式加权并没有破坏可对角化结构。这解释了为何一定要选线性-$\theta$（或其 RK 推广）这类「循环友好」的积分器。

外层 Parareal 仍写成标准形式

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal F_\alpha^*(T_n,T_{n+1},\boldsymbol u_n^{k+1})
+\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
-\mathcal F_\alpha^*(T_n,T_{n+1},\boldsymbol u_n^k). \tag{4.26}
$$

注意 (4.26) 与标准 Parareal 的粗校正外形完全一样，只是把粗传播子换成了 $\mathcal F_\alpha^*$——这正是「保留 CGC 外形、只改粗传播子」的体现，与路线一「改 CGC 本身」形成互补。

### 线性系统、并行性与极限

$f(\boldsymbol u)=A\boldsymbol u$ 时，(4.23) 化为

$$
(C_\alpha\otimes I_x
-\widetilde C_{\theta,\alpha}\otimes\Delta tA)\boldsymbol V
=\boldsymbol b(\boldsymbol u_n), \tag{4.27}
$$

$$
\boldsymbol b(\boldsymbol u_n)
=([(I_x+\Delta t(1-\theta)A)(1-\alpha)\boldsymbol u_n]^\top,0,\ldots,0)^\top.
$$

粗传播只需取末点，$\mathcal F_\alpha^*=(H_J\otimes I_x)\boldsymbol V=\boldsymbol v_J$（$H_J=(0,\ldots,0,1)\in\mathbb R^{1\times J}$）。它等价于求解

$$
\left\{
\begin{aligned}
\boldsymbol v_{j+1}-\boldsymbol v_j
&=\Delta tA[\theta\boldsymbol v_{j+1}+(1-\theta)\boldsymbol v_j],\\
\boldsymbol v_0&=\alpha\boldsymbol v_J+(1-\alpha)\boldsymbol u_n^k.
\end{aligned}
\right. \tag{4.28}
$$

两个极端揭示了 $\alpha$ 的作用：

- $\alpha=0$ 时 (4.28) 的初值退化为 $\boldsymbol v_0=\boldsymbol u_n^k$，粗传播子就等于顺序细传播子 $\mathcal F$，于是 Parareal (4.26) **一轮即收敛**，但因为要顺序解 $J$ 步、**完全没有加速**。
- $0<\alpha<1$ 时，(4.28) 变成首尾耦合系统，可由对角化一次性求解，$J$ 个细点并行；若空间系统有足够计算资源，粗传播的墙钟成本约为顺序细传播的 $1/J$。

> [!tip] 本站洞见
> $\alpha$ 在这里是「精度—并行度」旋钮，含义与路线一不同。路线一的 $\alpha$ 主要权衡「循环耦合舍入 vs 保持标准收敛率」；路线二的 $\alpha$ 直接权衡「粗传播子对细传播子的逼真度（$\alpha\to0$ 完全逼真但无并行）vs 并行加速（$\alpha>0$ 换来 $1/J$ 成本但引入收敛因子）」。理解这一点，就能读懂 Theorem 4.8 里 $\rho$ 随 $\alpha$ 直接变化的原因。

### Theorem 4.8：抛物谱与双曲谱

**Theorem 4.8（Gander and Wu 2020）.** 对线性初值问题 $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$，$\boldsymbol u(0)=\boldsymbol u_0$，$A\in\mathbb C^{N_x\times N_x}$，设 $\{\boldsymbol u_n^k\}$ 为变体 (4.26) 第 $k$ 轮迭代、$\{\boldsymbol u_n\}$ 为收敛解。若 $\mathcal F$ 与 $\mathcal F_\alpha^*$ 都用稳定的单步 Runge–Kutta 方法，令

$$
e^k=\max_{1\le n\le N_t}\|\boldsymbol u_n-\boldsymbol u_n^k\|_\infty,
$$

则

$$
e^k\le\rho^ke^0,
\qquad
\rho=
\begin{cases}
\alpha,&\sigma(A)\subset\mathbb R_-,\\[4pt]
\dfrac{2\alpha N_t}{1+\alpha},&\sigma(A)\subset i\mathbb R.
\end{cases} \tag{4.29}
$$

两种谱给出截然不同的行为，这正是本节的核心结论：

- **负实谱（抛物型）**：$A$ 由半离散热方程 $A\approx\Delta$ 得到，$\sigma(A)\subset\mathbb R_-$，此时 $\rho=\alpha$，**与粗区间数 $N_t$ 无关**。选 $\alpha=O(10^{-2})$ 即可获得极快且鲁棒的收敛。
- **纯虚谱（双曲型）**：二阶波动方程 (2.7)、Schrödinger 方程等半离散后 $\sigma(A)\subset i\mathbb R$，此时 $\rho=2\alpha N_t/(1+\alpha)$，**随 $N_t$ 线性增长**。但这是一个上界，$\alpha N_t$ 较大时可能很松，并不必然意味着实际收敛变差——只要 $\alpha$ 足够小、$N_t$ 不太大，仍可保持快速收敛。

> [!tip] 本站洞见
> 这里出现了全节最关键的对照：为什么路线二能处理双曲、而路线一（Theorem 4.7，$\rho_{\mathrm{new}}=\rho$）只停留在抛物？根源在粗传播子的构造。路线一的粗传播子 $\mathcal G$ 是**大步长廉价积分器**，对高频/振荡分量有强耗散或相位误差，它继承了标准 Parareal 在双曲问题上「粗传播子无法忠实还原细传播子」的老毛病，所以其收敛因子只能等同于标准 Parareal（抛物友好、双曲困难）。路线二的粗传播子 $\mathcal F_\alpha^*$ 用**与细传播子相同的步长**，逐频精确搬运，唯一的近似来自 $\alpha$ 引入的首尾耦合扰动；因此双曲情形的收敛因子退化为温和的 $2\alpha N_t/(1+\alpha)$——退化程度可由 $\alpha$ 直接控制，而非像大步长粗传播那样不可控。一句话：路线二把「双曲失效」从结构性缺陷降级为一个可调参数问题。

![原论文 Figure 4.14：热方程上 rho=alpha 的锐利预测](assets/papers/time-parallelization/source-figures/figure-4-14.svg)

热方程采用齐次 Dirichlet 边界、$u_0=\sin^2(2\pi x)$、梯形规则、$\Delta T=1/12$、$J=10$、$\Delta x=1/100$。左、右面板分别取 $N_t=36$ 与 $72$，每幅都比较 $\alpha=10^{-1},10^{-2},10^{-3}$。实测虚线与理论点线（$\rho=\alpha$）几乎平行，$N_t$ 加倍没有改变由 $\rho=\alpha$ 决定的斜率——直接印证了负实谱因子与 $N_t$ 无关且上界锐利。

![原论文 Figure 4.15：波动方程上 alpha 与粗区间数的共同影响](assets/papers/time-parallelization/source-figures/figure-4-15.svg)

波动方程取周期边界、$u_0=\sin^2(2\pi x)$、$u_t(x,0)=0$，半离散后 $A\approx\Delta$ 的特征值全为纯虚，故按 (4.29) 收敛率随 $N_t$ 增长而恶化。(a) 固定 $\alpha=0.01$（相对较大），比较 $N_t=24,48,96$，区间数增加明显减慢，验证 $2\alpha N_t/(1+\alpha)$ 的线性趋势；(b) 固定 $\alpha=10^{-4}$（很小），比较 $N_t=24,48,96,960$，$N_t$ 从 24 增到 960 只多约两轮就达到离散误差量级 $\max\{\Delta t^2,\Delta x^2\}$。两个面板把 $\alpha N_t$ 的联合作用分离出来：$N_t$ 的负面影响可被小 $\alpha$ 有效压制。

![原论文 Figure 4.16：小 alpha Nt 时理论因子较锐利，大乘积时出现超线性](assets/papers/time-parallelization/source-figures/figure-4-16.svg)

与热方程不同，波动方程里 (4.29) 的因子是否锐利取决于乘积 $\alpha N_t$。Figure 4.16 考察三组 $(\alpha,N_t)$：只有 $\alpha=10^{-4},N_t=24$ 的小乘积组合紧贴 (4.29) 的点线上界；另两组实测曲线出现**超线性**下降，线性上界明显偏保守。这说明 $2\alpha N_t/(1+\alpha)$ 是安全的悲观估计，实际表现常常更好。

![原论文 Figure 4.17：Burgers 方程达到 1e-8 所需迭代数](assets/papers/time-parallelization/source-figures/figure-4-17.svg)

非线性收敛分析（Gander and Wu 2020, Section 4）在「精确求解全时间系统 (4.23)」与「$f$ 满足 Lipschitz 条件」两个假设下，给出小 $\alpha$ 时 $\rho=O(\alpha)$，与线性情形一致。Burgers 实验取周期边界、$u_0=\sin^2(2\pi x)$、$\Delta T=0.1$、$J=10$、$\Delta x=1/100$，三条曲线对应 $\nu=1,0.01,10^{-4}$。(a) 固定 $N_t=40$，展示达到全局误差 $10^{-8}$ 所需迭代数随 $\alpha$ 的变化：小 $\alpha$ 加快收敛；对 $\nu$ 的依赖显示小 $\alpha$ 下黏性影响很弱，大 $\alpha$ 下黏性减小才会拖慢收敛。(b) 固定 $\alpha=10^{-3}$，$N_t$ 从 10 到 160，迭代数只在 2–5 轮之间变化，说明收敛率对 $N_t$ 鲁棒。

## 两种路线的最终对照

本节两个变体 (4.15) 与 (4.26) 以不同方式把 ParaDiag 注入标准 Parareal：前者对 $N_t$ 个粗时间点对角化、改动 CGC；后者对每个大区间 $[T_n,T_{n+1}]$ 内的 $J$ 个细时间点对角化、定义特殊粗传播子而保持 CGC 不变。二者适用范围不同——第一种像标准 Parareal 一样主要适合抛物问题，第二种对抛物与双曲问题都有效。

| 问题         | 对角化 CGC (4.15)        | 对角化粗传播 (4.26)                                  |
| ------------ | ------------------------ | ---------------------------------------------------- |
| 对角化方向   | 全局 $N_t$ 个粗点        | 每个粗区间的 $J$ 个细点                              |
| 改动位置     | CGC                      | 粗传播子                                             |
| 粗细积分器   | 可不同                   | 相同积分器与步长                                     |
| 主要适用范围 | 抛物问题                 | 抛物与双曲问题                                       |
| 关键参数     | $\alpha\le\rho/(1+\rho)$ | 抛物因子 $\alpha$，双曲上界 $2\alpha N_t/(1+\alpha)$ |

## 公式、定理与图表覆盖核对

| 原文项目                               | 论文小节 | 覆盖状态                                                           |
| -------------------------------------- | -------- | ------------------------------------------------------------------ |
| (4.14)–(4.17)                          | 4.5.1    | 标准/首尾 CGC、线性全时间矩阵、三步并行解                          |
| Theorem 4.7, Figure 4.12               | 4.5.1    | $\alpha$ 阈值、舍入折中、负实谱已证 vs 复谱数值观察、热与 ADE 实验 |
| (4.18)–(4.19), Figure 4.13             | 4.5.1    | 非线性系统、平均 Jacobian 准 Newton、Burgers 实验                  |
| Remark 4.2, (4.20)                     | 4.5.1    | MGRiT 的一致首尾条件及收敛变体                                     |
| (4.21)–(4.26)                          | 4.5.2    | 同积分器细/粗传播、非线性全时间系统、准 Newton、外层更新           |
| (4.27)–(4.28)                          | 4.5.2    | 线性化、终点提取、$\alpha=0$ 极限与 $J$ 路并行                     |
| Theorem 4.8, (4.29), Figures 4.14–4.17 | 4.5.2    | 负实/纯虚谱界、热、波动与 Burgers 全部原图                         |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 4.5, pp. 460–472.
