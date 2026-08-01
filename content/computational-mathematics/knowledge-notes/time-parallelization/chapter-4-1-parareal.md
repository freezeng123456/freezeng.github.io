---
title: 4.1–4.2：历史脉络与 Parareal
description: 从多重打靶、粗细传播子到线性与非线性收敛界的逐式精读
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-1-parareal
tags:
  - 时间并行
  - Parareal
  - 抛物方程
---

> [!note] 阅读范围
> 本页对应论文 Sections 4、4.1 和 4.2（pp. 443–452），覆盖公式 (4.1)–(4.9)、Theorems 4.1–4.4、Remark 4.1 和 Figures 4.1–4.5。公式后的推理按论文证明顺序展开。

## 4.1.1 为何另设“面向抛物问题的方法”

第二章说明，抛物方程除极低频外具有较强的时间局部性；双曲方程的全部频率都能跨越长时间传播。第三章的方法在时间上同时处理长程耦合，因此也常能处理抛物问题；它们在非线性场景中仍有明显限制：OSWR 的优化 Robin 参数难以确定，ParaExp 与两类 ParaDiag 的外层 Newton 迭代会随时间窗增长而变慢，甚至失效。

第四章的方法主动利用耗散带来的时间局部性。Parareal、PFASST、MGRIT 和 STMG 在扩散充分时对线性、非线性问题都有效；把它们直接搬到弱扩散或双曲问题上，收敛会持续变慢，最后可能发散。

Parareal 的思想可追溯到多重打靶、波形松弛和 Nievergelt (1964) 的非迭代先驱工作。Lions et al. (2001) 独立提出现代算法。随后出现 PITA、PFASST、MGRIT、Parareal–ParaDiag 组合。另一条路线是时空多重网格：早期方法难以有效粗化时间；Gander and Neumüller (2016) 通过时间块 Jacobi 平滑重新建立了可扩展的 STMG。

## 4.2.1 更新公式与两层时间网格

将 $[0,T]$ 分成 $0=T_0<T_1<\cdots<T_{N_t}=T$。$\mathcal F$ 是昂贵而准确的细传播子，$\mathcal G$ 是便宜的粗传播子。Parareal 从界面初值猜测 $\boldsymbol u_n^0$ 出发：

$$
\boldsymbol u_{n+1}^{k+1}
=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^k). \tag{4.1}
$$

同一轮中，所有 $\mathcal F(T_n,T_{n+1},\boldsymbol u_n^k)$ 可并行；带新迭代指标 $k+1$ 的粗传播必须沿 $n$ 顺序推进。最后两项形成粗网格上的预测–校正，也可由多重打靶 Newton 的有限差分 Jacobian 解释。

![原论文 Figure 4.1：每个粗时间步包含 J 个细时间步](assets/papers/time-parallelization/source-figures/figure-4-1.svg)

论文主要讨论均匀网格，令 $\Delta T/\Delta t=J\ge2$；非均匀网格也可使用。目标解是 $\mathcal F$ 顺序运行得到的离散解，Parareal 并不在固定 $\mathcal F$ 之外额外改变离散目标。

## 4.2.2 Theorem 4.1：把误差化成逐模态 Toeplitz 迭代

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

### 证明链条

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

把所有 $n$ 叠成 $\boldsymbol\xi^k$，有

$$
M_g(z)\boldsymbol\xi^{k+1}
=[M_g(z)-M_f(z)]\boldsymbol\xi^k,
$$

于是 $\boldsymbol\xi^k=M^k(z)\boldsymbol\xi^0$，再取模态和时间节点的最大值便得 (4.2)。

> [!note] Remark 4.1：预条件视角
> $M(z)=I_t-M_g^{-1}(z)M_f(z)$。$M_fU=b$ 是细传播子的全时间系统，$M_g$ 是粗传播预条件器。残差中的每个细传播可以并行计算，粗预条件求解则沿时间顺序进行。这一视角直接导向 Section 4.5 的对角化粗校正。

严格下三角结构还说明：第 $k$ 轮后前 $k$ 个粗节点已与顺序细解一致，精确算术下至多 $N_t$ 轮终止。

## 4.2.3 Theorem 4.2：短时间超线性与长时间线性

$M_g^{-1}$ 的下三角元素为 $R_g^j(z)$。因此

$$
M(z)=[R_f^J(z/J)-R_g(z)]\widetilde M(R_g(z)),
$$

其中 $\widetilde M(\beta)$ 的第一条次对角线为 $1$，再往下依次为 $\beta,\beta^2,\ldots$。对其幂取无穷范数得到两类界。

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

$k=N_t$ 时乘积含零，显式体现有限步收敛。若 $|R_g(z)|<1$，还可用与 $N_t$ 无关的长时间界

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

![原论文 Figure 4.2：短时间超线性与长时间线性两种收敛阶段](assets/papers/time-parallelization/source-figures/figure-4-2.svg)

Figure 4.2 使用周期热方程、零源、$u_0(x)=\sin^2(2\pi x)$、$\Delta x=1/5$，粗细层都用后向 Euler，$J=10$。$T=0.02,N_t=6$ 时 $\varrho_s$ 准确描述超线性下降；更长区间中，误差按近似固定斜率下降，$\varrho_l$ 更合适。把空间网格细化到 $\Delta x=1/8$ 后，也会更早进入线性阶段。

## 4.2.4 Theorem 4.3：非线性超线性界

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

## 4.2.5 Theorem 4.4：约 0.3 的抛物长时间因子

若粗传播 $\mathcal G$ 用后向 Euler，细传播 $\mathcal F$ 用 L-稳定 Runge–Kutta 方法，则存在 $J_{\min}=O(1)$，使

$$
\max_{z\in\mathbb R_-}\varrho_l(J,z)\approx0.3,
\qquad J\ge J_{\min}. \tag{4.7}
$$

这个常数不随 $T,N_t$ 增长。它来自负实谱、粗层的高频耗散和有限的粗细传播差。用 Radau IIA 等组合时最坏因子还可降到约 $0.068$。

若细传播仅 A-稳定，例如梯形规则，高频不会随 $|z|\to\infty$ 消失。此时在有限谱区间 $[0,z_{\max}]$ 上仍有

$$
\max_{z\in[0,z_{\max}]}\varrho_l(J,z)\approx0.3,
\qquad
J\ge J_{\min}=O(\log_2 z_{\max}). \tag{4.8}
$$

细层需要更多小步，才能让其高频行为与物理耗散相符。

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

SDIRK22 的 (4.7) 在 $J_{\min}=2$ 成立，SDIRK23 需要 $J_{\min}=4$。

![原论文 Figure 4.3：不同细传播子与粗细比带来的 Parareal 收敛差异](assets/papers/time-parallelization/source-figures/figure-4-3.svg)

Figure 4.3 取周期热方程、$\Delta x=1/256$、$\Delta T=0.1$、$T=4$、扩散系数 $0.1$。$J$ 较小时，梯形规则与 SDIRK23 较慢；$J=50$ 时三种细方法都接近 0.3 的斜率。

## 4.2.6 扩散减弱后的退化

论文随后固定 $T=4$、$\Delta T=0.1$、$\Delta x=1/128$、$J=32$，粗层用后向 Euler，细层用 SDIRK22，考察周期对流扩散和 Burgers 方程。

![原论文 Figure 4.4：三种黏性下每个对流扩散特征值对应的长时间因子](assets/papers/time-parallelization/source-figures/figure-4-4.svg)

对流扩散的谱由负实轴附近逐渐向虚轴展开。$\nu$ 减小时，$\max\varrho_l$ 逼近 $1$，说明粗传播越来越难修正长期存在的传播模态。

![原论文 Figure 4.5：对流扩散与 Burgers 方程上黏性降低导致的 Parareal 退化](assets/papers/time-parallelization/source-figures/figure-4-5.svg)

Figure 4.5(a) 与谱因子预测一致。Burgers 方程缺少同样精确的模态分析，但 Figure 4.5(b) 呈现相同趋势；约在 $\nu\le10^{-3}$ 时常规迭代会发散。严格下三角结构保证的有限步性质仍在，只是需要的轮数失去实用价值。波动方程也通常不收敛，这正是第四章方法适用范围的边界。

## 公式、定理与图表覆盖核对

| 原文项目                             | 本页位置 | 覆盖状态                                        |
| ------------------------------------ | -------- | ----------------------------------------------- |
| Section 4 引言与 4.1                 | §4.1.1   | 抛物局部性、四类算法、历史路线                  |
| (4.1), Figure 4.1                    | §4.2.1   | 更新式、并行细传播、顺序粗校正、双层时间网格    |
| (4.2)–(4.4), Theorem 4.1, Remark 4.1 | §4.2.2   | 模态化、Toeplitz 误差矩阵、完整证明、预条件解释 |
| (4.5), Theorem 4.2, Figure 4.2       | §4.2.3   | 超线性与长时间线性界、两种实验阶段              |
| (4.6), Theorem 4.3                   | §4.2.4   | 非线性假设、误差界和有限步含义                  |
| (4.7)–(4.9), Theorem 4.4, Figure 4.3 | §4.2.5   | L/A 稳定细传播、0.3 因子、SDIRK 对比            |
| Figures 4.4–4.5                      | §4.2.6   | 对流扩散谱、Burgers 实验和弱扩散失效            |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 4–4.2, pp. 443–452.
