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

## 4.1 历史发展

### 为何另设“面向抛物问题的方法”

第二章说明，抛物方程除极低频外具有较强的时间局部性；双曲方程的全部频率都能跨越长时间传播。第三章的方法在时间上同时处理长程耦合，因此也常能处理抛物问题；它们在非线性场景中仍有明显限制：OSWR 的优化 Robin 参数难以确定，ParaExp 与两类 ParaDiag 的外层 Newton 迭代会随时间窗增长而变慢，甚至失效。

第四章的方法主动利用耗散带来的时间局部性。Parareal、PFASST、MGRiT 和 STMG 在扩散充分时对线性、非线性问题都有效；把它们直接搬到弱扩散或双曲问题上，收敛会持续变慢，最后可能发散。

Parareal 的思想可追溯到 Nievergelt（1964）的非迭代先驱工作和多重打靶（Bellen 与 Zennaro 1989；Chartier 与 Philippe 1993）。Saha、Stadel 与 Tremaine（1997）已经用粗模型给出过这一算法，并提到它与波形松弛的联系。Lions、Maday 与 Turinici（2001）在虚拟控制的背景下独立提出现代算法，并强调它是非侵入式的。收敛理论见 Gander 与 Vandewalle（2007）、Gander 与 Hairer（2008, 2014）以及 Gander 与 Lunet（2024）。随后出现 PITA（Farhat 与 Chandesris 2003；Farhat et al. 2006；Cortial 与 Farhat 2009）、PFASST（Minion 2011；Emmett 与 Minion 2012；Minion et al. 2015）、MGRiT（Falgout et al. 2014；Dobrev et al. 2017；Hessenthaler et al. 2020）以及 Parareal 与 ParaDiag 的组合（Wu 2018；Gander 与 Wu 2020）。另一条路线是时空多重网格，从 Hackbusch（1984）的抛物多重网格和 Lubich 与 Ostermann（1987）的多重网格波形松弛开始；早期方法难以有效粗化时间，Gander 与 Neumüller（2016）通过时间块 Jacobi 平滑重新建立了可扩展的 STMG。

在这一谱系中，Parareal 既是两层时间网格方法，也是许多后续 PinT
算法的模板：时间方向使用两层或更多网格，空间方向则仍只保留一层。

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

![原论文 Figure 4.1：每个粗时间步包含 J 个细时间步](assets/papers/time-parallelization/source-figures/figure-4-1.svg)

以下以均匀网格 $\Delta T/\Delta t=J\ge2$ 为主，非均匀网格同样
可用。目标解始终是 $\mathcal F$ 顺序运行得到的离散解；Parareal
只改变求解过程，不改变固定 $\mathcal F$ 后的离散目标。

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
> $M(z)=I_t-M_g^{-1}(z)M_f(z)$。若 $M_fU=b$ 是细传播子的
> 全时间系统，则一次预条件校正显式写成
>
> $$
> M_g(z)\Delta U^k=r^k:=b-M_f(z)U^k,
> \qquad U^{k+1}=U^k+\Delta U^k.
> $$
>
> 第 $n$ 个残差块是
>
> $$
> r_n^k=b_n-
> \left[u_n^k-\mathcal F(T_{n-1},T_n,u_{n-1}^k)\right]
> =b_n-u_n^k+R_f^J(z/J)u_{n-1}^k.
> $$
>
> 因而各残差块中的细传播可以并行计算，$M_g$ 的粗预条件求解仍沿
> 时间顺序进行。这一视角直接导向 Section 4.5 的对角化粗校正。

严格下三角结构还说明：第 $k$ 轮后前 $k$ 个粗节点已与顺序细解一致，精确算术下至多 $N_t$ 轮终止。

### Theorem 4.2：短时间超线性与长时间线性

本定理中的误差记号与 Theorem 4.1 略有不同：这里 $\boldsymbol e_n^k=V_A(\boldsymbol u_n^k-\boldsymbol u_n)$ 是模态误差。$M_g^{-1}$ 的下三角元素为 $R_g^j(z)$。因此

$$
M(z)=[R_f^J(z/J)-R_g(z)]\widetilde M(R_g(z)),
$$

其中 $\widetilde M(\beta)$ 的第一条次对角线为 $1$，再往下依次为 $\beta,\beta^2,\ldots$。Gander 与 Vandewalle（2007, Lemma 4.4）给出

$$
\|\widetilde M^k(R_g)\|_\infty\le
\begin{cases}
\min\left\{
\left(\dfrac{1-|R_g|^{N_t-1}}{1-|R_g|}\right)^{\!k},\
\dbinom{N_t-1}{k}
\right\},
&|R_g|<1,\\[10pt]
\dbinom{N_t-1}{k},
&|R_g|=1,
\end{cases}
$$

下面两类界正来自这两个分支：二项式系数给出 (4.5a) 中的 $\prod_{j=1}^k(N_t-j)/k!$，几何和给出与 $N_t$ 无关的 (4.5b)。

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

Figure 4.2 使用周期热方程、零源、
$u_0(x)=\sin^2(2\pi x)$、$\Delta x=1/5$，粗细层都用后向
Euler，$J=10$。左图 $(T,N_t)=(0.02,6)$，$\varrho_s$ 准确描述
超线性下降；右图 $(T,N_t)=(0.5,64)$，误差按近似固定斜率下降，
$\varrho_l$ 更合适。空间网格细化到 $\Delta x=1/8$ 后也观察到
线性收敛，但原文没有断言它“更早进入”同一两阶段过程。

### Theorem 4.3：非线性超线性界

Theorem 4.3 引自 Gander 与 Hairer（2008, Theorem 1；另见
Gander 与 Lunet 2024, Theorem 2.6）。令 $\mathcal F$ 为精确传播，
$\mathcal G$ 为 $p$ 阶方法，局部截断误差不超过
$C_3\Delta T^{p+1}$。在相关有界状态集上，假设

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
\|(\mathcal F-\mathcal G)(\boldsymbol v)
-(\mathcal F-\mathcal G)(\boldsymbol w)\|
\le C_1\Delta T^{p+1}\|\boldsymbol v-\boldsymbol w\|,
$$

其中 $C_1$ 来自展开系数导数的一致界。对充分小的 $\Delta T$，

$$
\|\boldsymbol u(T_n)-\boldsymbol u_n^k\|
\le
\frac{C_3}{C_1}
\frac{(C_1\Delta T^{p+1})^{k+1}}{(k+1)!}
(1+C_2\Delta T)^{n-k-1}
\prod_{j=0}^{k}(n-j). \tag{4.6}
$$

当 $k\ge n$ 时乘积出现零，前 $n$ 个界面已经精确。小 $\Delta T$ 下，每多一轮又带来一份 $\Delta T^{p+1}$，这就是非线性超线性阶段的来源。

> [!warning] 原文公式核对：(4.6) 多出的一份 $\Delta T^{p+1}$
> 期刊版与 arXiv 版把前因子印成
> $C_3\Delta T^{p+1}(C_1\Delta T^{p+1})^{k+1}$。
> 但被引用的 Gander–Hairer 定理是上面的
> $(C_3/C_1)(C_1\Delta T^{p+1})^{k+1}$。取 $n=1,k=0$ 也能看出：
> 一次粗步的误差应为 $O(\Delta T^{p+1})$，不是
> $O(\Delta T^{2p+2})$。

### Theorem 4.4：约 0.3 的抛物长时间因子

若粗传播 $\mathcal G$ 用后向 Euler，细传播 $\mathcal F$ 用 L-稳定 Runge–Kutta 方法，则存在 $J_{\min}=O(1)$，使

$$
\max_{z\in\mathbb R_-}\varrho_l(J,z)\approx0.3,
\qquad J\ge J_{\min}. \tag{4.7}
$$

定理还假设 $A$ 为负半定矩阵。这个常数不随 $T,N_t$ 增长。
$\mathcal F$ 取后向 Euler 的情形见 Mathew、Sarkis 与
Schaerer（2010），一般 L-稳定 Runge–Kutta 方法见 Yang、Yuan 与
Zhou（2023）。Wu（2015）和 Wu 与 Zhou（2015）还分析了 BDF2
以及下面的 SDIRK 方法；梯形规则不属于本定理的 L-稳定情形，而属于
(4.8) 的有限谱扩展。

这一结果的源头是 Gander 与 Vandewalle（2007, Table 5.1）在**连续层面**（即 $\mathcal F$ 取精确传播子 $\exp(\Delta TA)$）对**其他粗传播子**得到的结论，收缩还可以更好，例如把粗传播子换成 Radau IIA 时约为 $0.068$。注意这里改变的是粗传播子，不是细传播子。

若细传播仅 A-稳定，例如梯形规则，高频不会随
$|z|\to\infty$ 消失。由于 $z=\Delta T\lambda(A)\le0$，有限谱区间应
写成 $[-z_{\max},0]$：

$$
\max_{z\in[-z_{\max},0]}\varrho_l(J,z)\approx0.3,
\qquad
J\ge J_{\min}=O(\log^2 z_{\max}). \tag{4.8}
$$

> [!warning] 原文公式核对：(4.8) 的谱区间
> 正式版印成 $z\in[0,z_{\max}]$，但同一节定义
> $z=\Delta T\lambda(A)$ 且 $A$ 负半定，所以 $z\le0$。也可令
> $s=-z\ge0$，把左端写成
> $\max_{s\in[0,z_{\max}]}\varrho_l(J,-s)$。

细层需要更多小步，才能让其高频行为与物理耗散相符。这与 $\mathcal F$ 取精确传播子 $\exp(\Delta TA)$ 的情形差别很大：后者只要 $J\ge2$ 就有约 $0.3$ 的收敛速率。(4.8) 由 Wu 与 Zhou（2015）对梯形规则和一个四阶 Gauss Runge–Kutta 方法证明。

两个 SDIRK 细传播子把稳定性条件具体化：

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

SDIRK22 的 (4.7) 在 $J_{\min}=2$ 成立，SDIRK23 需要 $J_{\min}=4$（Wu 2015；Wu 与 Zhou 2015）。

> [!warning] 原文公式核对：SDIRK22 的 Butcher 表
> 期刊版与 arXiv 版在 (4.9) 中把 SDIRK22 的第二个节点印成 $c_2=1-\gamma$，权重印成 $b=(1-\gamma,1-\gamma)$。取 $\gamma=(2-\sqrt2)/2$ 时 $b$ 之和为 $\sqrt2\ne1$，且 $c_2$ 与该行系数之和 $a_{21}+a_{22}=1$ 不符。上表按标准的 L-稳定、刚性精确 SDIRK22（$c_2=1$，$b$ 取 $A$ 的最后一行）排印。相邻的 SDIRK23 与原文一致，无需修正。

![原论文 Figure 4.3：不同细传播子与粗细比带来的 Parareal 收敛差异](assets/papers/time-parallelization/source-figures/figure-4-3.svg)

Figure 4.3 取周期热方程、$\Delta x=1/256$、$\Delta T=0.1$、$T=4$、扩散系数 $0.1$。三个面板从左到右对应 $J=2,10,50$。$J=2$ 时，梯形规则在约 $10^{-4}$ 处停滞，SDIRK23 也明显慢于 SDIRK22；$J=10$ 时两种 SDIRK 曲线已经接近，而梯形规则仍保留较慢的尾段；$J=50$ 时三种细传播都贴近 $0.3^k$ 参考斜率。图中变化验证了 Theorem 4.4 的限定条件：细传播子的稳定性类型和粗细比需要一起判断。

### 扩散减弱后的退化

负实轴上的结论不会自动延伸到双曲极限。为检验边界，下面固定
$T=4$、$\Delta T=0.1$、$\Delta x=1/128$、$J=32$，粗层用
后向 Euler，细层用 SDIRK22，并逐步降低周期 ADE 与 Burgers 的黏性。

![原论文 Figure 4.4：三种黏性下每个对流扩散特征值对应的长时间因子](assets/papers/time-parallelization/source-figures/figure-4-4.svg)

实验取零源项和 $u(x,0)=\sin(2\pi x)$。$\nu$ 减小时，$\max\varrho_l$ 逼近 $1$：三个面板分别标出 $\nu=1$ 时 $\varrho_{l,\max}=0.23$、$\nu=0.1$ 时 $0.39$、$\nu=0.02$ 时 $0.79$，说明粗传播越来越难修正长期存在的传播模态。

![原论文 Figure 4.5：对流扩散与 Burgers 方程上黏性降低导致的 Parareal 退化](assets/papers/time-parallelization/source-figures/figure-4-5.svg)

Figure 4.5 使用与 Figure 4.4 相同的三个黏性值，其中 (a) 与谱因子预测一致。Burgers 方程缺少同样精确的模态分析，但 Figure 4.5(b) 呈现相同趋势；约在 $\nu\le10^{-3}$ 时常规迭代会发散。严格下三角结构保证的有限步性质仍在，只是需要的轮数失去实用价值。

双曲困难的机制在于：任意小的高频分量都能传播任意远，很难让
$\mathcal G$ 与 $\mathcal F$ 同时便宜又足够接近；一旦把粗传播做得
足够准确，粗校正本身又会吞掉加速收益。相关分析见 Gander 与
Vandewalle（2007）、Gander 与 Lunet（2020a,b）、Gander、Lunet 与
Pogoželskytė（2023a）以及 Gander、Lunet、Ruprecht 与 Speck
（2023b）。

MGRiT 把 Parareal 推广到多层。为使它适用于输运问题，一条路线使用
半 Lagrange 优化粗求解器（Howse et al. 2019；De Sterck et al.
2021, 2023a, 2023b），但非线性情形仍然困难；另一条路线对角化粗
求解器（Gander 与 Wu 2020），见第 4.5 节。

## 公式、定理与图表覆盖核对

| 原文项目                             | 论文小节 | 覆盖状态                                        |
| ------------------------------------ | -------- | ----------------------------------------------- |
| Section 4 引言与 4.1                 | 4、4.1   | 抛物局部性、四类算法、历史路线                  |
| (4.1), Figure 4.1                    | 4.2      | 更新式、并行细传播、顺序粗校正、双层时间网格    |
| (4.2)–(4.4), Theorem 4.1, Remark 4.1 | 4.2      | 模态化、Toeplitz 误差矩阵、完整证明、预条件解释 |
| (4.5), Theorem 4.2, Figure 4.2       | 4.2      | 超线性与长时间线性界、两种实验阶段              |
| (4.6), Theorem 4.3                   | 4.2      | 非线性假设、误差界和有限步含义                  |
| (4.7)–(4.9), Theorem 4.4, Figure 4.3 | 4.2      | L/A 稳定细传播、0.3 因子、SDIRK 对比            |
| Figures 4.4–4.5                      | 4.2      | 对流扩散谱、Burgers 实验和弱扩散失效            |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 4–4.2, pp. 443–452.
