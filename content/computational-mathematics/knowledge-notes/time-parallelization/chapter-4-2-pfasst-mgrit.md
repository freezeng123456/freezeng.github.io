---
title: 4.3–4.4：PFASST 与 MGRiT
description: 从配置方程、SDC 平滑与层间插值到 FCF 重叠 Parareal，完整推导两类算法的矩阵结构、收敛因子与等成本比较
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-4-2-pfasst-mgrit
tags:
  - 时间并行
  - PFASST
  - MGRiT
---

> [!note] 阅读范围
> 本页对应论文 Sections 4.3–4.4（pp. 452–460），覆盖公式 (4.10)–(4.13)、Theorems 4.5–4.6 和 Figures 4.6–4.11。PFASST 的配置方程、层间传递、SDC 近似以及 MGRiT 的重叠解释均展开到可实现的矩阵形式，并逐步补足“为什么这样构造”的推理链。

## 4.3 PFASST

### 来源与两层节点

PFASST（parallel full approximation scheme in space–time，时空并行全近似格式）由 Emmett and Minion (2012) 正式提出。它的核心想法其实早两年就已出现：Minion (2010) 为了压低 Parareal 中一次迭代的成本，用**一次 SDC（spectral deferred correction，谱延迟校正）迭代替代昂贵的完整细传播** $\mathcal F$。SDC 本身源自 Dutt, Greengard and Rokhlin (2000)，是一种把高阶配置解通过低阶扫掠（sweep）逐次逼近的框架。因此从一开始 PFASST 就带有“用廉价近似传播换取并行度”的基因。

PFASST 长期缺乏清晰的描述与理论分析，原因在于它同时叠加了三重结构：时间分块的 Parareal 迭代、每个块内的 SDC 扫掠、以及细/粗两层配置。直到 Bolten, Moser and Speck (2017) 借助 Minion et al. (2015) 给出的 **SDC 代数表示**，把 PFASST 重新解释为一种**时间多重网格**方法，并在 Bolten, Moser and Speck (2018) 中给出了收敛分析，这一算法才有了可分析的代数骨架。Gander et al. (2023b) 进一步用块迭代（block iteration）形式对模型问题给出精确刻画，本页采用的正是这一表述。

> [!tip] 本站洞见
> 把 PFASST 读成“多层 SDC”是理解它的钥匙。SDC 的一次扫掠本质上是对稠密配置系统 $\Phi$ 的一个**低阶预处理器求解** $\widetilde\Phi^{-1}$（见 (4.11)），扮演多重网格中的“光滑器”；而细/粗两层配置加上 Lagrange 传递算子扮演“粗网格校正”。于是 PFASST = 时间分块的 Parareal 外层 + 每块内的两层 SDC V-循环。三重循环之所以能并行，是因为最贵的细层扫掠对各时间块是解耦的，只有廉价的粗层沿时间串行传播。

把 $(0,T)$ 分成 $N_t$ 个大区间 $[T_n,T_{n+1}]$。每个区间内设置 $M_f$ 个细节点和 $M_c$ 个粗节点：

$$
t_{n,m}^{f}=T_n+\tau_m^f\Delta t,
\qquad
t_{n,m}^{c}=T_n+\tau_m^c\Delta t,
$$

其中 $\tau_0^{f,c}=0$、$\tau_{M_{f,c}}^{f,c}=1$，且 $M_f>M_c$。要求 $M_f>M_c$ 是为了让两层承担不同角色：细层节点多、配置阶数高，负责精度；粗层节点少、系统小，负责廉价地把校正沿时间传播。上标 $f$、$c$ 分别标记细、粗时间网格。

### 配置方程 (4.10)

对 $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$，在任一层的 $M$ 个节点上用数值求积（配置）写成

$$
\boldsymbol u_{n,m}
=\boldsymbol u_{n,0}
+\Delta t\sum_{j=1}^{M}q_{m,j}
[A\boldsymbol u_{n,j}+\boldsymbol g(t_{n,j})],
\qquad m=1,\ldots,M. \tag{4.10}
$$

其中 $\boldsymbol u_{n,j}$ 是 $\boldsymbol u$ 在 $t=t_{n,j}$ 的近似，$M=M_f$ 或 $M_c$、$t_{n,j}=t_{n,j}^f$ 或 $t_{n,j}^c$。(4.10) 是把区间上的积分 $\int_{T_n}^{t_{n,m}}$ 用配置权重 $q_{m,j}$ 逼近的结果——每个节点值都依赖**全部**节点，因此这是一个把 $M$ 个节点耦合在一起的稠密隐式系统，直接求解成本很高。这正是后面要用 SDC 扫掠来近似求解的原因。

令

$$
Q=(q_{m,j}),\qquad
\boldsymbol u_n=(\boldsymbol u_{n,1}^\top,\ldots,\boldsymbol u_{n,M}^\top)^\top,
$$

$$
\chi=
\begin{bmatrix}
0&\cdots&0&1\\
\vdots&&\vdots&\vdots\\
0&\cdots&0&1
\end{bmatrix},
\qquad \boldsymbol\chi=\chi\otimes I_x,
$$

则末节点被复制为下一区间所有配置节点的初值：

$$
\boldsymbol u_n
=\Delta t(Q\otimes A)\boldsymbol u_n
+\boldsymbol\chi\boldsymbol u_{n-1}
+\Delta t\boldsymbol b_n.
$$

这里 $\boldsymbol b_n=(Q\otimes I_x)\boldsymbol g_n$，而复制矩阵 $\chi$ 的每一行都是 $(0,\ldots,0,1)$：它把上一区间的**末节点** $\boldsymbol u_{n-1,M}$（即 $t=T_n$ 处的解，因为 $\tau_M=1$）抽出来，作为本区间所有配置节点的公共初值 $\boldsymbol u_{n,0}$。这一“复制”动作是块与块之间唯一的串行耦合，也正是 Parareal 式序列依赖在配置层面的体现。

细、粗层分别记为

$$
\Phi_f=I_f-\Delta tQ_f\otimes A,
\qquad
\Phi_c=I_c-\Delta tQ_c\otimes A,
$$

$$
\boldsymbol u_n^f
=\Phi_f^{-1}(\boldsymbol\chi_f\boldsymbol u_{n-1}^f+\Delta t\boldsymbol b_n^f),
\qquad
\boldsymbol u_n^c
=\Phi_c^{-1}(\boldsymbol\chi_c\boldsymbol u_{n-1}^c+\Delta t\boldsymbol b_n^c).
$$

其中 $I_c=I_{M_c}\otimes I_x$、$I_f=I_{M_f}\otimes I_x$。矩阵 $\Phi_f,\Phi_c$ 就是配置系统的系统矩阵：给定上一区间末值，求解 $\Phi^{-1}$ 即求解该区间的整块配置解。$\Phi_f$ 稠密且大（$M_f$ 个节点），是全流程中最贵的一步，也正是需要被 SDC 扫掠近似替代的对象。

### Lagrange 层间传递

要在两层配置之间搬运数据，需要延拓 $T_{c\to f}$（粗→细）与限制 $T_{f\to c}$（细→粗）算子。它们都由 Lagrange 插值定义：粗层节点上的数据张成插值多项式

$$
p^c(\tau;\boldsymbol u^c)
=\sum_{m=1}^{M_c}u_m^cL_m^c(\tau),
\qquad
L_m^c(\tau)=
\prod_{\substack{j=1\\j\ne m}}^{M_c}
\frac{\tau-\tau_j^c}{\tau_m^c-\tau_j^c}.
$$

在细节点 $\{\tau_m^f\}$ 处取值即得延拓矩阵 $T_{c\to f}$；反过来用细层 Lagrange 基 $L_m^f$ 在粗节点取值，得到限制矩阵 $T_{f\to c}$。两个矩阵都与 $I_x$ 做 Kronecker 积，因此它们只沿配置节点方向做多项式插值，**完全不触碰空间变量**——这与空间多重网格里限制/延拓算子会重采样空间自由度形成对照，是“时间方向多重网格”的特征。

> [!tip] 本站洞见
> 传递算子选 Lagrange 插值而非简单注入或平均，是因为配置节点通常取 Radau/Gauss 这类非等距的高斯型点。用节点定义的插值多项式在两层间转换，可保持配置解所隐含的多项式表示一致，避免在层间传递时引入额外的低阶截断误差污染高阶细解。

在 Gander et al. (2023b) 的块迭代表述中，PFASST 写成

$$
\boldsymbol u_{n+1}^{k+1}
=B_{01}\boldsymbol u_{n+1}^k
+B_{10}(\boldsymbol\chi\boldsymbol u_n^{k+1}+\Delta t\boldsymbol b_n)
+B_{00}(\boldsymbol\chi\boldsymbol u_n^k+\Delta t\boldsymbol b_n),
$$

$$
\begin{aligned}
B_{01}&=[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]
(I_f-\widetilde\Phi_f^{-1}\Phi_f),\\
B_{10}&=T_{c\to f}\Phi_c^{-1}T_{f\to c},\\
B_{00}&=[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]\widetilde\Phi_f^{-1}.
\end{aligned}
$$

三个下标遵循块迭代约定：第一个下标标记作用于本区间当前迭代（下标 1）还是相邻区间（下标 0），第二个下标区分新迭代 $k{+}1$（下标 0）与旧迭代 $k$（下标 1）的输入。逐块理解其结构，就能看清 PFASST 的多重网格骨架：

- 公共因子 $[I_f-T_{c\to f}\Phi_c^{-1}T_{f\to c}\Phi_f]$ 是**粗网格校正算子**。它把细层残差经 $T_{f\to c}$ 限制到粗层，用 $\Phi_c^{-1}$ 廉价求解粗配置系统，再经 $T_{c\to f}$ 延拓回细层并从恒等中扣除——这正是两层多重网格里“限制—粗解—延拓”的标准结构。
- 因子 $(I_f-\widetilde\Phi_f^{-1}\Phi_f)$ 是**细层 SDC 光滑器的迭代矩阵**：$\widetilde\Phi_f^{-1}$ 是一次廉价扫掠（预处理器），$(I_f-\widetilde\Phi_f^{-1}\Phi_f)$ 就是它作用一次后的误差传播算子。于是 $B_{01}$ = 先在旧迭代上做一次细层 SDC 光滑、再做粗网格校正，恰是一个两层 V-循环。
- $B_{10}=T_{c\to f}\Phi_c^{-1}T_{f\to c}$ 用**新迭代**的相邻区间数据（含复制初值 $\boldsymbol\chi\boldsymbol u_n^{k+1}$）在粗层向前传播一步——这条串行链正是让整体收敛在有限步内的粗层“信息高速公路”。
- $B_{00}$ 用**旧迭代**的相邻数据做校正，$\widetilde\Phi_f^{-1}$ 再次体现廉价细扫掠；它与 $B_{10}$ 的差别相当于 Parareal 校正中新旧粗传播之差的角色。

其中 $\widetilde\Phi_f$ 是易解的细层配置矩阵近似（预处理器）。

> [!tip] 本站洞见
> 把这三块与经典 Parareal 校正 $\mathcal G_{new}-\mathcal G_{old}+\mathcal F_{old}$ 对照，会发现 PFASST 只是把“完整细传播 $\mathcal F$”替换成了“细层 SDC 光滑 $\widetilde\Phi_f^{-1}$ + 粗网格校正”。这解释了它为什么比 Parareal 便宜：单次迭代不再要求解稠密的 $\Phi_f^{-1}$，只需一次三角化的扫掠 $\widetilde\Phi_f^{-1}$ 加一个小粗系统 $\Phi_c^{-1}$。

### SDC 近似 (4.11) 与 Figure 4.6

实际中用细节点间的隐式 Euler 构造预处理器 $\widetilde\Phi_f$：

$$
\frac{\boldsymbol u_{n,m+1}-\boldsymbol u_{n,m}}
{\Delta t(\tau_{m+1}^f-\tau_m^f)}
=A\boldsymbol u_{n,m+1}+\boldsymbol g(t_{n,m+1}^f),
\quad m=0,\ldots,M_f-1. \tag{4.11}
$$

即

$$
\widetilde\Phi_f=
\begin{bmatrix}
1&&&\\
-1&1&&\\
&\ddots&\ddots&\\
&&-1&1
\end{bmatrix}\otimes I_x
-\Delta t
\begin{bmatrix}
\tau_1^f-\tau_0^f&&&\\
&\tau_2^f-\tau_1^f&&\\
&&\ddots&\\
&&&\tau_{M_f}^f-\tau_{M_f-1}^f
\end{bmatrix}\otimes A.
$$

它由“节点差构成的对角块”乘 $I_x$、再减去“下双对角差分矩阵”与 $\Delta tA$ 的 Kronecker 项组成。关键在于 $\widetilde\Phi_f$ 是**下双对角（块下三角）**的：相邻节点用隐式 Euler 逐段推进，因此 $\widetilde\Phi_f^{-1}$ 可以一次前向扫掠、逐节点求解，每个节点只需解一个 $(I-\Delta t(\tau_{m+1}^f-\tau_m^f)A)$ 型的空间子系统。这与稠密的 $\Phi_f$ 形成鲜明对比：一次 $\widetilde\Phi_f^{-1}$ 就对应一轮低阶 SDC 扫掠，成本远低于完整配置求解，却在多重网格意义上足以充当高频光滑器。

> [!tip] 本站洞见
> 为什么隐式 Euler 差分能当好光滑器？因为 $\widetilde\Phi_f$ 与 $\Phi_f$ 共享同一套配置节点，只是把稠密求积权重 $Q_f$ 替换成节点差构成的下三角“近似积分”。二者之差集中在描述节点间高阶耦合的部分，恰是空间高频、时间高阶的误差分量；隐式 Euler 对这些分量强阻尼，因而 $(I_f-\widetilde\Phi_f^{-1}\Phi_f)$ 能高效衰减它们，把残余的低频误差留给粗网格校正处理。这正是“光滑器 + 粗校正”分工的由来。

数值实验取 $T=3$、周期边界、零初值，源项用 (2.4) 且 $\sigma=1000$，$\Delta x=1/128$、$\Delta t=1/64$。细层采用三节点 Radau IIA（$M_f=3$）：

$$
\left\{0,\frac{4-\sqrt6}{10},\frac{4+\sqrt6}{10},1\right\},
$$

粗层采用两节点 Radau IIA（$M_c=2$）$\{0,1/3,1\}$。选 Radau IIA 是因为它右端含节点 $\tau_M=1$（配合复制矩阵 $\chi$ 抽取区间末值）且 L-稳定，适合刚性问题。相应的求积矩阵

$$
Q_f=
\begin{bmatrix}
\frac{88-7\sqrt6}{360}&\frac{296-169\sqrt6}{1800}&\frac{-2+3\sqrt6}{225}\\[2pt]
\frac{296+169\sqrt6}{1800}&\frac{88+7\sqrt6}{360}&\frac{-2-3\sqrt6}{225}\\[2pt]
\frac{16-\sqrt6}{36}&\frac{16+\sqrt6}{36}&\frac19
\end{bmatrix},
\qquad
Q_c=
\begin{bmatrix}
\frac{5}{12}&-\frac{1}{12}\\[2pt]
\frac34&\frac14
\end{bmatrix},
$$

由此得到的两个传递矩阵为

$$
T_{c\to f}=
\begin{bmatrix}
1.2674&-0.2674\\
0.5325&0.4674\\
0&1
\end{bmatrix}\otimes I_x,
\qquad
T_{f\to c}=
\begin{bmatrix}
0.5018&0.6833&-0.1851\\
0&0&1
\end{bmatrix}\otimes I_x.
$$

这些矩阵完全由节点位置和 Lagrange 插值确定，不含任何可调参数。注意两个传递矩阵的末行都是 $(0,\ldots,0,1)$：末节点 $\tau=1$ 在两层都存在，插值退化为恒等，保证区间末值在层间传递时不失真。

![原论文 Figure 4.6：热方程与三种黏性对流扩散方程上的 PFASST 误差](assets/papers/time-parallelization/source-figures/figure-4-6.svg)

热方程收敛最快。对流扩散的黏性减小时，持续传播的高频越来越难由粗配置层表示，收敛按同一方向恶化。这与前面 Parareal、MGRiT 中观察到的粗细传播差机制一致：粗层配置（节点少）无法准确表示弱扩散下长寿命、几乎无衰减的高频模态，粗网格校正随之失效。

## 4.4 MGRiT

### 三种解释：作为带重叠的 Parareal

时间多重网格缩减（multigrid reduction in time，MGRiT）由 Falgout et al. (2014) 提出，是 Parareal 的又一变体。它至少有三种等价视角：

1. **代数多重网格 + FCF 松弛**：把时间层级视作 AMG 的粗/细自由度，F 松弛更新细点、C 松弛更新粗点，FCF 即 F–C–F 一轮复合松弛；
2. **块迭代**（Gander et al. 2023b；Gander and Lunet 2024, 第 4.6 章）；
3. **带重叠的 Parareal**（Gander et al. 2018b, Theorem 4 与 Corollary 1）。

本页采用第三种视角，直接对非线性 ODE 系统 (2.2) $\boldsymbol u'=\boldsymbol f(\boldsymbol u,t)$、$\boldsymbol u(0)=\boldsymbol u_0$ 给出两层 FCF 松弛的更新：

$$
\boldsymbol u_0^{k+1}=\boldsymbol u_0,
\qquad
\boldsymbol u_1^{k+1}=\mathcal F(T_0,T_1,\boldsymbol u_0),
$$

$$
\begin{aligned}
\boldsymbol u_{n+1}^{k+1}
={}&\mathcal F\!\left(T_n,T_{n+1},
\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)\right)\\
&+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})\\
&-\mathcal G\!\left(T_n,T_{n+1},
\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)\right),
\quad n=1,\ldots,N_t-1. \tag{4.12}
\end{aligned}
$$

其中 $\mathcal G,\mathcal F$ 是 Parareal 中的粗、细传播（见 (4.1)）。对照 Parareal 校正 $\mathcal F(\boldsymbol u_{n-1}^k)+\mathcal G(\boldsymbol u_n^{k+1})-\mathcal G(\boldsymbol u_n^k)$ 可以看出唯一区别：MGRiT 把被校正的“上一区间值”从 $\boldsymbol u_{n-1}^k$ 换成了**先做一次细传播** $\mathcal F(T_{n-1},T_n,\boldsymbol u_{n-1}^k)$。因此 (4.12) 每轮出现两个细传播（内层 $\mathcal F$ 加外层 $\mathcal F$），而 Parareal 每轮只有一个。

![原论文 Figure 4.7：FCF-MGRiT 是重叠宽度为一个粗时间步的 Parareal](assets/papers/time-parallelization/source-figures/figure-4-7.svg)

Figure 4.7 中深色圆是运行粗传播 $\mathcal G$ 的粗点。额外那次内层 F 松弛把旧迭代信息先向前推进一个粗区间 $\Delta T$，再进入常规的 Parareal 式校正——几何上等价于让相邻两个 Parareal“窗口”重叠了一个粗步。因此两层 FCF-MGRiT 就是**重叠宽度为 $\Delta T$ 的 Parareal**，同样在有限步内精确收敛：全局误差至多经 $k=\lceil N_t/2\rceil$ 轮归零（Gander et al. 2018b, Theorem 5）。更一般地，$F(CF)^\nu$ 松弛对应重叠宽度 $\nu\Delta T$ 的 Parareal（Gander et al. 2018b, Theorem 6 与 Corollary 1 给出了相应的超线性收敛结果）。

> [!tip] 本站洞见
> 为什么重叠能把精确收敛步数从 Parareal 的 $N_t$ 降到 $\lceil N_t/2\rceil$？Parareal 每轮把“已被精确求解”的时间前沿只推进一个区间；而 FCF 的额外 F 松弛让前沿每轮多跨一个粗区间，前进速度翻倍，于是达到全局精确所需轮数减半。代价是每轮多一次细求解——这恰好把“迭代次数减半”换成了“单轮成本翻倍”，两者在细求解总数上打平。这是后面等成本比较的直觉根源。

### Theorem 4.5：长时间模态因子

线性情形的收敛估计来自 Dobrev et al. (2017)。考虑 $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$，$A$ 可对角化且谱 $\sigma(A)\subset\mathbb C^-$。沿用 Theorem 4.2 的记号，若粗传播稳定 $|R_g(z)|<1$，两层 FCF-MGRiT 满足

$$
\max_n\|\boldsymbol e_n^k\|_\infty
\le
\max_{z\in\sigma(\Delta TA)}\varrho_l^k(J,z)
\max_n\|\boldsymbol e_n^0\|_\infty,
$$

$$
\varrho_l(J,z)=
\frac{|R_f^J(z/J)|\,|R_g(z)-R_f^J(z/J)|}
{1-|R_g(z)|}. \tag{4.13}
$$

其中 $R_g,R_f$ 是粗传播 $\mathcal G$、细传播 $\mathcal F$ 的稳定性函数，$\boldsymbol e_n^k:=V_A(\boldsymbol u_n^k-\boldsymbol u_n)$，$J=\Delta T/\Delta t$ 是每个粗区间内的细步数。$R_f^J(z/J)$ 表示以步长 $\Delta t=\Delta T/J$（对应参数 $z/J$）连续走 $J$ 个细步、即在一个粗区间上的整段细传播。与 Parareal 收敛因子 (4.5b)

$$
\varrho_{l,\mathrm{Parareal}}
=\frac{|R_g(z)-R_f^J(z/J)|}{1-|R_g(z)|}
$$

对照，立即得到

$$
\varrho_{l,\mathrm{MGRiT}}
=|R_f^J(z/J)|\,\varrho_{l,\mathrm{Parareal}}.
$$

> [!tip] 本站洞见
> 为什么 FCF 恰好在 Parareal 因子前乘上一个 $|R_f^J(z/J)|$？看 (4.12)：被送进“粗校正差 $\mathcal G_{new}-\mathcal G_{old}$”的不再是旧误差本身，而是先经一段完整细传播 $\mathcal F$ 之后的误差。线性化后，这段细传播把误差乘以 $R_f^J(z/J)$，随后再乘以 Parareal 那套校正因子 $|R_g-R_f^J|/(1-|R_g|)$。于是额外 F 松弛的效果精确地表现为“多乘一个细传播幅度”。因为细传播稳定，$|R_f^J(z/J)|\le 1$，所以这一步只会收缩、不会放大——但它也实打实地多花一次昂贵的并行细求解。换言之，FCF 用“更多细求解”买“更小收敛因子”，是否划算取决于按什么单位比较。

另外，(4.13) 的分母 $1-|R_g(z)|$ 揭示了两类方法共同的软肋：当 $z\in\sigma(\Delta TA)$ 逼近虚轴（弱扩散、对流主导）时 $|R_g(z)|\to 1$，分母趋零、因子爆炸。这正是随后实验中黏性变弱、收敛急剧恶化的解析原因。

![原论文 Figure 4.8：等细求解成本下 MGRiT 与两轮 Parareal 的复平面收敛域](assets/papers/time-parallelization/source-figures/figure-4-8.svg)

Figure 4.8 用后向 Euler 粗传播 $R_g(z)=1/(1-z)$ 与精确细传播 $R_f(z)=e^z$，比较一轮 FCF-MGRiT 和两轮 Parareal。上排 (a) 是 $\varrho_{l,\mathrm{MGRiT}}\le\widehat\varrho$ 的区域，下排 (b) 是按相同两次细求解成本绘制的 $\varrho_{l,\mathrm{Parareal}}^2\le\widehat\varrho$ 区域；三列从左到右对应 $\widehat\varrho=0.2,0.4,0.6$。每一列的着色区域轮廓都很接近，说明**公平比较的单位应是细传播次数**：一轮 FCF-MGRiT（两次细解）与两轮 Parareal（也是两次细解）在复平面上的收敛域大致相当。换用 (4.9) 中的两种 SDIRK 细传播，结论相同。

### Theorem 4.6：等成本常数比较

Wu and Zhou (2019) 对 $z\in\mathbb R^-$ 给出更定量的比较。若细传播 L-稳定且 $J=\Delta T/\Delta t=O(1)$，后向 Euler 粗传播下负实轴上的最坏因子约为

$$
\max_{z\ge0}\varrho_l\approx
\begin{cases}
0.2984,&\text{Parareal},\\
0.1115,&\text{FCF-MGRiT}.
\end{cases}
$$

改用二阶 Lobatto IIIC（LIIIC-2）粗传播则为

$$
\max_{z\ge0}\varrho_l\approx
\begin{cases}
0.0817,&\text{Parareal},\\
0.0197,&\text{FCF-MGRiT}.
\end{cases}
$$

按两次细求解配平后比较：一轮 MGRiT（两次细解）应与两轮 Parareal（两次细解）对齐，即拿 $\varrho_{l,\mathrm{Parareal}}^2$ 去比 $\varrho_{l,\mathrm{MGRiT}}$。结果是

$$
0.2984^2=0.0890<0.1115,\qquad
0.0817^2=0.0067<0.0197.
$$

因此在这两种 L-稳定粗传播下，一轮 FCF-MGRiT 反而**略慢于**两轮 Parareal。

> [!tip] 本站洞见
> 这一结论常被误读为“MGRiT 不如 Parareal”。真正的含义是：**单看每次细求解带来的收缩，FCF 的额外 F 松弛并不自动占优**。原因在于 FCF 单轮曲线更陡——它在一轮内就把误差降到 Parareal 两轮才达到的量级，收敛更“成块”；但当把成本严格摊到每次细解时，Parareal 连续两次独立校正的复合略优于 MGRiT 一次带重叠的校正。选择哪种，取决于实际瓶颈是迭代轮数（同步、通信开销大 → 偏向 MGRiT 的少轮数）还是细求解总量（→ 二者相当，Parareal 常数略好）。

### Figures 4.9–4.11：线性与非线性实验

线性实验对热方程取齐次 Dirichlet 边界、对 ADE 取周期边界，共用

$$
u_0(x)=\sin^2(8\pi(1-x)^2),\quad
T=5,\quad J=20,\quad \Delta T=1/8,\quad \Delta x=1/160,
$$

粗层用后向 Euler，细层用 (4.9) 的 SDIRK22。

![原论文 Figure 4.9：热方程与两种黏性 ADE 上的模态因子分布](assets/papers/time-parallelization/source-figures/figure-4-9.svg)

Figure 4.9 上排 (a) 是 MGRiT，下排 (b) 是 Parareal；三列依次对应热方程、$\nu=0.1$ 的 ADE 和 $\nu=0.01$ 的 ADE，画的是 $\varrho_l(J,z)$ 在 $z\in\sigma(\Delta TA)$ 上的分布。图中标注的最大因子 $\varrho_{l,\max}=\max_{z}\varrho_l$ 分别为 MGRiT 的 $0.08375,0.2718,0.9021$，以及 Parareal 的 $0.2822,0.4453,0.9986$。前两列中，MGRiT 因子大致等于 Parareal 因子的平方（$0.2822^2\approx0.0796\approx0.08375$，$0.4453^2\approx0.198\lesssim0.2718$），印证 Theorem 4.6 所说“一轮 MGRiT ≈ 两轮 Parareal 且成本翻倍”；第三列两者都接近 $1$，粗传播已难以表示弱扩散下长寿命的对流模态。

![原论文 Figure 4.10：按两次细求解配平后的 Parareal 与 MGRiT 实测误差](assets/papers/time-parallelization/source-figures/figure-4-10.svg)

Figure 4.10 的四个面板按“热方程、$\nu=0.1$、$\nu=0.01$、$\nu=0.002$”排列，Parareal 横坐标把每两轮当作一轮（对齐两次细解），点划线标出截断误差量级 $\max\{\Delta t^2,\Delta x^2\}$（实际不会迭代到其之下）。前两幅中两条曲线接近；$\nu=0.01$ 时都很慢，且 Parareal 退化更明显——因为它在**每次**细传播后都紧跟一次已经失真的粗传播（见 (4.1)），而 FCF-MGRiT 连做两次细解、中间不插粗传播（见 (4.12)），从而少受粗层失真拖累。$\nu=0.002$ 时两者都发散，最大模态因子分别为 Parareal 的 $1.4211$ 和 MGRiT 的 $1.2812$。因此 Figure 4.9 的逐模态谱指标与 Figure 4.10 的实测曲线在四种扩散强度上一一对应。

非线性 Burgers 实验使用齐次 Dirichlet 边界、同一初值、$T=5$、$\Delta T=1/16$、$\Delta x=1/160$、$J=10$，空间用中心差分，粗层后向 Euler、细层 SDIRK22。非线性情形的收敛分析见 Gander et al. (2018b)，在对 $\mathcal G,\mathcal F$ 及其差满足一定 Lipschitz 条件下证明：只要粗传播足够准确，一轮 FCF-MGRiT（两次细解）的收缩与两轮 Parareal（两次细解）相当。

![原论文 Figure 4.11：三种黏性 Burgers 方程上等细求解成本的比较](assets/papers/time-parallelization/source-figures/figure-4-11.svg)

三个面板从左到右取 $\nu=0.5,0.01,0.002$，同样把 Parareal 两轮画作一轮以对齐细解次数。前两幅中，一轮 FCF-MGRiT 的下降大致对应两轮 Parareal；$\nu=0.002$ 时两条曲线都先经历一段较长的缓慢阶段、随后才进入快速下降。黏性降低使两种方法同步恶化，这与非线性 Lipschitz 分析的结论一致，也与线性谱分析中“弱扩散、对流主导 → $|R_g|\to1$ → 因子逼近 1”的机制相呼应。

## 公式、定理与图表覆盖核对

| 原文项目                        | 论文小节 | 覆盖状态                                                              |
| ------------------------------- | -------- | --------------------------------------------------------------------- |
| PFASST 来源（Minion 2010 等）   | 4.3      | SDC 替代细传播的溯源与多层 SDC 视角                                   |
| (4.10)                          | 4.3      | 细/粗配置节点、复制矩阵 $\chi$、配置系统 $\Phi$                       |
| Lagrange 传递与 PFASST 块迭代   | 4.3      | $T_{c\to f},T_{f\to c}$ 与 $B_{01},B_{10},B_{00}$ 的多重网格解读      |
| (4.11), Figure 4.6              | 4.3      | 隐式 Euler SDC 预处理器、Radau 节点、$Q_f/Q_c$、传递矩阵、PFASST 实验 |
| (4.12), Figure 4.7              | 4.4      | FCF 更新、两次细求解、重叠 $\Delta T$ 与有限步性质                    |
| (4.13), Theorem 4.5, Figure 4.8 | 4.4      | MGRiT 因子、$R_f^J$ 额外收缩因子、等成本收敛域                        |
| Theorem 4.6                     | 4.4      | 两类粗传播的四个最坏因子及平方（等成本）比较                          |
| Figures 4.9–4.11                | 4.4      | 线性谱、实测误差、非线性 Burgers 与弱扩散失效                         |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 4.3–4.4, pp. 452–460.

### 参考文献补充

- P. F. Emmett and M. L. Minion, _Toward an efficient parallel in time method for partial differential equations_, Commun. Appl. Math. Comput. Sci. **7** (2012), 105–132.（PFASST 的提出）
- M. L. Minion, _A hybrid parareal spectral deferred corrections method_, Commun. Appl. Math. Comput. Sci. **5** (2010), 265–301.（用一次 SDC 迭代替代细传播的先声）
- A. Dutt, L. Greengard and V. Rokhlin, _Spectral deferred correction methods for ordinary differential equations_, BIT **40** (2000), 241–266.（SDC 框架的起源）
- M. Bolten, D. Moser and R. Speck, _A multigrid perspective on the parallel full approximation scheme in space and time_, Numer. Linear Algebra Appl. **24** (2017), e2110；以及其收敛分析 (2018)。（把 PFASST 解释为时间多重网格）
- R. D. Falgout, S. Friedhoff, Tz. V. Kolev, S. P. MacLachlan and J. B. Schroder, _Parallel time integration with multigrid_, SIAM J. Sci. Comput. **36** (2014), C635–C661.（MGRiT 的提出）
- V. A. Dobrev, Tz. Kolev, N. A. Petersson and J. B. Schroder, _Two-level convergence theory for multigrid reduction in time (MGRiT)_, SIAM J. Sci. Comput. **39** (2017), S501–S527.（Theorem 4.5 的线性收敛估计）
