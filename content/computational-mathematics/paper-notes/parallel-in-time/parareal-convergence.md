---
title: parareal 的收敛分析
description: 编号 12、20、30、77：收缩因子对细传播子稳定性与细网格结构的依赖
lang: zh
translation: en/computational-mathematics/paper-notes/parallel-in-time/parareal-convergence
tags:
  - 论文笔记
  - 时间并行
  - 收敛分析
---

> [!note] 本页覆盖
> 编号 **12**（_SIAM J. Sci. Comput._ 37(2), 2015）、**20**（_J. Comput. Phys._ 329, 2017）、**30**（_J. Comput. Phys._ 358, 2018）、**77**（_SIAM J. Numer. Anal._ 62(5), 2024）。这四篇均无预印本，编号 20、77 的正文与摘要未能通过公开渠道完整获取，相应小节明确标出可核实范围。

## 收敛分析的标准装置

考虑由抛物 PDE 半离散得到的对称正定系统 $\boldsymbol u'(t)+A\boldsymbol u(t)=g(t)$。parareal 迭代为

$$
\boldsymbol u_{n+1}^{k+1}=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^{k})
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k}),
$$

每个粗区间含 $J=\Delta T/\Delta t\ge2$ 个细步。若 $A=V_ADV_A^{-1}$，且粗、细传播子是稳定函数为 $R_g,R_f$ 的一步法，则误差满足 $\boldsymbol\xi^{k+1}(z)=M(z)\boldsymbol\xi^k(z)$，其中

$$
M(z)=M_g^{-1}(z)\bigl[M_g(z)-M_f(z)\bigr]=I_t-M_g^{-1}(z)M_f(z),
$$

$M_g(z)$ 与 $M_f(z)$ 是单位对角、次对角分别为 $-R_g(z)$ 与 $-R_f^{J}(z/J)$ 的下三角 Toeplitz 矩阵，$z=\Delta T\lambda$。把标量因子提出来得到

$$
M(z)=\bigl[R_f^J(z/J)-R_g(z)\bigr]\,\widetilde M\bigl(R_g(z)\bigr),
$$

其中 $\widetilde M(\beta)$ 是次对角依次为 $1,\beta,\beta^2,\dots$ 的严格下三角 Toeplitz 矩阵。线性收敛因子与超线性估计因此是

$$
\varrho_l(J,z)=\frac{\bigl|R_g(z)-R_f^{J}(z/J)\bigr|}{1-|R_g(z)|},
\qquad
\varrho_s(J,z,N_t,k)=\frac{\bigl|R_g(z)-R_f^J(z/J)\bigr|^k}{k!}\prod_{j=1}^k(N_t-j),
$$

后者在 $k=N_t$ 处为零，对应精确算术下的有限步终止。

**这套装置里有一个被默认的结构假设，后面编号 77 正是要拆掉它。** $\mathcal F$ 在一个粗区间上的作用被写成单个标量稳定函数在单个自变量上的 $J$ 次**幂** $R_f^{J}(z/J)$，这一步要求区间内 $J$ 个细步等长。

## 12：细传播子只是 A-稳定时会发生什么

### 此前理论覆盖的与未覆盖的

理解得最清楚的情形是 parareal-Euler：粗细都用后向 Euler，Mathew、Sarkis 与 Schaerer 证明了与网格和 $T$ 无关、对所有 $J\ge2$ 成立、约为 $0.298$ 的收缩因子。

未被理解的是 $\mathcal F$ 取更高阶但**只是 A-稳定而非 L-稳定**的积分器，最重要的是梯形规则与四阶两级 Gauss Runge-Kutta 方法。以 $\boldsymbol u'+A\boldsymbol u=g$ 的约定（$z=\Delta T\lambda$，$\lambda>0$），四个相关积分器的稳定函数为

$$
\text{后向 Euler：}\ R(z)=\frac{1}{1+z},
\qquad
\text{梯形：}\ R(z)=\frac{2-z}{2+z},
\qquad
\text{四阶 Gauss：}\ R(z)=\frac{12-6z+z^2}{12+6z+z^2}.
$$

后向 Euler 是 L-稳定的，$R(\infty)=0$；梯形规则 A-稳定但不 L-稳定，$R(\infty)=-1$；四阶 Gauss 方法 A-稳定且辛，$R(\infty)=+1$。

后两者的 $|R_f(z)|\to1$（$z\to\infty$），因此产生 $\approx0.3$ 常数的论证失效：$\varrho_l$ 的分子不再随高频衰减，$\lambda_{\max}$ 很大时收缩因子可以趋近 $1$。此前的分析要么假设 $\mathcal F$ 是精确传播子 $\exp(\Delta TA)$，要么假设 L-稳定，因此对实践中重要的仅 A-稳定情形无话可说，尤其在粗细比 $J$ 较小时。

这篇论文的贡献正是给出这一情形下的临界粗细比框架：需要多大的 $J$ 才能把高频行为压回可控范围。相关结论在[[computational-mathematics/knowledge-notes/time-parallelization/index|时间并行综述精读]]的第四章中以 $(4.8)$ 的形式出现，那里给出 $J_{\min}=O(\log^2 z_{\max})$ 的有限谱扩展。

## 20 与 30：两种「分数阶」带来两种不同的困难

这两篇常被并列提及，但它们处理的是**结构上不同**的困难。

### 20：空间分数阶——粗传播子成为瓶颈

设定是含分数阶 Laplacian $(-\Delta)^{\alpha}$ 的时间依赖 PDE，空间半离散后得到刚性系统 $\boldsymbol u'+A\boldsymbol u=g$。此时 $A$ 继承算子的非局部性，因此是**稠密**的（一维情形有 Toeplitz 结构），且 $\lambda_{\max}$ 极大。两个困难叠加：每次与 $A$ 的隐式求解远比经典 Laplacian 情形昂贵，因此顺序执行、直接位于关键路径上的**粗**传播子 $\mathcal G$ 主导运行时间；而极大的 $\lambda_{\max}$ 正是编号 12 中临界粗细比按 $O(\log^2(\Delta T\lambda_{\max}))$ 增长的区域。题名中的「快速」指的是加速 parareal **迭代本身**，而不只是给出收缩因子的界。

> [!warning] 可核实范围
> 该文摘要未能从 Crossref、OpenAlex、Semantic Scholar、NASA ADS 或出版社页面获取，正文亦不可得。可确认的是：第三方文献明确记载该文在**常步长**下分析空间分数阶扩散问题的 parareal 收敛性质。其具体的分数阶离散方式、定理、收缩因子与数值算例本站均未核实，因此不报告任何数值常数。

### 30：时间分数阶——负载不均衡

时间分数阶方程的困难是结构性的：Caputo 或 Riemann-Liouville 导数在时刻 $t$ 依赖整个 $[0,t]$ 上的解历史。而 parareal 的全部前提是 $[T_n,T_{n+1}]$ 上的细传播只需要局部初值 $\boldsymbol u_n^k$。算子在时间上非局部后，在第 $n$ 个子区间上推进还需要 $[0,T_n]$ 上的历史，因此朴素推广会让处理器 $n$ 的工作量正比于 $n$：**计算时间在进程间不均衡**，最后一个处理器的工作量是第一个的 $N_t$ 倍，即使迭代收敛，并行效率也被破坏。

论文的对策是用两个局部时间积分器把分数阶算子**局部化**：引入辅助变量，其演化由局部（常微分）动力学支配，从而把非局部历史编码进有限个额外的局部未知量，而不是一个不断增长的卷积和。写成关于（解变量，辅助变量）的增广**局部**系统后，每次细传播调用在各子区间上代价相同，负载因此均衡。

真正新的算法成分是一种**混合式粗网格校正**：辅助变量与解变量**分开**校正，而不是把单一的 parareal 更新式作用在整个增广状态向量上。论文报告该算法具有稳健的收敛速率——在这一系列工作中，「稳健」一贯指与空间算子的特征值、粗细比 $J$ 以及粗区间数 $N_t$ 无关。

> [!note] 可核实范围
> 上述问题设定、局部化思路与混合校正的**存在**均可从摘要确认。所采用的两个局部时间积分器的具体名称、混合校正的确切公式、收敛因子的显式值与假设，以及负载均衡是否被定量化，本站均未核实。

## 77：拆掉「细网格等长」这个假设

### 为什么这不是常规推广

该文研究**均匀粗网格**配**任意分布的非均匀细网格**时的 parareal 收敛行为。「任意分布」是实质所在：粗区间内的细网格既不假设均匀，也不假设分级，甚至不假设各粗区间之间相同。

回到本页开头的装置就能看出困难在哪。所有经典线性收敛结果——Gander 与 Vandewalle 的界、$0.3$ 附近的常数、编号 12 的临界粗细比公式、编号 39 的 MGRIT 因子——都依赖恒等式

$$
\mathcal F(T_n,T_{n+1},\cdot)\ \longleftrightarrow\ R_f^{J}\!\Bigl(\frac{z}{J}\Bigr),
$$

即粗区间上的细传播是单个标量稳定函数在单个自变量上的 $J$ 次幂。细网格非均匀时得到的是**互不相同因子的乘积**：

$$
\mathcal F(T_n,T_{n+1},\cdot)\ \longleftrightarrow\ \prod_{i=1}^{J_n}R_f(\theta_{n,i}z),
\qquad
\theta_{n,i}=\frac{\Delta t_{n,i}}{\Delta T},\quad \sum_i\theta_{n,i}=1 .
$$

两件事随之失效。第一，$\varrho_l(J,z)$ 不再有定义，因为没有单一的 $J$。第二，也更严重：若细网格**在各粗区间之间不同**，则 $M_f(z)$ 的次对角元素逐行变化，$M_f(z)$ **不再是 Toeplitz** 矩阵，因式分解 $M(z)=[R_f^J(z/J)-R_g(z)]\widetilde M(R_g(z))$ 失效，而 Gander 与 Vandewalle 关于 $\|\widetilde M^k\|_\infty$ 的引理——这条文献中每个收敛因子结果的引擎——无法应用。在这一设定下恢复稳健的、与网格无关的收缩因子因此是全新的分析，而不是常规推广。

想要非均匀细网格的理由是明确的：自适应时间步、初值不光滑或不相容时 $t=0$ 附近的分级网格、快速瞬态附近的局部细化，都会产生非均匀细网格。把 parareal 理论限制在均匀细网格上，等于排除了几乎全部自适应实践。

> [!warning] 可核实范围
> 该文摘要未能获取：SIAM 不向 Crossref 提交摘要，OpenAlex、Semantic Scholar 与 NASA ADS 对该 DOI 均返回空值。可确认的是上述问题设定（均匀粗网格加任意分布的非均匀细网格），以及该文在文献中的定位——其参考文献列表是经典 parareal 收敛理论的阅读清单（编号 12、Gander 与 Hairer、Bal、Gander 与 Vandewalle、Southworth 等），因此这是一篇纯分析工作，而不是新算法工作。本站不报告其定理、收敛因子或门槛。
>
> 用于对照的量级是：均匀细网格、L-稳定 $\mathcal F$、后向 Euler 粗传播下的 $\approx0.3$ 常数。该文在非均匀细网格上是否达到可比的常数，本站未核实。

## 一处值得指出的张力

非均匀时间步在这两条路线上的地位恰好相反。在 parareal 分支中，它是**理论的障碍**：它破坏 Toeplitz 结构，从而破坏全部收敛因子论证。而在直接对角化分支中（见[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|对角化技术]]），非均匀步长是**必要条件**：Maday 与 Rønquist 需要各 $\Delta t_n$ 互不相同，才能让时间矩阵可对角化。同一个建模自由度，在一条路线上是麻烦，在另一条路线上是前提。

## 覆盖核对

| 内容                             | 论文   | 覆盖状态                                     |
| -------------------------------- | ------ | -------------------------------------------- |
| parareal 迭代与误差传播矩阵      | 12     | 迭代式、$M(z)$、因式分解、两个收敛因子       |
| 四个积分器的稳定函数与其极限行为 | 12     | 后向 Euler、梯形、四阶 Gauss 的 $R(\infty)$  |
| 仅 A-稳定情形为何超出旧理论      | 12     | 分子不衰减与临界粗细比的必要性               |
| 空间分数阶下粗传播子成为瓶颈     | 20     | 稠密 $A$、极大 $\lambda_{\max}$（限定核实）  |
| 时间分数阶的历史效应与负载不均衡 | 30     | 前提失效、工作量正比于 $n$                   |
| 局部化与混合式粗网格校正         | 30     | 辅助变量、增广局部系统、分开校正（限定核实） |
| 非均匀细网格破坏 Toeplitz 结构   | 77     | 幂到乘积、两处失效、为何是新分析             |
| 非均匀步长在两条路线上的相反地位 | 12, 77 | 障碍与前提的对照                             |

## 本页原文

- S. Wu and T. Zhou, [_Convergence analysis for three parareal solvers_](https://doi.org/10.1137/140970756), SIAM J. Sci. Comput. 37(2) (2015), pp. A970-A992。
- S. Wu and T. Zhou, [_Fast parareal iterations for fractional diffusion equations_](https://doi.org/10.1016/j.jcp.2016.10.046), J. Comput. Phys. 329 (2017), pp. 210-226。
- S. Wu and T. Zhou, [_Parareal algorithms with local time-integrators for time fractional differential equations_](https://doi.org/10.1016/j.jcp.2017.12.029), J. Comput. Phys. 358 (2018), pp. 135-149。
- S.-L. Wu and T. Zhou, [_Convergence analysis of the parareal algorithm with nonuniform fine time grid_](https://doi.org/10.1137/23M1592481), SIAM J. Numer. Anal. 62(5) (2024), pp. 2308-2330。
