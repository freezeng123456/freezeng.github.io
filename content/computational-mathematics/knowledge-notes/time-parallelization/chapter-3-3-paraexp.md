---
title: 3.4：ParaExp
description: 完整推导线性 ParaExp 的精确叠加、矩阵指数作用工具箱、非线性迭代与 Parareal 等价关系，并补充参考脉络
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-3-paraexp
tags:
  - 时间并行
  - ParaExp
  - 矩阵指数
---

> [!note] 阅读范围
> 本页对应论文 Section 3.4（pp. 412–415），覆盖公式 (3.13)–(3.21)、Theorem 3.4 和 Figures 3.7–3.8。线性重构的归纳证明、矩阵指数作用的算法工具箱、非线性版本的有限步性质及其与 Parareal 的等价关系均逐点展开。原始方法见 Gander and Güttel (2013)，非线性推广见 Gander et al. (2018a)。

## 3.4 ParaExp

### 线性问题的两组并行子问题

ParaExp 面向

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol g(t),
\qquad \boldsymbol u(0)=\boldsymbol u_0.
$$

ParaExp 属于**直接**时间并行方法：它不是靠迭代逐步逼近，而是把解一次性地代数重构出来（对比之下，Parareal、波形松弛等属于迭代类）。这类“把长时间区间拆开、并行求解再拼合”的思路可以上溯到最早的时间并行构想（Nievergelt 1964，参考补充）；ParaExp 的贡献是找到一种拆分，使拼合是**精确**的且拼合所需的传播足够便宜。

将时间区间切成 $[T_{n-1},T_n]$，$n=1,\ldots,N_t$。第一组“红色”问题保留源项，并在每个子区间使用零初值：

$$
\boldsymbol v_n'(t)=A\boldsymbol v_n(t)+\boldsymbol g(t),
\quad t\in(T_{n-1},T_n],
\qquad
\boldsymbol v_n(T_{n-1})=0. \tag{3.13}
$$

这些问题相互独立，可以全部并行。第二组“蓝色”问题删除源项，将前一段的红色末值向后传播到全局终点：

$$
\boldsymbol w_n'(t)=A\boldsymbol w_n(t),
\quad t\in(T_{n-1},T],
\qquad
\boldsymbol w_n(T_{n-1})=\boldsymbol v_{n-1}(T_{n-1}), \tag{3.14}
$$

其中 $\boldsymbol v_0(T_0)=\boldsymbol u_0$。所有蓝色问题也彼此独立。

![原论文 Figure 3.7：ParaExp 的红色局部受迫问题和蓝色齐次尾部传播](assets/papers/time-parallelization/source-figures/figure-3-7.svg)

Figure 3.7 的红线只覆盖本地子区间，蓝色虚线从每个接口一直延伸到 $T$。在第 $n$ 段上，精确解由本地红色响应和已经启动的全部蓝色响应叠加：

$$
\boldsymbol u(t)=\boldsymbol v_n(t)
+\sum_{j=1}^{n}\boldsymbol w_j(t),
\qquad t\in[T_{n-1},T_n]. \tag{3.15}
$$

![ParaExp 将局部受迫响应与全局齐次传播分开](assets/diagrams/pint/zh/paraexp-decomposition.svg)

> [!tip] 本站洞见
> 这两组子问题正是线性 ODE 解的“变差—常数”结构的并行化：红色 $\boldsymbol v_n$ 承担**受迫的特解部分**（右端 $\boldsymbol g$，零初值），蓝色 $\boldsymbol w_n$ 承担**齐次的初值传播部分**（无右端，携带接口值）。红色的零初值让各段互不依赖，从而立即并行；代价是每段丢掉了“来自过去的历史”，而这段历史恰恰由蓝色尾部补回。理解这一分工，就能看出 (3.15) 不是近似而是恒等式（下一节证明），也能解释为什么真正的开销集中在蓝色的长程传播上（后文说明它为何仍便宜）。

### 公式 (3.15) 的归纳证明

在第一段上，将 (3.13) 与 (3.14) 相加：

$$
(\boldsymbol v_1+\boldsymbol w_1)'
=A(\boldsymbol v_1+\boldsymbol w_1)+\boldsymbol g,
$$

且

$$
\boldsymbol v_1(0)+\boldsymbol w_1(0)
=0+\boldsymbol u_0.
$$

解的唯一性给出 $\boldsymbol u=\boldsymbol v_1+\boldsymbol w_1$，所以 (3.15) 对 $n=1$ 成立。

假设它对第 $n$ 段成立，则接口值满足

$$
\boldsymbol u(T_n)=\boldsymbol v_n(T_n)
+\sum_{j=1}^{n}\boldsymbol w_j(T_n).
$$

由 (3.14)，$\boldsymbol w_{n+1}(T_n)=\boldsymbol v_n(T_n)$，故

$$
\boldsymbol u(T_n)=\sum_{j=1}^{n+1}\boldsymbol w_j(T_n).
$$

令 $\boldsymbol w=\sum_{j=1}^{n+1}\boldsymbol w_j$。在 $(T_n,T_{n+1}]$ 上，$\boldsymbol w'=A\boldsymbol w$，初值为 $\boldsymbol u(T_n)$。再加上零初值受迫解 $\boldsymbol v_{n+1}$，便得到原初值问题在下一段的唯一解。归纳完成。

> [!tip] 本站洞见
> 归纳的每一步只用到两件事：**线性**（把两个子问题相加仍满足同一个方程）和**解的唯一性**（满足方程且初值正确者即为真解）。这说明 (3.15) 是一条**代数恒等式**，而非收敛意义下的极限——只要各子问题精确求解，叠加就一次到位，不存在“迭代次数”这一维度。归纳中反复用到的关键接口衔接由 $\boldsymbol w_{n+1}(T_n)=\boldsymbol v_n(T_n)$ 提供：第 $n$ 段红色积累出的历史，恰好成为第 $n+1$ 条蓝色尾部的“出发点”，于是全部历史被逐段接力送到区间末端。正因如此，线性 ParaExp 的误差只可能来自两处：矩阵指数作用的近似误差，以及红色受迫子问题的局部离散误差，而**不来自任何拆分假设**。

### 为什么蓝色长尾仍然便宜

齐次问题有闭式表达

$$
\boldsymbol w_n(t)=
\exp\!\left((t-T_{n-1})A\right)
\boldsymbol v_{n-1}(T_{n-1}),
\qquad t\in[T_{n-1},T]. \tag{3.16}
$$

乍看之下，越靠前的蓝色问题要积分到越靠后的 $T$，似乎和原问题一样贵。关键在于 (3.16) 是齐次的：所需的只是**矩阵指数作用** $e^{\tau A}\boldsymbol b$（$\tau=t-T_{n-1}$，$\boldsymbol b=\boldsymbol v_{n-1}(T_{n-1})$），它可以**直接跳到**目标时刻，成本由 $A$ 的谱性质和所需精度决定，而**不与中间时间步数成正比**。相比之下，普通时间步进的开销随区间长度线性累积，正是长尾看似昂贵的来源。这一“跳跃”能力把长程传播的代价与区间长度解耦，是 ParaExp 高并行效率的根本。

论文给出计算 $e^{\tau A}\boldsymbol b$ 的成熟工具（综述见 Higham 2008；Moler and Van Loan 2003），可按矩阵规模与结构选择：

- **有理 Krylov 方法**：面向大型稀疏 $A$。在有理 Krylov 子空间中逼近 $e^{\tau A}\boldsymbol b$，每步需求解带位移的线性系统 $(A-\sigma I)\boldsymbol x=\boldsymbol y$；对谱靠近负实轴的刚性/抛物型算子收敛尤其快。
- **Chebyshev 展开**：当 $A$ 对称（或谱落在已知实区间，如经中心差分离散的波动/扩散算子）时，用切比雪夫多项式一致逼近标量函数 $e^{\tau\lambda}$，只需反复做矩阵–向量乘，无需线性求解。
- **scaling-and-squaring 加 Padé 逼近**：形成**完整**矩阵指数 $e^{\tau A}$ 的稠密算法，适合较小的稠密矩阵。
- **`expmv` 型作用算法**：直接计算作用 $e^{\tau A}\boldsymbol b$ 而不显式形成整个指数矩阵（其底层思想见 Al-Mohy and Higham 2011，参考补充），MATLAB R2023b 及之后版本内置 `expmv`。
- **REXI / 基于 Laplace 变换的早期 PinT**：把 $e^{\tau A}\boldsymbol b$ 写成若干带位移预解式 $(A-\sigma_\ell I)^{-1}\boldsymbol b$ 的加权和，各项彼此独立、可再叠一层并行；对振荡（双曲）问题尤为契合（REXI 见 Schreiber et al. 2018）。

> [!tip] 本站洞见
> 选择哪种工具本质上是在“是否形成整个 $e^{\tau A}$”与“$A$ 的谱几何”之间权衡：稠密小矩阵可以负担 scaling-and-squaring 一次算全；大型稀疏矩阵则应避免生成稠密指数，改用只需矩阵–向量乘或位移求解的作用型方法（Krylov / Chebyshev / expmv）。谱的位置进一步细分：谱偏负实轴（强扩散）时 Krylov/Chebyshev 的多项式逼近收敛快；谱贴近虚轴（强振荡、波动）时 REXI 这类预解式求和更稳健。值得注意的是，REXI 的“预解式求和”本身又是一层可并行结构，于是 ParaExp 内部还能再嵌套并行——这正是它在双曲问题上能逼近高并行效率的技术底座。

Gander and Güttel (2013) 的波动方程实验曾报告约 80% 的时间并行效率。这个数字依赖指数算法、矩阵结构、分区和硬件，表达的是该实现的结果，并非方法的普适上界；换一套指数作用实现或不同的谱结构，效率会随之变化。

### 非线性拆分

设非线性项能写成

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B(\boldsymbol u(t))+\boldsymbol g(t). \tag{3.17}
$$

线性情形的直接叠加在这里失效，因为 $B$ 会把各段重新耦合。沿用线性思路的初始构造是把 $\boldsymbol u=\boldsymbol w+\boldsymbol v$，令齐次部分 $\boldsymbol w'=A\boldsymbol w$ 携带初值 $\boldsymbol w(0)=\boldsymbol u_0$，而非线性部分满足 $\boldsymbol v'=B(\boldsymbol v+\boldsymbol w)+\boldsymbol g$、零初值 $\boldsymbol v(0)=0$；两者之和仍解 (3.17)。问题在于 $B$ 项里出现 $\boldsymbol v+\boldsymbol w$，使各子区间通过接口值 $\boldsymbol v(T_{n-1})$ 相互牵连，无法像线性那样一次性并行，只能改成迭代。

若在每个非线性局部问题中显式计算全部蓝色尾部 $\sum_j\boldsymbol w_j^k(t)$，则每段都要把线性问题积分到整个 $[T_{n-1},T]$，大型 $A$ 会造成冗余与浪费。论文改用 $\boldsymbol v_n^k=\boldsymbol u_n^k-\sum_{j=1}^n\boldsymbol w_j^k$ 消去这一冗余，得到只需局部积分的两步迭代。

先按 $n=1,\ldots,N_t$ 构造齐次传播：

$$
\begin{aligned}
(\boldsymbol w_n^k)'(t)&=A\boldsymbol w_n^k(t),
&&t\in[T_{n-1},T],\\
\boldsymbol w_n^k(T_{n-1})
&=\boldsymbol u_{n-1}^{k-1}(T_{n-1})
-\sum_{j=1}^{n-1}\boldsymbol w_j^{k-1}(T_{n-1}),
&&\boldsymbol w_1^k(T_0)=\boldsymbol u_0.
\end{aligned} \tag{3.18}
$$

随后在所有时间子区间上并行求解完整非线性问题：

$$
\begin{aligned}
(\boldsymbol u_n^k)'(t)
&=A\boldsymbol u_n^k(t)+B(\boldsymbol u_n^k(t))+\boldsymbol g(t),
&&t\in[T_{n-1},T_n],\\
\boldsymbol u_n^k(T_{n-1})
&=\sum_{j=1}^{n}\boldsymbol w_j^k(T_{n-1}).
\end{aligned} \tag{3.19}
$$

第 $k$ 轮的全局近似在第 $n$ 段上取 $\boldsymbol u^k(t)=\boldsymbol u_n^k(t)$。

> [!tip] 本站洞见
> 换元 $\boldsymbol v_n^k=\boldsymbol u_n^k-\sum_{j=1}^n\boldsymbol w_j^k$ 的巧妙之处在于把“昂贵且冗余”的量替换掉：原始形式里 $B$ 依赖 $\sum_j\boldsymbol w_j^k$，迫使线性尾部在**每个**局部非线性问题中都被重复积分到区间末端；换元后，非线性问题 (3.19) 只在本地 $[T_{n-1},T_n]$ 上求解完整方程 $\boldsymbol u'=A\boldsymbol u+B(\boldsymbol u)+\boldsymbol g$，而线性尾部 (3.18) 独立地做一次性传播。注意 (3.18) 里蓝色初值取上一轮的 $\boldsymbol u_{n-1}^{k-1}-\sum_{j}\boldsymbol w_j^{k-1}$，正是被消掉的 $\boldsymbol v_{n-1}^{k-1}(T_{n-1})$；(3.19) 里非线性初值取当前轮的 $\sum_j\boldsymbol w_j^k(T_{n-1})$。于是线性传播与非线性求解被彻底解耦到两个可分别并行的步骤，这也是它下节能写成 Parareal 形式的结构基础。

### Theorem 3.4：有限步收敛与 Parareal 等价

**有限步结论。** 第 $k$ 轮后，$\boldsymbol u^k(t)$ 在 $[0,T_k]$ 上与精确解相同，即迭代 ParaExp 在有限步内收敛。理由可按时间段归纳：第一段总从真实初值 $\boldsymbol u_0$ 出发，故第 1 轮起就在 $[0,T_1]$ 上精确；若前 $k-1$ 轮已在 $[0,T_{k-1}]$ 上精确，则 (3.18) 在下一个接口构造出精确初值，(3.19) 把正确性再向前推进一段，于是第 $k$ 轮覆盖到 $T_k$。信息每轮沿时间方向前进一个粗区间，因此至多 $N_t$ 轮即得到精确解——这与迭代法的“渐近收敛”本质不同，是**有限步终止**。

> [!tip] 本站洞见
> 有限步性质其实是线性精确叠加 (3.15) 在非线性情形下的“残影”：非线性 $B$ 破坏了一次性叠加，但每轮迭代仍能把一个粗区间的解“钉死”为精确值，因此收敛是逐段推进而非全局衰减。这也解释了收敛曲线的形状——误差通常呈台阶式下降，直到第 $N_t$ 轮触底。实践中并不需要跑满 $N_t$ 轮：一旦误差降到离散截断误差量级即可停止（见下文 Figure 3.8 的横线）。

在粗节点 $T_n$ 上，该迭代等价于

$$
\boldsymbol U_n^k
=\mathcal G(T_{n-1},T_n,\boldsymbol U_{n-1}^k)
+\mathcal F(T_{n-1},T_n,\boldsymbol U_{n-1}^{k-1})
-\mathcal G(T_{n-1},T_n,\boldsymbol U_{n-1}^{k-1}), \tag{3.20a}
$$

其中粗传播子只解线性齐次问题

$$
\boldsymbol u'=A\boldsymbol u,
\qquad
\boldsymbol u(T_{n-1})=\boldsymbol U,
\qquad t\in[T_{n-1},T_n], \tag{3.20b}
$$

细传播子解完整非线性问题

$$
\boldsymbol u'=A\boldsymbol u+B(\boldsymbol u)+\boldsymbol g,
\qquad
\boldsymbol u(T_{n-1})=\boldsymbol U,
\qquad t\in[T_{n-1},T_n]. \tag{3.20c}
$$

这是本文首次出现 Parareal（第 4 节详述）。标准 Parareal 的粗传播子 $\mathcal G$ 通常也**近似求解完整非线性问题 (3.20c)**（只是用更粗的步长），因此廉价而仍抓住主要动力学。这里的 $\mathcal G$ 只保留线性部分 $A$、丢掉 $B$（即求解 (3.20b) 而非 (3.20c)），是一个**简化版本**。换言之，ParaExp 的“粗模型”是原问题的线性化，而非粗离散化。

> [!tip] 本站洞见
> 把非线性 ParaExp 认成“$\mathcal G$ 只解线性、$\mathcal F$ 解完整非线性”的 Parareal，能一眼看清它的成败关键：Parareal 的收敛依赖粗传播子对细传播子的**逼近质量**。标准 Parareal 用粗步长离散**同一个**非线性算子，误差主要来自时间分辨率；而 ParaExp 的 $\mathcal G$ 干脆扔掉 $B$，误差来自**模型失配**——只要 $B$ 携带了主要传播机制（如强对流），线性 $\mathcal G$ 就与 $\mathcal F$ 相去甚远，修正项 $\mathcal F-\mathcal G$ 不再是小量，迭代随之变慢乃至发散。因此“$A+B$ 的拆分是否把主导动力学放进 $A$”直接决定非线性 ParaExp 的表现，这一点在下面的 Burgers 实验里表现得非常直观。论文也提前提示：标准 Parareal 本就不擅长双曲问题，简化版更难指望在该情形奏效。

### Figure 3.8：Burgers 方程上的拆分失效

论文使用

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B\boldsymbol u^2(t),
\qquad t\in(0,2), \tag{3.21}
$$

它来自周期 Burgers 方程的中心差分，$\Delta x=1/100$，$A=\nu A_{xx}/\Delta x^2$，$B=-A_x/(2\Delta x)$，$A_{xx},A_x$ 见 (3.12)。这里线性部分 $A$ 对应黏性扩散项、由黏度 $\nu$ 标定，非线性部分 $B\boldsymbol u^2$ 对应对流项。ParaExp 和标准 Parareal 的细传播子都使用后向 Euler，细步长为 $0.01/20$。标准 Parareal 的粗传播子仍用后向 Euler，粗步长为 $0.01$；ParaExp 的线性粗传播调用 MATLAB `expmv`。

![原论文 Figure 3.8：三组黏性下非线性 ParaExp 与标准 Parareal 的误差](assets/papers/time-parallelization/source-figures/figure-3-8.svg)

三个面板从左到右取 $\nu=1,0.1,0.02$，横线表示离散截断误差 $\max\{\Delta t,\Delta x^2\}$，实际迭代达到这条线后即可停止。$\nu=1$（强扩散）时，$A$ 捕获主要动力学，ParaExp 明显快于标准 Parareal；$\nu=0.1$ 时，标准 Parareal 反而更快，ParaExp 仍缓慢下降；$\nu=0.02$（对流主导）时，**ParaExp 率先失效**——其误差连续增大即发散，而标准 Parareal 仍能越过截断误差线。继续减小黏性后，标准 Parareal 最终也会失效（详见第 4 节）。这三个面板展示的是拆分主导项随黏性变化发生转移，不能只概括为“扩散减弱后收敛变慢”。

> [!tip] 本站洞见
> 三个 $\nu$ 恰好扫过“$A$ 主导 → 势均力敌 → $B$ 主导”的相变：$\nu=1$ 时扩散（在 $A$ 里）主导，ParaExp 的线性粗模型几乎抓住全部动力学，因而最快；$\nu=0.02$ 时对流（在 $B$ 里）主导，被 ParaExp 的 $\mathcal G$ 完全忽略，于是它**比标准 Parareal 更早发散**——因为标准 Parareal 的粗传播子毕竟离散了含对流的完整方程，只是分辨率粗。这印证了上一节的判断：谁把主导项放进了粗模型，谁就在对应区间胜出。工程含义是，若要把 ParaExp 用于对流—扩散型问题，应设法把对流也纳入 $A$（例如围绕线性化算子做 $A+B$ 的再拆分），否则再精确的矩阵指数也补不回模型失配。

> [!important] 适用边界
> 线性 ParaExp 的 (3.15) 是精确代数分解，只受指数作用近似和局部受迫求解误差影响。非线性 ParaExp 已变成迭代法，收敛取决于 $A+B$ 的拆分。线性部分若没有覆盖主导传播机制，矩阵指数再精确也无法补偿模型失配。

## 公式与图表覆盖核对

| 原文项目                     | 论文小节 | 覆盖状态                                                 |
| ---------------------------- | -------- | -------------------------------------------------------- |
| (3.13)–(3.16), Figure 3.7    | 3.4      | 两组子问题、精确叠加、归纳证明、指数作用工具箱与完整原图 |
| (3.17)–(3.19)                | 3.4      | 非线性拆分、初始解耦、冗余来源与换元后的两步迭代         |
| (3.20a)–(3.20c), Theorem 3.4 | 3.4      | 有限步结论、归纳理由及与标准 Parareal 粗传播子的差异     |
| (3.21), Figure 3.8           | 3.4      | Burgers 离散、全部时间参数、截断误差线与三阶段相变结论   |

## 延伸文献（参考补充）

- 直接时间并行的最早构想：J. Nievergelt, _Parallel methods for integrating ordinary differential equations_, _Comm. ACM_ 7 (1964)。
- 线性 ParaExp 原始方法与波动方程效率：M. J. Gander and S. Güttel, _ParaExp: A parallel integrator for linear initial-value problems_, _SIAM J. Sci. Comput._ 35 (2013)。
- 矩阵指数作用综述与算法：N. J. Higham, _Functions of Matrices_ (SIAM 2008)；C. Moler and C. Van Loan, _Nineteen dubious ways to compute the exponential of a matrix, twenty-five years later_, _SIAM Review_ 45 (2003)；A. Al-Mohy and N. J. Higham, _Computing the action of the matrix exponential_, _SIAM J. Sci. Comput._ 33 (2011)。
- 指数型 PinT 与预解式求和：M. Schreiber, P. Peixoto, T. Haut and B. Wingate (2018)（REXI）。

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3.4, pp. 412–415.
