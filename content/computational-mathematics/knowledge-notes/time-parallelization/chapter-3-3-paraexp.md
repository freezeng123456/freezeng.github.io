---
title: 3.4：ParaExp
description: 完整推导线性 ParaExp 的精确叠加、非线性迭代和 Parareal 等价关系
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-3-paraexp
tags:
  - 时间并行
  - ParaExp
  - 矩阵指数
---

> [!note] 阅读范围
> 本页对应论文 Section 3.4（pp. 412–415），覆盖公式 (3.13)–(3.21)、Theorem 3.4 和 Figures 3.7–3.8。线性重构的归纳证明和非线性版本的有限步性质均展开说明。

## 3.4 ParaExp

### 线性问题的两组并行子问题

ParaExp 面向

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol g(t),
\qquad \boldsymbol u(0)=\boldsymbol u_0.
$$

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

### 为什么蓝色长尾仍然便宜

齐次问题有闭式表达

$$
\boldsymbol w_n(t)=
\exp\!\left((t-T_{n-1})A\right)
\boldsymbol v_{n-1}(T_{n-1}),
\qquad t\in[T_{n-1},T]. \tag{3.16}
$$

计算 $e^{\tau A}\boldsymbol b$ 可以直接跳到目标时刻，成本不必与中间步数
成正比。大型稀疏矩阵适合有理 Krylov 或 Chebyshev 展开；较小矩阵可用
scaling-and-squaring 加 Padé（Higham 2008；Moler 与 Van Loan 2003）。
MATLAB R2023b 及之后的 `expmv` 提供矩阵指数作用实现。
Schreiber、Peixoto、Haut 与 Wingate（2018）的 REXI 以及早期
Laplace 变换 PinT 技术也属于这类指数近似。

Gander 与 Güttel（2013）的波动方程 (2.7) 算例报告了高达 80% 的
时间并行效率，支持了一个明确判断：**ParaExp 很适合线性双曲问题**。
这个具体数字仍取决于指数算法、矩阵结构、分区和硬件。

### 非线性拆分

设非线性项能写成

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B(\boldsymbol u(t))+\boldsymbol g(t). \tag{3.17}
$$

非线性扩展来自 Gander、Güttel 与 Petcu（2018a）。仍令
$\boldsymbol u=\boldsymbol w+\boldsymbol v$，则保持原方程所需的精确拆分是

$$
\begin{aligned}
\boldsymbol w'(t)&=A\boldsymbol w(t),
&\boldsymbol w(0)&=\boldsymbol u_0,\\
\boldsymbol v'(t)&=A\boldsymbol v(t)
+B\!\left(\boldsymbol v(t)+\boldsymbol w(t)\right)+\boldsymbol g(t),
&\boldsymbol v(0)&=\boldsymbol 0.
\end{aligned}
$$

因为 $(\boldsymbol w+\boldsymbol v)'=
A(\boldsymbol w+\boldsymbol v)+B(\boldsymbol w+\boldsymbol v)+\boldsymbol g$。
与 (3.14) 不同，各时间区间不再解耦：$\boldsymbol w$ 在
$T_{n-1}$ 的初值依赖同一时刻的 $\boldsymbol v(T_{n-1})$。

> [!warning] 原文公式核对：非线性拆分
> 期刊版与 arXiv 版先漏掉了 $\boldsymbol v'$ 中的
> $A\boldsymbol v$，随后在未编号的区间迭代中又把
> $A\boldsymbol v_n^k$、$\boldsymbol v_n^k(T_{n-1})=0$
> 误排成 $A\boldsymbol u_n^k$、
> $\boldsymbol u_n^k(T_{n-1})=0$。这些写法与
> $\boldsymbol u=\boldsymbol v+\sum_j\boldsymbol w_j$ 不自洽；
> 本页按变量定义修正。

取 $\boldsymbol v_n^0(T_n)=0$（或其他接口近似），并从 $k=1$
开始，用上一轮接口值启动齐次问题。第一种并行迭代写成

$$
\begin{aligned}
(\boldsymbol w_n^k)'&=A\boldsymbol w_n^k,
&t&\in[T_{n-1},T],\\
\boldsymbol w_1^k(T_0)&=\boldsymbol u_0,\\
\boldsymbol w_n^k(T_{n-1})
&=\boldsymbol v_{n-1}^{k-1}(T_{n-1}),
&n&=2,\ldots,N_t,\\
(\boldsymbol v_n^k)'&=A\boldsymbol v_n^k\\
&\quad+B\!\left(
\boldsymbol v_n^k+\sum_{j=1}^n\boldsymbol w_j^k
\right)+\boldsymbol g,
&t&\in[T_{n-1},T_n],\\
\boldsymbol v_n^k(T_{n-1})&=\boldsymbol 0,\\
\boldsymbol u_n^k&=\boldsymbol v_n^k+\sum_{j=1}^n\boldsymbol w_j^k,
&n&=1,\ldots,N_t.
\end{aligned}
$$

若在每个局部问题中显式计算全部蓝色尾部
$\sum_j\boldsymbol w_j^k(t)$，大型 $A$ 会造成冗余。代入
$\boldsymbol v_n^k=\boldsymbol u_n^k-\sum_{j=1}^n\boldsymbol w_j^k$
可消去这项显式依赖。改写后的初始化取
$\boldsymbol u_n^0(T_n)=\boldsymbol w_j^0(T_n)=0$（所有相关
$j,n$），并从 $k=1$ 开始执行 (3.18)–(3.19)。其他初始化也可以，
但与 Parareal 比较时必须在粗节点上保持一致。

先对所有 $n=1,\ldots,N_t$ 构造齐次传播：

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

这些初值只含第 $k-1$ 轮数据，因此各个 $n$ 在第 $k$ 轮彼此独立，
可以并行求解。

> [!warning] 原文公式核对：`sequentially`
> 正式版在 (3.18) 前写成 “sequentially”，但 (3.18) 的右端只依赖
> 上一轮数据，与它前后两次 “in parallel” 的说明及公式依赖都冲突。
> 本页按算法的数据依赖把这一阶段解释为并行。

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

第 $k$ 轮的全局近似采用半开区间约定
$\boldsymbol u^k(t)=\boldsymbol u_n^k(t)$，
$t\in[T_{n-1},T_n)$；最终端点单独由粗节点值定义。

### Theorem 3.4：有限步收敛与 Parareal 等价

**有限步结论。** Theorem 3.4 引自 Gander 等（2018a）：第 $k$ 轮后，$\boldsymbol u^k(t)$ 在 $[0,T_k]$ 上与精确解相同，因此迭代在有限步内收敛。

> [!note] 本站补充：归纳理由
> 原文只给出结论。按迭代次数 $k$ 归纳可以看清机制：$k=1$ 时第一段从真实初值出发，故 $[0,T_1]$ 已精确；若第 $k-1$ 轮后 $[0,T_{k-1}]$ 精确，则 (3.18) 在 $T_{k-1}$ 处构造出精确初值，(3.19) 把精确区间再推进一段到 $T_k$。

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

标准 Parareal 的粗传播子通常也近似完整非线性问题 (3.20c)。
这里的 $\mathcal G$ 只保留 $A$，因此是一个简化版本。等价关系只在
粗节点成立：

$$
\boldsymbol u^k(T_n)=\boldsymbol u_{n+1}^k(T_n)=\boldsymbol U_n^k,
\qquad n=0,\ldots,N_t-1,
$$

且最终端点满足

$$
\boldsymbol u^k(T_{N_t})
:=\boldsymbol U_{N_t}^k,
\qquad \boldsymbol U_0^k=\boldsymbol u_0.
$$

一般而言，$\boldsymbol u_{N_t}^k(T_{N_t})$ 是最后一次局部细传播的
末值，并不等于校正后的 $\boldsymbol U_{N_t}^k$。整段局部轨迹也不与
Parareal 轨迹逐点相同。第四章再详细讨论标准 Parareal。

### Figure 3.8：Burgers 方程上的拆分失效

用周期 Burgers 方程检验这项非线性拆分：

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B\boldsymbol u^2(t),
\qquad t\in(0,2), \tag{3.21}
$$

这里 $\boldsymbol u^2$ 按分量计算，$A_{xx},A_x$ 是 (3.12) 的
无量纲 stencil。取 $\Delta x=1/100$ 后，一致的半离散系数为

$$
A=\frac{\nu}{\Delta x^2}A_{xx},
\qquad
B=-\frac{1}{4\Delta x}A_x.
$$

第二个系数来自
$-\frac12\partial_x(u^2)\approx-\frac1{4\Delta x}A_x\boldsymbol u^2$。
ParaExp 和标准 Parareal 的细传播子都使用后向 Euler，细步长为
$0.01/20$。标准 Parareal 的粗传播仍用后向 Euler，粗步长为
$0.01$；ParaExp 的线性粗传播调用 MATLAB `expmv`。

> [!warning] 原文公式核对：Figure 3.8 的空间缩放
> 期刊版与 arXiv 版写成 $A=A_{xx}$、$B=-A_x/2$，遗漏了
> $\nu$ 与 $\Delta x$ 的缩放。若照字面实现，$A$ 不含 $\nu$，
> 三组黏性实验便无法成立。上式由模型方程 (2.6) 和 (3.12) 的
> stencil 定义直接得到。

![原论文 Figure 3.8：三组黏性下非线性 ParaExp 与标准 Parareal 的误差](assets/papers/time-parallelization/source-figures/figure-3-8.svg)

三个面板从左到右取 $\nu=1,0.1,0.02$，横线表示离散截断误差
$\max\{\Delta t,\Delta x^2\}$，迭代达到这条线后即可停止。
$\nu=1$ 时，线性部分覆盖主导动力学，ParaExp 明显快于标准
Parareal；$\nu=0.1$ 时次序反转；$\nu=0.02$ 时 ParaExp 误差持续
增长，而标准 Parareal 仍能越过截断误差线。继续减小黏性后，
标准 Parareal 也会失效。关键不是统一的“扩散越弱、收敛越慢”，
而是 $A+B$ 拆分中主导动力学发生了转移；简化后的粗传播比标准
Parareal 更早失去代表性。

> [!important] 适用边界
> 线性 ParaExp 的 (3.15) 是精确代数分解，只受指数作用近似和局部受迫求解误差影响。非线性 ParaExp 已变成迭代法，收敛取决于 $A+B$ 的拆分。线性部分若没有覆盖主导传播机制，矩阵指数再精确也无法补偿模型失配。

## 公式与图表覆盖核对

| 原文项目                     | 论文小节 | 覆盖状态                                           |
| ---------------------------- | -------- | -------------------------------------------------- |
| (3.13)–(3.16), Figure 3.7    | 3.4      | 两组子问题、精确叠加、归纳证明、指数作用与完整原图 |
| (3.17)–(3.19)                | 3.4      | 非线性拆分、冗余来源和改写后的迭代                 |
| (3.20a)–(3.20c), Theorem 3.4 | 3.4      | 有限步结论、归纳理由及 Parareal 传播子             |
| (3.21), Figure 3.8           | 3.4      | Burgers 离散、全部时间参数、截断误差线与三阶段结论 |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3.4, pp. 412–415.
