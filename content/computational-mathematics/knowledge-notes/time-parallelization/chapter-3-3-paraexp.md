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

计算 $e^{\tau A}\boldsymbol b$ 可以直接跳到目标时刻，成本不必与中间步数成正比。论文列出以下工具：大型稀疏矩阵适合有理 Krylov 和 Chebyshev 展开；较小矩阵可以用 scaling-and-squaring 加 Padé；MATLAB R2023b 及之后的 `expmv` 提供矩阵指数作用实现。REXI 和基于 Laplace 变换的早期 PinT 技术也属于同一类指数近似思路。

Gander and Güttel (2013) 的波动方程实验曾报告约 80% 的时间并行效率。这个数字依赖指数算法、矩阵结构、分区和硬件，表达的是该实现的结果。

### 非线性拆分

设非线性项能写成

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B(\boldsymbol u(t))+\boldsymbol g(t). \tag{3.17}
$$

线性情形的直接叠加在这里失效，因为 $B$ 会把各段重新耦合。初始构造把 $\boldsymbol u=\boldsymbol w+\boldsymbol v$，令 $\boldsymbol w'=A\boldsymbol w$，并在 $\boldsymbol v$ 方程中使用 $B(\boldsymbol v+\boldsymbol w)$。为了并行，需要用上一轮的接口值启动蓝色问题。

若在每个非线性局部问题中显式计算全部蓝色尾部 $\sum_j\boldsymbol w_j^k(t)$，大型 $A$ 会造成冗余。论文改用 $\boldsymbol v_n^k=\boldsymbol u_n^k-\sum_{j=1}^n\boldsymbol w_j^k$，得到两步迭代。

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

### Theorem 3.4：有限步收敛与 Parareal 等价

**有限步结论。** 第 $k$ 轮后，$\boldsymbol u^k(t)$ 在 $[0,T_k]$ 上与精确解相同。理由可按时间段归纳：第一段总从真实初值出发；若前 $k-1$ 段已经精确，(3.18) 在下一个接口构造出精确初值，(3.19) 就把正确性再推进一段。

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

标准 Parareal 的粗传播子通常也近似完整非线性问题。这里的 $\mathcal G$ 只保留 $A$，因此是一个简化版本。拆分是否抓住主要动力学，决定非线性 ParaExp 的表现。

### Figure 3.8：Burgers 方程上的拆分失效

论文使用

$$
\boldsymbol f(\boldsymbol u(t),t)
=A\boldsymbol u(t)+B\boldsymbol u^2(t),
\qquad t\in(0,2), \tag{3.21}
$$

它来自周期 Burgers 方程的中心差分，$\Delta x=1/100$，$A=\nu A_{xx}/\Delta x^2$，$B=-A_x/(2\Delta x)$，$A_{xx},A_x$ 见 (3.12)。ParaExp 和标准 Parareal 的细传播子都使用后向 Euler，细步长为 $0.01/20$。标准 Parareal 的粗传播子仍用后向 Euler，粗步长为 $0.01$；ParaExp 的线性粗传播调用 MATLAB `expmv`。

![原论文 Figure 3.8：三组黏性下非线性 ParaExp 与标准 Parareal 的误差](assets/papers/time-parallelization/source-figures/figure-3-8.svg)

三个面板从左到右取 $\nu=1,0.1,0.02$，横线表示离散截断误差 $\max\{\Delta t,\Delta x^2\}$，实际迭代达到这条线后即可停止。$\nu=1$ 时，$A$ 捕获主要动力学，ParaExp 明显快于标准 Parareal；$\nu=0.1$ 时，标准 Parareal 反而更快，ParaExp 仍缓慢下降；$\nu=0.02$ 时，ParaExp 的误差连续增大，标准 Parareal 仍能越过截断误差线。继续减小黏性后，标准 Parareal 最终也会失效。这三个面板展示的是拆分主导项随黏性变化发生转移，不能只概括为“扩散减弱后收敛变慢”。

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
