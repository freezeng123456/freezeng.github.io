---
title: 随机 Galerkin 与随机配点
description: 编号 1、2、3、4、5、7、38：耦合系统的代数性质与随机变量方向的正则性
lang: zh
translation: en/computational-mathematics/paper-notes/stochastic-approximation/stochastic-galerkin-and-collocation
tags:
  - 论文笔记
  - 不确定性量化
  - 随机 Galerkin
---

> [!note] 本页覆盖
> 编号 **1**（_J. Comput. Phys._ 229, 2010）、**2**（_Commun. Comput. Phys._ 8, 2010）、**3**（_J. Comput. Appl. Math._ 236, 2011）、**4**（_J. Comput. Math._ 30, 2012）、**5**（_J. Sci. Comput._ 51, 2012）、**7**（_Adv. Appl. Math. Mech._ 6, 2014）、**38**（_J. Sci. Comput._ 80, 2019）。

这七篇属于**侵入式**路线：把解按广义多项式混沌展开、再做 Galerkin 投影，或在参数点上解耦地求解。共同的技术问题有两个：耦合系统的代数性质，以及解在**随机变量方向**的正则性。

## 1：耦合系统的系数矩阵何时对角占优

设随机场仿射依赖参数 $\kappa(x,y)=\kappa_0(x)+\sum_{i=1}^{N}\kappa_i(x)y_i$，$\kappa_0>0$，且一致椭圆 $\kappa\ge\kappa_{\min}>0$。广义多项式混沌展开 $u\approx\sum_{m=0}^{M}v_m(x,t)\Phi_m(y)$（$M=\binom{N+P}{N}$，$\mathbb E[\Phi_m\Phi_n]=\delta_{mn}$）配 Galerkin 投影给出 $\partial_t v=\nabla\cdot(A\nabla_x v)+f$，其中

$$
a_{jk}=\sum_{i=0}^{N}\kappa_i(x)\,e_{ijk},
\qquad
e_{ijk}=\int y_i\,\Phi_j(y)\Phi_k(y)\rho(y)\,\mathrm dy,
$$

按构造 $A=A^{T}$。

$A(x)$ 是否严格对角占优，决定对角与非对角分裂求解器——混合显隐时间推进、Jacobi 迭代、解耦预条件共轭梯度——是否稳定且以与网格无关的速率收敛。Xiu 与 Shen 已证明 $A$ 逐点正定且每行至多 $2N+1$ 个非零元（因为 $y_i\Phi_j$ 在每个方向上只与最近邻耦合），并对**对称** Beta 密度与 Gamma 密度证明了强化形式的对角占优 $a_{jj}\ge\kappa_{\min}+\sum_{k\ne j}|a_{jk}|$，但把一般非对称 Beta 的情形列为公开问题。

### 归结为三项递推系数的不等式

本文的技术路径是把矩阵问题化为 Jacobi 递推系数的不等式。把方向 $i$ 上的归一化三项递推写成 $xP_j=a^i_jP_{j+1}+b^i_jP_j+c^i_jP_{j-1}$，则

$$
a_{jj}=\kappa_0(x)+\sum_{i=1}^{N}\kappa_i(x)\,b^i_j,
\qquad
\sum_{k\neq j}|a_{jk}|=\sum_{i=1}^{N}|\kappa_i(x)|\,|a^i_j+c^i_j| .
$$

由于各 $y_i$ 同分布，上标 $i$ 可以去掉。**这一步是全篇的技术核心**：它把一个关于 $M\times M$ 矩阵的问题换成一个关于两三个递推系数的初等不等式。

### Theorem 1

若各 $y_i$ 同分布于 $(-1,1)$ 上的 Beta 密度 $\rho(y_i)=(1-y_i)^{\alpha}(1+y_i)^{\beta}$（$\alpha,\beta>-1$）且

$$
|\alpha|\ \ge\ \tfrac12,\qquad |\beta|\ \ge\ \tfrac12,
$$

或各 $y_i$ 同分布于 $(0,+\infty)$ 上的 Gamma 密度 $\rho(y_i)=y_i^{\alpha}e^{-y_i}$（$\alpha>-1$），则相应的 $A(x)$ 对所有 $x\in\Omega$ 严格对角占优：

$$
a_{jj}\ \ge\ \kappa_{\min}+\sum_{k\ne j}|a_{jk}|,
\qquad 1\le j\le M .
$$

公开问题因此被**有条件地**肯定回答：允许非对称的 Beta 参数，代价是显式阈值 $|\alpha|,|\beta|\ge1/2$。论文没有断言对所有 $\alpha\ne\beta$ 都成立；阈值能否去掉在该文中未解决。

一处指数约定值得注意：在 $\rho(y_i)=(1-y_i)^\alpha(1+y_i)^\beta$、$\alpha,\beta>-1$ 下，条件 $|\alpha|\ge1/2$ 同时容纳 $\alpha\ge1/2$ 与 $-1<\alpha\le-1/2$；**均匀（Legendre）情形 $\alpha=\beta=0$ 不被 Theorem 1 覆盖**，但被对称情形的已有结论覆盖。

本文没有数值实验：六页内容就是问题复述、递推系数代数与证明。

## 2 与 4：随机波速使解在随机变量方向正则性很低

### 2：配点的收敛性由数据的假设推出

标量双曲方程 $\partial_t u=c(y)\partial_x u$（$x\in(-1,1)$）在随机波速 $c(y)$ 可**变号**时，流入边界随之在区域两端切换：$c(y)<0$ 时在 $x=-1$ 给数据，$c(y)>0$ 时在 $x=1$ 给数据。因此即使数据光滑，解关于随机变量的正则性也很低。

Gottlieb 与 Xiu 对 Galerkin 情形的收敛证明假设展开系数快速渐近衰减——这是对**答案**的假设，不是对**数据**的假设。本文分析配点方案（边界与初值条件在 $y$ 上逐点施加因此平凡），并从数据的假设推出收敛率。

技术装置是一组带权数据假设，把密度除以波速：

$$
\int_\Gamma\!\!\int_D \rho(y)\bigl(\partial_y u_0(x;y)\bigr)^2\mathrm dx\,\mathrm dy<\infty,
$$

$$
\int_0^T\!\!\int_{\Gamma^+}\frac{\rho(y)}{c(y)}\bigl(\partial_y u_R(t;y)\bigr)^2\mathrm dy\,\mathrm dt<\infty,
\qquad
\int_0^T\!\!\int_{\Gamma^-}\frac{\rho(y)}{|c(y)|}\bigl(\partial_y u_L(t;y)\bigr)^2\mathrm dy\,\mathrm dt<\infty,
$$

其中 $\Gamma^{\pm}$ 是 $c$ 取正负的子集。**这个 $1/|c|$ 权正是把波速与可容许边界数据耦合起来的地方**：波速接近零的参数区域要求边界数据在该处更光滑。在 $|c'(y)|\le C$ 加上上述有限性下解具有 $y$ 方向的 $H^1$ 正则性（能量估计加 Gronwall），再加二阶假设得 $H^2$。

### 4：为什么能观察到指数收敛

编号 2 与 5 给出有限阶正则性结果（$H^k$、$BV$）与相应的**代数**收敛率，但数值上观察到的是指数收敛，理论无法解释。编号 4 补上这个缺口：它建立随机变量方向的**解析**正则性——到参数区间复邻域的解析延拓——并把它转成谱（指数）收敛。论文明确说明它有意忽略确定性求解器的误差，只研究随机空间的离散误差。

关键假设是对数据的全部 $y$ 导数施加**几何增长**界：

$$
\max_{y\in\Gamma}|\partial_y^k c(y)|\le\gamma^k,
\qquad
\max|\partial_y^k u_R|\le\delta_R^k,
\qquad
\max_{y\in\Gamma}\|\partial_y^k u_0\|_V^2\le\eta^k,
$$

常数被抬到 $k$ 次幂——正是这个关于 $k$ 的几何结构使导数界可以求和成一个解析延拓。相应的定理给出

$$
\max_{\Gamma}\|\partial_y^k u(\cdot,t,\cdot)\|_V^2
\ \le\ C_k(T)\bigl(\delta_R^k+\delta_L^k+\eta^k\bigr)<+\infty,
$$

以及在 $\Sigma(\Gamma,\tau)=\{z\in\mathbb C:\mathrm{dist}(z,\Gamma)\le\tau\}$（$0<\tau<1/\sqrt{\zeta}$）内的解析延拓。由此对随机配点得到

$$
\min_{v}\|u-v\|_{L^\infty[T,C^0(\Gamma,V)]}
\ \le\ \frac{2}{\varrho-1}\,e^{-p\log\varrho}
\max_{z\in\Sigma(\Gamma,\tau)}\|u(z)\|_{L^\infty(T,V)},
\qquad
\varrho=\tau+\sqrt{1+\tau^2}>1 ,
$$

即关于多项式次数 $p$ 的指数衰减，速率 $\log\varrho$ 由解析带宽 $\tau$ 决定。

> [!warning] 一处重要的限定
> 论文自己指出：与随机椭圆或抛物问题不同，随机**双曲**方程的解关于随机参数**一般不解析**。因此这条复解析的加强只在特殊情形下可用，例如周期边界条件配解析数据。论文还把它给出的随机 Galerkin 代数估计
> $$\|u-u^N_{SG}\|_{L^2(\Gamma,V)}\le C_\Gamma\sqrt{C(T)}\bigl(\sqrt{\zeta}N\bigr)^{-m}$$
> 自评为「相当粗」。引用指数收敛结论时应带上这两条限定。

数值算例取 $u_t=y\,u_x$、$u(x,0;y)=\cos(y)$，精确解 $u=\cos(x-yt)$。它对每个正整数 $m$ 都属于 $H^{(m)}_y(-1,1)$，因此预期指数收敛，实验对 Galerkin 与配点都观察到。论文同时报告误差随 $t$ 增长——多项式混沌与随机配点长时间退化的已知现象。

## 5 与 3：用双正交多项式解耦 Galerkin 系统

编号 5 面向随机双曲问题的 Galerkin 方法，指出三重困难：双曲方程正则性最差，参数解可能只有 $BV$ 或低阶 $H^k$；Galerkin 投影耦合各模态，需要解耦技术，而与配点不同，边界条件无法逐点施加；若随机场是 $N$ 变量的截断 Karhunen-Loève 展开，各向同性张量基的规模关于 $N$ 指数增长。

论文的对策是**双正交（双重正交）多项式基**，其定义条件为一对：

$$
\int_{\Gamma_i}\rho(y_i)\phi_{j,i}\phi_{k,i}\,\mathrm dy_i=\delta_{jk},
\qquad
\int_{\Gamma_i} y_i\,\rho(y_i)\phi_{j,i}\phi_{k,i}\,\mathrm dy_i=C_{k,i}\,\delta_{jk} .
$$

第二条带权 $y_i$ 的正交性把「乘以 $y_i$」这个算子对角化，从而把 Galerkin 系统解耦成互相独立的确定性求解。这正是编号 1 中系数矩阵 $A$ 的结构被「用掉」的方式：与其研究 $A$ 何时对角占优，不如换一组基使它对角。

针对张量基规模问题，论文用一个敏感性估计给出各向异性的多项式阶选择：形如

$$
\int\!\!\int\rho(y)u_{y_i}^2\,\mathrm dx\,\mathrm dy\le C(T)\bigl(\sqrt{\lambda_i}+\sqrt{\mu_i}+\nu_i\bigr)
$$

的界说明对小 $i$（Karhunen-Loève 谱较大的方向）应取更大的阶 $r_i$。

编号 3 把同一装置用于带随机输入的椭圆**界面**问题。带界面的椭圆问题在确定性情形已经困难：标准有限元除非网格贴合界面否则精度下降。加上随机输入后，朴素随机 Galerkin 路线双重昂贵：投影给出耦合的确定性界面问题系统，每个都需要贴合界面的网格。论文一次去掉两重耦合：随机方向用双正交基解耦，物理方向用**浸入有限元**空间，使网格不必与界面对齐。

## 7 与 38：两个应用方向

- **7（带随机输入的延迟微分方程）** 处理延迟结构。延迟使解对参数的依赖带有额外的分段特征：延迟项把 $[t-\tau,t]$ 上的历史引入，而历史本身依赖参数。
- **38（带随机输入的 Maxwell 方程组）** 给出高效随机 Galerkin 方法。Maxwell 系统是一阶双曲型系统，因此继承编号 2、4、5 关于双曲问题正则性的全部困难，而其散度约束又带来额外的结构要求。

> [!note] 覆盖进度
> 编号 1、2、4、5 的定理与关键构造已按原文核对。编号 7、38 只给出问题设定与定位（限于摘要与元数据层面）。编号 3 是唯一无法核实的一篇：Crossref、OpenAlex、zbMATH 与 Semantic Scholar 均无其摘要，出版社页面返回 403，因此其定理与数值结果本站不报告，上文只给出可从 zbMATH 索引关键词与同组独立来源确认的构造思路。

## 覆盖核对

| 内容                                | 论文 | 覆盖状态                                   |
| ----------------------------------- | ---- | ------------------------------------------ |
| Galerkin 系数矩阵与其稀疏性         | 1    | $a_{jk}$、$e_{ijk}$、对称性、$2N+1$ 非零元 |
| 归结为三项递推系数                  | 1    | 两条恒等式与其技术意义                     |
| Theorem 1 与阈值 $\|\alpha\|\ge1/2$ | 1    | 两类密度、结论、限定与指数约定             |
| 变号波速与流入边界切换              | 2    | 问题设定与正则性后果                       |
| 带 $1/\|c\|$ 权的数据假设           | 2    | 三条积分条件与其耦合含义                   |
| 几何增长假设与解析延拓              | 4    | 假设形式、导数界、解析区域                 |
| 配点的指数收敛与其速率              | 4    | 估计式、$\varrho$ 与 $\tau$ 的关系         |
| 双曲问题一般不解析的限定            | 4    | 论文自陈的两条限定                         |
| 双正交多项式的两条定义条件          | 5    | 定义、解耦机制、与编号 1 的关系            |
| 各向异性阶选择的敏感性依据          | 5    | 界的形式与其结论                           |
| 界面问题的双重解耦                  | 3    | 随机方向与物理方向各自的手段               |

## 本页原文

- T. Zhou and T. Tang, [_Note on coefficient matrices from stochastic Galerkin methods for random diffusion equations_](https://doi.org/10.1016/j.jcp.2010.07.016), J. Comput. Phys. 229 (2010), pp. 8225-8230。
- T. Tang and T. Zhou, [_Convergence analysis for stochastic collocation methods to scalar hyperbolic equations with a random wave speed_](https://doi.org/10.4208/cicp.060109.130110a), Commun. Comput. Phys. 8 (2010), pp. 226-248。
- T. Zhou, [_Stochastic Galerkin methods for elliptic interface problems with random input_](https://doi.org/10.1016/j.cam.2011.05.033), J. Comput. Appl. Math. 236 (2011), pp. 782-792。
- T. Zhou and T. Tang, [_Convergence analysis for spectral approximation to a scalar transport equation with a random wave speed_](https://doi.org/10.4208/jcm.1206-m4012), J. Comput. Math. 30 (2012), pp. 643-656。
- T. Zhou and T. Tang, [_Galerkin methods for stochastic hyperbolic problems using bi-orthogonal polynomials_](https://doi.org/10.1007/s10915-011-9508-0), J. Sci. Comput. 51 (2012), pp. 274-292。
- T. Zhou, [_A stochastic collocation method for delay differential equations with random input_](https://doi.org/10.4208/aamm.2012.m38), Adv. Appl. Math. Mech. 6 (2014), pp. 403-418。
- Z. Feng, J. Li, T. Tang, and T. Zhou, [_Efficient stochastic Galerkin methods for Maxwell's equations with random inputs_](https://doi.org/10.1007/s10915-019-00936-z), J. Sci. Comput. 80 (2019), pp. 248-267。
