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

## 先说机制：随机波速为什么让精度随时间流失

本页大半篇幅围绕同一个模型问题：标量输运方程 $\partial_t u=c(y)\partial_x u$，波速 $c$ 依赖随机参数 $y$。在相容的情形下解沿特征线平移，

$$
u(x,t;y)=u_0\bigl(x+c(y)t\bigr),
$$

这正是编号 5 在其式 (3.14) 中用来估计参数敏感度的那个显式解。对它关于 $y$ 求导，链式法则每次都吐出一个因子 $t$：

$$
\partial_y u=t\,c'(y)\,u_0'\bigl(x+c(y)t\bigr),
\qquad
\partial_y^k u=\bigl(t\,c'(y)\bigr)^k u_0^{(k)}\bigl(x+c(y)t\bigr)+\cdots,
$$

省略项带更低次的 $t$ 与 $c$ 的高阶导数。这个式子就是全页的物理直觉：**两个参数值 $y_1\ne y_2$ 对应的特征线以速度差 $c(y_1)-c(y_2)$ 分离，到时刻 $t$ 已相距 $|c(y_1)-c(y_2)|t$**。因此在固定的 $(x,t)$ 上看，映射 $y\mapsto u(x,t;y)$ 的振荡尺度约为 $1/t$：积分得越久，多项式在随机变量方向要分辨的结构越细。一个 $p$ 次多项式在 $p\gtrsim t$ 之前分辨不了尺度 $1/t$ 的结构。

这条机制解释了本页三处反复出现的现象：

1. 下文每一条定理的常数都带 $T$（$C(T)$、$C_k(T)$、$C_\Sigma C(T)$），**没有一条对时间一致**——收敛性分析必须是含时的；
2. 编号 4 与编号 5 的数值实验都报告固定离散下误差随 $t$ 增长，编号 5 把这个依赖描述为通常呈线性；
3. 同一套机器用在双曲问题上比用在椭圆或抛物问题上差，因为后者没有一个时间变量去逐次放大参数方向的导数。编号 4 给出这一点的尖锐版本：随机**双曲**方程的解关于随机参数**一般不解析**。

> [!note] 这一段的性质
> 上面的求导是对这几篇共用的模型问题做的初等演算，不是某一篇的定理。论文把时间增长作为观察结果报告，并把它编码进 $T$ 依赖的常数里；「参数尺度约为 $1/t$」是本站为串起这些结果给出的解释。编号 4 的数值例子上可以把它算成定量形式，见下文。

## 1：随机 Galerkin 系数矩阵何时严格对角占优

### 直觉

广义多项式混沌 Galerkin 投影把一个随机偏微分方程换成 $M+1$ 个确定性偏微分方程组成的**耦合**系统，全部耦合都装在一个矩阵 $A(x)$ 里。要便宜地推进这个系统，通常的做法是把 $A$ 劈成对角与非对角两块：对角部分隐式、非对角部分显式的混合时间推进，Jacobi 迭代，或者以对角块为预条件的共轭梯度。这类分裂的收敛速度由 $D^{-1}(A-D)$ 的大小控制，而**严格对角占优正是让这个压缩因子小于一的条件**。

关键在于占优是逐点成立的：结论对每个 $x\in\Omega$ 都成立，于是空间离散之后仍然成立，分裂求解器的收敛率与网格无关。而且论文追求的是加强形式

$$
a_{jj}\ \ge\ \kappa_{\min}+\sum_{k\ne j}|a_{jk}|,
$$

余量是椭圆常数 $\kappa_{\min}$ 而不仅仅是「大于零」——这个不随 $j$ 与 $x$ 退化的间隙才是网格无关性的来源。若没有它，分裂迭代可能发散，或者收敛率随网格加密、随多项式阶 $P$ 增大而变差。

Xiu 与 Shen 只对**对称** Beta 密度（$\alpha=\beta$）与 Gamma 密度证明了这一点，并明确把一般非对称 Beta（$\alpha\ne\beta$）列为公开问题。这篇六页短文的全部内容就是回答这个问题：没有新的离散格式，也没有新的求解器。

### 问题设定

模型问题为

$$
\partial_t u(x,y,t)=\nabla\cdot\bigl(\kappa(x,y)\nabla_x u(x,y,t)\bigr)+f(x,y,t),
\qquad x\in\Omega\subset\mathbb R^d\ (d=1,2,3),\ t\in(0,T],
$$

配 $u(x,y,0)=u_0(x,y)$ 与 $u(\cdot,y,t)|_{\partial\Omega}=0$；定常对应物是 $\nabla\cdot(\kappa\nabla_x u)=f$。随机场对参数**仿射**依赖，

$$
\kappa(x,y)=\kappa_0(x)+\sum_{i=1}^{N}\kappa_i(x)\,y_i,
\qquad \kappa_0(x)>0,
\qquad \kappa(x,y)\ge\kappa_{\min}>0,
$$

$y=(y_1,\dots,y_N)$ 的各分量独立同分布。广义多项式混沌展开取

$$
u\approx\sum_{m=0}^{M}v_m(x,t)\Phi_m(y),
\qquad
f\approx\sum_{m=0}^{M}f_m(x,t)\Phi_m(y),
\qquad
M=\binom{N+P}{N},
$$

其中 $\Phi_m(y)=\phi_{m_1}(y_1)\cdots\phi_{m_N}(y_N)$，$m_1+\cdots+m_N\le P$，单变量因子正交归一 $\int\phi_j(y_i)\phi_k(y_i)\rho_i(y_i)\,\mathrm dy_i=\delta_{jk}$，故 $\mathbb E[\Phi_m\Phi_n]=\delta_{mn}$，$\rho(y)=\prod_{i=1}^{N}\rho_i(y_i)$。

### 推导

把展开代入方程、乘以 $\Phi_k$ 取期望，得到分量形式的耦合系统

$$
\partial_t v_k=\sum_{j=1}^{M}\nabla\cdot\bigl(a_{jk}(x)\nabla v_j\bigr)+f_k(x,t),
\qquad k=1,\dots,M,
$$

即矩阵形式 $\partial_t v=\nabla\cdot(A\nabla_x v)+f$，$A(x)=(a_{jk})_{1\le j,k\le M}$，

$$
a_{jk}=\sum_{i=0}^{N}\kappa_i(x)\,e_{ijk},
\qquad
e_{ijk}=\int y_i\,\Phi_j(y)\Phi_k(y)\rho(y)\,\mathrm dy .
$$

按构造 $A=A^{T}$。约定 $y_0\equiv1$，则 $e_{0jk}=\int\Phi_j\Phi_k\rho\,\mathrm dy=\delta_{jk}$，于是求和中 $i=0$ 的那一项只往对角上加 $\kappa_0(x)$；所有真正的耦合都来自 $i\ge1$ 的项，也就是**乘以 $y_i$** 这个算子。

稀疏性由三项递推读出。把方向 $i$ 上的归一化递推写成

$$
y_i\,\phi_{m}(y_i)=a^i_{m}\,\phi_{m+1}(y_i)+b^i_{m}\,\phi_{m}(y_i)+c^i_{m}\,\phi_{m-1}(y_i),
$$

则对多重指标 $j=(j_1,\dots,j_N)$，

$$
y_i\Phi_j=\bigl(a^i_{j_i}\phi_{j_i+1}+b^i_{j_i}\phi_{j_i}+c^i_{j_i}\phi_{j_i-1}\bigr)\prod_{l\ne i}\phi_{j_l},
$$

因此 $e_{ijk}$ 只在 $k$ 与 $j$ 在除 $i$ 外的所有方向都相同、且 $k_i\in\{j_i-1,j_i,j_i+1\}$ 时非零。对 $i=1,\dots,N$ 求并，每行至多 $2N+1$ 个非零元（对角只算一次）——这就是 Xiu–Shen 稀疏性结论的来源。

同一展开给出全篇的技术核心，两条恒等式：

$$
a_{jj}=\kappa_0(x)+\sum_{i=1}^{N}\kappa_i(x)\,b^i_j,
\qquad
\sum_{k\neq j}|a_{jk}|=\sum_{i=1}^{N}|\kappa_i(x)|\,|a^i_j+c^i_j| .
$$

第一条来自 $e_{ijj}=b^i_{j_i}$；第二条把方向 $i$ 上两个近邻的贡献 $\kappa_i a^i_j$ 与 $\kappa_i c^i_j$ 合并，合成单个绝对值 $|a^i_j+c^i_j|$ 用到了归一化递推系数的符号结构。由于各 $y_i$ 同分布，上标 $i$ 可以去掉。**这一步是全篇的技术核心**：它把一个关于 $M\times M$ 矩阵的问题换成一个关于两三个递推系数的初等不等式。

余下的就是 Jacobi 多项式的公式。记 $\varpi:=2n+\alpha+\beta$，正交性给出

$$
\int_{-1}^{1}(1-x)^{\alpha}(1+x)^{\beta}P^{\alpha,\beta}_mP^{\alpha,\beta}_n\,\mathrm dx
=\frac{2^{\alpha+\beta+1}}{2n+\alpha+\beta+1}\frac{\Gamma(n+\alpha+1)\Gamma(n+\beta+1)}{\Gamma(n+\alpha+\beta+1)\,n!}\,\delta_{nm}
=:h^{\alpha,\beta}_n\delta_{nm},
$$

未归一化的递推系数为

$$
a_n=\frac{2(n+1)(n+\alpha+\beta+1)}{(\varpi+1)(\varpi+2)},
\qquad
b_n=\frac{\beta^2-\alpha^2}{\varpi(\varpi+2)},
\qquad
c_n=\frac{2(n+\alpha)(n+\beta)}{\varpi(\varpi+1)},
$$

归一化取 $\widetilde P^{\alpha,\beta}_n=P^{\alpha,\beta}_n/\sqrt{h^{\alpha,\beta}_n}$，得到 $x\widetilde P_n=\tilde a_n\widetilde P_{n+1}+\tilde b_n\widetilde P_n+\tilde c_n\widetilde P_{n-1}$，其中 $\tilde a_n$ 的分子为 $2\sqrt{(n+1)(n+\alpha+1)(n+\beta+1)(n+\alpha+\beta+1)}$。

> [!warning] 一处未能核实的公式
> $\tilde a_n$ 的**分母**（若干 $\varpi$ 平移量之积）在本站可获得的文本中被抽取过程破坏，无法确认其确切形式，因此不在此转录。分子已核实。

$b_n$ 的表达式解释了为什么对称情形先被解决：$b_n=0$ 当且仅当 $\alpha^2=\beta^2$。在 $\alpha=\beta$ 时对角元恰好等于 $\kappa_0(x)$，占优只剩下「用 $\kappa_0-\kappa_{\min}$ 压住非对角行和」这一件事；而 $\alpha\ne\beta$ 时 $b_n\ne0$，在 $\kappa_i b_j<0$ 的地方它会把对角元**拉低**，这正是新增的困难。阈值 $|\alpha|,|\beta|\ge1/2$ 就是使该文的估计能够闭合的条件；文中没有主张这个阈值是最优的。

### 定理

**引理 1（引自 Xiu–Shen）**：$A(x)$ 对每个 $x$ 正定，且每行至多 $2N+1$ 个非零元。

**引理 2（引自 Xiu–Shen）**：若各 $y_i$ 同分布于 $(-1,1)$ 上的**对称** Beta 密度 $\rho(y_i)=(1-y_i)^{\alpha}(1+y_i)^{\alpha}$，或 $(0,+\infty)$ 上的 Gamma 密度 $\rho(y_i)=y_i^{\alpha}e^{-y_i}$（$\alpha>-1$，尺度常数略去），则 $A(x)$ 对所有 $x\in\Omega$ 严格对角占优，且为加强形式 $a_{jj}\ge\kappa_{\min}+\sum_{k\ne j}|a_{jk}|$。

**定理 1（本文的新结论）**：设（I）各 $y_i$ 同分布于 $(-1,1)$ 上的 Beta 密度 $\rho(y_i)=(1-y_i)^{\alpha}(1+y_i)^{\beta}$，$\alpha,\beta>-1$ 且

$$
|\alpha|\ \ge\ \tfrac12,\qquad |\beta|\ \ge\ \tfrac12,
$$

或（II）各 $y_i$ 同分布于 $(0,+\infty)$ 上的 Gamma 密度 $\rho(y_i)=y_i^{\alpha}e^{-y_i}$，$\alpha>-1$（尺度常数略去）。则由相应广义多项式混沌基构造的 $A(x)$ 对所有 $x\in\Omega$ 严格对角占优：

$$
a_{jj}\ \ge\ \kappa_{\min}+\sum_{k\ne j}|a_{jk}|,
\qquad 1\le j\le M,\ \forall x\in\Omega .
$$

公开问题因此被**有条件地**肯定回答：允许非对称的 Beta 参数，代价是显式阈值 $|\alpha|,|\beta|\ge1/2$。论文没有断言对所有 $\alpha\ne\beta$（$\alpha,\beta>-1$）都成立；阈值能否去掉在该文中未解决。

一处指数约定值得注意：在 $\rho(y_i)=(1-y_i)^\alpha(1+y_i)^\beta$、$\alpha,\beta>-1$ 下，条件 $|\alpha|\ge1/2$ 同时容纳 $\alpha\ge1/2$ 与 $-1<\alpha\le-1/2$；**均匀（Legendre）情形 $\alpha=\beta=0$ 不被定理 1 覆盖**，但被引理 2 的对称情形覆盖。

### 数值实验

没有。这是一篇纯分析性的短文：六页内容是问题复述、递推系数代数与证明，全文不含数值算例。

### 与其他论文的关系

这是本组里最早、也最纯代数的一篇——它研究的是随机 Galerkin 投影**产生的线性代数**，而不是采样或逼近论。它与编号 3、5、38 共享侵入式 Galerkin 的主题，与分析非侵入式配点的编号 2、4 配对。与后来的采样设计路线（[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|最小二乘一页]]与[[computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce|稀疏恢复一页]]）基本不相交：那条路线放弃 Galerkin 投影，改从点值做最小二乘与 $\ell_1$ 恢复。

## 2：随机波速下配点法的收敛率由数据假设推出

### 直觉

标量双曲方程 $\partial_t u=c(y)\partial_x u$ 在随机波速 $c(y)$ 可**变号**时，流入边界随之在区域两端切换：$c(y)<0$ 时在 $x=-1$ 给数据，$c(y)>0$ 时在 $x=1$ 给数据。在使 $c$ 过零的那个参数值附近，问题的边界条件类型发生改变，于是**即使数据在 $x$ 与 $y$ 上都光滑，映射 $y\mapsto u(x,t;y)$ 也可能只有很低的正则性**。

这一点对两族方法的影响不同。配点法在每个节点 $y_j$ 上 $c(y_j)$ 有确定的符号，于是照常施加那个标量问题的适定边界条件即可，边界与初值条件在 $y$ 上逐点施加因而是平凡的。Galerkin 投影得到的却是一个特征值有正有负的对称双曲系统，边界条件必须在特征变量上一致地施加——这正是 Gottlieb 与 Xiu 处理的困难。

第二个动机在假设的位置上。Gottlieb 与 Xiu 对 Galerkin 情形的收敛证明假设展开系数快速渐近衰减——这是对**答案**的假设，不是对**数据**的假设。本文要的是从 $c$、$u_0$、$u_L$、$u_R$ 出发推出收敛率。

### 问题设定

$$
\partial_t u(x,t;y)=c(y)\,\partial_x u(x,t;y),
\qquad x\in D\equiv(-1,1),\ t>0,
$$

配 $u(x,0;y)=u_0(x;y)$ 与适定的边界条件：$c(y)<0$ 时 $u(-1,t;y)=u_L(t;y)$，$c(y)>0$ 时 $u(1,t;y)=u_R(t;y)$。随机变量 $y$ 在 $\Gamma\equiv[-1,1]$ 上有密度 $\rho(y)$。

配点格式取 Gauss 节点 $\{y_i\}_{i=0}^{N}$（$\Phi_{N+1}$ 的零点；Hermite 对应 Gauss 分布、Legendre 对应均匀分布、Laguerre 对应 Gamma 分布），求解 $N+1$ 个**互不耦合**的确定性问题，再插值：

$$
u^N(x,t;y)=\mathcal I^y_N u:=\sum_{k=0}^{N}u(x,t;y_k)\,F_k(y),
\qquad F_k\in\mathbb P_N,\ F_i(y_k)=\delta_{ik}.
$$

误差用两个泛函衡量：均方误差 $e_{ms}(u-u^N):=M[u-u^N]$ 与均值误差 $e_{mean}(u-u^N):=\mathbb E[|u-u^N|]$。随机方向的正则性在 $H^1$、$H^2$、$BV$ 中度量，用到张量空间同构 $L^2\otimes H^k(D)\simeq L^2(\Gamma;H^k(D))\simeq H^k(D;L^2(\Gamma))$。

技术装置是一组带权数据假设（式 (2.8)），把密度除以波速：

$$
\int_\Gamma\!\!\int_D \rho(y)\bigl(\partial_y u_0(x;y)\bigr)^2\mathrm dx\,\mathrm dy<\infty,
$$

$$
\int_0^T\!\!\int_{\Gamma^+}\frac{\rho(y)}{c(y)}\bigl(\partial_y u_R(t;y)\bigr)^2\mathrm dy\,\mathrm dt<\infty,
\qquad
\int_0^T\!\!\int_{\Gamma^-}\frac{\rho(y)}{|c(y)|}\bigl(\partial_y u_L(t;y)\bigr)^2\mathrm dy\,\mathrm dt<\infty,
$$

其中 $\Gamma^{\pm}$ 是 $c$ 取正负的子集。**这个 $1/|c|$ 权正是把波速与可容许边界数据耦合起来的地方**：波速接近零的参数区域要求边界数据在该处更光滑。

### 推导

正则性由对 $y$ 求导后的能量估计得到。令 $w=\partial_y u$，对方程求 $y$ 导数：

$$
\partial_t w=c(y)\,\partial_x w+c'(y)\,\partial_x u .
$$

在 $D$ 上乘 $w$ 积分，输运项分部积分成边界通量：

$$
\frac{\mathrm d}{\mathrm dt}\frac12\int_D w^2\,\mathrm dx
=\frac{c(y)}{2}\Bigl[w^2\Bigr]_{x=-1}^{x=1}
+\int_D c'(y)\,\partial_x u\,w\,\mathrm dx .
$$

$c(y)>0$ 时流入端是 $x=1$，那里 $w(1,t;y)=\partial_y u_R(t;y)$ 由数据给定，而流出端以有利的符号进入；$c(y)<0$ 时两端互换。再乘 $\rho(y)$ 在 $\Gamma$ 与 $[0,T]$ 上积分，边界项就落到 (2.8b–c) 要控制的那两个量上；交叉项 $\int c'\,\partial_x u\,w$ 用 Cauchy–Schwarz 处理，需要 $|c'(y)|\le C$ 与 $\partial_x u$ 的先验界（定理的假设中除 (2.8a–c) 外还引用了更早的有限性条件 (2.1)），最后由 Gronwall 闭合，$C(T)$ 即由此而来。二阶导数重复同一过程。

> [!note] 关于 $1/|c|$ 权的来源
> (2.8b–c) 的权是 $\rho/|c|$，而上面这条能量恒等式的边界通量带的是 $c$ 本身，二者相差 $c^2$。一个自然的对账是：在流入端由方程本身有 $\partial_x u=\partial_t u/c$，把边界上的时间导数换成空间导数正好付出一个 $1/c$；编号 4 的假设里直接出现 $d_R:=\partial_t u_R/c(y)$、$d_L:=\partial_t u_L/c(y)$，与这一换算一致。本站没有拿到该处的完整推导，这里只指出权的形式与这一换算相符，(2.8) 本身是按原文转录的。

收敛率来自一条标准的插值估计与上述正则性的组合：**引理 3.1**（引自 Canuto 等，第 289 页）说，对 $w^{(m)}\in L^2(-1,1)$ 与其在 $N+1$ 个 Gauss / Gauss–Radau / Gauss–Lobatto 点上的插值 $\mathcal I_N w$，

$$
\|w-\mathcal I_N w\|_{L^2(D)}\le C\,N^{-m}\,\|w^{(m)}\|_{L^2(-1,1)},
\qquad m\le N .
$$

把 $m$ 取成随机方向的正则性指标，就得到定理 3.1 与 3.2：**收敛阶等于正则性阶**。

$BV$ 情形走的是另一条路。此时没有可用的插值估计，改用 Gauss 求积对有界变差函数的余项界（**引理 3.2**，引自文献）：若 $f$ 在 $[-1,1]$ 上的全变差为 $V(f)$，则 $N$ 点 Gauss 求积余项满足

$$
|R_N(f)|\le\frac{\pi}{2N+1}V(f).
$$

均值误差 $\mathbb E[|u-u^N|]$ 本身就是一个对 $\rho$ 的积分，可以直接用这条余项界；而均方误差需要的是逐点控制，$BV$ 给不出来。**这就是为什么定理 3.3 只控制均值、不控制均方**。

### 定理

**定理 2.1（$H^1$ 正则性）**：设 $|c'(y)|\le C$ 在 $\Gamma$ 上几乎处处成立（分布意义下有界），且有限性假设 (2.1) 与 (2.8a–c) 成立，则解在 $y$ 方向有有界的 $H^1$ 正则性，对 $0<t\le T$ 一致。

**定理 2.2（$H^2$ 正则性）**：若另设 $c''(y)$ 在分布意义下有界，且相应的二阶导数数据假设成立，则解在 $y$ 方向有 $H^2$ 正则性。

**定理 2.3（$BV$ 正则性）**：设 $|c'(y)|\le C$，且数据满足有界变差型假设（式 (2.30a) 及其后，例如 $\int_\Gamma\int_D\rho(y)|\partial_x u_0(x;y)|\,\mathrm dx\,\mathrm dy<+\infty$），则解在 $y$ 方向有界变差。

**定理 3.1**：在定理 2.1 的假设下（$y$ 方向 $H^1$），

$$
e_{ms}(u-u^N)\le C(T)\,N^{-1},
\qquad
e_{mean}(u-u^N)\le C(T)\,N^{-1},
\qquad 0<t\le T,
$$

$C(T)$ 依赖 $T$ 但与 $N$ 无关。

**定理 3.2**：在定理 2.2 的假设下（$y$ 方向 $H^2$），$e_{ms}\le C(T)N^{-2}$ 且 $e_{mean}\le C(T)N^{-2}$，$0<t\le T$。

**定理 3.3**：在定理 2.3 的假设下（$y$ 方向 $BV$），**只有均值误差**被控制：$e_{mean}=\mathbb E[|u-u^N|]\le C(T)N^{-1}$，$0<t\le T$。

总的判断是：**收敛率由解在 $y$ 方向的正则性支配，而这个正则性由随机波速与初边值数据共同决定，不是单靠数据光滑就能得到的。**

### 数值实验

实验取 $y$ 在 $\Gamma$ 上均匀分布，用 Legendre–Gauss 配点。4.1 节求解 $u_t=y\,u_x$，取三组初值 $u(x,0;y)=\sin(x)+4\,\mathrm{sgn}(y)\,y^{k}$，随 $k$ 增大抬高 $y$ 方向的正则性：

| 初值中的 $4\,\mathrm{sgn}(y)y^{k}$ | 该项的形状          | $y$ 方向正则性 | 依据                | 预测阶   |
| ---------------------------------- | ------------------- | -------------- | ------------------- | -------- |
| $k=1$                              | $\lvert y\rvert$    | $H^1$          | 定理 2.1 + 定理 3.1 | $N^{-1}$ |
| $k=2$                              | $y\lvert y\rvert$   | $H^2$          | 定理 2.2 + 定理 3.2 | $N^{-2}$ |
| $k=3$                              | $y^2\lvert y\rvert$ | $H^3$          | 引理 3.1 取 $m=3$   | $N^{-3}$ |

论文报告观察到的阶随正则性指标同步上升。（注意论文列出的定理只到 $H^2$；$H^3$ 这一行的预测阶来自引理 3.1 对一般 $m$ 的形式，不是文中单列的定理。本站未获得逐点误差数值，只能报告论文陈述的阶。）

另一个算例取满足 $BV$ 假设 (2.32a) 但**不**满足 $H^1$ 假设 (2.8a) 的数据：

| 误差       | 理论预测                      | 论文观察               |
| ---------- | ----------------------------- | ---------------------- |
| $e_{mean}$ | $O(N^{-1})$（定理 3.3）       | 一阶，与定理相符       |
| $e_{ms}$   | 定理 3.1 的假设不满足，无预测 | **不是一阶**，约为半阶 |

这一组才是实验的要点所在：它说明**定理 3.1 的假设是尖锐的**——去掉 $H^1$ 假设之后，定理 3.3 保住的均值误差仍是一阶，而均方误差立刻掉下来。论文另外考察了一种精度增强技术，沿用 Foo、Wan 与 Karniadakis 的多元素配点方法（_J. Comput. Phys._ 227 (2008) 9572–9595）。

### 与其他论文的关系

这是与编号 4 配对的配点一半：编号 4 对同一输运模型在谱／Galerkin 框架下重做分析。编号 5 则用双正交多项式构造 Galerkin 格式，绕开这里让 Galerkin 变得别扭的耦合。它与编号 1 共享 Gottlieb–Xiu / Xiu–Shen 的侵入式方法背景。方法论上它属于「结构化节点、一维参数」的阶段，后面的[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|最小二乘一页]]与[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|最优采样一页]]转向随机的或经过设计的多元样本。

## 4：解析正则性与谱收敛，以及它为什么随时间变差

### 直觉

编号 2 与 5 给出有限阶正则性结果（$H^k$、$BV$）与相应的**代数**收敛率，但数值上观察到的是指数收敛，理论无法解释。缺口在哪里很清楚：指数收敛要的不是「导数很多」，而是「导数的增长至多是几何式的」，这恰好就是解析性。于是这篇论文的路线是：对数据假设导数的几何增长，把它沿方程传递到解上，把 Taylor 级数求和成一个到复邻域的解析延拓，再套用经典的最佳逼近估计。论文明确说明它有意忽略确定性求解器的误差，只研究随机空间的离散误差。

### 问题设定

模型与编号 2 相同，取 $V=L^2(D)$。**假设 2.1**（式 (2.7a)–(2.7d)）对数据的全部 $y$ 导数施加几何增长界。记 $d_R:=\partial_t u_R/c(y)$，$d_L:=\partial_t u_L/c(y)$，对每个 $k\in\mathbb N$：

$$
\max_{y\in\Gamma}|\partial_y^k c(y)|\le\gamma^k;
\qquad
\max_{\Gamma\otimes T}|\partial_y^k d_R|\le\delta_R^k,\quad
\max_{\Gamma\otimes T}|\partial_y^k u_R|\le\delta_R^k;
$$

$$
\max_{\Gamma\otimes T}|\partial_y^k d_L|\le\delta_L^k,\quad
\max_{\Gamma\otimes T}|\partial_y^k u_L|\le\delta_L^k;
\qquad
\max_{y\in\Gamma}\|\partial_y^k u_0'\|_V^2\le\eta^k,
\quad
\max_{y\in\Gamma}\|\partial_y^k u_0\|_V^2\le\eta^k,
$$

其中 $u'=\partial_x u$，不妨设 $\gamma\ge\max\{\delta_R,\delta_L,\eta\}$。**常数被抬到 $k$ 次幂**——正是这个关于 $k$ 的几何结构使导数界可以求和成一个解析延拓。

量 $d_R$ 值得单独说一句：在流入端 $x=1$ 上由方程本身有 $\partial_x u=\partial_t u/c$，因此 $d_R=\partial_t u_R/c$ 就是解在那个边界上的空间导数。这也解释了为什么 $c$ 出现在分母上，并与编号 2 中 $1/|c|$ 权的位置一致。

### 推导

对方程求 $k$ 次 $y$ 导数，Leibniz 法则给出

$$
\partial_t\bigl(\partial_y^k u\bigr)
=c(y)\,\partial_x\partial_y^k u
+\sum_{l=1}^{k}\binom{k}{l}\bigl(\partial_y^l c\bigr)\,\partial_x\partial_y^{k-l}u .
$$

首项是原来的输运算子，作用在 $\partial_y^k u$ 上；求和里的每一项都含更低阶的 $y$ 导数，可用归纳假设控制，而 $|\partial_y^l c|\le\gamma^l$ 给出系数的界。像编号 2 那样做能量估计：边界通量由 $\partial_y^k u_R$、$\partial_y^k u_L$ 给出，分别被 $\delta_R^k$、$\delta_L^k$ 压住，初值项被 $\eta^k$ 压住，Gronwall 在时间上闭合并产生依赖 $k$、$\gamma$ 与 $T$ 的常数 $C_k(T)$。论文的证明就是对 $k$ 的归纳。

有了这族导数界，解析延拓由幂级数

$$
u(z,x,t)=\sum_{k=0}^{\infty}(z-y)^k\,\partial_y^k u(y,x,t)
$$

给出，收敛区域为

$$
\Sigma(\Gamma,\tau)\equiv\{z\in\mathbb C:\ \mathrm{dist}(z,\Gamma)\le\tau\},
\qquad 0<\tau<1/\sqrt{\zeta}.
$$

两件事在这里对上了。定理 2.1 控制的是导数的**平方** $\|\partial_y^k u\|_V^2\lesssim\zeta^{k}$，于是 $\|\partial_y^k u\|_V\lesssim\zeta^{k/2}$，级数的公比是 $\sqrt{\zeta}\,|z-y|$，收敛半径正好是 $1/\sqrt{\zeta}$——**这就是定理 2.2 中 $\tau<1/\sqrt\zeta$ 的来历**。这也说明上面那个级数写成没有 $1/k!$ 的形式并非笔误：若带 $1/k!$，几何增长的导数界会给出处处收敛，$1/\sqrt\zeta$ 这个有限半径就无从谈起。

> [!note] 两处限定
> 级数按原文转录，其中不出现 $1/k!$，与通常的 Taylor 写法不同。常数 $\zeta$ 在论文前文引入，本站无法核实其确切定义，只能确认它是由假设 2.1 决定的、依赖数据的常数。

### 定理

**定理 2.1**：在假设 2.1 下，对每个 $k\in\mathbb N$，

$$
\max_{\Gamma}\|\partial_y^k u(\cdot,t,\cdot)\|_V^2
\ \le\ C_k(T)\bigl(\delta_R^k+\delta_L^k+\eta^k\bigr)\ <\ +\infty,
$$

$C_k(T)$ 依赖 $k$、$\gamma$ 与 $T$（对 $k$ 归纳证明）。论文由此指出 $u\in L^\infty[T,C^0(\Gamma,V)]$。

**定理 2.2**：解作为 $y$ 的函数 $u:\Gamma\to L^\infty(T,V)$ 在 $\Sigma(\Gamma,\tau)$（$0<\tau<1/\sqrt\zeta$）内存在解析延拓 $u(x,t;z)$，$z\in\mathbb C$。

**定理 3.1（随机配点；论文把证明归于文献并略去）**：若 $u\in L^\infty[T,C^0(\Gamma,V)]$ 且对某个 $\tau>0$ 在 $\Sigma(\Gamma,\tau)$ 内解析，则

$$
\min_{v\in L^1[T,\mathbb P_p(\Gamma)\otimes V]}\|u-v\|_{L^\infty[T,C^0(\Gamma,V)]}
\ \le\ \frac{2}{\varrho-1}\,e^{-p\log\varrho}
\max_{z\in\Sigma(\Gamma,\tau)}\|u(z)\|_{L^\infty(T,V)},
\qquad
\varrho=\tau+\sqrt{1+\tau^2}>1 ,
$$

即关于多项式次数 $p$ 的指数衰减，速率 $\log\varrho$ 由解析带宽 $\tau$ 决定。

**引理 3.3（投影误差，引自谱方法文献）**：对 $u\in H^m(I)$、$m\ge1$ 与充分大的 $N$，$\|u-\Pi_N u\|_{L^2(I)}\le C_I N^{-m}|u|_{H^m(I)}$，其中 $|u|_{H^m(I)}=\|\partial^m u\|_{L^2(I)}$，$C_I$ 与 $N$ 无关。

**定理 3.2（随机 Galerkin，代数形式）**：

$$
\|u-u^N_{SG}(\cdot,t,\cdot)\|_{L^2(\Gamma,V)}
\ \le\ C_\Gamma\sqrt{C(T)}\,\bigl(\sqrt{\zeta}\,N\bigr)^{-m},
$$

$m$ 是与解在随机空间中正则性相联系的整数指标。

**引理 3.4（引自编号 2）**：对任意有限时刻 $t$，广义多项式混沌误差满足 $\mathbb E\|u-u^p\|_2^2\le C(T)\sum_{k=p+1}^{\infty}\|\tilde u_k\|_1^2$，其中 $u^p=\sum_{k=0}^{p}\tilde u_kP_k(y)$，$\|u\|_1^2=\int_D(u^2+u_x^2)\,\mathrm dx$。

**引理 3.5（周期特例，$c(z)$ 与 $u_0$ 解析）**：

$$
\|\tilde u_n\|_1\ \le\ \frac{C_\Sigma}{2^n}\sqrt{\frac{2n+1}{2}}\int_{-1}^{1}\left(\frac{1-y^2}{1-|y|+\tau}\right)^{n}\mathrm dy .
$$

**定理 3.3（广义多项式混沌的谱收敛，按原文所印）**：

$$
\bigl(\mathbb E\|u-u^p\|_2^2\bigr)^{1/2}
\ \le\ \sqrt{C_\Sigma C(T)}\left(\frac{\sqrt{\pi}}{\sqrt{1-r^2}}+O\!\left(\frac{1}{p^{1/3}}\right)\right)\frac{r^{p+1}}{\sqrt{1-r^2}},
\qquad \xi=-1-\tau<-1 ,
$$

几何因子 $r^{p+1}$（$0<r<1$）就是谱收敛率。$r$ 的闭式在本站可获得的文本中被抽取过程破坏，只能确认形如 $r\equiv1/(|\xi|+\cdots)$ 且 $\xi=-1-\tau$；其余部分**未核实**，但它的角色是清楚的：一个由解析带宽 $\tau$ 决定的 Bernstein 椭圆型参数。

> [!warning] 一处重要的限定
> 论文自己指出：与随机椭圆或抛物问题不同，随机**双曲**方程的解关于随机参数**一般不解析**。因此这条复解析的加强只在特殊情形下可用，例如周期边界条件配解析数据。论文还把它给出的随机 Galerkin 代数估计
> $$\|u-u^N_{SG}\|_{L^2(\Gamma,V)}\le C_\Gamma\sqrt{C(T)}\bigl(\sqrt{\zeta}N\bigr)^{-m}$$
> 自评为「相当粗」。引用指数收敛结论时应带上这两条限定。

### 数值实验

第 4 节只有一个算例，沿用编号 2 的设定：

| 项目     | 设定                                                |
| -------- | --------------------------------------------------- |
| 方程     | $u_t(x,t;y)=y\,u_x(x,t;y)$，$0<x<2\pi$，$t>0$       |
| 初值     | $u(x,0;y)=\cos(y)$（按原文转录，见下面的注）        |
| 边界     | 取成使精确解为 $u(x,t;y)=\cos(x-yt)$                |
| 扫描量   | 投影项数（随机 Galerkin）／配点节点数（随机配点）   |
| 记录     | 均方误差对扫描量的曲线（图 4.1、4.2），若干个时间层 |
| 观察结果 | 两族方法都观察到指数收敛；误差随 $t$ 增长           |

这个解对每个正整数 $m$ 都属于 $H^{(m)}_y(-1,1)$，因此预期指数收敛，实验对 Galerkin 与配点都观察到。误差随 $t$ 增长是多项式混沌与随机配点长时间退化的已知现象。

> [!warning] 算例数据的一处不自洽
> 上表转录的三条信息互不相容：$u(x,0;y)=\cos(y)$ 与 $u=\cos(x-yt)$ 在 $t=0$ 处对不上（后者给出 $\cos(x)$）；而且 $\cos(x-yt)$ 满足的是 $u_t=-y\,u_x$，不是 $u_t=y\,u_x$。合理的修正是取 $u_0=\cos(x)$ 并相应统一符号约定（$u_t=y\,u_x$ 对应 $u=\cos(x+yt)$）。本站按可获得的文本转录并在此标出，不擅改；无论采用哪一种修正，例子的要点——解在 $y$ 方向属于任意阶 $H^{(m)}$、因而预期并观察到指数收敛——都不受影响。

这个例子还能把本页开头的时间机制算成定量形式。对 $u=\cos(x\mp yt)$ 有 $|\partial_y^k u|\le t^k$，于是幂级数 $\sum(z-y)^k\partial_y^k u$ 在 $|z-y|<1/t$ 内收敛：**可用的解析带宽 $\tau$ 随时间像 $1/t$ 一样收缩**。小 $\tau$ 时 $\varrho=\tau+\sqrt{1+\tau^2}\approx1+\tau$，故 $\log\varrho\approx\tau\approx1/t$，定理 3.1 的因子 $e^{-p\log\varrho}$ 变成大致的 $e^{-p/t}$：要把误差保持在同一水平，多项式次数必须随 $t$ 成比例增长。等价的读法是把 $t$ 的影响放进前因子——延拓到宽度 $\tau$ 的带上时 $|\cos(x\mp zt)|$ 可达 $\cosh(\tau t)$ 量级，于是速率 $\log\varrho$ 与时间无关而前因子 $\max_{z\in\Sigma}\|u(z)\|$ 随 $t$ 增长。两种读法都与论文报告的「误差随 $t$ 增长」相符。这段定量化是本站在论文自己的例子上做的演算，论文只给出定性的观察。

### 与其他论文的关系

这是输运方程三部曲 2、4、5 的理论收口：编号 2 给出 $H^1$／$H^2$／$BV$ 正则性与配点的 $N^{-1}$／$N^{-2}$ 代数率，编号 5 给出解耦的 Galerkin 格式，本文给出解析正则性与**两族方法共用**的谱收敛率，论文自己把它描述为那些数值论文的理论补充。它「可达速率由参数正则性决定」的判断，在后面的[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|最小二乘]]与[[computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce|稀疏恢复]]论文里变成另一个问题：需要多少样本，而不是能达到什么阶。

## 5：双正交多项式把 Galerkin 系统精确解耦

### 直觉

编号 1 问的是「系数矩阵 $A$ 什么时候对角占优」。编号 5 问了一个更好的问题：**能不能换一组基，让它直接是对角的**。

这个想法之所以可行，是因为耦合的来源非常单一。波速对参数仿射依赖，所以 Galerkin 系统里唯一的耦合算子就是「乘以 $y_i$」。普通的正交基只让 $L^2_\rho$ 内积对角化（即质量矩阵是单位阵），对「乘以 $y_i$」什么也没做；**双正交基同时对角化这两个二次型**，于是耦合项也变成对角，整个 Galerkin 系统散成一组互不相干的标量输运问题。

值得点明的是：编号 5 要对角化的那个矩阵，元素是 $\int y_i\phi_j\phi_k\rho\,\mathrm dy_i$——与编号 1 里 $e_{ijk}$ 完全是同一个对象，也就是三项递推的 Jacobi 矩阵，其元素正是编号 1 的 $a^i_j$、$b^i_j$、$c^i_j$。**两篇论文处理的是同一个矩阵：编号 1 界定它的非对角部分，编号 5 把它对角化掉。**

不这样做的代价是实打实的：Galerkin 投影耦合各模态；与配点不同，边界条件无法逐点施加（耦合系统是特征值有正有负的对称双曲系统）；而若随机场是 $N$ 变量的截断 Karhunen–Loève 展开，各向同性张量基的规模关于 $N$ 指数增长。加上双曲方程本身正则性最差（参数解可能只有 $BV$ 或低阶 $H^k$），这就是论文正面处理的三重困难。

### 问题设定

$$
\partial_t u(x,t,y)=c(x,y)\,\partial_x u(x,t,y),\qquad x\in D\equiv[-1,1],
$$

配 $u(1,t;y)=u_R(t,y)$（在 $c(1,y)>0$ 处）、$u(-1,t;y)=u_L(t,y)$（在 $c(-1,y)<0$ 处）与 $u(x,0;y)=u_0(x,y)$。关键的推广是：波场 $c$（等价地记作 $\kappa$）同时依赖 $x$ 与**向量** $y=(y_1,\dots,y_N)$，是真正的随机场，而不是编号 2 里的一个随机变量。

（记号上文中随机变量个数在模型处记作 $N$、在基的构造处记作 $M$，指的是同一个计数。）

**双正交（双重正交）基**是全篇的中心对象。对每个方向 $i$，要求 $\mathbb P_{r_i}$ 的单变量基 $\{\phi_{j,i}\}_{j=0}^{r_i}$ 同时满足**两条**正交关系：

$$
\int_{\Gamma_i}\rho(y_i)\phi_{j,i}\phi_{k,i}\,\mathrm dy_i=\delta_{jk},
\qquad
\int_{\Gamma_i} y_i\,\rho(y_i)\phi_{j,i}\phi_{k,i}\,\mathrm dy_i=C_{k,i}\,\delta_{jk},
$$

其中 $\{C_{k,i}\}_{k=0}^{r_i}$ 是非零常数。论文把这个空间归于文献，并指出构造它等价于求解**特征值问题**，其代价在 $r_j$ 不大时相比求解耦合系统可以忽略。$C_{k,i}$ 是「乘以 $y_i$」算子的特征值，等价地就是 Gauss 求积节点——这正是论文注 3.2 所说的、该 Galerkin 方法**等价于某种配点方法**的原因，只是保留了更便于分析的 Galerkin 框架。

多变量基取张量积：$\mathbb P_{\mathbf r}=\mathbb P_{r_1}\otimes\cdots\otimes\mathbb P_{r_M}\subset L^2(\Gamma,\rho)$，多重指标 $\mathbf i=(i_1,\dots,i_M)\le\mathbf r$，

$$
\Phi_{\mathbf i}(y)=\prod_{k=1}^{M}\phi_{i_k,k}(y_k),
\qquad
\int_\Gamma y_k\,\rho(y)\,\Phi_{\mathbf i}(y)\Phi_{\mathbf j}(y)\,\mathrm dy=C_{i_k,k}\,\delta_{\mathbf i\mathbf j},
$$

共 $N_y=\prod_{i=1}^{M}(r_i+1)$ 个基函数。

### 推导

写 $u=\sum_{\mathbf i\le\mathbf r}u_{\mathbf i}\Phi_{\mathbf i}(y)$，代入仿射波场并对 $\Phi_{\mathbf j}$ 投影。$\kappa_0$ 那一项由第一条正交性给出 $\delta_{\mathbf i\mathbf j}$，$\kappa_k y_k$ 那一项由第二条正交性给出 $C_{i_k,k}\delta_{\mathbf i\mathbf j}$——**两项都是对角的**，于是得到 $N_y$ 个彼此独立的确定性方程：

$$
\frac{\partial u_{\mathbf i}}{\partial t}=\kappa_{\mathbf i}(x)\frac{\partial u_{\mathbf i}}{\partial x},
\qquad
u_{\mathbf i}(x,0)=u_{0\mathbf i}(x),
\qquad
\kappa_{\mathbf i}(x)=\kappa_0(x)+\sum_{k=1}^{M}\kappa_k(x)\,C_{i_k}.
$$

耦合的 Galerkin 系统于是坍缩成一族标量输运问题，各自带一个**有效波速** $\kappa_{\mathbf i}(x)$。边界条件可以逐模态施加：$\kappa_{\mathbf i}(1)>0$ 时取 $u_{\mathbf i}(1,t)=u_{R\mathbf i}(t)$，$\kappa_{\mathbf i}(-1)<0$ 时取 $u_{\mathbf i}(-1,t)=u_{L\mathbf i}(t)$。

这一步需要有效波速的符号判断是有意义的，论文给出的理由很干净：在第二条定义式中取 $j=k$ 得

$$
C_{k,i}=\int_{\Gamma_i}y_i\,\rho(y_i)\,\phi_{k,i}^2\,\mathrm dy_i,
$$

而由第一条定义式，$\rho\phi_{k,i}^2$ 是一个积分为 $1$ 的概率密度。**所以 $C_{k,i}$ 是 $y_i$ 在某个概率密度下的平均值，必然落在 $y_i$ 的上下界之间、即落在其支撑的凸包内。**有效波速因此仍在原波速的取值范围内，符号判断与原问题一致。

统计量也是白拿的：

$$
\mathbb E(u^M)=\sum_{\mathbf i\le\mathbf r}u_{\mathbf i}(x,t)\int_\Gamma\rho(y)\Phi_{\mathbf i}\,\mathrm dy,
\qquad
\int_\Gamma\rho(y)\Bigl(\sum_{\mathbf i\le\mathbf r}u_{\mathbf i}\Phi_{\mathbf i}\Bigr)^2\mathrm dy=\sum_{\mathbf i\le\mathbf r}\bigl(u_{\mathbf i}(x,t)\bigr)^2 .
$$

剩下的是张量基规模问题，论文用敏感性估计换成各向异性的阶选择。若 Karhunen–Loève 特征值按 $\lambda_i\sim i^{-2m}$（$m\ge1$，衰减指标由协方差函数决定）衰减，且对系数 $\kappa_i$ 为常数的相容问题精确解 $u(x,y,t)=u_0(x+\kappa(y)t)$ 给出 $|\partial u/\partial y_i|\sim i^{-m}$、$|\partial^2u/\partial y_i\partial y_j|\sim i^{-m}j^{-m}$（式 (3.14)），则主要误差来自低指标方向的一阶导数。**规则因此是：小 $i$ 用较大的 $r_i$，大 $i$ 用较小的 $r_i$，并服从总次数约束 $\sum_{i=1}^{M}r_i\le P$。**这把 $N_y$ 压到远低于各向同性的 $r^M$。

注 3.1 还指出构造可直接推广到多个空间维度，例如 $\partial_t u=\kappa_1(x_1,x_2,y)\partial_{x_1}u+\kappa_2(x_1,x_2,y)\partial_{x_2}u$ 配相互独立的随机波场，只需在每个方向建基后相乘；推广到双曲**方程组**也被称为是直接的。

### 定理

**引理 3.1**：设 $\kappa(x,y)=\kappa_0(x)+\sum_{i=1}^{N}\sqrt{\lambda_i}\,\kappa_i(x)y_i$ 满足 $0<\kappa_{\min}\le\kappa(x,y)<\kappa_{\max}$ 与 $|\partial_x\kappa(x,y)|<\bar\kappa_{\max}$。若

$$
\int_\Gamma\!\!\int_D\rho(y)\bigl(\partial_x u_0(x;y)\bigr)^2\mathrm dx\,\mathrm dy<\infty,
\qquad
\int_0^T\!\!\int_\Gamma\rho(y)\bigl(\partial_t u_R(t;y)\bigr)^2\mathrm dy\,\mathrm dt<\infty
$$

以及 $u_L$ 的同类条件成立，则 $\int_\Gamma\int_D\rho(y)u_x^2\,\mathrm dx\,\mathrm dy\le C(T)<\infty$，$0<t\le T$。

**定理 3.1（敏感性估计，式 (3.17)）**：设 $u_0$ 确定，边界数据为 Karhunen–Loève 形式

$$
u_R(y,t)=u_{R0}(t)+\sum_{i=1}^{N}\sqrt{\mu_i}\,u_{Ri}(t)y_i,
\qquad
u_L(y,t)=u_{L0}(t)+\sum_{i=1}^{N}\sqrt{\nu_i}\,u_{Li}(t)y_i,
$$

其中 $\int_0^T u_{Ri}^2\,\mathrm dt<\infty$，且假设 (3.15) 成立，则

$$
\int_\Gamma\!\!\int_D \rho(y)\,u_{y_i}^2\,\mathrm dx\,\mathrm dy
\ \le\ C(T)\bigl(\sqrt{\lambda_i}+\sqrt{\mu_i}+\nu_i\bigr),
\qquad 0<t\le T .
$$

配套的混合导数界（式 (3.18)）在进一步假设下给出 $\int_\Gamma\int_D\rho(y)(\partial^2u/\partial y_i\partial y_j)^2\mathrm dx\,\mathrm dy\lesssim C(T)(\lambda_i+\mu_i+\lambda_i\mu_i+\lambda_i\nu_i)$。

> [!warning] 一处指数不对称
> (3.17) 按原文转录：前两族特征值带根号而第三族 $\nu_i$ 不带。这很可能是文中的排印不一致，因此 $\nu_i$ 的指数**按印刷形式核实、按意图未核实**。用这条估计做各向异性阶选择时，结论（低指标方向权重更大）不依赖这个指数，但若要拿它做定量比较则需回查原文。

论文另有一项被列为贡献的结果：给出保证解落在合适随机空间（$BV$ 与 $H^k$）的**充分条件**，把编号 2 从随机*变量*波速推广到随机*场*波速。

### 数值实验

例 4.1 的设定与测量如下：

| 项目 | 设定                                                                |
| ---- | ------------------------------------------------------------------- |
| 方程 | $\partial_t u=\kappa(y)\partial_x u$，$x\in[-1,1]$                  |
| 波速 | $\kappa(y)=0.5\sum_{n=1}^{4}n^{-2}y_n$，各 $y_n$ 在 $(-1,1)$ 上均匀 |
| 初值 | $u(x,0,y)=\sin(x)$                                                  |
| 边界 | $u(\pm1,t,y)=\sin(\pm1+\kappa(y)t)$，按 $\kappa$ 的符号施加         |
| 误差 | $e_{mean}$、$e_{std}$、$e_2$（定义见下）                            |
| 表 1 | $t=2$ 处的非均匀（各向异性）收敛结果                                |

三个误差量分别是

$$
e_{mean}=\max_x\bigl|\mathbb E(u_{num})-\mathbb E(u_{exa})\bigr|,
\quad
e_{std}=\max_x\bigl|\sigma_{u_{num}}-\sigma_{u_{exa}}\bigr|,
\quad
e_2=\max_x\bigl(\mathbb E(u_{num}-u_{exa})^2\bigr)^{1/2}.
$$

全部实验用均匀分布，理由是敏感性结果 (3.17)–(3.18) 与分布无关。

这个波速的取值范围可以直接算出：$0.5\sum_{n=1}^{4}n^{-2}\approx0.5\times1.4236\approx0.712$，故 $\kappa$ 在参数域上取遍 $[-0.712,\,0.712]$，**确实变号**，上文的逐模态符号判断因此是被真正用到的、而不是形式上的。

把 $\kappa_n$ 视为常数，则 $\sqrt{\lambda_n}\propto n^{-2}$，敏感性估计 (3.17) 中各方向的权重按 $n^{-2}$ 递减，与论文所取的各向异性阶一致：

| $i$ | 系数 $0.5\,i^{-2}$ | 相对权重 | 所取阶 $r_i$ | 该方向基数 $r_i+1$ |
| --- | ------------------ | -------- | ------------ | ------------------ |
| 1   | $0.5$              | $1$      | 3            | 4                  |
| 2   | $0.125$            | $1/4$    | 2            | 3                  |
| 3   | $\approx0.0556$    | $1/9$    | 1            | 2                  |
| 4   | $0.03125$          | $1/16$   | 1            | 2                  |

论文的结论是：即使取如此低的各向异性阶 $(r_1,r_2,r_3,r_4)=(3,2,1,1)$，数值解与精确解仍吻合良好；各向同性加密收敛很快；误差随 $t$ 增长，论文把这个对时间的依赖描述为通常呈线性，并指出这种「长时间积分」退化对广义多项式混沌与配点法同样是已知问题。

按论文自己的公式 $N_y=\prod_i(r_i+1)$，这组阶给出 $4\times3\times2\times2=48$ 个基函数，而在同样的最高阶 $r=3$ 上各向同性张量基需要 $4^4=256$ 个——**约五分之一的自由度**。（这两个数是用文中公式与文中阶组算出的，不是论文报告的数值。）本站未获得表 1 的具体误差数值，因此不报告。

### 与其他论文的关系

编号 5 是编号 2 的构造性 Galerkin 对应物，共享同一个模型问题；编号 4 随后给出覆盖两者的解析正则性与谱收敛率。编号 3 把完全相同的双正交解耦装置用到椭圆界面问题上。注 3.2 中「双重正交等于隐藏的配点」这一观察，是通往后面基于采样的论文的概念桥梁；各向异性阶选择规则则预示了[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|最小二乘]]、[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|最优采样]]与[[computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce|稀疏恢复]]诸篇中的各向异性与稀疏指标集选择。

## 3：把同一装置用于椭圆界面问题（无法核实）

> [!warning] 本篇无法核实
> Crossref、OpenAlex、MaRDI/zbMATH 与 Semantic Scholar 均不载该 DOI 的摘要（Semantic Scholar 的记录明确写明摘要字段已被出版社隐去），ScienceDirect 对自动检索返回 403。本站因此**不报告该文的任何定理、假设、常数、收敛阶或数值结果**，也不引用任何文字充当其摘要。下面只写可由 zbMATH 索引关键词、该文参考文献表与同组同期论文独立确认的内容。

可核实的题录是：Tao Zhou 独著，_J. Comput. Appl. Math._ 236(5)（2011 年 10 月 1 日），pp. 782–792，DOI `10.1016/j.cam.2011.05.033`。zbMATH/MaRDI 记录（Q651914）给出的索引关键词为「immersed finite element」「random elliptic interface problems」「bi-orthogonal polynomials」「Galerkin method」。参考文献表中包含 Li–Lin–Lin–Rogers（_Numer. Methods Partial Differential Equations_ **20** (2004) 338–367）与 He–Lin–Lin（同刊 **24** (2008) 1265–1300）两篇关于浸入有限元空间及其逼近能力的工作。Semantic Scholar 还挂有一段由模型生成的摘要式文字——**那不是作者的措辞，也不能替代摘要**，本站不引用它，只记录它与其他来源一致地佐证了两条结构性判断：双正交多项式是解耦装置，输出是一族**互不耦合**的确定性界面问题。

问题背景本身是标准的：带界面的椭圆问题在确定性情形已经困难，标准有限元除非网格贴合界面否则精度下降；加上随机输入后，朴素的随机 Galerkin 路线双重昂贵，因为投影给出耦合的确定性界面问题系统，而每一个都需要贴合界面的网格。上述两条可核实的判断合起来说明该文一次去掉两重耦合：随机方向用双正交基解耦（机制见编号 5），物理方向用浸入有限元空间，使网格不必与界面对齐。

不报告的部分需要说清楚：该文是否给出误差估计、给出什么形式的估计、常数与收敛阶如何，本站一概不知；二手记录显示文中附有数值实验，但测试问题与观察到的结果同样无法确认；文中用于选取各向异性多项式阶的敏感性估计的确切形式也未核实——编号 5 中已核实的定理 3.1／式 (3.17)–(3.18) 只能说明这一类装置的形态，不能代替该文自己的陈述。

## 7：带随机输入的延迟微分方程（摘要层面）

> [!note] 核实层级
> 本节内容限于摘要、关键词与文章前置信息。方程级细节、假设、常数与实验数据本站**未核实**，因此不报告。

### 直觉

带随机系数的延迟微分方程把两类困难叠在一起。延迟项使解在**时间**方向的光滑性在滞后量的各个整数倍处依次断裂，这就是经典的导数间断传播现象；随机系数则提出通常的问题：解在**参数**方向是否足够光滑，使随机空间里的谱方法划得来。这两者是有张力的：时间方向的 Legendre 谱配点要快，靠的正是解在区间上的光滑性，而延迟恰好在特定时刻破坏它。因此确立参数正则性在这里是前提而不是形式，论文也确实先做正则性、再提格式。

### 构造

随机空间用**随机配点**：在参数空间的每个配点上独立求解确定性延迟微分方程，再插值。MaRDI 关键词表中出现「sparse grid」，说明多个随机输入时配点取 Smolyak 型稀疏网格。每个确定性延迟微分方程在时间上用 **Legendre 谱配点方法**求解。论文先建立随机空间中的正则性，条件是「给定数据满足某些合理的假设」——与编号 2、4、5 同一模式：把数据假设转成参数正则性，再转成收敛率。该文这些假设的确切形式与正则性估计的形式本站未核实。

### 主要结论（摘要层面）

在对数据的合理假设下精确解在随机空间中具有良好的正则性；论文给出所提方法的收敛性分析；方法在**随机空间与时间空间上都达到常见的指数收敛阶**。本站核实过的来源里恢复不出任何显式常数、显式指数或样本复杂度关系，因此不给出具体速率常数。论文如何处理时间方向在滞后倍数处的导数间断，本站同样无法报告。

### 数值实验

摘要写明给出了数值算例以印证理论结果。具体求解了哪些延迟问题、取了什么滞后量、观察到的误差量级如何，本站未核实，不报告。

### 与其他论文的关系

方法论上这是 Zhou「先正则性、后谱收敛率」这条线的最后一篇，模式与编号 2、4、5 相同，只是对象从输运偏微分方程换成延迟常微分方程，也是本组中唯一以稀疏网格为采样装置的一篇。它处在整份名单的转折点上：到编号 7 为止问的是**给定**的结构化节点集能达到什么精度，从编号 9 起问的是如何**设计**节点集——见[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|最小二乘一页]]，其中同年发表、同样把稀疏网格列为候选的编号 6 是那篇过渡性的比较研究，其结论是结构化稀疏网格在很低维之外并不占优。

## 38：带随机输入的 Maxwell 方程组（摘要层面）

> [!note] 核实层级
> 只有摘要、元数据与参考文献表可获取，全文在付费墙后且未找到预印本。本节按摘要陈述转述其构造与结论；假设、常数、CFL 条件与收敛阶本站**未核实**，因此不报告。

### 直觉

含时 Maxwell 方程组是一阶双曲型系统，因此继承编号 2、4、5 关于双曲问题正则性的全部困难。它另有一件东西值得保住：**能量守恒律**。一个 Galerkin 截断能否继承这样的二次不变量，取决于投影产生的耦合系统的结构——也就是编号 1 所研究的那一类系数矩阵。论文关心的正是这两件事：投影是否保能量，以及耦合的 Galerkin 系统如何高效推进。

与编号 5 的对照很干净：两篇都是「把双曲系统用广义多项式混沌投影，再设法让得到的耦合确定性系统解耦」，但**编号 5 靠换基（双正交多项式）实现精确解耦，编号 38 靠改造时间推进格式实现解耦**。

### 构造（摘要层面）

先用广义多项式混沌把随机 Maxwell 系统转成关于展开系数的确定性 **Galerkin 系统**；物理空间用**有限元方法**离散；时间方向构造三种格式：**Crank–Nicolson** 格式、**经典蛙跳**格式，以及一种**改造的蛙跳型**格式，后者的设计目的就是让耦合的 Galerkin 系统能以解耦的方式推进。

### 主要结论（摘要层面）

- 随机 Galerkin 途径在广义多项式混沌投影这一层**保持能量守恒律**。
- 论文给出求解 Galerkin 系统所用有限元方法的误差估计。**收敛阶未核实。**
- 对 **Crank–Nicolson** 格式，**全离散**格式被证明是保能量的。
- 对**经典蛙跳**格式，给出**条件**能量稳定性。**显式的 CFL 型条件未核实。**
- **改造的蛙跳型**格式把 Galerkin 系统解耦，摘要称其「得到一个非常高效的数值途径」。**任何定量加速比未核实。**

### 数值实验

摘要写明给出数值算例支持理论结论，文章含三幅插图。具体测试问题与定量结果本站未核实，不报告。

### 与其他论文的关系

这是侵入式 Galerkin 一脉的最后一篇，直接引用了本组中的四篇：编号 3（椭圆界面问题的随机 Galerkin）、编号 5（双正交多项式的随机双曲 Galerkin）、编号 2（标量双曲方程随机配点的收敛分析）与编号 14（非结构网格上的随机配点，见[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|最小二乘一页]]）。其参考文献表中还有 Xiu–Shen 的《Efficient stochastic Galerkin methods for random diffusion equations》，即编号 1 所分析的那些系数矩阵的出处。

> [!note] 覆盖进度
> 编号 **1、2、4、5** 已按原文全文核对，本页给出它们的直觉、推导链条、带假设的定理与数值实验设定（编号 1 本身没有数值实验，全文是分析）。编号 **7、38** 只到摘要与元数据层面：构造与结论按摘要陈述转述，其假设、常数、CFL 条件、收敛阶与实验数据本站不报告。编号 **3** 是唯一完全无法核实的一篇：Crossref、OpenAlex、MaRDI/zbMATH 与 Semantic Scholar 均无其摘要（Semantic Scholar 明确记录摘要字段被出版社隐去），ScienceDirect 对自动检索返回 403，因此本页不给出它的任何定理、常数、收敛阶或数值结果，只报告可由 zbMATH 索引关键词、参考文献表与同组同期论文独立确认的构造思路。
> 若干处按原文转录但存在问题的公式已就地标出：编号 1 归一化递推系数 $\tilde a_n$ 的分母、编号 4 的解析延拓级数与常数 $r$ 的闭式、编号 4 数值算例三条信息之间的不自洽、编号 5 式 (3.17) 中 $\nu_i$ 的指数。

## 覆盖核对

| 内容                                     | 论文  | 覆盖状态                                         |
| ---------------------------------------- | ----- | ------------------------------------------------ |
| 特征线分离与精度随时间流失               | 2/4/5 | 显式解、$\partial_y^k u$ 的 $t^k$ 增长与后果     |
| Galerkin 系数矩阵与其稀疏性              | 1     | $a_{jk}$、$e_{ijk}$、对称性、$2N+1$ 非零元的推导 |
| 归结为三项递推系数                       | 1     | 两条恒等式、其推导与技术意义                     |
| Jacobi 递推系数与 $b_n=0$ 的判据         | 1     | $h_n$、$a_n$、$b_n$、$c_n$，分母未核实的标注     |
| 定理 1 与阈值 $\lvert\alpha\rvert\ge1/2$ | 1     | 两类密度、结论、限定与指数约定                   |
| 变号波速与流入边界切换                   | 2     | 问题设定与正则性后果                             |
| 带 $1/\lvert c\rvert$ 权的数据假设       | 2     | 三条积分条件、能量推导与权的对账                 |
| 三条正则性定理与三条收敛定理             | 2     | $H^1$/$H^2$/$BV$ 与 $N^{-1}$/$N^{-2}$/均值       |
| 假设尖锐性的数值验证                     | 2     | 三组初值的阶、$BV$ 例中均方误差掉到约半阶        |
| 几何增长假设与解析延拓                   | 4     | 假设形式、Leibniz 归纳、$1/\sqrt\zeta$ 的来历    |
| 配点的指数收敛与其速率                   | 4     | 估计式、$\varrho$ 与 $\tau$ 的关系               |
| 双曲问题一般不解析的限定                 | 4     | 论文自陈的两条限定与「相当粗」的自评             |
| 数值算例与时间退化的定量读法             | 4     | 设定表、算例不自洽的标注、$\tau\sim1/t$          |
| 双正交多项式的两条定义条件               | 5     | 定义、特征值构造、与编号 1 同一矩阵              |
| 精确解耦与有效波速的符号判断             | 5     | (3.8)–(3.10) 的推导与 $C_{k,i}$ 落在支撑内       |
| 各向异性阶选择的敏感性依据               | 5     | (3.17)、$\nu_i$ 指数的标注、(3.18)               |
| 例 4.1 的设定与自由度节省                | 5     | 波速变号区间、权重表、$48$ 对 $256$              |
| 界面问题的双重解耦                       | 3     | 仅关键词与参考文献可确认的构造思路               |
| 延迟微分方程的配点与谱时间离散           | 7     | 摘要层面的构造与「两方向指数收敛」的结论         |
| Maxwell 的保能量与解耦时间推进           | 38    | 摘要层面的三种格式与各自的结构性质               |

## 本页原文

- T. Zhou and T. Tang, [_Note on coefficient matrices from stochastic Galerkin methods for random diffusion equations_](https://doi.org/10.1016/j.jcp.2010.07.016), J. Comput. Phys. 229 (2010), pp. 8225-8230。
- T. Tang and T. Zhou, [_Convergence analysis for stochastic collocation methods to scalar hyperbolic equations with a random wave speed_](https://doi.org/10.4208/cicp.060109.130110a), Commun. Comput. Phys. 8 (2010), pp. 226-248。
- T. Zhou, [_Stochastic Galerkin methods for elliptic interface problems with random input_](https://doi.org/10.1016/j.cam.2011.05.033), J. Comput. Appl. Math. 236 (2011), pp. 782-792。
- T. Zhou and T. Tang, [_Convergence analysis for spectral approximation to a scalar transport equation with a random wave speed_](https://doi.org/10.4208/jcm.1206-m4012), J. Comput. Math. 30 (2012), pp. 643-656。
- T. Zhou and T. Tang, [_Galerkin methods for stochastic hyperbolic problems using bi-orthogonal polynomials_](https://doi.org/10.1007/s10915-011-9508-0), J. Sci. Comput. 51 (2012), pp. 274-292。
- T. Zhou, [_A stochastic collocation method for delay differential equations with random input_](https://doi.org/10.4208/aamm.2012.m38), Adv. Appl. Math. Mech. 6 (2014), pp. 403-418。
- Z. Feng, J. Li, T. Tang, and T. Zhou, [_Efficient stochastic Galerkin methods for Maxwell's equations with random inputs_](https://doi.org/10.1007/s10915-019-00936-z), J. Sci. Comput. 80 (2019), pp. 248-267。
