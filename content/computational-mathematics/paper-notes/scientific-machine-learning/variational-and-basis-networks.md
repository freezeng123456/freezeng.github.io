---
title: 约束怎么进损失：增广 Lagrange 与神经网络基
description: 编号 60、90、102：把边界条件、局部低正则性与散度约束交给结构而不是罚项
lang: zh
translation: en/computational-mathematics/paper-notes/scientific-machine-learning/variational-and-basis-networks
tags:
  - 论文笔记
  - 科学机器学习
  - 神经网络基
---

> [!note] 本页覆盖
> 编号 **60**（_Commun. Comput. Phys._ 31(3), 2022）、**90**（_Commun. Comput. Phys._ 39(2), 2026）、**102**（投稿 _Commun. Comput. Phys._，[arXiv:2603.17906](https://arxiv.org/abs/2603.17906)）。三篇的共同作者线索是黄建国与周涛，共同技术取向是**把网络训练还原成（一列）最小二乘问题**，并尽量让约束由结构而不是罚项承担。

## 60：用增广 Lagrange 取代罚项

### 问题

带本质边界条件（Dirichlet 型）的变分问题对网络求解器不友好：网络函数无法精确满足约束，连在插值节点上也不行。论文指出两种通行做法各有缺陷。一种是把约束写进拟设 $\varphi(x;\theta)=\ell(x)\psi(x;\theta)+\bar g(x)$，这需要已知一个类水平集函数 $\ell$ 与一个已知延拓 $\bar g$；另一种是在边界上加二次罚项，精度对罚参数敏感。本文换成真正的鞍点表述，使罚参数不再需要激进调节。

### 增广 Lagrange 与鞍点

设 $V$ 为 Hilbert 空间，$B:V\to W$ 有界线性（通常 $W=L^2(\Gamma)$），$V_g=\{v\in V:Bv=g\text{ on }\Gamma\}$，原问题为 $\min_{v\in V_g}J(v)$。增广 Lagrange 函数取

$$
\mathcal L_\beta(v,\mu)=J(v)-\langle\mu,\,Bv-g\rangle_W+\frac{\beta}{2}\|Bv-g\|_W^2 .
$$

注意乘子项前是**减号**，这是论文实际采用的符号约定。相应的极小极大问题 $\min_{v\in V}\max_{\mu\in W}\mathcal L_\beta(v,\mu)$ 与原问题等价，因为

$$
\max_{\mu\in W}\mathcal L_\beta(v,\mu)=
\begin{cases}
J(v), & v\in V_g,\\
+\infty, & \text{否则},
\end{cases}
$$

对偶函数为 $F_\beta(\mu)=\min_{v\in V}\mathcal L_\beta(v,\mu)$。

对二阶椭圆问题 $-\operatorname{div}(A\nabla u)+cu=f$、$u|_\Gamma=g$，

$$
J(v)=\tfrac12\int_\Omega\bigl[A\nabla v\cdot\nabla v+cv^2-2fv\bigr]\mathrm dx,
$$

$$
\mathcal L_\beta(v,\mu)=J(v)-\int_\Gamma\mu(x)\bigl(v(x)-g(x)\bigr)\mathrm dx
+\frac{\beta}{2}\int_\Gamma\bigl(v(x)-g(x)\bigr)^2\mathrm dx .
$$

线性特征值问题 $-\nabla\cdot(p\nabla u)+qu=\rho u$ 通过归一化 $\tilde v=v/\|v\|_{0,\Omega}$ 化为极小化 $\tilde J(v)=\int_\Omega[p\nabla\tilde v\cdot\nabla\tilde v+q\tilde v^2]\mathrm dx$，非线性（Gross-Pitaevskii 型）情形加入 $\tilde v^4$ 项，边界项形式不变。

### 乘子更新在参数空间的实现

原始与对偶变量各用一个 ResNet：$u\approx\varphi^u(x;\theta_u)$，$\lambda\approx\varphi^\lambda(x;\theta_\lambda)$，网络定义为 $h_0=Vx$、$h_\ell=h_{\ell-1}+\sigma(W_\ell h_{\ell-1}+b_\ell)$、$\varphi=a^{\mathsf T}h_L$。

这里出现了本文的核心技术困难：无穷维乘子更新

$$
\mu_{k+1}=\mu_k-\beta_k(Bv_k-g)
$$

无法直接作用在网络**参数**上——右端是一个函数，而参数空间没有对应的加法。论文的解决办法是在参数空间做一次最小二乘投影：

$$
\theta^\mu_{k+1}=\arg\min_{\theta_\nu}
\bigl\|\varphi^\nu-\varphi^\mu_k-\beta_k\bigl(B\varphi^v_k-g\bigr)\bigr\|^2_{L^2(\Gamma)},
$$

其中 $\varphi^\nu=\varphi^\lambda(x;\theta_\nu)$、$\varphi^\mu_k=\varphi^\lambda(x;\theta^\mu_k)$、$\varphi^v_k=\varphi^u(x;\theta^v_k)$。论文称之为投影技术。它把「函数空间的显式更新」翻译成「参数空间的拟合问题」，这一步是把经典增广 Lagrange 算法搬到网络上的关键。

两个损失的 Monte Carlo 形式为

$$
\mathcal L_\beta(v,\mu)=|\Omega|\,\mathbb E_\xi\Bigl[\tfrac12 A(\xi)|\nabla v(\xi)|^2
+\tfrac12 c(\xi)v^2(\xi)-f(\xi)v(\xi)\Bigr]
-|\Gamma|\,\mathbb E_\eta\Bigl[\mu(\eta)\bigl(v(\eta)-g(\eta)\bigr)
-\tfrac{\beta}{2}\bigl(v(\eta)-g(\eta)\bigr)^2\Bigr],
$$

$$
J_\lambda(\nu;\mu_k,v_k)=|\Gamma|\,\mathbb E_\eta
\Bigl[\bigl(\nu(\eta)-\mu_k(\eta)-\beta(v_k(\eta)-g(\eta))\bigr)^2\Bigr],
$$

$\xi\sim\mathrm{Unif}(\Omega)$、$\eta\sim\mathrm{Unif}(\Gamma)$。

外层循环是：固定 $v_k,\mu_k$ 用 Adam 更新乘子网络；固定 $\beta_k,\mu_{k+1}$ 用 Adam 更新解网络；再令 $\beta_{k+1}=\alpha\beta_k$（$\alpha\ge1$）。**罚参数是递增的**，但因为乘子承担了大部分约束力，起始 $\beta_0$ 不必很大，这正是「罚参数选择灵活」的机制。

论文给出的是无穷维层面的经典等价性与鞍点刻画，并未对离散（网络）方法证明收敛或逼近定理。实验覆盖二维与三维的二阶椭圆问题、线性特征值问题与非线性特征值问题，基线是罚方法与直接用随机梯度下降上升解极小极大问题；论文在一条注记中明确指出，为该极小极大问题设计高效的下降上升法仍是开放问题。

## 90：把低正则性局部化，再给局部换一个尺度

### 随机基方法的失效模式

随机基 PDE 求解器（极限学习机、随机特征方法、可迁移神经网络）把隐层权重与偏置随机取定后冻结，只对输出层系数做最小二乘。解光滑时极其精确；论文直接指出，当精确解正则性低时精度显著恶化。本文针对的是**局部**正则性丢失：尖峰与凹角奇性，即特征长度尺度比区域尺度低若干量级的情形，此时全局统一尺度的随机基无法分辨。论文把「低正则」定义为不属于 $H^2$，或导数在某些点附近非常大。

浅层网络写成基的形式：

$$
u_{NN}(x)=\sum_{m=1}^{M}\alpha_m\,\sigma\bigl(w_m^{\top}x+b_m\bigr)+\alpha_0,
$$

由于 $\{w_m,b_m\}$ 预设并冻结，$u_{NN}$ 落在 $V_M=\mathrm{span}\{\psi_0,\dots,\psi_M\}$，其中 $\psi_0\equiv1$、$\psi_m(x)=\sigma(\omega_m^{\top}x+b_m)$。

预设隐层有两种方式。第一种是在 $[-R,R]^d$ 上均匀抽 $w_m$、在 $[-R,R]$ 上抽 $b_m$。第二种是可迁移神经网络的重参数化：把每个神经元写成

$$
\sigma\bigl(w_m^{\top}x+b_m\bigr)=\sigma\bigl(\gamma_m(a_m^{\top}x+r_m)\bigr),
\qquad \|a_m\|_2=1,\ \gamma_m\ge0,
$$

对应 $w_m=\gamma_ma_m$、$b_m=\gamma_mr_m$。**位置**参数取

$$
a_m=\frac{X_m}{\|X_m\|_2},\qquad r_m=U_m,
$$

$X_m$ 是独立 $d$ 维标准高斯、$U_m$ 是 $[0,1]$ 上独立均匀，这使分割超平面在单位球内均匀分布；**形状**参数共享，$\gamma_m\equiv\gamma$，并用高斯随机场的实现作为辅助函数调节，不使用 PDE 信息。这个分解把「随机基好不好」拆成了两个可分别控制的部分：超平面位置的覆盖性与激活的有效波长。

### 非重叠区域分解与无权最小二乘

把 $\bar\Omega=\cup_{k=0}^{K}\bar\Omega_k$ 取为内部互不相交，每个 $\Omega_k$（$k\ge1$）只与 $\Omega_0$ 沿 $\Gamma_k=\partial\Omega_k\cap\partial\Omega_0$ 相接。原边值问题化为带 $C^1$ 传输条件的子域系统：

$$
Lu_k=f\ \text{in}\ \Omega_k,
\qquad u_k=g\ \text{on}\ \partial\Omega_k\cap\partial\Omega,
$$

$$
u_k=u_0\ \text{on}\ \Gamma_k,
\qquad
\frac{\partial u_k}{\partial n_k}=\frac{\partial u_0}{\partial n_k}\ \text{on}\ \Gamma_k .
$$

配点最小二乘损失把全部残差以**系数 1** 相加，没有罚权：

$$
\min_{\alpha_0,\dots,\alpha_K}\
\sum_{k=0}^{K}\Bigl(\sum_{x\in X_{f_k}}\bigl|L\tilde u_k(x)-f(x)\bigr|^2
+\sum_{x\in X_{g_k}}\bigl|\tilde u_k(x)-g(x)\bigr|^2\Bigr)
+\sum_{k=1}^{K}\sum_{x\in X_{\Gamma_k}}
\Bigl(\bigl|\tilde u_k(x)-\tilde u_0(x)\bigr|^2
+\Bigl|\tfrac{\partial\tilde u_k}{\partial n_k}(x)-\tfrac{\partial\tilde u_0}{\partial n_k}(x)\Bigr|^2\Bigr).
$$

这一点与带权 PINN 损失形成对照：因为基已冻结，线性算子下这就是一个线性最小二乘问题 $\min_\alpha\|F\alpha-T\|_2^2$，直接求解即可，权重不再是自由参数。

半线性算子用 Gauss-Newton：令 $u_k^{n+1}=u_k^{n}+v_k^{n}$，增量解线性化传输问题

$$
DL(u_k^n;v_k^n)(x)=f(x)-Lu_k^n(x)\ \text{in}\ \Omega_k,
$$

$$
v_k^n=g-u_k^n\ \text{on}\ \partial\Omega_k\cap\partial\Omega,
\qquad
v_k^n-v_0^n=u_0^n-u_k^n\ \text{on}\ \Gamma_k,
\qquad
\frac{\partial v_k^n}{\partial n_k}-\frac{\partial v_0^n}{\partial n_k}
=\frac{\partial u_k^n}{\partial n_k}-\frac{\partial u_0^n}{\partial n_k}\ \text{on}\ \Gamma_k,
$$

其中 $DL(u;v)$ 是 $L$ 在 $u$ 处沿 $v$ 的 Gâteaux 导数。把增量在同一冻结基上展开，每个 Newton 步同样是线性最小二乘。

### 自适应的三个部件

**指标。** 在 $\Omega_0$ 上取平均残差

$$
\mathcal L_{\Omega_0}(\alpha_0^*)=\frac{1}{|X_{f_0}|}\sum_{x\in X_{f_0}}
\bigl[L\tilde u_0(\alpha_0^*,x)-f(x)\bigr]^2 ,
$$

新峰中心取逐点残差的最大值点

$$
x_{K+1}=\arg\max_{x\in X_{f_0}}
\Bigl|L\Bigl(\sum_{m=0}^{M_0}\alpha^*_{m,0}\psi_{m,0}\Bigr)(x)-f(x)\Bigr| ,
$$

随后把区域切成 $\Omega^{*}=\Omega\setminus B_r(x_{K+1})$ 与 $\Omega^{**}=B_r(x_{K+1})\cap\Omega$。

**局部基的重定心与重缩放。** 在新建子域 $\Omega_K=B_{r_K}(x_K)\cap\Omega$ 上

$$
\psi_{0,K}(x)=1,
\qquad
\psi_{m,K}(x,c_K)=\sigma\bigl(c_K\,w_{m,K}^{\top}(x-x_K)+b_{m,K}\bigr).
$$

这是多尺度网络思想的局部化版本：$c_K$ 越大，每个基函数的有效波长越短，尖峰因此可分辨。

**尺度选择。** $c_K$ 由整数 $1,\dots,10$ 上的穷举搜索确定：对每个候选构造缩放基，只用 $\Omega_K$ 上的四条局部方程（内部残差、边界数据、以及与**冻结的**当前 $\tilde u_0$ 的两个界面匹配条件）解一个局部最小二乘问题，比较结果取最优。穷举在这里是可行的，因为每次求解都是小规模线性最小二乘，而不是一次网络训练。

## 102：让散度约束由算子恒等式承担

### 问题

不可压流的网络求解器几乎都把 $\operatorname{div}u=0$ 作为罚项，论文直接指出两个后果：罚参数的选择是关键问题；无论用速度-压力还是速度-涡量表述，得到的都是耦合 PDE 系统，未知量必须同时求解，复杂度高。二维流函数技巧能让约束精确成立，但论文指出把它推广到三维是主要困难。目标因此是：在二维与三维都**构造性**无散，并且速度与压力可以**顺序**求解。

### 一般原则与两个恒等式

对线性约束 $\mathcal A(u)=0$，寻找满足 $\mathcal A\circ\mathcal G=0$ 的 $\mathcal G$，再取拟设 $u=\mathcal G(v)$。这里 $\mathcal A=\operatorname{div}$，$\mathcal G=\mathbf{curl}$。二维取

$$
\mathbf{curl}\,\psi:=\Bigl[\tfrac{\partial\psi}{\partial y},\,-\tfrac{\partial\psi}{\partial x}\Bigr]^{\mathrm T},
\qquad
\operatorname{curl}\,v:=\tfrac{\partial v_2}{\partial x}-\tfrac{\partial v_1}{\partial y}.
$$

论文的两条定理把对流项在势变量下重写。二维：对 $\phi\in C^3(\Omega)$、$u=\mathbf{curl}\,\phi$，

$$
\operatorname{curl}\bigl((u\cdot\nabla)u\bigr)=-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi .
$$

三维：对 $\phi\in[C^3]^3$、$u=\mathbf{curl}\,\phi$ 且 $\operatorname{div}\phi=0$，

$$
\mathbf{curl}\bigl((u\cdot\nabla)u\bigr)
=(\Delta\phi\cdot\nabla)\mathbf{curl}\,\phi-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi .
$$

证明路径是 $\mathbf{curl}((u\cdot\nabla)u)=(u\cdot\nabla)\omega-(\omega\cdot\nabla)u$ 加上
$\omega=\mathbf{curl}\,\mathbf{curl}\,\phi=\nabla(\nabla\cdot\phi)-\Delta\phi=-\Delta\phi$。第二个等式用到了 $\operatorname{div}\phi=0$，这也说明三维情形为何必须把无散条件加在势上。

### 解耦后的子问题

对二维 Stokes 方程取 $\operatorname{curl}$ 消去压力（$\operatorname{curl}\nabla p=0$），得到单变量四阶问题

$$
\nu\Delta^{2}\phi=\operatorname{curl}\,f\ \text{in}\ \Omega,
\qquad
\phi=\frac{\partial\phi}{\partial n}=0\ \text{on}\ \Gamma .
$$

夹紧边界条件来自无滑移：$u\cdot n_0=\nabla\phi\cdot\tau_0=0$ 与 $u\cdot\tau_0=-\nabla\phi\cdot n_0=0$，故 $\phi$ 在 $\Gamma$ 上为常数，可归一化为零。二维 Navier-Stokes 相应为

$$
\nu\Delta^{2}\phi-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi=\operatorname{curl}\,f .
$$

三维 Stokes 为

$$
\nu\Delta^{2}\phi=\mathbf{curl}\,f,
\qquad
\operatorname{div}\phi=0\ \text{in}\ \Omega,
\qquad
\phi\cdot n=0,\ \ \mathbf{curl}\,\phi=0\ \text{on}\ \Gamma,
$$

三维 Navier-Stokes 是本文真正新的表述：

$$
\nu\Delta^{2}\phi+(\Delta\phi\cdot\nabla)\mathbf{curl}\,\phi
-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi=\mathbf{curl}\,f,
$$

配同样的三条约束。论文把它与涡量-矢势的有限差分/有限元工作对比：后者需要同时逼近势与涡量，而这里只有单一未知量 $\phi$，因此称为纯流函数表述在三维的推广。

压力恢复是一个**一阶梯度系统**而不是 Poisson 方程：

$$
\nabla p=f+\nu\Delta u\ \ (\text{Stokes}),
\qquad
\nabla p=f+\nu\Delta u-(u\cdot\nabla)u\ \ (\text{Navier-Stokes}),
\qquad p(x_0)=0,
$$

其中 $x_0$ 是任一内点，用于固定可加常数。论文指出这允许恢复经典解，而不是先前工作常用的弱解框架。

### 离散化与代价

基与编号 90 完全相同，向量值情形下**各分量共享同一组隐层基**，只有输出系数独立。二维 Stokes 速度步在内点与边界点上配置三组方程，得到超定线性系统并用最小二乘求解；压力步解

$$
\min_{\beta}\Bigl\|
\begin{bmatrix}\partial_x\Psi_{in}\\ \partial_y\Psi_{in}\\ \Psi_0^{\mathrm T}\end{bmatrix}\beta
-\begin{bmatrix}\tilde F_1\\ \tilde F_2\\ 0\end{bmatrix}\Bigr\|_2^2,
$$

最后一行就是 $p(x_0)=0$。Navier-Stokes 用 Gauss-Newton 线性化

$$
(u\cdot\nabla)u\ \longrightarrow\
u^{(k)}\!\cdot\!\nabla u^{(k+1)}+u^{(k+1)}\!\cdot\!\nabla u^{(k)}-u^{(k)}\!\cdot\!\nabla u^{(k)},
$$

论文称之为 Newton 加线性最小二乘框架，并指出它等价于先离散出非线性系统再用 Gauss-Newton；附录给出的 Picard/Oseen 型线性化用作初始化。

代价核算是解耦带来收益的定量说明。设 $\min_x\|Ax-b\|_2^2$、$A\in\mathbb R^{m\times n}$、$m\gg n$ 的代价为 $\mathcal O(mn^2)$，则

| 维数 | 解耦无散网络                                            | 耦合基线                      |
| ---- | ------------------------------------------------------- | ----------------------------- |
| 二维 | $\mathcal O((I+2J)(M+1)^2)+\mathcal O((2I+1)(M+1)^2)$   | $\mathcal O((3I+2J)(3M+3)^2)$ |
| 三维 | $\mathcal O((4I+4J)(3M+3)^2)+\mathcal O((3I+1)(M+1)^2)$ | $\mathcal O((4I+3J)(4M+4)^2)$ |

节省来自 $(M+1)^2$ 与 $(3M+3)^2$ 或 $(4M+4)^2$ 的差别：顺序解两个小系统，而不是一次解含 3 到 4 个耦合场的大系统。

> [!warning] 预印本中的一处符号
> 论文式 (4.11) 在线性化后的方程里印成 $\mu$，而上下文一贯使用黏性系数 $\nu$。按语境应为 $\nu$。

## 三篇的共同判断

| 编号 | 被移出罚项的约束   | 承担者                       | 训练变成什么           |
| ---- | ------------------ | ---------------------------- | ---------------------- |
| 60   | 本质边界条件       | 乘子网络 + 增广 Lagrange     | 交替的两个随机优化问题 |
| 90   | 界面匹配与局部尺度 | 区域分解 + 局部重缩放基      | 一列线性最小二乘       |
| 102  | 无散条件           | $\mathbf{curl}$ 拟设与恒等式 | 顺序的两个最小二乘     |

三条都在回答同一个问题：**当约束以罚项形式进入损失时，权重成为一个必须调的自由参数；把约束交给结构，这个自由参数就消失了。** 代价各不相同：编号 60 多一个网络与一次投影，编号 90 多一层区域分解与一次尺度搜索，编号 102 把二阶方程升成四阶方程。这些代价是可预测的，而罚权的调节不是。

## 覆盖核对

| 内容                            | 论文 | 覆盖状态                                    |
| ------------------------------- | ---- | ------------------------------------------- |
| 增广 Lagrange 与极小极大等价性  | 60   | 函数形式、符号约定、等价性与对偶函数        |
| 三类实例的具体能量              | 60   | 椭圆、线性特征值、非线性特征值              |
| 乘子更新的参数空间投影          | 60   | 困难来源、投影目标、外层循环与 $\beta$ 递增 |
| 冻结基与两种预设方式            | 90   | 基的定义、均匀抽样、重参数化与位置/形状分离 |
| 区域分解与无权最小二乘          | 90   | 传输条件、损失、线性与半线性情形            |
| 残差指标、局部重缩放与尺度搜索  | 90   | 三个部件及其可行性理由                      |
| 无散拟设与两条对流恒等式        | 102  | 一般原则、二维与三维恒等式及证明路径        |
| 四个解耦子问题与压力恢复        | 102  | 四阶速度方程、边界条件来源、一阶梯度系统    |
| 离散化、Gauss-Newton 与代价核算 | 102  | 共享基、最小二乘形式、线性化、复杂度对照    |

## 本页原文

- J. Huang, H. Wang, and T. Zhou, [_An augmented Lagrangian deep learning method for variational problems with essential boundary conditions_](https://doi.org/10.4208/cicp.OA-2021-0176), Commun. Comput. Phys. 31(3) (2022), pp. 966-986（预印本 [arXiv:2106.14348](https://arxiv.org/abs/2106.14348)）。
- J. Huang, H. Wu, and T. Zhou, [_Adaptive neural network basis methods for partial differential equations with low-regular solutions_](https://doi.org/10.4208/cicp.OA-2024-0310), Commun. Comput. Phys. 39(2) (2026), pp. 553-577（预印本 [arXiv:2411.01998](https://arxiv.org/abs/2411.01998)）。
- J. Cheng, J. Huang, H. Wang, and T. Zhou, _Decoupled divergence-free neural networks basis method for incompressible fluid problems_, [arXiv:2603.17906](https://arxiv.org/abs/2603.17906)，投稿 Commun. Comput. Phys.
