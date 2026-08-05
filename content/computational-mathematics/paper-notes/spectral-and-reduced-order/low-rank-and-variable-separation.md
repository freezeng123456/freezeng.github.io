---
title: 低秩流形与变量分离
description: 编号 15、17、83：让基与系数同时随时间演化，代价写在最小奇异值上
lang: zh
translation: en/computational-mathematics/paper-notes/spectral-and-reduced-order/low-rank-and-variable-separation
tags:
  - 论文笔记
  - 降阶模型
  - 低秩逼近
---

![在低秩流形上直接演化随机解](assets/diagrams/tao-zhou-papers/zh/low-rank-dynamics.svg)

## 低秩流形是降维，不是截断

三篇论文共享一个容易被误读的出发点，值得先说清楚：**把解限制在秩 $S$ 的表示上，与把一个固定展开截到前 $S$ 项，是两件不同的事。**

截断的做法是先把基定下来——随机方向上的广义多项式混沌多项式，或者由快照、由初值算出的本征正交分解模态与 Karhunen-Loève 模态——保留 $S$ 项，丢掉尾巴。此后未知量只剩下**固定线性子空间**里的 $S$ 个系数，精度完全由尾项在这组固定基下衰减得多快决定。当解对随机参数的依赖随时间明显改变时，那个固定子空间就不再是对的子空间，唯一的补救是加项：这正是编号 15 开篇所说的项数增长，也是编号 83 开篇所说的 Kolmogorov 障碍。

低秩流形的做法把状态空间换掉了。取 $\mathcal M_S$ 为全体秩 $S$ 场的集合，它不是线性子空间而是一张光滑流形；未知量是**在这张流形上运动的一个点**，表示的两个因子都是未知量。它的动力学由一个 Galerkin（Dirac-Frenkel）条件给出：要求残差 $\partial_t u_S-\mathcal L(u_S)$ 与当前点处的切空间 $T_{u_S}\mathcal M_S$ 正交，即

$$
\frac{\partial u_S}{\partial t}=P_{u_S}\bigl(\mathcal L(u_S)\bigr).
$$

由此有三个后果，恰好是这三篇论文的全部主题。

第一，解所占据的子空间可以在**秩不变**的前提下自行转动与变形，因为基本身是未知量的一部分。不必为了跟上一个移动的子空间而增大 $S$，Kolmogorov 障碍就是这样被绕开的。

第二，方程可以完全写在因子上：每个基函数一条偏微分方程，加上 $S$ 条系数方程。代价随 $S$ 与一次确定性场求解的代价增长，而不随乘积空间 $L^2(D)\otimes L^2(\Omega)$ 的维数增长；统计量也直接从因子读出，

$$
\mathbb E[u(x,t,\cdot)]\approx\bar u_S(x,t),
\qquad
\mathrm{Var}_T[u](t)\approx\sum_{i=1}^{S}\mathbb E\bigl[Y_i^2(t)\bigr],
$$

不需要再抽样。**这才是「降维」一词的实义：环境空间的维数不出现在演化的代价里。**

第三是代价。$\mathcal M_S$ 是弯的，它在一点处的曲率约为 $1/\sigma_S$，即最小保留奇异值的倒数；切空间投影里会显式地出现协方差矩阵的逆。后面的每一件事——编号 15 的 $e^{C/\rho}$ 常数、它对奇异协方差的 $C^{\dagger}$ 变通、编号 83 的整个设计动机——都从这一个逆矩阵长出来。

编号 17 属于这一专题的理由也在这里。逐时刻的最佳秩 $S$ 对象是截断 Karhunen-Loève 展开，它由协方差算子的 Fredholm 积分特征值问题 $\mathcal T_u Z_i=\mu_i Z_i$ 定义。它是编号 15 定理里的参照物，但要在每个 $t$ 上把它算出来，就得先知道 $u(t)$；而且 $t\mapsto z_S(t)$ 未必连续可微。于是只有两条路：把那个稠密特征问题算得更便宜（编号 17），或者用流形上的一条流代替逐时刻的特征分解（编号 15 与 83）。即便走后一条路，$t=0$ 处仍要做一次特征求解，因为动态正交流是从 $u_0$ 的 Karhunen-Loève 展开出发的。

## 15：动态正交逼近的误差理论

### 直觉：固定基为什么失效，动力低秩给出什么替代

时间依赖随机 PDE 通常用**固定**基展开：随机变量方向用广义多项式混沌，空间方向用本征正交分解模态。困难在于解对随机参数的依赖会随时间显著改变，因此固定基需要不断增加项数才能维持精度；自适应或贪心修补（时间依赖多项式混沌、推广的本征正交分解）只能部分缓解。

Sapsis 与 Lermusiaux 的动态正交逼近让空间基与随机系数**同时**演化，从而规避这一增长，但此前几乎没有误差理论：已有结果针对相近的多构型时间依赖 Hartree 方法与动力低秩矩阵情形，而不是随机 PDE。本文的第一个目标因此是**建立动态正交逼近与动力低秩逼近之间的精确对应**，把矩阵情形的理论整段搬过来。对应关系可以列成一张字典：

| 动力低秩（Koch-Lubich，矩阵）             | 动态正交（随机场）                                                |
| ----------------------------------------- | ----------------------------------------------------------------- |
| $Y=USV^{T}\in\mathcal M_r$                | $u_S=\bar u_S+\mathbf U^{T}\mathbf Y\in\mathcal M_S$              |
| 左因子 $U$，列正交                        | 空间基 $\{U_i(\cdot,t)\}$，$\langle U_i,U_j\rangle=\delta_{ij}$   |
| 右因子 $V$ 与核 $S$                       | 随机系数 $\{Y_i(t,\cdot)\}$，$\mathbb E[Y_i]=0$                   |
| 规范条件 $U^{T}\dot U=0$                  | 规范条件 $\langle\partial_t U_i,U_j\rangle=0$                     |
| $S^{-1}$ 出现在 $\dot U$、$\dot V$ 方程中 | $C^{-1}$ 出现在切空间投影中，$\sigma(u_S)=\sqrt{\mathrm{eig}(C)}$ |

最后一行是全篇的枢纽：矩阵情形里让流在 $\sigma_r\to0$ 时变刚的那个 $S^{-1}$，在随机场情形里就是 $C^{-1}$。

### 问题设定：模型、拟设与三条规范条件

在开区域 $D\subset\mathbb R^d$（$1\le d\le3$）与完备概率空间 $(\Omega,\mathcal A,\mathbb P)$ 上考虑

$$
\frac{\partial u(x,t,\omega)}{\partial t}=\mathcal L\bigl(u(x,t,\omega),\omega\bigr),
\qquad x\in D,\ t\in[0,T],\ \omega\in\Omega,
$$

$$
u(x,0,\omega)=u_0(x,\omega),
\qquad
\mathcal B\bigl(u(\sigma,t,\omega)\bigr)=h(\sigma,t),\ \sigma\in\partial D .
$$

误差分析专门针对抛物情形，此时算子取

$$
\mathcal L(u):=\nabla\!\cdot\!(a\nabla u)+f,
\qquad
\mathcal L^{*}(\cdot):=\mathcal L(\cdot)-\mathbb E[\mathcal L(\cdot)] .
$$

后面定理里的 $a_{\min}$ 与 $a_{\max}$ 属于这里的扩散系数 $a$。拟设为

$$
u_S(x,t,\omega)=\bar u_S(x,t)+\sum_{i=1}^{S}U_i(x,t)\,Y_i(t,\omega),
$$

表示的唯一性由三条条件固定：

$$
\mathbb E[Y_i(t,\cdot)]=0,
\qquad
\langle U_i(\cdot,t),U_j(\cdot,t)\rangle=\delta_{ij},
\qquad
\Bigl\langle \frac{\partial U_i(\cdot,t)}{\partial t},\,U_j(\cdot,t)\Bigr\rangle=0,
$$

对 $1\le i,j\le S$ 与所有 $t$ 成立，$\langle u,v\rangle=\int_D uv\,\mathrm dx$。第三条是**规范条件**：基的时间导数被要求与基自身的张成空间正交，也就是不允许基在自己的张成空间内旋转。它消去的正是分解的多余自由度——同一个 $\mathcal M_S$ 上的点有无穷多种因子写法，若不固定这项自由度，因子方程就是欠定的。

值得记下这组条件的不对称：只有空间基被约束，随机侧没有。论文自己指出这一点，并给出把随机侧也正交化的**双重动态正交**写法 $u_S^{*}=\mathbf U^{T}A\tilde{\mathbf Y}$，分析中用的正是它。

边界与初值条件同样由拟设决定：

$$
\mathcal B\bigl(\bar u_S(\sigma,t)\bigr)=h(\sigma,t),
\qquad
\sum_{i=1}^{S}C_{ij}(t)\,\mathcal B\bigl(U_i(\sigma,t)\bigr)=0,
$$

$$
\bar u_S(x,0)=\mathbb E[u_0(x,\cdot)],
\qquad
U_i(x,0)=Z_{i0}(x),
\qquad
Y_i(0,\omega)=\bigl\langle u_0(\cdot,\omega)-\bar u_0,\,Z_{i0}\bigr\rangle,
$$

其中 $\{Z_{i0}\}$ 是初值的 Karhunen-Loève 空间模态：**动态正交流是从 $u_0$ 的 Karhunen-Loève 展开点火的。**

### 推导：三个演化方程与它们不同的角色

把 Galerkin 投影与规范条件结合，得到

$$
\frac{\partial\bar u_S}{\partial t}=\mathbb E\bigl[\mathcal L(u_S)\bigr],
$$

$$
\sum_{i=1}^{S}C_{ij}(t)\,\frac{\partial U_i(x,t)}{\partial t}
=\Pi^{\perp}_{U}\,\mathbb E\bigl[\mathcal L(u_S)\,Y_j(t,\cdot)\bigr],
\qquad j=1,\dots,S,
$$

$$
\frac{\partial Y_i(t,\omega)}{\partial t}
=\bigl\langle \mathcal L^{*}(u_S(\cdot,t,\omega),\omega),\,U_i(\cdot,t)\bigr\rangle,
\qquad i=1,\dots,S,
$$

其中 $C_{ij}(t)=\mathbb E[Y_iY_j]$ 是协方差矩阵，$\mathcal L^{*}(u,\omega)=\mathcal L(u,\omega)-\mathbb E[\mathcal L(u)]$，而

$$
\Pi^{\perp}_{U}[v]=v-\Pi_U[v]=v-\sum_{i=1}^{S}\langle v,U_i\rangle U_i
$$

是到 $\mathcal U=\mathrm{span}\{U_1,\dots,U_S\}$ 的正交补的 $L^2(D)$ 投影。

三个方程的角色不同：均值方程就是平均后的 PDE；基方程由算子与随机系数的相关性中**被投影出张成空间的那部分**驱动，各模态之间仅通过 $C$ 耦合；系数方程是把中心化算子对每个基函数作检验得到的、以 $\omega$ 为参数的常微分方程。三者合起来仍然只有 $S+1$ 个空间场与 $S$ 条以 $\omega$ 为参数的常微分方程，这就是前面所说的「代价随 $S$ 而不随环境维数」。

### 推导：切空间投影与 $C^{-1}$ 的来源

把中心化部分写成 $u_S^{*}=\mathbf U^{T}\mathbf Y$，它落在秩 $S$ 流形 $\mathcal M_S$ 上。切空间由 $\delta u_S^{*}=\mathbf U^{T}\delta\mathbf Y+\delta\mathbf U^{T}\mathbf Y$ 构成，在规范条件的类比 $\langle\delta\mathbf U,\mathbf U^{T}\rangle=\mathbf 0$ 与 $\mathbb E[\delta\mathbf Y]=\mathbf 0$ 下分解唯一。切空间上的正交投影为

$$
P_{u_S^{*}}(v)=P_{u_S^{*}}(v^{*})=\mathbf U^{T}\langle v^{*},\mathbf U^{T}\rangle
+\Bigl(\Pi^{\perp}_{U}\bigl\{\mathbb E[v^{*}\mathbf Y^{T}]\bigr\}\,C^{-1}\Bigr)^{T}\mathbf Y,
$$

其中 $C=\mathbb E[\mathbf Y\mathbf Y^{T}]$ 按秩 $S$ 的定义满秩。**这就是后来出现在误差界里的 $\rho^{-1}$ 的来源。** 互补投影可以因子化为 $P^{\perp}_{u_S^{*}}v=\Pi^{\perp}_{U}\otimes\Pi^{\perp}_{\tilde{\mathcal Y}}v$，其中 $\tilde{\mathcal Y}=\mathrm{span}\{\tilde Y_1,\dots,\tilde Y_S\}$ 来自双重动态正交写法——这正是矩阵情形中 $P^{\perp}(Y)B=P^{\perp}_UBP^{\perp}_V$ 的对应物。

$\rho$ 不是抽象量。秩 $S$ 随机场的奇异值定义为 $\sigma(u_S):=\sigma(A)=\sqrt{\mathrm{eig}(C)}$，$A$ 是 $C$ 的平方根，因此 **$\rho$ 就是 $\sqrt{\lambda_{\min}(C(t))}$ 的一个下界**，是计算中看得见摸得着的东西。

用 Dirac-Frenkel 形式写出来，动态正交逼近就是 Galerkin 条件

$$
\mathbb E\Bigl[\Bigl\langle \frac{\partial u_S}{\partial t}-\mathcal L(u_S),\,v\Bigr\rangle\Bigr]=0
\quad \forall v=\bar v+v^{*},\ (\bar v,v^{*})\in H\times T_{u_S^{*}(t)}\mathcal M_S,
$$

等价地

$$
\frac{\partial u_S}{\partial t}=\mathbb E[\mathcal L(u_S)]
+P_{u_S^{*}(t)}\bigl(\mathcal L^{*}(u_S)\bigr).
$$

这个写法把「同时演化基与系数」翻译成一句话：**把右端投影到流形的切空间上。**

### 定理：曲率引理

误差分析的第一块砖是曲率引理。矩阵原型（Koch-Lubich）是：设 $X\in\mathcal M_r$ 的最小非零奇异值满足 $\sigma_r(X)\ge\rho>0$，$Y\in\mathcal M_r$ 且 $\|Y-X\|\le\tfrac18\rho$，则对一切 $B$，

$$
\bigl\|\bigl(P(Y)-P(X)\bigr)B\bigr\|\ \le\ 8\rho^{-1}\,\|Y-X\|\cdot\|B\|_2,
\qquad
\bigl\|P^{\perp}(Y)(Y-X)\bigr\|\ \le\ 4\rho^{-1}\,\|Y-X\|^{2},
$$

范数为 Frobenius 范数（$\|B\|_2$ 为谱范数）。本文的 Lemma 3.1 是**同一个陈述、同样的常数 $8\rho^{-1}$ 与 $4\rho^{-1}$**，搬到 $\mathcal M_S\subset L^2(D)\otimes L^2(\Omega)$ 上；论文把它归于 Conte 与 Lubich，只作记号上的小改动。两个不等式说的是同一件事：切空间投影随点变化的速度、以及流形偏离其切平面的程度，都被 $\rho^{-1}$ 控制——**曲率与最小奇异值成反比，就在这里量化。**

### 定理 4.1：拟最优性与它的代价

参照对象是每个时刻的最佳 $S$ 秩逼近，即截断的 Karhunen-Loève 展开

$$
z_S(x,\omega,t)=\bar u(x,t)+\sum_{i=1}^{S}\sqrt{\mu_i(t)}\,\gamma_i(t,\omega)\,Z_i(x,t),
$$

其中 $\{\mu_i,Z_i\}$ 是协方差算子的特征对。定理的假设有两条：该最佳逼近在 $0\le t\le\bar t$ 上连续可微地存在于 $(H^2(D)\cap H^1_0(D))\otimes L^2(\Omega)$ 内，且其最小奇异值一致有下界

$$
\sigma\bigl(z_S(t)\bigr)\ \ge\ \rho\ >\ 0,\qquad \forall t\in[0,\bar t].
$$

则存在 $0<\hat t\le\bar t$，使动态正交解 $u_S=\bar u_S+u_S^{*}$ 在 $u_S(0)=z_S(0)$ 下对一切 $0<t\le\hat t$ 满足

$$
\|u_S(t)-z_S(t)\|_0^2+a_{\min}\int_0^t |u_S(\tau)-z_S(\tau)|_1^2\,\mathrm d\tau
\ \le\ 2\alpha\,e^{2\beta(t)}\int_0^t \|z_S(\tau)-u(\tau)\|_1^2\,\mathrm d\tau,
$$

$$
\beta(t)=4\rho^{-1}\int_0^t\Bigl(4\|\mathcal L^{*}(z_S)\|_0+\|\mathcal L^{*}(u)\|_0
+\|\mathcal L^{*}(u_S)\|_0+\|\dot z_S^{*}\|_0^2\Bigr)\mathrm d\tau,
\qquad
\alpha=\max\Bigl\{\frac{a_{\max}^2}{2a_{\min}},\ 4\rho^{-1}\Bigr\},
$$

其中 $\|\cdot\|_1$ 与 $|\cdot|_1$ 是 $H^1(D)\otimes L^2(\Omega)$ 的范数与半范数，且要求界中各项都有定义。证明沿用 Lubich 等人对时间依赖数据矩阵的动力低秩论证。

结论的读法有两层。**上层是好消息**：动态正交误差被最佳秩 $S$ 误差控制，方法是拟最优的。**下层是代价**：常数按 $e^{C/\rho}$ 增长，即随最小奇异值倒数指数增长。这就是曲率障碍的定量形式——曲率引理里的 $4\rho^{-1}$ 原封不动地进了 $\beta(t)$ 与 $\alpha$，当第 $S$ 个奇异值塌缩时切空间投影变得病态，界失去意义。

假设进入的位置也值得逐一点名。**正则性**方面，分析需要 $u\in L^2\bigl(\mathcal T,(H^2(D)\cap H^1_0(D))\otimes L^2(\Omega)\bigr)$；论文证明了动态正交解也满足相应的能量估计，前提是 $\nabla a\in L^{\infty}(D\times\Omega)$ 且 $u(0),\dot u(0)\in H^1(D)\otimes L^2(\Omega)$。**弱形式的可用性**方面有一处关键的技术点（Prop. 3.5）：$-\Delta u_S^{*}\in T_{u_S^{*}}\mathcal M_S$，因此它是动态正交弱形式中的合法检验函数——没有这一条，抛物问题的能量论证就无法在流形上进行。**$\rho$** 只经由曲率引理进入，**连续可微性**则决定了结论成立的时间区间，下一节说明它为什么不是技术性假设。

### 假设为什么不是技术性的：特征值交叉

论文构造了一个特征值交叉的例子：确定性 Laplacian 加随机初值，前两个 Karhunen-Loève 特征值在某时刻 $t^{*}$ 交叉。此时最佳秩 1 误差与动态正交秩 1 误差都能闭式算出：

$$
\epsilon_{KL}(t)=\min\bigl\{\mathbb E[\alpha_1^2]e^{-2t},\ \mathbb E[\alpha_2^2]e^{-8t}\bigr\},
\qquad
\epsilon_{KL}(t)=\mathbb E[\alpha_2^2]e^{-8t}\ \ (t>t^{*}),
$$

$$
\epsilon_{DO}(t)=\mathbb E[\alpha_1^2]e^{-2t}=
\begin{cases}
\epsilon_{KL}(t), & t\le t^{*},\\[4pt]
\dfrac{\mathbb E[\alpha_1^2]}{\mathbb E[\alpha_2^2]}\,e^{6t}\,\epsilon_{KL}(t), & t>t^{*}.
\end{cases}
$$

读法很直接：动态正交误差按**最小**的那个 Laplacian 特征值定的速率衰减，而 Karhunen-Loève 误差在交叉后改按第二个特征值定的速率衰减，两者之比按 $e^{6t}$ 增长。**因此动态正交误差不可能被 Karhunen-Loève 误差一致地控制住。** 论文强调这与 Theorem 4.1 不矛盾：在 $t^{*}$ 处，秩 1 截断 Karhunen-Loève 展开在时间上不连续可微，定理的假设本来就不成立。换句话说，那条可微性假设不是技术副产品，**它就是「不许交叉」这个要求本身。**

> [!note] 障碍是界的障碍，未必总是方法的障碍
> 矩阵情形的原始文献对这一障碍给出两条限定。其一，若真实的 $\varepsilon$ 伪秩为 $q<r$ 而按秩 $r$ 积分，尽管核矩阵病态、$S^{-1}$ 出现在因子方程中，逼近质量并未受到严重影响。其二，$Y(t)$ 的奇异值合并处并不产生奇性，这一点与光滑奇异值分解的常微分方程不同。真正的失效模式是另一种：秩取得太小时，一个**没有被纳入**逼近的奇异值可以越过正在跟踪的那些而方法察觉不到——这恰好就是本文构造的特征值交叉例子。

### 数值实验：三个测试问题

论文给出三个测试问题，覆盖三种不同的随机性来源。

| 算例                    | 方程与区域                  | 随机来源与初值               | 用途                                                                        |
| ----------------------- | --------------------------- | ---------------------------- | --------------------------------------------------------------------------- |
| §6.2 线性抛物方程       | 一维区间                    | 随机初值                     | 对比动态正交空间模态 $U_i$ 与 Karhunen-Loève 模态 $Z_i$，跟踪随机系数的方差 |
| §6.3 线性抛物方程       | $x\in[0,1]$，齐次 Dirichlet | 随机扩散系数，确定性初值     | 初始协方差为零，专门检验 $C^{\dagger}$ 的处理                               |
| §6.4 非线性反应扩散方程 | $D=[0,1]^2$，齐次 Neumann   | 随机阈值电位，确定性阶跃初值 | 行波解，波速正比于激发率 $\beta$                                            |

第一个算例中还包含一次把正则化阈值取在 $10^{-16}$ 量级的运行。第二个算例的具体设置是

$$
\partial_t u-\mathrm{div}\bigl(a(x,\omega)\partial_x u\bigr)=0,
\qquad
a(x,\omega)=\bar a(x)+\sum_{i=1}^{2}\bigl(\eta_{2i-1}(\omega)\cos(i\pi x)+\eta_{2i}(\omega)\sin(i\pi x)\bigr),
$$

$\bar a=1.45$，$\eta_1,\dots,\eta_4$ 为独立零均值均匀随机变量且 $\mathbb E[\eta_i^2]=\tfrac13\cdot10^{-i+1}$，确定性初值 $u_0(x)=10\sin(\pi x)$。第三个算例为

$$
\partial_t u-\mu\Delta u=f(u),
\qquad
f(u)=\beta\,u(u-1)\bigl(\alpha(\omega)-u\bigr),
$$

初值为阶跃（$x_1\le0.5$ 处取 $1$，其余取 $0$），阈值电位 $\alpha(\omega)$ 服从均匀分布，激发率 $\beta$ 为常数。

实验的定性结论是：测试确认了理论界，并划出了方法可用与不可用的边界——只要奇异值保持分离，动态正交解就紧跟截断 Karhunen-Loève 展开；一旦进入分析所预言的交叉与秩不足区域，它就按预言的方式退化。

### 实现陷阱：奇异协方差与伪逆

上表第二个算例针对的是一个容易忽略的实现问题。协方差矩阵 $C(t^n)$ 可能奇异或极度病态——例如任何带随机系数但**确定性初值**的系统在起始时刻都有 $C\equiv0$。直接用 Moore-Penrose 伪逆 $C^{\dagger}$ 会把「非活跃」的基函数置零，从而**阻止秩增加**。论文改为把基方程直接写成

$$
\frac{\partial\mathbf U}{\partial t}=C^{\dagger}\,\Pi^{\perp}_{U}\,
\mathbb E\bigl[\mathbf Y\,\mathcal L(u_S)\bigr],
$$

在 $C$ 满秩时它与原方程一致；实际计算中每步对协方差矩阵做对角化以解耦系统，因为即便 $C(0)$ 对角，动态正交流也不保持 $Y_i$ 的不相关性。

## 17：让特征求解只在最粗网格上进行

Karhunen-Loève 展开需要协方差核对应的 Fredholm 积分算子的特征对：给定 $\mathrm{Cov}(x,y)$，求

$$
\int_D \mathrm{Cov}(x,y)\,u(y)\,\mathrm dy=\lambda\,u(x),
\qquad \|u\|_{L^2(D)}=1 .
$$

与微分算子的特征问题相比，这里有一处结构性差别：积分算子 $\mathcal T$ 是**紧的、光滑化的**，特征值在 $0$ 处聚集（而不是趋于 $+\infty$，因此感兴趣的是**最大**的那些特征值），而离散化给出**稠密**矩阵。因此既没有稀疏性可用，装配与求解的代价又随分辨协方差核所需的网格分辨率急剧上升。

先前工作从离散化一侧攻这个代价——小波 Galerkin、推广的快速多极、谱元、带张量结构的 Legendre-Galerkin、紧积分算子的多层增广方法——但特征求解本身仍必须在最细网格上完成。

编号 17 取消的正是这一点。它把 Fredholm 特征值问题改造成一列在网格层次上逐级上行的**积分迭代**，每级只作用一次积分算子，而真正的特征求解始终被限制在**最粗**的空间里（最粗空间再加一个由上一级带来的方向）；细网格上只做积分，从不解特征问题。它给出误差估计与复杂度分析，总工作量与在最细网格上作一次积分相当，因此任何高效的积分方案都可以直接插进来。做法改造自多层修正方法（Lin 与 Xie）在微分算子特征问题上的框架，其中「细空间源问题加粗空间小特征问题」这一交替是共有的骨架。

## 83：贪心地逐项加基，使每一步都解耦

### 直觉：解耦买来什么，「动态」又比静态分离多出什么

论文把定位写得很清楚。广义多项式混沌写成 $u\approx\sum_i\zeta_i(\xi)g_i(x,t)$，**参数系数在时间上被冻结**、空间模态随时间变化，并且需要来自全模型的大量先验信息。本征正交分解写成 $u\approx\sum_i\zeta_i(t;\xi)g_i(x)$，**空间模态在时间与参数两个方向上都被冻结**。

本征正交分解能否成功取决于解流形 $\mathcal M=\{u(\cdot,t;\xi):t\in[0,T],\xi\in\Omega\}$ 的 Kolmogorov 宽度是否快速衰减，而这**对一阶线性输运问题与双曲问题失效**——即 Kolmogorov 障碍。编号 15 分析的动态正交方法通过让两个因子都动来突破这一障碍，其变体（DyBO、对偶 DO、DBO，后者与前两者相差一个可逆矩阵变换）也都如此；但它们都对全部 $N$ 个模态求解一个**耦合**系统，而这正是协方差求逆及其条件数问题的来源。

「动态」相对静态变量分离多出的东西可以说得很具体。静态变量分离取 $u\approx\sum_i\zeta_i(\xi)g_i(x,t)$：参数系数 $\zeta_i(\xi)$ **不含时间**，也就是说「某个参数值下各模态如何配比」这件事一次决定、全程不变，模态权重只能靠空间-时间场自己去承担。本文把它改成 $\zeta_i(t;\xi)$，配比在每个时刻重新决定，于是每一项的两个因子都随时间变化，整个表示成为低秩流形上的一条流而不是一个固定的分离展开。这不是一处修饰：下面 §5.1.2 的数值对比显示，静态版本随 $N$ 增大**发散**，而动态版本单调下降。

本文的角度因此是：保留双向时间依赖的拟设，但用**贪心加基逐项构造**模态，使每一步解耦成两个标量系数子问题，从而全程不出现协方差求逆。

### 问题设定：模型与仿射假设

$$
\frac{\partial u}{\partial t}(x,t;\xi)=F\bigl(u(x,t;\xi);\xi\bigr),
\qquad u(x,0;\xi)=\mu(x;\xi),
\qquad
\mathcal B\bigl(u(x,t;\xi)\bigr)=g(x,t;\xi),
$$

变分形式为 $\langle\partial_t u(\cdot,t;\xi),v\rangle=\langle F(u(\cdot,t;\xi);\xi),v\rangle$（$\forall v\in V$）。结构性假设是 $F(u;\xi)=\mathcal C(\xi)+\mathcal A(u;\xi)+\mathcal H(u;\xi)$（$\mathcal A$ 线性、$\mathcal H$ 非线性）且**参数依赖是仿射的**：

$$
\mathcal C(\xi)=\sum_{i=1}^{N_C}\kappa^i_C(\xi)\mathcal C^i,
\qquad
\mathcal A(u;\xi)=\sum_{i=1}^{N_A}\kappa^i_A(\xi)\mathcal A^i(u),
\qquad
\mathcal H(u;\xi)=\sum_{i=1}^{N_H}\kappa^i_H(\xi)\mathcal H^i(u),
$$

并要求初值可分离 $\mu(x;\xi)=\sum_{i=1}^{N_{t_0}}p^i(\xi)q^i(x)$。**仿射性正是离线在线划分得以成立的前提**；论文指出若它不成立，可以先用变量分离方法构造一个精度损失可忽略的仿射逼近。

### 问题设定：拟设既没有均值场，也没有正交规范

$$
u(x,t;\xi)\approx u_N(x,t;\xi):=\sum_{i=1}^{N}\zeta_i(t;\xi)\,g_i(x,t),
$$

$\{\zeta_i\}$ 依赖参数、$\{g_i\}$ 不依赖参数，**两者都随时间变化**。与动态正交方法的 $u\approx\bar u(x,t)+\sum_i\zeta_i(t;\xi)g_i(x,t)$ 对比，后者保留一个统计均值场并要求 $\{g_i(\cdot,t)\}$ 在每个 $t$ 上正交归一。**本文把均值场与正交规范同时去掉**，唯一性改由顺序贪心构造提供。这也解释了为什么编号 15 的定理无法搬过来：那套理论正是建立在均值场与规范条件之上的。

### 推导：贪心准则与两个解耦子问题

第一步任取 $\xi_1$，令 $g_1(x,t)$ 为方程在 $\xi=\xi_1$ 处的**全解**，再用 $g_1$ 作检验得到 $\zeta_1(t;\xi)$ 的参数依赖常微分方程

$$
\Bigl\langle \frac{\partial\bigl(g_1(x,t)\zeta_1(t;\xi)\bigr)}{\partial t},\,g_1(x,t)\Bigr\rangle
=\bigl\langle F\bigl(g_1(x,t)\zeta_1(t;\xi);\xi\bigr),\,g_1(x,t)\bigr\rangle .
$$

第 $k\ge2$ 步以 $e(x,t;\xi):=u-u_{k-1}$ 为误差，选点准则是

$$
\xi_k\in\arg\max_{\xi\in\Xi}\triangle_k(\xi),
$$

其中 $\triangle_k$ 取 $\|e\|_{L^2([0,T];V)}$ 本身（若代价允许）或下面的后验界；当 $\triangle_k(\xi_k)<\varepsilon$ 时停止，否则把 $\xi_k$ 从候选集中移除。把方程按误差改写为

$$
\Bigl\langle \frac{\partial(e+u_{k-1})}{\partial t},v\Bigr\rangle
=\bigl\langle F\bigl((e+u_{k-1});\xi\bigr),v\bigr\rangle,
\qquad \forall v\in V,
$$

则两个子问题是：$g_k(x,t)$ 取该式在**单个**参数 $\xi=\xi_k$ 处的解（与参数无关的偏微分方程），而令 $\tilde e=g_k\zeta_k$、$v=g_k$ 得到 $\zeta_k(t;\xi)$ 的标量常微分方程

$$
\Bigl\langle \frac{\partial\bigl(g_k\zeta_k+u_{k-1}\bigr)}{\partial t},\,g_k\Bigr\rangle
=\bigl\langle F\bigl(g_k\zeta_k+u_{k-1};\xi\bigr),\,g_k\bigr\rangle .
$$

**两者之间没有耦合：空间侧只解一个参数处的问题，参数侧只解一条标量方程。** 这就是「每步解耦」的全部含义。

### 推导：初值构造在 $t=0$ 处补回规范条件

论文称初值是「最关键的成分之一」，其构造用与当前模态的 $L^2$ 匹配。记 $g_{k,0}(x)=g_k(x,0)$，则 $g_1(x,0)=\mu(x;\xi_1)$，并由匹配条件 $\langle u_1(\cdot,0;\xi),g_{1,0}\rangle=\langle\mu(\cdot;\xi),g_{1,0}\rangle$ 得

$$
\zeta_{1,0}(\xi)=\sum_{i=1}^{N_{t_0}}\frac{\langle q^i,g_{1,0}\rangle}{\langle g_{1,0},g_{1,0}\rangle}p^i(\xi);
$$

对 $k\ge2$，初始误差 $e_0(x;\xi)=\mu(x;\xi)-\sum_{j=1}^{k-1}g_{j,0}(x)\zeta_{j,0}(\xi)$ 给出

$$
g_{k,0}(x)=e_0(x;\xi_k)=\mu(x;\xi_k)-\sum_{j=1}^{k-1}g_{j,0}(x)\,\zeta_{j,0}(\xi_k),
$$

$$
\zeta_{k,0}(\xi)=\sum_{i=1}^{N_{t_0}}\frac{\langle q^i,g_{k,0}\rangle}{\langle g_{k,0},g_{k,0}\rangle}p^i(\xi)
-\sum_{j=1}^{k-1}\frac{\langle g_{j,0},g_{k,0}\rangle}{\langle g_{k,0},g_{k,0}\rangle}\zeta_{j,0}(\xi).
$$

**第二个和式是一个类 Gram-Schmidt 的修正项。** 这正是本方法在 $t=0$ 处隐式地、局部地补回动态正交方法用规范条件显式提供的那部分结构。一个边界情形（论文的 Remark 3.1）：若 $\mu\equiv0$，则每一步 $g_{k,0}=0$ 且 $\zeta_{k,0}=0$。

### 推导：线性情形中新模态由上一步的残差驱动

取 $F(u;\xi)=\mathcal C(\xi)+\mathcal A(u;\xi)$，误差方程成为

$$
\Bigl\langle\frac{\partial e}{\partial t},v\Bigr\rangle-\bigl\langle\mathcal A(e;\xi),v\bigr\rangle
=\bigl\langle r_k,v\bigr\rangle,
\qquad
r_k:=\mathcal A(u_{k-1};\xi)-\frac{\partial u_{k-1}}{\partial t}+\mathcal C(\xi),
$$

即**新的空间模态由上一步逼近的残差驱动**。系数方程是标量常微分方程

$$
\bigl\langle (g_k)_t,g_k\bigr\rangle\zeta_k
+\bigl\langle g_k,g_k\bigr\rangle(\zeta_k)_t
-\bigl\langle\mathcal A(g_k;\xi),g_k\bigr\rangle\zeta_k
=\bigl\langle r_k,g_k\bigr\rangle .
$$

取 $u_0\equiv0$、$r_1\equiv\mathcal C(\xi)$ 时它退化为第一步的方程。注意 $\langle(g_k)_t,g_k\rangle$ 这一项：在动态正交方法里它被规范条件**强制为零**，而这里它只是被一路带着。这句话把两条路线的差别落到了一个具体的项上。

### 推导：离散递推的除数是标量而不是矩阵

把 $[0,T]$ 分成 $N_t$ 步、步长 $\tau=T/N_t$，记 $g_{k,n}=g_k(t_n)$、$\zeta_{k,n}(\xi)=\zeta_k(t_n;\xi)$，用一阶差分与后向 Euler，并代入仿射展开与 $u_{k-1}$ 的分离形式，得到闭式递推

$$
\zeta_{k,n+1}(\xi)=\frac{c_{n+1}\,\zeta_{k,n}(\xi)+s_{n+1}(\xi)}{l_{n+1}(\xi)},
\qquad n=0,\dots,N_t-1,
$$

$$
c_{n+1}=\frac{\langle g_{k,n+1},g_{k,n+1}\rangle}{\tau},
\qquad
l_{n+1}(\xi)=2c_{n+1}-\frac{\langle g_{k,n},g_{k,n+1}\rangle}{\tau}
-\sum_{i=1}^{N_A}\kappa^i_A(\xi)\bigl\langle\mathcal A^i(g_{k,n+1}),g_{k,n+1}\bigr\rangle,
$$

$$
\begin{aligned}
s_{n+1}(\xi)=&\sum_{i=1}^{N_A}\sum_{j=1}^{k-1}\kappa^i_A(\xi)\,\zeta_{j,n+1}(\xi)\bigl\langle\mathcal A^i(g_{j,n+1}),g_{k,n+1}\bigr\rangle
+\sum_{i=1}^{N_C}\kappa^i_C(\xi)\bigl\langle \mathcal C^i,g_{k,n+1}\bigr\rangle\\
&-\sum_{j=1}^{k-1}\frac{\zeta_{j,n+1}(\xi)-\zeta_{j,n}(\xi)}{\tau}\bigl\langle g_{j,n+1},g_{k,n+1}\bigr\rangle
-\sum_{j=1}^{k-1}\zeta_{j,n+1}(\xi)\Bigl\langle \frac{g_{j,n+1}-g_{j,n}}{\tau},\,g_{k,n+1}\Bigr\rangle .
\end{aligned}
$$

$s_{n+1}$ 完全由仿射项与前 $k-1$ 个已算出的模态组装。**这就是本文对动态正交方法协方差求逆的回答：$l_{n+1}(\xi)$ 是一个标量而不是矩阵，因此全程不出现 $C^{-1}$ 或 $C^{\dagger}$。**

离线阶段不存 $\zeta_{k,n}(\xi)$，只存与参数无关的标量内积

$$
\langle g_{k,n+1},g_{k,n+1}\rangle,\quad
\langle g_{k,n},g_{k,n+1}\rangle,\quad
\langle\mathcal A^i(g_{k,n+1}),g_{k,n+1}\rangle,\quad
\langle\mathcal A^i(g_{j,n+1}),g_{k,n+1}\rangle,
$$

$$
\langle g_{j,n+1},g_{k,n+1}\rangle,\quad
\langle g_{j,n},g_{k,n+1}\rangle,\quad
\langle\mathcal C^i,g_{k,n+1}\rangle .
$$

在线阶段对新参数 $\bar\xi$ 只需重解一条廉价的标量常微分方程

$$
\alpha(t)\,(\zeta_k)_t(t;\bar\xi)+\beta(t;\bar\xi)\,\zeta_k(t;\bar\xi)=\gamma(t;\bar\xi),
$$

$$
\alpha(t)=\langle g_k(t),g_k(t)\rangle,
\qquad
\beta(t;\bar\xi)=\langle (g_k)_t(t),g_k(t)\rangle-\sum_{i=1}^{N_A}\kappa^i_A(\bar\xi)\langle\mathcal A^i(g_k(t)),g_k(t)\rangle,
$$

$\gamma$ 由所存标量组装。$\beta$ 与 $\gamma$ 都关于参数仿射，因此**在线代价与原问题的空间离散无关**。

### 非线性情形：Burgers 的半隐处理

取 $\mathcal A(u;\xi)=\kappa(\xi)\partial_x^2u$、$\mathcal H(u;\xi)=-u\partial_xu$。离散后的残差方程含 $\langle e_{n+1}\partial_xe_{n+1},v\rangle$，论文用**半隐格式**把它近似为 $e_n\partial_xe_{n+1}$ 以控制代价，并说明同一手法可搬到 Allen-Cahn。

### 定理：后验估计与「没有收敛定理」

后验界基于 $F$ 在 $u$ 处的**局部对数 Lipschitz 常数**

$$
L_V[F](u):=\sup_{v\in V,\ v\ne u}\frac{\langle v-u,\;F(v;\xi)-F(u;\xi)\rangle}{\|v-u\|_V^2},
$$

对误差方程 $\partial_te=F(u;\xi)-F(u_{k-1};\xi)+r_k$（其中 $r_k=F(u_{k-1};\xi)-\partial_tu_{k-1}$）用 $e$ 作检验并配比较引理，得到 Proposition A.1：$\|e(x,t;\xi)\|_V\le\delta_k(t;\xi)$，

$$
\delta_k(t;\xi)=\int_0^{t}\alpha(s;\xi)e^{\int_s^{t}\beta(\tau;\xi)\mathrm d\tau}\mathrm ds
+e^{\int_0^{t}\beta(\tau;\xi)\mathrm d\tau}\|e(x,0;\xi)\|_V,
$$

其中 $\alpha(t;\xi)=\|r_k(x,t;\xi)\|_V$、$\beta(t;\xi)=L_V[F](u_{k-1}(x,t;\xi))$，选点准则取 $\triangle_k(\xi)=\int_0^T\delta_k^2\,\mathrm dt$。**对耗散型 $F$，对数 Lipschitz 常数可以为负**，此时指数因子起衰减而不是放大作用——这正是该估计在长时间区间上可用的原因。

论文的结论部分明确把「在合理假设下对本方法作严格收敛分析」列为未来工作：严格的内容就是上面的后验界与那条精确的代数递推，而「相对许多已有低秩分离技术降低复杂度、提高效率」由数值实验而非定理支持。论文也自陈一处限制：存储与在线计算仍依赖时空离散的规模，如何降低这一依赖仍是开放问题。

### 数值实验：四个算例

精度指标是 $M$ 个测试样本上的平均相对误差 $\epsilon$，参考解由空间有限元加时间后向 Euler（下称 FEM-BE）给出。

**算例一：参数依赖边界条件的反应扩散方程**（$D=[0,1]$）

$$
\partial_t u+\xi_1u=2\xi_2\,\partial_x^2u+\xi_3,
\qquad
u(x,0;\xi)=u|_{\partial D}=2(x+1)\xi_4,
$$

$T=1$，$\xi\in[1,3]^4$，网格 $h=0.02$，$\tau=10^{-3}$，训练集 $|\Xi|=11$，测试 $M=10^3$。论文同时比较了两种贪心策略：按真实误差 $e$ 选点与按估计量 $\delta_k$ 选点。对照方法是 Billaud-Friess 与 Nouy 的时间依赖约化基降阶方法（下称 MTD）。

| $N$ | MTD 的 $\epsilon$   | 本方法的 $\epsilon$ |
| --- | ------------------- | ------------------- |
| 2   | $2.48\times10^{-4}$ | $3.43\times10^{-4}$ |
| 4   | $9.19\times10^{-6}$ | $1.46\times10^{-4}$ |
| 7   | $2.80\times10^{-5}$ | $4.66\times10^{-5}$ |

| $N$ | 方法   | 离线 (s) | 在线 (s) | 合计 (s)  | 每样本在线 (s)      |
| --- | ------ | -------- | -------- | --------- | ------------------- |
| 2   | MTD    | $8.43$   | $19.93$  | $28.36$   | $1.99\times10^{-2}$ |
| 2   | 本方法 | $8.01$   | $1.98$   | $9.99$    | $1.98\times10^{-3}$ |
| 4   | MTD    | $18.78$  | $22.48$  | $41.26$   | $2.25\times10^{-2}$ |
| 4   | 本方法 | $16.51$  | $3.96$   | $20.47$   | $3.96\times10^{-3}$ |
| 7   | MTD    | $39.92$  | $118.52$ | $158.44$  | $1.19\times10^{-1}$ |
| 7   | 本方法 | $29.86$  | $7.20$   | $37.06$   | $7.20\times10^{-3}$ |
| —   | FEM-BE | —        | —        | $3752.14$ | $3.75$              |

这组数字给出一个坦率的权衡。**MTD 的误差下降更快**：在 $N=4$ 处它比本方法准 $15.9$ 倍。论文自己给出的解释是 MTD 在每次迭代中更新全部参数系数，而本方法不修改先前步已算出的时间参数基函数。但 MTD 的误差在 $N=4$ 处触底（$9.19\times10^{-6}$），到 $N=7$ 已回升到 $2.80\times10^{-5}$，此时优势只剩 $1.7$ 倍；论文指出在 $N=7$ 之后其误差快速上升，归因于其参数系数线性系统的病态。本方法的误差在三个 $N$ 上单调下降，且在线时间小一个数量级以上（$N=7$ 时每样本 $7.20\times10^{-3}$ 秒对 $1.19\times10^{-1}$ 秒）。论文的总结是：对照方法在某些条件下精度更高，本方法在计算效率与误差稳定性上更好。

**算例二：参数依赖源项的二维热方程**（$D=[0,\pi]^2$，$\partial_tu=\kappa(\xi)\Delta u+f$，$\kappa(\xi)=\xi_1$，$T=1$，$u(x,0;\xi)=\sin x_1\sin x_2+1$，$u|_{\partial D}=1$），对照对象是**静态基**变量分离方法 $u\approx\sum_i\zeta_i(\xi)g_i(x,t)$。下表为 $t=1$、$M=10^3$ 的结果，时间列为论文所报告的耗时。

| $N$ | 本方法的 $\epsilon$ | 时间 (s)           | 静态基的 $\epsilon$ | 时间 (s) |
| --- | ------------------- | ------------------ | ------------------- | -------- |
| 2   | $9.81\times10^{-4}$ | $0.71$             | $1.24\times10^{-1}$ | $0.47$   |
| 4   | $3.66\times10^{-4}$ | $2.37$             | $4.91\times10^{-2}$ | $1.52$   |
| 6   | $1.77\times10^{-4}$ | $4.98$             | $6.37\times10^{-2}$ | $3.12$   |
| 8   | $4.52\times10^{-5}$ | $8.52$             | $2.18$              | $5.28$   |
| 10  | $4.27\times10^{-5}$ | $13.01$            | $9.64\times10^{-1}$ | $8.00$   |
| —   | FEM-BE              | $7.06\times10^{2}$ | —                   | —        |

**这是全篇最干净的一处证据。** 静态基方法在 $N=4$ 处达到最好（$4.91\times10^{-2}$）之后掉头，$N=8$ 时平均相对误差 $2.18$ 已经超过 $1$，即逼近不再携带任何信息；本方法则单调下降到 $4.27\times10^{-5}$。论文据此说明「对参数与空间两个变量都使用时间依赖基函数」的必要性——这正是前面所说「动态」二字多出的那部分。

**算例三：参数依赖初值的 Burgers 方程**（$x\in[0,1]$，$T=2$）

$$
\partial_tu+u\,\partial_xu=\frac{\xi_1}{50}\,\partial_x^2u,
\qquad
u(x,0;\xi)=x(1-x)^2\xi_2,
\qquad
u(0,t)=u(1,t)=0,
$$

$\xi\in[1,3]^2$，$h=0.01$，$\tau=10^{-4}$，$|\Xi|=12$，$M=10^3$。

| $N$ | $\epsilon$（$t=1$） | 在线 (s)            | $\epsilon$（$t=2$） | 在线 (s)            |
| --- | ------------------- | ------------------- | ------------------- | ------------------- |
| 2   | $1.17\times10^{-2}$ | $1.63\times10^{-2}$ | $0.87\times10^{-2}$ | $2.92\times10^{-2}$ |
| 4   | $1.09\times10^{-3}$ | $5.48\times10^{-2}$ | $2.62\times10^{-3}$ | $9.87\times10^{-2}$ |
| 6   | $2.17\times10^{-4}$ | $1.23\times10^{-1}$ | $4.10\times10^{-4}$ | $2.17\times10^{-1}$ |
| 8   | $5.68\times10^{-5}$ | $2.18\times10^{-1}$ | $1.07\times10^{-4}$ | $3.82\times10^{-1}$ |
| 10  | $2.76\times10^{-5}$ | $3.44\times10^{-1}$ | $7.12\times10^{-5}$ | $5.92\times10^{-1}$ |
| —   | FEM-BE              | $10.72$             | FEM-BE              | $21.47$             |

两个时刻上误差都随 $N$ 单调下降。只有 $N=2$ 时 $t=2$ 的误差比 $t=1$ 小，$N\ge4$ 之后 $t=2$ 的误差一律高于 $t=1$，即项数够多以后长时间积分的误差累积才显出来。基场呈清楚的幅度层级：

| 基场 | $g_1$              | $g_2$            | $g_3$            | $g_5$            | $g_7$              | $g_9$            |
| ---- | ------------------ | ---------------- | ---------------- | ---------------- | ------------------ | ---------------- |
| 幅度 | $3.5\times10^{-1}$ | $2\times10^{-2}$ | $6\times10^{-3}$ | $1\times10^{-3}$ | $1.2\times10^{-4}$ | $3\times10^{-5}$ |

即首个场承载解的核心信息、末几个场承载细尺度信息。

**算例四：Allen-Cahn 方程**（$D=[0,1]^2$，$T=1$）

$$
\partial_tu=\xi^2\Delta u-f(u),
\qquad
f=F',\quad F(u)=\tfrac14(u^2-1)^2,
$$

初值 $u(x,0;\xi)=\sqrt5\,(x_1^2-x_1)(x_2^2-x_2)$，齐次 Dirichlet 边界，$\xi\in[0.1,0.2]$；$h_{x_1}=h_{x_2}=0.05$，$\tau=10^{-4}$，$|\Xi|=8$，$M=10^3$。结果与其余三例不同：平均相对误差先随项数下降，随后**趋于平台**，即出现精度地板。

## 三篇的关系

| 编号 | 演化对象                 | 模态之间的关系        | 主要风险                 |
| ---- | ------------------------ | --------------------- | ------------------------ |
| 15   | 空间基与随机系数         | 同时演化，经 $C$ 耦合 | 最小奇异值塌缩           |
| 17   | 不演化（一次性特征求解） | 不适用                | 稠密矩阵的装配与求解代价 |
| 83   | 逐项贪心加入的分离表示   | 每步解耦成两个方程    | 贪心序列不保证最优       |

编号 17 与编号 15 的关系值得指出：编号 15 的误差界以截断 Karhunen-Loève 逼近为参照，而编号 17 解决的正是**如何算出那个参照对象**——控制 Theorem 4.1 常数的那些奇异值 $\sqrt{\mu_i}$，恰好就是编号 17 要算的量。两篇一起给出了这条路线的完整代价结构：要么承担计算 Karhunen-Loève 基的稠密特征问题，要么让基自己演化并承担最小奇异值的风险。

编号 83 与编号 15 之间是一次有意识的交换。编号 15 证明动态正交方法对截断 Karhunen-Loève 展开拟最优，但常数形如 $e^{C/\rho}$，根源是切空间投影中的 $C^{-1}=(\mathbb E[\mathbf Y\mathbf Y^T])^{-1}$，在编号 15 自己的 $C^{\dagger}$ 变通中也看得见。编号 83 干脆不组装 $C$：模态逐项构造，系数更新除以的是标量 $l_{n+1}(\xi)$。代价则由编号 83 自己的对照实验坦白说出——先前模态永不回访，因此每加一项误差下降得比重解全部系数的方法慢，它买到的是稳定性与代价而不是速率。编号 83 同时丢掉了编号 15 的理论所依赖的两样东西（均值场 $\bar u_S$ 与规范条件 $\langle\partial_tU_i,U_j\rangle=0$），这既解释了为什么编号 15 的定理搬不过去，也解释了编号 83 为什么没有收敛定理。

编号 83 与编号 17 之间还有一条次要的联系：两者都把一个大的耦合问题换成一列「每次只加一个方向」的小空间增广（这里是时间上的贪心加基，那里是网格层次上的修正），并且都由不确定性量化驱动。三篇合起来把这条线索的选择空间划清了——**耦合演化换来对最佳秩 $S$ 逼近的拟最优性，解耦贪心换来不出现协方差求逆，一次性特征求解换来最优基但要付稠密矩阵的代价。**

## 本页原文

- E. Musharbash, F. Nobile, and T. Zhou, [_Error analysis of the dynamically orthogonal approximation of time dependent random PDEs_](https://doi.org/10.1137/140967787), SIAM J. Sci. Comput. 37(2) (2015), pp. A776-A810。
- H. Xie and T. Zhou, [_A multilevel finite element method for Fredholm integral eigenvalue problems_](https://doi.org/10.1016/j.jcp.2015.09.043), J. Comput. Phys. 303 (2015), pp. 173-184。
- L. Chen, Y. Chen, Q. Li, and T. Zhou, [_A dynamical variable-separation method for parameter-dependent dynamical systems_](https://doi.org/10.1137/24M168427X), SIAM J. Sci. Comput. 47(3) (2025), pp. A1783-A1808（预印本 [arXiv:2502.08464](https://arxiv.org/abs/2502.08464)）。
- 曲率引理与动力低秩背景取自 O. Koch and C. Lubich, _Dynamical low-rank approximation_, SIAM J. Matrix Anal. Appl. 29(2) (2007), pp. 434-454。
- 编号 17 所改造的多层修正框架见 Q. Lin and H. Xie, [_A multi-level correction scheme for eigenvalue problems_](https://doi.org/10.1090/s0025-5718-2014-02825-1), Math. Comp. 84(291) (2014), pp. 71-88（预印本 [arXiv:1107.0223](https://arxiv.org/abs/1107.0223)）。
