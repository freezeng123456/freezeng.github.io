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

> [!note] 本页覆盖
> 编号 **15**（_SIAM J. Sci. Comput._ 37(2), 2015）、**17**（_J. Comput. Phys._ 303, 2015）、**83**（_SIAM J. Sci. Comput._ 47(3), 2025）。编号 17 的期刊正文需订阅访问且无预印本，该节区分「本文可确认的内容」与「所改造的通用框架」；编号 15 与 83 的内容依据可公开获取的全文逐式核对。

![在低秩流形上直接演化随机解](assets/diagrams/tao-zhou-papers/zh/low-rank-dynamics.svg)

## 15：动态正交逼近的误差理论

### 为什么需要演化的基

时间依赖随机 PDE 通常用**固定**基展开：随机变量方向用广义多项式混沌，空间方向用本征正交分解模态。困难在于解对随机参数的依赖会随时间显著改变，因此固定基需要不断增加项数才能维持精度；自适应或贪心修补（时间依赖多项式混沌、推广的本征正交分解）只能部分缓解。

Sapsis 与 Lermusiaux 的动态正交逼近让空间基与随机系数**同时**演化，从而规避这一增长，但此前几乎没有误差理论：已有结果针对相近的多构型时间依赖 Hartree 方法与动力低秩矩阵情形，而非随机 PDE。本文的第一个目标就是在动态正交逼近与动力低秩逼近之间建立精确对应，从而把已有理论搬过来。

### 拟设与规范条件

在有界 $D\subset\mathbb R^d$（$1\le d\le3$）与完备概率空间上考虑

$$
\frac{\partial u(x,t,\omega)}{\partial t}=\mathcal L\bigl(u(x,t,\omega),\omega\bigr),
$$

拟设为

$$
u_S(x,t,\omega)=\bar u_S(x,t)+\sum_{i=1}^{S}U_i(x,t)\,Y_i(t,\omega).
$$

表示的唯一性由三条条件固定：

$$
\mathbb E[Y_i(t,\cdot)]=0,
\qquad
\langle U_i(\cdot,t),U_j(\cdot,t)\rangle=\delta_{ij},
\qquad
\Bigl\langle \frac{\partial U_i(\cdot,t)}{\partial t},\,U_j(\cdot,t)\Bigr\rangle=0 .
$$

第三条是**规范条件**：基的时间导数被要求与基自身的张成空间正交，也就是不允许基在自己的张成空间内旋转。这一条消去的正是分解的多余自由度。

### 演化方程与三处不同的角色

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
\Pi^{\perp}_{U}[v]=v-\sum_{i=1}^{S}\langle v,U_i\rangle U_i
$$

是到 $\mathcal U=\mathrm{span}\{U_1,\dots,U_S\}$ 的正交补的 $L^2(D)$ 投影。

三个方程的角色不同：均值方程就是平均后的 PDE；基方程由算子与随机系数的相关性中**被投影出张成空间的那部分**驱动，各模态之间仅通过 $C$ 耦合；系数方程是把中心化算子对每个基函数作检验得到的、以 $\omega$ 为参数的常微分方程。

### 切空间与 $C^{-1}$ 的来源

把中心化部分写成 $u_S^{*}=\mathbf U^{T}\mathbf Y$，它落在秩 $S$ 流形 $\mathcal M_S$ 上。切空间由 $\delta u_S^{*}=\mathbf U^{T}\delta\mathbf Y+\delta\mathbf U^{T}\mathbf Y$ 构成，在规范条件的类比 $\langle\delta\mathbf U,\mathbf U^{T}\rangle=\mathbf 0$ 与 $\mathbb E[\delta\mathbf Y]=\mathbf 0$ 下分解唯一。切空间上的正交投影为

$$
P_{u_S^{*}}(v)=\mathbf U^{T}\langle v^{*},\mathbf U^{T}\rangle
+\Bigl(\Pi^{\perp}_{U}\bigl\{\mathbb E[v^{*}\mathbf Y^{T}]\bigr\}\,C^{-1}\Bigr)^{T}\mathbf Y,
$$

其中 $C=\mathbb E[\mathbf Y\mathbf Y^{T}]$ 按秩 $S$ 的定义满秩。**这就是后来出现在误差界里的 $\rho^{-1}$ 的来源。** 互补投影可以因子化为 $P^{\perp}_{u_S^{*}}v=\Pi^{\perp}_{U}\otimes\Pi^{\perp}_{\tilde{\mathcal Y}}v$。

用 Dirac-Frenkel 形式写出来，动态正交逼近就是 Galerkin 条件

$$
\mathbb E\Bigl[\Bigl\langle \frac{\partial u_S}{\partial t}-\mathcal L(u_S),\,v\Bigr\rangle\Bigr]=0
\quad \forall v=\bar v+v^{*},\ v^{*}\in T_{u_S^{*}(t)}\mathcal M_S,
$$

等价地

$$
\frac{\partial u_S}{\partial t}=\mathbb E[\mathcal L(u_S)]
+P_{u_S^{*}(t)}\bigl(\mathcal L^{*}(u_S)\bigr).
$$

这个写法把「同时演化基与系数」这件事翻译成一句话：**把右端投影到流形的切空间上。**

### Theorem 4.1：拟最优性与它的代价

参照对象是每个时刻的最佳 $S$ 秩逼近，即截断的 Karhunen-Loève 展开

$$
z_S(x,\omega,t)=\bar u(x,t)+\sum_{i=1}^{S}\sqrt{\mu_i(t)}\,\gamma_i(t,\omega)\,Z_i(x,t),
$$

其中 $\{\mu_i,Z_i\}$ 是协方差算子的特征对。定理假设该最佳逼近连续可微存在于 $(H^2(D)\cap H^1_0(D))\otimes L^2(\Omega)$ 内，且其最小奇异值一致有下界

$$
\sigma\bigl(z_S(t)\bigr)\ \ge\ \rho\ >\ 0,\qquad \forall t\in[0,\bar t].
$$

则存在 $0<\hat t\le\bar t$，使动态正交解在 $u_S(0)=z_S(0)$ 下满足

$$
\|u_S(t)-z_S(t)\|_0^2+a_{\min}\int_0^t |u_S(\tau)-z_S(\tau)|_1^2\,\mathrm d\tau
\ \le\ 2\alpha\,e^{2\beta(t)}\int_0^t \|z_S(\tau)-u(\tau)\|_1^2\,\mathrm d\tau,
$$

$$
\beta(t)=4\rho^{-1}\int_0^t\Bigl(4\|\mathcal L^{*}(z_S)\|_0+\|\mathcal L^{*}(u)\|_0
+\|\mathcal L^{*}(u_S)\|_0+\|\dot z_S^{*}\|_0^2\Bigr)\mathrm d\tau,
\qquad
\alpha=\max\Bigl\{\frac{a_{\max}^2}{2a_{\min}},\ 4\rho^{-1}\Bigr\}.
$$

结论的读法有两层。**上层是好消息**：动态正交误差被最佳秩 $S$ 误差控制，方法是拟最优的。**下层是代价**：常数按 $e^{C/\rho}$ 增长，即随最小奇异值倒数指数增长。这就是曲率障碍的定量形式——流形 $\mathcal M_S$ 的曲率约为 $1/\sigma_S$，当第 $S$ 个奇异值塌缩时切空间投影变得病态，界失去意义。

### 一个具体的失效场景与一个实现细节

论文构造了一个特征值交叉的例子：确定性 Laplacian 加随机初值，前两个 Karhunen-Loève 特征值在某时刻 $t^{*}$ 交叉。此时最佳秩 1 误差与动态正交秩 1 误差都可以闭式算出，交叉前后由不同分支主导。这个例子说明上面的假设不是技术性的：秩不足时最小奇异值确实会塌缩。

实现层面还有一个容易忽略的问题。协方差矩阵 $C(t^n)$ 可能奇异或极度病态——例如任何带随机系数但**确定性初值**的系统在起始时刻都有 $C\equiv0$。直接用 Moore-Penrose 伪逆 $C^{\dagger}$ 会把「非活跃」的基函数置零，从而**阻止秩增加**。论文改为把基方程直接写成

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

与微分算子的特征问题相比，这里有一处结构性差别：积分算子 $\mathcal T$ 是**紧的、光滑化的**，特征值在 $0$ 处聚集（而不是趋于 $+\infty$），离散化给出**稠密**矩阵。因此既没有稀疏性可用，装配与求解的代价又随分辨协方差核所需的网格分辨率急剧上升。

先前工作从离散化一侧攻这个代价——小波 Galerkin、推广的快速多极、谱元、带张量结构的 Legendre-Galerkin、紧积分算子的多层增广方法——但特征求解本身仍必须在最细网格上完成。本文的贡献正是取消这一点：把特征值问题转化为一列**积分迭代**加上只在**最粗网格**上的特征求解，因此任何高效积分方案都可以插进来，总工作量与最细网格上的一次积分相当。

其所改造的多层修正框架（Lin 与 Xie）的通用一步是：给定最粗空间 $V_H$ 与嵌套层次 $V_H\subset V_{h_2}\subset\cdots$，先在细空间上解一个**源问题**

$$
a(\tilde u_{h_{k+1}},v)=\lambda_{h_k}\,b(u_{h_k},v),\qquad \forall v\in V_{h_{k+1}},
$$

再在小的增广空间 $V_{H,h_{k+1}}=V_H+\mathrm{span}\{\tilde u_{h_{k+1}}\}$ 内解特征值问题

$$
a(u_{h_{k+1}},v)=\lambda_{h_{k+1}}\,b(u_{h_{k+1}},v),\qquad \forall v\in V_{H,h_{k+1}} .
$$

第二步的意义在于：特征问题只在维数为 $\dim V_H+1$ 的空间里解，即最粗空间加一个方向，从不在细网格上解。对积分算子而言，第一步的「源问题」不是线性求解而是**一次求积**（即 $\mathcal T$ 的一次作用），这正是摘要所说的积分迭代。

> [!warning] 可核实范围
> 该文正文需订阅访问且无预印本。可确认的是问题设定、关键词（不确定性量化、Karhunen-Loève 展开、Fredholm 特征值问题、多重网格有限元）与摘要中的贡献表述。上面给出的通用修正步来自被改造的 Lin-Xie 框架；**该框架对积分算子的具体改造形式、常数与数值结果本站均未核实**，引用时应回查期刊版。

## 83：贪心地逐项加基，使每一步都解耦

### 三种降阶各自冻结了错的东西

论文把定位写得很清楚。广义多项式混沌写成 $u\approx\sum_i\zeta_i(\xi)g_i(x,t)$，**参数系数在时间上被冻结**、空间模态随时间变化，并且需要来自全模型的大量先验信息。本征正交分解写成 $u\approx\sum_i\zeta_i(t;\xi)g_i(x)$，**空间模态在时间与参数两个方向上都被冻结**。

本征正交分解能否成功取决于解流形 $\mathcal M=\{u(\cdot,t;\xi)\}$ 的 Kolmogorov 宽度是否快速衰减，而这**对一阶线性输运问题与双曲问题失效**——即 Kolmogorov 障碍。编号 15 分析的动态正交方法通过让两个因子都动来突破这一障碍，但它（以及 DyBO、对偶 DO、DBO 等变体）都对全部 $N$ 个模态求解一个**耦合**系统，而这正是协方差求逆及其条件数问题的来源。

本文的角度是：保留双向时间依赖的拟设，但用**贪心加基逐项构造**模态，使每一步解耦成两个标量系数子问题。

### 模型与仿射假设

$$
\frac{\partial u}{\partial t}(x,t;\xi)=F\bigl(u(x,t;\xi);\xi\bigr),
\qquad u(x,0;\xi)=\mu(x;\xi),
$$

结构性假设是 $F(u;\xi)=\mathcal C(\xi)+\mathcal A(u;\xi)+\mathcal H(u;\xi)$（$\mathcal A$ 线性、$\mathcal H$ 非线性）且**参数依赖是仿射的**：

$$
\mathcal C(\xi)=\sum_{i}\kappa^i_C(\xi)\mathcal C^i,
\qquad
\mathcal A(u;\xi)=\sum_{i}\kappa^i_A(\xi)\mathcal A^i(u),
\qquad
\mathcal H(u;\xi)=\sum_{i}\kappa^i_H(\xi)\mathcal H^i(u),
$$

并要求初值可分离 $\mu(x;\xi)=\sum_i p^i(\xi)q^i(x)$。**仿射性正是离线在线划分得以成立的前提**；论文指出若它不成立，可以先用变量分离方法构造一个精度损失可忽略的仿射逼近。

### 拟设：既没有均值场，也没有正交规范

$$
u(x,t;\xi)\approx u_N(x,t;\xi):=\sum_{i=1}^{N}\zeta_i(t;\xi)\,g_i(x,t),
$$

$\{\zeta_i\}$ 依赖参数、$\{g_i\}$ 不依赖参数，**两者都随时间变化**。与动态正交方法的 $u\approx\bar u(x,t)+\sum_i\zeta_i(t;\xi)g_i(x,t)$ 对比，后者保留一个统计均值场并要求 $\{g_i(\cdot,t)\}$ 在每个 $t$ 上正交归一。**本文把均值场与正交规范同时去掉**，唯一性改由顺序贪心构造提供。

### 贪心与两个解耦子问题

第一步任取 $\xi_1$，令 $g_1(x,t)$ 为方程在 $\xi=\xi_1$ 处的**全解**，再用 $g_1$ 作检验得到 $\zeta_1(t;\xi)$ 的参数依赖常微分方程。第 $k\ge2$ 步以 $e:=u-u_{k-1}$ 为误差，选点准则是

$$
\xi_k\in\arg\max_{\xi\in\Xi}\triangle_k(\xi),
$$

其中 $\triangle_k$ 取 $\|e\|_{L^2([0,T];V)}$ 本身（若代价允许）或下面的后验界；当 $\triangle_k(\xi_k)<\varepsilon$ 时停止，否则把 $\xi_k$ 从候选集中移除。把方程按误差改写为

$$
\Bigl\langle \frac{\partial(e+u_{k-1})}{\partial t},v\Bigr\rangle
=\bigl\langle F\bigl((e+u_{k-1});\xi\bigr),v\bigr\rangle,
\qquad \forall v\in V,
$$

则两个子问题是：$g_k(x,t)$ 取该式在**单个**参数 $\xi=\xi_k$ 处的解（与参数无关的偏微分方程），而令 $e=g_k\zeta_k$、$v=g_k$ 得到 $\zeta_k(t;\xi)$ 的标量常微分方程。

### 初值构造：贪心方法在 $t=0$ 处补回规范条件

论文称初值是「最关键的成分之一」，其构造用与当前模态的 $L^2$ 匹配。记 $g_{k,0}(x)=g_k(x,0)$，则 $g_1(x,0)=\mu(x;\xi_1)$ 且

$$
\zeta_{1,0}(\xi)=\sum_{i}\frac{\langle q^i,g_{1,0}\rangle}{\langle g_{1,0},g_{1,0}\rangle}p^i(\xi);
$$

对 $k\ge2$，初始误差 $e_0(x;\xi)=\mu(x;\xi)-\sum_{j<k}g_{j,0}(x)\zeta_{j,0}(\xi)$ 给出 $g_{k,0}=e_0(\cdot;\xi_k)$ 与

$$
\zeta_{k,0}(\xi)=\sum_{i}\frac{\langle q^i,g_{k,0}\rangle}{\langle g_{k,0},g_{k,0}\rangle}p^i(\xi)
-\sum_{j=1}^{k-1}\frac{\langle g_{j,0},g_{k,0}\rangle}{\langle g_{k,0},g_{k,0}\rangle}\zeta_{j,0}(\xi).
$$

**第二个和式是一个类 Gram-Schmidt 的修正项。** 这正是本方法在 $t=0$ 处隐式地、局部地补回动态正交方法用规范条件显式提供的那部分结构。一个边界情形：若 $\mu\equiv0$，则每一步 $g_{k,0}=0$ 且 $\zeta_{k,0}=0$。

### 线性情形：新模态由上一步的残差驱动

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

注意 $\langle(g_k)_t,g_k\rangle$ 这一项：在动态正交方法里它被规范条件**强制为零**，而这里它只是被一路带着。这句话把两条路线的差别落到了一个具体的项上。

### 关键之处：除数是标量而不是矩阵

把 $[0,T]$ 分成 $N_t$ 步、用一阶差分与后向 Euler，并代入仿射展开与 $u_{k-1}$ 的分离形式，得到闭式递推

$$
\zeta_{k,n+1}(\xi)=\frac{c_{n+1}\,\zeta_{k,n}(\xi)+s_{n+1}(\xi)}{l_{n+1}(\xi)},
\qquad n=0,\dots,N_t-1,
$$

$$
c_{n+1}=\frac{\langle g_{k,n+1},g_{k,n+1}\rangle}{\tau},
\qquad
l_{n+1}(\xi)=2c_{n+1}-\frac{\langle g_{k,n},g_{k,n+1}\rangle}{\tau}
-\sum_{i}\kappa^i_A(\xi)\bigl\langle\mathcal A^i(g_{k,n+1}),g_{k,n+1}\bigr\rangle,
$$

而 $s_{n+1}(\xi)$ 由仿射项与前 $k-1$ 个已算出的模态组装。**这就是本文对动态正交方法协方差求逆的回答：$l_{n+1}(\xi)$ 是一个标量而不是矩阵，因此全程不出现 $C^{-1}$ 或 $C^{\dagger}$。**

离线阶段只存与参数无关的标量内积（如 $\langle g_{k,n+1},g_{k,n+1}\rangle$、$\langle\mathcal A^i(g_{j,n+1}),g_{k,n+1}\rangle$ 等）；在线阶段对新参数 $\bar\xi$ 只需重解一个廉价的标量常微分方程，其系数关于参数仿射，因此**在线代价与原问题的空间离散无关**。

非线性情形（Burgers）中残差方程含 $\langle e_{n+1}\partial_xe_{n+1},v\rangle$，论文用半隐格式把它近似为 $e_n\partial_xe_{n+1}$ 以控制代价，并说明同一手法可搬到 Allen-Cahn。

### 后验估计与「没有收敛定理」

后验界基于 $F$ 在 $u$ 处的**局部对数 Lipschitz 常数**

$$
L_V[F](u):=\sup_{v\ne u}\frac{\langle v-u,\;F(v;\xi)-F(u;\xi)\rangle}{\|v-u\|_V^2},
$$

对误差方程 $\partial_te=F(u;\xi)-F(u_{k-1};\xi)+r_k$ 用 $e$ 作检验并配比较引理，得到 $\|e\|_V\le\delta_k(t;\xi)$，

$$
\delta_k(t;\xi)=\int_0^{t}\alpha(s;\xi)e^{\int_s^{t}\beta(\tau;\xi)\mathrm d\tau}\mathrm ds
+e^{\int_0^{t}\beta(\tau;\xi)\mathrm d\tau}\|e(\cdot,0;\xi)\|_V,
$$

其中 $\alpha=\|r_k\|_V$、$\beta=L_V[F](u_{k-1})$，选点准则取 $\triangle_k(\xi)=\int_0^T\delta_k^2\mathrm dt$。**对耗散型 $F$，对数 Lipschitz 常数可以为负**，此时指数因子起衰减而不是放大作用——这正是该估计在长时间区间上可用的原因。

> [!warning] 没有收敛定理
> 论文的结论部分明确把「在合理假设下对本方法作严格收敛分析」列为未来工作。严格的内容是上面的后验界与那条精确的代数递推；「相对许多已有低秩分离技术降低复杂度、提高效率」是摘要的表述，由数值实验而非定理支持。论文也自陈一处限制：存储与在线计算仍依赖时空离散的规模，如何降低这一依赖仍是开放问题。

### 两组最有信息量的数值对比

**与时间依赖约化基方法的正面比较**（一维反应扩散，参数依赖边界条件，$\xi\in[1,3]^4$，测试 $M=10^3$）给出一个坦率的权衡。对照方法在 $N=4$ 时误差 $9.19\times10^{-6}$、本方法 $1.46\times10^{-4}$，即**对照方法的误差下降更快**；论文自己给出的解释是对照方法在每次迭代中更新全部参数系数，而本方法不修改先前步已算出的时间参数基函数。但**在 $N=7$ 之后对照方法的误差快速上升**，归因于其参数系数线性系统的病态；同时本方法的在线时间小一个数量级以上（$N=7$ 时每样本 $7.20\times10^{-3}$ 秒对 $1.19\times10^{-1}$ 秒，而全阶有限元加后向 Euler 为每样本 $3.75$ 秒）。论文的总结是：对照方法在某些条件下精度更高，本方法在计算效率与误差稳定性上更好。

**与静态基变量分离方法的比较**（二维热方程）是全篇最干净的一处证据。随 $N$ 增大，静态基方法**发散**（$N=8$ 时误差 $2.18$、$N=10$ 时 $9.64\times10^{-1}$），而本方法单调下降（$N=8$ 时 $4.52\times10^{-5}$）。论文据此说明「对参数与空间两个变量都使用时间依赖基函数」的必要性。

Burgers 算例中基场呈清楚的幅度层级（$g_1\sim3.5\times10^{-1}$ 递降到 $g_9\sim3\times10^{-5}$），即首个场承载核心信息、末几个场承载细尺度信息。Allen-Cahn 算例的误差先随项数下降而后**趋于平台**，与其余算例不同。

> [!note] 题名差异
> 主页把该文列为 _A dynamical variable-separation method for dynamical systems with random input_，而出版版本、预印本与第三方记录均为 _A Dynamical Variable-Separation Method for Parameter-Dependent Dynamical Systems_。本站按出版版本记录；预印本为 [arXiv:2502.08464](https://arxiv.org/abs/2502.08464)。

## 三篇的关系

| 编号 | 演化对象                 | 模态之间的关系        | 主要风险                 |
| ---- | ------------------------ | --------------------- | ------------------------ |
| 15   | 空间基与随机系数         | 同时演化，经 $C$ 耦合 | 最小奇异值塌缩           |
| 17   | 不演化（一次性特征求解） | 不适用                | 稠密矩阵的装配与求解代价 |
| 83   | 逐项贪心加入的分离表示   | 每步解耦成两个方程    | 贪心序列不保证最优       |

编号 17 与编号 15 的关系值得指出：编号 15 的误差界以截断 Karhunen-Loève 逼近为参照，而编号 17 解决的正是**如何算出那个参照对象**。两篇一起给出了这条路线的完整代价结构：要么承担计算 Karhunen-Loève 基的稠密特征问题，要么让基自己演化并承担最小奇异值的风险。

编号 83 提供了第三种取法：既不预先算出最优基，也不让全部模态同时演化，而是逐项贪心地加入并让每步解耦。三篇合起来把这条线索的选择空间划清了——**耦合演化换来对最佳秩 $S$ 逼近的拟最优性，解耦贪心换来不出现协方差求逆，一次性特征求解换来最优基但要付稠密矩阵的代价。**

## 覆盖核对

| 内容                                           | 论文 | 覆盖状态                                  |
| ---------------------------------------------- | ---- | ----------------------------------------- |
| 拟设与三条规范条件                             | 15   | 形式、条件与规范条件的含义                |
| 三个演化方程与各自的角色                       | 15   | 均值、基、系数方程及投影算子              |
| 切空间投影与 $C^{-1}$ 的来源                   | 15   | 切空间、投影公式、Dirac-Frenkel 形式      |
| Theorem 4.1 与 $e^{C/\rho}$ 常数               | 15   | 假设、结论、两层读法与曲率解释            |
| 特征值交叉的失效场景                           | 15   | 例子构造与其说明的问题                    |
| 奇异协方差的处理                               | 15   | 伪逆的陷阱、改写形式、逐步对角化          |
| Fredholm 特征问题的结构差别                    | 17   | 紧性、谱聚集、稠密矩阵                    |
| 多层修正步与积分迭代                           | 17   | 通用两步、维数论证、积分与线性求解的差别  |
| 三种降阶各自冻结的对象                         | 83   | 多项式混沌、本征正交分解、Kolmogorov 障碍 |
| 模型、仿射假设与拟设                           | 83   | 结构分解、仿射依赖、去掉均值场与规范      |
| 贪心准则与两个解耦子问题                       | 83   | 选点准则、停机、两个方程的来源            |
| 初值构造与类 Gram-Schmidt 修正                 | 83   | 两个初值公式及其与规范条件的关系          |
| 残差驱动与被带着的 $\langle(g_k)_t,g_k\rangle$ | 83   | 线性情形的两个方程与差别落点              |
| 标量除数的闭式递推                             | 83   | 递推、$c_{n+1}$、$l_{n+1}$、无 $C^{-1}$   |
| 离线在线划分与在线代价                         | 83   | 所存标量、在线标量方程、与空间离散无关    |
| 后验界与对数 Lipschitz 常数                    | 83   | 定义、$\delta_k$、可为负及其后果          |
| 没有收敛定理与自陈限制                         | 83   | 结论部分的表述与存储依赖                  |
| 两组数值对比                                   | 83   | 约化基方法的权衡、静态基方法的发散        |

## 本页原文

- E. Musharbash, F. Nobile, and T. Zhou, [_Error analysis of the dynamically orthogonal approximation of time dependent random PDEs_](https://doi.org/10.1137/140967787), SIAM J. Sci. Comput. 37(2) (2015), pp. A776-A810。
- H. Xie and T. Zhou, [_A multilevel finite element method for Fredholm integral eigenvalue problems_](https://doi.org/10.1016/j.jcp.2015.09.043), J. Comput. Phys. 303 (2015), pp. 173-184。
- L. Chen, Y. Chen, Q. Li, and T. Zhou, [_A dynamical variable-separation method for parameter-dependent dynamical systems_](https://doi.org/10.1137/24M168427X), SIAM J. Sci. Comput. 47(3) (2025), pp. A1783-A1808（预印本 [arXiv:2502.08464](https://arxiv.org/abs/2502.08464)）。
