---
title: 鞅型深度学习
description: 编号 86、93、96、97、100、108：把「残差为零」换成「某个过程是鞅」
lang: zh
translation: en/computational-mathematics/paper-notes/fbsde-and-control/martingale-deep-learning
tags:
  - 论文笔记
  - 随机最优控制
  - 深度学习
---

> [!note] 本页覆盖
> 编号 **86**（_SIAM J. Sci. Comput._ 47(4), 2025）、**93**（_J. Comput. Phys._ 555, 2026）、**96**（投稿 _SIAM Rev._，[arXiv:2408.14395](https://arxiv.org/abs/2408.14395)）、**97**（投稿 SIAM/ASA JUQ）、**100**（投稿 _Numer. Math._）、**108**（预印本）。

![把方程残差改写成鞅性质](assets/diagrams/tao-zhou-papers/zh/martingale-training.svg)

## 86：不显式求 $\inf_u H$

### 三处约束

高维 Hamilton-Jacobi-Bellman 方程的求解需要在每个时空点求

$$
\inf_{\kappa\in U}H\bigl(t,x,\kappa,\partial_xv,\partial^2_{xx}v\bigr),
\qquad
H(t,x,\kappa,z,p)=\tfrac12\mathrm{Tr}\bigl(p\,\bar\sigma\bar\sigma^\top(t,x,\kappa)\bigr)
+z^\top\bar\mu(t,x,\kappa)+c(t,x,\kappa).
$$

当 $U\subset\mathbb R^m$ 高维、或 $H$ 没有闭式极小点时，这个内层极小化本身受维数灾困扰。而大多数已有深度 PDE 求解器要么需要显式的 $\inf_u H$，要么需要沿时间递推地训练网络。此外，标准 deep BSDE 架构依赖 Euler-Maruyama 的**强**收敛（阶 $1/2$），因为它用到 $X$ 的逐轨道性质。

本文一次去掉三处约束：不需要显式 $\inf_u H$、不需要时间递推训练、只需要**弱**收敛（阶 $1$）。

### 第一步：把逐点极小原理换成积分极小原理

沿一条**无控制**扩散 $X_t=X_0+\int_0^t\mu\,\mathrm ds+\int_0^t\sigma\,\mathrm dB_s$ 定义 Hamilton 过程与其累积：

$$
H^{u,v}_t:=H\bigl(t,X_t,u(t,X_t),\partial_xv(t,X_t),\partial^2_{xx}v(t,X_t)\bigr),
\qquad
\mathcal M^{u,v}_t:=v(t,X_t)+\int_0^tH^{u,v}_s\,\mathrm ds .
$$

论文的引理指出：在适当可积性下，最优控制可以从**积分**条件

$$
\int_0^T\mathbb E\bigl[H^{u,v}_t\bigr]\mathrm dt
=\inf_{\bar u\in\mathcal U_{\rm ad}}\int_0^T\mathbb E\bigl[H^{\bar u,v}_t\bigr]\mathrm dt
$$

中得到。证明只有一行：$\varepsilon_t:=H^{u,v}_t-\inf_\kappa H\ge0$ 逐点成立，而上式迫使 $\int_0^T\mathbb E[\varepsilon_t]\mathrm dt\le0$，故 $\varepsilon_t=0$ 几乎处处。

**这一步的收益说得很清楚**：极小化对象从「在每个 $(t,X_t)$ 上对 $U$ 求极小」换成「对泛函 $\bar u\mapsto\int_0^T\mathbb E[H^{\bar u,v}_t]\mathrm dt$ 求极小」，后者是一个可用 Monte Carlo 估计的二重积分，且 $t$ 与 $x$ 可独立采样，因此天然可并行。

### 第二步：鞅刻画

在可积性条件下，对 $(u,v)\in\mathcal U_{\rm ad}\times C^{1,2}$，

$$
(\partial_t+\mathcal L)v(t,X_t)=-H^{u,v}_t\ \ \text{a.e.}
\qquad\Longleftrightarrow\qquad
\mathcal M^{u,v}_t=\mathbb E\bigl[\mathcal M^{u,v}_T\mid\mathcal F_t\bigr].
$$

正向由 Itô 公式给出：代入方程后 $\mathcal M^{u,v}_t=v(0,X_0)+\int_0^t\partial_xv\,\sigma\,\mathrm dB_s$，是一个鞅。反向由鞅表示定理：$\mathcal M^{u,v}_t=\mathcal M^{u,v}_0+\int_0^tZ_s\mathrm dB_s$，与 Itô 展开比较后

$$
Q_t=\int_0^t\bigl\{(\partial_t+\mathcal L)v+H^{u,v}\bigr\}(s,X_s)\,\mathrm ds
=\int_0^t\bigl\{\partial_xv\,\sigma-Z_s\bigr\}\mathrm dB_s
$$

同时是有限变差过程与连续鞅且 $Q_0=0$，故 $Q\equiv0$。

**这就是本文的中心陈述：PDE 残差为零当且仅当该过程是鞅。**

一条几何假设值得注意：要求无控制扩散的支撑包含最优受控扩散的支撑（$\Gamma(X_t)\supset\Gamma(X^*_t)$）。论文给出两种安排方式：取 $X_0\sim N(x_0,rI_d)$，$r>0$ 为超参数；或取 $X=X^{u_0}$，$u_0$ 是一个预先的控制。

另一条评注同样重要：整个表述只用到 $X$ 的**期望与条件期望**，不用逐轨道性质。因此 Euler-Maruyama 只通过它的**弱**阶（$1$）进入，而 deep BSDE 型方法依赖**强**阶（$1/2$）。这是论文观察到时间方向一阶收敛的理论原因。

### 第三步：从鞅到对抗式 min-max

由 $X$ 的 Markov 性，$\mathbb E[\mathcal M^{u,v}_T\mid\mathcal F_t]=\mathbb E[\mathcal M^{u,v}_T\mid X_t]$。为避免直接计算条件期望，论文把它换成对检验函数族的**弱**条件：

$$
\sup_{\rho\in\mathcal T}\Bigl|\int_0^{T-\Delta t}\mathbb E
\Bigl[\rho(t,X_t)\bigl(\mathcal M^{u,v}_{t+\Delta t}-\mathcal M^{u,v}_t\bigr)\Bigr]\mathrm dt\Bigr|^2=0 .
$$

论证依据是塔性质：$\mathbb E[\rho(t,X_t)(\mathcal M_{t+\Delta t}-\mathcal M_t)]=\mathbb E[\rho(t,X_t)\mathbb E[(\mathcal M_{t+\Delta t}-\mathcal M_t)\mid X_t]]$，因此对**所有** $\rho$ 成立就迫使 $\mathbb E[\mathcal M_{t+\Delta t}-\mathcal M_t\mid X_t]=0$——这正是条件期望的投影性质。这个结构就是对抗式学习：$\rho$ 是判别器。

增广 Lagrange 函数与相应的 min-max 问题为

$$
L(u,v,\rho,\lambda)=\int_0^T\mathbb E\bigl[H^{u,v}_t\bigr]\mathrm dt
+\lambda\Bigl|\int_0^{T-\Delta t}\mathbb E\Bigl[\rho(t,X_t)
\bigl(\mathcal M^{u,v}_{t+\Delta t}-\mathcal M^{u,v}_t\bigr)\Bigr]\mathrm dt\Bigr|^2,
$$

$$
(u,v)=\lim_{\lambda\to+\infty}\
\arg\min_{(\bar u,\bar v)}\Bigl\{\sup_{\rho\in\mathcal T}L(\bar u,\bar v,\rho,\lambda)\Bigr\},
\qquad
\mathcal V=\{v\in C^{1,2}:v(T,x)=g(x)\} .
$$

参数化后由控制网络 $u_\alpha$、值网络 $v_\theta$ 与对抗网络 $\rho_\eta$ 三者训练。

### 两个构造性的架构约束

控制网络把取值范围硬约束在 $U=\prod_i[a_i,b_i]$ 内：

$$
u_\alpha(t,x)=a+\frac{b-a}{6}\,\mathrm{ReLU6}\bigl(\psi_\alpha(t,x)\bigr),
\qquad
\mathrm{ReLU6}(y)=\min\{\max\{0,y\},6\},
$$

一般 $U$ 则用一个到 $U$ 的距离罚项。值网络把终值条件硬编码：$v_\theta(T,x)=g(x)$，$t<T$ 时 $v_\theta(t,x)=\phi_\theta(t,x)$。

**这两处与本站[[computational-mathematics/paper-notes/scientific-machine-learning/index|科学机器学习专题]]的取向一致**：能由架构保证的约束就不放进罚项。这里唯一留在罚项里的是鞅条件，而它带一个 $\lambda\to\infty$ 的极限。

代码见 [sx-fang/MartNet](https://github.com/sx-fang/MartNet)。

## 93、96、100：同一框架的推广方向

三篇的共同出发点是把编号 86 的鞅型判据**在时间上局部化**。对 $t\mapsto v(t,X_t^s)$ 用 Itô 公式，得

$$
\mathcal M_t^s:=v(t,X_t^s)-v(s,X_s^s)+\int_s^tf\bigl(r,X_r^s,v(r,X_r^s)\bigr)\,\mathrm dr
=\int_s^tR(r,X_r^s;v)\,\mathrm dr+\int_s^t(\partial_xv)^{\top}\sigma\,\mathrm dB_r,
$$

其中

$$
R(t,x;v):=(\partial_t+\mathcal L)v(t,x)+f\bigl(t,x,v(t,x)\bigr)
$$

**恰是 PDE 残差**。取条件期望消掉 Itô 积分，再取 $t=s+h$，得 $\mathbb E[\mathcal M_{s+h}^s\mid\hat X_s]=h\,R(s,\hat X_s;v)+O(h^2)$，于是鞅型判据可写成

$$
\mathbb E\bigl[\mathcal M_{t+h}^t\,\big|\,\hat X_t\bigr]=0,
\qquad 0\le t\le T-h .
$$

**这个写法的要点是残差从不被显式计算，却被这个条件所刻画。** 但它也带来一条明确的限制，论文自己在注记中指出：该条件只保证残差在**引导过程所探索到的区域内**为零，因此引导过程 $\hat X$ 必须以高概率覆盖关心的区域。这是这一族方法的实际风险所在。

- **93（深度随机差分方法）** 面向高维**准**线性抛物方程 $\mathcal Dv=f(t,x,v)$，其中 $\mathcal D=\partial_t+\mu^{\top}(t,x,v)\partial_x+\frac12\mathrm{Tr}[\sigma\sigma^{\top}(t,x,v)\partial_{xx}]$——注意 $\mu,\sigma$ 可以依赖 $v$ 本身，比编号 86 的固定 $\mathcal L$ 更一般。名称中的「随机差分」指用随机差分替代导数，因此进一步降低对自动微分的依赖：论文引用的对比是 PINN 型方法需要用自动微分组装 $d\times d$ 的 Hessian，在 $d\ge10^4$ 时会耗尽内存。它同时补上了此前鞅型方法缺失的收敛率估计，并且刻意只用 Taylor 展开与初等矩恒等式重新推导整套机制，以降低对随机分析背景的要求。
- **96（极高维准线性方程与随机最优控制的鞅型深度学习）** 是投向 _SIAM Review_ 的综合性工作，把上面的局部化鞅判据系统化并推向极高维度，其中引导过程与系统过程的分离是关键的结构装置。
- **100（无导数局部化随机方法）** 面向**半**线性抛物方程，用解耦的 FBSDE 表示 $Y_t=u(t,X_t)$、$Z_t=\sigma^{\top}\nabla u(t,X_t)$。在全局 Lipschitz 与线性增长假设下，再加上一致非退化 $\xi^{\top}\sigma\sigma^{\top}\xi\ge\lambda\|\xi\|^2$，该 FBSDE 有唯一适应解。「局部化」指把全局问题分解为局部子问题以控制方差，而「无导数」指整套离散只用条件期望，不含导数项。

## 97 与 108：均场方向

- **97（DeepSPoC）** 用深度学习实现**混沌的顺序传播**。混沌传播是把均场极限与有限粒子系统联系起来的经典结果；「顺序」意味着不是同时模拟大量粒子，而是顺序地更新一个表示，从而降低内存需求。
- **108（高维均场博弈的深度策略迭代）** 用带再生式重构的深度策略迭代处理均场博弈。均场博弈在 Hamilton-Jacobi-Bellman 方程之外还耦合一个 Fokker-Planck 方程，因此需要同时演化值函数与分布。

> [!note] 覆盖进度
> 本页六篇均已按全文逐式核对：编号 86 含出版版本与预印本，编号 93、96、97、100、108 按各自的 arXiv 全文（编号 96 覆盖其三个改题版本）。上文对编号 97 与 108 只给出问题设定与定位，是取舍而非核实缺口；两者的格式细节未在本页展开。

## 这条路线的一般判断

鞅型方法改变的是**残差的检验方式**。在高维空间中「逐点检验残差」不可行，因为没有网格；而「检验一个过程是否是鞅」只需要沿模拟轨道的期望，因此维数不再直接进入代价。

编号 86 把这一转换做得最完整，可以概括成三次替换：

1. 逐点极小原理 → 积分极小原理（去掉 $(t,x)$ 上的维数灾）；
2. PDE 残差 → 鞅性质（去掉网格）；
3. 条件期望 → 对检验函数族的弱条件（去掉条件期望的计算）。

每一次替换都把一个难以直接处理的条件换成一个可以用采样验证的等价条件。**第三次替换尤其值得注意**：它把「计算条件期望」变成「训练一个判别网络」，代价是引入对抗训练的全部不稳定性，收益是完全避开了条件期望的估计。

## 覆盖核对

| 内容                     | 论文 | 覆盖状态                                      |
| ------------------------ | ---- | --------------------------------------------- |
| HJB 方程与 Hamilton 函数 | 86   | 状态方程、代价、值函数、$H$ 的显式形式        |
| 三处被去掉的约束         | 86   | 显式 $\inf_u H$、时间递推、强收敛依赖         |
| 积分极小原理与其证明     | 86   | Hamilton 过程、积分条件、一行证明、并行性收益 |
| 鞅刻画与双向证明         | 86   | Itô 方向、鞅表示方向、有限变差与鞅的同时性    |
| 几何假设与两种安排       | 86   | 支撑包含条件、随机初值或预先控制              |
| 弱阶而非强阶             | 86   | 只用期望的后果与观察到的一阶收敛              |
| 对抗式弱条件与塔性质     | 86   | 检验函数族、投影性质、判别器的角色            |
| 增广 Lagrange 与三个网络 | 86   | 目标函数、min-max、$\lambda\to\infty$         |
| 两处硬约束的架构         | 86   | ReLU6 值域约束、终值条件硬编码                |

## 本页原文

- W. Cai, S. Fang, and T. Zhou, [_SOC-MartNet: a martingale neural network for the Hamilton-Jacobi-Bellman equation without explicit inf H in stochastic optimal controls_](https://doi.org/10.1137/24M1681033), SIAM J. Sci. Comput. 47(4) (2025), pp. C795-C819（预印本 [arXiv:2405.03169](https://arxiv.org/abs/2405.03169)；代码 [sx-fang/MartNet](https://github.com/sx-fang/MartNet)）。
- W. Cai, S. Fang, and T. Zhou, [_Deep random difference method for high-dimensional quasilinear parabolic partial differential equations_](https://doi.org/10.1016/j.jcp.2026.114767), J. Comput. Phys. 555 (2026), 114767。
- W. Cai, S. Fang, W. Zhang, and T. Zhou, _Martingale deep learning for very high dimensional quasi-linear partial differential equations and stochastic optimal controls_, [arXiv:2408.14395](https://arxiv.org/abs/2408.14395)，投稿 SIAM Rev.
- K. Du, Y. Xie, T. Zhou, and Y. Zhou, _DeepSPoC: a deep learning based sequential propagation of chaos method_，投稿 SIAM/ASA J. Uncertain. Quantif.
- S. Fang, C. Sheng, B. Su, and T. Zhou, _A derivative-free localized stochastic method for very high dimensional semi-linear parabolic PDEs_，投稿 Numer. Math.
- S. Fang, S. Wang, Z. Wu, H. Zhang, and T. Zhou, _Deep policy iteration for high-dimensional mean-field games with regenerative reformulation_，预印本。
