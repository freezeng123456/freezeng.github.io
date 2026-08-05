---
title: 二阶 FBSDE 与随机控制
description: 编号 16、19、25、26、41、50、51：全非线性方程与控制迭代
lang: zh
translation: en/computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control
tags:
  - 论文笔记
  - 随机最优控制
  - 全非线性方程
---

> [!note] 本页覆盖
> 编号 **16**（_Commun. Comput. Phys._ 18, 2015）、**19**（_Commun. Comput. Phys._ 21, 2017）、**25**（_Discrete Contin. Dyn. Syst. Ser. B_ 22(9), 2017）、**26**（_SIAM J. Numer. Anal._ 55(6), 2017）、**41**（_Numer. Math. Theor. Meth. Appl._ 13, 2020）、**50**（_J. Sci. Comput._ 85(2), 2020）、**51**（_SIAM J. Control Optim._ 58(6), 2020）。

## 为什么需要「二阶」FBSDE

一阶正倒向随机微分方程通过非线性 Feynman-Kac 关系对应**半线性**抛物方程：$Y_t$ 给出解，$Z_t$ 给出梯度。要覆盖**全非线性**方程（Hessian 以非线性方式出现），需要在概率表示中再引入一个过程。

Cheridito、Soner、Touzi 与 Victoir 的二阶正倒向随机微分方程正是为此设计：除 $Y_t$ 与 $Z_t$ 外还有 $\Gamma_t$，对应解的 Hessian。因此离散化的对象从两个变为三个，而每一个都需要一个条件期望的逼近。

编号 **16（全非线性抛物方程的概率型高阶格式）** 与 **19（二阶 FBSDE 的高阶格式）** 处理的正是这一设定。编号 19 的题名把应用写明：随机最优控制。二者与[[computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|多步格式一页]]中的编号 8 属同一技术传统，差别在于多一个过程要处理。编号 **35** 则用延迟校正处理同一设定。

编号 19 处理的耦合二阶系统写出来是三个方程：

$$
\begin{cases}
X_t=x+\displaystyle\int_0^tb(s,\Theta_s)\,\mathrm ds+\int_0^t\sigma(s,\Theta_s)\,\mathrm dW_s,\\[4pt]
Y_t=g(X_T)+\displaystyle\int_t^Tf(s,\Theta_s)\,\mathrm ds-\int_t^TZ_s\,\mathrm dW_s,\\[4pt]
Z_t=Z_0+\displaystyle\int_0^tA_s\,\mathrm ds+\int_0^t\Gamma_s\,\mathrm dW_s,
\end{cases}
\qquad
\Theta_t=(X_t,Y_t,Z_t,A_t,\Gamma_t)\in\mathbb R^m\times\mathbb R\times\mathbb R^d\times\mathbb S^d,
$$

$\mathbb S^d$ 为 $d\times d$ 实对称矩阵。**第三个方程是二阶设定的全部新内容**：它把 $Z$ 本身当作一个 Itô 过程展开，其扩散系数 $\Gamma$ 正是对应 Hessian 的那个过程。所谓「解耦」即指 $b,\sigma$ 不依赖 $(Y,Z,A,\Gamma)$。

论文对当时状况的判断是：已有的高阶 FBSDE 格式只适用于低维，而已有的高维格式只有低阶（它举出一个 12 维耦合 FBSDE 的例子，收敛阶为 1；而编号 25 的稀疏网格格式达到 6 维、阶数至多 3）。它声明的缺口是：「据我们所知，对二阶 FBSDE 尚无高阶数值方法的相关研究。」

## 25：条件期望求值点的数量

多维情形下，条件期望是多维高斯积分，用张量积求积会使求值点数随维数指数增长。编号 25 用**谱稀疏网格**处理这一增长，其出版题名为 _Efficient spectral sparse grid approximations for solving multi-dimensional forward backward SDEs_。

这一篇的位置值得指出：这条线索里的时间精度问题由多步格式与延迟校正解决，而**空间**表示的代价问题由稀疏网格与后来的 Sinc 逼近（[[computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes|编号 63]]）解决。两个方向的改进是独立的。

## 26：把控制迭代组织成梯度投影

### 问题与既有路线的缺口

目标是解带约束的随机最优控制问题

$$
\min_{u\in K}J(u)=\mathbb E\Bigl[\int_0^T\bigl(h(x_t^u)+j(u(t))\bigr)\mathrm dt+k(x_T^u)\Bigr],
\qquad
\mathrm dx_t^u=b(x_t^u,u(t))\,\mathrm dt+\sigma(x_t^u,u(t))\,\mathrm dW_t .
$$

论文把已有数值路线分成四类：化为有限维随机规划、动态规划（即解 HJB 方程，「最广泛使用的方法之一」）、鞍度型方法，以及基于**随机最大值原理**的方法。它指出的缺口很具体：随机最大值原理在理论研究中是常用工具，但「在数值设定下并未被广泛使用」。本文要补的正是这一处。

一个容易忽略的设定：这里的控制空间 $U=L^2([0,T];\mathbb R)$ 是**确定性**平方可积控制，论文的理由是工程与金融应用中的事前规划需要这种形式；自适应（反馈）情形另在第 5 节处理。

### 不动点刻画

一阶最优性条件是变分不等式 $(J'(u^\ast),v-u^\ast)\ge0$ 对一切 $v\in K$ 成立。记投影 $P_K\omega=\arg\min_{u\in K}\|u-\omega\|$，它等价于 $(P_K\omega-\omega,v-P_K\omega)\ge0$。把两者对照，即得对任意 $\rho>0$ 成立的**不动点刻画**

$$
u^\ast=P_K\bigl(u^\ast-\rho J'(u^\ast)\bigr).
$$

于是算法就是梯度投影迭代 $u^{i+1,N}=P_{K_N}\bigl(u^{i,N}-\rho_iJ_N'(u^{i,N})\bigr)$，其中控制空间被离散为分段常数 $U_N=\{\sum_n\alpha_n\chi_{I_n^N}\}$，而 $J_N'$ 是 $J'$ 的数值逼近。**这个写法把「带约束的最优控制」变成「反复求梯度并投影」，因此全部计算困难集中到一件事上：怎么算 $J'$。**

### 梯度由一个伴随 BSDE 给出

直接求 $J'$ 需要 Gâteaux 导数 $Dx_t^u(v)$，代价高。论文的做法是引入伴随过程，对 $p_t^uDx_t^u(v)$ 用 Itô 公式，$Dx_t^u(v)$ 的项相消，剩下

$$
J'(u)\big|_t=\mathbb E\bigl[p_t^u\,b_u'(x_t^u,u(t))+q_t^u\,\sigma_u'(x_t^u,u(t))\bigr]+j'(u(t)),
$$

其中 $(p^u,q^u)$ 解伴随 BSDE

$$
-\mathrm dp_t^u=f(x_t^u,p_t^u,q_t^u,u(t))\,\mathrm dt-q_t^u\,\mathrm dW_t,
\qquad p_T^u=g(x_T^u),
$$

生成元为 $f(x,p,q,u)=h'(x)+p\,b_x'(x,u)+q\,\sigma_x'(x,u)$。**注意它在 $p,q$ 上是线性的**——这一点决定了后面的格式虽名义上隐式却无需迭代。

论文特意与一条早期路线作对比（其 Remark 1）：那里的伴随方程是一个**预期型**随机微分方程，要求解是向后适应而非经典的向前适应，而论文指出「这样的要求一般并不成立」，即其适定性不清楚；上面的 BSDE 则由标准理论直接适定。**这是一处方法论上的实质差别，不只是技术偏好。**

### 离散与误差平衡

正倒向系统按 Euler 离散，倒向方程在 $[t_n,t_{n+1}]$ 上积分后取条件期望，用**左端点矩形公式**；条件期望在 Euler 状态下化为高斯积分，用 $L$ 点 Gauss-Hermite 公式配线性插值求值。

收敛性有两层。外层：若 $J'$ 在最优点附近 Lipschitz 且一致单调，且 $\epsilon_N=\sup_i\|J'(u^{i,N})-J_N'(u^{i,N})\|\to0$，则迭代收敛；进一步若 $u^\ast$ 与 $J'(u^\ast)$ Lipschitz，则 $\epsilon_N\sim O(\Delta t)$ 给出 $\|u^\ast-u^{i,N}\|\sim O(\Delta t)$。内层则是伴随系统的误差

$$
\hat{\mathbb E}\bigl[(\mu_n)^2\bigr]+\Delta t\sum_{n=0}^{N-1}\hat{\mathbb E}\bigl[(\nu_n)^2\bigr]
=O\bigl((\Delta t)^2\bigr)+O\bigl((\Delta x)^4/(\Delta t)^2\bigr),
$$

$(\mu_n,\nu_n)$ 为 $(p_n,q_n)$ 的误差。**第二项的形式值得记住：空间插值误差以 $(\Delta x)^4/(\Delta t)^2$ 进入，因此 $\Delta x$ 与 $\Delta t$ 不能独立地取小**——要让两项平衡，空间网格必须随时间步一同细化，这是插值型条件期望逼近的典型代价。

编号 **41（通过 FBSDE 求解随机最优控制的高精度格式）** 从另一侧改进同一问题：把高阶 FBSDE 格式用进控制问题，使每次正倒向求解本身更准。**这两篇的分工是这条线索的一个典型模式**：一篇改进外层迭代，一篇改进内层求解。

## 50 与 51：这一族的边缘

- **50（数据驱动反馈控制的高效数值算法）** 的设定是控制与**观测**同时存在：状态满足 $\mathrm dX_t=b(t,X_t,u_t)\mathrm dt+\sigma(t,X_t,u_t)\mathrm dW_t$，代价为 $J(u)=\mathbb E[\int_0^Tf(t,X_t,u_t)\mathrm dt+h(X_T)]$，而另有一个观测过程 $\mathrm dM_t=g(X_t)\mathrm dt+\mathrm dB_t$。与前面几篇的差别就在这第三个方程：控制律不能只从已知模型的最优性条件推出，还要用观测去更新对状态的估计，因此**滤波与控制耦合在一起**。
- **51（动态非线性互补问题的 Gauss-Seidel 型方法）** 处理的问题是

  $$
  \dot x(t)=F(t,x(t),y(t)),
  \qquad
  0\le y(t)\ \perp\ G(t,x(t),y(t))\ge0,
  $$

  $x(t)\in\mathbb R^m$、$y(t)\in\mathbb R^n_+$，并含差分半仿射系统与动态线性互补问题两个子类。按后向 Euler 离散后，每个时间点上都要解一个耦合的非线性系统，而论文指出已有两条主流做法各有确定的缺陷：**直接消元法**只在线性情形有效，消去 $x_j$ 后得到的约化矩阵 $M_h=hN(I-hA)^{-1}B+M$ **即使 $M$ 是 P-矩阵也可能不是**，而且形成它需要 $n$ 次大规模线性求解；**半光滑 Newton 法**只有局部收敛性，且需要 Clarke 广义 Jacobi 矩阵，对大规模问题代价高且不便。Gauss-Seidel 型分块迭代正是针对这两处。这一篇与 FBSDE 的联系较弱，但与编号 26 共享「带约束的动态最优化」这一主题。

> [!note] 覆盖进度
> 编号 19、25、26、50、51 已按预印本或作者接受稿全文核对。编号 16 与 41 未能获取正文——前者出版社对 PDF 直接下载返回 403，后者的 PDF 端点不可达，且两篇均无预印本——因此这两篇只给出问题设定与在本专题中的位置，其格式、定理与收敛阶本站不报告。编号 16 的构造思路按它所引用的姊妹篇（编号 19）反推，并在此标明。

## 七篇的定位

| 编号 | 处理的对象               | 相对其他篇的位置            |
| ---- | ------------------------ | --------------------------- |
| 16   | 全非线性抛物方程         | 二阶设定的起点              |
| 19   | 二阶 FBSDE 与随机控制    | 把二阶设定接到控制问题      |
| 25   | 多维条件期归的求值代价   | 空间表示方向的改进          |
| 26   | 带约束的随机最优控制     | 外层迭代（梯度投影）        |
| 41   | 随机最优控制的高精度格式 | 内层求解（高阶 FBSDE 格式） |
| 50   | 数据驱动反馈控制         | 控制律来自数据而非模型      |
| 51   | 动态非线性互补问题       | 非光滑约束下的分块迭代      |

## 本页原文

- T. Kong, W. Zhao, and T. Zhou, [_Probabilistic high order numerical schemes for fully nonlinear parabolic PDEs_](https://doi.org/10.4208/cicp.240515.280815a), Commun. Comput. Phys. 18 (2015), pp. 1482-1503。
- T. Kong, W. Zhao, and T. Zhou, [_High order numerical schemes for second-order FBSDEs with applications to stochastic optimal control_](https://doi.org/10.4208/cicp.OA-2016-0056), Commun. Comput. Phys. 21 (2017), pp. 808-834。
- Y. Fu, W. Zhao, and T. Zhou, [_Efficient spectral sparse grid approximations for solving multi-dimensional forward backward SDEs_](https://doi.org/10.3934/dcdsb.2017174), Discrete Contin. Dyn. Syst. Ser. B 22(9) (2017), pp. 3439-3458。
- B. Gong, W. Liu, T. Tang, W. Zhao, and T. Zhou, [_An efficient gradient projection method for stochastic optimal control problems_](https://doi.org/10.1137/17M1123559), SIAM J. Numer. Anal. 55(6) (2017), pp. 2982-3005。
- Y. Fu, W. Zhao, and T. Zhou, [_Highly accurate numerical schemes for stochastic optimal control via FBSDEs_](https://doi.org/10.4208/nmtma.OA-2019-0137), Numer. Math. Theor. Meth. Appl. 13 (2020), pp. 296-319。
- R. Archibald, F. Bao, J. Yong, and T. Zhou, [_An efficient numerical algorithm for solving data driven feedback control problems_](https://doi.org/10.1007/s10915-020-01358-y), J. Sci. Comput. 85(2) (2020), 58。
- S. Wu, T. Zhou, and X. Chen, [_A Gauss-Seidel type method for dynamic nonlinear complementarity problems_](https://doi.org/10.1137/19M1268884), SIAM J. Control Optim. 58(6) (2020), pp. 3389-3412。
