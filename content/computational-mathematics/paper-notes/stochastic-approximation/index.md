---
title: 随机逼近与配点设计
description: 二十三项工作把不确定性量化中的构造问题统一成一个加权采样问题
lang: zh
translation: en/computational-mathematics/paper-notes/stochastic-approximation
tags:
  - 计算数学
  - 论文笔记
  - 不确定性量化
---

这个专题包含 23 项工作，是七个专题中最早开始的一组（2010 年起），也是后面许多方向的技术来源。它们回答的是不确定性量化中最基础的问题：**给定一个昂贵的参数到解映射，用尽可能少的正问题求解，稳定地构造它的多项式逼近。**

![配点设计的统一流程](assets/diagrams/tao-zhou-papers/zh/sampling-design.svg)

## 从侵入式到非侵入式

前五篇（2010–2012）处理的是**侵入式**方法：把解按广义多项式混沌展开，再做 Galerkin 投影。这条路线把随机 PDE 变成一个耦合的确定性系统，问题随之变成这个系统的代数性质与收敛性。

例如编号 1 研究的正是这个耦合系统的系数矩阵。设随机场仿射依赖参数 $\kappa(x,y)=\kappa_0(x)+\sum_{i}\kappa_i(x)y_i$，Galerkin 投影给出 $\partial_t v=\nabla\cdot(A\nabla_x v)+f$，其中

$$
a_{jk}=\sum_{i=0}^{N}\kappa_i(x)\,e_{ijk},
\qquad
e_{ijk}=\int y_i\,\Phi_j(y)\Phi_k(y)\rho(y)\,\mathrm dy .
$$

矩阵 $A(x)$ 是否严格对角占优，决定了对角/非对角分裂求解器（混合显隐时间推进、Jacobi 迭代、解耦预条件共轭梯度）是否稳定且以与网格无关的速率收敛。Xiu 与 Shen 只对**对称** Beta 密度与 Gamma 密度证明了对角占优，并把一般非对称情形列为公开问题。编号 1 通过把问题归结为 Jacobi 三项递推系数的不等式，给出了有条件的肯定回答：对 $\rho(y_i)=(1-y_i)^{\alpha}(1+y_i)^{\beta}$，只要 $|\alpha|\ge\tfrac12$ 且 $|\beta|\ge\tfrac12$，则

$$
a_{jj}\ \ge\ \kappa_{\min}+\sum_{k\ne j}|a_{jk}|,
\qquad 1\le j\le M,\ \forall x\in\Omega .
$$

论文没有断言对所有 $\alpha\ne\beta$ 都成立，阈值能否去掉在该文中未解决。

从编号 6、9 起路线转向**非侵入式**：只在若干参数点上求解，再用离散最小二乘或稀疏恢复拟合。这一转向把问题从「代数性质」换成了「采样设计」，而后者是本专题真正的技术核心。

## 四组精读

| 精读页                                                                                                                           | 论文                   | 技术核心                 |
| -------------------------------------------------------------------------------------------------------------------------------- | ---------------------- | ------------------------ |
| [[computational-mathematics/paper-notes/stochastic-approximation/stochastic-galerkin-and-collocation\|随机 Galerkin 与随机配点]] | 1, 2, 3, 4, 5, 7, 38   | 耦合系统的代数与收敛性   |
| [[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares\|离散最小二乘逼近]]                      | 6, 9, 11, 13, 14       | 配点网格与随机求积       |
| [[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning\|最优采样与预条件]]        | 22, 24, 28, 45, 54     | Christoffel 权与贪心选点 |
| [[computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce\|稀疏恢复与数据驱动混沌]]   | 10, 21, 29, 32, 36, 44 | ℓ1 恢复与经验测度        |

## 核心量：Christoffel 函数

这一族工作的中心对象是 Christoffel 函数。给定正交基 $\{\phi_n\}_{n=1}^{N}$，它的离散形式给出加权最小二乘的权

$$
w(y)=\frac{N}{\sum_{n=1}^{N}\phi_n^{2}(y)} ,
$$

而它的倒数给出采样密度。这个量的作用可以一句话说清：**它衡量多项式空间在一点上的「集中程度」，按它加权采样使设计矩阵各行的范数接近一致。** 行范数一致，加权 Gram 矩阵 $G=A^{\mathsf T}W^2A$ 才接近单位矩阵，最小二乘才稳定。

这一观察把采样预算从「按目标测度取 $M\gtrsim N^2$ 量级」改善到「按 Christoffel 加权密度取 $M\gtrsim N\log N$ 量级」，而贪心选点（近似 Fekete 点，用列主元 QR 分解实现）进一步把 $M$ 压到接近 $N$。三条路线共享同一条流水线，差别只在「从哪个密度取点」。

## 逐篇核心思想

### 随机 Galerkin 与随机配点

- **1** 回答 Xiu 与 Shen 关于 Galerkin 系数矩阵对角占优的公开问题，给出非对称 Beta 密度下的显式阈值条件。
- **2** 分析随机波速标量双曲方程的随机配点收敛性。困难在于波速 $c(y)$ 可变号，流入边界随之在区域两端切换，因此即使数据光滑，解**关于随机变量**的正则性也很低。论文的技术装置是带权数据假设，把密度除以波速，例如

  $$
  \int_0^T\!\!\int_{\Gamma^-}\frac{\rho(y)}{|c(y)|}\bigl(\partial_y u_L(t;y)\bigr)^2\,\mathrm dy\,\mathrm dt<\infty,
  $$

  这个 $1/|c|$ 权正是把波速与可容许边界数据耦合起来的地方。与 Gottlieb 和 Xiu 的 Galerkin 分析相比，本文的收敛率由**数据**的假设推出，而不是由展开系数的渐近衰减假设推出。

- **3** 处理带随机输入的椭圆界面问题的随机 Galerkin 方法。
- **4** 给出随机波速标量输运方程谱逼近的收敛性分析。
- **5** 用双正交多项式做随机双曲问题的 Galerkin 方法。
- **7** 处理带随机输入的延迟微分方程，延迟结构使解对参数的依赖具有额外的分段特征。
- **38** 给出带随机输入的 Maxwell 方程组的高效随机 Galerkin 方法。

### 离散最小二乘逼近

- **6** 直接研究最小二乘多项式逼近的设计点选择，这是本专题从侵入式转向采样设计的起点。
- **9** 提出一类新的配点网格用于多元离散最小二乘逼近。
- **11** 处理无界区域上带随机求值的离散最小二乘投影，并将其用于参数不确定性量化。无界区域使「按什么密度采样」这个问题从技术细节变成决定性因素。
- **13** 用随机求积构造加权离散最小二乘多项式逼近。
- **14** 在非结构多元网格上做随机配点，处理的是参数域不是张量积结构时的情形。

### 最优采样与预条件

- **22** 提出 Christoffel 函数加权最小二乘算法，是这一支的理论中心。
- **24** 给出用于多项式混沌稀疏逼近的推广采样与预条件方案，把采样密度与预条件矩阵作为一对来设计。
- **28** 提出加权近似 Fekete 点，用贪心选点把样本数压到接近多项式维数。
- **45** 是发表在 _SIAM Review_ 上的综述，把最小二乘多项式逼近的构造方法统一组织，是进入这一支文献的推荐入口。
- **54** 处理核插值的最优设计，把同一套思路从多项式空间搬到再生核空间。

### 稀疏恢复与数据驱动混沌

- **10** 研究稀疏插值与确定性插值点的设计。
- **21** 用随机求积做基于 ℓ1 极小化的随机配点。
- **29** 与 **32** 引入梯度增强的 ℓ1 恢复：一次正问题求解同时提供函数值与全部偏导数，测量行数因此放大 $(1+d)$ 倍而未知量个数不变。这把「每个昂贵样本的信息量」变成一个可以设计的对象。
- **36** 与 **44** 处理数据驱动的多项式混沌：当输入分布只通过样本给出、没有解析形式时，正交基必须从经验测度构造，采样设计也随之改变。

![把梯度信息与采样密度一起写进恢复问题](assets/diagrams/tao-zhou-papers/zh/sparse-recovery.svg)

> [!note] 覆盖进度
> 编号 1、2、3、4、5、6、7、9、10、11 已有逐式核对的精读内容。编号 13、14、21、22、24、28、29、32、36、38、44、45、54 的精读页仍在补充：本页对这些工作只给出可从题名、摘要与相邻论文交叉引用确认的核心思想。

## 与其他专题的关系

这个专题是后面几个方向的技术来源。[[computational-mathematics/paper-notes/bayesian-inference/index|贝叶斯专题]]中多精度多项式混沌代理用的加权最小二乘与权 $w_i=M/\sum_m\Phi_m^2(z_i)$，正是这里的 Christoffel 权；[[computational-mathematics/paper-notes/scientific-machine-learning/index|科学机器学习专题]]中「用一个可采样的密度决定配点放在哪里」的做法，与这里「用 Christoffel 加权密度取样本」是同一条思路在不同函数类上的实现。

## 本专题的原文

编号与题录见[[computational-mathematics/paper-notes/catalog|论文总目录]]，各篇具体出处列在对应精读页末尾。
