---
title: 多步格式
description: 编号 8、18、23、33、35、61、68：用多个未来层提高倒向方程的时间阶
lang: zh
translation: en/computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes
tags:
  - 论文笔记
  - 随机微分方程
  - 多步格式
---

> [!note] 本页覆盖
> 编号 **8**（_SIAM J. Sci. Comput._ 36(4), 2014）、**18**（_J. Sci. Comput._ 69(2), 2016）、**23**（_Numer. Math. Theor. Meth. Appl._ 10(2), 2017）、**33**（_SIAM J. Numer. Anal._ 56(4), 2018）、**35**（_J. Sci. Comput._ 79, 2019）、**61**（_J. Comput. Math._ 40, 2022）、**68**（_J. Sci. Comput._ 94:53, 2023）。多数为 SIAM 或 Springer 期刊且无预印本，因此若干技术细节只能从摘要与参考文献确认；相应处均已标注。

![多步格式如何提高倒向方程的时间精度](assets/diagrams/tao-zhou-papers/zh/fbsde-multistep.svg)

## 8：把高阶要求全部转移到倒向方向

### 论文提出的问题

耦合的 Markov 型正倒向随机微分方程为

$$
X_t=X_0+\int_0^t b(s,X_s,Y_s,Z_s)\,\mathrm ds+\int_0^t\sigma(s,X_s,Y_s,Z_s)\,\mathrm dW_s,
$$

$$
Y_t=\xi+\int_t^T f(s,X_s,Y_s,Z_s)\,\mathrm ds-\int_t^T Z_s\,\mathrm dW_s,
\qquad \xi=\varphi(X_T).
$$

论文对当时状况的描述是：多数已有格式是 Euler 型、收敛率 $1/2$；对解耦系统的高阶方法「依赖同时对正向随机微分方程与倒向方程使用高阶方案」，而正向的高阶方案「计算量大且往往难以实施」。耦合情形下正向系数依赖倒向未知量，因此「似乎不容易设计高阶且高效的数值格式」。

论文明确提出并肯定回答的问题是：**如果正向随机微分方程只用 Euler 方法求解，倒向方程还能不能达到高阶精度？**

### 允许冻结正向系数的定理

设 $X_s$ 解 $\mathrm dX_s=b\,\mathrm ds+\sigma\,\mathrm dW_s$，其生成元为

$$
\mathcal A_t^x g(t,x)=\lim_{s\downarrow t}\frac{\mathbb E_t^x[g(s,X_s)]-g(t,x)}{s-t},
\qquad
\mathcal L^0_{t,x}=\frac{\partial}{\partial t}+\sum_i b_i\frac{\partial}{\partial x_i}
+\frac12\sum_{i,j}(\sigma\sigma^{\top})_{i,j}\frac{\partial^2}{\partial x_i\partial x_j}.
$$

论文的关键定理是：若 $f\in C^{1,2}$ 且 $\mathbb E^{x_0}_{t_0}[\mathcal L^0_{t,X_t}f(t,X_t)]<\infty$，则

$$
\frac{\mathrm d\,\mathbb E^{x_0}_{t_0}[f(t,X_t)]}{\mathrm dt}
=\mathbb E^{x_0}_{t_0}\bigl[\mathcal A_t^{X_t}f(t,X_t)\bigr],
$$

并且**在 $t=t_0$ 处**该导数与把 $X_t$ 换成**任何**系数仅在左端点匹配的扩散 $\bar X_t$ 后的导数相同：

$$
\bar b(t_0,\bar X_{t_0};t_0,x_0)=b(t_0,x_0),
\qquad
\bar\sigma(t_0,\bar X_{t_0};t_0,x_0)=\sigma(t_0,x_0)
\quad\Longrightarrow\quad
\left.\frac{\mathrm d\,\mathbb E[f(t,X_t)]}{\mathrm dt}\right|_{t_0}
=\left.\frac{\mathrm d\,\mathbb E[f(t,\bar X_t)]}{\mathrm dt}\right|_{t_0}.
$$

这条结论就是把真实正向扩散换成**系数冻结的**（Euler）扩散的理论许可：论文明确指出可以简单取 $\bar b(s,\bar X_s;t_0,x_0)=b(t_0,x_0)$、$\bar\sigma(s,\bar X_s;t_0,x_0)=\sigma(t_0,x_0)$ 对所有 $s\in[t_0,t]$ 成立。

**这条定理的一般意义值得单独指出：** 当一个耦合系统的两个方向精度要求不同时，先确认低精度方向的误差是否真的进入目标量。这里的答案是不会——在一步导数逼近的层面上，正向扩散只需在左端点匹配。

### 导数逼近的权由矩条件确定

对 $u\in C_b^{k+1}$ 与节点 $t_0<t_1<\cdots<t_k$、$\Delta t_i=t_i-t_0$，权由

$$
\sum_{i=0}^{k}\alpha_{k,i}\frac{(\Delta t_i)^j}{j!}=\delta_{j1},
\qquad j=0,1,\dots,k
$$

确定。这组条件把「用 $k+1$ 个层的值逼近一阶导数」写成一个线性系统，与常微分方程的多步法系数条件同形。

### 阶数的上限来自根条件

多步格式的稳定性由一个特征多项式的根条件决定，论文观察到的稳定窗口是 $1\le k\le6$。这个窗口在编号 8 中是经验性的；它的理论解释要等到编号 47（见[[computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes|稳定性理论一页]]），那里的根条件正是「稳定性」这一半。

## 18、23、35：三种不同的加阶路线

- **18（带跳的 FBSDE 多步格式）** 把同一思路推广到带跳过程。跳的出现使倒向方程多一项补偿泊松测度积分，因此需要额外处理与跳相关的条件期望。
- **23（延迟校正方法）** 用延迟校正而不是多步插值加阶：先用低阶格式算一个近似解，再解一列校正方程逐次提高阶。与多步法的差别是它不需要额外的启动值，代价是每加一阶多一遍扫描。
- **35（二阶 FBSDE 的显式延迟校正）** 把延迟校正用到二阶（全非线性）情形，那里除 $Y,Z$ 外还有二阶过程 $\Gamma$（见[[computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control|二阶 FBSDE 一页]]）。

## 33 与 61：均场情形

均场倒向随机微分方程是 McKean-Vlasov 型的，生成元依赖解的**分布**：

$$
Y_t=\xi+\int_t^T \mathbb E'\bigl[f(s,X'_s,Y'_s,Z'_s,X_s,Y_s,Z_s)\bigr]\mathrm ds
-\int_t^T Z_s\,\mathrm dW_s,
$$

其中 $\mathbb E'$ 是对 $(X,Y,Z)$ 的独立副本 $(X',Y',Z')$ 取期望。

编号 33 的摘要指出这「似乎是**首次尝试为均场倒向随机微分方程设计高阶数值格式**」——此前的均场格式都是一阶。所用格式是一族**显式** $\theta$ 格式，「显式」意味着均场生成元在已知的未来层上求值，因此不需要对 $Y^n$ 做非线性求解。编号 61 给出均场正倒向系统的显式多步格式，把编号 8 的路线搬到均场情形。

> [!note] 可核实范围
> 编号 33 与 61 的正文未能获取（无预印本，出版社阻止自动访问）。上述问题设定与「首次高阶」「显式 $\theta$ 格式」的表述可从摘要确认；$\theta$ 的具体放置、分布或期望 $\mathbb E'$ 的离散方式（粒子系统、求积或嵌套期望）以及 $Z$ 的处理，本站未核实。

## 68：一旦稳定性成为枢纽，就设计格式去最大化它

编号 68 的位置很清楚。到 2022 年，这一组已有（一）一族高阶多步格式，其稳定性只通过经验观察到的根条件窗口（例如编号 8 的 $1\le k\le6$）刻画；（二）一个一般的稳定性、相容性与收敛性框架（编号 47），证明了均方意义下的 Lax 等价定理，但没有告诉人如何**构造**稳定性好的格式。本文把这个环闭合：把**强稳定保持**的设计思想——来自双曲守恒律数值中的 Gottlieb-Shu-Tadmor 与 Gottlieb-Ketcheson-Shu 技术，其中人们最大化使凸组合（收缩性）性质得以保留的类 CFL 系数——引入 FBSDE。

摘要的表述是：作者先对一般类型的 FBSDE 多步格式做全面分析，据此给出**关于系数的新充分条件**使相应格式稳定且具有一定的相容阶，再据此提出**构造高阶强稳定保持多步格式的实用方法**。论文附录给出一张优化系数表，覆盖**一阶到五阶**，适用于均匀时间划分。

被分析的格式模板与编号 47 统一的一族同形：对 $k$ 步方法与 $\mathbb E_n[\cdot]=\mathbb E[\cdot\mid\mathcal F_{t_n}]$，

$$
\sum_{i=0}^{k}\alpha_i\,\mathbb E_n\bigl[Y^{n+i}\bigr]
=\Delta t\sum_{i=0}^{k}\beta_i\,\mathbb E_n\bigl[f(t_{n+i},X^{n+i},Y^{n+i},Z^{n+i})\bigr],
$$

配以由 $\mathbb E_n[Y^{n+i}\Delta W]/\Delta t$ 项构成的 $Z^n$ 递推。强稳定保持的问题因此是：对哪些系数向量 $(\alpha_i),(\beta_i)$，格式可以改写成若干类后向 Euler 步的**凸组合**，从而基础步的任何单调性或收缩性被继承，并且步长系数尽可能大。这与 Lenferink 以及 Spiteri 与 Ruuth 的强稳定保持线性多步理论完全对应，那些工作也在本文参考文献中。

> [!note] 可核实范围
> 上述定位、摘要中的三步表述与附录中一至五阶系数表的**存在**均可确认。充分条件的确切形式、所用的稳定性函数、以及得到「最优」系数的优化问题，本站未核实；系数的数值也未核实。

## 七篇的关系

| 编号 | 加阶手段           | 问题类型       | 稳定性依据             |
| ---- | ------------------ | -------------- | ---------------------- |
| 8    | 多层 Lagrange 插值 | 耦合 FBSDE     | 经验根条件窗口 $k\le6$ |
| 18   | 多步 + 跳处理      | 带跳 FBSDE     | 同族                   |
| 23   | 延迟校正           | FBSDE          | 同族                   |
| 33   | 显式 $\theta$ 格式 | 均场 BSDE      | 论文自证               |
| 35   | 延迟校正           | 二阶 FBSDE     | 同族                   |
| 61   | 显式多步           | 均场 FBSDE     | 同族                   |
| 68   | 按稳定性反向设计   | 一般多步 FBSDE | 新充分条件 + 优化      |

这条线索的形状值得总结：**先构造格式（编号 8 至 35），再统一分析（编号 47），最后按分析结果反向设计（编号 68）。** 编号 68 之所以可能，正因为编号 47 已把稳定性确立为收敛的充要一半，因此「设计稳定的格式」成为一个有明确目标的优化问题，而不是试错。

## 覆盖核对

| 内容                           | 论文       | 覆盖状态                                     |
| ------------------------------ | ---------- | -------------------------------------------- |
| 耦合 FBSDE 与论文提出的问题    | 8          | 系统、此前状况、明确的问题表述               |
| 生成元定理与左端点匹配         | 8          | 生成元、定理、冻结系数的许可与其一般意义     |
| 导数逼近权的矩条件             | 8          | 线性条件组                                   |
| 根条件窗口                     | 8          | 经验窗口与其后续理论解释的位置               |
| 跳、延迟校正、二阶延迟校正     | 18, 23, 35 | 各自的加阶路线与代价                         |
| 均场 BSDE 的形式与「首次高阶」 | 33         | McKean-Vlasov 形式、显式含义（限定核实）     |
| 强稳定保持的设计思想与格式模板 | 68         | 来源、三步表述、模板、凸组合问题（限定核实） |

## 本页原文

- W. Zhao, Y. Fu, and T. Zhou, [_New kinds of high-order multistep schemes for coupled forward backward stochastic differential equations_](https://doi.org/10.1137/130941274), SIAM J. Sci. Comput. 36(4) (2014), pp. A1731-A1751（预印本 [arXiv:1310.5307](https://arxiv.org/abs/1310.5307)）。
- Y. Fu, W. Zhao, and T. Zhou, [_Multistep schemes for forward backward stochastic differential equations with jumps_](https://doi.org/10.1007/s10915-016-0212-y), J. Sci. Comput. 69(2) (2016), pp. 651-672。
- T. Tang, W. Zhao, and T. Zhou, [_Deferred correction methods for forward backward stochastic differential equations_](https://doi.org/10.4208/nmtma.2017.s02), Numer. Math. Theor. Meth. Appl. 10(2) (2017), pp. 222-242。
- Y. Sun, W. Zhao, and T. Zhou, [_Explicit theta-schemes for mean-field backward stochastic differential equations_](https://doi.org/10.1137/17M1161944), SIAM J. Numer. Anal. 56(4) (2018), pp. 2672-2697。
- J. Yang, W. Zhao, and T. Zhou, [_Explicit deferred correction methods for second-order forward backward stochastic differential equations_](https://doi.org/10.1007/s10915-018-00896-w), J. Sci. Comput. 79 (2019), pp. 1409-1432。
- Y. Sun, J. Yang, W. Zhao, and T. Zhou, [_An explicit multistep scheme for mean-field forward-backward stochastic differential equations_](https://doi.org/10.4208/jcm.2011-m2019-0205), J. Comput. Math. 40 (2022), pp. 519-543。
- S. Fang, W. Zhao, and T. Zhou, [_Strong stability preserving multistep schemes for forward backward stochastic differential equations_](https://doi.org/10.1007/s10915-023-02111-x), J. Sci. Comput. 94 (2023), 53。
