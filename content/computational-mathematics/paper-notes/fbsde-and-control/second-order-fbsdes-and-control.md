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

## 25：条件期望求值点的数量

多维情形下，条件期望是多维高斯积分，用张量积求积会使求值点数随维数指数增长。编号 25 用**谱稀疏网格**处理这一增长，其出版题名为 _Efficient spectral sparse grid approximations for solving multi-dimensional forward backward SDEs_。

这一篇的位置值得指出：这条线索里的时间精度问题由多步格式与延迟校正解决，而**空间**表示的代价问题由稀疏网格与后来的 Sinc 逼近（[[computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes|编号 63]]）解决。两个方向的改进是独立的。

## 26：把控制迭代组织成梯度投影

编号 26 处理随机最优控制问题。随机最大值原理把最优性条件写成一个正倒向系统（Pontryagin 系统），而控制变量通常带约束（例如取值在一个凸集内）。此时最自然的迭代是**梯度投影**：沿代价对控制的 Gâteaux 导数下降，再投影回可行集。

论文的技术焦点是让这一迭代高效：每次梯度求值需要解一次正倒向系统，因此迭代次数与每次迭代的代价都要控制。编号 **41（通过 FBSDE 求解随机最优控制的高精度格式）** 从另一侧改进同一问题：把高阶 FBSDE 格式用进控制问题，使每次正倒向求解本身更准。

**这两篇的分工是这条线索的一个典型模式**：一篇改进外层迭代，一篇改进内层求解。

## 50 与 51：这一族的边缘

- **50（数据驱动反馈控制的高效数值算法）** 处理的是反馈控制在数据驱动设定下的计算。与前面几篇的差别在于：控制律不再从已知模型的最优性条件推出，而要从数据中获得，因此滤波与控制耦合在一起。
- **51（动态非线性互补问题的 Gauss-Seidel 型方法）** 处理带互补约束的动态问题。互补约束使问题非光滑，因此梯度型迭代不再直接适用，需要 Gauss-Seidel 型的分块迭代。这一篇与 FBSDE 的联系较弱，但与编号 26 共享「带约束的动态最优化」这一主题。

> [!note] 覆盖进度
> 本页尚未对编号 16、19、25、26、41、50、51 做逐式核对。上述内容给出的是各篇的问题设定与在本专题中的位置，未展开各自的格式、定理与收敛阶。这些工作多数发表在无预印本的期刊上，其正文需订阅访问。

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
