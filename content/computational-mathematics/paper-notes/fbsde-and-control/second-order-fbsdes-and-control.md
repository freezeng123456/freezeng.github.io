---
title: 二阶 FBSDE 与随机控制
description: 编号 16、19、25、26、41、50、51：全非线性方程、稀疏网格与控制迭代
lang: zh
translation: en/computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control
tags:
  - 论文笔记
  - 随机最优控制
  - 全非线性方程
---

> [!note] 本页覆盖
> 编号 **16**（_Commun. Comput. Phys._ 18(5), 2015）、**19**（_Commun. Comput. Phys._ 21(3), 2017）、**25**（_Discrete Contin. Dyn. Syst. Ser. B_ 22(9), 2017）、**26**（_SIAM J. Numer. Anal._ 55(6), 2017）、**41**（_Numer. Math. Theor. Meth. Appl._ 13(2), 2020）、**50**（_J. Sci. Comput._ 85(2), 2020）、**51**（_SIAM J. Control Optim._ 58(6), 2020）。
>
> 其中编号 **19**（arXiv 预印本全文）、**25**（arXiv 预印本全文）、**26**（作者接受稿全文）、**50**（arXiv 预印本全文）、**51**（已出版正文）已逐式核对，本页给出它们的完整推导、定理与数值实验。编号 **16** 与 **41** 的正文未能获取——两篇均无预印本，Global Science Press 的 PDF 端点不可达或对直接下载返回 403——因此只写出摘要与参考文献可支撑的内容；两篇的格式细节按其姊妹篇反推之处，均在文中标明为反推。

## 为什么需要「二阶」FBSDE

一阶正倒向随机微分方程通过非线性 Feynman-Kac 关系对应**半线性**抛物方程：$Y_t$ 给出解，$Z_t$ 给出（$\sigma$ 加权的）梯度。二阶导数只能**线性**地进入——它来自 Itô 公式里的二次变差项，没有别的入口。因此要覆盖 Hessian 以非线性方式出现的**全非线性**方程，必须在概率表示里再放进一个过程。

Cheridito、Soner、Touzi 与 Victoir 的二阶倒向随机微分方程正是为此设计。其定义是：对 $\mathrm dX_t=\mu(X_t)\mathrm dt+\sigma(X_t)\mathrm dW_t$ 与取值于 $\mathbb R\times\mathbb R^d\times\mathbb S^d\times\mathbb R^d$（$\mathbb S^d$ 为 $d\times d$ 实对称矩阵）的四元组 $(Y,Z,\Gamma,A)$，

$$
\mathrm dY_t=f(t,X_t,Y_t,Z_t,\Gamma_t)\,\mathrm dt+Z_t^{\top}\circ\mathrm dX_t,
\qquad
\mathrm dZ_t=A_t\,\mathrm dt+\Gamma_t\,\mathrm dX_t,
\qquad
Y_T=g(X_T),
$$

其中 $Z_t^\top\circ\mathrm dX_t$ 是 Fisk-Stratonovich 积分，与 Itô 积分的关系为 $Z_t^\top\circ\mathrm dX_t=Z_t^\top\mathrm dX_t+\frac12\operatorname{Tr}[\Gamma_t\sigma\sigma^\top]\mathrm dt$。相应的偏微分方程是

$$
-v_t(t,x)+f\bigl(t,x,v(t,x),Dv(t,x),D^2v(t,x)\bigr)=0,\qquad v(T,x)=g(x),
$$

而表示为 $Y_t=v(t,X_t)$、$Z_t=Dv(t,X_t)$、$\Gamma_t=D^2v(t,X_t)$、$A_t=\mathcal L\,Dv(t,X_t)$，$\mathcal L\varphi=\varphi_t+\frac12\operatorname{Tr}[D^2\varphi\,\sigma\sigma^\top]$ 是去掉漂移项后的 Dynkin 算子。理论要求 $f$ 关于 $\Gamma$ **单调下降**（这就是椭圆性／抛物性），并需要 PDE 在粘性意义下的比较原理才能反过来由 2BSDE 唯一确定 $v$。

**这里有一处对数值方法极重要的设计自由，值得单独指出。** CSTV 在给出 PDE 之后特别注明：PDE 的形式**不依赖**于 $\mu$ 与 $\sigma$，所以原则上可以只取 $\mu\equiv0$、$\sigma\equiv I_d$；但保留选择 $\mu,\sigma$ 的自由「为 Monte Carlo 格式的设计提供了额外的灵活性」。**这句话正是编号 16 头条论断的理论依据**，也与编号 8 的定理 2（一侧导数只看 $\bar b,\bar\sigma$ 的左端点值）严丝合缝地对上：在半线性情形，正向 SDE 由 PDE 指定；在全非线性情形，它是一个可以自由挑选的设计参数。

编号 16、19、35 使用的写法与 CSTV 略有不同，跨文献对公式时必须注意这两处差别：

1. **驱动噪声。** CSTV 用 $\mathrm dX$ 驱动 $Z$，编号 19 用 $\mathrm dW$ 驱动 $Z$。因为 $\mathrm dX=\mu\,\mathrm dt+\sigma\,\mathrm dW$，两边的 $\Gamma$ 差一个 $\sigma$ 因子，$A$ 则吸收了漂移项。
2. **$\sigma$ 加权。** CSTV 取 $Z=Dv$、$\Gamma=D^2v$；编号 19 的表示定理取 $Z_t=(\nabla_xu\,\sigma)(t,X_t)$、$\Gamma_t=(\nabla_x(\nabla_xu\,\sigma)\sigma)(t,X_t)$、$A_t=(\mathcal L(\nabla_xu\,\sigma))(t,X_t)$。**这保持了一阶情形的同一套 $\sigma$ 加权约定，也正是编号 8 的多步模板能逐字搬过来的原因。**

## 16：为全非线性抛物方程造概率型高阶格式

### 直觉

全非线性抛物方程的概率求解此前只有低阶方法（Fahim-Touzi-Warin、Guo-Zhang-Zhuo、Tan 等），而编号 8 的多步机器只为一阶 FBSDE 造过。把两者接起来是自然的下一步，但接起来之后还多出一件半线性情形没有的东西：**正向 SDE 不再由方程指定**。半线性方程的 Feynman-Kac 表示里，$b$ 与 $\sigma$ 就是方程的系数；全非线性方程里它们是自由的。这一篇把这条自由当作可优化的设计参数来用。

### 问题设定

全非线性抛物方程的 Cauchy 问题，通过 CSTV 的结果化为一个 2FBSDE：

$$
\begin{cases}
\mathrm dX_t=b(t,X_t)\,\mathrm dt+\sigma(t,X_t)\,\mathrm dW_t,\\
-\mathrm dY_t=f(t,X_t,Y_t,Z_t,\Gamma_t)\,\mathrm dt-Z_t\,\mathrm dW_t,\\
\mathrm dZ_t=A_t\,\mathrm dt+\Gamma_t\,\mathrm dW_t,
\end{cases}
\qquad Y_T=g(X_T).
$$

**这一组方程本站是按其姊妹篇编号 19 的式 (2.9) 反推写出的**，本文正文中的实际编号与写法未核实。

### 可核实的结论

摘要逐字确认两点。其一，方法是把编号 8 的结果加以推广而得（原话为「by extending our previous results [W. Zhao, Y. Fu and T. Zhou, SIAM J. Sci. Comput., 36 (2014), pp. A1731-A1751]」）。其二，也是本文特别标出的贡献：「在我们的数值格式中，**可以自由选择相关联的正向随机微分方程**，而恰当的选择能显著降低计算复杂度。」

这条自由的来源上一节已经说清：CSTV 的 PDE 不依赖 $\mu,\sigma$，而编号 8 的定理 2 说一步导数逼近只看左端点系数。两者叠起来，就允许挑一个方便的辅助扩散。**一个具体的好处（本站的推理，非论文原话）是：若把系数取成常数，则每个网格点上可以复用同一组 Gauss-Hermite 节点，省掉逐点重算。**

> [!note] 可核实范围
> 出版社对 PDF 直接下载返回 403，`doc.global-sci.org` 不可达，且本文无预印本，因此本站只读到落地页（含题名、作者、卷期页、DOI、日期、摘要与完整 28 项参考文献）与 OpenAlex 记录。四条参考常微分方程、$\alpha_{k,i}$ 权、正向的 Euler 步、局部 Lagrange 插值与 Gauss-Hermite 求积都可在姊妹篇编号 19 中读到，但**本文实际印出的格式编号与写法、以及其空间离散是否与编号 19 不同，本站未核实**。定理方面本站未读到任何定理陈述；按该族的惯例，$k$ 步格式对 $(Y,Z,\Gamma,A)$ 应为 $k$ 阶，但这一点对本文而言未核实。数值方面，摘要说「给出了包括 HJB 方程在内的多个数值例子」，**「包括 HJB 方程」是可核实的**，但具体测试问题与观察到的阶未核实。

### 与其他论文的关系

编号 16 与 19 是近乎孪生的两篇：**16 从 PDE 出发（全非线性抛物方程 → 2FBSDE → 格式），19 从 2FBSDE 出发并接上随机控制应用。** 从参考文献可以看出这批工作是同时进行的：16 的参考文献把 19 列为「arXiv:1502.03206, 2015」、把 25 列为「submitted, 2015」、把 23 列为「submitted, 2015」。编号 35 后来把这里的多步模板换成显式延迟校正，处理同一个 2FBSDE 类。

## 19：把多步机器搬到二阶方程上，并接到控制问题

### 直觉

二阶设定里要离散的量从两个变成四个：$Y,Z,\Gamma,A$。多出来的两个不是附赠品，$\Gamma$ 对应 Hessian，是全非线性项唯一的入口。表面上看这意味着四倍的条件期望，代价难以承受。

但编号 8 的结构在这里几乎原样重现，而且重现得比预想的更整齐。关键是**同一个「乘 $\Delta W$ 再取条件期望」的操作在每一层都降一次阶**：对 $Y$ 用它得到 $Z$，对 $Z$ 用它得到 $\Gamma$；而不乘 $\Delta W$ 直接对 $Z$ 的参考方程求导则得到 $A$。于是四个量由四条同构的参考常微分方程给出，用的是同一组 $\alpha_{k,j}$ 权。更关键的是，在 $t=t_n$ 处 $\Delta W_{t_n,t_n}=0$ 让所有含 $\Delta W$ 的项消失，于是 $Z^n,\Gamma^n,A^n$ **三个都是显式的**，只有 $Y^n$ 隐式。**代价并没有翻两番，因为多出来的三个量各只需一次求和。**

### 问题设定

耦合二阶正倒向随机微分方程（论文式 (1.1)）为

$$
\begin{cases}
X_t=x+\displaystyle\int_0^tb(s,\Theta_s)\,\mathrm ds+\int_0^t\sigma(s,\Theta_s)\,\mathrm dW_s,\\[4pt]
Y_t=g(X_T)+\displaystyle\int_t^Tf(s,\Theta_s)\,\mathrm ds-\int_t^TZ_s\,\mathrm dW_s,\\[4pt]
Z_t=Z_0+\displaystyle\int_0^tA_s\,\mathrm ds+\int_0^t\Gamma_s\,\mathrm dW_s,
\end{cases}
\qquad
\Theta_t=(X_t,Y_t,Z_t,A_t,\Gamma_t)\in\mathbb R^m\times\mathbb R\times\mathbb R^d\times\mathbb S^d,
$$

其中 $b:[0,T]\times\mathbb R^m\times\mathbb R\times\mathbb R^d\times\mathbb S^d\to\mathbb R^m$、$\sigma:\cdots\to\mathbb R^{m\times d}$、$f:\cdots\to\mathbb R$、$g:\mathbb R^m\to\mathbb R$。**第三个方程是二阶设定的全部新内容**：它把 $Z$ 本身当作 Itô 过程展开，其扩散系数 $\Gamma$ 正是对应 Hessian 的那个过程。所谓「解耦」即 $b,\sigma$ 不依赖 $(Y_t,Z_t,A_t,\Gamma_t)$。

论文对当时状况的判断是：已有的高阶 FBSDE 格式只适用于低维，而已有的高维格式只有低阶——它举出一个 12 维耦合 FBSDE 的例子，收敛阶为 1；而编号 25 的稀疏网格格式达到 6 维、阶数至多 3。它声明的缺口是：「据我们所知，对二阶 FBSDE 尚无高阶数值方法的相关研究。」

### 推导

**第一步：表示定理（定理 2.2）。** 设 $u(t,x)$ 解全非线性方程

$$
\mathcal Lu+f\bigl(t,x,u,\nabla_xu\,\sigma,\ \nabla_x(\nabla_xu\,\sigma)\sigma\bigr)=0,\qquad u(T,x)=g(x),
$$

$(X_t,Y_t,Z_t,\Gamma_t,A_t)$ 解解耦的 2FBSDE，则

$$
Y_t=u(t,X_t),\quad
Z_t=(\nabla_xu\,\sigma)(t,X_t),\quad
\Gamma_t=\bigl(\nabla_x(\nabla_xu\,\sigma)\sigma\bigr)(t,X_t),\quad
A_t=\bigl(\mathcal L(\nabla_xu\,\sigma)\bigr)(t,X_t),
$$

$\mathcal L$ 为生成元。这是非线性 Feynman-Kac 公式的全非线性推广。

**第二步：四条参考常微分方程。** 条件化给出两个积分恒等式

$$
\mathbb E^x_{t_n}[Y_t]=\mathbb E^x_{t_n}[g(X_T)]+\int_t^T\mathbb E^x_{t_n}[f(s,\Theta_s)]\,\mathrm ds,
\qquad
\mathbb E^x_{t_n}[Z_t]=\mathbb E^x_{t_n}[Z_{t_n}]+\int_{t_n}^t\mathbb E^x_{t_n}[A_s]\,\mathrm ds,
$$

各自乘以 $\Delta W^{\top}_{t_n,t}$ 再取条件期望又给出两个。对 $t$ 求导得四条参考方程：

$$
\frac{\mathrm d\,\mathbb E^x_{t_n}[Y_t]}{\mathrm dt}=-\mathbb E^x_{t_n}[f(t,\Theta_t)],
\qquad
\frac{\mathrm d\,\mathbb E^x_{t_n}[Z_t]}{\mathrm dt}=\mathbb E^x_{t_n}[A_t],
$$

$$
\frac{\mathrm d\,\mathbb E^x_{t_n}\bigl[Y_t\Delta W^{\top}_{t_n,t}\bigr]}{\mathrm dt}
=-\mathbb E^x_{t_n}\bigl[f(t,\Theta_t)\Delta W^{\top}_{t_n,t}\bigr]+\mathbb E^x_{t_n}[Z_t],
\qquad
\frac{\mathrm d\,\mathbb E^x_{t_n}\bigl[Z^{\top}_t\Delta W^{\top}_{t_n,t}\bigr]}{\mathrm dt}
=\mathbb E^x_{t_n}\bigl[A^{\top}_t\Delta W^{\top}_{t_n,t}\bigr]+\mathbb E^x_{t_n}[\Gamma_t].
$$

**这四条的模式值得看清楚**：$(Y,\Delta W)$ 恢复 $Z$，$(Z,\Delta W)$ 恢复 $\Gamma$，而不带 $\Delta W$ 的 $Z$ 方程恢复 $A$。在 $t=t_n$ 处每个 $(\cdot)\Delta W^{\top}_{t_n,t_n}$ 项都消失，因此 $Z^n,\Gamma^n,A^n$ 显式，只有 $Y^n$ 隐式。

**第三步：半离散格式 2。** 取 Euler 选择 $\bar b(s,\cdot)=b(t_n,x)$、$\bar\sigma(s,\cdot)=\sigma(t_n,x)$ 对一切 $s\in[t_n,T]$：

$$
X^{n,j}=X^n+b(t_n,X^n)\Delta t_{n,j}+\sigma(t_n,X^n)\Delta W_{n,j},\qquad j=1,\dots,k,
$$

$$
Z^n=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}\Delta W^{\top}_{n,j}\bigr],
\qquad
A^n=\sum_{j=0}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Z^{n+j}\bigr],
$$

$$
\Gamma^n=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[(\bar Z^{n+j})^{\top}\Delta W^{\top}_{n,j}\bigr],
\qquad
-\alpha_{k,0}Y^n=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}\bigr]+f(t_n,X^n,Y^n,Z^n,\Gamma^n).
$$

注意 **$A^n$ 的求和从 $j=0$ 起**（它用的是不带 $\Delta W$ 的 $Z$ 导数模板，要用到当前层），其余三式从 $j=1$ 起。

**第四步：全离散格式 3。** 同上，但 $x\in\mathcal D^n_h$，插值 $\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}$，求积 $\hat{\mathbb E}^{n,x}$：

$$
X^{n,j}=x+b(t_n,x)\Delta t_{n,j}+\sigma(t_n,x)\Delta W_{n,j},
\qquad
Z^n=\sum_{j=1}^{k}\alpha_{k,j}\hat{\mathbb E}^{n,x}\bigl[\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}Y^{n+j}\Delta W^{\top}_{n,j}\bigr],
$$

$$
\Gamma^n=\sum_{j=1}^{k}\alpha_{k,j}\hat{\mathbb E}^{n,x}\bigl[\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}(Z^{n+j})^{\top}\Delta W^{\top}_{n,j}\bigr],
\qquad
A^n=\sum_{j=0}^{k}\alpha_{k,j}\hat{\mathbb E}^{n,x}\bigl[\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}Z^{n+j}\bigr],
$$

$$
-\alpha_{k,0}Y^n=\sum_{j=1}^{k}\alpha_{k,j}\hat{\mathbb E}^{n,x}\bigl[\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}Y^{n+j}\bigr]+f(t_n,x,Y^n,Z^n,\Gamma^n).
$$

每个 $(t_n,x)$ 上的顺序是：Euler 步 → 显式算 $Z^n,\Gamma^n,A^n$ → 隐式解 $Y^n$。**因为各网格点互不依赖，论文明确指出「整个过程可以完全并行」。** 隐式的 $Y^n$ 可用 Picard 迭代，也可在 $f$ 关于 $y$ 可微时用 Newton 迭代

$$
Y^{n,l+1}=Y^{n,l}-\frac{\alpha_{k,0}Y^{n,l}+\sum_{j=1}^k\alpha_{k,j}\hat{\mathbb E}^{n,x}\bigl[\mathcal I^{n+j}_{\mathcal D,\bar X^{n,j}}Y^{n+j}\bigr]+f(t_n,x,Y^{n,l},Z^n,\Gamma^n)}{\alpha_{k,0}+f_y(t_n,x,Y^{n,l},Z^n,\Gamma^n)} .
$$

**第五步：耦合情形与空间离散。** 格式 4 把 Euler 步换成 $X^{n,j}=x+b(t_n,x,Y^n,Z^n,\Gamma^n)\Delta t_{n,j}+\sigma(\cdots)\Delta W_{n,j}$，格式 5 再套一层 Picard 迭代，初值取 $(Y^{n,0},Z^{n,0},\Gamma^{n,0})=(Y^{n+1},Z^{n+1},\Gamma^{n+1})$；若 $b,\sigma$ 与 $(Y,Z,\Gamma)$ 无关则退化为格式 3。空间上用均匀张量网格 $\mathcal D_h=\{x_0+j\cdot h\}$、局部 Lagrange 插值、每维 10 点 Gauss-Hermite 求积，误差平衡取 $h=(\Delta t)^{(k+1)/(r+1)}$，$r$ 为插值次数。

### 定理

- **局部截断误差（式 (3.22)）。** 假设 $\mathbb E^x_{t_n}[Y_t]$、$\mathbb E^x_{t_n}[Y_t\Delta W^{\top}_{t_n,t}]$、$\mathbb E^x_{t_n}[Z_t]$、$\mathbb E^x_{t_n}[Z^{\top}_t\Delta W^{\top}_{t_n,t}]$ 及其至 $k+1$ 阶的 $t$ 导数有界，则
  $$
  \bar R^k_{y,n}=O((\Delta t)^k),\quad
  \bar R^k_{z,n}=O((\Delta t)^k),\quad
  \bar R^k_{A,n}=O((\Delta t)^k),\quad
  \bar R^k_{\Gamma,n}=O((\Delta t)^k).
  $$
  全离散的截断误差分成 **12** 项：对 $\bullet\in\{y,z,A,\Gamma\}$ 各有导数逼近项 $R^k_{\bullet,n}$、求积项 $R^{k,E}_{\bullet,n}$ 与插值项 $R^{k,I}_{\bullet,n}$。
- **稳定范围。** 论文第 5.1 节的原话是：方法「具有 $k$ 阶收敛率，且在 $1\le k\le6$ 时保持稳定，这与经典的数值常微分方程理论以及我们先前的结果一致」——与编号 8 是同一个 $k\le6$ 壁垒。

> [!warning] 这篇同样没有收敛定理
> arXiv 版本**没有**证明带显式常数的收敛定理；高阶性由上面的截断误差估计加数值支撑。已出版的 CiCP 版本是否补上定理，本站未核实。

### 数值实验

统一设定：$T=1$，均匀网格，每维 10 个 Gauss-Hermite 点，$>6$ 个 Lagrange 插值点，**四精度**（`real(16)`，34 位有效数字），FORTRAN 95 + OpenMP，Intel Xeon E5-2620 v2。四精度这一项值得注意——编号 8 在双精度下 $k=6$ 的阶被舍入吃掉，这里换成四精度正好绕开同一个陷阱。

**例 (5.1)，解耦 2FBSDE。** $\mathrm dX_t=\sin(t+X_t)\mathrm dt+c\cos(t+X_t)\mathrm dW_t$，$-\mathrm dY_t=\bigl(-\cos(t+X_t)\frac1cZ_t-\cos(t+X_t)(Y_t^2+Y_t)-\frac14\Gamma_t\bigr)\mathrm dt-Z_t\mathrm dW_t$，$\mathrm dZ_t=A_t\mathrm dt+\Gamma_t\mathrm dW_t$，$Y_T=\sin(T+X_T)$，$x=0.5$。精确解为

$$
Y_t=\sin(t+X_t),\quad
Z_t=c\cos^2(t+X_t),\quad
\Gamma_t=-2c^2\sin(t+X_t)\cos^2(t+X_t),
$$

$$
A_t=-c\sin(2t+2X_t)\bigl(1+\sin(t+X_t)\bigr)-c^3\cos(2t+2X_t)\cos^2(t+X_t).
$$

报出的结论是 $k$ 阶、在 $1\le k\le6$ 内稳定。效率对比是这一例最有说服力的一处：

| 格式  | 达到 $\lvert Y^0-Y_0\rvert=5.172\times10^{-5}$ 所需 $N$ | 耗时     |
| ----- | ------------------------------------------------------- | -------- |
| $k=1$ | 8192                                                    | 135.0 秒 |
| $k=2$ | 2048                                                    | 10.82 秒 |

**从一阶换到二阶，时间省了一个数量级还多**，而两者用的是同一份代码、同一个 Euler 正向步。

**例 4（表 5.4），四个分量的收敛率。**

| $K$ | $Y$  | $Z$  | $\Gamma$ | $A$  | 耗时      |
| --- | ---- | ---- | -------- | ---- | --------- |
| 1   | 1.00 | 0.99 | 0.97     | 0.97 | 21.54 秒  |
| 2   | 1.98 | 1.95 | 1.99     | 1.97 | 70.10 秒  |
| 3   | 2.98 | 2.98 | 3.11     | 2.99 | 259.11 秒 |

**这张表是这一篇最实质的证据**：不只是 $Y$ 与 $Z$，连二阶过程 $\Gamma$ 与 $A$ 都按 $k$ 阶收敛。它同时给出代价的斜率——从 $k=1$ 到 $k=3$，同一网格上的耗时约翻了 12 倍，符合 $k$ 步格式每层访问 $k$ 个未来层的结构。

**第 5.2 节的随机控制例子，「显微镜下追踪一个粒子」。** 状态 $\mathrm dX_t=\beta\alpha_t\mathrm dt+\sigma\mathrm dW_t$，代价 $J(\alpha)=\mathbb E\bigl[p\int_0^TX_t^2\mathrm dt+q\int_0^T\alpha^2\mathrm dt\bigr]$。HJB 方程为

$$
0=\partial_tV+\inf_{\alpha\in\mathbb R}\Bigl\{\tfrac{\sigma^2}{2}\partial^2_{xx}V+\beta\alpha\partial_xV+px^2+q\alpha^2\Bigr\}
=\partial_tV+\tfrac{\sigma^2}{2}\partial^2_{xx}V-\tfrac{\beta}{4q}(\partial_xV)^2+px^2,
\qquad V(T,x)=0,
$$

最优控制 $\alpha^\ast_t=-\frac{\beta}{2q}\partial_xV(t,x)$。所解的 2FBSDE 为

$$
\mathrm dX_t=\beta c\,\mathrm dt+\sigma\,\mathrm dW_t,
\qquad
-\mathrm dY_t=\Bigl(-\tfrac{\beta^2}{4q\sigma^2}Z_t^2-\tfrac{\beta c}{\sigma}Z_t+pX_t^2\Bigr)\mathrm dt-Z_t\,\mathrm dW_t,
\qquad
\mathrm dZ_t=A_t\,\mathrm dt+\Gamma_t\,\mathrm dW_t,
$$

控制由 $\alpha^n=-\frac{\beta}{2q\sigma}Z^n$ 恢复，参数 $\mu=0.1$、$r=0.03$、$\sigma=0.5$、$c=0.1$。论文的结论是「该方法对 2FBSDE 的解与最优控制 $\alpha$ 都具有高阶精度」；**本站确认这一论断存在，但未转录其数据表。**

论文还给出一条通用配方：对 $\mathrm dX_t=b(t,X_t,\alpha_t)\mathrm dt+\sigma(t,X_t,\alpha_t)\mathrm dW_t$ 与 $J(\alpha)=\mathbb E[\int_0^Tf(t,X_t,\alpha_t)\mathrm dt+g(X_t)]$，HJB 为 $\partial_tV+\inf_{\alpha\in U}\{\frac{\sigma^2}{2}\partial^2_{xx}V+b\,\partial_xV-f\}=0$，$\alpha^\ast(t,x)=\arg\inf_\alpha\{\cdots\}$；代入后化为 $\partial_tV+G(t,x,\partial_xV,\partial^2_{xx}V)=0$，$G(t,x,p,P)=\frac{\sigma(t,x,\alpha^\ast)^2}{2}P+b(t,x,\alpha^\ast)p-f(t,x,\alpha^\ast)$；解相应的 2FBSDE，再由 $\alpha^\ast_t=g(X_t,Y_t,Z_t,\Gamma_t)$ 取回控制。**这条配方的关键前提是 $\inf_\alpha$ 能解析地求出来**——这正是后来编号 86 与 96 的深度学习路线要绕开的东西。

**这些实验建立了什么，又差在哪里。** 建立的是：二阶设定下四个分量同时达到 $k$ 阶，这是本文的核心主张，而表 5.4 直接支持它；效率对比也表明高阶在实测中确实划算。差在三处。第一，例子仍是低维光滑问题，条件期望用张量 Gauss-Hermite；论文自己在引言里把「高阶只在低维、高维只有低阶」当作研究缺口，而这一篇解决的是前半句，没有解决后半句。第二，四精度算术是一项实质的实现选择，说明这一族在双精度下能否稳定跑到 $k=6$ 仍是问题（编号 8 在双精度下的 $k=6$ 就被舍入污染了）。第三，控制例子的高阶只有文字论断，本站未转录数据表，因此本页不报告其阶。

### 与其他论文的关系

编号 8 的直接推广（其摘要原话为「我们把多步格式加以非平凡的更新后推广至此」）。与编号 16 是近乎孪生的两篇（同作者、同年、同一套 2FBSDE 机器，只是 PDE 优先与 2FBSDE 优先的框架差别）。编号 35 后来把这里的多步模板换成显式延迟校正。**它走 SOC → HJB → 2FBSDE 的动态规划路线，与编号 26、41 走的 Pontryagin／伴随路线是同一问题的两条不同攻法**；编号 86 与 96 后来在高维上攻同一个 HJB 问题，用神经网络，而且明确避开这里解析求出的 $\inf_\alpha$。它的引言引用编号 25 说明稀疏网格达到 6 维、阶数至多 3，正是用来框定「低阶／高维」这一权衡的基准。

## 25：把维数的代价从指数压到多项式

### 直觉

编号 8 的时间离散没有维数问题：$\alpha_{k,j}$ 权与维数无关。有维数问题的是**空间**那一半——用均匀张量网格存函数、用张量 Gauss-Hermite 算条件期望、用局部 Lagrange 插值取值，三样东西的代价都随维数指数增长。论文的诊断说得很直白：这个组合「对高维 FBSDE 效率较低，因为所需计算量随维数指数增长」。

因此这一篇**完全不动时间离散**，只换空间那一半，换成 Smolyak 稀疏网格。稀疏网格的想法是不取张量积的全部节点，只取多重指标之和受限的那些子张量积；对具有混合光滑性的函数，这样丢掉的精度很少，而节点数从指数降为「多项式乘对数」。论文再加一层：用**分层基**把稀疏插值写成单个谱展开，使系数可以用逐维快速变换求出。

### 问题设定

时间离散与编号 8 逐字相同：同样两条参考常微分方程、同样的 $\alpha_{k,i}$ Vandermonde 系统、同样的局部生成元定理，以及耦合情形下的迭代半离散算法

$$
X^{n,j}=X^n+b(t_n,X^n,Y^{n,l},Z^{n,l})\Delta t_{n,j}+\sigma(t_n,X^n,Y^{n,l},Z^{n,l})\Delta W_{n,j},\quad j=1,\dots,k,
$$

$$
Z^{n,l+1}(X^n)=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}(\Delta W_{n,j})^{\top}\bigr],
\qquad
\alpha_{k,0}Y^{n,l+1}(X^n)=-\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}\bigr]-f(t_n,X^n,Y^{n,l+1},Z^{n,l+1}),
$$

迭代至 $\max\{|Y^{n,l+1}-Y^{n,l}|,|Z^{n,l+1}-Z^{n,l}|\}<\epsilon_0$。论文重申已知性质：「$k$ 步格式具有 $k$ 阶收敛率，前提是 $1\le k\le6$。」本文自述为「关于求解耦合 FBSDE 的多步格式系列论文的第二部分」。

### 推导

**第一步：稀疏网格。** 由一维序列 $\chi_i=\{x^i_0,\dots,x^i_{N_i-1}\}$（通常 $N_i=2^i+1$）构造 $q$ 维 Smolyak 稀疏网格

$$
\chi^p_q=\bigcup_{q\le|\mathbf i|_1\le p}\chi_{i_1}\otimes\chi_{i_2}\otimes\cdots\otimes\chi_{i_q},\qquad p\ge q .
$$

若序列**嵌套**（$\chi_1\subset\chi_2\subset\cdots$），则可写成 $\chi^p_q=\bigcup_{q\le|\mathbf i|_1\le p}\tilde\chi_{i_1}\otimes\cdots\otimes\tilde\chi_{i_q}$，其中 $\tilde\chi_1=\chi_1$、$\tilde\chi_i=\chi_i\setminus\chi_{i-1}$（$i>1$）。论文用两族一维节点，分工明确：

- **Chebyshev-Gauss-Lobatto（CGL）**，$C_i=\{x^i_j=\cos(j\pi/2^i),\ j=0,\dots,2^i\}$，**嵌套**，用于函数逼近与插值；
- **Gauss-Hermite（GH）**，$G_i=\{x^i_j,\ j=1,\dots,2^i-1\}$（$x^i_j$ 为 $2^i-1$ 次 Hermite 多项式的根），论文明确指出它**不嵌套**，用于条件期望的求积。

**这一处分工不是随意的：插值要靠嵌套才能写成分层增量，而求积不需要，所以可以继续用不嵌套但精度更高的 Gauss-Hermite。**

**第二步：稀疏插值。** 记 $\Delta_1=\mathcal I_1$、$\Delta_i=\mathcal I_i-\mathcal I_{i-1}$（$i>1$），Smolyak 插值算子为

$$
\mathcal I^p_q[f]=\sum_{q\le|\mathbf i|_1\le p}\Delta_{i_1}\otimes\cdots\otimes\Delta_{i_q}[f]
=\sum_{p-q<|\mathbf i|_1\le p}(-1)^{p-|\mathbf i|_1}\binom{q-1}{p-|\mathbf i|_1}\ \mathcal I_{i_1}\otimes\cdots\otimes\mathcal I_{i_q}[f].
$$

**第三步：「谱」在何处——分层基加快速变换。** 若基 $\{\tilde\phi_k\}$ 满足：对一切 $j\in I_i$ 与 $k\notin I_i$ 有 $\tilde\phi_k(x_j)=0$，则称之为**分层的**（定义 3.1）。在嵌套网格与分层基下，展开系数**与层级指标无关**，于是插值退化为单个谱展开

$$
\mathcal I^p_q[f](x)=\sum_{\mathbf k\in I^p_q}b_{\mathbf k}\,\tilde\phi_{\mathbf k}(x),
\qquad
\tilde\phi_{\mathbf k}=\prod_{i=1}^{q}\tilde\phi_{k_i},
$$

系数由 $f(x_{\mathbf j})=\sum_{\mathbf k\in I^p_q}b_{\mathbf k}\tilde\phi_{\mathbf k}(x_{\mathbf j})$（对一切 $\mathbf j\in I^p_q$）确定。这个系统用**逐维快速变换**（算法 2，`FastTran`）求解：沿一个维度施加一维逆矩阵 $T=(\tilde\phi_k(x_j))^{-1}_{k,j}$，在分层块上逐维推进——这就是摘要里说的 FFT 型加速。所用基是经变换的 Chebyshev 基 $T_k(x)$。

**第四步：稀疏网格 Gauss-Hermite 求积。** 条件期望里的张量 GH 规则换成一维 GH 求积算子 $Q_k$ 的 Smolyak 组合 $\sum\Delta_{k_1}\otimes\cdots\otimes\Delta_{k_q}$。

**第五步：整体算法。** 输入 $i=0,\dots,k-1$ 时定义在 $C^{p_i}_q[a_i,b_i]$ 上的 $Y^{N-i}(x),Z^{N-i}(x)$；每个时间层先跑 `FastTran` 得谱系数 $\{\beta^{N-i}_{\mathbf j}\},\{\gamma^{N-i}_{\mathbf j}\}$，再用稀疏网格 GH 求积代入多步公式。论文用一张表把两种做法对照：

| 方法            | 网格         | 条件期望         | 逼近与插值   |
| --------------- | ------------ | ---------------- | ------------ |
| SSG（本文）     | 稀疏网格     | 稀疏网格 GH 求积 | 稀疏网格插值 |
| LTG（即编号 8） | 张量均匀网格 | 张量 GH 求积     | Lagrange     |

### 定理

**本站未在文中找到本文自证的收敛定理。** 论文把编号 8 的 $k$ 阶结论直接搬过来，并引用文献中已知的稀疏网格求积误差结果（原话是它「得到如下关于高维立方体上函数的稀疏网格求积结果」，但**具体常数本站未转录**）。它的实质贡献是复杂度层面的：运行时间随维数「按某种多项式（而非指数）水平」增长——**这是一条数值观察（图 2），不是定理。**

### 数值实验

FORTRAN 95，16 颗 Intel Xeon E5620 CPU（2.40 GHz），3.0 GB 内存。

**例 1（二维耦合、周期）。** $b_i=\cos(4(x_i+t))/4-1$，$\sigma$ 为对角阵、$\sigma_{ii}=\cos(4(x_i+t))\sin(4(x_i+t))/4$，$f$ 取成使精确解为

$$
Y_t=\sin\bigl(4(X_{t,1}+t)\bigr)\sin\bigl(4(X_{t,2}+t)\bigr),
\qquad
Z_{t,i}=\Bigl(\prod_{k=1}^{2}\sin\bigl(4(X_{t,k}+t)\bigr)\Bigr)\cos^2\bigl(4(t+X_{t,i})\bigr).
$$

周期性使问题可以放在 $[-\pi,\pi]^2$ 上。SSG 空间用 CGL 稀疏网格 $C^7_2$、求积用 GH 稀疏网格 $G^3_2$；LTG 用 $\Delta x=(\Delta t)^{(k+1)/(n+1)}$ 平衡误差。$N=8,\dots,128$ 上的收敛率：

| 格式     | $E_Y$ 收敛率 | $E_Z$ 收敛率 |
| -------- | ------------ | ------------ |
| SSG 一步 | 0.982        | 0.986        |
| SSG 二步 | 1.987        | 1.955        |
| SSG 三步 | 2.632        | 2.955        |
| LTG 一步 | 1.174        | —            |

运行时间显示在可比精度下 SSG 明显更便宜。

**例 2（$q$ 维解耦）。** $b_i=\frac1qx_ie^{-x_i^2}$、$\sigma_{ii}=\frac1qe^{-x_i^2}$，精确解 $Y_t=\frac1q\sum_{j=1}^q\bigl(X_{t,j}^2\prod_{k\ne j}(X_{t,k}+t)\bigr)$ 与相应的 $Z_{t,i}$。对 $q=3,4,5,6$ 求解，结论是「所提多步格式即使对六维问题也具有高阶收敛率」，且运行时间随 $q$「按某种多项式（非指数）水平」增长。

**例 3（$q$ 维耦合）。** $b_i=\frac t2\cos^2(y+x_i)$、$\sigma_{ii}=\frac t2\sin^2(y+x_i)$。

**这些实验建立了什么，又差在哪里。** 建立的是：稀疏网格把同一套时间格式推到六维仍保持高阶，且时间增长曲线不再是指数——这正是本文要证明的那一件事。差在三处。第一，**六维是这条路线的实测天花板，不是它的下限**；后来的深度学习几篇（编号 86、93、96、97、100）跑到上千维，靠的是彻底放弃网格。第二，例 1 的三步格式在 $E_Y$ 上只观察到 2.632，明显低于名义的 3，而 $E_Z$ 是 2.955——**误差表里 $Y$ 与 $Z$ 的阶不同步，论文没有解释这一处，本页也不代为解释。** 第三，例 1 依赖周期性把定义域压到有界方体上，这是稀疏网格能用的前提；一般无界问题上如何截断，实验没有覆盖。

### 与其他论文的关系

自述为编号 8 之后的「系列第二部分」，逐字保留其时间离散，只换空间机器——**这条线索里时间与空间两个方向的改进是完全独立的，这一篇是最清楚的例证。** 编号 19 引用它，说明稀疏网格达到 6 维、阶数至多 3，用作「高阶／低维 对 低阶／高维」这一权衡的基准。编号 16 的参考文献把它列为「submitted, 2015」。编号 63 从另一条路解决同一问题（对 BSDE 用 Sinc 求积与插值），好处是**连空间插值都不需要**。而编号 86、93、96、97、100 完全放弃网格，跑到上千维；**编号 25 标记的是确定性网格路线的实用上限。**

## 26：把控制迭代组织成梯度投影

### 直觉

带约束的最优控制难在两处：约束与梯度。约束这一处有一个干净的处理——把「在凸集 $K$ 内极小化」写成一个不动点方程 $u^\ast=P_K(u^\ast-\rho J'(u^\ast))$，于是算法就是「求梯度、走一步、投影回去」。这一步把全部困难集中到梯度上。

而梯度这一处的困难是：直接求 $J'$ 要对状态过程求 Gâteaux 导数 $Dx^u_t(v)$，也就是每个方向 $v$ 都要解一个变分 SDE，代价不可接受。伴随方法的作用正是把这个方向依赖**一次性消掉**：引入伴随过程 $(p,q)$，对乘积 $p^u_tDx^u_t(v)$ 用 Itô 公式，所有含 $Dx^u_t(v)$ 的项互相抵消，剩下的表达式里 $v$ 只以显式的线性方式出现。于是一次伴随求解给出全部方向的梯度。

**这一篇还有一个容易被略过的技术红利：伴随 BSDE 的生成元在 $(p,q)$ 上是线性的。** 因此它的左端点矩形格式虽然名义上隐式，却只需一次移项就能解出，无需任何迭代——这与一般 BSDE 格式必须跑 Picard 迭代形成对比。

### 问题设定

要解的是带约束的随机最优控制问题

$$
\min_{u\in K}J(u)=\mathbb E\Bigl[\int_0^T\bigl(h(x_t^u)+j(u(t))\bigr)\mathrm dt+k(x_T^u)\Bigr],
\qquad
\mathrm dx_t^u=b(x_t^u,u(t))\,\mathrm dt+\sigma(x_t^u,u(t))\,\mathrm dW_t,\quad x|_{t=0}=x_0 .
$$

论文把已有数值路线分成四类：化为有限维随机规划；动态规划，即解 HJB 方程（「最广泛使用的数值方法之一」）；鞅型方法；以及基于**随机最大值原理**的方法。它指出的缺口很具体：随机最大值原理在理论研究中是常用工具，但「在数值设定下并未被广泛使用」。本文要补的正是这一处。

一个容易忽略的设定：这里的控制空间 $U=L^2([0,T];\mathbb R)$ 是**确定性**平方可积控制，论文的理由是工程与金融应用中的事前规划需要这种形式；自适应（反馈）情形 $U_{\mathbb F}=L^2_{\mathbb F}([0,T]\times\Omega;\mathbb R)$ 另在第 5 节处理。论文按一维写，并声明「整套框架容易推广到多维」。

### 推导

**第一步：投影给出不动点刻画。** 一阶最优性条件是变分不等式 $(J'(u^\ast),v-u^\ast)\ge0$ 对一切 $v\in K$ 成立，其中 $(J'(u),v)=\lim_{\rho\downarrow0}\frac{J(u+\rho v)-J(u)}{\rho}$。记投影 $P_K\omega=\arg\min_{u\in K}\|u-\omega\|$，它等价地由 $(P_K\omega-\omega,v-P_K\omega)\ge0$（对一切 $v\in K$）刻画。把 $\omega=u^\ast-\rho J'(u^\ast)$ 代入并与变分不等式对照，即得对任意 $\rho>0$ 成立的**不动点刻画**

$$
u^\ast=P_K\bigl(u^\ast-\rho J'(u^\ast)\bigr).
$$

**第二步：离散控制空间与迭代。** 取 $0=t_0<\cdots<t_N=T$、$\Delta t=T/N$、$I^N_n=[t_{n-1},t_n)$，令 $U_N=\{u\in U\mid u=\sum_{n=1}^N\alpha_n\chi_{I^N_n}\ \text{a.e.},\ \alpha_n\in\mathbb R\}$、$K_N=K\cap U_N$（凸闭）。则 $u^{\ast,N}=P_{K_N}(u^{\ast,N}-\rho J'(u^{\ast,N}))$，算法为**梯度投影迭代**

$$
u^{i+1,N}=P_{K_N}\bigl(u^{i,N}-\rho_iJ'_N(u^{i,N})\bigr),\qquad i=1,2,\dots,
$$

$J'_N$ 是 $J'$ 的数值逼近。

**第三步：伴随 BSDE 给出梯度。** 直接写出

$$
(J'(u),v)=\mathbb E\Bigl[\int_0^Th'(x^u_t)Dx^u_t(v)\,\mathrm dt+\int_0^Tj'(u(t))v(t)\,\mathrm dt+k'(x^u_T)Dx^u_T(v)\Bigr],
\qquad
Dx^u_t(v)=\lim_{\rho\downarrow0}\frac{x^{u+\rho v}_t-x^u_t}{\rho},
$$

其中 $Dx^u_t(v)$ 解变分 SDE $\mathrm dDx^u_t(v)=(b'_xDx^u_t(v)+b'_uv(t))\mathrm dt+(\sigma'_xDx^u_t(v)+\sigma'_uv(t))\mathrm dW_t$。为消去 $Dx^u_t(v)$，引入**伴随 BSDE**

$$
-\mathrm dp_t^u=f(x_t^u,p_t^u,q_t^u,u(t))\,\mathrm dt-q_t^u\,\mathrm dW_t,
\qquad p_T^u=g(x_T^u)=k'(x_T^u),
$$

$$
f(x,p,q,u)=h'(x)+p\,b_x'(x,u)+q\,\sigma_x'(x,u).
$$

对 $p_t^uDx_t^u(v)$ 用 Itô 公式，含 $Dx^u_t(v)$ 的项相消，得**梯度表示**

$$
J'(u)\big|_t=\mathbb E\bigl[p_t^u\,b_u'(x_t^u,u(t))+q_t^u\,\sigma_u'(x_t^u,u(t))\bigr]+j'(u(t)).
$$

**注意生成元在 $p,q$ 上是线性的**——这一点决定了后面的格式虽名义上隐式却无需迭代。

> [!warning] 与一条早期路线的实质差别
> 论文的 Remark 1 特意作了对比：其参考文献 [12] 中的伴随方程是一个**预期型**随机微分方程，要求解是向后适应而非经典的向前适应，而论文指出「这样的要求一般并不成立」，即其适定性不清楚；上面的 BSDE 则由标准理论直接适定。**这是一处方法论上的实质差别，不只是技术偏好。**

**第四步：FBSDE 与其 Feynman-Kac 形式。** 合起来的系统与其非线性 Feynman-Kac 表示为

$$
\begin{cases}
\mathrm dx^u_t=b(x^u_t,u(t))\,\mathrm dt+\sigma(x^u_t,u(t))\,\mathrm dW_t,\quad x|_{t=0}=x_0,\\
-\mathrm dp^u_t=f(x^u_t,p^u_t,q^u_t,u(t))\,\mathrm dt-q^u_t\,\mathrm dW_t,\quad p^u_T=g(x^u_T),
\end{cases}
$$

$$
p_t=\eta(t,x_t),\qquad q_t=\sigma(x_t,u(t))\,\partial_x\eta(t,x_t),
$$

$\eta$ 解 $\mathcal L^0\eta(t,x)=-f(x,\eta,\sigma\partial_x\eta,u(t))$、$\eta(T,x)=g(x)$，$\mathcal L^0\eta=\partial_t\eta+b\,\partial_x\eta+\frac12\sigma^2\partial_{xx}\eta$。

**第五步：Euler 格式。** 在 $[t_n,t_{n+1}]$ 上积分倒向方程、取 $\mathbb E^x_{t_n}[\cdot]$，并对生成元用**左端点矩形公式**：

$$
p^x_{t_n}=\mathbb E^x_{t_n}[p_{t_{n+1}}]+\Delta t\,f(x,p^x_{t_n},q^x_{t_n},u(t_n))+\bar R^x_{p,n},
\qquad
q^x_{t_n}=\frac{1}{\Delta t}\Bigl(\mathbb E^x_{t_n}[p_{t_{n+1}}\Delta W_{n+1}]+\bar R^x_{q,n}\Bigr),
$$

丢掉截断项得半离散格式

$$
p^x_n=\mathbb E^x_{t_n}[p_{n+1}]+\Delta t\,f(x,p^x_n,q^x_n,u(t_n)),
\qquad
q^x_n=\frac{1}{\Delta t}\mathbb E^x_{t_n}[p_{n+1}\Delta W_{n+1}].
$$

用 $\theta$ 格式族的语言，这是 $\theta_1=\theta_2=1$ 的成员：生成元在**当前**层 $t_n$ 上用未知量 $(p^x_n,q^x_n)$ 求值，$q^x_n$ 由单项 $\mathbb E_n[p_{n+1}\Delta W_{n+1}]/\Delta t$ 给出。**它形式上隐式，但因为伴随生成元在 $p$ 上线性，一次移项即得 $p^x_n$ 的闭式，无需任何迭代。**

**第六步：条件期望。** Euler 状态 $\tilde x^{t_n,x}_{t_{n+1}}=x+b(x,u(t_n))\Delta t+\sigma(x,u(t_n))\Delta W_{n+1}$ 使条件期望化为对 $\rho(\xi)=\frac{1}{\sqrt{2\pi}}e^{-\xi^2/2}$ 的高斯积分，用 $L$ 点 Gauss-Hermite 公式：

$$
\hat{\mathbb E}^x_{t_n}[p_{t_{n+1}}]=\sum_{\ell=1}^{L}\mathcal I_hp_{t_{n+1}}\bigl(x+b(x,u(t_n))\Delta t+\sigma(x,u(t_n))\sqrt{\Delta t}\,\xi_\ell\bigr)\,\omega_\ell,
$$

$$
\hat{\mathbb E}^x_{t_n}[p_{t_{n+1}}\Delta W_{n+1}]=\sum_{\ell=1}^{L}\mathcal I_hp_{t_{n+1}}\bigl(x+b(x,u(t_n))\Delta t+\sigma(x,u(t_n))\sqrt{\Delta t}\,\xi_\ell\bigr)\sqrt{\Delta t}\,\xi_\ell\,\omega_\ell,
$$

$\mathcal I_h$ 为到均匀空间网格 $\mathcal R_h=\{x_k\}$（$|k|\le P$）上的**线性插值**。全离散格式即把上式代入，对 $n=N-1,\dots,0$ 与每个 $x_k$ 求 $(p^k_n,q^k_n)$，再由 $J'_N(u)|_{t_n}=\hat{\mathbb E}[p_nb'_u+q_n\sigma'_u]+j'(u(t_n))$ 得梯度。算法 1 的循环是：给初值 $u^0\in U_N$ 与容差 $\epsilon_0$ →（一）置终值 $p^k_N=g(x_k)$ →（二）倒推解 $(p_n,q_n)$ →（三）算 $J'_N$ →（四）投影更新 $u$ → 重复。

**第七步：反馈控制（算法 2）。** 若允许控制取反馈形式 $\bar u(t,x)$ 且约束集 $K$ 是**逐点**的（在时间与空间上都逐点），则梯度里的期望消失：

$$
J'(u)^k_n=p^k_n\,b'_u(x_k,\bar u(t_n,x_k))+q^k_n\,\sigma'_u(x_k,\bar u(t_n,x_k))+j'(\bar u(t_n,x_k)),
$$

投影退化为逐点的标量投影 $\bar u^\ast(t_n,x)=P_C(\cdots)$。**关键后果是**：论文指出「我们不再需要 $t$ 之前的历史信息来计算 $J'(u)_t$，只需要 $t$ 这一时刻的信息」，因此整个算法可以在**单次倒向扫描**中完成，控制更新嵌在时间循环里，存储量大幅下降。

### 定理

- **定理 1（迭代收敛）。** 假设 $J'$ 在 $u^\ast$ 与 $u^{\ast,N}$ 附近 Lipschitz 且一致单调：存在 $c,C>0$ 使对一切 $v\in K$
  $$
  \|J'(u^\ast)-J'(v)\|\le C\|u^\ast-v\|,
  \qquad
  \bigl(J'(u^\ast)-J'(v),\,u^\ast-v\bigr)\ge c\|u^\ast-v\|^2,
  $$
  且把 $(u^\ast,K)$ 换成 $(u^{\ast,N},K_N)$ 后同样两式成立。再设 $\epsilon_N=\sup_i\|J'(u^{i,N})-J'_N(u^{i,N})\|\to0$（$N\to\infty$）。若步长 $\rho_i$ 满足
  $$
  0<1-2c\rho_i+(1+2C)\rho_i^2\le\delta^2\quad\text{对某个 }0<\delta<1,
  $$
  则 $\|u^\ast-u^{i,N}\|\to0$（$i,N\to\infty$）。
- **推论 1。** 在定理 1 的条件下，若另设 $u^\ast$ 与 $J'(u^\ast)$ 在 $U$ 中 Lipschitz，则 $\epsilon_N\sim O(\Delta t)$ 蕴含 $\|u^\ast-u^{i,N}\|\sim O(\Delta t)$（$i\to\infty$）。
- **定理 2（主误差估计）。** 在假设 1 与引理 1 至 3 的条件下，记 $(\mu_n,\nu_n)$ 为 $(p^k_n,q^k_n)$ 的误差，则
  $$
  \hat{\mathbb E}\bigl[(\mu_n)^2\bigr]+\Delta t\sum_{n=0}^{N-1}\hat{\mathbb E}\bigl[(\nu_n)^2\bigr]
  =O\bigl((\Delta t)^2\bigr)+O\bigl((\Delta x)^4/(\Delta t)^2\bigr),
  $$
  从而
  $$
  \epsilon_N=\sup_i\|J'(u^{N,i})-J'_N(u^{N,i})\|=O(\Delta t)+O\bigl((\Delta x)^2/\Delta t\bigr),
  $$
  **特别地取 $\Delta x=\Delta t$ 给出 $\epsilon_N=O(\Delta t)$**，再由推论 1 得 $\|u^\ast-u^{N,i}\|=O(\Delta t)$（$i\to\infty$）。也就是说，**整个算法——控制迭代、BSDE 求解器、求积与插值合在一起——是一阶的。** 配套结果：命题 2 要求 $m\ge2$、$L\ge2$ 与 $\Delta x=O(\sqrt{\Delta t})$；引理 1 假设 $b,\sigma\in C^{0,4}_b$；局部截断误差满足 $\tilde R^k_{q,n}=O((\Delta t)^2)$。

> [!warning] $\Delta x$ 与 $\Delta t$ 不能独立细化
> 误差估计的第二项 $O((\Delta x)^4/(\Delta t)^2)$ 的形式值得记住：**空间插值误差被 $(\Delta t)^2$ 除，因此固定 $\Delta t$ 一味细化空间网格并不会一直变好，反过来固定 $\Delta x$ 细化时间步会让这一项爆掉。** 两项要平衡就必须让空间网格随时间步一同细化，这是插值型条件期望逼近的典型代价。$\Delta x=\Delta t$ 的平衡条件是编号 8、19、25 中 $h=(\Delta t)^{(k+1)/(r+1)}$ 规则的低阶特例。

### 数值实验

步长取 $\rho_i=1/\sqrt i$，论文的解释是：$\rho$ 小有助于估计收敛，但太小会拖慢迭代。

**例 1。** $J(u)=\frac12\int_0^T\mathbb E[(x_t-x^\ast(t))^2]\mathrm dt+\frac12\int_0^Tu^2(t)\mathrm dt$，$K=U$，状态 $\mathrm dx_t=u(t)x_t\mathrm dt+\sigma x_t\mathrm dW_t$。设定 $x_0=1$、$T=1$、$\sigma=0.1$、$M=10^5$ 个样本、容差 $\epsilon_0=10^{-5}$，观察到一阶收敛。**本站只部分转录了该例的精确控制表达式，故不在此写出。** 另一组测试取 $\sigma=0.1$、$M=10^5$、$\epsilon_0=10^{-5}$、$N=40,50,\dots,100$，结论是「数值解与精确解吻合很好，观察到一阶收敛率」。

**例 3——反馈控制。** 状态与代价同前，但 $K$ 取为**随机**控制集。算法 1（确定性控制）与算法 2（反馈控制）所得 $J(u)$ 的对比是这一篇最直观的一张表：

| $N$ | 算法 1（确定性控制） | 算法 2（反馈控制） |
| --- | -------------------- | ------------------ |
| 100 | 0.84833              | 0.62535            |
| 200 | 0.84797              | 0.64507            |
| 400 | 0.84777              | 0.65509            |
| 800 | 0.84770              | 0.66013            |

论文的结论是：「使用反馈控制确实能改进结果（给出更小的目标泛函值），这是合理的，因为我们是在一个更大的控制集上做极小化。」**顺带注意两列的走向相反**：算法 1 的值随 $N$ 单调下降并稳定在 $0.8477$ 附近，算法 2 的值随 $N$ 单调上升到 $0.660$ 附近——两者都在收敛，只是从不同侧逼近各自的极限。

**例 4——投资组合问题。** 参考最优值为 $J(u)=15023$；取 $N=1000,2000,4000,8000$ 与 $M=N^2/10$，结论是「方法具有一阶收敛率」。

**这些实验建立了什么，又差在哪里。** 建立的是：定理 2 预言的一阶在四个例子上都被观察到，而反馈控制的对比表给出了一条独立于阶数的、关于**控制集大小**的结论。差在三处。第一，一阶是这一篇的设计目标而非局限的证据——它有意只用 Euler 格式，配一阶分析；真正的问题是能否更高阶，那是编号 41 的题目。第二，本站只部分转录了例 1 的精确解与例 2 的内容，因此本页给不出这两例的误差表。第三，全部例子都是一维的，而论文对多维只声明「框架容易推广」，没有给出多维实验。

### 与其他论文的关系

**这是编号 19 所走的动态规划路线的 Pontryagin／随机最大值原理对应物。** 这里的 FBSDE 是伴随／Hamilton 系统，不是值函数的 Feynman-Kac 表示。它有意只用一阶 Euler 格式并配以一阶分析；编号 41 是其直接续篇，把这一处升级为编号 8 的高阶多步格式，同时保留梯度投影外层。它自述在 Euler 格式上「紧随」赵卫东等人的 BSDE 格式文献，即编号 8 与 47 所属的那一支。编号 50 也走伴随／SMP 路线，但设定是部分观测的数据驱动反馈控制。

## 41：把内层求解器换成高阶的

### 直觉

编号 26 的整体精度是一阶，瓶颈很清楚：伴随 BSDE 用的是 Euler 格式，外层用的是固定步长梯度投影。两处都可以换。内层换成编号 8 那一族的高阶成员，外层换成利用曲率信息的拟 Newton 更新。

**这一篇的头条论断继承了编号 8 的核心现象**：即使状态方程仍用 Euler 离散，整体仍达到二阶。理由与编号 8 相同——控制问题里真正被逼近的是伴随量与梯度，它们只通过条件期望看见状态过程。

### 问题设定

与编号 26 相同的随机最优控制问题。本文把它「转化为一个等价的 FBSDE 随机最优性系统」，即：正向状态 SDE、倒向伴随 BSDE，以及控制方向上的变分不等式／稳定性条件三者合成的耦合系统。

### 可核实的结论

摘要逐字给出三点：（一）先把随机最优控制问题转化为「等价的 FBSDE 随机最优性系统」；（二）为所得系统设计「一个高效的**二阶 FBSDE 求解器**」与「一个**拟 Newton 型优化求解器**」；（三）「值得注意的是，**即使状态方程用 Euler 格式逼近，我们的方法仍具有二阶收敛率**」。关键词为「FBSDE；随机最优控制；随机最大值原理；投影拟 Newton 方法」，AMS 分类为 60H35、93E20、93E25、49M29、65C20、65K15。

参考文献的构成同样是可核实的证据：它引用编号 8、19、25、26、23，以及 Zhao-Zhang-Ju（SINUM 2010）、Zhao-Chen-Peng（SISC 2006）、彭实戈的随机最大值原理、Pontryagin 等人、雍炯敏与周迅宇、Powell-Yuan 的信赖域、戴彧虹与袁亚湘的共轭梯度、何炳生的投影收缩方法。**这份清单直接说明它的两个部件各自站在哪里：高阶 FBSDE 求解器来自多步／延迟校正工具箱，优化器来自非线性规划文献。** 它还引用杨杰与赵卫东的《Convergence of recent multistep schemes for a forward-backward stochastic differential equation》（EAJAM 2015），即编号 47 后来推广的那套收敛理论。

> [!note] 可核实范围
> `doc.global-sci.org` 不可达、`global-sci.org` 对直接下载返回 403，且本文无预印本，因此本站只读到落地页（含完整 46 项参考文献）与 OpenAlex 记录。具体而言：**二阶 FBSDE 求解器究竟是编号 8／47 族的 $k=2$ 成员、还是某个等价的 $\theta$ 或 Crank-Nicolson 型格式，本站未核实**；拟 Newton 更新的具体形式（BFGS？有限内存？）与全局化策略（线搜索？信赖域？）也未核实。二阶这一论断是作为定理证明的还是数值演示的、其假设与常数如何，本站同样未核实——从期刊与参考文献看有证明是合理的推测，但本站未读到。数值方面，摘要说「给出若干数值例子以说明所提格式的有效性与精度」，测试问题与观察到的阶未核实。

### 与其他论文的关系

**编号 26 的直接续篇**（并引用它）：同一问题、同一条 SMP／伴随 BSDE 路线，但从一阶升到二阶，从定步长梯度投影升到拟 Newton。它同时引用编号 8、19、25、23，即整套多步与延迟校正工具箱，那正是高阶 FBSDE 求解器的来源。**编号 26 与 41 的分工是这条线索的一个典型模式：一篇改进外层迭代，一篇改进内层求解。** 它与编号 19 的 HJB／动态规划路线、以及编号 86／96 的神经网络路线形成对照；值得一提的是，它的参考文献里已经列入韩劼群与鄂维南的《Deep learning approximation for stochastic control problems》与 Pereira 等人的《Learning deep stochastic optimal control policies using forward-backward SDEs》，说明 2019 年时作者已经注意到深度学习这条替代路线。

## 50：观测进来之后，滤波与控制就分不开了

### 直觉

前面几篇都默认状态可以完全观测：知道 $X_t$，就能按最优性条件算控制。一旦只能看到含噪的间接观测 $M_t$，控制必须对**观测滤链**可测，问题的结构就变了。经典的补救是**分离原理**——先滤波估计状态，再把估计当作真值做控制——但它严格成立只限于线性时不变系统配线性观测；一般非线性问题仍是开放的。

即便退而求其次地把滤波与控制拼起来，代价也难以承受：最大值原理路线要求每一步梯度都在整个状态空间上解一遍伴随 FBSDE，再配一个 Zakai 方程滤波器，一维以上就跑不动。

**这一篇的做法是把「精确解 FBSDE」这个目标整个放弃。** 它注意到 FBSDE 在这里只是**梯度的载体**，而随机梯度下降本来就只需要无偏的含噪梯度。于是每一步梯度只抽**一条**轨道：从粒子云里抽一个粒子作起点，抽一列高斯增量，倒推一次。两层期望（粒子与 Monte Carlo 路径）被压成一次抽样。**这样 $(Y,Z)$ 作为 $X$ 的函数的空间表示被彻底消掉了**——这是它能在实时反馈里用的根本原因，也是它明确不再声称解准 FBSDE 的原因。

### 问题设定

状态、代价与观测分别为

$$
\mathrm dX_t=b(t,X_t,u_t)\,\mathrm dt+\sigma(t,X_t,u_t)\,\mathrm dW_t,\quad X_0=\xi,
\qquad
J(u)=\mathbb E\Bigl[\int_0^Tf(t,X_t,u_t)\,\mathrm dt+h(X_T)\Bigr],
$$

$$
\mathrm dM_t=g(X_t)\,\mathrm dt+\mathrm dB_t,\qquad M_0=0,
$$

$X$ 取值于 $\mathbb R^d$，$u$ 取值于 $U\subseteq\mathbb R^m$，$g:\mathbb R^d\to\mathbb R^\ell$，$B$ 为与 $W$ 独立的 $\ell$ 维 Brown 运动。**控制同时进入漂移与扩散。** 记 $\mathbb F^M$ 为 $M$ 生成的（增广）滤链，容许集为 $\mathcal U_{\rm ad}[0,T]=\{u:\ \mathbb F^M\text{-循序可测、取值于 }U\}$，问题 (C\*) 为 $J^\ast(u^\ast)=\inf_{u^M\in\mathcal U_{\rm ad}[0,T]}J^\ast(u^M)$。作者称之为**数据驱动反馈控制**问题。他们还指出，在 Girsanov 变换 $\mathrm d\mathbb P^M=\Theta^T_t\,\mathrm d\mathbb P$（$\Theta^T_t=\exp(-\int_t^Tg(X_s)\mathrm dB_s-\int_t^T\frac12|g(X_s)|^2\mathrm ds)$）下观测 $M$ 变成标准 Brown 运动。

### 推导

**第一步：由随机最大值原理给出梯度。** 当 $u^\ast$ 是 $\mathcal U_{\rm ad}$ 的内点时，

$$
(J^\ast)'_u(u^\ast_t)=\mathbb E\Bigl[b_u(t,X^\ast_t,u^\ast_t)^\top Y_t+\sigma_u(t,X^\ast_t,u^\ast_t)^\top Z_t+f_u(t,X^\ast_t,u^\ast_t)^\top\ \Big|\ \mathcal F^M_t\Bigr],
$$

其中 $(Y,Z,\zeta)$ 解**伴随（Pontryagin）FBSDE 系统**

$$
\begin{cases}
\mathrm dX^\ast_t=b(t,X^\ast_t,u^\ast_t)\,\mathrm dt+\sigma(t,X^\ast_t,u^\ast_t)\,\mathrm dW_t, & X_0=\xi,\\[2pt]
\mathrm dM^\ast_t=g(X^\ast_t)\,\mathrm dt+\mathrm dB_t, & M_0=0,\\[2pt]
\mathrm dY_t=\bigl(-b_x^\top Y_t-\sigma_x^\top Z_t-f_x^\top\bigr)\mathrm dt+Z_t\,\mathrm dW_t+\zeta_t\,\mathrm dB_t, & Y_T=h_x(X^\ast_T)^\top,
\end{cases}
$$

$Z$ 是 $Y$ 对 $W$ 的鞅表示被积项，$\zeta$ 是对 $B$ 的。**注意倒向方程有两个被积项，因为噪声有两个来源。**

**第二步：条件梯度下降。** 朴素迭代是 $u^{l+1,M}_t=u^{l,M}_t-\rho\,(J^\ast)'_u(u^{l,M}_t)$。但在当前时刻 $t$，未来的观测 $\{\mathcal F^M_s\}_{s>t}$ 不可得，因此把控制投影到 $\mathcal F^M_t$ 上，记 $u^{l,M}_s|_t:=\mathbb E[u^{l+1,M}_s|\mathcal F^M_t]$，迭代

$$
u^{l+1,M}_s\big|_t=u^{l,M}_s\big|_t-\rho\,\mathbb E\bigl[(J^\ast)'_u\bigl(u^{l,M}_s|_t\bigr)\ \big|\ \mathcal F^M_t\bigr],
\qquad l=0,1,\dots,\ s\in[t,T].
$$

观测方程被**有意排除**在驱动 FBSDE 之外，因为 $M$ 不出现在梯度表达式里，而未来的数据也不可得。

**第三步：时间离散。** 记 $\mathbb E_i[\cdot]=\mathbb E[\cdot\mid\mathcal F^{X,B}_{t_i}]$、$\Delta t_i=t_{i+1}-t_i$、$\Delta W_{t_i}=W_{t_{i+1}}-W_{t_i}$，对 $i=N_T-1,\dots,n$：

$$
X_{i+1}=X_i+b(t_i,X_i,u^{l,M}_{t_i}|_{t_n})\Delta t_i+\sigma(t_i,X_i)\Delta W_{t_i},
$$

$$
Y_i=\mathbb E_i[Y_{i+1}]+\mathbb E_i\Bigl[b_x^\top Y_{i+1}+\sigma_x^\top Z_{i+1}+f_x^\top\Bigr]\Delta t_i,
\qquad
Z_i=\mathbb E_i\bigl[Y_{i+1}\Delta W_{t_i}\bigr]\cdot(\Delta t_i)^{-1},
$$

其中 $b_x,\sigma_x,f_x$ 在 $t_{i+1}$ 层求值。也就是说：正向用**左端点（Euler-Maruyama）**规则；$Y$ 在条件化后对漂移积分用**右端点**规则，因此倒向一步是**显式**的（$Y_{i+1},Z_{i+1}$ 已知）；$Z$ 用左端点规则与经典的 $\mathbb E_i[Y_{i+1}\Delta W_{t_i}]/\Delta t_i$ 表示。**Remark 3.1 指出 $\zeta$ 既不出现在格式里也不出现在 $(J^\ast)'_u$ 里，因此根本不需要为它设计格式。**

**第四步：条件期望用 Monte Carlo，不用求积。** 与编号 8、19、25、26 不同，这里用纯 Monte Carlo，取 $K$ 个标准高斯样本 $\{\omega^k_i\}$、$\Delta W_{t_i}\approx\sqrt{\Delta t_i}\,\omega^k_i$：

$$
Y_i=\frac1K\sum_{k=1}^KY^k_{i+1}+\frac{\Delta t_i}{K}\sum_{k=1}^K\Bigl[b_x(t_{i+1},X^k_{i+1},\cdot)^\top Y^k_{i+1}+\sigma_x(t_{i+1},X^k_{i+1})^\top Z^k_{i+1}+f_x(t_{i+1},X^k_{i+1},\cdot)^\top\Bigr],
$$

$$
Z_{t_i}=\frac{1}{\Delta t_i}\sum_{k=1}^K\frac{Y^k_{i+1}\sqrt{\Delta t_i}\,\omega^k_i}{K}.
$$

作者给出的理由是在高维中有效。

**第五步：粒子滤波与「单实现」技巧。** 条件分布 $p(X_t|\mathcal F^M_t)$ 用粒子云 $\{x^{(s)}_n\}_{s=1}^S$ 表示，由粒子滤波器推进。本文的核心计算想法是把**两层**期望（$S$ 个粒子与 $\Lambda$ 条 Monte Carlo 路径）压成每个梯度步**一次**随机抽样，也就是随机梯度下降：

$$
Y^{(\hat l,\hat s)}_i=Y^{(\hat l,\hat s)}_{i+1}+\Bigl[b_x^\top Y^{(\hat l,\hat s)}_{i+1}+\sigma_x^\top Z^{(\hat l,\hat s)}_{i+1}+f_x^\top\Bigr]\Delta t_i,
\qquad
Z^{(\hat l,\hat s)}_i=Y^{(\hat l,\hat s)}_{i+1}\,\omega^{(\hat l,\hat s)}_i\,(\Delta t_i)^{-1/2},
$$

起点 $X^{(\hat l,\hat s)}_{t_n}=x^{(\hat s)}_n$ 自粒子云中抽取，$\omega^{(\hat l,\hat s)}_i\sim N(0,1)$。

> [!warning] 论文自己承认这不解准 FBSDE
> 作者明确表示这个格式**并不精确求解 FBSDE**；FBSDE 在这里只是梯度过程的载体，他们诉诸的是随机梯度下降的常规理由（无偏的含噪梯度）。**这正是它彻底消掉 $(Y,Z)$ 关于 $X$ 的空间逼近的地方，也是它与编号 26 在目标上的分野——26 买的是精度与理论，50 买的是速度与维数。**

整体算法（PF-SGD）的结构是：初始化粒子云与迭代次数 $L$；对 $n=0,1,\dots,N_T$，初始化估计控制过程与步长，对 $l=0,\dots,L$ 模拟**一条**轨道、倒推算 $\{Y^{(l,s)}_i\}$、更新控制；取 $\hat u^\ast(t_n)=u^{L,M}_{t_n}|_{t_n}$；用它推进粒子滤波器。

### 定理

**预印本没有给出 PF-SGD 的收敛定理。** 它转而引用外部分析（对上述倒向格式及其推广的数值分析见其参考文献，即 Bao-Cao-Webster 型与赵卫东学派的分析），以及标准的随机梯度下降收敛文献。**因此这是一篇算法／计算的论文，不是误差分析的论文。** 唯一明确陈述的理论部件是附录中由随机最大值原理与 Gâteaux 导数推出的梯度表示。

### 数值实验

**例 1——带非线性观测的 LQ 基准。** $\mathrm dX_t=A(t)X_t\mathrm dt+BU_t\mathrm dt+C\,\mathrm dW_t$，代价 $J(U)=\mathbb E\bigl[\frac12\int_0^T(\langle QX_t,X_t\rangle+\langle RU_t,U_t\rangle)\mathrm dt+\frac12\langle FX_T,X_T\rangle\bigr]$，其**完全观测**下的最优控制是 Riccati 反馈 $\bar U_t=-R^{-1}B^\top P(t)X_t$，$\dot P=-PA-A^\top P+PBR^{-1}B^\top P-Q$、$P(T)=F$。观测是**非线性**的：$M_t=\sin(X_t)+\eta_t$，高斯噪声协方差 $\Gamma$。用途是检验数据驱动控制能否跟上解析控制。**没有报告任何收敛阶。**

**例 2——一维非线性非二次控制，效率研究。** $\mathrm dX_t=\arctan(X_t+u_t)\mathrm dt+\sigma X_t\mathrm dW_t$，$\sigma=0.05$；代价 $J^\ast(u^M)=\mathbb E\bigl[\frac12\int_0^T\sin^2(X_t+u^M_t)\mathrm dt\bigr]$；观测 $M_t=X_t+\eta_t$，噪声标准差 $0.05$。基准是「完整求解法」（Zakai 方程滤波器加上带网格插值的完整 FBSDE 求解），在三套网格上运行。$T=1$；PF-SGD 用 $\Delta t=0.02$（$N_T=50$）、500 个粒子、1000 次 SGD 迭代；完整求解法的三套网格为 $\Delta t=0.1/\Delta x=0.1$、$\Delta t=0.05/\Delta x=\frac{\sqrt2}{2}\cdot0.1$、$\Delta t=0.025/\Delta x=0.05$，定义域 $[3,6]$。

| 方法            | 粗网格 | 较细网格 | 最细网格 | PF-SGD      |
| --------------- | ------ | -------- | -------- | ----------- |
| 总代价 $J^\ast$ | 0.0481 | 0.0318   | 0.0076   | **0.00095** |
| CPU 时间（秒）  | 29.78  | 220.47   | 1560.15  | **0.93**    |

**这张表是全文最强的一处证据**：PF-SGD 的代价比最细网格的完整求解法还低一个数量级，而耗时只有它的约 $1/1700$。另有一组重复实验的平均代价 $\hat J^\ast_t(\hat u)=\frac{1}{M_{\rm rept}}\sum_m\frac12\int_0^t\sin^2(\hat X^{(m)}_s+\hat u^{(m)}_s)\mathrm ds$（$M_{\rm rept}=50$ 次）在时间上确认了同样的次序。

**例 3——只有方位观测的 Dubins 车。** $\mathrm dX_t=\sin\theta_t\mathrm dt+\sigma\mathrm dW_t$、$\mathrm dY_t=\cos\theta_t\mathrm dt+\sigma\mathrm dW_t$、$\mathrm d\theta_t=u_t\mathrm dt+\sigma_2\mathrm dW_t$，$\sigma=0.2$；代价 $J^\ast(u^M)=\mathbb E\bigl[\int_0^T\frac12(u^M_t)^2\mathrm dt+\delta((X_T-X_P)^2+(Y_T-Y_P)^2)\bigr]$，$\delta=10$，目标点 $(X_P,Y_P)=(5,3)$。两个方位探测器位于 $(6,1)$ 与 $(-1,4)$，给出 $M_t=[\arctan\frac{X_t-6}{Y_t-1},\ \arctan\frac{X_t+1}{Y_t-4}]^\top+\eta_t$，噪声标准差 $0.01I_2$。$T=1$、$\Delta t=0.02$（$N_T=50$）、1000 个粒子、1000 次 SGD 步。只给出定性的轨迹图。

**这些实验建立了什么，又差在哪里。** 建立的是：在部分观测且观测非线性的问题上，PF-SGD 能给出可用的反馈控制，而且在例 2 这个可以与完整求解法对照的场合，它在代价与时间两项上同时占优。差在三处。第一，**没有任何收敛阶被报告，也没有收敛定理**，因此「更低的代价」是这一组问题上的观察，不是一般保证。第二，例 3 只有轨迹图，没有量化指标，因此三维状态上的表现无法与例 2 的结论对接。第三，例 2 的对照是与一套特定实现的完整求解法比，而后者用的是网格插值——**这张表说明的是「在这个问题上，放弃空间表示比细化空间表示更划算」，不是「PF-SGD 在一般意义下更准」。**

### 与其他论文的关系

**编号 26 的部分观测对应物。** 两篇都从同一个伴随 FBSDE 梯度表示出发；26 在完全观测下做梯度投影，用 Gauss-Hermite 求积加插值，并给出被证明的收敛率；50 把全部条件期望换成单样本随机梯度，并加一个粒子滤波器处理观测滤链。在条件期望的处理上它与编号 8、19、25、41 形成鲜明对照：那几篇用求积加（稀疏）网格插值，50 有意避开对 $(Y,Z)$ 的任何空间表示。它的倒向一步是编号 8 多步族的 $k=1$ 成员，而论文**不追求**高阶时间精度——它的目标是实时反馈。Feng Bao 与 Richard Archibald 不属于山东／中科院的多步一支；**这一篇是周涛的 FBSDE 控制机器与非线性滤波／资料同化社群的交汇点。** 它也指向编号 86 与 96：同样放弃基于网格的空间表示，但把 Monte Carlo／粒子表示换成神经网络。

## 51：非光滑约束下的分块迭代

### 直觉

这一篇处理的不是 FBSDE，而是把常微分方程与互补条件耦合在一起的动态非线性互补问题。用后向 Euler 离散后，每个时间点上都要解一个耦合的非线性系统，而它的两块结构截然不同：一块是光滑的微分系统，一块是非光滑的互补系统。**已有做法要么把两块合成一块（直接消元、半光滑 Newton），要么就没有收敛保证；这一篇的想法是不要合并，交替解。**

Gauss-Seidel 型的交替是最朴素的：用旧的 $x$ 解互补系统得新的 $y$，再用新的 $y$ 解微分系统得新的 $x$。这样两个子系统各自都可以直接调用现成的求解器。论文真正的贡献是它一直缺失的那一半：**为什么会收敛，以及收敛多快。** 而它给出的答案有一个漂亮的结构——固定时间区间时收敛是**超线性**的，且速率**与步长 $h$ 无关**；固定时间点数时则是 $O(h)$ 线性的。

### 问题设定

动态非线性互补问题为

$$
\dot x(t)=F(t,x(t),y(t)),
\qquad
0\le y(t)\ \perp\ G(t,x(t),y(t))\ge0,
\qquad t\in(0,T),\ x(0)=x_0,
$$

$x(t)\in\mathbb R^m$、$y(t)\in\mathbb R^n_+$。两个子类分别是**微分半仿射系统** $0\le y\perp Nx+My+g(t)\ge0$ 与**动态线性互补问题** $\dot x=Ax+By+f(t)$、$0\le y\perp Nx+My+g\ge0$。后向 Euler 离散给出

$$
0\le y_j\ \perp\ G(t_j,x_j,y_j)\ge0,
\qquad
x_j=x_{j-1}+hF(t_j,x_j,y_j),
\qquad j=1,\dots,N_t,\ h=T/N_t .
$$

论文指出已有两条主流做法各有确定的缺陷。**直接消元法**只在线性情形有效：消去 $x_j$ 后得到的约化线性互补系统的矩阵 $M_h:=hN(I-hA)^{-1}B+M$ **即使 $M$ 是 P-矩阵也可能不是**，而且形成它需要 $n$ 次大规模线性求解。**半光滑 Newton 法**只有局部收敛性，且需要 Clarke 广义 Jacobi 矩阵，对大规模问题代价高且不便。

### 推导

**第一步：单点 Gauss-Seidel 迭代。**

$$
0\le y^{k+1}_j\ \perp\ G\bigl(t_j,\,x^{k}_j,\,y^{k+1}_j\bigr)\ge0,
\qquad
x^{k+1}_j=x_{j-1}+hF\bigl(t_j,\,x^{k+1}_j,\,y^{k+1}_j\bigr).
$$

互补系统用**旧的** $x^k_j$，微分系统用**新的** $y^{k+1}_j$——这就是「Gauss-Seidel 式」的含义。两个子系统解耦，因此任何现成的 NCP 求解器（线性化方法、PATH 求解器）与任何针对光滑 $F$ 系统的 Newton 求解器都可以原样插入。

**第二步：多点（窗口）版本，也就是可并行的那个。** 把 $\{t_1,\dots,t_{N_t}\}$ 分成 $P$ 组、每组 $J=N_t/P$ 个点，在一组之内

$$
0\le y^{k+1}_j\ \perp\ G\bigl(t_j,x^k_j,y^{k+1}_j\bigr)\ge0,
\qquad
x^{k+1}_j=x^{k+1}_{j-1}+hF\bigl(t_j,x^{k+1}_j,y^{k+1}_j\bigr),
\qquad j=1,\dots,J .
$$

**关键在于 $\{y^{k+1}_j\}_{j=1}^J$ 互相独立**，因此 $J$ 个互补系统可以**并行**求解；$x$ 的递推随后按 $j$ 顺序进行。论文还写明几个特例：半仿射版本；动态线性互补版本，其每次迭代恰好只需一次线性求解 $(I-hA)x^{k+1}_j=x_{j-1}+hBy^{k+1}_j+hf_j$；以及 $M$ 为 Z-矩阵时，$y^{k+1}_j$ 可由线性规划 $\min\|y\|_1$（约束 $y\ge0$、$My+Nx^k_j+g_j\ge0$）得到。

**第三步：假设。** $G$ 关于 $y$ 是**一致 P-函数**：

$$
\max_{1\le l\le n}(\bar y_l-\tilde y_l)\bigl(G_l(t,x,\bar y)-G_l(t,x,\tilde y)\bigr)\ \ge\ L_0\|\bar y-\tilde y\|_2^2,\qquad L_0>0,
$$

并满足 $\|G(t,x,\bar y)-G(t,x,\tilde y)\|_2\le L_G\|\bar y-\tilde y\|_2$；$F$ 关于 $x$ 满足**单边 Lipschitz 条件**、关于 $y$ 满足普通 Lipschitz 条件：

$$
\langle F(t,\bar x,y)-F(t,\tilde x,y),\bar x-\tilde x\rangle\le L_1\|\bar x-\tilde x\|_2^2,
\qquad
\|F(t,x,\bar y)-F(t,x,\tilde y)\|_2\le L_2\|\bar y-\tilde y\|_2,
$$

其中 $L_1\in(-\infty,\infty)$ **可以为负**，$L_2>0$。引理 2.1 给出唯一的 $\mathcal Y(x)$ 解静态 NCP，且关于 $x$ Lipschitz，常数为 $\eta_t$。

### 定理

- **引理 2.2（组合恒等式）。** 记 $\psi(r,J,k)=\sum_{j_1=1}^{J}\sum_{j_2=1}^{j_1}\cdots\sum_{j_k=1}^{j_{k-1}}r^{J-j_k}$，则
  $$
  \psi(r,J,k)=
  \begin{cases}
  \dfrac{1}{(1-r)^k}-\displaystyle\sum_{l=1}^{k}\binom{J+k-l-1}{k-l}\dfrac{r^{J}}{(1-r)^l}, & r\neq1,\\[10pt]
  \dbinom{J+k-1}{k}, & r=1 .
  \end{cases}
  $$
  **这条恒等式是全部收敛结果的技术枢纽**：多点迭代的误差递推是一个 $k$ 重嵌套求和，$\psi$ 正是把它求成闭式的东西。
- **定理 2.3（基本误差界）。** 在上述假设与 $hL_1<1$ 下，多点方法的误差 $e^k_j=x_j-x^k_j$ 满足
  $$
  \max_{0\le j\le J}\|e^k_j\|_2\le
  \begin{cases}
  (h\tilde\eta)^k\,\psi(1,J,k)\,\max_{0\le j\le J}\|e^0_j\|_2, & L_1=0,\\[6pt]
  \psi\bigl((1-hL_1)^{-1},J,k\bigr)\Bigl(\dfrac{h\tilde\eta}{1-hL_1}\Bigr)^{k}\max_{1\le j\le J}\|e^0_j\|_2, & L_1\neq0,
  \end{cases}
  $$
  其中 $\tilde\eta=L_2\eta$，$x_j$ 是收敛后的（精确后向 Euler）解。
- **定理 2.4（固定时间点数 $\Rightarrow$ $O(h)$ 线性）。** 在定理 2.3 的假设下，当 $k\gg1$、$J$ 固定且 $h(L_1+\tilde\eta)<1$ 时，
  $$
  \max_{0\le j\le J}\|e^k_j\|_2\le\rho^k\max_{0\le j\le J}\|e^0_j\|_2,
  \qquad
  \rho=\frac{h\tilde\eta}{1-hL_1}=O(h),\quad\tilde\eta=\eta L_2 .
  $$
- **定理 2.5（固定时间区间 $\Rightarrow$ 超线性、与 $h$ 无关）。**
  $$
  \|x^k(t)-x(t)\|_2\le\max\{1,e^{L_1t}\}\,\frac{(t\tilde\eta)^k}{k!}\,\sup_{t\in[0,T]}\|x^0(t)-x(t)\|_2,
  \qquad t\in(0,T).
  $$
  **分母上的 $k!$ 使收敛超线性，且速率与 $h$ 无关**（论文的 Remark 2.2）。
- **Remark 2.3。** $L_1$ 为**负**会加速收敛；对动态线性互补问题，这发生在 $A$ 正交相似于其 Jordan 形且特征值实部为负时——**这正是常微分方程来自抛物型 PDE 半离散的典型情形**（例如抛物 Signorini 问题）。
- **$P_0$ 情形。** 若 $G$ 只是 $P_0$-函数（$L_0=0$），上述结果在 Tikhonov 正则化 $\widehat G(t,x,y)=G(t,x,y)+\varepsilon y$ 之后仍适用，$\varepsilon$ 取得与时间误差 $O(h)$ 同量级。
- **线性互补情形（第 3 节，$G(t,x,y)=My(t)+\widetilde G(t,x(t))$）。** 定理 3.2 给出**最小范数解的唯一存在性**：对半正定的 $M$，在 $(0,T^\ast)$ 上存在唯一的 $(x,y)\in C^1(0,t)\times C(0,t)$，且时间视界为
  $$
  T^\ast=
  \begin{cases}
  T, & L_1\le-\dfrac{\tilde\eta\beta+C_0}{\beta},\\[8pt]
  \min\Bigl\{T,\ \dfrac1{L_1}\log\Bigl(1+\dfrac{L_1\beta}{\tilde\eta\beta+C_0}\Bigr)\Bigr\}, & L_1>-\dfrac{\tilde\eta\beta+C_0}{\beta},
  \end{cases}
  $$
  其中 $\tilde\eta=L_2\eta_0\eta_1$、$C_0=\max_{t\in[0,T]}\|F(x_0,\mathcal Y(\widetilde G(t,x_0)))\|_2$、$\mathcal B(x_0,\beta)=\{v:\|v-x_0\|_2\le\beta\}$，并假设在该球上 $\mathrm{FEA}(M,\widetilde G(t,v))\neq\emptyset$。定理 3.3 表明最小范数迭代法是**良定的**（可行性沿迭代保持），条件是 $hL_1<1$ 与对可用时间点数 $J^\ast$ 的一个明确限制。$M$ 为 Z-矩阵的情形另在 Remark 3.1 中处理。

### 数值实验

停止准则为 $\max_{0\le j\le J}\|x^k_j-x_j\|_2\le10^{-8}$，比较对象是收敛后的解 $\{x_j\}$。两个例子分别是：**四二极管桥式整流器**（含非线性电阻与一个**取随机值**的电容的非光滑电路——周涛的不确定性量化兴趣在这里进来），以及**由空间价格均衡问题导出的投影动力系统**。

报告的结论是该方法在稳健性、复杂度与计算时间三方面都优于已有方法，且两个理论区间——固定 $T$ 时对 $h$ 稳健的超线性收敛，与固定 $J$ 时的 $O(h)$ 线性收敛——都被观察到。**本站未转录各实验的具体迭代次数与耗时表，因此本页不给出量化对比。**

### 与其他论文的关系

**它与 FBSDE 一线基本不相交**：没有 BSDE，没有条件期望，没有概率型格式。唯一的结构性亲缘是**分裂／解耦**的哲学——交替求解两个耦合子系统，正如编号 8、23、26 的耦合 FBSDE 求解器在正向 SDE 与倒向方程之间迭代。编号 26 的梯度投影外层同样是对耦合系统作不动点迭代并投影到凸容许集上，**在这份列表里它是数学上最近的亲戚**。吴树林是并行时间（波形松弛／parareal）方向的专家，多点窗口版本正是波形松弛型的想法，这也是「$J$ 个时间点上的并行」为什么重要的原因。第 4.1 节的随机电路参数连的是周涛的不确定性量化工作，不是 FBSDE 工作。

## 七篇的定位

| 编号 | 处理的对象               | 相对其他篇的位置             | 本站核实程度   |
| ---- | ------------------------ | ---------------------------- | -------------- |
| 16   | 全非线性抛物方程         | 二阶设定的起点，PDE 优先框架 | 摘要与参考文献 |
| 19   | 二阶 FBSDE 与随机控制    | 把二阶设定接到控制问题       | 全文逐式       |
| 25   | 多维条件期望的求值代价   | 空间表示方向的改进           | 全文逐式       |
| 26   | 带约束的随机最优控制     | 外层迭代（梯度投影），一阶   | 全文逐式       |
| 41   | 随机最优控制的高精度格式 | 内层求解（高阶 FBSDE），二阶 | 摘要与参考文献 |
| 50   | 数据驱动反馈控制         | 部分观测，滤波与控制耦合     | 全文逐式       |
| 51   | 动态非线性互补问题       | 非光滑约束下的分块迭代       | 全文逐式       |

这一页的形状有两条。第一条是**同一个控制问题的两条攻法**：编号 19 走动态规划（HJB → 2FBSDE → 解析求 $\inf_\alpha$），编号 26 与 41 走 Pontryagin（伴随 BSDE → 梯度 → 投影或拟 Newton）。前者需要 $\inf_\alpha$ 可解析求出，后者需要梯度可算；两者的困难位置完全不同。第二条是**代价从哪里省**：编号 25 在网格内部省（稀疏网格把指数压成多项式），编号 50 则整个放弃网格（单样本随机梯度），并因此放弃了精度保证。

## 覆盖核对

| 内容                                          | 论文 | 覆盖状态                                 |
| --------------------------------------------- | ---- | ---------------------------------------- |
| CSTV 二阶设定、$\Gamma$ 的作用、两处记号差别  | 背景 | 完整                                     |
| 正向 SDE 可自由选择这条设计自由               | 16   | 仅摘要可支撑的内容；方程按编号 19 反推   |
| 表示定理、四条参考 ODE、五个格式、Newton 迭代 | 19   | 完整推导                                 |
| 截断误差、12 项分解、$k\le6$ 稳定范围         | 19   | 完整；并注明缺收敛定理                   |
| 四分量收敛率表与效率对比                      | 19   | 表 5.4 与例 5.1 完整；控制例子的表未转录 |
| 稀疏网格、分层基、快速变换、稀疏求积          | 25   | 完整推导                                 |
| 二维与 $q$ 维实验、多项式增长                 | 25   | 例 1 完整；例 2、3 仅设定与结论          |
| 不动点刻画、伴随 BSDE、左端点矩形格式         | 26   | 完整推导，含 Remark 1 的批评             |
| 定理 1、推论 1、定理 2 与误差平衡             | 26   | 完整，含假设                             |
| 反馈控制对比表与一阶收敛                      | 26   | 例 3 完整；例 1、2、4 部分转录           |
| 二阶求解器与拟 Newton 优化器                  | 41   | 仅摘要与参考文献可支撑的内容             |
| 部分观测设定、伴随系统、PF-SGD 单实现技巧     | 50   | 完整推导，含论文自陈的局限               |
| 效率对照表与三个例子                          | 50   | 例 2 完整；例 1、3 仅设定与定性结论      |
| Gauss-Seidel 迭代、多点并行版本、四条定理     | 51   | 完整                                     |
| 两个数值例子                                  | 51   | 仅设定与结论；具体数据表未转录           |

## 本页原文

- T. Kong, W. Zhao, and T. Zhou, [_Probabilistic high order numerical schemes for fully nonlinear parabolic PDEs_](https://doi.org/10.4208/cicp.240515.280815a), Commun. Comput. Phys. 18(5) (2015), pp. 1482-1503。
- T. Kong, W. Zhao, and T. Zhou, [_High order numerical schemes for second-order FBSDEs with applications to stochastic optimal control_](https://doi.org/10.4208/cicp.OA-2016-0056), Commun. Comput. Phys. 21(3) (2017), pp. 808-834（预印本 [arXiv:1502.03206](https://arxiv.org/abs/1502.03206)）。
- Y. Fu, W. Zhao, and T. Zhou, [_Efficient spectral sparse grid approximations for solving multi-dimensional forward backward SDEs_](https://doi.org/10.3934/dcdsb.2017174), Discrete Contin. Dyn. Syst. Ser. B 22(9) (2017), pp. 3439-3458（预印本 [arXiv:1607.06897](https://arxiv.org/abs/1607.06897)）。
- B. Gong, W. Liu, T. Tang, W. Zhao, and T. Zhou, [_An efficient gradient projection method for stochastic optimal control problems_](https://doi.org/10.1137/17M1123559), SIAM J. Numer. Anal. 55(6) (2017), pp. 2982-3005。
- Y. Fu, W. Zhao, and T. Zhou, [_Highly accurate numerical schemes for stochastic optimal control via FBSDEs_](https://doi.org/10.4208/nmtma.OA-2019-0137), Numer. Math. Theor. Meth. Appl. 13(2) (2020), pp. 296-319。
- R. Archibald, F. Bao, J. Yong, and T. Zhou, [_An efficient numerical algorithm for solving data driven feedback control problems_](https://doi.org/10.1007/s10915-020-01358-y), J. Sci. Comput. 85(2) (2020), 58（预印本 [arXiv:2006.03047](https://arxiv.org/abs/2006.03047)）。
- S. Wu, T. Zhou, and X. Chen, [_A Gauss-Seidel type method for dynamic nonlinear complementarity problems_](https://doi.org/10.1137/19M1268884), SIAM J. Control Optim. 58(6) (2020), pp. 3389-3412。
- 用于交叉印证的外部来源：P. Cheridito, H. M. Soner, N. Touzi, and N. Victoir, [_Second-order backward stochastic differential equations and fully nonlinear parabolic PDEs_](https://doi.org/10.1002/cpa.20168), Comm. Pure Appl. Math. 60(7) (2007), pp. 1081-1110（本站据此转录二阶设定的定义、PDE 对应与「正向扩散可自由选择」的原始表述）。
