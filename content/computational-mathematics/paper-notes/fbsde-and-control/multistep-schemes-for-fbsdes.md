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
> 编号 **8**（_SIAM J. Sci. Comput._ 36(4), 2014）、**18**（_J. Sci. Comput._ 69(2), 2016）、**23**（_Numer. Math. Theor. Meth. Appl._ 10(2), 2017）、**33**（_SIAM J. Numer. Anal._ 56(4), 2018）、**35**（_J. Sci. Comput._ 79(3), 2019）、**61**（_J. Comput. Math._ 40(4), 2022）、**68**（_J. Sci. Comput._ 94:53, 2023）。
>
> 其中编号 **8**（据 arXiv 预印本全文）与 **23**（据已出版正文）已逐式核对，本页给出它们的完整推导、定理与数值实验表。编号 **18、33、35、61、68** 的正文未能获取——五篇均无预印本，Springer 与 SIAM 需订阅，Global Science Press 的 PDF 端点返回 403 或不可达——因此这五篇只写出摘要与参考文献可支撑的内容，其格式、定理、收敛阶与实验数据本站不报告。

![多步格式如何提高倒向方程的时间精度](assets/diagrams/tao-zhou-papers/zh/fbsde-multistep.svg)

## 8：把高阶要求全部转移到倒向方向

### 直觉

把一个正倒向系统数值求解，直觉上应当两个方向都精确：正向的 $X$ 算不准，倒向的 $Y,Z$ 就无从谈起。但倒向方程只通过**条件期望**看见正向过程——格式里出现的从来不是某一条 $X$ 的轨道，而是 $\mathbb E^x_{t_n}[\,\cdot\,]$ 这样的量。条件期望只需要**转移律**，不需要轨道。而在一个时间步的尺度上，从 $(t_n,x)$ 出发的扩散的转移律，其一阶变化只由**该点处的系数**决定：那正是生成元的含义。

这条观察是编号 8 的全部机制。它意味着可以把真实的正向扩散换成一条系数被冻结在 $(t_n,x)$ 的扩散——也就是 Euler 步——而不损失倒向方向上的任何阶。少了这条观察会发生什么，论文说得很清楚：为了让 $Y,Z$ 达到 $k$ 阶，就必须为正向 SDE 也构造 $k$ 阶强收敛格式，而这类格式（Itô-Taylor 型）需要重复的多重随机积分，计算量大且在耦合情形下几乎无法实施——耦合时正向系数依赖尚未算出的 $(Y,Z)$。

剩下的问题就纯粹是倒向方向的了：怎样用多个未来时间层的信息把时间阶提上去。论文的答案不是对被积函数插值，而是**先把积分恒等式对 $t$ 求导变成两条常微分方程，再用一个前向的 BDF 型模板逼近导数**。这个选择带来一个意外的红利：在 $t=t_n$ 处，所有形如 $f\cdot\Delta W_{t_n,t}$ 的项因 $\Delta W_{t_n,t_n}=0$ 而消失，于是 $Z^n$ 是**显式**的，只有 $Y^n$ 需要迭代。

### 问题设定

耦合的 Markov 型正倒向随机微分方程为

$$
X_t=X_0+\int_0^t b(s,X_s,Y_s,Z_s)\,\mathrm ds+\int_0^t\sigma(s,X_s,Y_s,Z_s)\,\mathrm dW_s,
$$

$$
Y_t=\xi+\int_t^T f(s,X_s,Y_s,Z_s)\,\mathrm ds-\int_t^T Z_s\,\mathrm dW_s,
\qquad \xi=\varphi(X_T),
$$

其中 $b:\Omega\times[0,T]\times\mathbb R^q\times\mathbb R^p\times\mathbb R^{p\times d}\to\mathbb R^q$、$\sigma:\cdots\to\mathbb R^{q\times d}$、生成元 $f:\cdots\to\mathbb R^p$，未知三元组 $(X_t,Y_t,Z_t)$ 取值于 $\mathbb R^q\times\mathbb R^p\times\mathbb R^{p\times d}$。当 $b,\sigma$ 不依赖 $(Y_t,Z_t)$ 时系统称为解耦。

论文对当时状况的描述是：多数已有格式是 Euler 型、收敛率 $1/2$；对解耦系统的高阶方法「依赖同时对正向随机微分方程与倒向方程使用高阶方案」，而正向的高阶方案「计算量大且往往难以实施」。耦合情形下正向系数依赖倒向未知量，因此「似乎不容易设计高阶且高效的数值格式」。论文明确提出并肯定回答的问题是：**如果正向随机微分方程只用 Euler 方法求解，倒向方程还能不能达到高阶精度？**

### 推导

**第一步：生成元只看左端点。** 设 $X_s$ 解 $\mathrm dX_s=b\,\mathrm ds+\sigma\,\mathrm dW_s$，其生成元为

$$
\mathcal A_t^x g(t,x)=\lim_{s\downarrow t}\frac{\mathbb E_t^x[g(s,X_s)]-g(t,x)}{s-t},
\qquad
\mathcal L^0_{t,x}=\frac{\partial}{\partial t}+\sum_i b_i\frac{\partial}{\partial x_i}
+\frac12\sum_{i,j}(\sigma\sigma^{\top})_{i,j}\frac{\partial^2}{\partial x_i\partial x_j},
$$

且在 $C^{1,2}$ 函数上 $\mathcal A_t^x f=\mathcal L^0_{t,x}f$。论文的定理 2 是：若 $f\in C^{1,2}$ 且 $\mathbb E^{x_0}_{t_0}[\mathcal L^0_{t,X_t}f(t,X_t)]<\infty$，则对 $t\ge t_0$

$$
\frac{\mathrm d\,\mathbb E^{x_0}_{t_0}[f(t,X_t)]}{\mathrm dt}
=\mathbb E^{x_0}_{t_0}\bigl[\mathcal A_t^{X_t}f(t,X_t)\bigr],
$$

并且**在 $t=t_0$ 处**该导数与把 $X_t$ 换成**任何**系数仅在左端点匹配的扩散 $\bar X_t=x+\int_{t_0}^t\bar b_s\,\mathrm ds+\int_{t_0}^t\bar\sigma_s\,\mathrm dW_s$ 后的导数相同：

$$
\bar b(t_0,\bar X_{t_0};t_0,x_0)=b(t_0,x_0),
\qquad
\bar\sigma(t_0,\bar X_{t_0};t_0,x_0)=\sigma(t_0,x_0)
\quad\Longrightarrow\quad
\left.\frac{\mathrm d\,\mathbb E[f(t,X_t)]}{\mathrm dt}\right|_{t_0}
=\left.\frac{\mathrm d\,\mathbb E[f(t,\bar X_t)]}{\mathrm dt}\right|_{t_0}.
$$

论文明确指出可以简单取 $\bar b(s,\bar X_s;t_0,x_0)=b(t_0,x_0)$、$\bar\sigma(s,\bar X_s;t_0,x_0)=\sigma(t_0,x_0)$ 对所有 $s\in[t_0,t]$ 成立——这就是 Euler 步。**这条定理的一般意义值得单独指出：** 当一个耦合系统的两个方向精度要求不同时，先确认低精度方向的误差是否真的进入目标量。这里的答案是不会，因为进入目标量的只是一步导数逼近，而它只看左端点。

**第二步：导数逼近的权。** 对 $u\in C_b^{k+1}$ 与节点 $t_0<t_1<\cdots<t_k$、$\Delta t_i=t_i-t_0$，权由矩条件

$$
\sum_{i=0}^{k}\alpha_{k,i}\frac{(\Delta t_i)^j}{j!}=\delta_{j1}
=\begin{cases}1,&j=1\\0,&j\ne1\end{cases},
\qquad j=0,1,\dots,k
$$

确定，给出 $\frac{\mathrm du}{\mathrm dt}(t_0)=\sum_{i=0}^k\alpha_{k,i}u(t_i)+\mathrm{Err}$，$\mathrm{Err}=O\bigl(\sum_i\alpha_{k,i}(\Delta t_i)^{k+1}\bigr)$。在均匀网格 $\Delta t_i=i\Delta t$ 上这是一个 Vandermonde 系统

$$
\begin{pmatrix}
1&1&\cdots&1\\ 0&1&\cdots&k\\ 0&1^2&\cdots&k^2\\ \vdots&&&\vdots\\ 0&1^k&\cdots&k^k
\end{pmatrix}
\begin{pmatrix}\alpha_{k,0}\Delta t\\ \alpha_{k,1}\Delta t\\ \vdots\\ \alpha_{k,k}\Delta t\end{pmatrix}
=\begin{pmatrix}0\\1\\0\\\vdots\\0\end{pmatrix}.
$$

解出的 $\alpha_{k,i}\Delta t$ 恰是**反向读的 BDF 系数**（即时间朝前的向后差分公式）：

| $k$ | $i=0$     | $1$ | $2$     | $3$    | $4$     | $5$   | $6$    |
| --- | --------- | --- | ------- | ------ | ------- | ----- | ------ |
| 1   | $-1$      | $1$ |         |        |         |       |        |
| 2   | $-3/2$    | $2$ | $-1/2$  |        |         |       |        |
| 3   | $-11/6$   | $3$ | $-3/2$  | $1/3$  |         |       |        |
| 4   | $-25/12$  | $4$ | $-3$    | $4/3$  | $-1/4$  |       |        |
| 5   | $-137/60$ | $5$ | $-5$    | $10/3$ | $-5/4$  | $1/5$ |        |
| 6   | $-49/20$  | $6$ | $-15/2$ | $20/3$ | $-15/4$ | $6/5$ | $-1/6$ |

**第三步：两条参考常微分方程。** 这是论文的核心结构。对倒向方程取 $\mathbb E^x_{t_n}[\cdot]$ 得积分恒等式 $\mathbb E^x_{t_n}[Y_t]=\mathbb E^x_{t_n}[\xi]+\int_t^T\mathbb E^x_{t_n}[f(s,X_s,Y_s,Z_s)]\,\mathrm ds$，对 $t$ 求导：

$$
\frac{\mathrm d\,\mathbb E^x_{t_n}[Y_t]}{\mathrm dt}
=-\,\mathbb E^x_{t_n}\bigl[f(t,X_t,Y_t,Z_t)\bigr],\qquad t\in[t_n,T].
$$

把 $Y_{t_n}=Y_t+\int_{t_n}^tf\,\mathrm ds-\int_{t_n}^tZ_s\,\mathrm dW_s$ 乘以 $(\Delta W_{t_n,t})^{\top}$ 再取 $\mathbb E^x_{t_n}[\cdot]$（Itô 等距使 $\mathbb E^x_{t_n}[Y_{t_n}(\Delta W_{t_n,t})^\top]=0$ 与随机积分项化为 $\int_{t_n}^t\mathbb E^x_{t_n}[Z_s]\mathrm ds$），得

$$
0=\mathbb E^x_{t_n}\bigl[Y_t(\Delta W_{t_n,t})^{\top}\bigr]
+\int_{t_n}^t\mathbb E^x_{t_n}\bigl[f(s,X_s,Y_s,Z_s)(\Delta W_{t_n,s})^{\top}\bigr]\mathrm ds
-\int_{t_n}^t\mathbb E^x_{t_n}[Z_s]\,\mathrm ds,
$$

再对 $t$ 求导给出第二条参考方程：

$$
\frac{\mathrm d\,\mathbb E^x_{t_n}\bigl[Y_t(\Delta W_{t_n,t})^{\top}\bigr]}{\mathrm dt}
=-\,\mathbb E^x_{t_n}\bigl[f(t,X_t,Y_t,Z_t)(\Delta W_{t_n,t})^{\top}\bigr]+\mathbb E^x_{t_n}[Z_t].
$$

**在 $t=t_n$ 处取值是关键一步**：$\Delta W_{t_n,t_n}=0$ 使含 $f\cdot\Delta W$ 的项整体消失，两条方程的右端只剩 $-f(t_n,x,Y_{t_n},Z_{t_n})$ 与 $Z_{t_n}$。这就是下面的格式里未来层上**没有 $f$ 项**的原因，也是 $Z$ 方程完全显式的原因。

**第四步：半离散格式。** 把导数逼近权用到两条参考方程的左端并借定理 2 换成冻结扩散，得

$$
\sum_{i=0}^k\alpha_{k,i}\,\mathbb E^x_{t_n}\bigl[\bar Y_{t_{n+i}}\bigr]=-f(t_n,x,Y_{t_n},Z_{t_n})+R^k_{y,n},
\qquad
\sum_{i=1}^k\alpha_{k,i}\,\mathbb E^x_{t_n}\bigl[\bar Y_{t_{n+i}}(\Delta W_{n,i})^{\top}\bigr]=Z_{t_n}+R^k_{z,n}.
$$

丢掉截断项并取 Euler 选择 $\bar b(s,\cdot)=b(t_n,x)$、$\bar\sigma(s,\cdot)=\sigma(t_n,x)$，得到论文的格式 2：

$$
X^{n+j}=X^n+b(t_n,X^n)\Delta t_{n,j}+\sigma(t_n,X^n)\Delta W_{n,j},\qquad j=1,\dots,k,
$$

$$
Z^n=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}(\Delta W_{n,j})^{\top}\bigr],
\qquad
\alpha_{k,0}Y^n=-\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{X^n}_{t_n}\bigl[\bar Y^{n+j}\bigr]-f(t_n,X^n,Y^n,Z^n).
$$

一个容易读漏的细节：$X^{n+j}$ 对**所有** $j=1,\dots,k$ 都由从 $t_n$ 出发的**单个** Euler 步给出，步长 $\Delta t_{n,j}=t_{n+j}-t_n$、增量 $\Delta W_{n,j}=W_{t_{n+j}}-W_{t_n}$，而不是 $j$ 次连续 Euler 步。这正是定理 2 允许的事情，也是格式便宜的原因。

**第五步：全离散格式。** 引入空间网格 $\mathcal D^n_h\subset\mathbb R^q$，密度 $h_n=\max_x\mathrm{dist}(x,\mathcal D^n_h)$，局部邻点集 $\mathcal D^n_{h,x}$（$\#\mathcal D^n_{h,x}\le N_e$），局部插值算子 $\mathcal I^n_{h,x}$，以及条件期望的求积逼近 $\mathbb E^{x,h}_{t_n}[\cdot]$。给定 $i=0,\dots,k-1$ 时的 $Y^{N-i},Z^{N-i}$，对 $n=N-k,\dots,0$ 与每个 $x\in\mathcal D^n_h$：

$$
Z^n=\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{x,h}_{t_n}\bigl[\mathcal I^{n+j}_{h,X^{n+j}}Y^{n+j}(\Delta W_{n,j})^{\top}\bigr],
$$

$$
\alpha_{k,0}Y^{n}=-\sum_{j=1}^{k}\alpha_{k,j}\,\mathbb E^{x,h}_{t_n}\bigl[\mathcal I^{n+j}_{h,X^{n+j}}Y^{n+j}\bigr]-f(t_n,x,Y^n,Z^n).
$$

每个时间层的操作顺序是：先 Euler 步得 $X^{n+j}$，再显式算 $Z^n$，最后隐式解 $Y^n$，用 Picard 迭代 $\alpha_{k,0}Y^{n,l+1}=-\sum_j\alpha_{k,j}\mathbb E^{x,h}_{t_n}[\mathcal I Y^{n+j}]-f(t_n,x,Y^{n,l},Z^n)$ 直到 $|Y^{n,l+1}-Y^{n,l}|\le\epsilon_0$。**代价在这里也要记清：$k$ 步格式需要 $k$ 个起始层，且每层要访问 $k$ 个未来层。** 这是它与延迟校正（编号 23）在结构上的分野。

**第六步：耦合情形。** 格式 4 把 Euler 步换成 $X^{n+j}=X^n+b(t_n,X^n,Y^n,Z^n)\Delta t_{n,j}+\sigma(t_n,X^n,Y^n,Z^n)\Delta W_{n,j}$，此时 $X^{n+j}$ 与 $(Y^n,Z^n)$ 互相依赖。格式 5 用一层外 Picard 迭代解开：取 $Y^{n,0}=Y^{n+1}$、$Z^{n,0}=Z^{n+1}$，对 $l=0,1,\dots$ 先用 $(Y^{n,l},Z^{n,l})$ 走 Euler 步，再显式算 $Z^{n,l+1}$、隐式算 $Y^{n,l+1}$，直到 $\max(|Y^{n,l+1}-Y^{n,l}|,|Z^{n,l+1}-Z^{n,l}|)<\epsilon_0$。若 $b,\sigma$ 与 $(Y,Z)$ 无关，格式 5 退化为格式 3。

**第七步：条件期望怎么算。** Euler 步使 $\Delta W_{n,j}\sim\sqrt{\Delta t_{n,j}}\,N(0,I_d)$，于是所有条件期望都是对已知高斯密度的积分。论文用 Gauss-Hermite 求积

$$
\int_{-\infty}^{+\infty}e^{-x^2}g(x)\,\mathrm dx\approx\sum_{j=1}^{L}\omega_jg(a_j),
\qquad
\omega_j=\frac{2^{L+1}L!\sqrt\pi}{\bigl(H_L'(a_j)\bigr)^2},
$$

$\{a_j\}$ 为 $L$ 次 Hermite 多项式 $H_L$ 的根，截断误差 $R(g,L)=\frac{L!\sqrt\pi}{2^L(2L)!}g^{(2L)}(\eta)$，因此对次数 $\le2L-1$ 的多项式精确。$d$ 维用张量规则 $\omega_{\mathbf j}=\prod_i\omega_{j_i}$，并借 $\mathbb E[g(N)]=\pi^{-d/2}\int g(\sqrt2x)e^{-x^\top x}\mathrm dx$ 转换。空间上 $\mathcal I^n_{h,x}$ 取局部 Lagrange 插值。**求积节点是 Hermite 多项式的根，一般不落在空间网格上，所以插值这一步不可省**——这正是编号 63 后来用 Sinc 求积去掉的东西。

### 定理

- **生成元定理（定理 2）。** 假设 $f\in C^{1,2}$、$\mathbb E^{x_0}_{t_0}[\mathcal L^0_{t,X_t}f(t,X_t)]<\infty$。结论如上：条件期望的时间导数等于生成元的条件期望，且在左端点处只依赖左端点系数。这是全部构造的许可证。
- **局部截断误差。** 半离散层面 $\bar R^k_{y,n}=\bar R^k_{z,n}=O((\Delta t)^k)$。全离散层面截断误差分成六项：$R^k_{y,n},R^k_{z,n}$（导数逼近）、$R^{k,I_h}_{y,n},R^{k,I_h}_{z,n}$（插值）、$R^{k,E}_{y,n},R^{k,E}_{z,n}$（条件期望求积）。在数据光滑且插值多项式次数为 $r$ 的条件下，
  $$
  R^k_{y,n}=O\bigl((\Delta t_n)^k\bigr),\quad R^k_{z,n}=O\bigl((\Delta t_n)^k\bigr),\quad
  R^{k,I_h}_{y,n}=O\bigl(h^{r+1}\bigr),\quad R^{k,I_h}_{z,n}=O\bigl(h^{r+1}\bigr).
  $$
- **误差平衡。** 为使时间与空间两项同阶，论文取 $h=(\Delta t)^{(k+1)/(r+1)}$。
- **零稳定性给出阶数上限。** 窗口 $1\le k\le6$ 不是含糊的经验观察，而是一次可以复算的检验。把同一组权用到确定性常微分方程 $\mathrm dY/\mathrm dt=f(t,Y)$ 上得 $\alpha_{k,0}Y^n+\sum_{j=1}^k\alpha_{k,j}Y^{n+j}=f(t_n,Y^n)$，其特征多项式为

  $$
  P(\lambda)=\alpha_{k,0}\lambda^{k}+\sum_{j=1}^{k}\lambda^{k-j}=0,
  $$

  要求根条件 $|\lambda_{k,j}|\le1$（等号时须为单根）。论文报出的最大根模（除去公共根 $1.0$）是

  | $k$      | 2      | 3      | 4      | 5      | 6      | 7          | 8          |
  | -------- | ------ | ------ | ------ | ------ | ------ | ---------- | ---------- |
  | 最大根模 | 0.3333 | 0.4264 | 0.5608 | 0.7087 | 0.8633 | **1.0222** | **1.1839** |

  **根模在 $k=6$ 与 $k=7$ 之间越过 $1$，因此格式对 $k\ge7$ 不稳定**，这就是只列到 $k=6$ 的原因。

> [!warning] 这篇没有收敛定理
> arXiv 版本**没有**给出带显式常数的严格收敛定理。高阶性由上面的截断误差估计加数值证据支撑；论文谈到格式 4 时的原话只是「我们可以期望格式 4 是一个高阶数值格式……我们将用数值验证这一点」。已出版的 SIAM 版本是否补上定理，本站未核实。这一族的严格稳定性与收敛性分析要到编号 47 才出现，见[[computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes|稳定性理论一页]]。

### 数值实验

统一设定：$T=1.0$，均匀网格，8 个 Gauss-Hermite 节点（使求积误差可忽略），FORTRAN 95 实现。下表报的是观察到的收敛率，不是误差绝对值。

**例 1（解耦，有解析解）。** 精确解 $Y_t=\dfrac{e^{t+X_t}}{1+e^{t+X_t}}$、$Z_t=\dfrac{(e^{t+X_t})^2}{(1+e^{t+X_t})^3}$，$x=1.0$，$N=16,\dots,256$：

| $k$        | 1     | 2     | 3     | 4     | 5     | 6              | 7     | 8        |
| ---------- | ----- | ----- | ----- | ----- | ----- | -------------- | ----- | -------- |
| $Y$ 收敛率 | 1.000 | 1.973 | 3.002 | 3.922 | 5.196 | 5.116（6.273） | 4.382 | $-5.487$ |
| $Z$ 收敛率 | 1.000 | 2.021 | 2.893 | 3.919 | 5.017 | 5.687（6.256） | 4.759 | $-7.170$ |

括号内是 $N=16,\dots,128$ 上的收敛率，即双精度舍入误差主导之前的值——**$k=6$ 的名义阶只有在更粗的网格上才看得出来，细网格上舍入已经吃掉了它**。$k=7$ 的率退化到 $4$ 到 $5$ 之间而非接近 $7$，$k=8$ 出现负率（发散）。这与上面的根条件计算一致。

**例 2（Black-Scholes 下的欧式看涨期权）。** 参数 $b=0.05$、$\sigma=0.2$、$r=0.03$、$d=0.04$、$T=1.0$、$K=S_0=100$：

| $k$        | 1     | 2     | 3     | 4     |
| ---------- | ----- | ----- | ----- | ----- |
| $Y$ 收敛率 | 1.002 | 1.964 | 2.935 | 3.957 |
| $Z$ 收敛率 | 1.002 | 1.998 | 2.819 | 3.991 |

**例 3（耦合，两个系统）。** 系统 (5.10) 的精确解为 $Y_t=\sin(t+X_t)$、$Z_t=\sqrt2\cos(t+X_t)\sin^2(t+X_t)$，其扩散系数 $\sigma=\sqrt2\,Y_s\sin(s+X_s)$ **有意违反**论文的一致椭圆性假设 (2.11)；系统 (5.11) 的 $Z_t=\sqrt2\cos(t+X_t)(\sin^2(t+X_t)+1)$，满足该假设。

| $k$               | 1     | 2     | 3     | 4     |
| ----------------- | ----- | ----- | ----- | ----- |
| (5.10) $Y$ 收敛率 | 0.984 | 1.984 | 2.935 | 4.006 |
| (5.10) $Z$ 收敛率 | 0.977 | 2.064 | 2.976 | 4.055 |
| (5.11) $Y$ 收敛率 | 1.081 | 2.273 | 2.916 | 3.822 |
| (5.11) $Z$ 收敛率 | 1.023 | 1.948 | 2.970 | 4.213 |

例 4 把 $\sigma$ 推广到依赖 $(X_s,Y_s,Z_s)$；本站确认该例存在，但未转录其数据表。

**这些实验建立了什么，又差在哪里。** 建立的是：在光滑的一维问题上、包括一个耦合且**违反椭圆性假设**的问题上，$k$ 步格式确实呈现 $k$ 阶，且不稳定的起点与确定性根条件预言的位置吻合。差在三处。第一，全部例子的空间维数都很低，条件期望用张量 Gauss-Hermite 求积算，因此这些表说明不了维数增长时会发生什么——那要等编号 25 的稀疏网格。第二，$k=6$ 的率被舍入污染，说明**双精度下这一族的实际可用上限低于理论窗口的上限**。第三，例 3 中椭圆性被违反仍观察到 $k$ 阶，这是一条数值观察而非定理；论文没有、也不能据此放宽假设。

### 与其他论文的关系

这是整条线索的**根论文**：编号 16、18、19、23、25、41、61、68 都把自己描述为它的推广。编号 16 与 19 推到二阶（全非线性）情形，18 加跳，25 换空间离散为稀疏网格，23 与 35 换成延迟校正，47 补上严格的稳定性与收敛性框架，63 换成 Sinc 插值与求积，68 从稳定性反向设计格式，61 与 33 搬到均场情形。

### 与另一条多步路线的区别

这里需要区分两种都叫「多步」的构造，否则容易把结论互相错配。

- **插值再积分**（Zhao、Zhang 与 Ju，_SIAM J. Numer. Anal._ 48(4) 2010）：把参考积分恒等式中的被积函数用若干未来时间层作 Lagrange 插值，再积分，得到形如 $h\sum_j b_j\mathbb E_i[f(t_{i+j},Y_{i+j},Z_{i+j})]$ 的 Newton-Cotes 型权。这条路线的稳定窗口是**两个方向不同**的：$Y$ 的参考方程只在 $K_y\in\{1,\dots,7,9\}$ 时稳定（注意 $K_y=8$ 被排除在外），而 $Z$ 的参考方程只在 $K_z\in\{1,2,3\}$ 时稳定。
- **微分成参考 ODE**（编号 8）：**不**对被积函数插值再积分，而是把参考积分恒等式对 $t$ 求导，化成两个参考常微分方程，再用上面的导数逼近权离散。窗口因此是单一的 $1\le k\le6$。

用常微分方程数值的语言说：Zhao-Zhang-Ju 是 Adams 型（插值被积函数再积分），编号 8 是 BDF 型（差分导数）。两者都被编号 47 的统一框架覆盖。

顺带纠正一处常见的归属错误：多步（插值再积分）构造出自 Zhao-Zhang-Ju（2010），而 Zhao、Chen 与 Peng（_SIAM J. Sci. Comput._ 28(4) 2006）是 **$\theta$-格式**的来源，多步格式是对后者的推广。两者不宜混称。

## 18：带跳的正倒向系统

### 直觉

把编号 8 的机制搬到跳扩散上，障碍是组合性的：一个 $k$ 步格式跨越区间 $[t_n,t_{n+k}]$，而在这段时间里跳的次数原则上没有上限，逐一枚举跳的组合会使代价爆炸。论文用两条来化解：其一，跳扩散的生成元同样是**局部**的（现在是一个积分微分算子 $\mathcal L^0+\int_E[\cdot]\,\nu(\mathrm de)$），所以正向仍可以只用 Euler；其二，每个时间步只计入**一次跳**，把多跳事件当作更高阶的 $\Delta t$ 贡献丢掉。

### 问题设定

带跳的解耦正倒向系统的一般形式为

$$
X_t=X_0+\int_0^tb(s,X_{s^-})\,\mathrm ds+\int_0^t\sigma(s,X_{s^-})\,\mathrm dW_s
+\int_0^t\!\!\int_E\gamma(s,X_{s^-},e)\,\tilde\mu(\mathrm ds,\mathrm de),
$$

$$
Y_t=\varphi(X_T)+\int_t^Tf(s,X_s,Y_s,Z_s,\Gamma_s)\,\mathrm ds
-\int_t^TZ_s\,\mathrm dW_s-\int_t^T\!\!\int_EU_s(e)\,\tilde\mu(\mathrm ds,\mathrm de),
$$

$\tilde\mu$ 为补偿跳测度，$U$ 是与跳部分对应的额外未知量。这是 Barles-Buckdahn-Pardoux 的标准设定，广义 Feynman-Kac 公式把它与一个**偏积分微分方程**联系起来，非局部项来自跳测度。

### 可核实的结论

摘要逐字确认三点：（一）「受跳扩散过程局部性质的启发，用 Euler 方法求解相应的带跳正向随机微分方程，这大幅降低了整体计算复杂度，然而带跳倒向方程中关心的量仍保持高阶收敛率」；（二）「在每个时间步中，计算过程只涉及**一次跳**，这再次大幅降低了计算复杂度」；（三）「借助广义 Feynman-Kac 公式，方法容易应用于偏积分微分方程（以及某些非局部 PDE 模型）」。

> [!note] 可核实范围
> 论文正文付费且无预印本，本站只读到摘要与参考文献。上面的方程组按 Barles-Buckdahn-Pardoux 的标准设定写出，**论文实际使用的记号未核实**；多步模板中与跳分量 $U$ 对应的参考方程、「只涉及一次跳」在公式层面的确切实现、以及所声称高阶的具体阶数、假设与常数，本站均不报告。摘要只说「高阶收敛率」而未给出阶。数值实验方面，摘要提到「若干数值实验」，但测试问题与观察到的阶本站未核实。

### 与其他论文的关系

编号 8 的直系兄弟：同样的「正向用 Euler」技巧、同样的 $\alpha_{k,j}$ 多步模板，推广到跳的设定与非局部方程。它与编号 16、19（走向二阶 FBSDE）和编号 25（走向高维）是正交的三个方向。第一作者付宇也是编号 8、25、41 的作者。

## 23：用重复的便宜格式换一次昂贵格式

### 直觉

多步法加阶的办法是**造一个更宽的模板**：要 $k$ 阶就用 $k$ 个未来层，付出的是 $k$ 个起始值与一个只在 $k\le6$ 内成立的根条件。延迟校正走另一条路：**只用最便宜的 Euler 格式，但用它解一列残差方程**。关键观察是，若把已算出的低阶解的插值 $\mathcal I u$ 代回原方程，误差 $\delta=y-\mathcal Iu$ 满足一个与原方程**同形**的方程，因此可以用同一个 Euler 格式去解它；解完加回去，阶就升一级。重复 $J$ 次，阶升 $J$ 级（直到子网格能支撑的上限为止）。

这条路线的好处是结构上的：没有 $\alpha_{k,i}$ 模板，没有起始值问题，也就没有 $k\le6$ 的稳定性壁垒。论文自己的诊断是：「由于随机性的介入以及正向与倒向方程的耦合，很难为 FBSDE 设计高阶且相对『干净』的数值格式。」延迟校正的卖点正是「简单与稳健」。

### 问题设定与常微分方程情形的推导

先看常微分方程 $y'(t)=f(t,y(t))$、$y(0)=y_0$。取划分 $0=t_0<\cdots<t_N=T$，在每个 $I_n=[t_n,t_{n+1}]$ 内再加密为 $\mathcal G^n_K=\{t_{n,k}\}_{k=0}^K$，子步长 $\delta t=(t_{n+1}-t_n)/K$。给定低阶值 $\{u_{n,k}\}$ 及其连续插值 $\mathcal Iu(t)$，误差 $\delta(t)=y(t)-\mathcal Iu(t)$ 满足**残差方程**

$$
\delta'(t)=f\bigl(t,\delta(t)+\mathcal Iu(t)\bigr)-\frac{\mathrm d}{\mathrm dt}\mathcal Iu(t),
\qquad\delta(0)=0 .
$$

右端第一项与原方程同形，第二项是已知函数，因此**同一个低阶格式**可以直接解它，然后更新 $u_{n,k}^{\rm new}=u_{n,k}+\delta_k$。重复 $J$ 次给出（Hairer 1978）的收敛率

$$
O\bigl((\delta t)^{\min(J,K)+1}\bigr).
$$

$\min(J,K)$ 这个形式说明两件事：校正次数 $J$ 不能超过子网格分辨率 $K$ 所能支撑的阶，反之亦然。

### 推导：搬到 FBSDE 上

在 $I_n$ 上，从 $Y_t=Y_{t_{n+1}}+\int_t^{t_{n+1}}f(s,X_s,Y_s,Z_s)\mathrm ds-\int_t^{t_{n+1}}Z_s\mathrm dW_s$ 出发，定义误差过程 $\delta Y_t=Y_t-\mathcal I_hY_t$、$\delta Z_t=Z_t-\mathcal I_hZ_t$。它们满足**残差 BSDE**

$$
\delta Y_t=\delta Y_{t_{n+1}}+\int_t^{t_{n+1}}F(s,X_s,\delta Y_s,\delta Z_s)\,\mathrm ds
-\int_t^{t_{n+1}}\delta Z_s\,\mathrm dW_s+E(t),
$$

$$
F(s,X_s,\delta Y_s,\delta Z_s)=f\bigl(s,X_s,\delta Y_s+\mathcal I_hY_s,\ \delta Z_s+\mathcal I_hZ_s\bigr),
\qquad
E(t)=\mathcal I_hY_{t_{n+1}}-\int_t^{t_{n+1}}\mathcal I_hZ_s\,\mathrm dW_s-\mathcal I_hY_t .
$$

低阶格式取冻结系数的一步方法（正是 Euler，也正是编号 8 那一族的 $k=1$ 成员，论文明确这样说）：对 $k=K-1,\dots,0$，

$$
X^{k+1}=X^k+b(\tau_k,X^k)\delta t+\sigma(\tau_k,X^k)\Delta W_{\tau_k,\tau_{k+1}},
$$

$$
Z^k=\mathbb E^{X^k}_{\tau_k}\bigl[\bar Y^{k+1}(\Delta W_{\tau_k,\tau_{k+1}})^{\top}\bigr]\big/\delta t,
\qquad
Y^k=\mathbb E^{X^k}_{\tau_k}\bigl[\bar Y^{k+1}\bigr]+\delta t\cdot f(\tau_k,X^k,Y^k,Z^k).
$$

校正格式则是：置 $\delta Y^K=\delta Z^K=0$，对 $k=K-1,\dots,0$，

$$
\delta Z^k=\mathbb E^{X^k}_{\tau_k}\bigl[\delta\bar Y^{k+1}(\Delta W_{\tau_k,\tau_{k+1}})^{\top}\bigr]\big/\delta t
\;-\;Z^k\;+\;\nabla(\mathcal I_hY_{\tau_k})\,\sigma(\tau_k,X^k),
$$

$$
\delta Y^k=\mathbb E^{X^k}_{\tau_k}\bigl[\delta\bar Y^{k+1}\bigr]
+\delta t\Bigl(f\bigl(\tau_k,X^k,\delta Y^k+Y^k,\ \delta Z^k+Z^k\bigr)
+\mathcal L^0_{\tau_k,X^k}(\mathcal I_hY_{\tau_k})\Bigr).
$$

末尾那两个多出来的项来自残差 $E(t)$，它们由论文的引理 2.1（即编号 8 的生成元定理）给出：

$$
\left.\frac{\mathrm d\,\mathbb E^{X^k}_{\tau_k}[\mathcal I_hY_t]}{\mathrm dt}\right|_{t=\tau_k}
=\mathcal L^0_{\tau_k,X^k}(\mathcal I_hY_{\tau_k}),
\qquad
\left.\frac{\mathrm d\,\mathbb E^{X^k}_{\tau_k}[\mathcal I_hY_t(\Delta W_{\tau_k,\tau_{k+1}})^{\top}]}{\mathrm dt}\right|_{t=\tau_k}
=\nabla(\mathcal I_hY_{\tau_k})\,\sigma(\tau_k,X^k).
$$

> [!warning] 延迟校正的真实代价在插值算子上
> 论文自己点明：「延迟校正格式的高阶精度严重依赖 $\partial(\mathcal I_hY_t)/\partial t$、$\partial(\mathcal I_hY_t)/\partial x$ 与 $\partial^2(\mathcal I_hY_t)/\partial x^2$ 的逼近质量」，因为 $\mathcal L^0$ 里含一个二阶空间导数。也就是说，多步法把负担放在时间模板上，延迟校正把负担放在**插值算子必须可微两次且微分后仍准确**上。这不是免费的简化，而是代价的转移。

整体算法：给定 $Y^N_i,Z^N_i$；对 $n=N-1,\dots,0$，置 $Y^{n,K}_i=Y^{n+1}_i$、$Z^{n,K}_i=Z^{n+1}_i$；对 $j=1,\dots,J$ 依次（一）用低阶格式从 $k=K-1$ 倒推得 $Y^{n,k,[j]}_i,Z^{n,k,[j]}_i$，（二）置 $\delta Y^{K,[j]}_i=\delta Z^{K,[j]}_i=0$ 并用同一低阶格式解 $\delta Y^{k,[j]}_i,\delta Z^{k,[j]}_i$，（三）更新 $Y^{n,k,[j+1]}_i=Y^{n,k,[j]}_i+\delta Y^{k,[j]}_i$、$Z^{n,k,[j+1]}_i=Z^{n,k,[j]}_i+\delta Z^{k,[j]}_i$。最后取 $Y^n_i=Y^{n,0,[J]}_i$、$Z^n_i=Z^{n,0,[J]}_i$。空间框架与编号 8 相同：网格 $\mathcal D_h=\{x_i\}$、密度 $h=\max_x\mathrm{dist}(x,\mathcal D_h)$、邻点集 $\#\mathcal D_{h,x}\le N_e$。

### 定理

- 常微分方程情形的率 $O((\delta t)^{\min(J,K)+1})$ 是转述 Hairer 的已知结果。
- 对 FBSDE，论文**声称** $K$ 阶并用数值验证；本站在正文中**未找到**带显式常数的 FBSDE 收敛定理。
- 稳定性是数值断言：「稳定且有效，至少对 $K=1,\dots,4$ 是 $K$ 阶方法」。关键在于延迟校正**没有 $k\le6$ 壁垒**：测试 2 一直做到 $K=12$。

### 数值实验

**测试 1，解耦系统 (5.1)。** $\mathrm dX_t=\frac{1}{1+2e^{t+X_t}}\mathrm dt+\frac{e^{t+X_t}}{1+e^{t+X_t}}\mathrm dW_t$，$-\mathrm dY_t=\bigl(-\frac{2Y_t}{1+2e^{t+X_t}}-\frac12(\frac{Y_tZ_t}{1+e^{t+X_t}}-Y_t^2Z_t)\bigr)\mathrm dt-Z_t\mathrm dW_t$，$Y_T=\frac{e^{T+X_T}}{1+e^{T+X_T}}$，$x=1$；精确解 $Y_t=\frac{e^{t+X_t}}{1+e^{t+X_t}}$、$Z_t=\frac{(e^{t+X_t})^2}{(1+e^{t+X_t})^3}$。$N=4,6,8,10,12$ 上的收敛率：

| $K$        | 1     | 2     | 3     | 4     |
| ---------- | ----- | ----- | ----- | ----- |
| $Y$ 收敛率 | 0.995 | 1.993 | 3.109 | 4.097 |
| $Z$ 收敛率 | 0.994 | 1.980 | 2.982 | 4.024 |

**注意 $N$ 有多小。** 这些是 4 到 12 个大时间步，每步内再分 $K$ 个子步——延迟校正把精度从「时间步数」搬到「每步的校正次数」上。

**测试 1，耦合系统 (5.3)。** $\mathrm dX_t=\frac{1}{1+e^{t+X_t}}\cdot\frac{1}{1+Y_t}\mathrm dt+Y_t\,\mathrm dW_t$，倒向方程与终值同上，$x=0$。$K=1$ 的收敛率为 $0.937/1.008$；更高的 $K$ 沿同样的 $K$ 阶模式（$K=2$ 的误差自 $3.666\times10^{-4}$ 起），本站只部分转录，故不列表。

**测试 2，系统 (5.4)。** $\mathrm dX_t=\sin(t+X_t)\mathrm dt+\frac{3}{10}\cos(t+X_t)\mathrm dW_t$，$-\mathrm dY_t=\bigl(\frac{3}{20}Y_tZ_t-\cos(t+X_t)(1+Y_t)\bigr)\mathrm dt-Z_t\mathrm dW_t$，$Y_T=\sin(T+X_T)$，$x=0.5$；精确解 $Y_t=\sin(t+X_t)$、$Z_t=\frac{3}{10}\cos^2(t+X_t)$。论文的结论是：「我们的延迟校正方法是 $K$ 阶方法（$K=1,2,\dots,12$），稳定、有效，收敛率非常高（可达 12）。」抽样的率为 $K=1$：$0.906/1.001$；$K=2$：$2.193/2.031$。效率上的说法是：$N=4$ 且 $K=2$ 得到的误差**远小于** $N=12$ 且 $K=1$（即 Euler 格式）得到的误差。

**这些实验建立了什么，又差在哪里。** 建立的是：延迟校正在一维光滑问题上确实按 $K$ 阶收敛，包括耦合情形，而且**做到了 $K=12$——远超多步族的 $k\le6$ 壁垒**，这是这条路线相对多步法最实在的证据。差在两处。第一，这一切都是数值观察，论文没有给出 FBSDE 的收敛定理，因此「为什么没有壁垒」在这一篇里仍是空白（编号 47 的框架给出了讨论它的语言，但本站未核实其对延迟校正族的具体结论）。第二，本站只转录到部分表格，耦合测试与高 $K$ 的完整误差表未纳入，因此上表不足以支持任何跨 $K$ 的效率比较；论文自己给出的效率论断只有 $N=4,K=2$ 对 $N=12,K=1$ 这一处。

### 与其他论文的关系

它是编号 8 的**方法论替代**：同样的「正向用 Euler」技巧与同样的空间框架，但阶来自延迟校正迭代而非 $k$ 步模板。$\min(J,K)+1$ 的率结构与没有 $k\le6$ 壁垒这两点，是延迟校正区别于多步族的地方，也是编号 68 后来做强稳定保持工作的动机之一。编号 35 是它的直接续篇：把延迟校正用到编号 19 的二阶方程类上，并做成显式。

## 33：均场倒向方程的显式 $\theta$ 格式

### 直觉

均场（McKean-Vlasov 型）倒向方程的生成元依赖解自身的**分布**，因此每算一次生成元都要对另一层期望求值。这一层期望是新的代价来源，也是把已有格式搬过来时的主要障碍：任何隐式处理都会把「解一个非线性方程」变成「解一个含未知分布的非线性方程」。「显式」$\theta$ 格式的用意就在这里——把均场生成元放在**已知的未来层**上求值，于是每步不需要非线性求解。

### 问题设定

均场倒向随机微分方程的标准形式（Buckdahn-Djehiche-Li-Peng 与 Buckdahn-Li-Peng）是

$$
Y_t=\xi+\int_t^T \mathbb E'\bigl[f(s,X'_s,Y'_s,Z'_s,X_s,Y_s,Z_s)\bigr]\mathrm ds
-\int_t^T Z_s\,\mathrm dW_s,
$$

其中 $\mathbb E'$ 是对 $(X,Y,Z)$ 的独立副本 $(X',Y',Z')$ 取期望。

### 可核实的结论

摘要给出三条，本站逐字确认：其一，所提出的是一族**显式 $\theta$ 格式**；其二，「我们首先证明一个严格的**稳定性结果**」，然后「在此基础上给出**尖锐的误差估计**，表明所提出的 $\theta$ 格式具有**二阶收敛率**」；其三，「这似乎是**首次尝试为均场倒向随机微分方程设计高阶数值格式**」——此前的均场格式都是一阶。

值得单独指出的是分析的**架构**：先稳定性，再误差估计。这与编号 47 后来抽象成一般框架的正是同一套次序，也与编号 63 的做法一致。三篇都是 SIAM J. Numer. Anal. 上赵卫东与周涛的文章。

> [!note] 可核实范围
> 论文无预印本，SIAM 阻止本环境的自动访问，因此除摘要外的一切本站未核实。具体而言：$\theta$ 的确切放置位置、分布或期望 $\mathbb E'$ 的离散方式（粒子系统？求积？对高斯的嵌套期望？）、$Z$ 的处理、稳定性不等式的确切形式、误差估计的假设与常数、哪些 $\theta$ 值给出二阶、以及该二阶是仅对 $Y$ 还是对 $(Y,Z)$ 同时成立——这些本站都不报告。数值方面，摘要说「进行了若干数值实验以验证理论结果」，测试问题与观察到的阶未核实。

### 与其他论文的关系

它开启了这个研究计划的均场分支；编号 61 是其显式续篇（自述为「我们关于均场正倒向随机微分方程数值方法系列工作之一」），从均场 **B**SDE 扩到均场 **FB**SDE，从 $\theta$ 格式扩到多步格式，并直接引用编号 33。$\theta$ 格式族后来被编号 47 证明是其统一格式的特例。编号 97（DeepSPoC）与 108 从深度学习与混沌传播的完全不同方向攻同类问题。

## 35：把延迟校正用到二阶方程上

### 直觉

编号 23 的校正步是**隐式**的（$\delta Y^k$ 出现在自身方程的右端），因此每次校正都要解一个非线性方程。若能把校正做成显式，则每轮迭代的代价就与一次 Euler 扫描相当，而阶仍逐轮上升。这正是这一篇宣称的东西，而它作用的对象是编号 16 与 19 的二阶（全非线性）方程类。

### 问题设定

目标方程类与编号 19 相同（见[[computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control|二阶 FBSDE 一页]]）：$\mathrm dX=b\,\mathrm dt+\sigma\,\mathrm dW$、$-\mathrm dY=f(t,X,Y,Z,\Gamma)\mathrm dt-Z\,\mathrm dW$、$\mathrm dZ=A\,\mathrm dt+\Gamma\,\mathrm dW$，其解表示一个全非线性抛物方程。

### 可核实的结论

摘要逐字给出四点：（一）「这是我们关于正倒向随机微分方程延迟校正方法系列论文的第二部分」，直接续编号 23；（二）「我们提出一类用于 2FBSDE 的**显式**延迟校正格式」——「显式」是与编号 23 的关键区别；（三）「关键特征是用**简单的 Euler 格式作为初始化**」，随后「通过一个简单的延迟校正迭代格式，就能得到**精度非常高**的近似解」；（四）「而在每次迭代中，**计算复杂度始终与 Euler 求解器相当**」。

论文还给出一句范围上的自陈，值得原样保留：「我们相信本文提出的格式在处理**中等维数**的 2FBSDE 时是有前景的」——这是对方法不能达到高维的明确承认。

> [!note] 可核实范围
> 论文付费且无预印本，本站只读到摘要与参考文献。残差 2FBSDE 的确切形式、校正扫描次数 $J$ 与所达阶的关系、$\Gamma$ 与 $A$ 在校正器中的处理方式，本站均未核实。摘要声称「非常高的精度」与每轮常数代价，但**未给出阶、假设或常数**，本站也未读到任何定理。数值方面，摘要说「给出数值例子以说明所提格式的有效性」，测试问题与观察到的阶未核实。

### 与其他论文的关系

它位于编号 23（延迟校正方法论）与编号 19、16（二阶方程类）的交点。第一作者杨杰同时是编号 47 与 61 的作者；编号 47 提供的一般稳定性与相容性框架，正是能同时容纳多步族与延迟校正／$\theta$ 族的语言。「每次迭代复杂度与 Euler 求解器相当」这条声明，是延迟校正族对多步族的结构性卖点——后者的 $k$ 步格式每层大约花 $k$ 倍 Euler 步的代价，还要 $k$ 个起始值。

## 61：均场正倒向系统的显式多步格式

### 直觉

编号 33 用 $\theta$ 格式把均场 BSDE 推到二阶。要再往上，自然的做法是换成多步模板——这正是编号 8 在非均场情形做的事。难点仍然是那一层对分布的期望 $\mathbb E'$：它必须被离散成某种可计算的东西，且这一离散不能破坏时间方向的阶。

### 问题设定

McKean-Vlasov 型的正倒向系统为

$$
X_t=X_0+\int_0^t\mathbb E'\bigl[b(s,X'_s,X_s)\bigr]\mathrm ds+\int_0^t\mathbb E'\bigl[\sigma(s,X'_s,X_s)\bigr]\mathrm dW_s,
$$

$$
Y_t=\mathbb E'\bigl[\varphi(X'_T,X_T)\bigr]+\int_t^T\mathbb E'\bigl[f(s,X'_s,Y'_s,Z'_s,X_s,Y_s,Z_s)\bigr]\mathrm ds-\int_t^TZ_s\,\mathrm dW_s .
$$

### 可核实的结论

摘要与关键词逐字确认：格式是**显式多步**的、「易于实现」、「具有高阶收敛率」，且「给出了所提多步格式的**严格误差估计**」。关键词为「均场正倒向随机微分方程；显式多步格式；误差估计」，MSC 为 60H35、65C20、60H10。

参考文献的构成也是可核实的证据，它说明这篇站在哪些工作上：Zhao-Chen-Peng（SISC 28 (2006) 1563-1581，Lagrange 插值多步思想的来源）、Zhao-Zhang-Ju（SINUM 48 (2010) 1369-1394，稳定多步格式；以及 NMTMA 9 (2016) 262-288，解耦 FBSDE 的多步格式）、编号 8、33、35、47、25，Sun-Yang-Zhao（NMTMA 10 (2017) 798-828，均场 SDE 的 Itô-Taylor 格式）、Kloeden-Shardlow（SISC 39 (2017) A2784-A2807，一维均场 SDE 的高斯求积方法），以及 Buckdahn-Djehiche-Li-Peng 与 Buckdahn-Li-Peng。**Kloeden-Shardlow 出现在参考文献里，使得对 $\mathbb E'$ 采用求积处理成为一种合理猜测，但这只是猜测，本站不作断言。**

> [!note] 可核实范围
> 正文付费且无预印本，本站只读到落地页（含完整 38 项参考文献）、OpenAlex 与 Semantic Scholar 记录。多步系数、所用的步数 $k$、以及任何均场格式最关键的一环——$\mathbb E'$ 究竟如何离散——本站均未核实。所声称高阶的具体阶数、范数、假设，以及是否出现类似编号 8 的 $1\le k\le6$ 根条件窗口，也未核实。数值方面，摘要说「进行了数值实验以说明所提格式的效率与精度」，测试问题与观察到的阶未核实。
>
> 另有一处著录更正：常见的条目把页码写作 40 (2022) 519-543，而 Global Science Press 落地页与 OpenAlex 都给出 **517-540**，本页采用后者。

### 与其他论文的关系

编号 33 的直接续篇（并引用它）：33 对均场 BSDE 用 $\theta$ 格式达到二阶，61 对均场 FBSDE 用多步格式声称高阶。方法论上的父辈是编号 8。第二作者杨杰与编号 35、47 相同，而编号 47 的统一框架被本文引用，这提示其误差分析大概率沿用同一套「先稳定性」的架构——同样只是提示。

## 68：一旦稳定性成为枢纽，就设计格式去最大化它

### 直觉

到 2022 年，这一组手里有两样东西：一族高阶多步格式，其稳定性只由经验观察到的根条件窗口刻画（例如编号 8 的 $1\le k\le6$）；以及一个一般的稳定性、相容性与收敛性框架（编号 47），它证明了均方 Lax 等价定理，却没有告诉人如何**构造**稳定性好的格式。把两者接起来的自然一步是：既然稳定性是收敛的充要一半，那就把稳定性当作**设计目标**，在满足相容阶的系数集合里挑稳定性最好的那些。这正是双曲守恒律数值中强稳定保持方法的思路——那里人们最大化使凸组合（收缩性）性质得以保留的类 CFL 系数。

### 可核实的结论

摘要的表述是：作者先对一般类型的 FBSDE 多步格式做全面分析，据此给出**关于系数的新充分条件**使相应格式稳定且具有一定的相容阶，再据此提出**构造高阶强稳定保持多步格式的实用方法**。Springer 页面上的附录标题也可核实：附录 A「Additional Optimal SSPM Schemes」含表 7「均匀时间划分下的 SSPM 格式（第 2 部分）」，给出「**阶数至多 5** 的最优 SSPM 格式的系数」。所以论文的实用产出是一张一到五阶的最优系数表。

被分析的格式模板与编号 47 统一的一族同形：对 $k$ 步方法与 $\mathbb E_n[\cdot]=\mathbb E[\cdot\mid\mathcal F_{t_n}]$，

$$
\sum_{i=0}^{k}\alpha_i\,\mathbb E_n\bigl[Y^{n+i}\bigr]
=\Delta t\sum_{i=0}^{k}\beta_i\,\mathbb E_n\bigl[f(t_{n+i},X^{n+i},Y^{n+i},Z^{n+i})\bigr],
$$

配以由 $\mathbb E_n[Y^{n+i}\Delta W]/\Delta t$ 项构成的 $Z^n$ 递推。强稳定保持的问题因此是：对哪些系数向量 $(\alpha_i),(\beta_i)$，格式可以改写成若干类后向 Euler 步的**凸组合**，从而基础步的任何单调性或收缩性被继承，并且步长系数尽可能大。这与 Lenferink（_Numer. Math._ 55 (1989) 213-223；_Math. Comp._ 56 (1991) 177-199）以及 Spiteri 与 Ruuth（_SIAM J. Numer. Anal._ 40 (2002) 469-491）的强稳定保持线性多步理论完全对应，那三篇都在本文参考文献中。

> [!note] 可核实范围
> 上述定位、摘要中的三步表述与附录中一至五阶系数表的**存在**均可确认。上面的格式模板是按编号 47 与编号 8 的同族格式反推的，**并非转录自本文**。充分条件的确切形式、FBSDE 意义下「强稳定」的精确定义（这是本文最关键的技术定义，从摘要无法恢复）、所用的稳定性泛函、得到「最优」系数的优化问题、以及系数的数值，本站均未核实；空间离散与条件期望的处理方式也无法从摘要推出。
>
> 数值方面，摘要说「进行了数值实验以展示我们 SSPM 格式的**强稳定性**」，Springer 页面显示至少有一张图。测试问题与观察到的阶本站未核实；从摘要的侧重看，这些实验更可能是稳定性演示而非阶数验证表，但这也只是推测。

### 与其他论文的关系

编号 8 与 47 的直接后继，两者都在其参考文献中。它同时引用编号 63（同组另一篇 2022 年的 SINUM 文章），也引用 Tang 与 Xiong（_IMA J. Numer. Anal._ 42 (2022) 1789-1805，Markov 型 BSDE 一般多步方法的稳定性分析）与 Chassagneux（_SIAM J. Numer. Anal._ 52 (2014) 2815-2836，BSDE 的线性多步格式）——这两篇外部工作也是编号 47 对话的对象。第一作者方水新在此首次出现，他后来成为深度学习几篇（编号 86、93、96、100、108）的主要合作者；**这篇是他的经典数值分析入口，而鞅型深度学习几篇继承了他的稳定性直觉。**

## 七篇的关系

| 编号 | 加阶手段           | 问题类型       | 稳定性依据         | 本站核实程度   |
| ---- | ------------------ | -------------- | ------------------ | -------------- |
| 8    | 微分成参考 ODE     | 耦合 FBSDE     | 根条件窗口 $k\le6$ | 全文逐式       |
| 18   | 多步 + 跳处理      | 带跳 FBSDE     | 未核实             | 摘要与参考文献 |
| 23   | 延迟校正           | FBSDE          | 数值断言，无壁垒   | 全文逐式       |
| 33   | 显式 $\theta$ 格式 | 均场 BSDE      | 论文自证（未核实） | 摘要           |
| 35   | 显式延迟校正       | 二阶 FBSDE     | 未核实             | 摘要与参考文献 |
| 61   | 显式多步           | 均场 FBSDE     | 未核实             | 摘要与参考文献 |
| 68   | 按稳定性反向设计   | 一般多步 FBSDE | 新充分条件 + 优化  | 摘要与附录标题 |

这条线索的形状值得总结：**先构造格式（编号 8 至 35），再统一分析（编号 47），最后按分析结果反向设计（编号 68）。** 编号 68 之所以可能，正因为编号 47 已把稳定性确立为收敛的充要一半，因此「设计稳定的格式」成为一个有明确目标的优化问题，而不是试错。

另一条形状是两种加阶手段的分工。多步法把负担放在**时间模板**上：$k$ 个起始值、$k$ 个未来层、一个根条件窗口。延迟校正把负担放在**插值算子**上：它必须光滑到可以微分两次，因为残差方程里出现 $\mathcal L^0(\mathcal I_hY)$。两者的实测上限也不同——编号 8 到 $k=6$（双精度下实际更低），编号 23 到 $K=12$。

## 覆盖核对

| 内容                                       | 论文 | 覆盖状态                                 |
| ------------------------------------------ | ---- | ---------------------------------------- |
| 耦合 FBSDE、既有状况与论文提出的问题       | 8    | 完整                                     |
| 生成元定理、导数权、两条参考 ODE、五个格式 | 8    | 完整推导                                 |
| 截断误差分解、误差平衡、根条件窗口         | 8    | 完整；并注明缺收敛定理                   |
| 三组数值实验与其局限                       | 8    | 例 1 至例 3 的收敛率表；例 4 未转录      |
| 带跳系统与三条摘要论断                     | 18   | 仅摘要可支撑的内容                       |
| 延迟校正：残差方程、低阶格式、校正格式     | 23   | 完整推导，含两个生成元恒等式             |
| 延迟校正的两组实验与 $K=12$                | 23   | 测试 1 解耦完整；耦合与测试 2 部分转录   |
| 均场 BSDE 形式、「首次高阶」与二阶声明     | 33   | 仅摘要可支撑的内容                       |
| 显式延迟校正的四条摘要论断与维数自陈       | 35   | 仅摘要可支撑的内容                       |
| 均场 FBSDE 形式、误差估计声明与参考构成    | 61   | 仅摘要与参考文献可支撑的内容；含页码更正 |
| 强稳定保持的设计思想、模板与系数表存在性   | 68   | 仅摘要与附录标题可支撑的内容             |

## 本页原文

- W. Zhao, Y. Fu, and T. Zhou, [_New kinds of high-order multistep schemes for coupled forward backward stochastic differential equations_](https://doi.org/10.1137/130941274), SIAM J. Sci. Comput. 36(4) (2014), pp. A1731-A1751（预印本 [arXiv:1310.5307](https://arxiv.org/abs/1310.5307)）。
- Y. Fu, W. Zhao, and T. Zhou, [_Multistep schemes for forward backward stochastic differential equations with jumps_](https://doi.org/10.1007/s10915-016-0212-y), J. Sci. Comput. 69(2) (2016), pp. 651-672。
- T. Tang, W. Zhao, and T. Zhou, [_Deferred correction methods for forward backward stochastic differential equations_](https://doi.org/10.4208/nmtma.2017.s02), Numer. Math. Theor. Meth. Appl. 10(2) (2017), pp. 222-242。
- Y. Sun, W. Zhao, and T. Zhou, [_Explicit theta-schemes for mean-field backward stochastic differential equations_](https://doi.org/10.1137/17M1161944), SIAM J. Numer. Anal. 56(4) (2018), pp. 2672-2697。
- J. Yang, W. Zhao, and T. Zhou, [_Explicit deferred correction methods for second-order forward backward stochastic differential equations_](https://doi.org/10.1007/s10915-018-00896-w), J. Sci. Comput. 79(3) (2019), pp. 1409-1432。
- Y. Sun, J. Yang, W. Zhao, and T. Zhou, [_An explicit multistep scheme for mean-field forward-backward stochastic differential equations_](https://doi.org/10.4208/jcm.2011-m2019-0205), J. Comput. Math. 40(4) (2022), pp. 517-540。
- S. Fang, W. Zhao, and T. Zhou, [_Strong stability preserving multistep schemes for forward backward stochastic differential equations_](https://doi.org/10.1007/s10915-023-02111-x), J. Sci. Comput. 94(3) (2023), 53。
- 用于交叉印证的外部来源：W. Zhao, G. Zhang, and L. Ju, [_A stable multistep scheme for solving backward stochastic differential equations_](https://doi.org/10.1137/09076979X), SIAM J. Numer. Anal. 48(4) (2010), pp. 1369-1394；W. Zhao, L. Chen, and S. Peng, [_A new kind of accurate numerical method for backward stochastic differential equations_](https://doi.org/10.1137/05063341X), SIAM J. Sci. Comput. 28(4) (2006), pp. 1563-1581；L. Teng, A. Lapitckii, and M. Günther, _A multi-step scheme based on cubic spline for solving backward stochastic differential equations_, [arXiv:1809.00324](https://arxiv.org/abs/1809.00324)（本站据此转录 Zhao-Zhang-Ju 的参考方程与稳定窗口）。
