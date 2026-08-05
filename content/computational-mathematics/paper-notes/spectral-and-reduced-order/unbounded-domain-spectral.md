---
title: 无界区域与分数阶算子
description: 编号 27、42、92：让基函数的远场行为与分数阶算子的作用相容
lang: zh
translation: en/computational-mathematics/paper-notes/spectral-and-reduced-order/unbounded-domain-spectral
tags:
  - 论文笔记
  - 谱方法
  - 分数阶算子
---

> [!note] 本页覆盖
> 编号 **27**（_Commun. Comput. Phys._ 24(4), 2018）、**42**（_SIAM J. Sci. Comput._ 42(2), 2020）、**92**（_SIAM J. Sci. Comput._ 48(1), 2026）。编号 92 的正文未能通过公开渠道获取，该节只报告摘要、关键词与参考文献可确认的内容。

模型问题是

$$
(-\Delta)^{\alpha/2}u(x)+\rho\,u(x)=f(x),\quad x\in\mathbb R^d,
\qquad u(x)\to0,\ |x|\to\infty,
$$

其中 $\rho>0$、$\alpha\in(0,2)$。分数阶 Laplacian 有两个等价定义，超奇异积分

$$
(-\Delta)^{\alpha/2}u(x)=C_{d,\alpha}\,\mathrm{p.v.}\!\int_{\mathbb R^d}
\frac{u(x)-u(y)}{|x-y|^{\,d+\alpha}}\,\mathrm dy,
\qquad
C_{d,\alpha}=\alpha\,2^{\alpha-1}\,
\frac{\Gamma\bigl(\tfrac{d+\alpha}{2}\bigr)}{\pi^{d/2}\,\Gamma\bigl(\tfrac{2-\alpha}{2}\bigr)},
$$

与拟微分算子形式 $(-\Delta)^{\alpha/2}u=\mathcal F^{-1}[|\xi|^{\alpha}\mathcal F[u]]$。

## 无界区域的三条标准路线都不适用

编号 42 把困难说得很清楚：分数阶算子的解在无穷远**只是代数衰减**，服从幂律。这一条使三条标准策略同时失效。

区域截断适合快速衰减的解，这里不适合；而且朴素截断会在无界区域被终止的界面处引入非物理奇性。透明边界条件与吸收层因算子非局部而变得极难构造。无界区域上的正交函数——Hermite 与 Laguerre——是为**指数衰减**调好的，与幂律尾部不匹配。

映射 Jacobi（有理）基对**整数阶**且解代数衰减的问题早已知道优于 Hermite 与 Laguerre。但把映射技巧推广到分数阶远非平凡，原因是结构性的：普通 Laplacian 把有理基函数映到同类函数，而分数阶 Laplacian 把它映到**性质完全不同的一类函数**。

## 27：Hermite 基与显式微分矩阵

有界区域上成功的做法是选一个与分数阶算子相容的基——Jacobi 多重分数阶函数（一个分数阶 Sturm-Liouville 问题的特征函数）与广义 Jacobi 函数——它们的分数阶导数仍是同族中参数移位后的函数，因此分数阶算子实际上变成局部的。无界区域上没有这样的装置。

同时期的替代方案（Mao 与 Shen 的 Hermite 谱方法）转到频率空间做配点，因此每次求值都要一次正向与一次反向 Hermite 变换，对非线性问题很不方便。编号 27 的贡献是在**物理空间**直接配点，依据是显式闭式微分矩阵。

论文给出两套 Hermite 基。归一化 Hermite 函数

$$
\widehat H_n(x)=\frac{1}{\sqrt{2^n n!}}\,e^{-x^2/2}H_n(x),
\qquad
\int_{\mathbb R}\widehat H_m\widehat H_n\,\mathrm dx=\sqrt{\pi}\,\delta_{mn}
$$

关于权 $\omega\equiv1$ 正交；过缩放（Brinkman）基

$$
\widetilde H_n(x)=e^{-x^2/2}\widehat H_n(x)=\frac{1}{\sqrt{2^n n!}}e^{-x^2}H_n(x),
\qquad
\int_{\mathbb R}\widetilde H_m\widetilde H_n\,e^{x^2}\mathrm dx=\sqrt{\pi}\,\delta_{mn}
$$

由 Brinkman 为 Fokker-Planck 方程引入，其中速度方向的概率密度按这组函数展开，此后成为该方程的标准基之一。两套基的差别在于把 $e^{-x^2/2}$ 放在基里还是放在权里，而这决定了离散内积与微分矩阵的形式。

## 42：把有理基分解成两种初等形状

### 映射与基

一对一映射与其恒等式为

$$
x=\frac{t}{\sqrt{1-t^2}}
\ \Longleftrightarrow\
t=\frac{x}{\sqrt{1+x^2}},
\qquad
1-t^2=\frac{1}{1+x^2},
\qquad
\frac{\mathrm dx}{\mathrm dt}=\frac{1}{(1-t^2)^{3/2}} .
$$

这是**奇异**的（无界区间）代数映射，正是它产生基函数的幂律尾部。修正映射 Gegenbauer 函数定义为

$$
R^{\lambda}_n(x):=\bigl(1+x^2\bigr)^{-\frac{\lambda+1}{2}}
C^{\lambda}_n\!\Bigl(\frac{x}{\sqrt{1+x^2}}\Bigr),
\qquad x\in\mathbb R,\ \lambda>-\tfrac12 .
$$

「修正」指的是基**把权函数的平方根吸收进来**：写成 $R^{\lambda}_n(x)=S(t)C^{\lambda}_n(t)$，其中

$$
S(t)=\sqrt{\omega_{\lambda}(t)\frac{\mathrm dt}{\mathrm dx}}=(1-t^2)^{\frac{\lambda+1}{2}} ,
$$

因此这组函数关于 $\mathbb R$ 上的**均匀**权正交：

$$
\int_{-\infty}^{\infty}R^{\lambda}_n(x)R^{\lambda}_m(x)\,\mathrm dx=\gamma^{\lambda}_n\delta_{nm},
\qquad
\gamma^{\lambda}_n=\frac{\pi\,2^{1-2\lambda}\Gamma(n+2\lambda)}{n!\,(n+\lambda)\,\Gamma^2(\lambda)} .
$$

基自身满足三项递推

$$
n\,R^{\lambda}_n(x)=\frac{2x}{\sqrt{1+x^2}}(n+\lambda-1)R^{\lambda}_{n-1}(x)
-(n+2\lambda-2)R^{\lambda}_{n-2}(x),\qquad n\ge2,
$$

$$
R^{\lambda}_0(x)=\frac{1}{(1+x^2)^{\frac{\lambda+1}{2}}},
\qquad
R^{\lambda}_1(x)=\frac{2\lambda x}{(1+x^2)^{1+\frac{\lambda}{2}}} .
$$

远场行为是关键：

$$
\lim_{x\to+\infty}(1+x^2)^{\frac{\lambda+1}{2}}R^{\lambda}_n(x)=\frac{(2\lambda)_n}{n!},
$$

因此每个基函数都按 $|x|^{-(\lambda+1)}$ 衰减，$\lambda$ 就是把基调到解的幂律的旋钮。

### 两种初等形状

用 Pfaff 变换把基展开后，**每个基函数都是两种初等形状的有限线性组合**：

$$
\frac{1}{(1+x^2)^{\gamma}}\ \bigl(\gamma=k+\tfrac{\lambda+1}{2}\bigr),
\qquad
\frac{x}{(1+x^2)^{\gamma}}\ \bigl(\gamma=k+\tfrac{\lambda}{2}+1\bigr).
$$

于是只需对这两种形状求分数阶 Laplacian。对实 $s>0$：

$$
(-\Delta)^{s}\Bigl\{\frac{1}{(1+x^2)^{\gamma}}\Bigr\}
=A^{\gamma}_s\ {}_2F_1\!\Bigl(s+\gamma,\ s+\tfrac12;\ \tfrac12;\ -x^2\Bigr),
\qquad \gamma>0,
$$

$$
(-\Delta)^{s}\Bigl\{\frac{x}{(1+x^2)^{\gamma}}\Bigr\}
=(2s+1)A^{\gamma}_s\,x\ {}_2F_1\!\Bigl(s+\gamma,\ s+\tfrac32;\ \tfrac32;\ -x^2\Bigr),
\qquad \gamma>\tfrac12,
$$

$$
A^{\gamma}_s=\frac{2^{2s}\,\Gamma(s+\gamma)\,\Gamma\bigl(s+\tfrac12\bigr)}{\sqrt{\pi}\,\Gamma(\gamma)} .
$$

证明路径值得记录，因为它说明了修正 Bessel 函数从哪里进来。先把 $v(x)=(1+x^2)^{-\gamma}$ 写成 ${}_2F_1(\gamma,\tfrac12;\tfrac12;-x^2)$，用余弦变换恒等式得到第二类修正 Bessel 函数

$$
\mathcal F[v](\xi)=\frac{2^{1-\gamma}}{\Gamma(\gamma)}\,\xi^{\gamma-\frac12}
K_{\gamma-\frac12}(\xi),\qquad \xi>0,
$$

再由 $(-\Delta)^s v=\mathcal F^{-1}[|\xi|^{2s}\mathcal F[v]]$ 化为一个 Bessel 积分，用相应的积分恒等式完成。奇次形状经 $\mathcal F[xv](\xi)=-i\frac{2^{1-\gamma}}{\Gamma(\gamma)}\xi^{\gamma-\frac12}K_{\gamma-\frac32}(\xi)$ 同法处理。

### 整数阶与分数阶的分岔

这是全篇最有信息量的结构观察。若 $s$ 为正整数，两个 ${}_2F_1$ 终止，算子**增加 $s$ 阶衰减**：

$$
(-\Delta)^{s}\Bigl\{\tfrac{1}{(1+x^2)^{\gamma}}\Bigr\}\sim\frac{1}{(1+x^2)^{s+\gamma}} .
$$

若 $s$ 非整数，超几何函数可能在 $|x|\to\infty$ 时发散，图景改变。取 $\gamma>1/2$ 时

$$
(-\Delta)^{s}\bigl\{(1+x^2)^{-\gamma}\bigr\}\sim(1+x^2)^{-(s+1/2)},
$$

**衰减率与 $\gamma$ 无关**；$\gamma=1/2$ 时出现对数因子 $\ln(1+x^2)/(1+x^2)^{s+1/2}$；只有 $0<\gamma<1/2$ 保留整数阶的行为。在基的层面上因此有

$$
(-\Delta)^{s}R^{\lambda}_{2n}(x)\sim
\begin{cases}
(1+x^2)^{-\left(s+\frac{\lambda+1}{2}\right)}, & -\tfrac12<\lambda<0,\\[3pt]
\ln(1+x^2)\,(1+x^2)^{-\left(s+\frac12\right)}, & \lambda=0,\\[3pt]
(1+x^2)^{-\left(s+\frac12\right)}, & \lambda>0 .
\end{cases}
$$

也就是说，分数阶 Laplacian 作用在这组基上**并不总是带来 $1/(1+x^2)^s$ 的衰减增益**。这正是引言中「性质完全不同」这句话的定量版本。

### 两个可实现性细节

分数阶 Laplacian 作用在基上的公式不如 Hermite 情形紧凑，但可用递推高效求值。令 $F_k(x)={}_2F_1(a,b;c;-x^2)$，$a=s+k+\tfrac{\lambda+1}{2}$、$b=s+\tfrac12$、$c=\tfrac12$，则

$$
F_{k+1}(x)=\frac{c-a}{a(1+x^2)}F_{k-1}(x)
+\frac{(2a-c)+(a-b)x^2}{a(1+x^2)}F_k(x),\qquad k\ge1 .
$$

另有一个尺度参数 $\mu>0$ 用于提高分辨率：映射改为 $x=\mu t/\sqrt{1-t^2}$，基改为

$$
R^{\lambda}_{n,\mu}(x)=\frac{\mu^{\lambda+\frac12}}{(\mu^2+x^2)^{\frac{\lambda+1}{2}}}
C^{\lambda}_n\!\Bigl(\frac{x}{\sqrt{\mu^2+x^2}}\Bigr)
=\mu^{-\frac12}R^{\lambda}_n\!\Bigl(\frac{x}{\mu}\Bigr).
$$

它是编号 27 中 Hermite 尺度因子的有理对应物，重要差别在于：这里 $\mu$ 控制代数型剖面的**宽度**，而 $\lambda$ 独立控制尾部的**指数**。两个参数分离，是有理基相对 Hermite 基的实际优势。

### 逼近空间与两种格式

$\mathcal V^{\lambda}_N=\mathrm{span}\{R^{\lambda}_n\}=\{S(t)P(t):P\in\mathcal P_N\}$。换元 $\breve u(x)=u(x)/s(x)=\breve U(t)$（$s(x)=(1+x^2)^{-(\lambda+1)/2}=S(t)$）给出 $\hat u^{\lambda}_n=\widehat{\breve U}_n$，即 $\mathbb R$ 上的修正映射 Gegenbauer 逼近**就是** $(-1,1)$ 上的 Gegenbauer 逼近。映射后的 Gauss 节点与权为

$$
x^{\lambda}_j=\frac{t^{\lambda}_j}{\sqrt{1-(t^{\lambda}_j)^2}},
\qquad
\omega^{\lambda}_j=\bigl(1+(t^{\lambda}_j)^2\bigr)^{-\lambda}\rho^{\lambda}_j .
$$

Galerkin 格式取 $\tilde a_s(u,v)=((-\Delta)^{s/2}u,(-\Delta)^{s/2}v)+\rho(u,v)$，在 $\mathcal V^{\lambda}_N$ 内求解，源项被插值，因此误差界另带一个 $f$ 项。配点格式用映射 Gauss 网格上的 Lagrange 基给出分数阶微分矩阵

$$
\mathcal D_{i,j}=(-\Delta)^{\alpha/2}l_j(x^{\lambda}_i)
=\sum_{k=0}^{N-1}b^j_k\,(-\Delta)^{\alpha/2}R^{\lambda}_k(x^{\lambda}_i),
\qquad
b^j_k=\frac{R^{\lambda}_k(x^{\lambda}_j)\,\omega^{\lambda}_j}{\gamma^{\lambda}_k} .
$$

论文明确指出配点格式的收敛分析「看起来不平凡且基本仍未解决」，这一自陈的限定值得保留。

多维情形转到 Fourier 域，那里问题是对角的：

$$
(|\xi|^{\alpha}+\rho)\hat u(\xi)=\hat f(\xi),
\qquad
\hat u(\xi)=a(\xi)\hat f(\xi),\quad a(\xi)=\frac{1}{|\xi|^{\alpha}+\rho} .
$$

算法是：在张量基网格上插值 $f$，变换，乘 $a(\xi)$，再投影回来。这里有一个自对偶结构：$R^{\lambda}_n$ 的逆变换由**同一组**展开公式给出。选择有理基在这里有额外的理由——函数的 Fourier 变换通常比函数本身衰减更慢，因此在变换后的变量里，代数衰减的基正是合适的工具。

## 92：一致性而非速度

编号 92 处理的是**一致性**缺口。摘要指出，已有方法在 $\alpha$ 接近 $0$ 或 $1$ 时表现不佳。结构原因是清楚的：标准积分表示带前因子 $\sin(\alpha\pi)$，它在整数 $\alpha$ 处为零，而被积函数的端点行为强度由 $\alpha$ 支配。例如 Balakrishnan 型与其负幂形式分别为

$$
A^{\alpha}=\frac{\sin(\alpha\pi)}{\alpha\pi}\,A\int_0^{\infty}\bigl(t^{1/\alpha}I+A\bigr)^{-1}\mathrm dt,
\qquad
\mathcal L^{-\alpha}=\frac{2\sin(\alpha\pi)}{\pi}\int_0^{\infty}t^{2\alpha-1}
\bigl(\mathcal I+t^{2}\mathcal L\bigr)^{-1}\mathrm dt .
$$

在最佳一致有理逼近一族中，误差渐近也带同一个因子：$t^{\beta-\alpha}$ 在 $[0,1]$ 上的最佳 $(k,k)$ 阶有理逼近满足

$$
\lim_{k\to\infty}e^{2\pi\sqrt{(\beta-\alpha)k}}E_{\alpha}(k,k;\beta)
=4^{\,1+\beta-\alpha}\,\bigl|\sin\pi(\beta-\alpha)\bigr| ,
$$

即根指数收敛，而 $|\sin\pi(\beta-\alpha)|$ 正是 $\alpha\to0,1$ 处的退化来源。

论文的目标是让代价与精度在 $\alpha\in(0,1)$ 上一致，处理两个对象：$A^{-\alpha}$（$A$ 正定），以及 $(q\mathcal I+A^{\alpha})^{-1}$（$q>0$）。后者正是空间分数阶演化方程隐式时间步所需的解算子，这解释了论文列出的应用。

> [!warning] 可核实范围
> 该文正文未能通过公开渠道获取，且无预印本。可确认的信息是：期刊关键词为梯形规则、Laguerre-Gauss 求积、分数阶 Laplacian、分数阶 Poisson 方程、空间分数阶 Allen-Cahn 方程；参考文献中包含 AAA 有理逼近算法、分数阶扩散的双指数求积、矩阵分数次幂的双指数公式、指数收敛的梯形规则等工作。据此可判断方法属于**基于求积的方案**，并与双指数求积与有理逼近方案作比较。本站不报告其定理、收敛率或「一致快速」的确切含义，这些内容需回查期刊版。

## 三篇的对照

| 编号 | 基 / 表示                | 远场衰减 | 算子的处理方式             |
| ---- | ------------------------ | -------- | -------------------------- |
| 27   | Hermite 函数与过缩放基   | 高斯型   | 物理空间的显式微分矩阵     |
| 42   | 修正映射 Gegenbauer 函数 | $        | x                          | ^{-(\lambda+1)}$ | 分解成两种初等有理形状 |
| 92   | 不涉及基（矩阵函数）     | 不适用   | 积分表示加求积（限定核实） |

## 覆盖核对

| 内容                           | 论文   | 覆盖状态                                     |
| ------------------------------ | ------ | -------------------------------------------- |
| 模型问题与算子的两个定义       | 27, 42 | 超奇异积分、常数、拟微分形式                 |
| 三条标准路线为何失效           | 42     | 截断、透明边界、指数衰减基                   |
| 两套 Hermite 基                | 27     | 归一化与过缩放形式及其正交关系               |
| 映射、修正基与正交性           | 42     | 映射恒等式、权吸收、三项递推、远场极限       |
| 两种初等形状的分数阶 Laplacian | 42     | 闭式、常数 $A^{\gamma}_s$、Bessel 证明路径   |
| 整数阶与分数阶的衰减分岔       | 42     | 三分支结论与其含义                           |
| 递推求值与尺度参数             | 42     | 连带关系、$\mu$ 与 $\lambda$ 的分工          |
| 逼近空间、节点权、两种格式     | 42     | 等价性、映射 Gauss 数据、微分矩阵、开放问题  |
| 多维 Fourier 域算法与自对偶    | 42     | 对角化、算法步骤、选基理由                   |
| 一致性缺口与 $\sin(\alpha\pi)$ | 92     | 两个积分表示、有理逼近误差渐近、两个目标对象 |

## 本页原文

- T. Tang, H. Yuan, and T. Zhou, [_Hermite spectral collocation methods for fractional PDEs in unbounded domains_](https://doi.org/10.4208/cicp.2018.hh80.12), Commun. Comput. Phys. 24(4) (2018), pp. 1143-1168（预印本 [arXiv:1801.09073](https://arxiv.org/abs/1801.09073)）。
- T. Tang, L.-L. Wang, H. Yuan, and T. Zhou, [_Rational spectral methods for PDEs involving fractional Laplacian in unbounded domains_](https://doi.org/10.1137/19M1244299), SIAM J. Sci. Comput. 42(2) (2020), pp. A585-A611（预印本 [arXiv:1905.02476](https://arxiv.org/abs/1905.02476)）。
- Y. Duan, F. Zeng, H. Zhang, and T. Zhou, [_Fast computation of the fractional power of a matrix_](https://doi.org/10.1137/25M1757411), SIAM J. Sci. Comput. 48(1) (2026), pp. A309-A334。
