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

## 在无界区域上，要匹配的是衰减率

有界区域上选基，选的是光滑性；无界区域上选基，选的是**远场衰减率**。基函数在 $|x|\to\infty$ 处以什么方式归零，决定了有限展开能以多快的速度吃掉解的尾部；两者的衰减型一旦不同，精度先在尾部失效，再把整体收敛率拖下来。

两族候选基的尾部属于不同类型。Hermite 函数带一个高斯因子 $e^{-x^2/2}$，因此与**指数型衰减**的解相配；映射有理（Gegenbauer）基按 $|x|^{-(\lambda+1)}$ 代数衰减，而且指数由参数 $\lambda$ 调，因此与**代数衰减**的解相配。这一句话就是编号 27 与编号 42 的全部分歧所在。

强迫人作出这个选择的是分数阶算子自己。超奇异核 $|x-y|^{-(d+\alpha)}$ 衰减很慢，算子是非局部的：远处的取值会经核持续反馈到近处，因此**分数阶问题的解一般不按高斯衰减，哪怕数据按高斯衰减**。编号 42 的数值实验把这件事做成了实证：源项分别取指数衰减的 $e^{-x^2/2}(1+x)$ 与代数衰减的 $(1+x^2)^{-2}$，两者的解都按同一条幂律

$$
\lvert u(x)\rvert\sim\lvert x\rvert^{-\alpha-1}
$$

衰减。**尾指数是算子定的，不是数据定的**——普通 Laplacian 没有这个现象，这也是三条标准路线在下一节同时失效的根源。

两篇论文的成绩与限制因此可以预先读出来：编号 27 的 Hermite 基对指数衰减的解给出谱收敛，而在它自己的特征值算例上——那里的特征函数只代数衰减——退化为代数收敛；编号 42 把基换成尾指数可调的有理函数，正是为了让基的尾巴对上 $|x|^{-\alpha-1}$。

## 无界区域的三条标准路线都不适用

编号 42 把困难说得很清楚：分数阶算子的解在无穷远**只是代数衰减**，服从幂律。这一条使三条标准策略同时失效。

区域截断适合快速衰减的解，这里不适合；而且朴素截断会在无界区域被终止的界面处引入非物理奇性。透明边界条件与吸收层因算子非局部而变得极难构造。无界区域上的正交函数——Hermite 与 Laguerre——是为**指数衰减**调好的，与幂律尾部不匹配。

映射 Jacobi（有理）基对**整数阶**且解代数衰减的问题早已知道优于 Hermite 与 Laguerre。但把映射技巧推广到分数阶远非平凡，原因是结构性的：普通 Laplacian 把有理基函数映到同类函数，而分数阶 Laplacian 把它映到**性质完全不同的一类函数**。

## 27：Hermite 基与显式微分矩阵

### 直觉：把超奇异积分一次算完，换来一张离线矩阵

有界区域上成功的做法是选一个与分数阶算子相容的基——Jacobi 多重分数阶函数（一个分数阶 Sturm-Liouville 问题的特征函数）与广义 Jacobi 函数——它们的分数阶导数仍是同族中参数移位后的函数，因此分数阶算子实际上变成局部的。无界区域上没有这样的装置。

没有这个装置，代价立刻现形：分数阶 Laplacian 非局部，离散化本来就给出稠密矩阵，而每个矩阵元还是一个超奇异积分，要在求解过程中反复求值。同时期的替代方案（Mao 与 Shen 的 Hermite 谱方法）转到频率空间做配点，因此每次求值都要一次正向与一次反向 Hermite 变换，对非线性问题很不方便。编号 27 的贡献是在**物理空间**直接配点，依据是显式闭式微分矩阵。

这条路之所以走得通，靠的是一个很朴素的性质：Hermite 函数的 Fourier 像仍是 $\xi$ 的单项式乘一个高斯因子。于是「乘 $|\xi|^{\alpha}$ 再变回去」这一步落在一族已知的积分里，结果是合流超几何函数 ${}_1F_1$ 的闭式。**超奇异积分因此只算一次，而且是解析地算完的**：矩阵与数据无关，可以完全离线预计算，非线性项在配点上逐点求值并用 Newton 迭代处理。与频率空间方案的实质差别就在这里——那边每次残差求值付两次 Hermite 变换，这边只付一次矩阵-向量乘。

### 问题设定

本文的模型问题带一个（线性或非线性的）反应项：

$$
(-\Delta)^{\alpha/2}u(x)+\gamma f(u)=g(x),\quad x\in\mathbb R^d,
\qquad u(x)=0,\ |x|\to\infty .
$$

分数阶 Laplacian 取开头给出的超奇异形式，也等价地取拟微分形式 $\mathcal F[(-\Delta)^{\alpha/2}u](\xi)=|\xi|^{\alpha}\mathcal F[u](\xi)$；$0<\alpha<2$，$\alpha\to2$ 时回到普通 Laplacian。**全篇的推导用的是拟微分定义**，超奇异形式只用来说明算子的非局部性。

Hermite 多项式由三项递推 $H_0=1$、$H_1(x)=2x$、$H_{n+1}(x)=2xH_n(x)-2nH_{n-1}(x)$ 给出，满足 $\int_{\mathbb R}H_mH_n e^{-x^2}\mathrm dx=\gamma_n\delta_{mn}$，$\gamma_n=\sqrt{\pi}\,2^n n!$。由它构造两套基。归一化 Hermite 函数

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

由 Brinkman 为 Fokker-Planck 方程引入，其中速度方向的概率密度按这组函数展开，此后成为该方程的标准基之一。两套基的差别在于把 $e^{-x^2/2}$ 放在基里还是放在权里，而这决定了离散内积与微分矩阵的形式——后面会看到，它还决定了条件数的增长方式，这是全篇第二个主题。

### 推导：过缩放基的分数阶像

配点格式先把解按过缩放基展开，$u(x)\approx u_N(x)=\sum_{n=0}^{N-1}c_n\widetilde H_n(x)$，再在第 $N$ 个 Hermite 多项式的根 $\{x_i\}_{i=0}^{N-1}$ 上要求方程成立：

$$
\sum_{n=0}^{N-1}c_n\,(-\Delta)^{\alpha/2}\widetilde H_n(x_i)+\gamma f\bigl(u_N(x_i)\bigr)=g(x_i),
$$

即 $\widetilde{\mathcal D}^{\alpha}c+\gamma F(c)=g$，其中 $\widetilde{\mathcal D}^{\alpha}_{i,j}=(-\Delta)^{\alpha/2}\widetilde H_j(x_i)$。**全部工作因此归结为一件事：把 $(-\Delta)^{\alpha/2}\widetilde H_n$ 算出来。**

推导链只有三步。第一步写下过缩放基的 Fourier 像，它是单项式乘高斯：

$$
\mathcal F\bigl[\widetilde H_{2n}\bigr](\xi)
=\frac{(-1)^n}{\sqrt2\,\sqrt{2^{2n}(2n)!}}\,\xi^{2n}e^{-\xi^2/4}.
$$

第二步按拟微分定义乘上 $|\xi|^{\alpha}$ 再作逆变换；由偶性，逆变换塌缩成一个余弦积分

$$
\int_{\mathbb R^+}\xi^{2n+\alpha}e^{-\xi^2/4}\cos(x\xi)\,\mathrm d\xi .
$$

第三步用这个积分的已知求值收尾，得到 ${}_1F_1$。结果是定理 3.1 与 3.2：对 $0<\alpha<2$，

$$
(-\Delta)^{\alpha/2}\widetilde H_{2n}(x)=\frac{(-1)^n\,2^{\,n+\alpha}}{\sqrt{\pi}\,\sqrt{(2n)!}}\;
\Gamma\!\Bigl(n+\tfrac{\alpha+1}{2}\Bigr)\;
{}_1F_1\!\Bigl(n+\tfrac{\alpha+1}{2},\ \tfrac12;\ -x^2\Bigr),
$$

$$
(-\Delta)^{\alpha/2}\widetilde H_{2n+1}(x)=\frac{(-1)^n\,2^{\,n+\alpha+\frac32}}{\sqrt{\pi}\,\sqrt{(2n+1)!}}\;
\Gamma\!\Bigl(n+\tfrac{\alpha+3}{2}\Bigr)\;
x\,{}_1F_1\!\Bigl(n+\tfrac{\alpha+3}{2},\ \tfrac32;\ -x^2\Bigr).
$$

值得注意的是奇偶两式的结构差别只在 ${}_1F_1$ 的第二个参数（$\tfrac12$ 对 $\tfrac32$）与前面多出的一个 $x$，这正是分数阶算子保奇偶性的痕迹。用 $\Gamma\bigl(n+\tfrac12\bigr)=\sqrt{\pi}\,(2n)!/(4^n n!)$ 可以把上面的前因子换成 $\Gamma$ 比值的写法，两者等价。取 $n=0$、$\alpha\to2$ 即回到 $-\bigl(e^{-x^2}\bigr)''=(2-4x^2)e^{-x^2}$，可作校验。

实现上还有一步节省。逐个矩阵元调用 ${}_1F_1$ 太贵，论文只对 $j=0,1,2,3$ 用闭式（这几个需要一个快而准的 ${}_1F_1$ 求值程序），其余 $4<j\le N-1$ 用 ${}_1F_1$ 的连带关系

$$
(2a-b+x)\,{}_1F_1(a,b;x)=a\,{}_1F_1(a+1,b;x)-(b-a)\,{}_1F_1(a-1,b;x)
$$

递推填充。矩阵与数据无关，因此整张表离线算一次即可。

### 推导：归一化基走 Fourier 特征函数这条路

归一化基的推导换一个入口，因为它是 Fourier 变换的特征函数：

$$
\mathcal F\bigl[\widehat H_n\bigr](\xi)
=\frac{1}{\sqrt{2\pi}}\int_{\mathbb R}\widehat H_n(x)e^{-i\xi x}\mathrm dx
=(-i)^n\widehat H_n(\xi).
$$

于是把 $\widehat H_n$ 自己按 $\widehat H_n(\xi)=\sum_{k=0}^{n}\hat a_{n,k}e^{-\xi^2/2}\xi^k$ 展开，$\hat a_{n,k}=a_{n,k}/\sqrt{2^n n!}$，多项式系数由

$$
a_{0,0}=1,\quad a_{1,0}=0,\quad a_{1,1}=2,
\qquad
a_{n+1,0}=-a_{n,1},
\qquad
a_{n+1,k}=2a_{n,k-1}-(k+1)a_{n,k+1}\ (k>0)
$$

生成。这样分数阶算子只落在 $e^{-\xi^2/2}\xi^k$ 这一族上，令 $F_k(x)=\mathcal F^{-1}\bigl[e^{-\xi^2/2}\xi^k|\xi|^{\alpha}\bigr](x)$，两个奇偶情形都有闭式：

$$
F_{2m}(x)=\frac{2^{\,m+\frac{\alpha}{2}}}{\sqrt{\pi}}\,
\Gamma\!\Bigl(\frac{2m+1+\alpha}{2}\Bigr)\,
{}_1F_1\!\Bigl(\frac{2m+1+\alpha}{2},\ \frac12;\ -\frac{x^2}{2}\Bigr),
$$

$$
F_{2m+1}(x)=\frac{2^{\,m+1+\frac{\alpha}{2}}\,i}{\sqrt{\pi}}\,
\Gamma\!\Bigl(\frac{2m+3+\alpha}{2}\Bigr)\,
x\,{}_1F_1\!\Bigl(\frac{2m+3+\alpha}{2},\ \frac32;\ -\frac{x^2}{2}\Bigr),
$$

微分矩阵随之写成一个和：

$$
\widehat{\mathcal D}^{\alpha}_{mn}=(-\Delta)^{\alpha/2}\widehat H_n(x_m)
=(-i)^n\sum_{k=0}^{N-1}\hat a_{n,k}F_k(x_m),
\qquad 0\le n,m\le N-1,
$$

约定 $k>n$ 时 $a_{n,k}=0$。注意 $F_k$ 的自变量是 $-x^2/2$ 而过缩放基那边是 $-x^2$，两套公式不能互相搬用。

同一个逼近空间还可以配 Lagrange 型的节点基 $h_j(x_k)=\delta_{jk}$，这给出**另一张**微分矩阵，但它由同一批 $(-\Delta)^{\alpha/2}\widehat H_k(x_i)$ 值组装——两张矩阵的条件数不同，这是下一小节要用到的事实。

### 推导：尺度因子、多项模型与二维

**尺度因子。** Hermite 型谱收敛在解的衰减率与基的衰减率不匹配时会变差。若 $|u(x)|\approx0$ 对一切 $|x|>M$ 成立，就改为展开

$$
u(x)=\sum_{n=0}^{N-1}c_n\widetilde H_n(rx)
\iff
u(x/r)=\sum_{n=0}^{N-1}c_n\widetilde H_n(x),\qquad r>0,
$$

使缩放后的配点 $\{x_k/r\}$ 落在 $u$ 的有效支集内。论文建议的取法是

$$
\max_{0\le k\le N-1}|x_k|/r\le M
\ \Longrightarrow\
r=\max_{0\le k\le N-1}|x_k|/M ,
$$

并明确指出实践中 $M$ 本身不易确定，所以最优 $r$ 一般难以拿到。带尺度时，微分矩阵的元就是同样的公式在 $z_i=rx_i$ 处求值，没有新的推导。

**多项与分布阶模型。** 对 $\sum_{j=1}^{J}(-\Delta)^{\alpha_j/2}u+\gamma f(u)=g$（其动机是分布阶模型的求积近似），配点系统只是

$$
\mathcal Jc+\gamma F(c)=g,
\qquad
\mathcal J=\sum_{j=1}^{J}\widetilde{\mathcal D}^{\alpha_j},
$$

即微分矩阵直接相加。**这是显式闭式路线的一个实惠：多加一个分数阶项不带来任何结构代价。**

**二维。** 取张量积基 $\{\widetilde H_n(x)\widetilde H_m(y)\}$，由 $\mathcal F$ 的可分性，

$$
(-\Delta)^{\alpha/2}\bigl[\widetilde H_n(x)\widetilde H_m(y)\bigr]
=\frac{1}{2\pi}\int_{\mathbb R}\int_{\mathbb R}(\xi^2+\eta^2)^{\alpha/2}
\mathcal F[\widetilde H_n](\xi)\,\mathcal F[\widetilde H_m](\eta)\,
e^{ix\xi}e^{iy\eta}\,\mathrm d\xi\,\mathrm d\eta .
$$

角向积分用 Bessel 恒等式化掉——关键是 $J_{-1/2}(x)=\sqrt{2/(\pi x)}\cos x$ 与一个 Sonine-Gegenbauer 型公式——最后每个矩阵元仍是单个 ${}_1F_1$ 值。双偶块为

$$
\widetilde{\mathcal D}^{\alpha}_{(i-1)N+j,\,(p-1)N+q}
=\frac{2^{\alpha}\,\Gamma\bigl(n+m+\tfrac{\alpha}{2}+1\bigr)}
{\sqrt{2^{2n}(2n)!}\,\sqrt{2^{2m}(2m)!}\;\Gamma(2n+2m+1)}\;
{}_1F_1\!\Bigl(n+m+\tfrac{\alpha}{2}+1;\ 2n+2m+1;\ -(x_i^2+y_j^2)\Bigr)
$$

（$p=2n$、$q=2m$，另有三个类似情形）。这里有一个漂亮的结构：**二维矩阵元只通过 $x_i^2+y_j^2$ 依赖节点**，因此二维并没有把公式复杂度按维数放大。

### 定理：有闭式公式，但没有收敛定理

这是一篇算法论文，**它不含收敛定理**。严格的内容就是上面那几组闭式求值（过缩放基的定理 3.1-3.2 与归一化基的 $F_k$ 公式）。谱收敛是经验结论：引言的措辞是「两种方法对在无穷远指数衰减的解都具有谱收敛」，这是一句断言而不是被证明的定理。

论文真正下了功夫的定量结论在条件数上。对过缩放基 $\{\widetilde H_n\}$，$\widetilde{\mathcal D}^{\alpha}$ 的条件数「随 $N$ 增长很快」（图 1 的纵轴跨 $10^0$ 到 $10^{10}$，$\alpha=0.4,1,1.6$），论文把原因归给这组基自身的病态，这也是引入归一化基的动机。对归一化基，增长是**代数**的，原文给出拟合斜率：

| 微分矩阵                                | $\alpha=0.4$ | $\alpha=1$ | $\alpha=1.6$ |
| --------------------------------------- | ------------ | ---------- | ------------ |
| $\widehat{\mathcal D}^{\alpha}$（图 2） | $0.446$      | $1.065$    | $1.684$      |
| Lagrange 节点基矩阵（图 3）             | $0.428$      | $1.054$    | $1.659$      |

三个斜率与三个 $\alpha$ 几乎相等，看起来条件数按 $N^{\alpha}$ 增长，不过论文没有把指数认定为 $\alpha$。论文明确说，实际应用中应当为此设计有效的预条件子。

### 数值实验：六个算例

论文 §6 给出六组算例，误差量级与收敛斜率以图呈现。

| 算例                      | 方程与精确解                                                                                                | 尺度因子                                        | 结论                                                                                   |
| ------------------------- | ----------------------------------------------------------------------------------------------------------- | ----------------------------------------------- | -------------------------------------------------------------------------------------- |
| §6.1 一维基准             | 精确解 $u(x)=e^{-x^2}\sin x$                                                                                | 过缩放基 $r=1$（无需缩放）；归一化基 $r=\sqrt2$ | $\alpha=0.4,1,1.6$ 在加权范数与最大范数下均谱收敛；过缩放基在大 $N$ 处被条件数增长污染 |
| §6.2 尺度因子（过缩放基） | 精确解 $u(x)=e^{-x^2/2}x^2\cos x$                                                                           | $r=1/\sqrt2$ 与 $r=1$ 对比                      | 取对的尺度因子给出明显更快的速率                                                       |
| §6.2 尺度因子（归一化基） | 精确解 $u(x)=e^{-2x^2}x^2\cos x$                                                                            | 最优 $r=2$ 与 $r=1$ 对比                        | 同上                                                                                   |
| §6.3 二维                 | $(-\Delta)^{\alpha/2}u+2u=g$，精确解 $u=e^{-(x^2+y^2)}\sin(x+y)$                                            | 张量积基                                        | —                                                                                      |
| §6.4 多项模型             | $\sum_{j=1}^{4}(-\Delta)^{\alpha_j/2}u=g$，精确解 $u=e^{-3x^2/2}(\sin x+x^6+x^2\cos x)$                     | $r=\sqrt{1.5},\sqrt{1.3},1$ 对比                | $r=\sqrt{1.5}$ 明显最好                                                                |
| §6 非线性                 | $(-\Delta)^{\alpha/2}u+u^2=g$，精确解 $e^{-x^2}(\sin x+x^2)$（过缩放）与 $e^{-x^2/2}(\sin x+x^2)$（归一化） | 未报告                                          | Newton 迭代容差 $10^{-16}$，观察到谱收敛                                               |
| §6 特征值                 | $\bigl[(-\Delta)^{\alpha/2}+x^2\bigr]u=\lambda u$                                                           | 未报告                                          | 只观察到误差的**代数**衰减                                                             |

多项模型的四个阶取自变换后的 Legendre-Gauss 点：

| $j$        | 1       | 2       | 3       | 4       |
| ---------- | ------- | ------- | ------- | ------- |
| $\alpha_j$ | $0.139$ | $0.660$ | $1.340$ | $1.861$ |

特征值算例是全篇最有诊断价值的一个，因为它有独立参照：$\alpha=1$ 时谱已知，$\lambda_{2k-1}=-a'_k$、$\lambda_{2k}=-a_k$，其中 $a_k,a'_k$ 是 Airy 函数

$$
A(x)=\frac1\pi\int_0^{\infty}\cos\Bigl(\frac{t^3}{3}+xt\Bigr)\mathrm dt
$$

及其导数按降序排列的零点。所用参考值为

| 特征值      | 参考值（$\alpha=1$） |
| ----------- | -------------------- |
| $\lambda_1$ | $1.01879297164747$   |
| $\lambda_2$ | $2.33810741045976$   |
| $\lambda_3$ | $3.24819758217983$   |

而这里只观察到代数收敛，论文把原因归给特征函数只**代数**衰减。

这组实验建立的是：物理空间配点加离线预计算矩阵确实对指数衰减的解给出谱收敛，且在三个 $\alpha$ 上一致；尺度因子不是修饰而是决定速率的量；多项模型确实只需把矩阵相加；非线性问题可以直接上 Newton，而这正是频率空间方案处理起来别扭的情形。

它在三处不及理论。第一，谱收敛没有定理支撑，条件数的结论也只是拟合斜率，而过缩放基在大 $N$ 处「收敛被污染」意味着**实际可达精度由条件数而非逼近性质决定**，论文只能把预条件列为待办。第二，最优尺度因子依赖有效支集 $M$，而论文自陈 $M$ 难求，所以 §6.2 里那些「取对了的 $r$」在真实问题上未必拿得到。第三，特征值算例把方法的适用边界画了出来：解一旦只代数衰减，谱收敛立刻退化成代数收敛。**这一条正是编号 42 要解决的问题，编号 27 自己已经把它演示了一遍。**

## 42：把有理基分解成两种初等形状

### 直觉：把基的尾巴对到算子产生的尾巴上

编号 42 的全部设计都从开头那条经验事实出发：解的尾巴是 $|x|^{-\alpha-1}$。要让谱逼近吃掉这样的尾巴，基函数自己就得按幂律衰减，而且幂次要可调——固定成高斯型或指数型都不行，因为 $\alpha$ 是问题的参数，尾指数会随它变。修正映射 Gegenbauer 函数正是这样一族：它们按 $|x|^{-(\lambda+1)}$ 衰减，$\lambda$ 是尾指数的旋钮，另有一个尺度参数 $\mu$ 独立控制剖面的宽度。基与解的两个衰减率因此可以真正对齐，而不是靠加大 $N$ 去硬顶。

代价出现在另一头，而且是结构性的。普通 Laplacian 把有理基函数映成同类的有理函数，所以整套映射谱方法在整数阶下是自封闭的；分数阶 Laplacian 不封闭，它把有理基映到超几何函数，甚至连衰减率都不一定跟着基的参数走。于是本文要做的事很具体：**把基分解成两种初等有理形状，对这两种形状把 $(-\Delta)^{s}$ 一次算清楚，再线性组合回去。** 能不能拿到这两个闭式，是整篇论文成立与否的分界——拿到了，就有与编号 27 同构的离线微分矩阵；拿不到，映射技巧在分数阶下就无法落地。

### 问题设定：映射与基

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

### 推导：分解成两种初等形状

用 Pfaff 变换把 Gegenbauer 多项式改写成以 $(1+x^2)^{-1}$ 为自变量的 ${}_2F_1$，基函数就展成有限和：

$$
R^{\lambda}_{2n}(x)=a^{\lambda}_n\sum_{k=0}^{n}
\frac{(-n)_k(n+\lambda)_k}{(\lambda+\tfrac12)_k\,k!}\cdot
\frac{1}{(1+x^2)^{\,k+\frac{\lambda+1}{2}}},
$$

$$
R^{\lambda}_{2n+1}(x)=b^{\lambda}_n\sum_{k=0}^{n}
\frac{(-n)_k(n+\lambda+1)_k}{(\lambda+\tfrac12)_k\,k!}\cdot
\frac{x}{(1+x^2)^{\,k+\frac{\lambda}{2}+1}} .
$$

两个常数 $a^{\lambda}_n$、$b^{\lambda}_n$ 只依赖 $n$ 与 $\lambda$，结构性结论与它们的具体形式无关：无论两个常数是什么，展开只含两种形状，所以

**每个基函数都是两种初等形状的有限线性组合**：

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

再由 $(-\Delta)^s v=\mathcal F^{-1}[|\xi|^{2s}\mathcal F[v]]$ 化为一个 Bessel 积分

$$
(-\Delta)^{s}v=\frac{2^{1-\gamma}}{\sqrt{2\pi}\,\Gamma(\gamma)}
\int_0^{\infty}\cos(x\xi)\,\xi^{2s+\gamma-\frac12}K_{\gamma-\frac12}(\xi)\,\mathrm d\xi,
$$

再用相应的积分恒等式（取 $\lambda=2s+\gamma-\tfrac12$、$\mu=\gamma-\tfrac12$、$b=x$）完成。奇次形状经 $\mathcal F[xv](\xi)=-i\frac{2^{1-\gamma}}{\Gamma(\gamma)}\xi^{\gamma-\frac12}K_{\gamma-\frac32}(\xi)$ 与 $\Gamma(z+1)=z\Gamma(z)$ 同法处理。**修正 Bessel 函数是这条推导的必经之处**，论文自己也点名「与修正 Bessel 函数有关的积分恒等式」。

同一结果用 Pfaff/Euler 变换还有一种等价写法（推论 3.3），论文说前一种便于计算、后一种便于分析：

$$
(-\Delta)^{s}\Bigl\{\frac{1}{(1+x^2)^{\gamma}}\Bigr\}
=\frac{A^{\gamma}_s}{(1+x^2)^{s+\gamma}}\
{}_2F_1\!\Bigl(-s,\ s+\gamma;\ \tfrac12;\ \tfrac{x^2}{1+x^2}\Bigr),
$$

$$
(-\Delta)^{s}\Bigl\{\frac{x}{(1+x^2)^{\gamma}}\Bigr\}
=\frac{(2s+1)A^{\gamma}_s\,x}{(1+x^2)^{s+\gamma}}\
{}_2F_1\!\Bigl(-s,\ s+\gamma;\ \tfrac32;\ \tfrac{x^2}{1+x^2}\Bigr).
$$

第二种写法把衰减因子 $(1+x^2)^{-(s+\gamma)}$ 显式提到前面，${}_2F_1$ 的自变量则被压在 $[0,1)$ 内——**这正是下一小节判断衰减率的入口**：整数阶与分数阶的差别，全在这个 ${}_2F_1$ 是否终止。

### 推导：整数阶与分数阶的分岔

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

### 基函数的分数阶像与两个可实现性细节

把两个初等形状的闭式代回 Pfaff 展开，就得到论文所称的「主要构件」，即基函数自身的分数阶像（定理 3.4，对实 $s>0$ 与 $\lambda>-1/2$）：

$$
(-\Delta)^{s}R^{\lambda}_{2n}(x)=a^{\lambda}_n\sum_{k=0}^{n}
\frac{(-n)_k(n+\lambda)_k}{(\lambda+\tfrac12)_k\,k!}\;
A^{\,k+\frac{\lambda+1}{2}}_{s}\
{}_2F_1\!\Bigl(s+k+\tfrac{\lambda+1}{2},\ s+\tfrac12;\ \tfrac12;\ -x^2\Bigr),
$$

$$
(-\Delta)^{s}R^{\lambda}_{2n+1}(x)=(2s+1)\,b^{\lambda}_n\,x\sum_{k=0}^{n}
\frac{(-n)_k(n+\lambda+1)_k}{(\lambda+\tfrac12)_k\,k!}\;
A^{\,k+\frac{\lambda}{2}+1}_{s}\
{}_2F_1\!\Bigl(s+k+\tfrac{\lambda}{2}+1,\ s+\tfrac32;\ \tfrac32;\ -x^2\Bigr).
$$

分数阶 Laplacian 作用在基上的公式不如 Hermite 情形紧凑（论文自己也这么说），但可用递推高效求值；论文报告用 Maple 或 Mathematica 可以把有理基的分数阶 Laplacian 准确算到 $10^3$ 量级的次数。令 $F_k(x)={}_2F_1(a,b;c;-x^2)$，$a=s+k+\tfrac{\lambda+1}{2}$、$b=s+\tfrac12$、$c=\tfrac12$，则

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

论文明确指出配点格式的收敛分析「看起来不平凡且基本仍未解决」。

多维情形转到 Fourier 域，那里问题是对角的：

$$
(|\xi|^{\alpha}+\rho)\hat u(\xi)=\hat f(\xi),
\qquad
\hat u(\xi)=a(\xi)\hat f(\xi),\quad a(\xi)=\frac{1}{|\xi|^{\alpha}+\rho} .
$$

算法是：在张量基网格上插值 $f$，变换，乘 $a(\xi)$，再投影回来。这里有一个自对偶结构：$R^{\lambda}_n$ 的逆变换由**同一组**展开公式给出。选择有理基在这里有额外的理由——函数的 Fourier 变换通常比函数本身衰减更慢，因此在变换后的变量里，代数衰减的基正是合适的工具。

### 定理：投影、插值与 Galerkin 的最优估计

误差理论用两个空间。一个是通常的分数阶 Sobolev 空间

$$
H^r(\mathbb R)=\Bigl\{u\in L^2:\int_{\mathbb R}(1+|\xi|^2)^r\bigl|\mathcal F[u](\xi)\bigr|^2\mathrm d\xi<\infty\Bigr\};
$$

另一个是为有理基定制的空间，它用**映射导数**而不是普通导数来度量正则性：

$$
D_xu:=a(x)\frac{\mathrm d\breve u}{\mathrm dx}=\frac{\mathrm d\breve U}{\mathrm dt},
\qquad
D^2_xu=a(x)\frac{\mathrm d}{\mathrm dx}\Bigl\{a(x)\frac{\mathrm d\breve u}{\mathrm dx}\Bigr\}=\frac{\mathrm d^2\breve U}{\mathrm dt^2},\ \dots
$$

其中 $a(x)$ 是映射带来的因子，第二个等号说明 $D_x$ 就是在参考变量 $t$ 上求导。半范为

$$
\lvert u\rvert_{\mathbb B^m_{\lambda}(\mathbb R)}
=\bigl\|(1+x^2)^{-\frac{\lambda+m+1}{2}}D^m_xu\bigr\|_{L^2(\mathbb R)} .
$$

权的指数不是凑出来的，它由换元定死：$(-1,1)$ 上的 Gegenbauer 估计以 $\omega_{\lambda+m}(t)=(1-t^2)^{\lambda+m-\frac12}$ 为权，而 $1-t^2=(1+x^2)^{-1}$、$\mathrm dt=(1+x^2)^{-3/2}\mathrm dx$，两个因子相乘正好给出 $(1+x^2)^{-(\lambda+m+1)}$，开方即为上式。**这个选择是全部估计能取到最优阶的原因**：$\mathbb B^m_{\lambda}$ 中的正则性不是解在 $\mathbb R$ 上的普通光滑性，而是它拉回到 $(-1,1)$ 之后的光滑性，而幂律尾部拉回去之后正好是光滑的。

下面几条估计沿同一条路线得到。基本工具是空间插值引理（引理 4.1）：对 $r=(1-\theta)r_0+\theta r_1$、$\theta\in[0,1]$ 有 $\|u\|_{H^r}\le\|u\|^{1-\theta}_{H^{r_0}}\|u\|^{\theta}_{H^{r_1}}$，在 Fourier 侧用指数 $p=1/(1-\theta)$、$q=1/\theta$ 的 Hölder 不等式即得。

**定理 4.2（$L^2$ 投影的最优误差）。** 设 $u\in H^s(\mathbb R)\cap\mathbb B^m_{\lambda}(\mathbb R)$，整数 $1\le m\le N+1$，$s\in(0,1)$，$\lambda>-1/2$，则

$$
\|\pi^{\lambda}_Nu-u\|_{H^s(\mathbb R)}\le c\,N^{\,s-m}\,|u|_{\mathbb B^m_{\lambda}(\mathbb R)},
$$

$c$ 与 $N$、$u$ 无关。证明路线正是上面那条等价性：先把问题搬到 $(-1,1)$ 上的 Gegenbauer 逼近，对 $l=0,1$ 用已有估计 $\|(\Pi^{\lambda}_N\Phi-\Phi)^{(l)}\|_{L^2_{\omega_{\lambda+l}}}\le cN^{\,l-m}\|\Phi^{(m)}\|_{L^2_{\omega_{\lambda+m}}}$，再用插值引理在 $L^2$ 与 $H^1$ 之间插到 $H^s$。

**定理 4.3（$H^s$ 正交投影）。** 取 $a_s(u,v)=\bigl((-\Delta)^{s/2}u,(-\Delta)^{s/2}v\bigr)+(u,v)$ 及其诱导的投影 $\pi^{s}_{N,\lambda}$，则它是 $H^s$ 意义下的最佳逼近，且同阶：

$$
\|\pi^{s}_{N,\lambda}u-u\|_{H^s}=\inf_{\phi\in\mathcal V^{\lambda}_N}\|\phi-u\|_{H^s}
\le cN^{\,s-m}|u|_{\mathbb B^m_{\lambda}} .
$$

**定理 4.4（插值）。** 映射 Gegenbauer-Gauss 插值同阶：$\|I^{\lambda}_Nu-u\|_{H^s(\mathbb R)}\le cN^{\,s-m}|u|_{\mathbb B^m_{\lambda}(\mathbb R)}$。

**定理 5.1（Galerkin 格式的最优收敛）。** 设 $u\in H^s(\mathbb R)\cap\mathbb B^m_{\lambda}(\mathbb R)$，$f\in\mathbb B^{k}_{\lambda}(\mathbb R)$，整数 $1\le m,k\le N+1$，$s=\alpha/2\in(0,1)$，$\lambda>-1/2$，则

$$
\|u-u_N\|_{H^s(\mathbb R)}\le c\,N^{\,s-m}\,|u|_{\mathbb B^m_{\lambda}(\mathbb R)}
+c\,N^{-k}\,|f|_{\mathbb B^{k}_{\lambda}(\mathbb R)},
$$

$c$ 与 $N$、$u$、$f$ 均无关。两项的来历不同：第一项是解自身的逼近误差，第二项来自源项被插值这一步。**论文 §5.2 指出实践中通常是第二项占主导**，下一小节的数值结果就是这句话的直接后果。

多维情形的界形状相同（式 (7.21)）：若 $a\hat f\in\mathbb B^m_{\lambda}(\mathbb R^d)$、$f\in\mathbb B^{m'}_{\lambda}(\mathbb R^d)$，$m\ge0$、$m'\ge d$，则

$$
\|u_{NM}-u\|_{L^2(\mathbb R^d)}\le c\,M^{-m}\,|a\hat f|_{\mathbb B^m_{\lambda}(\mathbb R^d)}
+c\,N^{-m'}\,|f|_{\mathbb B^{m'}_{\lambda}(\mathbb R^d)} ,
$$

其中 $M$ 是频率侧的自由度、$N$ 是插值 $f$ 的自由度，两项分别由 $\|\pi^{\lambda}_Mu-u\|_{L^2}\le cM^{-m}|u|_{\mathbb B^m_{\lambda}}$ 与 $\|I^{\lambda}_Nu-u\|_{L^2}\le cN^{-m}|u|_{\mathbb B^m_{\lambda}}$（$m\ge d$）支撑。备注 7.2 记下了一处未追究的地方：若把 $\pi^{\lambda}_M$ 换成 $I^{\lambda}_M$，就需要插值算子在含 $d$ 阶偏导的范数下稳定，论文没有做这件事。

上述估计全部只覆盖 Galerkin 格式，**配点格式没有收敛理论**。理由也很清楚：分数阶 Laplacian 把有理基送出了逼近空间的自然类，前面的分岔结论就是它的定量形式。

### 数值实验：两类源项、两个 $\lambda$ 与一个二维算例

共同设置：误差在离散 $H^{\alpha/2}$ 范数下报告，$\rho=1$；只测两个 $\lambda$——$\lambda=0$（修正映射 Chebyshev 有理函数）与 $\lambda=0.5$（修正映射 Legendre 函数）；$\alpha$ 一律取 $0.4,1,1.6$。

| 算例                   | 数据                                                                  | 有理基参数                     | 对照方法                                 |
| ---------------------- | --------------------------------------------------------------------- | ------------------------------ | ---------------------------------------- |
| 算例 1（指数衰减真解） | 精确解 $u_e(x)=e^{-x^2}$，源项闭式                                    | $\mu=5$                        | 无（构造解）                             |
| 算例 1（代数衰减真解） | 精确解 $u_a(x)=(1+x^2)^{-r}$，$r=2.3$，源项闭式                       | $\mu=3$                        | 无（构造解）                             |
| 算例 2（指数衰减源项） | $f(x)=e^{-x^2/2}(1+x)$，无闭式解                                      | $\mu=5$，$\lambda\in\{0,0.5\}$ | Hermite 函数法（Mao-Shen），尺度 $1/0.4$ |
| 算例 3（代数衰减源项） | $f(x)=(1+x^2)^{-2}$                                                   | $\mu=3$，两个 $\lambda$ 都测   | Hermite 函数法，尺度 $1/0.7$             |
| 配点：多项模型         | 两个源项 $e^{-x^2/2}(1+x)$ 与 $(1+x^2)^{-1.8}$                        | 未报告                         | Hermite 函数法                           |
| 二维（§7）             | $f(x,y)=e^{-\sqrt{x^2+y^2}}$，$\mathcal F[f]=(1+\xi^2+\eta^2)^{-3/2}$ | 未报告                         | Hermite 配点                             |

算例 1 的两个源项都有闭式，这是它能当基准的原因：

$$
f_e(x)=\rho e^{-x^2}+\frac{2^{\alpha}\Gamma\bigl(\tfrac{\alpha+1}{2}\bigr)}{\Gamma\bigl(\tfrac12\bigr)}\
{}_1F_1\!\Bigl(\tfrac{\alpha+1}{2};\tfrac12;-x^2\Bigr),
$$

$$
f_a(x)=\rho(1+x^2)^{-r}
+\frac{2^{\alpha}\Gamma\bigl(\tfrac{\alpha}{2}+r\bigr)\Gamma\bigl(\tfrac{\alpha+1}{2}\bigr)}
{\Gamma(r)\,\Gamma\bigl(\tfrac12\bigr)}\
{}_2F_1\!\Bigl(\tfrac{\alpha}{2}+r,\ \tfrac{\alpha+1}{2};\ \tfrac12;\ -x^2\Bigr).
$$

结果本身是这一节最有说明性的地方：**即便真解是 $e^{-x^2}$ 这种指数衰减的函数，误差也只按代数速率下降**。原因由定理 5.1 的第二项给出——源项 $f_e$ 只按 $(1+x^2)^{-\frac{\alpha+1}{2}}$ 衰减，因此误差被**源项的插值误差**支配。定理 5.1 加一次直接计算给出预测速率 $O\bigl(N^{-(\alpha+\frac12)}\bigr)$，论文称与数值结果吻合良好。把三个被测 $\alpha$ 代进去：

| $\alpha$ | 预测速率   |
| -------- | ---------- |
| $0.4$    | $N^{-0.9}$ |
| $1$      | $N^{-1.5}$ |
| $1.6$    | $N^{-2.1}$ |

这三行是把 $\alpha$ 代入 $N^{-(\alpha+1/2)}$ 得到的预测值，不是实测值。取 $u_a$（$r=2.3$）时速率相同，因为 $f_a\sim(1+x^2)^{-\min(r,\frac{\alpha+1}{2})}$，最小值仍落在 $\frac{\alpha+1}{2}$ 一侧。

算例 2 与 3 是与 Hermite 方法的正面对照。算例 2 没有闭式解，参考解取 $N=600$ 的数值解。两个算例的结论一致：**有理基在所有情形下都优于 Hermite 逼近，收敛率明显更高**，速率数值列在原文的表 1（$\alpha=1$）与表 2 中。

原文图 4 报告的是计算所得解的渐近行为，它是解释性的那一组，也是开头那条事实的出处：两个衰减方式完全不同的源项给出的解都按 $|u(x)|\sim|x|^{-\alpha-1}$ 衰减。论文用这一句解释「为什么修正映射 Gegenbauer 函数比 Hermite 函数表现更好」。

配点格式用一个多项（分布阶型）模型测试，$\sum_{j=1}^{4}\rho_j(-\Delta)^{\alpha_j/2}u=f$，参数为

| $j$        | 1       | 2       | 3       | 4       |
| ---------- | ------- | ------- | ------- | ------- |
| $\alpha_j$ | $0$     | $0.5$   | $1.5$   | $2$     |
| $\rho_j$   | $\pi/6$ | $\pi/3$ | $\pi/3$ | $\pi/6$ |

四个权的比例为 $1:2:2:1$，总和为 $\pi$。结论与 Galerkin 相同：有理基在所有情形下都明显优于 Hermite 函数法。论文同时指出配点格式对变系数与非线性问题更实用。

二维算例 $f(x,y)=e^{-\sqrt{x^2+y^2}}$ 挑得很有心思：它在 $x$ 上指数衰减，而它的 Fourier 变换 $(1+\xi^2+\eta^2)^{-3/2}$ 只代数衰减，正是「在变换后的变量里该用有理基」这条论断的干净示例。结论是有理基配点比 Hermite 配点更准且收敛更快。

这组实验建立的是三件事：有理基在两类源项、Galerkin 与配点两种格式、以及一维与二维上都胜过 Hermite 方法；预测速率 $N^{-(\alpha+1/2)}$ 被实现，也就是说定理 5.1 中起作用的确实是源项那一项；幂律尾部 $|x|^{-\alpha-1}$ 由数值证实，尾指数由算子而非数据决定。

它在四处不及理论。第一，**收敛定理只覆盖 Galerkin，而多项模型与二维算例都用配点**，被测得最多的格式恰好是无理论的那个。第二，所有算例的瓶颈都是源项正则性，因此定理 5.1 的第一项（解自身的逼近误差，含 $N^{s-m}$）在实验里从未被真正检验；换句话说，实验证实的是理论中较弱的那一半。第三，$\lambda$ 只测两个值、$\alpha$ 只测三个值，而尺度 $\mu$ 是手工给的（$5$ 或 $3$），论文没有给出选取规则——这与编号 27 里「$M$ 难求所以最优 $r$ 难得」是同一处空白的两种形式。第四，算例 2 的参考解是自家 $N=600$ 的数值解，不是独立的精确解。

## 92：一致性而非速度

编号 92 处理的不是速度而是**一致性**：已有方法在 $\alpha$ 接近 $0$ 或 $1$ 时表现不佳。退化的来源是结构性的——矩阵分数次幂的标准积分表示带一个在整数 $\alpha$ 处为零的正弦前因子，而最佳一致有理逼近的误差渐近里出现的也是同一个正弦因子，所以两条主流路线在区间端点附近同时变差。

论文要的是代价与精度在整个 $\alpha\in(0,1)$ 上一致，处理两个对象：$A^{-\alpha}$（$A$ 正定），以及 $(q\mathcal I+A^{\alpha})^{-1}$（$q>0$）。后者正是空间分数阶演化方程隐式时间步所需的解算子。方法走求积路线（梯形规则与 Laguerre-Gauss 求积），测试问题为分数阶 Laplacian、分数阶 Poisson 方程与空间分数阶 Allen-Cahn 方程。

## 三篇的对照

| 编号 | 基 / 表示                | 远场衰减                        | 算子的处理方式             |
| ---- | ------------------------ | ------------------------------- | -------------------------- |
| 27   | Hermite 函数与过缩放基   | 高斯型                          | 物理空间的显式微分矩阵     |
| 42   | 修正映射 Gegenbauer 函数 | $\lvert x\rvert^{-(\lambda+1)}$ | 分解成两种初等有理形状     |
| 92   | 不涉及基（矩阵函数）     | 不适用                          | 积分表示加求积             |

编号 42 明确说自己是沿编号 27 的思路做的，两篇共享同一套设计模式：先求出基函数的分数阶像的闭式，从而把分数阶微分矩阵预计算下来，再在物理空间直接配点；连多项/分布阶的应用、用超几何函数的连带关系廉价填矩阵、以及引入一个尺度参数这三件事都是继承来的。**升级只在基上**：$e^{-x^2/2}H_n$ 换成 $(1+x^2)^{-(\lambda+1)/2}C^{\lambda}_n(x/\sqrt{1+x^2})$，指数尾换成幂律尾 $\lvert x\rvert^{-(\lambda+1)}$，以便对上分数阶 Laplacian 真正产生的 $\lvert x\rvert^{-\alpha-1}$。

两篇的理论状态也刚好互补。编号 27 没有收敛定理，编号 42 为 Galerkin 格式证了最优估计——但**配点格式在两篇里都没有分析**，而它恰恰是两篇都主推的格式。两篇的限制则呈镜像：编号 27 在一个代数衰减的特征值问题上退化为代数收敛，编号 42 在源项正则性成为瓶颈时以 $N^{-(\alpha+1/2)}$ 代数收敛。

编号 92 是同一个算子的离散代数一侧。编号 27 与 42 靠挑基函数让 $(-\Delta)^{\alpha/2}$ 有闭式像，编号 92 接受一个已经离散好的 $A$，直接算 $A^{-\alpha}$ 与 $(q\mathcal I+A^{\alpha})^{-1}$；后者正是编号 42 多维格式里那个 Fourier 乘子 $1/(|\xi|^{\alpha}+\rho)$ 的矩阵版本，它列出的应用（分数阶 Laplacian、分数阶 Poisson 方程）也就是编号 42 模型问题的离散形式。两条路线互补而非竞争：**能自己选离散化时选编号 42 的路，$A$ 是别人交过来的时候才轮到编号 92 的路。**

## 本页原文

- T. Tang, H. Yuan, and T. Zhou, [_Hermite spectral collocation methods for fractional PDEs in unbounded domains_](https://doi.org/10.4208/cicp.2018.hh80.12), Commun. Comput. Phys. 24(4) (2018), pp. 1143-1168（预印本 [arXiv:1801.09073](https://arxiv.org/abs/1801.09073)）。
- T. Tang, L.-L. Wang, H. Yuan, and T. Zhou, [_Rational spectral methods for PDEs involving fractional Laplacian in unbounded domains_](https://doi.org/10.1137/19M1244299), SIAM J. Sci. Comput. 42(2) (2020), pp. A585-A611（预印本 [arXiv:1905.02476](https://arxiv.org/abs/1905.02476)）。
- Y. Duan, F. Zeng, H. Zhang, and T. Zhou, [_Fast computation of the fractional power of a matrix_](https://doi.org/10.1137/25M1757411), SIAM J. Sci. Comput. 48(1) (2026), pp. A309-A334。
