---
title: 离散化的稳定性理论
description: 编号 47、63：均方意义下的 Lax 等价定理，以及第一个全离散误差分析
lang: zh
translation: en/computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes
tags:
  - 论文笔记
  - 随机微分方程
  - 稳定性理论
---

> [!note] 本页覆盖
> 编号 **47**（_SIAM J. Numer. Anal._ 58(4), 2020）、**63**（_SIAM J. Numer. Anal._ 60(4), 2022）。
>
> 两篇**均无预印本**，SIAM 的正文在本环境下不可访问（`epubs.siam.org` 被 Cloudflare 拦截），因此两篇的正文本站都未读到。本页的写法是：摘要逐字可确认的内容照写；构成其技术部件的**外部独立来源**（Chessari-Kawai-Shinozaki-Yamada 的 BSDE 数值综述、孙亚兵与赵卫东关于 Sinc 求积的 CSIAM-AM 文章、王旭与赵卫东的 Sinc 多步续篇、Tang-Xiong 与 Chassagneux 的稳定性文章）另行核实并注明出处；其余按上下文重构的形式一律标为重构。**两篇的数值实验数据本站均未取得，因此本页没有这两篇的实验表。**

## 47：把逐格式证明换成一条等价定理

### 直觉

到 2020 年，这一组手里的 FBSDE 时间离散已经堆成一个小型动物园：后向 Euler、$\theta$ 格式族、Zhao-Zhang-Ju 多步格式、[[computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|编号 8]] 的多步族、延迟校正。每一个都配一份自己的收敛证明，每一个都带自己的稳定性条件——编号 8 的 $1\le k\le6$ 是数值算出的根条件窗口，编号 23 的「没有壁垒」只是数值断言，编号 33 则先证稳定性再给误差估计。

这种局面本身就是一个信号：**当一族方法的收敛证明逐个特设时，说明还没有人把「稳定性」这个概念单独提出来定义。** 确定性数值分析里早就有现成的模板——线性发展方程的 Lax 等价定理说，相容的格式收敛当且仅当稳定。它的价值不在于让证明变短，而在于把「收敛」这个含混的目标换成「稳定」这个可以逐格式检验、甚至可以优化的性质。

这一篇做的就是把这条模板搬到 FBSDE 上。搬过来有两处非平凡：范数必须是均方（$L^2(\Omega)$）的，因为解是随机过程；问题是非线性的（生成元只假定 Lipschitz），而经典 Lax 定理讲的是线性方程。

### 问题设定

考虑解耦 Markov 型倒向方程

$$
Y_t=\varphi(X_T)+\int_t^Tf(s,X_s,Y_s,Z_s)\,\mathrm ds-\int_t^TZ_s\,\mathrm dW_s .
$$

摘要说得很清楚：所提出的是「一个用于 FBSDE 数值求解的**一般离散化框架**」，且「该框架涵盖文献中若干已有的时间离散概率型格式」。Chessari、Kawai、Shinozaki 与 Yamada 的 BSDE 数值综述对这一篇的描述与之独立吻合，并把覆盖范围写得更具体：「文献 [311] 构造了一个一般框架，以统一的方式研究 FBSDE 离散化格式的稳定性、相容性与收敛性，**包括后向 Euler 方法、$\theta$ 格式与各类多步方法**。」

### 推导

**第一步：被统一的两族格式长什么样。** 一族是 $\theta$ 格式。下面这两条递推是从综述逐字转录的（综述把它归于赵卫东、王小捷与彭实戈，以及赵卫东、李颖与张启峰）：对 $k=n-1,\dots,0$，

$$
Y^\pi_{t_k}=\mathbb E\bigl[Y^\pi_{t_{k+1}}\mid\mathcal F_{t_k}\bigr]
+\theta_1\Delta_n f\bigl(t_k,Y^\pi_{t_k},Z^\pi_{t_k}\bigr)
+(1-\theta_1)\Delta_n\,\mathbb E\bigl[f(t_{k+1},Y^\pi_{t_{k+1}},Z^\pi_{t_{k+1}})\mid\mathcal F_{t_k}\bigr],
$$

$$
Z^\pi_{t_k}=-\frac{1-\theta_2}{\theta_2}\mathbb E\bigl[Z^\pi_{t_{k+1}}\mid\mathcal F_{t_k}\bigr]
+\frac{1}{\theta_2\Delta_n}\mathbb E\bigl[Y^\pi_{t_{k+1}}\Delta W_k\mid\mathcal F_{t_k}\bigr]
+\frac{1-\theta_2}{\theta_2}\mathbb E\bigl[f(t_{k+1},Y^\pi_{t_{k+1}},Z^\pi_{t_{k+1}})\Delta W_k\mid\mathcal F_{t_k}\bigr],
$$

参数 $\theta_1\in[0,1]$、$\theta_2\in(0,1]$，初值取 $Y^\pi_{t_n}=\Phi(X^\pi_{t_n})$、$Z^\pi_{t_n}=(\nabla\Phi(X^\pi_{t_n}))\sigma(t_n,X^\pi_{t_n})$。它由精确恒等式 $Y_{t_k}=\mathbb E[Y_{t_{k+1}}|\mathcal F_{t_k}]+\int_{t_k}^{t_{k+1}}\mathbb E[f(s,Y_s,Z_s)|\mathcal F_{t_k}]\mathrm ds$ 出发，把积分用隐式（$t_k$）与显式（$t_{k+1}$）两端值的凸组合逼近而来。

**记号上有一处容易读反：$\theta_1$ 乘的是当前层 $t_k$ 的生成元，而在倒向递推里当前层才是未知的那一层。** 因此：

| $(\theta_1,\theta_2)$ | 名称与性质                                                            |
| --------------------- | --------------------------------------------------------------------- |
| $\theta_1=1$          | 左端点矩形公式，生成元隐式（类后向 Euler）                            |
| $(1,1)$               | 上式再令 $Z$ 递推塌缩为单项 $\mathbb E_n[Y_{n+1}\Delta W_n]/\Delta_n$ |
| $\theta_1=0$          | 右端点公式，完全显式                                                  |
| $(1/2,1/2)$           | Crank-Nicolson 成员，唯一的二阶成员                                   |

$(1,1)$ 这个成员正是编号 26 用在其伴随 BSDE 上的格式（见[[computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control|二阶 FBSDE 与随机控制一页]]），它名义上隐式却无需迭代，因为那里的伴随生成元在 $p$ 上线性。另一族是编号 8 的多步格式，其形式见[[computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|多步格式一页]]。

**第二步：统一格式的形状。** 一个能同时容纳这两族的 $K$ 步概率型离散化，形状大致是一对递推

$$
\sum_{j=0}^{K}a_j\,\mathbb E_{t_n}\bigl[Y^{n+j}\bigr]
=\Delta t\sum_{j=0}^{K}b_j\,\mathbb E_{t_n}\bigl[f(t_{n+j},X^{n+j},Y^{n+j},Z^{n+j})\bigr],
$$

$$
\sum_{j=0}^{K}c_j\,\mathbb E_{t_n}\bigl[Z^{n+j}\bigr]
=\frac{1}{\Delta t}\sum_{j=0}^{K}d_j\,\mathbb E_{t_n}\bigl[Y^{n+j}\Delta W_{n,j}\bigr]
+\sum_{j=0}^{K}e_j\,\mathbb E_{t_n}\bigl[f(t_{n+j},\cdots)\Delta W_{n,j}\bigr],
$$

系数向量为 $(a_j),(b_j),(c_j),(d_j),(e_j)$。**这一对递推是本站按被覆盖的两族反推的形状，不是从本文转录的**；作者实际使用的指标范围、归一化与记号本站未核实。可核实的是**覆盖声明**本身：框架涵盖后向 Euler、$\theta$ 格式与各类多步方法。

**第三步：三个概念在均方意义下的样子。** 记 $\|\cdot\|=\|\cdot\|_{L^2(\Omega)}$。下面三条定义是本站按「均方版 Lax 等价定理」这一提法重构的，**不要据此引用常数或精确范数**：

- **$p$ 阶相容性。** 把**精确解** $(Y_{t_n},Z_{t_n})$ 代入离散算子，留下局部截断残量 $(R^Y_n,R^Z_n)$，且
  $$
  \Bigl(\sum_n\|R^Y_n\|^2+\Delta t\sum_n\|R^Z_n\|^2\Bigr)^{1/2}=O\bigl((\Delta t)^p\bigr)\xrightarrow[\Delta t\to0]{}0 .
  $$
  **$Z$ 那一项多带一个 $\Delta t$ 的权重是 BSDE 里的标准做法**，它反映 $Z$ 比 $Y$ 少「半个导数」的正则性。
- **稳定性。** 把离散数据扰动 $(\delta^Y_n,\delta^Z_n)$（$k$ 步格式还包括它需要的起始值），离散解的扰动被一个**与 $\Delta t$ 无关**的常数乘以累积扰动所控制，且对 $N$ 一致：
  $$
  \max_n\|\delta Y^n\|^2+\Delta t\sum_n\|\delta Z^n\|^2\le C\bigl(\|\delta Y^N\|^2+\text{累积的数据扰动}\bigr).
  $$
- **收敛性。** $\max_n\|Y_{t_n}-Y^n\|\to0$ 与相应的 $Z$ 量趋于零（$\Delta t\to0$）；若为 $O((\Delta t)^p)$ 则称 $p$ 阶收敛。

### 定理

**主定理（内容可核实，精确陈述未核实）：一个相容的 FBSDE 离散化格式收敛，当且仅当它稳定。** 摘要的原话是「特别地，我们证明了一个**随机均方版本的 Lax 等价定理**——表明 FBSDE 的相容离散化格式收敛当且仅当稳定」。这与适定线性发展方程的经典 Lax 等价定理完全对应，差别在于范数是均方的，而问题是非线性的（生成元 Lipschitz）。

摘要还说「分析结果对已有数值格式的应用也做了讨论」，即抽象定理被实例化，用来恢复或加强后向 Euler、$\theta$ 格式与多步格式的已知收敛结果。**这正是编号 8 的经验窗口获得理论归宿的地方**：那里对 $P(\lambda)=\alpha_{k,0}\lambda^k+\sum_{j=1}^k\lambda^{k-j}$ 的根条件就是等价定理里「稳定性」的那一半，相容性由构造保证，两者合起来给出收敛。

> [!note] 「稳定」在这一族里具体指什么——三条独立核实的旁证
> 本文正文未读到，但同一时期有三项外部工作可以从侧面确定这个概念的边界，它们也都是本文对话的对象。
>
> **其一，根条件。** 编号 8 把其权用到确定性常微分方程上得 $\alpha_{k,0}Y^n+\sum_{j=1}^k\alpha_{k,j}Y^{n+j}=f(t_n,Y^n)$，特征多项式 $P(\lambda)=\alpha_{k,0}\lambda^k+\sum_{j=1}^k\lambda^{k-j}=0$ 须满足经典根条件（$|\lambda|\le1$，模为 1 的根须单）。其报出的最大根模（除去公共根 $1.0$）为 $k=2$：0.3333，3：0.4264，4：0.5608，5：0.7087，6：0.8633，7：**1.0222**，8：**1.1839**，故 $k\ge7$ 不稳定。
>
> **其二，根条件确实蕴含均方稳定性，这一点由外部工作证明。** Tang 与 Xiong（_IMA J. Numer. Anal._ 42(2) (2022) 1789-1805）的摘要逐字说：「在经典根条件下，我们证明一般线性多步方法对生成元同时依赖 $y$ 与 $z$ 的解耦 FBSDE 是**均方（零）稳定**的。基于该稳定性结果，我们进一步建立了一条基本收敛定理。」**「均方（零）稳定」几乎肯定就是本文所形式化的概念。**
>
> **其三，同样的架构在别处也出现过。** Chassagneux（_SIAM J. Numer. Anal._ 52(6) (2014) 2815-2836）的摘要逐字说：「我们证明，在关于系数的一个充分条件下，这些格式享有一条**基本的稳定性性质**。把这一结果与截断误差分析结合，使我们能够设计具有任意收敛阶的逼近。」这正是「稳定性加截断误差推出阶」的同一套架构，而编号 47 把它公理化了。
>
> 另外两处已经记录过的具体稳定窗口可以作为标定：Zhao-Zhang-Ju 的 $Y$ 方程只在 $K_y\in\{1,\dots,7,9\}$ 稳定、$Z$ 方程只在 $K_z\in\{1,2,3\}$ 稳定；编号 8 与 19 都是 $1\le k\le6$。

> [!note] 可核实范围
> 主定理的**内容**（相容加稳定等价于收敛）与覆盖声明（后向 Euler、$\theta$ 格式、各类多步格式）可从摘要确认，并由 Chessari 等人的综述独立印证。一般离散化族的确切写法、稳定性与相容性的精确范数定义、显式常数、精确假设（大概率是 $f$ 关于 $(y,z)$ 一致 Lipschitz、终值平方可积、正向扩散的正则性），以及框架是否覆盖**耦合**而非仅解耦系统，本站均未核实。上面给出的 $\theta$ 格式是从独立的第三方综述转录的标准形式，统一格式的形状则是本站反推的。

### 数值实验

**摘要没有提到任何数值实验**，只提到「分析结果对已有数值格式的应用」。因此本页不报告这一篇的任何实验数据。从摘要判断，这很可能是一篇纯分析的论文，但这一点本站未核实——正文未读到，无法排除文中另有算例。

### 与其他论文的关系

**这是整个系列的理论拱心石。** 编号 8 与 19 提出多步格式并给出截断误差估计，但没有收敛定理；编号 23 与 35 提出延迟校正格式，同样没有 FBSDE 收敛定理；编号 33 对一个特定族先证稳定性再给误差估计。编号 47 把这个反复出现的模式抽象成框架，并证明稳定性加相容性等价于收敛性。它直接解释了编号 8 与 19 中经验观察到、并在编号 25 中被数值断言的 $k\le6$ 壁垒。

**它使编号 68 的反向设计成为可能**：一旦稳定性被确立为枢纽性质，「设计一个收敛的格式」就变成「设计一个稳定的格式」，而后者是一个有明确目标的优化问题——编号 68 正是据此给出关于系数的新充分条件，并优化出一到五阶的强稳定保持格式。编号 63 则把分析从**半离散**（只算时间，假设条件期望精确）推进到**全离散**（时间加空间），两篇共用「先稳定性」的架构。第一作者杨杰与编号 35、61 相同；杨杰与赵卫东更早的《Convergence of recent multistep schemes for a forward-backward stochastic differential equation》（EAJAM 5 (2015) 387-404，见于编号 41 的参考文献）是它的前身。

## 63：第一个全离散误差分析

### 直觉

每个概率型 BSDE 格式有两个误差来源：倒向方程的**时间**离散，以及条件期望 $\mathbb E_{t_n}[\cdot]$ 的**空间**逼近——在 Markov 情形下它们是 $d$ 维高斯积分，而被积函数只在网格点上已知。

这一组此前的流水线（编号 8、19、25、41）用 Gauss-Hermite 求积算这些积分。它的节点是 $x_n+\sqrt{\Delta t}\,\lambda_j$，其中 $\lambda_j$ 是 Hermite 多项式的根，**一般是无理数，不会落在空间网格上**，因此每个节点、每个时间步都必须做一次局部多项式插值。这一步有三重代价：它本身花时间；它把空间精度封顶在插值阶；它使**全离散**误差分析极其困难。后果是，在本文之前，几乎所有严格的 BSDE 误差分析都是**半离散**的——假设条件期望精确，只界定时间误差。

**这一篇的机制可以用一句话说清：换一个求积公式，让它的节点自己落到网格上。** Sinc 求积的节点是等距的 $kh$，而且 $h$ 是一个**自由参数**。于是只要把 $h$ 调成让节点间距等于空间步长，插值这一步就整个消失了。附带的好处是 Sinc 求积对适当类的函数是**指数收敛**的，比插值阶高得多。

### 问题设定

标准（一阶）BSDE，**没有** $\Gamma$ 过程——这一篇处理的不是 2FBSDE。论文的自述是：这似乎是**首次尝试分析 BSDE 的全离散格式**，达到时间二阶收敛与空间指数收敛。zbMATH 关键词为「倒向随机微分方程；误差估计；条件数学期望；Sinc-$\theta$ 格式」，MSC 为 60H10、60H35、65C30。

### 推导

**第一步：时间方向用 $\theta$ 格式族。** 摘要逐字确认「我们对时间离散考虑 $\theta$ 格式」。该族的两条递推与上一节编号 47 处所引的相同（同样转录自 Chessari 等人的综述）。综述同时记录：该格式在 $\theta_1=\theta_2=1/2$ 时二阶收敛，否则一阶——**这与本文声称的「时间二阶收敛率」正好对上，即所用的是 Crank-Nicolson 成员。**

阶数这一处有一个值得记清的细节。赵卫东、王小捷与彭实戈（_DCDS-B_ 12(4) (2009) 905-924）研究的是**单参数** $\theta$ 格式与**不依赖 $z$** 的生成元，其摘要说：一般 $\theta$ 时对 $y$ 一阶收敛；$\theta=\frac12$ 时对 $y$ 二阶、**对 $z$ 一阶**。要在生成元依赖 $(y,z)$ 时对 $Y$ 与 $Z$ **同时**达到二阶，需要赵卫东、李颖与张启峰（_DCDS-B_ 17(5) (2012) 1585-1603）「引入更多参数」的推广版本。**这两篇都在编号 63 的引用清单里**，这正是它的「时间二阶」与 $\theta_1=\theta_2=1/2$ 相容的原因。

**第二步：空间方向用 Sinc 逼近。** 摘要逐字确认「随后采用 Sinc 逼近来逼近相关的条件数学期望」。下面的构造转录自**同一组人**的另一篇文章（孙亚兵与赵卫东，_CSIAM Trans. Appl. Math._ 6(1) (2025) 176-206，可免费下载，其 4.1 节）：

$$
\mathrm{sinc}(x)=
\begin{cases}
\dfrac{\sin(\pi x)}{\pi x},&x\neq0,\\[6pt]
1,&x=0 .
\end{cases}
$$

设 $B(h)$ 是这样一类整函数 $g$：在实轴上 $g\in L^2(\mathbb R)$，且对一切 $z\in\mathbb C$ 有 $|g(z)|\le K\exp(\pi|z|/h)$。**Whittaker 基数展开（Stenger）**：若 $g\in B(h)$，则 $g(z)=\sum_{k=-\infty}^{\infty}g(kh)\,\mathrm{sinc}\bigl(\frac{z-kh}{h}\bigr)$；若 $\sum_kg(kh)$ 收敛，则对充分小的 $h$ 有 $\int_{\mathbb R}g(x)\mathrm dx=h\sum_{k=-\infty}^{\infty}g(kh)$。据此定义 **Sinc 求积公式**

$$
T_M(g,h)=h\sum_{k=-M}^{M}g(kh),
\qquad
\eta_M(g,h)=\int_{\mathbb R}g(x)\,\mathrm dx-T_M(g,h),
$$

其误差定理是：若 $g$ 有界，则对充分小的 $h$ 与任意满足 $\gamma_0\le Mh^2$ 的 $\gamma_0>0$，

$$
|\eta_M(g,h)|\ \le\ C_{\gamma_0,g}\,h\,\exp\Bigl(-\frac{M^2h^2}{2}\Bigr),
$$

常数 $C_{\gamma_0,g}$ 只依赖 $\gamma_0$ 与 $\|g\|_{L^\infty}$。**这就是「空间指数收敛」的来源。**

**第三步：用到高斯条件期望上。** 对 $X^{t,x}_r=x+\sigma_0(W_r-W_t)$ 与光滑的 $v:\mathbb R^d\to\mathbb R$，

$$
\mathbb E\bigl[v(X^{t,x}_r)\bigr]
=\int_{\mathbb R^d}v\bigl(x+\sigma_0\sqrt{r-t}\,p\bigr)\Bigl(\tfrac{1}{\sqrt{2\pi}}\Bigr)^{d}e^{-p^\top p/2}\,\mathrm dp
\ \approx\ \sum_{\mathbf k=-M}^{M}v\bigl(x+\sigma_0\sqrt{r-t}\,h\mathbf k\bigr)\,\alpha_{\mathbf k},
$$

$$
\alpha_{\mathbf k}=\prod_{i=1}^{d}\alpha_{k_i},
\qquad
\alpha_{k_i}=\frac{h}{\sqrt{2\pi}}\exp\Bigl(-\frac{k_i^2h^2}{2}\Bigr),
\qquad
\beta^M_{\mathbf k}=\frac{\alpha_{\mathbf k}}{\sum_{\mathbf k=-M}^{M}\alpha_{\mathbf k}} .
$$

$\beta^M_{\mathbf k}$ 是**重新归一化**后的权，它们和为 1——这使所得矩阵的行和为 1，从而保持有界性。**（上式转录自孙亚兵与赵卫东的文章；编号 63 中相应的一步本站未核实，是按此类推的。）**

**第四步：为什么 Sinc 消掉了插值。** 摘要逐字确认「通过恰当选择 Sinc 求积公式中的参数，我们的格式被证明非常高效，因为**不需要任何空间插值**」。机制在孙亚兵与赵卫东的 Remark 4.2 里说得最直白（逐字转录）：

> 「除了谱精度之外，我们选择 Sinc 求积公式而非其他求积公式的主要原因是：它的求积节点是**均匀的**，且含有一个**自由参数 $h$**。因此，对不同中间时刻之间的不同时间增量设定不同的 $h$ 值，我们就能得到不同的均匀求积节点，从而在逼近条件期望时**避免使用插值**。」

具体地说：对应增量 $r-t$ 的 Sinc 节点位于 $x+\sigma_0\sqrt{r-t}\,h\mathbf k$，是一个间距为 $\sigma_0\sqrt{r-t}\,h$ 的**等差数列**。把 $h$ 选成使 $\sigma_0\sqrt{r-t}\,h=\Delta x$（空间步长），每个求积节点就**恰好落在网格点上**。**Gauss-Hermite 做不到这件事，因为它的节点是 Hermite 多项式的根。** 同一批作者的 Sinc 多步续篇（王旭与赵卫东，_Adv. Appl. Math. Mech._ 15(3) (2023)）把同一机制写得更明确：「通过使用**积分变量变换**并在 Sinc 求积公式中**恰当选择空间步长参数**，不需要任何空间插值。」

**第五步：合起来是什么。** 把 Sinc 公式代入 $\theta$ 格式中的四个条件期望，每一个都变成有限的、与网格对齐的加权和，于是格式成为网格值 $\{Y^n_{\mathbf i}\},\{Z^n_{\mathbf i}\}$ 上的**矩阵递推**，转移矩阵带状／Toeplitz 且行和为 1。**这一段是本站的重构**：本文是否写成矩阵形式、空间定义域的截断如何处理（Sinc 和在 $\pm M$ 处截断，这本身就是定义域截断）、以及 $h,M,\Delta t,\Delta x$ 之间被强加的确切关系，本站均未核实。

$Z$ 的处理值得单独一提（同样是重构）：如上面的 $\theta$ 格式所示，$Z^n$ 来自 $\mathbb E_n[Y^{n+1}\Delta W_n]/(\theta_2\Delta t)$ 加上含 $\mathbb E_n[Z^{n+1}]$ 与 $\mathbb E_n[f^{n+1}\Delta W_n]$ 的修正项，**全程不对插值函数求导**。这一点重要，因为对插值函数求导正是 $Z$ 的精度在别处的常见瓶颈——延迟校正（编号 23）付的代价就在这里。

### 定理

- **主结论（摘要逐字）：** 「我们进行了**严格的稳定性分析与误差估计**，这似乎是首次尝试分析 BSDE 的**全离散**格式，达到**时间二阶收敛率**与**空间指数收敛**。」
- 因此头条估计的形状大约是
  $$
  \max_n\bigl\|Y_{t_n}-Y^{n,M,h}\bigr\|\ \le\ C\Bigl(\underbrace{(\Delta t)^2}_{\theta_1=\theta_2=1/2}+\underbrace{C'\,h\,e^{-M^2h^2/2}}_{\text{Sinc}}\Bigr),
  $$
  可能还带一项定义域截断误差。**这个形状是本站由两个已核实的速率拼出来的，不是转录**；确切陈述、范数、常数、所需正则性、同样的速率对 $Z$ 是否成立，以及 $M,h,\Delta t$ 之间的精确耦合条件，本站均未核实。
- **稳定性分析确实做了**（摘要可确认），但「稳定」在这里的精确含义本站未核实。鉴于本文引用编号 47，它几乎肯定就是那个框架里的均方稳定性概念，只是现在扩展到包含**求积引入的空间扰动**——这一扩展正是本文的技术新意所在。

> [!note] 一条独立核实的姊妹结果，可用作标定
> Sinc **多步**续篇（王旭与赵卫东，AAMM 2023）的摘要逐字说：「$K$ 步 Sinc 多步格式的稳定性与**时间 $K$ 阶误差估计**得到了理论证明（$1\le K\le6$）。这似乎是首次分析 FBSDE 的**全时空离散多步格式**。」其贡献清单还补充：最优 $K$ 阶是「对漂移项 $b=0$、扩散项 $\sigma$ 为非奇异矩阵的 FBSDE (1.1)」，通过「运用**数值代数理论**与逼近算子的性质」得到。**这里的 $1\le K\le6$ 正是编号 8 的根条件窗口**——同一个壁垒在全离散分析里再次出现，说明它来自时间模板而不是空间逼近。

> [!note] 可核实范围
> 「首次全离散分析」「时间二阶、空间指数收敛」「用 $\theta$ 格式做时间离散、用 Sinc 逼近处理条件期望」「不需要空间插值」均可从摘要逐字确认。$\theta$ 格式的具体形式转录自第三方综述；Sinc 求积的定义、误差定理与用于高斯条件期望的写法转录自同组的 CSIAM-AM 文章；「均匀节点加自由参数因而免插值」这条机制转录自同一篇的 Remark 4.2 与 Sinc 多步续篇的摘要。**本文自身的正文本站未读到**，因此其定理的精确陈述、常数、假设，以及全离散格式的具体写法都不报告。

### 数值实验

摘要说「给出了若干数值实验以验证理论结果并展示我们格式的效率与精度」。**具体的测试问题与观察到的阶本站未核实，因此本页不列表。** 按理论结果可以预期的是：一张时间加密表显示 $\theta_1=\theta_2=1/2$ 时率约为 2、其余情形约为 1；一张关于 $M$ 或 $h$ 的半对数图显示空间误差呈直线（指数）下降；以及与「Gauss-Hermite 加插值」格式的 CPU 时间对比。**但这三项本站一项也没有核实，不应当作事实陈述。**

### 与其他论文的关系

**这是整个研究计划在空间离散方向上的对应物。** 编号 8、19、25、41 都用 Gauss-Hermite 求积加局部 Lagrange 插值逼近条件期望；编号 63 用单一的 Sinc 求积同时替掉两者，并让节点与网格重合。它是这份清单里第一篇证明**全离散**（时间加空间）误差估计的论文。

它引用编号 47 并沿用其「先稳定性」的架构，也引用编号 8、18、25 与 33，即把自己放在整条多步／$\theta$ 线索的对面。编号 68 引用它（作为其参考文献 39），两篇 2022 至 2023 年的文章在稳定性主题上是同伴。编号 93（深度随机差分方法）也引用它，说明经典的 Sinc 工作反过来喂给了深度学习那一支。它的直接续篇是王旭与赵卫东的《Sinc-Multistep Schemes for FBSDEs》（AAMM 2023），把 Sinc 空间规则与编号 8 式的多步时间规则结合，并证明 $1\le K\le6$ 时的 $K$ 阶。

**值得记一笔的外溢：** 同一套 Sinc 机器后来被孙亚兵与赵卫东用于 Allen-Cahn 方程的**保最大值原理**随机 Runge-Kutta 格式（CSIAM-AM 2024/2025，第一与第二部分），这是 BSDE 数值向结构保持型抛物 PDE 求解器的一次跨越。周涛也是该线索所引用的 BDF2／Allen-Cahn 文章（廖洪林、汤涛与周涛，_SIAM J. Numer. Anal._ 58 (2020) 2294-2314）的作者之一。

## 两篇的关系

| 编号 | 覆盖的误差来源 | 主要结论                      | 对其他篇的作用                                            | 本站核实程度         |
| ---- | -------------- | ----------------------------- | --------------------------------------------------------- | -------------------- |
| 47   | 时间           | 均方 Lax 等价定理             | 给编号 8 的经验窗口理论归宿；使编号 68 的反向设计成为可能 | 摘要 + 综述独立印证  |
| 63   | 时间 + 空间    | 首个全离散分析（二阶 + 指数） | 补上编号 25 之外的另一条空间方向改进                      | 摘要 + 部件独立核实  |

一条一般判断：**当一族方法的收敛证明逐个特设时，值得停下来把「稳定性」定义清楚。** 编号 47 之后，「设计一个收敛的格式」变成「设计一个稳定的格式」，而后者是一个可以优化的目标——[[computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|编号 68]] 正是这样做的。

另一条判断关于编号 63：**它改的不是精度而是节点的位置。** Gauss-Hermite 的精度并不差，问题在于它的节点由多项式的根决定，无法与网格对齐，于是逼出一个插值步，而这个插值步同时吃掉时间、精度与可分析性。换成节点等距且间距可调的求积公式，三样代价一起消失。**这是一个「换工具而非加工具」的改进，在数值分析里比看上去更少见。**

## 覆盖核对

| 内容                                              | 论文 | 覆盖状态                                     |
| ------------------------------------------------- | ---- | -------------------------------------------- |
| $\theta$ 格式族的两条递推与四个具名成员           | 47   | 转录自第三方综述                             |
| 统一格式的形状                                    | 47   | 本站反推，已标明                             |
| 相容性、稳定性、收敛性的均方定义                  | 47   | 本站重构，已标明                             |
| 均方 Lax 等价定理与覆盖声明                       | 47   | 摘要逐字，并由综述独立印证                   |
| 「稳定」的三条外部旁证（根条件、Tang-Xiong、Chassagneux） | 47 | 逐条独立核实                             |
| 无数值实验                                        | 47   | 摘要未提；已注明                             |
| $\theta$ 格式的阶与两篇 DCDS-B 的分工             | 63   | 独立核实                                     |
| Sinc 求积的定义、误差定理与高斯条件期望的写法     | 63   | 转录自同组的 CSIAM-AM 文章                   |
| 「均匀节点加自由参数因而免插值」的机制            | 63   | 逐字转录 Remark 4.2 与续篇摘要               |
| 全离散估计的形状与稳定性分析                      | 63   | 两个速率可核实；形状为本站拼出，已标明       |
| Sinc 多步续篇的 $1\le K\le6$ 与 $b=0$ 条件        | 63   | 独立核实                                     |
| 数值实验                                          | 63   | 摘要提到存在；数据未核实，未列表             |

## 本页原文

- J. Yang, W. Zhao, and T. Zhou, [_A unified probabilistic discretization scheme for FBSDEs: stability, consistency, and convergence analysis_](https://doi.org/10.1137/19M1260177), SIAM J. Numer. Anal. 58(4) (2020), pp. 2351-2375。
- X. Wang, W. Zhao, and T. Zhou, [_Sinc-theta schemes for backward stochastic differential equations_](https://doi.org/10.1137/21M1444679), SIAM J. Numer. Anal. 60(4) (2022), pp. 1799-1823。
- 用于交叉印证的外部来源：C. Chessari, R. Kawai, Y. Shinozaki, and T. Yamada, [_Numerical methods for backward stochastic differential equations: a survey_](https://doi.org/10.1214/23-PS18), Probab. Surveys 20 (2023), pp. 486-567（[arXiv:2101.08936](https://arxiv.org/abs/2101.08936)；本站据此转录 $\theta$ 格式与对编号 47 的独立描述）；Y. Sun and W. Zhao, _Stochastic Runge-Kutta methods for preserving maximum bound principle of semilinear parabolic equations, part II: Sinc quadrature rule_, [CSIAM Trans. Appl. Math. 6(1) (2025), pp. 176-206](https://doi.org/10.4208/csiam-am.SO-2024-0012)（本站据此转录 Sinc 求积的定义、误差定理与免插值机制）；X. Wang and W. Zhao, [_Sinc-multistep schemes for forward backward stochastic differential equations_](https://doi.org/10.4208/aamm.OA-2022-0073), Adv. Appl. Math. Mech. 15(3) (2023)（编号 63 的多步续篇）；X. Tang and J. Xiong, [_Stability analysis of general multistep methods for Markovian backward stochastic differential equations_](https://doi.org/10.1093/imanum/drab023), IMA J. Numer. Anal. 42(2) (2022), pp. 1789-1805；J.-F. Chassagneux, [_Linear multistep schemes for BSDEs_](https://doi.org/10.1137/120902951), SIAM J. Numer. Anal. 52(6) (2014), pp. 2815-2836；W. Zhao, J. Wang, and S. Peng, [_Error estimates of the theta-scheme for backward stochastic differential equations_](https://doi.org/10.3934/dcdsb.2009.12.905), Discrete Contin. Dyn. Syst. Ser. B 12(4) (2009), pp. 905-924；W. Zhao, Y. Li, and G. Zhang, _A generalized theta-scheme for solving backward stochastic differential equations_, Discrete Contin. Dyn. Syst. Ser. B 17(5) (2012), pp. 1585-1603。
