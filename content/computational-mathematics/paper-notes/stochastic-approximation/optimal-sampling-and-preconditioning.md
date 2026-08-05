---
title: 最优采样与预条件
description: 编号 22、24、28、45、54：把「采样密度」与「正交性密度」解耦
lang: zh
translation: en/computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning
tags:
  - 论文笔记
  - 不确定性量化
  - 最优采样
---

> [!note] 本页覆盖
> 编号 **22**（_Math. Comput._ 86, 2017）、**24**（_SIAM J. Sci. Comput._ 39(3), 2017）、**28**（_SIAM J. Sci. Comput._ 40(1), 2018）、**45**（_SIAM Rev._ 62(2), 2020）、**54**（_J. Comput. Phys._ 430, 2021）。

![配点设计的统一流程](assets/diagrams/tao-zhou-papers/zh/sampling-design.svg)

## 贯穿整页的那一个念头

这一族论文全部建立在同一句话上：**你所采样的密度，不必是定义正交性的那个密度。**

把这句话摊开来说。设 $V$ 是 $L^2_w$ 中的 $N$ 维空间，$\{v_n\}$ 是它的 $w$-正交归一基，最小二乘的法方程矩阵是秩一矩阵的平均

$$
G=\frac1M\sum_{m=1}^{M}v_mv_m^{T},
\qquad
v_m=\frac{1}{\sqrt M}\bigl(v_1(x_m),\dots,v_N(x_m)\bigr)^{T} .
$$

无论怎样采样，$\mathbb E G=I$。真正决定「要多少个样本 $G$ 才接近 $I$」的，是这些秩一项**大小是否均匀**：若某些点上 $\sum_n v_n^2(x)$ 比别处大出几个量级，那么绝大多数样本对 $G$ 几乎没有贡献，而少数样本主宰一切，于是 $M$ 必须很大才能让平均稳定下来。这个最大值

$$
\bigl\|K\bigr\|_\infty=\sup_{x\in D}\sum_{n=1}^{N}v_n^2(x)
$$

就是**稳定性因子**。从正交性密度 $w$ 采样时，它常常随多项式次数爆炸，这正是标准 Monte Carlo 需要超线性样本数的原因。

**要把各行拉平，有两个互为倒数的动作**：一是在行范数大的地方多放样本（换采样密度），二是把大的行按其大小缩小（加权）。两者以恰好互逆的幅度同时施行时，$M\to\infty$ 的极限完全不变——这只是一次重要性采样的变量替换——改变的只有**前渐近**行为，而前渐近行为就是全部争论所在。做这件事所需的那把尺子，就是**在各点度量该函数空间集中程度**的量，即（倒数）Christoffel 函数。

后面五篇论文的分野只在于把这把尺子用在哪一侧：编号 22、24、28 从平衡测度采样、用 Christoffel 函数加权（在次数趋于无穷的意义下最优）；编号 45 干脆把倒数 Christoffel 函数本身当作采样密度（在每个有限 $N$ 上都最优）；编号 28 更进一步，把随机性整个换成贪心选点。

### 本页的记号约定

本页统一取

$$
K(z)=K_\Lambda(z)=\sum_{\alpha\in\Lambda}\varphi_\alpha^2(z),
\qquad
\lambda_\Lambda(z)=\frac{1}{K_\Lambda(z)},
\qquad
N\lambda_\Lambda(z)=\frac{N}{K_\Lambda(z)},
$$

即 $K$ 是**正交归一基平方和**，而（归一化的）**Christoffel 函数**是 $N/K(z)$，所以权总以 $1/K$ 的形式出现。$K$ 与基的选取无关：写 $\varphi$ 为基向量，则 $K(z)=\varphi^T\varphi$，任意正交变换 $\psi\leftarrow U\varphi$ 都不改变它；等价地说，$\sum_n v_n^2$ 是子空间 $V$ 的性质而非某组基的性质。

原文之间的记号并不统一，这是阅读这一族论文时最常见的混淆来源：

| 论文 | 记号              | 含义                                                                                                  |
| ---- | ----------------- | ----------------------------------------------------------------------------------------------------- |
| 22   | $K$、$K_k$        | $\sum_\alpha\varphi_\alpha^2$，即 Christoffel 函数的倒数；称 $N/K_k$ 为「（归一化）Christoffel 函数」 |
| 24   | $\lambda_\Lambda$ | $1/\sum_{i\in\Lambda}\varphi_i^2$，即 Christoffel 函数本身；预条件用 $N\lambda_\Lambda$               |
| 28   | $K_\Lambda$       | $\sum_\alpha\psi_\alpha^2$，倒数；所加的权是 $1/\sqrt{K_\Lambda}$                                     |
| 45   | $q^2$             | $\frac1N\sum_n v_n^2$，倒数且归一化；直接用作采样偏置                                                 |
| 54   | $K$、$\Lambda_N$  | $K$ 在这里指**核函数本身**（矩阵写作 $\mathbf A=K(\Xi,\Xi)$），$\Lambda_N$ 指 **Lebesgue 常数**——两个记号都与上面各行相撞，读到这一篇须整体切换 |

## 22：稳定性因子是一个可以设计掉的东西

### 直觉

标准做法是从 $w$ 采样、并在 $L^2_w$ 中做最小二乘，看上去天经地义，因为「采样」与「正交性」用了同一个 $w$。但这两件事其实各司其职：$w$ 决定谁与谁正交，采样密度决定信息从哪里来。把它们绑在一起，就等于放弃了唯一一个可以自由调节的旋钮。

**论文的洞察是：稳定性因子是这种错配造成的假象。** 只要从多势论平衡测度采样、再用 Christoffel 函数重新加权，有效的稳定性因子就降到它的最优值 $N$。用一句更贴近实现的话说：**按 Christoffel 函数加权等价于把设计矩阵的各行归一化**——这是论文自己给出的等价解读，它把「最优采样」翻译成一个纯线性代数的操作，也解释了为什么行范数一致是关键。

不这样做会怎样：对许多常用权而言，$\|K\|_\infty/N$ 随多项式次数爆炸，样本数被迫超线性增长，常常远差于二次。

### 问题设定

设 $K_k$ 是多项式空间 $\mathbb P_k$ 在 $L^2_w$ 中的再生核对角线，$N=\dim\mathbb P_k$。量 $N/K_k(z)$ 即（归一化的）**Christoffel 函数**，方法由此得名。

CLS 实际上是在带权 $\tilde w$ 的 $L^2$ 中做逼近，

$$
\tilde w(z)\triangleq\frac{N}{K(z)}\,v(z),
$$

其中 $v$ 是平衡测度的 Lebesgue 密度。因此全部理论都在度量 $\tilde w$ 离 $w$ 有多远，度量它的有两个**差异对象**：与函数无关的，是 $w$-正交归一基在 $\tilde w$ 下的 Gram 矩阵

$$
(R)_{m,n}=\int_D\varphi_m(z)\varphi_n(z)\,\tilde w(z)\,\mathrm dz;
$$

与函数有关的，是 $d(f)=\|\tilde\Pi f-\Pi f\|_w$，其中 $\Pi$、$\tilde\Pi$ 分别是到 $P$ 上的 $L^2_w$-、$L^2_{\tilde w}$-正交投影。$f\in P$ 时 $d(f)=0$；$\tilde w=w$ 时对一切 $f$ 都为零。

### 推导

第一步是写下标准 Monte Carlo 的判据（定理 4.1）：样本要求是 $\frac{S}{N\log S}$ 不小于常数乘以 $\frac{\|K(z)\|_\infty}{N}$。**一切都压在这个因子上。**

第二步是一条渐近恒等式（论文式 (2)）：若 $D$ 紧、内部非空且具有正的 $d$ 维 Lebesgue 测度，$w$ 在内部连续并容许一族正交多项式，则

$$
\lim_{k\to\infty}\frac{N}{K_k(z)}=\frac{w(z)}{v(z)}\quad\text{在 }D\text{ 中几乎处处成立},
$$

$v$ 为 $D$ 的多势论平衡测度的 Lebesgue 密度。一维 $D=[-1,1]$ 上，$v$ 就是反正弦（Chebyshev）密度。

第三步是把这条恒等式直接用成一个构造。定义非多项式函数

$$
\psi_n(z)=\sqrt{\frac{N}{K_k(z)}}\,\varphi_n(z),
$$

它们是 $\frac{1}{\sqrt{K_k}}\mathbb P_k$ 的一组基，在修正权 $\frac{wK_k}{N}$ 下正交；由上面的极限，它们对 $v$ **近似正交归一**。而它们的「近似再生核对角线」是

$$
\tilde K_k=\sum_n\psi_n^2=\frac{N K_k}{K_k}=N,
$$

**恰好取到上确界所能取到的最小值 $N$**。于是：从 $v$ 采样、用这组函数做最小二乘，就正好命中最好的样本数判据。论文自己的话是，这就是 CLS 算法的全部内容。

第四步是把算法写出来。紧域 $D$ 上（论文的算法 2）：

1. 从平衡测度 $\mu_D$ 独立同分布抽 $S$ 个样本 $\{z_s\}$；
2. 组装数据向量 $u$，$(u)_s=u(z_s)$；
3. 计算最小二乘权 $K$，$(K)_{s,s}=N/K(z_s)$；
4. 组装 $S\times N$ 的类 Vandermonde 矩阵 $V$，$(V)_{s,n}=\varphi_n(z_s)$；
5. 解 $c=\arg\min_{g}\bigl\|\sqrt{K}Vg-\sqrt{K}u\bigr\|$。

算法 3 是无界域的变体。这里有一个值得记住的细节：对有界 $D$，采样密度 $v(z)=v_D(z)=\frac{\mathrm d\mu_D}{\mathrm dz}$ **与正交性密度 $w$ 无关**——在区间上，无论 $w$ 是什么，CLS 都规定 Chebyshev 采样。论文把这一点与「区间上 Chebyshev 测度是普适的」这一经验说法联系起来。$D=[-1,1]^d$ 时 $v_D$ 是一维反正弦测度的乘积；单位球等特殊域有显式公式；无界锥形域上 $v$ 是 $\sqrt w$-加权多势论平衡测度的一个缩放。

第五步是代价核算：既然实际逼近发生在 $L^2_{\tilde w}$ 而非 $L^2_w$，误差里就必然留下 $R$ 与 $d(f)$ 的痕迹。这正是下面主定理最后一项的来源。

### 定理

**定理 4.1（Cohen–Davenport–Leviatan 的重述；标准 MC 基准）。** 设 $P$ 是 $L^2_w$ 的任意 $N$ 维子空间，$\{z_s\}_{s=1}^S$ 从 $w$ 独立同分布。若

$$
\frac{S}{N\log S}\ \ge\ \Bigl[\frac{1+r}{c_\delta}\Bigr]\frac{\|K(z)\|_\infty}{N},
\qquad c_\delta\triangleq\delta+(1-\delta)\log(1-\delta),\ \delta\in(0,1),\ r>0,
$$

则离散 Gram 矩阵满足 $\Pr\bigl[|\!|\!|G-I|\!|\!|>\delta\bigr]\le\frac{2}{S^r}$，并附带一条对满足 $\|f\|_\infty\le L$ 的 $f$ 的精度估计。

**定理 4.2（Berman、Bloom–Levenberg 等的重述；Christoffel 渐近）。** 设 $D\subset\mathbb R^d$ 势论可容许，$\mathrm dV(z)=q(z)\mathrm dz$ 光滑，$\rho$ 有界连续且使 $\rho\,\mathrm dV$ 在 $L^2_{q\rho}(D)$ 中定义一族正交归一函数，$K^{(k)}_k$ 为 $\mathbb P_k$ 在 $L^2_{q\rho^{2k}}$ 中的再生核对角线，则

$$
\lim_{k\to\infty}\frac{1}{N}\rho^{2k}(z)K^{(k)}_k(z)\,\mathrm dV(z)=\mathrm d\mu_{D,Q}(z)\quad\text{（弱收敛）} .
$$

**推论 4.1。** (1) 有界情形：$D$ 紧连通、$w$ 连续且容许正交基，取 $\rho\equiv1$、$\mathrm dV=w\,\mathrm dz$，则 $\lim_{k\to\infty}\frac1NK_k(z)=\frac{\mathrm d\mu_D}{\mathrm dV}=\frac{v_D(z)}{w(z)}$。(2) 无界凸锥 $D$、$w=\exp(-2Q)$、$\rho=\sqrt w$、$\mathrm dV=\mathrm dz$，则 $\lim_{k\to\infty}\frac1N\rho^{2k}(z)K^{(k)}_k(z)=v_{D,Q}(z)$。

**定理 4.3（最优测度，取自 Bos 等）。** 若 $\mu_k$ 是 $\mathbb P_k$ 在 $D$ 上带权 $\rho=\exp(-Q)$ 的最优测度，则 (i) $\kappa_{k,\rho}=N$ 在 $\mu_{D,Q}$ 意义下几乎处处成立；(ii) $k\to\infty$ 时 $\mu_k$ 弱收敛到 $\mu_{D,Q}$。由定理 4.1，这意味着**渐近地取到简单的对数线性标度 $S\log S\gtrsim N$，即可能的最好判据**。

> [!warning] 定理 4.3 的三条限制，论文自己列出
> 其一，该结果在**固定的 $k$ 上不给任何最优性**，因此 $k$ 较小时按平衡测度采样可能相当次优；其二，$k$ 要相对 $d$ 多大才能让渐近生效并不清楚；其三，高维下 $k$ 大在计算上不可行，因为 $\dim\mathbb P_k\sim k^d$。这三条正是编号 45 的诱导采样要消除的缺口。

**定理 5.1（CLS 稳定性）。** 对紧域 $D$ 与可容许的 $w$、指标集 $\Lambda$，若

$$
\frac{S}{N\log S}\ \ge\ \Bigl[\frac{1+r}{c_\delta}\Bigr]\frac{1}{\lambda_{\min}(R)},
$$

则 CLS 的离散 Gram 矩阵 $G$ 以高概率接近 $R$。**与定理 4.1 的对比是这条线索的核心：$\|K\|_\infty/N$ 被换成了 $1/\lambda_{\min}(R)$。**

**定理 5.2（CLS 精度，本文主定理）。** 设 $D$ 紧，CLS 即从平衡测度 $v$ 独立同分布采样并取权 $N/K(z)$。对任意 $r>0$，只要

$$
\frac{S}{N\log S}\ \ge\ C\,\frac{1+r}{\lambda_{\min}(R)},
$$

该过程即以高概率稳定（$C$ 为绝对常数）。进一步，若 $|f|\le L$，以 $T_L(x)=\mathrm{sgn}(x)\min\{|x|,L\}$ 为截断函数、$\tilde\Pi_Sf$ 为 $S$ 样本的 CLS 估计，则

$$
\mathbb E\Bigl[\bigl\|f-T_L(\tilde\Pi_Sf)\bigr\|_w^2\Bigr]
\le\|f-\Pi f\|_w^2
+\frac{\varepsilon(S)}{\lambda_{\min}(R)}\|f-\Pi f\|_{\tilde w}^2
+\frac{8L^2}{S^r}
+4\kappa^2(R)\,d^2(f),
$$

其中 $\varepsilon(S)=\frac{2-2\log2}{(1+r)\log S}\to0$，$\kappa(R)=\lambda_{\max}(R)/\lambda_{\min}(R)$。

无界域有对应的定理 5.3 与 5.4。论文说分析本身并无大的困难，**困难在实现**：对所关心的那些权，加权多势论平衡测度没有已知的显式表达式。论文为此**猜测**了若干形式，并用模拟支持这些猜测（第 6 节、表 2）。这些形式应当按猜测对待，不是定理。

> [!warning] 论文对自己理论的保留意见
> 这一点值得如实记下，因为它容易在转述中被抹掉。上式最后一项 $4\kappa^2(R)d^2(f)$ **不随 $S\to\infty$ 消失**，因此就表面形式而言，这套理论比已有的标准 Monte Carlo 理论要**弱**。论文明确表示目前无法严格证明两者可比，只能给出两条较弱的支持：一是论证该界与 Monte Carlo 界在量级上相同（第 5.1.1 节），二是数值上表明加权最小二乘**经常**（而非总是）更优（第 6 节）。论文另外指出，一维情形下 $\lambda_{\min}(R)$ 与 $\kappa(R)$ 在一般权下数值表现都很好（图 3）。
>
> 换言之，这条路线的说服力此处更多来自实验而非定理，而下文编号 45 的诱导采样判据才把「线性区域」变成一条干净的定理。

### 数值实验

本站核实到的是实验的**设计与论文明述的定性结论**，图中的具体数值（样本数、条件数、误差量级）未逐项核实，因此下面不引用任何数字。

| 实验   | 设定                                                                                          | 论文陈述的结论                                     |
| ------ | --------------------------------------------------------------------------------------------- | -------------------------------------------------- |
| 图 1   | 一维对称 Jacobi 族（$\alpha=\beta$），$N=\dim\mathbb P_k=k+1$；分别以 $\beta$ 和以 $k$ 为横轴 | 标准 MC 的稳定性因子随次数爆炸，这就是方法的出发点 |
| §6.1.1 | 有界域的矩阵稳定性：均匀／Legendre 及其他 Jacobi 情形                                         | CLS 稳定性优于标准 MC                              |
| §6.1.2 | 无界域（高斯）的矩阵稳定性                                                                    | 支持论文对加权平衡测度所作的猜测形式（表 2）       |
| §6.1.3 | 非全次数的 $\ell_p$ 多项式空间                                                                | CLS 仍表现良好，尽管理论是按全次数空间叙述的       |
| §6.2.1 | 代数函数                                                                                      | CLS 与标准 MC 的精度对比                           |
| §6.2.2 | 一维空间中的非均匀扩散方程                                                                    | 同上                                               |
| §6.2.3 | 电阻网络                                                                                      | 同上                                               |

总体结论按摘要的措辞是：CLS 在许多值得关心的情形下优于标准 Monte Carlo；论文在别处补上了「许多（但不是全部）」。

**这些实验建立了什么。** 第一，图 1 证实稳定性因子的爆炸不是纸面上的担忧，而是常用权族的真实行为；第二，§6.1.3 表明方法的适用范围比它的定理更宽（理论按全次数空间写，实验在 $\ell_p$ 空间也成立）；第三，§6.1.2 是无界域那几条猜测形式**唯一**的支撑。

**它们没有建立什么。** 定理 5.2 里那个不随 $S$ 消失的 $4\kappa^2(R)d^2(f)$ 项，实验并没有补上：数值上「经常更优」不等于理论上可比，论文自己也没有这样宣称。此外无界域部分是在猜测的采样密度上做的实验，因此它验证的是「若猜测正确则方法有效」，而不是猜测本身。

### 与其他论文的关系

这一篇是整个采样设计纲领的理论重心。它把[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 14]] 一类工作中编目的势论词汇（平衡测度、加权平衡测度、Fekete 阵列）变成一个算法，并把 Cohen–Davenport–Leviatan 界中的 $\|K\|_\infty/N$ 换成 $1/\lambda_{\min}(R)$。它的权 $N/K(z)$ 与随机化求积路线中隐含出现的 Gauss 权是同一个对象。把**倒数** Christoffel 函数直接用作采样密度、而不是用作平衡测度样本上的权，就是编号 45 的诱导采样改良，它消掉了定理 4.3 留下的「渐近于 $k$」的缺口。用同一个 Christoffel 权、但把 Monte Carlo 换成列主元 QR 的确定性对应物，是编号 28。

## 45：把这条线索写成综述

### 直觉

综述的论点只有一句：**最「显然」的采样密度经常是糟糕的选择。** 它的例 5.2 就是为此设的——一个 $d=2$ 的高斯问题，即便取 $M=10N$ 个样本，从 $w$ 采样得到的估计仍然比最佳逼近差上几个数量级。

修法在上文已经说过：把采样密度偏置成 $\rho=q^2w$，再用 $1/q^2$ 把偏置除回去。**关键在于这一对操作不改变任何 $M\to\infty$ 的极限**——Gram 矩阵仍趋于 $I$，右端仍趋于 $\hat f_n$，目标泛函仍趋于 $\|g-f\|^2_{L^2_w}$——所以 $q$ 是一个完全自由的设计变量，可以纯粹为了「让前渐近的 Gram 矩阵尽快集中」而选。而这个选择问题有精确解。

### 问题设定

$D\subset\mathbb R^d$，$w:D\to[0,\infty)$ 是**概率密度**，$L^2_w$ 配内积 $\langle u,v\rangle=\int_Duvw\,\mathrm dx$。$V\subset L^2_w$，$\dim V=N$，$L^2_w$-正交归一基 $v_1,\dots,v_N$。最佳逼近

$$
f_N(x)=\sum_{n=1}^N\hat f_nv_n(x),\quad \hat f_n=\langle f,v_n\rangle,
\qquad f_N=\arg\min_{v\in V}\|f-v\|.
$$

由于 $g_N-f_N\in V$ 而 $f-f_N\perp V$，勾股关系给出一个干净的质量指标

$$
\eta_N:=\frac{\|g_N-f_N\|}{\|f-f_N\|},
\qquad \|f-g_N\|^2=(1+\eta_N^2)\|f-f_N\|^2 .
$$

$\eta_N\approx1$ 就表示算出来的逼近与最佳逼近一样好。**综述里所有实验报告的都是这一个数。**

多指标空间：$V(\Lambda):=\mathrm{span}\{x^\lambda:\lambda\in\Lambda\}$，$\ell^p$ 球 $\Lambda_p(k)$，双曲十字 $\Lambda^{HC}(k):=\{\lambda:\|\log(\lambda+1)\|_1\le\log(k+1)\}$，以及 $\Lambda^{TD}(k)=\Lambda_1(k)$（全次数）、$\Lambda^{ED}(k)=\Lambda_2(k)$（Euclid 次数）、$\Lambda^{TP}(k)=\Lambda_\infty(k)$（张量积）。四者的维数有固定的序：

| 指标集      | 记号              | 维数渐近                                                     |
| ----------- | ----------------- | ------------------------------------------------------------ |
| 双曲十字    | $\Lambda^{HC}(k)$ | $\sim(k+1)\log(k+1)^{d-1}$                                   |
| 全次数      | $\Lambda^{TD}(k)$ | $(k+1)^d/\Gamma(d+1)$                                        |
| Euclid 次数 | $\Lambda^{ED}(k)$ | $\bigl[\tfrac{\sqrt\pi}{2}(k+1)\bigr]^d/\Gamma(\tfrac d2+1)$ |
| 张量积      | $\Lambda^{TP}(k)$ | $(k+1)^d$                                                    |

设计矩阵：无偏置时 $(A)_{m,n}=\frac{1}{\sqrt M}v_n(x_m)$、$(f)_m=\frac{1}{\sqrt M}f(x_m)$，法方程 $Gc=g$，$G=A^TA$，$g=A^Tf$。论文提醒实际求解要用 QR 而不要显式组装法方程，法方程只用于分析。

> [!note] 一处印刷不一致
> 例 2.2 把高斯权印作 $w(x)=(2\pi)^{-d}\exp(-\|x\|_2^2)$，这与后面例 8.1 实际使用的 $w=\exp(-\|x\|_2^2)/\pi^{d/2}$ 不一致。数值部分用的是后者，本页也按后者叙述并在此标出。

### 推导

**偏置采样。** 取任意正的 $q\in L^2_w$ 且 $\|q\|=1$，于是 $\rho(x):=q^2(x)w(x)$ 仍是 $D$ 上的概率密度。从 $\rho$ 独立同分布抽 $x_1,\dots,x_M$，令

$$
(A)_{m,n}=\frac{1}{\sqrt{Mq^2(x_m)}}v_n(x_m),
\qquad
(f)_m=\frac{1}{\sqrt{Mq^2(x_m)}}f(x_m),
$$

等价于 $g_N=\arg\min_{g\in V}\frac1M\sum_m\frac{(g(x_m)-f(x_m))^2}{q(x_m)^2}$。取 $q\equiv1$ 就退回无偏置的情形。前面强调的不变性正是在这里：从 $\rho=q^2w$ 采样**并且**按 $1/q^2$ 加权，$(G)_{m,n}\to\delta_{m,n}$、$(g)_n\to\hat f_n$、目标 $\to\|g-f\|^2_{L^2_w}$ 全都不变。

**稳定性的证明是一个矩阵 Chernoff 论证**，值得完整复述，因为常数就是从这里来的。写 $G=\sum_mV_m$，$V_m=v_mv_m^T$ 为秩一半正定矩阵，$v_m^T=\frac{1}{\sqrt Mq(x_m)}(v_1(x_m),\dots,v_N(x_m))^T$。则 $\mathbb EV_m=\frac1MI$，故 $\mathbb EG=I$，$\tau_{\min}=\tau_{\max}=1$。单项的范数界是

$$
\lambda_{\max}(V_m)=\|v_m\|^2
=\frac1M\sum_{n=1}^N\Bigl(\frac{v_n(x_m)}{q(x_m)}\Bigr)^2
\le\frac1M\sup_{x\in D}\sum_{n=1}^N\Bigl(\frac{v_n(x)}{q(x)}\Bigr)^2=:Q .
$$

两侧的 Chernoff 界分别是 $\Pr[\mathcal E_{\min}]\le N(2/e)^{\tau_{\min}/2Q}$ 与 $\Pr[\mathcal E_{\max}]\le N(8e/27)^{\tau_{\max}/2Q}$，合起来

$$
\Pr[\mathcal E]\le 2N\exp\Bigl(-\frac{1}{2Q}\log\frac{27}{8e}\Bigr).
$$

要它不超过 $2M^{-r}$，只需 $\frac{\log(27/8e)}{2Q}\ge(r+1)\log M$，整理即得下面的判据，且 $\Pr[\mathcal E]\le2NM^{-r}M^{-1}\le2M^{-r}$。**常数 $C=2/\log(27/8e)$ 就是这样出现的**，它不是拟合出来的。

**最优性论证。** 判据里唯一依赖 $w$、$D$、$d$ 的量是 $\sup_x\sum_n(v_n/q)^2$，而它永远不可能小于 $N$：

$$
\sup_{x\in D}\sum_{n=1}^N\Bigl(\frac{v_n(x)}{q(x)}\Bigr)^2
\ \ge\ \int_D\sum_{n=1}^N\Bigl(\frac{v_n(x)}{q(x)}\Bigr)^2\rho(x)\,\mathrm dx
\ =\ \sum_{n=1}^N\int_Dv_n^2(x)w(x)\,\mathrm dx\ =\ N .
$$

第一步只是「上确界不小于平均」，第二步把 $\rho=q^2w$ 代入后 $q^2$ 恰好约掉——**这个约分就是整条线索的技术核心**。而

$$
q^2(x)=\frac1N\sum_{n=1}^Nv_n^2(x)
$$

**恰好取到这个下界**（归于 Cohen–Migliorati）。它就是**诱导分布**

$$
\rho(x)=\frac1N\sum_{n=1}^Nv_n^2(x)\,w(x),
$$

即 $q^2=K_\Lambda/N=1/(N\lambda_\Lambda)$：倒数 Christoffel 函数本身。它与基无关，是子空间 $V$ 的性质。「诱导」一词借自 Gautschi 与 Li 关于由给定正交多项式诱导出的一族正交多项式的工作。

代回判据就得到 $M/\log M\ge C(r+1)N$。**至此，问题的全部依赖被从样本复杂度搬进了采样密度**：复杂度只剩 $N$，而所有困难集中到「怎么从 $\rho$ 抽样」这一件事上。

### 定理

**定理 6.1（Cohen–Davenport–Leviatan 定理 1 的特例；法方程稳定性）。** 按上式定义 $A$ 与 $G$。若对某个 $r>0$，

$$
\frac{M}{\log M}\ \ge\ C(r+1)\,
\sup_{x\in D}\sum_{n=1}^{N}\Bigl(\frac{v_n(x)}{q(x)}\Bigr)^{2},
\qquad
C=\frac{2}{\log(27/8e)}\approx9.24,
$$

则以至少 $1-2M^{-r}$ 的概率有 $\|G-I\|_2\le\frac12$。

**样本复杂度判据（诱导采样）。** 代入最优的 $q$ 后，

$$
\frac{M}{\log M}\ \ge\ C(r+1)\,N,
$$

即 $M\sim N\log N$，在 $N$ 上除对数因子外已是最优。**这条判据只依赖 $N=\dim V$**：不依赖维数 $d$、不依赖域 $D$、不依赖权 $w$，甚至不依赖具体取了哪个 $N$ 维子空间。代价很明确：必须从非标准密度 $\rho$ 抽样，而 $\rho$ **确实**依赖 $(V,w,D)$。

**引理 5.1（渐近相合）。** 对任意可容许的 $q$，$M\to\infty$ 时 $g_N\to f_N$ 在 $L^2_w$ 中几乎必然成立（强大数律，加上 $G^{-1}\to I$、$g\to\hat f$ 与连续映射定理）。**所以 $q$ 的选择只在前渐近区间起作用**，这一点值得反复强调。

**定理 7.2（取自 Cohen–Migliorati；带截断的精度）。** 设 $f\in L^2_w$ 有界，$\sup_D|f|=L<\infty$，$T_L(y)=\mathrm{sign}(y)\min\{|y|,L\}$。设 $x_m$ 从诱导密度 $\rho$ 独立同分布，$M$ 满足上面的判据，则

$$
\mathbb E\,\|f-T_L\circ g_N\|^2\ \le\ \Bigl(1+\frac{4}{C(1+r)\log M}\Bigr)\|f-f_N\|^2\ +\ 8L^2M^{-r} .
$$

也就是说，在期望意义下，截断后的最小二乘逼近所犯的误差与最好可能的误差相当。论文指出对**未截断**的 $g_N$ 也有相应的高概率型结论。

**样本复杂度对照。** 综述最有信息量的部分是把三条路线并置：

- **从正交性（均匀）测度独立同分布采样**（Cohen-Davenport-Leviatan）：对任意 $r>0$，若 $M/\log M\ge C_rN^2$，则 $\Pr[|\!|\!|\hat A-I|\!|\!|\ge\frac12]\le2M^{-r}$。这是**二次**要求 $M\gtrsim N^2\log N$，且在任意**下**指标集与多维情形下同样成立。
- **从 Chebyshev 测度 Monte Carlo 采样**：要求降到 $M\sim N^{\log3/\log2}$，严格优于 $N^2$ 但仍超线性。
- **确定性 Weil 点集**（来自编号 9）：$M\ge C(d)N^2$ 给出唯一解与近最佳逼近；这个二次要求比 Chebyshev Monte Carlo 更强，其补偿是确定性。

论文同时明确指出实践与理论之间的落差：实践者通常取 $M\simeq cN$，$c$ 在 $2$ 到 $3$ 之间，即**线性**，而线性区域的确定性理论「尚未确定地可得」。**这一段是整条线索的坐标系**：编号 22 与 45 指出的线性区域，正是编号 28 用完全不同的手段去够的目标。

**渐近诱导测度（第 8 节）。** 一维、$w$ 在 $[-1,1]$ 上均匀、$V$ 为次数 $\le N-1$ 的多项式时，$\lim_{N\to\infty}\rho(x)=\rho_\infty(x)=\frac{1}{\pi\sqrt{1-x^2}}$（弱收敛），即反正弦（Chebyshev）密度——这解释了为什么 Chebyshev 网格是区间上多项式逼近的经典答案，也正是编号 22 所识别的大 $N$ 最优策略。多维、$D=[-1,1]^d$、$w$ 均匀、$V=V(\Lambda^{TD}(k))$ 时，极限是张量化的 Chebyshev 密度 $\frac{1}{\pi^d\prod_j\sqrt{1-(x^{(j)})^2}}$。

> [!warning] 高斯情形只是猜测
> 对 $D=\mathbb R^d$ 上的高斯权，编号 22 **猜测**
> $$\lim_{k\to\infty}\rho\bigl(x/\sqrt k\bigr)=C\bigl(2-\|x\|_2^2\bigr)^{d/2},$$
> $C$ 为归一化常数，并注意到输入要按 $1/\sqrt k$ 缩放。**这是猜测，不是定理**，本页凡出现此式处都按猜测对待。综述本身也说得很直白：多维情形「更多的东西是未知的，很多时候我们只有猜测」，而 $\rho^{TP}$、$\rho^{ED}$、$\rho^{TD}$、$\rho^{HC}$ 四者的大 $k$ 渐近「基本上没有被研究过」。

**从 $\rho$ 抽样是便宜的（第 8.1 节）。** 一般的多维密度抽样代价高昂，但诱导密度是「张量积密度的加性混合」，因此可以用**关于维数 $d$ 线性**的复杂度抽样。实现细节见 Narayan 关于诱导正交多项式分布计算的文章，软件在 `https://github.com/akilnarayan/induced-distributions`。

### 数值实验

综述的三组实验都有完整设定，是本页数值证据最扎实的一部分。

**例 5.2——标准采样失效的样板。**

| 项目     | 设定                                                                                                            |
| -------- | --------------------------------------------------------------------------------------------------------------- |
| 域与权   | $D=\mathbb R^2$，$w(x)=\exp(-\lVert x\rVert^2)/\pi$                                                             |
| 空间     | $\Lambda=\Lambda^{TD}(k)$，$N=\binom{k+2}{k}=(k+1)(k+2)/2$                                                      |
| 次数     | $k=1,\dots,25$                                                                                                  |
| 测试函数 | $f(x)=B\bigl(\lVert x/4-(0.2,-0.1)\rVert_2\bigr)$，$B$ 为一维 bump 函数                                         |
| 样本数   | $M=10N$                                                                                                         |
| 重复     | 100 次                                                                                                          |
| 结果     | 从 $w$ 采样得到的 $g_N$「极不准确」，$\eta_N$ 极大；诱导采样在 $k=20$（$N=231$）时 $\eta_N\sim1$，且 $M/N$ 适中 |

**例 8.1——维数的影响。**

| 项目     | 设定                                                                                                                             |
| -------- | -------------------------------------------------------------------------------------------------------------------------------- |
| 域与权   | $D=\mathbb R^d$，$w=\exp(-\lVert x\rVert_2^2)/\pi^{d/2}$                                                                         |
| 测试函数 | $f(x)=\prod_{j=1}^d\exp\bigl([x^{(j)}]^2/j\bigr)$，取乘积形式是为了高维下 $f_N$ 仍可计算                                         |
| 空间     | $\Lambda=\Lambda^{HC}(k)$                                                                                                        |
| $(d,k)$  | $(4,20)$、$(8,10)$、$(20,5)$                                                                                                     |
| 重复     | 100 次                                                                                                                           |
| 结果     | 诱导采样一致优于标准采样，**但优势随 $d$ 增大而减小**；诱导采样的 $\eta_N$ 的定性行为随 $d$ 基本不变，与定理 6.1、7.2 的预期一致 |

优势变小的原因论文说得很清楚：$d$ 大时只能取小 $k$，而低次空间使 $\rho$ 与 $w$ 本来就接近。

**第 9 节——参数化热扩散方程。**

| 项目     | 设定                                                                        |
| -------- | --------------------------------------------------------------------------- |
| 方程     | 单位正方形上 $-\nabla\cdot(a(y,x)\nabla u)=S$，$u\rvert_{\partial\Omega}=0$ |
| 源项     | $S(y,x)=100\chi_F(y)$，$F$ 为边长 $0.2$ 的居中方形子域                      |
| 几何     | 四个半径 $r=0.13$ 的圆形夹杂，关于中心对称                                  |
| 系数     | $a(y,x)=1+\sum_{i=1}^4x^{(i)}\chi^i(y)$，故 $d=4$                           |
| 参数分布 | $x^{(i)}\sim U(-0.99,0.2)$ 独立同分布，对应张量 Legendre 基                 |
| 感兴趣量 | $f(x)=u\bigl((0.25,0.375),x\bigr)$                                          |
| 离散     | 空间用 $P_1$ 有限元；$V=V(\Lambda^{HC}(25))$                                |
| 误差估计 | $Q=10\,000$ 个参考样本，100 次重复                                          |
| 结果     | 诱导分布加权最小二乘「也许略好一点」，但「差别并不显著」                    |

另有两组图：图 4 在 $N-1=19$ 处对三个一维权（$[-1,1]$ 上均匀、$[0,\infty)$ 上指数、$\mathbb R$ 上高斯）比较 $w$、$\rho$ 与 $\rho_\infty$，显示 $\rho$ 与 $\rho_\infty$ 视觉上接近，而 $w$ 与 $\rho$ 可以差别很大——$\rho$ 要求在 $w$ 几乎不要样本的区域大量取样；图 5 画出 $d=2$ 高斯、$k=8$ 时 TD／ED／TP／HC 四种诱导密度的等值线，显示它们彼此之间、以及与 $w$ 之间都相差明显。

**这些实验建立了什么。** 例 5.2 把「标准采样有时会坏得离谱」从警告变成了可复现的事实；例 8.1 验证了定理最要紧的那条性质——诱导采样的表现**不随维数退化**；图 4、图 5 则说明「换密度」不是微调，而是量级上的改变。

**它们没有建立什么。** 第 9 节的真实 PDE 例子里，诱导采样只是「略好」。论文没有回避这一点，本页也照录：**能站住的结论是诱导采样始终位列最好的一档、并且是唯一有最小样本数定理支撑的那一个，而不是它总能在数值上取胜。** 此外定理 7.2 保证的是**截断后**估计在**期望**意义下近最佳，而实验报告的是 $\eta_N$ 的分布，两者不能直接互推。

### 与其他论文的关系

编号 45 是整条线索的教学总结。它把编号 22 呈现为**渐近**版的最优采样结果（平衡测度采样加 Christoffel 加权，只在次数趋于无穷时最优），而把诱导采样（Cohen–Migliorati）呈现为**非渐近的精确**版本——这正是它的判据没有渐近附加条件、而编号 22 的定理 4.3 有的原因。它的 $q^2=\frac1N\sum v_n^2$ 就是编号 24 记号里的 $1/(N\lambda_\Lambda)$、编号 28 记号里的 $K_\Lambda/N$：**三篇论文操作的是同一个函数，区别只在于它是被用作平衡测度样本上的权（22、24、28），还是被用作采样密度本身（45）。** 它的 Chebyshev 极限 $\rho_\infty$ 恰是编号 24 的 CSA-a 所用的 Chebyshev 预条件；它所引的高斯猜测与编号 28 §4.2 的候选密度是同一个函数形式。

## 28：用贪心选点代替随机采样

### 直觉

随机设计只能给出「以高概率成立」的保证；而经典的**近似 Fekete 点**（贪心地最大化 $|\det V|$）只对**紧** $\Gamma$ 有定义，因此对高斯／Hermite 问题**根本不存在**。此外可以优化的目标其实有四个——全局行列式、全局条件数，以及两者各自的贪心版本——先验地没有理由认为它们一致。

**把 Christoffel 权塞进行列式目标里，一举解决了这两件事。** 机制可以说得很具体：加权空间 $Q$ 的 Vandermonde 矩阵，其**每一行的范数恒等于 1**，因为第 $j$ 行的平方范数是 $\sum_\alpha\psi_\alpha^2(y_j)/K_\Lambda(y_j)=1$。对行范数为 1 的矩阵，Hadamard 不等式给出 $|\det V|\le1$，且等号成立当且仅当各行两两正交，也就是 $\kappa(V)=1$。**于是「行列式最大」与「条件数最小」在最优点上必然重合**，这是加权空间独有的性质，未加权情形没有对应物。（这一段是本页对定理 3.1 的直观解读，论文自身的证明本站未逐行核对。）注意它与编号 22 的自我评述完全对上：Christoffel 加权就是把设计矩阵的行归一化。

至于无界域：$1/\sqrt{K_\Lambda}$ 在远处压住了多项式的增长，因此加权表述对非紧 $\Gamma$ 有定义，而未加权的 AFP 没有。

### 问题设定

$\Gamma=\prod_{i=1}^d\Gamma_i\subset\mathbb R^d$ 为张量域，$\rho(y)=\prod_i\rho_i(y_i)$ 为张量积概率密度（各分量独立），$\phi^i_n$ 为 $\rho_i$ 的 $n$ 次正交归一多项式，多元基 $\psi_\alpha(y)=\prod_j\phi^j_{\alpha_j}(y_j)$ 满足 $\langle\psi_\alpha,\psi_\beta\rangle=\delta_{\alpha,\beta}$。指标集取全次数 $\Lambda^{TD}_k=\{\alpha:|\alpha|\le k\}$ 与双曲十字 $\Lambda^{HC}_k=\{\alpha:\prod_j(\alpha_j+1)\le k+1\}$；$P=\mathrm{span}\{\psi_\alpha:\alpha\in\Lambda\}$，$N=|\Lambda|$。

倒数 Christoffel 函数与加权空间：

$$
K_\Lambda(y)=\sum_{\alpha\in\Lambda}\psi_\alpha^2(y),
\qquad
Q=\mathrm{span}\Bigl\{\tfrac{\psi_\alpha}{\sqrt{K_\Lambda}}\ \Big|\ \alpha\in\Lambda\Bigr\}.
$$

类 Vandermonde 矩阵 $(V(A_m,P))_{j,k}=\psi_{\alpha(k)}(y_j)$，加权版本 $(V(A_m,Q))_{j,k}=\psi_{\alpha(k)}(y_j)/\sqrt{K_\Lambda(y_j)}$。所解的加权最小二乘问题是

$$
\hat f=\arg\min_{v\in\mathbb R^N}\bigl\|V(A_M,Q)\,v-Wf\bigr\|_2,
\qquad (W)_{m,m}=1/\sqrt{K_\Lambda(y_m)} .
$$

矩形矩阵的行列式模按 $|\det V|=\sqrt{|\det(VV^T)|}$ 定义（$1\le m\le N$，$m=N$ 时与通常的行列式模一致）。未加权的 Fekete 点与条件数最优点集分别是

$$
A^F_N(P):=\arg\max_{A_N\in\Gamma^N}|\det V(A_N,P)|,
\qquad
A^C_N(P):=\arg\min_{A_N\in\Gamma^N}\kappa\bigl(V(A_N,P)\bigr),
$$

一般来说这是两个不同的集合。Fekete 点给出经典的 Lebesgue 常数界 $\|I_N\|_{C(\Gamma)\to C(\Gamma)}\le N$（论文注明实践中观察到的是对数增长）。把上面两式中的 $P$ 全部换成 $Q$，就得到加权的全局目标；相应的贪心迭代是

$$
y^{F*}_{n+1}=\arg\max_{y\in\Gamma}\bigl|\det V\bigl(A^{F*}_n\cup y,\ Q\bigr)\bigr|,
\qquad
y^{C*}_{n+1}=\arg\min_{y\in\Gamma}\kappa\bigl(V(A^{C*}_n\cup y,\ Q)\bigr),
$$

前者就是 **CFP 方法**（Christoffel 加权近似 Fekete 点），简写为 $y_{n+1}=\arg\max_y\det|WVV^TW|$。起点 $y^{F*}_1$ 在 $\mathbb R$ 中可任取，但**确实影响结果**；出现平局时会有多个「分支」。

### 推导

一维情形有完整的结构，这也是全篇最漂亮的部分。取 $\Gamma=\mathbb R$ 上任意概率密度 $\rho$、$\Lambda=\{0,\dots,N-1\}$，记 $\phi_N$ 为 $N$ 次正交归一多项式，并定义亚纯函数

$$
r_N(y)=\frac{\phi_N(y)}{\phi_{N-1}(y)} .
$$

则对任意 $y\notin\phi_{N-1}^{-1}(0)$，集合 $A_N(y)=r_N^{-1}\bigl(r_N(y)\bigr)$（取集值逆）**唯一确定**，且它承载一个 $N$ 点正权求积公式，对 $2N-2$ 次以下精确，其权恰是 **Christoffel 权**：

$$
\int_\Gamma p(z)\rho(z)\,\mathrm dz=\sum_{z\in A_N(y)}\frac{1}{K_\Lambda(z)}\,p(z),
\qquad \deg p\le2N-2 .
$$

特别地，$y\in\phi_N^{-1}(0)$ 时 $A_N(y)$ 就是 $N$ 点 Gauss 公式。**这里出现的 $1/K_\Lambda(z)$ 与编号 22 的权 $N/K(z)$ 是同一个对象**（相差归一化因子 $N$），两条看起来不同的路线在这一点上合流；这也是 Gauss 求积权本来就是倒数 Christoffel 权这一事实的来源。

由此，一维情形的证明思路是：先证明满足 $|\det V(A_N,Q)|=1$ 的构型存在且由 $r_N$ 的等值集刻画，再证明贪心迭代从这样的构型中任一点出发都会重建整个构型。

论文实现上还有三处不能略过的工程细节：

- **过采样靠扩充指标集（§4.1）。** CFP 本身只产生 $N=\dim Q$ 个点，而最小二乘要 $M>N$。做法是：取下封闭的 $\Lambda$，令 $n=\max\{|\alpha|:\alpha\in\Lambda\}$，$S:=\Lambda^{TD}_n\setminus\Lambda$（若 $\Lambda=\Lambda^{TD}_n$ 则先令 $n\leftarrow n+1$），把 $S$ 按总次数排序、同次数以逆字典序断平，取前 $\Delta N$ 个并入得到大小 $M=N+\Delta N$ 的 $\tilde\Lambda$，定义 $\tilde Q:=Q_{\tilde\Lambda}$。
- **连续极大化换成候选池（§4.2）。** 对 $\Gamma$ 的极大化换成在有限候选集 $\tilde A$ 上的极大化。
- **贪心一步就是列主元 QR（§4.3）。** 组装 $V(\tilde A,\tilde Q)$ 后贪心地挑出张成体积最大的那些行，「这件事只需对 $V^T$ 作列主元 QR 分解即可完成」，主元的顺序就给出 $A_M$。注意这个转置：在 $V^T$ 的**列**上选主元，就是在候选 Vandermonde 矩阵的**行**上做选择，也就是在选点。
- **最终求解（§4.4）。** 点集 $A_M$ 定下来之后，在**原来**的大小为 $N$ 的指标集 $\Lambda$ 上求解加权问题，而不是在扩充后的指标集上。

### 定理

**定理 3.1（行列式最优 $\Leftrightarrow$ 条件数最优，且贪心达到全局最优）。** 设 $\rho:\Gamma\to[0,\infty)$ 是 $\mathbb R^d$ 上的概率密度，$\Lambda$ 是大小为 $N$ 的任意多指标集并由此定义 $Q$。构型 $A_N$ 满足 $|\det V(A_N,Q)|=1$ **当且仅当** $\kappa(V(A_N,Q))=1$；因此两个全局问题取到最优值的解相互重合。进而，若 $A_N$ 达到其中之一（从而两者皆达到），则：$y^{F*}_1\in A_N$ 时行列式贪心迭代有一个分支给出 $A^{F*}_N=A_N$；$y^{C*}_1\in A_N$ 时条件数贪心迭代有一个分支给出 $A^{C*}_N=A_N$。

> [!warning] 定理 3.1 的前提
> 该定理**预设了单位条件数构型的存在性**。作者明确指出，$d>1$ 时这一存在性难以验证（他们引用了若干已知的多维非平凡例子）。也就是说，多维情形下的最优性结论是**有条件的**。

**引理 3.1（一维结构）。** 上文推导中的五条结论：$y\notin\phi_{N-1}^{-1}(0)$ 时存在满足最优条件的 $A_N=A_N(y)$；$A_N(y)$ 作为 $y$ 的函数唯一；$A_N(y)=r_N^{-1}(r_N(y))$；$A_N(y)$ 承载 $2N-2$ 次精确的 $N$ 点正权求积公式，权为 Christoffel 权；$y\in\phi_N^{-1}(0)$ 时 $A_N(y)$ 是 Gauss 公式。论文提醒即便 $\rho$ 的支集是紧的，也必须把它当作整个 $\mathbb R$ 上的密度来看待。

**定理 3.2（一维：贪心＝全局＝最优，对几乎所有起点成立）。** 取 $\Gamma=\mathbb R$ 上任意密度 $\rho$、$\Lambda=\{0,\dots,N-1\}$，固定 $y\notin\phi_{N-1}^{-1}(0)$，由引理 3.1 取 $A_N(y)$，令 $A^F_N(Q)=A^C_N(Q)=A_N(y)$ 并以 $y^{F*}_1=y^{C*}_1=y$ 初始化。则四个点集完全重合，且

$$
\bigl|\det V\bigr|=1=\kappa(V),
$$

即所得 Vandermonde 矩阵是**完美条件**的。也就是说一维情形下贪心不是近似，而恰好取到最优。更要紧的是这条结论对**任意**一维密度成立，**包括非紧支撑的密度**——而这正是近似 Fekete 点所不具备的性质，因此它也是对无界域困难的一个正面回答。

> [!note] 一处印刷错误
> 该定理正文把起点条件印作 $y\in\phi_{N-1}^{-1}(0)$，这与它所依赖的引理 3.1 及推论 3.1 相矛盾（那里要求 $y$ **不在**零集中），应为 $y\notin\phi_{N-1}^{-1}(0)$。本页按后者叙述并在此标出。

**推论 3.1。** 对 $y\notin\phi_{N-1}^{-1}(0)$，从 $y^{F*}_1=y$ 出发的贪心迭代给出唯一的、具有最优多项式精度的正权 $L^2_\rho$ 求积节点；若 $y\in\phi_N^{-1}(0)$，给出的就是 $\rho$-加权的 Gauss 节点。

**没有被证明的部分（论文自述的限制）。** 这里**没有样本复杂度定理**，即没有形如「$M\gtrsim N\log N$」的界。第 4 节的过采样办法被作者明确称为「基本上是临时性的」，并未宣称其最优性。$d>1$ 时的最优性结论都以定理 3.1 的存在性假设为前提。

### 数值实验

全部实验共用同一套协议：

| 项目   | 取值                                                                                                   |
| ------ | ------------------------------------------------------------------------------------------------------ |
| 候选池 | $\tilde M=10^4$，一半从 $\rho$ 独立同分布抽，一半从受编号 22 启发的次数渐近密度抽                      |
| 过采样 | $\Delta N=\lfloor0.05N\rfloor$，即 $M=1.05N$——**线性、5% 过采样**                                      |
| 重复   | 每种配置 50 次；报告条件数均值以及 20%、80% 分位                                                       |
| 对照   | MC（从 $\rho$ 独立同分布、未加权 $V(A_M,P)$）、Fekete（AFP、未加权）、C-Fekete（CFP、加权 $V(A_M,Q)$） |

候选池的第二个系综：$[-1,1]^d$ 上均匀 $\rho$ 时取张量积 Chebyshev；$\mathbb R^d$ 上高斯 $\rho$ 时取支于半径 $\sqrt{2n}$ 的球、密度为

$$
C_d\Bigl(1-\tfrac{1}{2n}\sum_{k=1}^ds_k^2\Bigr)^{d/2},\qquad s\in\mathbb R^d,
$$

$C_d$ 为归一化常数。**这与编号 45 所引的高斯渐近诱导测度是同一个函数形式，而那个形式只是猜测**，因此这里的候选密度也应看作由猜测启发的启发式选择。作者说明本可以用弱可容许网格作候选，但已知构造在维数上增长太快，故未采用。

| 实验    | 设定                                                                                         | 结果                                               |
| ------- | -------------------------------------------------------------------------------------------- | -------------------------------------------------- |
| §5.1.1  | 矩阵稳定性，有界域：$[-1,1]^d$ 上 Legendre，$d=2,6,10$，TD 与 HC 两种指标集                  | CFP 给出的系统明显比 AFP 与 MC 更稳定              |
| §5.1.2  | 矩阵稳定性，无界域：高斯 $\rho\propto\exp(-\lVert y\rVert_2^2)$，张量 Hermite，$d=2,6,10,25$ | CFP 在每一例都更稳定，但作者自述**高维下改善有限** |
| §5.2(a) | 精度：$f(y)=\exp\bigl(-\sum_{j=1}^dy_j^2\bigr)$，Legendre 逼近，$d=2$                        | CFP 与 AFP 相当且不更差，远好于 MC                 |
| §5.2(b) | 精度：一维空间的随机椭圆方程，$d=2,6$（TD）与 $d=10,25$（HC）                                | 同上                                               |

§5.2(b) 的方程是 $(0,1)\times\mathbb R^d$ 上的

$$
-\frac{\mathrm d}{\mathrm dx}\Bigl[\kappa(x,y)\frac{\mathrm du}{\mathrm dx}(x,y)\Bigr]=f,
\qquad u(0,y)=u(1,y)=0,\quad f=2,
$$

扩散系数取 KL 型展开

$$
\kappa(x,y)=1+\sigma\sum_{k=1}^d\frac{1}{k^2\pi^2}\cos(2\pi kx)\,y_k,
$$

感兴趣量为 $u(0.5,y)$，$\rho$ 在 $[-1,1]^d$ 上均匀。（$\sigma$ 的取值本站未核实，故不引用。）

**这些实验建立了什么。** 第一，在仅 5% 过采样、即 $M=1.05N$ 这个极其接近 $N$ 的样本预算下，CFP 仍给出可用的条件数——这正是编号 45 所说「实践者取 $M\simeq cN$」而理论够不到的那个区域；第二，无界域一栏本身就是结论，因为未加权的 AFP 在那里**根本没有定义**；第三，两组精度实验表明加权没有以精度为代价换稳定性。

**它们没有建立什么。** 报告的量是**条件数**，而不是当初激励 Fekete 点的 Lebesgue 常数；实验没有、也无法弥补缺失的样本复杂度定理——$M=1.05N$ 是一个被观察到可行的选择，不是被证明可行的选择；高维下相对 AFP 的改善作者自己都说有限；而 $d>1$ 的最优性理论仍系于定理 3.1 那个难以验证的存在性假设。

### 与其他论文的关系

CFP 是编号 22 那套随机 Christoffel 加权最小二乘的确定性、贪心对应物——同一个权 $1/K_\Lambda$，甚至候选网格就是从编号 22 的平衡测度型密度抽出来的，区别只在于选点用列主元 QR 而不是 Monte Carlo。编号 45 与诱导采样一路是靠**精确诱导采样**取得非渐近最优性，编号 28 则是靠**确定性优化**取得（有条件的）精确最优性；论文的历史回顾明确把这两者摆成对编号 22 遗留的「渐近于次数」缺口的两种竞争性解答。

**编号 22、45 与 28 因此构成一个完整的谱**：

| 路线                 | 样本要求             | 结论类型    | 需要的额外结构     |
| -------------------- | -------------------- | ----------- | ------------------ |
| 从正交性测度采样     | $M\gtrsim N^2\log N$ | 概率        | 无                 |
| Christoffel 加权采样 | 对数线性             | 概率        | 平衡测度 + $K(z)$  |
| 贪心（近似 Fekete）  | 接近 $N$             | 确定性/势论 | 候选池 + 列主元 QR |

## 24：把采样与预条件作为一对来设计

### 直觉

编号 24 是编号 22 的 $\ell_1$ 版本。稀疏恢复的样本复杂度由**有界正交系**的一致上界 $L$ 决定；而高次或无界域上的多项式基**不是**一致有界的，所以从 $w$ 做 Monte Carlo 采样在次数高时恢复得很差。

已有的两种补救各有短板：Rauhut–Ward 的「渐近采样」（从 Chebyshev 抽样、对 Vandermonde 矩阵作预条件）只适用于**有界**随机变量，且精度随参数维数增长而退化；Hampton–Doostan 的相干最优采样虽然通用，却需要 MCMC 才能从它的测度抽样。论文要的是一个显式、对有界域上任意可容许权与无界域上一大类指数权都适用、并且带有证明过的样本数判据的统一方案。

**核心认识是那个互为倒数的结构**：采样密度取 $\propto w/(N\lambda)$，即 $w$ 乘以 Christoffel 函数的**倒数**；预条件权取 Christoffel 函数**本身**。采样与加权互逆，复合之后的系统才一致有界。因此「采样密度」与「预条件矩阵」应当作为一对来设计，而不是先定采样再补预条件。

### 问题设定

Christoffel 函数与预条件（论文式 (9)）：

$$
(W)_{m,m}=N\lambda_\Lambda(Z^{(m)}),
\qquad
\lambda_\Lambda(Z)=\frac{1}{\sum_{i\in\Lambda}\varphi_i^2(Z)} .
$$

一元情形 $\Lambda=\{0,1,\dots,n\}$、$N=n+1$ 时写作 $\lambda_{n+1}(z)=1/\sum_{k=0}^n\varphi_k^2(z)$。

> [!note] 一处印刷不一致
> 论文用 $\lambda$ 表示 Christoffel 函数本身、$N\lambda$ 表示其缩放版本，但算法 1 第 4 行印的是 $(W)_{m,m}=N/\lambda_\Lambda(Z^{(m)})$，与式 (9) 的 $N\lambda_\Lambda$ 相矛盾，两者必有其一是笔误。本页按原文照录并在此标出，不作静默更正。

**算法 1（Christoffel 稀疏逼近，CSA）** 五步：(1) 从平衡密度 $v=\frac{\mathrm d\mu}{\mathrm dZ}$ 独立同分布抽 $M$ 个样本 $\{Z^{(m)}\}$；(2) 组装 $f$，$f_m=f(Z^{(m)})$；(3) 组装 $M\times N(\Lambda)$ 的类 Vandermonde 矩阵 $\Phi$，$\Phi_{m,i}=\varphi_i(Z^{(m)})$；(4) 由（缩放的）Christoffel 函数取值算出对角预条件；(5) 解预条件后的基追踪去噪问题 $\alpha^\star=\arg\min_\alpha\|\alpha\|_1$ 受约束于 $\|\sqrt W\Phi\alpha-\sqrt Wf\|_2\le\varepsilon$。

**三种具体情形**：

- **CSA-a**：$Z$ 在 $[-1,1]$ 上服从 Beta 分布，形状参数 $\beta+1,\alpha+1\ge\frac12$（即 Jacobi 参数 $\alpha,\beta\ge-\frac12$）——从 **Chebyshev 密度**采样，与 $n$ 无关，支集 $S_n\equiv[-1,1]$。
- **CSA-b**：$Z$ 在 $\mathbb R$ 上服从双边指数型分布，$w(z)=\exp(-|z|^\alpha)$，$\alpha>1$——从一个**扩张的**平衡测度采样，密度 $v_n$ 与支集 $S_n$ 都依赖 $n$。
- **CSA-c**：$Z$ 在 $[0,\infty)$ 上服从单边指数型分布，$w(z)=\exp(-|z|^\alpha)$，$\alpha>\frac12$——类似的 $n$-依赖扩张平衡测度。

多维无界构造的抽样过程以一次缩放 $Z=4W_n$ 收尾，生成的样本落在 $\mathbb R^d$ 中 $\ell_1$ 范数不超过 $4n$ 的点集上。

Christoffel 加权基的 Gram 矩阵（式 (21)）：

$$
R_{k,\ell}=\int_{S_n}\varphi_{k-1}(z)\varphi_{\ell-1}(z)\,\bigl(N\lambda_n(z)\bigr)\,v_n(z)\,\mathrm dz,
\qquad 1\le k,\ell\le N=n+1 .
$$

$R$ 正定，且每个固定位置的元素随 $n\to\infty$ 收敛到单位阵的对应元素；$R^{1/2}$ 是它唯一的对称正定平方根；$\|A\|_1$ 指诱导的 $\ell_1$ 矩阵范数（最大列 $\ell_1$ 范数）。

### 推导

设计原则可以写成一句等式：采样密度须满足 $v\approx w/(N\lambda_N)$，其中 $N\to\infty$ 时的等号由加权多势论平衡测度给出——也就是说，Christoffel 函数相对 $w$ 的渐近行为由平衡测度支配。这与编号 22 的推导是同一条渐近恒等式，只是用途从「加权」换成了「配对」。

复合之后为什么一致有界，一行就能看出来：加权基 $\sqrt{N\lambda_\Lambda}\,\varphi_i$ 满足

$$
\sum_{i\in\Lambda}\bigl(\sqrt{N\lambda_\Lambda}\,\varphi_i\bigr)^2
=N\lambda_\Lambda\sum_{i\in\Lambda}\varphi_i^2=N
$$

**恒等地成立**——与编号 22 中 $\tilde K_k\equiv N$ 是同一个约分。有了一致有界性，有界正交系的限制等距理论就可以接手，样本要求形如 $M\gtrsim L\,s\log^3(s)\log(N)$，其中 $L$ 是该系统的一致上界。论文的三种情形之间的全部差别，就落在这个 $L$ 上。

代价则在 $R$：加权后的基并不在 $N\lambda_n v_n$ 下正交归一，必须再经 $R$ 变换才是。这就是为什么下面的定理恢复的是 $R^{1/2}\alpha$ 而不是 $\alpha$，也是为什么误差界的两项分别以 $\lambda_{\min}(R)$ 与 $\|R^{-1/2}\|_1$ 为系数——**与编号 22 一样，$R$ 的谱是最终出现在界里的量。**

### 定理

**定理 4.1（主恢复定理，一元）。** 从平衡密度 $v_n$ 独立同分布抽 $M$ 个点，按 $\Phi_{ij}=\varphi_j(Z^{(i)})$ 组装 $\Phi$、按式 (9) 组装对角阵 $W$。设

$$
M\ \ge\ L(n)\,\bigl\|R^{-1/2}\bigr\|_1^2\;s\log^3(s)\log(N),
$$

其中 $L(n)$ 按三种情形分别为：

| 情形  | $L(n)$                                          | 分支展开                                                                               |
| ----- | ----------------------------------------------- | -------------------------------------------------------------------------------------- |
| CSA-a | $C(\alpha,\beta)$，关于 $n\ge1$ 一致            | 与 $n$ 无关                                                                            |
| CSA-b | $Cn^{\max\{1/\alpha,\,2/3\}}$，$C=C(\alpha)$    | $\alpha\ge\tfrac32$ 时为 $Cn^{2/3}$；$1<\alpha<\tfrac32$ 时为 $Cn^{1/\alpha}$          |
| CSA-c | $Cn^{\max\{1/(2\alpha),\,2/3\}}$，$C=C(\alpha)$ | $\alpha\ge\tfrac34$ 时为 $Cn^{2/3}$；$\tfrac12<\alpha<\tfrac34$ 时印作 $Cn^{1/\alpha}$ |

则以超过 $1-N^{-\gamma\log^3(s)}$ 的概率，对一切多项式 $p(x)=\sum_j\alpha_j\varphi_j(x)$，在带噪样本 $f=\Phi\alpha+\eta$、$\|W\eta\|_\infty\le\varepsilon$ 下，解

$$
R^{1/2}\alpha^\star=\arg\min_\alpha\|R^{1/2}\alpha\|_1
\quad\text{s.t.}\quad
\|\sqrt W\Phi\alpha-\sqrt Wf\|_2\le\varepsilon
$$

给出

$$
\|\alpha-\alpha^\star\|_2\le C_1\frac{\sigma_s\bigl(R^{1/2}\alpha\bigr)_1}{\sqrt s\,\lambda_{\min}(R)}+C_2\frac{\varepsilon}{\sqrt{\lambda_{\min}(R)}},
$$

$$
\|\alpha-\alpha^\star\|_1\le D_1\sigma_s\bigl(R^{1/2}\alpha\bigr)_1\bigl\|R^{-1/2}\bigr\|_1+D_2\sqrt s\,\bigl\|R^{-1/2}\bigr\|_1\varepsilon .
$$

> [!note] 一处内部不一致
> CSA-c 的 $\max$ 表达式写作 $\max\{1/(2\alpha),2/3\}$，而紧随其后分支展开时印的是 $Cn^{1/\alpha}$（而非 $Cn^{1/(2\alpha)}$）。本页按原文照录并在此标出，不作静默更正。

**有界与无界的差别全部落在 $L(n)$ 上。** 有界（Chebyshev 型）情形 $L(n)$ 关于 $n$ 一致有界，因此要求就是 $M\gtrsim s\log^3(s)\log(N)$——在稀疏度 $s$ 上近乎线性，即**最优**；无界情形要求退化为 $M\gtrsim s\,n^{2/3}$，多出一个随多项式次数增长的因子。论文的注 4.2 指出，在它所用的有界正交系分析框架内，$\alpha\ge\frac32$ 时的 $n^{2/3}$ 与 $\alpha<\frac32$ 时的 $n^{1/\alpha}$ 都是**紧的**，CSA-c 有类似陈述。

**推论 4.1（干净的 Legendre 情形）。** CSA-a 中取 $\alpha=\beta=0$（$Z$ 均匀、Legendre 基）时，对一切 $n$ 有 $R=I$，于是在 $M\ge Cs\log^3(s)\log(N)$ 下，**不作变换的**问题 $\alpha^\star=\arg\min\|\alpha\|_1$ s.t. $\|\sqrt W\Phi\alpha-\sqrt Wf\|_2\le\varepsilon$ 就满足 $\|\alpha-\alpha^\star\|_2\le C_1\sigma_s(\alpha)_1/\sqrt s+C_2\varepsilon$ 与 $\|\alpha-\alpha^\star\|_1\le D_1\sigma_s(\alpha)_1+D_2\sqrt s\varepsilon$。

> [!warning] 定理与所实现算法之间的三条落差，论文自述
> 其一，**理论只有一元**。注 4.3 说结论可以推广到张量积域与张量积权，只要按各一元密度作张量式采样；此时判据形式不变，但 $L(n)$ 变成 $d$ 个一元因子的**乘积**——维数依赖是指数的。注 4.1 预期 Beta／Jacobi 情形可以推广到紧区间上几乎任意有界权。
> 其二，定理恢复的是 $R^{1/2}\alpha$ 而**不是** $\alpha$，因为 PCE 基只有在经 $R$ 变换之后才在 $N\lambda_N(z)v_n(z)$ 下正交归一；而**实际的算法 1 直接恢复 $\alpha$**，所以定理并没有字面上分析所实现的那个算法。
> 其三，无界情形的采样密度依赖 $n$，实现上需要那族扩张平衡测度，这与编号 22 在无界域遇到的是同一类困难。

论文明确把框架的适用范围写成「有界或无界」域，因此它同时回应了[[computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares|编号 11]] 在最小二乘情形遇到的无界域困难。

### 数值实验

本站核实到的是四组实验的**主题与论文明述的定性结论**；图中给出的是各参数随 $N$ 的变化曲线，具体的维数、稀疏度、样本数与误差量级未逐项核实，故不引用数字。

| 实验 | 设定                                       | 论文陈述的结论                      |
| ---- | ------------------------------------------ | ----------------------------------- |
| §6.1 | Beta 随机变量上的渐近（Chebyshev）采样     | CSA-a 情形的验证                    |
| §6.2 | 高斯随机变量，Chebyshev 方法在此**不适用** | 这是无界情形只能靠 CSA-b 的直接证据 |
| §6.3 | 人为构造的稀疏解                           | 恢复质量的受控检验                  |
| §6.4 | 带随机输入的椭圆型偏微分方程               | 应用层面的检验                      |

摘要中给出的总体结论有两条：CSA 在许多值得关心的情形下优于标准 Monte Carlo；并且**即便与专门针对 Legendre 或 Hermite 基设计的算法相比，精度也相当甚至更好**——对一个通用方案而言这是个值得注意的说法。论文另外指出，虽然理论是一元的，数值上在高维多元设定中表现良好。

**这些实验建立了什么。** 最关键的一条是：**多元情形的全部证据都在实验里**，因为定理只覆盖一元。§6.2 则说明这套方案在 Chebyshev 预条件根本不适用的地方仍然工作，这正是它相对 Rauhut–Ward 的卖点。

**它们没有建立什么。** 实验不能弥补注 4.3 指出的指数维数依赖：多元推广的常数是 $d$ 个一元因子的乘积，而实验只能表明在所试的若干维数上表现良好。定理与算法之间那条「恢复 $R^{1/2}\alpha$ 还是 $\alpha$」的落差，实验同样没有触及。

### 与其他论文的关系

编号 24 之于 $\ell_1$，正如编号 22 之于最小二乘——同一批作者、同一套 Christoffel 权与平衡测度的配对，论文自己也明说 CSA「基于我们在 [22] 中引入的一个用于离散最小二乘逼近的类似算法」。它关于 Christoffel 加权多项式的上界被其他工作用于量化「子采样张量 Gauss 求积网格」情形的样本数；它那个足以同时覆盖有界与无界域的预条件框架，后来被推广到包含梯度观测的情形；而它的 CSA-a（Chebyshev 采样）情形与 Rauhut–Ward 的预条件是一回事。它显式写出的那个互逆结构——从 $\propto w/(N\lambda)$ 采样、用 $N\lambda$ 加权——正是编号 45 的诱导分布所**精确**实现（而非渐近实现）的同一件事。

## 54：把同一思路搬到再生核空间

### 直觉

前面四篇都在多项式空间里工作，但回头看整页开头的那个念头，多项式其实从未真正参与：起作用的只是「一个有限维函数空间在各点有多集中」。凡是有再生核的地方，这个量就存在，同样的设计原则就该适用。**编号 54 说明的正是这条线索的抽象层次**：多项式空间被换成再生核 Hilbert 空间后，Christoffel 函数的角色由核的对角线承担，而「把点放到空间集中的地方去」这一设计原则不变。

这句话要落到实处，得说清楚是**哪一条**对角线，否则类比会走样。对平移不变的径向核，$\mathbf A_{ii}=\Phi(0)$ 对一切 $i$ 都相同，所以**核本身的对角线是常数**，没有任何东西可以拉平。真正与多项式情形的 $K_\Lambda=\sum_\alpha\psi_\alpha^2$ 对应的，是**已选点所张成的那个子空间**的再生核对角线，即 $k(z)^{T}\mathbf A_N^{-1}k(z)$（$k(z)$ 为 $z$ 与已选点之间的核向量）；它的补

$$
P_N^2(z)=\Phi(0)-k(z)^{T}\mathbf A_N^{-1}k(z)
$$

就是**幂函数**的平方，也就是以已有数据为条件的高斯过程后验方差。论文选下一个点的方式正是把它取到最大。**于是「把点放到空间集中的地方」在核的语言里读作「把下一个点放到当前张成空间够不到的地方」——同一条规则的正反两面。** 论文自己也是这样定位的：多项式文献里的贪心叫近似 Fekete 点，「对核插值而言，它是一个幂函数极大化方法」。

**为什么点一多就不稳。** 论文给的机制很具体（注 2）：两个中心靠得太近时，以它们为心、形状相同的基函数在全部节点上取值几乎相同，插值矩阵就出现两列近乎相同的值。而在有界的参数域上，$N$ 增大意味着填充距离必然缩小，成簇是躲不掉的——所以病态不是坏运气，是点数增长的必然后果。论文观察到，Sobol 点、低差异点、随机点这些常用选择，在 $N\sim\mathcal O(10^2)$ 时都出现这个问题。形状参数 $\epsilon$ 也帮不上忙，它两头受限：$\epsilon$ 太大插值不准，$\epsilon$ 太小矩阵病态到线性求解器失效，而怎么选最好的 $\epsilon$ 论文明说「仍是一个开问题」。梯度增强的情形更糟：把导数观测接进来之后，系统矩阵 $\mathbf B$ 的条件数比普通插值矩阵 $\mathbf A$ 增长得**快得多**。

**为什么从大候选池里贪心选点能对付它。** 论文把无从下手的问题换成了一个可解的问题：不去优化 $\epsilon$、也不去在连续域上找最优中心，而是问「给定一大堆候选点，怎样挑出一个『最优』子集」。这样问的好处有两条。其一，挑选只用到线性代数，本质上就是设计矩阵的 Cholesky 型分解；其二——这是实用上更要紧的一条——挑选**与数据无关**：$\mathbf A$ 与 $\mathbf B$ 只依赖点的位置，不依赖 $u$ 的取值，所以整个设计可以在任何一次昂贵的 PDE 求解之前离线完成。

### 问题设定

参数化问题写成一般形式：$Z=(Z^{(1)},\dots,Z^{(d)})\in I_Z\subseteq\mathbb R^d$ 为不确定输入，

$$
\begin{cases}
u_t(x,t,Z)=\mathcal L(u), & \mathcal D\times(0,T]\times I_Z,\\
\mathcal B(u)=0, & \partial\mathcal D\times(0,T]\times I_Z,\\
u=u_0, & \mathcal D\times\{t=0\}\times I_Z,
\end{cases}
$$

$\mathcal D\subset\mathbb R^l$（$l=1,2,3$）为物理域。固定 $x$ 与 $t$ 后记 $u_j=u(x,t,z_j)$，目标是由 $(z_j,u_j)$、$j=1,\dots,N$ 构造 $u_N(Z)\approx u(Z)$。

**核插值。** 取平移不变的径向核 $K=\Phi(\epsilon\lVert\cdot-\cdot\rVert)$，$\Phi$ 为径向基函数，$\epsilon$ 为形状参数。以 $\Xi=\{z_j\}_{j=1}^N$ 为中心，试探空间为

$$
\mathcal U_\Xi=\mathrm{span}\bigl\{\Phi(\epsilon\lVert\cdot-z_j\rVert)\ \big|\ z_j\in\Xi\bigr\},
\qquad
u_N(Z)=\sum_{j=1}^Nc_j\,\Phi(\epsilon\lVert Z-z_j\rVert).
$$

论文只考虑**正定**核，正定性保证逼近的存在唯一性；所用的三种是高斯 $\Phi(r)=\exp(-r^2)$、逆多重二次（IMQ）$\Phi(r)=1/\sqrt{1+r^2}$ 与紧支撑（CS）族。插值条件 $u_N(z_i)=u_i$ 给出

$$
\mathbf{Ac}=\mathbf u,\qquad \mathbf A=K(\Xi,\Xi),\quad \mathbf A_{ij}=\Phi(\epsilon\lVert z_i-z_j\rVert),
$$

$\mathbf A$ 对称，$\mathbf A$ 可逆时 $\mathbf c$ 唯一。注 1 记下这套过程与高斯过程的关系：以 $K$ 为协方差、由 $\Xi$ 上的实现 $\mathbf u$ 求最优无偏估计，得到的就是同一个函数（附录 A 给出简单克里金的推导，结论是 $u_N(Z)=K(Z,\Xi)^TK(\Xi,\Xi)^{-1}\mathbf u=\hat f(Z)$，两个估计**恒等**）。注 2 就是上文那条病态机制。

**LOOCV 选形状参数。** 论文不打算系统分析 $\epsilon$ 的选取，而是用留一交叉验证：误差向量的分量为

$$
e_i(\epsilon)=\frac{c_i}{\mathbf A^{-1}_{ii}},
\qquad
\lVert e(\epsilon^*)\rVert=\min_\epsilon\lVert e(\epsilon)\rVert,
$$

其中 $\mathbf A^{-1}_{ii}$ 是 $\mathbf A^{-1}$ 的第 $i$ 个对角元。这个形式的用处在于：留一误差可以从**已经算出的** $\mathbf c$ 与 $\mathbf A^{-1}$ 的对角线直接读出，不必真的解 $N$ 个去掉一点的系统。

**梯度增强表述。** 若除 $\{z_i,u(z_i)\}_{i=1}^N$ 外还有 $\{z_i,u'_m(z_i)\}_{i=1}^N$（$m=1,\dots,d$，$u'_m=\partial u/\partial z^{(m)}$），则插值函数取

$$
u_N(Z)=\sum_{j=1}^Nc_j\Phi(\epsilon\lVert Z-z_j\rVert)-\sum_{m=1}^d\sum_{j=1}^N\beta_{m,j}\Phi'_m(\epsilon\lVert Z-z_j\rVert),
$$

并要求满足 $N(d+1)$ 个广义插值条件 $\lambda_iu=\lambda_iu_N$，其中 $\lambda_i$ 或为在 $z_i$ 处取值、或为在 $z_i$ 处取某个导数。系统矩阵是分块的

$$
\mathbf B=\begin{pmatrix}
\mathbf A_{0,0}&\mathbf A_{0,1}&\cdots&\mathbf A_{0,d}\\
\mathbf A_{1,0}&\mathbf A_{1,1}&\cdots&\mathbf A_{1,d}\\
\vdots&\vdots&\ddots&\vdots\\
\mathbf A_{d,0}&\mathbf A_{d,1}&\cdots&\mathbf A_{d,d}
\end{pmatrix},
\qquad
(\mathbf A_{m,n})_{i,j}=\begin{cases}
\Phi(\epsilon\lVert z_i-z_j\rVert), & m=n=0,\\
-\Phi'_n(\epsilon\lVert z_i-z_j\rVert), & m=0,\ n\neq0,\\
\Phi'_m(\epsilon\lVert z_i-z_j\rVert), & m\neq0,\ n=0,\\
-\Phi''_{m,n}(\epsilon\lVert z_i-z_j\rVert), & m\neq0,\ n\neq0,
\end{cases}
$$

$\mathbf B$ 是对称的 $(d+1)N\times(d+1)N$ 矩阵。附录 B 指出：当 RBF 核与协方差核取同一个时，$\mathbf B$ 与梯度增强高斯过程的联合协方差矩阵**恒等**。这条恒等正是下一节贪心准则的意义所在——被极大化的量因此有一个概率解释。

### 推导

**从 Fekete 点出发。** 把插值写成基数形式

$$
u_N(Z)=\sum_{j=1}^Nc_j\Phi(\epsilon\lVert Z-z_j\rVert)=\sum_{j=1}^Nd_j\,\ell_j(Z),
\qquad
\ell_j(Z)=\frac{\det\mathbf A(z_1,\dots,z_{j-1},Z,z_{j+1},\dots,z_N)}{\det\mathbf A(z_1,\dots,z_N)},
$$

Lebesgue 常数为 $\Lambda_N:=\max_{Z\in I_Z}\sum_{j=1}^N\lvert\ell_j(Z)\rvert$。**Fekete 点就是使 Vandermonde 型行列式最大的构型**

$$
\Xi^*=\arg\max_{\Xi=\{z_1,\dots,z_N\}\subset I_Z}\bigl\lvert\det\mathbf A(z_1,\dots,z_N)\bigr\rvert .
$$

由基数公式，$\Xi^*$ 上任一 $\ell_j$ 的分子都是把某个 $z_j$ 换成 $Z$ 之后的行列式，其绝对值不超过分母，故 $\lvert\ell_j(Z)\rvert\le1$，于是 $\Lambda_N\le N$。论文补充说实践中观察到的 Lebesgue 常数增长常常是次线性的。困难在于：除若干特殊的多项式情形外，这些点没有已知的显式刻画，而直接求解上面的优化问题「是一件令人生畏的事」。标准的松弛办法是贪心：一次加一个点，每次让行列式最大。

> [!note] §3.1 贪心式的一处印刷错误
> 论文把贪心迭代印作 $z_{N+1}=\arg\max_{z\in I_Z}\lvert\det\mathbf A(z_1,\dots,z_N)\rvert$，而被极大化的表达式里**根本没有出现 $z$**，因此按字面读它是一个与自变量无关的常数。应为 $\det\mathbf A(z_1,\dots,z_N,z)$，即把新点并入之后的行列式——§3.2 对 $\mathbf B$ 写的对应式子正是这个形式。本页按后者叙述并在此标出。（另有一处较轻的：图 4 的图题印作「条件数关于形状参数关于采样点数 $N$」，而该图的横轴是 $N$、$\epsilon=3,5$ 是固定参数。）

论文明说，多项式文献里把这套贪心称作**近似 Fekete 点**，而**对核插值而言它是一个幂函数极大化方法**——这就是上文那条等价关系的出处。

**贪心一步到底在最大化什么。** §3.2 把贪心写在梯度增强矩阵上：

$$
z_{N+1}:=\arg\max_{z\in I_Z}\det\boldsymbol B(z_1,\dots,z_N,z).
$$

直接按定义算这个行列式要重新组装整个 $\boldsymbol B$，代价太高。论文的做法是用 Schur 补把它拆开。取置换矩阵 $\boldsymbol P_{N+1}$，把新点带来的 $d+1$ 个元素搬到末尾（其余 $N(d+1)$ 个保持原序），则

$$
\boldsymbol P_{N+1}\boldsymbol B\boldsymbol P_{N+1}^{-T}
=\begin{pmatrix}\boldsymbol B_N&\boldsymbol W(z)\\\boldsymbol W^T(z)&\boldsymbol B(z)\end{pmatrix},
\qquad
\boldsymbol B_N:=\boldsymbol B(z_1,\dots,z_N),
$$

于是

$$
\det\boldsymbol B(z_1,\dots,z_N,z)
=\det\boldsymbol B_N\cdot\det\bigl(\boldsymbol B(z)-\boldsymbol W^T\boldsymbol B_N^{-1}\boldsymbol W\bigr).
$$

**因子 $\det\boldsymbol B_N$ 与 $z$ 无关，可以整个丢掉**，贪心因此等价于

$$
z_{N+1}=\arg\max_{z\in I_Z}F(z),
\qquad
F(z):=\det\bigl(\boldsymbol B(z)-\boldsymbol W^T\boldsymbol B_N^{-1}\boldsymbol W\bigr).
$$

这就回答了「每个主元最大化什么」：**是 Schur 补的行列式，即以已经放好的全部函数值与梯度观测为条件、候选点 $z$ 上那 $(d+1)\times(d+1)$ 条件协方差块的行列式。** 由附录 B 的恒等关系，这个块正是梯度增强高斯过程在 $z$ 处的后验协方差，所以贪心的意思是「下一个点放在当前模型最不确定的地方」。把它退化到无梯度的标量情形，$\boldsymbol B(z)$ 变成 $\Phi(0)$、$\boldsymbol W$ 变成 $k(z)$，Schur 补就是 $\Phi(0)-k(z)^T\mathbf A_N^{-1}k(z)=P_N^2(z)$，即上文那个幂函数的平方。

**为什么这是带主元的 Cholesky 分解。** 评估 $F$ 的主要代价是解 $\boldsymbol B_N\boldsymbol X=\boldsymbol W$（$\boldsymbol X\in\mathbb R^{N(d+1)\times(d+1)}$），其余不过是一个 $(d+1)\times(d+1)$ 行列式，可以忽略。论文的办法是保留置换后 $\boldsymbol B_N$ 的 Cholesky 因子

$$
\boldsymbol P_N\boldsymbol B_N\boldsymbol P_N^T=\boldsymbol L_N\boldsymbol L_N^T,
$$

存下 $\boldsymbol L_N^{-1}$，则每次评估 $F$ 只需作用一次 $\boldsymbol L_N^{-1}$，代价 $\mathcal O(d^3N^2)$，优化写成

$$
z_{N+1}=\arg\max_{z\in I_Z}\det\bigl(\boldsymbol B(z)-\boldsymbol V^T\boldsymbol V\bigr),
\qquad
\boldsymbol V:=\boldsymbol L_N^{-1}\boldsymbol W(z).
$$

选定 $z_{N+1}$ 之后，因子按分块更新，无需重新分解：

$$
\boldsymbol L_{N+1}=\begin{pmatrix}\boldsymbol L_N&\boldsymbol 0\\\boldsymbol V^T&\widetilde{\boldsymbol L}\end{pmatrix},
\qquad
\widetilde{\boldsymbol L}\widetilde{\boldsymbol L}^T=\boldsymbol B(z_{N+1})-\boldsymbol V^T\boldsymbol V,
$$

$$
\boldsymbol L_{N+1}^{-1}=\begin{pmatrix}\boldsymbol L_N^{-1}&\boldsymbol 0\\-\widetilde{\boldsymbol L}^{-1}\boldsymbol V^T\boldsymbol L_N^{-1}&\widetilde{\boldsymbol L}^{-1}\end{pmatrix}.
$$

一步的全部动作因此是三件事：(i) 由上式算 $\boldsymbol V(z_{N+1})$；(ii) 求 Schur 补的 $(d+1)\times(d+1)$ Cholesky 因子 $\widetilde{\boldsymbol L}$ 之逆；(iii) 拼出 $\boldsymbol L_{N+1}^{-1}$，再作用一次 $\boldsymbol L_N^{-1}$。**这就是带主元的 Cholesky 分解本身**：$\boldsymbol L_N$ 是已选点核矩阵的 Cholesky 因子，每一步把剩余 Schur 补里「最大」的那一块提上来作为下一个主元；论文的总结一节也正是这样描述自己的方法的（「对 Vandermonde 型插值矩阵使用高效的带主元 Cholesky 分解来挑选样本点」）。最后一步是实现上的让步：**连续域 $I_Z$ 上的极大化换成有限候选集上的极大化**，论文明说实际计算的从来是后者。

### 定理

**这篇论文不含任何定理。** 全文没有 theorem、lemma、proposition、corollary 或 proof 环境，一条也没有。文中出现的一般性论断只有两类：其一是引自 Fekete 点经典结论的 $\lvert\ell_j\rvert\le1\Rightarrow\Lambda_N\le N$，其二是附录 A、B 中核插值与（梯度增强）高斯过程估计的恒等关系；两者都以陈述形式给出、归于文献，本文不作证明。

这一点必须说清楚，因为标题里的「optimal design」与正文里的「quasi-optimal」都容易被读成有定量含义的词。**它们没有。** 论文没有给出贪心解与真正极大化解之间的任何比较，没有条件数的界，没有所选点的 Lebesgue 常数的界，也没有收敛率。「拟最优」在这里只是「对一个无从直接求解的组合优化所作的贪心松弛」的名字。全部定量证据都在图里。

> [!note] 与编号 28 的对照就在这里最尖锐
> 编号 28 在多项式空间里对结构完全相同的贪心证明了两条实质结论：行列式最优与条件数最优在加权空间中重合（定理 3.1，附带一个在 $d>1$ 难以验证的存在性前提），以及一维情形下贪心恰好取到最优（定理 3.2）。编号 54 把同一策略搬到核空间，**但没有带上任何一条定理**。核空间里对应的结论是否成立，本文既未证明也未讨论。

### 数值实验

**共用协议。** 全部实验用同一套设置。

| 项目   | 取值                                                                                       |
| ------ | ------------------------------------------------------------------------------------------ |
| 候选池 | $M=10^4$ 个在 $I_Z$ 上随机取的中心；另比较「uniform」（$I_Z$ 上等距点）等候选              |
| 选点   | 由 §3.2 的贪心从候选池中选出 $N$ 个作为 RBF 中心（文中记作「Cholesky」）                    |
| 对照   | random（随机点）、sobol（Sobol 序列）、halton（Halton 序列）；PDE 一例另加稀疏网格         |
| 误差   | $E_{\ell_2}=\bigl(\tfrac1Q\sum_{i=1}^Q\lvert s(z_i)-u(z_i)\rvert^2\bigr)^{1/2}$，$Q=1000$ 个随机测试点 |
| 重复   | 每种配置 50 次；报告均值与 20%、80% 分位                                                   |
| 核     | 高斯 $\Phi(r)=e^{-r^2}$、IMQ $\Phi(r)=1/\sqrt{1+r^2}$、Wendland 紧支撑（取 $l=\lfloor d/2\rfloor+4$） |

> [!warning] 本节不能给出条件数与误差的具体数值
> 该文的**全部数值结果只以图给出，没有一张数据表**：图 1–图 14 都是曲线，正文里没有任何条件数或误差以数字形式出现。（该文的 HTML 版本里能数出四十余个表格元素，但它们无一例外是公式的排版容器，不是数据表。）因此本节报告的是可以逐项核实的**实验设定与论文明述的结论**，不引用任何条件数量级或误差数值。

**条件数实验。** 前两张图在方法提出之前就摆出了问题本身，后三张检验方法。

| 图   | 对象                                              | 横轴          | 固定参数                    |
| ---- | ------------------------------------------------- | ------------- | --------------------------- |
| 图 1 | $\mathbf A$（Lagrange，实线）与 $\mathbf B$（Hermite，虚线） | 形状参数 $\epsilon$ | $N=50,100,300$；$d=2$ |
| 图 2 | 同上                                              | 点数 $N$      | $\epsilon=0.1,1,3$；$d=2$   |
| 图 3 | $\mathbf A$，四种选点法对照                        | $\epsilon$    | $N=100,300$；$d=2$          |
| 图 4 | $\mathbf A$，四种选点法对照                        | $N$           | $\epsilon=3,5$；$d=2$       |
| 图 5 | $\mathbf A$，不同候选池                            | $N$           | $\epsilon=3,5$              |

论文对图 1、图 2 的结论是：形状参数趋近于零或点数变大时，设计矩阵都趋于奇异，而 $\mathbf B$ 的恶化远快于 $\mathbf A$。对图 3、图 4 的结论是：**高斯与 IMQ 两种核下，本方法都比其他采样方法稳定得多**——原文这句话只点了这两种核，紧支撑核虽然画在同一批图里，却没有被这句结论覆盖。图 5 的结论是候选点的选取（随机或等距）不显著影响本方法的表现。

**精度实验。** 第一个测试函数是 $[0,1]^2$ 上的 Franke 基准函数

$$
\begin{aligned}
u(z)=&\tfrac34e^{-((9z^{(1)}-2)^2+(9z^{(2)}-2)^2)/4}+\tfrac34e^{-(9z^{(1)}+1)^2/49-(9z^{(2)}+1)^2/10}\\
&+\tfrac12e^{-((9z^{(1)}-7)^2+(9z^{(2)}-3)^2)/4}-\tfrac15e^{-(9z^{(1)}-4)^2-(9z^{(2)}-7)^2}.
\end{aligned}
$$

> [!note] 第二项的指数印刷有误
> 原文该项印作 $-(9z^{(2)}+1)^2)/10$，括号不配对；且标准的 Franke 函数在这一项上是**一次**的 $-(9z^{(2)}+1)/10$，不是平方。上式按原文照录（并补全括号）并在此标出差异。

第二个是一维空间中的随机椭圆方程，与编号 28 的 §5.2(b) 是**同一个基准问题**：

$$
-\frac{\mathrm d}{\mathrm dx}\Bigl[\kappa(x,z)\frac{\mathrm du}{\mathrm dx}(x,z)\Bigr]=f,
\quad(x,z)\in(0,1)\times\mathbb R^d,
\qquad u(0,z)=u(1,z)=0,\quad f=2,
$$

$$
\kappa(x,z)=1+\sigma\sum_{k=1}^d\frac{1}{k^2\pi^2}\cos(2\pi kx)\,z^{(k)},
$$

感兴趣量为 $u(z)=u(0.5,z)$，$z^{(i)}\sim U[-1,1]$ 独立同分布。（$\sigma$ 的取值原文未给出，故不引用。）

| 图    | 测试问题                    | 维数  | 对照                                    | 备注                    |
| ----- | --------------------------- | ----- | --------------------------------------- | ----------------------- |
| 图 6  | Franke 函数，$[0,1]^2$      | $d=2$ | random／sobol／halton                   | $\epsilon=3,5$          |
| 图 7  | 同上，换候选池              | $d=2$ | uniform 等候选                          | $\epsilon=3,5$          |
| 图 8  | 同上                        | $d=2$ | 不同 $\epsilon$，以及 LOOCV 选出的 $\epsilon$ | 左条件数、右误差   |
| 图 9  | 随机椭圆方程                | $d=3$ | **加稀疏网格**（Legendre 全次数 $k=8$） | 左高斯／中 IMQ／右 CS   |
| 图 10 | 随机椭圆方程                | $d=6$ | **加稀疏网格**（Legendre 全次数 $k=4$） | 左高斯／中 IMQ／右 CS   |

论文对图 6 的结论有两层，第二层比第一层有意思：本方法精度优于其余采样方法；而且**随 $N$ 增大呈现清晰的收敛型态，其他采样方法的误差曲线则明显表现出「点加得越多不一定越准」**。图 8 的结论是 $\epsilon$ 越小矩阵越病态，而用 LOOCV 选 $\epsilon$ 可以得到很好的结果。

**梯度增强实验。** 三个测试函数：

| 图    | 测试函数                                                                              | 维数  | 域          |
| ----- | ------------------------------------------------------------------------------------- | ----- | ----------- |
| 图 12 | corner peak：$u(z)=\bigl(1+\sum_{i=1}^d\omega_iz^{(i)}\bigr)^{-(d+1)}$，$\omega_i=1/i^2$ | $d=2$ | $[0,1]^2$   |
| 图 13 | Rastrigin：$u(z)=20+\sum_{i=1}^2\bigl((z^{(i)})^2-10\cos(2\pi z^{(i)})\bigr)$          | $d=2$ | $[-4,4]^2$  |
| 图 14 | Friedman：$u(z)=10\sin(\pi z^{(1)}z^{(2)})+20(z^{(3)}-0.5)^2+10z^{(4)}+5z^{(5)}$       | $d=5$ | 原文未标明  |

（图 11 画的是 Rastrigin 函数本身，不含结果。）论文的结论是：在本算法之下梯度增强设计矩阵 $\mathbf B$ 可以是良态的，且精度仍优于其余采样方法。

**关于稀疏网格的比较，实际发生了什么。** 摘要说本方法「在许多值得关心的情形下可以优于稀疏网格方法」。逐节核对全文之后，与稀疏网格的比较**只出现一次**：§4.1.2 的随机椭圆方程，$d=3$（Legendre 全次数 $k=8$）与 $d=6$（$k=4$），结果见图 9、图 10，每图三种核各一栏。论文对此的原话是 RBF 逼近方法「明显优于稀疏网格方法」。**全文没有任何一处报告过本方法输给稀疏网格，也没有任何一处报告过它输给 random、sobol、halton 中的任何一种。**

所以这句对冲措辞的含义与它听上去的不一样：它不是「我们试了很多情形，其中许多赢了（言下之意另一些没赢）」，而是「我们做了一次比较，赢了，并谨慎地不把它宣布成普遍结论」。**支撑它的证据是一个基准问题、两个参数维数、三种核。** 论文也没有把两者的**代价**并置：图 9、图 10 的横轴是点数 $N$，而贪心选点本身要在 $M=10^4$ 个候选上、每步每个候选付出 $\mathcal O(d^3N^2)$，这部分未计入；反过来稀疏网格的节点数也不是可以任意指定的连续量，两条曲线在横轴上并不真正对齐。

**这些实验建立了什么。** 第一，图 1、图 2 把「病态随点数增长」从担忧变成了可复现的事实，并且指明梯度增强让它严重得多——这是全篇的出发点，也解释了为什么论文要在 $\mathbf B$ 而不是 $\mathbf A$ 上做贪心。第二，图 3、图 4 表明在同样的点数下，选点方式本身就能把条件数拉开可见的差距，而这不需要任何关于 $u$ 的信息。第三，图 5、图 7 表明结果对候选池的具体取法不敏感，这对一个「先挑候选、再挑点」的两段式方法是必要的稳健性检验。第四，图 6 的那条观察——其他采样方法加点不一定更准——本身就是稳定性论证的经验形式。

**它们没有建立什么。** 最要紧的一条上面已经说过：没有定理，所以「拟最优」没有任何量化含义，全部结论都是这些图上的曲线。其次，形状参数除图 8 外一律固定在 $\epsilon\in\{3,5\}$，而论文自己说最优 $\epsilon$ 是开问题，LOOCV 只在一张图里被检验过。其三，全部测试的参数维数都在 $d\le6$，而立论动机是高维不确定性量化；候选池固定为 $M=10^4$ 且不随 $d$ 调整，$d=5,6$ 时这是相当稀的覆盖，而论文除「随机对等距」之外没有研究 $M$ 的影响。其四，文中显式出现的最大点数是 $N=300$（图 1、图 3），以 $N$ 为横轴的各图其范围在正文与图题中都未给出，因此「推迟不稳定性」这一说法能覆盖到多大的 $N$，从文本无法判断。

### 与其他论文的关系

**最该被看见的是它与编号 28 的结构同一性。** 两篇论文做的是同一件事：不在连续域上优化，而是先备一个大候选池，再用一次带主元的矩阵分解把点一个一个贪心地挑出来。

| 项目       | 编号 28                        | 编号 54                              |
| ---------- | ------------------------------ | ------------------------------------ |
| 函数空间   | 多项式空间（加权后的 $Q$）     | 再生核 Hilbert 空间                  |
| 被最大化量 | Vandermonde 型行列式           | 核矩阵行列式（Schur 补的行列式）     |
| 分解       | 对 $V^T$ 作**列主元 QR**       | 对 $\boldsymbol B$ 作**带主元 Cholesky** |
| 候选池     | $\tilde M=10^4$                | $M=10^4$                             |
| 统计       | 50 次重复，20%／80% 分位       | 50 次重复，20%／80% 分位             |
| 定理       | 定理 3.1、3.2（有条件的最优性）| 无                                   |

两处分解之所以能各自实现贪心，理由是同一个：**两种分解都把行列式写成主元的连乘**（QR 给出 $\prod\lvert r_{ii}\rvert$，Cholesky 给出对角块的连乘），所以「加一个点使行列式最大」自动等于「取下一个主元」。差别只在函数空间，以及由此决定的分解类型——非对称的设计矩阵用 QR，对称正定的核矩阵用 Cholesky。

两篇论文对「权」的处理为何不同，也顺着这条线说得通。编号 28 必须插入 Christoffel 权 $1/\sqrt{K_\Lambda}$，因为多项式基在各点的集中程度千差万别，加权之后 Vandermonde 矩阵的行范数才恒为 1，Hadamard 不等式才咬得住；而在核这一侧，平移不变性使 $\mathbf A_{ii}=\Phi(0)$ 恒为常数，**这是平移不变性白送的**，所以没有可插的权。剩下的全部空间变化转移到了条件化之后的残差对角线上，也就是幂函数——这正是编号 54 的贪心所极大化的东西。两篇论文一个在**加权的**行列式上贪心、一个在**残差的**对角线上贪心，是同一条设计原则在两种几何下的形态。

两篇甚至共用同一个 PDE 基准：一维空间的随机椭圆方程，KL 型扩散系数 $\kappa=1+\sigma\sum_k\frac{1}{k^2\pi^2}\cos(2\pi kx)z^{(k)}$，感兴趣量 $u(0.5,\cdot)$，$z$ 在 $[-1,1]^d$ 上均匀。作者群重叠，实验协议也重叠。

**梯度增强这一侧则通向另一页。** 一次伴随求解同时给出函数值与全部 $d$ 个偏导，测量矩阵的行数因此乘以 $d+1$ 而代价增加甚少——这正是[[computational-mathematics/paper-notes/stochastic-approximation/sparse-recovery-and-data-driven-pce|稀疏恢复与数据驱动混沌]]一页上编号 29 与编号 32 所利用的同一件事，只是那里额外信息喂给的是 $\ell_1$ 稀疏恢复，这里喂给的是高斯过程仿真器。**两条线索付出的代价也是同一类。** 编号 32 发现朴素堆叠梯度行会破坏均值各向同性，必须用行预条件与列归一化修回来；编号 54 遇到的是它在条件数上的版本——图 1、图 2 记录了 $\mathbf B$ 的条件数比 $\mathbf A$ 增长得快得多。两者的修法各在一端：编号 32 修的是**矩阵**（预条件），编号 54 修的是**点**（选址）。

**在本页内部，它是这条线索在函数空间层面的推广。** 编号 22、24、28 用的是多项式空间的 $K_\Lambda=\sum_\alpha\psi_\alpha^2$，编号 45 用同一个量的归一化倒数作采样密度，而核方法里的对应物是已选子空间的再生核对角线及其补（幂函数）。四者是同一个设计原则在不同空间上的实例——**唯一的差别是，只有编号 54 这一例没有配上定理。**

> [!note] 覆盖进度
> 编号 **22** 与 **45** 已按可核实的材料完整展开：构造、推导链条（包括编号 45 的矩阵 Chernoff 证明与常数 $C=2/\log(27/8e)$ 的来源）、带假设的定理与实验设定。其中编号 45 的三组实验（例 5.2、例 8.1、第 9 节热扩散 PDE）有完整的数值设定与结论；编号 **22** 的实验只有设计与定性结论，图中的具体数值本站未核实。编号 **28** 的定理 3.1、引理 3.1、定理 3.2、推论 3.1 与全部实验协议（候选池 $10^4$、$M=1.05N$、50 次重复、三方对照）已核实。编号 **24** 的主恢复定理、三种情形的 $L(n)$ 与推论 4.1 已核实，其四组实验只有主题与定性结论。编号 **54** 现已按其预印本（arXiv:2104.06291）全文核实：问题设定、Schur 补与带主元 Cholesky 的完整推导链条、以及全部实验设定（候选池 $M=10^4$、$Q=1000$ 个测试点、50 次重复、三种核、四类测试问题）均取自原文。但该文**不含任何定理**，且其全部数值结果只以图给出、没有一张数据表，因此本页不引用它的条件数与误差数值。
>
> 本页保留了六处原文层面的问题：编号 22 定理 4.3 的三条限制与主定理中不消失的 $4\kappa^2(R)d^2(f)$ 项、编号 45 例 2.2 与例 8.1 之间的高斯归一化常数不一致、编号 28 定理 3.2 的起点条件印刷错误、编号 24 算法 1 第 4 行与式 (9) 的不一致、编号 24 CSA-c 分支展开的内部不一致、以及编号 54 §3.1 贪心式中遗漏新点 $z$ 的印刷错误。高斯渐近诱导测度在本页每一处出现时都按**猜测**对待。

## 覆盖核对

| 内容                          | 论文 | 覆盖状态                                                     |
| ----------------------------- | ---- | ------------------------------------------------------------ |
| 稳定性因子与其基无关性        | 22   | $K(z)=\varphi^T\varphi$、$\lVert K\rVert_\infty/N$、不变性   |
| 解耦采样密度与正交性密度      | 22   | 洞察表述、$\tilde K_k\equiv N$ 的约分与 Christoffel 函数定义 |
| 算法五步与行归一化的等价说法  | 22   | 完整步骤与等价解读                                           |
| Christoffel 渐近与平衡测度    | 22   | 定理 4.2、推论 4.1、定理 4.3 及其三条限制                    |
| CLS 稳定性与精度定理          | 22   | 定理 5.1、5.2 与不消失的 $4\kappa^2(R)d^2(f)$ 项             |
| 实验设计与定性结论            | 22   | 图 1、§6.1.1–6.1.3、§6.2.1–6.2.3；具体数值未核实             |
| 三种采样的样本复杂度对照      | 45   | 二次、$N^{\log3/\log2}$、确定性二次                          |
| 矩阵 Chernoff 证明与常数来源  | 45   | 定理 6.1 的完整论证与 $C=2/\log(27/8e)$                      |
| 诱导分布与最优性论证          | 45   | 定义 7.1、下界 $N$ 的取等、判据 $M/\log M\ge C(r+1)N$        |
| 带截断的精度定理              | 45   | 定理 7.2 与引理 5.1                                          |
| 三组数值实验                  | 45   | 例 5.2、例 8.1、第 9 节 PDE 的完整设定与结论                 |
| 渐近诱导测度                  | 45   | 一维与多维 Chebyshev 极限；高斯情形按猜测处理                |
| 加权近似 Fekete 点与列主元 QR | 28   | 选点机制、加权含义、过采样与候选池的工程细节                 |
| 行范数归一与两个目标的等价    | 28   | 定理 3.1、存在性前提与本页给出的直观解读                     |
| 一维贪心即最优                | 28   | 引理 3.1、定理 3.2、推论 3.1 与印刷错误                      |
| 三方对照实验                  | 28   | $\tilde M=10^4$、$M=1.05N$、50 次重复、四组设定与结论        |
| 采样与预条件成对设计          | 24   | 互逆结构、算法 1 五步、三种情形与 Gram 矩阵 $R$              |
| 主恢复定理与有界／无界的差别  | 24   | 定理 4.1、$L(n)$ 三种情形、推论 4.1 与两处不一致             |
| 理论一元而实验多元            | 24   | 注 4.3 的指数维数依赖与四组实验的定性结论                    |
| 从多项式空间到再生核空间      | 54   | 幂函数取代 Christoffel 函数；Schur 补贪心与带主元 Cholesky   |
| 与编号 28 的结构同一性        | 54   | 同为候选池贪心；QR 主元对 Cholesky 主元，两种行列式连乘      |
| 四类实验与稀疏网格比较        | 54   | $M=10^4$、$Q=1000$、50 次重复；稀疏网格比较只出现一次        |
| 原文不含任何定理              | 54   | 无 theorem／lemma／proof 环境，「拟最优」不带量化含义        |

## 本页原文

- A. Narayan, J. Jakeman, and T. Zhou, [_A Christoffel function weighted least squares algorithm for collocation approximations_](https://doi.org/10.1090/mcom/3192), Math. Comput. 86 (2017), pp. 1913-1947。
- J. Jakeman, A. Narayan, and T. Zhou, [_A generalized sampling and preconditioning scheme for sparse approximation of polynomial chaos expansions_](https://doi.org/10.1137/16M1063885), SIAM J. Sci. Comput. 39(3) (2017), pp. A1114-A1144。
- L. Guo, A. Narayan, L. Yan, and T. Zhou, [_Weighted approximate Fekete points: sampling for least-squares polynomial approximation_](https://doi.org/10.1137/17M1140960), SIAM J. Sci. Comput. 40(1) (2018), pp. A366-A387。
- L. Guo, A. Narayan, and T. Zhou, [_Constructing least-squares polynomial approximations_](https://doi.org/10.1137/18M1234151), SIAM Rev. 62(2) (2020), pp. 483-508。
- A. Narayan, L. Yan, and T. Zhou, [_Optimal design for kernel interpolation: applications to uncertainty quantification_](https://doi.org/10.1016/j.jcp.2020.110094), J. Comput. Phys. 430 (2021), 110094。
