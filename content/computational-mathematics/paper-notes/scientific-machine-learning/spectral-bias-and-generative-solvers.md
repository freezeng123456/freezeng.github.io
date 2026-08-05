---
title: 频谱偏差与生成式求解
description: 编号 81、89、94、101、103、105：频率内容应当被测量，而不是被假设
lang: zh
translation: en/computational-mathematics/paper-notes/scientific-machine-learning/spectral-bias-and-generative-solvers
tags:
  - 论文笔记
  - 科学机器学习
  - 频谱偏差
---

> [!note] 本页覆盖
> 编号 **81**（_Comput. Methods Appl. Mech. Engrg._ 437, 2025）、**89**（_Neural Networks_ 194, 2026，[arXiv:2401.02080](https://arxiv.org/abs/2401.02080)）、**94**（_J. Comput. Phys._ 558, 2026）、**101**（[arXiv:2512.18586](https://arxiv.org/abs/2512.18586)）、**103**（[arXiv:2606.22514](https://arxiv.org/abs/2606.22514)）、**105**（[arXiv:2604.07169](https://arxiv.org/abs/2604.07169)）。

## 频谱偏差是一个可以被测量的量

### 两个独立的发现，同一个现象

频谱偏差（spectral bias）与频率原则（F-Principle）由两组人各自发现、各自命名，两篇论文互相承认对方。

Rahaman 等人的说法是：深度网络对**低频函数**有学习偏好，即「全局变化而无局部起伏的函数」，表现为一个**依赖频率的学习速度**。他们的机制并不是关于优化器的陈述，而是关于**网络自身 Fourier 谱的衰减速率**。其定理 1 说：对 ReLU 网络 $f_\theta$（取 $f(\bm x)=(2\pi)^{d/2}\int\tilde f(\bm k)e^{\mathrm i\bm k\cdot\bm x}\mathrm d\bm k$ 的 Fourier 约定），

$$
\tilde f_{\theta}(\bm k)=\sum_{n=0}^{d}
\frac{C_{n}(\theta,\bm k)\,\mathbf 1_{H^{\theta}_{n}}(\bm k)}{k^{n+1}},
\qquad k=\|\bm k\| ,
$$

其中 $H^{\theta}_{n}$ 是与网络某个线性区域 $P_\epsilon$ 的某个 $n$ 余维面正交的 $n$ 维子空间之并，且 $k\to\infty$ 时 $C_n=\Theta(1)$。两条推论值得记：衰减是**各向异性**的——「在 $\mathbb R^d$ 的几乎所有方向上衰减为 $k^{-d-1}$，但在与界定线性区域的 $d-1$ 维面正交的特定方向上可以慢到 $k^{-2}$」；而分子被 $N_fL_f$ 控制（$N_f$ 为线性区域数，$L_f=\max_\epsilon\|W_\epsilon\|$ 为 Lipschitz 常数），

$$
L_{f}\le\prod_{k=1}^{L+1}\|W^{(k)}\|
\le\|\theta\|_{\infty}^{L+1}\sqrt d\prod_{k=1}^{L}d_{k} .
$$

动力学上的读法是原文自己给的：$L_f$ 被参数范数控制，而参数范数在梯度下降中只能逐渐增长，**因此高频只能在优化的后期被学到**。

许志钦等人的频率原则则直接陈述为：「深度网络在训练过程中常常按响应频率从低到高地拟合目标函数。」这里的「频率」指**响应**频率，即输入到输出映射的频率内容，而不是图像意义上的频率；度量方式是离散 Fourier 变换 $\hat f_k=\frac1n\sum_if(x_i)e^{-\mathrm i2\pi ik/n}$ 与逐频相对误差

$$
\Delta_F(k)=\frac{|\hat h_k-\hat f_k|}{|\hat f_k|} .
$$

这一支与经典迭代解法的反差是本页各篇最相关的一点，原文写得很直接：「深度网络的这一频率原则与 Jacobi 方法的行为相反——后者作为一种常规迭代数值格式，在各类科学计算问题上对**更高**频率收敛更快。」

其定理 1（理想化设定：单隐层、$\sigma=\tanh$）说：若目标满足 $|\hat f(k_1)|>0$、$|\hat f(k_2)|>0$、$|k_2|>|k_1|>0$，则存在 $c,C>0$，使得对充分小的 $\delta$，

$$
\frac{\mu\bigl(\{W:|\partial L(k_{1})/\partial\theta_{lj}|>|\partial L(k_{2})/\partial\theta_{lj}|\ \forall\,l,j\}\cap B_{\delta}\bigr)}{\mu(B_{\delta})}
\ \ge\ 1-C\exp(-c/\delta),
$$

$B_\delta$ 是权空间原点处半径 $\delta$ 的球、$\mu$ 是 Lebesgue 测度。定理 2 把结论从逐坐标梯度升级到损失的实际下降速率。**注意其作用域：小权重（初始化附近）、单隐层、$\tanh$——这不是全局陈述。** 论文把机制归给激活函数的光滑性：对 $\tanh$ 有 $|\hat h(k)|\le C\sum_j\frac{|a_j|}{|w_j|}\exp\bigl(-\frac{\pi k}{2|w_j|}\bigr)$，即固定权重下网络自身的谱在 $k$ 上**指数**衰减。

### 多尺度网络：把高频搬到低频再学

MscaleDNN 的构造是本页多篇的共同底座。取带限目标 $\operatorname{supp}\hat f\subset B(K_{\max})$，把 $k$ 空间分成 $M$ 个同心环带

$$
A_{i}=\{\bm k:(i-1)K_{0}\le|\bm k|\le iK_{0}\},\qquad K_{0}=K_{\max}/M,
$$

于是 $\hat f=\sum_i\chi_{A_i}\hat f=:\sum_i\hat f_i$。再**按半径降尺度**：令 $\hat f_i^{\rm scale}(\bm k):=\hat f_i(\alpha_i\bm k)$、$\alpha_i>1$，则 $f_i^{\rm scale}(\bm r)=f_i(\bm r/\alpha_i)$ 且

$$
\operatorname{supp}\hat f_{i}^{\rm scale}\subset
\Bigl\{\tfrac{(i-1)K_{0}}{\alpha_{i}}\le|\bm k|\le\tfrac{iK_{0}}{\alpha_{i}}\Bigr\},
$$

只要 $\alpha_i$ 够大这就是一个**低**频带。按频率原则，普通网络学 $f_i^{\rm scale}$ 很快；于是 $f(\bm r)\sim\sum_{i=1}^{M}h_i(\alpha_i\bm r,\theta^{n_i})$。原文明确说这段推导**本身不是算法**——$d$ 维的频率选择卷积仍受维数诅咒——它只是「给出了函数空间的一个可信形式」。**按半径而非按坐标缩放，正是 MscaleDNN 区别于更早的 PhaseDNN 之处**：后者沿每个坐标做相移，因而继承了维数诅咒。

具体架构有两种。MscaleDNN-1 把第一隐层分成 $N$ 段，第 $i$ 段接收输入 $a_i\bm x$：

$$
f_{\bm\theta}(\bm x)=\bm W^{[L-1]}\sigma\circ\bigl(\cdots(\bm W^{[0]}(\bm K\odot\bm x)+\bm b^{[0]})\cdots\bigr)+\bm b^{[L-1]},
$$

$$
\bm K=(\underbrace{a_{1},\dots,a_{1}}_{\text{第 1 段}},\dots,\underbrace{a_{N},\dots,a_{N}}_{\text{第 }N\text{ 段}})^{\top},
\qquad a_i=i\ \text{ 或 }\ a_i=2^{i-1} .
$$

「与普通全连接网络唯一的差别就是第一隐层的输入。」MscaleDNN-2 则是 $N$ 个**子网络之和**，即 $W^{[1]}$ 到 $W^{[L-1]}$ 为块对角。原文的经验判断是 MscaleDNN-2 更好且连接数远少，还便于动态增删尺度，因此后续实验一律用它。激活函数取紧支撑的形式（动机来自小波理论中的紧支撑尺度函数，使得把激活按 $\alpha$ 缩放就把其频率内容按 $\alpha$ 缩放）：

$$
\mathrm{sReLU}(x)=\mathrm{ReLU}(-(x-1))\times\mathrm{ReLU}(x)\quad(\text{支撑 }[0,1]),
$$

并用 $(\mathrm{sReLU})^2$、$(\mathrm{sReLU})^3$ 得到 $C^1$、$C^2$ 光滑性，以及二次样条型的

$$
\phi(x)=\mathrm{ReLU}(x)^{2}-3\,\mathrm{ReLU}(x-1)^{2}+3\,\mathrm{ReLU}(x-2)^{2}-\mathrm{ReLU}(x-3)^{2}.
$$

### 一个方向相反的提醒

上面这些都是在**回归**损失下说的。**但这个方向并非绝对**：编号 101 指出在 PDE 残差型损失下它会整体反转，详见下文 101 一节。这是本页最重要的一条判断，读其余各篇时都应当带着它。

## 81：先给出误差理论，再据它设计

### 直觉

多尺度网络好用，但它的好用建立在降尺度因子 $a_i$ 选对上，而 $a_i$ 通常被先验固定为 $2^{i-1}$。对 $f(x)=\sin(\pi x)+\sin(100\pi x)$，理想尺度是 $a_1=\pi$、$a_2=100\pi$——实践中不可得，因为要知道它就得先知道解。

带（随机）Fourier 特征的网络看上去能解决这个问题：把 $\sin(\bm k\cdot\bm x)$、$\cos(\bm k\cdot\bm x)$ 直接喂进网络，只要 $\bm k$ 对了，网络等于只需拟合一个常数。问题是**对不上就更糟**：论文用 $k=38\pi$ 对目标频率 $40\pi$ 的例子表明，此时带 Fourier 特征的网络比普通网络还差。

这篇论文的处理分三步。第一步是给出误差理论，说清楚降尺度为什么有用、以及用真实主模态作特征时误差界为什么不含频率。第二步是**混合嵌入**：把线性降尺度项和正弦项拼在一起，让频率对不上时还有线性项兜底。第三步是**后验捕获**：先训一个粗解，对它做离散 Fourier 变换，把能量占比最大的那几个模态取出来当作新的特征，再据此重建网络。整个循环用「捕获集合不再变化」当停机准则。

### 问题设定

多尺度网络的拟设为

$$
f_{\theta}(\bm x)=\sum_{i=1}^{W} a_i^{d}\, f_{\theta_i}(a_i\bm x),
$$

由把紧支撑 Fourier 变换按同心环带 $\mathbb K_i=\{\bm k:(i-1)K_0\le|\bm k|\le iK_0\}$（$K_0=K_{\max}/W$）分解、每带按 $a_i>1$ 降尺度得到，而 $a_i$ 通常先验固定为 $2^{i-1}$。随机 Fourier 特征映射为 $\beta[\bm A](\bm x)=[\sin(\bm A\bm x);\cos(\bm A\bm x)]$，$\bm A\in\mathbb R^{m\times d}$ 的元素抽自 $\mathcal N(0,\sigma^2)$。PDE 损失取

$$
\min_{\theta\in\Theta}\ w_r\mathcal L_r(\theta)+w_b\mathcal L_b(\theta),
\qquad
\mathcal L_r=\sum_{i=1}^{N_r}\bigl|\mathcal N(\bm x^i_r;u_{\text{net}})\bigr|^2,
\quad
\mathcal L_b=\sum_{i=1}^{N_b}\bigl|\mathcal B(\bm x^i_b;u_{\text{net}})\bigr|^2 .
$$

注意这两项按原文是**求和**而非取均值。

### 推导

**混合特征嵌入的动机是一处具体计算。** 把 $f(x)=\sin(40\pi x)$ 经 $F(y)=\sin(\hat k\arcsin y)$（$\hat k=40\pi/k$）逼近时，

$$
F'(y)=\frac{\hat k\cos(\hat k\arcsin y)}{\sqrt{1-y^2}},
\qquad
\lim_{y\to1}|F'(y)|=\infty ,
$$

因此只要 $|\hat k|\ne1$ 就有 $\|F\|_{C^1[-1,1]}=\infty$，定理 3.3 的假设失效。**这就是纯 Fourier 特征在频率不匹配时退化的精确原因**，而线性项的加入正是为了补上这一处。于是新的输入映射把多尺度网络的线性降尺度与正弦 Fourier 特征合起来：

$$
\Phi[\bm k](\bm x)=\begin{bmatrix}\bm k\cdot\bm x\\ \cos(\bm k\cdot\bm x)\\ \sin(\bm k\cdot\bm x)\end{bmatrix}.
$$

**后验捕获主模态。** 对当前网络解做离散 Fourier 变换

$$
\hat f_{\text{net},\bm k}=\int_{[0,1]^d}f_{\text{net}}(\bm x)e^{-\mathrm i2\pi\bm k\cdot\bm x}\,\mathrm d\bm x,
$$

保留模最大的 $N_0$ 个系数，其中 $N_0$ 由能量占比条件确定：

$$
\sum_{j=1}^{N_0}\bigl|\hat f_{\text{net},\bm k_j}\bigr|^2
\ \ge\ (1-\delta)\,\|f_{\text{net}}\|^2_{L^2([0,1]^d)},
\qquad 0\le\delta<1 .
$$

这一步的工作假设是 $\|f_{\text{net}}-f\|_{L^2([0,1]^d)}\le\epsilon$ 且 $0<\epsilon\ll\|f\|_{L^2}$——也就是说粗解已经足够好到它的谱能代表真解的谱。

**据此按两条准则重建网络。** 准则 A：若捕获模态数 $N_0\le M_0$（子网络数），建 $N_0$ 个子网络，第 $j$ 个的输入为 $\Phi[\bm k_j](\bm x)$，输出

$$
y=\sum_{j=1}^{N_0} h_j\,y_j,\qquad h_j=\hat u_{\text{net},0,\bm k_j},
$$

系数被复用以加速收敛。准则 B：若 $N_0>M_0$，则**不扩网络**，而把 $\mathbb B_1=\{\bm k_1,\dots,\bm k_{N_0}\}$ 分成 $M_0$ 块（通常 $\mathbb B_1^1=\{\bm k_1,\dots,\bm k_{\lfloor N_0/M_0\rfloor}\}$，依此类推），第 $j$ 个子网络接收该块内各 $\bm k$ 的 $\Phi[\bm k](\bm x)$，输出

$$
y=\sum_{j=1}^{M_0}h_j\,y_j,\qquad h_j=\sum_{\bm k\in\mathbb B^j_1}\bigl|\hat u_{\text{net},0,\bm k}\bigr| .
$$

一处实现细节值得记录（其注记 4.1）：最后一层写成 $u_{\text{net}}=WG+b$，其中 $W=(h_1,h_2,\dots)$ 与 $b=0$ 是**固定的、不可学习的**，论文报告这样既更准也更省。

**算法。** 给定子网络数 $M_0$、初始特征集 $\mathbb B_0=\{2^0,\dots,2^{M_0-1}\}$、自适应总步数 $I$、阈值参数 $\lambda$：反复地训练当前特征集下的多尺度网络、对解做离散 Fourier 变换、把满足 $|\hat u_{\text{net},It,\bm k}|>\lambda\max_{\bm k}|\hat u_{\text{net},It,\bm k}|$ 的 $\bm k$ 放进下一轮特征集、按准则 A 或 B 重建网络。两条停机准则：固定的最大轮数 $I$，以及捕获特征集不再变化（$\mathbb B_{It}=\mathbb B_{It+1}$）。

### 定理

四条定理，都在函数逼近的意义上陈述；本次核对了陈述，未审阅证明。

**定理 3.1** 给出**标准**网络（无缩放映射）的拟合误差界，作为对照基准。

**定理 3.2（降尺度为何有效）。** 设 $f:\mathbb R^d\to\mathbb R$ 是 $(p,C_0)$-光滑的，$p=q+s$、$q\in\mathbb N_0$、$s\in(0,1]$，且 $M$ 充分大。则在特征映射函数类 $\widetilde{\mathcal F}(1,L,\iota,\delta)$ 中存在带**降尺度**映射 $\Phi(\bm x)=k\bm x$ 的网络 $f_{\text{net}}$，使

$$
\|f_{\text{net}}-f\|_{\infty,[-h,h]^d}
\ \le\ \frac{C_2\bigl(\max\{kh,\,C_1\}\bigr)^{5q+3}}{M^{2p}} .
$$

与定理 3.1 相比，要点是**缩小 $kh$ 就缩小了分子**。

**定理 3.3** 处理 $f=F(\sin(\bm k\cdot\bm x))$ 型的复合，其关键假设是 $\|F\|_{C^1}<\infty$——上面那处 $\lim_{y\to1}|F'(y)|=\infty$ 的计算正是在说这个假设何时失效。

**定理 3.4（带限目标的频率无关界）。** 设 $f\in L^2([-\pi,\pi]^d)$ 带限，

$$
f(\bm x)=\sum_{\bm k\in\mathbb B}\bigl[b_{\bm k}\cos(\bm k\cdot\bm x)+c_{\bm k}\sin(\bm k\cdot\bm x)\bigr],
\qquad |\mathbb B|=N,
$$

系数有界。取 Fourier 特征 $\Phi(\bm x)=(\dots,\sin(\bm k\cdot\bm x),\cos(\bm k\cdot\bm x),\dots)_{\bm k\in\mathbb B}$ 以及 $L=C_4(1,q)$、$\iota=C_5(1,q,M)$、$\delta=C_6(1,M,1,1)$，则存在 $f_{\text{net}}\in\widetilde{\mathcal F}(2N,L,\iota,\delta)$ 使

$$
\|f_{\text{net}}-f\|_{\infty,[-\pi,\pi]^d}\ \le\ \frac{C_2\sqrt N}{M^{2p}},
\qquad M^{2p}\ge\max\{2c_2,c_3,12\} .
$$

**这个界不显含频率。** 这正是「用实际主模态构造嵌入」的理论依据。

### 数值实验

第 5 节覆盖四类方程：

| 小节 | 方程                              | 记录到的自适应轨迹                                                  |
| ---- | --------------------------------- | ------------------------------------------------------------------- |
| 5.1  | Poisson 方程                      | 捕获特征集 $\{2\pi,4\pi,200\pi,202\pi\}$，按准则 A 重建为四个子网络 |
| 5.2  | 热方程                            | 按准则 B 重建                                                       |
| 5.3  | 波方程                            | —                                                                   |
| 5.4  | 半经典极限附近的 Schrödinger 方程 | 按准则 B 重建                                                       |

5.1 节另外报告了一处消融：把最后一层的 $W$、$b$ **固定**比让它们可学更好，既更准也更省。

论文的头条结论是：频率自适应版本相对标准多尺度网络把精度提高**两到三个数量级**。它的参照点也值得记：原始 MscaleDNN 文献（Cai-Xu）报告的误差典型量级是 $10^{-2}\sim10^{-1}$，这正是本文要做后验精化的动机。

这组实验建立的是：**从解自身的谱里读出频率，比先验假设频率更可靠**，并且这个循环在四类方程上都收敛到稳定的特征集。它没有建立的是理论与实验的对接——定理 3.2 与 3.4 是逼近论意义上的存在性界，并不保证优化能找到那个网络，论文也没有测量实际误差与界之间的差距。

> [!warning] 逐例误差数值未转录
> 各例的具体误差数值位于第 5 节的图表中，本次核对未逐项转录，因此本页不列这些数字，只记录已核对的结构、消融与头条量级。

### 与其他论文的关系

这是该组「频谱偏差」一支的起点。编号 101 从另一个角度攻同一问题：不按离散 Fourier 变换捕获的模态重建网络，而是保留一个固定的多尺度随机 Fourier 基库、学一组**依输入而变的注意力权重**，同时保留「从中间解的谱里增量地丰富特征」这一思路。编号 94 是直接应用：它把多尺度网络（而非 PINN）用在半经典 Schrödinger 方程导出的高斯波包常微分方程系统上，也正是本文 5.4 节的例子；两文的前两位作者相同。编号 90 的子域局部缩放基是同一降尺度思想在随机基方法中的对应物，见[[computational-mathematics/paper-notes/scientific-machine-learning/variational-and-basis-networks|变分与基网络一页]]。

## 94：半经典极限——把空间维度整体消掉

### 直觉

$\varepsilon\ll1$ 时半经典 Schrödinger 方程的解在空间与时间上都以 $O(\varepsilon^{-1})$ 的尺度振荡。任何直接离散——有限差分、谱方法、物理信息网络——都要在这个尺度上分辨解，代价随 $\varepsilon$ 变小而爆炸。

但这些振荡的**形状**是已知的：解近似是一个高斯包络乘以一个相位因子 $\exp[\mathrm i(\cdot)/\varepsilon]$，包络中心、宽度、动量、相位都是随时间演化的少数几个标量。把这个形式当作拟设代入方程，就得到关于这几个标量的常微分方程系统，**空间变量 $x$ 被完全消掉**。网络要学的于是不再是一个高频振荡的函数，而是四条随时间演化的曲线；高频结构被解析地放进了重建公式。

代价也很具体：因为这几个标量坐在 $\exp[\mathrm i(\cdot)/\varepsilon]$ 里面，常微分方程解的误差在重建时被放大 $1/\varepsilon$ 倍。所以常微分方程必须解得极准——这正是论文改用多尺度网络而非普通 PINN 的理由。

第二个独立的问题是：更换初值会强制重训，因此需要一个算子学习的表述。论文的答案是把算子分解成三段，只学中间那段有限维的映射。

### 问题设定

目标方程为

$$
\psi_t=\frac{\mathrm i\varepsilon}{2}\Delta\psi-\frac{\mathrm i}{\varepsilon}V(\bm x)\psi,
\qquad
\psi(\bm x,0)=\varphi(\bm x)\exp\bigl(\mathrm i\phi(\bm x)/\varepsilon\bigr),
$$

$\varepsilon$ 是无量纲 Planck 常数。$\varepsilon\ll1$ 时解「在空间与时间上都以 $O(\varepsilon^{-1})$ 的尺度高频振荡」。经典分裂求解器（Strang 分裂在 $\Delta x/\varepsilon$ 上谱精度、时间上 $O(\Delta t^2/\varepsilon)$；Chin-Chen 空间谱精度、时间四阶）的全局截断误差都带一个 $1/\varepsilon$ 因子，因此网格与时间步必须随 $\varepsilon$ 一同缩小。

普通物理信息网络更差。PINN 基线损失为 $\mathcal L(\theta)=w_r\mathcal L_r+w_b\mathcal L_b+w_i\mathcal L_i$，

$$
\mathcal L_r=\Bigl\|\psi_t-\tfrac{\mathrm i\varepsilon}{2}\psi_{xx}+\tfrac{\mathrm i}{\varepsilon}V\psi\Bigr\|^2_{2,\Omega\times[0,T]},
\quad
\mathcal L_b=\|\mathcal P(\psi)-\psi\|^2_{2,\partial\Omega\times[0,T]},
\quad
\mathcal L_i=\bigl\|\psi(\cdot,0)-\varphi e^{\mathrm i\phi/\varepsilon}\bigr\|^2_{2,\Omega},
$$

$\mathcal P$ 为周期边界算子。论文直接给出证据：把朴素网络用于 $(-2\pi,2\pi)\times[0,1]$ 上的该方程，在 $\varepsilon=0.1$ 时相对 $L^2$ 误差**已经超过 1**。

### 推导

**第一步：Heller 高斯波包拟设。** 一维情形为

$$
\psi(x,t)=\exp\left[\frac{\mathrm i}{\varepsilon}
\Bigl(\alpha(t)\bigl(x-q(t)\bigr)^2+p(t)\bigl(x-q(t)\bigr)+\gamma(t)\Bigr)\right].
$$

包络中心在 $q(t)$、标准差正比于 $\sqrt{\varepsilon/\alpha_{\rm im}(t)}$，$x=q(t)$ 处的振荡波长为 $2\pi\varepsilon/p(t)$，而 $\mathrm{Re}\,\alpha$ 在尾部产生更细的振荡。

**第二步：关键约化。** 把拟设代入方程，对二次势与高斯初值它是**精确解**，条件是

$$
\dot q=p,\qquad
\dot p=-V'(q),\qquad
\dot\alpha=-2\alpha^{2}-\tfrac12V''(q),\qquad
\dot\gamma=\tfrac12p^{2}-V(q)+\mathrm i\,\alpha\varepsilon .
$$

**空间变量 $x$ 被完全消掉。** 对非调和势，高斯波包的建模误差是 $O(\sqrt\varepsilon)$，即它随 $\varepsilon\to0$ **变好**（这一条引自 Heller 的工作，不是本文新证）。

**第三步：误差放大这一关键事实。** 因为常微分方程的未知量坐在 $\exp[\mathrm i(\cdot)/\varepsilon]$ 里面，「若 PINN 解常微分方程系统的总误差为 $\mathcal E_t$，则重建出的解的总误差为 $\mathcal E_t/\varepsilon$」。等价地说，用 $k$ 阶时间积分器解常微分方程会给出 $(\Delta t)^k/\varepsilon$ 量级的总截断误差。**这个 $1/\varepsilon$ 放大就是常微分方程必须解得极准的原因**，也是论文改用多尺度网络的理由。

**第四步：常微分方程系统的损失。** 网络输出 $\bm y_\theta=(q_{\theta_1},p_{\theta_2},\alpha_{\theta_3},\gamma_{\theta_4})$ 只以 $t\in[0,T]$ 为输入，损失配点四条残差，例如 $\gamma$ 那一项为

$$
\bigl|\dot\gamma_\theta(t_i)-\tfrac12p_\theta^2(t_i)+V(q_\theta(t_i))-\mathrm i\alpha_\theta(t_i)\varepsilon\bigr|^2 ,
$$

重建场为

$$
\psi(x,t;\theta)=\exp\left[\frac{\mathrm i}{\varepsilon}
\Bigl(\alpha_{\theta_3}(t)(x-q_{\theta_1}(t))^2+p_{\theta_2}(t)(x-q_{\theta_1}(t))+\gamma_{\theta_4}(t)\Bigr)\right].
$$

**第五步：多尺度网络的具体取法。** 论文重述了 MscaleDNN 的降尺度构造并写成

$$
f_\theta(\bm z)=\sum_{i=1}^{M}a_i^{d+1}f_{\theta_i}(a_i\bm z),
$$

但**本文实验并不用通常的几何尺度 $a_i=2^{i-1}$，而用 100 个线性等距尺度 $\{a_i\}=\{0.1,0.2,\dots,10.0\}$**。理由是观察到 $q(t)$ 的 Fourier 谱有长尾，几何尺度覆盖不够密。

**第六步：算子分解。** 与其直接学 $\mathcal G:\psi(\cdot,0)\mapsto\psi(\cdot,t)$，不如把它穿过波包参数分解：

$$
\mathcal G:\ \psi(x,0)\ \longrightarrow\ \bm y(0)
\ \overset{\mathcal G'}{\longrightarrow}\ \bm y(t)\ \longrightarrow\ \psi(x,t),
\qquad \bm y=(q,p,\alpha,\gamma),
$$

只学**有限维**的 $\mathcal G':\bm y(0)\mapsto\bm y(t)$。**关键在于分支输入是 $\bm y(0)$、主干输入只有 $t$，两者都与 $x$ 无关**，这正是把振荡的 $x$ 依赖从学习问题中移除的地方。架构为

$$
\mathcal G'^{(i)}_\theta[\bm y(0)](t)=\sum_{j=J_{i-1}+1}^{J_i}b_j[\bm y(0)]\,c_j(t),
\quad i=1,\dots,I,
$$

即把 $J$ 个基指标划成 $I$ 块、每块对应一个输出分量。损失完全**无数据**：

$$
\hat{\mathcal L}(\theta)=
\underbrace{\frac1N\sum_{i=1}^{N}\sum_{j=1}^{I}
\bigl|\mathcal G'^{(j)}_\theta[\bm y^{(i)}(0)](0)-y^{(i)}_j(0)\bigr|^2}_{\hat{\mathcal L}_{\rm boundary}}
+\underbrace{\frac{1}{NQ}\sum_{i=1}^{N}\sum_{j=1}^{Q}\sum_{l=1}^{I}
\Bigl|\mathcal N^{(l)}\bigl(\mathcal G'_\theta[\bm y^{(i)}(0)](t_j^{(i)})\bigr)\Bigr|^2}_{\hat{\mathcal L}_{\rm physics}},
$$

其中 $(\mathcal N^{(1)},\dots,\mathcal N^{(I)})$ 是常微分方程系统的残差算子。

一般（非高斯）初值另有两条出路：附录 A 推广到 **Hagedorn 波包**（输出 $q,p,Q,P,S$），附录 B 用**高斯束分解**把 WKB 型初值拆成 $\mathfrak N$ 束。这两处附录的存在与覆盖范围已核对，其内部公式本次未细读。

### 定理

**本文没有证明定理。** 高斯波包的 $O(\sqrt\varepsilon)$ 建模误差与 $\mathcal E_t/\varepsilon$ 放大都引自前人工作（Heller 及相关文献），不是本文新证。论文的小结明说：「一个解释多尺度网络为何在求解常微分方程系统上优于 PINN 的稳健理论框架仍然缺失。」

另有一条明确的局限（其注记 3.1）：单个高斯波包只在 **Ehrenfest 时间**内准确，因此「本文提出的神经网络求解器不适用于长时间模拟」，即便常微分方程本身解得仍准。

### 数值实验

参考解：非调和势用四阶格式，调和势用波包常微分方程上的 RK4。度量为 $\Omega\times[0,T]$ 上的相对 $L^2$。激活取 **Softened Fourier Mapping** $\sigma(x)=0.5\sin x+0.5\cos x$，Glorot 正态初始化，Adam、学习率 $10^{-3}$、10 万轮、批量 1000。

**例 1：自由粒子的常微分方程系统。** $V\equiv0$、$\varepsilon=0.01$，精确解 $q=1+2t$、$p=2$、$\alpha=(2t+\mathrm i)/(4t^2+1)$、$\gamma=2t-0.005\arctan(2t)+\tfrac14(\log(4t^2+1)-\log(200/\pi))\mathrm i$。PINN 取 $[1,100,400,400,400,400,6]$，多尺度网络取 $[100,400,400,400,400,6]$ 配 100 个线性尺度。

| 分量     | PINN                 | MscaleDNN            |
| -------- | -------------------- | -------------------- |
| $q$      | $3.287\times10^{-5}$ | $2.928\times10^{-7}$ |
| $p$      | $2.585\times10^{-4}$ | $4.141\times10^{-6}$ |
| $\alpha$ | $3.105\times10^{-3}$ | $1.941\times10^{-5}$ |
| $\gamma$ | $2.864\times10^{-4}$ | $4.161\times10^{-6}$ |

诊断也很清楚：$q(t)$ 的离散 Fourier 变换有长尾谱，PINN 的**误差谱与之同形**，而多尺度网络的误差谱要均匀得多。嵌入数从 10 增到 100 时精度持续改善。

**例 2：一维调和势 $V=\tfrac12x^2$。** $\varepsilon\in\{4/25,1/25,1/100,1/400,1/1600,1/6400\}$，PINN 误差从 $1.522\times10^{-3}$ 增到 $6.877\times10^{-1}$，多尺度网络从 $3.192\times10^{-5}$ 增到 $7.956\times10^{-3}$。误差随 $\varepsilon$ 变小而**线性**增长，与 $1/\varepsilon$ 放大以及调和情形下建模误差为零这两点一致。采样密度研究表明 6400 个点已足够，此后一律用 10000 个。**把标准多尺度网络直接用在原方程上，在 $\varepsilon=1/32$ 时就已失败**——这是把波包约化单独隔离出来的干净消融。

**例 3：一维扭转势 $V=1-\cos x$。** 大 $\varepsilon$ 时 PINN 与多尺度网络几乎一样，因为 $O(\sqrt\varepsilon)$ 的**建模**误差占主导；$\varepsilon<1/100$ 后多尺度网络才拉开差距（$\varepsilon=1/1600$ 时 $3.408\times10^{-3}$ 对 $1.048\times10^{-2}$）。这是「建模误差加求解误差」这一两项预算的好例证。

**例 4：二维与四维。** 二维时 $q,p$ 成向量、$\alpha$ 成复对称矩阵 $A$，输出维数 12。四维时 $V(\bm x)=\tfrac12\bm x^{T}A_p\bm x$，$A_p$ 对角为 $(1.0,2.0,2.0,3.0)$、所有非对角元为 $0.2$，$\bm q(0)=(1.3,0,-1,0)^T$、$\bm p(0)=(0,1.3,0,1)^T$；$A$ 需要 20 个输出分量，因此用四个独立网络分别处理 $(\bm q,\bm p)$、$A$、$\gamma_{\rm re}$、$\gamma_{\rm im}$。

| $\varepsilon$（四维） | PINN                | MscaleDNN           |
| --------------------- | ------------------- | ------------------- |
| $0.01$                | $2.08\times10^{-4}$ | $7.48\times10^{-5}$ |
| $0.001$               | $2.07\times10^{-3}$ | $6.47\times10^{-4}$ |

**例 5：DeepONet，一维扭转势。** $\varepsilon=0.01$、$T=1$；$q(0)\sim\mathcal U[0.8,1.8]$，$p(0),\alpha_{\rm re}(0)\sim\mathcal U[-0.5,0.5]$，$\alpha_{\rm im}(0)\sim\mathcal U[0.5,1.5]$，$\gamma_{\rm re}(0)=0$，$\gamma_{\rm im}(0)$ 由归一化定；$I=6$、$J_1=\dots=J_I=100$；分支 $[6,100,100,100,100,600]$、主干 $[1,100,100,100,100,100]$；$N=2000$ 个初值函数、$Q=500$ 个时间点、20 万次迭代。在 $q(0)=\pi/2,1,1.5$ 处的误差为 $2.69\times10^{-2}$、$2.39\times10^{-2}$、$1.67\times10^{-2}$；在 $N_{\rm test}=100$ 上均值 $2.44\times10^{-2}$、标准差 $1.54\times10^{-2}$。

**例 6：DeepONet 与 RK4 的计时。** 调和情形下：

| $N_{\rm test}$ | RK4     | DeepONet |
| -------------- | ------- | -------- |
| 10000          | 11.36 s | 4.15 s   |
| 20000          | 22.61 s | 4.72 s   |
| 40000          | 45.63 s | 5.44 s   |

**DeepONet 的推断代价几乎不随 $N_{\rm test}$ 增长。** 另有一组：为 10000 个算例各生成 100 个网格点，DeepONet 用 4.58 s，而 $\Delta t=0.01$ 的 RK4 用 67.85 s。计时对比是诚实校准过的——RK4 的 $\Delta t=0.1$ 是选来让其误差与 DeepONet 相当的。分量误差（调和一维）的均值为 $q$ $5.33\times10^{-5}$、$p$ $6.95\times10^{-5}$、$\alpha$ $4.04\times10^{-4}$、$\gamma$ $1.85\times10^{-4}$，而 **$\psi$ 为 $1.93\times10^{-2}$**——从常微分方程参数到 $\psi$ 的这两个数量级跳跃，恰好就是 $1/\varepsilon=100$ 的放大。

二维 DeepONet 用 $V=(x_1^2+x_2^2)/2$；附录 A 用 Hagedorn 波包重做同一研究；附录 B 研究 WKB 型初值的高斯束分解并报告运行时间与精度随束数 $\mathfrak N$ 的变化。

这组实验建立的是：**把振荡结构解析地移出学习问题之后，剩下的低维问题可以被解得非常准**，且多尺度网络在这个低维问题上仍显著优于 PINN。它没有建立的是长时间行为——Ehrenfest 时间之外单个波包本身就失效，与网络无关。

### 与其他论文的关系

多尺度网络与编号 101 的多尺度随机 Fourier 基库攻的是同一个频谱偏差：编号 94 用**先验固定**的 100 个线性尺度，编号 101 则通过交叉注意力学基库上的幅度、还能把离散 Fourier 变换找到的频率**加进去**。两文都以目标或误差的离散 Fourier 变换作为诊断工具；编号 101 是更通用的机器，编号 94 是特定领域的应用。编号 90 的子域局部缩放基是同一降尺度思想在随机基方法中的对应物。

本文用的物理信息 DeepONet 正是编号 95 推广为置换不变、带不确定性版本的那个算子学习底座，也是编号 103 用算子分裂架构替换掉的那个底座；三者都是面向演化问题的、物理信息的、无数据算子学习。而「把 PDE 约化成小的常微分方程系统、再学这个系统的流映射」这一策略，在结构上与编号 107 把 Fokker-Planck 算子学习约化成带线性化随机微分方程基分布的转移密度是同一个动作，见[[computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning|算子学习与不确定性一页]]。

## 101：固定基库，可学权重

### 直觉

这篇论文最值得记录的一点不是架构，而是它 3.1 节的一个观察：**在物理信息训练下，频谱偏差的方向会反过来。**

理由在 Fourier 端的缩放。由 $\widehat{\partial_xu}(k)=\mathrm ik\,\hat u(k)$，Deep Ritz 型能量项 $\|\partial_xu\|_{L^2}^2$ 给模态 $k$ 的权重是 $k^2$；由 $\widehat{u_{xx}}(k)=-k^2\hat u(k)$，平方 PDE 残差给模态 $k$ 的权重则**像 $k^4$**。也就是说微分算子本身**放大**高频，于是在残差损失下被欠解的反而是**低频**部分。

论文把这一点测了出来：对 $u=\sin(\pi x)+\sin(5\pi x)+\sin(20\pi x)$ 在 $k\in\{1,5,20\}$ 上看逐频相对误差 $\Delta_F(k)=|\hat u_k^{\rm pred}-\hat u_k|/|\hat u_k|$，纯回归损失下确实低频先收敛，而 Deep Ritz 与 PINN 损失下高频衰减得更快。

架构那一侧的直觉则是：既然事先不知道解有哪些频率，就准备一个**足够宽**的固定频率基库，让模型自己决定在定义域的哪个部分强调哪一段频带。做法是交叉注意力——查询来自网络的隐状态，键与值来自频率基库，softmax 权重依输入而变。基库本身不训练，只训练一个作用在幅度上的衰减包络。若事后发现基库里缺了关键频率，就把离散 Fourier 变换找到的模态**追加**成新的词元，而不重建基库。

> [!note] 附录 A 未核对
> 正文的 Fourier 缩放论证（$k^2$ 与 $k^4$）已核对。附录 A 中「微分算子放大高频」的那段分析本次未读完，其细节按未核对处理。

### 问题设定

频谱偏差指网络「快速收敛到低频成分……却难以表示高频或强振荡特征」。论文对已有修法的批评很具体：随机 Fourier 特征、MscaleDNN、FMMNN、相移法、随机特征法「通常依赖预先指定的频率基或规定好的多尺度变换，这可能限制它们对实例相关的谱需求的适应性」；FG-PINN 则「强烈依赖源项或初边值数据中是否存在有信息量的高频内容」。论文的论点是：「克服频谱偏差不仅需要更丰富的频率字典，还需要能够依据输入与演化中的解结构动态地路由与重加权谱成分的架构机制。」

### 推导

**第一步：多尺度随机 Fourier 基库。** 基频

$$
\omega_m\sim\mathcal N\bigl(\bm 0,\ \sigma^{-2}I_{d_{\rm in}}\bigr),\qquad m=1,\dots,M_{\rm base},
$$

按二进尺度 $\widetilde\omega_{m,k}=2^k\omega_m$（$k=0,\dots,K$）展成 $\overline\Omega\in\mathbb R^{M\times d_{\rm in}}$、$M=M_{\rm base}(K+1)$，相位 $b_{m,k}\sim\mathrm{Uniform}(0,2\pi)$ 一次抽定。特征映射带一个**可学的幅度包络**：

$$
\phi(x)=\sqrt{\tfrac1M}\Bigl[a_{m,k}\cos\bigl(\widetilde\omega_{m,k}^{\top}x+b_{m,k}\bigr)\Bigr]_{(m,k)}\in\mathbb R^M,
\qquad
a_{m,k}=\exp\bigl(-\beta\|\widetilde\omega_{m,k}\|_2\bigr).
$$

**频率本身固定且不可学，唯一可学的是衰减率 $\beta\ge0$**（经 softplus 参数化保正）。这与编号 81 的取舍正好互补：81 让特征去对齐真实频率，101 用一个足够宽的固定基库配一个可学的衰减包络。

**第二步：把特征重排成词元。** 取词元宽度 $d_q\mid M$，把 $\phi(x)\in\mathbb R^{M}$ 重排成 $H(x)\in\mathbb R^{N_{\rm tok}\times d_q}$、$N_{\rm tok}=M/d_q$。「$d_q$ 的一个简单取法是 $d_q=M_{\rm base}$，使 $H(x)$ 的每一行对应某一个尺度」——**即一个尺度一个词元，这就使得对词元的注意力字面上就是对频带的注意力。**

**第三步：交叉注意力残差栈。** 初始隐状态

$$
Q^{(0)}(x)=\sigma\bigl(W^{(0)}\psi(x)+b^{(0)}\bigr)\in\mathbb R^{d_q},
$$

其中 $\psi(x)=x$ 给出 **NN-CA** 变体、$\psi(x)=\phi(x)$ 给出 **RFF-CA** 变体。记 $Q_l=Q^{(l)}W_Q^{(l)}$、$K_l=H(x)W_K^{(l)}$、$V_l=H(x)W_V^{(l)}$，$W_\bullet^{(l)}\in\mathbb R^{d_q\times d_q}$，则对 $l=0,\dots,L-1$，

$$
\mathrm{CA}\bigl(Q^{(l)},H\bigr)=\mathrm{softmax}\Bigl(\frac{Q_lK_l^{\top}}{\sqrt{d_q}}\Bigr)V_l,
$$

$$
\widetilde Q^{(l)}=Q^{(l)}+\mathrm{CA}\bigl(Q^{(l)},H\bigr),
\qquad
Q^{(l+1)}=\widetilde Q^{(l)}+\sigma\bigl(W^{(l)}\widetilde Q^{(l)}+b^{(l)}\bigr),
$$

输出 $u_\theta(x)=W_{\rm out}Q^{(L)}(x)+b_{\rm out}$。查询来自**隐状态**，键与值来自**频率基库**——所以是交叉注意力而非自注意力。要点在于 softmax 权重依输入而变，因此定义域的不同区域可以强调不同频带。

> [!warning] 预印本的一处笔误
> 描述 NN-CA 与 RFF-CA 两种情形的那句话把 $\psi(x)=\phi(x)$ 写了两遍；按上下文，第一种情形应为 $\psi(x)=x$。

**第四步：自适应频率增强（AFE）。** 先训一个初步模型 $u_\theta^{(0)}$，在周期区域 $\Omega$ 的均匀网格上做离散 Fourier 变换得 $\hat u_{\theta,k}$（指标集 $B$），取 $\zeta=\max_{k\in B}|\hat u_{\theta,k}|$ 并按相对阈值挑出

$$
\mathcal K_{\rm post}=\bigl\{k\in B:\ |\hat u_{\theta,k}|>\lambda\zeta\bigr\},\qquad 0<\lambda<1 .
$$

据此造确定性频率 $\omega_k^{\rm post}=2k\pi$ 与后验特征

$$
\phi_{\rm post}(x)=\sqrt{\tfrac{2}{M_{\rm post}}}\cos\bigl(\Omega_{\rm post}x+b_{\rm post}\bigr),
\qquad
\Omega_{\rm post}=[\omega_k^{\rm post}]_{k\in\mathcal K_{\rm post}},\ M_{\rm post}=|\mathcal K_{\rm post}| ,
$$

重排成 $H_{\rm post}$（若 $d_q\nmid M_{\rm post}$ 则补零），在**词元维度**上拼接成 $H_{\rm aug}=[H_{\rm base};H_{\rm post}]$。新词元不是硬接入，而是经一个加性对数掩码平滑放开：

$$
A^{(l)}=\frac{Q_lK_l^{\top}}{\sqrt{d_q}}+\mathcal M^{(l)},
\qquad
\mathcal M^{(l)}=[\,\bm 0;\ \eta_l\bm 1\,],\quad \eta_l\le0,
$$

零块作用在 $H_{\rm base}$ 上、常数块 $\eta_l$ 作用在 $H_{\rm post}$ 上，再让 $\eta_l\uparrow0$ 逐步释放注入的词元。因此它**增广**基库而不重建基库，主干保持不动——这正是与编号 81 的分野所在。注记 2.1 记下两处实践上的注意事项：区域复杂时先嵌入超立方体再做掩码；高维时从一维分量函数提取模态，而不做完整的 $d$ 维离散 Fourier 变换。

**第五步：两网络的 PDE 表示。** 对 $\mathcal N[u]=f$ 于 $\Omega$、$u=g$ 于 $\partial\Omega$，取

$$
u(\bm x;\theta)=u_h(\bm x;\theta_h)+\alpha\,u_\ell(\bm x;\theta_\ell),
$$

$u_h$ 是上述交叉注意力网络（管高频），$u_\ell$ 是普通全连接网络（管低频），边界条件按 $u_h=g$、$u_\ell=0$ 分开承担。损失为

$$
L=\int_\Omega\bigl(\mathcal N[u_h+\alpha u_\ell]-f\bigr)^2\rho_r\,\mathrm d\bm x
+\gamma\int_{\partial\Omega}\bigl((u_h-g)^2+\alpha^2u_\ell^2\bigr)\rho_b\,\mathrm ds .
$$

混合系数 $\alpha$ 有两种取法：可训练标量，或**闭式最优**。后者由 $\partial L/\partial\alpha=0$ 在 $\mathcal N$ 线性时解出：

$$
\alpha_{\rm opt}=-\frac{\displaystyle\int_\Omega\bigl(\mathcal N[u_h]-f\bigr)\,\mathcal N[u_\ell]\,\rho_r\,\mathrm d\bm x}
{\displaystyle\int_\Omega\bigl(\mathcal N[u_\ell]\bigr)^2\rho_r\,\mathrm d\bm x+\gamma\int_{\partial\Omega}u_\ell^2\,\rho_b\,\mathrm ds},
$$

由配点上的 Monte Carlo 估计。实践中为稳定估计，把 $\alpha_{\rm opt}$ 从边界罚项里去掉，用

$$
\widehat L\approx\frac{1}{N_r}\sum_i\bigl(\mathcal N[u_h+\alpha_{\rm opt}u_\ell]-f\bigr)^2(x_r^{(i)})
+\frac{\gamma}{N_b}\sum_i\Bigl((u_h-g)^2+u_\ell^2\Bigr)(x_b^{(i)}).
$$

### 定理

**没有定理。** 附录 A 被描述为「对微分算子放大高频的一个简单分析」，而不是收敛定理；其细节本次未核对。结论一节说得很直白：「基于交叉注意力的多尺度 Fourier 表示，其逼近与优化动力学的严格分析仍然缺失。」另外自陈的开放问题包括：离散 Fourier 变换引导的增强「对周期或网格友好的设定最自然」，以及向物理信息**算子**学习的推广。

### 数值实验

共享默认超参数：$m_{\rm base}=128$、$K=3$、宽度 $d_q=64$、$n_{\rm heads}=4$、$L=4$ 个块、$\beta_0=0.1$，双精度，Adam 并把梯度裁剪在 1.0。度量为相对 $L^2$，图像任务另加 PSNR 与 HFEN。

**例 1：二维合成回归。** 在 $[-1,1]^2$ 上取三个刻意刁难的目标：$f_1$ 是一个复合函数（用 logistic 角度门做扇区各向异性、用带通门做一个窄的高频环、一个局部化的 Gabor 型螺旋波包、一个弱的星形间断 $0.12\,\mathrm{sign}(r-r_*(\theta))$ 其中 $r_*(\theta)=0.55+0.10\cos5\theta$，再加一个平稳的交叉项）；$f_2=\cos(2\pi[(w_0+w_1r)(x_1\cos\kappa\theta+x_2\sin\kappa\theta)])$、$\kappa=5$、$w_0=4$、$w_1=3$（一个瞬时频率沿半径变化的旋涡）；$f_3=\mathrm{sign}(\sin2\pi f_xx_1\sin2\pi f_yx_2)$、$f_x=f_y=1$（一个不连续的棋盘）。$500\times500$ 网格，批量 4000，Adam 学习率 $2\times10^{-3}$，StepLR 参数对 $f_1,f_2$ 取 $(100,0.5)$、对 $f_3$ 取 $(50,0.5)$，分别训 1000/1000/500 轮。结论：RFF-CA 在三者上的损失与误差都优于 RFF-NN；NN-CA 与 RFF-CA 孰优则依问题而定。论文给的是曲线而非数值表。

**例 2：把图像当函数回归（DIV2K）。** 四张验证图（最大 $2040\times1536$），坐标按 $x_j=\tfrac{j+1/2}{W}\cdot2-1$、$y_i=\tfrac{i+1/2}{H}\cdot2-1$ 映到 $[-1,1]^2$，四倍降采样，全批量，5000 轮。度量含 HFEN，其定义经一个 $15\times15$、$\sigma\approx1.5$ 像素的 Laplacian-of-Gaussian 滤波器：

$$
\mathrm{HFEN}_{\rm rel}=\frac{\|\hat I^{\rm HP}-I^{\rm HP}\|_2}{\|I^{\rm HP}\|_2} .
$$

结论：NN-CA 在四张图上的 HFEN 与相对 $L^2$ 都低于 RFF-NN，PSNR 都更高。

**例 3：自适应频率增强的演示。** 目标 $u(x)=\sin(2\pi\cdot2x)+0.5\sin(2\pi\cdot20x+0.3)+0.5\cos(2\pi\cdot40x-0.2)$ 于 $(0,1)$——刻意取稀疏谱。$N_{\rm train}=2048$、$N_{\rm test}=4096$、$m_{\rm base}=128$、$n_{\rm scales}=1$、$d_q=64$、$L=3$、4 个头，全批量 Adam、学习率 $10^{-3}$ 每 500 轮乘 $0.9$。第一阶段 5000 轮；在 $N_{\rm fft}=4096$ 上做离散 Fourier 变换、取 $\lambda=0.02$，**恰好恢复出主模态 $k=2,20,40$**。第二阶段再训 5000 轮，$\eta_{\rm start}=-6$ 在前 70% 保持不变、后 30% 按余弦释放到 0。结果是词元释放之后误差衰减明显加速。

**例 4：一维 Poisson。** 于 $[-1,1]$，

$$
u=\sin(0.1\pi x)+0.2\sin(\pi x)+0.4\sin(\tfrac\nu3\pi x)+0.6\sin(\tfrac{2\nu}{3}\pi x)+\sin(\nu\pi x),
\qquad \nu=100 .
$$

对 $\alpha\in\{0,1,\text{可学},\text{最优}\}$ 做消融：灵活性增加会单调降低**训练损失**，且可学的 $\alpha$ 给出最低的损失——**但最优缩放才给出最低的相对 $L^2$ 误差**。这是一处显式的「损失与误差背离」，作者自己标出并作了解释；训练早期 $\alpha_{\rm opt}$ 取一个大的负值，注入了一个强的低频修正。架构消融：RFF-CA 与 NN-CA 都比 RFF-NN 衰减更快，NN-CA 总体最好。幅度消融刻意取失配的 $\sigma=0.02$（基频**高于**解的频率）：可学的 $\beta$ 在两个边界罚权 $\lambda=10^3,10^4$ 下都优于固定的 $\beta\equiv0$；而 $\beta=0$ 时更大的边界权反而有帮助，作者的解释是「相当于向损失注入了更多低频的边界信息」。

**例 5：二维 Poisson。** $u=\sin(\mu x_1^2)+\sin(\mu x_2^2)$，每次迭代取 $N_r=10^4$ 个内部点、每边 $N_b=1000$ 个边界点，$2\times10^4$ 步 AdamW，StepLR $(2000,0.5)$，边界权 $\lambda=10^4$，用最优 $\alpha$，同一随机种子。$\mu=50$ 时两者都行、RFF-CA 略好且更稳；**$\mu=100$ 时 RFF-NN 饱和不再下降，而 RFF-CA 继续下降。** 这是「好处随频率增大而增大」最清楚的证据。

**例 6：三维 Poisson-Boltzmann（Deep Ritz）。** $-\nabla\cdot(\epsilon\nabla u)+\kappa u=f$，界面上带跳跃条件 $[u]=0$、$[\epsilon\,\partial u/\partial n]=0$；精确解

$$
u=\frac{e^{\sin\mu x_1+\sin\mu x_2+\sin\mu x_3}}{|x|^2+1}\bigl(|x|^2-1\bigr),
\qquad \mu=15,
$$

$\epsilon\equiv1$，$\kappa$ 在 $\Omega_1$ 取 1、在 $\Omega_2$ 取 5。**区域按构造带几何奇性**：原点处半径 0.5 的球，并上 20 个半径抽自 $[0.1,0.2]$、球心随机落在其表面的小球，再被单位球截断——这些折角很难剖分网格。Ritz 损失

$$
L_{\rm Ritz}=\tfrac12\int_\Omega\bigl(|\epsilon\nabla u_\theta|^2+\kappa u_\theta^2\bigr)
-\int_\Omega fu_\theta+\gamma\int_{\partial\Omega}|u_\theta-g|^2,
$$

$g\equiv0$、$\gamma=10^4$；$L=3$、$\sigma=1$、$\alpha$ 可学、10000 轮、Adam 学习率 $10^{-3}$、StepLR $(1000,0.6)$，每轮 5000 个内部点加 4000 个边界点。结论：RFF-CA 给出明显更准的重建，而「RFF-NN 距真解仍然很远」。

这一组实验的对照是有控制的——RFF-NN 与 RFF-CA「共享同一个多尺度随机 Fourier 词元化器与同样的模型容量」，因此差别可以归给交叉注意力本身。它们建立的是：**在一个足够宽的固定基库上学依输入而变的权重，比把频率写死更稳健**，并且这一好处随频率升高而扩大。它们没有建立的是逼近或优化的理论保证——论文自己也这么说；此外离散 Fourier 变换那一步的适用面被限制在周期或网格友好的设定上。

### 与其他论文的关系

它是编号 94 所用多尺度网络的直接方法学后继：两者都把谱分解成尺度，但编号 94 固定 100 个线性尺度并把子网络输出相加，编号 101 保留一个二进基库并通过注意力**学在哪里用哪些尺度**，还能追加离散 Fourier 变换发现的模态。编号 101 明确把 MscaleDNN 列为「依赖预先指定频率基」的方法。

「频率随机抽定后冻结」这一点与编号 90、102 共享（那里的隐层参数也是抽定后不训练），但读出层不同：这里是训练出来的深注意力栈，那里是最小二乘求解。Xiaodong Feng、Xiaoliang Wan、Tao Zhou 三人还共同署名编号 98 与 105；在这个子群里编号 101 是确定性、逼近论那一支，另两篇是概率那一支。而「用两个网络分管低频与高频、再配一个解析最优的混合系数」这一做法，在概念上与编号 90 把光滑与低正则子域分开是相邻的：都把难的那一部分隔离到专门的逼近空间里，一个在频率上、一个在空间上。

## 103：把算子分裂做成网络，并把时间放回输入

### 直觉

算子分裂是解演化方程的经典手段：把 $u_t=\mathcal Lu+\mathcal Nu$ 中的线性与非线性部分交替推进，每一段都用最适合它自己的方法。DOSnet 的想法是把这个乘积**直接写成网络**——可学的卷积层扮演 $e^{\tau\mathcal L}$，而**激活函数就是非线性子问题的精确流 $e^{\tau\mathcal N}$**，不是 ReLU 或 tanh。于是中间层的激活值是物理上有意义的中间时刻状态，网络不再是黑箱。

但 DOSnet 有两处缺陷。其一，它是数据驱动的，而演化问题的每一对训练样本都要跑一次求解器，代价高得离谱。其二，学到的算子并未被要求满足方程，因此**不能在任意时刻取值**——它只输出终端时刻的解，中间块的输出只有靠与参考数据比对才对应得上中间时刻，且可取回的中间时刻数还受块数限制。

这篇论文的改动很小但很关键：取等步长 $dt=t/K$，并把线性部分的指数换成 $dt$ 的显式二阶 Taylor 展开。这样整块就成了 $t$ 的**显式函数**，于是 $t$ 回到输入里，解可在任意时刻求值，$\partial\hat u/\partial t$ 可以自动微分算出，训练就能改用 PDE 残差、完全不要成对数据。还有一个附带的好处：$t=0$ 时 $dt=0$，每一块都退化为恒等映射，**初始条件自动精确成立**，不需要初值损失项。

### 问题设定

设自治演化方程

$$
u_t=\mathcal Lu+\mathcal Nu=:\mathcal Fu\ \text{ 于 }\Omega\times(0,t^\star],
\qquad u(0,\cdot)=u_0,\qquad \mathcal Bu=0\ \text{ 于 }\partial\Omega,
$$

$\mathcal L$ 线性、$\mathcal N$ 非线性、$\mathcal B$ 为边界算子，均不含 $t$。把 $[0,t^\star]$ 在 $0<t_1<\cdots<t_K=t^\star$ 处分段、步长为 $\tau_1,\dots,\tau_K$，Lie-Trotter（一阶）与 Strang（二阶）分裂分别给出

$$
u(T,\bm x)\approx e^{\tau_K\mathcal N}e^{\tau_K\mathcal L}\cdots e^{\tau_1\mathcal N}e^{\tau_1\mathcal L}u(0,\bm x),
$$

$$
u(T,\bm x)\approx e^{\frac{\tau_K}2\mathcal L}e^{\tau_K\mathcal N}e^{\frac{\tau_K}2\mathcal L}
\cdots e^{\frac{\tau_1}2\mathcal L}e^{\tau_1\mathcal N}e^{\frac{\tau_1}2\mathcal L}u(0,\bm x).
$$

分裂误差来自 $\mathcal L$ 与 $\mathcal N$ 不可交换；两者可交换时它消失。

DOSnet 把这个乘积直接写成网络：$\psi_{\bm\theta_T}=\psi_{\bm\theta_K}\circ\cdots\circ\psi_{\bm\theta_1}$，每个分裂块交替放可学线性层（卷积 $\phi_{\mathcal L_{\bm\theta_{l,i}}}$）与非线性层 $\phi_{\mathcal N_{l,i}}=e^{\tau_{l,i}\mathcal N}$，并约束 $\sum_{l,i}\tau_{l,i}=T$，按数据训练：

$$
\mathfrak L(\bm\theta_T)=\frac1N\sum_n\bigl\|\psi_{\bm\theta_T}(u_0^{(n)})-u^{(n)}(T,\cdot)\bigr\|_{L^2(\Omega)} .
$$

### 推导

**第一步：把块写成 $t$ 的显式函数。** 取等步长 $\tau_{1,1}=\dots=\tau_{1,K}=dt=t/K$，并把线性部分的指数换成 $dt$ 的显式二阶 Taylor 展开：

$$
e^{dt\,\mathcal L}\approx\mathcal I+dt\,\mathcal L_{\bm\theta_i}+\frac{dt^2}{2}\mathcal L_{\bm\theta_i}^2,
$$

$$
u_i(t_{i-1}+dt,\bm x)=\psi_{\bm\theta_i,dt}(u_{i-1})
=e^{dt\,\mathcal N}\circ\Bigl(\mathcal I+dt\,\mathcal L_{\bm\theta_i}+\tfrac{dt^2}{2}\mathcal L_{\bm\theta_i}^2\Bigr)u_{i-1}(t_{i-1},\bm x),
$$

$$
\hat u(t,\bm x)=\Phi_{\bm\theta_T}(u_0,t)=\psi_{\bm\theta_K,dt}\circ\cdots\circ\psi_{\bm\theta_1,dt}\,u_0(\bm x).
$$

**一个漂亮的结构后果：$t=0$ 时 $dt=0$，故 $\psi_{i,dt}=\mathcal I$，初始条件精确且自动成立**——不需要初值损失项。

**第二步：两种导数各走各的路。** 时间导数 $\partial\hat u/\partial t$ 用**自动微分**，这只有在 $t$ 经 $dt=t/K$ 显式进入时才可能；空间导数用**有限差分**，因为 $\bm x$ 不是网络输入（网络把固定网格上的函数映到函数）。

**第三步：物理信息损失。** 残差 $f:=u_t-\mathcal Fu$，

$$
\mathfrak L(\bm\theta)=\lambda_r\mathfrak L_r+\lambda_b\mathfrak L_b,
\qquad
\mathfrak L_r=\frac{1}{N_uN_r}\sum_{n=1}^{N_u}\sum_{j=1}^{N_r}\bigl|f^{(n)}_{r,j}\bigr|^2,
\qquad
\mathfrak L_b=\frac{1}{N_uN_b}\sum_{n=1}^{N_u}\sum_{j=1}^{N_b}\bigl|\mathcal Bu^{(n)}_{b,j}\bigr|^2 .
$$

空间配点取自固定的等距点集 $\{\bm x_i\}_{i=1}^{N_x}$，时间点从 $(0,T]$ 随机抽。边界条件若由结构保证（如周期），$\mathfrak L_b$ 就去掉。**完全没有数据项。**

**第四步：长时间推断靠迭代。** 直接令 $t=t^\star>T$ 精度很差，所以对 $t^\star\in(mT,(m+1)T]$ 改用复合：$\hat u(T)=\Phi_{\bm\theta_T}(u_0,T)$，再 $\hat u(2T)=\Phi_{\bm\theta_T}(\hat u(T),T)$，如此 $m+1$ 次。

**第五步：用残差自诊断、决定何时重训。** 因为没有数据，论文拿**物理残差自身**当有效性监视器：

$$
\mathcal R_t=\frac{1}{|\mathfrak I_0|N_x}\sum_{n}\sum_{j}
\Bigl|\partial_t\hat u^{(n)}(t,\bm x_j)-\mathcal F\hat u^{(n)}(t,\bm x_j)\Bigr|^2
+\Bigl|\mathcal B\hat u^{(n)}(t,\bm x_j)\Bigr|^2,
$$

只要 $\mathcal R_t\le\epsilon\cdot\mathcal R_T$（$\epsilon\ge1$ 由用户设定）就接受预测。若在 $t=mT$ 处被违反，就从 $\{\hat u^{(n)}((m-1)T,\cdot)\}$ 中随机抽 $N_{\rm add}$ 个函数追加到初值训练集 $\mathfrak I$ 里并重训。动机说得很明确：振荡出现在「$\mathfrak I$ 中训练初值函数张成的子空间不足以刻画已显著演化的状态」之时。

### 定理

**这是本页少数几篇有真定理的论文之一。**

**线性稳定性分析。** 在试验方程 $u_t=\mathcal Lu+\lambda u$ 上取 $\mathcal Lu=-q_{\bm\theta}u$、$q_{\bm\theta}>0$ 为实数、$\lambda$ 为复数，一块给出

$$
u_{i+1}=e^{\lambda dt}\Bigl(1-q_{\bm\theta}dt+\tfrac12q_{\bm\theta}^2dt^2\Bigr)u_i,
$$

增长因子

$$
\xi=\Bigl(1-q_{\bm\theta}dt+\tfrac12q_{\bm\theta}^2dt^2\Bigr)
e^{\mathrm{Re}(\lambda dt)+\mathrm i\,\mathrm{Im}(\lambda dt)},
$$

稳定边界 $|\xi|=1$ 等价于

$$
\mathrm{Re}(\lambda dt)=\ln\!\left(\frac{1}{1-q_{\bm\theta}dt+\tfrac12q_{\bm\theta}^2dt^2}\right).
$$

由此得到的结论是：$0<q_{\bm\theta}dt<1$ 时稳定域随 $q_{\bm\theta}dt$ 增大而扩大，所有这些稳定域的交是左半平面；而 $q_{\bm\theta}dt\ge2$ 时二阶展开对任何 $\lambda>0$ 都**无条件不稳定**。若 $\max_u\partial(\mathcal Nu)/\partial u\le0$（即 $\mathrm{Re}\,\lambda\le0$）则 PI-DOSnet 无条件稳定；否则需要

$$
dt\le\frac{1}{\mathrm{Re}\lambda}\ln\Bigl(1\big/\bigl(1-q_{\bm\theta}dt+\tfrac12q_{\bm\theta}^2dt^2\bigr)\Bigr),
$$

论文称之为「与 CFL 条件类似」，等价地写作 $q_{\bm\theta}dt\le C$。**实践上最关键的一条是：对一个学出来的 $\mathcal L_{\bm\theta}$，其特征值 $q_{\bm\theta}$「对网格尺寸 $h$ 的敏感性远低于」有限差分离散**，这正是它能在细网格上取大 $dt$ 的原因。

**定理 3.1。** 设 $\mathcal X$ 是 Hilbert 空间、$\mathcal L:\mathcal X\to\mathcal X$ 线性、紧、自伴。则对任何满足 $\|u\|_{\mathcal X}\le C_0$ 的 $u\in\mathcal X$，

$$
\Bigl\|e^{dt\,\mathcal L}u-\bigl(\mathcal I+dt\,\mathcal L+\tfrac{dt^2}{2}\mathcal L^2\bigr)u\Bigr\|_{\mathcal X}\le C\,dt^3,
$$

$C$ 依赖 $\|u\|_{\mathcal X}$ 与 $\|\mathcal L\|_{\mathcal X\to\mathcal X}$。证明在 $\mathcal L$ 的特征基上做谱分解。

**误差分解。** 单块的总误差分成 $\epsilon_1+\epsilon_2$，$\epsilon_1$ 是 Taylor 截断、$\epsilon_2$ 是 $\mathcal L$ 与 $\mathcal L_{\bm\theta}$ 之差。定义 3.1 引入平移等变性 $\mathcal L(f(x-a))=(\mathcal Lf)(x-a)$；对一维 $\mathcal L=\Delta$ 配核 $\tfrac{1}{h^2}[1,-2,1]$，论文给出

$$
\bigl|[e^{dt\Delta}u](\bm x_j)-[\tilde{\mathcal L}_{\bm\theta}u](\bm x_j)\bigr|
\le C_1\,dt^3+C_2\,dt\,h^2 .
$$

定理后的注记说：$dt^3$ 这一项与网络参数无关且随 $dt$ 变小而变小；$dt$ 较大时（例如 $dt=1/4$）可以靠**加大卷积核**来缓解误差。论文还引 DOSnet 指出算子分裂块的逼近误差本身是 $dt^3$ 阶，因此不再重推。

> [!warning] 摘要与正文的一处强弱差别
> 摘要称对 Allen-Cahn 方程「即使用大时间步也得到能量稳定的解」，但正文说预测能量呈现「近乎单调的下降」，且「能量耗散律并未被严格施加」。因此「能量稳定」是一个经验观察而非已证性质；作者也说设计严格满足能量耗散的网络「仍是有挑战的任务，留待将来」。

### 数值实验

度量是给定时刻的相对 $L^2$ 误差与 $[0,\tilde t]$ 上的累计相对 $L^2$ 误差；实现用 PyTorch 与 Adam。

**例 1：一维对流。** $u_t+\beta u_x=0$ 于 $(-\pi,\pi)$、周期边界；初值 $u_0(x)=\sum_{i=1}^{10}(c_i\sin ix+q_i\cos ix)$、$c_i,q_i\sim\mathcal N(0,1)$；300 个等距空间点；4 个块、通道 1-4-1、核尺寸 61、$dt=0.075$；以 $T=0.3$ 训 20000 轮。一次前向就能预测 $2T,\dots,10T$；$10T$ 处绝对误差低于 0.05，时空平均相对 $L^2$ 误差 $3.149\times10^{-3}$。

**Taylor 阶数消融**（一阶对二阶）：

| $t$ | 一阶                 | 二阶                 |
| --- | -------------------- | -------------------- |
| 3   | $1.640\times10^{-2}$ | $5.318\times10^{-3}$ |
| 6   | $3.283\times10^{-2}$ | $1.062\times10^{-2}$ |
| 9   | $4.927\times10^{-2}$ | $1.592\times10^{-2}$ |
| 15  | $8.220\times10^{-2}$ | $2.650\times10^{-2}$ |
| 21  | $1.151\times10^{-1}$ | $3.704\times10^{-2}$ |
| 30  | $2.350\times10^{-1}$ | $5.278\times10^{-2}$ |

大致是 3 到 4 倍的差距，且随时间拉大。**这就是二阶展开的实验依据。**

**例 2：扩散反应。** $u_t=Du_{xx}+ku^2$ 于 $(0,1)$、$D=k=0.001$、齐次 Dirichlet、$T_{\rm end}=50$，训练窗口 $T\in\{0.5,1,2\}$。4 个块、每块 3 个卷积层、通道 1-8-8-1、核 21，**总共只有 6720 个可训练参数**，101 个空间点。

**例 3：Allen-Cahn。** $u_t-0.0004\Delta u+(u^3-u)=0$ 于 $(-1,1)^d$、周期边界、$d=1,2$；初值取随机 Fourier 型（一维 7 个模态，二维一个 $8\times8$ 的双重和再乘 0.02）；参考解用 Chebfun。一维：$T=1$、$T_{\rm end}=10$、$200\times20$ 个配点、4 个块、通道 1-4-1、核 15、8000 轮、1000 个初值函数。推断在 $t=3$ 处终止；追加 100 个 $t=2$ 处的解并重训后可达 $t=6$；再来一轮覆盖 $[0,10T]$。两次重训后均方误差降到约 $4\times10^{-4}$，绝对误差小于 0.04。能量 $E=\tfrac{0.0004}{2}\int u_x^2+\int\tfrac{(u^2-1)^2}{4}$ 近乎单调下降并贴合 Chebfun 参考。

**网格加密与稳定性对照——这是全篇最强的结果：**

| $N_x$ | PI-DOSnet（学出的卷积） | 二阶中心差分核 | 四阶中心差分核 |
| ----- | ----------------------- | -------------- | -------------- |
| 200   | 4 块，0.07 s            | 8 块，0.08 s   | 11 块，0.10 s  |
| 400   | 4 块，0.07 s            | 32 块，0.26 s  | 42 块，0.34 s  |
| 800   | 4 块，0.07 s            | 128 块，1.00 s | 170 块，1.35 s |

PI-DOSnet 在三个网格上都只用 **4 块**，$t=1$ 处的 $L^2$ 误差为 $1.81\times10^{-3}$、$1.97\times10^{-3}$、$1.11\times10^{-3}$，$t=10$ 处为 $9.86\times10^{-3}$、$1.19\times10^{-2}$、$9.57\times10^{-3}$。机制由谱半径证实：$-\mathcal L_{\bm\theta}$ 的谱半径随 $N_x$「只缓慢增长」，而二阶中心差分 $-\mathcal L$ 的谱半径「二次增长」，从而强加了 CFL 限制。二维情形：$T=1$、$T_{\rm end}=10$、20 个时刻、$200\times200$ 网格、4 个块、4 个通道、$15\times15$ 核、400 个初值函数；推断同样停在 $t=3$，追加 40 个 $t=2$ 处的预测、再追加 $t=5$ 处的预测。

**例 4：Gross-Pitaevskii 方程。** $\mathrm i\partial_t\psi=[-\tfrac{\partial_{xx}}{2m}+V(x)+g(x)|\psi|^2]\psi$（复值）。孤立波由

$$
F=\frac{1}{N^2}\Bigl|\int\psi^*(T_{\rm end},x)\psi(0,x)\,\mathrm dx\Bigr|^2,
\qquad N=\int|\psi(0,x)|^2\mathrm dx
$$

最小化 $|F-1|$ 来识别。$1024\times10$ 个配点，训练区间 $[0,0.1]$，$\tilde T_{\rm end}=1$，8 个块、深度 2、核 21、15000 轮、1000 个初值函数；参考解用时间分裂谱方法。在 100 个随机测试初值上，$\tilde T=0.1$ 处平均相对 $L^2$ 误差 $5.325\times10^{-3}$，$\tilde T_{\rm end}=1$ 处 $4.087\times10^{-2}$。

这些实验建立的是：**把分裂结构写进架构、再让 $t$ 显式进入输入，就同时买到了三样东西**——无数据训练、任意时刻取值、以及不随网格加密而恶化的时间步。它们没有建立的是能量耗散的严格性（见上面的警告），也没有覆盖论文自己列出的三个开放问题：复杂几何、高维问题、低正则解。

### 与其他论文的关系

小结中列为开放问题的「低正则解」，恰恰是编号 90 的目标——两文是同一组研究计划的互补两半，见[[computational-mathematics/paper-notes/scientific-machine-learning/variational-and-basis-networks|变分与基网络一页]]。它与编号 94（配高斯波包的物理信息 DeepONet）和编号 95（UQ-SONet）共享物理信息算子学习的设定：三者都学从初值到解的映射，其中 94 与 103 无成对数据、95 带不确定性。架构上的差别是 PI-DOSnet 的归纳偏置**就是分裂格式本身**，以 $e^{dt\mathcal N}$ 当激活函数，而不是分支/主干的因子分解。

Jizu Huang 与 Tao Zhou 也共同署名编号 94；PI-DOSnet 用的算子分裂，正是编号 94 中因 $O(\Delta t^2/\varepsilon)$ 失效而促成波包重述的那件经典工具。把耦合系统顺序解耦成更便宜的子问题，在时间上是编号 102 空间解耦的对应物。而残差触发的重训循环，是编号 90 残差触发的区域加密与编号 101 离散 Fourier 变换触发的词元注入的自诊断版本：三者都用当前逼近上的一个可计算指标来决定在哪里增加容量。

## 105：一个共享的摘要网络，把滤波与平滑绑在一起

### 直觉

在状态空间模型里，滤波要的是 $p(u_t\mid y_{1:t})$，平滑要的是 $p(u_{1:t}\mid y_{1:t})$。两者的条件变量都是**整段观测历史**，其维数随 $t$ 增长，所以固定输入的正规化流用不上。办法是先用一个循环网络把历史压成定长摘要 $s_t$，再让流以 $s_t$ 为条件。

到这一步为止都算常规。这篇论文真正的动作是：**让滤波流与平滑流共用同一个摘要网络。** 理由有两层。第一层是代数的：平滑有一个精确的因果分解，其中后向核的条件是 $y_{1:k}$ 而不是 $y_{1:t}$，因此同一个 $s_k$ 可以直接复用；在一个特定的权重取值上，两项损失重新组合成「逐时刻滤波似然」加「整条轨迹的平滑似然」，即这**一个**损失同时是两个任务的极大似然。第二层是充分性的：若存在一个充分摘要，它同时是两个目标各自的最优摘要。

第二层只是充分条件，不是必然。论文自己的引理指出两个目标的最优摘要一般**不重合**，所以共享不是免费的。而它的消融实验给出了另一个方向的答案：不共享时，两个流单独看都几乎不受影响，但后向**迭代采样**整体崩掉。也就是说共享在这里不是有益，而是必需。

### 问题设定

状态空间模型

$$
u_{t}=f(u_{t-1},\epsilon_{u,t}),\qquad y_{t}=h(u_{t},\epsilon_{y,t}),
$$

$u_t\in\mathbb R^{n_u}$、$y_t\in\mathbb R^{n_y}$，$\epsilon_{u,t}$ 与 $\epsilon_{y,t}$ 相互独立且沿时间独立同分布，诱导 Markov 分解

$$
p(u_{1:t},y_{1:t})=p(u_{1})p(y_{1}\mid u_{1})\prod_{k=2}^{t}p(u_{k}\mid u_{k-1})p(y_{k}\mid u_{k}).
$$

训练数据是 $N$ 条模拟轨迹 $\{u^i_{1:T},y^i_{1:T}\}_{i=1}^{N}$。

论文对已有工作的指控分三条：高斯滤波器（EKF、EnKF）依赖「只到前两阶矩的矩封闭拟设，对一般非线性系统有系统性偏差」；序贯 Monte Carlo 与粒子滤波「随状态维数与时间跨度增大通常出现严重的权退化」，平滑时更糟；近年的深度生成式滤波器「通常把滤波与平滑分开处理，且往往依赖代价高昂的逐实例优化」，因而不适合在线与大规模部署。**论文对自己的约束比通常更强**：只假设有一个**模拟器**，$f$、$h$ 与噪声分布的函数形式都不假设已知，转移密度与观测密度也不必可算。

条件正规化流的定义为：对条件变量 $s$ 与参照 $Z\sim p_Z=\mathcal N(0,I)$，

$$
z=f_{\theta}(u;s),\qquad
p_{\theta}(u\mid s)=p_{Z}\bigl(f_{\theta}(u;s)\bigr)\bigl|\det\nabla_{u}f_{\theta}(u;s)\bigr|,
$$

$f_\theta=f_{[L]}\circ\cdots\circ f_{[1]}$，行列式按层分解，按极大似然训练（等价于极小化对 $s$ 平均的条件 KL）。

### 推导

**第一步：把历史压成定长摘要。** 用多层 LSTM，单层单元写全为

$$
i_{t}=\sigma(W_{yi}y_{t}+W_{hi}h_{t-1}+b_{i}),\qquad
f_{t}=\sigma(W_{yf}y_{t}+W_{hf}h_{t-1}+b_{f}),
$$

$$
o_{t}=\sigma(W_{yo}y_{t}+W_{ho}h_{t-1}+b_{o}),\qquad
g_{t}=\tanh(W_{yg}y_{t}+W_{hg}h_{t-1}+b_{g}),
$$

$$
c_{t}=f_{t}\odot c_{t-1}+i_{t}\odot g_{t},\qquad
h_{t}=o_{t}\odot\tanh(c_{t}),\qquad (h_0,c_0)=(\bm 0,\bm 0),
$$

叠 $L$ 层（$z^{(1)}_t:=y_t$、$z^{(\ell)}_t:=h^{(\ell-1)}_t$），摘要由顶层经仿射映射读出：

$$
s_{t}=W_{s}h_{t}^{(L)}+b_{s},\qquad W_{s}\in\mathbb R^{d_{s}\times d_{h}} .
$$

注记 2.1 给出选 LSTM 而非普通 RNN 的理由：把加性的细胞更新展开，

$$
c_{t}=\Bigl(\prod_{j=1}^{t}f_{j}\Bigr)\odot c_{0}
+\sum_{s=1}^{t}\Bigl(\prod_{j=s+1}^{t}f_{j}\Bigr)\odot(i_{s}\odot g_{s}),
$$

故 $\partial c_t/\partial c_{t-k}=\prod_{j=t-k+1}^{t}f_j$，「沿细胞状态的梯度由门的乘积支配，而不是由同一个线性映射与非线性的反复作用支配」。

**第二步：两个流。** 前向流直接以摘要为条件：

$$
p(u_{t}\mid y_{1:t})\approx p_{\theta_1,\psi}(u_{t}\mid s_{t}),
\qquad s_t=\mathrm{Enc}(y_{1:t};\psi)\in\mathbb R^{h} .
$$

平滑一侧用的是一个**精确分解**

$$
p(u_{1:t}\mid y_{1:t})=p(u_{t}\mid y_{1:t})\prod_{k=1}^{t-1}p(u_{k}\mid u_{k+1},y_{1:k}),
$$

注意第二个因子的条件是 $y_{1:k}$ 而**不是** $y_{1:t}$——正是这一点让递推保持因果，于是同一个摘要 $s_k$ 可以直接复用。第二个流学

$$
p(u_{t}\mid u_{t+1},y_{1:t})\approx p_{\theta_2,\psi}(u_{t}\mid u_{t+1},s_{t}).
$$

**两个流共用同一个摘要网络 $\psi$，这是论文冠名的贡献。**

**第三步：联合目标与那个代数恒等式。**

$$
\min_{\theta_{1},\theta_{2},\psi}\
-\frac1{NT}\sum_{i=1}^N\sum_{t=1}^T\log p_{\theta_1,\psi}(u_t^i\mid s_t^i)
-\frac{\lambda}{N(T-1)}\sum_{i=1}^N\sum_{t=1}^{T-1}\log p_{\theta_2,\psi}(u_t^i\mid u_{t+1}^i,s_t^i).
$$

在 $\lambda=(T-1)/T$ 处两项重新组合。定义诱导的轨迹密度

$$
p^{\mathrm{smoothing}}_{\theta_{1},\theta_{2},\psi}(u_{1:T}\mid s_{1:T})
:=p_{\theta_{1},\psi}(u_{T}\mid s_{T})\prod_{t=1}^{T-1}p_{\theta_{2},\psi}(u_{t}\mid u_{t+1},s_{t}),
$$

则损失恰好成为

$$
\mathcal L=-\frac{1}{NT}\sum_{i=1}^{N}\sum_{t=1}^{T-1}\log p_{\theta_{1},\psi}(u^{i}_{t}\mid s^{i}_{t})
-\frac{1}{NT}\sum_{i=1}^{N}\log p^{\mathrm{smoothing}}_{\theta_{1},\theta_{2},\psi}(u^{i}_{1:T}\mid s^{i}_{1:T}),
$$

即「中间时刻滤波的边缘项」加「诱导轨迹模型的路径项」。因为两者都以同一组 $\{s_t\}$ 为条件，「学到的后向转移必须在一个共同的观测历史表示下与终端时刻的滤波分布保持一致」。**这就是摘要里说的隐式一致性正则化，而且它是一个精确改写而非启发式。**

**第四步：条件 KRnet 的细节。** $z=T_\Theta(u;c)$，$\log p_\Theta(u\mid c)=\log p_Z(T_\Theta(u;c))+\log|\det\nabla_uT_\Theta(u;c)|$。模块有三类：条件尺度-偏置层 $T^{\mathrm s}_{\psi}(u;c)=\exp(\eta_{\psi}(c))\odot u+\xi_{\psi}(c)$，对数行列式为 $\sum_j\eta_{\psi,j}(c)$；条件仿射耦合层，按 $u=[u^{(1)},u^{(2)}]^\top$、$u^{(1)}\in\mathbb R^{k}$、$k=\lfloor n_u/2\rfloor$ 划分，

$$
\tilde u^{(1)}=u^{(1)},
\qquad
\tilde u^{(2)}=\bigl(\bm 1+\alpha\tanh s_{\omega}(u^{(1)},c)\bigr)\odot u^{(2)}
+\gamma\odot\tanh t_{\omega}(u^{(1)},c),
$$

其中 **$\alpha=0.6$ 固定**、$\gamma$ 为正的可学向量；由 $\tanh\in(-1,1)$ 知尺度落在 $(1-\alpha,1+\alpha)$ 内，可逆性由此保证，对数行列式为 $\sum_j\log(1+\alpha\tanh s_{\omega,j})$。耦合网络是**随机 Fourier 特征网络**：以 $h_0=[u^{(1)},c]^\top$ 起步，

$$
h_{1}=\begin{bmatrix}\sin(e^{-\sigma}Fh_{0}+b_{0})\\ \cos(e^{-\sigma}Fh_{0}+b_{0})\\ h_{0}\end{bmatrix},
\qquad
h_{i}=\mathrm{SiLU}(W_{i-1}h_{i-1}+b_{i-1}),
$$

$(F,b_0)$ 为**固定的随机特征**，带宽型标量 $\sigma$ 与 $\{W_i,b_i\}$ 可训练。整体复合为 $T_{\Theta}=\Pi_{K}\circ T^{\mathrm{coup}}_{\omega_{K}}\circ\cdots\circ\Pi_{1}\circ T^{\mathrm{coup}}_{\omega_{1}}\circ T^{\mathrm{s}}_{\psi}$，$\Pi_k$ 是行列式为 1 的固定坐标置换。

**第五步：基于流的粒子滤波变体。** 不用完全摊销的边缘，而用**完全适配提议**的分解

$$
p(u_{k},y_{k}\mid u_{k-1})=p(y_{k}\mid u_{k-1})\,p(u_{k}\mid y_{k},u_{k-1}),
$$

把两个因子都学成条件流：$p_{\theta_3}(y_k\mid u_{k-1})$ 与 $p_{\theta_4}(u_k\mid y_k,u_{k-1})$，各按自己在模拟三元组上的负对数似然训练。一步滤波为：预测加权 $\tilde\omega^{(j)}_{k-1}\propto\omega^{(j)}_{k-1}p_{\theta_3}(y_k\mid u^{(j)}_{k-1})$；重采样 $a^{(j)}\sim\mathrm{Cat}(\tilde\omega^{(1:N)}_{k-1})$；传播 $u^{(j)}_{k}\sim p_{\theta_4}(\cdot\mid y_{k},u^{(a^{(j)})}_{k-1})$。注记 3.2 与 bootstrap 型替代（$p_{\theta_5}(u_k\mid u_{k-1})$、$p_{\theta_6}(y_k\mid u_k)$）对比：适配提议看得到 $y_k$，因此「预期能生成信息量更大的粒子」。

**第六步：有效样本量诊断（仅在真转移密度与似然可得时可用）。** 比较

$$
p(u_t,y_t,u_{t-1}\mid y_{1:t-1})=p(u_t\mid u_{t-1})p(y_t\mid u_t)p(u_{t-1}\mid y_{1:t-1})
$$

与学到的 $q=p_{\theta_4}(u_t\mid y_t,u_{t-1})p_{\theta_3}(y_t\mid u_{t-1})p(u_{t-1}\mid y_{1:t-1})$，权

$$
\omega_{i}=\frac{p(u^{(i)}_{t}\mid u^{(i)}_{t-1})p(y^{(i)}_{t}\mid u^{(i)}_{t})}
{p_{\theta_{4}}(u^{(i)}_{t}\mid y^{(i)}_{t},u^{(i)}_{t-1})p_{\theta_{3}}(y^{(i)}_{t}\mid u^{(i)}_{t-1})},
\qquad
\mathrm{ESS}=\frac{(\sum_i\omega_i)^2}{\sum_i\omega_i^2},
\qquad
\mathrm{RESS}=\frac{\mathrm{ESS}}{N},
$$

并由 $\chi^{2}(p\|q)=\int p^2/q\,\mathrm dx-1\approx N/\mathrm{ESS}-1$ 与卡方散度相连。注记 3.3 给出 $q$ 的祖先采样顺序，并说明取 $N=10^{6}$ 才有稳定的 Monte Carlo 估计。

**推断。** 滤波（算法 2）只需算 $s_t$ 再从 $p_{\theta_1,\psi}(\cdot\mid s_t)$ 抽样，是**摊销**的，测试时不做逐实例优化。平滑（算法 3）先算 $s_1,\dots,s_t$、从终端滤波分布抽 $\{u^j_t\}$，再对 $k=t-1,\dots,1$ 依次抽 $u^j_k\sim p_{\theta_2,\psi}(\cdot\mid u^j_{k+1},s_k)$ 并前接到路径上；在任意 $k$ 处提前终止就得到 $p(u_k\mid y_{1:t})$ 的边缘平滑样本。

### 定理

**理论只有附录 B 中的一条引理与一条命题，都是关于共享摘要统计量的。** 记 $X:=u_k$、$Y:=y_{1:k}$、$Z:=u_{k+1}$，$S=f(Y)$ 是任意可测摘要。两个目标 $\min\mathbb E[-\log p_{\theta_1}(X\mid f(Y))]$ 与 $\min\mathbb E[-\log p_{\theta_2}(X\mid Z,f(Y))]$ 在对流参数取最优时分别达到条件熵 $H(X\mid S)$ 与 $H(X\mid Z,S)$。记 $S^*_F=\arg\min_{S=f(Y)}H(X\mid S)$、$S^*_S=\arg\min_{S=f(Y)}H(X\mid Z,S)$。

**引理 B.1。** $I(X,S^{*}_{F})\ge I(X,S^{*}_{S})$ 且 $I(X;Z\mid S^{*}_{S})\ge I(X;Z\mid S^{*}_{F})$。（证明用 $H(X\mid S)=H(X\mid Z,S)+I(X;Z\mid S)$ 与两处极小性。）读法：**滤波与平滑各自的最优摘要不必重合——所以共享不是免费的。**

**定义 B.1 与命题 B.1。** 称摘要 $S^{\dagger}=f(Y)$ **充分**，若 $p(X\mid Y)=p(X\mid S^{\dagger})$ 几乎必然成立。若充分摘要存在，则 $S^{\dagger}=S^{*}_{F}=S^{*}_{S}$，即一个共享的摘要网络对两个目标同时最优。

**对 FLUID 本身没有收敛、相合、误差或稳定性定理。** 结论一节把「建立严格的理论保证」列为将来的工作。

### 数值实验

共享配置：两个流都是一个条件尺度-偏置层加 **6 个条件仿射耦合层**，每个耦合网络是随机 Fourier 特征嵌入接一个**深度 6、宽度 64**的多层感知机；摘要网是 **4 层 LSTM**（线性对流扩散那一例减为 **1 层**）；$\dim(s_t)=3\dim(y_t)$，随机波动率一例取 $5\dim(y_t)$；Adam、初始学习率 $0.001$；PyTorch，单张 **NVIDIA V100（32 GB）**。度量为 KL（只在精确后验可得时）、RMSE、MMD（高斯核 $\sigma=2$，对真状态处的 delta 测度）与 CRPS。

**例 1：一维线性对流扩散。** $\partial_tu=a\partial_xu+\kappa\partial_{xx}u$ 于 $[0,1]$、周期边界、$u(0,x)=\sin(2\pi x)$，离散成 $u_k=Mu_{k-1}+\epsilon_{u,k}$、$y_k=Hu_k+\epsilon_{y,k}$、$u_0\sim\mathcal N(\mu,\sigma^2I)$——一个**滤波与平滑律都有闭式的线性高斯模型**，这正是 KL 那一列能存在的原因。

情形 1（$a=-1$、$\kappa=0$ 的纯对流，迎风 $M_{\delta t}=I-\nu A$、$\nu=\delta t/\Delta x$，$Q=qI_n$，$H$ 隔点取样即 $n_y=n/2$，$R=rI$）：$\Delta t=0.05$、$q=0.01$、$r=0.1$、$\sigma=0.05$，随 $n$ 增大而**保持固定**，因此问题确实越来越难；$n=10,20,30,40,50$；$N_{\rm train}=2000$ 条长度 $T_{\rm train}=500$ 的轨迹，$N_{\rm test}=200$。

| $n$ | 滤波 KL | 后向核 KL | 滤波 RMSE | 平滑 RMSE | MMD    | CRPS   |
| --- | ------- | --------- | --------- | --------- | ------ | ------ |
| 10  | 0.0191  | 0.0146    | 0.1457    | 0.1284    | 0.0514 | 0.0822 |
| 50  | 0.0597  | 0.0534    | 0.1228    | 0.1134    | 0.1710 | 0.0693 |

注意滤波 RMSE 与 CRPS 随 $n$ 增大反而**下降**，只有 MMD 上升；论文自己标出了 MMD 的增长，并论证其余度量保持稳定。$n=50$ 的算例被跑到 $T=1000$（物理时间 $t=50$），即**训练时长的两倍**，一致性在虚线之后仍然维持。

情形 2（$a=1$、$\kappa=0.01$ 的离散一致区，Lax-Wendroff 型 $M_{\delta t}$，$Q_{\delta t}=\frac{\delta t}{n}I_n$，块平均观测 $y^i_k=\frac{1}{|\mathcal I_i|}\sum_{j\in\mathcal I_i}u^j_k+\epsilon^i_{y,k}$、$n_{\mathcal I}=8$ 固定）：$\Delta t=0.01$、$r=0.01$、$\sigma=0.05/n$；$n=16,32,48,64$。滤波 KL 停在一个窄带里（$0.1355$、$0.1363$、$0.1135$、$0.1268$），而 RMSE、MMD、CRPS 都随 $n$ **下降**（RMSE 从 $0.0558$ 到 $0.0214$）——这正是网格一致行为该有的样子。

**例 2：两因子随机波动率模型**（全篇唯一涉及真实数据的一例）。$\bm u_t=\bm\alpha+A(\bm u_{t-1}-\bm\alpha)+D_\sigma\epsilon_{u,t}$、$\bm y_t=\beta\exp(\tfrac12\bm u_t)\odot\epsilon_{y,t}$、$\bm u_0\sim\mathcal N(0,\operatorname{diag}(\tau_1^2,\tau_2^2))$，参数 $\bm\alpha=0$、$\gamma_1=\gamma_2=0.97$、$\sigma_1=\sigma_2=0.3$、$\beta=0.835$（取自文献中标普 500 的 MAP 估计），$\tau_i^2=\sigma_i^2/(1-\gamma_i^2)$。$N=2000$ 条长度 $T_{\rm train}=1000$ 的轨迹，$N_{\rm test}=200$。

| 方法             | RMSE   | MMD    | CRPS   |
| ---------------- | ------ | ------ | ------ |
| FLUID            | 0.6117 | 0.1571 | 0.3435 |
| 基于流的粒子滤波 | 0.6163 | 0.1591 | 0.3462 |
| FBF              | 0.7481 | 0.2163 | 0.4188 |

后向核为 $0.2746/0.0363/0.1551$、平滑为 $0.4805/0.1038/0.2710$，相对滤波都有大幅改善。测试轨迹跑到 $T=2000$，是训练时长的两倍。**真实数据应用**：把训练好的粒子滤波用在 2018 年 12 月 31 日到 2022 年 12 月 29 日的标普 500 日收益上（二元序列由处理后的标普 500 收益加一条合成的单因子轨迹拼成），恢复出的隐波动率在 $k\approx300$ 附近有一个峰、90% 可信区间同时变宽（2020 年初的疫情冲击），并从 $k\approx800$ 起持续偏高（2022 年的动荡）。这是一个**解释性**而非定量的结论。

**例 3：随机 Burgers 方程。** $du+(u\partial_xu-\nu\partial^2_{xx}u)dt=\sigma\,dW(t)$，$x\in[-1,1]$、$t\in[0,1]$、$u(0,x)=-\sin(\pi x)$、$u(t,\pm1)=0$，$\sigma=1.0$、$\nu=0.05$（附录 E.1 另做 $\nu=0.01$）。混合有限差分在 $201\times50$ 网格上，状态是 50 个空间值、观测隔点取（$n_y=25$）并加 $\mathcal N(0,r^2I)$ 噪声，$\Delta t=0.005$；$N=3000$ 条长度 $T_{\rm train}=200$ 的训练轨迹，$N_{\rm test}=200$。

| $r^2$ | FLUID 滤波 RMSE | FBF 滤波 RMSE |
| ----- | --------------- | ------------- |
| 0.01  | 0.0751          | 0.0752        |
| 0.04  | 0.0898          | 0.0917        |
| 0.09  | 0.1003          | 0.1015        |
| 0.16  | 0.1084          | 0.1113        |
| 0.25  | 0.1149          | 0.1179        |

论文自己的读法是两者「极其接近」而 FLUID 有「微小但一致的优势」，且两者都接近这个问题的精度上限。后向核比滤波锐利得多（$r^2=0.25$ 处 RMSE $0.0584$ 对 $0.1149$）。附录 E.1 在 $\nu=0.01$ 上重做同一扫描，序不变（例如 $r^2=0.25$ 处 FLUID $0.1546$ 对 FBF $0.1640$）。

**例 4：Lorenz-96 与两尺度 Lorenz。** 参数 $J=32$、$F\in\{5,8,16\}$、$h=1$、$b=10$、$c=4$、$\sigma_u=0.1$、$\sigma_v=0.01$；观测是**立方**的，$y_{k,i}=u_{k,i}^{3}+\epsilon_{y,i}$、$\epsilon_{y,i}\sim\mathcal N(0,1)$，强非线性。

单尺度（$c=0$、$\sigma_u=1$、$F=8$、全部指标可观测、$u_{0,j}=\sin(2\pi j/n)$、$\Delta t=0.05$、$K=10,\dots,50$）：

| $K$ | FLUID 滤波 RMSE | FBF 滤波 RMSE |
| --- | --------------- | ------------- |
| 10  | 0.1632          | 0.2044        |
| 20  | 0.1945          | 0.2439        |
| 30  | 0.2081          | 0.2604        |
| 40  | 0.2255          | 0.2742        |
| 50  | 0.2605          | 0.3106        |

FLUID 在每个维数上都更好，MMD 差距还在拉大（$K=50$ 处 $0.5508$ 对 $0.6918$）。

两尺度（$K=16$、$J=32$、只观测奇数指标、$u_{0,i}=F+\sigma_u\epsilon_u$、$v_{0,j}=\sigma_v\epsilon_v$）：此处**把快变量边缘化之后，慢变量的有效转移密度不可算**，因而模型违反 FBF 所要求的标准形式。

| $F$ | FLUID 滤波 RMSE | FBF 滤波 RMSE | 平滑 RMSE |
| --- | --------------- | ------------- | --------- |
| 5   | 0.4544          | 0.5817        | 0.4063    |
| 8   | 0.4397          | 0.8301        | 0.3795    |
| 16  | 0.5218          | **3.7536**    | 0.4837    |

FBF 在 $F=16$ 处的失败（CRPS $1.3578$ 对 $0.2350$）是全篇最强的对比证据。

### 共享不是有益，而是必需

最值得记录的是消融实验。状态维数 $K=20$ 时，共享与独立摘要的对比是：

| 指标        | 共享摘要 | 独立摘要    |
| ----------- | -------- | ----------- |
| 滤波 RMSE   | 0.1945   | 0.1958      |
| 后向核 RMSE | 0.1240   | 0.1596      |
| 平滑 RMSE   | 0.1525   | **60.3744** |
| 平滑 CRPS   | 0.0742   | **33.2893** |

两个流**单独**看几乎不受影响，而后向**迭代采样**整体崩掉。论文的解释是独立摘要下误差在平滑递推中迅速累积，$K\ge20$ 时结果已不可用。这正是上面那个损失重组恒等式的实验对应物，也是本文的实际贡献所在。

> [!warning] 两处正文与表格不一致
> 论文正文称平滑一致优于滤波，但 Table 9 在 $K=50$ 处给出平滑 RMSE 0.5423 而滤波 0.2605，Burgers 算例的 Table 6 在 $r^2=0.25$ 处有同样的反转（平滑 0.1185 对滤波 0.1149）。本页按表格数值记录。$K=50$ 处那次明显退化，正文并未评论。

> [!warning] 投稿状态未核实
> arXiv 记录上没有期刊引用字段，因此「投稿 CMAME」这一条**未由预印本本身核实**。

### 与其他论文的关系

技术上最近的是编号 107：两者都用**条件 KRnet**、耦合网络都用**随机 Fourier 特征**，都学一族由条件变量索引的条件密度，也都引同一条时间正规化流的谱系。差别在于什么在给流做条件、以及训练信号从哪来——编号 107 以初始状态与时间为条件、按 **PDE 残差**训练（物理信息，无样本），编号 105 以观测历史的 LSTM 摘要为条件、按模拟轨迹上的**极大似然**训练（基于模拟，无 PDE）。两者恰好框出了拟合条件正规化流的两条路。

Xiaodong Feng、Xiaoliang Wan、Tao Zhou 三人还共同署名编号 101 与 98。耦合网络里的随机 Fourier 特征嵌入，正是编号 101 用来搭多尺度基库的同一个原语——在编号 105 里它是一个带单个**可训练**带宽标量 $\sigma$ 的固定特征化器，在编号 101 里则是一个刻意多尺度、带自适应放大的基库。

它与编号 98、95 同属该组的不确定性量化一支，但概率工具不同：那两篇用带高斯解码器的隐变量模型与证据下界/信息瓶颈目标（给出的是似然的**界**），编号 105 用可逆流，给出的是**精确**似然，因而是朴素的极大似然目标，且模型各因子已知时重要性权可以精确计算。它与编号 89 同为该组的采样/生成模型论文，都用基于 KL 的训练与带有效样本量诊断的重要性重加权——编号 89 报告 Boltzmann 采样的相对有效样本量，编号 105 报告粒子滤波的 RESS；但编号 89 面对的是**已知的未归一化能量**且无数据，编号 105 面对的是**未知的后验**且只能通过模拟接触它。

它与确定性最小二乘 PDE 求解器（编号 90、102、94、103）在机器上无关：那些方法把一个正问题解到高精度而没有后验的概念，编号 105 则根本不解 PDE——对流扩散与 Burgers 系统在这里只作为数据生成器。

## 89：从已知未归一化能量采样

### 直觉

目标是从 $\pi(x)=Z^{-1}\exp(-U(x))$ 采样：能量 $U$ 可以逐点求值，归一化常数 $Z$ 不可得，而且**没有来自 $\pi$ 的训练数据**。最后这一条把普通生成模型的目标函数全部排除——Jensen-Shannon 散度、MMD、Wasserstein 距离都需要真实样本。

只剩下反向 KL 可用。但反向 KL 需要生成模型的密度 $\log p_D(x)$，这正是麻烦所在：要闭式的密度，解码器就得是**双射**（Boltzmann 生成器那一路），而双射限制了有效表达力。另一路是扩散型采样器，它不要求双射，但需要在训练回路里跑数值 SDE 或 ODE 求解器来算时间积分，代价高。

这篇论文同时去掉这两条限制。做法是：解码器用一个普通的（**非双射**的）高斯网络，然后在**隐空间**上再跑一段前向扩散。这段扩散有两个好处。第一，它的转移密度 $p_D(z_t\mid z_0)$ 解析可得，所以训练时可以直接精确采样 $z_t$，**根本不用解微分方程**。第二，它的时间反向过程给出一个可以当作「编码器」的对象，用它把反向 KL 放大成一个带隐路径的增广 KL——增广之后就不再需要 $\log p_D(x)$ 本身了。$\log Z$ 在这个界里只是一个与参数无关的常数，因此不知道它也能优化。

### 问题设定

目标为

$$
\pi(x)=\frac1Z\exp(-U(x)),\qquad Z=\int\exp(-U(x))\,\mathrm dx,\quad U:\mathbb R^d\to\mathbb R .
$$

解码器取 $x\mid z_0\sim p_D(x\mid z_0;\phi)$、$z_0\sim p_D(z_0)$ 为已知先验（标准多元正态），条件分布取高斯

$$
p_D(x\mid z_0;\phi)=\mathcal N\bigl(x\mid\mu(z_0;\phi),\,\Sigma(z_0;\phi)\bigr),
$$

$\mu$、$\Sigma$ 是网络，**不施加可逆性**。

### 推导

**第一步：增广的变分界。** 由于 $\pi$ 的样本不可得，只能用反向 KL，并把它放大成一个带隐变量的增广 KL：

$$
D_{\rm KL}\bigl(p_D(x)\,\|\,\pi(x)\bigr)
\le
D_{\rm KL}\bigl(p_{D}(z_{0})p_{D}(x\mid z_{0};\phi)\,\big\|\,\pi(x)p_{E}(z_{0}\mid x;\theta)\bigr)
$$

$$
=\mathbb E_{p_D(z_0)p_D(x\mid z_0;\phi)}
\left[\log\frac{p_D(z_0)p_D(x\mid z_0;\phi)}{p_E(z_0\mid x;\theta)}+U(x)\right]+\log Z,
$$

等号在 $p_E(z_0\mid x;\theta)$ 匹配解码器诱导的 $z_0\mid x$ 条件分布时取到。

**第二步：解码过程 = 解码器 + 隐空间前向扩散。**

$$
z_{0}\sim p_{D}(z_{0})\triangleq\mathcal N(\cdot\mid0,I),
\qquad x\mid z_{0}\sim p_{D}(x\mid z_{0};\phi),
$$

$$
\mathrm dz_{t}=f(z_{t},t)\,\mathrm dt+g(t)\,\mathrm dW_{t},\qquad t\in[0,T].
$$

两条结构性假设：（a）转移密度 $p_D(z_t\mid z_0)$ 解析可得，因此不需要解 Fokker-Planck 方程；（b）$z_T$ 近似无信息，$p_D(z_T)\approx p_D(z_T\mid z_0)$。**关键的结构观察是：边缘地看这段隐扩散只是噪声到噪声，但在 $x$ 的条件下，它把复杂的 $p_D(z_0\mid x)\propto p_D(z_0)p_D(x\mid z_0)$ 输运成可处理的 $p_D(z_T\mid x)=p_D(z_T)$。**

**第三步：反向 SDE 与编码过程。** 记 $\tilde t=T-t$，精确的反向时间 SDE 为

$$
\mathrm dz_{\tilde t}=-\Bigl(f(z_{\tilde t},\tilde t)-g(\tilde t)^{2}\nabla_{z_{\tilde t}}\log p_{D}(z_{\tilde t}\mid x)\Bigr)\mathrm d\tilde t
+g(\tilde t)\,\mathrm d\widetilde W_{\tilde t} .
$$

把其中不可得的分数换成网络 $s(z_{\tilde t},x,\tilde t;\theta)$，并把它与目标配对：$x\sim\pi(x)$、$z_T\sim p_E(z_T)\triangleq p_D(z_T)$。**编码过程从不被模拟**，它只提供损失里用到的密度 $p_E$；训练中 $x$ 与 $z_t$ 的样本全部来自解码过程。

**第四步：免模拟的无偏估计。** 由定理 1（见下）得到的可训练目标 $\mathcal L(\theta,\phi)$ 含一个散度项，用 Hutchinson 迹估计加对 $t$ 的 Monte Carlo 处理：

$$
\mathcal L(\theta,\phi)=\mathbb E_{p_D}\bigl[\log p_D(x\mid z_0;\phi)+U(x)\bigr]
+\mathbb E_{t\sim\mathcal U[0,T],\ \epsilon\sim p(\epsilon),\ (x,z_t)\sim p_D}\bigl[\mathcal L_t(x,z_t,\epsilon;\theta)\bigr],
$$

$$
\mathcal L_{t}(x,z_{t},\epsilon;\theta)=\frac{Tg(t)^{2}}{2}
\left(\|s(z_{t},x,t;\theta)\|^{2}
+2\,\frac{\partial\bigl[\epsilon^{\top}s(z_{t},x,t;\theta)\bigr]}{\partial z_{t}}\,\epsilon\right),
$$

$p(\epsilon)\sim\mathrm{Rademacher}^{d}$（$\mathbb E[\epsilon]=0$、$\mathrm{Cov}(\epsilon)=I$）。**「免模拟」就来自这里**：标准扩散 SDE 的条件分布 $p_D(z_t\mid z_0)$ 可以精确采样，所以训练时不需要任何数值 SDE/ODE 积分。

**第五步：把边界条件写进分数网络。** 精确分数在两端满足 $\nabla_{z_0}\log p_D(z_0\mid x)=\nabla_{z_0}[\log p_D(x\mid z_0)+\log p_D(z_0)]$ 与 $\nabla_{z_T}\log p_D(z_T\mid x)=\nabla_{z_T}\log p_D(z_T)$。于是网络写成一个恰好满足两端的插值：

$$
s(z,x,t;\theta)=\Bigl(1-\tfrac tT\Bigr)\nabla_{z_{0}}\bigl[\log p_{D}(x\mid z_{0}=z)+\log p_{D}(z_{0}=z)\bigr]
+\tfrac tT\,\nabla_{z_{T}}\log p_{D}(z_{T}=z)
+\tfrac tT\Bigl(1-\tfrac tT\Bigr)s'(z,x,t;\theta),
$$

$s'$ 是可训练网络。这使 $s$ 的误差在 $t=0$ 与 $t=T$ 处为零。

**第六步：广义 Hamilton 动力学（GHD）解码器。** 从 $z_0$ 生成一对初始位置与速度 $(y,v)$，再作蛙跳型更新（$l$ 为迭代指标、$\odot$ 为逐元素积）：

$$
v:=v-\frac{\epsilon_{0}e^{\epsilon_{0}\epsilon(l;\phi)}}{2}
\Bigl(\nabla U(y)\odot e^{\frac{\epsilon_{0}}{2}Q_{v}(y,\nabla U(y),l;\phi)}+T_{v}(y,\nabla U(y),l;\phi)\Bigr),
$$

$$
y:=y+\epsilon_{0}e^{\epsilon_{0}\epsilon(l;\phi)}
\Bigl(v_{k}\odot e^{\epsilon_{0}Q_{y}(v_{k},l;\phi)}+T_{y}(v_{k},l;\phi)\Bigr),
$$

再重复第一步。最终输出与其密度为

$$
x=y-\epsilon_{0}e^{\epsilon_{0}\eta(y;\phi)}\nabla U(y)+\sqrt{2\epsilon_{0}e^{\epsilon_{0}\eta(y;\phi)}}\,\xi,
\qquad \xi\sim\mathcal N(0,I),
$$

$$
p_{D}(x\mid z_{0};\phi)=\mathcal N\Bigl(x\Bigm|
y-\epsilon_{0}e^{\epsilon_{0}\eta(y;\phi)}\nabla U(y),\ 2\epsilon_{0}e^{\epsilon_{0}\eta(y;\phi)}I\Bigr).
$$

论文把最后一步解释为 Brownian 动力学 $\mathrm dy=-\nabla U(y)\mathrm dt+\mathrm dW_t$ 的有限步近似。$Q_v,T_v,Q_y,T_y$ 是可训练网络，步长参数化为可训练的正函数。声明的两点好处是：解码器直接用上 $\nabla U$，对多峰目标有帮助；以及只需要少数几次迭代。

**第七步：事后重要性重加权。**

$$
w(x,z_{0})=\frac{\exp(-U(x))\,p_{E}(z_{0}\mid x)}{p_{D}(z_{0})\,p_{D}(x\mid z_{0})}
\ \propto\ \frac{\pi(x)p_{E}(z_{0}\mid x)}{p_{D}(z_{0})p_{D}(x\mid z_{0})},
\qquad
\mathbb E_{\pi}[O(x)]\approx\frac{\sum_{n}w(x^{n},z_{0}^{n})O(x^{n})}{\sum_{n}w(x^{n},z_{0}^{n})} .
$$

障碍在于 $p_E(z_0\mid x)$ 是一个路径边缘。它经**概率流 ODE** 恢复：

$$
\mathrm dz_{t}=\Bigl(f(z_{t},t)-\tfrac12g(t)^{2}s(z_{t},x,t;\theta)\Bigr)\mathrm dt,
\qquad z_{T}\sim p_{E}(z_{T}),
$$

$$
\log p_{E}(z_{0}\mid x)=\log p_{E}(z_{T})
+\int_{0}^{T}\nabla\cdot\Bigl(f(z_{t},t)-\tfrac12g(t)^{2}s(z_{t},x,t;\theta)\Bigr)\mathrm dt .
$$

注意这里**确实**要解 ODE，但只在评价与重加权时解，训练时从不解。由此还得到一个自由能下界

$$
\log Z\ \ge\ \mathbb E_{p_{D}(x,z_{0})}\bigl[\log w(x,z_{0})\bigr],
$$

在 $p_D(x,z_0)=p_E(x,z_0)$ 时取紧。论文明确建议把这个下界当作**训练质量指标**——越大越好。相对有效样本量为

$$
\mathrm{rESS}=\frac{\bigl(\sum_{n=1}^{N}w(x^{n},z_{0}^{n})\bigr)^{2}}{N\sum_{n=1}^{N}w(x^{n},z_{0}^{n})^{2}}\in(0,1].
$$

### 定理

**定理 1（全篇唯一带编号的定理）。** 对上述解码与编码两个过程，

$$
D_{\rm KL}\bigl(p_D(x)\|\pi(x)\bigr)
\le D_{\rm KL}\bigl(p_D(x,z_{[\cdot]})\,\|\,p_E(x,z_{[\cdot]})\bigr)
=\mathcal L(\theta,\phi)
+\int_{0}^{T}\frac{g(t)^{2}}{2}\mathbb E_{p_D}\bigl[\|\nabla_{z_t}\log p_D(z_t)\|^{2}\bigr]\mathrm dt
+\log Z,
$$

其中 $z_{[\cdot]}=\{z_t\}_{t\in[0,T]}$ 是整条隐路径，而

$$
\mathcal L(\theta,\phi)=\mathbb E_{p_D}\bigl[\log p_D(x\mid z_0;\phi)+U(x)\bigr]
+\int_{0}^{T}\frac{g(t)^{2}}{2}\,\mathbb E_{p_D}
\Bigl[\|s(z_{t},x,t;\theta)\|^{2}+2\,\nabla_{z_{t}}\!\cdot s(z_{t},x,t;\theta)\Bigr]\mathrm dt .
$$

**取等条件是：解码过程中 $z_0\perp z_T$，且 $s(z,x,t;\theta)\equiv\nabla_{z_t}\log p_D(z_t\mid x;\phi)$。** 界中的后两项与参数无关，因此 $\mathcal L$ 单独就是可训练的目标。证明在附录 C，配套的命题在附录 A、B（$z_T$ 与 $x$ 的独立性）、D、E、F。

### 数值实验

**例 1：二维合成能量。** MoG2（两个各向同性高斯，$\sigma^2=0.5$，中心相距 10）、MoG2(i)（方差不等，$\sigma_1^2=1.5$、$\sigma_2^2=0.3$，中心相距 10）、MoG6（$\sigma^2=0.1$）、MoG9（$\sigma^2=0.3$）、Ring、Ring5。度量是对参考样本的 MMD，每种方法 50 万个样本。

| 目标    | EDG      | BG   | PIS  | V-HMC | L2HMC |
| ------- | -------- | ---- | ---- | ----- | ----- |
| MoG2    | 0.01     | 1.90 | —    | —     | —     |
| MoG2(i) | **0.50** | 1.63 | —    | 1.56  | 0.94  |
| MoG6    | 0.01     | 2.64 | —    | —     | —     |
| MoG9    | 0.02     | 0.07 | 0.42 | —     | —     |
| Ring    | 0.01     | 0.05 | —    | —     | —     |
| Ring5   | 0.02     | 0.18 | 0.78 | —     | —     |

**方差不等的双峰 MoG2(i) 是那个区分性的算例**：所有基线都 $\ge0.94$，只有 EDG 到 0.50。BG 在多峰情形上普遍很差，PIS 在 MoG9 与 Ring5 上崩掉，V-HMC 与 L2HMC 除 MoG2(i) 外都有竞争力。

**例 2：贝叶斯 logistic 回归。** Australian（15 个协变量）、German（25 个）、Heart（14 个）；准确率与 AUC 在 32 次独立运行上取平均，HMC 步长 0.01。EDG 在全部六个数字上最好：AU $84.96\pm1.67$ / $92.82\pm0.69$，GE $79.40\pm1.74$ / $82.79\pm1.46$，HE $88.02\pm3.90$ / $95.10\pm1.23$。

**例 3：Covertype。** 二分类，581012 个样本、54 个特征，层次贝叶斯后验；BG 与 EDG 用小批量无偏对数密度

$$
\log\pi(x)\approx\log p(x)+\frac{|\mathcal D_{\rm train}|}{|\mathcal B|}\sum_{(L,D)\in\mathcal B}\log p(L\mid D,x).
$$

EDG 准确率 $70.13\pm2.13$，而四个基线全部落在 $49.88$–$51.51$——**全篇差距最大的一处**。

**例 4：Lennard-Jones 系统。** LJ13（$d=39$）与 LJ55（$d=165$），非周期边界，参考 MCMC 数据取自文献。定性上 EDG 给出的原子间距直方图贴合测试数据，而 L2HMC 与 PIS 偏离到被放进单独的附录图里。定量上（256 个样本、10 个随机种子的相对有效样本量）：

| 体系 | EDG             | BG              | PIS             |
| ---- | --------------- | --------------- | --------------- |
| LJ13 | $0.132\pm0.048$ | $0.006\pm0.002$ | $0.005\pm0.001$ |
| LJ55 | $0.098\pm0.014$ | $0.004\pm0.000$ | $0.004\pm0.000$ |

约 20 倍的有效样本量改进，**但仍远低于 1**。论文明说不处理等变性，留待将来。

**例 5：二维 Ising 模型。** $16\times16$ 即 $d=256$，取连续松弛

$$
\pi(x)=\exp\Bigl(-\tfrac12x^{\top}(K(T)+\alpha I)^{-1}x\Bigr)\prod_{i=1}^{N}\cosh(x_{i}),
$$

$K$ 为依赖温度的对称矩阵、$\alpha$ 选得使 $K+\alpha I\succ0$；离散自旋由 $\pi(s\mid x)=\prod_i(1+e^{-2s_ix_i})^{-1}$ 恢复；并有

$$
\log Z=\log Z_{\text{Ising}}+\tfrac12\ln\det(K+\alpha I)-\tfrac N2\bigl[\ln(2/\pi)-\alpha\bigr].
$$

在 $T=2.0,\dots,2.7$ 上报告 $\log Z_{\rm Ising}$ 的下界估计（批量 $n=256$，均值的标准差为 $\mathrm{std}/\sqrt n$）：EDG 在 $T=2.0$ 处最大即最紧（$270.32\pm0.18$，NeuralRG $260.24\pm0.13$，PIS $210.17\pm0.43$），在 $2.1$、$2.2$、$2.3$、$2.5$、$2.7$ 处也最大，而 **$T=2.4$ 与 $T=2.6$ 处 NeuralRG 更好**。因此这条结论是「多数温度区间」而不是全部。

**例 6：消融（5000 个样本上的 MMD）。**

| 变体       | MoG2 | MoG2(i) | MoG6 | MoG9 | Ring | Ring5 |
| ---------- | ---- | ------- | ---- | ---- | ---- | ----- |
| VAE 无 GHD | 1.86 | 1.62    | 2.57 | 2.10 | 0.12 | 0.23  |
| VAE 有 GHD | 0.01 | 2.43    | —    | —    | 1.68 | —     |
| EDG 无 GHD | 0.06 | 1.01    | 0.04 | 0.05 | 0.02 | 0.04  |
| 完整 EDG   | 0.01 | 0.50    | 0.01 | 0.02 | 0.01 | 0.02  |

读法：**扩散编码器承担了大部分收益，GHD 再加上一层一致的改进——两个成分单独都不够。**

**代价。** 免模拟训练是一条结构性质而非实验结论：因为 $p_D(z_t\mid z_0)$ 解析，训练算法从不调用 SDE 或 ODE 求解器。可测的后果是训练显存 EDG **16.75 MB**，对 L2HMC 17.19、BG 22.41、PIS 55.83（在单张 NVIDIA RTX 2080 Ti / 11 GB 上、Mog2 任务、按 `torch.cuda.memory_allocated` 计）。论文承认 BG 的**单轮**时间更短，但主张 EDG 用更少轮次就达到更好的性能。另一条自陈的代价：重加权时 EDG 比 PIS **略慢**，因为神经 ODE 的对数密度需要算散度，而 PIS 不需要。

### 与其他论文的关系

它与[[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|密度流一支]]共享同一取向，但设定恰好相反：那里有样本、要估密度；这里有密度、要造样本。

它与编号 95、98 共享 Ling Guo + Hao Wu + Tao Zhou 这一作者群，三者都建立在隐变量与变分机器上：编号 89 用一个变分自编码器式的证据下界而把编码器换成扩散过程，编号 95 用带集合 transformer 编码器的条件变分自编码器与负证据下界，编号 98 用一个耦合到高斯过程的隐变量模型。贯穿三者的想法是**把手工设计的变分后验换成结构上更丰富的东西**——编号 89 是一条反向 SDE，编号 95 是一个置换不变的集合编码器，编号 98 是一个高斯过程插值的隐场。

在这个三人组之外，方法学上最近的是编号 105 与 107：都是把流或扩散模型用在动力学上的概率模型，而编号 107 同样把基分布选成能吸收目标中最难的那一部分。编号 89 用概率流 ODE 算对数密度所用的瞬时换元机器，正是连续正规化流用的那一套。它与本页的确定性一半（编号 90、102、94、101、103）形成鲜明对照：那些方法用最小二乘或配点求解一个**给定**的 PDE，而编号 89 求解的是一个**采样**问题，对目标的唯一抓手是 $U$ 与 $\nabla U$ 的逐点求值。两半共同的一点是：目标的支配对象（一个微分算子，或一个能量）都是解析已知的，因此都不需要带标签的数据。Ising 与 Lennard-Jones 这两组实验使编号 89 成为本清单中最接近计算统计物理的一篇。

> [!note] 覆盖进度
> 六篇的构造、损失、算法与主要结果均已按预印本或期刊全文逐式核对。定量结果的可核对程度不一：编号 94、103、105、89 的表格数值为逐项转录；编号 101 主要报告曲线而非数值表，本页因此只记录其设定、消融结构与定性结论；编号 81 的逐例误差位于第 5 节图表中，本次未逐项转录，只保留已核对的头条量级（相对标准多尺度网络提高两到三个数量级）。仍未核对的是编号 101 附录 A 中对高频放大的那段分析（正文的 Fourier 缩放论证已核对），以及编号 101、103、105 的期刊状态——三者在 arXiv 记录上均无期刊信息。背景文献一节的 Rahaman 定理 1、频率原则定理 1 与 2、MscaleDNN 构造均已按原文核对，其作用域限制也一并记下。

## 六篇的定位

| 编号 | 频率内容如何被处理               | 结构是否变动         |
| ---- | -------------------------------- | -------------------- |
| 81   | 从当前解的 DFT 后验捕获主模态    | 重建（两条准则）     |
| 94   | 解析地吸收进高斯波包拟设         | 不适用               |
| 101  | 固定基库 + 可学包络 + 交叉注意力 | 只增广词元，主干不动 |
| 103  | 不涉及（分裂结构 + $t$ 入输入）  | 不变动               |
| 105  | 不涉及（条件流 + 共享摘要）      | 不变动               |
| 89   | 不涉及（Boltzmann 采样）         | 不变动               |

理论分量的差别同样值得并排记下：

| 编号 | 定理内容                                                                                                                 |
| ---- | ------------------------------------------------------------------------------------------------------------------------ |
| 81   | 四条逼近定理：标准网络的界、降尺度网络的界（分子含 $(\max\{kh,C_1\})^{5q+3}$）、$C^1$ 复合的界、带限目标的**频率无关**界 |
| 94   | 无新证定理；$O(\sqrt\varepsilon)$ 建模误差与 $\mathcal E_t/\varepsilon$ 放大均引自前人                                   |
| 101  | 无定理；附录 A 的高频放大分析未核对，论文自陈缺少严格分析                                                                |
| 103  | 线性稳定性分析（含 $q_{\bm\theta}dt\ge2$ 时无条件不稳定）、定理 3.1 的 $dt^3$ Taylor 界、$C_1dt^3+C_2dt\,h^2$ 的误差分解 |
| 105  | 引理 B.1（两个最优摘要不必重合）与命题 B.1（充分摘要存在时二者相同）；对 FLUID 本身无收敛或误差定理                      |
| 89   | 定理 1：增广路径 KL 界与 $\mathcal L(\theta,\phi)$ 的显式分解，附取等条件                                                |

一条贯穿的判断：**频率内容要么被测量（编号 81 的 DFT 捕获、编号 101 的后验阈值），要么被解析地消掉（编号 94 的波包拟设），要么让模型自己在一个足够宽的固定基库上加权（编号 101 的交叉注意力）。凡是把频率内容当作先验固定下来的做法，都要承担频率不匹配的代价——编号 81 用 $\|F\|_{C^1[-1,1]}=\infty$ 把这个代价算成了一个具体的量。**

而编号 101 又给这条判断加了一个限定：**「偏向低频」并不是网络的固有属性，而是损失的属性。** 回归损失下网络偏低频，$k^4$ 加权的残差损失下它偏高频。所以真正要避免的不是某一个方向的偏差，而是在不测量的情况下假设方向。

另有一条与频率无关、但同样贯穿的判断：**这六篇里凡是有效的自适应，都建立在一个可以在当前近似上计算出来的指标上。** 编号 81 用解的离散 Fourier 谱决定重建哪些子网络，编号 101 用同一个谱决定注入哪些词元，编号 103 用物理残差 $\mathcal R_t$ 决定何时重训，编号 89 用自由能下界与相对有效样本量判断训练质量。这些指标都不需要参考解——这正是它们在实际问题上可用的原因。

## 本页原文

- J. Huang, R. You, and T. Zhou, [_Frequency-adaptive multi-scale deep neural networks_](https://doi.org/10.1016/j.cma.2025.117751), Comput. Methods Appl. Mech. Engrg. 437 (2025), 117751（预印本 [arXiv:2410.00053](https://arxiv.org/abs/2410.00053)）。
- Y. Wang, L. Guo, H. Wu, and T. Zhou, [_Energy-based diffusion generator for efficient sampling of Boltzmann distributions_](https://doi.org/10.1016/j.neunet.2025.108126), Neural Networks 194 (2026), 108126（预印本 [arXiv:2401.02080](https://arxiv.org/abs/2401.02080)）。
- J. Huang, R. You, and T. Zhou, [_Deep learning for the semi-classical limit of the Schrödinger equation_](https://doi.org/10.1016/j.jcp.2026.114869), J. Comput. Phys. 558 (2026), 114869（预印本 [arXiv:2509.04453](https://arxiv.org/abs/2509.04453)）。
- X. Feng, T. Tang, X. Wan, and T. Zhou, _Overcoming spectral bias via cross-attention_, [arXiv:2512.18586](https://arxiv.org/abs/2512.18586)，投稿 J. Comput. Phys.（arXiv 记录上无期刊信息）。
- J. Huang, Y. Qian, and T. Zhou, _PI-DOSnet: a physics-informed deep operator-splitting network for evolution partial differential equations_, [arXiv:2606.22514](https://arxiv.org/abs/2606.22514)，投稿 J. Comput. Phys.（arXiv 记录上无期刊信息）。
- T. Cui, X. Feng, C. Pei, X. Wan, and T. Zhou, _FLUID: flow-based unified inference for dynamics_, [arXiv:2604.07169](https://arxiv.org/abs/2604.07169)（arXiv 记录上无期刊引用字段，故「投稿 Comput. Methods Appl. Mech. Engrg.」未由预印本核实）。

背景文献：

- N. Rahaman, A. Baratin, D. Arpit, F. Draxler, M. Lin, F. Hamprecht, Y. Bengio, and A. Courville, _On the spectral bias of neural networks_, ICML 2019（[arXiv:1806.08734](https://arxiv.org/abs/1806.08734)）。
- Z.-Q. J. Xu, Y. Zhang, T. Luo, Y. Xiao, and Z. Ma, _Frequency principle: Fourier analysis sheds light on deep neural networks_, Commun. Comput. Phys.（[arXiv:1901.06523](https://arxiv.org/abs/1901.06523)）。
- Z. Liu, W. Cai, and Z.-Q. J. Xu, _Multi-scale deep neural network (MscaleDNN) for solving Poisson-Boltzmann equation in complex domains_, Commun. Comput. Phys. 28 (2020), 1970–2001（[arXiv:2007.11207](https://arxiv.org/abs/2007.11207)）。
