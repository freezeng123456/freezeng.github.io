---
title: 频谱偏差与生成式求解
description: 编号 81、89、94、101、103、105：网络偏向低频，而目标常常不是低频
lang: zh
translation: en/computational-mathematics/paper-notes/scientific-machine-learning/spectral-bias-and-generative-solvers
tags:
  - 论文笔记
  - 科学机器学习
  - 频谱偏差
---

> [!note] 本页覆盖
> 编号 **81**（_Comput. Methods Appl. Mech. Engrg._ 437, 2025）、**89**（_Neural Networks_ 194, 2026）、**94**（_J. Comput. Phys._ 558, 2026）、**101**（投稿 _J. Comput. Phys._，[arXiv:2512.18586](https://arxiv.org/abs/2512.18586)）、**103**（投稿 _J. Comput. Phys._）、**105**（投稿 _Comput. Methods Appl. Mech. Engrg._）。

## 频谱偏差是一个可以被测量的量

频谱偏差（频率原则）指网络快速收敛到低频成分而难以表示高频或强振荡特征。多尺度网络用降尺度映射改善这一点：把紧支撑 Fourier 变换按同心环带 $\mathbb K_i=\{\bm k:(i-1)K_0\le|\bm k|\le iK_0\}$ 分解，每带按 $a_i>1$ 降尺度，逼近

$$
f_{\theta}(\bm x)=\sum_{i=1}^{W} a_i^{d}\, f_{\theta_i}(a_i\bm x),
$$

而 $a_i$ 通常先验固定为 $2^{i-1}$。

## 81：先给出误差理论，再据它设计

### 两个具体的失效例子

论文的引导例子是 $f(x)=\sin(\pi x)+\sin(100\pi x)$，其理想尺度是 $a_1=\pi$、$a_2=100\pi$——实践中不可得。第二个问题是：带（随机）Fourier 特征的网络在嵌入频率**恰好正确**时比多尺度网络更准，而在**稍有偏差**时比普通网络更差。论文用 $k=38\pi$ 对目标频率 $40\pi$ 的例子把这一点具体化。

### 两条误差界

第一条解释降尺度为何有效。对 $(p,C_0)$-光滑的 $f:\mathbb R^d\to\mathbb R$（$p=q+s$，$q\in\mathbb N_0$，$s\in(0,1]$）与充分大的 $M$，存在带**降尺度**映射 $\Phi(\bm x)=k\bm x$ 的网络使

$$
\|f_{\text{net}}-f\|_{\infty,[-h,h]^d}
\ \le\ \frac{C_2\bigl(\max\{kh,\,C_1\}\bigr)^{5q+3}}{M^{2p}} .
$$

与不带缩放映射的标准网络界相比，要点是**缩小 $kh$ 就缩小了分子**。

第二条解释为何应当用真实主模态作特征。若 $f\in L^2([-\pi,\pi]^d)$ 带限，$f(\bm x)=\sum_{\bm k\in\mathbb B}[b_{\bm k}\cos(\bm k\cdot\bm x)+c_{\bm k}\sin(\bm k\cdot\bm x)]$、$|\mathbb B|=N$ 且系数有界，则用 $\Phi(\bm x)=(\dots,\sin(\bm k\cdot\bm x),\cos(\bm k\cdot\bm x),\dots)_{\bm k\in\mathbb B}$ 作 Fourier 特征时存在网络使

$$
\|f_{\text{net}}-f\|_{\infty,[-\pi,\pi]^d}\ \le\ \frac{C_2\sqrt N}{M^{2p}} .
$$

**这个界不显含频率。** 这正是「用实际主模态构造嵌入」的理论依据。

### 混合特征嵌入

论文的新输入映射把多尺度网络的线性降尺度与正弦 Fourier 特征合起来：

$$
\Phi[\bm k](\bm x)=\begin{bmatrix}\bm k\cdot\bm x\\ \cos(\bm k\cdot\bm x)\\ \sin(\bm k\cdot\bm x)\end{bmatrix}.
$$

动机是一个具体计算：把 $f(x)=\sin(40\pi x)$ 经 $F(y)=\sin(\hat k\arcsin y)$（$\hat k=40\pi/k$）逼近时，$F'(y)=\hat k\cos(\hat k\arcsin y)/\sqrt{1-y^2}$，故 $\lim_{y\to1}|F'(y)|=\infty$，因此只要 $|\hat k|\ne1$ 就有 $\|F\|_{C^1[-1,1]}=\infty$，相应定理的假设失效。**这就是纯 Fourier 特征在频率不匹配时退化的精确原因**，而线性项的加入正是为了补上这一处。

### 后验捕获主模态

对当前网络解做离散 Fourier 变换 $\hat f_{\text{net},\bm k}=\int_{[0,1]^d}f_{\text{net}}(\bm x)e^{-\mathrm i2\pi\bm k\cdot\bm x}\mathrm d\bm x$，保留模最大的 $N_0$ 个系数，其中 $N_0$ 由能量占比条件确定：

$$
\sum_{j=1}^{N_0}\bigl|\hat f_{\text{net},\bm k_j}\bigr|^2
\ \ge\ (1-\delta)\,\|f_{\text{net}}\|^2_{L^2([0,1]^d)},
\qquad 0\le\delta<1 .
$$

据此按两条准则重建网络。若捕获模态数 $N_0\le M_0$（子网络数），建 $N_0$ 个子网络，第 $j$ 个的输入为 $\Phi[\bm k_j](\bm x)$，输出

$$
y=\sum_{j=1}^{N_0} h_j\,y_j,\qquad h_j=\hat u_{\text{net},0,\bm k_j},
$$

系数被复用以加速收敛。若 $N_0>M_0$，则**不扩网络**，而把 $\mathbb B_1=\{\bm k_1,\dots,\bm k_{N_0}\}$ 分成 $M_0$ 块，第 $j$ 个子网络接收该块内各 $\bm k$ 的 $\Phi[\bm k](\bm x)$，输出

$$
y=\sum_{j=1}^{M_0}h_j\,y_j,\qquad h_j=\sum_{\bm k\in\mathbb B^j_1}\bigl|\hat u_{\text{net},0,\bm k}\bigr| .
$$

一处实现细节值得记录：最后一层写成 $u_{\text{net}}=WG+b$，其中 $W=(h_1,h_2,\dots)$ 与 $b=0$ 是**固定的、不可学习的**，论文报告这样既更准也更省。

算法有两条停机准则：固定的最大轮数，以及捕获特征集不再变化。论文报告的量级是频率自适应版本相对标准多尺度网络把精度提高两到三个数量级。实验覆盖 Poisson 方程、热方程、波方程与半经典极限附近的 Schrödinger 方程。

## 94：半经典极限——把空间维度整体消掉

### 问题的量级

目标方程为

$$
\psi_t=\frac{\mathrm i\varepsilon}{2}\Delta\psi-\frac{\mathrm i}{\varepsilon}V(\bm x)\psi,
\qquad
\psi(\bm x,0)=\varphi(\bm x)\exp\bigl(\mathrm i\phi(\bm x)/\varepsilon\bigr),
$$

$\varepsilon$ 是无量纲 Planck 常数。$\varepsilon\ll1$ 时解在空间与时间上都以 $O(\varepsilon^{-1})$ 的尺度高频振荡。经典分裂求解器（Strang 分裂在空间上谱精度、时间上 $O(\Delta t^2/\varepsilon)$；Chin-Chen 空间谱精度、时间四阶）的全局截断误差都带一个 $1/\varepsilon$ 因子，因此网格与时间步必须随 $\varepsilon$ 一同缩小。

普通物理信息网络更差。论文直接给出证据：把朴素网络用于 $(-2\pi,2\pi)\times[0,1]$ 上的该方程，在 $\varepsilon=0.1$ 时相对 $L^2$ 误差**已经超过 1**。

### Heller 高斯波包与关键约化

论文采用 Heller 的高斯波包拟设，一维情形为

$$
\psi(x,t)=\exp\left[\frac{\mathrm i}{\varepsilon}
\Bigl(\alpha(t)\bigl(x-q(t)\bigr)^2+p(t)\bigl(x-q(t)\bigr)+\gamma(t)\Bigr)\right].
$$

包络中心在 $q(t)$、标准差正比于 $\sqrt{\varepsilon/\alpha_{\rm im}(t)}$，$x=q(t)$ 处的振荡波长为 $2\pi\varepsilon/p(t)$，而 $\mathrm{Re}\,\alpha$ 在尾部产生更细的振荡。

**关键在于：把这个拟设代入方程后得到关于 $(\alpha,q,p,\gamma)$ 的常微分方程系统，空间变量 $x$ 被完全消掉。** 于是网络要学的不再是一个在空间上高频振荡的函数，而是若干个随时间演化的标量——高频结构被解析地放进了拟设。这是[[computational-mathematics/paper-notes/scientific-machine-learning/index|本专题]]中「把结构写进架构」这一取向在振荡问题上的实现。

论文处理的第二个独立问题是：更换初值会强制重训，因此需要一个算子学习的表述。

## 101：保持基库固定，学注意力权重

编号 101 与编号 81 处理同一个频谱偏差问题，但角度不同。论文对已有修补（随机 Fourier 特征、多尺度网络等）的批评在于它们把频率内容的选择固定在架构里。它的替代方案是保持一组**固定的**多尺度随机 Fourier 基，并学习依输入变化的**交叉注意力**权重。

两条路线的对比值得写清楚：编号 81 按捕获的主模态**重建网络**（离散的、需要重训的），编号 101 让权重**连续可学**（不重建结构）。前者的优点是特征恰好对齐目标频率，后者的优点是不需要离散的重建步骤。

## 103 与 105：把结构写进网络的另两种方式

- **103（PI-DOSnet）** 把**算子分裂**的结构写进网络，用于演化型方程。Lie-Trotter 与 Strang 分裂把一个复杂演化拆成若干可分别处理的子步；把这个结构显式写进网络意味着网络的层次对应分裂的子步，而不是一个无结构的映射。
- **105（FLUID）** 把编号 62 的随机场流路线推进为面向动力系统的统一推断框架，用条件流替代带 Karhunen-Loève 结构的参考场。

## 89：从已知未归一化能量采样

编号 89 处理的是**从 Boltzmann 分布采样**：给定未归一化能量 $E(x)$，目标分布是 $\pi(x)\propto e^{-E(x)}$，归一化常数不可得。这属于生成模型而不是 PDE 求解，但它与[[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|密度流一支]]共享同一取向：用一个可逆或可采样的模型表示分布，从而把「采样」变成「训练一个模型」。

论文的方法是能量型扩散生成器。扩散或基于分数的方法与正规化流的差别在于：流给出显式密度但受可逆性限制，扩散模型表达力更强但密度只能通过一个随机微分方程隐式得到。在只有能量已知（没有样本）的设定下，训练信号只能来自能量本身，因此需要一个反向 KL 或类似的变分目标。

> [!note] 覆盖进度
> 编号 81 与 94 的构造、误差界与算法已核对。编号 89、101、103、105 的具体损失、架构与定理本页未逐式核对；上述内容限于可从题名、摘要与相邻论文交叉引用确认的部分。

## 六篇的定位

| 编号 | 频率内容如何被处理        | 是否重建结构   |
| ---- | ------------------------- | -------------- |
| 81   | 从当前解的 DFT 后验捕获   | 是（两条准则） |
| 94   | 解析地放进高斯波包拟设    | 不适用         |
| 101  | 固定基库 + 可学交叉注意力 | 否             |
| 103  | 不涉及（算子分裂结构）    | 否             |
| 105  | 不涉及（条件流表示）      | 否             |
| 89   | 不涉及（Boltzmann 采样）  | 否             |

一条贯穿的判断：**频率内容要么被测量（编号 81），要么被解析地消掉（编号 94），要么让模型自己在一个足够宽的基库上加权（编号 101）。假设它是低频的做法在这三条之外，而那正是频谱偏差的定义。**

## 本页原文

- J. Huang, R. You, and T. Zhou, [_Frequency-adaptive multi-scale deep neural networks_](https://doi.org/10.1016/j.cma.2025.117751), Comput. Methods Appl. Mech. Engrg. 437 (2025), 117751（预印本 [arXiv:2410.00053](https://arxiv.org/abs/2410.00053)）。
- Y. Wang, L. Guo, H. Wu, and T. Zhou, [_Energy-based diffusion generator for efficient sampling of Boltzmann distributions_](https://doi.org/10.1016/j.neunet.2025.108126), Neural Networks 194 (2026), 108126。
- J. Huang, R. You, and T. Zhou, [_Deep learning for the semi-classical limit of the Schrödinger equation_](https://doi.org/10.1016/j.jcp.2026.114869), J. Comput. Phys. 558 (2026), 114869（预印本 [arXiv:2509.04453](https://arxiv.org/abs/2509.04453)）。
- X. Feng, T. Tang, X. Wan, and T. Zhou, _Overcoming spectral bias via cross-attention_, [arXiv:2512.18586](https://arxiv.org/abs/2512.18586)，投稿 J. Comput. Phys.
- J. Huang, Y. Qian, and T. Zhou, _PI-DOSnet: a physics-informed deep operator-splitting network for evolution partial differential equations_，投稿 J. Comput. Phys.
- T. Cui, X. Feng, C. Pei, X. Wan, and T. Zhou, _FLUID: flow-based unified inference for dynamics_，投稿 Comput. Methods Appl. Mech. Engrg.
