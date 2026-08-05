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

频谱偏差（频率原则）指网络快速收敛到低频成分而难以表示高频或强振荡特征。**但这个方向并非绝对**：编号 101 指出在 PDE 残差型损失下它会整体反转，详见下文。先按通常的方向看，多尺度网络用降尺度映射改善这一点：把紧支撑 Fourier 变换按同心环带 $\mathbb K_i=\{\bm k:(i-1)K_0\le|\bm k|\le iK_0\}$ 分解，每带按 $a_i>1$ 降尺度，逼近

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

## 101：固定基库，可学权重

### 一个方向相反的观察

这篇论文最值得记录的一点不是架构，而是它 3.1 节的一个观察：**在物理信息训练下，频谱偏差的方向会反过来。** 理由在 Fourier 端的缩放。由 $\widehat{\partial_xu}(k)=\mathrm ik\,\hat u(k)$，Deep Ritz 型能量项 $\|\partial_xu\|_{L^2}^2$ 给模态 $k$ 的权重是 $k^2$；由 $\widehat{u_{xx}}(k)=-k^2\hat u(k)$，平方 PDE 残差给模态 $k$ 的权重则是 $k^4$。也就是说微分算子本身**放大**高频，于是在残差损失下被欠解的反而是**低频**部分。

论文把这一点测了出来：对 $u=\sin(\pi x)+\sin(5\pi x)+\sin(20\pi x)$ 在 $k\in\{1,5,20\}$ 上看逐频相对误差 $\Delta_F(k)=|\hat u_k^{\rm pred}-\hat u_k|/|\hat u_k|$，纯回归损失下确实低频先收敛，而 Deep Ritz 与 PINN 损失下高频衰减得更快。

这直接决定了解的表示形式——两个网络相加

$$
u(\bm x;\theta)=u_h(\bm x;\theta_h)+\alpha\,u_\ell(\bm x;\theta_\ell),
$$

$u_h$ 是下述交叉注意力网络（管高频），$u_\ell$ 是普通全连接网络（管低频），边界条件按 $u_h=g$、$u_\ell=0$ 分开承担。

### 频率固定，衰减可学

基频 $\omega_m\sim\mathcal N(\bm 0,\sigma^{-2}I_{d_{\rm in}})$，按二进尺度 $\widetilde\omega_{m,k}=2^k\omega_m$（$k=0,\dots,K$）展成 $M=M_{\rm base}(K+1)$ 个频率，相位 $b_{m,k}\sim\mathrm{Uniform}(0,2\pi)$ 一次抽定。特征映射带一个**可学的幅度包络**：

$$
\phi(\bm x)=\sqrt{\tfrac1M}\Bigl[a_{m,k}\cos\bigl(\widetilde\omega_{m,k}^{\top}\bm x+b_{m,k}\bigr)\Bigr]_{(m,k)},
\qquad
a_{m,k}=\exp\bigl(-\beta\|\widetilde\omega_{m,k}\|_2\bigr).
$$

**频率本身固定且不可学，唯一可学的是衰减率 $\beta$**（经 softplus 参数化保正）。这与编号 81 的取舍正好互补：81 让特征去对齐真实频率，101 用一个足够宽的固定基库配一个可学的衰减包络。

### 交叉注意力：查询来自解，键值来自基库

记隐状态 $Q^{(l)}(\bm x)$、基库特征 $H(\bm x)$，并令 $Q_l=Q^{(l)}W_Q^{(l)}$、$K_l=H(\bm x)W_K^{(l)}$、$V_l=H(\bm x)W_V^{(l)}$，则

$$
\mathrm{CA}\bigl(Q^{(l)},H\bigr)=\mathrm{softmax}\Bigl(\frac{Q_lK_l^{\top}}{\sqrt{d_q}}\Bigr)V_l,
$$

$$
\widetilde Q^{(l)}=Q^{(l)}+\mathrm{CA}\bigl(Q^{(l)},H\bigr),
\qquad
Q^{(l+1)}=\widetilde Q^{(l)}+\sigma\bigl(W^{(l)}\widetilde Q^{(l)}+b^{(l)}\bigr).
$$

查询来自**隐状态**，键与值来自**频率基库**——所以是交叉注意力而非自注意力。要点在于 softmax 权重依输入而变，因此定义域的不同区域可以强调不同频带。

### 后验频率增强

先训一个初步模型，在均匀网格上做 DFT 得 $\hat u_{\theta,k}$，取 $\zeta=\max_{k\in B}|\hat u_{\theta,k}|$ 并按相对阈值挑出

$$
\mathcal K_{\rm post}=\bigl\{k\in B:\ |\hat u_{\theta,k}|>\lambda\zeta\bigr\},\qquad 0<\lambda<1,
$$

据此造确定性频率 $\omega_k^{\rm post}=2k\pi$ 与相应特征，在**词元维度**上拼接成 $H_{\rm aug}=[H_{\rm base};H_{\rm post}]$。新词元不是硬接入，而是经一个加性对数掩码平滑放开：

$$
A^{(l)}=\frac{Q_lK_l^{\top}}{\sqrt{d_q}}+\mathcal M^{(l)},
\qquad
\mathcal M^{(l)}=[\,\bm 0;\ \eta_l\bm 1\,],\quad \eta_l\le0,
$$

零块作用在 $H_{\rm base}$ 上、常数块 $\eta_l$ 作用在 $H_{\rm post}$ 上，再让 $\eta_l\uparrow0$ 逐步释放注入的词元。因此它**增广**基库而不重建基库，主干保持不动——这正是与编号 81 的分野所在。

## 103：把算子分裂做成网络，并把时间放回输入

设自治演化方程

$$
u_t=\mathcal Lu+\mathcal Nu=:\mathcal Fu\ \text{ 于 }\Omega\times(0,t^\star],
\qquad u(0,\cdot)=u_0,\qquad \mathcal Bu=0\ \text{ 于 }\partial\Omega,
$$

$\mathcal L$ 线性、$\mathcal N$ 非线性、$\mathcal B$ 为边界算子，均不含 $t$。Lie-Trotter 与 Strang 分裂分别给出

$$
u(T,\bm x)\approx e^{\tau_K\mathcal N}e^{\tau_K\mathcal L}\cdots e^{\tau_1\mathcal N}e^{\tau_1\mathcal L}u(0,\bm x),
$$

$$
u(T,\bm x)\approx e^{\frac{\tau_K}2\mathcal L}e^{\tau_K\mathcal N}e^{\frac{\tau_K}2\mathcal L}\cdots e^{\frac{\tau_1}2\mathcal L}e^{\tau_1\mathcal N}e^{\frac{\tau_1}2\mathcal L}u(0,\bm x).
$$

本文所基于的 DOSnet 把这个乘积直接写成网络：$\psi_{\bm\theta_T}=\psi_{\bm\theta_K}\circ\cdots\circ\psi_{\bm\theta_1}$，每个分裂块交替放可学线性层（卷积）与非线性层 $\phi_{\mathcal N_{l,i}}=e^{\tau_{l,i}\mathcal N}$，并约束 $\sum_{l,i}\tau_{l,i}=T$。**其中最值得注意的是：激活函数就是非线性子问题的精确流，而不是 ReLU 或 tanh。**

两处缺陷是本文的出发点。其一，数据驱动的算子学习需要大量成对样本，而演化问题的每一对都要跑一次求解器。其二，学到的算子并未被要求满足方程，因此**不能在任意时刻取值**——DOSnet 只输出终端时刻的解，中间块的输出只有靠与参考数据比对才对应得上中间时刻，且可取回的中间时刻数还受块数限制。

PI-DOSnet 的中心改动是取等步长 $\tau_{1,1}=\dots=\tau_{1,K}=dt=t/K$，并把线性部分的指数换成 $dt$ 的显式二阶 Taylor 展开，从而使整块成为 $t$ 的**显式函数**。于是 $t$ 回到输入里，解可在任意时刻求值，训练也不再依赖成对数据。

## 105：一个共享的摘要网络，把滤波与平滑绑在一起

### 设定

状态空间模型

$$
u_t=f(u_{t-1},\epsilon_{u,t}),\qquad y_t=h(u_t,\epsilon_{y,t}),
$$

要同时做**滤波** $p(u_t\mid y_{1:t})$ 与**平滑** $p(u_{1:t}\mid y_{1:t})$。论文对自己的约束比通常更强：只假设有一个**模拟器**，$f$、$h$ 与噪声分布的函数形式都不假设已知，转移密度与观测密度也不必可算。

核心困难是滤波的条件变量为整段历史 $y_{1:t}$，其维数随 $t$ 增长，固定输入的流用不上。FLUID 用多层 LSTM 把它压成定长摘要 $s_t=\mathrm{Enc}(y_{1:t};\psi)\in\mathbb R^h$，再用条件 KRnet 流写

$$
p(u_t\mid y_{1:t})\approx p_{\theta_1,\psi}(u_t\mid s_t).
$$

### 因果分解与共享

平滑一侧用的是一个精确分解

$$
p(u_{1:t}\mid y_{1:t})=p(u_t\mid y_{1:t})\prod_{k=1}^{t-1}p(u_k\mid u_{k+1},y_{1:k}),
$$

注意第二个因子的条件是 $y_{1:k}$ 而**不是** $y_{1:t}$——正是这一点让递推保持因果，于是同一个摘要 $s_k$ 可以直接复用。第二个流学 $p(u_t\mid u_{t+1},y_{1:t})\approx p_{\theta_2,\psi}(u_t\mid u_{t+1},s_t)$。

**两个流共用同一个摘要网络 $\psi$，这是论文冠名的贡献。** 联合目标为

$$
\min_{\theta_1,\theta_2,\psi}\
-\frac1{NT}\sum_{i=1}^N\sum_{t=1}^T\log p_{\theta_1,\psi}(u_t^i\mid s_t^i)
-\frac{\lambda}{N(T-1)}\sum_{i=1}^N\sum_{t=1}^{T-1}\log p_{\theta_2,\psi}(u_t^i\mid u_{t+1}^i,s_t^i).
$$

在 $\lambda=(T-1)/T$ 处有一个代数恒等式：两项重新组合成「逐时刻滤波似然」加「整条轨迹的平滑似然」，即这**一个**损失同时是两个任务的极大似然。共享还有一条充分性理由：若存在充分摘要 $S^{\dagger}$（即 $p(X\mid Y)=p(X\mid S^{\dagger})$ 几乎必然成立），则它同时是两个目标各自的最优摘要。

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
> 论文正文称平滑一致优于滤波，但 Table 9 在 $K=50$ 处给出平滑 RMSE 0.5423 而滤波 0.2605，Burgers 算例的 Table 6 在 $r^2=0.25$ 处有同样的反转。本页按表格数值记录。

## 89：从已知未归一化能量采样

目标是从 $\pi(x)=Z^{-1}\exp(-U(x))$ 采样，其中 $U$ 可逐点求值、$Z$ 不可得，且**没有来自 $\pi$ 的训练数据**。后一条排除了普通生成模型的目标函数——Jensen-Shannon 散度、MMD、Wasserstein 距离都需要真实样本。

论文针对的是两族方法各自的具体缺陷。基于正规化流的 Boltzmann 生成器必须用**双射**解码器才能得到闭式 $\log p_D(x)$，而双射限制了有效表达力；扩散型采样器（PIS 及其变体）需要在训练回路里跑数值 SDE/ODE 求解器来算时间积分，代价高。能量型扩散生成器同时去掉这两条限制：**非双射**解码器加**免模拟**的损失。

起点是一个变分界。由于 $\pi$ 的样本不可得，只能用反向 KL，并把它放大成一个带隐变量的增广 KL：

$$
D_{\rm KL}\bigl(p_D(x)\,\|\,\pi(x)\bigr)
\le
\mathbb E_{p_D(z_0)p_D(x\mid z_0;\phi)}
\left[\log\frac{p_D(z_0)p_D(x\mid z_0;\phi)}{p_E(z_0\mid x;\theta)}+U(x)\right]+\log Z,
$$

等号在 $p_E(z_0\mid x;\theta)$ 匹配解码器诱导的 $z_0\mid x$ 条件分布时取到。解码器取高斯

$$
p_D(x\mid z_0;\phi)=\mathcal N\bigl(x\mid\mu(z_0;\phi),\,\Sigma(z_0;\phi)\bigr),
$$

$\mu$、$\Sigma$ 是网络，**不施加可逆性**；解码之后在隐空间再跑一段前向扩散 $\mathrm dz_t=f(z_t,t)\mathrm dt+g(t)\mathrm dW_t$。要点是 $\log Z$ 只是与参数无关的常数，因此无须知道它也能优化上界。

这与[[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|密度流一支]]共享同一取向，但设定恰好相反：那里有样本、要估密度；这里有密度、要造样本。

> [!note] 覆盖进度
> 六篇的构造、损失与主要结果均已按预印本或期刊全文逐式核对。仍未核对的是编号 101 附录 A 中对高频放大的那段分析（正文的 Fourier 缩放论证已核对），以及编号 101、103、105 的期刊状态——三者在 arXiv 记录上均无期刊信息。

## 六篇的定位

| 编号 | 频率内容如何被处理               | 结构是否变动         |
| ---- | -------------------------------- | -------------------- |
| 81   | 从当前解的 DFT 后验捕获主模态    | 重建（两条准则）     |
| 94   | 解析地吸收进高斯波包拟设         | 不适用               |
| 101  | 固定基库 + 可学包络 + 交叉注意力 | 只增广词元，主干不动 |
| 103  | 不涉及（分裂结构 + $t$ 入输入）  | 不变动               |
| 105  | 不涉及（条件流 + 共享摘要）      | 不变动               |
| 89   | 不涉及（Boltzmann 采样）         | 不变动               |

一条贯穿的判断：**频率内容要么被测量（编号 81 的 DFT 捕获、编号 101 的后验阈值），要么被解析地消掉（编号 94 的波包拟设），要么让模型自己在一个足够宽的固定基库上加权（编号 101 的交叉注意力）。凡是把频率内容当作先验固定下来的做法，都要承担频率不匹配的代价——编号 81 用 $\|F\|_{C^1[-1,1]}=\infty$ 把这个代价算成了一个具体的量。**

而编号 101 又给这条判断加了一个限定：**「偏向低频」并不是网络的固有属性，而是损失的属性。** 回归损失下网络偏低频，$k^4$ 加权的残差损失下它偏高频。所以真正要避免的不是某一个方向的偏差，而是在不测量的情况下假设方向。

## 本页原文

- J. Huang, R. You, and T. Zhou, [_Frequency-adaptive multi-scale deep neural networks_](https://doi.org/10.1016/j.cma.2025.117751), Comput. Methods Appl. Mech. Engrg. 437 (2025), 117751（预印本 [arXiv:2410.00053](https://arxiv.org/abs/2410.00053)）。
- Y. Wang, L. Guo, H. Wu, and T. Zhou, [_Energy-based diffusion generator for efficient sampling of Boltzmann distributions_](https://doi.org/10.1016/j.neunet.2025.108126), Neural Networks 194 (2026), 108126（预印本 [arXiv:2401.02080](https://arxiv.org/abs/2401.02080)）。
- J. Huang, R. You, and T. Zhou, [_Deep learning for the semi-classical limit of the Schrödinger equation_](https://doi.org/10.1016/j.jcp.2026.114869), J. Comput. Phys. 558 (2026), 114869（预印本 [arXiv:2509.04453](https://arxiv.org/abs/2509.04453)）。
- X. Feng, T. Tang, X. Wan, and T. Zhou, _Overcoming spectral bias via cross-attention_, [arXiv:2512.18586](https://arxiv.org/abs/2512.18586)，投稿 J. Comput. Phys.
- J. Huang, Y. Qian, and T. Zhou, _PI-DOSnet: a physics-informed deep operator-splitting network for evolution partial differential equations_, [arXiv:2606.22514](https://arxiv.org/abs/2606.22514)，投稿 J. Comput. Phys.
- T. Cui, X. Feng, C. Pei, X. Wan, and T. Zhou, _FLUID: flow-based unified inference for dynamics_, [arXiv:2604.07169](https://arxiv.org/abs/2604.07169)，投稿 Comput. Methods Appl. Mech. Engrg.
