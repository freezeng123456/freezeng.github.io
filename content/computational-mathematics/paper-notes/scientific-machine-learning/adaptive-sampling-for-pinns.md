---
title: 采样点放在哪里：失效概率与重要采样
description: 编号 66、70、73、76、80：把配点选择变成可靠性分析与方差缩减问题
lang: zh
translation: en/computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns
tags:
  - 论文笔记
  - 科学机器学习
  - 自适应采样
---

> [!note] 本页覆盖
> 编号 **66**（_Comput. Methods Appl. Mech. Engrg._ 400, 2022）、**70**（_SIAM J. Sci. Comput._ 45(4), 2023）、**73**（_Commun. Appl. Math. Comput._ 6, 2024）、**76**（_CSIAM Trans. Appl. Math._ 5(3), 2024）、**80**（_Commun. Appl. Math. Comput._ 7(3), 2025）。

![用可靠性分析决定下一批配点](assets/diagrams/tao-zhou-papers/zh/failure-informed-sampling.svg)

物理信息神经网络的精度对配点位置高度敏感。这组工作的共同判断是：**配点选择不是实现细节，而是可以用一个明确的数学目标来驱动的对象**。四种目标各不相同——残差过大的概率、变分损失的方差、条件期望的无偏性——但都归结为「先定义一个可估计的标量，再据它生成新点」。

## 66：先把算子本身写成期望

### 直觉

分数阶 Laplacian 是非局部的：要算它在一个点上的值，原则上需要解在**整个空间**的信息。早先的 fPINN 用一张辅助网格离散这个积分，于是每个配点都要访问网格上的许多邻点，存储与计算随维数指数增长。

编号 66 的出发点是一个很朴素的观察：那个奇异积分

$$
\mathrm{P.V.}\!\int_{\mathbb R^d}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\,\mathrm dy
$$

**本来就是一个加权平均**——权重 $\|x-y\|^{-d-\alpha}$ 只是没有归一化。把它归一化成一个概率密度，积分就变成期望，而期望可以用少量随机样本无偏地估计。于是分数阶算子既不需要组装，也不需要存储：每个配点只要随机抽几个方向和几个半径，评估网络在几个偏移点上的值即可。维数只出现在「抽一个 $d$ 维方向」这一步，而这一步的代价是线性的。

这个思路有一个陷阱，而它恰恰是本文摘要点名的关键成分。物理信息损失要的是残差的**平方**，而带噪估计的平方是有偏的：$\mathbb E[\widehat L^2]=(\mathbb E\widehat L)^2+\mathrm{Var}(\widehat L)$。如果直接平方，优化器最小化的就不是真实残差，而是残差加上估计量的方差——后者与 PDE 无关，却可以通过让网络的二阶差分本身变小来降低。论文的解决办法是抽**两组**独立随机数，把两个独立估计相乘而不是平方一个估计。

### 问题设定

模型方程为分数阶对流扩散：

$$
L[u(x,t)]:=\frac{\partial^\gamma u(x,t)}{\partial t^\gamma}
+c\,(-\Delta)^{\alpha/2}u(x,t)+v\cdot\nabla u(x,t)=f(x,t),
$$

其中 Caputo 时间导数与积分型分数阶 Laplacian 分别为

$$
\frac{\partial^\gamma u(x,t)}{\partial t^\gamma}
\triangleq\frac{1}{\Gamma(1-\gamma)}\int_0^t (t-\tau)^{-\gamma}
\frac{\partial u(x,\tau)}{\partial\tau}\,\mathrm d\tau,
\qquad 0<\gamma<1,
$$

$$
(-\Delta)^{\alpha/2}u(x)\triangleq C_{d,\alpha}\,\mathrm{P.V.}\!
\int_{\mathbb R^d}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\,\mathrm dy .
$$

训练集写成三部分 $\mathcal D=(\mathcal D_f,\mathcal D_g,\mathcal D_u)$，分别是方程点、边界／初值点与观测点，总损失

$$
\mathcal{LOSS}(\theta)=w_{equ}L_{equ}(\theta)+w_g L_g(\theta)+w_u L_u(\theta).
$$

### 推导

**分数阶 Laplacian 的内外分裂。** 以 $x$ 为心、半径 $r_0$ 把 $\mathbb R^d$ 分成球内与球外两块，两块分别处理，因为奇性只在球内。

球内部分。取极坐标 $y=x+r\xi$，$\mathrm dy=r^{d-1}\mathrm dr\,\mathrm dS(\xi)$，被积函数的 $\|x-y\|^{-d-\alpha}$ 变成 $r^{-d-\alpha}$，两个 $r$ 的幂次相消后剩下 $r^{-1-\alpha}$。再对 $\xi\mapsto-\xi$ 作对称化（$S^{d-1}$ 在这个对合下不变），主值意义下的一阶差分被替换成二阶差分：

$$
u(x)-u(x+r\xi)\ \longrightarrow\ \tfrac12\bigl[2u(x)-u(x+r\xi)-u(x-r\xi)\bigr].
$$

这一步是全篇的技术核心：**对称化把 $O(r)$ 的一阶项消掉，剩下的二阶差分是 $O(r^2)$，恰好抵消 $r^{-1-\alpha}$ 里 $\alpha<2$ 时残余的奇性。** 提出 $r^2$ 后半径方向的权重是 $r^{1-\alpha}$，它在 $[0,r_0]$ 上可积，归一化常数为 $\int_0^{r_0}r^{1-\alpha}\mathrm dr=r_0^{2-\alpha}/(2-\alpha)$。于是

$$
\int_{y\in B_{r_0}(x)}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\mathrm dy
=\frac{|S^{d-1}|\,r_0^{2-\alpha}}{2(2-\alpha)}\;
\mathbb E_{\xi,\,r\sim f_I}\!\left[\frac{2u(x)-u(x-r\xi)-u(x+r\xi)}{r^2}\right],
$$

$$
f_I(r)=\frac{2-\alpha}{r_0^{2-\alpha}}\,r^{1-\alpha}\mathbf 1_{r\in[0,r_0]},
\qquad r/r_0\sim\mathrm{Beta}(2-\alpha,1).
$$

Beta 分布是从累积分布函数读出来的：$f_I$ 的原函数是 $(r/r_0)^{2-\alpha}$，而 $\mathrm{Beta}(a,1)$ 的累积分布函数正是 $x^a$，因此 $r/r_0$ 可以由一次幂变换直接采样，不需要拒绝法。

球外部分不含奇性，因此不需要除以 $r^2$，半径权重是 $r^{-1-\alpha}$，归一化常数 $\int_{r_0}^\infty r^{-1-\alpha}\mathrm dr=r_0^{-\alpha}/\alpha$：

$$
\int_{y\notin B_{r_0}(x)}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\mathrm dy
=\frac{|S^{d-1}|\,r_0^{-\alpha}}{2\alpha}\;
\mathbb E_{\xi,\,r\sim f_O}\bigl[2u(x)-u(x-r\xi)-u(x+r\xi)\bigr],
$$

$$
f_O(r)=\alpha\,r_0^{\alpha}\,r^{-1-\alpha}\mathbf 1_{r\in[r_0,\infty)},
\qquad r_0/r\sim\mathrm{Beta}(\alpha,1).
$$

**数值保护。** 二阶差分商在 $r\to0$ 时是良定的，

$$
\lim_{r\to0}\frac{2u(x)-u(x-r\xi)-u(x+r\xi)}{r^2}
=\partial^2_r u(x+r\xi)\big|_{r=0},
$$

但浮点运算下分子是两个几乎相等的量相减，$1/r^2$ 会把相消误差放大。因此实现中把半径换成 $r_\epsilon=\max\{\epsilon,r_I\}$。合起来给出网络上的估计式

$$
(-\Delta)^{\alpha/2}u_{NN}(x)=C_{d,\alpha}\frac{|S^{d-1}|r_0^{2-\alpha}}{2(2-\alpha)}
\mathbb E_{\xi,r_I\sim f_I}\!\left[
\frac{2u_{NN}(x)-u_{NN}(x-r_\epsilon\xi)-u_{NN}(x+r_\epsilon\xi)}{r_\epsilon^2}\right]
$$

$$
\qquad+\;C_{d,\alpha}\frac{|S^{d-1}|r_0^{-\alpha}}{2\alpha}\,
\mathbb E_{\xi,r_o\sim f_O}\bigl[2u_{NN}(x)-u_{NN}(x-r_o\xi)-u_{NN}(x+r_o\xi)\bigr].
$$

**Caputo 时间导数**用同一手法处理，结果是一个期望加一个显式项：

$$
\frac{\partial^\gamma u_{NN}(x,t)}{\partial t^\gamma}
=\frac{\gamma}{1-\gamma}\,t^{1-\gamma}\,
\mathbb E_{\tau\sim f_{I,t}}\!\left[\frac{u_{NN}(x,t)-u_{NN}(x,t-\tau)}{\tau}\right]
+\frac{u_{NN}(x,t)-u_{NN}(x,0)}{t^{\gamma}},
$$

$$
f_{I,t}(\tau)=(1-\gamma)\tau^{-\gamma}\mathbf 1_{\tau\in[0,1]},
\qquad \tau\sim\mathrm{Beta}(1-\gamma,1),
\qquad \tau_\epsilon=\max\{\tau,\epsilon_t t^{-1}\}.
$$

第二项显式含初值 $u_{NN}(x,0)$，这与 Caputo 导数「先减去初值再求分数阶导」的定义是一致的。三项（分数阶 Laplacian、Caputo 导数、由自动微分给出的对流项 $v\cdot\nabla u_{NN}$）合成一个随机算子 $\widehat L$。

**无偏性构造。** 把 $\widehat L$ 的平方换成两个独立副本的乘积：

$$
\hat L_{equ}(\theta)=\frac{1}{mN_u}\sum_{i,j}
\widehat L\bigl[u_{NN}(x_i,t_i;\theta);\epsilon,\epsilon_t,\tau_j,\xi_j,r_{Ij},r_{oj}\bigr]\cdot
\widehat L\bigl[u_{NN}(x_i,t_i;\theta);\epsilon,\epsilon_t,\tau'_j,\xi'_j,r'_{Ij},r'_{oj}\bigr].
$$

两组随机参数独立，因此乘积的期望等于期望的乘积，即真实残差的平方；方差项不再出现。

### 定理

本文不给收敛定理。其理论内容是上面两个期望恒等式的**精确性**，加上一条无偏性陈述：在 $\epsilon=\epsilon_t=0$（即不加数值保护）且忽略舍入误差的前提下，

$$
\mathbb E\bigl[\hat L_{equ}(\theta)\bigr]=L_{equ}(\theta).
$$

条件值得记住：实际实现里 $\epsilon,\epsilon_t>0$，因此严格意义上的无偏性只在极限下成立，$\epsilon$ 同时是稳定性与偏差之间的一个权衡旋钮。

### 数值实验

算法本身很短（论文的 Algorithm 1）：给定 $\mathcal D_f,\mathcal D_g,\mathcal D_u$，每步从训练数据抽 $N$ 个快照；对每个残差点从 $\Omega\subset\mathbb R^d$ 与时间区间均匀抽 $x_i,t_i$，再抽两组独立的 $\{\tau_i,\xi_i,r_{Ii},r_{oi}\}_{i=1}^m$ 与 $\{\tau'_i,\xi'_i,r'_{Ii},r'_{oi}\}_{i=1}^m$；按上式组装损失，用 Adam 更新，直至收敛。注意随机性有两层：配点的随机性与算子估计的随机性，后者每步重抽。

实验分三组：

| 实验组 | 目标                                             |
| ------ | ------------------------------------------------ |
| 一     | 高维积分型分数阶 Laplacian 方程（正问题）        |
| 二     | 时空分数阶 PDE 的参数辨识，即从观测反推分数阶阶数 |
| 三     | 带随机输入的分数阶扩散方程                       |

论文的定量主张是**总体计算代价低于 fPINN**，从而打开高维分数阶 PDE 的可能性。

> [!warning] 数值结果的可核实范围
> 本页依据的核对材料确认了实验的分组与设定，但未能从可访问的渲染中逐个转录误差表：**具体的维数上限与误差量级未经核实**，因此本页不给数字。需要引用具体误差时应回查期刊版第 4 节。

### 与其他论文的关系

这个内外分裂估计被编号 72 原样引入，作为其 Lemma 3.1（同样的 $f_I$、$f_O$、Beta 采样与 $r_\epsilon$ 保护），差别在于施加的对象：编号 66 作用在一个通用 PINN 代理 $u_{NN}$ 上，编号 72 作用在一个**正规化流密度**上，因而还要顾及非负性与归一化。编号 72 同时给出一条完全避开随机性的替代路线（高斯径向基辅助模型），见[[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|密度流一页]]。

> [!note] 题名与原文核对
> 预印本题为 "Monte Carlo PINNs: deep learning approach for forward and inverse problems involving high dimensional fractional partial differential equations"，期刊版为 "Monte Carlo fPINNs: Deep learning method for ..."。另外，Caputo 导数定义的前因子在预印本渲染中写作 $\Gamma(1-\alpha)$，而指数用的是 $\gamma$；按 $0<\gamma<1$ 的设定应为 $\Gamma(1-\gamma)$。上面的期望恒等式链在渲染版中是从不带 $1/\Gamma(1-\gamma)$ 前因子的积分出发的，引用前因子时建议回查期刊版。

## 70：把「残差过大」定义成失效事件

### 直觉

标准的残差自适应细化在均匀候选池中取残差最大的 $m$ 个点。它失效的场景很具体：当高残差区只占定义域很小一部分时——尖峰、内层、无界区域上支撑集中的解——**均匀候选几乎落不进去**，于是「取最大的 $m$ 个」取到的仍是一堆无关紧要的点。

编号 70 换了一个学科来借工具。结构可靠性分析里有一个标准问题：给定一个荷载分布与一个「结构失效」判据，估计失效概率，而失效概率往往极小，因此不能用朴素 Monte Carlo，要用重要采样把提议密度推向失效区。把「结构失效」换成「残差超过容差」，整套机器可以逐字搬过来：极限状态函数、失效集、自适应重要采样。**副产品是一个后验误差指标**——失效概率本身，其地位与自适应有限元中的后验误差估计完全对应，而且（见下面的定理）它真的会出现在误差界里。

### 问题设定

物理信息神经网络的损失是

$$
\mathcal L(\theta)=\mathcal L_c(\theta)+\lambda\mathcal L_b(\theta),
\qquad
\mathcal L_c(\theta)=\frac{1}{N_c}\sum_{i=1}^{N_c}\bigl|r(x^c_i;\theta)\bigr|^2,
\qquad
\mathcal L_b(\theta)=\frac{1}{N_b}\sum_{i=1}^{N_b}\bigl|b(x^b_i;\theta)\bigr|^2 .
$$

一般的极限状态函数写成 $g(x)=\mathcal Q(x)-\epsilon_r$，其中 $\mathcal Q$ 把定义域映到某个关心的量；对 PINN 取 $\mathcal Q(x)=|r(x;\theta)|$，于是

$$
g(x)=\bigl|r(x;\theta)\bigr|-\epsilon_r,
$$

零水平集把定义域分成安全集 $\Omega_{\mathcal S}=\{g<0\}$ 与失效集 $\Omega_{\mathcal F}=\{g>0\}$。在先验 $\omega(x)$ 下的失效概率

$$
P_{\mathcal F}=\int_\Omega \omega(x)\,\mathbb I_{\Omega_{\mathcal F}}(x)\,\mathrm dx
$$

充当后验误差指标：当 $P_{\mathcal F}<\epsilon_p$ 时宣布网络可靠。

与残差自适应细化的一个具体差别值得指出：基于 Monte Carlo 的加点会把**所有**落入 $\Omega_{\mathcal F}$ 的候选点都加进去，因此每轮加点个数 $m$ 是变化的；残差细化每轮固定加 $m$ 个。

### 推导

朴素估计 $\hat P^{MC}_{\mathcal F}=\frac{1}{|\mathcal S|}\sum_{x\in\mathcal S}\mathbb I_{\Omega_{\mathcal F}}(x)$（$\mathcal S\sim\omega$）在 $P_{\mathcal F}$ 很小时几乎全是零。论文用自适应重要采样逐步把提议密度推向失效区：从 $h_1=\omega$ 出发，第 $k$ 步抽 $N_1$ 个样本，按极限状态函数**降序**排列得到 $\widetilde x^k_1,\dots,\widetilde x^k_{N_1}$，令

$$
N_\eta=\max_{1\le i\le N_1}\{i:\ g(\widetilde x_i)>0\},
\qquad
N_p=\lfloor p_0N_1\rfloor .
$$

$N_\eta$ 是落在失效集里的样本数，$N_p$ 是一个固定分位数。若 $N_\eta<N_p$（失效样本还太少，说明提议还不够集中），用前 $N_p$ 个样本的矩更新截断高斯提议：

$$
\mu_{k+1}=\frac{1}{N_p}\sum_{i=1}^{N_p}\widetilde x^k_i,
\qquad
\Sigma_{k+1}=\frac{1}{N_p-1}\sum_{i=1}^{N_p}
\bigl(\widetilde x^k_i-\mu_{k+1}\bigr)\otimes\bigl(\widetilde x^k_i-\mu_{k+1}\bigr),
$$

并令 $h_{k+1}=\mathcal N_T(\mu_{k+1},\Sigma_{k+1})$，即限制在 $\Omega$ 上的截断高斯。这是一个典型的「交叉熵法」步骤：用当前提议下最极端的一小撮样本去定义下一轮提议。反过来，当 $N_\eta\ge N_p$ 时说明失效样本已经足够多，循环终止。

终止时最终提议改用**先验加权**均值：

$$
\mu_{opt}=\frac{\sum_{i=1}^{N_p}\widetilde x_i\,\omega(\widetilde x_i)}{\sum_{i=1}^{N_p}\omega(\widetilde x_i)},
\qquad
\Sigma_{opt}=\frac{1}{N_p-1}\sum_{i=1}^{N_p}
\bigl(\widetilde x_i-\mu_{opt}\bigr)\otimes\bigl(\widetilde x_i-\mu_{opt}\bigr).
$$

失效概率的重要采样估计为

$$
\hat P^{SAIS}_{\mathcal F}=\frac{1}{N_2}\sum_{i=1}^{N_2}
\frac{\omega(x_i)}{\hat h_{opt}(x_i)}\,\mathbb I_{\Omega_{\mathcal F}}(x_i),
\qquad x_i\sim\hat h_{opt}.
$$

同一批样本里落在 $\Omega_{\mathcal F}$ 的那些直接充当新配点，因此估计指标与生成点是**同一次采样的两个输出**。论文报告 $p_0=0.1$ 使这一过程快速自终止且精度良好。论文的 Remark 1 补充说，若 $\Omega$ 无界，中间提议就直接取（不截断的）高斯，与先验 $\omega$ 的形式一致；本文只处理截断高斯的情形。

> [!note] 中间更新与最终更新不对称
> 中间步的 $\mu_{k+1}$ 是**不加权**平均，最终步的 $\mu_{opt}$ 是 $\omega$ **加权**平均。中间步只需把提议推向失效区，最终步则要让提议逼近零方差最优密度 $\mathbb I_{\Omega_{\mathcal F}}\omega/P_{\mathcal F}$，后者带先验因子，因此加权是自然的。

### 定理

**假设 4.1（稳定性）.** 存在与 $v$ 无关的常数 $C_1,C_2>0$，使 $\|v\|$ 与 $\|\mathcal Av\|_{2,\Omega}+\|\mathcal Bv\|_{2,\partial\Omega}$ 之间成立双边界。

**假设 4.2.** 边界残差满足 $\|\mathcal B(u-u(\cdot;\theta^\ast))\|_{2,\partial\Omega}\le\epsilon_b$。

**假设 4.3.** 残差在闭域上有界，$M:=\max_{x\in\Omega}|r(x;\theta^\ast)|<\infty$。

**Theorem 4.4.** 设 $\Omega$ 有界，FI-PINNs 的解 $u(x;\theta^\ast)$ 满足假设 4.1–4.3，则

$$
\bigl\|u(x)-u(x;\theta^\ast)\bigr\|_{2,\Omega}
\le\sqrt2\,C_1^{-1}\Bigl(S_\Omega\bigl(M^2\epsilon_p+\epsilon_r^2\bigr)+\epsilon_b^2\Bigr)^{1/2},
$$

其中 $S_\Omega$ 是 $\Omega$ 的面积，$\epsilon_r,\epsilon_p$ 是两个**预设**容差。

证明路径很短且值得记住：把残差的 $L^2$ 范数按安全集与失效集分开，

$$
\|r\|^2_{2,\Omega}=\int_{\Omega_{\mathcal F}}r^2+\int_{\Omega_{\mathcal S}}r^2 .
$$

失效集上残差可能大但面积小，$S_{\Omega_{\mathcal F}}<S_\Omega\epsilon_p$ 给出 $\int_{\Omega_{\mathcal F}}r^2\le M^2S_\Omega\epsilon_p$；安全集上按极限状态函数的定义有 $|r|<\epsilon_r$，给出 $\int_{\Omega_{\mathcal S}}r^2\le S_\Omega\epsilon_r^2$。于是 $\|r\|^2_{2,\Omega}\le S_\Omega(M^2\epsilon_p+\epsilon_r^2)$，再用假设 4.1 的下界把残差范数换成解的误差范数，用假设 4.2 处理边界项。两个预设容差因此直接进入误差界，这正是「失效概率是后验误差指标」这句话的定量含义，而不是一个类比。

> [!warning] 证明中的一处先验因子
> 证明步骤把 $P_{\mathcal F}=\int_\Omega\mathbb I_{\Omega_{\mathcal F}}(x)\,\mathrm dx=S_{\Omega_{\mathcal F}}/S_\Omega$，省去了定义式中的 $\omega(x)$。两者仅在均匀先验 $\omega\equiv1/S_\Omega$ 时一致。对非均匀先验，$S_{\Omega_{\mathcal F}}<S_\Omega\epsilon_p$ 这一步需要额外条件。

### 数值实验

外层循环（Algorithm 1）的输入是：网络解 $u(x;\theta)$、边界点集 $\mathcal D_b$、配点集 $\mathcal D_c$、最大轮数 $M$、残差容差 $\epsilon_r$、失效概率容差 $\epsilon_p$。每轮先在当前集合上训练，再用自适应重要采样（Algorithm 2）同时得到 $\hat P_{\mathcal F}$ 与新点集 $\mathcal D_{adaptive}$；若 $\hat P_{\mathcal F}<\epsilon_p$ 则停止，否则 $\mathcal D_c\leftarrow\mathcal D_c\cup\mathcal D_{adaptive}$ 并继续。**训练集单调增长**，这是第二部分要改的地方。

第 5 节的算例覆盖五类问题，都是挑选出来让均匀采样吃亏的：

| 算例                          | 为什么难                         |
| ----------------------------- | -------------------------------- |
| 二维 Poisson，解带尖峰／奇性  | 高残差区只占定义域极小一部分     |
| Burgers 方程                  | 解含随时间变陡的内层             |
| 高维 Poisson 问题             | 均匀候选在高维几乎不覆盖有效区域 |
| 无界二维区域上的 Poisson 问题 | 无界域上均匀采样无从定义         |
| 无界域上的时间依赖问题        | 前两条同时出现                   |

基线为均匀采样与残差自适应细化。定性结论是自适应重要采样把点集中到高残差区，误差衰减快于均匀采样与残差细化。作者所在小组的公开仓库为 [SEU-YL-UQ/FI-PINNs](https://github.com/SEU-YL-UQ/FI-PINNs)；论文的 arXiv 渲染本身没有印出这个地址，因此宜按「作者小组仓库」引用。

> [!warning] 数值结果的可核实范围
> 算例的类型、基线与定性结论已核实；**逐例的误差量级，以及「高维 Poisson」到底是几维，未能从可访问的渲染中确认**，因此本页不给数字，也不给维数。

### 与其他论文的关系

这是三部曲的第一部，由编号 73（第二部分：重采样加子集模拟）与编号 76（第三部分：反问题，截断高斯**混合**）接续。它与密度驱动的自适应（编号 72、80、87）构成对照：FI-PINNs 在**残差**大的地方加点，而 ADDA/B-KRnet 一支从**解密度本身**采样。论文自己把基于 KRnet 的生成式采样列为第三类自适应策略，与本文并列。

## 73：定量训练集与子集模拟

### 直觉

第一部分留下两个可见的缺口。其一，训练集单调增长，每一轮都比上一轮贵；到后期，早期加进来的点大多已经不再提供信息，却仍然参与每一次梯度计算。其二，失效概率用**单个**截断高斯估计——一个单峰、椭球形的提议。当失效区是多峰的，或者形状复杂（比如沿一条曲线分布），单峰提议要么漏掉一部分模态，要么为了覆盖全部而变得过宽，两种情况下重要采样权重的方差都会爆炸。

第二部分分别针对这两点：把「累加」换成「定量重采样」，把「单个截断高斯」换成「子集模拟」。

### 推导

**扩展一：余弦退火下的定量重采样。** 训练集规模保持不变，配点组成按余弦退火从均匀逐步过渡到自适应。这个安排的动机是训练早期网络还很差，此时的残差场不可信，按它加点等于放大噪声；随着训练推进，残差场才逐渐变成一个有意义的指标。退火因此把「什么时候可以相信自适应」写成了一条随训练进度变化的规则。

**扩展二：子集模拟作为后验模型。** 子集模拟把一个小概率写成一列**嵌套中间失效水平**上较大条件概率的乘积：取一列阈值 $\epsilon_r^{(1)}>\epsilon_r^{(2)}>\cdots$，每一层的条件概率都不太小，可以用适中的样本量估计，层内的样本由 MCMC 从上一层的样本出发生成。这样既得到 $P_{\mathcal F}$ 的估计，又得到失效区中的样本，而且不预设失效区的形状——这正是它相对单个截断高斯的优势所在。

### 数值实验

论文的实证主张是在若干「有挑战性的问题」上相对原始算法有显著改进；已发表的期刊摘要补充说，子集模拟后验模型「能更有效地估计失效概率并在失效区生成新的有效训练点」。

> [!warning] 可核实范围
> 余弦退火的具体公式、中间失效水平的构造与所用 MCMC 核，以及基准问题的清单与误差表，本页依据的核对材料未能逐式确认，因此这一节只报告构造的**作用与位置**，不给公式也不给数字。多峰失效区的说法应当读作动机，而不是论文报告的实验结果。

### 与其他论文的关系

三部曲的中间一部。与第一部分相比改的是**后验模型**（截断高斯 → 子集模拟）与**集合管理策略**（累加 → 重采样）；第三部分再一次更换后验模型（→ 截断高斯混合）并转向反问题。余弦退火式的「从均匀过渡到自适应」与 ADDA 一支（编号 64、72、80、87）中「先均匀、再从当前模型采样」的时间表在概念上是同一件事。

> [!note] 作者姓名
> 部分文献聚合记录把第三作者列作 "Yan Liang"，正确形式是 **Liang Yan**（东南大学），与编号 70、76 一致。

## 76：截断高斯混合与反问题

### 直觉

反问题把第一部分的困难放大了一档。此时有两个网络在联合训练——状态 $u(x;\theta_u)$ 与系数 $\gamma(x;\theta_\gamma)$——损失里除了方程残差与边界残差，还有一项观测失配。两个网络互相牵制，残差场因此更容易分裂成**若干互不连通的高残差区**：系数在某处估错，会让状态在另一处偏离。单峰提议在这种几何下是不够的。

第三部分的选择是把提议换成截断高斯**混合**，并用 EM 拟合。混合模型天然处理多模态，而 EM 是拟合它的标准工具，代价可控。

### 问题设定

反问题写成

$$
\mathcal A[x;u(x),\gamma(x)]=0\ \text{在}\ \Omega,
\qquad
\mathcal B[x;u(x),\gamma(x)]=0\ \text{在}\ \partial\Omega,
$$

带间接观测 $y(x)=\mathcal G[x;u(x),\gamma(x)]$，$x\in D_{\mathrm{indirect}}$。两个网络在

$$
\mathcal L(\theta)=\mathcal L_c(\theta)+\lambda\mathcal L_b(\theta)+\mu\mathcal L_d(\theta),
\qquad
\mathcal L_c=\|r(x;\theta)\|_R^2,\quad
\mathcal L_b=\|b(x;\theta)\|_B^2,\quad
\mathcal L_d=|d(\theta)|_D^2
$$

上训练，其中 $d(x_i;\theta)=y(x_i)-\mathcal G[x_i;u(x_i,\theta_u),\gamma(x_i,\theta_\gamma)]$。极限状态函数相应地写成 $g(x)=|r(x;\theta)|_{\tilde R}-\varepsilon_r$，$\tilde R$ 是进入 $R$ 范数的各项绝对值之和；对下面的电阻抗成像算例它退化为 $|r(x;\theta)|$。

### 推导

本文明确写出零方差最优提议

$$
h_{\mathrm{opt}}(x)=\frac{\mathbb I_{\Omega_{\mathcal F}}(x)\,\omega(x)}{P_{\mathcal F}}
=\frac{\mathbb I_{\{g(x)>0\}}(x)\,\omega(x)}{\int_\Omega\mathbb I_{\{g(x)>0\}}(x)\omega(x)\,\mathrm dx},
$$

并用截断高斯混合去逼近它：

$$
h_{\mathrm{opt}}\approx h(x;\eta)=\sum_{k=1}^{K}\pi_k\,\mathcal N(x;\mu_k,\Sigma_k),
\qquad \pi_k\ge0,\ \sum_k\pi_k=1 .
$$

拟合准则是 $\min_\eta D_{\mathrm{KL}}(h_{\mathrm{opt}}\|h(\cdot;\eta))$。展开 KL 散度，只有交叉熵一项依赖 $\eta$，因此离散化后化为在 $x_j\sim h_{\mathrm{opt}}$ 上的极大似然

$$
\max_\eta\ \log\prod_{j=1}^{N_c}h(x_j,\eta),
$$

而 $x_j\sim h_{\mathrm{opt}}$ 只需从当前提议中取落在失效集里的样本。参数由 EM 迭代给出，E 步计算责任度

$$
q^{(t)}_{k,j}=\frac{\pi^{(t)}_k\,\mathcal N\bigl(x_j;\mu^{(t)}_k,\Sigma^{(t)}_k\bigr)}
{\sum_{k=1}^{K}\pi^{(t)}_k\,\mathcal N\bigl(x_j;\mu^{(t)}_k,\Sigma^{(t)}_k\bigr)},
$$

M 步更新三组参数

$$
\pi^{(t+1)}_k=\frac1{N_c}\sum_{j}q^{(t)}_{k,j},
\qquad
\mu^{(t+1)}_k=\frac{\sum_j q^{(t)}_{k,j}x_j}{\sum_j q^{(t)}_{k,j}},
\qquad
\Sigma^{(t+1)}_k=\frac{\sum_j q^{(t)}_{k,j}\bigl(x_j-\mu^{(t+1)}_k\bigr)\bigl(x_j-\mu^{(t+1)}_k\bigr)^{\!\top}}{\sum_j q^{(t)}_{k,j}} .
$$

实现直接调用 `sklearn.mixture`（论文 Remark 3.1）。

截断由投影实现，而不是由带截断的密度公式：$\mathrm{Proj}_\Omega(x)=\arg\min_{y\in\bar\Omega}\|x-y\|_2$。第 $l$ 轮的 EM 更新写成 $\eta^{\star,l+1}=\mathrm{EM}(\mathrm{Proj}_\Omega(x^l),K^l)$，最终拟合用极限状态函数值最大的 $N_p$ 个失效点。重要采样估计相应写成

$$
\hat P^{\mathrm{SAIS}}_{\mathcal F}=\frac1{N_2}\sum_{i=1}^{N_2}
\frac{\omega\bigl(\mathrm{Proj}_\Omega(x_i)\bigr)}{\hat h_{\mathrm{opt}}\bigl(\mathrm{Proj}_\Omega(x_i)\bigr)}
\,\mathbb I_{\Omega_{\mathcal F}}\bigl(\mathrm{Proj}_\Omega(x_i)\bigr).
$$

### 定理

本文不给新的收敛定理。它的贡献是提议模型加拟合路线：与另一类做法（用风险函数极大化拟合混合、在高维下被迫把 $\Sigma_k$ 对角化）相比，EM 允许完整的协方差。Remark 3.3 给出代价判断：EM 增加的时间相对网络训练可以忽略。

### 数值实验

两类反问题。第一类是**电阻抗成像中的反导率问题**：

$$
-\nabla\!\cdot\!(\gamma\nabla u)-f=0\ \text{在}\ \Omega,
\qquad
u=u_b,\ \ \gamma=\gamma_b,\ \ \partial_{\vec n}u=u_n\ \text{在}\ \partial\Omega,
$$

损失为 $\mathcal L_c+\lambda_{u_b}\mathcal L_{u_b}+\lambda_{\gamma_b}\mathcal L_{\gamma_b}+\lambda_{\partial u}\mathcal L_{\partial u}$，其中 $\mathcal L_c=\|-\nabla\cdot(\gamma\nabla u)-f\|^2_{L^2(\Omega)}$，三个边界项都是 $L^2(\partial\Omega)$ 范数。加权范数定义为 $\|f\|^2_{L^2(\Omega)}=\int_\Omega f^2(x)\omega(x)\,\mathrm dx$，本例取 $\omega\equiv1$——这一点值得留意，因为上面那条关于先验因子的告诫在 $\omega\equiv1$ 时不构成问题。第二类是**抛物系统中的反源问题**。

> [!warning] 数值结果的可核实范围
> 两类算例的设定与损失结构已按期刊 PDF 核实；**逐例的误差表未转录**，因此本页不给数字。

### 与其他论文的关系

三部曲的第三部，也是唯一面向反问题（两个耦合网络）的一部。三部的差别几乎完全集中在估计 $P_{\mathcal F}$ 并生成新点的**后验模型**：单个截断高斯（70）→ 子集模拟（73）→ EM 拟合的截断高斯混合（76）。「用 KL 散度逼近最优密度」这一变分原理，与编号 80 和 87 用正规化流代替混合模型时用的是同一条。

## 80：变分损失没有残差，于是改用方差

### 直觉

深度 Ritz 方法最小化的是变分能量而不是逐点残差，于是前面整套机器都用不上：没有 $r(x;\theta)$，就没有极限状态函数，也没有失效集。

论文换了一个角度看问题。深度 Ritz 的损失是一个积分，而这个积分是用 Monte Carlo 求积近似的，因此离散误差天然分成两块：网络能不能表示真解（逼近误差），以及有限个样本能不能把积分算准（**统计误差**）。第二块在被积函数正则性低时会主导。论文用一个可以手算的例子把这句话变成数字：取 $G(u(x))=\frac{1}{\sqrt{2\pi}}e^{-x^2/2\sigma^2}$ 放在 $[-1,1]$ 上，相对 Monte Carlo 误差是 $C(\sigma N)^{-1/2}$，即要达到 $O(1)$ 的相对精度需要 $O(1/\sigma)$ 个均匀样本。当 $\sigma$ 很小时——也就是被积函数集中在一个窄峰上时——均匀采样是灾难性的。

因此自适应的目标应当是**变分损失的方差缩减**，而重要采样正是为此设计的工具：把被积函数当成未归一化密度，学它，然后按它采样。

### 问题设定

带罚的深度 Ritz 问题为

$$
\min_{u\in V} J(u),
\qquad
J(u)=\int_\Omega G\bigl(u(x)\bigr)\,\mathrm dx
+\beta\,\bigl\|B\bigl(x,u(x)\bigr)\bigr\|^2_{\partial\Omega,2},
$$

在均匀配点集 $S_\Omega=\{x_i\}_{i=1}^{N_v}$、$S_{\partial\Omega}=\{x_j\}_{j=1}^{N_b}$ 上离散为

$$
J_N\bigl(u(x,\theta)\bigr)=\frac1{N_v}\sum_{i=1}^{N_v}G\bigl(u(x_i,\theta)\bigr)
+\frac{\beta}{N_b}\sum_{j=1}^{N_b}B\bigl(x_j,u(x_j,\theta)\bigr)^2 .
$$

重要采样估计为

$$
I(u)=\int_\Omega G\bigl(u(x,\theta)\bigr)\mathrm dx
=\mathbb E_p\Bigl[\frac{G(u(X,\theta))}{p(X)}\Bigr]
\approx\frac1{N_v}\sum_{i=1}^{N_v}\frac{G(u(x_i,\theta))}{p(x_i)},
\qquad x_i\sim p .
$$

### 推导

**变号被积函数：自适应无法把方差降到零。** 若 $G\ge0$，最优提议 $p^\star=G/\mu$（$\mu=\int_\Omega G$）使被积量 $G/p$ 恒等于常数 $\mu$，方差为零。但变分能量的被积函数不必非负——Dirichlet 能量减去源项那一部分本来就会变号。此时最优选择是

$$
p^\star(x)=\frac{|G(u(x,\theta))|}{\mu},
\qquad \mu=\int_\Omega|G(u(x,\theta))|\,\mathrm dx,
$$

由 Cauchy-Schwarz 它在所有密度中方差最小，但残余方差严格为正：

$$
\sigma_{p^\star}=\Bigl(\int_\Omega|G(u(x,\theta))|\mathrm dx\Bigr)^2
-\Bigl(\int_\Omega G(u(x,\theta))\mathrm dx\Bigr)^2>0 .
$$

右端是 $|G|$ 的积分平方减去 $G$ 的积分平方，只有在 $G$ 不变号时才为零。这是一个结构性结论：**被积函数变号时，仅靠自适应不能把统计误差消掉，还必须增加样本量** $N_v$，因为 Monte Carlo 误差是 $\sigma_p/\sqrt{N_v}\ge\sigma_{p^\star}/\sqrt{N_v}>0$。这一点把「自适应采样」的作用界定得很清楚。

**两个网络交替更新。** 密度模型取有界 KRnet（见[[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|密度流一页]]）：

$$
p_{\text{bKRnet}}(x,\theta_f)=p_Z\bigl(f_{\text{bKRnet}}(x,\theta_f)\bigr)\,
\bigl|\det\nabla_x f_{\text{bKRnet}}\bigr|,
\qquad Z\sim\mathrm{Unif}([-1,1]^d),
$$

采样由逆流 $X=f^{-1}_{\text{bKRnet}}(Z)$ 给出。训练准则 $\min_{\theta_f}D_{\mathrm{KL}}(p^\star\|p_{\text{bKRnet}})$ 中只有交叉熵项依赖 $\theta_f$，因此归结为极小化 $H(p^\star,p_{\text{bKRnet}})$；由于 $p^\star$ 只能按未归一化形式 $|G|$ 求值，交叉熵本身要用重要采样估计：

$$
H(p^\star,p_{\text{bKRnet}})\approx
-\frac1{N_v}\sum_{i=1}^{N_v}
\frac{|G(u(x_i,\theta))|\,\log p_{\text{bKRnet}}(x_i,\theta_f)}{\mu\,\tilde p(x_i)},
\qquad x_i\sim\tilde p .
$$

两个网络在第 $k$ 轮交替更新，解网络在当前密度下做重要采样加权的 Ritz 最小化，密度网络在新解下做交叉熵最小化：

$$
\theta^{k+1}=\arg\min_\theta\ \frac1{N_v}\sum_{i}
\frac{G\bigl(u(x^k_{\Omega,i},\theta)\bigr)}{p_{\text{bKRnet}}(x^k_{\Omega,i},\theta_f^k)}
+\frac{\beta}{N_b}\sum_{j}B^2\bigl(x^k_{\partial\Omega,j},u(x^k_{\partial\Omega,j},\theta)\bigr),
$$

$$
\theta_f^{k+1}=\arg\min_{\theta_f}\
-\frac1{N_B}\sum_{l}
\frac{\bigl|G\bigl(u(x^k_{B,l},\theta^{k+1})\bigr)\bigr|\,
\log p_{\text{bKRnet}}(x^k_{B,l},\theta_f)}
{\mu^{k+1}\,p_{\text{bKRnet}}(x^k_{B,l},\theta_f^k)} .
$$

初始密度 $p_{\text{bKRnet}}(\cdot,\theta_f^0)$ 取均匀，$\mu^{k+1}$ 作为常数可以省去。

**两个混合模型给出密度下界。** 比值 $G/p$ 在学到的密度很小处会爆炸，因此论文给密度加一个下界。模型一与均匀分布混合：

$$
p_{\mathrm{mixture}}(x,\theta_f^{k+1})
=\epsilon\,p_{\text{bKRnet}}(x,\theta_f^{k+1})+(1-\epsilon)\,p_{\mathrm{uniform}}(x)
\ \ge\ \frac{1-\epsilon}{|\Omega|} .
$$

模型二与**历史上所有**流递归混合：

$$
p_{\mathrm{mixture}}\bigl(x,\{\theta_f^t\}_{t\le k+1}\bigr)
=\epsilon\,p_{\text{bKRnet}}(x,\theta_f^{k+1})
+(1-\epsilon)\,p_{\mathrm{mixture}}\bigl(x,\{\theta_f^t\}_{t\le k}\bigr)
=\sum_{t=1}^{k+1}\epsilon(1-\epsilon)^{k+1-t}p_{\text{bKRnet}}(x,\theta_f^t)
+\frac{(1-\epsilon)^{k+1}}{|\Omega|},
$$

下界为 $(1-\epsilon)^{k+1}/|\Omega|$，随轮数指数变小——这是模型二的代价，它换来的是提议覆盖历史上所有被认为重要的区域。

### 定理

本文没有关于耦合自适应迭代的收敛定理。它的理论内容是方差序：由 Cauchy-Schwarz，对任意密度 $p$ 有 $\sigma_{p^\star}\le\sigma_p$；配上变号情形下 $\sigma_{p^\star}>0$ 的显式表达，得到 Monte Carlo 误差的下界

$$
\mathrm{Var}_p^{1/2}[\cdot]=\frac{\sigma_p}{\sqrt{N_v}}\ \ge\ \frac{\sigma_{p^\star}}{\sqrt{N_v}}>0 .
$$

### 数值实验

Algorithm 1（基本版）的结构是：输入初始解网络 $u(x,\theta^0)$、初始密度 $p_{\text{bKRnet}}(x,\theta_f^0)$、最大轮数 $N_e$、批量 $m$，以及三个初始训练集 $S^0_\Omega$、$S^0_{\partial\Omega}$、$S^0_B$；外层对 $k=0,\dots,N_{\mathrm{adaptive}}-1$ 循环，内层先用随机梯度法解 PDE，再用随机梯度法训练有界 KRnet，最后从 $p_{\text{bKRnet}}(x,\theta_f^{k+1})$ 采出新的 $S^{k+1}_\Omega$ 与 $S^{k+1}_B$。Algorithm 2 把新集合改成均匀点与流样本的并集并用混合模型一算密度值；Algorithm 3 保存**每一轮**的 $\theta_f^t$，从全部流中重新采点，并用混合模型二算密度值。

四个算例与四种采样策略并列比较：

| 算例      | 问题                     |
| --------- | ------------------------ |
| 第 4.1 节 | 二维单峰问题             |
| 第 4.2 节 | 二维双峰问题             |
| 第 4.3 节 | 二维带奇性的问题         |
| 第 4.4 节 | 高维 Poisson 问题        |

四种策略是：常规深度 Ritz（每轮重新均匀采点）、Algorithm 1、Algorithm 2、Algorithm 3。共同设定为：误差用张量网格上的相对离散 $L_2$ 误差（高维改用均匀样本），激活取 $\sin^3(x)$ 的类 ResNet 结构，两个网络都用 Adam，有界 KRnet 的架构与编号 87 相同。摘要的定量主张是：相对原始深度 Ritz，自适应方法提高精度，**尤其在低正则性与高维问题上**。

> [!warning] 数值结果的可核实范围
> 算例清单、误差度量、网络与优化设定、四种策略的对照设计均已核实；**逐例误差数字未能确认**，因此本页不给数字。

### 与其他论文的关系

编号 80 是 FI-PINNs 一支（编号 70、73、76）在深度 Ritz 上的对应物——同样的「学一个密度、重采样、再训练」循环，但目标密度是 $|G(u)|/\mu$（变分被积函数）而不是基于残差的失效指标。它的密度模型正是编号 87 的有界 KRnet，而编号 87 的结论部分明确把这个耦合列为后续工作。与均匀分布混合（模型一）这个保险措施，与编号 87 里的更新率 $\gamma$ 和编号 73 里余弦退火的均匀到自适应过渡是同一类安全阀。它还与编号 60 互补：后者修的是同一个深度 Ritz 目标里**边界条件的处理方式**。

## 五篇的共同结构

| 编号 | 被驱动的目标           | 提议 / 密度模型         | 训练集更新方式         |
| ---- | ---------------------- | ----------------------- | ---------------------- |
| 66   | 算子估计的无偏性       | 两组独立随机参数        | 不涉及（采样在算子内） |
| 70   | 残差失效概率           | 单个截断高斯            | 单调累加               |
| 73   | 残差失效概率           | 子集模拟                | 定量重采样 + 余弦退火  |
| 76   | 残差失效概率（含观测） | 截断高斯混合（EM 拟合） | 单调累加               |
| 80   | 变分损失的方差         | 有界 KRnet + 均匀混合   | 每轮从当前密度重采     |

三条可迁移的判断：

- **指标要能直接进入误差界。** 编号 70 的 Theorem 4.4 让两个预设容差出现在最终估计里，这是「失效概率是后验误差指标」的定量版本，而不是类比。
- **提议模型的表达力决定上限。** 从单个截断高斯到子集模拟再到混合模型，改进的对象一直是同一个量的估计方式。
- **自适应的作用有上界。** 编号 80 的残余方差说明，被积函数变号时样本量仍是独立的必要资源。

## 覆盖核对

| 内容                            | 论文 | 覆盖状态                                          |
| ------------------------------- | ---- | ------------------------------------------------- |
| 分数阶 Laplacian 的内外分裂估计 | 66   | 完整推导链、两个 Beta 采样、$r_\epsilon$ 保护     |
| Caputo 导数的期望表示           | 66   | 公式与采样分布（含前因子核对说明）                |
| 平方偏差与两组独立随机数        | 66   | 无偏性构造与其成立条件                            |
| 66 的数值实验                   | 66   | 算法与三组实验的设定；**误差量级未核实**          |
| 极限状态函数与失效概率          | 70   | 定义、安全/失效集、与残差细化的差别               |
| 自适应重要采样                  | 70   | 迭代更新、终止条件、加权最终提议、估计式          |
| Theorem 4.4 与其证明骨架        | 70   | 三条假设、误差界、分集估计、先验因子核对          |
| 70 的数值实验                   | 70   | 五类算例与基线；**误差量级与维数未核实**          |
| 定量重采样与子集模拟            | 73   | 作用与位置（限定可核实范围）                      |
| 零方差最优提议与截断高斯混合    | 76   | 最优密度、KL 准则、EM 更新、投影截断              |
| 76 的数值实验                   | 76   | 电阻抗成像方程与损失、反源问题；**误差表未核实**  |
| 变分损失的统计误差              | 80   | 误差分解、算例结论、$O(1/\sigma)$ 样本需求        |
| 变号被积函数的残余方差          | 80   | 最优密度、方差序与严格正下界的结构含义            |
| 两网络交替更新与两个混合模型    | 80   | 两个目标函数、两个密度下界                        |
| 80 的数值实验                   | 80   | 三个算法版本、四个算例、四种策略；**数字未核实**  |

## 本页原文

- L. Guo, H. Wu, X. Yu, and T. Zhou, _Monte Carlo fPINNs: deep learning method for forward and inverse problems involving high dimensional fractional partial differential equations_, Comput. Methods Appl. Mech. Engrg. 400 (2022), 115523（预印本 [arXiv:2203.08501](https://arxiv.org/abs/2203.08501)）。
- Z. Gao, L. Yan, and T. Zhou, [_Failure-informed adaptive sampling for PINNs_](https://doi.org/10.1137/22M1527763), SIAM J. Sci. Comput. 45(4) (2023), pp. A1971-A1994（预印本 [arXiv:2210.00279](https://arxiv.org/abs/2210.00279)）。
- Z. Gao, T. Tang, L. Yan, and T. Zhou, [_Failure-informed adaptive sampling for PINNs, Part II: combining with re-sampling and subset simulation_](https://doi.org/10.1007/s42967-023-00312-7), Commun. Appl. Math. Comput. 6 (2024), pp. 1720-1741（预印本 [arXiv:2302.01529](https://arxiv.org/abs/2302.01529)）。
- W. Liu, L. Yan, T. Zhou, and Y. Zhou, [_Failure-informed adaptive sampling for PINNs, Part III: applications to inverse problems_](https://doi.org/10.4208/csiam-am.SO-2023-0059), CSIAM Trans. Appl. Math. 5(3) (2024), pp. 636-670。
- X. Wan, T. Zhou, and Y. Zhou, _Adaptive importance sampling for deep Ritz_, Commun. Appl. Math. Comput. 7(3) (2025), pp. 929-953（预印本 [arXiv:2310.17185](https://arxiv.org/abs/2310.17185)）。
