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

在讨论采样点之前，编号 66 先解决另一个问题：分数阶算子的组装本身太贵。早先的 fPINN 在辅助网格上离散分数阶导数，这在低维之外无法扩展。本文的选择是**根本不离散化分数阶算子**，而把它写成期望，用 Monte Carlo 在损失内部估计。

模型方程为分数阶对流扩散：

$$
L[u(x,t)]:=\frac{\partial^\gamma u(x,t)}{\partial t^\gamma}
+c\,(-\Delta)^{\alpha/2}u(x,t)+v\cdot\nabla u(x,t)=f(x,t),
$$

其中积分型分数阶 Laplacian 为

$$
(-\Delta)^{\alpha/2}u(x)\triangleq C_{d,\alpha}\,\mathrm{P.V.}\!
\int_{\mathbb R^d}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\,\mathrm dy .
$$

### 内外分裂与两个 Beta 分布

在 $x$ 周围以半径 $r_0$ 分裂积分。球内部分，取方向 $\xi\sim\mathrm{Unif}(S^{d-1})$、半径 $r\sim f_I$：

$$
\int_{y\in B_{r_0}(x)}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\mathrm dy
=\frac{|S^{d-1}|\,r_0^{2-\alpha}}{2(2-\alpha)}\;
\mathbb E_{\xi,\,r\sim f_I}\!\left[\frac{2u(x)-u(x-r\xi)-u(x+r\xi)}{r^2}\right],
$$

$$
f_I(r)=\frac{2-\alpha}{r_0^{2-\alpha}}\,r^{1-\alpha}\mathbf 1_{r\in[0,r_0]},
\qquad r/r_0\sim\mathrm{Beta}(2-\alpha,1).
$$

球外部分，$r\sim f_O$：

$$
\int_{y\notin B_{r_0}(x)}\frac{u(x)-u(y)}{\|x-y\|_2^{d+\alpha}}\mathrm dy
=\frac{|S^{d-1}|\,r_0^{-\alpha}}{2\alpha}\;
\mathbb E_{\xi,\,r\sim f_O}\bigl[2u(x)-u(x-r\xi)-u(x+r\xi)\bigr],
$$

$$
f_O(r)=\alpha\,r_0^{\alpha}\,r^{-1-\alpha}\mathbf 1_{r\in[r_0,\infty)},
\qquad r_0/r\sim\mathrm{Beta}(\alpha,1).
$$

这两个恒等式的意义在于，被积函数的奇性被吸收进采样密度，剩下的期望里只有二阶差分。由于

$$
\lim_{r\to0}\frac{2u(x)-u(x-r\xi)-u(x+r\xi)}{r^2}
=\partial^2_r u(x+r\xi)\big|_{r=0},
$$

数值上仍需一个下界保护 $r_\epsilon=\max\{\epsilon,r_I\}$，防止 $1/r^2$ 放大舍入误差。

Caputo 时间导数用同一手法处理：

$$
\frac{\partial^\gamma u(x,t)}{\partial t^\gamma}
=\frac{\gamma}{1-\gamma}\,t^{1-\gamma}\,
\mathbb E_{\tau\sim f_{I,t}}\!\left[\frac{u(x,t)-u(x,t-\tau)}{\tau}\right]
+\frac{u(x,t)-u(x,0)}{t^{\gamma}},
\qquad
f_{I,t}(\tau)=(1-\gamma)\tau^{-\gamma}\mathbf 1_{\tau\in[0,1]} .
$$

### 无偏性：两组独立随机数

这是本文摘要点明的关键成分。若把带噪的算子估计直接平方，残差损失会被估计量的方差**系统性抬高**：$\mathbb E[\widehat L^2]=(\mathbb E\widehat L)^2+\mathrm{Var}(\widehat L)$。解决办法是对每个配点抽**两组独立**随机参数，把两个独立估计相乘而不是平方一个估计：

$$
\hat L_{equ}(\theta)=\frac{1}{mN_u}\sum_{i,j}
\widehat L\bigl[u_{NN}(x_i,t_i;\theta);\tau_j,\xi_j,r_{Ij},r_{oj}\bigr]\cdot
\widehat L\bigl[u_{NN}(x_i,t_i;\theta);\tau'_j,\xi'_j,r'_{Ij},r'_{oj}\bigr].
$$

在 $\epsilon=\epsilon_t=0$ 且忽略舍入的前提下，$\mathbb E[\hat L_{equ}(\theta)]=L_{equ}(\theta)$。

实验分三组：高维积分型分数阶 Laplacian 方程、时空分数阶 PDE 的参数辨识（反问题）、以及带随机输入的分数阶扩散方程。

> [!note] 题名与原文核对
> 预印本题为 "Monte Carlo PINNs: deep learning approach for forward and inverse problems involving high dimensional fractional partial differential equations"，期刊版为 "Monte Carlo fPINNs: Deep learning method for ..."。另外，Caputo 导数定义的前因子在预印本渲染中写作 $\Gamma(1-\alpha)$，而指数用的是 $\gamma$；按 $0<\gamma<1$ 的设定应为 $\Gamma(1-\gamma)$。引用前因子时建议回查期刊版。

## 70：把「残差过大」定义成失效事件

### 为什么残差最大的点不够用

标准的残差自适应细化在均匀候选池中取残差最大的 $m$ 个点。当高残差区只占定义域很小一部分时，均匀候选几乎落不进去，这个策略就失效了。编号 70 换了框架：把配点放置视为**可靠性分析**问题。

物理信息神经网络的损失是

$$
\mathcal L(\theta)=\mathcal L_c(\theta)+\lambda\mathcal L_b(\theta),
\qquad
\mathcal L_c(\theta)=\frac{1}{N_c}\sum_{i=1}^{N_c}\bigl|r(x^c_i;\theta)\bigr|^2,
\qquad
\mathcal L_b(\theta)=\frac{1}{N_b}\sum_{i=1}^{N_b}\bigl|b(x^b_i;\theta)\bigr|^2 .
$$

定义极限状态函数

$$
g(x)=\bigl|r(x;\theta)\bigr|-\epsilon_r,
$$

它的零水平集把定义域分成安全集 $\Omega_{\mathcal S}=\{g<0\}$ 与失效集 $\Omega_{\mathcal F}=\{g>0\}$。在先验 $\omega(x)$ 下的失效概率

$$
P_{\mathcal F}=\int_\Omega \omega(x)\,\mathbb I_{\Omega_{\mathcal F}}(x)\,\mathrm dx
$$

充当**后验误差指标**，其地位与自适应有限元中的后验误差估计完全对应：当 $P_{\mathcal F}<\epsilon_p$ 时宣布网络可靠。

与残差自适应细化的一个具体差别值得指出：基于 Monte Carlo 的加点会把**所有**落入 $\Omega_{\mathcal F}$ 的候选点都加进去，因此每轮加点个数是变化的；残差细化每轮固定加 $m$ 个。

### 自适应重要采样

失效概率小时朴素 Monte Carlo 估计不可用。论文用自适应重要采样逐步把提议密度推向失效区：从 $h_1=\omega$ 出发，第 $k$ 步抽 $N_1$ 个样本，按极限状态函数**降序**排列得到 $\widetilde x^k_1,\dots,\widetilde x^k_{N_1}$，令

$$
N_\eta=\max_{1\le i\le N_1}\{i:\ g(\widetilde x_i)>0\},
\qquad
N_p=\lfloor p_0N_1\rfloor .
$$

若 $N_\eta<N_p$（失效样本还太少），用前 $N_p$ 个样本更新截断高斯提议：

$$
\mu_{k+1}=\frac{1}{N_p}\sum_{i=1}^{N_p}\widetilde x^k_i,
\qquad
\Sigma_{k+1}=\frac{1}{N_p-1}\sum_{i=1}^{N_p}
\bigl(\widetilde x^k_i-\mu_{k+1}\bigr)\otimes\bigl(\widetilde x^k_i-\mu_{k+1}\bigr),
$$

并令 $h_{k+1}=\mathcal N_T(\mu_{k+1},\Sigma_{k+1})$。终止时最终提议改用**先验加权**均值：

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

论文报告 $p_0=0.1$ 使这一过程快速自终止且精度良好。

> [!note] 中间更新与最终更新不对称
> 中间步的 $\mu_{k+1}$ 是**不加权**平均，最终步的 $\mu_{opt}$ 是 $\omega$ **加权**平均。中间步只需把提议推向失效区，最终步则要让提议逼近零方差最优密度 $\mathbb I_{\Omega_{\mathcal F}}\omega/P_{\mathcal F}$，后者带先验因子，因此加权是自然的。

### Theorem 4.4：把两个容差直接写进误差界

在三条假设下——(4.1) 存在与 $v$ 无关的 $C_1,C_2>0$ 给出 $\|v\|$ 与 $\|\mathcal Av\|_{2,\Omega}+\|\mathcal Bv\|_{2,\partial\Omega}$ 的双边界；(4.2) 边界残差 $\|\mathcal B(u-u(\cdot;\theta^\ast))\|_{2,\partial\Omega}\le\epsilon_b$；(4.3) $M:=\max_{x\in\Omega}|r(x;\theta^\ast)|<\infty$——对有界 $\Omega$ 有

$$
\bigl\|u(x)-u(x;\theta^\ast)\bigr\|_{2,\Omega}
\le\sqrt2\,C_1^{-1}\Bigl(S_\Omega\bigl(M^2\epsilon_p+\epsilon_r^2\bigr)+\epsilon_b^2\Bigr)^{1/2},
$$

其中 $S_\Omega$ 是 $\Omega$ 的面积。

证明路径很短且值得记住：把残差的 $L^2$ 范数按安全集与失效集分开，

$$
\|r\|^2_{2,\Omega}=\int_{\Omega_{\mathcal F}}r^2+\int_{\Omega_{\mathcal S}}r^2 .
$$

失效集上残差可能大但面积小，$S_{\Omega_{\mathcal F}}<S_\Omega\epsilon_p$ 给出 $\int_{\Omega_{\mathcal F}}r^2\le M^2S_\Omega\epsilon_p$；安全集上 $|r|<\epsilon_r$ 给出 $\int_{\Omega_{\mathcal S}}r^2\le S_\Omega\epsilon_r^2$。两个预设容差 $\epsilon_p,\epsilon_r$ 因此直接进入误差界，这正是「失效概率是后验误差指标」这句话的定量含义。

> [!warning] 证明中的一处先验因子
> 证明步骤把 $P_{\mathcal F}=\int_\Omega\mathbb I_{\Omega_{\mathcal F}}(x)\,\mathrm dx=S_{\Omega_{\mathcal F}}/S_\Omega$，省去了定义式中的 $\omega(x)$。两者仅在均匀先验 $\omega\equiv1/S_\Omega$ 时一致。对非均匀先验，$S_{\Omega_{\mathcal F}}<S_\Omega\epsilon_p$ 这一步需要额外条件。

实验覆盖带尖峰或奇性解的二维 Poisson 方程、Burgers 方程、高维 Poisson 问题、无界二维区域上的 Poisson 问题以及一个无界区域上的时间依赖问题，基线为均匀采样与残差自适应细化。作者所在小组的公开仓库为 [SEU-YL-UQ/FI-PINNs](https://github.com/SEU-YL-UQ/FI-PINNs)。

## 73 与 76：换掉后验模型

三部曲的差别几乎完全集中在**用什么模型估计失效概率并生成新点**。

### 73：定量训练集与子集模拟

编号 70 的训练集单调增长（$\mathcal D_c\leftarrow\mathcal D_c\cup\mathcal D_{adaptive}$），训练代价随轮数上升；而单个截断高斯在失效区多峰或形状复杂时是很粗的提议。第二部分同时处理这两点：训练集规模保持不变，配点组成按余弦退火从均匀逐步过渡到自适应；失效概率估计改用子集模拟，即把一个小概率写成一列嵌套中间失效水平上较大条件概率的乘积，每层内部用 MCMC 采样。

> [!note] 可核实范围
> 余弦退火的具体公式、中间失效水平的构造与所用 MCMC 核，本页依据的公开材料未能逐式确认，因此只报告其作用与位置。

### 76：截断高斯混合与反问题

第三部分面向**反问题**：状态与系数两个网络联合训练，残差可能集中在若干互不连通的区域，单峰提议不够。损失为

$$
\mathcal L(\theta)=\mathcal L_c(\theta)+\lambda\mathcal L_b(\theta)+\mu\mathcal L_d(\theta),
$$

三项分别是方程残差、边界残差与观测残差 $d(x_i;\theta)=y(x_i)-\mathcal G[x_i;u(x_i,\theta_u),\gamma(x_i,\theta_\gamma)]$。

本文明确写出零方差最优提议

$$
h_{\mathrm{opt}}(x)=\frac{\mathbb I_{\Omega_{\mathcal F}}(x)\,\omega(x)}{P_{\mathcal F}}
=\frac{\mathbb I_{\{g(x)>0\}}(x)\,\omega(x)}{\int_\Omega\mathbb I_{\{g(x)>0\}}(x)\omega(x)\,\mathrm dx},
$$

并用截断高斯混合去逼近它：

$$
h_{\mathrm{opt}}\approx h(x;\eta)=\sum_{k=1}^{K}\pi_k\,\mathcal N(x;\mu_k,\Sigma_k),
\qquad \pi_k\ge0,\ \sum_k\pi_k=1,
$$

拟合准则是 $\min_\eta D_{\mathrm{KL}}(h_{\mathrm{opt}}\|h(\cdot;\eta))$，离散化后化为在 $x_j\sim h_{\mathrm{opt}}$ 上的极大似然 $\max_\eta\log\prod_{j=1}^{N_c}h(x_j,\eta)$。参数由 EM 迭代给出：

$$
q^{(t)}_{k,j}=\frac{\pi^{(t)}_k\,\mathcal N\bigl(x_j;\mu^{(t)}_k,\Sigma^{(t)}_k\bigr)}
{\sum_{k=1}^{K}\pi^{(t)}_k\,\mathcal N\bigl(x_j;\mu^{(t)}_k,\Sigma^{(t)}_k\bigr)},
$$

$$
\pi^{(t+1)}_k=\frac1{N_c}\sum_{j}q^{(t)}_{k,j},
\qquad
\mu^{(t+1)}_k=\frac{\sum_j q^{(t)}_{k,j}x_j}{\sum_j q^{(t)}_{k,j}},
\qquad
\Sigma^{(t+1)}_k=\frac{\sum_j q^{(t)}_{k,j}\bigl(x_j-\mu^{(t+1)}_k\bigr)\bigl(x_j-\mu^{(t+1)}_k\bigr)^{\!\top}}{\sum_j q^{(t)}_{k,j}} .
$$

截断由投影实现：$\mathrm{Proj}_\Omega(x)=\arg\min_{y\in\bar\Omega}\|x-y\|_2$，重要采样估计相应写成

$$
\hat P^{\mathrm{SAIS}}_{\mathcal F}=\frac1{N_2}\sum_{i=1}^{N_2}
\frac{\omega\bigl(\mathrm{Proj}_\Omega(x_i)\bigr)}{\hat h_{\mathrm{opt}}\bigl(\mathrm{Proj}_\Omega(x_i)\bigr)}
\,\mathbb I_{\Omega_{\mathcal F}}\bigl(\mathrm{Proj}_\Omega(x_i)\bigr).
$$

论文把这条路线与另一类做法（用风险函数极大化拟合混合、在高维下被迫把 $\Sigma_k$ 对角化）作了对比，并指出 EM 增加的时间相对网络训练可以忽略。实验为电阻抗成像中的反导率问题与抛物系统中的反源问题。

## 80：变分损失没有残差，于是改用方差

### 为什么不能照搬

深度 Ritz 方法最小化的是变分能量而不是逐点残差，残差型自适应无从下手。论文指出，变分损失的离散误差分成逼近误差与**统计（Monte Carlo 求积）误差**两部分，而当被积函数正则性低时后者主导。论文的算例把 $G(u(x))=\frac{1}{\sqrt{2\pi}}e^{-x^2/2\sigma^2}$ 放在 $[-1,1]$ 上，相对 Monte Carlo 误差为 $C(\sigma N)^{-1/2}$，即要达到 $O(1)$ 的相对精度需要 $O(1/\sigma)$ 个均匀样本。因此自适应的目标应当是**变分损失的方差缩减**。

带罚的深度 Ritz 问题为

$$
J(u)=\int_\Omega G\bigl(u(x)\bigr)\,\mathrm dx
+\beta\,\bigl\|B\bigl(x,u(x)\bigr)\bigr\|^2_{\partial\Omega,2},
$$

重要采样估计为

$$
I(u)=\int_\Omega G\bigl(u(x,\theta)\bigr)\mathrm dx
=\mathbb E_p\Bigl[\frac{G(u(X,\theta))}{p(X)}\Bigr]
\approx\frac1{N_v}\sum_{i=1}^{N_v}\frac{G(u(x_i,\theta))}{p(x_i)},
\qquad x_i\sim p .
$$

### 变号被积函数：自适应无法把方差降到零

若 $G\ge0$，最优提议 $p^\star=G/\mu$（$\mu=\int_\Omega G$）给出零方差。但变分能量的被积函数不必非负。此时最优选择是

$$
p^\star(x)=\frac{|G(u(x,\theta))|}{\mu},
\qquad \mu=\int_\Omega|G(u(x,\theta))|\,\mathrm dx,
$$

由 Cauchy-Schwarz 它在所有密度中方差最小，但残余方差严格为正：

$$
\sigma_{p^\star}=\Bigl(\int_\Omega|G(u(x,\theta))|\mathrm dx\Bigr)^2
-\Bigl(\int_\Omega G(u(x,\theta))\mathrm dx\Bigr)^2>0 .
$$

这是一个结构性结论：**被积函数变号时，仅靠自适应不能把统计误差消掉，还必须增加样本量** $N_v$。这一点把「自适应采样」的作用界定得很清楚。

### 两个网络交替更新

密度模型取有界 KRnet（见[[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|密度流一页]]）：

$$
p_{\text{bKRnet}}(x,\theta_f)=p_Z\bigl(f_{\text{bKRnet}}(x,\theta_f)\bigr)\,
\bigl|\det\nabla_x f_{\text{bKRnet}}\bigr|,
\qquad Z\sim\mathrm{Unif}([-1,1]^d),
$$

训练准则 $\min_{\theta_f}D_{\mathrm{KL}}(p^\star\|p_{\text{bKRnet}})$ 归结为极小化交叉熵，在重要采样下写成

$$
H(p^\star,p_{\text{bKRnet}})\approx
-\frac1{N_v}\sum_{i=1}^{N_v}
\frac{|G(u(x_i,\theta))|\,\log p_{\text{bKRnet}}(x_i,\theta_f)}{\mu\,\tilde p(x_i)},
\qquad x_i\sim\tilde p .
$$

两个网络在第 $k$ 轮交替更新：

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

### 两个混合模型给出密度下界

比值 $G/p$ 在学到的密度很小处会爆炸，因此论文给密度加一个下界。模型一与均匀分布混合：

$$
p_{\mathrm{mixture}}(x,\theta_f^{k+1})
=\epsilon\,p_{\text{bKRnet}}(x,\theta_f^{k+1})+(1-\epsilon)\,p_{\mathrm{uniform}}(x)
\ \ge\ \frac{1-\epsilon}{|\Omega|} .
$$

模型二与**历史上所有**流递归混合：

$$
p_{\mathrm{mixture}}\bigl(x,\{\theta_f^t\}_{t\le k+1}\bigr)
=\sum_{t=1}^{k+1}\epsilon(1-\epsilon)^{k+1-t}p_{\text{bKRnet}}(x,\theta_f^t)
+\frac{(1-\epsilon)^{k+1}}{|\Omega|},
$$

下界为 $(1-\epsilon)^{k+1}/|\Omega|$。三个算法版本分别对应不混合、混合模型一、混合模型二，四种采样策略（含常规深度 Ritz）在实验中并列比较：二维单峰问题、二维双峰问题、二维奇性问题与高维 Poisson 问题；误差用张量网格（高维用均匀样本）上的相对离散 $L_2$ 误差，激活取 $\sin^3(x)$ 的类 ResNet 结构。

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

| 内容                            | 论文 | 覆盖状态                                      |
| ------------------------------- | ---- | --------------------------------------------- |
| 分数阶 Laplacian 的内外分裂估计 | 66   | 两个恒等式、两个 Beta 采样、$r_\epsilon$ 保护 |
| Caputo 导数的期望表示           | 66   | 公式与采样分布（含前因子核对说明）            |
| 平方偏差与两组独立随机数        | 66   | 无偏性构造与其成立条件                        |
| 极限状态函数与失效概率          | 70   | 定义、安全/失效集、与残差细化的差别           |
| 自适应重要采样                  | 70   | 迭代更新、终止条件、加权最终提议、估计式      |
| Theorem 4.4 与其证明骨架        | 70   | 三条假设、误差界、分集估计、先验因子核对      |
| 定量重采样与子集模拟            | 73   | 作用与位置（限定可核实范围）                  |
| 零方差最优提议与截断高斯混合    | 76   | 最优密度、KL 准则、EM 更新、投影截断          |
| 变分损失的统计误差              | 80   | 误差分解、算例结论、$O(1/\sigma)$ 样本需求    |
| 变号被积函数的残余方差          | 80   | 最优密度与严格正下界的结构含义                |
| 两网络交替更新与两个混合模型    | 80   | 两个目标函数、两个密度下界                    |

## 本页原文

- L. Guo, H. Wu, X. Yu, and T. Zhou, _Monte Carlo fPINNs: deep learning method for forward and inverse problems involving high dimensional fractional partial differential equations_, Comput. Methods Appl. Mech. Engrg. 400 (2022), 115523（预印本 [arXiv:2203.08501](https://arxiv.org/abs/2203.08501)）。
- Z. Gao, L. Yan, and T. Zhou, [_Failure-informed adaptive sampling for PINNs_](https://doi.org/10.1137/22M1527763), SIAM J. Sci. Comput. 45(4) (2023), pp. A1971-A1994（预印本 [arXiv:2210.00279](https://arxiv.org/abs/2210.00279)）。
- Z. Gao, T. Tang, L. Yan, and T. Zhou, [_Failure-informed adaptive sampling for PINNs, Part II: combining with re-sampling and subset simulation_](https://doi.org/10.1007/s42967-023-00312-7), Commun. Appl. Math. Comput. 6 (2024), pp. 1720-1741（预印本 [arXiv:2302.01529](https://arxiv.org/abs/2302.01529)）。
- W. Liu, L. Yan, T. Zhou, and Y. Zhou, [_Failure-informed adaptive sampling for PINNs, Part III: applications to inverse problems_](https://doi.org/10.4208/csiam-am.SO-2023-0059), CSIAM Trans. Appl. Math. 5(3) (2024), pp. 636-670。
- X. Wan, T. Zhou, and Y. Zhou, _Adaptive importance sampling for deep Ritz_, Commun. Appl. Math. Comput. 7(3) (2025), pp. 929-953（预印本 [arXiv:2310.17185](https://arxiv.org/abs/2310.17185)）。
