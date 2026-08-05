---
title: 算子学习与不确定性
description: 编号 75、95、98、107：观测稀疏、分布外与初值可变时，预测应当给出分布
lang: zh
translation: en/computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning
tags:
  - 论文笔记
  - 科学机器学习
  - 算子学习
---

> [!note] 本页覆盖
> 编号 **75**（_J. Comput. Phys._ 510, 2024）、**95**（_J. Comput. Phys._ 2026，[arXiv:2509.25646](https://arxiv.org/abs/2509.25646)）、**98**（投稿 _J. Comput. Phys._，[arXiv:2507.22493](https://arxiv.org/abs/2507.22493)）、**107**（预印本 [arXiv:2606.09434](https://arxiv.org/abs/2606.09434)）。

![算子学习中的三处不确定性](assets/diagrams/tao-zhou-papers/zh/operator-learning-uq.svg)

## 75：用信息瓶颈处理分布外行为

### 问题定位

确定性网络回归与确定性算子网络只给点预测，没有校准的不确定性。标准补救（带 Hamiltonian Monte Carlo 后验采样的贝叶斯神经网络、高斯过程、深度集成）要么随数据量扩展性差，要么在远离训练数据处**不会**放大不确定性。论文的具体目标正是**外推**行为：一个不确定性量化方法应当在输入远离训练分布时报告高不确定性。作者特别指出深度集成在分布外测试输入上并不增加其离散度。

### 三处构造

**信息瓶颈目标带数据增强。** 基本信息瓶颈把输入 $X$ 编码为隐变量 $Z$、再从 $Z$ 解码 $Y$，并极大化 $I(Z;Y)-\beta I(Z;X)$。本文改成

$$
\mathcal L_{\mathrm{IB}}=I(Z;Y)-\beta\,I(\tilde Z;\tilde X),
$$

其中 $\tilde X\sim\tilde p(x)\propto p(x)^{1/\tau}$（$\tau>1$）是训练输入分布的**压平**版本，由一个保体积的正规化流生成，$\tilde Z$ 是其编码。**把 $I(Z;X)$ 换成 $I(\tilde Z;\tilde X)$ 是针对分布外问题的关键改动**：它迫使隐变量对分布外输入近似独立。

**置信感知编码器。**

$$
z=\operatorname{diag}\bigl(m(x)\bigr)\,\bar z(x)
+\operatorname{diag}\bigl(\mathbf 1-m(x)\bigr)\,z_0,
\qquad z_0\sim\mathcal N(0,I),
$$

其中 $\bar z(x)$ 是确定性特征向量，$m(x)\in[0,1]^{\dim z}$ 的每个分量是一个逐特征的置信度。$m(x)=\mathbf 1$ 意味着「$x$ 像训练数据，确定性预测」；$m(x)=\mathbf 0$ 意味着「$x$ 在分布外，$z$ 退化为纯噪声」。

**实际优化的变分界。** 解码器取高斯 $y\mid z\sim\mathcal N(\mu_D(z),\Sigma_D(z))$（$\Sigma_D$ 对角），则

$$
\mathcal L_{\mathrm{VIB}}
=\mathbb E\bigl[\log q_D(y\mid z)\bigr]
-\beta\,\mathbb E\Bigl[\log\frac{q_E(\tilde z\mid\tilde x)}{e(\tilde z)}\Bigr]
\ \le\ \mathcal L_{\mathrm{IB}}-\mathcal H(Y),
$$

其中 $e(\tilde z)$ 是 $\tilde Z$ 边缘的**正规化流**模型，$\mathcal H(Y)$ 是与模型无关的 $Y$ 的熵。

置信门的实现有一处值得记录的细节：$\log m(x)$ 由一个带 LogSigmoid 输出的网络给出，再重新门控使其在生成模型密度低的区域衰减到接近零，具体做法用到 $\tilde p$ 的第 1 与第 5 百分位作为阈值。**这一步是让 $m(x)$ 在极低密度输入处降到 $10^{-6}$ 量级的机制**，也就是分布外行为的实际来源。

### 算子学习版本

DeepONet 基线为 $\mathcal F_\theta(u)(y)=\sum_{i=1}^{n}b_i(o)t_i(y)$，$o=\{u(x_j)\}_{j=1}^m$ 是传感器读数。信息瓶颈被放在**分支网络**里：

$$
\mathcal L_{\mathrm{IBONet}}=I\bigl((Y,Z);s(Y)\bigr)-\beta\,I(\tilde Z;\tilde O),
$$

编码器同上；解码器让 DeepONet **同时**输出均值与对数标准差：

$$
\begin{pmatrix}\mu_D(y,z)\\ \log\sigma_D(y,z)\end{pmatrix}
=\sum_{i=1}^{n} b_i(z)\,t_i(y),
\qquad
\Sigma_D(y,z)=\operatorname{diag}\bigl(\sigma_D(y,z)\bigr)^2 .
$$

实验包括一维不连续函数回归（对比高斯过程、Hamiltonian Monte Carlo 贝叶斯网络、深度集成）、扩散反应方程与对流方程的算子学习、加州房价这一真实表格数据，以及一个大规模气候模型。论文报告的对比结论是：与高斯过程和 Hamiltonian Monte Carlo 一样（而与深度集成不同），本方法在数据稀疏区域给出更大的标准差。

## 95：传感器数量与位置可以变化

### 问题

DeepONet 要求编码输入函数的传感器位置 $\{x_i\}_{i=1}^m$ 在**每个训练与测试样本上数量与位置都相同**，这对实验数据或不规则模拟数据不现实。可变输入算子网络用轻量 transformer 让分支网络置换不变，修好了第一重限制，但仍假设观测足够密集且信息充分，从而存在一个**确定性**的输入到输出映射。当传感器稀疏、不规则或带噪，或者算子本身是随机的（例如一个随机微分方程）时，输出是真正不确定的，任何确定性代理都不可能正确。论文指出的空白是：很少有研究在置换不变的算子学习框架下明确处理不确定性量化。

### 随机算子模型

$$
\hat{\mathcal G}\bigl(h(\mathcal O),z\bigr)(y)=\sum_{n=1}^{p}b^{n}\bigl(h(\mathcal O),z\bigr)\,t^{n}(y),
$$

其中 $h(\mathcal O)$ 是一个**集合 transformer**，给出可变规模观测集 $\mathcal O$ 的置换不变、固定维数表示；$z$ 是带简单先验（标准高斯）的隐随机变量。**隐变量 $z$ 同时承担两种不确定性**：观测不完整带来的认知不确定性，与算子本身的随机性带来的偶然不确定性。

集合 transformer 分两步。输入嵌入把坐标与值分别嵌入再**相加**（而不是拼接）：

$$
\Lambda_{i}=\Lambda_{x}(x_{i})+\Lambda_{\kappa}\bigl(\kappa(x_{i})\bigr)\in\mathbb R^{d_{\mathrm{emb}}} .
$$

注意力池化用 $H$ 个头，每个头有可训练的注意力权网络 $w_l$ 与值投影网络 $v_l$：

$$
h^{(l)}(\mathcal O)=\sum_{i=1}^{m}
\frac{\exp\bigl(w_{l}(\Lambda_{i})/\sqrt{d_{\mathrm{emb}}}\bigr)}
{\sum_{k=1}^{m}\exp\bigl(w_{l}(\Lambda_{k})/\sqrt{d_{\mathrm{emb}}}\bigr)}\,
v_{l}(\Lambda_{i}),
\qquad
h(\mathcal O)=\bigl[h^{(1)},\dots,h^{(H)}\bigr]^{\top} .
$$

**对 $i$ 求和并做 softmax 归一化，正是同时给出顺序不变性与规模无关性的原因。**

### 一处刻意的尺度设计

解码器取网格值 $\bar u=(u(y_1),\dots,u(y_M))^{\top}$ 与

$$
p_{D}(\bar u\mid\mathcal O,z)=\prod_{i=1}^{M}\mathcal N
\Bigl(u(y_{i})\Bigm| \hat{\mathcal G}(h(\mathcal O),z)(y_{i}),\ M\sigma_{u}^{2}\Bigr),
$$

注意人为噪声方差被乘以 $M$。这是**刻意**的：$M\to\infty$ 时离散损失才收敛到定义在**函数空间**上的条件变分自编码器损失而不退化。论文附录把这一点写清楚：对均匀网格与大 $M$，$\sum_i(\cdot)^2/M\approx\|\cdot\|_{L^2}^2/|\Omega_u|$，取 $\sigma=\sqrt{|\Omega_u|}\,\sigma_u$ 就使离散损失成为函数型损失的逼近。

训练损失是负证据下界，闭式为

$$
\mathcal L=\tfrac12\mathbb E\Bigl[-\log|\Sigma_{z}|+\mathrm{tr}\,\Sigma_{z}+\|\mu_{z}\|^{2}-d_{z}\Bigr]
+\frac{1}{2\sigma_{u}^{2}}\mathbb E\Bigl[\sum_{i=1}^{M}
\frac{\bigl(\hat{\mathcal G}(h(\mathcal O),z)(y_{i})-u(y_{i})\bigr)^{2}}{M}\Bigr],
$$

第一项是高斯到标准高斯的解析 KL，第二项是重构项。齐次 Dirichlet 条件由架构硬约束：一维 $[-1,1]$ 上取 $u(x)=(1-x)(1+x)u_\theta(x)$，二维 $[0,1]^2$ 上取 $u(x,y)=x(1-x)y(1-y)u_\theta(x,y)$。

### 让不确定性量化的结论可被检验

评价用 Wasserstein-2 距离与均值、标准差的相对误差。**参考分布不是对训练好的模型做 Monte Carlo，而是精确的高斯过程条件分布**：由于 $\log\kappa\sim\mathcal{GP}(\mu,\mathcal K)$，给定传感器值后的后验仍是高斯过程，

$$
\mu_{\text{post}}(x)=\mu(x)+\mathcal K(x,x_{\text{obs}})
\mathcal K(x_{\text{obs}},x_{\text{obs}})^{-1}(y_{\text{obs}}-\mu(x_{\text{obs}})),
$$

$$
\mathcal K_{\text{post}}(x,x')=\mathcal K(x,x')-\mathcal K(x,x_{\text{obs}})
\mathcal K(x_{\text{obs}},x_{\text{obs}})^{-1}\mathcal K(x_{\text{obs}},x') ,
$$

带噪时把逆换成 $[\mathcal K+\sigma^2I]^{-1}$。每个观测集抽 1000 个后验样本并推过有限差分求解器，构成参考条件分布。**这一步是让不确定性量化的结论可被否证而不是自我循环的原因**，值得作为方法论借鉴。

推断阶段的 $z$ 从**先验**而非编码器采样，因为测试时 $\bar u$ 未知。

## 98：把噪声项换成高斯过程先验

编号 98 沿用编号 75 的置信感知编码器思路，但做两处替换：把噪声项 $z_0\sim\mathcal N(0,I)$ 换成**高斯过程先验**，把高斯解码器的均值换成**神经算子**。前一处替换的意义是给出带**相关结构**的不确定性：独立高斯噪声给出的是逐点不确定性，而高斯过程先验让不同位置的不确定性相关，这更符合 PDE 解的实际。

> [!note] 覆盖进度
> 编号 98 的正文已在预印本中公开，但本页尚未逐式核对其损失函数与训练算法，因此只给出与编号 75 的差别与其含义。

## 107：把学习对象从初值依赖换成转移密度

### 问题

在集合预报与数据同化中，随机动力学是固定的，而**初始分布随每个场景或后验更新而变**，因此经典 Fokker-Planck 求解器或物理信息求解器每次都要从头重跑。把它写成算子学习（初始分布 $\mapsto$ 时间依赖解）是自然的，但论文指出直接把 DeepONet、Fourier 神经算子或 Green 函数型算子学习用于 Fokker-Planck 方程「需要对边界条件与概率密度约束（非负性与单位质量）做额外处理」。

### 核心一步：学转移密度

论文的替代方案是学**转移密度**而不是初值依赖：

$$
\frac{\partial p(\bm x,t\mid\bm x_{0})}{\partial t}
=\mathcal L^{*}_{\bm f,\bm g}\,p(\bm x,t\mid\bm x_{0}),
\qquad
p(\bm x,0\mid\bm x_{0})=\delta(\bm x-\bm x_{0}),
$$

然后由 Chapman-Kolmogorov 方程恢复任意初值：

$$
p(\bm x,t)=\int_{\Omega_{0}}p(\bm x,t\mid\bm x_{0})\,p_{0}(\bm x_{0})\,\mathrm d\bm x_{0} .
$$

**这是全篇的核心动作：被学习的对象与 $p_0$ 无关，因此新初值完全不需要重训。** 代价是转移密度满足一个带 **Dirac 初值**的方程，而 Dirac 是奇异的，任何可逆映射都无法从光滑基分布到达它。

### 条件流与线性化基过程

密度模型取条件正规化流：令变换 $T$ 与基分布同时依赖 $(\bm x_0,t)$，

$$
p_{\bm X_{t}\mid\bm X_{0}=\bm x_{0}}(\bm x)
=p_{\bm Z_{t}(\bm x_{0})}\bigl(T(\bm x;\theta(\bm x_{0},t))\bigr)\,
\bigl|\det\nabla_{\bm x}T(\bm x;\theta(\bm x_{0},t))\bigr| .
$$

由于是流，非负性与单位质量**由构造成立**——这正是论文对「概率密度约束」这一反对意见的回答，与[[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|密度流一页]]的取向一致。

处理 Dirac 奇性的办法是给基过程一个**线性化的随机微分方程**：对漂移做一阶 Taylor、对扩散做零阶 Taylor，都在 $\bm x_0$ 处展开，

$$
\mathrm d\widehat{\bm X}_{t}
=\bigl(\bm f(\bm x_{0})+\nabla\bm f(\bm x_{0})(\widehat{\bm X}_{t}-\bm x_{0})\bigr)\mathrm dt
+\bm g(\bm x_{0})\mathrm d\bm W_{t},
\qquad \widehat{\bm X}_{0}=\bm x_{0},
$$

这是一个高斯过程，其均值与协方差有闭式：

$$
\mathbb E\bigl[\widehat{\bm X}_{t}\bigr]
=\bm x_{0}+\int_{0}^{t}e^{\nabla\bm f(\bm x_{0})(t-s)}\bm f(\bm x_{0})\,\mathrm ds,
$$

$$
\bm\Sigma\bigl[\widehat{\bm X}_{t}\bigr]
=\int_{0}^{t}e^{\nabla\bm f(\bm x_{0})(t-s)}\bm g(\bm x_{0})\bm g(\bm x_{0})^{\top}
e^{\nabla\bm f(\bm x_{0})^{\top}(t-s)}\,\mathrm ds .
$$

**这个基过程在 $t\to0$ 时自动退化为以 $\bm x_0$ 为中心、协方差趋零的高斯**，因此流不需要从光滑分布跨越到 Dirac：奇性被基分布本身吸收了。这是把架构与解的解析结构对齐的一个清晰例子。

> [!note] 题名与投稿状态
> 该文实际题名为 _A transition-density-based operator learning method for Fokker-Planck equations with various initial conditions_，与主页所列略有差别；主页题名接近其页眉，可能是早期工作题名。预印本采用 SIAM 格式，但未标注期刊、DOI 或期刊引用信息。

## 四篇的对照

| 编号 | 不确定性来自       | 表示方式                  | 约束由谁保证         |
| ---- | ------------------ | ------------------------- | -------------------- |
| 75   | 输入落在分布外     | 置信门 + 高斯解码器       | 无（软约束）         |
| 95   | 观测稀疏或算子随机 | 隐变量 + 集合 transformer | 边界条件由架构硬约束 |
| 98   | 同 75，带相关结构  | 高斯过程先验 + 神经算子   | 待补                 |
| 107  | 不涉及（初值可变） | 条件流 + 线性化基过程     | 密度约束由流保证     |

一条贯穿的判断：**「预测应当给出分布」这件事，在算子学习里可以落在三个不同的位置。** 编号 75 把它放在编码器的置信门上，编号 95 放在隐变量上，编号 98 放在先验的相关结构上。而编号 107 说明另一件事：有时正确的做法不是给预测加不确定性，而是换一个与参数无关的学习对象。

## 本页原文

- L. Guo, H. Wu, W. Zhou, Y. Wang, and T. Zhou, [_IB-UQ: information bottleneck based uncertainty quantification for neural function regression and neural operator learning_](https://doi.org/10.1016/j.jcp.2024.113089), J. Comput. Phys. 510 (2024), 113089（预印本 [arXiv:2302.03271](https://arxiv.org/abs/2302.03271)）。
- L. Ma, L. Guo, H. Wu, and T. Zhou, [_Deep set based operator learning with uncertainty quantification_](https://doi.org/10.1016/j.jcp.2026.115011), J. Comput. Phys. (2026)（预印本 [arXiv:2509.25646](https://arxiv.org/abs/2509.25646)）。
- X. Feng, L. Guo, X. Wan, H. Wu, T. Zhou, and W. Zhou, _LVM-GP: uncertainty-aware PDE solver via coupling latent variable model and Gaussian process_, [arXiv:2507.22493](https://arxiv.org/abs/2507.22493)，投稿 J. Comput. Phys.
- L. Zeng, X. Wan, Y. Wang, F. Nobile, and T. Zhou, _A transition-density-based operator learning method for Fokker-Planck equations with various initial conditions_, [arXiv:2606.09434](https://arxiv.org/abs/2606.09434)。
