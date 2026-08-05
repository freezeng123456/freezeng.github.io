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

### 直觉

确定性网络在训练数据附近与远离训练数据处给出的输出，从形式上看没有任何区别：都是一串数字，都不带「我不知道」的标记。深度集成之所以不能补上这一点，是因为集成中的每个成员看的是同一批数据，在外推区域它们往往**一致地**给出同一个错误答案，于是集成方差反而小。

这篇论文的办法是在网络中间插一道**闸门**。输入先被编码成隐表示 $z$，再由 $z$ 解码出 $y$；闸门 $m(x)$ 决定这个 $z$ 有多少来自输入、有多少来自纯噪声。$x$ 落在训练数据密集的地方，闸门全开，$z$ 完全由 $x$ 决定，预测是确定的；$x$ 落在训练数据看不到的地方，闸门关闭，$z$ 退化成标准高斯噪声，解码器只能吐出一个宽分布。**因此「不知道」不是事后估出来的，而是架构在低密度区强制发生的行为。**

剩下的问题是怎么训练这个闸门：训练集里没有分布外样本，闸门没有理由学会关闭。论文用信息瓶颈目标的一个改写来提供这个理由——把压缩项从「压缩训练输入」改成「压缩一个**压平后的**输入分布」，压平分布的支撑比训练分布宽，于是闸门在训练时就被迫在训练分布的外围练习关闭。

### 问题设定

基本信息瓶颈把输入 $X$ 编码为隐变量 $Z$、再从 $Z$ 解码 $Y$，并极大化

$$
I(Z;Y)-\beta\,I(Z;X),
$$

$I(\cdot;\cdot)$ 为互信息，$\beta\in[0,1]$ 在预测与压缩之间取舍。论文的具体目标是**外推**行为：一个不确定性量化方法应当在输入远离训练分布时报告高不确定性。标准补救（带 Hamiltonian Monte Carlo 后验采样的贝叶斯神经网络、高斯过程、深度集成）要么随数据量扩展性差，要么在远离训练数据处**不会**放大不确定性；作者特别指出深度集成在分布外测试输入上并不增加其离散度。

### 推导

**第一步：把压缩项换成压平分布上的压缩项。** 论文的目标函数（其式 6）为

$$
\mathcal L_{\mathrm{IB}}=I(Z;Y)-\beta\,I(\tilde Z;\tilde X),
$$

其中 $\tilde X\sim\tilde p(x)\propto p(x)^{1/\tau}$（$\tau>1$）是训练输入分布的**压平**版本，由一个保体积的正规化流（general incompressible-flow network，GIN）生成，$\tilde Z$ 是其编码。**把 $I(Z;X)$ 换成 $I(\tilde Z;\tilde X)$ 是针对分布外问题的关键改动**：压平分布把质量推向训练分布的外围，压缩项在那里惩罚 $\tilde Z$ 对 $\tilde X$ 的依赖，于是隐变量对分布外输入近似独立。

**第二步：置信感知编码器。** 编码器（其式 7）写成一个逐分量的插值

$$
z=\operatorname{diag}\bigl(m(x)\bigr)\,\bar z(x)
+\operatorname{diag}\bigl(\mathbf 1-m(x)\bigr)\,z_0,
\qquad z_0\sim\mathcal N(0,I),
$$

其中 $\bar z(x)$ 是确定性特征向量，$m(x)\in[0,1]^{\dim z}$ 的每个分量是一个逐特征的置信度。$m(x)=\mathbf 1$ 意味着「$x$ 像训练数据，确定性预测」；$m(x)=\mathbf 0$ 意味着「$x$ 在分布外，$z$ 退化为纯噪声」。

**第三步：高斯解码器与变分界。** 解码器取 $y\mid z\sim\mathcal N(\mu_D(z),\Sigma_D(z))$，$\Sigma_D$ 对角，$\mu_D$ 与 $\log\operatorname{diag}\Sigma_D$ 都是多层感知机的输出。实际优化的是（其式 8–9）

$$
\mathcal L_{\mathrm{VIB}}
=\mathbb E_{(x,y,z)\sim p(x,y)q_E(z\mid x)}\bigl[\log q_D(y\mid z)\bigr]
-\beta\,\mathbb E_{(\tilde x,\tilde z)\sim\tilde p(\tilde x)q_E(\tilde z\mid\tilde x)}
\Bigl[\log\frac{q_E(\tilde z\mid\tilde x)}{e(\tilde z)}\Bigr]
\ \le\ \mathcal L_{\mathrm{IB}}-\mathcal H(Y),
$$

其中 $e(\tilde z)$ 是 $\tilde Z$ 边缘的**正规化流**模型，$\mathcal H(Y)$ 是与模型无关的 $Y$ 的熵——正因为它与参数无关，极大化 $\mathcal L_{\mathrm{VIB}}$ 就是在提升 $\mathcal L_{\mathrm{IB}}$ 的一个下界。$(q_E,q_D,e)$ 都训练到位且 $\mathcal L_{\mathrm{VIB}}$ 取极大时取等。

小批量估计（补充材料 S4）为

$$
\widehat{I(Z;Y)}\approx\frac1B\sum_{b=1}^{B}\log q_D(y_b\mid z_b),
\qquad
\widehat{I(\tilde Z;\tilde X)}\approx\frac1B\sum_{b=1}^{B}
\bigl[\log q_E(\tilde z_b\mid\tilde x_b)-\log e(\tilde z_b)\bigr].
$$

**第四步：闸门真正关上的机制。** 置信门的实现有一处值得记录的细节：$\log m(x)$ 由一个带 LogSigmoid 输出的网络给出，再重新门控使其在生成模型密度低的区域衰减到接近零，

$$
\log m(x)^{\mathrm{new}}=\mathrm{LogSigmoid}\!\left(
\log\frac{m(x)^-}{1-m(x)^-}
+\log\!\left(\tanh\!\left(
\frac{\mathrm{relu}\bigl(\log\tilde p(x)-\log\tilde p'\bigr)}
{\log\tilde p''-\log\tilde p'}\right)+10^{-12}\right)\right),
$$

其中 $m(x)^-=\min\{m(x),\mathrm{Sigmoid}(10^6)\}$，$\tilde p'$、$\tilde p''$ 分别是 GIN 密度 $p_G(x)$ 的第 1 与第 5 百分位。$\log\tilde p(x)$ 低于第 1 百分位时 relu 归零、$\tanh$ 项趋于 $10^{-12}$，$\log m$ 被推向 $-\infty$。**这一步是让 $m(x)$ 在极低密度输入处降到 $10^{-6}$ 量级的机制**，也就是分布外行为的实际来源。

### 算子学习版本

DeepONet 基线为 $\mathcal F_\theta(u)(y)=\sum_{i=1}^{n}b_i(o)t_i(y)$，$o=\{u(x_j)\}_{j=1}^m$ 是传感器读数，训练损失为

$$
\mathcal L(\theta)=\frac{1}{NP}\sum_{l}\sum_{j}
\bigl|\mathcal F_\theta(u^l)(y^l_j)-s^l(y^l_j)\bigr|^2 .
$$

信息瓶颈被放在**分支网络**里：

$$
\mathcal L_{\mathrm{IBONet}}=I\bigl((Y,Z);s(Y)\bigr)-\beta\,I(\tilde Z;\tilde O),
$$

编码器同上，只是作用在传感器读数上：$z(o)=\operatorname{diag}(m(o))\bar z(o)+\operatorname{diag}(\mathbf 1-m(o))z_0$。解码器让 DeepONet **同时**输出均值与对数标准差：

$$
\begin{pmatrix}\mu_D(y,z)\\ \log\sigma_D(y,z)\end{pmatrix}
=\sum_{i=1}^{n} b_i(z)\,t_i(y),
\qquad
\Sigma_D(y,z)=\operatorname{diag}\bigl(\sigma_D(y,z)\bigr)^2 ,
$$

相应的界为 $\mathcal L_{\mathrm{VIBONet}}\le\mathcal L_{\mathrm{IBONet}}-\mathcal H(s(Y))$。

> [!warning] 归一化常数的一处印刷错误
> 上式的归一化常数印为 $NP$，而内层求和的上界写作 $M$。按维数只能有一个是对的，这是原文的一处笔误。

训练算法（补充材料 S4，算法 1）依次为：先把 GIN 拟合到 $\{x_1,\dots,x_N\}$ 得到 $\tilde p$；每步抽一个小批量，可选地做 mixup（抽随机置换 $(I_1,\dots,I_B)$ 与 $\lambda_b\sim\mathrm{Beta}(\alpha,\alpha)$，置 $x_b:=\lambda_bx_b+(1-\lambda_b)x_{I_b}$，$y_b$ 同理）；从 GIN 抽 $\tilde x_1,\dots,\tilde x_B$；按置信门构造 $z_b$ 与 $\tilde z_b$；估计 $\widehat{\mathcal L_{\mathrm{VIB}}}$ 并对 $q_D,q_E,e$ 的全部参数做梯度**上升**。预测时对新的 $x$ 反复抽 $(Y,Z)\sim q_E(z\mid x)q_D(y\mid z)$，报告样本均值与方差，即用 Monte Carlo 求 $\hat p(y\mid x)=\int q_E(z\mid x)q_D(y\mid z)\mathrm dz$。

### 定理

**没有收敛定理。** 全篇的理论内容是两处：补充材料 S1 中变分下界的推导，以及 S5 中「GIN 数据增强为什么使隐变量对分布外输入近似独立」的分析。论文没有给出误差率或一致性结论。

### 数值实验

五组实验，覆盖回归、算子学习、真实表格数据与大规模气候模型：

| 例子           | 设定                                                                                   | 对比基线                           |
| -------------- | -------------------------------------------------------------------------------------- | ---------------------------------- |
| 一维不连续回归 | $u(x)=\tfrac12[\sin^3(2\pi x)-1]$ 于 $[-1,0)$，$\tfrac12[\sin^3(3\pi x)+1]$ 于 $[0,1]$ | 高斯过程、HMC 贝叶斯网络、深度集成 |
| 扩散反应方程   | 算子学习                                                                               | —                                  |
| 对流方程       | 算子学习                                                                               | —                                  |
| 加州房价       | 真实表格数据回归                                                                       | —                                  |
| 大规模气候模型 | 大数据量下的效率论证                                                                   | HMC 贝叶斯网络                     |

不连续回归的设定完整可查：$N=32$ 个等距含噪样本，但**只落在 $[-0.8,-0.2]\cup[0.2,0.8]$**，噪声 $\epsilon_u\sim\mathcal N(0,\sigma_u^2)$，$\beta=0.3$。留白区间 $[-1,-0.8]\cup[-0.2,0.2]\cup[0.8,1]$ 正是用来检验外推行为的。论文报告的对比结论是：与高斯过程和 Hamiltonian Monte Carlo 一样（而与深度集成不同），本方法在这三段数据稀疏区间给出更大的标准差。

这组实验建立的是：**置信门确实在无数据区打开了不确定性**，而且这一行为在低维回归、算子学习与真实数据上都出现。它没有建立的是校准质量——论文没有给出覆盖率、CRPS 之类的标定指标，「更大的标准差」只是与参考区间的定性比较。

> [!warning] 定量结果未核对
> 具体误差数值与运行时间位于图与补充材料 S7–S8 中，本次核对未逐项转录，因此本页不列这些数字。上表只记录已核对的实验清单与设定。

### 与其他论文的关系

这是该组三篇「用隐变量做不确定性量化」论文中的第一篇。编号 98 沿用**同一个置信感知编码器**，但把噪声项 $z_0\sim\mathcal N(0,I)$ 换成**高斯过程先验**、把高斯解码器的均值换成神经算子；编号 95 保留条件变分自编码器与证据下界这套机器，但把置信门换成置换不变的集合 transformer 编码器。边缘模型 $e(\tilde z)$ 用的正规化流把本页接到[[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|密度流一页]]，GIN 本身就是一个保体积的正规化流。

## 95：传感器数量与位置可以变化

### 直觉

DeepONet 的分支网络吃的是一个**固定长度的向量**——第一个分量是第一个传感器的读数，第二个是第二个传感器的读数，如此等等。这意味着传感器一旦少一个、挪一个位置，网络就完全用不了。实验数据和不规则模拟数据几乎总是这种情形。

修法是把分支网络的输入从「向量」改成「**集合**」：每个传感器贡献一个 $(x_i,\kappa(x_i))$ 对，网络先把每个对嵌入到同一个特征空间，再用 softmax 加权求和把它们汇总成一个定长向量。求和是对集合元素做的，所以顺序无关；softmax 的分母跑遍同一个集合，所以数量也无关。

但只做到这里还不够。传感器稀疏时，同一组读数其实对应许多条不同的真实系数曲线，因而对应许多个不同的解——**输入到输出根本不是一个函数，而是一个条件分布**。任何确定性代理在原理上就是错的，它最多学到条件均值。所以论文再加一个隐变量 $z$：给定观测集，不同的 $z$ 给出不同的候选解，$z$ 的分布就是预测分布。这个 $z$ 同时兼管两件事——观测不完整造成的**认知**不确定性，以及算子本身随机（例如系数来自一个随机场，或右端项是第二个未观测的随机源）造成的**偶然**不确定性。

### 问题设定

DeepONet 要求编码输入函数的传感器位置 $\{x_i\}_{i=1}^m$ 在**每个训练与测试样本上数量与位置都相同**。可变输入算子网络（VIDON）用轻量 transformer 让分支网络置换不变，修好了第一重限制，但仍假设观测足够密集且信息充分，从而存在一个**确定性**的输入到输出映射。论文指出的空白是：很少有研究在置换不变的算子学习框架下明确处理不确定性量化。论文还记下一条前人工作的次要限制：原始 DeepONet 处理随机算子时用的是 Karhunen–Loève 展开，因而被限制在**低维**随机输入上。

基线 DeepONet（其式 2.1）为

$$
\hat{\mathcal G}(\kappa)(y)=\sum_{n=1}^{p}
\underbrace{b^{n}(\mathcal O)}_{\text{分支}}\,
\underbrace{t^{n}(y)}_{\text{主干}}
=\sum_{n=1}^{p}b^{n}\bigl(\kappa(x_1),\dots,\kappa(x_m)\bigr)\,t^{n}(y),
$$

$\mathcal O=\{(x_i,\kappa(x_i))\}_{i=1}^m$ 是观测集，$p$ 是截断秩。本文的中心对象（其式 2.2）为

$$
\hat{\mathcal G}\bigl(h(\mathcal O),z\bigr)(y)=\sum_{n=1}^{p}b^{n}\bigl(h(\mathcal O),z\bigr)\,t^{n}(y),
$$

其中 $h(\mathcal O)$ 是一个**集合 transformer**，给出可变规模观测集 $\mathcal O$ 的置换不变、固定维数表示；$z$ 是带简单先验（标准高斯）的隐随机变量。

### 推导

**集合 transformer。** 分两步。输入嵌入（其式 3.1）把坐标与值分别嵌入再**相加**（而不是拼接）：

$$
\Lambda_{i}=\Lambda_{x}(x_{i})+\Lambda_{\kappa}\bigl(\kappa(x_{i})\bigr)\in\mathbb R^{d_{\mathrm{emb}}},
\qquad i=1,\dots,m .
$$

注意力池化（其式 3.2–3.3）用 $H$ 个头，每个头有可训练的注意力权网络 $w_l$ 与值投影网络 $v_l$：

$$
h^{(l)}(\mathcal O)=\sum_{i=1}^{m}
\underbrace{\frac{\exp\bigl(w_{l}(\Lambda_{i})/\sqrt{d_{\mathrm{emb}}}\bigr)}
{\sum_{k=1}^{m}\exp\bigl(w_{l}(\Lambda_{k})/\sqrt{d_{\mathrm{emb}}}\bigr)}}_{\text{注意力权}}\,
\underbrace{v_{l}(\Lambda_{i})}_{\text{值投影}},
\qquad
h(\mathcal O)=\bigl[h^{(1)},\dots,h^{(H)}\bigr]^{\top} .
$$

**对 $i$ 求和并做 softmax 归一化，正是同时给出顺序不变性与规模无关性的原因**：分子与分母跑遍同一个指标集，置换 $\{1,\dots,m\}$ 不改变结果，而 $m$ 变化时结构仍然良定义。该模块声明沿用 VIDON 的架构。

**解码器与一处刻意的尺度设计。** $u$ 用网格值 $\bar u=(u(y_1),\dots,u(y_M))^{\top}$ 表示，解码器（其式 3.4）为

$$
p_{D}(\bar u\mid\mathcal O,z)=\prod_{i=1}^{M}\mathcal N
\Bigl(u(y_{i})\Bigm| \hat{\mathcal G}(h(\mathcal O),z)(y_{i}),\ M\sigma_{u}^{2}\Bigr),
$$

先验 $p(z)=\mathcal N(z\mid0,I)$，$\sigma_u^2$ 是一个小的正超参数（「人为噪声」）。注意人为噪声方差被乘以 $M$。这是**刻意**的（其注记 3.1）：$M\to\infty$ 时离散损失才收敛到定义在**函数空间**上的条件变分自编码器损失而不退化。附录 A.3 把这一点写清楚：对均匀网格与大 $M$，$\sum_i(\cdot)^2/M\approx\|\cdot\|_{L^2}^2/|\Omega_u|$，取 $\sigma=\sqrt{|\Omega_u|}\,\sigma_u$ 就使离散损失成为函数型损失 $\mathcal L_f$ 的逼近。

**编码器与证据下界。** 变分后验取

$$
q_{E}(z\mid\mathcal O,\bar u)=\mathcal N
\Bigl(z\Bigm|\mu_{z}\bigl(h(\mathcal O),\bar u\bigr),\ \Sigma_{z}\bigl(h(\mathcal O),\bar u\bigr)\Bigr),
$$

$\Sigma_z$ **对角**，且编码器复用**同一个** $h(\mathcal O)$ 嵌入。对边缘似然用 Jensen 不等式，

$$
\log p_{D}(\bar u\mid\mathcal O)=\log\!\int p_{D}(\bar u\mid\mathcal O,z)p(z)\,\mathrm dz
\ \ge\ -D_{\mathrm{KL}}\bigl(q_{E}(z\mid\bar u,\mathcal O)\,\|\,p(z)\bigr)
+\mathbb E_{q_{E}}\bigl[\log p_{D}(\bar u\mid\mathcal O,z)\bigr],
$$

$q_E$ 等于真后验 $p_D(z\mid\bar u,\mathcal O)$ 时取等。训练损失是负证据下界，因先验与后验都是高斯而有闭式（其式 3.5）：

$$
\mathcal L=\tfrac12\mathbb E\Bigl[-\log|\Sigma_{z}|+\mathrm{tr}\,\Sigma_{z}+\|\mu_{z}\|^{2}-d_{z}\Bigr]
+\frac{1}{2\sigma_{u}^{2}}\mathbb E\Bigl[\sum_{i=1}^{M}
\frac{\bigl(\hat{\mathcal G}(h(\mathcal O),z)(y_{i})-u(y_{i})\bigr)^{2}}{M}\Bigr],
$$

第一项是高斯到标准高斯的解析 KL，第二项是重构项，对均匀网格与大 $M$ 正比于均方 $L^2$ 范数。$\mathcal L$ 上界控制 $\mathbb E[-\log p_D(\bar u\mid\mathcal O)]$。

**边界条件硬约束。** 齐次 Dirichlet 条件由架构承担：一维 $[-1,1]$ 上取 $u(x)=(1-x)(1+x)u_\theta(x)$，二维 $[0,1]^2$ 上取 $u(x,y)=x(1-x)y(1-y)u_\theta(x,y)$，$u_\theta$ 是网络原始输出。

**训练与推断。** 训练时每个小批量随机选一个传感器数 $m\in\mathcal M$，为每个样本随机抽 $m$ 个位置生成 $\mathcal O_i$，算 $h(\mathcal O_i)$，从 $q_E$ 抽 $z_i$，形成随机损失

$$
\hat{\mathcal L}=\frac{1}{2B}\sum_{i=1}^{B}
\Bigl(-\log|\Sigma_{z}|+\mathrm{tr}\,\Sigma_{z}+\|\mu_{z}\|^{2}-d_{z}\Bigr)
+\frac{1}{2BM\sigma_{u}^{2}}\sum_{i=1}^{B}\sum_{j=1}^{M}
\Bigl(\hat{\mathcal G}(h(\mathcal O_{i}),z_{i})(y_{j})-u(y_{j})\Bigr)^{2},
$$

再用 Adam 更新。**推断阶段的 $z$ 从先验而非编码器采样**，因为测试时 $\bar u$ 未知；反复抽 $z\sim\mathcal N(0,I)$ 就累积出 $u$ 的条件分布。

### 定理

**论文没有编号定理或命题。** 理论内容是 3.3 节的证据下界推导、注记 3.1 关于 $\sigma_u^2$ 乘 $M$ 的说明，以及附录 A.3 的函数型变分自编码器极限论证。

置换不变性这一条的理论依据来自集合学习的背景文献而非本文。可援引的是 Deep Sets 的定理 2：**取值于可数全集的集合函数是置换不变的，当且仅当它能分解为 $\rho\bigl(\sum_{x\in X}\phi(x)\bigr)$**。本文的注意力池化正是这种形式的 softmax 加权版本。

> [!warning] 背景文献的核对边界
> Deep Sets 的定理 2 及其作用域限制已核对：原文自己说明「$\mathfrak X$ 不可数（如 $\mathfrak X=\mathbb R$）时，只能对**固定规模**的集合证明该分解」，因此那条干净的充要条件是可数全集下的结论。Set Transformer 论文的 MAB/SAB/ISAB/PMA 具体公式**未核对**——本次只读到其补充材料。因此本页对置换不变性的论证只依赖 Deep Sets 定理与本文自己的式 3.2，不依赖那组未核对的公式。

### 数值实验

全局设定：Tanh 激活，截断秩 $p=100$，隐维数 $d_z=10$（因为 $W_2$ 距离在此处趋稳，其图 2(a)），一维算例用 10000 组训练样本。**所有误差都是相对 $L^2$ 误差并乘以 $10^{-2}$**，即表中「5.35」表示 $5.35\times10^{-2}$。

评价用 Wasserstein-2 距离

$$
W_{2}(\mu,\nu)=\Bigl(\inf_{\gamma\in\Gamma(\mu,\nu)}\int_{M\times M}d^{2}(x_{1},x_{2})\,\mathrm d\gamma\Bigr)^{1/2}
$$

与均值、标准差的相对误差

$$
err_{E[u]}=\frac{\|E[u]-E_{\theta}[u]\|_{2}}{\|E[u]\|_{2}},
\qquad
err_{\sigma[u]}=\frac{\|\sigma[u]-\sigma_{\theta}[u]\|_{2}}{\|\sigma[u]\|_{2}} .
$$

**参考分布不是对训练好的模型做 Monte Carlo，而是精确的高斯过程条件分布**（附录 A.2）：由于 $\log\kappa\sim\mathcal{GP}(\mu,\mathcal K)$，给定传感器值后的后验仍是高斯过程，

$$
\mu_{\text{post}}(x)=\mu(x)+\mathcal K(x,x_{\text{obs}})
\mathcal K(x_{\text{obs}},x_{\text{obs}})^{-1}(y_{\text{obs}}-\mu(x_{\text{obs}})),
$$

$$
\mathcal K_{\text{post}}(x,x')=\mathcal K(x,x')-\mathcal K(x,x_{\text{obs}})
\mathcal K(x_{\text{obs}},x_{\text{obs}})^{-1}\mathcal K(x_{\text{obs}},x') ,
$$

带噪时把逆换成 $[\mathcal K(x_{\text{obs}},x_{\text{obs}})+\sigma^2I]^{-1}$。**每个观测集抽 1000 个后验样本并推过有限差分求解器**，构成参考条件分布。**这一步是让不确定性量化的结论可被否证而不是自我循环的原因**，值得作为方法论借鉴。

**例 1：一维扩散。** $-\tfrac1{10}\tfrac{\mathrm d}{\mathrm dx}\bigl(k(x)\tfrac{\mathrm du}{\mathrm dx}\bigr)=f(x)$ 于 $[-1,1]$，$u(\pm1)=0$，$f(x)=2\sin(2\pi x)$，学 $\mathcal G:k\mapsto u$。输入取 $\log k\sim\mathcal{GP}\bigl(\sin(2\pi x),\sigma^2\exp(-(x-x')^2/l^2)\bigr)$，$\sigma=0.5$、$l=0.1$；401 点传感网格，二阶有限差分求解器，$N=10^4$ 组数据对，训练输出用 101 个点。可变传感器数 $\mathcal M=\{1,\dots,10\}$，10 万次迭代，测试在 401 点网格上，5 次重复。

| $m$ | UQ-SONet $err_{E[u]}$ | VIDON $err_{E[u]}$ |
| --- | --------------------- | ------------------ |
| 1   | 4.47 ± 0.83           | 5.35 ± 0.15        |
| 7   | 5.88 ± 0.74           | 8.73 ± 1.03        |

UQ-SONet 在 $m=1$ 到 $10$ 的**每一个**传感器数上都优于 VIDON；$err_{\sigma[u]}$ 落在 7.68–8.95 之间。预测带随 $m$ 增大而收窄。噪声鲁棒性另测（其表 3）：对 $k$ 施乘性噪声 $k(x)\exp(\epsilon(x))$、$\epsilon\sim\mathcal N(0,\sigma^2)$、$\sigma\in\{0.1,0.3,0.5,0.7,1.0\}$，误差基本持平——$m=7$ 时 $err_{E[u]}$ 从 7.49 变到 6.50，$\sigma$ 从 0.1 变到 1.0。**退化对噪声水平出人意料地不敏感。**

**例 2：二维 Poisson。** $-\tfrac1{10}\Delta u=f$ 于 $[0,1]^2$，齐次 Dirichlet，学 $\mathcal G:f\mapsto u$。$f\sim\mathcal{GP}$，均值 $4(\sin2\pi x+\sin2\pi y)$、$l_1=l_2=0.1$，在 $101\times101$ 网格上生成；$N=80000$ 组数据对；$\mathcal M=\{1,2,3,4\}$；传感器位置用**规则空间聚类**（`pyemma.coordinates.cluster_regspace`）在 $[0.1,0.9]^2$ 内的置乱 $81\times81$ 网格上选取，$m=2$ 时最小间距 $d_{\min}=0.8$，$m=3,4$ 时 $d_{\min}=0.5$；训练输出网格 $51\times51$，测试 $101\times101$；2 万次迭代，3 次重复。

| $m$ | UQ-SONet（无噪） | VIDON（无噪） | UQ-SONet（噪声 $\mathcal N(0,1.0^2)$） |
| --- | ---------------- | ------------- | -------------------------------------- |
| 1   | 4.00             | 3.26          | 3.40                                   |
| 2   | 4.83             | 4.67          | 3.40                                   |
| 3   | 5.92             | 5.99          | 4.16                                   |
| 4   | 5.48             | 7.01          | 5.54                                   |

**这一组是全篇最诚实的一张表：$m=1,2$ 时 VIDON 反而更准**，UQ-SONet 只在 $m=3,4$ 时占优。也就是说「加了不确定性还顺带更准」这条并非无条件成立。

**例 3：一维椭圆随机微分方程。** $-\tfrac1{10}\tfrac{\mathrm d}{\mathrm dx}\bigl(k(x;\omega)\tfrac{\mathrm d}{\mathrm dx}u(x;\omega)\bigr)=f(x;\omega)$；$\log k$ 取 $l=0.05$、$\sigma=0.3$、$\mu=\sin(\pi x+1)$，$f$ 取 $l=0.1$、$\sigma=0.1$、$\mu=\sin(2\pi x)+0.1$；学 $\mathcal G:k\mapsto u$，此处 **$f$ 是第二个、未被观测的随机源**。$m=1,\dots,10$、5 次重复下，$err_{E[u]}$ 落在 1.76–2.45，$err_{\sigma[u]}$ 落在 3.61–5.66——全篇最好的精度。定性上更重要的观察是：**预测离散度在 $m\to10$ 时并不收缩到零**，因为 $f$ 的内在随机性始终存在，被消解的只是认知部分。这是全篇对「两类不确定性被分开」这一主张最干净的证据。

**例 4：二维椭圆随机微分方程。** $-\nabla\cdot(k(x,y;\omega)\nabla u)=f(x,y;\omega)$ 于 $[0,1]^2$；$\log k$ 零均值、$l_1=l_2=0.1$；$f$ 均值 $4(\sin2\pi x+\sin2\pi y)$；学 $\mathcal G:f\mapsto u$；$\mathcal M=\{1,2,3,4\}$，5 万次迭代，训练网格 $51\times51$、测试 $101\times101$。在 $x=0.5,0.7$ 与 $y=0.5,0.7$ 四条切片上评价，3 次重复：$err_{E[u]}$ 落在 3.23–7.57，$err_{\sigma[u]}$ 落在 6.52–8.50。

**例 5：二维不可压 Navier-Stokes。** 涡量-速度形式、周期边界，$\partial_tw+u\cdot\nabla w=\nu\Delta w+f$、$\nabla\cdot u=0$，$\nu=0.001$，$f(x,y)=0.1\sin(2\pi(x+y))+0.1\cos(2\pi(x+y))$；学 $w_0\mapsto w|_{t=10}$。初始涡量取 $g(\bm x;\omega)=x^{1/3}(1-x)^{1/3}y^{1/3}(1-y)^{1/3}h(\bm x;\omega)$，$h$ 为零均值高斯过程、$l_1=l_2=0.1$；$100\times100$ 网格，伪谱流函数求解器；训练观测在 $50\times50$ 上；4 批各 2 万，5 万次迭代，3 次重复。

| $m$ | UQ-SONet：$w$ 在 $t=10$ 的 $err_{E}$ | VIDON |
| --- | ------------------------------------ | ----- |
| 1   | 2.35                                 | 6.01  |
| 2   | 1.89                                 | 5.98  |
| 3   | 2.55                                 | 8.74  |
| 4   | 2.69                                 | 8.08  |

$err_\sigma$ 落在 6.25–7.32。**均值误差约为 VIDON 的三分之一，这是全篇最强的结果**，且正是在最难的算例上。

**深度集成基线（附录 A.4）。** 用不同随机种子与数据洗牌重训 VIDON 构成集成。报告的结论是：集成的均值是准的，但其不确定性带**一致地窄于参考**，即低估了真实条件不确定性。论文的解释是集成离散度衡量的是确定性预测器之间的差异，而不是条件律 $p(u\mid\mathcal O)$。这是全篇关于不确定性质量的主要论证，且它是**定性的**——没有给出覆盖率之类的定量标定表；除图 A.1 外是否另有定量比较，本次未能核实。

**代价。** 论文自己记下 UQ-SONet 在训练与推断上都比 VIDON 贵（多一个编码器网络；推断时要采样隐变量并统计）。表 A.2 在 NVIDIA RTX 3090 上给出：

| 算例   | UQ-SONet 训练 | VIDON 训练 |
| ------ | ------------- | ---------- |
| 4.1 节 | 6.03 h        | 5.24 h     |
| 4.2 节 | 4.37 h        | 4.04 h     |
| 4.3 节 | 6.01 h        | —          |
| 4.4 节 | 11.34 h       | —          |
| 4.5 节 | 10.32 h       | —          |

推断 0.026–0.031 s，VIDON 0.012–0.025 s。

**论文自陈的局限。** 结论一节明说三条：条件高斯解码器「对强非高斯或多峰后验可能过于受限」；长时间滚动预测下不确定性的累积未被处理；**输出观测也稀缺**的设定未被覆盖。此外传感器布点的主动学习被提出但未做。

### 与其他论文的关系

与编号 89 是直接的方法学同胞（同为 Ling Guo / Hao Wu / Tao Zhou 组合）：都是带高斯解码器与对角高斯变分后验、按证据下界训练的变分自编码器式隐变量模型；差别在隐变量买到了什么——编号 89 用扩散过程替掉编码器以**收紧**一个反向 KL 界从而采样，这里则是常规的摊销高斯编码器，隐变量承载的是**预测**不确定性。

与编号 98 同样是直接同胞：两者都给神经算子解场附上不确定性，都用编码器产生随机隐变量再解码成解的条件高斯。编号 98 把标准高斯隐先验换成由学得的置信函数混合的高斯过程先验；编号 95 保留 $\mathcal N(0,I)$ 而把**条件化**做成置换不变。合起来正好覆盖两条自然的轴：更丰富的先验（98）对更丰富的条件（95）。

编号 95 的引言明确引用了「基于信息瓶颈的不确定性量化框架……用于函数回归与算子学习」，即编号 75——那是重叠作者群的前作，而不是本文的一部分。它与编号 94、107 共享 DeepONet 底座；与编号 107 意图最近：两者都学一个**分布值**对象，编号 107 用条件正规化流学转移密度，编号 95 用条件变分自编码器学解场。它与[[computational-mathematics/paper-notes/scientific-machine-learning/variational-and-basis-networks|最小二乘配点求解器]]（编号 90、102）及编号 101、103 成对照：那些方法在求解时需要方程且给出一个确定性场，编号 95 纯数据驱动、推断时不需要方程、给出的是分布。

## 98：把噪声项换成高斯过程先验

### 直觉

这篇论文最值得记录的是它 3.1 节给出的理由，而不是它换了哪个模块。

设想用一个隐变量模型给 PDE 解配上不确定性。最省事的做法是让每个点 $\bm x$ 独立地带一份噪声。问题在于解并不是这样生成的。考虑 $\lambda u_{xx}=f$：解是 $f$ 的**二重积分**，因此（论文明确标注为启发式的关系）

$$
\mathrm{Var}\bigl(u(x)\bigr)\ \sim\ \iint_{\Omega_x}\mathrm{Var}\bigl(f(\xi)\bigr)\,\mathrm d\xi .
$$

也就是说 $u$ 的不确定性是 $f$ 的不确定性在整个区域上的**全局累积**，而不是它的逐点映射。一个低维独立同分布的隐变量表达不了这种结构，而一个带空间相关性的隐**过程**可以。**这是把 $\mathcal N(0,I)$ 换成高斯过程的实际论据**，而不是「相关结构更符合实际」这样的一般说法。

论文对已有路线的批评也各自具体：普通高斯过程回归对非线性方程吃力且规模化差；带 Hamiltonian Monte Carlo 的 Bayesian PINN 准确但后验计算昂贵；深度集成便宜，但在数据稀疏或含噪时会产生伪振荡与不可靠的置信度。论文的总结是，多数已有方法「要么缺少表示隐函数不确定性的结构化机制，要么计算负担仍然很重」。

### 问题设定

对 $\mathcal N_{\bm x}(u;\lambda)=f$ 于 $\Omega$、$\mathcal B_{\bm x}(u;\lambda)=b$ 于 $\partial\Omega$，PINN 基线取残差 $r(\bm x;\theta)=\mathcal N_{\bm x}(u_{\mathrm{NN}}(\bm x;\theta);\lambda_{\mathrm{NN}}(\bm x;\theta))-f(\bm x)$，损失 $\mathcal L(\theta)=w_r\mathcal L_r+w_b\mathcal L_b+w_d\mathcal L_d$，三项都是半平方和。观测是含噪的：$\bar u^{(i)}=u(\bm x_u^{(i)})+\epsilon_u^{(i)}$，$\bar f^{(i)}=f(\bm x_f^{(i)})+\epsilon_f^{(i)}$，$\bar b^{(i)}=b(\bm x_b^{(i)})+\epsilon_b^{(i)}$。

### 推导

**第一步：置信感知编码器（其式 6）。**

$$
\bm z(\bm x;\theta_E)=\mathrm{diag}\bigl(m(\bm x;\theta_m)\bigr)\,\bar{\bm z}(\bm x;\theta_{\bar{\bm z}})
+\mathrm{diag}\bigl(1-m(\bm x;\theta_m)\bigr)\,\bm z_0,
$$

其中 $\bar{\bm z}$ 是确定性深网络、$m$ 的每个分量取值于 $[0,1]$，而 $\bm z_0:\mathbb R^d\to\mathbb R^{d_z}$ 是**向量值高斯过程**

$$
\bm z_0\sim\mathcal{GP}\bigl(0,K(\bm x,\bm x')\bigr),
\qquad
K(\bm x,\bm x')=\sigma_K^2\exp\Bigl(-\frac{\|\bm x-\bm x'\|^2}{2\ell^2}\Bigr)\bm I_{d_z},
$$

各输出维独立并共用一个标量平方指数核（全程取 $\sigma_K=1$）。$m(\bm x)=1$ 表示 $\bm x$ 靠近训练数据、$\bm z$ 是确定的；$m(\bm x)=0$ 表示远离数据、$\bm z$ 承接全部高斯过程方差。由 $\mathrm{Var}(\bm z(\bm x))=(1-m(\bm x))^2$，论文指出 $m$ 也可以读作**用非平稳核参数化一个高斯过程**的另一种写法，诱导的编码器分布为

$$
q_E(\bm z\mid\bm x)\sim\mathcal N\Bigl(m(\bm x)\bar{\bm z}(\bm x),\ \operatorname{diag}\bigl(1-m(\bm x)\bigr)^2\Bigr).
$$

与编号 95、89 的差别正在这里：隐先验是定义在输入域上的过程而非 $\mathcal N(0,I)$，因此隐随机性带空间相关；而混合权重逐点可学，所以模型能在见过数据处自信、在别处不自信。

**第二步：神经算子解码器（其式 7–8、10）。** 隐场经一个 FNO 式积分算子栈传播，$\bm z_1(\bm x,\bm\omega_E)=\bm z(\bm x,\bm\omega_E)$、$\bm\omega_E\sim\mathcal N(0,\bm I_M)$：

$$
\bm z_i=\sigma\Bigl(W_i\bm z_{i-1}+b_i+\int_{\Omega}k_i(\bm x,\bm x')\bm z_{i-1}(\bm x',\bm\omega_E)\,\mathrm d\bm x'\Bigr),
\quad i=2,\dots,L-1,
\qquad
\bm z_L=W_L\bm z_{L-1}+b_L,
$$

核取「位置 Transformer」形式

$$
k_i(\bm x,\bm x')=\frac{\exp(-\alpha_i\|\bm x-\bm x'\|^2)}{\int_{\Omega}\exp(-\alpha_i\|\bm x-\bm y\|^2)\mathrm d\bm y}\,V_i,
$$

$\alpha_i>0$ 为长度尺度、$V_i$ 可学，积分用求积、FFT 或 Monte Carlo 计算。输出 $\mu_u(\bm x,\bm\omega_E):=\bm z_L(\bm x,\bm\omega_E)$ 是 $u$ 的预测均值。**这一层积分是「全局累积」在架构上的对应物**：$\mu_u$ 在 $\bm x$ 处的值依赖整个区域上的隐场，而不只是 $\bm z(\bm x)$。

**第三步：概率解码器与两类不确定性的分离（其式 9–10）。**

$$
q_D^u(u\mid\bm x,\bm\omega_E)\sim\mathcal N\bigl(\mu_u(\bm x,\bm\omega_E),\ \sigma_u^2(\bm x;\theta_\sigma^u)\bigr),
\qquad
u(\bm x,\bm\omega_u)=\bm z_L(\bm x,\bm\omega_E)+\sigma_u(\bm x;\theta_\sigma^u)\cdot\omega_D^u,
$$

$\omega_D^u\sim\mathcal N(0,1)$。于是 $\sigma_u$ 建模**偶然**不确定性、$\bm\omega_E$（经 $m$ 与高斯过程）建模**认知**不确定性，两者在结构上被分开。

**第四步：物理作为硬性结构恒等式（其式 11–12）。** 论文假定 $\sigma_u(\bm x;\theta_\sigma^u)$ 不参与对 $\bm x$ 的求导，于是

$$
\mathcal N_{\bm x}[u(\bm x,\bm\omega_u)]=\mathcal N_{\bm x}[\mu_u(\bm x,\bm\omega_E)],
\qquad
\mu_f=\mathcal N_{\bm x}[\mu_u],
\qquad
\mu_b=\mathcal B_{\bm x}[\mu_u],
$$

$$
f(\bm x,\bm\omega_f)=\mathcal N_{\bm x}[\mu_u]+\sigma_f(\bm x;\theta_\sigma^f)\omega_D^f,
\qquad
b(\bm x,\bm\omega_b)=\mathcal B_{\bm x}[\mu_u]+\sigma_b(\bm x;\theta_\sigma^b)\omega_D^b .
$$

也就是说物理约束是一个**硬性的结构恒等式**，而不是另加一个残差网络去拟合——这与编号 75 把物理放在损失项里的做法不同。链式法则要穿过积分层，而对 $\bm z_1$ 求导又要对高斯过程样本 $\bm z_0$ 求导；附录 A 给出两条路：**精确的 Karhunen–Loève 展开**，与一个**有限维联合高斯采样恒等式**。

**第五步：数据似然。** $\mathcal L_{\mathrm{data}}=\mathcal L_{\mathrm{data},u}+\mathcal L_{\mathrm{data},f}+\mathcal L_{\mathrm{data},b}$，每项形如

$$
\mathcal L_{\mathrm{data},u}=\frac{1}{N_uN_{\bm\omega}}\sum_{i=1}^{N_u}\sum_{j=1}^{N_{\bm\omega}}
\log q_D^u\bigl(\bar u^{(i)}\bigm|\bm x_u^{(i)},\bm\omega_E^{(j)};\theta_u\bigr),
\qquad \bm\omega_E^{(j)}\sim\mathcal N(0,\bm I_M).
$$

**第六步：信息瓶颈正则项——为阻止隐变量塌缩（其式 13–15）。** 论文观察到只用数据损失时模型会把观测点处的 $m(\bm x_i)$ 推到 $1$，导致 $\bm z$ 在整个区域退化为确定映射。仿照编号 75 加上

$$
\mathcal L_{\mathrm{reg}}(Z,X)=\mathbb E_{\bm x\sim\mathcal U(\Omega)}\,
D_{\mathrm{KL}}\bigl[q_E(\bm z\mid\bm x)\,\|\,e_1(\bm z)\bigr],
\qquad e_1(\bm z)\sim\mathcal N(0,\bm I_{d_z}),
$$

经验估计为 $\mathcal L_{\mathrm{reg}}=\frac1N\sum_i[\log q_E(\bm z_i\mid\bm x_i)-\log e_1(\bm z_i)]$、$\bm x_i\sim\mathcal U(\Omega)$，总目标为

$$
\mathcal L=\mathcal L_{\mathrm{data}}-\beta\cdot\mathcal L_{\mathrm{reg}} .
$$

> [!warning] 符号约定与叙述不一致
> 论文**极大化** $\mathcal L$（算法 1 用 $W=W+lr\cdot\partial\mathcal L/\partial W$），且 $\mathcal L_{\mathrm{data}}$ 写作 $+\log q_D$，所以它是对数似然而非负对数似然，尽管周边文字称之为「负对数似然损失」。另外，全程「非贝叶斯地」训练，不做后验推断。

**第七步：带相关性的正则项（其注记 3.3，式 16–17）——实验中实际使用的版本。** 因为 $\mathcal L_{\mathrm{reg}}$ 忽略了输入之间的相互作用，论文定义一个对整批的版本，参照分布取高斯过程 $e_2(\bm z_{1:B}\mid\bm x_{1:B})\sim\mathcal N(0,\bm K)$、$\bm K=(K(\bm x_i,\bm x_j))_{i,j\le B}$，闭式为

$$
\widetilde{\mathcal L}_{\mathrm{reg}}=\mathbb E_{\bm x_{1:B}\sim\mathcal U(\Omega^B)}\frac{1}{2B}
\Bigl[-2\log\det(\bm D_m)-B+\mathrm{tr}(\bm K^{-1}\bm D_m\bm K\bm D_m)
+(\bm{m\bar z})^{\top}\bm K^{-1}(\bm{m\bar z})\Bigr],
$$

其中 $\bm D_m=\operatorname{diag}(1-m(\bm x_1),\dots,1-m(\bm x_B))$、$\bm{m\bar z}=(m(\bm x_1)\bar z(\bm x_1),\dots)$。极小化它会**同时**优化 $m$ 与高斯过程的超参数。论文明说「$B\to\infty$ 时 $\widetilde{\mathcal L}_{\mathrm{reg}}$ 的统计性质尚待理论分析」——这是它自陈的一处空白。

**第八步：反问题扩展（其式 18–20）。** 加一个解码器 $q_D^\lambda(\lambda\mid\bm x,\bm\omega_E;\theta_\lambda)\sim\mathcal N(\mu_\lambda,\sigma_\lambda^2)$，**与 $u$ 共用同一份编码器随机性 $\bm\omega_E$**；$\lambda$ 是标量或向量而非场时取 $\mu_\lambda(\bm\omega_E)=\int_\Omega\mu_\lambda(\bm x,\bm\omega_E)\mathrm d\bm x$。于是 $f=\mathcal N_{\bm x}[\mu_u;\mu_\lambda]+\sigma_f\omega_D^f$、$b=\mathcal B_{\bm x}[\mu_u;\mu_\lambda]+\sigma_b\omega_D^b$，总目标加一项 $\mathcal L_{\mathrm{data},\lambda}$。

**解码器的替代方案（其注记 3.1）。** 也可用 DeepONet 式解码器 $\mu_u(\bm x,\bm\omega_E)=\mathrm{Branch}(\bm z(\bm x_{\mathrm{grid},1},\bm\omega_E),\dots,\bm z(\bm x_{\mathrm{grid},p},\bm\omega_E))\cdot\mathrm{Trunk}(\bm z(\bm x,\bm\omega_E))$，但论文指出分支输入维数 $p\times d_z$ 可能极大。FNO 式积分朴素实现是 $\mathcal O(N^2)$，用 FFT 是 $\mathcal O(N\log N)$。该变体的结果在附录 B。

### 定理

**没有定理、没有命题、没有误差分析。** 论文在结论一节自己写明：「稳定性与精度的严格分析仍然缺失，修改后的正则项的数学适定性也尚未完全建立。」自陈的后续工作包括理论、时间依赖问题，以及把该框架并入物理信息**算子**学习。

### 数值实验

几乎全部算例都在两个噪声水平 $\mathcal N(0,0.01^2)$ 与 $\mathcal N(0,0.1^2)$ 上测，对比对象始终是 B-PINN-HMC 与深度集成，且**网络架构完全相同**。

**例 1：一维 Poisson（正问题）。** $\lambda\partial_x^2u=f$ 于 $[-0.7,0.7]$，$\lambda=0.01$，精确解 $u(x)=\sin^3(6x)$；32 个等距 $f$ 传感器加 2 个边界 $u$ 传感器；$\sigma_f,\sigma_u$ 是可学标量，**初值取经验噪声水平的两倍**。低噪声下三种方法都好；$\sigma=0.1$ 时深度集成对 $f$ 给不出可靠的不确定性，而 LVM-GP 的均值更好且不确定性带比 B-PINN-HMC 更**紧**。启用带相关性的正则项与可学核长度尺度 $L_c$ 后：**$L_c$ 从约 $1.0$ 出发，在约 8000 次迭代后稳定在 $0.58$ 附近**——这是模型自行调整先验相关尺度的一个具体证据。

**例 2：带边界层的一维多孔介质流（正问题）。** $-\frac{\nu_e}{\phi}\partial_x^2u+\frac{\nu u}{K}=g$ 于 $[0,1]$，$u(0)=u(1)=0$，解析解

$$
u=\frac{gK}{\nu}\Bigl[1-\frac{\cosh\bigl(r(x-H/2)\bigr)}{\cosh(rH/2)}\Bigr],
\qquad r=\sqrt{\nu\phi/(\nu_eK)},
$$

参数 $\nu_e=\nu=10^{-3}$、$\phi=0.4$、$K=10^{-3}$、$g=1$；16 个 $g$ 传感器加 2 个边界传感器；编码器与解码器都是 3 层全连接、每隐层 20 个神经元。结论：LVM-GP 与 HMC 相当；深度集成在高噪声下既不准、不确定性也不可靠。训练用两段式：**先只优化预测均值 2000 步，再优化标准差并微调均值 8000 步**。

**例 3：一维非线性 Poisson（反问题）。** $k\partial_x^2u+\lambda\tanh(u)=f$ 于 $[-0.7,0.7]$，$k=0.01$，精确解 $u=\sin^3(6x)$，未知反应率 $\lambda=0.7$；32 个 $f$ 传感器、2 个边界 $u$ 传感器、6 个内部 $u$ 传感器；$\Lambda$ 网络是单神经元隐层加 $\tanh$。

| 噪声 | LVM-GP 均值（标准差）          | B-PINN-HMC                     | 深度集成                       |
| ---- | ------------------------------ | ------------------------------ | ------------------------------ |
| 0.01 | 0.6976（$9.816\times10^{-3}$） | 0.6967（$4.225\times10^{-3}$） | 0.6966（$2.493\times10^{-4}$） |
| 0.1  | 0.6965（$6.954\times10^{-2}$） | 0.6787（$4.166\times10^{-2}$） | 0.6959（$3.691\times10^{-2}$） |

读法：高噪声下 LVM-GP 的均值最准，但报告的标准差也**最大**——与论文「竞争方法过度自信」的论点一致，尽管论文没有显式计算标定。

**外推测试（其图 9）。** 在**缺失单侧边界条件**的设定下，只给 $x\le0$ 处的 4 个 $u$ 测量值加 40 个均匀 $f$ 测量值（噪声标准差 0.01），模型靠 PDE 约束把 $u$ 外推到 $x>0$，参数恢复为均值 $0.7040$、标准差 $1.361\times10^{-2}$。

> [!warning] 原文的一处符号笔误
> 论文文字称这是「$k$ 的预测均值」，但该方程中 $k=0.01$ 而 $\lambda=0.7$，因此 $0.7040$ 只能指 $\lambda$。此处按原样记录并标出，而不做无声更正。

**例 4：二维非线性扩散反应（反问题）。** $k(\partial_{x_1}^2u+\partial_{x_2}^2u)+\lambda u^2=f$ 于 $[-1,1]^2$，$k=0.01$，精确解 $u=\sin(\pi x_1)\sin(\pi x_2)$，未知 $\lambda=1$；100 个内部 $u$ 传感器、484 个 $f$ 传感器、每边 25 个边界 $u$ 传感器；3 层全连接、每隐层 128 个神经元。

| 噪声 | LVM-GP 均值（标准差）         | B-PINN-HMC                    | 深度集成                      |
| ---- | ----------------------------- | ----------------------------- | ----------------------------- |
| 0.01 | 1.0003（$4.58\times10^{-3}$） | 1.0005（$5.75\times10^{-3}$） | 1.0047（$4.12\times10^{-3}$） |
| 0.1  | 0.9916（$5.70\times10^{-3}$） | 0.9781（$4.98\times10^{-2}$） | 0.9302（$2.60\times10^{-2}$） |

**这是全篇最强的定量结果**：高噪声下 LVM-GP 既最准，标准差又比两个基线小一个数量级。

**例 5：六维污染源反演。** $-\lambda(\partial_{x_1}^2+\partial_{x_2}^2)u-f_2=f_1$ 于 $[0,1]^2$、零边界条件，$\lambda=0.02$，已知 $f_1=0.1\sin(\pi x_1)\sin(\pi x_2)$，未知源

$$
f_2=\sum_{i=1}^{3}k_i\exp\Bigl[-0.5\frac{\|\bm x-\bm x_{c,i}\|^2}{0.15^2}\Bigr],
\qquad \bm k=(2,-3,0.5)\ \text{已知},
$$

要反演三个中心（共 6 个未知量）。真值为 $(0.3,0.3)$、$(0.75,0.75)$、$(0.2,0.7)$；噪声 $\epsilon_u\sim\mathcal N(0,0.1^2)$、$\epsilon_{f_1}\sim\mathcal N(0,0.01^2)$；参考 $u$ 用 FEniCS 有限元算；1000 个随机 $u$ 样本、200 个 $f_1$ 样本；3 层 128 神经元的编码器与解码器、激活取 **Mish**，2 层 128 神经元的 $\Lambda$ 网络、激活取 $\tanh$；$\sigma_u$ 初值 0.1；2 万次 Adam 迭代、学习率固定 $0.001$，按 1 万（只优化均值）加 1 万（均值与标准差）划分。

| 中心 | 真值          | LVM-GP 均值       | LVM-GP 标准差（$\times10^{2}$） | B-PINN-HMC 均值   | HMC 标准差（$\times10^{3}$） |
| ---- | ------------- | ----------------- | ------------------------------- | ----------------- | ---------------------------- |
| 1    | $(0.3,0.3)$   | $(0.2927,0.3022)$ | $(2.79,6.25)$                   | $(0.3014,0.2883)$ | $(3.08,3.45)$                |
| 2    | $(0.75,0.75)$ | $(0.7433,0.7542)$ | $(2.87,2.36)$                   | $(0.7473,0.7496)$ | $(3.51,2.52)$                |
| 3    | $(0.2,0.7)$   | $(0.2065,0.7569)$ | $(5.49,6.34)$                   | $(0.2268,0.6519)$ | $(18.97,11.47)$              |

**注意第三个源**：B-PINN-HMC 给出 $0.6519$（真值 $0.7$）却报告很小的标准差，而 LVM-GP 给出 $0.7569$ 并配一个更大的标准差——又一次是过度自信的模式。两行标准差在原文中用的是**不同的**比例因子（$\times10^2$ 与 $\times10^3$），此处照录。

这五组实验合起来建立的是：在噪声与稀疏数据下，LVM-GP 的均值精度可与 B-PINN-HMC 比肩而不需后验采样，且其不确定性宽窄的**方向**更合理。它们没有建立的是标定质量本身——论文全程用「标准差大小是否与误差匹配」做定性判断，没有算覆盖率、区间得分或 CRPS。

### 与其他论文的关系

与编号 95 是同胞（都属 Ling Guo / Hao Wu / Tao Zhou 一系）：两者都是「编码器–隐变量–解码器」结构，解码器给出解场上的条件高斯，均值由神经算子给出。分工很干净——编号 98 丰富**先验**（高斯过程取代 $\mathcal N(0,I)$，再配一个学得的置信混合），并以信息瓶颈式 KL 正则项非贝叶斯地训练；编号 95 保留标准高斯先验而丰富**条件**（对传感器做置换不变的集合 transformer），并以正规的证据下界训练。

与编号 89 在机器上同源（隐变量、KL 散度、重参数化采样）而目的相反：编号 89 从已知能量采样，编号 98 从含噪数据回归 PDE 解。信息瓶颈目标 $I(Z;Y)-\beta I(Z;X)$ 在本文 2.2.2 节被显式写出并引到编号 75，本文的 $\mathcal L=\mathcal L_{\mathrm{data}}-\beta\mathcal L_{\mathrm{reg}}$ 就是它的一个变分实例。

与[[computational-mathematics/paper-notes/scientific-machine-learning/variational-and-basis-networks|确定性 PDE 求解器一页]]（编号 90、102）以及编号 101、103 构成对照：那些方法假定 PDE 数据精确，用最小二乘把残差压到机器精度；编号 98 面对的是**含噪、稀疏**的数据，此时把残差压到零反而是过拟合，所以物理被放进似然里当作软约束。

## 107：把学习对象从初值依赖换成转移密度

### 直觉

在集合预报与数据同化中，随机动力学是固定的，而**初始分布随每个场景或后验更新而变**。若把「初始分布 $\mapsto$ 时间依赖解」直接当作要学的算子，那么每换一族初值就要重新准备一批训练数据。

这篇论文换了一个学习对象：不学初值依赖，而学**转移密度** $p(\bm x,t\mid\bm x_0)$——从一个确定的初始点出发、经过时间 $t$ 之后落在 $\bm x$ 的密度。它只依赖动力学本身，与初始分布无关。有了它，任何初始分布的解都由 Chapman-Kolmogorov 方程一次积分得到，完全不需要重训。

代价立刻出现：转移密度满足一个带 **Dirac 初值**的方程。正规化流的输出密度总是光滑的，任何可逆映射都不可能把一个光滑基分布送到 $\delta(\bm x-\bm x_0)$。论文的解法是把奇性交给**基分布**：把漂移在 $\bm x_0$ 处线性化、扩散取常值，得到的线性随机微分方程有闭式高斯解，且这个高斯在 $t\to0$ 时自动塌成 $\delta_{\bm x_0}$。于是流只需要在这个已经近似正确的基分布上学一个**接近恒等**的修正。

同一思路再走一步就得到第二个技巧：让耦合层的每个非线性项都显式带一个因子 $t$，于是层在 $t=0$ 时**恰好**是恒等映射，初始条件不用任何惩罚项就成立。第三个技巧是时间加权：即使近似再好，$t\to0$ 时 PDE 残差本身也会发散，所以损失要按发散速率反向加权。

### 问题设定

设非退化随机微分方程

$$
\mathrm d\bm X_t=\bm f(\bm X_t)\mathrm dt+\bm g(\bm X_t)\mathrm d\bm W_t,
\qquad \bm f:\mathbb R^d\to\mathbb R^d,\quad \bm g:\mathbb R^d\to\mathbb R^{d\times m},\ m\ge d,
$$

其密度满足 Fokker-Planck 方程

$$
\frac{\partial p}{\partial t}=\mathcal L^{*}_{\bm f,\bm g}p
:=-\nabla\cdot(p\bm f)+\frac12\nabla\cdot\nabla\cdot(\bm g\bm g^{\mathrm T}p),
$$

并带 $p(\bm x)\to0$（$|\bm x|_2\to\infty$）、$\int_{\mathbb R^d}p(\bm x,t)\mathrm d\bm x\equiv1$ 与 $p\ge0$。论文指出直接把 DeepONet、Fourier 神经算子或 Green 函数型算子学习用于 Fokker-Planck 方程「需要对边界条件与概率密度约束（非负性与单位质量）做额外处理」。

### 推导

**第一步：改学转移密度。** 不再对一整类初值 $p_0\in\mathcal G$ 逐个求解，而是一次性求

$$
\frac{\partial p(\bm x,t\mid\bm x_{0})}{\partial t}
=\mathcal L^{*}_{\bm f,\bm g}\,p(\bm x,t\mid\bm x_{0}),
\qquad
p(\bm x,0\mid\bm x_{0})=\delta(\bm x-\bm x_{0}),
$$

再由 Chapman-Kolmogorov 方程恢复任意初值：

$$
p(\bm x,t)=\int_{\Omega_{0}}p(\bm x,t\mid\bm x_{0})\,p_{0}(\bm x_{0})\,\mathrm d\bm x_{0} .
$$

**这是全篇的核心动作：被学习的对象与 $p_0$ 无关，因此新初值完全不需要重训。**

**第二步：条件正规化流。** 由换元公式 $p_{\bm X}(\bm x)=p_{\bm Z}(T(\bm x))|\det\nabla_{\bm x}T(\bm x)|$，令变换 $T$ 与基分布同时依赖 $(\bm x_0,t)$：

$$
\bm X_t|_{\bm X_0=\bm x_0}=T^{-1}\bigl(\bm Z_t(\bm x_0);\theta(\bm x_0,t)\bigr),
\qquad
p_{\bm X_{t}\mid\bm X_{0}=\bm x_{0}}(\bm x)
=p_{\bm Z_{t}(\bm x_{0})}\bigl(T(\bm x;\theta(\bm x_{0},t))\bigr)\,
\bigl|\det\nabla_{\bm x}T(\bm x;\theta(\bm x_{0},t))\bigr| .
$$

由于是流，非负性与单位质量**由构造成立**——这正是论文对「概率密度约束」这一反对意见的回答，与[[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|密度流一页]]的取向一致。

**第三步：线性化基过程。** 对漂移做一阶 Taylor、对扩散做零阶 Taylor，都在 $\bm x_0$ 处展开，

$$
\mathrm d\widehat{\bm X}_{t}
=\bigl(\bm f(\bm x_{0})+\nabla\bm f(\bm x_{0})(\widehat{\bm X}_{t}-\bm x_{0})\bigr)\mathrm dt
+\bm g(\bm x_{0})\mathrm d\bm W_{t},
\qquad \widehat{\bm X}_{0}=\bm x_{0} .
$$

这是一个 Ornstein-Uhlenbeck 型线性方程，解可以显式写出，

$$
\widehat{\bm X}_{t}=\bm x_{0}
+\int_{0}^{t}e^{\nabla\bm f(\bm x_{0})(t-s)}\bm f(\bm x_{0})\,\mathrm ds
+\int_{0}^{t}e^{\nabla\bm f(\bm x_{0})(t-s)}\bm g(\bm x_{0})\,\mathrm d\bm W_{s},
$$

因而是高斯过程，其均值与协方差有闭式：

$$
\mathbb E\bigl[\widehat{\bm X}_{t}\bigm|\widehat{\bm X}_0=\bm x_0\bigr]
=\bm x_{0}+\int_{0}^{t}e^{\nabla\bm f(\bm x_{0})(t-s)}\bm f(\bm x_{0})\,\mathrm ds,
$$

$$
\bm\Sigma\bigl[\widehat{\bm X}_{t}\bigm|\widehat{\bm X}_0=\bm x_0\bigr]
=\int_{0}^{t}e^{\nabla\bm f(\bm x_{0})(t-s)}\bm g(\bm x_{0})\bm g(\bm x_{0})^{\top}
e^{\nabla\bm f(\bm x_{0})^{\top}(t-s)}\,\mathrm ds .
$$

两个积分都用 **Gauss-Legendre 求积**计算：取 $[-1,1]$ 上的标准节点与权 $\{(\hat w_i,\hat s_i)\}$，映射为 $w_i=\tfrac t2\hat w_i$、$s_i=\tfrac t2(\hat s_i+1)$。$\nabla\bm f(\bm x_0)$ 可逆时均值还有闭式 $\bm x_0+(\nabla\bm f(\bm x_0))^{-1}(e^{\nabla\bm f(\bm x_0)t}-\bm I)\bm f(\bm x_0)$。**这个基过程在 $t\to0$ 时协方差趋零、均值趋 $\bm x_0$，自动退化为 $\delta_{\bm x_0}$**，因此流不需要从光滑分布跨越到 Dirac：奇性被基分布本身吸收了。

**第四步：C-KRnet 架构。** 流写成复合 $T=T_{[K]}\circ\cdots\circ T_{[1]}$，行列式按 $|\det\nabla_{\bm x}T|=\prod_i|\det\nabla_{\bm x_{[i-1]}}T_{[i]}|$ 分解，采用**主动维数递减**的 KRnet 结构（Knothe-Rosenblatt 重排加仿射耦合层）：

$$
T_{[1]}=\tilde T_{K},
\qquad
T_{[k]}=\begin{pmatrix}\tilde T_{K+1-k}\\ \mathrm{Id}_{K+2-k:K}\end{pmatrix},
\quad k=2,\dots,K,
$$

对应分块 $\bm x^{\top}=((\bm x^{(1)})^{\top},\dots,(\bm x^{(K)})^{\top})$、$\sum_id_i=d$。外层是 $K$ 个阶段，内层 $\tilde T_k$ 是 $l_k$ 个耦合层；$\tilde T_k$ 之后第 $k$ 块被冻结，所以后面的阶段可以用更小的 $l_k$。

**第五步：在 $t=0$ 处内建恒等的条件耦合层（其式 15–16）。**

$$
\hat{\bm y}_{1}=\bm y_{1},
\qquad
\hat{\bm y}_{2}=\bm y_{2}\odot\bigl(\bm 1_{m-m_{1}}+\beta\tanh(t\cdot\bm s)\bigr)
+e^{\bm\zeta}\odot\tanh(t\cdot\bm q),
\qquad
(\bm s,\bm q)=\mathrm{NN}(\bm y_{1},\bm\xi,t),
$$

$\beta\in(0,1)$ 由用户设定、$\bm\zeta$ 可训练。因为两个 $\tanh$ 的宗量都显式带因子 $t$，**该层在 $t=0$ 时恰好是恒等映射**，从而与基分布对齐、无需任何惩罚项就满足初始条件。$\bm y_1,\bm y_2$ 的角色在层间交替，使所有坐标都被更新。

**第六步：随机 Fourier 特征条件网络（其式 17）。**

$$
\bm h_{0}=[\bm y_{1},\bm\xi,t]^{\top},
\qquad
\bm h_{1}=\Bigl[\sin\bigl(\tfrac{1}{e^{\gamma}}\bm F\bm h_{0}+\bm b_{0}\bigr),\
\cos\bigl(\tfrac{1}{e^{\gamma}}\bm F\bm h_{0}+\bm b_{0}\bigr),\ \bm h_{0}\Bigr]^{\top},
$$

$$
\bm h_{j}=\mathrm{SiLU}(\bm W_{j-1}\bm h_{j-1}+\bm b_{j-1}),\ j=2,\dots,M,
\qquad
(\bm s,\bm q)=\bm W_{M}\bm h_{M}+\bm b_{M},
$$

其中 $\bm F\in\mathbb R^{r_h/2\times\dim(h_0)}$ 抽自标准正态、$\bm b_0\in\mathbb R^{r_h/2}$ 抽自 $[0,2\pi]^{r_h/2}$ 上的均匀分布，二者**初始化后固定**；$\gamma$ 与 $\{\bm W_j,\bm b_j\}$ 可训练。声明的动机是加强初始状态 $\bm x_0$ 的影响。注意频率尺度 $1/e^{\gamma}$ 是可学的。

**第七步：残差损失与自适应采样。** 损失只有残差一项（其式 18）：

$$
\mathcal L(\bm\theta)=\int_{\mathbb R^{d}\times(0,t_{f}]\times\Omega_{0}}
|r(\bm x,\bm x_{0},t;p_{\bm\theta})|^{2}\,\mathrm d\rho(\bm x,\bm x_{0},t),
\qquad
r=\frac{\partial p_{\bm\theta}}{\partial t}-\mathcal L^{*}_{\bm f,\bm g}[p_{\bm\theta}],
$$

采样测度分解为 $\rho(\bm x,\bm x_0,t)=\rho(\bm x\mid\bm x_0,t)\rho(t\mid\bm x_0)\rho(\bm x_0)$，其中 $\rho(t\mid\bm x_0)\rho(\bm x_0)$ 在 $\Omega_0\times[0,t_f]$ 上均匀。**没有初始条件损失项**——它由结构满足。空间部分的测度按阶段更新（其式 19）：

$$
\rho_{0}(\bm x\mid\bm x_{0},t)=\frac{\gamma_{1}}{\gamma_{1}+\gamma_{3}}\mu(\bm x)
+\frac{\gamma_{3}}{\gamma_{1}+\gamma_{3}}\rho_{\bm Z}(\bm x\mid\bm x_{0},t),
$$

$$
\rho_{k+1}(\bm x\mid\bm x_{0},t)=\gamma_{1}\mu(\bm x)+\gamma_{2}\rho_{k}(\bm x\mid\bm x_{0},t)
+\gamma_{3}\rho_{\text{C-KRnet},\bm\theta}(\bm x\mid\bm x_{0},t),
\qquad \gamma_{1}+\gamma_{2}+\gamma_{3}=1,
$$

混合一个用户给定盒上的均匀测度 $\mu$、上一阶段的测度，以及**当前**流诱导的测度。第 0 阶段从**基过程** $\rho_{\bm Z}$ 而非均匀起步，正因为基分布在小 $t$ 处已经准确。

**第八步：时间加权。** 命题 4.1 说残差在 $t\to0$ 时按两个不同速率发散，而这两个速率对应**两种不同的采样测度**。于是经验损失用两个不同的指数（其式 22）：均匀采样部分配权 $w_1(t)=t^{\frac d2+2}$，两个解采样部分配权 $w_2(t)=t^{d+2}$，

$$
\widehat{\mathcal L}_{w}^{k}(\bm\theta;S)
:=\frac{1}{N}\sum_{i=1}^{N}\bigl|r(\bm x^{i},\bm x^{i}_{0},t^{i};p_{\bm\theta})\bigr|^{2}
\Bigl(\eta^{i}t^{d+2}+(1-\eta^{i})t^{\frac d2+2}\Bigr),
$$

其中 $N=N_1+N_2+N_3$、$N_1:N_2:N_3=\gamma_1:\gamma_2:\gamma_3$，指示量 $\eta^i=0$ 对应均匀样本 $S_1$、$\eta^i=1$ 对应 $S_2,S_3$。论文把这一加权描述为「在保持因果性与缓解优化困难之间的折中」——指数恰好抵消命题 4.1 的发散速率，使加权残差为 $\mathcal O(1)$。

**第九步：评价阶段的重要性采样（其式 24–27）。** 朴素 Monte Carlo $p(\bm x,t)\approx\frac1M\sum_ip(\bm x,t\mid\bm x_0^i)$、$\bm x_0^i\sim p_0$ 方差很大，因为小 $t$ 时转移密度几乎是 Dirac，绝大多数 $\bm x_0^i$ 对给定的 $\bm x$ 贡献近乎为零。改用

$$
p(\bm x,t)\approx\hat p_{\bm\theta}(\bm x,t)
:=\frac1M\sum_{i=1}^{M}\frac{p(\bm x,t\mid\bm x_{0}^{i})p_{0}(\bm x_{0}^{i})}{q(\bm x_{0}^{i}\mid\bm x,t)},
\qquad \bm x_{0}^{i}\sim q(\cdot\mid\bm x,t),
$$

提议分布取「时间反向线性化」的高斯

$$
q_{1}(\bm x_{0}\mid\bm x,t)=\mathcal N\bigl(\bm x_{0};m(\bm x,t),\tilde{\bm\Sigma}(\bm x,t)\bigr),
$$

$$
m(\bm x,t)=\bm x-\int_{0}^{t}e^{\nabla\bm f(\bm x)s}\bm f(\bm x)\,\mathrm ds,
\qquad
\tilde{\bm\Sigma}(\bm x,t)=\int_{0}^{t}e^{\nabla\bm f(\bm x)s}\bm g(\bm x)\bm g(\bm x)^{\top}
e^{\nabla\bm f(\bm x)^{\top}s}\,\mathrm ds .
$$

$q_1$ 随 $t$ 增大而变差，于是再取**按时间退火的混合**

$$
q_{2}(\bm x_{0}\mid\bm x,t)=\alpha(t)q_{1}(\bm x_{0}\mid\bm x,t)+(1-\alpha(t))p_{0}(\bm x_{0}),
\qquad \alpha(t)=\exp(-at),\ a>0,
$$

采样时从 $q_1$ 抽 $\alpha(t)M$ 个点、从 $p_0$ 抽 $(1-\alpha(t))M$ 个点。

**算法。** 外层是 $k=0,\dots,N_{\mathrm{adaptive}}$ 的自适应循环。$k=0$ 时取 $N_1=\lfloor\frac{\gamma_1}{\gamma_1+\gamma_3}N\rfloor$、$N_3=N-N_1$，$S_1$ 按 $\mu(\bm x)\rho(t\mid\bm x_0)\rho(\bm x_0)$ 生成并置 $\eta^i=0$，$S_3$ 从当前流生成并置 $\eta^i=1$，令 $S=S_1\cup S_3$；$k\ge1$ 时取 $N_1=\lfloor\gamma_1N\rfloor$、$N_2=\lfloor\gamma_2N\rfloor$、$N_3=N-N_1-N_2$，$S_1$ 重新均匀生成，$S_2$ 从**已有的 $S$** 中随机取 $N_2$ 个，$S_3$ 从当前流抽，更新 $S=S_1\cup S_2\cup S_3$。内层跑 $N_e$ 轮，每 $n_s$ 步把学习率乘以衰减率；每轮把 $S$ 随机切成小批量，按式 22 计损失并用 Adam 更新。评价时对任意新的 $p_0$ 从 $q_2$ 抽点、按重要性采样公式组装 $\hat p_{\bm\theta}$，不重训。

> [!warning] 符号复用
> 原文中 $m$ 同时表示扩散矩阵的列数与 Gauss-Legendre 求积的节点数，$\eta$ 同时表示学习率衰减率与损失中的样本指示量。按上下文可以区分，但阅读时需留意。

### 定理

四条命题，证明分别在附录 A–D。

**命题 3.1（线性化基过程的短时精度）。** 设 $\bm f=\bm f_\nu$ 可微、$\bm g$ **为常数**、$\bm G=\tfrac12\bm g\bm g^{\top}$ 非退化，且 $\|\nabla\bm f_\nu(\bm x)\|\le C_{\bm f}$ 对所有 $\bm x$ 成立。记 $\nu_t$、$\sigma_t$ 分别为真实过程与线性化过程从 $\delta_{\bm x_0}$ 出发的分布。则对 $t\in(0,t_f)$，

$$
H(\sigma_{t}\mid\nu_{t})\le\frac12\int_{0}^{t}\!\!\int_{\mathbb R^{d}}
\bigl|\bm G^{-1/2}(\bm f_{\sigma}(\bm x)-\bm f_{\nu}(\bm x))\bigr|^{2}\mathrm d\sigma_{s}\mathrm ds,
$$

$$
\|\nu_{t}-\sigma_{t}\|_{\mathrm{TV}}^{2}\le\frac14\int_{0}^{t}\!\!\int_{\mathbb R^{d}}
\bigl|\bm G^{-1/2}(\bm f_{\sigma}-\bm f_{\nu})\bigr|^{2}\mathrm d\sigma_{s}\mathrm ds,
$$

且 $t\to0$ 时 $H(\sigma_t\mid\nu_t)$ 与 $\|\nu_t-\sigma_t\|^2_{\mathrm{TV}}$ 都是 $\mathcal O(t^2)$（故全变差距离本身是 $\mathcal O(t)$）。**若再设 $\|\nabla^2\bm f_\nu(\bm x)\|\le H_{\bm f}$，二者都改进到 $\mathcal O(t^3)$。** 证明用 Bogachev 等人的相对熵界加 Pinsker-Csiszár-Kullback 不等式，中间估计为

$$
\mathrm{tr}(\bm\Sigma_{s})\le\mathrm{tr}(\bm G)\frac{e^{2C_{\bm f}s}-1}{C_{\bm f}}=\mathcal O(s),
\qquad
\bigl|\mathbb E_{\sigma_s}[\widehat{\bm X}_{s}]-\bm x_{0}\bigr|^{2}
\le|\bm f(\bm x_{0})|^{2}\Bigl(\frac{e^{C_{\bm f}s}-1}{C_{\bm f}}\Bigr)^{2}=\mathcal O(s^{2}).
$$

**命题 3.2（全变差下期望的传递）。** 对 $\mathbb R^d$ 上的概率测度 $\nu,\sigma$ 与 $\bm f\in L^2(\nu)\cap L^2(\sigma)$，

$$
\Bigl|\int\bm f\,\mathrm d\nu-\int\bm f\,\mathrm d\sigma\Bigr|
\le\sqrt{2\|\nu-\sigma\|_{\mathrm{TV}}}
\Bigl(\|\bm f\|_{L^{2}(\mathbb R^{d},\nu;\mathbb R^{m})}+\|\bm f\|_{L^{2}(\mathbb R^{d},\sigma;\mathbb R^{m})}\Bigr).
$$

**命题 3.3（Knothe-Rosenblatt 映射趋于恒等）。** 若 $\{\nu_i\}$ 绝对连续且二阶矩有限、$\mu$ 绝对连续、$\|\nu_i-\mu\|_{\mathrm{TV}}\to0$ 且 $\|\bm x\|_{L^2(\nu_i)}\to\|\bm x\|_{L^2(\mu)}$，则满足 $(T_i)_{\#}\mu=\nu_i$ 的 Knothe-Rosenblatt 映射 $T_i$ 在 $L^2(\mu)$ 中收敛到 $\mathrm{Id}$。**这一条是「耦合层在 $t=0$ 处取恒等」这一设计的理论依据**：命题 3.1 说基分布在小 $t$ 处全变差意义下接近目标，命题 3.3 于是说所需的输运接近恒等，所以流只用学一个小扰动。

**命题 4.1（$t\to0$ 时残差的发散——时间加权的依据）。** 在与命题 3.1 相同的假设下（$\bm g$ 为常数、$\bm G$ 非退化、$\|\nabla\bm f\|\le C_{\bm f}$），对 $\tilde p(\bm x,t\mid\bm x_0)=\mathcal N(\bm x;\bm m(t),\bm\Sigma(t))$、$\bm\Sigma(t)\sim t\bm G$ 且 $\tfrac{\mathrm d}{\mathrm dt}\bm m(t)$ 有界，

$$
\int_{\mathbb R^{d}}|r(\bm x,\bm x_{0},t;\tilde p)|^{2}\mathrm d\bm x=\mathcal O\bigl(t^{-(\frac d2+2)}\bigr),
\qquad
\mathbb E_{\bm x\sim\tilde p}\bigl[|r(\bm x,\bm x_{0},t;\tilde p)|^{2}\bigr]=\mathcal O\bigl(t^{-(d+2)}\bigr).
$$

要点在于两个速率对应**两种不同的采样测度**，这正是经验损失用两个不同指数的原因。**即便逼近任意好，残差在 $t\to0$ 处仍然发散**；加权后的残差才是 $\mathcal O(1)$。

### 数值实验

实现用 PyTorch，初始学习率 $0.001$，式 13 的 Gauss-Legendre 求积取 10 个节点。度量有三种：转移密度的相对 $L^2$ 误差、$\hat p_{\bm\theta}(\cdot,t)$ 的相对 $L^2$ 误差，以及一个混合带宽核的最大均值差异（核取 $K(\bm x,\bm y)=\tfrac13\sum$ 三个高斯，带宽为 $\hat\sigma^2/4$、$\hat\sigma^2$、$4\hat\sigma^2$，$\hat\sigma$ 是中位距离）。参考解在转移密度解析可得时用 Gauss-Kronrod 求积，否则用 **ADI 有限差分格式**（$\Delta t=0.001$、$\Delta x=0.01$）；最大均值差异的样本来自步长 $0.001$ 的 Euler-Maruyama。

**例 1：二维 Beneš 方程。** $\mathrm d\bm X_t=\tanh\bm X_t\,\mathrm dt+\bm I_d\mathrm d\bm W_t$，$\bm x_0\in[-1,1]^d$，转移密度有精确解

$$
p(\bm x,t\mid\bm x_{0})=\frac{1}{(2\pi t)^{d/2}}
\exp\Bigl(-\frac d2t\Bigr)\exp\Bigl(-\frac{1}{2t}\|\bm x-\bm x_{0}\|^{2}\Bigr)
\prod_{i=1}^{d}\frac{\cosh x_{i}}{\cosh x_{0,i}} .
$$

设定：8 个耦合层，每个条件网络用 32 个随机 Fourier 特征与两个 32 神经元隐层；10 次自适应迭代 × 300 轮；Adam，学习率 $0.001$ 每 2000 轮减半；$2\times10^5$ 个训练点、批量 $10^4$；验证用 $10^6$ 个样本。

三种损失的消融（都取 $\gamma_1=0.2$、$\gamma_2=0.6$、$\gamma_3=0.2$）：

| 损失      | 时间区间    | 结果                                                                         |
| --------- | ----------- | ---------------------------------------------------------------------------- |
| 朴素 PINN | $(0,1.5)$   | 误差曲线与**基分布本身重合**——$t\to0$ 附近的残差压倒一切，那里的训练毫无贡献 |
| 朴素 PINN | $(0.1,1.5)$ | 去掉 $t\to0$ 附近的 PDE 约束，$t>0.1$ 的精度改善，但原点附近只能靠流的连续性 |
| 时间加权  | $(0,1.5)$   | 全区间准确                                                                   |

从 $\bm x_0=(0,0)$ 出发在 $t=0.01,0.5,1.5$ 三个时刻取快照；$t=1.5$ 时流捕捉到一个**双峰**转移密度。误差随自适应迭代单调下降。

**例 2：Beneš 的初值迁移（二维）。** 初值取 $[-1,1]^d$ 上的乘积 Beta 密度族

$$
\textbf{Beta}(\bm x;\bm\zeta,\bm\eta)=\prod_{i=1}^{d}\frac{1}{\bm B(\zeta_{i},\eta_{i})}
\Bigl(\frac{1+x_{i}}{2}\Bigr)^{\zeta_{i}-1}\Bigl(\frac{1-x_{i}}{2}\Bigr)^{\eta_{i}-1},
$$

从均匀情形 $\bm\zeta=\bm\eta=(1,1)$ 开始；验证在 $[-5,5]^2$ 的 $10^4$ 个网格点上，$M\in\{10^4,10^5,10^6\}$。结论：相对误差在 $M\le10^5$ 之前遵循 Monte Carlo 的 $\mathcal O(M^{-1/2})$ 标度，**再往上则被转移密度本身的逼近误差主导**——这是对精度上限的一次干净识别。提议分布的对比：$q_1$ 在 $t\to0$ 时方差远低，但方差随 $t$ 增长；$q_2$（混合，$\alpha(t)=e^{-6t}$）在全区间平衡方差。作者还指出**均匀初始分布在 $\mathbb R^2$ 上是不连续的、会拖累传统求解器**，而本方法不受影响；最大均值差异这一度量正是为这类不连续 $p_0$ 引入的。

**例 3：四维 Beneš。** 三个 KRnet 阶段，$l_1=10$、$l_2=8$、$l_3=6$；每层网络同上；20 次自适应迭代 × 300 轮；$\gamma_1=0.2$、$\gamma_2=0.4$、$\gamma_3=0.4$；学习率 $0.001$ 每 4000 轮减半；$2\times10^5$ 个训练样本、批量 $2\times10^4$；验证用 $2\times10^6$ 个样本，域为 $[-1,1]^4\times[-5,5]^4$；评价用 $q_2$、$\alpha(t)=e^{-6t}$、$M=10^4$，测试点 $10^5$ 个落在 $[-5,5]^4$。**一处值得记录的工程细节：残差在被放大 1000 倍的解上计算，以避免数值下溢**——高维时 Fokker-Planck 解的数值会变得非常小。

**例 4：非线性漂移、常扩散（二维）。**

$$
\mathrm d(X_1,X_2)^{\top}=(2X_2,\ 2X_1-0.8X_2-0.2X_1^{3})^{\top}\mathrm dt
+\mathrm{diag}(\sqrt{0.4},\sqrt{0.8})\,\mathrm d\bm W_{t},
$$

架构与训练同二维 Beneš；评价用 $q_2$、$\alpha(t)=e^{-6t}$、$M=10^4$；验证在 $[-6,6]^2$ 的 $10^4$ 个网格点上，最大均值差异用 $5\times10^4$ 个样本，参考解来自 ADI。结论：小 $t$ 时高度准确（基分布接近真值），$t$ 增大后仍然可靠。

**例 5：非线性漂移、状态依赖扩散（二维）。**

$$
\mathrm d(X_1,X_2)^{\top}=(X_2,\ -0.5X_1-0.3X_1^{3}+\sin X_2)^{\top}\mathrm dt
+\mathrm{diag}(0.5+0.3X_1,\ 0.4+0.1\sin X_2)\,\mathrm d\bm W_{t},
$$

$\bm x_0\in[-1,1]^2$；4 次自适应迭代 × 1000 轮；$4\times10^5$ 个训练点、批量 $5\times10^4$；验证在 $[-5,5]^2$ 的 $10^4$ 个网格点上对 ADI，最大均值差异用 $5\times10^4$ 个样本。**注意理论（命题 3.1 与 4.1）假设 $\bm g$ 为常数，所以这个例子被刻意放在分析的适用范围之外**，作为鲁棒性检验。

这五组实验建立的是：线性化基分布加时间加权损失能在整个时间区间上给出准确的转移密度，而单用基分布或朴素 PINN 只在小 $t$ 处准、随 $t$ 增长而退化；并且**未观察到误差随时间增长**，作者归因于时间加权。它们没有建立的是状态依赖扩散下的理论保证，也没有与 DeepONet、Fourier 神经算子等算子学习基线做同题对比。

> [!warning] 定量结果的可核对性
> 论文正文没有数值表格，全部定量结果都在图 1 与图 5–11 中，因此**具体误差量级无法仅凭文本核对**。上文只记录已核对的设定、消融结构与定性结论。

> [!note] 题名与投稿状态
> 该文实际题名为 _A transition-density-based operator learning method for Fokker-Planck equations with various initial conditions_，与主页所列略有差别；主页题名接近其页眉，可能是早期工作题名。预印本采用 SIAM 格式并带「Submitted to the editor's DATE」脚注，但未标注期刊、DOI 或期刊引用信息。

### 与其他论文的关系

与编号 105 最近（两边共享 Xiaoliang Wan 与 Tao Zhou）：都是面向随机动力学的正规化流方法，都以**摊销**为目标——编号 107 靠学转移核摊销掉初始**分布**，编号 105 靠学以共享摘要为条件的前向与后向流摊销掉观测历史。两者都依赖 Markov 结构（这里是 Chapman-Kolmogorov，那里是后向递推）。

与编号 89 共享「让基分布承担难的那部分」这一设计：编号 89 的隐扩散有解析可得的 $p_D(z_t\mid z_0)$ 因而无需解随机微分方程，这里的基过程是解析可解的线性化随机微分方程因而流只学一个近恒等修正。两者也都用换元与对数行列式那套密度机器。

与编号 95 共享算子学习的取景：两者都学一个分布值映射，使新输入无需重训——编号 95 摊销的是对系数的传感器观测，编号 107 摊销的是初始律。与编号 103 及 PINN 文献共享**残差加权与因果性**的关切：编号 107 的 $t^{d+2}$ 加权是 PDE 残差的时间重加权，与因果训练权同属一类，但它是从一个显式的渐近速率（命题 4.1）导出的，而不是启发式的。随机 Fourier 特征这一成分与编号 101 共享，只是编号 101 拿它当多尺度基库，编号 107 拿它锐化对 $\bm x_0$ 的敏感性。KRnet、Fokker-Planck 的自适应采样，以及 $\gamma_1/\gamma_2/\gamma_3$ 混合测度都来自作者们自己的前作，见[[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|密度流一页]]与[[computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|自适应采样一页]]。

> [!note] 覆盖进度
> 四篇的构造、损失、算法与主要结果均已按预印本或期刊全文逐式核对。定量结果的可核对程度不一：编号 95 与 98 的表格数值为逐项转录；编号 75 的误差与运行时间位于图与补充材料 S7–S8，本次未逐项转录，因此本页不列其数字；编号 107 正文没有数值表格，全部量级都在图 1 与图 5–11 中，本页因此只记录设定、消融结构与定性结论。书目方面：编号 95 与 98 的 arXiv 记录上都没有 DOI 与期刊引用，所以其期刊状态未由 arXiv 核实；编号 107 只带 SIAM 格式的投稿脚注，无期刊信息。背景文献方面，Set Transformer 的 MAB/SAB/ISAB/PMA 公式未核对（只读到其补充材料），本页对置换不变性的论证因此只依赖 Deep Sets 定理与编号 95 自己的注意力池化公式。

## 四篇的对照

| 编号 | 不确定性来自       | 表示方式                  | 约束由谁保证         |
| ---- | ------------------ | ------------------------- | -------------------- |
| 75   | 输入落在分布外     | 置信门 + 高斯解码器       | 无（软约束）         |
| 95   | 观测稀疏或算子随机 | 隐变量 + 集合 transformer | 边界条件由架构硬约束 |
| 98   | 认知与偶然分开建模 | 高斯过程先验 + 神经算子   | 结构恒等式（硬约束） |
| 107  | 不涉及（初值可变） | 条件流 + 线性化基过程     | 密度约束由流保证     |

四篇的理论分量也差得很远，值得一并记下：

| 编号 | 定理内容                                                                                                                                                                   |
| ---- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 75   | 无定理；理论内容是变分下界推导与「压平分布使隐变量对分布外输入近似独立」的分析                                                                                             |
| 95   | 无定理；理论内容是证据下界推导与 $\sigma_u^2$ 乘 $M$ 的函数空间极限论证                                                                                                    |
| 98   | 无定理、无命题、无误差分析（论文自陈）                                                                                                                                     |
| 107  | 四条命题：基分布的短时全变差为 $\mathcal O(t)$（$\bm f$ 二阶导有界时其平方改进到 $\mathcal O(t^{3})$）、全变差下的期望传递、KR 映射趋于恒等、残差在 $t\to0$ 的两个发散速率 |

一条贯穿的判断：**「预测应当给出分布」这件事，在算子学习里可以落在三个不同的位置。** 编号 75 把它放在编码器的置信门上，编号 95 放在隐变量上，编号 98 放在先验的相关结构上。而编号 107 说明另一件事：有时正确的做法不是给预测加不确定性，而是换一个与参数无关的学习对象。

编号 98 还补上了一条前三者都没有明说的判断依据：**不确定性的表示形式应当由算子本身的结构决定。** $u$ 是 $f$ 的二重积分，所以 $u$ 的方差是 $f$ 的方差的全局累积，逐点独立的隐噪声在原理上就表达不了它——这不是精度问题，而是表示能力问题。

还有一条共同的方法论教训，来自编号 95 与编号 98 各自对基线的处理：**不确定性量化的结论只有在参考分布独立于被评价的模型时才可否证。** 编号 95 用精确的高斯过程后验加有限差分求解器造参考分布，因此「深度集成的带子太窄」是一个可检验的陈述；编号 98 缺少这样的参考，只能靠「均值误差大而标准差小」这种间接迹象论证过度自信，说服力相应弱一档。

## 本页原文

- L. Guo, H. Wu, W. Zhou, Y. Wang, and T. Zhou, [_IB-UQ: information bottleneck based uncertainty quantification for neural function regression and neural operator learning_](https://doi.org/10.1016/j.jcp.2024.113089), J. Comput. Phys. 510 (2024), 113089（预印本 [arXiv:2302.03271](https://arxiv.org/abs/2302.03271)）。
- L. Ma, L. Guo, H. Wu, and T. Zhou, [_Deep set based operator learning with uncertainty quantification_](https://doi.org/10.1016/j.jcp.2026.115011), J. Comput. Phys. (2026)（预印本 [arXiv:2509.25646](https://arxiv.org/abs/2509.25646)；arXiv v2 的署名顺序为 Lei Ma, Ling Guo, Hao Wu, Tao Zhou，arXiv 记录上无 DOI 与期刊引用，故期刊状态未由 arXiv 核实）。
- X. Feng, L. Guo, X. Wan, H. Wu, T. Zhou, and W. Zhou, _LVM-GP: uncertainty-aware PDE solver via coupling latent variable model and Gaussian process_, [arXiv:2507.22493](https://arxiv.org/abs/2507.22493)，投稿 J. Comput. Phys.（arXiv 记录上无 DOI 与期刊引用）。
- L. Zeng, X. Wan, Y. Wang, F. Nobile, and T. Zhou, _A transition-density-based operator learning method for Fokker-Planck equations with various initial conditions_, [arXiv:2606.09434](https://arxiv.org/abs/2606.09434)。
