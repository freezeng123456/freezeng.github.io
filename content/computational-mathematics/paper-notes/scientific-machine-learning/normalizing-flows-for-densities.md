---
title: 用可逆映射表示密度：从随机场到分数阶方程
description: 编号 62、64、72、87：把未知密度写成参考分布的推前，并用方程残差训练
lang: zh
translation: en/computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities
tags:
  - 论文笔记
  - 科学机器学习
  - 正规化流
---

![用可逆映射把密度方程写成参数化问题](assets/diagrams/tao-zhou-papers/zh/density-flow-solvers.svg)

这四篇的共同判断是：当方程的解**本身就是一个概率密度**时，用网格上的数值去逼近它是在浪费结构。密度必须非负且积分为一，这两条在网格方法里是约束，在可逆映射里是恒等式。

推前公式

$$
p_{X}(x)=p_{Z}\bigl(f(x)\bigr)\,\bigl|\det\nabla_x f(x)\bigr|
$$

自动给出非负性（右端是密度乘正的行列式）与归一化（变量替换）。剩下的问题只有两个：$f$ 怎么参数化才能既表达力强又让 $\det\nabla_x f$ 可算，以及在没有标签数据时用什么损失训练它。四篇论文分别在**目标对象**（随机场 / 时间依赖密度 / 分数阶方程 / 有界支撑）上推进。

## 62：把随机场写成带 Karhunen-Loève 结构的流

### 直觉

一个随机场是一族按空间位置编号的相关随机变量。要给它建模，必须同时管住两件事：**空间上的相关结构**与**单点边缘分布的非高斯性**。这两件事各有成熟工具，但各自都不够。

Karhunen-Loève 展开管相关结构：把场写成一组正交模态乘以互不相关的随机系数。它是高斯场的标准表示，但要用它得先知道协方差核，而在数据驱动的场景下协方差恰恰是未知的；而且截断阶数一旦变高，后续在参数空间建代理就撞上维数灾——论文点名的任意多项式混沌深度网络方法，其项数随有效维数指数增长。

正规化流管非高斯性：一个可逆网络把简单分布推成复杂分布。但流本身没有任何空间结构，它作用在一个固定维数的向量上，对「场」这个概念一无所知。

编号 62 的做法是把两者串起来：**先造一个带 Karhunen-Loève 展开结构的高斯参考场，但把展开系数交给网络去学；再用一个流把这个高斯场推成目标场。** 前半段提供空间相关性且不需要协方差先验，后半段提供非高斯性。

串起来之后还多得到一样东西，而这正是相对物理信息 GAN 的实质优势：**流是双射，因此似然可以显式写出。** 训练准则可以是老老实实的极大似然，而不是对抗训练。另一个容易被忽略的收益是，模型给出的是任意位置上场值的分布，而不是一个固定长度向量的分布，因此**传感器位置不必在各次观测之间保持固定**——这在实际测量中是常态而非例外。

### 问题设定

问题类为

$$
\mathcal N_x\bigl[u(x;\omega);k(x;\omega)\bigr]=f(x;\omega),
\quad x\in\mathcal D,\ \omega\in\Omega,
\qquad
\mathcal B_x\bigl[u(x;\omega)\bigr]=0,\quad x\in\Gamma,
$$

其中 $k$ 是随机输入场（如扩散系数），$f$ 是随机源项，$u$ 是随机解场。观测只有散点传感器读数，随机性的协方差结构未知。

### 推导

论文把模型的构造明确写成三步：

1. 构造一个带**截断 Karhunen-Loève 展开结构**的参考高斯随机场 $z(x,\omega)$，其中展开系数由深度网络参数化；
2. 构造参考场与目标随机场之间的双射变换，即一个正规化流；
3. 在散点观测上最大化对数似然之和，同时训练全部参数。

对随机微分方程与随机 PDE，已知物理由一个残差损失加入总损失，得到**物理信息**版本。也就是说，同一套表示既可以在纯数据模式下工作（只有第 3 步的似然），也可以在物理约束模式下工作（似然加残差），而正问题、反问题与混合问题的区别只在于哪些量被观测、哪些量被残差约束。这就是论文所说的「统一框架」的具体含义。

论文自陈的优势有三条：同一框架处理正问题、反问题与混合问题；不要求传感器位置在各次观测间固定；缓解多项式混沌的维数灾。

### 定理

本文不给收敛或逼近结果。这一点值得明说：它的贡献是表示与训练框架，不是分析。

### 数值实验

第 4 节分三组：

| 实验组 | 内容                                       |
| ------ | ------------------------------------------ |
| 一     | 学习随机过程，含非高斯场与**混合**非高斯场 |
| 二     | 正向随机椭圆方程                           |
| 三     | 反向随机椭圆方程                           |

「混合非高斯场」这一组是三组里最能说明问题的：它检验的正是「Karhunen-Loève 结构提供相关性、流提供非高斯性」这个分工是否真的成立——如果流不够表达，混合分布会被抹成单峰。

定性结论是模型能够学习非高斯过程，并求解若干类随机 PDE。

### 与其他论文的关系

这是本主题里最早的正规化流论文，为编号 64、72、87 以及[[computational-mathematics/paper-notes/scientific-machine-learning/spectral-bias-and-generative-solvers|编号 105]]、[[computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning|编号 107]] 立下模板：**用可逆网络表示一个未知的分布或场，并用物理残差而不是标签数据去训练它。** 它是编号 105（FLUID）的直接前身，后者把带 Karhunen-Loève 结构的参考场换成一个循环的摘要统计量加条件流。编号 64 在正文中把它作为基于流的不确定性量化代理引用。

## 64：时间不是一个额外维度

### 直觉

时间依赖 Fokker-Planck 方程的解不是一个密度，而是**一族按时间编号的密度**。要让流表示这一族，有三条看起来都合理的路。

第一条是每个时间片训一个流。这在计算上是灾难，而且时间方向上没有任何连续性保证——两个相邻时刻的解可能被学成完全不同的东西。

第二条是把时间当作第 $d+1$ 个维度，让流去学时空联合密度。这条路**在数学上就是错的**，理由很干脆：密度在每个固定时刻对空间积分为一，但对时空一起积分就等于时间区间长度而不是一，$\int p(x,t)\,\mathrm dx\,\mathrm dt\neq1$。流的推前公式保证的是「联合归一化」，而方程要的是「逐时刻归一化」，两者不兼容。

第三条是论文的选择：**条件化而非增广。** 一个流，但它的耦合层把 $t$ 当作额外输入；潜变量的时间坐标钉在真实时间上，$t^\ast=t$。这样一来 Jacobian 在时间方向退化成 $1$，剩下的只有空间块，逐时刻归一化自动成立，而一套权重覆盖整个时间区间。

还有一层收益与流本身无关：训练信号是方程残差，因此**既不需要标签数据，也不需要生成随机微分方程的样本轨道**。基于路径的方法要靠大量 Monte Carlo 轨道去估计密度，代价和噪声都在那里；残差路线把这一整块省掉了。

### 问题设定

底层随机微分方程与时间依赖 Fokker-Planck 方程为

$$
\mathrm d\bm X_t=\bm\mu(\bm X_t,t)\,\mathrm dt+\bm\sigma(\bm X_t,t)\,\mathrm d\bm W_t,
$$

$$
\frac{\partial p(\bm x,t)}{\partial t}
=-\sum_{i=1}^d\frac{\partial}{\partial x_i}\bigl[\mu_i(\bm x,t)p(\bm x,t)\bigr]
+\sum_{i=1}^d\sum_{j=1}^d\frac{\partial^2}{\partial x_i\partial x_j}
\bigl[D_{ij}(\bm x,t)p(\bm x,t)\bigr],
\qquad D=\tfrac12\sigma\sigma^{\mathsf T}.
$$

稳态版本去掉时间导数：

$$
\sum_i\frac{\partial}{\partial x_i}\bigl[\mu_i p\bigr]
+\sum_{i,j}\frac{\partial^2}{\partial x_i\partial x_j}\bigl[D_{ij}p\bigr]=0 .
$$

### 推导

**从增广到条件化。** 先按最朴素的想法写，令 $\widehat{\bm x}=(\bm x,t)$、$\widehat{\bm z}=(\bm z,t^\ast)$，

$$
p_{\widehat{\bm X}}(\widehat{\bm x})=p_{\widehat{\bm Z}}(\widehat{\bm z})\,|\det J|,
\qquad \widehat{\bm z}=f(\widehat{\bm x}).
$$

再把潜时间钉住，$t^\ast=t$。这一钉，Jacobian 的最后一行变成 $(0,\dots,0,1)$，行列式按块展开：

$$
\det J=\begin{vmatrix}\partial\bm z/\partial\bm x & \partial\bm z/\partial t\\ 0 & 1\end{vmatrix}
=\Bigl|\frac{\partial \bm z(\bm x,t)}{\partial\bm x}\Bigr| .
$$

注意右上块 $\partial\bm z/\partial t$ 一般不为零——流确实依赖时间——但它落在行列式的展开之外，不进入结果。于是

$$
p_{\widehat{\bm X}}(\bm x,t)=p_{\widehat{\bm Z}}(\bm z,t)\,
\Bigl|\frac{\partial\bm z}{\partial\bm x}\Bigr|,
\qquad \bm z=f(\bm x,t).
$$

这句话值得单列：**流是时间条件化的，不是时间增广的。** 一个字之差决定了归一化是否成立。

**复合与 Jacobian。** 流按层复合，Jacobian 按层连乘：

$$
\bm z=f(\bm x,t)=f_{[L]}\circ\cdots\circ f_{[1]}(\bm x,t),
\qquad
\bm x=f^{-1}(\bm z,t),
$$

$$
\bigl|\det\nabla_{\bm x}f(\cdot,t)\bigr|
=\prod_{i=1}^L\bigl|\det\nabla_{\bm x_{[i-1]}}f_{[i]}(\cdot,t)\bigr| .
$$

**架构。** 论文把自己的模型描述为「KRnet 从空间域到时空域的简化推广」。每个 $f_{[i]}$ 是一个 Actnorm 层接一个修改过的**时间依赖**仿射耦合层；最后一层用多项式样条变换提高表达力。

Actnorm（尺度加偏置）层及其逆：

$$
\bm y_{[i]}=\bm a_i\odot \bm x_{[i]}+\bm b_i,
\qquad
\bm x_{[i]}=\frac{\bm y_{[i]}-\bm b_i}{\bm a_i},
$$

其中 $\bm a_i,\bm b_i$ 由小批统计量作数据依赖初始化。

**时间依赖仿射耦合层**是与 real NVP 和普通 KRnet 对照时应当引用的对象：

$$
\bm x_{[i],1}=\bm x_{[i-1],1},
$$

$$
\bm x_{[i],2}=\bm x_{[i-1],2}\odot\Bigl(\bm 1_{d-m}
+\beta\tanh\bigl(\bm s_i(\bm x_{[i-1],1},t)\bigr)\Bigr)
+e^{\bm\zeta_i}\odot\tanh\bigl(\bm q_i(\bm x_{[i-1],1},t)\bigr),
$$

其中 $|\beta|<1$ 由用户指定（例如 $0.6$），$\bm s_i,\bm q_i:\mathbb R^{m+1}\to\mathbb R^{d-m}$，$\bm\zeta_i\in\mathbb R^{d-m}$ 可训练。

与标准 real NVP 的差别有两处，都有明确目的。其一，$\tanh$ 加上 $\beta$ 把尺度因子夹在 $(1-\beta,1+\beta)$ 内，而 real NVP 的尺度是 $\exp(\cdot)$，可以任意大或任意接近零；夹住之后 Jacobian 的对数行列式不会爆炸，训练更稳。其二，$\bm s_i,\bm q_i$ 的输入维数是 $m+1$ 而不是 $m$，多出来的那一维正是 $t$，这就是「一套权重覆盖整个时间区间」的实现方式。

**损失。** 训练信号是时间依赖 Fokker-Planck 残差，没有标签数据：在时空配点 $(\bm x^{(i)},t^{(i)})$ 上取均方残差，与 ADDA 一支的做法一致。

### 定理

本文不宣称定理。它的主张是方法论层面的：格式无网格、不需要标签数据、易于推广到高维。

### 数值实验

论文自己的算法特征列表可以概括为四步：用潜时间钉住真实时间的时间正规化流建模 $p(\bm x,t)$；由方程构造物理信息残差，不用样本轨道也不用标签；在时空配点上极小化均方残差；再从当前流自适应地刷新配点集并重复。

实验覆盖三类：

| 类别       | 内容                                      |
| ---------- | ----------------------------------------- |
| 线性漂移   | 时间依赖 Fokker-Planck 方程，漂移项为线性 |
| 非线性漂移 | 同上，漂移项含非线性                      |
| 高维       | 高维时间依赖问题                          |

线性漂移那一组的作用是校准：这类问题常有解析解或高精度参考解，因此可以判断误差究竟来自流的表达力还是来自残差的离散。非线性漂移与高维两组才是方法真正要打的目标。

### 与其他论文的关系

编号 64 是稳态 ADDA-KRnet 框架与后续流论文之间的时间依赖桥梁。编号 72 随后处理**分数阶**时间依赖 Fokker-Planck 方程，用的是非常相似的时间依赖耦合层（对比编号 72 的 MCTNF 层，它把 $\tanh$ 的自变量额外乘了一个 $t$）。编号 87 把无界 KRnet 换成有界版本，而[[computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning|编号 107]] 把「一个方程配一个流」换成一个**跨初值的算子**。时间正规化流这个想法本身，论文归功于 Both 与 Kusters。

## 72：分数阶算子的两种处理

### 直觉

编号 72 把这条路线推到分数阶 Fokker-Planck 方程，也就是同时受 Lévy 噪声与高斯噪声驱动的粒子密度。这里**三个困难同时出现**：解定义在无界区域上、维数可能很高、分数阶 Laplacian 是非局部的。网格方法在三条上同时失败。

流一次解决前两条：推前密度天然定义在整个 $\mathbb R^d$ 上且归一，维数只影响网络宽度而不引发网格爆炸。第三条需要单独处理，而论文给出两条路，它们的取舍是这篇文章最值得记住的部分。

**路线一是接受随机性。** 把分数阶 Laplacian 写成期望，用少量样本估计。代价与维数无关，但每次求值都带噪声。

**路线二是回避随机性。** 高斯函数的分数阶 Laplacian 有闭式（是一个合流超几何函数），因此换一个由高斯径向基混合构成的辅助模型，在它上面精确地算分数阶 Laplacian，再用一致性罚把辅助模型与流拉在一起。求值不带噪声，但径向基中心数随维数增长。

一条是维数友好但有噪声，另一条是无噪声但维数不友好，两者恰好互补。这种「把难算的算子换到一个算子已知的模型上」的手法在本主题的其他论文里没有对应物，是本文最独特的想法。

### 问题设定

带 Lévy 噪声的随机微分方程与相应的分数阶 Fokker-Planck 方程：

$$
\mathrm d X_t=\mu(X_t,t)\,\mathrm dt+\sigma(X_t,t)\,\mathrm dW_t+\mathrm dL^\alpha_t,
\qquad
\frac{\partial p}{\partial t}=\mathcal L p-(-\Delta)^{\alpha/2}p,
$$

$$
\mathcal L p=-\nabla\cdot(p\mu)+\tfrac12\nabla\cdot\nabla\cdot(\sigma\sigma^{\mathsf T}p),
\qquad
(-\Delta)^{\alpha/2}p=C_{d,\alpha}\,\mathrm{P.V.}\!
\int_{\mathbb R^d\setminus\{0\}}\frac{p(x)-p(y)}{|x-y|_2^{d+\alpha}}\,\mathrm dy .
$$

密度模型是流：$p_X(x)=p_Z(f(x))\,|\det\nabla_x f(x)|$，$f$ 取「简化 KRnet」，由 Actnorm 层

$$
y_{[i]}=a_i\odot x_{[i]}+b_i
$$

与仿射耦合层

$$
x_{[i],1}=x_{[i-1],1},
\qquad
x_{[i],2}=x_{[i-1],2}\odot\Bigl(\bm 1_{d-m}+\beta\tanh\bigl(s_i(x_{[i-1],1})\bigr)\Bigr)
+e^{\zeta_i}\odot\tanh\bigl(q_i(x_{[i-1],1})\bigr)
$$

复合而成。

> [!warning] 一处记号冲突
> 原文在这个耦合层公式里把 $\tanh$ 的尺度界印作 $\alpha$，而 $\alpha$ 在全篇是分数阶指数。应读作 $\beta$，与编号 64 的对应公式以及原始 KRnet 层一致。

### 推导

**路线一：Monte Carlo（MCNF）。** 损失是分数阶 Fokker-Planck 残差的均方：

$$
L(p_{\mathrm{KRnet},\theta})=\frac{1}{N_S}\sum_{i=1}^{N_S}\bigl|R_\theta(x^i)\bigr|^2,
\qquad
R_\theta(x)=\bigl(\mathcal L-(-\Delta)^{\alpha/2}\bigr)p_{\mathrm{KRnet},\theta}(x).
$$

其中的分数阶 Laplacian 由[[computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|编号 66]]的内外分裂估计给出，在本文中作为 Lemma 3.1：

$$
(-\Delta)^{\alpha/2}u(\bm x)=C_{d,\alpha}\frac{|S^{d-1}|r_0^{2-\alpha}}{2(2-\alpha)}\,
\mathbb E_{\bm\xi\sim U(S^{d-1}),\,r_1\sim f_{\mathrm I}}\!
\left[\frac{2u(\bm x)-u(\bm x-r_1\bm\xi)-u(\bm x+r_1\bm\xi)}{r_1^2}\right]
$$

$$
\qquad+\;C_{d,\alpha}\frac{|S^{d-1}|r_0^{-\alpha}}{2\alpha}\,
\mathbb E_{\bm\eta\sim U(S^{d-1}),\,r_2\sim f_{\mathrm O}}
\bigl[2u(\bm x)-u(\bm x-r_2\bm\eta)-u(\bm x+r_2\bm\eta)\bigr],
$$

$$
f_{\mathrm I}(r)=\frac{2-\alpha}{r_0^{2-\alpha}}r^{1-\alpha}\mathbf 1_{r\in[0,r_0]},
\qquad
f_{\mathrm O}(r)=\alpha r_0^\alpha r^{-1-\alpha}\mathbf 1_{r\in[r_0,\infty)},
$$

$$
r_1/r_0\sim\mathrm{Beta}(2-\alpha,1),
\qquad
r_0/r_2\sim\mathrm{Beta}(\alpha,1),
\qquad
r_\epsilon=\max\{\epsilon,r_1\}.
$$

施加的对象与编号 66 不同：那里是通用 PINN 代理，这里是一个**密度**，因此非负性与归一化是由流保证的，而不需要额外约束。公式本身也有一处小差别：编号 72 在内层与外层分别用两个**独立**的方向变量 $\bm\xi$ 与 $\bm\eta$，编号 66 的对应公式复用同一个 $\xi$。

一个容易忽略但实际重要的观察：$\mathrm{Beta}(a,1)$ 在 $a\to0$ 时向原点集中，而内层半径服从 $r_1/r_0\sim\mathrm{Beta}(2-\alpha,1)$，因此 $\alpha$ 越接近 $2$，第一个参数 $2-\alpha$ 越接近零，样本 $r_1$ 越堆积在零附近，$1/r_1^2$ 越容易放大相消误差，需要**更大**的下界保护 $r_\epsilon$ 才稳定。这把「参数取值影响数值稳定性」写成了可操作的规则，而不是一句泛泛的告诫。

**路线二：解析辅助模型（GRBFNF）。** 出发点是 Lemma 3.2：对高斯 $u(\bm x)=\exp(-\sigma^{-2}|\bm x-\bm x_0|_2^2)$，

$$
(-\Delta)^{\alpha/2}u(\bm x)=c_{\alpha,d}\,|\sigma|^{-\alpha}\;
{}_1F_1\!\Bigl(\frac{d+\alpha}{2};\frac{d}{2};-\sigma^{-2}|\bm x-\bm x_0|^2_2\Bigr),
\qquad
c_{\alpha,d}=\frac{2^\alpha\Gamma\bigl(\frac{d+\alpha}{2}\bigr)}{\Gamma\bigl(\frac d2\bigr)},
$$

${}_1F_1$ 是合流超几何函数。于是引入一个高斯径向基混合作为辅助模型，中心集为 $S_{\text{center}}=\{\tilde{\bm x}_i\}_{i=1}^M$：

$$
p_{\mathrm{GRBF},\tilde\theta}(\bm x)=\sum_{i=1}^M w_i\,\mathcal N(\tilde{\bm x}_i,\sigma_i^2\mathbf I)(\bm x),
\qquad 0\le w_i\le1,\ \sum_i w_i=1,
$$

$$
\mathcal N(\tilde{\bm x}_i,\sigma_i^2\mathbf I)(\bm x)
=(2\pi)^{-d/2}\sigma_i^{-d}\exp\!\Bigl(-\frac{|\bm x-\tilde{\bm x}_i|_2^2}{2\sigma_i^2}\Bigr),
$$

其中 $w_i$ 与 $\sigma_i$ 都可训练。由 Lemma 3.2 逐项作用并整理常数，它的分数阶 Laplacian 有闭式：

$$
(-\Delta)^{\alpha/2}p_{\mathrm{GRBF},\tilde\theta}(\bm x)
=c_{\alpha,d}\,\pi^{-d/2}2^{-\frac{d+\alpha}{2}}\sum_{i=1}^M w_i\,|\sigma_i|^{-(d+\alpha)}\,
{}_1F_1\!\Bigl(\frac{d+\alpha}{2};\frac d2;-\frac{|\bm x-\tilde{\bm x}_i|_2^2}{2\sigma_i^2}\Bigr).
$$

耦合损失把流与辅助模型绑在一起。第一项把**流的局部部分**与**辅助模型的非局部部分**放在同一个残差里，第二项是一致性罚：

$$
\tilde L\bigl(p_{\mathrm{KRnet},\theta},p_{\mathrm{GRBF},\tilde\theta}\bigr)
=\frac{1}{N_S}\sum_{i=1}^{N_S}\Bigl(\mathcal L p_{\mathrm{KRnet},\theta}(x^i)
-(-\Delta)^{\alpha/2}p_{\mathrm{GRBF},\tilde\theta}(x^i)\Bigr)^2
+\frac{\beta_m}{N_S}\sum_{i=1}^{N_S}\Bigl(p_{\mathrm{KRnet},\theta}(x^i)
-p_{\mathrm{GRBF},\tilde\theta}(x^i)\Bigr)^2 .
$$

第一项的构造值得停一下：如果没有第二项，两个模型可以各自漂移到无关的函数上而第一项仍然很小；一致性罚 $\beta_m$ 正是把「辅助模型必须代表同一个密度」这句话变成一个可优化的量。这是全篇最独特的想法：**把难算的算子换到一个算子已知的模型上，再用一致性罚把两个表示拉在一起。**

**时间依赖版本（MCTNF）中的一个细节。** 耦合层是

$$
x_{[i],2}=x_{[i-1],2}\odot\Bigl(\bm 1_{d-m}
+\beta\tanh\bigl(t\,s_{i,t}(x_{[i-1],1},t)\bigr)\Bigr)
+e^{\zeta_i}\odot\tanh\bigl(t\,q_{i,t}(x_{[i-1],1},t)\bigr).
$$

注意网络输出前显式乘了一个 $t$。$\tanh(0)=0$，因此 $t=0$ 时尺度因子变成 $\bm 1$、平移项变成 $\bm 0$，整层退化为恒等映射，流在初始时刻就是参考分布本身。这意味着**初值条件被构造性地满足**，而不是通过一个罚项去逼近。这类「把约束写进架构」的做法与推前公式保证归一化是同一思路。

### 定理

本文的形式结果是两条引理，都不是关于网络的收敛性陈述，而是关于算子求值的恒等式：

**Lemma 3.1（分数阶 Laplacian 的 Monte Carlo 表示）.** 见上，条件是分裂半径 $r_0>0$、$0<\alpha<2$，以及两个 Beta 采样律；$r_\epsilon$ 只是数值保护，不进入恒等式本身。

**Lemma 3.2（高斯的解析分数阶 Laplacian）.** 见上，对形如 $\exp(-\sigma^{-2}|\bm x-\bm x_0|^2)$ 的函数成立。论文注明这条引自其参考文献而非自证。

### 数值实验

自适应策略本身很短（对 MCNF 陈述，GRBFNF 类似）：

1. 从一个简单分布（均匀或高斯）抽初始训练集；
2. 在当前训练集上极小化残差损失训练流（MCNF 用 Monte Carlo 分数阶 Laplacian；GRBFNF 用辅助模型的解析分数阶 Laplacian 加一致性项）；
3. 把参考样本推过逆流生成**新的**训练集，也就是从当前近似密度采样；
4. 在刷新后的集合上继续训练，交替第 2、3 步直到停机。

这就是「交替细化训练集与近似解」的具体含义。它与[[computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|编号 70]]那一支的区别在于所依据的量：那里是残差失效概率，这里是**解密度本身**。

算例覆盖四类：

| 算例                                     | 检验的是什么                                                |
| ---------------------------------------- | ----------------------------------------------------------- |
| 二维分数阶 FP，仅由分数阶 Laplacian 驱动 | 非局部项的处理，隔离掉漂移与扩散的影响                      |
| 双峰目标分布                             | 流的表达力与自适应是否会塌到单峰                            |
| 更高维的稳态 FP 方程                     | 两条路线在维数上的分野                                      |
| 时间依赖分数阶 FP，解为 Cauchy 分布      | 时间条件化与重尾解；Cauchy 是 $\alpha$ 稳定噪声的自然稳态律 |

最后一例的选择很讲究：Cauchy 分布是重尾的，二阶矩不存在，正是有限区域网格方法最难处理的对象，而流在整个 $\mathbb R^d$ 上定义，重尾只是参考分布经过映射之后的形状问题。

### 与其他论文的关系

编号 72 等于[[computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|编号 66]]的 Monte Carlo 分数阶算子乘以编号 64 与 Tang-Wan-Liao 的 ADDA/KRnet 密度框架。它的耦合层与编号 64 同族，自适应刷新与 ADDA 同理。编号 87 随后把无界 KRnet 换成适用于有界与混合区域的 B-KRnet，而[[computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning|编号 107]] 把「解一个 FP 方程」的目标换成「学一个跨初值的 FP 解算子」。

## 87：有界支撑需要另一种耦合层

### 直觉

KRnet 与一般正规化流映到 $\mathbb R^d$ 上的高斯参考，诱导密度的支撑因此是整个空间。对真正支撑在有界集上的密度这是错的——论文举的日常例子是人的年龄，模型不该给 $-3$ 岁和 $700$ 岁分配正概率——对定义在超矩形上的 PDE 解同样是错的。

直接的补救是把流硬性限制在盒子里，但这会破坏精确可逆性：一个映到盒外的点没有像，逆映射就不再处处有定义，推前公式随之失效。

论文的思路是**从要求出发反推构件**。需要的是一个 $[-1,1]\to[-1,1]$ 的单调双射，而且它的逆要能显式算。哪一类函数天然满足？**累积分布函数。** 任何在 $[-1,1]$ 上严格正的密度，其累积分布函数把 $[-1,1]$ 单调地映满 $[0,1]$，再作仿射拉伸就落回 $[-1,1]$。剩下的只是选一个既有表达力又能显式求逆的密度族。

论文选的是**分段线性**密度。它的累积分布函数是分段二次的，求逆就是解一元二次方程，有根式解。这个选择把「可逆」从一个需要小心维护的性质变成了一行代数。

代价也很明确，而且论文自己写在结论里：模型密度是分段线性的，因此**只有一阶可微**。二阶 PDE 无法直接作用在它上面，必须先改写成一阶系统。

### 问题设定

流密度与其复合结构：

$$
p_{\bm X}(\bm x)=p_{\bm Z}\bigl(f(\bm x)\bigr)\,\bigl|\det\nabla_{\bm x}f(\bm x)\bigr|,
\qquad
f=f_{[L]}\circ\cdots\circ f_{[1]},
$$

$$
f^{-1}=f^{-1}_{[1]}\circ\cdots\circ f^{-1}_{[L]},
\qquad
\bigl|\det\nabla_{\bm x}f\bigr|=\prod_i\bigl|\det\nabla_{\bm x_{[i-1]}}f_{[i]}\bigr| .
$$

关键差别在参考分布：$\bm Z$ 取 $[-1,1]^d$ 上的**均匀分布**，而不是 $\mathbb R^d$ 上的高斯。

结构上的依据是 Knothe-Rosenblatt 重排（论文的 Prop. 2.2）：存在三角映射 $T$ 使 $T_\#\rho=\pi$，

$$
\bm z=T(\bm x)=\bigl(T_1(x_1),\,T_2(x_1,x_2),\dots,T_d(x_1,\dots,x_d)\bigr)^{\mathsf T},
$$

它可以分解成 $d$ 个每次只更新一个坐标的映射的复合。

### 推导

**伪三角（块）结构。** 把 $\bm x$ 分块为 $\bm x=\bigl((\bm x^{(1)})^{\mathsf T},\dots,(\bm x^{(K)})^{\mathsf T}\bigr)^{\mathsf T}$，$\sum_i d_i=d$，把逐坐标的三角结构放松成逐块的：

$$
\bm z=f_{\text{KR}}(\bm x)=
\begin{pmatrix}\tilde f_1\\ \bm x^{(2:K)}\end{pmatrix}\circ
\begin{pmatrix}\tilde f_2\\ \bm x^{(3:K)}\end{pmatrix}\circ\cdots\circ
\begin{pmatrix}\tilde f_{K-1}\\ \bm x^{(K)}\end{pmatrix}\circ \tilde f_K(\bm x),
$$

其中 $\tilde f_k:[-1,1]^{\sum_{i\le k}d_i}\to[-1,1]^{\sum_{i\le k}d_i}$。外层循环有 $K$ 个阶段，$f^{outer}_{[k]}=\tilde f_{K-k+1}$，每完成一个阶段就冻结一个块（论文称之为使若干维度失活的压缩操作）；内层再由若干 CDF 耦合层复合，$f^{outer}_{[k]}=L_{\text{CDF},[k,l_k]}\circ\cdots\circ L_{\text{CDF},[k,1]}$。

**CDF 耦合层。** 先在网格 $-1=s_0<s_1<\cdots<s_n=1$ 上定义一个分段线性密度

$$
p(s)=\frac{w_{i+1}-w_i}{h_i}(s-s_i)+w_i,
\quad s\in[s_i,s_{i+1}],
\quad p(s_i)=w_i\ge0,\ h_i=s_{i+1}-s_i .
$$

逐段积分（一次项给出二次项，常数项给出一次项，前面各段的贡献按梯形面积累加）得到**分段二次**的累积分布函数：

$$
F(s)=\frac{w_{i+1}-w_i}{2h_i}(s-s_i)^2+w_i(s-s_i)
+\sum_{k=0}^{i-1}\frac{w_k+w_{k+1}}{2}(s_{k+1}-s_k),
\quad s\in[s_i,s_{i+1}].
$$

求逆就是在每段上解一个一元二次方程。记 $q_0=0$、$q_i=\sum_{k=0}^{i-1}\frac{w_k+w_{k+1}}{2}(s_{k+1}-s_k)$，则

$$
F^{-1}(q)=s_i+\frac{-w_i+\sqrt{w_i^2+2(w_{i+1}-w_i)(q-q_i)/h_i}}{(w_{i+1}-w_i)/h_i}
=s_i+\frac{2(q-q_i)}{w_i+\sqrt{w_i^2+2(w_{i+1}-w_i)(q-q_i)/h_i}} .
$$

两种写法在数学上等价，第二种由分子有理化得到。它在 $w_{i+1}\approx w_i$（该段几乎是常密度，第一种写法的分母趋于零）时数值稳定，这是实现层面必要的细节而非风格偏好。

令 $\tilde F(s)=2F(s)-1$ 得到 $[-1,1]\to[-1,1]$ 的双射，逐分量作用为 $F_i(y_i)=\tilde F(y_i;\bm\theta_i)$。耦合层对 $\bm y=(\bm y_1^{\mathsf T},\bm y_2^{\mathsf T})^{\mathsf T}\in[-1,1]^l$、$\bm y_1\in\mathbb R^m$ 定义为

$$
\begin{pmatrix}\hat{\bm y}_1\\ \hat{\bm y}_2\end{pmatrix}
=L_{\text{CDF},\bm\theta}\begin{pmatrix}\bm y_1\\ \bm y_2\end{pmatrix}:
\qquad
\hat{\bm y}_1=\bm y_1,
\qquad
\hat{\bm y}_2=\bm F\bigl(\bm y_2;\bm\theta(\bm y_1)\bigr),
$$

参数 $\bm\theta=(\bm s_1^{\mathsf T},\dots,\bm s_{n-1}^{\mathsf T},\bm w_0^{\mathsf T},\dots,\bm w_n^{\mathsf T})^{\mathsf T}\in\mathbb R^{2n(l-m)}$ 由一个网络 $\mathrm{NN}(\bm y_1)$ 产生，并重参数化以保证网格严格递增与权非负：

$$
\bm s_1=-1+\frac2n\bigl(1+\beta_1\tanh(\hat{\bm s}_1)\bigr),
\qquad
\bm s_{i+1}=\bm s_i+\Bigl(\frac{1-\bm s_i}{n-i}\Bigr)\bigl(1+\beta_{i+1}\tanh(\hat{\bm s}_{i+1})\bigr),
\qquad
\bm w_i=\frac{1+\gamma_i\tanh(\hat{\bm w}_i)}{C},
$$

$C$ 是归一化常数，实现中取 $\beta_1=\tfrac{65}{66}$、$\beta_i=0.97$（$i\ge2$）、$\gamma_i=0.99$。这些常数都严格小于 $1$，作用是把 $1+\beta\tanh(\cdot)$ 夹在 $(1-\beta,1+\beta)\subset(0,2)$ 内：网格增量恒为正因而 $s_i$ 严格递增，权 $w_i$ 恒为正因而密度不退化。之后交换 $\bm y_1,\bm y_2$ 的角色再作用一层，按层号奇偶交替，使每个坐标都被更新到。

自由度按层累加，$\mathrm{DOFs}(p_{\text{B-KRnet},\bm\theta})=\sum_{k=1}^{K}\sum_{i=1}^{l_k}\mathrm{DOFs}(L_{\text{CDF},[k,i]})$，其中当 $d_1=1$ 时 $\mathrm{DOFs}(f^{outer}_{[K]})=0$。

**一个边界情形（Remark 2.3）。** CDF 耦合层需要输入维数至少为 $2$（否则没有 $\bm y_1$ 去产生参数），因此当 $\bm x^{(1)}=x_1\in\mathbb R$ 时，$f_{\text{KR}}$ 只由 $K-1$ 个双射组成，最外层 $f^{outer}_{[K]}$ 取恒等映射。论文给的理由是此时 $\bm x^{(1)}_{[K-1]}$ 已经接近均匀，$2F(\cdot)-1$ 也接近恒等，因此损失不大。

### 与 KRnet 的四点差别

| 项目       | KRnet                  | B-KRnet                      |
| ---------- | ---------------------- | ---------------------------- |
| 参考分布   | $\mathbb R^d$ 上的高斯 | $[-1,1]^d$ 上的均匀分布      |
| 耦合层     | 仿射（尺度加平移）     | 单调分段二次累积分布函数     |
| 尺度偏置层 | 有                     | **去掉**（区域有界，不需要） |
| 光滑性     | 由激活决定             | 密度分段线性，只有一阶可微   |

第二条是本质：仿射映射把 $[-1,1]^m$ 映到别处，CDF 映射保持它不变。第三条是第二条的推论——尺度偏置层的用处是把数据规整到合适的量程，而在有界域上这件事由域本身完成。第四条有直接后果：**二阶 PDE 必须改写成一阶系统。**

一般乘积区域 $\prod_i[a_i,b_i]$ 由线性映射 $\bm y=\hat{\bm a}\odot\bm x+\hat{\bm b}$ 归一到 $[-1,1]^d$，其中 $\hat a_i=2/(b_i-a_i)$、$\hat b_i=-(b_i+a_i)/(b_i-a_i)$。

### 定理

**Prop. 2.4（推前的弱收敛）.** 若 $f_n\to T$ 在 $L^p([-1,1]^d)$ 中成立且 $T_\#\rho=\pi$，则 $(f_n)_\#\rho\rightharpoonup T_\#\rho$ 且 $(f_n^{-1})_\#\pi\rightharpoonup(T^{-1})_\#\pi$。论文给出证明梗概。

**Prop. 2.6（有界域上的普适性）.** 对有界 $\Omega$ 上任何绝对连续的 $\rho$，存在一列 B-KRnet 变换 $f_k$ 使 $(f_k)_\#\rho$ 弱收敛到 $\Omega$ 上的均匀测度。论文陈述了这条结论，但把证明留给另一篇文章。

论文另外指出，在 $l_k=1$ 且 $n\to\infty$ 的极限下，CDF 耦合层就是 Knothe-Rosenblatt 重排的直接逼近——这把架构与经典构造之间的关系说清楚了：不是启发式类比，而是有极限意义的逼近。

### 数值实验

**用 B-KRnet 解 PDE。** 对二阶算子 $\mathcal N[\bm x;\{\partial^{\bm\alpha}_{\bm x}p,|\bm\alpha|\le2\}]=0$ 配 $\int_\Omega p=1$、$p\ge0$（若质量是 $M$ 而非 $1$，把 $M$ 吸收进一个可训练标量 $\zeta$，使 $\zeta\hat p\approx p$），引入 $\bm g(\bm x)=\nabla p$ 由单独网络 $\bm g_{\text{NN}}$ 逼近，降为一阶系统，极小化

$$
\mathcal L\bigl(p_{\text{B-KRnet},\bm\theta},\bm g_{\text{NN}}\bigr)
=\lambda_{pde}\mathcal L_{pde}+\lambda_b\mathcal L_b+\lambda_{\bm g}\mathcal L_{\bm g},
$$

$$
\mathcal L_{pde}=\mathbb E_{\bm x\sim\rho}\bigl[\bigl|\mathcal N[\bm x;p_{\text{B-KRnet},\bm\theta},\bm g_{\text{NN}}]\bigr|^2\bigr],
\qquad
\mathcal L_b=\mathbb E_{\bm x\sim\rho_b}\bigl[|\mathcal B[\cdots]|^2\bigr],
$$

$$
\mathcal L_{\bm g}=\mathbb E_{\bm x\sim\rho}
\bigl[\|\bm g_{\text{NN}}(\bm x)-\nabla p_{\text{B-KRnet},\bm\theta}(\bm x)\|_2^2\bigr].
$$

第三项是把 $\bm g_{\text{NN}}$ 钉在真梯度上的一致性罚，与编号 72 里把 GRBF 钉在流上的那一项是同一类构造。非负性与质量守恒**由构造成立**，不进入罚项——这正是用流表示密度的全部理由。

**混合区域。** 动理学 Fokker-Planck 方程的解是位置 $\bm x\in\Omega$（有界）与速度 $\bm v\in\mathbb R^d$（无界）的联合密度。论文按 $p(\bm x,\bm v)=h(\bm v|\bm x)q(\bm x)$ 分解并建模

$$
p_{\bm\theta}(\bm x,\bm v)=h_{\text{KRnet},\bm\theta_1}(\bm v|\bm x)\cdot q_{\text{B-KRnet},\bm\theta_2}(\bm x),
$$

即无界 KRnet 作为速度的条件流乘以 B-KRnet 作为位置的密度。这里的条件化机制正是编号 64 的时间条件化，只是把条件变量从时间换成位置——同一个技术零件在不同论文里承担不同角色。

**自适应策略（第 4.2 节四步加 Algorithm 1）：**

1. 取 $\rho(\bm x)=1/|\Omega|$，抽 $C^0_{pde}=\{\bm x^{i,0}\}_{i=1}^{N_{pde}}\sim\mathrm{Unif}(\Omega)$ 与 $C_b=\{\bm x^i_b\}_{i=1}^{N_b}\sim\mathrm{Unif}(\partial\Omega)$；
2. 在 $C^0_{pde}$、$C_b$ 上训练 B-KRnet 得 $\bm\theta^{*,0}$；
3. 更新 $\rho(\bm x)\leftarrow(1-\gamma)\rho(\bm x)+\gamma\,p_{\text{B-KRnet},\bm\theta^{*,0}}(\bm x)$，取 $N_{new}=\gamma N_{pde}$，新点由逆流生成，$\bm z^{i,1}\sim\mathrm{Unif}(-1,1)^d$、$\bm x^{i,1}=(f_{\bm\theta^{*,0}})^{-1}(\bm z^{i,1})$，其余点保留；
4. 重复第 2、3 步共 $N_{\mathrm{adaptive}}$ 次。

Algorithm 1 在此之上加了小批量、Adam 与阶梯衰减学习率（每 $n_s$ 步 $l_r\leftarrow\eta\,l_r$）。第 3 步的混合形式值得注意：新的采样密度不是直接换成当前流，而是与旧密度按 $\gamma$ 加权混合——这与[[computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|编号 80]]给密度加下界、编号 73 用余弦退火从均匀过渡到自适应，属于同一类安全阀。

**实验清单。** 密度估计（第 3 节）与 PDE 逼近（第 6 节）各三例：

| 节  | 问题                                          | 检验的是什么                     |
| --- | --------------------------------------------- | -------------------------------- |
| 3.1 | 环形区域上的密度估计                          | 非凸支撑                         |
| 3.2 | 高斯混合                                      | 多峰表达力                       |
| 3.3 | 带孔洞的 logistic 分布                        | 支撑不连通                       |
| 6.1 | 四维 $-\Delta p+p=f$                          | 一阶系统改写在中等维数下是否可用 |
| 6.2 | 稳态 Keller-Segel 系统（两个耦合密度）        | 多个密度联立                     |
| 6.3 | 二维稳态动理学 Fokker-Planck，KRnet ⊗ B-KRnet | 有界与无界维度混合               |

共同设定为：`Tanh` 激活；**每个 CDF 耦合层三个子区间**；$\mathrm{NN}(\bm y_1)$ 是两个隐层的全连接网络；PyTorch 中的 Adam；**更新率 $\gamma=0.8$**。

定性结论是：B-KRnet 生成的样本与真值高度吻合，自适应密度逼近格式在四维问题、Keller-Segel 方程与动理学 Fokker-Planck 方程上都有效。论文自陈的限制有两条：B-KRnet 的光滑性只到一阶，要提高需要更高次多项式；在深度 Ritz 框架内做自适应点更新留作开放问题——后者由[[computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|编号 80]]接手。

### 与其他论文的关系

B-KRnet 是贯穿本主题的 KRnet 家族里的有界成员。它是[[computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|编号 80]]（深度 Ritz 的自适应重要采样）所用的密度模型。它处理混合区域的条件流技巧来自编号 64，只是把条件变量从时间换成位置。它是编号 72 与 Tang-Wan-Liao 原始 ADDA-FPE 工作中所用 KRnet 的有界对应物，而[[computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning|编号 107]] 把条件流用于 Fokker-Planck 的想法推广成一整个跨初值的算子。

## 四篇的推进关系

| 编号 | 目标对象               | 参考分布           | 训练信号                   | 新增技术零件               |
| ---- | ---------------------- | ------------------ | -------------------------- | -------------------------- |
| 62   | 散点观测下的随机场     | 带 KL 结构的高斯场 | 对数似然（可加物理残差）   | 系数由网络参数化的 KL 展开 |
| 64   | 时间依赖 Fokker-Planck | 无界高斯           | 方程残差                   | 潜时间钉住真实时间         |
| 72   | 分数阶 Fokker-Planck   | 无界高斯           | 方程残差                   | Monte Carlo 或解析辅助模型 |
| 87   | 有界支撑的密度         | $[-1,1]^d$ 均匀    | 方程残差（改写成一阶系统） | 分段二次 CDF 耦合层        |

一条贯穿的判断：**结构性约束应当由架构保证，而不是由罚项逼近。** 归一化由变量替换保证，非负性由推前公式保证，初值条件由耦合层乘 $t$ 保证，有界支撑由 CDF 耦合层保证。每加一条这样的保证，损失函数里就少一个需要调权重的罚项。

第二条判断关于代价：每一条构造性保证都要用表达力或光滑性去换。乘 $t$ 让该层在 $t=0$ 附近的变化受限；CDF 耦合层换来精确可逆，但把模型密度压成分段线性，于是二阶方程必须降阶。**没有免费的结构。**

## 本页原文

- L. Guo, H. Wu, and T. Zhou, _Normalizing field flows: solving forward and inverse stochastic differential equations using physics-informed flow models_, J. Comput. Phys. 461 (2022), 111202（预印本 [arXiv:2108.12956](https://arxiv.org/abs/2108.12956)）。
- X. Feng, L. Zeng, and T. Zhou, [_Solving time dependent Fokker-Planck equations via temporal normalizing flow_](https://doi.org/10.4208/cicp.OA-2022-0090), Commun. Comput. Phys. 32(2) (2022), pp. 401-423（预印本 [arXiv:2112.14012](https://arxiv.org/abs/2112.14012)）。
- L. Zeng, X. Wan, and T. Zhou, [_Adaptive deep density approximation for fractional Fokker-Planck equations_](https://doi.org/10.1007/s10915-023-02379-z), J. Sci. Comput. 97 (2023), 68（预印本 [arXiv:2210.14402](https://arxiv.org/abs/2210.14402)）。
- L. Zeng, X. Wan, and T. Zhou, _Bounded KRnet and its applications to density estimation and approximation_, SIAM J. Sci. Comput. 47(6) (2025), pp. C1294-C1318（预印本 [arXiv:2305.09063](https://arxiv.org/abs/2305.09063)）。
