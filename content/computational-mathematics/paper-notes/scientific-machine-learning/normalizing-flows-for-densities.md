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

> [!note] 本页覆盖
> 编号 **62**（_J. Comput. Phys._ 461, 2022）、**64**（_Commun. Comput. Phys._ 32(2), 2022）、**72**（_J. Sci. Comput._ 97:68, 2023）、**87**（_SIAM J. Sci. Comput._ 47(6), 2025）。

![用可逆映射把密度方程写成参数化问题](assets/diagrams/tao-zhou-papers/zh/density-flow-solvers.svg)

这四篇的共同判断是：当方程的解**本身就是一个概率密度**时，用网格上的数值去逼近它是在浪费结构。密度必须非负且积分为一，这两条在网格方法里是约束，在可逆映射里是恒等式。

推前公式

$$
p_{X}(x)=p_{Z}\bigl(f(x)\bigr)\,\bigl|\det\nabla_x f(x)\bigr|
$$

自动给出非负性（右端是密度乘正的行列式）与归一化（变量替换）。剩下的问题只有两个：$f$ 怎么参数化才能既表达力强又让 $\det\nabla_x f$ 可算，以及在没有标签数据时用什么损失训练它。四篇论文分别在**目标对象**（随机场 / 时间依赖密度 / 分数阶方程 / 有界支撑）上推进。

## 62：把随机场写成带 Karhunen-Loève 结构的流

### 问题

数据驱动的正反随机 PDE 问题中，随机输入场（扩散系数、源项）或解只能通过散点传感器观测，而随机性的协方差结构事先并不知道。经典路线先用 Karhunen-Loève 展开降维、再在参数空间建代理，这需要协方差的先验知识，且受维数灾困扰；论文引用的任意多项式混沌深度网络方法的项数随有效维数指数增长，而物理信息 GAN 不给出可计算的似然。

### 三步构造

论文把模型的构造明确写成三步：

1. 构造一个带**截断 Karhunen-Loève 展开结构**的参考高斯随机场 $z(x,\omega)$，其中展开系数由深度网络参数化；
2. 构造参考场与目标随机场之间的双射变换，即一个正规化流；
3. 在散点观测上最大化对数似然之和，同时训练全部参数。

对随机微分方程与随机 PDE，已知物理由一个残差损失加入总损失，得到物理信息版本。

这一构造的价值在于它同时解决两件事：Karhunen-Loève 结构给出场的空间相关性表示，而流给出非高斯性；由于流是双射，似然可以显式写出，因此训练准则是极大似然而不是对抗。论文自己列出的优势包括：同一框架处理正问题、反问题与混合问题；不要求传感器位置在各次观测间固定；缓解多项式混沌的维数灾。

实验分三组：学习随机过程（含非高斯与混合非高斯场）、正向随机椭圆方程、反向随机椭圆方程。

> [!note] 可核实范围
> Karhunen-Loève 结构参考场的代数形式、流的耦合层形式与似然加物理的组合损失，本页依据的公开材料未能逐式确认。需要这些公式时应直接阅读原文第 3 节。

## 64：时间不是一个额外维度

### 关键困难

时间依赖 Fokker-Planck 方程的解是密度，定义在无界甚至高维区域上。要把流推广到时间方向，最自然的想法是把时间当作额外维度，即令 $\widehat x=(x,t)$、$\widehat z=(z,t^\ast)$，写

$$
p_{\widehat X}(\widehat x)=p_{\widehat Z}(\widehat z)\,|\det J|,
\qquad \widehat z=f(\widehat x).
$$

这样做是错的：质量沿时间轴不守恒，$\int p(x,t)\,\mathrm dx\,\mathrm dt\neq1$。论文的处理是把**潜时间钉在真实时间上**，$t^\ast=t$。Jacobian 随之退化为

$$
\det J=\begin{vmatrix}\partial z/\partial x & \partial z/\partial t\\ 0 & 1\end{vmatrix}
=\Bigl|\frac{\partial z(x,t)}{\partial x}\Bigr| ,
$$

于是

$$
p_{\widehat X}(x,t)=p_{\widehat Z}(z,t)\,\Bigl|\frac{\partial z}{\partial x}\Bigr|,
\qquad z=f(x,t).
$$

这句话值得单列：**流是时间条件化的，不是时间增广的。** 一个字之差决定了归一化是否成立。

### 时间依赖仿射耦合层

结构是 KRnet 从空间域到时空域的简化推广。每层由 Actnorm 层

$$
y_{[i]}=a_i\odot x_{[i]}+b_i
$$

接一个时间依赖仿射耦合层构成：

$$
x_{[i],1}=x_{[i-1],1},
$$

$$
x_{[i],2}=x_{[i-1],2}\odot\Bigl(\mathbf 1_{d-m}
+\beta\tanh\bigl(s_i(x_{[i-1],1},t)\bigr)\Bigr)
+e^{\zeta_i}\odot\tanh\bigl(q_i(x_{[i-1],1},t)\bigr),
$$

其中 $|\beta|<1$ 由用户指定（例如 $0.6$），$s_i,q_i:\mathbb R^{m+1}\to\mathbb R^{d-m}$，$\zeta_i\in\mathbb R^{d-m}$ 可训练。最后一层用多项式样条变换提高表达力。

与标准 real NVP 的差别有两处：$\tanh$ 把尺度因子限制在 $(1-\beta,1+\beta)$ 内，避免尺度爆炸；$s_i,q_i$ 以 $t$ 为额外输入，使同一组权重覆盖整个时间区间。

损失是时间依赖 Fokker-Planck 残差，无标签数据，也不生成随机微分方程的样本轨道。实验覆盖线性漂移、非线性漂移与高维问题。

## 72：分数阶算子的两种处理

编号 72 把这条路线推到分数阶 Fokker-Planck 方程。随机微分方程含 Lévy 噪声：

$$
\mathrm d X_t=\mu(X_t,t)\,\mathrm dt+\sigma(X_t,t)\,\mathrm dW_t+\mathrm dL^\alpha_t,
\qquad
\frac{\partial p}{\partial t}=\mathcal L p-(-\Delta)^{\alpha/2}p,
$$

$$
\mathcal L p=-\nabla\cdot(p\mu)+\tfrac12\nabla\cdot\nabla\cdot(\sigma\sigma^{\mathsf T}p).
$$

三个困难同时出现：无界区域、维数、非局部性。流解决前两个，第三个有两条路。

### 路线一：Monte Carlo（MCNF）

直接搬用[[computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|编号 66]]的内外分裂估计，作为本文的 Lemma 3.1，损失为

$$
L(p_{\mathrm{KRnet},\theta})=\frac{1}{N_S}\sum_{i=1}^{N_S}\bigl|R_\theta(x^i)\bigr|^2,
\qquad
R_\theta(x)=\bigl(\mathcal L-(-\Delta)^{\alpha/2}\bigr)p_{\mathrm{KRnet},\theta}(x).
$$

一个容易忽略但实际重要的观察：$\mathrm{Beta}(a,1)$ 在 $a\to0$ 时向原点集中，因此 $\alpha$ 越接近 $2$，内层半径样本 $r_1$ 越堆积在零附近，需要**更大**的下界保护 $r_\epsilon$ 才稳定。这把「参数取值影响数值稳定性」写成了可操作的规则。

> [!note] 与编号 66 的一处细微差别
> 编号 72 在内层与外层分别使用两个**独立**的方向变量 $\xi$ 与 $\eta$，而编号 66 的对应公式复用同一个 $\xi$。差别很小但确实存在。

### 路线二：解析辅助模型（GRBFNF）

第二条路避开随机性：换一个分数阶 Laplacian 已知闭式的模型。对高斯 $u(x)=\exp(-\sigma^{-2}|x-x_0|_2^2)$，

$$
(-\Delta)^{\alpha/2}u(x)=c_{\alpha,d}\,|\sigma|^{-\alpha}\;
{}_1F_1\!\Bigl(\frac{d+\alpha}{2};\frac{d}{2};-\sigma^{-2}|x-x_0|^2_2\Bigr),
\qquad
c_{\alpha,d}=\frac{2^\alpha\Gamma\bigl(\frac{d+\alpha}{2}\bigr)}{\Gamma\bigl(\frac d2\bigr)},
$$

${}_1F_1$ 是合流超几何函数。于是引入一个高斯径向基混合作为辅助模型

$$
p_{\mathrm{GRBF},\tilde\theta}(x)=\sum_{i=1}^M w_i\,\mathcal N(\tilde x_i,\sigma_i^2\mathbf I)(x),
\qquad 0\le w_i\le1,\ \sum_i w_i=1,
$$

其中 $w_i$ 与 $\sigma_i$ 都可训练，它的分数阶 Laplacian 有闭式：

$$
(-\Delta)^{\alpha/2}p_{\mathrm{GRBF},\tilde\theta}(x)
=c_{\alpha,d}\,\pi^{-d/2}2^{-\frac{d+\alpha}{2}}\sum_{i=1}^M w_i\,|\sigma_i|^{-(d+\alpha)}\,
{}_1F_1\!\Bigl(\frac{d+\alpha}{2};\frac d2;-\frac{|x-\tilde x_i|_2^2}{2\sigma_i^2}\Bigr).
$$

耦合损失把流与辅助模型绑在一起，第二项是一致性罚：

$$
\tilde L\bigl(p_{\mathrm{KRnet},\theta},p_{\mathrm{GRBF},\tilde\theta}\bigr)
=\frac{1}{N_S}\sum_{i}\Bigl(\mathcal L p_{\mathrm{KRnet},\theta}(x^i)
-(-\Delta)^{\alpha/2}p_{\mathrm{GRBF},\tilde\theta}(x^i)\Bigr)^2
+\frac{\beta_m}{N_S}\sum_{i}\Bigl(p_{\mathrm{KRnet},\theta}(x^i)
-p_{\mathrm{GRBF},\tilde\theta}(x^i)\Bigr)^2 .
$$

这是全篇最独特的想法：**把难算的算子换到一个算子已知的模型上，再用一致性罚把两个表示拉在一起。** 代价是径向基中心数随维数增长，因此这条路在低维更合适，而 Monte Carlo 路线在高维更合适。

### 时间依赖版本中的一个细节

时间依赖模型（MCTNF）的耦合层是

$$
x_{[i],2}=x_{[i-1],2}\odot\Bigl(\mathbf 1_{d-m}
+\beta\tanh\bigl(t\,s_{i,t}(x_{[i-1],1},t)\bigr)\Bigr)
+e^{\zeta_i}\odot\tanh\bigl(t\,q_{i,t}(x_{[i-1],1},t)\bigr).
$$

注意网络输出前显式乘了一个 $t$。这让该层在 $t=0$ 时退化为恒等映射，从而**初值条件被构造性地满足**，而不是通过一个罚项去逼近。这类「把约束写进架构」的做法与推前公式保证归一化是同一思路。

训练采用交替细化：用当前流采新的训练集，再在新集合上继续训练。实验含仅由分数阶 Laplacian 驱动的二维方程、双峰目标分布、更高维稳态方程，以及一个解为 Cauchy 分布的时间依赖问题（$\alpha$ 稳定噪声的自然稳态律）。

> [!note] 作者顺序
> 主页列出的作者顺序为 Xiaoliang Wan、Li Zeng、Tao Zhou；预印本与期刊版均为 **Li Zeng、Xiaoliang Wan、Tao Zhou**。本站按出版版本记录。

## 87：有界支撑需要另一种耦合层

### 为什么不能只是截断

KRnet 与一般正规化流映到 $\mathbb R^d$ 上的高斯参考，诱导密度支撑无界。对真正支撑在有界集上的密度（论文举的日常例子是人的年龄），以及定义在超矩形上的 PDE 解，这是错的。而把流硬性限制在盒子里会破坏精确可逆性。

### CDF 耦合层

本文的核心新对象是一个把 $[-1,1]$ 严格映到 $[-1,1]$ 的单调映射。先在网格 $-1=s_0<s_1<\cdots<s_n=1$ 上定义一个分段线性密度

$$
p(s)=\frac{w_{i+1}-w_i}{h_i}(s-s_i)+w_i,
\quad s\in[s_i,s_{i+1}],
\quad p(s_i)=w_i\ge0,\ h_i=s_{i+1}-s_i,
$$

其累积分布函数是**分段二次**的：

$$
F(s)=\frac{w_{i+1}-w_i}{2h_i}(s-s_i)^2+w_i(s-s_i)
+\sum_{k=0}^{i-1}\frac{w_k+w_{k+1}}{2}(s_{k+1}-s_k).
$$

它的逆是一个二次方程的根。记 $q_0=0$、$q_i=\sum_{k=0}^{i-1}\frac{w_k+w_{k+1}}{2}(s_{k+1}-s_k)$，则

$$
F^{-1}(q)=s_i+\frac{-w_i+\sqrt{w_i^2+2(w_{i+1}-w_i)(q-q_i)/h_i}}{(w_{i+1}-w_i)/h_i}
=s_i+\frac{2(q-q_i)}{w_i+\sqrt{w_i^2+2(w_{i+1}-w_i)(q-q_i)/h_i}} .
$$

第二种写法在分母接近零时数值稳定，这是实现层面的必要细节。令 $\tilde F(s)=2F(s)-1$ 得到 $[-1,1]\to[-1,1]$ 的双射，逐分量作用。耦合层对 $y=(y_1^{\mathsf T},y_2^{\mathsf T})^{\mathsf T}\in[-1,1]^l$、$y_1\in\mathbb R^m$ 为

$$
\hat y_1=y_1,
\qquad
\hat y_2=F\bigl(y_2;\theta(y_1)\bigr),
$$

参数 $\theta=(s_1^{\mathsf T},\dots,s_{n-1}^{\mathsf T},w_0^{\mathsf T},\dots,w_n^{\mathsf T})^{\mathsf T}\in\mathbb R^{2n(l-m)}$ 由一个网络 $\mathrm{NN}(y_1)$ 产生，并重参数化以保证网格递增与权非负：

$$
s_1=-1+\frac2n\bigl(1+\beta_1\tanh(\hat s_1)\bigr),
\qquad
s_{i+1}=s_i+\Bigl(\frac{1-s_i}{n-i}\Bigr)\bigl(1+\beta_{i+1}\tanh(\hat s_{i+1})\bigr),
\qquad
w_i=\frac{1+\gamma_i\tanh(\hat w_i)}{C},
$$

$C$ 是归一化常数，实现中取 $\beta_1=\tfrac{65}{66}$、$\beta_i=0.97$（$i\ge2$）、$\gamma_i=0.99$。这些常数把变换限制在严格单调的内部，避免退化。之后交换 $y_1,y_2$ 的角色再作用一层，按层号奇偶交替。

### 与 KRnet 的四点差别

| 项目       | KRnet                  | B-KRnet                      |
| ---------- | ---------------------- | ---------------------------- |
| 参考分布   | $\mathbb R^d$ 上的高斯 | $[-1,1]^d$ 上的均匀分布      |
| 耦合层     | 仿射（尺度加平移）     | 单调分段二次累积分布函数     |
| 尺度偏置层 | 有                     | **去掉**（区域有界，不需要） |
| 光滑性     | 由激活决定             | 密度分段线性，只有一阶可微   |

最后一条有直接后果：**二阶 PDE 必须改写成一阶系统。** 论文正是这样做的：引入 $g(x)=\nabla p$ 由单独网络 $g_{\mathrm{NN}}$ 逼近，损失为

$$
\mathcal L=\lambda_{pde}\mathcal L_{pde}+\lambda_b\mathcal L_b+\lambda_{g}\mathcal L_{g},
$$

$$
\mathcal L_{pde}=\mathbb E_{x\sim\rho}\bigl[|\mathcal N[x;p_{\text{B-KRnet},\theta},g_{\mathrm{NN}}]|^2\bigr],
\qquad
\mathcal L_{g}=\mathbb E_{x\sim\rho}\bigl[\|g_{\mathrm{NN}}(x)-\nabla p_{\text{B-KRnet},\theta}(x)\|_2^2\bigr].
$$

非负性与质量守恒**由构造成立**，不进入罚项。一般乘积区域 $\prod_i[a_i,b_i]$ 由线性映射 $y=\hat a\odot x+\hat b$（$\hat a_i=2/(b_i-a_i)$、$\hat b_i=-(b_i+a_i)/(b_i-a_i)$）归一到 $[-1,1]^d$。

### 伪三角结构与外层冻结

把 $x$ 分块为 $x=((x^{(1)})^{\mathsf T},\dots,(x^{(K)})^{\mathsf T})^{\mathsf T}$，$\sum_i d_i=d$，Knothe-Rosenblatt 重排的三角结构放松为块三角：

$$
z=f_{\text{KR}}(x)=
\begin{pmatrix}\tilde f_1\\ x^{(2:K)}\end{pmatrix}\circ
\begin{pmatrix}\tilde f_2\\ x^{(3:K)}\end{pmatrix}\circ\cdots\circ
\begin{pmatrix}\tilde f_{K-1}\\ x^{(K)}\end{pmatrix}\circ \tilde f_K(x),
$$

外层每完成一个阶段就冻结一个块（论文称为使若干维度失活的压缩操作），内层再由若干 CDF 耦合层复合。一个边界情形值得注意：CDF 耦合层需要输入维数至少为 $2$，因此当 $x^{(1)}=x_1\in\mathbb R$ 时最外层取恒等映射，理由是此时 $x^{(1)}$ 已接近均匀，$2F(\cdot)-1$ 也接近恒等。

### 混合区域

动理学 Fokker-Planck 方程的解是位置 $x\in\Omega$（有界）与速度 $v\in\mathbb R^d$（无界）的联合密度。论文按 $p(x,v)=h(v|x)q(x)$ 分解并建模

$$
p_{\theta}(x,v)=h_{\text{KRnet},\theta_1}(v|x)\cdot q_{\text{B-KRnet},\theta_2}(x),
$$

即无界 KRnet 作为速度的条件流乘以 B-KRnet 作为位置的密度。这里的条件化机制正是编号 64 的时间条件化，只是把条件变量从时间换成位置——同一个技术零件在不同论文里承担不同角色。

### 自适应策略

$$
\rho(x)\leftarrow(1-\gamma)\rho(x)+\gamma\,p_{\text{B-KRnet},\theta^{*,k}}(x),
\qquad N_{new}=\gamma N_{pde},
$$

新点由逆流生成：$z^{i}\sim\mathrm{Unif}(-1,1)^d$、$x^{i}=(f_{\theta^{*,k}})^{-1}(z^{i})$，其余点保留。实现中取更新率 $\gamma=0.8$、每个 CDF 耦合层三个子区间、`Tanh` 激活、Adam 加阶梯衰减。

密度估计实验为环形区域、高斯混合、带孔洞的 logistic 分布；PDE 实验为四维 $-\Delta p+p=f$、稳态 Keller-Segel 系统（两个耦合密度）、二维稳态动理学 Fokker-Planck 方程。论文自陈的限制是：B-KRnet 的光滑性只到一阶，要提高需要更高次多项式；在深度 Ritz 框架内做自适应点更新留作开放问题——后者由[[computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|编号 80]]接手。

> [!note] 版本差异
> 预印本 v1 的摘要说 B-KRnet「由一系列活跃变换维数逐步减少的耦合层组成，受 Knothe-Rosenblatt 重排的三角结构启发」，并把参考测度称为 base 分布；v3 改为「把伪三角结构纳入正规化流模型」，并称之为 prior 分布。引用时应注明版本。

## 四篇的推进关系

| 编号 | 目标对象               | 参考分布           | 训练信号                   | 新增技术零件               |
| ---- | ---------------------- | ------------------ | -------------------------- | -------------------------- |
| 62   | 散点观测下的随机场     | 带 KL 结构的高斯场 | 对数似然（可加物理残差）   | 系数由网络参数化的 KL 展开 |
| 64   | 时间依赖 Fokker-Planck | 无界高斯           | 方程残差                   | 潜时间钉住真实时间         |
| 72   | 分数阶 Fokker-Planck   | 无界高斯           | 方程残差                   | Monte Carlo 或解析辅助模型 |
| 87   | 有界支撑的密度         | $[-1,1]^d$ 均匀    | 方程残差（改写成一阶系统） | 分段二次 CDF 耦合层        |

一条贯穿的判断：**结构性约束应当由架构保证，而不是由罚项逼近。** 归一化由变量替换保证，非负性由推前公式保证，初值条件由耦合层乘 $t$ 保证，有界支撑由 CDF 耦合层保证。每加一条这样的保证，损失函数里就少一个需要调权重的罚项。

## 覆盖核对

| 内容                                | 论文 | 覆盖状态                                   |
| ----------------------------------- | ---- | ------------------------------------------ |
| 带 KL 结构的参考场与三步构造        | 62   | 构造步骤与优势（限定可核实范围）           |
| 时间条件化与 Jacobian 退化          | 64   | 错误做法、正确做法与后果                   |
| 时间依赖仿射耦合层                  | 64   | 公式、$\beta$ 的作用、与 real NVP 的差别   |
| Lévy 噪声下的分数阶 FP 方程         | 72   | 方程、两条算子处理路线                     |
| Beta 分布集中与 $r_\epsilon$ 的关系 | 72   | 参数依赖的稳定性规则                       |
| 高斯的解析分数阶 Laplacian 与 GRBF  | 72   | 闭式公式、混合模型、耦合损失与一致性罚     |
| 乘 $t$ 使初值条件构造性成立         | 72   | 耦合层形式与其含义                         |
| 分段线性密度与分段二次 CDF          | 87   | 密度、CDF、稳定的逆、重参数化与数值常数    |
| 与 KRnet 的四点差别                 | 87   | 参考分布、耦合层、去掉尺度偏置、一阶光滑性 |
| 一阶系统改写与三项损失              | 87   | 引入 $g$、三项损失、构造性满足的两个性质   |
| 伪三角结构、边界情形与混合区域      | 87   | 块分解、外层冻结、$d_1=1$ 情形、条件流乘积 |

## 本页原文

- L. Guo, H. Wu, and T. Zhou, _Normalizing field flows: solving forward and inverse stochastic differential equations using physics-informed flow models_, J. Comput. Phys. 461 (2022), 111202（预印本 [arXiv:2108.12956](https://arxiv.org/abs/2108.12956)）。
- X. Feng, L. Zeng, and T. Zhou, [_Solving time dependent Fokker-Planck equations via temporal normalizing flow_](https://doi.org/10.4208/cicp.OA-2022-0090), Commun. Comput. Phys. 32(2) (2022), pp. 401-423（预印本 [arXiv:2112.14012](https://arxiv.org/abs/2112.14012)）。
- L. Zeng, X. Wan, and T. Zhou, [_Adaptive deep density approximation for fractional Fokker-Planck equations_](https://doi.org/10.1007/s10915-023-02379-z), J. Sci. Comput. 97 (2023), 68（预印本 [arXiv:2210.14402](https://arxiv.org/abs/2210.14402)）。
- L. Zeng, X. Wan, and T. Zhou, _Bounded KRnet and its applications to density estimation and approximation_, SIAM J. Sci. Comput. 47(6) (2025), pp. C1294-C1318（预印本 [arXiv:2305.09063](https://arxiv.org/abs/2305.09063)）。
