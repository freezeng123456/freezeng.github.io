---
title: 采样器、滤波器与函数空间流
description: 编号 55、56、82、88、99、106：换掉采样器而不是代理
lang: zh
translation: en/computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering
tags:
  - 论文笔记
  - 贝叶斯反问题
  - 变分推断
---

> [!note] 本页覆盖
> 编号 **55**（_J. Comput. Math._ 39(6), 2021）、**56**（_Comput. Methods Appl. Mech. Engrg._ 386, 2021）、**82**（_Comput. Phys. Commun._ 311, 2025）、**88**（_Int. J. Mech. Sci._ 313, 2026）、**99**（投稿 _Math. Comput._，[arXiv:2411.13277](https://arxiv.org/abs/2411.13277)）、**106**（投稿 _SIAM J. Sci. Comput._，[arXiv:2605.29373](https://arxiv.org/abs/2605.29373)）。
>
> 编号 82 的正文（第 2 至 4 节）需订阅访问，本页对该篇只报告摘要、引言与节首片段可确认的内容，并对无法核实处明确标注。

[[computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|上一页]]的四篇固定采样器、替换代理。本页六篇反过来：代理的角色相对稳定，变的是**如何在后验上前进**。

## 55：把代理放进优化型提议

### randomize-then-optimize 的结构

RTO 把后验采样变成重复求解一个被随机扰动的非线性最小二乘问题。在白化坐标下（先验 $u\sim\mathcal N(u_{\mathrm{pr}},\Gamma_{\mathrm{pr}})$，$u=S_{\mathrm{pr}}v+u_{\mathrm{pr}}$），线性化点取

$$
v_{\mathrm{ref}}=\arg\min_{v}\tfrac12\|H(v)\|^{2},
$$

提议由

$$
v^{(i)}_{\mathrm{prop}}=\arg\min_{v}\ \tfrac12\bigl\|\bar Q^{T}H(v)-\xi^{(i)}\bigr\|^{2},
\qquad \xi^{(i)}\sim\mathcal N(0,I_n)
$$

给出。这个提议的密度是显式的：

$$
\pi_{\mathrm{RTO}}(v)=(2\pi)^{-n/2}\bigl|{\det}\bigl(Q^{T}\nabla H(v)\bigr)\bigr|
\exp\Bigl(-\tfrac12\bigl\|Q^{T}H(v)\bigr\|^{2}\Bigr),
$$

因此可以做独立提议的 Metropolis 修正。接受比化简为一个权的比值：

$$
\frac{\pi_{\mathrm{tar}}(v^{(i)}_{\mathrm{prop}})\,\pi_{\mathrm{RTO}}(v^{(i-1)})}
{\pi_{\mathrm{tar}}(v^{(i-1)})\,\pi_{\mathrm{RTO}}(v^{(i)}_{\mathrm{prop}})}
=\frac{w\bigl(v^{(i)}_{\mathrm{prop}}\bigr)}{w\bigl(v^{(i-1)}\bigr)},
$$

$$
w(v)=\bigl|{\det}\bigl(Q^{T}\nabla H(v)\bigr)\bigr|^{-1}
\exp\Bigl(-\tfrac12\|H(v)\|^{2}+\tfrac12\|Q^{T}H(v)\|^{2}\Bigr).
$$

RTO 的优点是提议独立，自相关远低于随机游走；代价是每个样本都需要多次正模型与**其 Jacobian** 的求值，既在优化器内部，也在权中。

论文同时复述了 Bardsley 等人的可扩展 RTO：对 $\nabla f(v_{\mathrm{ref}})=\Psi\Lambda\Phi^{T}$ 取秩 $r$ 的简化奇异值分解，拆分 $v=\Phi v_r+v^{\perp}$，则

$$
v^{\perp}=(I_n-\Phi\Phi^{T})\xi,
\qquad
v_r=\arg\min_{z}\bigl\|(\Lambda^{2}+I_r)^{-1/2}z+\Lambda\Psi^{T}f(v^{\perp}+\Phi z)-\Phi^{T}\xi\bigr\|^{2},
$$

把 $n$ 维优化换成 $r$ 维优化，行列式相应简化为

$$
\bigl|{\det}(\tilde Q^{T}\nabla H(v))\bigr|
=\bigl|{\det}(\Lambda^{2}+I_r)^{-1/2}\bigr|\cdot
\bigl|{\det}\bigl(I_r+\Lambda\Psi^{T}\nabla f(v)\Phi\bigr)\bigr| .
$$

### 本文的实质贡献：训练点取自近似后验

代理换成神经网络后，$\nabla_v\mathcal{NN}(v;\theta)$ 由反向传播给出，优化与权都变便宜。但真正决定成败的是训练点的位置：论文从一个**近似后验** $\widetilde\pi_{\mathrm{pos}}$ 抽训练点，而不是从先验抽。消融基线 "NN-RTO-pr" 就是同一算法改用先验抽点。

与[[computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates|上一页]]四篇的区别在于：这里的代理**离线训练一次，不做在线细化**。因此它是这一族里的静态成员，也是讨论在线细化价值时的自然对照。

效率用按 CPU 时间折算的有效样本量衡量：

$$
\mathrm{ESS}=\frac{n_{\mathrm{samps}}}{1+2\sum_{k=1}^{K}\rho(k)} .
$$

### 数值证据

基准 Darcy 型椭圆反问题：从带噪的逐点压力观测恢复渗透率 $\kappa(x)$，源项 $f(x)=100\sin(\pi x_1)\sin(\pi x_2)$，双线性 Galerkin 有限元在 $40\times40$ 均匀网格上求解。三种方法对比：真模型 RTO、先验训练的 NN-RTO-pr、后验训练的 NN-RTO。网络 3 隐层 × 40 神经元，Adam 学习率 $5\times10^{-4}$，内层优化用信赖域反射 Newton 法。结论是 NN-RTO 在 $N=50$ 个训练点时就复现真模型 RTO 的一维与二维后验边缘分布，而 NN-RTO-pr 不能——这个消融正好隔离出「训练点放在后验附近」的价值。

论文不含定理；摘要中「收敛到直接 RTO」是论证性与经验性的表述，不是已证明的命题。

## 56：粒子流与在线细化的结合

### 为什么 SVGD 需要梯度

Stein 变分梯度下降用核化的函数梯度流把粒子集合推向目标。粒子更新是 $x_i\leftarrow x_i+\varepsilon\varphi(x_i)$，最优方向来自泛函优化

$$
\varphi^{*}=\arg\max_{\varphi\in\mathcal S}
\Bigl\{-\tfrac{d}{d\varepsilon}\mathrm{KL}\bigl(q_{[\varepsilon\varphi]}\,\|\,\pi\bigr)\Bigr\}.
$$

Stein 算子恒等式把这个导数写成期望：

$$
-\frac{d}{d\varepsilon}\mathrm{KL}\bigl(q_{[\varepsilon\varphi]}\|\pi\bigr)\Big|_{\varepsilon=0}
=\mathbb E_{x\sim q}\bigl[\mathcal A_{\pi}^{\top}\varphi(x)\bigr],
\qquad
\mathcal A_{\pi}^{\top}\varphi(x)=\nabla_x\log\pi(x)^{\top}\varphi(x)+\nabla_x^{\top}\varphi(x).
$$

在再生核 Hilbert 空间单位球 $\mathcal S=\{\varphi\in\mathcal H^{d}:\|\varphi\|_{\mathcal H^{d}}\le1\}$ 上取极大，最优方向有闭式：

$$
\varphi^{*}(\cdot)\ \propto\ \mathbb E_{x\sim q}
\bigl[\nabla_x\log\pi(x)\,\kappa(x,\cdot)+\nabla_x\kappa(x,\cdot)\bigr],
$$

对应的经验更新是

$$
x^{(l+1)}_{i}\leftarrow x^{(l)}_{i}+\varepsilon_l Q_l\bigl(x^{(l)}_i\bigr),
\qquad
Q_l(x)=\frac1N\sum_{j=1}^{N}
\Bigl[\nabla_{x^{(l)}_j}\log\pi\bigl(x^{(l)}_j\bigr)\kappa\bigl(x^{(l)}_j,x\bigr)
+\nabla_{x^{(l)}_j}\kappa\bigl(x^{(l)}_j,x\bigr)\Bigr].
$$

两项有清楚的分工：第一项是 $\log\pi$ 的核加权上升，起吸引作用；第二项是核的梯度，起排斥作用，防止粒子塌缩到一点。核取径向基 $\kappa(x,x')=\exp(-h\|x-x'\|^{2})$。步长用 AdaGrad，论文明确指出高维下常数步长会发散。

问题在于 $\nabla_x\log\pi(x)$：对 PDE 约束后验，这意味着每个粒子每步一次伴随求解。若梯度不可得（遗留代码、不可微求解器）或太贵，原始 SVGD 无法使用。

### 检验点与选点约束

论文用粒子均值作为设计点：

$$
x^{*}=\frac1M\sum_{i=1}^{M}x^{(t+1)}_{i},
$$

指标是该点处的相对 $\ell^2$ 误差

$$
\mathrm{err}(x^{*})=\frac{\bigl\|f(x^{*})-\tilde f(x^{*})\bigr\|_{2}}{\bigl\|f(x^{*})\bigr\|_{2}} .
$$

新训练点的选择方式值得单独指出：**不从新的随机抽样中取，而从现有粒子中取**，并带一个分离约束

$$
x_i=\arg\min_{x'\in X^{(t+1)}}\|x'-x^{*}\|_{2}
\quad\text{s.t.}\quad \|x'-x\|_{2}\ge R\ \ \forall x\in\mathcal D .
$$

约束是针对退化的：若不要求与已有设计 $\mathcal D$ 分离，点会聚集，训练问题条件变差。若无可行解则跳出循环并收缩 $R\leftarrow\rho R$（$0<\rho<1$）。在线重训用迁移学习，从上一轮网络出发。

外层循环的结构是：跑 $T$ 步内层 SVGD，再判定是否细化，重复至多 $I_{\max}$ 轮。在线真求解总量是 $N_{\mathrm{eval}}=\sum_t q_t$，其中 $q_t$ 是第 $t$ 轮实际接受的点数。

### 数值证据

三个问题。第一个是二维「双香蕉」后验：固定的先验训练网络（10 个训练点，3 隐层 × 20 神经元）把粒子推向**错误**的高概率区，而带局部细化的版本在 100 轮后铺开到真支撑上。第二个是二维热源反演，第三个是时间分数阶方程扩散系数估计。默认参数为 $Q=5$、$R=0.2$、容差 $10^{-2}$、$\rho=0.8$、$I_{\max}=30$ 外层、$T=10$ 内层、动量 $0.9$；样本质量用共享中位数带宽的最大均值差异评分。第三个例子中，带局部细化的版本在 $n_t=100$ 时用 50 次在线真求解达到误差 $0.1732$，而同等条件下固定网络为 $0.3968$；$n_t=500$ 时分别是 80 次求解、误差 $0.1155$ 与 $0.2941$。

论文摘要中「在不破坏收敛性的前提下在线细化」是设计目标的表述，本文未给出相应定理。

## 82：时间演化场景下的两级学习

> [!warning] 可核实范围
> 该文第 2 至 4 节需订阅访问。以下内容来自摘要、完整引言与节首片段，**不含可核实的公式**；算法轮廓由摘要与引言重构，未见于算法框。

### 问题

对时间依赖的参数化 PDE 做状态与参数的联合估计时，集合 Kalman 滤波每个成员每个同化步都要一次 PDE 求解。缩小集合会引入协方差的采样误差（局部化只能部分缓解）。经典降阶模型依赖线性模态叠加，而参数化动力系统的解流形常有缓慢衰减的 Kolmogorov 宽度，降维维数随之增长，收益消失。数据驱动的深度降阶模型不受线性子空间限制，但离线快照生成本身很贵。本文针对的正是最后这项代价。

### 两级结构

第一级在**粗网格**上做算子推断：给约化模型设定一个与物理方程结构对应的多项式形式，通过最小二乘回归从模拟数据学出约化算子，而不是做侵入式投影。快照在粗网格 $H$ 上生成（$H\gg h$，$h$ 为全阶模型的细网格），这既是离线便宜的原因，也是约化模型有偏的原因。

第二级用神经网络学**模型误差**：从时间与参数映到约化误差，依据是这类误差在整个时间域上相对稳定。修正后的对象嵌入集合 Kalman 滤波，修正作用在观测输出一侧，使 Kalman 更新看到去偏后的预测观测。

论文的新颖性表述是：这是第一个基于算子推断的集合 Kalman 方法，用于非线性时间依赖 PDE 的状态与参数联合估计。误差度量为 $\mathrm{rel}_{L^2}(\mu)=\|\mu-\mu^{*}\|_2/\|\mu^{*}\|_2$ 及状态的相应相对 $L^2$ 误差。

### 与本专题其他工作的位置关系

结构上这是编号 34 的成分替换版：廉价模型放进 Kalman 方法，再显式修正代理误差。区别在于**修正在何时学**：编号 34、37、49、56、79 都在线细化并由指标驱动；本文把差异作为时间与参数的函数**离线**学好，在线只应用不再细化，这一点与编号 55 同类。基准问题为 Burgers 方程、FitzHugh-Nagumo 模型与对流扩散反应系统，三方对比为本文方法、未修正的降阶集合 Kalman 滤波、以及全阶模型集合 Kalman 滤波。摘要给出的结论是在不损失精度的前提下取得可观加速；具体加速倍数在可核实材料中未见。

## 88：分层与分阶段处理不可辨识性

这一篇是本专题里唯一不属于 PDE 反问题的工作，恰因此有参照价值。

### 问题结构

把介观红细胞模型与实验对齐，是一个带两重结构困难的反问题。第一，数据来自不同实验平台（干涉显微给几何、光镊给拉伸、波动与松弛实验给黏弹性），跨平台差异互不一致。第二，参数无法从任何单一数据集同时辨识。此外正模型是耗散粒子动力学模拟，远达不到直接 MCMC 的速度。

待推参数为 $\vartheta_{\mathrm{in}}=(A_0,v,\mu_{sh},k_b,\eta_m)$，即目标表面积、约化体积、剪切刚度、弯曲刚度与膜黏度；模拟输出为 $y=(D_{eq},h_{\max},h_{\min},D_{ax},D_{tr},t_c,W_{fl})$，即平衡直径、最大与最小厚度、拉伸下的轴向与横向变形、松弛特征时间与膜波动统计量。观测模型把全部模型误差与不确定性归入一个零均值高斯项 $y_j=\text{model}_j(x;\vartheta)+\sigma\epsilon_j$。

### 三个结构性决定

**动态退火确定无应力基态。** 不同约化体积 $v$ 的细胞由 $v=0.64$ 的双凹形变形得到；$v$ 较大时固定三角网格残留较大键力，因此当 $|v-v_{\mathrm{desired}}|$ 超过容差时，把平衡键长逐步调向当前键长。这一步定义了四组实验共同参照的无应力基线。

**八个代理网络。** 四类实验 × 两种细胞（健康红细胞与恶性疟原虫感染细胞），每个三隐层、`tanh` 激活、均方误差损失、Adam 加阶梯学习率，各用 10000 组模拟结果训练；网络宽度按敏感性分析逐实验选择。摘要报告预测误差低于 $10^{-2}$。

**两阶段分层推断。** 第一阶段用平衡形态与拉伸实验的单层模型做分层推断，得到几何参数与剪切模量的稳定分布；第二阶段带入第一阶段信息，加上波动与松弛实验，辨识包含黏弹性参数的完整参数集。超参数编码跨数据集的信息共享。实验锚点取自公开文献：健康红细胞的 $D_{eq}=7.82\pm0.62\,\mu\mathrm m$、$h_{\max}=2.58\pm0.27\,\mu\mathrm m$、$h_{\min}=0.81\pm0.35\,\mu\mathrm m$；感染细胞由表面积与体积测量推得 $D_{eq}=6.9\,\mu\mathrm m$、$h_{\max}=h_{\min}=3.2\,\mu\mathrm m$。

结论是后验分布统计稳健，病理细胞被推断为更硬更黏，且跨平台数据融合缓解了单层推断无法处理的多源不确定性。本文无数学定理。

> [!note] 题名差异
> 预印本题为 _An RBC-MsUQ Framework for Red Blood Cell Morpho-Mechanics_，作者与摘要相同；期刊版题名为主页所列，本站按期刊版记录。

## 99：直接在函数空间做变分推断

### 核心困难是测度论的

离散后再做变分推断，会让算法行为随网格细化而退化。直接在函数空间做能修正这一点，但引出新障碍：无穷维空间上的概率测度一般**互相奇异**，先验经任意流变换后可能与先验之间不存在 Radon-Nikodym 导数，此时 KL 损失甚至没有定义。

无穷维 Bayes 公式给出 $\dfrac{d\mu}{d\mu_0}(u)=\dfrac{1}{Z_\mu}\exp(-\Phi(u))$，逼近族 $\mathcal M(\mathcal H_u)$ 中每个测度都被要求与先验 $\mu_0$ 等价，才能使

$$
D_{\mathrm{KL}}(\nu\|\mu)=\int_{\mathcal H_u}
\Bigl[\ln\frac{d\nu}{d\mu_0}(u)-\ln\frac{d\mu}{d\mu_0}(u)\Bigr]\nu(du)
$$

有意义。

### Theorem 2.3：等价性与显式 Radon-Nikodym 导数

对复合 $f_\theta=f^{(N)}_{\theta_N}\circ\cdots\circ f^{(1)}_{\theta_1}$，若每层双射并满足论文给出的条件，则 $\mu_{f_\theta}\sim\mu_0$，且

$$
\frac{d\mu_{f_\theta}}{d\mu_0}\bigl(f_\theta(u)\bigr)
=\prod_{n=1}^{N}\bigl|{\det}_{1}\bigl(Df^{(n)}_{\theta_n}(u_{n-1})\bigr)\bigr|^{-1}
\exp\Bigl(\tfrac12\bigl\langle f_\theta(u)-u,\ f_\theta(u)-u\bigr\rangle_{\mathcal H}
+\bigl\langle u,\ u-f_\theta(u)\bigr\rangle_{\mathcal H}\Bigr),
$$

其中 ${\det}_1$ 是 Carleman-Fredholm 行列式，这是无穷维下的正确对象。实际使用的损失分解把所有期望推回**先验**上：

$$
D_{\mathrm{KL}}(\mu_{f_\theta}\|\mu)
=\mathbb E_{\mu_0}\ln\frac{d\mu_{f_\theta}}{d\mu_0}\bigl(f_\theta(u)\bigr)
-\mathbb E_{\mu_{f_\theta}}\ln\frac{d\mu}{d\mu_0}(u),
$$

先验是可直接采样的，因此 Monte Carlo 梯度可算：每步从 $\mu_0$ 抽 $N$ 个样本、推过流、按

$$
\nabla_{\theta_k}L(\theta_k)\approx
\frac1N\sum_{i=1}^{N}\nabla_{\theta_k}\ln\frac{d\mu_{f_{\theta_k}}}{d\mu_0}\bigl(f_{\theta_k}(u_i)\bigr)
-\frac1N\sum_{i=1}^{N}\nabla_{\theta_k}\ln\frac{d\mu}{d\mu_0}\bigl(f_{\theta_k}(u_i)\bigr)
$$

更新，实验用 Adam。

### 四类流层与低秩的必要性

论文给出四种满足条件的具体变换：函数 Householder 流（线性，像空间一维）、函数投影变换流（线性，像空间 $M$ 维，因此比 Householder 更具表达力）、函数平面流（非线性）、函数 Sylvester 流（非线性）。线性与非线性族的可逆性分别由两条引理保证。

函数平面流的形式是 $f_n(u)=u+u_n h\bigl(\langle u,w_n\rangle_{\mathcal H_u}+b\bigr)$，其低秩参数化

$$
u_n=\sum_{i=1}^{r}\lambda_i\alpha_i\phi_i,
\qquad
w_n=\sum_{i=1}^{r}\lambda_i\beta_i\phi_i
$$

建立在先验协方差 $\mathcal C_0$ 的特征对 $\{\lambda_i,\phi_i\}$ 上，$\{\alpha_i,\beta_i\}$ 可训练。论文强调这里低秩不是计算上的便利而是**理论上的必要**：没有它，测度等价性失效，KL 损失不再有限。这是本文与神经算子里的低秩技巧的实质区别。

论文另给出一条离散不变性命题（若 $\mathcal H_u$ 连续嵌入 $C(D)$，则四类流的各层都是离散不变的），以及一个条件变体：固定正模型时不同观测数据诱导不同后验，无条件版本需按数据重训，条件版本把流条件在数据上并能处理不同维数的观测向量。

数值实验为三个反问题——一维光滑方程、二维稳态 Darcy 流、电阻抗成像——并以 pCN 作为 MCMC 基线，重点验证与理论一致、相对 pCN 的效率，以及经验上的离散不变性（同一个流在不同离散层次上行为一致）。

> [!note] 版本差异
> 预印本 v1 只报告两个反问题（光滑方程与稳态 Darcy 流），v2 与 v3 增加了电阻抗成像。本页依据 v3。

## 106：隐变量流与自适应先验

### 三个耦合的失效模式

高维 PDE 反问题中，非高斯甚至多峰后验让基于高斯近似的采样器（无迹与集合 Kalman 反演）乃至 pCN 都吃力；先验样本预训练的神经代理在后验集中到别处后落到分布外；先验均值通常给错，而手工调它正是希望消除的人工干预。标准正规化流在维数上帮不上忙，因为它是双射，保持维数。

### Variational Flow

架构是 VAE 式非线性降维（隐变量 $z\in\mathbb R^{k}$、数据 $x\in\mathbb R^{d}$、$k<d$）加**双流**：一个作用在隐变量先验上，一个把编码器从对角高斯换成条件正规化流。采样映射为

$$
z=f^{-1}_{\mathrm{pr},\beta}(v),\ v\sim\mathcal N(0,I)
\quad\Longrightarrow\quad
\xi=\mu_{\mathrm{de},\theta}(z)+\sigma_{\mathrm{de},\theta}(z)\odot\epsilon,\ \epsilon\sim\mathcal N(0,I),
$$

目标是未归一化后验 $\hat p(x)=\exp(-\Phi(\xi,y))\pi_0(\xi)$。

论文给出一条证据下界严格改进定理：在两个条件之一成立时，该模型的证据下界严格高于标准 VAE；两者同时成立时，流先验与条件流编码器各自贡献严格正的改进。证明用

$$
\mathrm{ELBO}=\text{const}-\mathbb E_{p_x}D_{\mathrm{KL}}\bigl(q_{z|x,\alpha}\,\|\,p_{z|x,\theta^{*},\tilde\beta}\bigr)
$$

分两步进行：先固定最优 VAE 编解码器、把隐变量先验换成流先验，再固定解码器与流先验、把编码器从对角高斯扩展为条件流。论文明确指出自适应循环本身缺少收敛保证。

### 先验均值的动量式更新

$$
\mu^{(k,i)}_{\mathrm{prior}}=\alpha\,\mu^{(k,i)}_{\mathrm{post}}+(1-\alpha)\,\mu^{(k-1)}_{\mathrm{prior}},
\qquad
\mu^{(k,i)}_{\mathrm{post}}=\frac1M\sum_{j=1}^{M}\xi^{(j)},
\quad \xi^{(j)}\sim p^{(k,i-1)}_{\mathrm{VF}}(\xi).
$$

两个时间尺度：后验均值每个 epoch 重估，锚点 $\mu^{(k-1)}_{\mathrm{prior}}$ 在一个阶段内固定、阶段结束时更新。协方差**刻意保持** $\Sigma_0$ 不变，以防模式塌缩；$\alpha$ 越小正则化越强，$\alpha=1$ 去掉正则化，论文引用的先前工作发现 $\alpha\approx0.5$ 接近最优，初值取 $\mu_0=0$、$\Sigma_0=I$。

### 代理微调与激进的数据替换

后验样本在用作训练数据前先加扰动：

$$
\hat\xi^{(j)}=\xi^{(j)}_{\mathrm{post}}+\gamma\,\nu^{(j)},
\qquad \nu^{(j)}\sim\mathcal N(0,I),\ \gamma>1,
$$

扰动施加在 Karhunen-Loève 系数空间，那里先验是标准正态，因此 $\gamma$ 的单位是先验标准差。论文明确与编号 79 对比：不累积历史数据、也不做代价较高的贪心筛选，而是**丢弃**旧数据，只在新采的局部集合上微调 Fourier 神经算子。

停机准则是先验均值处数据残差的相对变化：

$$
\frac{\bigl|\Phi\bigl(\mu^{(k-1)}_{\mathrm{prior}},y\bigr)-\Phi\bigl(\mu^{(k)}_{\mathrm{prior}},y\bigr)\bigr|}
{\Phi\bigl(\mu^{(k-1)}_{\mathrm{prior}},y\bigr)}<\epsilon .
$$

这与编号 79 的 $e_D$ 是同一类目标导向量，区别只在于检验点：那里是真模型选出的锚点，这里是演化中的先验均值。

### 数值证据

四个问题。第一个是 100 维 Rosenbrock 反问题，用于单独检验后验逼近质量，因此**刻意关闭**先验更新模块，与 VAE、MCMC、SVGD、无迹 Kalman 反演比较二维边缘；报告结论是本文模型最准，且显著优于原始 VAE，作者据此认为证据下界定理得到经验支持。其余三个是一维 Darcy 流、二维 Darcy 流与二维 Navier-Stokes（从 $T=1$ 的涡度场恢复初始涡度），均在 Karhunen-Loève 截断 $d\in\{32,64\}$ 与噪声 $\delta\in\{1\%,5\%,10\%\}$ 下运行、3 次取平均，与 pCN、SVGD、全阶与算子代理版无迹 Kalman 反演比较。报告结论是一维 Darcy 全面占优、高噪声最明显；二维 Darcy 在中高噪声多数情形最优，且真值取自均匀分布（出分布）时，先验更新是它找到真后验的原因，而缺少先验更新的 pCN 与 SVGD 停留在初始 $\mathcal N(0,I)$ 附近，把场的量级搞错；Navier-Stokes 在所有截断维数与噪声水平下反演误差最低。

> [!note] 定量表格
> 该文的具体误差数值在可用材料中不易可靠转录，本页只报告可确认的定性结论与实验配置。

## 六篇的对照表

| 编号 | 后验逼近方式                | 代理             | 细化时机         | 理论结果                 |
| ---- | --------------------------- | ---------------- | ---------------- | ------------------------ |
| 55   | 优化型独立提议 + Metropolis | 前馈网络         | 无（离线一次）   | 无                       |
| 56   | 粒子流                      | 前馈网络         | 在线，按粒子均值 | 无                       |
| 82   | 集合 Kalman 滤波            | 算子推断降阶模型 | 无（离线学误差） | 无                       |
| 88   | 分层分阶段 MCMC             | 八个前馈网络     | 无               | 无                       |
| 99   | 函数空间变分流              | 不用代理         | 不适用           | 测度等价性与显式 RN 导数 |
| 106  | 隐变量变分流                | Fourier 神经算子 | 在线，替换式     | 证据下界严格改进         |

## 覆盖核对

| 内容                         | 论文 | 覆盖状态                                       |
| ---------------------------- | ---- | ---------------------------------------------- |
| RTO 提议、密度与接受权       | 55   | 线性化点、提议、密度、权与可扩展变体           |
| 训练点取自近似后验的消融     | 55   | 与先验训练基线的对照及其结论                   |
| Stein 算子与 SVGD 闭式方向   | 56   | 泛函优化、恒等式、闭式解与吸引/排斥读法        |
| 分离约束的选点规则           | 56   | 约束形式、退化动机与半径收缩                   |
| 两级降阶 + 误差网络          | 82   | 粗网格算子推断、误差网络、输出侧修正（含限定） |
| 分层两阶段架构与动态退火     | 88   | 参数与输出、退火、八个代理、两阶段分工         |
| 测度等价性、RN 导数与四类流  | 99   | 定理结论、损失分解、低秩必要性、离散不变性     |
| 隐变量流、先验更新与激进替换 | 106  | 采样映射、动量更新、扰动策略、停机准则         |

## 本页原文

- L. Yan and T. Zhou, [_An acceleration strategy for randomize-then-optimize sampling via deep neural networks_](https://doi.org/10.4208/jcm.2102-m2020-0339), J. Comput. Math. 39(6) (2021), pp. 848-864（预印本 [arXiv:2104.06285](https://arxiv.org/abs/2104.06285)）。
- L. Yan and T. Zhou, [_Stein variational gradient descent with local approximations_](https://doi.org/10.1016/j.cma.2021.114087), Comput. Methods Appl. Mech. Engrg. 386 (2021), 114087（预印本 [arXiv:2104.06276](https://arxiv.org/abs/2104.06276)）。
- Y. Wang, L. Yan, and T. Zhou, [_Deep learning-enhanced reduced-order ensemble Kalman filter for efficient Bayesian data assimilation of parametric PDEs_](https://doi.org/10.1016/j.cpc.2025.109544), Comput. Phys. Commun. 311 (2025), 109544.
- S. Wang, L. Ma, L. Guo, X. Li, and T. Zhou, [_Multi-stage uncertainty quantification framework for red blood cell morpho-mechanics_](https://doi.org/10.1016/j.ijmecsci.2026.111352), Int. J. Mech. Sci. 313 (2026), 111352（预印本 [arXiv:2508.06852](https://arxiv.org/abs/2508.06852)）。
- Y. Zhao, H. Lu, J. Jia, and T. Zhou, _Functional normalizing flow for statistical inverse problems of partial differential equations_, [arXiv:2411.13277](https://arxiv.org/abs/2411.13277)，投稿 Math. Comput.；参考实现 [jjx323/FunctionalNormalizingFlow](https://github.com/jjx323/FunctionalNormalizingFlow)。
- Y. Wang, X. Wang, K. Tang, X. Wan, T. Zhou, and C. Yang, _Deep adaptive dimension reduction for Bayesian inference in inverse problems_, [arXiv:2605.29373](https://arxiv.org/abs/2605.29373)，投稿 SIAM J. Sci. Comput.
