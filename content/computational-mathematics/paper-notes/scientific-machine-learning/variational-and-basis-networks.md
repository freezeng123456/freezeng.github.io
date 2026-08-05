---
title: 约束怎么进损失：增广 Lagrange 与神经网络基
description: 编号 60、90、102：把边界条件、局部低正则性与散度约束交给结构而不是罚项
lang: zh
translation: en/computational-mathematics/paper-notes/scientific-machine-learning/variational-and-basis-networks
tags:
  - 论文笔记
  - 科学机器学习
  - 神经网络基
---

> [!note] 本页覆盖
> 编号 **60**（_Commun. Comput. Phys._ 31(3), 2022）、**90**（_Commun. Comput. Phys._ 39(2), 2026）、**102**（投稿 _Commun. Comput. Phys._，[arXiv:2603.17906](https://arxiv.org/abs/2603.17906)）。三篇的共同作者线索是黄建国与周涛，共同技术取向是**把网络训练还原成（一列）最小二乘问题**，并尽量让约束由结构而不是罚项承担。

## 60：用增广 Lagrange 取代罚项

### 直觉

带本质边界条件的变分问题对网络求解器不友好，原因很具体：网络函数无法精确满足 Dirichlet 条件，连在插值节点上也不行。论文点出两种通行做法各有缺陷。

一种是把约束写进拟设，$\varphi(x;\theta)=\ell(x)\psi(x;\theta)+\bar g(x)$。这在纸面上完美——约束恒等成立——但要求已知一个在边界上取零、内部为正的类水平集函数 $\ell$，以及边界值的一个已知延拓 $\bar g$。区域一旦复杂，这两样东西本身就是难题。

另一种是在边界上加二次罚项。这总是可行的，但它把精度换成了一个调参问题：罚参数太小则约束不被满足，太大则损失病态、优化器寸步难行。**罚参数没有一个自然的量级**，只能试。

编号 60 的判断是：这个两难在经典优化里早有标准答案，就是增广 Lagrange 法。引入一个乘子把约束的「一阶信息」显式携带起来，罚项只负责保证局部凸性。这样一来罚参数不需要趋于无穷，起始值也不必很大——**乘子承担了主要的约束力，罚项退居辅助**。这就是摘要所说「罚参数的选择灵活而稳健」的机制。

把这套经典方法搬到网络上，只有一处真正的障碍，也正是论文的技术核心：经典算法的乘子更新是**函数空间里的显式加法**，而网络的乘子是由参数决定的，参数空间没有对应的加法。论文的解决办法是把这一步改写成一次最小二乘拟合。

### 问题设定

设 $V$ 为 Hilbert 空间，$B:V\to W$ 有界线性（通常 $W=L^2(\Gamma)$），$V_g=\{v\in V:Bv=g\text{ on }\Gamma\}$，原问题为

$$
\min_{v\in V_g}J(v).
$$

增广 Lagrange 函数取

$$
\mathcal L_\beta(v,\mu)=J(v)-\langle\mu,\,Bv-g\rangle_W+\frac{\beta}{2}\|Bv-g\|_W^2 .
$$

注意乘子项前是**减号**，这是论文实际采用的符号约定。

三类实例的能量分别是：

**二阶椭圆问题** $-\operatorname{div}(A\nabla u)+cu=f$、$u|_\Gamma=g$，取 $V=H^1(\Omega)$、$W=L^2(\Gamma)$，

$$
J(v)=\tfrac12\int_\Omega\bigl[A\nabla v\cdot\nabla v+cv^2-2fv\bigr]\mathrm dx,
$$

$$
\mathcal L_\beta(v,\mu)=J(v)-\int_\Gamma\mu(x)\bigl(v(x)-g(x)\bigr)\mathrm dx
+\frac{\beta}{2}\int_\Gamma\bigl(v(x)-g(x)\bigr)^2\mathrm dx .
$$

论文援引其参考文献指出这一问题存在唯一鞍点。

**线性特征值问题** $-\nabla\cdot(p\nabla u)+qu=\rho u$、$u|_\Gamma=0$：通过归一化 $\tilde v=v/\|v\|_{0,\Omega}$ 把 Rayleigh 商化为极小化

$$
\tilde J(v)=\int_\Omega\bigl[p\nabla\tilde v\cdot\nabla\tilde v+q\tilde v^2\bigr]\mathrm dx,
\qquad
\mathcal L_\beta(v,\mu)=\tilde J(v)-\int_\Gamma\mu(x)\tilde v(x)\mathrm dx
+\frac{\beta}{2}\int_\Gamma\tilde v^2(x)\mathrm dx .
$$

**非线性（Gross-Pitaevskii 型）特征值问题** $-\nabla\cdot(A\nabla u)+Vu+u^3=\rho u$、$u|_{\partial\Omega}=0$、$\|u\|_{0,\Omega}=1$，能量加一个四次项：

$$
\tilde J(v)=\tfrac12\int_\Omega\bigl[A\nabla\tilde v\cdot\nabla\tilde v+V\tilde v^2+\tilde v^4\bigr]\mathrm dx,
$$

边界项形式不变。论文援引其参考文献指出，在 $1\le d\le3$ 时非负基态解唯一。

### 推导

**从约束问题到极小极大问题。** 对固定的 $v$ 作内层极大化。若 $v\in V_g$，则 $Bv-g=0$，两个含 $\mu$ 的项都消失，$\mathcal L_\beta(v,\mu)=J(v)$ 与 $\mu$ 无关；若 $v\notin V_g$，则 $-\langle\mu,Bv-g\rangle_W$ 关于 $\mu$ 线性且非平凡，沿 $\mu=-t(Bv-g)$、$t\to+\infty$ 可以让它无界增长。于是

$$
\max_{\mu\in W}\mathcal L_\beta(v,\mu)=
\begin{cases}
J(v), & v\in V_g,\\
+\infty, & \text{否则},
\end{cases}
$$

因此 $\min_{v\in V}\max_{\mu\in W}\mathcal L_\beta(v,\mu)$ 与原约束问题等价。对偶函数为 $F_\beta(\mu)=\min_{v\in V}\mathcal L_\beta(v,\mu)$，对偶问题为 $\max_{\mu\in W}F_\beta(\mu)$。

**网络参数化。** 原始与对偶变量各用一个 ResNet：$u\approx\varphi^u(x;\theta_u)$、$\lambda\approx\varphi^\lambda(x;\theta_\lambda)$，网络定义为

$$
h_0=Vx,
\qquad
h_\ell=h_{\ell-1}+\sigma(W_\ell h_{\ell-1}+b_\ell),
\qquad
\varphi(x;\theta)=a^{\mathsf T}h_L .
$$

**乘子更新的翻译。** 经典算法的乘子更新是

$$
\mu_{k+1}=\mu_k-\beta_k(Bv_k-g),
$$

右端是一个函数，而 $\mu$ 现在由 $\theta_\lambda$ 决定，参数空间里没有与这个函数加法对应的操作。论文把它改写成一次参数空间的最小二乘投影：

$$
\theta^\mu_{k+1}=\arg\min_{\theta_\nu}J_\lambda(\varphi^\nu;\varphi^\mu_k,\varphi^v_k),
\qquad
J_\lambda=\bigl\|\varphi^\nu-\varphi^\mu_k-\beta_k\bigl(B\varphi^v_k-g\bigr)\bigr\|^2_{L^2(\Gamma)},
$$

其中 $\varphi^\nu=\varphi^\lambda(x;\theta_\nu)$、$\varphi^\mu_k=\varphi^\lambda(x;\theta^\mu_k)$、$\varphi^v_k=\varphi^u(x;\theta^v_k)$。论文称之为**投影技术**。它把「函数空间的显式更新」翻译成「参数空间的拟合问题」：目标函数里的 $\varphi^\mu_k-\beta_k(B\varphi^v_k-g)$ 是一个可以逐点求值的固定函数，于是问题变成让乘子网络去逼近它。这一步是把经典增广 Lagrange 算法搬到网络上的关键。

**两个损失的 Monte Carlo 形式。** 以椭圆情形为例，

$$
\mathcal L_\beta(v,\mu)=|\Omega|\,\mathbb E_\xi\Bigl[\tfrac12 A(\xi)|\nabla v(\xi)|^2
+\tfrac12 c(\xi)v^2(\xi)-f(\xi)v(\xi)\Bigr]
-|\Gamma|\,\mathbb E_\eta\Bigl[\mu(\eta)\bigl(v(\eta)-g(\eta)\bigr)
-\tfrac{\beta}{2}\bigl(v(\eta)-g(\eta)\bigr)^2\Bigr],
$$

$$
J_\lambda(\nu;\mu_k,v_k)=|\Gamma|\,\mathbb E_\eta
\Bigl[\bigl(\nu(\eta)-\mu_k(\eta)-\beta(v_k(\eta)-g(\eta))\bigr)^2\Bigr],
$$

$\xi\sim\mathrm{Unif}(\Omega)$、$\eta\sim\mathrm{Unif}(\Gamma)$。

> [!note] 一处括号分组
> Monte Carlo 形式在可访问渲染中的 $-|\Gamma|\mathbb E_\eta[\cdot]$ 括号分组有歧义。**符号结构 $-\mu(v-g)+\tfrac{\beta}{2}(v-g)^2$ 由上面的连续形式确定**，引用时应以连续形式为准。

### 定理

论文给出的是无穷维层面的经典结论：约束变分问题与极小极大问题的等价性，以及相应的鞍点刻画。**这是被引用到无穷维设定中的经典增广 Lagrange 理论，不是关于神经网络的新定理。** 论文没有对离散（网络）方法证明收敛或逼近定理，这一点值得明说，因为它划定了本文主张的范围：方法是有原理依据的，但离散层面的保证不在其中。

### 数值实验

**算法 2（ALDL）的结构：**

1. 输入 $\beta_0>0$、增长系数 $\alpha\ge1$、外层轮数 `Epoch`、内层轮数 `Epoch_u` 与 `Epoch_λ`；
2. 用 PyTorch 默认随机初始化 $\theta^v_0$、$\theta^\mu_0$；
3. 对 $k=0,1,\dots$：固定 $\varphi^v_k,\varphi^\mu_k$，从 $\theta^\mu_k$ 出发用 Adam 迭代 `Epoch_λ` 步极小化 $J_\lambda$ 得 $\theta^\mu_{k+1}$；固定 $\beta_k,\varphi^\mu_{k+1}$，从 $\theta^v_k$ 出发用 Adam 迭代 `Epoch_u` 步极小化 $\mathcal L_\beta$ 得 $\theta^v_{k+1}$；令 $\beta_{k+1}=\alpha\beta_k$；
4. 输出 $u=\varphi^v_{\texttt{Epoch}+1}$、$\lambda=\varphi^\mu_{\texttt{Epoch}+1}$。

它离散的是一个无穷维算法（论文的 Algorithm 1）：先更新 $\mu_{k+1}=\mu_k-\beta_k(Bv_k-g)$，再求 $\mathcal L_\beta(v,\mu_{k+1})$ 的近似极小点 $v_{k+1}$ 使 $\|\partial_1\mathcal L_\beta(v_{k+1},\mu_{k+1})\|\le\tau_k$，然后 $\beta_{k+1}=\alpha\beta_k$ 并选新的容差 $\tau_{k+1}$。

**罚参数是递增的**，但因为乘子承担了大部分约束力，起始 $\beta_0$ 不必很大——这正是「罚参数选择灵活」的实现方式，而不只是一句宣称。

第 4 节的算例与基线：

| 算例                     | 维数    | 内容                       |
| ------------------------ | ------- | -------------------------- |
| 二阶椭圆（Poisson 型）   | 2D、3D  | 带本质边界条件             |
| 线性特征值问题           | 2D、3D  | 最小特征值与相应特征函数   |
| 非线性特征值问题         | 2D、3D  | Gross-Pitaevskii 型        |

基线有两个：**罚方法的深度学习版本**，以及**直接用随机梯度下降上升（SGDA）解极小极大问题**。第二个基线的存在很关键——它检验的是「用增广 Lagrange 加交替更新」是否真的比「照着极小极大的字面意思做」更好。论文在 Remark 3.1 中明确指出，为这个极小极大问题设计高效的下降上升法**仍是开放问题**，这既解释了为何要用交替方案，也承认了这条路线的理论基础尚不完整。

定性结论是：ALDL 对 $\beta$ 的敏感性显著低于罚方法，并在相当的计算代价下达到更好的精度。

> [!warning] 数值结果的可核实范围
> 算例的构成、两个基线与定性结论已核实；**逐例误差表未能从可访问渲染中确认**。流传的具体倍数（诸如特征函数改进 2 至 20 倍、特征值改进至多 100 倍、快 25% 至 30%）**均未核实，不应引用**。

### 与其他论文的关系

这是本主题里面向约束优化的一篇，与以采样为中心的几篇构成对照。它与[[computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|编号 80]]共享同一个深度 Ritz 变分损失设定：编号 60 修的是**边界约束如何进入目标函数**，编号 80 修的是**求积点放在哪里**。它与编号 90、102 共享作者线索与「把训练变成最小二乘」的取向，但手段相反——编号 60 引入乘子去执行约束，编号 90 与 102 则彻底取消罚权。

## 90：把低正则性局部化，再给局部换一个尺度

### 直觉

随机基 PDE 求解器（极限学习机、随机特征方法、可迁移神经网络）的做法是：把隐层权重与偏置随机取定后**冻结**，只对输出层系数做最小二乘。这一步把 PDE 求解从非凸优化变成了线性最小二乘，因此解光滑时精度极高，能达到离散层面而非训练层面的误差水平。

论文直接指出它的失效模式：精确解正则性低时精度显著恶化。而且失效的方式很有针对性——**低正则性通常是局部的**。尖峰、凹角奇性：解在绝大部分区域上光滑，只在很小的邻域里有极大的导数。此时全局统一尺度的随机基是两头不讨好的：为了分辨尖峰要把激活的有效波长压到峰宽的量级，可这样一来同一组基在光滑区域就变成了高频噪声，条件数崩坏。

论文的对策是拒绝在「全局尺度」这一个自由度上妥协，转而**让不同区域用不同尺度的基**。三个部件把这句话变成算法：用残差找到需要细化的位置；在那里切出一个子域并给它一组重定心、重缩放的基；用非重叠区域分解把各子域的解拼起来，界面上按 $C^1$ 匹配。

值得一提的是这与自适应采样的对偶关系：**残差大的地方，自适应采样加点，这里加的是一组重缩放的局部基。** 同一个指标，两种响应。

论文把「低正则」定义为不属于 $H^2$，或导数在某些点附近非常大。

### 问题设定

浅层网络写成基的形式：

$$
u_{NN}(x)=\sum_{m=1}^{M}\alpha_m\,\sigma\bigl(w_m^{\top}x+b_m\bigr)+\alpha_0,
$$

由于 $\{w_m,b_m\}$ 预设并冻结，$u_{NN}$ 落在 $V_M=\mathrm{span}\{\psi_0,\dots,\psi_M\}$，其中

$$
\psi_0(x)=1,
\qquad
\psi_m(x)=\sigma\bigl(\omega_m^{\top}x+b_m\bigr),
\quad 1\le m\le M .
$$

预设隐层有两种方式。第一种是在 $\Omega=[-1,1]^d$ 上从 $[-R,R]^d$ 均匀抽 $w_m$、从 $[-R,R]$ 抽 $b_m$，$R$ 是用户参数。第二种是可迁移神经网络的重参数化：把每个神经元写成

$$
\sigma\bigl(w_m^{\top}x+b_m\bigr)=\sigma\bigl(\gamma_m(a_m^{\top}x+r_m)\bigr),
\qquad \|a_m\|_2=1,\ \gamma_m\ge0,
$$

对应关系是 $w_m=\gamma_ma_m$、$b_m=\gamma_mr_m$，反之 $a_m=w_m/\|w_m\|_2$、$r_m=b_m/\|w_m\|_2$、$\gamma_m=\|w_m\|_2$。**位置**参数取

$$
a_m=\frac{X_m}{\|X_m\|_2},\qquad r_m=U_m,\quad m=1,\dots,M,
$$

$X_m$ 是独立 $d$ 维标准高斯、$U_m$ 是 $[0,1]$ 上独立均匀，这使分割超平面在单位球 $B_1(0)$ 内均匀分布；**形状**参数共享，$\gamma_m\equiv\gamma$，并用高斯随机场的实现作为辅助函数调节，**不使用任何 PDE 信息**。

这个分解是本文的先决条件：它把「随机基好不好」拆成两个可分别控制的部分——超平面位置的覆盖性（由 $a_m,r_m$ 决定）与激活的有效波长（由 $\gamma_m$ 决定）。有了这个分离，「给局部换一个尺度」才有明确含义，就是只动 $\gamma$。

### 推导

**非重叠区域分解与无权最小二乘。** 取 $\bar\Omega=\cup_{k=0}^{K}\bar\Omega_k$ 内部互不相交，每个 $\Omega_k$（$k\ge1$）只与 $\Omega_0$ 沿 $\Gamma_k=\partial\Omega_k\cap\partial\Omega_0$ 相接。原边值问题 $Lu=f$、$u|_{\partial\Omega}=g$ 化为带 $C^1$ 传输条件的 $(K{+}1)$ 子域系统：

$$
Lu_k=f\ \text{in}\ \Omega_k,
\qquad u_k=g\ \text{on}\ \partial\Omega_k\cap\partial\Omega,
$$

$$
u_k=u_0\ \text{on}\ \Gamma_k,
\qquad
\frac{\partial u_k}{\partial n_k}=\frac{\partial u_0}{\partial n_k}\ \text{on}\ \Gamma_k,
\quad 1\le k\le K .
$$

以 $\tilde u_k(\alpha_k,x)=\sum_{m=0}^{M_k}\alpha_{m,k}\psi_{m,k}(x)$ 展开，并在 $3K{+}2$ 组配点上取残差，配点最小二乘损失为

$$
\min_{\alpha_0,\dots,\alpha_K}\
\sum_{k=0}^{K}\Bigl(\sum_{x\in X_{f_k}}\bigl|L\tilde u_k(x)-f(x)\bigr|^2
+\sum_{x\in X_{g_k}}\bigl|\tilde u_k(x)-g(x)\bigr|^2\Bigr)
+\sum_{k=1}^{K}\sum_{x\in X_{\Gamma_k}}
\Bigl(\bigl|\tilde u_k(x)-\tilde u_0(x)\bigr|^2
+\Bigl|\tfrac{\partial\tilde u_k}{\partial n_k}(x)-\tfrac{\partial\tilde u_0}{\partial n_k}(x)\Bigr|^2\Bigr).
$$

**全部残差项以系数 $1$ 相加，没有罚权。** 这与带权 PINN 损失形成鲜明对照，而且不是风格差异：因为基已冻结，线性 $L$ 下这是一个线性最小二乘问题 $\min_\alpha\|F\alpha-T\|_2^2$，可以直接求解；最小二乘的解不依赖于对残差整体乘一个常数，各项之间的相对权重才是唯一自由度，而论文选择让它们都为 $1$。换句话说，**罚权之所以能取消，是因为求解器换了。**

**半线性算子的 Gauss-Newton。** 令 $u_k^{n+1}=u_k^{n}+v_k^{n}$，增量 $v^n$ 解线性化的传输问题

$$
DL(u_k^n;v_k^n)(x)=f(x)-Lu_k^n(x)\ \text{in}\ \Omega_k,
$$

$$
v_k^n=g-u_k^n\ \text{on}\ \partial\Omega_k\cap\partial\Omega,
\qquad
v_k^n-v_0^n=u_0^n-u_k^n\ \text{on}\ \Gamma_k,
$$

$$
\frac{\partial v_k^n}{\partial n_k}-\frac{\partial v_0^n}{\partial n_k}
=\frac{\partial u_k^n}{\partial n_k}-\frac{\partial u_0^n}{\partial n_k}\ \text{on}\ \Gamma_k,
$$

其中 $DL(u;v)$ 是 $L$ 在 $u$ 处沿 $v$ 的 Gâteaux 导数。把增量在同一冻结基上展开，每个 Newton 步同样是线性最小二乘 $\min_{a^n}\|F^na^n-T^n\|_2^2$。停机准则用残差范数的相对变化，

$$
R_{\text{emse}}=\frac{\bigl|\|F^na^{n,*}-T^n\|_2^2-\|F^{n-1}a^{n-1,*}-T^{n-1}\|_2^2\bigr|}
{\|F^{n-1}a^{n-1,*}-T^{n-1}\|_2^2}<\text{tol},
$$

全部实验取 $\text{tol}=10^{-5}$。

**自适应的三个部件。** 其一是**指标**：在 $\Omega_0$ 上取平均残差

$$
\mathcal L_{\Omega_0}(\alpha_0^*)=\frac{1}{|X_{f_0}|}\sum_{x\in X_{f_0}}
\bigl[L\tilde u_0(\alpha_0^*,x)-f(x)\bigr]^2 ,
$$

新峰中心取逐点残差的最大值点

$$
x_{K+1}=\arg\max_{x\in X_{f_0}}
\Bigl|L\Bigl(\sum_{m=0}^{M_0}\alpha^*_{m,0}\psi_{m,0}\Bigr)(x)-f(x)\Bigr| ,
$$

随后把区域切成 $\Omega^{*}=\Omega\setminus B_r(x_{K+1})$ 与 $\Omega^{**}=B_r(x_{K+1})\cap\Omega$。

其二是**局部基的重定心与重缩放**。在新建子域 $\Omega_K=B_{r_K}(x_K)\cap\Omega$ 上

$$
\psi_{0,K}(x)=1,
\qquad
\psi_{m,K}(x,c_K)=\sigma\bigl(c_K\,w_{m,K}^{\top}(x-x_K)+b_{m,K}\bigr),
\quad 1\le m\le M_K .
$$

两处改动各有分工：减去 $x_K$ 把基的关注点移到峰上，乘以 $c_K$ 压缩每个基函数的有效波长。这是多尺度网络思想的局部化版本。

其三是**尺度选择**。$c_K$ 由整数 $1,\dots,10$ 上的穷举搜索确定：对每个候选 $s$ 构造缩放基，解一个只涉及 $\Omega_K$ 上四条局部方程（内部残差、边界数据、以及与**冻结的**当前 $\tilde u_0$ 的两个界面匹配条件）的局部最小二乘问题 $\arg\min_{\alpha_K}\mathcal L_K(\alpha_K,s)$，记录 $\mathrm{Loss}_s$，取 $s^*=\arg\min_{1\le s\le10}\mathrm{Loss}_s$。穷举在这里是可行的，因为每次求解都是小规模线性最小二乘，而不是一次网络训练——这是冻结基路线的一个具体红利。论文指出也可以改用二分法。

最终逼近是分片展开，$u(x)\approx\tilde u_k(\alpha_k^*,x)$ 于 $x\in\Omega_k$。因此表示**在 $\Gamma_k$ 上是构造性间断的**，$C^1$ 匹配只在界面配点上弱地施加。

### 定理

**本文没有定理，没有收敛证明，也没有误差估计。** 这是一篇纯算法与数值的文章，结论部分自称是「攻克这类问题的第一次尝试」。把这一点说清楚是必要的：论文的说服力全部来自实验，因此下一节的细节就是它的全部证据。

### 数值实验

全部实验使用第二种预设方式（可迁移神经网络的重参数化）。误差度量为相对 $L^2$：

$$
\mathrm{err}_{L^2}=\sqrt{\sum_i|\hat u(x_i)-u(x_i)|^2\Big/\sum_i|u(x_i)|^2}.
$$

**算例一：二维 Poisson，$P$ 个高斯尖峰（第 4.2 节）。** $\Omega=[-1,1]^2$，精确解

$$
u=\sum_{p=1}^{P}\exp\bigl(-1000\bigl[(x-x_p)^2+(y-y_p)^2\bigr]\bigr),
$$

三个情形分别是一个峰在 $(0.5,0.5)$、两个峰在对角的 $(\pm0.5,\pm0.5)$、四个峰在 $(\pm0.5,\pm0.5)$。指数里的 $1000$ 是关键：它把峰宽压到 $O(0.03)$，比区域尺度小两个量级左右，正是全局统一尺度的基分辨不了的对象。

| 设定项           | 取值                              |
| ---------------- | --------------------------------- |
| 容差 $\epsilon$  | $10^{-4}$                         |
| 初始基数 $M_0$   | 200                               |
| 细化半径 $r$     | 0.15                              |
| 内部配点 $X_{f_0}$ | $50\times50$ 网格，$J_{f_0}=2500$ |
| 边界配点 $X_{g_0}$ | 400                               |
| 误差评估网格     | $256\times256$                    |

结果：三个情形分别收敛到 2、3、5 个子域，即**检测到的子域数恰好等于峰数加一个光滑主域**；检测到的中心是 $(\pm0.5102,\pm0.5102)$，与真实峰位相符，且事先没有提供任何峰位信息。误差随每峰子域基数 $M^*\in\{700,800,900,1000\}$ 单调下降，在 $M^*=1000$ 时 $\mathrm{err}_{L^2}\approx10^{-4}$，论文称其「显著小于」两篇自适应采样 PINN 工作的结果。每个峰子域选出的整数尺度都是 $c_k=5$。

**算例二：二维半线性尖峰问题（第 4.3 节）。** 方程改为 $-\Delta u+u^2=f$，三个峰情形与全部设定不变。这一例的作用是检验 Gauss-Newton 分支。结果：分区与算例一相同，选出的尺度同样是 $c_k=5$，而误差**比线性情形更小**。

**算例三：二维凹角奇性（第 4.4 节）。** L 形区域 $\Omega=[-1,1]^2\setminus[0,1]^2$，精确解 $u=(x^2+y^2)^{1/3}$。这是与尖峰不同的一类低正则性：奇性来自区域几何而不是源项。

| 设定项           | 取值        |
| ---------------- | ----------- |
| 容差 $\epsilon$  | $10^{-3}$   |
| 初始基数 $M_0$   | 600         |
| 细化半径 $r$     | 0.32        |
| 内部配点数       | $J_{f_0}=1875$ |
| 边界配点数       | $J_{g_0}=400$  |

结果：只需一次细化，$\Omega_1=B_r(x_1)$，$x_1=(-0.0345,-0.0345)$——非常接近凹角顶点。尺度在 $M^*=800,900$ 时为 $c_1=5$，在 $M^*=1000$ 时为 $c_1=4$；$M^*=1000$ 时 $\mathrm{err}_{L^2}\approx10^{-3}$。误差比尖峰算例大一个量级，而且尺度不随 $M^*$ 单调，这两点都提示**代数型奇性比高斯尖峰更难**：尖峰有一个明确的特征波长可以匹配，而 $r^{2/3}$ 型奇性在所有尺度上都不光滑。

**算例四：三维 Poisson，单个尖峰（第 4.5 节）。** $\Omega=[-1,1]^3$，$u=\exp(-1000[(x-0.5)^2+(y-0.5)^2+(z-0.5)^2])$。

| 设定项           | 取值                        |
| ---------------- | --------------------------- |
| 容差 $\epsilon$  | $10^{-4}$                   |
| 初始基数 $M_0$   | 2000                        |
| 细化半径 $r$     | 0.11                        |
| 内部配点数       | $J_{f_0}=10000$             |
| 边界配点数       | $J_{g_0}=2400$（每面 400）  |
| 子域内部配点     | 8500                        |
| 界面配点         | 600                         |

结果：检测中心 $(0.5170,0.5050,0.5000)$；尺度随基数增大，$M^*=3000,4000,5000$ 分别给出 $c_1=6,7,8$；$M^*=5000$ 时 $\mathrm{err}_{L^2}\approx10^{-3}$。

把四例放在一起看，可以读出方法的边界：**二维尖峰是最有利的情形（$10^{-4}$），而凹角奇性与三维尖峰都停在 $10^{-3}$**。三维那一例还需要 $M_0=2000$ 的初始基与 $5000$ 的子域基，比二维高出一个量级——冻结基方法在维数上的代价没有被消除，只是被推后了。

论文自陈的开放问题有三条：一次迭代找出所有峰子域（现在是一峰一轮）；确定缩放系数的更好办法（现在是整数穷举）；处理时间依赖的低正则性问题。

> [!note] 出版信息
> 期刊卷期与 DOI 已核实；**页码 553-577 未独立复核**。

### 与其他论文的关系

与编号 102 同属随机基一支：同样的冻结随机基、同样的线性最小二乘加 Gauss-Newton，区别在于所要修的东西——编号 90 调基的**尺度**以对付低正则性，编号 102 换**表述**以让约束精确成立。

缩放基 $\sigma(c_Kw^\top(x-x_K)+b)$ 是多尺度网络构造在随机基上的对应物，也是[[computational-mathematics/paper-notes/scientific-machine-learning/spectral-bias-and-generative-solvers|编号 101]] 多尺度随机 Fourier 特征库的近亲。三者是对同一个谱偏差问题的三种回答：编号 94 用多尺度网络，编号 90 在子域上离散搜索尺度，编号 101 在固定的多尺度库上用交叉注意力学习振幅。

残差最大值点驱动的细化，是同一批作者在别处所用深度自适应采样的「硬」版本：不是在残差大的地方**加样本**，而是在那里**加一组重缩放的局部基**。

与编号 60 的关系是同一个问题的两种答案：两者都针对普通罚项处理不好的约束，但编号 90 完全不用罚权而依赖最小二乘求解器，编号 60 则为边界条件引入乘子。

> [!note] 引用关系
> 编号 90 是否引用编号 60，本页依据的核对材料未能确认。

## 102：让散度约束由算子恒等式承担

### 直觉

不可压流的网络求解器几乎都把 $\operatorname{div}u=0$ 作为罚项。论文直接指出两个后果，而且这两个后果是叠加的。第一，罚参数的选择是关键问题，与编号 60 面对的困难同源。第二，无论用速度-压力还是速度-涡量表述，得到的都是**耦合** PDE 系统，几个未知场必须同时求解，最小二乘矩阵的列数随场数成倍增长，而最小二乘的代价对列数是平方的。

二维有一个经典出路：流函数。取 $u=\mathbf{curl}\,\psi$，无散条件因为 $\operatorname{div}\mathbf{curl}\equiv0$ 而恒等成立，一个标量未知量代替两个分量。但论文指出，把它推广到三维是主要困难——三维的势是向量，本身还需要一个规范条件。

本文的做法可以概括成一句一般原则：**对线性约束 $\mathcal A(u)=0$，找一个满足 $\mathcal A\circ\mathcal G=0$ 的算子 $\mathcal G$，然后用拟设 $u=\mathcal G(v)$。** 这里 $\mathcal A=\operatorname{div}$、$\mathcal G=\mathbf{curl}$。约束不再是需要执行的东西，而是恒等式。

真正的技术工作在第二步。取了 $\mathbf{curl}$ 拟设之后，方程要在势变量下重写，而对流项 $(u\cdot\nabla)u$ 在势变量下长什么样并不显然。论文的两条定理算出了它，三维那一条是本文新表述的关键。作为副产品，对整个动量方程取 $\operatorname{curl}$ 会消掉压力梯度（$\operatorname{curl}\nabla p=0$），于是**速度子问题里根本没有压力**，可以先解速度、再由一个独立的方程恢复压力。耦合就此解除。

### 问题设定

二维的两个旋度算子分别是

$$
\mathbf{curl}\,\psi:=\Bigl[\tfrac{\partial\psi}{\partial y},\,-\tfrac{\partial\psi}{\partial x}\Bigr]^{\mathrm T},
\qquad
\operatorname{curl}\,\bm v:=\tfrac{\partial v_2}{\partial x}-\tfrac{\partial v_1}{\partial y}.
$$

两条表示引理给出拟设的合法性。**Lemma 2.1（二维）：** 无散的 $\bm v$ 可写成 $\bm v=\mathbf{curl}\,\psi$，$\psi\in H^1(\Omega)$。**Lemma 2.2（三维，单连通）：** $\bm v=\mathbf{curl}\,\bm\psi$ 且 $\operatorname{div}\bm\psi=0$，可取 $\bm\psi\in H(\mathbf{curl};\Omega)$。三维的额外条件 $\operatorname{div}\bm\psi=0$ 就是规范固定，它在下面的推导里还要起第二个作用。

### 推导

**两条对流恒等式。**

**Theorem 2.1（二维）.** 对 $\phi\in C^3(\Omega)$、$\bm u=\mathbf{curl}\,\phi$，

$$
\operatorname{curl}\bigl((\bm u\cdot\nabla)\bm u\bigr)=-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi .
$$

**Theorem 2.2（三维）.** 对 $\bm\phi\in[C^3]^3$、$\bm u=\mathbf{curl}\,\bm\phi$ 且 $\operatorname{div}\bm\phi=0$，

$$
\mathbf{curl}\bigl((\bm u\cdot\nabla)\bm u\bigr)
=(\Delta\bm\phi\cdot\nabla)\mathbf{curl}\,\bm\phi-(\mathbf{curl}\,\bm\phi\cdot\nabla)\Delta\bm\phi .
$$

证明路径是两步。先用向量恒等式

$$
\mathbf{curl}\bigl((\bm u\cdot\nabla)\bm u\bigr)=(\bm u\cdot\nabla)\bm\omega-(\bm\omega\cdot\nabla)\bm u,
$$

再把涡量用势表示，

$$
\bm\omega=\mathbf{curl}\,\mathbf{curl}\,\bm\phi=\nabla(\nabla\cdot\bm\phi)-\Delta\bm\phi=-\Delta\bm\phi .
$$

第二个等号用掉了 $\operatorname{div}\bm\phi=0$。这说明规范条件在三维不是可有可无的技术附件：**没有它，涡量就不能干净地写成 $-\Delta\bm\phi$，整条恒等式随之失效。** 这是三维情形必须把无散条件加在势上的原因。

**四个解耦子问题。** 对二维 Stokes 方程 $-\nu\Delta\bm u+\nabla p=\bm f$ 取 $\operatorname{curl}$：压力项因 $\operatorname{curl}\nabla p=0$ 消失，粘性项因 $-\nu\Delta\operatorname{curl}\mathbf{curl}\,\phi=\nu\Delta^2\phi$ 变成双调和算子，得到单变量四阶问题

$$
\nu\Delta^{2}\phi=\operatorname{curl}\,\bm f\ \text{in}\ \Omega,
\qquad
\phi=\frac{\partial\phi}{\partial n}=0\ \text{on}\ \Gamma .
$$

夹紧边界条件不是假定，而是从无滑移推出来的：$\bm u\cdot\bm n_0=\nabla\phi\cdot\bm\tau_0=0$ 说明 $\phi$ 沿边界的切向导数为零，故 $\phi$ 在 $\Gamma$ 上为常数，可归一化为零；$\bm u\cdot\bm\tau_0=-\nabla\phi\cdot\bm n_0=0$ 给出法向导数为零。对滑移条件 $\bm u=\bm g$，两式改为 $\nabla\phi\cdot\bm\tau_0=\bm g\cdot\bm n_0$ 与 $-\nabla\phi\cdot\bm n_0=\bm g\cdot\bm\tau_0$。

二维 Navier-Stokes 把 Theorem 2.1 加进去：

$$
\nu\Delta^{2}\phi-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi=\operatorname{curl}\,\bm f,
\qquad
\phi=\frac{\partial\phi}{\partial n}=0\ \text{on}\ \Gamma .
$$

三维 Stokes 为

$$
\nu\Delta^{2}\bm\phi=\mathbf{curl}\,\bm f,
\qquad
\operatorname{div}\bm\phi=0\ \text{in}\ \Omega,
\qquad
\bm\phi\cdot\bm n=0,\ \ \mathbf{curl}\,\bm\phi=\bm 0\ \text{on}\ \Gamma
$$

（滑移条件把最后一条换成 $\mathbf{curl}\,\bm\phi=\bm g$），而三维 Navier-Stokes 是本文真正新的表述：

$$
\nu\Delta^{2}\bm\phi+(\Delta\bm\phi\cdot\nabla)\mathbf{curl}\,\bm\phi
-(\mathbf{curl}\,\bm\phi\cdot\nabla)\Delta\bm\phi=\mathbf{curl}\,\bm f,
$$

配同样的三条约束。论文的 Remark 3.3 把它与涡量-矢势的有限差分／有限元工作对比：后者需要**同时**逼近势与涡量，而这里只有单一未知量 $\bm\phi$，因此称之为纯流函数表述在三维的推广。

**压力恢复。** 关键在于它是一个**一阶梯度系统**而不是 Poisson 方程：

$$
\text{Stokes：}\ \nabla p=\bm f+\nu\Delta\bm u,
\qquad
\text{Navier-Stokes：}\ \nabla p=\bm f+\nu\Delta\bm u-(\bm u\cdot\nabla)\bm u,
\qquad p(\bm x_0)=0,
$$

其中 $\bm x_0$ 是任一内点，用于固定可加常数。论文指出这允许恢复**经典解**，而不是先前工作常用的弱解框架——取散度化成 Poisson 方程会丢掉一阶信息并引入额外的边界条件问题，直接配点一阶系统则不会。

### 定理

Theorem 2.1 与 Theorem 2.2 是本文仅有的两条定理，而且它们是**代数恒等式，不是收敛结果**。论文没有误差分析，也没有关于方法的收敛定理。这一点与编号 90 一致：随机基一支的说服力来自数值证据。

### 数值实验

**基与离散。** 基与编号 90 完全相同（可迁移神经网络的三步重参数化，冻结），向量值情形下**各分量共享同一组隐层基** $\{\psi_m\}$，只有输出系数 $\bm\alpha^{(i)}\in\mathbb R^{M+1}$ 独立——这正是下面块矩阵对角线上重复出现 $\nu\Delta^2\bm\Psi_{in}$ 的原因。

以二维 Stokes 速度步为例，在 $\bm X^{in}=\{x^{in}_i\}_{i=1}^{I}\subset\Omega$ 与 $\bm X^{bd}=\{x^{bd}_j\}_{j=1}^{J}\subset\Gamma$ 上配置三组方程

$$
\nu\Delta^{2}\!\sum_{m=0}^{M}\!\alpha_m\psi_m(x^{in}_i)=\operatorname{curl}\,\bm f(x^{in}_i),
\qquad
\sum_m\alpha_m\psi_m(x^{bd}_j)=0,
\qquad
\frac{\partial}{\partial n}\sum_m\alpha_m\psi_m(x^{bd}_j)=0,
$$

得到超定线性系统，用 NumPy 的线性求解器按最小二乘解出。压力步解

$$
\min_{\bm\beta}\Bigl\|
\begin{bmatrix}\partial_x\bm\Psi_{in}\\ \partial_y\bm\Psi_{in}\\ \bm\Psi_0^{\mathrm T}\end{bmatrix}\bm\beta
-\begin{bmatrix}\tilde{\bm F}_1\\ \tilde{\bm F}_2\\ 0\end{bmatrix}\Bigr\|_2^2,
$$

最后一行就是 $p(\bm x_0)=0$。

Navier-Stokes 用 Gauss-Newton 线性化

$$
(\bm u\cdot\nabla)\bm u\ \longrightarrow\
\bm u^{(k)}\!\cdot\!\nabla\bm u^{(k+1)}+\bm u^{(k+1)}\!\cdot\!\nabla\bm u^{(k)}
-\bm u^{(k)}\!\cdot\!\nabla\bm u^{(k)},
$$

二维下每步的线性 PDE 为

$$
\nu\Delta^{2}\phi^{(k+1)}-(\mathbf{curl}\,\phi^{(k)}\!\cdot\!\nabla)\Delta\phi^{(k+1)}
-(\mathbf{curl}\,\phi^{(k+1)}\!\cdot\!\nabla)\Delta\phi^{(k)}
=\operatorname{curl}\,\bm f-(\mathbf{curl}\,\phi^{(k)}\!\cdot\!\nabla)\Delta\phi^{(k)} .
$$

论文称之为 Newton 加线性最小二乘框架，并指出它等价于先离散出非线性系统再用 Gauss-Newton。附录给出的 Picard/Oseen 型线性化（「Scheme I」，$\nu\Delta^2\phi^{(k+1)}-\mathbf{curl}\,\phi^{(k)}\cdot\nabla\Delta\phi^{(k+1)}=\mathbf{curl}\,\bm f$）用作初始化，因为 Gauss-Newton 只是局部收敛的。

> [!warning] 预印本中的一处符号
> 论文式 (4.11) 在线性化后的方程里印成 $\mu$，而上下文一贯使用黏性系数 $\nu$。按语境应为 $\nu$。

**代价核算。** 设 $\min_x\|Ax-b\|_2^2$、$A\in\mathbb R^{m\times n}$、$m\gg n$ 的代价为 $\mathcal O(mn^2)$，则

| 维数 | 解耦无散网络                                            | 耦合基线                      |
| ---- | ------------------------------------------------------- | ----------------------------- |
| 二维 | $\mathcal O((I+2J)(M+1)^2)+\mathcal O((2I+1)(M+1)^2)$   | $\mathcal O((3I+2J)(3M+3)^2)$ |
| 三维 | $\mathcal O((4I+4J)(3M+3)^2)+\mathcal O((3I+1)(M+1)^2)$ | $\mathcal O((4I+3J)(4M+4)^2)$ |

节省来自 $(M+1)^2$ 与 $(3M+3)^2$ 或 $(4M+4)^2$ 的差别：顺序解两个小系统，而不是一次解含 3 到 4 个耦合场的大系统。因为最小二乘代价对列数是平方的，把 $3M$ 换成 $M$ 就是九倍。

**误差度量与基线。** 相对 $L^2$ 误差与绝对散度误差分别为

$$
\mathrm{error}\_g=\frac{\sqrt{\tfrac1N\sum_i[g-g_{\rm NN}]^2}}{\sqrt{\tfrac1N\sum_i g^2}},
\qquad
\mathrm{error}\_{\rm div}=\sqrt{\tfrac1N\sum_i\bigl[\operatorname{div}\bm g_{\rm NN}(x_i)\bigr]^2}.
$$

基线有两个：标准的**耦合** TransNet，以及一个 PINN（ResNet，宽 $M=30$、深 $L=4$、$\tanh$、Adam、512 内部点加 256 边界点、10000 次迭代）。硬件为 Intel Core i9-12900H 笔记本、16 GB 内存、RTX 3060；Python 3.12 / PyTorch 2.1。

**算例一：二维 Stokes，Kovasznay 型解（第 5.1 节）。** $\Omega=(0,2)\times(-0.5,1.5)$，

$$
u=1-e^{\zeta x}\cos(2\pi y),
\quad
v=\tfrac{\zeta}{2\pi}e^{\zeta x}\sin(2\pi y),
\quad
p=\tfrac12 e^{2\zeta x},
\quad
\zeta=\tfrac{1}{2\nu}-\sqrt{\tfrac{1}{4\nu^2}+4\pi^2}.
$$

训练用 $50\times50=2500$ 内部点加 200 边界点的均匀网格，测试用 $111\times111=12321$ 点。

| 对照                       | 结果                                                                 |
| -------------------------- | -------------------------------------------------------------------- |
| 对耦合 TransNet，$\nu=10^{-4}$ | 速度与压力均更优；散度误差至少 $\mathcal O(10^{-12})$              |
| 扫描 $\nu$，$M=1000$       | 大 $\nu$ 时 TransNet 略好；$\nu\le10^{-2}$ 时本方法优近两个量级       |
| 运行时间                   | 约 2 秒，对 TransNet 的略超 5 秒                                      |
| 对 PINN                    | 速度与压力误差至少小四个量级；PINN 散度误差仅 $\mathcal O(10^{-3})$，且在低黏性下完全失去精度 |

低黏性下的分野是这一例最有信息量的部分：$\nu$ 越小，Kovasznay 解的指数层越陡，耦合表述的条件数越差，而解耦表述把速度问题变成一个纯双调和问题，不受压力耦合牵连。

**算例二：二维 Navier-Stokes（第 5.2 节）。** 精确解

$$
u=16y(y-1)(2y-1)\sin^2(\pi x),
\quad
v=-8\pi y^2(y-1)^2\sin(2\pi x),
\quad
p=\sin(\pi x)\cos(\pi y),
$$

网格与算例一相同；Gauss-Newton 的初值取 PINN 训练 10000 次迭代的结果，最多 40 次 Gauss-Newton 迭代。$M=1000$ 时的结果是：$\nu\le10^{-2}$ 下精度与 TransNet **相当**——这是预期之中的，因为两者用同一个初值和同一类 Newton 求解器——但代价低 60% 以上，且散度保持严格更优。这一例因此不是精度声明，而是**代价声明加约束声明**。

**算例三：三维 Stokes（第 5.3 节）。** $(0,1)^3$ 上，

$$
u=e^{\cos(\pi y)}\sin(\pi z),
\quad
v=e^{\cos(\pi z)}\sin(\pi x),
\quad
w=e^{\cos(\pi x)}\sin(\pi y),
$$

$$
p=e^{\cos(\pi x)+\sin(\pi y)}+e^{\cos(\pi z)+\sin(\pi x)} .
$$

采样用 Halton 序列：10000 内部点加 2400 边界点（$20\times20\times6$），2000 个 Halton 测试点。

| 对照           | 结果                                                     |
| -------------- | -------------------------------------------------------- |
| 速度           | 在所有 $M$ 上都比 TransNet 更准                          |
| 压力           | 大 $M$ 时 TransNet 约好一个量级，但耗时接近两倍          |
| 散度           | $10^{-14}$，比 TransNet 好六个量级                       |
| 扫描 $\nu$     | 速度在 $\nu\le10^{-4}$、压力在 $\nu\le10^{-3}$ 时本方法更优 |
| 代价           | 约为 TransNet 的一半                                     |

压力那一行值得注意，因为它是论文自己报告的一处劣势：**顺序求解意味着压力步继承了速度步的误差**，而耦合求解让两者互相校正。论文没有回避这一点。

> [!note] 预印本内部不一致
> 图 5.7 的说明写 $M=1500$，而正文写 $M=2500$；图 5.9 的说明误重复为「3D Stokes」。

**算例四：三维 Navier-Stokes（第 5.4 节）。** $(0,1)^3$ 上的三角多项式解

$$
u=2(y-1)(z-1)\sin y\sin z-2(y-1)\sin y\cos z-2(z-1)\cos y\sin z
$$

（$v,w$ 由轮换给出），$p=xyz+x^3y^3z-\tfrac{5}{32}$。初值取 Scheme I 迭代 10 次的结果，最多 15 次 Gauss-Newton 迭代。结果：$\nu<10^{-3}$ 时速度与压力精度更高，加上一贯的散度与代价优势。

**总体判断。** 论文自陈：散度误差接近机器精度（$10^{-14}$）；在基函数数量较大时执行时间比 TransNet 一类框架少约 50%，二维 Navier-Stokes 情形超过 60%；低黏性／高 Reynolds 数下速度精度占优。它同时诚实地指出一条限制：**Stokes 的表现普遍好于 Navier-Stokes，原因归结为 Gauss-Newton 初值的质量。** 这与算例二、四都要靠外部初始化（PINN 或 Scheme I）是一致的——解耦解决了耦合与约束的问题，但没有解决非线性的全局收敛问题。

代码「可向通讯作者合理请求获取」，无公开仓库。

### 与其他论文的关系

编号 102 是编号 90 的直系兄弟：同一个学派、同一套冻结随机基（逐字相同的三步重参数化）、同一套线性最小二乘加 Gauss-Newton 求解。编号 90 调基的尺度以对付低正则性，编号 102 换表述以让约束精确成立；两者合起来构成本主题的随机基一支。

「用 $u=\mathcal G(v)$ 硬性消去约束」这一取向，与编号 60 为本质边界条件引入乘子的路线正好相反：一个是结构性消元，一个是乘子执行。哪一条可行取决于约束是否有现成的势表示——散度约束有 $\mathbf{curl}$，Dirichlet 条件在复杂区域上没有同样干净的对应物。

「把耦合系统拆成一列更便宜的子问题」这个动作，与[[computational-mathematics/paper-notes/scientific-machine-learning/spectral-bias-and-generative-solvers|编号 103]] 在时间方向做算子分裂是同一个思路，只是拆分的方向不同：这里拆的是场之间的耦合，那里拆的是算子之间的耦合。

## 三篇的共同判断

| 编号 | 被移出罚项的约束   | 承担者                       | 训练变成什么           |
| ---- | ------------------ | ---------------------------- | ---------------------- |
| 60   | 本质边界条件       | 乘子网络 + 增广 Lagrange     | 交替的两个随机优化问题 |
| 90   | 界面匹配与局部尺度 | 区域分解 + 局部重缩放基      | 一列线性最小二乘       |
| 102  | 无散条件           | $\mathbf{curl}$ 拟设与恒等式 | 顺序的两个最小二乘     |

三条都在回答同一个问题：**当约束以罚项形式进入损失时，权重成为一个必须调的自由参数；把约束交给结构，这个自由参数就消失了。** 代价各不相同：编号 60 多一个网络与一次投影，编号 90 多一层区域分解与一次尺度搜索，编号 102 把二阶方程升成四阶方程。这些代价是可预测的，而罚权的调节不是。

第二条判断关于求解器：编号 90 与 102 之所以能完全取消罚权，是因为它们把训练换成了最小二乘。最小二乘对各残差块的相对权重仍然敏感，但它不需要一个能让梯度下降收敛的权重——这两件事的难度不在一个量级上。编号 60 保留了随机优化，因此仍需要 $\beta$，只是让乘子把 $\beta$ 的量级要求降了下来。

第三条判断关于证据：**这三篇都没有关于其离散方法的收敛定理。** 编号 60 引用的是无穷维层面的经典理论，编号 90 与 102 明说没有分析。因此三篇的说服力都落在实验设计上，而其中编号 90 与 102 的实验报告得足够细致，可以被复核。

## 覆盖核对

| 内容                            | 论文 | 覆盖状态                                    |
| ------------------------------- | ---- | ------------------------------------------- |
| 两种通行做法的缺陷与乘子的作用  | 60   | 直觉、拟设法与罚法各自的问题                |
| 增广 Lagrange 与极小极大等价性  | 60   | 函数形式、符号约定、内层极大化的推导、对偶函数 |
| 三类实例的具体能量              | 60   | 椭圆、线性特征值、非线性特征值及其唯一性引用 |
| 乘子更新的参数空间投影          | 60   | 困难来源、投影目标、Monte Carlo 形式与括号说明 |
| 60 的数值实验                   | 60   | 算法 2 与无穷维算法 1、三类算例、两个基线、Remark 3.1；**误差表与倍数未核实** |
| 冻结基与两种预设方式            | 90   | 基的定义、均匀抽样、重参数化与位置/形状分离 |
| 区域分解与无权最小二乘          | 90   | 传输条件、损失、罚权可取消的原因、Gauss-Newton 与停机准则 |
| 残差指标、局部重缩放与尺度搜索  | 90   | 三个部件、穷举可行的理由、构造性间断        |
| 90 的数值实验                   | 90   | 四个算例的完整设定与结果，含检测中心、尺度与误差量级、三条开放问题 |
| 无散拟设与两条对流恒等式        | 102  | 一般原则、两条引理、两条定理及证明路径、规范条件的作用 |
| 四个解耦子问题与压力恢复        | 102  | 四阶速度方程、边界条件的推导、一阶梯度系统与经典解 |
| 离散化、Gauss-Newton 与代价核算 | 102  | 共享基、两步最小二乘、线性化与初始化、复杂度对照 |
| 102 的数值实验                  | 102  | 四个算例的解、网格、基线、时间与散度量级，以及自陈的两处劣势 |

## 本页原文

- J. Huang, H. Wang, and T. Zhou, [_An augmented Lagrangian deep learning method for variational problems with essential boundary conditions_](https://doi.org/10.4208/cicp.OA-2021-0176), Commun. Comput. Phys. 31(3) (2022), pp. 966-986（预印本 [arXiv:2106.14348](https://arxiv.org/abs/2106.14348)）。
- J. Huang, H. Wu, and T. Zhou, [_Adaptive neural network basis methods for partial differential equations with low-regular solutions_](https://doi.org/10.4208/cicp.OA-2024-0310), Commun. Comput. Phys. 39(2) (2026), pp. 553-577（预印本 [arXiv:2411.01998](https://arxiv.org/abs/2411.01998)）。
- J. Cheng, J. Huang, H. Wang, and T. Zhou, _Decoupled divergence-free neural networks basis method for incompressible fluid problems_, [arXiv:2603.17906](https://arxiv.org/abs/2603.17906)，投稿 Commun. Comput. Phys.
