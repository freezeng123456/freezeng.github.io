---
title: 多精度代理与在线细化
description: 编号 34、37、49、79：一个闭环的四次改装，含完整推导、定理假设与数值实验
lang: zh
translation: en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates
tags:
  - 论文笔记
  - 贝叶斯反问题
  - 代理模型
---

编号 **37**、**34**、**49**、**79** 是同一个闭环的四次改装，每篇替换其中一个部件。

## 四篇共享的骨架

把四篇放在一起看，它们回答的是同一组三个问题，只是每篇换一个答案：

1. **用什么代理**——低阶多项式混沌加高精度修正（37、34）、复合神经网络（49）、算子网络（79）。
2. **什么时候承认代理不够用**——在一个精心挑选的单点上比较代理与真模型（37、34、49），或者比较拟合本身有没有变坏（79）。
3. **触发之后在哪里补点**——局部球内随机取点（37、49）、球内取点（34）、在代理输出空间中贪心选点（79）。

循环本身是固定的：用当前代理跑一段采样，在一个点上算指标，指标超阈值就花几次真求解重训代理，回到采样。**代价控制全部集中在第二步**：指标每轮只允许一次（或几次）真求解，因此它必须是单点量，而单点量要有意义，这个点就必须落在后验主体附近。四篇的差别归根结底是「怎么找到那个点」。

## 37：闭环的原型

### 直觉

先验多项式混沌代理是这样构造的：从先验抽一批参数，在每个参数上求解真问题，然后把输入到输出的映射拟合成一个多项式。这个代理在先验抽到样本的地方准确——问题在于后验不在那儿。观测数据一旦有信息量，后验就会收缩到先验支撑的一个薄子集上；数据信息越强，这个子集越薄，甚至可能落到先验的尾部。于是代理的精度分布和后验的质量分布**几乎正交**：算力花在了后验永远不会访问的区域。

有两条常规出路，论文都拒绝了。一是提高多项式阶数直到全局精确——总次数空间的维数 $\binom{N+n_z}{n_z}$ 随参数维数组合式增长，中等维数就不可行。二是换一个更靠近后验的先验——但后验位置事先未知，这是循环论证。

论文的出路是第三条：**保持低阶，但让代理跟着链走**。链自己会走到后验主体，那就在链走到的地方补一次真求解，把代理在那一带修正好。这样代理永远不是全局准确的，它只在链关心的地方准确——而下面的 Theorem 2 恰好说明这就够了。

还剩一个问题：链是用代理跑的，它走到的地方是**代理的后验**主体，不一定是真后验主体。如果完全信任代理的接受概率，链就会在代理自己的高概率区里打转，指标在那里当然合格，细化永远不被触发，误差被自己确认成了没有误差。论文的补丁是在每一轮外循环里插入一步用**真模型接受概率**的 Metropolis 判定，用它产生的被接受点作为检验点。这一步的成本是每轮两次真求解（在 $z^-$ 和 $z^+$ 上各一次），换来的是检验点由真后验而不是代理后验挑选。**这是「代理放进 MCMC 需要接受步纠正」这句话在本族里的确切含义**——注意它纠正的是细化点的位置，不是样本本身，后一点见下文。

### 问题设定

参数 $z\in\Gamma\subset\mathbb R^{n_z}$，先验密度 $\pi(z)$，观测数据 $d\in\mathbb R^{n_d}$，正算子记 $u^{H}:\Gamma\to\mathbb R^{n_d}$（高精度，即真模型），代理记 $u^{L}$（低精度）。观测噪声为独立同分布高斯。真后验记 $\pi^{d}$，代理诱导的后验记 $\widetilde\pi^{d}_{N}$。

低精度模型是先验多项式混沌展开：以 $\{\Phi_m\}_{m=1}^{M}$ 记关于先验测度正交的多项式基，$M=\dim\Lambda_N$ 为总次数指标集 $\Lambda_N$ 的大小。

### 推导

**第一步：系数的加权离散最小二乘。** 在 $Q$ 个样本点 $\{z_i\}$ 上，

$$
\{c_m\}_{m=1}^{M}=\arg\min_{c}\sum_{i=1}^{Q}
\Bigl[G(z_i)-\sum_{m=1}^{M}c_m\Phi_m(z_i)\Bigr]^{2},
\qquad
w_i=\frac{M}{\sum_{m=1}^{M}\Phi_m^{2}(z_i)} .
$$

写成矩阵形式：记 $\Psi\in\mathbb R^{Q\times M}$、$\Psi_{im}=\Phi_m(z_i)$、$W=\mathrm{diag}(w_i)$，则求解的是 $c^{\#}=\arg\min_c\|\sqrt W\Psi c-\sqrt W b\|_2^2$。

这个权不是随手取的。记再生核对角线 $K(z)=\sum_{m=1}^{M}\Phi_m^{2}(z)$，则 $w_i=M/K(z_i)$，于是 $\sqrt W\Psi$ 第 $i$ 行的平方范数为

$$
\sum_{m=1}^{M}w_i\Phi_m^{2}(z_i)=\frac{M}{K(z_i)}\cdot K(z_i)=M,
$$

**与 $z_i$ 落在哪里无关**。也就是说加权把设计矩阵的每一行都归一化到同一长度 $\sqrt M$，从而 $\frac1Q\Psi^{T}W\Psi$ 的迹恒等于 $M$；采样密度剩下的唯一作用，是决定这份固定的迹如何分摊到各个特征值上。

$w_i$ 正是离散形式的 Christoffel 函数，与[[computational-mathematics/paper-notes/stochastic-approximation/optimal-sampling-and-preconditioning|最优采样与预条件]]一页研究的是同一个对象：那里的稳定性因子是 $\|K\|_\infty/N$，即上式中权的倒数取上确界，而诱导采样 $q\propto\sum_n v_n^{2}$ 就是让这个上确界降到它的下界。采样密度按度数渐近设计：均匀测度用张量积 Chebyshev，高斯测度用支撑在半径 $\sqrt{2N}$ 球上的密度（Narcowich–Ward 型设计）。

**第二步：修正项只在低阶指标上展开。** 定义两个精度之差

$$
C(z)=u^{H}(z)-u^{L}(z)\ \approx\ \sum_{\mathbf m\in\Lambda_{N_C}}u^{C}_{\mathbf m}\Phi_{\mathbf m}(z),
\qquad N_C\ll N .
$$

**第三步：合并只是一次重新分组。** 把 $u^{L}=\sum_{\alpha\in\Lambda^{d}_{N}}u^{L}_{\alpha}\Psi_{\alpha}$ 与 $C$ 相加，因为 $\Lambda^{d}_{N_C}\subset\Lambda^{d}_{N}$，按指标归并即得

$$
f_{M}(\theta)=\sum_{\alpha\in\Lambda^{d}_{N_C}}\bigl(u^{L}_{\alpha}+u^{C}_{\alpha}\bigr)\Psi_{\alpha}
+\sum_{\alpha\in\Lambda^{d}_{N}\setminus\Lambda^{d}_{N_C}}u^{L}_{\alpha}\Psi_{\alpha}.
$$

这一步没有近似，只是把和式重排。它的价值在于把结构显示出来：**只有低阶系数被高精度数据校正，高阶系数保留先验拟合的值**。这不是疏漏而是刻意的——每轮只有 $Q=2\binom{N_C+n_z}{n_z}$ 个新的真求解，即修正空间维数的两倍，恰好是加权最小二乘的常用过采样量；把这点样本摊到高阶指标上必然过拟合。

**第四步：触发准则。** 细化在被接受点上判定，用的是**绝对** $\ell^\infty$ 误差：

$$
\mathrm{err}(y)=\bigl\|u^{H}(y)-u^{L}(y)\bigr\|_{\infty}.
$$

这里的 $y$ 不是任意链状态，而是用高精度接受概率接受的点，理由已在直觉一节说明。

**第五步：局部细化设计。** 若 $\mathrm{err}(y)>\epsilon$，在 $\ell^\infty$ 球 $B(y,R)=\{x:\|x-y\|_{\infty}\le R\}$ 内取 $Q$ 个随机点重拟合。半径按输入常数 $\rho$ 逐轮收缩，理由是早期被接受点可能仍离后验主体较远：初期需要大球以覆盖不确定的位置，后期需要小球以换取局部精度。

**算法的整体节奏。** 外循环每轮先用当前低精度模型跑 $m-1$ 步标准 Metropolis-Hastings 子链（子链的作用是让末状态与起点去相关），再提议 $z^{*}\sim q(\cdot\mid z_{m-1})$，用真模型接受概率在 $z^{-}=z_{m-1}$ 与 $z^{+}=z^{*}$ 之间判定得到 $y$，算指标，必要时细化，最多 $I_{\max}$ 轮。**最终的后验样本是所有子链样本的并集。**

> [!warning] 这个算法采的是谁的后验
> 子链是用当前低精度模型跑的，因此并集里的样本服从的是**代理诱导的后验** $\widetilde\pi^{d}_{N}$，不是真后验 $\pi^{d}$。真模型只出现在两个位置：挑选检验点，以及细化时的 $Q$ 次求值。所以这套算法的正确性依据不是 Metropolis 修正把样本纠回真后验（那是 Conrad–Marzouk–Pillai–Smith 的渐近精确性路线），而是下面的 Theorem 2：只要 $\widetilde\pi^{d}_{N}$ 与 $\pi^{d}$ 在 KL 意义下足够近，用前者代替后者就够了。两条路线的保证强度不同，不应混为一谈。

### 定理

**Assumption 1.** 正算子一致有界：$\sup_{z\in\Gamma}\|u^{H}(z)\|=:C_H<\infty$。

定义 $\epsilon$-可行集及其后验测度

$$
\Gamma_{N}(\epsilon)=\bigl\{y\in\Gamma:\ \|u^{H}(y)-u^{L}(y)\|_{\infty}\le\epsilon\bigr\},
\qquad
\mu\bigl(\Gamma_{N}(\epsilon)\bigr)=\int_{\Gamma_{N}(\epsilon)}\pi^{d}(z)\,dz,
$$

补集记 $\Gamma^{\perp}_{N}(\epsilon)=\Gamma\setminus\Gamma_{N}(\epsilon)$。

**Theorem 2（KL 上界）.** 假设 $u^{H}$ 与 $u^{L}$ 关于 $N$ 一致地满足 Assumption 1，观测噪声独立同分布高斯，则对给定的 $\epsilon>0$ 存在常数 $K_1,K_2>0$ 使

$$
D_{\mathrm{KL}}\bigl(\widetilde\pi^{d}_{N}\,\|\,\pi^{d}\bigr)
\le\Bigl(K_1\epsilon+K_2\,\mu\bigl(\Gamma^{\perp}_{N}(\epsilon)\bigr)\Bigr)^{2}.
$$

**Corollary 3（Hellinger 上界，论文归于其参考文献 [8]）.** 同样假设下

$$
D_{\mathrm{Hell}}\bigl(\widetilde\pi^{d}_{N}\,\|\,\pi^{d}\bigr)
\le K_1\epsilon+K_2\,\mu\bigl(\Gamma^{\perp}_{N}(\epsilon)\bigr),
$$

即去掉平方。

**推论（论文式 21，直接给出算法目标）.** 若采样充分好，使得 $\mu(\Gamma^{\perp}_{N}(\epsilon))\le\epsilon$，则 KL 距离完全由 $\epsilon^{2}$ 刻画。

界的形状本身说明了它的来源：两项对应把 KL 积分按 $\Gamma_{N}(\epsilon)$ 与其补集拆开——在可行集上势函数之差被 $\epsilon$ 控制，给出 $K_1\epsilon$；在补集上代理可以任意坏，只能用 Assumption 1 的一致界压住被积函数，再乘以那块区域的后验质量，给出 $K_2\mu(\Gamma^{\perp}_{N}(\epsilon))$。这解释了为什么需要一个看上去多余的一致有界假设。论文本身只叙述结论，没有给出这个拆解。

> [!warning] 已证明的与只是论证的
> 论文对算法与定理之间的联系是明确表述为**机制**而非定理的：只要候选点落在 $\Gamma^{\perp}_{N}(\epsilon)$ 内，算法就在其邻域细化，因此 $\mu(\Gamma^{\perp}_{N}(\epsilon))$ 随细化渐近衰减。**这是一个论证，不是收敛速率，也没有给出 $\mu(\Gamma^{\perp}_{N}(\epsilon))$ 何时或以多快速度趋于零。** 这一族的论文在这一点上处境相同。

### 数值实验

两个非线性 PDE 反问题：

| 例  | 问题                 | 参数维数 | 设计意图                                                     |
| --- | -------------------- | -------- | ------------------------------------------------------------ |
| 1   | 二维热源反演         | $n_z=2$  | 维数小到可以负担真正精确的高阶先验代理，因此这里买到的是精度 |
| 2   | 椭圆方程扩散系数反演 | $n_z=9$  | 全局精确的先验代理很贵，因此这里同时买到精度与代价           |

两个例子的分工是刻意的：例 1 排除了「自适应方法只是因为对手太弱才赢」这一解释，因为在 $n_z=2$ 上对手可以做到全局精确；例 2 才是方法的目标场景。论文报告的量级是相对于纯真模型 MCMC 有若干数量级的效率提升。

**实验建立了什么，又差在哪里。** 建立的是：低阶代理加局部修正能在两类非线性 PDE 反问题上恢复精度。差的是三点。第一，两个例子的参数维数是 2 和 9，而论文自己诉诸的动机是高维下阶数爆炸，这个动机没有在实验里被直接检验（编号 49 正是为此而写）。第二，Theorem 2 的界含未定常数 $K_1,K_2$，实验没有也无法据此校准 $\epsilon$，容差是当作调参给的。第三，「$\mu(\Gamma^{\perp}_{N}(\epsilon))$ 随细化衰减」这个机制没有被直接测量——实验测的是后验边缘分布是否对上，不是那块坏集的测度是否在缩小。

### 与其他论文的关系

这是整族的起点。编号 34 是集合 Kalman 版的移植，编号 49 是神经网络版的移植（并明确批评本篇的两个弱点：多项式对低正则性映射逼近差，且受维数诅咒），编号 [[computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|55 与 56]] 把「训练点放在后验附近」这一思想分别移植到 randomize-then-optimize 与 Stein 变分梯度下降。本篇的 Theorem 2 就是编号 79 的 Theorem 2.1（两篇引同一来源），也是[[computational-mathematics/paper-notes/bayesian-inference/index|专题首页]]中那条 $\epsilon$-可行集界。

## 34：同一个修正装进集合 Kalman 反演

### 直觉

集合 Kalman 反演不需要正算子的导数，这在 PDE 约束问题里是决定性的优点：伴随方程可能没有实现，或者求解器根本不可微。代价是它用一个粒子集合的经验协方差代替导数信息，而经验协方差要可信，集合就不能小；$N_e$ 个成员乘 $J$ 轮迭代就是约 $N_eJ$ 次正求解。

把正模型换成先验多项式混沌代理，在线代价几乎归零。但集合 Kalman 反演有一个和 MCMC 不同的失效方式，本篇的实验设计正是冲着它来的：**集合会主动往数据指的方向走，如果数据指向先验之外，集合就会走出代理的可靠区域，而代理不会报警，它只会继续给出光滑的、错误的输出**。实验里真值从 $\log\theta_i\sim U(-4,4)$ 抽取而先验是 $\log\theta_i\sim\mathcal N(0,1)$，就是把这个失效方式放到最大。

检验点的选择在这里比在 MCMC 里更自然。MCMC 需要额外一步真模型接受判定才能得到一个可信的检验点；集合 Kalman 反演直接就有一个现成的候选——**集合均值**，它就是当前集合对真值位置的最佳猜测。这省掉了编号 37 里那次额外的接受判定。

### 问题设定

正问题与观测（论文式 2.1）：

$$
y=f(\theta)+\xi,\qquad \xi\sim\mathcal N(0,\Gamma),
\qquad \theta\in\mathbb R^{d},\ y\in\mathbb R^{m}.
$$

集合规模 $N_e$，第 $n$ 轮第 $j$ 个成员记 $\theta^{(j)}_n$，其输出 $\omega^{(j)}_n=f(\theta^{(j)}_n)$，扰动数据 $y^{(j)}=y+\xi^{(j)}$，$\xi^{(j)}\sim\mathcal N(0,\Gamma)$。

### 推导

**第一步：更新式从哪里来。** Iglesias–Law–Stuart 的集合 Kalman 反演在**扩张状态空间** $Z=X\times Y$ 上工作，$z=(u,p)$、$p=\mathcal G(u)$，投影 $H=[0,I]$。预测步把每个粒子按人工动力学推进 $\hat z^{(j)}_{n+1}=(u^{(j)}_n,\mathcal G(u^{(j)}_n))$，分析步用 Kalman 增益

$$
K_n=C_nH^{T}\bigl(HC_nH^{T}+\Gamma\bigr)^{-1},
\qquad
z^{(j)}_{n+1}=\hat z^{(j)}_{n+1}+K_{n+1}\bigl(y^{(j)}_{n+1}-H\hat z^{(j)}_{n+1}\bigr).
$$

把 $C_n$ 按 $z=(u,p)$ 分块，代入 $H=[0,I]$ 并用 $H^{\perp}=[I,0]$ 投影出参数分量，$HC_nH^{T}$ 就是 $C^{pp}_n$、$C_nH^{T}$ 的参数块就是 $C^{up}_n$，于是

$$
u^{(j)}_{n+1}=u^{(j)}_{n}+C^{up}_{n+1}\bigl(C^{pp}_{n+1}+\Gamma\bigr)^{-1}
\bigl(y^{(j)}_{n+1}-\mathcal G(u^{(j)}_{n})\bigr).
$$

（这一步归约是标准的；Iglesias 等人的原文写的是分块前的形式。）编号 34 采用的正是这个式子，只在 $\Gamma$ 前多一个正则化参数 $\alpha_n$：

$$
\theta^{(j)}_{n+1}=\theta^{(j)}_{n}
+C^{\theta\omega}_{n}\bigl(C^{\omega\omega}_{n}+\alpha_{n}\Gamma\bigr)^{-1}\bigl(y^{(j)}-\omega^{(j)}_{n}\bigr),
\qquad j=1,\dots,N_e,
$$

经验协方差为

$$
C^{\theta\omega}_{n}=\frac{1}{N_e-1}\sum_{j=1}^{N_e}\bigl(\theta^{(j)}_{n}-\bar\theta_{n}\bigr)\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)^{\!T},
\qquad
C^{\omega\omega}_{n}=\frac{1}{N_e-1}\sum_{j=1}^{N_e}\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)^{\!T}.
$$

**第二步：为什么这叫「迭代正则化」。** 三件事一起起作用。其一是结构性的：集合始终停留在初始集合的线性张成 $\mathcal A=\mathrm{span}\{\psi^{(j)}\}_{j=1}^{J}$ 之内，因此在 $X$ 上病态的残差极小化，实际上是在紧集 $\mathcal A$ 上做的。其二是可比对的：在线性情形 $\mathcal G(u)=Gu$ 下，$J\to\infty$ 时一步集合更新收敛到 Tikhonov–Phillips 解

$$
u_{TP}=\bar u+CG^{*}\bigl(GCG^{*}+\Gamma\bigr)^{-1}(y-G\bar u),
$$

而这**恰好是线性高斯后验均值**（Stuart 的 (3.4) 在有限维数据下的形式）。也就是说集合 Kalman 反演在线性情形下不是近似而是精确，代价只是它不需要 $\mathcal G$ 的导数。其三是算法性的：迭代必须按偏差原理停机。

**第三步：正则化参数与停机。** $\alpha_n$ 取满足

$$
\alpha^{N}_{n}\bigl\|\Gamma^{1/2}\bigl(C^{\omega\omega}_{n}+\alpha^{N}_{n}\Gamma\bigr)^{-1}\bigl(y^{(j)}-\bar\omega_{n}\bigr)\bigr\|
\ \ge\ \rho\,\bigl\|\Gamma^{-1/2}\bigl(y^{(j)}-\bar\omega_{n}\bigr)\bigr\|
$$

的第一个整数 $N$ 所对应的值；停机用偏差原理

$$
\bigl\|\Gamma^{-1/2}(y-\bar\omega_{n})\bigr\|\le\tau\eta,
\qquad \rho<1,\ \tau\ge1/\rho .
$$

两条都来自 Iglesias 的正则化集合 Kalman 平滑器。**Iglesias 等人自己明确说明：完整的收敛与正则化分析超出该文范围，偏差原理是数值上支持的，不是被证明的。** 因此这条停机准则在本族里始终只是经验规则。Iglesias 等人还记录了一个反直觉的经验事实：把扰动数据换成不扰动（$\eta^{(j)}_{n+1}=0$）结果反而**更差**，作者猜测噪声帮助算法探索 $\mathcal A$。

**第四步：指标改成相对形式，检验点改成集合均值。**

$$
\mathrm{err}=\frac{\bigl\|f(\bar\theta_{n+1})-f_{M}(\bar\theta_{n+1})\bigr\|_{\infty}}{\bigl\|f(\bar\theta_{n+1})\bigr\|_{\infty}},
\qquad
\bar\theta_{n+1}=\frac{1}{N_e}\sum_{j}\theta^{(j)}_{n+1}.
$$

若 $\mathrm{err}\le\mathrm{tol}$ 就继续用当前代理，否则按多精度最小二乘重建 $f_M$。样本预算是显式的：初始化用 $Q_1=2\binom{N+d}{d}$ 个先验样本，每次细化用 $Q_2=2\binom{N_C+d}{d}$ 个新点。

### 定理

**本篇不给出自适应格式的收敛定理，全部论断都是计算性的。** 上一节引用的线性情形 Tikhonov–Phillips 极限与偏差原理都来自 Iglesias 等人，不是本篇的结果。

### 数值实验

正问题是二维**时间分数阶**反扩散问题（Caputo 导数，阶 $0<\alpha<1$）。

| 项目        | 设置                                |
| ----------- | ----------------------------------- |
| 时间离散    | 有限差分，$\Delta t=0.01$           |
| 空间离散    | 谱方法，多项式次数 $P=6$            |
| 数据生成    | $P=10$，以避免反演犯罪              |
| 例 1 未知量 | 渗透率中的 9 个径向基权重（$d=9$）  |
| 例 2 未知量 | Karhunen-Loève 参数化的对数渗透率场 |

**例 1：真值刻意落在先验之外。** 先验取 $\log\theta_i\sim\mathcal N(0,1)$，真值从 $\log\theta_i\sim U(-4,4)$ 抽取；$N_e=100$，噪声水平 $\sigma=10^{-3}$。

| 方法                           | 多项式阶 | 相对误差 | 真模型求值数  | CPU 时间   |
| ------------------------------ | -------- | -------- | ------------- | ---------- |
| 常规集合 Kalman 反演           | —        | 0.0461   | 2000（在线）  | 约 56.71 s |
| 固定先验多项式混沌代理         | $N=4$    | 0.7921   | 1430（离线）  | 约 0.82 s  |
| 固定先验多项式混沌代理         | $N=6$    | 0.2892   | 10010（离线） | 未报告     |
| 自适应多精度（容差 $10^{-2}$） | $N=2$    | 未报告   | 250           | 未报告     |
| 自适应多精度（容差 $10^{-3}$） | $N=2$    | 未报告   | 575           | 未报告     |

几点值得单独指出。

- 离线求值数与预算公式**对得上**：$d=9$ 时 $Q_1=2\binom{N+d}{d}$ 在 $N=4$ 给 $2\binom{13}{9}=1430$，在 $N=6$ 给 $2\binom{15}{9}=10010$，与表中数字完全一致。这反过来确认了参数维数确实是 9，也确认了预算公式在实现中是照字面执行的。
- 常规方法的 2000 次求值在 $N_e=100$ 下对应 20 轮迭代，与 $N_eJ$ 的代价模型一致。
- 自适应版本的相对误差论文没有给出数值，只说局部修正基本恢复了全模型集合 Kalman 反演的精度。
- 求解次数的对比是清楚的：$250$ 对 $2000$ 是八分之一，$575$ 对 $2000$ 是不到三分之一，且都远低于把固定代理提到 $N=6$ 所需的 10010。

**例 2：高维随机场。** 保留 95% 先验能量对应 $d=22$ 个 Karhunen-Loève 模态，$N_e=300$，噪声 $\mathcal N(0,0.01^2)$。

| 方法                   | 多项式阶 | 相对误差 |
| ---------------------- | -------- | -------- |
| 固定先验多项式混沌代理 | $N=2$    | 0.3430   |
| 固定先验多项式混沌代理 | $N=3$    | 0.2146   |
| 自适应多精度           | $N=1$    | 0.0889   |

**实验建立了什么，又差在哪里。** 建立的是一个相当锋利的结论：**当真值出先验时，提高全局代理的阶数是低效的修补方式**。从 $N=4$ 到 $N=6$ 把离线求值数从 1430 提到 10010（七倍），误差只从 0.7921 降到 0.2892，仍然比常规方法差六倍；而自适应方法用 250 次求解、基于阶数只有 2 的展开，就够了。例 2 更极端：自适应方法在**一阶**展开上做到 0.0889，比二阶和三阶的固定代理都好一个量级。

差的地方有四处。第一，自适应版本的 CPU 时间没有报告，因此「便宜」这个结论建立在真求解次数上，而不是实测墙钟时间上——多精度重拟合本身的开销没有被计入对照。第二，自适应版本的相对误差没有数值，只有定性判断，因此「基本恢复全模型精度」无法量化。第三，两个例子都在同一个 PDE 族（时间分数阶扩散）里，跨族的普适性未被检验。第四，本篇没有定理，出先验这一情形恰恰是 Theorem 2 那类界最难用的场合——界里的 $\mu(\Gamma^{\perp}_{N}(\epsilon))$ 在真值出先验时初值很大，而它如何衰减仍然只是论证。

### 与其他论文的关系

这是编号 37 的集合 Kalman 版：同一个多精度修正、同一族 $\ell^\infty$ 触发准则，驱动者从 MCMC 链换成集合迭代，指标从绝对改成相对，检验点从被接受点换成集合均值。编号 49 把代理换成复合神经网络并保留触发准则，编号 79 再换一次代理（DeepONet）、再换一次采样器（无迹 Kalman 反演）、并把指标换成数据残差。编号 [[computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|82]] 把「廉价模型加学到的修正装进 Kalman 方法」这一思路搬到时间演化的数据同化。

## 49：把多项式换成复合神经网络

### 直觉

编号 37 的循环没有问题，问题在于它的代理。多项式混沌有两个此处致命的弱点：对正则性有限的映射逼近能力差，以及基函数个数随参数维数组合式爆炸。两者都直接卡住向高维的推广。

换成神经网络之后，出现一个新的具体问题：**每轮只买得起十次左右真求解，用十个样本训练一个网络必然过拟合。** 论文的解法不是把网络做小到无用，而是改变要学的对象——不学 $z\mapsto f^{H}(z)$，而学「已有代理错在哪里」。而且不是像编号 37 那样学一个可加的修正项，而是把已训练的低精度网络当作**输入变量**送进一个新的浅网络。

这个区别值得说清楚。可加修正 $f^{H}\approx f^{L}+f_{\mathrm{CORR}}$ 假设两个精度之差本身是一个容易拟合的函数；复合形式 $f^{H}\approx\mathcal{NN}(z,\mathcal{NN}^{L}(z))$ 不作这个假设，它允许高精度输出以任意方式依赖低精度输出，例如乘性缩放或者分段的行为切换。因为两个模型高度相关，这个映射本身很简单——**这正是浅网络够用的原因，而浅网络够用又正是十个样本不至于过拟合的原因**。容量约束在这里是链条中的一环，不是事后的正则化。

### 问题设定

参数 $z$，数据 $d$，似然 $\mathcal L(d,\cdot)$，先验 $\pi$。低精度模型 $f^{L}=\mathcal{NN}^{L}$ 是按先验训练的深网络（实验用 4 隐层 × 40 神经元，50 个训练点），高精度模型 $f^{H}$ 是真求解。

### 推导

**第一步：复合多精度网络。** 设两个精度之间存在一个未知的非线性关系

$$
f^{H}(z)=\mathcal F\bigl(z,f^{L}(z)\bigr)=\mathcal F\bigl(z,\mathcal{NN}^{L}(z)\bigr),
$$

并用一个网络逼近这个关系：

$$
f^{H}(z)\approx\mathcal{NN}^{H}(z;\theta):=\mathcal{NN}\bigl(z,\mathcal{NN}^{L}(z);\theta\bigr).
$$

论文自己给出与编号 37 的对照：多项式混沌的修正 $f^{H}\approx f^{L}_{\mathrm{PCE}}+f_{\mathrm{CORR}}$ 是**线性**叠加，上式学的是两个精度之间的**非线性相关**。

**第二步：训练集的构造。** 选 $Q$ 个点 $\{z_k\}_{k=1}^{Q}$，在每个点上跑一次 $\mathcal{NN}^{L}$（便宜）与一次 $f^{H}$（贵），组成

$$
\mathcal D=\Bigl\{\bigl((z_k,\mathcal{NN}^{L}(z_k)),\ f^{H}(z_k)\bigr)\Bigr\}_{k=1}^{Q},
$$

在 $\mathcal D$ 上训练 $\mathcal{NN}^{H}(z;\theta)$。注意输入是**二元组**，低精度预测值作为一个额外的输入通道进入网络——这就是复合的实现方式。

**第三步：容量的约束链。** $Q$ 必须小（真求解贵，实验取 $Q=10$）$\Rightarrow$ 修正网络必须浅（实验用单隐层 50 神经元，论文表述为至多两个隐层）$\Rightarrow$ 只有在两个精度高度相关、待学映射本身简单时这才可行。**论文明确把浅网络说成防过拟合的约束，而不是性能优化。**

**第四步：指标与细化。**

$$
\mathrm{err}(\tilde z)=\frac{\bigl\|f^{H}(\tilde z)-f^{L}(\tilde z)\bigr\|_{\infty}}{\bigl\|f^{H}(\tilde z)\bigr\|_{\infty}} .
$$

相比编号 37 的绝对误差改成相对误差。若 $\mathrm{err}(\tilde z)>\mathrm{tol}$，在 $B(\tilde z,R)=\{z:\|z-\tilde z\|_\infty\le R\}$ 内均匀取 $Q$ 点重训 $\mathcal{NN}^{H}$，然后令 $f^{L}\leftarrow\mathcal{NN}^{H}$。

**第五步：复合关系逐轮嵌套。** 上一步的赋值有一个不易察觉的后果：下一轮的「低精度模型」是上一轮的复合网络，因此第 $k$ 轮的网络把第 $k-1$ 轮的网络当输入，而后者又把第 $k-2$ 轮的当输入。复合深度随外循环增长，尽管每一层修正都是浅的。论文没有讨论这个嵌套的稳定性。

**外循环。** 输入初始先验训练代理 $f^{L}=\mathcal{NN}^{L}$、提议密度 $q$、子链长度 $m$（实验取 $m=1000$）、最大修正次数 $I_{\max}$。每轮跑 $m-1$ 步子链、提议 $z^{*}$、按需细化、算接受概率并接受或拒绝，最终返回汇集的后验样本。

> [!warning] 印出的接受概率分子分母是反的
> 论文式 (10) 把高精度接受概率印成
> $\beta=\min\{1,\ \mathcal L(d,f^{H}(z^{-}))\pi(z^{-})/\mathcal L(d,f^{H}(z^{+}))\pi(z^{+})\}$。
> 在对称提议的 Metropolis-Hastings 中，被提议状态 $z^{+}$ 应在分子上。按印出的形式，$\beta$ 会是被提议状态后验密度的**减函数**：数据拟合得越好的提议越不容易被接受，链会朝低概率区漂移——这与论文自己给该点设定的用途（让检验点落在后验主体附近）直接矛盾，所以印出的形式不可能是本意。按上下文这里应为
> $\beta=\min\{1,\ \mathcal L(d,f^{H}(z^{+}))\pi(z^{+})/\mathcal L(d,f^{H}(z^{-}))\pi(z^{-})\}$。

### 定理

**本篇不证明任何定理**，并明确说明沿用编号 37 的分析框架。设计论断只有一条，而它是前提性的：$Q$ 必须小，因此修正网络的容量必须受限——这是防过拟合的约束。

### 数值实验

基准椭圆方程反问题的两个配置。

| 项目         | 例 1                                                     |
| ------------ | -------------------------------------------------------- |
| 未知量       | 9 参数渗透率                                             |
| 先验代理     | 4 隐层 × 40 神经元，50 个训练点                          |
| 修正网络     | 单隐层 × 50 神经元                                       |
| 每轮细化点数 | $Q=10$                                                   |
| 容差         | $\mathrm{tol}\in\{0.1,\ 0.05\}$                          |
| 正则化系数   | $\lambda=0$                                              |
| 子链长度     | $m=1000$                                                 |
| 基线         | Direct（真模型 MCMC）、DNN（固定先验网络）、ADNN（本文） |

例 2 换成由 Karhunen-Loève 展开参数化的高维随机场渗透率。

**结果。** 固定的先验训练网络给出**可见错误**的后验边缘分布；自适应版本恢复真模型 MCMC 的边缘分布；容差从 $0.1$ 收紧到 $0.05$ 后一致性进一步提高。

**实验建立了什么，又差在哪里。** 建立的是消融意义上的结论：同一条链、同一个初始代理，唯一的差别是在线细化开还是关，开着的那条能对上真模型 MCMC 的边缘分布，关着的那条对不上。这直接把改进归因到细化机制，而不是网络架构或训练技巧。差的地方在于：论文的动机是高维（多项式的基函数爆炸），而例 1 只有 9 个参数；没有任何量化指标，因此无法比较不同容差、不同 $Q$ 之间的取舍；复合嵌套随外循环加深，其影响完全未被测量。

### 与其他论文的关系

编号 37 的直接后继：同一个外层 MCMC 骨架、同一套「在被接受点周围的 $\ell^\infty$ 球内重训」设计，把多项式混沌换成复合神经网络，把绝对误差指标换成相对误差指标。复合思想（把上一个代理当作下一个的输入）是本篇的标志性贡献，编号 79 与 [[computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|106]] 都不再采用，它们改用迁移学习微调同一个网络。论文自己指出向编号 34（集合 Kalman）与 [[computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|55]]（randomize-then-optimize）的推广「也是可能的」。

## 79：目标导向指标与贪心设计

### 直觉

前三篇的指标都在问同一个问题：代理在这一点上准不准。这个问题有一个隐蔽的缺陷——**代理不准不一定要紧**。反演真正在乎的是数据能不能被拟合上；代理在某个方向上有误差，如果这个方向对数据残差没有影响，那么修正它是浪费。反过来，代理在一个方向上误差很小，但如果这个方向恰好主导残差，那么这点误差就要紧。前三篇通过「把检验点选在后验主体附近」来间接对齐这两件事，而本篇直接换掉被测量的对象：**测数据残差本身**。

这一换带来两个后果。好的一面：指标现在和反演目标同单位、同尺度，阈值可以按残差的相对变化设定，不必猜代理误差的量级。代价的一面：残差是一个标量，它不告诉你代理在哪个方向上错了，因此选点规则必须自己补上这个信息——本篇的贪心设计就是干这个的，它在**代理输出空间**里挑彼此差异最大的点，也就是代理自己认为最不一样、因而最可能出错的地方。

第二处直觉是关于采样器的。无迹 Kalman 反演每轮只要 $2N_m+1$ 次正求值，通常十来轮收敛，因此整个反演的正求解总量是几百次量级——比 MCMC 少两三个数量级。这使得「每轮花一次真求解算指标」的相对开销变得可观，也就更需要指标本身廉价。贪心打分只用代理在候选池上的预测值，因此选点这一步几乎不花代价，真求解只花在最终选中的 $Q$ 个点上。

### 问题设定

无穷维 Bayes 表述：未知 $m\in\mathcal M$，正算子 $\mathcal G$，观测 $y=\mathcal G(m)+\eta$，$\eta\sim\mathcal N(0,\Sigma_\eta)$，势函数 $\Phi(m;y)=\tfrac12\|y-\mathcal G(m)\|^{2}_{\Sigma_\eta}$。代理 $\widehat{\mathcal G}_t=\mathcal O\circ\mathcal F_\theta$ 由 DeepONet $\mathcal F_\theta$ 与观测算子 $\mathcal O$ 复合而成，$\nu_t$ 是第 $t$ 轮代理诱导的后验近似。

### 推导

**第一步：先写出诚实的指标，再拒绝它。** 局部模型误差应当是

$$
e_{M}(t):=\mathbb E_{\nu_t}\bigl\|\mathcal G-\widehat{\mathcal G}_t\bigr\|
=\Bigl(\int_{\mathcal M}\bigl|\mathcal G(m)-\widehat{\mathcal G}_t(m)\bigr|^{2}\,\nu_t(dm)\Bigr)^{1/2}.
$$

论文明确拒绝它作为**可实现**的指标：它需要一个高维积分。这一步定下了后面所有单点指标的地位：它们是这个量的廉价替身，不是它的估计。

**第二步：锚点。** 从当前代理后验 $\nu_t$ 抽 $T$ 个样本 $\mathcal M^{(t)}=\{m^{(t)}_k\}_{k=1}^{T}$，用**真模型**在其中挑数据拟合最好的一个：

$$
r_t=\arg\min_{m\in\mathcal M^{(t)}}\ \tfrac12\bigl\|y-\mathcal G(m)\bigr\|^{2}_{\Sigma_\eta}.
$$

**第三步：指标与触发。** 指标就是锚点处的数据残差，触发用它的**相对变化**：

$$
e_{D}(t):=\Phi(r_t;y)=\tfrac12\bigl\|y-\mathcal G(r_t)\bigr\|^{2}_{\Sigma_\eta},
\qquad
\frac{\bigl|e_D(t)-e_D(t-1)\bigr|}{e_D(t)}\ >\ \epsilon .
$$

不满足就停止细化。用相对变化而不是绝对值有一层意思：$e_D$ 本身的下界由噪声水平决定，不可能降到零，因此以「还在改善吗」而不是「够小了吗」作为判据。实验取 $\epsilon=0.01$、$I_{\max}=10$。

**第四步：贪心设计。** 从 $\nu_t$ 抽一个大候选池 $\Gamma=\{m_1,\dots,m_K\}$，逐点生长子集 $\gamma_Q=\{\hat m_1,\dots,\hat m_Q\}\subset\Gamma$：

$$
\hat m_{j+1}=\arg\max_{m\in\Gamma\setminus\gamma_j}
\Bigl\{d\bigl(\widehat{\mathcal G}_t(m),\widehat{\mathcal G}^{j}_t\bigr)-\lambda\|m-r_t\|_{2}\Bigr\},
\qquad
d\bigl(\widehat{\mathcal G}_t(\cdot),\widehat{\mathcal G}^{j}_t\bigr)
=\max_{\hat m\in\gamma_j}\bigl\|\widehat{\mathcal G}_t(\cdot)-\widehat{\mathcal G}_t(\hat m)\bigr\|_{2},
$$

其中 $\widehat{\mathcal G}^{j}_t:=\{\widehat{\mathcal G}_t(\hat m_i)\}_{i=1}^{j}$，全部实验取 $\lambda=1$。两项刻意对抗：第一项要求在**代理输出空间**中彼此远离，服务于泛化；第二项把选点拉回锚点附近，服务于局部精度。打分只用代理预测值，因此选点相对于随后的全阶求解几乎免费。

**第五步：重训用迁移学习。** DeepONet 参数从上一轮的权重出发，不从头训练。

**第六步：无迹 Kalman 反演。** 后验近似由随机动力系统给出：

$$
m_{n+1}=r_0+\alpha(m_n-r_0)+\omega_{n+1},\quad \omega_{n+1}\sim\mathcal N(0,\Sigma_\omega),
\qquad
y_{n+1}=\mathcal G(m_{n+1})+\eta_{n+1},\quad\eta_{n+1}\sim\mathcal N(0,\Sigma_\eta),
$$

$\alpha\in(0,1]$ 是正则化参数。预测步

$$
\hat r_{n+1}=\alpha r_n+(1-\alpha)r_0,
\qquad
\hat C_{n+1}=\alpha^{2}C_n+\Sigma_\omega ,
$$

分析步用高斯联合分布给出

$$
r_{n+1}=\hat r_{n+1}+\hat C^{my}_{n+1}\bigl(\hat C^{yy}_{n+1}\bigr)^{-1}\bigl(y_{n+1}-\hat y_{n+1}\bigr),
\qquad
C_{n+1}=\hat C_{n+1}-\hat C^{my}_{n+1}\bigl(\hat C^{yy}_{n+1}\bigr)^{-1}\bigl(\hat C^{my}_{n+1}\bigr)^{T},
$$

其中 $\hat C^{yy}_{n+1}=\mathrm{Cov}[\mathcal G(m_{n+1})\mid Y_n]+\Sigma_\eta$。期望用修正无迹变换在 $2N_m+1$ 个对称 $\sigma$ 点上求值：

$$
m^{0}=r,\qquad m^{j}=r+c_j[\sqrt C]_j,\qquad m^{j+N_m}=r-c_j[\sqrt C]_j,\qquad 1\le j\le N_m .
$$

每轮代价是 $2N_m+1$ 次前向求值，论文报告通常在 $O(10)$ 轮收敛。

**闭环。** 离线用少量先验样本训练 $\mathcal F_\theta$；用当前代理跑无迹 Kalman 反演得到 $\nu_t$ 与 $T$ 个样本（这一步只碰代理与数据）；算锚点与 $e_D(t)$ 并检验相对变化；若触发则抽候选池、贪心选 $\gamma_Q$、在其上跑全阶模型、从当前权重微调；重复至多 $I_{\max}$ 轮。

### 定理

**Theorem 2.1（引自论文参考文献 [43]）.** 与编号 37 的 Theorem 2 同形。记 $\nu$ 为真后验、$\widehat\nu$ 为代理诱导后验、$\mathcal M(\epsilon)$ 为 $\epsilon$-可行集、$\mathcal M^{\perp}(\epsilon)$ 为其补集，则存在 $K_1,K_2>0$ 使

$$
D_{\mathrm{KL}}(\widehat\nu\,\|\,\nu)\ \le\ \bigl(K_1\epsilon+K_2\,\nu(\mathcal M^{\perp}(\epsilon))\bigr)^{2}.
$$

这是「只要求局部精度」的理论依据，本篇把它作为动机引用。

**Assumptions 3.2–3.4.**

- (3.2) 对任意 $\epsilon$，**线性**神经算子 $\widehat G:\mathbb R^{N_m}\to\mathbb R^{N_y}$ 可训练到 $\|\widehat G-G\|_2<\epsilon$；
- (3.3) $\|G\|_2<H$；
- (3.4) $G^{T}\Sigma_\eta^{-1}G\succ0$ 且 $\|G^{T}\Sigma_\eta^{-1}G\|_2>C_1$。

**Lemma 3.5.** 在 3.2–3.4 下 $\|\widehat G^{T}\Sigma_\eta^{-1}\widehat G\|_2>C_2$，途径是

$$
\bigl\|G^{T}\Sigma_\eta^{-1}G-\widehat G^{T}\Sigma_\eta^{-1}\widehat G\bigr\|_2\le 2\epsilon H\bigl\|\Sigma_\eta^{-1}\bigr\|_2 .
$$

这一步靠的是一个初等恒等式：

$$
A^{T}\Sigma^{-1}A-B^{T}\Sigma^{-1}B=A^{T}\Sigma^{-1}(A-B)+(A-B)^{T}\Sigma^{-1}B,
$$

取 $A=G$、$B=\widehat G$，两项各以 $H\|\Sigma_\eta^{-1}\|_2\epsilon$ 为界即得。随后由反向三角不等式得到 $C_2$ 的存在性。（论文只叙述结论并抽象地断言某个 $C_2>0$ 存在，这一步的恒等式是常规补充。）

**Theorem 3.6（线性情形的收敛）.** 假设：3.2–3.4，$\mathrm{Range}(G^{T})=\mathrm{Range}(\widehat G^{T})=\mathbb R^{N_m}$，$\Sigma_\omega\succ0$，$\Sigma_\eta\succ0$。结论：代理驱动的无迹 Kalman 反演不动点 $(\widehat r_\infty,\widehat C_\infty^{-1})$ 收敛到全模型不动点 $(r_\infty,C_\infty^{-1})$，且

$$
\bigl\|\widehat C^{-1}_\infty-C^{-1}_\infty\bigr\|_2\le\frac{2\epsilon H H_\eta}{1-\beta},
\qquad
\bigl\|\widehat r_\infty-r_\infty\bigr\|_2
\le\frac{K_1H_\eta H_y}{C_1}\Bigl(1+\frac{2(1+\alpha\beta)K_2H_\eta H^{2}}{(1-\beta)C_2}\Bigr)\epsilon ,
$$

其中 $\beta,C_1,C_2,K_1,K_2,H_\eta,H_y,H$ 都是正的有界常数。两个界都是 $O(\epsilon)$，即关于代理的算子范数误差线性。

**这是本族里唯一一条关于采样器输出的收敛结果。** 编号 37 与本篇的 Theorem 2.1 界的是后验之间的距离，而这一条界的是算法实际算出的均值与协方差。

> [!warning] 定理适用范围的限定
> 论文的 Remark 1 说明，为满足 Theorem 3.6 可以让 $\widehat G_\theta$ 线性化，做法是去掉分支网络中的非线性激活、保留主干网络的激活。也就是说，这条定理覆盖的是一个受限的 DeepONet，不是非线性实验中实际使用的网络。引用该定理时应带上这一限定。

### 数值实验

| 项目           | 设置                                                                               |
| -------------- | ---------------------------------------------------------------------------------- |
| 基准问题       | Darcy 流、热源反演（Case I / Case II）、反应扩散                                   |
| 代理           | DeepONet，分支与主干均为 5 隐层 × 100 神经元，`tanh` 激活                          |
| 离线训练       | $1\times10^{5}$ 步，$N_{\mathrm{prior}}=1000$ 个高斯随机场先验样本                 |
| 观测数据       | $y_{\mathrm{obs}}=y_{\mathrm{ref}}+\max\{\lvert y_{\mathrm{ref}}\rvert\}\delta\xi$ |
| 指标阈值       | $\epsilon=0.01$                                                                    |
| 重训上限       | $I_{\max}=10$                                                                      |
| 贪心权重       | $\lambda=1$                                                                        |
| 每轮采样器代价 | $2N_m+1$ 次前向求值，通常 $O(10)$ 轮                                               |
| 基线           | FEM-UKI（全模型）、DeepOnet-UKI-Direct（固定代理）、DeepOnet-UKI-Adaptive（本文）  |

热源反演 Case I 是两参数源位置问题：DeepONet 在 $[0.5,1]\times[0.5,1]$ 上用 500 个均匀样本训练，无迹 Kalman 反演从 $[0.6,0.6]$ 出发；Case II 是更高维的变体。

**关键的实验设计是分布内与分布外两组真值。** 分布外的真值就是把「代理在后验所在处不可靠」这一失效模式显式制造出来，也正是固定代理失效的场合。

**结果。** 固定代理只给出粗略估计，并产生可见错误的反演轨迹；自适应版本随细化单调降低模型误差，并接近全阶模型精度。第三个例子中细化在六轮后由停机准则终止——这是停机准则确实会触发而不只是被 $I_{\max}$ 截断的证据。

**实验建立了什么，又差在哪里。** 建立的是三点：目标导向指标可以自行停机；贪心选点在三个不同类型的 PDE 上都有效；分布外真值下固定算子代理确实失败，而自适应版本恢复精度。差的是：全部实验用的是非线性 DeepONet，而 Theorem 3.6 只覆盖线性化的分支网络，因此**定理与实验之间存在一条本篇自己指出的缝**；模型误差单调下降是定性描述的，没有配套的数值序列；无迹 Kalman 反演给出的是高斯近似，因此对多峰后验的适用性未被检验（这正是编号 106 接手的地方）。

### 与其他论文的关系

本族里理论最完整的一篇：既把编号 37 的 $\epsilon$-可行集 KL 界作为动机引用，又补上了一条真正的（线性情形的）采样器收敛定理，而编号 34、37、49 与 [[computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|55、56]] 都没有。它把编号 34、37、49 的逐点 $\ell^\infty$ 指标换成由数据残差本身构造的目标导向指标，把球内随机取点换成多样性与邻近性的显式权衡。它也是 [[computational-mathematics/paper-notes/bayesian-inference/sampling-and-filtering|106]] 的直接方法论前身：后者保留自适应微调闭环，把无迹 Kalman 反演换成变分流、DeepONet 换成 Fourier 神经算子、贪心筛选换成激进的数据替换。

## 四篇的对照表

| 部件     | 37                     | 34                     | 49                     | 79                        |
| -------- | ---------------------- | ---------------------- | ---------------------- | ------------------------- |
| 代理     | 多精度多项式混沌       | 多精度多项式混沌       | 复合多精度网络         | DeepONet                  |
| 采样器   | Metropolis-Hastings    | 正则化集合 Kalman 反演 | Metropolis-Hastings    | 无迹 Kalman 反演          |
| 指标     | 绝对 $\ell^\infty$     | 相对 $\ell^\infty$     | 相对 $\ell^\infty$     | 数据残差的相对变化        |
| 检验点   | 高精度接受概率接受的点 | 集合均值               | 高精度接受概率接受的点 | 真模型选出的锚点          |
| 新点选择 | 球内随机，半径收缩     | 球内随机               | 球内均匀               | 代理输出空间贪心 + 邻近项 |
| 重训方式 | 重新最小二乘           | 重新最小二乘           | 训练浅修正网络         | 从上轮权重迁移学习        |
| 理论     | KL / Hellinger 上界    | 无定理                 | 无定理                 | 线性情形的不动点收敛      |

## 本页原文

- L. Yan and T. Zhou, [_Adaptive multi-fidelity polynomial chaos approach to Bayesian inference in inverse problems_](https://doi.org/10.1016/j.jcp.2018.12.025), J. Comput. Phys. 381 (2019), pp. 110-128（预印本 [arXiv:1807.00618](https://arxiv.org/abs/1807.00618)）。
- L. Yan and T. Zhou, [_An adaptive multifidelity PC-based ensemble Kalman inversion for inverse problems_](https://doi.org/10.1615/Int.J.UncertaintyQuantification.2019029059), Int. J. Uncertain. Quantif. 9(3) (2019), pp. 205-220（预印本 [arXiv:1809.08931](https://arxiv.org/abs/1809.08931)）。
- L. Yan and T. Zhou, [_An adaptive surrogate modeling based on deep neural networks for large-scale Bayesian inverse problems_](https://doi.org/10.4208/cicp.OA-2020-0186), Commun. Comput. Phys. 28 (2020), pp. 2180-2205（预印本 [arXiv:1911.08926](https://arxiv.org/abs/1911.08926)）。
- Z. Gao, L. Yan, and T. Zhou, [_Adaptive operator learning for infinite-dimensional Bayesian inverse problems_](https://doi.org/10.1137/24M1643815), SIAM/ASA J. Uncertain. Quantif. 12(4) (2024), pp. 1389-1423（预印本 [arXiv:2310.17844](https://arxiv.org/abs/2310.17844)）。
- 背景来源：M. A. Iglesias, K. J. H. Law, and A. M. Stuart, [_Ensemble Kalman methods for inverse problems_](https://doi.org/10.1088/0266-5611/29/4/045001), Inverse Problems 29 (2013), 045001；A. M. Stuart, [_Inverse problems: a Bayesian perspective_](https://doi.org/10.1017/S0962492910000061), Acta Numerica 19 (2010), pp. 451-559。
