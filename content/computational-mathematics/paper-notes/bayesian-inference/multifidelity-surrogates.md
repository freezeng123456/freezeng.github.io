---
title: 多精度代理与在线细化
description: 编号 34、37、49、79：同一个闭环在四种代理与四种采样器上的实现
lang: zh
translation: en/computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates
tags:
  - 论文笔记
  - 贝叶斯反问题
  - 代理模型
---

> [!note] 本页覆盖
> 编号 **37**（_J. Comput. Phys._ 381, 2019）、**34**（_Int. J. Uncertain. Quantif._ 9(3), 2019）、**49**（_Commun. Comput. Phys._ 28, 2020）、**79**（_SIAM/ASA J. Uncertain. Quantif._ 12(4), 2024）。四篇共享同一个闭环骨架，逐篇替换其中一个部件。

## 37：闭环的原型

### 问题定位

先验多项式混沌代理让 MCMC 变便宜，但它在先验支撑上拟合，而后验集中在先验支撑的薄子集上。提高阶数可以消除误差，代价却随参数维数组合式增长。论文的选择是**保持低阶，在线修正**。

### 多精度修正的具体形式

设低精度模型是先验多项式混沌展开 $u^L$，高精度模型 $u^H$ 就是真正问题。修正项

$$
C(z)=u^{H}(z)-u^{L}(z)\approx\sum_{\mathbf m\in\Lambda_{N_C}}u^{C}_{\mathbf m}\Phi_{\mathbf m}(z)
$$

只在一个**低阶**指标集 $\Lambda_{N_C}$ 上展开，$N_C\ll N$。合并后的多精度模型把两组系数按指标相加：

$$
f_{M}(\theta)=\sum_{\alpha\in\Lambda^{d}_{N_C}}\bigl(u^{L}_{\alpha}+u^{C}_{\alpha}\bigr)\Psi_{\alpha}
+\sum_{\alpha\in\Lambda^{d}_{N}\setminus\Lambda^{d}_{N_C}}u^{L}_{\alpha}\Psi_{\alpha}.
$$

这个写法说明了修正的结构性质：**只有低阶系数被高精度数据校正**，高阶系数保留先验拟合的值。这一点是刻意的——高阶指标上的样本量不足，强行校正会过拟合。

系数由加权离散最小二乘确定：

$$
\{c_m\}_{m=1}^{M}=\arg\min_{c}\sum_{i=1}^{Q}\Bigl[G(z_i)-\sum_{m=1}^{M}c_m\Phi_m(z_i)\Bigr]^{2},
\qquad
w_i=\frac{M}{\sum_{m=1}^{M}\Phi_m^{2}(z_i)} .
$$

这里的权 $w_i$ 正是 Christoffel 函数的离散形式，与[[computational-mathematics/paper-notes/stochastic-approximation/index|随机逼近专题]]中的最优采样理论同源；采样密度按度数渐近设计，均匀测度用张量 Chebyshev，高斯测度用支撑在半径 $\sqrt{2N}$ 球上的密度。

### 触发准则与局部细化

细化在被接受点上判定，且用的是**绝对** $\ell^\infty$ 误差：

$$
\mathrm{err}(y)=\bigl\|u^{H}(y)-u^{L}(y)\bigr\|_{\infty}.
$$

这里的 $y$ 不是任意链状态，而是用**高精度接受概率**接受的点。这一细节是整个设计的关键：用真模型接受一个点，代价是一次真求解，换来的是被检验的点更接近后验主体。

若 $\mathrm{err}(y)>\epsilon$，在 $\ell^\infty$ 球 $B(y,R)=\{x:\|x-y\|_\infty\le R\}$ 内取 $Q$ 个随机点重拟合，总次数空间取 $Q=2\binom{N_C+n_z}{n_z}$。半径按一个输入常数 $\rho$ 逐轮收缩，理由是早期被接受点可能仍离后验主体较远。

### Theorem 2：为什么局部精度够用

定义 $\epsilon$-可行集及其后验测度

$$
\Gamma_{N}(\epsilon)=\bigl\{y\in\Gamma:\ \|u^{H}(y)-u^{L}(y)\|_{\infty}\le\epsilon\bigr\},
\qquad
\mu\bigl(\Gamma_{N}(\epsilon)\bigr)=\int_{\Gamma_{N}(\epsilon)}\pi^{d}(z)\,dz .
$$

在 Assumption 1（$\sup_{z\in\Gamma}\|u^H(z)\|=:C_H<\infty$）与独立同分布高斯观测噪声下，存在 $K_1,K_2>0$ 使

$$
D_{\mathrm{KL}}\bigl(\widetilde\pi^{d}_{N}\,\|\,\pi^{d}\bigr)
\le\Bigl(K_1\epsilon+K_2\,\mu\bigl(\Gamma^{\perp}_{N}(\epsilon)\bigr)\Bigr)^{2},
$$

相应的 Hellinger 界为 $K_1\epsilon+K_2\mu(\Gamma^{\perp}_{N}(\epsilon))$。推论给出算法目标：若采样充分使 $\mu(\Gamma^{\perp}_{N}(\epsilon))\le\epsilon$，则 KL 距离完全由 $\epsilon^2$ 刻画。

论文对算法与定理之间的联系是明确表述为**机制**而非定理的：只要候选点落在 $\Gamma^{\perp}_{N}(\epsilon)$ 内，算法就在其邻域细化，因此 $\mu(\Gamma^{\perp}_{N}(\epsilon))$ 随细化渐近衰减。这是一个论证，不是收敛速率。

### 数值证据

两个非线性 PDE 反问题。第一个是二维热源反演，参数维数 $n_z=2$，小到可以负担一个真正精确的高阶先验代理——在这个例子里自适应方法买到的是精度而不是时间。第二个是九参数椭圆方程扩散系数反演，全局精确的先验代理很贵，自适应方法同时改进精度与代价。摘要给出的量级判断是相对纯真模型 MCMC 有若干数量级的效率提升。

## 34：同一个修正装进集合 Kalman 反演

### 换掉采样器

集合 Kalman 反演不需要导数，适合 PDE 约束的参数估计，但它是 Monte Carlo 方法：$N_e$ 个成员、$J$ 轮迭代意味着约 $N_e J$ 次正求解。本文把 37 的多精度修正搬进正则化迭代集合 Kalman 平滑器。更新式为

$$
\theta^{(j)}_{n+1}=\theta^{(j)}_{n}
+C^{\theta\omega}_{n}\bigl(C^{\omega\omega}_{n}+\alpha_{n}\Gamma\bigr)^{-1}\bigl(y^{(j)}-\omega^{(j)}_{n}\bigr),
\qquad j=1,\dots,N_e,
$$

其中 $\omega^{(j)}_n=f(\theta^{(j)}_n)$，$y^{(j)}=y+\xi^{(j)}$ 是扰动数据，经验协方差为

$$
C^{\theta\omega}_{n}=\frac{1}{N_e-1}\sum_{j}\bigl(\theta^{(j)}_{n}-\bar\theta_{n}\bigr)\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)^{\!T},
\qquad
C^{\omega\omega}_{n}=\frac{1}{N_e-1}\sum_{j}\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)^{\!T}.
$$

正则化参数 $\alpha_n$ 取满足

$$
\alpha^{N}_{n}\bigl\|\Gamma^{1/2}\bigl(C^{\omega\omega}_{n}+\alpha^{N}_{n}\Gamma\bigr)^{-1}\bigl(y^{(j)}-\bar\omega_{n}\bigr)\bigr\|
\ \ge\ \rho\,\bigl\|\Gamma^{-1/2}\bigl(y^{(j)}-\bar\omega_{n}\bigr)\bigr\|
$$

的第一个整数 $N$ 对应的值，停机用偏差原理 $\|\Gamma^{-1/2}(y-\bar\omega_{n})\|\le\tau\eta$，其中 $\rho<1$、$\tau\ge1/\rho$。这两条来自 Iglesias 的正则化集合 Kalman 平滑器。

### 指标改成相对形式，检验点改成集合均值

$$
\mathrm{err}=\frac{\bigl\|f(\bar\theta_{n+1})-f_{M}(\bar\theta_{n+1})\bigr\|_{\infty}}{\bigl\|f(\bar\theta_{n+1})\bigr\|_{\infty}},
\qquad
\bar\theta_{n+1}=\frac{1}{N_e}\sum_{j}\theta^{(j)}_{n+1}.
$$

集合均值在这里承担 37 中「被接受点」的角色：它是当前集合认为最可能的位置。样本预算是显式的：初始化用 $Q_1=2\binom{N+d}{d}$ 个先验样本，每次细化用 $Q_2=2\binom{N_C+d}{d}$ 个新点。

### 数值证据中最有信息量的一点

实验刻意让真值落在先验之外：先验取 $\log\theta_i\sim\mathcal N(0,1)$，真值从 $\log\theta_i\sim U(-4,4)$ 抽取。在二维时间分数阶反扩散问题上（Caputo 导数，9 个径向基权重，$N_e=100$，噪声水平 $\sigma=10^{-3}$），固定的先验多项式混沌代理给出相对误差 $0.7921$（阶数 $N=4$）与 $0.2892$（$N=6$），而常规集合 Kalman 反演是 $0.0461$。代价方面，常规方法约 $56.71$ 秒，多项式混沌方法约 $0.82$ 秒。自适应方法在 $N=2$ 时只用 250 次（容差 $10^{-2}$）或 575 次（容差 $10^{-3}$）真模型求解，而常规方法用 2000 次。

第二个例子恢复由 Karhunen–Loève 展开参数化的对数渗透率场，保留 95% 先验能量对应 $d=22$ 个模态，$N_e=300$。固定代理在 $N=2$ 时相对误差 $0.3430$、$N=3$ 时 $0.2146$，而自适应方法即使在 $N=1$ 也达到 $0.0889$。

这组数字说明了一件事：**当真值出先验时，提高全局代理的阶数是低效的修补方式**，而局部修正即使基于极低阶展开也足够。

## 49：把多项式换成复合神经网络

### 为什么换

论文明确指出多项式代理的两个弱点：对正则性有限的映射逼近能力差，基函数个数随参数维数爆炸。神经网络在这两点上更合适。

### 复合多精度网络

关键构造不是「再训练一个网络」，而是把已训练的低精度网络当作**输入变量**送进高精度网络：

$$
f^{H}(z)=\mathcal F\bigl(z,f^{L}(z)\bigr)=\mathcal F\bigl(z,\mathcal{NN}^{L}(z);\theta\bigr),
\qquad
f^{H}(z)\approx\mathcal{NN}^{H}(z;\theta):=\mathcal{NN}\bigl(z,\mathcal{NN}^{L}(z);\theta\bigr).
$$

论文自己给出了与 37 的对比：多项式混沌的修正 $f^{H}\approx f^{L}_{\mathrm{PCE}}+f_{\mathrm{CORR}}$ 是**线性**叠加，而上式学习的是两个精度之间的**非线性相关**。因为两个模型高度相关，$\mathcal{NN}^{H}$ 可以很浅——数值实验用单隐层 50 个神经元。

浅网络在这里不是性能优化而是**防过拟合的约束**：每轮只能负担 $Q$ 次真求解（实验取 $Q=10$），网络容量必须与之匹配。

### 指标

$$
\mathrm{err}(\tilde z)=\frac{\bigl\|f^{H}(\tilde z)-f^{L}(\tilde z)\bigr\|_{\infty}}{\bigl\|f^{H}(\tilde z)\bigr\|_{\infty}} .
$$

相比 37 的绝对误差改为相对误差，检验点 $\tilde z$ 仍由高精度接受概率产生。触发后在 $B(\tilde z,R)$ 内均匀取 $Q$ 点重训 $\mathcal{NN}^{H}$，然后令 $f^{L}\leftarrow\mathcal{NN}^{H}$——注意这一步使得下一轮的「低精度模型」是上一轮的复合网络，复合关系逐轮嵌套。

> [!warning] 原文公式核对：高精度接受概率的分子分母
> 论文式 (10) 把高精度接受概率印成
> $\beta=\min\{1,\ \mathcal L(d,f^{H}(z^{-}))\pi(z^{-})/\mathcal L(d,f^{H}(z^{+}))\pi(z^{+})\}$。
> 在对称提议的 Metropolis-Hastings 中，被提议状态 $z^{+}$ 应在分子上。按上下文这里应为
> $\beta=\min\{1,\ \mathcal L(d,f^{H}(z^{+}))\pi(z^{+})/\mathcal L(d,f^{H}(z^{-}))\pi(z^{-})\}$。
> 本条为阅读者的核对结论，引用时建议回查期刊版。

### 数值证据

基准椭圆方程反问题的两个配置。第一个是九参数渗透率反演，比较三条基线：真模型 MCMC、固定的先验训练网络、以及本文的自适应版本。先验网络用 50 个训练点、4 隐层 × 40 神经元；修正网络单隐层 × 50 神经元；$Q=10$；容差取 $0.1$ 与 $0.05$。第二个例子是 Karhunen–Loève 参数化的高维随机场。结论是固定的先验网络给出可见错误的后验边缘分布，自适应版本恢复真模型 MCMC 的边缘分布，且容差从 $0.1$ 收紧到 $0.05$ 后一致性进一步提高。

本文不含定理，明确说明沿用 37 的分析框架。

## 79：目标导向指标与贪心设计

### 两处实质改动

编号 79 把代理换成算子网络（DeepONet），把采样器换成无迹 Kalman 反演，并且——这是方法论上更重要的一步——把指标从「代理准不准」换成「拟合坏没坏」。

论文先写出诚实的局部模型误差

$$
e_{M}(t):=\mathbb E_{\nu_t}\bigl\|\mathcal G-\widehat{\mathcal G}_t\bigr\|
=\Bigl(\int_{\mathcal M}\bigl|\mathcal G(m)-\widehat{\mathcal G}_t(m)\bigr|^{2}\,\nu_t(dm)\Bigr)^{1/2},
$$

然后明确拒绝它：这需要一个高维积分。替代方案分两步。先在当前代理后验 $\nu_t$ 的 $T$ 个样本 $\mathcal M^{(t)}=\{m^{(t)}_k\}$ 中，用**真模型**挑出数据拟合最好的锚点

$$
r_t=\arg\min_{m\in\mathcal M^{(t)}}\tfrac12\bigl\|y-\mathcal G(m)\bigr\|^{2}_{\Sigma_\eta},
$$

再把指标定义为锚点处的数据残差本身，并用它的**相对变化**触发细化：

$$
e_{D}(t):=\Phi(r_t;y)=\tfrac12\bigl\|y-\mathcal G(r_t)\bigr\|^{2}_{\Sigma_\eta},
\qquad
\frac{\bigl|e_D(t)-e_D(t-1)\bigr|}{e_D(t)}>\epsilon .
$$

数值实验取 $\epsilon=0.01$，重训次数上限 $I_{\max}=10$。

### 贪心设计：多样性与局部性的显式权衡

从 $\nu_t$ 抽一个大候选池 $\Gamma=\{m_1,\dots,m_K\}$，逐点生长子集 $\gamma_Q$：

$$
\hat m_{j+1}=\arg\max_{m\in\Gamma\setminus\gamma_j}
\Bigl\{d\bigl(\widehat{\mathcal G}_t(m),\widehat{\mathcal G}^{j}_t\bigr)-\lambda\|m-r_t\|_{2}\Bigr\},
\qquad
d\bigl(\widehat{\mathcal G}_t(\cdot),\widehat{\mathcal G}^{j}_t\bigr)
=\max_{\hat m\in\gamma_j}\bigl\|\widehat{\mathcal G}_t(\cdot)-\widehat{\mathcal G}_t(\hat m)\bigr\|_{2},
$$

全部实验取 $\lambda=1$。两项刻意对抗：第一项要求在**代理输出空间**中彼此远离，服务于泛化；第二项把选点拉回锚点附近，服务于局部精度。这个打分只用代理在候选池上的预测值，因此几乎不花代价。重训采用迁移学习，从上一轮权重出发而非从头开始。

### 无迹 Kalman 反演

后验近似由随机动力系统给出：

$$
m_{n+1}=r_0+\alpha(m_n-r_0)+\omega_{n+1},\quad \omega_{n+1}\sim\mathcal N(0,\Sigma_\omega),
\qquad
y_{n+1}=\mathcal G(m_{n+1})+\eta_{n+1},
$$

其中 $\alpha\in(0,1]$ 是正则化参数。预测步 $\hat r_{n+1}=\alpha r_n+(1-\alpha)r_0$、$\hat C_{n+1}=\alpha^{2}C_n+\Sigma_\omega$，分析步

$$
r_{n+1}=\hat r_{n+1}+\hat C^{my}_{n+1}\bigl(\hat C^{yy}_{n+1}\bigr)^{-1}\bigl(y_{n+1}-\hat y_{n+1}\bigr),
\qquad
C_{n+1}=\hat C_{n+1}-\hat C^{my}_{n+1}\bigl(\hat C^{yy}_{n+1}\bigr)^{-1}\bigl(\hat C^{my}_{n+1}\bigr)^{T},
$$

$\hat C^{yy}_{n+1}=\mathrm{Cov}[\mathcal G(m_{n+1})\mid Y_n]+\Sigma_\eta$。期望用修正无迹变换在 $2N_m+1$ 个对称 $\sigma$ 点上求值：$m^{0}=r$、$m^{j}=r+c_j[\sqrt C]_j$、$m^{j+N_m}=r-c_j[\sqrt C]_j$。每轮代价是 $2N_m+1$ 次前向求值，论文报告通常在 $O(10)$ 轮收敛。

### Theorem 3.6：这一族中唯一的收敛结果

在三条假设下——(3.2) 对任意 $\epsilon$ 可训练**线性**神经算子 $\widehat G$ 使 $\|\widehat G-G\|_2<\epsilon$；(3.3) $\|G\|_2<H$；(3.4) $G^{T}\Sigma_\eta^{-1}G\succ0$ 且 $\|G^{T}\Sigma_\eta^{-1}G\|_2>C_1$——再加上 $\mathrm{Range}(G^{T})=\mathrm{Range}(\widehat G^{T})=\mathbb R^{N_m}$、$\Sigma_\omega\succ0$、$\Sigma_\eta\succ0$，代理驱动的无迹 Kalman 反演不动点收敛到真模型不动点，且

$$
\bigl\|\widehat C^{-1}_\infty-C^{-1}_\infty\bigr\|_2\le\frac{2\epsilon H H_\eta}{1-\beta},
\qquad
\bigl\|\widehat r_\infty-r_\infty\bigr\|_2
\le\frac{K_1H_\eta H_y}{C_1}\Bigl(1+\frac{2(1+\alpha\beta)K_2H_\eta H^{2}}{(1-\beta)C_2}\Bigr)\epsilon .
$$

两个界都是 $O(\epsilon)$，即关于代理的算子范数误差线性。中间引理 3.5 给出 $\|\widehat G^{T}\Sigma_\eta^{-1}\widehat G\|_2>C_2$，途径是
$\|G^{T}\Sigma_\eta^{-1}G-\widehat G^{T}\Sigma_\eta^{-1}\widehat G\|_2\le 2\epsilon H\|\Sigma_\eta^{-1}\|_2$。

> [!warning] 定理适用范围的限定
> 论文的 Remark 1 说明，为满足 Theorem 3.6 可以让 $\widehat G_\theta$ 线性化，做法是去掉分支网络中的非线性激活、保留主干网络的激活。也就是说，这条定理覆盖的是一个受限的 DeepONet，不是非线性实验中实际使用的网络。引用该定理时应带上这一限定。

### 数值证据

三个基准：Darcy 流、热源反演（Case I 是两参数源位置问题，DeepONet 在 $[0.5,1]^2$ 上用 500 个均匀样本训练，无迹 Kalman 反演从 $[0.6,0.6]$ 出发；Case II 是更高维变体）、以及反应扩散问题。分支与主干网络均为 5 隐层 × 100 神经元、`tanh` 激活；离线训练 $1\times10^{5}$ 步、$N_{\mathrm{prior}}=1000$ 个高斯随机场先验样本；观测数据 $y_{\mathrm{obs}}=y_{\mathrm{ref}}+\max\{|y_{\mathrm{ref}}|\}\delta\xi$。同时测试分布内与**分布外**真值，后者正是固定代理失效的场合。结论是固定代理只给出粗略估计并产生可见错误的反演轨迹，自适应版本随细化单调降低模型误差并接近全阶模型精度；第三个例子在六轮后由停机准则终止。

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

## 覆盖核对

| 内容                                  | 论文   | 覆盖状态                                               |
| ------------------------------------- | ------ | ------------------------------------------------------ |
| 多精度修正的指标结构                  | 37     | 修正项、系数合并、加权最小二乘与权                     |
| $\epsilon$-可行集与 KL / Hellinger 界 | 37, 79 | 假设、结论、算法与定理的连接方式                       |
| 正则化集合 Kalman 反演的完整更新式    | 34     | 更新、协方差、正则化选择、停机准则                     |
| 出先验真值下的定量对比                | 34     | 两个例子的相对误差与真求解次数                         |
| 复合多精度网络                        | 49     | 网络耦合形式、与线性修正的对比、容量约束               |
| 目标导向指标与贪心设计                | 79     | 被拒绝的诚实指标、锚点、相对变化、贪心打分             |
| 无迹 Kalman 反演与线性收敛定理        | 79     | 预测分析步、$\sigma$ 点、两个 $O(\epsilon)$ 界与其限定 |

## 本页原文

- L. Yan and T. Zhou, [_Adaptive multi-fidelity polynomial chaos approach to Bayesian inference in inverse problems_](https://doi.org/10.1016/j.jcp.2018.12.025), J. Comput. Phys. 381 (2019), pp. 110-128（预印本 [arXiv:1807.00618](https://arxiv.org/abs/1807.00618)）。
- L. Yan and T. Zhou, [_An adaptive multifidelity PC-based ensemble Kalman inversion for inverse problems_](https://doi.org/10.1615/Int.J.UncertaintyQuantification.2019029059), Int. J. Uncertain. Quantif. 9(3) (2019), pp. 205-220（预印本 [arXiv:1809.08931](https://arxiv.org/abs/1809.08931)）。
- L. Yan and T. Zhou, [_An adaptive surrogate modeling based on deep neural networks for large-scale Bayesian inverse problems_](https://doi.org/10.4208/cicp.OA-2020-0186), Commun. Comput. Phys. 28 (2020), pp. 2180-2205（预印本 [arXiv:1911.08926](https://arxiv.org/abs/1911.08926)）。
- Z. Gao, L. Yan, and T. Zhou, [_Adaptive operator learning for infinite-dimensional Bayesian inverse problems_](https://doi.org/10.1137/24M1643815), SIAM/ASA J. Uncertain. Quantif. 12(4) (2024), pp. 1389-1423（预印本 [arXiv:2310.17844](https://arxiv.org/abs/2310.17844)）。
