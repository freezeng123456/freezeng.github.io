---
title: 鞅型深度学习
description: 编号 86、93、96、97、100、108：把「残差为零」换成「某个过程是鞅」
lang: zh
translation: en/computational-mathematics/paper-notes/fbsde-and-control/martingale-deep-learning
tags:
  - 论文笔记
  - 随机最优控制
  - 深度学习
---

> [!note] 本页覆盖
> 编号 **86**（_SIAM J. Sci. Comput._ 47(4), 2025）、**93**（_J. Comput. Phys._ 555, 2026）、**96**（投稿 _SIAM Rev._，[arXiv:2408.14395](https://arxiv.org/abs/2408.14395)）、**97**（投稿 SIAM/ASA JUQ，[arXiv:2408.16403](https://arxiv.org/abs/2408.16403)）、**100**（投稿 _Numer. Math._，[arXiv:2510.02635](https://arxiv.org/abs/2510.02635)）、**108**（预印本 arXiv:2604.26782）。
>
> **六篇全部按全文逐式核对**：编号 86 同时核对了 SIAM 出版版与 arXiv 预印本；编号 93、97、100、108 按各自的 arXiv 全文（含附录）；编号 96 覆盖其三个改题版本。因此本页给出完整的推导、定理与**数值实验表**。定理与实验数据的编号沿用各文自身的编号。

![把方程残差改写成鞅性质](assets/diagrams/tao-zhou-papers/zh/martingale-training.svg)

这六篇共享一个转换：高维空间里没有网格，「在每个点检验 PDE 残差」不可行；但「检验一个过程是否是鞅」只需要沿模拟轨道的期望，维数不再直接进入代价。编号 86 首次把这一转换做完整，其余五篇分别沿**去掉导数**（96）、**补上收敛率**（93）、**换掉神经网络**（100）、**换到前向均场方程**（97）、**换到均场博弈**（108）五个方向推进。

## 86：不显式求 $\inf_u H$

### 直觉

高维 Hamilton-Jacobi-Bellman 方程的求解需要在每个时空点求

$$
\inf_{\kappa\in U}H\bigl(t,x,\kappa,\partial_xv,\partial^2_{xx}v\bigr),
\qquad
H(t,x,\kappa,z,p)=\tfrac12\mathrm{Tr}\bigl(p\,\bar\sigma\bar\sigma^\top(t,x,\kappa)\bigr)
+z^\top\bar\mu(t,x,\kappa)+c(t,x,\kappa).
$$

当 $U\subset\mathbb R^m$ 高维、或 $H$ 没有闭式极小点时，这个内层极小化本身受维数灾困扰。而大多数已有深度 PDE 求解器要么需要显式的 $\inf_u H$，要么需要沿时间递推地训练网络。此外，标准 deep BSDE 架构依赖 Euler-Maruyama 的**强**收敛（阶 $1/2$），因为它用到 $X$ 的逐轨道性质。

本文一次去掉三处约束：不需要显式 $\inf_u H$、不需要时间递推训练、只需要**弱**收敛（阶 $1$）。做法是三次替换，每一次都把一个难以直接处理的条件换成一个可以用采样验证的等价条件：逐点极小原理换成积分极小原理，PDE 残差换成鞅性质，条件期望换成对检验函数族的弱条件。

### 问题设定

反馈控制类 $\mathcal U_{\rm ad}:=\{u:[0,T]\times\mathbb R^d\to U\mid u\ \text{Borel 可测}\}$，$U\subset\mathbb R^m$。状态与代价为

$$
X^u_t=x_0+\int_0^t\bar\mu\bigl(s,X^u_s,u(s,X^u_s)\bigr)\mathrm ds
+\int_0^t\bar\sigma\bigl(s,X^u_s,u(s,X^u_s)\bigr)\mathrm dB_s,
$$

$$
J(u):=\mathbb E\Bigl[\int_0^Tc\bigl(s,X^u_s,u(s,X^u_s)\bigr)\mathrm ds+g(X^u_T)\Bigr],
\qquad
\text{求 }u^*\ \text{使}\ J(u^*)=\inf_{u\in\mathcal U_{\rm ad}}J(u).
$$

值函数 $v(t,x)=\inf_uJ(t,x,u)$ 是上述全非线性 HJB 方程的粘性解。本文的实际目标略广一些，是**HJB 型方程** $\partial_tv+\mathcal Lv+\inf_{\kappa\in U}H(t,x,\kappa,\partial_xv,\partial^2_{xx}v)=0$、$v(T,x)=g(x)$，其中 $\mathcal L$ 是一个固定的、**不受控**的二阶算子。第 3 节列出它覆盖的四种情形：控制只出现在漂移里；扩散受控但含不受控部分 $\bar\sigma=\bar\sigma_0+\bar\sigma_1$；已知一个预先控制 $u_0$（取 $\mathcal L=\mathcal L^{u_0}$、$H=(\mathcal L^\kappa-\mathcal L^{u_0})v+c$）；以及 $\bar H(t,x,z,p)=\inf_uH$ 显式可写，此时退化为半线性或全非线性抛物方程。

全节假定 $v\in C^{1,2}$，其依据是非退化条件：$\sigma\sigma^\top(t,x)$ 在 $[0,T]\times\mathbb R^d$ 上一致正定。

### 推导

**第一步：把逐点极小原理换成积分极小原理。** 验证给出的是逐点极小原理 $H(t,X^*_t,u^*(t,X^*_t),\partial_xv,\partial^2_{xx}v)=\inf_{\kappa\in U}H(\cdots)$，而它在 $(t,x)$ 上正是一处维数灾来源。沿一条**无控制**扩散 $X_t=X_0+\int_0^t\mu\,\mathrm ds+\int_0^t\sigma\,\mathrm dB_s$ 定义 Hamilton 过程与代价过程：

$$
H^{u,v}_t:=H\bigl(t,X_t,u(t,X_t),\partial_xv(t,X_t),\partial^2_{xx}v(t,X_t)\bigr),
\qquad
\mathcal M^{u,v}_t:=v(t,X_t)+\int_0^tH^{u,v}_s\,\mathrm ds .
$$

**引理 3.2。** 若 $\int_0^T\mathbb E\bigl[\bigl|\inf_\kappa H(t,X_t,\kappa,\partial_xv,\partial^2_{xx}v)\bigr|\bigr]\mathrm dt<\infty$，且逐点问题在 $\mathcal U_{\rm ad}$ 中有解，则最优控制可以从**积分**条件

$$
\int_0^T\mathbb E\bigl[H^{u,v}_t\bigr]\mathrm dt
=\inf_{\bar u\in\mathcal U_{\rm ad}}\int_0^T\mathbb E\bigl[H^{\bar u,v}_t\bigr]\mathrm dt
$$

中得到。证明只有一行：$\varepsilon_t:=H^{u,v}_t-\inf_\kappa H\ge0$ 逐点成立，而上式迫使 $\int_0^T\mathbb E[\varepsilon_t]\mathrm dt\le0$，故 $\varepsilon_t=0$ 对 $\mathrm dt\times\mathbb P$ 几乎处处成立。

**这一步的收益（注记 3.3）说得很清楚**：极小化对象从「在每个 $(t,X_t)$ 上对 $U$ 求极小」换成「对泛函 $\bar u\mapsto\int_0^T\mathbb E[H^{\bar u,v}_t]\mathrm dt$ 求极小」，后者是一个可用 Monte Carlo 估计的二重积分，且 $t$ 与 $x$ 可独立采样，因此天然可并行。

**第二步：鞅刻画。** **引理 3.4。** 在可积性条件 $\int_0^T\mathbb E[|\partial_xv\,\sigma(t,X_t)|^2]\mathrm dt<\infty$、$\int_0^T\mathbb E[|H^{u,v}_t|^2]\mathrm dt<\infty$、$\mathbb E[|v(T,X_T)|^2]<\infty$ 下，对 $(u,v)\in\mathcal U_{\rm ad}\times C^{1,2}$，

$$
(\partial_t+\mathcal L)v(t,X_t)=-H^{u,v}_t\ \ \text{a.e.}
\qquad\Longleftrightarrow\qquad
\mathcal M^{u,v}_t=\mathbb E\bigl[\mathcal M^{u,v}_T\mid\mathcal F_t\bigr].
$$

正向由 Itô 公式给出：$v(t,X_t)=v(0,X_0)+\int_0^t(\partial_t+\mathcal L)v\,\mathrm ds+\int_0^t\partial_xv\,\sigma\,\mathrm dB_s$，代入方程后 $\mathcal M^{u,v}_t=v(0,X_0)+\int_0^t\partial_xv\,\sigma\,\mathrm dB_s$，是一个鞅。反向由鞅表示定理：$\mathcal M^{u,v}_t=\mathcal M^{u,v}_0+\int_0^tZ_s\mathrm dB_s$，与 Itô 展开比较后

$$
Q_t=\int_0^t\bigl\{(\partial_t+\mathcal L)v+H^{u,v}\bigr\}(s,X_s)\,\mathrm ds
=\int_0^t\bigl\{\partial_xv\,\sigma-Z_s\bigr\}\mathrm dB_s
$$

同时是有限变差过程与连续鞅且 $Q_0=0$，故 $Q\equiv0$。

**这就是本文的中心陈述：PDE 残差为零当且仅当该过程是鞅。**

一条几何假设值得注意：要求无控制扩散的支撑包含最优受控扩散的支撑（$\Gamma_t=\Gamma(X_t)\supset\Gamma(X^*_t)$）。注记 3.1 给出两种安排方式：取 $X_0\sim N(x_0,rI_d)$，$r>0$ 为超参数；或取 $X=X^{u_0}$，$u_0$ 是一个预先的控制。

**注记 3.6 同样重要**：整个表述只用到 $X$ 的**期望与条件期望**，不用逐轨道性质。因此 Euler-Maruyama 只通过它的**弱**阶（$1$）进入，而 deep BSDE 型方法依赖**强**阶（$1/2$）。这是第 4.3 节观察到时间方向一阶收敛的理论原因。

**第三步：从鞅到对抗式 min-max。** 由 $X$ 的 Markov 性，$\mathbb E[\mathcal M^{u,v}_T\mid\mathcal F_t]=\mathbb E[\mathcal M^{u,v}_T\mid X_t]$。为避免直接计算条件期望（原始 DeepMartNet 就是直接算的），论文把它换成对检验函数族 $\mathcal T:=\{\rho:[0,T]\times\mathbb R^d\to\mathbb R\mid\rho\ \text{光滑有界}\}$ 的**弱**条件：

$$
\sup_{\rho\in\mathcal T}\Bigl|\int_0^{T-\Delta t}\mathbb E
\Bigl[\rho(t,X_t)\bigl(\mathcal M^{u,v}_{t+\Delta t}-\mathcal M^{u,v}_t\bigr)\Bigr]\mathrm dt\Bigr|^2=0 .
$$

论证依据是塔性质：$\mathbb E[\rho(t,X_t)(\mathcal M_{t+\Delta t}-\mathcal M_t)]=\mathbb E[\rho(t,X_t)\mathbb E[(\mathcal M_{t+\Delta t}-\mathcal M_t)\mid X_t]]$，因此对**所有** $\rho$ 成立就迫使 $\mathbb E[\mathcal M_{t+\Delta t}-\mathcal M_t\mid X_t]=0$——这正是条件期望的投影性质。这个结构就是对抗式学习：$\rho$ 是判别器。

增广 Lagrange 函数与相应的 min-max 问题为

$$
L(u,v,\rho,\lambda)=\int_0^T\mathbb E\bigl[H^{u,v}_t\bigr]\mathrm dt
+\lambda\Bigl|\int_0^{T-\Delta t}\mathbb E\Bigl[\rho(t,X_t)
\bigl(\mathcal M^{u,v}_{t+\Delta t}-\mathcal M^{u,v}_t\bigr)\Bigr]\mathrm dt\Bigr|^2,
$$

$$
(u,v)=\lim_{\lambda\to+\infty}\
\arg\min_{(\bar u,\bar v)}\Bigl\{\sup_{\rho\in\mathcal T}L(\bar u,\bar v,\rho,\lambda)\Bigr\},
\qquad
\mathcal V=\{v\in C^{1,2}:v(T,x)=g(x)\} .
$$

参数化后由控制网络 $u_\alpha$、值网络 $v_\theta$ 与对抗网络 $\rho_\eta$ 三者训练。

**第四步：两个构造性的架构约束。** 控制网络把取值范围硬约束在 $U=\prod_i[a_i,b_i]$ 内：

$$
u_\alpha(t,x)=a+\frac{b-a}{6}\,\mathrm{ReLU6}\bigl(\psi_\alpha(t,x)\bigr),
\qquad
\mathrm{ReLU6}(y)=\min\{\max\{0,y\},6\},
$$

一般 $U$ 则用一个到 $U$ 的距离罚项 $\bar\lambda\int_0^T\mathbb E[\mathrm{dist}(u_\alpha(t,X_t),U)]\mathrm dt$（注记 3.7）。值网络把终值条件硬编码：$v_\theta(T,x)=g(x)$，$t<T$ 时 $v_\theta(t,x)=\phi_\theta(t,x)$。对抗网络则**刻意做得很浅但取值为向量**：

$$
\rho_\eta(t,x)=\sin\bigl(W_1t+W_2x+b\bigr)\in\mathbb R^r,
\qquad
\eta=(W_1,W_2,b)\in\mathbb R^r\times\mathbb R^{r\times d}\times\mathbb R^r,
$$

论文的说法是「取为向量值以增强对抗训练的稳定性」。

**这两处硬约束与本站[[computational-mathematics/paper-notes/scientific-machine-learning/index|科学机器学习专题]]的取向一致**：能由架构保证的约束就不放进罚项。这里唯一留在罚项里的是鞅条件，而它带一个 $\lambda\to\infty$ 的极限。

**第五步：离散与训练。** 在 $\pi_N=\{0=t_0<\cdots<t_N=T\}$ 上，无控制扩散用 Euler 推进，代价过程增量用**梯形公式**：

$$
\Delta\mathcal M^{\alpha,\theta}_{n+1}
:=v_\theta(t_{n+1},X_{n+1})-v_\theta(t_n,X_n)
-\tfrac12\bigl(H^{\alpha,\theta}_n+H^{\alpha,\theta}_{n+1}\bigr)\Delta t_n .
$$

小批量损失在随机指标集 $A\subset\{0,\dots,N\}\times\{1,\dots,M\}$ 上取

$$
L(\alpha,\theta,\eta,\lambda;A)=\frac1{|A|}\sum_{(n,m)\in A}H^{\alpha,\theta,(m)}_n\Delta t_n
+\lambda\bigl|G(\alpha,\theta,\eta;A)\bigr|^2,
\qquad
G:=\frac1{|A|}\sum_{(n,m)\in A}\rho_\eta\bigl(t_n,X^{(m)}_n\bigr)\Delta\mathcal M^{\alpha,\theta,(m)}_{n+1}\Delta t_n .
$$

算法 3.1 的一轮是：抽 $A_i$；做 $J$ 步下降更新 $\alpha,\theta$；再做 $K$ 步上升更新 $\eta$，并按 $\lambda\leftarrow\min\{\lambda,\bar\lambda+\delta_4|G|^2\}$ 调整罚参数。**注记 3.8 是本文速度的来源**：扩散 $X$ 与控制、值函数都无关，所有 $M$ 条 Euler 轨道**一次性离线生成**；而且损失与训练都不含**时间递推**——这是与 deep BSDE 的结构性差别。

第 3.4 节给出抛物方程的特化：对 $\partial_tv+\mathcal Lv+f(t,x,v,\partial_xv,\partial^2_{xx}v)=0$，令 $\widetilde{\mathcal M}^v_t:=v(t,X_t)+\int_0^tf\,\mathrm ds$，**不需要控制网络**，训练目标退化为 $\theta^*=\arg\min_\theta\{\max_\eta|\widetilde G(\theta,\eta)|^2\}$（算法 3.2）。注意此处 $f$ 可以依赖 $v$ 自身，这一点比 $\bar H$ 更一般。

### 定理

**定理 3.5（鞅型表述）。** 设 $(u,v)\in\mathcal U_{\rm ad}\times C^{1,2}$ 满足上述可积性条件（3.7）与（3.12）。则满足逐点极小原理的最优反馈控制 $u$ 与满足 HJB 方程（对 $t\in[0,T]$、$x\in\Gamma(X_t)$）的值函数 $v$，可以从**两个**条件

$$
\int_0^T\mathbb E\bigl[H^{u,v}_t\bigr]\mathrm dt=\inf_{\bar u\in\mathcal U_{\rm ad}}\int_0^T\mathbb E\bigl[H^{\bar u,v}_t\bigr]\mathrm dt,
\qquad
\mathcal M^{u,v}_t=\mathbb E\bigl[\mathcal M^{u,v}_T\mid\mathcal F_t\bigr],\quad t\in[0,T]
$$

中得到。引理 3.2 与引理 3.4 是它的两半。

> [!warning] 这是一条刻画定理，不是收敛定理
> 本文**没有**给出神经网络逼近的收敛率定理，**没有**泛化界，也**没有**证明 min-max 训练收敛。$v\in C^{1,2}$ 的正则性是全节**假定**的，其依据是 $\sigma\sigma^\top$ 一致正定这一非退化条件。唯一的定量速率主张——$\Delta t$ 方向的一阶弱收敛——是注记 3.6 用 Euler-Maruyama 的弱阶**论证**的，并在第 4.3 节**数值验证**，而非证明。这一缺口后来由编号 93 的定理 4 补上。

### 数值实验

实现：Python 3.12 + PyTorch 2.51，在一个装有 8 块 NVIDIA A100-SXM4-80GB 的节点上用分布式数据并行运行。

**（一）线性抛物（第 4.1 节）。** $(\partial_t+\tfrac12\Delta_x)v-f=0$，精确解 $v(t,x)=1+\tfrac1d\sum_i\sin(t+x_i)$，在线段 $S_1=\{se_1\}$、$S_2=\{s\mathbf 1_d\}$（$s\in[-1,1]$）上求解。维数 $d=100,1000,2000$，5 次独立运行的**运行时间分别为 37、112、363 秒**。

**（二）半线性抛物（第 4.2 节）。** 取自 deep BSDE 文献的算例 $(\partial_t+\Delta_x)v-|\partial_xv|^2=0$，$v(T,x)=1+g(x)$，精确解 $v(t,x)=1-\ln\mathbb E[\exp(-g(x+\sqrt2B_{T-t}))]$，终值**振荡**：$g(x)=\frac1d\sum_i\{\sin(x_i-\frac\pi2)+\sin((\epsilon_0+x_i^2)^{-1})\}$。参考解由 $10^6$ 样本的 Monte Carlo 给出。

**（三）时间收敛率（第 4.3 节）。** 变系数算子 $\mathcal L=\sum_i\sin(2x_i)\partial_{x_i}+\tfrac12\sum_i(1+0.5\sin(5t+x_i))^2\partial^2_{x_i}$，配 **Allen-Cahn 型**源项 $f=v-v^3+\bar f$，$d=100$。**观察到的阶为 $\mathrm{RE}_1=10^{-0.75}N^{-1.01}$，即一阶 $\mathcal O(N^{-1.01})$**，与注记 3.6 的弱阶论证一致。

**（四）不需要显式 $\inf_uH$ 的非退化 HJB（第 4.4 节）。** $(\partial_t+b^\top\partial_x+\epsilon_1\Delta_x)v+\inf_{\kappa\in\mathbb R^d}(2\kappa^\top\partial_xv+c_1|\kappa|^2)=0$，对应的 SOCP 为 $J(u)=1+\mathbb E[\int_0^Tc_1|u_s|^2\mathrm ds+g(X^u_T)]$、$X^u_t=X_0+\int_0^t(b+2u_s)\mathrm ds+\int_0^t\sqrt{2\epsilon_1}\mathrm dB_s$。取 $c_1=\epsilon_1^{-2}$ 时精确解为 $v(t,x)=-\ln\mathbb E[\exp(-g(X^{t,x}_T))]$。三个实例：**HJB-1**（$b=0$，$\epsilon_0=0.1\pi$，$\epsilon_1=1$）、**HJB-2**（$b=\mathbf 1_d$，$\epsilon_0=0.3\pi$，$\epsilon_1=0.2$）、**HJB-3**。维数最高做到 $d=10\,000$，全程不使用显式的 $\inf_\kappa H$。第 4.5 节测试在一个时空区域上的泛化（HJB-2，$d=1000$）；第 4.6 节测试终值代价极小点偏离原点的 SOCP；第 4.7 节则用一个**真正没有闭式极小点**的 Hamilton 量 $\inf_\kappa(2\kappa^\top\partial_xv+c_1|\kappa|^2+\varepsilon\sin(\mathbf 1_d^\top\kappa))$，以 $\varepsilon=0$ 时 HJB-2 的解作基准。

**（五）与 deep BSDE 的效率对比（第 4.8 节，表 1）。** HJB-2，单块 A100；SOC-MartNet 迭代 $I=1000$ 步，deep BSDE 迭代 $I=2000$ 步且**被直接告知显式的 $\inf_uH$**。RE 为相对误差，SD 为标准差，RT 为运行秒数。

| $W_h$  | $d$   | RE (deep BSDE) | RE (SOC-MartNet) | RT (deep BSDE) | RT (SOC-MartNet) |
| ------ | ----- | -------------- | ---------------- | -------------- | ---------------- |
| $256$  | 100   | 3.23E-03       | 1.24E-03         | 73             | **48**           |
| $256$  | 300   | 1.18E-03       | 1.14E-03         | 90             | **53**           |
| $256$  | 500   | 1.05E-03       | 2.89E-03         | 116            | **58**           |
| $256$  | 800   | 2.54E-03       | 4.35E-03         | 145            | **73**           |
| $256$  | 1,000 | 1.86E-03       | 5.83E-03         | 170            | **118**          |
| $d+10$ | 100   | 2.86E-03       | 2.94E-03         | 53             | **20**           |
| $d+10$ | 300   | 3.27E-04       | 8.99E-04         | 103            | **51**           |
| $d+10$ | 500   | 6.41E-04       | 6.92E-04         | 184            | **103**          |
| $d+10$ | 800   | 1.28E-03       | 4.00E-03         | 386            | **255**          |
| $d+10$ | 1,000 | 3.77E-03       | 1.11E-03         | 615            | **360**          |

同表的标准差（SD）为：$W_h=256$ 时 deep BSDE 依次 8.90E-04、8.19E-04、5.71E-04、1.90E-03、2.08E-03，SOC-MartNet 依次 2.43E-03、1.20E-03、1.22E-03、1.41E-03、2.57E-03；$W_h=d+10$ 时 deep BSDE 依次 1.09E-03、1.27E-04、3.58E-04、1.14E-03、5.87E-03，SOC-MartNet 依次 1.67E-03、1.02E-03、4.99E-04、3.78E-03、4.26E-04。

**这张表建立的结论是：精度相当，墙钟时间约为一半，且不需要 Hamilton 量的显式极小点。** 它没有建立的是精度优势——十组里 SOC-MartNet 有六组不如 deep BSDE，$d=1000$、$W_h=256$ 一组差了三倍。

**（六）多 GPU 扩展性（第 4.7 节，表 2）。** HJB-3，$I=6000$，$W_h=d+10$，单位为秒。

| GPU 数        | $d=100$ | $d=500$ | $d=800$ | $d=1{,}000$ | $d=2{,}000$ |
| ------------- | ------- | ------- | ------- | ----------- | ----------- |
| $1\times$A100 | 153     | 775     | 1,350   | 1,909       | 5,032       |
| $2\times$A100 | 151     | 430     | 721     | 1,001       | 2,582       |
| $4\times$A100 | 142     | 233     | 393     | 536         | 1,387       |
| $8\times$A100 | 148     | 153     | 231     | 302         | 773         |

$d=100$ 一行几乎不随 GPU 数变化（153 → 148），说明该规模下瓶颈是通信与启动开销而非计算；$d\ge500$ 起接近线性加速（$d=2000$ 从 5032 秒降到 773 秒，加速比 6.5）。**这条曲线正是「损失不含时间递推」的直接后果。**

> [!warning] 论文自陈的局限（第 5 节）
> 「仍需进一步工作以避免计算 (3.4) 的鞅中值函数 $v(t,x)$ 的那个极大的 $d\times d$ Hessian 矩阵，从而得到一个无导数的方法。」**这正是编号 96 与 100 所填的缺口。**

### 与其他论文的关系

它建立在 **DeepMartNet**（Cai 及其合作者用于边值与特征值问题）之上，加入了对抗学习与控制网络；SOC-MartNet 是其随机控制版本。

与整条经典线索（编号 8、19、23、25、26、35、41、47、63、68）互补：那些方法在低维上以高阶精度求解 FBSDE，条件期望用求积公式算，并有证明了的误差估计；编号 86 把 $d$ 做到 $10^4$，但没有证明的速率。**两者的桥梁是同一个对象**——鞅 $\mathcal M_t=v(t,X_t)+\int_0^tH_s\mathrm ds$ 就是伪装的 BSDE $Y_t=v(t,X_t)$，引理 3.4 的证明就是非线性 Feynman-Kac 论证。

与编号 26、41、50 的路线对照（见[[computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control|二阶 FBSDE 与随机控制一页]]）：那三篇走 Pontryagin 极大值原理与伴随 FBSDE，编号 86 走**动态规划／HJB**，因此从不构造伴随 BSDE，且得到的是一个区域上的**反馈**控制 $u_\alpha(t,x)$ 而非沿单一初值的开环控制。编号 68 是方水鑫（Shuixin Fang）在经典分析方向上的前身工作。注记 3.6 的弱阶对强阶论证，是与 deep BSDE 家族（E-Han-Jentzen；Han-Jentzen-E）以及一切逐轨道概率型格式最锋利的单点概念对照。

代码见 [sx-fang/MartNet](https://github.com/sx-fang/MartNet)。

## 93：把二阶算子换成一阶随机差分

### 直觉

高维抛物方程的深度求解器分成两族，缺陷正好互补。**直接／PINN 型**方法在采样点上极小化强残差 $\|\mathcal Dv-f\|^2$，跨样本点完美并行，但需要用自动微分组装 $d\times d$ 的 Hessian $\partial_{xx}v$，这是内存与时间的瓶颈——作者引用的对比数据是 PINN 在 $d\ge10^4$ 时**内存耗尽**。**SDE／FBSDE 型**方法完全避开 $\mathcal D\hat v$，代价是**丢掉时间方向的并行性**：准线性问题的样本轨道依赖未知的 $Y_t=v(t,X_t)$，训练中必须顺序重新模拟，而倒向格式本身就是顺序的。已有的鞅型方法（编号 86、96、DeepMartNet）同时修好了这两处，但它们经由 Itô 微积分推导——作者认为这对没有随机分析背景的读者是一道门槛——而且**没有证明过收敛率**。

**本文的机制可以用一句话说清：一个二阶微分算子可以由 $v$ 在恰好两点上的一阶差商重现，代价是对一个随机跳跃取期望。** 没有 Hessian，没有梯度，没有自动微分。附带的收获是，这套重推只用 Taylor 展开与初等矩恒等式，并顺带补上了整族方法缺失的误差估计。

### 问题设定

终值型准线性抛物方程

$$
\mathcal Dv(t,x)=f\bigl(t,x,v(t,x)\bigr),\quad(t,x)\in[0,T]\times\mathbb R^d,
\qquad v(T,x)=g(x),
$$

$$
\mathcal D:=\partial_t+\mu^\top\bigl(t,x,v(t,x)\bigr)\partial_x
+\tfrac12\mathrm{Tr}\Bigl[\sigma\sigma^\top\bigl(t,x,v(t,x)\bigr)\partial_{xx}\Bigr],
$$

$\mu\in\mathbb R^d$，$\sigma\in\mathbb R^{d\times q}$。**准线性**指 $\mu,\sigma$ 可以依赖 $v$ 本身。对应的 FBSDE 是一个**生成元不依赖 $Z$** 的系统，非线性 Feynman-Kac 关系为 $Y_t=v(t,X_t)$、$Z_t=\partial_xv(t,X_t)\sigma(t,X_t,Y_t)$。

### 推导

**第一步：基本随机差分展开。** 对四阶导有界的 $F:\mathbb R^q\to\mathbb R$ 与随机向量 $\xi=(\xi_1,\dots,\xi_q)^\top$，在 $z=\sqrt h\,\xi$ 处 Taylor 展开并取期望：

$$
\mathbb E\bigl[F(\sqrt h\,\xi)\bigr]
=F(0)+\sum_{k=1}^{3}\frac{h^{k/2}}{k!}\mathbb E\bigl[(\xi^\top\partial_z)^kF(0)\bigr]
+\frac{h^2}{4!}\mathbb E\bigl[(\xi^\top\partial_z)^4F(c\sqrt h\xi)\bigr].
$$

施加**矩条件**

$$
\mathbb E[\xi_i]=0,\qquad\mathbb E[\xi_i\xi_j]=\delta_{ij},\qquad
\mathbb E[\xi_i\xi_j\xi_k]=0,\qquad\mathbb E\bigl[|\xi_i\xi_j\xi_k\xi_l|\bigr]<\infty,
$$

一阶项与三阶项被杀掉，剩下 $\mathbb E[F(\sqrt h\xi)]=F(0)+\frac h2\sum_i\partial^2_{z_i}F(0)+O(h^2)$。**关键在于左边只含 $F$ 的取值，右边却出现了二阶导。** 可用的 $\xi$（注记 1）：$\xi\sim N(0,I_q)$；或各分量独立、服从三点分布 $\mathbb P(\xi_i=\pm c)=\frac1{2c^2}$、$\mathbb P(\xi_i=0)=1-c^{-2}$（$c\ge1$）。

**第二步：对流扩散算子的随机差分算子（RDO）。** 固定 $(t,x)$，令 $V(s,z):=v(t+s,\,x+\mu s+\sigma z)$（系数冻结在 $(t,x,v(t,x))$），对 $F=V(h,\cdot)$ 用上式，再把 $V(h,0)$ 与 $\sum_i\partial^2_{z_i}V(h,0)$ 在 $h$ 中展开，得

$$
\mathcal D_hv(t,x):=\mathbb E\!\left[\frac{v\bigl(t+h,\;x+\mu h+\sigma\sqrt h\,\xi\bigr)-v(t,x)}{h}\right]
=\mathcal Dv(t,x)+O(h).
$$

注记 3 把它与**随机光滑化 PINN**（RS-PINN）作了对照：后者光滑的是**网络**，$v_\theta(x)=\mathbb E[\phi_\theta(x+\sqrt h\xi)]$，再用 Stein 恒等式 $\partial_xv_\theta=\mathbb E[\frac{\xi}{\sqrt h}\phi_\theta]$、$\partial_{xx}v_\theta=\mathbb E[\frac{\xi\xi^\top-I_d}{h}\phi_\theta]$ 求导。**RS-PINN 的行为像中心差分，而 RDO 像迎风差分**，类似于一个物质导数，因此更适合对流占优的问题。

**第三步：RDM 表述与强形式损失。** 把 $\mathcal D\to\mathcal D_h$ 代入方程，得 $\mathbb E[R(t,x,\xi;v)]=O(h)$，其中

$$
R(t,x,\xi;v):=\frac{v(t+h,\,x+\xi_h)-v(t,x)}{h}-f\bigl(t,x,v(t,x)\bigr),
\qquad
\xi_h:=\mu\bigl(t,x,v(t,x)\bigr)h+\sigma\bigl(t,x,v(t,x)\bigr)\sqrt h\,\xi,
$$

对应的强形式损失是 $\mathcal L_{\rm rdm}(\hat v)=\int_0^{T-h}\int_{\mathbb R^d}|\mathbb E[R(t,x,\xi;\hat v)]|^2p(t,x)\,\mathrm dx\,\mathrm dt$。

**第四步：采样密度就是 Fokker-Planck 方程的解。** 这一步是本文对「该在哪里采样」这一普遍启发式的理论交代。写 $\mathcal D=\partial_t+\mathcal L$，令误差 $\epsilon=\hat v-v$、线性化残差 $\hat r:=(\partial_t+\mathcal L)\hat v-f(t,x,v)$（$\mathcal L,f$ 在**精确** $v$ 处取值），则 $(\partial_t+\mathcal L)\epsilon=\hat r$。与权 $p$ 配对并用 $L^2$ 伴随 $\mathcal L^*$ 得 $\partial_t\int\epsilon p=\int\hat rp+\int\epsilon(\partial_tp-\mathcal L^*p)$。**取 $p$ 解伴随（Fokker-Planck）问题** $(\partial_t-\mathcal L^*)p=0$、$p(0,\cdot)=\delta_{x_0}$，最后一项消失，留下精确的误差表示

$$
\epsilon(0,x_0)=\int_{\mathbb R^d}\epsilon(T,x)p(T,x)\,\mathrm dx
-\int_0^T\!\!\int_{\mathbb R^d}\hat r(s,x)p(s,x)\,\mathrm dx\,\mathrm ds .
$$

于是 $p(s,x)$ **就是** $(0,x_0)$ 处误差对 $(s,x)$ 处残差的敏感度；而 $p(t,\cdot)$ 正是正向 SDE 的 $X_t$ 的密度，所以最优采样就是用弱 Euler-Maruyama 模拟那条 SDE，且用与 RDO 中**同一个** $\xi$ 分布。

**第五步：两条等值关系。** 两处识别都由直接计算给出，RDM 方向**完全不用随机微积分**。

- **与鞅型方法等值。** 令 $M_t:=v(t,X_t)+\int_0^tf(s,X_s,v(s,X_s))\mathrm ds$，鞅条件 $\mathbb E^x_t[M_{t+h}-M_t]=0$ 化为 $\mathbb E^x_t[v(t+h,X_{t+h})]-v(t,x)+\int_t^{t+h}\mathbb E^x_t[f]\mathrm ds=0$。用弱 Euler-Maruyama 逼近条件律、用左矩形公式逼近积分，两处各差 $O(h^2)$，合起来给出 $h\,\mathbb E[R(t,x,\xi;v)]=O(h^2)$。**即 RDM 表述就是编号 86、96 与 DeepMartNet 的离散鞅条件。**
- **与 FBSDE 的隐式 Euler 格式等值。** 在 $[t_n,t_{n+1}]$ 上对倒向方程取 $\mathbb E^x_{t_n}$ 杀掉 Itô 积分，代入 $Y=v(t,X_t)$ 并用同样两处逼近，得 $v(t_n,x)=\mathbb E[v(t_{n+1},x+\xi_h)]-h\,f(t_n,x,v(t_n,x))$，论文称之为赵卫东-陈艳萍-彭实戈 $\theta$ 格式的特例。**注记 5 写出了一般 $Z$ 情形**：$Y$ 部分变成 $v(t_n,x)=\mathbb E[v(t_{n+1},x+\xi_h)]-hf(t_n,x,v,z)$，需要配一条 $Z$ 部分，最简单的一例是
  $$
  z(t_n,x)=\frac1h\mathbb E\Bigl[v(t_{n+1},x+\xi_h)\,\sqrt h\,\xi^\top\Bigr].
  $$
  **这是整份清单里深度学习与经典 FBSDE 格式之间最干净的一座桥**：它就是「乘以 Brown 增量再除以 $h$」这条 $Z$ 的标准表示。

**第六步：为控制方差改用 Galerkin 弱形式。** 强形式损失的期望在平方**里面**，因此用 $M$ 个点、$2K$ 次 $\xi$ 抽样的无偏小批量估计，方差是

$$
\mathrm{Var}\bigl[\hat{\mathcal L}_{\rm rdm}(\hat v)\bigr]=O\!\Bigl(\tfrac1M\bigl(1+\tfrac1{K^2}\bigr)\Bigr),
$$

它**不按标准的 $1/(MK)$ 衰减**。改为对 $\rho\in\mathcal T$ 检验，并把 $p(t,x)$ 吸收进 $X_t$ 的分布，二重积分塌缩成单个期望：

$$
\min_{\hat v\in\mathcal V}\max_{\rho\in\mathcal T}\bigl|\mathcal L(\hat v,\rho)\bigr|^2,
\qquad
\mathcal L(\hat v,\rho):=\int_0^{T-h}\mathbb E\bigl[\rho(t,X_t)R(t,X_t,\xi;\hat v)\bigr]\mathrm dt,
$$

在同样 $2MK$ 次残差求值下，相应估计量达到 $\mathrm{Var}[|\hat{\mathcal L}|^2]=O(1/(MK))$。**无偏小批量乘积**取

$$
|\mathcal L(\hat v,\rho)|^2\approx\mathcal L^\top(\hat v,\rho;A_1)\,\mathcal L(\hat v,\rho;A_2),
\qquad
A_i=N_i\times M_i,\quad M_1\cap M_2=\varnothing,
$$

**使乘积无偏的正是「轨道」指标集的不交性**。

**第七步：网络。** 终值条件硬编码：$t\le t_{N-1}$ 时 $v_\theta=\phi_\theta(t,x)$，$t>t_{N-1}$ 时 $v_\theta=g(x)$。对抗网络是一个浅层**多尺度（MscaleDNN 式）正弦网络**，输出维数很宽：

$$
\rho_\eta(t,x)=\sin\bigl(\Lambda(W_1t+W_2x+b)\bigr)\in\mathbb R^r,
\qquad
\Lambda(y_1,\dots,y_r)=(c_1y_1,\dots,c_ry_r)^\top,\quad c_i=1+(i-1)c .
$$

算法 1 交替 $J$ 步 $\theta$ 下降与 $K$ 步 $\eta$ 上升，并每 $I_0$ 次迭代刷新 $r\%$ 的存储轨道（注记 6：轨道模拟在时间上是顺序的，因此离线做、少刷新）。

**第八步：HJB 扩展用策略改进。** 对 $\inf_{\kappa\in U}\{\mathcal D^\kappa v+c(t,x,\kappa)\}=0$，拆成 $u(t,x)=\arg\min_\kappa\{\mathcal D^\kappa v+c\}$ 与 $\mathcal D^uv+c(t,x,u)=0$ 两步，对 $\mathcal D^u$ 用 RDO，两步都写成弱形式：

$$
\min_{\hat v\in\mathcal V}\ \sup_{\rho\in\mathcal T}\bigl|\mathcal L(\hat u,\hat v,\rho)\bigr|^2,
\qquad
\min_{\hat u\in\mathcal U_{\rm ad}}\ \mathcal L(\hat u,\hat v,\mathbf 1).
$$

**控制那一步取常检验函数 $\mathbf 1$**，这正是把 Hamilton 量的逐点极小化换成对采样 $(t,x)$ 的**平均**极小化——与编号 86、96 是同一个装置。HJB 情形的轨道由**受控** SDE 生成。

### 定理

这是 Cai-Fang-Zhou 鞅型深度学习线索里**第一篇给出收敛率证明**的论文。

- **假设 1。** $v\in C^{2,4}$；$\hat v$ Borel 可测；多项式增长界 $|\mu(t,x,v)|+|\mu(t,x,\hat v)|\le C_g(1+|x|^{p_\mu})$，$\sigma$ 同理带 $p_\sigma$，以及 $\sum_{\alpha\in M_{2,4}}|D^\alpha v(t,x)|\le C_g(1+|x|^{p_v})$；$\xi$ 满足矩条件且 $\mathbb E[|\xi_i|^{\bar m}]<\infty$，$\bar m:=\max\{2p_v,8\}$。
- **定理 1（局部、冻结系数版）。** 对一切 $(t,x)$ 与 $0<h<\min\{1,T-t\}$，
  $$
  \bigl|\hat{\mathcal D}_hv(t,x)-\hat{\mathcal D}v(t,x)\bigr|\le C_{\rm loc}\,h\,\bigl(1+|x|^{\bar p}\bigr),
  \qquad
  \bar p:=p_v+3p_\mu\max\{1,p_v\}+3p_\sigma\max\{2,p_v\},
  $$
  其中 $\hat{\mathcal D},\hat{\mathcal D}_h$ 由 $\hat\mu=\mu(t,x,\hat v(t,x))$、$\hat\sigma=\sigma(t,x,\hat v(t,x))$ 构造，$C_{\rm loc}$ 与 $\hat v,t,h,T,x$ 无关。
- **推论 2（局部截断误差，取 $\hat v=v$）。** $|\mathcal D_hv-\mathcal Dv|\le C_{\rm loc}h(1+|x|^{\bar p})$，即 **RDO 关于 $h$ 是一阶精度**。
- **假设 2（稳定性假设）。** 第三个变元上的 Lipschitz 条件 $|\mu(t,x,y_1)-\mu(t,x,y_2)|\le C_\mu|y_1-y_2|$（$\sigma$ 带 $C_\sigma$，$f$ 带 $C_f$），再加**非线性界** $C_\mu|\partial_xv|+C_\sigma|\sigma||\partial_{xx}v|\le C_{\rm nl}$。$\mathcal D$ 线性时可取 $C_\mu=C_\sigma=0$、$C_{\rm nl}=0$。**假设 3：** Euler-Maruyama 轨道的矩界 $\max_n\mathbb E[|X^m_n|^{2\bar p}]\le C_{\rm EM}$。
- **引理 3（零稳定性）。** 记 $M^p_n[\psi]:=\int|\psi(t_n,x)|^pP_{X_n}(x)\mathrm dx$（对**采样轨道的分布**取的平均 $p$ 次矩），$\Delta v:=v-\hat v$、$\Delta R:=\mathbb E[R(\cdot;v)]-\mathbb E[R(\cdot;\hat v)]$，则当 $0<h\le\min\{1,(24C_{\rm nl}^2+12C_f^2)^{-1}\}$ 时
  $$
  \max_{0\le n\le N-1}M^2_n[\Delta v]\le C_{\rm st1}\exp(C_{\rm st2}T)
  \Bigl\{M^2_N[\Delta v]+\frac TN\sum_{n=0}^{N-1}M^2_n[\Delta R]+Th^2\Bigr\}.
  $$
  证明用的是**倒向离散 Grönwall 不等式**。
- **定理 4（全局误差）。** 在假设 1-3 与同样的步长限制下，
  $$
  \max_{0\le n\le N-1}M^2_n[v-\hat v]\ \le\ C_1\exp(C_2T)
  \Bigl\{M^2_N[v-\hat v]+\mathcal L_{{\rm rdm},\pi}(\hat v)+Th^2\Bigr\},
  $$
  $$
  \mathcal L_{{\rm rdm},\pi}(\hat v):=\frac TN\sum_{n=0}^{N-1}\int_{\mathbb R^d}
  \bigl|\mathbb E[R(t_n,x,\xi;\hat v)]\bigr|^2P_{X_n}(x)\,\mathrm dx,
  $$
  $C_1=C_{\rm st1}\max\{1,2+4C_{\rm loc}^2(1+C_{\rm EM})\}$、$C_2=C_{\rm st2}$，均与 $h,N,\hat v$ 无关。三项各有归宿：第一项对网络恒为零；第二项由对抗训练驱向零；第三项给出**时间一阶精度**，即 $M^2_n[v-\hat v]=O(h^2)$、$L^2$ 范数下 $O(h)$。

> [!warning] 精度只在采样轨道探索到的区域内成立
> 定理 4 的范数 $M^2_n$ 是按 $X^m_n$ 的分布加权的，因此**它只保证在样本轨道所探索的区域内精确**。这是整族鞅型方法（编号 86、93、96、100、108）共有的限制，编号 96 的注记 1 把同一件事表述为「引导过程必须以高概率覆盖关心的区域」，编号 108 的 $d=1$ 可视化里也直接看得到。

- **注记 9（本文最锋利的比较性主张）。** 由于 RDM 表述与离散鞅条件等值，定理 4 **同时证明了编号 86、96 与 DeepMartNet 的时间一阶收敛**。这高于 deep BSDE 型方法在一般系数下典型的 $O(h^{1/2})$。机制是：鞅／RDM 方法只用 Euler-Maruyama 逼近 $X_{t+h}$ 在 $X_t=x$ 下的**条件律**，误差由 EM 的**弱**阶 $O(h)$ 支配；deep BSDE 用 EM 逼近**样本轨道**，因此被**强**阶 $O(h^{1/2})$ 封顶。
- **缺席的部分：** 没有按网络规模或训练过程给出的误差界；没有证明 HJB 的 min-max 问题能控制 $v$ 的误差——**注记 7 明确把这一点列为未决**；没有关于 $d$ 的速率。

### 数值实验

误差报告为在采样点集 $D_n=\{X^m_n\}$ 上的相对 $L^1$ 与 $L^\infty$ 误差 $\mathrm{RE}_1(t_n)$、$\mathrm{RE}_\infty(t_n)$，除非另有说明均在 $t=0$ 处评估；测试曲线为 $S_2=\{s\mathbf 1_d:s\in[-1,1]\}$（直线，长度随 $d$ 增长）与 $S_3=\{l(s)\}$，$l_i(s)=s\,\mathrm{sgn}(\sin i)+\cos(i+\pi s)$（一条在 $\mathbb R^d$ 中盘绕的曲线）。

**（一）带陡梯度的对流扩散，$d=10^3$（第 4.1 节）。** $(\partial_t+\mu^\top\partial_x+\frac{\bar\sigma^2}2\sum_i\partial^2_{x_i})v=0$，$T=2$，$\bar\sigma^2=0.1$，$\mu_i=c\tanh(10x_i)$，终值 $v(T,x)=\frac1d\sum_i\{\tanh(x_i)+\cos(10x_i)\}$（强振荡）；参考解由 $10^6$ 条轨道、EM 步长 $T/100$ 的 Monte Carlo 给出。取 $c=1$ 与 $c=5$，$x=0$ 处的陡梯度都被捕捉到，作者归因于基于 Fokker-Planck 的采样。

**（二）准线性抛物，$d=10^4$ 与 $10^5$（第 4.2 节，表 2）。** 精确解 $v(t,x)=V((t-0.5)\mathbf 1_d+x)$，$V(x)=\sum_{i=1}^{d-1}c_iK(x_i,x_{i+1})+c_dK(x_d,x_1)$，$c_i=(1.5-\cos(i\pi/d))/d$，$K(x_i,x_j)=\sin(x_i+\cos(x_j)+x_j\cos(x_i))$，$T=1$。三个算例：**QLP-1** $\mathcal D=\partial_t+\frac{v^2}2\sum_i\partial^2_{x_i}$、$f=v-v^3+Q$；**QLP-2a** $\mathcal D=\partial_t+(\frac v2-1)\sum_i\partial_{x_i}+\frac{v^2}2\sum_i\partial^2_{x_i}$、$f=v^2+Q$；**QLP-2b** 漂移按 $1/d$ 缩放但扩散是**稠密**的 $\frac1{2d^2}\sum_{i,j,k}\sigma_{ik}\sigma_{jk}\partial_{x_i}\partial_{x_j}$，$\sigma_{ij}=\cos(x_i)+v\sin(x_j)$——作者指出这一例对 PINN 实际上不可行，因为需要完整的 Hessian。

| 方程   | $d$    | 平均 $\mathrm{RE}_1$ | SD      | 平均 $\mathrm{RE}_\infty$ | SD      | 运行时间 (s) |
| ------ | ------ | -------------------- | ------- | ------------------------- | ------- | ------------ |
| QLP-1  | $10^4$ | 1.98E-2              | 7.24E-3 | 3.91E-2                   | 8.75E-3 | 1585         |
| QLP-1  | $10^5$ | 1.29E-2              | 1.40E-3 | 3.20E-2                   | 3.66E-3 | 3761         |
| QLP-2a | $10^4$ | 2.82E-2              | 1.28E-2 | 9.08E-2                   | 1.94E-2 | 1589         |
| QLP-2a | $10^5$ | 4.06E-2              | 1.14E-3 | 1.37E-1                   | 1.33E-2 | 3773         |
| QLP-2b | $10^4$ | 5.77E-2              | 2.00E-3 | 1.12E-1                   | 1.18E-2 | 1602         |
| QLP-2b | $10^5$ | 5.19E-2              | 7.69E-4 | 1.05E-1                   | 1.09E-2 | 3822         |

**值得注意的是维数从 $10^4$ 涨到 $10^5$ 时误差基本不变、时间只涨约 2.4 倍。** 但精度停在 $10^{-2}$ 量级，这与定理 4 的一阶速率一致，也与经典格式六阶精度的量级完全不同。

**（三）HJB，$d=10^4$（第 4.3 节，表 3）。** $\partial_tv+\inf_{\kappa\in\mathbb R^d}\{(b+c\sigma\kappa)^\top\partial_xv+\frac12|\kappa|^2\}+\frac12\mathrm{Tr}[\sigma\sigma^\top\partial_{xx}v]=0$，$c=2$，$T=1$，$v(T,x)=\ln(1+\frac1d\sum_ix_i^2+0.5\sin(10x_i))$。参考解经 **Cole-Hopf 变换** $v(t,x)=-c^{-2}\ln\mathbb E[\exp(-c^2v(T,X^0_T))\mid X^0_t=x]$ 用 $10^6$ 样本、EM 步长 $T/100$ 算出。实例：**HJB-1a** $b_i=\sin(t+i+x_{i+1}-1)$（循环，$x_{d+1}:=x_1$，耦合各分量），$\sigma=0.5I_d$；**HJB-1b** 同上但 $\sigma=0.025I_d$（解更不光滑）；**HJB-2** $b_i=\sin(x_{i+1})$，$\sigma_{ij}=0.5\delta_{ij}\tanh((t-0.5)^2+d^{-1}\sum_kx_k^2)$（变扩散）。网络为全连接，宽度 $W=d+10$。

| 方程   | 平均 $\mathrm{RE}_1$ | SD      | 平均 $\mathrm{RE}_\infty$ | SD      | 运行时间 (s) |
| ------ | -------------------- | ------- | ------------------------- | ------- | ------------ |
| HJB-1a | 8.76E-3              | 8.09E-4 | 3.76E-2                   | 5.74E-3 | 2709         |
| HJB-1b | 1.73E-2              | 8.07E-4 | 4.74E-2                   | 1.94E-3 | 2707         |

**（四）采样策略消融，$d=1$（第 4.4 节）。** 刻意取低维以便可视化，在 $c=1$ 与 $c=5$ 的对流扩散问题上比较「动力学」采样（按 Fokker-Planck 分析所得的 SDE）与「朴素」采样 $X^m_n=X^m_0+B^m_{t_n}$、$X^m_0\sim U[-1,1]$。前者正是第 2.4 节误差传播分析所支持的那一个。

**（五）与 PINN、RS-PINN 的正面对比（第 4.5 节，表 4）。** 在 Hu 等人的椭圆基准上：Allen-Cahn $\Delta v+v-v^3=f$ 与 Sine-Gordon $\Delta v+\sin v=f$，定义域为单位球 $B_d$，$v|_{\partial B_d}=0$，精确解 $v(x)=(1-|x|^2)\sum_{i=1}^{d-1}c_i\sin(x_i+\cos(x_{i+1})+x_{i+1}\cos(x_i))$，$c_i\sim N(0,1)$ 独立同分布。PINN 与 RS-PINN 的数字逐字引自该文献。表中还含**强形式 DRDM**，作者自己指出它与 Xu-Zhang 的 shotgun 方法在离线轨道采样与对偶变量之外基本相同。指标为相对 $L^2$ 误差／运行时间（分钟）／内存（MB）。

| 方法                    | 指标          | $d=10^2$ | $10^3$   | $5\times10^3$ | $10^4$   | $10^5$   |
| ----------------------- | ------------- | -------- | -------- | ------------- | -------- | -------- |
| PINN                    | AC 相对 $L^2$ | 7.187E-3 | 5.617E-4 | 1.773E-3      | 不可行   | 不可行   |
| PINN                    | SG 相对 $L^2$ | 7.192E-3 | 5.642E-4 | 1.782E-3      | 不可行   | 不可行   |
| PINN                    | 时间 (min)    | 3        | 285      | 1832.4        | 不可行   | 不可行   |
| PINN                    | 内存 (MB)     | 1328     | 4425     | 56563         | >81252   | >81252   |
| RS-PINN                 | AC 相对 $L^2$ | 7.923E-3 | 5.504E-4 | 1.802E-3      | 1.860E-3 | 2.192E-3 |
| RS-PINN                 | SG 相对 $L^2$ | 7.835E-3 | 6.744E-4 | 1.795E-3      | 1.854E-3 | 2.176E-3 |
| RS-PINN                 | 时间 (min)    | 1.8      | 7.2      | 31.8          | 66       | 720      |
| RS-PINN                 | 内存 (MB)     | 1413     | 1815     | 3593          | 5789     | 45599    |
| 强形式 DRDM             | AC 相对 $L^2$ | 1.311E-2 | 6.171E-3 | 4.714E-3      | 3.112E-3 | 7.423E-4 |
| 强形式 DRDM             | SG 相对 $L^2$ | 2.009E-2 | 6.191E-3 | 4.715E-3      | 3.112E-3 | 2.232E-2 |
| 强形式 DRDM             | 时间 (min)    | 1.04     | 1.11     | 1.81          | 2.70     | 48.64    |
| 强形式 DRDM             | 内存 (MB)     | 75       | 254      | 1204          | 2389     | 23739    |
| **弱形式 DRDM（本文）** | AC 相对 $L^2$ | 3.963E-2 | 6.211E-3 | 4.699E-3      | 3.118E-3 | 7.381E-4 |
| **弱形式 DRDM**         | SG 相对 $L^2$ | 4.172E-2 | 6.306E-3 | 4.701E-3      | 3.118E-3 | 4.769E-3 |
| **弱形式 DRDM**         | 时间 (min)    | 1.22     | 1.26     | 1.43          | 1.59     | 13.26    |
| **弱形式 DRDM**         | 内存 (MB)     | 57       | 185      | 859           | 1698     | 16829    |

**作者自己的解读异常坦率：DRDM 的精度一般不如 PINN 与 RS-PINN**，他们归因于 RDO 那项 $O(h)$ 截断误差——基于自动微分的方法不产生这一项。收益完全在运行时间与内存上，而且收益随 $d$ 增长：$d=10^5$ 时弱形式 DRDM 用 13.26 分钟、16.8 GB，RS-PINN 用 720 分钟、45.6 GB，PINN 则完全跑不动。弱形式与强形式精度相当，但每步梯度、小批量 $n_b$ 下，强形式要求值 $v_\theta$ 共 $n_b+n_b\ell$ 次（他们的运行中每点抽 $\ell=128$ 个 $\xi$），弱形式只要 $2n_b$ 次，因此 $d$ 越大弱形式越占优。

### 与其他论文的关系

**它是编号 86 与 96 的理论收口。** 那两篇分别引入了鞅型损失与无导数变体，但都没有证明速率；编号 93 指出它们的离散鞅条件与 RDM 表述在代数上是同一个东西，因此定理 4 同时给三者背书了**时间一阶收敛**。作者的代码库把这一识别写得很直白：`DfSocMartNet` 类标注为「DRDM／无导数 MartNet……这两种方法是等价的」。

**它也是深度学习支线与经典 FBSDE 支线的正式重接点。** 第 2.6 节推出隐式 Euler 的 $Y$ 部分，注记 5 给出 $Z$ 部分，并引用 Zhang (2004) 与赵卫东-陈艳萍-彭实戈 (2006)，指向编号 8、47、63 以及赵-李-张／赵-张-鞠的一般理论。换句话说，**RDO 不过是那一族里的一步弱 Euler 成员**，区别只在于条件期望被留作期望，而不是像编号 8、16、19、25 那样用 Gauss-Hermite 求积算出来。相应地，经典诸篇证明高阶（$O(h^k)$，$k$ 最高到 6）但只在 $d\lesssim10$；编号 93 只证一阶，却做到 $d=10^5$。

它是编号 96 的直接后继、编号 100 的同胞：随机差分算子最早出现在编号 96 的 v1 里，编号 93 才把它单独提出、命名、分析并做基准测试；编号 100 是这一无导数三重奏的第三位，走的是局部化而非全局的路线。

## 96：无导数，且时间与空间同时并行

### 直觉

编号 86 仍然需要 $\partial_xv$，更麻烦的是需要 Hamilton 量里那个 $d\times d$ 的 **Hessian $\partial^2_{xx}v$**，每个样本都要靠自动微分算一次。$d\sim10^4$ 时这就是决定性开销，编号 86 自己的结论也把它列为未决问题。第二处瓶颈是：一切基于 SDE 模型的深度方法（deep BSDE、DeepMartNet，以及**受控**情形下的编号 86）都必须模拟依赖未知 $v$ 或 $u$ 的轨道，因此每次网络更新后都要**沿时间顺序重新模拟**。

本文一次去掉两处。**无导数的机制一句话就够：把鞅增量写成只含网络在两个邻近点上的两次取值，中间不出现任何导数。** 时间并行的机制则来自一个结构装置：把鞅性质从**一条**贯穿 $[0,T]$ 的长轨道，换成**一族**短时程过程，每个时间步一个，于是时间步之间解耦。

### 问题设定

准线性抛物方程

$$
(\partial_t+\mathcal L)v(t,x)+f\bigl(t,x,v(t,x)\bigr)=0,\quad(t,x)\in[0,T)\times\mathbb R^d,
\qquad v(T,x)=g(x),
$$

$$
\mathcal L:=\mu^\top\bigl(t,x,v(t,x)\bigr)\partial_x
+\tfrac12\mathrm{Tr}\Bigl\{\sigma\sigma^\top\bigl(t,x,v(t,x)\bigr)\partial^2_{xx}\Bigr\}.
$$

注意这里的**准线性**：$\mu,\sigma$ 依赖 $v$ 本身，比编号 86 的固定 $\mathcal L$ 更一般。全文的常设假设是「$\mu$、$\sigma$ 与 $v$ 足够光滑，以使所涉的截断误差估计成立」。

### 推导

**第一步：引导过程与系统过程。** 这是本文的关键结构装置。

- **引导过程** $\hat X$，只用来探索 $\mathbb R^d$，由一个**初始猜测** $\hat v$ 构造：
  $$
  \hat X_t=\hat X_0+\int_0^t\hat\mu(s,\hat X_s)\mathrm ds+\int_0^t\hat\sigma(s,\hat X_s)\mathrm dB_s,
  \qquad
  \hat\mu(t,x):=\mu(t,x,\hat v(t,x)),\ \ \hat\sigma(t,x):=\sigma(t,x,\hat v(t,x)).
  $$
- **系统过程** $X^s_t$（$0\le s\le t\le T$），**从 $\hat X_s$ 出发**、由真正的 $\mathcal L$ 生成：
  $$
  X^s_t=\hat X_s+\int_s^t\mu(r,X^s_r,v(r,X^s_r))\,\mathrm dr+\int_s^t\sigma(r,X^s_r,v(r,X^s_r))\,\mathrm dB_r .
  $$

**第二步：时间局部化的鞅表述。** 对 $t\mapsto v(t,X^s_t)$ 用 Itô 公式：

$$
\mathcal M^s_t:=v(t,X^s_t)-v(s,X^s_s)+\int_s^tf\bigl(r,X^s_r,v(r,X^s_r)\bigr)\mathrm dr
=\int_s^tR(r,X^s_r;v)\,\mathrm dr+\int_s^t(\partial_xv)^\top\sigma\,\mathrm dB_r,
$$

其中

$$
R(t,x;v):=(\partial_t+\mathcal L)v(t,x)+f\bigl(t,x,v(t,x)\bigr)
$$

**恰是 PDE 残差**。取条件期望消掉 Itô 积分，再取 $t=s+h$，得 $\mathbb E[\mathcal M^s_{s+h}\mid\hat X_s]=h\,R(s,\hat X_s;v)+O(h^2)$，于是鞅型判据为

$$
\mathbb E\bigl[\mathcal M^t_{t+h}\,\big|\,\hat X_t\bigr]=0,
\qquad 0\le t\le T-h .
$$

**这个写法的要点是残差从不被显式计算，却被这个条件所刻画。**

> [!warning] 注记 1：判据只在引导过程探索到的区域内起作用
> 上式只保证残差 $R(t,\hat X_t;v)$ 在**引导过程所探索到的区域内**为零，因此引导过程必须以高概率覆盖关心的区域。**这是这一族方法（编号 86、93、96、100、108）共同的实际风险所在**，在编号 93 中表现为加权范数 $M^2_n[\cdot]$，在编号 108 的可视化里可以直接看到。

**与 DeepMartNet／编号 86 的关键差别**在于：那里鞅性质加在**一条**贯穿 $[0,T]$ 的长轨道 $t\mapsto X^0_t$ 上；这里加在**一族**短时程系统过程 $X^t_{t+h}$ 上，每个时间步一个，从而**把时间步解耦**，于是可以在时间维度上做小批量抽样并完全并行。

**第三步：去掉条件期望——Galerkin 加对抗。** 由塔性质 $\mathbb E[\rho(t,\hat X_t)\mathbb E[\mathcal M^t_{t+h}|\hat X_t]]=\mathbb E[\rho(t,\hat X_t)\mathcal M^t_{t+h}]$，再用一步 Euler 把 $\mathcal M^t_{t+h}$ 换成完全显式、**无导数**的量

$$
\mathcal M(t,x,w;v):=v\bigl(t+h,\ x+\mu(t,x)h+\sigma(t,x)\sqrt h\,w\bigr)-v(t,x)+h\,f\bigl(t,x,v(t,x)\bigr),
\qquad \xi\sim N(0,I_q),
$$

得到 min-max 问题 $\min_{v\in\mathcal V}\sup_{\rho\in\mathcal T}|G(v,\rho)|^2=0$，$G(v,\rho):=\int_0^{T-h}\mathbb E[\rho(t,\hat X_t)\mathcal M(t,\hat X_t,\xi;v)]\mathrm dt$，$\mathcal V=\{v\in C^{1,2}:v(T,x)=g(x)\}$。**注记 2** 给出这一步的逼近误差为 $O(h^2)$：Euler-Maruyama 的弱二阶**局部**截断误差，加上积分的左矩形公式同样局部二阶。

**注意 $\mathcal M$ 里完全没有 $v$ 的导数**，只有网络在两个邻近点上的两次取值。这就是无导数机制的全部。

**第四步：平方期望的无偏小批量估计。** 引导轨道由 Euler 生成，经验损失取 $G(v,\rho;A):=\frac h{|A|}\sum_{(n,m)\in A}\rho(t_n,\hat X^m_n)\mathcal M(t_n,\hat X^m_n,\xi^m_n;v)$，再取

$$
|G(v,\rho)|^2\ \approx\ G(v,\rho;A_1)\,G(v,\rho;A_2),
\qquad A_1\cap A_2=\varnothing .
$$

**用两个不交指标集才使得对平方的估计无偏**；朴素的 $|G(\cdot;A)|^2$ 会被小批量方差系统性地抬高。

**第五步：网络。** 值网络硬编码终值条件；**多尺度对抗网络**比编号 86 的更丰富：

$$
\rho_\eta(t,x)=\sin\bigl(\Lambda(W_1t+W_2x)+b\bigr)\in\mathbb R^r,
\qquad
\Lambda(y_1,\dots,y_r)=(cy_1,2cy_2,\dots,rcy_r)^\top,
$$

尺度层 $\Lambda$ 沿用多尺度神经网络的想法，给判别器一段频率谱。控制网络与编号 86 相同，用 ReLU6 硬约束到 $U=[a,b]$。

**第六步：HJB 扩展用策略改进（PIA）。** 目标为 $\partial_tv+\inf_{\kappa\in U}\{\mathcal L^\kappa v+c(t,x,\kappa)\}=0$，拆成策略改进对：给定 $u$ 时第一式**线性**，用第二至四步处理；控制步则把逐点 $\arg\min$ 换成积分版本，由于 $\mathbb E[\mathcal M(t,\hat X_t,\xi;u,v)]=h\{(\partial_t+\mathcal L^u)v+c\}+O(h^2)$，而 $h^{-1}$ 与 $\partial_tv$ 都不影响极小点，控制步塌缩成

$$
\min_{u\in\mathcal U_{\rm ad}}\ G(u,v,\mathbf 1),
$$

即取常检验函数 $\rho\equiv1$。**这是极小原理的无导数实现**：用「极小化值函数的均值」代替「逐点极小化 Hamilton 量」。

> [!note] 注记 4：一处诚实的保留
> 相对编号 86，损失里确实不再有 $\partial_xv$ 与 $\partial^2_{xx}v$；**但在受控情形下**，随机跳跃 $\xi^{t,x,u}_h:=\mu(t,x,u(t,x))h+\sigma(t,x,u(t,x))\sqrt h\,\xi$ 依赖 $u$，因此**这些跳跃无法在训练前预先算好**，能离线预生成的只有引导轨道 $\hat X$。

第一版（v1）从**随机有限差分算子**出发推出同一表述，并观察到：半线性抛物方程的无导数鞅网络等价于一个以转移密度 $p(t,x)$ 为权的**加权 Galerkin 型弱对抗网络**（WAN），也就是弱形式 PINN。这正是通向编号 93 的概念桥梁。

### 定理

本文是**构造性／算法性**的。它的「结果」是鞅表述、无导数弱形式、无导数极小原理，以及两个算法。

- 唯一的定量估计是**注记 2 的 $O(h^2)$ 局部逼近误差**，来自 Euler-Maruyama 的弱二阶局部截断误差加左矩形求积的局部二阶。在 $O(1/h)$ 步上累积，这就是全局一阶弱精度的常规机制，但**本文没有陈述全局收敛定理**。
- 缺席的部分：没有神经逼近的收敛定理，没有按网络规模的误差界，没有对抗训练收敛的证明，没有关于 $d$ 的速率。**这些缺口中的时间阶部分由编号 93 的定理 4 补上。**

### 数值实验

$u_\alpha$ 与 $v_\theta$ 都是全连接网络，**4 个隐层**、宽度 $W$；引导轨道 $M=10^4$ 条；小批量 $|M_i|$ 在 $d\le1000$、$d=2000$、$d=10^4$ 时分别取 256、128、64。

**（一）Allen-Cahn，$d=100$（第 4.1 节）。** $\partial_tv+\Delta_xv+v-v^3=0$，$v(T,x)=1/(2+0.4|x|^2)$，$T=0.3$，在 $x_0=0$ 处评估。参考值 $v(0,x_0)\approx0.0528$ 来自分支扩散方法。SIAM 排版版本对这一例报告：迭代 500 步时平均相对误差 $3.2\times10^{-3}$、标准差 $2.1\times10^{-3}$，**运行时间不到 6.8 秒**。

**（二）极高维扩散，$d=10^4$（第 4.2 节，表 1）。** 加一个源项 $Q$ 使精确解为 $v(t,x)=V((t-0.5)\mathbf 1_d+x)$，$V$ 同编号 93 第 4.2 节所用的成对耦合函数——刻意做成不可分离、含成对相互作用。迭代 9000 步。

| 设置 | 求解集合 $S$ | 宽度 $W$ | 平均 RE            | RE 的 SD           | 运行时间 (s) |
| ---- | ------------ | -------- | ------------------ | ------------------ | ------------ |
| 1    | $\{-0.5\}$   | 1024     | $5.5\times10^{-3}$ | $2.0\times10^{-3}$ | 295          |
| 2    | $[-1.5,1.5]$ | 1024     | $1.8\times10^{-2}$ | $3.1\times10^{-3}$ | 296          |
| 3    | $[-1.5,1.5]$ | 10240    | $5.4\times10^{-3}$ | $5.3\times10^{-4}$ | 5410         |

**这张表的结论很明确：固定宽度下把求解区域扩大会使精度退化，加宽网络能恢复，代价是约 18 倍的时间。也就是说，瓶颈是网络的表达能力，而不是鞅型机制本身。**

**（三）准线性方程，$d=10^4$（第 4.3 节）。** 三个算例 QLP-1、QLP-2a、QLP-2b 的形式与编号 93 第 4.2 节所列相同（含那个稠密扩散矩阵的算例）。

**（四）HJB 方程（第 4.4 节，表 2）。** $\mathcal L^\kappa=(b+2\kappa)^\top\partial_x+\delta^2\mathrm{Tr}\{\partial^2_{xx}\}$，$c(t,x,\kappa)=\delta^{-2}|\kappa|^2$，$U=\mathbb R^d$，$T=1$；精确解 $v(t,x)=-\ln\mathbb E[\exp(-g(X^{t,x}_T))]$，$X^{t,x}_T=x+(T-t)b+\sqrt2\delta B_{T-t}$，参考解由 $10^6$ 样本的 Monte Carlo 给出。实例：**HJB-1**（$b=0$，$\delta=1$，$g=\ln(0.5(1+|x|^2))$）；**HJB-2a**（$b=\mathbf 1_d$，$\delta=0.1$）；**HJB-2b**（$\delta=0.05$）；**HJB-3a**（$b=\mathbf 1_d$，$\delta=0.2$，$\bar g(x)=\frac1d\sum_i\{\sin(x_i-\frac\pi2)+\sin((0.1\pi+x_i^2)^{-1})\}$，在 0 附近强振荡）；**HJB-3b**（$\delta=0.1$）。$d=10^4$，求 $v(0,s\mathbf 1_d)$（$s\in[-1,1]$），迭代 9000 步。

| 方程   | 隐层数 $H$ | 平均 RE            | RE 的 SD           | 运行时间 (s) |
| ------ | ---------- | ------------------ | ------------------ | ------------ |
| HJB-1  | 4          | $2.2\times10^{-3}$ | $3.4\times10^{-4}$ | 9432         |
| HJB-2a | 4          | $7.5\times10^{-3}$ | $4.5\times10^{-4}$ | 9423         |
| HJB-2b | 4          | $2.1\times10^{-2}$ | $4.3\times10^{-4}$ | 9425         |
| HJB-3a | 4          | $2.4\times10^{-2}$ | $1.7\times10^{-3}$ | 9422         |
| HJB-3b | 6          | $2.3\times10^{-2}$ | $5.4\times10^{-4}$ | 13996        |

**$\delta$ 变小（HJB-2a → 2b）时误差涨了近三倍**，即扩散越弱、解越不光滑，方法越吃力——这与「引导过程要覆盖关心区域」的限制是同一件事的两面。$d=2000$ 的宽度研究（HJB-3b，6000 次迭代）进一步印证表 1 的结论：$W=d+10$ 给出平均 RE $6.9\times10^{-2}$、SD $3.6\times10^{-3}$、540 秒；$W=5d+10$ 给出 $2.0\times10^{-2}$、$1.2\times10^{-3}$、9050 秒。

### 与其他论文的关系

**它是编号 86 的直接后继。** 编号 86 的结论要求一个避开 $d\times d$ Hessian 的「无导数方法」，编号 96 正是交付了这一点，并把编号 86 列为它所改进的对象。两者共享鞅性质的对抗／Galerkin 实施与 ReLU6 控制网络，但编号 96 用短时程系统过程替换了单条长鞅，并从损失中彻底消掉了自动微分。

**它是编号 93 的直接前身：** 编号 96 的 v1 里引入的随机有限差分算子——一个一元、一阶的随机差分却逼近二阶算子——正是编号 93 所发展的对象；v1 关于「结果格式是一个以转移密度为权的弱对抗网络」的观察是二者的概念纽带。编号 100 是这一无导数家族的第三位。

与经典线索（编号 8、19、25、26、41、47、63、68）的对照仍然是那句话：那些方法在 $d\lesssim10$ 上证明阶与稳定性，编号 96 做到 $d=10^4$ 但没有理论。有意思的是**两条线共享同一个结构恒等式**——抛物 PDE 的残差等于过程 $\mathcal M_t=v(t,X_t)+\int f$ 的漂移，这正是 FBSDE 格式所离散的非线性 Feynman-Kac 关系。这里的准线性设定（$\mu,\sigma$ 依赖 $v$）则是编号 8 与 23 所处理的**耦合** FBSDE 在深度学习一侧的对应物。

## 97：DeepSPoC，把粒子的记忆交给网络

### 直觉

非线性 Fokker-Planck 方程（多孔介质、分数阶多孔介质、Keller-Segel、Curie-Weiss）是**均场 SDE** $\mathrm dX_t=b(t,X_{t-},\mu_t)\mathrm dt+\sigma(t,X_{t-},\mu_t)\mathrm dZ_t$（$\mu_t=\mathrm{Law}(X_t)$，$Z$ 为 Lévy 过程）的前向方程。它们无法用朴素 Monte Carlo 模拟，因为均场项破坏了 Markov 性。经典的**混沌传播（PoC）粒子方法**把 $\mu_t$ 换成 $N$ 个相互作用粒子的经验测度，但有两处代价：$N$ 个粒子必须同时推进，很贵；而且重建解需要存储**所有粒子的轨道**，高维下不可承受。

**本文的机制是：让神经网络充当粒子的记忆。** 它建立在 Du-Jiang-Li 的**顺序混沌传播（SPoC）**理论上——每个新粒子只与**在先的**粒子相互作用，于是系统有了迭代形式；再把存储的测度整个换成一个神经网络密度。轨道于是从不被存储。

**还有一处观念上的反常值得先说清楚**：这里的「损失」不是一个被极小化的目标。粒子批次每个 epoch 都在变，损失本身也随之改变；梯度下降只是网络**吸收**新批次的机制，用来模仿加权平均 $\mu^n_t=(1-\alpha_n)\mu^{n-1}_t+\alpha_n\hat\mu^n_t$。这就是为什么每个 epoch 只做**恰好一步**下降，也是为什么**训练中损失值并不下降**（图 2：损失走平，而相对 $L^2$ 误差持续下落）。损失因此不能用作停机准则——这正是第 3.2 节那个后验误差估计的动机。

### 问题设定

经典 PoC 粒子系统与 SPoC 粒子系统分别为

$$
\mathrm dX^{n,N}_t=b(t,X^{n,N}_{t-},\mu^N_t)\mathrm dt+\sigma(t,X^{n,N}_{t-},\mu^N_t)\mathrm dZ^n_t,
\qquad \mu^N_t=\frac1N\sum_{i=1}^N\delta_{X^{i,N}_t},
$$

$$
\mathrm dX^n_t=b(t,X^n_{t-},\mu^{n-1}_t)\mathrm dt+\sigma(t,X^n_{t-},\mu^{n-1}_t)\mathrm dZ^n_t,
\qquad \mu^n_t=\mu^{n-1}_t+\alpha_n\bigl(\delta_{X^n_t}-\mu^{n-1}_t\bigr),
$$

其中 $\alpha_n\downarrow0$ 是更新率。批量版本取 $\hat\mu^n_t=\frac1K\sum_{i=1}^K\delta_{X^{i,n}_t}$ 与 $\mu^n_t=\mu^{n-1}_t+\alpha_n(\hat\mu^n_t-\mu^{n-1}_t)$。

### 推导

**第一步：deepSPoC 系统。** 把存储的测度 $\mu^{n-1}_t$ 换成神经网络密度 $\rho_{{\rm NN},\theta_{n-1}}(t,\cdot)$：

$$
\begin{cases}
\mathrm dX^{i,n}_t=b\bigl(t,X^{i,n}_{t-},\rho_{{\rm NN},\theta_{n-1}}(t,\cdot)\bigr)\mathrm dt
+\sigma\bigl(t,X^{i,n}_{t-},\rho_{{\rm NN},\theta_{n-1}}(t,\cdot)\bigr)\mathrm dZ^{i,n}_t,
& i=1,\dots,K,\\[4pt]
\rho_{{\rm NN},\theta_n}=\mathcal F\bigl(\rho_{{\rm NN},\theta_{n-1}},\hat\mu^n,\alpha_n\bigr).
\end{cases}
$$

算子 $\mathcal F$ 就是：由当前网络与新批次的经验测度构造损失，求 $\nabla_\theta$，走**恰好一步**学习率为 $\alpha_n$ 的梯度下降。学习率扮演 SPoC 更新率的角色，按 $\alpha_n=\alpha_0\gamma^{\lfloor n/\Gamma\rfloor}$（$\gamma<1$）衰减。

**第二步：时间离散与磨光。** Euler-Maruyama 在 $0=t_0<\cdots<t_M=T$、$\Delta t=T/M$ 上推进，$X^i_0\overset{\rm iid}\sim\mu_0$，给出 $\hat\mu_{t_m}=\frac1K\sum_i\delta_{X^i_{t_m}}$。磨光后的经验密度为 $\hat\rho_{t_m}:=\hat\mu_{t_m}*f_\epsilon$，用 Gauss 磨光核 $f_\epsilon(x)=(2\pi\epsilon^2)^{-d/2}e^{-|x|^2/(2\epsilon^2)}$。

**第三步：三种损失。**

1. **$L^2$ 距离损失 $\mathcal L_{\rm sq}$**，配普通全连接网络 $\rho_{FC,\theta}$：
   $$
   \mathcal L_{\rm sq}=\sum_{m=0}^M\bigl\|\rho_{FC,\theta}(t_m,\cdot)-\hat\rho_{t_m}(\cdot)\bigr\|^2_{L^2}
   \approx\frac1N\sum_{m=0}^M\sum_{x\in S}\Bigl|\rho_{FC,\theta}(t_m,x)
   -\frac1K\sum_{i=1}^K\frac{e^{-|x-X^i_{t_m}|^2/(2\epsilon^2)}}{(2\pi\epsilon^2)^{d/2}}\Bigr|^2,
   $$
   $S=\{x_j\}_{j=1}^N$ 从截断盒 $\Omega_0$ 中均匀抽取。**注记 2**：反向传播时 $\hat\rho_{t_m}(x_j)$ 被**分离**（detach），即视作与 $\theta$ 无关的常数，尽管粒子本身是用网络生成的。
2. **KL 损失 $\mathcal L_{\rm kl}$**，配**时间型正规化流**（KRnet）——它自动非负、质量为 1，且可直接采样：
   $$
   \mathcal L_{\rm kl}=-\sum_{m=0}^M\int_{\mathbb R^d}\hat\rho_{t_m}(x)\log\rho_{NF,\theta}(t_m,x)\,\mathrm dx
   \approx-\sum_{m=0}^M\frac1N\sum_{x\in S}\hat\rho_{t_m}(x)\log\rho_{NF,\theta}(t_m,x).
   $$
3. **轨道／负对数似然损失 $\mathcal L_{\rm path}$**，唯一**不需要磨光**的选项，因此也是最高维算例所用的那个：
   $$
   \mathcal L_{\rm path}=-\sum_{m=0}^M\frac1K\sum_{i=1}^K\log\rho_{NF,\theta}\bigl(t_m,X^i_{t_m}\bigr).
   $$

**第四步：自适应空间采样。** 高维下均匀采样 $\Omega_0$ 会失效，因为密度的支撑只占盒子的极小一部分。修法是复用算法自己的粒子：$S_m=S_{\rm uniform}\cup S^m_{\rm adaptive}$，其中 $S^m_{\rm adaptive}=\{x_j=\hat X^j_{t_m}+\sigma\delta_j\}$，$\delta_j\overset{\rm iid}\sim N(0,1)$，$\hat X^j_{t_m}$ 取自上一个 epoch 的粒子，Gauss 抖动 $\sigma$ 保留探索性。KL 损失的自适应策略不同（对积分做重要性采样）；$\mathcal L_{\rm path}$ 根本不需要训练点集，因此不设自适应。三个算法分别是：全连接网络加 $\mathcal L_{\rm sq}$（算法 1）、正规化流加 $\mathcal L_{\rm kl}$（算法 2）、正规化流加 $\mathcal L_{\rm path}$（算法 3）。

### 定理

分析刻意放在一个**简化的替身设定**里：密度用截断 Fourier 基表示，配 $L^2$ 投影 $P_N$，于是对 $\mathcal L_{\rm sq}$ 走一步梯度**恰好**是一个仿射更新。

- **精确更新恒等式。** 在 Fourier 替身里 $\mathcal F$ 的作用是 $\rho^{\theta_t}_n=(1-2\alpha_n)\rho^{\theta_t}_{n-1}+2\alpha_nP_N(\hat\rho_t)$，展开得 $\rho^{\theta_t}_{n-1}=P_N(\sum_{l=0}^{n-1}\frac{\beta_l}K\sum_i\delta^\epsilon_{\tilde X^{i,l}_t})$，$\beta_l=(1-2\alpha_{n-1})\cdots(1-2\alpha_{l+1})2\alpha_l$，且 $\sum_{l}\beta_l=1$。**这就是「deepSPoC 就是 SPoC」的精确含义：网络更新是所有历史批次的加权平均。**
- **基规模条件。** $\hat\rho_t$ 的 Lipschitz 常数为 $\tilde C/(K\epsilon^{d+1})$，故 $\|\hat\rho_t-P_N\hat\rho_t\|_\infty\le C(\epsilon,K,L,d)(\log N)^d/N$；取 $N$ 使 $\|\rho-P_N\rho\|_\infty\le\epsilon/(2L)^{d+1}$，即得 $W_1(\rho,P_N\rho)\le\epsilon$ 与 $|\|\rho\|_1-\|P_N\rho\|_1|\le\epsilon$。**选 $W_1$ 而非一般 $W_r$ 是有意的**：由 $\|\cdot\|_\infty$ 控制的齐次性只在 $r=1$ 时成立。
- **假设 3.1。** $|b(t,x,\mu)-b(t,y,\nu)|+\|\sigma(t,x,\mu)-\sigma(t,x,\nu)\|\le C(|x-y|+W_1(\mu,\nu))$，且 $|b(t,0,\mu)|+\|\sigma(t,0,\mu)\|\le C(1+\|\mu\|_1)$。
- **命题 3.2（粒子矩的一致界，为截断辩护）。** 设 $\mu_0\in\mathcal P_3$ 且 $N$ 满足基规模条件，则存在与 $L$ 无关的 $C_0$ 使 $\mathbb E[\sup_{0\le s\le T}|X^{i,n}_s|^3]\le C_0$ 且 $\mathbb E[\mathbf 1_{\{\sup_s|X^{i,n}_s|^2\ge L_0\}}\sup_s|X^{i,n}_s|^2]\le C_0/L_0$。
- **定理 3.3（收敛性）。** 在假设 3.1、$\mu_0\in\mathcal P_3$、$N$ 满足基规模条件、**批权相等** $\beta_l=1/n$，以及截断半径 $L$ 取得使 $\mathbb E[\sup_s|X^{i,n}_s-\tilde X^{i,n}_s|^2]\le\epsilon^2$ 之下，存在与 $n,K,L,\epsilon$ 无关的 $C$ 使
  $$
  \mathbb E\,W_1^2\bigl(\rho^{\theta_t}_{n-1},\mu_t\bigr)\ \le\ C\Bigl(\epsilon^2+(Kn)^{-\frac{1}{1+d/2}}\Bigr).
  $$
  第二项是经验测度／PoC 的常规速率，$Kn$ 就是模拟过的粒子总数；第一项是磨光与截断造成的下限。**作者自己写明的保留：这条结论成立于 Fourier 投影替身，而非真正的神经网络。**
- **后验误差估计（第 3.2 节）。** 限制到 Brown 驱动的情形，在**假设 3.2**（同样的 Lipschitz 条件但用 $W_2$ 代替 $W_1$；论文指出它比 3.1 **更弱**，因为 $W_1\le W_2$）下，在 $\mathcal P_{2,\infty}([0,T])$ 上定义度量 $H_\alpha(\mu_\cdot,\nu_\cdot):=(\int_0^Te^{-\alpha t}W_2^2(\mu_t,\nu_t)\mathrm dt)^{1/2}$ 与解映射 $\Phi(\mu_\cdot):=t\mapsto\mathrm{Law}(X_t)$，其中 $X$ 服从**冻结了均场项**、因而是 Markov 且可直接模拟的 SDE。
  - **命题 3.5：** 取 $C_0=2(T+1)C^2$ 与 $\alpha>C_0$，则 $H_\alpha(\Phi(\mu_\cdot),\Phi(\nu_\cdot))\le\sqrt{C_0/(\alpha-C_0)}\,H_\alpha(\mu_\cdot,\nu_\cdot)$，即 $\Phi$ 是压缩映射。
  - **定理 3.6：** 设 $\mu^*_\cdot$ 是不动点（真解），则对 $\alpha>2C_0=4(T+1)C^2$ 与**任意** $\mu_\cdot\in\mathcal P_{2,\infty}([0,T])$，
    $$
    H_\alpha(\mu_\cdot,\mu^*_\cdot)\ \le\ \Bigl(1-\sqrt{\tfrac{C_0}{\alpha-C_0}}\Bigr)^{-1}H_\alpha\bigl(\mu_\cdot,\Phi(\mu_\cdot)\bigr).
    $$
  计算上的读法是：把学到的密度喂进那条已成为 Markov 的 SDE，模拟，再把所得的律与喂进去的比较——这个差就是真误差的一个**可计算上界**。**这是本页六篇里唯一的可计算误差界**（编号 93 给的是先验速率）。
- 缺席的部分：对真实神经网络参数化、对 Lévy 驱动或奇异核情形、对 KL 与轨道损失都没有分析；作者也明确指出若干测试方程「超出了 SPoC 理论甚至 PoC 理论的范围」。

### 数值实验

网络：全连接为 6 个隐层 $\times$ 512 神经元、ReLU；时间型 KRnet 为 10 个仿射耦合块，每块两个 512 神经元隐层；优化器 Adam；PyTorch。**测试集刻意选的是退化性与非局部性，而不是阶的验证——全文没有报告任何 $\Delta t$ 方向的收敛阶。**

**（一）多孔介质方程（第 4.1 节）。** $\partial_t\rho=\Delta\rho^m$，准线性且在 $\rho=0$ 处退化，基准为 **Barenblatt 解** $U_{m,C}(t,x)=t^{-\alpha}\{C-\frac{m-1}{2m}\beta\frac{|x|^2}{t^{2\beta}}\}_+^{1/(m-1)}$，$\alpha=\frac d{d(m-1)+2}$、$\beta=\alpha/d$——一个在自由边界处失去经典可微性的弱解。各维的配置为：

| 维数   | 算法                     | 参数                                                                                                                         |
| ------ | ------------------------ | ---------------------------------------------------------------------------------------------------------------------------- |
| 1D     | 算法 1（均匀）           | $t_0=1$，$T=1$，$\Delta t=0.01$，$N=1000$                                                                                    |
| 3D     | 算法 1（自适应）         | $t_0=0.1$，$T=0.2$，$\Delta t=0.005$，噪声强度 $0.2$，$\Omega_0=[-2,2]^3$                                                    |
| 5D     | 算法 1（自适应）与算法 3 | 2000 均匀点 + 4000 自适应点，$\sigma=0.3$，$\epsilon=0.05$，$K=8000$，$\Omega_0=[-3,3]^5$                                    |
| 6D、8D | 仅算法 3                 | $t_0=1$，$T=1.5$，$\delta t=0.025$，$10^4$ 样本，$\Omega_0=[-3,3]^6$／$[-3,3]^8$，每 $\Gamma=500$ epoch 按 $\gamma=0.5$ 衰减 |

**一条明确的发现：算法 1 在 6D／8D 上很难用，因为高维磨光不准；正是不需要磨光的 $\mathcal L_{\rm path}$（算法 3）使这些维数变得可达。**

**（二）确定性粒子变体（第 4.1.5 节）。** 同一方程改从 ODE $\mathrm dX_t=\nu m\rho^{m-2}(X_t)\nabla\rho(X_t)\mathrm dt$（Carrillo-Craig-Patacchini 的 blob 方法／Wasserstein 梯度流形式）入手。因为需要 $\nabla\rho$，激活函数从 ReLU 换成 **SoftPlus** $\frac1\beta\log(1+e^{\beta x})$（$\beta=20$），导数由 PyTorch 自动微分提供。**这说明 deepSPoC 对底层粒子方法是随机还是确定性并不敏感。**

**（三）Keller-Segel（第 4.2 节）。** $\partial_t\mu=\nabla\cdot((\nabla W*\mu)\mu)+\Delta\mu$，$d=2$，核 $W(x)=\frac1{2\pi}\ln|x|$ **奇异**——一个依赖分布且相互作用奇异的情形，落在 PoC 理论之外。算法 1 与算法 3 都做了。

**（四）Curie-Weiss 均场方程（第 4.3 节）。** $\mathrm dX_t=\{-\beta(X_t^3-X_t)+\beta K\mathbb EX_t\}\mathrm dt+\mathrm dB_t$，$\beta=1$，$K=-0.1$——依赖**期望**的均场 SDE。长时行为与显式不变密度 $p^*(x)=C^{-1}\exp\{-2\beta(\frac{x^4}4-\frac{x^2}2)\}$ 比较，也与一个用 **500 万粒子**算出的 PoC 参考解比较；deepSPoC 用 5000 个 epoch $\times$ 每 epoch 1000 个粒子，**粒子总数相同，但从不同时驻留内存**——这就是内存论证的具体化。

**（五）分数阶多孔介质方程（第 4.4 节）。** $\partial_t\rho=-(-\Delta)^{\alpha/2}(|\rho|^{m-1}\rho)$，$m>1$，$\alpha\in(0,2)$，分数阶 Laplace 算子写成主值积分——一个**非局部**算子，驱动过程是 $\alpha$ 稳定 Lévy 过程而非 Brown 运动。

### 与其他论文的关系

**它在方法论上是本页的异类**，也是全清单里唯一处理**前向**（Fokker-Planck／McKean-Vlasov）问题而非后向（FBSDE／HJB）问题的一篇。编号 8、47、63、86、93、96 都用条件期望把信息从终值条件**向后**传播，而 deepSPoC 把一个**律**向前传播并用网络拟合它。

工程上的血缘却是共享的：网络取代存储的轨道；训练由新鲜模拟的 Euler-Maruyama 轨道驱动；高维下自适应／基于动力学的空间采样不可或缺。尤其是编号 93 从 Fokker-Planck 伴随方程推出的「在过程真正去的地方采样」原则，在这里字面上就是被学习的那个对象。

**与编号 33、61 的均场联系**（见[[computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control|二阶 FBSDE 一页]]与[[computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|多步格式一页]]）：那两篇在低维上用高阶 $\theta$ 格式与多步格式离散**均场 BSDE／FBSDE**，即生成元依赖解的律的后向方程；编号 97 处理同一圈问题的前向 McKean-Vlasov 一侧，没有阶的保证但把 $d$ 做到 8。两者是「分布依赖 SDE 的数值」这一课题互补的两半。与编号 108 的联系是：这里的 McKean-Vlasov 前向动力学恰是均场博弈中人口演化的那一半。

**后验估计是它最独特的地方。** 定理 3.6 是一条真正**可计算**的误差界，本页其余各篇（编号 86、93、96、100、108）都没有；编号 93 给的是先验速率。其机制——冻结非线性、解所得的 Markov 问题、量出失配——是不动点／压缩论证，而不是截断误差论证。

## 100：把网络换成局部线性回归

### 直觉

在 $d\gg1$ 的半线性抛物方程上，两大族方法各有缺口。**深度学习求解器**（deep BSDE、SOC-MartNet、DRDM）能上极高维，但本文的诊断很直接：「优化过程的稳定性有限，对超参数高度敏感，且缺乏严格的先验误差估计。」**经典概率型／回归 Monte Carlo BSDE 求解器**（Gobet-Lemor-Warin 一系）可以分析，但它们用**全局**基函数回归条件期望，基的个数随 $d$ 组合爆炸，$d\approx200$ 以上就病态；而且它们用 Picard 迭代解耦合的 $(Y,Z)$ 系统，代价在高维下失控。稀疏网格停在 $d\approx10$ 左右；分支扩散表示在长时程上方差爆炸且代价 $O(d^2)$。

**本文的机制是：保留鞅型时间离散与无导数精神，但把神经网络整个换成一小群粒子上的局部线性回归。** 关键的两处判断是：条件期望用**全体粒子的系综平均**近似（而粒子数 $M$ 只通过 $e^{-cM}$ 进入误差，因此 $M\le100$ 就够）；$Z$ 由**加权最小二乘拟合 $\nabla u$** 得到，从而把顺序改成 $X\to Z\to Y$，不再需要 Picard 迭代。换来的是可解释性、严格的先验界，以及笔记本量级的成本。

### 问题设定

终值型半线性抛物方程

$$
(\partial_t+\mathcal L)u(t,x)+f\bigl(t,x,u(t,x),\sigma^\top\nabla u(t,x)\bigr)=0,\qquad u(T,x)=g(x),
\qquad
\mathcal Lu=\tfrac12\mathrm{Tr}\bigl[\sigma\sigma^\top\mathrm{Hess}_xu\bigr]+\langle\mu,\nabla u\rangle,
$$

配**解耦** FBSDE 表示 $Y_t=u(t,X_t)$、$Z_t=\sigma^\top(t,X_t)\nabla u(t,X_t)$。**假设 2.1：** 各变元上全局 Lipschitz（常数 $L$），并有线性增长 $\|\mu\|+\|\sigma\|+|f|+|g|\le C(1+\|x\|+|y|+\|z\|)$。**引理 2.1：** 再加一致非退化 $\xi^\top\sigma\sigma^\top\xi\ge\lambda\|\xi\|^2$，该 FBSDE 在 $S^2(\mathbb R^d)\times S^2(\mathbb R)\times H^2(\mathbb R^d)$ 中有唯一适应解。

### 推导

**第一步：鞅／条件期望型时间离散。** 在 $[t_k,t_{k+1}]$ 上积分倒向方程并取 $\mathbb E_k[\cdot]$ 杀掉 Itô 积分，得 $Y_k=\mathbb E_k[Y_{k+1}+\int_{t_k}^{t_{k+1}}f\,\mathrm ds]$；把系数冻结在 $t_k$（左矩形，即 Euler）给出**对 $Y$ 隐式**的半离散格式

$$
\tilde Y_k=\mathbb E_k\bigl[\tilde Y_{k+1}+f(t_k,X_k,\tilde Y_k,\tilde Z_k)\Delta t\bigr],\qquad 0\le k\le N .
$$

**这就是 $\theta=1$ 的那个成员**，与编号 93 注记 5 写下的 $Y$ 部分相同，也是编号 8、47、63 各格式的祖先。

**第二步：小规模随机粒子方法。** 模拟 $M$ 条**独立**的 Euler-Maruyama 轨道，条件期望用**全体 $M$ 个粒子的系综平均**近似：

$$
\tilde Y^m_k\approx\frac1M\sum_{j=1}^M\tilde Y^j_{k+1}+f\bigl(t_k,\tilde X^m_k,\tilde Y^m_k,\tilde Z^m_k\bigr)\Delta t,
\qquad 1\le m\le M,
\qquad \tilde Y_0=\frac1M\sum_m\tilde Y^m_0 .
$$

**关键在于 $M$ 很小**——一般 $M\le100$——因为按定理 3.1，$M$ 只通过一个被指数压制的坏事件概率 $e^{-cM}$ 进入误差，而不是通过 $1/\sqrt M$ 的 Monte Carlo 项。

**第三步：$Z$ 由加权局部线性回归得到（本文最有特色的一步）。** 不去 Picard 迭代那个耦合的 $(d+1)$ 维系统，而是**先解耦、先算 $Z$**。在锚点 $\tilde X^m_k$ 附近局部展开 $u(t_{k+1},\cdot)$：

$$
u(t_{k+1},\cdot)\approx
\underbrace{u(t_k,\tilde X^m_k)+\partial_tu(t_k,\tilde X^m_k)\Delta t}_{=:\alpha}
+\underbrace{\nabla u(t_k,\tilde X^m_k)}_{=:\alpha_x}{}^\top(\cdot-\tilde X^m_k),
$$

再把 $(\alpha,\alpha_x)\in\mathbb R^{d+1}$ 用加权最小二乘拟合到**已经算出的**值 $\tilde Y^j_{k+1}$ 对**当前**位置 $\tilde X^j_k$ 的回归上：

$$
J(\boldsymbol\alpha)=\sum_{j=1}^Mw_j\bigl(\tilde Y^j_{k+1}-\alpha-\alpha_x^\top D_j\bigr)^2,
\qquad D_j:=\tilde X^j_k-\tilde X^m_k,
\qquad w_j=\frac{K(\|D_j\|/\varepsilon_k)}{\sum_{i=1}^MK(\|D_i\|/\varepsilon_k)},
$$

$K$ 为核函数（实现中用 Gauss 核）。于是 $\tilde Z^m_k=\sigma^\top(t_k,\tilde X^m_k)\alpha_x$，**截距 $\alpha$ 被丢弃**。三处实现细节值得记下：

- **免矩阵求解。** 用 LSQR／PCG，只需要两个乘积 $\beta_j=w_j(\alpha+D_j^\top\alpha_x)$ 与 $(\mathbf D^\top\beta)_0=\sum_j\beta_j$、$(\mathbf D^\top\beta)_{1:d}=\sum_j\beta_jD_j$。存储 $O(d)$，每时间步代价 $O(Md)$，$\mathbf D^\top W\mathbf D$ 从不显式形成。
- **岭正则化（注记 2.2）。** $M\ll d$ 是常规工况，法方程欠定，因此实现极小化 $J_\lambda(\boldsymbol\alpha)=\sum_jw_j(\cdots)^2+\lambda\|\boldsymbol\alpha\|^2$。
- **为什么用核权而不用 $k$ 近邻（注记 2.3）。** 高维下距离集中会使 $k$ 近邻的半径被迫膨胀；但对 Gauss 核，权的**比值** $w_j/w_i=\exp\{-(\|D_j\|^2-\|D_i\|^2)/\varepsilon_k^2\}$ 仍然按相对差 $|\|D_j\|-\|D_i\||$ 区分，即使绝对距离已经集中。用**全部**粒子还顺带去掉了常规局部线性回归里那个半径超参数。

**第四步：$Y$ 由标量 Newton 迭代得到。** 固定 $\tilde X^m_k$ 与 $\tilde Z^m_k$ 后，上式是一个**一维**求根问题 $F(\tilde Y^m_k)=0$，用 Newton 法解，报告每步只需 2 至 3 次迭代。

**算法 2.1** 的整体流水线是：并行地置终值 $Y^j_N=g(X^j_N)$；做**一次**前向扫描生成全部 $\tilde X^j_k$（轨道被存储而非重算）；再对 $k=N-1,\dots,0$ 倒推，每个粒子上依次做免矩阵局部线性回归求 $\alpha_x$、置 $\tilde Z_k\leftarrow\sigma^\top\nabla u$、Newton 解 $\tilde Y^m_k$，**这一层对 $m$ 完全并行**；最后取 $\tilde Y_0=\frac1M\sum_m\tilde Y^m_0$。**整条流水线的顺序是 $X\to Z\to Y$，与回归 Monte Carlo 耦合求解 $(Y,Z)$ 正相反。**

### 定理

- **引理 3.1。** Euler-Maruyama 的标准强阶 $1/2$：$\max_{0\le t\le T}\mathbb E\|X_t-\tilde X_t\|^2\le C\Delta t$，且 $\mathbb E[\sup_t\|X_t-\tilde X_t\|^2]\le C\Delta t$。
- **引理 3.3（倒向格式的局部截断误差）。** 若 $f\in C^{1,2}$ 且满足 Lipschitz 条件，则
  $$
  |\mathcal E_k|:=\Bigl|\mathbb E_k\Bigl[\int_{t_k}^{t_{k+1}}f(s,X_s,Y_s,Z_s)\mathrm ds\Bigr]-f(t_k,X_k,Y_k,Z_k)\Delta t\Bigr|\le C(\Delta t)^2 .
  $$
- **引理 3.8（技术核心：局部线性回归梯度估计的误差）。** 设 $\alpha_x$ 是带扰动数据 $\tilde Y^j_{k+1}=Y^j_{k+1}+\delta^j_{k+1}$ 的有限样本加权最小二乘极小点，误差 $\{\delta^j_{k+1}\}$ 在 $\mathcal F_{t_k}$ 下独立同分布，则
  $$
  \mathbb E_k\bigl[\|\alpha_x-\nabla u(t_k,x)\|^2\bigr]\ \le\ C\varepsilon_k^2
  \;+\;C\varepsilon_k^{-2}\,\mathbb E_k\bigl[|\delta_{k+1}|^2\bigr]
  \;+\;C\,e^{-C_{A_k}M\varepsilon_k^d}.
  $$
  三项分别是**偏差**（一阶 Taylor 截断，$O(\varepsilon_k^2)$）、**方差放大**（噪声除以带宽的平方，局部线性回归的经典偏差-方差权衡），以及「$\varepsilon_k$ 球内采样不足」这一坏事件的概率，它按 $e^{-CM\varepsilon_k^d}$ 衰减。
- **定理 3.1（全局误差）。** 记 $\delta_k:=\tilde Y_k-Y_k$。在假设 2.1 与引理 3.3 的前提下，对一切充分大的 $M$，
  $$
  \mathbb E\bigl[|\delta_0|^2\bigr]\ \le\ C\,\Delta t\;+\;C\,\Delta t\,e^{-c_1M},
  $$
  其中 $C>0$ 只依赖 $T$ 与 Lipschitz 常数 $L$、与 $\Delta t$ 无关，$c_1$ 依赖 $\varepsilon_k$。证明结构是标准的：精确解减离散解，拆成 Monte Carlo 偏差 $\frac1M\sum_j\tilde Y^j_{k+1}-\mathbb E_k[\tilde Y_{k+1}]$、传播误差 $\mathbb E_k[\delta_{k+1}]$、Lipschitz 项 $|f(t_k,X_k,\tilde Y_k,\tilde Z_k)-f(t_k,X_k,Y_k,Z_k)|\le L(|\delta_k|+\|\sigma\|\|\alpha_{x,k}-\nabla u_k\|)$、截断误差 $\mathcal E_k$ 四部分，再用条件 Jensen、带参数 $\eta$ 的 Young 不等式与离散 Grönwall。
- **这条界说的是什么。** 估计「对粒子数 $M$ 与时间步长 $\Delta t$ 都是完全显式的」；因为 $M$ 只出现在 $e^{-c_1M}$ 里而不是以 $M^{-1/2}$ 出现，**中等的 $M$（经验上约 100）就已足够**；速率为 $O(\Delta t)$，即**时间一阶**，与编号 93 对鞅／RDM 家族证明的一致，并优于 deep BSDE 分析的 $O(\Delta t^{1/2})$。
- 缺席的部分：常数对 $d$ 的依赖被断言是良性的，但没有量化为显式的与维数无关的陈述；$\varepsilon_k$、$M$ 与 $d$ 三者的相互作用（界需要 $M\varepsilon_k^d$ 大，任何残余的维数灾都会藏在这里）没有优化；代码里实际使用的岭参数 $\lambda$ 与 Newton 迭代都没有分析。

### 数值实验

全部实验在一台 **MacBook Pro（Apple M1 Pro，10 核，32 GB）** 上运行——论文特意强调不需要任何专用硬件。

**（一）Allen-Cahn（第 4.1 节）。** $\partial_tu+\Delta u+f(u)=0$。

- **算例 1**，双阱 $f(u)=u-u^3$，$u(T,x)=1/(2+0.4\|x\|^2)$，$T=0.3$，$d=100$，在 $x_0=0$ 处评估；参考值 $u(0,x_0)\approx0.0528$ 来自分支扩散方法（与 deep BSDE 文献和编号 96 用的是同一个基准）。**观察到的阶：绝对误差与相对误差的对数-对数斜率都接近 1，即时间一阶**，与定理 3.1 相符；改变 $M$ 只改变精度水平而不改变速率。在 $N=10^4$（$\Delta t=3\times10^{-5}$）、$M=100$ 时绝对误差约为 $1.2\times10^{-5}$。
- **算例 2**，对数势 $f(u)=\frac\theta2\ln\frac{1+u}{1-u}-\theta_cu$（$\theta<\theta_c$），人造解 $u(t,x)=\cos(\prod_{j=1}^dx_j)e^{\cos t-\|x\|^2}$ 配相应源项；$T=1$，在 $x=0$ 与 $x=(0.1,\dots,0.1)$ 处评估，$d=100$ 与 $d=1000$。两个维数都观察到一阶收敛；$d=100$ 时相对误差在 $\Delta t<6.25\times10^{-5}$（$N\ge16000$）后达到约 $10^{-3}$，之后线性下降。**值得注意的是这个 $f$ 并非全局 Lipschitz，Newton 迭代却仍然收敛得很快**。在 $N$ 与 $M$ 都大时，主导误差变成引理 3.8 的**局部回归偏差**，而不是时间步进或 Monte Carlo 噪声。运行时间（秒，$\Delta t$ 依次减半四次）：

| 维数     | 起始 $\Delta t$ | $M$ | 五次加密的运行时间                       |
| -------- | --------------- | --- | ---------------------------------------- |
| $d=100$  | 0.002           | 50  | 1.24, 2.28, 4.57, 8.83, 17.27            |
| $d=100$  | 0.002           | 100 | 3.81, 7.45, 14.99, 29.55, 59.74          |
| $d=1000$ | 0.0005          | 50  | 36.52, 72.32, 145.91, 288.43, 583.69     |
| $d=1000$ | 0.0005          | 100 | 129.24, 257.59, 521.38, 1039.55, 2082.74 |

运行时间基本随 $N\cdot M$ 线性增长。

**（二）Burgers 方程，$d=10^4$（第 4.2 节）。** $\partial_tu+(u-\frac{2+d}2)\sum_i\partial_{x_i}u+\frac{d^2}2\nu\Delta u=0$，终值 $u(T,x)=\frac{\exp(T+\sum_ix_i/d)}{1+\exp(T+\sum_ix_i/d)}$，使得 $x_0=0\in\mathbb R^{10000}$ 处 $u(0,x_0)=0.5$，$T=0.3$。**接近一阶的时间收敛**。不需要人工黏性，局部线性回归代理直接分辨出陡梯度结构——作者认为这正是显式估计 $\nabla u$ 相对深度求解器占便宜的地方。运行时间（秒，$\Delta t=0.003$ 依次减半四次）：$M=100$ 时 1040.15、2160.11、4757.37、10593.83、22174.49；$M=200$ 时 2189.28、4633.74、9674.18、21194.52、45724.85。

**（三）带梯度依赖汇项的 Hamilton-Jacobi 型方程，$d=500$ 与 $d=2000$（第 4.3 节）。** $\partial_tu+u+f(t,x,u,\nabla u)=0$，$R(u,\nabla u)=\kappa u\|\nabla u\|^2$、$\kappa=0.1$，$f=\frac{4d}{(1+4t)^{(d+2)/2}}\frac{e^{-\|x\|^2}}{1+4t}-R(u,\nabla u)$，精确解 $u(t,x)=(1+4t)^{-d/2}\exp\{-\frac{\|x\|^2}{1+4t}\}$；$T=0.5$，$N$ 最高 $3\times10^4$，在 $x=0$ 处评估。两个维数都是**一阶收敛**，Newton 每步 2 至 3 次迭代。一条报告出来的观察是：**用大约全局多项式基点数的 10% 经局部线性回归就能达到同样精度**，而 $d=2000$ 时全局多项式拟合会严重过拟合或病态。运行时间（秒，$d=500$，$\Delta t=0.0005$ 依次减半四次）：$M=50$ 时 89.07、183.69、363.86、740.39、1512.62；$M=100$ 时 185.06、376.91、765.88、1517.73、3016.35。

**论文提出的比较性主张**（作为主张核实）：在 Allen-Cahn 算例 2 上以约**少 40% 的总样本量**达到与 deep BSDE 相当的精度；代价关于 $d$ 线性，而分支扩散是 $O(d^2)$；不需要繁重的训练。

### 与其他论文的关系

**它是这个无导数三重奏里「经典数值反打一枪」的那一篇。** 编号 96 与 93 用神经网络把鞅型／无导数的想法做到 $d=10^4$–$10^5$；编号 100 保留 FBSDE／鞅型时间离散与无导数精神，但把网络整个换成**小粒子系综上的局部线性回归**，换回了可解释性、严格的先验界与笔记本量级的成本。它的引言点名引用编号 86、93、96 并把上述诊断写在明处。

**它与经典格式共享同一个 $Y$ 部分，只是条件期望的算法不同。** 那条 $\theta=1$ 的递推正是编号 8 的多步格式、编号 23 的延迟校正、编号 47 的统一格式、编号 63 的 $\theta$ 格式与编号 68 的强稳定保持格式所推广的对象。区别在于：那些论文用**张量／稀疏网格上的 Gauss-Hermite 求积**算 $\mathbb E_k$（因而止步于 $d\le6$ 上下），编号 100 用**粒子系综平均**；那些论文用第二条求积恒等式 $Z\approx\frac1{\Delta t}\mathbb E_k[Y_{k+1}\Delta W^\top]$ 得到 $Z$，编号 100 用 **$\nabla u$ 的加权最小二乘拟合**。**这一处替换就是「$d\le10$ 上的六阶」与「$d=10^4$ 上的一阶」之间的全部差别。**

与编号 25 的对照也值得一记：编号 25 用**谱稀疏网格**攻同一个条件期望瓶颈，那是确定性的答案，能上到中等维数；编号 100 是随机的答案。与 Gobet-Lemor-Warin 回归 Monte Carlo 的对照则是：全局基加耦合 $(Y,Z)$ 的 Picard 迭代，对局部基加解耦的 $X\to Z\to Y$。它与编号 93 在时间一阶这一点上殊途同归——一个经由 Euler-Maruyama 的弱阶，一个经由局部截断加 Grönwall——而两篇都明确把它与逐轨道 deep BSDE 的 $O(\Delta t^{1/2})$ 对照。

## 108：均场博弈的深度策略迭代

### 直觉

有限时程均场博弈由一个**前后向耦合的 PDE 系统**刻画：代表性个体的值函数满足向后的 HJB 方程，人口密度满足向前的 Fokker-Planck 方程。高维下这有双重困难。第一，耦合迫使反复计算对状态分布的高维积分。第二——这是本文更锋利的观察——**同一套空间采样策略必须同时**逼近那些分布积分、施加 HJB 方程、分辨 FP 方程；设计一个三件事都做得好的采样器「极不平凡」，而这一困难在已有的机器学习方法里依然存在。此外，轨道型（虚拟对局式）方法每轮外迭代都要在 $[0,T)$ 上模拟完整轨道，天然在时间上顺序；HJB 型方法则需要逐点极小化 Hamilton 量并计算二阶空间导数。

**本文的机制是「再生式改写」：把 $[0,T)$ 上的有限时程博弈重写成一个由长度为 $T$ 的循环反复拼成的无限时程再生过程。** 每个循环末尾触发一次重置。由于每个循环上的律相同，有限时程代价与均衡条件就变成**长程平均代价**与**不变占据测度**。收益是：循环下标可以与**策略迭代下标**认同，于是策略评估、策略改进与测度估计都在一个短时程上「逐循环」完成，而不是去解耦合的 HJB-FP 系统。

### 问题设定

状态取为 $x=(t,z)\in[0,T)\times\mathbb R^d$，**时间作为状态的第一个分量**，于是

$$
b(x,\mu,u)=\bigl(1,\,b_z^\top(x,\mu,u)\bigr)^\top,
\qquad
\sigma(x,\mu,u)=\bigl(0_q,\,\sigma_z^\top(x,\mu,u)\bigr)^\top,
$$

$$
J(u,\mu)=\mathbb E\Bigl[\int_0^Tf\bigl(X^{\mu,u}_s,\mu,u(X^{\mu,u}_s)\bigr)\mathrm ds
+g\bigl(X^{\mu,u}_{T-},\mu\bigr)\Bigr].
$$

**均衡**是一对 $(u^*,\mu^*)$，满足 $J(u^*,\mu^*)=\inf_uJ(u,\mu^*)$ 且 $\mu^*(A)=\frac1T\int_0^T\mathbb E[\mathbf 1(X^{\mu^*,u^*}_t\in A)]\mathrm dt$ 对一切 $A\in\mathcal B([0,T)\times\mathbb R^d)$ 成立——即 $\mu^*$ 是**时空上的占据测度**，而不是一族时间边缘。**注记 1** 指出这比标准表述略广：标准表述中 $b,\sigma,f,g$ 只通过时刻 $t$ 的边缘 $\mu_t$ 依赖 $\mu$；正是时空表述使再生式改写得以成立。

### 推导

**第一步：再生式改写。** 在第 $i$ 个循环 $t\in[iT,(i+1)T)$ 内，状态服从同一条从 $X^{\mu,u}_{iT}$ 出发的均场 SDE；每个循环末尾触发**重置机制** $X^{\mu,u}_{(i+1)T}=\xi_{i+1}$，$\xi_i\overset{\rm iid}\sim\mu_{X_0}$。由于每个循环上 $X^{\mu,u}$ 的律相同，

$$
J(u,\mu)=\lim_{t\to\infty}\frac Tt\mathbb E\Bigl[\int_0^tf\,\mathrm ds
+\sum_{i\ge1}g\bigl(X^{\mu,u}_{iT-},\mu\bigr)\mathbf 1(iT\le t)\Bigr],
\qquad
\mu^*(A)=\lim_{t\to\infty}\frac1t\int_0^t\mathbb E\bigl[\mathbf 1(X^{\mu^*,u^*}_s\in A)\bigr]\mathrm ds .
$$

**第二步：以动态规划残差做策略迭代。** 冻结 $\mu^*$ 后，值函数为 $v(x)=\inf_u\mathbb E[\int_t^Tf\,\mathrm ds+g(X_{T-},\mu^*)\mid X_t=x]$，迭代为

$$
\text{(PE)}\ \ \text{求 }v_{i+1}\in\mathcal V^{\mu^*}\ \text{使}\ \mathcal R_h(x;u_i,v_{i+1},\mu^*)=0,
\qquad
\text{(PI)}\ \ \text{求 }u_{i+1}\in\arg\min_u\mathcal R_h(x;u,v_{i+1},\mu^*),
$$

$$
\mathcal R_h(x;u,v,\mu):=\mathbb E\Bigl[\int_t^{t+h}f\bigl(X^{\mu,u}_s,\mu,u(X^{\mu,u}_s)\bigr)\mathrm ds
+v\bigl(X^{\mu,u}_{(t+h)-}\bigr)-v(x)\ \Big|\ X^{\mu,u}_t=x\Bigr].
$$

$\mathcal R_h$ 正是**鞅／动态规划残差**；令它为零就是编号 86、93、96 的离散鞅条件，只是写在一个受控、依赖均场的过程上。

**第三步：用一步随机映射代替轨道模拟。** 在均匀网格 $\Pi_h=\{nh:n=0,\dots,N-1\}$ 上（注意因为有重置，$T\notin\Pi_h$）定义**随机映射**

$$
\Phi(x,\mu,u)\overset{\mathcal L}=
\begin{cases}
x+b(x,\mu,u(x))h+\sigma(x,\mu,u(x))\sqrt h\,\zeta, & x\in\Pi_h\times\mathbb R^d,\\
X_0, & x\in\{T\}\times\mathbb R^d,
\end{cases}
\qquad \zeta\sim N(0,I_q),
$$

即循环内是 Euler-Maruyama、循环边界上是**从初始律重新抽样**。**注记 2** 给出 $\Phi$ 的弱局部截断误差为 $O(h^2)$。人口测度由一个粒子系综承载，每轮外迭代**只在一个随机小批量上推进一步**：

$$
X^m_i:=\begin{cases}\Phi(X^m_{i-1},\mu_{i-1},u_\alpha), & m\in A_i,\\ X^m_{i-1}, & m\in A_i^c,\end{cases}
\qquad
\mu_i:=\frac1M\sum_{m=1}^M\delta_{X^m_i}.
$$

**注记 3 对读懂算法很关键：$i$ 是策略迭代下标，不是时间下标。** 对每个固定的 $i$，系综 $\{X^m_i\}$ 是 $[0,T)\times\mathbb R^d$ 的一组**时空**样本，逼近占据测度 $\mu^*$，而不是某一个时间层上的样本。**这正是时间顺序瓶颈被去掉的原因：从来没有哪条轨道被端到端地模拟过。**

**第四步：弱形式加对抗训练。** 每个当前状态只有**一个**下一状态样本可用，因此条件期望靠检验去掉：

$$
\text{(PE)}\ \min_\theta\sup_{\rho\in\mathcal T}
\bigl|\mathbb E\bigl[\rho(X^m_i)\,\mathcal M(X^m_i;\mu_i,u_\alpha,v_\theta)\bigr]\bigr|^2,
\qquad
\text{(PI)}\ \min_\alpha\ \mathbb E\bigl[\mathcal M(X^m_i;\mu_i,u_\alpha,v_\theta)\bigr],
$$

其中 $\mathcal M(x;\mu,u,v):=\{hf(x,\mu,u)+v\circ\Phi(x,\mu,u)-v(x)\}\mathbf 1(t\in\Pi_h)$ 且 $\mathcal R_h=\mathbb E[\mathcal M]+O(h^2)$，$\mathcal T:=\{\rho:[0,T)\times\mathbb R^d\to[-1,1]^r\ \text{光滑}\}$。依据是塔性质：$\mathbb E[\rho\mathcal M]=\mathbb E[\rho\,\mathbb E[\mathcal M|X^m_i]]$，所以 PE 就是残差方程按 $X^m_i$ 的律加权的**Galerkin 表述**；而由全期望公式，只要 $u_\alpha$ 足够有表达力，极小化**平均**目标就使 $u_\alpha(X^m_i)$ 对几乎每个实现都极小化内层条件期望——即**平均化的优化保持逐点最优性**，论文对严格版本引用编号 86 的引理 3.2。**注记 5**：取向量值检验函数且 $r\ge600$ 相比 $r=1$ 能显著稳定对抗训练。

**第五步：网络与小批量估计。** 控制被裁剪进箱 $U=\prod_j[a_j,b_j]$：$u_{\alpha,j}(x)=\min\{\max\{a_j,\psi_{\alpha,j}(x)\},b_j\}$。值网络把终值条件硬编码，而且**用的是当前经验测度**：$[0,T)$ 上 $v_\theta=\varphi_\theta(x)$，$\{T\}\times\mathbb R^d$ 上 $v_\theta=g(x,\mu_i)$。判别器沿用编号 93 的浅层多尺度正弦网络。小批量估计同样用**不交**的 $A_{i,1},A_{i,2}$ 使平方项（近似）无偏。

**四条自陈的优点：**（一）没有高阶导数——弱形式从不对 $v_\theta$ 求 $x$ 的导数，因此没有 Hessian；（二）没有完整轨道模拟——只有一步转移，在时间维度上可并行；（三）没有逐点 Hamilton 极小化——换成 $u_\alpha$ 的平均化优化；（四）**计算共享**——同一批转移 $\{X^m_i\}_{m\in A_i}\to\{X^m_{i+1}\}$ 同时服务于策略评估、策略改进与测度更新，这恰好回答了引言里「一个采样器要干三件事」的困难。

### 定理

**本文是构造性的，没有收敛定理。** 结论一节写明：「未来工作包括建立严格的收敛性与误差估计、把框架推广到更一般的均场相互作用与公共噪声，以及发展自适应采样与网络架构。」

唯一的定量估计是局部的：**注记 2** 给出 $\Phi$ 的弱局部截断误差 $O(h^2)$；以及 $\mathcal R_h=\mathbb E[\mathcal M]+O(h^2)$，其中 $O(h^2)$ 来自左矩形求积加弱 Euler 误差。没有全局速率，没有按 $M$、$|A_i|$、网络规模或策略迭代次数给出的界，也没有证明这个类虚拟对局的外循环收敛到均衡。**注记 4** 给出定位：该迭代是一种**虚拟对局**——第 $i$ 轮里 $u_\alpha$ 对上一轮继承来的人口测度做最优应对——与已有的虚拟对局型均场博弈求解器的区别在于，测度是用**一步 Markov 转移**推进的，而不是重新模拟整条轨道。

### 数值实验

配置：$\psi_\alpha,\varphi_\theta$ 为全连接，深度 $H=6$，$d=1$ 时宽 $W=104$、$d=1000$ 时宽 $W=1008$，ReLU；判别器 $r=1200$、$c=10$；学习率 $\delta_1=\delta_2=\delta_0\times10^{-3}\times0.01^{i/I}$、$\delta_3=\delta_0\times10^{-2}\times0.01^{i/I}$，$\delta_0=3d^{-0.5}$；$I=9000$ 次迭代；$J=2K=2$；全批 $M=1024\times10^3$，小批量 $|A_{i,1}|=|A_{i,2}|=M/20$；$h=T/100$；RMSProp；PyTorch 2.6.0，float32 加自动混合精度，分布式数据并行运行在 **8 块 NVIDIA RTX 4090** 上。$Z_0\sim\mathrm{Uniform}(\{s\mathbf 1_d:s\in[-c,c]\})$，第 4.1 节取 $c=1$、第 4.2 节取 $c=3$。指标为相对代价 $\mathrm{RC}=|\hat J-J^*|/|J^*|$ 以及值函数在 1000 个固定测试点上的 $\mathrm{RE}_1$ 与 $\mathrm{RE}_\infty$。

**（一）有显式解的线性二次均场博弈（第 4.1 节，表 1）。** $b_z=u+c_0(\bar m_t-z)+c_1(z^*(t)-\bar m_t)$，$\sigma_z=c_\sigma I_d$，$f=\frac12\{c_2|u|^2+c_3|z-\bar m_t|^2+c_4|z-z^*(t)|^2\}$，$g=\frac{c_5}2|z-z^*(T)|^2+\frac12$，$T=1$，$\bar m_t$ 为人口均值。参考的 $v$ 与 $u^*$ 由 Riccati 型 ODE 用 `scipy.integrate.solve_ivp` 的 RK45、rtol $10^{-8}$ 解出（附录给出推导）。三个变体：**LQ-1**（$c_0=c_1=c_4=0$，$c_2=1$，$c_3=c_5=1/d$，$c_\sigma=0.5/\sqrt d$，无目标跟踪）；**LQ-2**（$c_0=1$，$c_1=0$，$c_2=c_3=c_4=c_5=1/d$，目标为单位圆 $z^*(t)=y^*(t)/|y^*(t)|$，$y^*_i(t)=\sin(2\pi t+i\pi/2)$）；**LQ-3**（同 LQ-2 但 $c_1=-0.5$，目标为螺旋线 $z^*(t)=2t\,y^*(t)/|y^*(t)|$）。初始猜测刻意取得很差：$\hat\mu$ 取 $(t,Z_{{\rm init},t})$ 的律，$Z_{{\rm init},t}=Z_0+5t\mathbf 1_d+c_\sigma B_t$，「有意地远离真正的 $\mu^*$」。下表为 5 次运行的均值（括号内为标准差），MEM 为每块 GPU 的峰值内存（另有 $\times8$）。

| 方程 | $d$  | $\mathrm{RE}_1$   | $\mathrm{RE}_\infty$ | RC                | MEM (MB) | RT (s) |
| ---- | ---- | ----------------- | -------------------- | ----------------- | -------- | ------ |
| LQ-1 | 1    | 4.05E-3 (2.09E-3) | 1.26E-2 (5.20E-3)    | 9.69E-3 (6.18E-3) | 535      | 294    |
| LQ-1 | 1000 | 1.68E-2 (2.13E-3) | 4.33E-2 (5.99E-3)    | 1.55E-2 (1.26E-3) | 3964     | 579    |
| LQ-2 | 1    | 9.92E-3 (2.41E-3) | 2.11E-2 (3.91E-3)    | 9.84E-3 (3.24E-3) | 535      | 308    |
| LQ-2 | 1000 | 1.97E-2 (4.57E-3) | 2.86E-2 (5.18E-3)    | 1.95E-2 (3.12E-3) | 3995     | 602    |
| LQ-3 | 1    | 6.67E-3 (5.62E-3) | 1.67E-2 (7.65E-3)    | 6.54E-3 (5.23E-3) | 535      | 317    |
| LQ-3 | 1000 | 1.75E-2 (2.15E-3) | 2.88E-2 (5.22E-3)    | 1.79E-2 (3.01E-3) | 3997     | 615    |

**这张表建立的是维数依赖的温和性：$d$ 从 1 涨到 1000，$\mathrm{RE}_1$、$\mathrm{RE}_\infty$、RC 只分别涨了约 2.0-4.1、1.4-3.4、1.6-2.7 倍，而运行时间只涨约 2 倍、内存涨约 7.5 倍。** 另有一次 **$d=10\,000$ 的 LQ-1 运行**：$v_\theta$ 沿对角线与 $v$ 吻合良好，RE 与 RC 稳定下降，**运行时间 5258 秒，每块 GPU 峰值内存 55,645 MB**（$\times8$）。$d=1$ 的可视化还显示：值网络的误差在模拟轨道访问过的区域很小、在未探索区域较大——**与编号 93 用加权范数 $M^2_n[\cdot]$ 证明的是同一个现象。**

**（二）银行间系统性风险均场博弈（第 4.2 节）。** Carmona-Fouque-Sun 式模型，为与本格式相容而**去掉公共噪声**：$b_z=c_1(\bar m_t-z)+u$，$\sigma_z=c_2$，运行代价 $f=\frac12u^2-c_3u(\bar m_t-z)+\frac{c_4}2(\bar m_t-z)^2$，终值 $g=\frac{c_5}2(\bar m_T-z)^2$，$T=1$，$c_1=1$、$c_2=0.5$、$c_3=1$、$c_4=2$、$c_5=1$，标量状态。参考解同样来自 Riccati ODE。报告的诊断量是 RE 与 RC 随迭代的曲线，以及学到的经验平均储备 $\hat m_t$ 与精确 $m^*_t$ 的对比。

**（三）带非线性相互作用的目标跟踪均场博弈，$d=2$（第 4.3 节）。** 为便于可视化而取低维：$b_z=u$，$\sigma_z=0.5I_d$，$f=c_1|z-z^*(t)|^2+c_2|u|^2+F(x,\mu)$，$g=0.5c_4|z-z^*(T)|^2$，圆形目标 $z^*_i(t)=\sin(2\pi t+i\pi/2)$，以及**拥堵／反聚集项**

$$
F(x,\mu):=\int_{[0,T)\times\mathbb R^d}\delta_t(s)\exp\bigl(-c_F|z-y|^2\bigr)\,\mu(\mathrm ds\times\mathrm dy),
$$

$c_1=1$、$c_2=0.1$、$c_4=10$、$c_F=1$。没有显式解，因此只做定性验证：经验代价在训练中下降约 0.1，个体轨迹描出近似圆形并终止于目标附近。

**（四）带非线性动力学与相互作用的均场博弈，$d=2$（第 4.4 节）。** $b_z=u$，$\sigma_z=\frac{c_1}{\sqrt d}\sum_{i=1}^d\sin(i+z_i)I_d$（状态依赖的非常数扩散），$f=c_0/(1+c_2|z-0.5\mathbf 1_d|^2)+c_3|u|^2+F(x,\mu)$，$g=0.5c_4|z-z^*(T)|^2$，$T=1$，$c_1=c_3=0.1$、$c_2=c_4=1$、$c_0=5$、$c_F=1$，$Z_0$ 为标准 Gauss。$f$ 的第一项在 $0.5\mathbf 1_d$ 处造出一道**壁垒**，观察到个体散开、绕过壁垒、汇聚到目标。

### 与其他论文的关系

**它是 Cai-Fang-Zhou 鞅／对抗家族在均场博弈方向的成员。** 机器几乎原样继承自编号 86、93、96：鞅／动态规划残差 $\mathcal M$、带对抗检验网络的 Galerkin 弱形式、输出维数很大的多尺度正弦判别器、用于平方期望无偏化的不交小批量 $A_{i,1},A_{i,2}$、箱裁剪的控制网络，以及硬编码进 $v_\theta$ 的终值条件。注记 6 归功于弱对抗网络的想法，平均化策略改进步的逐点最优性论证则引用编号 86 的引理 3.2。

**相对编号 86、93、96 真正新的东西**是**再生式改写**（循环加重置，把有限时程均场博弈变成带不变占据测度的无限时程问题）与粒子系综的**一步随机映射**更新。两者合起来使外循环下标兼作循环下标，于是没有任何完整轨道被模拟过——这是编号 96 与 93 的「离线生成轨道」技巧的均场博弈版本，且推得更远。

**均场方向上的联系：** 编号 33（均场 BSDE 的显式 $\theta$ 格式）与编号 61（均场 FBSDE 的显式多步格式）是均场**后向**方程的高阶、低维、可证收敛的处理；编号 97（DeepSPoC）处理均场**前向**（McKean-Vlasov／Fokker-Planck）一侧；编号 108 把两侧耦合起来求均衡，$d$ 做到 $10^4$ 但没有理论。**编号 33／61／97／108 张成的正是编号 8／47／63 与编号 86／93／96 在非均场情形下张成的同一张权衡表：阶与证明，对维数。**

与随机控制线索的联系：编号 26、41、50 走 Pontryagin／FBSDE 路线（梯度投影、高阶 FBSDE 格式、带粒子滤波的条件梯度），编号 108 走动态规划与策略迭代，并额外加上均衡不动点。人员上，方水鑫与周涛同时是编号 86、93、96、100 的作者；吴臻是山东大学的控制论学者，而山东大学也是赵卫东的单位——赵卫东是编号 8、16、18、19、23、25、26、33、35、41、47、61、63、68 的经典 FBSDE 合作者。

## 这条路线的一般判断

鞅型方法改变的是**残差的检验方式**。在高维空间中「逐点检验残差」不可行，因为没有网格；而「检验一个过程是否是鞅」只需要沿模拟轨道的期望，因此维数不再直接进入代价。

编号 86 把这一转换做得最完整，可以概括成三次替换：

1. 逐点极小原理 → 积分极小原理（去掉 $(t,x)$ 上的维数灾）；
2. PDE 残差 → 鞅性质（去掉网格）；
3. 条件期望 → 对检验函数族的弱条件（去掉条件期望的计算）。

每一次替换都把一个难以直接处理的条件换成一个可以用采样验证的等价条件。**第三次替换尤其值得注意**：它把「计算条件期望」变成「训练一个判别网络」，代价是引入对抗训练的全部不稳定性，收益是完全避开了条件期望的估计。

后五篇各自去掉了这套机器的一件零件，而这个「去掉一件」的谱系本身很能说明问题：

| 论文 | 去掉的是                           | 换来的                                      | 代价                                |
| ---- | ---------------------------------- | ------------------------------------------- | ----------------------------------- |
| 96   | 损失里的所有导数与轨道重模拟       | $d=10^4$，时空双向并行                      | 仍无全局定理；受控情形跳跃不能预算  |
| 93   | 随机分析（改用 Taylor 与矩恒等式） | 第一条收敛率定理（时间一阶），并反证 86、96 | 精度低于 PINN；RDO 自带 $O(h)$ 截断 |
| 100  | 神经网络本身                       | 严格先验界，笔记本可跑，可解释              | 只到 $d=10^4$，只有一阶             |
| 97   | 存储的粒子轨道                     | 前向均场方程，且给出唯一的可计算后验界      | 收敛定理只在 Fourier 替身里成立     |
| 108  | 完整轨道模拟与耦合 HJB-FP 求解     | 均场博弈均衡，$d=10^4$                      | 没有任何收敛定理                    |

**一条横贯全表的限制值得单独记住**：这一族方法的判据只把 PDE 残差压到零**于引导过程／采样系综所探索到的区域内**。编号 96 的注记 1 明说了这一点，编号 93 的定理 4 把它编码在加权范数 $M^2_n[\cdot]$ 里，编号 108 在 $d=1$ 的可视化里直接展示了它。**「解得准」在这里始终是一句有定义域的话。**

另一条判断关于精度的量级。经典多步／$\theta$ 格式在 $d\lesssim10$ 上做到六阶，本页各篇在 $d=10^4$ 上做到一阶、相对误差 $10^{-2}$ 到 $10^{-3}$。**这不是同一件事做得好坏的差别，而是两件不同的事。** 编号 93 的表 4 尤其诚实：它的方法比 PINN 精度更低，全部收益在时间与内存上，而且收益随 $d$ 增长——这正是一个「用精度换可行性」的交易被明码标价的例子。

## 覆盖核对

| 内容                                      | 论文 | 覆盖状态                       |
| ----------------------------------------- | ---- | ------------------------------ |
| HJB 方程、Hamilton 函数与四种覆盖情形     | 86   | 全文核实                       |
| 积分极小原理（引理 3.2）与其一行证明      | 86   | 全文核实                       |
| 鞅刻画（引理 3.4）与双向证明              | 86   | 全文核实                       |
| 定理 3.5、几何假设、弱阶对强阶            | 86   | 全文核实；已注明它不是收敛定理 |
| 对抗弱形式、增广 Lagrange、三个网络       | 86   | 全文核实                       |
| 实验：运行时间、一阶率、表 1、表 2        | 86   | 全文核实，已列表               |
| 随机差分展开、RDO、Fokker-Planck 采样     | 93   | 全文核实                       |
| 与鞅方法及隐式 Euler 的两条等值关系       | 93   | 全文核实                       |
| 假设 1-3、定理 1、推论 2、引理 3、定理 4  | 93   | 全文核实，含步长限制与常数     |
| 注记 9（同时证明 86、96 的一阶收敛）      | 93   | 全文核实                       |
| 实验：表 2、表 3、表 4 与作者自评         | 93   | 全文核实，已列表               |
| 引导过程／系统过程、局部化鞅、注记 1      | 96   | 全文核实；区域限制已保留       |
| 无导数 $\mathcal M$、不交小批量、PIA      | 96   | 全文核实                       |
| 只有 $O(h^2)$ 局部估计，无全局定理        | 96   | 全文核实                       |
| 实验：表 1、表 2、宽度研究                | 96   | 全文核实，已列表               |
| SPoC、deepSPoC 系统、三种损失、自适应采样 | 97   | 全文核实                       |
| 定理 3.3 及其替身设定的保留               | 97   | 全文核实                       |
| 后验估计（命题 3.5、定理 3.6）            | 97   | 全文核实                       |
| 实验：五类方程与各自参数                  | 97   | 全文核实；无阶的报告，已注明   |
| 局部线性回归、免矩阵求解、Newton 迭代     | 100  | 全文核实                       |
| 引理 3.1、3.3、3.8 与定理 3.1             | 100  | 全文核实                       |
| 实验：三类方程、运行时间表、硬件          | 100  | 全文核实，已列表               |
| 再生式改写、一步随机映射、注记 3          | 108  | 全文核实                       |
| 无收敛定理，只有 $O(h^2)$ 局部估计        | 108  | 全文核实                       |
| 实验：表 1、$d=10^4$ 运行、四类算例       | 108  | 全文核实，已列表               |

## 本页原文

- W. Cai, S. Fang, and T. Zhou, [_SOC-MartNet: a martingale neural network for the Hamilton-Jacobi-Bellman equation without explicit inf H in stochastic optimal controls_](https://doi.org/10.1137/24M1681033), SIAM J. Sci. Comput. 47(4) (2025), pp. C795-C819（预印本 [arXiv:2405.03169](https://arxiv.org/abs/2405.03169)；代码 [sx-fang/MartNet](https://github.com/sx-fang/MartNet)）。
- W. Cai, S. Fang, and T. Zhou, [_Deep random difference method for high-dimensional quasilinear parabolic partial differential equations_](https://doi.org/10.1016/j.jcp.2026.114767), J. Comput. Phys. 555 (2026), 114767（预印本 [arXiv:2506.20308](https://arxiv.org/abs/2506.20308)；代码 [sx-fang/DRDM](https://github.com/sx-fang/DRDM)）。
- W. Cai, S. Fang, W. Zhang, and T. Zhou, _Martingale deep learning for very high dimensional quasi-linear partial differential equations and stochastic optimal controls_, [arXiv:2408.14395](https://arxiv.org/abs/2408.14395)，投稿 SIAM Rev.
- K. Du, Y. Xie, T. Zhou, and Y. Zhou, _DeepSPoC: a deep learning based sequential propagation of chaos_，[arXiv:2408.16403](https://arxiv.org/abs/2408.16403)，投稿 SIAM/ASA J. Uncertain. Quantif.
- S. Fang, C. Sheng, B. Su, and T. Zhou, _A derivative-free localized stochastic method for very high dimensional semi-linear parabolic PDEs_，[arXiv:2510.02635](https://arxiv.org/abs/2510.02635)，投稿 Numer. Math.
- S. Fang, S. Wang, Z. Wu, H. Zhang, and T. Zhou, _Deep policy iteration for high-dimensional mean-field games with regenerative reformulation_，预印本 arXiv:2604.26782。
