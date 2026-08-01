---
title: 3.5.2：迭代 ParaDiag（ParaDiag-II）
description: 从循环预条件、首尾耦合波形松弛到谱界、舍入误差和非线性 Newton–Krylov 的完整推导
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-5-paradiag-ii
tags:
  - 时间并行
  - ParaDiag
  - 循环预条件
---

> [!note] 阅读范围
> 本页对应论文 Section 3.5.2（pp. 431–443），覆盖公式 (3.49)–(3.68)、Theorems 3.8–3.9 和 Figures 3.15–3.18。文中从离散循环预条件和连续首尾耦合两条路线出发，随后证明两条路线落到同一个 ParaDiag-II 迭代。

## 3.5.2.1 从多步公式到循环预条件

ParaDiag-I 很难直接推广到多级 Runge–Kutta 等高阶积分器。ParaDiag-II 改为近似时间步进矩阵，并在定常迭代或 Krylov 方法中求解全时间系统。McDonald et al. (2018) 从离散系统出发；Gander and Wu (2019) 同期从连续波形松弛出发。两种表述在离散后等价。

先用 $m$ 步线性多步法离散 (2.1)：

$$
\sum_{l=0}^{m}a_l\boldsymbol u_{n-l}
=\Delta t\sum_{l=0}^{m}b_lA\boldsymbol u_{n-l}
+\bar{\boldsymbol g}_n,
\qquad n=1,\ldots,N_t.
$$

已知初始历史 $\{\boldsymbol u_{-(m-1)},\ldots,\boldsymbol u_0\}$。令

$$
\boldsymbol U=(\boldsymbol u_1^\top,\ldots,\boldsymbol u_{N_t}^\top)^\top,
\qquad
K=B_1\otimes I_x-B_2\otimes(\Delta t A),
$$

其中 $B_1,B_2$ 是由 $\{a_l\}$ 和 $\{b_l\}$ 生成的下三角带状 Toeplitz 矩阵。把两者换成 Strang 型循环矩阵 $C_1,C_2$，得到

$$
P=C_1\otimes I_x-C_2\otimes(\Delta t A).
$$

若 $P$ 已足够接近 $K$，可直接使用定常迭代

$$
P\Delta\boldsymbol U^k=\boldsymbol r^k
:=\boldsymbol b-K\boldsymbol U^k,
\qquad
\boldsymbol U^{k+1}=\boldsymbol U^k+\Delta\boldsymbol U^k,
\qquad k=0,1,\ldots. \tag{3.49}
$$

其渐近速度由 $\rho(P^{-1}K)$ 控制。即使谱半径不小，$P^{-1}K$ 的谱若高度聚集，预条件 GMRES 仍可能很快。

## 3.5.2.2 Fourier 对角化与三步并行求解

循环矩阵彼此交换，可由同一个离散 Fourier 矩阵同时对角化：

$$
C_l=F^*D_lF,\qquad l=1,2,
$$

$$
F=\frac1{\sqrt{N_t}}
\begin{bmatrix}
1&1&\cdots&1\\
1&\omega&\cdots&\omega^{N_t-1}\\
\vdots&\vdots&&\vdots\\
1&\omega^{N_t-1}&\cdots&\omega^{(N_t-1)^2}
\end{bmatrix},
\qquad
\omega=\exp\!\left(\frac{2\pi i}{N_t}\right). \tag{3.50}
$$

若 $C_l(:,1)$ 表示第一列，则

$$
D_l=\operatorname{diag}\!\left(\sqrt{N_t},F C_l(:,1)\right),
\qquad l=1,2. \tag{3.51}
$$

于是

$$
P=(F^*\otimes I_x)
\left(D_1\otimes I_x-D_2\otimes(\Delta t A)\right)
(F\otimes I_x),
$$

一次预条件作用分为

$$
\left\{
\begin{aligned}
\boldsymbol U^a&=(F\otimes I_x)\boldsymbol r^k,\\
(\lambda_{1,n}I_x-\lambda_{2,n}\Delta t A)\boldsymbol u_n^b
&=\boldsymbol u_n^a,
&&n=1,\ldots,N_t,\\
\boldsymbol U&=(F^*\otimes I_x)\boldsymbol U^b.
\end{aligned}
\right. \tag{3.52}
$$

首尾两步由 FFT 完成，复杂度为 $O(N_xN_t\log N_t)$；中间的 $N_t$ 个空间系统相互独立。

### Theorem 3.8：精确聚集的上界

若 $A\in\mathbb R^{N_x\times N_x}$ 对称负定，则 $P^{-1}K$ 至多有 $mN_x$ 个特征值不等于 $1$。因此精确算术下，预条件 GMRES 至多 $mN_x+1$ 步终止。该结论对大 $N_x$ 并未给出很小的迭代上界；$A$ 失去对称性后，谱聚集也会变弱。

## 3.5.2.3 Figure 3.15：三类 PDE 上的谱与 GMRES

论文比较热方程、两种黏性的对流扩散方程，以及波动方程。全部使用齐次 Dirichlet 边界和 $u(x,0)=\sin(2\pi x)$；波动方程另取 $u_t(x,0)=0$。空间中心差分后，波动方程为

$$
\boldsymbol u''(t)=A\boldsymbol u(t),
\qquad
\boldsymbol u(0)=\boldsymbol u_0,
\quad \boldsymbol u'(0)=0,
\quad t\in(0,T], \tag{3.53}
$$

其中 $A=\operatorname{Tri}[1,-2,1]/\Delta x^2$。一阶问题使用梯形规则，二阶问题使用带参数的 Numerov 型方法：

$$
\left\{
\begin{aligned}
\widetilde{\boldsymbol u}_n-\boldsymbol u_n
+\gamma\Delta t^2A
(\boldsymbol u_{n+1}-2\boldsymbol u_n+\boldsymbol u_{n-1})&=0,\\
\boldsymbol u_{n+1}-2\boldsymbol u_n+\boldsymbol u_{n-1}
-\frac{\Delta t^2A}{12}
(\boldsymbol u_{n+1}+10\widetilde{\boldsymbol u}_n+\boldsymbol u_n)&=0.
\end{aligned}
\right. \tag{3.54}
$$

$\gamma=0$ 给出经典 Numerov 四阶方法，但仅条件稳定；$\gamma\ge1/120$ 时保持四阶并达到无条件稳定。

![原论文 Figure 3.15：三类 PDE 上循环预条件后的谱和 GMRES 残差](assets/papers/time-parallelization/source-figures/figure-3-15.svg)

Figure 3.15 使用 $T=2$、$\Delta t=1/50$、$\Delta x=1/100$、$\gamma=1/100$。热方程的特征值几乎全聚在 $1$；对流增强后聚集逐渐变松，GMRES 需要更多轮；切到波动方程后预条件效果明显下降。黏性 $\nu=10^{-3}$ 时仍可工作，再减小黏性便会持续恶化。

## 3.5.2.4 连续首尾耦合与 $\alpha$-循环矩阵

Gander and Wu (2019) 在连续层面引入

$$
\boldsymbol u_t^k(t)=A\boldsymbol u^k(t)+\boldsymbol g(t),
\qquad
\boldsymbol u^k(0)=
\alpha[\boldsymbol u^k(T)-\boldsymbol u^{k-1}(T)]+\boldsymbol u_0, \tag{3.55}
$$

其中 $\alpha\in\mathbb C$。收敛后括号内的差消失，初值回到 $\boldsymbol u_0$。

用一类单步公式离散：

$$
\left\{
\begin{aligned}
r_1(\Delta t A)\boldsymbol u_n^k
&=r_2(\Delta t A)\boldsymbol u_{n-1}^k+\widetilde{\boldsymbol g}_n,
&&n=1,\ldots,N_t,\\
\boldsymbol u_0^k
&=\alpha(\boldsymbol u_{N_t}^k-\boldsymbol u_{N_t}^{k-1})+\boldsymbol u_0.
\end{aligned}
\right. \tag{3.56}
$$

后向 Euler 与梯形规则分别对应

$$
\left\{
\begin{aligned}
r_1&=I_x-\Delta tA,&r_2&=I_x,
&&\text{后向 Euler},\\
r_1&=I_x-\tfrac12\Delta tA,&r_2&=I_x+\tfrac12\Delta tA,
&&\text{梯形规则}.
\end{aligned}
\right. \tag{3.57}
$$

把首尾条件代入第一步，可写成

$$
P_\alpha\boldsymbol U^k=\boldsymbol b^k, \tag{3.58a}
$$

$$
\boldsymbol b^k=\boldsymbol b-
\alpha
\begin{bmatrix}
r_2(\Delta tA)\boldsymbol u_{N_t}^{k-1}\\0\\\vdots\\0
\end{bmatrix},
\qquad
\boldsymbol b=
\begin{bmatrix}
r_2(\Delta tA)\boldsymbol u_0+\widetilde{\boldsymbol g}_1\\
\widetilde{\boldsymbol g}_2\\\vdots\\
\widetilde{\boldsymbol g}_{N_t}
\end{bmatrix}. \tag{3.58b}
$$

$$
P_\alpha=
\begin{bmatrix}
r_1& & &-\alpha r_2\\
-r_2&r_1\\
&\ddots&\ddots\\
&&-r_2&r_1
\end{bmatrix}
=I_t\otimes r_1(\Delta tA)-C_\alpha\otimes r_2(\Delta tA), \tag{3.58c}
$$

$$
C_\alpha=
\begin{bmatrix}
0&&&\alpha\\
1&0\\
&\ddots&\ddots\\
&&1&0
\end{bmatrix}.
$$

对 Strang 型 $\alpha$-循环矩阵，

$$
C_\alpha=V_\alpha D_\alpha V_\alpha^{-1}, \tag{3.59a}
$$

$$
D_\alpha=\operatorname{diag}\!\left(\sqrt{N_t},F\Lambda_\alpha C_\alpha(:,1)\right),
\quad
V_\alpha=\Lambda_\alpha F^*,
\quad
\Lambda_\alpha=\operatorname{diag}
(1,\alpha^{-1/N_t},\ldots,\alpha^{-(N_t-1)/N_t}). \tag{3.59b}
$$

因此求解 (3.58a) 仍是三步：

$$
\left\{
\begin{aligned}
\boldsymbol U^a&=(V_\alpha^{-1}\otimes I_x)\boldsymbol b^k,\\
\left(r_1(\Delta tA)-\lambda_n r_2(\Delta tA)\right)\boldsymbol u_n^b
&=\boldsymbol u_n^a,
&&n=1,\ldots,N_t,\\
\boldsymbol U^k&=(V_\alpha\otimes I_x)\boldsymbol U^b.
\end{aligned}
\right. \tag{3.60}
$$

$\alpha=1$ 时 $V_\alpha=F^*$，该算法与 (3.52) 完全一致。$\alpha\ne1$ 时，$V_\alpha$ 只是 Fourier 变换前后多乘一个对角缩放，FFT 仍可使用。

![原论文 Figure 3.16：两种 alpha 下首尾耦合波形松弛的误差衰减](assets/papers/time-parallelization/source-figures/figure-3-16.svg)

Figure 3.16 沿用 Figure 3.15 的 PDE 参数并改为周期边界。此时 $\alpha=1$ 的 $P$ 奇异，需要取 $\alpha<1$。较小的 $\alpha$ 使强对流和波动问题也能在没有 Krylov 加速的情况下迅速收敛。

## 3.5.2.5 两条路线的离散等价

首尾耦合离散可以写成预条件定常迭代

$$
P_\alpha\Delta\boldsymbol U^{k-1}
=\boldsymbol r^{k-1}:=\boldsymbol b-K\boldsymbol U^{k-1},
\qquad
\boldsymbol U^k=\boldsymbol U^{k-1}+\Delta\boldsymbol U^{k-1}. \tag{3.61}
$$

这里

$$
K=I_t\otimes r_1(\Delta tA)-B\otimes r_2(\Delta tA), \tag{3.62a}
$$

$$
B=
\begin{bmatrix}
0\\1&0\\&\ddots&\ddots\\&&1&0
\end{bmatrix}. \tag{3.62b}
$$

由 $\boldsymbol b^k=(P_\alpha-K)\boldsymbol U^{k-1}+\boldsymbol b$，(3.58a) 立即化为 (3.61)。$\alpha=1$ 时，它就是 McDonald et al. 的 (3.49)；$0<\alpha<1$ 时，对应 Banjai and Peterseim (2012) 的并行方法。

一般单步公式

$$
r_1(\Delta tA)\boldsymbol u_n
=r_2(\Delta tA)\boldsymbol u_{n-1}+\widetilde{\boldsymbol g}_n,
\qquad n=1,\ldots,N_t, \tag{3.63}
$$

生成全时间方程 $K\boldsymbol U=\boldsymbol b$，$P_\alpha$ 正是其广义块 $\alpha$-循环预条件器。

## 3.5.2.6 直接离散二阶系统

把二阶系统改写为一阶会让每个时间点的存储翻倍。论文因此考虑对称两步格式

$$
r_1(\Delta t^2A)\boldsymbol u_{n+1}
-r_2(\Delta t^2A)\boldsymbol u_n
+r_1(\Delta t^2A)\boldsymbol u_{n-1}
=\widetilde{\boldsymbol g}_n,
\quad n=1,\ldots,N_t-1, \tag{3.64}
$$

并假设第二个初值 $\boldsymbol u_1$ 已知。Numerov 型方法对应

$$
r_1(z)=I_x-\frac{z}{12}+\frac{10\gamma z^2}{12},
\qquad
r_2(z)=2I_x+\frac{10z}{12}+\frac{20\gamma z^2}{12},
\quad z=\Delta t^2A.
$$

全时间矩阵及预条件器为

$$
K=\widetilde B\otimes r_1(\Delta t^2A)-B\otimes r_2(\Delta t^2A),
\qquad
P_\alpha=\widetilde C_\alpha\otimes r_1(\Delta t^2A)
-C_\alpha\otimes r_2(\Delta t^2A), \tag{3.65a}
$$

其中 $B,C_\alpha$ 如上，

$$
\widetilde B=
\begin{bmatrix}
1\\0&1\\1&0&1\\&\ddots&\ddots&\ddots\\&&1&0&1
\end{bmatrix},
\qquad
\widetilde C_\alpha=
\begin{bmatrix}
1&&\alpha\\0&1&&\alpha\\1&0&1\\&\ddots&\ddots&\ddots\\&&1&0&1
\end{bmatrix}. \tag{3.65b}
$$

$C_\alpha$ 与 $\widetilde C_\alpha$ 可同时对角化，所以 (3.60) 的三步求解也适用于二阶预条件器。

## 3.5.2.7 Theorem 3.9：稳定性、谱界与收敛因子

对于一阶问题，如果 (3.63) 稳定，即

$$
|r_1^{-1}(z)r_2(z)|\le1,
\qquad z\in\sigma(\Delta tA)\subset\mathbb C_-,
$$

论文印刷版给出

$$
\frac1{1-\alpha}
\le |\lambda(P_\alpha^{-1}K)|
\le\frac1{1+\alpha},
\qquad 0<\alpha<1. \tag{3.66, as printed}
$$

二阶两步格式若满足 $|r_1^{-1}(z)r_2(z)|\le2$，且等号只在 $z=0$ 取得，也声称满足同一界。

> [!warning] 原文公式核对
> 对 $0<\alpha<1$，印刷式左端大于右端，区间为空。这里完整保留论文排版，并把疑似校正记为 $1/(1+\alpha)\le|\lambda|\le1/(1-\alpha)$。后文给出的迭代矩阵界也与“分母交换后的次序”相容。此处是读者核对后的注释，不改写原论文陈述。

定常迭代矩阵为

$$
\mathcal M=I-P_\alpha^{-1}K, \tag{3.67}
$$

论文据此得到 $\rho(\mathcal M)\le\alpha/(1-\alpha)$。小 $\alpha$ 加快收敛，Figure 3.16 的实验与此吻合。底层积分器稳定是谱界成立的充分条件；论文的数值结果也表明它接近必要条件。

![原论文 Figure 3.17：稳定阈值两侧的 Numerov 型方法迭代谱](assets/papers/time-parallelization/source-figures/figure-3-17.svg)

Figure 3.17 取 $\Delta t=1/16$、$\Delta x=1/128$、$\alpha=0.02$，并比较 $T=0.5,10,20$。上排 $\gamma=1/120$ 位于无条件稳定阈值，所有特征值都在半径 $\alpha/(1-\alpha)$ 的虚线圆内；下排 $\gamma=1/120.01$ 略低于阈值，长时间窗下谱界失效。

## 3.5.2.8 小 $\alpha$ 的舍入限制与实现方式

$\alpha$ 不能无限减小。由 $V_\alpha=\Lambda_\alpha F^*$、$\operatorname{Cond}_2(F^*)=1$，

$$
\operatorname{Cond}_2(V_\alpha)=\frac1\alpha,
\qquad
\mathrm{err}_{\mathrm{ro}}
=O\!\left(\epsilon\operatorname{Cond}_2(V_\alpha)\right)
=O\!\left(\frac\epsilon\alpha\right).
$$

这一矩阵分解误差并不必然等量传入最终迭代。直接由 (3.58) 求 $\boldsymbol U^k$ 会显露 $O(\epsilon/\alpha)$ 的放大；先由 (3.61) 求误差增量 $\Delta\boldsymbol U^{k-1}$，再更新解，能显著减轻舍入误差。

![原论文 Figure 3.18：三种 alpha 下两种实现方式的舍入误差](assets/papers/time-parallelization/source-figures/figure-3-18.svg)

Figure 3.18 的左图对应直接解当前迭代，误差随 $\alpha$ 减小更早停滞；右图对应增量更新，能够下降到接近双精度极限。论文把更系统的舍入分析指向 Wu, Yang and Zhou (2025)。

对于一般多步法，(3.66) 未必成立。Volterra 偏积分微分方程会产生稠密下三角 Toeplitz 时间矩阵；若积分权重正且单调，已有结果给出 $|\lambda(P_\alpha^{-1}K)|=1+O(\alpha)$。

## 3.5.2.9 非线性 Newton–Krylov

对 $\boldsymbol u'=f(t,\boldsymbol u)$ 使用后向 Euler，先在非线性全时间方程上做 Newton：

$$
J\Delta\boldsymbol U^l=\boldsymbol b-F(\boldsymbol U^l),
\qquad
\boldsymbol U^{l+1}=\boldsymbol U^l+\Delta\boldsymbol U^l. \tag{3.68}
$$

其中

$$
J=B\otimes I_x-\nabla F_l,
\qquad
\nabla F_l=\operatorname{blkdiag}\!\left(
\nabla f(\boldsymbol u_1^l,t_1),\ldots,
\nabla f(\boldsymbol u_{N_t}^l,t_{N_t})
\right),
$$

$$
B=\frac1{\Delta t}
\begin{bmatrix}
1\\-1&1\\&\ddots&\ddots\\&&-1&1
\end{bmatrix}.
$$

内层 GMRES 使用

$$
P_\alpha=C_\alpha\otimes I_x-I_t\otimes A_l,
$$

$A_l$ 为所有 $\nabla f(\boldsymbol u_n^l,t_n)$ 的平均。一般有 $\rho(P_\alpha^{-1}J)>1$，所以定常迭代 (3.61) 不适合该 Jacobian 系统；预条件谱仍会聚集，GMRES 因而有效。缩短时间窗会增强聚集并降低内层迭代数，Section 3.5.1 的最近 Kronecker 近似还可进一步改善预条件。

## 公式与图表覆盖核对

| 原文项目                                             | 本页位置    | 覆盖状态                                             |
| ---------------------------------------------------- | ----------- | ---------------------------------------------------- |
| 多步全时间系统、循环替换、(3.49)–(3.52), Theorem 3.8 | §§3.5.2.1–2 | 离散起点、FFT 对角化、三步求解、谱聚集               |
| (3.53)–(3.54), Figure 3.15                           | §3.5.2.3    | 三类 PDE、Numerov 稳定阈值、谱与 GMRES               |
| (3.55)–(3.60), Figure 3.16                           | §3.5.2.4    | 连续首尾条件、$\alpha$-循环矩阵、FFT 实现与实验      |
| (3.61)–(3.65)                                        | §§3.5.2.5–6 | 两条路线等价、一阶与直接二阶全时间系统               |
| Theorem 3.9, (3.66)–(3.67), Figure 3.17              | §3.5.2.7    | 稳定性条件、原文矛盾标注、谱半径和阈值实验           |
| Figure 3.18                                          | §3.5.2.8    | $O(\epsilon/\alpha)$、直接解与增量更新               |
| (3.68)                                               | §3.5.2.9    | 非线性 Newton、$\alpha$-循环预条件 GMRES、时间窗影响 |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3.5.2, pp. 431–443.
