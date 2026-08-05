---
title: IMEX 与松弛型格式
description: 编号 78、91、104：把非线性项显式处理，同时保住能量论证
lang: zh
translation: en/computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation
tags:
  - 论文笔记
  - 相场模型
  - 隐显格式
---

> [!note] 本页覆盖
> 编号 **78**（_J. Comput. Phys._ 515, 2024）、**91**（_Math. Comput._ 95(359), 2026）、**104**（投稿 _SIAM J. Numer. Anal._，[arXiv:2605.05619](https://arxiv.org/abs/2605.05619)）。
>
> 三篇都已逐式核对：**91** 与 **104** 依作者自己的 arXiv 源文件核对，**78** 依已出版的 PDF 全文核对，因此下面的公式、定理假设、参数窗口与实验数值是转录而非转述。三篇都给出完整的推导链与带假设的定理。
>
> 数值实验方面：**78** 与 **91** 的实验设置与实测数据在本页照实报告；**104** 的这一版**没有偏微分方程数值实验**——摘要里那句「给出数值实验以支持理论」在源文件中被注释掉了，结论部分说明将由后续报告给出，本页照实说明，不作补写。
>
> 两处出版信息本站未能独立核实：编号 91 的卷期页码取自阅读清单，未与 AMS 页面核对；编号 78 的 DOI 字符串未在 PDF 提取文本中出现（卷号、年份与文章号 113225 已核实）。

## 贯穿本页的一个问题

隐显格式把线性刚性部分隐式处理、非线性部分显式处理，因此每级（或每步）不需要内层非线性迭代，这是它在相场计算里受欢迎的全部理由。代价集中在一处。

**全隐格式的能量论证之所以走得通，是因为非线性项与被检验的差分同号。** 用 $\nabla_\tau u^n$ 去测试 $-F'(u^n)$，得到的是 $F$ 在两层之间的差，凸性或者一条代数恒等式就把它变成一个电报式相消加一个非负余项。一旦把 $F'$ 挪到显式一侧，它求值的位置与被检验的差分错开一层，这个同号性就没有了：显式项只提供一个**符号不定**的贡献，必须靠别处的耗散把它压住。

于是问题变成：**要压住它，需要多少耗散？这个「多少」由什么对象决定？** 本页三篇给出三种回答，而它们回答的其实是同一个结构问题的三种投影。

- **编号 91** 处理隐显 Runge-Kutta。承载困难的对象是一个**微分矩阵** $D(z)$，索引跑遍级而不是时间层；它就是[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|BDF 一族]]里离散正交卷积核的级版本。论文的做法是设计 Butcher 表，使 $D$ 与网格参数**完全无关**。
- **编号 104** 处理隐显多步法。承载困难的对象是三串卷积核（$\vec a$ 差分、$\vec b$ 隐式、$\vec c$ 显式）的复合。论文把它们各自的**半生成函数**在单位圆上取极值，三个数一算，稳定性判据就出来了。
- **编号 78** 不动时间离散，改动**能量**。把非线性势用一个辅助变量二次化，格式就线性了，能量论证也就重新变成一件平凡的事——代价是耗散的对象换成了修正能量。

三条路线的共同点是：**都不试图直接证明显式项贡献符号确定，而是各自换一个可控的代数对象。** 这与[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|变步长 BDF 那一页]]的取向完全一致——先找到承载全部困难的那个对象，再为它建立独立的判据。

## 91：让平均耗散率与时空参数无关

### 直觉

隐显 Runge-Kutta 用在相场模型上有两个缺陷，一个是理论的，一个是实践的，本文一次解决两个。

**理论缺陷。** 要证明格式继承**原始**能量耗散律，必须先知道**各级解在最大范数下一致有界**——因为稳定化分裂里的参数 $\kappa$ 要压住 $\|F''\|_\infty$，而 $F''(\Phi)=3\Phi^2-1$ 只有在 $\Phi$ 有界时才有界。此前的统一框架绕开这一步的办法是直接**假设**非线性主体 $F'$ 全局 Lipschitz，而 Cahn-Hilliard 的四次势 $F(\Phi)=\frac14(\Phi^2-1)^2$ 根本不满足。这不是技术上的小疏漏：Cahn-Hilliard 方程**没有**最大值原理，所以也不能像 Allen-Cahn 那样从不变区间白拿一个 $\ell^\infty$ 界。

**实践缺陷。** 一般隐显 Runge-Kutta 的**平均耗散率**依赖 $\tau_n\overline\lambda_{\mathrm{ML}}$，即时间步乘以离散算子的平均特征值。这意味着自适应算法一旦取大步长，方法的**有效耗散就变了**，算出来的能量曲线随 $\tau_{\max}$ 漂移。对一个以自适应为卖点的方法来说，这是釜底抽薪的毛病：你无法判断能量曲线上的变化是物理的还是数值参数造成的。

**本文的关键判断是：这两件事由同一个矩阵决定。** 把格式写成差分形式后，出现一个下三角矩阵 $D(z)=D_{\mathrm E}-zD_{\mathrm{EI}}$，其中 $z=-\tau_n\overline\lambda_{\mathrm{ML}}$ 一类的量。它在级指标上扮演的角色，与 DOC 核在时间层指标上扮演的角色**完全相同**：把多级算子逆回一个逐级差分，能量论证才走得动。平均耗散率就是 $D$ 的迹的组合，因此

$$
\mathcal R\ \text{与网格无关}
\quad\Longleftrightarrow\quad
D\ \text{与}\ z\ \text{无关}
\quad\Longleftrightarrow\quad
D_{\mathrm{EI}}=\mathbf 0 .
$$

**最后一条是 Butcher 表上的一个代数方程，可以直接解。** 这就是「精化」（refined）的全部含义。而 $D_{\mathrm R}$ 一旦是纯数字矩阵，它的最小特征值也就与网格无关，正好可以用来做那个一致最大范数界的自举——**理论缺陷与实践缺陷被同一个条件消掉。**

### 问题设定

周期区域 $\Omega\subseteq\mathbb R^2$ 上，Ginzburg-Landau 自由能与其 $H^{-1}$ 梯度流：

$$
E[\Phi]=\int_\Omega\Bigl[\tfrac{\epsilon^2}{2}|\nabla\Phi|^2+F(\Phi)\Bigr]\mathrm d\mathbf x,
\qquad F(\Phi):=\tfrac14(\Phi^2-1)^2,
\qquad 0<\epsilon<1,
$$

$$
\partial_t\Phi=\Delta\bigl[F'(\Phi)-\epsilon^2\Delta\Phi\bigr],
\qquad
(\Phi(t),1)=(\Phi(t_0),1),
$$

$$
\frac{\mathrm dE}{\mathrm dt}
=\Bigl(\frac{\delta E}{\delta\Phi},\partial_t\Phi\Bigr)_{L^2}
=-\bigl((-\Delta)^{-1}\partial_t\Phi,\partial_t\Phi\bigr)_{L^2}\le0 .
$$

以稳定化参数 $\kappa\ge0$ 改写：

$$
L_\kappa\Phi:=-\epsilon^2\Delta\Phi+\kappa\Phi,
\qquad
f_\kappa(\Phi):=\kappa\Phi-F'(\Phi),
\qquad
\partial_t\Phi=\Delta\bigl[L_\kappa\Phi-f_\kappa(\Phi)\bigr].
$$

非均匀网格 $0=t_0<\cdots<t_N=T$、$\tau_n=t_n-t_{n-1}$ 上的 $s$ 级隐显 Runge-Kutta 格式为

$$
u_h^{n,i}=u_h^{n,1}
+\tau_n\sum_{j=1}^{i}a_{ij}\Delta_hL_{\kappa,h}u_h^{n,j}
-\tau_n\sum_{j=1}^{i-1}\hat a_{ij}\Delta_h f_\kappa(u_h^{n,j}),
\qquad
u_h^{n,1}:=\phi_h^{n-1},\ \ \phi_h^{n}:=u_h^{n,s},
$$

其中 $c_1=0$、$c_s=1$；隐式部分 $A$ 是**刚性精确、第一级显式**（首尾相同）的对角隐式 Runge-Kutta，$\widehat A$ 严格下三角（显式）；并要求**节点条件** $\hat{\mathbf c}=\mathbf c$（等价于 $A\mathbf 1=\widehat A\mathbf 1$），它使方法在所有级上相容，并保持平衡态 $L_{\kappa,h}\phi_h^{*}=f_\kappa(\phi_h^{*})$。空间用 **Fourier 拟谱**离散，给出 $\Delta_h$ 与 $L_{\kappa,h}$。

三阶以内的阶条件为：一阶 $\mathbf b^T\mathbf 1=\hat{\mathbf b}^T\mathbf 1=1$；二阶 $\mathbf b^T\mathbf c=\hat{\mathbf b}^T\mathbf c=\frac12$；三阶 $\mathbf b^T\mathbf c^{.2}=\hat{\mathbf b}^T\mathbf c^{.2}=\frac13$、$\mathbf b^TA\mathbf c=\hat{\mathbf b}^T\widehat A\mathbf c=\frac16$，再加两条**耦合**条件 $\mathbf b^T\widehat A\mathbf c=\hat{\mathbf b}^TA\mathbf c=\frac16$。

**步长比限制：没有。** 网格任意非均匀，$\tau_n$ 自由，**不施加任何形式的步长比条件**——这正是「稳健时间自适应」的含义。

> [!warning] 本文里出现的 $\frac{1+\sqrt2}{4}$ 不是步长比
> 全文检索可见数值 $\frac{1+\sqrt2}{4}$，但它是**对照方法** IERK(2,3)（取自 Liao-Wang-Wen）的 Butcher 系数 $a_{33}$，是一个 Runge-Kutta 系数，**不是步长比**。$1+\sqrt2$、$3.561$、$4.8645$ 这三个门槛在本文中**不出现**，本文也不需要它们。把这个 $\frac{1+\sqrt2}{4}$ 与[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 48]] 的条件 S0 联系起来是误读。

### 推导

**第一步：差分形式与微分矩阵。** 记 $s_{\mathrm I}:=s-1$、$\delta_\tau u^{n,\ell+1}:=u^{n,\ell+1}-u^{n,\ell}$、$u^{n,\ell+\frac12}:=(u^{n,\ell+1}+u^{n,\ell})/2$，并令

$$
A_{\mathrm I}:=(a_{i+1,j+1})_{i,j=1}^{s_{\mathrm I}},
\qquad
A_{\mathrm E}:=(\hat a_{i+1,j})_{i,j=1}^{s_{\mathrm I}},
\qquad
E_{s_{\mathrm I}}:=(1_{i\ge j})\ \text{（下三角全一矩阵）},
$$

则格式等价于

$$
\sum_{\ell=1}^{i}d_{i\ell}\bigl(\tau_n\Delta_hL_{\kappa,h}\bigr)\,\delta_\tau u_h^{n,\ell+1}
=\tau_n\Delta_h\Bigl[L_{\kappa,h}u_h^{n,i+\frac12}-f_\kappa(u_h^{n,i})\Bigr],
\qquad 1\le i\le s_{\mathrm I},
$$

$$
D(z):=D_{\mathrm E}-zD_{\mathrm{EI}},
\qquad
D_{\mathrm E}:=A_{\mathrm E}^{-1}E_{s_{\mathrm I}},
\qquad
D_{\mathrm{EI}}:=A_{\mathrm E}^{-1}A_{\mathrm I}E_{s_{\mathrm I}}-E_{s_{\mathrm I}}+\tfrac12 I_{s_{\mathrm I}} .
$$

下三角矩阵 $D$ 称为正（半）定，指其对称部分 $\mathcal S(D)=(D+D^T)/2$ 正（半）定。**这个 $D$ 承担的正是 DOC 核在 BDF 各篇中承担的角色**，只是索引从时间层换成级；相应的正交性恒等式为

$$
\sum_{i=j}^{k}d^{(R)}_{k,i}\,\underline{\hat a}_{i+1,j}\equiv\delta_{kj},
\qquad\text{于是}\qquad
\sum_{i=1}^{k}d^{(R)}_{k,i}\sum_{j=1}^{i}\underline{\hat a}_{i+1,j}v^j\equiv v^k,
\qquad
(\underline{\hat a}_{i+1,j}):=E_{s_{\mathrm I}}^{-1}A_{\mathrm E}.
$$

**第二步：平均耗散率。** 在 $D_{\mathrm E}$ 与 $D_{\mathrm{EI}}$ 半正定的前提下，各级能量律带有

$$
\mathcal R=\frac1{s_{\mathrm I}}\mathrm{tr}(D_{\mathrm E})
+\frac1{s_{\mathrm I}}\mathrm{tr}(D_{\mathrm{EI}})\,\tau_n\overline\lambda_{\mathrm{ML}}
=\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}\frac1{\hat a_{k+1,k}}
+\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}
\Bigl(\frac{a_{k+1,k+1}}{\hat a_{k+1,k}}-\frac12\Bigr)\tau_n\overline\lambda_{\mathrm{ML}}
\ \ge 0,
$$

其中 $\overline\lambda_{\mathrm{ML}}>0$ 是对称正定矩阵 $-\Delta_hL_{\kappa,h}$ 的平均特征值。一个方法「好」的标准是：$\mathcal R$ 在很大的 $\tau_n\overline\lambda_{\mathrm{ML}}$ 范围内都尽可能接近 $1$。

**第三步：精化条件。** $\mathcal R$ 与 $\tau_n\overline\lambda_{\mathrm{ML}}$ 无关**当且仅当** $D_{\mathrm{EI}}=\mathbf 0$，即

$$
A_{\mathrm E}^{-1}A_{\mathrm I}E_{s_{\mathrm I}}-E_{s_{\mathrm I}}+\tfrac12I_{s_{\mathrm I}}=\mathbf 0
\qquad\Longleftrightarrow\qquad
A_{\mathrm I}=A_{\mathrm E}P_{s_{\mathrm I}},
\qquad
P_{s_{\mathrm I}}:=I_{s_{\mathrm I}}-\tfrac12E_{s_{\mathrm I}}^{-1} .
$$

此时 $D_{\mathrm R}:=(d^{(R)}_{ij})=D_{\mathrm E}=A_{\mathrm E}^{-1}E_{s_{\mathrm I}}$ **与 $z$ 无关，因而与两个网格参数都无关**，且

$$
\mathcal R_{\mathrm R}=\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}\frac1{\hat a_{k+1,k}} .
$$

格式塌缩为紧凑的精化形式

$$
u_h^{n,i+1}=u_h^{n,1}+\tau_n\sum_{j=1}^{i}\hat a_{i+1,j}\Delta_h
\Bigl[L_{\kappa,h}u_h^{n,j+\frac12}-f_\kappa(u_h^{n,j})\Bigr],
$$

等价地 $\delta_\tau u_h^{n,i+1}=\tau_n\sum_{j=1}^{i}\underline{\hat a}_{i+1,j}\Delta_h[L_{\kappa,h}u_h^{n,j+\frac12}-f_\kappa(u_h^{n,j})]$。

**第四步：一条结构性的排除。** 节点条件迫使隐式部分的第一列

$$
\mathbf a_1=\bigl(\tfrac12\hat a_{21},\dots,\tfrac12\hat a_{s1}\bigr)^T\ne\mathbf 0,
$$

因此这类方法必然是 **Lobatto 型，绝不可能是 Radau 型或 ARS 型**。论文把这一点写成命题：不存在平均耗散率与 $\tau_n\overline\lambda_{\mathrm{ML}}$ 无关的 Radau 型或 ARS 型隐显 Runge-Kutta 方法。这类结论在方法设计中很有价值——它不是「我们选了 Lobatto 型」，而是「要这条性质就只能是 Lobatto 型」。相应地，隐式部分的级阶为二，且**不必代数稳定**。

**第五步：三族具体方法。**

| 方法                      | 构造                                                                                                                                                                                                                         | 参数窗口与耗散率                                                                                                 |
| ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| R-IERK(1,2)               | 在两级族 IERK(1,2;$\theta$) 中取 $\theta=\frac12$（该族 $D^{(1,2)}(z)=1-z(\theta-\frac12)$，$z\le0$ 时正定当且仅当 $\theta\ge\frac12$，$\mathcal R^{(1,2)}(\theta)=1+(\theta-\frac12)\tau_n\overline\lambda_{\mathrm{ML}}$） | $\mathcal R^{(1,2)}_{\mathrm R}=1$ 精确成立                                                                      |
| R-IERK(2,4;$c_2$)         | 二阶无三级方法；四级单参数族取 $\hat a_{32}=\hat a_{43}=c_2$、$\hat a_{42}=\frac1{2c_2}-c_3$                                                                                                                                 | $\frac{\sqrt3-1}{2}\le c_2<1$ 或 $1<c_2\le\frac{2+\sqrt6}{2}$，$\mathcal R^{(2,4)}_{\mathrm R}(c_2)=\frac1{c_2}$ |
| R-IERK(3,6;$\hat a_{52}$) | 三阶无四级方法，五级也未找到正（半）定者；六级族固定 $c_2=1$、$c_3=\frac45$、$c_4=\frac7{10}$、$c_5=\frac{12}{25}$、$\hat a_{32}=\frac65$，$\hat a_{51}=\frac{17\hat a_{52}}{103}-\frac{86756827361}{181963162950}$          | $0.664767<\hat a_{52}<0.751947$                                                                                  |

R-IERK(1,2) 就是 Crank-Nicolson 型的

$$
\delta_\tau\phi_h^{n}=\tau_n\Delta_h\Bigl[\tfrac12L_{\kappa,h}(\phi_h^{n}+\phi_h^{n-1})
-f_\kappa(\phi_h^{n-1})\Bigr].
$$

四级族中 $c_3$ 是

$$
c_3^2-\Bigl(\frac1{2c_2}+c_2\Bigr)c_3+(c_2-1)^2=0
$$

的根，显式地 $c_3^{*}=\frac32$（当 $c_2=1$），否则

$$
c_3^{*}=\frac1{4c_2}+\frac{c_2}{2}
-\frac{\sqrt{(-2c_2^2+4c_2+1)(6c_2^2-4c_2+1)}}{4c_2}>0,
\qquad 0<c_2\le\tfrac{2+\sqrt6}{2}.
$$

正定性由两个显式行列式判定：

$$
\mathrm{Det}\,\mathcal S\bigl(D^{(2,4)}_{\mathrm R,2};c_2\bigr)=\frac{c_3(4c_2-c_3)}{4c_2^4}>0,
\qquad
\mathrm{Det}\,\mathcal S\bigl(D^{(2,4)}_{\mathrm R};c_2\bigr)
=\frac{1}{8c_2^7}(c_2-1)^2(2c_2^2+2c_2-1)>0 .
$$

六级族的耗散率为

$$
\mathcal R^{(3,6)}_{\mathrm R}(\hat a_{52})
=\frac{36392632590}{123664285500\,\hat a_{52}+89552314349}
+\frac{219055887768899}{156795586342500}
\approx\frac{1.39708\,\hat a_{52}+1.30599}{\hat a_{52}+0.724157}.
$$

> [!note] 论文自己划定的两条边界
> 第一，$c_2=1$ 给出最优耗散率 $\mathcal R=1$，但此时 $\mathcal S(D_{\mathrm R}^{(2,4)};1)$ **奇异**（$\lambda_{\min}=0$），因此被理论**排除在外**（注 3.1）——尽管它在数值上表现最好。这是一处真实的理论-实践落差，论文没有回避。第二，六级方法能否达到 $\mathcal R^{(3,6)}_{\mathrm R}=1$ 仍是公开问题；四阶需要至少九级（28 条阶条件），论文未作尝试。

**第六步：怎样替代最大值原理。** Cahn-Hilliard 方程**没有**最大值原理，所以**本文没有最大值原理定理**。取而代之的是各级解的一致最大范数界，由一次**更新过的时空误差分裂**得到：

1. 先用能量论证配合各级亏量的粗略设定，得到**时间半离散**各级解 $U^{n,i}$ 的 $H^{m+4}$ 正则性；
2. 把全离散各级误差分裂为 $U^{n,\ell}-u_h^{n,\ell}=(U^{n,\ell}-U_M^{n,\ell})+e_h^{n,\ell}$，其中 $U_M=P_MU$ 是 $L^2$ 投影；
3. 自举出一个粗略的最大范数界 $\hat c_4^{*}/\epsilon^2$。

两条引理让这套论证成立。**引理 4.3（DOC 矩阵不等式）**：若 $D_{\mathrm R}$ 与 $D_{\mathrm R}^{-1}$ 都正定，记 $\lambda_{\min}$、$\sigma_{\min}$ 分别为 $\mathcal S(D_{\mathrm R})$ 与 $\mathcal S(D_{\mathrm R}^{-1})$ 的最小特征值，则对任意序列 $\{v^j\},\{u^j\}$，

$$
\sum_{i=1}^{k}\sum_{j=1}^{i}d^{(R)}_{i,j}v^jv^i\ge\lambda_{\min}\sum_i(v^i)^2,
\qquad
\sum_{i}\sum_{j\le i}\underline{\hat a}_{i+1,j}v^jv^i\ge\sigma_{\min}\sum_i(v^i)^2,
$$

$$
\sum_i\sum_{j\le i}d^{(R)}_{i,j}v^ju^i\le\frac1{\sigma_{\min}}\sum_i|v^i||u^i|,
\qquad
\sum_i\sum_{j\le i}\underline{\hat a}_{i+1,j}v^ju^i\le\frac1{\lambda_{\min}}\sum_i|v^i||u^i| .
$$

**这里的 $\lambda_{\min},\sigma_{\min}$ 是纯数字，因为 $D_{\mathrm R}$ 与网格无关——这正是精化条件在收敛性证明里的第二次回报。**

**引理 4.4（非线性项的局部处理）**：用

$$
F'(v)-F'(w)=(v-w)\int_0^1F''\bigl[\gamma v+(1-\gamma)w\bigr]\,\mathrm d\gamma
$$

以及 $\delta_\tau F'$ 的相应恒等式，把 Lipschitz 常数局部化到一个球上。**这一步就是全局 Lipschitz 假设被去掉的地方。**

各级能量不等式本身为（设各级解最大范数被 $c_0$ 界住，且 $\kappa\ge\max_{\|\xi\|_\infty\le c_0}\|F''(\xi)\|_\infty$）

$$
\bigl(L_{\kappa,h}u^{n,i+\frac12}-f_\kappa(u^{n,i}),\delta_\tau u^{n,i+1}\bigr)
\le E[u^{n,i}]-E[u^{n,i+1}]
-\tfrac12\|\delta_\tau u^{n,i+1}\|^2
\bigl(\kappa-\max_{\xi_h\in\mathcal B_{n,i}}\|F''(\xi)\|_\infty\bigr).
$$

另外，**逐级体积守恒是精确的**：$(u^{n,i+1},1)=(u^{n,1},1)=(\phi^0,1)$ 对一切 $n,i$ 成立。

### 定理

**（标准正则性假设）** 存在整数 $m\ge1$、$p\ge1$ 与 $c_\phi>0$ 使

$$
\|\Phi^0\|_{H^{m+4}}+\sum_{k=0}^{2}\|\partial_t^{(k)}\Phi(t)\|_{H^{m+4-k}}
+\sum_{k=3}^{p+1}\|\partial_t^{(k)}\Phi(t)\|_{L^2}\le c_\phi,
\qquad 0<t<T .
$$

**（引理 2.1，一般隐显 Runge-Kutta 的各级能量律）** 若 $D_{\mathrm E},D_{\mathrm{EI}}$ 正（半）定，各级解最大范数被 $c_0$ 界住，且 $\kappa\ge\max_{\|\xi\|_\infty\le c_0}\|F''(\xi)\|_\infty$，则对 $n\ge1$、$1\le i\le s_{\mathrm I}$，

$$
E[u^{n,i+1}]-E[u^{n,1}]\le\frac1{\tau_n}\sum_{k=1}^{i}
\Bigl(\Delta_h^{-1}\delta_\tau u^{n,k+1},\ \sum_{\ell=1}^{k}
d_{k\ell}(\tau_n\Delta_hL_{\kappa,h})\,\delta_\tau u^{n,\ell+1}\Bigr),
$$

右端即上面那个非负的平均耗散率 $\mathcal R$ 所度量的量。

**（定理 4.1，时间半离散各级解的正则性）** 在上述正则性假设（$m\ge1$）下，并设 **$D_{\mathrm R}=A_{\mathrm E}^{-1}E_{s_{\mathrm I}}$ 与 $D_{\mathrm R}^{-1}$ 都正定**，若最大时间步 $\tau$ 充分小，则存在 $C_\phi>0$、**与步长 $\tau_n$ 无关**，使

$$
\|U^{n,i}\|_{H^{m+4}}+\bigl\|(U^{n,i}-U^{n,1})/\tau_n\bigr\|_{H^{m}}
\le C_\phi/\epsilon^{2},
\qquad 1\le n\le N,\ 2\le i\le s .
$$

**（定理 4.2，各级都满足的原始能量耗散律）** 在同样假设下，再设 $h$ 与 $\tau$ 充分小、$\kappa\ge\max_{\|\xi\|_\infty\le\hat c_4^{*}/\epsilon^2}\|F''(\xi)\|_\infty$，则各级解 $u_h^{n,i}$ 最大范数有界，且

$$
E[u^{n,j+1}]-E[u^{n,1}]\le\frac1{\tau}\sum_{k=1}^{j}
\Bigl(\Delta_h^{-1}\delta_\tau u^{n,k+1},\ \sum_{\ell=1}^{k}
d^{(R)}_{k\ell}\,\delta_\tau u^{n,\ell+1}\Bigr),
\qquad 1\le n\le N,\ 1\le j\le s_{\mathrm I} .
$$

**衰减的是原始能量 $E$，不是修正能量，也不含辅助变量；而且没有任何步长比限制。** 这里的 $d^{(R)}_{k\ell}$ 是纯数字，与网格无关。

**（定理 4.3，全精度 $L^2$ 收敛）** 同样假设下，再设各级亏量满足 $\zeta^{n,i+1}_{\mathrm R}=0$（$1\le i\le s_{\mathrm I}$）与 $\|\zeta^{n,s}_{\mathrm R}\|\le c_2\tau^{p}$，则 $u_h^n$ 在 $L^2$ 范数下以 $\mathcal O(\tau^{p}+h^{m})$ 收敛。这里的**无条件**指：$\tau$ 与 $h$ 之间不需要任何耦合，也不需要步长比限制。

论文声明这是**首次**在**不假设非线性主体全局 Lipschitz** 的情况下，为 Cahn-Hilliard 模型的隐显 Runge-Kutta 方法建立原始能量耗散律与无条件 $L^2$ 收敛。

### 数值实验

空间一律用 Fourier 拟谱；自适应步长规则为

$$
\tau_{\mathrm{ada}}=\max\{\tau_{\min},\ \tau_{\max}/\Pi_\eta(\phi)\},
\qquad
\Pi_\eta(\phi):=\sqrt{1+\eta\|\partial_\tau\phi^n\|^2},
$$

默认 $\eta=1000$、$\tau_{\min}=10^{-4}$、$\tau_1=\tau_{\min}$。

**算例 5.1（精度检验）。** 制造解 $\Phi(x,y;t)=e^{-t}\sin(\pi x)\sin(\pi y)$，$\Omega=(0,2)^2$，$\epsilon=0.2$，$64\times64$ 网格，$T=1$，$\kappa=4$，步长 $\tau=2^{-k}/10$（$0\le k\le9$），误差 $e(\tau)=\max_n\|\Phi_h^n-\Phi(t_n)\|_\infty$。结果：R-IERK(2,4;$c_2$) 二阶，R-IERK(3,6;$\hat a_{52}$) 三阶。**参数的影响很不对称**：不同的 $c_2$ 给出肉眼可辨的精度差异（$\tau<10^{-2}$ 时 $c_2=1$ 最准），而不同的 $\hat a_{52}$ 给出几乎相同的解。

**算例 2.1（引出自适应问题的算例）。** $\Omega=(-\pi,\pi)^2$，$\epsilon=0.1$，$\kappa=2.5$，$T=1000$，初值

$$
\Phi^0=\tfrac12\tanh(|x|+|y|+1)-e^{-5(|x|+|y|-2)^2}
+\tfrac12e^{-2(|x|-1)^2}+\tfrac1{10}\sin\bigl(e^{|y|-1}\bigr),
$$

参考解由 Lobatto 型 IERK(2,3)（$a_{33}=\frac{1+\sqrt2}{4}$）在 $\tau=10^{-4}$ 下给出。取 $\tau_{\max}\in\{0.01,0.05,0.1\}$：**R-IERK(1,2)（$\mathcal R=1$）的能量曲线始终贴住参考解，而 IERK(1,2;$\theta$)（$\theta=1,\frac32$）的曲线随 $\tau_{\max}$ 增大而漂移。** 这就是引言里那个实践缺陷的直接演示。

**算例（R-IERK(2,4;$c_2$) 的能量与效率）。** 取 $c_2=\frac12,1,\frac32$，$\tau_{\max}\in\{0.1,0.2,0.5\}$；所有离散能量都单调下降；IERK(2,3) 的曲线随 $\tau_{\max}$ 显著变化，R-IERK 的曲线稳健。效率表（R-IERK(2,4;1)，$T=450$）：

| 步长设定            | CPU 时间    | 时间层数          |
| ------------------- | ----------- | ----------------- |
| 固定 $\tau=10^{-4}$ | $6724.63$ s | $4.5\times10^{6}$ |
| $\tau_{\max}=0.1$   | $11.90$ s   | $11714$           |
| $\tau_{\max}=0.2$   | $6.42$ s    | $5840$            |
| $\tau_{\max}=0.5$   | $2.62$ s    | $2304$            |

**这张表是本页最直接的一处「理论换来性能」的记录**：约三个数量级的加速，而这只有在能量曲线不随 $\tau_{\max}$ 漂移的前提下才敢用——也就是说，$\mathcal R_{\mathrm R}$ 的网格无关性不是一条美学性质，它换来的是能不能真的把步长放大。

**算例（R-IERK(3,6;$\hat a_{52}$) 的能量）。** 取 $\hat a_{52}=\frac23,\frac7{10},\frac34$，自适应参数 $\eta=500$，$\tau_{\max}\in\{0.2,0.5,0.8\}$；参考方法为 Lobatto 型 IERK(3,5)（$a_{43}=-\frac35$，$\mathcal R^{(3,5)}_{\mathrm L}=\frac54+\frac25\tau\overline\lambda_{\mathrm{ML}}$）。三个 $\tau_{\max}$ 下 R-IERK 的能量曲线**互相无法分辨**，而 IERK(3,5) 的曲线显著变化。注意 $\mathcal R^{(3,5)}_{\mathrm L}$ 显式含 $\tau\overline\lambda_{\mathrm{ML}}$，这就是它漂移的原因，对照非常干净。

**算例（粗化动力学）。** 用 R-IERK(3,6;$\frac34$)、$\eta=500$、$\tau_{\max}=0.5$，给出 $t=0,2,20,200,450,1000$ 处的解形貌。

**这批实验建立了什么、没建立什么。** 建立的是：精化方法的能量曲线确实与 $\tau_{\max}$ 无关，因而自适应可以放心取大步长，效率提升约三个数量级；观测阶与定理一致。没有建立的是理论边界的尖锐性——$c_2=1$ 在理论上被排除（$\mathcal S(D_{\mathrm R})$ 奇异）却在实验里最准，说明正定性要求对二阶方法而言可能是充分而非必要的。论文如实记录了这个矛盾，未加掩饰。

### 与其他论文的关系

微分矩阵 $D$ 是 DOC 核的**级指标化身**：编号 [[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|48、52、58、67]] 在时间层指标上用它，本文在级指标上用它；论文明确把 DOC 技术归于 Liao-Ji-Wang-Zhang（2022）、Liao-Ji-Zhang（2022，相场晶体）与 Liao-Zhang（2021）。

**方法论上最锋利的对照在这里**：编号 48、52 为变步长 BDF2 争取稳健性的办法是**限制步长比**（$r_k<3.561$），而本文为一个**多级**方法争取稳健性的办法是**消掉耗散率对网格的依赖**，因此**完全不需要步长比条件**。两条路线针对的是同一个毛病的两种表现。

与编号 78 的分歧在能量的类型：本文衰减的是**原始**能量 $E[\Phi]$，不是含辅助变量的修正能量——这与[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 40、57]] 及编号 48、52 的立场一致，而本文引言正是批评基于 SAV 的高阶格式只能建立关于修正能量的稳定性。编号 104 处理的是同一问题的**多步**版本，工具从微分矩阵换成生成函数。

最后，本文解决的那个「各级解一致有界」的困难，正是此前 Fu-Tang-Yang 与 Shin-Lee-Lee 的 Radau 型／ARS 型构造被卡住的地方——两者都需要全局 Lipschitz 连续性。

## 104：用半生成函数取代逐族寻找乘子

### 直觉

对非线性抛物方程，隐显多步方法把线性刚性算子隐式、非线性显式处理。已有的族很多——加权 BDF（Li-Xie）、修正 BDF（Akrivis-Karakatsani）、推广 BDF（Huang-Shen）、NIMEX（Rosales-Seibold-Shirokoff-Zhou）——都以扩大绝对稳定域为设计目标。但三阶及更高阶的变体**不 $A$-稳定**，因此严格的**离散能量**稳定性此前是长期开放的问题。

两条现有路线都要求**逐族**寻找特设对象，这是本文真正针对的痛点。第一条是基于 Dahlquist $G$-稳定的 Nevanlinna-Odeh 型乘子：加权 BDF、修正 BDF、推广 BDF、NIMEX 各需要一族不同的乘子，而乘子难找。第二条是 Huang-Shen 型的隐式部分分解

$$
\sum_{j=0}^{k-1}b_{\mathrm G,j}^{(k,\beta)}v^{n-j}
=\eta_k(\beta)\sum_{j=0}^{k-1}c_{\mathrm G,j}^{(k,\beta)}v^{n-j}
+\sum_{j=0}^{k-1}d_j^{(k,\beta)}v^{n-j},
\qquad
\eta_2=\tfrac{\beta-1}{\beta},\
\eta_3=\tfrac{\beta-1}{\beta+1},\
\eta_4=\tfrac{\beta-1}{\beta+3},
$$

稳定性条件为 $\eta_k(\beta)>\mu_0/\varpi$。这条路线的局限很具体：它在 $\beta=1$ 处退化（因此**从不覆盖普通 BDF-$k$**），只对 $2\le k\le4$ 存在，其三项精化版本（Huang-Shen 2025）只在固定的 $\beta_k=3,6,9$ 处有效；对推广 BDF5 以及加权 BDF、修正 BDF、NIMEX 都没有对应物。

**本文把上述全部替换为一次计算：三个显式有理函数在单位圆上的三个极值。**

支撑这一步的观察出奇简单。多步格式在被 DOC 核作用之后，隐式与显式部分各变成一个**下三角 Toeplitz** 矩阵，而下三角 Toeplitz 矩阵的**对称部分**恰好有经典的 Grenander-Szegő 生成函数——并且这个生成函数就等于**单边**和 $a(\theta)=\sum_{k\ge0}a_ke^{\imath k\theta}$ 的**实部**。也就是说，虽然核是因果的（单边的），Toeplitz-Carathéodory 那套判据仍然可用，只要换成取实部。这就是「半生成函数」这个名字的含义。

有了它，三串核的所有问题都变成三个复有理多项式在单位圆上的极值问题，而**下三角 Toeplitz 矩阵的乘积仍是下三角 Toeplitz 且可交换**，所以复合核的生成函数就是各自生成函数的乘积与商。**逐族寻找乘子的工作被一次性取消了。**

### 问题设定

在 Hilbert 三元组 $V\subset H=H'\subset V'$ 上考虑

$$
u_t+\varpi\mathcal Lu=\mathcal F(u),\qquad 0<t<T,\qquad u(0)=u^0\in H,\ \varpi>0,
$$

$\mathcal L:V\to V'$ 正定、自伴、有界、线性，$\mathcal F$ 可非线性。范数取 $\|v\|_H=\langle v,v\rangle^{1/2}$、$\|v\|_V=\langle\mathcal Lv,v\rangle^{1/2}$，对偶范数 $\|v\|_{\star}=\sup_{\|w\|_V=1}|\langle v,w\rangle|$。标准假设是球 $\mathcal B_{u(t)}=\{v\in V:\|v-u(t)\|_V\le1\}$ 上的**局部** Lipschitz 条件

$$
\|\mathcal F(v)-\mathcal F(w)\|_{\star}\le\mu_0\|v-w\|_V+\mu_1\|v-w\|_H,
\qquad \mu_0\in(0,\varpi),\ \mu_1\ \text{任意} .
$$

该框架也覆盖非自伴的 $\mathcal L=\mathcal L_s+\mathcal L_a$：把低阶的反自伴部分 $\mathcal L_a$ 移入显式项即可。

**均匀网格** $\tau=t_j-t_{j-1}$ 上的一般 $k$ 步隐显多步方法为

$$
\sum_{j=0}^{k-1}a_j^{(k)}\partial_\tau u^{n-j}
+\varpi\sum_{j=0}^{k}b_j^{(k)}\mathcal Lu^{n-j}
=\sum_{j=0}^{k-1}c_j^{(k)}\mathcal F(u^{n-j-1})+\mathfrak C^{(k)}_n(u^0),
\qquad n\ge1,
$$

其中 $a_0^{(k)},b_0^{(k)},c_0^{(k)}>0$，启动修正项 $\mathfrak C^{(k)}_n(u^0)$（$n\ge k$ 时为零）假定可得，使格式从第一步起就 $k$ 阶相容。**因此一个方法就是一个三元组 $(\vec a^{(k)},\vec b^{(k)},\vec c^{(k)})$。** 把核按零延拓（$j\ge k+1$ 时 $b_j^{(k)}=0$，$j\ge k$ 时 $a_j^{(k)}=c_j^{(k)}=0$），得到等价的全局卷积形式

$$
\sum_{j=1}^{n}a^{(k)}_{n-j}\partial_\tau u^{j}
+\varpi\sum_{j=1}^{n}b^{(k)}_{n-j}\mathcal Lu^{j}
=\sum_{j=1}^{n}c^{(k)}_{n-j}\mathcal F(u^{j-1})+\mathfrak C^{(k)}_n(u^0).
$$

**步长比限制：没有——网格是均匀的。最大值原理论证：没有。**

### 推导

**第一步：DOC 核与全局离散能量方法。** $\vec a^{(k)}$ 的 DOC 核为

$$
a_0^{(-1,k)}:=\frac{1}{a_0^{(k)}},
\qquad
a_j^{(-1,k)}:=-\frac{1}{a_0^{(k)}}\sum_{i=1}^{j}a^{(-1,k)}_{j-i}a_i^{(k)}\ (j\ge1),
$$

满足双向正交性

$$
\sum_{\ell=j}^{n}a^{(-1,k)}_{n-\ell}a^{(k)}_{\ell-j}\equiv\delta_{nj}
=\sum_{\ell=j}^{n}a^{(k)}_{n-\ell}a^{(-1,k)}_{\ell-j},
\qquad 1\le j\le n .
$$

作用之后格式化为「差分」形式，这是能量论证的起点：

$$
\partial_\tau u^{n}+\varpi\sum_{\ell=1}^{n}\hat b^{(k)}_{n-\ell}\mathcal Lu^{\ell}
=\sum_{\ell=1}^{n}\hat c^{(k)}_{n-\ell}\mathcal F(u^{\ell-1})
+\sum_{\ell=1}^{n}a^{(-1,k)}_{n-\ell}\mathfrak C^{(k)}_\ell(u^0),
$$

**复合核**为

$$
\hat b_j^{(k)}:=\sum_{i=0}^{j}a^{(-1,k)}_{j-i}b_i^{(k)},
\qquad
\hat c_j^{(k)}:=\sum_{i=0}^{j}a^{(-1,k)}_{j-i}c_i^{(k)} .
$$

矩阵形式下，记 $A_{L,k},B_{L,k},C_{L,k}$ 为三串系数对应的下三角 Toeplitz 矩阵，则 $A^{(-1)}_{L,k}=A_{L,k}^{-1}$，且

$$
\widehat B_{L,k}=A_{L,k}^{-1}B_{L,k},
\qquad
\widehat C_{L,k}=A_{L,k}^{-1}C_{L,k},
$$

**下三角 Toeplitz 矩阵之积仍是下三角 Toeplitz 且可交换**——这一条平凡的事实是后面一切化归的技术基础。由于整个过程没有丢弃任何信息，作者称之为**全局离散能量方法**。

**第二步：半生成函数（本文的装置，引理 2.1）。** 对实序列 $\{a_0,a_1,\dots\}$（$j<0$ 时 $a_j=0$）定义

$$
a(\theta):=\sum_{k=0}^{\infty}a_ke^{\imath k\theta}\in L^2([0,2\pi)) .
$$

对二次型 $Q_n:=\sum_{k=1}^{n}w_k\sum_{j=1}^{k}a_{k-j}w_j$（对应下三角 Toeplitz 矩阵 $P_{L,n}$）：

1. $Q_n$ 正定**当且仅当** $\Re[a(\theta)]>0$ 于 $[0,2\pi)$；
2. $\min_\theta\Re[a(\theta)]\le\lambda_j(Q_n)\le\max_\theta\Re[a(\theta)]$；
3. 特征值与 $\Re[a(2\pi j/n)]$ **等分布**：$\lim_{n\to\infty}\frac1n\sum_{j=0}^{n-1}\bigl(\lambda_j(Q_n)-\Re[a(2\pi j/n)]\bigr)=0$。

证明的关键观察只有一句：**对称部分** $\mathcal S(P_{L,n})$ 的经典 Grenander-Szegő 生成函数是

$$
\mathrm g(\theta)=a_0+\sum_{k\ge1}a_k\cos k\theta=\Re[a(\theta)],
$$

也就是说**单边（因果）和的实部恰好等于经典的双边生成函数**。这正是把 Toeplitz-Carathéodory 判据推广到单边序列的那一步。

**第三步：复合规则（引理 2.2）。** (i) 若 $\hat b_j=\sum_{k=0}^{j}a_{j-k}b_k$ 则 $\hat b(\theta)=a(\theta)b(\theta)$；(ii) 若 $\{\xi_j\}$ 是 $\{a_j\}$ 的 DOC 核则 $\xi(\theta)=1/a(\theta)$。于是，记

$$
a^{(k)}(\theta)=\sum_{j=0}^{k-1}a_j^{(k)}e^{\imath j\theta},
\qquad
b^{(k)}(\theta)=\sum_{j=0}^{k}b_j^{(k)}e^{\imath j\theta},
\qquad
c^{(k)}(\theta)=\sum_{j=0}^{k-1}c_j^{(k)}e^{\imath j\theta},
$$

就有 $a^{(-1,k)}(\theta)=1/a^{(k)}(\theta)$、$\hat b^{(k)}(\theta)=b^{(k)}(\theta)/a^{(k)}(\theta)$、$\hat c^{(k)}(\theta)=c^{(k)}(\theta)/a^{(k)}(\theta)$。**一切都归结为单位圆上的三个复有理多项式。**

**第四步：三个极值常数（引理 2.4）。**

$$
\sigma_{\mathrm F}^{(k)}=\max_{\theta\in[0,2\pi)}\Bigl|\frac{1}{a^{(k)}(\theta)}\Bigr|,
\qquad
\sigma_{\mathrm E}^{(k)}=\max_{\theta\in[0,2\pi)}\Bigl|\frac{c^{(k)}(\theta)}{a^{(k)}(\theta)}\Bigr|,
\qquad
\lambda_{\mathrm I}^{(k)}=\min_{\theta\in[0,2\pi)}\Re\Bigl[\frac{b^{(k)}(\theta)}{a^{(k)}(\theta)}\Bigr],
$$

分别称为**扰动放大因子**、**非线性放大因子**与**耗散保持因子**。三者的作用是

$$
\|A_{L,k}^{-1}\|_{\ell^2}\le\sigma_{\mathrm F}^{(k)},
\qquad
\|\widehat C_{L,k}\|_{\ell^2}\le\sigma_{\mathrm E}^{(k)},
\qquad
\lambda\bigl(\mathcal S(\widehat B_{L,k})\bigr)>\lambda_{\mathrm I}^{(k)},
$$

对任意序列 $\{v^i\},\{u^i\}$ 展开则是

$$
\sum_{i=1}^{n}\sum_{j=1}^{i}\hat b^{(k)}_{i-j}v^jv^i\ge\lambda_{\mathrm I}^{(k)}\sum_{i=1}^{n}|v^i|^2,
$$

$$
\sum_{i}\sum_{j\le i}a^{(-1,k)}_{i-j}v^ju^i
\le\sigma_{\mathrm F}^{(k)}\Bigl(\sum_i|v^i|^2\Bigr)^{1/2}\Bigl(\sum_i|u^i|^2\Bigr)^{1/2},
$$

以及把 $a^{(-1,k)}$ 换成 $\hat c^{(k)}$、$\sigma_{\mathrm F}^{(k)}$ 换成 $\sigma_{\mathrm E}^{(k)}$ 的第三条。证明的成分是：第一条用 Grenander-Szegő 加 Cauchy 交错定理；谱范数界 $\|P_L\|_{\ell^2}\le\max_\theta|a(\theta)|$ 用 $\ell^2(\mathbb N)$ 上的 Parseval 等式，取 $Y(\theta)=X(\theta)\overline{a(\theta)}$。

**第五步：尖锐的归一化（引理 2.5 与推论 3.4）。** 相容性意味着 $a^{(k)}(0)=b^{(k)}(0)=c^{(k)}(0)=1$，于是

$$
\sigma_{\mathrm F}^{(k)}\ge1,
\qquad
\sigma_{\mathrm E}^{(k)}\ge1,
\qquad
\lambda_{\mathrm I}^{(k)}\le1,
$$

且三者同时取等号的**唯一**情形是隐显 Euler 格式（$k=1$）。**这条归一化让三个因子成为可比较的量**：它们都以隐显 Euler 为基准，偏离得越远，方法就越「难控」。

### 定理

**（定理 3.2，统一的无条件稳定性与收敛性）** 在 $\mathcal F$ 的局部 Lipschitz 条件、$u$ 的充分正则性、引理 2.4 的假设与 $k$ 阶相容条件之下，**若**

$$
\frac{\lambda_{\mathrm I}^{(k)}}{\sigma_{\mathrm E}^{(k)}}>\frac{\mu_0}{\varpi}
$$

且步长 $\tau$（依赖于 $\lambda_{\mathrm I}^{(k)}/\sigma_{\mathrm F}^{(k)}$）充分小，则该 $k$ 步隐显多步方法稳定且以 $\mathcal O(\tau^{k})$ 阶收敛。证明是对误差 $\tilde u^j=U^j-u^j$ 的界 $\|\tilde u^{\ell}\|_V\le1$ 作完全数学归纳，截断误差满足 $\|R_n^{(k)}\|_{\star}\le c_u\tau^{k}$。

**这个判据的形状值得停一下**：左端 $\lambda_{\mathrm I}^{(k)}/\sigma_{\mathrm E}^{(k)}$ 只依赖**方法**，右端 $\mu_0/\varpi$ 只依赖**模型**。判据把两者干净地分开了，而这正是它能一次覆盖所有族的原因。

**（隐显可控强度）** 由此定义

$$
\mathfrak I_{\mathrm{IE}}^{(k)}:=\frac{\lambda_{\mathrm I}^{(k)}}{\sigma_{\mathrm E}^{(k)}}
=\frac{\min_{\theta\in[0,2\pi)}\Re\bigl[b^{(k)}(\theta)/a^{(k)}(\theta)\bigr]}
{\max_{\theta\in[0,2\pi)}\bigl|c^{(k)}(\theta)/a^{(k)}(\theta)\bigr|}\ \le\ 1,
$$

最优值 $1$ 由隐显 Euler 格式达到。**这是本文的实用产出：每个方法一个数，据此可以给格式排序、给参数定值。**

**（新的 $\gamma$ 参数化 SIEMS 方法）** 由三个特征多项式定义：

$$
\tilde\varrho^{(k)}_{a,\mathrm S}(\zeta):=\sum_{j=1}^{k}\frac{f_{\mathrm S}^{(j)}(1)}{j!}(\zeta-1)^{j-1},
\qquad f_{\mathrm S}(z)=(\gamma z-\gamma+1)^{k-1}z\ln z,
$$

$$
\varrho^{(k)}_{b,\mathrm S}(\zeta):=\zeta(\gamma\zeta-\gamma+1)^{k-1},
\qquad
\varrho^{(k)}_{c,\mathrm S}(\zeta):=\zeta(\gamma\zeta-\gamma+1)^{k-1}-\gamma^{k-1}(\zeta-1)^{k}.
$$

$\gamma>\frac12$ 时 $\varrho_{b,\mathrm S}^{(k)}$ 的根全在 $|\zeta|<1$ 内；由 Routh-Hurwitz 判据，$\varrho^{(k)}_{c,\mathrm S}$ 的根全在 $|\zeta|<1$ 内的条件对 $k=2,\dots,8$ 依次为

| $k$       | $2$         | $3$                   | $4$        | $5$        | $6$ | $7$       | $8$      |
| --------- | ----------- | --------------------- | ---------- | ---------- | --- | --------- | -------- |
| $\gamma>$ | $-\tfrac12$ | $\tfrac{\sqrt2-1}{2}$ | $\tfrac38$ | $0.658691$ | $1$ | $1.37957$ | $1.7863$ |

SIEMS-2 与 WBDF2／GBDF2 重合。与推广 BDF（只到 $k=5$ 零稳定）不同，**SIEMS-$k$ 在适当的 $\gamma$ 下零稳定且满足引理 2.4 的前提，一直到 $k=8$**，因此理论给出**八阶**的无条件稳定性；作者说未见过阶数高于七的已知无条件稳定隐显多步格式。

**（比较结论）** 论文评估了五个族：$\alpha$ 参数化的加权 BDF、$s$ 参数化的修正 BDF、$\beta$ 参数化的推广 BDF、$\delta$ 参数化的 NIMEX，以及新的 $\gamma$ 参数化 SIEMS。按 $\mathfrak I_{\mathrm{IE}}^{(k)}$ 的理论取值范围，**推广 BDF-$k$（$2\le k\le5$）与 SIEMS-$k$（$2\le k\le6$）对非线性抛物模型的适应性优于其余已有的隐显多步格式。**

### 数值实验

**这一版没有偏微分方程数值实验。** 摘要里那句「给出数值实验以支持理论」在源文件中被注释掉了，结论部分写明「后续报告将说明隐显可控强度的用法」。

论文实际给出的是**三个因子的计算评估**：图形文件是 $\lambda_{\mathrm I}$、$\sigma_{\mathrm E}$、$\sigma_{\mathrm F}$ 随族参数变化的曲线，覆盖 GBDF4、GBDF5、WBDF5 与 SIEMS4 至 SIEMS8。

**因此本文与[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 58 与 74]] 属于同一类**：结论全是定理，没有格式层面的实测数据可以对照。要注意区分——那些因子曲线检验的是**引理的可计算性与各族的相对优劣**，不是格式在具体方程上的收敛阶。上面那句「$2\le k\le5$ 与 $2\le k\le6$ 更优」也是基于 $\mathfrak I_{\mathrm{IE}}^{(k)}$ 的理论取值范围，而非数值算例。

### 与其他论文的关系

它是[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 58]] 的自然续篇：两篇都用基于 DOC 的能量论证取代 Nevanlinna-Odeh 乘子技术，都在对称化的 Toeplitz 型上用 Grenander-Szegő。差别在于编号 58 处理**均匀网格上的全隐 BDF-$k$**，本文处理**真正非线性抛物问题上的隐显多步法**，并把生成函数推广到单边（「半」）版本，好让 $b/a$ 与 $c/a$ 这样的**比值**也能被处理。

它是编号 91 的多步对应物：编号 91 用级指标的微分矩阵 $D_{\mathrm R}$ 与其最小特征值 $\lambda_{\min}$，本文用复合 Toeplitz 矩阵 $\widehat B_{L,k}$ 与其最小特征值 $\lambda_{\mathrm I}^{(k)}$。**两篇的结构完全一样：一个对称化卷积／微分矩阵的特征值，压住显式部分的一个范数界。**

用生成函数判定正定性，是[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 74]] 那套代数判据 C1-C4 的**均匀网格版本**；两者都上溯到 Toeplitz-Carathéodory 与 Grenander-Szegő，但编号 74 逃到了任意非均匀网格上——那里根本没有生成函数可用。

**本文不出现任何步长比常数**，因为它和编号 58 一样是均匀网格的论文。这里相关的「小数字」换成了 $\mu_0/\varpi$（模型给出的门槛）与 $\mathfrak I_{\mathrm{IE}}^{(k)}\le1$（方法给出的强度）。

## 78：线性松弛与正则化能量重构

编号 78 处理的是另一类思路：**不改进时间离散，而改写能量本身。** 论文把自己的方法命名为 RRER（relaxation with regularized energy reformulation，带正则化能量重构的松弛法）。

### 直觉

能量稳定格式的常见做法有凸分裂（精度受限且常常非线性）、指数时间差分、稳定化（通常一阶）、不变能量二次化（IEQ）、标量辅助变量（SAV）、Lagrange 乘子法与补充变量法。IEQ 与 SAV 的共同想法是把非线性能量写成带辅助变量的**二次**形式，于是每步只解线性系统。IEQ 是 Lagrange 乘子法的推广，通常给出**耦合**且系数**依赖时间**的系统；SAV 保留了它的优点，但给出**解耦**且系数为**常数**的系统。

两者也有一处共同的做法，而这正是本文针对的地方：**辅助变量的演化方程是靠对该变量求时间导数得到的。** 那一步引入了它自己的截断误差，也让离散系统对原方程的忠实度下降。论文的原话是它「不需要对辅助变量求时间导数」。

**本文的办法是让辅助变量只由一个代数关系定义，从不参与时间求导。** 光有这一条还不够——如果 $q$ 与 $\phi$ 在同一层上，(c) 式代回去仍会让格式非线性。第二个成分是 Jiang 等人的**松弛**想法：把 $q$ 放到**半整数层**上。这样定义 $q^{n+\frac12}$ 的代数关系右端只含**已知的** $\phi^n$，而 $q^{n+\frac12}$ 出现在动量方程里时又只与 $\frac{\phi^{n+1}+\phi^n}{2}$ 相乘——于是整个格式对未知量是线性的，而且还是二阶的。**交错网格是让「代数定义」与「二阶精度」同时成立的那个技巧。**

### 问题设定

$\Omega\subset\mathbb R^d$（$d=2,3$），周期边界，$(f,g)=\int_\Omega fg\,\mathrm d\mathbf x$、$\|f\|=\sqrt{(f,f)}$。从简化的自由能与其梯度流出发（$\mathcal G\ge0$ 为半正定迁移算子，$\mathcal L$ 线性）：

$$
E(\phi)=\tfrac12(\mathcal L\phi,\phi)+\bigl(F(\phi),1\bigr),
\qquad
\frac{\partial\phi}{\partial t}=-\mathcal G\bigl(\mathcal L\phi+F'(\phi)\bigr),
$$

$$
\frac{\mathrm d}{\mathrm dt}E(\phi)
=\Bigl(\frac{\delta E}{\delta\phi},\frac{\partial\phi}{\partial t}\Bigr)
=-\bigl(\mathcal L\phi+F'(\phi),\ \mathcal G(\mathcal L\phi+F'(\phi))\bigr)\le0 .
$$

**测试模型之一：带斜率选择的分子束外延模型。**

$$
E(\phi)=\int_\Omega\Bigl(\frac{\epsilon^2}{2}(\Delta\phi)^2
+\frac14\bigl(|\nabla\phi|^2-1\bigr)^2\Bigr)\mathrm d\mathbf x,
\qquad \mathcal G=I,
$$

$$
\phi_t=-\epsilon^2\Delta^2\phi+\nabla\cdot\bigl((|\nabla\phi|^2-1)\nabla\phi\bigr),
$$

在 $\partial\phi/\partial\mathbf n|_{\partial\Omega}=0$、$\partial\Delta\phi/\partial\mathbf n|_{\partial\Omega}=0$ 下守质量：$\frac{\mathrm d}{\mathrm dt}\int_\Omega\phi\,\mathrm d\mathbf x=0$。

> [!note] 与编号 52 处理的是**不同**的分子束外延模型
> 本文是**带斜率选择**的模型，自由能里是双阱 $\frac14(|\nabla\phi|^2-1)^2$；[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 52]] 处理的是**无斜率选择**的模型，自由能里是对数项 $-\frac12\ln(1+|\nabla\phi|^2)$。**无斜率选择的情形对能量稳定性而言更难**，因为那里的非线性项不是多项式有界的，自由能也没有下界。编号 52 的结论部分明说其技术「不适用于带斜率选择的分子束外延模型」，而本文用完全不同的方法处理后者。两者不可混为一谈。

**测试模型之二：相场晶体模型（PFC）。**

$$
E(\phi)=\int_\Omega\Bigl(\tfrac12\phi(a_0+\Delta)^2\phi+\tfrac14\phi^4
-\tfrac{b_0}{2}\phi^2\Bigr)\mathrm d\mathbf x,
\qquad
\mathcal G=-\lambda\Delta,\quad 0<b_0<a_0,\ b_0\ll1,\ \lambda>0,
$$

$$
\phi_t=\lambda\Delta\mu,\qquad \mu=(a_0+\Delta)^2\phi+\phi^3-b_0\phi .
$$

**步长比限制：没有**（步长 $\delta t$ 均匀，稳定性无条件）。**最大值原理论证：没有。DOC／DCC 核：没有**——本文与 Liao-Tang-Zhou 的卷积核纲领无关。

### 推导

**第一步：正则化辅助变量。** 取 $C_0$ 使 $F(\phi)+C_0\ge0$，并引入**稳定化参数** $\gamma$：

$$
q=\sqrt{4\bigl(F(\phi)+C_0\bigr)}-\gamma
\qquad\Longrightarrow\qquad
F(\phi)=\tfrac14q^2+\tfrac12\gamma q+\tfrac{\gamma^2}{4}-C_0,
\qquad
F'(\phi)=\tfrac12qq'+\tfrac{\gamma}{2}q',\ \ q'=\tfrac{\mathrm dq}{\mathrm d\phi} .
$$

等价系统与相应能量为

$$
\frac{\partial\phi}{\partial t}=-\mathcal G\mu,
\quad
\mu=\mathcal L\phi+\tfrac12qq'+\tfrac{\gamma}{2}q',
\quad
q=\sqrt{4\bigl(F(\phi)+C_0\bigr)}-\gamma,
$$

$$
\widehat E(\phi,q)=\tfrac12(\mathcal L\phi,\phi)
+\Bigl(\tfrac14q^2+\tfrac12\gamma q+\tfrac{\gamma^2}{4}-C_0,\,1\Bigr),
\qquad
\frac{\mathrm d}{\mathrm dt}\widehat E(\phi,q)=-(\mu,\mathcal G\mu)\le0 .
$$

注 2.1 把整件事说得最清楚：这个重构「与 IEQ 和 SAV 的做法不同，因为我们不对辅助变量求时间导数，而那会在数值计算中引入截断误差」。

**第二步：对分子束外延模型的具体重构。** 取 $q=|\nabla\phi|^2-1-\gamma$（$\gamma>0$），得

$$
\begin{cases}
\phi_t=-\epsilon^2\Delta g+\nabla\cdot(q\nabla\phi)+\gamma g,\\
g=\Delta\phi,\\
q=|\nabla\phi|^2-1-\gamma,
\end{cases}
$$

$$
\widehat E(\phi,q)=\frac{\epsilon^2}{2}\|\Delta\phi\|^2+\frac{\gamma}{2}\|\nabla\phi\|^2
+\frac12\bigl(q(|\nabla\phi|^2-1-\gamma),1\bigr)-\frac14\|q\|^2
-\frac{2\gamma+\gamma^2}{4}|\Omega| .
$$

> [!warning] $\widehat E$ 与 $E$ 的关系在连续层面与离散层面**不同**
> **连续层面：$\widehat E(\phi,q)$ 与 $E(\phi)$ 精确相等，不只是近似**（注 2.2，论文把这个恒等式逐行推了出来）。这与 SAV、IEQ 的修正能量不同，后者与原能量之间只有辅助变量相容意义下的近似关系。
>
> **离散层面：不再精确相等。** 由于 $q$ 住在交错网格的半整数层上，注 2.6 与 2.11 指出离散修正能量只是原能量的**二阶近似**，论文自己说这一点「与 IEQ 或 SAV 的修正能量是一致的」。
>
> 引用时必须区分这两个层面。下面「为什么是本专题的例外」一节的全部张力都系于此。

**第三步：交错时间网格使格式线性（算法一）。** 已知 $(\phi^n,q^{n-\frac12})$，求 $(\phi^{n+1},g^{n+\frac12},q^{n+\frac12})$：

$$
\begin{aligned}
&\text{(a)}\ \ \frac{\phi^{n+1}-\phi^{n}}{\delta t}
=-\epsilon^2\Delta g^{n+\frac12}
+\nabla\cdot\Bigl(q^{n+\frac12}\nabla\frac{\phi^{n+1}+\phi^{n}}{2}\Bigr)
+\gamma g^{n+\frac12},\\
&\text{(b)}\ \ g^{n+\frac12}=\Delta\frac{\phi^{n+1}+\phi^{n}}{2},
\qquad
\text{(c)}\ \ \frac{q^{n+\frac12}+q^{n-\frac12}}{2}=|\nabla\phi^{n}|^2-1-\gamma .
\end{aligned}
$$

**格式线性的原因就在 (c)：它是代数的，且右端只含已知的 $\phi^n$。** 交错（$q$ 在半整数层、$\phi$ 在整数层）加上这一点，使每步只解一个线性代数系统。启动值取 $q^{\frac12}=|\nabla\phi^0|^2-1-\gamma$、$g^{\frac12}=\Delta\frac{\phi^1+\phi^0}{2}$，再走 (a)。时间二阶精度（注 2.3）。

**第四步：相场晶体模型（算法二）。** 取 $q=\phi^2-b_0-\gamma$（$\gamma>0$），得

$$
\text{(a)}\ \phi_t=\lambda\Delta\mu,
\quad
\text{(b)}\ \mu=(a_0+\Delta)g+q\phi+\gamma\phi,
\quad
\text{(c)}\ g=(a_0+\Delta)\phi,
\quad
\text{(d)}\ q=\phi^2-b_0-\gamma,
$$

$$
\widehat E(\phi,q)=\tfrac12\|\Delta\phi\|^2-a_0\|\nabla\phi\|^2
+\tfrac12(\gamma+a_0^2)\|\phi\|^2-\tfrac14\|q\|^2
+\tfrac12\bigl(q(\phi^2-b_0-\gamma),1\bigr)-\tfrac{(\gamma+b_0)^2}{4}|\Omega| .
$$

同样 $\widehat E\equiv E$ 在连续层面精确成立（注 2.7）。算法二是相应的交错 Crank-Nicolson 格式，二阶（注 2.8）且守质量（注 2.12）。

**第五步：耦合模型。** 论文第 2.4 节把 RRER 推广到**三元相场模型**（2.4.1）与**晶粒生长相场模型**（2.4.2）。**这两个模型的确切能量与格式本站未核实**，因此本页只报告推广的存在，不写出其形式。

### 定理

**（定理 2.4，质量守恒，分子束外延）** 算法一保总质量：

$$
\int_\Omega\phi^{n+1}\,\mathrm d\mathbf x=\int_\Omega\phi^{n}\,\mathrm d\mathbf x .
$$

证明是把 (a) 在 $\Omega$ 上积分、代入积分后的 (b)、再用 Green 公式。

**（定理 2.5，无条件能量稳定性，分子束外延）** 对任意 $n\ge1$，

$$
\widehat E^{n+1}\bigl(\phi^{n+1},q^{n+\frac12}\bigr)
-\widehat E^{n}\bigl(\phi^{n},q^{n-\frac12}\bigr)
=-\delta t\,\bigl\|\mu^{n+\frac12}\bigr\|^2\ \le\ 0,
$$

$$
\mu^{n+\frac12}=\epsilon^2\Delta g^{n+\frac12}
-\nabla\cdot\Bigl(q^{n+\frac12}\nabla\frac{\phi^{n+1}+\phi^{n}}{2}\Bigr)
-\gamma g^{n+\frac12},
$$

$$
\widehat E^{n+1}=\frac{\epsilon^2}{2}\|\Delta\phi^{n+1}\|^2
+\frac{\gamma}{2}\|\nabla\phi^{n+1}\|^2
+\frac12\bigl(q^{n+\frac12}(|\nabla\phi^{n+1}|^2-1-\gamma),1\bigr)
-\frac14\|q^{n+\frac12}\|^2-\frac{2\gamma+\gamma^2}{4}|\Omega| .
$$

**注意这是一个等式，不是不等式**——格式精确耗散 $\delta t\|\mu^{n+\frac12}\|^2$。证明把 (a) 与 $\delta t\,\mu^{n+\frac12}$ 配对、把 (c) 与 $q^{n+\frac12}-q^{n-\frac12}$ 配对，后者给出关键恒等式

$$
\tfrac12\bigl(\|q^{n+\frac12}\|^2-\|q^{n-\frac12}\|^2\bigr)
=\bigl((q^{n+\frac12}-q^{n-\frac12})(|\nabla\phi^n|^2-1-\gamma),1\bigr).
$$

**对 $\delta t$ 没有任何限制。**

**（定理 2.10，无条件能量稳定性，相场晶体）** 对任意 $n\ge1$，

$$
\widehat E^{n+1}\bigl(\phi^{n+1},q^{n+\frac12}\bigr)
-\widehat E^{n}\bigl(\phi^{n},q^{n-\frac12}\bigr)
=-\lambda\,\delta t\,\bigl\|\nabla\mu^{n+\frac12}\bigr\|^2\ \le\ 0,
$$

$$
\widehat E^{n+1}=\tfrac12\|\Delta\phi^{n+1}\|^2-a_0\|\nabla\phi^{n+1}\|^2
+\tfrac12(\gamma+a_0^2)\|\phi^{n+1}\|^2-\tfrac14\|q^{n+\frac12}\|^2
+\tfrac12\bigl(q^{n+\frac12}[(\phi^{n+1})^2-b_0-\gamma],1\bigr)
-\tfrac{(\gamma+b_0)^2}{4}|\Omega| .
$$

### 数值实验

空间用标准**有限元**（$P_1$ 元），在 **FreeFEM** 中实现，周期边界。六组算例。

**算例一（收敛性检验，3.1 节）。** $\gamma=2.0$，$\Omega=[0,1]^2$，$T=1$，制造解 $\phi(x,y,t)=e^{-t}\cos(\pi x)\cos(\pi y)$，对分子束外延（$\epsilon=1$）与相场晶体（$a_0=1.0$、$b_0=0.01$、$\lambda=1.0$）分别检验，取 $\delta t=h$。与 IEQ 及指数型 SAV（ESAV）对照。分子束外延一组的实测为：

| $\delta t$  | $\frac18$    | $\frac1{16}$ | $\frac1{32}$ | $\frac1{64}$ |
| ----------- | ------------ | ------------ | ------------ | ------------ |
| RRER 误差   | $3.2923$e-03 | $8.2307$e-04 | $2.0577$e-04 | $5.1442$e-05 |
| RRER 观测阶 | —            | $2.00$       | $2.00$       | $2.00$       |
| IEQ 观测阶  | —            | $1.94$       | $1.97$       | $1.99$       |

ESAV 与 RRER 相当。**三个方法都是二阶**；在相场晶体一组里 RRER 比 IEQ 与 ESAV 都更准，且 CPU 时间表显示 RRER 最省。

**算例二（能量耗散与质量守恒，3.2 节）。** 初值 $\phi(x,y,0)=\cos\pi x\cos\pi y$；在若干 $\delta t$ 与 $\epsilon$ 下同时画出原始能量与修正能量。两个模型在所测的一切 $\delta t$ 下都单调下降（支持无条件稳定性），质量守恒（分子束外延取 $\gamma=20$、$\epsilon=1.0$；相场晶体取 $\gamma=1.0$、$a_0=1.0$、$b_0=0.325$、$\delta t=0.1$）。

**算例三（二维分子束外延粗化，3.3 节）。** $\Omega=[0,2\pi]^2$，$\epsilon^2=0.1$，$\gamma=20$，$\phi(x,y,0)=0.1(\sin3x\sin2y+\sin5x\sin5y)$，$h=2\pi/128$，$\delta t=10^{-4}$；$t=0,0.05,2.5,8,15,30$ 处的形貌重现了**带斜率选择**的分子束外延已知相图。

**算例四（二维相场晶体，3.4 节）。** $\Omega=[0,100]^2$，$\phi(x,y,0)=\hat\phi_0+0.01\,\mathrm{rand}(x,y)$，$a_0=1.0$，$b_0=0.35$，$T=200$，$\gamma=2.0$，$h=100/128$，$\delta t=0.1$：$\hat\phi_0=0$ 给出**条纹**图样，$\hat\phi_0=0.2$ 给出**三角**图样，与文献相图一致。

**算例五（曲面上的相场晶体，3.5 节）。** 半径 $R=64$ 的球面（依 $\hat\phi_0$ 不同，$t=50,100,200$ 处给出条纹或六角图样）；外半径 $50$、内半径 $20$ 的环面（$h=1$、$\delta t=0.1$，条纹对六角）；以及 $[0,100]^3$ 上的立方体表面情形（$h=100/64$、$\delta t=0.1$）。

**算例六（三维相场晶体与耦合系统，3.6 节）。** 另有一个三元系统算例，取 $m=3$、$L_1=L_2=L_3=1$、$\alpha=\beta=\gamma=1$、$k_1=k_2=k_3=2$，区域 $[0,1]^2$。

**这批实验建立了什么、没建立什么。** 建立的是：格式确实二阶（且时间精度不逊于 IEQ、ESAV），在所测的一切步长下能量单调下降与质量守恒，形貌与文献相图一致，且 CPU 成本最低。没有建立的是**大步长下的定量精度**——无条件稳定意味着不会爆，不意味着大步长下准确；算例二只报告了能量的单调性，没有报告不同 $\delta t$ 之间解的差异。另外，$\gamma$ 作为稳定化参数在不同算例中取值差别很大（$1.0$、$2.0$、$20$），论文没有给出选取 $\gamma$ 的准则。

### 为什么这一篇是本专题的例外

值得把这一点写清楚，因为它关系到整个专题的内在张力。编号 78 没有 Liao 的合作、没有 DOC/DCC 核、没有变步长、也没有步长比分析——它属于 IEQ/SAV/松弛这一支，而不是卷积核那一支。周涛是唯一的交集。

更要紧的是**能量陈述的类型不同**：编号 78 证明的是**修正**能量 $\widehat E(\phi,q)$ 的耗散，而[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 40、57]] 与编号 48、52、91 恰恰在设法**避免**这一点，它们要的是原能量或变分能量的耗散。这个张力在文献里是明写着的：**编号 91 的引言正是批评基于 SAV 的高阶格式只能建立「关于含辅助变量的修正能量」的稳定性。**

编号 78 的缓解之处在于 $\widehat E\equiv E$ 在**连续层面精确成立**，因此它的修正能量不是一个新对象，而是同一个能量的另一种写法。但上面那条警告也必须一并记住：在**离散层面**，由于 $q$ 住在交错层上，$\widehat E$ 退化为原能量的二阶近似，论文自承这与 IEQ／SAV 的情形一致。**因此这不是一个把张力取消掉的论证，而是一个把它缩小到二阶的论证。**

此外，本文与编号 91 都处理相场晶体／Cahn-Hilliard 型的四阶动力学，但稳定化哲学完全不同：这里是辅助变量二次化，那里是显式稳定化加 Runge-Kutta。

## 三篇的对策对照

| 编号 | 被显式处理的对象      | 承载困难的代数对象            | 保住能量论证的手段                               | 衰减的能量                                       | 数值实验                       |
| ---- | --------------------- | ----------------------------- | ------------------------------------------------ | ------------------------------------------------ | ------------------------------ |
| 78   | 非线性势（经二次化）  | 无（改的是能量而非离散）      | 重构能量，交错网格使格式线性                     | **修正**能量 $\widehat E$（连续层面 $\equiv E$） | 六组，二阶，相图一致           |
| 91   | 非线性主体 $f_\kappa$ | 微分矩阵 $D(z)$（级指标 DOC） | 令 $D_{\mathrm{EI}}=\mathbf 0$，耗散率与网格无关 | **原始**能量 $E$                                 | 五组，含三个数量级的加速表     |
| 104  | 非线性项 $\mathcal F$ | 三串卷积核的复合              | 三个半生成函数在单位圆上的极值                   | 抽象框架的 $V$ 范数稳定性                        | 无偏微分方程算例；仅三因子曲线 |

**读这张表的方式**：三条路线都不试图直接证明「显式项贡献符号确定」，而是各自换一个可控对象——编号 78 换能量，编号 91 换 Butcher 表的结构条件，编号 104 换判据的形式。真正可比的是「衰减的能量」一列：**编号 91 与 104 留在原能量／原范数一侧，编号 78 走到修正能量一侧**，这条分界线也正是本专题内部唯一一处方法论上的真分歧。

## 覆盖核对

| 内容                                 | 论文 | 覆盖状态                                                                                                 |
| ------------------------------------ | ---- | -------------------------------------------------------------------------------------------------------- |
| 显式处理为何破坏能量论证             | 全页 | 开篇一节：同号性失效与三种替代对象                                                                       |
| Cahn-Hilliard 能量与 $H^{-1}$ 梯度流 | 91   | 能量、方程、耗散律、体积守恒                                                                             |
| 两个缺陷                             | 91   | 全局 Lipschitz 假设失效（且 CH 无最大值原理）；耗散率随步长漂移                                          |
| 稳定化、节点条件与阶条件             | 91   | $L_\kappa$、$\hat{\mathbf c}=\mathbf c$、三阶以内含耦合条件                                              |
| 差分形式与微分矩阵                   | 91   | $D(z)=D_{\mathrm E}-zD_{\mathrm{EI}}$ 及其正交性恒等式                                                   |
| $D$ 与 DOC 核的对应                  | 91   | 索引由时间层换成级                                                                                       |
| 平均耗散率与精化条件                 | 91   | $\mathcal R$、$D_{\mathrm{EI}}=\mathbf 0$、$A_{\mathrm I}=A_{\mathrm E}P$、$\mathcal R_{\mathrm R}$      |
| Lobatto 型的必然性                   | 91   | 第一列非零、排除性命题、级阶二且不必代数稳定                                                             |
| 三族具体方法与参数窗口               | 91   | R-IERK(1,2)/(2,4;$c_2$)/(3,6;$\hat a_{52}$)、两个行列式、$c_2=1$ 被排除                                  |
| $\frac{1+\sqrt2}{4}$ 不是步长比      | 91   | 它是对照方法 IERK(2,3) 的 Butcher 系数 $a_{33}$                                                          |
| 替代最大值原理的时空分裂             | 91   | 三步自举、引理 4.3 与 4.4、各级能量不等式                                                                |
| 四条主要结果                         | 91   | 各级能量律、时间半离散正则性、原始能量耗散、无条件 $L^2$ 收敛                                            |
| 数值实验与效率                       | 91   | 自适应规则、五组算例、CPU 表（$6724.63$ s 对 $2.62$ s）                                                  |
| 现有两条路线的具体局限               | 104  | 乘子逐族寻找；分解在 $\beta=1$ 退化、$k$ 受限、精化版仅 $\beta_k=3,6,9$                                  |
| 抽象设定与局部 Lipschitz 条件        | 104  | Hilbert 三元组、$\mu_0\in(0,\varpi)$、非自伴处理                                                         |
| 三组核与全局离散能量方法             | 104  | 三元组、DOC 核、复合核、下三角 Toeplitz 可交换                                                           |
| 半生成函数及其证明观察               | 104  | 三条结论、$\mathrm g(\theta)=\Re[a(\theta)]$、复合规则                                                   |
| 三个极值常数与归一化                 | 104  | $\sigma_{\mathrm F},\sigma_{\mathrm E},\lambda_{\mathrm I}$、三条不等式、隐显 Euler 取等                 |
| 主定理与可控强度                     | 104  | $\lambda_{\mathrm I}/\sigma_{\mathrm E}>\mu_0/\varpi$、方法与模型的分离、$\mathfrak I_{\mathrm{IE}}\le1$ |
| SIEMS 族与八阶                       | 104  | 三个特征多项式、$k=2,\dots,8$ 的 $\gamma$ 门槛表、比较结论                                               |
| 「本版无数值实验」                   | 104  | 摘要那句被注释掉；仅有三因子曲线                                                                         |
| IEQ／SAV 的共同做法与 RRER 的分歧    | 78   | 耦合对解耦、对辅助变量求时间导数、注 2.1                                                                 |
| 正则化辅助变量与等价系统             | 78   | $q$ 的定义、$F$ 的二次化、$\widehat E$ 与其耗散律                                                        |
| $\widehat E$ 与 $E$ 在两个层面的关系 | 78   | 连续层面精确相等；离散层面退为二阶近似                                                                   |
| 交错网格与格式的线性性               | 78   | 算法一三式、(c) 为何是关键、启动值、二阶                                                                 |
| 相场晶体与耦合模型                   | 78   | 算法二与其修正能量；三元／晶粒生长仅报告存在                                                             |
| 三条定理                             | 78   | 质量守恒、两条无条件能量稳定（等式形式）与证明中的关键恒等式                                             |
| 六组数值实验                         | 78   | 收敛表、能量与质量、粗化、条纹／三角、曲面、三维与三元                                                   |
| 与编号 52 的模型差别                 | 78   | 带斜率选择对无斜率选择，及后者更难的原因                                                                 |
| 本专题的例外与那条张力               | 78   | 无 Liao／无 DOC／无变步长；修正能量对原能量之争                                                          |

## 本页原文

- J. Zhang, X. Guo, M. Jiang, T. Zhou, and J. Zhao, [_Linear relaxation method with regularized energy reformulation for phase field models_](https://doi.org/10.1016/j.jcp.2024.113225), J. Comput. Phys. 515 (2024), 113225。（卷号、年份与文章号已核实；**DOI 字符串未在 PDF 提取文本中出现**，此处链接按通用格式给出。）
- H.-l. Liao, T. Tang, X. Wang, and T. Zhou, [_A class of refined implicit-explicit Runge-Kutta methods with robust time adaptability and unconditional convergence for the Cahn-Hilliard model_](https://doi.org/10.1090/mcom/4090), Math. Comput. 95(359) (2026), pp. 1293-1325（预印本 [arXiv:2412.07321](https://arxiv.org/abs/2412.07321)）。**卷期页码取自阅读清单，本站未与 AMS 页面核对**；正文内容依 arXiv 源文件核实。
- H.-l. Liao, C. Quan, T. Tang, and T. Zhou, _A semi-generating function approach to the stability of implicit-explicit multistep methods for nonlinear parabolic equations_, [arXiv:2605.05619](https://arxiv.org/abs/2605.05619)，投稿 SIAM J. Numer. Anal.
