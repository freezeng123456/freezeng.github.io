---
title: 时间分数阶相场
description: 编号 40、43、57：分数阶梯度流的能量律从积分型走向微分型
lang: zh
translation: en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field
tags:
  - 论文笔记
  - 相场模型
  - 分数阶导数
---

> [!note] 本页覆盖
> 编号 **40**（_SIAM J. Sci. Comput._ 41(6), 2019）、**43**（_J. Comput. Phys._ 414, 2020）、**57**（_SIAM J. Sci. Comput._ 43(5), 2021）。三篇构成一条清晰的推进线：编号 40 证明积分型能量律，编号 57 把它换成微分型的变分能量律。

## 40：能量耗散在分数阶情形是积分型的

### 标准论证在哪里失效

相场模型（Allen-Cahn、Cahn-Hilliard、分子束外延）都来自一个随时间耗散的自由能，而这条耗散律是长时间数值稳定性的支柱。把时间导数换成阶为 $\alpha\in(0,1)$ 的 Caputo 导数

$$
\frac{\partial^{\alpha}}{\partial t^{\alpha}}\phi
={}^{C}_{0}D_{t}^{\alpha}\phi(t)
:=\frac{1}{\Gamma(1-\alpha)}\int_{0}^{t}\frac{\phi'(s)}{(t-s)^{\alpha}}\,\mathrm ds
$$

后，标准论证失效：乘以 $\phi_t$ 不再逐点给出符号确定的项，因为 $\partial_t^\alpha\phi$ 是整个历史的**时间非局部**卷积。在本文之前，时间分数阶相场模型的能量耗散行为只在数值上被观察到，无论连续层面还是离散层面都没有证明。

本文补上这个缺口，并且指出了这条律的**正确形式**：它是**积分型**的（在 $[0,T]$ 上累积），而不是逐点的微分型。这一点值得记住，因为它正是编号 57 的出发点。

三个模型的具体形式为

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}
=\gamma\Bigl(\varepsilon\Delta\phi-\frac{1}{\varepsilon}F'(\phi)\Bigr)
\qquad\text{（Allen-Cahn）},
$$

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}=\gamma\Delta\mu,
\quad \mu=-\varepsilon\Delta\phi+\frac{1}{\varepsilon}F'(\phi)
\qquad\text{（Cahn-Hilliard）},
$$

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}
=\gamma\Bigl(-\varepsilon\Delta^{2}\phi
+\frac{1}{\varepsilon}\nabla\cdot\mathbf f_m(\nabla\phi)\Bigr)
\qquad\text{（分子束外延）},
$$

对应的自由能为

$$
E[\phi]=\frac{\varepsilon}{2}\|\nabla\phi\|^{2}
+\frac{1}{\varepsilon}\bigl\langle F(\phi),1\bigr\rangle,
\qquad F(\phi)=\tfrac14(1-\phi^{2})^{2},
$$

$$
E_{m}[\phi]=\frac{\varepsilon}{2}\|\Delta\phi\|^{2}
+\frac{1}{\varepsilon}\bigl\langle F_m(\nabla\phi),1\bigr\rangle,
$$

其中分子束外延模型的 $F_m$ 按是否带斜率选择分两种：$\frac14(|\mathbf v|^{2}-1)^{2}$ 与 $-\frac12\ln|1+|\mathbf v|^{2}|$，相应地 $\mathbf f_m(\nabla\phi)=(|\nabla\phi|^{2}-1)\nabla\phi$ 与 $-\nabla\phi/(1+|\nabla\phi|^{2})$。

### 唯一的分析装置：分数阶卷积核的正定性

替代逐点符号论证的是一条关于卷积核的正定性结论。对 $h,g\in L^{p}(0,T)$ 定义

$$
I_\alpha(h,g):=\frac{1}{\Gamma(\alpha)}\int_{0}^{T}\!\!\int_{0}^{t}
\frac{h(s)g(t)}{(t-s)^{1-\alpha}}\,\mathrm ds\,\mathrm dt .
$$

论文证明：对 $\alpha\in(0,1)$ 与 $p\ge\frac{2}{1+\alpha}$，核 $I_\alpha(h,h)$ 具有正定性。这一条把「历史的加权累积」变成一个有符号的量，从而使能量论证可以在**积分**层面完成。

这个装置的意义超出本文：它是这一系列工作中反复出现的模式的连续版本，离散版本就是[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|变步长 BDF 一页]]中的二次型正定性。**连续与离散层面的困难与对策是同构的：都要把非局部的历史项汇成一个正定的量。**

## 43：非均匀步长下的二阶保最大值原理格式

编号 43 处理时间分数阶 Allen-Cahn 方程的数值格式，目标是同时具备二阶精度、非均匀步长与离散最大值原理。非均匀步长在这里不是可选项：分数阶问题在 $t=0$ 附近通常有初始奇性，需要分级网格；而随后的粗化过程极慢，需要放大步长。

这一篇在整条线索中的位置是承接的：它把编号 40 的连续层面结论带到离散层面，并处理最大值原理；而它留下的问题——能量律在离散层面是什么形式——由编号 57 回答。

## 57：把提问方式改掉

### 为什么编号 40 的结论不够用

编号 40 证明了 $E[u](t)\le E[u](0)$ 与最大值原理，但这个能量表述弱于经典耗散律

$$
\frac{\mathrm dE}{\mathrm dt}+\Bigl\|\frac{\delta E}{\delta u}\Bigr\|^2=0,
$$

而且不是**微分**形式的律，因此不能像经典律那样使用。同时实践需要**变**步长（初始奇性加粗化），而当时没有格式能在变步长下**同时**保住能量稳定性与最大值原理。

本文的解决方式是改变提问：不问**原**能量是否耗散，而定义一个确实耗散的**变分能量**。

### 关键一步：把非局部性移到变分导数上

模型为 $\partial_t^{\alpha}u=\varepsilon^{2}\Delta u-f(u)$，$F(u)=\tfrac14(1-u^2)^2$，$f=F'$，读作**分数阶梯度流** $\partial_t^{\alpha}u=-\delta E/\delta u$。用 Riemann-Liouville 导数 ${}^{R}\!\partial_t^{\alpha}v:=\partial_t\mathcal I_t^{1-\alpha}v$ 与半群恒等式

$$
{}^{R}\!\partial_t^{1-\alpha}\bigl(\partial_t^{\alpha}v\bigr)
=\partial_t\mathcal I_t^{1}v'=v',
$$

方程被**等价改写**为

$$
\partial_t u=-{}^{R}\!\partial_t^{1-\alpha}\Bigl(\frac{\delta E}{\delta u}\Bigr).
$$

这一步是全篇的核心：左端的时间导数现在是**局部**的，因此标准能量测试函数可用；非局部性被移到变分导数上，而那里有一条可用的正定性。

### 变分能量与其微分律

利用 Riemann-Liouville 不等式

$$
v(t)\bigl({}^{R}\!\partial_t^{1-\alpha}v\bigr)(t)
\ \ge\ \tfrac12\bigl({}^{R}\!\partial_t^{1-\alpha}v^{2}\bigr)(t)
+\tfrac12\omega_{\alpha}(t)v^{2}(t),
\qquad \omega_{\mu}(t)=\frac{t^{\mu-1}}{\Gamma(\mu)},
$$

定义变分能量并证明其律：

$$
\mathcal{E}_{\alpha}[u]:=E[u]+\frac12\,\mathcal{I}_t^{\alpha}
\Bigl\|\frac{\delta E}{\delta u}\Bigr\|^{2},
\qquad
\frac{\mathrm d\mathcal{E}_{\alpha}}{\mathrm dt}
+\frac12\,\omega_{\alpha}(t)\Bigl\|\frac{\delta E}{\delta u}\Bigr\|^{2}\le0,
\qquad \forall t>0 .
$$

$\alpha\to1$ 时 $\mathcal I_t^\alpha\to\mathcal I_t^1$ 且 $\omega_\alpha(t)\to1$，于是恢复经典耗散律。这条性质论文称为**渐近保能量耗散**：新律不是对经典律的类比，而是它的推广，且在极限处退化回去。

### L1$_R$ 公式与不需要步长比限制的正定性

变步长下把 Riemann-Liouville 导数离散为

$$
\bigl({}^{R}\!\partial_{\tau}^{1-\alpha}v\bigr)^{n-\frac12}
:=\frac{1}{\tau_{n}}\int_{t_{n-1}}^{t_{n}}\frac{\partial}{\partial t}
\int_{0}^{t}\omega_{\alpha}(t-s)(\Pi_{0}v)(s)\,\mathrm ds\,\mathrm dt
\ \triangleq\ \frac{1}{\tau_{n}}\sum_{k=1}^{n}a_{n-k}^{(n)}v^{k-\frac12},
$$

其中 $\Pi_0v$ 是在 $(t_{k-1},t_k]$ 上取值 $v^{k-\frac12}$ 的分段常数插值。辅助序列与核为

$$
q_{n-k}^{(n)}:=\int_{t_{k-1}}^{t_{k}}\omega_{\alpha}(t_{n}-s)\,\mathrm ds
=\sum_{j=k}^{n}a_{j-k}^{(j)}>0,
$$

$$
a_{0}^{(n)}:=q_{0}^{(n)}>0\ (n\ge1),
\qquad
a_{n-k}^{(n)}:=q_{n-k}^{(n)}-q_{n-k-1}^{(n-1)}<0\ (n\ge k+1\ge2).
$$

符号模式值得注意：第一个核为正，之后全为负。显式地 $a^{(n)}_0=\omega_{1+\alpha}(\tau_n)=\tau_n^{\alpha}/\Gamma(1+\alpha)$。

核的正定性在离散层面单独证明（不是从连续核继承）：对任意实序列 $\{w_k\}$，

$$
2\sum_{k=1}^{n}w_{k}\sum_{j=1}^{k}a_{k-j}^{(k)}w_{j}
\ \ge\ \sum_{k=1}^{n}\Bigl(q_{n-k}^{(n)}+\sum_{j=1}^{k}a_{k-j}^{(k)}\Bigr)w_{k}^{2}
\ >\ 0,
\qquad n\ge1,\ w\not\equiv0 .
$$

证明由两条恒等式承担：$\omega_\alpha$ 的完全单调性给出 $q_{k-j-1}^{(k-1)}-q_{k-j}^{(k)}>0$，以及

$$
\sum_{j=1}^{k}a_{k-j}^{(k)}
=\int_{t_{k-1}}^{t_{k}}\omega_{\alpha}(s)\,\mathrm ds>0 .
$$

**这条正定性不需要任何步长比限制。** 这一点与[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|变步长 BDF]] 形成鲜明对照：那里的二次型正定性给出 $r_k<(3+\sqrt{17})/2$ 之类的门槛，而这里的分数阶核由完全单调性直接保证正定。步长上界只在最大值原理的论证中出现，与能量律无关。

### 格式与两个构造性设计

设 $v:=-\delta E/\delta u$，把方程分裂为

$$
\partial_t u={}^{R}\!\partial_t^{1-\alpha}v,
\qquad
v=\varepsilon^{2}\Delta u-f(u),
$$

离散为 Crank-Nicolson 型

$$
\partial_{\tau}u^{n-\frac12}=\bigl({}^{R}\!\partial_{\tau}^{1-\alpha}v\bigr)^{n-\frac12},
\qquad
v^{n-\frac12}=\varepsilon^{2}D_{h}u^{n-\frac12}-H(u^{n},u^{n-1}),
$$

其中非线性项的二阶逼近取

$$
H(u^{n},u^{n-1}):=\tfrac13(u^{n})^{.3}+\tfrac12(u^{n-1})^{.2}\!\circ u^{n}
+\tfrac16(u^{n-1})^{.3}-\tfrac12\bigl(u^{n}+u^{n-1}\bigr)
$$

（$\circ$ 与幂为逐元素运算）。这个特定的 $H$ 是**设计出来的**，使 $H(a,b)(a-b)\ge F(a)-F(b)$ 逐点成立——这是链式法则的离散对应物，也是能量论证得以进行的原因。

离散变分能量为

$$
\mathcal{E}_{\alpha}[u^{n}]:=E[u^{n}]+\frac12h^{2}\sum_{i,j}\sum_{k=1}^{n}
q_{n-k}^{(n)}\bigl(v_{ij}^{k-\frac12}\bigr)^{2},
$$

其中 $q$ 核本身构成一个数值分数阶积分：$(\mathcal I_\tau^{\alpha}v)^n=\sum_{k}q_{n-k}^{(n)}v^{k-\frac12}$。

### DOC 核与两种分数阶导数之间的可逆变换

论文的另一个「首次」是给出两种分数阶导数离散之间的显式等价。定义 L1$_R$ 核的离散正交卷积核 $\theta$ 由正交性恒等式确定：

$$
\sum_{j=k}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}\equiv\delta_{nk},
\qquad
\theta_{0}^{(n)}=\frac{1}{a^{(n)}_{0}} .
$$

**双向**正交性同时成立（这就是可逆性）：

$$
\sum_{j=k}^{n}a^{(n)}_{n-j}\theta_{j-k}^{(j)}\equiv\delta_{nk},
\qquad
\sum_{j=k}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}\equiv\delta_{nk},
$$

并且 $\theta$ 与 $q$ 之间有互补性 $\sum_{j=k}^{n}q_{n-j}^{(n)}\theta_{j-k}^{(j)}\equiv1$。论文证明 $\theta^{(n)}_0=\Gamma(1+\alpha)\tau_n^{-\alpha}$、全部 $\theta^{(n)}_j>0$，且对 $n\ge2$ 有单调递减 $\theta_{0}^{(n)}>\theta_{1}^{(n)}>\cdots>\theta_{n-1}^{(n)}>0$。

把格式的第一个方程与 $\theta^{(n)}_{n-j}$ 卷积并用正交性，得到

$$
\sum_{j=1}^{n}\theta_{n-j}^{(n)}\nabla_{\tau}u^{j}
=\sum_{k=1}^{n}v^{k-\frac12}\sum_{j=k}^{n}\theta_{n-j}^{(n)}a_{j-k}^{(j)}
=v^{n-\frac12},
$$

因此 Riemann-Liouville 形式的 Crank-Nicolson 格式**等价于** Caputo 型形式

$$
\sum_{j=1}^{n}\theta_{n-j}^{(n)}\nabla_{\tau}u^{j}
=\varepsilon^{2}D_{h}u^{n-\frac12}-H(u^{n},u^{n-1}) .
$$

也就是说 DOC 核定义了一个**新的**离散 Caputo 导数，其核在非均匀网格上正且单调递减，与经典 L1 核同性质。论文自己提醒这是一个**间接**逼近，其精度与直接 L1 公式（误差阶 $2-\alpha$）不同。

最大值原理的论证正是在这个 Caputo 等价形式上、按归纳法进行的，用到 $D_h$ 的对称负半定性与一条 $\ell^\infty$ 引理。

## 三篇的推进关系

| 编号 | 能量律的形式             | 步长限制       | 主要分析装置               |
| ---- | ------------------------ | -------------- | -------------------------- |
| 40   | 积分型（$[0,T]$ 上累积） | 连续层面不涉及 | 分数阶卷积核的正定性       |
| 43   | 离散层面的最大值原理     | 非均匀步长     | 二阶非均匀离散             |
| 57   | 微分型（变分能量）       | 能量律无限制   | 方程改写 + L1$_R$ 核正定性 |

编号 40 到编号 57 的推进值得单独总结：**当一条定律的形式不对时，不要削弱结论，而要换一个对象。** 编号 40 得到的是原能量的积分型不等式；编号 57 不去加强它，而是构造一个新的能量 $\mathcal E_\alpha$，使其满足微分型律并在 $\alpha\to1$ 时退化回经典律。代价是这个能量含一个分数阶积分项，因此不是原能量本身。

## 覆盖核对

| 内容                         | 论文 | 覆盖状态                                                         |
| ---------------------------- | ---- | ---------------------------------------------------------------- |
| 三个分数阶相场模型与其自由能 | 40   | Allen-Cahn、Cahn-Hilliard、分子束外延两型                        |
| 标准能量论证为何失效         | 40   | 逐点符号性丢失与非局部卷积                                       |
| 分数阶卷积核正定性           | 40   | 定义、条件与其在论证中的作用                                     |
| 积分型而非微分型             | 40   | 结论形式及其对后续工作的影响                                     |
| 二阶非均匀保最大值原理格式   | 43   | 目标与在线索中的位置                                             |
| 方程改写与非局部性的移位     | 57   | 半群恒等式、改写式、为何标准测试函数可用                         |
| 变分能量与微分型律           | 57   | Riemann-Liouville 不等式、$\mathcal E_\alpha$、$\alpha\to1$ 极限 |
| L1$_R$ 核与其正定性          | 57   | 核定义、符号模式、离散正定性、无步长比限制                       |
| 非线性项的构造性设计         | 57   | $H$ 的形式与 $H(a,b)(a-b)\ge F(a)-F(b)$                          |
| DOC 核与可逆变换             | 57   | 正交性、互补性、单调性、等价的 Caputo 形式                       |

## 本页原文

- T. Tang, H. Yu, and T. Zhou, [_On energy dissipation theory and numerical stability for time-fractional phase-field equations_](https://doi.org/10.1137/18M1203560), SIAM J. Sci. Comput. 41(6) (2019), pp. A3757-A3778（预印本 [arXiv:1808.01471](https://arxiv.org/abs/1808.01471)）。
- H.-l. Liao, T. Tang, and T. Zhou, [_A second-order and nonuniform time-stepping maximum-principle preserving scheme for time-fractional Allen-Cahn equations_](https://doi.org/10.1016/j.jcp.2020.109473), J. Comput. Phys. 414 (2020), 109473。
- H.-l. Liao, T. Tang, and T. Zhou, [_An energy stable and maximum bound preserving scheme with variable time steps for time fractional Allen-Cahn equation_](https://doi.org/10.1137/20M1384105), SIAM J. Sci. Comput. 43(5) (2021), pp. A3503-A3526（预印本 [arXiv:2012.10740](https://arxiv.org/abs/2012.10740)）。
