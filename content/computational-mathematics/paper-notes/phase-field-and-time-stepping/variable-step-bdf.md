---
title: 变步长 BDF 与卷积核
description: 编号 48、52、58、67、69、74：把符号不定的多步核换成非负递减核，能量论证才成为可能
lang: zh
translation: en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf
tags:
  - 论文笔记
  - 相场模型
  - 离散能量
---

## 贯穿这一线索的一个机制

在读任何一篇之前，值得先把这条线索的共同困难与共同解法说清楚，因为六篇论文都是它的变奏。

**困难。** 均匀步长下的 BDF2 是一个 Toeplitz 卷积：核 $b_j$ 只依赖滞后指标 $j$，于是生成函数 $\hat b(z)=\sum_jb_jz^j$、Grenander-Szegő 定理、Toeplitz-Carathéodory 判据这一整套经典工具都可以用。变步长一来，核变成 $b^{(n)}_j$——多了一个上标，记录「当前在第几层」。它不再是 Toeplitz 的，算子也不再自伴，生成函数机器整个失效。

**更具体的困难在于符号。** 变步长 BDF2 的两个核是

$$
b^{(n)}_0=\frac{1+2r_n}{\tau_n(1+r_n)}>0,
\qquad
b^{(n)}_1=-\frac{r_n^2}{\tau_n(1+r_n)}<0 .
$$

第二个是**负的**。这一个负号同时挡住了两条路：

- **能量论证**要的是二次型 $\sum_kw_k\sum_jb^{(k)}_{k-j}w_j$ 的正定性。核有负项，正定性不再自动成立，只在负项被正项压住时才成立——这就是门槛 $r_k<(3+\sqrt{17})/2$ 的由来。
- **最大值原理**要的是把格式写成「当前层的正定算子」$=$「历史层的**非负**组合」。有了非负权重，归纳假设 $\|u^k\|_\infty\le1$ 才能把右端压住，因为非负权重的加权平均不会超过被平均量的最大值。**权重一旦带负号，这个论证就彻底失效**：负系数把一个上界变成没有上界。

**解法。** 有三种：

1. **核重排（KRC）**——换一个变量 $\bar v^k=v^k-\eta v^{k-1}$，使新核 $d^{(n)}_j$ **非负且单调递减**。这是编号 48 为最大值原理走的路，可行范围恰好是 $r_k<1+\sqrt2$。
2. **离散正交卷积核（DOC）**——构造 $\theta$ 使 $\Theta B=I$，把多步算子从格式上「剥掉」，剩下一个单纯的一阶差分。这是编号 52、58、67、74 走的路，把多步稳定性问题化归为二次型正定性问题。
3. **离散互补卷积核（DCC）**——构造 $Q_d$ 使行和恒为 $1$，用来造离散 Grönwall 不等式。编号 48 用它，[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 43]] 在分数阶一侧先用它。

一句话对照：**DOC 求逆（$\equiv\delta_{nk}$，产出能量与 $L^2$ 估计），DCC 求补（$\equiv1$，产出 Grönwall 与 $\ell^\infty$ 估计），KRC 换基（把核变成非负递减，产出最大值原理）。**

## 48：两个步长比常数，两条不同的性质

### 直觉

Allen-Cahn 方程在连续层面同时具有能量耗散与最大值原理，而解跨越两个差别极大的时间尺度——初期演化快、随后粗化极慢——因此实践上必须用变步长。但在本文写作时，即便对线性与半线性抛物方程，变步长格式的数值分析用作者的话说仍「远未完成」。具体的缺口有两处：没有人对 Allen-Cahn 方程的二阶变步长格式证明过**离散最大值原理**；而已有的非均匀 BDF2 能量稳定性结果所依赖的步长比条件，其尖锐性也不清楚。

本文的关键判断是：**能量稳定性与最大值原理需要核的两种不同性质，因此得到两个不同的门槛，二者不能互换。**

能量那一侧要的是**正定性**。用 $\nabla_\tau u^n$ 去测试格式，左端出现 $\sum_k\nabla_\tau u^k\,D_2u^k$，这是一个以 BDF2 核为系数的二次型。核里有负项，所以要问的是负项什么时候被正项压住。答案是 $r_k<(3+\sqrt{17})/2\approx3.561$。

最大值原理那一侧要的是**非负性加单调性**，这是更强的要求，理由如上：$\ell^\infty$ 的归纳只对非负权重有效。BDF2 的核不满足，于是本文做一次**变量替换**：令 $\bar v^k=v^k-\eta v^{k-1}$。这个替换的效果是把一个符号不定的两项核换成一个几何衰减的非负核族 $d^{(n)}_j$，代价是 $\eta$ 必须落在一个区间里，而该区间非空当且仅当 $r_k<1+\sqrt2$。

**这就是为什么 $1+\sqrt2$ 与 $3.561$ 不是「同一个结论的两个版本」，而是两个不同问题的答案。** 更妙的是，$1+\sqrt2$ 恰好就是 Grigorieff（1983）对常微分方程给出的变步长 BDF2 经典**零稳定性**界，而本文是在这个经典极限上**恰好**拿到离散最大值原理与最大范数二阶收敛的——此前的 $L^2$ 分析都需要比 $1+\sqrt2$ 更苛刻的限制。

### 问题设定

$\Omega=(0,L)^2$，周期边界，

$$
\partial_{t}u=\varepsilon^{2}\Delta u-f(u),\quad f(u)=u^{3}-u,
\qquad u(\mathbf x,0)=u_{0}(\mathbf x),
$$

它是 Ginzburg-Landau 自由能的 $L^2$ 梯度流，

$$
E[u](t):=\int_{\Omega}\Bigl(\tfrac12\varepsilon^{2}|\nabla u|^{2}+F[u]\Bigr)\mathrm d\mathbf x,
\qquad F[u]=\tfrac14(1-u^{2})^{2},
$$

因而 $\mathrm dE/\mathrm dt\le0$，且 $|u(\mathbf x,0)|\le1$ 蕴含 $|u(\mathbf x,t)|\le1$。

网格与步长比（这是整条文献的标准记号）：

$$
0=t_0<t_1<\cdots<t_N=T,
\qquad \tau_k:=t_k-t_{k-1},
\qquad \tau:=\max_{1\le k\le N}\tau_k,
$$

$$
r_{k}:=\frac{\tau_{k}}{\tau_{k-1}}\ (2\le k\le N),
\qquad r_1\equiv0\ \text{（约定）} .
$$

以 $\Pi_{n,2}v$ 记 $v$ 在 $t_{n-2},t_{n-1},t_n$ 上的二次插值，变步长 BDF2 公式为

$$
D_{2}v^{n}:=(\Pi_{n,2}v)'(t_{n})
=\frac{1+2r_{n}}{\tau_{n}(1+r_{n})}\nabla_{\tau}v^{n}
-\frac{r_{n}^{2}}{\tau_{n}(1+r_{n})}\nabla_{\tau}v^{n-1},
\qquad n\ge2,
$$

第一层用 $D_2v^1:=D_1v^1=\nabla_\tau v^1/\tau_1$。$r_n=1$ 时它退化为经典的 $(3v^n-4v^{n-1}+v^{n-2})/(2\tau)$。写成离散卷积 $D_{2}v^{n}=\sum_{k=1}^{n}b^{(n)}_{n-k}\nabla_{\tau}v^{k}$ 时，

$$
b^{(1)}_0:=\frac{1}{\tau_1};
\qquad
b^{(n)}_0:=\frac{1+2r_n}{\tau_n(1+r_n)},
\qquad
b^{(n)}_1:=-\frac{r_n^2}{\tau_n(1+r_n)},
\qquad
b^{(n)}_j:=0\ (2\le j\le n).
$$

全离散格式（$\Lambda_h$ 为周期中心差分离散 Laplace 算子）是**完全隐式**的，既无稳定化也无凸分裂：

$$
D_{2}u^{n}=\varepsilon^{2}\Lambda_{h}u^{n}-f(u^{n}),
\quad n\ge1,
\qquad f(u^n):=(u^n)^{.3}-u^n .
$$

论文用两个命名条件（下面的措辞与原文一致）：

- **S1**：$0<r_{k}<\dfrac{3+\sqrt{17}}{2}\approx3.561$（$2\le k\le N$），用于**能量稳定性**。这个常数并非本文首创：它出自 Liao 与 Zhang（_Math. Comp._ 90 (2021) 1207-1226）的 Lemma 2.1，作为 BDF2 核正定性的条件，编号 48 与 52 都是把它作为 S1 引入的。编号 52 对它的评价很直接——称其为「由条件 S1 带来的人为常数」。
- **S0**：$0<r_{k}<1+\sqrt{2}\approx2.414$（$2\le k\le N$），用于**离散最大值原理**与**最大范数**收敛；论文明确指出它与 Grigorieff（1983）对常微分方程的零稳定性条件一致。

> [!warning] 两个常数对应两条不同性质
> **S1（$3.561$）给的是能量稳定性，S0（$1+\sqrt2$）给的是最大值原理与最大范数收敛。** 能量稳定性来自二次型的正定性，最大值原理来自核重排后的非负单调性，二者不能互换，也不是一个「更好」一个「更差」。引用时把它们混起来是这条文献里最常见的错误之一。

### 推导

**第一步：修正能量与 $3.561$ 的精确来源。** 耗散的对象不是普通离散能量，而是修正能量

$$
\widehat{E}[u^{k}]:=E[u^{k}]+\frac{r_{k+1}\tau_{k}}{2(1+r_{k+1})}
\sum_{i=1}^{M}\bigl(\partial_{\tau}u_{i}^{k}\bigr)^{2},
\quad k\ge1,
\qquad
\widehat{E}[u^{0}]:=E[u^{0}],
$$

$$
E[u^{k}]:=-\frac{\varepsilon^{2}}{2}(u^{k})^{T}\Lambda_{h}u^{k}
+\frac{1}{4}\sum_{i=1}^{M}\bigl(1-(u_{i}^{k})^{2}\bigr)^{2},
\quad k\ge0 .
$$

修正项是 $O(\tau)$，故 $\tau\to0$ 时 $\widehat E\to E$。注意它用的是 $r_{k+1}$，即**下一个**步长比；这个前瞻正是让电报式求和成立的原因。

用 $2a(a-b)=a^2-b^2+(a-b)^2$ 与核的定义可得逐点的不等式

$$
D_{2}u_{i}^{n}\,(\nabla_\tau u_{i}^{n})
\ \ge\
\frac{r_{n+1}\tau_{n}}{2(1+r_{n+1})}(\partial_{\tau}u_{i}^{n})^{2}
-\frac{r_{n}\tau_{n-1}}{2(1+r_{n})}(\partial_{\tau}u_{i}^{n-1})^{2}
+\Bigl(\frac{2+4r_{n}-r_{n}^{2}}{1+r_{n}}-\frac{r_{n+1}}{1+r_{n+1}}\Bigr)
\frac{\tau_{n}}{2}(\partial_{\tau}u_{i}^{n})^{2}.
$$

**这一行是全篇能量论证的全部内容**：前两项电报式抵消进 $\widehat E$，最后一项在吸收非线性项贡献的 $-\frac{\tau_n^2}{2}\sum_i(\partial_\tau u_i^n)^2$ 之后必须非负。把括号在最坏情形下的正性写出来，就是

$$
\frac{r_{k+1}}{1+r_{k+1}}<\frac{r_{s}}{1+r_{s}}=\frac{\sqrt{17}-1}{4}\approx0.78,
\qquad
r_{s}=\frac{3+\sqrt{17}}{2}\ \text{是}\ 2+3r-r^{2}=0\ \text{的正根}.
$$

数值上 $r_s=3.5616\ldots$，且 $r_s/(1+r_s)=3.5616/4.5616=0.7808=(\sqrt{17}-1)/4$，自洽。

论文另给出对应的步长上界讨论：令 $h(x):=\dfrac{2+4x-x^{2}}{1+x}$，则

$$
h'(x)=\frac{x+1+\sqrt3}{(1+x)^{2}}\bigl(\sqrt3-1-x\bigr),
$$

因此 $h$ 在 $(0,\sqrt3-1]$ 上递增、之后递减，$h(0)=2$、$h(\sqrt2+1)=1+\frac{\sqrt2}{2}$。据此分三段：

1. $0<r_k\le\sqrt3-1$：$\tau_k\le\min\{1,\frac{9-\sqrt{17}}{4}\}=1$ 即可；
2. $\sqrt3-1<r_k\le\sqrt2+1$：$\tau_k\le1+\frac{\sqrt2}{2}-\frac{\sqrt{17}-1}{4}\approx0.93$ 即可；
3. $\sqrt2+1<r_k<r_s$：必须控制**下一个**步长比，例如 $r_{k+1}\le\frac{2h(r_s)-1}{3-2h(r_s)}\approx0.39$ 时 $\tau_k\le\frac12$ 即可。

第三段值得留意：**步长比越接近 $r_s$，对下一步的约束越紧**，这正是修正能量里那个前瞻项在起作用。

> [!warning] 关于常被并列引用的 $4.8645$
> BDF2 核正定性的门槛后来被改进到 $r_*\approx4.864$（$1+2r-r^{3/2}=0$ 的正根，精确值 $4.864536512317583$），但这个改进**不属于本专题的任何一篇**。它出自 Liao、Ji、Wang 与 Zhang（_J. Sci. Comput._ 92 (2022) 52），那里印的是 $4.864$；流传更广的 $4.8645$ 写法出自 Zhang 与 Zhao（_J. of Math._ (PRC) 41(6) (2021) 471-488）。编号 69 与 91 都引用了 Liao-Ji-Wang-Zhang，这大概就是这个常数与这批工作产生关联的路径；完整的辨析见[[computational-mathematics/paper-notes/phase-field-and-time-stepping/index|本专题首页]]。此外要注意：$4.8645$ 改进的是 $3.561$（同一条正定性引理），而**不是**取代 $1+\sqrt2$——后者管的是另一件事。

**第二步：核重排与 $1+\sqrt2$ 的来源。** 这是本文的技术招牌，推导干净得出奇。引入参数 $\eta\in\mathbb R$ 与

$$
\bar{v}^{0}:=v^{0},
\qquad
\bar{v}^{k}:=v^{k}-\eta\,v^{k-1}\ (k\ge1),
\qquad
v^{k}=\sum_{\ell=0}^{k}\eta^{k-\ell}\bar{v}^{\ell} .
$$

代入并交换求和顺序，得到**重排后的** BDF2 公式

$$
D_{2}v^{n}\equiv\sum_{j=1}^{n}d^{(n)}_{n-j}\nabla_{\tau}\bar{v}^{j}+d^{(n)}_{n}\bar{v}^{0},
\qquad
d^{(n)}_{n-j}:=\sum_{k=j}^{n}b^{(n)}_{n-k}\eta^{k-j},
\qquad
d^{(n)}_{n}:=\eta\,d^{(n)}_{n-1},
$$

其闭式极其简单：

$$
d^{(n)}_{0}=b_{0}^{(n)},
\qquad
d^{(n)}_{j}=\eta^{j-1}\bigl(b_{0}^{(n)}\eta+b_{1}^{(n)}\bigr)\quad(1\le j\le n).
$$

目标是让新核**非负且单调递减**，$d^{(n)}_{0}\ge d^{(n)}_{1}\ge\cdots\ge d^{(n)}_{n}\ge0$。代入 BDF2 核，这恰好等价于

$$
\frac{r_{k}^{2}}{1+2r_{k}}\ \le\ \eta\ <\ 1\qquad(k\ge2).
$$

**这样的 $\eta$ 存在当且仅当 $\dfrac{r_k^{2}}{1+2r_k}<1$，即 $r_{k}^{2}-2r_{k}-1<0$，即 $r_{k}<1+\sqrt{2}$。** 这正是条件 S0，也正是 Grigorieff 的经典零稳定性界。

$\eta$ 的最优取法也给出了：令

$$
K(\eta):=\frac{1-\eta}{\eta^{2}}\cdot\frac{(1+2r_{n})\eta-r_{n}^{2}}{1+r_{n}},
\qquad
K'(\eta)=\frac{1+r_{n}}{\eta^{3}}\Bigl(\frac{2r_{n}^{2}}{(1+r_{n})^{2}}-\eta\Bigr),
$$

故 $K$ 在 $\eta\to\frac{2r_n^2}{(1+r_n)^2}$ 处取最大，论文取

$$
\eta:=\frac{2r_{s}^{2}}{(1+r_{s})^{2}},
\qquad r_{s}\in[1,\,1+\sqrt2)\ \text{为最大步长比} .
$$

两个算例：$r_s=1$ 给出 $\eta=\frac12$ 与 $\tau\le\frac{1}{2(S_n+4\varepsilon^{2}h^{-2})}$；$r_s=2$ 给出 $\eta=\frac89$ 与 $\tau_n\le\frac{1}{48}\cdot\frac{1}{S_n+4\varepsilon^{2}h^{-2}}$。**注意 $\eta$ 越接近 $1$，$1-\eta$ 越小，允许的步长就越小**——这就是「步长比越大、步长上界越紧」的代数原因。

**第三步：最大值原理的论证。** 在重排变量下格式变为

$$
\bigl(d^{(n)}_{0}+S_{n}-\varepsilon^{2}\Lambda_{h}\bigr)\bar{u}^{n}
=\sum_{j=0}^{n-1}Q_{n-j}^{(n)}\bar{u}^{j}+(S_{n}+1)u^{n}-(u^{n})^{.3},
$$

$$
Q_j^{(n)}:=\bigl(d^{(n)}_{j-1}-d^{(n)}_{j}-S_{n}\eta^{j}\bigr)I+\eta^{j}\varepsilon^{2}\Lambda_{h},
$$

其中 $S_n$ 是第二个自由参数。在步长限制

$$
\tau_{n}\le\frac{(1+2r_{n})\eta-r_{n}^{2}}{\eta^{2}(1+r_{n})}
\cdot\frac{1-\eta}{S_{n}+4\varepsilon^{2}h^{-2}}
$$

之下，$Q_j^{(n)}$ 的每个元素都非负，且 $\|Q_j^{(n)}\|_{\infty}\le d^{(n)}_{j-1}-d^{(n)}_{j}-S_{n}\eta^{j}$。**「每个元素都非负」正是重排的全部目的**：没有它，下面的归纳一步也走不动。

配上 $\ell^\infty$ 立方引理（若 $B=(b_{ij})$ 满足 $b_{ii}=-\max_i\sum_{j\ne i}|b_{ij}|$ 且 $A=aI-B$、$a>0$，则对 $c>0$ 有 $\|AV+cV^{3}\|_{\infty}\ge a\|V\|_{\infty}+c\|V\|_{\infty}^{3}$），完全归纳给出承载全篇的**辅助断言**

$$
\|\bar{u}^{k}\|_{\infty}\le 1-\eta\qquad(1\le k\le N),\ \text{只要}\ \|\bar u^0\|_\infty\le1,
$$

再由 $u^{k}=\sum_{\ell}\eta^{k-\ell}\bar u^{\ell}$ 与 $\sum_{j\ge0}\eta^j(1-\eta)=1$ 得到 $\|u^{k}\|_{\infty}\le1$。**这里的 $1-\eta$ 不是巧合**：它就是几何级数 $\sum\eta^j$ 的倒数，重排变量的界必须缩小这个因子，还原后才恰好是 $1$。论证中用到的两条初等事实是：$c\ge1$ 时 $g_c(z):=(c-1)z+z^{3}-c$ 在 $z>0$ 上递增；$|z|\le1$ 时 $|(c+1)z-z^{3}|\le c$。

**第四步：DCC 核与新的 Grönwall 不等式。** 建立在**重排后**的核 $d^{(n)}_j$ 之上：

$$
(Q_{d})_{0}^{(n)}:=\frac{1}{d_{0}^{(n)}},
\qquad
(Q_{d})_{n-j}^{(n)}:=\sum_{k=j+1}^{n}\frac{d_{k-j-1}^{(k)}-d_{k-j}^{(k)}}{d_{0}^{(j)}}(Q_{d})_{n-k}^{(n)},
\quad 1\le j\le n-1,
$$

$$
\sum_{j=k}^{n}(Q_{d})^{(n)}_{n-j}\,d_{j-k}^{(j)}\equiv 1\quad(1\le k\le n),
\qquad
0<(Q_{d})^{(n)}_{n-j}\le\frac{1}{d^{(j)}_{0}} .
$$

实用界为 $(Q_{d})^{(n)}_{n-1}\le1/b^{(1)}_{0}=\tau_1$ 与 $(Q_{d})^{(n)}_{n-j}\le1/b^{(j)}_{0}=\frac{1+r_j}{1+2r_j}\tau_j\le\tau_j$（$2\le j\le n$）。这类核归于 Liao-Li-Zhang（分数阶 Caputo）与 Liao-McLean-Zhang（一般核）。**注意 DCC 核的非负性同样依赖 $d^{(n)}_j$ 非负递减，所以它也是重排的下游产物。**

### 定理

**（唯一可解性）** 若 $\tau_{n}<\dfrac{1+2r_{n}}{1+r_{n}}$（$n\ge1$），格式唯一可解；论文注明取 $\tau_n<1$ 就够。证明先得 $b_0^{(n)}>1$，于是 $G_h:=b_0^{(n)}-1-\varepsilon^{2}\Lambda_h$ 正定，$u^n$ 是严格凸泛函 $\frac12w^TG_hw+\frac14\sum_kw_k^4-w^Tg(u^{n-1})$ 的唯一极小点，其中 $g(u^{n-1}):=b_0^{(n)}u^{n-1}-b_1^{(n)}\nabla_\tau u^{n-1}$。

**（能量稳定性，$3.561$ 的结论）** 设 **S1** 成立且

$$
\tau_{k}\le\min\Bigl\{\frac{1+2r_{k}}{1+r_{k}},\
\frac{2+4r_{k}-r_{k}^{2}}{1+r_{k}}-\frac{r_{k+1}}{1+r_{k+1}}\Bigr\},
\qquad k\ge1,
$$

则变步长 BDF2 解满足离散能量耗散律

$$
\widehat{E}[u^{k}]\le\widehat{E}[u^{k-1}],\qquad k\ge1 .
$$

**注意衰减的是修正能量 $\widehat E$ 而非 $E$。** 证明另用到代数恒等式

$$
4(a^{3}-a)(a-b)+2(1-a^{2})(a-b)^{2}
=(1-a^{2})^{2}-(1-b^{2})^{2}+(a^{2}-b^{2})^{2}.
$$

**（离散最大值原理，$1+\sqrt2$ 的结论）** 设 **S0** 成立，取 $\eta:=\frac{2r_s^2}{(1+r_s)^2}$ 且

$$
\tau_{n}\le\frac{(1+2r_{n})\eta-r_{n}^{2}}{\eta^{2}(1+r_{n})}
\cdot\frac{1-\eta}{2+4\varepsilon^{2}h^{-2}},\qquad n\ge1
$$

（即上面的步长条件取 $S_n=2$），则 $\|u^{0}\|_{\infty}\le1$ 蕴含 $\|u^{k}\|_{\infty}\le1$（$1\le k\le N$）。**这个步长限制是真实的，它把 $\tau$ 与 $h$、$\varepsilon$ 通过 $4\varepsilon^2h^{-2}$ 耦合起来**；论文自己称之为「一个实用的时间步长约束」。

**（借互补核建立的离散 Grönwall 不等式）** 设 $\kappa>0$、$\lambda\in(0,1)$，非负序列 $\{g^k\},\{w^k\}$ 满足

$$
\sum_{k=1}^{n}d^{(n)}_{n-k}\nabla_{\tau}w^{k}
\le\kappa\sum_{k=1}^{n}\lambda^{n-k}w^{k}+g^{n},
\qquad 1\le n\le N .
$$

若 $b^{(n)}_{0}\ge2\kappa$，则

$$
w^{n}\le2\exp\Bigl(\frac{2\kappa t_{n}}{1-\lambda}\Bigr)
\Bigl(w^{0}+\sum_{j=1}^{n}\frac{g^{j}}{b^{(j)}_{0}}\Bigr),
\qquad 1\le n\le N .
$$

证明就是「乘 $(Q_d)^{(n)}_{n-j}$、交换求和、用互补恒等式」这一招：

$$
\sum_{j=1}^n(Q_d)^{(n)}_{n-j}\sum_{k=1}^jd^{(j)}_{j-k}\nabla_\tau w^k
=\sum_{k=1}^n\nabla_\tau w^k\sum_{j=k}^n(Q_d)^{(n)}_{n-j}d^{(j)}_{j-k}
=w^n-w^0 .
$$

**（最大范数二阶收敛）** 设 $u_0$ 光滑、$|u_0|\le1$，精确解充分光滑。在 **S0** 与最大值原理定理同样的步长限制下，

$$
\|u(\mathbf x_{h},t_{n})-u_{h}^{n}\|_{\infty}
\le\frac{C_{u}t_{n}}{1-\eta}\exp\Bigl(\frac{4t_{n}}{1-\eta}\Bigr)\bigl(\tau^{2}+h^{2}\bigr),
\qquad 1\le n\le N,
$$

其中 **$C_u$ 与步长大小、步长比都无关**。非线性项不需要任何 Lipschitz 假设：由离散最大值原理知解落在 $[-1,1]$ 内，再用 $|(a^{3}-a)-(b^{3}-b)|\le2|a-b|$（$a,b\in[-1,1]$）得到 $\|f(U^n)-f(u^n)\|_\infty\le2\|e^n\|_\infty$，于是 $\kappa=2$，与前面取 $S_n=2$ 相合。论文声明这是在 Grigorieff 零稳定性条件下非均匀 BDF2 的**第一个**此类收敛结果。

### 数值实验

三个算例。

| 算例 | 设置                                                                                                                                                                                                                      | 检验对象                       |
| ---- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------ |
| 1    | 带外力 Allen-Cahn：$\partial_tu=\frac{1}{8\pi^2}\Delta u-f(u)+g$，$(0,1)^2\times(0,1)$，制造解 $u=\sin(2\pi x)\sin(2\pi y)\sin t$，**随机网格** $\tau_k=T\epsilon_k/S$（$\epsilon_k\in(0,1)$ 随机，$S=\sum_k\epsilon_k$） | 时间二阶精度与步长比稳健性     |
| 2    | 四气泡合并：$\varepsilon=0.02$，$\Omega=(-1,1)^2$，每方向 $128$ 个网格，初值为中心在 $(\pm0.3,0)$ 与 $(0,\pm0.3)$、半径 $0.2$ 的四个 $\tanh$ 剖面之积                                                                     | 界面演化的定性正确性           |
| 3    | 粗化动力学：$\varepsilon=0.01$，$\Omega=(0,1)^2$，$128\times128$ 均匀空间网格，随机初值                                                                                                                                   | 多个 $\tau$ 下的最大范数与能量 |

自适应策略为 $\tau_{\mathrm{ada}}(e,\tau)=\rho\bigl(\frac{tol}{e}\bigr)^{1/2}\tau_{\mathrm{cur}}$，默认参数 $\rho=0.6$、$tol=10^{-4}$、$\tau_{\max}=0.1$、$\tau_{\min}=10^{-3}$；一阶/二阶配对取后向 Euler 与自适应 BDF2。

**算例一是三者中最有信息量的。** 它用的不是精心设计的网格，而是**纯随机**步长序列，因此步长比是不受控的；在这样的网格上仍观测到时间二阶精度，说明分析对步长比的稳健性不只是纸面上的。算例三跟踪最大范数与能量，二者分别验证最大值原理与能量耗散。

粗化算例的初值在原文中描述得自相矛盾：一处写作 $u_{0}=0.95+\mathrm{rand}(\mathbf x)\times0.05$，另一处又说随机数「从 $-0.05$ 变化到 $0.05$」。两者不能同时成立——前者给出均值 $0.95$ 的初值，后者暗示以 $0$ 为中心的扰动。

理论与实验之间的落差在于：定理给的是 $\tau^2+h^2$ 的**上界**，而算例一在随机网格上（步长比可以很大）依然二阶，说明 **S0/S1 是充分条件而非必要条件**。这一点在编号 52 的实验里被量化得更彻底。

### 与其他论文的关系

这是整个系列的**枢纽**。它把互补核的想法从分数阶设定（[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 43]] 以及 Liao-Li-Zhang、Liao-McLean-Zhang）移植到整数阶变步长 BDF2，办法就是那个把符号不定的 $\{b^{(n)}_0>0,b^{(n)}_1<0\}$ 换成非负递减 $\{d^{(n)}_j\}$ 的重排。编号 52 把同一个 $3.561$ 能量门槛搬到无斜率选择的分子束外延模型，并把机器从 KRC 升级到 **DOC 核**，在同一门槛下额外拿到 $L^2$ 稳定性与误差估计；编号 58 把 DOC 推广到 BDF-$k$（$3\le k\le5$）；编号 67 做变步长 BDF3 并得到更紧的 $1.4877$；编号 74 则把「二次型是否正定」这个问题一般地解决。本文结论部分提出的两个问题后来都有了回应：S1/S0 是否最优（部分地由 Liao-Ji-Wang-Zhang 的 $4.864$ 回答），以及如何为**时间分数阶**相场方程构造非均匀 BDF2 型格式（由[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 57]] 回答）。论文还与 Chen-Wang-Yan-Zhang（2019）作了比较：后者用广义离散 Grönwall 不等式分析了带凸分裂的非均匀 BDF2 用于 Cahn-Hilliard，所需条件为 $r_k<1.53$。

## 52：同一套分析用于分子束外延模型

### 直觉

无斜率选择的分子束外延模型比 Allen-Cahn 难在两处。第一，自由能里有一个**对数项**，它上有界（被 $0$ 界住）但**下无界且没有相对极小**——所以没有能量上偏好的 $\nabla\Phi$ 值，也就没有「斜率选择」，物理上表现为真正的粗-平-粗多尺度形貌。第二，方程是**四阶**的。

这两点各自封死了一条路。四阶加上高度函数的性质意味着**这里没有最大值原理可用**：$\Phi$ 是薄膜高度，不是取值受限的序参量，论文在连续层面建立的是体积守恒、能量耗散与一条 $L^2$ 增长界，其中没有一条能给出 $\ell^\infty$ 的不变区域。于是编号 48 那套 KRC$+\ell^\infty$ 的机器在这里**无处着力**，能量与 $L^2$ 成为仅有的工具。

本文的一步是换到卷积的另一侧。KRC 是把核换掉；**DOC 是把算子求逆**。构造 $\theta$ 使得下三角矩阵 $\Theta_2=B_2^{-1}$，于是把格式在第 $j$ 层乘上 $\theta^{(n)}_{n-j}$ 再对 $j$ 求和，左端的多步算子就**塌成一个一阶差分** $\nabla_\tau u^n$，右端变成空间项与非线性项的一个卷积。**多步稳定性问题就此化归为一个二次型的正定性问题**，而后者恰好就是编号 48 里给出 $3.561$ 的那个问题。这就是为什么两篇拿到同一个门槛，尽管模型完全不同。

### 问题设定

$\Omega\subset\mathbb R^2$ 有界、周期边界，

$$
\Phi_t=-\varepsilon\Delta^{2}\Phi-\nabla\cdot\mathbf f(\nabla\Phi),
\qquad
\mathbf f(\mathbf v):=\frac{\mathbf v}{1+|\mathbf v|^{2}},
\qquad \Phi(\mathbf x,0)=\Phi_0(\mathbf x),
$$

它是自由能

$$
E[\Phi]=\int_{\Omega}\Bigl[\frac{\varepsilon}{2}(\Delta\Phi)^{2}
-\frac{1}{2}\ln\bigl(1+|\nabla\Phi|^{2}\bigr)\Bigr]\mathrm d\mathbf x
$$

的 $L^2(\Omega)$ 梯度流。这里 $\Phi$ 是共动坐标系中薄膜的无量纲高度，非线性二阶项刻画 Ehrlich-Schwoebel 效应，线性四阶项刻画表面扩散，$\varepsilon>0$ 表示原本刻面化的晶体薄膜上圆角的宽度。连续性质有三条：体积守恒 $(\Phi(t),1)=(\Phi_0,1)$、能量耗散

$$
\frac{\mathrm d}{\mathrm dt}E[\Phi]=-\|\Phi_t\|_{L^2(\Omega)}^{2}\le0,
$$

以及 $L^2$ 估计 $\|\Phi\|_{L^{2}(\Omega)}\le e^{t/(4\varepsilon)}\|\Phi_0\|_{L^{2}(\Omega)}$。适定性归于 Li 与 Liu（2003）：$\Phi_0\in H^m_{per}(\Omega)$（$m\ge2$）时存在唯一弱解，$\Phi\in L^{\infty}(0,T;H^m)\cap L^2(0,T;H^{m+2})$、$\partial_t\Phi\in L^2(0,T;H^{m-2})$。

时间离散与编号 48 完全相同（同样的 $b^{(n)}_j$，同样以 $r_1\equiv0$ 把第一层嵌成 BDF1），空间用中心差分，格式完全隐式：

$$
D_{2}\phi_{h}^{n}+\varepsilon\Delta_{h}^{2}\phi_{h}^{n}
+\nabla_{h}\cdot\mathbf f(\nabla_{h}\phi_{h}^{n})=0,
\qquad 1\le n\le N,\quad \phi_h^0=\Phi_0(\mathbf x_h).
$$

> [!note] 与编号 78 处理的是**不同**的分子束外延模型
> 本文是**无斜率选择**的模型，自由能里是对数项 $-\frac12\ln(1+|\mathbf v|^2)$；[[computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation|编号 78]] 处理的是**带斜率选择**的模型，自由能里是双阱 $\frac14(|\mathbf v|^2-1)^2$。**无斜率选择的情形更难**，因为非线性项不是多项式有界的，自由能也没有下界。本文结论部分明说，这里的技术「不适用于带斜率选择的分子束外延模型」，把它留作未来工作。两者不可混为一谈。

### 推导

**第一步：DOC 核的正式定义。** 本文是把它明确写出来的那一篇（并把线性扩散情形归于 Liao-Zhang）：

$$
\theta_{0}^{(n)}:=\frac{1}{b_{0}^{(n)}},
\qquad
\theta_{n-k}^{(n)}:=-\frac{1}{b_{0}^{(k)}}\sum_{j=k+1}^{n}\theta_{n-j}^{(n)}b_{j-k}^{(j)}
\qquad(1\le k\le n-1),
$$

**正交恒等式**为

$$
\sum_{j=k}^{n}\theta_{n-j}^{(n)}\,b_{j-k}^{(j)}\;\equiv\;\delta_{nk}
\qquad(1\le k\le n),
$$

$\delta_{nk}$ 是 Kronecker 记号。真正驱动一切的推论是

$$
\sum_{j=1}^{n}\theta_{n-j}^{(n)}\,D_{2}v^{j}\;=\;\nabla_{\tau}v^{n}
\qquad(n\ge1).
$$

也就是说，**把格式与 DOC 核作卷积，就把多步算子剥掉，只剩一个一阶差分**，于是标准的一步能量论证重新可用。与编号 43、48 的 DCC 核对照：后者的恒等式是 $\sum_{j=k}^{n}P^{(n)}_{n-j}A^{(j)}_{j-k}\equiv1$（互补/求和型，给出 $\sum_jP\sum_kA\nabla w^k=w^n-w^0$）。**DOC 是正交（Kronecker $\delta$），DCC 是互补（常数 $1$）。**

矩阵形式把这件事说得最清楚：设 $\mathbf B_2$ 是 BDF2 核构成的下双对角矩阵（$\otimes\,\mathbf I_2$），$\mathbf\Theta_2$ 是 DOC 核构成的下三角矩阵，则正交恒等式就是

$$
\mathbf\Theta_{2}=\mathbf B_{2}^{-1},
$$

而正定性引理说的是 $\mathbf B:=\mathbf B_{2}+\mathbf B_{2}^{T}$ 正定。**DOC 核字面上就是 BDF2 矩阵之逆的元素**——这是解释「DOC 核是什么」最干净的说法。

**第二步：DOC 核的显式公式与三条性质。** 在 **S1** 下，

1. DOC 核 $\theta^{(n)}_{n-j}$ 正定；
2. 它们为正，且有闭式

$$
\theta_{n-j}^{(n)}=\frac{1}{b^{(j)}_{0}}\prod_{i=j+1}^{n}\frac{r_i^{2}}{1+2r_i},
\qquad 1\le j\le n;
$$

3. $\displaystyle\sum_{j=1}^{n}\theta_{n-j}^{(n)}=\tau_n$，因而 $\displaystyle\sum_{k=1}^{n}\sum_{j=1}^{k}\theta_{k-j}^{(k)}=t_n$（$n\ge1$）。

性质 2 把与 $1+\sqrt2$ 的关系显示得一清二楚：**乘积 $\prod_i\frac{r_i^2}{1+2r_i}$ 衰减当且仅当每个因子小于 $1$，即 $r_i<1+\sqrt2$**。性质 3 则是把空间相容误差转成 $t_nh^2$ 的那一步。

**第三步：二次型正定性与 $3.561$。** 在 **S1** 下，对任意非零序列 $\{w_k\}_{k=1}^n$，

$$
2w_{k}\sum_{j=1}^{k}b_{k-j}^{(k)}w_{j}
\ \ge\ \frac{r_{k+1}}{1+r_{k+1}}\frac{w_{k}^{2}}{\tau_{k}}
-\frac{r_{k}}{1+r_{k}}\frac{w_{k-1}^{2}}{\tau_{k-1}}
+\Bigl(\frac{2+4r_{k}-r_{k}^{2}}{1+r_{k}}-\frac{r_{k+1}}{1+r_{k+1}}\Bigr)\frac{w_{k}^{2}}{\tau_{k}},
\quad k\ge2,
$$

求和后

$$
\sum_{k=1}^{n}w_{k}\sum_{j=1}^{k}b_{k-j}^{(k)}w_{j}
\ \ge\ \frac{1}{2}\sum_{k=1}^{n}\Bigl(\frac{2+4r_{k}-r_{k}^{2}}{1+r_{k}}
-\frac{r_{k+1}}{1+r_{k+1}}\Bigr)\frac{w_{k}^{2}}{\tau_{k}}\ >\ 0,
\qquad n\ge2 .
$$

括号对一切容许步长比为正，恰好就是 **S1**：$0<r_{k}<r_{s}=\frac{3+\sqrt{17}}{2}\approx3.561$。**这与编号 48 里那条不等式是同一个代数事实的两种写法**（那里以 $\partial_\tau u$ 为变量，这里以一般序列 $w$ 为变量），这正是两篇拿到同一个常数的原因。论文注明它「与线性情形的结果一致」（Liao-Zhang，扩散方程），并说「就目前而言这似乎是文献中非线性问题的最好结果」。

**第四步：修正能量。**

$$
\mathcal{E}[\phi^{n}]:=E[\phi^{n}]
+\frac{r_{n+1}}{2(1+r_{n+1})\tau_{n}}\|\nabla_{\tau}\phi^{n}\|^{2},
\qquad 0\le n\le N,
\qquad \mathcal{E}[\phi^{0}]=E[\phi^{0}],
$$

$$
E[\phi^{n}]:=\frac{\varepsilon}{2}\|\Delta_{h}\phi^{n}\|^{2}
-\frac{1}{2}\bigl\langle\ln(1+|\nabla_{h}\phi^{n}|^{2}),1\bigr\rangle .
$$

修正项 $\frac{r_{n+1}}{2(1+r_{n+1})\tau_n}\|\nabla_\tau\phi^n\|^2$ 与编号 48 的 $\frac{r_{n+1}\tau_n}{2(1+r_{n+1})}\|\partial_\tau\phi^n\|^2$ 是同一个量的两种写法。

**第五步：与网格无关的常数 $\mathcal M_r$。** 这是本文误差估计能称为「稳健」的关键。用步长重标度后的矩阵定义

$$
\mathcal{M}_r:=\max_{n\ge1}\|\widetilde{\mathbf B}_{2}\|^{2}\,\|\mathbf L^{-1}\|^{4}
=\max_{n\ge1}\frac{\lambda_{\max}(\widetilde{\mathbf B}_{2}^{T}\widetilde{\mathbf B}_{2})}
{\lambda_{\min}^{2}(\widetilde{\mathbf B})} .
$$

在 **S1** 下有粗略估计 $\mathcal M_r<39$；而若实际计算不**持续**使用接近稳定性极限 $r_s=3.561$ 的大步长比，则 $\mathcal M_r\le4$。这两个数值得原样记住。

**第六步：处理对数非线性的两条不等式。** 凸性用

$$
g(\lambda):=\tfrac12\ln(1+|\mathbf u+\lambda\mathbf v|^2),
\qquad
g''(0)=\frac{1-|\mathbf u|^{2}}{(1+|\mathbf u|^{2})^{2}}\mathbf v^T\mathbf v\le\mathbf v^T\mathbf v,
$$

能量论证用向量不等式

$$
\frac{2(\mathbf u-\mathbf v)^{T}\mathbf u}{1+|\mathbf u|^{2}}
\le\ln\frac{1+|\mathbf u|^{2}}{1+|\mathbf v|^{2}}+|\mathbf u-\mathbf v|^{2} .
$$

后者是对数势的「链式法则替代品」，作用与[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 57]] 里那个特制的 $H$ 完全一致。

**第七步：两条 DOC 二次型不等式（本文的技术新意）。** (a) DOC 型 Cauchy-Young 不等式：在 **S1** 下，对任意 $\epsilon\ge0$，

$$
\sum_{k=1}^n\sum_{j=1}^k\theta^{(k)}_{k-j}(\mathbf z^k)^T\mathbf w^j
\le\frac{\epsilon}{2}\mathbf z^T\mathbf\Theta\mathbf z
+\frac{1}{2\epsilon}\mathbf w^T\mathbf B^{-1}\mathbf w ;
$$

(b) 针对非线性项的 DOC 不等式：

$$
\sum_{k,j}\theta^{(k)}_{k-j}(\mathbf z^k)^T
\bigl[\mathbf f(\mathbf v^j+\mathbf z^j)-\mathbf f(\mathbf v^j)\bigr]
\le2\sqrt{\mathcal M_r}\sum_{k,j}\theta^{(k)}_{k-j}(\mathbf z^k)^T\mathbf z^j,
$$

其中 **$\mathcal M_r$ 与 $t_n$、步长 $\tau_n$、步长比 $r_n$ 都无关**。这个常数的网格无关性正是全篇的要点所在。

### 定理

**（唯一可解性／凸性）** 若 $\tau_{n}\le4\varepsilon$，BDF2 格式是凸的，因而唯一可解。关键一步是 $\tau_n\le4\varepsilon$ 给出 $b_0^{(n)}>\frac{1}{4\varepsilon}$，足以压住上面 $g''(0)\le\mathbf v^T\mathbf v$ 带来的负贡献。

**（离散能量耗散律，分子束外延的 $3.561$ 结论）** 设 **S1** 成立且

$$
\tau_{n}\le4\varepsilon\min\Bigl\{1,\
\frac{2+4r_{n}-r_{n}^{2}}{1+r_{n}}-\frac{r_{n+1}}{1+r_{n+1}}\Bigr\},
\qquad n\ge1,
$$

则

$$
\mathcal{E}[\phi^{n}]\le\mathcal{E}[\phi^{n-1}]\le\mathcal{E}[\phi^{0}]=E[\phi^{0}],
\qquad n\ge1 .
$$

**（$L^2$ 稳定性）** 设 **S1** 成立且 $\tau_{n}\le\varepsilon/(16\mathcal M_r^{2})$。对扰动初值 $\bar\phi_h^0$，

$$
\|\bar{\phi}^{n}-\phi^{n}\|\le2\exp\bigl(16\mathcal M_r^{2}\,t_{n-1}/\varepsilon\bigr)\,
\|\bar{\phi}^{0}-\phi^{0}\|,
\qquad 1\le n\le N .
$$

证明就是那套 DOC 招式：把扰动方程乘 $\theta^{(k)}_{k-j}$ 并对 $j=1,\dots,k$ 求和，用 $\sum_j\theta^{(k)}_{k-j}D_2z^j=\nabla_\tau z^k$，再与 $2z^k$ 作内积并对 $k$ 求和。

**（$L^2$ 收敛性）** 设 $\Phi\in C^{(6,3)}_{\mathbf x,t}(\Omega\times(0,T])$，**S1** 成立且 $\tau_{n}\le\varepsilon/(16\mathcal M_r^{2})$，则

$$
\|\Phi^{n}-\phi^{n}\|\le C_{\phi}\exp\bigl(16\mathcal M_r^{2}t_{n-1}/\varepsilon\bigr)
\Bigl[\tau_{1}^{2}\sum_{k=1}^{n}\prod_{i=2}^{k}\frac{r_{i}^{2}}{1+2r_{i}}
+t_{n}(\tau^{2}+h^{2})\Bigr],
\qquad 1\le n\le N .
$$

**第一项值得单独看**：起步步长的误差 $\tau_1^2$ 被乘积 $\prod_i\frac{r_i^2}{1+2r_i}$ 阻尼，而这个乘积正是上面 DOC 核性质 2 里的那个因子。也就是说，**起步误差随时间被 DOC 核自动衰减掉**，这是把 DOC 核的显式闭式算出来才能看到的结论。

### 数值实验

**算例一（随机时间网格）。** 带外力的分子束外延模型，$\varepsilon=0.1$，$\Omega=(0,2\pi)^2$，制造解 $\Phi(x,t)=\cos(t)\sin(x)\sin(y)$，$T=1$，空间取 $3000$ 个网格点，时间步长 $\tau_k:=T\sigma_k/S$，$\sigma_k\in(0,1)$ 均匀随机。论文的表里除了 $L^2$ 误差与阶，还专门记录 $\max r_k$ 与 $N_1$——**步长比达到或超过 $(3+\sqrt{17})/2$ 的时间层数**。逐次加密时记录的数值为：

| 记录量                        | 逐次加密时的取值                                |
| ----------------------------- | ----------------------------------------------- |
| $\max_k r_k$                  | $2.94,\ 11.98,\ 34.82,\ 37.72,\ 71.89,\ 850.80$ |
| $N_1$（$r_k\ge3.561$ 的层数） | $0,\ 3,\ 7,\ 13,\ 24,\ 49$                      |
| 观测时间阶（相邻两行之间）    | $1.84,\ 2.29,\ 2.35,\ 2.42,\ 2.00$              |

**这张表是整条线索里最有说服力的一处理论-实验落差记录。** 结论很直接：即便有几十个时间层的步长比**大幅**违反 $3.561$（最大处 $850.80$，是门槛的两百多倍），二阶收敛依然保持。也就是说 **S1 是充分条件，远非必要条件**；变步长 BDF2 在经验上比理论保证的稳健得多。论文自己也把 $r_s=3.561$ 称为「由条件 S1 带来的人为常数」，这张表就是这句自评的证据。

**算例二（自适应步长）** 用于加速趋向定态解的计算。

**论文自己划定的边界**：结论部分明确指出，本文的技术「不适用于带斜率选择的分子束外延模型」，留作未来工作。带斜率选择的那一支由[[computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation|编号 78]] 用完全不同的方法（线性松弛加正则化能量重构，无变步长分析）处理。

### 与其他论文的关系

它是编号 48（KRC/DCC、$\ell^\infty$、Allen-Cahn）与一般 DOC 理论之间的桥。它从 Liao-Zhang 关于线性扩散的分析中引入 DOC 核，并证明这套工具经得起一个真正非线性、四阶、对数能量的模型，在**同一个** $3.561$ 门槛下同时给出能量稳定性与 $L^2$ 稳定性／收敛性。编号 58 随后把 DOC 推广到 BDF-$k$（$3\le k\le5$）；编号 67 用同样的哲学做变步长 BDF3 并得到更小的 $1.4877$；编号 69、[[computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation|91 与 104]] 则分别把 DOC 用到时间滤波后向 Euler、隐显 Runge-Kutta 与隐显多步法上。后续的 Liao-Ji-Wang-Zhang（_J. Sci. Comput._ 92 (2022) 52）明确把本文列为三个 $3.561$ 的 DOC 结果之一，并通过改变加权指数把 Cahn-Hilliard 模型的门槛放宽到 $\approx4.864$。

## 58：把不 A-稳定的 BDF-$k$ 拉回教科书式论证

### 直觉

$A$-稳定的 BDF1 与 BDF2 允许直接的教科书式离散能量证明：在 $L^2$ 中用 $u^n$ 作检验函数就够了。对 $3\le k\le5$，BDF-$k$ **不** $A$-稳定，这条直接论证失效。自 Lubich、Mansour 与 Venkataraman（_IMA J. Numer. Anal._ 33 (2013) 1365-1385）以来的标准对策是 **Nevanlinna-Odeh 乘子技术**（_Numer. Funct. Anal. Optim._ 3 (1981) 377-423），它基于 Dahlquist 关于 $A$-稳定与 $G$-稳定等价的结论，改用 $u^n-\sum_i\eta_iu^{n-i}$ 而不是 $u^n$ 作检验；另一条路是 Liu 的伸缩公式。两者都引入人为乘子，而且关键地**对启动数据要求更强的范数**——乘子型稳定性估计里会出现 $H^1$ 型的量。

本文提出并回答的问题是：对 $3\le k\le5$ 的 BDF-$k$，能不能恢复那个直接的离散能量分析？

**它的办法不是换检验函数，而是换方程。** 既然 $\Theta=B^{-1}$，就把整个格式左乘 $\Theta$：多步算子被剥掉，左端回到一个单纯的一阶差分，右端变成一个卷积。因为 $B\Theta=\Theta B=I$ 两侧都成立，这个变换是**可逆**的，没有信息损失。变换之后再用 $u^n$ 作检验，教科书论证就照常走完，而全部困难被集中到一个问题上：**DOC 核构成的二次型是否正定？**

而这个问题在均匀网格上有完整答案。$\{\theta_j\}$ 正定当且仅当 $\{b_j\}$ 正定，后者由对称化 Toeplitz 矩阵的生成函数 $\mathrm g^{(k)}(\varphi)=2\sum_jb_j^{(k)}\cos(j\varphi)$ 与 Grenander-Szegő 定理判定。**关键的事实是：$3\le k\le5$ 时 $\mathrm g^{(k)}_{\min}>0$，尽管 BDF-$k$ 并不 $A$-稳定。** 正是这条比 $A$-稳定弱的性质，替换了论证中原本需要 $A$-稳定的位置。

### 问题设定

有界凸区域 $\Omega$ 上的线性反应扩散方程，齐次 Dirichlet 边界：

$$
\partial_t u-\varepsilon\Delta u=\beta(x,t)\,u+f(t,x),
\qquad x\in\Omega,\ 0<t<T,
\qquad u|_{\partial\Omega}=0,\ u(0,x)=u_0(x),
$$

$\varepsilon>0$ 为常数，$|\beta(x,t)|\le\beta^{*}$。**注意本文的时间网格是均匀的**，$t_k=k\tau$，$\tau:=T/N$；$k\ge3$ 的变步长情形在结论部分被明确列为公开问题。

把 BDF-$k$ 写成离散卷积

$$
D_{k}v^n:=\frac1{\tau}\sum_{j=1}^n b_{n-j}^{(k)}\,\nabla_\tau v^{j},
\qquad n\ge k,
$$

其核对 $j\ge k$ 为零，由生成函数

$$
\sum_{\ell=1}^{k}\frac{1}{\ell}(1-\zeta)^{\ell-1}
=\sum_{\ell=0}^{k-1}b_{\ell}^{(k)}\zeta^{\ell},
\qquad 3\le k\le 5
$$

给出，显式为

| BDF-$k$ | $b_0^{(k)}$ | $b_1^{(k)}$ | $b_2^{(k)}$ | $b_3^{(k)}$ | $b_4^{(k)}$ |
| ------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| $k=2$   | $3/2$       | $-1/2$      |             |             |             |
| $k=3$   | $11/6$      | $-7/6$      | $1/3$       |             |             |
| $k=4$   | $25/12$     | $-23/12$    | $13/12$     | $-1/4$      |             |
| $k=5$   | $137/60$    | $-163/60$   | $137/60$    | $-21/20$    | $1/5$       |

**符号交替，负项还不止一个**——这就是为什么这里需要的不是「非负核」（那不可能），而是「二次型正定」。时间离散格式为（$u^1,\dots,u^{k-1}$ 假定已由启动过程给出）

$$
D_{k}u^{j}=\varepsilon\Delta u^{j}+\beta^j u^{j}+f^j,
\qquad k\le j\le N,
$$

弱形式为 $\langle D_k u^j,w\rangle+\varepsilon\langle\nabla u^j,\nabla w\rangle=\langle\beta^j u^j,w\rangle+\langle f^j,w\rangle$（$\forall w\in H^1_0(\Omega)$）。

### 推导

**第一步：DOC-$k$ 核与可逆变换。**

$$
\theta_{0}^{(k)}:=\frac{1}{b_{0}^{(k)}},
\qquad
\theta_{n-j}^{(k)}:=-\frac{1}{b_{0}^{(k)}}\sum_{\ell=j+1}^{n}\theta_{n-\ell}^{(k)}\,b_{\ell-j}^{(k)},
\qquad j=n-1,n-2,\dots,k,
$$

（约定 $i>j$ 时 $\sum_{k=i}^{j}\cdot=0$），满足正交恒等式

$$
\sum_{\ell=j}^{n}\theta_{n-\ell}^{(k)}\,b_{\ell-j}^{(k)}\equiv\delta_{nj}
\qquad(k\le j\le n),
$$

以及相互正交性 $\sum_{\ell=j}^{n}b_{n-\ell}^{(k)}\theta_{\ell-j}^{(k)}\equiv\delta_{nj}$（原文此处把 $\delta_{nj}$ 排成了 $\delta_{mk}$，是明显的排印疏漏）。于是

$$
\sum_{j=k}^{n}\theta_{n-j}^{(k)}\sum_{\ell=k}^{j}b_{j-\ell}^{(k)}\nabla_\tau u^{\ell}
=\sum_{\ell=k}^{n}\nabla_\tau u^{\ell}\sum_{j=\ell}^{n}\theta_{n-j}^{(k)}b_{j-\ell}^{(k)}
=\nabla_\tau u^{n},
$$

即**用 DOC 核作用把 BDF-$k$ 卷积逆回单个一阶差分**，只差一个启动值余项：

$$
\sum_{j=k}^{n}\theta_{n-j}^{(k)}D_{k}u^{j}
=\frac1{\tau}u_{\mathrm I}^{(k,n)}+\partial_{\tau}u^{n},
\qquad
u_{\mathrm I}^{(k,n)}:=\sum_{\ell=1}^{k-1}\nabla_\tau u^{\ell}
\sum_{j=k}^{n}\theta_{n-j}^{(k)}b_{j-\ell}^{(k)} .
$$

变换后的弱格式为：对一切 $w\in H^1_0(\Omega)$ 与 $k\le j\le N$，

$$
\langle\partial_\tau u^{j},w\rangle
+\varepsilon\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle\nabla u^{\ell},\nabla w\rangle
=-\frac1\tau\langle u_{\mathrm I}^{(k,j)},w\rangle
+\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle\beta^{\ell}u^{\ell},w\rangle
+\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle f^{\ell},w\rangle .
$$

取 $w=2\tau u^j$ 并对 $j=k,\dots,n$ 求和（丢掉 $\sum_j\|u^j-u^{j-1}\|^2$ 一项），得到**经典形式**的能量不等式

$$
\|u^n\|^2-\|u^{k-1}\|^2
\le -2\sum_{j=k}^{n}\langle u_{\mathrm I}^{(k,j)},u^j\rangle
-2\varepsilon\tau\sum_{j=k}^{n}\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle\nabla u^{\ell},\nabla u^{j}\rangle
+2\tau\sum_{j=k}^{n}\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle\beta^{\ell}u^{\ell},u^{j}\rangle
+2\tau\sum_{j=k}^{n}\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle f^{\ell},u^{j}\rangle .
$$

余下三个技术成分把论证收口：(i) $\theta_j^{(k)}$ 的正定性；(ii) $\theta_j^{(k)}$ 的衰减；(iii) 启动项 $u_{\mathrm I}^{(k,j)}$ 的衰减。

**第二步（i）：正定性由生成函数判定。** 引理 2.1 建立等价：$\{b_j^{(k)}\}$ 正（半）定当且仅当 $\{\theta_j^{(k)}\}$ 正（半）定。而 $\{b_j^{(k)}\}$ 的正定性由对称化 Toeplitz 矩阵 $B_k=B_{k,l}+B_{k,l}^{T}$ 的生成函数

$$
\mathrm g^{(k)}(\varphi)=2\sum_{j=0}^{k-1}b_j^{(k)}\cos(j\varphi)
$$

与 Grenander-Szegő 定理判定：$\mathrm g^{(k)}_{\min}\le\lambda_{\min}(B_k)\le\lambda_{\max}(B_k)\le\mathrm g^{(k)}_{\max}$。显式地，

$$
\mathrm g^{(3)}(\varphi)=\tfrac13\bigl(11-7\cos\varphi+2\cos2\varphi\bigr)
=\tfrac43\bigl(\cos\varphi-\tfrac78\bigr)^2+\tfrac{95}{48},
$$

$$
\mathrm g^{(4)}(\varphi)=\tfrac16\bigl(25-23\cos\varphi+13\cos2\varphi-3\cos3\varphi\bigr),
$$

$$
\mathrm g^{(5)}(\varphi)=\tfrac1{30}\bigl(137-163\cos\varphi+137\cos2\varphi-63\cos3\varphi+12\cos4\varphi\bigr).
$$

$k=3$ 的配方一眼可见下界为正。由此得到二次型下界（引理 2.4）：对 $3\le k\le5$ 与任意实序列 $\{w_j\}$，

$$
2\sum_{m=k}^{n}w_m\sum_{j=k}^{m}b_{m-j}^{(k)}w_j
\ \ge\ \sigma_k\sum_{j=k}^{n}w_j^2,
\qquad n\ge k,
$$

带**显式常数**

$$
\sigma_3=\frac{95}{48}\approx1.97919,
\qquad
\sigma_4=\frac{2656-43\sqrt{43}}{1458}\approx1.62828,
\qquad
\sigma_5\approx0.477683 .
$$

$k=4$ 时 $Z_4(x)=12-14x+26x^2-12x^3$ 在 $[-1,1]$ 上的极小点为 $x_{*}=(13-\sqrt{43})/18$；$k=5$ 时 $Z_5(x)=12+26x+178x^2-252x^3+96x^4$ 的极小点为

$$
x^{*}=\tfrac{1}{96}\Bigl(63-\sqrt[3]{49041-16\sqrt{3891895}}
-\tfrac{1121}{\sqrt[3]{49041-16\sqrt{3891895}}}\Bigr)\approx-0.064041,
\qquad Z_5(x^{*})\approx14.3305 .
$$

**$\sigma_k$ 随 $k$ 迅速变小（$1.979\to1.628\to0.478$）是有含义的**：它量化了「BDF-$k$ 离失去正定性还有多远」，$k=5$ 已经只剩不到 $k=3$ 的四分之一。

**第三步（ii）：DOC-$k$ 核的几何衰减。** 虽然 $b_j^{(k)}=0$（$j\ge k$），DOC 核 $\theta_j^{(k)}$ 却**永不消失**——求逆把一个有限支撑的核变成了一个无限长的核。好在它几何衰减：

$$
\bigl|\theta_j^{(k)}\bigr|\le\frac{\rho_k}{4}\Bigl(\frac{k}{7}\Bigr)^{j},
\qquad 3\le k\le5,\ j\ge0,
\qquad
\rho_3=\frac{10}{3},\quad \rho_4=6,\quad \rho_5=\frac{96}{5}.
$$

公比 $k/7$（即 $3/7$、$4/7$、$5/7$）通过 $\sum_{j\ge\ell}(k/7)^{j-\ell}\le\frac{7}{7-k}$ 产生了最终常数里那个无处不在的因子 $\frac{7}{7-k}$，对 $k=3,4,5$ 分别是 $\frac74$、$\frac73$、$\frac72$。

**第四步（iii）：启动项的衰减。** 存在有限的 $c_{\mathrm I,k}>1$ 使

$$
\bigl|u_{\mathrm I}^{(k,j)}\bigr|\le\frac{c_{\mathrm I,k}\rho_k}{8}
\Bigl(\frac{k}{7}\Bigr)^{j-k}\sum_{\ell=1}^{k-1}\bigl|\nabla_\tau u^{\ell}\bigr|,
\qquad
\sum_{j=k}^{n}\bigl|u_{\mathrm I}^{(k,j)}\bigr|
\le\frac{7c_{\mathrm I,k}\rho_k}{8(7-k)}\sum_{\ell=1}^{k-1}\bigl|\nabla_\tau u^{\ell}\bigr| .
$$

论文只显式给出 $c_{\mathrm I,3}=11/7$；$c_{\mathrm I,4}$ 与 $c_{\mathrm I,5}$ 只断言有限，做法是把 $n=4,5$ 与 $n=5,6,7$ 这有限几个低指标情形单独吸收进常数，**它们的具体数值没有给出**。

### 定理

**（耗散情形 $\beta=\beta(x)\le0$）** 对 $3\le k\le5$ 与 $n\ge k$，

$$
\|u^n\|\le\|u^{k-1}\|
+\frac{7c_{\mathrm I,k}\rho_k}{4(7-k)}\sum_{\ell=1}^{k-1}\|\nabla_\tau u^{\ell}\|
+\frac{7\rho_k}{2(7-k)}\sum_{\ell=k}^{n}\tau\|f^{\ell}\|
\ \le\ \frac{7\rho_k}{2(7-k)}\Bigl(c_{\mathrm I,k}\sum_{\ell=0}^{k-1}\|u^{\ell}\|
+\sum_{\ell=k}^{n}\tau\|f^{\ell}\|\Bigr).
$$

**这一条完全没有步长限制。**

**（一般有界 $\beta$）** 若 $|\beta(x,t)|\le\beta^{*}$ 且均匀步长满足

$$
\tau\le\frac{7-k}{7\rho_k\beta^{*}},
$$

则对 $k\le n\le N$，

$$
\|u^n\|\le\frac{7\rho_k}{7-k}\exp\Bigl(\frac{7\rho_k}{7-k}\beta^{*}t_{n-k}\Bigr)
\Bigl(c_{\mathrm I,k}\sum_{\ell=0}^{k-1}\|u^{\ell}\|+\sum_{\ell=k}^{n}\tau\|f^{\ell}\|\Bigr).
$$

论文只印出符号形式 $\tau\le(7-k)/(7\rho_k\beta^{*})$；代入 $\rho_k$ 的值，这个限制读作 $\tau\le6/(35\beta^{*})$（$k=3$）、$\tau\le1/(14\beta^{*})$（$k=4$）、$\tau\le2/(134.4\,\beta^{*})$（$k=5$）。

**（$L^2$ 收敛，$k$ 阶）** 在同样的步长限制与相容性界 $|\eta^j|=|D_k u(t_j)-\partial_t u(t_j)|\le C_u\tau^{k}\max_{t_k\le t\le T}|\partial_t^{(k+1)}u(t)|\le C_u\tau^{k}$（$j\ge k$）下，

$$
\|u(t_n)-u^{n}\|\le\frac{7\rho_k c_{\mathrm I,k}}{7-k}
\exp\Bigl(\frac{7\rho_k\beta^{*}t_{n-k}}{7-k}\Bigr)
\Bigl(\sum_{\ell=0}^{k-1}\|u(t_{\ell})-u^{\ell}\|+C_u\,t_{n-k}\,\tau^{k}\Bigr),
\qquad k\le n\le N .
$$

时间 $k$ 阶，且**启动值误差只以 $L^2$ 范数进入**——这正是相对乘子技术所宣称的正则性优势，后者需要启动数据的更强范数（论文引 Akrivis-Katsoprinakis，_Math. Comp._ 85 (2016) 2195-2216 的 Proposition 5.1 与 Theorem 5.1 作对照）。

### 数值实验

**本文没有数值实验。** 它是一篇纯分析论文，全文仅有的两组图形是：(a) $3\le k\le5$ 时生成函数 $\mathrm g^{(k)}(\varphi)$ 在 $[-\pi,\pi]$ 上的曲线，用来显示 $\mathrm g^{(k)}_{\min}>0$；(b) 计算出的 DOC-$k$ 核与界 $\frac{\rho_k}{4}(k/7)^{x}$ 的对照图，用来显示几何衰减。

这两组图与上面的两条引理一一对应，可以看作对**引理**而非对**格式**的数值验证。**本文因此没有关于收敛阶的实测数据**，$k$ 阶只是定理的结论。

### 与其他论文的关系

DOC 技术本身来自 Liao-Zhang（_Math. Comp._ 90 (2021) 1207-1226，那里用于**变步长 BDF2**）与编号 52 的配套分析；编号 48 与 52 把 DOC/DCC 用于非线性相场问题的变步长 BDF2，而本文是同一纲领的**均匀网格、高阶（$k=3,4,5$）**分支。它是编号 67（变步长 BDF3）的直接方法学祖先：编号 67 接手的正是本文结论部分第三条（「研究 BDF-$k$（$3\le k\le5$）变步长情形的离散能量技术」），并在那里得到门槛 $r_k<1.4877$。它的正定性机器（对卷积核构成的对称化 Toeplitz 型用 Grenander-Szegő）与编号 74 推广到**变步长 L1 型核**的那套是同一件事的两个版本。它所取代的那条路线——Akrivis 等人的 Nevanlinna-Odeh 乘子／$G$-稳定性——正是[[computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation|编号 104]] 用半生成函数论证重新审视的对象。

## 67：变步长 BDF3

### 直觉

变步长 BDF3 此前基本只有一条经典结果：Calvo 与 Grigorieff（2002）在步长比条件 $r_k<1.199$ 下证明 $L^2$ 稳定性，估计形如

$$
\|u^n\|\le C\exp(C\Gamma_n)\Bigl(\|u_0\|+\sum_{j=1}^{n}\tau_j\|f^j\|\Bigr),
\qquad \Gamma_n:=\sum_{k=2}^{n}|r_k-r_{k-1}| .
$$

**问题不在门槛，而在前因子 $\exp(C\Gamma_n)$ 不是网格稳健的。** 取交替步长 $\{\tau_1,\mu\tau_1,\tau_1,\mu\tau_1,\dots\}$（$\mu\ne1$）并固定 $T=\frac M2(1+\mu)\tau_1$，则

$$
\Gamma_M=(M-1)\bigl|\mu-\mu^{-1}\bigr|\longrightarrow\infty
\qquad(\tau_1\to0),
$$

界恰好在它本应覆盖的自适应区域内退化。**一个只在网格几乎均匀时才有意义的变步长定理，价值有限。**

本文的取舍很明确：用一个**常数与步长比完全无关**的分析取代它，代价是一个步长比门槛。做法有两层。能量那一层是构造一个显式的 Lyapunov 泛函 $G$，使得

$$
2v_n\tau_n\sum_{j}d^{(n)}_{n-j}v_j
=G[v_n,v_{n-1}]-G[v_{n-1},v_{n-2}]+F[v_n,v_{n-1},v_{n-2}],
\qquad F\ge\frac{\tau_n}{50}v_n^2 .
$$

这与编号 48、52 的「电报式相消加一个非负余项」是同一个模式，只是 BDF3 是三项算子，所以 $G$ 必须**同时携带两层**。$L^2$ 那一层则用 DOC 核，关键性质是 $\sum_j|\vartheta|\le K_3$ 且 $K_3$ 与步长比无关——这就是标题里「网格稳健」的来源。

### 问题设定

有界凸区域 $\Omega$ 上 $\partial_t u-\varepsilon\Delta u=\kappa(x)u+f(t,x)$，$u|_{\partial\Omega}=0$，$\varepsilon>0$ 常数，$|\kappa(x)|\le\kappa^{*}$。网格与步长比同前，$r_k:=\tau_k/\tau_{k-1}$（$2\le k\le N$）。

变步长 BDF3 公式（Calvo-Grigorieff 形式）为

$$
D_3v^n=d_0(r_n,r_{n-1})\partial_\tau v^{n}
+d_1(r_n,r_{n-1})\partial_\tau v^{n-1}
+d_2(r_n,r_{n-1})\partial_\tau v^{n-2},
$$

$$
d_0(x,y):=\frac{1+2x}{1+x}+\frac{xy}{1+y+xy},
$$

$$
d_1(x,y):=-\frac{x}{1+x}-\frac{xy}{1+y+xy}
-\frac{xy^2}{1+y+xy}\cdot\frac{1+x}{1+y},
\qquad
d_2(x,y):=\frac{xy^2}{1+y+xy}\cdot\frac{1+x}{1+y},
$$

对 $x,y\ge0$。写成卷积 $D_3v^n=\sum_{j=1}^{n}d^{(n)}_{n-j}\partial_\tau v^{j}$，其中 $d^{(n)}_j:=d_j(r_n,r_{n-1})$（$j=0,1,2$）、$d^{(n)}_j:=0$（$n\ge j+1\ge4$）。格式为 $D_3u^k=\varepsilon\Delta u^k+\kappa u^k+f^k$（$3\le k\le N$），$u^1,u^2$ 由启动过程给出。

### 推导

**第一步：变步长 DOC 核。**

$$
\vartheta_0^{(n)}:=\frac{1}{d_0^{(n)}},
\qquad
\vartheta_{n-j}^{(n)}:=-\frac{1}{d_0^{(j)}}\sum_{i=j+1}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-j},
\qquad 3\le j\le n-1,
$$

正交恒等式为 $\sum_{i=j}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-j}\equiv\delta_{nj}$（$3\le j\le n$），等价于矩阵形式 $\Theta_3D_3=I_{m\times m}$（$m=n-2$）；由 $D_3\Theta_3=I$ 又得相互正交性 $\sum_{i=j}^{n}d_{n-i}^{(n)}\vartheta^{(i)}_{i-j}\equiv\delta_{nj}$。

**注意变步长的 DOC 递归用的是 $1/d_0^{(j)}$（下标随运行指标 $j$ 变），与编号 58 均匀网格版本里固定的 $1/b_0^{(k)}$ 不同。** 作用 DOC 核后的等价卷积格式为

$$
\partial_\tau u^n=-\mathcal I_3^n[u]
+\sum_{k=3}^{n}\vartheta_{n-k}^{(n)}\bigl(\varepsilon\Delta u^k+\kappa u^k\bigr)
+\sum_{k=3}^{n}\vartheta_{n-k}^{(n)}f^k,
\qquad 3\le n\le N,
$$

$$
\mathcal I_3^n[v]:=\sum_{j=1}^{2}\partial_\tau v^{j}\sum_{i=3}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-j}
=\partial_\tau v^2\sum_{i=3}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-2}
+\vartheta_{n-3}^{(n)}d^{(3)}_2\,\partial_\tau v^1 .
$$

论文的分工是：**能量（梯度结构）论证用原始形式，$L^2$ 论证用 DOC 变换后的形式。**

**第二步：门槛 $R_e\approx1.4877$ 的来历。** 这个门槛不是一个天然常数，而是一次**参数取舍**的结果，值得写清楚。取 $\gamma=7/10$，则 $R_e$ 是

$$
d_1(R_e,0)+\tfrac{7}{10}\sqrt{R_e}\,d_2(R_e,R_e)=0
\qquad\Longleftrightarrow\qquad
\frac{10}{7(R_e+1)}-\frac{R_e^2\sqrt{R_e}}{R_e^2+R_e+1}=0
$$

的**唯一正根**，数值上 $R_e\approx1.4877$。

$\gamma=7/10$ 的来历是这样：离散梯度分解需要两个条件 $q_{n+1}\ge0$ 与 $p_{n+1}>0$ 同时成立，而它们含 $r_{n+1},r_n,r_{n-1},\gamma,R_e$ 五个变量，无法精确求解。论文的做法是退到等比网格上处理：由 $q_{n+1}\ge0$ 在 $r_{n-1}=0$、$r_{n+1}=r_n=r$ 时得到 $\gamma\le-d_1(r,0)/(\sqrt r\,d_2(r,r))$；把第二个条件也放到等比网格上，得 $\bar R_e\approx1.4965$ 与 $\bar\gamma\approx0.6924$。作者随后**固定 $\gamma=7/10$**（与 $0.6924$ 很接近）以换取可处理性，由此得到 $R_e\approx1.4877$。这个取法的依据是：$q_{n+1}\ge0$ 是「必要且尖锐」的，而 $p_{n+1}>0$ 可以放松。

**第三步：门槛的尖锐程度被量化了。** $\{\tau_nd^{(n)}_{n-k}\}$ 的正定性由五对角对称矩阵 $B_3=B_L+B_L^T$ 支配（$B_L$ 是元素为 $\tau_nd^{(n)}_j$ 的下三角矩阵）。对步长重标度后的 $\widetilde B_3=\Lambda_\tau^{-1}(B_L+B_L^T)\Lambda_\tau^{-1}$（$\Lambda_\tau=\mathrm{diag}(\sqrt{\tau_3},\dots,\sqrt{\tau_n})$），在 $r_k\sim U(0,R_e)$ 的随机网格上取 $200$ 次运行的最小值，得到最小特征值：

| $n$   | $R_e=1.20$ | $R_e=1.50$ | $R_e=1.69$  | $R_e=1.70$  |
| ----- | ---------- | ---------- | ----------- | ----------- |
| $50$  | $1.12$     | $5.08$e-01 | $6.12$e-02  | $-4.55$e-02 |
| $100$ | $1.07$     | $4.35$e-01 | $4.58$e-02  | $-5.29$e-02 |
| $200$ | $1.08$     | $4.18$e-01 | $-2.06$e-02 | $-8.49$e-02 |

**因此数值上 $R_e<1.69$ 是必要的，而理论给出 $R_e<1.4877$ 是充分的，两者之间的间隙很小。** 这与编号 52 对 $3.561$ 的自我评价（「由条件 S1 带来的人为常数」）恰好构成对照——同样是人为选定的常数，这里的距离是**被量化了**的。表里 $n=200$、$R_e=1.69$ 一格已经翻负，说明门槛的必要值本身也随 $n$ 缓慢下移。

**第四步：离散梯度结构（引理 3.3，本文的分析装置）。** 对 $0<x,y,z<R_e$ 定义

$$
d_*(x,y):=-\tfrac{10}{7}\sqrt{x}\,d_1(x,y)-\sqrt{xy}\,d_2(x,y),
$$

$$
p(x,y,z):=2d_0(y,z)-\sqrt{yz}\,d_2(y,z)-\tfrac{49}{100}d_*(y,z)-d_*(x,y),
\qquad
q(x,y,z):=d_*(y,z)-\sqrt{xy}\,d_2(x,y).
$$

若 $0<r_k<R_e$，则存在非负泛函 $G$ 与 $F$ 使

$$
2v_n\tau_n\sum_{j=3}^{n}d^{(n)}_{n-j}v_j
=G[v_n,v_{n-1}]-G[v_{n-1},v_{n-2}]+F[v_n,v_{n-1},v_{n-2}],
\qquad n\ge3,
$$

$$
G[v_n,v_{n-1}]:=d_*(r_{n+1},r_n)\,\tau_nv_n^2
+\sqrt{r_{n+1}r_n}\,d_2(r_{n+1},r_n)
\Bigl(\tfrac{7}{10}\sqrt{\tau_n}\,v_n-\sqrt{\tau_{n-1}}\,v_{n-1}\Bigr)^2,
$$

$$
F[v_n,v_{n-1},v_{n-2}]:=p(r_{n+1},r_n,r_{n-1})\tau_nv_n^2
+q(r_{n+1},r_n,r_{n-1})\Bigl(\tfrac{7}{10}\sqrt{\tau_n}v_n-\sqrt{\tau_{n-1}}v_{n-1}\Bigr)^2
$$

$$
\qquad\qquad
+\sqrt{r_nr_{n-1}}\,d_2(r_n,r_{n-1})
\Bigl(\sqrt{\tau_n}v_n-\tfrac{7}{10}\sqrt{\tau_{n-1}}v_{n-1}+\sqrt{\tau_{n-2}}v_{n-2}\Bigr)^2
\ \ge\ \frac{\tau_n}{50}v_n^2 .
$$

$G$ 与 $F$ 里那些 $7/10$ 就是上面的 $\gamma$，而常数 $1/50$ 是显式正定性常数的出处。由此（引理 3.4）：若 $0<r_k<R_e\approx1.4877$，则

$$
2\sum_{k=3}^{n}\xi_k\sum_{j=3}^{k}\tau_kd^{(k)}_{k-j}\xi_j
\ \ge\ \frac{1}{50}\sum_{k=3}^{n}\tau_k\xi_k^2,
\qquad n\ge3 .
$$

**第五步：DOC 核的正定性与一致可和性。** 在同一条件下（引理 4.1）$\{\tau_n\vartheta^{(n)}_{n-k}\}$ 正定。证明是标准的 DOC 招式：令 $\eta_k=\sum_{j=3}^{k}\vartheta^{(k)}_{k-j}\xi_j$，用相互正交性得 $\sum_kd^{(n)}_{n-k}\eta_k=\xi_n$，再对 $\{\eta_k\}$ 用引理 3.4。

另有一条（引自 Li-Liao 的引理 3.1）：若 $r_k\le R_e$，存在 $K_3>0$，**与 $t_n$ 及步长比 $r_k\in(0,R_e]$ 都无关**，使

$$
\sum_{j=3}^{n}\bigl|\vartheta^{(n)}_{n-j}\bigr|\le K_3,
\qquad
\sum_{j=i}^{n}\bigl|\vartheta^{(j)}_{j-i}\bigr|\le K_3 .
$$

**DOC 核不一定为正，但指数衰减**，这正是最终估计得以「网格稳健」的原因。

（本文对象是线性扩散方程，因此**没有**最大值原理的结论。）

### 定理

**（离散能量耗散律）** 设 $\kappa\le0$、$f\equiv0$。若 $0<r_k<R_e\approx1.4877$（$k\ge2$），则变步长 BDF3 格式**无条件能量稳定**：$E^n\le E^{n-1}$（$n\ge3$），其中**修正离散能量**为

$$
E^n:=\varepsilon\|\nabla u^n\|^2+\bigl(-\kappa u^n,u^n\bigr)
+\bigl(1,\,G[\partial_\tau u^n,\partial_\tau u^{n-1}]\bigr),
\qquad n\ge2 .
$$

注意这个能量是 $H^1$ 半范型的，并带一个非负的 $G$ 修正项（Lyapunov 型修正），与编号 48、52 的修正能量同构。

**（$L^2$ 稳定性，耗散情形 $\kappa<0$）** 若 $0<r_k<R_e$，

$$
\|u^n\|\le\|u^2\|+K_3\tau\|\partial_\tau u^1\|+4K_3\tau\|\partial_\tau u^2\|
+2\sum_{k=3}^{n}\tau_k\sum_{j=3}^{k}\bigl|\vartheta^{(k)}_{k-j}\bigr|\,\|f^j\|
\ \le\ \|u^2\|+K_3\tau\|\partial_\tau u^1\|+4K_3\tau\|\partial_\tau u^2\|
+2K_3t_n\max_{3\le k\le n}\|f^k\| .
$$

**（$L^2$ 稳定性，一般有界 $\kappa$）** 若 $|\kappa|\le\kappa^{*}$、$0<r_k<R_e$ 且 $\tau\le 1/(4K_3\kappa^{*})$，则对 $3\le n\le N$，

$$
\|u^n\|\le2\exp\bigl(4K_3\kappa^{*}t_{n-1}\bigr)
\Bigl(\|u^2\|+K_3\tau\|\partial_\tau u^1\|+4K_3\tau\|\partial_\tau u^2\|
+2K_3t_n\max_{3\le k\le n}\|f^k\|\Bigr).
$$

**（三阶 $L^2$ 收敛）** 在同样假设下，误差 $\tilde u^n=u(t_n)-u^n$ 满足

$$
\|\tilde u^n\|\le2\exp\bigl(4K_3\kappa^{*}t_{n-1}\bigr)
\Bigl(\|\tilde u^2\|+K_3\tau\|\partial_\tau\tilde u^1\|
+4K_3\tau\|\partial_\tau\tilde u^2\|+2K_3K_ut_n\tau^{3}\Bigr),
\qquad 3\le n\le N,
$$

其中 **$K_3$ 与 $K_u$ 与 $t_n$、步长 $\tau_n$、步长比 $r_n$ 都无关——即便 $r_n\to R_e$ 也如此**。这个无关性就是标题里的「网格稳健」，也正是 Calvo-Grigorieff 的前因子 $\exp(C\Gamma_n)$ 所缺的东西。

### 数值实验

热方程 $\partial_t u-\varepsilon\Delta u=f$，$\Omega=(0,2\pi)^2$，周期边界，$\varepsilon=0.1$，制造解 $u=\cos(t)\sin(x)\sin(y)$，$T=1$；误差取 $e(N)=\max_{1\le n\le N}\|v(t_n)-v^n\|$。两族网格：(a) 周期步长 $\{\tau_1,\mu\tau_1,\tau_1,\mu\tau_1,\dots\}$，$\tau_1=2/(N(1+\mu))$，$r_{\max}=\mu$；(b) 随机步长 $\tau_k=\epsilon_k/\sum\epsilon_k$，$\epsilon_k\sim U(0,1)$。启动值取两级三阶 SDIRK 或变步长 BDF2。

| 网格与启动方式               | $N=160,\dots,1280$ 上的观测阶 | 违反 $R_e$ 的层数 |
| ---------------------------- | ----------------------------- | ----------------- |
| $\mu=2R_e$，Runge-Kutta 启动 | $2.98,\ 2.99,\ 3.00,\ 3.00$   | 约 $N_1=N/2$      |
| $\mu=4R_e$，BDF2 启动        | $2.98,\ 2.99,\ 3.00,\ 2.99$   | 约 $N_1=N/2$      |

三条结论。第一，**大约一半的时间层步长比超过理论限制，格式依然稳定且三阶**——门槛充分但保守。第二，三阶 SDIRK 与二阶 BDF2 都足以作为启动过程达到三阶精度，与收敛定理的预言一致（定理里启动值只以 $\|\tilde u^2\|$、$\tau\|\partial_\tau\tilde u^1\|$、$\tau\|\partial_\tau\tilde u^2\|$ 的形式进入，因此二阶启动就够）。第三，在随机网格上方法同样网格稳健、同样三阶，即便许多步长比远超 $R_e$。

上面第三步里那张最小特征值表其实也是一个数值实验，只不过测的是引理而不是格式。两张表合起来把这个门槛的地位说清楚了：**$1.4877$ 充分，$1.69$ 必要，实际算例里超过它也往往没事。**

### 与其他论文的关系

论文自称是「我们变步长 BDF 离散能量分析系列工作之一」，前作是编号 48（_SIAM J. Numer. Anal._ 58:2294-2314）与 Liao-Zhang（_Math. Comp._ 90:1207-1226），两者都在摘要里按编号被引。它回答编号 58 结论部分的第三个公开问题（「研究 BDF-$k$（$3\le k\le5$）变步长情形的离散能量技术」）中 $k=3$ 的情形。DOC 装置与编号 48、52、58 在结构上完全相同；本文额外的装置是那个实现**离散梯度结构**的显式 Lyapunov 泛函 $G$，与编号 48、52 里修正能量所用的是同一件工具。

步长比的次序值得单独列清楚：

$$
1.199\ <\ 1.4877\ <\ 2.4142\ <\ 2.553\ <\ 3.5616\ <\ 4.8645 .
$$

其中 $1.199$ 是 Calvo-Grigorieff（2002）的经典 BDF3 $L^2$ 门槛（带不稳健前因子，被本文取代），$1.4877$ 是本文的 BDF3 能量与 $L^2$ 门槛，$2.4142=1+\sqrt2$ 是编号 48 的 S0，$2.553$ 是 Li-Liao 对常微分方程 BDF3 的 $R_3$（本文引用并使用），$3.5616$ 是 S1，$4.8645$ 属于 Liao-Ji-Wang-Zhang（不在本专题内）。**数值大小不代表强弱可比**：它们分别是不同格式、不同结论的门槛。

## 69：变步长时间滤波后向 Euler

### 直觉

时间滤波是数值天气预报里的一种廉价后处理：先算一步后向 Euler，再加上若干已算出层的一个线性组合，就把精度从一阶提到二阶，代价只是在遗留求解器里多写一行代码。Guzel 与 Layton（_BIT_ 58 (2018) 1-15）证明常步长下滤波后向 Euler（FiBE）是二阶且 $A$-稳定的，DeCaria、Guzel、Layton 与 Li（_SIAM J. Sci. Comput._ 43(3) (2021) A2130-A2160）在此基础上建了一个变步长变阶族。

变步长理论的障碍是结构性的，后续文献把它说得很明白：**滤波算法把滤波解与未滤波解耦合在一起，所以它是一个预估-校正格式而不是一步法，而且变步长 FiBE 不是 $A$-稳定的**，于是常步长情形所用的 $A$-稳定／$G$-稳定论证全部失效。

本文的桥是**单腿改写**。消去中间量之后，FiBE 等价于一个单腿多步法（OLM），而这个 OLM 的**左端恰好就是变步长 BDF2 差分算子**，右端把 $f$ 在 $u^{n},u^{n-1},u^{n-2}$ 的一个二阶精度组合处求值。左端一旦是 BDF2，编号 48、52、67 那套梯度结构机器就可以搬过来——摘要里「我们为时间滤波后向 Euler 格式**所对应的单腿多步格式**建立离散梯度结构」这句话说的就是这件事。

### 格式

变步长 FiBE 的一步出自 DeCaria、Guzel、Layton 与 Li。以 $\tau=k_{n+1}/k_n$ 记当前步与前一步之比：

$$
\text{后向 Euler 预估：}\qquad
\frac{y^{1}_{n+2}-y_{n+1}}{k_{n+1}}=f\bigl(t_{n+2},y^{1}_{n+2}\bigr),
$$

$$
\text{时间滤波：}\qquad
y_{n+2}=y^{1}_{n+2}-\frac{\tau(1+\tau)}{1+2\tau}
\Bigl(\frac{1}{1+\tau}\,y^{1}_{n+2}-y_{n+1}+\frac{\tau}{1+\tau}\,y_{n}\Bigr).
$$

改写到 Liao 的记号（$u^n$ 在 $t_n$，$\tau_n=t_n-t_{n-1}$，$r_n=\tau_n/\tau_{n-1}$，带波浪号者为未滤波的预估值）只是指标重标：

$$
\frac{\tilde u^{\,n}-u^{n-1}}{\tau_n}=f(t_n,\tilde u^{\,n}),
\qquad
u^{n}=\tilde u^{\,n}-\frac{r_n(1+r_n)}{1+2r_n}
\Bigl(\frac{1}{1+r_n}\tilde u^{\,n}-u^{n-1}+\frac{r_n}{1+r_n}u^{n-2}\Bigr).
$$

一致性检验：$r_n=1$ 时滤波退化为 Guzel-Layton 的经典形式

$$
u^{n}=\tilde u^{\,n}-\tfrac13\bigl(\tilde u^{\,n}-2u^{n-1}+u^{n-2}\bigr),
$$

直接计算即可验证。

滤波项其实是一个重标度的**二阶差商**：

$$
u^{n}=\tilde u^{\,n}-\eta^{(2)}\,\delta^{2}\tilde u,
\qquad
\eta^{(2)}=\frac{k_{n+1}}{\dfrac{1}{k_{n+1}}+\dfrac{1}{k_{n+1}+k_{n}}},
\qquad
\delta^{2}\tilde u=\frac{\dfrac{\tilde u^{\,n}-u^{n-1}}{k_{n+1}}
-\dfrac{u^{n-1}-u^{n-2}}{k_{n}}}{k_{n+1}+k_{n}} .
$$

消去中间值 $\tilde u$ 把 FiBE 变成等价的单腿多步法

$$
\sum_{j=1}^{2}\Bigl[\prod_{i=1}^{j-1}(t_{n+m}-t_{n+m-i})\Bigr]\delta^{j}y
\;=\;f\Bigl(t_{n+m},\;y_{n+m}
+\frac{\eta^{(2)}}{1-\eta^{(2)}c^{(2)}_{m}}\,\delta^{2}y\Bigr),
\qquad m=2 .
$$

其左端正是变步长 BDF2 差分算子——这正是编号 69 所利用的结构事实。

### 结果

本文为变步长 FiBE 所对应的单腿多步格式建立**离散梯度结构**，据此在耗散情形下给出离散能量耗散律，并在下述步长比条件下给出 $L^2$ 稳定性与 $L^2$ 误差估计。论文声明这「似乎是变步长时间滤波刚性求解器的**第一个**能量稳定性与 $L^2$ 范数误差估计」。分析装置是**两类新的离散正交卷积核**——「两类」与 OLM 的结构相合：左端的 BDF2 卷积与右端的单腿／插值卷积各需要自己的正交对偶。

处理的对象是**线性抛物方程**而非 Allen-Cahn，因此这一篇不含最大值原理的结论。

### 步长比条件

稳定性与 $L^2$ 误差估计在**双边的「实用」步长比约束**

$$
\tfrac12\ \le\ \frac{\tau_k}{\tau_{k-1}}\ \le\ 2
\qquad(k\ge2)
$$

下成立。

> [!warning] $[1/2,\,2]$ 与其余门槛不是同一类对象
> **第一，它是双边的。** 其他各篇的结果都只限制步长比的**上界**（$r_k<1+\sqrt2$、$r_k<3.561$、$r_k<1.4877$），而这里步长比还被 $1/2$ 从下方限制——也就是说步长既不能放大太快，也不能**收缩**太快。这反映的是变步长 FiBE 失去了 $A$-稳定性。
>
> **第二，$[1/2,2]$ 不是一个尖锐的解析门槛。** $1+\sqrt2$ 是 $r^2-2r-1=0$ 的根，$3.561$ 是 $r^2-3r-2=0$ 的根，$1.4877$ 是一条显式方程的唯一正根；而 $[1/2,2]$ 是自适应程序里常用的经验保护区间——DeCaria-Guzel-Layton-Li 就把自适应实验限制在「最大二、最小二分之一，这是变步长方法中常见的经验做法」。论文自己称之为「实用」约束，正是出于这个原因。
>
> **因此不应把 $1/2\le r_k\le2$ 与 $1+\sqrt2$、$3.561$、$4.8645$、$1.4877$ 并列引用，好像它们是同一种量。**

### 与其他论文的关系

它是同一纲领里的**线性抛物、时间滤波**分支：与编号 48、52、58、67 用同一套离散梯度结构加 DOC 核工具箱，只是作用对象从 BDF 公式本身换成了一个**预估-校正**方法。参考文献直接印证这条谱系：它引用了编号 48（_SIAM J. Numer. Anal._ 58:2294-2314）、编号 58（_CSIAM Trans. Appl. Math._ 3:318-334）与编号 67（_J. Comput. Math._，DOI `10.4208/jcm.2207-m2022-0020`），以及 Liao-Zhang（_Math. Comp._ 90:1207-1226）、Li-Liao（_SIAM J. Numer. Anal._ 60:2253-2272）、Liao-Ji-Wang-Zhang（_J. Sci. Comput._ 92:52，即 $4.8645$ 的出处）与 Liao-Ji-Zhang（_IMA J. Numer. Anal._ 42:649-679）。

**单腿改写是通往 BDF2 各篇的桥**：FiBE 所对应的 OLM 左端恰好是变步长 BDF2，因此编号 48、52 的梯度结构机器可以直接转移，只是非线性项与外力在一个滤波过的宗量处求值。它的步长比结论则是这一组里的**异类**：双边且经验，而非单侧的尖锐代数门槛，这正反映了变步长 FiBE 失去 $A$-稳定性这一事实。

下游方面，后续论文《Energy dissipation laws of time filtered BDF methods up to fourth-order for the molecular beam epitaxial equation》明确在本文基础上把理论推广到 FiBDF-$k$ 与分子束外延模型；Wang-Liao-Zhao（_Numer. Math. Theor. Meth. Appl._ 16(1) (2023) 165-181）则处理**常步长**滤波后向 Euler 用于带斜率选择的分子束外延。

## 74：把工具本身作为研究对象

### 直觉

前面五篇每一篇都要回答同一个问题：某族卷积核构成的实二次型

$$
\sum_{k=1}^{n}w_k\sum_{j=1}^{k}a^{(k)}_{k-j}w_j
$$

什么时候正定？每一篇都是就自己的核族现证一遍。本文把这个问题**作为独立对象**来解决。

在**均匀**网格上有完整的经典答案：由 Toeplitz-Carathéodory 定理，若 $\hat a(z)=\sum_{k\ge0}a_kz^k$ 在开单位圆盘 $\mathcal D_z$ 上解析，则该二次型半正定**当且仅当** $\mathrm{Re}[\hat a(z)]\ge0$ 于 $\mathcal D_z$；López-Marcos（1990，Prop. 5.2）从中提炼出易用的充分条件

$$
a_j\ge0,\qquad a_{j-1}\ge a_j,\qquad a_{j-1}-a_j\ge a_j-a_{j+1} .
$$

在**变**网格上核多了一个指标 $a^{(n)}_{n-k}$（取值依赖于当前时刻 $t_n$），生成函数机器整个失效，此前只有针对特殊核的「保半正定」构造。本文补上了缺失的一般判据。

**判据是纯代数的，与背后的连续核没有任何显式联系——这正是要点。** 证明的策略也很干净：对一般核族同时造出 DOC 与 DCC 两套核，先从条件推出 DCC 核的单调性，再用「DOC 与原核同正定」这一等价性收口。

### 问题设定

对变步长核 $\{a^{(n)}_{n-k}\}_{k=1}^{n}$，考察上式对任意 $\{w_1,\dots,w_n\}$ 的符号。结构性假设：$a^{(n)}_j\ne0$（$0\le j\le n-1$），且若对某个 $2\le n_0\le n-1$ 有 $a^{(n)}_{n_0}=0$，则 $a^{(n)}_j=0$ 对一切 $j\ge n_0$ 成立。

### 推导

**主定理（定理 1.1）。** 对固定的 $n\ge2$，若

$$
\textbf{C1}:\ a^{(n)}_{j}>0\quad(0\le j\le n-1);
\qquad
\textbf{C2}:\ a^{(n-1)}_{j-1}>a^{(n)}_{j}\quad(1\le j\le n-1);
$$

$$
\textbf{C3}:\ a^{(n-1)}_{j-1}a^{(n)}_{j+1}\ \ge\ a^{(n-1)}_{j}a^{(n)}_{j}\quad(1\le j\le n-2);
\qquad
\textbf{C4}:\ a^{(n)}_{j-1}\ \ge\ a^{(n)}_{j}\quad(1\le j\le n-1),
$$

则离散卷积核 $a^{(n)}_{n-k}$ **正定**。

四条的含义分别是：C1 是正性；C4 是固定层内对滞后指标的单调性；C2 是**跨层**的单调性（这一条在均匀网格上是平凡的，正是变步长带来的新要求）；C3 是一条带对数凸性／全正性味道的跨层不等式，扮演的是均匀网格上 López-Marcos 那条凸性条件 $a_{j-1}-a_j\ge a_j-a_{j+1}$ 的角色。

**工具一：DOC 核。** 本文是**一般变步长**定义的规范出处：

$$
\theta_0^{(n)}:=\frac{1}{a^{(n)}_0},
\qquad
\theta_{n-k}^{(n)}:=-\frac{1}{a^{(k)}_0}\sum_{j=k+1}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k},
\qquad k=n-1,n-2,\dots,1,
$$

满足离散正交恒等式 $\sum_{j=k}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}\equiv\delta_{nk}$（$1\le k\le n$）。引理 2.1 给出相互正交性

$$
\sum_{j=k}^{n}a^{(n)}_{n-j}\theta_{j-k}^{(j)}=\delta_{nk}
\qquad\text{与}\qquad
\sum_{j=k}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}=\delta_{nk},
\qquad 1\le k\le n,
$$

从而 **$\theta_{n-j}^{(n)}$ 正（半）定当且仅当 $a^{(n)}_{n-k}$ 正（半）定**。DOC 的原始构造归于 Liao-Zhang（那里用于变步长 BDF2 的 $L^2$ 稳定性）。

**工具二：DCC 核。** 本文的做法是把 DCC 核**从 DOC 核沿层指标求和**定义出来：

$$
p_{n-k}^{(n)}:=\sum_{j=k}^{n}\theta_{j-k}^{(j)},
\qquad 1\le k\le n,
$$

因而

$$
\theta_0^{(n)}=p_0^{(n)},
\qquad
\theta_{n-k}^{(n)}=p_{n-k}^{(n)}-p_{n-k-1}^{(n-1)}\quad(1\le k\le n-1),
$$

并满足离散互补恒等式

$$
\sum_{j=k}^{n}p_{n-j}^{(n)}a^{(j)}_{j-k}\equiv 1,
\qquad 1\le k\le n .
$$

证明是令 $\Xi_k^{(n)}:=\sum_{j=k}^{n}p^{(n)}_{n-j}a^{(j)}_{j-k}$，验证 $\Xi_k^{(k)}=1$ 与 $\Xi_k^{(n)}=\Xi_k^{(n-1)}$，再归纳。**这两个名字因此是名副其实的：DOC 把卷积求逆（给出 $\delta_{nk}$），DCC 把卷积求补（给出 $1$）。** 两者的闭式（引理 2.3 与 2.6）用辅助序列 $\psi^{(m)}_j$、$\chi^{(k)}_\ell$ 表出：

$$
\theta_{j-k}^{(j)}=-\frac{1}{a_0^{(j)}}\psi_1^{(k+1)}
\prod_{\ell=k+2}^{j}\bigl(\chi_2^{(j-\ell)}\psi_2^{(\ell)}-\psi_1^{(\ell)}\bigr)
\quad(j\ge k+1),
$$

$$
p_{n-k}^{(n)}=\frac{1}{a_0^{(k)}}-\psi_1^{(k+1)}\sum_{j=k+1}^{n}\frac{1}{a_0^{(j)}}
\prod_{\ell=k+2}^{j}\bigl(\chi_2^{(j-\ell)}\psi_2^{(\ell)}-\psi_1^{(\ell)}\bigr).
$$

**证明的引擎是 DCC 单调性（引理 2.7）。** 对 $n\ge2$：若 $a^{(n)}_{n-k}$ 满足 **C1-C3**，则

$$
p_0^{(n)}>0,
\qquad
p_0^{(n-1)}>p_1^{(n)},
\qquad
p_{j-1}^{(n-1)}\ge p_j^{(n)}\ \ (2\le j\le n-1);
$$

再由 **C4** 补上余下的单调性，即可断定 $\theta$ 正定。定理 1.1 就是把 DCC 单调性与引理 2.1 的 DOC 等价性合起来。

**步长比限制：没有。** C1-C4 是对**核**的条件，在**任意非均匀网格上**成立，对 $r_k=\tau_k/\tau_{k-1}$ 不加任何限制。这是本文的卖点之一，也是它与编号 58（Grenander-Szegő，均匀网格 Toeplitz）与编号 67（显式离散梯度结构，变步长 BDF3）的分水岭：三条路里只有这一条经得起任意网格。

### 定理

**定理 1.1**（C1-C4 $\Rightarrow$ 正定）如上。论文称它是「关于变步长卷积系数正定性的第一个具有简单代数条件的结果」。四个应用如下。

**应用一：变步长 L1 公式（命题 4.1）。** 记 $\omega_\gamma(t):=t^{\gamma-1}/\Gamma(\gamma)$，Riemann-Liouville 积分 $(\mathcal I^\gamma v)(t)=\int_0^t\omega_\gamma(t-s)v(s)\mathrm ds$，Caputo 导数 $(\partial_t^\alpha v)(t)=(\mathcal I^{1-\alpha}v')(t)$（$0<\alpha<1$）。非均匀网格上的 **L1 公式**把 $v'$ 在 $(t_{k-1},t_k)$ 上换成常数 $\nabla_\tau v^k/\tau_k$：

$$
(\partial_\tau^\alpha v)^n:=\sum_{k=1}^{n}c^{(n,\alpha)}_{n-k}\nabla_\tau v^{k},
\qquad
c^{(n,\alpha)}_{n-k}:=\frac1{\tau_k}\int_{t_{k-1}}^{t_k}\omega_{1-\alpha}(t_n-s)\,\mathrm ds .
$$

则

$$
c^{(n,\alpha)}_{j}>0,
\quad
c^{(n,\alpha)}_{j-1}>c^{(n,\alpha)}_{j},
\quad
c^{(n-1,\alpha)}_{j-1}>c^{(n,\alpha)}_{j},
\quad
c^{(n-1,\alpha)}_{j-1}c^{(n,\alpha)}_{j+1}>c^{(n-1,\alpha)}_{j}c^{(n,\alpha)}_{j},
$$

即 C1-C4 全部成立，因而**变步长 L1 核在任意非均匀网格上正定**。证明很漂亮，值得记下：积分中值定理给出 $c^{(n,\alpha)}_{n-k}=\omega_{1-\alpha}(t_n-s_{nk})$（某个 $s_{nk}\in[t_{k-1},t_k]$），这直接给出 C1 与 C4——**离散核字面上就是在采样连续核 $\omega_{1-\alpha}$**；再对 $c_{n,k}(\mu):=\frac1{\tau_k}\int_{t_{k-1}}^{t_{k-1}+\mu\tau_k}\omega_{1-\alpha}(t_n-s)\mathrm ds$ 用 Cauchy 中值定理，得

$$
\psi^{(n,\alpha)}_{n-k}:=\frac{c^{(n,\alpha)}_{n-k}}{c^{(n-1,\alpha)}_{n-1-k}}
=\Bigl(\frac{t_{n-1}-t_{k-1}-\xi_{1k}\tau_k}{t_n-t_{k-1}-\xi_{1k}\tau_k}\Bigr)^{\alpha},
\qquad \xi_{1k}\in(0,1),
$$

从而

$$
\Bigl(\frac{t_{n-1}-t_k}{t_n-t_k}\Bigr)^{\alpha}<\psi^{(n,\alpha)}_{n-k}
<\Bigl(\frac{t_{n-1}-t_{k-1}}{t_n-t_{k-1}}\Bigr)^{\alpha},
\qquad
0<\psi_1^{(n,\alpha)}<\psi_2^{(n,\alpha)}<\cdots<\psi_{n-1}^{(n,\alpha)}<(t_{n-1}/t_n)^{\alpha}<1 .
$$

$\psi^{(n,\alpha)}_j$ 关于 $j$ 单调递增就正好是 C2 与 C3。

**（尖锐性，注 4.1）** C1-C4 充分但**不必要**。Ji-Liao-Gong-Zhang 的二阶 **L1$^{+}$** 公式，核为

$$
\bar c_{n-k}^{(n)}:=\frac{1}{\tau_n\tau_k}\int_{t_{n-1}}^{t_n}
\int_{t_{k-1}}^{\min\{t,t_k\}}\omega_{1-\alpha}(t-s)\,\mathrm ds\,\mathrm dt,
$$

**确实**半正定（从连续核继承），但条件 $\bar c_0^{(n)}\ge\bar c_1^{(n)}$ 对某些 $n$ 会失败。寻找充要的代数条件仍是公开问题。

**应用二：时间分数阶 Allen-Cahn 的能量稳定性（命题 4.2）。** 对 $\partial_t^\alpha u=\varepsilon^2\Delta u-F'(u)$、$F(u)=\frac14(1-u^2)^2$，连续能量律 $E[u(t)]\le E[u(0)]$ 归于 Tang 与 Tang-Yu-Zhou，即[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 40]]。L1 核的正定性给出

$$
\sum_{k=1}^{n}\nabla_\tau v^{k}\,(\partial_\tau^\alpha v)^{k}
=\sum_{k=1}^{n}\nabla_\tau v^{k}\sum_{j=1}^{k}c^{(k,\alpha)}_{k-j}\nabla_\tau v^{j}>0
\qquad(\nabla_\tau v^k\not\equiv0).
$$

用到 Ji-Liao-Zhang 的一阶稳定化显隐格式

$$
(\partial_\tau^\alpha u)^n=\varepsilon^2D_hu^n-F'(u^{n-1})-S(u^n-u^{n-1}),
\qquad n\ge1,
$$

结论是：**若稳定化参数 $S\ge2$**，格式保离散最大值原理，且 $E_h^n\le E_h^0$（$n\ge1$），其中

$$
E_h^n=-\frac{\varepsilon^2}{2}(u^n)^TD_hu^n+\sum_{x_h\in\Omega_h}F(u_h^n).
$$

论文注明均匀网格上的对应结果已知（Ji-Liao-Zhang；Tang-Yu-Zhou），而这个**非均匀网格**版本似乎是新的。

> [!note] 三条能量律是三个不同的陈述
> 编号 40 证的是**分数阶（非局部、积分型）**能量律；[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 57]] 证的是**变分**能量的微分型律；这里证的是**离散**能量 $E_h^n\le E_h^0$。三者的对象与形式都不同，不能互相替代着引用。

**应用三：Riemann-Liouville 积分与分数阶波方程。** 对 $\partial_t u=\mathcal I^\gamma\Delta u+f$（解在 $t=0$ 附近有 $\partial_{tt}u\sim\mathcal O(t^{\gamma-1})$，故分级网格是自然的对策），后向 Euler 格式 $\partial_\tau u^{n-\frac12}=(\mathcal I^\gamma_\tau\Delta u)^n+f(x,t_n)$ 配中点求积 $(\mathcal I^\gamma_\tau v)^n=\sum_{k=1}^{n}c^{(n,1-\gamma)}_{n-k}\tau_kv^{k-\frac12}$，在**一般**非均匀网格上满足

$$
\|u^n\|_{L^2}\le\|u^0\|_{L^2}+\sum_{k=1}^{n}\tau_k\|f(t_k)\|_{L^2}.
$$

**应用四：弱奇异 Volterra 方程（命题 4.3）。** 对 $\partial_t u=\mathcal K_t^{(\beta)}\Delta u+f$，$(\mathcal K^{(\beta)}_tv)(t)=\int_0^t\kappa_\beta(t-s)v(s)\mathrm ds$，中点核 $\kappa^{(n,\beta)}_{n-k}:=\frac1{\tau_k}\int_{t_{k-1}}^{t_k}\kappa_\beta(t_n-s)\mathrm ds$：**若 $\kappa_\beta>0$、$\kappa_\beta'<0$、$\kappa_\beta''\ge0$**，则 $\kappa^{(n,\beta)}_{n-k}$ 满足 C1-C4（论证用到 $y=\kappa_\beta(t_n-x)/\kappa_\beta(t_{n-1}-x)$ 在 $t_{n-1}<x<t_n$ 上递减），因而正定。

**这四个应用一起说明了把工具单列出来的价值**：同一条判据一次性覆盖了分数阶 Caputo 导数、Riemann-Liouville 积分与一般弱奇异核，而且**对网格不加任何限制**。

### 数值实验

**本文没有数值实验。** 它是纯分析论文，全文唯一的图形是核的示意图。因此它给出的所有结论都是定理，没有实测数据可以对照——这一点与编号 58 相同。

### 与其他论文的关系

它是这批论文分数阶那一半**底下的一般理论**。这里证明的变步长 L1 正定性，正是[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 43]]（非均匀网格上时间分数阶 Allen-Cahn 的 Alikhanov／L1 型核）与[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 57]]（Riemann-Liouville 改写所用的 L1$_R$ 核）以各种面貌用到的性质——编号 57 在证明其 DOC 核单调性时直接调用本文的三条代数判据，因此本文在逻辑上先于编号 57，尽管出版更晚。它所形式化的 DOC 构造，与编号 48、52 在 BDF2 层面、编号 58 在 BDF-$k$ 层面、编号 67 在 BDF3 层面用的是同一个；它形式化的 DCC 构造则是编号 43 互补核界的出处。**本文是把一般变步长定义与那两条恒等式（DOC 的 $\equiv\delta_{nk}$ 与 DCC 的 $\equiv1$）对任意核族干净地陈述出来的地方。**

**把编号 74 单列出来的意义在于**：这一系列工作的技术核心不是某个格式，而是「历史项汇成的实二次型何时正定」这一代数问题。一旦这个问题有了独立的判据，同一套论证就可以搬到三阶 BDF、时间滤波 Euler、分数阶 L1 逼近与隐显 Runge-Kutta 上。

## 六篇的关系

| 编号 | 对象                           | 网格   | 门槛                             | 核心工具                     | 数值实验                             |
| ---- | ------------------------------ | ------ | -------------------------------- | ---------------------------- | ------------------------------------ |
| 48   | Allen-Cahn，BDF2               | 变步长 | S1 $3.561$；S0 $1+\sqrt2$        | 修正能量；核重排 KRC；DCC    | 随机网格上二阶；气泡合并；粗化       |
| 52   | 分子束外延（无斜率选择），BDF2 | 变步长 | $3.561$（能量与 $L^2$ 同门槛）   | DOC 核；$\mathcal M_r$       | 随机网格：$\max r_k$ 到 $850$ 仍二阶 |
| 58   | 线性反应扩散，BDF-$k$          | 均匀   | 无（$3\le k\le5$）               | DOC 核；Grenander-Szegő      | 无；仅有引理的示意图                 |
| 67   | 扩散方程，BDF3                 | 变步长 | $1.4877$（充分）／$1.69$（必要） | 离散梯度结构 $G$；变步长 DOC | 三阶；半数层违反门槛仍稳定           |
| 69   | 线性抛物，滤波 Euler           | 变步长 | $[1/2,2]$（双边，经验）          | 单腿改写；两类新 DOC 核      | —                                    |
| 74   | 二次型本身                     | 任意   | 无（给出代数判据 C1-C4）         | DOC 与 DCC 的一般理论        | 无                                   |

**读这张表的方式**：门槛一列的数值大小不可横向比较——它们分别是不同格式、不同结论的条件。真正可比的是最后一列：**每当论文自己做了随机网格实验，观测到的稳健范围都远宽于定理保证的范围**（编号 52 的 $850.80$、编号 67 的半数层违反），这是六篇里最一致的经验事实。六篇中没有一篇给出打破结论的反面算例；唯一那个能把最大值界真的打破的实验在[[computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|编号 43]]，而它破坏的是**步长-空间耦合条件**而非步长比条件。两相对照说明：耦合条件比步长比条件更接近必要。

## 原文

- H.-l. Liao, T. Tang, and T. Zhou, [_On energy stable, maximum-principle preserving, second-order BDF scheme with variable steps for the Allen-Cahn equation_](https://doi.org/10.1137/19M1289157), SIAM J. Numer. Anal. 58(4) (2020), pp. 2294-2314（预印本 [arXiv:2003.00421](https://arxiv.org/abs/2003.00421)）。
- H.-l. Liao, X. Song, T. Tang, and T. Zhou, [_Analysis of the second-order BDF scheme with variable steps for the molecular beam epitaxial model without slope selection_](https://doi.org/10.1007/s11425-020-1817-4), Sci. China Math. 64 (2021), pp. 887-902（预印本 [arXiv:2008.03185](https://arxiv.org/abs/2008.03185)）。
- H.-l. Liao, T. Tang, and T. Zhou, [_A new discrete energy technique for multi-step backward difference formulas_](https://doi.org/10.4208/csiam-am.SO-2021-0032), CSIAM Trans. Appl. Math. 3 (2022), pp. 318-334（预印本 [arXiv:2102.04644](https://arxiv.org/abs/2102.04644)）。
- H.-l. Liao, T. Tang, and T. Zhou, [_Discrete energy analysis of the third-order variable-step BDF time-stepping for diffusion equations_](https://doi.org/10.4208/jcm.2207-m2022-0020), J. Comput. Math. 41 (2023), pp. 325-344（预印本 [arXiv:2204.12742](https://arxiv.org/abs/2204.12742)）。
- H.-l. Liao, T. Tang, and T. Zhou, [_Stability and convergence of the variable-step time filtered backward Euler scheme for parabolic equations_](https://doi.org/10.1007/s10543-023-00982-y), BIT Numer. Math. 63 (2023), 39。
- H.-l. Liao, T. Tang, and T. Zhou, [_Positive definiteness of real quadratic forms resulting from the variable-step L1-type approximations of convolution operators_](https://doi.org/10.1007/s11425-022-2229-5), Sci. China Math. 67 (2024), pp. 237-252（预印本 [arXiv:2011.13383](https://arxiv.org/abs/2011.13383)）。

方法出处与背景（上文引用但不属于本专题）：

- V. DeCaria, S. Guzel, W. Layton, and Y. Li, [_A variable stepsize, variable order family of low complexity_](https://doi.org/10.1137/19M1291666), SIAM J. Sci. Comput. 43(3) (2021), pp. A2130-A2160（预印本 [arXiv:1810.06670](https://arxiv.org/abs/1810.06670)）——编号 69 所分析的滤波格式的出处。
- H.-l. Liao, B. Ji, L. Wang, and Z. Zhang, [_Mesh-robustness of an energy stable BDF2 scheme with variable steps for the Cahn-Hilliard model_](https://doi.org/10.1007/s10915-022-01923-7), J. Sci. Comput. 92 (2022), 52（预印本 [arXiv:2102.03731](https://arxiv.org/abs/2102.03731)）——$4.864$ 的出处。
