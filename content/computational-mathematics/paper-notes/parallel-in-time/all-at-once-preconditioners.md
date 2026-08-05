---
title: 全时间预条件与谱分析
description: 编号 59、65、71、84、85：把「构造一个预条件子」变成「证明它的谱受控」
lang: zh
translation: en/computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners
tags:
  - 论文笔记
  - 时间并行
  - 预条件
---

> [!note] 本页覆盖
> 编号 **59**（_Adv. Comput. Math._ 48:16, 2022）、**65**（_SIAM J. Matrix Anal. Appl._ 43(3), 2022）、**71**（_SIAM J. Matrix Anal. Appl._ 44(4), 2023）、**84**（_J. Sci. Comput._ 103:82, 2025）、**85**（_Acta Numer._ 34, 2025）。
>
> 核实深度差别很大，本页据此分层书写。编号 **59** 有 arXiv 预印本且其 LaTeX 源码已逐式核对，因此定理、证明技术、条件数结论与摘要报告的加速比都可以完整给出。编号 **65** 的摘要已核实，其主定理被同一作者群的 _Acta Numerica_ 综述完整转述，本页据此给出定理与假设，并用一段自己的标量通道计算确定界的方向（综述把两个端点印反了）。编号 **71** 的摘要在所有公开来源中数学符号均被剥离，**$\alpha$ 的定标律因此只能记为形式而不给指数**。编号 **84** 的 Springer 摘要已核实，构造思路清楚，但聚集半径与实验参数未核实。编号 **85** 的逐节精读见[[computational-mathematics/knowledge-notes/time-parallelization/index|时间并行专题]]，本页只给出它在这条线索中的位置。

## 59：换一个时间离散，使 $V$ 天生条件良好

### 直觉

时间推进是串行的，因为第 $n$ 层依赖第 $n-1$ 层。把所有时间层叠成一个长向量，这条依赖表现为一个块下三角矩阵，而时间推进恰好就是对它做块前代——前代无法并行，它本质上就是一条链。

对角化改变的是**观察这条链所用的坐标**。若时间矩阵可分解为 $B=VDV^{-1}$，在 $V^{-1}$ 给出的坐标里各时间自由度之间的耦合完全散开：「第 $n$ 层等第 $n-1$ 层」变成「第 $n$ 个**模态**独立求解」。每个模态只留下一个复标量 $\lambda_n$ 与一个复移位的空间问题 $(\lambda_nI_x+\cdots)x=y$，这 $N_t$ 个问题互不相干，可以同时算。

代价全部集中在一句话里：$V^{-1}$ 与 $V$ 是在浮点下作用的，因此计算结果携带 $O(\epsilon\,\mathrm{Cond}_2(V))$ 量级的舍入污染，$\epsilon$ 是机器精度。**理论上可对角化但 $\mathrm{Cond}_2(V)=10^{14}$ 的时间矩阵在双精度下毫无用处。**

这就是 Maday 与 Rønquist 路线的困境。标准等步长离散给出的 $B$ 根本不可对角化（是一个 Jordan 块），他们的对策是把步长取成互不相同的，使特征值分开。但步长越接近，$B$ 越接近亏损，$\mathrm{Cond}_2(V)$ 越大；步长差别越大，离散精度越差。两头夹击的结果是方法只能用到约 $20$ 到 $25$ 个时间点。

**这一篇把设计问题反过来问。** 前面的工作都是「给定时间离散，设法安排步长使 $V$ 可用」；这一篇问的是「什么样的时间离散会**天生**给出条件良好的 $V$」。答案是一个**边值方法**：一个不能作为时间推进格式运行、只在全时间意义下成立的格式。它在等步长下就可对角化，而且 $\mathrm{Cond}_2(V)$ 只按 $n^2$ 增长。付出的代价是放弃时间推进的解释——但既然目标本来就是一次性全时间求解，这个代价是免费的。

### 问题设定

对 $u'(t)+Au(t)=g(t)$ 形成全时间系统 $\mathcal K\boldsymbol u:=(B\otimes I_x+I_t\otimes A)\boldsymbol u=\boldsymbol b$，若 $B=VDV^{-1}$ 则

$$
\mathcal K=(V\otimes I_x)\bigl(D\otimes I_x+I_t\otimes A\bigr)(V^{-1}\otimes I_x),
$$

于是三步求解：$(V^{-1}\otimes I_x)$ 作用于右端、$n$ 个完全解耦的空间求解、$(V\otimes I_x)$ 作用回来。**障碍在于标准等步长离散给出的 $B$ 根本不可对角化**：等步长后向 Euler 的 $B$ 是下二对角 Toeplitz 矩阵，即单个 Jordan 块；多步法的 $B$ 是下三角 Toeplitz 矩阵，同样亏损。

Maday 与 Rønquist（2008）的对策是取**互不相同**的步长 $\{\Delta t_j\}$。此时 $B$ 的对角元 $1/\Delta t_j$ 两两不同，$B$ 有 $n$ 个互异特征值因而可对角化，$D=\mathrm{diag}(1/\Delta t_1,\dots,1/\Delta t_n)$。但舍入误差满足

$$
\texttt{舍入误差}=\mathcal O\bigl(\epsilon\,\mathrm{Cond}_2(V)\bigr),
\qquad \epsilon=\text{机器精度}\ (\text{双精度 } 2.22\times10^{-16}),
$$

而对几何递增步长 $\Delta t_j=\Delta t_1\tau^{\,j-1}$，比值 $\tau>1$ 被夹在两边：$\tau\to1$ 使 $B$ 接近亏损、$\mathrm{Cond}_2(V)$ 爆炸；$\tau\gg1$ 使步长指数增长、离散精度被破坏。

这不是定性描述，两边都有闭形式。几何步长配线性 $\theta$-方法时 $V$ 与 $V^{-1}$ **都是单位下三角 Toeplitz 矩阵**，

$$
V=\mathbb T(p_1,\dots,p_{n-1}),
\qquad
V^{-1}=\mathbb T(q_1,\dots,q_{n-1}),
\qquad
\mathbb T(a_1,\dots,a_{n-1}):=
\begin{bmatrix}1\\a_1&1\\\vdots&\ddots&\ddots\\a_{n-1}&\cdots&a_1&1\end{bmatrix},
$$

后向 Euler（$\theta=1$）情形

$$
p_j=\frac{1}{\prod_{i=1}^{j}\bigl(1-\tau^{\,i}\bigr)},
\qquad
q_j=(-1)^j\,\tau^{\,j(j-1)/2}\,p_j .
$$

> [!warning] 这里的参数是几何**比值**，不是增量
> 综述在转述这条闭形式时写成 $p_j=1/\prod_i(1-\varrho^{\,i})$，其中 $\varrho=\tau-1$ 是**增量**。这不可能对：按增量写时 $\varrho\to0$ 会让 $p_j\to1$、$V$ 趋于良态，而这恰好是所有人都同意 $V$ 变得灾难性病态的极限。由 $BV=VD$ 的 $(2,1)$ 元素可直接读出 $p_1=1/(1-\tau)$，确认应当用比值。本站[[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i|ParaDiag-I 一章]]对此有完整核对。

把两侧的误差量化后可以解出最优拉伸。写 $\varrho=\tau-1$，对 $\boldsymbol u'=A\boldsymbol u+\boldsymbol g$（$\sigma(A)\subset\mathbb R^-$，$|\lambda(A)|\le\lambda_{\max}$），记 $\boldsymbol u_{N_t}(\varrho)$ 为几何网格上的后向 Euler 解、$\boldsymbol u_{N_t}(0)$ 为等步长解、$\tilde{\boldsymbol u}_n(\varrho)$ 为三步对角化在浮点下实际算出的量，则

$$
\underbrace{\|\boldsymbol u_{N_t}(\varrho)-\boldsymbol u_{N_t}(0)\|\lesssim C(\lambda_*T,N_t)\,\varrho^{2}}_{\text{拉伸的截断代价，随 }\varrho\text{ 增大}},
\qquad
\underbrace{\|\tilde{\boldsymbol u}_n(\varrho)-\boldsymbol u_n(\varrho)\|\lesssim\epsilon\,\frac{N_t^2(2N_t+1)(N_t+\lambda_{\max}T)}{\phi(N_t)}\,\varrho^{-(N_t-1)}}_{\text{对角化的舍入代价，随 }\varrho\to0\text{ 爆炸}},
$$

其中 $C(x,N_t):=\tfrac{N_t(N_t^2-1)}{24}r(x/N_t,N_t)$、$r(\tilde x,N_t):=\bigl(\tfrac{\tilde x}{1+\tilde x}\bigr)^2(1+\tilde x)^{-N_t}$、$\lambda_*:=N_t\tilde x_*/T$（$\tilde x_*$ 是 $r(\cdot,N_t)$ 在 $[0,\infty)$ 上的最大值点），而

$$
\phi(N_t):=
\begin{cases}
\bigl(\tfrac{N_t}{2}\bigr)!\bigl(\tfrac{N_t}{2}-1\bigr)!, & N_t\ \text{为偶},\\[4pt]
\Bigl(\bigl(\tfrac{N_t-1}{2}\bigr)!\Bigr)^2, & N_t\ \text{为奇}.
\end{cases}
$$

两者平衡给出

$$
\varrho_{\rm opt}=\left(\epsilon\,\frac{N_t^{2}(2N_t+1)(N_t+\lambda_{\max}T)}{\phi(N_t)\,C(\lambda_*T,N_t)}\right)^{\frac{1}{N_t+1}} .
$$

**关键是那个 $\varrho^{-(N_t-1)}$：舍入界按时间步数指数退化。** 所以 Maday-Rønquist 路线本质上是一个**短窗口**方法，实用范围约 $N_t\approx20$ 到 $30$，长时间问题只能逐窗口做。取消这个上限就是编号 59 的目标。

### 推导

**第一步：换一个只在全时间层面成立的格式。** 论文的时间离散是混合的：前 $n-1$ 步用中心（跳蛙）差分，只有最后一步用隐式 Euler：

$$
\begin{cases}
\dfrac{u_{j+1}-u_{j-1}}{2\Delta t}+Au_j=g_j, & j=1,2,\dots,n-1,\\[6pt]
\dfrac{u_n-u_{n-1}}{\Delta t}+Au_n=g_n. &
\end{cases}
$$

这是一个**边值方法**：它**不能**作为时间推进格式运行，因为中心格式那样用是不稳定的；只有一次性全时间求解才有意义。这一构造出自 Axelsson 与 Verwer（1985），他们证明即使最后一步只有一阶，同时得到的全部解仍具有**一致二阶精度**；更早 Fox（1954）与 Fox-Mitchell（1957）用过 BDF2 作末步。Brugnano、Mazzia 与 Trigiante（1993）用迭代方法求解相应的全时间系统，本文改为用对角化**直接**求解。

关键收益是步长可以**等长**，因此完全没有步长比参数：

$$
B=\frac{1}{\Delta t}
\begin{bmatrix}
0&\tfrac12&&&\\
-\tfrac12&0&\tfrac12&&\\
&\ddots&\ddots&\ddots&\\
&&-\tfrac12&0&\tfrac12\\
&&&-1&1
\end{bmatrix},
\qquad
\boldsymbol b=\begin{bmatrix}\tfrac{u_0}{2\Delta t}+g_1\\ g_2\\ \vdots\\ g_n\end{bmatrix}.
$$

只有 $u_0$ 进入右端，全部时间步一次解出。注意 $B$ 除最后一行外是反对称 Toeplitz 的——**正是最后一行破坏了对称性，也正是它使谱可以闭形式算出**。二阶（波型）方程 $u''+Au=g$ 用同一离散给出 $\bigl(B^2\otimes I_x+I_t\otimes A\bigr)\boldsymbol u=\boldsymbol b$，而 $B^2$ 与 $B$ **共用同一个** $V$，因此同一套条件数分析直接适用。这就是一篇论文能同时覆盖一阶与二阶问题的原因。

**第二步：把特征问题化归到 Chebyshev 多项式。** 结论本身容易记，机制才是这篇的实质：**特征问题可以化归到 Chebyshev 多项式，于是 $V$ 与 $V^{-1}$ 都有显式表达式。** 取缩放后的 $\mathcal B:=\Delta t\,B$，并记第一类与第二类 Chebyshev 多项式

$$
T_n(x)=\cos(n\arccos x),
\qquad
U_n(x)=\frac{\sin\bigl[(n+1)\arccos x\bigr]}{\sin(\arccos x)} .
$$

则 $\mathcal B$ 的 $n$ 个特征值为 $\lambda_j=\mathrm ix_j$，其中 $\{x_j\}_{j=1}^n$ 是**特征方程**

$$
U_{n-1}(x)-\mathrm i\,T_n(x)=0
$$

的 $n$ 个根，对应特征向量的分量为

$$
p_{j,k}=\mathrm i^{\,k}\,U_k(x_j),\qquad k=0,1,\dots,n-1,\qquad p_{j,0}=1 .
$$

**第三步：证明根的结构。** 这些根的结构恰好提供了所需的三件事：全部根**单重**（因此 $\mathcal B$ 可对角化）、都是虚部为负的复数（因此特征值落在正确的半平面）、且模小于 $1+1/\sqrt{2n}$（这一条正是控制 $\mathrm{Cond}_2(V)$ 的关键）；此外 $x$ 是根则 $-\bar x$ 也是根。证明的做法是代入

$$
y=x+\sqrt{x^2-1}
\qquad\Longrightarrow\qquad
T_n=\tfrac12\bigl(y^n+y^{-n}\bigr),
\qquad
U_{n-1}=\frac{y^n-y^{-n}}{y-y^{-1}},
$$

于是特征方程 $U_{n-1}-\mathrm iT_n=0$ 化为

$$
y^{2n}=-\frac{(y-\mathrm i)^2}{(y+\mathrm i)^2},
$$

再用 $|y|>1$ 与 Chebyshev 的 Pythagoras 恒等式 $T_n^2(x)+(1-x^2)U_{n-1}^2(x)=1$ 定出模的上界。**这一步是全文的技术核心**：$y$-变量把一个含两类 Chebyshev 多项式的超越方程变成一个纯代数的 $2n$ 次方程，右端是一个模为 $1$ 的 Möbius 型量的平方，因此 $|y|^{2n}$ 被夹住，根的位置随之被夹住。

**第四步：把条件数归约到 Vandermonde 型矩阵。** 特征向量矩阵分解为一个**酉**对角因子乘一个 Chebyshev-Vandermonde 型矩阵：

$$
V=\underbrace{\mathrm{diag}\bigl(\mathrm i^0,\mathrm i^1,\dots,\mathrm i^{\,n-1}\bigr)}_{=:\,\Theta,\ \text{酉}}\;
\underbrace{\begin{bmatrix}
U_0(x_1)&\cdots&U_0(x_n)\\
\vdots&&\vdots\\
U_{n-1}(x_1)&\cdots&U_{n-1}(x_n)
\end{bmatrix}}_{=:\,\Phi},
$$

因此 $\mathrm{Cond}_2(V)=\mathrm{Cond}_2(\Theta\Phi)=\mathrm{Cond}_2(\Phi)$——酉因子是免费的，条件数问题被完全归约到 $\Phi$。最后一步用 **Christoffel-Darboux 公式**与相关正交多项式的性质估计 $\|\Phi\|_2$ 与 $\|\Phi^{-1}\|_2$，得到 $\mathrm{Cond}_2(\Phi)=\mathcal O(n^2)$。

### 定理

**（主定理：条件数的多项式增长。）** 设 $B$ 来自上述边值方法的等步长离散，$V$ 是 $\mathcal B=\Delta tB$ 的特征向量矩阵（按上式规范化，$p_{j,0}=1$），则

$$
\mathrm{Cond}_2(V)=\mathcal O(n^2) .
$$

配合 $\texttt{舍入误差}=\mathcal O(\epsilon\,\mathrm{Cond}_2(V))$，这意味着舍入误差只随时间点数**多项式**增长，而不是几何步长路线的指数增长。摘要给出的直接推论是：与其他直接时间并行算法相比，可以用**大得多的 $n$** 来获得令人满意的并行度。（本站[[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i|ParaDiag-I 一章]]转述的版本把这条界限定在 $n\ge8$。）

**（界是保守的。）** ParaDiag 综述与本文 arXiv v2 源码中一段被注释掉的说明都记录：$\mathcal O(n^2)$ 是**保守**的上界，数值上观察到的是 $\mathrm{Cond}_2(V)=\mathcal O(n^{1.75})$，作者未能证明。

**（快速谱分解算法。）** 论文另设计了一个 $\mathcal O(n^2)$ 的结构利用型算法来计算 $\mathcal B$ 的谱分解，特别是 $V^{-1}=\Phi^{-1}\Theta^{*}$，其基础是第二类 Chebyshev 多项式的三项递推 $2yU_j(y)=U_{j+1}(y)+U_{j-1}(y)$。论文报告它比 MATLAB 的 `eig` 快得多。

**（二阶问题免费覆盖。）** 二阶情形的全时间矩阵是 $B^2\otimes I_x+I_t\otimes A$，$B^2$ 与 $B$ 同一个 $V$，因此不需要新的分析。

非线性问题在论文 §2 处理（该节存在已核实），**具体处理方式本站未核实**。

把两条 ParaDiag-I 路线并排看，这一篇改动了什么就很清楚：

| 路线                      | 网格       | 可对角化的来源     | $\mathrm{Cond}_2(V)$ 的增长          | 实用 $n$     |
| ------------------------- | ---------- | ------------------ | ------------------------------------ | ------------ |
| Maday-Rønquist / 几何步长 | 变步长     | 步长互异使谱分开   | 随 $\tau\to1$ 指数爆炸               | 约 $20$–$25$ |
| 编号 59 / 边值方法        | **等步长** | 格式本身给出单重根 | $\mathcal O(n^2)$（实测 $n^{1.75}$） | 无原理性上限 |

### 数值实验

**本文自身报告的（摘要，已核实）：** 在并行机上给出了数值结果，**256 核下取得超过 60 倍的加速**；论文 §4.1 与 Gander 等的几何步长直接算法作对比，并记录了后者约 $n\approx20$ 到 $25$ 的限制。

| 项目                       | 数值 / 内容                 | 核实状态                            |
| -------------------------- | --------------------------- | ----------------------------------- |
| 并行规模                   | 256 核                      | 摘要，已核实                        |
| 报告加速比                 | 超过 $60$ 倍                | 摘要，已核实                        |
| 对比基线                   | Gander 等的几何步长直接算法 | 论文 §4.1，已核实                   |
| 基线的可用时间点数         | 约 $20$–$25$                | 论文 §1，已核实                     |
| $\mathrm{Cond}_2(V)$ 实测  | $\mathcal O(n^{1.75})$      | ParaDiag 综述作图，已核实（无证明） |
| 测试 PDE、网格、墙钟时间表 | —                           | **未核实**（未读正文）              |

$60/256\approx23\%$ 的并行效率对一个直接方法是合理的量级，但**在没有墙钟表的情况下无法把损失归因**：三步中的 (a)、(c) 是时间方向的全局变换（通信密集），而 (b) 的 $N_t$ 个复移位空间系统的移位各不相同，若用迭代法求解则迭代次数不同、天然负载不均。这两项各占多少，摘要无法回答。还有一点值得明确：**256 核是否全部用于时间方向、还是时空混合划分，本站未核实**——若是纯时间并行，则 $n\ge256$ 本身就是「上限被取消」最直接的证据；若是混合划分，这一点便无法从加速比读出。

**基线一侧的定量证据（来自 GWZ 综述，即本站[[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-4-paradiag-i|ParaDiag-I 一章]]，非本文实验）。** 编号 59 要取消的那个上限，在综述里有完整的实验刻画，正好可以作为本文动机的量化背景：

| 项目       | 设定                                                                 |
| ---------- | -------------------------------------------------------------------- |
| 方程       | 一维热方程与对流扩散方程（黏性 $10^{-2}$）                           |
| 边界与初值 | 齐次 Dirichlet，$u_0(x)=\sin(2\pi x)$                                |
| 空间网格   | $\Delta x=1/50$                                                      |
| 时间窗口   | $T=0.2$（扫描 $\varrho$）；$T=0.5$（扫描 $N_t$）                     |
| 扫描范围   | $\varrho\in[10^{-2},1]$，五组 $N_t$；以及 $N_t=2^4,2^5,\dots,2^{10}$ |
| 误差度量   | 全部时间节点上的最大 $L^\infty$ 误差                                 |

| 观察对象                            | 结果                                      |
| ----------------------------------- | ----------------------------------------- |
| 误差随 $\varrho$ 的变化             | 每条曲线都有唯一极小点                    |
| $\varrho_{\rm opt}$ 的预测能力      | 对流扩散很准；热方程在小 $N_t$ 处略有偏差 |
| 取数值最优 $\varrho$ 后误差随 $N_t$ | 先降，在不足 $100$ 步处越过阈值并迅速上升 |
| 对照：等步长后向 Euler              | 误差随 $N_t$ 单调下降                     |

最后两行是这组实验的要害：**上限是舍入上限，不是收敛上限。** 同一个时间离散在等步长下误差随 $N_t$ 一路下降，只有在几何步长加对角化之后才在约 $100$ 步处折返，可见损失完全来自对角化本身。这正是编号 59 通过更换时间离散所消除的那个机制。

**这组实验没有确立什么。** 综述的实验用的是一维小规模问题，只测精度不测墙钟；本文的实验报告墙钟加速比但具体设置未核实。因此「$\mathrm{Cond}_2(V)=\mathcal O(n^2)$ 在实际问题上真的允许 $n$ 取到几百」这一条，在本站可核实的范围内只有间接支持。

### 与其他论文的关系

这一篇是 **ParaDiag-I**（直接、Maday-Rønquist 分支）的决定性修补。[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|编号 31 与 46]] 利用的是时间**周期**问题所免费提供的对角化（$V$ 酉、$\mathrm{Cond}_2(V)=1$）；编号 59 则为真正的**初值**问题拿到 $\mathcal O(n^2)$，途径是换**时间积分器**而不是换步长。

与 **ParaDiag-II** 分支（编号 39、53、65、71、84）的对比是这条线索的主轴：后者保留等步长与标准积分器，改为用 $\alpha$-循环矩阵作**预条件子**，代价是 $\mathrm{Cond}_2(V)\le1/\alpha$ 并且需要外层迭代。编号 59 的卖点是**直接**——不迭代、无 $\alpha$，而条件数只多项式增长。

编号 65 与本篇共享「一阶与二阶演化问题在同一框架内」的抱负，且作者群重叠：编号 65 在预条件迭代分支上做，编号 59 在直接分支上做。最后，用 Chebyshev 与 Christoffel-Darboux 分析一个 Vandermonde 型特征向量矩阵，在这一组论文里没有第二例，是其中最「特殊函数」的一篇。

## 65：用「稳定」这一条经典性质换掉全部结构性假设

### 直觉

到 2022 年，块 $\alpha$-循环预条件已成为主流的时间并行方法，尤其对双曲问题，因为 $\mathcal P_\alpha^{-1}$ 的作用经 FFT 归约为各时间层上的独立求解。但它的理论是**逐情形堆积**起来的：每篇论文选定一个时间积分器，利用该积分器的特殊结构算出 $\mathcal P_\alpha^{-1}\mathcal K$ 的谱。综述对此的评价很直接：这些分析繁复，并且严重依赖时间推进矩阵的特殊性质，例如**稀疏性、Toeplitz 结构与对角占优**。

缺少的是一条覆盖全部积分器的定理，而且其假设应当是一条**经典**性质而不是某个特设的结构性质。本文给出的那条性质是**稳定性**。

为什么稳定性恰好是对的假设？可以先看清预条件子到底改了什么。$\mathcal P_\alpha$ 与全时间矩阵 $\mathcal K$ 的差别只在右上角一个 $\alpha$：它把最后一个时间层以权重 $\alpha$ 绕回第一个时间层。所以预条件所犯的错误，正是**让传播子绕整个时间窗口走一圈**。若方法稳定，走完一整圈不会放大——用标量语言说 $|r|^{N_t}\le1$——那么这个绕回扰动的大小就被 $\alpha$ 本身封顶，**与步数、网格尺度、$A$ 的刚性统统无关**。下面的推导会把这句话写成 $\lambda=1/(1-\alpha r^{N_t})$，网格无关性就是这个式子里 $r^{N_t}$ 被 $|r|\le1$ 控制住的直接后果。（这一段是本页对定理的读法，不是论文的表述。）

第二个缺口是内存。二阶（波型）问题通常通过把 $u''=Au+g$ 改写成一阶系统来处理，这会把每个时间层的内存翻倍——对精细空间网格与高维问题不利——因此需要对二阶形式的**两步**离散做直接分析。

### 问题设定

一阶情形把一般一步法写成两矩阵形式

$$
r_1(\Delta t A)\,\boldsymbol u_n=r_2(\Delta t A)\,\boldsymbol u_{n-1}+\tilde{\boldsymbol g}_n,
\qquad n=1,\dots,N_t,
$$

因此方法的稳定函数是 $r_1^{-1}(z)r_2(z)$（例如后向 Euler 取 $r_1(z)=I-z$、$r_2(z)=I$）。堆叠后得到全时间矩阵与其块 $\alpha$-循环预条件子

$$
\mathcal K=I_t\otimes r_1(\Delta tA)-B\otimes r_2(\Delta tA),
\qquad
\mathcal P_\alpha=I_t\otimes r_1(\Delta tA)-C_\alpha\otimes r_2(\Delta tA),
$$

写开就是

$$
\mathcal P_\alpha=
\begin{bmatrix}
r_1(\Delta tA)&&&-\alpha\,r_2(\Delta tA)\\
-r_2(\Delta tA)&r_1(\Delta tA)&&\\
&\ddots&\ddots&\\
&&-r_2(\Delta tA)&r_1(\Delta tA)
\end{bmatrix},
\qquad
C_\alpha=
\begin{bmatrix}0&&&\alpha\\1&0&&\\&\ddots&\ddots&\\&&1&0\end{bmatrix},
$$

其中 $B$ 是严格下移矩阵、$C_\alpha=B+\alpha e_1e_{N_t}^\top$。也就是说，**预条件子是把全时间矩阵的右上角放上 $\alpha$ 得到的**，即用 $\alpha$-循环矩阵替换（块）Toeplitz 矩阵。$\alpha=1$ 给出普通循环矩阵（Strang 型），$\alpha\to0$ 回到 $\mathcal K$ 本身。

二阶情形不做一阶化，而用对称两步法

$$
r_1(\Delta t^2A)\,\boldsymbol u_{n+1}-r_2(\Delta t^2A)\,\boldsymbol u_n
+r_1(\Delta t^2A)\,\boldsymbol u_{n-1}=\tilde{\boldsymbol g}_n,
\qquad n=1,\dots,N_t-1,
$$

例如 Numerov 型方法取 $r_1=I_x-\tfrac{\Delta t^2A}{12}+\tfrac{10\gamma(\Delta t^2A)^2}{12}$、$r_2=2I_x+\tfrac{10\Delta t^2A}{12}+\tfrac{20\gamma(\Delta t^2A)^2}{12}$。相应的全时间矩阵与预条件子为

$$
\mathcal K=\tilde B\otimes r_1(\Delta t^2A)-B\otimes r_2(\Delta t^2A),
\qquad
\mathcal P_\alpha=\tilde C_\alpha\otimes r_1(\Delta t^2A)-C_\alpha\otimes r_2(\Delta t^2A),
$$

$$
\tilde B=\begin{bmatrix}1&&&&\\0&1&&&\\1&0&1&&\\&\ddots&\ddots&\ddots&\\&&1&0&1\end{bmatrix},
\qquad
\tilde C_\alpha=\begin{bmatrix}1&&&\alpha&\\0&1&&&\alpha\\1&0&1&&\\&\ddots&\ddots&\ddots&\\&&1&0&1\end{bmatrix},
$$

其中两步模板使 $\alpha$ 出现在**两个**角落位置。

$\mathcal P_\alpha^{-1}$ 的作用方式与 $\alpha$-循环族其余成员相同：$C_\alpha=VDV^{-1}$，$V=\Gamma_\alpha^{-1}\mathbb F^*$，$\Gamma_\alpha=\mathrm{diag}(1,\alpha^{1/N_t},\dots,\alpha^{(N_t-1)/N_t})$，因此步骤 (a) 与 (c) 是（带缩放的）FFT，步骤 (b) 是 $N_t$ 个独立的复空间求解，而 $\mathrm{Cond}_2(V)\le1/\alpha$ 是小 $\alpha$ 的舍入代价。这套装置在[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|对角化技术一页]]有完整推导。

### 推导

方向可以自己算出来，不必依赖转述。逐个取 $z\in\sigma(\Delta tA)$，记 $r:=r_1^{-1}(z)r_2(z)$，稳定性即 $|r|\le1$。在该标量通道内 $\mathcal K=I_t-rB$、$\mathcal P_\alpha=I_t-rC_\alpha$，而 $C_\alpha=B+\alpha e_1e_{N_t}^{\top}$，于是

$$
\mathcal K=\mathcal P_\alpha+\alpha r\,e_1e_{N_t}^{\top}
\qquad\Longrightarrow\qquad
\mathcal P_\alpha^{-1}\mathcal K=I_t+\alpha r\,\bigl(\mathcal P_\alpha^{-1}e_1\bigr)e_{N_t}^{\top},
$$

这是单位矩阵的**秩一修正**。因此它的特征值是重数 $N_t-1$ 的 $1$，加上唯一一个 $1+\alpha r\,e_{N_t}^{\top}\mathcal P_\alpha^{-1}e_1$。用 $\mathcal P_\alpha^{-1}=\sum_{j\ge0}r^jC_\alpha^{\,j}$ 以及 $C_\alpha^{\,j}e_1=e_{1+j}$（$j<N_t$）、$C_\alpha^{N_t}=\alpha I_t$，得

$$
e_{N_t}^{\top}\mathcal P_\alpha^{-1}e_1=\sum_{m\ge0}r^{\,N_t-1+mN_t}\alpha^m=\frac{r^{\,N_t-1}}{1-\alpha r^{\,N_t}},
\qquad
\lambda=1+\frac{\alpha r^{\,N_t}}{1-\alpha r^{\,N_t}}=\frac{1}{1-\alpha r^{\,N_t}} .
$$

由 $|r|\le1$ 得 $|\alpha r^{N_t}|\le\alpha$，故 $|1-\alpha r^{N_t}|\in[1-\alpha,1+\alpha]$，正好给出 $\frac1{1+\alpha}\le|\lambda|\le\frac1{1-\alpha}$；同一个计算还给出综述自己的推论 $\rho(\mathcal M)=\bigl|\alpha r^{N_t}/(1-\alpha r^{N_t})\bigr|\le\alpha/(1-\alpha)$。两个结论出自同一步，方向由此确定。

这个计算顺带说明了两件事。其一，每个标量通道内除**一个**特征值外全部精确等于 $1$，所以整个系统最多有 $N_x$ 个特征值偏离 $1$——这正是 McDonald、Pestana 与 Wathen 在 $\alpha=1$ 情形下聚集性定理的对应物。其二，偏离量带因子 $r^{N_t}$，因此对**严格**压缩的通道（$|r|<1$，即有耗散）聚集性远好于最坏情形的界。**这解释了为什么该方法在抛物问题上的表现明显优于双曲问题**：双曲情形 $|r|\approx1$，$r^{N_t}$ 不衰减，最坏界几乎是紧的。

> [!note] 关于这段计算的性质
> 以上是本页为确定界的方向而做的推导，只覆盖标量通道，并非论文本身的证明。论文发表于 _SIAM J. Matrix Anal. Appl._，其定理可能以更精细的、依赖 $z$ 的端点形式给出；若要引用论文原始表述，仍应核对期刊正文。

### 定理

**（主定理，一阶情形。）** 设 $\mathcal K$ 来自上述一步积分器，$\mathcal P_\alpha$ 是 $\alpha\in(0,1)$ 的块 $\alpha$-循环预条件子。若该积分器**稳定**，即

$$
\bigl|r_1^{-1}(z)\,r_2(z)\bigr|\le1
\qquad \forall z\in\sigma(\Delta t A)\subset\mathbb C^-,
$$

则预条件矩阵的每个特征值满足**与网格无关**的界

$$
\frac{1}{1+\alpha}\ \le\ \bigl|\lambda(\mathcal P_\alpha^{-1}\mathcal K)\bigr|\ \le\ \frac{1}{1-\alpha} .
$$

> [!warning] 界的方向
> 综述中转述这条定理时把两个端点印反了（写成 $\frac{1}{1-\alpha}\le|\lambda|\le\frac{1}{1+\alpha}$），这对 $\alpha\in(0,1)$ 不可能成立，因为 $1/(1-\alpha)>1/(1+\alpha)$。上式给出的是正确方向，理由见上面的直接计算。

**（主定理，二阶情形。）** 对上述对称两步法，若

$$
\bigl|r_1^{-1}(z)r_2(z)\bigr|\le2
\qquad\forall z\in\sigma(\Delta t^2A)\subset\mathbb R^-,\ \text{且仅在 }z=0\text{ 取等},
$$

则同样的两条界成立。其中的「$\le2$」是对称两步格式的自然稳定性条件，与中间项 $-r_2\boldsymbol u_n$ 的系数相匹配。

**（推论：定常迭代的收缩。）** 对迭代矩阵 $\mathcal M=\mathcal I-\mathcal P_\alpha^{-1}\mathcal K$，

$$
\rho(\mathcal M)\le\frac{\alpha}{1-\alpha},
$$

即一个 $\mathcal O(\alpha)$ 的收缩，**与 $\Delta t$、$\Delta x$、$N_t$、$T$ 及 $\sigma(A)$ 全都无关**。这恰好就是[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|编号 53]] 所猜测的收敛率，现在对全部稳定一步法得证。

把这三条写成数字，$\alpha$ 该怎么取一目了然：

| $\alpha$  | $1/(1+\alpha)$ | $1/(1-\alpha)$ | $\rho(\mathcal M)\le\alpha/(1-\alpha)$ | $\mathrm{Cond}_2(V)\le1/\alpha$ |
| --------- | -------------- | -------------- | -------------------------------------- | ------------------------------- |
| $10^{-1}$ | $0.9091$       | $1.1111$       | $0.1111$                               | $10$                            |
| $10^{-2}$ | $0.9901$       | $1.0101$       | $0.0101$                               | $10^{2}$                        |
| $10^{-3}$ | $0.9990$       | $1.0010$       | $0.0010$                               | $10^{3}$                        |

（前四列由定理与推论直接算出；最后一列是 $\alpha$-循环对角化的通用界。）再叠上推导给出的第二条信息——每个标量通道只有**一个**特征值偏离 $1$——就可以看出这条定理对 Krylov 方法有多友好：整个 $N_xN_t$ 维系统里至多 $N_x$ 个特征值不等于 $1$，其余全部精确聚在 $1$。

**（假设的锐性。）** 稳定性是**充分**条件；综述报告数值上它也是**必要**的，示例是 Numerov 型方法在参数 $\gamma=1/120$（无条件稳定、四阶）与 $\gamma=1/120.01$（不稳定）之间的对比，后者的 $\sigma(\mathcal M)$ 逃出半径 $\alpha/(1-\alpha)$ 的圆。

**（适用范围与其边界。）** 摘要的措辞是：一阶问题的分析适用于**全部稳定单步积分器**，二阶问题的分析适用于**一大类可任意高阶的对称两步方法**。但对**一般多步**方法（超出一步与对称两步之外）界不必成立：例如 Volterra 偏积分微分方程给出的 $B$ 是**稠密**下三角 Toeplitz 矩阵，此时只能得到 $|\lambda(\mathcal P_\alpha^{-1}\mathcal K)|=1+\mathcal O(\alpha)$，且还需要求积权满足正性或单调性条件。**这条边界正是编号 84 的出发点。**

> [!note] 未核实项
> SIMAX 正文中的定理编号、以及**奇异值**是否与特征值一样被夹住，本站未核实。

### 数值实验

摘要只确认「给出了示例性数值实验以配合理论」。**具体 PDE、网格参数与迭代次数未核实。** 但综述在两处引用了本文的实验，这两处都值得记录，因为它们说明的是理论**没有**覆盖的情形：

| 引用的实验                   | 观察到的现象                                                                                    |
| ---------------------------- | ----------------------------------------------------------------------------------------------- |
| 非线性问题，Newton + GMRES   | $\mathcal P_\alpha^{-1}\mathcal J$ 的谱聚集性；可能出现 $\rho>1$，定常迭代失效而 GMRES 仍然好用 |
| 时间窗口 $T$ 的影响          | 较短的 $T$ 给出更聚集的谱与更快的 GMRES 收敛                                                    |
| Numerov 参数（综述自身实验） | $\gamma=1/120$ 落在半径 $\alpha/(1-\alpha)$ 的圆内；$\gamma=1/120.01$ 逃出                      |

第一行的意义在于把定理与实践分开：定理保证的是 $\rho(\mathcal M)\le\alpha/(1-\alpha)$，这是**定常迭代**的收缩；一旦把 $\mathcal K$ 换成非线性问题的 Jacobi 矩阵 $\mathcal J$，定理的假设不再成立，$\rho$ 可以超过 $1$，但 GMRES 只需要谱**聚集**而不需要 $\rho<1$，所以仍然收敛得快。这是选 Krylov 而不是定常迭代的具体理由。

第二行需要小心解读，因为它与线性理论看似相反。线性情形的界 $\frac1{1+\alpha}\le|\lambda|\le\frac1{1-\alpha}$ 完全不含 $T$；而 Gander 与 Wu 对后向 Euler 的加细界——按 $\boldsymbol u'+A\boldsymbol u=\boldsymbol g$ 的符号约定（与本节相反），设 $\Re\lambda(A)\ge\kappa\ge0$——是

$$
\rho\le\frac{\alpha e^{-T\kappa}}{1-\alpha e^{-T\kappa}},
$$

它**随 $T$ 增大而变好**。所以「短 $T$ 更好」不可能是线性图景里的现象，只能来自非线性设定：$\mathcal P_\alpha$ 是用某个冻结的 Jacobi 矩阵造出来的，窗口越短，Jacobi 矩阵在窗口内变化越小，冻结近似越准。**这一解释是本页的读法，不是综述或论文的表述**，但两条界的 $T$-依赖是已核实的，因此至少可以确定：那条观察不能用线性理论去解释。

**这些实验没有确立什么。** 全部三行都是定性的：没有一处给出迭代次数、网格规模或墙钟时间。本文的贡献本来就是一条定理而不是一个新算法，所以实验的作用是佐证而非度量；但这也意味着「$\alpha=10^{-2}$ 在真实规模问题上需要几步 GMRES」这个使用者最关心的数字，在本站可核实的范围内取不到。

### 与其他论文的关系

编号 53 提出了猜想——「隐式 Runge-Kutta 方法的 A-稳定性足以保证 $\mathcal O(\alpha)$ 的稳健收敛率」——并对一族方法（两级 SDIRK，条件 $\gamma\ge1/4$）作了验证；编号 65 对全部稳定一步法证明了对应的命题。两者合起来是**个案与定理**的关系。

[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|编号 39 与 31]] 建造了 $\alpha$-循环机制，编号 65 是它的谱理论。编号 46 在**前向-后向**（KKT）设定下对两个具体积分器证明了特征值与奇异值的聚集性；本页的编号 71 则是编号 65 的前向-后向对应物——同样的「一条定理覆盖全部稳定一步法」的抱负，用在耦合系统上，代价是一条 $\alpha$ 随 $N_t$ 的定标律。编号 59 与本篇共享「一阶与二阶统一处理」的范围，但在直接分支上：编号 59 界的是精确对角化的 $\mathrm{Cond}_2(V)=\mathcal O(n^2)$，编号 65 界的是廉价近似所得预条件系统的谱。编号 84 把分析推到时间谱离散，那里时间块稠密、$\mathcal K$ 根本不是块二对角的——正是本篇假设之外的情形之一。

## 71：前向-后向情形有多个 Toeplitz 块

### 直觉

前面所有方法处理的都是**一个**方向的演化。最优控制与反问题给出的却是一对方向相反且互相耦合的演化：状态方程从 $t=0$ 向前，伴随方程从 $t=T$ 向后，而两者通过控制变量绑在一起。此时连「串行算法」都不存在——没有哪个方向可以单独推进，必须一次性求解耦合系统。**时间并行在这里不是加速手段，而是唯一的组织方式。**

机制上的想法与编号 65 完全一样：给时间方向加一个权重为 $\alpha$ 的绕回，把 Toeplitz 换成 $\alpha$-循环，从而可 FFT 对角化。差别在于绕回要做**几次**：耦合系统里出现的不是一个 Toeplitz 矩阵，而是若干个（前向矩阵及其转置至少各一个），每一个都要循环化。摘要的用词是复数的「那些 Toeplitz 矩阵」，这个细节本身就说明了结构。

这样得到的东西并不是全新的：论文作者自己在综述中把它认定为 Pearson、Stoll 与 Wathen 的**匹配 Schur 补**预条件子的**并行版本**。也就是说，这条路线的贡献不是发明一个预条件子，而是把一个已经成熟的串行预条件子改造成可对角化、从而可时间并行的形式。

代价是这一篇最值得记住的技术区别：**$\alpha$ 不再能自由固定，必须随 $N_t$ 缩小。** 编号 65 对单个前向演化给出的界对任意固定 $\alpha\in(0,1)$ 都成立；前向-后向耦合夺走了这个自由。

### 问题设定

论文的两类应用是 PDE 约束最优控制问题（抛物 KKT 系统）与抛物源识别问题（从观测恢复未知源，同样给出前向-后向最优性系统）。两者的连续形式是同一类：

$$
\begin{cases}
\partial_t y-\Delta y=f+\tfrac1\gamma p, & y(\cdot,0)=y_0 \quad(\text{前向}),\\[4pt]
-\partial_t p-\Delta p=y_d-y, & p(\cdot,T)=0 \quad(\text{后向}).
\end{cases}
$$

> [!note] 这组方程是重建的
> 上式是这类问题的标准形式，用于说明结构；**论文自身的记号与具体形式本站未核实**。

摘要（散文部分）确认的是：主要计算负担是求解一个大型线性系统，其中心对象是**前向子问题**经时空离散后的全时间矩阵，因此高效求解器需要一个好的预条件子。前向-后向情形特有的困难是：两个演化耦合且方向相反，都不能独立时间推进，得到的是鞍点型或 Schur 补型系统而不是单个块下三角 Toeplitz 系统。编号 65 的 $\alpha$-循环理论只覆盖单一前向演化。

### 推导

**预条件子的构造（摘要散文已核实，符号被删）。** 「受 $[\mathcal K]$ 的结构启发，我们用 $[\mathcal P_\alpha]$ 作预条件，其中 $[\mathcal P_\alpha]$ 是把 $[\mathcal K]$ 中的**那些 Toeplitz 矩阵**替换为 $\alpha$-循环矩阵所得的块 $\alpha$-循环矩阵。」

**重建的显式形式。** 记前向全时间矩阵与其 $\alpha$-循环替身为

$$
\mathcal K=I_t\otimes r_1(\Delta tA)-B\otimes r_2(\Delta tA),
\qquad
\mathcal K_\alpha=I_t\otimes r_1(\Delta tA)-C_\alpha\otimes r_2(\Delta tA),
\qquad
C_\alpha=B+\alpha e_1e_{N_t}^\top,
$$

耦合系统的算子同时含 $\mathcal K$ 与 $\mathcal K^\top$，预条件子则把两处都换成 $\mathcal K_\alpha$ 与 $\mathcal K_\alpha^\top$。**这一显式排布是按编号 65 与 ParaDiag 综述中波动方程控制问题的模板重建的，论文的实际写法未核实。**

**与匹配 Schur 补的关系（这一条已核实）。** 综述明确把本文的构造描述为 Pearson-Stoll-Wathen（2012）匹配 Schur 补预条件子的并行版本。记 $\mathcal M_h$ 为时空质量矩阵、$\gamma$ 为正则化参数，MSC 的想法是把精确 Schur 补

$$
S=\mathcal K\mathcal M_h^{-1}\mathcal K^\top+\tfrac1\gamma\mathcal M_h
$$

近似成一个**配对乘积**

$$
\widehat S=\Bigl(\mathcal K+\tfrac{1}{\sqrt\gamma}\mathcal M_h\Bigr)\mathcal M_h^{-1}\Bigl(\mathcal K+\tfrac{1}{\sqrt\gamma}\mathcal M_h\Bigr)^\top,
$$

展开后首项与末项与 $S$ 完全一致，而且 $\widehat S$ 对称正定。把其中的 $\mathcal K$ 换成 $\mathcal K_\alpha$，整个乘积就变成可 FFT 对角化的，于是每一次预条件作用都能在全部时间步上并行。**$\widehat S$ 的公式是 Pearson-Stoll-Wathen 的标准构造，不是从本文读到的**；已核实的是「本文的预条件子是 MSC 的并行版本」这一认定。$\widehat S$ 的对称正定性正是外层可以用共轭梯度法的原因，这一点与实验部分互相印证。

**作用方式（摘要已核实）。** 「通过对 $[\mathcal P_\alpha]$ 作块 Fourier 对角化，预条件步骤的计算对全部时间步是可并行的。」即通常的 FFT → 解耦复空间求解 → 逆 FFT 三步，$V=\Gamma_\alpha^{-1}\mathbb F^*$，$\mathrm{Cond}_2(V)\le1/\alpha$。

### 定理

**（主定理，散文形式已核实，公式未核实。）** 「我们给出预条件矩阵 $[\mathcal P_\alpha^{-1}\mathcal K]$ 的谱分析，并证明**对任意一步稳定时间积分器**，$[\mathcal P_\alpha^{-1}\mathcal K]$ 的特征值**分布在一个与网格无关的区间内**，只要参数 $[\alpha]$ **随时间步数弱定标为** $[\text{公式被删}]$，其中 $[C]$ 是一个自由常数。」

定理的**结构**完全可核实：假设是「任意一步稳定时间积分器」，与编号 65 相同；结论是「与网格无关的特征值区间」；而与编号 65 关键的不同在于，$\alpha$ 必须**随 $N_t$ 定标**而不能自由固定。本站把这条定标律记为

$$
\alpha=C\,N_t^{-\theta},\qquad \theta>0,\ C\ \text{为自由常数},
$$

**其中指数 $\theta$ 未核实，区间端点同样未核实。**

> [!warning] 可核实范围
> 该文摘要在所有公开来源中（OpenAlex 倒排索引、Crossref 衍生聚合站、经检索中继获取的出版社页面）**内联数学符号均被删除**，产生形如「求解线性系统 ___ 通常是主要计算负担……其中 ___ 是所谓的全时间矩阵」这样的句子。散文部分可靠，符号无法从公开来源恢复。检索合成曾把该公式渲染为 $\alpha=C/N_t$，但那是对一个符号已被剥离的页面所作的转述，不构成来源；而「弱定标」这个用词本身也可能指向更慢的衰减（如 $N_t^{-1/2}$）。**本页因此不给出 $\theta$ 的任何数值。**

不过定标律的**形式**已经足以推出一个后果。把它代入 $\alpha$-循环对角化的通用界 $\mathrm{Cond}_2(V)\le1/\alpha$：

$$
\mathrm{Cond}_2(V)\ \lesssim\ \frac{N_t^{\theta}}{C}
\qquad\Longrightarrow\qquad
\texttt{舍入地板}\ \approx\ \epsilon\,\frac{N_t^{\theta}}{C}.
$$

也就是说，**耦合把编号 65 里那个与 $N_t$ 无关的舍入地板变成了随 $N_t$ 多项式增长的地板。** 这个退化比 ParaDiag-I 的 $\varrho^{-(N_t-1)}$ 温和得多（多项式对指数），但它不再像编号 65 那样完全不存在。$\theta$ 的具体数值决定了这个方法在多长的时间窗口上仍然可用，而这恰好就是公开来源删掉的那个符号。（此推论由两条已核实的事实相乘得到，指数本身仍未核实。）

论文的其余可核实结论：两类应用（PDE 约束最优控制、抛物源识别）各自被完整做出；构造被作者自己认定为并行 MSC 预条件子，因而属于 Pearson-Stoll-Wathen 一脉而非孤立发明。**定理编号、奇异值是否同样受控、以及对正则化参数 $\gamma$ 的依赖，本站均未核实。**

### 数值实验

**摘要已核实：** 两类应用的数值结果都表明，谱分析**很好地预测了预条件共轭梯度法的收敛率**。这句话里有两处值得单独拿出来。

其一，外层求解器是 **CG**，不是 GMRES 或 MINRES。这意味着被迭代的预条件系统是**对称正定**的，与上面 MSC/Schur 补的读法一致，也是与相邻工作的一处实质差别：

| 论文    | 外层求解器       | 隐含的谱性质 |
| ------- | ---------------- | ------------ |
| 编号 46 | GMRES / BiCGStab | 非对称       |
| 编号 65 | GMRES            | 非对称       |
| 编号 71 | **CG**           | **对称正定** |

其二，论文的主张是预测区间与实测 CG 收敛率之间的**定量吻合**，而不只是「比不预条件快」。这是一个比加速比更强的主张，因为它把定理里的区间端点直接对上了实验曲线的斜率。

**未核实：** 网格规模、$\gamma$ 的取值、$N_t$、所选的自由常数 $C$、以及迭代次数，全部无法从公开来源取得。

**这留下的缺口是实用性的，不只是学术性的。** 使用这个方法必须先决定 $\alpha$，而定理说 $\alpha$ 要随 $N_t$ 定标——不知道 $\theta$ 与 $C$，就没有配方。编号 65 没有这个问题（任意固定的 $\alpha$ 都行，实践取 $10^{-2}$ 到 $10^{-3}$），编号 71 恰恰在这一点上把参数选择变成了必须查阅正文才能完成的一步。

### 与其他论文的关系

[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|编号 46]] 是直接前身：同一类物理问题（抛物 PDE 约束最优化，「演化方向相反的耦合 PDE 系统」），同一个策略（把不可对角化的时间矩阵换成循环型替身并用作预条件子）。编号 46 对**两个**具体积分器（后向 Euler 与梯形规则）证明了特征值与奇异值的聚集性并用 GMRES/BiCGStab；编号 71 对**全部**一步稳定积分器证明了与网格无关的区间并改用 CG。**编号 71 之于编号 46，正如编号 65 之于编号 53：用一条一般定理取代若干个案。**

编号 65 是解耦前向问题上的姊妹结论，两者的对照是这一对论文最锋利的一句话：单个前向演化时，$|\lambda(\mathcal P_\alpha^{-1}\mathcal K)|\in[\tfrac{1}{1+\alpha},\tfrac{1}{1-\alpha}]$ 对**任意固定** $\alpha\in(0,1)$ 成立；前向-后向系统则要求 $\alpha$ **随 $N_t$ 收缩**。耦合的代价就是失去固定 $\alpha$ 的自由。编号 53 贡献了「积分器稳定即足以给出 $\mathcal O(\alpha)$ 收敛率」的猜想，编号 71 在耦合设定中确认了同一模式（「对任意一步稳定时间积分器」）；编号 31 与 39 提供底层的 $\alpha$-循环与 FFT 机制；编号 84 是同一纲领的最新成员，走向稠密的时间谱块。

## 84：没有 Toeplitz 结构可循环化时怎么办

### 直觉

时间谱方法用基函数（例如多项式）的组合逼近解，是空间谱离散的自然搭配，时间方向精度很高。但全部组合系数必须**一次性**通过求解一个全时间系统得到，完全没有时间推进的选项——这一点与编号 31 的时间周期问题结构相同：离散本身就不提供串行算法，全时间求解是必需而非可选。

问题出在这条路线赖以运作的那句配方上。ParaDiag-II 一族此前的每篇论文都依赖时间矩阵是（块）**Toeplitz** 的——通常还是块二对角 Toeplitz——从而预条件子可以由一句话得到：「把 Toeplitz 矩阵换成它的 $\alpha$-循环对应物」。时间谱方法给出的矩阵 $M$ 稠密且不带 Toeplitz 结构，**因此没有可循环化的对象，配方直接失效**。

本文的解法可以用五个字概括：**先因式分解，再循环化。** 把稠密无结构的 $M$ 先分解成若干因子，使其中的一部分**确实是** Toeplitz 的，然后只对这些 Toeplitz 因子作 $\alpha$-循环替换。$M$ 本身不需要有结构，只要它能被分解出结构来。

### 问题设定

全时间系统仍是 Kronecker 张量形式

$$
\bigl(M\otimes I_x+I_t\otimes A\bigr)\boldsymbol u=\boldsymbol b,
\qquad
M\in\mathbb R^{N_t\times N_t}\ \text{稠密、无结构},
$$

（这一显式写法是按摘要的描述重建的，论文的实际记号未核实）。与前面各篇对照，这里唯一变化的就是 $M$：一步法给出块二对角 Toeplitz 的 $B$，对称两步法给出带两条对角线的 $\tilde B$，而时间谱给出一个满矩阵。

时间矩阵为什么会变稠密，在配置型方法上看得最清楚。在节点 $t_0<t_1<\dots<t_N$ 上做积分延迟校正（IDC），其更新式里含一项

$$
\int_{t_m}^{t_{m+1}}f\bigl(\boldsymbol u^k(\tau),\tau\bigr)\,d\tau
\approx\sum_{j=1}^{N}\omega_{m,j}\,f(\boldsymbol u_j^k,t_j),
\qquad
\omega_{m,j}=\int_{t_m}^{t_{m+1}}\Bigl(\prod_{i\ne j}\frac{\tau-t_i}{t_j-t_i}\Bigr)d\tau,
$$

即用穿过**全部**节点的 Lagrange 插值多项式来积分。于是矩阵 $\Omega=[\omega_{m,j}]$ 是满的：第 $m$ 个子区间上的更新依赖**每一个**节点，包括它右边的节点，下三角结构被彻底破坏。时间谱方法是同一现象的极端形式——基函数在整个区间上有支集，所以所有系数从一开始就完全耦合。

### 推导

**核心构造（摘要已核实）。** 预条件子「由对 $M$ 的一个**新颖因式分解**、再把该分解中的 **Toeplitz 矩阵**替换为相应的 **Strang 型 $\alpha$-循环矩阵**而得到」。

两个词需要解释。**Strang 型**指经典的 Strang 循环预条件子：从 Toeplitz 矩阵的中央若干条对角线构造一个循环矩阵；$\alpha$-循环推广把绕回位置的 $1$ 换成 $\alpha$。这一支的来源是 Strang（1986），此后 $\sigma(C^{-1}B)$ 的谱分析已有三十年文献（Chan-Ng 1996；Ng 2004；Bini-Latouche-Meini 2005）。**新颖之处**在于：ParaDiag-II 里的块 $r_1(\Delta tA)$、$r_2(\Delta tA)$ 本身**不是** Toeplitz 的，所以经典的块 Toeplitz 理论不能直接搬用——这正是编号 65、71、84 三篇要填的缺口，而编号 84 面对的是其中最极端的一种：连时间方向的 Toeplitz 结构都没有了。

**范围（摘要已核实）。** 论文对 **Legendre 对偶 Petrov-Galerkin 方法**具体给出了这个因式分解，并声称对其他常用时间谱方法同样成立。参考文献确认了相关背景：Shen 的对偶 Petrov-Galerkin 方法（SINUM 2003）、Shen-Wang 对双曲方程的 Legendre 与 Chebyshev 对偶 Petrov-Galerkin 方法、Kong-Shen-Wang-Xiang 关于 Legendre 对偶 Petrov-Galerkin 初值问题的特征值分析（Adv. Comput. Math. 2024）、Tang-Ma 的单区间与多区间 Legendre $\tau$-方法，以及 Yang-Wang 的 Chebyshev-Gauss 谱配置方法。

**设计判据（摘要已核实）。** 预条件子「允许**良态的对角化**，因此每一步预条件都能以高效的时间并行方式求解」。注意「良态对角化」这个判据与编号 59 完全相同，只是这里施加在**预条件子**上而不是直接求解器上：编号 59 要求精确对角化本身条件良好，编号 84 要求近似对角化条件良好，而近似的误差交给 GMRES 去吸收。

> [!note] 未核实的关键细节
> **因式分解的显式形式本站未核实**：分成几个因子、各是什么、非 Toeplitz 的余项是对角的、三角的还是低秩的，都无法从摘要判断。这是本篇最实质的技术内容，也是本页最大的缺口。

### 定理

摘要可确认三条：（i）预条件矩阵的谱分析显示特征值**高度聚集**，从而促进 **GMRES** 的快速收敛；（ii）该因式分解对 Legendre 对偶 Petrov-Galerkin 方法具体给出，并声称对其他常用时间谱方法成立；（iii）预条件子允许良态对角化，故 $\mathcal P^{-1}$ 的作用在全部时间层上并行。

> [!warning] 不要替它补一个收敛因子
> 摘要中**没有**给出任何收敛因子或谱区间。聚集半径、它对 $\alpha$、$N_t$ 与多项式次数的依赖、以及定理编号，本站均未核实。特别地，不能把编号 65 的 $[\tfrac{1}{1+\alpha},\tfrac{1}{1-\alpha}]$ 挪用到这里——编号 65 的假设（块二对角 Toeplitz 的 $\mathcal K$）在本文的设定下恰好不成立，这正是本文存在的理由。

### 数值实验

**部分可核实。** Springer 页面显示文章含**六幅图**（Fig. 1–Fig. 6），与一个分量不轻的数值部分相符，外层求解器为 **GMRES**；数据可用性声明为「数据集可向通讯作者索取」。**具体测试问题、谱方法次数、网格规模与迭代次数未核实。**

这个缺口的位置很不巧。本文最核心的量化主张是「高度聚集」，而这是一个形容词。编号 65 的同类主张是一条带端点的区间，因此可以直接换算成 GMRES 的迭代次数上界；本文的主张若没有聚集半径，就无法与之比较，也无法回答「稠密的 $M$ 被因式分解之后，付出的代价相对于块二对角情形是多少」这个自然问题。那六幅图很可能正是回答它的地方。

### 与其他论文的关系

编号 65 是被推广的理论骨干，也被本文引用。编号 65 覆盖全部稳定**一步**法（$\mathcal K$ 为块二对角 Toeplitz）与对称**两步**法；综述明确指出，对一般多步与非 Toeplitz 的时间结构，编号 65 的界不必成立。时间谱的 $M$ 是那个缺口的极端情形：完全稠密，毫无 Toeplitz 结构。

方法论上最接近的先例是[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|编号 53]]：它早一步遇到同一类困难——**多级**积分器的差分方程无法直接堆成 Toeplitz 全时间系统——解决办法是构造一个「结构完全不同、实现细节也不同」的 $\alpha$-循环预条件子。编号 84 对**稠密**时间矩阵重复了这个动作，但走的是因式分解而不是重新推导全时间形式。编号 31 确立的模式在这里再次出现：一种**强制**全时间求解的离散（那里是时间周期，这里是谱系数）不是障碍而是时间并行的机会。编号 59 共享「良态对角化」这一判据并被引用，区别是直接求解器（编号 59，$\mathrm{Cond}_2(V)=\mathcal O(n^2)$）与预条件 GMRES（编号 84）。编号 71 是同组紧邻的前一篇，也被引用。

按时间顺序这是列表中最新的一篇，读起来像整个纲领的当前前沿：一步、两步、多级与前向-后向的情形都已覆盖之后，剩下的边疆就是无结构、稠密的时间算子。

## 85：把上述工作放进一个统一框架

编号 85 是 2025 年的 _Acta Numerica_ 综述，把上述工作与更广的文献统一组织为两类：对**传播**型问题仍然有效的方法（Schwarz 波形松弛、积分延迟校正、ParaExp、ParaDiag），与主要为**耗散**问题设计的方法（parareal、PFASST、MGRIT、对角化 parareal、时空多重网格）。这个二分法本身就是一条判断：本页的 $\alpha$-循环路线之所以被归入第一类，是因为它的收缩率 $\alpha/(1-\alpha)$ 不依赖 $\sigma(A)$，因而不像 parareal 那样在双曲问题上失效——尽管由编号 65 的推导可知，它在双曲情形只是**不失效**，聚集性仍然明显弱于抛物情形。

本站对它的逐节精读——含原论文全部 48 个图表资产与可复现的 Python 实验——是一个独立专题：[[computational-mathematics/knowledge-notes/time-parallelization/index|双曲与抛物问题的时间并行方法]]。

## 五篇的位置关系

| 编号 | 被解决的障碍                                 | 手段                                                                |
| ---- | -------------------------------------------- | ------------------------------------------------------------------- |
| 59   | 等步长下 $B$ 不可对角化，几何步长下 $V$ 病态 | 换用边值方法，使 $\mathrm{Cond}_2(V)=\mathcal O(n^2)$               |
| 65   | 逐情形的谱分析依赖特设结构假设               | 只假设稳定性，得到 $[\frac{1}{1+\alpha},\frac{1}{1-\alpha}]$ 的模界 |
| 71   | 前向-后向系统不是单个 Toeplitz 系统          | 循环化系统中的多个 Toeplitz 块                                      |
| 84   | 时间谱矩阵没有 Toeplitz 结构                 | 先因式分解出 Toeplitz 因子，再循环化                                |
| 85   | 文献分散                                     | 按动力学类型统一组织                                                |

一条贯穿的判断：**这条路线的每一步进展都表现为「把一个假设换成一个更弱或更经典的假设」。** 编号 59 把「步长互异」换成「等步长但换格式」；编号 65 把「Toeplitz 加对角占优」换成「稳定」；编号 71 把「单个前向演化」换成「前向-后向耦合」，代价是失去固定 $\alpha$ 的自由；编号 84 指出剩下的那个假设（Toeplitz 结构）在什么时候不再成立，并给出绕过它的办法。

把五篇按「对角化要付多少舍入代价」排成一列，整条线索的经济学就清楚了：

| 路线                       | 时间矩阵                               | 对角化的舍入放大                                 | 对 $N_t$ 的依赖           | 代表           |
| -------------------------- | -------------------------------------- | ------------------------------------------------ | ------------------------- | -------------- |
| ParaDiag-I，几何步长       | 变步长后向 Euler                       | 界含因子 $\varrho^{-(N_t-1)}$                    | **指数**                  | Maday-Rønquist |
| ParaDiag-I，边值方法       | 等步长 BVM                             | $\mathcal O(n^2)$（实测 $\mathcal O(n^{1.75})$） | 多项式                    | 编号 59        |
| ParaDiag-II，固定 $\alpha$ | $\alpha$-循环                          | $\le1/\alpha$                                    | **无**                    | 编号 65        |
| ParaDiag-II，前向-后向     | $\alpha$-循环，$\alpha=CN_t^{-\theta}$ | $\lesssim N_t^{\theta}/C$                        | 多项式（$\theta$ 未核实） | 编号 71        |
| ParaDiag-II，时间谱        | 因式分解后循环化                       | 「良态」，具体界未核实                           | 未核实                    | 编号 84        |

第三行是这张表的中心：**固定 $\alpha$ 的 $\alpha$-循环预条件是这一族里唯一一个舍入代价完全不随时间窗口增长的方案**，代价是它只是近似、需要外层 Krylov 迭代。其余四行都在解释在什么情况下这个理想状态无法达到，以及退化得有多快。

## 覆盖核对

| 内容                                 | 论文 | 覆盖状态                                                                  |
| ------------------------------------ | ---- | ------------------------------------------------------------------------- |
| 三步分解与 $B$ 不可对角化的原因      | 59   | 分解式、Jordan 块、下三角 Toeplitz 的亏损性                               |
| 几何步长的两难与 $n\approx20$ 天花板 | 59   | 闭形式 $V$、比值/增量核对、截断与舍入两条界、$\varrho_{\rm opt}$          |
| 边值方法及其一致二阶精度             | 59   | 混合格式、不可作为时间推进、$B$ 与 $\boldsymbol b$、二阶情形用 $B^2$      |
| Chebyshev 机制与条件数证明           | 59   | 特征方程、根结构、$y$-代换、$V=\Theta\Phi$、Christoffel-Darboux           |
| 加速比与基线上限的实验               | 59   | 256 核 60 倍（已核实）、基线扫描表（来自综述）、未核实项已标              |
| 一阶与二阶设定的统一两矩阵写法       | 65   | $r_1,r_2$、$\mathcal K$ 与 $\mathcal P_\alpha$、两个角落的 $\alpha$       |
| 主定理与其假设的替换                 | 65   | 稳定性条件、模界、方向核对、二阶情形的 $\le2$、$\rho\le\alpha/(1-\alpha)$ |
| 秩一计算及其两个副产品               | 65   | 特征值 $1/(1-\alpha r^{N_t})$、至多 $N_x$ 个偏离、抛物优于双曲            |
| 假设的锐性与失效边界                 | 65   | Numerov $\gamma$ 对比、一般多步只得 $1+\mathcal O(\alpha)$                |
| 非线性与 $T$-依赖的实验读法          | 65   | 三行实验表；$T$-依赖与线性理论的张力已辨明（本页读法）                    |
| 前向-后向系统的多块循环化            | 71   | 障碍、配方的复数形式、MSC 及其对称正定性、CG 的理由                       |
| $\alpha$ 定标律及其后果              | 71   | 形式 $\alpha=CN_t^{-\theta}$、$\theta$ 未核实、舍入地板随 $N_t$ 增长      |
| 时间谱矩阵的非结构性                 | 84   | 一次性求解、稠密机制（Lagrange 求积）、配方失效的精确原因                 |
| 先因式分解再循环化                   | 84   | Strang 型 $\alpha$-循环、Legendre 对偶 Petrov-Galerkin、分解形式未核实    |

## 本页原文

- J. Liu, X.-S. Wang, S.-L. Wu, and T. Zhou, [_A well-conditioned direct PinT algorithm for first- and second-order evolutionary equations_](https://doi.org/10.1007/s10444-022-09928-4), Adv. Comput. Math. 48 (2022), 16（预印本 [arXiv:2108.01716](https://arxiv.org/abs/2108.01716)）。
- S.-L. Wu, T. Zhou, and Z. Zhou, [_A uniform spectral analysis for a preconditioned all-at-once system from first-order and second-order evolutionary problems_](https://doi.org/10.1137/21M145358X), SIAM J. Matrix Anal. Appl. 43(3) (2022), pp. 1331-1353。
- S.-L. Wu, Z. Wang, and T. Zhou, [_PinT preconditioner for forward-backward evolutionary equations_](https://doi.org/10.1137/22M1516476), SIAM J. Matrix Anal. Appl. 44(4) (2023), pp. 1771-1798。
- C. Tang, S.-L. Wu, T. Zhou, and Y. Zhou, [_Parallel-in-time preconditioner for the time spectral methods_](https://doi.org/10.1007/s10915-025-02899-w), J. Sci. Comput. 103 (2025), 82。
- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time parallelization for hyperbolic and parabolic problems_](https://doi.org/10.1017/S0962492924000072), Acta Numer. 34 (2025), pp. 385-489。
