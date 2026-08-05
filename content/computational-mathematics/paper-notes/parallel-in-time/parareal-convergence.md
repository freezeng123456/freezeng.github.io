---
title: parareal 的收敛分析
description: 编号 12、20、30、77：收缩因子对细传播子稳定性与细网格结构的依赖，含可核实的数值证据
lang: zh
translation: en/computational-mathematics/paper-notes/parallel-in-time/parareal-convergence
tags:
  - 论文笔记
  - 时间并行
  - 收敛分析
---

> [!note] 本页覆盖
> 编号 **12**（_SIAM J. Sci. Comput._ 37(2), 2015）、**20**（_J. Comput. Phys._ 329, 2017）、**30**（_J. Comput. Phys._ 358, 2018）、**77**（_SIAM J. Numer. Anal._ 62(5), 2024）。四篇均无预印本。编号 12 的收敛装置与门槛结论已逐式核对，其数值算例只能确认到应用类别；编号 30 的问题设定与算法思想可从摘要确认，定理常数不可得；编号 20 与 77 的摘要在所有公开渠道均未获取，本页不报告其定理、常数或实验。

## 收敛分析的标准装置

### 迭代与它读作什么

考虑由抛物 PDE 半离散得到的对称正定系统 $\boldsymbol u'(t)+A\boldsymbol u(t)=g(t)$。取 $0=T_0<T_1<\dots<T_{N_t}=T$、$T_n=n\Delta T$，parareal 迭代为

$$
\boldsymbol u_{n+1}^{k+1}=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^{k})
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k}),
$$

每个粗区间含 $J=\Delta T/\Delta t\ge2$ 个细步。昂贵的 $\mathcal F$ 只用第 $k$ 轮的数据，因此在 $n$ 上完全独立、可并行；两个 $\mathcal G$ 项构成一次廉价的顺序扫描。它也可以读成一次**预条件的定常迭代**：记 $M_f$、$M_g$ 分别为细、粗求解器的全时间矩阵，则上式等价于

$$
M_g(z)\,\Delta U^k=b-M_f(z)U^k,
\qquad U^{k+1}=U^k+\Delta U^k,
$$

即用块下二对角的 $M_g$ 预条件同为块下二对角的 $M_f$。并行性全部落在残量 $b-M_f U^k$ 的形成上。这个读法是本专题另外两页的桥梁：把 $M_g$ 换成一个可 FFT 对角化的矩阵，就得到[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|对角化技术]]与[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|全时间预条件]]两条路线。

### 误差传播与两个收敛因子

若 $A=V_ADV_A^{-1}$，且粗、细传播子是稳定函数为 $R_g,R_f$ 的一步法，则误差按模态解耦，满足 $\boldsymbol\xi^{k+1}(z)=M(z)\boldsymbol\xi^k(z)$，其中

$$
M(z)=M_g^{-1}(z)\bigl[M_g(z)-M_f(z)\bigr]=I_t-M_g^{-1}(z)M_f(z),
$$

$M_g(z)$ 与 $M_f(z)$ 是单位对角、次对角分别为 $-R_g(z)$ 与 $-R_f^{J}(z/J)$ 的下三角 Toeplitz 矩阵，$z=\Delta T\lambda$。注意这里的定标：细求解器在一个粗区间内走 $J$ 步、每步长 $\Delta T/J$，所以出现的是 $R_f^{J}(z/J)$ 而不是 $R_f(z)$。把标量因子提出来得到

$$
M(z)=\bigl[R_f^J(z/J)-R_g(z)\bigr]\,\widetilde M\bigl(R_g(z)\bigr),
$$

其中 $\widetilde M(\beta)$ 是次对角依次为 $1,\beta,\beta^2,\dots$ 的严格下三角 Toeplitz 矩阵。**这一步是整条文献的引擎**：它把一个矩阵幂的估计化归为一个标量因子乘以一个只依赖 $R_g$ 的固定 Toeplitz 矩阵，而后者的 $\|\widetilde M^k\|_\infty$ 由 Gander 与 Vandewalle 的引理控制。由此得到线性收敛因子与超线性估计

$$
\varrho_l(J,z)=\frac{\bigl|R_g(z)-R_f^{J}(z/J)\bigr|}{1-|R_g(z)|},
\qquad
\varrho_s(J,z,N_t,k)=\frac{\bigl|R_g(z)-R_f^J(z/J)\bigr|^k}{k!}\prod_{j=1}^k(N_t-j),
$$

并有 $\max_n\|\boldsymbol e_n^k\|_\infty\le\bigl(\max_{z\in\sigma(\Delta TA)}\varrho_l\bigr)^k\max_n\|\boldsymbol e_n^0\|_\infty$。$\varrho_s$ 在 $k=N_t$ 处为零，对应精确算术下的有限步终止；$\varrho_l$ 要求 $|R_g(z)|<1$，是长时间区间上真正起作用的那一个。

$\varrho_l$ 的两部分各自有明确含义：**分子是一个粗步内粗、细传播子的失配，分母是粗传播子的耗散裕量。** 前者说明校正量有多大，后者说明这个校正沿顺序粗扫描累积时会被压缩多少。

### 非线性情形

细传播子取精确、粗传播子为 $p$ 阶且局部截断误差不超过 $C_3\Delta T^{p+1}$，并设 Lipschitz 条件 $\|\mathcal G(\cdot,\boldsymbol v)-\mathcal G(\cdot,\boldsymbol w)\|\le(1+C_2\Delta T)\|\boldsymbol v-\boldsymbol w\|$ 与展开 $\mathcal F-\mathcal G=c_{p+1}(\boldsymbol v)\Delta T^{p+1}+\cdots$（$c_j$ 连续可微），则（Gander 与 Hairer，2008）

$$
\|\boldsymbol u(T_n)-\boldsymbol u_n^k\|\le
\frac{C_3\Delta T^{p+1}\bigl(C_1\Delta T^{p+1}\bigr)^{k+1}}{(k+1)!}
\,(1+C_2\Delta T)^{n-k-1}\prod_{j=0}^{k}(n-j).
$$

这是超线性型的界，与 $\varrho_s$ 同构：阶乘压过多项式，因此短时间窗口上收敛极快，长时间窗口上要靠 $\varrho_l$ 型的线性界。

### 这套装置在实践中的边界

线性理论给出的常数在耗散问题上很准，在传播型问题上会失效，而失效是连续发生的。综述报告的一组直接算例可以作为标尺：对流扩散方程，$T=4$、$\Delta T=0.1$、$\Delta x=1/128$、$J=32$，$\mathcal G$ 取后向 Euler、$\mathcal F$ 取两级二阶 SDIRK，随黏性 $\nu$ 减小 parareal 的表现单调变差，**大约在 $\nu\le10^{-3}$ 处发散**；二阶波动方程上则根本不收敛。作为另一端的标尺，Gander 与 Vandewalle 在连续层面（$\mathcal F$ 精确）给出的收缩可以远好于 $0.3$，例如 Radau IIA 约为 $0.068$。

> [!note] 这些算例的出处
> 上述 $\nu$ 扫描与 Radau IIA 的 $0.068$ 出自 Gander、Wu 与 Zhou 的 _Acta Numerica_ 综述（编号 85）与 Gander–Vandewalle（2007，Table 5.1），**不是**本页四篇论文自己的实验。本页引用它们是为了给下面各篇的常数一个可比的量级，而不是把它们归给这些论文。

**这套装置里还有一个被默认的结构假设，编号 77 正是要拆掉它。** $\mathcal F$ 在一个粗区间上的作用被写成单个标量稳定函数在单个自变量上的 $J$ 次**幂** $R_f^{J}(z/J)$，这一步要求区间内 $J$ 个细步等长。

## 12：细传播子只是 A-稳定时会发生什么

### 直觉

parareal 为什么会收敛？因为在**刚性**问题上，粗、细传播子在最危险的那些模态上其实是一致的：高频模态被两者同样地杀死，所以它们的差很小；低频模态被两者同样准确地积分，所以它们的差也很小；只有中间那一段频率上粗传播子明显不如细传播子，而那一段的宽度有限。这就是 $\approx0.3$ 这个与网格无关、与 $T$ 无关的常数的来源——它是 $\varrho_l(J,\cdot)$ 在中频段的一个峰值，两端都被压下去了。

一旦细传播子只是 A-稳定而不是 L-稳定，**高频这一端就不再被压下去**。梯形规则在 $z\to\infty$ 时 $R_f\to-1$，四阶 Gauss 方法 $R_f\to+1$：它们不衰减高频，只是把高频的振幅原样保留（甚至逐步反号）。而粗传播子后向 Euler 仍然衰减高频。于是二者在高频上的差趋于 $1$，$\varrho_l$ 也趋于 $1$，parareal 在**恰恰使问题成为刚性问题的那些模态上停止收缩**。

论文的救援方式是承认离散问题的谱是有界的：只需要 $z\in[0,z_{\max}]$ 上一致小的收缩，而不是全轴上的。给定 $z_{\max}$，只要粗细比 $J$ 足够大，细传播子在每个细步上看到的自变量 $z/J$ 就落回它仍然衰减的区域，$J$ 次幂重新提供衰减。剩下的问题就只是**这个「足够大」有多大**——这正是论文的临界粗细比 $J_{\rm cri}$。

### 问题设定

$A$ 对称正定，来自抛物 PDE 的空间半离散，允许含分数阶 Laplacian（此时 $A$ 稠密、$\lambda_{\max}$ 极大）。粗传播子固定为后向 Euler。细传播子取三种：梯形规则、三阶两级 DIRK、四阶两级 Gauss RK。以 $\boldsymbol u'+A\boldsymbol u=g$ 的约定（$z=\Delta T\lambda$，$\lambda>0$），相关稳定函数为

$$
\text{后向 Euler：}\ R(z)=\frac{1}{1+z},
\qquad
\text{梯形：}\ R(z)=\frac{2-z}{2+z},
\qquad
\text{四阶 Gauss：}\ R(z)=\frac{12-6z+z^2}{12+6z+z^2}.
$$

后向 Euler 是 L-稳定的，$R(\infty)=0$；梯形规则 A-稳定但不 L-稳定，$R(\infty)=-1$；四阶 Gauss 方法 A-稳定且辛，$R(\infty)=+1$。三阶两级 DIRK 在这条线索中通常取的 Butcher 表是

$$
\begin{array}{c|cc}
\gamma & \gamma & 0\\
1-\gamma & -\tfrac{1}{\sqrt3} & \gamma\\ \hline
 & \tfrac12 & \tfrac12
\end{array},
\qquad \gamma=\frac{3+\sqrt3}{6},
$$

其稳定函数是 $(2,2)$ 型有理函数，A-稳定但不 L-稳定。（该表出自作者自己的综述，编号 12 是否逐字采用同一张表本站未核实。）

### 推导：高频极限直接说明旧论证为什么失效

把 $z\to\infty$ 代入 $\varrho_l$ 即可。$\mathcal G$ 为后向 Euler 时 $R_g(z)=1/(1+z)\to0$，分母 $\to1$；而 $R_f^J(z/J)\to R_f(\infty)^J$，因此

$$
\lim_{z\to\infty}\varrho_l(J,z)=\bigl|R_f(\infty)\bigr|^{J}.
$$

四种细传播子的结果一目了然：

| 细传播子 $\mathcal F$      | $R_f(\infty)$ | $\lim_{z\to\infty}\varrho_l(J,z)$ |
| -------------------------- | ------------- | --------------------------------- |
| 后向 Euler（L-稳定）       | $0$           | $0$                               |
| 精确传播子 $e^{-z}$        | $0$           | $0$                               |
| 梯形规则（仅 A-稳定）      | $-1$          | $1$                               |
| 四阶 Gauss RK（仅 A-稳定） | $+1$          | $1$                               |

（由上面 $\varrho_l$ 的表达式与三个稳定函数直接计算，$\mathcal G$ 取后向 Euler。）

上半表就是经典结论成立的原因：L-稳定或精确的 $\mathcal F$ 使 $\varrho_l$ 在全实轴上的上确界由中频段决定，从而对所有 $J\ge J_{\min}=O(1)$（$\mathcal F$ 精确时甚至对所有 $J\ge2$）都有 $\max_{z\in\mathbb R^-}\varrho_l\approx0.3$。下半表说明为什么同一论证对梯形规则与四阶 Gauss 方法完全失效：$\varrho_l$ 的上确界是 $1$，与 $J$ 无关，无论把 $J$ 取多大都不会改变全轴上的上确界。

因此唯一的出路是把 $z$ 限制在实际出现的有限区间 $[0,z_{\max}]$ 上，$z_{\max}=\Delta T\lambda_{\max}$。这时 $J$ 重新起作用：$J$ 越大，$R_f$ 被求值的点 $z/J$ 越靠近细传播子仍有衰减的区域，$J$ 次幂把衰减放大。论文给出的正是这个门槛的定量形式。

### 定理

**（A-稳定但不 L-稳定的细传播子。）** 设 $\mathcal G$ 为后向 Euler、$\mathcal F$ 为梯形规则或四阶 Gauss RK，则

$$
\max_{z\in[0,z_{\max}]}\varrho_l(J,z)\approx0.3
\qquad\forall\,J\ge J_{\min}=O\bigl(\log^2 (z_{\max})\bigr).
$$

门槛在 $z_{\max}=\Delta T\lambda_{\max}$ 上的**对数平方**增长，是摘要中「$J_{\rm cri}$ 依赖 $\Delta T$、$\Delta t$ 与 $\lambda_{\max}$」这句话的精确形式，也是全文的核心技术内容。

**（L-稳定细传播子的对照。）** 若 $\mathcal F$ 是 L-稳定的，则 $\max_{z\in\mathbb R^-}\varrho_l\approx0.3$ 对所有 $J\ge J_{\min}=O(1)$ 成立；若 $\mathcal F=\exp(\Delta TA)$ 精确，该速率对所有 $J\ge2$ 已经成立。

**（三阶 DIRK 的特殊性。）** 对三阶两级 DIRK 细传播子，$J_{\rm cri}=4$，**与 $\Delta T$、$\Delta t$、$\lambda_{\max}$ 全部无关**。这是摘要中三个结论里唯一一个不带参数依赖的，也是这条线索里最实用的一条：一个固定的、可以直接写进代码的粗细比。

> [!note] 未核实的部分
> 论文声明给出计算 $J_{\rm cri}$ 的**简明闭式公式**（对梯形规则与四阶 Gauss RK），这些公式本身在任何可获取的来源中都读不到，本站不给出。定理编号、对 $g$ 的确切假设、以及是否有非线性推广，同样未核实。

### 数值实验

摘要确认的部分只有算例的**类别**：「含分数阶 PDE 与不确定性量化应用的数值算例被用来支持理论预测」。这两类不是随手选的：

- **分数阶 Laplacian $(-\Delta)^{\alpha}$ 的时间依赖 PDE。** 空间半离散给出稠密的对称正定 $A$（一维为 Toeplitz 结构），$\lambda_{\max}$ 极大，因此 $z_{\max}=\Delta T\lambda_{\max}$ 极大——这正是 $J_{\rm cri}$ 的公式产生实际差别的区域。若门槛只有 $O(1)$，这一类算例就没有信息量。
- **不确定性量化。** 同一个 ODE 系统的大量参数实现各自需要长时间积分，是时间并行真正能兑现收益的场景，同时也检验收缩因子对系数矩阵谱的**稳健性**：不同实现给出不同的 $\sigma(A)$，而理论宣称的界不依赖它。

具体的网格规模、$J$ 的取值、迭代次数、实测收缩因子与加速比，本站**均未核实**（论文无预印本，出版社全文对自动化访问不可得）。

最接近「被报告的测量值」的，是这条文献中围绕本文流传的一组常数。它们出自第三方论文的引言（Zhou、Liu 与 Wu，_Parareal-CG_，arXiv:2304.10152 §1）对本文的转述，可核实为**转述**而非原文排版：

| 粗/细组合                      | 收缩因子       | 门槛                                | 归属                   |
| ------------------------------ | -------------- | ----------------------------------- | ---------------------- |
| Parareal-Euler                 | $\approx0.298$ | $J\ge2$                             | Mathew–Sarkis–Schaerer |
| Parareal-2sDIRK                | $0.316$        | $J\ge2$                             | Wu                     |
| Parareal-TR/BDF2               | $0.333$        | $J\ge2$                             | Wu                     |
| 三阶 DIRK 细传播子             | $\approx0.333$ | $J\ge4$                             | **编号 12**            |
| 梯形规则 / 四阶 Gauss 细传播子 | $\approx0.333$ | $J\ge J^*_{\min}(\rho(A),\Delta t)$ | **编号 12**            |

> [!warning] 关于这组常数的舍入
> 论文本身是否排印 $0.333$、还是别的舍入，本站未核实；作者自己的综述把同一结论四舍五入为 $\approx0.3$。综述另外记录了两个具体门槛：两级二阶 SDIRK 为 $J_{\min}=2$（Wu, IMA J. Numer. Anal. 2015），两级三阶 SDIRK 为 $J_{\min}=4$（本文）——后者与上表第四行一致。

**这些实验建立了什么、又缺什么。** 能建立的是：门槛现象真实存在，且在极大 $\lambda_{\max}$ 的分数阶算例上仍与预测一致；不能建立的是门槛公式的**紧性**——没有可核实的数据说明 $J$ 略小于 $J_{\rm cri}$ 时收缩因子退化得有多快，也没有可核实的墙钟时间说明在真实并行机上这套选择带来多少加速。

### 与其他论文的关系

这是本系列 parareal 分析的**奠基之作**，后续工作都以它为基准。本页的编号 20 与编号 30 把同一个「哪一对 $\mathcal F/\mathcal G$ 稳健收缩」的问题带进分数阶问题；编号 77 拆掉 $R_f^J(z/J)$ 所依赖的均匀细网格假设。换一个分析对象则得到[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|编号 39]]：那里分析的是两级 MGRIT，并给出与本文可直接比较的常数 $0.2984$（parareal）与 $0.1115$（MGRIT-FCF）。相关结论在[[computational-mathematics/knowledge-notes/time-parallelization/index|时间并行综述精读]]第四章以 $(4.8)$ 的形式出现。

## 20：空间分数阶让粗传播子成为瓶颈

### 直觉

在经典抛物问题上，parareal 的成本模型很简单：细传播贵但并行，粗传播便宜但顺序，只要粗传播的总代价远小于细传播的一轮，加速就成立。**空间分数阶把这个成本模型打破了。** 分数阶 Laplacian 是非局部算子，半离散后的 $A$ 稠密，因此「解一次隐式方程」这件事本身变得昂贵——而粗传播子每轮都要在 $N_t$ 个粗区间上顺序地做这件事，直接落在关键路径上。于是即使收敛只需要三四轮，运行时间也可能被那条顺序链条吃掉。

同时，稠密的 $A$ 带来极大的 $\lambda_{\max}$，正好落在编号 12 中门槛按 $O(\log^2(\Delta T\lambda_{\max}))$ 增长的区域。两个困难是叠加的：粗传播越贵，越希望取大的 $\Delta T$（少的粗区间），而大的 $\Delta T$ 又抬高 $z_{\max}$、抬高需要的 $J$。题名中的「快速」指的是加速 parareal **迭代本身**，而不只是给出收缩因子的界。

### 问题设定

含分数阶 Laplacian $(-\Delta)^{\alpha}$ 的时间依赖 PDE，空间半离散后得到刚性系统 $\boldsymbol u'+A\boldsymbol u=g$，$A$ 稠密（一维情形有 Toeplitz 结构）、对称正定、$\lambda_{\max}$ 极大。第三方文献明确记载该文在**常步长**下分析这一设定下的 parareal 收敛性质。

> [!warning] 可核实范围
> 该文摘要未能从 Crossref、OpenAlex、Semantic Scholar、NASA ADS 或出版社页面获取（Elsevier 未向任何聚合站提交摘要，ScienceDirect 与 NASA ADS 阻断自动化检索），正文亦不可得。可确认的只有上述问题设定与「常步长」这一条第三方记载。**其具体的分数阶离散方式（Riesz / 移位 Grünwald / 矩阵转移 $A=Q^\alpha$ 或其他）、算法构造、定理、收缩因子与全部数值算例，本站均未核实，因此不报告任何数值常数、表格或实验结论。**
>
> 需要特别指出的一点：同一作者在相邻时间发表的几篇工作中出现过「用 IMEX-Euler 显式处理稠密分数阶算子作为廉价粗传播子」与「一级复 Rosenbrock 粗传播子」等机制（分别见 Wu, _Appl. Math. Comput._ 2017，DOI `10.1016/j.amc.2017.02.019`；Wu 与 Zhou, _Math. Methods Appl. Sci._ 2017，DOI `10.1002/mma.4273`）。**这些机制不能在未核对原文的情况下归给编号 20**，本站只把它们记为邻近文献中的候选思路。

### 数值实验

**未核实。** 由题名与引用它的文献可以推断测试类别几乎必然是含分数阶 Laplacian 的时间依赖 PDE，但这是推断而非核实，因此本站不列任何设置或结果。

### 与其他论文的关系

编号 12 是直接前身：它建立了对称正定 $A$ 上的 $J_{\rm cri}$ 框架，并已把「分数阶 PDE 应用」作为动机算例。编号 30 是另一种分数阶的对应物——编号 20 处理**空间**分数阶（空间非局部、时间局部，parareal 的结构得以保留），编号 30 处理**时间**分数阶（时间非局部，破坏子区间之间的独立性）。第三方文献明确作出这一区分。至于[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|编号 31]]，它同样面对分数阶 Laplacian，却完全放弃迭代，转向时间周期情形的直接对角化。

## 30：时间分数阶让负载不再均衡

### 直觉

parareal 的全部前提是：$[T_n,T_{n+1}]$ 上的细传播只需要局部初值 $\boldsymbol u_n^k$。时间分数阶导数使这一条不成立——Caputo 或 Riemann-Liouville 导数在时刻 $t$ 依赖整个 $[0,t]$ 上的解历史。

值得注意的是，**这里的失效方式不是「不收敛」，而是「并行度消失」。** 朴素推广仍然可以写出来：在第 $n$ 个子区间上推进时，把 $[0,T_n]$ 上的历史也一并带上。迭代照样收敛，但处理器 $n$ 的工作量正比于 $n$，最后一个处理器的工作量是第一个的 $N_t$ 倍。总工作量从 $O(N_t)$ 变成

$$
\sum_{n=1}^{N_t}n=\frac{N_t(N_t+1)}{2},
$$

而并行时间由最慢的那个处理器决定，因此关键路径长度也是 $O(N_t)$ 个细步——**与串行时间推进同阶**。这就是摘要说的「计算时间在各进程间不均衡」：算法在纸面上并行，在机器上不并行。

修好它的思路不是改迭代，而是改**方程的写法**：把非局部的历史压缩进有限个满足局部常微分方程的辅助变量。一旦如此，每次细传播调用在每个子区间上代价相同，负载自动均衡，parareal 的原始成本模型恢复。

### 问题设定

时间分数阶微分方程，分数阶算子具有历史效应。目标是设计一个**代价在各进程间均衡**且收敛速率稳健的 parareal 算法。摘要给出的路线是采用两个「近期发展出来的时间分数阶算子的局部时间积分器」，两者都通过引入辅助变量来局部化分数阶算子。

### 推导

**第一步：局部化。** 设增广状态为 $(\boldsymbol u_n^k,\boldsymbol z_n^k)$，其中 $\boldsymbol z$ 是辅助变量，其演化由局部动力学支配。分数阶卷积被替换为有限个额外的局部未知量，因此第 $n$ 个子区间上的细传播不再需要 $[0,T_n]$ 上的历史，只需要 $(\boldsymbol u_n^k,\boldsymbol z_n^k)$。这一步把「工作量正比于 $n$」变回「工作量与 $n$ 无关」。

**第二步：不能照搬 parareal 更新式。** 把单一的 parareal 更新作用在整个增广状态上，即

$$
\begin{pmatrix}\boldsymbol u_{n+1}^{k+1}\\ \boldsymbol z_{n+1}^{k+1}\end{pmatrix}
=\mathcal F\!\begin{pmatrix}\boldsymbol u_{n}^{k}\\ \boldsymbol z_{n}^{k}\end{pmatrix}
+\mathcal G\!\begin{pmatrix}\boldsymbol u_{n}^{k+1}\\ \boldsymbol z_{n}^{k+1}\end{pmatrix}
-\mathcal G\!\begin{pmatrix}\boldsymbol u_{n}^{k}\\ \boldsymbol z_{n}^{k}\end{pmatrix},
$$

正是论文说明**不采用**的做法。真正新的算法成分是一种**混合式粗网格校正**：辅助变量与解变量**分开**校正，两块使用不同的校正规则。

**为什么必须分开。** 解变量与辅助变量在动力学上的角色不同：辅助变量是为逼近记忆核而引入的，其模态的时间尺度分布很宽（这正是核压缩的目的），把它们和解变量一起用同一个粗传播子校正，等于用同一个 $\Delta T$ 去处理跨若干数量级的时间尺度。分开校正使两类变量各自获得合适的处理。（这一句是对「分开校正」这一已核实事实的机理解读，不是论文的论证。）

> [!note] 可核实范围
> 上式（即论文**不用**的那个朴素增广更新）是按摘要描述重建的；**混合校正的确切公式本站未核实**。所采用的两个局部时间积分器的名称同样未核实。与摘要描述相符的标准候选有两类：其一是把分数阶核用指数和 / 有理逼近压缩，从而把卷积换成少量辅助 ODE 模态；其二是 Yuan–Agrawal 型扩散表示，把分数阶导数写成对连续辅助模态的积分，再用求积化为有限个辅助变量。**这两类只是文献中的标准候选，论文实际引用哪两篇未核实。**

### 定理

摘要确认的表述是：所提出的 parareal 算法**具有稳健的收敛速率**。在这一系列工作中，「稳健」一贯指与空间算子的特征值、粗细比 $J$ 以及粗区间数 $N_t$ 无关。

**收敛因子的显式值、分数阶阶数 $\alpha$ 的允许范围、$\mathcal G$ 与 $\mathcal F$ 的取法、辅助变量动力学所需的稳定性条件、定理编号，本站均未核实，因此不给出任何数字。** 负载均衡的说法是否被定量化（例如「每处理器代价与 $n$ 无关」是否作为命题写出）同样未核实。

### 数值实验

摘要只确认「给出了数值算例以支持结论」。NASA ADS 记录的关键词是 _Parareal; Time-fractional differential equations; Local time-integrators_。**具体的测试方程、分数阶阶数、网格规模、进程数、实测迭代次数与加速比，本站均未核实。**

值得说明的是这里缺失的信息为什么重要：本文的核心主张是**并行效率**（负载均衡）而不只是收敛速率，而并行效率恰恰是唯一必须靠实测才能确认的量。收敛速率还可以从理论侧间接判断，负载均衡不能。因此在拿到原文的实验表之前，本站对该文效率主张的支持程度只能停留在「机制可信」这一层。

### 与其他论文的关系

编号 20 是空间分数阶的姊妹篇；两者的技术困难没有交集——历史效应导致的负载不均衡在编号 20 中不存在。编号 12 提供了分析风格（稳健、与参数无关的收缩因子），编号 30 宣称在分数阶设定下复现它。[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|编号 31]] 从相反方向进攻同类问题：不修复迭代法，而用直接对角化绕开收敛因子问题。此外，「用辅助变量把非局部算子局部化」这一手法在结构上会在[[computational-mathematics/paper-notes/parallel-in-time/all-at-once-preconditioners|编号 71]] 再次出现——那里的非局部性来自前向与后向演化的耦合，而不是分数阶核，但对策同样是扩大状态并预条件扩大后的系统。

## 77：拆掉「细网格等长」这个假设

### 直觉

自适应时间步、初值不光滑或不相容时 $t=0$ 附近的分级网格、快速瞬态附近的局部细化，都会产生非均匀细网格。把 parareal 理论限制在均匀细网格上，等于排除了几乎全部自适应实践。所以这个推广的动机是显然的。

不显然的是它为什么困难。直觉上，非均匀只是把「$J$ 个等长细步」换成「$J_n$ 个不等长细步」，看起来像是记号问题。**但收敛分析依赖的不是细步的数量，而是细传播子在整个时间轴上的「同一性」。** 均匀网格下，每个粗区间上的细传播都是同一个标量 $R_f^J(z/J)$，因此 $M_f(z)$ 是 Toeplitz 的；非均匀网格下，每个粗区间上的细传播是各自不同的乘积，$M_f(z)$ 逐行变化。收敛分析的引擎——把误差矩阵分解成「一个标量因子乘一个固定的 Toeplitz 矩阵」——就没有东西可以提出来了。

### 问题设定

该文研究**均匀粗网格**配**任意分布的非均匀细网格**时的 parareal 收敛行为。「任意分布」是实质所在：粗区间内的细网格既不假设均匀，也不假设分级，甚至不假设各粗区间之间相同。迭代式本身不变：

$$
\boldsymbol u_{n+1}^{k+1}=\mathcal F(T_n,T_{n+1},\boldsymbol u_n^{k})
+\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k+1})
-\mathcal G(T_n,T_{n+1},\boldsymbol u_n^{k}),
\qquad T_n=n\Delta T .
$$

### 推导：旧装置在哪一步断掉

所有经典线性收敛结果——Gander 与 Vandewalle 的界、$0.3$ 附近的常数、编号 12 的临界粗细比公式、编号 39 的 MGRIT 因子——都依赖恒等式

$$
\mathcal F(T_n,T_{n+1},\cdot)\ \longleftrightarrow\ R_f^{J}\!\Bigl(\frac{z}{J}\Bigr),
$$

即粗区间上的细传播是单个标量稳定函数在单个自变量上的 $J$ 次幂。细网格非均匀时得到的是**互不相同因子的乘积**：

$$
\mathcal F(T_n,T_{n+1},\cdot)\ \longleftrightarrow\ \prod_{i=1}^{J_n}R_f(\theta_{n,i}z),
\qquad
\theta_{n,i}=\frac{\Delta t_{n,i}}{\Delta T},\quad \sum_i\theta_{n,i}=1 .
$$

两件事随之失效。**第一，$\varrho_l(J,z)$ 不再有定义**，因为没有单一的 $J$：分子里的 $R_f^J(z/J)$ 要换成上面那个乘积，而这个乘积依赖 $n$。**第二，也更严重：若细网格在各粗区间之间不同**，则 $M_f(z)$ 的次对角元素逐行变化，$M_f(z)$ **不再是 Toeplitz** 矩阵，因式分解

$$
M(z)=\bigl[R_f^J(z/J)-R_g(z)\bigr]\,\widetilde M\bigl(R_g(z)\bigr)
$$

失效——右端那个标量因子根本提不出来。而 Gander 与 Vandewalle 关于 $\|\widetilde M^k\|_\infty$ 的引理（本页开头那个「引擎」）正是作用在 $\widetilde M$ 上的，没有这个分解就无法应用。在这一设定下恢复稳健的、与网格无关的收缩因子因此是全新的分析，而不是常规推广。

> [!warning] 这一节的性质
> 以上「为什么断掉」的论证是本页从编号 12 与 39 的**已核实公式**出发做的推理，逻辑上是这些公式的严格推论，但**不是**从编号 77 原文读到的说法。该文自己如何刻画这一困难、以及它用什么替代 $\varrho_l$，本站未核实。

### 定理

**未核实。** 无法从任何可获取来源确认该文的定理陈述、收敛因子或门槛，因此本站**不把任何数值常数归给它**。

可核实并且有助于定位的周边信息有两项。其一，MaRDI/zbMATH 记录的参考文献列表显示该文对标的是经典 parareal 收敛理论的完整阅读清单：编号 12、Gander 与 Hairer 的非线性收敛分析、Bal 的收敛与稳定性、Gander 与 Vandewalle 的分析、Wu 关于二阶 parareal 算法的收敛分析（IMA J. Numer. Anal. 2014）、Wu 关于并行粗网格校正的工作（SISC 2018）、Southworth 的必要条件与紧的两级界等。这确认它是一篇**纯分析**工作，而不是新算法工作。

其二，用于对照的量级来自最接近的均匀网格结果（Yang、Yuan 与 Zhou，_Robust Convergence of Parareal Algorithms with Arbitrarily High-Order Fine Propagators_，CSIAM Trans. Appl. Math.，arXiv:2109.05203），其摘要可核实：

| 假设                                                                  | 结论                                                     |
| --------------------------------------------------------------------- | -------------------------------------------------------- |
| $\mathcal G$ 后向 Euler，单步细传播子满足 $\lvert r(-\infty)\rvert<1$ | 存在临界 $J_*$，使所有 $J\ge J_*$ 上线性收敛率接近 $0.3$ |
| 同上，且初值不光滑 / 不相容                                           | 该速率仍然稳健                                           |
| $\mathcal F$ 取二、三、四级 Lobatto IIIC                              | 收缩因子 $<0.31$，且 $J_*=2$                             |

**编号 77 在非均匀细网格上是否达到可比的常数，本站未核实。**

### 数值实验

**未核实。** 无法确认任何测试问题或结果。本站不列表、不给数字。

### 与其他论文的关系

编号 12 是被推广的对象：它固定 $\mathcal G$ 为后向 Euler、允许三种 $\mathcal F$，在**均匀细网格**上找临界比 $J_{\rm cri}$；编号 77 保留均匀粗网格但完全放开细网格。编号 12 的依赖对象是标量 $J$，编号 77 里根本没有一个标量 $J$ 可以依赖。[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|编号 39]] 通过 MGRIT 收敛因子中的 $|R_f^J(z/J)|$ 继承了同一条均匀细网格限制，因此编号 77 的技术是通向非均匀网格 MGRIT 理论的自然路径。编号 20 与编号 30 被第三方明确记载为在「常步长」下做分析，编号 77 正是把这条长期假设从作者自己的工作中移除的那一篇。

## 一处值得指出的张力

非均匀时间步在这两条路线上的地位恰好相反。在 parareal 分支中，它是**理论的障碍**：它破坏 Toeplitz 结构，从而破坏全部收敛因子论证。而在直接对角化分支中（见[[computational-mathematics/paper-notes/parallel-in-time/diagonalization-technique|对角化技术]]），非均匀步长是**必要条件**：Maday 与 Rønquist 需要各 $\Delta t_n$ 互不相同，才能让时间矩阵可对角化。同一个建模自由度，在一条路线上是麻烦，在另一条路线上是前提。

顺着这条对照还能看出两条路线的代价结构不同。parareal 分支的代价是一个**收缩因子**，它由 $\mathcal F$ 与 $\mathcal G$ 的稳定函数决定，因此随问题类型（耗散还是传播）剧烈变化；对角化分支的代价是一个**条件数**，它由时间矩阵的结构决定，与问题的物理性质基本无关。这解释了为什么随黏性减小 parareal 会在 $\nu\approx10^{-3}$ 处崩溃，而 $\alpha$-循环预条件对双曲问题仍然可用。

## 覆盖核对

| 内容                             | 论文   | 覆盖状态                                                   |
| -------------------------------- | ------ | ---------------------------------------------------------- |
| parareal 迭代、预条件读法        | 12     | 迭代式、$M_g\Delta U^k=b-M_fU^k$、并行性所在               |
| 误差传播矩阵与两个收敛因子       | 12     | $M(z)$、因式分解、$\varrho_l$ 与 $\varrho_s$、有限终止     |
| 非线性超线性界                   | —      | Gander–Hairer 界（背景，非本页论文）                       |
| 装置的实践边界                   | —      | $\nu$ 扫描与 Radau IIA（出自综述，非本页论文）             |
| 四个积分器的稳定函数与高频极限   | 12     | $R(\infty)$ 表、$\lim\varrho_l=\lvert R_f(\infty)\rvert^J$ |
| A-稳定情形的门槛定理             | 12     | $J_{\min}=O(\log^2 z_{\max})$、L-稳定与精确的对照          |
| 三阶 DIRK 的参数无关门槛         | 12     | $J_{\rm cri}=4$                                            |
| 流传的常数表与其舍入警告         | 12     | 五行常数（第三方转述）、$0.333$ 与 $0.3$ 的差异            |
| 实验类别与它们各自检验什么       | 12     | 分数阶 PDE、UQ；网格与加速比未核实                         |
| 空间分数阶下粗传播子成为瓶颈     | 20     | 稠密 $A$、极大 $\lambda_{\max}$（限定核实）                |
| 邻近文献机制不得归属             | 20     | 明确列出不可归属的候选机制                                 |
| 时间分数阶的历史效应与负载不均衡 | 30     | 前提失效、$\sum n=N_t(N_t+1)/2$、关键路径                  |
| 局部化与混合式粗网格校正         | 30     | 辅助变量、被弃用的朴素更新式、分开校正（限定核实）         |
| 非均匀细网格破坏 Toeplitz 结构   | 77     | 幂到乘积、两处失效、为何是新分析                           |
| 非均匀网格结果的对照量级         | 77     | Yang–Yuan–Zhou 的三行结论（外部可核实）                    |
| 非均匀步长在两条路线上的相反地位 | 12, 77 | 障碍与前提的对照、两种代价结构                             |

## 本页原文

- S. Wu and T. Zhou, [_Convergence analysis for three parareal solvers_](https://doi.org/10.1137/140970756), SIAM J. Sci. Comput. 37(2) (2015), pp. A970-A992。
- S. Wu and T. Zhou, [_Fast parareal iterations for fractional diffusion equations_](https://doi.org/10.1016/j.jcp.2016.10.046), J. Comput. Phys. 329 (2017), pp. 210-226。
- S. Wu and T. Zhou, [_Parareal algorithms with local time-integrators for time fractional differential equations_](https://doi.org/10.1016/j.jcp.2017.12.029), J. Comput. Phys. 358 (2018), pp. 135-149。
- S.-L. Wu and T. Zhou, [_Convergence analysis of the parareal algorithm with nonuniform fine time grid_](https://doi.org/10.1137/23M1592481), SIAM J. Numer. Anal. 62(5) (2024), pp. 2308-2330。
