---
title: 3.1–3.2：历史脉络与 Schwarz 波形松弛
description: 从波形松弛与区域分解出发，完整推导 OSWR、Robin minimax 参数、有限步收敛与 tent pitching，并补足 Dirichlet 多子域理论
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-3-1-history-and-swr
tags:
  - 时间并行
  - SWR
  - 双曲-PDE
---

> [!note] 阅读范围
> 本页对应论文 Sections 3.1–3.2（pp. 396–405），覆盖历史发展、公式 (3.1)–(3.4)、Theorems 3.1–3.2 和 Figures 3.1–3.3。原论文图像保持原样；正文、推导说明与 intuition 均重新撰写。

> [!info] 图表许可
> 原论文以 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) 发布。本页复制 Figure 3.1–3.3 的完整图形，并在相邻文字中给出中文图解和出处。

## 3.1 历史发展

### 四条方法线索怎样形成

论文把四类对双曲问题有效的 PinT 方法放在一起讨论，并用第二章的四个 PDE 检查它们的理论性质。之所以要先梳理源流，是因为每一类方法的并行性都来自一个不同的数学结构，理解来源才能理解它们各自擅长和失效的场景。

1. **Schwarz 波形松弛（SWR）**来自区域分解与波形松弛的结合。连续时空子问题的做法由 Gander (1999) 首次针对抛物问题提出，Giladi and Keller (2002) 独立给出相近构造。它的两条根系分别是：区域分解（DD），一种可追溯到 Schwarz (1870) 的并行求解 PDE 的经典技术；以及波形松弛（WR），源自 Lelarasmee et al. (1982) 的电路仿真。Gander (1997) 同时对抛物与双曲情形做了发展与分析，SWR 这一名称由 Gander et al. (1999) 正式提出。非线性抛物问题的后续结果见 Gander (1999) 与 Gander and Rohde (2005)。使用更有效传输条件的优化版本（OSWR）在抛物情形由 Gander and Halpern (2007)、Bennequin, Gander and Halpern (2009) 和 Bennequin, Gander, Gouarin and Halpern (2016) 发展，在双曲情形由 Gander, Halpern and Nataf (2003) 和 Gander and Halpern (2004) 发展；非线性对流扩散见 Gander, Lunowa and Rohde (2023)。此外还有 Dirichlet–Neumann 和 Neumann–Neumann 波形松弛变体（Gander, Kwok and Mandal 2016b, 2021b）。Ciaramella, Gander and Mazzieri (2023) 的 unmapped tent pitching（UTP）也建立在 SWR 上。
2. **并行积分延迟校正（IDC）**源于 Böhmer and Stetter (1984) 的演化问题校正思想。Dutt, Greengard and Rokhlin (2000) 将其识别为一种能通过精确处理积分项、逐轮提高阶数的专用时间积分器。Guibert and Tromeur-Dervout (2007) 提出流水线版 PIDC，Christlieb et al. (2010) 提出修订版 RIDC，二者都能在时间方向并行。
3. **ParaExp**由 Gander and Güttel (2013) 提出。它把初值传播和源项响应分开处理，随后产生了线性实现（Merkel et al. 2017；Kooij et al. 2017）和非线性迭代扩展（Gander, Güttel and Petcu 2018a）。
4. **ParaDiag**从 Maday and Rønquist (2008) 的直接对角化方法发展而来，作为不含迭代的直接时间并行求解器。后续研究扩展到抛物（Gander et al. 2016a）、双曲（Gander et al. 2019）和非线性问题（Gander and Halpern 2017），也形成波形松弛（Gander and Wu 2019）、Parareal（Gander and Wu 2020）、Krylov 预条件（McDonald et al. 2018；Liu and Wu 2020）、插值（Kressner, Massei and Zhu 2023）以及 Sherman–Morrison–Woodbury 加 Krylov（Gander and Palitta 2024）等实现路径。

这段历史的作用是标明并行性的来源。SWR 利用空间接口上的整段波形，IDC 利用校正层之间的流水线，ParaExp 利用线性叠加，ParaDiag 利用全时间矩阵的谱结构。四者所依赖的结构不同，因此在双曲与抛物问题上的表现也系统性地不同——这正是后续几节要逐一验证的。

> [!tip] 本站洞见：为什么 SWR 排在第一位
> 四类方法里，只有 SWR 把"并行"直接建立在 PDE 的物理传播机制上：它先按空间切开区域，再让每个子域独立求解整段时间窗，并在人工界面上交换时间函数。对双曲问题，有限传播速度让这种界面交换在有限步内精确收敛（见 3.2.2）；对抛物问题，无限传播速度使得必须用优化传输条件来加速。方法与 PDE 结构的这种"贴合"是它能贯穿抛物/双曲两类问题的根本原因。

## 3.2 Schwarz 波形松弛（SWR）

### 从 WR 到 SWR

SWR 融合了经典 Schwarz DD 与 WR 的长处，同时克服了两者各自的局限。要看清这一点，需要先分别理解这两种前身方法为什么不够用。

### 经典空间区域分解的限制

在演化 PDE 的语境里，经典 Schwarz DD 通常先做统一的隐式时间离散，再在每个时间步上用 DD 技术顺序求解由此得到的椭圆问题（Cai 1991；Meurant 1991；Cai 1994）。这带来两个结构性代价：其一，所有子域必须在当前时间步的 DD 迭代收敛后，才能一起进入下一步，时间方向被完全串行化；其二，各子域必须共用同一套时间离散，这恰恰抹掉了 DD 最有价值的优势——为每个子域单独裁剪数值处理（不同的空间网格、时间步长或积分器）。

### 波形松弛的分裂困难

经典波形松弛换一个起点：从演化 PDE 空间离散得到的一阶 ODE 系统

$$
\boldsymbol u'(t)=A\boldsymbol u(t)+\boldsymbol f(t),
\qquad A=M+N
$$

出发，用一种类似 Picard 迭代但带有系统分块的动态迭代求解：

$$
\frac{d\boldsymbol u^k}{dt}-M\boldsymbol u^k
=N\boldsymbol u^{k-1}+\boldsymbol f,
\qquad \boldsymbol u^k(0)=\boldsymbol u_0.
$$

这里 $k\geq1$ 是迭代指标，对所有 $k$ 都有 $\boldsymbol u^k(0)=\boldsymbol u_0$，而 $(M,N)$ 是 $A=M+N$ 的一个相容分裂。分裂方式直接决定并行结构：

- **Jacobi（对角）型 $M$**：求解 $\boldsymbol u^k$ 退化为一组彼此独立的标量 ODE，可以完全并行。它之所以成为 PinT 方法，是因为所有未知量的"未来"都在其相邻未知量的未来已知之前就被近似出来——时间轴被提前推进。
- **Gauss–Seidel（三角）型分裂**：同样退化为标量 ODE，并可借助红黑或其他着色获得并行性。
- **循环约化（cyclic reduction）**：还能进一步增加并行层次（Worley 1991；Horton, Vandewalle and Worley 1995；Simoens and Vandewalle 2000）。

WR 的真正难点不在于并行，而在于**找到一个能保证快速收敛的有效分裂**。收敛速度完全取决于 $A=M+N$ 是否把强耦合保留在 $M$ 中、只把弱耦合留给 $N$。正如 Nevanlinna (1989) 所强调的：人们真正关心的是"什么样的子系统划分能让迭代快速收敛……如何分裂才能让耦合保持微弱，是一个重要问题"。这句话点出了症结：分裂通常被当作"给定的"，但一个糟糕的分裂会把强耦合项丢进 $N$，使收敛任意地慢，甚至发散，从而让 WR 在实践中不可用。

> [!tip] 本站洞见：SWR 如何"绕过"分裂难题
> SWR 的关键一步是改变解耦的对象。它不先做空间离散、再去寻找一个好的代数分裂，而是**先在连续层面切开空间区域**，然后像 WR 那样在各子域上独立求解连续时空 PDE。这样做的好处是：既然我们知道特定 PDE 下相邻子域到底通过什么机制耦合，就可以直接设计**传输条件**去逼近这种耦合，让迭代快速收敛。于是 Nevanlinna 指出的"如何分裂使耦合微弱"这个难题，被"如何设计传输条件使子问题近似解耦"所取代——后者有 PDE 的物理结构可依据，因而是可解的。这正是 OSWR 的出发点。

SWR 在空间离散前切分连续区域。每个子域独立求解完整时间窗上的 PDE，并交换人工边界处的时间函数。这样既能为各子域选择不同的时空离散，也能根据 PDE 的传播机制设计传输条件。优化后的 Robin、Ventcel 或卷积条件试图逼近连续问题的 Dirichlet-to-Neumann（DtN）映射：DtN 正是把一个子域外部对界面的全部影响浓缩成的算子，若传输条件与它完全一致，迭代一步即精确收敛，因此传输条件越接近 DtN，收敛越快。OSWR 已被用于众多 PDE，例如浅水方程（Martin 2009）、时域 Maxwell（Courvoisier and Gander 2013）、Schrödinger（Halpern and Szeftel 2010；Besse and Xing 2017；Antoine and Lorin 2017）、海洋原始方程（Audusse, Dreyfuss and Merlet 2010）、量子波问题（Antoine and Lorin 2016）、分数扩散（Wu 2017）以及耦合 Ekman 边界层（Thery et al. 2022）等。

![Schwarz 波形松弛在完整时间窗上交换界面波形](assets/diagrams/pint/zh/schwarz-waveform-relaxation.svg)

OSWR 对一阶抛物问题（如对流扩散方程 (2.5) 和非线性 Burgers 方程 (2.6)）与二阶双曲问题（如波动方程 (2.7)）的收敛特性截然不同，下面分别讨论。

## 3.2.1 一阶抛物问题

### OSWR 的完整迭代

考虑区间 $(0,L)$ 上带齐次 Dirichlet 边界 $u(0,t)=u(L,t)=0$ 和初值 $u(x,0)=u_0(x)$ 的对流扩散方程，记

$$
\mathcal L=\partial_x-\nu\partial_{xx}.
$$

取重叠子域 $\Omega_1=(0,\beta L)$、$\Omega_2=(\alpha L,L)$，其中 $0<\alpha\leq\beta<1$。这里 $\alpha,\beta$ 是相对于区间长度 $L$ 的无量纲坐标，因此 $\beta-\alpha$ 是**以区间长度为单位的重叠比例**，而真正的**物理重叠宽度**为 $l=(\beta-\alpha)L$。论文的 Robin OSWR 写成

$$
\left\{
\begin{aligned}
\partial_tu_1^k+\mathcal Lu_1^k&=0,
&& (x,t)\in\Omega_1\times(0,T],\\
u_1^k(0,t)&=0,\\
\frac1p\partial_xu_1^k(\beta L,t)+u_1^k(\beta L,t)
&=\frac1p\partial_xu_2^{k-1}(\beta L,t)+u_2^{k-1}(\beta L,t),
\end{aligned}
\right.
$$

$$
\left\{
\begin{aligned}
\partial_tu_2^k+\mathcal Lu_2^k&=0,
&& (x,t)\in\Omega_2\times(0,T],\\
\frac1p\partial_xu_2^k(\alpha L,t)-u_2^k(\alpha L,t)
&=\frac1p\partial_xu_1^{k-1}(\alpha L,t)-u_1^{k-1}(\alpha L,t),\\
u_2^k(L,t)&=0.
\end{aligned}
\right. \tag{3.1}
$$

每一轮都恢复物理初值：$u_i^k(x,0)=u_0(x)$。$k=0$ 的两条界面波形 $\{u_1^0(\alpha L,t),u_2^0(\beta L,t)\}$ 可以任意给定。参数 $p>0$ 决定 Robin 传输强度：左子域在其右界面 $x=\beta L$、右子域在其左界面 $x=\alpha L$ 各用一条 Robin 条件把邻域上一轮的解（值加导数的组合）传进来。当 $p\to\infty$ 时，$\tfrac1p\partial_x$ 项消失，Robin 条件退化为经典 Dirichlet 交换 $u_1^k(\beta L,t)=u_2^{k-1}(\beta L,t)$、$u_2^k(\alpha L,t)=u_1^{k-1}(\alpha L,t)$。多子域推广沿每个人工边界重复相同构造，直接而无需额外理论；非线性版本保留迭代结构，只把 $\mathcal L$ 换成相应的非线性算子。

### Theorem 3.1：Robin 参数的 minimax 选择

Theorem 3.1（Gander and Halpern 2007）在**无界空间、两个子域、连续时空**的简化假设下分析 (3.1)。令 $l>0$ 为物理重叠宽度，

$$
y_0=\frac{l}{\nu},
\qquad
y=\frac{l\omega}{\nu},
\qquad
\omega\in\left[\frac\pi T,\frac\pi{\Delta t}\right].
$$

优化参数具有 $p^*=\widetilde p^*/\nu$ 的缩放（论文原式 $p^*=\widetilde p^*\nu/s$ 的等价写法）。定义单频收敛因子

$$
R_0(y,\widetilde p,y_0)
=\frac{(y-\widetilde p)^2+y^2-y_0^2}
{(y+\widetilde p)^2+y^2-y_0^2}e^{-y},
$$

以及内部极值位置

$$
\bar y(y_0,\widetilde p)
=\sqrt{\frac{
y_0^2+2\widetilde p
+\sqrt{\widetilde p\left(-\widetilde p^3-4\widetilde p^2
+(4+2y_0^2)\widetilde p+8y_0^2\right)}}{2}}.
$$

$R_0$ 是在 Fourier 空间得到的 OSWR 收敛因子，其中 $y$ 对应单个 Fourier 模 $\omega\in[\pi/T,\pi/\Delta t]$，即 $y=l\omega/\nu$。优化 $p$ 的目标是让 $R_0$ 在整个相关频率区间上的最大值尽量小，因此分成两种情形：

若 $y_0<y_c$，其中 $y_c=1.618386576\ldots$，则 $\widetilde p^*$ 是下式的唯一解：

$$
R_0(y_0,\widetilde p^*,y_0)
=R_0\!\left(\bar y(y_0,\widetilde p^*),\widetilde p^*,y_0\right). \tag{3.2a}
$$

若 $y_0\geq y_c$，则 $\widetilde p^*$ 由

$$
y_0=\widetilde p^*\sqrt{\frac{\widetilde p^*}{4+\widetilde p^*}} \tag{3.2b}
$$

唯一确定。令

$$
y_{\min}=\frac{l\pi}{\nu T},
\qquad
y_{\max}=\frac{l\pi}{\nu\Delta t},
$$

则用优化 Robin 参数 $p^*$ 时，所有相关频率上的最坏收敛因子满足

$$
\rho:=\max_{y\in[y_{\min},y_{\max}]}
R_0(y,\widetilde p^*,y_0)
\leq
R_0\!\left(\bar y(y_0,\widetilde p^*),\widetilde p^*,y_0\right). \tag{3.2c}
$$

经典 Dirichlet 交换对应 $p=\infty$。把 $p=\infty$ 代入 $R_0$（有理系数趋于 $1$），只剩下穿越重叠的指数衰减，于是

$$
\rho\leq e^{-y_{\min}}
=\exp\!\left(-\frac{l\pi}{\nu T}\right). \tag{3.3}
$$

> [!tip] 推导 intuition：为什么出现等峰条件与两个分支
> 对时间作 Fourier 变换后，每个频率 $\omega$ 都变成一个独立的空间界面误差传播问题。穿越重叠区给出 $e^{-y}$ 的衰减，Robin 条件给出前面那个有理反射系数 $\big[(y-\widetilde p)^2+y^2-y_0^2\big]/\big[(y+\widetilde p)^2+y^2-y_0^2\big]$。优化 $p$ 等价于压低整个频率区间上最高的峰（一个 minimax 问题）。当两个候选峰一样高时，继续降低其中一个必然抬高另一个，公式 (3.2a) 正是这个"等峰"（equioscillation）平衡条件：区间端点 $y_0$ 处的峰与内部极值 $\bar y$ 处的峰相等。$y_0=y_c\approx1.618$ 是一个结构分界点：当重叠相对扩散足够大（$y_0\geq y_c$）时，内部极值不再是有效约束，最优 $p$ 转由更简单的代数关系 (3.2b) 决定。式 (3.3) 则说明 Dirichlet 只能靠重叠衰减，因子随 $l/(\nu T)$ 变化，明显弱于 Robin 情形。

### Figure 3.1：连续理论与离散实验的距离

![原论文 Figure 3.1：OSWR 理论收敛因子和四子域迭代数](assets/papers/time-parallelization/source-figures/figure-3-1.svg)

论文使用 $L=8.2$、$T=5$、$\Delta t=0.01$、$\Delta x=0.02$、$l=2\Delta x$，空间用中心差分，时间用后向 Euler，初值为

$$
u_0(x)=e^{-10(x-L/2)^2}.
$$

Figure 3.1(a) 画出 Dirichlet 和优化 Robin 条件的理论收敛因子随扩散参数 $\nu$ 的变化。$\nu$ 越小，对流项越占主导，因子越小，方法收敛越快。Figure 3.1(b) 把 $(0,L)$ 分成四个子域，从随机界面猜测出发，当迭代解与收敛解之差小于 $10^{-8}$ 时停止。实测迭代数同样随 $\nu$ 减小而下降，与 (a) 的理论预测吻合良好。

在 $\nu=0.1$ 时，Dirichlet 与优化 Robin 实测分别需要 92 和 28 轮；而二子域连续理论给出的 $\rho$ 预测为 32 和 4 轮，明显更小。差异来自三项设置变化：理论使用**无界域**，实验使用有界域；理论只有**两个子域**，实验有四个；实验还离散了空间和时间。论文据此提醒，Figure 3.1(a) 描述的是趋势和理想因子，不能直接充当多子域离散实现的迭代数公式。

关于多子域的理论现状需要分两层说清楚：**Dirichlet 传输条件的多子域收敛分析是存在的**，见 Gander and Stuart (1998) 与 Wu, Huang and Huang (2012)；**而 Robin 传输条件在多子域情形下的完整收敛分析目前仍然缺失**。对 Robin 条件，两个子域的半离散收敛分析可见 Wu and Al-Khaleel (2014)；稳态情形下连续与离散、有界与无界域之间的详细对比见 Gander, Halpern, Hubert and Krell (2020, 2021a)。

> [!tip] 本站洞见：把 Figure 3.1 读对
> Figure 3.1(a) 是"两子域、无界、连续"的理想因子，Figure 3.1(b) 是"四子域、有界、离散"的真实迭代数，两者之间隔着三重简化。因此正确的读法是：理论曲线用于判断**趋势与相对优劣**（Robin 系统性地优于 Dirichlet、小 $\nu$ 更快），而**绝对迭代数必须以离散实验为准**。这也解释了为什么补上 Dirichlet 多子域理论（Gander and Stuart 1998；Wu, Huang and Huang 2012）之外，仍有必要指出 Robin 多子域理论的缺口——后者正是理论与 Figure 3.1(b) 之间最大的一段未闭合距离。

### 更精确的传输条件

除 Dirichlet 与 Robin 外，还可以用 Ventcel 条件进一步加速 OSWR（Bennequin et al. 2016）。Ventcel 条件本质上是对最优传输条件的局部高阶近似——后者在 Fourier（或 Laplace）空间由 Gander and Halpern (2007) 给出，形如

$$
\partial_x-\frac{1}{2\nu}
\mathcal F^{-1}\!\left(1+\sqrt{1+4i\nu\omega}\right),
$$

其中 $i^2=-1$，$\mathcal F^{-1}$ 是关于 Fourier 模 $\omega$ 的逆变换。这个算子非局部（含平方根、耦合所有频率），无法直接实现，故 Ventcel 用局部微分算子去逼近它。在渐近意义下，若取 $l=C_1\Delta x$、$\Delta t=C_1\Delta x^\beta$ 且 $\Delta x$ 很小，收敛因子满足 $\rho=1-O(\Delta x^\gamma)$，指数 $\gamma>0$ 由 $\beta$ 决定（Gander and Halpern 2007；Bennequin et al. 2009, 2016）。这意味着网格越细因子越接近 $1$，收敛随之变慢。Wu and Xu (2017) 的卷积传输条件则给出**与网格无关**的常数因子 $\rho=1-C$，$C\in(0,1)$，从而避免了这种随网格加密的退化；它尤其适合处理带非局部项的演化方程，如 Volterra 型偏积分微分方程。代价是卷积接口算子需要保存时间历史。

> [!tip] 本站洞见：三档传输条件的取舍
> 从 Dirichlet 到 Robin 到 Ventcel/卷积，是一条"越来越逼近 DtN、也越来越贵"的谱系。Dirichlet 最简单但只靠重叠衰减；Robin 用一个可优化标量 $p$ 换取显著加速；Ventcel 用局部高阶算子进一步逼近最优算子，但在网格加密时因子仍趋于 $1$；卷积条件以保存时间历史为代价换来网格无关的常数因子。选择哪一档，取决于能承受多少接口复杂度与内存，以及网格是否会显著加密。

## 3.2.2 二阶双曲问题

### Dirichlet 条件也能有限步收敛

与一阶抛物问题不同，对二阶双曲问题（如波动方程 (2.7)），即便只用简单的 Dirichlet 传输条件，SWR 也能在**有限步内**收敛。两个重叠子域上的 SWR 为

$$
\left\{
\begin{aligned}
\partial_{tt}u_1^k&=c^2\partial_{xx}u_1^k+g,
&& (x,t)\in\Omega_1\times(0,T],\\
u_1^k(x,0)&=u_0(x),
&\partial_tu_1^k(x,0)&=\widetilde u_0(x),\\
u_1^k(0,t)&=0,
&u_1^k(\beta L,t)&=u_2^{k-1}(\beta L,t),
\end{aligned}
\right.
$$

$$
\left\{
\begin{aligned}
\partial_{tt}u_2^k&=c^2\partial_{xx}u_2^k+g,
&& (x,t)\in\Omega_2\times(0,T],\\
u_2^k(x,0)&=u_0(x),
&\partial_tu_2^k(x,0)&=\widetilde u_0(x),\\
u_2^k(\alpha L,t)&=u_1^{k-1}(\alpha L,t),
&u_2^k(L,t)&=0.
\end{aligned}
\right. \tag{3.4}
$$

这里 $c>0$ 是波速，$0<\alpha<\beta<1$。两子域各自恢复物理初值 $u_0,\widetilde u_0$，仅在人工界面 $x=\beta L$、$x=\alpha L$ 上用 Dirichlet 方式交换上一轮邻域的界面值。

**Theorem 3.2（Gander 1997, Theorem 6.3.3）.** 对 $0<\alpha<\beta<1$，只要

$$
k>\frac{Tc}{\beta-\alpha},
$$

第 $k$ 轮在两个人工界面上的误差便在整个 $(0,T)$ 内为零：

$$
u_1^k(\alpha L,t)-u(\alpha L,t)=0,
\qquad
u_2^k(\beta L,t)-u(\beta L,t)=0.
$$

该结论来自双曲问题固有的有限传播速度。误差本身满足零源、零初值的波动方程，故误差只能以速度 $c$ 沿特征线传播。每一轮界面交换，只能把"已知正确"的区域沿特征锥向内推进一段物理宽度 $l=(\beta-\alpha)L$；这段宽度对应的时间高度约为 $(\beta-\alpha)L/c$。经过 $k$ 轮，正确区域累计推进的时间高度约为 $k(\beta-\alpha)L/c$。注意此处的 $\beta-\alpha$ 是重叠占区间长度的比例，$T$ 已按同一无量纲化处理，于是当累计高度超过整个时间窗时，条件恰为 $k>Tc/(\beta-\alpha)$，此后整段界面波形都已精确。多子域、高维分解（Gander and Halpern 2004）以及一维非线性守恒律（Gander and Rohde 2005）也有相应结果。

> [!tip] 本站洞见：抛物与双曲的分水岭就在传播速度
> Theorem 3.1 与 Theorem 3.2 的对照是全章的枢纽。抛物问题传播速度无限，任何界面误差瞬间影响整个子域，故只能得到"每轮乘一个 $<1$ 因子"的线性收敛，且因子依赖 $\nu,l,T$。双曲问题传播速度有限，误差被锁在特征锥内，故每轮"精确"地扩大一块正确区域，得到有限步收敛，且步数只依赖几何量 $Tc/(\beta-\alpha)$。重叠比例 $\beta-\alpha$ 越大，每轮推进越多、所需步数越少——这直接启发了下面用"大重叠"换取并行的红黑构造。

### Figure 3.2：红黑 SWR 如何推进正确区域

![原论文 Figure 3.2：带大重叠的红黑 SWR 四阶段几何过程](assets/papers/time-parallelization/source-figures/figure-3-2.svg)

Theorem 3.2 揭示了一条可利用的性质（早在 Gander et al. 2003, Figure 3.1 已指出）：每个子域在"仅受初值影响、尚未被可能错误的界面数据污染"的特征锥内，算出的就是精确解。据此可以刻意选择时空子域，用**直接并行的子域求解**代替反复迭代。Figure 3.2 采用五个交错子域、并配以慷慨的重叠，计算顺序如下。

1. 在 $\Omega_1,\Omega_3,\Omega_5$ 三个红色子域上并行求解到 $T_1$。位于 $x_2,x_4$ 的内部界面此时未知，可先给任意值。由有限传播速度，每个子域底部由特征线界定的三角形 tent 已经精确；两端子域因外边界已知，还各多出一小块精确区域。此步在正确 tent 上方也算出了"尚不正确"的近似——这正是 Nievergelt 为换取并行而主张的冗余计算。
2. 在 $\Omega_2,\Omega_4$ 两个黑色子域上并行求解到 $T_2$。由于上一轮蓝色区域已正确，它们据此获得正确界面数据，从而生成两个正确的菱形 tent，同样附带一些冗余计算。
3. 再回到红色子域，把时间推进到区间 $(T_1,T_3)$。
4. 红黑交替继续进行，精确区域沿时间不断上推，最终覆盖到 $T_4$。

每个子域在正确 tent 上方还计算了暂时不可靠的区域，这些冗余工作换来了并行性。Nievergelt 的思想在这里以特征锥的几何方式重现：宁可多算一些注定要被丢弃的区域，也要让各子域在同一步内彼此独立地推进。

### MTP、UTP 与 Figure 3.3

上述红黑 SWR 其实是实现当前最强大的双曲时空求解器之一——mapped tent pitching（MTP，Gopalakrishnan, Schöberl and Wintersteiger 2017；时域 Maxwell 应用见 Gopalakrishnan et al. 2020）——的一种简单方式。MTP 把红黑 SWR 中出现的倾斜 tent 映射成时空柱体，在柱体内用经典时间步进求解，再映回原几何，从而**避免冗余计算**。代价有二：一是需要计算映射本身；二是映射后计算域尺寸与红黑 SWR 中的时空子域相当，计算量可比。更麻烦的是，映射会引入**阶数下降**（order reduction），需要专门设计的时间积分器来补偿。

相反，红黑 SWR（现也称 unmapped tent pitching，UTP）保留原坐标，因而**不发生阶数下降**，而且可以非常简单地实现：它等价于把 DD 里的限制加性 Schwarz（RAS）技术直接施加到全时间"all-at-once"时空系统上（解释见 Gander 2008），这一等价关系让高维实现可以复用成熟的 RAS 基础设施。

![原论文 Figure 3.3：UTP 对二阶波动方程误差区域的逐轮消除](assets/papers/time-parallelization/source-figures/figure-3-3.svg)

Figure 3.3 展示用红黑 SWR / UTP 求解波动方程模型问题 (2.7)。(a) 是随机界面猜测产生的初始误差；(b)、(c)、(d) 分别显示第 4 次红更新、第 8 次黑更新和第 12 次红更新。可以看到，UTP 无需预先知道 tent 结构，就在红、黑 tent 内构造出精确解，推进方式与 MTP 完全一致。蓝色零误差区沿 tent 逐层上升。UTP 同样可轻松用于非线性双曲问题；若事先不知道 tent 高度，只需观察计算解的残差——残差在时间方向降为零的高度自然指示了当前 tent 高度，据此即可自适应地选取时间窗长度 $T_i-T_{i-1}$。

原始 MTP 依赖有限传播速度，因此**无法直接用于抛物问题**：抛物问题传播速度无限，不存在"解在其中精确"的 tent。不过 SWR 以及 UTP 仍可用于抛物方程，其中优化 SWR 变体（Gander and Halpern 2007；Bennequin et al. 2009）尤其有效。对弱扩散的对流主导问题（如本文的对流扩散模型），可以考虑照样使用 UTP，只是在每个时间 slab 里多加一两轮迭代来修正跨 tent 的影响。

> [!tip] 本站洞见：UTP=对全时间系统做 RAS
> UTP 最实用的一点，是它把"tent pitching"这种看似特殊的时空推进，还原为一个人们已经非常熟悉的对象——对 all-at-once 时空系统做限制加性 Schwarz。这带来三重好处：无需显式构造映射、不产生阶数下降、可直接借用现成 RAS 代码扩展到高维。再配合"用残差降零高度自适应 tent 高度"的技巧，UTP 既保留了 MTP 的精确推进结构，又回避了它的映射成本与阶数问题；对弱扩散抛物问题，只需每 slab 多迭代一两轮即可迁移使用。

## 公式与图表覆盖核对

| 原文项目                          | 论文小节 | 覆盖状态                                                                    |
| --------------------------------- | -------- | --------------------------------------------------------------------------- |
| 历史发展，pp. 396–398             | 3.1      | 四类方法的来源、主要分支与关键文献均已列出，含 DN/NN 变体                   |
| WR 分裂与 SWR 动机                | 3.2 导论 | 连续区域分解限制、$A=M+N$ 分裂并行性、坏分裂困境、传输条件绕过难题          |
| (3.1)                             | 3.2.1    | 两个子域的 PDE、物理边界、Robin 交换、$p\to\infty$ 极限与初值               |
| (3.2a)–(3.2c), (3.3), Theorem 3.1 | 3.2.1    | 参数缩放、单频因子、$y_c$ 两分支、minimax 等峰、最坏频率界与 Dirichlet 极限 |
| Figure 3.1                        | 3.2.1    | 完整原图、参数、停止准则、三重差异，及 Dirichlet/Robin 多子域理论现状       |
| Ventcel 与卷积条件                | 3.2.1    | 最优 DtN 算子形式、渐近因子、网格无关卷积因子与非局部代价                   |
| (3.4), Theorem 3.2                | 3.2.2    | 双子域迭代、有限步条件 $k>Tc/(\beta-\alpha)$、有限传播速度论证              |
| Figures 3.2–3.3                   | 3.2.2    | 完整原图、红黑推进、Nievergelt 冗余、MTP/UTP、RAS 等价与残差自适应          |

## 本页原文

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Sections 3.1–3.2, pp. 396–405.
