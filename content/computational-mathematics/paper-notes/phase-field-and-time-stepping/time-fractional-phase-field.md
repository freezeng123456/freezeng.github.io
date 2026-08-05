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

编号 40、43、57 构成一条清晰的推进线：编号 40 证明积分型能量律并明确指出它不能被加强，编号 43 在非均匀网格上补上二阶保最大值原理的格式，编号 57 换掉能量对象、换来一条微分型的律。

## 40：分数阶能量耗散律是积分型的

### 直觉

相场模型（Allen-Cahn、Cahn-Hilliard、分子束外延）都来自一个随时间耗散的自由能，这条耗散律是长时间数值稳定性的支柱。整数阶情形的论证只有一行：用 $\phi_t$ 去测试方程，右端变成 $-\|\phi_t\|^2$，**逐点**为负，于是 $\mathrm dE/\mathrm dt\le0$。

把时间导数换成 Caputo 导数以后，这一行立刻断掉。$\partial_t^\alpha\phi$ 是整个历史的加权卷积，$\int_\Omega\phi_t\,\partial_t^\alpha\phi$ 在某个时刻可以是任意符号的：当前时刻的变化率与过去的变化率通过核 $(t-s)^{-\alpha}$ 耦合，没有理由逐点同号。在本文之前，时间分数阶相场模型的能量耗散只在数值上被观察到，连续与离散层面都没有证明。

本文的机制是**把符号性从逐点搬到累积**。核 $|t-s|^{-\alpha}$ 作为 $[0,T]^2$ 上的二次型是正定的——这是 Riesz 位势正定性的一个初等实例，可以由 Fourier 变换看出来。因此虽然 $\phi_t\cdot\partial_t^\alpha\phi$ 在单个时刻无符号，把它在 $[0,T]$ 上积一次，得到的双重积分是一个正量。**代价是所得的律只能是积分型的：它给出 $E[\phi(T)]\le E[\phi(0)]$，但不给出任何逐点的下降速率。** 这不是技术上的遗憾，而是分数阶梯度流的真实面貌，编号 57 的数值实验会直接把它显示出来。

> [!warning] 这条律不能被加强着引用
> 论文明确指出：所得的积分型能量律**不蕴含** $\frac{\mathrm d}{\mathrm dt}E\le0$，也**不蕴含** $\frac{\mathrm d^{\alpha}}{\mathrm dt^{\alpha}}E\le0$。这是下游文献中最常被误述的一点——把这里的结论转述成逐点的（整数阶或分数阶）微分耗散律并不成立，论文在结论部分把这两条列为公开问题。正因为加强不成立，编号 57 才必须**更换能量对象**，而不是把这条律推得更强。

### 问题设定

Caputo 导数取标准定义

$$
\frac{\partial^{\alpha}}{\partial t^{\alpha}}\phi
={}^{C}_{0}D_{t}^{\alpha}\phi(t)
:=\frac{1}{\Gamma(1-\alpha)}\int_{0}^{t}\frac{\phi'(s)}{(t-s)^{\alpha}}\,\mathrm ds,
\qquad t>0,\ \alpha\in(0,1).
$$

三个模型分别为（$\gamma$ 为迁移率，$\varepsilon$ 为界面厚度，边界条件取齐次 Dirichlet、Neumann 或周期）

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}
=\gamma\Bigl(\varepsilon\Delta\phi-\frac{1}{\varepsilon}F'(\phi)\Bigr)
\qquad\text{（Allen-Cahn）},
$$

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}=\gamma\Delta\mu,
\quad \mu=-\varepsilon\Delta\phi+\frac{1}{\varepsilon}F'(\phi)
\qquad\text{（Cahn-Hilliard，取周期或无流条件 }\partial_n\mu=\partial_n\phi=0\text{）},
$$

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}
=\gamma\Bigl(-\varepsilon\Delta^{2}\phi
+\frac{1}{\varepsilon}\nabla\cdot\mathbf f_m(\nabla\phi)\Bigr)
\qquad\text{（分子束外延，}\mathbf f_m=\partial F_m/\partial\mathbf v\text{）},
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

其中 $\langle\cdot,\cdot\rangle$ 与 $\|\cdot\|$ 是 $L^2(\Omega)$ 的内积与范数。分子束外延模型的 $F_m$ 按是否带斜率选择分两种：

$$
F_{m}(\mathbf v)=
\begin{cases}
\tfrac14\bigl(|\mathbf v|^{2}-1\bigr)^{2}, & \text{带斜率选择},\\[2pt]
-\tfrac12\ln\bigl|1+|\mathbf v|^{2}\bigr|, & \text{无斜率选择},
\end{cases}
$$

相应地 $\mathbf f_m(\nabla\phi)=(|\nabla\phi|^{2}-1)\nabla\phi$ 与 $-\nabla\phi/(1+|\nabla\phi|^{2})$。这两种模型在本专题里分属不同的页：无斜率选择的那一支由[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 52]] 处理，带斜率选择的那一支由[[computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation|编号 78]] 处理。

### 推导

**第一步：连续核的正定性。** 对 $h,g\in L^{p}(0,T)$ 定义

$$
I_\alpha(h,g):=\frac{1}{\Gamma(\alpha)}\int_{0}^{T}\!\!\int_{0}^{t}
\frac{h(s)g(t)}{(t-s)^{1-\alpha}}\,\mathrm ds\,\mathrm dt .
$$

论文证明的引理是：对 $\alpha\in(0,1)$ 与 $p\ge\frac{2}{1+\alpha}$，

$$
I_\alpha(h,h)=\int_0^T\bigl(I_{0+}^{\alpha}h\bigr)(t)\,h(t)\,\mathrm dt
\ \ge\ \cos\frac{\alpha\pi}{2}\,
\bigl\|I_{0+}^{\alpha/2}h\bigr\|^{2}_{L^2(0,T)}\ \ge\ 0 .
$$

真正被用到的是它在 $\alpha\mapsto1-\alpha$ 处的推论。令 $A_\alpha(h,g):=I_{1-\alpha}(h,g)$，则对 $p\ge\frac{2}{2-\alpha}$，

$$
A_\alpha(h,h)
=\frac{1}{2}\cdot\frac{1}{\Gamma(1-\alpha)}\int_{0}^{T}\!\!\int_{0}^{T}
\frac{h(s)h(t)}{|t-s|^{\alpha}}\,\mathrm ds\,\mathrm dt
\ \ge\ \sin\frac{\alpha\pi}{2}\,
\bigl\|I_{0+}^{(1-\alpha)/2}h\bigr\|^{2}_{L^2(0,T)}\ \ge\ 0 .
$$

**第一行的对称化恒等式是关键的代数一步**：Volterra 型的三角形区域二重积分被换成整个正方形上带 $|t-s|^{-\alpha}$ 的对称（Riesz 型）二重积分，这样才谈得上「二次型」。

证明路线是把 $h$ 从 $[0,T]$ 零延拓到 $\mathbb R$，用 Liouville 半群性质 $I_{+}^{\alpha}I_{+}^{\beta}=I_{+}^{\alpha+\beta}$ 把被积函数写成 $(I_{+}^{\alpha/2}h)(I_{-}^{\alpha/2}h)$，再用

$$
\bigl(\mathcal F I_{\pm}^{\alpha}f\bigr)(\xi)=\frac{(\mathcal Ff)(\xi)}{(\mp i\xi)^{\alpha}},
\qquad
(\mp i\xi)^{\alpha}=|\xi|^{\alpha}e^{\mp i\alpha\pi\,\mathrm{sgn}(\xi)/2}
$$

与 Parseval 等式，最后用 Hardy-Littlewood 不等式 $\|I_{0+}^\alpha f\|_{L^q}\le K\|f\|_{L^p}$（$q=p/(1-\alpha p)$）保证可积性。**常数 $\cos\frac{\alpha\pi}{2}$ 与 $\sin\frac{\alpha\pi}{2}$ 就是 $e^{\mp i\alpha\pi/2}$ 的实部**，这也解释了为什么 $\alpha\to1$ 时前一条退化（$\cos\frac\pi2=0$）而后一条趋于 $1$。

**第二步：连续能量律。** 用 $\phi_t$ 测试方程并在 $[0,T]$ 上积分，右端的椭圆项与势能项按通常方式合成 $E[\phi(T)]-E[\phi(0)]$，左端 $\frac1\gamma\int_\Omega\int_0^T\phi_t\,\partial_t^\alpha\phi$ 恰好是 $\frac1\gamma\int_\Omega A_\alpha(\phi_t,\phi_t)\,\mathrm dx$，由上一步非负。

**第三步：L1 离散与离散核的完全单调性。** 取均匀步长 $\tau=T/n$，L1（分段线性 Caputo）核为

$$
b_{j}=\frac{1}{\Gamma(1-\alpha)}\int_{j\tau}^{(j+1)\tau}\frac{\mathrm dt}{t^{\alpha}}
=\frac{\tau^{1-\alpha}}{\Gamma(2-\alpha)}
\bigl[(j+1)^{1-\alpha}-j^{1-\alpha}\bigr],\qquad j\ge0 .
$$

（积分限处的 $\tau$ 在原文中误排为 $\pi$；右端的闭式无歧义，实际使用的是它。）这些核满足

$$
b_{k}>0,\qquad b_{k}-b_{k+1}>0,\qquad
\sum_{j=0}^{k-1}(b_{j}-b_{j+1})+b_{k}=b_{0},
$$

以及重排恒等式

$$
\sum_{j=0}^{k}b_{j}\frac{u^{k+1-j}-u^{k-j}}{\tau}
=\frac{1}{\tau}\Bigl[b_{0}u^{k+1}
-\sum_{j=0}^{k-1}(b_{j}-b_{j+1})u^{k-j}-b_{k}u^{0}\Bigr].
$$

**这条重排恒等式是全篇离散论证的枢纽**：它把「当前层的系数」与「历史层的非负权重」分开，历史部分的权重 $b_j-b_{j+1}$ 与 $b_k$ 全为正，且总和恰为 $b_0$。稳定化 L1 格式（以 Allen-Cahn 为例，$f=F'$，$S>0$ 为稳定化常数）为

$$
\frac{1}{\gamma}\sum_{j=0}^{k}b_{j}\frac{\phi^{k+1-j}-\phi^{k-j}}{\tau}
=\varepsilon\Delta\phi^{k+1}-\frac{1}{\varepsilon}f(\phi^{k})
-\frac{S}{\gamma}\bigl(\phi^{k+1}-\phi^{k}\bigr),
$$

Cahn-Hilliard 与分子束外延格式同构，只是把稳定化项分别放进化学势与作用在 $\Delta$ 上。

**第四步：离散二次型的正定性。** 这是上面连续推论的离散对应物。对任意 $(u_1,\dots,u_n)^T\in\mathbb R^{n}$ 令

$$
B:=2\sum_{k=1}^{n}\sum_{j=1}^{k}b_{|k-j|}u_{j}u_{k}
=\sum_{k=1}^{n}b_{0}u_{k}^{2}+\sum_{k=1}^{n}\sum_{j=1}^{n}b_{|k-j|}u_{j}u_{k}
\ \ge\ \sum_{k=1}^{n}b_{0}u_{k}^{2},
$$

更精细地

$$
B\ \ge\ \frac{2}{\tau}\,\sin\frac{\alpha\pi}{2}\,
\bigl\|I_{0+}^{(1-\alpha)/2}u^{n}(t)\bigr\|^{2}_{L^{2}(0,T)}
+s_{n}\sum_{k=1}^{n}u_{k}^{2},
\qquad
s_{n}=\Bigl(\frac{n+1}{2}\Bigr)^{-\alpha}\frac{\tau^{1-\alpha}}{\Gamma(1-\alpha)}>0,
$$

其中 $u^{n}(t)=u_{\lfloor t/\tau\rfloor+1}$。证明的办法是把 $b_{|k|}$ 与精确积分核

$$
\tilde b_{|k|}=\frac{\tau^{1-\alpha}}{\Gamma(3-\alpha)}
\bigl((k+1)^{2-\alpha}-2k^{2-\alpha}+(k-1)^{2-\alpha}\bigr)\ (k\ge1),
\qquad
\tilde b_{0}=\frac{2\tau^{1-\alpha}}{\Gamma(3-\alpha)}
$$

作比较：$\tilde b_{|k|}-b_{|k|}\ge0$（$k\ge1$）、$2b_{0}-\tilde b_{0}\ge0$，于是差矩阵 $C=\{b_{|k-j|}-\tilde b_{|k-j|}\}$ 是对称正定 $M$-矩阵，对角占优给出上述下界，且 $s_n=2b_0-\tilde b_0-c_0$。

**第五步：最大值原理的论证是比较/归纳，不是能量论证。** 值得强调，因为后来的 Liao-Tang-Zhou 系列都改用 DOC 核，而本文没有。假设势满足 $F\in C^2(\mathbb R)$ 且存在 $M_1<0<M_2$ 使

$$
F'(M_1)=F'(M_2)=0;\qquad F'(u)>0\ (u>M_2);\qquad F'(u)<0\ (u<M_1),
$$

四次双阱以 $M_1=-1$、$M_2=1$ 满足之。用重排恒等式把格式写成

$$
\Bigl(\frac{b_0}{\gamma\tau}+\frac{S}{\gamma}\Bigr)\phi^{k+1}-\varepsilon\Delta\phi^{k+1}
=\frac{1}{\gamma\tau}\Bigl[\sum_{j=0}^{k-1}(b_{j}-b_{j+1})\phi^{k-j}+b_{k}\phi^{0}\Bigr]
+\frac{S}{\gamma}\phi^{k}-\frac{1}{\varepsilon}f(\phi^{k}),
$$

再写 $f(\phi^k)=f(\phi^k)-f(M_2)=f'(\xi)(\phi^k-M_2)$，对算子 $\bigl(\frac{b_0}{\gamma\tau}+\frac{S}{\gamma}\bigr)I-\varepsilon\Delta$ 用椭圆最大值原理即可。**历史项系数非负是这一步能走通的唯一原因**，而这正是核的完全单调性。

### 定理

**（连续，时间分数阶 Allen-Cahn）** 若 $E[\phi(0)]$ 有限，则

$$
E[\phi(T)]-E[\phi(0)]
=-\frac{1}{\gamma}\int_{\Omega}A_{\alpha}(\phi_{t},\phi_{t})\,\mathrm dx\ \le\ 0 .
$$

注意这是一个**等式**加一个符号判断。它同时给出耗散量 $\frac1\gamma\int_\Omega A_\alpha(\phi_t,\phi_t)\mathrm dx$ 被 $E[\phi(0)]$ 控制，因而 $\|I_{0+}^{(1-\alpha)/2}\phi_t\|^2_{L^2(0,T)}$ 有界——这是这条律给出的、比「能量不增」更强的信息。

**（连续，时间分数阶 Cahn-Hilliard）** 先证总质量守恒，再得

$$
E[\phi(T)]-E[\phi(0)]
=-\frac{1}{\gamma}\int_{\Omega}A_{\alpha}(\nabla\psi,\nabla\psi)\,\mathrm dx\le0,
\qquad -\Delta\psi=\phi_{t}.
$$

**（连续，时间分数阶分子束外延）** 印刷形式为 $E_{m}[\phi(T)]-E_{m}[\phi(0)]\le-\frac1\gamma\int_\Omega A_\alpha(\phi,\phi)\,\mathrm dx\le0$。

> [!warning] 一处疑似排印
> 分子束外延这一条里 $A_\alpha$ 的两个宗量印作 $(\phi,\phi)$；按 Allen-Cahn 与 Cahn-Hilliard 两条的类比以及证明的结构，应为 $(\phi_t,\phi_t)$。引用这条显式时按后者理解。

**（离散，分数阶 Allen-Cahn）** 设势被改造成二阶导全局有界，$\max_u|F''(u)|\le L$；论文用的是二次增长截断

$$
F(\phi)=\begin{cases}
\tfrac{11}{2}(\phi-2)^{2}+6(\phi-2)+\tfrac94, & \phi>2,\\
\tfrac14(\phi^{2}-1)^{2}, & \phi\in[-2,2],\\
\tfrac{11}{2}(\phi+2)^{2}-6(\phi+2)+\tfrac94, & \phi<-2 .
\end{cases}
$$

若

$$
S+\frac{b_{0}}{2\tau}\ \ge\ \frac{\gamma L}{2\varepsilon},
$$

则以 $\delta_{t}\phi^{k+1}:=\phi^{k+1}-\phi^{k}$ 记，

$$
E[\phi^{n}]-E[\phi^{0}]
\le-\frac{b_{0}}{2\gamma\tau}\sum_{k=0}^{n-1}\|\delta_{t}\phi^{k+1}\|^{2}
-\sum_{k=0}^{n-1}\Bigl\{\frac{\varepsilon}{2}\|\nabla\delta_{t}\phi^{k+1}\|^{2}
+\Bigl\langle\frac{S}{\gamma}-\frac{1}{2\varepsilon}f'(\xi^{k}),
(\delta_{t}\phi^{k+1})^{2}\Bigr\rangle\Bigr\}.
$$

特别地，**取 $S\ge\gamma L/(2\varepsilon)$ 时格式无条件能量稳定**，对一切 $\tau>0$ 与 $n>0$ 有 $E[\phi^{n}]\le E[\phi^{0}]$。

**（离散，分数阶 Cahn-Hilliard）** 同样的结构，但耗散项取 $H^{-1}$ 范数（反映 $H^{-1}$ 梯度流），条件为 $\sqrt{b_{0}\varepsilon/(\gamma\tau)}+S/\gamma\ge L/(2\varepsilon)$，$S\ge\gamma L/(2\varepsilon)$ 时无条件。

**（离散，分数阶分子束外延）** 条件为 $\sqrt{b_{0}\varepsilon/(\gamma\tau)}+S/\gamma\ge\frac{1}{2\varepsilon}\lambda_{\max}(\mathbf f_m'(\xi^{k}))$。对**无斜率选择**的模型可以把这个特征值算清楚：

$$
\mathbf f_m'(\mathbf v)=\frac{2\mathbf v^{2}-(|\mathbf v|^{2}+1)I}{(1+|\mathbf v|^{2})^{2}},
\qquad
\lambda_{\max}\bigl(\mathbf f_m'\bigr)\le\frac18,
$$

因此 **$S\ge\dfrac{\gamma}{16\varepsilon}$ 就足以保证无条件能量稳定**。这个 $1/8$ 与 $\gamma/(16\varepsilon)$ 是原文印出的精确常数。

**（离散最大值原理，分数阶 Allen-Cahn）** 设 $\phi_0\in C^0$ 且 $M_1\le\phi_0\le M_2$，用**未截断**的标准双阱势。若

$$
\frac{b_{0}-b_{1}}{\tau}+S\ \ge\ \gamma\,\frac{\max_{M_1\le u\le M_2}|f'(u)|}{\varepsilon},
$$

则 $M_{1}\le\phi^{k}(x)\le M_{2}$ 对一切 $k\ge1$ 成立。由于

$$
b_0-b_1=\frac{\tau^{1-\alpha}}{\Gamma(2-\alpha)}\bigl(2-2^{1-\alpha}\bigr),
$$

这是一个**温和但真实的步长限制**，不是无条件结论。它的价值在于：一旦最大值原理成立，$f$ 的全局 Lipschitz 假设就可以丢掉，上面那个截断势也就不再需要。

### 数值实验

空间离散是 $\Omega=[0,L_x]\times[0,L_y]$ 上的周期 Fourier-Galerkin。历史项若直接求和，第 $n$ 层要用到全部 $n$ 个历史值，总代价 $O(N^2)$；论文用 Jiang 等人的**指数和（SOE）快速算法**来算分数阶导数的历史部分，这是让长时间粗化模拟可行的实现前提。

| 算例 | 模型                               | 分数阶参数             | 初值与观察对象               |
| ---- | ---------------------------------- | ---------------------- | ---------------------------- |
| 1    | 时间分数阶 Allen-Cahn              | $\alpha=1,\,0.5,\,0.3$ | 随机初值；解的快照与离散能量 |
| 2    | 时间分数阶 Cahn-Hilliard           | $\alpha=1,\,0.5,\,0.3$ | 随机初值；离散能量与粗化速率 |
| 3    | 时间分数阶分子束外延（带斜率选择） | $\alpha=1,\,0.7,\,0.4$ | 随机初值；离散能量与粗化速率 |

两条结论。第一，三个模型的计算能量都**单调下降**，与上面的离散定理一致。第二——这是本文数值部分真正的产出——在粗化阶段，能量耗散率服从幂律，其渐近指数为

$$
-\frac{\alpha}{3}
$$

对时间分数阶 Cahn-Hilliard 与时间分数阶分子束外延**都**成立。$\alpha=1$ 时它退化为经典的 $-1/3$ 粗化律，这是一致性检验。

**这条 $-\alpha/3$ 律超出了论文自己的理论。** 作者在结论中直说它是经验观察，「需要严格的理论论证」。它也是本文数值部分与理论部分之间最大的落差：定理保证的是能量不增，而实验测到的是能量下降的**速率**，后者完全在定理的射程之外。

另有两处理论未覆盖之处值得记下：格式在时间上只有**一阶**、网格是**均匀**的，因此分数阶解在 $t=0$ 附近的初始奇性没有被处理，形式阶实际上拿不到；这两点正是编号 43 的出发点。

### 与其他论文的关系

这是分数阶这一支的**奠基论文**，文献中所称「Tang-Yu-Zhou 分数阶能量耗散律」即此。它的核正定性引理与其离散对应物（L1 二次型 $B$ 的正定性）是[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 74]] 变步长正定性理论的直接祖先。编号 43 与 57 攻击的正是本文自己承认的两处局限：编号 43 给出非均匀网格上的二阶保最大值原理格式，编号 57 给出变步长下同时保能量稳定性与最大值原理的格式，并把这里的积分型律升级为一条 $\alpha\to1$ 时退化回经典律的微分型律。整数阶那一支（编号 48、52、58、67、69、91、104）则是把同一套「卷积正定性」哲学搬到变步长 BDF 与隐显方法上。

## 43：非均匀步长上的二阶保最大值原理格式

### 直觉

编号 40 的格式在时间上是一阶的，网格是均匀的。这在分数阶问题上有两个具体的坏处。第一，时间分数阶 Allen-Cahn 的解有内禀的初始奇性 $u_t\sim\mathcal O(t^{\alpha-1})$，均匀网格会把形式阶打掉，必须在 $t=0$ 附近加密（分级网格）。第二，动力学跨越多个时间尺度——初期演化快、随后粗化极慢——长时间模拟必须放大步长。

但非均匀网格会破坏分数阶卷积核的经典正定性与单调性论证，而当时没有二阶的保最大值原理格式。本文的机制有两层。

外层是**选一个在非均匀网格上仍能保住核结构的公式**：Alikhanov（L2-$1_\sigma$）公式，把方程配置在离网点 $t_{n-\theta}$（取 $\theta=\alpha/2$）上，于是核 $A^{(n)}_{n-k}$ 依然正、依然对滞后指标单调递减，只是它带**两个**指标——上标记录当前层。两个指标意味着它不再是 Toeplitz 卷积，这正是变步长之所以困难的地方。

内层是**离散互补卷积（DCC）核**：构造一族 $P^{(n)}_{n-j}\ge0$ 使得

$$
\sum_{j=k}^{n}P^{(n)}_{n-j}A_{j-k}^{(j)}\equiv1 .
$$

它们像一个「反算子」，把 $\sum_kA^{(n)}_{n-k}\nabla_\tau v^k\le\cdots$ 这样的不等式两端同时卷一次，左端就塌回 $v^n-v^0$。**由此得到的收敛估计之所以尖锐，是因为局部截断误差本身也有卷积结构**：用 $P$ 去卷一列截断误差，得到的量比这列误差的最大值小得多。这一点在初始层上最明显——$|\Upsilon^{1}|$ 在 $\sigma=\alpha$ 时是 $\mathcal O(1)$ 的，而 $P_0^{(1)}|\Upsilon^1|$ 却是 $\mathcal O(\tau_1^{\sigma})$。

### 问题设定

二维时间分数阶 Allen-Cahn 方程，$\Omega=(0,L)^2$，周期边界：

$$
\partial_{t}^{\alpha}u=\varepsilon^{2}\Delta u-f(u),\qquad f(u)=u^{3}-u,
\qquad u(\mathbf x,0)=u_{0}(\mathbf x),
$$

Caputo 导数写成卷积形式 $(\partial_{t}^{\alpha}v)(t)=(\mathcal I_{t}^{1-\alpha}v')(t)$，其中 $\omega_{\mu}(t):=t^{\mu-1}/\Gamma(\mu)$。这套 $\omega_\mu$ 记号贯穿整个 Liao-Tang-Zhou 系列。从编号 40 继承的连续事实是能量律 $E(t)\le E(0)$ 与最大值原理 $|u|\le1$。

> [!warning] 步长比约定与 BDF 各篇互为倒数
> 本文用的是
>
> $$
> \rho_k:=\frac{\tau_k}{\tau_{k+1}},\qquad \rho:=\max_{k\ge1}\rho_k,
> $$
>
> 而[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 48、52]] 用的是 $r_k:=\tau_k/\tau_{k-1}$，因此 $\rho_k=1/r_{k+1}$。**跨篇引用门槛时混用这两套约定，是这条文献里最容易出错的地方。**

两条网格假设：

- **M1**：最大步长比 $\rho=7/4$。作者的解读是「可以使用一串递减步长，缩减因子低至 $4/7$」，而对**放大**步长**不加任何限制**。
- **M2**（只用于收敛性）：存在与网格无关的 $C_{1\gamma},C_{2\gamma}>0$ 使 $\tau_k\le\tau\min\{1,C_{1\gamma}t_k^{1-1/\gamma}\}$ 且 $t_k\le C_{2\gamma}t_{k-1}$。准均匀网格对应 $\gamma=1$，模型情形是分级网格 $t_k=T(k/N)^{\gamma}$。

### 推导

**第一步：Alikhanov 公式与它的两指标核。** 以 $\Pi_{1,k}$、$\Pi_{2,k}$ 记 $\{t_{k-1},t_k\}$ 与 $\{t_{k-1},t_k,t_{k+1}\}$ 上的线性、二次插值，

$$
(\partial_{\tau}^{\alpha}v)^{n-\theta}
:=\int_{t_{n-1}}^{t_{n-\theta}}\omega_{1-\alpha}(t_{n-\theta}-s)(\Pi_{1,n}v)'(s)\,\mathrm ds
+\sum_{k=1}^{n-1}\int_{t_{k-1}}^{t_k}\omega_{1-\alpha}(t_{n-\theta}-s)(\Pi_{2,k}v)'(s)\,\mathrm ds,
$$

展开后由两族系数

$$
a_{n-k}^{(n)}:=\frac{1}{\tau_{k}}\int_{t_{k-1}}^{\min\{t_{k},t_{n-\theta}\}}
\omega_{1-\alpha}(t_{n-\theta}-s)\,\mathrm ds,
\qquad
b_{n-k}^{(n)}:=\frac{2}{\tau_{k}(\tau_{k}+\tau_{k+1})}
\int_{t_{k-1}}^{t_{k}}\bigl(s-t_{k-\frac12}\bigr)\omega_{1-\alpha}(t_{n-\theta}-s)\,\mathrm ds
$$

组装成紧凑的卷积形式 $(\partial_{\tau}^{\alpha}v)^{n-\theta}=\sum_{k=1}^{n}A_{n-k}^{(n)}\nabla_{\tau}v^{k}$，

$$
A_{n-k}^{(n)}:=
\begin{cases}
a_{0}^{(n)}+\rho_{n-1}b_{1}^{(n)}, & k=n,\\
a_{n-k}^{(n)}+\rho_{k-1}b_{n-k+1}^{(n)}-b_{n-k}^{(n)}, & 2\le k\le n-1,\\
a_{n-1}^{(n)}-b_{n-1}^{(n)}, & k=1 .
\end{cases}
$$

全离散格式（空间中心差分 $D_h$，周期边界）为

$$
(\partial_{\tau}^{\alpha}u)^{n-\theta}=\varepsilon^{2}D_{h}u^{n-\theta}-f(u)^{n-\theta},
\qquad
f(u)^{n-\theta}:=\theta f(u^{n-1})+(1-\theta)f(u^{n}) .
$$

**第二步：三条核估计。** 在 M1 下，

1. $A^{(n)}_{0}\le\dfrac{24}{11\tau_{n}}\displaystyle\int_{t_{n-1}}^{t_{n}}\omega_{1-\alpha}(t_{n}-s)\,\mathrm ds$，并且各核有形如 $\dfrac{4}{11}$ 倍的下界；
2. 单调性：$A^{(n)}_{n-k-1}-A^{(n)}_{n-k}>0$（$1\le k\le n-1$）；
3. 首核支配次核：$\dfrac{1-2\theta}{1-\theta}A^{(n)}_{0}-A^{(n)}_{1}>0$（$n\ge2$），取 $\theta=\alpha/2$ 时前面的系数是 $\dfrac{2-2\alpha}{2-\alpha}$。

这三条比 Alikhanov 原来的均匀网格估计更强，它们蕴含分数阶 Grönwall 引理的两条假设成立，且带常数

$$
\pi_A=\frac{11}{4}.
$$

（第 1 条的下界常被印成 $\frac{4}{11\tau_n}\int_{t_{n-1}}^{t_n}$，与核的指标 $n-k$ 不匹配；按 Grönwall 引理的假设应为 $\frac{4}{11\tau_k}\int_{t_{k-1}}^{t_k}\omega_{1-\alpha}(t_n-s)\mathrm ds$。这是指标上的排印疏漏，常数 $4/11$ 本身与 $\pi_A=11/4$ 一致。）

**第三步：DCC 核与互补恒等式。**

$$
P_{0}^{(n)}:=\frac{1}{A_{0}^{(n)}},
\qquad
P_{n-j}^{(n)}:=\frac{1}{A_{0}^{(j)}}\sum_{k=j+1}^{n}
\bigl(A_{k-j-1}^{(k)}-A_{k-j}^{(k)}\bigr)P_{n-k}^{(n)},\quad 1\le j\le n-1 .
$$

（前因子在原文中印作 $1/p_0^{(j)}$，应为 $1/A_0^{(j)}$。）由核的单调性 $P^{(n)}_{n-j}\ge0$，并且

$$
\sum_{j=k}^{n}P^{(n)}_{n-j}A_{j-k}^{(j)}\equiv1,
\qquad
\sum_{j=1}^{n}P^{(n)}_{n-j}\,\omega_{1+m\alpha-\alpha}(t_{j})
\le\pi_{A}\,\omega_{1+m\alpha}(t_{n}),\quad m=0,1 .
$$

第二条是把 $P$ 作用在幂型函数上的定量估计，是 Grönwall 引理的燃料。

**第四步：离散分数阶 Grönwall 不等式。** 以 $E_\alpha$ 记 Mittag-Leffler 函数，$\lambda=\lambda_0+\lambda_1$，若最大步长满足

$$
\tau\le\frac{1}{\sqrt[\alpha]{2\Gamma(2-\alpha)\lambda\pi_{A}}},
$$

且非负序列满足 $\sum_{k=1}^{n}A_{n-k}^{(n)}\nabla_{\tau}v^{k}\le\lambda_{0}v^{n}+\lambda_{1}v^{n-1}+\xi^{n}+\eta^{n}$，则

$$
v^{n}\le2E_\alpha\bigl(2\max\{1,\rho\}\lambda\pi_{A}t_{n}^{\alpha}\bigr)
\Bigl(v^{0}+\Gamma(1-\alpha)\pi_{A}\max_{k\le n}\{t_{k}^{\alpha}\xi^{k}\}
+\pi_{A}\,\omega_{1+\alpha}(t_{n})\max_{k\le n}\eta^{k}\Bigr).
$$

**第五步：最大值原理是 $\ell^\infty$ 矩阵论证加归纳。** 三个成分：(a) $D_h$ 对称、负半定，且 $d_{ii}=-\max_i\sum_{j\ne i}|d_{ij}|$；(b) 一条立方引理——若 $B$ 满足同样的对角条件且 $A=aI-B$（$a>0$），则对 $c>0$，

$$
\|AV\|_{\infty}\ge a\|V\|_{\infty},
\qquad
\|AV+cV^{3}\|_{\infty}\ge a\|V\|_{\infty}+c\|V\|_{\infty}^{3};
$$

(c) 把格式整理成

$$
(A^{(n)}_{0}-1+\theta)u^{n}-(1-\theta)\varepsilon^{2}D_{h}u^{n}+(1-\theta)(u^n)^{.3}
=\mathcal{L}^{n-2}(u)+\cdots,
\qquad
\mathcal{L}^{n-2}(u):=\sum_{k=1}^{n-2}\bigl(A_{n-k-1}^{(n)}-A_{n-k}^{(n)}\bigr)u^{k}+A_{n-1}^{(n)}u^{0},
$$

历史算子的系数由单调性全为非负，因此归纳假设 $\|u^k\|_\infty\le1$（$k\le n-1$）能把右端压住，而左端由立方引理有下界。

### 定理

**（唯一可解性）** 若 M1 成立且 $\tau\le\sqrt[\alpha]{\omega_{2-\alpha}(1-\theta)/(1-\theta)}$，则非线性格式唯一可解。证明用 $A_0^{(n)}\ge a_0^{(n)}=\omega_{2-\alpha}(1-\theta)/\tau_n^{\alpha}\ge1-\theta$ 得到 $G_h:=A_0^{(n)}-1+\theta-(1-\theta)\varepsilon^{2}D_h$ 正定，$u^n$ 是严格凸泛函 $\frac12w^{T}G_hw+\frac{1-\theta}{4}\sum_kw_k^{4}-w^{T}g(u^{n-1})$ 的唯一极小点。

**（离散最大值原理）** 设 M1（$\rho=7/4$）且

$$
\tau\le\min\Bigl\{
\sqrt[\alpha]{\frac{\theta\,\omega_{2-\alpha}(1-\theta)}{2(1-\theta)}},\
\sqrt[\alpha]{\frac{h^{2}\,\omega_{2-\alpha}(1-\theta)}{4\varepsilon^{2}}}\Bigr\},
$$

则 $\|u^{0}\|_{\infty}\le1$ 蕴含 $\|u^{k}\|_{\infty}\le1$（$1\le k\le N$）。这是本文的招牌结论：「时间分数阶 Allen-Cahn 方程的第一个二阶保最大值原理格式」。**第二个约束把 $\tau$ 与空间网格 $h$ 及 $\varepsilon$ 耦合起来，是真实（虽然温和）的限制，不是无条件结果。**

**（尖锐的最大范数收敛）** 设 $\|u^0\|_{L^\infty}\le1$，正则性假设

$$
\|u(t)\|_{W^{4,\infty}(\Omega)}\le C_u,
\qquad
\|u^{(\ell)}(t)\|_{W^{2,\infty}(\Omega)}\le C_u(1+t^{\sigma-\ell})\ (\ell=1,2,3),\ \sigma\in(0,1),
$$

M1 成立，且 $\tau$ 满足上面两条再加 $\tau\le\sqrt[\alpha]{\omega_{2-\alpha}(1)/11}$，则

$$
\|u(t_{n})-u^{n}\|_{\infty}\le C_{u}\Bigl(\frac{\tau_{1}^{\sigma}}{\sigma}
+\frac{1}{1-\alpha}\max_{2\le k\le n}t_{k}^{\alpha}t_{k-1}^{\sigma-3}\tau_{k}^{3-\alpha}
+h^{2}\Bigr);
$$

若再补上 M2，则

$$
\|u(t_{n})-u^{n}\|_{\infty}\le\frac{C_{u}}{\sigma(1-\alpha)}\tau^{\min\{\gamma\sigma,\,2\}}+C_{u}h^{2},
$$

于是**最优的 $\mathcal O(\tau^{2})$ 阶在分级参数 $\gamma\ge\max\{1,2/\sigma\}$ 时达到**。证明用最大值原理来避免假设 $f$ 全局 Lipschitz。

**（误差的卷积结构）** 这是尖锐性的技术来源。在 M1 下，对 $v\in C^3((0,T])$ 且 $\int_0^Ts^2|v'''(s)|\mathrm ds<\infty$，局部相容误差 $\Upsilon^n[v]:=(\partial_t^\alpha v)(t_{n-\theta})-(\partial_\tau^\alpha v)^{n-\theta}$ 自身具有卷积结构

$$
|\Upsilon^{n}[v]|\le A_{0}^{(n)}G_{\mathrm{loc}}^{n}
+\sum_{k=1}^{n-1}\bigl(A_{n-k-1}^{(n)}-A_{n-k}^{(n)}\bigr)G_{\mathrm{his}}^{k},
$$

$$
G_{\mathrm{loc}}^{k}:=\frac{3}{2}\int_{t_{k-1}}^{t_{k-1/2}}(s-t_{k-1})^{2}|v'''(s)|\mathrm ds
+\frac{3\tau_{k}}{2}\int_{t_{k-1/2}}^{t_{k}}(t_{k}-s)|v'''(s)|\mathrm ds,
$$

$$
G_{\mathrm{his}}^{k}:=\frac{5}{2}\int_{t_{k-1}}^{t_{k}}(s-t_{k-1})^{2}|v'''(s)|\mathrm ds
+\frac{5}{2}\int_{t_{k}}^{t_{k+1}}(t_{k+1}-s)^{2}|v'''(s)|\mathrm ds .
$$

（$G_{\mathrm{loc}}$ 第二项的积分限在原文中写反，此处按显然的意图写正。）由于 $P$ 与 $A$ 互补，全局误差是 $\sum_jP_{n-j}^{(n)}|\Upsilon^j[v]|$，它远小于 $\max_j|\Upsilon^j[v]|$：在 $n=1$ 且 $\sigma=\alpha$ 时 $|\Upsilon^{1}|\le C_u\tau_1^{\sigma-\alpha}/\sigma$ 只是 $\mathcal O(1)$，而 $P_0^{(1)}|\Upsilon^1|\le G_{\mathrm{loc}}^1\le C_u\tau_1^{\sigma}/\sigma$。**这就是「借互补核实现超收敛」，也是估计能够反映时间正则性而不是最坏截断误差的原因。**

### 数值实验

实现上用基于指数和逼近的快速 Alikhanov 公式（SOE 绝对容差 $\epsilon=10^{-12}$），每层作非线性迭代（终止误差 $\eta=10^{-12}$），自适应步长选择器为

$$
\tau_{\mathrm{ada}}(e,\tau)=S_a\Bigl(\frac{tol}{e}\Bigr)^{1/2}\tau,
$$

其中 $e$ 是一阶（后向 Euler/L1）解与本文二阶解之间的相对差；并设一个**永不让步长缩到 $\frac23\tau_{n-1}$ 以下**的守卫。这个守卫与 M1 是配套的：$\rho\le7/4$ 允许的最大缩减是 $4/7\approx0.571$，$2/3$ 落在其内。

**算例一（精度检验，$\alpha=0.8$）。** 分别取时间正则性 $\sigma=0.8$ 与 $\sigma=0.4$，理论最优分级为 $\gamma_{\mathrm{opt}}=\max\{1,2/\sigma\}$：

| 时间正则性 $\sigma$ | 最优分级 $\gamma_{\mathrm{opt}}$ | 准均匀网格（$\gamma=1$）观测阶 | 分级网格观测阶      |
| ------------------- | -------------------------------- | ------------------------------ | ------------------- |
| $0.8$               | $2.5$                            | $\approx0.80$                  | $\approx2.1$–$2.34$ |
| $0.4$               | $5$                              | —                              | $\approx2.1$–$2.34$ |

$\sigma=0.8$、$\gamma=1$ 时观测到的 $0.80$ 与理论预言 $\min\{\gamma\sigma,2\}=0.8$ 精确吻合——**这不只是「格式收敛」，而是估计的尖锐性被验证了**：定理不仅给出上界，还预言了在正则性不足时会掉到多少阶。换成 $\gamma=\gamma_{\mathrm{opt}}$ 的分级网格后观测阶升到 $2$ 以上，符合 $\min\{\gamma\sigma,2\}=2$。

**算例二（液滴合并，长时间）。** $\alpha=0.4,0.7,0.9$，跑到 $T=100$，快照取 $t=1,10,50,100$，同时跟踪最大范数（始终 $\le1$，验证离散最大值原理）与离散能量（下降）。

**算例三（条件被违反时会发生什么）。** 在 $\varepsilon=0.02$ 与 $\varepsilon=0.08$ 下用不同的时间网格作对比，显示**忽略步长条件时最大值界确实会被破坏**。这类反面算例比正面算例更有价值：它说明定理里的 $\tau\le\sqrt[\alpha]{h^{2}\omega_{2-\alpha}(1-\theta)/(4\varepsilon^{2})}$ 不是证明技术的产物，而是真有内容的。

理论与实验之间仍有落差：本文证的是最大值原理与最大范数收敛，**没有**给出离散能量律的证明，能量下降只是算例二里被观察到的现象；这个缺口由编号 57 补上。

### 与其他论文的关系

它直接回答编号 40 结论部分提出的「高阶格式」公开问题，并从编号 40 继承连续能量律与最大值原理。它把 **DCC（离散互补卷积）核**及其互补恒等式与分数阶 Grönwall 机器引入本系列——这套工具随后被反复使用，并在 [[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 48、52、58、67]] 中被**对偶化**成 **DOC（离散正交卷积）核**：DCC 站在卷积的一侧、给出常数 $1$ 的互补恒等式、产出 $\ell^\infty$ 与 Grönwall 型估计；DOC 站在另一侧、给出 Kronecker $\delta$ 的正交恒等式、产出能量与 $L^2$ 估计。编号 57 是分数阶这一侧的直接续作，编号 74 则把这里隐含的「变步长 L1 型核何时正定」问题作为独立对象来解决。

## 57：把能量换掉，换来一条微分型的律

### 直觉

编号 40 的困难可以这样描述：**分数阶算子在方程的「错误一侧」。** 只要时间导数是 $\partial_t^\alpha$，用 $u_t$ 去测试就会得到一个非局部的二次型，它只在积分后才有符号。

本文的一步是把它搬到另一侧。分数阶梯度流 $\partial_t^{\alpha}u=-\delta E/\delta u$ 两端同时作用 ${}^{R}\!\partial_t^{1-\alpha}$，由半群性质左端塌回 $u_t$，于是

$$
\partial_t u=-{}^{R}\!\partial_t^{1-\alpha}\Bigl(\frac{\delta E}{\delta u}\Bigr).
$$

**现在时间导数是局部的，标准的能量测试函数可以照常使用；非局部性落在变分导数上，而那里恰好有一条可用的正定性不等式。** 这条不等式（Riemann-Liouville 型）比正定性更强一点：它把 $v\cdot{}^{R}\!\partial_t^{1-\alpha}v$ 下界化成「一个全导数」加「一个非负项」。全导数的部分可以吸收进能量的定义，非负项就是耗散。于是**新能量 $\mathcal E_\alpha$ 不是凭空造的，它就是那个全导数**。

这个构造还有一个漂亮的副产品：新律在 $\alpha\to1$ 时精确退化为经典律。也就是说 $\mathcal E_\alpha$ 不是经典能量的类比，而是它的推广。

### 问题设定

$\Omega=(0,L)^2$，周期边界，

$$
\partial_t^{\alpha}u=\varepsilon^{2}\Delta u-f(u),
\qquad F(u)=\tfrac14(1-u^2)^2,\ f=F',
$$

$$
\partial_t^{\alpha}v:=\mathcal I_t^{1-\alpha}v',
\qquad
(\mathcal I_t^{\mu}v)(t):=\int_0^t\omega_{\mu}(t-s)v(s)\,\mathrm ds,
\qquad
\omega_{\mu}(t):=\frac{t^{\mu-1}}{\Gamma(\mu)} .
$$

方程被读作**分数阶梯度流** $\partial_t^{\alpha}u=-\delta E/\delta u$，其中 $E[u]=\int_\Omega(\frac{\varepsilon^2}{2}|\nabla u|^2+F(u))\mathrm dx$。

### 推导

**第一步：改写。** 用 Riemann-Liouville 导数 ${}^{R}\!\partial_t^{\alpha}v:=\partial_t\mathcal I_t^{1-\alpha}v$ 与半群恒等式

$$
{}^{R}\!\partial_t^{1-\alpha}\bigl(\partial_t^{\alpha}v\bigr)
=\partial_t\mathcal I_t^{1}v'=v',
$$

方程等价改写为 $\partial_t u=-{}^{R}\!\partial_t^{1-\alpha}(\delta E/\delta u)$。

**第二步：变分能量与其微分律。** 用 Riemann-Liouville 不等式（归于 Alsaedi-Ahmad-Kirane）

$$
v(t)\bigl({}^{R}\!\partial_t^{1-\alpha}v\bigr)(t)
\ \ge\ \tfrac12\bigl({}^{R}\!\partial_t^{1-\alpha}v^{2}\bigr)(t)
+\tfrac12\omega_{\alpha}(t)v^{2}(t),
\qquad \forall v\in C[0,T],
$$

以及 $\frac{\mathrm dE}{\mathrm dt}=-\bigl(\frac{\delta E}{\delta u},{}^{R}\!\partial_t^{1-\alpha}\frac{\delta E}{\delta u}\bigr)$，定义

$$
\mathcal{E}_{\alpha}[u]:=E[u]+\frac12\,\mathcal{I}_t^{\alpha}
\Bigl\|\frac{\delta E}{\delta u}\Bigr\|^{2},
$$

即得

$$
\frac{\mathrm d\mathcal{E}_{\alpha}}{\mathrm dt}
+\frac12\,\omega_{\alpha}(t)\Bigl\|\frac{\delta E}{\delta u}\Bigr\|^{2}\le0,
\qquad \forall t>0 .
$$

两处 $\frac12$ 都是原文印出的。$\alpha\to1$ 时 $\mathcal I_t^\alpha\to\mathcal I_t^1$ 且 $\omega_\alpha(t)\to1$，恢复 $\frac{\mathrm dE}{\mathrm dt}+\|\frac{\delta E}{\delta u}\|^{2}\le0$。论文称这条性质为**渐近保能量耗散**。

**第三步：L1$_R$ 公式。** 变步长下把 Riemann-Liouville 导数离散为

$$
\bigl({}^{R}\!\partial_{\tau}^{1-\alpha}v\bigr)^{n-\frac12}
:=\frac{1}{\tau_{n}}\int_{t_{n-1}}^{t_{n}}\frac{\partial}{\partial t}
\int_{0}^{t}\omega_{\alpha}(t-s)(\Pi_{0}v)(s)\,\mathrm ds\,\mathrm dt
\ \triangleq\ \frac{1}{\tau_{n}}\sum_{k=1}^{n}a_{n-k}^{(n)}v^{k-\frac12},
$$

其中 $\Pi_0v$ 是在 $(t_{k-1},t_k]$ 上取值 $v^{k-\frac12}$ 的分段常数插值。辅助（DCO）序列与核为

$$
q_{n-k}^{(n)}:=\int_{t_{k-1}}^{t_{k}}\omega_{\alpha}(t_{n}-s)\,\mathrm ds
=\sum_{j=k}^{n}a_{j-k}^{(j)}>0,
$$

$$
a_{0}^{(n)}:=q_{0}^{(n)}>0\ (n\ge1),
\qquad
a_{n-k}^{(n)}:=q_{n-k}^{(n)}-q_{n-k-1}^{(n-1)}<0\ (n\ge k+1\ge2).
$$

符号模式值得注意：第一个核为正，之后全为负。显式地 $a^{(n)}_0=\omega_{1+\alpha}(\tau_n)=\tau_n^{\alpha}/\Gamma(1+\alpha)$，而 $a^{(j)}_{j-k}=\int_{t_{j-1}}^{t_j}\!\int_{t_{k-1}}^{t_k}\omega_{\alpha-1}(t-s)\,\mathrm ds\,\mathrm dt<0$。这个公式归于 Mustapha（线性次扩散情形，阶为 $1+\alpha$），作者把它命名为 L1$_R$ 以区别于 Caputo 导数的 L1 公式。

**第四步：核的正定性——在离散层面单独证明，不是从连续核继承。** 对任意实序列 $\{w_k\}$，

$$
2w_{k}\sum_{j=1}^{k}a_{k-j}^{(k)}w_{j}
\ \ge\ w_{k}^{2}\sum_{j=1}^{k}a_{k-j}^{(k)}
+\sum_{j=1}^{k}q_{k-j}^{(k)}w_{j}^{2}-\sum_{j=1}^{k-1}q_{k-j-1}^{(k-1)}w_{j}^{2},
$$

求和后

$$
2\sum_{k=1}^{n}w_{k}\sum_{j=1}^{k}a_{k-j}^{(k)}w_{j}
\ \ge\ \sum_{k=1}^{n}\Bigl(q_{n-k}^{(n)}+\sum_{j=1}^{k}a_{k-j}^{(k)}\Bigr)w_{k}^{2}
\ >\ 0,
\qquad n\ge1,\ w\not\equiv0 .
$$

证明由两条恒等式承担：$\omega_\alpha$ 的完全单调性给出

$$
q_{k-j-1}^{(k-1)}-q_{k-j}^{(k)}
=\int_{t_{j-1}}^{t_j}\bigl[\omega_{\alpha}(t_{k-1}-s)-\omega_{\alpha}(t_k-s)\bigr]\mathrm ds>0,
$$

以及

$$
\sum_{j=1}^{k}a_{k-j}^{(k)}
=\sum_{j=1}^{k}q_{k-j}^{(k)}-\sum_{j=1}^{k-1}q_{k-j-1}^{(k-1)}
=\int_{t_{k-1}}^{t_{k}}\omega_{\alpha}(s)\,\mathrm ds>0 .
$$

**这条正定性不需要任何步长比限制。** 这一点与[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|变步长 BDF]] 形成鲜明对照：那里的二次型正定性给出 $r_k<(3+\sqrt{17})/2$ 之类的门槛，而这里的分数阶核由完全单调性直接保证正定。步长上界只在最大值原理的论证中出现，与能量律无关。

**第五步：格式与两个构造性设计。** 设 $v:=-\delta E/\delta u$，把方程分裂为 $\partial_t u={}^{R}\!\partial_t^{1-\alpha}v$、$v=\varepsilon^{2}\Delta u-f(u)$，离散为 Crank-Nicolson 型

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

（$\circ$ 与幂为逐元素运算）。这个特定的 $H$ 是**设计出来的**，使 $H(a,b)(a-b)\ge F(a)-F(b)$ 逐点成立——这是链式法则的离散对应物，也是能量论证得以进行的原因。离散变分能量为

$$
\mathcal{E}_{\alpha}[u^{n}]:=E[u^{n}]+\frac12h^{2}\sum_{i,j=1}^{M_1}\sum_{k=1}^{n}
q_{n-k}^{(n)}\bigl(v_{ij}^{k-\frac12}\bigr)^{2},
\qquad
E[u^{n}]:=h^{2}\sum_{i,j}F(u_{ij}^{n})-\tfrac12\varepsilon^{2}h^{2}(u^{n})^{T}D_{h}u^{n},
$$

其中 $q$ 核本身构成一个数值分数阶积分：$(\mathcal I_\tau^{\alpha}v)^n=\sum_{k}q_{n-k}^{(n)}v^{k-\frac12}$，并且 $({}^{R}\!\partial_{\tau}^{1-\alpha}v)^{n-\frac12}=\partial_\tau(\mathcal I_\tau^{\alpha}v)^{n-\frac12}$。

**第六步：DOC 核与两种分数阶导数之间的可逆变换。** 定义 L1$_R$ 核的离散正交卷积核 $\theta$：

$$
\theta_{0}^{(n)}:=\frac{1}{a^{(n)}_{0}},
\qquad
\theta_{n-k}^{(n)}:=-\frac{1}{a^{(k)}_{0}}\sum_{j=k+1}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}
\ \ (1\le k\le n-1).
$$

**双向**正交性同时成立（这就是可逆性），并与 DCO 核 $q$ 互补：

$$
\sum_{j=k}^{n}a^{(n)}_{n-j}\theta_{j-k}^{(j)}\equiv\delta_{nk},
\qquad
\sum_{j=k}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}\equiv\delta_{nk},
\qquad
\sum_{j=k}^{n}q_{n-j}^{(n)}\theta_{j-k}^{(j)}\equiv1 .
$$

论文证明 $\theta^{(n)}_0=\Gamma(1+\alpha)\tau_n^{-\alpha}$、全部 $\theta^{(n)}_j>0$、尖锐的首项间隙

$$
\theta_{0}^{(n)}-\theta_{1}^{(n)}
=\frac{\omega_{1+\alpha}(r_{n}+1)-\omega_{1+\alpha}(r_{n})}
{\omega_{1+\alpha}(\tau_{n})\,\omega_{1+\alpha}(1)}
\ >\ \frac{\omega_{\alpha}(r_{n}+1)}{\omega_{1+\alpha}(\tau_{n})\,\omega_{1+\alpha}(1)},
$$

以及对 $n\ge2$ 的**单调递减** $\theta_{0}^{(n)}>\theta_{1}^{(n)}>\cdots>\theta_{n-1}^{(n)}>0$。单调性的证明引入辅助核 $\zeta^{(n)}_0:=\theta^{(n)}_0$、$\zeta^{(n)}_{n-j}:=\theta^{(n)}_{n-j}-\theta^{(n)}_{n-j-1}$（满足 $\sum_{j=k}^{n}\zeta^{(n)}_{n-j}q^{(j)}_{j-k}=\delta_{nk}$），并**直接调用[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 74]] 的代数判据** $q^{(n)}_j>0$、$q^{(n-1)}_{j-1}>q^{(n)}_j$、$q^{(n-1)}_{j-1}q^{(n)}_{j+1}>q^{(n-1)}_{j}q^{(n)}_{j}$——因此编号 74 在逻辑上先于本文，尽管出版更晚。

把格式的第一个方程与 $\theta^{(n)}_{n-j}$ 卷积并用正交性，

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

也就是说 DOC 核定义了一个**新的**离散 Caputo 导数，其核在非均匀网格上正且单调递减，与经典 L1 核同性质，而变换由双向正交性可逆。论文自己提醒这是一个**间接**逼近，其精度与直接 L1 公式（误差阶 $2-\alpha$）不同。

最大值原理的论证正是在这个 Caputo 等价形式上、按归纳法进行的：用 $D_h$ 的对称负半定性、一条 $\ell^\infty$ 引理

$$
\bigl\|(aI-B)V+U^{.2}\!\circ V+cV^{.3}\bigr\|_{\infty}
\ \ge\ a\|V\|_{\infty}+\|U\|_{\infty}^{2}\|V\|_{\infty}+c\|V\|_{\infty}^{3},
$$

以及整理后的格式

$$
\Bigl(\theta^{(n)}_{0}-\tfrac12+\tfrac12(u^{n-1})^{.2}-\tfrac{\varepsilon^{2}}{2}D_{h}\Bigr)u^{n}
+\tfrac13(u^n)^{.3}=\cdots+\mathcal{L}^{n-2}(u),
\qquad
\mathcal{L}^{n-2}(u):=\sum_{k=1}^{n-2}\bigl(\theta_{n-k-1}^{(n)}-\theta_{n-k}^{(n)}\bigr)u^{k}
+\theta_{n-1}^{(n)}u^{0},
$$

其历史系数非负，**恰恰因为 DOC 核正且递减**。

### 定理

**（离散变分能量耗散律，无条件）** Crank-Nicolson 格式**无条件**（无步长大小限制、无步长比限制）满足

$$
\partial_{\tau}\bigl(\mathcal{E}_{\alpha}[u]\bigr)^{n-\frac12}
+\frac{1}{2\tau_{n}}\int_{t_{n-1}}^{t_{n}}\omega_{\alpha}(s)\,\mathrm ds
\sum_{i,j=1}^{M_1}h^{2}\bigl(v_{ij}^{n-\frac12}\bigr)^{2}\ \le\ 0,
\qquad n\ge1 .
$$

两个成分：$H$ 的链式法则性质，以及 L1$_R$ 正定性引理的逐层形式 $v^{n-\frac12}\sum_k a^{(n)}_{n-k}v^{k-\frac12}\ge\frac12(\mathcal I_\tau^{\alpha}v^2)^n-\frac12(\mathcal I_\tau^{\alpha}v^2)^{n-1}+\frac12\int_{t_{n-1}}^{t_n}\omega_\alpha(s)\mathrm ds\,(v^{n-\frac12})^2$。

**（$\alpha\to1$ 的渐近保持）** 当 $\alpha\to1$ 时 $q^{(n)}_{n-k}\to\tau_k$，于是 $(\mathcal I_\tau^\alpha v)^n\to\sum_k\tau_kv^{k-\frac12}$、$({}^{R}\!\partial_\tau^{1-\alpha}v)^{n-\frac12}\to v^{n-\frac12}$，离散律变为

$$
\partial_{\tau}\bigl(E[u]\bigr)^{n-\frac12}
+\sum_{i,j}h^{2}\bigl(v_{ij}^{n-\frac12}\bigr)^{2}\ \le\ 0,
$$

即经典 Allen-Cahn 的标准离散能量耗散律。

**（唯一可解性）** 若 $\tau<\sqrt[\alpha]{2\Gamma(1+\alpha)}$，则 $\theta^{(n)}_0=\Gamma(1+\alpha)\tau_n^{-\alpha}>\frac12$，$G_h:=(\theta^{(n)}_0-\frac12+\frac12(u^{n-1})^{.2})I-\frac{\varepsilon^2}{2}D_h$ 正定，目标泛函 $\frac12w^TG_hw+\frac{1}{12}\sum_kw_k^4-w^T\mathrm G_0(u^{n-1})$ 严格凸，格式唯一可解。

**（离散最大值原理）** 若步长满足

$$
\tau_{n}\le\sqrt[\alpha]{\min\Bigl\{\frac12,\ \frac{h^{2}}{2\varepsilon^{2}}\Bigr\}
\cdot\frac{\alpha\,\Gamma(1+\alpha)}{(1+r_{n})^{1-\alpha}}},
$$

则 $\|u^{0}\|_{\infty}\le1$ 蕴含 $\|u^{k}\|_{\infty}\le1$。**注意 $r_n=\tau_n/\tau_{n-1}$ 只通过温和的因子 $(1+r_n)^{\alpha-1}<1$ 出现，对 $r_n$ 没有任何上界要求。** 这正是本文能够宣称「第一个在变步长下同时保住能量稳定性与最大值原理的工作」的原因。

### 数值实验

实现同样用基于指数和逼近的快速 L1$_R$ 算法（绝对容差 $\epsilon=10^{-12}$，截断时间 $\Delta t=10^{-12}$）。

| 算例 | 设置                                                                                                                                                                                                       | 检验对象                                         |
| ---- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------ |
| 1    | 带外力模型 $\partial_tu=-{}^{R}\!\partial_t^{1-\alpha}(\delta E/\delta u)+g$，$(0,1)^2\times(0,1]$，$\varepsilon=0.1$，制造解 $u=\omega_{1+\sigma}(t)\sin(2\pi x)\sin(2\pi y)$，如 $\alpha=0.6,\sigma=0.4$ | 时间收敛阶（期望 $1+\alpha$）                    |
| 2    | 最大值原理：$\alpha=0.7,\,0.9$，$\tau=0.1,\,0.8,\,1.0$                                                                                                                                                     | 离散最大范数是否 $\le1$                          |
| 3    | 粗化：$(0,2\pi)^2$，$\varepsilon=0.05$，$128\times128$ 空间网格，随机初值均匀分布于 $[-0.001,0.001]$                                                                                                       | 初始层 $u_t=\mathcal O(t^{\alpha-1})$ 与分级网格 |
| 4    | 自适应步长                                                                                                                                                                                                 | $E(t)$ 与 $\mathcal E_\alpha(t)$ 的对照          |

算例一的期望阶是 L1$_R$ 公式本身的阶 $1+\alpha$；在存在 $t^{\sigma}$ 初始奇性时，分级网格把它恢复回来。算例二用的 $\tau=0.8$ 与 $\tau=1.0$ 是**很大**的步长，最大范数依然不超过 $1$，与定理里步长条件的宽松程度相符。

**算例四是全篇最值得看的一张图，因为它把编号 40 的警告变成了可见的现象**：在自适应步长下同时画出原能量 $E(t)$ 与变分能量 $\mathcal E_\alpha(t)$，$\mathcal E_\alpha$ 单调下降，而 $E$ **不必**单调。这正是编号 40 所说「积分型律不蕴含逐点下降」的数值体现，也说明本文换掉能量对象不是修辞上的方便，而是必需的。

理论与实验的落差在于：本文的定理覆盖能量稳定性、可解性与最大值原理，但**没有给出误差估计**（那是编号 43 在 Alikhanov 公式下做的事），算例一测到的 $1+\alpha$ 阶因此是经验的。此外论文自己提醒，DOC 核诱导的那个「新 Caputo 导数」是间接逼近，其精度与直接 L1 公式不同。

### 与其他论文的关系

这是编号 40 与 43 在分数阶一侧的直接续作，也是[[computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|编号 48]] 结论中「为时间分数阶相场方程发展非均匀 BDF2 型格式」这一公开问题的答案。它把编号 40 的积分型律 $E[u](T)\le E[u](0)$ 换成一条对重新定义的变分能量 $\mathcal E_\alpha$ 成立的**微分型**律，并在 $\alpha\to1$ 时退化回经典律。它使用编号 52 与 Liao-Zhang 为 BDF2 引入的 DOC 核，但用途是新的：不是用来剥掉多步算子，而是用来**在两种离散分数阶导数之间作变换**。它的 DOC 单调性证明显式调用编号 74 的代数判据。与编号 43 合起来看，两者构成一对：「只保最大值原理、二阶、Alikhanov」（43）对「能量加最大值原理、步长比无限制、L1$_R$/Crank-Nicolson」（57）。

## 三篇的推进关系

| 编号 | 能量律的形式             | 步长限制                               | 主要分析装置                          | 数值实验的主产出                             |
| ---- | ------------------------ | -------------------------------------- | ------------------------------------- | -------------------------------------------- |
| 40   | 积分型（$[0,T]$ 上累积） | 连续层面不涉及；离散有稳定化条件       | 分数阶卷积核的正定性                  | $-\alpha/3$ 粗化幂律（经验）                 |
| 43   | 未证；只证最大值原理     | M1：$\rho=7/4$，加 $\tau$-$h$ 耦合     | Alikhanov 核估计 + DCC 核 + Grönwall  | 分级网格上 $\min\{\gamma\sigma,2\}$ 阶被验证 |
| 57   | 微分型（变分能量）       | 能量律无限制；最大值原理有 $\tau$ 上界 | 方程改写 + L1$_R$ 核正定性 + DOC 变换 | $\mathcal E_\alpha$ 单调而 $E$ 不必单调      |

编号 40 到编号 57 的推进值得单独总结：**当一条定律的形式不对时，不要削弱结论，而要换一个对象。** 编号 40 得到的是原能量的积分型不等式；编号 57 不去加强它，而是构造一个新的能量 $\mathcal E_\alpha$，使其满足微分型律并在 $\alpha\to1$ 时退化回经典律。代价是这个能量含一个分数阶积分项，因此不是原能量本身——而编号 57 的算例四恰恰显示，原能量确实可以不单调，所以这个代价是不可避免的。

## 原文

- T. Tang, H. Yu, and T. Zhou, [_On energy dissipation theory and numerical stability for time-fractional phase-field equations_](https://doi.org/10.1137/18M1203560), SIAM J. Sci. Comput. 41(6) (2019), pp. A3757-A3778（预印本 [arXiv:1808.01471](https://arxiv.org/abs/1808.01471)）。
- H.-l. Liao, T. Tang, and T. Zhou, [_A second-order and nonuniform time-stepping maximum-principle preserving scheme for time-fractional Allen-Cahn equations_](https://doi.org/10.1016/j.jcp.2020.109473), J. Comput. Phys. 414 (2020), 109473（预印本 [arXiv:1909.10216](https://arxiv.org/abs/1909.10216)）。
- H.-l. Liao, T. Tang, and T. Zhou, [_An energy stable and maximum bound preserving scheme with variable time steps for time fractional Allen-Cahn equation_](https://doi.org/10.1137/20M1384105), SIAM J. Sci. Comput. 43(5) (2021), pp. A3503-A3526（预印本 [arXiv:2012.10740](https://arxiv.org/abs/2012.10740)）。
