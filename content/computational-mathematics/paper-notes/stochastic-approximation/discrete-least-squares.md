---
title: 离散最小二乘逼近
description: 编号 6、9、11、13、14：从「取哪些点」到「取多少个点才稳定」
lang: zh
translation: en/computational-mathematics/paper-notes/stochastic-approximation/discrete-least-squares
tags:
  - 论文笔记
  - 不确定性量化
  - 最小二乘
---

> [!note] 本页覆盖
> 编号 **6**（_Commun. Comput. Phys._ 16, 2014）、**9**（_SIAM J. Sci. Comput._ 36(5), 2014）、**11**（_SIAM J. Sci. Comput._ 36(5), 2014）、**13**（_J. Comput. Phys._ 298, 2015）、**14**（_Commun. Comput. Phys._ 18, 2015）。

![配点设计的统一流程](assets/diagrams/tao-zhou-papers/zh/sampling-design.svg)

## 6：先用实验回答「取哪些点」

编号 6 是整个采样设计纲领的经验序幕。它直接问最小二乘多项式逼近的设计点应当怎么选，并用实验回答。它想消除的具体担忧是「随机点会破坏收敛率」——而这一担忧后来被编号 9、22、45 的稳定性理论解释掉了：真正的问题不是收敛率，而是**要让设计矩阵条件良好需要多少个样本**。

## 9：确定性点集与二次样本复杂度

### 稳定性定理

编号 9 是采样设计一族中第一篇给出理论的。它构造一类新的配点网格（基于 Weil 和的点集），并证明**确定性**的稳定性结论。设 $N=\#\Lambda$ 为多项式空间维数，$A$ 为相应的设计（Gram）矩阵，则若

$$
M\ \ge\ 4^{\,d+1}\cdot d^2\cdot N^2
$$

且 $M$ 为素数，则在谱范数下

$$
\Bigl|\!\Bigl|\!\Bigl|\frac{2^{d+1}}{M}A-I\Bigr|\!\Bigr|\!\Bigr|\le\frac12 .
$$

证明路径很短：取 $\delta=2^{d}((d-1)\sqrt{M}+1)/M$，用 Gershgorin 定理得 $|\lambda_i-1|\le N\delta$，再验证上述条件下 $N\delta\le\frac12$。

**这就是二次样本复杂度的陈述：$M\gtrsim N^2$，并带一个显式的、依赖维数的前因子 $4^{d+1}d^2$。** 它是后面几篇要攻击的基准。

由此立即得到唯一性（$A$ 严格对角占优故非奇异）与最优收敛：设 $P^\Lambda f$ 为 Chebyshev 测度下的最佳 $L^2$ 逼近、$P^\Lambda_m f$ 为该点集上的离散最小二乘解，则

$$
\|f-P^\Lambda_m f\|_{L^2_{\rho_c}}\ \le\ \Bigl(1+\frac{4}{d^2 N}\Bigr)\,\|f-P^\Lambda f\|_{L^\infty} .
$$

因子 $1+4/(d^2N)$ 随 $N$ 增大趋于 $1$，因此离散投影渐近上与最佳 $L^\infty$ 逼近同样好。**这个估计是确定性的，不带「以高概率」。** 这是它相对随机采样的代价所换来的收益。

### 覆盖到其他测度

论文引入一个支配条件：若存在与 $Y$ 无关的常数 $C$ 使 $0<\rho(Y)\le C\rho_c(Y)$（$\rho_c$ 为 Chebyshev 密度），则

$$
\|f-P^\Lambda_m f\|_{L^2_\rho}\ \le\ \sqrt{C}\,\Bigl(1+\frac{4}{d^2 N}\Bigr)\|f-P^\Lambda f\|_{L^\infty} .
$$

论文指出这覆盖均匀测度以及所有满足 $0<\rho_{\min}\le\rho\le\rho_{\max}$ 的测度，并进一步指出它对**认知不确定性**有用：当 $Y$ 的密度未知时，只要该未知密度满足支配条件，基于 Chebyshev 的逼近仍然有效。

**这篇论文确立了这一族的模板**：一个显式的点集或密度，一条量化所需样本数的稳定性定理，以及一条近最佳逼近的推论。

## 11：无界区域上情况在质上更差

### 问题

已有的最小二乘稳定性理论都针对**有界**参数域——$[-1,1]^d$ 上的均匀或 Chebyshev 测度，那里 $M\sim N^2$ 个样本够用（Chebyshev 更少）。而高斯与 Gamma 参数生活在无界域上，论文指出那里的情况在**质**上更差而不只是量上更差：用自然的 Hermite 或 Laguerre 多项式混沌基、并从高斯或 Gamma 测度采样时，设计矩阵只有在样本数关于逼近空间维数**指数增长**时才条件良好。

第二个独立的问题是 Hermite 展开的分辨率很差：论文直接引用 Gottlieb 与 Orszag 的说法——分辨 $\sin(x)$ 的 $M$ 个波长需要接近 $M^2$ 个 Hermite 多项式。

### 对策与它在这一族中的位置

论文的解决路线是同时更换**基**（从多项式换成函数）与**采样密度**（从高斯换成映射后的均匀分布），由此得到本专题第一个**对数线性**而非二次的样本要求。

**这条路线与后来的 Christoffel 加权是两个不同的选择**：编号 11 改基与密度；编号 22 与 45 保留多项式基，改密度并加权。两者都达到对数线性复杂度，但代价不同——前者放弃了多项式基的正交结构，后者需要计算 Christoffel 函数。

## 13 与 14：随机求积与非结构网格

### 13：Gauss 权就是 Christoffel 函数值

编号 13 对张量 Gauss 网格做随机子采样。它的关键观察是 Gauss 权恰好是 Christoffel 函数值：

$$
w_{\mathbf k}=\lambda_{\mathbf n}(z_{\mathbf k})
=\prod_{i=1}^{d}\lambda^i_{n_i}(z^i_{k_i})
=\prod_{i=1}^{d}\frac{1}{\sum_{k=0}^{n_i-1}\bigl[\phi^i_k(z^i_{k_i})\bigr]^2} .
$$

**因此「对 Gauss 网格做子采样」隐含地就是 Christoffel 加权采样。** 这一观察是编号 22 与 45 把采样密度直接取成 Christoffel 函数倒数的直接前身。

按第三方文献的表述，该方法在样本数关于多项式维数**线性**增长时即稳定——强于诱导采样的 $M\gtrsim N\log N$ 与朴素 Monte Carlo 的 $M\gtrsim N^2$。但这是在「张量 Gauss 网格子集」这个受限设计空间内的结论，其精确假设与常数本站未核实。

### 14：非结构多元网格

编号 14 处理参数域不是张量积结构的情形，在非结构多元网格上做随机配点。这一篇是综述性质较强的一篇（篇幅 36 页），在这一族中提供的是把前面几篇的构造放在一起比较的框架。

## 五篇的样本复杂度对照

| 编号 | 点集或密度                 | 样本要求                    | 结论类型         |
| ---- | -------------------------- | --------------------------- | ---------------- |
| 6    | 若干候选设计点（实验比较） | 不适用                      | 经验             |
| 9    | Weil 和确定性点集          | $M\ge4^{d+1}d^2N^2$（素数） | 确定性           |
| 11   | 映射后的均匀密度 + 函数基  | 对数线性                    | 概率（限定核实） |
| 13   | 张量 Gauss 网格的随机子集  | 线性（受限设计空间）        | 概率（限定核实） |
| 14   | 非结构多元网格             | 见原文                      | 综述与构造       |

一条贯穿的判断：**「取哪些点」这个问题的正确形式是「取多少个点，从什么密度取」。** 编号 6 问的是前者，编号 9 起把它换成后者，而后者有定理可证。这一转换是这条线索得以推进的原因。

> [!note] 覆盖进度
> 编号 9 与 11 的定理与样本复杂度陈述已核对。编号 6、13、14 的具体实验、假设与常数本站未逐项核实；上述内容限于可从摘要、第三方文献与参考关系确认的部分。

## 本页原文

- Z. Gao and T. Zhou, [_On the choice of design points for least square polynomial approximations with application to uncertainty quantification_](https://doi.org/10.4208/cicp.130813.060214a), Commun. Comput. Phys. 16 (2014), pp. 365-381。
- T. Zhou, A. Narayan, and Z. Xu, [_Multivariate discrete least-squares approximations with a new type of collocation grid_](https://doi.org/10.1137/130950434), SIAM J. Sci. Comput. 36(5) (2014), pp. A2401-A2422。
- T. Tang and T. Zhou, [_On discrete least-squares projection in unbounded domain with random evaluations and its application to parametric uncertainty quantification_](https://doi.org/10.1137/140961894), SIAM J. Sci. Comput. 36(5) (2014), pp. A2272-A2295。
- T. Zhou, A. Narayan, and D. Xiu, [_Weighted discrete least-squares polynomial approximation using randomized quadratures_](https://doi.org/10.1016/j.jcp.2015.06.042), J. Comput. Phys. 298 (2015), pp. 787-800。
- A. Narayan and T. Zhou, [_Stochastic collocation on unstructured multivariate meshes_](https://doi.org/10.4208/cicp.020215.070515a), Commun. Comput. Phys. 18 (2015), pp. 1-36。
