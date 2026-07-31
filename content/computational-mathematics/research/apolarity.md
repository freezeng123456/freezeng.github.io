---
title: Apolarity 引导的 Taylor Jet
description: 面向指定混合偏导数的 Waring 最优方向计划
lang: zh
translation: en/computational-mathematics/research/apolarity
tags:
  - 计算数学
  - 自动微分
  - 代数几何
---

代码仓库：[freezeng123456/apolarity](https://github.com/freezeng123456/apolarity)

## 问题

为了计算一个指定的 $p$ 阶混合偏导数，嵌套自动微分会构造越来越深的导数计算图。Taylor 模态自动微分可以沿一个方向传播截断 jet，但这又引出一个组合问题：哪些方向能够精确恢复目标混合偏导数？所需方向的最小数量是多少？

## 核心 insight

方向 Taylor 系数

$$
T_p(x;v)=\frac1{p!}\frac{d^p}{d\tau^p}u(x+\tau v)\bigg|_{\tau=0}
$$

关于方向 $v$ 是一个 $p$ 次齐次多项式。把指定混合偏导数写成若干 $T_p(x;v_r)$ 的线性组合，等价于把关联单项式 $z^\alpha$ 写成若干线性形式 $p$ 次幂之和，即 Waring 分解。

若非零指数按 $a_0\le\cdots\le a_n$ 排列，复数域 Waring 秩给出最少方向数：

$$
R_{\mathbb C}(z^\alpha)=\prod_{j=1}^{n}(a_j+1).
$$

单位根构造能够达到该秩，并产生可执行的复方向计划。

## 计算结构

```text
多重指标
  -> 非零指数
  -> 单位根 Waring 计划
  -> 方向 Taylor jet
  -> 加权得到精确导数
  -> 自定义反向规则
  -> 参数梯度训练
```

前向传播在每一层线性映射和 $\sinh$ 递推中传播至 $p$ 阶的 Taylor 系数；自定义反向规则再对导数计算结果关于网络参数求导。

## 代表性结果

- 在论文报告的验证实验中，计算得到的导数值及参数梯度与嵌套自动微分一致。
- 对含重复指标的混合偏导数，复 Waring 计划可少于实数极化所需的方向数。对 square-free 导数，两种构造都需要 $2^{p-1}$ 个方向，因此复数构造在这种情况下没有方向数量优势。
- 在统一网络宽度、统一 1200 秒预算和 5 个随机种子下，复数-$\sinh$ 模型在报告的 12 个 PDE 配置中取得最低的平均相对 $L^2$ 误差。测试覆盖多调和、径向 chirp 和有耗 Maxwell 问题。

![[assets/research/apolarity-chirp.png]]

## 适用边界

该构造是一个针对指定单项式偏导数的精确后端，并非对 Laplace 算子或 $\Delta^m$ 等算子和的通用加速。更短的代数计划也不意味着墙钟时间按比例下降，因为复数算术、批处理、数值精度和内存布局都会影响实测成本。

PDE 结果为论文设置下的复数-$\sinh$ 表示提供了证据，但不能把全部精度收益仅归因于 Waring 计划。若要做出这一归因，仍需匹配参数量、运行时间和峰值内存的对照实验。
