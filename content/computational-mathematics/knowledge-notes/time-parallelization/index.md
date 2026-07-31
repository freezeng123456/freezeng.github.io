---
title: 双曲与抛物问题的时间并行方法
description: 时间并行方法的分章笔记与可复现实验
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization
tags:
  - 计算数学
  - 时间并行
---

本专题依据 M. J. Gander、S.-L. Wu 和 T. Zhou 的综述 _Time Parallelization for Hyperbolic and Parabolic Problems_（Acta Numerica 34, 2025, pp. 385-489）整理。原文区分了对传播型问题仍然有效的方法，以及主要为耗散问题设计的方法。

## 阅读顺序

1. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-1-why-parallelize-in-time|第一章：为什么要做时间并行？]]介绍因果性、全时间耦合形式与性能判断标准。
2. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-2-model-problems|第二章：模型问题]]比较热方程、对流扩散方程、Burgers 方程和波动方程，并收录三个重新计算的解实验。
3. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-hyperbolic-methods|第三章：适用于双曲问题的方法]]讨论 SWR、PIDC/RIDC、ParaExp 和 ParaDiag，并收录 Heat、ADE 与 Wave 三类 ParaDiag-II 实验。
4. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-4-parabolic-methods|第四章：为抛物问题设计的方法]]介绍 Parareal、PFASST、MGRiT、基于对角化的 Parareal 和 STMG，并给出重新计算的收敛实验。
5. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view|第五章：统一视角与方法选择]]比较各类算法，记录完整实验清单和复现流程。

## 方法图谱

| 方法      | 并行单元           | 机制             | 自然适用区域     |
| --------- | ------------------ | ---------------- | ---------------- |
| SWR       | 重叠时空子域       | 波形传输         | 输运与波动问题   |
| PIDC/RIDC | 校正层与时间节点   | 延迟校正与流水线 | 初值问题         |
| ParaExp   | 非齐次与齐次子问题 | 指数传播         | 线性系统         |
| ParaDiag  | 全时间耦合矩阵     | 循环近似与 FFT   | 线性或线性化系统 |
| Parareal  | 粗时间区间         | 粗预测加细校正   | 中等到强耗散     |
| PFASST    | 跨时间步的配置节点 | SDC 与多层校正   | 高阶时间积分     |
| MGRiT     | 时间网格层次       | 松弛与粗网格校正 | 长时间区间       |
| STMG      | 完整时空网格       | 时空平滑与粗化   | 抛物系统         |

## 三条组织原则

### 因果性被重新表述，而不是被消除

单步方法 $u_{n+1}=\Phi_{\Delta t}(u_n)$ 本质上是串行的。时间并行方法转而构造全时间耦合系统逆算子的并行近似，使并行工作发生在迭代、分解或变换内部。

### 耗散决定粗表示是否包含有效信息

抛物动力学会衰减高频误差，因此粗传播子仍能表示长时间尺度上保留下来的慢变分量。双曲动力学保留相位信息，粗问题中的微小相位误差可能沿多个时间区间累积。

### 迭代次数不等于并行效率

一个简化成本模型是

$$
T_{\mathrm{parallel}}
\approx K(C_G+C_F/P)+C_{\mathrm{comm}}+C_{\mathrm{setup}},
$$

其中 $K$ 是迭代次数，$C_G$ 和 $C_F$ 分别是粗、细传播成本，$P$ 是时间并发度。本专题的数值实验测量算法收敛性，不测量端到端并行加速比。

> [!note] 数值来源
> 所有展示的实验图均于 2026 年 7 月 31 日在新的实验服务器上，由 Python 复现项目重新生成。初始正式结果使用 SciPy CPU 路径；随后加入 CuPy/T4 混合后端，将独立 Burgers 细传播批量放到 GPU。Figure 4.5 的停止迭代保持一致，完整论文套件由 263.57 秒降至 67.92 秒。

## 第一手资料

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), Acta Numerica 34 (2025), pp. 385-489.
- 原始 MATLAB 算例：[wushulin/ActaPinT](https://github.com/wushulin/ActaPinT)。
- Python 转换、扩充与正式结果：[freezeng123456/ActaPinT-Python](https://github.com/freezeng123456/ActaPinT-Python)。
