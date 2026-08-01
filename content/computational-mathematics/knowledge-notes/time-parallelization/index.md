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

1. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-1-why-parallelize-in-time|第一章：为什么要做时间并行？]]逐段梳理论文摘要与引言，介绍硬件背景、因果链、四条方法谱系、双重分类和全时间系统。
2. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-2-model-problems|第二章：模型问题]]逐图比较热方程、对流扩散方程、Burgers 方程和波动方程的边界条件、时间记忆，并收录三个重新计算的解实验。
3. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-3-hyperbolic-methods|第三章：对双曲问题仍然有效的方法]]完整推导 SWR、PIDC/RIDC、ParaExp 和 ParaDiag，并收录 Heat、ADE 与 Wave 三类 ParaDiag-II 实验。
4. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-4-parabolic-methods|第四章：主要为抛物问题设计的方法]]分析 Parareal、PFASST、MGRiT、两种对角化 Parareal 和 STMG，并给出重新计算的收敛实验。
5. [[computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view|第五章：结论、统一视角与复现边界]]区分论文结论与本站综合，记录方法选择、完整实验清单、GPU 性能和报告规范。

## 论文内容覆盖进度

“段落级完成”表示正文论点、公式、图示、历史线索、限定条件和章节关系都已与原文逐项核对。已有实验不会因此删除，它们会放在对应原文解释之后。

| 原文范围           | 网页对应 | 当前状态                                                                                             |
| ------------------ | -------- | ---------------------------------------------------------------------------------------------------- |
| 摘要与 Section 1   | 第一章   | **段落级完成**：覆盖 pp. 385–388 的全部论点，并重绘 Figure 1.1                                       |
| Sections 2.1–2.4   | 第二章   | **段落级完成**：覆盖四类模型、全部边界设置、Figures 2.1–2.4 的各组观察，并保留三个补充实验           |
| Sections 3.1–3.5.2 | 第三章   | **段落级完成**：覆盖历史、SWR、PIDC/RIDC、ParaExp、ParaDiag-I/II 的推导、定理与全部论文数值讨论      |
| Sections 4.1–4.6   | 第四章   | **段落级完成**：覆盖 Parareal、PFASST、MGRiT、对角化变体、STMG、Theorems 4.1–4.9 与 Figures 4.1–4.22 |
| Section 5          | 第五章   | **段落级完成**：覆盖论文结论，并将统一视角、GPU 分析和实验清单明确标为本站补充                       |

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

### 因果约束如何进入并行计算

单步方法 $u_{n+1}=\Phi_{\Delta t}(u_n)$ 形成串行递推。时间并行方法把所有时间点的耦合显式写出，再构造其逆算子的并行近似。并发工作由此进入迭代、分解或变换过程。

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
