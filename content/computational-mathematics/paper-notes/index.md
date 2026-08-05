---
title: 周涛研究工作精读
description: 按七个专题整理周涛主页发表列表的双语论文笔记
lang: zh
translation: en/computational-mathematics/paper-notes
tags:
  - 计算数学
  - 论文笔记
  - 不确定性量化
---

本专题系统整理周涛（中国科学院数学与系统科学研究院计算数学与科学工程计算研究所）在个人主页上列出的全部 108 项工作，其中 95 项已发表、13 项处于投稿或预印本状态。目录、题录与主题划分依据 2026 年 8 月 5 日抓取的[主页发表列表](https://lsec.cc.ac.cn/~tzhou/)，DOI 通过 Crossref 逐条核对。

> [!info] 本专题的定位与边界
> 这些页面是**第三方阅读笔记**，不是原文复述，也不代表作者本人的表述。每页给出问题背景、方法构造、主要结论与限定条件，并在结尾列出原文出处；证明细节请回到原文。已发表工作的技术内容依据公开摘要、公开预印本与公开正文整理；对无法从公开材料核实的细节，页面会明确标注为待核实，而不做推测性补全。

![七个专题共享同一条技术主线](assets/diagrams/tao-zhou-papers/zh/research-map.svg)

## 为什么按专题读，而不按年份读

主页的发表列表按时间倒序排列，这种排列适合检索，不适合理解。按时间读会反复在互不相邻的问题之间跳转：2014 年的离散最小二乘、2015 年的 parareal 收敛分析、2015 年的正倒向随机微分方程多步格式，各自属于不同的技术传统。

按专题读能看到另一件事：同一个技术动作在不同问题上反复出现。把串行的时间递推写成一个全时间算子，再为它构造并行的近似逆；把随机输入的高维积分写成加权最小二乘，再为它设计采样密度；把方程残差写成沿轨道的鞅性质，再用两个网络对抗训练。这三件事表面上分属时间并行、不确定性量化和科学机器学习，实际共享同一条思路：**先把难以直接处理的结构显式写成一个算子方程，再设计带权、带自适应的近似求逆过程，并用一个误差指标决定在哪里补算真实信息。**

## 七个专题

| 专题                                                                                                | 论文数 | 核心问题                                               | 主要技术量                     |
| --------------------------------------------------------------------------------------------------- | ------ | ------------------------------------------------------ | ------------------------------ |
| [[computational-mathematics/paper-notes/stochastic-approximation/index\|随机逼近与配点设计]]        | 23     | 用尽可能少的正问题求解，稳定地构造参数化解的多项式逼近 | Christoffel 函数、诱导采样密度 |
| [[computational-mathematics/paper-notes/spectral-and-reduced-order/index\|谱方法与降阶表示]]        | 6      | 无界区域与低秩结构下如何保持谱精度                     | 基函数远场衰减率、最小奇异值   |
| [[computational-mathematics/paper-notes/bayesian-inference/index\|贝叶斯反问题与数据同化]]          | 10     | 在昂贵正问题下完成后验采样而不引入代理偏差             | 后验误差指标、接受率修正       |
| [[computational-mathematics/paper-notes/scientific-machine-learning/index\|科学机器学习]]           | 22     | 神经网络求解 PDE 时，采样点与频谱结构决定成败          | 失效概率、频率自适应尺度       |
| [[computational-mathematics/paper-notes/fbsde-and-control/index\|FBSDE 与随机最优控制]]             | 22     | 用概率表示绕开高维网格，同时保持高阶时间精度           | 条件期望、鞅性质               |
| [[computational-mathematics/paper-notes/phase-field-and-time-stepping/index\|相场模型与变步长离散]] | 12     | 自适应步长下如何保住离散能量律与最大值原理             | 卷积核重排、步长比门槛         |
| [[computational-mathematics/paper-notes/parallel-in-time/index\|时间并行算法]]                      | 13     | 把串行时间推进改写成可并发的代数问题                   | 对角化条件数、全时间预条件     |

完整题录见[[computational-mathematics/paper-notes/catalog|论文总目录]]，其中按专题列出全部 108 项工作的编号、题名、期刊与对应精读页。

## 时间分布

![2010 至 2026 年的专题分布](assets/diagrams/tao-zhou-papers/zh/research-timeline.svg)

时间轴显示了三次重心转移。2010 至 2015 年集中在随机 Galerkin、随机配点与离散最小二乘，问题是「如何用有限的正问题求解刻画参数依赖」。2016 至 2021 年出现两条并行主线：一条把采样设计推向最优采样与稀疏恢复，并进入贝叶斯反问题；另一条转向时间方向，包括变步长相场格式的能量分析与基于对角化的时间并行算法。2022 年之后科学机器学习成为主要方向，但技术内核仍来自前两个阶段：自适应采样继承了配点设计，密度流继承了变量替换，鞅方法继承了正倒向随机微分方程的概率表示。

## 与站内其他专题的关系

[[computational-mathematics/knowledge-notes/time-parallelization/index|双曲与抛物问题的时间并行方法]]是对第 85 项工作（Gander、Wu 与 Zhou 的 2025 年 _Acta Numerica_ 综述）的逐节精读，含原论文图表与可复现实验。本专题的[[computational-mathematics/paper-notes/parallel-in-time/index|时间并行算法]]页面不重复那份精读，而是说明构成该综述基础的十二项原始工作各自解决了什么问题。

## 阅读建议

- 想了解整体思路：先看本页的三个反复出现的动作，再看任一专题的索引页。
- 想跟踪某条技术线：从专题索引页的「论文脉络」表进入，它按时间列出该专题内部的依赖关系。
- 想查某一篇具体论文：直接用[[computational-mathematics/paper-notes/catalog|论文总目录]]按编号检索。

## 数据与可复现性

- 题录数据：`scripts/data/tao-zhou-publications.json`，含编号、作者、题名、期刊、年份、状态、专题与 DOI。
- DOI 解析与核对：`scripts/enrich-publication-dois.mjs` 与 `scripts/verify-publication-dois.mjs`，均通过 Crossref 查询并逐条比对题名。
- 总目录页的表格由 `scripts/generate-publication-catalog.mjs` 从上述数据生成，正文部分手工撰写。
- 全部插图由 `scripts/generate-diagram-svgs.mjs` 生成双语 SVG 对。
