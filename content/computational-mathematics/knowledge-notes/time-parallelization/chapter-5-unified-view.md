---
title: 第五章：统一视角与方法选择
description: 共同代数解释、实验清单与可复现性边界
lang: zh
translation: en/computational-mathematics/knowledge-notes/time-parallelization/chapter-5-unified-view
tags:
  - 时间并行
  - 方法论
---

## 共同代数形式

许多 PinT 迭代都可写成

$$
U^{k+1}=U^k+M^{-1}(b-AU^k),
$$

其中 $A$ 是全时间耦合时空算子，$M^{-1}$ 是其逆算子的并行近似。

| 方法     | $M^{-1}$ 的来源              |
| -------- | ---------------------------- |
| Parareal | 由粗传播子定义的分块下三角逆 |
| MGRiT    | 时间多层 cycle               |
| STMG     | 完整时空网格上的多层 cycle   |
| ParaDiag | 对循环时间近似执行 FFT 求逆  |
| SWR      | 由波形传输耦合的局部时空逆   |

这一形式明确提出三个核心问题：$M^{-1}$ 能削减哪些误差模态？保留哪些物理信息？哪些操作真正并发？

## 方法选择

![时间并行方法选择地图](assets/diagrams/pint/zh/method-selection.svg)

- **强耗散、低到中阶时间积分：** MGRiT 与 STMG 是自然候选，Parareal 可作为简单原型；
- **高阶时间精度：** PFASST 引入配置与 SDC，但调度更复杂；
- **大型线性系统：** 当移位空间求解可扩展时，ParaDiag 有吸引力；当矩阵指数作用高效时，可考虑 ParaExp；
- **输运或波动占优问题：** 特征线传输、相位校正或 SWR 通常比强耗散粗求解器更忠实地保留传播。

## 完整实验清单

Python 复现项目提供 8 个基线实验和 1 个组合论文验证入口；组合入口生成 6 张论文对应图。由此得到的 14 组 SVG/PNG 图及数值结论均已归入第二至第四章。

| Python 输出              | 网页位置         | 机器可读结果                                         |
| ------------------------ | ---------------- | ---------------------------------------------------- |
| `solution_heat_ade`      | 第二章，对流扩散 | [[assets/pint/data/solution_heat_ade.json            | JSON]] |
| `solution_burgers`       | 第二章，Burgers  | [[assets/pint/data/solution_burgers.json             | JSON]] |
| `solution_wave`          | 第二章，波动方程 | [[assets/pint/data/solution_wave.json                | JSON]] |
| `parareal_heat_ade`      | 第四章，Parareal | [[assets/pint/data/parareal_heat_ade.json            | JSON]] |
| `parareal_burgers`       | 第四章，Parareal | [[assets/pint/data/parareal_burgers.json             | JSON]] |
| `mgrit_heat_ade`         | 第四章，MGRiT    | [[assets/pint/data/mgrit_heat_ade.json               | JSON]] |
| `iterative_paradiag_ade` | 第三章，ParaDiag | [[assets/pint/data/iterative_paradiag_ade.json       | JSON]] |
| `stmg_heat_ade`          | 第四章，STMG     | [[assets/pint/data/stmg_heat_ade.json                | JSON]] |
| Figure 3.15 验证         | 第三章，ParaDiag | [[assets/pint/data/figure_3_15_validation.json       | JSON]] |
| Figure 4.5 验证          | 第四章，Parareal | [[assets/pint/data/figure_4_5_validation.json        | JSON]] |
| Figures 4.9-4.10 验证    | 第四章，MGRiT    | [[assets/pint/data/figure_4_10_validation.json       | JSON]] |
| Figure 4.19 验证         | 第四章，STMG     | [[assets/pint/data/figures_4_19_4_20_validation.json | JSON]] |
| Figure 4.20 验证         | 第四章，STMG     | [[assets/pint/data/figures_4_19_4_20_validation.json | JSON]] |

跨实验简表见 [[assets/pint/data/paper_validation_summary.json|paper_validation_summary.json]]。

上游 MATLAB 仓库还包含直接 ParaDiag、对角化 Parareal、ParaExp、SWR、IDC/PIDC 与波动区域分解脚本。这些脚本已经登记，但尚未全部移植到当前 Python 实验接口；Figure 3.15 的波动 ParaDiag 已经完成。清单声称覆盖所有已经生成的 Python 结果产物，而不是复现每一份上游 MATLAB 脚本。

## Formal 复现

```bash
python3.11 -m pip install -e ".[test]"
actapint all --quick --output-dir results/quick
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  actapint paper_validation --output-dir results/paper-full
```

2026 年 7 月 31 日的 formal 运行在新实验服务器上使用 Python 3.11、NumPy 2.4.6、SciPy 1.17.1 与 Matplotlib 3.11.1。论文套件墙钟约 4 分 24 秒，峰值常驻内存低于 280 MiB。代码路径使用 CPU 稀疏分解、Krylov 方法和 FFT。

每个实验写出：

1. 由同一组分析数组生成的可编辑 SVG 与高分辨率 PNG；
2. 一份记录网格、物理参数、容差和指标的 JSON；
3. 当全时间耦合初值需要随机量时，使用确定性随机初始化。

## 解释边界

- 实验测量数值收敛性并复现论文中的部分趋势与数值，不测量时间维强、弱扩展；
- 没有图报告 MPI 进程数、通信量、初始化成本或墙钟加速比；
- 原 MATLAB 随机生成器与 NumPy 不会产生相同初始数组，应比较收敛因子、迭代次数和最终状态；
- 稀疏分解、FFT 排序与 GMRES 归约会在 $10^{-14}$ 至 $10^{-16}$ 附近产生正常差异；
- `MGRiT_Heat_ADE.m` 中无效表达式 `nu=0.002max;` 被解释为 $\nu=0.002$，与上下文分支及论文一致；
- STMG 论文验证保留后向 Euler 与原 MATLAB 残差约定，同时在 JSON 中保存一致的后平滑残差用于诊断。

## 后续实验的最低报告标准

面向性能的 PinT 研究应报告空间与时间离散、细网格参考解、时间区间数、硬件分配、粗细成本比、停止准则、通信与初始化时间、参数扫描以及强或弱扩展。只有“误差随迭代次数变化”足以支持算法收敛结论，但不足以支持并行效率结论。
