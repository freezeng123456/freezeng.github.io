---
title: Freezeng 知识库
description: 计算数学、计算机科学与推荐系统理论的个人知识库
tags:
  - 导航
---

# 从问题出发，建立可复用的知识

这里不是按时间堆叠的博客，而是一座持续生长的知识库。每个专题尽量回答三件事：**问题为什么难、方法靠什么成立、实现时最容易忽略什么**。

## 三个入口

| 领域                                          | 研究内容                                                                                                            | 知识整理                                                                                                     |
| --------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| [[计算数学/index\|Computational Mathematics]] | [[计算数学/个人科研/index\|Research]]: [[计算数学/个人科研/RSSE\|RSSE]]、[[计算数学/个人科研/Apolarity\|Apolarity]] | [[计算数学/知识整理/index\|Knowledge Notes]]: [[计算数学/知识整理/时间并行计算/index\|Time Parallelization]] |
| [[计算机科学/index\|计算机科学]]              | [[计算机科学/个人科研/MCoTTA\|MCoTTA]]                                                                              | 暂留空，后续扩展                                                                                             |
| [[推荐系统理论/index\|推荐系统理论]]          | —                                                                                                                   | [[推荐系统理论/OneRec/index\|OneRec 工业实现链路]]                                                           |

## 阅读约定

- **研究页**只写 insight、算法骨架、核心结果与适用边界，不公开尚未发表的推导细节和实现诀窍。
- **知识页**强调从公式到程序的映射，数值结果尽量附带参数、误差指标与复现实验入口。
- **OneRec 页**是临时内部评审稿；仅保留架构与工程逻辑，剔除凭据、内网地址、人员信息和真实业务样本。

> [!tip] 推荐阅读路径
> 如果第一次来，可先看 [[计算数学/知识整理/时间并行计算/index\|Time Parallelization method map]]，再看 [[推荐系统理论/OneRec/端到端链路\|OneRec 端到端链路]]。前者展示如何组织理论，后者展示如何组织大型系统实现。

## 当前建设状态

- Computational Mathematics: research and knowledge notes now use separate directories; the time-parallelization notes and their numerical evidence are written in academic English.
- Time parallelization: all formal Python result artifacts are assigned to the corresponding Acta Numerica chapter.
- 计算机科学科研：已整理 MCoTTA 的最小算法改动与实验边界。
- OneRec：已形成数据、训练、语义 ID、在线召回、参数与稳定性闭环，等待内部评审后迁移或关闭。
