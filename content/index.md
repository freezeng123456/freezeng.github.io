---
title: Freezeng 知识库
description: 面向计算数学、机器学习与推荐系统理论的个人知识库
lang: zh
translation: en
tags:
  - 导航
---

# 从问题出发，构建可复用的知识

这里不是按时间堆叠的博客，而是一座持续生长的知识库。每个专题尽量回答三件事：**问题为什么难、方法靠什么成立、实现时最容易忽略什么**。

## 三个入口

| 领域                                          | 研究内容                                                                                                                                                               | 知识整理                                                                                                                                                                                                                 |
| --------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [[computational-mathematics/index\|计算数学]] | [[computational-mathematics/research/index\|个人科研]]：[[computational-mathematics/research/rsse\|RSSE]]、[[computational-mathematics/research/apolarity\|Apolarity]] | [[computational-mathematics/knowledge-notes/index\|知识整理]]：[[computational-mathematics/knowledge-notes/time-parallelization/index\|时间并行方法]]；[[computational-mathematics/paper-notes/index\|周涛研究工作精读]] |
| [[计算机科学/index\|机器学习]]                | [[计算机科学/个人科研/MCoTTA\|MCoTTA]]                                                                                                                                 | [[计算机科学/知识整理/index\|知识整理]]：[[计算机科学/知识整理/efficient-llm-inference\|高效大模型推理]]、[[计算机科学/知识整理/how-gpus-work\|GPU 如何工作]]                                                            |
| [[推荐系统理论/index\|推荐系统理论]]          | —                                                                                                                                                                      | [[computational-mathematics/paper-notes/推荐系统/长序列建模\|长序列建模]]；工程专题：[[推荐系统理论/OneRec/index\|OneRec 工业实现链路]]                                                                                  |

## 阅读约定

- **研究页**只写 insight、算法骨架、核心结果与适用边界，不公开尚未发表的推导细节和实现诀窍。
- **知识页**强调从公式到程序的映射，数值结果尽量附带参数、误差指标与复现实验入口。
- **OneRec 页**是临时内部评审稿；仅保留架构与工程逻辑，剔除凭据、内网地址、人员信息和真实业务样本。

> [!tip] 推荐阅读路径
> 如果第一次来，可先看 [[computational-mathematics/knowledge-notes/time-parallelization/index\|时间并行方法地图]]，再看 [[computational-mathematics/paper-notes/推荐系统/长序列建模\|推荐系统长序列建模]]。前者展示如何组织数学理论，后者展示如何把算法路线、工程预算与实验验证放入同一套框架。

## 当前建设状态

- 计算数学：个人科研与知识整理已采用独立目录；时间并行笔记及其数值证据已按学术写作规范整理。
- 时间并行：所有正式 Python 结果产物均已归入对应的 Acta Numerica 章节。
- 机器学习：已整理 MCoTTA 的最小算法改动与实验边界，并新增高效 LLM 推理和 GPU 工作机制两篇可视化知识笔记。
- 推荐系统：新增长序列建模的双语技术地图、混合架构、参数手册和选型框架；OneRec 已形成数据、训练、语义 ID、在线召回、参数与稳定性闭环，等待内部评审后迁移或关闭。
