---
title: Freezeng Knowledge Base
description: A personal knowledge base for computational mathematics, machine learning, and recommender systems
lang: en
translation: index
tags:
  - Navigation
---

# Building Reusable Knowledge from Problems

This is not a chronological blog. It is a growing knowledge base in which each topic aims to answer three questions: **Why is the problem difficult? What makes the method work? What is most easily overlooked in implementation?**

## Three Starting Points

| Area                                                              | Research                                                                                                                                                                        | Knowledge Notes                                                                                                                                                                                                                |
| ----------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [[en/computational-mathematics/index\|Computational Mathematics]] | [[en/computational-mathematics/research/index\|Research]]: [[en/computational-mathematics/research/rsse\|RSSE]], [[en/computational-mathematics/research/apolarity\|Apolarity]] | [[en/computational-mathematics/knowledge-notes/index\|Knowledge Notes]]: [[en/computational-mathematics/knowledge-notes/time-parallelization/index\|Time Parallelization]]                                                     |
| [[en/computer-science/index\|Machine Learning]]                   | [[en/computer-science/research/mcotta\|MCoTTA]]                                                                                                                                 | [[en/computer-science/knowledge-notes/index\|Knowledge Notes]]: [[en/computer-science/knowledge-notes/efficient-llm-inference\|Efficient LLM Inference]], [[en/computer-science/knowledge-notes/how-gpus-work\|How GPUs Work]] |
| [[en/recommender-systems/index\|Recommender Systems]]             | —                                                                                                                                                                               | [[en/recommender-systems/long-sequence-modeling\|Long-Sequence Modeling]]; production study: [[en/recommender-systems/onerec/index\|The OneRec Production Pipeline]]                                                           |

## Reading Conventions

- **Research pages** present the insight, algorithmic skeleton, core results, and limits of applicability. They intentionally withhold unpublished derivations and implementation know-how.
- **Knowledge pages** emphasize the mapping from formulas to programs. Numerical results include parameters, error metrics, and reproducibility entry points whenever possible.
- **OneRec pages** are temporary internal-review drafts. They retain architecture and engineering logic while excluding credentials, internal addresses, identities, and real business samples.

> [!tip] Suggested route
> For a first visit, start with the [[en/computational-mathematics/knowledge-notes/time-parallelization/index\|time-parallelization method map]], then read [[en/recommender-systems/long-sequence-modeling\|long-sequence modeling for recommender systems]]. The former shows how mathematical theory is organized; the latter places algorithmic routes, engineering budgets, and experimental evidence in one framework.

## Current Status

- Computational Mathematics: research and knowledge notes use separate directories; the time-parallelization notes and their numerical evidence follow an academic presentation.
- Time parallelization: all formal Python result artifacts are assigned to the corresponding Acta Numerica chapter.
- Machine Learning: the minimal MCoTTA algorithmic change and its experimental limits are documented, together with visual notes on efficient LLM inference and GPU execution.
- Recommender Systems: a bilingual map of long-sequence methods now includes a hybrid architecture, parameter reference, and selection framework. The OneRec loop connecting data, training, semantic IDs, online retrieval, parameters, and reliability remains documented for internal review, after which it should be migrated or closed.
