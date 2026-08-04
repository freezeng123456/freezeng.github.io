---
title: GFlowGR Method
description: How trajectory sampling, behavior-aware rewards, and GFlowNet losses fine-tune generative recommendation
lang: en
translation: 推荐系统理论/GFlowGR/方法框架
tags:
  - GFlowGR
  - GFlowNet
  - Generative-Recommendation
---

## 1. Objects in Generative Recommendation

A standard GR stack has three parts:

1. an **item tokenizer** that maps each item to a discrete identifier $\{t_l\}_{l\le L}$;
2. a **user prompt** that assembles task text, history, and profile features;
3. a **generative LLM** that autoregressively emits the next-item identifier.

SFT usually minimizes next-token prediction:

$$
\mathcal{L}_{\mathrm{GR}}(U,v)=-\sum_{l=1}^{L}\log P_{\mathrm{GR}}(t_l\mid U,\{t_i\}_{i<l}).
$$

The paper’s running examples use $L=3$. This objective watches only one positive item. It neither distinguishes interaction value nor directly optimizes the quality of a candidate set.

## 2. Why GFlowNets Fit

Write identifier generation as a state trajectory

$$
\tau=(s_0\to s_1\to\cdots\to s_L),\qquad s_0=\emptyset,\quad s_l=\{t_i\}_{i\le l}.
$$

![Item-identifier generation as a GFlowNet trajectory](assets/diagrams/gflowgr/en/trajectory-flow.svg)

GFlowNets require the probability of any complete trajectory to be proportional to its terminal reward: $P(\tau)\propto R(s_L)$. For GR this simplifies in two ways:

- the forward transition is the LLM token probability, $P_F(s_{l+1}\mid s_l)=P_{\mathrm{GR}}(t_{l+1}\mid U,s_l)$;
- the backward transition is deterministic, $P_B(s_l\mid s_{l+1})=1$, because prefix concatenation is unique.

The remaining design work is therefore **where trajectories come from** and **how rewards are defined**.

## 3. Three Modules

![The three modules used during GFlowGR training](assets/diagrams/gflowgr/en/framework.svg)

For a training pair $(U,v)$:

1. the trajectory sampler expands the positive into a set $\{\tau_n\}_{n\le N}$;
2. the reward model scores each terminal state with $R(s_L)$;
3. a DB or TB loss supplies token-level supervision and is added to the SFT loss.

GFlowGR changes **only the fine-tuning objective**. Inference keeps the original GR decoder. A collaborative model (CM) is optional and can serve both sampling and reward estimation.

## 4. Trajectory Sampler

The sampler builds the $N$ trajectories needed for set-wise learning. $\tau_1$ is the positive item $v$; the rest are augmentations.

| Strategy         | Procedure                                        | Best when                            | Risk                                     |
| ---------------- | ------------------------------------------------ | ------------------------------------ | ---------------------------------------- |
| Interaction logs | Take $N-1$ other items from the same session     | Industrial logs are rich             | Exposure bias enters training directly   |
| Random negatives | Draw $N-1$ unobserved items uniformly            | Starting on academic data            | Negatives become too easy later          |
| CM curriculum    | Start from easy low-score negatives, then harden | A reliable CM exists                 | CM noise can distort the curriculum      |
| LLM on-policy    | Use current LLM beam search for top-$(N-1)$      | You want model-error-aware negatives | Costly; early distributions are unstable |

The default public setup uses **CM-based curriculum sampling**: early epochs prefer easy low-score negatives; later epochs introduce harder items closer to user preference.

In practice, score a fixed candidate pool (for example $2N-1$ items) and then select $N-1$ augmentations from that pool, instead of scoring the full catalog every step.

## 5. Behavior-Aware Reward

![Fusion of three reward signals](assets/diagrams/gflowgr/en/reward-model.svg)

For a positive trajectory $\tau$ and an augmented trajectory $\tau'$, the paper uses three composable signals:

1. **Interaction signal $r_a$**: behavior-level values. The paper’s example uses liked $=10$, clicked $=1$, and augmented / unpresented $=0$.
2. **Estimated score $r_{\hat{y}}$**: a collaborative score $\mathrm{CM}(U,\cdot)$, or the LLM’s own normalized generation probability.
3. **Format similarity $r_{\mathrm{sim}}$**: the number of shared tokens with the positive, measuring partial correctness:

$$
r'_{\mathrm{sim}}=\sum_{l=1}^{L}\mathbb{I}(t_l=t'_l),\qquad r_{\mathrm{sim}}=L.
$$

Default fusion:

$$
R(s_L)=r_a+r_{\hat{y}}+r_{\mathrm{sim}}.
$$

A learnable Weighted Sum or a deeper MLP is optional. Ablations show Weighted Sum helps top-5, Deep Integration does not, and removing any signal hurts; the interaction signal matters most.

Assigning zero interaction reward to unobserved items does not discard them. Collaborative scores and token similarity can still express latent relevance, which is central to mitigating exposure bias.

## 6. GFlowNet Training Objective

After substituting $P_F=P_{\mathrm{GR}}$ and $P_B=1$, detailed balance (DB) and trajectory balance (TB) become:

$$
\mathcal{L}_{\mathrm{DB}}(\tau)=\sum_{l=0}^{L-1}\left(\log\frac{F(s_l)\,P_{\mathrm{GR}}(t_{l+1}\mid U,s_l)}{F(s_{l+1})}\right)^2,
$$

$$
\mathcal{L}_{\mathrm{TB}}(\tau)=\left(\log\frac{Z\prod_{l=0}^{L-1}P_{\mathrm{GR}}(t_{l+1}\mid U,s_l)}{R(s_L)}\right)^2.
$$

$F$ is a learnable flow estimator with $F(s_L)=R(s_L)$ under DB; $Z$ approximates $F(s_0)$. The full loss is

$$
\mathcal{L}(U,v,\{\tau_n\}_{n\le N})=\mathcal{L}_{\mathrm{GR}}(U,v)+\lambda\sum_{n=1}^{N}\mathcal{L}_{\mathrm{GFN}}(\tau_n).
$$

Two implementation points matter:

- **Keep $\mathcal{L}_{\mathrm{GR}}$**. Otherwise flow matching can erase basic recommendation competence; large $\lambda$ clearly regresses in the paper.
- **TB usually beats DB**. Intermediate states do not uniquely identify items, so whole-trajectory balance is a better match.

## 7. Train and Inference Flow

Following Algorithm 1 in the paper, the training loop is:

```text
tokenize items → build prompt U
→ sample N trajectories (positive + augmented)
→ compute P_GR, flow F, reward R
→ optimize L_GR + λ Σ L_GFN
```

Inference is unchanged:

```text
build prompt → generate identifier tokens → decode item
```

GFlowGR is therefore a **training-time plug-in**, not a new online decoder. That makes it easy to attach to TIGER, LETTER, or other GR backbones.
