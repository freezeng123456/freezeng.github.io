---
title: FBSDEs and Stochastic Optimal Control
description: Twenty-two papers that avoid high-dimensional grids with a probabilistic representation and concentrate accuracy in conditional expectations
lang: en
translation: computational-mathematics/paper-notes/fbsde-and-control
tags:
  - computational-mathematics
  - paper-notes
  - stochastic-optimal-control
---

This topic holds 22 papers spanning 2014 to 2026, in roughly two phases. From 2014 to 2023, mostly with Weidong Zhao, the subject is **constructing high-order schemes for forward-backward stochastic differential equations and building a stability theory for them**. From 2025, mostly with Wei Cai, the subject is **using martingale properties to push the same class of problems into very high dimension**.

![How multistep schemes raise the backward accuracy](assets/diagrams/tao-zhou-papers/en/fbsde-multistep.svg)

## Why a probabilistic representation

Consider the coupled Markovian forward-backward system

$$
X_t=X_0+\int_0^t b(s,X_s,Y_s,Z_s)\,\mathrm ds+\int_0^t\sigma(s,X_s,Y_s,Z_s)\,\mathrm dW_s,
$$

$$
Y_t=\xi+\int_t^T f(s,X_s,Y_s,Z_s)\,\mathrm ds-\int_t^T Z_s\,\mathrm dW_s,
$$

with $\xi=\varphi(X_T)$. The system is decoupled when $b$ and $\sigma$ do not depend on $(Y_t,Z_t)$. The nonlinear Feynman-Kac relation ties $(Y_t,Z_t)$ to the solution of a semilinear parabolic equation and its gradient, so this framework supplies a solution route that needs **no spatial grid**: every quantity is a conditional expectation along simulated paths. The cost moves accordingly, and accuracy becomes a question about how those conditional expectations are discretised.

## Four close readings

| Close reading                                                                                                                    | Papers                     | Technical core                                   |
| -------------------------------------------------------------------------------------------------------------------------------- | -------------------------- | ------------------------------------------------ |
| [[en/computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes\|Multistep schemes]]                   | 8, 18, 23, 33, 35, 61, 68  | raising temporal order with future levels        |
| [[en/computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control\|Second-order FBSDEs and control]]  | 16, 19, 25, 26, 41, 50, 51 | fully nonlinear equations and control iterations |
| [[en/computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes\|Stability theory for discretisations]] | 47, 63                     | stability, consistency and convergence           |
| [[en/computational-mathematics/paper-notes/fbsde-and-control/martingale-deep-learning\|Martingale deep learning]]                | 86, 93, 96, 97, 100, 108   | rewriting the residual as a martingale property  |

![Rewrite the equation residual as a martingale property](assets/diagrams/tao-zhou-papers/en/martingale-training.svg)

## Core idea of each paper

### Multistep schemes

- **8 (high-order multistep schemes for coupled FBSDEs)** poses and answers a concrete question: **if the forward stochastic differential equation is solved by the Euler method alone, can the backward equation still be high-order accurate?** Earlier high-order methods for decoupled systems relied on high-order schemes for both directions, and high-order schemes for the forward equation are computationally heavy and often awkward to implement; in the coupled case the forward coefficients depend on the backward unknowns, which is harder still. The paper's licence is a theorem about the generator of the diffusion: if another diffusion's coefficients merely match at the **left endpoint**,

  $$
  \bar b(t_0,\bar X_{t_0};t_0,x_0)=b(t_0,x_0),
  \qquad
  \bar\sigma(t_0,\bar X_{t_0};t_0,x_0)=\sigma(t_0,x_0),
  $$

  then the first derivative of the conditional expectation at $t_0$ is the same. That permits replacing the true forward diffusion by a **frozen-coefficient** (Euler) one inside the one-step derivative approximation, moving the entire high-order requirement onto the backward direction. The derivative weights are fixed by the moment conditions

  $$
  \sum_{i=0}^{k}\alpha_{k,i}\frac{(\Delta t_i)^j}{j!}=\delta_{j1},\qquad j=0,1,\dots,k .
  $$

- **18 (multistep schemes for FBSDEs with jumps)** carries the idea to jump processes.
- **23 (deferred correction methods)** raises accuracy by deferred correction, an alternative to multistep interpolation.
- **33 (explicit theta-schemes for mean-field backward equations)** and **61 (an explicit multistep scheme for mean-field FBSDEs)** handle the mean-field case, where the coefficients depend on the law of the solution itself.
- **35 (explicit deferred correction for second-order FBSDEs)** applies deferred correction to the second-order setting.
- **68 (strong stability preserving multistep schemes)** imports the strong-stability-preserving concept from ODE numerics into FBSDE multistep schemes.

### Second-order FBSDEs and control

- **16 (probabilistic high-order schemes for fully nonlinear parabolic PDEs)** and **19 (high-order schemes for second-order FBSDEs)** treat second-order (fully nonlinear) systems in the sense of Cheridito, Soner, Touzi and Victoir, which correspond to fully nonlinear parabolic equations. A second-order process $\Gamma$ appears alongside $Y$ and $Z$, so there is one more object to discretise.
- **25 (spectral sparse grid approximations for multi-dimensional FBSDEs)** uses sparse grids for the spatial representation, targeting the growth in the number of evaluation points for conditional expectations in several dimensions.
- **26 (an efficient gradient projection method for stochastic optimal control)** and **41 (highly accurate schemes for stochastic optimal control via FBSDEs)** treat stochastic optimal control, the first organising the control iteration as a gradient projection and the second bringing high-order FBSDE schemes to control problems.
- **50 (data-driven feedback control)** and **51 (a Gauss-Seidel type method for dynamic nonlinear complementarity problems)** sit at the edge of the family, handling data-driven feedback control and dynamic problems with complementarity constraints respectively.

### Stability theory for discretisations

- **47 (a unified probabilistic discretisation scheme for FBSDEs)** is the theoretical core of this line. It sets up notions of stability and consistency for FBSDE discretisations so that convergence follows from those two, as in ODE numerical analysis, rather than being proved scheme by scheme.
- **63 (Sinc-theta schemes for backward stochastic differential equations)** uses Sinc approximation for the conditional expectations, an alternative to Gauss-Hermite quadrature.

### Martingale deep learning

- **86 (SOC-MartNet)** treats the Hamilton-Jacobi-Bellman equation with the specific goal of **not computing $\inf_u H$ explicitly**. The residual is rewritten as a martingale property of a process along paths: the residual vanishes if and only if that process is a martingale. A test network probes the property with a family of measurable functions and picks the most violated direction, a value network minimises the detected violation, and a control network updates the feedback control inside the same loop.
- **93 (deep random difference method for high-dimensional quasilinear parabolic equations)** and **96 (martingale deep learning for very high dimensional quasi-linear equations and stochastic optimal controls)** push the same framework to higher dimension and a broader equation class.
- **97 (DeepSPoC)** implements sequential propagation of chaos with deep learning.
- **100 (a derivative-free localised stochastic method for very high dimensional semi-linear parabolic PDEs)** removes the dependence on derivatives.
- **108 (deep policy iteration for high-dimensional mean-field games)** handles mean-field games through a regenerative reformulation.

> [!note] Coverage status
> Papers 8, 16, 18, 19, 23, 25, 26, 33, 35, 41, 47 and 50 have close-reading content checked equation by equation. The close readings for 51, 61, 63, 68, 86, 93, 96, 97, 100 and 108 are still being filled in.

## Two judgements that run through the topic

### High-order accuracy can be added in one direction only

The conclusion of paper 8 generalises: when the two directions of a coupled system carry different accuracy requirements, first check whether the error from the cheaper direction actually reaches the quantity of interest. The generator theorem there is exactly such a check — it shows that at the level of the one-step derivative approximation the forward diffusion need only match at the left endpoint, so using Euler forward does not pollute the backward order.

### Replace "the residual vanishes" with a testable property

The value of the martingale route is that it changes how the residual is tested. Pointwise residual testing is unavailable in high dimension because there is no grid, whereas testing whether a process is a martingale needs only expectations along simulated paths, so dimension no longer enters the cost directly. That is the same conversion as "represent a distribution by a sampleable model" in the [[en/computational-mathematics/paper-notes/scientific-machine-learning/index|scientific machine learning topic]]: swap a condition that cannot be verified directly for an equivalent one that sampling can verify.

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]; per-paper references appear at the end of each close-reading page.
