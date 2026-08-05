---
title: Martingale Deep Learning
description: Papers 86, 93, 96, 97, 100 and 108 - replacing "the residual vanishes" with "a process is a martingale"
lang: en
translation: computational-mathematics/paper-notes/fbsde-and-control/martingale-deep-learning
tags:
  - paper-notes
  - stochastic-optimal-control
  - deep-learning
---

> [!note] Coverage of this page
> Papers **86** (_SIAM J. Sci. Comput._ 47(4), 2025), **93** (_J. Comput. Phys._ 555, 2026), **96** (submitted to _SIAM Rev._, [arXiv:2408.14395](https://arxiv.org/abs/2408.14395)), **97** (submitted to SIAM/ASA JUQ), **100** (submitted to _Numer. Math._) and **108** (preprint).

![Rewrite the equation residual as a martingale property](assets/diagrams/tao-zhou-papers/en/martingale-training.svg)

## 86: not computing $\inf_u H$ explicitly

### Three constraints removed

Solving a high-dimensional Hamilton-Jacobi-Bellman equation requires evaluating

$$
\inf_{\kappa\in U}H\bigl(t,x,\kappa,\partial_xv,\partial^2_{xx}v\bigr),
\qquad
H(t,x,\kappa,z,p)=\tfrac12\mathrm{Tr}\bigl(p\,\bar\sigma\bar\sigma^\top(t,x,\kappa)\bigr)
+z^\top\bar\mu(t,x,\kappa)+c(t,x,\kappa)
$$

at every space-time point. When $U\subset\mathbb R^m$ is high-dimensional or $H$ has no closed-form minimiser, that inner minimisation is itself subject to the curse of dimensionality. Most existing deep PDE solvers either need $\inf_u H$ in explicit form or need pathwise time-recursive training. Moreover the standard deep BSDE architecture relies on the **strong** convergence of Euler-Maruyama (order $1/2$) because it uses pathwise properties of $X$.

This paper removes all three constraints at once: no explicit $\inf_u H$, no time-recursive training, and only **weak** convergence (order 1) is needed.

### Step one: replace the pointwise minimum principle by an integral one

Along an **uncontrolled** diffusion $X_t=X_0+\int_0^t\mu\,\mathrm ds+\int_0^t\sigma\,\mathrm dB_s$, define the Hamiltonian process and its accumulation

$$
H^{u,v}_t:=H\bigl(t,X_t,u(t,X_t),\partial_xv(t,X_t),\partial^2_{xx}v(t,X_t)\bigr),
\qquad
\mathcal M^{u,v}_t:=v(t,X_t)+\int_0^tH^{u,v}_s\,\mathrm ds .
$$

Under suitable integrability the paper's lemma states that an optimal control can be obtained from the **integral** condition

$$
\int_0^T\mathbb E\bigl[H^{u,v}_t\bigr]\mathrm dt
=\inf_{\bar u\in\mathcal U_{\rm ad}}\int_0^T\mathbb E\bigl[H^{\bar u,v}_t\bigr]\mathrm dt .
$$

The proof is one line: $\varepsilon_t:=H^{u,v}_t-\inf_\kappa H\ge0$ pointwise, and the display forces $\int_0^T\mathbb E[\varepsilon_t]\mathrm dt\le0$, so $\varepsilon_t=0$ almost everywhere.

**The payoff is stated plainly**: minimisation moves from "minimise over $U$ at every $(t,X_t)$" to "minimise the functional $\bar u\mapsto\int_0^T\mathbb E[H^{\bar u,v}_t]\mathrm dt$", a double integral estimable by Monte Carlo in which $t$ and $x$ can be sampled independently, so it parallelises trivially.

### Step two: the martingale characterisation

Under integrability conditions, for $(u,v)\in\mathcal U_{\rm ad}\times C^{1,2}$,

$$
(\partial_t+\mathcal L)v(t,X_t)=-H^{u,v}_t\ \ \text{a.e.}
\qquad\Longleftrightarrow\qquad
\mathcal M^{u,v}_t=\mathbb E\bigl[\mathcal M^{u,v}_T\mid\mathcal F_t\bigr].
$$

The forward direction is Itô's formula: substituting the equation gives $\mathcal M^{u,v}_t=v(0,X_0)+\int_0^t\partial_xv\,\sigma\,\mathrm dB_s$, a martingale. The converse uses the martingale representation theorem: with $\mathcal M^{u,v}_t=\mathcal M^{u,v}_0+\int_0^tZ_s\mathrm dB_s$, comparison against the Itô expansion gives

$$
Q_t=\int_0^t\bigl\{(\partial_t+\mathcal L)v+H^{u,v}\bigr\}(s,X_s)\,\mathrm ds
=\int_0^t\bigl\{\partial_xv\,\sigma-Z_s\bigr\}\mathrm dB_s,
$$

which is simultaneously a finite-variation process and a continuous martingale with $Q_0=0$, hence $Q\equiv0$.

**This is the paper's central statement: the PDE residual vanishes if and only if that process is a martingale.**

One geometric assumption deserves note: the support of the uncontrolled diffusion must contain that of the optimally controlled one, $\Gamma(X_t)\supset\Gamma(X^*_t)$. The paper gives two ways to arrange it: sample $X_0\sim N(x_0,rI_d)$ with $r>0$ a hyperparameter, or take $X=X^{u_0}$ for a preliminary control $u_0$.

Another remark matters just as much: the whole formulation uses only **expectations and conditional expectations** of $X$, not pathwise properties. Euler-Maruyama therefore enters only through its **weak** rate (order 1), whereas deep-BSDE-type methods depend on the **strong** rate (order $1/2$). That is the theoretical reason for the first-order temporal convergence observed in the experiments.

### Step three: from the martingale to an adversarial min-max

Since $X$ is Markov, $\mathbb E[\mathcal M^{u,v}_T\mid\mathcal F_t]=\mathbb E[\mathcal M^{u,v}_T\mid X_t]$. To avoid computing a conditional expectation, the paper replaces it by a **weak** condition against a family of test functions:

$$
\sup_{\rho\in\mathcal T}\Bigl|\int_0^{T-\Delta t}\mathbb E
\Bigl[\rho(t,X_t)\bigl(\mathcal M^{u,v}_{t+\Delta t}-\mathcal M^{u,v}_t\bigr)\Bigr]\mathrm dt\Bigr|^2=0 .
$$

The justification is the tower property: $\mathbb E[\rho(t,X_t)(\mathcal M_{t+\Delta t}-\mathcal M_t)]=\mathbb E[\rho(t,X_t)\mathbb E[(\mathcal M_{t+\Delta t}-\mathcal M_t)\mid X_t]]$, so holding for **all** $\rho$ forces $\mathbb E[\mathcal M_{t+\Delta t}-\mathcal M_t\mid X_t]=0$, which is exactly the projection property of conditional expectation. That structure is adversarial learning, with $\rho$ as the discriminator.

The augmented Lagrangian and the min-max problem are

$$
L(u,v,\rho,\lambda)=\int_0^T\mathbb E\bigl[H^{u,v}_t\bigr]\mathrm dt
+\lambda\Bigl|\int_0^{T-\Delta t}\mathbb E\Bigl[\rho(t,X_t)
\bigl(\mathcal M^{u,v}_{t+\Delta t}-\mathcal M^{u,v}_t\bigr)\Bigr]\mathrm dt\Bigr|^2,
$$

$$
(u,v)=\lim_{\lambda\to+\infty}\
\arg\min_{(\bar u,\bar v)}\Bigl\{\sup_{\rho\in\mathcal T}L(\bar u,\bar v,\rho,\lambda)\Bigr\},
\qquad
\mathcal V=\{v\in C^{1,2}:v(T,x)=g(x)\},
$$

and after parameterisation the training involves a control network $u_\alpha$, a value network $v_\theta$ and an adversarial network $\rho_\eta$.

### Two constraints built into the architecture

The control network hard-constrains its range to $U=\prod_i[a_i,b_i]$:

$$
u_\alpha(t,x)=a+\frac{b-a}{6}\,\mathrm{ReLU6}\bigl(\psi_\alpha(t,x)\bigr),
\qquad
\mathrm{ReLU6}(y)=\min\{\max\{0,y\},6\},
$$

with a distance penalty to $U$ for general sets. The value network hard-encodes the terminal condition: $v_\theta(T,x)=g(x)$ and $v_\theta(t,x)=\phi_\theta(t,x)$ for $t<T$.

**Both moves match the stance of the [[en/computational-mathematics/paper-notes/scientific-machine-learning/index|scientific machine learning topic]]**: a constraint the architecture can guarantee does not belong in a penalty. The only thing left in a penalty here is the martingale condition, and it carries a $\lambda\to\infty$ limit.

Code is at [sx-fang/MartNet](https://github.com/sx-fang/MartNet).

## 93, 96 and 100: directions in which the framework extends

- **93 (deep random difference method)** targets high-dimensional quasilinear parabolic equations. The "random difference" in the name refers to replacing derivatives by random differences, reducing the dependence on automatic differentiation.
- **96 (martingale deep learning for very high dimensional quasi-linear equations and stochastic optimal controls)** is the comprehensive submission to _SIAM Review_ that systematises the martingale framework and pushes it to very high dimension.
- **100 (a derivative-free localised stochastic method)** removes the dependence on derivatives for very high dimensional semilinear parabolic equations. "Localised" refers to decomposing the global problem into local subproblems to control variance.

## 97 and 108: the mean-field direction

- **97 (DeepSPoC)** implements **sequential propagation of chaos** with deep learning. Propagation of chaos is the classical result connecting a mean-field limit to a finite particle system; "sequential" means updating one representation in sequence rather than simulating many particles simultaneously, which lowers the memory requirement.
- **108 (deep policy iteration for high-dimensional mean-field games)** treats mean-field games by deep policy iteration with a regenerative reformulation. A mean-field game couples a Fokker-Planck equation to the Hamilton-Jacobi-Bellman equation, so the value function and the distribution must evolve together.

> [!note] Coverage status
> Paper 86 has been checked equation by equation, including the published version. The close readings for 93, 96, 97, 100 and 108 are still being filled in: this page gives only the positioning confirmable from titles, abstracts and their relation to the framework of paper 86, without expanding their formulas and theorems.

## The general judgement behind this route

The martingale route changes **how the residual is tested**. Pointwise residual testing is unavailable in high dimension because there is no grid, whereas testing whether a process is a martingale needs only expectations along simulated paths, so dimension no longer enters the cost directly.

Paper 86 carries the conversion furthest, and it can be summarised as three substitutions:

1. pointwise minimum principle to integral minimum principle (removing the curse of dimensionality in $(t,x)$);
2. PDE residual to martingale property (removing the grid);
3. conditional expectation to a weak condition against test functions (removing the conditional-expectation computation).

Each substitution swaps a condition that cannot be handled directly for an equivalent one that sampling can verify. **The third is the most noteworthy**: it converts "compute a conditional expectation" into "train a discriminator network", at the price of inheriting the instability of adversarial training and with the benefit of avoiding conditional-expectation estimation entirely.

## Coverage check

| Item                                          | Paper | Status                                                                   |
| --------------------------------------------- | ----- | ------------------------------------------------------------------------ |
| HJB equation and the Hamiltonian              | 86    | state equation, cost, value function, explicit $H$                       |
| The three constraints removed                 | 86    | explicit $\inf_u H$, time recursion, strong-rate dependence              |
| Integral minimum principle and its proof      | 86    | Hamiltonian process, integral condition, one-line proof, parallelism     |
| Martingale characterisation, both directions  | 86    | Itô direction, representation direction, the finite-variation argument   |
| Geometric assumption and two remedies         | 86    | support containment, random initial data or preliminary control          |
| Weak rather than strong rate                  | 86    | consequence of using only expectations, observed first-order convergence |
| Adversarial weak condition and tower property | 86    | test function family, projection property, role of the discriminator     |
| Augmented Lagrangian and three networks       | 86    | objective, min-max, $\lambda\to\infty$                                   |
| Two architectural hard constraints            | 86    | ReLU6 range constraint, hard-coded terminal condition                    |

## Sources for this page

- W. Cai, S. Fang, and T. Zhou, [_SOC-MartNet: a martingale neural network for the Hamilton-Jacobi-Bellman equation without explicit inf H in stochastic optimal controls_](https://doi.org/10.1137/24M1681033), SIAM J. Sci. Comput. 47(4) (2025), pp. C795-C819 (preprint [arXiv:2405.03169](https://arxiv.org/abs/2405.03169); code [sx-fang/MartNet](https://github.com/sx-fang/MartNet)).
- W. Cai, S. Fang, and T. Zhou, [_Deep random difference method for high-dimensional quasilinear parabolic partial differential equations_](https://doi.org/10.1016/j.jcp.2026.114767), J. Comput. Phys. 555 (2026), 114767.
- W. Cai, S. Fang, W. Zhang, and T. Zhou, _Martingale deep learning for very high dimensional quasi-linear partial differential equations and stochastic optimal controls_, [arXiv:2408.14395](https://arxiv.org/abs/2408.14395), submitted to SIAM Rev.
- K. Du, Y. Xie, T. Zhou, and Y. Zhou, _DeepSPoC: a deep learning based sequential propagation of chaos method_, submitted to SIAM/ASA J. Uncertain. Quantif.
- S. Fang, C. Sheng, B. Su, and T. Zhou, _A derivative-free localized stochastic method for very high dimensional semi-linear parabolic PDEs_, submitted to Numer. Math.
- S. Fang, S. Wang, Z. Wu, H. Zhang, and T. Zhou, _Deep policy iteration for high-dimensional mean-field games with regenerative reformulation_, preprint.
