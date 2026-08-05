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
> Papers **86** (_SIAM J. Sci. Comput._ 47(4), 2025), **93** (_J. Comput. Phys._ 555, 2026), **96** (submitted to _SIAM Rev._, [arXiv:2408.14395](https://arxiv.org/abs/2408.14395)), **97** (submitted to SIAM/ASA JUQ, [arXiv:2408.16403](https://arxiv.org/abs/2408.16403)), **100** (submitted to _Numer. Math._, [arXiv:2510.02635](https://arxiv.org/abs/2510.02635)) and **108** (preprint arXiv:2604.26782).
>
> **All six were checked equation by equation against their full texts**: paper 86 against both the published SIAM version and the arXiv preprint; papers 93, 97, 100 and 108 against their arXiv full texts including appendices; paper 96 across its three retitled versions. This page therefore carries complete derivations, theorems and **experiment tables**. Theorem and experiment numbering follows each paper's own.

![Rewrite the equation residual as a martingale property](assets/diagrams/tao-zhou-papers/en/martingale-training.svg)

All six share one substitution. In high dimension there is no grid, so "check the PDE residual at every point" is unavailable; but "check whether a process is a martingale" needs only expectations along simulated paths, so dimension no longer enters the cost directly. Paper 86 is the first to carry the substitution out completely, and the other five push it along five separate axes: **remove the derivatives** (96), **supply the missing rate** (93), **remove the neural network** (100), **move to forward mean-field equations** (97), and **move to mean-field games** (108).

## 86: not computing $\inf_u H$ explicitly

### The idea

Solving a high-dimensional Hamilton-Jacobi-Bellman equation requires evaluating

$$
\inf_{\kappa\in U}H\bigl(t,x,\kappa,\partial_xv,\partial^2_{xx}v\bigr),
\qquad
H(t,x,\kappa,z,p)=\tfrac12\mathrm{Tr}\bigl(p\,\bar\sigma\bar\sigma^\top(t,x,\kappa)\bigr)
+z^\top\bar\mu(t,x,\kappa)+c(t,x,\kappa)
$$

at every space-time point. When $U\subset\mathbb R^m$ is high-dimensional or $H$ has no closed-form minimiser, that inner minimisation is itself subject to the curse of dimensionality. Most existing deep PDE solvers either need $\inf_u H$ in explicit form or need pathwise time-recursive training. Moreover the standard deep BSDE architecture relies on the **strong** convergence of Euler-Maruyama (order $1/2$) because it uses pathwise properties of $X$.

This paper removes all three constraints at once: no explicit $\inf_u H$, no time-recursive training, and only **weak** convergence (order 1). It does so by three substitutions, each replacing a condition that cannot be checked directly by an equivalent one that can be checked by sampling: the pointwise minimum principle becomes an integral minimum principle, the PDE residual becomes a martingale property, and the conditional expectation becomes a weak condition against a family of test functions.

### Setting

Feedback controls $\mathcal U_{\rm ad}:=\{u:[0,T]\times\mathbb R^d\to U\mid u\ \text{Borel measurable}\}$ with $U\subset\mathbb R^m$. The state and cost are

$$
X^u_t=x_0+\int_0^t\bar\mu\bigl(s,X^u_s,u(s,X^u_s)\bigr)\mathrm ds
+\int_0^t\bar\sigma\bigl(s,X^u_s,u(s,X^u_s)\bigr)\mathrm dB_s,
$$

$$
J(u):=\mathbb E\Bigl[\int_0^Tc\bigl(s,X^u_s,u(s,X^u_s)\bigr)\mathrm ds+g(X^u_T)\Bigr],
\qquad
\text{find }u^*\text{ with }J(u^*)=\inf_{u\in\mathcal U_{\rm ad}}J(u).
$$

The value function $v(t,x)=\inf_uJ(t,x,u)$ is the viscosity solution of the fully nonlinear HJB equation above. The paper's actual target is slightly more general, an **HJB-type equation** $\partial_tv+\mathcal Lv+\inf_{\kappa\in U}H(t,x,\kappa,\partial_xv,\partial^2_{xx}v)=0$ with $v(T,x)=g(x)$, where $\mathcal L$ is a fixed, **uncontrolled** second-order operator. Section 3 lists the four situations it covers: controls only in the drift; controlled volatility with an uncontrolled part $\bar\sigma=\bar\sigma_0+\bar\sigma_1$; a known preliminary control $u_0$ (take $\mathcal L=\mathcal L^{u_0}$ and $H=(\mathcal L^\kappa-\mathcal L^{u_0})v+c$); and the case where $\bar H(t,x,z,p)=\inf_uH$ is explicit, which degenerates to a semilinear or fully nonlinear parabolic equation.

Throughout, $v\in C^{1,2}$ is assumed, justified by a nondegeneracy condition: $\sigma\sigma^\top(t,x)$ is uniformly positive definite on $[0,T]\times\mathbb R^d$.

### Derivation

**Step one: replace the pointwise minimum principle by an integral one.** Verification gives the pointwise minimum principle $H(t,X^*_t,u^*(t,X^*_t),\partial_xv,\partial^2_{xx}v)=\inf_{\kappa\in U}H(\cdots)$, and that condition on $(t,x)$ is precisely one source of the curse of dimensionality. Along an _uncontrolled_ diffusion $X_t=X_0+\int_0^t\mu\,\mathrm ds+\int_0^t\sigma\,\mathrm dB_s$, define the Hamiltonian process and the cost process:

$$
H^{u,v}_t:=H\bigl(t,X_t,u(t,X_t),\partial_xv(t,X_t),\partial^2_{xx}v(t,X_t)\bigr),
\qquad
\mathcal M^{u,v}_t:=v(t,X_t)+\int_0^tH^{u,v}_s\,\mathrm ds .
$$

**Lemma 3.2.** If $\int_0^T\mathbb E\bigl[\bigl|\inf_\kappa H(t,X_t,\kappa,\partial_xv,\partial^2_{xx}v)\bigr|\bigr]\mathrm dt<\infty$ and the pointwise problem admits a solution in $\mathcal U_{\rm ad}$, then an optimal control can be found from the **integral** condition

$$
\int_0^T\mathbb E\bigl[H^{u,v}_t\bigr]\mathrm dt
=\inf_{\bar u\in\mathcal U_{\rm ad}}\int_0^T\mathbb E\bigl[H^{\bar u,v}_t\bigr]\mathrm dt .
$$

The proof is one line: $\varepsilon_t:=H^{u,v}_t-\inf_\kappa H\ge0$ holds pointwise, while the display forces $\int_0^T\mathbb E[\varepsilon_t]\mathrm dt\le0$, hence $\varepsilon_t=0$ almost everywhere for $\mathrm dt\times\mathbb P$.

**The payoff (Remark 3.3) is stated clearly**: the object being minimised changes from "minimise over $U$ at each $(t,X_t)$" to "minimise the functional $\bar u\mapsto\int_0^T\mathbb E[H^{\bar u,v}_t]\mathrm dt$", and the latter is a double integral evaluable by Monte Carlo in which $t$ and $x$ can be sampled independently, so it parallelises trivially.

**Step two: the martingale characterisation.** **Lemma 3.4.** Under the integrability conditions $\int_0^T\mathbb E[|\partial_xv\,\sigma(t,X_t)|^2]\mathrm dt<\infty$, $\int_0^T\mathbb E[|H^{u,v}_t|^2]\mathrm dt<\infty$ and $\mathbb E[|v(T,X_T)|^2]<\infty$, for $(u,v)\in\mathcal U_{\rm ad}\times C^{1,2}$,

$$
(\partial_t+\mathcal L)v(t,X_t)=-H^{u,v}_t\ \ \text{a.e.}
\qquad\Longleftrightarrow\qquad
\mathcal M^{u,v}_t=\mathbb E\bigl[\mathcal M^{u,v}_T\mid\mathcal F_t\bigr].
$$

The forward direction is Itô's formula: $v(t,X_t)=v(0,X_0)+\int_0^t(\partial_t+\mathcal L)v\,\mathrm ds+\int_0^t\partial_xv\,\sigma\,\mathrm dB_s$, so substituting the equation gives $\mathcal M^{u,v}_t=v(0,X_0)+\int_0^t\partial_xv\,\sigma\,\mathrm dB_s$, a martingale. The converse uses the martingale representation theorem: $\mathcal M^{u,v}_t=\mathcal M^{u,v}_0+\int_0^tZ_s\mathrm dB_s$, and comparing with the Itô expansion,

$$
Q_t=\int_0^t\bigl\{(\partial_t+\mathcal L)v+H^{u,v}\bigr\}(s,X_s)\,\mathrm ds
=\int_0^t\bigl\{\partial_xv\,\sigma-Z_s\bigr\}\mathrm dB_s
$$

is simultaneously a finite-variation process and a continuous martingale with $Q_0=0$, hence $Q\equiv0$.

**This is the paper's central statement: the PDE residual vanishes if and only if the process is a martingale.**

One geometric assumption is worth noting: the support of the uncontrolled diffusion must contain that of the optimally controlled one ($\Gamma_t=\Gamma(X_t)\supset\Gamma(X^*_t)$). Remark 3.1 gives two ways to arrange it: sample $X_0\sim N(x_0,rI_d)$ with $r>0$ a hyperparameter, or take $X=X^{u_0}$ for a preliminary control $u_0$.

**Remark 3.6 matters just as much**: the whole formulation uses only **expectations and conditional expectations** of $X$, never its pathwise properties. Hence Euler-Maruyama enters only through its **weak** order ($1$), whereas deep-BSDE-type methods depend on the **strong** order ($1/2$). This is the theoretical reason for the first-order time convergence observed in Section 4.3.

**Step three: from the martingale to an adversarial min-max.** Since $X$ is Markov, $\mathbb E[\mathcal M^{u,v}_T\mid\mathcal F_t]=\mathbb E[\mathcal M^{u,v}_T\mid X_t]$. To avoid computing that conditional expectation (the original DeepMartNet computed it directly), the paper replaces it by a **weak** condition against test functions $\mathcal T:=\{\rho:[0,T]\times\mathbb R^d\to\mathbb R\mid\rho\ \text{smooth and bounded}\}$:

$$
\sup_{\rho\in\mathcal T}\Bigl|\int_0^{T-\Delta t}\mathbb E
\Bigl[\rho(t,X_t)\bigl(\mathcal M^{u,v}_{t+\Delta t}-\mathcal M^{u,v}_t\bigr)\Bigr]\mathrm dt\Bigr|^2=0 .
$$

The justification is the tower property: $\mathbb E[\rho(t,X_t)(\mathcal M_{t+\Delta t}-\mathcal M_t)]=\mathbb E[\rho(t,X_t)\mathbb E[(\mathcal M_{t+\Delta t}-\mathcal M_t)\mid X_t]]$, so holding for **all** $\rho$ forces $\mathbb E[\mathcal M_{t+\Delta t}-\mathcal M_t\mid X_t]=0$ — exactly the projection property of conditional expectation. That structure is adversarial learning: $\rho$ is the discriminator.

The augmented Lagrangian and the resulting min-max problem are

$$
L(u,v,\rho,\lambda)=\int_0^T\mathbb E\bigl[H^{u,v}_t\bigr]\mathrm dt
+\lambda\Bigl|\int_0^{T-\Delta t}\mathbb E\Bigl[\rho(t,X_t)
\bigl(\mathcal M^{u,v}_{t+\Delta t}-\mathcal M^{u,v}_t\bigr)\Bigr]\mathrm dt\Bigr|^2,
$$

$$
(u,v)=\lim_{\lambda\to+\infty}\
\arg\min_{(\bar u,\bar v)}\Bigl\{\sup_{\rho\in\mathcal T}L(\bar u,\bar v,\rho,\lambda)\Bigr\},
\qquad
\mathcal V=\{v\in C^{1,2}:v(T,x)=g(x)\} .
$$

After parameterisation the training involves three networks: the control network $u_\alpha$, the value network $v_\theta$, and the adversarial network $\rho_\eta$.

**Step four: two constraints built into the architecture.** The control network hard-constrains its range to $U=\prod_i[a_i,b_i]$:

$$
u_\alpha(t,x)=a+\frac{b-a}{6}\,\mathrm{ReLU6}\bigl(\psi_\alpha(t,x)\bigr),
\qquad
\mathrm{ReLU6}(y)=\min\{\max\{0,y\},6\},
$$

with a distance penalty $\bar\lambda\int_0^T\mathbb E[\mathrm{dist}(u_\alpha(t,X_t),U)]\mathrm dt$ for general $U$ (Remark 3.7). The value network hard-wires the terminal condition: $v_\theta(T,x)=g(x)$, and $v_\theta(t,x)=\phi_\theta(t,x)$ for $t<T$. The adversarial network is **deliberately shallow but vector-valued**:

$$
\rho_\eta(t,x)=\sin\bigl(W_1t+W_2x+b\bigr)\in\mathbb R^r,
\qquad
\eta=(W_1,W_2,b)\in\mathbb R^r\times\mathbb R^{r\times d}\times\mathbb R^r,
$$

"chosen to be vector-valued to enhance the stability of adversarial training", in the paper's words.

**Both hard constraints match the disposition of this site's [[en/computational-mathematics/paper-notes/scientific-machine-learning/index|scientific machine learning topic]]**: a constraint the architecture can guarantee should not be put into a penalty. The only thing left in a penalty here is the martingale condition, and it comes with a $\lambda\to\infty$ limit.

**Step five: discretisation and training.** On $\pi_N=\{0=t_0<\cdots<t_N=T\}$ the uncontrolled diffusion is advanced by Euler and the cost-process increment by the **trapezoid rule**:

$$
\Delta\mathcal M^{\alpha,\theta}_{n+1}
:=v_\theta(t_{n+1},X_{n+1})-v_\theta(t_n,X_n)
-\tfrac12\bigl(H^{\alpha,\theta}_n+H^{\alpha,\theta}_{n+1}\bigr)\Delta t_n .
$$

The minibatch loss over a random index set $A\subset\{0,\dots,N\}\times\{1,\dots,M\}$ is

$$
L(\alpha,\theta,\eta,\lambda;A)=\frac1{|A|}\sum_{(n,m)\in A}H^{\alpha,\theta,(m)}_n\Delta t_n
+\lambda\bigl|G(\alpha,\theta,\eta;A)\bigr|^2,
\qquad
G:=\frac1{|A|}\sum_{(n,m)\in A}\rho_\eta\bigl(t_n,X^{(m)}_n\bigr)\Delta\mathcal M^{\alpha,\theta,(m)}_{n+1}\Delta t_n .
$$

One round of Algorithm 3.1 draws $A_i$, takes $J$ descent steps in $\alpha,\theta$, then $K$ ascent steps in $\eta$, updating the penalty by $\lambda\leftarrow\min\{\lambda,\bar\lambda+\delta_4|G|^2\}$. **Remark 3.8 is where the method's speed comes from**: the diffusion $X$ is independent of both the control and the value function, so all $M$ Euler paths are generated **once, offline**; and neither the loss nor the training involves **recursion in time**, which is a structural difference from deep BSDE.

Section 3.4 specialises to parabolic equations: for $\partial_tv+\mathcal Lv+f(t,x,v,\partial_xv,\partial^2_{xx}v)=0$, set $\widetilde{\mathcal M}^v_t:=v(t,X_t)+\int_0^tf\,\mathrm ds$; **no control network is needed** and the objective collapses to $\theta^*=\arg\min_\theta\{\max_\eta|\widetilde G(\theta,\eta)|^2\}$ (Algorithm 3.2). Note that $f$ here may depend on $v$ itself, which is more general than $\bar H$.

### Theorems

**Theorem 3.5 (the martingale formulation).** Assume $(u,v)\in\mathcal U_{\rm ad}\times C^{1,2}$ satisfies the integrability conditions (3.7) and (3.12) above. Then an optimal feedback control $u$ satisfying the pointwise minimum principle, together with the value function $v$ satisfying the HJB equation for $t\in[0,T]$ and $x\in\Gamma(X_t)$, can be found from the **two** conditions

$$
\int_0^T\mathbb E\bigl[H^{u,v}_t\bigr]\mathrm dt=\inf_{\bar u\in\mathcal U_{\rm ad}}\int_0^T\mathbb E\bigl[H^{\bar u,v}_t\bigr]\mathrm dt,
\qquad
\mathcal M^{u,v}_t=\mathbb E\bigl[\mathcal M^{u,v}_T\mid\mathcal F_t\bigr],\quad t\in[0,T].
$$

Lemmas 3.2 and 3.4 are its two halves.

> [!warning] This is a characterisation theorem, not a convergence theorem
> The paper gives **no** convergence-rate theorem for the neural-network approximation, **no** generalisation bound, and **no** proof that the min-max training converges. The regularity $v\in C^{1,2}$ is _assumed_ throughout, justified by the nondegeneracy of $\sigma\sigma^\top$. The one quantitative rate claim — first-order weak convergence in $\Delta t$ — is _argued_ in Remark 3.6 from the weak order of Euler-Maruyama and _verified numerically_ in Section 4.3, not proved. That gap is filled later by Theorem 4 of paper 93.

### Numerical experiments

Implementation: Python 3.12 with PyTorch 2.51, distributed data parallel on a node with 8 NVIDIA A100-SXM4-80GB GPUs.

**(i) Linear parabolic (Section 4.1).** $(\partial_t+\tfrac12\Delta_x)v-f=0$ with exact solution $v(t,x)=1+\tfrac1d\sum_i\sin(t+x_i)$, solved on the line segments $S_1=\{se_1\}$ and $S_2=\{s\mathbf 1_d\}$ for $s\in[-1,1]$. Dimensions $d=100,1000,2000$; **running times over five independent runs were 37, 112 and 363 seconds** respectively.

**(ii) Semilinear parabolic (Section 4.2).** A benchmark from the deep BSDE literature: $(\partial_t+\Delta_x)v-|\partial_xv|^2=0$, $v(T,x)=1+g(x)$, with exact solution $v(t,x)=1-\ln\mathbb E[\exp(-g(x+\sqrt2B_{T-t}))]$ and an **oscillatory** terminal condition $g(x)=\frac1d\sum_i\{\sin(x_i-\frac\pi2)+\sin((\epsilon_0+x_i^2)^{-1})\}$. The reference is Monte Carlo with $10^6$ samples.

**(iii) Time convergence rate (Section 4.3).** The variable-coefficient operator $\mathcal L=\sum_i\sin(2x_i)\partial_{x_i}+\tfrac12\sum_i(1+0.5\sin(5t+x_i))^2\partial^2_{x_i}$ with an **Allen-Cahn-type** source $f=v-v^3+\bar f$, $d=100$. **The observed order is $\mathrm{RE}_1=10^{-0.75}N^{-1.01}$, that is first order $\mathcal O(N^{-1.01})$**, consistent with the weak-order argument of Remark 3.6.

**(iv) Nondegenerate HJB without an explicit $\inf_uH$ (Section 4.4).** $(\partial_t+b^\top\partial_x+\epsilon_1\Delta_x)v+\inf_{\kappa\in\mathbb R^d}(2\kappa^\top\partial_xv+c_1|\kappa|^2)=0$, with the corresponding control problem $J(u)=1+\mathbb E[\int_0^Tc_1|u_s|^2\mathrm ds+g(X^u_T)]$ and $X^u_t=X_0+\int_0^t(b+2u_s)\mathrm ds+\int_0^t\sqrt{2\epsilon_1}\mathrm dB_s$. With $c_1=\epsilon_1^{-2}$ the exact solution is $v(t,x)=-\ln\mathbb E[\exp(-g(X^{t,x}_T))]$. Three instances: **HJB-1** ($b=0$, $\epsilon_0=0.1\pi$, $\epsilon_1=1$), **HJB-2** ($b=\mathbf 1_d$, $\epsilon_0=0.3\pi$, $\epsilon_1=0.2$) and **HJB-3**. Dimensions up to $d=10\,000$, with the explicit $\inf_\kappa H$ never used. Section 4.5 tests generalisation over a space-time region (HJB-2, $d=1000$); Section 4.6 a control problem whose terminal-cost minimum is away from the origin; and Section 4.7 a genuinely non-explicit Hamiltonian $\inf_\kappa(2\kappa^\top\partial_xv+c_1|\kappa|^2+\varepsilon\sin(\mathbf 1_d^\top\kappa))$, benchmarked against the $\varepsilon=0$ solution of HJB-2.

**(v) Efficiency against deep BSDE (Section 4.8, Table 1).** HJB-2 on a single A100; SOC-MartNet runs $I=1000$ iterations against deep BSDE's $I=2000$, and deep BSDE **is given the explicit $\inf_uH$**. RE is relative error and RT is running time in seconds.

| $W_h$  | $d$   | RE (deep BSDE) | RE (SOC-MartNet) | RT (deep BSDE) | RT (SOC-MartNet) |
| ------ | ----- | -------------- | ---------------- | -------------- | ---------------- |
| $256$  | 100   | 3.23E-03       | 1.24E-03         | 73             | **48**           |
| $256$  | 300   | 1.18E-03       | 1.14E-03         | 90             | **53**           |
| $256$  | 500   | 1.05E-03       | 2.89E-03         | 116            | **58**           |
| $256$  | 800   | 2.54E-03       | 4.35E-03         | 145            | **73**           |
| $256$  | 1,000 | 1.86E-03       | 5.83E-03         | 170            | **118**          |
| $d+10$ | 100   | 2.86E-03       | 2.94E-03         | 53             | **20**           |
| $d+10$ | 300   | 3.27E-04       | 8.99E-04         | 103            | **51**           |
| $d+10$ | 500   | 6.41E-04       | 6.92E-04         | 184            | **103**          |
| $d+10$ | 800   | 1.28E-03       | 4.00E-03         | 386            | **255**          |
| $d+10$ | 1,000 | 3.77E-03       | 1.11E-03         | 615            | **360**          |

The standard deviations in the same table are: for $W_h=256$, deep BSDE 8.90E-04, 8.19E-04, 5.71E-04, 1.90E-03, 2.08E-03 and SOC-MartNet 2.43E-03, 1.20E-03, 1.22E-03, 1.41E-03, 2.57E-03; for $W_h=d+10$, deep BSDE 1.09E-03, 1.27E-04, 3.58E-04, 1.14E-03, 5.87E-03 and SOC-MartNet 1.67E-03, 1.02E-03, 4.99E-04, 3.78E-03, 4.26E-04.

**What this table establishes is comparable accuracy at roughly half the wall-clock time, without the explicit Hamiltonian minimiser.** What it does not establish is an accuracy advantage: SOC-MartNet is worse in six of the ten rows, and by a factor of three at $d=1000$ with $W_h=256$.

**(vi) Multi-GPU scaling (Section 4.7, Table 2).** HJB-3, $I=6000$, $W_h=d+10$, times in seconds.

| GPUs          | $d=100$ | $d=500$ | $d=800$ | $d=1{,}000$ | $d=2{,}000$ |
| ------------- | ------- | ------- | ------- | ----------- | ----------- |
| $1\times$A100 | 153     | 775     | 1,350   | 1,909       | 5,032       |
| $2\times$A100 | 151     | 430     | 721     | 1,001       | 2,582       |
| $4\times$A100 | 142     | 233     | 393     | 536         | 1,387       |
| $8\times$A100 | 148     | 153     | 231     | 302         | 773         |

The $d=100$ row barely moves with the GPU count (153 to 148), showing that at that size the bottleneck is communication and startup rather than computation; from $d\ge500$ the scaling is close to linear ($d=2000$ drops from 5032 to 773 seconds, a speed-up of 6.5). **That curve is the direct consequence of a loss with no recursion in time.**

> [!warning] The paper's own stated limitation (Section 5)
> "Further work is needed to avoid computing the very large $d\times d$ Hessian matrix of the value function $v(t,x)$ in the martingale in (3.4) for a derivative-free method." **This is exactly the gap that papers 96 and 100 fill.**

### Relation to the others

It builds on **DeepMartNet** (Cai and co-authors, for boundary-value and eigenvalue problems), adding adversarial learning and the control network; SOC-MartNet is its stochastic-control version.

It is complementary to the entire classical line (papers 8, 19, 23, 25, 26, 35, 41, 47, 63, 68): those solve low-dimensional FBSDEs to high order with quadrature-based conditional expectations and proven error estimates, while paper 86 reaches $d=10^4$ with no proven rates. **The bridge is the same underlying object**: the martingale $\mathcal M_t=v(t,X_t)+\int_0^tH_s\,\mathrm ds$ _is_ the BSDE $Y_t=v(t,X_t)$ in disguise, and the proof of Lemma 3.4 is exactly the nonlinear Feynman-Kac argument.

Contrast with the route of papers 26, 41 and 50 (see the [[en/computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control|second-order FBSDEs and control page]]): those go through the Pontryagin maximum principle and an adjoint FBSDE, whereas paper 86 goes through **dynamic programming and the HJB equation**, so it never forms the adjoint BSDE, and it obtains a _feedback_ control $u_\alpha(t,x)$ on a region rather than an open-loop control along a single initial condition. Paper 68 is Shuixin Fang's classical-analysis precursor. Remark 3.6's weak-versus-strong-order argument is the sharpest single conceptual contrast with the deep BSDE family (E-Han-Jentzen; Han-Jentzen-E) and with all pathwise probabilistic schemes.

Code is at [sx-fang/MartNet](https://github.com/sx-fang/MartNet).

## 93: replacing a second-order operator by a first-order random difference

### The idea

Deep solvers for high-dimensional parabolic PDEs fall into two families with complementary defects. **Direct/PINN-type** methods minimise the strong residual $\|\mathcal Dv-f\|^2$ over sampled points and parallelise perfectly across those points, but they need automatic differentiation to build the $d\times d$ Hessian $\partial_{xx}v$, which is the memory and time bottleneck — the comparison the authors quote is that PINNs **run out of memory** for $d\ge10^4$. **SDE/FBSDE-type** methods avoid $\mathcal D\hat v$ altogether but pay by **losing parallelism in time**: for quasilinear problems the sample paths depend on the unknown $Y_t=v(t,X_t)$ and must be re-simulated sequentially during training, and backward schemes are sequential by construction. Existing martingale methods (papers 86, 96, DeepMartNet) fix both defects at once, but they are derived through Itô calculus, which the authors regard as a barrier for readers without a stochastic-analysis background, and they had **no proven convergence rate**.

**The mechanism of this paper fits in one sentence: a second-order differential operator can be reproduced by a first-order difference quotient of $v$ evaluated at exactly two points, at the price of an expectation over a random jump.** No Hessian, no gradient, no automatic differentiation. As a bonus, the whole re-derivation uses only Taylor expansion and elementary moment identities, and it supplies the error estimate the family had been missing.

### Setting

Terminal-value quasilinear parabolic PDE

$$
\mathcal Dv(t,x)=f\bigl(t,x,v(t,x)\bigr),\quad(t,x)\in[0,T]\times\mathbb R^d,
\qquad v(T,x)=g(x),
$$

$$
\mathcal D:=\partial_t+\mu^\top\bigl(t,x,v(t,x)\bigr)\partial_x
+\tfrac12\mathrm{Tr}\Bigl[\sigma\sigma^\top\bigl(t,x,v(t,x)\bigr)\partial_{xx}\Bigr],
$$

with $\mu\in\mathbb R^d$ and $\sigma\in\mathbb R^{d\times q}$. **Quasilinear** means $\mu,\sigma$ may depend on $v$ itself. The associated FBSDE has a **generator not depending on $Z$**, with nonlinear Feynman-Kac relation $Y_t=v(t,X_t)$, $Z_t=\partial_xv(t,X_t)\sigma(t,X_t,Y_t)$.

### Derivation

**Step one: the basic random-difference expansion.** For $F:\mathbb R^q\to\mathbb R$ with bounded fourth derivatives and a random vector $\xi=(\xi_1,\dots,\xi_q)^\top$, Taylor expansion at $z=\sqrt h\,\xi$ followed by taking the expectation gives

$$
\mathbb E\bigl[F(\sqrt h\,\xi)\bigr]
=F(0)+\sum_{k=1}^{3}\frac{h^{k/2}}{k!}\mathbb E\bigl[(\xi^\top\partial_z)^kF(0)\bigr]
+\frac{h^2}{4!}\mathbb E\bigl[(\xi^\top\partial_z)^4F(c\sqrt h\xi)\bigr].
$$

Imposing the **moment conditions**

$$
\mathbb E[\xi_i]=0,\qquad\mathbb E[\xi_i\xi_j]=\delta_{ij},\qquad
\mathbb E[\xi_i\xi_j\xi_k]=0,\qquad\mathbb E\bigl[|\xi_i\xi_j\xi_k\xi_l|\bigr]<\infty
$$

kills the first- and third-order terms and leaves $\mathbb E[F(\sqrt h\xi)]=F(0)+\frac h2\sum_i\partial^2_{z_i}F(0)+O(h^2)$. **The point is that the left-hand side contains only values of $F$, while the right-hand side has produced a second derivative.** Admissible $\xi$ (Remark 1): $\xi\sim N(0,I_q)$, or independent components with the three-point law $\mathbb P(\xi_i=\pm c)=\frac1{2c^2}$, $\mathbb P(\xi_i=0)=1-c^{-2}$ for $c\ge1$.

**Step two: the random difference operator (RDO) for the convection-diffusion operator.** Fix $(t,x)$, set $V(s,z):=v(t+s,\,x+\mu s+\sigma z)$ with coefficients frozen at $(t,x,v(t,x))$, apply the expansion with $F=V(h,\cdot)$, then Taylor-expand $V(h,0)$ and $\sum_i\partial^2_{z_i}V(h,0)$ in $h$:

$$
\mathcal D_hv(t,x):=\mathbb E\!\left[\frac{v\bigl(t+h,\;x+\mu h+\sigma\sqrt h\,\xi\bigr)-v(t,x)}{h}\right]
=\mathcal Dv(t,x)+O(h).
$$

Remark 3 contrasts this with **randomized-smoothing PINNs** (RS-PINNs), which smooth the _network_, $v_\theta(x)=\mathbb E[\phi_\theta(x+\sqrt h\xi)]$, and differentiate via Stein's identity $\partial_xv_\theta=\mathbb E[\frac{\xi}{\sqrt h}\phi_\theta]$, $\partial_{xx}v_\theta=\mathbb E[\frac{\xi\xi^\top-I_d}{h}\phi_\theta]$. **RS-PINNs behave like a central difference while the RDO behaves like an upwind one**, analogous to a material derivative, so the RDO is better suited to convection-dominated problems.

**Step three: the RDM formulation and the strong-form loss.** Substituting $\mathcal D\to\mathcal D_h$ turns the equation into $\mathbb E[R(t,x,\xi;v)]=O(h)$ where

$$
R(t,x,\xi;v):=\frac{v(t+h,\,x+\xi_h)-v(t,x)}{h}-f\bigl(t,x,v(t,x)\bigr),
\qquad
\xi_h:=\mu\bigl(t,x,v(t,x)\bigr)h+\sigma\bigl(t,x,v(t,x)\bigr)\sqrt h\,\xi,
$$

with the strong-form loss $\mathcal L_{\rm rdm}(\hat v)=\int_0^{T-h}\int_{\mathbb R^d}|\mathbb E[R(t,x,\xi;\hat v)]|^2p(t,x)\,\mathrm dx\,\mathrm dt$.

**Step four: the sampling density is the solution of the Fokker-Planck equation.** This step supplies the theoretical account of the "where should one sample" heuristic. Write $\mathcal D=\partial_t+\mathcal L$, let $\epsilon=\hat v-v$ and let $\hat r:=(\partial_t+\mathcal L)\hat v-f(t,x,v)$ be the _linearised_ residual, with $\mathcal L,f$ evaluated at the **exact** $v$. Then $(\partial_t+\mathcal L)\epsilon=\hat r$; pairing with a weight $p$ and using the $L^2$-adjoint $\mathcal L^*$ gives $\partial_t\int\epsilon p=\int\hat rp+\int\epsilon(\partial_tp-\mathcal L^*p)$. **Choosing $p$ to solve the adjoint (Fokker-Planck) problem** $(\partial_t-\mathcal L^*)p=0$, $p(0,\cdot)=\delta_{x_0}$ eliminates the last term and leaves the exact error representation

$$
\epsilon(0,x_0)=\int_{\mathbb R^d}\epsilon(T,x)p(T,x)\,\mathrm dx
-\int_0^T\!\!\int_{\mathbb R^d}\hat r(s,x)p(s,x)\,\mathrm dx\,\mathrm ds .
$$

So $p(s,x)$ _is_ the sensitivity of the error at $(0,x_0)$ to the residual at $(s,x)$; and since $p(t,\cdot)$ is the density of $X_t$ under the forward SDE, optimal sampling means simulating that SDE with weak Euler-Maruyama, using the **same** law of $\xi$ as in the RDO.

**Step five: two equivalences.** Both identifications are proved by direct computation, with **no stochastic calculus** in the RDM direction.

- **Equivalent to the martingale methods.** With $M_t:=v(t,X_t)+\int_0^tf(s,X_s,v(s,X_s))\mathrm ds$, the martingale condition $\mathbb E^x_t[M_{t+h}-M_t]=0$ becomes $\mathbb E^x_t[v(t+h,X_{t+h})]-v(t,x)+\int_t^{t+h}\mathbb E^x_t[f]\mathrm ds=0$. Approximating the conditional law by weak Euler-Maruyama and the integral by the left-rectangle rule, each to $O(h^2)$, gives $h\,\mathbb E[R(t,x,\xi;v)]=O(h^2)$. **In other words the RDM formulation _is_ the discrete martingale condition of papers 86, 96 and DeepMartNet.**
- **Equivalent to the implicit Euler scheme for the FBSDE.** Taking $\mathbb E^x_{t_n}$ of the backward equation on $[t_n,t_{n+1}]$ kills the Itô integral; substituting $Y=v(t,X_t)$ and using the same two approximations gives $v(t_n,x)=\mathbb E[v(t_{n+1},x+\xi_h)]-h\,f(t_n,x,v(t_n,x))$, described in the paper as a special case of the Zhao-Chen-Peng $\theta$-schemes. **Remark 5 writes out the general-$Z$ case**: the $Y$-part becomes $v(t_n,x)=\mathbb E[v(t_{n+1},x+\xi_h)]-hf(t_n,x,v,z)$ and needs a companion $Z$-part, the simplest instance being
  $$
  z(t_n,x)=\frac1h\mathbb E\Bigl[v(t_{n+1},x+\xi_h)\,\sqrt h\,\xi^\top\Bigr].
  $$
  **This is the cleanest bridge in the whole list between the deep-learning papers and the classical FBSDE schemes**: it is the standard "multiply by the Brownian increment and divide by $h$" representation of $Z$.

**Step six: a Galerkin weak form to control the variance.** The strong loss has an expectation _inside_ a square, so an unbiased minibatch estimator with $M$ points and $2K$ draws of $\xi$ has variance

$$
\mathrm{Var}\bigl[\hat{\mathcal L}_{\rm rdm}(\hat v)\bigr]=O\!\Bigl(\tfrac1M\bigl(1+\tfrac1{K^2}\bigr)\Bigr),
$$

which does **not** decay like the standard $1/(MK)$. Testing against $\rho\in\mathcal T$ and absorbing $p(t,x)$ into the law of $X_t$ collapses the double integral into a single expectation:

$$
\min_{\hat v\in\mathcal V}\max_{\rho\in\mathcal T}\bigl|\mathcal L(\hat v,\rho)\bigr|^2,
\qquad
\mathcal L(\hat v,\rho):=\int_0^{T-h}\mathbb E\bigl[\rho(t,X_t)R(t,X_t,\xi;\hat v)\bigr]\mathrm dt,
$$

for which the analogous estimator attains $\mathrm{Var}[|\hat{\mathcal L}|^2]=O(1/(MK))$ at the same $2MK$ residual evaluations. The **unbiased minibatch product** is

$$
|\mathcal L(\hat v,\rho)|^2\approx\mathcal L^\top(\hat v,\rho;A_1)\,\mathcal L(\hat v,\rho;A_2),
\qquad
A_i=N_i\times M_i,\quad M_1\cap M_2=\varnothing,
$$

**and it is the disjointness of the _path_ index sets that makes the product unbiased.**

**Step seven: networks.** The terminal condition is hard-wired: $v_\theta=\phi_\theta(t,x)$ for $t\le t_{N-1}$ and $v_\theta=g(x)$ for $t>t_{N-1}$. The adversary is a shallow **multiscale (MscaleDNN-style) sine network** with a wide output:

$$
\rho_\eta(t,x)=\sin\bigl(\Lambda(W_1t+W_2x+b)\bigr)\in\mathbb R^r,
\qquad
\Lambda(y_1,\dots,y_r)=(c_1y_1,\dots,c_ry_r)^\top,\quad c_i=1+(i-1)c .
$$

Algorithm 1 alternates $J$ descent steps in $\theta$ with $K$ ascent steps in $\eta$, refreshing $r\%$ of the stored paths every $I_0$ iterations (Remark 6: path simulation is sequential in time, so it is done offline and refreshed rarely).

**Step eight: the HJB extension by policy improvement.** For $\inf_{\kappa\in U}\{\mathcal D^\kappa v+c(t,x,\kappa)\}=0$, split into $u(t,x)=\arg\min_\kappa\{\mathcal D^\kappa v+c\}$ and $\mathcal D^uv+c(t,x,u)=0$, apply the RDO to $\mathcal D^u$, and put both into weak form:

$$
\min_{\hat v\in\mathcal V}\ \sup_{\rho\in\mathcal T}\bigl|\mathcal L(\hat u,\hat v,\rho)\bigr|^2,
\qquad
\min_{\hat u\in\mathcal U_{\rm ad}}\ \mathcal L(\hat u,\hat v,\mathbf 1).
$$

**The constant test function $\mathbf 1$ in the control step** is what replaces pointwise minimisation of the Hamiltonian by an **averaged** minimisation over sampled $(t,x)$ — the same device as in papers 86 and 96. Paths in the HJB case are generated by the _controlled_ SDE.

### Theorems

This is the first paper in the Cai-Fang-Zhou martingale line with a **proved convergence rate**.

- **Assumption 1.** $v\in C^{2,4}$; $\hat v$ Borel measurable; polynomial-growth bounds $|\mu(t,x,v)|+|\mu(t,x,\hat v)|\le C_g(1+|x|^{p_\mu})$, similarly for $\sigma$ with $p_\sigma$, and $\sum_{\alpha\in M_{2,4}}|D^\alpha v(t,x)|\le C_g(1+|x|^{p_v})$; $\xi$ satisfies the moment conditions with $\mathbb E[|\xi_i|^{\bar m}]<\infty$, $\bar m:=\max\{2p_v,8\}$.
- **Theorem 1 (local, frozen-coefficient version).** For all $(t,x)$ and $0<h<\min\{1,T-t\}$,
  $$
  \bigl|\hat{\mathcal D}_hv(t,x)-\hat{\mathcal D}v(t,x)\bigr|\le C_{\rm loc}\,h\,\bigl(1+|x|^{\bar p}\bigr),
  \qquad
  \bar p:=p_v+3p_\mu\max\{1,p_v\}+3p_\sigma\max\{2,p_v\},
  $$
  where $\hat{\mathcal D},\hat{\mathcal D}_h$ are built from $\hat\mu=\mu(t,x,\hat v(t,x))$ and $\hat\sigma=\sigma(t,x,\hat v(t,x))$, and $C_{\rm loc}$ is independent of $\hat v,t,h,T,x$.
- **Corollary 2 (local truncation error, set $\hat v=v$).** $|\mathcal D_hv-\mathcal Dv|\le C_{\rm loc}h(1+|x|^{\bar p})$, that is, **the RDO is first-order accurate in $h$**.
- **Assumption 2 (stability hypotheses).** Lipschitz continuity in the third argument, $|\mu(t,x,y_1)-\mu(t,x,y_2)|\le C_\mu|y_1-y_2|$ (with $C_\sigma$ for $\sigma$ and $C_f$ for $f$), plus the **nonlinearity bound** $C_\mu|\partial_xv|+C_\sigma|\sigma||\partial_{xx}v|\le C_{\rm nl}$. If $\mathcal D$ is linear one may take $C_\mu=C_\sigma=0$ and $C_{\rm nl}=0$. **Assumption 3:** the moment bound $\max_n\mathbb E[|X^m_n|^{2\bar p}]\le C_{\rm EM}$ on the Euler-Maruyama paths.
- **Lemma 3 (zero-stability).** With $M^p_n[\psi]:=\int|\psi(t_n,x)|^pP_{X_n}(x)\mathrm dx$ (the mean $p$-th moment against **the law of the sampled paths**), $\Delta v:=v-\hat v$ and $\Delta R:=\mathbb E[R(\cdot;v)]-\mathbb E[R(\cdot;\hat v)]$, for $0<h\le\min\{1,(24C_{\rm nl}^2+12C_f^2)^{-1}\}$,
  $$
  \max_{0\le n\le N-1}M^2_n[\Delta v]\le C_{\rm st1}\exp(C_{\rm st2}T)
  \Bigl\{M^2_N[\Delta v]+\frac TN\sum_{n=0}^{N-1}M^2_n[\Delta R]+Th^2\Bigr\}.
  $$
  The proof uses a **backward discrete Grönwall inequality**.
- **Theorem 4 (global error).** Under Assumptions 1-3 and the same step restriction,
  $$
  \max_{0\le n\le N-1}M^2_n[v-\hat v]\ \le\ C_1\exp(C_2T)
  \Bigl\{M^2_N[v-\hat v]+\mathcal L_{{\rm rdm},\pi}(\hat v)+Th^2\Bigr\},
  $$
  $$
  \mathcal L_{{\rm rdm},\pi}(\hat v):=\frac TN\sum_{n=0}^{N-1}\int_{\mathbb R^d}
  \bigl|\mathbb E[R(t_n,x,\xi;\hat v)]\bigr|^2P_{X_n}(x)\,\mathrm dx,
  $$
  with $C_1=C_{\rm st1}\max\{1,2+4C_{\rm loc}^2(1+C_{\rm EM})\}$ and $C_2=C_{\rm st2}$, both independent of $h,N,\hat v$. Each term has its own home: the first vanishes identically for the network; the second is driven to zero by adversarial training; the third gives **first-order accuracy in time**, $M^2_n[v-\hat v]=O(h^2)$, that is $O(h)$ in the $L^2$ norm.

> [!warning] Accuracy holds only where the sampled paths go
> The norm $M^2_n$ in Theorem 4 is weighted by the law of $X^m_n$, so **it guarantees accuracy only inside the region the sample paths explore**. This limitation is common to the whole family (papers 86, 93, 96, 100, 108). Remark 1 of paper 96 states the same thing as "a good pilot process must cover the region of interest with high probability", and the $d=1$ visualisation of paper 108 shows it directly.

- **Remark 9, the sharpest comparative claim in the paper.** Because the RDM formulation is equivalent to the discrete martingale condition, Theorem 4 **also proves first-order-in-time convergence for papers 86, 96 and DeepMartNet**. This is higher than the $O(h^{1/2})$ typical of deep-BSDE methods for general coefficients. The mechanism: martingale/RDM methods use Euler-Maruyama only to approximate the _conditional law_ of $X_{t+h}$ given $X_t=x$, so the error is governed by the **weak** order $O(h)$; deep BSDE methods use it to approximate _sample paths_, so they are capped by the **strong** order $O(h^{1/2})$.
- What is absent: no error bound in terms of network size or training; no proof that the min-max problem for HJB controls the error of $v$ — **Remark 7 explicitly flags this as open**; no rate in $d$.

### Numerical experiments

Errors are reported as relative $L^1$ and $L^\infty$ errors $\mathrm{RE}_1(t_n)$, $\mathrm{RE}_\infty(t_n)$ over sampled point sets $D_n=\{X^m_n\}$, evaluated at $t=0$ unless stated otherwise. Test curves are $S_2=\{s\mathbf 1_d:s\in[-1,1]\}$ (a straight line whose length grows with $d$) and $S_3=\{l(s)\}$ with $l_i(s)=s\,\mathrm{sgn}(\sin i)+\cos(i+\pi s)$ (a curve winding through $\mathbb R^d$).

**(i) Convection-diffusion with a steep gradient, $d=10^3$ (Section 4.1).** $(\partial_t+\mu^\top\partial_x+\frac{\bar\sigma^2}2\sum_i\partial^2_{x_i})v=0$, $T=2$, $\bar\sigma^2=0.1$, $\mu_i=c\tanh(10x_i)$, terminal $v(T,x)=\frac1d\sum_i\{\tanh(x_i)+\cos(10x_i)\}$ (highly oscillatory); reference by Monte Carlo with $10^6$ paths at EM step $T/100$. Solved for $c=1$ and $c=5$; the steep gradient at $x=0$ is captured, which the authors attribute to the Fokker-Planck-based sampling.

**(ii) Quasilinear parabolic, $d=10^4$ and $10^5$ (Section 4.2, Table 2).** Exact solution $v(t,x)=V((t-0.5)\mathbf 1_d+x)$ with $V(x)=\sum_{i=1}^{d-1}c_iK(x_i,x_{i+1})+c_dK(x_d,x_1)$, $c_i=(1.5-\cos(i\pi/d))/d$, $K(x_i,x_j)=\sin(x_i+\cos(x_j)+x_j\cos(x_i))$, $T=1$. Three cases: **QLP-1** with $\mathcal D=\partial_t+\frac{v^2}2\sum_i\partial^2_{x_i}$, $f=v-v^3+Q$; **QLP-2a** with $\mathcal D=\partial_t+(\frac v2-1)\sum_i\partial_{x_i}+\frac{v^2}2\sum_i\partial^2_{x_i}$, $f=v^2+Q$; **QLP-2b** with the drift scaled by $1/d$ but a **dense** diffusion $\frac1{2d^2}\sum_{i,j,k}\sigma_{ik}\sigma_{jk}\partial_{x_i}\partial_{x_j}$, $\sigma_{ij}=\cos(x_i)+v\sin(x_j)$ — noted as effectively impossible for PINNs, which would need the full Hessian.

| Equation | $d$    | Mean $\mathrm{RE}_1$ | SD      | Mean $\mathrm{RE}_\infty$ | SD      | RT (s) |
| -------- | ------ | -------------------- | ------- | ------------------------- | ------- | ------ |
| QLP-1    | $10^4$ | 1.98E-2              | 7.24E-3 | 3.91E-2                   | 8.75E-3 | 1585   |
| QLP-1    | $10^5$ | 1.29E-2              | 1.40E-3 | 3.20E-2                   | 3.66E-3 | 3761   |
| QLP-2a   | $10^4$ | 2.82E-2              | 1.28E-2 | 9.08E-2                   | 1.94E-2 | 1589   |
| QLP-2a   | $10^5$ | 4.06E-2              | 1.14E-3 | 1.37E-1                   | 1.33E-2 | 3773   |
| QLP-2b   | $10^4$ | 5.77E-2              | 2.00E-3 | 1.12E-1                   | 1.18E-2 | 1602   |
| QLP-2b   | $10^5$ | 5.19E-2              | 7.69E-4 | 1.05E-1                   | 1.09E-2 | 3822   |

**The notable feature is that going from $d=10^4$ to $10^5$ barely changes the error while the time rises only by a factor of about 2.4.** But the accuracy sits at the $10^{-2}$ level, which is consistent with Theorem 4's first-order rate and is an entirely different scale from the sixth-order accuracy of the classical schemes.

**(iii) HJB, $d=10^4$ (Section 4.3, Table 3).** $\partial_tv+\inf_{\kappa\in\mathbb R^d}\{(b+c\sigma\kappa)^\top\partial_xv+\frac12|\kappa|^2\}+\frac12\mathrm{Tr}[\sigma\sigma^\top\partial_{xx}v]=0$, $c=2$, $T=1$, $v(T,x)=\ln(1+\frac1d\sum_ix_i^2+0.5\sin(10x_i))$. The reference comes from the **Cole-Hopf transform** $v(t,x)=-c^{-2}\ln\mathbb E[\exp(-c^2v(T,X^0_T))\mid X^0_t=x]$ by Monte Carlo with $10^6$ samples at EM step $T/100$. Instances: **HJB-1a** with $b_i=\sin(t+i+x_{i+1}-1)$ (cyclic, $x_{d+1}:=x_1$, coupling the components) and $\sigma=0.5I_d$; **HJB-1b** the same with $\sigma=0.025I_d$ (a less smooth solution); **HJB-2** with $b_i=\sin(x_{i+1})$ and $\sigma_{ij}=0.5\delta_{ij}\tanh((t-0.5)^2+d^{-1}\sum_kx_k^2)$ (variable diffusion). Networks are fully connected of width $W=d+10$.

| Equation | Mean $\mathrm{RE}_1$ | SD      | Mean $\mathrm{RE}_\infty$ | SD      | RT (s) |
| -------- | -------------------- | ------- | ------------------------- | ------- | ------ |
| HJB-1a   | 8.76E-3              | 8.09E-4 | 3.76E-2                   | 5.74E-3 | 2709   |
| HJB-1b   | 1.73E-2              | 8.07E-4 | 4.74E-2                   | 1.94E-3 | 2707   |

**(iv) Sampling-strategy ablation, $d=1$ (Section 4.4).** Deliberately low-dimensional for visualisation, on the convection-diffusion problem with $c=1$ and $c=5$: "Dynamics" sampling (by the SDE the Fokker-Planck analysis prescribes) versus "Plain" sampling $X^m_n=X^m_0+B^m_{t_n}$ with $X^m_0\sim U[-1,1]$. The former is the one supported by the error-propagation analysis of Section 2.4.

**(v) Head-to-head against PINNs and RS-PINNs (Section 4.5, Table 4).** On the elliptic benchmarks of Hu et al.: Allen-Cahn $\Delta v+v-v^3=f$ and Sine-Gordon $\Delta v+\sin v=f$ on the unit ball $B_d$ with $v|_{\partial B_d}=0$ and exact solution $v(x)=(1-|x|^2)\sum_{i=1}^{d-1}c_i\sin(x_i+\cos(x_{i+1})+x_{i+1}\cos(x_i))$, $c_i\sim N(0,1)$ i.i.d. The PINN and RS-PINN numbers are quoted verbatim from that literature. The table also includes **strong-form DRDM**, which the authors themselves identify as essentially the shotgun method of Xu-Zhang up to offline path sampling and antithetic variates. Metrics are relative $L^2$ error, runtime in minutes, memory in MB.

| Method                      | Metric        | $d=10^2$ | $10^3$   | $5\times10^3$ | $10^4$   | $10^5$   |
| --------------------------- | ------------- | -------- | -------- | ------------- | -------- | -------- |
| PINNs                       | AC rel. $L^2$ | 7.187E-3 | 5.617E-4 | 1.773E-3      | N.A.     | N.A.     |
| PINNs                       | SG rel. $L^2$ | 7.192E-3 | 5.642E-4 | 1.782E-3      | N.A.     | N.A.     |
| PINNs                       | RT (min)      | 3        | 285      | 1832.4        | N.A.     | N.A.     |
| PINNs                       | Memory (MB)   | 1328     | 4425     | 56563         | >81252   | >81252   |
| RS-PINNs                    | AC rel. $L^2$ | 7.923E-3 | 5.504E-4 | 1.802E-3      | 1.860E-3 | 2.192E-3 |
| RS-PINNs                    | SG rel. $L^2$ | 7.835E-3 | 6.744E-4 | 1.795E-3      | 1.854E-3 | 2.176E-3 |
| RS-PINNs                    | RT (min)      | 1.8      | 7.2      | 31.8          | 66       | 720      |
| RS-PINNs                    | Memory (MB)   | 1413     | 1815     | 3593          | 5789     | 45599    |
| Strong-form DRDM            | AC rel. $L^2$ | 1.311E-2 | 6.171E-3 | 4.714E-3      | 3.112E-3 | 7.423E-4 |
| Strong-form DRDM            | SG rel. $L^2$ | 2.009E-2 | 6.191E-3 | 4.715E-3      | 3.112E-3 | 2.232E-2 |
| Strong-form DRDM            | RT (min)      | 1.04     | 1.11     | 1.81          | 2.70     | 48.64    |
| Strong-form DRDM            | Memory (MB)   | 75       | 254      | 1204          | 2389     | 23739    |
| **Weak-form DRDM (theirs)** | AC rel. $L^2$ | 3.963E-2 | 6.211E-3 | 4.699E-3      | 3.118E-3 | 7.381E-4 |
| **Weak-form DRDM**          | SG rel. $L^2$ | 4.172E-2 | 6.306E-3 | 4.701E-3      | 3.118E-3 | 4.769E-3 |
| **Weak-form DRDM**          | RT (min)      | 1.22     | 1.26     | 1.43          | 1.59     | 13.26    |
| **Weak-form DRDM**          | Memory (MB)   | 57       | 185      | 859           | 1698     | 16829    |

**The authors' own reading is unusually candid: DRDM is generally _less_ accurate than PINNs and RS-PINNs**, which they attribute to the $O(h)$ truncation error of the RDO that the AD-based methods do not incur. The gain is entirely in runtime and memory, and it grows with $d$: at $d=10^5$ the weak-form DRDM needs 13.26 minutes and 16.8 GB, RS-PINNs need 720 minutes and 45.6 GB, and PINNs cannot run at all. Weak and strong forms have comparable accuracy, but per gradient step with minibatch $n_b$ the strong form evaluates $v_\theta$ $n_b+n_b\ell$ times ($\ell=128$ draws of $\xi$ per point in their runs) whereas the weak form evaluates it $2n_b$ times, so the weak form wins as $d$ grows.

### Relation to the others

**It is the theoretical capstone of papers 86 and 96.** Those two introduced the martingale loss and the derivative-free variant respectively, and neither proved a rate; paper 93 shows their discrete martingale condition is algebraically identical to the RDM formulation, so Theorem 4 certifies **first-order-in-time convergence for all three**. The authors' code repository makes the identification explicit: the `DfSocMartNet` class is labelled "DRDM / derivative-free MartNet ... these two methods are equivalent".

**It is also where the deep-learning branch is formally reconnected to the classical FBSDE branch.** Section 2.6 derives the implicit Euler $Y$-part and Remark 5 the $Z$-part, citing Zhang (2004) and Zhao-Chen-Peng (2006) and pointing to papers 8, 47, 63 and to Zhao-Li-Zhang / Zhao-Zhang-Ju for the general theory. In other words **the RDO is nothing but the one-step, weak-Euler member of that family**, differing only in that the conditional expectation is left as an expectation instead of being evaluated by Gauss-Hermite quadrature as in papers 8, 16, 19 and 25. Correspondingly, the classical papers prove high orders ($O(h^k)$ with $k$ up to 6) but only for $d\lesssim10$, whereas paper 93 proves order 1 only, at $d=10^5$.

It is the direct successor of paper 96 and the sibling of paper 100: the random difference operator first appeared in paper 96's v1, and paper 93 is the one that isolates, names, analyses and benchmarks it; paper 100 is the third member of the derivative-free trio, localised rather than global.

## 96: derivative-free, and parallel in time as well as space

### The idea

Paper 86 still needs $\partial_xv$ and, worse, the $d\times d$ **Hessian $\partial^2_{xx}v$** inside the Hamiltonian, computed by automatic differentiation at every sample. For $d\sim10^4$ that is the binding cost, and paper 86's own conclusion flags it as the open problem. A second bottleneck is that all SDE-model-based deep methods (deep BSDE, DeepMartNet, and paper 86 in the _controlled_ case) must simulate paths that depend on the unknown $v$ or $u$, so the paths must be **re-simulated sequentially in time** after every network update.

This paper removes both. **The derivative-free mechanism fits in one sentence: write the martingale increment so that it contains only two evaluations of the network at two nearby points, with no derivative anywhere in between.** The time-parallelism mechanism comes from a structural device: impose the martingale property not on **one** long path spanning $[0,T]$ but on a **family** of short-horizon processes, one per time step, so the time steps decouple.

### Setting

Quasilinear parabolic equation

$$
(\partial_t+\mathcal L)v(t,x)+f\bigl(t,x,v(t,x)\bigr)=0,\quad(t,x)\in[0,T)\times\mathbb R^d,
\qquad v(T,x)=g(x),
$$

$$
\mathcal L:=\mu^\top\bigl(t,x,v(t,x)\bigr)\partial_x
+\tfrac12\mathrm{Tr}\Bigl\{\sigma\sigma^\top\bigl(t,x,v(t,x)\bigr)\partial^2_{xx}\Bigr\}.
$$

Note the **quasilinearity**: $\mu,\sigma$ depend on $v$ itself, which is more general than paper 86's fixed $\mathcal L$. The standing assumption throughout is that "the functions $\mu$, $\sigma$ and $v$ are smooth enough to validate the involved truncation error estimates".

### Derivation

**Step one: pilot process and system process.** This is the paper's key structural device.

- The **pilot process** $\hat X$, used only to explore $\mathbb R^d$ and built from an _initial guess_ $\hat v$:
  $$
  \hat X_t=\hat X_0+\int_0^t\hat\mu(s,\hat X_s)\mathrm ds+\int_0^t\hat\sigma(s,\hat X_s)\mathrm dB_s,
  \qquad
  \hat\mu(t,x):=\mu(t,x,\hat v(t,x)),\ \ \hat\sigma(t,x):=\sigma(t,x,\hat v(t,x)).
  $$
- The **system process** $X^s_t$ for $0\le s\le t\le T$, **started at $\hat X_s$** and generated by the true $\mathcal L$:
  $$
  X^s_t=\hat X_s+\int_s^t\mu(r,X^s_r,v(r,X^s_r))\,\mathrm dr+\int_s^t\sigma(r,X^s_r,v(r,X^s_r))\,\mathrm dB_r .
  $$

**Step two: the martingale formulation, localised in time.** Itô applied to $t\mapsto v(t,X^s_t)$ gives

$$
\mathcal M^s_t:=v(t,X^s_t)-v(s,X^s_s)+\int_s^tf\bigl(r,X^s_r,v(r,X^s_r)\bigr)\mathrm dr
=\int_s^tR(r,X^s_r;v)\,\mathrm dr+\int_s^t(\partial_xv)^\top\sigma\,\mathrm dB_r,
$$

where

$$
R(t,x;v):=(\partial_t+\mathcal L)v(t,x)+f\bigl(t,x,v(t,x)\bigr)
$$

is **exactly the PDE residual**. Conditioning kills the Itô integral, and setting $t=s+h$ gives $\mathbb E[\mathcal M^s_{s+h}\mid\hat X_s]=h\,R(s,\hat X_s;v)+O(h^2)$, so the martingale formulation reads

$$
\mathbb E\bigl[\mathcal M^t_{t+h}\,\big|\,\hat X_t\bigr]=0,
\qquad 0\le t\le T-h .
$$

**The point of this formulation is that the residual is never computed explicitly, yet it is exactly what the condition characterises.**

> [!warning] Remark 1: the criterion only acts where the pilot process goes
> The condition only guarantees that the residual $R(t,\hat X_t;v)$ vanishes **in the region explored by the pilot process**, so the pilot must cover the region of interest with high probability. **This is the practical risk shared by the whole family (papers 86, 93, 96, 100, 108)**: it appears in paper 93 as the weighted norm $M^2_n[\cdot]$ and is directly visible in paper 108's visualisation.

**The crucial difference from DeepMartNet and paper 86** is that there, the martingale property is imposed on **one** long path $t\mapsto X^0_t$ over all of $[0,T]$; here it is imposed on a **family** of short-horizon system processes $X^t_{t+h}$, one per time step, which **decouples the time steps** and permits minibatch sampling over time and full parallelism.

**Step three: removing the conditional expectation with Galerkin plus adversarial training.** The tower property gives $\mathbb E[\rho(t,\hat X_t)\mathbb E[\mathcal M^t_{t+h}|\hat X_t]]=\mathbb E[\rho(t,\hat X_t)\mathcal M^t_{t+h}]$, and one Euler step approximates $\mathcal M^t_{t+h}$ by the fully explicit, **derivative-free** quantity

$$
\mathcal M(t,x,w;v):=v\bigl(t+h,\ x+\mu(t,x)h+\sigma(t,x)\sqrt h\,w\bigr)-v(t,x)+h\,f\bigl(t,x,v(t,x)\bigr),
\qquad \xi\sim N(0,I_q).
$$

This yields $\min_{v\in\mathcal V}\sup_{\rho\in\mathcal T}|G(v,\rho)|^2=0$ with $G(v,\rho):=\int_0^{T-h}\mathbb E[\rho(t,\hat X_t)\mathcal M(t,\hat X_t,\xi;v)]\mathrm dt$ and $\mathcal V=\{v\in C^{1,2}:v(T,x)=g(x)\}$. **Remark 2** gives the approximation error of this step as $O(h^2)$, from the weak-second-order _local_ truncation error of Euler-Maruyama plus the locally second-order left-rectangle quadrature.

**Note that $\mathcal M$ contains no derivative of $v$ at all**, only two evaluations of the network at two nearby points. That is the entire derivative-free mechanism.

**Step four: an unbiased minibatch estimate of a squared expectation.** Pilot paths come from Euler, the empirical loss is $G(v,\rho;A):=\frac h{|A|}\sum_{(n,m)\in A}\rho(t_n,\hat X^m_n)\mathcal M(t_n,\hat X^m_n,\xi^m_n;v)$, and then

$$
|G(v,\rho)|^2\ \approx\ G(v,\rho;A_1)\,G(v,\rho;A_2),
\qquad A_1\cap A_2=\varnothing .
$$

**Using two disjoint index sets is what makes the estimate of the _square_ unbiased**; the naive $|G(\cdot;A)|^2$ is biased upward by the minibatch variance.

**Step five: networks.** The value network hard-wires the terminal condition. The **multiscale adversarial network** is richer than paper 86's:

$$
\rho_\eta(t,x)=\sin\bigl(\Lambda(W_1t+W_2x)+b\bigr)\in\mathbb R^r,
\qquad
\Lambda(y_1,\dots,y_r)=(cy_1,2cy_2,\dots,rcy_r)^\top,
$$

the scale layer $\Lambda$ following the multiscale-neural-network idea and giving the adversary a spread of frequencies. The control network is the same ReLU6 construction as in paper 86.

**Step six: the HJB extension via policy improvement.** The target is $\partial_tv+\inf_{\kappa\in U}\{\mathcal L^\kappa v+c(t,x,\kappa)\}=0$, split into a policy-improvement pair. Given $u$ the first equation is **linear**, so steps two to four apply. For the control step, the pointwise $\arg\min$ is replaced by its integral version; since $\mathbb E[\mathcal M(t,\hat X_t,\xi;u,v)]=h\{(\partial_t+\mathcal L^u)v+c\}+O(h^2)$ and neither $h^{-1}$ nor $\partial_tv$ affects the minimiser, the control step collapses to

$$
\min_{u\in\mathcal U_{\rm ad}}\ G(u,v,\mathbf 1),
$$

that is, the constant test function $\rho\equiv1$. **This is a derivative-free implementation of the minimum principle**: find the optimal feedback control by minimising the mean of the value function instead of minimising the Hamiltonian pointwise.

> [!note] Remark 4: an honest caveat
> Relative to paper 86 the loss is genuinely free of $\partial_xv$ and $\partial^2_{xx}v$. **But in the controlled case** the random jump $\xi^{t,x,u}_h:=\mu(t,x,u(t,x))h+\sigma(t,x,u(t,x))\sqrt h\,\xi$ depends on $u$, so **those jumps cannot be pre-computed before training**; only the pilot paths $\hat X$ can.

Version v1 derives the same formulation from a **random finite difference operator** and observes that the derivative-free martingale network for a semilinear parabolic equation is equivalent to a **weighted Galerkin-type weak adversarial network** with weight equal to the transition density $p(t,x)$ of $X_t$, that is, a weak-form PINN. This is the conceptual bridge to paper 93.

### Theorems

The paper is **constructive and algorithmic**. Its "results" are the martingale formulation, the derivative-free weak form, the derivative-free minimum principle, and the two algorithms.

- The only quantitative estimate is **Remark 2's $O(h^2)$ local approximation error**, from the weak-second-order local truncation error of Euler-Maruyama plus second-order left-rectangle quadrature. Over $O(1/h)$ steps this is the usual mechanism behind first-order global weak accuracy, but **the paper states no global convergence theorem**.
- Absent: no convergence theorem for the neural approximation, no error bound in terms of network size, no proof that the adversarial training converges, no rate in $d$. **The temporal-order part of that gap is filled by Theorem 4 of paper 93.**

### Numerical experiments

Both $u_\alpha$ and $v_\theta$ are fully connected with **4 hidden layers** of width $W$; $M=10^4$ pilot paths; minibatch $|M_i|=256$, $128$, $64$ for $d\le1000$, $d=2000$ and $d=10^4$ respectively.

**(i) Allen-Cahn, $d=100$ (Section 4.1).** $\partial_tv+\Delta_xv+v-v^3=0$, $v(T,x)=1/(2+0.4|x|^2)$, $T=0.3$, evaluated at $x_0=0$. The reference $v(0,x_0)\approx0.0528$ comes from the branching-diffusion method. The SIAM-formatted version reports for this test a mean relative error of $3.2\times10^{-3}$ with standard deviation $2.1\times10^{-3}$ at iteration 500, and a **runtime under 6.8 seconds**.

**(ii) Very-high-dimensional diffusion, $d=10^4$ (Section 4.2, Table 1).** A source $Q$ is added so that the exact solution is $v(t,x)=V((t-0.5)\mathbf 1_d+x)$ with the same pairwise-coupled $V$ used in paper 93's Section 4.2, deliberately non-separable. 9000 iterations.

| Setting | Region $S$   | Width $W$ | Mean RE            | SD of RE           | RT (s) |
| ------- | ------------ | --------- | ------------------ | ------------------ | ------ |
| 1       | $\{-0.5\}$   | 1024      | $5.5\times10^{-3}$ | $2.0\times10^{-3}$ | 295    |
| 2       | $[-1.5,1.5]$ | 1024      | $1.8\times10^{-2}$ | $3.1\times10^{-3}$ | 296    |
| 3       | $[-1.5,1.5]$ | 10240     | $5.4\times10^{-3}$ | $5.3\times10^{-4}$ | 5410   |

**The conclusion is clean: accuracy degrades when the solved region is enlarged at fixed width, and is recovered by widening the network at about $18\times$ the cost. In other words the bottleneck is network expressivity, not the martingale machinery.**

**(iii) Quasilinear PDEs at $d=10^4$ (Section 4.3).** The three cases QLP-1, QLP-2a and QLP-2b have the same form as those listed under paper 93's Section 4.2, including the dense-diffusion case.

**(iv) HJB equations (Section 4.4, Table 2).** $\mathcal L^\kappa=(b+2\kappa)^\top\partial_x+\delta^2\mathrm{Tr}\{\partial^2_{xx}\}$, $c(t,x,\kappa)=\delta^{-2}|\kappa|^2$, $U=\mathbb R^d$, $T=1$; exact solution $v(t,x)=-\ln\mathbb E[\exp(-g(X^{t,x}_T))]$ with $X^{t,x}_T=x+(T-t)b+\sqrt2\delta B_{T-t}$, reference by Monte Carlo with $10^6$ samples. Instances: **HJB-1** ($b=0$, $\delta=1$, $g=\ln(0.5(1+|x|^2))$); **HJB-2a** ($b=\mathbf 1_d$, $\delta=0.1$); **HJB-2b** ($\delta=0.05$); **HJB-3a** ($b=\mathbf 1_d$, $\delta=0.2$, $\bar g(x)=\frac1d\sum_i\{\sin(x_i-\frac\pi2)+\sin((0.1\pi+x_i^2)^{-1})\}$, highly oscillatory near 0); **HJB-3b** ($\delta=0.1$). At $d=10^4$, solving $v(0,s\mathbf 1_d)$ for $s\in[-1,1]$, 9000 iterations.

| Equation | Hidden layers $H$ | Mean RE            | SD of RE           | RT (s) |
| -------- | ----------------- | ------------------ | ------------------ | ------ |
| HJB-1    | 4                 | $2.2\times10^{-3}$ | $3.4\times10^{-4}$ | 9432   |
| HJB-2a   | 4                 | $7.5\times10^{-3}$ | $4.5\times10^{-4}$ | 9423   |
| HJB-2b   | 4                 | $2.1\times10^{-2}$ | $4.3\times10^{-4}$ | 9425   |
| HJB-3a   | 4                 | $2.4\times10^{-2}$ | $1.7\times10^{-3}$ | 9422   |
| HJB-3b   | 6                 | $2.3\times10^{-2}$ | $5.4\times10^{-4}$ | 13996  |

**Shrinking $\delta$ (HJB-2a to 2b) raises the error by nearly a factor of three**: the weaker the diffusion and the less smooth the solution, the harder the method works — the other face of the "the pilot must cover the region of interest" limitation. The $d=2000$ width study (HJB-3b, 6000 iterations) confirms the reading of Table 1: $W=d+10$ gives mean RE $6.9\times10^{-2}$, SD $3.6\times10^{-3}$ and 540 seconds, while $W=5d+10$ gives $2.0\times10^{-2}$, $1.2\times10^{-3}$ and 9050 seconds.

### Relation to the others

**It is the direct successor of paper 86.** That paper's conclusion asks for a derivative-free method avoiding the $d\times d$ Hessian, and paper 96 delivers exactly that, citing paper 86 as what it improves. The two share the adversarial/Galerkin enforcement of the martingale property and the ReLU6 control network, but paper 96 replaces the single long martingale with short-horizon system processes and eliminates all automatic differentiation from the loss.

**It is the direct predecessor of paper 93:** the random finite difference operator introduced in paper 96's v1 — a univariate first-order random difference that nonetheless approximates a second-order operator — is precisely the object developed in paper 93, and v1's observation that the resulting scheme is a transition-density-weighted weak adversarial network is the conceptual link. Paper 100 is the third member of the derivative-free family.

The contrast with the classical line (papers 8, 19, 25, 26, 41, 47, 63, 68) is the same sentence as before: those prove orders and stability for $d\lesssim10$, while paper 96 reaches $d=10^4$ with no theory. What is interesting is that **both lines share the same structural identity** — the residual of the parabolic PDE equals the drift of the process $\mathcal M_t=v(t,X_t)+\int f$, which is the nonlinear Feynman-Kac relation the FBSDE schemes discretise. The quasilinear setting here ($\mu,\sigma$ depending on $v$) is the deep-learning analogue of the **coupled** FBSDEs treated in papers 8 and 23.

## 97: DeepSPoC, handing the particles' memory to a network

### The idea

Nonlinear Fokker-Planck equations (porous medium, fractional porous medium, Keller-Segel, Curie-Weiss) are the forward equations of **mean-field SDEs** $\mathrm dX_t=b(t,X_{t-},\mu_t)\mathrm dt+\sigma(t,X_{t-},\mu_t)\mathrm dZ_t$ with $\mu_t=\mathrm{Law}(X_t)$ and $Z$ a Lévy process. They cannot be simulated by plain Monte Carlo because the mean-field term destroys the Markov property. Classical **propagation-of-chaos (PoC) particle methods** replace $\mu_t$ by the empirical measure of $N$ interacting particles, at two costs: all $N$ particles must be advanced simultaneously, which is expensive; and reconstructing the solution requires storing **all particle trajectories**, which is prohibitive in high dimension.

**The mechanism of this paper is to let the neural network serve as the particles' memory.** It builds on the **sequential propagation of chaos (SPoC)** theory of Du-Jiang-Li, in which each new particle interacts only with the _preceding_ ones so that the system takes an iterative form; the stored measure is then replaced wholesale by a neural network density, and trajectories are never stored.

**One conceptual oddity is worth settling first**: the "loss" here is not an objective being minimised. The particle batch, and hence the loss itself, changes every epoch; gradient descent is only the mechanism by which the network _absorbs_ the new batch, mimicking the weighted average $\mu^n_t=(1-\alpha_n)\mu^{n-1}_t+\alpha_n\hat\mu^n_t$. That is why exactly one descent step is taken per epoch, and why **the loss value does not decay during training** (Figure 2: the loss plateaus while the relative $L^2$ error keeps dropping). The loss therefore cannot serve as a stopping criterion, which is precisely what motivates the posterior error estimate of Section 3.2.

### Setting

The classical PoC particle system and the SPoC particle system are

$$
\mathrm dX^{n,N}_t=b(t,X^{n,N}_{t-},\mu^N_t)\mathrm dt+\sigma(t,X^{n,N}_{t-},\mu^N_t)\mathrm dZ^n_t,
\qquad \mu^N_t=\frac1N\sum_{i=1}^N\delta_{X^{i,N}_t},
$$

$$
\mathrm dX^n_t=b(t,X^n_{t-},\mu^{n-1}_t)\mathrm dt+\sigma(t,X^n_{t-},\mu^{n-1}_t)\mathrm dZ^n_t,
\qquad \mu^n_t=\mu^{n-1}_t+\alpha_n\bigl(\delta_{X^n_t}-\mu^{n-1}_t\bigr),
$$

with $\alpha_n\downarrow0$ the update rate. The batch version uses $\hat\mu^n_t=\frac1K\sum_{i=1}^K\delta_{X^{i,n}_t}$ and $\mu^n_t=\mu^{n-1}_t+\alpha_n(\hat\mu^n_t-\mu^{n-1}_t)$.

### Derivation

**Step one: the deepSPoC system.** Replace the stored measure $\mu^{n-1}_t$ by a neural network density $\rho_{{\rm NN},\theta_{n-1}}(t,\cdot)$:

$$
\begin{cases}
\mathrm dX^{i,n}_t=b\bigl(t,X^{i,n}_{t-},\rho_{{\rm NN},\theta_{n-1}}(t,\cdot)\bigr)\mathrm dt
+\sigma\bigl(t,X^{i,n}_{t-},\rho_{{\rm NN},\theta_{n-1}}(t,\cdot)\bigr)\mathrm dZ^{i,n}_t,
& i=1,\dots,K,\\[4pt]
\rho_{{\rm NN},\theta_n}=\mathcal F\bigl(\rho_{{\rm NN},\theta_{n-1}},\hat\mu^n,\alpha_n\bigr).
\end{cases}
$$

The operator $\mathcal F$ is: build a loss from the current network and the fresh batch's empirical measure, compute $\nabla_\theta$, and take **exactly one** gradient-descent step with learning rate $\alpha_n$. The learning rate plays the role of the SPoC update rate and is decayed as $\alpha_n=\alpha_0\gamma^{\lfloor n/\Gamma\rfloor}$ with $\gamma<1$.

**Step two: time discretisation and mollification.** Euler-Maruyama on $0=t_0<\cdots<t_M=T$ with $\Delta t=T/M$ and $X^i_0\overset{\rm iid}\sim\mu_0$ gives $\hat\mu_{t_m}=\frac1K\sum_i\delta_{X^i_{t_m}}$. The mollified empirical density is $\hat\rho_{t_m}:=\hat\mu_{t_m}*f_\epsilon$ with the Gaussian mollifier $f_\epsilon(x)=(2\pi\epsilon^2)^{-d/2}e^{-|x|^2/(2\epsilon^2)}$.

**Step three: three loss functions.**

1. **$L^2$-distance loss $\mathcal L_{\rm sq}$**, used with a plain fully connected network $\rho_{FC,\theta}$:
   $$
   \mathcal L_{\rm sq}=\sum_{m=0}^M\bigl\|\rho_{FC,\theta}(t_m,\cdot)-\hat\rho_{t_m}(\cdot)\bigr\|^2_{L^2}
   \approx\frac1N\sum_{m=0}^M\sum_{x\in S}\Bigl|\rho_{FC,\theta}(t_m,x)
   -\frac1K\sum_{i=1}^K\frac{e^{-|x-X^i_{t_m}|^2/(2\epsilon^2)}}{(2\pi\epsilon^2)^{d/2}}\Bigr|^2,
   $$
   with $S=\{x_j\}_{j=1}^N$ sampled uniformly from a truncation box $\Omega_0$. **Remark 2**: in back-propagation $\hat\rho_{t_m}(x_j)$ is _detached_, treated as a $\theta$-independent constant even though the particles were generated using the network.
2. **KL loss $\mathcal L_{\rm kl}$**, used with a **temporal normalizing flow** (KRnet), which is automatically nonnegative with unit mass and can be sampled directly:
   $$
   \mathcal L_{\rm kl}=-\sum_{m=0}^M\int_{\mathbb R^d}\hat\rho_{t_m}(x)\log\rho_{NF,\theta}(t_m,x)\,\mathrm dx
   \approx-\sum_{m=0}^M\frac1N\sum_{x\in S}\hat\rho_{t_m}(x)\log\rho_{NF,\theta}(t_m,x).
   $$
3. **Path / negative-log-likelihood loss $\mathcal L_{\rm path}$**, the only mollification-free option and hence the one used in the highest dimensions:
   $$
   \mathcal L_{\rm path}=-\sum_{m=0}^M\frac1K\sum_{i=1}^K\log\rho_{NF,\theta}\bigl(t_m,X^i_{t_m}\bigr).
   $$

**Step four: adaptive spatial sampling.** Uniform sampling of $\Omega_0$ fails in high $d$ because the support of the density occupies a vanishing fraction of the box. The fix reuses the algorithm's own particles: $S_m=S_{\rm uniform}\cup S^m_{\rm adaptive}$ with $S^m_{\rm adaptive}=\{x_j=\hat X^j_{t_m}+\sigma\delta_j\}$ and $\delta_j\overset{\rm iid}\sim N(0,1)$, where the $\hat X^j_{t_m}$ come from the previous epoch and the Gaussian jitter $\sigma$ preserves exploration. For the KL loss the adaptive strategy differs (importance sampling of the integrals); for $\mathcal L_{\rm path}$ no training set is needed at all, so no adaptivity is designed. The three algorithms are: fully connected network with $\mathcal L_{\rm sq}$ (Algorithm 1), normalizing flow with $\mathcal L_{\rm kl}$ (Algorithm 2), and normalizing flow with $\mathcal L_{\rm path}$ (Algorithm 3).

### Theorems

The analysis is deliberately set in a **simplified surrogate**: the density is represented by a truncated Fourier basis with $L^2$-projection $P_N$, so that one gradient step on $\mathcal L_{\rm sq}$ is _exactly_ an affine update.

- **Exact-update identity.** In the Fourier surrogate $\mathcal F$ acts as $\rho^{\theta_t}_n=(1-2\alpha_n)\rho^{\theta_t}_{n-1}+2\alpha_nP_N(\hat\rho_t)$, so unrolling gives $\rho^{\theta_t}_{n-1}=P_N(\sum_{l=0}^{n-1}\frac{\beta_l}K\sum_i\delta^\epsilon_{\tilde X^{i,l}_t})$ with $\beta_l=(1-2\alpha_{n-1})\cdots(1-2\alpha_{l+1})2\alpha_l$ and $\sum_l\beta_l=1$. **This is the precise sense in which deepSPoC "is" SPoC: the network update is a weighted average of all past batches.**
- **Basis-size condition.** $\hat\rho_t$ is Lipschitz with constant $\tilde C/(K\epsilon^{d+1})$, so $\|\hat\rho_t-P_N\hat\rho_t\|_\infty\le C(\epsilon,K,L,d)(\log N)^d/N$; choosing $N$ so that $\|\rho-P_N\rho\|_\infty\le\epsilon/(2L)^{d+1}$ gives $W_1(\rho,P_N\rho)\le\epsilon$ and $|\|\rho\|_1-\|P_N\rho\|_1|\le\epsilon$. **The choice of $W_1$ rather than a general $W_r$ is deliberate**: control by $\|\cdot\|_\infty$ is homogeneous only for $r=1$.
- **Assumption 3.1.** $|b(t,x,\mu)-b(t,y,\nu)|+\|\sigma(t,x,\mu)-\sigma(t,x,\nu)\|\le C(|x-y|+W_1(\mu,\nu))$ and $|b(t,0,\mu)|+\|\sigma(t,0,\mu)\|\le C(1+\|\mu\|_1)$.
- **Proposition 3.2 (uniform particle moments, justifying the truncation).** With $\mu_0\in\mathcal P_3$ and $N$ as above, there is $C_0$ independent of $L$ with $\mathbb E[\sup_{0\le s\le T}|X^{i,n}_s|^3]\le C_0$ and $\mathbb E[\mathbf 1_{\{\sup_s|X^{i,n}_s|^2\ge L_0\}}\sup_s|X^{i,n}_s|^2]\le C_0/L_0$.
- **Theorem 3.3 (convergence).** Under Assumption 3.1, $\mu_0\in\mathcal P_3$, $N$ satisfying the basis-size condition, **equal batch weights $\beta_l=1/n$**, and a truncation radius $L$ chosen so that $\mathbb E[\sup_s|X^{i,n}_s-\tilde X^{i,n}_s|^2]\le\epsilon^2$, there exists $C$ independent of $n,K,L,\epsilon$ with
  $$
  \mathbb E\,W_1^2\bigl(\rho^{\theta_t}_{n-1},\mu_t\bigr)\ \le\ C\Bigl(\epsilon^2+(Kn)^{-\frac{1}{1+d/2}}\Bigr).
  $$
  The second term is the usual PoC/empirical-measure rate in $Kn$, the total number of particles simulated; the first is the mollification and truncation floor. **The caveat the authors state themselves: this holds in the Fourier-projection surrogate, not for the actual neural network.**
- **Posterior error estimate (Section 3.2).** Restricted to Brownian-driven dynamics, under **Assumption 3.2** (the same Lipschitz conditions with $W_2$ in place of $W_1$, noted to be _weaker_ than 3.1 since $W_1\le W_2$), define on $\mathcal P_{2,\infty}([0,T])$ the metric $H_\alpha(\mu_\cdot,\nu_\cdot):=(\int_0^Te^{-\alpha t}W_2^2(\mu_t,\nu_t)\mathrm dt)^{1/2}$ and the solution map $\Phi(\mu_\cdot):=t\mapsto\mathrm{Law}(X_t)$ for the _frozen_, hence Markovian and directly simulable, SDE.
  - **Proposition 3.5:** with $C_0=2(T+1)C^2$ and $\alpha>C_0$, $H_\alpha(\Phi(\mu_\cdot),\Phi(\nu_\cdot))\le\sqrt{C_0/(\alpha-C_0)}\,H_\alpha(\mu_\cdot,\nu_\cdot)$, so $\Phi$ is a contraction.
  - **Theorem 3.6:** with $\mu^*_\cdot$ the fixed point (the true solution), for $\alpha>2C_0=4(T+1)C^2$ and **any** $\mu_\cdot\in\mathcal P_{2,\infty}([0,T])$,
    $$
    H_\alpha(\mu_\cdot,\mu^*_\cdot)\ \le\ \Bigl(1-\sqrt{\tfrac{C_0}{\alpha-C_0}}\Bigr)^{-1}H_\alpha\bigl(\mu_\cdot,\Phi(\mu_\cdot)\bigr).
    $$
  Computationally: feed the learned density into the now-Markovian SDE, simulate, and compare the resulting law with what you fed in; the discrepancy is a **computable upper bound** on the true error. **This is the only computable error bound among the six papers on this page** — paper 93 gives an a priori rate instead. The mechanism, freeze the nonlinearity and measure the mismatch, is a fixed-point argument rather than a truncation-error argument.
- Absent: no analysis for the actual neural-network parameterisation, for the Lévy-driven or singular-kernel cases, or for the KL and path losses; the authors explicitly note that several of their test equations "go beyond the scope of the SPoC theory or even PoC theory".

### Numerical experiments

Networks: fully connected with 6 hidden layers of 512 neurons and ReLU; temporal KRnet with 10 affine coupling blocks, each with two 512-neuron hidden layers; Adam; PyTorch. **The test suite is chosen for degeneracy and nonlocality rather than for order verification — no convergence orders in $\Delta t$ are reported anywhere.**

**(i) Porous medium equation (Section 4.1).** $\partial_t\rho=\Delta\rho^m$, quasilinear and degenerate at $\rho=0$, benchmarked against the **Barenblatt solution** $U_{m,C}(t,x)=t^{-\alpha}\{C-\frac{m-1}{2m}\beta\frac{|x|^2}{t^{2\beta}}\}_+^{1/(m-1)}$ with $\alpha=\frac d{d(m-1)+2}$ and $\beta=\alpha/d$ — a weak solution that loses classical differentiability at the free boundary. The runs are:

| Dimension | Algorithm                 | Parameters                                                                                                                   |
| --------- | ------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| 1D        | Alg. 1 (uniform)          | $t_0=1$, $T=1$, $\Delta t=0.01$, $N=1000$                                                                                    |
| 3D        | Alg. 1 (adaptive)         | $t_0=0.1$, $T=0.2$, $\Delta t=0.005$, noise intensity $0.2$, $\Omega_0=[-2,2]^3$                                             |
| 5D        | Alg. 1 (adaptive), Alg. 3 | 2000 uniform + 4000 adaptive points, $\sigma=0.3$, $\epsilon=0.05$, $K=8000$, $\Omega_0=[-3,3]^5$                            |
| 6D, 8D    | Alg. 3 only               | $t_0=1$, $T=1.5$, $\delta t=0.025$, $10^4$ samples, $\Omega_0=[-3,3]^6$ / $[-3,3]^8$, $\gamma=0.5$ every $\Gamma=500$ epochs |

**One explicit finding: Algorithm 1 is hard to use at 6D and 8D because high-dimensional mollification is inaccurate; it is the mollification-free $\mathcal L_{\rm path}$ (Algorithm 3) that makes those dimensions reachable.**

**(ii) Deterministic-particle variant (Section 4.1.5).** The same PDE is attacked through the ODE $\mathrm dX_t=\nu m\rho^{m-2}(X_t)\nabla\rho(X_t)\mathrm dt$, the blob-method / Wasserstein-gradient-flow form of Carrillo-Craig-Patacchini. Because this needs $\nabla\rho$, the activation is switched from ReLU to **SoftPlus** $\frac1\beta\log(1+e^{\beta x})$ with $\beta=20$, and PyTorch AD supplies the derivative. **This demonstrates that deepSPoC is agnostic to whether the underlying particle method is stochastic or deterministic.**

**(iii) Keller-Segel (Section 4.2).** $\partial_t\mu=\nabla\cdot((\nabla W*\mu)\mu)+\Delta\mu$ in $d=2$ with the **singular** kernel $W(x)=\frac1{2\pi}\ln|x|$ — a distribution-dependent, singular-interaction case outside the PoC theory. Solved with both Algorithm 1 and Algorithm 3.

**(iv) Curie-Weiss mean-field equation (Section 4.3).** $\mathrm dX_t=\{-\beta(X_t^3-X_t)+\beta K\mathbb EX_t\}\mathrm dt+\mathrm dB_t$ with $\beta=1$ and $K=-0.1$, an _expectation_-dependent mean-field SDE. Long-time behaviour is compared with the explicit invariant density $p^*(x)=C^{-1}\exp\{-2\beta(\frac{x^4}4-\frac{x^2}2)\}$ and with a PoC reference computed using **5 million particles**; deepSPoC uses 5000 epochs of 1000 particles each — **the same total particle count, but never held in memory at once, which is the memory argument made concrete.**

**(v) Fractional porous medium equation (Section 4.4).** $\partial_t\rho=-(-\Delta)^{\alpha/2}(|\rho|^{m-1}\rho)$ with $m>1$ and $\alpha\in(0,2)$, the fractional Laplacian written as a principal-value integral — a **nonlocal** operator driven by an $\alpha$-stable Lévy process rather than Brownian motion.

### Relation to the others

**It is the methodological odd one out on this page**, and the only paper in the list about **forward** (Fokker-Planck / McKean-Vlasov) problems rather than backward (FBSDE / HJB) ones. Where papers 8, 47, 63, 86, 93 and 96 propagate information _backward_ from a terminal condition using conditional expectations, deepSPoC propagates a _law_ forward and fits it with a network.

The engineering DNA is nevertheless shared: the network replaces stored trajectories; training is driven by freshly simulated Euler-Maruyama paths; adaptive, dynamics-based spatial sampling is essential in high $d$. In particular the "sample where the process actually goes" principle that paper 93 derives from the Fokker-Planck adjoint equation is here literally the object being learned.

**The mean-field link to papers 33 and 61** (see the [[en/computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control|second-order FBSDEs page]] and the [[en/computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|multistep schemes page]]): those discretise **mean-field BSDEs and FBSDEs**, backward equations whose generator depends on the law of the solution, with high-order $\theta$- and multistep schemes at low $d$; paper 97 treats the forward McKean-Vlasov side of the same circle of problems with no order guarantees but at $d$ up to 8. The two are complementary halves of "numerics for distribution-dependent SDEs". The connection to paper 108 is that the McKean-Vlasov forward dynamics here are exactly the population-evolution half of a mean-field game.

## 100: replacing the network with local linear regression

### The idea

For semilinear parabolic PDEs in $d\gg1$, the two dominant families both have gaps. **Deep learning solvers** (deep BSDE, SOC-MartNet, DRDM) reach very high $d$, but the paper's diagnosis is blunt: "limited stability of the optimization procedure, pronounced sensitivity to hyperparameters, and a lack of rigorous a priori error estimates." **Classical probabilistic and regression Monte Carlo BSDE solvers** (Gobet-Lemor-Warin and descendants) are analysable, but they approximate conditional expectations by regression onto _global_ basis functions whose number grows combinatorially with $d$ and which become ill-conditioned beyond $d\approx200$, and they solve the coupled $(Y,Z)$ system by Picard iteration, whose cost explodes in high dimension. Sparse grids stall around $d\approx10$; branching-diffusion representations suffer variance blow-up over long horizons and cost $O(d^2)$.

**The mechanism here is to keep the martingale time discretisation and the derivative-free spirit but replace the neural network entirely by local linear regression on a small particle ensemble.** Two judgements carry the design: the conditional expectation is approximated by an **ensemble average over all particles** (and the particle number $M$ enters the error only through $e^{-cM}$, so $M\le100$ suffices), and $Z$ comes from a **weighted least-squares fit of $\nabla u$**, which reorders the pipeline to $X\to Z\to Y$ and removes the Picard iteration. What is bought back is interpretability, a rigorous a priori bound, and laptop-scale cost.

### Setting

Terminal-value semilinear parabolic PDE

$$
(\partial_t+\mathcal L)u(t,x)+f\bigl(t,x,u(t,x),\sigma^\top\nabla u(t,x)\bigr)=0,\qquad u(T,x)=g(x),
\qquad
\mathcal Lu=\tfrac12\mathrm{Tr}\bigl[\sigma\sigma^\top\mathrm{Hess}_xu\bigr]+\langle\mu,\nabla u\rangle,
$$

with the **decoupled** FBSDE representation $Y_t=u(t,X_t)$, $Z_t=\sigma^\top(t,X_t)\nabla u(t,X_t)$. **Assumption 2.1:** global Lipschitz continuity in all arguments with constant $L$, plus linear growth $\|\mu\|+\|\sigma\|+|f|+|g|\le C(1+\|x\|+|y|+\|z\|)$. **Lemma 2.1:** with the additional uniform nondegeneracy $\xi^\top\sigma\sigma^\top\xi\ge\lambda\|\xi\|^2$, the FBSDE has a unique adapted solution in $S^2(\mathbb R^d)\times S^2(\mathbb R)\times H^2(\mathbb R^d)$.

### Derivation

**Step one: the martingale / conditional-expectation time discretisation.** Integrating the backward equation over $[t_k,t_{k+1}]$ and taking $\mathbb E_k[\cdot]$ kills the Itô integral, giving $Y_k=\mathbb E_k[Y_{k+1}+\int_{t_k}^{t_{k+1}}f\,\mathrm ds]$; freezing the coefficients at $t_k$ (the left-rectangle, that is Euler, rule) gives the semidiscrete scheme, **implicit in $Y$**:

$$
\tilde Y_k=\mathbb E_k\bigl[\tilde Y_{k+1}+f(t_k,X_k,\tilde Y_k,\tilde Z_k)\Delta t\bigr],\qquad 0\le k\le N .
$$

**This is the $\theta=1$ member of the family**, the same $Y$-part that paper 93's Remark 5 writes down, and the ancestor of the schemes in papers 8, 47 and 63.

**Step two: a small-scale stochastic particle method.** Simulate $M$ _independent_ Euler-Maruyama trajectories and approximate the conditional expectation by the **ensemble average over all $M$ particles**:

$$
\tilde Y^m_k\approx\frac1M\sum_{j=1}^M\tilde Y^j_{k+1}+f\bigl(t_k,\tilde X^m_k,\tilde Y^m_k,\tilde Z^m_k\bigr)\Delta t,
\qquad 1\le m\le M,
\qquad \tilde Y_0=\frac1M\sum_m\tilde Y^m_0 .
$$

**Crucially $M$ is small** — typically $M\le100$ — because by Theorem 3.1 it enters the error only through an exponentially suppressed bad-event probability $e^{-cM}$, not through a $1/\sqrt M$ Monte Carlo term.

**Step three: $Z$ by weighted local linear regression, the distinctive piece.** Instead of Picard-iterating the coupled $(d+1)$-dimensional system, the authors **decouple and compute $Z$ first**. Locally expand $u(t_{k+1},\cdot)$ around the anchor $\tilde X^m_k$,

$$
u(t_{k+1},\cdot)\approx
\underbrace{u(t_k,\tilde X^m_k)+\partial_tu(t_k,\tilde X^m_k)\Delta t}_{=:\alpha}
+\underbrace{\nabla u(t_k,\tilde X^m_k)}_{=:\alpha_x}{}^\top(\cdot-\tilde X^m_k),
$$

and fit $(\alpha,\alpha_x)\in\mathbb R^{d+1}$ by weighted least squares against the _already computed_ values $\tilde Y^j_{k+1}$ regressed on the _current_ positions $\tilde X^j_k$:

$$
J(\boldsymbol\alpha)=\sum_{j=1}^Mw_j\bigl(\tilde Y^j_{k+1}-\alpha-\alpha_x^\top D_j\bigr)^2,
\qquad D_j:=\tilde X^j_k-\tilde X^m_k,
\qquad w_j=\frac{K(\|D_j\|/\varepsilon_k)}{\sum_{i=1}^MK(\|D_i\|/\varepsilon_k)},
$$

with $K$ a kernel (Gaussian in practice). Then $\tilde Z^m_k=\sigma^\top(t_k,\tilde X^m_k)\alpha_x$ and **the intercept $\alpha$ is discarded**. Three implementation details are worth recording:

- **Matrix-free solve.** LSQR or PCG using only the two products $\beta_j=w_j(\alpha+D_j^\top\alpha_x)$ and $(\mathbf D^\top\beta)_0=\sum_j\beta_j$, $(\mathbf D^\top\beta)_{1:d}=\sum_j\beta_jD_j$. Storage is $O(d)$, cost $O(Md)$ per time step, and $\mathbf D^\top W\mathbf D$ is never formed.
- **Ridge regularisation (Remark 2.2).** Since $M\ll d$ is the normal regime, the normal equations are underdetermined, so the implementation minimises $J_\lambda(\boldsymbol\alpha)=\sum_jw_j(\cdots)^2+\lambda\|\boldsymbol\alpha\|^2$.
- **Why kernel weights instead of $k$-nearest-neighbours (Remark 2.3).** Under distance concentration in high $d$, $k$-NN inflates the radii; but for the Gaussian kernel the _ratio_ $w_j/w_i=\exp\{-(\|D_j\|^2-\|D_i\|^2)/\varepsilon_k^2\}$ still discriminates on the relative gap $|\|D_j\|-\|D_i\||$ even when absolute distances concentrate. Using _all_ particles also removes the radius-tuning hyperparameter of conventional local linear regression.

**Step four: $Y$ by a scalar Newton iteration.** With $\tilde X^m_k$ and $\tilde Z^m_k$ fixed, the scheme is a **one-dimensional** root problem $F(\tilde Y^m_k)=0$, solved by Newton's method and reported to need only 2 to 3 iterations per step.

**Algorithm 2.1** runs: set the terminal condition $Y^j_N=g(X^j_N)$ in parallel; make a **single forward pass** generating all $\tilde X^j_k$ (paths are stored, not recomputed); then loop backward over $k=N-1,\dots,0$ doing, per particle, the matrix-free regression for $\alpha_x$, the assignment $\tilde Z_k\leftarrow\sigma^\top\nabla u$, and Newton for $\tilde Y^m_k$, **fully parallel over $m$**; finally take $\tilde Y_0=\frac1M\sum_m\tilde Y^m_0$. **The whole pipeline is the ordering $X\to Z\to Y$, in contrast to the coupled $(Y,Z)$ solves of regression Monte Carlo.**

### Theorems

- **Lemma 3.1.** The standard strong order $1/2$ for Euler-Maruyama: $\max_{0\le t\le T}\mathbb E\|X_t-\tilde X_t\|^2\le C\Delta t$ and $\mathbb E[\sup_t\|X_t-\tilde X_t\|^2]\le C\Delta t$.
- **Lemma 3.3 (local truncation error of the backward scheme).** If $f\in C^{1,2}$ and satisfies the Lipschitz condition, then
  $$
  |\mathcal E_k|:=\Bigl|\mathbb E_k\Bigl[\int_{t_k}^{t_{k+1}}f(s,X_s,Y_s,Z_s)\mathrm ds\Bigr]-f(t_k,X_k,Y_k,Z_k)\Delta t\Bigr|\le C(\Delta t)^2 .
  $$
- **Lemma 3.8 (the technical heart: error of the local-linear-regression gradient estimator).** Let $\alpha_x$ be the finite-sample minimiser of the weighted least-squares problem with perturbed data $\tilde Y^j_{k+1}=Y^j_{k+1}+\delta^j_{k+1}$, the errors being i.i.d. conditional on $\mathcal F_{t_k}$. Then
  $$
  \mathbb E_k\bigl[\|\alpha_x-\nabla u(t_k,x)\|^2\bigr]\ \le\ C\varepsilon_k^2
  \;+\;C\varepsilon_k^{-2}\,\mathbb E_k\bigl[|\delta_{k+1}|^2\bigr]
  \;+\;C\,e^{-C_{A_k}M\varepsilon_k^d}.
  $$
  The three terms are exactly **bias** (the first-order Taylor truncation, $O(\varepsilon_k^2)$), **variance amplification** (noise divided by the bandwidth squared, the classic bias-variance trade-off), and the probability of the "insufficient sampling in the $\varepsilon_k$-ball" event, which decays like $e^{-CM\varepsilon_k^d}$.
- **Theorem 3.1 (global error).** With $\delta_k:=\tilde Y_k-Y_k$, under Assumption 2.1 and the hypotheses of Lemma 3.3, for all sufficiently large $M$,
  $$
  \mathbb E\bigl[|\delta_0|^2\bigr]\ \le\ C\,\Delta t\;+\;C\,\Delta t\,e^{-c_1M},
  $$
  where $C>0$ depends only on $T$ and the Lipschitz constant $L$ and is independent of $\Delta t$, and $c_1$ depends on $\varepsilon_k$. The proof is the standard structure: subtract exact from discrete, split into the Monte Carlo deviation $\frac1M\sum_j\tilde Y^j_{k+1}-\mathbb E_k[\tilde Y_{k+1}]$, the propagated error $\mathbb E_k[\delta_{k+1}]$, the Lipschitz term $|f(t_k,X_k,\tilde Y_k,\tilde Z_k)-f(t_k,X_k,Y_k,Z_k)|\le L(|\delta_k|+\|\sigma\|\|\alpha_{x,k}-\nabla u_k\|)$ and the truncation error $\mathcal E_k$; then conditional Jensen, Young's inequality with a parameter $\eta$, and discrete Grönwall.
- **What the bound says.** The estimate is "fully explicit in both the particle number $M$ and the time step size $\Delta t$"; because $M$ appears only inside $e^{-c_1M}$ rather than as $M^{-1/2}$, **a moderate $M$ (empirically about 100) already suffices**; and the rate is $O(\Delta t)$, **first order in time**, matching what paper 93 proves for the martingale/RDM family and beating the $O(\Delta t^{1/2})$ of deep-BSDE analyses.
- Absent: the constants' dependence on $d$ is asserted to be benign but is not quantified as an explicit dimension-free statement; the interplay between $\varepsilon_k$, $M$ and $d$ (the bound needs $M\varepsilon_k^d$ large, which is where any residual curse could hide) is not optimised; there is no analysis of the ridge parameter $\lambda$ actually used in the code, nor of the Newton iteration.

### Numerical experiments

All experiments run on a **MacBook Pro (Apple M1 Pro, 10 cores, 32 GB)** — the paper makes a point of requiring no specialised hardware.

**(i) Allen-Cahn (Section 4.1).** $\partial_tu+\Delta u+f(u)=0$.

- **Case 1**, double-well $f(u)=u-u^3$, $u(T,x)=1/(2+0.4\|x\|^2)$, $T=0.3$, $d=100$, evaluated at $x_0=0$; the reference $u(0,x_0)\approx0.0528$ comes from the branching-diffusion method (the same benchmark used in the deep BSDE literature and by paper 96). **The observed order: log-log slopes close to 1 in both absolute and relative error, that is first order in $\Delta t$**, matching Theorem 3.1; varying $M$ changes the accuracy level but not the rate. At $N=10^4$ ($\Delta t=3\times10^{-5}$) with $M=100$ the absolute error is about $1.2\times10^{-5}$.
- **Case 2**, logarithmic potential $f(u)=\frac\theta2\ln\frac{1+u}{1-u}-\theta_cu$ with $\theta<\theta_c$, with the manufactured solution $u(t,x)=\cos(\prod_{j=1}^dx_j)e^{\cos t-\|x\|^2}$ and a compensating source; $T=1$, evaluated at $x=0$ and $x=(0.1,\dots,0.1)$, at $d=100$ and $d=1000$. First-order convergence is observed at both dimensions; at $d=100$ the relative error reaches about $10^{-3}$ once $\Delta t<6.25\times10^{-5}$ ($N\ge16000$) and decreases linearly thereafter. **Notably this $f$ is not globally Lipschitz, yet Newton still converges quickly.** At large $N$ and $M$ the dominant error becomes the **local regression bias** of Lemma 3.8 rather than time stepping or Monte Carlo noise. Runtimes (seconds, with $\Delta t$ halved four times):

| Dimension | Initial $\Delta t$ | $M$ | Runtimes over five refinements           |
| --------- | ------------------ | --- | ---------------------------------------- |
| $d=100$   | 0.002              | 50  | 1.24, 2.28, 4.57, 8.83, 17.27            |
| $d=100$   | 0.002              | 100 | 3.81, 7.45, 14.99, 29.55, 59.74          |
| $d=1000$  | 0.0005             | 50  | 36.52, 72.32, 145.91, 288.43, 583.69     |
| $d=1000$  | 0.0005             | 100 | 129.24, 257.59, 521.38, 1039.55, 2082.74 |

Runtime grows essentially linearly in $N\cdot M$.

**(ii) Burgers equation, $d=10^4$ (Section 4.2).** $\partial_tu+(u-\frac{2+d}2)\sum_i\partial_{x_i}u+\frac{d^2}2\nu\Delta u=0$ with terminal condition $u(T,x)=\frac{\exp(T+\sum_ix_i/d)}{1+\exp(T+\sum_ix_i/d)}$, so that $u(0,x_0)=0.5$ at $x_0=0\in\mathbb R^{10000}$, $T=0.3$. **Near first-order convergence in time.** No artificial viscosity is needed and the local-regression surrogate resolves the sharp-gradient structure directly, which the authors argue is where explicitly estimating $\nabla u$ pays off relative to deep solvers. Runtimes (seconds, $\Delta t=0.003$ halved four times): $M=100$ gives 1040.15, 2160.11, 4757.37, 10593.83, 22174.49; $M=200$ gives 2189.28, 4633.74, 9674.18, 21194.52, 45724.85.

**(iii) Hamilton-Jacobi type with a gradient-dependent sink, $d=500$ and $d=2000$ (Section 4.3).** $\partial_tu+u+f(t,x,u,\nabla u)=0$ with $R(u,\nabla u)=\kappa u\|\nabla u\|^2$, $\kappa=0.1$, and $f=\frac{4d}{(1+4t)^{(d+2)/2}}\frac{e^{-\|x\|^2}}{1+4t}-R(u,\nabla u)$, exact solution $u(t,x)=(1+4t)^{-d/2}\exp\{-\frac{\|x\|^2}{1+4t}\}$; $T=0.5$, $N$ up to $3\times10^4$, evaluated at $x=0$. **First-order convergence** at both dimensions, with Newton needing 2 to 3 iterations per step. One reported observation: **using roughly 10% of the global polynomial basis points via local linear regression gives the same accuracy**, whereas a global polynomial fit at $d=2000$ would be hopelessly overfitted or ill-conditioned. Runtimes (seconds, $d=500$, $\Delta t=0.0005$ halved four times): $M=50$ gives 89.07, 183.69, 363.86, 740.39, 1512.62; $M=100$ gives 185.06, 376.91, 765.88, 1517.73, 3016.35.

**Comparative claims made** (recorded as claims): comparable accuracy to deep BSDE with roughly **40% fewer total samples** on Allen-Cahn Case 2; cost linear in $d$ versus $O(d^2)$ for branching diffusion; and no heavy training.

### Relation to the others

**It is the "classical numerics strikes back" paper of the derivative-free trio.** Papers 96 and 93 made the martingale and derivative-free idea work at $d=10^4$ to $10^5$ using neural networks; paper 100 keeps the FBSDE/martingale time discretisation and the derivative-free spirit but replaces the network entirely by **local linear regression on a small particle ensemble**, recovering interpretability, a rigorous a priori bound and laptop-scale cost. Its introduction cites papers 86, 93 and 96 by name and states exactly this motivation.

**It shares the same $Y$-part as the classical schemes, differing only in how the conditional expectation is evaluated.** The $\theta=1$ recursion is precisely what paper 8's multistep scheme, paper 23's deferred corrections, paper 47's unified scheme, paper 63's $\theta$-schemes and paper 68's strong-stability-preserving schemes generalise. Where those papers evaluate $\mathbb E_k$ by **Gauss-Hermite quadrature on a tensor or sparse grid** (and therefore stop around $d\le6$), paper 100 evaluates it by a **particle ensemble average**; and where those papers get $Z$ from a second quadrature identity $Z\approx\frac1{\Delta t}\mathbb E_k[Y_{k+1}\Delta W^\top]$, paper 100 gets it from a **weighted least-squares fit of $\nabla u$**. **That substitution is the whole difference between order six at $d\le10$ and order one at $d=10^4$.**

The contrast with paper 25 is worth recording too: paper 25 attacks the same conditional-expectation bottleneck with **spectral sparse grids**, the deterministic answer, reaching moderate $d$, while paper 100 is the stochastic answer. Against Gobet-Lemor-Warin regression Monte Carlo the contrast is global bases plus Picard iteration on the coupled $(Y,Z)$ system, versus local bases plus the decoupled $X\to Z\to Y$. It agrees with paper 93 on first order in time by a different argument — one via the weak order of Euler-Maruyama, the other via local truncation plus Grönwall — and both explicitly contrast it with the $O(\Delta t^{1/2})$ of path-based deep BSDE.

## 108: deep policy iteration for mean-field games

### The idea

A finite-horizon mean-field game is characterised by a **coupled forward-backward PDE system**: an HJB equation running backward for the representative agent's value function, coupled to a Fokker-Planck equation running forward for the population density. In high dimension this is doubly hard. First, the coupling forces repeated evaluation of high-dimensional integrals against the state distribution. Second — the paper's sharper observation — **a single spatial sampling strategy must simultaneously** approximate those distributional integrals, enforce the HJB equation and resolve the Fokker-Planck equation; designing one sampler that does all three well is "highly nontrivial", and the difficulty persists in existing machine-learning approaches. On top of this, trajectory-based (fictitious-play style) methods must simulate full paths over $[0,T)$ at every outer iteration, which is inherently sequential in time, while HJB-based methods need pointwise Hamiltonian minimisation and second-order spatial derivatives.

**The mechanism is a regenerative reformulation: recast the finite-horizon game on $[0,T)$ as an infinite-horizon regenerative process on $[0,\infty)$ built from repeated cycles of length $T$**, with a resetting mechanism at the end of each cycle. Since the law is identical over every cycle, the finite-horizon cost and equilibrium condition become a **long-run average cost** and an **invariant occupation measure**. The payoff is that the cycle index can be identified with the **policy-iteration index**, so policy evaluation, policy improvement and measure estimation all happen cycle by cycle on a single short horizon rather than on the coupled HJB-Fokker-Planck system.

### Setting

The state is $x=(t,z)\in[0,T)\times\mathbb R^d$ with **the time variable carried as the first component of the state**, so

$$
b(x,\mu,u)=\bigl(1,\,b_z^\top(x,\mu,u)\bigr)^\top,
\qquad
\sigma(x,\mu,u)=\bigl(0_q,\,\sigma_z^\top(x,\mu,u)\bigr)^\top,
$$

$$
J(u,\mu)=\mathbb E\Bigl[\int_0^Tf\bigl(X^{\mu,u}_s,\mu,u(X^{\mu,u}_s)\bigr)\mathrm ds
+g\bigl(X^{\mu,u}_{T-},\mu\bigr)\Bigr].
$$

An **equilibrium** is a pair $(u^*,\mu^*)$ with $J(u^*,\mu^*)=\inf_uJ(u,\mu^*)$ and $\mu^*(A)=\frac1T\int_0^T\mathbb E[\mathbf 1(X^{\mu^*,u^*}_t\in A)]\mathrm dt$ for all $A\in\mathcal B([0,T)\times\mathbb R^d)$ — so $\mu^*$ is the **occupation measure on space-time**, not a family of time marginals. **Remark 1** notes this is slightly more general than the standard formulation, in which $b,\sigma,f,g$ depend on $\mu$ only through the time-$t$ marginal $\mu_t$; the space-time formulation is what makes the regenerative reformulation work.

### Derivation

**Step one: the regenerative reformulation.** Within cycle $i$, $t\in[iT,(i+1)T)$, the state obeys the same mean-field SDE started at $X^{\mu,u}_{iT}$; at the end of each cycle a **resetting mechanism** fires, $X^{\mu,u}_{(i+1)T}=\xi_{i+1}$ with $\xi_i\overset{\rm iid}\sim\mu_{X_0}$. Since the law of $X^{\mu,u}$ is identical over every cycle,

$$
J(u,\mu)=\lim_{t\to\infty}\frac Tt\mathbb E\Bigl[\int_0^tf\,\mathrm ds
+\sum_{i\ge1}g\bigl(X^{\mu,u}_{iT-},\mu\bigr)\mathbf 1(iT\le t)\Bigr],
\qquad
\mu^*(A)=\lim_{t\to\infty}\frac1t\int_0^t\mathbb E\bigl[\mathbf 1(X^{\mu^*,u^*}_s\in A)\bigr]\mathrm ds .
$$

**Step two: policy iteration through a dynamic-programming residual.** With $\mu^*$ frozen, the value function is $v(x)=\inf_u\mathbb E[\int_t^Tf\,\mathrm ds+g(X_{T-},\mu^*)\mid X_t=x]$ and the iteration is

$$
\text{(PE)}\ \ \text{find }v_{i+1}\in\mathcal V^{\mu^*}\ \text{with}\ \mathcal R_h(x;u_i,v_{i+1},\mu^*)=0,
\qquad
\text{(PI)}\ \ \text{find }u_{i+1}\in\arg\min_u\mathcal R_h(x;u,v_{i+1},\mu^*),
$$

$$
\mathcal R_h(x;u,v,\mu):=\mathbb E\Bigl[\int_t^{t+h}f\bigl(X^{\mu,u}_s,\mu,u(X^{\mu,u}_s)\bigr)\mathrm ds
+v\bigl(X^{\mu,u}_{(t+h)-}\bigr)-v(x)\ \Big|\ X^{\mu,u}_t=x\Bigr].
$$

$\mathcal R_h$ is precisely the **martingale / dynamic-programming residual**; setting it to zero is the discrete martingale condition of papers 86, 93 and 96 written for a controlled, mean-field-dependent process.

**Step three: a one-step random mapping replaces path simulation.** On the uniform grid $\Pi_h=\{nh:n=0,\dots,N-1\}$ (note $T\notin\Pi_h$ because of the reset), define the **random mapping**

$$
\Phi(x,\mu,u)\overset{\mathcal L}=
\begin{cases}
x+b(x,\mu,u(x))h+\sigma(x,\mu,u(x))\sqrt h\,\zeta, & x\in\Pi_h\times\mathbb R^d,\\
X_0, & x\in\{T\}\times\mathbb R^d,
\end{cases}
\qquad \zeta\sim N(0,I_q),
$$

that is, Euler-Maruyama inside a cycle and **resampling from the initial law** at the cycle boundary. **Remark 2** gives $\Phi$ a weak local truncation error of $O(h^2)$. The population measure is carried by a particle ensemble updated **one step per outer iteration, on a random mini-batch only**:

$$
X^m_i:=\begin{cases}\Phi(X^m_{i-1},\mu_{i-1},u_\alpha), & m\in A_i,\\ X^m_{i-1}, & m\in A_i^c,\end{cases}
\qquad
\mu_i:=\frac1M\sum_{m=1}^M\delta_{X^m_i}.
$$

**Remark 3 is essential for reading the algorithm: $i$ is the policy-iteration index, not a time index.** For each fixed $i$ the ensemble $\{X^m_i\}$ is a set of **space-time** samples of $[0,T)\times\mathbb R^d$ approximating the occupation measure $\mu^*$, not samples at one time level. **This is what removes the sequential-in-time bottleneck: no trajectory is ever simulated end to end.**

**Step four: weak form plus adversarial training.** Only _one_ sample of the next state is available per current state, so the conditional expectations are removed by testing:

$$
\text{(PE)}\ \min_\theta\sup_{\rho\in\mathcal T}
\bigl|\mathbb E\bigl[\rho(X^m_i)\,\mathcal M(X^m_i;\mu_i,u_\alpha,v_\theta)\bigr]\bigr|^2,
\qquad
\text{(PI)}\ \min_\alpha\ \mathbb E\bigl[\mathcal M(X^m_i;\mu_i,u_\alpha,v_\theta)\bigr],
$$

where $\mathcal M(x;\mu,u,v):=\{hf(x,\mu,u)+v\circ\Phi(x,\mu,u)-v(x)\}\mathbf 1(t\in\Pi_h)$ with $\mathcal R_h=\mathbb E[\mathcal M]+O(h^2)$, and $\mathcal T:=\{\rho:[0,T)\times\mathbb R^d\to[-1,1]^r\ \text{smooth}\}$. The justification is the tower property, $\mathbb E[\rho\mathcal M]=\mathbb E[\rho\,\mathbb E[\mathcal M|X^m_i]]$, so the policy-evaluation step is a **weighted Galerkin formulation** of the residual equation, weighted by the law of $X^m_i$. And by the law of total expectation, if $u_\alpha$ is expressive enough then minimising the _averaged_ objective makes $u_\alpha(X^m_i)$ minimise the inner conditional expectation for almost every realisation, so **the averaged optimisation preserves pointwise optimality** — the paper cites paper 86's Lemma 3.2 for the rigorous version. **Remark 5**: a vector-valued test function with $r\ge600$ substantially stabilises adversarial training compared with $r=1$.

**Step five: networks and minibatch estimators.** The control is clipped into the box $U=\prod_j[a_j,b_j]$ by $u_{\alpha,j}(x)=\min\{\max\{a_j,\psi_{\alpha,j}(x)\},b_j\}$. The value network hard-wires the terminal condition **using the current empirical measure**: $v_\theta=\varphi_\theta(x)$ on $[0,T)\times\mathbb R^d$ and $v_\theta=g(x,\mu_i)$ on $\{T\}\times\mathbb R^d$. The adversary is the same shallow multiscale sine network as in paper 93. The minibatch estimators again use **disjoint** sets $A_{i,1},A_{i,2}$ so that the squared term is (nearly) unbiased.

**Four stated advantages:** (i) no higher-order derivatives — the weak form never differentiates $v_\theta$ in $x$, so no Hessians; (ii) no full path simulation — one-step transitions, parallelisable across the time dimension; (iii) no pointwise Hamiltonian minimisation — replaced by averaged optimisation of $u_\alpha$; (iv) **shared computation** — the same transitions $\{X^m_i\}_{m\in A_i}\to\{X^m_{i+1}\}$ serve the policy-evaluation step, the policy-improvement step and the measure update, which is precisely the "one sampler for three jobs" problem identified in the introduction.

### Theorems

**The paper is constructive; there is no convergence theorem.** The conclusion says so: "Future work includes establishing rigorous convergence and error estimates, extending the framework to more general mean-field interactions and common noise, and developing adaptive sampling and network architectures."

The only quantitative estimates are local: **Remark 2**, the weak local truncation error $O(h^2)$ for $\Phi$; and $\mathcal R_h=\mathbb E[\mathcal M]+O(h^2)$, where the $O(h^2)$ comes from the left-rectangle quadrature plus the weak Euler error. There is no global rate, no bound in terms of $M$, $|A_i|$, network size or the number of policy iterations, and no proof that the fictitious-play-like outer loop converges to the equilibrium. **Remark 4** supplies the positioning: the iteration is a form of **fictitious play** — at iteration $i$, $u_\alpha$ best-responds to the population measure inherited from the previous iteration — differing from prior fictitious-play mean-field-game solvers in that the measure is advanced by a **one-step Markovian transition** instead of by re-simulating whole trajectories.

### Numerical experiments

Setup: $\psi_\alpha,\varphi_\theta$ fully connected, depth $H=6$, width $W=104$ for $d=1$ and $W=1008$ for $d=1000$, ReLU; adversary with $r=1200$ and $c=10$; learning rates $\delta_1=\delta_2=\delta_0\times10^{-3}\times0.01^{i/I}$ and $\delta_3=\delta_0\times10^{-2}\times0.01^{i/I}$ with $\delta_0=3d^{-0.5}$; $I=9000$ iterations; $J=2K=2$; full batch $M=1024\times10^3$ with mini-batches $|A_{i,1}|=|A_{i,2}|=M/20$; $h=T/100$; RMSProp; PyTorch 2.6.0 with float32 and automatic mixed precision, distributed data parallel on **8 NVIDIA RTX 4090 GPUs**. $Z_0\sim\mathrm{Uniform}(\{s\mathbf 1_d:s\in[-c,c]\})$ with $c=1$ in Section 4.1 and $c=3$ in Section 4.2. Metrics: the relative cost $\mathrm{RC}=|\hat J-J^*|/|J^*|$ plus $\mathrm{RE}_1$ and $\mathrm{RE}_\infty$ for the value function on 1000 fixed test points.

**(i) Linear-quadratic mean-field game with an explicit solution (Section 4.1, Table 1).** $b_z=u+c_0(\bar m_t-z)+c_1(z^*(t)-\bar m_t)$, $\sigma_z=c_\sigma I_d$, $f=\frac12\{c_2|u|^2+c_3|z-\bar m_t|^2+c_4|z-z^*(t)|^2\}$, $g=\frac{c_5}2|z-z^*(T)|^2+\frac12$, $T=1$, with $\bar m_t$ the population mean. The reference $v$ and $u^*$ come from Riccati-type ODEs solved with `scipy.integrate.solve_ivp`, RK45 at rtol $10^{-8}$; the appendix derives them. Variants: **LQ-1** ($c_0=c_1=c_4=0$, $c_2=1$, $c_3=c_5=1/d$, $c_\sigma=0.5/\sqrt d$, no target tracking); **LQ-2** ($c_0=1$, $c_1=0$, $c_2=c_3=c_4=c_5=1/d$, target the unit circle $z^*(t)=y^*(t)/|y^*(t)|$ with $y^*_i(t)=\sin(2\pi t+i\pi/2)$); **LQ-3** (as LQ-2 but with $c_1=-0.5$ and a helix target $z^*(t)=2t\,y^*(t)/|y^*(t)|$). The initial guess is deliberately adversarial: $\hat\mu$ is the law of $(t,Z_{{\rm init},t})$ with $Z_{{\rm init},t}=Z_0+5t\mathbf 1_d+c_\sigma B_t$, "intentionally far from the true $\mu^*$". Means over five runs with standard deviations in parentheses; MEM is peak per-GPU memory (times 8).

| Equation | $d$  | $\mathrm{RE}_1$   | $\mathrm{RE}_\infty$ | RC                | MEM (MB) | RT (s) |
| -------- | ---- | ----------------- | -------------------- | ----------------- | -------- | ------ |
| LQ-1     | 1    | 4.05E-3 (2.09E-3) | 1.26E-2 (5.20E-3)    | 9.69E-3 (6.18E-3) | 535      | 294    |
| LQ-1     | 1000 | 1.68E-2 (2.13E-3) | 4.33E-2 (5.99E-3)    | 1.55E-2 (1.26E-3) | 3964     | 579    |
| LQ-2     | 1    | 9.92E-3 (2.41E-3) | 2.11E-2 (3.91E-3)    | 9.84E-3 (3.24E-3) | 535      | 308    |
| LQ-2     | 1000 | 1.97E-2 (4.57E-3) | 2.86E-2 (5.18E-3)    | 1.95E-2 (3.12E-3) | 3995     | 602    |
| LQ-3     | 1    | 6.67E-3 (5.62E-3) | 1.67E-2 (7.65E-3)    | 6.54E-3 (5.23E-3) | 535      | 317    |
| LQ-3     | 1000 | 1.75E-2 (2.15E-3) | 2.88E-2 (5.22E-3)    | 1.79E-2 (3.01E-3) | 3997     | 615    |

**What this table establishes is the mildness of the dimension dependence: from $d=1$ to $d=1000$, $\mathrm{RE}_1$, $\mathrm{RE}_\infty$ and RC grow only by factors of about 2.0-4.1, 1.4-3.4 and 1.6-2.7 respectively, while runtime rises by about $2\times$ and memory by about $7.5\times$.** A separate **$d=10\,000$ run for LQ-1** shows close agreement of $v_\theta$ with $v$ along the diagonal and steady decay of the errors, with a **runtime of 5258 seconds and peak per-GPU memory of 55,645 MB** (times 8). The $d=1$ visualisation also shows that the value network's error is small in regions visited by the simulated paths and larger in unexplored regions — **the same phenomenon paper 93 proves via the weighted norm $M^2_n[\cdot]$.**

**(ii) Interbank systemic-risk mean-field game (Section 4.2).** A Carmona-Fouque-Sun style model, **without common noise** for compatibility with the scheme: $b_z=c_1(\bar m_t-z)+u$, $\sigma_z=c_2$, running cost $f=\frac12u^2-c_3u(\bar m_t-z)+\frac{c_4}2(\bar m_t-z)^2$, terminal $g=\frac{c_5}2(\bar m_T-z)^2$, $T=1$, $c_1=1$, $c_2=0.5$, $c_3=1$, $c_4=2$, $c_5=1$, scalar state. The reference again comes from Riccati ODEs. Reported diagnostics are the errors versus iteration and the learned empirical mean reserve $\hat m_t$ against the exact $m^*_t$.

**(iii) Target-tracking mean-field game with nonlinear interaction, $d=2$ (Section 4.3).** Chosen low-dimensional for visualisation: $b_z=u$, $\sigma_z=0.5I_d$, $f=c_1|z-z^*(t)|^2+c_2|u|^2+F(x,\mu)$, $g=0.5c_4|z-z^*(T)|^2$, circular target $z^*_i(t)=\sin(2\pi t+i\pi/2)$, and the **congestion / anti-aggregation term**

$$
F(x,\mu):=\int_{[0,T)\times\mathbb R^d}\delta_t(s)\exp\bigl(-c_F|z-y|^2\bigr)\,\mu(\mathrm ds\times\mathrm dy),
$$

with $c_1=1$, $c_2=0.1$, $c_4=10$, $c_F=1$. There is no explicit solution, so validation is qualitative: the empirical cost drops by about 0.1 over training and the agents trace an approximately circular pattern ending near the target.

**(iv) Mean-field game with nonlinear dynamics and interactions, $d=2$ (Section 4.4).** $b_z=u$, $\sigma_z=\frac{c_1}{\sqrt d}\sum_{i=1}^d\sin(i+z_i)I_d$ (a state-dependent, non-constant diffusion), $f=c_0/(1+c_2|z-0.5\mathbf 1_d|^2)+c_3|u|^2+F(x,\mu)$, $g=0.5c_4|z-z^*(T)|^2$, $T=1$, $c_1=c_3=0.1$, $c_2=c_4=1$, $c_0=5$, $c_F=1$, with $Z_0$ standard Gaussian. The first term of $f$ creates a **barrier** at $0.5\mathbf 1_d$; the agents are observed to spread, route around the barrier and converge on the target.

### Relation to the others

**It is the mean-field-game member of the Cai-Fang-Zhou martingale and adversarial family.** The machinery is inherited almost intact from papers 86, 93 and 96: the martingale / dynamic-programming residual $\mathcal M$, the Galerkin weak form with an adversarial test network, the multiscale sine adversary with large output dimension, the disjoint mini-batches for an unbiased squared expectation, the box-clipped control network, and the terminal condition hard-wired into $v_\theta$. Remark 6 credits the weak-adversarial-network idea, and the pointwise-optimality argument for the averaged improvement step cites paper 86's Lemma 3.2.

**What is genuinely new relative to papers 86, 93 and 96** is the **regenerative reformulation** with cycles and resetting, which converts a finite-horizon mean-field game into an infinite-horizon problem with an invariant occupation measure, and the **one-step random mapping** update of the particle ensemble. Together these make the outer loop index double as the cycle index, so no full trajectory is ever simulated — the mean-field-game analogue of the offline-path-generation trick of papers 96 and 93, pushed further.

**Mean-field connections:** paper 33 (explicit $\theta$-schemes for mean-field BSDEs) and paper 61 (explicit multistep schemes for mean-field FBSDEs) are the high-order, low-dimensional, provably convergent treatments of mean-field _backward_ equations; paper 97 (DeepSPoC) handles the mean-field _forward_ side; paper 108 couples both sides and solves the equilibrium problem at $d=10^4$ but with no theory. **Papers 33, 61, 97 and 108 span exactly the same grid of trade-offs — order and proof versus dimension — that papers 8, 47 and 63 versus 86, 93 and 96 span for the non-mean-field case.**

The link to the stochastic-control thread: papers 26, 41 and 50 solve control problems through the Pontryagin/FBSDE route (gradient projection, high-order FBSDE schemes, conditional gradient with particle filters), whereas paper 108 goes through dynamic programming and policy iteration and adds the equilibrium fixed point. On authorship, Shuixin Fang and Tao Zhou also co-author papers 86, 93, 96 and 100; Zhen Wu is a control theorist at Shandong University, which is also Weidong Zhao's institution — the classical FBSDE co-author on papers 8, 16, 18, 19, 23, 25, 26, 33, 35, 41, 47, 61, 63 and 68.

## A general judgement about this line

Martingale methods change **how the residual is tested**. Checking a residual pointwise is unavailable in high dimension because there is no grid; checking whether a process is a martingale needs only expectations along simulated paths, so dimension no longer enters the cost directly.

Paper 86 carries the conversion out most completely, and it can be summarised as three substitutions:

1. pointwise minimum principle to integral minimum principle (removing the curse of dimensionality in $(t,x)$);
2. PDE residual to martingale property (removing the grid);
3. conditional expectation to a weak condition against test functions (removing the computation of the conditional expectation).

Each substitution replaces a condition that cannot be handled directly by an equivalent one that can be verified by sampling. **The third is the most noteworthy**: it converts "compute a conditional expectation" into "train a discriminator network", at the cost of importing all the instability of adversarial training and with the benefit of avoiding conditional-expectation estimation entirely.

The other five each remove one component of that machinery, and the genealogy of removals is itself informative:

| Paper | What it removes                                         | What it buys                                                                                   | What it costs                                                        |
| ----- | ------------------------------------------------------- | ---------------------------------------------------------------------------------------------- | -------------------------------------------------------------------- |
| 96    | all derivatives in the loss, and path re-simulation     | $d=10^4$, parallelism in time and space                                                        | still no global theorem; controlled jumps cannot be pre-computed     |
| 93    | stochastic analysis (Taylor and moments instead)        | the first convergence-rate theorem (first order in time), which retroactively covers 86 and 96 | lower accuracy than PINNs; the RDO carries its own $O(h)$ truncation |
| 100   | the neural network itself                               | a rigorous a priori bound, laptop-scale cost, interpretability                                 | only $d=10^4$, only first order                                      |
| 97    | stored particle trajectories                            | forward mean-field equations, and the only computable posterior bound                          | the convergence theorem holds only in a Fourier surrogate            |
| 108   | full trajectory simulation and the coupled HJB-FP solve | mean-field game equilibria at $d=10^4$                                                         | no convergence theorem at all                                        |

**One limitation runs across the whole table and deserves to be remembered separately**: the criterion of this family only forces the PDE residual to vanish **in the region the pilot process or the sampled ensemble explores**. Remark 1 of paper 96 states it, Theorem 4 of paper 93 encodes it in the weighted norm $M^2_n[\cdot]$, and paper 108 displays it in its $d=1$ visualisation. **"Solved accurately" is always a statement with a domain attached here.**

A second judgement concerns the scale of accuracy. The classical multistep and $\theta$-schemes reach sixth order at $d\lesssim10$, while the papers on this page reach first order and relative errors of $10^{-2}$ to $10^{-3}$ at $d=10^4$. **These are not the same task done well or badly; they are two different tasks.** Paper 93's Table 4 is especially honest about it: its method is less accurate than PINNs, all the gain is in time and memory, and the gain grows with $d$ — a case of trading accuracy for feasibility with the price written on the label.

## Coverage checklist

| Item                                                                   | Paper | Coverage status                                |
| ---------------------------------------------------------------------- | ----- | ---------------------------------------------- |
| HJB equation, Hamiltonian, and the four covered situations             | 86    | verified from the full text                    |
| Integral minimum principle (Lemma 3.2) with its one-line proof         | 86    | verified from the full text                    |
| Martingale characterisation (Lemma 3.4), both directions               | 86    | verified from the full text                    |
| Theorem 3.5, the geometric assumption, weak versus strong order        | 86    | verified; noted as not a convergence theorem   |
| Adversarial weak form, augmented Lagrangian, three networks            | 86    | verified from the full text                    |
| Experiments: runtimes, first-order rate, Tables 1 and 2                | 86    | verified, tabulated here                       |
| Random-difference expansion, RDO, Fokker-Planck sampling               | 93    | verified from the full text                    |
| The two equivalences (martingale methods, implicit Euler)              | 93    | verified from the full text                    |
| Assumptions 1-3, Theorem 1, Corollary 2, Lemma 3, Theorem 4            | 93    | verified, with step restrictions and constants |
| Remark 9 (proving first order for papers 86 and 96 as well)            | 93    | verified from the full text                    |
| Experiments: Tables 2, 3, 4 and the authors' self-assessment           | 93    | verified, tabulated here                       |
| Pilot/system process, localised martingale, Remark 1                   | 96    | verified; the region limitation preserved      |
| Derivative-free $\mathcal M$, disjoint minibatches, policy improvement | 96    | verified from the full text                    |
| Only an $O(h^2)$ local estimate, no global theorem                     | 96    | verified from the full text                    |
| Experiments: Tables 1 and 2, and the width study                       | 96    | verified, tabulated here                       |
| SPoC, the deepSPoC system, three losses, adaptive sampling             | 97    | verified from the full text                    |
| Theorem 3.3 and its surrogate-setting caveat                           | 97    | verified from the full text                    |
| Posterior estimate (Proposition 3.5, Theorem 3.6)                      | 97    | verified from the full text                    |
| Experiments: five equation families with their parameters              | 97    | verified; no orders reported, noted            |
| Local linear regression, matrix-free solve, Newton iteration           | 100   | verified from the full text                    |
| Lemmas 3.1, 3.3, 3.8 and Theorem 3.1                                   | 100   | verified from the full text                    |
| Experiments: three problems, runtime tables, hardware                  | 100   | verified, tabulated here                       |
| Regenerative reformulation, one-step random mapping, Remark 3          | 108   | verified from the full text                    |
| No convergence theorem, only an $O(h^2)$ local estimate                | 108   | verified from the full text                    |
| Experiments: Table 1, the $d=10^4$ run, four test problems             | 108   | verified, tabulated here                       |

## Sources for this page

- W. Cai, S. Fang, and T. Zhou, [_SOC-MartNet: a martingale neural network for the Hamilton-Jacobi-Bellman equation without explicit inf H in stochastic optimal controls_](https://doi.org/10.1137/24M1681033), SIAM J. Sci. Comput. 47(4) (2025), pp. C795-C819 (preprint [arXiv:2405.03169](https://arxiv.org/abs/2405.03169); code [sx-fang/MartNet](https://github.com/sx-fang/MartNet)).
- W. Cai, S. Fang, and T. Zhou, [_Deep random difference method for high-dimensional quasilinear parabolic partial differential equations_](https://doi.org/10.1016/j.jcp.2026.114767), J. Comput. Phys. 555 (2026), 114767 (preprint [arXiv:2506.20308](https://arxiv.org/abs/2506.20308); code [sx-fang/DRDM](https://github.com/sx-fang/DRDM)).
- W. Cai, S. Fang, W. Zhang, and T. Zhou, _Martingale deep learning for very high dimensional quasi-linear partial differential equations and stochastic optimal controls_, [arXiv:2408.14395](https://arxiv.org/abs/2408.14395), submitted to SIAM Rev.
- K. Du, Y. Xie, T. Zhou, and Y. Zhou, _DeepSPoC: a deep learning based sequential propagation of chaos_, [arXiv:2408.16403](https://arxiv.org/abs/2408.16403), submitted to SIAM/ASA J. Uncertain. Quantif.
- S. Fang, C. Sheng, B. Su, and T. Zhou, _A derivative-free localized stochastic method for very high dimensional semi-linear parabolic PDEs_, [arXiv:2510.02635](https://arxiv.org/abs/2510.02635), submitted to Numer. Math.
- S. Fang, S. Wang, Z. Wu, H. Zhang, and T. Zhou, _Deep policy iteration for high-dimensional mean-field games with regenerative reformulation_, preprint arXiv:2604.26782.
