---
title: Spectral Bias and Generative Solvers
description: Papers 81, 89, 94, 101, 103 and 105 - frequency content should be measured, not assumed
lang: en
translation: computational-mathematics/paper-notes/scientific-machine-learning/spectral-bias-and-generative-solvers
tags:
  - paper-notes
  - scientific-machine-learning
  - spectral-bias
---

> [!note] Coverage of this page
> Papers **81** (_Comput. Methods Appl. Mech. Engrg._ 437, 2025), **89** (_Neural Networks_ 194, 2026, [arXiv:2401.02080](https://arxiv.org/abs/2401.02080)), **94** (_J. Comput. Phys._ 558, 2026), **101** ([arXiv:2512.18586](https://arxiv.org/abs/2512.18586)), **103** ([arXiv:2606.22514](https://arxiv.org/abs/2606.22514)) and **105** ([arXiv:2604.07169](https://arxiv.org/abs/2604.07169)).

## Spectral bias is a measurable quantity

Spectral bias, or the frequency principle, is the observation that networks converge rapidly to low-frequency components and struggle to represent high-frequency or highly oscillatory features. **That direction is not absolute**: paper 101 shows it reverses wholesale under a PDE residual loss, as set out below. Taking the usual direction first, multi-scale networks improve on it with a down-scaling map: decompose the compactly supported Fourier transform over concentric annuli $\mathbb K_i=\{\bm k:(i-1)K_0\le|\bm k|\le iK_0\}$, scale each band down by $a_i>1$, and approximate

$$
f_{\theta}(\bm x)=\sum_{i=1}^{W} a_i^{d}\, f_{\theta_i}(a_i\bm x),
$$

with $a_i$ normally fixed a priori at $2^{i-1}$.

## 81: prove an error theory first, then design from it

### Two concrete failure examples

The paper's motivating example is $f(x)=\sin(\pi x)+\sin(100\pi x)$, whose ideal scales are $a_1=\pi$ and $a_2=100\pi$ — unavailable in practice. The second problem is that a network with (random) Fourier features is more accurate than a multi-scale network when the embedding frequency is **exactly right** and **worse than a plain network** when it is slightly off, which the paper demonstrates concretely with $k=38\pi$ against a target frequency of $40\pi$.

### Two error bounds

The first explains why down-scaling helps. For a $(p,C_0)$-smooth $f:\mathbb R^d\to\mathbb R$ with $p=q+s$, $q\in\mathbb N_0$, $s\in(0,1]$, and $M$ sufficiently large, there is a network with the **down-scaling** map $\Phi(\bm x)=k\bm x$ such that

$$
\|f_{\text{net}}-f\|_{\infty,[-h,h]^d}
\ \le\ \frac{C_2\bigl(\max\{kh,\,C_1\}\bigr)^{5q+3}}{M^{2p}} .
$$

Compared with the bound for a standard network without a scaling map, the point is that **shrinking $kh$ shrinks the numerator**.

The second explains why the true dominant modes should be used as features. If $f\in L^2([-\pi,\pi]^d)$ is band-limited, $f(\bm x)=\sum_{\bm k\in\mathbb B}[b_{\bm k}\cos(\bm k\cdot\bm x)+c_{\bm k}\sin(\bm k\cdot\bm x)]$ with $|\mathbb B|=N$ and bounded coefficients, then with Fourier features $\Phi(\bm x)=(\dots,\sin(\bm k\cdot\bm x),\cos(\bm k\cdot\bm x),\dots)_{\bm k\in\mathbb B}$ there is a network with

$$
\|f_{\text{net}}-f\|_{\infty,[-\pi,\pi]^d}\ \le\ \frac{C_2\sqrt N}{M^{2p}} .
$$

**This bound has no explicit frequency dependence.** That is the theoretical justification for building the embedding from the actual dominant modes.

### A hybrid feature embedding

The paper's new input map combines the linear down-scaling of multi-scale networks with a sinusoidal Fourier feature:

$$
\Phi[\bm k](\bm x)=\begin{bmatrix}\bm k\cdot\bm x\\ \cos(\bm k\cdot\bm x)\\ \sin(\bm k\cdot\bm x)\end{bmatrix}.
$$

The motivation is a concrete computation: approximating $f(x)=\sin(40\pi x)$ through $F(y)=\sin(\hat k\arcsin y)$ with $\hat k=40\pi/k$ gives $F'(y)=\hat k\cos(\hat k\arcsin y)/\sqrt{1-y^2}$, so $\lim_{y\to1}|F'(y)|=\infty$ and therefore $\|F\|_{C^1[-1,1]}=\infty$ whenever $|\hat k|\ne1$, at which point the relevant theorem's hypothesis fails. **That is the precise reason pure Fourier features degrade under frequency mismatch**, and adding the linear component is what repairs it.

### Capturing the dominant modes a posteriori

Take the discrete Fourier transform of the current network solution, $\hat f_{\text{net},\bm k}=\int_{[0,1]^d}f_{\text{net}}(\bm x)e^{-\mathrm i2\pi\bm k\cdot\bm x}\mathrm d\bm x$, and keep the $N_0$ coefficients of largest modulus with $N_0$ fixed by an energy-fraction condition:

$$
\sum_{j=1}^{N_0}\bigl|\hat f_{\text{net},\bm k_j}\bigr|^2
\ \ge\ (1-\delta)\,\|f_{\text{net}}\|^2_{L^2([0,1]^d)},
\qquad 0\le\delta<1 .
$$

The network is then rebuilt by one of two criteria. If the number of captured modes satisfies $N_0\le M_0$ (the number of sub-networks), build $N_0$ sub-networks whose $j$-th input is $\Phi[\bm k_j](\bm x)$ with output

$$
y=\sum_{j=1}^{N_0} h_j\,y_j,\qquad h_j=\hat u_{\text{net},0,\bm k_j},
$$

reusing the coefficients to accelerate convergence. If $N_0>M_0$, do **not** grow the network: partition $\mathbb B_1=\{\bm k_1,\dots,\bm k_{N_0}\}$ into $M_0$ blocks, feed sub-network $j$ the inputs $\Phi[\bm k](\bm x)$ for $\bm k$ in its block, and output

$$
y=\sum_{j=1}^{M_0}h_j\,y_j,\qquad h_j=\sum_{\bm k\in\mathbb B^j_1}\bigl|\hat u_{\text{net},0,\bm k}\bigr| .
$$

One implementation detail is worth recording: the final layer is written $u_{\text{net}}=WG+b$ with $W=(h_1,h_2,\dots)$ and $b=0$ **fixed and not learnable**, which the paper reports is both more accurate and cheaper.

The algorithm has two stopping criteria: a fixed maximum number of rounds, and stagnation of the captured feature set. The reported improvement is two to three orders of magnitude in accuracy over standard multi-scale networks, on a Poisson equation, a heat equation, a wave equation and a Schrödinger equation near the semi-classical limit.

## 94: the semi-classical limit, removing the spatial variable entirely

### The scale of the problem

The target equation is

$$
\psi_t=\frac{\mathrm i\varepsilon}{2}\Delta\psi-\frac{\mathrm i}{\varepsilon}V(\bm x)\psi,
\qquad
\psi(\bm x,0)=\varphi(\bm x)\exp\bigl(\mathrm i\phi(\bm x)/\varepsilon\bigr),
$$

with $\varepsilon$ the non-dimensional Planck constant. For $\varepsilon\ll1$ the solution oscillates at scale $O(\varepsilon^{-1})$ in both space and time. Classical splitting solvers (Strang splitting with spectral accuracy in space and $O(\Delta t^2/\varepsilon)$ in time; Chin-Chen with spectral space and fourth-order time) all have global truncation errors scaled by a factor $1/\varepsilon$, so both mesh and time step must shrink with $\varepsilon$.

Plain physics-informed networks are worse still. The paper supplies the evidence directly: applying a vanilla network to this equation on $(-2\pi,2\pi)\times[0,1]$ gives a relative $L^2$ error already **exceeding 1** at $\varepsilon=0.1$.

### Heller's Gaussian wave packet and the decisive reduction

The paper adopts Heller's Gaussian wave packet ansatz, in one dimension

$$
\psi(x,t)=\exp\left[\frac{\mathrm i}{\varepsilon}
\Bigl(\alpha(t)\bigl(x-q(t)\bigr)^2+p(t)\bigl(x-q(t)\bigr)+\gamma(t)\Bigr)\right].
$$

The envelope is centred at $q(t)$ with standard deviation proportional to $\sqrt{\varepsilon/\alpha_{\rm im}(t)}$, the oscillation wavelength at $x=q(t)$ is $2\pi\varepsilon/p(t)$, and $\mathrm{Re}\,\alpha$ produces finer oscillations in the tails.

**The decisive point is that substituting this ansatz produces a system of ordinary differential equations for $(\alpha,q,p,\gamma)$, with the spatial variable $x$ removed entirely.** What the network must learn is therefore no longer a function oscillating rapidly in space but a handful of scalars evolving in time: the high-frequency structure has been placed analytically inside the ansatz. This is the [[computational-mathematics/paper-notes/scientific-machine-learning/index|topic's]] "write the structure into the architecture" stance realised for an oscillatory problem.

A second, independent problem the paper addresses is that changing the initial condition forces retraining, which motivates an operator-learning formulation.

## 101: a fixed feature bank with learnable weights

### An observation that points the other way

The most notable thing in this paper is not its architecture but an observation in its Section 3.1: **under physics-based training the direction of the spectral bias reverses.** The reason is scaling on the Fourier side. Since $\widehat{\partial_xu}(k)=\mathrm ik\,\hat u(k)$, a Deep Ritz energy term $\|\partial_xu\|_{L^2}^2$ weights mode $k$ by $k^2$; since $\widehat{u_{xx}}(k)=-k^2\hat u(k)$, a squared PDE residual weights mode $k$ by $k^4$. The differential operator therefore **amplifies** high frequencies on its own, and what ends up under-resolved under a residual loss is the **low**-frequency part.

The paper measures this. For $u=\sin(\pi x)+\sin(5\pi x)+\sin(20\pi x)$ it tracks the frequency-wise relative error $\Delta_F(k)=|\hat u_k^{\rm pred}-\hat u_k|/|\hat u_k|$ at $k\in\{1,5,20\}$: under pure regression the low frequencies do converge first, but under Deep Ritz and PINN losses the high frequencies decay faster.

That fixes the form of the solution representation — a sum of two networks,

$$
u(\bm x;\theta)=u_h(\bm x;\theta_h)+\alpha\,u_\ell(\bm x;\theta_\ell),
$$

where $u_h$ is the cross-attention network below (responsible for high frequencies) and $u_\ell$ is a plain fully connected network (low frequencies), with the boundary condition split as $u_h=g$ and $u_\ell=0$.

### Fixed frequencies, learnable decay

Base frequencies $\omega_m\sim\mathcal N(\bm 0,\sigma^{-2}I_{d_{\rm in}})$ are expanded over dyadic scales $\widetilde\omega_{m,k}=2^k\omega_m$ for $k=0,\dots,K$, giving $M=M_{\rm base}(K+1)$ frequencies, with phases $b_{m,k}\sim\mathrm{Uniform}(0,2\pi)$ drawn once. The feature map carries a **learnable amplitude envelope**:

$$
\phi(\bm x)=\sqrt{\tfrac1M}\Bigl[a_{m,k}\cos\bigl(\widetilde\omega_{m,k}^{\top}\bm x+b_{m,k}\bigr)\Bigr]_{(m,k)},
\qquad
a_{m,k}=\exp\bigl(-\beta\|\widetilde\omega_{m,k}\|_2\bigr).
$$

**The frequencies themselves are fixed and non-trainable; the only learnable quantity is the decay rate $\beta$**, parameterised through a softplus to keep it non-negative. This is exactly complementary to paper 81's trade-off: 81 moves the features onto the true frequencies, 101 pairs a sufficiently wide fixed bank with a learnable decay envelope.

### Cross-attention: queries from the solution, keys and values from the bank

Write $Q^{(l)}(\bm x)$ for the latent state and $H(\bm x)$ for the bank features, and set $Q_l=Q^{(l)}W_Q^{(l)}$, $K_l=H(\bm x)W_K^{(l)}$, $V_l=H(\bm x)W_V^{(l)}$. Then

$$
\mathrm{CA}\bigl(Q^{(l)},H\bigr)=\mathrm{softmax}\Bigl(\frac{Q_lK_l^{\top}}{\sqrt{d_q}}\Bigr)V_l,
$$

$$
\widetilde Q^{(l)}=Q^{(l)}+\mathrm{CA}\bigl(Q^{(l)},H\bigr),
\qquad
Q^{(l+1)}=\widetilde Q^{(l)}+\sigma\bigl(W^{(l)}\widetilde Q^{(l)}+b^{(l)}\bigr).
$$

The query comes from the **latent state** and the keys and values from the **frequency bank** — hence cross-attention rather than self-attention. The point is that the softmax weights are input-dependent, so different regions of the domain can emphasise different frequency bands.

### A posteriori frequency enhancement

Train a preliminary model, take a DFT on a uniform grid to get $\hat u_{\theta,k}$, put $\zeta=\max_{k\in B}|\hat u_{\theta,k}|$ and select by relative threshold

$$
\mathcal K_{\rm post}=\bigl\{k\in B:\ |\hat u_{\theta,k}|>\lambda\zeta\bigr\},\qquad 0<\lambda<1 .
$$

Deterministic frequencies $\omega_k^{\rm post}=2k\pi$ and their features are built from this set and concatenated in the **token** dimension as $H_{\rm aug}=[H_{\rm base};H_{\rm post}]$. The new tokens are not switched on abruptly but released smoothly through an additive logit mask,

$$
A^{(l)}=\frac{Q_lK_l^{\top}}{\sqrt{d_q}}+\mathcal M^{(l)},
\qquad
\mathcal M^{(l)}=[\,\bm 0;\ \eta_l\bm 1\,],\quad \eta_l\le0,
$$

with the zero block acting on $H_{\rm base}$ and the constant block $\eta_l$ on $H_{\rm post}$, and a schedule $\eta_l\uparrow0$ releasing the injected tokens gradually. So the bank is **augmented** rather than rebuilt and the backbone is left untouched — which is exactly where this parts company with paper 81.

## 103: operator splitting as a network, with time back in the input

Take an autonomous evolution equation

$$
u_t=\mathcal Lu+\mathcal Nu=:\mathcal Fu\ \text{ in }\Omega\times(0,t^\star],
\qquad u(0,\cdot)=u_0,\qquad \mathcal Bu=0\ \text{ on }\partial\Omega,
$$

with $\mathcal L$ linear, $\mathcal N$ nonlinear and $\mathcal B$ the boundary operator, none of them depending on $t$. Lie-Trotter and Strang splitting give

$$
u(T,\bm x)\approx e^{\tau_K\mathcal N}e^{\tau_K\mathcal L}\cdots e^{\tau_1\mathcal N}e^{\tau_1\mathcal L}u(0,\bm x),
$$

$$
u(T,\bm x)\approx e^{\frac{\tau_K}2\mathcal L}e^{\tau_K\mathcal N}e^{\frac{\tau_K}2\mathcal L}\cdots e^{\frac{\tau_1}2\mathcal L}e^{\tau_1\mathcal N}e^{\frac{\tau_1}2\mathcal L}u(0,\bm x).
$$

DOSnet, the architecture this paper builds on, writes that product directly as a network: $\psi_{\bm\theta_T}=\psi_{\bm\theta_K}\circ\cdots\circ\psi_{\bm\theta_1}$, each splitting block alternating learnable linear layers (convolutions) with nonlinear layers $\phi_{\mathcal N_{l,i}}=e^{\tau_{l,i}\mathcal N}$ subject to $\sum_{l,i}\tau_{l,i}=T$. **The striking part is that the activation function is the exact flow of the nonlinear subproblem, not ReLU or tanh.**

Two defects motivate the paper. First, data-driven operator learning needs many paired samples, and for an evolution problem every pair costs a solver run. Second, the learned operator is never required to satisfy the equation, so it **cannot be evaluated at an arbitrary time**: DOSnet outputs only the terminal-time solution, its intermediate block outputs correspond to intermediate times only by comparison against reference data, and the number of retrievable intermediate times is tied to the number of blocks.

PI-DOSnet's central change is to take equal steps $\tau_{1,1}=\dots=\tau_{1,K}=dt=t/K$ and replace the exponential of the linear part by an explicit second-order Taylor expansion in $dt$, which makes the whole block an **explicit function of $t$**. Time returns to the input, the solution can be evaluated at any instant, and training no longer depends on paired data.

## 105: one shared summary network tying filtering to smoothing

### Setting

For a state-space model

$$
u_t=f(u_{t-1},\epsilon_{u,t}),\qquad y_t=h(u_t,\epsilon_{y,t}),
$$

the goal is to do **filtering** $p(u_t\mid y_{1:t})$ and **smoothing** $p(u_{1:t}\mid y_{1:t})$ at once. The paper constrains itself more than usual: it assumes only a **simulator**, so the functional forms of $f$, $h$ and the noise distributions are not assumed known and the transition and observation densities need not be evaluable.

The central difficulty is that the conditioning variable for filtering is the entire history $y_{1:t}$, whose dimension grows with $t$, so a fixed-input flow cannot be used. FLUID compresses it with a multi-layer LSTM into a fixed-length summary $s_t=\mathrm{Enc}(y_{1:t};\psi)\in\mathbb R^h$ and writes, with a conditional KRnet flow,

$$
p(u_t\mid y_{1:t})\approx p_{\theta_1,\psi}(u_t\mid s_t).
$$

### A causal factorisation, and sharing

The smoothing side rests on an exact factorisation

$$
p(u_{1:t}\mid y_{1:t})=p(u_t\mid y_{1:t})\prod_{k=1}^{t-1}p(u_k\mid u_{k+1},y_{1:k}),
$$

where the second factor conditions on $y_{1:k}$ and **not** on $y_{1:t}$ — that is what keeps the recursion causal, so the same summary $s_k$ can be reused directly. A second flow learns $p(u_t\mid u_{t+1},y_{1:t})\approx p_{\theta_2,\psi}(u_t\mid u_{t+1},s_t)$.

**Both flows share the same summary network $\psi$, and this is the paper's named contribution.** The joint objective is

$$
\min_{\theta_1,\theta_2,\psi}\
-\frac1{NT}\sum_{i=1}^N\sum_{t=1}^T\log p_{\theta_1,\psi}(u_t^i\mid s_t^i)
-\frac{\lambda}{N(T-1)}\sum_{i=1}^N\sum_{t=1}^{T-1}\log p_{\theta_2,\psi}(u_t^i\mid u_{t+1}^i,s_t^i).
$$

At $\lambda=(T-1)/T$ there is an algebraic identity: the two terms recombine into a per-time filtering likelihood plus a whole-trajectory smoothing likelihood, so this **single** loss is the maximum-likelihood objective for both tasks simultaneously. Sharing has a sufficiency justification as well: if a sufficient summary $S^{\dagger}$ exists, meaning $p(X\mid Y)=p(X\mid S^{\dagger})$ almost surely, then it is simultaneously the optimal summary for each of the two objectives.

### Sharing is not helpful but necessary

The ablation is the result most worth recording. At state dimension $K=20$, shared against independent summaries:

| Metric               | Shared | Independent |
| -------------------- | ------ | ----------- |
| filtering RMSE       | 0.1945 | 0.1958      |
| backward-kernel RMSE | 0.1240 | 0.1596      |
| smoothing RMSE       | 0.1525 | **60.3744** |
| smoothing CRPS       | 0.0742 | **33.2893** |

Each flow **on its own** is barely affected, while the backward **iterative sampling** collapses entirely. The paper's explanation is that with independent summaries the errors accumulate rapidly through the smoothing recursion and the estimates become unusable for $K\ge20$. This is the empirical counterpart of the loss-recombination identity above, and it is where the paper's real contribution lies.

> [!warning] Two places where the text and the tables disagree
> The prose claims smoothing improves uniformly on filtering, but Table 9 at $K=50$ reports smoothing RMSE 0.5423 against filtering 0.2605, and the Burgers Table 6 has the same inversion at $r^2=0.25$. This page records the table values.

## 89: sampling from a known unnormalised energy

The goal is to sample from $\pi(x)=Z^{-1}\exp(-U(x))$, where $U$ can be evaluated pointwise, $Z$ is unavailable, and **no training data from $\pi$ exist**. The last condition rules out the objectives used for ordinary generative models: Jensen-Shannon divergence, MMD and Wasserstein distance all need real samples.

The paper targets a specific defect in each of two families. Normalizing-flow Boltzmann generators must use a **bijective** decoder to have $\log p_D(x)$ in closed form, and bijectivity limits effective capacity; diffusion-based samplers (PIS and its relatives) need a numerical SDE or ODE solver inside the training loop to evaluate a time integral, which is expensive. The energy-based diffusion generator removes both restrictions at once: a **non-invertible** decoder together with a **simulation-free** loss.

The starting point is a variational bound. Since samples from $\pi$ are unavailable the reverse KL must be used, and it is bounded by an augmented KL carrying a latent variable:

$$
D_{\rm KL}\bigl(p_D(x)\,\|\,\pi(x)\bigr)
\le
\mathbb E_{p_D(z_0)p_D(x\mid z_0;\phi)}
\left[\log\frac{p_D(z_0)p_D(x\mid z_0;\phi)}{p_E(z_0\mid x;\theta)}+U(x)\right]+\log Z,
$$

with equality when $p_E(z_0\mid x;\theta)$ matches the decoder-induced conditional of $z_0$ given $x$. The decoder is Gaussian,

$$
p_D(x\mid z_0;\phi)=\mathcal N\bigl(x\mid\mu(z_0;\phi),\,\Sigma(z_0;\phi)\bigr),
$$

with $\mu$ and $\Sigma$ networks and **no invertibility imposed**; after decoding, a forward diffusion $\mathrm dz_t=f(z_t,t)\mathrm dt+g(t)\mathrm dW_t$ runs in latent space. The point is that $\log Z$ is a constant independent of the parameters, so the bound can be optimised without knowing it.

This shares the stance of the [[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|density-flow line]] with the setting exactly reversed: there one has samples and wants a density, here one has a density and wants samples.

> [!note] Coverage status
> The constructions, losses and principal results of all six have been checked equation by equation against the preprint or journal full text. What remains unchecked is the analysis of high-frequency amplification in paper 101's Appendix A — the Fourier scaling argument in its main text has been checked — together with the journal status of papers 101, 103 and 105, none of which carries journal information on its arXiv record.

## Where the six sit

| No. | How frequency content is handled                       | Change to the structure         |
| --- | ------------------------------------------------------ | ------------------------------- |
| 81  | dominant modes captured a posteriori from a DFT        | rebuilt (two criteria)          |
| 94  | absorbed analytically into a Gaussian wave packet      | not applicable                  |
| 101 | fixed bank, learnable envelope, cross-attention        | tokens augmented, backbone kept |
| 103 | not applicable (splitting structure, $t$ in the input) | unchanged                       |
| 105 | not applicable (conditional flows, shared summary)     | unchanged                       |
| 89  | not applicable (Boltzmann sampling)                    | unchanged                       |

One judgement runs through all of them: **frequency content is either measured (paper 81's DFT capture, paper 101's a posteriori threshold), removed analytically (paper 94's wave-packet ansatz), or weighted by the model itself over a sufficiently wide fixed bank (paper 101's cross-attention). Anything that fixes the frequency content a priori pays for frequency mismatch — and paper 81 turns that cost into a concrete quantity via $\|F\|_{C^1[-1,1]}=\infty$.**

Paper 101 then adds a qualification to that judgement: **"biased towards low frequencies" is a property of the loss, not an intrinsic property of the network.** Under a regression loss the network favours low frequencies; under a residual loss weighted like $k^4$ it favours high ones. What has to be avoided is therefore not a bias in one particular direction but assuming the direction without measuring it.

## Sources for this page

- J. Huang, R. You, and T. Zhou, [_Frequency-adaptive multi-scale deep neural networks_](https://doi.org/10.1016/j.cma.2025.117751), Comput. Methods Appl. Mech. Engrg. 437 (2025), 117751 (preprint [arXiv:2410.00053](https://arxiv.org/abs/2410.00053)).
- Y. Wang, L. Guo, H. Wu, and T. Zhou, [_Energy-based diffusion generator for efficient sampling of Boltzmann distributions_](https://doi.org/10.1016/j.neunet.2025.108126), Neural Networks 194 (2026), 108126 (preprint [arXiv:2401.02080](https://arxiv.org/abs/2401.02080)).
- J. Huang, R. You, and T. Zhou, [_Deep learning for the semi-classical limit of the Schrödinger equation_](https://doi.org/10.1016/j.jcp.2026.114869), J. Comput. Phys. 558 (2026), 114869 (preprint [arXiv:2509.04453](https://arxiv.org/abs/2509.04453)).
- X. Feng, T. Tang, X. Wan, and T. Zhou, _Overcoming spectral bias via cross-attention_, [arXiv:2512.18586](https://arxiv.org/abs/2512.18586), submitted to J. Comput. Phys.
- J. Huang, Y. Qian, and T. Zhou, _PI-DOSnet: a physics-informed deep operator-splitting network for evolution partial differential equations_, [arXiv:2606.22514](https://arxiv.org/abs/2606.22514), submitted to J. Comput. Phys.
- T. Cui, X. Feng, C. Pei, X. Wan, and T. Zhou, _FLUID: flow-based unified inference for dynamics_, [arXiv:2604.07169](https://arxiv.org/abs/2604.07169), submitted to Comput. Methods Appl. Mech. Engrg.
