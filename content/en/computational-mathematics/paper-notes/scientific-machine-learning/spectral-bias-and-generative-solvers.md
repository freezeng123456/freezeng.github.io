---
title: Spectral Bias and Generative Solvers
description: Papers 81, 89, 94, 101, 103 and 105 - networks favour low frequencies, and targets often are not low frequency
lang: en
translation: computational-mathematics/paper-notes/scientific-machine-learning/spectral-bias-and-generative-solvers
tags:
  - paper-notes
  - scientific-machine-learning
  - spectral-bias
---

> [!note] Coverage of this page
> Papers **81** (_Comput. Methods Appl. Mech. Engrg._ 437, 2025), **89** (_Neural Networks_ 194, 2026), **94** (_J. Comput. Phys._ 558, 2026), **101** (submitted to _J. Comput. Phys._, [arXiv:2512.18586](https://arxiv.org/abs/2512.18586)), **103** (submitted to _J. Comput. Phys._) and **105** (submitted to _Comput. Methods Appl. Mech. Engrg._).

## Spectral bias is a measurable quantity

Spectral bias, or the frequency principle, is the observation that networks converge rapidly to low-frequency components and struggle to represent high-frequency or highly oscillatory features. Multi-scale networks improve on that with a down-scaling map: decompose the compactly supported Fourier transform over concentric annuli $\mathbb K_i=\{\bm k:(i-1)K_0\le|\bm k|\le iK_0\}$, scale each band down by $a_i>1$, and approximate

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

## 101: keep the feature bank fixed and learn attention weights

Paper 101 attacks the same spectral-bias problem as paper 81 from a different angle. Its criticism of the existing fixes (random Fourier features, multi-scale networks and others) is that they freeze the choice of frequency content into the architecture. Its alternative keeps a **fixed** bank of multiscale random Fourier features and learns input-dependent **cross-attention** weights over it.

The contrast between the two routes is worth stating: paper 81 **rebuilds the network** around the captured dominant modes, which is discrete and requires retraining; paper 101 makes the weights **continuously learnable**, with no structural rebuild. The first has features aligned exactly with the target frequencies; the second needs no discrete rebuild step.

## 103 and 105: two other ways of writing structure into a network

- **103 (PI-DOSnet)** writes the structure of **operator splitting** into the network for evolution equations. Lie-Trotter and Strang splitting decompose a complex evolution into substeps that can be handled separately; writing that structure explicitly into a network means the network's layers correspond to the substeps rather than to an unstructured map.
- **105 (FLUID)** carries the random-field flow route of paper 62 into a unified inference framework for dynamics, replacing the Karhunen-Loève-structured reference field with conditional flows.

## 89: sampling from a known unnormalised energy

Paper 89 treats **sampling from a Boltzmann distribution**: given an unnormalised energy $E(x)$, the target is $\pi(x)\propto e^{-E(x)}$ with an unavailable normalising constant. This is a generative-model paper rather than a PDE solver, but it shares the stance of the [[computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|density-flow line]]: represent a distribution by an invertible or sampleable model, so "sampling" becomes "training a model".

Its approach is an energy-based diffusion generator. Diffusion or score-based methods differ from normalizing flows in a specific way: flows give an explicit density but are constrained by invertibility, while diffusion models are more expressive but supply the density only implicitly through a stochastic differential equation. When only the energy is known and no samples are available, the training signal has to come from the energy itself, so a reverse-KL or comparable variational objective is needed.

> [!note] Coverage status
> The constructions, error bounds and algorithm of papers 81 and 94 have been checked. The specific losses, architectures and theorems of papers 89, 101, 103 and 105 have not been checked equation by equation here; the content above is limited to what titles, abstracts and cross-references in neighbouring papers confirm.

## Where the six sit

| No. | How frequency content is handled                  | Structural rebuild |
| --- | ------------------------------------------------- | ------------------ |
| 81  | captured a posteriori from the solution's DFT     | yes (two criteria) |
| 94  | placed analytically inside a Gaussian wave packet | not applicable     |
| 101 | fixed feature bank plus learnable cross-attention | no                 |
| 103 | not applicable (operator-splitting structure)     | no                 |
| 105 | not applicable (conditional flow representation)  | no                 |
| 89  | not applicable (Boltzmann sampling)               | no                 |

One judgement runs through all of them: **frequency content is either measured (paper 81), removed analytically (paper 94), or weighted by the model itself over a sufficiently wide feature bank (paper 101). Assuming it is low-frequency is outside all three, and that assumption is precisely what spectral bias means.**

## Sources for this page

- J. Huang, R. You, and T. Zhou, [_Frequency-adaptive multi-scale deep neural networks_](https://doi.org/10.1016/j.cma.2025.117751), Comput. Methods Appl. Mech. Engrg. 437 (2025), 117751 (preprint [arXiv:2410.00053](https://arxiv.org/abs/2410.00053)).
- Y. Wang, L. Guo, H. Wu, and T. Zhou, [_Energy-based diffusion generator for efficient sampling of Boltzmann distributions_](https://doi.org/10.1016/j.neunet.2025.108126), Neural Networks 194 (2026), 108126.
- J. Huang, R. You, and T. Zhou, [_Deep learning for the semi-classical limit of the Schrödinger equation_](https://doi.org/10.1016/j.jcp.2026.114869), J. Comput. Phys. 558 (2026), 114869 (preprint [arXiv:2509.04453](https://arxiv.org/abs/2509.04453)).
- X. Feng, T. Tang, X. Wan, and T. Zhou, _Overcoming spectral bias via cross-attention_, [arXiv:2512.18586](https://arxiv.org/abs/2512.18586), submitted to J. Comput. Phys.
- J. Huang, Y. Qian, and T. Zhou, _PI-DOSnet: a physics-informed deep operator-splitting network for evolution partial differential equations_, submitted to J. Comput. Phys.
- T. Cui, X. Feng, C. Pei, X. Wan, and T. Zhou, _FLUID: flow-based unified inference for dynamics_, submitted to Comput. Methods Appl. Mech. Engrg.
