---
title: Densities as Invertible Maps - From Random Fields to Fractional Equations
description: Papers 62, 64, 72 and 87 - write an unknown density as the push-forward of a reference and train it on the equation residual
lang: en
translation: computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities
tags:
  - paper-notes
  - scientific-machine-learning
  - normalizing-flows
---

> [!note] Coverage of this page
> Papers **62** (_J. Comput. Phys._ 461, 2022), **64** (_Commun. Comput. Phys._ 32(2), 2022), **72** (_J. Sci. Comput._ 97:68, 2023) and **87** (_SIAM J. Sci. Comput._ 47(6), 2025).

![Write a density equation as a parametrised transport map](assets/diagrams/tao-zhou-papers/en/density-flow-solvers.svg)

The shared judgement across these four papers is that when the solution of an equation **is itself a probability density**, approximating it on a grid wastes structure. A density must be non-negative and integrate to one; on a grid those are constraints, but under an invertible map they are identities.

The push-forward formula

$$
p_{X}(x)=p_{Z}\bigl(f(x)\bigr)\,\bigl|\det\nabla_x f(x)\bigr|
$$

delivers non-negativity for free (a density times a positive determinant) and normalisation for free (change of variables). Only two questions remain: how to parameterise $f$ so that it is expressive and $\det\nabla_x f$ stays computable, and what loss trains it when no labelled data exist. The four papers push on the **target object**: random fields, time-dependent densities, fractional operators, bounded supports.

## 62: a random field as a flow with Karhunen-Loève structure

### The problem

In data-driven forward and inverse stochastic PDE problems, the random input field (diffusivity, forcing) or the solution is observed only through scattered sensors, and the covariance structure of the randomness is not known in advance. The classical route reduces the field by a Karhunen-Loève expansion and then builds a surrogate in parameter space, which requires prior knowledge of the covariance and suffers from the curse of dimensionality. The arbitrary-polynomial-chaos network methods the paper cites have term counts growing exponentially in the effective dimension, and physics-informed GANs do not give a tractable likelihood.

### A three-step construction

The paper states its construction in three steps:

1. Build a reference Gaussian random field $z(x,\omega)$ with a **truncated Karhunen-Loève expansion structure**, whose expansion coefficients are parameterised by deep networks.
2. Build a bijective transformation, a normalizing flow, between the reference field and the target stochastic field.
3. Train all parameters by maximising the sum of the log-likelihood over the scattered measurements.

For stochastic differential and partial differential equations, known physics enters through an added residual loss, giving the physics-informed version.

The value of this construction is that it settles two things at once: the Karhunen-Loève structure represents spatial correlation while the flow supplies non-Gaussianity, and because the flow is a bijection the likelihood is explicit, so training is maximum likelihood rather than adversarial. The advantages the paper lists are that one framework handles forward, inverse and mixed problems; sensor locations need not be fixed across snapshots; and the curse of dimensionality afflicting polynomial chaos is alleviated.

The experiments fall into three groups: learning stochastic processes, including non-Gaussian and mixed non-Gaussian fields; forward stochastic elliptic equations; and inverse stochastic elliptic equations.

> [!note] What could be verified
> The algebraic form of the Karhunen-Loève-structured reference field, the coupling layers of the flow, and the combined likelihood-plus-physics loss could not be confirmed equation by equation from the public material used here. Read Section 3 of the source directly if those formulas are needed.

## 64: time is not an extra dimension

### The key obstruction

The solution of a time-dependent Fokker-Planck equation is a density on an unbounded, possibly high-dimensional domain. The natural way to extend a flow to time is to treat time as an extra dimension, setting $\widehat x=(x,t)$ and $\widehat z=(z,t^\ast)$ and writing

$$
p_{\widehat X}(\widehat x)=p_{\widehat Z}(\widehat z)\,|\det J|,
\qquad \widehat z=f(\widehat x).
$$

That is wrong: mass is not conserved along the time axis, so $\int p(x,t)\,\mathrm dx\,\mathrm dt\neq1$. The paper's fix is to **pin the latent time to the real time**, $t^\ast=t$, which collapses the Jacobian to

$$
\det J=\begin{vmatrix}\partial z/\partial x & \partial z/\partial t\\ 0 & 1\end{vmatrix}
=\Bigl|\frac{\partial z(x,t)}{\partial x}\Bigr| ,
$$

giving

$$
p_{\widehat X}(x,t)=p_{\widehat Z}(z,t)\,\Bigl|\frac{\partial z}{\partial x}\Bigr|,
\qquad z=f(x,t).
$$

This deserves a sentence of its own: **the flow is time-conditioned, not time-augmented.** The one-word difference decides whether normalisation holds.

### The time-dependent affine coupling layer

The architecture is a simplified extension of KRnet from the spatial to the space-time domain. Each layer is an Actnorm layer

$$
y_{[i]}=a_i\odot x_{[i]}+b_i
$$

followed by a time-dependent affine coupling layer:

$$
x_{[i],1}=x_{[i-1],1},
$$

$$
x_{[i],2}=x_{[i-1],2}\odot\Bigl(\mathbf 1_{d-m}
+\beta\tanh\bigl(s_i(x_{[i-1],1},t)\bigr)\Bigr)
+e^{\zeta_i}\odot\tanh\bigl(q_i(x_{[i-1],1},t)\bigr),
$$

where $|\beta|<1$ is user-specified (for example $0.6$), $s_i,q_i:\mathbb R^{m+1}\to\mathbb R^{d-m}$, and $\zeta_i\in\mathbb R^{d-m}$ is trainable. A final polynomial spline transformation increases modelling power.

Two things differ from standard real NVP: the $\tanh$ confines the scale factor to $(1-\beta,1+\beta)$ and prevents scale blow-up, and $s_i,q_i$ take $t$ as an additional input so one set of weights covers the whole time interval.

The loss is the time-dependent Fokker-Planck residual, with no labelled data and no simulated sample paths of the underlying stochastic differential equation. Experiments cover linear drift, nonlinear drift and high-dimensional problems.

## 72: two ways to handle a fractional operator

Paper 72 pushes this line to fractional Fokker-Planck equations. The stochastic differential equation carries Lévy noise:

$$
\mathrm d X_t=\mu(X_t,t)\,\mathrm dt+\sigma(X_t,t)\,\mathrm dW_t+\mathrm dL^\alpha_t,
\qquad
\frac{\partial p}{\partial t}=\mathcal L p-(-\Delta)^{\alpha/2}p,
$$

$$
\mathcal L p=-\nabla\cdot(p\mu)+\tfrac12\nabla\cdot\nabla\cdot(\sigma\sigma^{\mathsf T}p).
$$

Three difficulties appear together: an unbounded domain, dimension, and nonlocality. The flow handles the first two; the third has two routes.

### Route one: Monte Carlo (MCNF)

Import the inner-outer split estimator of [[en/computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|paper 66]] as this paper's Lemma 3.1, with loss

$$
L(p_{\mathrm{KRnet},\theta})=\frac{1}{N_S}\sum_{i=1}^{N_S}\bigl|R_\theta(x^i)\bigr|^2,
\qquad
R_\theta(x)=\bigl(\mathcal L-(-\Delta)^{\alpha/2}\bigr)p_{\mathrm{KRnet},\theta}(x).
$$

One easily missed but practically important observation: $\mathrm{Beta}(a,1)$ concentrates at the origin as $a\to0$, so the closer $\alpha$ gets to $2$, the more the inner radii $r_1$ pile up near zero and the **larger** the floor $r_\epsilon$ must be for stability. That turns "parameter values affect numerical stability" into an actionable rule.

> [!note] A small difference from paper 66
> Paper 72 uses two **independent** direction variables $\xi$ (inner) and $\eta$ (outer), whereas the corresponding formula in paper 66 reuses a single $\xi$. The difference is small but real.

### Route two: an analytic auxiliary model (GRBFNF)

The second route avoids randomness by switching to a model whose fractional Laplacian is known in closed form. For a Gaussian $u(x)=\exp(-\sigma^{-2}|x-x_0|_2^2)$,

$$
(-\Delta)^{\alpha/2}u(x)=c_{\alpha,d}\,|\sigma|^{-\alpha}\;
{}_1F_1\!\Bigl(\frac{d+\alpha}{2};\frac{d}{2};-\sigma^{-2}|x-x_0|^2_2\Bigr),
\qquad
c_{\alpha,d}=\frac{2^\alpha\Gamma\bigl(\frac{d+\alpha}{2}\bigr)}{\Gamma\bigl(\frac d2\bigr)},
$$

with ${}_1F_1$ the confluent hypergeometric function. A Gaussian radial basis mixture then serves as auxiliary model,

$$
p_{\mathrm{GRBF},\tilde\theta}(x)=\sum_{i=1}^M w_i\,\mathcal N(\tilde x_i,\sigma_i^2\mathbf I)(x),
\qquad 0\le w_i\le1,\ \sum_i w_i=1,
$$

with both $w_i$ and $\sigma_i$ trainable, whose fractional Laplacian is explicit:

$$
(-\Delta)^{\alpha/2}p_{\mathrm{GRBF},\tilde\theta}(x)
=c_{\alpha,d}\,\pi^{-d/2}2^{-\frac{d+\alpha}{2}}\sum_{i=1}^M w_i\,|\sigma_i|^{-(d+\alpha)}\,
{}_1F_1\!\Bigl(\frac{d+\alpha}{2};\frac d2;-\frac{|x-\tilde x_i|_2^2}{2\sigma_i^2}\Bigr).
$$

A coupled loss ties the flow to the auxiliary model, with the second term a consistency penalty:

$$
\tilde L\bigl(p_{\mathrm{KRnet},\theta},p_{\mathrm{GRBF},\tilde\theta}\bigr)
=\frac{1}{N_S}\sum_{i}\Bigl(\mathcal L p_{\mathrm{KRnet},\theta}(x^i)
-(-\Delta)^{\alpha/2}p_{\mathrm{GRBF},\tilde\theta}(x^i)\Bigr)^2
+\frac{\beta_m}{N_S}\sum_{i}\Bigl(p_{\mathrm{KRnet},\theta}(x^i)
-p_{\mathrm{GRBF},\tilde\theta}(x^i)\Bigr)^2 .
$$

This is the most distinctive idea in the paper: **move the hard operator onto a model where it is analytic, then tie the two representations together with a consistency penalty.** The price is that the number of radial basis centres grows with dimension, so this route suits low dimension while the Monte Carlo route suits high dimension.

### One detail in the time-dependent version

The time-dependent model (MCTNF) uses the coupling layer

$$
x_{[i],2}=x_{[i-1],2}\odot\Bigl(\mathbf 1_{d-m}
+\beta\tanh\bigl(t\,s_{i,t}(x_{[i-1],1},t)\bigr)\Bigr)
+e^{\zeta_i}\odot\tanh\bigl(t\,q_{i,t}(x_{[i-1],1},t)\bigr).
$$

Note the explicit factor $t$ multiplying the network outputs. It makes the layer the identity at $t=0$, so the **initial condition holds by construction** rather than being approximated by a penalty. Writing constraints into the architecture is the same move as letting the push-forward guarantee normalisation.

Training alternates refinement: draw a new training set from the current flow, then keep training on it. The experiments include a two-dimensional equation driven by the fractional Laplacian alone, a bimodal target, higher-dimensional stationary equations, and a time-dependent problem whose solution is a Cauchy distribution, the natural stationary law for $\alpha$-stable noise.

> [!note] Author order
> The homepage lists the authors as Xiaoliang Wan, Li Zeng, Tao Zhou; both the preprint and the journal version give **Li Zeng, Xiaoliang Wan, Tao Zhou**. This site records the published order.

## 87: a bounded support needs a different coupling layer

### Why truncation is not enough

KRnet and normalizing flows in general map to a Gaussian reference on $\mathbb R^d$, so the induced density has unbounded support. That is wrong for densities genuinely supported on a bounded set (the paper's everyday example is human ages) and wrong for PDE solutions posed on a hyperrectangle. Forcibly restricting a flow to a box destroys exact invertibility.

### The CDF coupling layer

The paper's core new object is a monotone map taking $[-1,1]$ exactly onto $[-1,1]$. Define a piecewise-linear density on a mesh $-1=s_0<s_1<\cdots<s_n=1$,

$$
p(s)=\frac{w_{i+1}-w_i}{h_i}(s-s_i)+w_i,
\quad s\in[s_i,s_{i+1}],
\quad p(s_i)=w_i\ge0,\ h_i=s_{i+1}-s_i,
$$

whose cumulative distribution function is **piecewise quadratic**:

$$
F(s)=\frac{w_{i+1}-w_i}{2h_i}(s-s_i)^2+w_i(s-s_i)
+\sum_{k=0}^{i-1}\frac{w_k+w_{k+1}}{2}(s_{k+1}-s_k).
$$

Its inverse is the root of a quadratic. With $q_0=0$ and $q_i=\sum_{k=0}^{i-1}\frac{w_k+w_{k+1}}{2}(s_{k+1}-s_k)$,

$$
F^{-1}(q)=s_i+\frac{-w_i+\sqrt{w_i^2+2(w_{i+1}-w_i)(q-q_i)/h_i}}{(w_{i+1}-w_i)/h_i}
=s_i+\frac{2(q-q_i)}{w_i+\sqrt{w_i^2+2(w_{i+1}-w_i)(q-q_i)/h_i}} .
$$

The second form is the numerically stable one when the denominator is small, an implementation detail that matters. Setting $\tilde F(s)=2F(s)-1$ gives a bijection of $[-1,1]$ onto itself, applied componentwise. For $y=(y_1^{\mathsf T},y_2^{\mathsf T})^{\mathsf T}\in[-1,1]^l$ with $y_1\in\mathbb R^m$, the coupling layer is

$$
\hat y_1=y_1,
\qquad
\hat y_2=F\bigl(y_2;\theta(y_1)\bigr),
$$

where $\theta=(s_1^{\mathsf T},\dots,s_{n-1}^{\mathsf T},w_0^{\mathsf T},\dots,w_n^{\mathsf T})^{\mathsf T}\in\mathbb R^{2n(l-m)}$ is produced by a network $\mathrm{NN}(y_1)$ and reparameterised so the mesh increases and the weights stay non-negative:

$$
s_1=-1+\frac2n\bigl(1+\beta_1\tanh(\hat s_1)\bigr),
\qquad
s_{i+1}=s_i+\Bigl(\frac{1-s_i}{n-i}\Bigr)\bigl(1+\beta_{i+1}\tanh(\hat s_{i+1})\bigr),
\qquad
w_i=\frac{1+\gamma_i\tanh(\hat w_i)}{C},
$$

with $C$ a normalising constant and, in practice, $\beta_1=\tfrac{65}{66}$, $\beta_i=0.97$ for $i\ge2$, and $\gamma_i=0.99$. Those constants keep the transformation strictly inside the monotone regime and prevent degeneracy. A second layer with the roles of $y_1$ and $y_2$ swapped completes the update, and layers alternate by parity.

### Four differences from KRnet

| Item             | KRnet                     | B-KRnet                                    |
| ---------------- | ------------------------- | ------------------------------------------ |
| Reference        | Gaussian on $\mathbb R^d$ | uniform on $[-1,1]^d$                      |
| Coupling layer   | affine (scale and shift)  | monotone piecewise-quadratic CDF           |
| Scale-bias layer | present                   | **removed** (unnecessary on a bounded box) |
| Smoothness       | set by the activation     | density piecewise linear, only $C^1$       |

The last row has an immediate consequence: **a second-order PDE must be recast as a first-order system.** That is what the paper does, introducing $g(x)=\nabla p$ approximated by a separate network $g_{\mathrm{NN}}$ and minimising

$$
\mathcal L=\lambda_{pde}\mathcal L_{pde}+\lambda_b\mathcal L_b+\lambda_{g}\mathcal L_{g},
$$

$$
\mathcal L_{pde}=\mathbb E_{x\sim\rho}\bigl[|\mathcal N[x;p_{\text{B-KRnet},\theta},g_{\mathrm{NN}}]|^2\bigr],
\qquad
\mathcal L_{g}=\mathbb E_{x\sim\rho}\bigl[\|g_{\mathrm{NN}}(x)-\nabla p_{\text{B-KRnet},\theta}(x)\|_2^2\bigr].
$$

Non-negativity and mass conservation **hold by construction** and never enter as penalties. A general product domain $\prod_i[a_i,b_i]$ is mapped to $[-1,1]^d$ by the linear map $y=\hat a\odot x+\hat b$ with $\hat a_i=2/(b_i-a_i)$ and $\hat b_i=-(b_i+a_i)/(b_i-a_i)$.

### Pseudo-triangular structure and outer freezing

Partition $x=((x^{(1)})^{\mathsf T},\dots,(x^{(K)})^{\mathsf T})^{\mathsf T}$ with $\sum_i d_i=d$; the triangular structure of the Knothe-Rosenblatt rearrangement relaxes to a block-triangular one,

$$
z=f_{\text{KR}}(x)=
\begin{pmatrix}\tilde f_1\\ x^{(2:K)}\end{pmatrix}\circ
\begin{pmatrix}\tilde f_2\\ x^{(3:K)}\end{pmatrix}\circ\cdots\circ
\begin{pmatrix}\tilde f_{K-1}\\ x^{(K)}\end{pmatrix}\circ \tilde f_K(x),
$$

where the outer loop freezes one block per stage (the paper calls it a squeezing operation that deactivates dimensions) and the inner loop composes several CDF coupling layers. One boundary case deserves note: a CDF coupling layer needs input dimension at least $2$, so when $x^{(1)}=x_1\in\mathbb R$ the outermost map is taken to be the identity, justified because $x^{(1)}$ is already close to uniform by then, making $2F(\cdot)-1$ close to the identity.

### Mixed domains

The solution of a kinetic Fokker-Planck equation is a joint density of position $x\in\Omega$ (bounded) and velocity $v\in\mathbb R^d$ (unbounded). The paper factors $p(x,v)=h(v|x)q(x)$ and models

$$
p_{\theta}(x,v)=h_{\text{KRnet},\theta_1}(v|x)\cdot q_{\text{B-KRnet},\theta_2}(x),
$$

the unbounded KRnet as a conditional flow for velocity times B-KRnet for position. The conditioning mechanism is exactly the time conditioning of paper 64 with position in place of time: the same technical component serving a different role in a different paper.

### The adaptive strategy

$$
\rho(x)\leftarrow(1-\gamma)\rho(x)+\gamma\,p_{\text{B-KRnet},\theta^{*,k}}(x),
\qquad N_{new}=\gamma N_{pde},
$$

with new points generated through the inverse flow, $z^{i}\sim\mathrm{Unif}(-1,1)^d$ and $x^{i}=(f_{\theta^{*,k}})^{-1}(z^{i})$, keeping the rest. The implementation uses update rate $\gamma=0.8$, three subintervals per CDF coupling layer, `Tanh` activations, and Adam with a step-decay schedule.

Density estimation is tested on an annulus, a Gaussian mixture and a logistic distribution with holes; PDE approximation on a four-dimensional $-\Delta p+p=f$, the stationary Keller-Segel system (two coupled densities) and a two-dimensional stationary kinetic Fokker-Planck equation. The paper's own stated limitations are that B-KRnet is smooth only to first order, so improving it would need higher-order polynomials, and that adaptive point updating inside a deep Ritz framework is left open — which is what [[en/computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|paper 80]] takes up.

> [!note] Version drift
> The v1 abstract says B-KRnet "consists of a series of coupling layers with progressively fewer active transformation dimensions, inspired by the triangular structure of the Knothe-Rosenblatt rearrangement" and calls the reference measure a base distribution; v3 says the structure "adapts the pseudo-triangular structure into a normalizing flow model" and calls it a prior distribution. Cite the version you quote.

## How the four papers build on each other

| No. | Target object                    | Reference                  | Training signal                  | New component                           |
| --- | -------------------------------- | -------------------------- | -------------------------------- | --------------------------------------- |
| 62  | random field from scattered data | Gaussian with KL structure | log-likelihood (plus physics)    | KL expansion with network coefficients  |
| 64  | time-dependent Fokker-Planck     | unbounded Gaussian         | equation residual                | latent time pinned to real time         |
| 72  | fractional Fokker-Planck         | unbounded Gaussian         | equation residual                | Monte Carlo or analytic auxiliary model |
| 87  | densities with bounded support   | uniform on $[-1,1]^d$      | residual of a first-order system | piecewise-quadratic CDF coupling layer  |

One judgement runs through all four: **structural constraints belong in the architecture, not in a penalty.** Normalisation comes from the change of variables, positivity from the push-forward, the initial condition from the factor $t$ in the coupling layer, and bounded support from the CDF coupling layer. Each such guarantee removes one penalty term whose weight would otherwise need tuning.

## Coverage check

| Item                                                 | Paper | Status                                                           |
| ---------------------------------------------------- | ----- | ---------------------------------------------------------------- |
| KL-structured reference field and three-step build   | 62    | construction steps and advantages, with verification limits      |
| Time conditioning and the collapsed Jacobian         | 64    | the wrong route, the right route, and the consequence            |
| Time-dependent affine coupling layer                 | 64    | formula, role of $\beta$, differences from real NVP              |
| Fractional Fokker-Planck under Lévy noise            | 72    | the equation and both operator routes                            |
| Beta concentration and the $r_\epsilon$ floor        | 72    | the parameter-dependent stability rule                           |
| Analytic fractional Laplacian of a Gaussian and GRBF | 72    | closed form, mixture model, coupled loss, consistency term       |
| The factor $t$ enforcing the initial condition       | 72    | coupling-layer form and its meaning                              |
| Piecewise-linear density and piecewise-quadratic CDF | 87    | density, CDF, stable inverse, reparameterisation, constants      |
| Four differences from KRnet                          | 87    | reference, coupling layer, removed scale-bias, $C^1$ limit       |
| First-order recasting and the three-term loss        | 87    | introducing $g$, three terms, the two built-in properties        |
| Pseudo-triangular structure, edge case, mixed domain | 87    | block factorisation, freezing, $d_1=1$ case, conditional product |

## Sources for this page

- L. Guo, H. Wu, and T. Zhou, _Normalizing field flows: solving forward and inverse stochastic differential equations using physics-informed flow models_, J. Comput. Phys. 461 (2022), 111202 (preprint [arXiv:2108.12956](https://arxiv.org/abs/2108.12956)).
- X. Feng, L. Zeng, and T. Zhou, [_Solving time dependent Fokker-Planck equations via temporal normalizing flow_](https://doi.org/10.4208/cicp.OA-2022-0090), Commun. Comput. Phys. 32(2) (2022), pp. 401-423 (preprint [arXiv:2112.14012](https://arxiv.org/abs/2112.14012)).
- L. Zeng, X. Wan, and T. Zhou, [_Adaptive deep density approximation for fractional Fokker-Planck equations_](https://doi.org/10.1007/s10915-023-02379-z), J. Sci. Comput. 97 (2023), 68 (preprint [arXiv:2210.14402](https://arxiv.org/abs/2210.14402)).
- L. Zeng, X. Wan, and T. Zhou, _Bounded KRnet and its applications to density estimation and approximation_, SIAM J. Sci. Comput. 47(6) (2025), pp. C1294-C1318 (preprint [arXiv:2305.09063](https://arxiv.org/abs/2305.09063)).
