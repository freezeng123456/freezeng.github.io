---
title: Multi-Fidelity Surrogates Refined Online
description: Papers 34, 37, 49 and 79 - one closed loop realised with four surrogates and four samplers
lang: en
translation: computational-mathematics/paper-notes/bayesian-inference/multifidelity-surrogates
tags:
  - paper-notes
  - bayesian-inverse-problems
  - surrogate-modelling
---

> [!note] Coverage of this page
> Papers **37** (_J. Comput. Phys._ 381, 2019), **34** (_Int. J. Uncertain. Quantif._ 9(3), 2019), **49** (_Commun. Comput. Phys._ 28, 2020) and **79** (_SIAM/ASA J. Uncertain. Quantif._ 12(4), 2024). All four share one closed-loop skeleton and replace one component at a time.

## 37: the prototype of the loop

### What the paper targets

A prior-based polynomial chaos surrogate makes MCMC cheap, but it is fitted where the prior lives while the posterior concentrates on a thin subset of the prior support. Raising the polynomial order removes the error at a cost that grows combinatorially in the parameter dimension. The paper's choice is to **keep the order low and correct online**.

### The exact form of the multi-fidelity correction

Let the low-fidelity model be the prior-based polynomial chaos expansion $u^L$ and let the high-fidelity model $u^H$ be the true forward solve. The correction

$$
C(z)=u^{H}(z)-u^{L}(z)\approx\sum_{\mathbf m\in\Lambda_{N_C}}u^{C}_{\mathbf m}\Phi_{\mathbf m}(z)
$$

is expanded on a **low-order** index set $\Lambda_{N_C}$ with $N_C\ll N$. Merging adds the two coefficient sets index by index:

$$
f_{M}(\theta)=\sum_{\alpha\in\Lambda^{d}_{N_C}}\bigl(u^{L}_{\alpha}+u^{C}_{\alpha}\bigr)\Psi_{\alpha}
+\sum_{\alpha\in\Lambda^{d}_{N}\setminus\Lambda^{d}_{N_C}}u^{L}_{\alpha}\Psi_{\alpha}.
$$

Written this way, the structural property of the correction becomes visible: **only the low-order coefficients are corrected by high-fidelity data**, while the high-order coefficients keep their prior-fitted values. This is deliberate. There are not enough samples to identify the high-order indices, and correcting them would overfit.

The coefficients come from weighted discrete least squares:

$$
\{c_m\}_{m=1}^{M}=\arg\min_{c}\sum_{i=1}^{Q}\Bigl[G(z_i)-\sum_{m=1}^{M}c_m\Phi_m(z_i)\Bigr]^{2},
\qquad
w_i=\frac{M}{\sum_{m=1}^{M}\Phi_m^{2}(z_i)} .
$$

The weight $w_i$ is the discrete Christoffel function, the same object that governs optimal sampling in the [[en/computational-mathematics/paper-notes/stochastic-approximation/index|stochastic approximation topic]]. The sampling density follows degree-asymptotic designs: tensor-product Chebyshev for a uniform measure, and for a Gaussian measure a density supported on the ball of radius $\sqrt{2N}$.

### Trigger and local refinement

Refinement is decided at an accepted point, using an **absolute** $\ell^\infty$ error:

$$
\mathrm{err}(y)=\bigl\|u^{H}(y)-u^{L}(y)\bigr\|_{\infty}.
$$

Here $y$ is not an arbitrary chain state but a point accepted with the **high-fidelity acceptance probability**. That detail carries the design: accepting a point with the true model costs one true solve and buys a test point that sits closer to the posterior bulk.

If $\mathrm{err}(y)>\epsilon$, the method draws $Q$ random points in the $\ell^\infty$ ball $B(y,R)=\{x:\|x-y\|_\infty\le R\}$ and refits, with $Q=2\binom{N_C+n_z}{n_z}$ in the total-degree space. The radius shrinks each round by an input constant $\rho$, because early accepted points may still be far from the posterior bulk.

### Theorem 2: why local accuracy suffices

Define the $\epsilon$-feasible set and its posterior measure

$$
\Gamma_{N}(\epsilon)=\bigl\{y\in\Gamma:\ \|u^{H}(y)-u^{L}(y)\|_{\infty}\le\epsilon\bigr\},
\qquad
\mu\bigl(\Gamma_{N}(\epsilon)\bigr)=\int_{\Gamma_{N}(\epsilon)}\pi^{d}(z)\,dz .
$$

Under Assumption 1 ($\sup_{z\in\Gamma}\|u^H(z)\|=:C_H<\infty$) and i.i.d. Gaussian observational noise, there are constants $K_1,K_2>0$ with

$$
D_{\mathrm{KL}}\bigl(\widetilde\pi^{d}_{N}\,\|\,\pi^{d}\bigr)
\le\Bigl(K_1\epsilon+K_2\,\mu\bigl(\Gamma^{\perp}_{N}(\epsilon)\bigr)\Bigr)^{2},
$$

and a Hellinger bound $K_1\epsilon+K_2\mu(\Gamma^{\perp}_{N}(\epsilon))$. The corollary states the algorithmic goal: if sampling is good enough that $\mu(\Gamma^{\perp}_{N}(\epsilon))\le\epsilon$, the KL distance is characterised entirely by $\epsilon^2$.

The paper is careful to present the link between algorithm and theorem as a **mechanism**, not a theorem: whenever a candidate falls in $\Gamma^{\perp}_{N}(\epsilon)$, the algorithm refines near it, so $\mu(\Gamma^{\perp}_{N}(\epsilon))$ decays asymptotically with refinement. That is an argument, not a rate.

### Numerical evidence

Two nonlinear PDE inverse problems. The first is a two-dimensional heat source inversion with $n_z=2$, small enough that a genuinely accurate high-order prior surrogate is affordable; here the adaptive method buys accuracy rather than time. The second infers the diffusion coefficient of an elliptic PDE with nine random parameters, where a globally accurate prior surrogate is expensive and the adaptive method improves both accuracy and cost. The abstract's order-of-magnitude claim is several orders of efficiency gain over MCMC on the true model alone.

## 34: the same correction inside ensemble Kalman inversion

### Replacing the sampler

Ensemble Kalman inversion is derivative-free, which suits PDE-constrained parameter estimation, but it is a Monte Carlo method: $N_e$ members over $J$ iterations means roughly $N_e J$ forward solves. This paper puts the multi-fidelity correction of 37 inside a regularising iterative ensemble Kalman smoother. The update reads

$$
\theta^{(j)}_{n+1}=\theta^{(j)}_{n}
+C^{\theta\omega}_{n}\bigl(C^{\omega\omega}_{n}+\alpha_{n}\Gamma\bigr)^{-1}\bigl(y^{(j)}-\omega^{(j)}_{n}\bigr),
\qquad j=1,\dots,N_e,
$$

with $\omega^{(j)}_n=f(\theta^{(j)}_n)$ and perturbed data $y^{(j)}=y+\xi^{(j)}$, and empirical covariances

$$
C^{\theta\omega}_{n}=\frac{1}{N_e-1}\sum_{j}\bigl(\theta^{(j)}_{n}-\bar\theta_{n}\bigr)\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)^{\!T},
\qquad
C^{\omega\omega}_{n}=\frac{1}{N_e-1}\sum_{j}\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)\bigl(\omega^{(j)}_{n}-\bar\omega_{n}\bigr)^{\!T}.
$$

The regularisation parameter $\alpha_n$ is the value attached to the first integer $N$ satisfying

$$
\alpha^{N}_{n}\bigl\|\Gamma^{1/2}\bigl(C^{\omega\omega}_{n}+\alpha^{N}_{n}\Gamma\bigr)^{-1}\bigl(y^{(j)}-\bar\omega_{n}\bigr)\bigr\|
\ \ge\ \rho\,\bigl\|\Gamma^{-1/2}\bigl(y^{(j)}-\bar\omega_{n}\bigr)\bigr\| ,
$$

and the iteration stops by the discrepancy principle $\|\Gamma^{-1/2}(y-\bar\omega_{n})\|\le\tau\eta$ with $\rho<1$ and $\tau\ge1/\rho$. Both rules come from Iglesias' regularising ensemble Kalman smoother.

### A relative indicator, tested at the ensemble mean

$$
\mathrm{err}=\frac{\bigl\|f(\bar\theta_{n+1})-f_{M}(\bar\theta_{n+1})\bigr\|_{\infty}}{\bigl\|f(\bar\theta_{n+1})\bigr\|_{\infty}},
\qquad
\bar\theta_{n+1}=\frac{1}{N_e}\sum_{j}\theta^{(j)}_{n+1}.
$$

The ensemble mean plays the role that the accepted state plays in 37: it is where the current ensemble believes the truth lies. The sample budget is explicit: initialisation uses $Q_1=2\binom{N+d}{d}$ prior samples and each refinement uses $Q_2=2\binom{N_C+d}{d}$ new points.

### The most informative part of the experiments

The truth is deliberately placed outside the prior: the prior is $\log\theta_i\sim\mathcal N(0,1)$ while the truth is drawn from $\log\theta_i\sim U(-4,4)$. On a two-dimensional time-fractional inverse diffusion problem (Caputo derivative, nine radial-basis weights, $N_e=100$, noise level $\sigma=10^{-3}$), the fixed prior-based polynomial chaos surrogate gives relative error $0.7921$ at order $N=4$ and $0.2892$ at $N=6$, against $0.0461$ for conventional ensemble Kalman inversion. On cost, the conventional method takes about $56.71$ seconds and the polynomial chaos version about $0.82$ seconds. The adaptive method at order $N=2$ needs only 250 true solves (tolerance $10^{-2}$) or 575 (tolerance $10^{-3}$), against 2000 for the conventional method.

The second example recovers a log-permeability field parameterised by a Karhunen-Loève expansion, retaining 95% of the prior energy in $d=22$ modes with $N_e=300$. The fixed surrogate gives relative error $0.3430$ at $N=2$ and $0.2146$ at $N=3$, while the adaptive method reaches $0.0889$ even at $N=1$.

These numbers make one point sharply: **when the truth is out of prior, raising the order of a global surrogate is an inefficient repair**, while a local correction built on a very low-order expansion is enough.

## 49: replacing polynomials with a composite network

### Why replace them

The paper names two weaknesses of polynomial surrogates explicitly: poor approximation of maps with limited regularity, and a basis count that explodes with parameter dimension. Networks are better suited on both counts.

### The composite multi-fidelity network

The key construction is not "train another network" but feeding the already-trained low-fidelity network into the high-fidelity network **as an input variable**:

$$
f^{H}(z)=\mathcal F\bigl(z,f^{L}(z)\bigr)=\mathcal F\bigl(z,\mathcal{NN}^{L}(z);\theta\bigr),
\qquad
f^{H}(z)\approx\mathcal{NN}^{H}(z;\theta):=\mathcal{NN}\bigl(z,\mathcal{NN}^{L}(z);\theta\bigr).
$$

The paper draws the contrast with 37 itself: the polynomial chaos correction $f^{H}\approx f^{L}_{\mathrm{PCE}}+f_{\mathrm{CORR}}$ is a **linear** superposition, while the display above learns a **nonlinear** correlation between fidelities. Because the two models are highly correlated, $\mathcal{NN}^{H}$ can be shallow; the experiments use one hidden layer with 50 neurons.

The shallow network is not a performance optimisation but an **anti-overfitting constraint**: only $Q$ true solves are affordable per round ($Q=10$ in the experiments), and the network capacity has to match that.

### The indicator

$$
\mathrm{err}(\tilde z)=\frac{\bigl\|f^{H}(\tilde z)-f^{L}(\tilde z)\bigr\|_{\infty}}{\bigl\|f^{H}(\tilde z)\bigr\|_{\infty}} .
$$

The absolute error of 37 becomes relative; the test point $\tilde z$ is still produced by the high-fidelity acceptance probability. Once triggered, the method draws $Q$ uniform points in $B(\tilde z,R)$, retrains $\mathcal{NN}^{H}$, and sets $f^{L}\leftarrow\mathcal{NN}^{H}$. Note the consequence: the next round's "low-fidelity model" is the previous composite network, so the composition nests round by round.

> [!warning] Checking the source: numerator and denominator of the acceptance probability
> Equation (10) of the paper prints the high-fidelity acceptance probability as
> $\beta=\min\{1,\ \mathcal L(d,f^{H}(z^{-}))\pi(z^{-})/\mathcal L(d,f^{H}(z^{+}))\pi(z^{+})\}$.
> In Metropolis-Hastings with a symmetric proposal the proposed state $z^{+}$ belongs in the numerator, so in context this should read
> $\beta=\min\{1,\ \mathcal L(d,f^{H}(z^{+}))\pi(z^{+})/\mathcal L(d,f^{H}(z^{-}))\pi(z^{-})\}$.
> This is a reader's reconciliation; check the journal version before quoting.

### Numerical evidence

A benchmark elliptic inverse problem in two configurations. The first is a nine-parameter permeability inference compared against three baselines: MCMC on the true model, a fixed prior-trained network, and the adaptive version of this paper. The prior network uses 50 training points with 4 hidden layers of 40 neurons; the correction network has one hidden layer of 50 neurons; $Q=10$; tolerances $0.1$ and $0.05$. The second example uses a Karhunen-Loève parameterised high-dimensional random field. The fixed prior-trained network gives visibly wrong posterior marginals, the adaptive version recovers the true-model MCMC marginals, and tightening the tolerance from $0.1$ to $0.05$ tightens the agreement further.

The paper proves no theorem and states explicitly that it inherits the analytical setting of 37.

## 79: a goal-oriented indicator and a greedy design

### Two substantive changes

Paper 79 replaces the surrogate with an operator network (DeepONet) and the sampler with unscented Kalman inversion, and — methodologically the more important step — replaces the question "is the surrogate accurate?" with "has the fit been damaged?".

The paper first writes the honest local model error

$$
e_{M}(t):=\mathbb E_{\nu_t}\bigl\|\mathcal G-\widehat{\mathcal G}_t\bigr\|
=\Bigl(\int_{\mathcal M}\bigl|\mathcal G(m)-\widehat{\mathcal G}_t(m)\bigr|^{2}\,\nu_t(dm)\Bigr)^{1/2},
$$

then rejects it explicitly: it needs a high-dimensional integral. The replacement has two steps. Among $T$ samples $\mathcal M^{(t)}=\{m^{(t)}_k\}$ from the current surrogate-induced posterior $\nu_t$, use the **true model** to pick the anchor with the best data fit,

$$
r_t=\arg\min_{m\in\mathcal M^{(t)}}\tfrac12\bigl\|y-\mathcal G(m)\bigr\|^{2}_{\Sigma_\eta},
$$

then define the indicator as the data misfit at the anchor and trigger on its **relative change**:

$$
e_{D}(t):=\Phi(r_t;y)=\tfrac12\bigl\|y-\mathcal G(r_t)\bigr\|^{2}_{\Sigma_\eta},
\qquad
\frac{\bigl|e_D(t)-e_D(t-1)\bigr|}{e_D(t)}>\epsilon .
$$

The experiments use $\epsilon=0.01$ and cap retraining at $I_{\max}=10$.

### The greedy design: diversity against locality

Draw a large candidate pool $\Gamma=\{m_1,\dots,m_K\}$ from $\nu_t$ and grow a subset $\gamma_Q$ one point at a time:

$$
\hat m_{j+1}=\arg\max_{m\in\Gamma\setminus\gamma_j}
\Bigl\{d\bigl(\widehat{\mathcal G}_t(m),\widehat{\mathcal G}^{j}_t\bigr)-\lambda\|m-r_t\|_{2}\Bigr\},
\qquad
d\bigl(\widehat{\mathcal G}_t(\cdot),\widehat{\mathcal G}^{j}_t\bigr)
=\max_{\hat m\in\gamma_j}\bigl\|\widehat{\mathcal G}_t(\cdot)-\widehat{\mathcal G}_t(\hat m)\bigr\|_{2},
$$

with $\lambda=1$ in every experiment. The two terms are deliberately adversarial: the first demands separation in the **surrogate's output space**, which serves generalisation, and the second pulls selections back toward the anchor, which serves local accuracy. The score uses only surrogate predictions on the pool, so it costs essentially nothing. Retraining is transfer learning from the previous weights rather than from scratch.

### Unscented Kalman inversion

The posterior approximation comes from a stochastic dynamical system

$$
m_{n+1}=r_0+\alpha(m_n-r_0)+\omega_{n+1},\quad \omega_{n+1}\sim\mathcal N(0,\Sigma_\omega),
\qquad
y_{n+1}=\mathcal G(m_{n+1})+\eta_{n+1},
$$

with regularisation parameter $\alpha\in(0,1]$. The prediction step is $\hat r_{n+1}=\alpha r_n+(1-\alpha)r_0$ and $\hat C_{n+1}=\alpha^{2}C_n+\Sigma_\omega$, and the analysis step is

$$
r_{n+1}=\hat r_{n+1}+\hat C^{my}_{n+1}\bigl(\hat C^{yy}_{n+1}\bigr)^{-1}\bigl(y_{n+1}-\hat y_{n+1}\bigr),
\qquad
C_{n+1}=\hat C_{n+1}-\hat C^{my}_{n+1}\bigl(\hat C^{yy}_{n+1}\bigr)^{-1}\bigl(\hat C^{my}_{n+1}\bigr)^{T},
$$

with $\hat C^{yy}_{n+1}=\mathrm{Cov}[\mathcal G(m_{n+1})\mid Y_n]+\Sigma_\eta$. Expectations are evaluated by a modified unscented transform on $2N_m+1$ symmetric $\sigma$-points $m^{0}=r$, $m^{j}=r+c_j[\sqrt C]_j$, $m^{j+N_m}=r-c_j[\sqrt C]_j$. Each iteration costs $2N_m+1$ forward evaluations, and the paper reports convergence typically in $O(10)$ iterations.

### Theorem 3.6: the only convergence result in this family

Under three assumptions — (3.2) for any $\epsilon$ a **linear** neural operator $\widehat G$ can be trained with $\|\widehat G-G\|_2<\epsilon$; (3.3) $\|G\|_2<H$; (3.4) $G^{T}\Sigma_\eta^{-1}G\succ0$ with $\|G^{T}\Sigma_\eta^{-1}G\|_2>C_1$ — together with $\mathrm{Range}(G^{T})=\mathrm{Range}(\widehat G^{T})=\mathbb R^{N_m}$, $\Sigma_\omega\succ0$ and $\Sigma_\eta\succ0$, the surrogate-driven unscented Kalman inversion fixed point converges to the full-model fixed point, with

$$
\bigl\|\widehat C^{-1}_\infty-C^{-1}_\infty\bigr\|_2\le\frac{2\epsilon H H_\eta}{1-\beta},
\qquad
\bigl\|\widehat r_\infty-r_\infty\bigr\|_2
\le\frac{K_1H_\eta H_y}{C_1}\Bigl(1+\frac{2(1+\alpha\beta)K_2H_\eta H^{2}}{(1-\beta)C_2}\Bigr)\epsilon .
$$

Both bounds are $O(\epsilon)$, linear in the surrogate's operator-norm error. The intermediate Lemma 3.5 gives $\|\widehat G^{T}\Sigma_\eta^{-1}\widehat G\|_2>C_2$ through
$\|G^{T}\Sigma_\eta^{-1}G-\widehat G^{T}\Sigma_\eta^{-1}\widehat G\|_2\le 2\epsilon H\|\Sigma_\eta^{-1}\|_2$.

> [!warning] The scope of the theorem
> Remark 1 of the paper explains that $\widehat G_\theta$ can be made linear, as Theorem 3.6 requires, by dropping the nonlinear activations in the branch net while keeping them in the trunk net. The theorem therefore covers a restricted DeepONet, not the network used in the nonlinear experiments. Quote it with that qualification attached.

### Numerical evidence

Three benchmarks: Darcy flow, heat source inversion (Case I is a two-parameter source-location problem with DeepONet trained on $[0.5,1]^2$ from 500 uniform samples and unscented Kalman inversion started at $[0.6,0.6]$; Case II is a higher-dimensional variant), and a reaction-diffusion problem. Branch and trunk nets are fully connected with 5 hidden layers of 100 neurons and `tanh` activations; offline training runs $1\times10^{5}$ iterations on $N_{\mathrm{prior}}=1000$ Gaussian-random-field prior samples; the data are $y_{\mathrm{obs}}=y_{\mathrm{ref}}+\max\{|y_{\mathrm{ref}}|\}\delta\xi$. Both in-distribution and **out-of-distribution** truths are tested, the latter being exactly where the fixed surrogate fails. The fixed surrogate yields only a rough estimate and a visibly wrong inversion trajectory; the adaptive version drives the model error down monotonically with refinements and approaches full-order accuracy. In the third example refinement terminated after six iterations by the stopping criterion.

## Side-by-side comparison

| Component  | 37                              | 34                              | 49                               | 79                                    |
| ---------- | ------------------------------- | ------------------------------- | -------------------------------- | ------------------------------------- |
| Surrogate  | multi-fidelity polynomial chaos | multi-fidelity polynomial chaos | composite multi-fidelity network | DeepONet                              |
| Sampler    | Metropolis-Hastings             | regularising ensemble Kalman    | Metropolis-Hastings              | unscented Kalman inversion            |
| Indicator  | absolute $\ell^\infty$          | relative $\ell^\infty$          | relative $\ell^\infty$           | relative change in the data misfit    |
| Test point | accepted via the true model     | ensemble mean                   | accepted via the true model      | anchor chosen by the true model       |
| New points | random in a shrinking ball      | random in a ball                | uniform in a ball                | greedy in output space plus proximity |
| Retraining | refit least squares             | refit least squares             | train a shallow correction net   | transfer learning from prior weights  |
| Theory     | KL and Hellinger bounds         | none                            | none                             | fixed-point convergence, linear case  |

## Coverage check

| Item                                                 | Paper  | Status                                                                                    |
| ---------------------------------------------------- | ------ | ----------------------------------------------------------------------------------------- |
| Index structure of the multi-fidelity correction     | 37     | correction term, coefficient merging, weighted least squares                              |
| $\epsilon$-feasible set with KL and Hellinger bounds | 37, 79 | assumptions, conclusions, how algorithm and theorem connect                               |
| Full regularising ensemble Kalman update             | 34     | update, covariances, regularisation choice, stopping rule                                 |
| Quantitative comparison with out-of-prior truth      | 34     | relative errors and true-solve counts for both examples                                   |
| Composite multi-fidelity network                     | 49     | coupling form, contrast with a linear correction, capacity limit                          |
| Goal-oriented indicator and greedy design            | 79     | the rejected honest indicator, anchor, relative change, greedy score                      |
| Unscented Kalman inversion and linear convergence    | 79     | prediction and analysis steps, $\sigma$-points, both $O(\epsilon)$ bounds and their scope |

## Sources for this page

- L. Yan and T. Zhou, [_Adaptive multi-fidelity polynomial chaos approach to Bayesian inference in inverse problems_](https://doi.org/10.1016/j.jcp.2018.12.025), J. Comput. Phys. 381 (2019), pp. 110-128 (preprint [arXiv:1807.00618](https://arxiv.org/abs/1807.00618)).
- L. Yan and T. Zhou, [_An adaptive multifidelity PC-based ensemble Kalman inversion for inverse problems_](https://doi.org/10.1615/Int.J.UncertaintyQuantification.2019029059), Int. J. Uncertain. Quantif. 9(3) (2019), pp. 205-220 (preprint [arXiv:1809.08931](https://arxiv.org/abs/1809.08931)).
- L. Yan and T. Zhou, [_An adaptive surrogate modeling based on deep neural networks for large-scale Bayesian inverse problems_](https://doi.org/10.4208/cicp.OA-2020-0186), Commun. Comput. Phys. 28 (2020), pp. 2180-2205 (preprint [arXiv:1911.08926](https://arxiv.org/abs/1911.08926)).
- Z. Gao, L. Yan, and T. Zhou, [_Adaptive operator learning for infinite-dimensional Bayesian inverse problems_](https://doi.org/10.1137/24M1643815), SIAM/ASA J. Uncertain. Quantif. 12(4) (2024), pp. 1389-1423 (preprint [arXiv:2310.17844](https://arxiv.org/abs/2310.17844)).
