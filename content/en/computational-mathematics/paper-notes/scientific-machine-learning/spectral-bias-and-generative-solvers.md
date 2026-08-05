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

## Spectral bias is a quantity that can be measured

### Two independent discoveries of the same phenomenon

Spectral bias and the frequency principle (F-Principle) were discovered and named independently by two groups, and the two papers acknowledge each other.

Rahaman et al. state it as a learning preference for **low-frequency functions**, meaning "functions that vary globally without local fluctuations", manifesting as a **frequency-dependent learning speed**. Their mechanism is not a statement about the optimiser but about the **decay rate of the network's own Fourier spectrum**. Their Theorem 1 says: for a ReLU network $f_\theta$ (with the Fourier convention $f(\bm x)=(2\pi)^{d/2}\int\tilde f(\bm k)e^{\mathrm i\bm k\cdot\bm x}\mathrm d\bm k$),

$$
\tilde f_{\theta}(\bm k)=\sum_{n=0}^{d}
\frac{C_{n}(\theta,\bm k)\,\mathbf 1_{H^{\theta}_{n}}(\bm k)}{k^{n+1}},
\qquad k=\|\bm k\| ,
$$

where $H^{\theta}_{n}$ is a union of $n$-dimensional subspaces orthogonal to some $n$-codimensional face of a linear region $P_\epsilon$ of the network, and $C_n=\Theta(1)$ as $k\to\infty$. Two corollaries are worth recording. The decay is **anisotropic** — "$k^{-d-1}$ in almost every direction of $\mathbb R^d$, but as slow as $k^{-2}$ in the particular directions orthogonal to the $(d-1)$-dimensional faces that bound the linear regions". And the numerator is controlled by $N_fL_f$ ($N_f$ the number of linear regions, $L_f=\max_\epsilon\|W_\epsilon\|$ the Lipschitz constant),

$$
L_{f}\le\prod_{k=1}^{L+1}\|W^{(k)}\|
\le\|\theta\|_{\infty}^{L+1}\sqrt d\prod_{k=1}^{L}d_{k} .
$$

The dynamical reading is the paper's own: $L_f$ is controlled by the parameter norm, and the parameter norm can only grow gradually under gradient descent, **so high frequencies can only be learned late in the optimisation**.

Xu et al.'s frequency principle is stated directly: "deep networks often fit target functions from low to high frequencies during training." Here "frequency" means the **response** frequency — the frequency content of the input-to-output map, not frequency in the image sense; it is measured through the discrete Fourier transform $\hat f_k=\frac1n\sum_if(x_i)e^{-\mathrm i2\pi ik/n}$ and the frequency-wise relative error

$$
\Delta_F(k)=\frac{|\hat h_k-\hat f_k|}{|\hat f_k|} .
$$

The contrast with classical iterative solvers is the point most relevant to this page, and the original states it plainly: this frequency principle of deep networks is the opposite of the behaviour of the Jacobi method, which as a conventional iterative numerical scheme converges faster on **higher** frequencies across scientific-computing problems.

Their Theorem 1 (an idealised setting: one hidden layer, $\sigma=\tanh$) says: if the target satisfies $|\hat f(k_1)|>0$, $|\hat f(k_2)|>0$ and $|k_2|>|k_1|>0$, then there exist $c,C>0$ such that for sufficiently small $\delta$,

$$
\frac{\mu\bigl(\{W:|\partial L(k_{1})/\partial\theta_{lj}|>|\partial L(k_{2})/\partial\theta_{lj}|\ \forall\,l,j\}\cap B_{\delta}\bigr)}{\mu(B_{\delta})}
\ \ge\ 1-C\exp(-c/\delta),
$$

with $B_\delta$ the ball of radius $\delta$ at the origin of weight space and $\mu$ Lebesgue measure. Theorem 2 upgrades the conclusion from coordinatewise gradients to the actual rate of decrease of the loss. **Note the scope: small weights (near initialisation), one hidden layer, $\tanh$ — this is not a global statement.** The paper attributes the mechanism to the smoothness of the activation: for $\tanh$ one has $|\hat h(k)|\le C\sum_j\frac{|a_j|}{|w_j|}\exp\bigl(-\frac{\pi k}{2|w_j|}\bigr)$, so at fixed weights the network's own spectrum decays **exponentially** in $k$.

### Multi-scale networks: move high frequencies down and then learn them

The MscaleDNN construction is the shared substrate for several papers on this page. Take a band-limited target $\operatorname{supp}\hat f\subset B(K_{\max})$ and partition $k$-space into $M$ concentric annuli

$$
A_{i}=\{\bm k:(i-1)K_{0}\le|\bm k|\le iK_{0}\},\qquad K_{0}=K_{\max}/M,
$$

so that $\hat f=\sum_i\chi_{A_i}\hat f=:\sum_i\hat f_i$. Then **scale down radially**: set $\hat f_i^{\rm scale}(\bm k):=\hat f_i(\alpha_i\bm k)$ with $\alpha_i>1$, so $f_i^{\rm scale}(\bm r)=f_i(\bm r/\alpha_i)$ and

$$
\operatorname{supp}\hat f_{i}^{\rm scale}\subset
\Bigl\{\tfrac{(i-1)K_{0}}{\alpha_{i}}\le|\bm k|\le\tfrac{iK_{0}}{\alpha_{i}}\Bigr\},
$$

which is a **low**-frequency band as soon as $\alpha_i$ is large enough. By the frequency principle an ordinary network learns $f_i^{\rm scale}$ quickly, so $f(\bm r)\sim\sum_{i=1}^{M}h_i(\alpha_i\bm r,\theta^{n_i})$. The original is explicit that this derivation **is not itself an algorithm** — frequency-selective convolution in $d$ dimensions still suffers the curse of dimensionality — it only "gives a plausible form for the function space". **Scaling radially rather than coordinatewise is precisely what distinguishes MscaleDNN from the earlier PhaseDNN**, which phase-shifts along each coordinate and therefore inherits the curse of dimensionality.

There are two concrete architectures. MscaleDNN-1 splits the first hidden layer into $N$ segments, with the $i$-th segment receiving input $a_i\bm x$:

$$
f_{\bm\theta}(\bm x)=\bm W^{[L-1]}\sigma\circ\bigl(\cdots(\bm W^{[0]}(\bm K\odot\bm x)+\bm b^{[0]})\cdots\bigr)+\bm b^{[L-1]},
$$

$$
\bm K=(\underbrace{a_{1},\dots,a_{1}}_{\text{segment }1},\dots,\underbrace{a_{N},\dots,a_{N}}_{\text{segment }N})^{\top},
\qquad a_i=i\ \text{ or }\ a_i=2^{i-1} .
$$

"The only difference from a normal fully connected network is the input to the first hidden layer." MscaleDNN-2 is instead a **sum of $N$ subnetworks**, so $W^{[1]}$ through $W^{[L-1]}$ are block-diagonal. The original's empirical judgement is that MscaleDNN-2 is better and has far fewer connections, and it makes adding or removing scales easy, so all subsequent experiments use it. The activations are compactly supported (the motivation comes from compactly supported scaling functions in wavelet theory, so that scaling the activation by $\alpha$ scales its frequency content by $\alpha$):

$$
\mathrm{sReLU}(x)=\mathrm{ReLU}(-(x-1))\times\mathrm{ReLU}(x)\quad(\text{support }[0,1]),
$$

with $(\mathrm{sReLU})^2$ and $(\mathrm{sReLU})^3$ giving $C^1$ and $C^2$ smoothness, plus the quadratic-spline form

$$
\phi(x)=\mathrm{ReLU}(x)^{2}-3\,\mathrm{ReLU}(x-1)^{2}+3\,\mathrm{ReLU}(x-2)^{2}-\mathrm{ReLU}(x-3)^{2}.
$$

### A warning that points the other way

Everything above is stated under a **regression** loss. **That direction is not absolute**: paper 101 shows it reverses wholesale under a PDE residual loss, as set out in the section on 101 below. This is the most important judgement on the page, and the remaining papers should be read with it in mind.

## 81: prove an error theory first, then design from it

### Intuition

Multi-scale networks work well, but their usefulness rests on choosing the down-scaling factors $a_i$ correctly, and $a_i$ is normally fixed a priori at $2^{i-1}$. For $f(x)=\sin(\pi x)+\sin(100\pi x)$ the ideal scales are $a_1=\pi$ and $a_2=100\pi$ — unavailable in practice, because knowing them requires knowing the solution first.

A network with (random) Fourier features looks like it solves this: feed $\sin(\bm k\cdot\bm x)$ and $\cos(\bm k\cdot\bm x)$ straight into the network, and if $\bm k$ is right the network only has to fit a constant. The problem is that **being slightly off is worse**: the paper's example with $k=38\pi$ against a target frequency of $40\pi$ shows the Fourier-feature network is then worse than a plain network.

The paper's treatment has three steps. First, prove an error theory that says why down-scaling helps and why the error bound carries no frequency when the true dominant modes are used as features. Second, a **hybrid embedding**: concatenate the linear down-scaling term with the sinusoidal terms, so that the linear term catches the fall when the frequency is mismatched. Third, **posterior capture**: train a coarse solution, take its discrete Fourier transform, extract the modes carrying the most energy, use them as new features, and rebuild the network accordingly. The loop uses "the captured set stops changing" as its stopping criterion.

### Problem setup

The multi-scale ansatz is

$$
f_{\theta}(\bm x)=\sum_{i=1}^{W} a_i^{d}\, f_{\theta_i}(a_i\bm x),
$$

obtained by decomposing a compactly supported Fourier transform over concentric annuli $\mathbb K_i=\{\bm k:(i-1)K_0\le|\bm k|\le iK_0\}$ ($K_0=K_{\max}/W$) and scaling each band down by $a_i>1$, with $a_i$ normally fixed a priori at $2^{i-1}$. The random Fourier feature map is $\beta[\bm A](\bm x)=[\sin(\bm A\bm x);\cos(\bm A\bm x)]$ with the entries of $\bm A\in\mathbb R^{m\times d}$ drawn from $\mathcal N(0,\sigma^2)$. The PDE loss is

$$
\min_{\theta\in\Theta}\ w_r\mathcal L_r(\theta)+w_b\mathcal L_b(\theta),
\qquad
\mathcal L_r=\sum_{i=1}^{N_r}\bigl|\mathcal N(\bm x^i_r;u_{\text{net}})\bigr|^2,
\quad
\mathcal L_b=\sum_{i=1}^{N_b}\bigl|\mathcal B(\bm x^i_b;u_{\text{net}})\bigr|^2 .
$$

Note that per the original these two terms are **sums** rather than means.

### Derivation

**The hybrid feature embedding is motivated by one concrete computation.** Approximating $f(x)=\sin(40\pi x)$ through $F(y)=\sin(\hat k\arcsin y)$ with $\hat k=40\pi/k$ gives

$$
F'(y)=\frac{\hat k\cos(\hat k\arcsin y)}{\sqrt{1-y^2}},
\qquad
\lim_{y\to1}|F'(y)|=\infty ,
$$

so as soon as $|\hat k|\ne1$ one has $\|F\|_{C^1[-1,1]}=\infty$ and the hypothesis of Theorem 3.3 fails. **This is the precise reason pure Fourier features degrade under frequency mismatch**, and the linear term is added exactly to cover this case. The new input map therefore combines the multi-scale linear down-scaling with sinusoidal Fourier features:

$$
\Phi[\bm k](\bm x)=\begin{bmatrix}\bm k\cdot\bm x\\ \cos(\bm k\cdot\bm x)\\ \sin(\bm k\cdot\bm x)\end{bmatrix}.
$$

**Posterior capture of the dominant modes.** Take the discrete Fourier transform of the current network solution,

$$
\hat f_{\text{net},\bm k}=\int_{[0,1]^d}f_{\text{net}}(\bm x)e^{-\mathrm i2\pi\bm k\cdot\bm x}\,\mathrm d\bm x,
$$

and keep the $N_0$ coefficients of largest modulus, where $N_0$ is fixed by an energy-fraction condition:

$$
\sum_{j=1}^{N_0}\bigl|\hat f_{\text{net},\bm k_j}\bigr|^2
\ \ge\ (1-\delta)\,\|f_{\text{net}}\|^2_{L^2([0,1]^d)},
\qquad 0\le\delta<1 .
$$

The working assumption for this step is $\|f_{\text{net}}-f\|_{L^2([0,1]^d)}\le\epsilon$ with $0<\epsilon\ll\|f\|_{L^2}$ — that is, the coarse solution is already good enough for its spectrum to stand in for the true one.

**Rebuild the network by two criteria.** Criterion A: if the number of captured modes satisfies $N_0\le M_0$ (the number of subnetworks), build $N_0$ subnetworks, with the $j$-th taking input $\Phi[\bm k_j](\bm x)$ and output

$$
y=\sum_{j=1}^{N_0} h_j\,y_j,\qquad h_j=\hat u_{\text{net},0,\bm k_j},
$$

reusing the coefficients to accelerate convergence. Criterion B: if $N_0>M_0$, **do not enlarge the network**; instead split $\mathbb B_1=\{\bm k_1,\dots,\bm k_{N_0}\}$ into $M_0$ blocks (typically $\mathbb B_1^1=\{\bm k_1,\dots,\bm k_{\lfloor N_0/M_0\rfloor}\}$ and so on), let the $j$-th subnetwork receive $\Phi[\bm k](\bm x)$ for each $\bm k$ in that block, and output

$$
y=\sum_{j=1}^{M_0}h_j\,y_j,\qquad h_j=\sum_{\bm k\in\mathbb B^j_1}\bigl|\hat u_{\text{net},0,\bm k}\bigr| .
$$

One implementation detail is worth recording (their Remark 4.1): the last layer is written $u_{\text{net}}=WG+b$ with $W=(h_1,h_2,\dots)$ and $b=0$ **fixed and non-trainable**, and the paper reports this is both more accurate and cheaper.

**Algorithm.** Given the number of subnetworks $M_0$, an initial feature set $\mathbb B_0=\{2^0,\dots,2^{M_0-1}\}$, a total number of adaptive steps $I$ and a threshold parameter $\lambda$: repeatedly train the multi-scale network under the current feature set, take the discrete Fourier transform of the solution, add every $\bm k$ satisfying $|\hat u_{\text{net},It,\bm k}|>\lambda\max_{\bm k}|\hat u_{\text{net},It,\bm k}|$ to the next feature set, and rebuild the network by criterion A or B. Two stopping criteria: a fixed maximum number of rounds $I$, and the captured feature set ceasing to change ($\mathbb B_{It}=\mathbb B_{It+1}$).

### Theorems

Four theorems, all stated in the sense of function approximation; the statements were checked in this pass, the proofs were not reviewed.

**Theorem 3.1** gives the fitting error bound for a **standard** network (no scaling map), as a reference baseline.

**Theorem 3.2 (why down-scaling works).** Let $f:\mathbb R^d\to\mathbb R$ be $(p,C_0)$-smooth with $p=q+s$, $q\in\mathbb N_0$, $s\in(0,1]$, and let $M$ be sufficiently large. Then within the feature-map function class $\widetilde{\mathcal F}(1,L,\iota,\delta)$ there is a network $f_{\text{net}}$ with the **down-scaling** map $\Phi(\bm x)=k\bm x$ such that

$$
\|f_{\text{net}}-f\|_{\infty,[-h,h]^d}
\ \le\ \frac{C_2\bigl(\max\{kh,\,C_1\}\bigr)^{5q+3}}{M^{2p}} .
$$

Compared with Theorem 3.1, the point is that **shrinking $kh$ shrinks the numerator**.

**Theorem 3.3** handles compositions of the form $f=F(\sin(\bm k\cdot\bm x))$, and its key hypothesis is $\|F\|_{C^1}<\infty$ — the computation $\lim_{y\to1}|F'(y)|=\infty$ above says exactly when this hypothesis fails.

**Theorem 3.4 (a frequency-free bound for band-limited targets).** Let $f\in L^2([-\pi,\pi]^d)$ be band-limited,

$$
f(\bm x)=\sum_{\bm k\in\mathbb B}\bigl[b_{\bm k}\cos(\bm k\cdot\bm x)+c_{\bm k}\sin(\bm k\cdot\bm x)\bigr],
\qquad |\mathbb B|=N,
$$

with bounded coefficients. Take Fourier features $\Phi(\bm x)=(\dots,\sin(\bm k\cdot\bm x),\cos(\bm k\cdot\bm x),\dots)_{\bm k\in\mathbb B}$ together with $L=C_4(1,q)$, $\iota=C_5(1,q,M)$ and $\delta=C_6(1,M,1,1)$. Then there is an $f_{\text{net}}\in\widetilde{\mathcal F}(2N,L,\iota,\delta)$ with

$$
\|f_{\text{net}}-f\|_{\infty,[-\pi,\pi]^d}\ \le\ \frac{C_2\sqrt N}{M^{2p}},
\qquad M^{2p}\ge\max\{2c_2,c_3,12\} .
$$

**This bound does not contain the frequency.** That is exactly the theoretical justification for building the embedding out of the actual dominant modes.

### Numerical experiments

Section 5 covers four classes of equation:

| Subsection | Equation                                           | Recorded adaptive trajectory                                                                   |
| ---------- | -------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| 5.1        | Poisson equation                                   | captured feature set $\{2\pi,4\pi,200\pi,202\pi\}$, rebuilt as four subnetworks by criterion A |
| 5.2        | heat equation                                      | rebuilt by criterion B                                                                         |
| 5.3        | wave equation                                      | —                                                                                              |
| 5.4        | Schrödinger equation near the semi-classical limit | rebuilt by criterion B                                                                         |

Section 5.1 additionally reports one ablation: **fixing** $W$ and $b$ in the last layer beats making them learnable, on both accuracy and cost.

The paper's headline conclusion is that the frequency-adaptive version improves accuracy over a standard multi-scale network by **two to three orders of magnitude**. Its reference point is worth recording too: the original MscaleDNN literature (Cai–Xu) reports errors typically of order $10^{-2}$ to $10^{-1}$, which is precisely the motivation for a posterior refinement.

What this set of experiments establishes is that **reading the frequencies off the solution's own spectrum is more reliable than assuming them a priori**, and that the loop converges to a stable feature set on all four classes of equation. What it does not establish is a link between theory and experiment — Theorems 3.2 and 3.4 are existence bounds in the approximation-theoretic sense and do not guarantee that optimisation finds that network, and the paper does not measure the gap between the realised error and the bound.

> [!warning] Per-example error values not transcribed
> The specific error values for each example live in the figures and tables of Section 5 and were not transcribed item by item in this pass, so this page does not list those numbers and records only the verified structure, ablations and headline order of magnitude.

### Relation to the other papers

This is the starting point of the group's "spectral bias" strand. Paper 101 attacks the same problem from the other side: rather than rebuilding the network around modes captured by a discrete Fourier transform, it keeps a fixed multi-scale random Fourier basis bank and learns a set of **input-dependent attention weights**, while retaining the idea of incrementally enriching features from the spectrum of an intermediate solution. Paper 94 is a direct application: it uses a multi-scale network (not a PINN) on the Gaussian-wave-packet ODE system derived from the semi-classical Schrödinger equation, which is exactly the example in Section 5.4 here; the first two authors are the same on both. Paper 90's subdomain-local scaled bases are the counterpart of the same down-scaling idea inside random-basis methods, discussed on the [[en/computational-mathematics/paper-notes/scientific-machine-learning/variational-and-basis-networks|variational and basis networks page]].

## 94: the semi-classical limit — eliminate the spatial dimension outright

### Intuition

When $\varepsilon\ll1$ the solution of the semi-classical Schrödinger equation oscillates on the scale $O(\varepsilon^{-1})$ in both space and time. Any direct discretisation — finite differences, spectral methods, physics-informed networks — has to resolve the solution at that scale, and the cost explodes as $\varepsilon$ shrinks.

But the **shape** of those oscillations is known: the solution is approximately a Gaussian envelope times a phase factor $\exp[\mathrm i(\cdot)/\varepsilon]$, and the envelope's centre, width, momentum and phase are a handful of scalars evolving in time. Substituting that form into the equation as an ansatz gives a system of ordinary differential equations for those few scalars, with **the spatial variable $x$ eliminated entirely**. What the network has to learn is no longer a highly oscillatory function but four curves in time; the high-frequency structure has been put analytically into the reconstruction formula.

The cost is equally concrete: because those scalars sit inside $\exp[\mathrm i(\cdot)/\varepsilon]$, the error in the ODE solution is amplified by $1/\varepsilon$ on reconstruction. So the ODEs must be solved extremely accurately — which is exactly why the paper switches from a plain PINN to a multi-scale network.

A second, independent problem is that changing the initial condition forces retraining, so an operator-learning formulation is needed. The paper's answer is to factor the operator into three stages and learn only the finite-dimensional middle one.

### Problem setup

The target equation is

$$
\psi_t=\frac{\mathrm i\varepsilon}{2}\Delta\psi-\frac{\mathrm i}{\varepsilon}V(\bm x)\psi,
\qquad
\psi(\bm x,0)=\varphi(\bm x)\exp\bigl(\mathrm i\phi(\bm x)/\varepsilon\bigr),
$$

with $\varepsilon$ the dimensionless Planck constant. When $\varepsilon\ll1$ the solution "oscillates at high frequency on the scale $O(\varepsilon^{-1})$ in both space and time". Classical splitting solvers (Strang splitting: spectral accuracy in $\Delta x/\varepsilon$, $O(\Delta t^2/\varepsilon)$ in time; Chin–Chen: spectral in space, fourth order in time) all carry a $1/\varepsilon$ factor in their global truncation error, so mesh and time step must shrink together with $\varepsilon$.

Plain physics-informed networks do worse. The PINN baseline loss is $\mathcal L(\theta)=w_r\mathcal L_r+w_b\mathcal L_b+w_i\mathcal L_i$ with

$$
\mathcal L_r=\Bigl\|\psi_t-\tfrac{\mathrm i\varepsilon}{2}\psi_{xx}+\tfrac{\mathrm i}{\varepsilon}V\psi\Bigr\|^2_{2,\Omega\times[0,T]},
\quad
\mathcal L_b=\|\mathcal P(\psi)-\psi\|^2_{2,\partial\Omega\times[0,T]},
\quad
\mathcal L_i=\bigl\|\psi(\cdot,0)-\varphi e^{\mathrm i\phi/\varepsilon}\bigr\|^2_{2,\Omega},
$$

where $\mathcal P$ is the periodic boundary operator. The paper gives the evidence directly: applying a naive network to this equation on $(-2\pi,2\pi)\times[0,1]$ gives a relative $L^2$ error **already above 1** at $\varepsilon=0.1$.

### Derivation

**Step 1: the Heller Gaussian wave-packet ansatz.** In one dimension,

$$
\psi(x,t)=\exp\left[\frac{\mathrm i}{\varepsilon}
\Bigl(\alpha(t)\bigl(x-q(t)\bigr)^2+p(t)\bigl(x-q(t)\bigr)+\gamma(t)\Bigr)\right].
$$

The envelope is centred at $q(t)$ with standard deviation proportional to $\sqrt{\varepsilon/\alpha_{\rm im}(t)}$; the oscillation wavelength at $x=q(t)$ is $2\pi\varepsilon/p(t)$, and $\mathrm{Re}\,\alpha$ produces still finer oscillations in the tails.

**Step 2: the key reduction.** Substituting the ansatz into the equation makes it an **exact** solution for a quadratic potential and Gaussian initial data, provided

$$
\dot q=p,\qquad
\dot p=-V'(q),\qquad
\dot\alpha=-2\alpha^{2}-\tfrac12V''(q),\qquad
\dot\gamma=\tfrac12p^{2}-V(q)+\mathrm i\,\alpha\varepsilon .
$$

**The spatial variable $x$ is eliminated entirely.** For non-harmonic potentials the Gaussian wave packet has modelling error $O(\sqrt\varepsilon)$, so it **improves** as $\varepsilon\to0$ (this is cited from Heller's work, not proved here).

**Step 3: the crucial fact about error amplification.** Because the ODE unknowns sit inside $\exp[\mathrm i(\cdot)/\varepsilon]$, "if the total error of the PINN solution of the ODE system is $\mathcal E_t$, the total error of the reconstructed solution is $\mathcal E_t/\varepsilon$". Equivalently, solving the ODEs with a $k$-th order time integrator gives a total truncation error of order $(\Delta t)^k/\varepsilon$. **This $1/\varepsilon$ amplification is why the ODEs have to be solved extremely accurately**, and why the paper switches to a multi-scale network.

**Step 4: the loss for the ODE system.** The network output $\bm y_\theta=(q_{\theta_1},p_{\theta_2},\alpha_{\theta_3},\gamma_{\theta_4})$ takes only $t\in[0,T]$ as input, and the loss collocates the four residuals; the $\gamma$ term, for example, is

$$
\bigl|\dot\gamma_\theta(t_i)-\tfrac12p_\theta^2(t_i)+V(q_\theta(t_i))-\mathrm i\alpha_\theta(t_i)\varepsilon\bigr|^2 ,
$$

with the reconstructed field

$$
\psi(x,t;\theta)=\exp\left[\frac{\mathrm i}{\varepsilon}
\Bigl(\alpha_{\theta_3}(t)(x-q_{\theta_1}(t))^2+p_{\theta_2}(t)(x-q_{\theta_1}(t))+\gamma_{\theta_4}(t)\Bigr)\right].
$$

**Step 5: the specific choice of multi-scale network.** The paper restates the MscaleDNN down-scaling construction as

$$
f_\theta(\bm z)=\sum_{i=1}^{M}a_i^{d+1}f_{\theta_i}(a_i\bm z),
$$

but **the experiments here do not use the usual geometric scales $a_i=2^{i-1}$; they use 100 linearly spaced scales $\{a_i\}=\{0.1,0.2,\dots,10.0\}$**. The reason is the observation that the Fourier spectrum of $q(t)$ has a long tail that geometric scales do not cover densely enough.

**Step 6: operator decomposition.** Rather than learning $\mathcal G:\psi(\cdot,0)\mapsto\psi(\cdot,t)$ directly, factor it through the wave-packet parameters:

$$
\mathcal G:\ \psi(x,0)\ \longrightarrow\ \bm y(0)
\ \overset{\mathcal G'}{\longrightarrow}\ \bm y(t)\ \longrightarrow\ \psi(x,t),
\qquad \bm y=(q,p,\alpha,\gamma),
$$

and learn only the **finite-dimensional** $\mathcal G':\bm y(0)\mapsto\bm y(t)$. **The key point is that the branch input is $\bm y(0)$ and the trunk input is only $t$, neither of which involves $x$** — this is where the oscillatory $x$-dependence is removed from the learning problem. The architecture is

$$
\mathcal G'^{(i)}_\theta[\bm y(0)](t)=\sum_{j=J_{i-1}+1}^{J_i}b_j[\bm y(0)]\,c_j(t),
\quad i=1,\dots,I,
$$

that is, the $J$ basis indices are split into $I$ blocks, one per output component. The loss is entirely **data-free**:

$$
\hat{\mathcal L}(\theta)=
\underbrace{\frac1N\sum_{i=1}^{N}\sum_{j=1}^{I}
\bigl|\mathcal G'^{(j)}_\theta[\bm y^{(i)}(0)](0)-y^{(i)}_j(0)\bigr|^2}_{\hat{\mathcal L}_{\rm boundary}}
+\underbrace{\frac{1}{NQ}\sum_{i=1}^{N}\sum_{j=1}^{Q}\sum_{l=1}^{I}
\Bigl|\mathcal N^{(l)}\bigl(\mathcal G'_\theta[\bm y^{(i)}(0)](t_j^{(i)})\bigr)\Bigr|^2}_{\hat{\mathcal L}_{\rm physics}},
$$

where $(\mathcal N^{(1)},\dots,\mathcal N^{(I)})$ are the residual operators of the ODE system.

There are two further routes for general (non-Gaussian) initial data: Appendix A generalises to **Hagedorn wave packets** (outputs $q,p,Q,P,S$), and Appendix B uses a **Gaussian-beam decomposition** to split WKB-type initial data into $\mathfrak N$ beams. The existence and scope of these two appendices were verified; their internal formulas were not read closely in this pass.

### Theorems

**This paper proves no theorems.** Both the $O(\sqrt\varepsilon)$ modelling error of the Gaussian wave packet and the $\mathcal E_t/\varepsilon$ amplification are cited from prior work (Heller and related literature), not proved here. The summary says so explicitly: a robust theoretical framework explaining why multi-scale networks outperform PINNs on solving the ODE system is still missing.

There is also one explicit limitation (their Remark 3.1): a single Gaussian wave packet is only accurate within the **Ehrenfest time**, so "the neural network solver proposed here is not suitable for long-time simulations", even when the ODEs themselves are still solved accurately.

### Numerical experiments

Reference solutions: a fourth-order scheme for non-harmonic potentials, RK4 on the wave-packet ODEs for harmonic ones. The metric is relative $L^2$ over $\Omega\times[0,T]$. The activation is a **softened Fourier mapping** $\sigma(x)=0.5\sin x+0.5\cos x$, with Glorot normal initialisation, Adam, learning rate $10^{-3}$, 100,000 epochs and batch size 1000.

**Example 1: the free-particle ODE system.** $V\equiv0$, $\varepsilon=0.01$, with exact solution $q=1+2t$, $p=2$, $\alpha=(2t+\mathrm i)/(4t^2+1)$ and $\gamma=2t-0.005\arctan(2t)+\tfrac14(\log(4t^2+1)-\log(200/\pi))\mathrm i$. The PINN is $[1,100,400,400,400,400,6]$ and the multi-scale network is $[100,400,400,400,400,6]$ with 100 linear scales.

| Component | PINN                 | MscaleDNN            |
| --------- | -------------------- | -------------------- |
| $q$       | $3.287\times10^{-5}$ | $2.928\times10^{-7}$ |
| $p$       | $2.585\times10^{-4}$ | $4.141\times10^{-6}$ |
| $\alpha$  | $3.105\times10^{-3}$ | $1.941\times10^{-5}$ |
| $\gamma$  | $2.864\times10^{-4}$ | $4.161\times10^{-6}$ |

The diagnostic is clear too: the discrete Fourier transform of $q(t)$ has a long-tailed spectrum, the PINN's **error spectrum has the same shape**, and the multi-scale network's error spectrum is far more uniform. Accuracy keeps improving as the number of embeddings grows from 10 to 100.

**Example 2: the 1D harmonic potential $V=\tfrac12x^2$.** For $\varepsilon\in\{4/25,1/25,1/100,1/400,1/1600,1/6400\}$ the PINN error grows from $1.522\times10^{-3}$ to $6.877\times10^{-1}$, and the multi-scale network from $3.192\times10^{-5}$ to $7.956\times10^{-3}$. The error grows **linearly** as $\varepsilon$ shrinks, consistent with both the $1/\varepsilon$ amplification and the vanishing modelling error in the harmonic case. A sampling-density study shows 6400 points suffice, and 10,000 are used thereafter. **Applying a standard multi-scale network directly to the original equation already fails at $\varepsilon=1/32$** — a clean ablation isolating the wave-packet reduction on its own.

**Example 3: the 1D torsional potential $V=1-\cos x$.** At large $\varepsilon$ the PINN and the multi-scale network are nearly identical because the $O(\sqrt\varepsilon)$ **modelling** error dominates; the multi-scale network only pulls ahead below $\varepsilon=1/100$ ($3.408\times10^{-3}$ against $1.048\times10^{-2}$ at $\varepsilon=1/1600$). This is a good illustration of the two-term budget of modelling error plus solution error.

**Example 4: two and four dimensions.** In two dimensions $q$ and $p$ become vectors and $\alpha$ becomes a complex symmetric matrix $A$, for output dimension 12. In four dimensions $V(\bm x)=\tfrac12\bm x^{T}A_p\bm x$ with $A_p$ having diagonal $(1.0,2.0,2.0,3.0)$ and all off-diagonal entries $0.2$, $\bm q(0)=(1.3,0,-1,0)^T$ and $\bm p(0)=(0,1.3,0,1)^T$; $A$ needs 20 output components, so four independent networks handle $(\bm q,\bm p)$, $A$, $\gamma_{\rm re}$ and $\gamma_{\rm im}$ separately.

| $\varepsilon$ (4D) | PINN                | MscaleDNN           |
| ------------------ | ------------------- | ------------------- |
| $0.01$             | $2.08\times10^{-4}$ | $7.48\times10^{-5}$ |
| $0.001$            | $2.07\times10^{-3}$ | $6.47\times10^{-4}$ |

**Example 5: DeepONet on the 1D torsional potential.** $\varepsilon=0.01$, $T=1$; $q(0)\sim\mathcal U[0.8,1.8]$, $p(0),\alpha_{\rm re}(0)\sim\mathcal U[-0.5,0.5]$, $\alpha_{\rm im}(0)\sim\mathcal U[0.5,1.5]$, $\gamma_{\rm re}(0)=0$, and $\gamma_{\rm im}(0)$ fixed by normalisation; $I=6$ and $J_1=\dots=J_I=100$; branch $[6,100,100,100,100,600]$ and trunk $[1,100,100,100,100,100]$; $N=2000$ initial functions, $Q=500$ time points, 200,000 iterations. The errors at $q(0)=\pi/2,1,1.5$ are $2.69\times10^{-2}$, $2.39\times10^{-2}$ and $1.67\times10^{-2}$; over $N_{\rm test}=100$ the mean is $2.44\times10^{-2}$ with standard deviation $1.54\times10^{-2}$.

**Example 6: timing DeepONet against RK4.** In the harmonic case:

| $N_{\rm test}$ | RK4     | DeepONet |
| -------------- | ------- | -------- |
| 10000          | 11.36 s | 4.15 s   |
| 20000          | 22.61 s | 4.72 s   |
| 40000          | 45.63 s | 5.44 s   |

**DeepONet's inference cost barely grows with $N_{\rm test}$.** In a further run generating 100 grid points for each of 10,000 cases, DeepONet takes 4.58 s against 67.85 s for RK4 at $\Delta t=0.01$. The timing comparison is honestly calibrated — RK4's $\Delta t=0.1$ was chosen so that its error matches DeepONet's. The mean component errors (harmonic, 1D) are $5.33\times10^{-5}$ for $q$, $6.95\times10^{-5}$ for $p$, $4.04\times10^{-4}$ for $\alpha$ and $1.85\times10^{-4}$ for $\gamma$, while **$\psi$ is $1.93\times10^{-2}$** — that two-order jump from ODE parameters to $\psi$ is exactly the $1/\varepsilon=100$ amplification.

The 2D DeepONet study uses $V=(x_1^2+x_2^2)/2$; Appendix A redoes the same study with Hagedorn wave packets; Appendix B studies the Gaussian-beam decomposition of WKB-type initial data and reports runtime and accuracy against the number of beams $\mathfrak N$.

What these experiments establish is that **once the oscillatory structure is moved analytically out of the learning problem, the remaining low-dimensional problem can be solved very accurately**, and that the multi-scale network still beats the PINN substantially on that low-dimensional problem. What they do not establish is long-time behaviour — beyond the Ehrenfest time a single wave packet fails on its own terms, independently of the network.

### Relation to the other papers

The multi-scale network and paper 101's multi-scale random Fourier basis bank attack the same spectral bias: paper 94 uses 100 **a priori fixed** linear scales, while paper 101 learns amplitudes over a basis bank through cross-attention and can additionally **inject** frequencies found by a discrete Fourier transform. Both use the discrete Fourier transform of the target or the error as a diagnostic; paper 101 is the more general machine and paper 94 the domain-specific application. Paper 90's subdomain-local scaled bases are the counterpart of the same down-scaling idea inside random-basis methods.

The physics-informed DeepONet used here is exactly the operator-learning substrate that paper 95 generalises to a permutation-invariant, uncertainty-aware version, and the one that paper 103 replaces with an operator-splitting architecture; all three are physics-informed, data-free operator learning aimed at evolution problems. And the strategy of reducing a PDE to a small ODE system and then learning that system's flow map is structurally the same move as paper 107's reduction of Fokker–Planck operator learning to a transition density with a linearised-SDE base distribution, discussed on the [[en/computational-mathematics/paper-notes/scientific-machine-learning/uncertainty-aware-operator-learning|operator learning and uncertainty page]].

## 101: a fixed basis bank with learnable weights

### Intuition

The most worthwhile thing to record about this paper is not its architecture but an observation in its Section 3.1: **under physics-informed training the direction of spectral bias reverses.**

The reason is a scaling on the Fourier side. Since $\widehat{\partial_xu}(k)=\mathrm ik\,\hat u(k)$, a Deep Ritz energy term $\|\partial_xu\|_{L^2}^2$ weights mode $k$ by $k^2$; since $\widehat{u_{xx}}(k)=-k^2\hat u(k)$, a squared PDE residual weights mode $k$ **like $k^4$**. In other words, the differential operator itself **amplifies** high frequencies, so what ends up under-resolved under a residual loss is the **low**-frequency part.

The paper measures this: for $u=\sin(\pi x)+\sin(5\pi x)+\sin(20\pi x)$ it looks at the frequency-wise relative error $\Delta_F(k)=|\hat u_k^{\rm pred}-\hat u_k|/|\hat u_k|$ at $k\in\{1,5,20\}$, and finds that under a pure regression loss the low frequencies do converge first, while under Deep Ritz and PINN losses the high frequencies decay faster.

The intuition on the architecture side is: since one does not know in advance which frequencies the solution contains, prepare a **sufficiently wide** fixed frequency basis bank and let the model decide which band to emphasise in which part of the domain. The mechanism is cross-attention — queries come from the network's hidden state, keys and values from the frequency basis bank, and the softmax weights vary with the input. The bank itself is not trained; only a decay envelope acting on the amplitudes is. If a key frequency turns out to be missing from the bank, modes found by a discrete Fourier transform are **appended** as new tokens rather than the bank being rebuilt.

> [!note] Appendix A not verified
> The main-text Fourier scaling argument ($k^2$ and $k^4$) has been verified. The analysis in Appendix A of "the differential operator amplifies high frequencies" was not read through in this pass and is treated as unverified.

### Problem setup

Spectral bias means networks "converge rapidly to low-frequency components ... yet struggle to represent high-frequency or highly oscillatory features". The paper's criticism of existing fixes is specific: random Fourier features, MscaleDNN, FMMNN, phase-shift methods and random feature methods "typically rely on pre-specified frequency bases or prescribed multi-scale transformations, which may limit their adaptivity to instance-dependent spectral demands"; FG-PINN "depends strongly on the presence of informative high-frequency content in the source term or in the initial-boundary data". The paper's thesis is that overcoming spectral bias needs not just a richer frequency dictionary but an architectural mechanism able to route and reweight spectral components dynamically according to the input and the evolving solution structure.

### Derivation

**Step 1: the multi-scale random Fourier basis bank.** Base frequencies

$$
\omega_m\sim\mathcal N\bigl(\bm 0,\ \sigma^{-2}I_{d_{\rm in}}\bigr),\qquad m=1,\dots,M_{\rm base},
$$

are expanded over dyadic scales $\widetilde\omega_{m,k}=2^k\omega_m$ ($k=0,\dots,K$) into $\overline\Omega\in\mathbb R^{M\times d_{\rm in}}$ with $M=M_{\rm base}(K+1)$, and phases $b_{m,k}\sim\mathrm{Uniform}(0,2\pi)$ are drawn once. The feature map carries a **learnable amplitude envelope**:

$$
\phi(x)=\sqrt{\tfrac1M}\Bigl[a_{m,k}\cos\bigl(\widetilde\omega_{m,k}^{\top}x+b_{m,k}\bigr)\Bigr]_{(m,k)}\in\mathbb R^M,
\qquad
a_{m,k}=\exp\bigl(-\beta\|\widetilde\omega_{m,k}\|_2\bigr).
$$

**The frequencies themselves are fixed and not learnable; the only learnable quantity is the decay rate $\beta\ge0$** (kept positive through a softplus parameterisation). This is exactly complementary to paper 81's trade-off: 81 makes the features align with the true frequencies, 101 uses a sufficiently wide fixed bank plus a learnable decay envelope.

**Step 2: rearrange the features into tokens.** Choose a token width $d_q\mid M$ and reshape $\phi(x)\in\mathbb R^{M}$ into $H(x)\in\mathbb R^{N_{\rm tok}\times d_q}$ with $N_{\rm tok}=M/d_q$. "A simple choice of $d_q$ is $d_q=M_{\rm base}$, so that each row of $H(x)$ corresponds to one scale" — **one token per scale, which makes attention over tokens literally attention over frequency bands.**

**Step 3: the cross-attention residual stack.** The initial hidden state is

$$
Q^{(0)}(x)=\sigma\bigl(W^{(0)}\psi(x)+b^{(0)}\bigr)\in\mathbb R^{d_q},
$$

where $\psi(x)=x$ gives the **NN-CA** variant and $\psi(x)=\phi(x)$ gives the **RFF-CA** variant. Writing $Q_l=Q^{(l)}W_Q^{(l)}$, $K_l=H(x)W_K^{(l)}$ and $V_l=H(x)W_V^{(l)}$ with $W_\bullet^{(l)}\in\mathbb R^{d_q\times d_q}$, for $l=0,\dots,L-1$,

$$
\mathrm{CA}\bigl(Q^{(l)},H\bigr)=\mathrm{softmax}\Bigl(\frac{Q_lK_l^{\top}}{\sqrt{d_q}}\Bigr)V_l,
$$

$$
\widetilde Q^{(l)}=Q^{(l)}+\mathrm{CA}\bigl(Q^{(l)},H\bigr),
\qquad
Q^{(l+1)}=\widetilde Q^{(l)}+\sigma\bigl(W^{(l)}\widetilde Q^{(l)}+b^{(l)}\bigr),
$$

with output $u_\theta(x)=W_{\rm out}Q^{(L)}(x)+b_{\rm out}$. The queries come from the **hidden state** and the keys and values from the **frequency basis bank** — hence cross-attention rather than self-attention. The point is that the softmax weights vary with the input, so different regions of the domain can emphasise different frequency bands.

> [!warning] A typo in the preprint
> The sentence describing the NN-CA and RFF-CA cases writes $\psi(x)=\phi(x)$ twice; from the context the first case should be $\psi(x)=x$.

**Step 4: adaptive frequency enhancement (AFE).** Train a preliminary model $u_\theta^{(0)}$, take a discrete Fourier transform on a uniform grid over the periodic domain $\Omega$ to get $\hat u_{\theta,k}$ (index set $B$), set $\zeta=\max_{k\in B}|\hat u_{\theta,k}|$ and select by a relative threshold

$$
\mathcal K_{\rm post}=\bigl\{k\in B:\ |\hat u_{\theta,k}|>\lambda\zeta\bigr\},\qquad 0<\lambda<1 .
$$

From these, build deterministic frequencies $\omega_k^{\rm post}=2k\pi$ and posterior features

$$
\phi_{\rm post}(x)=\sqrt{\tfrac{2}{M_{\rm post}}}\cos\bigl(\Omega_{\rm post}x+b_{\rm post}\bigr),
\qquad
\Omega_{\rm post}=[\omega_k^{\rm post}]_{k\in\mathcal K_{\rm post}},\ M_{\rm post}=|\mathcal K_{\rm post}| ,
$$

reshape them into $H_{\rm post}$ (zero-padding if $d_q\nmid M_{\rm post}$) and concatenate along the **token dimension** into $H_{\rm aug}=[H_{\rm base};H_{\rm post}]$. The new tokens are not switched on abruptly but released smoothly through an additive log-mask:

$$
A^{(l)}=\frac{Q_lK_l^{\top}}{\sqrt{d_q}}+\mathcal M^{(l)},
\qquad
\mathcal M^{(l)}=[\,\bm 0;\ \eta_l\bm 1\,],\quad \eta_l\le0,
$$

with the zero block acting on $H_{\rm base}$ and the constant block $\eta_l$ on $H_{\rm post}$, and $\eta_l\uparrow0$ gradually releasing the injected tokens. So it **augments** the bank rather than rebuilding it, and the backbone is untouched — which is precisely the dividing line with paper 81. Remark 2.1 records two practical caveats: for complex domains, embed in a hypercube first and then mask; in high dimensions, extract modes from one-dimensional component functions rather than taking a full $d$-dimensional discrete Fourier transform.

**Step 5: the two-network PDE representation.** For $\mathcal N[u]=f$ in $\Omega$ with $u=g$ on $\partial\Omega$, take

$$
u(\bm x;\theta)=u_h(\bm x;\theta_h)+\alpha\,u_\ell(\bm x;\theta_\ell),
$$

with $u_h$ the cross-attention network above (handling high frequencies) and $u_\ell$ a plain fully connected network (handling low frequencies), the boundary condition being split as $u_h=g$ and $u_\ell=0$. The loss is

$$
L=\int_\Omega\bigl(\mathcal N[u_h+\alpha u_\ell]-f\bigr)^2\rho_r\,\mathrm d\bm x
+\gamma\int_{\partial\Omega}\bigl((u_h-g)^2+\alpha^2u_\ell^2\bigr)\rho_b\,\mathrm ds .
$$

There are two choices for the mixing coefficient $\alpha$: a trainable scalar, or a **closed-form optimum**. The latter comes from solving $\partial L/\partial\alpha=0$ when $\mathcal N$ is linear:

$$
\alpha_{\rm opt}=-\frac{\displaystyle\int_\Omega\bigl(\mathcal N[u_h]-f\bigr)\,\mathcal N[u_\ell]\,\rho_r\,\mathrm d\bm x}
{\displaystyle\int_\Omega\bigl(\mathcal N[u_\ell]\bigr)^2\rho_r\,\mathrm d\bm x+\gamma\int_{\partial\Omega}u_\ell^2\,\rho_b\,\mathrm ds},
$$

estimated by Monte Carlo on the collocation points. In practice, to stabilise the estimate, $\alpha_{\rm opt}$ is dropped from the boundary penalty and

$$
\widehat L\approx\frac{1}{N_r}\sum_i\bigl(\mathcal N[u_h+\alpha_{\rm opt}u_\ell]-f\bigr)^2(x_r^{(i)})
+\frac{\gamma}{N_b}\sum_i\Bigl((u_h-g)^2+u_\ell^2\Bigr)(x_b^{(i)})
$$

is used instead.

### Theorems

**There are no theorems.** Appendix A is described as "a simple analysis of the amplification of high frequencies by the differential operator", not a convergence theorem; its details were not verified in this pass. The conclusion says so plainly: a rigorous analysis of the approximation and optimisation dynamics of cross-attention-based multi-scale Fourier representations is still missing. Other self-declared open problems include that the discrete-Fourier-transform-guided enhancement "is most natural for periodic or grid-friendly settings", and the extension to physics-informed **operator** learning.

### Numerical experiments

Shared default hyperparameters: $m_{\rm base}=128$, $K=3$, width $d_q=64$, $n_{\rm heads}=4$, $L=4$ blocks, $\beta_0=0.1$, double precision, Adam with gradients clipped at 1.0. The metric is relative $L^2$, with PSNR and HFEN added for the image task.

**Example 1: 2D synthetic regression.** Three deliberately hostile targets on $[-1,1]^2$: $f_1$ is a composite (a logistic angular gate producing sector anisotropy, a band-pass gate producing a narrow high-frequency ring, a localised Gabor-type spiral wave packet, a weak star-shaped discontinuity $0.12\,\mathrm{sign}(r-r_*(\theta))$ with $r_*(\theta)=0.55+0.10\cos5\theta$, plus a smooth cross term); $f_2=\cos(2\pi[(w_0+w_1r)(x_1\cos\kappa\theta+x_2\sin\kappa\theta)])$ with $\kappa=5$, $w_0=4$, $w_1=3$ (a vortex whose instantaneous frequency varies with radius); and $f_3=\mathrm{sign}(\sin2\pi f_xx_1\sin2\pi f_yx_2)$ with $f_x=f_y=1$ (a discontinuous checkerboard). A $500\times500$ grid, batch size 4000, Adam at learning rate $2\times10^{-3}$, StepLR parameters $(100,0.5)$ for $f_1,f_2$ and $(50,0.5)$ for $f_3$, trained for 1000/1000/500 epochs. Conclusion: RFF-CA beats RFF-NN on loss and error for all three; which of NN-CA and RFF-CA is better is problem-dependent. The paper gives curves rather than numerical tables.

**Example 2: image regression as function fitting (DIV2K).** Four validation images (largest $2040\times1536$), coordinates mapped into $[-1,1]^2$ via $x_j=\tfrac{j+1/2}{W}\cdot2-1$ and $y_i=\tfrac{i+1/2}{H}\cdot2-1$, downsampled by a factor of four, full batch, 5000 epochs. The metrics include HFEN, defined through a $15\times15$ Laplacian-of-Gaussian filter with $\sigma\approx1.5$ pixels:

$$
\mathrm{HFEN}_{\rm rel}=\frac{\|\hat I^{\rm HP}-I^{\rm HP}\|_2}{\|I^{\rm HP}\|_2} .
$$

Conclusion: NN-CA has lower HFEN and lower relative $L^2$ than RFF-NN on all four images, and higher PSNR.

**Example 3: a demonstration of adaptive frequency enhancement.** The target is $u(x)=\sin(2\pi\cdot2x)+0.5\sin(2\pi\cdot20x+0.3)+0.5\cos(2\pi\cdot40x-0.2)$ on $(0,1)$ — a deliberately sparse spectrum. $N_{\rm train}=2048$, $N_{\rm test}=4096$, $m_{\rm base}=128$, $n_{\rm scales}=1$, $d_q=64$, $L=3$, 4 heads, full-batch Adam at learning rate $10^{-3}$ multiplied by $0.9$ every 500 epochs. Stage one runs 5000 epochs; a discrete Fourier transform at $N_{\rm fft}=4096$ with $\lambda=0.02$ **recovers exactly the dominant modes $k=2,20,40$**. Stage two runs another 5000 epochs with $\eta_{\rm start}=-6$ held constant for the first 70% and released to 0 on a cosine schedule over the last 30%. The result is a visibly faster error decay once the tokens are released.

**Example 4: 1D Poisson.** On $[-1,1]$,

$$
u=\sin(0.1\pi x)+0.2\sin(\pi x)+0.4\sin(\tfrac\nu3\pi x)+0.6\sin(\tfrac{2\nu}{3}\pi x)+\sin(\nu\pi x),
\qquad \nu=100 .
$$

The ablation over $\alpha\in\{0,1,\text{learnable},\text{optimal}\}$ shows that increasing flexibility monotonically lowers the **training loss**, and that a learnable $\alpha$ gives the lowest loss — **but the optimal scaling gives the lowest relative $L^2$ error**. This is an explicit divergence between loss and error, which the authors flag and explain themselves: early in training $\alpha_{\rm opt}$ takes a large negative value, injecting a strong low-frequency correction. Architecture ablation: RFF-CA and NN-CA both decay faster than RFF-NN, with NN-CA best overall. The amplitude ablation deliberately takes a mismatched $\sigma=0.02$ (base frequencies **above** the solution's frequencies): a learnable $\beta$ beats a fixed $\beta\equiv0$ at both boundary-penalty weights $\lambda=10^3,10^4$; and with $\beta=0$ a larger boundary weight actually helps, which the authors explain as "injecting more low-frequency boundary information into the loss".

**Example 5: 2D Poisson.** $u=\sin(\mu x_1^2)+\sin(\mu x_2^2)$, with $N_r=10^4$ interior points and $N_b=1000$ boundary points per side at each iteration, $2\times10^4$ AdamW steps, StepLR $(2000,0.5)$, boundary weight $\lambda=10^4$, the optimal $\alpha$, and the same random seed. At $\mu=50$ both work, with RFF-CA slightly better and more stable; **at $\mu=100$ RFF-NN saturates and stops decreasing while RFF-CA keeps going**. This is the clearest evidence that the benefit grows with frequency.

**Example 6: 3D Poisson–Boltzmann (Deep Ritz).** $-\nabla\cdot(\epsilon\nabla u)+\kappa u=f$ with jump conditions $[u]=0$ and $[\epsilon\,\partial u/\partial n]=0$ across the interface; the exact solution is

$$
u=\frac{e^{\sin\mu x_1+\sin\mu x_2+\sin\mu x_3}}{|x|^2+1}\bigl(|x|^2-1\bigr),
\qquad \mu=15,
$$

with $\epsilon\equiv1$ and $\kappa$ equal to 1 on $\Omega_1$ and 5 on $\Omega_2$. **The domain is constructed to be geometrically singular**: a ball of radius 0.5 at the origin, unioned with 20 small balls of radii drawn from $[0.1,0.2]$ whose centres land randomly on its surface, then truncated by the unit ball — the resulting creases are hard to mesh. The Ritz loss is

$$
L_{\rm Ritz}=\tfrac12\int_\Omega\bigl(|\epsilon\nabla u_\theta|^2+\kappa u_\theta^2\bigr)
-\int_\Omega fu_\theta+\gamma\int_{\partial\Omega}|u_\theta-g|^2,
$$

with $g\equiv0$ and $\gamma=10^4$; $L=3$, $\sigma=1$, learnable $\alpha$, 10,000 epochs, Adam at learning rate $10^{-3}$, StepLR $(1000,0.6)$, and 5000 interior plus 4000 boundary points per epoch. Conclusion: RFF-CA gives a visibly more accurate reconstruction, while "RFF-NN remains far from the true solution".

This set of experiments is controlled — RFF-NN and RFF-CA "share the same multi-scale random Fourier tokeniser and the same model capacity", so the difference can be attributed to cross-attention itself. What they establish is that **learning input-dependent weights over a sufficiently wide fixed basis bank is more robust than hard-coding the frequencies**, and that the benefit widens as frequency grows. What they do not establish is any approximation or optimisation guarantee — the paper says as much — and the applicability of the discrete-Fourier-transform step is limited to periodic or grid-friendly settings.

### Relation to the other papers

It is the direct methodological successor to the multi-scale network used in paper 94: both decompose the spectrum into scales, but paper 94 fixes 100 linear scales and sums the subnetwork outputs, while paper 101 keeps a dyadic bank and **learns where to use which scales** through attention, and can append modes discovered by a discrete Fourier transform. Paper 101 explicitly lists MscaleDNN among the methods that "rely on pre-specified frequency bases".

Drawing frequencies at random and then freezing them is shared with papers 90 and 102 (where the hidden-layer parameters are also drawn and left untrained), but the read-out differs: here it is a trained deep attention stack, there a least-squares solve. Xiaodong Feng, Xiaoliang Wan and Tao Zhou are also joint authors on papers 98 and 105; within this subgroup paper 101 is the deterministic, approximation-theoretic strand and the other two the probabilistic one. And using two networks to handle low and high frequencies with an analytically optimal mixing coefficient is conceptually adjacent to paper 90's separation of smooth from low-regularity subdomains: both isolate the hard part into a dedicated approximation space, one in frequency and one in space.

## 103: make operator splitting the architecture, and put time back into the input

### Intuition

Operator splitting is a classical way to solve evolution equations: advance the linear and nonlinear parts of $u_t=\mathcal Lu+\mathcal Nu$ alternately, each with the method best suited to it. DOSnet's idea is to write that product **directly as a network** — learnable convolution layers play the role of $e^{\tau\mathcal L}$, and **the activation function is the exact flow of the nonlinear subproblem, $e^{\tau\mathcal N}$**, not ReLU or tanh. The activations of the intermediate layers are then physically meaningful intermediate-time states, and the network is no longer a black box.

But DOSnet has two shortcomings. First, it is data-driven, and every training pair for an evolution problem requires a solver run, which is prohibitively expensive. Second, the learned operator is not required to satisfy the equation, so **it cannot be evaluated at arbitrary times** — it outputs only the terminal-time solution, the intermediate blocks' outputs correspond to intermediate times only by comparison against reference data, and the number of recoverable intermediate times is limited by the number of blocks.

This paper's change is small but decisive: take a uniform step $dt=t/K$ and replace the exponential of the linear part with an explicit second-order Taylor expansion in $dt$. The whole block then becomes an **explicit function of $t$**, so $t$ returns to the input, the solution can be evaluated at any time, $\partial\hat u/\partial t$ can be computed by automatic differentiation, and training can switch to a PDE residual with no paired data at all. There is a bonus: at $t=0$ one has $dt=0$, so every block degenerates to the identity and **the initial condition is satisfied exactly and automatically**, with no initial-condition loss term needed.

### Problem setup

Consider the autonomous evolution equation

$$
u_t=\mathcal Lu+\mathcal Nu=:\mathcal Fu\ \text{ in }\Omega\times(0,t^\star],
\qquad u(0,\cdot)=u_0,\qquad \mathcal Bu=0\ \text{ on }\partial\Omega,
$$

with $\mathcal L$ linear, $\mathcal N$ nonlinear and $\mathcal B$ a boundary operator, none of which involves $t$. Partition $[0,t^\star]$ at $0<t_1<\cdots<t_K=t^\star$ with steps $\tau_1,\dots,\tau_K$; Lie–Trotter (first order) and Strang (second order) splitting give respectively

$$
u(T,\bm x)\approx e^{\tau_K\mathcal N}e^{\tau_K\mathcal L}\cdots e^{\tau_1\mathcal N}e^{\tau_1\mathcal L}u(0,\bm x),
$$

$$
u(T,\bm x)\approx e^{\frac{\tau_K}2\mathcal L}e^{\tau_K\mathcal N}e^{\frac{\tau_K}2\mathcal L}
\cdots e^{\frac{\tau_1}2\mathcal L}e^{\tau_1\mathcal N}e^{\frac{\tau_1}2\mathcal L}u(0,\bm x).
$$

The splitting error comes from $\mathcal L$ and $\mathcal N$ not commuting; it vanishes when they do.

DOSnet writes that product directly as a network: $\psi_{\bm\theta_T}=\psi_{\bm\theta_K}\circ\cdots\circ\psi_{\bm\theta_1}$, with each splitting block alternating learnable linear layers (convolutions $\phi_{\mathcal L_{\bm\theta_{l,i}}}$) with nonlinear layers $\phi_{\mathcal N_{l,i}}=e^{\tau_{l,i}\mathcal N}$, subject to $\sum_{l,i}\tau_{l,i}=T$, and trained on data:

$$
\mathfrak L(\bm\theta_T)=\frac1N\sum_n\bigl\|\psi_{\bm\theta_T}(u_0^{(n)})-u^{(n)}(T,\cdot)\bigr\|_{L^2(\Omega)} .
$$

### Derivation

**Step 1: write the block as an explicit function of $t$.** Take uniform steps $\tau_{1,1}=\dots=\tau_{1,K}=dt=t/K$ and replace the exponential of the linear part with an explicit second-order Taylor expansion in $dt$:

$$
e^{dt\,\mathcal L}\approx\mathcal I+dt\,\mathcal L_{\bm\theta_i}+\frac{dt^2}{2}\mathcal L_{\bm\theta_i}^2,
$$

$$
u_i(t_{i-1}+dt,\bm x)=\psi_{\bm\theta_i,dt}(u_{i-1})
=e^{dt\,\mathcal N}\circ\Bigl(\mathcal I+dt\,\mathcal L_{\bm\theta_i}+\tfrac{dt^2}{2}\mathcal L_{\bm\theta_i}^2\Bigr)u_{i-1}(t_{i-1},\bm x),
$$

$$
\hat u(t,\bm x)=\Phi_{\bm\theta_T}(u_0,t)=\psi_{\bm\theta_K,dt}\circ\cdots\circ\psi_{\bm\theta_1,dt}\,u_0(\bm x).
$$

**A neat structural consequence: at $t=0$ one has $dt=0$, so $\psi_{i,dt}=\mathcal I$ and the initial condition holds exactly and automatically** — no initial-condition loss term is needed.

**Step 2: the two kinds of derivative take different routes.** The time derivative $\partial\hat u/\partial t$ is taken by **automatic differentiation**, which is only possible because $t$ enters explicitly through $dt=t/K$; the spatial derivatives are taken by **finite differences**, because $\bm x$ is not a network input (the network maps functions on a fixed grid to functions).

**Step 3: the physics-informed loss.** With the residual $f:=u_t-\mathcal Fu$,

$$
\mathfrak L(\bm\theta)=\lambda_r\mathfrak L_r+\lambda_b\mathfrak L_b,
\qquad
\mathfrak L_r=\frac{1}{N_uN_r}\sum_{n=1}^{N_u}\sum_{j=1}^{N_r}\bigl|f^{(n)}_{r,j}\bigr|^2,
\qquad
\mathfrak L_b=\frac{1}{N_uN_b}\sum_{n=1}^{N_u}\sum_{j=1}^{N_b}\bigl|\mathcal Bu^{(n)}_{b,j}\bigr|^2 .
$$

The spatial collocation points come from a fixed uniform set $\{\bm x_i\}_{i=1}^{N_x}$, and the time points are drawn at random from $(0,T]$. If the boundary condition is structurally guaranteed (periodic, say), $\mathfrak L_b$ is dropped. **There is no data term at all.**

**Step 4: long-time inference by iteration.** Setting $t=t^\star>T$ directly is inaccurate, so for $t^\star\in(mT,(m+1)T]$ the composition is used instead: $\hat u(T)=\Phi_{\bm\theta_T}(u_0,T)$, then $\hat u(2T)=\Phi_{\bm\theta_T}(\hat u(T),T)$, and so on $m+1$ times.

**Step 5: self-diagnosis by residual, to decide when to retrain.** Since there is no data, the paper uses the **physics residual itself** as a validity monitor:

$$
\mathcal R_t=\frac{1}{|\mathfrak I_0|N_x}\sum_{n}\sum_{j}
\Bigl|\partial_t\hat u^{(n)}(t,\bm x_j)-\mathcal F\hat u^{(n)}(t,\bm x_j)\Bigr|^2
+\Bigl|\mathcal B\hat u^{(n)}(t,\bm x_j)\Bigr|^2,
$$

accepting the prediction as long as $\mathcal R_t\le\epsilon\cdot\mathcal R_T$ with $\epsilon\ge1$ set by the user. If this is violated at $t=mT$, draw $N_{\rm add}$ functions at random from $\{\hat u^{(n)}((m-1)T,\cdot)\}$, append them to the initial-condition training set $\mathfrak I$ and retrain. The motivation is stated plainly: oscillations appear when "the subspace spanned by the training initial functions in $\mathfrak I$ is insufficient to characterise a state that has evolved substantially".

### Theorems

**This is one of the few papers on this page with real theorems.**

**Linear stability analysis.** On the test equation $u_t=\mathcal Lu+\lambda u$ with $\mathcal Lu=-q_{\bm\theta}u$, $q_{\bm\theta}>0$ real and $\lambda$ complex, one block gives

$$
u_{i+1}=e^{\lambda dt}\Bigl(1-q_{\bm\theta}dt+\tfrac12q_{\bm\theta}^2dt^2\Bigr)u_i,
$$

with amplification factor

$$
\xi=\Bigl(1-q_{\bm\theta}dt+\tfrac12q_{\bm\theta}^2dt^2\Bigr)
e^{\mathrm{Re}(\lambda dt)+\mathrm i\,\mathrm{Im}(\lambda dt)},
$$

so the stability boundary $|\xi|=1$ is equivalent to

$$
\mathrm{Re}(\lambda dt)=\ln\!\left(\frac{1}{1-q_{\bm\theta}dt+\tfrac12q_{\bm\theta}^2dt^2}\right).
$$

The conclusions drawn from this are: for $0<q_{\bm\theta}dt<1$ the stability region grows with $q_{\bm\theta}dt$, and the intersection of all these regions is the left half-plane; while for $q_{\bm\theta}dt\ge2$ the second-order expansion is **unconditionally unstable** for any $\lambda>0$. If $\max_u\partial(\mathcal Nu)/\partial u\le0$ (that is, $\mathrm{Re}\,\lambda\le0$) then PI-DOSnet is unconditionally stable; otherwise one needs

$$
dt\le\frac{1}{\mathrm{Re}\lambda}\ln\Bigl(1\big/\bigl(1-q_{\bm\theta}dt+\tfrac12q_{\bm\theta}^2dt^2\bigr)\Bigr),
$$

which the paper calls "analogous to a CFL condition", equivalently written $q_{\bm\theta}dt\le C$. **The practically decisive point is that for a learned $\mathcal L_{\bm\theta}$, the eigenvalue $q_{\bm\theta}$ is "far less sensitive to the mesh size $h$" than a finite-difference discretisation** — which is exactly why it can take large $dt$ on fine meshes.

**Theorem 3.1.** Let $\mathcal X$ be a Hilbert space and $\mathcal L:\mathcal X\to\mathcal X$ linear, compact and self-adjoint. Then for any $u\in\mathcal X$ with $\|u\|_{\mathcal X}\le C_0$,

$$
\Bigl\|e^{dt\,\mathcal L}u-\bigl(\mathcal I+dt\,\mathcal L+\tfrac{dt^2}{2}\mathcal L^2\bigr)u\Bigr\|_{\mathcal X}\le C\,dt^3,
$$

with $C$ depending on $\|u\|_{\mathcal X}$ and $\|\mathcal L\|_{\mathcal X\to\mathcal X}$. The proof works by spectral decomposition in the eigenbasis of $\mathcal L$.

**Error decomposition.** The total error of a single block splits into $\epsilon_1+\epsilon_2$, with $\epsilon_1$ the Taylor truncation and $\epsilon_2$ the difference between $\mathcal L$ and $\mathcal L_{\bm\theta}$. Definition 3.1 introduces translation equivariance $\mathcal L(f(x-a))=(\mathcal Lf)(x-a)$; for $\mathcal L=\Delta$ in one dimension with kernel $\tfrac{1}{h^2}[1,-2,1]$, the paper gives

$$
\bigl|[e^{dt\Delta}u](\bm x_j)-[\tilde{\mathcal L}_{\bm\theta}u](\bm x_j)\bigr|
\le C_1\,dt^3+C_2\,dt\,h^2 .
$$

The remark following the theorem says: the $dt^3$ term is independent of the network parameters and shrinks with $dt$; at larger $dt$ (say $dt=1/4$) the error can be mitigated by **enlarging the convolution kernel**. The paper also cites DOSnet for the fact that the approximation error of an operator-splitting block is itself of order $dt^3$, and so does not rederive it.

> [!warning] A gap in strength between abstract and main text
> The abstract says that for the Allen–Cahn equation the method gives "energy-stable solutions even with large time steps", but the main text says the predicted energy shows a "nearly monotone" decrease and that "the energy dissipation law is not strictly enforced". So "energy stability" is an empirical observation, not a proved property; the authors also say that designing a network satisfying energy dissipation strictly "remains a challenging task, left for the future".

### Numerical experiments

The metrics are relative $L^2$ error at a given time and cumulative relative $L^2$ error over $[0,\tilde t]$; the implementation uses PyTorch and Adam.

**Example 1: 1D convection.** $u_t+\beta u_x=0$ on $(-\pi,\pi)$ with periodic boundaries; initial data $u_0(x)=\sum_{i=1}^{10}(c_i\sin ix+q_i\cos ix)$ with $c_i,q_i\sim\mathcal N(0,1)$; 300 uniform spatial points; 4 blocks, channels 1-4-1, kernel size 61, $dt=0.075$; trained at $T=0.3$ for 20,000 epochs. A single forward pass predicts $2T,\dots,10T$; the absolute error at $10T$ is below 0.05, and the space-time-averaged relative $L^2$ error is $3.149\times10^{-3}$.

**Taylor order ablation** (first against second order):

| $t$ | first order          | second order         |
| --- | -------------------- | -------------------- |
| 3   | $1.640\times10^{-2}$ | $5.318\times10^{-3}$ |
| 6   | $3.283\times10^{-2}$ | $1.062\times10^{-2}$ |
| 9   | $4.927\times10^{-2}$ | $1.592\times10^{-2}$ |
| 15  | $8.220\times10^{-2}$ | $2.650\times10^{-2}$ |
| 21  | $1.151\times10^{-1}$ | $3.704\times10^{-2}$ |
| 30  | $2.350\times10^{-1}$ | $5.278\times10^{-2}$ |

Roughly a factor of 3 to 4, widening with time. **This is the experimental justification for the second-order expansion.**

**Example 2: diffusion-reaction.** $u_t=Du_{xx}+ku^2$ on $(0,1)$ with $D=k=0.001$, homogeneous Dirichlet conditions and $T_{\rm end}=50$, with training windows $T\in\{0.5,1,2\}$. Four blocks with 3 convolution layers each, channels 1-8-8-1, kernel 21, for **only 6720 trainable parameters in total**, and 101 spatial points.

**Example 3: Allen–Cahn.** $u_t-0.0004\Delta u+(u^3-u)=0$ on $(-1,1)^d$ with periodic boundaries, $d=1,2$; initial data of random Fourier type (7 modes in 1D; an $8\times8$ double sum times 0.02 in 2D); reference solutions from Chebfun. In 1D: $T=1$, $T_{\rm end}=10$, $200\times20$ collocation points, 4 blocks, channels 1-4-1, kernel 15, 8000 epochs, 1000 initial functions. Inference breaks down at $t=3$; appending 100 solutions at $t=2$ and retraining reaches $t=6$; another round covers $[0,10T]$. After two retrainings the mean squared error drops to about $4\times10^{-4}$ with absolute error below 0.04. The energy $E=\tfrac{0.0004}{2}\int u_x^2+\int\tfrac{(u^2-1)^2}{4}$ decreases nearly monotonically and tracks the Chebfun reference.

**Mesh refinement and stability — the strongest result in the paper:**

| $N_x$ | PI-DOSnet (learned convolution) | second-order central-difference kernel | fourth-order central-difference kernel |
| ----- | ------------------------------- | -------------------------------------- | -------------------------------------- |
| 200   | 4 blocks, 0.07 s                | 8 blocks, 0.08 s                       | 11 blocks, 0.10 s                      |
| 400   | 4 blocks, 0.07 s                | 32 blocks, 0.26 s                      | 42 blocks, 0.34 s                      |
| 800   | 4 blocks, 0.07 s                | 128 blocks, 1.00 s                     | 170 blocks, 1.35 s                     |

PI-DOSnet uses **4 blocks** on all three meshes, with $L^2$ error at $t=1$ of $1.81\times10^{-3}$, $1.97\times10^{-3}$ and $1.11\times10^{-3}$, and at $t=10$ of $9.86\times10^{-3}$, $1.19\times10^{-2}$ and $9.57\times10^{-3}$. The mechanism is confirmed by the spectral radius: that of $-\mathcal L_{\bm\theta}$ grows "only slowly" with $N_x$, while that of the second-order central-difference $-\mathcal L$ grows "quadratically", imposing a CFL restriction. In 2D: $T=1$, $T_{\rm end}=10$, 20 time instants, a $200\times200$ grid, 4 blocks, 4 channels, $15\times15$ kernels, 400 initial functions; inference again stops at $t=3$, with 40 predictions at $t=2$ appended and then predictions at $t=5$ appended.

**Example 4: the Gross–Pitaevskii equation.** $\mathrm i\partial_t\psi=[-\tfrac{\partial_{xx}}{2m}+V(x)+g(x)|\psi|^2]\psi$ (complex-valued). Solitary waves are identified by minimising $|F-1|$ where

$$
F=\frac{1}{N^2}\Bigl|\int\psi^*(T_{\rm end},x)\psi(0,x)\,\mathrm dx\Bigr|^2,
\qquad N=\int|\psi(0,x)|^2\mathrm dx .
$$

$1024\times10$ collocation points, training interval $[0,0.1]$, $\tilde T_{\rm end}=1$, 8 blocks of depth 2, kernel 21, 15,000 epochs, 1000 initial functions; the reference solution comes from a time-splitting spectral method. Over 100 random test initial conditions, the mean relative $L^2$ error is $5.325\times10^{-3}$ at $\tilde T=0.1$ and $4.087\times10^{-2}$ at $\tilde T_{\rm end}=1$.

What these experiments establish is that **writing the splitting structure into the architecture and letting $t$ enter the input explicitly buys three things at once** — data-free training, evaluation at arbitrary times, and a time step that does not degrade under mesh refinement. What they do not establish is the rigour of energy dissipation (see the warning above), and they do not cover the three open problems the paper itself lists: complex geometries, high-dimensional problems and low-regularity solutions.

### Relation to the other papers

"Low-regularity solutions", listed as an open problem in the summary, is exactly paper 90's target — the two are complementary halves of the same research programme, discussed on the [[en/computational-mathematics/paper-notes/scientific-machine-learning/variational-and-basis-networks|variational and basis networks page]]. It shares the physics-informed operator-learning setting with paper 94 (physics-informed DeepONet with Gaussian wave packets) and paper 95 (UQ-SONet): all three learn the map from initial data to solution, with 94 and 103 data-free and 95 uncertainty-aware. The architectural difference is that PI-DOSnet's inductive bias **is the splitting scheme itself**, with $e^{dt\mathcal N}$ as the activation, rather than a branch/trunk factorisation.

Jizu Huang and Tao Zhou are also joint authors on paper 94; the operator splitting PI-DOSnet uses is precisely the classical tool whose $O(\Delta t^2/\varepsilon)$ failure in paper 94 motivated the wave-packet reformulation. Decoupling a coupled system sequentially into cheaper subproblems is the temporal counterpart of paper 102's spatial decoupling. And the residual-triggered retraining loop is the self-diagnostic version of paper 90's residual-triggered subdomain refinement and paper 101's discrete-Fourier-transform-triggered token injection: all three use a computable indicator on the current approximation to decide where to add capacity.

## 105: a shared summary network that ties filtering and smoothing together

### Intuition

In a state-space model, filtering wants $p(u_t\mid y_{1:t})$ and smoothing wants $p(u_{1:t}\mid y_{1:t})$. Both condition on the **entire observation history**, whose dimension grows with $t$, so a fixed-input normalizing flow cannot be used. The fix is to compress the history into a fixed-length summary $s_t$ with a recurrent network and then condition the flow on $s_t$.

Up to here this is routine. The paper's real move is to **let the filtering flow and the smoothing flow share the same summary network**. There are two reasons. The first is algebraic: smoothing has an exact causal factorisation in which the backward kernel conditions on $y_{1:k}$ rather than $y_{1:t}$, so the same $s_k$ can be reused directly; and at one particular weight the two loss terms recombine into "per-time-step filtering likelihood" plus "whole-trajectory smoothing likelihood", so this **single** loss is simultaneously the maximum likelihood for both tasks. The second is sufficiency: if a sufficient summary exists, it is simultaneously optimal for both objectives.

The second reason is only a sufficient condition, not a necessity. The paper's own lemma notes that the optimal summaries for the two objectives generally **do not coincide**, so sharing is not free. And its ablation gives the answer from the other direction: without sharing, each flow individually is nearly unaffected, but backward **iterative sampling** collapses entirely. That is, sharing here is not beneficial but necessary.

### Problem setup

The state-space model

$$
u_{t}=f(u_{t-1},\epsilon_{u,t}),\qquad y_{t}=h(u_{t},\epsilon_{y,t}),
$$

with $u_t\in\mathbb R^{n_u}$, $y_t\in\mathbb R^{n_y}$ and $\epsilon_{u,t}$, $\epsilon_{y,t}$ mutually independent and i.i.d. in time, induces the Markov factorisation

$$
p(u_{1:t},y_{1:t})=p(u_{1})p(y_{1}\mid u_{1})\prod_{k=2}^{t}p(u_{k}\mid u_{k-1})p(y_{k}\mid u_{k}).
$$

The training data are $N$ simulated trajectories $\{u^i_{1:T},y^i_{1:T}\}_{i=1}^{N}$.

The paper's charge against existing work has three parts: Gaussian filters (EKF, EnKF) rely on "moment-closure assumptions up to the first two moments, which are systematically biased for general nonlinear systems"; sequential Monte Carlo and particle filters "typically suffer severe weight degeneracy as state dimension and time horizon grow", and worse for smoothing; and recent deep generative filters "usually treat filtering and smoothing separately and often rely on expensive per-instance optimisation", making them unsuitable for online and large-scale deployment. **The paper constrains itself more tightly than usual**: it assumes only a **simulator**, with no assumption that the functional forms of $f$, $h$ and the noise distributions are known, and no requirement that the transition and observation densities be computable.

A conditional normalizing flow is defined by: for a conditioning variable $s$ and reference $Z\sim p_Z=\mathcal N(0,I)$,

$$
z=f_{\theta}(u;s),\qquad
p_{\theta}(u\mid s)=p_{Z}\bigl(f_{\theta}(u;s)\bigr)\bigl|\det\nabla_{u}f_{\theta}(u;s)\bigr|,
$$

with $f_\theta=f_{[L]}\circ\cdots\circ f_{[1]}$, the determinant factorising across layers, and training by maximum likelihood (equivalently, minimising the conditional KL averaged over $s$).

### Derivation

**Step 1: compress the history into a fixed-length summary.** A multi-layer LSTM is used; a single layer written out in full is

$$
i_{t}=\sigma(W_{yi}y_{t}+W_{hi}h_{t-1}+b_{i}),\qquad
f_{t}=\sigma(W_{yf}y_{t}+W_{hf}h_{t-1}+b_{f}),
$$

$$
o_{t}=\sigma(W_{yo}y_{t}+W_{ho}h_{t-1}+b_{o}),\qquad
g_{t}=\tanh(W_{yg}y_{t}+W_{hg}h_{t-1}+b_{g}),
$$

$$
c_{t}=f_{t}\odot c_{t-1}+i_{t}\odot g_{t},\qquad
h_{t}=o_{t}\odot\tanh(c_{t}),\qquad (h_0,c_0)=(\bm 0,\bm 0),
$$

stacked $L$ deep ($z^{(1)}_t:=y_t$, $z^{(\ell)}_t:=h^{(\ell-1)}_t$), with the summary read off the top layer through an affine map:

$$
s_{t}=W_{s}h_{t}^{(L)}+b_{s},\qquad W_{s}\in\mathbb R^{d_{s}\times d_{h}} .
$$

Remark 2.1 gives the reason for choosing an LSTM over a plain RNN: unrolling the additive cell update,

$$
c_{t}=\Bigl(\prod_{j=1}^{t}f_{j}\Bigr)\odot c_{0}
+\sum_{s=1}^{t}\Bigl(\prod_{j=s+1}^{t}f_{j}\Bigr)\odot(i_{s}\odot g_{s}),
$$

so $\partial c_t/\partial c_{t-k}=\prod_{j=t-k+1}^{t}f_j$, and "the gradient along the cell state is governed by a product of gates rather than by repeated application of the same linear map and nonlinearity".

**Step 2: the two flows.** The forward flow conditions directly on the summary:

$$
p(u_{t}\mid y_{1:t})\approx p_{\theta_1,\psi}(u_{t}\mid s_{t}),
\qquad s_t=\mathrm{Enc}(y_{1:t};\psi)\in\mathbb R^{h} .
$$

The smoothing side uses an **exact factorisation**

$$
p(u_{1:t}\mid y_{1:t})=p(u_{t}\mid y_{1:t})\prod_{k=1}^{t-1}p(u_{k}\mid u_{k+1},y_{1:k}),
$$

noting that the second factor conditions on $y_{1:k}$ and **not** on $y_{1:t}$ — this is exactly what keeps the recursion causal, so the same summary $s_k$ can be reused directly. The second flow learns

$$
p(u_{t}\mid u_{t+1},y_{1:t})\approx p_{\theta_2,\psi}(u_{t}\mid u_{t+1},s_{t}).
$$

**Both flows share the same summary network $\psi$, which is the paper's namesake contribution.**

**Step 3: the joint objective and that algebraic identity.**

$$
\min_{\theta_{1},\theta_{2},\psi}\
-\frac1{NT}\sum_{i=1}^N\sum_{t=1}^T\log p_{\theta_1,\psi}(u_t^i\mid s_t^i)
-\frac{\lambda}{N(T-1)}\sum_{i=1}^N\sum_{t=1}^{T-1}\log p_{\theta_2,\psi}(u_t^i\mid u_{t+1}^i,s_t^i).
$$

At $\lambda=(T-1)/T$ the two terms recombine. Defining the induced trajectory density

$$
p^{\mathrm{smoothing}}_{\theta_{1},\theta_{2},\psi}(u_{1:T}\mid s_{1:T})
:=p_{\theta_{1},\psi}(u_{T}\mid s_{T})\prod_{t=1}^{T-1}p_{\theta_{2},\psi}(u_{t}\mid u_{t+1},s_{t}),
$$

the loss becomes exactly

$$
\mathcal L=-\frac{1}{NT}\sum_{i=1}^{N}\sum_{t=1}^{T-1}\log p_{\theta_{1},\psi}(u^{i}_{t}\mid s^{i}_{t})
-\frac{1}{NT}\sum_{i=1}^{N}\log p^{\mathrm{smoothing}}_{\theta_{1},\theta_{2},\psi}(u^{i}_{1:T}\mid s^{i}_{1:T}),
$$

that is, "marginal terms for filtering at intermediate times" plus "a path term for the induced trajectory model". Because both condition on the same set $\{s_t\}$, "the learned backward transitions must remain consistent with the terminal filtering distribution under a common representation of the observation history". **This is the implicit consistency regularisation mentioned in the abstract, and it is an exact rewriting rather than a heuristic.**

**Step 4: details of the conditional KRnet.** With $z=T_\Theta(u;c)$ and $\log p_\Theta(u\mid c)=\log p_Z(T_\Theta(u;c))+\log|\det\nabla_uT_\Theta(u;c)|$, there are three kinds of module. A conditional scale-bias layer $T^{\mathrm s}_{\psi}(u;c)=\exp(\eta_{\psi}(c))\odot u+\xi_{\psi}(c)$ with log-determinant $\sum_j\eta_{\psi,j}(c)$; and a conditional affine coupling layer that partitions $u=[u^{(1)},u^{(2)}]^\top$ with $u^{(1)}\in\mathbb R^{k}$, $k=\lfloor n_u/2\rfloor$,

$$
\tilde u^{(1)}=u^{(1)},
\qquad
\tilde u^{(2)}=\bigl(\bm 1+\alpha\tanh s_{\omega}(u^{(1)},c)\bigr)\odot u^{(2)}
+\gamma\odot\tanh t_{\omega}(u^{(1)},c),
$$

where **$\alpha=0.6$ is fixed** and $\gamma$ is a positive learnable vector; since $\tanh\in(-1,1)$ the scale lies in $(1-\alpha,1+\alpha)$, which guarantees invertibility, and the log-determinant is $\sum_j\log(1+\alpha\tanh s_{\omega,j})$. The coupling networks are **random Fourier feature networks**: starting from $h_0=[u^{(1)},c]^\top$,

$$
h_{1}=\begin{bmatrix}\sin(e^{-\sigma}Fh_{0}+b_{0})\\ \cos(e^{-\sigma}Fh_{0}+b_{0})\\ h_{0}\end{bmatrix},
\qquad
h_{i}=\mathrm{SiLU}(W_{i-1}h_{i-1}+b_{i-1}),
$$

with $(F,b_0)$ **fixed random features** and the bandwidth-like scalar $\sigma$ and $\{W_i,b_i\}$ trainable. The whole composition is $T_{\Theta}=\Pi_{K}\circ T^{\mathrm{coup}}_{\omega_{K}}\circ\cdots\circ\Pi_{1}\circ T^{\mathrm{coup}}_{\omega_{1}}\circ T^{\mathrm{s}}_{\psi}$, with $\Pi_k$ fixed coordinate permutations of determinant 1.

**Step 5: a flow-based particle filter variant.** Instead of a fully amortised marginal, use the **fully adapted proposal** factorisation

$$
p(u_{k},y_{k}\mid u_{k-1})=p(y_{k}\mid u_{k-1})\,p(u_{k}\mid y_{k},u_{k-1}),
$$

and learn both factors as conditional flows: $p_{\theta_3}(y_k\mid u_{k-1})$ and $p_{\theta_4}(u_k\mid y_k,u_{k-1})$, each trained by its own negative log-likelihood on simulated triples. One filtering step is: predictive weighting $\tilde\omega^{(j)}_{k-1}\propto\omega^{(j)}_{k-1}p_{\theta_3}(y_k\mid u^{(j)}_{k-1})$; resampling $a^{(j)}\sim\mathrm{Cat}(\tilde\omega^{(1:N)}_{k-1})$; and propagation $u^{(j)}_{k}\sim p_{\theta_4}(\cdot\mid y_{k},u^{(a^{(j)})}_{k-1})$. Remark 3.2 contrasts this with a bootstrap-type alternative ($p_{\theta_5}(u_k\mid u_{k-1})$, $p_{\theta_6}(y_k\mid u_k)$): the adapted proposal sees $y_k$ and is therefore "expected to generate more informative particles".

**Step 6: an effective-sample-size diagnostic (available only when the true transition density and likelihood are).** Comparing

$$
p(u_t,y_t,u_{t-1}\mid y_{1:t-1})=p(u_t\mid u_{t-1})p(y_t\mid u_t)p(u_{t-1}\mid y_{1:t-1})
$$

against the learned $q=p_{\theta_4}(u_t\mid y_t,u_{t-1})p_{\theta_3}(y_t\mid u_{t-1})p(u_{t-1}\mid y_{1:t-1})$, the weights are

$$
\omega_{i}=\frac{p(u^{(i)}_{t}\mid u^{(i)}_{t-1})p(y^{(i)}_{t}\mid u^{(i)}_{t})}
{p_{\theta_{4}}(u^{(i)}_{t}\mid y^{(i)}_{t},u^{(i)}_{t-1})p_{\theta_{3}}(y^{(i)}_{t}\mid u^{(i)}_{t-1})},
\qquad
\mathrm{ESS}=\frac{(\sum_i\omega_i)^2}{\sum_i\omega_i^2},
\qquad
\mathrm{RESS}=\frac{\mathrm{ESS}}{N},
$$

connected to the chi-squared divergence through $\chi^{2}(p\|q)=\int p^2/q\,\mathrm dx-1\approx N/\mathrm{ESS}-1$. Remark 3.3 gives the ancestral sampling order for $q$ and notes that $N=10^{6}$ is needed for a stable Monte Carlo estimate.

**Inference.** Filtering (Algorithm 2) only needs $s_t$ and then samples from $p_{\theta_1,\psi}(\cdot\mid s_t)$; it is **amortised**, with no per-instance optimisation at test time. Smoothing (Algorithm 3) first computes $s_1,\dots,s_t$ and draws $\{u^j_t\}$ from the terminal filtering distribution, then for $k=t-1,\dots,1$ draws $u^j_k\sim p_{\theta_2,\psi}(\cdot\mid u^j_{k+1},s_k)$ and prepends it to the path; stopping early at any $k$ gives marginal smoothing samples from $p(u_k\mid y_{1:t})$.

### Theorems

**The theory consists of one lemma and one proposition in Appendix B, both about the shared summary statistic.** Write $X:=u_k$, $Y:=y_{1:k}$, $Z:=u_{k+1}$, and let $S=f(Y)$ be any measurable summary. The two objectives $\min\mathbb E[-\log p_{\theta_1}(X\mid f(Y))]$ and $\min\mathbb E[-\log p_{\theta_2}(X\mid Z,f(Y))]$ attain, at the optimal flow parameters, the conditional entropies $H(X\mid S)$ and $H(X\mid Z,S)$ respectively. Write $S^*_F=\arg\min_{S=f(Y)}H(X\mid S)$ and $S^*_S=\arg\min_{S=f(Y)}H(X\mid Z,S)$.

**Lemma B.1.** $I(X,S^{*}_{F})\ge I(X,S^{*}_{S})$ and $I(X;Z\mid S^{*}_{S})\ge I(X;Z\mid S^{*}_{F})$. (The proof uses $H(X\mid S)=H(X\mid Z,S)+I(X;Z\mid S)$ and the two minimality properties.) The reading: **the optimal summaries for filtering and for smoothing need not coincide — so sharing is not free.**

**Definition B.1 and Proposition B.1.** A summary $S^{\dagger}=f(Y)$ is called **sufficient** if $p(X\mid Y)=p(X\mid S^{\dagger})$ almost surely. If a sufficient summary exists, then $S^{\dagger}=S^{*}_{F}=S^{*}_{S}$, so a single shared summary network is simultaneously optimal for both objectives.

**There is no convergence, consistency, error or stability theorem for FLUID itself.** The conclusion lists "establishing rigorous theoretical guarantees" as future work.

### Numerical experiments

Shared configuration: both flows are one conditional scale-bias layer plus **6 conditional affine coupling layers**, with each coupling network a random Fourier feature embedding followed by a multilayer perceptron of **depth 6 and width 64**; the summary network is a **4-layer LSTM** (reduced to **1 layer** for the linear advection-diffusion example); $\dim(s_t)=3\dim(y_t)$, or $5\dim(y_t)$ for the stochastic volatility example; Adam with initial learning rate $0.001$; PyTorch on a single **NVIDIA V100 (32 GB)**. The metrics are KL (only where the exact posterior is available), RMSE, MMD (Gaussian kernel with $\sigma=2$, against the delta measure at the true state) and CRPS.

**Example 1: 1D linear advection-diffusion.** $\partial_tu=a\partial_xu+\kappa\partial_{xx}u$ on $[0,1]$ with periodic boundaries and $u(0,x)=\sin(2\pi x)$, discretised as $u_k=Mu_{k-1}+\epsilon_{u,k}$, $y_k=Hu_k+\epsilon_{y,k}$, $u_0\sim\mathcal N(\mu,\sigma^2I)$ — a **linear Gaussian model whose filtering and smoothing laws are both available in closed form**, which is exactly why the KL column can exist.

Case 1 (pure advection with $a=-1$, $\kappa=0$, upwind $M_{\delta t}=I-\nu A$ with $\nu=\delta t/\Delta x$, $Q=qI_n$, $H$ subsampling every other point so $n_y=n/2$, $R=rI$): $\Delta t=0.05$, $q=0.01$, $r=0.1$, $\sigma=0.05$, all **held fixed** as $n$ grows, so the problem genuinely gets harder; $n=10,20,30,40,50$; $N_{\rm train}=2000$ trajectories of length $T_{\rm train}=500$, $N_{\rm test}=200$.

| $n$ | filtering KL | backward-kernel KL | filtering RMSE | smoothing RMSE | MMD    | CRPS   |
| --- | ------------ | ------------------ | -------------- | -------------- | ------ | ------ |
| 10  | 0.0191       | 0.0146             | 0.1457         | 0.1284         | 0.0514 | 0.0822 |
| 50  | 0.0597       | 0.0534             | 0.1228         | 0.1134         | 0.1710 | 0.0693 |

Note that filtering RMSE and CRPS actually **decrease** as $n$ grows, with only MMD rising; the paper flags the MMD growth itself and argues the remaining metrics stay stable. The $n=50$ case is run to $T=1000$ (physical time $t=50$), **twice the training horizon**, and consistency is maintained beyond the dashed line.

Case 2 (the discretisation-consistent regime with $a=1$, $\kappa=0.01$, a Lax–Wendroff-type $M_{\delta t}$, $Q_{\delta t}=\frac{\delta t}{n}I_n$, and block-averaged observations $y^i_k=\frac{1}{|\mathcal I_i|}\sum_{j\in\mathcal I_i}u^j_k+\epsilon^i_{y,k}$ with $n_{\mathcal I}=8$ held fixed): $\Delta t=0.01$, $r=0.01$, $\sigma=0.05/n$; $n=16,32,48,64$. The filtering KL stays in a narrow band ($0.1355$, $0.1363$, $0.1135$, $0.1268$), while RMSE, MMD and CRPS all **decrease** with $n$ (RMSE from $0.0558$ to $0.0214$) — exactly what mesh-consistent behaviour should look like.

**Example 2: a two-factor stochastic volatility model** (the only example touching real data). $\bm u_t=\bm\alpha+A(\bm u_{t-1}-\bm\alpha)+D_\sigma\epsilon_{u,t}$, $\bm y_t=\beta\exp(\tfrac12\bm u_t)\odot\epsilon_{y,t}$, $\bm u_0\sim\mathcal N(0,\operatorname{diag}(\tau_1^2,\tau_2^2))$, with $\bm\alpha=0$, $\gamma_1=\gamma_2=0.97$, $\sigma_1=\sigma_2=0.3$, $\beta=0.835$ (MAP estimates for the S&P 500 taken from the literature) and $\tau_i^2=\sigma_i^2/(1-\gamma_i^2)$. $N=2000$ trajectories of length $T_{\rm train}=1000$, $N_{\rm test}=200$.

| Method                     | RMSE   | MMD    | CRPS   |
| -------------------------- | ------ | ------ | ------ |
| FLUID                      | 0.6117 | 0.1571 | 0.3435 |
| flow-based particle filter | 0.6163 | 0.1591 | 0.3462 |
| FBF                        | 0.7481 | 0.2163 | 0.4188 |

The backward kernel gives $0.2746/0.0363/0.1551$ and smoothing gives $0.4805/0.1038/0.2710$, both large improvements over filtering. The test trajectories run to $T=2000$, twice the training horizon. **Real-data application**: the trained particle filter is applied to daily S&P 500 returns from 31 December 2018 to 29 December 2022 (the bivariate series is the processed S&P 500 returns paired with a synthetic single-factor trajectory), and the recovered latent volatility has a peak near $k\approx300$ with a simultaneously widening 90% credible interval (the early-2020 pandemic shock), and stays persistently elevated from $k\approx800$ onward (the turbulence of 2022). This is an **interpretive** rather than a quantitative conclusion.

**Example 3: the stochastic Burgers equation.** $du+(u\partial_xu-\nu\partial^2_{xx}u)dt=\sigma\,dW(t)$ with $x\in[-1,1]$, $t\in[0,1]$, $u(0,x)=-\sin(\pi x)$, $u(t,\pm1)=0$, $\sigma=1.0$ and $\nu=0.05$ (Appendix E.1 repeats it at $\nu=0.01$). A mixed finite-difference scheme on a $201\times50$ grid; the state is 50 spatial values, observations subsample every other point ($n_y=25$) with $\mathcal N(0,r^2I)$ noise, $\Delta t=0.005$; $N=3000$ training trajectories of length $T_{\rm train}=200$, $N_{\rm test}=200$.

| $r^2$ | FLUID filtering RMSE | FBF filtering RMSE |
| ----- | -------------------- | ------------------ |
| 0.01  | 0.0751               | 0.0752             |
| 0.04  | 0.0898               | 0.0917             |
| 0.09  | 0.1003               | 0.1015             |
| 0.16  | 0.1084               | 0.1113             |
| 0.25  | 0.1149               | 0.1179             |

The paper's own reading is that the two are "extremely close" with a "small but consistent" advantage to FLUID, and that both are near the accuracy ceiling for this problem. The backward kernel is far sharper than filtering (RMSE $0.0584$ against $0.1149$ at $r^2=0.25$). Appendix E.1 repeats the same sweep at $\nu=0.01$ with the ordering unchanged (for example FLUID $0.1546$ against FBF $0.1640$ at $r^2=0.25$).

**Example 4: Lorenz-96 and two-scale Lorenz.** Parameters $J=32$, $F\in\{5,8,16\}$, $h=1$, $b=10$, $c=4$, $\sigma_u=0.1$, $\sigma_v=0.01$; the observations are **cubic**, $y_{k,i}=u_{k,i}^{3}+\epsilon_{y,i}$ with $\epsilon_{y,i}\sim\mathcal N(0,1)$, which is strongly nonlinear.

Single-scale ($c=0$, $\sigma_u=1$, $F=8$, all indices observed, $u_{0,j}=\sin(2\pi j/n)$, $\Delta t=0.05$, $K=10,\dots,50$):

| $K$ | FLUID filtering RMSE | FBF filtering RMSE |
| --- | -------------------- | ------------------ |
| 10  | 0.1632               | 0.2044             |
| 20  | 0.1945               | 0.2439             |
| 30  | 0.2081               | 0.2604             |
| 40  | 0.2255               | 0.2742             |
| 50  | 0.2605               | 0.3106             |

FLUID is better at every dimension, and the MMD gap widens ($0.5508$ against $0.6918$ at $K=50$).

Two-scale ($K=16$, $J=32$, only odd indices observed, $u_{0,i}=F+\sigma_u\epsilon_u$, $v_{0,j}=\sigma_v\epsilon_v$): here **marginalising out the fast variables makes the effective transition density of the slow variables intractable**, so the model violates the standard form FBF requires.

| $F$ | FLUID filtering RMSE | FBF filtering RMSE | smoothing RMSE |
| --- | -------------------- | ------------------ | -------------- |
| 5   | 0.4544               | 0.5817             | 0.4063         |
| 8   | 0.4397               | 0.8301             | 0.3795         |
| 16  | 0.5218               | **3.7536**         | 0.4837         |

FBF's failure at $F=16$ (CRPS $1.3578$ against $0.2350$) is the strongest comparative evidence in the paper.

### Sharing is not beneficial, it is necessary

The most worthwhile result to record is the ablation. At state dimension $K=20$, shared against independent summaries gives:

| Metric               | shared summary | independent summaries |
| -------------------- | -------------- | --------------------- |
| filtering RMSE       | 0.1945         | 0.1958                |
| backward-kernel RMSE | 0.1240         | 0.1596                |
| smoothing RMSE       | 0.1525         | **60.3744**           |
| smoothing CRPS       | 0.0742         | **33.2893**           |

Each flow **individually** is nearly unaffected, while backward **iterative sampling** collapses entirely. The paper explains this as errors accumulating rapidly through the smoothing recursion under independent summaries, with results already unusable at $K\ge20$. This is the experimental counterpart of the loss-recombination identity above, and it is where the paper's real contribution lies.

> [!warning] Two disagreements between prose and tables
> The main text claims smoothing is uniformly better than filtering, but Table 9 at $K=50$ gives smoothing RMSE 0.5423 against filtering 0.2605, and Table 6 for the Burgers example has the same inversion at $r^2=0.25$ (smoothing 0.1185 against filtering 0.1149). This page records the table values. The clear degradation at $K=50$ is not commented on in the text.

> [!warning] Submission status not verified
> The arXiv record carries no journal-reference field, so "submitted to CMAME" is **not verified by the preprint itself**.

### Relation to the other papers

The technically closest is paper 107: both use a **conditional KRnet** with **random Fourier feature** coupling networks, both learn a family of conditional densities indexed by a conditioning variable, and both cite the same temporal-normalizing-flow lineage. The difference is what conditions the flow and where the training signal comes from — paper 107 conditions on the initial state and time and trains on the **PDE residual** (physics-informed, no samples), while paper 105 conditions on an LSTM summary of the observation history and trains by **maximum likelihood** on simulated trajectories (simulation-based, no PDE). Together they bracket the two routes to fitting conditional normalizing flows.

Xiaodong Feng, Xiaoliang Wan and Tao Zhou are also joint authors on papers 101 and 98. The random Fourier feature embedding in the coupling networks is the same primitive paper 101 uses to build its multi-scale bank — in paper 105 it is a fixed featuriser with a single **trainable** bandwidth scalar $\sigma$, in paper 101 a deliberately multi-scale bank with adaptive amplification.

It belongs to the group's uncertainty-quantification strand along with papers 98 and 95, but with different probabilistic tools: those two use latent-variable models with Gaussian decoders and evidence-lower-bound or information-bottleneck objectives (giving a **bound** on the likelihood), while paper 105 uses invertible flows giving the **exact** likelihood, hence a plain maximum-likelihood objective, and importance weights that can be computed exactly when the model factors are known. It shares with paper 89 the status of a sampling/generative-model paper in this group, with both using KL-based training and importance reweighting with an effective-sample-size diagnostic — paper 89 reports relative effective sample size for Boltzmann sampling and paper 105 reports RESS for the particle filter; but paper 89 faces a **known unnormalised energy** with no data, while paper 105 faces an **unknown posterior** reachable only through simulation.

It is mechanically unrelated to the deterministic least-squares PDE solvers (papers 90, 102, 94, 103): those solve a forward problem to high accuracy with no notion of a posterior, while paper 105 does not solve a PDE at all — the advection-diffusion and Burgers systems appear here only as data generators.

## 89: sampling from a known unnormalised energy

### Intuition

The goal is to sample from $\pi(x)=Z^{-1}\exp(-U(x))$: the energy $U$ can be evaluated pointwise, the normalising constant $Z$ is unavailable, and **there is no training data from $\pi$**. That last condition rules out the usual generative-model objectives — Jensen–Shannon divergence, MMD and Wasserstein distance all need true samples.

That leaves the reverse KL. But the reverse KL needs the generative model's density $\log p_D(x)$, which is where the trouble starts: getting a closed-form density requires the decoder to be a **bijection** (the Boltzmann-generator route), and bijectivity limits effective expressiveness. The other route is diffusion-type samplers, which do not require bijectivity but need a numerical SDE or ODE solver inside the training loop, which is expensive.

This paper removes both restrictions at once. The decoder is an ordinary (**non-bijective**) Gaussian network, and a forward diffusion is then run in **latent space**. That diffusion has two benefits. First, its transition density $p_D(z_t\mid z_0)$ is available analytically, so $z_t$ can be sampled exactly during training with **no differential equation solved at all**. Second, its time-reversal provides an object that can serve as an "encoder", which is used to enlarge the reverse KL into an augmented KL over latent paths — and after augmentation $\log p_D(x)$ itself is no longer needed. Within this bound $\log Z$ is just a parameter-independent constant, so optimisation is possible without knowing it.

### Problem setup

The target is

$$
\pi(x)=\frac1Z\exp(-U(x)),\qquad Z=\int\exp(-U(x))\,\mathrm dx,\quad U:\mathbb R^d\to\mathbb R .
$$

The decoder takes $x\mid z_0\sim p_D(x\mid z_0;\phi)$ with $z_0\sim p_D(z_0)$ a known prior (standard multivariate normal), and the conditional distribution is Gaussian,

$$
p_D(x\mid z_0;\phi)=\mathcal N\bigl(x\mid\mu(z_0;\phi),\,\Sigma(z_0;\phi)\bigr),
$$

with $\mu$ and $\Sigma$ networks and **no invertibility imposed**.

### Derivation

**Step 1: the augmented variational bound.** Since samples from $\pi$ are unavailable, only the reverse KL is usable, and it is enlarged into an augmented KL with latent variables:

$$
D_{\rm KL}\bigl(p_D(x)\,\|\,\pi(x)\bigr)
\le
D_{\rm KL}\bigl(p_{D}(z_{0})p_{D}(x\mid z_{0};\phi)\,\big\|\,\pi(x)p_{E}(z_{0}\mid x;\theta)\bigr)
$$

$$
=\mathbb E_{p_D(z_0)p_D(x\mid z_0;\phi)}
\left[\log\frac{p_D(z_0)p_D(x\mid z_0;\phi)}{p_E(z_0\mid x;\theta)}+U(x)\right]+\log Z,
$$

with equality when $p_E(z_0\mid x;\theta)$ matches the conditional distribution of $z_0\mid x$ induced by the decoder.

**Step 2: the decoding process = decoder + forward diffusion in latent space.**

$$
z_{0}\sim p_{D}(z_{0})\triangleq\mathcal N(\cdot\mid0,I),
\qquad x\mid z_{0}\sim p_{D}(x\mid z_{0};\phi),
$$

$$
\mathrm dz_{t}=f(z_{t},t)\,\mathrm dt+g(t)\,\mathrm dW_{t},\qquad t\in[0,T].
$$

Two structural assumptions: (a) the transition density $p_D(z_t\mid z_0)$ is available analytically, so no Fokker–Planck equation need be solved; and (b) $z_T$ is nearly uninformative, $p_D(z_T)\approx p_D(z_T\mid z_0)$. **The key structural observation is that marginally this latent diffusion is just noise to noise, but conditioned on $x$ it transports the complicated $p_D(z_0\mid x)\propto p_D(z_0)p_D(x\mid z_0)$ into the tractable $p_D(z_T\mid x)=p_D(z_T)$.**

**Step 3: the reverse SDE and the encoding process.** Writing $\tilde t=T-t$, the exact reverse-time SDE is

$$
\mathrm dz_{\tilde t}=-\Bigl(f(z_{\tilde t},\tilde t)-g(\tilde t)^{2}\nabla_{z_{\tilde t}}\log p_{D}(z_{\tilde t}\mid x)\Bigr)\mathrm d\tilde t
+g(\tilde t)\,\mathrm d\widetilde W_{\tilde t} .
$$

The unavailable score in it is replaced by a network $s(z_{\tilde t},x,\tilde t;\theta)$ and paired with the target: $x\sim\pi(x)$, $z_T\sim p_E(z_T)\triangleq p_D(z_T)$. **The encoding process is never simulated**; it only supplies the density $p_E$ used in the loss, and during training all samples of $x$ and $z_t$ come from the decoding process.

**Step 4: a simulation-free unbiased estimator.** The trainable objective $\mathcal L(\theta,\phi)$ delivered by Theorem 1 (below) contains a divergence term, handled by a Hutchinson trace estimator plus Monte Carlo over $t$:

$$
\mathcal L(\theta,\phi)=\mathbb E_{p_D}\bigl[\log p_D(x\mid z_0;\phi)+U(x)\bigr]
+\mathbb E_{t\sim\mathcal U[0,T],\ \epsilon\sim p(\epsilon),\ (x,z_t)\sim p_D}\bigl[\mathcal L_t(x,z_t,\epsilon;\theta)\bigr],
$$

$$
\mathcal L_{t}(x,z_{t},\epsilon;\theta)=\frac{Tg(t)^{2}}{2}
\left(\|s(z_{t},x,t;\theta)\|^{2}
+2\,\frac{\partial\bigl[\epsilon^{\top}s(z_{t},x,t;\theta)\bigr]}{\partial z_{t}}\,\epsilon\right),
$$

with $p(\epsilon)\sim\mathrm{Rademacher}^{d}$ ($\mathbb E[\epsilon]=0$, $\mathrm{Cov}(\epsilon)=I$). **This is where "simulation-free" comes from**: the conditional distribution $p_D(z_t\mid z_0)$ of a standard diffusion SDE can be sampled exactly, so training needs no numerical SDE or ODE integration.

**Step 5: write the boundary conditions into the score network.** The exact score satisfies $\nabla_{z_0}\log p_D(z_0\mid x)=\nabla_{z_0}[\log p_D(x\mid z_0)+\log p_D(z_0)]$ at one end and $\nabla_{z_T}\log p_D(z_T\mid x)=\nabla_{z_T}\log p_D(z_T)$ at the other. The network is therefore written as an interpolation satisfying both ends exactly:

$$
s(z,x,t;\theta)=\Bigl(1-\tfrac tT\Bigr)\nabla_{z_{0}}\bigl[\log p_{D}(x\mid z_{0}=z)+\log p_{D}(z_{0}=z)\bigr]
+\tfrac tT\,\nabla_{z_{T}}\log p_{D}(z_{T}=z)
+\tfrac tT\Bigl(1-\tfrac tT\Bigr)s'(z,x,t;\theta),
$$

with $s'$ a trainable network. This makes the error in $s$ vanish at $t=0$ and $t=T$.

**Step 6: the generalised Hamiltonian dynamics (GHD) decoder.** From $z_0$, generate an initial position-velocity pair $(y,v)$ and apply leapfrog-type updates ($l$ the iteration index, $\odot$ the elementwise product):

$$
v:=v-\frac{\epsilon_{0}e^{\epsilon_{0}\epsilon(l;\phi)}}{2}
\Bigl(\nabla U(y)\odot e^{\frac{\epsilon_{0}}{2}Q_{v}(y,\nabla U(y),l;\phi)}+T_{v}(y,\nabla U(y),l;\phi)\Bigr),
$$

$$
y:=y+\epsilon_{0}e^{\epsilon_{0}\epsilon(l;\phi)}
\Bigl(v_{k}\odot e^{\epsilon_{0}Q_{y}(v_{k},l;\phi)}+T_{y}(v_{k},l;\phi)\Bigr),
$$

and then repeat the first step. The final output and its density are

$$
x=y-\epsilon_{0}e^{\epsilon_{0}\eta(y;\phi)}\nabla U(y)+\sqrt{2\epsilon_{0}e^{\epsilon_{0}\eta(y;\phi)}}\,\xi,
\qquad \xi\sim\mathcal N(0,I),
$$

$$
p_{D}(x\mid z_{0};\phi)=\mathcal N\Bigl(x\Bigm|
y-\epsilon_{0}e^{\epsilon_{0}\eta(y;\phi)}\nabla U(y),\ 2\epsilon_{0}e^{\epsilon_{0}\eta(y;\phi)}I\Bigr).
$$

The paper interprets the last step as a finite-step approximation of the Brownian dynamics $\mathrm dy=-\nabla U(y)\mathrm dt+\mathrm dW_t$. Here $Q_v,T_v,Q_y,T_y$ are trainable networks and the step sizes are parameterised as trainable positive functions. The two claimed benefits are that the decoder makes direct use of $\nabla U$, which helps on multimodal targets, and that only a few iterations are needed.

**Step 7: post-hoc importance reweighting.**

$$
w(x,z_{0})=\frac{\exp(-U(x))\,p_{E}(z_{0}\mid x)}{p_{D}(z_{0})\,p_{D}(x\mid z_{0})}
\ \propto\ \frac{\pi(x)p_{E}(z_{0}\mid x)}{p_{D}(z_{0})p_{D}(x\mid z_{0})},
\qquad
\mathbb E_{\pi}[O(x)]\approx\frac{\sum_{n}w(x^{n},z_{0}^{n})O(x^{n})}{\sum_{n}w(x^{n},z_{0}^{n})} .
$$

The obstacle is that $p_E(z_0\mid x)$ is a path marginal. It is recovered through the **probability flow ODE**:

$$
\mathrm dz_{t}=\Bigl(f(z_{t},t)-\tfrac12g(t)^{2}s(z_{t},x,t;\theta)\Bigr)\mathrm dt,
\qquad z_{T}\sim p_{E}(z_{T}),
$$

$$
\log p_{E}(z_{0}\mid x)=\log p_{E}(z_{T})
+\int_{0}^{T}\nabla\cdot\Bigl(f(z_{t},t)-\tfrac12g(t)^{2}s(z_{t},x,t;\theta)\Bigr)\mathrm dt .
$$

Note that an ODE **is** solved here, but only at evaluation and reweighting time, never during training. This also yields a free-energy lower bound

$$
\log Z\ \ge\ \mathbb E_{p_{D}(x,z_{0})}\bigl[\log w(x,z_{0})\bigr],
$$

tight when $p_D(x,z_0)=p_E(x,z_0)$. The paper explicitly recommends using this lower bound as a **training-quality indicator** — larger is better. The relative effective sample size is

$$
\mathrm{rESS}=\frac{\bigl(\sum_{n=1}^{N}w(x^{n},z_{0}^{n})\bigr)^{2}}{N\sum_{n=1}^{N}w(x^{n},z_{0}^{n})^{2}}\in(0,1].
$$

### Theorems

**Theorem 1 (the only numbered theorem in the paper).** For the decoding and encoding processes above,

$$
D_{\rm KL}\bigl(p_D(x)\|\pi(x)\bigr)
\le D_{\rm KL}\bigl(p_D(x,z_{[\cdot]})\,\|\,p_E(x,z_{[\cdot]})\bigr)
=\mathcal L(\theta,\phi)
+\int_{0}^{T}\frac{g(t)^{2}}{2}\mathbb E_{p_D}\bigl[\|\nabla_{z_t}\log p_D(z_t)\|^{2}\bigr]\mathrm dt
+\log Z,
$$

where $z_{[\cdot]}=\{z_t\}_{t\in[0,T]}$ is the whole latent path and

$$
\mathcal L(\theta,\phi)=\mathbb E_{p_D}\bigl[\log p_D(x\mid z_0;\phi)+U(x)\bigr]
+\int_{0}^{T}\frac{g(t)^{2}}{2}\,\mathbb E_{p_D}
\Bigl[\|s(z_{t},x,t;\theta)\|^{2}+2\,\nabla_{z_{t}}\!\cdot s(z_{t},x,t;\theta)\Bigr]\mathrm dt .
$$

**Equality holds when $z_0\perp z_T$ in the decoding process and $s(z,x,t;\theta)\equiv\nabla_{z_t}\log p_D(z_t\mid x;\phi)$.** The last two terms of the bound are parameter-independent, so $\mathcal L$ alone is the trainable objective. The proof is in Appendix C, with supporting propositions in Appendices A, B (independence of $z_T$ and $x$), D, E and F.

### Numerical experiments

**Example 1: 2D synthetic energies.** MoG2 (two isotropic Gaussians, $\sigma^2=0.5$, centres 10 apart), MoG2(i) (unequal variances, $\sigma_1^2=1.5$, $\sigma_2^2=0.3$, centres 10 apart), MoG6 ($\sigma^2=0.1$), MoG9 ($\sigma^2=0.3$), Ring and Ring5. The metric is MMD against reference samples, with 500,000 samples per method.

| Target  | EDG      | BG   | PIS  | V-HMC | L2HMC |
| ------- | -------- | ---- | ---- | ----- | ----- |
| MoG2    | 0.01     | 1.90 | —    | —     | —     |
| MoG2(i) | **0.50** | 1.63 | —    | 1.56  | 0.94  |
| MoG6    | 0.01     | 2.64 | —    | —     | —     |
| MoG9    | 0.02     | 0.07 | 0.42 | —     | —     |
| Ring    | 0.01     | 0.05 | —    | —     | —     |
| Ring5   | 0.02     | 0.18 | 0.78 | —     | —     |

**The unequal-variance bimodal MoG2(i) is the discriminating case**: every baseline is $\ge0.94$, and only EDG reaches 0.50. BG is poor across the multimodal cases, PIS collapses on MoG9 and Ring5, and V-HMC and L2HMC are competitive except on MoG2(i).

**Example 2: Bayesian logistic regression.** Australian (15 covariates), German (25) and Heart (14); accuracy and AUC averaged over 32 independent runs, with HMC step size 0.01. EDG is best on all six numbers: AU $84.96\pm1.67$ / $92.82\pm0.69$, GE $79.40\pm1.74$ / $82.79\pm1.46$, HE $88.02\pm3.90$ / $95.10\pm1.23$.

**Example 3: Covertype.** Binary classification with 581,012 samples and 54 features, under a hierarchical Bayesian posterior; BG and EDG use the minibatch unbiased log-density

$$
\log\pi(x)\approx\log p(x)+\frac{|\mathcal D_{\rm train}|}{|\mathcal B|}\sum_{(L,D)\in\mathcal B}\log p(L\mid D,x).
$$

EDG reaches accuracy $70.13\pm2.13$ while all four baselines land between $49.88$ and $51.51$ — **the widest gap in the paper**.

**Example 4: Lennard-Jones systems.** LJ13 ($d=39$) and LJ55 ($d=165$), non-periodic boundaries, with reference MCMC data taken from the literature. Qualitatively, EDG's interatomic-distance histograms match the test data while L2HMC and PIS deviate enough to be relegated to a separate appendix figure. Quantitatively (relative effective sample size over 256 samples and 10 random seeds):

| System | EDG             | BG              | PIS             |
| ------ | --------------- | --------------- | --------------- |
| LJ13   | $0.132\pm0.048$ | $0.006\pm0.002$ | $0.005\pm0.001$ |
| LJ55   | $0.098\pm0.014$ | $0.004\pm0.000$ | $0.004\pm0.000$ |

About a twentyfold improvement in effective sample size, **but still far below 1**. The paper says explicitly that it does not handle equivariance, leaving it for future work.

**Example 5: the 2D Ising model.** $16\times16$, so $d=256$, under the continuous relaxation

$$
\pi(x)=\exp\Bigl(-\tfrac12x^{\top}(K(T)+\alpha I)^{-1}x\Bigr)\prod_{i=1}^{N}\cosh(x_{i}),
$$

with $K$ a temperature-dependent symmetric matrix and $\alpha$ chosen so that $K+\alpha I\succ0$; discrete spins are recovered through $\pi(s\mid x)=\prod_i(1+e^{-2s_ix_i})^{-1}$, and

$$
\log Z=\log Z_{\text{Ising}}+\tfrac12\ln\det(K+\alpha I)-\tfrac N2\bigl[\ln(2/\pi)-\alpha\bigr].
$$

Lower-bound estimates of $\log Z_{\rm Ising}$ are reported at $T=2.0,\dots,2.7$ (batch $n=256$, standard deviation of the mean $\mathrm{std}/\sqrt n$): EDG is largest, that is tightest, at $T=2.0$ ($270.32\pm0.18$ against NeuralRG $260.24\pm0.13$ and PIS $210.17\pm0.43$), and also largest at $2.1$, $2.2$, $2.3$, $2.5$ and $2.7$, while **NeuralRG is better at $T=2.4$ and $T=2.6$**. So the conclusion is "over most of the temperature range" rather than all of it.

**Example 6: ablation (MMD over 5000 samples).**

| Variant         | MoG2 | MoG2(i) | MoG6 | MoG9 | Ring | Ring5 |
| --------------- | ---- | ------- | ---- | ---- | ---- | ----- |
| VAE without GHD | 1.86 | 1.62    | 2.57 | 2.10 | 0.12 | 0.23  |
| VAE with GHD    | 0.01 | 2.43    | —    | —    | 1.68 | —     |
| EDG without GHD | 0.06 | 1.01    | 0.04 | 0.05 | 0.02 | 0.04  |
| Full EDG        | 0.01 | 0.50    | 0.01 | 0.02 | 0.01 | 0.02  |

The reading: **the diffusion encoder carries most of the gain, and GHD adds a further consistent improvement — neither component suffices on its own.**

**Cost.** Simulation-free training is a structural property rather than an experimental finding: because $p_D(z_t\mid z_0)$ is analytic, the training algorithm never calls an SDE or ODE solver. The measurable consequence is training memory of **16.75 MB** for EDG against 17.19 for L2HMC, 22.41 for BG and 55.83 for PIS (on a single NVIDIA RTX 2080 Ti / 11 GB, on the MoG2 task, measured by `torch.cuda.memory_allocated`). The paper concedes that BG's **per-epoch** time is shorter, but argues EDG reaches better performance in fewer epochs. Another self-declared cost: EDG is **slightly slower** than PIS at reweighting, because the neural-ODE log-density requires a divergence computation that PIS does not.

### Relation to the other papers

It shares its orientation with the [[en/computational-mathematics/paper-notes/scientific-machine-learning/normalizing-flows-for-densities|density-flow strand]], but with the setting exactly reversed: there one has samples and wants a density, here one has a density and wants samples.

It shares the Ling Guo + Hao Wu + Tao Zhou author group with papers 95 and 98, and all three are built on latent-variable and variational machinery: paper 89 uses a variational-autoencoder-style evidence lower bound with the encoder replaced by a diffusion process, paper 95 uses a conditional variational autoencoder with a set-transformer encoder and a negative evidence lower bound, and paper 98 uses a latent-variable model coupled to a Gaussian process. The idea running through all three is **replacing a hand-designed variational posterior with something structurally richer** — a reverse SDE in paper 89, a permutation-invariant set encoder in paper 95, a Gaussian-process-interpolated latent field in paper 98.

Outside that trio, the methodologically closest are papers 105 and 107: both are probabilistic models applying flows or diffusion to dynamics, and paper 107 likewise chooses a base distribution able to absorb the hardest part of the target. The instantaneous-change-of-variables machinery paper 89 uses to compute log-densities along the probability flow ODE is the same machinery behind continuous normalizing flows. It contrasts sharply with the deterministic half of this page (papers 90, 102, 94, 101, 103): those solve a **given** PDE by least squares or collocation, while paper 89 solves a **sampling** problem whose only handle on the target is pointwise evaluation of $U$ and $\nabla U$. What both halves share is that the object governing the target — a differential operator, or an energy — is known analytically, so neither needs labelled data. The Ising and Lennard-Jones experiments make paper 89 the item on this list closest to computational statistical physics.

> [!note] Coverage status
> The constructions, losses, algorithms and main results of all six papers have been checked equation by equation against the preprint or journal full text. The quantitative results are verifiable to differing degrees: the table values for papers 94, 103, 105 and 89 are transcribed item by item; paper 101 mainly reports curves rather than numerical tables, so this page records only its setups, ablation structure and qualitative conclusions; and paper 81's per-example errors live in the figures and tables of Section 5 and were not transcribed here, leaving only the verified headline order of magnitude (two to three orders of magnitude over a standard multi-scale network). Still unverified are the analysis of high-frequency amplification in paper 101's Appendix A (its main-text Fourier scaling argument is verified) and the journal status of papers 101, 103 and 105, none of which carries journal information on its arXiv record. Rahaman's Theorem 1, Theorems 1 and 2 of the frequency principle, and the MscaleDNN construction in the background section have all been checked against the originals, with their scope restrictions recorded alongside.

## Where the six papers stand

| Paper | How frequency content is handled                                    | Does the structure change                 |
| ----- | ------------------------------------------------------------------- | ----------------------------------------- |
| 81    | posterior capture of dominant modes from the current solution's DFT | rebuilt (two criteria)                    |
| 94    | absorbed analytically into the Gaussian wave-packet ansatz          | not applicable                            |
| 101   | fixed bank + learnable envelope + cross-attention                   | tokens augmented only, backbone untouched |
| 103   | not addressed (splitting structure + $t$ in the input)              | unchanged                                 |
| 105   | not addressed (conditional flows + shared summary)                  | unchanged                                 |
| 89    | not addressed (Boltzmann sampling)                                  | unchanged                                 |

The difference in theoretical content is equally worth setting side by side:

| Paper | Theoretical content                                                                                                                                                                                                                        |
| ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 81    | four approximation theorems: a bound for the standard network, a bound for the down-scaled network (numerator carrying $(\max\{kh,C_1\})^{5q+3}$), a bound for $C^1$ compositions, and a **frequency-free** bound for band-limited targets |
| 94    | no new theorems; both the $O(\sqrt\varepsilon)$ modelling error and the $\mathcal E_t/\varepsilon$ amplification are cited from prior work                                                                                                 |
| 101   | no theorems; the high-frequency amplification analysis in Appendix A is unverified, and the paper concedes the lack of a rigorous analysis                                                                                                 |
| 103   | linear stability analysis (including unconditional instability for $q_{\bm\theta}dt\ge2$), the $dt^3$ Taylor bound of Theorem 3.1, and the $C_1dt^3+C_2dt\,h^2$ error decomposition                                                        |
| 105   | Lemma B.1 (the two optimal summaries need not coincide) and Proposition B.1 (they coincide when a sufficient summary exists); no convergence or error theorem for FLUID itself                                                             |
| 89    | Theorem 1: the augmented path-KL bound and the explicit decomposition of $\mathcal L(\theta,\phi)$, with equality conditions                                                                                                               |

One judgement runs through all of them: **frequency content is either measured (paper 81's DFT capture, paper 101's posterior threshold), or eliminated analytically (paper 94's wave-packet ansatz), or left to the model to weight over a sufficiently wide fixed bank (paper 101's cross-attention). Any approach that fixes frequency content a priori pays the price of frequency mismatch — and paper 81 turned that price into a concrete quantity through $\|F\|_{C^1[-1,1]}=\infty$.**

Paper 101 then adds a qualification to that judgement: **"biased toward low frequencies" is not an intrinsic property of the network but a property of the loss.** Under a regression loss the network favours low frequencies; under a $k^4$-weighted residual loss it favours high ones. So what really has to be avoided is not a bias in one particular direction but assuming a direction without measuring it.

There is a second judgement, unrelated to frequency but equally pervasive: **every adaptation that works in these six papers rests on an indicator computable on the current approximation.** Paper 81 uses the solution's discrete Fourier spectrum to decide which subnetworks to rebuild, paper 101 uses the same spectrum to decide which tokens to inject, paper 103 uses the physics residual $\mathcal R_t$ to decide when to retrain, and paper 89 uses the free-energy lower bound and relative effective sample size to judge training quality. None of these indicators needs a reference solution — which is exactly why they are usable on real problems.

## Sources for this page

- J. Huang, R. You, and T. Zhou, [_Frequency-adaptive multi-scale deep neural networks_](https://doi.org/10.1016/j.cma.2025.117751), Comput. Methods Appl. Mech. Engrg. 437 (2025), 117751 (preprint [arXiv:2410.00053](https://arxiv.org/abs/2410.00053)).
- Y. Wang, L. Guo, H. Wu, and T. Zhou, [_Energy-based diffusion generator for efficient sampling of Boltzmann distributions_](https://doi.org/10.1016/j.neunet.2025.108126), Neural Networks 194 (2026), 108126 (preprint [arXiv:2401.02080](https://arxiv.org/abs/2401.02080)).
- J. Huang, R. You, and T. Zhou, [_Deep learning for the semi-classical limit of the Schrödinger equation_](https://doi.org/10.1016/j.jcp.2026.114869), J. Comput. Phys. 558 (2026), 114869 (preprint [arXiv:2509.04453](https://arxiv.org/abs/2509.04453)).
- X. Feng, T. Tang, X. Wan, and T. Zhou, _Overcoming spectral bias via cross-attention_, [arXiv:2512.18586](https://arxiv.org/abs/2512.18586), submitted to J. Comput. Phys. (no journal information on the arXiv record).
- J. Huang, Y. Qian, and T. Zhou, _PI-DOSnet: a physics-informed deep operator-splitting network for evolution partial differential equations_, [arXiv:2606.22514](https://arxiv.org/abs/2606.22514), submitted to J. Comput. Phys. (no journal information on the arXiv record).
- T. Cui, X. Feng, C. Pei, X. Wan, and T. Zhou, _FLUID: flow-based unified inference for dynamics_, [arXiv:2604.07169](https://arxiv.org/abs/2604.07169) (the arXiv record carries no journal-reference field, so "submitted to Comput. Methods Appl. Mech. Engrg." is not verified by the preprint itself).

Background literature:

- N. Rahaman, A. Baratin, D. Arpit, F. Draxler, M. Lin, F. Hamprecht, Y. Bengio, and A. Courville, _On the spectral bias of neural networks_, ICML 2019 ([arXiv:1806.08734](https://arxiv.org/abs/1806.08734)).
- Z.-Q. J. Xu, Y. Zhang, T. Luo, Y. Xiao, and Z. Ma, _Frequency principle: Fourier analysis sheds light on deep neural networks_, Commun. Comput. Phys. ([arXiv:1901.06523](https://arxiv.org/abs/1901.06523)).
- Z. Liu, W. Cai, and Z.-Q. J. Xu, _Multi-scale deep neural network (MscaleDNN) for solving Poisson-Boltzmann equation in complex domains_, Commun. Comput. Phys. 28 (2020), 1970–2001 ([arXiv:2007.11207](https://arxiv.org/abs/2007.11207)).
