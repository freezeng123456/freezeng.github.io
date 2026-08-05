---
title: How Constraints Enter the Loss - Augmented Lagrangians and Neural Network Bases
description: Papers 60, 90 and 102 - handing boundary conditions, local low regularity and the divergence constraint to structure rather than to penalties
lang: en
translation: computational-mathematics/paper-notes/scientific-machine-learning/variational-and-basis-networks
tags:
  - paper-notes
  - scientific-machine-learning
  - neural-network-basis
---

> [!note] Coverage of this page
> Papers **60** (_Commun. Comput. Phys._ 31(3), 2022), **90** (_Commun. Comput. Phys._ 39(2), 2026) and **102** (submitted to _Commun. Comput. Phys._, [arXiv:2603.17906](https://arxiv.org/abs/2603.17906)). The common authorship thread is Jianguo Huang and Tao Zhou, and the common technical orientation is **reducing network training to a (sequence of) least-squares problems** while letting structure rather than penalties carry the constraints.

## 60: replace the penalty with an augmented Lagrangian

### The idea

Variational problems with essential boundary conditions are unfriendly to network solvers for a very concrete reason: a network function cannot satisfy a Dirichlet condition exactly, not even at interpolation nodes. The paper names two prevailing workarounds and the defect of each.

The first builds the constraint into the ansatz, $\varphi(x;\theta)=\ell(x)\psi(x;\theta)+\bar g(x)$. On paper this is perfect — the constraint holds identically — but it requires a known level-set-like function $\ell$ vanishing on the boundary and positive inside, plus a known extension $\bar g$ of the boundary data. Once the domain is complicated, producing those two objects is a problem in itself.

The second adds a quadratic penalty on the boundary. That always works, but it trades accuracy for a tuning problem: too small a penalty and the constraint is not enforced, too large and the loss becomes ill-conditioned so the optimiser cannot move. **The penalty parameter has no natural scale**, so it has to be found by trial.

Paper 60's judgement is that classical optimisation solved this dilemma long ago with the augmented Lagrangian method. Introducing a multiplier carries the first-order information of the constraint explicitly, leaving the penalty responsible only for local convexity. The penalty parameter then need not tend to infinity, and it need not start large — **the multiplier does most of the enforcing and the penalty is auxiliary**. That is the mechanism behind the abstract's claim that "the choice of the penalty parameter is flexible and robust".

Transplanting this classical method to networks meets exactly one real obstacle, and it is the paper's technical core: the classical multiplier update is an **explicit addition in function space**, while a network's multiplier is determined by parameters and parameter space has no corresponding addition. The paper's answer is to rewrite that step as a least-squares fit.

### Setting

Let $V$ be a Hilbert space, $B:V\to W$ bounded and linear (typically $W=L^2(\Gamma)$), and $V_g=\{v\in V:Bv=g \text{ on }\Gamma\}$. The original problem is

$$
\min_{v\in V_g}J(v).
$$

The augmented Lagrangian is

$$
\mathcal L_\beta(v,\mu)=J(v)-\langle\mu,\,Bv-g\rangle_W+\frac{\beta}{2}\|Bv-g\|_W^2 .
$$

Note the **minus** sign on the multiplier term; that is the sign convention the paper actually uses.

The three instances carry these energies.

**Second-order elliptic problem** $-\operatorname{div}(A\nabla u)+cu=f$ with $u|_\Gamma=g$, taking $V=H^1(\Omega)$ and $W=L^2(\Gamma)$:

$$
J(v)=\tfrac12\int_\Omega\bigl[A\nabla v\cdot\nabla v+cv^2-2fv\bigr]\mathrm dx,
$$

$$
\mathcal L_\beta(v,\mu)=J(v)-\int_\Gamma\mu(x)\bigl(v(x)-g(x)\bigr)\mathrm dx
+\frac{\beta}{2}\int_\Gamma\bigl(v(x)-g(x)\bigr)^2\mathrm dx .
$$

The paper cites a reference for the existence of a unique saddle point here.

**Linear eigenvalue problem** $-\nabla\cdot(p\nabla u)+qu=\rho u$ with $u|_\Gamma=0$: normalising $\tilde v=v/\|v\|_{0,\Omega}$ turns the Rayleigh quotient into a minimisation of

$$
\tilde J(v)=\int_\Omega\bigl[p\nabla\tilde v\cdot\nabla\tilde v+q\tilde v^2\bigr]\mathrm dx,
\qquad
\mathcal L_\beta(v,\mu)=\tilde J(v)-\int_\Gamma\mu(x)\tilde v(x)\mathrm dx
+\frac{\beta}{2}\int_\Gamma\tilde v^2(x)\mathrm dx .
$$

**Nonlinear (Gross-Pitaevskii-type) eigenvalue problem** $-\nabla\cdot(A\nabla u)+Vu+u^3=\rho u$ with $u|_{\partial\Omega}=0$ and $\|u\|_{0,\Omega}=1$ adds a quartic term:

$$
\tilde J(v)=\tfrac12\int_\Omega\bigl[A\nabla\tilde v\cdot\nabla\tilde v+V\tilde v^2+\tilde v^4\bigr]\mathrm dx,
$$

with the boundary terms unchanged. The paper cites a reference stating that the non-negative ground-state solution is unique for $1\le d\le3$.

### Derivation

**From the constrained problem to a minimax problem.** Maximise over $\mu$ for fixed $v$. If $v\in V_g$ then $Bv-g=0$, both $\mu$-dependent terms vanish, and $\mathcal L_\beta(v,\mu)=J(v)$ independently of $\mu$. If $v\notin V_g$ then $-\langle\mu,Bv-g\rangle_W$ is linear and non-trivial in $\mu$, and taking $\mu=-t(Bv-g)$ with $t\to+\infty$ makes it grow without bound. Hence

$$
\max_{\mu\in W}\mathcal L_\beta(v,\mu)=
\begin{cases}
J(v), & v\in V_g,\\
+\infty, & \text{otherwise},
\end{cases}
$$

so $\min_{v\in V}\max_{\mu\in W}\mathcal L_\beta(v,\mu)$ is equivalent to the original constrained problem. The dual functional is $F_\beta(\mu)=\min_{v\in V}\mathcal L_\beta(v,\mu)$ and the dual problem is $\max_{\mu\in W}F_\beta(\mu)$.

**Network parameterisation.** The primal and dual variables each get their own ResNet, $u\approx\varphi^u(x;\theta_u)$ and $\lambda\approx\varphi^\lambda(x;\theta_\lambda)$, with

$$
h_0=Vx,
\qquad
h_\ell=h_{\ell-1}+\sigma(W_\ell h_{\ell-1}+b_\ell),
\qquad
\varphi(x;\theta)=a^{\mathsf T}h_L .
$$

**Translating the multiplier update.** The classical update is

$$
\mu_{k+1}=\mu_k-\beta_k(Bv_k-g),
$$

whose right-hand side is a function, while $\mu$ is now determined by $\theta_\lambda$ and parameter space offers no operation corresponding to that addition. The paper rewrites the step as a least-squares projection in parameter space:

$$
\theta^\mu_{k+1}=\arg\min_{\theta_\nu}J_\lambda(\varphi^\nu;\varphi^\mu_k,\varphi^v_k),
\qquad
J_\lambda=\bigl\|\varphi^\nu-\varphi^\mu_k-\beta_k\bigl(B\varphi^v_k-g\bigr)\bigr\|^2_{L^2(\Gamma)},
$$

with $\varphi^\nu=\varphi^\lambda(x;\theta_\nu)$, $\varphi^\mu_k=\varphi^\lambda(x;\theta^\mu_k)$ and $\varphi^v_k=\varphi^u(x;\theta^v_k)$. The paper calls this its **projection technique**. It translates an explicit update in function space into a fitting problem in parameter space: the target $\varphi^\mu_k-\beta_k(B\varphi^v_k-g)$ is a fixed function that can be evaluated pointwise, so the problem becomes making the multiplier network approximate it. This is the step that carries the classical augmented Lagrangian algorithm onto networks.

**Monte Carlo form of both losses.** For the elliptic instance,

$$
\mathcal L_\beta(v,\mu)=|\Omega|\,\mathbb E_\xi\Bigl[\tfrac12 A(\xi)|\nabla v(\xi)|^2
+\tfrac12 c(\xi)v^2(\xi)-f(\xi)v(\xi)\Bigr]
-|\Gamma|\,\mathbb E_\eta\Bigl[\mu(\eta)\bigl(v(\eta)-g(\eta)\bigr)
-\tfrac{\beta}{2}\bigl(v(\eta)-g(\eta)\bigr)^2\Bigr],
$$

$$
J_\lambda(\nu;\mu_k,v_k)=|\Gamma|\,\mathbb E_\eta
\Bigl[\bigl(\nu(\eta)-\mu_k(\eta)-\beta(v_k(\eta)-g(\eta))\bigr)^2\Bigr],
$$

with $\xi\sim\mathrm{Unif}(\Omega)$ and $\eta\sim\mathrm{Unif}(\Gamma)$.

> [!note] A bracket grouping
> The grouping of the $-|\Gamma|\mathbb E_\eta[\cdot]$ bracket in the accessible rendering of the Monte Carlo form is ambiguous. **The sign structure $-\mu(v-g)+\tfrac{\beta}{2}(v-g)^2$ is fixed by the continuous form above**, which should be treated as the authority.

### Theorems

What the paper supplies is the classical infinite-dimensional result: the equivalence of the constrained variational problem with the minimax problem, and the corresponding saddle-point characterisation. **This is classical augmented Lagrangian theory quoted for the infinite-dimensional setting, not a new theorem about neural networks.** The paper proves no convergence or approximation theorem for the discretised (network) method, and saying so plainly marks the boundary of its claims: the method rests on a principle, but discrete-level guarantees are not part of it.

### Numerical experiments

**Algorithm 2 (ALDL) is structured as follows:**

1. Input $\beta_0>0$, an increase parameter $\alpha\ge1$, an outer count `Epoch`, and inner counts `Epoch_u` and `Epoch_λ`;
2. Initialise $\theta^v_0$, $\theta^\mu_0$ with PyTorch's default random initialisation;
3. For $k=0,1,\dots$: fix $\varphi^v_k,\varphi^\mu_k$ and minimise $J_\lambda$ with Adam for `Epoch_λ` steps starting from $\theta^\mu_k$, giving $\theta^\mu_{k+1}$; then fix $\beta_k,\varphi^\mu_{k+1}$ and minimise $\mathcal L_\beta$ with Adam for `Epoch_u` steps starting from $\theta^v_k$, giving $\theta^v_{k+1}$; then set $\beta_{k+1}=\alpha\beta_k$;
4. Output $u=\varphi^v_{\texttt{Epoch}+1}$ and $\lambda=\varphi^\mu_{\texttt{Epoch}+1}$.

It discretises an infinite-dimensional algorithm (the paper's Algorithm 1): update $\mu_{k+1}=\mu_k-\beta_k(Bv_k-g)$; find an approximate minimiser $v_{k+1}$ of $\mathcal L_\beta(v,\mu_{k+1})$ with $\|\partial_1\mathcal L_\beta(v_{k+1},\mu_{k+1})\|\le\tau_k$; then set $\beta_{k+1}=\alpha\beta_k$ and pick a new tolerance $\tau_{k+1}$.

**The penalty parameter does increase**, but since the multiplier carries most of the enforcement, $\beta_0$ need not be large — which is how "flexible choice of the penalty parameter" is actually implemented rather than merely asserted.

The examples of Section 4 and the baselines:

| Example                          | Dimension | Content                                    |
| -------------------------------- | --------- | ------------------------------------------ |
| second-order elliptic (Poisson type) | 2D, 3D | with essential boundary conditions          |
| linear eigenvalue problem        | 2D, 3D    | smallest eigenvalue and its eigenfunction   |
| nonlinear eigenvalue problem     | 2D, 3D    | Gross-Pitaevskii type                       |

There are two baselines: the **penalty-method deep-learning version**, and **solving the minimax problem directly by stochastic gradient descent-ascent (SGDA)**. The second matters, because it tests whether "augmented Lagrangian plus alternating updates" really beats taking the minimax formulation at face value. Remark 3.1 states plainly that designing an efficient descent-ascent method for this minimax problem **remains open**, which both explains the alternating scheme and concedes that its theoretical footing is incomplete.

The qualitative conclusion is that ALDL is markedly less sensitive to $\beta$ than the penalty method and reaches better accuracy at comparable cost.

> [!warning] How far the numerical results were verified
> The composition of the examples, both baselines and the qualitative conclusion are verified; **the per-example error tables could not be confirmed in an accessible rendering**. The specific improvement factors sometimes quoted (eigenfunctions better by 2 to 20 times, eigenvalues by up to 100 times, 25 to 30 per cent faster) are **all unverified and should not be cited**.

### Relation to the others

This is the constrained-optimisation entry in the topic, a counterpoint to the sampling-centred papers. It shares the deep Ritz variational-loss setting with [[en/computational-mathematics/paper-notes/scientific-machine-learning/adaptive-sampling-for-pinns|paper 80]]: paper 60 fixes **how the boundary constraint enters the objective**, paper 80 fixes **where the quadrature points go**. It shares authorship and the "turn training into least squares" orientation with papers 90 and 102, but by the opposite means — paper 60 introduces a multiplier to enforce the constraint, while 90 and 102 abolish penalty weights altogether.

## 90: localise the low regularity, then give the locality a different scale

### The idea

Randomised-basis PDE solvers — extreme learning machines, the random feature method, transferable neural networks — draw the hidden weights and biases at random and then **freeze** them, solving only for the output-layer coefficients by least squares. That single move converts PDE solving from a non-convex optimisation into a linear least-squares problem, which is why these methods are extremely accurate on smooth solutions: their error lives at discretisation level rather than at training level.

The paper states the failure mode plainly: accuracy deteriorates remarkably when the exact solution has low regularity. And the failure has a specific shape — **low regularity is usually local**. Sharp peaks, re-entrant-corner singularities: the solution is smooth over almost the entire domain and has enormous derivatives in a tiny neighbourhood. A globally uniform random basis then loses both ways. To resolve the peak, the effective wavelength of the activations must be pushed down to the peak width; but the same basis then behaves like high-frequency noise over the smooth region and the conditioning collapses.

The paper's response is to refuse the compromise on a single global scale and instead **let different regions carry bases at different scales**. Three parts turn that sentence into an algorithm: use the residual to locate where refinement is needed; carve out a subdomain there and give it a recentred, rescaled basis; and stitch the subdomains together by non-overlapping domain decomposition with $C^1$ matching at the interfaces.

The duality with adaptive sampling is worth naming: **where the residual is large, adaptive sampling adds points; here the method adds a rescaled local basis.** One indicator, two responses.

The paper defines "low-regular" as not $H^2$-smooth, or having very large derivatives around some points.

### Setting

A shallow network written as a basis expansion:

$$
u_{NN}(x)=\sum_{m=1}^{M}\alpha_m\,\sigma\bigl(w_m^{\top}x+b_m\bigr)+\alpha_0,
$$

and because $\{w_m,b_m\}$ are pre-set and frozen, $u_{NN}$ lies in $V_M=\mathrm{span}\{\psi_0,\dots,\psi_M\}$ with

$$
\psi_0(x)=1,
\qquad
\psi_m(x)=\sigma\bigl(\omega_m^{\top}x+b_m\bigr),
\quad 1\le m\le M .
$$

There are two presetting strategies. The first draws $w_m$ uniformly on $[-R,R]^d$ and $b_m$ on $[-R,R]$ over $\Omega=[-1,1]^d$, with $R$ a user parameter. The second is the TransNet re-parameterisation, writing each neuron as

$$
\sigma\bigl(w_m^{\top}x+b_m\bigr)=\sigma\bigl(\gamma_m(a_m^{\top}x+r_m)\bigr),
\qquad \|a_m\|_2=1,\ \gamma_m\ge0,
$$

with the correspondence $w_m=\gamma_ma_m$, $b_m=\gamma_mr_m$ and conversely $a_m=w_m/\|w_m\|_2$, $r_m=b_m/\|w_m\|_2$, $\gamma_m=\|w_m\|_2$. The **location** parameters are drawn as

$$
a_m=\frac{X_m}{\|X_m\|_2},\qquad r_m=U_m,\quad m=1,\dots,M,
$$

with $X_m$ i.i.d. $d$-dimensional standard Gaussians and $U_m$ i.i.d. uniform on $[0,1]$, which distributes the partition hyperplanes uniformly in the unit ball $B_1(0)$. The **shape** parameter is shared, $\gamma_m\equiv\gamma$, and tuned using realisations of Gaussian random fields as auxiliary functions, **with no PDE information at all**.

This split is the precondition for everything that follows: it separates "is the random basis any good" into two independently controllable parts — how well the hyperplanes cover the domain (through $a_m,r_m$) and the effective wavelength of the activations (through $\gamma_m$). Only with that separation does "give the locality a different scale" have a precise meaning, namely change $\gamma$ alone.

### Derivation

**Non-overlapping decomposition and an unweighted least-squares loss.** Take $\bar\Omega=\cup_{k=0}^{K}\bar\Omega_k$ with disjoint interiors, each $\Omega_k$ ($k\ge1$) touching only $\Omega_0$ along $\Gamma_k=\partial\Omega_k\cap\partial\Omega_0$. The boundary-value problem $Lu=f$, $u|_{\partial\Omega}=g$ becomes a $(K{+}1)$-subdomain system with $C^1$ transmission:

$$
Lu_k=f\ \text{in}\ \Omega_k,
\qquad u_k=g\ \text{on}\ \partial\Omega_k\cap\partial\Omega,
$$

$$
u_k=u_0\ \text{on}\ \Gamma_k,
\qquad
\frac{\partial u_k}{\partial n_k}=\frac{\partial u_0}{\partial n_k}\ \text{on}\ \Gamma_k,
\quad 1\le k\le K .
$$

Expanding $\tilde u_k(\alpha_k,x)=\sum_{m=0}^{M_k}\alpha_{m,k}\psi_{m,k}(x)$ and collocating on the $3K{+}2$ point sets gives

$$
\min_{\alpha_0,\dots,\alpha_K}\
\sum_{k=0}^{K}\Bigl(\sum_{x\in X_{f_k}}\bigl|L\tilde u_k(x)-f(x)\bigr|^2
+\sum_{x\in X_{g_k}}\bigl|\tilde u_k(x)-g(x)\bigr|^2\Bigr)
+\sum_{k=1}^{K}\sum_{x\in X_{\Gamma_k}}
\Bigl(\bigl|\tilde u_k(x)-\tilde u_0(x)\bigr|^2
+\Bigl|\tfrac{\partial\tilde u_k}{\partial n_k}(x)-\tfrac{\partial\tilde u_0}{\partial n_k}(x)\Bigr|^2\Bigr).
$$

**Every residual term enters with coefficient $1$; there are no penalty weights.** That contrasts sharply with weighted PINN losses, and the contrast is not stylistic: because the basis is frozen, for linear $L$ this is a linear least-squares problem $\min_\alpha\|F\alpha-T\|_2^2$ that can be solved directly, its solution is unchanged by scaling all residuals by a common constant, and only the relative weights between blocks are a genuine degree of freedom — which the paper sets to unity. In other words, **the penalty weights can be dropped because the solver changed.**

**Gauss-Newton for semilinear operators.** Set $u_k^{n+1}=u_k^{n}+v_k^{n}$; the increment solves the linearised transmission problem

$$
DL(u_k^n;v_k^n)(x)=f(x)-Lu_k^n(x)\ \text{in}\ \Omega_k,
$$

$$
v_k^n=g-u_k^n\ \text{on}\ \partial\Omega_k\cap\partial\Omega,
\qquad
v_k^n-v_0^n=u_0^n-u_k^n\ \text{on}\ \Gamma_k,
$$

$$
\frac{\partial v_k^n}{\partial n_k}-\frac{\partial v_0^n}{\partial n_k}
=\frac{\partial u_k^n}{\partial n_k}-\frac{\partial u_0^n}{\partial n_k}\ \text{on}\ \Gamma_k,
$$

where $DL(u;v)$ is the Gâteaux derivative of $L$ at $u$ in direction $v$. Expanding the increment in the same frozen basis makes each Newton step another linear least-squares solve $\min_{a^n}\|F^na^n-T^n\|_2^2$. The stopping rule is the relative change of the residual norm,

$$
R_{\text{emse}}=\frac{\bigl|\|F^na^{n,*}-T^n\|_2^2-\|F^{n-1}a^{n-1,*}-T^{n-1}\|_2^2\bigr|}
{\|F^{n-1}a^{n-1,*}-T^{n-1}\|_2^2}<\text{tol},
$$

with $\text{tol}=10^{-5}$ in all experiments.

**The three adaptive parts.** First, the **indicator**: the mean residual on $\Omega_0$,

$$
\mathcal L_{\Omega_0}(\alpha_0^*)=\frac{1}{|X_{f_0}|}\sum_{x\in X_{f_0}}
\bigl[L\tilde u_0(\alpha_0^*,x)-f(x)\bigr]^2 ,
$$

with the new peak centre taken at the argmax of the pointwise residual,

$$
x_{K+1}=\arg\max_{x\in X_{f_0}}
\Bigl|L\Bigl(\sum_{m=0}^{M_0}\alpha^*_{m,0}\psi_{m,0}\Bigr)(x)-f(x)\Bigr| ,
$$

after which the domain splits into $\Omega^{*}=\Omega\setminus B_r(x_{K+1})$ and $\Omega^{**}=B_r(x_{K+1})\cap\Omega$.

Second, **recentring and rescaling the local basis**. On the new subdomain $\Omega_K=B_{r_K}(x_K)\cap\Omega$,

$$
\psi_{0,K}(x)=1,
\qquad
\psi_{m,K}(x,c_K)=\sigma\bigl(c_K\,w_{m,K}^{\top}(x-x_K)+b_{m,K}\bigr),
\quad 1\le m\le M_K .
$$

The two modifications do different jobs: subtracting $x_K$ moves the basis's attention onto the peak, and multiplying by $c_K$ compresses the effective wavelength of every basis function. This is the MscaleDNN idea localised to a subdomain.

Third, **scale selection**. The value $c_K$ is chosen by brute-force search over the integers $1,\dots,10$: for each candidate $s$, build the scaled basis, solve a local least-squares problem $\arg\min_{\alpha_K}\mathcal L_K(\alpha_K,s)$ involving only the four local equations on $\Omega_K$ (interior residual, boundary data, and the two interface matching conditions against the **frozen** current $\tilde u_0$), record $\mathrm{Loss}_s$, and take $s^*=\arg\min_{1\le s\le10}\mathrm{Loss}_s$. Brute force is affordable here precisely because each trial is a small linear least-squares solve rather than a network training run — a concrete dividend of the frozen-basis route. The paper notes a bisection search could be used instead.

The final approximation is a piecewise expansion, $u(x)\approx\tilde u_k(\alpha_k^*,x)$ for $x\in\Omega_k$. The representation is therefore **discontinuous across $\Gamma_k$ by construction**, with $C^1$ matching imposed only weakly at interface collocation points.

### Theorems

**The paper contains no theorem, no convergence proof and no error estimate.** It is a purely algorithmic and numerical article; its conclusion calls it "our first attempt to attack such class of problems". Saying this plainly matters, because it means the paper's persuasiveness rests entirely on its experiments — so the details in the next section are the whole of its evidence.

### Numerical experiments

All experiments use the second presetting strategy (the TransNet re-parameterisation). The error metric is relative $L^2$:

$$
\mathrm{err}_{L^2}=\sqrt{\sum_i|\hat u(x_i)-u(x_i)|^2\Big/\sum_i|u(x_i)|^2}.
$$

**Example one: 2-D Poisson with $P$ Gaussian peaks (Section 4.2).** On $\Omega=[-1,1]^2$ with exact solution

$$
u=\sum_{p=1}^{P}\exp\bigl(-1000\bigl[(x-x_p)^2+(y-y_p)^2\bigr]\bigr),
$$

three cases: one peak at $(0.5,0.5)$; two peaks diagonally at $(\pm0.5,\pm0.5)$; four peaks at $(\pm0.5,\pm0.5)$. The $1000$ in the exponent is the point: it squeezes the peak width to $O(0.03)$, roughly two orders below the domain scale, which is exactly what a globally scaled basis cannot resolve.

| Setting               | Value                              |
| --------------------- | ---------------------------------- |
| tolerance $\epsilon$  | $10^{-4}$                          |
| initial basis $M_0$   | 200                                |
| refinement radius $r$ | 0.15                               |
| interior points $X_{f_0}$ | $50\times50$ grid, $J_{f_0}=2500$ |
| boundary points $X_{g_0}$ | 400                            |
| error evaluation grid | $256\times256$                     |

Outcome: the three cases converge to 2, 3 and 5 subdomains, meaning **the number of detected subdomains is exactly the number of peaks plus the one smooth master domain**; the detected centres are $(\pm0.5102,\pm0.5102)$, matching the true peak locations with no prior information supplied. The error decreases monotonically over the per-peak basis counts $M^*\in\{700,800,900,1000\}$ and reaches $\mathrm{err}_{L^2}\approx10^{-4}$ at $M^*=1000$, which the paper describes as significantly smaller than the results of two adaptive-sampling PINN papers. The selected integer scale is $c_k=5$ for every peak subdomain.

**Example two: 2-D semilinear peak problem (Section 4.3).** The equation becomes $-\Delta u+u^2=f$, with the same three peak cases and the same settings. This example exercises the Gauss-Newton branch. Outcome: the same partition, the same selected scales $c_k=5$, and errors **smaller** than in the linear case.

**Example three: 2-D corner singularity (Section 4.4).** L-shaped $\Omega=[-1,1]^2\setminus[0,1]^2$ with exact solution $u=(x^2+y^2)^{1/3}$. This is a different species of low regularity: the singularity comes from the geometry rather than from the source term.

| Setting               | Value          |
| --------------------- | -------------- |
| tolerance $\epsilon$  | $10^{-3}$      |
| initial basis $M_0$   | 600            |
| refinement radius $r$ | 0.32           |
| interior points       | $J_{f_0}=1875$ |
| boundary points       | $J_{g_0}=400$  |

Outcome: a single refinement suffices, $\Omega_1=B_r(x_1)$ with $x_1=(-0.0345,-0.0345)$, very close to the re-entrant corner. The scale is $c_1=5$ at $M^*=800,900$ and $c_1=4$ at $M^*=1000$, with $\mathrm{err}_{L^2}\approx10^{-3}$ at $M^*=1000$. The error is an order of magnitude worse than in the peak examples and the scale is not monotone in $M^*$; both hint that **an algebraic singularity is harder than a Gaussian peak**, because a peak has one clear characteristic wavelength to match while an $r^{2/3}$ singularity is non-smooth at every scale.

**Example four: 3-D Poisson with one peak (Section 4.5).** On $\Omega=[-1,1]^3$ with $u=\exp(-1000[(x-0.5)^2+(y-0.5)^2+(z-0.5)^2])$.

| Setting               | Value                          |
| --------------------- | ------------------------------ |
| tolerance $\epsilon$  | $10^{-4}$                      |
| initial basis $M_0$   | 2000                           |
| refinement radius $r$ | 0.11                           |
| interior points       | $J_{f_0}=10000$                |
| boundary points       | $J_{g_0}=2400$ (400 per face)  |
| subdomain interior points | 8500                       |
| interface points      | 600                            |

Outcome: the detected centre is $(0.5170,0.5050,0.5000)$; the scale grows with the basis count, $c_1=6,7,8$ at $M^*=3000,4000,5000$; and $\mathrm{err}_{L^2}\approx10^{-3}$ at $M^*=5000$.

Read together, the four examples mark the method's boundary: **2-D peaks are the favourable case ($10^{-4}$), while both the corner singularity and the 3-D peak stop at $10^{-3}$**. The 3-D case also needs $M_0=2000$ initial basis functions and $5000$ per subdomain, an order of magnitude above the 2-D runs — the dimensional cost of a frozen-basis method is not removed, only deferred.

The paper's stated open problems are three: find all peak subdomains in a single iteration (currently one peak per round); find a better way to determine the scaling coefficients (currently an integer brute-force search); and handle time-dependent low regularity.

> [!note] Publication details
> The journal volume, issue and DOI are verified; **the page range 553-577 has not been independently re-checked**.

### Relation to the others

Paper 90 belongs to the randomised-basis strand together with paper 102: the same frozen random basis, the same linear-least-squares plus Gauss-Newton machinery, differing in what they repair — paper 90 adapts the **scale** of the basis to handle low regularity, paper 102 adapts the **formulation** to make a constraint exact.

The scaled basis $\sigma(c_Kw^\top(x-x_K)+b)$ is the randomised-basis counterpart of the MscaleDNN construction and a close relative of the multiscale random Fourier feature bank of [[en/computational-mathematics/paper-notes/scientific-machine-learning/spectral-bias-and-generative-solvers|paper 101]]. The three are different answers to the same spectral-bias problem: paper 94 uses a multiscale network, paper 90 searches discretely for a scale on a subdomain, and paper 101 learns amplitudes over a fixed multiscale bank via cross-attention.

The residual-argmax refinement is the "hard" analogue of the deep adaptive sampling the same group uses elsewhere: instead of **adding samples** where the residual is large, it **adds a rescaled local basis** there.

Its relation to paper 60 is two answers to one question: both target constraints that plain penalty losses handle badly, but paper 90 uses no penalty weights at all and relies on the least-squares solver, while paper 60 introduces multipliers for the boundary condition.

> [!note] Citation relation
> Whether paper 90 cites paper 60 could not be confirmed from the verification material behind this page.

## 102: let an operator identity carry the divergence constraint

### The idea

Network solvers for incompressible flow almost always impose $\operatorname{div}u=0$ as a penalty. The paper names two consequences, and they compound. First, choosing the penalty parameter is a critical issue — the same difficulty paper 60 faces. Second, whether one uses a velocity-pressure or a velocity-vorticity formulation, the result is a **coupled** PDE system in which several unknown fields must be solved simultaneously, so the number of columns in the least-squares matrix grows in proportion to the number of fields, and least-squares cost is quadratic in the column count.

Two dimensions have a classical escape: the stream function. Taking $u=\mathbf{curl}\,\psi$ makes incompressibility hold identically because $\operatorname{div}\mathbf{curl}\equiv0$, replacing two components by one scalar. But the paper points out that extending it to three dimensions is a major challenge — in 3-D the potential is a vector and needs a gauge condition of its own.

The paper's move can be stated as a general principle: **for a linear constraint $\mathcal A(u)=0$, find an operator $\mathcal G$ with $\mathcal A\circ\mathcal G=0$ and use the ansatz $u=\mathcal G(v)$.** Here $\mathcal A=\operatorname{div}$ and $\mathcal G=\mathbf{curl}$. The constraint stops being something to enforce and becomes an identity.

The real technical work is the second step. Once the curl ansatz is adopted, the equations must be rewritten in the potential variable, and it is not obvious what the advection term $(u\cdot\nabla)u$ looks like there. The paper's two theorems compute it, and the 3-D one is the key to its new formulation. As a by-product, taking the curl of the whole momentum equation annihilates the pressure gradient (since $\operatorname{curl}\nabla p=0$), so **the velocity subproblem contains no pressure at all**: velocity can be solved first and the pressure recovered afterwards from a separate equation. The coupling is thereby dissolved.

### Setting

The two curl operators in 2-D are

$$
\mathbf{curl}\,\psi:=\Bigl[\tfrac{\partial\psi}{\partial y},\,-\tfrac{\partial\psi}{\partial x}\Bigr]^{\mathrm T},
\qquad
\operatorname{curl}\,\bm v:=\tfrac{\partial v_2}{\partial x}-\tfrac{\partial v_1}{\partial y}.
$$

Two representation lemmas license the ansatz. **Lemma 2.1 (2-D):** a divergence-free $\bm v$ can be written $\bm v=\mathbf{curl}\,\psi$ with $\psi\in H^1(\Omega)$. **Lemma 2.2 (3-D, simply connected):** $\bm v=\mathbf{curl}\,\bm\psi$ with $\operatorname{div}\bm\psi=0$, and one may take $\bm\psi\in H(\mathbf{curl};\Omega)$. The extra 3-D condition $\operatorname{div}\bm\psi=0$ is the gauge fixing, and it has a second job in the derivation below.

### Derivation

**Two advection identities.**

**Theorem 2.1 (2-D).** For $\phi\in C^3(\Omega)$ with $\bm u=\mathbf{curl}\,\phi$,

$$
\operatorname{curl}\bigl((\bm u\cdot\nabla)\bm u\bigr)=-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi .
$$

**Theorem 2.2 (3-D).** For $\bm\phi\in[C^3]^3$ with $\bm u=\mathbf{curl}\,\bm\phi$ and $\operatorname{div}\bm\phi=0$,

$$
\mathbf{curl}\bigl((\bm u\cdot\nabla)\bm u\bigr)
=(\Delta\bm\phi\cdot\nabla)\mathbf{curl}\,\bm\phi-(\mathbf{curl}\,\bm\phi\cdot\nabla)\Delta\bm\phi .
$$

The proof runs in two steps. First the vector identity

$$
\mathbf{curl}\bigl((\bm u\cdot\nabla)\bm u\bigr)=(\bm u\cdot\nabla)\bm\omega-(\bm\omega\cdot\nabla)\bm u,
$$

then the expression of vorticity through the potential,

$$
\bm\omega=\mathbf{curl}\,\mathbf{curl}\,\bm\phi=\nabla(\nabla\cdot\bm\phi)-\Delta\bm\phi=-\Delta\bm\phi .
$$

The second equality consumes $\operatorname{div}\bm\phi=0$. This shows that the gauge condition is not an optional technical appendix in 3-D: **without it the vorticity does not reduce cleanly to $-\Delta\bm\phi$ and the whole identity fails.** That is why the divergence-free condition must be imposed on the potential in three dimensions.

**Four decoupled subproblems.** Take the curl of the 2-D Stokes equation $-\nu\Delta\bm u+\nabla p=\bm f$: the pressure term vanishes because $\operatorname{curl}\nabla p=0$, and the viscous term becomes biharmonic because $-\nu\Delta\operatorname{curl}\mathbf{curl}\,\phi=\nu\Delta^2\phi$, leaving a single-variable fourth-order problem

$$
\nu\Delta^{2}\phi=\operatorname{curl}\,\bm f\ \text{in}\ \Omega,
\qquad
\phi=\frac{\partial\phi}{\partial n}=0\ \text{on}\ \Gamma .
$$

The clamped boundary conditions are derived, not assumed, from no-slip: $\bm u\cdot\bm n_0=\nabla\phi\cdot\bm\tau_0=0$ says the tangential derivative of $\phi$ along the boundary vanishes, so $\phi$ is constant on $\Gamma$ and may be normalised to zero; and $\bm u\cdot\bm\tau_0=-\nabla\phi\cdot\bm n_0=0$ gives the vanishing normal derivative. For a slip condition $\bm u=\bm g$ the two relations become $\nabla\phi\cdot\bm\tau_0=\bm g\cdot\bm n_0$ and $-\nabla\phi\cdot\bm n_0=\bm g\cdot\bm\tau_0$.

The 2-D Navier-Stokes velocity problem adds Theorem 2.1:

$$
\nu\Delta^{2}\phi-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi=\operatorname{curl}\,\bm f,
\qquad
\phi=\frac{\partial\phi}{\partial n}=0\ \text{on}\ \Gamma .
$$

3-D Stokes reads

$$
\nu\Delta^{2}\bm\phi=\mathbf{curl}\,\bm f,
\qquad
\operatorname{div}\bm\phi=0\ \text{in}\ \Omega,
\qquad
\bm\phi\cdot\bm n=0,\ \ \mathbf{curl}\,\bm\phi=\bm 0\ \text{on}\ \Gamma
$$

(with the last condition replaced by $\mathbf{curl}\,\bm\phi=\bm g$ for slip), and 3-D Navier-Stokes is the genuinely new formulation:

$$
\nu\Delta^{2}\bm\phi+(\Delta\bm\phi\cdot\nabla)\mathbf{curl}\,\bm\phi
-(\mathbf{curl}\,\bm\phi\cdot\nabla)\Delta\bm\phi=\mathbf{curl}\,\bm f,
$$

under the same three constraints. Remark 3.3 positions this against vorticity-vector-potential finite-difference and finite-element work: those require simultaneous approximation of potential **and** vorticity, whereas this has the single unknown $\bm\phi$, so the paper calls it the extension of the "pure" stream-function formulation to three dimensions.

**Pressure recovery.** The point is that it is a **first-order gradient system**, not a Poisson equation:

$$
\text{Stokes: }\ \nabla p=\bm f+\nu\Delta\bm u,
\qquad
\text{Navier-Stokes: }\ \nabla p=\bm f+\nu\Delta\bm u-(\bm u\cdot\nabla)\bm u,
\qquad p(\bm x_0)=0,
$$

with $\bm x_0$ any interior point fixing the additive constant. The paper notes this allows recovery of a **classical** solution rather than the weak-solution frameworks commonly adopted in earlier work — taking a divergence to form a Poisson equation would discard first-order information and introduce extra boundary-condition questions, whereas collocating the first-order system directly does not.

### Theorems

Theorems 2.1 and 2.2 are the paper's only theorems, and they are **algebraic identities, not convergence results**. There is no error analysis and no convergence theorem for the method. In this respect it matches paper 90: the randomised-basis strand rests on numerical evidence.

### Numerical experiments

**Basis and discretisation.** The basis is identical to paper 90's (the frozen three-step TransNet re-parameterisation), and for vector-valued fields **all components share the same hidden basis** $\{\psi_m\}$ with independent output coefficients $\bm\alpha^{(i)}\in\mathbb R^{M+1}$ — which is why $\nu\Delta^2\bm\Psi_{in}$ repeats along the diagonal of the block matrices.

For the 2-D Stokes velocity step, on $\bm X^{in}=\{x^{in}_i\}_{i=1}^{I}\subset\Omega$ and $\bm X^{bd}=\{x^{bd}_j\}_{j=1}^{J}\subset\Gamma$ the three collocated equations are

$$
\nu\Delta^{2}\!\sum_{m=0}^{M}\!\alpha_m\psi_m(x^{in}_i)=\operatorname{curl}\,\bm f(x^{in}_i),
\qquad
\sum_m\alpha_m\psi_m(x^{bd}_j)=0,
\qquad
\frac{\partial}{\partial n}\sum_m\alpha_m\psi_m(x^{bd}_j)=0,
$$

and the resulting over-determined system is solved as a linear least-squares problem using the linear solvers provided by NumPy. The pressure step solves

$$
\min_{\bm\beta}\Bigl\|
\begin{bmatrix}\partial_x\bm\Psi_{in}\\ \partial_y\bm\Psi_{in}\\ \bm\Psi_0^{\mathrm T}\end{bmatrix}\bm\beta
-\begin{bmatrix}\tilde{\bm F}_1\\ \tilde{\bm F}_2\\ 0\end{bmatrix}\Bigr\|_2^2,
$$

whose last row is exactly $p(\bm x_0)=0$.

Navier-Stokes uses the Gauss-Newton linearisation

$$
(\bm u\cdot\nabla)\bm u\ \longrightarrow\
\bm u^{(k)}\!\cdot\!\nabla\bm u^{(k+1)}+\bm u^{(k+1)}\!\cdot\!\nabla\bm u^{(k)}
-\bm u^{(k)}\!\cdot\!\nabla\bm u^{(k)},
$$

giving in 2-D the per-iteration linear PDE

$$
\nu\Delta^{2}\phi^{(k+1)}-(\mathbf{curl}\,\phi^{(k)}\!\cdot\!\nabla)\Delta\phi^{(k+1)}
-(\mathbf{curl}\,\phi^{(k+1)}\!\cdot\!\nabla)\Delta\phi^{(k)}
=\operatorname{curl}\,\bm f-(\mathbf{curl}\,\phi^{(k)}\!\cdot\!\nabla)\Delta\phi^{(k)} .
$$

The paper calls this the Newton-LLSQ framework and notes it is equivalent to first deriving a discretised nonlinear system and then applying Gauss-Newton. The appendix's Picard/Oseen-type linearisation ("Scheme I", $\nu\Delta^2\phi^{(k+1)}-\mathbf{curl}\,\phi^{(k)}\cdot\nabla\Delta\phi^{(k+1)}=\mathbf{curl}\,\bm f$) is used as an initialiser, because Gauss-Newton is only locally convergent.

> [!warning] A symbol in the preprint
> Equation (4.11) prints $\mu$ in the linearised equation while the surrounding text consistently uses the viscosity $\nu$. Read it as $\nu$.

**Complexity accounting.** Taking the cost of $\min_x\|Ax-b\|_2^2$ with $A\in\mathbb R^{m\times n}$, $m\gg n$, to be $\mathcal O(mn^2)$:

| Dimension | Decoupled divergence-free network                       | Coupled baseline              |
| --------- | ------------------------------------------------------- | ----------------------------- |
| 2-D       | $\mathcal O((I+2J)(M+1)^2)+\mathcal O((2I+1)(M+1)^2)$   | $\mathcal O((3I+2J)(3M+3)^2)$ |
| 3-D       | $\mathcal O((4I+4J)(3M+3)^2)+\mathcal O((3I+1)(M+1)^2)$ | $\mathcal O((4I+3J)(4M+4)^2)$ |

The saving comes from the $(M+1)^2$ versus $(3M+3)^2$ or $(4M+4)^2$ factor: two small systems solved in sequence rather than one system with three or four coupled fields. Because least-squares cost is quadratic in the column count, replacing $3M$ by $M$ is a factor of nine.

**Error metrics and baselines.** The relative $L^2$ error and the absolute divergence error are

$$
\mathrm{error}\_g=\frac{\sqrt{\tfrac1N\sum_i[g-g_{\rm NN}]^2}}{\sqrt{\tfrac1N\sum_i g^2}},
\qquad
\mathrm{error}\_{\rm div}=\sqrt{\tfrac1N\sum_i\bigl[\operatorname{div}\bm g_{\rm NN}(x_i)\bigr]^2}.
$$

There are two baselines: the standard **coupled** TransNet, and a PINN (ResNet of width $M=30$, depth $L=4$, $\tanh$, Adam, 512 interior plus 256 boundary points, 10000 iterations). Hardware: an Intel Core i9-12900H laptop, 16 GB RAM, RTX 3060; Python 3.12 with PyTorch 2.1.

**Example one: 2-D Stokes, Kovasznay-type solution (Section 5.1).** On $\Omega=(0,2)\times(-0.5,1.5)$,

$$
u=1-e^{\zeta x}\cos(2\pi y),
\quad
v=\tfrac{\zeta}{2\pi}e^{\zeta x}\sin(2\pi y),
\quad
p=\tfrac12 e^{2\zeta x},
\quad
\zeta=\tfrac{1}{2\nu}-\sqrt{\tfrac{1}{4\nu^2}+4\pi^2}.
$$

Training uses a uniform grid of $50\times50=2500$ interior points plus 200 boundary points; testing uses $111\times111=12321$ points.

| Comparison                        | Result                                                              |
| --------------------------------- | -------------------------------------------------------------------- |
| versus coupled TransNet at $\nu=10^{-4}$ | better on both velocity and pressure; divergence at least $\mathcal O(10^{-12})$ |
| sweeping $\nu$ at $M=1000$        | TransNet slightly better at large $\nu$; this method better by nearly two orders for $\nu\le10^{-2}$ |
| runtime                           | about 2 s versus slightly over 5 s for TransNet                      |
| versus PINN                       | velocity and pressure errors at least four orders smaller; PINN divergence error only $\mathcal O(10^{-3})$, and PINN loses all accuracy at low viscosity |

The divergence at low viscosity is the most informative part of this example: the smaller $\nu$, the steeper the exponential layer in the Kovasznay solution and the worse the conditioning of the coupled formulation, whereas the decoupled formulation turns the velocity problem into a pure biharmonic one, unentangled from the pressure.

**Example two: 2-D Navier-Stokes (Section 5.2).** Exact solution

$$
u=16y(y-1)(2y-1)\sin^2(\pi x),
\quad
v=-8\pi y^2(y-1)^2\sin(2\pi x),
\quad
p=\sin(\pi x)\cos(\pi y),
$$

with the same grids; the Gauss-Newton loop is initialised from a PINN trained for 10000 iterations, with at most 40 Gauss-Newton iterations. At $M=1000$ the accuracy for $\nu\le10^{-2}$ is **comparable** to TransNet — as expected, since both use the same initial guess and the same Newton-type solver — but at more than 60 per cent lower cost and with strictly better divergence preservation. This example is therefore not an accuracy claim but a **cost claim plus a constraint claim**.

**Example three: 3-D Stokes (Section 5.3).** On $(0,1)^3$,

$$
u=e^{\cos(\pi y)}\sin(\pi z),
\quad
v=e^{\cos(\pi z)}\sin(\pi x),
\quad
w=e^{\cos(\pi x)}\sin(\pi y),
$$

$$
p=e^{\cos(\pi x)+\sin(\pi y)}+e^{\cos(\pi z)+\sin(\pi x)} .
$$

Halton sampling gives 10000 interior points plus 2400 boundary points ($20\times20\times6$), with 2000 Halton test points.

| Comparison       | Result                                                       |
| ---------------- | ------------------------------------------------------------- |
| velocity         | more accurate than TransNet at every $M$                      |
| pressure         | TransNet about one order better at large $M$, but nearly twice as slow |
| divergence       | $10^{-14}$, six orders better than TransNet                   |
| sweeping $\nu$   | this method better for $\nu\le10^{-4}$ on velocity and $\nu\le10^{-3}$ on pressure |
| cost             | about half of TransNet                                        |

The pressure row deserves attention because it is a weakness the paper reports itself: **solving sequentially means the pressure step inherits the error of the velocity step**, whereas a coupled solve lets the two correct one another. The paper does not sidestep this.

> [!note] Internal inconsistencies in the preprint
> The caption of Fig. 5.7 says $M=1500$ while the surrounding text says $M=2500$; the caption of Fig. 5.9 mistakenly repeats "3D Stokes".

**Example four: 3-D Navier-Stokes (Section 5.4).** On $(0,1)^3$ with the trigonometric-polynomial solution

$$
u=2(y-1)(z-1)\sin y\sin z-2(y-1)\sin y\cos z-2(z-1)\cos y\sin z
$$

($v,w$ by cyclic permutation) and $p=xyz+x^3y^3z-\tfrac{5}{32}$. The initial guess comes from Scheme I after 10 iterations, with at most 15 Gauss-Newton iterations. Outcome: higher accuracy in both velocity and pressure for $\nu<10^{-3}$, plus the usual divergence and cost advantages.

**Overall judgement.** The paper's own summary: divergence error near machine precision ($10^{-14}$); roughly 50 per cent less execution time than comparable frameworks such as TransNet when using many basis functions, and more than 60 per cent for the 2-D Navier-Stokes case; superior velocity accuracy in the low-viscosity, high-Reynolds-number regime. It also concedes one limitation honestly: **Stokes performance is generally better than Navier-Stokes, which the authors attribute to the quality of the Gauss-Newton initial guess.** That is consistent with examples two and four both requiring external initialisation (a PINN or Scheme I): decoupling solved the coupling and the constraint, but not the global convergence of the nonlinear solve.

The code is "available from the corresponding author on reasonable request"; there is no public repository.

### Relation to the others

Paper 102 is the direct sibling of paper 90: the same school, the same frozen random basis (verbatim the same three-step re-parameterisation), the same linear-least-squares plus Gauss-Newton solve. Paper 90 adapts the scale of the basis to handle low regularity, paper 102 adapts the formulation to make a constraint exact; together they are the randomised-basis strand of this topic.

The "eliminate the constraint structurally via $u=\mathcal G(v)$" philosophy is the opposite of paper 60's multiplier route for essential boundary conditions. Which is available depends on whether the constraint has a ready-made potential representation — the divergence constraint has $\mathbf{curl}$, while Dirichlet conditions on a complicated domain have no equally clean counterpart.

The move of "splitting a coupled system into a sequence of cheaper subproblems" is conceptually the same as the temporal operator splitting of [[en/computational-mathematics/paper-notes/scientific-machine-learning/spectral-bias-and-generative-solvers|paper 103]], differing only in the direction of the split: coupling between fields here, coupling between operators there.

## What the three share

| No. | Constraint taken out of the penalty | Carried instead by                       | What training becomes         |
| --- | ----------------------------------- | ---------------------------------------- | ----------------------------- |
| 60  | essential boundary conditions       | a multiplier network plus augmented Lagrangian | two alternating stochastic optimisations |
| 90  | interface matching and local scale  | domain decomposition plus rescaled local bases | a sequence of linear least-squares problems |
| 102 | the divergence-free condition       | a $\mathbf{curl}$ ansatz and its identities | two least-squares solves in sequence |

All three answer the same question: **when a constraint enters the loss as a penalty, its weight becomes a free parameter that must be tuned; hand the constraint to structure and that free parameter disappears.** The prices differ — paper 60 pays an extra network and a projection, paper 90 pays a layer of domain decomposition and a scale search, paper 102 raises a second-order equation to fourth order. Those prices are predictable; tuning a penalty weight is not.

A second judgement concerns the solver: papers 90 and 102 can abolish penalty weights entirely because they replace training with least squares. Least squares is still sensitive to the relative weighting between residual blocks, but it does not need a weight that makes gradient descent converge, and those two difficulties are not of the same order. Paper 60 retains stochastic optimisation and therefore still needs $\beta$; the multiplier only lowers how large $\beta$ has to be.

A third judgement concerns evidence: **none of the three proves a convergence theorem for its discrete method.** Paper 60 quotes classical theory at the infinite-dimensional level; papers 90 and 102 state outright that no analysis is offered. All three therefore stand on their experimental design — and papers 90 and 102 report theirs in enough detail to be re-checked.

## Coverage check

| Item                                             | Paper | Status                                                    |
| ------------------------------------------------ | ----- | ----------------------------------------------------------- |
| Defects of the two workarounds and the multiplier's role | 60 | intuition, problems with the ansatz and the penalty        |
| Augmented Lagrangian and minimax equivalence      | 60   | functional form, sign convention, derivation of the inner maximisation, dual functional |
| Energies of the three instances                   | 60   | elliptic, linear eigenvalue, nonlinear eigenvalue with cited uniqueness |
| Parameter-space projection of the multiplier update | 60  | source of the difficulty, projection target, Monte Carlo forms and the bracket note |
| Numerical experiments of 60                       | 60   | Algorithm 2 and infinite-dimensional Algorithm 1, three example classes, two baselines, Remark 3.1; **error tables and quoted factors unverified** |
| Frozen basis and the two presetting strategies    | 90   | basis definition, uniform draws, re-parameterisation and the location/shape split |
| Domain decomposition and unweighted least squares | 90   | transmission conditions, loss, why penalty weights vanish, Gauss-Newton and its stopping rule |
| Residual indicator, local rescaling, scale search | 90   | all three parts, why brute force is affordable, the constructive discontinuity |
| Numerical experiments of 90                       | 90   | complete settings and results for all four examples, including detected centres, scales, error magnitudes and three open problems |
| Divergence-free ansatz and the two advection identities | 102 | general principle, two lemmas, two theorems with the proof route, the role of the gauge |
| Four decoupled subproblems and pressure recovery  | 102  | fourth-order velocity equation, derivation of the boundary conditions, first-order gradient system and classical solutions |
| Discretisation, Gauss-Newton and cost accounting  | 102  | shared basis, two least-squares steps, linearisation and initialisation, complexity comparison |
| Numerical experiments of 102                      | 102  | solutions, grids, baselines, timings and divergence magnitudes for all four examples, plus the two conceded weaknesses |

## Sources for this page

- J. Huang, H. Wang, and T. Zhou, [_An augmented Lagrangian deep learning method for variational problems with essential boundary conditions_](https://doi.org/10.4208/cicp.OA-2021-0176), Commun. Comput. Phys. 31(3) (2022), pp. 966-986 (preprint [arXiv:2106.14348](https://arxiv.org/abs/2106.14348)).
- J. Huang, H. Wu, and T. Zhou, [_Adaptive neural network basis methods for partial differential equations with low-regular solutions_](https://doi.org/10.4208/cicp.OA-2024-0310), Commun. Comput. Phys. 39(2) (2026), pp. 553-577 (preprint [arXiv:2411.01998](https://arxiv.org/abs/2411.01998)).
- J. Cheng, J. Huang, H. Wang, and T. Zhou, _Decoupled divergence-free neural networks basis method for incompressible fluid problems_, [arXiv:2603.17906](https://arxiv.org/abs/2603.17906), submitted to Commun. Comput. Phys.
