---
title: How Constraints Enter the Loss - Augmented Lagrangians and Network Bases
description: Papers 60, 90 and 102 - giving boundary conditions, local irregularity and the divergence constraint to structure instead of penalties
lang: en
translation: computational-mathematics/paper-notes/scientific-machine-learning/variational-and-basis-networks
tags:
  - paper-notes
  - scientific-machine-learning
  - network-bases
---

> [!note] Coverage of this page
> Papers **60** (_Commun. Comput. Phys._ 31(3), 2022), **90** (_Commun. Comput. Phys._ 39(2), 2026) and **102** (submitted to _Commun. Comput. Phys._, [arXiv:2603.17906](https://arxiv.org/abs/2603.17906)). The connecting thread is Jianguo Huang and Tao Zhou, and the shared technical stance is to **reduce network training to a (sequence of) least-squares problems** while letting structure rather than penalties carry the constraints.

## 60: replacing the penalty with an augmented Lagrangian

### The problem

Variational problems with essential (Dirichlet-type) boundary conditions are awkward for network solvers because a network function cannot satisfy the constraint exactly, not even at interpolation nodes. The paper names two standard workarounds and their defects. One writes the constraint into the ansatz as $\varphi(x;\theta)=\ell(x)\psi(x;\theta)+\bar g(x)$, which requires a known level-set-like function $\ell$ and a known extension $\bar g$. The other adds a quadratic boundary penalty whose accuracy is sensitive to the penalty weight. This paper switches to a genuine saddle-point formulation so the penalty parameter no longer needs aggressive tuning.

### The augmented Lagrangian and its saddle point

With $V$ a Hilbert space, $B:V\to W$ bounded linear (typically $W=L^2(\Gamma)$) and $V_g=\{v\in V:Bv=g\text{ on }\Gamma\}$, the problem is $\min_{v\in V_g}J(v)$ and the augmented Lagrangian is

$$
\mathcal L_\beta(v,\mu)=J(v)-\langle\mu,\,Bv-g\rangle_W+\frac{\beta}{2}\|Bv-g\|_W^2 .
$$

Note the **minus** sign on the multiplier term; that is the sign convention the paper uses. The minimax problem $\min_{v\in V}\max_{\mu\in W}\mathcal L_\beta(v,\mu)$ is equivalent to the original, because

$$
\max_{\mu\in W}\mathcal L_\beta(v,\mu)=
\begin{cases}
J(v), & v\in V_g,\\
+\infty, & \text{otherwise},
\end{cases}
$$

with dual functional $F_\beta(\mu)=\min_{v\in V}\mathcal L_\beta(v,\mu)$.

For the second-order elliptic problem $-\operatorname{div}(A\nabla u)+cu=f$ with $u|_\Gamma=g$,

$$
J(v)=\tfrac12\int_\Omega\bigl[A\nabla v\cdot\nabla v+cv^2-2fv\bigr]\mathrm dx,
$$

$$
\mathcal L_\beta(v,\mu)=J(v)-\int_\Gamma\mu(x)\bigl(v(x)-g(x)\bigr)\mathrm dx
+\frac{\beta}{2}\int_\Gamma\bigl(v(x)-g(x)\bigr)^2\mathrm dx .
$$

The linear eigenvalue problem $-\nabla\cdot(p\nabla u)+qu=\rho u$ is handled by normalising $\tilde v=v/\|v\|_{0,\Omega}$ and minimising $\tilde J(v)=\int_\Omega[p\nabla\tilde v\cdot\nabla\tilde v+q\tilde v^2]\mathrm dx$; the nonlinear (Gross-Pitaevskii-type) case adds a $\tilde v^4$ term, with the boundary terms unchanged.

### Realising the multiplier update in parameter space

The primal and dual variables get their own ResNet, $u\approx\varphi^u(x;\theta_u)$ and $\lambda\approx\varphi^\lambda(x;\theta_\lambda)$, with $h_0=Vx$, $h_\ell=h_{\ell-1}+\sigma(W_\ell h_{\ell-1}+b_\ell)$ and $\varphi=a^{\mathsf T}h_L$.

Here the paper's central technical difficulty appears: the infinite-dimensional multiplier update

$$
\mu_{k+1}=\mu_k-\beta_k(Bv_k-g)
$$

cannot be applied to network **parameters**, because the right-hand side is a function and parameter space has no matching addition. The fix is a least-squares projection in parameter space:

$$
\theta^\mu_{k+1}=\arg\min_{\theta_\nu}
\bigl\|\varphi^\nu-\varphi^\mu_k-\beta_k\bigl(B\varphi^v_k-g\bigr)\bigr\|^2_{L^2(\Gamma)},
$$

with $\varphi^\nu=\varphi^\lambda(x;\theta_\nu)$, $\varphi^\mu_k=\varphi^\lambda(x;\theta^\mu_k)$ and $\varphi^v_k=\varphi^u(x;\theta^v_k)$. The paper calls this its projection technique. It translates an explicit update in function space into a fitting problem in parameter space, and that step is what carries the classical augmented Lagrangian algorithm over to networks.

The Monte Carlo forms of the two losses are

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

The outer loop fixes $v_k,\mu_k$ and updates the multiplier network with Adam, then fixes $\beta_k,\mu_{k+1}$ and updates the solution network with Adam, then sets $\beta_{k+1}=\alpha\beta_k$ with $\alpha\ge1$. The penalty parameter **increases**, but because the multiplier carries most of the constraint force the initial $\beta_0$ need not be large, and that is the mechanism behind "the choice of the penalty parameter is flexible".

What the paper establishes is the classical equivalence and saddle-point characterisation at the infinite-dimensional level; it proves no convergence or approximation theorem for the discretised network method. The experiments cover second-order elliptic problems, linear eigenvalue problems and nonlinear eigenvalue problems in two and three dimensions, against the penalty method and against solving the minimax problem directly by stochastic gradient descent-ascent, and a remark states explicitly that an efficient descent-ascent method for this minimax problem remains open.

## 90: localise the irregularity, then give the locality a different scale

### How random-basis methods fail

Randomised-basis PDE solvers (extreme learning machines, random feature methods, transferable neural networks) fix hidden-layer weights and biases at random and solve only for the output coefficients by least squares. These are extremely accurate for smooth solutions, and the paper states plainly that accuracy deteriorates markedly when the exact solution has low regularity. The target here is **localised** loss of regularity — sharp peaks and re-entrant-corner singularities, where the feature length scale is orders of magnitude below the domain scale and a globally scaled random basis cannot resolve it. The paper defines "low-regular" as not $H^2$-smooth, or with derivatives that are very large near some points.

A shallow network is written as a basis expansion,

$$
u_{NN}(x)=\sum_{m=1}^{M}\alpha_m\,\sigma\bigl(w_m^{\top}x+b_m\bigr)+\alpha_0,
$$

and because $\{w_m,b_m\}$ are preset and frozen, $u_{NN}$ lies in $V_M=\mathrm{span}\{\psi_0,\dots,\psi_M\}$ with $\psi_0\equiv1$ and $\psi_m(x)=\sigma(\omega_m^{\top}x+b_m)$.

There are two ways to preset the hidden layer. The first draws $w_m$ uniformly on $[-R,R]^d$ and $b_m$ on $[-R,R]$. The second is the transferable-network reparameterisation, writing each neuron as

$$
\sigma\bigl(w_m^{\top}x+b_m\bigr)=\sigma\bigl(\gamma_m(a_m^{\top}x+r_m)\bigr),
\qquad \|a_m\|_2=1,\ \gamma_m\ge0,
$$

so that $w_m=\gamma_ma_m$ and $b_m=\gamma_mr_m$. The **location** parameters are drawn as

$$
a_m=\frac{X_m}{\|X_m\|_2},\qquad r_m=U_m,
$$

with $X_m$ i.i.d. $d$-dimensional standard Gaussians and $U_m$ i.i.d. uniform on $[0,1]$, giving uniformly distributed partition hyperplanes in the unit ball. The **shape** parameter is shared, $\gamma_m\equiv\gamma$, and tuned using realisations of Gaussian random fields as auxiliary functions with no PDE information. This split turns "is the random basis any good?" into two separately controllable questions: how well the hyperplanes cover the domain, and what effective wavelength the activations carry.

### Domain decomposition and an unweighted least-squares loss

Write $\bar\Omega=\cup_{k=0}^{K}\bar\Omega_k$ with disjoint interiors, where each $\Omega_k$ for $k\ge1$ touches only $\Omega_0$ along $\Gamma_k=\partial\Omega_k\cap\partial\Omega_0$. The boundary value problem becomes a subdomain system with $C^1$ transmission,

$$
Lu_k=f\ \text{in}\ \Omega_k,
\qquad u_k=g\ \text{on}\ \partial\Omega_k\cap\partial\Omega,
$$

$$
u_k=u_0\ \text{on}\ \Gamma_k,
\qquad
\frac{\partial u_k}{\partial n_k}=\frac{\partial u_0}{\partial n_k}\ \text{on}\ \Gamma_k .
$$

The collocation least-squares loss sums every residual with **coefficient one**, with no penalty weights:

$$
\min_{\alpha_0,\dots,\alpha_K}\
\sum_{k=0}^{K}\Bigl(\sum_{x\in X_{f_k}}\bigl|L\tilde u_k(x)-f(x)\bigr|^2
+\sum_{x\in X_{g_k}}\bigl|\tilde u_k(x)-g(x)\bigr|^2\Bigr)
+\sum_{k=1}^{K}\sum_{x\in X_{\Gamma_k}}
\Bigl(\bigl|\tilde u_k(x)-\tilde u_0(x)\bigr|^2
+\Bigl|\tfrac{\partial\tilde u_k}{\partial n_k}(x)-\tfrac{\partial\tilde u_0}{\partial n_k}(x)\Bigr|^2\Bigr).
$$

The contrast with a weighted PINN loss is the point: because the basis is frozen, for a linear operator this is a linear least-squares problem $\min_\alpha\|F\alpha-T\|_2^2$ solved directly, and the weights stop being free parameters.

Semilinear operators use Gauss-Newton. Setting $u_k^{n+1}=u_k^{n}+v_k^{n}$, the increment solves the linearised transmission problem

$$
DL(u_k^n;v_k^n)(x)=f(x)-Lu_k^n(x)\ \text{in}\ \Omega_k,
$$

$$
v_k^n=g-u_k^n\ \text{on}\ \partial\Omega_k\cap\partial\Omega,
\qquad
v_k^n-v_0^n=u_0^n-u_k^n\ \text{on}\ \Gamma_k,
\qquad
\frac{\partial v_k^n}{\partial n_k}-\frac{\partial v_0^n}{\partial n_k}
=\frac{\partial u_k^n}{\partial n_k}-\frac{\partial u_0^n}{\partial n_k}\ \text{on}\ \Gamma_k,
$$

where $DL(u;v)$ is the Gâteaux derivative of $L$ at $u$ in direction $v$. Expanding the increment in the same frozen basis makes each Newton step another linear least-squares solve.

### Three components of the adaptivity

**The indicator.** On $\Omega_0$, take the mean residual

$$
\mathcal L_{\Omega_0}(\alpha_0^*)=\frac{1}{|X_{f_0}|}\sum_{x\in X_{f_0}}
\bigl[L\tilde u_0(\alpha_0^*,x)-f(x)\bigr]^2 ,
$$

and place the new peak centre at the pointwise-residual maximiser

$$
x_{K+1}=\arg\max_{x\in X_{f_0}}
\Bigl|L\Bigl(\sum_{m=0}^{M_0}\alpha^*_{m,0}\psi_{m,0}\Bigr)(x)-f(x)\Bigr| ,
$$

then split into $\Omega^{*}=\Omega\setminus B_r(x_{K+1})$ and $\Omega^{**}=B_r(x_{K+1})\cap\Omega$.

**Recentring and rescaling the local basis.** On a new subdomain $\Omega_K=B_{r_K}(x_K)\cap\Omega$,

$$
\psi_{0,K}(x)=1,
\qquad
\psi_{m,K}(x,c_K)=\sigma\bigl(c_K\,w_{m,K}^{\top}(x-x_K)+b_{m,K}\bigr).
$$

This localises the multi-scale network idea: a larger $c_K$ shortens the effective wavelength of every basis function, so the peak becomes resolvable.

**Scale selection.** The coefficient $c_K$ is found by brute-force search over the integers $1,\dots,10$: for each candidate, build the scaled basis and solve a local least-squares problem using only the four local equations on $\Omega_K$ (interior residual, boundary data, and the two interface matching conditions against the **frozen** current $\tilde u_0$), then keep the best. Brute force is affordable here precisely because each trial is a small linear least-squares solve rather than a network training run.

## 102: let an operator identity carry the divergence constraint

### The problem

Network solvers for incompressible flow almost always impose $\operatorname{div}u=0$ as a penalty, and the paper names two consequences directly: choosing the penalty parameter is a critical issue, and whether one uses a velocity-pressure or velocity-vorticity formulation, the resulting system consists of coupled PDEs whose unknowns must be solved simultaneously at high complexity. The two-dimensional stream-function trick makes the constraint exact, but the paper notes that extending it to three dimensions is a major challenge. The target is therefore a formulation that is divergence-free **by construction** in both two and three dimensions and lets velocity and pressure be solved **sequentially**.

### A general principle and two identities

For a linear constraint $\mathcal A(u)=0$, find $\mathcal G$ with $\mathcal A\circ\mathcal G=0$ and use the ansatz $u=\mathcal G(v)$. Here $\mathcal A=\operatorname{div}$ and $\mathcal G=\mathbf{curl}$. In two dimensions,

$$
\mathbf{curl}\,\psi:=\Bigl[\tfrac{\partial\psi}{\partial y},\,-\tfrac{\partial\psi}{\partial x}\Bigr]^{\mathrm T},
\qquad
\operatorname{curl}\,v:=\tfrac{\partial v_2}{\partial x}-\tfrac{\partial v_1}{\partial y}.
$$

Two theorems rewrite the advection term in the potential variable. In two dimensions, for $\phi\in C^3(\Omega)$ with $u=\mathbf{curl}\,\phi$,

$$
\operatorname{curl}\bigl((u\cdot\nabla)u\bigr)=-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi .
$$

In three dimensions, for $\phi\in[C^3]^3$ with $u=\mathbf{curl}\,\phi$ and $\operatorname{div}\phi=0$,

$$
\mathbf{curl}\bigl((u\cdot\nabla)u\bigr)
=(\Delta\phi\cdot\nabla)\mathbf{curl}\,\phi-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi .
$$

The proof route is $\mathbf{curl}((u\cdot\nabla)u)=(u\cdot\nabla)\omega-(\omega\cdot\nabla)u$ together with $\omega=\mathbf{curl}\,\mathbf{curl}\,\phi=\nabla(\nabla\cdot\phi)-\Delta\phi=-\Delta\phi$. The second identity uses $\operatorname{div}\phi=0$, which is why the three-dimensional formulation must impose a divergence-free condition on the potential.

### The decoupled subproblems

Taking $\operatorname{curl}$ of the two-dimensional Stokes equations eliminates the pressure, since $\operatorname{curl}\nabla p=0$, leaving a single-variable fourth-order problem:

$$
\nu\Delta^{2}\phi=\operatorname{curl}\,f\ \text{in}\ \Omega,
\qquad
\phi=\frac{\partial\phi}{\partial n}=0\ \text{on}\ \Gamma .
$$

The clamped boundary condition comes from no-slip: $u\cdot n_0=\nabla\phi\cdot\tau_0=0$ and $u\cdot\tau_0=-\nabla\phi\cdot n_0=0$, so $\phi$ is constant on $\Gamma$ and can be normalised to zero. Two-dimensional Navier-Stokes becomes

$$
\nu\Delta^{2}\phi-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi=\operatorname{curl}\,f .
$$

Three-dimensional Stokes reads

$$
\nu\Delta^{2}\phi=\mathbf{curl}\,f,
\qquad
\operatorname{div}\phi=0\ \text{in}\ \Omega,
\qquad
\phi\cdot n=0,\ \ \mathbf{curl}\,\phi=0\ \text{on}\ \Gamma,
$$

and three-dimensional Navier-Stokes is the genuinely new formulation:

$$
\nu\Delta^{2}\phi+(\Delta\phi\cdot\nabla)\mathbf{curl}\,\phi
-(\mathbf{curl}\,\phi\cdot\nabla)\Delta\phi=\mathbf{curl}\,f,
$$

with the same three constraints. The paper contrasts this with vorticity-vector-potential finite-difference and finite-element work, which must approximate the potential and the vorticity simultaneously, whereas here the single unknown is $\phi$; it calls the result the extension of the pure stream-function formulation to three dimensions.

Pressure recovery is a **first-order gradient system**, not a Poisson equation:

$$
\nabla p=f+\nu\Delta u\ \ (\text{Stokes}),
\qquad
\nabla p=f+\nu\Delta u-(u\cdot\nabla)u\ \ (\text{Navier-Stokes}),
\qquad p(x_0)=0,
$$

with $x_0$ any interior point fixing the additive constant. The paper notes that this recovers a classical solution rather than the weak-solution framework common in earlier work.

### Discretisation and cost

The basis is identical to paper 90, and for vector-valued unknowns **all components share the same hidden basis** with independent output coefficients. The two-dimensional Stokes velocity step collocates three groups of equations at interior and boundary points and solves the resulting over-determined system by least squares; the pressure step solves

$$
\min_{\beta}\Bigl\|
\begin{bmatrix}\partial_x\Psi_{in}\\ \partial_y\Psi_{in}\\ \Psi_0^{\mathrm T}\end{bmatrix}\beta
-\begin{bmatrix}\tilde F_1\\ \tilde F_2\\ 0\end{bmatrix}\Bigr\|_2^2,
$$

where the last row is $p(x_0)=0$. Navier-Stokes uses the Gauss-Newton linearisation

$$
(u\cdot\nabla)u\ \longrightarrow\
u^{(k)}\!\cdot\!\nabla u^{(k+1)}+u^{(k+1)}\!\cdot\!\nabla u^{(k)}-u^{(k)}\!\cdot\!\nabla u^{(k)},
$$

which the paper calls a Newton plus linear least-squares framework and notes is equivalent to discretising the nonlinear system first and then applying Gauss-Newton; an appendix gives a Picard/Oseen-type linearisation used as an initialiser.

The cost accounting quantifies what decoupling buys. Taking $\mathcal O(mn^2)$ for $\min_x\|Ax-b\|_2^2$ with $A\in\mathbb R^{m\times n}$ and $m\gg n$:

| Dimension | Decoupled divergence-free network                       | Coupled baseline              |
| --------- | ------------------------------------------------------- | ----------------------------- |
| Two       | $\mathcal O((I+2J)(M+1)^2)+\mathcal O((2I+1)(M+1)^2)$   | $\mathcal O((3I+2J)(3M+3)^2)$ |
| Three     | $\mathcal O((4I+4J)(3M+3)^2)+\mathcal O((3I+1)(M+1)^2)$ | $\mathcal O((4I+3J)(4M+4)^2)$ |

The saving is the $(M+1)^2$ against $(3M+3)^2$ or $(4M+4)^2$ factor: two small systems solved in sequence instead of one large system with three or four coupled fields.

> [!warning] A symbol in the preprint
> Equation (4.11) prints $\mu$ in the linearised equation while the viscosity is written $\nu$ everywhere else; in context it should be $\nu$.

## The shared judgement

| No. | Constraint taken out of the penalty | Carried instead by                            | What training becomes                    |
| --- | ----------------------------------- | --------------------------------------------- | ---------------------------------------- |
| 60  | essential boundary condition        | a multiplier network and augmented Lagrangian | two alternating stochastic optimisations |
| 90  | interface matching and local scale  | domain decomposition and rescaled local bases | a sequence of linear least squares       |
| 102 | the divergence condition            | a $\mathbf{curl}$ ansatz and two identities   | two sequential least-squares solves      |

All three answer one question: **when a constraint enters the loss as a penalty, its weight becomes a free parameter that must be tuned; giving the constraint to structure makes that parameter disappear.** The prices differ. Paper 60 adds a network and a projection, paper 90 adds a decomposition layer and a scale search, and paper 102 raises a second-order equation to fourth order. These prices are predictable, whereas tuning a penalty weight is not.

## Coverage check

| Item                                              | Paper | Status                                                                          |
| ------------------------------------------------- | ----- | ------------------------------------------------------------------------------- |
| Augmented Lagrangian and minimax equivalence      | 60    | functional form, sign convention, equivalence, dual                             |
| Concrete energies for the three instances         | 60    | elliptic, linear eigenvalue, nonlinear eigenvalue                               |
| Parameter-space projection for the multiplier     | 60    | the obstruction, the projection target, the loop and increasing $\beta$         |
| Frozen basis and the two presetting routes        | 90    | basis definition, uniform draws, reparameterisation, location/shape split       |
| Domain decomposition and unweighted least squares | 90    | transmission conditions, the loss, linear and semilinear cases                  |
| Residual indicator, local rescaling, scale search | 90    | all three components and why brute force is affordable                          |
| Divergence-free ansatz and the two identities     | 102   | general principle, two- and three-dimensional identities, proof route           |
| Four decoupled subproblems and pressure recovery  | 102   | fourth-order velocity equations, origin of boundary conditions, gradient system |
| Discretisation, Gauss-Newton and cost accounting  | 102   | shared basis, least-squares form, linearisation, complexity table               |

## Sources for this page

- J. Huang, H. Wang, and T. Zhou, [_An augmented Lagrangian deep learning method for variational problems with essential boundary conditions_](https://doi.org/10.4208/cicp.OA-2021-0176), Commun. Comput. Phys. 31(3) (2022), pp. 966-986 (preprint [arXiv:2106.14348](https://arxiv.org/abs/2106.14348)).
- J. Huang, H. Wu, and T. Zhou, [_Adaptive neural network basis methods for partial differential equations with low-regular solutions_](https://doi.org/10.4208/cicp.OA-2024-0310), Commun. Comput. Phys. 39(2) (2026), pp. 553-577 (preprint [arXiv:2411.01998](https://arxiv.org/abs/2411.01998)).
- J. Cheng, J. Huang, H. Wang, and T. Zhou, _Decoupled divergence-free neural networks basis method for incompressible fluid problems_, [arXiv:2603.17906](https://arxiv.org/abs/2603.17906), submitted to Commun. Comput. Phys.
