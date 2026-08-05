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

The two phases want almost opposite things. The first reaches sixth order in time on one-dimensional smooth problems, computes conditional expectations by quadrature, and supplies error estimates with constants. The second works at $d=10^4$ to $10^5$ but attains only first order in time and relative errors of $10^{-2}$ to $10^{-3}$; three of its six papers (86, 96 and 108) prove no convergence theorem at all, and the theorems of the other three each carry a restricted scope — paper 93 gives a first-order rate in time, paper 100 gives rigorous a priori bounds, and the theorem of paper 97 holds only in a Fourier surrogate setting. **This is not one thing done well and done badly; it is two different things.**

![How multistep schemes raise the backward accuracy](assets/diagrams/tao-zhou-papers/en/fbsde-multistep.svg)

## Why a probabilistic representation

Consider the coupled Markovian forward-backward system

$$
X_t=X_0+\int_0^t b(s,X_s,Y_s,Z_s)\,\mathrm ds+\int_0^t\sigma(s,X_s,Y_s,Z_s)\,\mathrm dW_s,
$$

$$
Y_t=\xi+\int_t^T f(s,X_s,Y_s,Z_s)\,\mathrm ds-\int_t^T Z_s\,\mathrm dW_s,
$$

with $\xi=\varphi(X_T)$. The system is decoupled when $b$ and $\sigma$ do not depend on $(Y_t,Z_t)$. The nonlinear Feynman-Kac relation ties $(Y_t,Z_t)$ to the solution of a semilinear parabolic equation and its gradient, so this framework supplies a solution route that needs **no spatial grid**: every quantity is a conditional expectation along simulated paths.

The cost moves accordingly, and where it moves to is the key to the whole line of work. At the discrete level the backward equation always takes the shape

$$
Y_i=\mathbb E_i\bigl[Y_{i+1}\bigr]+\Delta t\,\mathbb E_i\bigl[f(\cdots)\bigr]+\cdots,
\qquad
\mathbb E_i[\,\cdot\,]:=\mathbb E[\,\cdot\mid\mathcal F_{t_i}],
$$

so **every difficulty of a method lands in two mutually independent directions**: how few time levels suffice for a given temporal order, and how the conditional expectation is evaluated at each node. The two can be improved separately — paper 25 describes itself as "part two" of the series begun by paper 8, keeping its time discretisation verbatim and replacing only the spatial machinery, which is the clearest evidence of that division of labour.

### The cost lands on the conditional expectation

The algorithm chosen for the conditional expectation decides how many dimensions a paper can actually reach. The table below reorders the whole topic by that algorithm; the last column records **the highest experimental dimension each paper reports**, not a theoretical ceiling.

| How the conditional expectation is computed         | Papers              | Highest tested dimension |
| --------------------------------------------------- | ------------------- | ------------------------ |
| Gauss-Hermite quadrature + Lagrange interpolation   | 8, 19, 23, 26       | 1                        |
| Sparse grids (CGL space + GH quadrature)            | 25                  | 6                        |
| Sinc quadrature (uniform nodes, no interpolation)   | 63                  | —                        |
| Particle filter + single-sample stochastic gradient | 50                  | 3                        |
| Small particle ensemble + local linear regression   | 100                 | $10^4$                   |
| Neural network + adversarial test functions         | 86, 93, 96, 97, 108 | $10^5$                   |

**This table is the shortest possible summary of the topic** — from the first row to the last, the reachable dimension grows by five orders of magnitude, paid for by retreating from "sixth order with a proof" to "first order with numerical evidence".

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

  Its stability window is the single range $1\le k\le6$, fixed by the modulus of the characteristic roots crossing 1 between $k=6$ and $k=7$; numerically the observed order at $k=7$ degrades to between 4 and 5, and $k=8$ gives negative rates. **A second construction of the same name must be kept distinct**: interpolating the integrand over several future levels and then integrating is due to Zhao, Zhang and Ju (_SIAM J. Numer. Anal._ 48(4) 2010), and its windows differ by direction, $K_y\in\{1,\dots,7,9\}$ and $K_z\in\{1,2,3\}$. Paper 8 instead differentiates the reference integral identities in $t$ into reference ordinary differential equations. Both extend the $\theta$-scheme of Zhao, Chen and Peng (_SIAM J. Sci. Comput._ 28(4) 2006).

- **18 (multistep schemes for FBSDEs with jumps)** carries the idea to jump processes, where the backward equation acquires a further unknown $U$ conjugate to the compensated jump measure. The combinatorial obstacle is that the number of jumps inside the span of a $k$-step scheme is unbounded; the remedy is to **count only one jump per time step**, while the forward equation still uses Euler and the backward quantities of interest retain high-order rates.
- **23 (deferred correction methods)** takes a different route to order: instead of adding time levels it repeatedly solves a **residual equation** inside each large step, each sweep raising the accuracy of the interpolant by one order. Its strongest evidence is reaching **$K=12$**, far past the $k\le6$ barrier of the multistep family, at the price of an interpolation operator smooth enough to differentiate twice.
- **33 (explicit theta-schemes for mean-field backward equations)** and **61 (an explicit multistep scheme for mean-field FBSDEs)** handle the mean-field case, where the coefficients depend on the law of the solution itself, so each time level must advance the solution and its law together.
- **35 (explicit deferred correction for second-order FBSDEs)** applies the deferred correction of paper 23 to the second-order class of paper 19, and makes it explicit.
- **68 (strong stability preserving multistep schemes)** imports the strong-stability-preserving concept from ODE numerics: **once paper 47 has made stability half of convergence, "design a convergent scheme" becomes an optimisation problem with a definite objective**, and this paper designs coefficients against that objective.

### Second-order FBSDEs and control

- **16 (probabilistic high-order schemes for fully nonlinear parabolic PDEs)** and **19 (high-order schemes for second-order FBSDEs)** treat second-order (fully nonlinear) systems in the sense of Cheridito, Soner, Touzi and Victoir, which correspond to fully nonlinear parabolic equations. A second-order process $\Gamma$ appears alongside $Y$ and $Z$ (as does the drift $A$ of $Z$ itself), so the number of objects to discretise rises from two to four. The design freedom that matters in this setting is that **the forward diffusion may be chosen freely**, since it is only a probe of the PDE rather than part of the problem. The experiments of paper 19 show all four components attaining order $k$ **simultaneously**, and connect the framework to stochastic optimal control — provided the Hamiltonian's $\inf_\alpha$ can be computed in closed form, which is exactly what papers 86 and 96 later set out to avoid.
- **25 (spectral sparse grid approximations for multi-dimensional FBSDEs)** uses sparse grids for the spatial representation, targeting the exponential growth in evaluation points for conditional expectations. It keeps high order in six dimensions with runtimes growing polynomially rather than exponentially; **six dimensions is the measured ceiling of the deterministic-grid route, not its floor**.
- **26 (an efficient gradient projection method for stochastic optimal control)** and **41 (highly accurate schemes for stochastic optimal control via FBSDEs)** take the Pontryagin route: an adjoint BSDE supplies the gradient of the cost functional, and an outer loop moves the control along the negative gradient and projects back onto the convex admissible set. The adjoint equation of paper 26 is **linear in $(p,q)$**, which is why its nominally implicit left-endpoint rectangle scheme needs no iteration; its error estimate is $O((\Delta t)^2)+O((\Delta x)^4/(\Delta t)^2)$, so **the two mesh parameters cannot be refined independently**. Paper 41 is its direct sequel, replacing the inner solver with the high-order scheme of paper 8.
- **50 (data-driven feedback control)** and **51 (a Gauss-Seidel type method for dynamic nonlinear complementarity problems)** sit at the edge of the family. Paper 50 treats **partial observation**: once observations enter, filtering and control cannot be separated, so it maintains the conditional distribution with a particle filter and replaces every conditional expectation by a single-sample stochastic gradient. Paper 51 contains no BSDE at all; it alternates between a smooth differential subsystem and a nonsmooth complementarity subsystem, and its contribution is the convergence proof that this naive alternation had always lacked.

### Stability theory for discretisations

- **47 (a unified probabilistic discretisation scheme for FBSDEs)** is the theoretical core of this line. It sets up mean-square notions of stability and consistency for FBSDE discretisations so that convergence follows from those two, as in ODE numerical analysis, rather than being proved scheme by scheme. **The empirical window of paper 8 finds its theoretical home here, and the reverse-engineering of paper 68 becomes possible because of it.**
- **63 (Sinc-theta schemes for backward stochastic differential equations)** uses Sinc approximation for the conditional expectations. What it changes is not the accuracy but **the location of the nodes**: Gauss-Hermite nodes are roots of a polynomial and cannot be aligned with the grid, which forces a spatial interpolation step that costs time, caps the spatial accuracy, and destroys analysability at once. Sinc nodes are uniform with a tunable spacing, so all three costs disappear together, and the paper therefore **describes itself as the first attempt to analyse a fully discrete BSDE scheme**, with second order in time and exponential accuracy in space.

### Martingale deep learning

- **86 (SOC-MartNet)** treats the Hamilton-Jacobi-Bellman equation with the specific goal of **not computing $\inf_u H$ explicitly**. It makes three substitutions: the pointwise minimum principle becomes an integral one, the PDE residual becomes a martingale property, and the conditional expectation becomes a weak condition against a family of test functions. A test network probes the property with measurable functions and picks the most violated direction, a value network minimises the detected violation, and a control network updates the feedback control in the same loop. **The loss contains no time recursion**, so training parallelises across time.
- **93 (deep random difference method for high-dimensional quasilinear parabolic equations)** replaces the second-order operator by a **first-order random difference**: a difference quotient between two nearby points, averaged over a random jump size, reproduces a second-order operator, so neither Hessians nor automatic differentiation are needed. It is also the theoretical closure of the line — it identifies its own discretisation with the martingale conditions of papers 86 and 96 as algebraically the same object, so its first-order-in-time convergence theorem underwrites all three.
- **96 (martingale deep learning for very high dimensional quasi-linear equations and stochastic optimal controls)** separates a pilot process from a system process, which localises the martingale structure in time and so **parallelises across time and space at once**, with no derivatives in the loss and no path re-simulation.
- **97 (DeepSPoC)** implements sequential propagation of chaos with deep learning: **the network itself serves as the particles' memory**, so a new particle only interacts with its predecessors and the $N$ particles never have to reside in memory simultaneously.
- **100 (a derivative-free localised stochastic method for very high dimensional semi-linear parabolic PDEs)** replaces the neural network entirely by **local linear regression over a small particle ensemble**, recovering rigorous a priori error bounds and interpretability, with every experiment run on a laptop.
- **108 (deep policy iteration for high-dimensional mean-field games)** uses a **regenerative reformulation** to rewrite the finite-horizon game as an infinite-horizon regenerative process, so a one-step random map replaces full path simulation in the coupled HJB and Fokker-Planck system.

## The numerical evidence at a glance

The two tables below place one representative measurement from each paper side by side; the details and the remaining examples are on the close-reading pages.

**Classical schemes.**

| Paper | Setup                                                      | Measured result                                                                                                                       |
| ----- | ---------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| 8     | $d=1$, $N=16$–$256$, 8 GH nodes, double precision          | orders $k$ observed for $k=1$–$5$ (1.000/1.973/3.002/3.922/5.196); $k=6$ polluted by rounding, $k=8$ diverges                         |
| 19    | $d=1$, 10 GH points per dimension, **quadruple precision** | all four of $Y,Z,\Gamma,A$ attain order $k$; at equal error $k=2$ beats $k=1$ by more than an order of magnitude (10.82 s vs 135.0 s) |
| 23    | $d=1$, $N=4$–$12$ large steps                              | order $K$ up to $K=12$, **no $k\le6$ barrier**                                                                                        |
| 25    | $q=2$–$6$, CGL/GH sparse grids                             | high order retained in six dimensions, runtime polynomial in $q$; the three-step scheme in example 1 reaches only 2.632 in $E_Y$      |
| 26    | $d=1$, $M=10^5$ samples, $\rho_i=1/\sqrt i$                | first order observed in all four examples; feedback control improves $J$ from 0.848 to 0.660                                          |
| 50    | $d=1$, 500 particles, 1000 SGD steps                       | cost $9.5\times10^{-4}$ against $7.6\times10^{-3}$ for the finest-grid full solver, in 0.93 s against 1560 s                          |
| 51    | rectifier circuit and a projected dynamical system         | both theoretical regimes observed                                                                                                     |

**Martingale deep learning.**

| Paper | Highest dimension and hardware         | Measured result                                                                                                                                          |
| ----- | -------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 86    | $d=10^4$ (HJB), 8 A100 GPUs            | first order in time, $\mathcal O(N^{-1.01})$; accuracy comparable to deep BSDE at about half the wall clock; eight GPUs cut 5032 s to 773 s at $d=2000$  |
| 93    | $d=10^5$                               | QLP-1 relative $L^1$ error $1.29\times10^{-2}$ in 3761 s; at $d=10^5$, 13.26 min and 16.8 GB against RS-PINN's 720 min and 45.6 GB, with PINN infeasible |
| 96    | $d=10^4$                               | relative error $5.5\times10^{-3}$ in 295 s; Allen-Cahn at $d=100$ reaches $3.2\times10^{-3}$ in **under 6.8 s**                                          |
| 97    | $d=8$ (porous medium)                  | handles Barenblatt weak solutions, a singular Keller-Segel kernel and Lévy-driven fractional equations; **no order is reported anywhere in the paper**   |
| 100   | $d=10^4$ (Burgers), MacBook Pro M1 Pro | first order in time on all three equation classes; Allen-Cahn at $d=100$ reaches absolute error $1.2\times10^{-5}$                                       |
| 108   | $d=10^4$ (LQ-1), 8 RTX 4090 GPUs       | going from $d=1$ to $d=1000$ inflates the three error measures by only 1.4–4.1 times; the $d=10^4$ run takes 5258 s                                      |

> [!warning] What limits all of these numbers
> Every convergence rate in the classical table comes from **smooth examples in one to at most six dimensions** with conditional expectations evaluated by tensor or sparse quadrature, so none of them says anything about behaviour as the dimension grows. Every error in the deep-learning table is measured **inside the region the pilot process or sampling ensemble actually explored** — remark 1 of paper 96 says so outright, paper 93 encodes it in a weighted norm, and paper 108 displays it directly in its $d=1$ visualisation. **"Accurate" in this family is always a statement with a domain attached.**

## Three judgements that run through the topic

### High-order accuracy can be added in one direction only

The conclusion of paper 8 generalises: when the two directions of a coupled system carry different accuracy requirements, first check whether the error from the cheaper direction actually reaches the quantity of interest. The generator theorem there is exactly such a check — it shows that at the level of the one-step derivative approximation the forward diffusion need only match at the left endpoint, so using Euler forward does not pollute the backward order. The same sentence is used again in paper 19, where the forward diffusion of a second-order system may be chosen freely, and in another form in paper 25, where the time discretisation is kept verbatim and only the spatial machinery is replaced.

### Replace "the residual vanishes" with a testable property

The value of the martingale route is that it changes how the residual is tested. Pointwise residual testing is unavailable in high dimension because there is no grid, whereas testing whether a process is a martingale needs only expectations along simulated paths, so dimension no longer enters the cost directly. That is the same conversion as "represent a distribution by a sampleable model" in the [[en/computational-mathematics/paper-notes/scientific-machine-learning/index|scientific machine learning topic]]: swap a condition that cannot be verified directly for an equivalent one that sampling can verify.

The third of paper 86's substitutions deserves separate notice, because it turns "compute a conditional expectation" into "train a test network". The price is every instability of adversarial training; the return is that conditional expectations are never estimated at all. Each of the five later papers then removes one more part of the machine: paper 96 removes the derivatives and the path re-simulation, paper 93 removes stochastic analysis in favour of Taylor expansions and moment identities and thereby earns the family's first convergence rate, which underwrites papers 86 and 96 along the way, paper 100 removes the neural network itself, paper 97 removes the stored particle trajectories, and paper 108 removes full path simulation.

### Accuracy against dimension is a priced trade

Read the two tables side by side and the point is unavoidable: the classical multistep and $\theta$-schemes reach sixth order with error estimates at $d\lesssim10$, while the martingale methods reach first order and relative errors of $10^{-2}$ to $10^{-3}$ at $d=10^4$ to $10^5$, mostly without a convergence theorem. **The comparison table in paper 93 is unusually honest about this**: its accuracy is generally below PINN and RS-PINN, which the authors attribute to the $O(h)$ truncation carried by the random difference operator, and the entire benefit sits in runtime and memory, growing with dimension.

The practical form of this judgement is: **first ask how many significant digits the quantity of interest needs, then choose an end of the spectrum.** For a one- or low-dimensional problem that needs high accuracy, the mature answer in this topic is the line of papers 8 and 19; for a thousand-dimensional problem that needs two or three digits, the answer is the line from 86 to 108. Paper 100 points to the middle ground: keep the martingale time discretisation, swap the network for an analysable local regression, and pay with a dimension ceiling back at $10^4$ and an order back at one, in exchange for rigorous a priori bounds and laptop-scale cost.

## Sources for this topic

Numbers and records are in the [[en/computational-mathematics/paper-notes/catalog|catalogue]]; per-paper references appear at the end of each close-reading page.
