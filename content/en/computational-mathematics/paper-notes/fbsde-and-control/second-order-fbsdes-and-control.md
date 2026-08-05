---
title: Second-Order FBSDEs and Stochastic Control
description: Papers 16, 19, 25, 26, 41, 50 and 51 - fully nonlinear equations and control iterations
lang: en
translation: computational-mathematics/paper-notes/fbsde-and-control/second-order-fbsdes-and-control
tags:
  - paper-notes
  - stochastic-optimal-control
  - fully-nonlinear-equations
---

> [!note] Coverage of this page
> Papers **16** (_Commun. Comput. Phys._ 18, 2015), **19** (_Commun. Comput. Phys._ 21, 2017), **25** (_Discrete Contin. Dyn. Syst. Ser. B_ 22(9), 2017), **26** (_SIAM J. Numer. Anal._ 55(6), 2017), **41** (_Numer. Math. Theor. Meth. Appl._ 13, 2020), **50** (_J. Sci. Comput._ 85(2), 2020) and **51** (_SIAM J. Control Optim._ 58(6), 2020).

## Why a "second-order" FBSDE is needed

A first-order forward-backward system corresponds through the nonlinear Feynman-Kac relation to a **semilinear** parabolic equation: $Y_t$ gives the solution and $Z_t$ gives the gradient. Covering **fully nonlinear** equations, where the Hessian enters nonlinearly, requires one more process in the probabilistic representation.

The second-order forward-backward system of Cheridito, Soner, Touzi and Victoir is designed for exactly that: alongside $Y_t$ and $Z_t$ there is $\Gamma_t$, corresponding to the Hessian of the solution. The number of objects to discretise therefore rises from two to three, and each one needs its own conditional-expectation approximation.

Papers **16 (probabilistic high-order schemes for fully nonlinear parabolic PDEs)** and **19 (high-order schemes for second-order FBSDEs)** treat this setting, with the title of paper 19 naming the application: stochastic optimal control. Both belong to the same technical tradition as paper 8 on the [[en/computational-mathematics/paper-notes/fbsde-and-control/multistep-schemes-for-fbsdes|multistep page]], with one extra process to handle. Paper **35** treats the same setting by deferred correction.

Written out, the coupled second-order system that paper 19 treats is three equations:

$$
\begin{cases}
X_t=x+\displaystyle\int_0^tb(s,\Theta_s)\,\mathrm ds+\int_0^t\sigma(s,\Theta_s)\,\mathrm dW_s,\\[4pt]
Y_t=g(X_T)+\displaystyle\int_t^Tf(s,\Theta_s)\,\mathrm ds-\int_t^TZ_s\,\mathrm dW_s,\\[4pt]
Z_t=Z_0+\displaystyle\int_0^tA_s\,\mathrm ds+\int_0^t\Gamma_s\,\mathrm dW_s,
\end{cases}
\qquad
\Theta_t=(X_t,Y_t,Z_t,A_t,\Gamma_t)\in\mathbb R^m\times\mathbb R\times\mathbb R^d\times\mathbb S^d,
$$

with $\mathbb S^d$ the $d\times d$ real symmetric matrices. **The third equation is the whole of what is new in the second-order setting**: it expands $Z$ itself as an Itô process, and its diffusion coefficient $\Gamma$ is the process corresponding to the Hessian. "Decoupled" means $b$ and $\sigma$ do not depend on $(Y,Z,A,\Gamma)$.

The paper's reading of the state of the art is that existing high-order FBSDE schemes work only in low dimension while existing high-dimensional schemes are low order — it cites a 12-dimensional coupled FBSDE example converging at order 1, against paper 25's sparse-grid schemes reaching dimension 6 with rates up to 3. Its stated gap: to the authors' knowledge there was no study of high-order numerical methods for second-order FBSDEs.

## 25: the number of evaluation points for a conditional expectation

In several dimensions a conditional expectation is a multidimensional Gaussian integral, and tensor-product quadrature makes the number of evaluation points grow exponentially with dimension. Paper 25 addresses that growth with **spectral sparse grids**, under the published title _Efficient spectral sparse grid approximations for solving multi-dimensional forward backward SDEs_.

Its position is worth naming: in this thread the temporal accuracy problem is solved by multistep schemes and deferred correction, while the cost of the **spatial** representation is solved by sparse grids and later by Sinc approximation ([[en/computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes|paper 63]]). The two directions of improvement are independent.

## 26: organising the control iteration as a gradient projection

### The problem and the gap in existing routes

The goal is the constrained stochastic optimal control problem

$$
\min_{u\in K}J(u)=\mathbb E\Bigl[\int_0^T\bigl(h(x_t^u)+j(u(t))\bigr)\mathrm dt+k(x_T^u)\Bigr],
\qquad
\mathrm dx_t^u=b(x_t^u,u(t))\,\mathrm dt+\sigma(x_t^u,u(t))\,\mathrm dW_t .
$$

The paper sorts existing numerical routes into four: reduction to finite-dimensional stochastic programming; dynamic programming, that is solving the HJB equation, which it calls one of the most widely used methods; martingale-based methods; and methods based on the **stochastic maximum principle**. The gap it identifies is specific: while the stochastic maximum principle is a popular tool for theoretical work, it has not been widely used in the numerical setting. That is what this paper addresses.

One easily missed part of the setup: the control space $U=L^2([0,T];\mathbb R)$ consists of **deterministic** square-integrable controls, which the paper defends on the grounds that engineering and financial applications need this kind of advance planning; the adapted (feedback) case is handled separately in its Section 5.

### A fixed-point characterisation

The first-order optimality condition is the variational inequality $(J'(u^\ast),v-u^\ast)\ge0$ for all $v\in K$. Writing $P_K\omega=\arg\min_{u\in K}\|u-\omega\|$ for the projection, equivalently characterised by $(P_K\omega-\omega,v-P_K\omega)\ge0$, and comparing the two gives, for any $\rho>0$, the **fixed-point characterisation**

$$
u^\ast=P_K\bigl(u^\ast-\rho J'(u^\ast)\bigr).
$$

The algorithm is then the gradient projection iteration $u^{i+1,N}=P_{K_N}\bigl(u^{i,N}-\rho_iJ_N'(u^{i,N})\bigr)$, with the control space discretised into piecewise constants $U_N=\{\sum_n\alpha_n\chi_{I_n^N}\}$ and $J_N'$ a numerical approximation of $J'$. **This turns constrained optimal control into repeated gradient evaluation and projection, so the entire computational difficulty concentrates in one question: how to compute $J'$.**

### The gradient comes from an adjoint BSDE

Computing $J'$ directly needs the Gâteaux derivative $Dx_t^u(v)$, which is expensive. The paper instead introduces an adjoint process and applies Itô's formula to $p_t^uDx_t^u(v)$; the $Dx_t^u(v)$ terms cancel, leaving

$$
J'(u)\big|_t=\mathbb E\bigl[p_t^u\,b_u'(x_t^u,u(t))+q_t^u\,\sigma_u'(x_t^u,u(t))\bigr]+j'(u(t)),
$$

where $(p^u,q^u)$ solves the adjoint BSDE

$$
-\mathrm dp_t^u=f(x_t^u,p_t^u,q_t^u,u(t))\,\mathrm dt-q_t^u\,\mathrm dW_t,
\qquad p_T^u=g(x_T^u),
$$

with generator $f(x,p,q,u)=h'(x)+p\,b_x'(x,u)+q\,\sigma_x'(x,u)$. **Note that it is linear in $p$ and $q$** — which is why the scheme below, nominally implicit, needs no iteration.

The paper deliberately contrasts this with an earlier route in its Remark 1: there the adjoint equation is an **anticipating** stochastic differential equation whose solution must be backward-adapted rather than forward-adapted in the classical sense, and the paper observes that such a requirement does not hold in general, so its well-posedness is unclear. The BSDE above is well posed by standard theory. **That is a substantive methodological difference, not a matter of preference.**

### Discretisation and an error balance

The forward-backward system is discretised by Euler, the backward equation integrated over $[t_n,t_{n+1}]$ and conditioned, using the **left-point rectangle rule**; under the Euler state the conditional expectations become Gaussian integrals, evaluated by an $L$-point Gauss-Hermite rule with linear interpolation.

Convergence has two layers. Outer: if $J'$ is Lipschitz and uniformly monotone near the optimum and $\epsilon_N=\sup_i\|J'(u^{i,N})-J_N'(u^{i,N})\|\to0$, the iteration converges; if in addition $u^\ast$ and $J'(u^\ast)$ are Lipschitz, then $\epsilon_N\sim O(\Delta t)$ gives $\|u^\ast-u^{i,N}\|\sim O(\Delta t)$. Inner is the adjoint system's error,

$$
\hat{\mathbb E}\bigl[(\mu_n)^2\bigr]+\Delta t\sum_{n=0}^{N-1}\hat{\mathbb E}\bigl[(\nu_n)^2\bigr]
=O\bigl((\Delta t)^2\bigr)+O\bigl((\Delta x)^4/(\Delta t)^2\bigr),
$$

with $(\mu_n,\nu_n)$ the errors in $(p_n,q_n)$. **The form of the second term is worth remembering: the spatial interpolation error enters as $(\Delta x)^4/(\Delta t)^2$, so $\Delta x$ and $\Delta t$ cannot be refined independently** — balancing the two forces the spatial grid to be refined along with the time step, which is the characteristic price of interpolation-based conditional expectations.

Paper **41 (highly accurate schemes for stochastic optimal control via FBSDEs)** improves the same problem from the other side, bringing high-order FBSDE schemes into the control setting so each forward-backward solve is itself more accurate. **The division of labour between these two is a recurring pattern in this thread**: one paper improves the outer iteration, another improves the inner solve.

## 50 and 51: the edge of the family

- **50 (an efficient numerical algorithm for data-driven feedback control)** has control and **observation** together: the state satisfies $\mathrm dX_t=b(t,X_t,u_t)\mathrm dt+\sigma(t,X_t,u_t)\mathrm dW_t$ with cost $J(u)=\mathbb E[\int_0^Tf(t,X_t,u_t)\mathrm dt+h(X_T)]$, and there is a separate observation process $\mathrm dM_t=g(X_t)\mathrm dt+\mathrm dB_t$. That third equation is the difference from the preceding papers: the control law cannot come from the optimality conditions of a known model alone but must use the observations to update the state estimate, so **filtering and control become coupled**.
- **51 (a Gauss-Seidel type method for dynamic nonlinear complementarity problems)** treats

  $$
  \dot x(t)=F(t,x(t),y(t)),
  \qquad
  0\le y(t)\ \perp\ G(t,x(t),y(t))\ge0,
  $$

  with $x(t)\in\mathbb R^m$ and $y(t)\in\mathbb R^n_+$, including the differential semiaffine system and the dynamic linear complementarity problem as subclasses. Backward Euler in time leaves a coupled nonlinear system to solve at every time point, and the paper identifies definite defects in both mainstream approaches: **direct elimination** is valid only in the linear case, and after eliminating $x_j$ the reduced matrix $M_h=hN(I-hA)^{-1}B+M$ **may fail to be a P-matrix even when $M$ is**, while forming it requires $n$ large linear solves; **semismooth Newton** is only locally convergent and needs the Clarke generalized Jacobian, which is expensive and awkward at scale. The Gauss-Seidel-type block iteration targets exactly those two. Its connection to FBSDEs is looser, but it shares the theme of constrained dynamic optimisation with paper 26.

> [!note] Coverage status
> Papers 19, 25, 26, 50 and 51 have been checked against full preprint or author-accepted texts. The full texts of papers 16 and 41 could not be obtained — the publisher returns 403 to direct PDF download for the former and the PDF endpoint is unreachable for the latter, and neither has a preprint — so those two receive only their problem setting and position in the topic, with no schemes, theorems or convergence orders reported. Paper 16's constructional idea is inferred from the sister paper it cites, paper 19, and is flagged as such here.

## Where the seven sit

| No. | What it treats                                               | Position relative to the others                |
| --- | ------------------------------------------------------------ | ---------------------------------------------- |
| 16  | fully nonlinear parabolic equations                          | the starting point of the second-order setting |
| 19  | second-order FBSDEs and stochastic control                   | connects the second-order setting to control   |
| 25  | evaluation cost of multidimensional conditional expectations | improvement in the spatial direction           |
| 26  | constrained stochastic optimal control                       | the outer iteration (gradient projection)      |
| 41  | high-accuracy schemes for stochastic control                 | the inner solve (high-order FBSDE schemes)     |
| 50  | data-driven feedback control                                 | the control law comes from data, not a model   |
| 51  | dynamic nonlinear complementarity                            | block iteration under nonsmooth constraints    |

## Sources for this page

- T. Kong, W. Zhao, and T. Zhou, [_Probabilistic high order numerical schemes for fully nonlinear parabolic PDEs_](https://doi.org/10.4208/cicp.240515.280815a), Commun. Comput. Phys. 18 (2015), pp. 1482-1503.
- T. Kong, W. Zhao, and T. Zhou, [_High order numerical schemes for second-order FBSDEs with applications to stochastic optimal control_](https://doi.org/10.4208/cicp.OA-2016-0056), Commun. Comput. Phys. 21 (2017), pp. 808-834.
- Y. Fu, W. Zhao, and T. Zhou, [_Efficient spectral sparse grid approximations for solving multi-dimensional forward backward SDEs_](https://doi.org/10.3934/dcdsb.2017174), Discrete Contin. Dyn. Syst. Ser. B 22(9) (2017), pp. 3439-3458.
- B. Gong, W. Liu, T. Tang, W. Zhao, and T. Zhou, [_An efficient gradient projection method for stochastic optimal control problems_](https://doi.org/10.1137/17M1123559), SIAM J. Numer. Anal. 55(6) (2017), pp. 2982-3005.
- Y. Fu, W. Zhao, and T. Zhou, [_Highly accurate numerical schemes for stochastic optimal control via FBSDEs_](https://doi.org/10.4208/nmtma.OA-2019-0137), Numer. Math. Theor. Meth. Appl. 13 (2020), pp. 296-319.
- R. Archibald, F. Bao, J. Yong, and T. Zhou, [_An efficient numerical algorithm for solving data driven feedback control problems_](https://doi.org/10.1007/s10915-020-01358-y), J. Sci. Comput. 85(2) (2020), 58.
- S. Wu, T. Zhou, and X. Chen, [_A Gauss-Seidel type method for dynamic nonlinear complementarity problems_](https://doi.org/10.1137/19M1268884), SIAM J. Control Optim. 58(6) (2020), pp. 3389-3412.
