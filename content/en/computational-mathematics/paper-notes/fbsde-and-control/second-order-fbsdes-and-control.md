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

## 25: the number of evaluation points for a conditional expectation

In several dimensions a conditional expectation is a multidimensional Gaussian integral, and tensor-product quadrature makes the number of evaluation points grow exponentially with dimension. Paper 25 addresses that growth with **spectral sparse grids**, under the published title _Efficient spectral sparse grid approximations for solving multi-dimensional forward backward SDEs_.

Its position is worth naming: in this thread the temporal accuracy problem is solved by multistep schemes and deferred correction, while the cost of the **spatial** representation is solved by sparse grids and later by Sinc approximation ([[en/computational-mathematics/paper-notes/fbsde-and-control/stability-theory-for-fbsdes|paper 63]]). The two directions of improvement are independent.

## 26: organising the control iteration as a gradient projection

Paper 26 treats stochastic optimal control. The stochastic maximum principle writes the optimality condition as a forward-backward (Pontryagin) system, and the control variable is usually constrained, for instance to a convex set. The natural iteration is then a **gradient projection**: descend along the Gâteaux derivative of the cost with respect to the control, then project back onto the feasible set.

The technical focus is making that iteration efficient: each gradient evaluation requires one solve of the forward-backward system, so both the iteration count and the per-iteration cost must be controlled. Paper **41 (highly accurate schemes for stochastic optimal control via FBSDEs)** improves the same problem from the other side, bringing high-order FBSDE schemes into the control setting so each forward-backward solve is itself more accurate.

**The division of labour between these two is a recurring pattern in this thread**: one paper improves the outer iteration, another improves the inner solve.

## 50 and 51: the edge of the family

- **50 (an efficient numerical algorithm for data-driven feedback control)** computes feedback control in a data-driven setting. The difference from the preceding papers is that the control law is no longer derived from the optimality conditions of a known model but must be obtained from data, so filtering and control become coupled.
- **51 (a Gauss-Seidel type method for dynamic nonlinear complementarity problems)** treats dynamic problems with complementarity constraints. Complementarity makes the problem nonsmooth, so gradient-type iterations no longer apply directly and a Gauss-Seidel-type block iteration is used instead. Its connection to FBSDEs is looser, but it shares the theme of constrained dynamic optimisation with paper 26.

> [!note] Coverage status
> This page has not checked papers 16, 19, 25, 26, 41, 50 and 51 equation by equation. What it gives is the problem setting of each and its position within the topic, without expanding the schemes, theorems and convergence orders. Most of these appeared in journals with no preprint, and their texts require a subscription.

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
