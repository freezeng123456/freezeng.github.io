---
title: Time-Fractional Phase Field
description: Papers 40, 43 and 57 - moving the energy law for a fractional gradient flow from integral to differential form
lang: en
translation: computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field
tags:
  - paper-notes
  - phase-field
  - fractional-derivatives
---

> [!note] Coverage of this page
> Papers **40** (_SIAM J. Sci. Comput._ 41(6), 2019), **43** (_J. Comput. Phys._ 414, 2020) and **57** (_SIAM J. Sci. Comput._ 43(5), 2021). The three form a clear progression: paper 40 proves an integral energy law and paper 57 replaces it by a differential law for a variational energy.

## 40: energy dissipation in the fractional case is an integral statement

### Where the standard argument fails

Phase-field models (Allen-Cahn, Cahn-Hilliard, molecular beam epitaxy) all come from a free energy that dissipates in time, and that dissipation law is the backbone of long-time numerical stability. Replacing the time derivative by a Caputo derivative of order $\alpha\in(0,1)$,

$$
\frac{\partial^{\alpha}}{\partial t^{\alpha}}\phi
={}^{C}_{0}D_{t}^{\alpha}\phi(t)
:=\frac{1}{\Gamma(1-\alpha)}\int_{0}^{t}\frac{\phi'(s)}{(t-s)^{\alpha}}\,\mathrm ds,
$$

breaks the standard argument: multiplying by $\phi_t$ no longer produces a sign-definite term pointwise in time, because $\partial_t^\alpha\phi$ is a **nonlocal-in-time** convolution of the whole history. Before this paper, the energy dissipation behaviour of time-fractional phase-field models had only been observed numerically, with no proof at either the continuous or the discrete level.

The paper closes that gap and identifies the **correct form** of the law: it is an **integral** statement accumulated over $[0,T]$, not a pointwise differential one. That is worth remembering, because it is exactly the starting point of paper 57.

The three models read

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}
=\gamma\Bigl(\varepsilon\Delta\phi-\frac{1}{\varepsilon}F'(\phi)\Bigr)
\qquad\text{(Allen-Cahn)},
$$

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}=\gamma\Delta\mu,
\quad \mu=-\varepsilon\Delta\phi+\frac{1}{\varepsilon}F'(\phi)
\qquad\text{(Cahn-Hilliard)},
$$

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}
=\gamma\Bigl(-\varepsilon\Delta^{2}\phi
+\frac{1}{\varepsilon}\nabla\cdot\mathbf f_m(\nabla\phi)\Bigr)
\qquad\text{(molecular beam epitaxy)},
$$

with free energies

$$
E[\phi]=\frac{\varepsilon}{2}\|\nabla\phi\|^{2}
+\frac{1}{\varepsilon}\bigl\langle F(\phi),1\bigr\rangle,
\qquad F(\phi)=\tfrac14(1-\phi^{2})^{2},
$$

$$
E_{m}[\phi]=\frac{\varepsilon}{2}\|\Delta\phi\|^{2}
+\frac{1}{\varepsilon}\bigl\langle F_m(\nabla\phi),1\bigr\rangle,
$$

where $F_m$ is $\frac14(|\mathbf v|^{2}-1)^{2}$ with slope selection and $-\frac12\ln|1+|\mathbf v|^{2}|$ without, giving $\mathbf f_m(\nabla\phi)=(|\nabla\phi|^{2}-1)\nabla\phi$ and $-\nabla\phi/(1+|\nabla\phi|^{2})$ respectively.

### The one analytical device: positivity of the fractional convolution kernel

What replaces the pointwise sign argument is a positivity result for the convolution kernel. For $h,g\in L^{p}(0,T)$ define

$$
I_\alpha(h,g):=\frac{1}{\Gamma(\alpha)}\int_{0}^{T}\!\!\int_{0}^{t}
\frac{h(s)g(t)}{(t-s)^{1-\alpha}}\,\mathrm ds\,\mathrm dt .
$$

The paper proves that for $\alpha\in(0,1)$ and $p\ge\frac{2}{1+\alpha}$, the kernel $I_\alpha(h,h)$ is positive. That converts "a weighted accumulation of history" into a signed quantity, so the energy argument can be completed at the **integral** level.

The device matters beyond this paper: it is the continuous version of a pattern that recurs throughout this series, whose discrete counterpart is the quadratic-form positivity on the [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|variable-step BDF page]]. **The difficulty and the remedy are isomorphic at the continuous and discrete levels: collect the nonlocal history into one positive-definite quantity.**

## 43: a second-order maximum-principle-preserving scheme on nonuniform steps

Paper 43 constructs a scheme for the time-fractional Allen-Cahn equation aiming simultaneously at second-order accuracy, nonuniform steps and a discrete maximum principle. Nonuniform steps are not optional here: fractional problems typically have an initial singularity near $t=0$ that needs a graded mesh, while the subsequent coarsening is extremely slow and needs growing steps.

Its position in the thread is as the bridge: it carries the continuous conclusion of paper 40 into the discrete setting and handles the maximum principle, while the question it leaves — what form the energy law takes discretely — is answered by paper 57.

## 57: changing the question

### Why the conclusion of paper 40 is not enough

Paper 40 proves $E[u](t)\le E[u](0)$ and the maximum bound principle, but that statement is weaker than the classical dissipation law

$$
\frac{\mathrm dE}{\mathrm dt}+\Bigl\|\frac{\delta E}{\delta u}\Bigr\|^2=0,
$$

and it is not a **differential** law, so it cannot be used the way the classical law is. Meanwhile practice needs **variable** steps (initial singularity plus coarsening), and no scheme was known that preserved energy stability and the maximum bound principle **simultaneously** under variable steps.

The resolution is to change the question: instead of asking whether the **original** energy dissipates, define a **variational energy** that does.

### The key step: move the nonlocality onto the variational derivative

The model is $\partial_t^{\alpha}u=\varepsilon^{2}\Delta u-f(u)$ with $F(u)=\tfrac14(1-u^2)^2$ and $f=F'$, read as a **fractional gradient flow** $\partial_t^{\alpha}u=-\delta E/\delta u$. Using the Riemann-Liouville derivative ${}^{R}\!\partial_t^{\alpha}v:=\partial_t\mathcal I_t^{1-\alpha}v$ and the semigroup identity

$$
{}^{R}\!\partial_t^{1-\alpha}\bigl(\partial_t^{\alpha}v\bigr)
=\partial_t\mathcal I_t^{1}v'=v',
$$

the equation is **equivalently rewritten** as

$$
\partial_t u=-{}^{R}\!\partial_t^{1-\alpha}\Bigl(\frac{\delta E}{\delta u}\Bigr).
$$

This is the core of the paper: the time derivative on the left is now **local**, so the standard energy test function applies, and the nonlocality has moved onto the variational derivative, where a positivity property is available.

### The variational energy and its differential law

Using the Riemann-Liouville inequality

$$
v(t)\bigl({}^{R}\!\partial_t^{1-\alpha}v\bigr)(t)
\ \ge\ \tfrac12\bigl({}^{R}\!\partial_t^{1-\alpha}v^{2}\bigr)(t)
+\tfrac12\omega_{\alpha}(t)v^{2}(t),
\qquad \omega_{\mu}(t)=\frac{t^{\mu-1}}{\Gamma(\mu)},
$$

the paper defines a variational energy and proves its law:

$$
\mathcal{E}_{\alpha}[u]:=E[u]+\frac12\,\mathcal{I}_t^{\alpha}
\Bigl\|\frac{\delta E}{\delta u}\Bigr\|^{2},
\qquad
\frac{\mathrm d\mathcal{E}_{\alpha}}{\mathrm dt}
+\frac12\,\omega_{\alpha}(t)\Bigl\|\frac{\delta E}{\delta u}\Bigr\|^{2}\le0,
\qquad \forall t>0 .
$$

As $\alpha\to1$, $\mathcal I_t^\alpha\to\mathcal I_t^1$ and $\omega_\alpha(t)\to1$, recovering the classical law. The paper calls this **asymptotically energy-dissipation preserving**: the new law is not an analogue of the classical one but a generalisation that degenerates back to it.

### The L1$_R$ formula and positivity without any step-ratio restriction

On variable steps the Riemann-Liouville derivative is discretised as

$$
\bigl({}^{R}\!\partial_{\tau}^{1-\alpha}v\bigr)^{n-\frac12}
:=\frac{1}{\tau_{n}}\int_{t_{n-1}}^{t_{n}}\frac{\partial}{\partial t}
\int_{0}^{t}\omega_{\alpha}(t-s)(\Pi_{0}v)(s)\,\mathrm ds\,\mathrm dt
\ \triangleq\ \frac{1}{\tau_{n}}\sum_{k=1}^{n}a_{n-k}^{(n)}v^{k-\frac12},
$$

with $\Pi_0v$ the piecewise constant interpolant equal to $v^{k-\frac12}$ on $(t_{k-1},t_k]$. The auxiliary sequence and kernels are

$$
q_{n-k}^{(n)}:=\int_{t_{k-1}}^{t_{k}}\omega_{\alpha}(t_{n}-s)\,\mathrm ds
=\sum_{j=k}^{n}a_{j-k}^{(j)}>0,
$$

$$
a_{0}^{(n)}:=q_{0}^{(n)}>0\ (n\ge1),
\qquad
a_{n-k}^{(n)}:=q_{n-k}^{(n)}-q_{n-k-1}^{(n-1)}<0\ (n\ge k+1\ge2),
$$

a sign pattern worth noting: the first kernel is positive and all later ones negative, with $a^{(n)}_0=\omega_{1+\alpha}(\tau_n)=\tau_n^{\alpha}/\Gamma(1+\alpha)$ explicitly.

Kernel positivity is proved discretely rather than inherited from the continuous kernel: for any real sequence $\{w_k\}$,

$$
2\sum_{k=1}^{n}w_{k}\sum_{j=1}^{k}a_{k-j}^{(k)}w_{j}
\ \ge\ \sum_{k=1}^{n}\Bigl(q_{n-k}^{(n)}+\sum_{j=1}^{k}a_{k-j}^{(k)}\Bigr)w_{k}^{2}
\ >\ 0,
\qquad n\ge1,\ w\not\equiv0 .
$$

Two identities carry the proof: complete monotonicity of $\omega_\alpha$ gives $q_{k-j-1}^{(k-1)}-q_{k-j}^{(k)}>0$, and

$$
\sum_{j=1}^{k}a_{k-j}^{(k)}
=\int_{t_{k-1}}^{t_{k}}\omega_{\alpha}(s)\,\mathrm ds>0 .
$$

**This positivity holds with no step-ratio restriction whatsoever.** The contrast with the [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|variable-step BDF]] case is sharp: there the quadratic form yields thresholds such as $r_k<(3+\sqrt{17})/2$, while here the fractional kernels are positive definite directly by complete monotonicity. A step-size bound appears only in the maximum-bound argument, not in the energy law.

### The scheme and two constructive designs

Setting $v:=-\delta E/\delta u$ splits the equation as

$$
\partial_t u={}^{R}\!\partial_t^{1-\alpha}v,
\qquad
v=\varepsilon^{2}\Delta u-f(u),
$$

discretised in Crank-Nicolson form,

$$
\partial_{\tau}u^{n-\frac12}=\bigl({}^{R}\!\partial_{\tau}^{1-\alpha}v\bigr)^{n-\frac12},
\qquad
v^{n-\frac12}=\varepsilon^{2}D_{h}u^{n-\frac12}-H(u^{n},u^{n-1}),
$$

with the second-order nonlinear approximation

$$
H(u^{n},u^{n-1}):=\tfrac13(u^{n})^{.3}+\tfrac12(u^{n-1})^{.2}\!\circ u^{n}
+\tfrac16(u^{n-1})^{.3}-\tfrac12\bigl(u^{n}+u^{n-1}\bigr)
$$

where $\circ$ and the powers are element-wise. This particular $H$ is **engineered** so that $H(a,b)(a-b)\ge F(a)-F(b)$ holds pointwise, the discrete counterpart of the chain rule and the reason the energy argument can proceed.

The discrete variational energy is

$$
\mathcal{E}_{\alpha}[u^{n}]:=E[u^{n}]+\frac12h^{2}\sum_{i,j}\sum_{k=1}^{n}
q_{n-k}^{(n)}\bigl(v_{ij}^{k-\frac12}\bigr)^{2},
$$

where the $q$ kernels themselves form a numerical fractional integral, $(\mathcal I_\tau^{\alpha}v)^n=\sum_{k}q_{n-k}^{(n)}v^{k-\frac12}$.

### DOC kernels and a reversible transformation between the two fractional derivatives

Another stated first in the paper is an explicit equivalence between the discretisations of the two fractional derivatives. Define the discrete orthogonal convolution kernels $\theta$ of the L1$_R$ kernels by the orthogonality identity

$$
\sum_{j=k}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}\equiv\delta_{nk},
\qquad
\theta_{0}^{(n)}=\frac{1}{a^{(n)}_{0}} .
$$

Orthogonality holds in **both** directions, which is the reversibility:

$$
\sum_{j=k}^{n}a^{(n)}_{n-j}\theta_{j-k}^{(j)}\equiv\delta_{nk},
\qquad
\sum_{j=k}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}\equiv\delta_{nk},
$$

with complementarity $\sum_{j=k}^{n}q_{n-j}^{(n)}\theta_{j-k}^{(j)}\equiv1$ between $\theta$ and $q$. The paper proves $\theta^{(n)}_0=\Gamma(1+\alpha)\tau_n^{-\alpha}$, that all $\theta^{(n)}_j>0$, and monotone decrease $\theta_{0}^{(n)}>\theta_{1}^{(n)}>\cdots>\theta_{n-1}^{(n)}>0$ for $n\ge2$.

Convolving the first scheme equation with $\theta^{(n)}_{n-j}$ and using orthogonality gives

$$
\sum_{j=1}^{n}\theta_{n-j}^{(n)}\nabla_{\tau}u^{j}
=\sum_{k=1}^{n}v^{k-\frac12}\sum_{j=k}^{n}\theta_{n-j}^{(n)}a_{j-k}^{(j)}
=v^{n-\frac12},
$$

so the Riemann-Liouville-form Crank-Nicolson scheme is **equivalent** to the Caputo-type form

$$
\sum_{j=1}^{n}\theta_{n-j}^{(n)}\nabla_{\tau}u^{j}
=\varepsilon^{2}D_{h}u^{n-\frac12}-H(u^{n},u^{n-1}) .
$$

The DOC kernels therefore define a **new** discrete Caputo derivative whose kernels are positive and monotonically decreasing on nonuniform meshes, just like the classical L1 kernels. The paper cautions that this is an **indirect** approximation whose accuracy differs from the direct L1 formula, which has error order $2-\alpha$.

The maximum-bound argument runs on this Caputo-equivalent form, by induction, using the symmetry and negative semi-definiteness of $D_h$ together with an $\ell^\infty$ lemma.

## The progression across the three papers

| No. | Form of the energy law                 | Step restriction             | Main analytical device                          |
| --- | -------------------------------------- | ---------------------------- | ----------------------------------------------- |
| 40  | integral, accumulated over $[0,T]$     | none at the continuous level | positivity of the fractional convolution kernel |
| 43  | discrete maximum principle             | nonuniform steps             | second-order nonuniform discretisation          |
| 57  | differential, for a variational energy | none for the energy law      | rewriting plus L1$_R$ kernel positivity         |

The move from paper 40 to paper 57 deserves its own summary: **when a law comes out in the wrong form, do not weaken the conclusion, change the object.** Paper 40 obtains an integral inequality for the original energy; paper 57 does not try to strengthen it but constructs a new energy $\mathcal E_\alpha$ satisfying a differential law that degenerates to the classical one as $\alpha\to1$. The price is that this energy contains a fractional integral term and is therefore not the original energy.

## Coverage check

| Item                                                   | Paper | Status                                                                    |
| ------------------------------------------------------ | ----- | ------------------------------------------------------------------------- |
| Three fractional phase-field models and their energies | 40    | Allen-Cahn, Cahn-Hilliard, both MBE variants                              |
| Why the standard energy argument fails                 | 40    | loss of pointwise sign, nonlocal convolution                              |
| Positivity of the fractional convolution kernel        | 40    | definition, conditions, role in the argument                              |
| Integral rather than differential form                 | 40    | the shape of the conclusion and its consequence                           |
| Second-order nonuniform maximum-principle scheme       | 43    | goal and position in the thread                                           |
| Rewriting and relocating the nonlocality               | 57    | semigroup identity, rewritten equation, why the test function works       |
| Variational energy and differential law                | 57    | Riemann-Liouville inequality, $\mathcal E_\alpha$, the $\alpha\to1$ limit |
| L1$_R$ kernels and their positivity                    | 57    | definitions, sign pattern, discrete positivity, no step-ratio bound       |
| Constructive design of the nonlinear term              | 57    | form of $H$ and the inequality it satisfies                               |
| DOC kernels and the reversible transformation          | 57    | orthogonality, complementarity, monotonicity, equivalent Caputo form      |

## Sources for this page

- T. Tang, H. Yu, and T. Zhou, [_On energy dissipation theory and numerical stability for time-fractional phase-field equations_](https://doi.org/10.1137/18M1203560), SIAM J. Sci. Comput. 41(6) (2019), pp. A3757-A3778 (preprint [arXiv:1808.01471](https://arxiv.org/abs/1808.01471)).
- H.-l. Liao, T. Tang, and T. Zhou, [_A second-order and nonuniform time-stepping maximum-principle preserving scheme for time-fractional Allen-Cahn equations_](https://doi.org/10.1016/j.jcp.2020.109473), J. Comput. Phys. 414 (2020), 109473.
- H.-l. Liao, T. Tang, and T. Zhou, [_An energy stable and maximum bound preserving scheme with variable time steps for time fractional Allen-Cahn equation_](https://doi.org/10.1137/20M1384105), SIAM J. Sci. Comput. 43(5) (2021), pp. A3503-A3526 (preprint [arXiv:2012.10740](https://arxiv.org/abs/2012.10740)).
