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
> Papers **40** (_SIAM J. Sci. Comput._ 41(6), 2019), **43** (_J. Comput. Phys._ 414, 2020) and **57** (_SIAM J. Sci. Comput._ 43(5), 2021). All three were checked equation by equation against the authors' own arXiv LaTeX sources, so the equations, theorem hypotheses and constants below are transcribed rather than paraphrased; all three carry a full derivation chain, theorems with hypotheses, and numerical experiments.
>
> The three form a clear progression: paper 40 proves an integral energy law and states plainly that it cannot be strengthened, paper 43 supplies the missing second-order maximum-principle-preserving scheme on nonuniform meshes, and paper 57 changes the energy object and gets a differential law in return.

## 40: energy dissipation in the fractional case is an integral statement

### The idea

Phase-field models (Allen-Cahn, Cahn-Hilliard, molecular beam epitaxy) all descend from a free energy that dissipates in time, and that dissipation law is the backbone of long-time numerical stability. In the integer-order case the argument is one line: test the equation with $\phi_t$, the right-hand side becomes $-\|\phi_t\|^2$, which is negative **pointwise**, so $\mathrm dE/\mathrm dt\le0$.

Replace the time derivative by a Caputo derivative and that line breaks immediately. Now $\partial_t^\alpha\phi$ is a weighted convolution of the entire history, and $\int_\Omega\phi_t\,\partial_t^\alpha\phi$ can carry either sign at any given instant: the present rate of change is coupled to past rates of change through the kernel $(t-s)^{-\alpha}$, and there is no reason for the two to agree in sign pointwise. Before this paper, energy dissipation for time-fractional phase-field models had only been observed numerically, with no proof at either the continuous or the discrete level.

The mechanism of the paper is to **move the sign property from pointwise to accumulated**. The kernel $|t-s|^{-\alpha}$ is positive definite as a quadratic form on $[0,T]^2$ — an elementary instance of Riesz-potential positivity, visible through the Fourier transform. So although $\phi_t\cdot\partial_t^\alpha\phi$ has no sign at a single instant, integrating it once over $[0,T]$ produces a double integral that is a positive quantity. **The price is that the resulting law can only be of integral type: it gives $E[\phi(T)]\le E[\phi(0)]$ but no pointwise rate of decrease.** That is not a technical shortfall but the true shape of a fractional gradient flow, and the numerical experiments of paper 57 display it directly.

> [!warning] The law must not be quoted in a strengthened form
> The paper states explicitly that the integral energy law it obtains **does not imply** $\frac{\mathrm d}{\mathrm dt}E\le0$, and **does not imply** $\frac{\mathrm d^{\alpha}}{\mathrm dt^{\alpha}}E\le0$ either. This is the single most frequently mis-stated point in the downstream literature: restating the conclusion as a pointwise dissipation law, whether in the integer-order or the fractional derivative, does not hold, and the paper lists both as open problems in its concluding remarks. Precisely because no such strengthening is available, paper 57 has to **change the energy object** rather than push this law further.

### Setting

The Caputo derivative is taken in its standard form,

$$
\frac{\partial^{\alpha}}{\partial t^{\alpha}}\phi
={}^{C}_{0}D_{t}^{\alpha}\phi(t)
:=\frac{1}{\Gamma(1-\alpha)}\int_{0}^{t}\frac{\phi'(s)}{(t-s)^{\alpha}}\,\mathrm ds,
\qquad t>0,\ \alpha\in(0,1).
$$

The three models read (with $\gamma$ the mobility, $\varepsilon$ the interfacial width, and homogeneous Dirichlet, Neumann or periodic boundary conditions)

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}
=\gamma\Bigl(\varepsilon\Delta\phi-\frac{1}{\varepsilon}F'(\phi)\Bigr)
\qquad\text{(Allen-Cahn)},
$$

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}=\gamma\Delta\mu,
\quad \mu=-\varepsilon\Delta\phi+\frac{1}{\varepsilon}F'(\phi)
\qquad\text{(Cahn-Hilliard, periodic or no-flux }\partial_n\mu=\partial_n\phi=0\text{)},
$$

$$
\frac{\partial^{\alpha}\phi}{\partial t^{\alpha}}
=\gamma\Bigl(-\varepsilon\Delta^{2}\phi
+\frac{1}{\varepsilon}\nabla\cdot\mathbf f_m(\nabla\phi)\Bigr)
\qquad\text{(molecular beam epitaxy, }\mathbf f_m=\partial F_m/\partial\mathbf v\text{)},
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

where $\langle\cdot,\cdot\rangle$ and $\|\cdot\|$ are the $L^2(\Omega)$ inner product and norm. The MBE potential $F_m$ comes in two variants according to whether slope selection is present:

$$
F_{m}(\mathbf v)=
\begin{cases}
\tfrac14\bigl(|\mathbf v|^{2}-1\bigr)^{2}, & \text{with slope selection},\\[2pt]
-\tfrac12\ln\bigl|1+|\mathbf v|^{2}\bigr|, & \text{without slope selection},
\end{cases}
$$

giving $\mathbf f_m(\nabla\phi)=(|\nabla\phi|^{2}-1)\nabla\phi$ and $-\nabla\phi/(1+|\nabla\phi|^{2})$ respectively. The two variants live on different pages of this topic: the case without slope selection is handled by [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|paper 52]], the case with slope selection by [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation|paper 78]].

### Derivation

**Step one: positivity of the continuous kernel.** For $h,g\in L^{p}(0,T)$ define

$$
I_\alpha(h,g):=\frac{1}{\Gamma(\alpha)}\int_{0}^{T}\!\!\int_{0}^{t}
\frac{h(s)g(t)}{(t-s)^{1-\alpha}}\,\mathrm ds\,\mathrm dt .
$$

The lemma proved in the paper is that for $\alpha\in(0,1)$ and $p\ge\frac{2}{1+\alpha}$,

$$
I_\alpha(h,h)=\int_0^T\bigl(I_{0+}^{\alpha}h\bigr)(t)\,h(t)\,\mathrm dt
\ \ge\ \cos\frac{\alpha\pi}{2}\,
\bigl\|I_{0+}^{\alpha/2}h\bigr\|^{2}_{L^2(0,T)}\ \ge\ 0 .
$$

What actually gets used is its corollary at $\alpha\mapsto1-\alpha$. Setting $A_\alpha(h,g):=I_{1-\alpha}(h,g)$, one has for $p\ge\frac{2}{2-\alpha}$

$$
A_\alpha(h,h)
=\frac{1}{2}\cdot\frac{1}{\Gamma(1-\alpha)}\int_{0}^{T}\!\!\int_{0}^{T}
\frac{h(s)h(t)}{|t-s|^{\alpha}}\,\mathrm ds\,\mathrm dt
\ \ge\ \sin\frac{\alpha\pi}{2}\,
\bigl\|I_{0+}^{(1-\alpha)/2}h\bigr\|^{2}_{L^2(0,T)}\ \ge\ 0 .
$$

**The symmetrisation identity on the first line is the crucial algebraic step**: a Volterra-type double integral over a triangle is replaced by a symmetric (Riesz-type) double integral over the whole square with kernel $|t-s|^{-\alpha}$, which is what makes the phrase "quadratic form" meaningful at all.

The proof extends $h$ from $[0,T]$ to $\mathbb R$ by zero, uses the Liouville semigroup property $I_{+}^{\alpha}I_{+}^{\beta}=I_{+}^{\alpha+\beta}$ to write the integrand as $(I_{+}^{\alpha/2}h)(I_{-}^{\alpha/2}h)$, then applies

$$
\bigl(\mathcal F I_{\pm}^{\alpha}f\bigr)(\xi)=\frac{(\mathcal Ff)(\xi)}{(\mp i\xi)^{\alpha}},
\qquad
(\mp i\xi)^{\alpha}=|\xi|^{\alpha}e^{\mp i\alpha\pi\,\mathrm{sgn}(\xi)/2}
$$

together with Parseval, and finally invokes the Hardy-Littlewood inequality $\|I_{0+}^\alpha f\|_{L^q}\le K\|f\|_{L^p}$ with $q=p/(1-\alpha p)$ to secure integrability. **The constants $\cos\frac{\alpha\pi}{2}$ and $\sin\frac{\alpha\pi}{2}$ are exactly the real parts of $e^{\mp i\alpha\pi/2}$**, which also explains why the first bound degenerates as $\alpha\to1$ (since $\cos\frac\pi2=0$) while the second tends to $1$.

**Step two: the continuous energy law.** Test the equation with $\phi_t$ and integrate over $[0,T]$. On the right the elliptic and potential terms assemble into $E[\phi(T)]-E[\phi(0)]$ in the usual way; on the left, $\frac1\gamma\int_\Omega\int_0^T\phi_t\,\partial_t^\alpha\phi$ is precisely $\frac1\gamma\int_\Omega A_\alpha(\phi_t,\phi_t)\,\mathrm dx$, nonnegative by the previous step.

**Step three: L1 discretisation and complete monotonicity of the discrete kernels.** On a uniform step $\tau=T/n$ the L1 (piecewise-linear Caputo) kernels are

$$
b_{j}=\frac{1}{\Gamma(1-\alpha)}\int_{j\tau}^{(j+1)\tau}\frac{\mathrm dt}{t^{\alpha}}
=\frac{\tau^{1-\alpha}}{\Gamma(2-\alpha)}
\bigl[(j+1)^{1-\alpha}-j^{1-\alpha}\bigr],\qquad j\ge0 .
$$

(The original LaTeX source misprints $\tau$ as $\pi$ in the limits of integration; the closed form on the right is unambiguous and is what is actually used.) These kernels satisfy

$$
b_{k}>0,\qquad b_{k}-b_{k+1}>0,\qquad
\sum_{j=0}^{k-1}(b_{j}-b_{j+1})+b_{k}=b_{0},
$$

together with the rearrangement identity

$$
\sum_{j=0}^{k}b_{j}\frac{u^{k+1-j}-u^{k-j}}{\tau}
=\frac{1}{\tau}\Bigl[b_{0}u^{k+1}
-\sum_{j=0}^{k-1}(b_{j}-b_{j+1})u^{k-j}-b_{k}u^{0}\Bigr].
$$

**This rearrangement identity is the hinge of every discrete argument in the paper**: it separates the coefficient of the current level from nonnegative weights on the history, the history weights $b_j-b_{j+1}$ and $b_k$ all being positive and summing exactly to $b_0$. The stabilised L1 scheme (for Allen-Cahn, with $f=F'$ and stabilisation constant $S>0$) is

$$
\frac{1}{\gamma}\sum_{j=0}^{k}b_{j}\frac{\phi^{k+1-j}-\phi^{k-j}}{\tau}
=\varepsilon\Delta\phi^{k+1}-\frac{1}{\varepsilon}f(\phi^{k})
-\frac{S}{\gamma}\bigl(\phi^{k+1}-\phi^{k}\bigr),
$$

and the Cahn-Hilliard and MBE schemes are structurally identical, with the stabilisation placed in the chemical potential and acting through $\Delta$ respectively.

**Step four: positive definiteness of the discrete quadratic form.** This is the discrete counterpart of the corollary above. For any $(u_1,\dots,u_n)^T\in\mathbb R^{n}$ set

$$
B:=2\sum_{k=1}^{n}\sum_{j=1}^{k}b_{|k-j|}u_{j}u_{k}
=\sum_{k=1}^{n}b_{0}u_{k}^{2}+\sum_{k=1}^{n}\sum_{j=1}^{n}b_{|k-j|}u_{j}u_{k}
\ \ge\ \sum_{k=1}^{n}b_{0}u_{k}^{2},
$$

and, more sharply,

$$
B\ \ge\ \frac{2}{\tau}\,\sin\frac{\alpha\pi}{2}\,
\bigl\|I_{0+}^{(1-\alpha)/2}u^{n}(t)\bigr\|^{2}_{L^{2}(0,T)}
+s_{n}\sum_{k=1}^{n}u_{k}^{2},
\qquad
s_{n}=\Bigl(\frac{n+1}{2}\Bigr)^{-\alpha}\frac{\tau^{1-\alpha}}{\Gamma(1-\alpha)}>0,
$$

where $u^{n}(t)=u_{\lfloor t/\tau\rfloor+1}$. The proof compares $b_{|k|}$ against the exactly-integrated kernel

$$
\tilde b_{|k|}=\frac{\tau^{1-\alpha}}{\Gamma(3-\alpha)}
\bigl((k+1)^{2-\alpha}-2k^{2-\alpha}+(k-1)^{2-\alpha}\bigr)\ (k\ge1),
\qquad
\tilde b_{0}=\frac{2\tau^{1-\alpha}}{\Gamma(3-\alpha)}
$$

via $\tilde b_{|k|}-b_{|k|}\ge0$ for $k\ge1$ and $2b_{0}-\tilde b_{0}\ge0$, so that the difference matrix $C=\{b_{|k-j|}-\tilde b_{|k-j|}\}$ is a symmetric positive definite $M$-matrix; diagonal dominance then gives the bound above with $s_n=2b_0-\tilde b_0-c_0$.

**Step five: the maximum principle argument is comparison plus induction, not an energy argument.** This is worth stressing, because the later Liao-Tang-Zhou papers all switch to DOC kernels and this one does not. Assume the potential satisfies $F\in C^2(\mathbb R)$ and that there exist $M_1<0<M_2$ with

$$
F'(M_1)=F'(M_2)=0;\qquad F'(u)>0\ (u>M_2);\qquad F'(u)<0\ (u<M_1),
$$

which the quartic double well meets with $M_1=-1$, $M_2=1$. Use the rearrangement identity to write the scheme as

$$
\Bigl(\frac{b_0}{\gamma\tau}+\frac{S}{\gamma}\Bigr)\phi^{k+1}-\varepsilon\Delta\phi^{k+1}
=\frac{1}{\gamma\tau}\Bigl[\sum_{j=0}^{k-1}(b_{j}-b_{j+1})\phi^{k-j}+b_{k}\phi^{0}\Bigr]
+\frac{S}{\gamma}\phi^{k}-\frac{1}{\varepsilon}f(\phi^{k}),
$$

then write $f(\phi^k)=f(\phi^k)-f(M_2)=f'(\xi)(\phi^k-M_2)$ and apply the elliptic maximum principle to the operator $\bigl(\frac{b_0}{\gamma\tau}+\frac{S}{\gamma}\bigr)I-\varepsilon\Delta$. **Nonnegativity of the history coefficients is the single reason this step goes through**, and that is exactly the complete monotonicity of the kernels.

### Theorems

**(Continuous, time-fractional Allen-Cahn.)** If $E[\phi(0)]$ is finite, then

$$
E[\phi(T)]-E[\phi(0)]
=-\frac{1}{\gamma}\int_{\Omega}A_{\alpha}(\phi_{t},\phi_{t})\,\mathrm dx\ \le\ 0 .
$$

Note that this is an **equality** plus a sign determination. It simultaneously shows that the dissipated quantity $\frac1\gamma\int_\Omega A_\alpha(\phi_t,\phi_t)\mathrm dx$ is controlled by $E[\phi(0)]$, hence that $\|I_{0+}^{(1-\alpha)/2}\phi_t\|^2_{L^2(0,T)}$ is bounded — information strictly stronger than "the energy does not increase".

**(Continuous, time-fractional Cahn-Hilliard.)** After proving conservation of total mass,

$$
E[\phi(T)]-E[\phi(0)]
=-\frac{1}{\gamma}\int_{\Omega}A_{\alpha}(\nabla\psi,\nabla\psi)\,\mathrm dx\le0,
\qquad -\Delta\psi=\phi_{t}.
$$

**(Continuous, time-fractional MBE.)** As printed, $E_{m}[\phi(T)]-E_{m}[\phi(0)]\le-\frac1\gamma\int_\Omega A_\alpha(\phi,\phi)\,\mathrm dx\le0$.

> [!warning] A suspected misprint in the source
> In the authors' accepted-version LaTeX source, the two arguments of $A_\alpha$ in the MBE theorem are written $(\phi,\phi)$, whereas by analogy with the Allen-Cahn and Cahn-Hilliard statements and by the structure of the proof they should be $(\phi_t,\phi_t)$. This site judges it a typographical slip in the source, but the **official SIAM published version was not checked**, so anyone quoting this particular display must go back to the published PDF.

**(Discrete, fractional Allen-Cahn.)** Suppose the potential is modified so that its second derivative is globally bounded, $\max_u|F''(u)|\le L$; the paper uses the quadratically-growing truncation

$$
F(\phi)=\begin{cases}
\tfrac{11}{2}(\phi-2)^{2}+6(\phi-2)+\tfrac94, & \phi>2,\\
\tfrac14(\phi^{2}-1)^{2}, & \phi\in[-2,2],\\
\tfrac{11}{2}(\phi+2)^{2}-6(\phi+2)+\tfrac94, & \phi<-2 .
\end{cases}
$$

If

$$
S+\frac{b_{0}}{2\tau}\ \ge\ \frac{\gamma L}{2\varepsilon},
$$

then, writing $\delta_{t}\phi^{k+1}:=\phi^{k+1}-\phi^{k}$,

$$
E[\phi^{n}]-E[\phi^{0}]
\le-\frac{b_{0}}{2\gamma\tau}\sum_{k=0}^{n-1}\|\delta_{t}\phi^{k+1}\|^{2}
-\sum_{k=0}^{n-1}\Bigl\{\frac{\varepsilon}{2}\|\nabla\delta_{t}\phi^{k+1}\|^{2}
+\Bigl\langle\frac{S}{\gamma}-\frac{1}{2\varepsilon}f'(\xi^{k}),
(\delta_{t}\phi^{k+1})^{2}\Bigr\rangle\Bigr\}.
$$

In particular, **taking $S\ge\gamma L/(2\varepsilon)$ makes the scheme unconditionally energy stable**, with $E[\phi^{n}]\le E[\phi^{0}]$ for all $\tau>0$ and all $n>0$.

**(Discrete, fractional Cahn-Hilliard.)** The same structure, but with the dissipation measured in the $H^{-1}$ norm (reflecting the $H^{-1}$ gradient flow); the condition is $\sqrt{b_{0}\varepsilon/(\gamma\tau)}+S/\gamma\ge L/(2\varepsilon)$, again unconditional once $S\ge\gamma L/(2\varepsilon)$.

**(Discrete, fractional MBE.)** The condition reads $\sqrt{b_{0}\varepsilon/(\gamma\tau)}+S/\gamma\ge\frac{1}{2\varepsilon}\lambda_{\max}(\mathbf f_m'(\xi^{k}))$. For the model **without slope selection** this eigenvalue can be computed outright:

$$
\mathbf f_m'(\mathbf v)=\frac{2\mathbf v^{2}-(|\mathbf v|^{2}+1)I}{(1+|\mathbf v|^{2})^{2}},
\qquad
\lambda_{\max}\bigl(\mathbf f_m'\bigr)\le\frac18,
$$

so **$S\ge\dfrac{\gamma}{16\varepsilon}$ already suffices for unconditional energy stability**. The $1/8$ and the $\gamma/(16\varepsilon)$ are the exact constants printed in the paper.

**(Discrete maximum principle, fractional Allen-Cahn.)** Let $\phi_0\in C^0$ with $M_1\le\phi_0\le M_2$, using the **untruncated** standard double-well potential. If

$$
\frac{b_{0}-b_{1}}{\tau}+S\ \ge\ \gamma\,\frac{\max_{M_1\le u\le M_2}|f'(u)|}{\varepsilon},
$$

then $M_{1}\le\phi^{k}(x)\le M_{2}$ for all $k\ge1$. Since

$$
b_0-b_1=\frac{\tau^{1-\alpha}}{\Gamma(2-\alpha)}\bigl(2-2^{1-\alpha}\bigr),
$$

this is a **mild but genuine step-size restriction**, not an unconditional conclusion. Its value is that once the maximum principle holds, the global Lipschitz assumption on $f$ can be dropped and the truncated potential above is no longer needed.

### Numerical experiments

The spatial discretisation is a periodic Fourier-Galerkin method on $\Omega=[0,L_x]\times[0,L_y]$. Evaluating the history term directly would require all $n$ previous values at level $n$, at total cost $O(N^2)$; the paper instead uses the **sum-of-exponentials (SOE) fast algorithm** of Jiang et al. for the history part of the fractional derivative, which is the implementation prerequisite that makes long-time coarsening runs feasible.

| Test | Model                             | Fractional parameters  | Initial data and observable            |
| ---- | --------------------------------- | ---------------------- | -------------------------------------- |
| 1    | time-fractional Allen-Cahn        | $\alpha=1,\,0.5,\,0.3$ | random data; snapshots, discrete energy |
| 2    | time-fractional Cahn-Hilliard     | $\alpha=1,\,0.5,\,0.3$ | random data; discrete energy, coarsening rate |
| 3    | time-fractional MBE (with slope selection) | $\alpha=1,\,0.7,\,0.4$ | random data; discrete energy, coarsening rate |

Two conclusions. First, the computed energies decay **monotonically** in all three models, in agreement with the discrete theorems above. Second — and this is what the numerical section really produces — in the coarsening stage the energy dissipation rate obeys a power law with asymptotic exponent

$$
-\frac{\alpha}{3}
$$

for **both** the time-fractional Cahn-Hilliard equation and the time-fractional MBE model. At $\alpha=1$ it degenerates to the classical $-1/3$ coarsening law, which serves as a consistency check.

**This $-\alpha/3$ law lies outside the paper's own theory.** The authors say so directly in the concluding remarks: it is an empirical observation for which "a rigorous theoretical justification is needed". It is also the largest gap between the numerical and the theoretical parts of this paper: the theorems guarantee that the energy does not increase, while the experiments measure the **rate** at which it decreases, which the theorems say nothing about.

Two further limitations are worth recording: the scheme is only **first order** in time and the mesh is **uniform**, so the initial singularity of the fractional solution near $t=0$ is not addressed and the formal order is not actually attained. Those two points are precisely the starting point of paper 43.

### Relation to the others

This is the **foundational** paper of the fractional strand, and it is the origin of what the literature calls the "fractional energy dissipation law of Tang-Yu-Zhou". Its kernel-positivity lemma and its discrete counterpart (positive definiteness of the L1 quadratic form $B$) are the direct ancestors of the variable-step positivity theory of [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|paper 74]]. Papers 43 and 57 attack exactly the two limitations the paper concedes about itself: paper 43 supplies a second-order maximum-principle-preserving scheme on nonuniform meshes, and paper 57 supplies a variable-step scheme preserving energy stability and the maximum bound principle simultaneously, upgrading the integral law here into a differential one that degenerates to the classical law as $\alpha\to1$. The integer-order branch (papers 48, 52, 58, 67, 69, 91, 104) carries the same "convolution positivity" philosophy over to variable-step BDF and implicit-explicit methods.

## 43: a second-order maximum-principle-preserving scheme on nonuniform steps

### The idea

The scheme of paper 40 is first order in time on a uniform mesh, and for fractional problems that has two concrete costs. First, the solution of a time-fractional Allen-Cahn equation carries an intrinsic initial singularity $u_t\sim\mathcal O(t^{\alpha-1})$, so a uniform mesh destroys the formal order and refinement near $t=0$ (a graded mesh) is mandatory. Second, the dynamics span several time scales — fast initial evolution followed by extremely slow coarsening — so long-time simulation has to grow the step.

But nonuniform meshes break the classical positivity and monotonicity arguments for fractional convolution kernels, and no second-order maximum-principle-preserving scheme was available. The mechanism here has two layers.

The outer layer is **choosing a formula whose kernel structure survives on a nonuniform mesh**: the Alikhanov (L2-$1_\sigma$) formula, which collocates the equation at the off-grid point $t_{n-\theta}$ (with $\theta=\alpha/2$), so that the kernels $A^{(n)}_{n-k}$ remain positive and remain monotonically decreasing in the lag index — except that they now carry **two** indices, the superscript recording the current level. Two indices means the kernel is no longer a Toeplitz convolution, and that is exactly what makes variable steps hard.

The inner layer is the **discrete complementary convolution (DCC) kernels**: a family $P^{(n)}_{n-j}\ge0$ constructed so that

$$
\sum_{j=k}^{n}P^{(n)}_{n-j}A_{j-k}^{(j)}\equiv1 .
$$

They act as an inverse operator: convolve an inequality of the form $\sum_kA^{(n)}_{n-k}\nabla_\tau v^k\le\cdots$ once on both sides and the left collapses to $v^n-v^0$. **The convergence estimate obtained this way is sharp because the local truncation error itself has convolution structure**: convolving a sequence of truncation errors with $P$ yields something far smaller than the maximum of that sequence. This is most visible at the first level — $|\Upsilon^{1}|$ is only $\mathcal O(1)$ when $\sigma=\alpha$, whereas $P_0^{(1)}|\Upsilon^1|$ is $\mathcal O(\tau_1^{\sigma})$.

### Setting

The two-dimensional time-fractional Allen-Cahn equation on $\Omega=(0,L)^2$ with periodic boundary conditions:

$$
\partial_{t}^{\alpha}u=\varepsilon^{2}\Delta u-f(u),\qquad f(u)=u^{3}-u,
\qquad u(\mathbf x,0)=u_{0}(\mathbf x),
$$

with the Caputo derivative written in convolution form $(\partial_{t}^{\alpha}v)(t)=(\mathcal I_{t}^{1-\alpha}v')(t)$, where $\omega_{\mu}(t):=t^{\mu-1}/\Gamma(\mu)$. This $\omega_\mu$ notation runs through the whole Liao-Tang-Zhou series. Inherited from paper 40 at the continuous level are the energy law $E(t)\le E(0)$ and the maximum principle $|u|\le1$.

> [!warning] The step-ratio convention here is the reciprocal of the BDF papers'
> This paper uses
>
> $$
> \rho_k:=\frac{\tau_k}{\tau_{k+1}},\qquad \rho:=\max_{k\ge1}\rho_k,
> $$
>
> whereas [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|papers 48 and 52]] use $r_k:=\tau_k/\tau_{k-1}$, so $\rho_k=1/r_{k+1}$. **Mixing the two conventions when quoting thresholds across papers is the easiest mistake to make in this literature.**

Two mesh assumptions:

- **M1**: the maximum step ratio is $\rho=7/4$. The authors read this as permitting "a sequence of decreasing steps with reduction factor as low as $4/7$", while placing **no restriction whatsoever** on step **increase**.
- **M2** (used only for convergence): there exist mesh-independent $C_{1\gamma},C_{2\gamma}>0$ with $\tau_k\le\tau\min\{1,C_{1\gamma}t_k^{1-1/\gamma}\}$ and $t_k\le C_{2\gamma}t_{k-1}$. Quasi-uniform meshes correspond to $\gamma=1$; the model case is the graded mesh $t_k=T(k/N)^{\gamma}$.

### Derivation

**Step one: the Alikhanov formula and its two-index kernels.** Writing $\Pi_{1,k}$ and $\Pi_{2,k}$ for linear and quadratic interpolation on $\{t_{k-1},t_k\}$ and $\{t_{k-1},t_k,t_{k+1}\}$,

$$
(\partial_{\tau}^{\alpha}v)^{n-\theta}
:=\int_{t_{n-1}}^{t_{n-\theta}}\omega_{1-\alpha}(t_{n-\theta}-s)(\Pi_{1,n}v)'(s)\,\mathrm ds
+\sum_{k=1}^{n-1}\int_{t_{k-1}}^{t_k}\omega_{1-\alpha}(t_{n-\theta}-s)(\Pi_{2,k}v)'(s)\,\mathrm ds,
$$

which expands via the two families of coefficients

$$
a_{n-k}^{(n)}:=\frac{1}{\tau_{k}}\int_{t_{k-1}}^{\min\{t_{k},t_{n-\theta}\}}
\omega_{1-\alpha}(t_{n-\theta}-s)\,\mathrm ds,
\qquad
b_{n-k}^{(n)}:=\frac{2}{\tau_{k}(\tau_{k}+\tau_{k+1})}
\int_{t_{k-1}}^{t_{k}}\bigl(s-t_{k-\frac12}\bigr)\omega_{1-\alpha}(t_{n-\theta}-s)\,\mathrm ds
$$

into the compact convolution form $(\partial_{\tau}^{\alpha}v)^{n-\theta}=\sum_{k=1}^{n}A_{n-k}^{(n)}\nabla_{\tau}v^{k}$ with

$$
A_{n-k}^{(n)}:=
\begin{cases}
a_{0}^{(n)}+\rho_{n-1}b_{1}^{(n)}, & k=n,\\
a_{n-k}^{(n)}+\rho_{k-1}b_{n-k+1}^{(n)}-b_{n-k}^{(n)}, & 2\le k\le n-1,\\
a_{n-1}^{(n)}-b_{n-1}^{(n)}, & k=1 .
\end{cases}
$$

The fully discrete scheme (central differences $D_h$ in space, periodic boundary) is

$$
(\partial_{\tau}^{\alpha}u)^{n-\theta}=\varepsilon^{2}D_{h}u^{n-\theta}-f(u)^{n-\theta},
\qquad
f(u)^{n-\theta}:=\theta f(u^{n-1})+(1-\theta)f(u^{n}) .
$$

**Step two: three kernel estimates.** Under M1,

1. $A^{(n)}_{0}\le\dfrac{24}{11\tau_{n}}\displaystyle\int_{t_{n-1}}^{t_{n}}\omega_{1-\alpha}(t_{n}-s)\,\mathrm ds$, and the kernels admit a matching lower bound with factor $\dfrac{4}{11}$;
2. monotonicity: $A^{(n)}_{n-k-1}-A^{(n)}_{n-k}>0$ for $1\le k\le n-1$;
3. the first kernel dominates the second: $\dfrac{1-2\theta}{1-\theta}A^{(n)}_{0}-A^{(n)}_{1}>0$ for $n\ge2$, the prefactor being $\dfrac{2-2\alpha}{2-\alpha}$ at $\theta=\alpha/2$.

These are stronger than Alikhanov's original uniform-mesh estimates; they imply the two hypotheses of the fractional Grönwall lemma, with constant

$$
\pi_A=\frac{11}{4}.
$$

(The lower bound in item 1 is printed in the arXiv source as $\frac{4}{11\tau_n}\int_{t_{n-1}}^{t_n}$, which does not match the kernel index $n-k$; the hypotheses of the Grönwall lemma require $\frac{4}{11\tau_k}\int_{t_{k-1}}^{t_k}\omega_{1-\alpha}(t_n-s)\mathrm ds$. This site judges it an index misprint; the constant $4/11$ itself is consistent with $\pi_A=11/4$.)

**Step three: DCC kernels and the complementarity identity.**

$$
P_{0}^{(n)}:=\frac{1}{A_{0}^{(n)}},
\qquad
P_{n-j}^{(n)}:=\frac{1}{A_{0}^{(j)}}\sum_{k=j+1}^{n}
\bigl(A_{k-j-1}^{(k)}-A_{k-j}^{(k)}\bigr)P_{n-k}^{(n)},\quad 1\le j\le n-1 .
$$

(The source writes the prefactor as $1/p_0^{(j)}$; it should be $1/A_0^{(j)}$.) Kernel monotonicity gives $P^{(n)}_{n-j}\ge0$, and

$$
\sum_{j=k}^{n}P^{(n)}_{n-j}A_{j-k}^{(j)}\equiv1,
\qquad
\sum_{j=1}^{n}P^{(n)}_{n-j}\,\omega_{1+m\alpha-\alpha}(t_{j})
\le\pi_{A}\,\omega_{1+m\alpha}(t_{n}),\quad m=0,1 .
$$

The second is a quantitative estimate for $P$ acting on power-type functions, and it is the fuel of the Grönwall lemma.

**Step four: the discrete fractional Grönwall inequality.** With $E_\alpha$ the Mittag-Leffler function and $\lambda=\lambda_0+\lambda_1$, if the maximum step satisfies

$$
\tau\le\frac{1}{\sqrt[\alpha]{2\Gamma(2-\alpha)\lambda\pi_{A}}},
$$

and a nonnegative sequence satisfies $\sum_{k=1}^{n}A_{n-k}^{(n)}\nabla_{\tau}v^{k}\le\lambda_{0}v^{n}+\lambda_{1}v^{n-1}+\xi^{n}+\eta^{n}$, then

$$
v^{n}\le2E_\alpha\bigl(2\max\{1,\rho\}\lambda\pi_{A}t_{n}^{\alpha}\bigr)
\Bigl(v^{0}+\Gamma(1-\alpha)\pi_{A}\max_{k\le n}\{t_{k}^{\alpha}\xi^{k}\}
+\pi_{A}\,\omega_{1+\alpha}(t_{n})\max_{k\le n}\eta^{k}\Bigr).
$$

**Step five: the maximum principle is an $\ell^\infty$ matrix argument plus induction.** Three ingredients: (a) $D_h$ is symmetric, negative semi-definite, with $d_{ii}=-\max_i\sum_{j\ne i}|d_{ij}|$; (b) a cubic lemma — if $B$ satisfies the same diagonal condition and $A=aI-B$ with $a>0$, then for $c>0$,

$$
\|AV\|_{\infty}\ge a\|V\|_{\infty},
\qquad
\|AV+cV^{3}\|_{\infty}\ge a\|V\|_{\infty}+c\|V\|_{\infty}^{3};
$$

(c) rearranging the scheme as

$$
(A^{(n)}_{0}-1+\theta)u^{n}-(1-\theta)\varepsilon^{2}D_{h}u^{n}+(1-\theta)(u^n)^{.3}
=\mathcal{L}^{n-2}(u)+\cdots,
\qquad
\mathcal{L}^{n-2}(u):=\sum_{k=1}^{n-2}\bigl(A_{n-k-1}^{(n)}-A_{n-k}^{(n)}\bigr)u^{k}+A_{n-1}^{(n)}u^{0},
$$

whose history coefficients are all nonnegative by monotonicity, so the induction hypothesis $\|u^k\|_\infty\le1$ for $k\le n-1$ controls the right-hand side while the cubic lemma bounds the left from below.

### Theorems

**(Unique solvability.)** Under M1 and $\tau\le\sqrt[\alpha]{\omega_{2-\alpha}(1-\theta)/(1-\theta)}$, the nonlinear scheme is uniquely solvable. The proof uses $A_0^{(n)}\ge a_0^{(n)}=\omega_{2-\alpha}(1-\theta)/\tau_n^{\alpha}\ge1-\theta$ to make $G_h:=A_0^{(n)}-1+\theta-(1-\theta)\varepsilon^{2}D_h$ positive definite, so that $u^n$ is the unique minimiser of the strictly convex functional $\frac12w^{T}G_hw+\frac{1-\theta}{4}\sum_kw_k^{4}-w^{T}g(u^{n-1})$.

**(Discrete maximum principle.)** Under M1 (with $\rho=7/4$) and

$$
\tau\le\min\Bigl\{
\sqrt[\alpha]{\frac{\theta\,\omega_{2-\alpha}(1-\theta)}{2(1-\theta)}},\
\sqrt[\alpha]{\frac{h^{2}\,\omega_{2-\alpha}(1-\theta)}{4\varepsilon^{2}}}\Bigr\},
$$

$\|u^{0}\|_{\infty}\le1$ implies $\|u^{k}\|_{\infty}\le1$ for $1\le k\le N$. This is the headline result: "the first second-order maximum-principle-preserving scheme for the time-fractional Allen-Cahn equation". **The second constraint couples $\tau$ to the spatial mesh $h$ and to $\varepsilon$; it is a real, if mild, restriction rather than an unconditional result.**

**(Sharp maximum-norm convergence.)** Assume $\|u^0\|_{L^\infty}\le1$, the regularity hypotheses

$$
\|u(t)\|_{W^{4,\infty}(\Omega)}\le C_u,
\qquad
\|u^{(\ell)}(t)\|_{W^{2,\infty}(\Omega)}\le C_u(1+t^{\sigma-\ell})\ (\ell=1,2,3),\ \sigma\in(0,1),
$$

M1, and the two step conditions above together with $\tau\le\sqrt[\alpha]{\omega_{2-\alpha}(1)/11}$. Then

$$
\|u(t_{n})-u^{n}\|_{\infty}\le C_{u}\Bigl(\frac{\tau_{1}^{\sigma}}{\sigma}
+\frac{1}{1-\alpha}\max_{2\le k\le n}t_{k}^{\alpha}t_{k-1}^{\sigma-3}\tau_{k}^{3-\alpha}
+h^{2}\Bigr);
$$

adding M2 gives

$$
\|u(t_{n})-u^{n}\|_{\infty}\le\frac{C_{u}}{\sigma(1-\alpha)}\tau^{\min\{\gamma\sigma,\,2\}}+C_{u}h^{2},
$$

so **the optimal $\mathcal O(\tau^{2})$ rate is attained once the grading parameter satisfies $\gamma\ge\max\{1,2/\sigma\}$**. The proof uses the maximum principle to avoid assuming that $f$ is globally Lipschitz.

**(Convolution structure of the error.)** This is the technical source of the sharpness. Under M1, for $v\in C^3((0,T])$ with $\int_0^Ts^2|v'''(s)|\mathrm ds<\infty$, the local consistency error $\Upsilon^n[v]:=(\partial_t^\alpha v)(t_{n-\theta})-(\partial_\tau^\alpha v)^{n-\theta}$ itself has convolution structure,

$$
|\Upsilon^{n}[v]|\le A_{0}^{(n)}G_{\mathrm{loc}}^{n}
+\sum_{k=1}^{n-1}\bigl(A_{n-k-1}^{(n)}-A_{n-k}^{(n)}\bigr)G_{\mathrm{his}}^{k},
$$

$$
G_{\mathrm{loc}}^{k}:=\frac{3}{2}\int_{t_{k-1}}^{t_{k-1/2}}(s-t_{k-1})^{2}|v'''(s)|\mathrm ds
+\frac{3\tau_{k}}{2}\int_{t_{k-1/2}}^{t_{k}}(t_{k}-s)|v'''(s)|\mathrm ds,
$$

$$
G_{\mathrm{his}}^{k}:=\frac{5}{2}\int_{t_{k-1}}^{t_{k}}(s-t_{k-1})^{2}|v'''(s)|\mathrm ds
+\frac{5}{2}\int_{t_{k}}^{t_{k+1}}(t_{k+1}-s)^{2}|v'''(s)|\mathrm ds .
$$

(The limits in the second term of $G_{\mathrm{loc}}$ are transposed in the source; they are written here in their evident intended order.) Because $P$ and $A$ are complementary, the global error is $\sum_jP_{n-j}^{(n)}|\Upsilon^j[v]|$, which is far smaller than $\max_j|\Upsilon^j[v]|$: at $n=1$ with $\sigma=\alpha$, $|\Upsilon^{1}|\le C_u\tau_1^{\sigma-\alpha}/\sigma$ is merely $\mathcal O(1)$, whereas $P_0^{(1)}|\Upsilon^1|\le G_{\mathrm{loc}}^1\le C_u\tau_1^{\sigma}/\sigma$. **This is superconvergence borrowed from the complementary kernels, and it is why the estimate reflects temporal regularity rather than the worst-case truncation error.**

### Numerical experiments

The implementation uses a fast Alikhanov formula built on a sum-of-exponentials approximation of the kernel (absolute SOE tolerance $\epsilon=10^{-12}$), a nonlinear iteration at each level with termination error $\eta=10^{-12}$, and the adaptive step selector

$$
\tau_{\mathrm{ada}}(e,\tau)=S_a\Bigl(\frac{tol}{e}\Bigr)^{1/2}\tau,
$$

where $e$ is the relative difference between a first-order (backward-Euler/L1) solution and the second-order solution of this paper, plus a guard that **never lets the step shrink below $\frac23\tau_{n-1}$**. That guard is matched to M1: $\rho\le7/4$ permits reduction by at most $4/7\approx0.571$, and $2/3$ sits inside that.

**Test one (accuracy, $\alpha=0.8$).** Prescribed temporal regularity $\sigma=0.8$ and $\sigma=0.4$, with theoretically optimal grading $\gamma_{\mathrm{opt}}=\max\{1,2/\sigma\}$:

| Temporal regularity $\sigma$ | Optimal grading $\gamma_{\mathrm{opt}}$ | Observed order, quasi-uniform ($\gamma=1$) | Observed order, graded |
| ---------------------------- | --------------------------------------- | ------------------------------------------ | ---------------------- |
| $0.8$                        | $2.5$                                   | $\approx0.80$                              | $\approx2.1$–$2.34$    |
| $0.4$                        | $5$                                     | not verified here                          | $\approx2.1$–$2.34$    |

The $0.80$ observed at $\sigma=0.8$, $\gamma=1$ matches the theoretical prediction $\min\{\gamma\sigma,2\}=0.8$ exactly — **this is not just "the scheme converges" but a verification that the estimate is sharp**: the theorem does not only give an upper bound, it predicts how far the order falls when regularity is insufficient. Switching to the graded mesh with $\gamma=\gamma_{\mathrm{opt}}$ raises the observed order above $2$, consistent with $\min\{\gamma\sigma,2\}=2$.

**Test two (merging drops, long time).** $\alpha=0.4,0.7,0.9$, run to $T=100$, with snapshots at $t=1,10,50,100$, tracking both the maximum norm (which stays $\le1$, confirming the discrete maximum principle) and the discrete energy (which decays).

**Test three (what happens when the conditions are ignored).** A comparison at $\varepsilon=0.02$ and $\varepsilon=0.08$ across different time meshes, showing that **the maximum bound really can be violated when the step-size condition is disregarded**. A negative test of this kind is worth more than a positive one: it establishes that the condition $\tau\le\sqrt[\alpha]{h^{2}\omega_{2-\alpha}(1-\theta)/(4\varepsilon^{2})}$ in the theorem is not an artefact of the proof technique but carries real content.

A gap between theory and experiment remains: what the paper proves is the maximum principle and maximum-norm convergence, **not** a discrete energy law; the energy decay in test two is an observation only. Paper 57 closes that gap.

### Relation to the others

It answers directly the "higher-order schemes" open problem stated in the concluding remarks of paper 40, and inherits from paper 40 the continuous energy law and maximum bound principle. It introduces into this series the **DCC (discrete complementary convolution) kernels**, their complementarity identity, and the fractional Grönwall machinery — tools re-used repeatedly afterwards and **dualised** in [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|papers 48, 52, 58 and 67]] into **DOC (discrete orthogonal convolution) kernels**: DCC sits on one side of the convolution, gives a complementarity identity with constant $1$, and yields $\ell^\infty$ and Grönwall-type estimates; DOC sits on the other side, gives an orthogonality identity with a Kronecker $\delta$, and yields energy and $L^2$ estimates. Paper 57 is the direct sequel on the fractional side, while paper 74 takes the question implicit here — when are variable-step L1-type kernels positive definite — and settles it as an object in its own right.

## 57: change the energy, get a differential law

### The idea

The difficulty in paper 40 can be described this way: **the fractional operator is on the wrong side of the equation.** As long as the time derivative is $\partial_t^\alpha$, testing with $u_t$ produces a nonlocal quadratic form that acquires a sign only after integration.

The move here is to put it on the other side. Apply ${}^{R}\!\partial_t^{1-\alpha}$ to both sides of the fractional gradient flow $\partial_t^{\alpha}u=-\delta E/\delta u$; by the semigroup property the left-hand side collapses to $u_t$, giving

$$
\partial_t u=-{}^{R}\!\partial_t^{1-\alpha}\Bigl(\frac{\delta E}{\delta u}\Bigr).
$$

**The time derivative is now local, so the standard energy test function applies as usual, and the nonlocality has landed on the variational derivative — where a usable positivity inequality happens to exist.** That inequality (of Riemann-Liouville type) is slightly stronger than positivity: it bounds $v\cdot{}^{R}\!\partial_t^{1-\alpha}v$ below by "a total derivative" plus "a nonnegative term". The total-derivative part can be absorbed into the definition of the energy, and the nonnegative part is the dissipation. So **the new energy $\mathcal E_\alpha$ is not invented out of thin air; it is that total derivative.**

The construction has an elegant by-product: the new law degenerates exactly to the classical one as $\alpha\to1$. That is, $\mathcal E_\alpha$ is not an analogue of the classical energy but a generalisation of it.

### Setting

On $\Omega=(0,L)^2$ with periodic boundary conditions,

$$
\partial_t^{\alpha}u=\varepsilon^{2}\Delta u-f(u),
\qquad F(u)=\tfrac14(1-u^2)^2,\ f=F',
$$

$$
\partial_t^{\alpha}v:=\mathcal I_t^{1-\alpha}v',
\qquad
(\mathcal I_t^{\mu}v)(t):=\int_0^t\omega_{\mu}(t-s)v(s)\,\mathrm ds,
\qquad
\omega_{\mu}(t):=\frac{t^{\mu-1}}{\Gamma(\mu)} .
$$

The equation is read as a **fractional gradient flow** $\partial_t^{\alpha}u=-\delta E/\delta u$ with $E[u]=\int_\Omega(\frac{\varepsilon^2}{2}|\nabla u|^2+F(u))\mathrm dx$.

### Derivation

**Step one: the rewriting.** With the Riemann-Liouville derivative ${}^{R}\!\partial_t^{\alpha}v:=\partial_t\mathcal I_t^{1-\alpha}v$ and the semigroup identity

$$
{}^{R}\!\partial_t^{1-\alpha}\bigl(\partial_t^{\alpha}v\bigr)
=\partial_t\mathcal I_t^{1}v'=v',
$$

the equation is equivalently rewritten as $\partial_t u=-{}^{R}\!\partial_t^{1-\alpha}(\delta E/\delta u)$.

**Step two: the variational energy and its differential law.** Using the Riemann-Liouville inequality (attributed to Alsaedi-Ahmad-Kirane)

$$
v(t)\bigl({}^{R}\!\partial_t^{1-\alpha}v\bigr)(t)
\ \ge\ \tfrac12\bigl({}^{R}\!\partial_t^{1-\alpha}v^{2}\bigr)(t)
+\tfrac12\omega_{\alpha}(t)v^{2}(t),
\qquad \forall v\in C[0,T],
$$

together with $\frac{\mathrm dE}{\mathrm dt}=-\bigl(\frac{\delta E}{\delta u},{}^{R}\!\partial_t^{1-\alpha}\frac{\delta E}{\delta u}\bigr)$, define

$$
\mathcal{E}_{\alpha}[u]:=E[u]+\frac12\,\mathcal{I}_t^{\alpha}
\Bigl\|\frac{\delta E}{\delta u}\Bigr\|^{2},
$$

whereupon

$$
\frac{\mathrm d\mathcal{E}_{\alpha}}{\mathrm dt}
+\frac12\,\omega_{\alpha}(t)\Bigl\|\frac{\delta E}{\delta u}\Bigr\|^{2}\le0,
\qquad \forall t>0 .
$$

Both factors of $\frac12$ are as printed in the paper. As $\alpha\to1$, $\mathcal I_t^\alpha\to\mathcal I_t^1$ and $\omega_\alpha(t)\to1$, recovering $\frac{\mathrm dE}{\mathrm dt}+\|\frac{\delta E}{\delta u}\|^{2}\le0$. The paper calls this property **asymptotically energy-dissipation preserving**.

**Step three: the L1$_R$ formula.** On variable steps the Riemann-Liouville derivative is discretised as

$$
\bigl({}^{R}\!\partial_{\tau}^{1-\alpha}v\bigr)^{n-\frac12}
:=\frac{1}{\tau_{n}}\int_{t_{n-1}}^{t_{n}}\frac{\partial}{\partial t}
\int_{0}^{t}\omega_{\alpha}(t-s)(\Pi_{0}v)(s)\,\mathrm ds\,\mathrm dt
\ \triangleq\ \frac{1}{\tau_{n}}\sum_{k=1}^{n}a_{n-k}^{(n)}v^{k-\frac12},
$$

with $\Pi_0v$ the piecewise constant interpolant equal to $v^{k-\frac12}$ on $(t_{k-1},t_k]$. The auxiliary (DCO) sequence and the kernels are

$$
q_{n-k}^{(n)}:=\int_{t_{k-1}}^{t_{k}}\omega_{\alpha}(t_{n}-s)\,\mathrm ds
=\sum_{j=k}^{n}a_{j-k}^{(j)}>0,
$$

$$
a_{0}^{(n)}:=q_{0}^{(n)}>0\ (n\ge1),
\qquad
a_{n-k}^{(n)}:=q_{n-k}^{(n)}-q_{n-k-1}^{(n-1)}<0\ (n\ge k+1\ge2).
$$

The sign pattern is worth noting: the first kernel is positive and all the rest are negative. Explicitly $a^{(n)}_0=\omega_{1+\alpha}(\tau_n)=\tau_n^{\alpha}/\Gamma(1+\alpha)$, while $a^{(j)}_{j-k}=\int_{t_{j-1}}^{t_j}\!\int_{t_{k-1}}^{t_k}\omega_{\alpha-1}(t-s)\,\mathrm ds\,\mathrm dt<0$. The formula is due to Mustapha (in the linear subdiffusion setting, with order $1+\alpha$); the authors name it L1$_R$ to distinguish it from the L1 formula for the Caputo derivative.

**Step four: kernel positivity — proved at the discrete level, not inherited from the continuous kernel.** For any real sequence $\{w_k\}$,

$$
2w_{k}\sum_{j=1}^{k}a_{k-j}^{(k)}w_{j}
\ \ge\ w_{k}^{2}\sum_{j=1}^{k}a_{k-j}^{(k)}
+\sum_{j=1}^{k}q_{k-j}^{(k)}w_{j}^{2}-\sum_{j=1}^{k-1}q_{k-j-1}^{(k-1)}w_{j}^{2},
$$

and after summation

$$
2\sum_{k=1}^{n}w_{k}\sum_{j=1}^{k}a_{k-j}^{(k)}w_{j}
\ \ge\ \sum_{k=1}^{n}\Bigl(q_{n-k}^{(n)}+\sum_{j=1}^{k}a_{k-j}^{(k)}\Bigr)w_{k}^{2}
\ >\ 0,
\qquad n\ge1,\ w\not\equiv0 .
$$

Two identities carry the proof. Complete monotonicity of $\omega_\alpha$ gives

$$
q_{k-j-1}^{(k-1)}-q_{k-j}^{(k)}
=\int_{t_{j-1}}^{t_j}\bigl[\omega_{\alpha}(t_{k-1}-s)-\omega_{\alpha}(t_k-s)\bigr]\mathrm ds>0,
$$

and

$$
\sum_{j=1}^{k}a_{k-j}^{(k)}
=\sum_{j=1}^{k}q_{k-j}^{(k)}-\sum_{j=1}^{k-1}q_{k-j-1}^{(k-1)}
=\int_{t_{k-1}}^{t_{k}}\omega_{\alpha}(s)\,\mathrm ds>0 .
$$

**This positivity requires no step-ratio restriction whatsoever.** The contrast with [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|variable-step BDF]] is sharp: there, positive definiteness of the quadratic form produces thresholds such as $r_k<(3+\sqrt{17})/2$, whereas here the fractional kernels are positive definite directly by complete monotonicity. A step-size bound appears only in the maximum-bound argument and has nothing to do with the energy law.

**Step five: the scheme and two constructive designs.** Setting $v:=-\delta E/\delta u$ splits the equation into $\partial_t u={}^{R}\!\partial_t^{1-\alpha}v$ and $v=\varepsilon^{2}\Delta u-f(u)$, discretised in Crank-Nicolson form as

$$
\partial_{\tau}u^{n-\frac12}=\bigl({}^{R}\!\partial_{\tau}^{1-\alpha}v\bigr)^{n-\frac12},
\qquad
v^{n-\frac12}=\varepsilon^{2}D_{h}u^{n-\frac12}-H(u^{n},u^{n-1}),
$$

with the second-order approximation of the nonlinear term taken to be

$$
H(u^{n},u^{n-1}):=\tfrac13(u^{n})^{.3}+\tfrac12(u^{n-1})^{.2}\!\circ u^{n}
+\tfrac16(u^{n-1})^{.3}-\tfrac12\bigl(u^{n}+u^{n-1}\bigr)
$$

($\circ$ and the powers being element-wise). This particular $H$ is **engineered** so that $H(a,b)(a-b)\ge F(a)-F(b)$ holds pointwise — the discrete counterpart of the chain rule, and the reason the energy argument can proceed. The discrete variational energy is

$$
\mathcal{E}_{\alpha}[u^{n}]:=E[u^{n}]+\frac12h^{2}\sum_{i,j=1}^{M_1}\sum_{k=1}^{n}
q_{n-k}^{(n)}\bigl(v_{ij}^{k-\frac12}\bigr)^{2},
\qquad
E[u^{n}]:=h^{2}\sum_{i,j}F(u_{ij}^{n})-\tfrac12\varepsilon^{2}h^{2}(u^{n})^{T}D_{h}u^{n},
$$

where the $q$ kernels themselves constitute a numerical fractional integral, $(\mathcal I_\tau^{\alpha}v)^n=\sum_{k}q_{n-k}^{(n)}v^{k-\frac12}$, with $({}^{R}\!\partial_{\tau}^{1-\alpha}v)^{n-\frac12}=\partial_\tau(\mathcal I_\tau^{\alpha}v)^{n-\frac12}$.

**Step six: DOC kernels and a reversible transformation between the two fractional derivatives.** Define the discrete orthogonal convolution kernels $\theta$ of the L1$_R$ kernels:

$$
\theta_{0}^{(n)}:=\frac{1}{a^{(n)}_{0}},
\qquad
\theta_{n-k}^{(n)}:=-\frac{1}{a^{(k)}_{0}}\sum_{j=k+1}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}
\ \ (1\le k\le n-1).
$$

Orthogonality holds in **both** directions — this is what reversibility means — and is complementary to the DCO kernels $q$:

$$
\sum_{j=k}^{n}a^{(n)}_{n-j}\theta_{j-k}^{(j)}\equiv\delta_{nk},
\qquad
\sum_{j=k}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}\equiv\delta_{nk},
\qquad
\sum_{j=k}^{n}q_{n-j}^{(n)}\theta_{j-k}^{(j)}\equiv1 .
$$

The paper proves $\theta^{(n)}_0=\Gamma(1+\alpha)\tau_n^{-\alpha}$, that every $\theta^{(n)}_j>0$, the sharp leading gap

$$
\theta_{0}^{(n)}-\theta_{1}^{(n)}
=\frac{\omega_{1+\alpha}(r_{n}+1)-\omega_{1+\alpha}(r_{n})}
{\omega_{1+\alpha}(\tau_{n})\,\omega_{1+\alpha}(1)}
\ >\ \frac{\omega_{\alpha}(r_{n}+1)}{\omega_{1+\alpha}(\tau_{n})\,\omega_{1+\alpha}(1)},
$$

and **monotone decrease** $\theta_{0}^{(n)}>\theta_{1}^{(n)}>\cdots>\theta_{n-1}^{(n)}>0$ for $n\ge2$. The monotonicity proof introduces the auxiliary kernels $\zeta^{(n)}_0:=\theta^{(n)}_0$ and $\zeta^{(n)}_{n-j}:=\theta^{(n)}_{n-j}-\theta^{(n)}_{n-j-1}$ (which satisfy $\sum_{j=k}^{n}\zeta^{(n)}_{n-j}q^{(j)}_{j-k}=\delta_{nk}$) and **invokes the algebraic criteria of [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|paper 74]] directly** — $q^{(n)}_j>0$, $q^{(n-1)}_{j-1}>q^{(n)}_j$, and $q^{(n-1)}_{j-1}q^{(n)}_{j+1}>q^{(n-1)}_{j}q^{(n)}_{j}$ — so paper 74 is logically prior to this one even though it was published later.

Convolving the first scheme equation with $\theta^{(n)}_{n-j}$ and using orthogonality,

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

In other words, the DOC kernels define a **new** discrete Caputo derivative whose kernels are positive and monotonically decreasing on nonuniform meshes, exactly as the classical L1 kernels are, with the transformation made reversible by two-sided orthogonality. The paper cautions that this is an **indirect** approximation whose accuracy differs from that of the direct L1 formula (error order $2-\alpha$).

The maximum-bound argument runs on this Caputo-equivalent form, by induction, using the symmetry and negative semi-definiteness of $D_h$, an $\ell^\infty$ lemma

$$
\bigl\|(aI-B)V+U^{.2}\!\circ V+cV^{.3}\bigr\|_{\infty}
\ \ge\ a\|V\|_{\infty}+\|U\|_{\infty}^{2}\|V\|_{\infty}+c\|V\|_{\infty}^{3},
$$

and the rearranged scheme

$$
\Bigl(\theta^{(n)}_{0}-\tfrac12+\tfrac12(u^{n-1})^{.2}-\tfrac{\varepsilon^{2}}{2}D_{h}\Bigr)u^{n}
+\tfrac13(u^n)^{.3}=\cdots+\mathcal{L}^{n-2}(u),
\qquad
\mathcal{L}^{n-2}(u):=\sum_{k=1}^{n-2}\bigl(\theta_{n-k-1}^{(n)}-\theta_{n-k}^{(n)}\bigr)u^{k}
+\theta_{n-1}^{(n)}u^{0},
$$

whose history coefficients are nonnegative **precisely because the DOC kernels are positive and decreasing**.

### Theorems

**(Discrete variational energy dissipation law, unconditional.)** The Crank-Nicolson scheme satisfies, **unconditionally** — no restriction on step size, none on step ratio —

$$
\partial_{\tau}\bigl(\mathcal{E}_{\alpha}[u]\bigr)^{n-\frac12}
+\frac{1}{2\tau_{n}}\int_{t_{n-1}}^{t_{n}}\omega_{\alpha}(s)\,\mathrm ds
\sum_{i,j=1}^{M_1}h^{2}\bigl(v_{ij}^{n-\frac12}\bigr)^{2}\ \le\ 0,
\qquad n\ge1 .
$$

Two ingredients: the chain-rule property of $H$, and the level-by-level form of the L1$_R$ positivity lemma, $v^{n-\frac12}\sum_k a^{(n)}_{n-k}v^{k-\frac12}\ge\frac12(\mathcal I_\tau^{\alpha}v^2)^n-\frac12(\mathcal I_\tau^{\alpha}v^2)^{n-1}+\frac12\int_{t_{n-1}}^{t_n}\omega_\alpha(s)\mathrm ds\,(v^{n-\frac12})^2$.

**(Asymptotic preservation as $\alpha\to1$.)** As $\alpha\to1$ one has $q^{(n)}_{n-k}\to\tau_k$, hence $(\mathcal I_\tau^\alpha v)^n\to\sum_k\tau_kv^{k-\frac12}$ and $({}^{R}\!\partial_\tau^{1-\alpha}v)^{n-\frac12}\to v^{n-\frac12}$, and the discrete law becomes

$$
\partial_{\tau}\bigl(E[u]\bigr)^{n-\frac12}
+\sum_{i,j}h^{2}\bigl(v_{ij}^{n-\frac12}\bigr)^{2}\ \le\ 0,
$$

the standard discrete energy dissipation law of the classical Allen-Cahn equation.

**(Unique solvability.)** If $\tau<\sqrt[\alpha]{2\Gamma(1+\alpha)}$, then $\theta^{(n)}_0=\Gamma(1+\alpha)\tau_n^{-\alpha}>\frac12$, the matrix $G_h:=(\theta^{(n)}_0-\frac12+\frac12(u^{n-1})^{.2})I-\frac{\varepsilon^2}{2}D_h$ is positive definite, the objective functional $\frac12w^TG_hw+\frac{1}{12}\sum_kw_k^4-w^T\mathrm G_0(u^{n-1})$ is strictly convex, and the scheme is uniquely solvable.

**(Discrete maximum bound principle.)** If the steps satisfy

$$
\tau_{n}\le\sqrt[\alpha]{\min\Bigl\{\frac12,\ \frac{h^{2}}{2\varepsilon^{2}}\Bigr\}
\cdot\frac{\alpha\,\Gamma(1+\alpha)}{(1+r_{n})^{1-\alpha}}},
$$

then $\|u^{0}\|_{\infty}\le1$ implies $\|u^{k}\|_{\infty}\le1$. **Note that $r_n=\tau_n/\tau_{n-1}$ enters only through the mild factor $(1+r_n)^{\alpha-1}<1$; there is no upper bound on $r_n$ at all.** That is what lets the paper claim "the first work on variable time-stepping schemes that can preserve both the energy stability and the maximum bound principle".

### Numerical experiments

The implementation again uses a fast L1$_R$ algorithm based on a sum-of-exponentials approximation (absolute tolerance $\epsilon=10^{-12}$, cut-off time $\Delta t=10^{-12}$).

| Test | Setup                                                                                     | What is checked                |
| ---- | ----------------------------------------------------------------------------------------- | ------------------------------ |
| 1    | forced model $\partial_tu=-{}^{R}\!\partial_t^{1-\alpha}(\delta E/\delta u)+g$ on $(0,1)^2\times(0,1]$, $\varepsilon=0.1$, manufactured solution $u=\omega_{1+\sigma}(t)\sin(2\pi x)\sin(2\pi y)$, e.g. $\alpha=0.6,\sigma=0.4$ | temporal order (expected $1+\alpha$) |
| 2    | maximum bound principle: $\alpha=0.7,\,0.9$ at $\tau=0.1,\,0.8,\,1.0$                       | whether the discrete maximum norm stays $\le1$ |
| 3    | coarsening: $(0,2\pi)^2$, $\varepsilon=0.05$, $128\times128$ spatial mesh, random data uniform on $[-0.001,0.001]$ | initial layer $u_t=\mathcal O(t^{\alpha-1})$ and graded meshes |
| 4    | adaptive time stepping                                                                      | $E(t)$ against $\mathcal E_\alpha(t)$ |

The expected rate in test one is $1+\alpha$, the order of the L1$_R$ formula itself; in the presence of a $t^{\sigma}$ initial singularity, graded meshes recover it. The steps $\tau=0.8$ and $\tau=1.0$ used in test two are **large**, and the maximum norm still does not exceed $1$, matching how permissive the step condition in the theorem is.

**Test four is the most instructive figure in the paper, because it turns the warning attached to paper 40 into a visible phenomenon**: under adaptive steps, the original energy $E(t)$ and the variational energy $\mathcal E_\alpha(t)$ are plotted together, and $\mathcal E_\alpha$ decays monotonically while $E$ **need not**. That is the numerical face of paper 40's statement that an integral law does not imply pointwise decay, and it shows that changing the energy object here is a necessity rather than a rhetorical convenience.

The gap between theory and experiment is this: the theorems cover energy stability, solvability and the maximum bound principle, but **no error estimate is given** (that is what paper 43 does under the Alikhanov formula), so the $1+\alpha$ rate measured in test one is empirical. The paper also cautions that the "new Caputo derivative" induced by the DOC kernels is an indirect approximation whose accuracy differs from the direct L1 formula.

### Relation to the others

This is the direct sequel to papers 40 and 43 on the fractional side, and the answer to the open question posed in the conclusions of [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|paper 48]]: "develop nonuniform BDF2-type schemes for the time-fractional phase field equations". It replaces paper 40's integral law $E[u](T)\le E[u](0)$ by a genuine **differential** law for a redefined variational energy $\mathcal E_\alpha$ that limits to the classical law as $\alpha\to1$. It uses the DOC kernels introduced for BDF2 by paper 52 and Liao-Zhang, but in a new role: not to strip a multistep operator, but to **transform between two discrete fractional derivatives**. Its DOC monotonicity proof explicitly invokes the algebraic criteria of paper 74. Taken together with paper 43 it forms a pair: "maximum principle only, second order, Alikhanov" (43) versus "energy plus maximum bound, no restriction on the ratio, L1$_R$/Crank-Nicolson" (57).

## The progression across the three papers

| No. | Form of the energy law              | Step restriction                                | Main analytical device                         | Main numerical output                         |
| --- | ----------------------------------- | ----------------------------------------------- | ---------------------------------------------- | --------------------------------------------- |
| 40  | integral, accumulated over $[0,T]$  | none at the continuous level; stabilisation condition discretely | positivity of the fractional convolution kernel | $-\alpha/3$ coarsening power law (empirical)  |
| 43  | not proved; only the maximum principle | M1: $\rho=7/4$, plus a $\tau$-$h$ coupling      | Alikhanov kernel estimates, DCC kernels, Grönwall | the order $\min\{\gamma\sigma,2\}$ verified on graded meshes |
| 57  | differential, for a variational energy | none for the energy law; a $\tau$ bound for the maximum principle | rewriting, L1$_R$ kernel positivity, DOC transformation | $\mathcal E_\alpha$ monotone while $E$ need not be |

The move from paper 40 to paper 57 deserves its own summary: **when a law comes out in the wrong form, do not weaken the conclusion, change the object.** Paper 40 obtains an integral inequality for the original energy; paper 57 does not try to strengthen it but constructs a new energy $\mathcal E_\alpha$ satisfying a differential law that degenerates to the classical one as $\alpha\to1$. The price is that this energy contains a fractional integral term and is therefore not the original energy — and test four of paper 57 shows exactly that the original energy really can fail to be monotone, so the price is unavoidable.

## Coverage check

| Item                                                    | Paper | Status                                                                    |
| ------------------------------------------------------- | ----- | ------------------------------------------------------------------------- |
| Three fractional phase-field models and their energies  | 40    | Allen-Cahn, Cahn-Hilliard, both MBE variants                              |
| Full proof route for continuous kernel positivity       | 40    | symmetrisation, Liouville semigroup, Fourier transform, origin of the $\cos/\sin$ constants |
| Integral rather than differential, and not strengthenable | 40  | form of the conclusion, the paper's own warning, consequence for paper 57 |
| L1 kernels, rearrangement identity, discrete positivity  | 40   | kernel formula, complete monotonicity, $M$-matrix comparison, $s_n$       |
| Continuous and discrete energy theorems with conditions  | 40   | three models, stabilisation conditions, $\lambda_{\max}\le1/8$, $S\ge\gamma/(16\varepsilon)$ |
| Discrete maximum principle and its step condition        | 40   | the $(b_0-b_1)/\tau+S$ condition and why it is not unconditional          |
| Numerical experiments and the $-\alpha/3$ power law      | 40   | Fourier-Galerkin, SOE, three tests, flagged as empirical                  |
| Reciprocal step-ratio conventions                        | 43   | $\rho_k=1/r_{k+1}$ and the citation hazard                                |
| Alikhanov formula, two-index kernels, three estimates    | 43   | definitions of $a,b,A$, $\pi_A=11/4$, misprints flagged                   |
| DCC kernels, complementarity, fractional Grönwall        | 43   | definitions, $\equiv1$, power-type estimate, full Grönwall statement      |
| Solvability, maximum principle, sharp convergence        | 43   | all three theorems with hypotheses                                        |
| Convolution structure of the error and superconvergence  | 43   | $G_{\mathrm{loc}}$, $G_{\mathrm{his}}$, the $n=1$ comparison              |
| Adaptive strategy and three tests                        | 43   | SOE tolerance, guard, order table, negative test                          |
| Rewriting and relocating the nonlocality                 | 57   | semigroup identity, rewritten equation, why the test function works       |
| Variational energy and differential law                  | 57   | Riemann-Liouville inequality, $\mathcal E_\alpha$, the $\alpha\to1$ limit |
| L1$_R$ kernels and their positivity                      | 57   | definitions, sign pattern, discrete positivity, no step-ratio bound       |
| Constructive design of the nonlinear term                | 57   | form of $H$ and $H(a,b)(a-b)\ge F(a)-F(b)$                                |
| DOC kernels and the reversible transformation            | 57   | orthogonality, complementarity, monotonicity (citing paper 74), equivalent Caputo form |
| Four theorems and four tests                             | 57   | unconditional energy law, asymptotic preservation, solvability, maximum bound; the figure where $E$ is not monotone |

## Sources for this page

- T. Tang, H. Yu, and T. Zhou, [_On energy dissipation theory and numerical stability for time-fractional phase-field equations_](https://doi.org/10.1137/18M1203560), SIAM J. Sci. Comput. 41(6) (2019), pp. A3757-A3778 (preprint [arXiv:1808.01471](https://arxiv.org/abs/1808.01471)).
- H.-l. Liao, T. Tang, and T. Zhou, [_A second-order and nonuniform time-stepping maximum-principle preserving scheme for time-fractional Allen-Cahn equations_](https://doi.org/10.1016/j.jcp.2020.109473), J. Comput. Phys. 414 (2020), 109473 (preprint [arXiv:1909.10216](https://arxiv.org/abs/1909.10216)).
- H.-l. Liao, T. Tang, and T. Zhou, [_An energy stable and maximum bound preserving scheme with variable time steps for time fractional Allen-Cahn equation_](https://doi.org/10.1137/20M1384105), SIAM J. Sci. Comput. 43(5) (2021), pp. A3503-A3526 (preprint [arXiv:2012.10740](https://arxiv.org/abs/2012.10740)).
