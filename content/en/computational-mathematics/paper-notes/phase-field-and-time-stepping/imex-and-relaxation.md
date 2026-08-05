---
title: IMEX and Relaxation Schemes
description: Papers 78, 91 and 104 - treating the nonlinearity explicitly while keeping the energy argument
lang: en
translation: computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation
tags:
  - paper-notes
  - phase-field
  - implicit-explicit-schemes
---

## One question running through the three papers

Implicit-explicit schemes treat the stiff linear part implicitly and the nonlinearity explicitly, so no inner nonlinear iteration is needed per stage or per step. That is the entire reason they are popular in phase-field computation. The cost is concentrated in one place.

**A fully implicit energy argument works because the nonlinear term and the difference being tested share a sign.** Testing $-F'(u^n)$ with $\nabla_\tau u^n$ produces the difference of $F$ between two levels, and either convexity or an algebraic identity converts that into a telescoping cancellation plus a non-negative remainder. Move $F'$ to the explicit side and the level at which it is evaluated no longer matches the difference being tested, so that shared sign is gone: the explicit term supplies only a **sign-indefinite** contribution that has to be dominated by dissipation from elsewhere.

So the question becomes: **how much dissipation is needed, and what object decides that "how much"?** The three papers here answer it three ways, and their answers are three projections of the same structural question.

- **Paper 91** treats implicit-explicit Runge-Kutta methods. The object carrying the difficulty is a **differentiation matrix** $D(z)$ indexed by stages rather than time levels; it is the stage-index version of the discrete orthogonal convolution kernels of the [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|BDF family]]. The paper designs the Butcher tableau so that $D$ becomes **completely independent of the mesh parameters**.
- **Paper 104** treats implicit-explicit multistep methods. The object carrying the difficulty is the composition of three kernel sequences ($\vec a$ for the difference, $\vec b$ implicit, $\vec c$ explicit). The paper takes extrema of their **semi-generating functions** on the unit circle: three numbers, and the stability criterion falls out.
- **Paper 78** leaves the time discretisation alone and changes the **energy**. Quadratising the nonlinear potential through an auxiliary variable makes the scheme linear and the energy argument nearly trivial — at the price that what dissipates is a modified energy.

What the three share is that **none of them tries to prove directly that the explicit term contributes with a definite sign; each substitutes a controllable algebraic object instead.** That is exactly the stance of the [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|variable-step BDF page]]: locate the object carrying the whole difficulty, then build independent criteria for it.

## 91: making the average dissipation rate independent of the discretisation parameters

### The idea

Implicit-explicit Runge-Kutta methods applied to phase-field models suffer from two defects, one theoretical and one practical, and this paper removes both at once.

**The theoretical defect.** Proving that the scheme inherits the **original** energy dissipation law requires knowing first that the **stage solutions are uniformly bounded in the maximum norm** — the stabilisation parameter $\kappa$ has to dominate $\|F''\|_\infty$, and $F''(\Phi)=3\Phi^2-1$ is bounded only when $\Phi$ is. The earlier unified framework sidestepped that by simply **assuming** the nonlinear bulk $F'$ is globally Lipschitz, which the Cahn-Hilliard quartic $F(\Phi)=\frac14(\Phi^2-1)^2$ is not. This is not a minor technical slip: the Cahn-Hilliard equation has **no** maximum bound principle, so unlike Allen-Cahn there is no invariant interval handing over an $\ell^\infty$ bound for free.

**The practical defect.** The **average dissipation rate** of a generic implicit-explicit Runge-Kutta method depends on $\tau_n\overline\lambda_{\mathrm{ML}}$, the time step times the average eigenvalue of the discrete operator. So as soon as an adaptive algorithm takes a large step, **the method's effective dissipation changes** and the computed energy curve drifts with $\tau_{\max}$. For a method sold on adaptivity that is fatal in practice: one cannot tell whether a feature of the energy curve is physical or an artefact of the step-size cap.

**The paper's key judgement is that both defects are governed by the same matrix.** Written in differential form the scheme produces a lower triangular matrix $D(z)=D_{\mathrm E}-zD_{\mathrm{EI}}$, with $z$ a quantity of the form $-\tau_n\overline\lambda_{\mathrm{ML}}$. Its role in the stage index is **exactly** that of the DOC kernels in the time-level index: it inverts the multi-stage operator back to stage-wise differences so that an energy argument can proceed. The average dissipation rate is a combination of traces of $D$, so

$$
\mathcal R\ \text{is mesh-independent}
\quad\Longleftrightarrow\quad
D\ \text{is independent of}\ z
\quad\Longleftrightarrow\quad
D_{\mathrm{EI}}=\mathbf 0 .
$$

**The last of these is an algebraic equation on the Butcher tableau and can simply be solved.** That is what "refined" means here. And once $D_{\mathrm R}$ is a matrix of pure numbers, its minimum eigenvalue is mesh-independent too, which is precisely what the bootstrap for the uniform maximum-norm bound needs — **one condition kills both defects.**

### Setting

On a periodic $\Omega\subseteq\mathbb R^2$, the Ginzburg-Landau free energy and its $H^{-1}$ gradient flow:

$$
E[\Phi]=\int_\Omega\Bigl[\tfrac{\epsilon^2}{2}|\nabla\Phi|^2+F(\Phi)\Bigr]\mathrm d\mathbf x,
\qquad F(\Phi):=\tfrac14(\Phi^2-1)^2,
\qquad 0<\epsilon<1,
$$

$$
\partial_t\Phi=\Delta\bigl[F'(\Phi)-\epsilon^2\Delta\Phi\bigr],
\qquad
(\Phi(t),1)=(\Phi(t_0),1),
$$

$$
\frac{\mathrm dE}{\mathrm dt}
=\Bigl(\frac{\delta E}{\delta\Phi},\partial_t\Phi\Bigr)_{L^2}
=-\bigl((-\Delta)^{-1}\partial_t\Phi,\partial_t\Phi\bigr)_{L^2}\le0 .
$$

With a stabilisation parameter $\kappa\ge0$,

$$
L_\kappa\Phi:=-\epsilon^2\Delta\Phi+\kappa\Phi,
\qquad
f_\kappa(\Phi):=\kappa\Phi-F'(\Phi),
\qquad
\partial_t\Phi=\Delta\bigl[L_\kappa\Phi-f_\kappa(\Phi)\bigr].
$$

The $s$-stage implicit-explicit Runge-Kutta scheme on a nonuniform mesh $0=t_0<\cdots<t_N=T$ with $\tau_n=t_n-t_{n-1}$ is

$$
u_h^{n,i}=u_h^{n,1}
+\tau_n\sum_{j=1}^{i}a_{ij}\Delta_hL_{\kappa,h}u_h^{n,j}
-\tau_n\sum_{j=1}^{i-1}\hat a_{ij}\Delta_h f_\kappa(u_h^{n,j}),
\qquad
u_h^{n,1}:=\phi_h^{n-1},\ \ \phi_h^{n}:=u_h^{n,s},
$$

with $c_1=0$, $c_s=1$; the implicit part $A$ a **stiffly accurate DIRK with explicit first stage** ("first same as last"), $\widehat A$ strictly lower triangular (explicit), and the **canopy node condition** $\hat{\mathbf c}=\mathbf c$ (equivalently $A\mathbf 1=\widehat A\mathbf 1$), which makes the method consistent at all stages and preserves the equilibria $L_{\kappa,h}\phi_h^{*}=f_\kappa(\phi_h^{*})$. Space is discretised by a **Fourier pseudo-spectral** method, giving $\Delta_h$ and $L_{\kappa,h}$.

The order conditions up to third order are: first order $\mathbf b^T\mathbf 1=\hat{\mathbf b}^T\mathbf 1=1$; second order $\mathbf b^T\mathbf c=\hat{\mathbf b}^T\mathbf c=\frac12$; third order $\mathbf b^T\mathbf c^{.2}=\hat{\mathbf b}^T\mathbf c^{.2}=\frac13$ and $\mathbf b^TA\mathbf c=\hat{\mathbf b}^T\widehat A\mathbf c=\frac16$, plus the two **coupling** conditions $\mathbf b^T\widehat A\mathbf c=\hat{\mathbf b}^TA\mathbf c=\frac16$.

**Step-ratio restriction: none.** The mesh is arbitrary and nonuniform, $\tau_n$ is free, and **no step-ratio condition of any kind is imposed** — which is what "robust time adaptability" means.

> [!warning] The $\frac{1+\sqrt2}{4}$ that appears in this paper is not a step ratio
> The numeral $\frac{1+\sqrt2}{4}$ does occur in this paper, but as the Butcher coefficient $a_{33}$ of the **comparison** method IERK(2,3) taken from Liao-Wang-Wen. It is a Runge-Kutta coefficient, **not a step ratio**. The thresholds $1+\sqrt2$, $3.561$ and $4.8645$ do **not** appear in this paper and it does not need them. Connecting this $\frac{1+\sqrt2}{4}$ to condition S0 of [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|paper 48]] is a misreading.

### Derivation

**Step one: differential form and the differentiation matrix.** Write $s_{\mathrm I}:=s-1$, $\delta_\tau u^{n,\ell+1}:=u^{n,\ell+1}-u^{n,\ell}$, $u^{n,\ell+\frac12}:=(u^{n,\ell+1}+u^{n,\ell})/2$, and

$$
A_{\mathrm I}:=(a_{i+1,j+1})_{i,j=1}^{s_{\mathrm I}},
\qquad
A_{\mathrm E}:=(\hat a_{i+1,j})_{i,j=1}^{s_{\mathrm I}},
\qquad
E_{s_{\mathrm I}}:=(1_{i\ge j})\ \text{(lower-triangular all-ones)} .
$$

The scheme is then equivalent to

$$
\sum_{\ell=1}^{i}d_{i\ell}\bigl(\tau_n\Delta_hL_{\kappa,h}\bigr)\,\delta_\tau u_h^{n,\ell+1}
=\tau_n\Delta_h\Bigl[L_{\kappa,h}u_h^{n,i+\frac12}-f_\kappa(u_h^{n,i})\Bigr],
\qquad 1\le i\le s_{\mathrm I},
$$

$$
D(z):=D_{\mathrm E}-zD_{\mathrm{EI}},
\qquad
D_{\mathrm E}:=A_{\mathrm E}^{-1}E_{s_{\mathrm I}},
\qquad
D_{\mathrm{EI}}:=A_{\mathrm E}^{-1}A_{\mathrm I}E_{s_{\mathrm I}}-E_{s_{\mathrm I}}+\tfrac12 I_{s_{\mathrm I}} .
$$

A lower triangular $D$ is called positive (semi-)definite when its symmetric part $\mathcal S(D)=(D+D^T)/2$ is. **This $D$ plays exactly the role the DOC kernels play in the BDF papers**, with the index changed from time levels to stages; the corresponding orthogonality identity is

$$
\sum_{i=j}^{k}d^{(R)}_{k,i}\,\underline{\hat a}_{i+1,j}\equiv\delta_{kj},
\qquad\text{so that}\qquad
\sum_{i=1}^{k}d^{(R)}_{k,i}\sum_{j=1}^{i}\underline{\hat a}_{i+1,j}v^j\equiv v^k,
\qquad
(\underline{\hat a}_{i+1,j}):=E_{s_{\mathrm I}}^{-1}A_{\mathrm E}.
$$

**Step two: the average dissipation rate.** Under positive semi-definiteness of $D_{\mathrm E}$ and $D_{\mathrm{EI}}$, the stage energy law comes with

$$
\mathcal R=\frac1{s_{\mathrm I}}\mathrm{tr}(D_{\mathrm E})
+\frac1{s_{\mathrm I}}\mathrm{tr}(D_{\mathrm{EI}})\,\tau_n\overline\lambda_{\mathrm{ML}}
=\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}\frac1{\hat a_{k+1,k}}
+\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}
\Bigl(\frac{a_{k+1,k+1}}{\hat a_{k+1,k}}-\frac12\Bigr)\tau_n\overline\lambda_{\mathrm{ML}}
\ \ge 0,
$$

where $\overline\lambda_{\mathrm{ML}}>0$ is the average eigenvalue of the symmetric positive definite matrix $-\Delta_hL_{\kappa,h}$. A method is "good" when $\mathcal R$ stays as close to $1$ as possible over a wide range of $\tau_n\overline\lambda_{\mathrm{ML}}$.

**Step three: the refinement condition.** $\mathcal R$ is independent of $\tau_n\overline\lambda_{\mathrm{ML}}$ **if and only if** $D_{\mathrm{EI}}=\mathbf 0$, that is

$$
A_{\mathrm E}^{-1}A_{\mathrm I}E_{s_{\mathrm I}}-E_{s_{\mathrm I}}+\tfrac12I_{s_{\mathrm I}}=\mathbf 0
\qquad\Longleftrightarrow\qquad
A_{\mathrm I}=A_{\mathrm E}P_{s_{\mathrm I}},
\qquad
P_{s_{\mathrm I}}:=I_{s_{\mathrm I}}-\tfrac12E_{s_{\mathrm I}}^{-1} .
$$

Then $D_{\mathrm R}:=(d^{(R)}_{ij})=D_{\mathrm E}=A_{\mathrm E}^{-1}E_{s_{\mathrm I}}$ is **independent of $z$ and hence of both mesh parameters**, with

$$
\mathcal R_{\mathrm R}=\frac1{s_{\mathrm I}}\sum_{k=1}^{s_{\mathrm I}}\frac1{\hat a_{k+1,k}} .
$$

The scheme collapses to the compact refined form

$$
u_h^{n,i+1}=u_h^{n,1}+\tau_n\sum_{j=1}^{i}\hat a_{i+1,j}\Delta_h
\Bigl[L_{\kappa,h}u_h^{n,j+\frac12}-f_\kappa(u_h^{n,j})\Bigr],
$$

equivalently $\delta_\tau u_h^{n,i+1}=\tau_n\sum_{j=1}^{i}\underline{\hat a}_{i+1,j}\Delta_h[L_{\kappa,h}u_h^{n,j+\frac12}-f_\kappa(u_h^{n,j})]$.

**Step four: a structural exclusion.** The canopy condition forces the first implicit column

$$
\mathbf a_1=\bigl(\tfrac12\hat a_{21},\dots,\tfrac12\hat a_{s1}\bigr)^T\ne\mathbf 0,
$$

so such methods are necessarily **Lobatto-type and never Radau- or ARS-type**. The paper states this as a proposition: there is no Radau-type or ARS-type implicit-explicit Runge-Kutta method whose average dissipation rate is independent of $\tau_n\overline\lambda_{\mathrm{ML}}$. Results of this kind are valuable in method design — the conclusion is not "we chose Lobatto-type" but "this property forces Lobatto-type". Correspondingly the implicit part has stage order two and is **not necessarily algebraically stable**.

**Step five: three concrete families.**

| Method                    | Construction                                                                                                                                                                                                                                                                                | Parameter window and rate                                                                                             |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| R-IERK(1,2)               | take $\theta=\frac12$ in the two-stage family IERK(1,2;$\theta$), for which $D^{(1,2)}(z)=1-z(\theta-\frac12)$ is positive definite for $z\le0$ exactly when $\theta\ge\frac12$, and $\mathcal R^{(1,2)}(\theta)=1+(\theta-\frac12)\tau_n\overline\lambda_{\mathrm{ML}}$                    | $\mathcal R^{(1,2)}_{\mathrm R}=1$ exactly                                                                            |
| R-IERK(2,4;$c_2$)         | no three-stage second-order method exists; the four-stage one-parameter family sets $\hat a_{32}=\hat a_{43}=c_2$ and $\hat a_{42}=\frac1{2c_2}-c_3$                                                                                                                                        | $\frac{\sqrt3-1}{2}\le c_2<1$ or $1<c_2\le\frac{2+\sqrt6}{2}$, with $\mathcal R^{(2,4)}_{\mathrm R}(c_2)=\frac1{c_2}$ |
| R-IERK(3,6;$\hat a_{52}$) | no four-stage third-order method exists and no positive-(semi-)definite five-stage one was found; the six-stage family fixes $c_2=1$, $c_3=\frac45$, $c_4=\frac7{10}$, $c_5=\frac{12}{25}$, $\hat a_{32}=\frac65$, $\hat a_{51}=\frac{17\hat a_{52}}{103}-\frac{86756827361}{181963162950}$ | $0.664767<\hat a_{52}<0.751947$                                                                                       |

R-IERK(1,2) is the Crank-Nicolson-type scheme

$$
\delta_\tau\phi_h^{n}=\tau_n\Delta_h\Bigl[\tfrac12L_{\kappa,h}(\phi_h^{n}+\phi_h^{n-1})
-f_\kappa(\phi_h^{n-1})\Bigr].
$$

In the four-stage family $c_3$ is a root of

$$
c_3^2-\Bigl(\frac1{2c_2}+c_2\Bigr)c_3+(c_2-1)^2=0,
$$

explicitly $c_3^{*}=\frac32$ when $c_2=1$ and otherwise

$$
c_3^{*}=\frac1{4c_2}+\frac{c_2}{2}
-\frac{\sqrt{(-2c_2^2+4c_2+1)(6c_2^2-4c_2+1)}}{4c_2}>0,
\qquad 0<c_2\le\tfrac{2+\sqrt6}{2}.
$$

Positive definiteness is decided by two explicit determinants:

$$
\mathrm{Det}\,\mathcal S\bigl(D^{(2,4)}_{\mathrm R,2};c_2\bigr)=\frac{c_3(4c_2-c_3)}{4c_2^4}>0,
\qquad
\mathrm{Det}\,\mathcal S\bigl(D^{(2,4)}_{\mathrm R};c_2\bigr)
=\frac{1}{8c_2^7}(c_2-1)^2(2c_2^2+2c_2-1)>0 .
$$

The six-stage family has rate

$$
\mathcal R^{(3,6)}_{\mathrm R}(\hat a_{52})
=\frac{36392632590}{123664285500\,\hat a_{52}+89552314349}
+\frac{219055887768899}{156795586342500}
\approx\frac{1.39708\,\hat a_{52}+1.30599}{\hat a_{52}+0.724157}.
$$

> [!note] Two boundaries the paper draws itself
> First, $c_2=1$ gives the optimal rate $\mathcal R=1$ but makes $\mathcal S(D_{\mathrm R}^{(2,4)};1)$ **singular** ($\lambda_{\min}=0$), so it is **excluded from the theory** (Remark 3.1) even though it performs best numerically. That is a genuine theory-versus-practice gap and the paper does not paper over it. Second, whether a six-stage method attaining $\mathcal R^{(3,6)}_{\mathrm R}=1$ exists is left open; fourth order would need at least nine stages (28 order conditions) and is not attempted.

**Step six: what replaces the maximum principle.** The Cahn-Hilliard equation has **no** maximum bound principle, so **this paper has no maximum-principle theorem**. What stands in for it is a uniform maximum-norm bound on the stage solutions, obtained by an **updated time-space error splitting**:

1. energy arguments with a rough setting of the stage defects give $H^{m+4}$ regularity of the **time-discrete** stage solutions $U^{n,i}$;
2. the fully discrete stage error is split as $U^{n,\ell}-u_h^{n,\ell}=(U^{n,\ell}-U_M^{n,\ell})+e_h^{n,\ell}$ with $U_M=P_MU$ the $L^2$ projection;
3. a rough maximum-norm bound $\hat c_4^{*}/\epsilon^2$ is bootstrapped.

Two lemmas make this work. **Lemma 4.3 (the DOC-matrix inequality)**: if $D_{\mathrm R}$ and $D_{\mathrm R}^{-1}$ are both positive definite, with $\lambda_{\min}$ and $\sigma_{\min}$ the minimum eigenvalues of $\mathcal S(D_{\mathrm R})$ and $\mathcal S(D_{\mathrm R}^{-1})$, then for arbitrary sequences $\{v^j\},\{u^j\}$,

$$
\sum_{i=1}^{k}\sum_{j=1}^{i}d^{(R)}_{i,j}v^jv^i\ge\lambda_{\min}\sum_i(v^i)^2,
\qquad
\sum_{i}\sum_{j\le i}\underline{\hat a}_{i+1,j}v^jv^i\ge\sigma_{\min}\sum_i(v^i)^2,
$$

$$
\sum_i\sum_{j\le i}d^{(R)}_{i,j}v^ju^i\le\frac1{\sigma_{\min}}\sum_i|v^i||u^i|,
\qquad
\sum_i\sum_{j\le i}\underline{\hat a}_{i+1,j}v^ju^i\le\frac1{\lambda_{\min}}\sum_i|v^i||u^i| .
$$

**Here $\lambda_{\min}$ and $\sigma_{\min}$ are pure numbers because $D_{\mathrm R}$ is mesh-independent — the second payoff of the refinement condition, this time inside the convergence proof.**

**Lemma 4.4 (local handling of the nonlinearity)**: using

$$
F'(v)-F'(w)=(v-w)\int_0^1F''\bigl[\gamma v+(1-\gamma)w\bigr]\,\mathrm d\gamma
$$

and the corresponding identity for $\delta_\tau F'$, the Lipschitz constant is localised to a ball. **This is the step at which the global Lipschitz assumption is dropped.**

The stage energy inequality itself reads (assuming the stage solutions bounded by $c_0$ in the maximum norm and $\kappa\ge\max_{\|\xi\|_\infty\le c_0}\|F''(\xi)\|_\infty$)

$$
\bigl(L_{\kappa,h}u^{n,i+\frac12}-f_\kappa(u^{n,i}),\delta_\tau u^{n,i+1}\bigr)
\le E[u^{n,i}]-E[u^{n,i+1}]
-\tfrac12\|\delta_\tau u^{n,i+1}\|^2
\bigl(\kappa-\max_{\xi_h\in\mathcal B_{n,i}}\|F''(\xi)\|_\infty\bigr).
$$

Additionally, **stage-wise volume conservation is exact**: $(u^{n,i+1},1)=(u^{n,1},1)=(\phi^0,1)$ for all $n$ and $i$.

### Theorems

**(Standing regularity assumption)** There are integers $m\ge1$, $p\ge1$ and a constant $c_\phi>0$ with

$$
\|\Phi^0\|_{H^{m+4}}+\sum_{k=0}^{2}\|\partial_t^{(k)}\Phi(t)\|_{H^{m+4-k}}
+\sum_{k=3}^{p+1}\|\partial_t^{(k)}\Phi(t)\|_{L^2}\le c_\phi,
\qquad 0<t<T .
$$

**(Lemma 2.1, stage energy law for a general IERK method)** If $D_{\mathrm E},D_{\mathrm{EI}}$ are positive (semi-)definite, the stage solutions are bounded by $c_0$ in the maximum norm, and $\kappa\ge\max_{\|\xi\|_\infty\le c_0}\|F''(\xi)\|_\infty$, then for $n\ge1$ and $1\le i\le s_{\mathrm I}$,

$$
E[u^{n,i+1}]-E[u^{n,1}]\le\frac1{\tau_n}\sum_{k=1}^{i}
\Bigl(\Delta_h^{-1}\delta_\tau u^{n,k+1},\ \sum_{\ell=1}^{k}
d_{k\ell}(\tau_n\Delta_hL_{\kappa,h})\,\delta_\tau u^{n,\ell+1}\Bigr),
$$

the right-hand side being measured by the non-negative average dissipation rate $\mathcal R$ displayed above.

**(Theorem 4.1, regularity of the time-discrete stage solutions)** Under the regularity assumption with $m\ge1$, and assuming **both $D_{\mathrm R}=A_{\mathrm E}^{-1}E_{s_{\mathrm I}}$ and $D_{\mathrm R}^{-1}$ are positive definite**, if the maximum time step $\tau$ is sufficiently small then there is $C_\phi>0$, **independent of the step sizes $\tau_n$**, with

$$
\|U^{n,i}\|_{H^{m+4}}+\bigl\|(U^{n,i}-U^{n,1})/\tau_n\bigr\|_{H^{m}}
\le C_\phi/\epsilon^{2},
\qquad 1\le n\le N,\ 2\le i\le s .
$$

**(Theorem 4.2, the original energy dissipation law at every stage)** Under the same hypotheses, with $h$ and $\tau$ sufficiently small and $\kappa\ge\max_{\|\xi\|_\infty\le\hat c_4^{*}/\epsilon^2}\|F''(\xi)\|_\infty$, the stage solutions $u_h^{n,i}$ are bounded in the maximum norm and

$$
E[u^{n,j+1}]-E[u^{n,1}]\le\frac1{\tau}\sum_{k=1}^{j}
\Bigl(\Delta_h^{-1}\delta_\tau u^{n,k+1},\ \sum_{\ell=1}^{k}
d^{(R)}_{k\ell}\,\delta_\tau u^{n,\ell+1}\Bigr),
\qquad 1\le n\le N,\ 1\le j\le s_{\mathrm I} .
$$

**What decays is the original energy $E$ — not a modified energy, and with no auxiliary variable — and there is no step-ratio restriction of any kind.** The $d^{(R)}_{k\ell}$ here are pure numbers, independent of the mesh.

**(Theorem 4.3, $L^2$ convergence at full accuracy)** Under the same hypotheses, with the stage-defect assumption $\zeta^{n,i+1}_{\mathrm R}=0$ for $1\le i\le s_{\mathrm I}$ and $\|\zeta^{n,s}_{\mathrm R}\|\le c_2\tau^{p}$, the solution $u_h^n$ converges in the $L^2$ norm with order $\mathcal O(\tau^{p}+h^{m})$. "Unconditional" here means that no coupling between $\tau$ and $h$ and no step-ratio restriction is required.

The paper states this is the **first** time the original energy dissipation law and unconditional $L^2$ convergence of implicit-explicit Runge-Kutta methods have been established for Cahn-Hilliard **without assuming global Lipschitz continuity of the nonlinear bulk**.

### Numerical experiments

Space is Fourier pseudo-spectral throughout, and the adaptive step-size rule used everywhere is

$$
\tau_{\mathrm{ada}}=\max\{\tau_{\min},\ \tau_{\max}/\Pi_\eta(\phi)\},
\qquad
\Pi_\eta(\phi):=\sqrt{1+\eta\|\partial_\tau\phi^n\|^2},
$$

with defaults $\eta=1000$, $\tau_{\min}=10^{-4}$ and $\tau_1=\tau_{\min}$.

**Example 5.1 (accuracy test).** Manufactured solution $\Phi(x,y;t)=e^{-t}\sin(\pi x)\sin(\pi y)$ on $\Omega=(0,2)^2$ with $\epsilon=0.2$, a $64\times64$ grid, $T=1$, $\kappa=4$, steps $\tau=2^{-k}/10$ for $0\le k\le9$, and error $e(\tau)=\max_n\|\Phi_h^n-\Phi(t_n)\|_\infty$. R-IERK(2,4;$c_2$) is second order and R-IERK(3,6;$\hat a_{52}$) is third order. **The parameter sensitivity is markedly asymmetric**: different $c_2$ give visibly different accuracy ($c_2=1$ is most accurate for $\tau<10^{-2}$), while different $\hat a_{52}$ give nearly identical solutions.

**Example 2.1 (the test that motivates the whole paper).** $\Omega=(-\pi,\pi)^2$, $\epsilon=0.1$, $\kappa=2.5$, $T=1000$, initial datum

$$
\Phi^0=\tfrac12\tanh(|x|+|y|+1)-e^{-5(|x|+|y|-2)^2}
+\tfrac12e^{-2(|x|-1)^2}+\tfrac1{10}\sin\bigl(e^{|y|-1}\bigr),
$$

with the reference produced by the Lobatto method IERK(2,3) ($a_{33}=\frac{1+\sqrt2}{4}$) at $\tau=10^{-4}$. For $\tau_{\max}\in\{0.01,0.05,0.1\}$: **the energy curve of R-IERK(1,2) ($\mathcal R=1$) stays close to the reference throughout, while IERK(1,2;$\theta$) with $\theta=1,\frac32$ drifts as $\tau_{\max}$ grows.** That is the practical defect from the introduction, demonstrated directly.

**Energy and efficiency tests for R-IERK(2,4;$c_2$).** With $c_2=\frac12,1,\frac32$ and $\tau_{\max}\in\{0.1,0.2,0.5\}$, all discrete energies decrease; the IERK(2,3) curves change significantly with $\tau_{\max}$ while the R-IERK curves are robust. Efficiency table for R-IERK(2,4;1) at $T=450$:

| Step-size setting    | CPU time    | time levels       |
| -------------------- | ----------- | ----------------- |
| fixed $\tau=10^{-4}$ | $6724.63$ s | $4.5\times10^{6}$ |
| $\tau_{\max}=0.1$    | $11.90$ s   | $11714$           |
| $\tau_{\max}=0.2$    | $6.42$ s    | $5840$            |
| $\tau_{\max}=0.5$    | $2.62$ s    | $2304$            |

**This table is the clearest record of theory bought back as performance among the three**: roughly three orders of magnitude, and it is only usable because the energy curve does not drift with $\tau_{\max}$. Mesh-independence of $\mathcal R_{\mathrm R}$ is not an aesthetic property; it is what makes large adaptive steps trustworthy.

**Energy tests for R-IERK(3,6;$\hat a_{52}$).** With $\hat a_{52}=\frac23,\frac7{10},\frac34$, adaptivity parameter $\eta=500$ and $\tau_{\max}\in\{0.2,0.5,0.8\}$; the reference is the Lobatto method IERK(3,5) with $a_{43}=-\frac35$ and $\mathcal R^{(3,5)}_{\mathrm L}=\frac54+\frac25\tau\overline\lambda_{\mathrm{ML}}$. The R-IERK energy curves are **indistinguishable** across all three values of $\tau_{\max}$, while the IERK(3,5) curves change significantly. Note that $\mathcal R^{(3,5)}_{\mathrm L}$ contains $\tau\overline\lambda_{\mathrm{ML}}$ explicitly, which is exactly why it drifts — a very clean contrast.

**Coarsening dynamics.** Solution profiles at $t=0,2,20,200,450,1000$ from R-IERK(3,6;$\frac34$) with $\eta=500$ and $\tau_{\max}=0.5$.

**What these experiments establish and what they do not.** They establish that the refined methods' energy curves really are independent of $\tau_{\max}$, so adaptivity can safely take large steps, with roughly three orders of magnitude of speed-up, and that the observed orders match the theorems. They do not establish sharpness of the theoretical boundary: $c_2=1$ is excluded by the theory (singular $\mathcal S(D_{\mathrm R})$) yet is the most accurate in the experiments, which suggests the positive-definiteness requirement may be sufficient rather than necessary at second order. The paper records the contradiction honestly.

### Relation to the others

The differentiation matrix $D$ is the **stage-index incarnation of the DOC kernels**: papers [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|48, 52, 58 and 67]] use them on the time-level index, this paper on the stage index; the paper credits the DOC technique to Liao-Ji-Wang-Zhang (2022), Liao-Ji-Zhang (2022, phase field crystal) and Liao-Zhang (2021).

**The sharpest methodological contrast in the collection sits here**: where papers 48 and 52 buy robustness for variable-step BDF2 by **restricting the step ratio** ($r_k<3.561$), this paper buys robustness for a **multi-stage** method by **removing the mesh dependence of the dissipation rate**, so **no step-ratio condition is needed at all**. The two routes address two symptoms of the same ailment.

Its disagreement with paper 78 is about the type of energy: what decays here is the **original** energy $E[\Phi]$, not a modified energy containing an auxiliary variable — the same stance as [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|papers 40 and 57]] and papers 48 and 52, and indeed this paper's introduction criticises SAV-based high-order schemes for establishing stability only with respect to a modified energy. Paper 104 treats the **multistep** version of the same problem, swapping the differentiation matrix for generating functions.

Finally, the "uniform boundedness of stage solutions" difficulty solved here is exactly what blocked the earlier Radau-type and ARS-type constructions of Fu-Tang-Yang and Shin-Lee-Lee, both of which required global Lipschitz continuity.

## 104: replacing a family-by-family search for multipliers

### The idea

For nonlinear parabolic equations, implicit-explicit multistep methods treat the stiff linear operator implicitly and the nonlinearity explicitly. Many families exist — weighted BDF (Li-Xie), modified BDF (Akrivis-Karakatsani), generalised BDF (Huang-Shen), NIMEX (Rosales-Seibold-Shirokoff-Zhou) — all designed to enlarge absolute stability regions. But third- and higher-order variants are **not A-stable**, so rigorous **discrete energy** stability had been a long-standing open problem.

Both available routes demand bespoke objects found **family by family**, and that is the real pain point this paper targets. The first is Nevanlinna-Odeh-type multipliers via Dahlquist G-stability: weighted, modified and generalised BDF and NIMEX would each need their own multiplier family, and multipliers are hard to find. The second is an implicit-part decomposition of Huang-Shen type,

$$
\sum_{j=0}^{k-1}b_{\mathrm G,j}^{(k,\beta)}v^{n-j}
=\eta_k(\beta)\sum_{j=0}^{k-1}c_{\mathrm G,j}^{(k,\beta)}v^{n-j}
+\sum_{j=0}^{k-1}d_j^{(k,\beta)}v^{n-j},
\qquad
\eta_2=\tfrac{\beta-1}{\beta},\
\eta_3=\tfrac{\beta-1}{\beta+1},\
\eta_4=\tfrac{\beta-1}{\beta+3},
$$

giving stability under $\eta_k(\beta)>\mu_0/\varpi$. Its limitations are concrete: it degenerates at $\beta=1$ and therefore **never covers plain BDF-$k$**; it exists only for $2\le k\le4$; its refined three-term version (Huang-Shen 2025) works only at the fixed $\beta_k=3,6,9$; and there is no analogue for generalised BDF5 or for weighted BDF, modified BDF and NIMEX.

**This paper replaces all of that with a single computation: three extremal values of three explicit rational functions on the unit circle.**

The observation behind it is surprisingly simple. After the multistep scheme is acted on by the DOC kernels, the implicit and explicit parts each become a **lower-triangular Toeplitz** matrix, and the **symmetric part** of a lower-triangular Toeplitz matrix has the classical Grenander-Szegő generating function — which turns out to be exactly the **real part** of the one-sided sum $a(\theta)=\sum_{k\ge0}a_ke^{\imath k\theta}$. So even though the kernels are causal (one-sided), the Toeplitz-Carathéodory criterion still applies, provided one takes real parts. That is what "semi-generating function" means.

With it in hand, every question about the three kernel sequences becomes an extremal problem for three complex rational polynomials on the unit circle, and since **products of lower-triangular Toeplitz matrices are again lower-triangular Toeplitz and commute**, the composite kernels' generating functions are just products and quotients of the individual ones. **The family-by-family search for multipliers is cancelled in one stroke.**

### Setting

On a Hilbert triple $V\subset H=H'\subset V'$ consider

$$
u_t+\varpi\mathcal Lu=\mathcal F(u),\qquad 0<t<T,\qquad u(0)=u^0\in H,\ \varpi>0,
$$

with $\mathcal L:V\to V'$ positive definite, self-adjoint, bounded and linear, and $\mathcal F$ possibly nonlinear. The norms are $\|v\|_H=\langle v,v\rangle^{1/2}$, $\|v\|_V=\langle\mathcal Lv,v\rangle^{1/2}$ and the dual $\|v\|_{\star}=\sup_{\|w\|_V=1}|\langle v,w\rangle|$. The standing hypothesis is a **local** Lipschitz condition on the ball $\mathcal B_{u(t)}=\{v\in V:\|v-u(t)\|_V\le1\}$:

$$
\|\mathcal F(v)-\mathcal F(w)\|_{\star}\le\mu_0\|v-w\|_V+\mu_1\|v-w\|_H,
\qquad \mu_0\in(0,\varpi),\ \mu_1\ \text{arbitrary} .
$$

The framework also covers a non-self-adjoint $\mathcal L=\mathcal L_s+\mathcal L_a$ by moving the low-order anti-self-adjoint part $\mathcal L_a$ into the explicit term.

A general $k$-step implicit-explicit multistep method on a **uniform mesh** $\tau=t_j-t_{j-1}$ reads

$$
\sum_{j=0}^{k-1}a_j^{(k)}\partial_\tau u^{n-j}
+\varpi\sum_{j=0}^{k}b_j^{(k)}\mathcal Lu^{n-j}
=\sum_{j=0}^{k-1}c_j^{(k)}\mathcal F(u^{n-j-1})+\mathfrak C^{(k)}_n(u^0),
\qquad n\ge1,
$$

with $a_0^{(k)},b_0^{(k)},c_0^{(k)}>0$ and starting corrections $\mathfrak C^{(k)}_n(u^0)$ (zero for $n\ge k$) assumed available, so that the scheme is $k$-th order consistent from the first step. **A method is therefore a triad $(\vec a^{(k)},\vec b^{(k)},\vec c^{(k)})$.** Extending the kernels by zero ($b_j^{(k)}=0$ for $j\ge k+1$, $a_j^{(k)}=c_j^{(k)}=0$ for $j\ge k$) gives the equivalent global convolution form

$$
\sum_{j=1}^{n}a^{(k)}_{n-j}\partial_\tau u^{j}
+\varpi\sum_{j=1}^{n}b^{(k)}_{n-j}\mathcal Lu^{j}
=\sum_{j=1}^{n}c^{(k)}_{n-j}\mathcal F(u^{j-1})+\mathfrak C^{(k)}_n(u^0).
$$

**Step-ratio restriction: none — the mesh is uniform. Maximum-principle argument: none.**

### Derivation

**Step one: DOC kernels and the global discrete energy method.** The DOC kernels of $\vec a^{(k)}$ are

$$
a_0^{(-1,k)}:=\frac{1}{a_0^{(k)}},
\qquad
a_j^{(-1,k)}:=-\frac{1}{a_0^{(k)}}\sum_{i=1}^{j}a^{(-1,k)}_{j-i}a_i^{(k)}\ (j\ge1),
$$

satisfying two-sided orthogonality

$$
\sum_{\ell=j}^{n}a^{(-1,k)}_{n-\ell}a^{(k)}_{\ell-j}\equiv\delta_{nj}
=\sum_{\ell=j}^{n}a^{(k)}_{n-\ell}a^{(-1,k)}_{\ell-j},
\qquad 1\le j\le n .
$$

Acting with them turns the scheme into the "differential" form that is the starting point of the energy argument:

$$
\partial_\tau u^{n}+\varpi\sum_{\ell=1}^{n}\hat b^{(k)}_{n-\ell}\mathcal Lu^{\ell}
=\sum_{\ell=1}^{n}\hat c^{(k)}_{n-\ell}\mathcal F(u^{\ell-1})
+\sum_{\ell=1}^{n}a^{(-1,k)}_{n-\ell}\mathfrak C^{(k)}_\ell(u^0),
$$

with **composite kernels**

$$
\hat b_j^{(k)}:=\sum_{i=0}^{j}a^{(-1,k)}_{j-i}b_i^{(k)},
\qquad
\hat c_j^{(k)}:=\sum_{i=0}^{j}a^{(-1,k)}_{j-i}c_i^{(k)} .
$$

In matrix form, with $A_{L,k},B_{L,k},C_{L,k}$ the lower-triangular Toeplitz matrices of the three coefficient sequences, $A^{(-1)}_{L,k}=A_{L,k}^{-1}$ and

$$
\widehat B_{L,k}=A_{L,k}^{-1}B_{L,k},
\qquad
\widehat C_{L,k}=A_{L,k}^{-1}C_{L,k},
$$

since **products of lower-triangular Toeplitz matrices are again lower-triangular Toeplitz and commute** — a trivial fact that is nonetheless the technical basis for everything that follows. Because no information is discarded, the authors call this the **global discrete energy method**.

**Step two: the semi-generating function (the paper's device, Lemma 2.1).** For a real sequence $\{a_0,a_1,\dots\}$ with $a_j=0$ for $j<0$, define

$$
a(\theta):=\sum_{k=0}^{\infty}a_ke^{\imath k\theta}\in L^2([0,2\pi)) .
$$

For the quadratic form $Q_n:=\sum_{k=1}^{n}w_k\sum_{j=1}^{k}a_{k-j}w_j$ with lower-triangular Toeplitz matrix $P_{L,n}$:

1. $Q_n$ is positive definite **if and only if** $\Re[a(\theta)]>0$ on $[0,2\pi)$;
2. $\min_\theta\Re[a(\theta)]\le\lambda_j(Q_n)\le\max_\theta\Re[a(\theta)]$;
3. the eigenvalues are **equally distributed** as $\Re[a(2\pi j/n)]$, that is $\lim_{n\to\infty}\frac1n\sum_{j=0}^{n-1}\bigl(\lambda_j(Q_n)-\Re[a(2\pi j/n)]\bigr)=0$.

The proof rests on a single observation: the classical Grenander-Szegő generating function of the **symmetric part** $\mathcal S(P_{L,n})$ is

$$
\mathrm g(\theta)=a_0+\sum_{k\ge1}a_k\cos k\theta=\Re[a(\theta)],
$$

so **the real part of the one-sided (causal) sum equals the classical two-sided generating function**. That is the extension of Toeplitz-Carathéodory to one-sided sequences.

**Step three: composition rules (Lemma 2.2).** (i) If $\hat b_j=\sum_{k=0}^{j}a_{j-k}b_k$ then $\hat b(\theta)=a(\theta)b(\theta)$. (ii) If $\{\xi_j\}$ are the DOC kernels of $\{a_j\}$ then $\xi(\theta)=1/a(\theta)$. Hence, writing

$$
a^{(k)}(\theta)=\sum_{j=0}^{k-1}a_j^{(k)}e^{\imath j\theta},
\qquad
b^{(k)}(\theta)=\sum_{j=0}^{k}b_j^{(k)}e^{\imath j\theta},
\qquad
c^{(k)}(\theta)=\sum_{j=0}^{k-1}c_j^{(k)}e^{\imath j\theta},
$$

one gets $a^{(-1,k)}(\theta)=1/a^{(k)}(\theta)$, $\hat b^{(k)}(\theta)=b^{(k)}(\theta)/a^{(k)}(\theta)$ and $\hat c^{(k)}(\theta)=c^{(k)}(\theta)/a^{(k)}(\theta)$. **Everything reduces to three complex rational polynomials on the unit circle.**

**Step four: the three extremal constants (Lemma 2.4).**

$$
\sigma_{\mathrm F}^{(k)}=\max_{\theta\in[0,2\pi)}\Bigl|\frac{1}{a^{(k)}(\theta)}\Bigr|,
\qquad
\sigma_{\mathrm E}^{(k)}=\max_{\theta\in[0,2\pi)}\Bigl|\frac{c^{(k)}(\theta)}{a^{(k)}(\theta)}\Bigr|,
\qquad
\lambda_{\mathrm I}^{(k)}=\min_{\theta\in[0,2\pi)}\Re\Bigl[\frac{b^{(k)}(\theta)}{a^{(k)}(\theta)}\Bigr],
$$

called respectively the **perturbation amplification factor**, the **nonlinear amplification factor** and the **dissipation preserving factor**. Their role is

$$
\|A_{L,k}^{-1}\|_{\ell^2}\le\sigma_{\mathrm F}^{(k)},
\qquad
\|\widehat C_{L,k}\|_{\ell^2}\le\sigma_{\mathrm E}^{(k)},
\qquad
\lambda\bigl(\mathcal S(\widehat B_{L,k})\bigr)>\lambda_{\mathrm I}^{(k)},
$$

which unfold, for arbitrary sequences $\{v^i\},\{u^i\}$, into

$$
\sum_{i=1}^{n}\sum_{j=1}^{i}\hat b^{(k)}_{i-j}v^jv^i\ge\lambda_{\mathrm I}^{(k)}\sum_{i=1}^{n}|v^i|^2,
$$

$$
\sum_{i}\sum_{j\le i}a^{(-1,k)}_{i-j}v^ju^i
\le\sigma_{\mathrm F}^{(k)}\Bigl(\sum_i|v^i|^2\Bigr)^{1/2}\Bigl(\sum_i|u^i|^2\Bigr)^{1/2},
$$

plus the third with $a^{(-1,k)}$ replaced by $\hat c^{(k)}$ and $\sigma_{\mathrm F}^{(k)}$ by $\sigma_{\mathrm E}^{(k)}$. The proof ingredients are Grenander-Szegő plus Cauchy interlacing for the first, and Parseval on $\ell^2(\mathbb N)$ with $Y(\theta)=X(\theta)\overline{a(\theta)}$ for the spectral-norm bound $\|P_L\|_{\ell^2}\le\max_\theta|a(\theta)|$.

**Step five: a sharp normalisation (Lemma 2.5 and Corollary 3.4).** Consistency means $a^{(k)}(0)=b^{(k)}(0)=c^{(k)}(0)=1$, whence

$$
\sigma_{\mathrm F}^{(k)}\ge1,
\qquad
\sigma_{\mathrm E}^{(k)}\ge1,
\qquad
\lambda_{\mathrm I}^{(k)}\le1,
$$

with equality throughout in the **unique** case of the implicit-explicit Euler scheme ($k=1$). **This normalisation is what makes the three factors comparable quantities**: they are all measured against implicit-explicit Euler, and the further they stray, the harder the method is to control.

### Theorems

**(Theorem 3.2, unified unconditional stability and convergence)** Under the local Lipschitz condition on $\mathcal F$, sufficient regularity of $u$, the assumptions of Lemma 2.4, and the $k$-th order conditions, **if**

$$
\frac{\lambda_{\mathrm I}^{(k)}}{\sigma_{\mathrm E}^{(k)}}>\frac{\mu_0}{\varpi}
$$

and the step size $\tau$ (depending on $\lambda_{\mathrm I}^{(k)}/\sigma_{\mathrm F}^{(k)}$) is sufficiently small, then the $k$-step implicit-explicit multistep method is stable and convergent of order $\mathcal O(\tau^{k})$. The proof is complete mathematical induction on the bound $\|\tilde u^{\ell}\|_V\le1$ for the errors $\tilde u^j=U^j-u^j$, with truncation error $\|R_n^{(k)}\|_{\star}\le c_u\tau^{k}$.

**The shape of this criterion is worth pausing on**: the left-hand side $\lambda_{\mathrm I}^{(k)}/\sigma_{\mathrm E}^{(k)}$ depends only on the **method**, the right-hand side $\mu_0/\varpi$ only on the **model**. The criterion separates the two cleanly, and that separation is precisely why one statement covers every family at once.

**(Implicit-explicit controllability intensity)** Accordingly,

$$
\mathfrak I_{\mathrm{IE}}^{(k)}:=\frac{\lambda_{\mathrm I}^{(k)}}{\sigma_{\mathrm E}^{(k)}}
=\frac{\min_{\theta\in[0,2\pi)}\Re\bigl[b^{(k)}(\theta)/a^{(k)}(\theta)\bigr]}
{\max_{\theta\in[0,2\pi)}\bigl|c^{(k)}(\theta)/a^{(k)}(\theta)\bigr|}\ \le\ 1,
$$

with the optimal value $1$ attained by the implicit-explicit Euler scheme. **This is the paper's practical output: one number per method, by which schemes can be ranked and parameters chosen.**

**(New $\gamma$-parameterised SIEMS methods)** Defined through three characteristic polynomials:

$$
\tilde\varrho^{(k)}_{a,\mathrm S}(\zeta):=\sum_{j=1}^{k}\frac{f_{\mathrm S}^{(j)}(1)}{j!}(\zeta-1)^{j-1},
\qquad f_{\mathrm S}(z)=(\gamma z-\gamma+1)^{k-1}z\ln z,
$$

$$
\varrho^{(k)}_{b,\mathrm S}(\zeta):=\zeta(\gamma\zeta-\gamma+1)^{k-1},
\qquad
\varrho^{(k)}_{c,\mathrm S}(\zeta):=\zeta(\gamma\zeta-\gamma+1)^{k-1}-\gamma^{k-1}(\zeta-1)^{k}.
$$

All roots of $\varrho_{b,\mathrm S}^{(k)}$ lie in $|\zeta|<1$ when $\gamma>\frac12$; by Routh-Hurwitz, all roots of $\varrho^{(k)}_{c,\mathrm S}$ lie in $|\zeta|<1$ under the following thresholds for $k=2,\dots,8$:

| $k$       | $2$         | $3$                   | $4$        | $5$        | $6$ | $7$       | $8$      |
| --------- | ----------- | --------------------- | ---------- | ---------- | --- | --------- | -------- |
| $\gamma>$ | $-\tfrac12$ | $\tfrac{\sqrt2-1}{2}$ | $\tfrac38$ | $0.658691$ | $1$ | $1.37957$ | $1.7863$ |

SIEMS-2 coincides with WBDF2 and GBDF2. Unlike generalised BDF (zero-stable only up to $k=5$), **SIEMS-$k$ is zero-stable and satisfies the hypotheses of Lemma 2.4 for suitable $\gamma$ all the way to $k=8$**, so the theory yields unconditional stability at **eighth order**; the authors note they are not aware of any previously known unconditionally stable implicit-explicit multistep scheme of order above seven.

**(Comparative conclusion)** Five families are evaluated: $\alpha$-parameterised weighted BDF, $s$-parameterised modified BDF, $\beta$-parameterised generalised BDF, $\delta$-parameterised NIMEX, and the new $\gamma$-parameterised SIEMS. Based on the theoretical range of $\mathfrak I_{\mathrm{IE}}^{(k)}$, **generalised BDF-$k$ ($2\le k\le5$) and SIEMS-$k$ ($2\le k\le6$) have better adaptability to the nonlinear parabolic model than the other existing implicit-explicit multistep schemes.**

### Numerical experiments

**This version has no PDE numerical experiments**, only plots of the three factors; the conclusions state that a forthcoming report will illustrate the usage of the implicit-explicit controllability intensity.

What the paper does present is a **computational evaluation of the three factors**: plots of $\lambda_{\mathrm I}$, $\sigma_{\mathrm E}$ and $\sigma_{\mathrm F}$ as functions of the family parameter, covering GBDF4, GBDF5, WBDF5 and SIEMS4 through SIEMS8.

**This puts the paper in the same category as [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|papers 58 and 74]]**: every conclusion is a theorem, with no scheme-level measured data to compare against. The distinction matters — those factor curves test the computability of the lemmas and the relative merits of the families, not the convergence order of any scheme on any equation. The comparative statement about $2\le k\le5$ and $2\le k\le6$ above likewise rests on the theoretical range of $\mathfrak I_{\mathrm{IE}}^{(k)}$, not on numerical examples.

### Relation to the others

It is the natural sequel to [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|paper 58]]: both replace the Nevanlinna-Odeh multiplier technique with a DOC-based energy argument, and both apply Grenander-Szegő to a symmetrised Toeplitz form. The difference is that paper 58 handles **fully implicit BDF-$k$ on a uniform grid** while this handles **implicit-explicit multistep methods for genuinely nonlinear parabolic problems**, generalising the generating function to the one-sided ("semi") version so that **ratios** like $b/a$ and $c/a$ can be treated.

It is the multistep counterpart of paper 91: paper 91 uses a stage-index differentiation matrix $D_{\mathrm R}$ and its minimum eigenvalue $\lambda_{\min}$, and here the analogous object is the composite Toeplitz matrix $\widehat B_{L,k}$ and its minimum eigenvalue $\lambda_{\mathrm I}^{(k)}$. **The structure is identical in both: an eigenvalue of a symmetrised convolution or differentiation matrix dominating a norm bound on the explicit part.**

Deciding positive definiteness through generating functions is the **uniform-mesh analogue** of the algebraic criteria C1-C4 of [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|paper 74]]; both trace back to Toeplitz-Carathéodory and Grenander-Szegő, but paper 74 escapes to arbitrary nonuniform meshes, where no generating function exists at all.

**No step-ratio constant appears here**, since like paper 58 this is a uniform-grid paper. The relevant "small numbers" are instead $\mu_0/\varpi$ (the threshold supplied by the model) and $\mathfrak I_{\mathrm{IE}}^{(k)}\le1$ (the intensity supplied by the method).

## 78: linear relaxation with a regularised energy reformulation

Paper 78 takes a different route: **do not improve the time discretisation, rewrite the energy.** The paper names its method RRER, for relaxation with regularized energy reformulation.

### The idea

Energy-stable schemes for phase-field models are usually built with convex splitting (limited in accuracy and often nonlinear), exponential time differencing, stabilisation (typically first order), invariant energy quadratisation (IEQ), the scalar auxiliary variable method (SAV), Lagrange multipliers or supplementary variables. What IEQ and SAV have in common is rewriting the nonlinear energy as a **quadratic** form in an auxiliary variable, so each step needs only a linear solve. IEQ generalises the method of Lagrange multipliers and typically produces a **coupled** system with **time-dependent** coefficients; SAV retains its advantages but produces a **decoupled** system with **constant** coefficients.

They also share one move, and that is what this paper targets: **the evolution equation for the auxiliary variable is obtained by differentiating that variable in time.** That step introduces its own truncation error and makes the discrete system less faithful to the original equation. The paper's own phrasing is that it does not need to take time derivatives of the auxiliary variables.

**Its answer is to let the auxiliary variable be defined by an algebraic relation only, never differentiated in time.** That alone is not enough: if $q$ and $\phi$ lived on the same level, substituting the algebraic relation back would still make the scheme nonlinear. The second ingredient is Jiang et al.'s **relaxation** idea — put $q$ on **half-integer** levels. Then the algebraic relation defining $q^{n+\frac12}$ has only the **known** $\phi^n$ on its right, while $q^{n+\frac12}$ appears in the evolution equation only multiplied by $\frac{\phi^{n+1}+\phi^n}{2}$, so the whole scheme is linear in the unknowns and still second order. **The staggered grid is the trick that lets "algebraic definition" and "second-order accuracy" hold at the same time.**

### Setting

On $\Omega\subset\mathbb R^d$ with $d=2,3$ and periodic boundary conditions, write $(f,g)=\int_\Omega fg\,\mathrm d\mathbf x$ and $\|f\|=\sqrt{(f,f)}$. Start from a simplified free energy and its gradient flow, with $\mathcal G\ge0$ a semi-positive-definite mobility operator and $\mathcal L$ linear:

$$
E(\phi)=\tfrac12(\mathcal L\phi,\phi)+\bigl(F(\phi),1\bigr),
\qquad
\frac{\partial\phi}{\partial t}=-\mathcal G\bigl(\mathcal L\phi+F'(\phi)\bigr),
$$

$$
\frac{\mathrm d}{\mathrm dt}E(\phi)
=\Bigl(\frac{\delta E}{\delta\phi},\frac{\partial\phi}{\partial t}\Bigr)
=-\bigl(\mathcal L\phi+F'(\phi),\ \mathcal G(\mathcal L\phi+F'(\phi))\bigr)\le0 .
$$

**First test bed: the molecular beam epitaxy model with slope selection.**

$$
E(\phi)=\int_\Omega\Bigl(\frac{\epsilon^2}{2}(\Delta\phi)^2
+\frac14\bigl(|\nabla\phi|^2-1\bigr)^2\Bigr)\mathrm d\mathbf x,
\qquad \mathcal G=I,
$$

$$
\phi_t=-\epsilon^2\Delta^2\phi+\nabla\cdot\bigl((|\nabla\phi|^2-1)\nabla\phi\bigr),
$$

with mass conservation $\frac{\mathrm d}{\mathrm dt}\int_\Omega\phi\,\mathrm d\mathbf x=0$ under $\partial\phi/\partial\mathbf n|_{\partial\Omega}=0$ and $\partial\Delta\phi/\partial\mathbf n|_{\partial\Omega}=0$.

> [!note] A **different** MBE model from the one in paper 52
> This paper treats the model **with slope selection**, whose free energy carries the double well $\frac14(|\nabla\phi|^2-1)^2$; [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf|paper 52]] treats the model **without slope selection**, whose free energy carries the logarithmic term $-\frac12\ln(1+|\nabla\phi|^2)$. **The case without slope selection is the harder one for energy stability**, because there the nonlinearity is not polynomially bounded and the free energy has no lower bound. Paper 52's conclusions say explicitly that its technique does not apply to the model with slope selection, which this paper handles by entirely different means. The two must not be conflated.

**Second test bed: the phase-field crystal (PFC) model.**

$$
E(\phi)=\int_\Omega\Bigl(\tfrac12\phi(a_0+\Delta)^2\phi+\tfrac14\phi^4
-\tfrac{b_0}{2}\phi^2\Bigr)\mathrm d\mathbf x,
\qquad
\mathcal G=-\lambda\Delta,\quad 0<b_0<a_0,\ b_0\ll1,\ \lambda>0,
$$

$$
\phi_t=\lambda\Delta\mu,\qquad \mu=(a_0+\Delta)^2\phi+\phi^3-b_0\phi .
$$

**Step-ratio restriction: none** (the step $\delta t$ is uniform and the stability unconditional). **Maximum-principle argument: none. DOC/DCC kernels: none** — this paper is unrelated to the Liao-Tang-Zhou convolution-kernel programme.

### Derivation

**Step one: the regularised auxiliary variable.** Choose $C_0$ so that $F(\phi)+C_0\ge0$ and introduce a **stabilisation parameter** $\gamma$:

$$
q=\sqrt{4\bigl(F(\phi)+C_0\bigr)}-\gamma
\qquad\Longrightarrow\qquad
F(\phi)=\tfrac14q^2+\tfrac12\gamma q+\tfrac{\gamma^2}{4}-C_0,
\qquad
F'(\phi)=\tfrac12qq'+\tfrac{\gamma}{2}q',\ \ q'=\tfrac{\mathrm dq}{\mathrm d\phi} .
$$

The equivalent system and its energy are

$$
\frac{\partial\phi}{\partial t}=-\mathcal G\mu,
\quad
\mu=\mathcal L\phi+\tfrac12qq'+\tfrac{\gamma}{2}q',
\quad
q=\sqrt{4\bigl(F(\phi)+C_0\bigr)}-\gamma,
$$

$$
\widehat E(\phi,q)=\tfrac12(\mathcal L\phi,\phi)
+\Bigl(\tfrac14q^2+\tfrac12\gamma q+\tfrac{\gamma^2}{4}-C_0,\,1\Bigr),
\qquad
\frac{\mathrm d}{\mathrm dt}\widehat E(\phi,q)=-(\mu,\mathcal G\mu)\le0 .
$$

Remark 2.1 puts the whole point plainly: this reformulation differs from the IEQ and SAV approaches because no time derivatives of the auxiliary variables are taken, and taking them would introduce truncation errors during the computation.

**Step two: the concrete reformulation for MBE.** With $q=|\nabla\phi|^2-1-\gamma$ for $\gamma>0$,

$$
\begin{cases}
\phi_t=-\epsilon^2\Delta g+\nabla\cdot(q\nabla\phi)+\gamma g,\\
g=\Delta\phi,\\
q=|\nabla\phi|^2-1-\gamma,
\end{cases}
$$

$$
\widehat E(\phi,q)=\frac{\epsilon^2}{2}\|\Delta\phi\|^2+\frac{\gamma}{2}\|\nabla\phi\|^2
+\frac12\bigl(q(|\nabla\phi|^2-1-\gamma),1\bigr)-\frac14\|q\|^2
-\frac{2\gamma+\gamma^2}{4}|\Omega| .
$$

> [!warning] The relation between $\widehat E$ and $E$ differs at the continuous and discrete levels
> **Continuous level: $\widehat E(\phi,q)$ is exactly equal to $E(\phi)$, not merely an approximation** (Remark 2.2; the paper derives the identity line by line). This is unlike the modified energies of SAV and IEQ, which agree with the original energy only up to the consistency of the auxiliary variable.
>
> **Discrete level: the exact equality is lost.** Because $q$ lives on the half-integer levels of the staggered grid, Remarks 2.6 and 2.11 note that the discrete modified energy is only a **second-order approximation** to the original energy, which the paper says is consistent with the modified energies of the IEQ and SAV approaches.
>
> Citations must keep the two levels apart. The whole tension in the section on why this paper is the exception hangs on this distinction.

**Step three: the staggered time grid makes the scheme linear (Algorithm 1).** Given $(\phi^n,q^{n-\frac12})$, compute $(\phi^{n+1},g^{n+\frac12},q^{n+\frac12})$ from

$$
\begin{aligned}
&\text{(a)}\ \ \frac{\phi^{n+1}-\phi^{n}}{\delta t}
=-\epsilon^2\Delta g^{n+\frac12}
+\nabla\cdot\Bigl(q^{n+\frac12}\nabla\frac{\phi^{n+1}+\phi^{n}}{2}\Bigr)
+\gamma g^{n+\frac12},\\
&\text{(b)}\ \ g^{n+\frac12}=\Delta\frac{\phi^{n+1}+\phi^{n}}{2},
\qquad
\text{(c)}\ \ \frac{q^{n+\frac12}+q^{n-\frac12}}{2}=|\nabla\phi^{n}|^2-1-\gamma .
\end{aligned}
$$

**What makes the scheme linear is (c): it is algebraic, and its right-hand side involves only the known $\phi^n$.** Together with the staggering ($q$ at half-integer levels, $\phi$ at integer levels), each step therefore solves a single linear algebraic system. The starting values are $q^{\frac12}=|\nabla\phi^0|^2-1-\gamma$ and $g^{\frac12}=\Delta\frac{\phi^1+\phi^0}{2}$, followed by (a). The scheme is second-order accurate in time (Remark 2.3).

**Step four: the phase-field crystal model (Algorithm 2).** With $q=\phi^2-b_0-\gamma$ for $\gamma>0$,

$$
\text{(a)}\ \phi_t=\lambda\Delta\mu,
\quad
\text{(b)}\ \mu=(a_0+\Delta)g+q\phi+\gamma\phi,
\quad
\text{(c)}\ g=(a_0+\Delta)\phi,
\quad
\text{(d)}\ q=\phi^2-b_0-\gamma,
$$

$$
\widehat E(\phi,q)=\tfrac12\|\Delta\phi\|^2-a_0\|\nabla\phi\|^2
+\tfrac12(\gamma+a_0^2)\|\phi\|^2-\tfrac14\|q\|^2
+\tfrac12\bigl(q(\phi^2-b_0-\gamma),1\bigr)-\tfrac{(\gamma+b_0)^2}{4}|\Omega| .
$$

Again $\widehat E\equiv E$ exactly at the continuous level (Remark 2.7). Algorithm 2 is the analogous staggered Crank-Nicolson scheme, second order (Remark 2.8) and mass conserving (Remark 2.12).

**Step five: coupled models.** Section 2.4 extends RRER to a **ternary phase-field model** (2.4.1) and a **phase-field model for grain growth** (2.4.2).

### Theorems

**(Theorem 2.4, mass conservation, MBE)** Algorithm 1 preserves the total mass:

$$
\int_\Omega\phi^{n+1}\,\mathrm d\mathbf x=\int_\Omega\phi^{n}\,\mathrm d\mathbf x .
$$

The proof integrates (a) over $\Omega$, substitutes the integrated (b) and applies Green's formula.

**(Theorem 2.5, unconditional energy stability, MBE)** For any $n\ge1$,

$$
\widehat E^{n+1}\bigl(\phi^{n+1},q^{n+\frac12}\bigr)
-\widehat E^{n}\bigl(\phi^{n},q^{n-\frac12}\bigr)
=-\delta t\,\bigl\|\mu^{n+\frac12}\bigr\|^2\ \le\ 0,
$$

$$
\mu^{n+\frac12}=\epsilon^2\Delta g^{n+\frac12}
-\nabla\cdot\Bigl(q^{n+\frac12}\nabla\frac{\phi^{n+1}+\phi^{n}}{2}\Bigr)
-\gamma g^{n+\frac12},
$$

$$
\widehat E^{n+1}=\frac{\epsilon^2}{2}\|\Delta\phi^{n+1}\|^2
+\frac{\gamma}{2}\|\nabla\phi^{n+1}\|^2
+\frac12\bigl(q^{n+\frac12}(|\nabla\phi^{n+1}|^2-1-\gamma),1\bigr)
-\frac14\|q^{n+\frac12}\|^2-\frac{2\gamma+\gamma^2}{4}|\Omega| .
$$

**Note this is an equality, not an inequality** — the scheme dissipates exactly $\delta t\|\mu^{n+\frac12}\|^2$. The proof pairs (a) with $\delta t\,\mu^{n+\frac12}$ and (c) with $q^{n+\frac12}-q^{n-\frac12}$, the latter giving the key identity

$$
\tfrac12\bigl(\|q^{n+\frac12}\|^2-\|q^{n-\frac12}\|^2\bigr)
=\bigl((q^{n+\frac12}-q^{n-\frac12})(|\nabla\phi^n|^2-1-\gamma),1\bigr).
$$

**There is no restriction on $\delta t$.**

**(Theorem 2.10, unconditional energy stability, PFC)** For any $n\ge1$,

$$
\widehat E^{n+1}\bigl(\phi^{n+1},q^{n+\frac12}\bigr)
-\widehat E^{n}\bigl(\phi^{n},q^{n-\frac12}\bigr)
=-\lambda\,\delta t\,\bigl\|\nabla\mu^{n+\frac12}\bigr\|^2\ \le\ 0,
$$

$$
\widehat E^{n+1}=\tfrac12\|\Delta\phi^{n+1}\|^2-a_0\|\nabla\phi^{n+1}\|^2
+\tfrac12(\gamma+a_0^2)\|\phi^{n+1}\|^2-\tfrac14\|q^{n+\frac12}\|^2
+\tfrac12\bigl(q^{n+\frac12}[(\phi^{n+1})^2-b_0-\gamma],1\bigr)
-\tfrac{(\gamma+b_0)^2}{4}|\Omega| .
$$

### Numerical experiments

Space is discretised by standard **finite elements** ($P_1$), implemented in **FreeFEM**, with periodic boundary conditions. Six groups of tests.

**Test 1 (convergence, §3.1).** $\gamma=2.0$, $\Omega=[0,1]^2$, $T=1$, manufactured solution $\phi(x,y,t)=e^{-t}\cos(\pi x)\cos(\pi y)$, for both MBE ($\epsilon=1$) and PFC ($a_0=1.0$, $b_0=0.01$, $\lambda=1.0$), with $\delta t=h$; compared against IEQ and exponential SAV (ESAV). The MBE numbers are:

| $\delta t$          | $\frac18$    | $\frac1{16}$ | $\frac1{32}$ | $\frac1{64}$ |
| ------------------- | ------------ | ------------ | ------------ | ------------ |
| RRER error          | $3.2923$e-03 | $8.2307$e-04 | $2.0577$e-04 | $5.1442$e-05 |
| RRER observed order | —            | $2.00$       | $2.00$       | $2.00$       |
| IEQ observed order  | —            | $1.94$       | $1.97$       | $1.99$       |

ESAV matches RRER. **All three methods are second order**; for PFC, RRER is more accurate than both IEQ and ESAV, and a CPU-time table shows RRER is the cheapest of the three.

**Test 2 (energy dissipation and mass conservation, §3.2).** Initial datum $\phi(x,y,0)=\cos\pi x\cos\pi y$; the original and modified energies are plotted for several $\delta t$ and $\epsilon$. Both models decay monotonically for every $\delta t$ tested, supporting the unconditional stability, and mass is conserved (MBE with $\gamma=20$, $\epsilon=1.0$; PFC with $\gamma=1.0$, $a_0=1.0$, $b_0=0.325$, $\delta t=0.1$).

**Test 3 (two-dimensional MBE coarsening, §3.3).** $\Omega=[0,2\pi]^2$, $\epsilon^2=0.1$, $\gamma=20$, $\phi(x,y,0)=0.1(\sin3x\sin2y+\sin5x\sin5y)$, $h=2\pi/128$, $\delta t=10^{-4}$; snapshots at $t=0,0.05,2.5,8,15,30$ reproduce the known phase diagram of MBE **with** slope selection.

**Test 4 (two-dimensional PFC, §3.4).** $\Omega=[0,100]^2$, $\phi(x,y,0)=\hat\phi_0+0.01\,\mathrm{rand}(x,y)$, $a_0=1.0$, $b_0=0.35$, $T=200$, $\gamma=2.0$, $h=100/128$, $\delta t=0.1$: a **stripe** pattern for $\hat\phi_0=0$ and a **triangle** pattern for $\hat\phi_0=0.2$, consistent with the literature phase diagram.

**Test 5 (PFC on curved surfaces, §3.5).** A sphere of radius $R=64$ (striped or hexagonal patterns at $t=50,100,200$ depending on $\hat\phi_0$); a ring torus of outer radius $50$ and inner radius $20$ with $h=1$, $\delta t=0.1$ (stripes versus hexagons); and a cube-surface case on $[0,100]^3$ with $h=100/64$, $\delta t=0.1$.

**Test 6 (three-dimensional PFC and a coupled system, §3.6).** Plus a ternary-system test with $m=3$, $L_1=L_2=L_3=1$, $\alpha=\beta=\gamma=1$, $k_1=k_2=k_3=2$ on $[0,1]^2$.

**What these tests establish and what they do not.** They establish that the scheme really is second order (with temporal accuracy no worse than IEQ or ESAV), that the energy decreases monotonically and mass is conserved at every step size tested, that the morphologies agree with the phase diagrams in the literature, and that the CPU cost is the lowest of the three methods. They do not establish **quantitative accuracy at large steps** — unconditional stability means the computation will not blow up, not that it is accurate with a large $\delta t$, and Test 2 reports only monotonicity of the energy, not the difference between solutions at different $\delta t$. Also, the stabilisation parameter $\gamma$ takes very different values across the tests ($1.0$, $2.0$, $20$) and the paper gives no criterion for choosing it.

### Why this paper is the exception in this topic

This is worth stating plainly, because it bears on a tension running through the whole topic. Paper 78 has no Liao co-authorship, no DOC or DCC kernels, no variable steps and no step-ratio analysis; it belongs to the IEQ/SAV/relaxation branch rather than the convolution-kernel branch, with Tao Zhou the only point of contact.

More importantly, **the type of its energy statement differs**: paper 78 proves dissipation of a **modified** energy $\widehat E(\phi,q)$, whereas [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|papers 40 and 57]] together with papers 48, 52 and 91 work hard precisely to **avoid** that, establishing dissipation of the original or variational energy instead. The tension is explicit in the literature: **paper 91's introduction criticises SAV-based high-order schemes for establishing stability only with respect to a modified energy involving the auxiliary variable.**

What mitigates it for paper 78 is that $\widehat E\equiv E$ holds **exactly at the continuous level**, so its modified energy is not a new object but another way of writing the same one. But the warning above must be carried along: at the **discrete** level, because $q$ lives on the staggered levels, $\widehat E$ degrades to a second-order approximation of the original energy, which the paper itself says is consistent with the IEQ and SAV situation. **So this is not an argument that dissolves the tension; it is an argument that shrinks it to second order.**

Finally, this paper and paper 91 both treat phase-field-crystal and Cahn-Hilliard-type fourth-order dynamics, but with entirely different stabilisation philosophies: auxiliary-variable quadratisation here, explicit stabilisation plus Runge-Kutta there.

## The three approaches side by side

| No. | Treated explicitly                    | Object carrying the difficulty                      | How the energy argument is kept                                    | Energy that decays                                                    | Numerical experiments                                        |
| --- | ------------------------------------- | --------------------------------------------------- | ------------------------------------------------------------------ | --------------------------------------------------------------------- | ------------------------------------------------------------ |
| 78  | the nonlinear potential (quadratised) | none — the energy changes, not the discretisation   | reformulate the energy; a staggered grid makes the scheme linear   | **modified** energy $\widehat E$ ($\equiv E$ at the continuous level) | six groups, second order, phase diagrams reproduced          |
| 91  | the nonlinear bulk $f_\kappa$         | the differentiation matrix $D(z)$ (stage-index DOC) | impose $D_{\mathrm{EI}}=\mathbf 0$ so the rate is mesh-independent | the **original** energy $E$                                           | five groups, with a three-orders-of-magnitude speed-up table |
| 104 | the nonlinear term $\mathcal F$       | the composition of three kernel sequences           | extrema of three semi-generating functions on the unit circle      | $V$-norm stability in an abstract framework                           | none for PDEs; only the three factor curves                  |

**How to read this table**: none of the three tries to prove directly that the explicit term contributes with a definite sign; each substitutes a controllable object — paper 78 the energy, paper 91 a structural condition on the Butcher tableau, paper 104 the form of the criterion. What is genuinely comparable is the "energy that decays" column: **papers 91 and 104 stay on the side of the original energy or norm, paper 78 moves to the modified-energy side**, and that line is the one real methodological disagreement inside this topic.

## Sources

- J. Zhang, X. Guo, M. Jiang, T. Zhou, and J. Zhao, [_Linear relaxation method with regularized energy reformulation for phase field models_](https://doi.org/10.1016/j.jcp.2024.113225), J. Comput. Phys. 515 (2024), 113225.
- H.-l. Liao, T. Tang, X. Wang, and T. Zhou, [_A class of refined implicit-explicit Runge-Kutta methods with robust time adaptability and unconditional convergence for the Cahn-Hilliard model_](https://doi.org/10.1090/mcom/4090), Math. Comput. 95(359) (2026), pp. 1293-1325 (preprint [arXiv:2412.07321](https://arxiv.org/abs/2412.07321)).
- H.-l. Liao, C. Quan, T. Tang, and T. Zhou, _A semi-generating function approach to the stability of implicit-explicit multistep methods for nonlinear parabolic equations_, [arXiv:2605.05619](https://arxiv.org/abs/2605.05619), submitted to SIAM J. Numer. Anal.
