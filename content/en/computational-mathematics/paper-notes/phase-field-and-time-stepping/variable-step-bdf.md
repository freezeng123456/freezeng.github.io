---
title: Variable-Step BDF and Convolution Kernels
description: Papers 48, 52, 58, 67, 69 and 74 - turning sign-indefinite multistep kernels into non-negative decreasing ones is what makes an energy argument possible
lang: en
translation: computational-mathematics/paper-notes/phase-field-and-time-stepping/variable-step-bdf
tags:
  - paper-notes
  - phase-field
  - discrete-energy
---

## One mechanism running through the line of work

Before any individual paper it pays to state the shared difficulty and the shared remedy, because all six papers are variations on it.

**The difficulty.** On a uniform grid, BDF2 is a Toeplitz convolution: the kernel $b_j$ depends only on the lag index $j$, so the whole classical apparatus — the generating function $\hat b(z)=\sum_jb_jz^j$, the Grenander-Szegő theorem, the Toeplitz-Carathéodory criterion — is available. Once the steps vary, the kernel becomes $b^{(n)}_j$, carrying an extra superscript that records which level we are currently on. It is no longer Toeplitz, the operator is no longer self-adjoint, and the generating-function machinery collapses.

**The difficulty is specifically about sign.** The two variable-step BDF2 kernels are

$$
b^{(n)}_0=\frac{1+2r_n}{\tau_n(1+r_n)}>0,
\qquad
b^{(n)}_1=-\frac{r_n^2}{\tau_n(1+r_n)}<0 .
$$

The second one is **negative**, and that single minus sign blocks two roads at once:

- **The energy argument** needs the quadratic form $\sum_kw_k\sum_jb^{(k)}_{k-j}w_j$ to be positive definite. With a negative entry in the kernel, positivity is no longer automatic; it holds only while the negative term is dominated by the positive one, and that is exactly where the threshold $r_k<(3+\sqrt{17})/2$ comes from.
- **The maximum principle** needs the scheme written as "a positive-definite operator on the current level" $=$ "a **non-negative** combination of past levels". With non-negative weights, the induction hypothesis $\|u^k\|_\infty\le1$ bounds the right-hand side, because a weighted average with non-negative weights cannot exceed the largest quantity being averaged. **The moment a weight goes negative, the argument dies outright**: a negative coefficient turns an upper bound into no bound at all.

**The remedies.** There are three:

1. **Kernel recombination (KRC)** — change variables to $\bar v^k=v^k-\eta v^{k-1}$ so that the new kernel $d^{(n)}_j$ is **non-negative and decreasing**. This is the route paper 48 takes for the maximum principle, and it is available exactly when $r_k<1+\sqrt2$.
2. **Discrete orthogonal convolution (DOC) kernels** — construct $\theta$ with $\Theta B=I$ and peel the multistep operator off the scheme, leaving a plain first difference behind. This is the route of papers 52, 58, 67 and 74, and it reduces a multistep stability question to a question about positivity of a quadratic form.
3. **Discrete complementary convolution (DCC) kernels** — construct $Q_d$ whose row sums are identically $1$, and use it to build a discrete Grönwall inequality. Paper 48 uses it; [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|paper 43]] used it first on the fractional side.

In one line: **DOC inverts ($\equiv\delta_{nk}$, producing energy and $L^2$ estimates), DCC complements ($\equiv1$, producing Grönwall and $\ell^\infty$ estimates), KRC changes basis (making the kernel non-negative and decreasing, producing the maximum principle).**

## 48: two step-ratio constants, two different properties

### The idea

The Allen-Cahn equation satisfies both energy dissipation and a maximum bound principle at the continuous level, and its solutions span two wildly separated time scales — fast early evolution, then extremely slow coarsening — so variable steps are a practical necessity. At the time of writing, though, the numerical analysis of variable-step schemes was, in the authors' words, far from complete even for linear and semilinear parabolic equations. Two specific gaps: nobody had proved a **discrete maximum principle** for a second-order variable-step scheme for Allen-Cahn, and the sharpness of the step-ratio conditions in the existing nonuniform BDF2 energy results was unclear.

The paper's key judgement is that **energy stability and the maximum principle need two different properties of the kernel, so they come with two different thresholds, and the two cannot be swapped.**

The energy side needs **positive definiteness**. Testing the scheme with $\nabla_\tau u^n$ produces $\sum_k\nabla_\tau u^k\,D_2u^k$ on the left, a quadratic form whose coefficients are the BDF2 kernels. The kernel has a negative entry, so the question is when the negative term is dominated. The answer is $r_k<(3+\sqrt{17})/2\approx3.561$.

The maximum-principle side needs **non-negativity plus monotonicity**, which is strictly more, for the reason given above: the $\ell^\infty$ induction only works with non-negative weights. BDF2's kernels do not qualify, so the paper performs a **change of variable**: set $\bar v^k=v^k-\eta v^{k-1}$. This trades a sign-indefinite two-term kernel for a geometrically decaying non-negative family $d^{(n)}_j$, at the cost of confining $\eta$ to an interval — and that interval is non-empty precisely when $r_k<1+\sqrt2$.

**That is why $1+\sqrt2$ and $3.561$ are not two versions of one conclusion but answers to two different questions.** Better still, $1+\sqrt2$ is exactly Grigorieff's (1983) classical **zero-stability** bound for variable-step BDF2 applied to ODEs, and this paper obtains the discrete maximum principle and second-order max-norm convergence right at that classical limit — whereas earlier $L^2$ analyses all required something stricter than $1+\sqrt2$.

### Setting

On $\Omega=(0,L)^2$ with periodic boundary conditions,

$$
\partial_{t}u=\varepsilon^{2}\Delta u-f(u),\quad f(u)=u^{3}-u,
\qquad u(\mathbf x,0)=u_{0}(\mathbf x),
$$

the $L^2$ gradient flow of the Ginzburg-Landau free energy

$$
E[u](t):=\int_{\Omega}\Bigl(\tfrac12\varepsilon^{2}|\nabla u|^{2}+F[u]\Bigr)\mathrm d\mathbf x,
\qquad F[u]=\tfrac14(1-u^{2})^{2},
$$

so that $\mathrm dE/\mathrm dt\le0$, and $|u(\mathbf x,0)|\le1$ implies $|u(\mathbf x,t)|\le1$.

The grid and step ratios (standard notation across this whole literature):

$$
0=t_0<t_1<\cdots<t_N=T,
\qquad \tau_k:=t_k-t_{k-1},
\qquad \tau:=\max_{1\le k\le N}\tau_k,
$$

$$
r_{k}:=\frac{\tau_{k}}{\tau_{k-1}}\ (2\le k\le N),
\qquad r_1\equiv0\ \text{(by convention)} .
$$

Writing $\Pi_{n,2}v$ for the quadratic interpolant of $v$ at $t_{n-2},t_{n-1},t_n$, the variable-step BDF2 formula is

$$
D_{2}v^{n}:=(\Pi_{n,2}v)'(t_{n})
=\frac{1+2r_{n}}{\tau_{n}(1+r_{n})}\nabla_{\tau}v^{n}
-\frac{r_{n}^{2}}{\tau_{n}(1+r_{n})}\nabla_{\tau}v^{n-1},
\qquad n\ge2,
$$

with $D_2v^1:=D_1v^1=\nabla_\tau v^1/\tau_1$ on the first level. At $r_n=1$ it degenerates to the classical $(3v^n-4v^{n-1}+v^{n-2})/(2\tau)$. As a discrete convolution $D_{2}v^{n}=\sum_{k=1}^{n}b^{(n)}_{n-k}\nabla_{\tau}v^{k}$,

$$
b^{(1)}_0:=\frac{1}{\tau_1};
\qquad
b^{(n)}_0:=\frac{1+2r_n}{\tau_n(1+r_n)},
\qquad
b^{(n)}_1:=-\frac{r_n^2}{\tau_n(1+r_n)},
\qquad
b^{(n)}_j:=0\ (2\le j\le n).
$$

The fully discrete scheme ($\Lambda_h$ the periodic central-difference Laplacian) is **fully implicit**, with neither stabilisation nor convex splitting:

$$
D_{2}u^{n}=\varepsilon^{2}\Lambda_{h}u^{n}-f(u^{n}),
\quad n\ge1,
\qquad f(u^n):=(u^n)^{.3}-u^n .
$$

The paper works with two named conditions (the wording below follows the original):

- **S1**: $0<r_{k}<\dfrac{3+\sqrt{17}}{2}\approx3.561$ for $2\le k\le N$, used for **energy stability**. The constant does not originate here: it is Lemma 2.1 of Liao and Zhang (_Math. Comp._ 90 (2021) 1207-1226), the condition for positive definiteness of the BDF2 kernels, and papers 48 and 52 both import it as S1. Paper 52's own verdict on it is blunt — it calls it an artificial constant that is due to the condition S1.
- **S0**: $0<r_{k}<1+\sqrt{2}\approx2.414$ for $2\le k\le N$, used for the **discrete maximum principle** and **max-norm** convergence; the paper notes explicitly that this coincides with Grigorieff's (1983) zero-stability condition for ODEs.

> [!warning] Two constants, two different properties
> **S1 ($3.561$) buys energy stability; S0 ($1+\sqrt2$) buys the maximum principle and max-norm convergence.** Energy stability comes from positivity of a quadratic form; the maximum principle comes from non-negativity and monotonicity after kernel recombination. They are not interchangeable, and neither is simply "better" than the other. Conflating them when citing is one of the most common errors in this literature.

### Derivation

**Step one: the modified energy and the exact origin of $3.561$.** What decays is not the plain discrete energy but a modified one,

$$
\widehat{E}[u^{k}]:=E[u^{k}]+\frac{r_{k+1}\tau_{k}}{2(1+r_{k+1})}
\sum_{i=1}^{M}\bigl(\partial_{\tau}u_{i}^{k}\bigr)^{2},
\quad k\ge1,
\qquad
\widehat{E}[u^{0}]:=E[u^{0}],
$$

$$
E[u^{k}]:=-\frac{\varepsilon^{2}}{2}(u^{k})^{T}\Lambda_{h}u^{k}
+\frac{1}{4}\sum_{i=1}^{M}\bigl(1-(u_{i}^{k})^{2}\bigr)^{2},
\quad k\ge0 .
$$

The correction is $O(\tau)$, so $\widehat E\to E$ as $\tau\to0$. Note that it uses $r_{k+1}$, the **next** step ratio; that look-ahead is exactly what makes the telescoping sum work.

From $2a(a-b)=a^2-b^2+(a-b)^2$ and the kernel definitions one gets the pointwise inequality

$$
D_{2}u_{i}^{n}\,(\nabla_\tau u_{i}^{n})
\ \ge\
\frac{r_{n+1}\tau_{n}}{2(1+r_{n+1})}(\partial_{\tau}u_{i}^{n})^{2}
-\frac{r_{n}\tau_{n-1}}{2(1+r_{n})}(\partial_{\tau}u_{i}^{n-1})^{2}
+\Bigl(\frac{2+4r_{n}-r_{n}^{2}}{1+r_{n}}-\frac{r_{n+1}}{1+r_{n+1}}\Bigr)
\frac{\tau_{n}}{2}(\partial_{\tau}u_{i}^{n})^{2}.
$$

**That single line is the entire energy argument**: the first two terms telescope into $\widehat E$, and the last must stay non-negative after absorbing the nonlinear contribution $-\frac{\tau_n^2}{2}\sum_i(\partial_\tau u_i^n)^2$. Writing out positivity of the bracket in the worst case gives

$$
\frac{r_{k+1}}{1+r_{k+1}}<\frac{r_{s}}{1+r_{s}}=\frac{\sqrt{17}-1}{4}\approx0.78,
\qquad
r_{s}=\frac{3+\sqrt{17}}{2}\ \text{the positive root of}\ 2+3r-r^{2}=0 .
$$

Numerically $r_s=3.5616\ldots$ and $r_s/(1+r_s)=3.5616/4.5616=0.7808=(\sqrt{17}-1)/4$, which is self-consistent.

The paper also derives the matching step-size bounds. With $h(x):=\dfrac{2+4x-x^{2}}{1+x}$,

$$
h'(x)=\frac{x+1+\sqrt3}{(1+x)^{2}}\bigl(\sqrt3-1-x\bigr),
$$

so $h$ increases on $(0,\sqrt3-1]$ and decreases afterwards, with $h(0)=2$ and $h(\sqrt2+1)=1+\frac{\sqrt2}{2}$. Three regimes follow:

1. $0<r_k\le\sqrt3-1$: it suffices that $\tau_k\le\min\{1,\frac{9-\sqrt{17}}{4}\}=1$;
2. $\sqrt3-1<r_k\le\sqrt2+1$: it suffices that $\tau_k\le1+\frac{\sqrt2}{2}-\frac{\sqrt{17}-1}{4}\approx0.93$;
3. $\sqrt2+1<r_k<r_s$: one must control the **next** ratio as well — for instance $\tau_k\le\frac12$ suffices when $r_{k+1}\le\frac{2h(r_s)-1}{3-2h(r_s)}\approx0.39$.

The third regime is worth noticing: **the closer the step ratio gets to $r_s$, the tighter the constraint on the following step**, which is the look-ahead term of the modified energy making itself felt.

> [!warning] About $4.8645$, so often quoted alongside these
> The positive-definiteness threshold for the BDF2 kernels was later improved to $r_*\approx4.864$ (the positive root of $1+2r-r^{3/2}=0$, exactly $4.864536512317583$), but that improvement **belongs to none of the papers in this topic**. It is due to Liao, Ji, Wang and Zhang (_J. Sci. Comput._ 92 (2022) 52), where it is printed as $4.864$; the more widely circulated form $4.8645$ comes from Zhang and Zhao (_J. of Math._ (PRC) 41(6) (2021) 471-488). Papers 69 and 91 both cite Liao-Ji-Wang-Zhang, which is probably how the constant came to be associated with this group of work. The full disentangling is on [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/index|the topic index]]. One more point: $4.8645$ improves $3.561$ (the same positive-definiteness lemma) and does **not** replace $1+\sqrt2$, which governs something else entirely.

**Step two: kernel recombination and the origin of $1+\sqrt2$.** This is the paper's signature technique, and the derivation is remarkably clean. Introduce a parameter $\eta\in\mathbb R$ and

$$
\bar{v}^{0}:=v^{0},
\qquad
\bar{v}^{k}:=v^{k}-\eta\,v^{k-1}\ (k\ge1),
\qquad
v^{k}=\sum_{\ell=0}^{k}\eta^{k-\ell}\bar{v}^{\ell} .
$$

Substituting and exchanging the order of summation gives the **recombined** BDF2 formula

$$
D_{2}v^{n}\equiv\sum_{j=1}^{n}d^{(n)}_{n-j}\nabla_{\tau}\bar{v}^{j}+d^{(n)}_{n}\bar{v}^{0},
\qquad
d^{(n)}_{n-j}:=\sum_{k=j}^{n}b^{(n)}_{n-k}\eta^{k-j},
\qquad
d^{(n)}_{n}:=\eta\,d^{(n)}_{n-1},
$$

whose closed form is extremely simple:

$$
d^{(n)}_{0}=b_{0}^{(n)},
\qquad
d^{(n)}_{j}=\eta^{j-1}\bigl(b_{0}^{(n)}\eta+b_{1}^{(n)}\bigr)\quad(1\le j\le n).
$$

The goal is to make the new kernel **non-negative and decreasing**, $d^{(n)}_{0}\ge d^{(n)}_{1}\ge\cdots\ge d^{(n)}_{n}\ge0$. Substituting the BDF2 kernels, that is exactly

$$
\frac{r_{k}^{2}}{1+2r_{k}}\ \le\ \eta\ <\ 1\qquad(k\ge2).
$$

**Such an $\eta$ exists if and only if $\dfrac{r_k^{2}}{1+2r_k}<1$, i.e. $r_{k}^{2}-2r_{k}-1<0$, i.e. $r_{k}<1+\sqrt{2}$.** That is condition S0, and it is Grigorieff's classical zero-stability bound.

The optimal choice of $\eta$ is also supplied. With

$$
K(\eta):=\frac{1-\eta}{\eta^{2}}\cdot\frac{(1+2r_{n})\eta-r_{n}^{2}}{1+r_{n}},
\qquad
K'(\eta)=\frac{1+r_{n}}{\eta^{3}}\Bigl(\frac{2r_{n}^{2}}{(1+r_{n})^{2}}-\eta\Bigr),
$$

$K$ is maximised as $\eta\to\frac{2r_n^2}{(1+r_n)^2}$, and the paper takes

$$
\eta:=\frac{2r_{s}^{2}}{(1+r_{s})^{2}},
\qquad r_{s}\in[1,\,1+\sqrt2)\ \text{the maximum step ratio} .
$$

Two worked cases: $r_s=1$ gives $\eta=\frac12$ and $\tau\le\frac{1}{2(S_n+4\varepsilon^{2}h^{-2})}$; $r_s=2$ gives $\eta=\frac89$ and $\tau_n\le\frac{1}{48}\cdot\frac{1}{S_n+4\varepsilon^{2}h^{-2}}$. **The closer $\eta$ gets to $1$, the smaller $1-\eta$ and hence the smaller the admissible step** — that is the algebraic reason behind "larger step ratios force tighter step bounds".

**Step three: the maximum-principle argument.** In the recombined variables the scheme becomes

$$
\bigl(d^{(n)}_{0}+S_{n}-\varepsilon^{2}\Lambda_{h}\bigr)\bar{u}^{n}
=\sum_{j=0}^{n-1}Q_{n-j}^{(n)}\bar{u}^{j}+(S_{n}+1)u^{n}-(u^{n})^{.3},
$$

$$
Q_j^{(n)}:=\bigl(d^{(n)}_{j-1}-d^{(n)}_{j}-S_{n}\eta^{j}\bigr)I+\eta^{j}\varepsilon^{2}\Lambda_{h},
$$

with $S_n$ a second free parameter. Under the step restriction

$$
\tau_{n}\le\frac{(1+2r_{n})\eta-r_{n}^{2}}{\eta^{2}(1+r_{n})}
\cdot\frac{1-\eta}{S_{n}+4\varepsilon^{2}h^{-2}}
$$

every entry of $Q_j^{(n)}$ is non-negative and $\|Q_j^{(n)}\|_{\infty}\le d^{(n)}_{j-1}-d^{(n)}_{j}-S_{n}\eta^{j}$. **"Every entry non-negative" is the whole point of the recombination**: without it the induction below cannot take a single step.

Combined with an $\ell^\infty$ cubic lemma (if $B=(b_{ij})$ satisfies $b_{ii}=-\max_i\sum_{j\ne i}|b_{ij}|$ and $A=aI-B$ with $a>0$, then $\|AV+cV^{3}\|_{\infty}\ge a\|V\|_{\infty}+c\|V\|_{\infty}^{3}$ for $c>0$), complete induction gives the **auxiliary claim** that carries the whole paper,

$$
\|\bar{u}^{k}\|_{\infty}\le 1-\eta\qquad(1\le k\le N),\ \text{provided}\ \|\bar u^0\|_\infty\le1,
$$

and then $u^{k}=\sum_{\ell}\eta^{k-\ell}\bar u^{\ell}$ together with $\sum_{j\ge0}\eta^j(1-\eta)=1$ yields $\|u^{k}\|_{\infty}\le1$. **The factor $1-\eta$ is no coincidence**: it is the reciprocal of the geometric series $\sum\eta^j$, so the bound in the recombined variable has to shrink by exactly that factor for the reconstruction to land on $1$. Two elementary facts are used along the way: for $c\ge1$, $g_c(z):=(c-1)z+z^{3}-c$ increases on $z>0$; and $|z|\le1$ implies $|(c+1)z-z^{3}|\le c$.

**Step four: DCC kernels and a new Grönwall inequality.** These are built on the **recombined** kernels $d^{(n)}_j$:

$$
(Q_{d})_{0}^{(n)}:=\frac{1}{d_{0}^{(n)}},
\qquad
(Q_{d})_{n-j}^{(n)}:=\sum_{k=j+1}^{n}\frac{d_{k-j-1}^{(k)}-d_{k-j}^{(k)}}{d_{0}^{(j)}}(Q_{d})_{n-k}^{(n)},
\quad 1\le j\le n-1,
$$

$$
\sum_{j=k}^{n}(Q_{d})^{(n)}_{n-j}\,d_{j-k}^{(j)}\equiv 1\quad(1\le k\le n),
\qquad
0<(Q_{d})^{(n)}_{n-j}\le\frac{1}{d^{(j)}_{0}} .
$$

The practical bounds are $(Q_{d})^{(n)}_{n-1}\le1/b^{(1)}_{0}=\tau_1$ and $(Q_{d})^{(n)}_{n-j}\le1/b^{(j)}_{0}=\frac{1+r_j}{1+2r_j}\tau_j\le\tau_j$ for $2\le j\le n$. Kernels of this kind are due to Liao-Li-Zhang (fractional Caputo) and Liao-McLean-Zhang (general kernels). **Note that non-negativity of the DCC kernels also depends on $d^{(n)}_j$ being non-negative and decreasing, so they too are downstream of the recombination.**

### Theorems

**(Unique solvability)** If $\tau_{n}<\dfrac{1+2r_{n}}{1+r_{n}}$ for $n\ge1$, the scheme is uniquely solvable; the paper notes that $\tau_n<1$ suffices. The proof first gets $b_0^{(n)}>1$, so $G_h:=b_0^{(n)}-1-\varepsilon^{2}\Lambda_h$ is positive definite and $u^n$ is the unique minimiser of the strictly convex functional $\frac12w^TG_hw+\frac14\sum_kw_k^4-w^Tg(u^{n-1})$ with $g(u^{n-1}):=b_0^{(n)}u^{n-1}-b_1^{(n)}\nabla_\tau u^{n-1}$.

**(Energy stability, the $3.561$ conclusion)** Assume **S1** and

$$
\tau_{k}\le\min\Bigl\{\frac{1+2r_{k}}{1+r_{k}},\
\frac{2+4r_{k}-r_{k}^{2}}{1+r_{k}}-\frac{r_{k+1}}{1+r_{k+1}}\Bigr\},
\qquad k\ge1 .
$$

Then the variable-step BDF2 solution obeys the discrete energy dissipation law

$$
\widehat{E}[u^{k}]\le\widehat{E}[u^{k-1}],\qquad k\ge1 .
$$

**What decays is the modified energy $\widehat E$, not $E$.** The proof also uses the algebraic identity

$$
4(a^{3}-a)(a-b)+2(1-a^{2})(a-b)^{2}
=(1-a^{2})^{2}-(1-b^{2})^{2}+(a^{2}-b^{2})^{2}.
$$

**(Discrete maximum principle, the $1+\sqrt2$ conclusion)** Assume **S0**, take $\eta:=\frac{2r_s^2}{(1+r_s)^2}$ and

$$
\tau_{n}\le\frac{(1+2r_{n})\eta-r_{n}^{2}}{\eta^{2}(1+r_{n})}
\cdot\frac{1-\eta}{2+4\varepsilon^{2}h^{-2}},\qquad n\ge1
$$

(that is, the step condition above with $S_n=2$). Then $\|u^{0}\|_{\infty}\le1$ implies $\|u^{k}\|_{\infty}\le1$ for $1\le k\le N$. **This step restriction is real, coupling $\tau$ to $h$ and $\varepsilon$ through $4\varepsilon^2h^{-2}$**; the paper calls it a practical time-step constraint.

**(A discrete Grönwall inequality via complementary kernels)** Let $\kappa>0$, $\lambda\in(0,1)$, and let non-negative sequences $\{g^k\},\{w^k\}$ satisfy

$$
\sum_{k=1}^{n}d^{(n)}_{n-k}\nabla_{\tau}w^{k}
\le\kappa\sum_{k=1}^{n}\lambda^{n-k}w^{k}+g^{n},
\qquad 1\le n\le N .
$$

If $b^{(n)}_{0}\ge2\kappa$, then

$$
w^{n}\le2\exp\Bigl(\frac{2\kappa t_{n}}{1-\lambda}\Bigr)
\Bigl(w^{0}+\sum_{j=1}^{n}\frac{g^{j}}{b^{(j)}_{0}}\Bigr),
\qquad 1\le n\le N .
$$

The proof is the one move — multiply by $(Q_d)^{(n)}_{n-j}$, exchange the summations, apply the complementarity identity:

$$
\sum_{j=1}^n(Q_d)^{(n)}_{n-j}\sum_{k=1}^jd^{(j)}_{j-k}\nabla_\tau w^k
=\sum_{k=1}^n\nabla_\tau w^k\sum_{j=k}^n(Q_d)^{(n)}_{n-j}d^{(j)}_{j-k}
=w^n-w^0 .
$$

**(Second-order max-norm convergence)** Let $u_0$ be smooth with $|u_0|\le1$ and let the exact solution be sufficiently smooth. Under **S0** and the same step restriction as in the maximum-principle theorem,

$$
\|u(\mathbf x_{h},t_{n})-u_{h}^{n}\|_{\infty}
\le\frac{C_{u}t_{n}}{1-\eta}\exp\Bigl(\frac{4t_{n}}{1-\eta}\Bigr)\bigl(\tau^{2}+h^{2}\bigr),
\qquad 1\le n\le N,
$$

where **$C_u$ is independent of both the step sizes and the step ratios**. The nonlinearity requires no Lipschitz assumption: the discrete maximum principle puts the solution in $[-1,1]$, and $|(a^{3}-a)-(b^{3}-b)|\le2|a-b|$ for $a,b\in[-1,1]$ gives $\|f(U^n)-f(u^n)\|_\infty\le2\|e^n\|_\infty$, so $\kappa=2$, consistent with taking $S_n=2$ earlier. The paper states this is the **first** such convergence result for nonuniform BDF2 under Grigorieff's zero-stability condition.

### Numerical experiments

Three examples.

| Example | Setup                                                                                                                                                                                                                                               | What it tests                                               |
| ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------- |
| 1       | Forced Allen-Cahn $\partial_tu=\frac{1}{8\pi^2}\Delta u-f(u)+g$ on $(0,1)^2\times(0,1)$, manufactured solution $u=\sin(2\pi x)\sin(2\pi y)\sin t$, **random grid** $\tau_k=T\epsilon_k/S$ with $\epsilon_k\in(0,1)$ random and $S=\sum_k\epsilon_k$ | second-order accuracy in time and robustness to step ratios |
| 2       | Four merging bubbles: $\varepsilon=0.02$, $\Omega=(-1,1)^2$, $128$ grid points per direction, initial data the product of four $\tanh$ profiles of radius $0.2$ centred at $(\pm0.3,0)$ and $(0,\pm0.3)$                                            | qualitative correctness of the interface evolution          |
| 3       | Coarsening dynamics: $\varepsilon=0.01$, $\Omega=(0,1)^2$, $128\times128$ uniform spatial grid, random initial data                                                                                                                                 | max norm and energy across several values of $\tau$         |

The adaptive strategy is $\tau_{\mathrm{ada}}(e,\tau)=\rho\bigl(\frac{tol}{e}\bigr)^{1/2}\tau_{\mathrm{cur}}$ with defaults $\rho=0.6$, $tol=10^{-4}$, $\tau_{\max}=0.1$, $\tau_{\min}=10^{-3}$; the first/second-order pair is backward Euler and adaptive BDF2.

**Example 1 is the most informative of the three.** It uses not a carefully designed grid but a **purely random** step sequence, so the step ratios are uncontrolled; second-order accuracy in time is still observed on such a grid, which says the analysis's robustness to step ratios is not merely on paper. Example 3 tracks max norm and energy, verifying the maximum principle and the energy dissipation respectively.

The coarsening example's initial data is described inconsistently: it is written as $u_{0}=0.95+\mathrm{rand}(\mathbf x)\times0.05$ in one place, while the random numbers are said to vary from $-0.05$ to $0.05$ in another. Those two statements cannot both hold — the first gives data with mean $0.95$, the second suggests a perturbation centred at $0$.

The gap between theory and experiment: the theorem gives an **upper bound** of the form $\tau^2+h^2$, yet example 1 remains second order on random grids where step ratios can be large, which says **S0 and S1 are sufficient rather than necessary**. Paper 52's experiments quantify this much more thoroughly.

### Relation to the others

This is the **hub** of the whole series. It transplants the complementary-kernel idea from the fractional setting ([[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|paper 43]], plus Liao-Li-Zhang and Liao-McLean-Zhang) to integer-order variable-step BDF2, via the recombination that turns the sign-indefinite pair $\{b^{(n)}_0>0,b^{(n)}_1<0\}$ into a non-negative decreasing family $\{d^{(n)}_j\}$. Paper 52 carries the same $3.561$ energy threshold to the molecular beam epitaxy model without slope selection and upgrades the machinery from KRC to **DOC kernels**, additionally obtaining $L^2$ stability and error estimates at the same threshold; paper 58 extends DOC to BDF-$k$ for $3\le k\le5$; paper 67 handles variable-step BDF3 with the tighter $1.4877$; and paper 74 solves the underlying "is this quadratic form positive definite" question in general. Both questions raised in this paper's conclusions were later answered: whether S1/S0 are optimal (partly, by Liao-Ji-Wang-Zhang's $4.864$), and how to construct nonuniform BDF2-type schemes for **time-fractional** phase-field equations (by [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|paper 57]]). The paper also compares itself with Chen-Wang-Yan-Zhang (2019), who analysed nonuniform BDF2 with convex splitting for Cahn-Hilliard using a generalised discrete Grönwall inequality, requiring $r_k<1.53$.

## 52: the same analysis for molecular beam epitaxy

### The idea

The molecular beam epitaxy model without slope selection is harder than Allen-Cahn in two ways. First, its free energy contains a **logarithmic term** that is bounded above (by $0$) but **unbounded below with no relative minimum** — so no value of $\nabla\Phi$ is energetically preferred, there is no slope selection, and the physics is a genuinely multiscale rough-flat-rough morphology. Second, the equation is **fourth order**.

Each of those closes off a route. Fourth order plus the nature of a height function means **there is no maximum principle to be had here**: $\Phi$ is film height, not an order parameter confined to an interval, and what the paper establishes at the continuous level is volume conservation, energy dissipation and an $L^2$ growth bound — none of which produces an $\ell^\infty$ invariant region. Paper 48's KRC-plus-$\ell^\infty$ machinery therefore has **nothing to grip on**, and energy and $L^2$ become the only tools.

The paper's move is to switch to the other side of the convolution. KRC replaces the kernel; **DOC inverts the operator**. Construct $\theta$ so that the lower triangular matrix $\Theta_2=B_2^{-1}$, then multiply the scheme at level $j$ by $\theta^{(n)}_{n-j}$ and sum over $j$: the multistep operator on the left **collapses to a single first difference** $\nabla_\tau u^n$, and the right-hand side becomes a convolution of the spatial and nonlinear terms. **The multistep stability question is thereby reduced to a positive-definiteness question about a quadratic form**, and that is precisely the question that produced $3.561$ in paper 48. This is why two papers on completely different models land on the same threshold.

### Setting

On a bounded $\Omega\subset\mathbb R^2$ with periodic boundary conditions,

$$
\Phi_t=-\varepsilon\Delta^{2}\Phi-\nabla\cdot\mathbf f(\nabla\Phi),
\qquad
\mathbf f(\mathbf v):=\frac{\mathbf v}{1+|\mathbf v|^{2}},
\qquad \Phi(\mathbf x,0)=\Phi_0(\mathbf x),
$$

the $L^2(\Omega)$ gradient flow of the free energy

$$
E[\Phi]=\int_{\Omega}\Bigl[\frac{\varepsilon}{2}(\Delta\Phi)^{2}
-\frac{1}{2}\ln\bigl(1+|\nabla\Phi|^{2}\bigr)\Bigr]\mathrm d\mathbf x .
$$

Here $\Phi$ is the dimensionless height of the film in a co-moving frame, the nonlinear second-order term models the Ehrlich-Schwoebel effect, the linear fourth-order term models surface diffusion, and $\varepsilon>0$ measures the width of the rounded corners on an otherwise faceted crystal film. Three continuous properties: volume conservation $(\Phi(t),1)=(\Phi_0,1)$, energy dissipation

$$
\frac{\mathrm d}{\mathrm dt}E[\Phi]=-\|\Phi_t\|_{L^2(\Omega)}^{2}\le0,
$$

and the $L^2$ estimate $\|\Phi\|_{L^{2}(\Omega)}\le e^{t/(4\varepsilon)}\|\Phi_0\|_{L^{2}(\Omega)}$. Well-posedness is due to Li and Liu (2003): for $\Phi_0\in H^m_{per}(\Omega)$ with $m\ge2$ there is a unique weak solution with $\Phi\in L^{\infty}(0,T;H^m)\cap L^2(0,T;H^{m+2})$ and $\partial_t\Phi\in L^2(0,T;H^{m-2})$.

The time discretisation is identical to paper 48's (same $b^{(n)}_j$, same $r_1\equiv0$ convention embedding BDF1 on the first level), space is central differences, and the scheme is fully implicit:

$$
D_{2}\phi_{h}^{n}+\varepsilon\Delta_{h}^{2}\phi_{h}^{n}
+\nabla_{h}\cdot\mathbf f(\nabla_{h}\phi_{h}^{n})=0,
\qquad 1\le n\le N,\quad \phi_h^0=\Phi_0(\mathbf x_h).
$$

> [!note] A **different** MBE model from the one in paper 78
> This paper treats the model **without slope selection**, whose free energy carries the logarithmic term $-\frac12\ln(1+|\mathbf v|^2)$; [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation|paper 78]] treats the model **with slope selection**, whose free energy carries the double well $\frac14(|\mathbf v|^2-1)^2$. **The case without slope selection is harder**, because the nonlinearity is not polynomially bounded and the free energy has no lower bound. This paper's conclusions say explicitly that its technique does not apply to the molecular beam epitaxy model with slope selection and leaves that as future work. The two must not be conflated.

### Derivation

**Step one: the formal definition of DOC kernels.** This is the paper that writes them out explicitly (attributing the linear-diffusion case to Liao-Zhang):

$$
\theta_{0}^{(n)}:=\frac{1}{b_{0}^{(n)}},
\qquad
\theta_{n-k}^{(n)}:=-\frac{1}{b_{0}^{(k)}}\sum_{j=k+1}^{n}\theta_{n-j}^{(n)}b_{j-k}^{(j)}
\qquad(1\le k\le n-1),
$$

with the **orthogonality identity**

$$
\sum_{j=k}^{n}\theta_{n-j}^{(n)}\,b_{j-k}^{(j)}\;\equiv\;\delta_{nk}
\qquad(1\le k\le n),
$$

$\delta_{nk}$ the Kronecker symbol. The consequence that drives everything is

$$
\sum_{j=1}^{n}\theta_{n-j}^{(n)}\,D_{2}v^{j}\;=\;\nabla_{\tau}v^{n}
\qquad(n\ge1).
$$

That is, **convolving the scheme with the DOC kernels peels off the multistep operator and leaves a single first difference**, so the standard one-step energy argument becomes available again. Compare the DCC kernels of papers 43 and 48: their identity is $\sum_{j=k}^{n}P^{(n)}_{n-j}A^{(j)}_{j-k}\equiv1$, of complementary/summation type, giving $\sum_jP\sum_kA\nabla w^k=w^n-w^0$. **DOC is orthogonal (Kronecker $\delta$); DCC is complementary (constant $1$).**

The matrix form states it most cleanly: with $\mathbf B_2$ the lower bidiagonal matrix of BDF2 kernels ($\otimes\,\mathbf I_2$) and $\mathbf\Theta_2$ the lower triangular matrix of DOC kernels, the orthogonality identity says

$$
\mathbf\Theta_{2}=\mathbf B_{2}^{-1},
$$

while the positivity lemma says $\mathbf B:=\mathbf B_{2}+\mathbf B_{2}^{T}$ is positive definite. **The DOC kernels are literally the entries of the inverse of the BDF2 matrix** — the cleanest possible answer to "what is a DOC kernel".

**Step two: the closed form and three properties.** Under **S1**:

1. the DOC kernels $\theta^{(n)}_{n-j}$ are positive definite;
2. they are positive, with the closed form

$$
\theta_{n-j}^{(n)}=\frac{1}{b^{(j)}_{0}}\prod_{i=j+1}^{n}\frac{r_i^{2}}{1+2r_i},
\qquad 1\le j\le n;
$$

3. $\displaystyle\sum_{j=1}^{n}\theta_{n-j}^{(n)}=\tau_n$, hence $\displaystyle\sum_{k=1}^{n}\sum_{j=1}^{k}\theta_{k-j}^{(k)}=t_n$ for $n\ge1$.

Property 2 exposes the connection to $1+\sqrt2$ plainly: **the product $\prod_i\frac{r_i^2}{1+2r_i}$ decays exactly when each factor is below $1$, i.e. when $r_i<1+\sqrt2$**. Property 3 is the step that converts the spatial consistency error into $t_nh^2$.

**Step three: positive definiteness of the quadratic form and $3.561$.** Under **S1**, for any non-zero sequence $\{w_k\}_{k=1}^n$,

$$
2w_{k}\sum_{j=1}^{k}b_{k-j}^{(k)}w_{j}
\ \ge\ \frac{r_{k+1}}{1+r_{k+1}}\frac{w_{k}^{2}}{\tau_{k}}
-\frac{r_{k}}{1+r_{k}}\frac{w_{k-1}^{2}}{\tau_{k-1}}
+\Bigl(\frac{2+4r_{k}-r_{k}^{2}}{1+r_{k}}-\frac{r_{k+1}}{1+r_{k+1}}\Bigr)\frac{w_{k}^{2}}{\tau_{k}},
\quad k\ge2,
$$

and after summation

$$
\sum_{k=1}^{n}w_{k}\sum_{j=1}^{k}b_{k-j}^{(k)}w_{j}
\ \ge\ \frac{1}{2}\sum_{k=1}^{n}\Bigl(\frac{2+4r_{k}-r_{k}^{2}}{1+r_{k}}
-\frac{r_{k+1}}{1+r_{k+1}}\Bigr)\frac{w_{k}^{2}}{\tau_{k}}\ >\ 0,
\qquad n\ge2 .
$$

The bracket is positive for all admissible step ratios exactly under **S1**: $0<r_{k}<r_{s}=\frac{3+\sqrt{17}}{2}\approx3.561$. **This is the same algebraic fact as the inequality in paper 48, written in different variables** (there with $\partial_\tau u$, here with a general sequence $w$), which is why the two papers arrive at the same constant. The paper notes it agrees with the linear-case result (Liao-Zhang, diffusion equation) and says this appears to be the best result in the literature so far for nonlinear problems.

**Step four: the modified energy.**

$$
\mathcal{E}[\phi^{n}]:=E[\phi^{n}]
+\frac{r_{n+1}}{2(1+r_{n+1})\tau_{n}}\|\nabla_{\tau}\phi^{n}\|^{2},
\qquad 0\le n\le N,
\qquad \mathcal{E}[\phi^{0}]=E[\phi^{0}],
$$

$$
E[\phi^{n}]:=\frac{\varepsilon}{2}\|\Delta_{h}\phi^{n}\|^{2}
-\frac{1}{2}\bigl\langle\ln(1+|\nabla_{h}\phi^{n}|^{2}),1\bigr\rangle .
$$

The correction $\frac{r_{n+1}}{2(1+r_{n+1})\tau_n}\|\nabla_\tau\phi^n\|^2$ and paper 48's $\frac{r_{n+1}\tau_n}{2(1+r_{n+1})}\|\partial_\tau\phi^n\|^2$ are the same quantity written two ways.

**Step five: the mesh-independent constant $\mathcal M_r$.** This is what lets the error estimate be called robust. Using the step-rescaled matrices,

$$
\mathcal{M}_r:=\max_{n\ge1}\|\widetilde{\mathbf B}_{2}\|^{2}\,\|\mathbf L^{-1}\|^{4}
=\max_{n\ge1}\frac{\lambda_{\max}(\widetilde{\mathbf B}_{2}^{T}\widetilde{\mathbf B}_{2})}
{\lambda_{\min}^{2}(\widetilde{\mathbf B})} .
$$

Under **S1** there is the rough estimate $\mathcal M_r<39$; and if the computation does not **persistently** use step ratios near the stability limit $r_s=3.561$, then $\mathcal M_r\le4$. Both numbers are worth remembering as printed.

**Step six: two inequalities for the logarithmic nonlinearity.** Convexity uses

$$
g(\lambda):=\tfrac12\ln(1+|\mathbf u+\lambda\mathbf v|^2),
\qquad
g''(0)=\frac{1-|\mathbf u|^{2}}{(1+|\mathbf u|^{2})^{2}}\mathbf v^T\mathbf v\le\mathbf v^T\mathbf v,
$$

and the energy argument uses the vector inequality

$$
\frac{2(\mathbf u-\mathbf v)^{T}\mathbf u}{1+|\mathbf u|^{2}}
\le\ln\frac{1+|\mathbf u|^{2}}{1+|\mathbf v|^{2}}+|\mathbf u-\mathbf v|^{2} .
$$

The latter is the substitute for a chain rule for the logarithmic potential, playing exactly the role of the bespoke $H$ in [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|paper 57]].

**Step seven: two DOC quadratic-form inequalities (the paper's technical novelty).** (a) A DOC-type Cauchy-Young inequality: under **S1**, for any $\epsilon\ge0$,

$$
\sum_{k=1}^n\sum_{j=1}^k\theta^{(k)}_{k-j}(\mathbf z^k)^T\mathbf w^j
\le\frac{\epsilon}{2}\mathbf z^T\mathbf\Theta\mathbf z
+\frac{1}{2\epsilon}\mathbf w^T\mathbf B^{-1}\mathbf w ;
$$

(b) a DOC inequality tailored to the nonlinear term:

$$
\sum_{k,j}\theta^{(k)}_{k-j}(\mathbf z^k)^T
\bigl[\mathbf f(\mathbf v^j+\mathbf z^j)-\mathbf f(\mathbf v^j)\bigr]
\le2\sqrt{\mathcal M_r}\sum_{k,j}\theta^{(k)}_{k-j}(\mathbf z^k)^T\mathbf z^j,
$$

where **$\mathcal M_r$ is independent of $t_n$, of the step sizes $\tau_n$ and of the step ratios $r_n$**. That mesh-independence is the whole point of the paper.

### Theorems

**(Unique solvability / convexity)** If $\tau_{n}\le4\varepsilon$, the BDF2 scheme is convex and hence uniquely solvable. The key step is that $\tau_n\le4\varepsilon$ gives $b_0^{(n)}>\frac{1}{4\varepsilon}$, enough to dominate the negative contribution coming from $g''(0)\le\mathbf v^T\mathbf v$ above.

**(Discrete energy dissipation law, the $3.561$ conclusion for MBE)** Assume **S1** and

$$
\tau_{n}\le4\varepsilon\min\Bigl\{1,\
\frac{2+4r_{n}-r_{n}^{2}}{1+r_{n}}-\frac{r_{n+1}}{1+r_{n+1}}\Bigr\},
\qquad n\ge1 .
$$

Then

$$
\mathcal{E}[\phi^{n}]\le\mathcal{E}[\phi^{n-1}]\le\mathcal{E}[\phi^{0}]=E[\phi^{0}],
\qquad n\ge1 .
$$

**($L^2$ stability)** Assume **S1** and $\tau_{n}\le\varepsilon/(16\mathcal M_r^{2})$. For perturbed initial data $\bar\phi_h^0$,

$$
\|\bar{\phi}^{n}-\phi^{n}\|\le2\exp\bigl(16\mathcal M_r^{2}\,t_{n-1}/\varepsilon\bigr)\,
\|\bar{\phi}^{0}-\phi^{0}\|,
\qquad 1\le n\le N .
$$

The proof is the DOC move: multiply the perturbation equation by $\theta^{(k)}_{k-j}$, sum over $j=1,\dots,k$, use $\sum_j\theta^{(k)}_{k-j}D_2z^j=\nabla_\tau z^k$, take the inner product with $2z^k$ and sum over $k$.

**($L^2$ convergence)** If $\Phi\in C^{(6,3)}_{\mathbf x,t}(\Omega\times(0,T])$, **S1** holds and $\tau_{n}\le\varepsilon/(16\mathcal M_r^{2})$, then

$$
\|\Phi^{n}-\phi^{n}\|\le C_{\phi}\exp\bigl(16\mathcal M_r^{2}t_{n-1}/\varepsilon\bigr)
\Bigl[\tau_{1}^{2}\sum_{k=1}^{n}\prod_{i=2}^{k}\frac{r_{i}^{2}}{1+2r_{i}}
+t_{n}(\tau^{2}+h^{2})\Bigr],
\qquad 1\le n\le N .
$$

**The first term deserves separate attention**: the starting-step error $\tau_1^2$ is damped by the product $\prod_i\frac{r_i^2}{1+2r_i}$, which is exactly the factor from DOC property 2 above. In other words, **the DOC kernels automatically damp the starting error over time**, a conclusion visible only once the closed form of the kernels has been computed.

### Numerical experiments

**Example 1 (random time grid).** Forced MBE, $\varepsilon=0.1$, $\Omega=(0,2\pi)^2$, manufactured solution $\Phi(x,t)=\cos(t)\sin(x)\sin(y)$, $T=1$, $3000$ spatial grid points, time steps $\tau_k:=T\sigma_k/S$ with $\sigma_k\in(0,1)$ uniformly random. Alongside the $L^2$ errors and rates, the paper's table records $\max r_k$ and $N_1$ — **the number of time levels whose step ratio reaches or exceeds $(3+\sqrt{17})/2$**. Under successive refinement:

| Quantity                                           | Values under successive refinement              |
| -------------------------------------------------- | ----------------------------------------------- |
| $\max_k r_k$                                       | $2.94,\ 11.98,\ 34.82,\ 37.72,\ 71.89,\ 850.80$ |
| $N_1$ (levels with $r_k\ge3.561$)                  | $0,\ 3,\ 7,\ 13,\ 24,\ 49$                      |
| observed temporal order (between consecutive rows) | $1.84,\ 2.29,\ 2.35,\ 2.42,\ 2.00$              |

**This table is the most convincing theory-versus-experiment record in the whole line of work.** The conclusion is direct: even with dozens of time levels violating $3.561$ **by a wide margin** (up to $850.80$, over two hundred times the threshold), second-order convergence persists. So **S1 is sufficient and far from necessary**; variable-step BDF2 is empirically much more robust than the theory guarantees. The paper itself calls $r_s=3.561$ an artificial constant that is due to the condition S1, and this table is the evidence behind that self-assessment.

**Example 2 (adaptive steps)** is used to accelerate the approach to the steady state.

**The boundary the paper draws itself**: the conclusions state explicitly that this technique does not apply to the molecular beam epitaxy model with slope selection, leaving it as future work. That branch is handled by [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation|paper 78]] with entirely different methods (linear relaxation plus a regularised energy reformulation, no variable-step analysis).

### Relation to the others

It is the bridge between paper 48 (KRC/DCC, $\ell^\infty$, Allen-Cahn) and the general DOC theory. It imports DOC kernels from Liao-Zhang's analysis of linear diffusion and shows the toolkit survives a genuinely nonlinear, fourth-order, logarithmic-energy model, delivering energy stability together with $L^2$ stability and convergence at the **same** $3.561$ threshold. Paper 58 then extends DOC to BDF-$k$ with $3\le k\le5$; paper 67 applies the same philosophy to variable-step BDF3 and obtains the smaller $1.4877$; and papers 69, [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation|91 and 104]] apply DOC to time-filtered backward Euler, implicit-explicit Runge-Kutta and implicit-explicit multistep methods respectively. The later Liao-Ji-Wang-Zhang (_J. Sci. Comput._ 92 (2022) 52) explicitly lists this paper as one of three $3.561$ DOC results and relaxes the threshold for the Cahn-Hilliard model to $\approx4.864$ by changing the weighting exponent.

## 58: bringing non-A-stable BDF-$k$ back to a textbook argument

### The idea

A-stable BDF1 and BDF2 admit a direct textbook discrete energy proof: testing with $u^n$ in $L^2$ suffices. For $3\le k\le5$ the BDF-$k$ formulas are **not** A-stable and that direct argument fails. Since Lubich, Mansour and Venkataraman (_IMA J. Numer. Anal._ 33 (2013) 1365-1385) the standard remedy has been the **Nevanlinna-Odeh multiplier technique** (_Numer. Funct. Anal. Optim._ 3 (1981) 377-423), which rests on Dahlquist's equivalence of A-stability and G-stability and tests with $u^n-\sum_i\eta_iu^{n-i}$ instead of $u^n$; another is Liu's telescope formulas. Both introduce artificial multipliers and, crucially, **demand a stronger norm on the starting data** — $H^1$-type quantities appear in the multiplier-based stability estimates.

The question the paper poses and answers is whether the direct discrete energy analysis can be recovered for BDF-$k$ with $3\le k\le5$.

**Its move is not to change the test function but to change the equation.** Since $\Theta=B^{-1}$, multiply the whole scheme on the left by $\Theta$: the multistep operator is peeled away, the left-hand side returns to a plain first difference, and the right-hand side becomes a convolution. Because $B\Theta=\Theta B=I$ holds on both sides, the transform is **reversible** with no loss of information. After that, testing with $u^n$ carries the textbook argument through, and the entire difficulty is concentrated into a single question: **is the quadratic form built from the DOC kernels positive definite?**

On a uniform grid that question has a complete answer. $\{\theta_j\}$ is positive definite if and only if $\{b_j\}$ is, and the latter is decided by the generating function $\mathrm g^{(k)}(\varphi)=2\sum_jb_j^{(k)}\cos(j\varphi)$ of the symmetrised Toeplitz matrix together with the Grenander-Szegő theorem. **The crucial fact is that $\mathrm g^{(k)}_{\min}>0$ for $3\le k\le5$, even though BDF-$k$ is not A-stable.** It is that weaker-than-A-stability property that replaces A-stability in the argument.

### Setting

A linear reaction-diffusion equation on a bounded convex domain $\Omega$ with homogeneous Dirichlet conditions:

$$
\partial_t u-\varepsilon\Delta u=\beta(x,t)\,u+f(t,x),
\qquad x\in\Omega,\ 0<t<T,
\qquad u|_{\partial\Omega}=0,\ u(0,x)=u_0(x),
$$

with $\varepsilon>0$ constant and $|\beta(x,t)|\le\beta^{*}$. **Note the time grid here is uniform**, $t_k=k\tau$ with $\tau:=T/N$; the variable-step case for $k\ge3$ is listed explicitly as an open problem in the concluding remarks.

Writing BDF-$k$ as a discrete convolution

$$
D_{k}v^n:=\frac1{\tau}\sum_{j=1}^n b_{n-j}^{(k)}\,\nabla_\tau v^{j},
\qquad n\ge k,
$$

its kernels vanish for $j\ge k$ and come from the generating function

$$
\sum_{\ell=1}^{k}\frac{1}{\ell}(1-\zeta)^{\ell-1}
=\sum_{\ell=0}^{k-1}b_{\ell}^{(k)}\zeta^{\ell},
\qquad 3\le k\le 5,
$$

explicitly

| BDF-$k$ | $b_0^{(k)}$ | $b_1^{(k)}$ | $b_2^{(k)}$ | $b_3^{(k)}$ | $b_4^{(k)}$ |
| ------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| $k=2$   | $3/2$       | $-1/2$      |             |             |             |
| $k=3$   | $11/6$      | $-7/6$      | $1/3$       |             |             |
| $k=4$   | $25/12$     | $-23/12$    | $13/12$     | $-1/4$      |             |
| $k=5$   | $137/60$    | $-163/60$   | $137/60$    | $-21/20$    | $1/5$       |

**The signs alternate and there is more than one negative entry** — which is why what is needed here is not a non-negative kernel (impossible) but a positive-definite quadratic form. The time-discrete scheme is (with $u^1,\dots,u^{k-1}$ assumed supplied by a starting procedure)

$$
D_{k}u^{j}=\varepsilon\Delta u^{j}+\beta^j u^{j}+f^j,
\qquad k\le j\le N,
$$

with weak form $\langle D_k u^j,w\rangle+\varepsilon\langle\nabla u^j,\nabla w\rangle=\langle\beta^j u^j,w\rangle+\langle f^j,w\rangle$ for all $w\in H^1_0(\Omega)$.

### Derivation

**Step one: DOC-$k$ kernels and the reversible transform.**

$$
\theta_{0}^{(k)}:=\frac{1}{b_{0}^{(k)}},
\qquad
\theta_{n-j}^{(k)}:=-\frac{1}{b_{0}^{(k)}}\sum_{\ell=j+1}^{n}\theta_{n-\ell}^{(k)}\,b_{\ell-j}^{(k)},
\qquad j=n-1,n-2,\dots,k,
$$

(with the convention $\sum_{k=i}^{j}\cdot=0$ when $i>j$), satisfying the orthogonality identity

$$
\sum_{\ell=j}^{n}\theta_{n-\ell}^{(k)}\,b_{\ell-j}^{(k)}\equiv\delta_{nj}
\qquad(k\le j\le n),
$$

and mutual orthogonality $\sum_{\ell=j}^{n}b_{n-\ell}^{(k)}\theta_{\ell-j}^{(k)}\equiv\delta_{nj}$ ($\delta_{nj}$ is typeset as $\delta_{mk}$ here, an evident typo). Hence

$$
\sum_{j=k}^{n}\theta_{n-j}^{(k)}\sum_{\ell=k}^{j}b_{j-\ell}^{(k)}\nabla_\tau u^{\ell}
=\sum_{\ell=k}^{n}\nabla_\tau u^{\ell}\sum_{j=\ell}^{n}\theta_{n-j}^{(k)}b_{j-\ell}^{(k)}
=\nabla_\tau u^{n},
$$

that is, **acting with the DOC kernels inverts the BDF-$k$ convolution back to a single first difference**, up to a starting-value remainder:

$$
\sum_{j=k}^{n}\theta_{n-j}^{(k)}D_{k}u^{j}
=\frac1{\tau}u_{\mathrm I}^{(k,n)}+\partial_{\tau}u^{n},
\qquad
u_{\mathrm I}^{(k,n)}:=\sum_{\ell=1}^{k-1}\nabla_\tau u^{\ell}
\sum_{j=k}^{n}\theta_{n-j}^{(k)}b_{j-\ell}^{(k)} .
$$

The transformed weak scheme reads: for all $w\in H^1_0(\Omega)$ and $k\le j\le N$,

$$
\langle\partial_\tau u^{j},w\rangle
+\varepsilon\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle\nabla u^{\ell},\nabla w\rangle
=-\frac1\tau\langle u_{\mathrm I}^{(k,j)},w\rangle
+\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle\beta^{\ell}u^{\ell},w\rangle
+\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle f^{\ell},w\rangle .
$$

Taking $w=2\tau u^j$ and summing over $j=k,\dots,n$ (discarding the term $\sum_j\|u^j-u^{j-1}\|^2$) gives an energy inequality in **classical form**:

$$
\|u^n\|^2-\|u^{k-1}\|^2
\le -2\sum_{j=k}^{n}\langle u_{\mathrm I}^{(k,j)},u^j\rangle
-2\varepsilon\tau\sum_{j=k}^{n}\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle\nabla u^{\ell},\nabla u^{j}\rangle
+2\tau\sum_{j=k}^{n}\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle\beta^{\ell}u^{\ell},u^{j}\rangle
+2\tau\sum_{j=k}^{n}\sum_{\ell=k}^{j}\theta_{j-\ell}^{(k)}\langle f^{\ell},u^{j}\rangle .
$$

Three technical ingredients close the argument: (i) positive definiteness of $\theta_j^{(k)}$; (ii) decay of $\theta_j^{(k)}$; (iii) decay of the starting term $u_{\mathrm I}^{(k,j)}$.

**Step two (i): positivity decided by the generating function.** Lemma 2.1 establishes the equivalence: $\{b_j^{(k)}\}$ is positive (semi-)definite if and only if $\{\theta_j^{(k)}\}$ is. And positivity of $\{b_j^{(k)}\}$ is decided by the generating function of the symmetrised Toeplitz matrix $B_k=B_{k,l}+B_{k,l}^{T}$,

$$
\mathrm g^{(k)}(\varphi)=2\sum_{j=0}^{k-1}b_j^{(k)}\cos(j\varphi),
$$

together with Grenander-Szegő: $\mathrm g^{(k)}_{\min}\le\lambda_{\min}(B_k)\le\lambda_{\max}(B_k)\le\mathrm g^{(k)}_{\max}$. Explicitly,

$$
\mathrm g^{(3)}(\varphi)=\tfrac13\bigl(11-7\cos\varphi+2\cos2\varphi\bigr)
=\tfrac43\bigl(\cos\varphi-\tfrac78\bigr)^2+\tfrac{95}{48},
$$

$$
\mathrm g^{(4)}(\varphi)=\tfrac16\bigl(25-23\cos\varphi+13\cos2\varphi-3\cos3\varphi\bigr),
$$

$$
\mathrm g^{(5)}(\varphi)=\tfrac1{30}\bigl(137-163\cos\varphi+137\cos2\varphi-63\cos3\varphi+12\cos4\varphi\bigr).
$$

For $k=3$ the completed square exhibits the positive lower bound at a glance. This yields the quadratic-form bound (Lemma 2.4): for $3\le k\le5$ and any real sequence $\{w_j\}$,

$$
2\sum_{m=k}^{n}w_m\sum_{j=k}^{m}b_{m-j}^{(k)}w_j
\ \ge\ \sigma_k\sum_{j=k}^{n}w_j^2,
\qquad n\ge k,
$$

with **explicit constants**

$$
\sigma_3=\frac{95}{48}\approx1.97919,
\qquad
\sigma_4=\frac{2656-43\sqrt{43}}{1458}\approx1.62828,
\qquad
\sigma_5\approx0.477683 .
$$

For $k=4$ the minimiser of $Z_4(x)=12-14x+26x^2-12x^3$ on $[-1,1]$ is $x_{*}=(13-\sqrt{43})/18$; for $k=5$ the minimiser of $Z_5(x)=12+26x+178x^2-252x^3+96x^4$ is

$$
x^{*}=\tfrac{1}{96}\Bigl(63-\sqrt[3]{49041-16\sqrt{3891895}}
-\tfrac{1121}{\sqrt[3]{49041-16\sqrt{3891895}}}\Bigr)\approx-0.064041,
\qquad Z_5(x^{*})\approx14.3305 .
$$

**That $\sigma_k$ shrinks rapidly with $k$ ($1.979\to1.628\to0.478$) carries meaning**: it quantifies how far BDF-$k$ is from losing positive definiteness, and by $k=5$ less than a quarter of the $k=3$ margin remains.

**Step three (ii): geometric decay of the DOC-$k$ kernels.** Although $b_j^{(k)}=0$ for $j\ge k$, the DOC kernels $\theta_j^{(k)}$ **never vanish** — inversion turns a finitely supported kernel into an infinitely long one. Fortunately they decay geometrically:

$$
\bigl|\theta_j^{(k)}\bigr|\le\frac{\rho_k}{4}\Bigl(\frac{k}{7}\Bigr)^{j},
\qquad 3\le k\le5,\ j\ge0,
\qquad
\rho_3=\frac{10}{3},\quad \rho_4=6,\quad \rho_5=\frac{96}{5}.
$$

The ratio $k/7$ (that is $3/7$, $4/7$, $5/7$) produces, through $\sum_{j\ge\ell}(k/7)^{j-\ell}\le\frac{7}{7-k}$, the factor $\frac{7}{7-k}$ that appears everywhere in the final constants — $\frac74$, $\frac73$, $\frac72$ for $k=3,4,5$.

**Step four (iii): decay of the starting term.** There exist finite $c_{\mathrm I,k}>1$ with

$$
\bigl|u_{\mathrm I}^{(k,j)}\bigr|\le\frac{c_{\mathrm I,k}\rho_k}{8}
\Bigl(\frac{k}{7}\Bigr)^{j-k}\sum_{\ell=1}^{k-1}\bigl|\nabla_\tau u^{\ell}\bigr|,
\qquad
\sum_{j=k}^{n}\bigl|u_{\mathrm I}^{(k,j)}\bigr|
\le\frac{7c_{\mathrm I,k}\rho_k}{8(7-k)}\sum_{\ell=1}^{k-1}\bigl|\nabla_\tau u^{\ell}\bigr| .
$$

The paper gives only $c_{\mathrm I,3}=11/7$ explicitly; $c_{\mathrm I,4}$ and $c_{\mathrm I,5}$ are asserted to be finite (the argument absorbs the finitely many low-index cases $n=4,5$ and $n=5,6,7$ separately into the constant). **Their numerical values are not given.**

### Theorems

**(Dissipative case $\beta=\beta(x)\le0$)** For $3\le k\le5$ and $n\ge k$,

$$
\|u^n\|\le\|u^{k-1}\|
+\frac{7c_{\mathrm I,k}\rho_k}{4(7-k)}\sum_{\ell=1}^{k-1}\|\nabla_\tau u^{\ell}\|
+\frac{7\rho_k}{2(7-k)}\sum_{\ell=k}^{n}\tau\|f^{\ell}\|
\ \le\ \frac{7\rho_k}{2(7-k)}\Bigl(c_{\mathrm I,k}\sum_{\ell=0}^{k-1}\|u^{\ell}\|
+\sum_{\ell=k}^{n}\tau\|f^{\ell}\|\Bigr).
$$

**This one carries no step restriction whatsoever.**

**(General bounded $\beta$)** If $|\beta(x,t)|\le\beta^{*}$ and the uniform step satisfies

$$
\tau\le\frac{7-k}{7\rho_k\beta^{*}},
$$

then for $k\le n\le N$,

$$
\|u^n\|\le\frac{7\rho_k}{7-k}\exp\Bigl(\frac{7\rho_k}{7-k}\beta^{*}t_{n-k}\Bigr)
\Bigl(c_{\mathrm I,k}\sum_{\ell=0}^{k-1}\|u^{\ell}\|+\sum_{\ell=k}^{n}\tau\|f^{\ell}\|\Bigr).
$$

The paper prints only the symbolic form $\tau\le(7-k)/(7\rho_k\beta^{*})$; substituting the values of $\rho_k$, this restriction reads $\tau\le6/(35\beta^{*})$ for $k=3$, $\tau\le1/(14\beta^{*})$ for $k=4$ and $\tau\le2/(134.4\,\beta^{*})$ for $k=5$.

**($L^2$ convergence, order $k$)** Under the same step restriction and the consistency bound $|\eta^j|=|D_k u(t_j)-\partial_t u(t_j)|\le C_u\tau^{k}\max_{t_k\le t\le T}|\partial_t^{(k+1)}u(t)|\le C_u\tau^{k}$ for $j\ge k$,

$$
\|u(t_n)-u^{n}\|\le\frac{7\rho_k c_{\mathrm I,k}}{7-k}
\exp\Bigl(\frac{7\rho_k\beta^{*}t_{n-k}}{7-k}\Bigr)
\Bigl(\sum_{\ell=0}^{k-1}\|u(t_{\ell})-u^{\ell}\|+C_u\,t_{n-k}\,\tau^{k}\Bigr),
\qquad k\le n\le N .
$$

Order $k$ in time, and **the starting-value errors enter in the $L^2$ norm only** — precisely the regularity advantage claimed relative to the multiplier technique, which needs stronger norms on the starting data (the paper cites Proposition 5.1 and Theorem 5.1 of Akrivis-Katsoprinakis, _Math. Comp._ 85 (2016) 2195-2216 for contrast).

### Numerical experiments

**This paper has no numerical experiments.** It is a pure analysis paper, and the only two sets of figures it contains are: (a) the generating functions $\mathrm g^{(k)}(\varphi)$ on $[-\pi,\pi]$ for $3\le k\le5$, showing $\mathrm g^{(k)}_{\min}>0$; and (b) the computed DOC-$k$ kernels plotted against the bound $\frac{\rho_k}{4}(k/7)^{x}$, showing the geometric decay.

Those two sets correspond one-to-one with the two lemmas above, so they are numerical verification of **lemmas**, not of the **scheme**. **The paper therefore contains no measured convergence rates**; order $k$ is a theorem, not an observation.

### Relation to the others

The DOC technique itself comes from Liao-Zhang (_Math. Comp._ 90 (2021) 1207-1226, where it serves **variable-step BDF2**) and the companion analysis in paper 52; papers 48 and 52 apply DOC/DCC to variable-step BDF2 for nonlinear phase-field problems, while this paper is the **uniform-grid, higher-order ($k=3,4,5$)** branch of the same programme. It is the direct methodological ancestor of paper 67 (variable-step BDF3): paper 67 takes up exactly the third item in this paper's concluding remarks (studying the discrete energy technique for BDF-$k$, $3\le k\le5$, with variable steps) and obtains the threshold $r_k<1.4877$ there. Its positivity machinery (Grenander-Szegő applied to symmetrised Toeplitz forms built from convolution kernels) is one version of the same thing that paper 74 generalises to **variable-step L1-type kernels**. And the route it displaces — Akrivis et al.'s Nevanlinna-Odeh multipliers and G-stability — is exactly what [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/imex-and-relaxation|paper 104]] revisits with semi-generating-function arguments.

## 67: variable-step BDF3

### The idea

Variable-step BDF3 previously had essentially one classical result: Calvo and Grigorieff (2002) proved $L^2$ stability under the step-ratio condition $r_k<1.199$, with an estimate of the form

$$
\|u^n\|\le C\exp(C\Gamma_n)\Bigl(\|u_0\|+\sum_{j=1}^{n}\tau_j\|f^j\|\Bigr),
\qquad \Gamma_n:=\sum_{k=2}^{n}|r_k-r_{k-1}| .
$$

**The problem is not the threshold but the prefactor $\exp(C\Gamma_n)$, which is not mesh robust.** Take the alternating steps $\{\tau_1,\mu\tau_1,\tau_1,\mu\tau_1,\dots\}$ with $\mu\ne1$ and fix $T=\frac M2(1+\mu)\tau_1$; then

$$
\Gamma_M=(M-1)\bigl|\mu-\mu^{-1}\bigr|\longrightarrow\infty
\qquad(\tau_1\to0),
$$

so the bound degenerates precisely in the adaptive regime it is supposed to cover. **A variable-step theorem that only says something when the mesh is nearly uniform is of limited value.**

The paper's trade is explicit: replace it with an analysis whose constants are **completely independent of the step ratios**, at the price of a step-ratio threshold. This happens on two levels. On the energy level it constructs an explicit Lyapunov functional $G$ with

$$
2v_n\tau_n\sum_{j}d^{(n)}_{n-j}v_j
=G[v_n,v_{n-1}]-G[v_{n-1},v_{n-2}]+F[v_n,v_{n-1},v_{n-2}],
\qquad F\ge\frac{\tau_n}{50}v_n^2 .
$$

That is the same pattern as papers 48 and 52 — telescoping cancellation plus a non-negative remainder — except BDF3 is a three-term operator, so $G$ must **carry two levels at once**. On the $L^2$ level it uses DOC kernels, whose key property is $\sum_j|\vartheta|\le K_3$ with $K_3$ independent of the step ratios, which is where the "mesh robust" in the title comes from.

### Setting

On a bounded convex $\Omega$, $\partial_t u-\varepsilon\Delta u=\kappa(x)u+f(t,x)$ with $u|_{\partial\Omega}=0$, $\varepsilon>0$ constant and $|\kappa(x)|\le\kappa^{*}$. Grid and step ratios as before, $r_k:=\tau_k/\tau_{k-1}$ for $2\le k\le N$.

The variable-step BDF3 formula (Calvo-Grigorieff form) is

$$
D_3v^n=d_0(r_n,r_{n-1})\partial_\tau v^{n}
+d_1(r_n,r_{n-1})\partial_\tau v^{n-1}
+d_2(r_n,r_{n-1})\partial_\tau v^{n-2},
$$

$$
d_0(x,y):=\frac{1+2x}{1+x}+\frac{xy}{1+y+xy},
$$

$$
d_1(x,y):=-\frac{x}{1+x}-\frac{xy}{1+y+xy}
-\frac{xy^2}{1+y+xy}\cdot\frac{1+x}{1+y},
\qquad
d_2(x,y):=\frac{xy^2}{1+y+xy}\cdot\frac{1+x}{1+y},
$$

for $x,y\ge0$. As a convolution $D_3v^n=\sum_{j=1}^{n}d^{(n)}_{n-j}\partial_\tau v^{j}$ with $d^{(n)}_j:=d_j(r_n,r_{n-1})$ for $j=0,1,2$ and $d^{(n)}_j:=0$ for $n\ge j+1\ge4$. The scheme is $D_3u^k=\varepsilon\Delta u^k+\kappa u^k+f^k$ for $3\le k\le N$, with $u^1,u^2$ from a starting procedure.

### Derivation

**Step one: variable-step DOC kernels.**

$$
\vartheta_0^{(n)}:=\frac{1}{d_0^{(n)}},
\qquad
\vartheta_{n-j}^{(n)}:=-\frac{1}{d_0^{(j)}}\sum_{i=j+1}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-j},
\qquad 3\le j\le n-1,
$$

with orthogonality identity $\sum_{i=j}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-j}\equiv\delta_{nj}$ for $3\le j\le n$, equivalently $\Theta_3D_3=I_{m\times m}$ with $m=n-2$; and $D_3\Theta_3=I$ gives mutual orthogonality $\sum_{i=j}^{n}d_{n-i}^{(n)}\vartheta^{(i)}_{i-j}\equiv\delta_{nj}$.

**Note that the variable-step DOC recursion uses $1/d_0^{(j)}$, whose index tracks the running index $j$, unlike the fixed $1/b_0^{(k)}$ of the uniform-grid version in paper 58.** After acting with the DOC kernels the equivalent convolution scheme is

$$
\partial_\tau u^n=-\mathcal I_3^n[u]
+\sum_{k=3}^{n}\vartheta_{n-k}^{(n)}\bigl(\varepsilon\Delta u^k+\kappa u^k\bigr)
+\sum_{k=3}^{n}\vartheta_{n-k}^{(n)}f^k,
\qquad 3\le n\le N,
$$

$$
\mathcal I_3^n[v]:=\sum_{j=1}^{2}\partial_\tau v^{j}\sum_{i=3}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-j}
=\partial_\tau v^2\sum_{i=3}^{n}\vartheta_{n-i}^{(n)}d^{(i)}_{i-2}
+\vartheta_{n-3}^{(n)}d^{(3)}_2\,\partial_\tau v^1 .
$$

The division of labour is: **the energy (gradient-structure) argument uses the original form, the $L^2$ argument uses the DOC-transformed form.**

**Step two: where $R_e\approx1.4877$ comes from.** The threshold is not a natural constant but the outcome of a **parameter trade-off**, worth spelling out. With $\gamma=7/10$, $R_e$ is the **unique positive root** of

$$
d_1(R_e,0)+\tfrac{7}{10}\sqrt{R_e}\,d_2(R_e,R_e)=0
\qquad\Longleftrightarrow\qquad
\frac{10}{7(R_e+1)}-\frac{R_e^2\sqrt{R_e}}{R_e^2+R_e+1}=0,
$$

numerically $R_e\approx1.4877$.

The value $\gamma=7/10$ arises as follows. The discrete gradient decomposition needs two conditions, $q_{n+1}\ge0$ and $p_{n+1}>0$, holding simultaneously in the five variables $r_{n+1},r_n,r_{n-1},\gamma,R_e$, which cannot be solved exactly. The paper retreats to a constant-ratio grid: $q_{n+1}\ge0$ with $r_{n-1}=0$ and $r_{n+1}=r_n=r$ forces $\gamma\le-d_1(r,0)/(\sqrt r\,d_2(r,r))$, and imposing the second condition on a constant-ratio grid as well gives $\bar R_e\approx1.4965$ with $\bar\gamma\approx0.6924$. The authors then **fix $\gamma=7/10$** (close to that $0.6924$) for tractability, yielding $R_e\approx1.4877$. The justification is that $q_{n+1}\ge0$ is necessary and sharp whereas $p_{n+1}>0$ can be relaxed.

**Step three: the sharpness of the threshold is quantified.** Positive definiteness of $\{\tau_nd^{(n)}_{n-k}\}$ is governed by the pentadiagonal symmetric matrix $B_3=B_L+B_L^T$ ($B_L$ the lower triangular matrix with entries $\tau_nd^{(n)}_j$). For the step-rescaled $\widetilde B_3=\Lambda_\tau^{-1}(B_L+B_L^T)\Lambda_\tau^{-1}$ with $\Lambda_\tau=\mathrm{diag}(\sqrt{\tau_3},\dots,\sqrt{\tau_n})$, taking the minimum over $200$ runs on random grids with $r_k\sim U(0,R_e)$ gives these smallest eigenvalues:

| $n$   | $R_e=1.20$ | $R_e=1.50$ | $R_e=1.69$  | $R_e=1.70$  |
| ----- | ---------- | ---------- | ----------- | ----------- |
| $50$  | $1.12$     | $5.08$e-01 | $6.12$e-02  | $-4.55$e-02 |
| $100$ | $1.07$     | $4.35$e-01 | $4.58$e-02  | $-5.29$e-02 |
| $200$ | $1.08$     | $4.18$e-01 | $-2.06$e-02 | $-8.49$e-02 |

**So numerically $R_e<1.69$ is necessary while the theory delivers $R_e<1.4877$ as sufficient, and the gap between them is small.** That contrasts instructively with paper 52's self-assessment of $3.561$ as an artificial constant due to condition S1 — both constants are chosen for convenience, but here the distance to necessity has been measured. Note that the $n=200$, $R_e=1.69$ entry has already gone negative, so the necessary value itself drifts slowly downward with $n$.

**Step four: the discrete gradient structure (Lemma 3.3, the paper's analytical device).** For $0<x,y,z<R_e$ define

$$
d_*(x,y):=-\tfrac{10}{7}\sqrt{x}\,d_1(x,y)-\sqrt{xy}\,d_2(x,y),
$$

$$
p(x,y,z):=2d_0(y,z)-\sqrt{yz}\,d_2(y,z)-\tfrac{49}{100}d_*(y,z)-d_*(x,y),
\qquad
q(x,y,z):=d_*(y,z)-\sqrt{xy}\,d_2(x,y).
$$

If $0<r_k<R_e$, there exist non-negative functionals $G$ and $F$ with

$$
2v_n\tau_n\sum_{j=3}^{n}d^{(n)}_{n-j}v_j
=G[v_n,v_{n-1}]-G[v_{n-1},v_{n-2}]+F[v_n,v_{n-1},v_{n-2}],
\qquad n\ge3,
$$

$$
G[v_n,v_{n-1}]:=d_*(r_{n+1},r_n)\,\tau_nv_n^2
+\sqrt{r_{n+1}r_n}\,d_2(r_{n+1},r_n)
\Bigl(\tfrac{7}{10}\sqrt{\tau_n}\,v_n-\sqrt{\tau_{n-1}}\,v_{n-1}\Bigr)^2,
$$

$$
F[v_n,v_{n-1},v_{n-2}]:=p(r_{n+1},r_n,r_{n-1})\tau_nv_n^2
+q(r_{n+1},r_n,r_{n-1})\Bigl(\tfrac{7}{10}\sqrt{\tau_n}v_n-\sqrt{\tau_{n-1}}v_{n-1}\Bigr)^2
$$

$$
\qquad\qquad
+\sqrt{r_nr_{n-1}}\,d_2(r_n,r_{n-1})
\Bigl(\sqrt{\tau_n}v_n-\tfrac{7}{10}\sqrt{\tau_{n-1}}v_{n-1}+\sqrt{\tau_{n-2}}v_{n-2}\Bigr)^2
\ \ge\ \frac{\tau_n}{50}v_n^2 .
$$

The $7/10$ appearing throughout $G$ and $F$ is the $\gamma$ above, and the constant $1/50$ is the source of the explicit positivity constant. Hence (Lemma 3.4) if $0<r_k<R_e\approx1.4877$ then

$$
2\sum_{k=3}^{n}\xi_k\sum_{j=3}^{k}\tau_kd^{(k)}_{k-j}\xi_j
\ \ge\ \frac{1}{50}\sum_{k=3}^{n}\tau_k\xi_k^2,
\qquad n\ge3 .
$$

**Step five: positivity and uniform summability of the DOC kernels.** Under the same condition (Lemma 4.1) $\{\tau_n\vartheta^{(n)}_{n-k}\}$ is positive definite. The proof is the standard DOC move: set $\eta_k=\sum_{j=3}^{k}\vartheta^{(k)}_{k-j}\xi_j$, use mutual orthogonality to get $\sum_kd^{(n)}_{n-k}\eta_k=\xi_n$, then apply Lemma 3.4 to $\{\eta_k\}$.

There is also (from Lemma 3.1 of Li-Liao): if $r_k\le R_e$, there exists $K_3>0$, **independent of both $t_n$ and the step ratios $r_k\in(0,R_e]$**, with

$$
\sum_{j=3}^{n}\bigl|\vartheta^{(n)}_{n-j}\bigr|\le K_3,
\qquad
\sum_{j=i}^{n}\bigl|\vartheta^{(j)}_{j-i}\bigr|\le K_3 .
$$

**The DOC kernels need not be positive, but they decay exponentially**, which is what makes the final estimate mesh robust.

(The subject here is the linear diffusion equation, so there is **no** maximum-principle result.)

### Theorems

**(Discrete energy dissipation law)** Let $\kappa\le0$ and $f\equiv0$. If $0<r_k<R_e\approx1.4877$ for $k\ge2$, the variable-step BDF3 scheme is **unconditionally energy stable**: $E^n\le E^{n-1}$ for $n\ge3$, where the **modified discrete energy** is

$$
E^n:=\varepsilon\|\nabla u^n\|^2+\bigl(-\kappa u^n,u^n\bigr)
+\bigl(1,\,G[\partial_\tau u^n,\partial_\tau u^{n-1}]\bigr),
\qquad n\ge2 .
$$

Note this energy is of $H^1$-seminorm type and carries a non-negative $G$ correction (a Lyapunov-type correction), structurally the same as the modified energies of papers 48 and 52.

**($L^2$ stability, dissipative case $\kappa<0$)** If $0<r_k<R_e$,

$$
\|u^n\|\le\|u^2\|+K_3\tau\|\partial_\tau u^1\|+4K_3\tau\|\partial_\tau u^2\|
+2\sum_{k=3}^{n}\tau_k\sum_{j=3}^{k}\bigl|\vartheta^{(k)}_{k-j}\bigr|\,\|f^j\|
\ \le\ \|u^2\|+K_3\tau\|\partial_\tau u^1\|+4K_3\tau\|\partial_\tau u^2\|
+2K_3t_n\max_{3\le k\le n}\|f^k\| .
$$

**($L^2$ stability, general bounded $\kappa$)** If $|\kappa|\le\kappa^{*}$, $0<r_k<R_e$ and $\tau\le 1/(4K_3\kappa^{*})$, then for $3\le n\le N$,

$$
\|u^n\|\le2\exp\bigl(4K_3\kappa^{*}t_{n-1}\bigr)
\Bigl(\|u^2\|+K_3\tau\|\partial_\tau u^1\|+4K_3\tau\|\partial_\tau u^2\|
+2K_3t_n\max_{3\le k\le n}\|f^k\|\Bigr).
$$

**(Third-order $L^2$ convergence)** Under the same hypotheses, the error $\tilde u^n=u(t_n)-u^n$ satisfies

$$
\|\tilde u^n\|\le2\exp\bigl(4K_3\kappa^{*}t_{n-1}\bigr)
\Bigl(\|\tilde u^2\|+K_3\tau\|\partial_\tau\tilde u^1\|
+4K_3\tau\|\partial_\tau\tilde u^2\|+2K_3K_ut_n\tau^{3}\Bigr),
\qquad 3\le n\le N,
$$

where **$K_3$ and $K_u$ are independent of $t_n$, of the step sizes $\tau_n$ and of the step ratios $r_n$ — even as $r_n\to R_e$**. That independence is the "mesh robustness" of the title, and it is exactly what the Calvo-Grigorieff prefactor $\exp(C\Gamma_n)$ lacks.

### Numerical experiments

The heat equation $\partial_t u-\varepsilon\Delta u=f$ on $\Omega=(0,2\pi)^2$ with periodic boundary conditions, $\varepsilon=0.1$, manufactured solution $u=\cos(t)\sin(x)\sin(y)$, $T=1$; the error is $e(N)=\max_{1\le n\le N}\|v(t_n)-v^n\|$. Two families of grids: (a) periodic steps $\{\tau_1,\mu\tau_1,\tau_1,\mu\tau_1,\dots\}$ with $\tau_1=2/(N(1+\mu))$ and $r_{\max}=\mu$; (b) random steps $\tau_k=\epsilon_k/\sum\epsilon_k$ with $\epsilon_k\sim U(0,1)$. Starting values come from a two-stage third-order SDIRK method or from variable-step BDF2.

| Grid and starting procedure   | observed order for $N=160,\dots,1280$ | levels violating $R_e$ |
| ----------------------------- | ------------------------------------- | ---------------------- |
| $\mu=2R_e$, Runge-Kutta start | $2.98,\ 2.99,\ 3.00,\ 3.00$           | about $N_1=N/2$        |
| $\mu=4R_e$, BDF2 start        | $2.98,\ 2.99,\ 3.00,\ 2.99$           | about $N_1=N/2$        |

Three conclusions. First, **roughly half the time levels exceed the theoretical ratio limit and the scheme remains stable and third order** — the threshold is sufficient but conservative. Second, both third-order SDIRK and second-order BDF2 suffice as starting procedures to reach third-order accuracy, matching the convergence theorem's prediction (starting values enter only through $\|\tilde u^2\|$, $\tau\|\partial_\tau\tilde u^1\|$ and $\tau\|\partial_\tau\tilde u^2\|$, so a second-order start is enough). Third, on random grids the method is likewise mesh robust and third order, even with many step ratios far above $R_e$.

The smallest-eigenvalue table in step three above is also a numerical experiment, only one that tests a lemma rather than the scheme. Together the two tables settle the status of this threshold: **$1.4877$ is sufficient, $1.69$ is necessary, and in practice exceeding it is often harmless.**

### Relation to the others

The paper describes itself as one of a series on discrete energy analysis of variable-step BDF discretisations, with paper 48 (_SIAM J. Numer. Anal._ 58:2294-2314) and Liao-Zhang (_Math. Comp._ 90:1207-1226) as predecessors, both cited by number in the abstract. It answers the $k=3$ case of the third open problem in paper 58's concluding remarks (studying the discrete energy technique for BDF-$k$, $3\le k\le5$, with variable steps). The DOC apparatus is structurally identical to that of papers 48, 52 and 58; the additional device here is the explicit Lyapunov functional $G$ realising the **discrete gradient structure**, the same tool used in the modified energies of papers 48 and 52.

The ordering of step-ratio constants is worth writing out:

$$
1.199\ <\ 1.4877\ <\ 2.4142\ <\ 2.553\ <\ 3.5616\ <\ 4.8645 .
$$

Here $1.199$ is Calvo-Grigorieff's (2002) classical BDF3 $L^2$ threshold (with the non-robust prefactor this paper removes), $1.4877$ is this paper's BDF3 energy and $L^2$ threshold, $2.4142=1+\sqrt2$ is paper 48's S0, $2.553$ is Li-Liao's $R_3$ for BDF3 applied to ODEs (cited and used here), $3.5616$ is S1, and $4.8645$ belongs to Liao-Ji-Wang-Zhang (outside this topic). **Numerical size does not translate into strength**: these are thresholds for different schemes and different conclusions.

## 69: variable-step time-filtered backward Euler

### The idea

Time filtering is a cheap post-processing step from numerical weather prediction: take one backward Euler step, then add a linear combination of already-computed levels, and the accuracy rises from first to second order at the cost of one extra line of code in a legacy solver. Guzel and Layton (_BIT_ 58 (2018) 1-15) proved that with constant steps the filtered backward Euler (FiBE) scheme is second order and A-stable, and DeCaria, Guzel, Layton and Li (_SIAM J. Sci. Comput._ 43(3) (2021) A2130-A2160) built a variable-step, variable-order family on that basis.

The obstruction to a variable-step theory is structural, and the follow-up literature states it plainly: **the filtering algorithm couples the filtered and unfiltered solutions, so it is a predictor-corrector scheme rather than a one-step method, and variable-step FiBE is not A-stable** — which invalidates the A-stability/G-stability arguments used in the constant-step case.

The paper's bridge is a **one-leg reformulation**. After eliminating the intermediate quantity, FiBE is equivalent to a one-leg multistep (OLM) method whose **left-hand side is exactly the variable-step BDF2 difference operator**, with the right-hand side evaluating $f$ at a second-order accurate combination of $u^{n},u^{n-1},u^{n-2}$. Once the left-hand side is BDF2, the gradient-structure machinery of papers 48, 52 and 67 can be carried over — which is what the abstract means when it says the discrete gradient structure is established for the one-leg multistep scheme associated with the time-filtered backward Euler scheme.

### The scheme

One step of variable-step FiBE is due to DeCaria, Guzel, Layton and Li. Writing $\tau=k_{n+1}/k_n$ for the ratio of the current step to the previous one:

$$
\text{backward Euler predictor:}\qquad
\frac{y^{1}_{n+2}-y_{n+1}}{k_{n+1}}=f\bigl(t_{n+2},y^{1}_{n+2}\bigr),
$$

$$
\text{time filter:}\qquad
y_{n+2}=y^{1}_{n+2}-\frac{\tau(1+\tau)}{1+2\tau}
\Bigl(\frac{1}{1+\tau}\,y^{1}_{n+2}-y_{n+1}+\frac{\tau}{1+\tau}\,y_{n}\Bigr).
$$

Rewriting in Liao's notation ($u^n$ at $t_n$, $\tau_n=t_n-t_{n-1}$, $r_n=\tau_n/\tau_{n-1}$, tildes for the unfiltered predictor) is only a reindexing:

$$
\frac{\tilde u^{\,n}-u^{n-1}}{\tau_n}=f(t_n,\tilde u^{\,n}),
\qquad
u^{n}=\tilde u^{\,n}-\frac{r_n(1+r_n)}{1+2r_n}
\Bigl(\frac{1}{1+r_n}\tilde u^{\,n}-u^{n-1}+\frac{r_n}{1+r_n}u^{n-2}\Bigr).
$$

Consistency check: at $r_n=1$ the filter degenerates to Guzel-Layton's classical form

$$
u^{n}=\tilde u^{\,n}-\tfrac13\bigl(\tilde u^{\,n}-2u^{n-1}+u^{n-2}\bigr),
$$

as direct computation confirms.

The filter term is in fact a rescaled **second divided difference**:

$$
u^{n}=\tilde u^{\,n}-\eta^{(2)}\,\delta^{2}\tilde u,
\qquad
\eta^{(2)}=\frac{k_{n+1}}{\dfrac{1}{k_{n+1}}+\dfrac{1}{k_{n+1}+k_{n}}},
\qquad
\delta^{2}\tilde u=\frac{\dfrac{\tilde u^{\,n}-u^{n-1}}{k_{n+1}}
-\dfrac{u^{n-1}-u^{n-2}}{k_{n}}}{k_{n+1}+k_{n}} .
$$

Eliminating the intermediate value $\tilde u$ turns FiBE into the equivalent one-leg multistep method

$$
\sum_{j=1}^{2}\Bigl[\prod_{i=1}^{j-1}(t_{n+m}-t_{n+m-i})\Bigr]\delta^{j}y
\;=\;f\Bigl(t_{n+m},\;y_{n+m}
+\frac{\eta^{(2)}}{1-\eta^{(2)}c^{(2)}_{m}}\,\delta^{2}y\Bigr),
\qquad m=2 .
$$

Its left-hand side is precisely the variable-step BDF2 difference operator, and that is the structural fact paper 69 exploits.

### The results

The paper establishes a **discrete gradient structure** for the one-leg multistep scheme corresponding to variable-step FiBE, and from it a discrete energy dissipation law in the dissipative case, together with $L^2$ stability and an $L^2$ error estimate under the step-ratio condition below. It claims priority: this appears to be the **first** energy stability and $L^2$-norm error estimate for a variable-step time-filtered stiff solver. The analytical apparatus is **two new classes of discrete orthogonal convolution kernel** — a plural that matches the OLM structure, since the BDF2 convolution on the left and the one-leg/interpolation convolution on the right each need their own orthogonal dual.

The equations treated are **linear parabolic**, not Allen-Cahn, so this paper carries no maximum-principle result.

### The step-ratio condition

The stability and $L^2$ error estimates hold under the **two-sided "practical" step-ratio constraint**

$$
\tfrac12\ \le\ \frac{\tau_k}{\tau_{k-1}}\ \le\ 2
\qquad(k\ge2).
$$

> [!warning] $[1/2,\,2]$ is not the same kind of object as the other thresholds
> **First, it is two-sided.** Every other result in this group restricts only the **upper** step ratio ($r_k<1+\sqrt2$, $r_k<3.561$, $r_k<1.4877$), whereas here the ratio is also bounded below by $1/2$ — the steps may neither grow nor **shrink** too abruptly. This reflects the loss of A-stability for variable-step FiBE.
>
> **Second, $[1/2,2]$ is not a sharp analytic threshold.** $1+\sqrt2$ is a root of $r^2-2r-1=0$, $3.561$ is a root of $r^2-3r-2=0$, and $1.4877$ is the unique positive root of an explicit equation; $[1/2,2]$ is the empirical safeguard interval standard in adaptive codes — DeCaria-Guzel-Layton-Li restrict their adaptive experiments to a maximum of two and a minimum of one half, a common empirical practice in variable-step methods. That is why the paper itself calls it a practical constraint.
>
> **It should therefore not be quoted alongside $1+\sqrt2$, $3.561$, $4.8645$ or $1.4877$ as though they were quantities of the same type.**

### Relation to the others

It is the **linear parabolic, time-filtering** branch of the same programme: the same discrete-gradient-structure-plus-DOC toolkit as papers 48, 52, 58 and 67, applied to a **predictor-corrector** method instead of to a BDF formula itself. The reference list confirms the lineage directly: it cites paper 48 (_SIAM J. Numer. Anal._ 58:2294-2314), paper 58 (_CSIAM Trans. Appl. Math._ 3:318-334) and paper 67 (_J. Comput. Math._, DOI `10.4208/jcm.2207-m2022-0020`), along with Liao-Zhang (_Math. Comp._ 90:1207-1226), Li-Liao (_SIAM J. Numer. Anal._ 60:2253-2272), Liao-Ji-Wang-Zhang (_J. Sci. Comput._ 92:52, the source of $4.8645$) and Liao-Ji-Zhang (_IMA J. Numer. Anal._ 42:649-679).

**The one-leg reformulation is the bridge to the BDF2 papers**: the OLM corresponding to FiBE has exactly the variable-step BDF2 operator on its left, so the gradient-structure machinery of papers 48 and 52 transfers directly, with the nonlinearity and forcing evaluated at a filtered argument. Its step-ratio conclusion is the **outlier** of the group: two-sided and empirical rather than a one-sided sharp algebraic threshold, reflecting exactly the loss of A-stability under variable steps.

Downstream, the follow-up paper _Energy dissipation laws of time filtered BDF methods up to fourth-order for the molecular beam epitaxial equation_ explicitly builds on this one to extend the theory to FiBDF-$k$ and to the MBE model; Wang-Liao-Zhao (_Numer. Math. Theor. Meth. Appl._ 16(1) (2023) 165-181) treats **constant-step** filtered backward Euler for MBE with slope selection.

## 74: making the tool itself the object of study

### The idea

Each of the five preceding papers has to answer the same question: when is the real quadratic form

$$
\sum_{k=1}^{n}w_k\sum_{j=1}^{k}a^{(k)}_{k-j}w_j
$$

built from some family of convolution kernels positive definite? Each proves it afresh for its own family. This paper solves the question **as an object in its own right**.

On a **uniform** grid there is a complete classical answer: by the Toeplitz-Carathéodory theorem, if $\hat a(z)=\sum_{k\ge0}a_kz^k$ is analytic on the open unit disc $\mathcal D_z$, then the form is positive semidefinite **if and only if** $\mathrm{Re}[\hat a(z)]\ge0$ on $\mathcal D_z$; López-Marcos (1990, Prop. 5.2) distilled from it the usable sufficient condition

$$
a_j\ge0,\qquad a_{j-1}\ge a_j,\qquad a_{j-1}-a_j\ge a_j-a_{j+1} .
$$

On a **variable** grid the kernel acquires the extra index $a^{(n)}_{n-k}$ (its value depends on the current time $t_n$), the generating-function machinery collapses, and until now there were only positivity-preserving constructions for special kernels. This paper supplies the missing general criterion.

**The criterion is purely algebraic, with no explicit connection to the underlying continuous kernel — that is the point.** The proof strategy is clean too: construct both DOC and DCC kernels for a general family, deduce monotonicity of the DCC kernels from the conditions, and close with the equivalence "DOC is positive definite exactly when the original kernel is".

### Setting

For variable-step kernels $\{a^{(n)}_{n-k}\}_{k=1}^{n}$, examine the sign of the form above for arbitrary $\{w_1,\dots,w_n\}$. The structural assumption is $a^{(n)}_j\ne0$ for $0\le j\le n-1$, and if $a^{(n)}_{n_0}=0$ for some $2\le n_0\le n-1$ then $a^{(n)}_j=0$ for all $j\ge n_0$.

### Derivation

**The main theorem (Theorem 1.1).** For fixed $n\ge2$, if

$$
\textbf{C1}:\ a^{(n)}_{j}>0\quad(0\le j\le n-1);
\qquad
\textbf{C2}:\ a^{(n-1)}_{j-1}>a^{(n)}_{j}\quad(1\le j\le n-1);
$$

$$
\textbf{C3}:\ a^{(n-1)}_{j-1}a^{(n)}_{j+1}\ \ge\ a^{(n-1)}_{j}a^{(n)}_{j}\quad(1\le j\le n-2);
\qquad
\textbf{C4}:\ a^{(n)}_{j-1}\ \ge\ a^{(n)}_{j}\quad(1\le j\le n-1),
$$

then the discrete convolution kernel $a^{(n)}_{n-k}$ is **positive definite**.

What the four say: C1 is positivity; C4 is monotonicity in the lag index within a fixed level; C2 is monotonicity **across levels** (trivial on a uniform grid, and precisely the new demand created by variable steps); and C3 is a cross-level inequality with a log-convexity or total-positivity flavour, playing the role that López-Marcos's convexity condition $a_{j-1}-a_j\ge a_j-a_{j+1}$ plays on a uniform grid.

**Tool one: DOC kernels.** This paper is the canonical source for the **general variable-step** definition:

$$
\theta_0^{(n)}:=\frac{1}{a^{(n)}_0},
\qquad
\theta_{n-k}^{(n)}:=-\frac{1}{a^{(k)}_0}\sum_{j=k+1}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k},
\qquad k=n-1,n-2,\dots,1,
$$

satisfying the discrete orthogonality identity $\sum_{j=k}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}\equiv\delta_{nk}$ for $1\le k\le n$. Lemma 2.1 supplies mutual orthogonality,

$$
\sum_{j=k}^{n}a^{(n)}_{n-j}\theta_{j-k}^{(j)}=\delta_{nk}
\qquad\text{and}\qquad
\sum_{j=k}^{n}\theta_{n-j}^{(n)}a^{(j)}_{j-k}=\delta_{nk},
\qquad 1\le k\le n,
$$

whence **$\theta_{n-j}^{(n)}$ is positive (semi-)definite if and only if $a^{(n)}_{n-k}$ is**. The original DOC construction is due to Liao-Zhang (there for the $L^2$ stability of variable-step BDF2).

**Tool two: DCC kernels.** This paper defines them by **summing the DOC kernels along the level index**:

$$
p_{n-k}^{(n)}:=\sum_{j=k}^{n}\theta_{j-k}^{(j)},
\qquad 1\le k\le n,
$$

so that

$$
\theta_0^{(n)}=p_0^{(n)},
\qquad
\theta_{n-k}^{(n)}=p_{n-k}^{(n)}-p_{n-k-1}^{(n-1)}\quad(1\le k\le n-1),
$$

and they satisfy the discrete complementarity identity

$$
\sum_{j=k}^{n}p_{n-j}^{(n)}a^{(j)}_{j-k}\equiv 1,
\qquad 1\le k\le n .
$$

The proof sets $\Xi_k^{(n)}:=\sum_{j=k}^{n}p^{(n)}_{n-j}a^{(j)}_{j-k}$, verifies $\Xi_k^{(k)}=1$ and $\Xi_k^{(n)}=\Xi_k^{(n-1)}$, then inducts. **The two names are therefore literal: DOC inverts the convolution (giving $\delta_{nk}$), DCC complements it (giving $1$).** Closed forms for both (Lemmas 2.3 and 2.6) are expressed through auxiliary sequences $\psi^{(m)}_j$ and $\chi^{(k)}_\ell$:

$$
\theta_{j-k}^{(j)}=-\frac{1}{a_0^{(j)}}\psi_1^{(k+1)}
\prod_{\ell=k+2}^{j}\bigl(\chi_2^{(j-\ell)}\psi_2^{(\ell)}-\psi_1^{(\ell)}\bigr)
\quad(j\ge k+1),
$$

$$
p_{n-k}^{(n)}=\frac{1}{a_0^{(k)}}-\psi_1^{(k+1)}\sum_{j=k+1}^{n}\frac{1}{a_0^{(j)}}
\prod_{\ell=k+2}^{j}\bigl(\chi_2^{(j-\ell)}\psi_2^{(\ell)}-\psi_1^{(\ell)}\bigr).
$$

**The engine of the proof is DCC monotonicity (Lemma 2.7).** For $n\ge2$: if $a^{(n)}_{n-k}$ satisfies **C1-C3**, then

$$
p_0^{(n)}>0,
\qquad
p_0^{(n-1)}>p_1^{(n)},
\qquad
p_{j-1}^{(n-1)}\ge p_j^{(n)}\ \ (2\le j\le n-1);
$$

**C4** then supplies the remaining monotonicity, and $\theta$ is positive definite. Theorem 1.1 is DCC monotonicity combined with the DOC equivalence of Lemma 2.1.

**Step-ratio restriction: none.** C1-C4 are conditions on the **kernel** and hold **on arbitrary nonuniform grids**, imposing nothing on $r_k=\tau_k/\tau_{k-1}$. That is one of the paper's selling points, and it is the watershed between this route, paper 58's (Grenander-Szegő, uniform-grid Toeplitz) and paper 67's (explicit discrete gradient structure, variable-step BDF3): of the three, only this one survives an arbitrary grid.

### Theorems

**Theorem 1.1** (C1-C4 $\Rightarrow$ positive definiteness) as above. The paper calls it the first result with simple algebraic conditions for positive definiteness of variable-step convolution coefficients. Four applications follow.

**Application one: the variable-step L1 formula (Proposition 4.1).** With $\omega_\gamma(t):=t^{\gamma-1}/\Gamma(\gamma)$, the Riemann-Liouville integral $(\mathcal I^\gamma v)(t)=\int_0^t\omega_\gamma(t-s)v(s)\mathrm ds$ and the Caputo derivative $(\partial_t^\alpha v)(t)=(\mathcal I^{1-\alpha}v')(t)$ for $0<\alpha<1$, the **L1 formula** on a nonuniform grid replaces $v'$ on $(t_{k-1},t_k)$ by the constant $\nabla_\tau v^k/\tau_k$:

$$
(\partial_\tau^\alpha v)^n:=\sum_{k=1}^{n}c^{(n,\alpha)}_{n-k}\nabla_\tau v^{k},
\qquad
c^{(n,\alpha)}_{n-k}:=\frac1{\tau_k}\int_{t_{k-1}}^{t_k}\omega_{1-\alpha}(t_n-s)\,\mathrm ds .
$$

Then

$$
c^{(n,\alpha)}_{j}>0,
\quad
c^{(n,\alpha)}_{j-1}>c^{(n,\alpha)}_{j},
\quad
c^{(n-1,\alpha)}_{j-1}>c^{(n,\alpha)}_{j},
\quad
c^{(n-1,\alpha)}_{j-1}c^{(n,\alpha)}_{j+1}>c^{(n-1,\alpha)}_{j}c^{(n,\alpha)}_{j},
$$

so C1-C4 all hold and **the variable-step L1 kernels are positive definite on any nonuniform grid**. The proof is elegant enough to record: the integral mean value theorem gives $c^{(n,\alpha)}_{n-k}=\omega_{1-\alpha}(t_n-s_{nk})$ for some $s_{nk}\in[t_{k-1},t_k]$, which yields C1 and C4 immediately — **the discrete kernel is literally a sampling of the continuous kernel $\omega_{1-\alpha}$**. Then applying the Cauchy mean value theorem to $c_{n,k}(\mu):=\frac1{\tau_k}\int_{t_{k-1}}^{t_{k-1}+\mu\tau_k}\omega_{1-\alpha}(t_n-s)\mathrm ds$ gives

$$
\psi^{(n,\alpha)}_{n-k}:=\frac{c^{(n,\alpha)}_{n-k}}{c^{(n-1,\alpha)}_{n-1-k}}
=\Bigl(\frac{t_{n-1}-t_{k-1}-\xi_{1k}\tau_k}{t_n-t_{k-1}-\xi_{1k}\tau_k}\Bigr)^{\alpha},
\qquad \xi_{1k}\in(0,1),
$$

hence

$$
\Bigl(\frac{t_{n-1}-t_k}{t_n-t_k}\Bigr)^{\alpha}<\psi^{(n,\alpha)}_{n-k}
<\Bigl(\frac{t_{n-1}-t_{k-1}}{t_n-t_{k-1}}\Bigr)^{\alpha},
\qquad
0<\psi_1^{(n,\alpha)}<\psi_2^{(n,\alpha)}<\cdots<\psi_{n-1}^{(n,\alpha)}<(t_{n-1}/t_n)^{\alpha}<1 .
$$

That $\psi^{(n,\alpha)}_j$ increases in $j$ is exactly C2 and C3.

**(Sharpness, Remark 4.1)** C1-C4 are sufficient but **not necessary**. The second-order **L1$^{+}$** formula of Ji-Liao-Gong-Zhang, with kernels

$$
\bar c_{n-k}^{(n)}:=\frac{1}{\tau_n\tau_k}\int_{t_{n-1}}^{t_n}
\int_{t_{k-1}}^{\min\{t,t_k\}}\omega_{1-\alpha}(t-s)\,\mathrm ds\,\mathrm dt,
$$

**is** positive semidefinite (inherited from the continuous kernel), yet the condition $\bar c_0^{(n)}\ge\bar c_1^{(n)}$ fails for some $n$. Finding necessary and sufficient algebraic conditions remains open.

**Application two: energy stability for the time-fractional Allen-Cahn equation (Proposition 4.2).** For $\partial_t^\alpha u=\varepsilon^2\Delta u-F'(u)$ with $F(u)=\frac14(1-u^2)^2$, the continuous energy law $E[u(t)]\le E[u(0)]$ is due to Tang and to Tang-Yu-Zhou, that is [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|paper 40]]. Positive definiteness of the L1 kernels gives

$$
\sum_{k=1}^{n}\nabla_\tau v^{k}\,(\partial_\tau^\alpha v)^{k}
=\sum_{k=1}^{n}\nabla_\tau v^{k}\sum_{j=1}^{k}c^{(k,\alpha)}_{k-j}\nabla_\tau v^{j}>0
\qquad(\nabla_\tau v^k\not\equiv0).
$$

Applied to Ji-Liao-Zhang's first-order stabilised semi-implicit scheme

$$
(\partial_\tau^\alpha u)^n=\varepsilon^2D_hu^n-F'(u^{n-1})-S(u^n-u^{n-1}),
\qquad n\ge1,
$$

the conclusion is: **if the stabilisation parameter satisfies $S\ge2$**, the scheme preserves the discrete maximum principle and $E_h^n\le E_h^0$ for $n\ge1$, where

$$
E_h^n=-\frac{\varepsilon^2}{2}(u^n)^TD_hu^n+\sum_{x_h\in\Omega_h}F(u_h^n).
$$

The paper notes that the uniform-grid counterpart was known (Ji-Liao-Zhang; Tang-Yu-Zhou) and that this **nonuniform-grid** version appears to be new.

> [!note] Three energy laws, three different statements
> Paper 40 proves a **fractional (nonlocal, integral-type)** energy law; [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|paper 57]] proves a differential-type law for a **variational** energy; the result here is about the **discrete** energy, $E_h^n\le E_h^0$. Their objects and their forms all differ, and citing one in place of another is a mistake.

**Application three: Riemann-Liouville integrals and fractional wave equations.** For $\partial_t u=\mathcal I^\gamma\Delta u+f$ (whose solution has $\partial_{tt}u\sim\mathcal O(t^{\gamma-1})$ near $t=0$, making graded meshes the natural remedy), the backward Euler scheme $\partial_\tau u^{n-\frac12}=(\mathcal I^\gamma_\tau\Delta u)^n+f(x,t_n)$ with the midpoint quadrature $(\mathcal I^\gamma_\tau v)^n=\sum_{k=1}^{n}c^{(n,1-\gamma)}_{n-k}\tau_kv^{k-\frac12}$ satisfies, on a **general** nonuniform grid,

$$
\|u^n\|_{L^2}\le\|u^0\|_{L^2}+\sum_{k=1}^{n}\tau_k\|f(t_k)\|_{L^2}.
$$

**Application four: weakly singular Volterra equations (Proposition 4.3).** For $\partial_t u=\mathcal K_t^{(\beta)}\Delta u+f$ with $(\mathcal K^{(\beta)}_tv)(t)=\int_0^t\kappa_\beta(t-s)v(s)\mathrm ds$ and midpoint kernels $\kappa^{(n,\beta)}_{n-k}:=\frac1{\tau_k}\int_{t_{k-1}}^{t_k}\kappa_\beta(t_n-s)\mathrm ds$: **if $\kappa_\beta>0$, $\kappa_\beta'<0$ and $\kappa_\beta''\ge0$**, then $\kappa^{(n,\beta)}_{n-k}$ satisfies C1-C4 (the argument uses that $y=\kappa_\beta(t_n-x)/\kappa_\beta(t_{n-1}-x)$ decreases on $t_{n-1}<x<t_n$) and is therefore positive definite.

**Together the four applications show the value of isolating the tool**: one criterion covers the fractional Caputo derivative, the Riemann-Liouville integral and general weakly singular kernels at once, and **without any restriction on the mesh**.

### Numerical experiments

**This paper has no numerical experiments.** It is pure analysis, and its only figure is a schematic of the kernels. Every conclusion it offers is a theorem with no measured data to compare against — the same situation as paper 58.

### Relation to the others

It is the **general theory underneath** the fractional half of this collection. The positive definiteness of variable-step L1 kernels proved here is the property used in various guises by [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|paper 43]] (Alikhanov/L1-type kernels for the time-fractional Allen-Cahn equation on nonuniform grids) and [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|paper 57]] (the L1$_R$ kernels of the Riemann-Liouville reformulation) — paper 57 invokes the three algebraic criteria from here directly when proving monotonicity of its DOC kernels, so this paper precedes paper 57 logically even though it was published later. The DOC construction it formalises is the same one papers 48 and 52 use at the BDF2 level, paper 58 at the BDF-$k$ level and paper 67 at the BDF3 level; the DCC construction it formalises is the source of paper 43's complementary-kernel bounds. **This is the paper that states the general variable-step definitions and the two identities (DOC's $\equiv\delta_{nk}$ and DCC's $\equiv1$) cleanly for an arbitrary kernel family.**

**The point of singling out paper 74** is that the technical core of this whole series is not any one scheme but the algebraic question of when the real quadratic form collected from the history terms is positive definite. Once that question has independent criteria, the same argument transfers to third-order BDF, time-filtered Euler, fractional L1 approximations and implicit-explicit Runge-Kutta methods.

## How the six relate

| No. | Object                             | Mesh      | Threshold                                  | Core tool                                          | Numerical experiments                                            |
| --- | ---------------------------------- | --------- | ------------------------------------------ | -------------------------------------------------- | ---------------------------------------------------------------- |
| 48  | Allen-Cahn, BDF2                   | variable  | S1 $3.561$; S0 $1+\sqrt2$                  | modified energy; KRC; DCC                          | second order on random grids; bubble merging; coarsening         |
| 52  | MBE without slope selection, BDF2  | variable  | $3.561$ (same for energy and $L^2$)        | DOC kernels; $\mathcal M_r$                        | random grid: still second order at $\max r_k=850$                |
| 58  | linear reaction-diffusion, BDF-$k$ | uniform   | none ($3\le k\le5$)                        | DOC kernels; Grenander-Szegő                       | none; only figures illustrating the lemmas                       |
| 67  | diffusion, BDF3                    | variable  | $1.4877$ (sufficient) / $1.69$ (necessary) | discrete gradient structure $G$; variable-step DOC | third order; stable with half the levels violating the threshold |
| 69  | linear parabolic, filtered Euler   | variable  | $[1/2,2]$ (two-sided, heuristic)           | one-leg reformulation; two new DOC classes         | —                                                                |
| 74  | the quadratic form itself          | arbitrary | none (supplies algebraic criteria C1-C4)   | general theory of DOC and DCC                      | none                                                             |

**How to read this table**: the numbers in the threshold column cannot be compared across rows, since each is a condition for a different scheme and a different conclusion. What is genuinely comparable is the last column: **wherever a paper ran random-grid experiments, the observed robust range was far wider than the theorems guarantee** (paper 52's $850.80$, paper 67's half the levels in violation). That is the most consistent empirical fact across the six. None of them reports a test that breaks its own conclusion; the one experiment that does break the maximum bound is in [[en/computational-mathematics/paper-notes/phase-field-and-time-stepping/time-fractional-phase-field|paper 43]], and what it violates there is the **step-size-to-mesh coupling condition**, not a step-ratio condition. The contrast is the point: the coupling condition comes far closer to being necessary than the step-ratio conditions do.

## Sources

- H.-l. Liao, T. Tang, and T. Zhou, [_On energy stable, maximum-principle preserving, second-order BDF scheme with variable steps for the Allen-Cahn equation_](https://doi.org/10.1137/19M1289157), SIAM J. Numer. Anal. 58(4) (2020), pp. 2294-2314 (preprint [arXiv:2003.00421](https://arxiv.org/abs/2003.00421)).
- H.-l. Liao, X. Song, T. Tang, and T. Zhou, [_Analysis of the second-order BDF scheme with variable steps for the molecular beam epitaxial model without slope selection_](https://doi.org/10.1007/s11425-020-1817-4), Sci. China Math. 64 (2021), pp. 887-902 (preprint [arXiv:2008.03185](https://arxiv.org/abs/2008.03185)).
- H.-l. Liao, T. Tang, and T. Zhou, [_A new discrete energy technique for multi-step backward difference formulas_](https://doi.org/10.4208/csiam-am.SO-2021-0032), CSIAM Trans. Appl. Math. 3 (2022), pp. 318-334 (preprint [arXiv:2102.04644](https://arxiv.org/abs/2102.04644)).
- H.-l. Liao, T. Tang, and T. Zhou, [_Discrete energy analysis of the third-order variable-step BDF time-stepping for diffusion equations_](https://doi.org/10.4208/jcm.2207-m2022-0020), J. Comput. Math. 41 (2023), pp. 325-344 (preprint [arXiv:2204.12742](https://arxiv.org/abs/2204.12742)).
- H.-l. Liao, T. Tang, and T. Zhou, [_Stability and convergence of the variable-step time filtered backward Euler scheme for parabolic equations_](https://doi.org/10.1007/s10543-023-00982-y), BIT Numer. Math. 63 (2023), 39.
- H.-l. Liao, T. Tang, and T. Zhou, [_Positive definiteness of real quadratic forms resulting from the variable-step L1-type approximations of convolution operators_](https://doi.org/10.1007/s11425-022-2229-5), Sci. China Math. 67 (2024), pp. 237-252 (preprint [arXiv:2011.13383](https://arxiv.org/abs/2011.13383)).

Method sources and background (cited here but outside this topic):

- V. DeCaria, S. Guzel, W. Layton, and Y. Li, [_A variable stepsize, variable order family of low complexity_](https://doi.org/10.1137/19M1291666), SIAM J. Sci. Comput. 43(3) (2021), pp. A2130-A2160 (preprint [arXiv:1810.06670](https://arxiv.org/abs/1810.06670)) — the source of the filtered scheme analysed in paper 69.
- H.-l. Liao, B. Ji, L. Wang, and Z. Zhang, [_Mesh-robustness of an energy stable BDF2 scheme with variable steps for the Cahn-Hilliard model_](https://doi.org/10.1007/s10915-022-01923-7), J. Sci. Comput. 92 (2022), 52 (preprint [arXiv:2102.03731](https://arxiv.org/abs/2102.03731)) — the source of $4.864$.
