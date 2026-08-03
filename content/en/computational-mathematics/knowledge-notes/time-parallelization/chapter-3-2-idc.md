---
title: "3.3: Parallel Integral Deferred Correction"
description: Deriving the IDC update from the integral residual, explaining why deferred correction raises the order, and giving a full account of the parallel scheduling and regularity limits of PIDC and RIDC
lang: en
translation: computational-mathematics/knowledge-notes/time-parallelization/chapter-3-2-idc
tags:
  - parallel-in-time
  - IDC
  - high-order-time-integration
---

> [!note] Reading scope
> This page corresponds to Section 3.3 of the paper (pp. 405–412) and covers equations (3.5)–(3.12), Theorem 3.3, and Figures 3.4–3.6. The derivation keeps every algebraic step, and the parallel schedules of PIDC and RIDC are explained separately. IDC was proposed by Dutt, Greengard, and Rokhlin (2000); PIDC was introduced by Guibert and Tromeur-Dervout (2007); RIDC was proposed by Christlieb, Macdonald, and Ong (2010).

## 3.3 Time-parallel IDC

Integral deferred correction (IDC), proposed by Dutt, Greengard, and Rokhlin (2000), uses round after round of iterative correction to gradually lift a low-order time integrator into a high-order method. The original IDC is sequential in time; this section presents two time-parallel reworkings of it: the pipeline IDC (PIDC) of Guibert and Tromeur-Dervout (2007) and the revisionist IDC (RIDC) of Christlieb et al. (2010). Both differ substantially from the original IDC in parallel structure, so we must first make clear how IDC itself works.

> [!tip] Insight
> The idea of deferred correction can be traced further back to Böhmer and Stetter (1984) and their iterative correction of discretization error; Dutt et al. (2000) recast it as residual correction on the **integral equation**, obtaining the numerically stable and easily order-raising spectral deferred correction (SDC). The key to understanding IDC is to view it as a Picard-type fixed-point iteration: each round uses the previous round's trajectory to estimate a residual once, then solves an error equation of the same structure but with a more accurate right-hand side.

### Where IDC begins

Consider a nonlinear initial-value problem. Rather than stepping the differential equation directly, IDC first writes the ODE in its equivalent integral form:

$$
\boldsymbol u(t)=\boldsymbol u_0+
\int_0^t\boldsymbol f(\boldsymbol u(\tau),\tau)\,d\tau,
\qquad t\in(0,T]. \tag{3.5}
$$

There are two reasons to choose the integral form over the differential form. First, the integral operator is **smoothing**: integrating $\boldsymbol f$ is more stable than differentiating it, and errors are not amplified. Second, once a trajectory approximation is available, the right-hand integral can be evaluated accurately with a high-order quadrature rule, which is precisely the source of the order increase later.

Start with a rough approximation $\widetilde{\boldsymbol u}(t)$. It may be the constant function $\boldsymbol u_0$ (i.e. $\widetilde{\boldsymbol u}(t)\equiv\boldsymbol u_0$) or the output of a low-accuracy integrator. Define the true error and the integral residual by

$$
\boldsymbol e(t)=\boldsymbol u(t)-\widetilde{\boldsymbol u}(t),
$$

$$
\boldsymbol r(t)=\boldsymbol u_0+
\int_0^t\boldsymbol f(\widetilde{\boldsymbol u}(\tau),\tau)\,d\tau
-\widetilde{\boldsymbol u}(t). \tag{3.6}
$$

The residual $\boldsymbol r$ measures the extent to which the current trajectory $\widetilde{\boldsymbol u}$ violates the integral equation (3.5): how much the result of substituting $\widetilde{\boldsymbol u}$ into the right-hand side of the integral equation differs from $\widetilde{\boldsymbol u}$ itself. When $\widetilde{\boldsymbol u}$ is exactly the true solution, $\boldsymbol r\equiv 0$; hence the smaller $\boldsymbol r$ is, the closer the trajectory is to the true solution. Note that $\boldsymbol r(0)=\boldsymbol u_0-\widetilde{\boldsymbol u}(0)=\boldsymbol 0$, and the error likewise satisfies $\boldsymbol e(0)=\boldsymbol 0$; the initial values of both are fixed.

Substituting $\boldsymbol u=\widetilde{\boldsymbol u}+\boldsymbol e$ into (3.5), then adding and subtracting $\int_0^t\boldsymbol f(\widetilde{\boldsymbol u},\tau)\,d\tau$, gives

$$
\begin{aligned}
\boldsymbol e(t)
&=\boldsymbol u_0+
\int_0^t\boldsymbol f(\widetilde{\boldsymbol u}(\tau)+\boldsymbol e(\tau),\tau)\,d\tau
-\widetilde{\boldsymbol u}(t)\\
&=\boldsymbol r(t)+
\int_0^t\left[
\boldsymbol f(\widetilde{\boldsymbol u}(\tau)+\boldsymbol e(\tau),\tau)
-\boldsymbol f(\widetilde{\boldsymbol u}(\tau),\tau)
\right]d\tau. \tag{3.7}
\end{aligned}
$$

The first line simply substitutes $\boldsymbol u=\widetilde{\boldsymbol u}+\boldsymbol e$ into the integral equation; the second line extracts $\int_0^t\boldsymbol f(\widetilde{\boldsymbol u},\tau)\,d\tau$ from the integral, which is exactly the definition (3.6) of the residual $\boldsymbol r(t)$, and what remains is the **increment** of the integrand with respect to $\boldsymbol e$. This step expresses the error $\boldsymbol e$ explicitly as "the known residual $\boldsymbol r$" plus "an integral equation in $\boldsymbol e$".

Differentiating in time gives the error equation

$$
\boldsymbol e'(t)-\boldsymbol r'(t)
=\boldsymbol f(\widetilde{\boldsymbol u}(t)+\boldsymbol e(t),t)
-\boldsymbol f(\widetilde{\boldsymbol u}(t),t). \tag{3.8}
$$

Equation (3.7) is the integral form of the error equation, and (3.8) is its differential form. Both forms are kept because the next discretization step **treats the two terms differently**: the residual term $\boldsymbol r'$ (equivalent to the integral increment of $\boldsymbol r$) is handled with a high-order quadrature rule, while the $\boldsymbol f$-increment term of the error equation can still use the simplest low-order time stepper.

> [!tip] Insight
> The right-hand side of (3.8) is $\boldsymbol f(\widetilde{\boldsymbol u}+\boldsymbol e,\cdot)-\boldsymbol f(\widetilde{\boldsymbol u},\cdot)$, and when $\boldsymbol e$ is already small it is approximately $\partial_{\boldsymbol u}\boldsymbol f\cdot\boldsymbol e$, a "linearized" equation in the error; its inhomogeneous term comes from the residual $\boldsymbol r$. This explains why the error equation can be solved with the same cheap stepper as the original equation—its solution is itself a small quantity, so the **absolute** error of the low-order stepper is small too. All the credit for raising the order goes to the high-order quadrature of the residual integral, not to the stepper.

### The discrete correction formula

Take nodes

$$
0=t_0<t_1<\cdots<t_M=T,
\qquad \boldsymbol u_m^k\approx\boldsymbol u(t_m),
$$

where the superscript $k$ is the correction round (sweep) and the subscript $m$ is the node index; $\boldsymbol u_m^k$ denotes the approximation at node $t_m$ in round $k$. Let $\Delta t_m=t_{m+1}-t_m$. Applying the linear $\theta$ method to (3.8) (i.e. weighting $\boldsymbol f$ at the two ends of the interval with a parameter $\theta\in[0,1]$) gives

$$
\begin{aligned}
\boldsymbol e_{m+1}-\boldsymbol e_m
={}&\boldsymbol r_{m+1}-\boldsymbol r_m\\
&+\Delta t_m(1-\theta)
\left[\boldsymbol f(\boldsymbol u_m^{k+1},t_m)
-\boldsymbol f(\boldsymbol u_m^k,t_m)\right]\\
&+\Delta t_m\theta
\left[\boldsymbol f(\boldsymbol u_{m+1}^{k+1},t_{m+1})
-\boldsymbol f(\boldsymbol u_{m+1}^k,t_{m+1})\right]. \tag{3.9}
\end{aligned}
$$

Here the new trajectory is written $\boldsymbol u^{k+1}=\boldsymbol u^k+\boldsymbol e$, and the $\widetilde{\boldsymbol u}+\boldsymbol e$ and $\widetilde{\boldsymbol u}$ on the right-hand side of (3.8) correspond to $\boldsymbol u^{k+1}$ and $\boldsymbol u^k$ respectively. The linear $\theta$ method uses the difference quotient $(\boldsymbol e_{m+1}-\boldsymbol e_m)/\Delta t_m$ for the left-hand $\boldsymbol e'$ and weights the right-hand side by $1-\theta$ at $t_m$ and by $\theta$ at $t_{m+1}$. Taking $\theta=1$ gives backward Euler, $\theta=0$ gives forward Euler, and $\theta=1/2$ gives the trapezoidal rule.

From the residual definition (3.6), the residual increment between two adjacent nodes is

$$
\boldsymbol r_{m+1}-\boldsymbol r_m
=\int_{t_m}^{t_{m+1}}\boldsymbol f(\boldsymbol u^k(\tau),\tau)\,d\tau
-\left(\boldsymbol u_{m+1}^k-\boldsymbol u_m^k\right).
$$

That is, the residual increment equals the **exact integral** of the old trajectory over the small interval $[t_m,t_{m+1}]$ minus the old trajectory's own endpoint jump. The former is approximated by high-order quadrature, and the latter is a known quantity, so the residual increment is fully computable. The local integral is evaluated by Lagrange interpolatory quadrature over the full set of nodes:

$$
\int_{t_m}^{t_{m+1}}\boldsymbol f(\boldsymbol u^k(\tau),\tau)\,d\tau
\approx\sum_{j=1}^{M}\omega_{m,j}
\boldsymbol f(\boldsymbol u_j^k,t_j), \tag{3.10a}
$$

$$
\omega_{m,j}=
\int_{t_m}^{t_{m+1}}
\prod_{\substack{i=1\\i\ne j}}^{M}
\frac{\tau-t_i}{t_j-t_i}\,d\tau. \tag{3.10b}
$$

The crucial point is this: the index $j$ in the sum (3.10a) ranges over **all** $M$ nodes, not just the two endpoints of the small interval $[t_m,t_{m+1}]$. The weights $\omega_{m,j}$ are obtained by integrating the Lagrange basis functions over the whole span $[0,T]$ on the small interval, so (3.10a) is a high-order quadrature rule that can reach order $M$. This is exactly the step that lets the cheap low-order stepper "see" the curvature information of the entire interval.

Then use $\boldsymbol u_m^{k+1}=\boldsymbol u_m^k+\boldsymbol e_m$ to eliminate the error variable $\boldsymbol e$: replace $\boldsymbol e_{m+1}-\boldsymbol e_m$ in (3.9) by $(\boldsymbol u_{m+1}^{k+1}-\boldsymbol u_{m+1}^k)-(\boldsymbol u_m^{k+1}-\boldsymbol u_m^k)$ and substitute the residual increment; the term $\boldsymbol u_{m+1}^k-\boldsymbol u_m^k$ cancels exactly with the identically named term in the residual. The node-by-node correction formula is

$$
\begin{aligned}
\boldsymbol u_{m+1}^{k+1}
={}&\boldsymbol u_m^{k+1}
+\Delta t_m(1-\theta)
\left[\boldsymbol f(\boldsymbol u_m^{k+1},t_m)
-\boldsymbol f(\boldsymbol u_m^k,t_m)\right]\\
&+\Delta t_m\theta
\left[\boldsymbol f(\boldsymbol u_{m+1}^{k+1},t_{m+1})
-\boldsymbol f(\boldsymbol u_{m+1}^k,t_{m+1})\right]\\
&+\sum_{j=1}^{M}\omega_{m,j}
\boldsymbol f(\boldsymbol u_j^k,t_j),
\qquad m=0,\ldots,M-1. \tag{3.11}
\end{aligned}
$$

Each correction level $k=0,1,\ldots,k_{\max}-1$ sweeps from left to right through $m=0,1,\ldots,M-1$. The choice $\theta=1$ gives a backward-Euler correction, and $\theta=1/2$ gives a trapezoidal correction. Note that in (3.11) only $\boldsymbol u_{m+1}^{k+1}$ (which appears inside $\boldsymbol f$ when $\theta\ne 0$) is unknown; everything else comes from the already-computed old trajectory or the already-advanced endpoints of the new trajectory, so each step only requires solving a (possibly nonlinear) algebraic equation of the same size as the base stepper.

> [!tip] Insight (the division of labor among the three terms in (3.11))
> The first line, $\boldsymbol u_m^{k+1}$, advances the new solution forward with the simple integrator; the two middle terms correct the dynamical difference $\boldsymbol f(\boldsymbol u^{k+1},\cdot)-\boldsymbol f(\boldsymbol u^k,\cdot)$ between the new and old trajectories at the local endpoints with weight $\theta$; the final quadrature term $\sum_j\omega_{m,j}\boldsymbol f(\boldsymbol u_j^k,t_j)$ injects high-order integral information from the old trajectory over the entire node set. The reason the correction can raise the order lies in this last term, which sees a more complete time interval than the base stepper—it uses an $M$-th-order-accurate integral as the "target", while the first two terms push the new solution toward that target. If the new and old trajectories coincide, the two middle terms vanish, and the formula reduces to pure deferred correction: "low-order stepping + high-order quadrature residual".

### Theorem 3.3: how much order is gained per round

If the base integrator has order $p$ and $M$ equally spaced nodes are used, the error order after correction $k$ is

$$
O\!\left(\Delta t^{\min\{M,(k+1)p\}}\right).
$$

Its meaning can be read in two parts. The $(k+1)p$ term shows that the initial guess (before $k=0$) reaches order $p$, and thereafter each completed correction round adds another $p$ orders—this comes from the fact that the error equation is itself solved with a $p$-th-order stepper, eliminating the current **leading error term** once more each round. The $\min\{M,\cdot\}$ is the ceiling of the quadrature accuracy: (3.10) uses Lagrange interpolation on $M$ nodes, whose quadrature error is $O(\Delta t^{M})$, so no matter how many rounds are iterated, the overall accuracy can never exceed the order of the quadrature rule itself. Hence backward Euler ($p=1$) raises the order by $1$ per round, the trapezoidal rule ($p=2$) raises it by $2$ per round, until saturation at $M$.

The original IDC of Dutt et al. (2000) uses Gauss-type nodes to raise the quadrature ceiling higher. For example, Gauss–Lobatto nodes can reach order $2J-1$ ($J$ being the number of nodes). This type of IDC with spectral-accuracy quadrature is usually called spectral deferred correction (SDC), and it is also the core component of the PFASST algorithm (see Section 4.3; PFASST was proposed by Emmett and Minion (2012), and its idea originates in Minion (2010), who replaced the fine solver of Parareal with a single step of SDC).

A long time interval is not suited to being approximated as a whole by a single high-order polynomial—a single high-degree polynomial on a large interval exhibits violent oscillations (the Runge phenomenon), and the quadrature accuracy actually degrades. The paper partitions $[0,T]$ into windows

$$
I_n=[T_{n-1},T_n],
\qquad n=1,\ldots,N_t,
$$

where $T_0=0$ and $T_{N_t}=T$. When the windows are small enough, a low-degree polynomial already gives accurate quadrature. Standard IDC first completes all corrections on $I_n$, then passes the endpoint value to $I_{n+1}$. This process is **completely sequential**: the initial value of $I_{n+1}$ is unknown until $I_n$ is finished, so the $(n+1)$-th window must wait for the $n$-th window to end; moreover, the node updates along $m$ within a window also proceed step by step. It is precisely this double sequential dependence that motivates the two parallel reworkings, PIDC and RIDC.

## 3.3.1 Pipeline IDC (PIDC)

### A pipeline laid out across windows

The first parallel IDC, PIDC, was proposed by Guibert and Tromeur-Dervout (2007). It organizes IDC into a pipeline, whose core is a simple observation applicable to any time-evolution computation, already made in Womble (1990): as soon as a **preliminary** endpoint initial value is computed on $I_n$, the downstream window $I_{n+1}=[T_n,T_{n+1}]$ can start immediately, without waiting for $I_n$ to fully converge.

Concretely, after $I_n$ completes its first sweep, a rough endpoint value $\boldsymbol u_{n,M}^{1}$ is obtained at $t=T_n$ (i.e. the rightmost solution after one sweep on window $I_n$). $I_{n+1}$ can immediately use $\boldsymbol u_{n,M}^{1}$ as its initial value to start its own first sweep, while $I_n$ continues with its second sweep. After this step, $I_{n+2}$ can start as well, while $I_{n+1}$ and $I_n$ each advance, giving three sweeps in parallel. In general, while $I_n$ performs sweep $k$, $I_{n+1}$ can simultaneously perform sweep $k-1$, $I_{n+2}$ sweep $k-2$, and so on, down to the first sweep on $I_{n+k-1}$. In this way the windows are staggered along the diagonal, forming a continually widening pipeline.

![Source Figure 3.4: PIDC pipeline startup and steady-state phases on four time windows](assets/papers/time-parallelization/source-figures/figure-3-4.svg)

Figure 3.4 takes $M=6$ and $k_{\max}=4$. Panels (a)–(d) launch the first four windows in turn, and the pipeline width grows from one sweep to four; this is the **fill (bootstrap) phase**: only after computing the preliminary endpoint value of an upstream window does the downstream window have an initial value to use. Panels (e) and (f) show two successive instants after the pipeline is full, when four different correction levels on $I_n$ through $I_{n+3}$ can run simultaneously; this is the **steady-state phase**, with concurrency width reaching $k_{\max}=4$. Black dashed lines record the history of completed sweeps, red lines with red dots indicate the sweeps currently running in parallel, and solid black lines show the exact solution. Thus these six panels contain both startup cost and steady-state concurrency width, and one should not treat (e) and (f) merely as a generic time-stepping sketch—the speedup of the pipeline is discounted by the overhead at the fill/drain ends.

> [!tip] Insight
> The parallelism of PIDC is not free: the steady-state concurrency width is at most $k_{\max}$, the number of correction rounds, not the number of windows $N_t$. When $N_t\gg k_{\max}$, most of the time is spent in steady state and the speedup approaches $k_{\max}$; but the initial value each window receives is a "rough initial value" that is continually rewritten by the upstream window, which brings the convergence hazard below.

Every PIDC window with $n\ge 1$ begins its sweep from a rough and **continually changing** initial value. Because the initial value of a downstream window comes from an upstream solution that has not yet converged, the accuracy of the resulting solution is **not guaranteed to decrease monotonically** as correction rounds advance. This is essentially different from sequential IDC (where each window receives the exact, converged initial value from upstream); it is a phenomenon specific to PIDC and needs regularity experiments to judge its impact.

### The periodic advection–diffusion experiment and the discrete matrix (3.12)

To test the phenomenon above, the paper applies IDC and PIDC to the advection–diffusion equation (2.5), taking two diffusion coefficients $\nu=1$ and $\nu=10^{-3}$. Using periodic boundaries, central-difference discretization, and grid $\Delta x=1/64$, one obtains the linear ODE system $\boldsymbol u'(t)=A\boldsymbol u(t)$, where the semidiscrete matrix is written

$$
A=\frac{\nu}{\Delta x^2}A_{xx}
-\frac{1}{2\Delta x}A_x,
$$

with the first term $\tfrac{\nu}{\Delta x^2}A_{xx}\approx\nu\partial_{xx}$ approximating diffusion and the second term $\tfrac{1}{2\Delta x}A_x\approx\partial_x$ approximating advection. The two difference matrices are

$$
A_{xx}=
\begin{bmatrix}
-2&1&&&1\\
1&-2&1&&\\
&\ddots&\ddots&\ddots&\\
&&1&-2&1\\
1&&&1&-2
\end{bmatrix},
\qquad
A_x=
\begin{bmatrix}
0&1&&&-1\\
-1&0&1&&\\
&\ddots&\ddots&\ddots&\\
&&-1&0&1\\
1&&&-1&0
\end{bmatrix}. \tag{3.12}
$$

The extra elements in the upper-right and lower-left corners of both matrices come from the **periodic boundary**: the first cell and the last cell are neighbors of each other, so the wrap-around terms are supplied cyclically. $A_{xx}$ is symmetric (diffusion is dissipative), and $A_x$ is antisymmetric (advection is conservative).

Take $T=3$, window width $\Delta T=1/10$ (hence $30$ windows in total), $M=5$, and backward Euler as the base integrator. The error on the $n$-th window after sweep $k$ is defined by

$$
\operatorname{err}_n^k=
\frac{\max_m\lVert\boldsymbol u_{\mathrm{ref}}^{n,m}
-\boldsymbol u_k^{n,m}\rVert_\infty}
{\max_{n,m}\lVert\boldsymbol u_{\mathrm{ref}}^{n,m}\rVert_\infty}.
$$

The numerator takes the maximum norm of the difference between the numerical solution and the reference solution over all nodes within window $n$, and the denominator normalizes by the maximum norm of the global reference solution, so this is a **windowwise relative error**. The reference solution $\boldsymbol u_{\mathrm{ref}}$ is computed by MATLAB's built-in ODE45 with both relative and absolute tolerances set to $10^{-13}$, and can be regarded as the "exact solution". The initial trajectory of the next window is fixed to the first-sweep endpoint of the previous window: $\boldsymbol u_{n+1,m}^0\equiv\boldsymbol u_{n,M}^1$ (for all $m=0,1,\ldots,M$)—this is exactly the rough initial value that downstream windows receive in the PIDC pipeline.

![Source Figure 3.5: windowwise errors of IDC and PIDC under different regularities and viscosities](assets/papers/time-parallelization/source-figures/figure-3-5.svg)

The four panels of Figure 3.5 form two pairs of comparisons. The horizontal axis is the odd-numbered time windows ($1,3,\ldots,29$), and the vertical axis is the relative error (log scale); each panel simultaneously shows the initial error and the IDC/PIDC errors after the first and second sweeps; the legend of (a) applies to the remaining panels.

- Panels (a) and (b) use the sharp source term $g(x,t)$ with $\sigma=1000$ (see (2.4)), approximating a $\delta$-function-type source, so the solution has insufficient regularity in time; (c) and (d) use the smooth source term with $\sigma=5$, and the solution is sufficiently regular.
- Panels (a) and (c) take large viscosity $\nu=1$; (b) and (d) take small viscosity $\nu=10^{-3}$.

Panel-by-panel conclusions:

- **(a) low regularity, strong diffusion**: IDC and PIDC perform similarly, but after the first correction the error **no longer decreases**—the solution is not regular enough, so the high-order approximation has nothing to work with.
- **(b) low regularity, weak diffusion**: with reduced viscosity, even the improvement of the first IDC sweep is far worse than in (a), and further iteration helps little; PIDC behaves the same way.
- **(c) high regularity, strong diffusion**: the solution is smooth enough that **both IDC and PIDC continue to reduce the error on the second sweep**, and PIDC performs comparably to sequential IDC—this is the ideal case in which the order-raising mechanism works normally.
- **(d) high regularity, weak diffusion**: small viscosity weakens the order-raising effect; the second PIDC sweep clearly lags behind sequential IDC, because under weak diffusion the operator is nearly hyperbolic, the error propagates along characteristics rather than being smeared out by dissipation, and the contamination from the rough initial value is harder to repair with subsequent corrections.

The paper concludes on this basis: for hyperbolic problems, if the solution's regularity is insufficient, PIDC is not suitable for PinT (parallel-in-time) computation. The benefit of order-raising methods depends essentially on the smoothness of the solution.

## 3.3.2 Revisionist IDC (RIDC)

### A sliding quadrature window

RIDC was proposed by Christlieb, Macdonald, and Ong (2010). It refines the parallel granularity from the "time window" down to the "correction level", using a **sliding IDC quadrature interval** to achieve finer-grained parallelism. For this it takes an equally spaced-node quadrature rule.

Its schedule unfolds processor by processor:

1. The **first processor** advances continuously forward with a low-order time stepper (such as backward Euler), just like ordinary IDC; but after computing the first $M$ steps it **does not stop**, continuing to advance steps $M+1$, $M+2$, …, producing a steady stream of low-order trajectory values.
2. The **second processor** starts the first level of IDC correction once the first processor has accumulated the first $M$ values. After computing the first IDC interval (correcting step $M$) it also does not stop, but slides its **quadrature interval and quadrature nodes together to the right by one fine time step**: where it previously used steps $1,2,\ldots,M$ of the first processor, it now uses steps $2,3,\ldots,M+1$ as the IDC interval and quadrature nodes, and from these computes the correction of step $M+1$; then slides to steps $3,4,\ldots,M+2$ to compute step $M+2$, and so on.
3. The **third processor** starts once enough second-level data is available, using the same sliding rule to apply another level of correction to the output of the second processor. Higher-level processors follow analogously.

In the stable phase, each processor is responsible for a fixed correction level, and the levels advance adjacent time steps simultaneously, forming a pipeline that translates along the time axis. Unlike the window-blocked PIDC, every step of RIDC uses a "centered" sliding quadrature window, so it must store a number of historical values (window width $M$) and, like all pipelines, handle the startup (fill) and drain phases.

![Source Figure 3.6: windowwise errors of IDC and RIDC on the same advection–diffusion data](assets/papers/time-parallelization/source-figures/figure-3-6.svg)

Figure 3.6 reuses the PDE, grid, source, and viscosity settings of Figure 3.5, replacing only PIDC with RIDC, and takes initial condition $u(x,0)=0$. The correspondence of the four panels is unchanged: the top row is $\sigma=1000$ (low regularity), the bottom row is $\sigma=5$ (high regularity); the left column is $\nu=1$ (strong diffusion), the right column is $\nu=10^{-3}$ (weak diffusion). Hence (a)–(d) still represent "low regularity/strong diffusion", "low regularity/weak diffusion", "high regularity/strong diffusion", and "high regularity/weak diffusion" respectively. The bottom row allows the correction levels to continue reducing the error, while the top row and right column expose the difficulties brought by insufficient regularity and weak diffusion.

> [!tip] Insight
> RIDC changes only the **scheduling and memory structure** (sliding window, processors assigned by level) and does not change the mathematical kernel of IDC, so it inherits IDC's dependence on regularity: RIDC is still a high-order approximation technique, and for low-regularity hyperbolic solutions the error likewise cannot keep decreasing by adding correction levels. In other words, PIDC and RIDC solve "how to parallelize", not "how to raise the order under low regularity"—the latter is a physical limit shared by the whole IDC family.

## Equation and figure coverage check

| Source item        | Paper section    | Coverage status                                                                                                                              |
| ------------------ | ---------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| (3.5)–(3.8)        | 3.3 introduction | integral equation, residual, integral error equation, and differential error equation, including initial values and the smoothing motivation |
| (3.9)–(3.11)       | 3.3 introduction | $\theta$ discretization, residual increment, Lagrange quadrature weights, and the final IDC update                                           |
| Theorem 3.3        | 3.3 introduction | base order, $p$ orders gained per round, quadrature order ceiling $M$, SDC/Gauss–Lobatto/PFASST                                              |
| Window partition   | 3.3 introduction | windowing over long times, source of the sequential dependence                                                                               |
| Figure 3.4         | 3.3.1            | complete PIDC startup and steady-state pipeline diagram, bootstrap, $M=6/k_{\max}=4$                                                         |
| (3.12), Figure 3.5 | 3.3.1            | periodic difference matrices, all parameters, error definition, and four-panel conclusions                                                   |
| Figure 3.6         | 3.3.2            | RIDC sliding-window per-processor mechanism and regularity experiment                                                                        |

## Source of this page

- M. J. Gander, S.-L. Wu, and T. Zhou, [_Time Parallelization for Hyperbolic and Parabolic Problems_](https://doi.org/10.1017/S0962492924000072), _Acta Numerica_ 34 (2025), Section 3.3, pp. 405–412.

Supplementary references (cited in the paper, by author and year):

- A. Dutt, L. Greengard and V. Rokhlin (2000), _Spectral deferred correction methods for ordinary differential equations_ — the original IDC/SDC reference.
- D. Guibert and D. Tromeur-Dervout (2007), _Parallel deferred correction method for CFD problems_ — PIDC proposed.
- D. E. Womble (1990), _A time-stepping algorithm for parallel computers_, SIAM J. Sci. — source of the PIDC pipeline idea.
- A. J. Christlieb, C. B. Macdonald and B. W. Ong (2010), _Parallel high-order integrators_ — RIDC proposed.
- M. L. Minion (2010), _A hybrid parareal spectral deferred corrections method_; M. Emmett and M. L. Minion (2012), PFASST — the combination of SDC and Parareal.
- K. Böhmer and H. J. Stetter (1984) — an earlier source of the deferred/defect correction idea.
