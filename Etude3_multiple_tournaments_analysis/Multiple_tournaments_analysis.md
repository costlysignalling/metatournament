# Known talents in multiple identical tournaments

This working analysis measures region sizes by the fraction of actors under the model's independent lognormal talent distribution, with median 1 and original-scale standard deviation 0.5. All numerical results below use deterministic quadrature of analytically derived regions, not simulated actors. Exactly equal investment and positive investment in every tournament are different events: exact equality has probability zero for independent continuous talents.

## Allocation and normalization

For identical tournaments with a symmetric population environment, the actor maximizes

$$F_T(B;\boldsymbol\alpha)=\max_{x_t\ge0,\,\sum_t x_t=B}\sum_{t=1}^T w(x_t,\alpha_t),\qquad w(x,a)=[1+a g(x)]^\varphi-1,\quad g(x)=(1+x)^\beta-1.$$

The common multiplier $B/(T M_{\varphi,T})$ changes utility levels but not the allocation. Write $K$ for the number of positive allocations, and $P_k^{(T)}=\Pr(K=k)$.

## The concave case gives exact global region boundaries

For $0<\varphi\le1$ and $0<\beta<1$, each summand is strictly concave in investment. There is a unique optimum. Positive allocations equalize marginal returns at a common value $\lambda$; an unused tournament satisfies $\varphi\beta\alpha_t\le\lambda$.

Complete specialization in a tournament with talent $a$ is optimal if every other talent is at most

$$d_B(a)=a(1+B)^{\beta-1}[1+a g(B)]^{\varphi-1}.$$

For talent density $f$ and cumulative distribution $F$, the exact specialist probability is therefore the one-dimensional integral

$$P_1^{(T)}=T\int_0^\infty f(a)F(d_B(a))^{T-1}\,da.$$

For three tournaments, condition on the largest two talents $a\ge b$. If the optimal two-tournament allocation invests in both, let $\lambda_2(a,b;B)$ be their common marginal return and set $c_2=\lambda_2/(\varphi\beta)$. The third tournament enters precisely when its talent exceeds $c_2$. Since its talent cannot exceed $b$ under this ordering,

$$P_3^{(3)}=6\int_0^\infty\int_b^\infty f(a)f(b)[F(b)-F(c_2(a,b;B))]_+\,da\,db.$$

If the two-tournament optimum is a specialist, set the bracket to zero: the third-ranked talent cannot enter while the second-ranked talent is unused. Finally, $P_2^{(3)}=1-P_1^{(3)}-P_3^{(3)}$. The R implementation evaluates these exact regions using scalar root-finding and adaptive integration.

## Linear spoils returns reduce all three probabilities to one-dimensional integrals

At $\varphi=1$, put $p=1/(1-\beta)$ and $z_t=\alpha_t^p$. For an active set $A$ of size $k$,

$$x_t^*=(B+k)\frac{z_t}{\sum_{j\in A}z_j}-1\quad(t\in A),\qquad x_t^*=0\quad(t\notin A).$$

The active set consists of the largest talents and is identified by checking positive allocations for its members and non-entry for the rest. Equivalently, $x_t^*=[(\alpha_t/c)^p-1]_+$ with $c$ chosen to exhaust the budget.

Let $s$ denote the standard deviation of log-talent, let $\Phi$ and $\phi_N$ be the standard-normal cumulative distribution and density, and define $L=(1-\beta)\log(1+B)$. The median talent cancels from the strategy probabilities. For two tournaments,

$$P_2^{(2)}=1-2\Phi\left(-\frac{L}{\sqrt2s}\right).$$

For three tournaments,

$$P_1^{(3)}=3\int_{-\infty}^{\infty}\phi_N(z)\Phi(z-L/s)^2\,dz.$$

To calculate full participation, condition on tournament 3 having the minimum talent. Then $D_i=\log(\alpha_i/\alpha_3)$ have variance $2s^2$ and correlation $1/2$. Full participation means $D_1,D_2\ge0$ and $e^{pD_1}+e^{pD_2}<B+2$. Thus, with $h(d)=(1-\beta)\log[B+2-e^{d/(1-\beta)}]$,

$$P_3^{(3)}=3\int_0^L\frac{\phi_N(d/(\sqrt2s))}{\sqrt2s}\left[\Phi\left(\frac{h(d)-d/2}{s\sqrt{3/2}}\right)-\Phi\left(\frac{-d/2}{s\sqrt{3/2}}\right)\right],dd.$$

These are analytical probability formulas evaluated numerically. They do not require simulated talents or a discretized allocation simplex. Two independent integral reductions agree within $10^{-8}$ at the checked parameter values; the general concave-region integral also agrees with this special case.

For $B=50$ and $\varphi=1$, two-tournament participation becomes more frequent than three-tournament participation at $\beta\approx0.83855$. This is a comparison between $P_2^{(3)}$ and $P_3^{(3)}$, not a claim that two tournaments are the majority strategy everywhere above that threshold; specialists eventually dominate.

## Does full participation become more common when a tournament is added?

With fixed $B$ and $0<\varphi\le1$, no. Couple the two- and three-tournament problems by retaining the first two talents and adding an independent third talent. At a fixed marginal-return level, the additional option weakly increases total demanded investment. The marginal-return level required to exhaust the unchanged budget must therefore weakly increase. Investment in every pre-existing option weakly decreases. An excluded option cannot become active just because another option was added.

Consequently, full participation in three implies full participation in the corresponding two-tournament problem, so $P_3^{(3)}\le P_2^{(2)}$. The same reasoning applies from $T$ to $T+1$. This proof concerns full participation, not the probability of a particular intermediate support size.

In the same concave regime, increasing $B$ weakly increases the active set for each fixed talent vector. Also, increasing $\varphi$ within $(0,1]$ weakly decreases the active set: for any fixed $x>0$, the marginal return divided by $\varphi\beta$ increases with $\varphi$, while its value at zero stays equal to talent. The normalized marginal-return cutoff must rise to maintain the fixed budget. At $\varphi=1$, increasing $\beta$ likewise weakly decreases the active set. The last statement should not be extended to arbitrary $\varphi$ without a separate proof.

## Proportional budget growth

Using $B_T=bT$ keeps available investment per tournament constant. It does not change the parameter $\beta$; it controls the investment scale at which the response function is used. With $B_2=50$, the proportional-budget comparison uses $B_3=75$.

For fixed $T$ in the concave regime, the larger budget increases full participation and reduces complete specialization. Intermediate support frequencies can rise or fall. Comparing different $T$ while also increasing $B$ combines two opposing changes, so the fixed-budget nesting proof no longer applies.

For example, with $\beta=0.8,\varphi=1$, full participation has probability 0.800034 in two tournaments at $B=50$, 0.557968 in three tournaments at $B=50$, and 0.634084 in three tournaments at $B=75$. A universal claim that full participation must fall even under proportional budgets would be false: at $\beta=0.2$, the corresponding two-tournament probability is 0.999999705, while the proportional-budget three-tournament probability is 0.999999812. The difference is small but resolved by the deterministic integrals.

For arbitrary $T$ at $\varphi=1$, the finite-dimensional formula above still solves allocations exactly. Under $B_T=bT$, the empirical cutoff converges to the positive solution of

$$\mathbb E\left[\left((\alpha/c)^p-1\right)_+\right]=b.$$

This expectation has a closed expression in normal cumulative distributions. If $\log\alpha\sim N(m,s^2)$, it equals

$$c^{-p}e^{pm+p^2s^2/2}\Phi\left(\frac{m+ps^2-\log c}{s}\right)-\Phi\left(\frac{m-\log c}{s}\right).$$

The limiting fraction of entered tournaments is $1-F(c)$, which is below 1. Thus proportional budget growth need not preserve investment in every available tournament, even though it can preserve a positive fraction of them.

## Arbitrary positive phi: global rather than local entry tests

The simple marginal-return boundaries need not identify global optima when $\varphi>1$. There is nevertheless an exact reduction using globally optimized two-tournament values. Define

$$V_2(r;a,b)=\max_{0\le x\le r}\{w(x,a)+w(r-x,b)\}.$$

Higher talent has a strictly larger marginal return at any given positive investment, because $w_{x\alpha}>0$ for all $\varphi>0$. Swapping allocations therefore shows that an optimum can assign investments in decreasing talent order. Its active set contains the largest $k$ talents.

For fixed largest talents $a\ge b$, define the third-tournament entry threshold

$$C_3(a,b;B)=\inf_{0<z\le B}\frac{\{1+V_2(B;a,b)-V_2(B-z;a,b)\}^{1/\varphi}-1}{g(z)}.$$

This follows by comparing $w(z,c)+V_2(B-z;a,b)$ with $V_2(B;a,b)$ and solving the inequality for the third talent $c$. The third tournament improves global utility precisely when $c>C_3$, apart from boundary ties. Thus the preceding double integral for $P_3^{(3)}$ remains valid with $c_2$ replaced by $C_3$. For $\varphi\le1$, the infimum is the limit at $z=0$, namely the marginal-entry threshold. For larger $\varphi$, entry can instead involve a finite jump in investment.

There is an exact global specialist test too. For competitor talents $b,c$, define $a_r$ as the unique solution of

$$w(B,a_r)-w(B-r,a_r)=V_2(r;b,c),\qquad 0<r\le B,$$

and $A_1(b,c;B)=\sup_{0<r\le B}a_r$. The left side is strictly increasing in $a_r$, so the inverse is well-defined. Complete specialization is globally optimal precisely when the focal talent is at least $A_1$. Consequently,

$$P_1^{(3)}=3\int_0^\infty\int_0^\infty f(b)f(c)\{1-F[A_1(b,c;B)]\}\,db\,dc,$$

and $P_2^{(3)}=1-P_1^{(3)}-P_3^{(3)}$. These equations characterize the full dependence on $\beta,\varphi,B$. They are not a closed-form evaluation: outside the concave regime, the nested maxima and infima require global one-dimensional optimization before deterministic integration. Those nonconcave probability integrals have not yet been numerically evaluated in the accompanying scripts.

## If the third tournament is excluded

Yes: at fixed realized talents and the same $B$, excluding tournament 3 leaves exactly the two-tournament allocation objective already analyzed. A possible common population-normalization factor does not change its maximizing allocation.

Two qualifications matter. First, exclusion must concern the global optimum, not only a locally unprofitable small entry. Second, the two retained talents are selected from three, so their conditional distribution differs from two fresh independent draws. Actor-level allocation formulas transfer exactly; population frequencies and absolute normalized utility levels need not.

## Reproduction and scope

Run `analytic_participation.R` for the exact linear-spoils calculation and its independent checks. Run `concave_participation.R` for $0<\varphi\le1$ region integration. Both use base R only. Their CSV files give unrounded probabilities; displayed zeros in logs can reflect rounding and should not be interpreted as empty regions. The quadrature errors control numerical integration, while the region formulas themselves are analytical. No Paper or Supplement edits have been made by this research step.
