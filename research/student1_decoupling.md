# Large Value Estimates and l^2 Decoupling: Finding Slack in Guth--Maynard

## Technical Report — Student 1

**Subject:** Analysis of the bottleneck in the Guth--Maynard zero density estimate and prospects for improvement via decoupling, restriction, and additive combinatorics.

**Date:** February 2026

---

## 1. Introduction and Statement of the Result

In May 2024, Larry Guth and James Maynard posted the paper "New large value estimates for Dirichlet polynomials" (arXiv:2405.20552), subsequently accepted for publication in the *Annals of Mathematics* (received July 2024, accepted April 2025). Their main consequence is the zero density estimate

$$
N(\sigma, T) \leq T^{\frac{30}{13}(1-\sigma) + o(1)},
$$

or equivalently, the uniform bound on the zero density exponent:

$$
\|A\|_\infty := \sup_{1/2 < \sigma < 1} A(\sigma) \leq \frac{30}{13} \approx 2.3077,
$$

where $A(\sigma)$ is the infimum of exponents $A$ such that $N(\sigma, T) \ll T^{A(1-\sigma)+o(1)}$. This represents the first improvement to the supremum norm of $A(\sigma)$ since Huxley's 1972 result giving $\|A\|_\infty \leq 12/5 = 2.4$, and by extension the first substantial improvement to the Ingham bound $N(3/4, T) \ll T^{3/5+o(1)}$ from 1940.

The purpose of this report is to:
1. Trace the origin of the exponent 30/13 through the chain of inequalities.
2. Identify the precise bottleneck step.
3. Evaluate whether newer decoupling results, the restriction conjecture, multilinear estimates, or additive combinatorics could yield further improvements.

---

## 2. Architecture of the Proof: Where Does 30/13 Come From?

### 2.1 The Large Value Framework

The standard approach to zero density estimates proceeds through *large value estimates* for Dirichlet polynomials. Define:

- $D_N(t) = \sum_{n \sim N} a_n n^{-it}$, where $|a_n| \leq 1$.
- $W \subset [0, T]$ is a 1-separated set of points where $|D_N(t)| \geq N^\sigma$.
- $C(\sigma, \alpha)$ is the large value exponent: the infimum of $C$ such that $|W| \leq T^{C+o(1)}$, where $\alpha = \log N / \log T$.

The connection to zero density is given by the fundamental inequality (Lemma 11.5 in the ANTEDB blueprint):

$$
A(\sigma)(1-\sigma) \leq \max\!\left(\sup_{\tau \geq 2} \frac{\mathrm{LV}_\zeta(\sigma, \tau)}{\tau},\; \limsup_{\tau \to \infty} \frac{\mathrm{LV}(\sigma, \tau)}{\tau}\right),
$$

where $\mathrm{LV}(\sigma, \tau)$ encodes the large value exponent in a normalized form, and $\mathrm{LV}_\zeta$ is the specialized version for partial sums of the zeta function.

### 2.2 Classical Inputs

Several classical large value estimates feed into the optimization:

**L^2 Mean Value Theorem (Equation 11 of Tao's blog):**
$$
C(\sigma, \alpha) \leq \min\!\big((2-2\sigma)\alpha,\; 1 + (1-2\sigma)\alpha\big).
$$

**Power Raising (Equation 12):** For natural $k \geq 1$,
$$
C(\sigma, \alpha) \leq C(\sigma, k\alpha).
$$

**Huxley's Extension (Equation 17):**
$$
C(\sigma, \alpha) \leq \min\!\big((2-2\sigma)\alpha,\; 1 + (4-6\sigma)\alpha\big),
$$
which improves on the L^2 bound for $\sigma > 3/4$, via a subdivision argument.

**Monotonicity (Equation 16):** The quantity $(1 - C(\sigma, \alpha))/\alpha$ is non-decreasing in $\alpha$.

### 2.3 The Halasz--Montgomery Barrier at sigma = 3/4

When $r = 2$ (i.e., using the second moment / Halasz--Montgomery method), the trace $\mathrm{tr}(M_W^* M_W)$ controls $s_1(M_W)^2$, and the resulting large value bound is essentially the L^2 mean value theorem. This yields no non-trivial information when $\sigma \leq 3/4$, because the bound from orthogonality alone gives $|W| \leq N^{2-2\sigma} \cdot T$, and for $\sigma = 3/4$ this is $|W| \leq N^{1/2} \cdot T \approx T^{1+\alpha/2}$, which after optimization over $\alpha$ reproduces the Ingham bound $A(3/4) \leq 3/(2 - 3/4) = 12/5$.

This is the fundamental barrier: the $r = 2$ method cannot see below $\sigma = 3/4$.

### 2.4 The Guth--Maynard Innovation: r = 3 (Sixth Power)

The key innovation is to use $r = 3$, i.e., bound the largest singular value $s_1(M_W)$ via the trace of $(M_W^* M_W)^3$:

$$
s_1(M_W) \leq \mathrm{tr}\!\big((M_W^* M_W)^3\big)^{1/6}.
$$

This involves six factors of the matrix and leads to what Tao described as "raising to the sixth (!) power." The resulting trace expansion decomposes into three sums:

- **$S_1$:** Negligible, bounded by $O_\epsilon(T^{-10})$ (Proposition 5.1).
- **$S_2$:** Controlled by Heath-Brown's theorem on double zeta sums (Proposition 6.1). This uses existing technology.
- **$S_3$:** The main term. Its estimation occupies the majority of the paper (roughly 34 of 48 pages).

**Proposition 8.1** ("$S_3$ controlled by energy") states:
$$
S_3 \ll T^2 \, |W|^{1/2} \, E(W)^{1/2},
$$
where $E(W) = |\{(w_1, w_2, w_3, w_4) \in W^4 : w_1 + w_2 = w_3 + w_4\}|$ is the additive energy of $W$.

### 2.5 The Parameter Optimization Yielding 30/13

The Guth--Maynard main technical estimate (Proposition 3.1) gives, in the notation of Tao's blog (Equation 20):

$$
C(\sigma, \alpha) \leq 1 + \left(\frac{12}{5} - 4\sigma\right)\alpha \qquad \text{for } 0.7 \leq \sigma \leq 0.8,\ \alpha = \frac{5}{6}.
$$

Feeding this into the zero density machinery (Equation 21) yields:

$$
A(\sigma) \leq \frac{15}{3 + 5\sigma} \qquad \text{for } 0.7 \leq \sigma \leq 0.8.
$$

The supremum of $A(\sigma)$ is attained by balancing the Guth--Maynard bound with the Ingham bound $A(\sigma) \leq 3/(2-\sigma)$:

$$
\frac{15}{3+5\sigma} = \frac{3}{2-\sigma} \implies 15(2-\sigma) = 3(3+5\sigma) \implies 30 - 15\sigma = 9 + 15\sigma \implies \sigma = \frac{21}{30} = \frac{7}{10}.
$$

At $\sigma = 7/10$:
$$
A(7/10) = \frac{15}{3 + 7/2} = \frac{15}{13/2} = \frac{30}{13}.
$$

**This is the origin of 30/13.** The supremum of $A(\sigma)$ has moved from $\sigma = 3/4$ (where it was stuck at $12/5$ for 80 years) to $\sigma = 7/10$, where the Guth--Maynard bound and the Ingham bound meet.

### 2.6 The Heuristic Principle (ANTEDB Section 11.3)

The ANTEDB blueprint codifies a useful heuristic (Corollary 11.9): if one can establish $\mathrm{LV}(\sigma, \tau_0) \leq 3 - 3\sigma$ for some parameter $\tau_0$ (a "Montgomery conjecture with 1.5x loss"), then one expects $A(\sigma) \leq 3/\tau_0$. For Guth--Maynard, the choice $\tau_0 = (3+5\sigma)/5$ at $\sigma = 7/10$ gives $\tau_0 = 13/10$, and indeed $3/\tau_0 = 30/13$.

---

## 3. Identifying the Precise Bottleneck

Having traced the exponent, we can now identify three distinct bottleneck layers:

### Bottleneck 1 (Primary): The limitation to r = 3

The authors state explicitly that they "do not know how to obtain good bounds when $r \geq 4$." If one could take $r = 4$ with comparable control on the resulting sums, the trace bound $s_1(M_W) \leq \mathrm{tr}((M_W^* M_W)^4)^{1/8}$ would provide strictly stronger information, potentially pushing the large value estimate below the current threshold. Indeed, the full Montgomery large value conjecture ($\mathrm{LV}(\sigma, \tau) \leq 2 - 2\sigma$ for all relevant $\tau$) would follow from sharp bounds on $\mathrm{tr}((M_W^* M_W)^r)$ for all $r$, and this would imply the density hypothesis $\|A\|_\infty \leq 2$.

**The inability to handle $r \geq 4$ is the single most important bottleneck.** The combinatorial complexity of the trace expansion grows rapidly with $r$: for $r = 3$, one must control sums over sextuples $(n_1, n_2, n_3, m_1, m_2, m_3)$ with constraints involving $\log(n_i/m_i)$. For $r = 4$, one would need to handle octuples, and the Poisson summation / stationary phase analysis becomes dramatically more intricate.

### Bottleneck 2 (Secondary): The $S_3$ bound via additive energy

The bound $S_3 \ll T^2 |W|^{1/2} E(W)^{1/2}$ (Proposition 8.1) is the heart of the new contribution. This estimate involves:

(a) A Poisson summation step exploiting the arithmetic structure of $\{\log n\}$.

(b) A "refusal to simplify" via stationary phase — Guth and Maynard keep a complicated Fourier integral in its natural form rather than applying standard asymptotic expansions, which would lose the crucial cancellation.

(c) Multilinear estimates in the spirit of Keil--Zhao, controlling the resulting oscillatory expressions.

The quality of the $S_3$ bound directly determines the strength of the large value estimate. If one could improve the exponent $1/2$ on $E(W)$ to something smaller (say by proving that the relevant oscillatory integrals exhibit additional cancellation), this would immediately improve the zero density estimate.

### Bottleneck 3 (Tertiary): The intersection point at sigma = 7/10

The exponent 30/13 arises as the maximum of two curves: the Ingham curve $3/(2-\sigma)$ (controlling the regime $\sigma < 7/10$) and the Guth--Maynard curve $15/(3+5\sigma)$ (controlling $7/10 < \sigma < 0.76$). Any improvement to the supremum must either:

- Improve the Guth--Maynard estimate in a neighborhood of $\sigma = 7/10$, or
- Improve the Ingham bound for $\sigma$ slightly below $7/10$, or
- Push the intersection point to a lower value of $A$.

Note that the Ingham bound $A(\sigma) \leq 3/(2-\sigma)$ is extremely robust — it follows from a one-page orthogonality argument that has resisted improvement for over a century except via the deep new ideas of Guth--Maynard. This suggests that Bottleneck 3 is essentially a consequence of Bottleneck 1.

**Verdict:** The primary bottleneck is the limitation to $r = 3$ in the singular value / trace method. This is not a decoupling limitation per se; it is a limitation in our ability to control higher-order multiplicative structure in the integers.

---

## 4. Can Newer Decoupling Results Help?

### 4.1 The Role of Decoupling in Guth--Maynard

It is important to clarify a common misconception: the Guth--Maynard paper does **not** directly invoke the Bourgain--Demeter l^2 decoupling theorem (2015) as a black box. Rather, the *ideas* from decoupling theory — wave packet decomposition, induction on scales, and the philosophy of exploiting curvature — inform the proof strategy, particularly in the analysis of $S_3$. The paper is largely Fourier-analytic in nature, and the key curvature being exploited is that of the frequency set $\Phi = \{\log n : n \sim N\}$, which has specific arithmetic properties quite different from the paraboloid.

The Bourgain--Demeter--Guth resolution of the Vinogradov mean value theorem (2015) — which *is* a direct application of l^2 decoupling — feeds into the zero density framework through the improved bounds on Weyl sums and the resulting large value estimates encoded in the L^2 mean value theorem. These are already incorporated into the Guth--Maynard optimization.

### 4.2 Recent Decoupling Advances (2024--2026)

Several notable refinements have appeared:

1. **Guth--Maldague--Oh (2024/2025):** Sharp $(l^2, L^p)$ decoupling for general surfaces in $\mathbb{R}^3$ with non-vanishing Gaussian curvature, for the full range $2 \leq p \leq 4$ (arXiv:2403.18431). This introduces a new box decomposition of the delta-neighborhood. However, the relevance to the Guth--Maynard setting is indirect, since the frequency set $\{\log n\}$ is one-dimensional, not a surface.

2. **Guth--Maldague--Wang (JEMS 2024):** Improved decoupling for the parabola, achieving the $(l^2, l^6)$ estimate with a logarithmic constant $(\log R)^c$. This is an improvement in the *constant* (or sub-polynomial factor) rather than in the *exponent*, and so would affect only the $o(1)$ terms in the zero density estimate — not the main exponent 30/13.

3. **Guth--Maldague (Analysis & PDE 2024):** Sharp small cap decoupling for the moment curve in $\mathbb{R}^3$, resolving a conjecture of Demeter--Guth--Wang. This is relevant to Vinogradov-type mean value theorems for cubic and higher-degree polynomials, and could potentially improve the large value inputs coming from the mean value theorem. However, the critical regime for Guth--Maynard is $\sigma \approx 3/4$, where the mean value theorem is already superseded by their direct approach.

### 4.3 Assessment

**Expected gain from newer decoupling: negligible for the main exponent.** The bottleneck is not the l^2 decoupling constant or the mean value theorem input; it is the limitation to the sixth moment ($r = 3$) in the singular value method. Improvements to decoupling constants would affect sub-polynomial factors only. Qualitatively new decoupling results (e.g., for the specific frequency set $\{\log n\}$) could potentially help with the $S_3$ bound, but no such results are currently known or conjectured.

---

## 5. The Restriction Conjecture Connection

### 5.1 What Would Full Restriction Give?

The restriction conjecture for the paraboloid in $\mathbb{R}^d$ asserts certain $L^p$ bounds for the Fourier extension operator. For $d = 2$ (the relevant case for one-dimensional exponential sums), the conjecture is already proved (Fefferman 1970, refined by Zygmund). The open cases in higher dimensions do not directly bear on the Guth--Maynard framework, which deals with one-dimensional sums $\sum a_n n^{-it}$.

However, there is an indirect connection via the **exponent pair conjecture**. The restriction conjecture for the paraboloid implies improved bounds on oscillatory integrals via Bourgain--Guth type arguments, which feed into improved exponent pairs $(k, \ell)$, which then feed into improved zero density estimates via the Ivic--Bourgain machinery.

### 5.2 Assuming the Full Restriction Conjecture

If one assumes the full restriction conjecture for the paraboloid in all dimensions, the ANTEDB machinery can optimize the resulting exponent pairs and zero density bounds. The expected effect is:

- Improvement to zero density estimates for $\sigma$ near 1 (say $\sigma > 0.85$), where the Ivic--Bourgain method already gives bounds close to the density hypothesis.
- Little or no improvement near the critical intersection point $\sigma = 7/10$, because that region is controlled by the Guth--Maynard estimate, not by exponent pair methods.

**Quantitative estimate:** Assuming full restriction, the Tao--Trudgian--Yang ANTEDB optimization would likely verify the density hypothesis $A(\sigma) \leq 2$ for $\sigma \geq \sigma_0$ for some $\sigma_0$ slightly below the current threshold of $25/32 \approx 0.781$ (Bourgain's result). But the supremum $\|A\|_\infty$ would remain at 30/13, governed by the $\sigma = 7/10$ intersection.

### 5.3 What Would Push the Density Estimate Further?

To improve $\|A\|_\infty$ below 30/13, one needs progress on the large value problem itself — specifically, a new bound on $\mathrm{LV}(\sigma, \tau)$ for $\sigma$ near $7/10$ that beats the current Guth--Maynard estimate. This would require either:

(a) Extending the $r = 3$ method with a sharper $S_3$ bound (a plausible but technically formidable target), or

(b) Making $r = 4$ work (a fundamental advance that currently seems out of reach), or

(c) An entirely new approach to the large value problem.

If one could prove Montgomery's full large value conjecture, one would obtain $\|A\|_\infty \leq 2$ (the density hypothesis). The Lindelof hypothesis would give $A(\sigma) \leq 0$ for $\sigma > 3/4$. The Riemann Hypothesis gives $A(\sigma) = -\infty$ for $\sigma > 1/2$.

---

## 6. Multilinear vs. Linear Estimates

### 6.1 The Bennett--Carbery--Tao Framework

Bennett, Carbery, and Tao (2006) proved $d$-linear restriction estimates that are stronger than the best known linear estimates in many parameter ranges. Guth proved the endpoint case of the multilinear Kakeya conjecture (2010) using the polynomial method. Bourgain and Guth (2011) applied these multilinear results to obtain improved linear restriction estimates via a "multilinear-to-linear" reduction.

The key question is: can the multilinear-to-linear reduction be sharpened in the specific setup that Guth--Maynard use?

### 6.2 Application to Guth--Maynard

The Guth--Maynard proof does not directly use multilinear restriction estimates in the Bennett--Carbery--Tao sense. The "multilinear" structure that appears is in the trace expansion $(M_W^* M_W)^3$, which involves products of six matrix entries. This is more naturally a question about multiplicative structure of integers than about restriction to curved surfaces.

However, the Bourgain--Guth philosophy — of decomposing into "broad" (multilinear) and "narrow" (concentrated) contributions — does inform certain steps. In particular, the distinction between the $S_2$ and $S_3$ terms can be viewed as a decomposition into contributions where the frequency variables are "transverse" versus "aligned."

### 6.3 Potential for Improvement

Tao's 2019 paper "Sharp bounds for multilinear curved Kakeya, restriction and oscillatory integral estimates away from the endpoint" removed epsilon losses in the multilinear estimates via a "ball inflation" induction-on-scales scheme inspired by decoupling proofs. If similar epsilon-removal techniques could be applied within the Guth--Maynard $S_3$ analysis, one might improve the sub-polynomial factors. But again, this would not affect the main exponent.

**Assessment:** The multilinear-to-linear reduction is not the binding constraint. The potential for improvement here is marginal — at best, refinements to sub-polynomial losses.

---

## 7. Additive Combinatorics Inputs

### 7.1 The Role of Additive Energy in Guth--Maynard

The additive energy $E(W)$ of the set $W$ plays a central role. The proof employs a dichotomy:

- **Low additive energy:** If $E(W) \ll |W|^{7/3-\delta}$ for some $\delta > 0$, then the $S_3$ bound (Proposition 8.1) directly gives an improvement over orthogonality.
- **High additive energy:** If $E(W) \approx |W|^{7/3}$, then $W$ has significant additive structure, and Heath-Brown's theorem (on double zeta sums for structured sets) provides the needed bound on $S_2$.

This is a "sum-product type" dichotomy: a set cannot simultaneously have high additive structure (large $E(W)$) and be "generic" with respect to multiplicative structure (which is what the oscillatory integrals in $S_3$ measure).

### 7.2 Kelley--Meka and Recent Advances

The Kelley--Meka breakthrough (2023) on 3-term arithmetic progressions establishes that any subset $A \subseteq \{1, \ldots, N\}$ with $|A| \geq \exp(-c(\log N)^{1/12}) \cdot N$ contains a non-trivial 3-AP. This uses Fourier-analytic techniques (density increment, Bohr sets, $L^p$ norm control) that have structural parallels with the additive energy arguments in Guth--Maynard.

However, the direct applicability is limited:

1. **Different setting:** Kelley--Meka work with additive structure in $\mathbb{Z}/N\mathbb{Z}$, while Guth--Maynard work with the set $W$ of large value points on the real line, with additive energy measured in $\mathbb{R}$.

2. **Different quantity:** The relevant quantity for Guth--Maynard is the additive energy $E(W)$, not the density of arithmetic progressions. While both involve four-variable additive relations ($w_1 + w_2 = w_3 + w_4$ vs. $a_1 - 2a_2 + a_3 = 0$), the structural implications are different.

3. **Sum-product estimates:** More directly relevant are recent improvements to sum-product estimates in $\mathbb{R}$, which quantify the extent to which a set can simultaneously have additive and multiplicative structure. Any improvement to the Erdos--Szemeredi conjecture (particularly the additive energy of multiplicatively structured sets) could in principle sharpen the dichotomy in Guth--Maynard.

### 7.3 Concrete Possibility

The most promising additive combinatorics input would be an improved *inverse theorem for additive energy*: a sharper structural characterization of sets $W \subset [0, T]$ with $E(W) \geq |W|^{7/3-\delta}$. If one could show that such sets must have very rigid arithmetic structure (e.g., concentration in generalized arithmetic progressions with specific parameters), this could allow the Heath-Brown type argument to be applied more efficiently, reducing the "overlap region" in the dichotomy.

**Assessment:** Plausible but speculative. A gain of $\epsilon$ in the additive energy threshold (from $|W|^{7/3}$ to $|W|^{7/3-\epsilon}$) would translate to a gain of roughly $O(\epsilon)$ in the zero density exponent. Major breakthroughs in additive combinatorics (like a full resolution of the Erdos--Szemeredi conjecture) would be needed for a substantial improvement.

---

## 8. The ANTEDB Perspective: Systematic Optimization

The Tao--Trudgian--Yang ANTEDB project (arXiv:2501.16779, January 2025) provides a systematic framework for optimizing zero density estimates from all available inputs. Key features:

### 8.1 Current State of the Art (ANTEDB Table 11.1)

| $\sigma$ range | Best bound on $A(\sigma)$ | Source |
|---|---|---|
| $1/2 \leq \sigma \leq 0.7$ | $3/(2-\sigma)$ | Ingham (1940) |
| $0.7 \leq \sigma < 0.76$ | $15/(3+5\sigma)$ | Guth--Maynard (2024) |
| $0.76 \leq \sigma < 0.89$ | Various piecewise | Ivic/Bourgain/Heath-Brown |
| $0.89 \leq \sigma < 0.9$ | $24/(30\sigma - 11)$ | Chen--Debruyne--Vidas |
| $\sigma \geq 0.9$ | Multiple piecewise | Pintz/Bourgain optimized |

The supremum $\|A\|_\infty = 30/13$ is attained at $\sigma = 7/10$.

### 8.2 What the ANTEDB Reveals About Slack

The ANTEDB optimization is *exhaustive* within the framework of known techniques. If there were unexploited slack in the chain of known inequalities — for instance, if combining two known results in a previously overlooked way could improve the bound — the computer search would have found it.

The fact that the ANTEDB reproduces $30/13$ as the optimal value from all known inputs confirms that **new mathematical input is required** for further improvement. The database does, however, identify specific targets:

- A new large value estimate at $(\sigma, \tau) = (0.7, 1.3)$ that beats the current bound of $\mathrm{LV}(0.7, 1.3) \leq 0.9$ would lower the supremum.
- Any improvement to the Ingham bound $A(\sigma) \leq 3/(2-\sigma)$ for $\sigma < 0.7$ would also lower the supremum (but this seems extremely hard).

### 8.3 Guth's 2025 Survey

In March 2025, Guth posted a 49-page survey (arXiv:2503.07410) titled "Large value estimates in number theory, harmonic analysis, and computer science," connecting the large value problem for Dirichlet polynomials with the restriction problem in harmonic analysis and problems from theoretical computer science. This survey makes explicit the observation that for a "large range of parameters where the operator norm method (based on orthogonality, taking less than a page, known for a hundred years) is still the best known," and that "it is really hard to improve on this one page argument." This underscores the depth of the bottleneck.

---

## 9. Concrete Proposals for Improvement

Based on the analysis above, I rank the following proposals by feasibility and expected gain:

### Proposal A: Sharpen the S_3 bound (Medium feasibility, moderate gain)

**Idea:** Improve Proposition 8.1 from $S_3 \ll T^2 |W|^{1/2} E(W)^{1/2}$ to $S_3 \ll T^2 |W|^{1/2} E(W)^{1/2-\delta}$ for some $\delta > 0$, by exploiting additional cancellation in the oscillatory integrals. This could come from a more careful Poisson summation argument or from new bounds on exponential sums over products of logarithms.

**Expected gain:** A gain of $\delta$ in the energy exponent would translate to a gain of roughly $\delta \cdot (1-\sigma)$ in the zero density exponent near $\sigma = 7/10$, potentially lowering $\|A\|_\infty$ by a constant times $\delta$.

**Feasibility:** This requires new ideas for the specific oscillatory integrals that arise in the $S_3$ analysis. The "refusal to use stationary phase" suggests that the current bound may already be close to optimal for the method, but this is not proven. A careful case analysis might reveal additional structure.

### Proposal B: Improve the additive energy dichotomy (Medium feasibility, small gain)

**Idea:** Use improved inverse theorems for additive energy (from recent work in additive combinatorics) to sharpen the threshold at which the Heath-Brown argument kicks in. Specifically, if one can characterize sets with $E(W) \geq |W|^{7/3-\delta}$ more precisely, the dichotomy can be tightened.

**Expected gain:** Small, perhaps lowering $30/13$ by $O(0.01)$.

**Feasibility:** Depends on progress in additive combinatorics. The Erdos--Szemeredi conjecture and its relatives are the relevant open problems.

### Proposal C: Extend to r = 4 (Low feasibility, large gain)

**Idea:** Develop techniques to bound $\mathrm{tr}((M_W^* M_W)^4)$ with sufficient precision. This would require controlling sums over octuples of integers with multiplicative constraints, a dramatic increase in complexity.

**Expected gain:** If achievable with comparable quality to the $r = 3$ case, this could lower $\|A\|_\infty$ significantly — potentially to around 2.1 or below, approaching the density hypothesis.

**Feasibility:** Very low with current techniques. The combinatorial explosion in the number of terms and the difficulty of the resulting oscillatory integral analysis make this a long-term goal rather than a near-term target. This would likely require a fundamentally new idea for organizing the trace expansion.

### Proposal D: New decoupling for log-frequency sets (Low feasibility, uncertain gain)

**Idea:** Develop an l^2 decoupling theorem specifically adapted to the frequency set $\Phi = \{\log n : n \sim N\}$. Unlike the paraboloid (which has non-vanishing curvature everywhere), $\Phi$ is a discrete set with specific arithmetic structure. A "decoupling" result for this set would directly bound the operator norm of $M_W$.

**Expected gain:** Uncertain. If such a result gave square-root cancellation in the appropriate sense, it could be transformative.

**Feasibility:** Very low. The arithmetic structure of $\{\log n\}$ is fundamentally different from the geometric structure (curvature) that drives standard decoupling results. No framework currently exists for this type of result.

### Proposal E: Computer-assisted optimization via ANTEDB (High feasibility, zero gain on main exponent)

**Idea:** Continue developing the ANTEDB to find previously overlooked combinations of known results.

**Expected gain:** The ANTEDB has already been run and confirms $30/13$ as optimal from known inputs. Further gains require new mathematical inputs.

**Feasibility:** High (the infrastructure exists), but the output has already been maximized.

---

## 10. Honest Assessment and Conclusions

### What is achievable in the near term (1--3 years)?

- **Refinements to sub-polynomial factors:** Improved decoupling constants (Guth--Maldague--Wang type) and epsilon-removal techniques (Tao 2019 type) could improve the $o(1)$ term in the zero density estimate. This is publishable but does not change the main exponent.
- **Minor improvements via Proposal B:** If progress is made on additive energy inverse theorems, a small improvement to $\|A\|_\infty$ (perhaps to $2.30$ or so) is conceivable.
- **ANTEDB integration of new inputs:** Any new large value estimate, exponent pair, or exponential sum bound can be immediately fed into the ANTEDB to check for downstream improvements.

### What would require a breakthrough?

- **Lowering $\|A\|_\infty$ below 2.2** likely requires either extending to $r = 4$ or finding a qualitatively new approach to the large value problem.
- **The density hypothesis** ($\|A\|_\infty \leq 2$) remains extremely distant. It would follow from the full Montgomery large value conjecture, which in turn would require sharp trace bounds for all $r$ — essentially a complete understanding of the multiplicative structure of the integers as seen through exponential sums.

### The fundamental tension

The Guth--Maynard result reveals a deep tension in the problem: the frequency set $\{\log n\}$ has both *arithmetic* structure (which makes Poisson summation effective) and *analytic* structure (which makes Fourier methods effective), but neither type of structure alone is sufficient to reach the density hypothesis. The l^2 decoupling theory is optimized for *curvature* — the non-vanishing of the second derivative of the phase function — while the arithmetic of logarithms is more subtle, involving the distribution of prime factors and the rigidity of multiplicative number theory.

The exponent $30/13$ is, in a precise sense, the best that can be achieved by combining the Ingham orthogonality argument (which is sharp for $\sigma < 7/10$) with a single application of the $r = 3$ trace method (which is sharp in the range $7/10 < \sigma < 0.76$). To go further, one must either improve the trace method itself or find an entirely different approach to the $\sigma < 7/10$ regime.

---

## References

1. L. Guth and J. Maynard, "New large value estimates for Dirichlet polynomials," [arXiv:2405.20552](https://arxiv.org/abs/2405.20552) (2024), accepted in *Annals of Mathematics* (2025).

2. T. Tao, "A computation-outsourced discussion of zero density theorems for the Riemann zeta function," [What's New blog](https://terrytao.wordpress.com/2024/07/07/a-computation-outsourced-discussion-of-zero-density-theorems-for-the-riemann-zeta-function/) (July 2024).

3. T. Tao, T. Trudgian, and A. Yang, "New exponent pairs, zero density estimates, and zero additive energy estimates: a systematic approach," [arXiv:2501.16779](https://arxiv.org/abs/2501.16779) (January 2025).

4. ANTEDB: [Analytic Number Theory Exponent Database](https://teorth.github.io/expdb/), GitHub: [teorth/expdb](https://github.com/teorth/expdb).

5. ANTEDB Blueprint, Zero Density Chapter: [teorth.github.io/expdb/blueprint/zero-density-chapter.html](https://teorth.github.io/expdb/blueprint/zero-density-chapter.html).

6. L. Guth, "Large value estimates in number theory, harmonic analysis, and computer science," [arXiv:2503.07410](https://arxiv.org/abs/2503.07410) (March 2025).

7. J. Bourgain and C. Demeter, "The proof of the l^2 decoupling conjecture," *Annals of Mathematics* 182(1), 351--389, 2015. [arXiv:1403.5335](https://arxiv.org/abs/1403.5335).

8. J. Bourgain, C. Demeter, and L. Guth, "Proof of the main conjecture in Vinogradov's mean value theorem for degrees higher than three," *Annals of Mathematics* 184(2), 633--682, 2016.

9. J. Bennett, A. Carbery, and T. Tao, "On the multilinear restriction and Kakeya conjectures," *Acta Mathematica* 196(2), 261--302, 2006.

10. L. Guth, "The endpoint case of the Bennett-Carbery-Tao multilinear Kakeya conjecture," *Acta Mathematica* 205, 263--286, 2010.

11. T. Tao, "Sharp bounds for multilinear curved Kakeya, restriction and oscillatory integral estimates away from the endpoint," [arXiv:1907.12348](https://arxiv.org/abs/1907.12348) (2019).

12. L. Guth, D. Maldague, and H. Wang, "Improved decoupling for the parabola," *JEMS* 26, 2024. [EMS Press](https://ems.press/journals/jems/articles/8736488).

13. L. Guth, D. Maldague, and C. Oh, "l^2 decoupling theorem for surfaces in R^3," [arXiv:2403.18431](https://arxiv.org/abs/2403.18431) (2024, revised 2025).

14. L. Guth and D. Maldague, "Small cap decoupling for the moment curve in R^3," *Analysis & PDE* 17(10), 2024. [MSP](https://msp.org/apde/2024/17-10/apde-v17-n10-p04-p.pdf).

15. J. Maynard and K. Pratt, "Half-isolated zeros and zero-density estimates," *IMRN* 2024(19), 12978--13014, 2024. [arXiv:2206.11729](https://arxiv.org/abs/2206.11729).

16. Z. Kelley and R. Meka, "Strong bounds for 3-progressions," [arXiv:2302.05537](https://arxiv.org/abs/2302.05537) (2023).

17. T. Bloom and O. Sisask, "The Kelley-Meka bounds for sets free of three-term arithmetic progressions," *Essential Number Theory* 2(1), 2023. [MSP](https://msp.org/ent/2023/2-1/ent-v2-n1-p02-p.pdf).

18. T. Tao, Mastodon announcement on Guth--Maynard, [mathstodon.xyz/@tao/112557248794707738](https://mathstodon.xyz/@tao/112557248794707738) (June 2024).
