# Mean Value Theorems for Dirichlet Polynomials: Sharpening the Halász–Montgomery Framework

## Technical Research Report

---

## 1. The Halász–Montgomery Zero-Detection Method

### 1.1 The Zero-Detecting Polynomial

The fundamental insight connecting Dirichlet polynomial theory to the distribution of zeta zeros is the construction of a **zero-detecting polynomial**. If $\rho = \sigma + it$ is a nontrivial zero of $\zeta(s)$ with $0 < \sigma < 1$ and $t \in [T, 2T]$, then the approximate functional equation gives

$$\zeta(\sigma + it) = \sum_{n} \frac{\eta(\log n - \log x)}{n^{\sigma+it}} + O_A(T^{-A}),$$

where $\eta$ is a suitable smooth cutoff and $x$ is a parameter to be chosen. To isolate the contribution of small $n$, one multiplies by a short Möbius partial sum $\sum_{d \leq T^\varepsilon} \mu(d) d^{-\sigma-it}$, producing coefficients

$$a_n := \sum_{\substack{d \mid n \\ d \leq T^\varepsilon}} \mu(d)\, \eta\!\left(\log \frac{n}{d} - \log x\right) - \mathbf{1}_{n=1}.$$

The key property is that if $\zeta(\sigma + it) = 0$, then the resulting Dirichlet polynomial satisfies

$$\left|\sum_{n} \frac{a_n}{n^{\sigma + it}}\right| > \frac{1}{2}$$

for sufficiently large $T$. This Dirichlet polynomial "detects" zeros: it must be large at every zero of $\zeta$.

### 1.2 The Chain of Reductions

The passage from zero density estimates to large value estimates proceeds through three steps:

**Step 1 (Zero detection).** If $\rho = \sigma + it$ is a zero of $\zeta$ with $t \in [T, 2T]$, then a suitable Dirichlet polynomial $D(s) = \sum_n a_n n^{-s}$ satisfies $|D(\sigma + it)| > 1/2$.

**Step 2 (Dyadic decomposition).** Decompose $D$ over dyadic blocks: for $N$ ranging over dyadic values with $T^\varepsilon \ll N \ll T^{1+\varepsilon}$, write

$$D(\sigma + it) = \sum_{\text{dyadic } N} D_N(\sigma + it), \quad D_N(s) = \sum_{N/2 \leq n \leq N} \frac{a_n}{n^s}.$$

By the pigeonhole principle, for at least one dyadic block $N$, we must have $|D_N(\sigma + it)| \gg (\log T)^{-1}$.

**Step 3 (Large value theorem application).** Apply a large value estimate to bound the number of well-separated points $t$ at which $|D_N(\sigma + it)|$ can exceed a given threshold $V$. If $W \subset [0, T]$ is a 1-separated set with $|D_N(t)| \geq V$ for all $t \in W$, then the large value theorem bounds $|W|$.

The zero density estimate $N(\sigma, T)$ — defined as the number of zeros $\rho = \beta + i\gamma$ of $\zeta$ with $\beta \geq \sigma$ and $|\gamma| \leq T$ — is then obtained by summing over the dyadic blocks. Specifically, if $A(\sigma)$ denotes the **zero density exponent**, i.e., the infimum of all $A$ such that

$$N(\sigma, T) \ll T^{A(1 - \sigma) + o(1)},$$

then the connection to large values is captured by:

$$A(\sigma)(1-\sigma) \leq \sup_{\tau \geq 2} \frac{\mathrm{LV}(\sigma, \tau)}{\tau},$$

where $\mathrm{LV}(\sigma, \tau)$ is the large value exponent controlling the cardinality of the set where a Dirichlet polynomial of length $N = T^\tau$ exceeds $N^\sigma$.

### 1.3 The Halász–Montgomery Inequality

The underlying analytic engine is the **Halász–Montgomery inequality**, which is a Hilbert-space bound: given vectors $\xi, \varphi_1, \ldots, \varphi_R$ in a Hilbert space,

$$\sum_{r=1}^{R} |\langle \xi, \varphi_r \rangle|^2 \leq \|\xi\|^2 \cdot \left(\sup_{r} \sum_{s} |\langle \varphi_r, \varphi_s \rangle|\right).$$

Applied with $\xi$ being the coefficient vector and $\varphi_r(n) = n^{-it_r}$ for well-separated frequencies $t_r$, this connects the number of large values to the additive structure of the frequency set.

---

## 2. Large Value Estimates: History and State of the Art

### 2.1 The Central Problem

Let $D_N(t) = \sum_{n \sim N} b_n n^{it}$ with $|b_n| \leq 1$, and let $W \subset [0, T]$ be a 1-separated set on which $|D_N(t)| \geq V = N^\sigma$ for all $t \in W$. The fundamental question is: **what is the largest possible $|W|$?**

The Montgomery–Vaughan mean value theorem (1973) provides the starting point. In its classical form:

$$\int_0^T \left|\sum_{n \leq N} a_n n^{-it}\right|^2 dt = (T + O(N)) \sum_{n \leq N} |a_n|^2.$$

By Chebyshev's inequality, if $|D_N(t)| \geq V$ on a 1-separated set $W$, then

$$|W| \cdot V^2 \leq \int_0^T |D_N(t)|^2 \, dt \ll (T + N) \sum |a_n|^2 \ll (T+N)N,$$

giving the **mean value theorem bound**:

$$|W| \ll \frac{(T + N)N}{V^2} = (T+N)N^{1-2\sigma}.$$

### 2.2 Halász (1968) and Montgomery (1969)

Halász's 1968 work on mean values of multiplicative functions introduced the bilinear form technique that Montgomery adapted in 1969 for Dirichlet polynomials. In "Mean and large values of Dirichlet polynomials" (*Inventiones Math.* **8**, 334–345, 1969), Montgomery proved the first large value theorem going beyond the naive mean value bound:

**Montgomery's Large Value Estimate.** After removing $O(T^\beta)$ unit intervals from $[0, T]$, the remaining $t$ satisfy

$$\left|\sum_{N/2 \leq n \leq N} \frac{a_n}{n^{\sigma + it}}\right| \ll N^{1-\sigma}\left(T^{-\beta/2} + T^{1/4} N^{-1/2}\right) \log^{O(1)}(NT).$$

The key improvement over the pure mean value bound comes from the $T^{1/4}N^{-1/2}$ term, which provides savings when $N$ is large relative to $T$.

### 2.3 Montgomery's Large Value Conjecture

Montgomery conjectured that the optimal large value exponent satisfies

$$\mathrm{LV}(\sigma, \tau) \leq 2 - 2\sigma$$

for all $\sigma \in (1/2, 1)$ and all $\tau \geq 0$. In other words, the $L^2$ mean value theorem already captures the true large-value behavior. If true, this would imply the **density hypothesis**: $A(\sigma) \leq 2$ for all $\sigma$. Bourgain provided a counterexample showing that some $\varepsilon$-loss is necessary (related to Besicovitch sets of measure zero), but the conjecture is expected to hold up to such losses.

### 2.4 Huxley (1972): The Reflection Argument

Huxley's decisive contribution in "On the difference between consecutive primes" (*Invent. Math.* **15**, 164–170, 1972) combined two ingredients:

1. **Ingham's estimate** (1940): $A(\sigma) \leq 3/(2 - \sigma)$ for $1/2 \leq \sigma \leq 3/4$.
2. **Huxley's own estimate**: $A(\sigma) \leq 3/(3\sigma - 1)$ for $3/4 \leq \sigma \leq 1$.

The second estimate uses a **reflection argument** through the functional equation of $\zeta$, relating large values near $\sigma$ to those near $1 - \sigma$, together with a subdivision technique: a large value estimate on intervals of length $N$ automatically implies estimates on longer intervals by covering.

The combination gives the uniform bound

$$A(\sigma) \leq \frac{12}{5} = 2.4 \quad \text{for all } \sigma \in (1/2, 1),$$

with the supremum attained at $\sigma = 3/4$. This remained the best uniform bound for over fifty years.

### 2.5 Jutila (1977) and Heath-Brown (1979): Restricted Improvements

**Jutila** (*Acta Arith.* **32**, 55–62, 1977) proved a parametric family of large value estimates: for any positive integer $k$,

$$\mathrm{LV}(\sigma, \tau) \leq \max\left(2 - 2\sigma,\; \tau + (4 - 2/k) - (6 - 2/k)\sigma,\; \tau + (6 - 8\sigma)k\right).$$

Montgomery's conjecture holds in Jutila's framework whenever $\tau \leq \min\bigl((4 - 2/k)\sigma - (2 - 2/k),\; (8k-2)\sigma - 6k + 2\bigr)$.

**Heath-Brown** (*J. London Math. Soc.* (2) **19**, 221–232, 1979) proved

$$\mathrm{LV}(\sigma, \tau) \leq \max(2 - 2\sigma,\; 10 + \tau - 13\sigma),$$

so that Montgomery's conjecture holds when $\tau \leq 11\sigma - 8$. Heath-Brown also introduced bounds on **double zeta sums** — bilinear expressions involving $n^{it}$ — which would later prove crucial for Guth–Maynard.

These estimates improved $A(\sigma)$ in restricted $\sigma$-ranges but did not lower the supremum $\sup_\sigma A(\sigma) = 12/5$.

### 2.6 Bourgain (2000): The Density Hypothesis for Large $\sigma$

Bourgain verified the density hypothesis $A(\sigma) \leq 2$ for $\sigma \geq 25/32$, improving on Jutila's $\sigma > 11/14$. His argument introduced a **dichotomy**: either the set $W$ of large values has small additive energy (amenable to standard methods) or it correlates with large values of $\zeta$ itself (amenable to moment bounds). This dichotomy foreshadowed the energy-based case splitting in Guth–Maynard.

### 2.7 Guth–Maynard (2024): The Breakthrough

In "New large value estimates for Dirichlet polynomials" (arXiv:2405.20552, accepted *Annals of Mathematics* 2025), Guth and Maynard proved:

**Theorem 1.1 (Guth–Maynard Large Value Estimate).** Let $D_N(t) = \sum_{n \sim N} b_n n^{it}$ with $|b_n| \leq 1$, and let $W \subset [0, T]$ be 1-separated with $|D_N(t)| \geq V$ for $t \in W$, where $N \in [T^{2/3}, T]$. When $V$ is close to $N^{3/4}$, specifically for $V \in [N^{7/10+o(1)}, N^{8/10-o(1)}]$, the bound improves upon the classical mean value theorem bound $|W| \leq N^2 V^{-2+o(1)}$.

This leads to:

**Theorem 1.2 (Zero Density Estimate).** For $1/2 < \sigma < 1$,

$$N(\sigma, T) \ll T^{15(1-\sigma)/(3+5\sigma) + o(1)}.$$

Combining with Ingham's estimate for $\sigma \leq 7/10$ yields the headline result:

$$A(\sigma) \leq \frac{15}{3 + 5\sigma}, \quad \text{and hence} \quad \sup_\sigma A(\sigma) \leq \frac{30}{13} \approx 2.307,$$

with the supremum now attained at $\sigma = 7/10$ rather than $\sigma = 3/4$.

---

## 3. The Guth–Maynard Method in Detail

### 3.1 The Matrix-Theoretic Framework

The core innovation is to reformulate the large value problem in terms of **singular values of a matrix**. Define the $|W| \times N$ matrix $M_W$ with entries

$$M_{t,n} = w(n/N)\, n^{it},$$

where $w$ is a smooth bump function supported on $[1,2]$ and $t$ ranges over $W$.

**Lemma 4.1 (Singular Value Control).** If $|D_N(t)| \geq N^\sigma$ for all $t \in W$ and $|b_n| \leq 1$, then

$$|W| \lesssim N^{1-2\sigma}\, s_1(M_W)^2,$$

where $s_1(M_W)$ is the largest singular value of $M_W$.

This is immediate: $|W| N^{2\sigma} \leq \sum_{t \in W} |D_N(t)|^2 = \|M_W b\|^2 \leq s_1(M_W)^2 \|b\|^2 \leq s_1(M_W)^2 N$.

### 3.2 Bounding the Largest Singular Value via Traces

For any positive integer $r$, the largest singular value satisfies

$$s_1(M_W)^{2r} \leq \mathrm{tr}\!\left((M_W^* M_W)^r\right).$$

The classical mean value theorem corresponds to $r = 1$: $\mathrm{tr}(M_W^* M_W) = \sum_{t,n} |w(n/N)|^2 \approx |W| N$, recovering $|W| \leq N^{1-2\sigma} \cdot |W| N / |W| = N^{2-2\sigma}$.

For $r = 2$, one encounters a sum that can be handled by the mean value theorem for Dirichlet polynomials, but this gives no improvement over $r = 1$ in the critical range near $\sigma = 3/4$.

Guth and Maynard work with $r = 3$ (the **cubic trace**), expanding $\mathrm{tr}\!\left((M_W M_W^*)^3\right)$ via Poisson summation.

### 3.3 The Cubic Trace Decomposition

**Lemma 4.5 (Expansion of the cubic trace).** For $W$ a $T^\varepsilon$-separated set, the trace $\mathrm{tr}\!\left((M_W M_W^*)^3\right)$ decomposes into three sums $S_1 + S_2 + S_3$.

The three components are handled differently:

- **$S_1$** (Proposition 5.1): This is $O_\varepsilon(T^{-10})$, hence negligible.
- **$S_2$**: Bounded using **Heath-Brown's theorem on double zeta sums** (Theorem 1.6 in their paper). This provides strong estimates when the set $W$ has large additive energy.
- **$S_3$**: This is the main term and occupies the bulk of the paper (~34 of 48 pages). It contains the novel content.

### 3.4 The Role of Additive Energy

The **additive energy** of $W$ is defined as

$$E(W) := \#\{(w_1, w_2, w_3, w_4) \in W^4 : |w_1 + w_2 - w_3 - w_4| < 1\}.$$

The proof of Theorem 1.1 divides into cases based on $E(W)$:

1. **Low energy** ($E(W) \approx |W|^2$): Standard analytic estimates on $S_3$ suffice.
2. **High energy** ($E(W)$ close to $|W|^3$): Heath-Brown's theorem (Theorem 1.6) gives strong bounds on $S_2$.
3. **Intermediate energy**: This requires the deepest part of the argument — a strengthened bound on $S_3$ valid for a wider range of energies.

**Lemma 1.7 (Energy bound):** Under the hypotheses of Theorem 1.1 with $N \in [T^{2/3}, T]$,

$$E(W) \leq |W|^3 N^{1-2\sigma+o(1)} + |W|^2 N^{2-2\sigma+o(1)}.$$

This energy bound is itself a new result, providing the first non-trivial constraint on the additive structure of the set of large values.

### 3.5 Key Technical Maneuvers

Several distinctive features of the argument deserve emphasis:

- **Raising to the sixth power.** By working with $(M_W M_W^*)^3$ rather than lower powers, they access information about three-fold correlations in $W$, which connects to the additive energy.
- **Refusing to use stationary phase.** At a critical juncture, Guth and Maynard deliberately avoid simplifying a Fourier integral via stationary phase, instead keeping the full oscillatory integral. This preserves cancellation that would be lost in the asymptotic expansion.
- **Fourier analytic character.** The arguments are largely Fourier-analytic, connecting to the decoupling and restriction philosophy in harmonic analysis.

---

## 4. The Role of Vinogradov's Mean Value Theorem

### 4.1 The Main Conjecture (Now a Theorem)

The **Vinogradov mean value** $J_{s,k}(X)$ counts the number of integer solutions to

$$x_1^j + \cdots + x_s^j = x_{s+1}^j + \cdots + x_{2s}^j, \quad j = 1, \ldots, k,$$

with $1 \leq x_i \leq X$. The **main conjecture** asserts

$$J_{s,k}(X) \ll_\varepsilon X^{s + \varepsilon} + X^{2s - k(k+1)/2 + \varepsilon}.$$

This was proved for $k = 3$ by Wooley (2012, efficient congruencing) and for all $k \geq 4$ by Bourgain–Demeter–Guth (2016, $\ell^2$ decoupling).

### 4.2 Connection to the Guth–Maynard Framework

The VMVT enters the story in several ways:

1. **Weyl-sum estimates.** Heath-Brown showed that the VMVT implies improved $k$-th derivative bounds for Weyl sums, which in turn give new bounds for $\zeta(\sigma + it)$ throughout the critical strip. These feed into the exponent pair machinery underlying zero density estimates.

2. **Exponential sum refinements.** The Tao–Trudgian–Yang systematic approach (arXiv:2501.16779, January 2025) shows that the optimal resolution of the VMVT, combined with exponent pair calculus, propagates to improved zero density estimates in $\sigma$-ranges away from $3/4$.

3. **The Guth–Maynard $S_3$ estimate.** The bound on $S_3$ involves multi-linear exponential sums where VMVT-type cancellation is relevant. Specifically, when expanding the cubic trace, one encounters sums of the form $\sum n^{it}$ over constrained sets of integers $n$, and the additive structure of these sums connects to the same circle of ideas as the VMVT.

### 4.3 Untapped Potential

There are several directions where the VMVT could yield further improvements:

**Restricted VMVT.** If one could prove VMVT-type bounds restricted to integers $n$ satisfying specific congruence or divisibility conditions (e.g., smooth numbers, numbers in arithmetic progressions), these could be applied within the Guth–Maynard $S_3$ decomposition where the summation is naturally restricted by the Poisson summation step. Formally, one would need bounds of the form

$$\sum_{\substack{n_1^j + \cdots + n_s^j = m_1^j + \cdots + m_s^j \\ n_i, m_i \in \mathcal{A}}} 1 \ll_\varepsilon |\mathcal{A}|^{s+\varepsilon}$$

for structured subsets $\mathcal{A}$ of $[1, X]$.

**Higher-degree decoupling within the matrix framework.** The Guth–Maynard method uses the cubic trace ($r = 3$) because they cannot control higher traces ($r \geq 4$). The VMVT, via the decoupling machinery of Bourgain–Demeter–Guth, gives optimal control of $r$-th power moments of exponential sums. If one could reformulate the matrix trace bounds in a way that interfaces directly with decoupling estimates at the level of $r = 4$ or $r = 5$, one might overcome the current barrier. The obstacle is that the matrix entries involve $\log n$ (from $n^{it} = e^{it \log n}$) rather than polynomial phases, so standard decoupling does not directly apply.

---

## 5. Moment vs. Large Value Duality: Where the Slack Lies

### 5.1 The Basic Mechanism

The passage from moments to large values goes through **Markov's inequality**: if $\int_0^T |D(t)|^{2k} dt \leq M_{2k}$, then

$$|\{t \in [0,T] : |D(t)| \geq V\}| \leq \frac{M_{2k}}{V^{2k}}.$$

The $L^2$ moment gives $M_2 \ll (T+N)N$ (Montgomery–Vaughan), yielding $|W| \ll (T+N)N/V^2$, which is the mean value theorem bound. Higher moments could in principle give tighter bounds, but computing $M_{2k}$ for $k \geq 2$ requires controlling higher correlations.

### 5.2 Where Slack Appears in Guth–Maynard

The Guth–Maynard argument does **not** proceed through a moment bound at all — it goes through **singular values**, which is a subtly different approach. The key slack appears at several points:

1. **From singular values to cardinality.** Lemma 4.1 uses $|W| \leq N^{1-2\sigma} s_1^2$. But $s_1$ is the **largest** singular value, while the Dirichlet polynomial $D_N(t) = M_W b$ for a specific coefficient vector $b$. If $b$ does not align with the top singular vector, the bound is wasteful. In principle, one could use spectral information about $b$ relative to $M_W$ to obtain sharper bounds, but this requires knowledge of the singular vectors, which seems difficult.

2. **From trace to largest singular value.** The bound $s_1^{2r} \leq \mathrm{tr}((M_W^* M_W)^r)$ is an instance of the power-mean inequality: the $\ell^\infty$ norm of the singular value sequence is bounded by the $\ell^r$ norm. This incurs a factor of (rank)$^{1/r}$. For $r = 3$, the trace captures some information about the singular value distribution, but not all. If the singular values are spread out (rather than concentrated at $s_1$), the trace method overestimates $s_1$.

3. **The energy dichotomy.** The case split based on $E(W)$ introduces slack at the boundary between regimes. The intermediate-energy case requires balancing two different bounds, and the optimal balance point determines the final exponent. Any improvement to either the low-energy or high-energy bound would shift this balance point.

### 5.3 Could a Direct Moment Approach Do Better?

A direct moment approach would compute $\int_0^T |D_N(t)|^{2k} dt$ for some $k > 1$ and apply Markov. The advantage would be avoiding the singular-value step entirely, which introduces the slack in point (1) above. The disadvantage is that for $k \geq 2$, the moment involves higher correlations:

$$\int_0^T |D_N(t)|^{2k} dt = \sum_{n_1 \cdots n_k = m_1 \cdots m_k} b_{n_1} \cdots b_{n_k} \overline{b_{m_1} \cdots b_{m_k}} \cdot (T + O(\cdots)),$$

and the constraint $n_1 \cdots n_k = m_1 \cdots m_k$ (multiplicative energy) is arithmetically very different from the additive energy $E(W)$ that appears in the Guth–Maynard framework.

The fourth moment ($k = 2$) involves the **divisor function** and is closely related to the **fourth moment of $\zeta$**. The sixth moment ($k = 3$) is connected to $\int |\zeta|^6$, which remains poorly understood. Thus a direct moment approach at $k = 3$ would require essentially solving the sixth moment problem for $\zeta$, which is far beyond current technology.

The Guth–Maynard method cleverly sidesteps this by working with the matrix trace (which involves **additive** rather than multiplicative correlations) and using the flexibility to choose $r = 3$ without needing the full sixth moment of $\zeta$.

---

## 6. The Twisted Fourth Moment and Potential Synergies

### 6.1 Classical Results

Heath-Brown (1979) established the asymptotic formula for the fourth moment $\int_0^T |\zeta(1/2 + it)|^4 dt$. The **twisted** fourth moment, where one multiplies by a Dirichlet polynomial, has been studied extensively:

- **Hughes–Young** (2010): Asymptotics for $\int_0^T |\zeta(1/2+it)|^4 |A(1/2+it)|^2 dt$ when $A$ is a Dirichlet polynomial of length $T^{1/11-\varepsilon}$.
- **Bettin–Bui–Li–Radziwiłł**: Extended the length to $T^{1/4-\varepsilon}$, using Watt's theorem on Kloosterman fractions.
- **Page (2025)**: The amplified fourth moment $\int_0^T |\zeta|^4 |A|^4 dt$ with $A$ of length $T^{1/4-\varepsilon}$, yielding unconditional lower bounds for the sixth and eighth moments.

### 6.2 Potential Combination with Guth–Maynard

The twisted fourth moment controls the **average** behavior of $\zeta$ weighted by a Dirichlet polynomial, while Guth–Maynard controls the **exceptional set** where a Dirichlet polynomial is large. There are two potential synergies:

**Approach A: Refined zero detection.** In the zero-detection step, one could use a more sophisticated polynomial $D(s)$ that incorporates information from the fourth moment. Specifically, if $A(s)$ is a "mollifier" for $\zeta$ (i.e., $A(s) \approx 1/\zeta(s)$ in a suitable sense), then $\zeta(s) A(s) \approx 1$, and a zero of $\zeta$ forces $|\zeta(s) A(s)|$ to be small rather than $|D(s)|$ to be large. This reverses the logic: instead of detecting where $D$ is large, one detects where $\zeta A$ deviates from its expected value, and the twisted fourth moment controls the average deviation.

**Approach B: Hybrid large value / moment bounds.** One could use the twisted fourth moment to bound the contribution from zeros with "typical" behavior (where $\zeta$ near $\rho$ behaves like a generic point on the critical line) and reserve the Guth–Maynard large value estimate for "exceptional" zeros. This hybrid approach would be most effective if the fourth moment gives savings in a $\sigma$-range complementary to where Guth–Maynard is strongest ($\sigma \approx 3/4$).

### 6.3 A Specific Hybrid Bound

Consider the mixed estimate: for a 1-separated set $W$ with $|D_N(t)| \geq V$ for $t \in W$,

$$|W| V^4 \leq \sum_{t \in W} |D_N(t)|^4 \leq \left(\sum_{t \in W} |\zeta(\sigma+it) A(\sigma+it)|^4\right) \cdot \sup_t \left|\frac{D_N(t)}{\zeta(\sigma+it) A(\sigma+it)}\right|^4.$$

If $D_N$ is related to $\zeta$ via the zero-detecting construction, the ratio $D_N/(\zeta A)$ is bounded, and the sum over $W$ is controlled by the (discrete) twisted fourth moment. When $|W|$ is large, the discrete sum over $W$ is comparable to the continuous integral, which is known. This gives

$$|W| \ll V^{-4} T^{1+\varepsilon} \cdot (\text{correction from polynomial length}).$$

The constraint is the length of the mollifier/amplifier $A$: current technology allows $T^{1/4-\varepsilon}$, which limits the applicability to $N \leq T^{1/4}$.

---

## 7. Recent Systematic Developments (2025)

### 7.1 The Tao–Trudgian–Yang Framework

In January 2025, Tao, Trudgian, and Yang (arXiv:2501.16779) introduced the **Analytic Number Theory Exponent Database (ANTEDB)**, a systematic framework encoding:

- Growth exponents $\mu(\sigma)$ of the Riemann zeta function
- Exponent pairs $(k, \ell)$ for bounding exponential sums
- Zero density exponents $A(\sigma)$
- Zero additive energy exponents

By systematically optimizing the relationships between these exponents using computer search, they obtained four new exponent pairs and several new zero density estimates. Their work confirms that the Guth–Maynard estimate is the binding constraint at $\sigma = 7/10$, while other classical estimates bind at other values of $\sigma$.

### 7.2 Extensions to $L$-functions

The Guth–Maynard method has been extended to Dirichlet $L$-functions by B. Chen and others (arXiv:2507.08296), who adapted the matrix framework to Dirichlet polynomials twisted by characters. The key innovation involves handling sums with GCD twists that arise from the character orthogonality relations.

---

## 8. Concrete Proposal: Improving the 30/13 Exponent

### 8.1 Identifying the Bottleneck

The exponent $30/13$ arises from $A(7/10) = 15/(3 + 5 \cdot 7/10) = 15/(6.5) = 30/13$. The supremum moved from $\sigma = 3/4$ (where Huxley's bound was 12/5) to $\sigma = 7/10$ (where Ingham's bound $3/(2-\sigma) = 3/1.3 = 30/13$ now matches the Guth–Maynard bound). Thus, **the bottleneck is now Ingham's classical estimate at $\sigma = 7/10$**, not the Guth–Maynard large value theorem.

### 8.2 Strategy: Improve Ingham's Bound Near $\sigma = 7/10$

**Ingham's bound** $A(\sigma) \leq 3/(2-\sigma)$ arises from the mean value theorem applied to the zero-detecting polynomial with $N \approx T^{2/(2-\sigma)}$, combined with the classical $L^2$ bound. To improve this near $\sigma = 7/10$, one needs a better large value estimate when $\tau = 2/(2-\sigma)$ is close to $20/13$.

**Specific Proposal.** Extend the Guth–Maynard singular value method to work in the range $N \approx T^{20/13}$, which requires $N > T$ (the "long Dirichlet polynomial" regime). The current Guth–Maynard result assumes $N \in [T^{2/3}, T]$, so the constraint $N \leq T$ is binding. In the regime $N > T$, the Montgomery–Vaughan mean value theorem gives $(T + N)N \approx N^2$ (the $N$ term dominates), and the mean value bound is $|W| \leq N^{2-2\sigma}$. The Guth–Maynard matrix $M_W$ has more columns than rows when $N > T > |W|$, and the singular value analysis needs modification.

**What would be needed:** A version of the cubic trace expansion (Lemma 4.5) valid when $N > T$, with appropriate Poisson summation estimates. The $S_3$ bound would need to account for the fact that the Poisson dual has longer support when $N > T$. One would also need a version of Heath-Brown's double zeta sum estimate in this range.

**Expected outcome:** If the Guth–Maynard bound could be extended to give $A(\sigma) < 3/(2-\sigma)$ in a neighborhood of $\sigma = 7/10$, the supremum would decrease below $30/13$. Even a modest improvement like $A(7/10) \leq 30/13 - \delta$ for small $\delta > 0$ would be significant.

### 8.3 Alternative Strategy: Exploit the Energy–Cardinality Tradeoff

An alternative approach targets the intermediate-energy case in the Guth–Maynard argument directly. Currently, the bound on $S_3$ in the intermediate-energy regime involves a balance between the additive energy $E(W)$ and $|W|$. One could attempt to:

1. **Upgrade the energy bound (Lemma 1.7)** using the VMVT. The additive energy of the set $W$ where $|D_N(t)| \geq V$ may satisfy stronger constraints when the coefficients $b_n$ have multiplicative structure (as they do for the zero-detecting polynomial). A VMVT-type estimate for "multiplicatively structured" additive energy could give $E(W) \ll |W|^{2+\varepsilon} N^{1-2\sigma+\varepsilon}$ with improved implicit constants.

2. **Use the twisted fourth moment to bound $S_2$.** Currently, $S_2$ is bounded by Heath-Brown's theorem. If instead one uses the twisted fourth moment of $\zeta$ to control the bilinear sums in $S_2$ when the Dirichlet polynomial has arithmetic structure (as it does in the zero-density application), one might obtain savings of the form $T^{-\delta}$ for some $\delta > 0$.

3. **Combine (1) and (2)** to shift the balance point in the energy dichotomy, obtaining a bound like

$$A(\sigma) \leq \frac{15}{3 + 5\sigma} - c(\sigma)$$

for a positive function $c(\sigma)$ in a neighborhood of $\sigma = 3/4$.

### 8.4 Quantitative Expectations

If strategy 8.3 succeeds with even modest savings in $S_2$ (say a factor $T^{-1/100}$), the balance in the intermediate-energy case would shift, and one could expect an improvement in the zero density exponent at $\sigma = 3/4$ from $30/13 \approx 2.307$ to roughly $2.30$. While numerically small, this would demonstrate that the Guth–Maynard framework has room for further optimization and would represent the first improvement to the $30/13$ exponent.

More ambitiously, if the VMVT could be brought to bear on the $r = 4$ trace (avoiding the current obstruction of uncontrolled quartic terms), one might achieve $A(\sigma) \leq 12/(2 + 4\sigma)$ near $\sigma = 3/4$, giving $A(3/4) = 12/5 \cdot 3/5 = 36/25 = 1.44$, far below $30/13$. This would require a genuinely new idea for the quartic trace and remains highly speculative.

---

## 9. Conclusion

The Guth–Maynard breakthrough represents a paradigm shift in the study of large values of Dirichlet polynomials: from moment-based methods to **singular-value methods**. The key innovation — using the cubic trace of a matrix of exponential phases, decomposed via Poisson summation into terms controlled by additive energy — opens several avenues for further progress.

The most promising near-term direction is to combine the Guth–Maynard framework with:
1. The twisted fourth moment technology (to handle the $S_2$ term more efficiently when the Dirichlet polynomial has arithmetic structure), and
2. Restricted VMVT estimates (to improve the $S_3$ bound in the intermediate-energy case).

The longer-term goal of improving the trace power from $r = 3$ to $r = 4$ would require a fundamentally new approach to the quartic trace expansion, but if achieved, would bring the zero density exponent significantly closer to the density hypothesis.

---

## References

1. G. Halász, "Über die Mittelwerte multiplikativer zahlentheoretischer Funktionen," *Acta Math. Acad. Sci. Hungaricae* **19** (1968), 365–403.
2. H. L. Montgomery, "Mean and large values of Dirichlet polynomials," *Invent. Math.* **8** (1969), 334–345.
3. H. L. Montgomery and R. C. Vaughan, "The large sieve," *Mathematika* **20** (1973), 119–134.
4. M. N. Huxley, "On the difference between consecutive primes," *Invent. Math.* **15** (1972), 164–170.
5. M. Jutila, "Zero-density estimates for $L$-functions," *Acta Arith.* **32** (1977), 55–62.
6. D. R. Heath-Brown, "Zero density estimates for the Riemann zeta-function and Dirichlet $L$-functions," *J. London Math. Soc.* (2) **19** (1979), 221–232.
7. D. R. Heath-Brown, "The fourth power moment of the Riemann zeta-function," *Proc. London Math. Soc.* (3) **38** (1979), 385–422.
8. J. Bourgain, "On large values estimates for Dirichlet polynomials and the density hypothesis for the Riemann zeta function," *Internat. Math. Res. Notices* (2000), 133–146.
9. T. D. Wooley, "Vinogradov's mean value theorem via efficient congruencing," *Ann. of Math.* **175** (2012), 1575–1632.
10. J. Bourgain, C. Demeter, and L. Guth, "Proof of the main conjecture in Vinogradov's Mean Value Theorem for degrees higher than three," *Ann. of Math.* **184** (2016), 633–682.
11. C. Hughes and M. Young, "The twisted fourth moment of the Riemann zeta function," *J. reine angew. Math.* **641** (2010), 203–236.
12. S. Bettin, H. Bui, X. Li, and M. Radziwiłł, "A quadratic divisor problem and moments of the Riemann zeta-function," *J. Eur. Math. Soc.* (to appear).
13. L. Guth and J. Maynard, "New large value estimates for Dirichlet polynomials," *Ann. of Math.* (to appear), arXiv:2405.20552.
14. T. Tao, T. Trudgian, and A. Yang, "New exponent pairs, zero density estimates, and zero additive energy estimates: a systematic approach," arXiv:2501.16779 (2025).
15. T. Tao, "A computation-outsourced discussion of zero density theorems for the Riemann zeta function," *What's New* blog, July 7, 2024.
16. T. Tao, "254A, Notes 6: Large values of Dirichlet polynomials, zero density estimates, and primes in short intervals," *What's New* blog, February 13, 2015.
