# Hybrid Moment-Density Methods: Random Matrix Theory Meets Zero Density Estimates

**Author:** Student 4 (RMT-Analytic Number Theory Interface)
**Date:** 2026-02-23
**Status:** Preliminary Research Report

---

## Abstract

This report investigates the largely unexplored interface between two powerful frameworks for understanding the zeros of the Riemann zeta function: **zero density estimates** (from the Ingham-to-Guth-Maynard tradition) and **random matrix theory** (from the Montgomery-Odlyzko-Keating-Snaith paradigm). Although these programs have developed largely independently over the past fifty years, several deep structural connections exist that could, if properly exploited, yield improved bounds. We examine the Gonek-Hughes-Keating hybrid formula, the Fyodorov-Hiary-Keating extreme value conjecture, the Conrey-Snaith ratios recipe, and the phenomenon of RMT-predicted zero repulsion, assessing each as a potential source of new constraints on off-line zero density. We also propose concrete numerical experiments using the existing computational infrastructure in the `/home/user/riemann/` codebase.

---

## 1. The Gonek-Hughes-Keating Hybrid Formula and the Zero-Detection Framework

### 1.1 The Hybrid Decomposition

The central innovation of Gonek, Hughes, and Keating (Duke Math. J. 136(3), 507-549, 2007; arXiv: math/0511182) is a representation of the Riemann zeta function that interpolates between its Euler product and its Hadamard product. For a parameter X with 2 <= X <= T, they write

    zeta(s) ~ P_X(s) * Z_X(s)

where P_X(s) is a truncated Euler product over primes p <= X:

    P_X(s) = prod_{p <= X} (1 - p^{-s})^{-1} * (smoothing corrections)

and Z_X(s) is a partial Hadamard product over zeros rho with |Im(rho) - t| <= X:

    Z_X(s) ~ prod_{|gamma - t| < X} (1 - s/rho) * (normalizing factors)

The parameter X mediates between two extremes: when X is small, P_X is close to 1 and Z_X carries essentially all the information (the "spectral" side); when X is large, P_X captures most of the arithmetic and Z_X becomes smooth. Near height T, the Hadamard component Z_X involves roughly (X log T)/(2 pi) zeros, while the Euler component involves primes up to X.

The key achievement of this decomposition for moment theory is that it **naturally produces the arithmetic factor a(k)** in the Keating-Snaith conjecture. In the original Keating-Snaith approach, the moment of zeta was modeled by the moment of a CUE characteristic polynomial, but the arithmetic factor (an Euler product over primes) had to be inserted by hand. In the GHK hybrid model, the 2k-th moment decomposes as:

    (1/T) int_0^T |zeta(1/2+it)|^{2k} dt ~ (moment of P_X) * (moment of Z_X)

The moment of P_X yields the arithmetic factor a(k) via standard probabilistic arguments on the Euler product, while the moment of Z_X is modeled by CUE characteristic polynomials, yielding the RMT factor g(k) = G(1+k)^2/G(1+2k). For k = 1 and k = 2, all steps in this approach have been made rigorous.

### 1.2 Connection to the Halasz-Montgomery Zero-Detection Method

The Halasz-Montgomery method for proving zero density estimates proceeds as follows. Suppose rho = beta + i*gamma is a zero of zeta with beta > 1/2. One constructs a **zero-detecting polynomial** D(s) = sum_{n <= N} a_n n^{-s} such that |D(rho)| is guaranteed to be large (typically >> N^{sigma - 1/2} for some range of sigma). Then a **large value theorem** for Dirichlet polynomials --- bounding the measure of the set where |D(s)| exceeds a threshold --- translates into an upper bound on the number of zeros.

The classical choices for D(s) involve partial sums of zeta or of its logarithmic derivative. The Guth-Maynard breakthrough (arXiv: 2405.20552, 2024) achieved the exponent 30/13 by proving new large value estimates for such Dirichlet polynomials via decoupling methods from harmonic analysis.

**Key Question: Can the GHK hybrid splitting improve the zero-detection polynomial?**

Consider replacing the standard zero-detecting polynomial by the hybrid decomposition evaluated at a zero rho:

    0 = zeta(rho) ~ P_X(rho) * Z_X(rho)

Since zeta vanishes at rho, either P_X(rho) is small or Z_X(rho) is small (or both). This creates a **dichotomy**:

- **Case 1 (P_X(rho) is small):** The partial Euler product is unusually small at a point 1/2 + epsilon + i*gamma off the critical line. By the probabilistic model for Euler products (log P_X is approximately Gaussian by Selberg's CLT), this is a rare event whose probability can be bounded explicitly. The probability that |P_X(sigma + it)| < delta for sigma > 1/2 is controlled by the left tail of a Gaussian with variance ~ (1/2) log log X. This gives:

      P(|P_X(sigma+it)| < delta) << delta^c * (log X)^{-A}

  for computable constants c, A depending on sigma.

- **Case 2 (Z_X(rho) is small):** The partial Hadamard product is small. But Z_X is modeled by CUE characteristic polynomials, and the value distribution of CUE characteristic polynomials near their zeros is well understood via RMT. The probability that a CUE polynomial is small at a point slightly off the unit circle (analogous to sigma > 1/2) decays in a specific way controlled by the GUE/CUE kernel.

The RMT prediction for Case 2 is particularly interesting. For an N x N CUE matrix U, the characteristic polynomial Z_N(e^{i*theta}) has zeros on the unit circle with mutual repulsion governed by the Vandermonde determinant. The probability that |Z_N(r * e^{i*theta})| is small for r = 1 - delta (slightly inside the unit circle) is controlled by the distance to the nearest eigenvalue. By the CUE correlation functions, the probability that the nearest eigenvalue is within distance epsilon of the test point is ~ epsilon * N, leading to:

    P(|Z_N(r * e^{i*theta})| < eta) << eta^{1+o(1)} * N^{-c(delta)}

for some exponent c(delta) > 0 that depends on the distance from the unit circle.

**Proposed approach:** Use the hybrid dichotomy to split the zeros of zeta into "P_X-small" and "Z_X-small" classes. Bound the first class using Selberg-type probabilistic arguments on Euler products. Bound the second class using RMT predictions for CUE characteristic polynomial values near (but not on) the unit circle. If the RMT predictions provide better control than the generic large value theorems for Dirichlet polynomials, this could improve zero density estimates in certain sigma-ranges.

### 1.3 The Critical Bottleneck

The main obstacle is that the GHK hybrid formula is currently proved only in an averaged sense (for moments) and as a pointwise approximation only under the Riemann Hypothesis and certain technical conditions on X relative to T. Making the dichotomy argument rigorous for zero detection would require either:

(a) An unconditional pointwise hybrid formula with sufficiently small error, or
(b) A way to work with the hybrid splitting at the level of mean values rather than pointwise.

Option (b) seems more tractable and connects to the "mean value approach" to zero density, where one bounds sum_{rho} |D(rho)|^2 rather than looking at individual zeros.

---

## 2. Selberg's Central Limit Theorem, GUE, and Off-Line Zeros

### 2.1 Selberg's CLT and Its RMT Enhancement

Selberg proved (1946, published 1991) that for t chosen uniformly from [T, 2T]:

    log |zeta(1/2 + it)| / sqrt((1/2) log log T) --> N(0, 1)

in distribution. That is, the logarithm of zeta on the critical line is approximately Gaussian with variance (1/2) log log T. Radziwill and Soundararajan (2017) gave a new, elegant proof of this result.

The GUE conjecture predicts more: the **full local statistics** of zeros, including all n-point correlation functions, match those of GUE random matrices. Keating and Snaith showed that this implies the determinant of a CUE matrix has a log-normal distribution matching Selberg's theorem, but the GUE/CUE model also predicts:

1. The complete moment generating function of log |zeta(1/2+it)|, not just the Gaussian approximation
2. Large deviation estimates beyond the Gaussian regime
3. Joint distributions of zeta at multiple points

Recent work by Fazzari et al. (arXiv: 2507.04150, 2025) establishes a **joint CLT** for log zeta weighted by linear statistics of the zeros, showing that the logarithm of zeta and the zero correlations behave approximately like two independent random variables.

### 2.2 From Value Distribution to Zero Density

The connection between value distribution on the critical line and off-line zeros proceeds via classical complex analysis tools:

**Jensen's Formula Route:** For a disk of radius R centered at 1/2 + it, Jensen's formula gives:

    sum_{|rho - (1/2+it)| < R} log(R/|rho - (1/2+it)|) = (1/2pi) int_0^{2pi} log|zeta(1/2 + it + R*e^{i*theta})| d*theta - log|zeta(1/2+it)|

The left side counts zeros (weighted by proximity). The right side involves the mean value of log|zeta| on a circle that extends into the region sigma > 1/2. If we have tight control on the distribution of log|zeta(1/2+it)| from Selberg's CLT (and especially its tails), this constrains the right side, and hence constrains the number of off-line zeros.

**The GUE Enhancement:** The key advantage of using the full GUE prediction (rather than just the Gaussian approximation) is in the **tails**. Selberg's CLT says log|zeta| is approximately Gaussian, but this breaks down for values larger than ~ (log log T)^{1/2+epsilon}. The GUE model, via the Keating-Snaith moment conjecture, predicts:

    P(|zeta(1/2+it)| > V) ~ exp(-c * (log V)^2 / log log T)

for V in the range exp(c_1 * sqrt(log log T)) < V < exp(c_2 * log log T). This is **sub-Gaussian** decay in the far tail, which is better than what the bare CLT gives.

**Hadamard Three-Lines Route:** If we know that |zeta(1/2+it)| <= M(t) on the critical line, and |zeta(1+it)| << log t on the 1-line, then the Hadamard three-lines theorem (or Phragmen-Lindelof) bounds |zeta(sigma+it)| for 1/2 < sigma < 1:

    log|zeta(sigma+it)| <= (1-sigma) * log M(t) + (sigma - 1/2) * log log t + lower order

The zero-free region of zeta near 1 + it then constrains how large zeta can be slightly to the right of the critical line. Combining this with a zero density framework: if sigma is such that zeta(sigma+it) = 0, then the function is oscillating, and Jensen's formula applied locally constrains the density of such oscillations.

### 2.3 Quantitative Estimate

Suppose the GUE conjecture gives us that on the critical line, for most t:

    |zeta(1/2+it)| <= exp((1+epsilon) * sqrt((1/2) * log log T * log 2))

(This is the "typical maximum" in the Gaussian model.) By three-lines interpolation with the bound on the 1-line, for sigma = 1/2 + delta:

    log|zeta(sigma+it)| <= (1 - 2*delta) * sqrt((1/2) * log log T * log 2) + 2*delta * log log T

For zeros at sigma = 1/2 + delta, the local zero density is constrained by how fast the function oscillates, which is controlled by its size. This suggests:

    N(1/2 + delta, T) << T * exp(-c * delta^2 * log log T)

for delta in a suitable range, which is much stronger than the polynomial bounds T^{A(1-sigma)} from classical zero density estimates --- but only for delta << 1/sqrt(log log T), i.e., very close to the critical line. This is exactly the regime where classical density estimates are weakest, suggesting complementary information.

---

## 3. The Fyodorov-Hiary-Keating Conjecture and Zero Density

### 3.1 The Conjecture and Its Resolution

Fyodorov, Hiary, and Keating (Phys. Rev. Lett. 2012; arXiv: 1202.4713) conjectured, based on the freezing transition in log-correlated random energy models, that:

    max_{t in [T, T+1]} |zeta(1/2+it)| ~ (log T) / (log log T)^{3/4}

More precisely, for tau uniform on [T, 2T]:

    max_{|h| <= 1} |zeta(1/2 + i*tau + ih)| * (log log T)^{3/4} / log T

converges in distribution to a random variable with tail decay ~ y * e^{-2y}.

This conjecture has now been resolved in full by Arguin, Bourgade, and Radziwill:
- **Part I (Upper bound):** arXiv: 2007.00988 (2020) proved the upper bound in strong form.
- **Part II (Lower bound):** arXiv: 2307.00982 (2023, updated 2024) proved the matching lower bound, confirming the tail decay is ~ y * e^{-2y}.

The proof identifies an **approximate branching random walk** (tree structure) in the Fourier decomposition of the zeta function, connecting it to the theory of extremes for log-correlated Gaussian fields.

### 3.2 Implications for Off-Line Zero Density

The FHK conjecture, now theorem, constrains off-line zeros through several mechanisms:

**Mechanism 1: Three-Lines Interpolation from Maximum Bounds.**

If we know that typically max_{[T,T+1]} |zeta(1/2+it)| ~ L := log T / (log log T)^{3/4}, and on the 1-line we have the standard bound |zeta(1+it)| << log T, then the Hadamard three-lines theorem gives for sigma in (1/2, 1):

    max_{[T,T+1]} |zeta(sigma+it)| << L^{2(1-sigma)} * (log T)^{2sigma-1}
                                     = (log T)^1 * (log log T)^{-(3/2)(1-sigma)}

This is essentially the Lindelof hypothesis on average (up to the (log log T) correction), and it is much stronger than the pointwise convexity bound. In particular, for sigma = 3/4:

    max_{[T,T+1]} |zeta(3/4+it)| << (log T) * (log log T)^{-3/8}

**Mechanism 2: Moment Constraints via Freezing.**

The FHK freezing transition predicts that the "free energy" of the zeta function:

    F_beta(T) = (1/beta) * log( (1/T) int_T^{2T} |zeta(1/2+it)|^{2*beta} dt )

satisfies F_beta = beta + 1/beta for beta < 1 (the "unfrozen" phase) and F_beta = 2 for beta >= 1 (the "frozen" phase). The critical point beta_c = 1 corresponds to the moment where the integral is dominated by the extreme values. The Keating-Snaith moment conjecture gives the unfrozen phase (2k-th moment ~ C_k * (log T)^{k^2}), while the frozen phase (k > 1) requires understanding the extreme value statistics. The exponent k^2 transitions to 2k - 1 in the frozen regime.

For zero density applications, the key insight is: the frozen phase moments are **weaker** than the naive Keating-Snaith extrapolation would suggest. This means that the contribution of extreme values of zeta to mean values is capped. When applied via Halasz-Montgomery, this means that the zero-detecting polynomial D(s), whose mean square is estimated using zeta moments, has a better-controlled contribution from extreme values than previously thought. This could sharpen the large value estimates that feed into zero density bounds.

**Mechanism 3: Local Maximum Controls Local Zero Count.**

The Riemann-von Mangoldt formula gives approximately (1/(2pi)) * log(T/(2pi)) zeros per unit interval near height T. If the maximum of |zeta| on [T, T+1] is M, then by Jensen's formula applied to the disk of radius 1 centered at 1/2 + iT, the number of zeros in the disk is << log M + log log T. The FHK result says M ~ log T / (log log T)^{3/4}, so:

    #{zeros in disk of radius 1} << log log T

This is consistent with the expected count from the smooth zero-counting function. But more importantly, for **off-line** zeros at sigma = 1/2 + delta, the contribution to Jensen's formula is weighted by log(1/delta), meaning that off-line zeros at distance delta from the critical line cost a factor log(1/delta) in the Jensen budget. Combined with the FHK bound on the maximum, this means:

    #{zeros with Re(rho) > 1/2 + delta in [T, T+1]} << log log T / log(1/delta)

for delta >> 1/T. This is a local (per unit interval) zero density estimate. Summing over T unit intervals and comparing with Guth-Maynard's global estimate N(sigma, T) << T^{(30/13)(1-sigma)+o(1)}, the FHK-derived local bound is stronger when delta is not too small, giving complementary information.

---

## 4. Farmer-Gonek-Hughes Zero Repulsion and RMT-Predicted Density Exponents

### 4.1 Zero Repulsion from the GUE/CUE Model

The phenomenon of zero repulsion refers to the observation, both from RMT and from numerical computation, that the zeros of zeta tend to **repel each other** and also to **stay close to the critical line**. In the RMT framework, this has multiple manifestations:

**Vertical repulsion (Montgomery pair correlation):** The probability of finding two zeros within normalized distance epsilon of each other is ~ epsilon^2 (from the sin^2(pi*x)/(pi*x)^2 pair correlation kernel). This quadratic vanishing is the hallmark of GUE-type repulsion.

**Horizontal repulsion (off-line zeros):** If rho = 1/2 + delta + i*gamma is a zero with delta > 0, then the nearby zeros on the critical line exert "repulsion" in the horizontal direction. The RMT model for this comes from studying where the zeros of the derivative zeta'(s) fall: Mezzadri (arXiv: math-ph/0207044, 2002) showed that for CUE characteristic polynomials, the fraction of roots of Z'(U,z) at distance ~ delta/N from the unit circle has a specific limiting distribution.

**The Farmer-Gonek-Hughes picture:** Building on the hybrid formula, Farmer, Gonek, Hughes, and collaborators studied mean values of products of logarithmic derivatives of zeta near the critical line and their connections to zero correlations. Their work with Lester established that correlations of zeros at different heights, as measured by mean values of zeta'/zeta, are controlled by a combination of RMT predictions (for the local spacing) and arithmetic (from the primes).

### 4.2 Quantifying Repulsion for Zero Density

The key question is: what zero density exponent does RMT predict?

Let us work within the CUE model with N ~ log T / (2 pi). An "off-line zero" at distance delta from the critical line corresponds, in the CUE model, to a zero of Z_N(z) at distance ~ delta * (2 pi / log T) from the unit circle. But the zeros of Z_N(z) (as a polynomial in z) all lie exactly on the unit circle! This is the CUE analog of the Riemann Hypothesis.

Thus, within the pure CUE model, there are **no off-line zeros at all**, which "predicts" N(sigma, T) = 0 for sigma > 1/2 --- i.e., the Riemann Hypothesis. This is not directly useful as a quantitative bound.

However, we can ask a more refined question: **what does the CUE model predict for near-misses?** That is, how close can Z_N(z) come to zero at a point z with |z| = 1 - delta/N? The answer, from the theory of characteristic polynomial values, is:

    P(|Z_N(r * e^{i*theta})| < epsilon) ~ epsilon^2 / delta^2   for epsilon << delta

The epsilon^2 (rather than epsilon) reflects the quadratic repulsion of eigenvalues. Translating back to zeta:

    P(|zeta(1/2 + delta + it)| < epsilon) ~ epsilon^2 * (log T)^{-2} / delta^2

This predicts that the density of "approximate zeros" (points where |zeta(sigma+it)| < epsilon) decreases like 1/delta^2 as delta increases from 0. Using the heuristic that a zero at sigma = 1/2 + delta can be detected as a point where |zeta| is very small in a delta-neighborhood, this would predict:

    N(1/2 + delta, T) << T / (delta^2 * (log T)^2)

For delta = 1 - sigma (so sigma close to 1/2), this gives:

    N(sigma, T) << T * (1-sigma)^{-2} * (log T)^{-2}

This is **much stronger** than the density hypothesis (which gives T^{2(1-sigma)}) for fixed sigma close to 1/2, but weaker for sigma close to 1. The two predictions cross near sigma ~ 1 - 1/log T, which is exactly the transition region where classical density estimates change character.

### 4.3 A Conditional Density Estimate from RMT

We can formulate a precise conditional result:

**Conjecture (RMT-Predicted Zero Density):** Assume the GUE conjecture for the local statistics of zeta zeros, and assume the Gonek-Hughes-Keating hybrid formula with parameter X = T^{alpha}. Then for 1/2 < sigma < 1:

    N(sigma, T) << T^{1+epsilon} * exp(-c(sigma) * log T)

where c(sigma) > 0 for all sigma > 1/2. That is, the number of off-line zeros is at most T^{1-c+epsilon}, which is sublinear in T for any fixed sigma > 1/2.

This would be much stronger than the density hypothesis but weaker than RH. Proving this rigorously (even conditionally on GUE) seems very difficult because the GUE conjecture is a statement about local statistics, while zero density is a global count. The gap between local and global is precisely where the difficulty lies.

---

## 5. The Conrey-Snaith Ratios Recipe and Joint Distributions

### 5.1 The CFKRS Recipe

Conrey, Farmer, Keating, Rubinstein, and Snaith ("Integral Moments of L-functions," Proc. London Math. Soc. 91(1), 33-104, 2005; arXiv: math/0206018) developed a systematic "recipe" for conjecturing the full asymptotic main terms of moments of L-functions. The recipe instructs one to:

1. Write the moment integral and replace zeta by its approximate functional equation
2. Identify the "diagonal" and "off-diagonal" contributions
3. Keep terms where the product of gamma factors is not oscillating rapidly
4. Factor the result into an arithmetic part (Euler product) and a combinatorial/RMT part

The result is a precise conjecture for all main terms (not just leading order) in the 2k-th moment. Conrey and Snaith (Proc. London Math. Soc. 94(3), 594-646, 2007; arXiv: math/0509480) further developed the **L-functions ratios conjecture** (building on Conrey-Farmer-Zirnbauer), which predicts averages of ratios of L-functions:

    < prod_{i} L(s_i, chi) / prod_{j} L(w_j, chi) >_{chi in Family}

This is incredibly powerful because it encodes, as special cases, moments, correlation functions of zeros, one-level density, and more.

### 5.2 Predicting Joint Distributions at Multiple sigma Values

The ratios conjecture naturally extends to predict the joint distribution of zeta(sigma_1 + it), zeta(sigma_2 + it), ..., zeta(sigma_k + it) for multiple values of sigma. Specifically, consider the joint moment:

    (1/T) int_0^T |zeta(sigma_1 + it)|^{2k_1} * |zeta(sigma_2 + it)|^{2k_2} dt

The CFKRS recipe can be applied to this integral by writing each zeta factor using its approximate functional equation and applying the standard recipe. The prediction should involve:

- A product of RMT factors (from the CUE model at each sigma value)
- An arithmetic factor (involving the primes, but now coupling the different sigma values)
- A **cross-correlation term** that encodes how the values at different sigma values are related

The cross-correlation term is the key new ingredient. In the RMT model, the values of a CUE characteristic polynomial Z_N(z) at z = e^{i*theta} and z = r * e^{i*theta} (with r < 1) are correlated: the characteristic polynomial at a point inside the unit circle is determined by the eigenvalues, and its absolute value at distance delta from the circle is approximately the product of distances from each eigenvalue, which is controlled by the nearest eigenvalue.

For zero density applications, the relevant quantity is:

    P(|zeta(sigma + it)| < epsilon for some sigma in [1/2 + delta, 1])

A zero at rho = beta + i*gamma with beta > 1/2 forces |zeta(beta + i*gamma)| = 0 and also constrains |zeta(1/2 + i*gamma)| via the functional equation and the hybrid model. The CFKRS joint distribution prediction would constrain how often this simultaneous condition can occur.

### 5.3 Potential Improvement Mechanism

The standard zero density argument bounds N(sigma, T) using only information about zeta at a single sigma value (the zero-detecting polynomial is evaluated at sigma). But the ratios conjecture gives joint information across multiple sigma values. One could potentially improve density estimates by using a **multi-point zero-detecting scheme**:

Construct a zero-detecting polynomial D(s) that is sensitive to zeros at sigma = sigma_0 but uses cancellation information from the behavior of zeta at sigma = 1/2 (where we have much better understanding via RMT) and sigma = 1 (where we have good pointwise bounds). The CFKRS recipe would predict the relevant mean values for such multi-point polynomials, and if the cross-correlations are favorable, this could sharpen the density bound.

This is analogous to how Levinson's method for critical-line zeros uses a mollifier that involves information from the 1-line to detect zeros on the critical line. The proposal here is to use the full sigma-profile, predicted by CFKRS, in the density estimate context.

---

## 6. Synthesis: A Roadmap for Hybrid Improvements

### 6.1 What We Have

Collecting the threads above, there are several distinct mechanisms by which RMT information could constrain zero density:

| Mechanism | Source | Type of Bound | Regime |
|-----------|--------|---------------|--------|
| Hybrid dichotomy | GHK formula | Case split: Euler vs Hadamard | All sigma |
| Selberg CLT + Jensen | Value distribution | Exponential decay in delta | sigma near 1/2 |
| FHK maximum + three-lines | Extreme values | Local density per interval | sigma not too close to 1/2 |
| GUE repulsion | Pair correlation | Quadratic repulsion ~delta^2 | sigma very near 1/2 |
| CFKRS joint distribution | Ratios conjecture | Cross-sigma constraints | All sigma |
| Freezing transition moments | FHK | Improved moment bounds for k > 1 | Affects Halasz-Montgomery |

### 6.2 The Most Promising Direction

In the author's assessment, the most promising short-term direction is **Mechanism 6: using the freezing transition to improve moment bounds that feed into Halasz-Montgomery**. Here is why:

The Guth-Maynard proof of the 30/13 density exponent ultimately relies on bounding the number of points where a Dirichlet polynomial takes large values. The key input is a **large value estimate**: if D(t) = sum a_n n^{-it} with |a_n| <= 1, and |D(t)| > N^{sigma} on a 1-separated set W in [0, T], then |W| <= T^{A(sigma)} * N^{B(sigma)} for certain exponents A, B.

The proof uses the singular value decomposition of a matrix whose entries involve the Dirichlet coefficients evaluated at the points of W. The RMT connection: this matrix is structurally similar to a random matrix (its entries are oscillatory sums), and the singular value distribution should be constrained by RMT predictions.

Concretely, if the FHK freezing transition implies that the 2k-th moment of the Dirichlet polynomial is bounded by (log N)^{min(k^2, 2k-1)} (the freezing cap at k = 1), then the contribution of extreme values to the large value count is reduced. For k > 1 (the frozen regime), the effective power saving is:

    k^2 - (2k - 1) = (k - 1)^2

which is quadratic in k. This means that the extreme-value tail contributes less than the naive Keating-Snaith extrapolation would predict, potentially allowing improved large value estimates for the relevant range of sigma.

### 6.3 Required Ingredients for a Proof

To make any of these approaches rigorous, one would need:

1. **Unconditional moment bounds capturing the freezing transition.** The best current unconditional bounds are Harper's theorem (2013) giving sharp bounds on the maximum of |zeta| on short intervals. The FHK conjecture is now a theorem (Arguin-Bourgade-Radziwill), but translating this into moment bounds of the specific form needed for Halasz-Montgomery requires additional work.

2. **A rigorous hybrid formula usable in zero detection.** Currently, the GHK hybrid formula has been made rigorous only for averaged quantities (moments). A pointwise version with controlled error would be needed for the dichotomy argument.

3. **RMT predictions for structured matrices.** The matrices arising in Guth-Maynard's proof have arithmetic structure (entries involve n^{it} for integers n). Proving that their singular value distribution matches RMT predictions would require new universality results.

---

## 7. Proposed Numerical Experiments

Using the existing computational infrastructure in `/home/user/riemann/`, we propose the following experiments to test whether RMT predictions constrain off-line zero density. The codebase provides: Riemann-Siegel evaluation (`riemann_siegel.py`), zero-finding (`zero_finding.py`, `zero_counting.py`), pair correlation statistics (`pair_correlation.py`), Keating-Snaith moments (`moments.py`), and zero density exponent comparisons (`zero_density.py`).

### Experiment 1: Hybrid Formula Verification

**Goal:** Numerically verify the GHK hybrid splitting and test its accuracy as a function of the parameter X.

**Method:**
1. Using `riemann_siegel.py::hardy_z(t)` and `zero_finding.py::find_zeros_in_range()`, compute zeta(1/2+it) and locate all zeros in [T, T+H] for T up to 10^5 and H = 100.
2. For each t in a fine grid, compute the partial Euler product P_X(1/2+it) = prod_{p <= X} (1 - p^{-(1/2+it)})^{-1} for several values of X (X = 10, 50, 100, T^{1/4}, T^{1/2}).
3. Compute Z_X(1/2+it) = zeta(1/2+it) / P_X(1/2+it).
4. Verify that the moments of Z_X match CUE predictions (using `moments.py::rmt_factor(k)`) and that the moments of P_X match the arithmetic factor (`moments.py::arithmetic_factor(k)`).
5. Test how the error in the hybrid formula depends on X and on the distance to the nearest zero.

**Implementation Note:** The function `moments.py::numerical_moment(k, T)` already computes zeta moments; we would add a `hybrid_moment(k, T, X)` function that separately computes the P_X and Z_X moments.

### Experiment 2: Off-Line Value Distribution from Three-Lines Interpolation

**Goal:** Test whether the known distribution of |zeta(1/2+it)| (Selberg CLT) correctly predicts the distribution of |zeta(sigma+it)| for sigma > 1/2 via three-lines interpolation.

**Method:**
1. For T in {1000, 5000, 10000, 50000}, compute |zeta(1/2+it)| and |zeta(sigma+it)| at a fine grid of t-values using `riemann_siegel.py`. (Note: evaluating zeta at sigma > 1/2 requires a modified Riemann-Siegel routine or direct Euler-Maclaurin summation.)
2. For sigma = 0.55, 0.60, 0.65, 0.70, 0.75, compute the empirical distribution of log|zeta(sigma+it)|.
3. Compare with: (a) the three-lines interpolation from the empirical sigma=1/2 distribution and the known sigma=1 distribution, and (b) the CFKRS joint distribution prediction.
4. Measure the discrepancy and assess whether the interpolation gives useful bounds.

### Experiment 3: Local Zero Density from FHK Maximum Bounds

**Goal:** Test the FHK-derived local zero density prediction against the known distribution of zeros.

**Method:**
1. Using `zero_finding.py::find_zeros_in_range()`, compute all zeros in [T, T+1] for many values of T (e.g., T = 100, 200, ..., 10000).
2. For each interval, compute max_{t in [T,T+1]} |Z(t)| using `riemann_siegel.py::hardy_z()`.
3. Verify the FHK prediction: max |Z| ~ log T / (log log T)^{3/4}.
4. Apply Jensen's formula numerically: for each unit interval, compute the Jensen integral over circles of various radii and compare with the zero count from step 1.
5. Vary the radius to probe the off-line direction: use circles extending to sigma = 1/2 + delta for delta = 0.01, 0.05, 0.1, and compute the Jensen-implied off-line zero bound.

### Experiment 4: Pair Correlation and Horizontal Distribution

**Goal:** Test whether the empirical pair correlation of zeros, as computed by `pair_correlation.py`, is consistent with the RMT-predicted horizontal zero repulsion.

**Method:**
1. Compute pair correlations using `pair_correlation.py::pair_correlation_empirical(zeros)` for zeros in [T, 2T].
2. From the pair correlation data, extract the small-spacing behavior: how does R_2(alpha) behave as alpha -> 0? Compare with the GUE prediction R_2(alpha) ~ alpha^2 * pi^2 / 3 (quadratic vanishing).
3. Using the connection between vertical repulsion and horizontal repulsion (via the GHK hybrid model), predict the off-line zero density and compare with `zero_density.py::estimate_N_sigma()` for various sigma values.
4. Test whether the empirical number variance (`pair_correlation.py::number_variance_empirical()`) is consistent with GUE rigidity, which would constrain how many zeros can "escape" to sigma > 1/2.

### Experiment 5: Moment Freezing and Large Value Statistics

**Goal:** Detect the signature of the FHK freezing transition in the moments of zeta and assess its implications for large value counts.

**Method:**
1. Using `moments.py::numerical_moment(k, T)`, compute the 2k-th moment for k = 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0 and T = 500, 1000, 5000.
2. Fit the growth rate: plot log(M_k(T)) vs log(log T) and extract the effective exponent alpha(k). Compare with the Keating-Snaith prediction alpha(k) = k^2 (unfrozen) and the FHK prediction alpha(k) = 2k-1 (frozen, for k > 1).
3. Compute the empirical large value statistics: for each T, count how often |Z(t)| > V for V = (log T)^a with a = 0.5, 1.0, 1.5, 2.0. Compare the empirical tail with the Gaussian (Selberg) and log-correlated (FHK) predictions.
4. Assess whether the freezing transition is visible at computationally accessible heights (T ~ 10^4 to 10^5) or whether much larger T is needed.

---

## 8. Conclusions and Open Problems

The interface between random matrix theory and zero density estimates is rich with unexploited connections. The key findings of this survey are:

1. **The GHK hybrid formula provides a natural framework** for splitting the zero-detection problem into an arithmetic part (Euler product, controlled by Selberg's CLT) and a spectral part (Hadamard product, modeled by CUE). Making this dichotomy rigorous in the zero-detection context would be a significant advance.

2. **The FHK conjecture (now theorem) has immediate consequences** for local zero density through the Jensen/three-lines mechanism. These consequences are strongest for sigma very close to 1/2, complementing classical density estimates that are strongest for sigma near 1.

3. **The freezing transition suggests improved moment bounds** for k > 1, which could sharpen the large value estimates feeding into Halasz-Montgomery. This is arguably the most promising route to an unconditional improvement.

4. **The CFKRS ratios conjecture encodes cross-sigma information** that is not used in classical density arguments. Developing a multi-point zero-detection scheme that exploits this joint distribution could yield new density bounds.

5. **RMT-predicted zero repulsion gives a heuristic density exponent** that is much stronger than the density hypothesis for sigma near 1/2 but becomes trivial for sigma near 1. The complementarity with classical estimates suggests that a hybrid approach could improve the overall bound.

### Key Open Problems

- **Problem 1:** Make the GHK hybrid dichotomy rigorous for zero detection. Even a conditional result (assuming GUE) would be very interesting.
- **Problem 2:** Prove that the FHK freezing transition implies improved large value estimates for Dirichlet polynomials.
- **Problem 3:** Use the CFKRS ratios conjecture to develop multi-sigma zero-detecting polynomials.
- **Problem 4:** Establish universality results for the singular value distribution of the matrices arising in Guth-Maynard's method, connecting to RMT.
- **Problem 5:** Determine whether the RMT-predicted density exponent can be rigorously established even in a restricted sigma-range.

---

## References

1. Gonek, S.M., Hughes, C.P., and Keating, J.P. "A hybrid Euler-Hadamard product for the Riemann zeta function." Duke Math. J. 136(3), 507-549 (2007). [arXiv:math/0511182](https://arxiv.org/abs/math/0511182)

2. Fyodorov, Y.V., Hiary, G.A., and Keating, J.P. "Freezing Transition, Characteristic Polynomials of Random Matrices, and the Riemann Zeta-Function." [arXiv:1202.4713](https://arxiv.org/abs/1202.4713)

3. Arguin, L.-P., Bourgade, P., and Radziwill, M. "The Fyodorov-Hiary-Keating Conjecture. I." [arXiv:2007.00988](https://arxiv.org/abs/2007.00988)

4. Arguin, L.-P., Bourgade, P., and Radziwill, M. "The Fyodorov-Hiary-Keating Conjecture. II." [arXiv:2307.00982](https://arxiv.org/abs/2307.00982)

5. Conrey, J.B., Farmer, D.W., Keating, J.P., Rubinstein, M.O., and Snaith, N.C. "Integral moments of L-functions." Proc. London Math. Soc. 91(1), 33-104 (2005). [arXiv:math/0206018](https://arxiv.org/abs/math/0206018)

6. Conrey, J.B. and Snaith, N.C. "Applications of the L-functions ratios conjectures." Proc. London Math. Soc. 94(3), 594-646 (2007). [arXiv:math/0509480](https://arxiv.org/abs/math/0509480)

7. Keating, J.P. and Snaith, N.C. "Random matrix theory and zeta(1/2+it)." Comm. Math. Phys. 214(1), 57-89 (2000).

8. Guth, L. and Maynard, J. "New large value estimates for Dirichlet polynomials." [arXiv:2405.20552](https://arxiv.org/abs/2405.20552) (2024).

9. Baluyot, S., Goldston, D., Suriajaya, A., and Turnage-Butterbaugh, C. "Pair correlation and the distribution of zeros." [arXiv:2503.15449](https://arxiv.org/abs/2503.15449) (2025).

10. Radziwill, M. and Soundararajan, K. "Selberg's central limit theorem for log|zeta(1/2+it)|." L'Enseignement Mathematique (2017). [arXiv:1509.06827](https://arxiv.org/abs/1509.06827)

11. Mezzadri, F. "Random matrix theory and the zeros of zeta'(s)." [arXiv:math-ph/0207044](https://arxiv.org/abs/math-ph/0207044) (2002).

12. Fazzari, A. et al. "Selberg's Central Limit Theorem weighted by Linear Statistics of Zeta Zeros." [arXiv:2507.04150](https://arxiv.org/abs/2507.04150) (2025).

13. Fyodorov, Y.V. and Keating, J.P. "Freezing transitions and extreme values: random matrix theory, and disordered landscapes." Phil. Trans. Roy. Soc. A 372(2007), 20120503 (2014).

14. Tao, T. "254A, Notes 6: Large values of Dirichlet polynomials, zero density estimates, and primes in short intervals." [Blog post](https://terrytao.wordpress.com/2015/02/13/254a-notes-6-large-values-of-dirichlet-polynomials-zero-density-estimates-and-primes-in-short-intervals/) (2015).

15. Arguin, L.-P. and Hamdan, R. "The Fyodorov-Hiary-Keating Conjecture on Mesoscopic Intervals." [arXiv:2405.06474](https://arxiv.org/abs/2405.06474) (2024).

16. Rodgers, B. "Arithmetic Consequences of the GUE Conjecture for Zeta Zeros." Michigan Math. J. 74(5) (2024).

17. "50 Years of Number Theory and Random Matrix Theory Conference." [IAS](https://www.ias.edu/math/events/50yntrmt).
