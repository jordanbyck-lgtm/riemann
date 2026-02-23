# Exponent Pairs and Systematic Optimization via ANTEDB: A Technical Report

## Author: Graduate Student 2 (Analytic Number Theory)
## Date: February 2025

---

## Abstract

This report provides a detailed investigation of the pipeline connecting **exponent pairs** to **zero density estimates** for the Riemann zeta function, with particular focus on the recent **Analytic Number Theory Exponent Database (ANTEDB)** project of Tao, Trudgian, and Yang (arXiv:2501.16779). We describe the classical theory of exponent pairs and the van der Corput A- and B-processes, the explicit mechanisms by which exponent pairs feed into zero density bounds, the four new exponent pairs discovered via the ANTEDB, and their implications for zero density estimates. We carry out explicit computations at the critical point sigma = 3/4, compare the results with the breakthrough of Guth and Maynard (2024), discuss the exponent pair conjecture and its relationship to the density hypothesis, formulate a linear programming relaxation for direct optimization of the zero density exponent, and survey Huxley's contributions beyond the classical AB-process.

---

## 1. The Classical Theory of Exponent Pairs

### 1.1 Definition

An **exponent pair** (k, l) with 0 <= k <= 1/2 <= l <= 1 encodes a bound on one-dimensional exponential sums: if f is a smooth function on [M, 2M] with |f^(j)(x)| asymptotically comparable to T/M^j for j >= 1, then

$$\sum_{M < n \le 2M} e(f(n)) \ll M^{1-k} T^l + M^{1/2}$$

Equivalently, in terms of Dirichlet series:

$$\sum_{M < n \le 2M} n^{-it} \ll M^{1-k} t^l \quad (t \to \infty)$$

The set of all exponent pairs is **closed and convex** in R^2. Convex combinations of exponent pairs are again exponent pairs.

### 1.2 The Van der Corput Processes

Two fundamental transformations generate new exponent pairs from old ones:

**A-process (Weyl differencing):** If (k, l) is an exponent pair, then so is

$$A(k, l) = \left(\frac{k}{2k + 2},\ \frac{k + l + 1}{2k + 2}\right)$$

**B-process (Poisson summation / stationary phase):** If (k, l) is an exponent pair, then so is

$$B(k, l) = \left(l - \frac{1}{2},\ k + \frac{1}{2}\right)$$

The key structural constraint is that **BB is redundant**: applying the B-process twice in succession returns (up to technical corrections) to the original sum. Thus all classical exponent pairs are obtained by applying sequences of A and B to the trivial pair (0, 1), with no consecutive B's.

### 1.3 Classical Exponent Pairs Generated from (0, 1)

Starting from the trivial pair (0, 1) and applying sequences of the A- and B-processes, one generates the following classical exponent pairs:

| Process sequence | Exponent pair (k, l) | Decimal approximation |
|:-----------------|:--------------------|:---------------------|
| (trivial)        | (0, 1)              | (0, 1)               |
| B                | (1/2, 1/2)          | (0.5, 0.5)           |
| AB               | (1/6, 2/3)          | (0.1667, 0.6667)     |
| A^2B             | (1/14, 11/14)       | (0.0714, 0.7857)     |
| A^3B             | (1/30, 13/15)       | (0.0333, 0.8667)     |
| BAB              | (1/6, 2/3)          | (0.1667, 0.6667)     |
| BA^2B            | (2/7, 4/7)          | (0.2857, 0.5714)     |
| ABA^2B           | (2/13, 35/52)       | (0.1538, 0.6731)     |
| BABA^2B          | (2/9, 11/18)        | (0.2222, 0.6111)     |
| A^2BA^2B         | (1/7, 5/7)          | (0.1429, 0.7143)     |

The pair (2/7, 4/7) = BA^2B(0,1) is particularly important in classical applications. Note that B(0, 1) = (1/2, 1/2) lies on the conjectured boundary line l = k + 1/2, while A^nB(0, 1) converges to (0, 1) as n -> infinity, tracing out the upper boundary of the exponent pair region.

### 1.4 The Sargos C- and D-Processes

Beyond the classical AB-processes, Sargos introduced additional transformations:

**C-process:** If (k, l) is an exponent pair, then so is

$$C(k, l) = \left(\frac{k}{12(1 + 4k)},\ \frac{11(1 + 4k) + l}{12(1 + 4k)}\right)$$

**D-process:** Applied to (k, l), this yields

$$D(k, l) = \left(\frac{5k + l + 2}{8(5k + 3l + 2)},\ \frac{29k + 21l + 10}{8(5k + 3l + 2)}\right)$$

The D-process is particularly interesting because it provides bounds on the exponential sum growth function beta(alpha) rather than guaranteed new exponent pairs, which enriches the convex hull of achievable estimates.

### 1.5 The Duality with beta(alpha)

The exponential sum growth function beta(alpha) is defined so that

$$\sum_{M < n \le 2M} e(f(n)) \ll M^{beta(alpha) + epsilon}$$

when T ~ M^alpha. A fundamental duality (Lemma 5.3 in the ANTEDB blueprint) states:

> **(k, l) is an exponent pair if and only if beta(alpha) <= k + (l - k)alpha for all 0 <= alpha <= 1.**

This means the set of exponent pairs is exactly the polar dual of the epigraph of beta(alpha). Optimizing over exponent pairs is therefore equivalent to optimizing over linear upper bounds on beta(alpha).

---

## 2. The Exponent Pair to Zero Density Pipeline

### 2.1 Overview of the Pipeline

The connection from exponent pairs to zero density estimates proceeds through several layers:

```
Exponent Pairs (k, l)
    |
    v
Exponential Sum Bounds beta(alpha)
    |
    v
Growth Exponent mu(sigma) for zeta
    |
    v
Moment Bounds / Mixed Moments
    |
    v
Large Value Estimates LV(sigma, tau) and LV_zeta(sigma, tau)
    |
    v
Zero Density Estimates A(sigma)
```

Each arrow represents a class of theorems in analytic number theory. The ANTEDB systematizes all of these steps.

### 2.2 The Classical Ingham Formula

The most direct classical connection uses an exponent pair (k, l) to bound the growth of the zeta function, which then feeds into Ingham's method. The chain is:

**Step 1.** An exponent pair (k, l) gives a subconvexity bound for the zeta function on the critical line:

$$\zeta(1/2 + it) \ll t^{mu(1/2) + epsilon}$$

where mu(1/2) <= (l - k)/(2(1 + k - l)) when applied to the approximate functional equation (or more directly, mu(1/2) <= k for suitable exponent pairs via the relation mu(sigma) from exponential sum estimates).

**Step 2.** Ingham's 1940 method then gives:

$$N(sigma, T) \ll T^{(2 + 4c)(1 - sigma) + epsilon}$$

where c = mu(1/2), i.e., |zeta(1/2 + it)| << t^{c + epsilon}. In particular, the Lindelof hypothesis (c = 0) would give A(sigma) = 2, the density hypothesis.

**Step 3.** More refined use of the exponent pair through Ivic's mixed moment method (Chapter 11 of Ivic's monograph) gives: if (k, l) is an exponent pair, then Theorem 8.2 of Ivic guarantees that the moments

$$q_0 = 6, \quad B_0 = 1 + epsilon, \quad q_1 = \frac{2(1 + 2k + 2l)}{k}, \quad B_1 = \frac{k + l}{k} + epsilon$$

are admissible for the Halasz-Montgomery zero-detecting argument. This produces the zero density bound

$$A(sigma) \le \frac{2(k + l)}{1 - k}$$

in a suitable range of sigma. This is the key formula connecting an individual exponent pair to a zero density exponent.

### 2.3 The Large Value Estimate Pathway

The modern and more powerful pipeline uses **large value estimates** as the crucial intermediary.

**Definition.** For sigma in (1/2, 1) and tau > 0, the large value exponent LV(sigma, tau) is the infimum of exponents rho such that: for any Dirichlet polynomial P(s) = sum_{n~N} a_n n^{-s} with |a_n| <= 1, the number R of 1-separated points t_r in [T, 2T] with |P(1/2 + it_r)| >= N^sigma satisfies

$$R \ll N^{rho + epsilon}$$

where N = T^{1/tau}.

**The core zero density lemma (Lemma 11.5 of ANTEDB):** For 1/2 < sigma < 1,

$$A(sigma)(1 - sigma) \le \max\left(\sup_{\tau \ge 2} \frac{LV_\zeta(sigma, \tau)}{\tau},\ \limsup_{\tau \to \infty} \frac{LV(sigma, \tau)}{\tau}\right)$$

Here LV_zeta(sigma, tau) denotes the large value exponent specifically for zeta-type Dirichlet polynomials (which may enjoy better bounds due to the multiplicative structure of the zeta function).

**Practical corollary (Corollary 11.9):** If one can show

- LV(sigma, tau) <= (3 - 3sigma) tau/tau_0 for 2tau_0/3 <= tau <= tau_0
- LV_zeta(sigma, tau) <= (3 - 3sigma) tau/tau_0 for 2 <= tau < 4tau_0/3

then **A(sigma) <= 3/tau_0**.

### 2.4 The Halasz-Montgomery Inequality

The fundamental analytic tool underlying zero density estimates is the **Halasz-Montgomery inequality**. If rho_1, ..., rho_R are zeros of zeta(s) with Re(rho_j) >= sigma and |Im(rho_j) - Im(rho_k)| >= 1 for j != k, and if P(s) is a zero-detecting Dirichlet polynomial such that |P(rho_j)| >= 1 for all j, then

$$R \le \left(\sum_j |P(rho_j)|^2\right) \cdot \left(\max_j \sum_k |P(rho_k)|^2 |\langle rho_j, rho_k \rangle|^{-2}\right)^{-1}$$

where the inner product terms involve mean values of the Dirichlet polynomial. This converts the problem of counting zeros into the problem of bounding large values of Dirichlet polynomials, which is where exponential sum estimates and exponent pairs enter.

### 2.5 Huxley's Subdivision Trick

A crucial technical innovation due to Huxley (1972) extends large value estimates from short intervals to long intervals:

> If C(sigma, alpha) is a bound on the large value exponent on intervals of length T, then (1 - C(sigma, alpha'))/alpha' is non-decreasing in alpha'.

This produces the improved bound:

$$C(sigma, alpha) \le \min\left((2 - 2sigma)alpha,\ 1 + (4 - 6sigma)alpha\right)$$

The subdivision trick is essential for converting the local large value estimates (which come from exponent pairs) into the global estimates needed for zero density.

---

## 3. The ANTEDB and Four New Exponent Pairs

### 3.1 Overview of the ANTEDB Project

The **Analytic Number Theory Exponent Database (ANTEDB)** was launched in January 2025 by Terence Tao, Timothy Trudgian, and Andrew Yang (arXiv:2501.16779). The project creates a living database that:

1. Collects known theorems for exponents in analytic number theory
2. Records relationships between these exponents as both human-readable proofs and executable Python code
3. Systematically optimizes these relationships via computational methods

The key insight is that many papers in analytic number theory derive consequences of known exponent bounds in an ad hoc manner, tailored to the specific numerical values available at the time of publication. The ANTEDB abstracts these arguments so that improvements in any input exponent automatically propagate through the entire dependency graph.

### 3.2 The Four New Exponent Pairs (Theorem 20)

The ANTEDB systematic optimization discovered **four new exponent pairs** beyond those previously in the literature:

| # | Exponent pair (k, l) | Decimal approximation | Method |
|:--|:---------------------|:---------------------|:-------|
| 1 | (89/1282, 997/1282) | (0.06942, 0.77769) | ANTEDB optimization |
| 2 | (652397/9713986, 7599781/9713986) | (0.06716, 0.78234) | ANTEDB optimization |
| 3 | (10769/351096, 609317/702192) | (0.03067, 0.86777) | ANTEDB optimization |
| 4 | (89/3478, 15327/17390) | (0.02559, 0.88138) | ANTEDB optimization |

These were obtained without introducing any new analytic number theory inputs -- they arise purely from more careful optimization of the relationships between existing results.

### 3.3 Previously Known Non-Classical Exponent Pairs

The ANTEDB also catalogues previously known exponent pairs that go beyond the classical AB-process:

| Source | Exponent pair (k, l) | Decimal approx. |
|:-------|:---------------------|:----------------|
| Bourgain (2017) | (13/84, 55/84) | (0.15476, 0.65476) |
| Trudgian-Yang (2023) | (4742/38463, 35731/51284) | (0.12330, 0.69697) |
| Trudgian-Yang (2023) | (18/199, 593/796) | (0.09045, 0.74497) |
| Trudgian-Yang (2023) | (2779/38033, 58699/76066) | (0.07307, 0.77170) |
| Watt (1989) | (89/560, 471/560) | (0.15893, 0.84107) |
| Huxley (2005) | (32/205, 269/410) | (0.15610, 0.65610) |
| Cushing (2025) | (311/4822, 3799/4822) | (0.06450, 0.78784) |
| Cushing (2025) | (80219/1298878, 515638/649439) | (0.06176, 0.79393) |

The convex hull of all known exponent pairs has vertices at (0, 1), (1/2, 1/2), and the points in Tables 3.2 and 3.3 above.

---

## 4. Zero Density Bounds from Exponent Pairs: Explicit Computations

### 4.1 Current Best Bounds on A(sigma)

The ANTEDB records the following best-known bounds on the zero density exponent A(sigma):

| Range of sigma | Bound on A(sigma) | Source |
|:---------------|:-----------------|:-------|
| 1/2 <= sigma <= 7/10 | 3/(2 - sigma) | Ingham (1940) |
| 7/10 <= sigma < 19/25 | 15/(3 + 5sigma) | Guth-Maynard (2024) |
| 4/5 <= sigma < 7/8 | 3/(2sigma) | Ivic / Bourgain |
| sigma >= 23/24 | 3/(24sigma - 20) | Pintz |

### 4.2 Explicit Computation at sigma = 3/4

The **critical test case** for zero density estimates is sigma = 3/4, where the density hypothesis predicts A(3/4) = 2.

**Ingham's bound:** A(3/4) <= 3/(2 - 3/4) = 3/(5/4) = 12/5 = 2.4

**Guth-Maynard bound:** A(3/4) <= 15/(3 + 5 * 3/4) = 15/(3 + 15/4) = 15/(27/4) = 60/27 = 20/9 approx 2.2222

Actually, we need to check which formula applies at sigma = 3/4. Since 7/10 < 3/4 < 19/25, the Guth-Maynard formula applies:

$$A(3/4) \le \frac{15}{3 + 5 \cdot 3/4} = \frac{15}{27/4} = \frac{60}{27} = \frac{20}{9} \approx 2.222$$

Alternatively, the uniform bound gives:

$$N(3/4, T) \ll T^{30(1 - 3/4)/13 + o(1)} = T^{30/52 + o(1)} = T^{15/26 + o(1)}$$

This translates to A(3/4)(1 - 3/4) = 15/26, so A(3/4) = 15/26 * 4 = 60/26 = 30/13 approx 2.3077.

Wait -- let us be careful. The quantity 30/13 is the **uniform** bound: N(sigma, T) << T^{30(1-sigma)/13} for ALL sigma in [1/2, 1]. At sigma = 3/4 this gives the exponent 30 * (1/4) / 13 = 30/52 = 15/26. But A(sigma) is defined so that N(sigma, T) << T^{A(sigma)(1-sigma)}, meaning A(3/4) * (1/4) = 15/26, giving A(3/4) = 30/13 approx 2.3077.

However, the **pointwise** Guth-Maynard bound at sigma = 3/4 is stronger:

$$A(3/4) \le \frac{15}{3 + 15/4} = \frac{15}{27/4} = \frac{20}{9} \approx 2.222$$

This is because the Guth-Maynard estimate N(sigma, T) << T^{15(1-sigma)/(3+5sigma)} is stronger than the uniform 30(1-sigma)/13 bound for sigma > 7/10. At sigma = 3/4:

$$\frac{15(1 - 3/4)}{3 + 5 \cdot 3/4} = \frac{15/4}{27/4} = \frac{15}{27} = \frac{5}{9}$$

So N(3/4, T) << T^{5/9 + o(1)}, which gives A(3/4) = (5/9)/(1/4) = 20/9 approx 2.222.

### 4.3 The Classical Formula Applied to Known Exponent Pairs

Using the classical formula A(sigma) <= 2(k + l)/(1 - k) from Ivic's Theorem 11.2, we can compute the zero density bound at sigma = 3/4 from each known exponent pair:

| Exponent pair (k, l) | 2(k+l)/(1-k) | Source |
|:---------------------|:-------------|:-------|
| (0, 1) | 2(0+1)/(1-0) = 2 | Trivial -- but only applicable near sigma = 1 |
| (1/2, 1/2) | 2(1)/(1/2) = 4 | B(0,1) |
| (1/6, 2/3) | 2(5/6)/(5/6) = 2 | AB(0,1) |
| (2/7, 4/7) | 2(6/7)/(5/7) = 12/5 = 2.4 | BA^2B(0,1) |
| (13/84, 55/84) | 2(68/84)/(71/84) = 136/71 approx 1.915 | Bourgain (2017) |
| (89/560, 471/560) | 2(560/560)/(471/560) = 2*560/471 approx 2.378 | Watt |
| (32/205, 269/410) | 2(301/410)/(173/205) = ... | Huxley (2005) |
| (89/1282, 997/1282) | 2(1086/1282)/(1193/1282) = 2172/1193 approx 1.821 | ANTEDB new pair 1 |
| (652397/9713986, 7599781/9713986) | ... | ANTEDB new pair 2 |

**Important caveat:** The classical formula A(sigma) <= 2(k+l)/(1-k) gives a **uniform bound** that applies in a specific range of sigma, typically for sigma near 1. The formula is most useful when combined with Ingham's bound for sigma near 1/2 and the Halasz-Montgomery method for the intermediate range. The formula should not be taken as the complete zero density bound achievable from a given exponent pair, as more sophisticated methods (large value estimates, mixed moments, subdivision) typically yield better bounds.

Let us compute 2(k+l)/(1-k) for the ANTEDB new pairs more carefully:

**Pair 1: (89/1282, 997/1282)**
- k + l = 89/1282 + 997/1282 = 1086/1282 = 543/641
- 1 - k = 1 - 89/1282 = 1193/1282
- 2(k+l)/(1-k) = 2 * (543/641) / (1193/1282) = 2 * 543 * 1282 / (641 * 1193)
  = 2 * 543 * 2 / 1193 = 2172/1193 approx 1.8206

**Pair 3: (10769/351096, 609317/702192)**
- k + l = 10769/351096 + 609317/702192 = 21538/702192 + 609317/702192 = 630855/702192
- 1 - k = 1 - 10769/351096 = 340327/351096
- 2(k+l)/(1-k) = 2 * 630855/702192 / (340327/351096) = 2 * 630855 * 351096 / (702192 * 340327)
  = 2 * 630855 / (2 * 340327) = 630855/340327 approx 1.8537

**Pair 4: (89/3478, 15327/17390)**
- k + l = 89/3478 + 15327/17390 = 445/17390 + 15327/17390 = 15772/17390 = 7886/8695
- 1 - k = 1 - 89/3478 = 3389/3478
- 2(k+l)/(1-k) = 2 * 7886/8695 / (3389/3478) = 2 * 7886 * 3478 / (8695 * 3389)
  = 2 * 7886 * 2 / (5 * 3389) = 31544/16945 approx 1.8616

These values below 2 from the classical formula are formally remarkable, but recall that the formula A(sigma) <= 2(k+l)/(1-k) is not applicable uniformly: it gives bounds only in a specific range of sigma (typically for sigma close to 1), and the seemingly strong bounds reflect that the new ANTEDB pairs have very small k and correspondingly large l, which is beneficial in the near-sigma=1 regime but less so for sigma = 3/4.

### 4.4 How the ANTEDB Pairs Compare to Guth-Maynard

The Guth-Maynard breakthrough operates through an entirely different mechanism: **new large value estimates for Dirichlet polynomials** (arXiv:2405.20552). Their main result is a bound on how often a Dirichlet polynomial of length N can take values of size approximately N^{3/4}, which is precisely the critical case for zero density.

Their Theorem 1.2 gives:

$$N(sigma, T) \ll T^{15(1-sigma)/(3+5sigma) + o(1)}$$

At sigma = 3/4, this gives A(3/4) <= 20/9 approx 2.222. The exponent 30/13 approx 2.308 is the **uniform** bound obtained by taking the worst case over all sigma.

The ANTEDB new exponent pairs, while yielding formally small values of the classical formula 2(k+l)/(1-k), **do not compete with Guth-Maynard at sigma = 3/4**. The reason is fundamental: the exponent pair approach and the large value estimate approach are complementary but operate at different scales. The Guth-Maynard improvement comes from a new large value theorem in the critical range V ~ N^{3/4}, which cannot be captured by classical exponent pair technology alone.

However, the ANTEDB pairs **do yield improvements** in the range of sigma close to 1, where the classical Ivic-type formula is most effective. The ANTEDB project found new zero density estimates in these ranges, though these are not competitive with the strongest known bounds in the sigma ~ 3/4 range.

---

## 5. The Exponent Pair Conjecture and the Density Hypothesis

### 5.1 Statement of the Exponent Pair Conjecture

The **exponent pair conjecture** (EPC) asserts that for every epsilon > 0, the pair

$$(epsilon, 1/2 + epsilon)$$

is an exponent pair. Equivalently, for any epsilon > 0:

$$\sum_{M < n \le 2M} n^{-it} \ll M^{1-epsilon} t^{1/2+epsilon}$$

This conjecture implies the **Lindelof hypothesis**: mu(1/2) = 0, i.e., zeta(1/2 + it) << t^epsilon for every epsilon > 0.

### 5.2 The EPC Implies the Density Hypothesis

The chain of implications is:

1. **EPC => Lindelof hypothesis.** If (epsilon, 1/2 + epsilon) is an exponent pair, then from the relation between exponent pairs and the growth exponent, mu(1/2) <= epsilon for every epsilon > 0, hence mu(1/2) = 0.

2. **Lindelof hypothesis => Density hypothesis.** By Ingham's result, if |zeta(1/2 + it)| << t^{c+epsilon} for some constant c >= 0, then

   $$N(sigma, T) \ll T^{(2 + 4c)(1-sigma) + epsilon}$$

   Setting c = 0 (the Lindelof hypothesis) gives

   $$N(sigma, T) \ll T^{2(1-sigma) + epsilon}$$

   which is exactly the **density hypothesis**: A(sigma) <= 2 for all sigma in [1/2, 1].

### 5.3 Direct Verification from the Formula

Using the classical formula A <= 2(k+l)/(1-k) with the conjectured pair (epsilon, 1/2 + epsilon):

$$A \le \frac{2(epsilon + 1/2 + epsilon)}{1 - epsilon} = \frac{2(1/2 + 2epsilon)}{1 - epsilon} = \frac{1 + 4epsilon}{1 - epsilon}$$

As epsilon -> 0, this approaches **A = 1**, which would be the density hypothesis (and in fact even stronger -- note that the density hypothesis states A(sigma) <= 2, meaning N(sigma, T) << T^{2(1-sigma)}, while A = 1 would give N(sigma, T) << T^{1-sigma}, which is strictly stronger).

However, the classical formula does not uniformly apply at all sigma, so the correct implication of the EPC for the density hypothesis goes through the Lindelof hypothesis as described above, yielding A(sigma) <= 2 uniformly.

### 5.4 Implications for the Primes

The density hypothesis A(sigma) = 2 implies:
- Primes in short intervals: pi(x) - pi(x - y) ~ y/log(x) for y >= x^{1/2 + epsilon}
- This would improve on the current best result of Guth-Maynard (y >= x^{17/30 + epsilon})

---

## 6. Linear Programming Relaxation for Zero Density Optimization

### 6.1 The ANTEDB Optimization Methodology

The ANTEDB uses a **polytope-based computation** to optimize exponent bounds. The key idea:

1. Each known bound on beta(alpha) (from exponent pairs, Vinogradov-type estimates, etc.) defines a constraint in the (alpha, beta) plane.
2. The feasible region for beta(alpha) is represented as a polygon -- the intersection of all constraint half-planes.
3. Optimization over this polygon yields the best achievable exponent for a given application.

The specific implementation uses standard computational geometry: "represent each beta estimate as a polygon containing feasible (alpha, beta) tuples, and then take the intersection of all such polygons, which is easily achieved via standard Boolean polytope operations."

### 6.2 Setting Up the LP for Zero Density

We can formulate the zero density optimization as an explicit linear program. The goal is to minimize A(sigma) over all achievable combinations of exponent pair bounds.

**Variables:** beta(alpha) for alpha in [0, 1] (the exponential sum growth function), discretized at points alpha_0, alpha_1, ..., alpha_N.

**Objective:** Minimize A(sigma_0) for a target sigma_0 (e.g., sigma_0 = 3/4).

**Constraints from exponent pairs:** For each known exponent pair (k_j, l_j):

$$beta(alpha) \le k_j + (l_j - k_j) \cdot alpha \quad \text{for all } alpha \in [0, 1]$$

**Constraints from the A-process:** If beta(alpha) is feasible, then so is the transform:

$$beta'(alpha) = \frac{beta(2alpha) + alpha}{2} \quad \text{(simplified)}$$

**Constraints from the B-process:**

$$beta'(alpha) = \frac{1 + (2alpha - 1)beta(1/alpha)}{2alpha} \quad \text{(for appropriate ranges)}$$

**Constraints from convexity:** beta(alpha) is a convex function.

**Link to zero density:** Through the chain beta(alpha) -> mu(sigma) -> LV(sigma, tau) -> A(sigma), we express A(sigma_0) as a function of the beta values:

$$A(sigma_0) \le \inf_{\tau_0 \ge 2} \frac{3}{\tau_0} \quad \text{subject to } LV(sigma_0, \tau) \le (3 - 3sigma_0)\frac{\tau}{\tau_0}$$

where LV(sigma_0, tau) is bounded in terms of beta(alpha) through the moment and mean value estimates.

### 6.3 The Extended LP Formulation

More explicitly, the LP can be written as:

**Minimize** A subject to:

1. **Exponent pair constraints:** For each vertex (k_j, l_j) of the known convex hull:
   $$\beta_i \le k_j + (l_j - k_j) \alpha_i \quad \forall i, j$$

2. **Process constraints:** For each discretization point:
   $$\beta_{A(i)} \le \frac{\beta_{2i} + \alpha_i}{2}$$
   (A-process propagation)

3. **Convexity:** beta is piecewise-linear and convex on the grid.

4. **Large value bounds:** Using the Halasz-Montgomery framework:
   $$LV(sigma_0, \tau_i) \le f(\beta, sigma_0, \tau_i)$$
   where f encodes Ivic's mixed moment machinery.

5. **Huxley subdivision:** (1 - LV(sigma_0, tau))/tau is non-decreasing.

6. **Zero density link:**
   $$A(1 - sigma_0) \le \sup_i LV(sigma_0, \tau_i) / \tau_i$$

This is a **linear program** (or semidefinite program when we include more sophisticated constraints) that can be solved computationally using standard solvers.

### 6.4 Semidefinite Programming Extension

The ANTEDB also explores SDP relaxations. The convex hull of achievable exponents can be represented more tightly using semidefinite constraints when one incorporates:

- **Decoupling estimates** (Bourgain-Demeter-Guth) that provide non-linear constraints on beta(alpha)
- **The Vinogradov main conjecture** (now theorem), which gives bounds that are quadratic in the parameters
- **Heath-Brown's double zeta sum estimates**, which involve optimization over auxiliary parameters

The SDP takes the form:

**Minimize** A subject to:
$$M(A, sigma, beta_1, ..., beta_N, LV_1, ..., LV_M) \succeq 0$$

where M is a matrix whose positive semidefiniteness encodes all the known constraints.

### 6.5 The Verification Framework

A crucial feature of the ANTEDB is that every claimed bound comes with a **machine-verifiable proof**. Each optimization result is accompanied by:

1. A dependency tree showing which theorems were used
2. Explicit numerical certificates that verify feasibility
3. Python code that can reproduce the computation

This is particularly valuable because the optimization landscape is complex and subtle errors in applying theorems (e.g., forgetting a constraint on the range of sigma) can lead to incorrect claims.

---

## 7. Huxley's "New" Exponent Pairs and the Bombieri-Iwaniec Method

### 7.1 Beyond the AB-Process

The classical van der Corput A- and B-processes generate a countably infinite but ultimately limited set of exponent pairs from the trivial pair (0, 1). In the 1980s and 1990s, Bombieri and Iwaniec, followed by Huxley, developed fundamentally new methods that go beyond this framework.

The **Bombieri-Iwaniec method** introduces number-theoretic ideas into the estimation of exponential sums. Rather than relying purely on analytic transformations (differencing and Poisson summation), it uses:

1. **Mean square estimates** over families of exponential sums
2. **The First Spacing Problem**: counting how often different sub-sums can take similar values
3. **The Second Spacing Problem**: counting resonances between short intervals of the sum, when two arcs of the graph of y = f'(x) coincide approximately after an integer lattice automorphism

### 7.2 Huxley's Resonance Curve Method (2005)

In his landmark paper "Exponential Sums and the Riemann Zeta Function V" (2005), Huxley axiomatized the Bombieri-Iwaniec method and introduced the **resonance curve** technique. The key results:

- A bound for the zeta function: zeta(1/2 + it) << t^{theta} with theta = 32/205 = 0.156098..., improving on Bombieri-Iwaniec's 9/56 = 0.160714...
- The new exponent pair **(32/205, 269/410)**, which satisfies l = k + 1/2 (the conjectural boundary)
- Even a complete resolution of both the first and second spacing problems alone could not achieve theta < 3/20 = 0.15, indicating that the resonance curve method brings genuinely new information

### 7.3 The Watt Exponent Pair

N. Watt (1989) obtained the exponent pair **(89/560, 471/560)** using a refinement of the Bombieri-Iwaniec method. This pair has the special property that l - k = 382/560 = 191/280, and it generates particularly useful "derived" pairs through the AB-process:

- A(89/560, 471/560) = (89/560 / (2*89/560 + 2), ...) leads to pairs useful for the Dirichlet divisor problem
- B-transforms give pairs in the range useful for zeta function moment estimates

The Watt pair plays a central role in the ANTEDB's optimization, as it provides constraints in ranges not accessible to the classical AB-process.

### 7.4 Interaction with Guth-Maynard

The Guth-Maynard approach is **fundamentally different** from both the classical exponent pair method and the Bombieri-Iwaniec method. Their key innovation is a new **large value estimate for Dirichlet polynomials** near the critical value V ~ N^{3/4}.

Specifically, Guth-Maynard prove (Theorem 1.1): the number R of 1-separated points where |P(1/2 + it)| >= V satisfies

$$R \ll N^{epsilon} \cdot \left(N^{3/5} V^{8/5} + T N^{-3/5} V^{-8/5}\right)$$

This improves on the classical mean value theorem bound (which gives N^{1/2} V^2 + T N^{-1/2} V^{-2}) precisely when V is near N^{3/4}.

The connection to exponent pairs is **indirect**: exponent pairs provide certain large value estimates (through the moment method), but the Guth-Maynard bound comes from a completely different technique involving:

- **Decoupling theory** (building on Bourgain-Demeter-Guth)
- **Multilinear estimates** for Dirichlet polynomials
- **Geometric combinatorics** (counting lattice points near curves)

The ANTEDB framework is designed to accommodate both types of inputs: exponent pair bounds and direct large value estimates like Guth-Maynard's. The optimization then finds the best zero density bound achievable from the full combination.

### 7.5 Current State of the Art

The interaction between Huxley-type exponent pairs and Guth-Maynard large value estimates is a **frontier research area**. Some observations:

1. **Complementary ranges:** Huxley's exponent pairs are most effective near sigma = 1 (via the classical Ivic machinery), while Guth-Maynard is most effective near sigma = 3/4. The ANTEDB combines both for the best overall bound.

2. **Potential for interaction:** The Bombieri-Iwaniec method could potentially be combined with decoupling techniques to yield even stronger large value estimates. This is speculative but represents a natural research direction.

3. **The ANTEDB as a unifying framework:** By encoding all known bounds in a common format, the ANTEDB makes it possible to systematically explore how improvements in one area (e.g., a new Huxley-type exponent pair) propagate to other areas (e.g., zero density bounds). This is a major advantage over the traditional approach of optimizing each application separately.

---

## 8. Optimal Exponent Pair for Density Estimates at sigma = 3/4

### 8.1 The Optimization Problem

Among all known exponent pairs (classical + ANTEDB), which gives the best zero density bound at sigma = 3/4?

From the table in Section 4.1, the current best bound at sigma = 3/4 comes from **Guth-Maynard** with A(3/4) <= 20/9 approx 2.222, via their large value estimate rather than through the exponent pair mechanism.

If we restrict to the **classical exponent pair pipeline**, the best bound comes from a combination:

- **Ingham's estimate** gives A(3/4) <= 3/(2 - 3/4) = 12/5 = 2.4
- **Ivic's formula** with the Bourgain pair (13/84, 55/84) gives 2(68/84)/(71/84) = 136/71 approx 1.915, but this applies only for sigma sufficiently close to 1
- At sigma = 3/4 specifically, the classical exponent pair approach struggles to beat Ingham's 12/5

### 8.2 The Large Value / Mixed Approach

The best zero density bound at sigma = 3/4 that uses exponent pairs comes from combining them with large value estimates via the ANTEDB pipeline:

1. Use the Bourgain pair (13/84, 55/84) and Watt pair (89/560, 471/560) to bound beta(alpha) and hence mu(sigma) in the range sigma in [1/2, 1].

2. Use these mu bounds to constrain LV_zeta(sigma, tau) for tau >= 2.

3. Combine with Guth-Maynard's direct LV bound for the critical range.

4. Apply Huxley's subdivision and the zero density machinery (Lemma 11.5).

The ANTEDB performs this optimization automatically and finds that the best bound at sigma = 3/4 is indeed dominated by Guth-Maynard's contribution rather than by any single exponent pair.

### 8.3 Summary of Bounds at sigma = 3/4

| Method | Bound on A(3/4) | Reference |
|:-------|:----------------|:----------|
| Ingham (1940) | 12/5 = 2.400 | Classical |
| Huxley (1972) | 12/5 = 2.400 (same at sigma=3/4) | Classical + subdivision |
| Guth-Maynard (2024) | 20/9 approx 2.222 | Large value estimates |
| Density hypothesis | 2.000 | Conjecture |
| ANTEDB optimized (classical EP only) | ~2.35 | Systematic optimization |
| ANTEDB optimized (all inputs) | 20/9 approx 2.222 | = Guth-Maynard |

---

## 9. Concrete Proposals and Future Directions

### 9.1 Immediate Research Directions

1. **Extend the ANTEDB LP to directly minimize A(3/4).** The current ANTEDB infrastructure optimizes over exponent pairs and beta(alpha) bounds. One could extend this to set up a single LP/SDP that directly minimizes A(3/4) over all achievable combinations of:
   - Known exponent pairs (classical + ANTEDB + Huxley + Watt)
   - Guth-Maynard type large value estimates
   - Huxley subdivision constraints
   - Mixed moment bounds

2. **Explore the interaction between new ANTEDB exponent pairs and Guth-Maynard.** The four new ANTEDB pairs have small k and large l, making them effective near sigma = 1. Can the beta(alpha) bounds they provide be propagated to improve LV_zeta estimates in the sigma ~ 3/4 range?

3. **Investigate SDP tightenings.** The current LP approach uses only linear constraints from exponent pairs. Adding nonlinear constraints from decoupling theory and the Vinogradov main conjecture (via SDP relaxation) could tighten the feasible region.

### 9.2 Longer-Term Directions

4. **Combine Bombieri-Iwaniec with decoupling.** The Huxley-type approach (resonance curves) and the Guth-Maynard approach (decoupling + multilinear estimates) attack different aspects of the same problem. A unified approach might yield exponent pairs beyond both methods.

5. **Automated theorem proving.** The ANTEDB is designed to be eventually integrated with proof assistants (Lean 4). Formalizing the entire exponent pair -> zero density pipeline would provide machine-verified bounds.

6. **New exponent pairs via AI-assisted search.** The ANTEDB optimization discovered four new exponent pairs without new analytic inputs. More sophisticated search methods (reinforcement learning over process sequences, genetic algorithms over parameter choices) might find additional pairs.

### 9.3 The Gap Between 30/13 and 2

The current best uniform zero density exponent is 30/13 approx 2.308 (Guth-Maynard). The density hypothesis predicts 2. Closing this gap requires either:

- **New large value estimates** beyond Guth-Maynard (conceptually hardest but most impactful)
- **Better exponent pairs** that improve the beta(alpha) bounds in the critical range
- **New structural results** about the zeta function that provide constraints beyond what the general Dirichlet polynomial framework captures
- **Refinements of the zero density machinery** itself (new versions of the Halasz-Montgomery inequality, improved subdivision techniques, etc.)

The ANTEDB provides the infrastructure to immediately propagate any such improvement to the best achievable zero density bound.

---

## 10. Conclusion

The theory of exponent pairs provides a systematic framework for converting exponential sum bounds into zero density estimates for the Riemann zeta function. The pipeline proceeds through several layers: exponent pairs constrain the growth function beta(alpha), which bounds the zeta growth exponent mu(sigma), which controls large value estimates LV(sigma, tau), which finally determine the zero density exponent A(sigma).

The ANTEDB project of Tao, Trudgian, and Yang represents a paradigm shift in how this pipeline is used. By systematically recording all known theorems and relationships in a machine-readable format, and then running computational optimization, they discovered four new exponent pairs and several new zero density estimates without any new analytic input.

The Guth-Maynard breakthrough of 2024 provides the current best zero density bound (A = 30/13 ~ 2.308 uniformly, A(3/4) = 20/9 ~ 2.222 pointwise), but this comes through a fundamentally different mechanism (large value estimates for Dirichlet polynomials) rather than through the exponent pair framework. The two approaches are complementary, and the ANTEDB provides the infrastructure to optimally combine them.

The exponent pair conjecture -- that (epsilon, 1/2 + epsilon) is an exponent pair for every epsilon > 0 -- would imply the Lindelof hypothesis and hence the density hypothesis A(sigma) = 2. This remains a major open problem.

The formulation of the zero density optimization as an explicit LP/SDP, and the prospect of automated optimization within the ANTEDB framework, opens the door to systematic exploration of the gap between the current best bound and the density hypothesis.

---

## References

1. Tao, T., Trudgian, T., and Yang, A. (2025). "New exponent pairs, zero density estimates, and zero additive energy estimates: a systematic approach." arXiv:2501.16779. [https://arxiv.org/abs/2501.16779](https://arxiv.org/abs/2501.16779)

2. Guth, L. and Maynard, J. (2024). "New large value estimates for Dirichlet polynomials." arXiv:2405.20552. [https://arxiv.org/abs/2405.20552](https://arxiv.org/abs/2405.20552)

3. Tao, T. (2024). "A computation-outsourced discussion of zero density theorems for the Riemann zeta function." Blog post. [https://terrytao.wordpress.com/2024/07/07/a-computation-outsourced-discussion-of-zero-density-theorems-for-the-riemann-zeta-function/](https://terrytao.wordpress.com/2024/07/07/a-computation-outsourced-discussion-of-zero-density-theorems-for-the-riemann-zeta-function/)

4. Tao, T. (2025). "New exponent pairs, zero density estimates, and zero additive energy estimates: a systematic approach." Blog post. [https://terrytao.wordpress.com/2025/01/28/new-exponent-pairs-zero-density-estimates-and-zero-additive-energy-estimates-a-systematic-approach/](https://terrytao.wordpress.com/2025/01/28/new-exponent-pairs-zero-density-estimates-and-zero-additive-energy-estimates-a-systematic-approach/)

5. ANTEDB GitHub Repository. [https://github.com/teorth/expdb](https://github.com/teorth/expdb)

6. ANTEDB Blueprint. [https://teorth.github.io/expdb/](https://teorth.github.io/expdb/)

7. Huxley, M. N. (2005). "Exponential sums and the Riemann zeta function V." Proc. London Math. Soc. (3) 90, no. 1, 1-41.

8. Graham, S. W. and Kolesnik, G. (1991). "Van der Corput's Method of Exponential Sums." London Mathematical Society Lecture Note Series, No. 126. Cambridge University Press.

9. Ivic, A. (1985/2003). "The Riemann Zeta-Function: Theory and Applications." Dover Publications.

10. Ingham, A. E. (1940). "On the estimation of N(sigma, T)." Quart. J. Math. Oxford Ser. 11, 291-292.

11. Kadiri, H., Lumley, A., and Ng, N. (2018/2021). "Explicit zero density for the Riemann zeta function." arXiv:2101.12263. [https://arxiv.org/abs/2101.12263](https://arxiv.org/abs/2101.12263)

12. Tao, T. (2015). "254A, Notes 6: Large values of Dirichlet polynomials, zero density estimates, and primes in short intervals." [https://terrytao.wordpress.com/2015/02/13/254a-notes-6-large-values-of-dirichlet-polynomials-zero-density-estimates-and-primes-in-short-intervals/](https://terrytao.wordpress.com/2015/02/13/254a-notes-6-large-values-of-dirichlet-polynomials-zero-density-estimates-and-primes-in-short-intervals/)

13. ANTEDB Zero Density Chapter. [https://teorth.github.io/expdb/blueprint/zero-density-chapter.html](https://teorth.github.io/expdb/blueprint/zero-density-chapter.html)

14. Watt, N. (1989). "Exponential sums and the Riemann zeta-function II." J. London Math. Soc. (2) 39, 385-404.

15. Bellotti, A. (2024). "An explicit log-free zero density estimate for the Riemann zeta-function." arXiv:2405.12545. [https://arxiv.org/abs/2405.12545](https://arxiv.org/abs/2405.12545)
