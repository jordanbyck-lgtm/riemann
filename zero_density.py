"""
Zero Density Estimates for the Riemann Zeta Function.

A zero density estimate bounds the number of zeros N(sigma, T) of
zeta(s) with Re(s) >= sigma and 0 < Im(s) <= T. The RH asserts
N(sigma, T) = 0 for sigma > 1/2, but unconditional bounds take the form:

    N(sigma, T) <= C * T^{A(sigma)*(1-sigma) + epsilon}

The exponent A(sigma) has been refined over the past century:

    A(sigma) = 3         Ingham (1940) for sigma near 3/4
    A(sigma) = 12/5      Huxley (1972) for sigma near 3/4
    A(sigma) = 30/13     Guth & Maynard (2024) for sigma near 3/4

Guth-Maynard Breakthrough (May 2024):
======================================
Larry Guth and James Maynard proved:

    N(sigma, T) <= T^{(30/13)(1-sigma) + o(1)}

This is the first improvement to Ingham's density bound at sigma = 3/4
in OVER 80 YEARS. Terence Tao called it "a remarkable breakthrough
towards the Riemann hypothesis."

Their method involves:
1. New bounds on large values of Dirichlet polynomials
2. A novel "efficient congruencing" / decoupling approach
3. Connections to additive combinatorics and restriction estimates

Consequences:
- Primes exist in all intervals [x, x + x^{17/30}] for large x
  (improving x^{7/12} from Huxley 1972)
- Better bounds on gaps between primes
- Improved error terms in the prime number theorem

This module also implements:
- Baluyot-Goldston-Suriajaya-Turnage-Butterbaugh (2025): Montgomery's
  pair correlation conjecture implies 100% simple zeros on the critical
  line WITHOUT assuming RH

Reference:
    Ingham, A.E. (1940). "On the estimation of N(sigma,T)"
    Huxley, M.N. (1972). "On the difference between consecutive primes"
    Guth, L. & Maynard, J. (2024). "New large value estimates for
        Dirichlet polynomials" arXiv:2405.20552
    Baluyot, S., Goldston, D., Suriajaya, A. & Turnage-Butterbaugh, C.
        (2025). "Pair correlation and the distribution of zeros"
        arXiv:2503.15449
"""

import numpy as np
from numpy import pi, log, exp


# ============================================================
# Zero Density Exponents A(sigma)
# ============================================================

def ingham_exponent(sigma):
    """
    Ingham's zero density exponent (1940).

    N(sigma, T) <= T^{3*(1-sigma) + epsilon}

    This held as the best bound near sigma = 3/4 for 80 years.

    Parameters
    ----------
    sigma : float or ndarray
        Real part threshold, 1/2 < sigma < 1.

    Returns
    -------
    float or ndarray
        Exponent 3*(1-sigma).
    """
    sigma = np.asarray(sigma, dtype=np.float64)
    return 3.0 * (1 - sigma)


def huxley_exponent(sigma):
    """
    Huxley's zero density exponent (1972).

    N(sigma, T) <= T^{(12/5)*(1-sigma) + epsilon}

    Improvement over Ingham for sigma near 3/4, via sieve methods.

    Parameters
    ----------
    sigma : float or ndarray
        Real part threshold.

    Returns
    -------
    float or ndarray
        Exponent (12/5)*(1-sigma).
    """
    sigma = np.asarray(sigma, dtype=np.float64)
    return (12.0 / 5) * (1 - sigma)


def guth_maynard_exponent(sigma):
    """
    Guth-Maynard zero density exponent (2024).

    N(sigma, T) <= T^{(30/13)*(1-sigma) + o(1)}

    First improvement at sigma = 3/4 since Ingham (1940).
    Uses new large value estimates for Dirichlet polynomials
    via decoupling methods from harmonic analysis.

    Parameters
    ----------
    sigma : float or ndarray
        Real part threshold.

    Returns
    -------
    float or ndarray
        Exponent (30/13)*(1-sigma).
    """
    sigma = np.asarray(sigma, dtype=np.float64)
    return (30.0 / 13) * (1 - sigma)


def density_hypothesis_exponent(sigma):
    """
    The Density Hypothesis exponent (conjectural).

    N(sigma, T) <= T^{2*(1-sigma) + epsilon}

    This is the strongest conjecture short of RH. It would follow
    from the Lindelöf Hypothesis and implies many of the same
    consequences for prime gaps.

    Parameters
    ----------
    sigma : float or ndarray
        Real part threshold.

    Returns
    -------
    float or ndarray
        Exponent 2*(1-sigma).
    """
    sigma = np.asarray(sigma, dtype=np.float64)
    return 2.0 * (1 - sigma)


def compare_density_bounds(sigma_values=None):
    """
    Compare all known zero density exponents.

    Displays how the bound N(sigma, T) <= T^{A*(1-sigma)} has
    improved from Ingham (1940) through Guth-Maynard (2024).

    Parameters
    ----------
    sigma_values : array-like or None
        Values of sigma to compare at.

    Returns
    -------
    dict
        Comparison data.
    """
    if sigma_values is None:
        sigma_values = np.linspace(0.55, 0.95, 20)

    sigma = np.asarray(sigma_values)

    return {
        'sigma': sigma,
        'ingham': ingham_exponent(sigma),
        'huxley': huxley_exponent(sigma),
        'guth_maynard': guth_maynard_exponent(sigma),
        'density_hypothesis': density_hypothesis_exponent(sigma),
    }


# ============================================================
# Consequences for Prime Gaps
# ============================================================

def prime_gap_exponent_from_density(A):
    """
    Convert a zero density exponent A to a prime gap exponent.

    If N(sigma, T) <= T^{A*(1-sigma)}, then for theta = 1 - 1/A,
    every interval [x, x + x^theta] contains a prime for large x.

    Ingham (1940), A=3:   theta = 2/3
    Huxley (1972), A=12/5: theta = 7/12
    Guth-Maynard (2024), A=30/13: theta = 17/30

    The Density Hypothesis (A=2): theta = 1/2
    RH: theta = 1/2 + epsilon

    Parameters
    ----------
    A : float
        Zero density exponent.

    Returns
    -------
    float
        Prime gap exponent theta.
    """
    return 1 - 1.0 / A


def demonstrate_prime_gap_improvement():
    """
    Show the progression of prime gap results from density estimates.

    Returns
    -------
    list of dict
        Historical progression.
    """
    results = [
        {'year': 1940, 'author': 'Ingham', 'A': 3.0,
         'theta': 1 - 1 / 3.0, 'note': 'First density estimate'},
        {'year': 1972, 'author': 'Huxley', 'A': 12 / 5,
         'theta': 1 - 5 / 12.0, 'note': 'Sieve methods'},
        {'year': 2024, 'author': 'Guth-Maynard', 'A': 30 / 13,
         'theta': 1 - 13 / 30.0,
         'note': 'Decoupling + large values (BREAKTHROUGH)'},
        {'year': None, 'author': 'Density Hypothesis', 'A': 2.0,
         'theta': 0.5, 'note': 'Conjectural'},
        {'year': None, 'author': 'Riemann Hypothesis', 'A': float('inf'),
         'theta': 0.5, 'note': 'Would give theta = 1/2 + epsilon'},
    ]

    return results


# ============================================================
# Pair Correlation and Simple Zeros (Baluyot et al. 2025)
# ============================================================

def pair_correlation_implies_simplicity():
    """
    Summarize the Baluyot-Goldston-Suriajaya-Turnage-Butterbaugh
    result (March 2025).

    They proved: Montgomery's Pair Correlation Conjecture (PCC)
    implies that asymptotically 100% of zeta zeros are:
      (a) Simple (no repeated zeros)
      (b) On the critical line Re(s) = 1/2

    This does NOT assume the Riemann Hypothesis.

    Previous results:
    - Gallagher & Mueller (1978): PCC + RH => simple zeros
    - The 2025 result removes the RH assumption using "symmetric
      diagonal terms" that extract horizontal distribution info
      from the pair correlation function.

    Returns
    -------
    dict
        Summary of the result.
    """
    return {
        'theorem': (
            "If Montgomery's Pair Correlation Conjecture holds, then "
            "asymptotically 100% of the non-trivial zeros of zeta(s) "
            "are simple and lie on the critical line Re(s) = 1/2."
        ),
        'authors': 'Baluyot, Goldston, Suriajaya, Turnage-Butterbaugh',
        'year': 2025,
        'arxiv': '2503.15449',
        'key_innovation': (
            "Symmetric diagonal terms in the pair correlation sum "
            "provide information about horizontal zero distribution "
            "without assuming RH."
        ),
        'previous': (
            "Gallagher & Mueller (1978) proved PCC + RH => simplicity. "
            "The 2025 result removes the RH assumption entirely."
        ),
    }


def estimate_N_sigma(sigma, T, method='guth_maynard'):
    """
    Estimate the number of zeros with Re(s) >= sigma, 0 < Im(s) <= T.

    Uses the chosen zero density bound. Note these are UPPER bounds;
    the true count could be much smaller (and is 0 if RH is true).

    Parameters
    ----------
    sigma : float
        Real part threshold (> 1/2).
    T : float
        Height.
    method : str
        'ingham', 'huxley', 'guth_maynard', or 'density_hypothesis'.

    Returns
    -------
    float
        Upper bound on N(sigma, T).
    """
    exponent_funcs = {
        'ingham': ingham_exponent,
        'huxley': huxley_exponent,
        'guth_maynard': guth_maynard_exponent,
        'density_hypothesis': density_hypothesis_exponent,
    }

    func = exponent_funcs.get(method, guth_maynard_exponent)
    exponent = func(sigma)

    return T ** exponent


if __name__ == "__main__":
    print("Zero Density Estimates for the Riemann Zeta Function")
    print("=" * 70)

    # Compare bounds
    print("\nZero density exponent A(sigma) such that N(sigma,T) <= T^{A(1-sigma)}:")
    print(f"\n{'sigma':>6s}  {'Ingham':>10s}  {'Huxley':>10s}  "
          f"{'Guth-May':>10s}  {'DensHyp':>10s}")
    print(f"{'':>6s}  {'(1940)':>10s}  {'(1972)':>10s}  "
          f"{'(2024)':>10s}  {'(conj)':>10s}")
    print("-" * 52)

    data = compare_density_bounds()
    for i, s in enumerate(data['sigma']):
        print(f"{s:6.3f}  {data['ingham'][i]:10.4f}  {data['huxley'][i]:10.4f}  "
              f"{data['guth_maynard'][i]:10.4f}  {data['density_hypothesis'][i]:10.4f}")

    # Prime gap consequences
    print("\n\nConsequences for primes in short intervals [x, x + x^theta]:")
    print("-" * 70)
    results = demonstrate_prime_gap_improvement()
    for r in results:
        year = str(r['year']) if r['year'] else '????'
        print(f"  {year}  {r['author']:<22s}  A={r['A']:<6.2f}  "
              f"theta={r['theta']:.4f}  {r['note']}")

    # The sigma = 3/4 case (where the breakthrough happened)
    print(f"\n\nAt sigma = 3/4 (the critical point):")
    print(f"  Ingham (1940):      exponent = {ingham_exponent(0.75):.4f}")
    print(f"  Huxley (1972):      exponent = {huxley_exponent(0.75):.4f}")
    print(f"  Guth-Maynard (2024): exponent = {guth_maynard_exponent(0.75):.4f}")
    print(f"  Improvement: {ingham_exponent(0.75) - guth_maynard_exponent(0.75):.4f} "
          f"({100*(1-guth_maynard_exponent(0.75)/ingham_exponent(0.75)):.1f}% reduction)")

    # Numerical comparison: how many off-line zeros are allowed?
    print(f"\n\nUpper bounds on N(3/4, T) for various T:")
    print(f"{'T':>12s}  {'Ingham':>14s}  {'Guth-Maynard':>14s}  {'Ratio':>10s}")
    print("-" * 56)
    for logT in [10, 15, 20, 25, 30]:
        T = 10.0 ** logT
        n_old = estimate_N_sigma(0.75, T, 'ingham')
        n_new = estimate_N_sigma(0.75, T, 'guth_maynard')
        ratio = n_old / n_new if n_new > 0 else float('inf')
        print(f"  10^{logT:<5d}  {n_old:14.2e}  {n_new:14.2e}  {ratio:10.1f}x")

    # Baluyot et al. result
    print("\n\nBaluyot-Goldston-Suriajaya-Turnage-Butterbaugh (2025):")
    print("-" * 70)
    result = pair_correlation_implies_simplicity()
    print(f"  {result['theorem']}")
    print(f"\n  Key innovation: {result['key_innovation']}")
    print(f"  Previous: {result['previous']}")
