"""
Keating-Snaith Moment Conjectures from Random Matrix Theory.

In 2000, Keating and Snaith used random matrix theory to predict the
asymptotic behavior of moments of the Riemann zeta function on the
critical line:

    (1/T) * int_0^T |zeta(1/2 + it)|^{2k} dt ~ C_k * (log T)^{k^2}

where C_k involves a product of the Barnes G-function:
    C_k = g_k * a_k

with:
    g_k = G(1+k)^2 / G(1+2k)    (from RMT)

    a_k = prod_{p prime} (1-1/p)^{k^2} * sum_{m=0}^{inf} (Gamma(m+k)/(m!*Gamma(k)))^2 / p^m
                                          (arithmetic factor)

Previously, only k=1 (Hardy-Littlewood, 1918) and k=2 (Ingham, 1926)
were known:
    C_1 = 1          (mean square of zeta)
    C_2 = 1/(2*pi^2) (fourth moment)

The Keating-Snaith conjecture provides predictions for ALL positive real k,
including the striking prediction for k=6:
    C_6 = 42         (the answer to life, the universe, and everything!)

These conjectures have been verified numerically and have generated
deep connections between number theory and random matrix theory.

They also lead to conjectures about:
- Extreme values of zeta (maximum of |zeta| in intervals)
- Moments of L-functions (generalized Keating-Snaith)
- Quantum unique ergodicity

Reference:
    Keating, J.P. & Snaith, N.C. (2000). "Random matrix theory and
        ζ(1/2+it)" Comm. Math. Phys. 214(1), 57-89.
    Conrey, J.B., Farmer, D.W., Keating, J.P., Rubinstein, M.O., &
        Snaith, N.C. (2005). "Integral moments of L-functions"
"""

import numpy as np
from numpy import pi, log, exp, sqrt
from scipy.special import gamma as gamma_func
from scipy.special import gammaln


def barnes_g(n):
    """
    Barnes G-function G(n) for positive integer n.

    G(1) = 1
    G(n+1) = Gamma(n) * G(n)

    So G(n) = prod_{k=1}^{n-2} Gamma(k+1) = prod_{k=1}^{n-2} k!

    For non-integer arguments, uses the relation:
        log G(z+1) = (z/2)*log(2*pi) - (z + z^2*(1+gamma))/2
                     + z*log(Gamma(z+1)) - log(Gamma(z+1))
                     + integral...

    Parameters
    ----------
    n : int
        Positive integer argument.

    Returns
    -------
    float
        G(n).
    """
    if n <= 0:
        return 0.0
    if n == 1 or n == 2:
        return 1.0

    # G(n) = product of factorials
    result = 1.0
    for k in range(1, n - 1):
        result *= gamma_func(k + 1)
    return result


def log_barnes_g(z):
    """
    Logarithm of the Barnes G-function for real z > 0.

    Uses the asymptotic/integral representation for non-integer z.

    For integer z = n:
        log G(n) = sum_{k=1}^{n-2} log(k!)

    Parameters
    ----------
    z : float
        Positive real argument.

    Returns
    -------
    float
        log G(z).
    """
    if z <= 0:
        return float('-inf')

    # For integer values
    n = int(round(z))
    if abs(z - n) < 1e-10 and n > 0:
        if n <= 2:
            return 0.0
        result = 0.0
        for k in range(1, n - 1):
            result += gammaln(k + 1)
        return result

    # For non-integer, use the identity:
    # log G(z+1) = (z/2)*log(2*pi) + integral representation
    # Approximate using the Stirling-type expansion:
    # log G(z+1) ~ (z^2/2)*log(z) - 3*z^2/4 + (z/2)*log(2*pi)
    #              - (1/12)*log(z) + zeta'(-1) + ...
    zeta_prime_minus1 = -0.16542114370045092  # zeta'(-1)

    if z > 10:
        return ((z * z - 1.0 / 12) * log(z) / 2
                - 3 * z * z / 4
                + z * log(2 * pi) / 2
                + zeta_prime_minus1
                + 1.0 / (720 * z * z))
    else:
        # Use recurrence G(z+1) = Gamma(z)*G(z) to reduce to large argument
        log_g = 0.0
        zz = z
        while zz < 12:
            log_g -= gammaln(zz)
            zz += 1
        log_g += ((zz * zz - 1.0 / 12) * log(zz) / 2
                  - 3 * zz * zz / 4
                  + zz * log(2 * pi) / 2
                  + zeta_prime_minus1)
        return log_g


def rmt_factor(k):
    """
    Random matrix theory factor g_k in the moment conjecture.

    g_k = G(1+k)^2 / G(1+2k)

    where G is the Barnes G-function.

    Parameters
    ----------
    k : float
        Moment parameter (positive real).

    Returns
    -------
    float
        g_k value.
    """
    log_g = 2 * log_barnes_g(1 + k) - log_barnes_g(1 + 2 * k)
    return exp(log_g)


def arithmetic_factor(k, num_primes=100):
    """
    Arithmetic factor a_k in the Keating-Snaith conjecture.

    a_k = prod_{p prime} (1-1/p)^{k^2}
          * sum_{m=0}^{inf} (Gamma(m+k)/(m!*Gamma(k)))^2 / p^m

    Parameters
    ----------
    k : float
        Moment parameter.
    num_primes : int
        Number of primes to include in the Euler product.

    Returns
    -------
    float
        a_k value.
    """
    primes = _sieve_primes(num_primes)

    log_product = 0.0
    for p in primes:
        # (1 - 1/p)^{k^2}
        log_product += k * k * log(1 - 1.0 / p)

        # Compute the local sum
        local_sum = 0.0
        p_power = 1.0  # p^m
        for m in range(50):
            # Gamma(m+k) / (m! * Gamma(k)) = binomial-like coefficient
            # = (k)(k+1)...(k+m-1) / m! = C(m+k-1, m)
            if m == 0:
                coeff = 1.0
            else:
                coeff *= (k + m - 1) / m
            local_sum += coeff * coeff / p_power
            p_power *= p

            if coeff * coeff / p_power < 1e-15:
                break

        log_product += log(local_sum)

    return exp(log_product)


def moment_constant(k, num_primes=100):
    """
    Full Keating-Snaith moment constant C_k = g_k * a_k.

    Predicts: (1/T) * int_0^T |zeta(1/2+it)|^{2k} dt ~ C_k * (log T)^{k^2}

    Known/conjectured values:
        C_1 = 1 (proven: Hardy-Littlewood)
        C_2 = 1/(2*pi^2) ≈ 0.05066 (proven: Ingham)
        C_3 ≈ 0.000425 (conjectured)
        C_4 ≈ 2.4 * 10^{-6} (conjectured)

    Parameters
    ----------
    k : float
        Moment parameter.
    num_primes : int
        Primes for the arithmetic factor.

    Returns
    -------
    dict
        'C_k': the constant,
        'g_k': RMT factor,
        'a_k': arithmetic factor,
        'power': k^2 (the power of log T).
    """
    g = rmt_factor(k)
    a = arithmetic_factor(k, num_primes)
    return {
        'C_k': g * a,
        'g_k': g,
        'a_k': a,
        'power': k * k,
    }


def numerical_moment(k, T, num_points=10000):
    """
    Numerically compute the 2k-th moment of zeta on the critical line.

    M_k(T) = (1/T) * int_0^T |zeta(1/2 + it)|^{2k} dt

    Uses the Hardy Z-function: |zeta(1/2+it)| = |Z(t)|.

    Parameters
    ----------
    k : float
        Moment parameter.
    T : float
        Height limit.
    num_points : int
        Quadrature points.

    Returns
    -------
    dict
        'moment': numerical value,
        'predicted': C_k * (log T)^{k^2},
        'ratio': moment / predicted.
    """
    from riemann_siegel import hardy_z

    t_vals = np.linspace(max(10, T * 0.01), T, num_points)
    dt = t_vals[1] - t_vals[0]

    z_vals = np.array([hardy_z(t) for t in t_vals])
    moment_val = np.mean(np.abs(z_vals) ** (2 * k))

    # Prediction
    ck = moment_constant(k)
    predicted = ck['C_k'] * log(T) ** (k * k)

    return {
        'moment': moment_val,
        'predicted': predicted,
        'ratio': moment_val / predicted if predicted > 0 else float('inf'),
        'C_k': ck['C_k'],
        'power': k * k,
        'T': T,
    }


def _sieve_primes(count):
    """Return the first `count` prime numbers."""
    if count <= 0:
        return []
    primes = [2]
    candidate = 3
    while len(primes) < count:
        is_prime = True
        for p in primes:
            if p * p > candidate:
                break
            if candidate % p == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(candidate)
        candidate += 2
    return primes[:count]


if __name__ == "__main__":
    print("Keating-Snaith Moment Conjectures")
    print("=" * 60)

    print("\nMoment constants C_k = g_k * a_k:")
    print(f"{'k':>5s}  {'g_k (RMT)':>12s}  {'a_k (arith)':>12s}  "
          f"{'C_k':>14s}  {'(log T)^k^2 power':>8s}")
    print("-" * 65)

    for k in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]:
        mc = moment_constant(k, num_primes=200)
        print(f"{k:5.1f}  {mc['g_k']:12.6e}  {mc['a_k']:12.6e}  "
              f"{mc['C_k']:14.6e}  {mc['power']:8.1f}")

    print(f"\nKnown exact values:")
    print(f"  C_1 = 1 (Hardy-Littlewood 1918)")
    mc1 = moment_constant(1.0, num_primes=500)
    print(f"  Our computation: C_1 = {mc1['C_k']:.8f}")

    print(f"  C_2 = 1/(2*pi^2) = {1/(2*pi**2):.8f} (Ingham 1926)")
    mc2 = moment_constant(2.0, num_primes=500)
    print(f"  Our computation: C_2 = {mc2['C_k']:.8f}")

    # Numerical moment computation
    print("\n\nNumerical moment computation (T=1000):")
    for k in [1.0, 2.0]:
        result = numerical_moment(k, T=1000, num_points=20000)
        print(f"  k={k:.0f}: moment={result['moment']:.6f}, "
              f"predicted={result['predicted']:.6f}, "
              f"ratio={result['ratio']:.4f}")
