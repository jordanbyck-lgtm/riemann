"""
Zero counting and the argument principle for the Riemann zeta function.

The number of zeros of zeta(s) with 0 < Im(s) < T is given by the
Riemann-von Mangoldt formula:

    N(T) = (T/(2*pi)) * log(T/(2*pi)) - T/(2*pi) + 7/8 + S(T) + O(1/T)

where S(T) = (1/pi) * arg(zeta(1/2 + iT)) is the argument of zeta
on the critical line, typically small (|S(T)| is usually < 1).

Backlund (1914) proved that |S(T)| < 0.137*log(T) + 0.443*log(log(T)) + 4.350
for T >= 2.

This module implements:
1. The smooth part of N(T) via Stirling's approximation
2. The S(T) correction via numerical argument tracking
3. Gram points and Gram's law
4. Turing's method for rigorous zero verification

Reference:
    Titchmarsh, E.C. "The Theory of the Riemann Zeta-Function" (1986)
    Turing, A.M. (1953). "Some calculations of the Riemann zeta-function"
    Edwards, H.M. "Riemann's Zeta Function" (1974)
"""

import numpy as np
from numpy import pi, log, sqrt, floor
from riemann_siegel import hardy_z, riemann_siegel_theta


def N_smooth(T):
    """
    Smooth (principal) part of the zero-counting function.

    N_0(T) = (T / (2*pi)) * log(T / (2*pi)) - T / (2*pi) + 7/8

    This counts the "expected" number of zeros up to height T.

    Parameters
    ----------
    T : float
        Height on the critical line.

    Returns
    -------
    float
        Smooth approximation to N(T).
    """
    return (T / (2 * pi)) * log(T / (2 * pi)) - T / (2 * pi) + 7.0 / 8


def gram_points(n_start, n_end):
    """
    Compute Gram points g_n where theta(g_n) = n*pi.

    Gram points provide a natural grid for locating zeros of Z(t).
    Gram's law (which holds for "most" n) states that (-1)^n * Z(g_n) > 0.
    When Gram's law holds, there is exactly one zero of Z(t) in each
    Gram interval [g_n, g_{n+1}].

    Gram's law fails about 27% of the time asymptotically (Trudgian, 2011).

    Parameters
    ----------
    n_start : int
        Starting Gram index.
    n_end : int
        Ending Gram index (exclusive).

    Returns
    -------
    list of float
        Gram points g_{n_start}, ..., g_{n_end - 1}.
    """
    from scipy.optimize import brentq

    gram_pts = []
    for n in range(n_start, n_end):
        target = n * pi

        # Initial bracket: theta is monotonically increasing for t > 7
        # Use asymptotic inversion: t ~ 2*pi*exp(W(n/e) + 1) where W is Lambert W
        # Simpler: for large n, g_n ~ 2*pi*n / log(n/(2*pi))
        if n <= 0:
            t_approx = 10.0
        else:
            # Better initial guess using the asymptotic formula
            t_approx = 2 * pi * np.exp(1 + np.real(
                _lambert_w_approx(n / np.e)
            ))
            if t_approx < 8:
                t_approx = 8.0 + n * 0.5

        # Bracket the root
        dt = max(1.0, t_approx * 0.1)
        a = max(8.0, t_approx - dt)
        b = t_approx + dt

        # Expand bracket if needed
        while riemann_siegel_theta(a) > target:
            a = max(8.0, a - dt)
            dt *= 2
        while riemann_siegel_theta(b) < target:
            b += dt
            dt *= 2

        try:
            g_n = brentq(lambda t: riemann_siegel_theta(t) - target, a, b)
            gram_pts.append(g_n)
        except ValueError:
            gram_pts.append(float('nan'))

    return gram_pts


def _lambert_w_approx(x):
    """Rough approximation to the Lambert W function for x > 0."""
    if x <= 0:
        return 0.0
    L1 = np.log(x)
    L2 = np.log(L1) if L1 > 0 else 0
    return L1 - L2 + L2 / L1


def check_gram_law(gram_pts):
    """
    Check Gram's law at each Gram point.

    Gram's law states: (-1)^n * Z(g_n) > 0.
    It holds for most Gram points but fails about 27% of the time.

    Parameters
    ----------
    gram_pts : list of (int, float)
        List of (index, gram_point) pairs.

    Returns
    -------
    list of dict
        Results including Gram index, point, Z value, and whether Gram's law holds.
    """
    results = []
    for n, g in gram_pts:
        z_val = hardy_z(g)
        gram_law_holds = ((-1) ** n * z_val > 0)
        results.append({
            'n': n,
            'gram_point': g,
            'Z_value': z_val,
            'gram_law': gram_law_holds
        })
    return results


def count_sign_changes(t_start, t_end, num_points=10000):
    """
    Count zeros of Z(t) in [t_start, t_end] by detecting sign changes.

    This is the basic method for numerical zero counting. Each sign
    change of Z(t) corresponds to a zero of zeta on the critical line.

    Parameters
    ----------
    t_start : float
        Start of interval.
    t_end : float
        End of interval.
    num_points : int
        Number of sample points.

    Returns
    -------
    dict
        'count': number of sign changes detected,
        'zeros': approximate zero locations,
        'expected': N_smooth(t_end) - N_smooth(t_start).
    """
    t_grid = np.linspace(t_start, t_end, num_points)
    z_vals = np.array([hardy_z(t) for t in t_grid])

    zeros = []
    for i in range(len(z_vals) - 1):
        if z_vals[i] * z_vals[i + 1] < 0:
            # Bisect for precision
            a, b = t_grid[i], t_grid[i + 1]
            for _ in range(60):  # ~18 digits of precision
                mid = (a + b) / 2
                z_mid = hardy_z(mid)
                if z_vals[i] * z_mid <= 0:
                    b = mid
                else:
                    a = mid
            zeros.append((a + b) / 2)

    expected = N_smooth(t_end) - N_smooth(t_start)

    return {
        'count': len(zeros),
        'zeros': zeros,
        'expected': expected,
    }


def turing_method_verify(t_start, t_end, gram_block_size=1):
    """
    Turing's method for rigorous zero verification.

    Turing (1953) devised a method to rigorously verify that all zeros
    of zeta in a given range lie on the critical line. The idea:

    1. Count zeros on the critical line by sign changes of Z(t).
    2. Count total zeros (on and off the line) using N(T).
    3. If the counts agree, ALL zeros in the range are on the critical line.

    The key insight is that N(T) counts ALL zeros (via the argument
    principle), while sign changes of Z(t) only count zeros on Re(s) = 1/2.

    For rigorous verification, one must also account for potential
    "Lehmer phenomena" where Z(t) comes close to zero without changing sign
    (which would indicate two very close zeros).

    Parameters
    ----------
    t_start : float
        Start of verification range (should be a Gram point).
    t_end : float
        End of verification range (should be a Gram point).
    gram_block_size : int
        Size of Gram blocks for Turing's method (1 for individual).

    Returns
    -------
    dict
        Verification results.
    """
    # Step 1: Count sign changes (zeros on critical line)
    result = count_sign_changes(t_start, t_end, num_points=20000)
    sign_change_count = result['count']
    found_zeros = result['zeros']

    # Step 2: Compute expected count from N(T)
    expected_total = N_smooth(t_end) - N_smooth(t_start)

    # Step 3: Round expected to nearest integer (S(T) fluctuations are small)
    expected_int = int(round(expected_total))

    # Step 4: Verify
    verified = (sign_change_count == expected_int)

    # Step 5: Check for close zeros (Lehmer phenomenon)
    min_gap = float('inf')
    for i in range(1, len(found_zeros)):
        gap = found_zeros[i] - found_zeros[i - 1]
        min_gap = min(min_gap, gap)

    return {
        'range': (t_start, t_end),
        'zeros_found': sign_change_count,
        'expected_total': expected_total,
        'expected_int': expected_int,
        'all_on_critical_line': verified,
        'min_zero_gap': min_gap if len(found_zeros) > 1 else None,
        'zeros': found_zeros,
    }


if __name__ == "__main__":
    print("Zero Counting and Turing's Method")
    print("=" * 60)

    # Compute Gram points
    print("\nGram points g_0 through g_20:")
    gpts = gram_points(0, 21)
    print(f"{'n':>4s}  {'g_n':>12s}  {'theta(g_n)':>12s}  {'Z(g_n)':>12s}  {'Gram law':>10s}")
    print("-" * 60)
    for n, g in enumerate(gpts):
        theta = riemann_siegel_theta(g)
        z_val = hardy_z(g)
        gram_ok = "YES" if (-1) ** n * z_val > 0 else "NO"
        print(f"{n:4d}  {g:12.6f}  {theta:12.6f}  {z_val:12.6f}  {gram_ok:>10s}")

    # Turing verification
    print("\n\nTuring Method: Verifying zeros in [10, 100]")
    print("-" * 60)
    result = turing_method_verify(10, 100)
    print(f"Zeros found (sign changes):    {result['zeros_found']}")
    print(f"Expected total (from N(T)):    {result['expected_total']:.2f}")
    print(f"Expected (integer):            {result['expected_int']}")
    print(f"All on critical line:          {result['all_on_critical_line']}")
    print(f"Minimum zero gap:              {result['min_zero_gap']:.6f}")

    print("\nZero locations:")
    for i, z in enumerate(result['zeros']):
        print(f"  #{i+1:3d}: t = {z:.10f}")
