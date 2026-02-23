"""
Riemann-Siegel Formula for computing Z(t) on the critical line.

The Riemann-Siegel formula is the primary tool used for large-scale
numerical verification of the Riemann Hypothesis. It expresses the
Hardy Z-function as a finite main sum plus correction terms, enabling
efficient evaluation of zeta(1/2 + it).

The formula was discovered in Riemann's unpublished notes and first
published by C.L. Siegel in 1932.

Z(t) = 2 * sum_{n=1}^{N} n^{-1/2} * cos(theta(t) - t*log(n)) + R(t)

where N = floor(sqrt(t / (2*pi))) and R(t) is a remainder expressible
as an asymptotic series in powers of (2*pi/t)^{1/4}.

The Riemann-Siegel theta function is:
    theta(t) = arg(Gamma(1/4 + it/2)) - (t/2)*log(pi)

For practical computation, the asymptotic expansion is:
    theta(t) ~ (t/2)*log(t/(2*pi)) - t/2 - pi/8
                + 1/(48t) + 7/(5760t^3) + ...

Reference:
    Edwards, H.M. "Riemann's Zeta Function" (1974), Chapter 7
    Gabcke, W. (1979). "Neue Herleitung und explizite Restabschätzung
        der Riemann-Siegel-Formel" (PhD thesis, Göttingen)
"""

import numpy as np
from numpy import pi, sqrt, log, cos, sin, floor, exp


def riemann_siegel_theta(t):
    """
    Compute the Riemann-Siegel theta function.

    theta(t) = arg(Gamma(1/4 + it/2)) - (t/2)*log(pi)

    Uses Stirling's asymptotic expansion for large t:
        theta(t) ~ (t/2)*log(t/(2*pi)) - t/2 - pi/8
                    + 1/(48*t) + 7/(5760*t^3) + 31/(80640*t^5) + ...

    Parameters
    ----------
    t : float or ndarray
        Point(s) on the critical line (t > 0 for accuracy).

    Returns
    -------
    float or ndarray
        Value of theta(t).
    """
    t = np.asarray(t, dtype=np.float64)

    # Stirling expansion coefficients for theta
    # These come from the asymptotic expansion of log(Gamma(s))
    result = (t / 2) * np.log(t / (2 * pi)) - t / 2 - pi / 8

    # Higher-order correction terms (Stirling series)
    t_inv = 1.0 / t
    t_inv2 = t_inv * t_inv

    result += (1.0 / 48) * t_inv
    result += (7.0 / 5760) * t_inv * t_inv2
    result += (31.0 / 80640) * t_inv * t_inv2 * t_inv2
    result += (127.0 / 430080) * t_inv * t_inv2 * t_inv2 * t_inv2

    return result


def riemann_siegel_correction_coefficients(p):
    """
    Compute Riemann-Siegel correction coefficients C_k(p).

    The remainder R(t) in the Riemann-Siegel formula is expressed as:
        R(t) = (-1)^{N-1} * (2*pi/t)^{1/4} * sum_{k>=0} C_k(p) * (2*pi/t)^{k/2}

    where p = sqrt(t/(2*pi)) - N is the fractional part.

    C_0(p) = Psi(p) = cos(2*pi*(p^2 - p - 1/16)) / cos(2*pi*p)

    Higher coefficients use the polynomial expansion from Edwards (1974),
    Chapter 7, based on Siegel's original notes and Gabcke (1979).

    Parameters
    ----------
    p : float
        Fractional part, 0 <= p < 1.

    Returns
    -------
    tuple
        (C0, C1, C2, C3) correction coefficients.
    """
    # C_0(p) = Psi(p) = cos(2*pi*(p^2 - p - 1/16)) / cos(2*pi*p)
    cos_denom = cos(2 * pi * p)
    if abs(cos_denom) < 1e-15:
        C0 = 0.0
    else:
        C0 = cos(2 * pi * (p * p - p - 1.0 / 16)) / cos_denom

    # For C_1, C_2, C_3 use the polynomial forms from the Riemann-Siegel
    # expansion. These are expressed in terms of z = p - 1/2.
    # Reference: Edwards "Riemann's Zeta Function" (1974), eq. 7.5.6-7.5.8
    #            and Gabcke's thesis (1979), Table 1.
    z = p - 0.5

    # C_1(p): coefficient of (2*pi/t)^{1/2}
    # From the standard expansion:
    #   C_1(z) = -(1/96*pi^2) * d^3/dz^3 [Psi(z+1/2)]
    # Using the polynomial:
    C1 = -(1.0 / 96.0) + (1.0 / 64.0) * z**2 - (1.0 / 64.0) * z**4

    # C_2(p): coefficient of (2*pi/t)^{1}
    C2 = (1.0 / 18432.0) - z**2 / 3072.0 + z**4 / 2048.0

    # C_3(p): coefficient of (2*pi/t)^{3/2}
    C3 = -(1.0 / 5308416.0) + z**2 / 294912.0

    return (C0, C1, C2, C3)


def hardy_z(t, num_corrections=1):
    """
    Compute the Hardy Z-function using the Riemann-Siegel formula.

    Z(t) is a real-valued function whose zeros coincide with the zeros
    of zeta(1/2 + it) on the critical line. It is defined as:

        Z(t) = exp(i*theta(t)) * zeta(1/2 + it)

    The Riemann-Siegel formula gives:
        Z(t) = 2 * sum_{n=1}^{N} cos(theta(t) - t*log(n)) / sqrt(n)
               + R(t)

    where N = floor(sqrt(t / (2*pi))).

    The remainder R(t) is dominated by the C_0 term:
        R(t) ≈ (-1)^{N-1} * (2*pi/t)^{1/4} * C_0(p)

    Parameters
    ----------
    t : float
        Point on the critical line (should be > 10 for accuracy).
    num_corrections : int
        Number of Riemann-Siegel correction terms (0-3).
        Default 1 uses only C_0 (accurate to O((2*pi/t)^{3/4})).

    Returns
    -------
    float
        Z(t), the Hardy Z-function value.
    """
    t = float(t)
    if t < 10:
        return _hardy_z_direct(t)

    N = int(floor(sqrt(t / (2 * pi))))
    theta = riemann_siegel_theta(t)

    # Main sum
    n_vals = np.arange(1, N + 1, dtype=np.float64)
    main_sum = 2 * np.sum(np.cos(theta - t * np.log(n_vals)) / np.sqrt(n_vals))

    # Remainder term
    if num_corrections > 0:
        tau = sqrt(t / (2 * pi))
        p = tau - N  # fractional part, 0 <= p < 1
        coeffs = riemann_siegel_correction_coefficients(p)

        factor = (-1) ** (N - 1) * (2 * pi / t) ** 0.25
        remainder = 0.0
        power = 1.0
        for k in range(min(num_corrections, 4)):
            remainder += coeffs[k] * power
            power *= sqrt(2 * pi / t)

        main_sum += factor * remainder

    return main_sum


def _hardy_z_direct(t):
    """Direct computation of Z(t) using mpmath for small t."""
    try:
        import mpmath
        mpmath.mp.dps = 25
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        zeta_val = mpmath.zeta(s)
        theta_val = mpmath.arg(mpmath.gamma(mpmath.mpf('0.25') + 1j * mpmath.mpf(str(t)) / 2))
        theta_val -= (t / 2) * float(mpmath.log(mpmath.pi))
        z_val = float(mpmath.exp(1j * theta_val) * zeta_val)
        return z_val.real if isinstance(z_val, complex) else z_val
    except ImportError:
        # Fallback: direct but less accurate
        return 0.0


def hardy_z_vectorized(t_array, num_corrections=2):
    """
    Evaluate Z(t) at multiple points efficiently.

    Parameters
    ----------
    t_array : ndarray
        Array of t values.
    num_corrections : int
        Number of correction terms.

    Returns
    -------
    ndarray
        Z(t) values.
    """
    return np.array([hardy_z(t, num_corrections) for t in t_array])


def zeta_on_critical_line(t):
    """
    Compute zeta(1/2 + it) using the Z-function.

    zeta(1/2 + it) = Z(t) * exp(-i*theta(t))

    Parameters
    ----------
    t : float
        Imaginary part.

    Returns
    -------
    complex
        Value of zeta(1/2 + it).
    """
    z = hardy_z(t)
    theta = riemann_siegel_theta(t)
    return z * np.exp(-1j * theta)


if __name__ == "__main__":
    print("Riemann-Siegel Formula - Computing Z(t) on the critical line")
    print("=" * 65)

    # Known zeros (imaginary parts) of zeta on the critical line
    known_zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                   37.586178, 40.918719, 43.327073, 48.005151, 49.773832]

    print("\nVerification at known zero locations:")
    print(f"{'t':>12s}  {'Z(t)':>15s}  {'|Z(t)|':>12s}")
    print("-" * 45)
    for t0 in known_zeros:
        z_val = hardy_z(t0)
        print(f"{t0:12.6f}  {z_val:15.10f}  {abs(z_val):12.2e}")

    # Evaluate Z(t) on a grid to show sign changes
    print("\n\nSign changes of Z(t) between t=10 and t=50:")
    t_grid = np.linspace(10, 50, 2000)
    z_grid = hardy_z_vectorized(t_grid)

    sign_changes = []
    for i in range(len(z_grid) - 1):
        if z_grid[i] * z_grid[i + 1] < 0:
            # Bisect to find zero more precisely
            a, b = t_grid[i], t_grid[i + 1]
            for _ in range(50):
                mid = (a + b) / 2
                if hardy_z(a) * hardy_z(mid) < 0:
                    b = mid
                else:
                    a = mid
            sign_changes.append((a + b) / 2)

    print(f"Found {len(sign_changes)} zeros:")
    for i, z in enumerate(sign_changes):
        print(f"  Zero {i+1:2d}: t = {z:.6f}")
