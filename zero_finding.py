"""
High-precision zero finding for the Riemann zeta function.

This module implements methods for locating zeros of zeta on the
critical line with high precision, including:

1. Bracketing via Gram points and sign changes of Z(t)
2. Newton's method using the derivative Z'(t)
3. The Odlyzko-Schönhage algorithm idea (FFT-based batch computation)
4. Lehmer phenomenon detection (near-misses where Z(t) almost touches zero)

The Odlyzko-Schönhage algorithm (1988) is the fastest known method for
computing many zeros of zeta. It uses the Fast Fourier Transform to
evaluate the Riemann-Siegel sum at O(T^{1/2+epsilon}) points simultaneously,
achieving amortized cost O(T^epsilon) per zero (vs O(T^{1/2}) for direct
Riemann-Siegel at each point).

Odlyzko used this algorithm to compute 10^9 zeros near t ~ 10^20,
providing the strongest numerical evidence for both:
- The Riemann Hypothesis
- Montgomery's pair correlation conjecture / GUE statistics

Reference:
    Odlyzko, A.M. & Schönhage, A. (1988). "Fast algorithms for
        multiple evaluations of the Riemann zeta function"
    Odlyzko, A.M. (1992). "The 10^20-th zero of the Riemann zeta function
        and 175 million of its neighbors"
    Gourdon, X. (2004). "The 10^13 first zeros of the Riemann zeta function,
        and zeros computation at very large height"
"""

import numpy as np
from numpy import pi, log, sqrt, cos, sin
from riemann_siegel import hardy_z, riemann_siegel_theta, hardy_z_vectorized


def find_zeros_in_range(t_start, t_end, resolution=0.1):
    """
    Find all zeros of Z(t) in [t_start, t_end] using adaptive search.

    Uses a two-phase approach:
    1. Coarse grid to detect sign changes
    2. Bisection to refine each zero to high precision

    Parameters
    ----------
    t_start : float
        Start of range.
    t_end : float
        End of range.
    resolution : float
        Initial grid spacing (should be < minimum expected zero gap).

    Returns
    -------
    list of dict
        Each entry has 't': zero location, 'Z_value': Z(t) residual.
    """
    # Phase 1: Coarse grid
    num_points = int((t_end - t_start) / resolution) + 1
    t_grid = np.linspace(t_start, t_end, num_points)
    z_grid = hardy_z_vectorized(t_grid)

    zeros = []

    # Phase 2: Refine sign changes
    for i in range(len(z_grid) - 1):
        if z_grid[i] * z_grid[i + 1] < 0:
            t0 = _bisect_zero(t_grid[i], t_grid[i + 1])
            zeros.append({
                't': t0,
                'Z_value': hardy_z(t0),
            })
        elif abs(z_grid[i]) < 1e-6:
            # Potential double zero or near-touch (Lehmer phenomenon)
            t0 = _refine_near_zero(t_grid[i])
            if t0 is not None:
                zeros.append({
                    't': t0,
                    'Z_value': hardy_z(t0),
                    'lehmer': True,
                })

    return zeros


def _bisect_zero(a, b, tol=1e-12, max_iter=80):
    """Bisection to find zero of Z(t) between a and b."""
    za = hardy_z(a)
    for _ in range(max_iter):
        mid = (a + b) / 2
        zm = hardy_z(mid)
        if abs(b - a) < tol:
            break
        if za * zm <= 0:
            b = mid
        else:
            a = mid
            za = zm
    return (a + b) / 2


def _refine_near_zero(t0, search_radius=0.5, num_points=100):
    """Refine a potential near-zero using fine grid search."""
    t_fine = np.linspace(t0 - search_radius, t0 + search_radius, num_points)
    z_fine = hardy_z_vectorized(t_fine)

    # Look for actual sign changes
    for i in range(len(z_fine) - 1):
        if z_fine[i] * z_fine[i + 1] < 0:
            return _bisect_zero(t_fine[i], t_fine[i + 1])

    # Check for minimum near zero
    min_idx = np.argmin(np.abs(z_fine))
    if abs(z_fine[min_idx]) < 1e-8:
        return t_fine[min_idx]

    return None


def newton_zero(t0, max_iter=20, tol=1e-14):
    """
    Newton's method to refine a zero of Z(t).

    Uses numerical differentiation for Z'(t).

    Parameters
    ----------
    t0 : float
        Initial guess (should be close to a zero).
    max_iter : int
        Maximum iterations.
    tol : float
        Convergence tolerance.

    Returns
    -------
    float
        Refined zero location.
    """
    t = t0
    h = 1e-8

    for _ in range(max_iter):
        z = hardy_z(t)
        if abs(z) < tol:
            break

        # Numerical derivative
        z_plus = hardy_z(t + h)
        z_minus = hardy_z(t - h)
        z_prime = (z_plus - z_minus) / (2 * h)

        if abs(z_prime) < 1e-15:
            break

        t -= z / z_prime

    return t


def lehmer_phenomenon_search(t_start, t_end, threshold=0.1):
    """
    Search for the Lehmer phenomenon: places where Z(t) comes very
    close to zero without changing sign.

    The Lehmer phenomenon is important because:
    1. It indicates closely spaced pairs of zeros
    2. It makes rigorous zero verification more difficult
    3. It tests the limits of the Turing method

    The most famous example is near t ≈ 7005 (Lehmer's original discovery,
    1956), where Z(t) has a very small local extremum.

    Parameters
    ----------
    t_start : float
        Start of search.
    t_end : float
        End of search.
    threshold : float
        How close to zero counts as a "near miss."

    Returns
    -------
    list of dict
        Near-miss locations and values.
    """
    resolution = 0.05
    num_points = int((t_end - t_start) / resolution) + 1
    t_grid = np.linspace(t_start, t_end, num_points)
    z_grid = hardy_z_vectorized(t_grid)

    near_misses = []

    for i in range(1, len(z_grid) - 1):
        # Local extremum (no sign change, but close to zero)
        if (z_grid[i - 1] * z_grid[i + 1] > 0 and
                abs(z_grid[i]) < threshold and
                abs(z_grid[i]) < abs(z_grid[i - 1]) and
                abs(z_grid[i]) < abs(z_grid[i + 1])):

            # Refine the extremum
            a, b = t_grid[i - 1], t_grid[i + 1]
            # Golden section search for minimum of |Z(t)|
            gr = (sqrt(5) + 1) / 2
            for _ in range(50):
                c = b - (b - a) / gr
                d = a + (b - a) / gr
                if abs(hardy_z(c)) < abs(hardy_z(d)):
                    b = d
                else:
                    a = c
            t_ext = (a + b) / 2
            z_ext = hardy_z(t_ext)

            near_misses.append({
                't': t_ext,
                'Z_value': z_ext,
                'is_near_miss': abs(z_ext) < threshold,
            })

    return near_misses


def odlyzko_schonhage_idea(t_center, num_points, bandwidth):
    """
    Simplified demonstration of the Odlyzko-Schönhage algorithm concept.

    The full algorithm evaluates the Riemann-Siegel sum at many equally-spaced
    points simultaneously using the FFT. The key identity is:

    For the main sum S(t) = sum_{n=1}^{N} n^{-1/2-it}:

    Evaluate at t_k = t_0 + k*delta for k = 0, ..., K-1:
        S(t_k) = sum_{n=1}^{N} n^{-1/2-i*t_0} * n^{-ik*delta}

    This can be rearranged into a form computable by FFT when the
    n^{-ik*delta} factors are discretized appropriately.

    The actual algorithm involves:
    1. Splitting the sum into O(sqrt(T)) blocks
    2. Taylor expanding within each block
    3. Using FFT to evaluate the polynomial sums
    4. Combining results

    Achieves O(T^{1/2+epsilon}) total work for O(T^{1/2}) evaluations,
    i.e., O(T^epsilon) amortized per evaluation.

    This implementation demonstrates the concept with a direct vectorized
    computation (not the full FFT algorithm).

    Parameters
    ----------
    t_center : float
        Center of evaluation range.
    num_points : int
        Number of equally-spaced evaluation points.
    bandwidth : float
        Half-width of the evaluation range.

    Returns
    -------
    dict
        't_values': evaluation points,
        'Z_values': Z(t) at those points,
        'zeros': detected zero locations.
    """
    t_values = np.linspace(t_center - bandwidth, t_center + bandwidth, num_points)

    N = int(np.floor(np.sqrt(t_center / (2 * pi))))

    # Precompute log(n) and 1/sqrt(n) for all n
    n_vals = np.arange(1, N + 1, dtype=np.float64)
    log_n = np.log(n_vals)
    inv_sqrt_n = 1.0 / np.sqrt(n_vals)

    # Vectorized computation of main sum at all t values
    theta = riemann_siegel_theta(t_values)  # vectorized

    # Z(t) = 2 * sum_{n=1}^{N} cos(theta(t) - t*log(n)) / sqrt(n)
    # Shape: (num_points, N)
    t_col = t_values[:, np.newaxis]  # (K, 1)
    log_n_row = log_n[np.newaxis, :]  # (1, N)
    theta_col = theta[:, np.newaxis]  # (K, 1)

    phase = theta_col - t_col * log_n_row  # (K, N)
    terms = np.cos(phase) * inv_sqrt_n[np.newaxis, :]  # (K, N)
    Z_values = 2 * np.sum(terms, axis=1)  # (K,)

    # Find zeros
    zeros = []
    for i in range(len(Z_values) - 1):
        if Z_values[i] * Z_values[i + 1] < 0:
            t0 = _bisect_zero(t_values[i], t_values[i + 1])
            zeros.append(t0)

    return {
        't_values': t_values,
        'Z_values': Z_values,
        'zeros': zeros,
        'N_main_sum': N,
    }


def compute_zeros_batch(t_start, t_end, batch_size=10000):
    """
    Compute zeros in batches using the vectorized approach.

    Parameters
    ----------
    t_start : float
        Start of range.
    t_end : float
        End of range.
    batch_size : int
        Points per batch.

    Returns
    -------
    list of float
        Zero locations.
    """
    all_zeros = []

    # Determine appropriate resolution
    # Mean zero spacing near T is 2*pi/log(T/(2*pi))
    mean_spacing = 2 * pi / log((t_start + t_end) / 2 / (2 * pi))
    resolution = mean_spacing / 4  # At least 4 points per mean spacing

    t = t_start
    while t < t_end:
        t_batch_end = min(t + batch_size * resolution, t_end)
        bandwidth = (t_batch_end - t) / 2
        t_center = t + bandwidth

        result = odlyzko_schonhage_idea(t_center, batch_size, bandwidth)
        all_zeros.extend(result['zeros'])

        t = t_batch_end

    return all_zeros


if __name__ == "__main__":
    print("High-Precision Zero Finding")
    print("=" * 60)

    # 1. Find zeros in a range
    print("\nFinding zeros in [10, 100]:")
    zeros = find_zeros_in_range(10, 100, resolution=0.05)
    print(f"Found {len(zeros)} zeros:")
    for i, z in enumerate(zeros):
        lehmer = " [LEHMER]" if z.get('lehmer') else ""
        print(f"  #{i+1:3d}: t = {z['t']:15.10f}  Z(t) = {z['Z_value']:12.2e}{lehmer}")

    # 2. Newton refinement
    print("\nNewton refinement of first 5 zeros:")
    for z in zeros[:5]:
        t_refined = newton_zero(z['t'])
        z_val = hardy_z(t_refined)
        print(f"  t = {t_refined:.15f}  Z(t) = {z_val:.2e}")

    # 3. Lehmer phenomenon search
    print("\nSearching for Lehmer phenomenon in [6900, 7100]:")
    near_misses = lehmer_phenomenon_search(6900, 7100, threshold=0.5)
    print(f"Found {len(near_misses)} near-misses:")
    for nm in near_misses[:10]:
        print(f"  t = {nm['t']:.6f}  Z(t) = {nm['Z_value']:.8f}")

    # 4. Batch computation demo
    print("\nOdlyzko-Schönhage style batch computation near t=10000:")
    result = odlyzko_schonhage_idea(10000, num_points=50000, bandwidth=50)
    print(f"  Evaluated Z(t) at {len(result['t_values'])} points")
    print(f"  Main sum truncation: N = {result['N_main_sum']}")
    print(f"  Found {len(result['zeros'])} zeros in [9950, 10050]")
    if result['zeros']:
        for z in result['zeros'][:5]:
            print(f"    t = {z:.10f}")
