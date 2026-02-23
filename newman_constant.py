"""
The De Bruijn-Newman Constant and the Rodgers-Tao Theorem.

The de Bruijn-Newman constant Λ is a real number defined through a
family of entire functions H_t(z) that deform the Riemann xi function.

Define:
    Phi(u) = sum_{n=1}^{inf} (2*pi^2*n^4*exp(9u) - 3*pi*n^2*exp(5u))
             * exp(-pi*n^2*exp(4u))

This is a super-exponentially decaying function related to the Jacobi
theta function. Then define:

    H_t(z) = int_0^{inf} exp(tu^2) * Phi(u) * cos(zu) du

Key facts:
- H_0(z) = (1/8) * xi(1/2 + iz/2)  (proportional to the Riemann xi function)
- The Riemann Hypothesis is equivalent to: all zeros of H_0 are real
- De Bruijn (1950) proved: H_t has only real zeros for t >= 1/2
- Newman (1976) proved: there exists Λ such that H_t has only real zeros
  iff t >= Λ, and conjectured Λ >= 0
- The Riemann Hypothesis is equivalent to Λ <= 0

Timeline of bounds on Λ:
    Λ <= 1/2     de Bruijn (1950)
    Λ >= -50     early bounds
    Λ >= -1.1e-11  Ki, Kim, Lee (2009)
    Λ >= 0       Rodgers & Tao (2018) [proved Newman's conjecture!]
    Λ <= 0.22    Platt & Trudgian (2021)

So currently: 0 <= Λ <= 0.22

The RH is equivalent to Λ = 0. Rodgers-Tao's proof that Λ >= 0 means
the RH, if true, is "just barely" true.

Rodgers-Tao approach (simplified):
- They studied how zeros of H_t evolve as t increases from 0
- Zeros satisfy an ODE (the "electrostatic" or "Dyson Brownian motion" equation)
- They showed that if H_0 had non-real zeros, they would persist for t slightly
  below 0, contradicting Newman's result. So Λ >= 0.

This module implements:
1. Numerical computation of H_t(z)
2. Visualization of zero evolution under the heat flow
3. The backward heat equation and zero dynamics
4. Upper bound verification methods

Reference:
    de Bruijn, N.G. (1950). "The roots of trigonometric integrals"
    Newman, C.M. (1976). "Fourier transforms with only real zeros"
    Rodgers, B. & Tao, T. (2020). "The de Bruijn-Newman constant is non-negative"
        Annals of Mathematics, 192(2), 599-642.
    Platt, D.J. & Trudgian, T.S. (2021). "The Riemann hypothesis is true
        up to 3·10^12"
"""

import numpy as np
from numpy import pi, exp, sqrt, log, cos, sin
from scipy import integrate


def phi_function(u, num_terms=30):
    """
    Compute the Phi function (kernel of the xi integral transform).

    Phi(u) = sum_{n=1}^{inf} (2*pi^2*n^4*exp(9u) - 3*pi*n^2*exp(5u))
             * exp(-pi*n^2*exp(4u))

    This function decays super-exponentially for large u.

    Parameters
    ----------
    u : float or ndarray
        Argument (u >= 0).
    num_terms : int
        Number of terms in the sum.

    Returns
    -------
    float or ndarray
        Phi(u) value.
    """
    u = np.asarray(u, dtype=np.float64)
    result = np.zeros_like(u, dtype=np.float64)

    for n in range(1, num_terms + 1):
        n2 = n * n
        n4 = n2 * n2
        e4u = np.exp(4 * u)
        exponential = np.exp(-pi * n2 * e4u)

        term = (2 * pi ** 2 * n4 * np.exp(9 * u) -
                3 * pi * n2 * np.exp(5 * u)) * exponential

        result += term

        # Early termination when terms are negligible
        if np.all(np.abs(term) < 1e-30 * (np.abs(result) + 1e-100)):
            break

    return result


def H_t(z, t, u_max=5.0, num_quad=500):
    """
    Compute H_t(z) = int_0^{inf} exp(t*u^2) * Phi(u) * cos(z*u) du.

    H_0(z) is proportional to the Riemann xi function:
        H_0(z) = (1/8) * xi(1/2 + iz/2)

    For t > 0, H_t is a "smoothed" version of xi (heat equation flow).
    For t = Λ, H_t first develops only real zeros.

    Parameters
    ----------
    z : float or complex
        Argument.
    t : float
        Heat flow parameter.
    u_max : float
        Upper limit of integration (Phi decays rapidly).
    num_quad : int
        Number of quadrature points.

    Returns
    -------
    complex
        H_t(z) value.
    """
    z = complex(z)

    def integrand_real(u):
        if u < 1e-15:
            return 0.0
        phi_val = float(phi_function(np.array([u]))[0])
        heat_kernel = exp(t * u * u)
        cos_part = np.real(np.cos(z * u))
        return heat_kernel * phi_val * cos_part

    def integrand_imag(u):
        if u < 1e-15:
            return 0.0
        phi_val = float(phi_function(np.array([u]))[0])
        heat_kernel = exp(t * u * u)
        sin_part = -np.imag(np.cos(z * u))
        return heat_kernel * phi_val * sin_part

    # Adaptive quadrature
    real_part, _ = integrate.quad(integrand_real, 0, u_max, limit=200)
    imag_part, _ = integrate.quad(integrand_imag, 0, u_max, limit=200)

    return complex(real_part, imag_part)


def find_H_t_zeros(t, z_max=50, num_search=500):
    """
    Find zeros of H_t on the real line by sign-change detection.

    For t >= Λ, all zeros should be real.
    For t < Λ, some zeros will be complex (off the real line).

    Parameters
    ----------
    t : float
        Heat flow parameter.
    z_max : float
        Search range [-z_max, z_max] (H_t is even, so only z > 0 needed).
    num_search : int
        Number of search points.

    Returns
    -------
    list of float
        Real zeros of H_t(z).
    """
    from scipy.optimize import brentq

    z_grid = np.linspace(0.1, z_max, num_search)
    h_vals = np.array([H_t(z, t).real for z in z_grid])

    zeros = []
    for i in range(len(h_vals) - 1):
        if h_vals[i] * h_vals[i + 1] < 0:
            try:
                z0 = brentq(lambda z: H_t(z, t).real, z_grid[i], z_grid[i + 1])
                zeros.append(z0)
            except (ValueError, RuntimeError):
                pass

    return zeros


def zero_dynamics_ode(zeros_init, t_span, dt=0.001):
    """
    Simulate the dynamics of zeros under the backward heat flow.

    As t varies, the zeros of H_t satisfy the ODE system
    (Dyson Brownian motion / electrostatic analogy):

        dx_j/dt = -2 * sum_{k != j} 1/(x_j - x_k)

    This is the "Coulomb gas" equation where zeros repel each other
    like charged particles on a line. The key insight of Rodgers-Tao
    is that this dynamics prevents complex zeros from forming as
    t decreases to 0.

    Parameters
    ----------
    zeros_init : array-like
        Initial zero positions at t = t_span[0].
    t_span : tuple
        (t_start, t_end) time interval.
    dt : float
        Time step.

    Returns
    -------
    dict
        't': time points,
        'zeros': zeros[i, j] = position of j-th zero at time t[i].
    """
    zeros = np.array(zeros_init, dtype=np.float64)
    N = len(zeros)

    t_values = [t_span[0]]
    zeros_history = [zeros.copy()]

    t = t_span[0]
    direction = 1 if t_span[1] > t_span[0] else -1

    while (direction > 0 and t < t_span[1]) or (direction < 0 and t > t_span[1]):
        # Compute the Coulomb force on each zero
        force = np.zeros(N)
        for j in range(N):
            for k in range(N):
                if k != j:
                    diff = zeros[j] - zeros[k]
                    if abs(diff) > 1e-10:
                        force[j] -= 2.0 / diff

        # Euler step (for illustration; RK4 would be more accurate)
        zeros = zeros + direction * dt * force
        t += direction * dt

        t_values.append(t)
        zeros_history.append(zeros.copy())

    return {
        't': np.array(t_values),
        'zeros': np.array(zeros_history),
    }


def verify_newman_conjecture_numerically(t_values=None, z_range=30):
    """
    Numerically verify that H_t has all real zeros for t >= 0 and
    check for complex zeros when t < 0.

    This demonstrates the Rodgers-Tao result: Λ >= 0.

    Parameters
    ----------
    t_values : array-like or None
        Values of t to check.
    z_range : float
        Range to search for zeros.

    Returns
    -------
    list of dict
        Results for each t value.
    """
    if t_values is None:
        t_values = [-0.1, -0.05, -0.01, 0.0, 0.01, 0.05, 0.1, 0.2, 0.5]

    results = []
    for t in t_values:
        # Find real zeros
        real_zeros = find_H_t_zeros(t, z_max=z_range, num_search=200)

        # Check for evidence of complex zeros by evaluating H_t
        # on a grid in the complex plane
        has_complex_evidence = False
        if t < 0:
            # Sample along lines Im(z) = epsilon
            for eps in [0.1, 0.5, 1.0]:
                z_test = np.linspace(1, z_range, 50) + 1j * eps
                h_vals = [abs(H_t(z, t)) for z in z_test]
                min_val = min(h_vals)
                if min_val < 0.01:
                    has_complex_evidence = True
                    break

        results.append({
            't': t,
            'real_zeros_count': len(real_zeros),
            'real_zeros': real_zeros[:10],  # First 10
            'complex_zero_evidence': has_complex_evidence,
        })

    return results


def compute_lambda_upper_bound(T_max=1000, num_zeros=50):
    """
    Compute an approximate upper bound for Λ using the method of
    Ki, Kim, and Lee.

    The idea: if gamma_1, gamma_2, ... are consecutive zeros of H_0
    (i.e., zeros of xi), then define:

        lambda_n = (1 / (gamma_{n+1} - gamma_n)^2) * (something involving
                   the gaps and derivatives)

    A simplified version: Λ <= 1/(2*D^2) where D is the minimum
    normalized gap between zeros.

    Platt & Trudgian (2021) proved Λ <= 0.22 using rigorous computation
    of zeros up to height 3*10^12.

    Parameters
    ----------
    T_max : float
        Height to search for zeros.
    num_zeros : int
        Expected number of zeros.

    Returns
    -------
    dict
        Upper bound estimate and supporting data.
    """
    from zero_counting import count_sign_changes

    # Find zeros
    result = count_sign_changes(10, T_max, num_points=max(20000, int(T_max * 50)))
    zeros = np.array(result['zeros'])

    if len(zeros) < 3:
        return {'upper_bound': float('inf'), 'num_zeros': len(zeros)}

    # Normalize spacings
    gaps = np.diff(zeros)
    midpoints = (zeros[:-1] + zeros[1:]) / 2
    density = np.log(midpoints / (2 * pi)) / (2 * pi)
    normalized_gaps = gaps * density

    # The lambda bound from normalized gaps
    # Simplified: Λ <= max over pairs of a functional involving inverse gaps
    min_normalized_gap = np.min(normalized_gaps)
    max_lambda_local = 1.0 / (2 * min_normalized_gap ** 2)

    # Lehmer-type bound: if all gaps are > threshold, then Λ <= bound
    mean_gap = np.mean(normalized_gaps)
    gap_variance = np.var(normalized_gaps)

    return {
        'upper_bound': max_lambda_local,
        'min_normalized_gap': min_normalized_gap,
        'mean_normalized_gap': mean_gap,
        'gap_variance': gap_variance,
        'num_zeros': len(zeros),
        'note': (f"Crude upper bound Λ <= {max_lambda_local:.4f}. "
                 f"Best known: 0 <= Λ <= 0.22 (Platt-Trudgian 2021)")
    }


if __name__ == "__main__":
    print("De Bruijn-Newman Constant & Rodgers-Tao Theorem")
    print("=" * 60)

    # 1. Show Phi function
    print("\nPhi function values:")
    u_vals = np.linspace(0, 3, 20)
    phi_vals = phi_function(u_vals)
    for u, p in zip(u_vals, phi_vals):
        print(f"  Phi({u:.2f}) = {p:15.8e}")

    # 2. Compute H_t zeros for various t
    print("\n\nFinding zeros of H_t for various t values:")
    print("(Demonstrating that zeros are real for t >= 0)")
    for t_val in [0.0, 0.1, 0.5]:
        zeros = find_H_t_zeros(t_val, z_max=30, num_search=300)
        print(f"\n  t = {t_val:.2f}: {len(zeros)} real zeros found")
        if zeros:
            for z in zeros[:5]:
                h_val = H_t(z, t_val)
                print(f"    z = {z:10.6f}, H_t(z) = {h_val.real:12.2e}")

    # 3. Compute upper bound
    print("\n\nComputing crude upper bound for Λ:")
    bound_result = compute_lambda_upper_bound(T_max=200)
    print(f"  Zeros used: {bound_result['num_zeros']}")
    print(f"  Min normalized gap: {bound_result['min_normalized_gap']:.6f}")
    print(f"  Crude bound: Λ <= {bound_result['upper_bound']:.4f}")
    print(f"  {bound_result['note']}")

    # 4. Zero dynamics
    print("\n\nSimulating zero dynamics (Dyson Brownian motion):")
    init_zeros = np.array([5.0, 10.0, 15.0, 20.0, 25.0])
    dynamics = zero_dynamics_ode(init_zeros, (0.0, 0.1), dt=0.005)
    print(f"  Initial zeros: {init_zeros}")
    print(f"  Final zeros (t=0.1): {dynamics['zeros'][-1]}")
    print("  (Zeros repel each other via the Coulomb interaction)")
