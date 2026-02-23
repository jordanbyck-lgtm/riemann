#!/usr/bin/env python3
"""
Riemann Hypothesis: Computational Demonstrations

This script runs all the major computational methods and produces
results demonstrating the strongest evidence for the Riemann Hypothesis.

It covers:
1. Riemann-Siegel computation and zero verification
2. Turing's method for rigorous verification
3. GUE random matrix statistics (Montgomery-Odlyzko)
4. Keating-Snaith moment conjectures
5. De Bruijn-Newman constant

Run: python demo.py [--plot] [--range T_MAX]
"""

import sys
import time
import numpy as np
from numpy import pi, log

# Local modules
from riemann_siegel import hardy_z, hardy_z_vectorized, riemann_siegel_theta
from zero_counting import N_smooth, gram_points, turing_method_verify, count_sign_changes
from zero_finding import (find_zeros_in_range, newton_zero,
                          lehmer_phenomenon_search, odlyzko_schonhage_idea)
from pair_correlation import (gue_pair_correlation, gue_nearest_neighbor_spacing,
                              poisson_spacing, normalize_zeros,
                              pair_correlation_empirical,
                              spacing_distribution_empirical,
                              number_variance_empirical, gue_number_variance)
from moments import moment_constant, numerical_moment, rmt_factor, arithmetic_factor
from newman_constant import (phi_function, H_t, find_H_t_zeros,
                             zero_dynamics_ode, compute_lambda_upper_bound)


def separator(title):
    """Print a section separator."""
    width = 70
    print("\n" + "=" * width)
    print(f"  {title}")
    print("=" * width)


def demo_riemann_siegel():
    """Demonstrate the Riemann-Siegel formula."""
    separator("1. RIEMANN-SIEGEL FORMULA")
    print("""
The Riemann-Siegel formula (1932) computes Z(t), the Hardy Z-function,
whose zeros on the real line correspond to zeros of zeta(1/2 + it).

Z(t) = 2 * sum_{n=1}^{N} cos(theta(t) - t*log(n))/sqrt(n) + R(t)
where N = floor(sqrt(t/(2*pi))).
""")

    # Known zeros
    known_zeros = [
        14.134725142, 21.022039639, 25.010857580, 30.424876126,
        32.935061588, 37.586178159, 40.918719012, 43.327073281,
        48.005150881, 49.773832478,
    ]

    print("Verification at the first 10 known zeros:")
    print(f"{'#':>3s}  {'t (known)':>15s}  {'Z(t)':>15s}  {'|Z(t)|':>12s}")
    print("-" * 50)
    for i, t0 in enumerate(known_zeros):
        z_val = hardy_z(t0)
        print(f"{i+1:3d}  {t0:15.9f}  {z_val:15.10f}  {abs(z_val):12.2e}")

    # Theta function values
    print("\nRiemann-Siegel theta function:")
    for t in [10, 50, 100, 500, 1000]:
        print(f"  theta({t}) = {riemann_siegel_theta(t):.8f}")


def demo_zero_counting():
    """Demonstrate zero counting and Turing's method."""
    separator("2. ZERO COUNTING & TURING'S METHOD")
    print("""
Turing's method (1953): Verify ALL zeros are on the critical line by
comparing sign-change count (zeros on Re=1/2 only) with N(T) formula
(total zeros from the argument principle). If they agree, RH holds
in that range.
""")

    # Count zeros in progressively larger ranges
    ranges = [(10, 100), (100, 500), (500, 1000)]
    for t_start, t_end in ranges:
        result = turing_method_verify(t_start, t_end)
        status = "PASS" if result['all_on_critical_line'] else "FAIL"
        print(f"  [{t_start:6.0f}, {t_end:6.0f}]: "
              f"found={result['zeros_found']:4d}, "
              f"expected={result['expected_int']:4d}, "
              f"status={status}")

    # Gram points
    print("\nGram points and Gram's law (first 15):")
    gpts = gram_points(0, 15)
    violations = 0
    for n, g in enumerate(gpts):
        z_val = hardy_z(g)
        holds = (-1) ** n * z_val > 0
        if not holds:
            violations += 1
        mark = " " if holds else "*"
        print(f"  g_{n:2d} = {g:10.6f}  Z(g_n) = {z_val:10.6f}  "
              f"Gram law: {'YES' if holds else 'NO ':>3s} {mark}")
    print(f"\n  Gram law violations: {violations}/{len(gpts)} "
          f"({100*violations/len(gpts):.0f}%, "
          f"expected ~27% asymptotically)")


def demo_pair_correlation():
    """Demonstrate Montgomery's pair correlation and GUE statistics."""
    separator("3. MONTGOMERY'S PAIR CORRELATION & GUE")
    print("""
Montgomery (1973) conjectured that normalized zero spacings of zeta
follow GUE (Gaussian Unitary Ensemble) statistics. Dyson immediately
recognized this as the eigenvalue statistics of random Hermitian matrices.

Odlyzko (1987) confirmed this spectacularly by computing millions of
zeros near t ~ 10^20.

Key signatures:
- Level repulsion: p(0) = 0 (zeros avoid each other)
- Pair correlation: R_2(x) = 1 - (sin(pi*x)/(pi*x))^2
- Logarithmic number variance (vs linear for Poisson)
""")

    # Compute a batch of zeros
    print("Computing zeros in [10, 2000]...")
    t0 = time.time()
    result = count_sign_changes(10, 2000, num_points=100000)
    zeros = np.array(result['zeros'])
    elapsed = time.time() - t0
    print(f"  Found {len(zeros)} zeros in {elapsed:.1f}s")

    # Spacing distribution
    spacing_result = spacing_distribution_empirical(zeros, num_bins=50)
    print(f"\nNearest-neighbor spacing distribution:")
    print(f"  Mean normalized spacing: {spacing_result['mean_spacing']:.4f} (should be ~1)")
    print(f"  Variance: {spacing_result['variance']:.4f}")
    gue_var = 1 - 3 / pi + 4 / (pi * pi) * (pi / 2 - 1)  # approximate
    print(f"  GUE Wigner surmise variance: ~0.178")

    print(f"\n  {'s':>6s}  {'Empirical':>10s}  {'GUE':>10s}  {'Poisson':>10s}")
    print("  " + "-" * 42)
    step = max(1, len(spacing_result['s']) // 15)
    for i in range(0, len(spacing_result['s']), step):
        print(f"  {spacing_result['s'][i]:6.3f}  "
              f"{spacing_result['p_empirical'][i]:10.4f}  "
              f"{spacing_result['p_gue'][i]:10.4f}  "
              f"{spacing_result['p_poisson'][i]:10.4f}")

    # Pair correlation
    pc_result = pair_correlation_empirical(zeros, num_bins=40)
    print(f"\nPair correlation R_2(alpha):")
    print(f"  (GUE: 1-(sin(pi*x)/(pi*x))^2, Poisson: 1)")
    print(f"  {'alpha':>8s}  {'Empirical':>10s}  {'GUE':>10s}")
    print("  " + "-" * 30)
    step = max(1, len(pc_result['alpha']) // 12)
    for i in range(0, len(pc_result['alpha']), step):
        print(f"  {pc_result['alpha'][i]:8.3f}  "
              f"{pc_result['R2_empirical'][i]:10.4f}  "
              f"{pc_result['R2_gue'][i]:10.4f}")


def demo_moments():
    """Demonstrate Keating-Snaith moment conjectures."""
    separator("4. KEATING-SNAITH MOMENT CONJECTURES")
    print("""
Keating & Snaith (2000) used Random Matrix Theory to predict:

  (1/T) * int_0^T |zeta(1/2+it)|^{2k} dt ~ C_k * (log T)^{k^2}

where C_k = g_k * a_k with g_k from RMT and a_k from arithmetic.

Known exact values:
  C_1 = 1          (Hardy-Littlewood, 1918)
  C_2 = 1/(2pi^2)  (Ingham, 1926)
""")

    print("Keating-Snaith constants C_k = g_k * a_k:")
    print(f"{'k':>5s}  {'g_k':>14s}  {'a_k':>14s}  {'C_k':>14s}  {'(logT)^':>7s}")
    print("-" * 62)
    for k in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0]:
        mc = moment_constant(k, num_primes=200)
        print(f"{k:5.1f}  {mc['g_k']:14.6e}  {mc['a_k']:14.6e}  "
              f"{mc['C_k']:14.6e}  k^2={mc['power']:.0f}")

    print(f"\nComparison with known values:")
    mc1 = moment_constant(1.0, num_primes=500)
    print(f"  C_1 computed: {mc1['C_k']:.8f}  exact: 1.00000000")
    mc2 = moment_constant(2.0, num_primes=500)
    print(f"  C_2 computed: {mc2['C_k']:.8f}  exact: {1/(2*pi**2):.8f}")

    # Numerical verification
    print("\nNumerical moment verification (comparing integral vs prediction):")
    for k_val in [1.0, 2.0]:
        for T in [200, 500, 1000]:
            result = numerical_moment(k_val, T, num_points=5000)
            print(f"  k={k_val:.0f}, T={T:5d}: "
                  f"integral={result['moment']:10.4f}, "
                  f"prediction={result['predicted']:10.4f}, "
                  f"ratio={result['ratio']:6.3f}")


def demo_newman_constant():
    """Demonstrate the de Bruijn-Newman constant."""
    separator("5. DE BRUIJN-NEWMAN CONSTANT (RODGERS-TAO)")
    print("""
The de Bruijn-Newman constant Lambda is defined by:
  H_t has only real zeros iff t >= Lambda

Key result: Rodgers & Tao (2020) proved Lambda >= 0.
Combined with Platt-Trudgian (2021): 0 <= Lambda <= 0.22.
The Riemann Hypothesis is equivalent to Lambda = 0.

Rodgers-Tao approach: zeros of H_t satisfy Dyson Brownian motion
(electrostatic repulsion). If H_0 had complex zeros, they would
persist for t < 0, contradicting de Bruijn's result.
""")

    # Phi function
    print("Phi function (kernel):")
    for u in [0.0, 0.5, 1.0, 1.5, 2.0]:
        val = phi_function(np.array([u]))[0]
        print(f"  Phi({u:.1f}) = {val:15.8e}")

    # H_t zeros
    print("\nReal zeros of H_t for various t:")
    for t_val in [0.0, 0.1, 0.5]:
        zeros = find_H_t_zeros(t_val, z_max=25, num_search=200)
        print(f"  t={t_val:.1f}: {len(zeros)} zeros", end="")
        if zeros:
            print(f" (first: {zeros[0]:.4f})", end="")
        print()

    # Upper bound
    print("\nCrude Lambda upper bound from computed zeros:")
    bound = compute_lambda_upper_bound(T_max=300)
    print(f"  Zeros used: {bound['num_zeros']}")
    print(f"  Min normalized gap: {bound['min_normalized_gap']:.6f}")
    print(f"  Crude bound: Lambda <= {bound['upper_bound']:.2f}")
    print(f"  Best known: 0 <= Lambda <= 0.22")

    # Zero dynamics
    print("\nZero dynamics (Dyson Brownian motion):")
    init = np.array([5.0, 10.0, 15.0, 20.0, 25.0])
    dyn = zero_dynamics_ode(init, (0.0, 0.05), dt=0.002)
    print(f"  t=0.00: {dyn['zeros'][0]}")
    mid = len(dyn['t']) // 2
    print(f"  t={dyn['t'][mid]:.2f}: {dyn['zeros'][mid]}")
    print(f"  t={dyn['t'][-1]:.2f}: {dyn['zeros'][-1]}")
    print("  (Zeros repel: gaps widen over time)")


def demo_lehmer():
    """Search for the Lehmer phenomenon."""
    separator("6. LEHMER PHENOMENON")
    print("""
Lehmer (1956) discovered a place near t ~ 7005 where Z(t) comes
extremely close to zero without changing sign, indicating two very
closely spaced zeros. This is significant because:
1. It challenges numerical verification methods
2. It tests the limits of Turing's method
3. Related to moments and extreme values of zeta
""")

    print("Searching for near-misses in [6900, 7100]...")
    near_misses = lehmer_phenomenon_search(6900, 7100, threshold=0.5)
    print(f"Found {len(near_misses)} near-misses:")
    for nm in near_misses[:10]:
        print(f"  t = {nm['t']:.8f}  Z(t) = {nm['Z_value']:12.8f}")

    # The famous Lehmer zero pair near 7005
    print("\nDetailed view near t = 7005 (Lehmer's discovery):")
    t_fine = np.linspace(7004, 7006, 500)
    z_fine = hardy_z_vectorized(t_fine)
    min_idx = np.argmin(np.abs(z_fine))
    print(f"  Minimum |Z(t)| at t = {t_fine[min_idx]:.6f}: "
          f"Z = {z_fine[min_idx]:.10f}")


def demo_summary():
    """Print a summary of the evidence."""
    separator("SUMMARY: EVIDENCE FOR THE RIEMANN HYPOTHESIS")
    print("""
The methods implemented here represent the major computational
approaches to the Riemann Hypothesis:

VERIFIED (proven results):
  * All zeros up to height T ~ 3*10^12 lie on the critical line
    (Platt & Trudgian, 2021) -- over 10^13 zeros verified
  * The de Bruijn-Newman constant Lambda >= 0 (Rodgers & Tao, 2020),
    meaning RH is "just barely true" if true at all
  * 0 <= Lambda <= 0.22 (Platt & Trudgian, 2021)

STRONG NUMERICAL EVIDENCE:
  * Zero statistics match GUE random matrix predictions to extraordinary
    precision (Odlyzko, zeros near t ~ 10^20)
  * Keating-Snaith moment conjectures consistent with all computed data
  * No Lehmer-type near-misses have led to actual counterexamples

OPEN:
  * A proof remains one of the great challenges of mathematics
  * Worth $1,000,000 (Clay Millennium Prize)
  * Would have profound consequences for prime number distribution

The RH connects: number theory <-> random matrices <-> quantum chaos
                  <-> nuclear physics <-> quantum gravity
""")


def main():
    """Run all demonstrations."""
    plot_mode = "--plot" in sys.argv

    print("RIEMANN HYPOTHESIS: COMPUTATIONAL METHODS")
    print("Implementations of the most significant advances")
    print(f"{'=' * 70}")

    t_start = time.time()

    demo_riemann_siegel()
    demo_zero_counting()
    demo_pair_correlation()
    demo_moments()
    demo_newman_constant()
    demo_lehmer()
    demo_summary()

    elapsed = time.time() - t_start
    print(f"\nTotal computation time: {elapsed:.1f}s")

    if plot_mode:
        create_plots()


def create_plots():
    """Generate matplotlib visualizations (optional)."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not available, skipping plots.")
        return

    print("\nGenerating plots...")

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # Plot 1: Z(t) and its zeros
    ax = axes[0, 0]
    t = np.linspace(10, 60, 2000)
    z = hardy_z_vectorized(t)
    ax.plot(t, z, 'b-', linewidth=0.8)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.set_title('Hardy Z-function Z(t)')
    ax.set_xlabel('t')
    ax.set_ylabel('Z(t)')

    # Plot 2: Spacing distribution
    ax = axes[0, 1]
    result = count_sign_changes(10, 2000, num_points=100000)
    zeros = np.array(result['zeros'])
    spacings = normalize_zeros(zeros)
    s = np.linspace(0, 3, 100)
    ax.hist(spacings, bins=40, density=True, alpha=0.6, label='Zeta zeros')
    ax.plot(s, gue_nearest_neighbor_spacing(s), 'r-', linewidth=2, label='GUE')
    ax.plot(s, poisson_spacing(s), 'g--', linewidth=2, label='Poisson')
    ax.set_title("Spacing Distribution\n(Montgomery-Odlyzko)")
    ax.set_xlabel('Normalized spacing s')
    ax.legend()

    # Plot 3: Pair correlation
    ax = axes[0, 2]
    pc = pair_correlation_empirical(zeros, num_bins=30)
    ax.bar(pc['alpha'], pc['R2_empirical'], width=pc['alpha'][1]-pc['alpha'][0],
           alpha=0.6, label='Zeta zeros')
    ax.plot(pc['alpha'], pc['R2_gue'], 'r-', linewidth=2, label='GUE')
    ax.axhline(y=1, color='g', linestyle='--', label='Poisson')
    ax.set_title("Pair Correlation R₂(α)")
    ax.set_xlabel('α')
    ax.legend()

    # Plot 4: Number variance
    ax = axes[1, 0]
    nv = number_variance_empirical(zeros)
    L = nv['L']
    ax.plot(L, nv['sigma2_empirical'], 'b.', markersize=6, label='Empirical')
    ax.plot(L, nv['sigma2_gue'], 'r-', linewidth=2, label='GUE')
    ax.plot(L, nv['sigma2_poisson'], 'g--', linewidth=2, label='Poisson')
    ax.set_title('Number Variance Σ²(L)')
    ax.set_xlabel('L')
    ax.legend()

    # Plot 5: Phi function (Newman constant)
    ax = axes[1, 1]
    u = np.linspace(0, 3, 200)
    phi = phi_function(u)
    ax.plot(u, phi, 'b-', linewidth=1.5)
    ax.set_title('Φ(u) kernel\n(de Bruijn-Newman)')
    ax.set_xlabel('u')
    ax.set_yscale('symlog')

    # Plot 6: Moment constants
    ax = axes[1, 2]
    k_vals = np.linspace(0.1, 4, 30)
    g_vals = [rmt_factor(k) for k in k_vals]
    a_vals = [arithmetic_factor(k, 100) for k in k_vals]
    c_vals = [g * a for g, a in zip(g_vals, a_vals)]
    ax.semilogy(k_vals, c_vals, 'b-', linewidth=2)
    ax.set_title('Keating-Snaith C_k')
    ax.set_xlabel('k')
    ax.set_ylabel('C_k')

    plt.tight_layout()
    plt.savefig('riemann_plots.png', dpi=150)
    print("Saved: riemann_plots.png")


if __name__ == "__main__":
    main()
