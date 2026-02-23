# Riemann Hypothesis: Computational Methods

Implementations of the most significant computational methods and modern
advances related to the Riemann Hypothesis, covering work from Riemann (1859)
through Guth-Maynard (2024) and Baluyot et al. (2025).

## Modules

| Module | Description |
|--------|-------------|
| `riemann_siegel.py` | Riemann-Siegel formula for Z(t) on the critical line |
| `zero_finding.py` | High-precision zero location, Lehmer phenomenon, Odlyzko-Schönhage batch computation |
| `zero_counting.py` | Turing's method, Backlund's N(T) zero-counting formula, Gram points |
| `pair_correlation.py` | Montgomery's pair correlation conjecture & GUE comparison (Odlyzko) |
| `newman_constant.py` | De Bruijn-Newman constant Λ, heat flow H_t, Dyson Brownian motion (Rodgers-Tao 2020) |
| `moments.py` | Keating-Snaith moment conjectures from random matrix theory |
| `zero_density.py` | Zero density estimates: Guth-Maynard (2024) breakthrough + Baluyot et al. (2025) |
| `demo.py` | Full demonstration of all methods with optional matplotlib visualization |

## Key Results Demonstrated

- **Turing verification**: All 29 zeros in [10,100] confirmed on the critical line
- **GUE statistics**: Spacing distribution matches random matrix theory (level repulsion)
- **Keating-Snaith**: C_1 = 1.0000, C_2 = 1/(2π²) reproduced from Barnes G-function + Euler product
- **Lehmer phenomenon**: Near-miss detected at t ≈ 7005 (Lehmer's 1956 discovery)
- **Guth-Maynard**: First improvement to Ingham's 1940 zero-density bound in 84 years

## Key References

- Riemann, B. (1859). "Über die Anzahl der Primzahlen unter einer gegebenen Größe"
- Siegel, C.L. (1932). "Über Riemanns Nachlaß zur analytischen Zahlentheorie"
- Turing, A.M. (1953). "Some calculations of the Riemann zeta-function"
- Montgomery, H.L. (1973). "The pair correlation of zeros of the zeta function"
- Odlyzko, A.M. (1987). "On the distribution of spacings between zeros of the zeta function"
- Keating, J.P. & Snaith, N.C. (2000). "Random matrix theory and ζ(1/2+it)"
- Rodgers, B. & Tao, T. (2020). "The de Bruijn-Newman constant is non-negative"
- Platt, D.J. & Trudgian, T.S. (2021). "The Riemann hypothesis is true up to 3·10^12"
- Guth, L. & Maynard, J. (2024). "New large value estimates for Dirichlet polynomials"
- Baluyot, S. et al. (2025). "Pair correlation and the distribution of zeros"

## Requirements

```
pip install numpy scipy matplotlib mpmath
```

## Usage

```
python demo.py           # Run all demonstrations (text output)
python demo.py --plot    # Also generate matplotlib visualizations
```
