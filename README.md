# Riemann Hypothesis: Computational Methods

Implementations of the most significant computational methods and modern
advances related to the Riemann Hypothesis.

## Modules

| Module | Description |
|--------|-------------|
| `riemann_siegel.py` | Riemann-Siegel formula for Z(t) on the critical line |
| `zero_finding.py` | Gram points, sign-change detection, Turing method verification |
| `zero_counting.py` | Hardy's Z-function, Backlund's N(T) zero-counting formula |
| `pair_correlation.py` | Montgomery's pair correlation conjecture & GUE comparison |
| `newman_constant.py` | De Bruijn-Newman constant Λ (Rodgers-Tao, 2018) |
| `moments.py` | Keating-Snaith moment conjectures from random matrix theory |
| `demo.py` | Demonstrations and visualizations |

## Key References

- Riemann, B. (1859). "Über die Anzahl der Primzahlen unter einer gegebenen Größe"
- Siegel, C.L. (1932). "Über Riemanns Nachlaß zur analytischen Zahlentheorie"
- Montgomery, H.L. (1973). "The pair correlation of zeros of the zeta function"
- Odlyzko, A.M. (1987). "On the distribution of spacings between zeros of the zeta function"
- Keating, J.P. & Snaith, N.C. (2000). "Random matrix theory and ζ(1/2+it)"
- Rodgers, B. & Tao, T. (2020). "The de Bruijn-Newman constant is non-negative"
- Platt, D.J. (2017). "Isolating some non-trivial zeros of zeta"

## Requirements

```
pip install numpy scipy matplotlib mpmath
```
