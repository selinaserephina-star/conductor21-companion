# Prime and zero statistics of a degree-8 Artin L-function of conductor 21¹⁰

**An empirical companion — Selina Stenberg and Ilya Balashov** · v1.3 (2026-09-19) · code MIT · paper/data CC BY 4.0

Companion paper to *The Conductor-21 L-Function: Certified Zeros, Exact Laws, and a
Density–Prime Operator* ([concept DOI 10.5281/zenodo.21982347](https://doi.org/10.5281/zenodo.21982347);
data [10.5281/zenodo.21958629](https://doi.org/10.5281/zenodo.21958629)).
Zenodo DOI for this paper: *to be added at deposit*.

## What this is

An empirical study of the degree-8 Artin L-function L(s, χ₈) of the Galois closure of
the Trinks field x⁷−7x+3 (Galois group PSL(2,7), conductor 21¹⁰), from 153 certified
zeros and exact Frobenius data for every prime to 10¹⁰:

- **The explicit formula as an identity** — zeros ⇄ primes at correlation 1.000,
  including the forced silence of the ramified prime 7.
- **The blindness suite** — Frobenius carries no pair information (twins, all gaps
  g ≤ 60, affine pairs 2p±1) at the 10⁻³ level, while the same data reproduce the
  Hardy–Littlewood singular series to 0.1%.
- **Chebotarev races** — biases μ_C = 1 − r₂(C) (−21 for the split class) and race
  shapes predicted by our own zero tables at r = 0.81–0.97.
- **Low-zero statistics** — GUE-consistent, symmetry type SO(even).
- **The bulk form factor from primes** — the conductor puts every height γ ≳ 2 in the
  correlation-sensitive regime at x = 10⁸; the Montgomery ramp F(α) = α holds at ~1%
  to heights γ ~ 10⁴; uncorrelated zeros excluded at χ² ≈ 5×10⁵.

All results are consistency measurements against classical predictions (Montgomery;
Goldston–Montgomery; Hardy–Littlewood; Rubinstein–Sarnak; Katz–Sarnak). No novelty is
claimed for the frameworks, and the Riemann Hypothesis for L(s, χ₈) — the open problem
of the main paper — is untouched.

## Repository map

```
paper/            companion.pdf (v1.3), LaTeX source, figures
reproducibility/  per-run folders (BK_OFFDIAG, CHEBRACE, GAPMAP, SIVAR, FINALE
                  + ext_1e10): briefs, findings notes, scripts, logs, result
                  JSONs, plots; data_seed/ (Frobenius classes to 1e8 + the four
                  zero tables); DAY_SUMMARY.md; SHA256SUMS.txt
```

Every number in the paper is produced by a named script on named data. The larger
caches regenerate from the included scripts: Frobenius classes to 10⁹ (~17 min,
`gapmap.py`), exact a_p on (10⁹, 10¹⁰] (~47 min, `sweep_1e10.py`).

## Reproducing

Python 3 with `numpy`, `python-flint`, `mpmath`, `matplotlib`; the (10⁹, 10¹⁰] sweep
also uses `numba` (on some Windows setups the parallel kernel needs
`NUMBA_THREADING_LAYER=workqueue`). Each `reproducibility/RUNS/*/` script runs
standalone against `data_seed/` and its own cached inputs.

## Verification note

Independent review of this work verified the theoretical inputs on which the
statistics rest — the bias-mean formula and its square-root counts, the
Frobenius–Schur indicators, and the framework identities, re-derived from first
principles — and did not re-run the numerical computations; the scripts and
data are released precisely so that any reader may perform that verification.

## License and citation

Dual-licensed: **scripts and code under MIT**; **paper, figures, and data tables
under CC BY 4.0** (see `LICENSE` and `LICENSE-CC-BY-4.0`). Cite via
`CITATION.cff`, or the Zenodo DOI once deposited.
Prepared with AI assistance (Claude); see the paper's Appendix A for the
model-attribution record.
