# FINDINGS — SIVAR: the genome's bulk form factor is the GUE ramp, measured from primes at the 1% level

**2026-09-11 · Stenberg + Claude (Fable 5) · `sivar.py` (46 s against the 10⁹ class
caches), log `run_sivar.log`, numbers `results/sivar_results.json`. Successor of
GAPMAP; executes menu item 2 (Goldston–Montgomery / Montgomery–Soundararajan
direction — framework cited, not claimed). Closes what BK_OFFDIAG Q2 could not
reach with 153 zeros.**

## §0. Scope lock (binding)

Exploration. NOT the wall (Op2/GRH(ψ)/R-209 — OPEN, claimed nowhere). Not RH/GRH:
the prediction formulas place zeros on the ½-line and use the stationary cluster
expansion (GUE Y₂ vs Poisson) — heuristic grade, stated; the measurement is exact
prime counting. "First bulk two-level statistics of the genome spectrum" is a
measurement claim only; GUE statistics are the *expected* behaviour
(consistency, not discovery — RUN 73 §0 wording). Refutations/nulls at equal
prominence.

## §1. The design (why this object can do what ζ cannot)

Var over x of [ψ_χ(xe^δ)−ψ_χ(x)]/√x equals, under the cluster expansion,
∫ρ̄(γ)·2|φ_γ(δ)|²·F(α)dγ with α = ln x/(2πρ̄(γ)) and F the Montgomery form
factor: **F = min(α,1) for GUE, F ≡ 1 for uncorrelated zeros.** For the genome,
2πρ̄ = ln(21¹⁰) − 8ln 2π + 8ln γ — the conductor makes the spectrum so dense that
**α < 1 for every γ ≳ 2 already at x = 10⁸**: the whole variance is
correlation-sensitive and GUE vs Poisson differ by ×2–5. For ζ at identical
(x, δ), α ≈ 3–4 ⇒ F = 1 ⇒ zero discrimination — the built-in control.
Derived en route (§4): in the deep-conductor regime F=α makes ρ̄ cancel exactly,
Var_GUE → δ·ln x = the independent-prime-signs benchmark.

## §2. Gates

- **Gate A — two independent prime engines agree:** increments from the 2026-07
  numba deep-sweep signal vs this month's flint class cache, overlapping window,
  δ=0.045: **corr = 1.000000, rms ratio = 1.0000.** (First cross-validation of
  the two engines; retro-certifies the deep-sweep signal data.)
- **Gate B — ζ control:** measured/predicted = **0.996 (W1), 0.994 (W2)**,
  χ² 10 and 8 on 26 points each. The pipeline and its normalization are
  calibrated at the half-percent level — including the α>1 physics (the
  prime-correlation suppression below the naive benchmark that GM predicts
  for ζ).

## §3. THE RESULT — the Montgomery ramp F(α)=α, measured in the genome's bulk

26 interval lengths δ ∈ [10⁻⁴, 0.3] × windows W1 = [10⁷,10⁸], W2 = [10⁸,10⁹]
(exact primes), W3 = [10⁹,10¹¹] (stored deep-sweep signal, coarse):

| window | regime | χ² vs GUE | χ² vs Poisson | mean R = Var/Var_Poisson |
|---|---|---|---|---|
| W1 | bulk δ ≤ 0.05 (n=20) | **5.5** | 585,788 | measured 0.327 · GUE 0.316 |
| W2 | bulk δ ≤ 0.05 (n=20) | **11.7** | 470,710 | measured/GUE = **1.009** at δ≤0.01 |
| W3 | all 8 δ | rides the curve within errors | — | e.g. 0.457±0.027 vs GUE 0.428 |

- The measured ratio R(δ) tracks **three separate GUE ramp curves** — the
  window-dependence (α ∝ ln x pushes W1 < W2 < W3) is itself a confirmed
  prediction (`results/sivar_ratio.png`).
- **Uncorrelated zeros (Poisson) are excluded at χ² ≈ 5×10⁵**; the variance is
  suppressed to ~20–50% of the uncorrelated value exactly as the ramp dictates.
- Effective coverage: **α ≈ 0.22–0.74, zero heights γ up to ~10⁴** — four
  hundred times beyond the certified spectrum (γ ≤ 27.9). This is the first
  two-level measurement of the genome's bulk, and it is **GUE to ~1%** where the
  data is strongest.
- Large-δ caveat (recorded, expected): at δ ≥ 0.15 the increments are dominated
  by the few lowest actual zeros — one realization, under half an oscillation
  per window — so ensemble predictions and block-bootstrap errors both degrade
  (W2's full-range χ² of 96 is entirely these points: 96 → 11.7 on δ ≤ 0.05).
  Not physics; noted as the boundary of the method.

**Grade: MEASURED — GUE CONFIRMED in the bulk at the 1% level (consistency with
the expected universality class; Poisson decisively excluded).**

## §4. The identity that fell out (worth keeping)

For α < 1, F = α makes the density cancel: **Var_GUE = δ·ln x exactly — the
independent-prime-signs benchmark.** Measured: bench/GUE = **1.001** at δ = 10⁻³.
Interpretation: in the conductor-dominated regime, GUE zero statistics ⟺ the
prime increments carrying no correlations at all — and the measured variance
lands on both. This ties SIVAR to GAPMAP: the integrated gap-correlation over
~xδ gaps (which GAPMAP's per-gap nulls could not constrain) is measured here to
be ≈ 0 at the percent level, consistent with the per-gap blindness *and* with
GUE zeros. The two runs are the same statement seen from the two ends of the
explicit formula.

## §5. Files

`README_BRIEF.md`, `sivar.py`, `run_sivar.log`, `results/sivar_results.json`,
`results/sivar_var.png`, `results/sivar_ratio.png`. Inputs: BK_OFFDIAG + GAPMAP
class caches (≤10⁹), `genome_deep_zeros_raw/signal_100000000000.npy` (read-only,
parent lane untouched).

— Stenberg + Claude (Fable 5), 2026-09-11. The day's arc closes: one-level =
identity (BK_OFFDIAG Q1), prime-pair channel = blind (GAPMAP), two-level bulk =
GUE ramp at 1% (SIVAR). Not the wall; not RH; the spectrum simply behaves.
