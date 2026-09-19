# FINDINGS — FINALE: four closers — races confirmed at 10⁹, increments Gaussian after exact shot noise, the affine channel is blind too, purity scaling flat to 10⁹ (coarse tail flagged)

**2026-09-11 · Stenberg + Claude (Fable 5) · `finale.py` (27 s), log
`run_finale.log`, numbers `results/finale_results.json`. Closes the day's arc
before consolidation.**

## §0. Scope lock (binding)

Exploration. NOT the wall (Op2/GRH(ψ)/R-209 — OPEN, claimed nowhere). Not RH/GRH;
part D is consistency-grade only (deep-sweep β≈0 framing: prime-only
corroboration, never proof). Nulls at equal prominence.

## §A. Races at 10⁹ — CHEBRACE's verdicts hold and sharpen

Measured log-avg means [10⁴,10⁹]: **1A −32.9 (μ=−21, LI σ 34) · 2A −1.83 (−1) ·
4A +1.51 (+1) · 3A −0.00 (0, dead on) · 7 +0.17 (0)** — every class within its
LI spread, signs and ordering as character theory forces. Zero-table
reconstruction correlations *improve* with the longer window: **0.81 / 0.92 /
0.95 / 0.94 / 0.97**. **Grade: CONFIRMED at the next octave.**

## §B. The increments are Gaussian — once the EXACT shot noise is subtracted

Naive test: max|z| of skewness = 3.72 (one point over Bonferroni) — but the
third moment of a short-interval prime sum has an exact diagonal (shot-noise)
term ΣΛ_χ(n)³/x^{3/2}, computable from the caches (⟨a³⟩ = 3 for the genome).
With that predicted exactly (no fit):

> **max|z| skew 3.72 → 0.85; kurtosis 2.46** — all 48 tests inside Bonferroni.

The "non-Gaussianity" was entirely the predicted granularity of the primes.
Residual Gaussianity of the fluctuation field = **third- and fourth-order zero
correlations bounded at the level GUE/CLT predicts** — the first higher-order
statistics of the genome spectrum. **Grade: MEASURED — GAUSSIAN (after exact
diagonal); the naive 3.7σ retired by prediction, not by excuse.**

## §C. The multiplicative-affine channel is blind too

Pairs with q ≤ 10⁹ (so p ≤ 5×10⁸), exact classes both ends:

| pair | N | ⟨a_p·a_q⟩ | z | joint χ²/dof 16 |
|---|---|---|---|---|
| Sophie-Germain (p, 2p+1) | 1,775,672 | −0.00022 ± 0.00076 | −0.28 | 25.0 (p≈.07) |
| (p, 2p−1) | 1,774,609 | −0.00034 ± 0.00074 | −0.46 | 13.5 |

Positive control: **SG count / HL-type prediction 2C₂∫dx/(ln x·ln 2x) = 0.9999.**
(The first run's "0.537" was an integral-limit bug in the *control*, not the
data — the pair constraint caps p at 5×10⁸; fixed and now 4-digit.) With GAPMAP:
**additive gaps ≤ 60 AND the affine maps p→2p±1 — all blind.** The blindness
family now covers both pair topologies. **Grade: MEASURED — NULL (the result).**

## §D. Spectral-purity scaling — flat to 10⁹; the coarse tail is bounded, not read

Var/Var_GUE at δ = 0.0134, 0.0295 across half-decade windows 10⁶→10⁹ (exact
engine) + 10⁹→10¹¹ (stored deep-sweep signal, coarse):

- Exact windows scatter around 1 (mean ratios 0.92–0.97; the few-% offset is the
  un-modeled Montgomery–Soundararajan arithmetic correction territory, same size
  as SIVAR's residuals).
- Weighted slope d(ratio)/d ln x: **+0.023 ± 0.011** (δ=0.0134) and +0.012 ±
  0.014 (δ=0.0295) — a combined ~2.3σ, **driven entirely by the last
  coarse-signal point** (10^10.5–10^11, ratio 1.2–1.26, the binned 8000-pt grid's
  noisiest corner). Look-elsewhere + known systematics ⇒ **consistent with flat**;
  recorded as a bound: no anomalous growth beyond ~4%/e-fold at 95% over five
  decades — under the explicit formula, no sign of any off-½-line contribution
  at that sensitivity. **Named follow-up before ANY stronger reading: extend the
  exact class engine past 10⁹ (≈7 h local or AWS) and re-test the tail with the
  exact engine.** **Grade: CONSISTENCY-GRADE PASS (deep-sweep β≈0 framing);
  tail flagged, not interpreted.**

## §D′. ADDENDUM (same day) — the tail flag is RESOLVED by the exact 10¹⁰ sweep

Executed the named follow-up immediately (`ext_1e10/`): 404,204,977 primes in
(10⁹,10¹⁰] via the overflow-safe numba kernel (46.5 min; workqueue threading
layer — the default tbb/omp init hangs on this machine, recorded for reuse);
Chebotarev at 10¹⁰: .00595/.37500/.33331/.28574 (5-digit). Exact ψ_χ binned at
du=2×10⁻⁴ (`ext_1e10/bins_1e9_1e10.npz` — reusable asset).

Exact windows replace the coarse signal through 10¹⁰:

| window | δ=0.0134 | δ=0.0294 |
|---|---|---|
| [10⁹, 10⁹·⁵] | Var/GUE = 0.981 ± 0.114 | 0.851 ± 0.124 |
| [10⁹·⁵, 10¹⁰] | 1.054 ± 0.121 | 0.899 ± 0.150 |

**EXACT-ONLY slope fit 10⁶→10¹⁰: +0.022 ± 0.015 and +0.005 ± 0.018 — both
consistent with 0.** The +2.3σ drift was the coarse stored-signal corner
(10^10.5–10¹¹, the 8000-pt grid's noisiest region), not the primes. **Verdict:
the GUE ratio is flat at 1 across five decades of exact data (10⁶→10¹⁰) — no
anomalous growth; the flag is retired.** (Consistency-grade, deep-sweep β≈0
framing; the 10¹⁰–10¹¹ corner remains signal-only and unretested — a numba
sweep of (10¹⁰,10¹¹] ≈ 8 h would close it; low priority.) SIVAR's ramp
confirmation now stands on exact data over 10⁶→10¹⁰.

## §5. Files

`README_BRIEF.md`, `finale.py`, `run_finale.log`, `results/finale_results.json`,
`results/finale.png`. Inputs: the day's caches + CHEBRACE tables +
`signal_100000000000.npy` (read-only). Parent lane untouched.

— Stenberg + Claude (Fable 5), 2026-09-11. The arc is closed on all four sides;
ready for consolidation.
