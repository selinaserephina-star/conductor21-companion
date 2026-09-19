# FINALE — four closing micro-probes on the day's caches

**2026-09-11 · Stenberg + Claude (Fable 5) · Merkabit_crystal, closes the
BK_OFFDIAG → CHEBRACE → GAPMAP → SIVAR arc before consolidation.**

## Scope lock (binding)

Exploration. NOT the wall (Op2/GRH(ψ)/R-209 — OPEN, claimed nowhere). Not RH/GRH.
Part D is a **consistency check in the deep-sweep tradition (β ≈ 0 framing:
prime-only corroboration, never proof, never wall progress)**. Nulls at equal
prominence. Framework credits: Rubinstein–Sarnak/Ng (A), Montgomery–Soundararajan
(B, D) — cited, not claimed.

## Parts

- **A — races at 10⁹** (CHEBRACE named successor): same five-class races on the
  full 10⁹ cache; one more octave; measured log-avg means vs μ = (−21,−1,+1,0,0);
  reconstruction correlation with the (unchanged) family zero tables.
- **B — Gaussianity of increments** (SIVAR extension): skewness + excess kurtosis
  of D(u) = Δψ_χ/√x across δ ∈ bulk regime (≤0.05), windows [10⁷,10⁸], [10⁸,10⁹],
  block-bootstrap errors. Prediction: Gaussian (CLT / MS-type). Any coherent
  non-Gaussianity = 3rd/4th-order zero correlations beyond GUE → bug-hunt first.
- **C — multiplicative pairs** (GAPMAP named successor): (p, 2p+1) Sophie-Germain
  and (p, 2p−1): ⟨a_p·a_q⟩ vs independence + joint 5×5 class χ²; HL-type count
  control 2C₂∫dx/(ln x · ln 2x). Prediction: blind (no congruence channel).
- **D — spectral-purity scaling**: Var(δ=0.0134, 0.03) across half-decade windows
  10⁶→10⁹ (exact) + [10⁹,10¹¹] (stored signal): measured/GUE ratio vs x. Under
  zeros-on-line: flat at 1. Anomalous growth ∝ x^(2β−1) would flag an off-line
  contribution — none expected; consistency-grade only.

## Inputs

10⁹ class caches (BK_OFFDIAG + GAPMAP), family zero tables (CHEBRACE `data/`),
`signal_100000000000.npy` (read-only). All reused, nothing recomputed.

## Files

`finale.py` → `FINDINGS_finale.md`, `results/finale_results.json`, `results/*.png`.
