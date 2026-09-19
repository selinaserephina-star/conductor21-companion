# DAY SUMMARY — 2026-09-11: six runs, one arc

**Stenberg + Claude (Fable 5) · Merkabit_crystal · all findings notes, scripts,
and result JSONs under `RUNS/`. Scope: exploration; NOT the wall; not RH/GRH;
conditional statements labelled; nulls at equal prominence.**

Object: L(s,χ₈), degree 8, Artin conductor 21¹⁰ = 3¹⁰7¹⁰, Γ_ℂ(s)⁴, ε = +1,
Trinks field x⁷−7x+3, Gal = PSL(2,7). Inputs: 153 certified zeros (γ ≤ 27.911);
family low-zero tables χ₃-pair (15), χ₆ (12), χ₇ (10); exact Frobenius classes
for all p ≤ 10⁹ and exact a_p to 10¹⁰ (455,052,511 primes; Chebotarev at 5
digits: .00595/.37500/.33331/.28574).

## 1 · BK_OFFDIAG — twin-blindness + the explicit-formula identity

- **Q1 (identity/calibration):** windowed resonance vs explicit-formula prime
  side: r = +1.000 (χ₈ and ζ control both); rms 0.904 = 0.904. Log-7 silence
  exact (+0.438 meas vs +0.436 tail-only pred; ζ sings −54.0). Class
  fingerprints (2A: tr=0, tr(Frob²)=8; 4A silent at both) as the local factors
  force. **Identity, not evidence — but an end-to-end coherence certificate**
  (zeros ⊗ conductor data ⊗ local factors).
- **Q2:** R₂ of the unfolded 153 fully GUE (χ² 23.2/24 bins); BK off-diagonal
  not detectable at N=153 (Δχ² = 2.0/1 dof) — bound only. (Later closed from
  the prime side by SIVAR.)
- **Q3 (headline):** over 440,310 twin pairs ≤ 10⁸: ⟨a_p·a_{p+2}⟩ = −0.00196 ±
  0.00148 (z = −1.32); joint 5×5 class χ² = 11.1/16 (p ≈ 0.80). Positive
  control: π₂(10⁸)/Hardy–Littlewood = 0.9996. **TWIN-BLIND.**

## 2 · CHEBRACE — races predicted from our own zero tables

- Bias means exact from character theory: μ_C = 1 − r₂(C), r₂ = 1+χ₆+χ₇+χ₈ →
  **(1A: −21, 2A: −1, 4A: +1, 3A: 0, 7: 0)**. Measured to 10⁸ (and 10⁹ in
  FINALE): consistent throughout; the 1A deficit is −71 primes at 10⁸ =
  −21·√x/ln x/168 on the nose.
- Zero-table reconstruction tracks the measured races at r = 0.80–0.96 (0.81–0.97
  at 10⁹) — **first dynamical validation of the χ₃-pair/χ₆/χ₇ tables**, and a
  measurement of Rubinstein–Sarnak low-zero dominance (10–15 zeros suffice).
- First bias densities for a PSL(2,7) field (conditional GRH+LI, labelled):
  P(1A trails each class) ≈ 0.71–0.73; δ(4A>3A) = 0.56; 3A vs 7 even.
  **Structural note: bias mean ∝ squares count but spread ∝ column norm
  (Σ deg χ² = 167) ⇒ large mean, weak density** — contrast classical mod-4
  (δ ≈ 0.996).

## 3 · GAPMAP — uniform gap-blindness, g ≤ 60, to 10⁹

- 30 even gaps × 3 windows (10⁸; 10⁹; disjoint decade), 3.0–6.9M pairs/gap at
  10⁹: **no gap over the Bonferroni threshold anywhere**; max |corr| = 0.00101;
  z-ensembles unit noise (std 1.13/0.95/0.96); Σχ² joint-class 475/480,
  480/480, 462/480.
- Same pairs reproduce the HL singular-series sawtooth ∏(p−1)/(p−2)
  (1.00→2.67): mean ratio **0.9989**, worst dev 0.0015. Blindness measured
  against a loudly structured background.
- **Cousin flag retired:** +2.76σ (g=4, 10⁸) → **−2.26σ (sign flip)** on the
  pre-registered disjoint decade. Textbook fluctuation death.

## 4 · SIVAR — the bulk form factor from primes

- Regime insight: 2πρ̄(γ) = ln Q − 8ln 2π + 8ln γ ≈ 15.7 + 8ln γ ⇒ α =
  ln x/(2πρ̄) < 1 for all γ ≳ 2 at x = 10⁸ — **the conductor is a form-factor
  microscope** (ζ at same (x,δ): α ≈ 3–4, diagonal only).
- **Montgomery ramp F(α)=α confirmed at ~1%**: bulk χ² vs GUE 5.5/20 and
  11.7/20 (windows 10⁷⁻⁸, 10⁸⁻⁹); vs uncorrelated zeros χ² ≈ 5×10⁵; W2
  measured/GUE = 1.009 at δ ≤ 0.01; three windows (to 10¹¹) ride three separate
  GUE curves (the ln x-dependence itself confirmed). Coverage α ≈ 0.22–0.74,
  heights to γ ~ 10⁴ (≈ 400× past the certified spectrum). **First two-level
  bulk statistics of this spectrum.**
- Deep-regime identity derived + measured: for α<1, F=α cancels ρ̄ ⇒
  **Var_GUE = δ·ln x**; measured bench/GUE = 1.001.
- ζ control through the identical pipeline: 0.996/0.994.
- Gate A: increments from your 2026-07 numba deep-sweep signal vs the flint
  cache: corr = 1.000000 — the two prime engines cross-certify.
- Method boundary recorded: δ ≥ 0.15 is realization-dominated (few lowest
  zeros; <½ oscillation per window) — excluded from fits by design.

## 5 · FINALE — four closers

- **Races at 10⁹ confirmed** (3A at −0.00; recon corr → 0.97).
- **Gaussianity:** naive 3.72σ skewness = exact shot noise (ΣΛ_χ³/x^{3/2},
  ⟨a³⟩ = 3); after predicting it: max|z| 0.85 (skew), 2.46 (kurt), 48 tests —
  3rd/4th-order zero correlations bounded at the Gaussian level.
- **Affine channel blind:** Sophie-Germain (p,2p+1): z = −0.28 over 1,775,672
  pairs; (p,2p−1): z = −0.46; SG count / HL-type prediction = **0.9999**.
- **Purity scaling:** Var/Var_GUE at fixed δ across x.

## 6 · ext_1e10 — the tail retest (exact engine to 10¹⁰)

- 404,204,977 primes in (10⁹,10¹⁰] in 46.5 min (overflow-safe numba kernel;
  note for reuse: NUMBA_THREADING_LAYER=workqueue needed on this machine).
- New exact windows: Var/GUE = 0.98/0.85/1.05/0.90 (±0.11–0.15).
- **Exact-only slope over 10⁶→10¹⁰: +0.022 ± 0.015 and +0.005 ± 0.018 — flat.**
  The earlier +2.3σ drift was the coarse legacy-signal corner; retired.
  Consistency-grade only (your deep-sweep β≈0 framing); the 10¹⁰–10¹¹ corner
  remains signal-only.

## The arc in one table

| level | run | verdict |
|---|---|---|
| one-point, zeros↔primes | BK_OFFDIAG Q1 | identity, exact (r = 1.000) |
| prime pairs, additive | BK_OFFDIAG Q3 + GAPMAP | uniformly blind (g ≤ 60, 10⁹) |
| prime pairs, affine | FINALE C | blind (SG control 0.9999) |
| one-point races | CHEBRACE (+FINALE A) | biases exact; zero tables validated |
| low zeros (RUN 72–74, thermo lane, cited) | | GUE-consistent; SO(even) |
| two-point bulk | SIVAR (+ext_1e10) | **GUE ramp @1% to γ~10⁴; flat to 10¹⁰** |
| 3rd/4th order | FINALE B | Gaussian after exact shot noise |

Flags raised during the day: 3. Flags surviving: 0 (sign-flip retest; exact
shot-noise prediction; exact-engine re-measurement).

— 2026-09-11. Not the wall; not RH; the spectrum simply behaves.
