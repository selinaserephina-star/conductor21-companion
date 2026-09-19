# FINDINGS — CHEBRACE: the Trinks-field Chebotarev races run as our zeros predict; first bias densities for a PSL(2,7) field

**2026-09-11 · Stenberg + Claude (Fable 5) · `chebrace.py` (3.2 s against the cached
10⁸ class data), log `run_chebrace.log`, numbers `results/chebrace_results.json`.
Framework: Rubinstein–Sarnak; non-abelian treatment per Ng (2000) — no framework
priority claimed. Our contribution: the object (Trinks x⁷−7x+3, PSL(2,7)) and the
fact that the predictions come from our own zero tables.**

## §0. Scope lock (binding)

Exploration. NOT the wall (Op2/GRH(ψ)/R-209 — OPEN, claimed nowhere). Not RH/GRH.
Prediction side is **conditional** (GRH for the family + LI, the standard R–S
hypotheses); the measurement side is **unconditional exact counting to 10⁸**.
χ₈ zeros certified; χ₃-pair/χ₆/χ₇ zeros numerical-anchored (PARI) — graded so.

## §1. Q1 — the squares bias is REAL and carries the predicted structure

Exact character theory gives the bias means of E′_C = (ln x/√x)[(168/|C|)π_C − π_unram]:
**μ = 1 − r₂(C)**, r₂ = 1+χ₆+χ₇+χ₈ ⇒ **1A: −21, 2A: −1, 4A: +1, 3A: 0, 7: 0.**

Measured log-averages over [10⁴, 10⁸]:

| class | predicted μ | measured | LI 1σ spread |
|---|---|---|---|
| 1A (split) | **−21** | **−31.6** | ±34 |
| 2A | −1 | −1.17 | ±5.6 |
| 4A | +1 | −0.20 | ±4.5 |
| 3A | 0 | +0.26 | ±3.5 |
| 7 | 0 | +1.04 | ±3.6 |

2A and 3A land on their predictions to ~0.2; 1A is deep in the predicted trench
(−31.6, i.e. currently below its −21 mean — and the reconstruction §2 shows *why*:
the window sits in a low phase of the slow modes). One honest caveat: ln(10⁸/10⁴)
≈ 9.2 is barely **one period of the slowest oscillations** (2π/γ₁ ≈ 6–13), so the
log-window average is ~one draw, not an ensemble; the ±spreads reflect that.
The 1A deficit is even visible naked-eye in the raw Chebotarev fractions
(.00594 vs .005952 = −71 primes at 10⁸ ≈ −21·√x/ln x/168 on the nose).

**Grade: MEASURED — CONSISTENT.** The fully-split primes trail, hard, exactly as
the squares-count r₂(1A)=22 forces; 4A (the only class nothing squares into — G
has no order-8 elements) is the only positively-biased class.

## §2. Q2 — the zero tables reproduce the races' SHAPE (first dynamical validation of the family tables)

Overlaying −Σ_χ w_χ(C)·Σ_γ 2Re[x^{iγ}/(½+iγ)] + μ_C from OUR tables (χ₃-pair 15,
χ₆ 12, χ₇ 10, χ₈ 153 zeros) on the measured races:

> correlation(measured, reconstruction) over [10⁴,10⁸]:
> **1A +0.796 · 2A +0.923 · 4A +0.934 · 3A +0.944 · 7 +0.956**

Five orders of magnitude of race wiggles — 2A's +14 excursion at 10⁶, 7's surge at
10⁸, 3A's late dive — reproduced from a few dozen numerically-anchored zeros
(`results/cheb_races.png`). Honest decomposition of what this is: the explicit
formula (identity) **plus** two nontrivial facts it certifies — (i) the R–S
low-zero dominance is quantitatively true here (10–15 zeros per character suffice
for r ≈ 0.95), and (ii) the χ₃-pair/χ₆/χ₇ tables are **correct and complete at low
height** — a single missed or misplaced low zero would visibly break the tracking.
These three tables had never been dynamically tested before.

**Grade: VERIFIED (tables validated dynamically; low-zero dominance measured).**

## §3. Q3 — first computed bias densities for a PSL(2,7) field (conditional numbers)

LI Monte Carlo (200k draws, phases shared across classes for correct covariance;
tail variance beyond each table from end-density estimate — chi3pair/6/7 tails
add ~13–32% to their variance factors, flagged):

- **P(E′_1A < 0) = 0.727**; pairwise: 1A trails each other class with lead-prob
  **δ(1A > X) ≈ 0.27–0.29**; 4A leads 3A and 7 at **δ ≈ 0.56**; 3A vs 7 dead even
  (0.498). Full 10-pair table in the JSON.
- Measured lead fractions in the (single-sample) 10⁴–10⁸ window agree in sign and
  ordering everywhere (1A trailed 93–97% of the log-window).

**The qualitative finding worth keeping:** despite the enormous bias *mean* (−21!),
the 1A bias *density* is only ≈ 0.73 — far weaker than the classical mod-4 race
(δ ≈ 0.996). Reason, now measured: the identity class couples to every character
with weight deg χ, so its oscillation variance is ∝ Σ deg²·S_χ (the |G|−1 = 167
column norm), and σ ≈ 34 ≫ |μ| = 21. **In a strongly non-abelian race the
split-prime bias is large in mean but weak in probability — mean grows like the
squares count, spread grows like the full column norm.** (Conditional on GRH+LI;
computed, no novelty claim on the mechanism — but this number set for a PSL(2,7)
field is, to our knowledge, first.)

**Grade: COMPUTED (conditional, labeled).**

## §4. Caveats of record

1. Prediction side conditional (GRH family + LI); measurement unconditional.
2. 7A/7B merged (the septic's cycle type cannot split them) — all class-averaged
   character weights real, attribution ambiguity of the χ₃-pair zeros vanishes.
3. Central zeros assumed absent (m_χ(½)=0) — supported by all four tables
   (lowest γ = 0.48); a hidden central zero would shift the means.
4. Tail variance = flat-density extrapolation past table ends (underestimates the
   log growth slightly); χ₈ tail negligible (T=27.9).
5. Benchmark is π_unram (not Li), so ζ's zeros cancel by construction and only the
   family drives the races.

## §5. Successors (named, bounded, wall-untouched)

- Extend χ₆/χ₇/χ₃-pair tables past t≈5 (PARI, cheap) → reconstruction residual
  rms (currently 1.0–2.6, 1A 11.6) should drop like the missing-tail power.
- 10⁹ class sweep (~25 min) → one more octave of race, sharpens §1 means.
- The gap-blindness map (menu item 3) shares the same cache.

## §6. Files

`README_BRIEF.md`, `chebrace.py`, `run_chebrace.log`,
`results/chebrace_results.json`, `results/race_curves.npz`,
`results/cheb_races.png`, `results/cheb_bias_bars.png`; inputs `data/` (read-only
copies; parent lane untouched). Class cache reused from `../BK_OFFDIAG_2026-09-11`.

— Stenberg + Claude (Fable 5), 2026-09-11. Exploration; Merkabit is self-contained,
no registry row; parent lane sealed and untouched.
