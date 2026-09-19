# CHEBRACE — the Chebotarev race for the Trinks field, predicted from our own zeros

**2026-09-11 · Stenberg + Claude (Fable 5) · Merkabit_crystal (successor thread to
BK_OFFDIAG; user-directed). Framework: Rubinstein–Sarnak 1994; non-abelian Chebotarev
version per Ng (2000). Our contribution is the *object*: first race analysis of a
PSL(2,7) septic (Trinks x⁷−7x+3) with the bias predicted from its own computed
zero tables — no priority claim on the framework.**

## Scope lock (binding)

Exploration. NOT the wall (Op2/GRH(ψ)/R-209 — OPEN, claimed nowhere). Not RH/GRH.
The *prediction* side assumes GRH-for-the-family + LI (standard R–S hypotheses —
stated, not asserted); the *measurement* side (races to 10⁸) is unconditional exact
counting. Computed/heuristic; refutations at equal prominence.

## Inputs

- Family zero tables (ours): χ₃-pair 15 zeros (γ≤5.31), χ₆ 12 (γ≤4.85), χ₇ 10
  (γ≤4.14) — `Zero_table_extensions_from_D_2026-07-11` — and χ₈ 153 certified
  (γ≤27.91, parent master). Read-only copies in `data/`.
- Exact Frobenius classes for all p ≤ 10⁸: cached
  `../BK_OFFDIAG_2026-09-11/results/classes_100000000.npz` (flint, verified
  Chebotarev to 4 digits).

## The mathematics being tested

Normalized race variables vs the unramified-π benchmark (so ζ's zeros cancel):
E′_C(x) = (ln x/√x)·[(|G|/|C|)·π_C(x) − π_unram(x)], observable classes
C ∈ {1A, 2A, 4A, 3A, 7=7A∪7B} of G = PSL(2,7), |G|=168.

- **Bias means (exact character theory, no central zeros — supported by the tables):**
  μ_C = 1 − r₂(C), r₂ = 1 + χ₆ + χ₇ + χ₈ (FS-positive characters)
  ⇒ **1A: −21 · 2A: −1 · 4A: +1 · 3A: 0 · 7: 0.**
  The squares bias: fully-split primes trail hugely; 4A (nothing squares into it —
  G has no order-8 elements) is the only positively-biased class.
- **Oscillation (the race's shape):** −Σ_χ w_χ(C)·Σ_{γ table} 2Re[x^{iγ}/(½+iγ)],
  w = class-averaged χ̄ (all real after merging 7A/7B, which the septic cycle type
  forces anyway). Reconstruction from OUR tables overlaid on the measured races —
  the first dynamical test of the χ₃-pair/χ₆/χ₇ tables.
- **LI lead-densities:** Monte Carlo over independent zero phases (shared per zero
  across classes ⇒ correct cross-class covariance), tail variance beyond each table
  estimated from end-of-table density (flat-density extrapolation, flagged).

## Questions

- **Q1:** do the measured races to 10⁸ carry the predicted means (−21,−1,+1,0,0)?
- **Q2:** does the zero-table reconstruction track the measured race shapes
  (per-class correlation over ln x)? [validates the family tables dynamically]
- **Q3:** LI lead-probabilities for the 10 pairwise races + P(E′_C<0) per class —
  the first computed bias densities for this field (conditional numbers, labeled).

## Files

`chebrace.py` → `FINDINGS_chebrace.md`, `results/*.json`, `results/*.png`.
