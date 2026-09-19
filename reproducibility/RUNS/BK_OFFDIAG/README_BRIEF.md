# BK_OFFDIAG — exploratory: the off-diagonal (Bogomolny–Keating) side of the genome zeros

**2026-09-11 · Stenberg + Claude (Fable 5) · Merkabit_crystal (user-directed home;
distinct thread from the PARKED conductor-21 scaffolding — this extends nothing there).**

## Scope lock (binding)

Exploration, instrument-class. **NOT** an attempt on the wall (Op2 deficiency /
GRH(ψ) / R-209 — OPEN, claimed nowhere). **NOT a twin-prime-conjecture claim of any
kind**: the pair-correlation ↔ prime-pair bridge (Montgomery; Bogomolny–Keating
1996; Conrey–Snaith) is a *heuristic*, and everything here is computed/heuristic,
no novelty asserted on the classical TPC. Not RH/GRH. Refutations recorded at
equal prominence (house rule).

**Motivating question (Selina, 2026-09-11):** can the genome program touch the twin
prime conjecture? Honest answer (recorded in session): not as a proof route — RH-type
statements don't imply TPC, and the sieve route needs the Dirichlet family, not one
L-function. The one genuine contact is the BK heuristic: zero *statistics* ↔ prime
*pair* arithmetic. This run probes exactly that contact, on our own object.

## Inputs

- `data/chi8_zeros_153.csv` — the 153 certified χ₈ zeros (parent lane master,
  `chi8_certified_zeros_MASTER_2026-07-14`, γ ≤ 27.911; read-only copy, lane untouched).
- `data/zeta_zeros_599.txt` — 599 classical ζ ordinates (control / instrument calibration).
- `data/gen_an_chi8_REFERENCE.py` — the class → local-factor tables (reference copy;
  the run imports its `frob_class` machinery for exact Frobenius classes of x⁷−7x+3).

Analytic data (paper, validated by RUN 73 gate A): degree 8, Q = 21¹⁰, Γ_ℂ⁴,
θ(T) = (T/2)ln Q − 4T ln 2π + 4 arg Γ(½+iT), N̄ = θ/π + 1 (slope 0.99895 on the 153 —
no fit). Symmetry type SO(even) (RUN 74). Bulk GUE-consistent, Var(s)=0.145 (RUN 73).
**None of that is redone here.**

## Predecessors (check-existing rule)

- `new_direction_exploration/zeros_to_primes_resonance.py` (June, 21 zeros):
  qualitative — troughs of Σcos(xγ) at x = log pᵏ, p=2 deepest. This run is the
  quantitative version at 153 zeros with exact per-class amplitude predictions.
- RUN 72/73/74 (thermo lane): unfold + one/two-level statistics — cited, not repeated.
- `genome_deep_zeros_raw/` (primes → zeros to 10¹¹): the dual direction, cited.

## Questions

- **Q1 (zeros → primes, amplitude-resolved).** Windowed resonance
  F(x) = Σⱼ w(γⱼ)cos(xγⱼ) minus its smooth-density main term. Prediction (explicit
  formula): dips at x = log pᵏ with amplitude ∝ Λ_χ(pᵏ)/pᵏ⸍² where
  Λ_χ(pᵏ) = (Σᵢ αᵢᵏ) ln p, αᵢ the exact local-factor roots of the class of Frob_p.
  Sharp discriminants: **silence at log 7** (ramified, local factor 1 — ζ sings there);
  p=3 sings (local factor 1−X); 2A-primes silent at log p but **loud at log p²**
  (Σα² = 8); 4A silent at log p AND log p² (Σα²=0), loud at log p⁴. Gate: fitted
  amplitudes vs predictions across prime-power targets + composite nulls;
  ζ control run through the identical pipeline calibrates the overall constant.
- **Q2 (BK off-diagonal in R₂).** Pair-difference residual (measured − GUE) of the
  unfolded 153 vs the arithmetic reconstruction Σₙ (Λ_χ(n)²/n)·cos(d·log n)-type
  term. Honest expectation: at N=153 this may only be bounded, not detected.
- **Q3 (twin-Frobenius correlation — the HL-analogue number).** Frobenius classes of
  all p ≤ 10⁸ (exact, flint). Over twin pairs (p, p+2): ⟨a_p·a_{p+2}⟩ and the joint
  5×5 class table vs the independence null. Structural prediction: PSL(2,7) is
  **simple** ⇒ no abelian quotient ⇒ no congruence coupling ⇒ the genome should be
  **twin-blind**: correlation 0 ± 1/√N_twin (a blindness result in the
  marking-identity family, if it holds). Controls: (p,p+4), (p,p+6), (p,next-prime).
  Either verdict is a result; a nonzero at >4σ would first be treated as a bug/bias.

## Files

`bk_explore.py` (the run) → `FINDINGS_bk_offdiag.md`, `results/*.json`, `results/*.png`.
