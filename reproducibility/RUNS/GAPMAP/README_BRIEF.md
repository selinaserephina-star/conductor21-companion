# GAPMAP — systematic gap-blindness map for the genome (g ≤ 60, to 10⁹)

**2026-09-11 · Stenberg + Claude (Fable 5) · Merkabit_crystal, successor of
BK_OFFDIAG (menu item 3 + the named cousin recheck). Reuses the 10⁸ Frobenius
cache; extends it to 10⁹.**

## Scope lock (binding)

Exploration. NOT the wall (Op2/GRH(ψ)/R-209 — OPEN, claimed nowhere). Not RH/GRH,
no TPC claim. The joint distribution of Frobenius classes at (p, p+g) for a
non-abelian extension is theorem-free territory; measurements here are exact
counting, predictions are the independence heuristic (forced structurally by
PSL(2,7) simplicity ⇒ no congruence channel). Any gap showing coherent structure
at both decades gets the full bug-hunt treatment BEFORE being called a signal
(house rule). Refutations/nulls at equal prominence.

## Questions

- **Q1 (the map):** for every even gap g = 2…60 at X = 10⁸ then 10⁹:
  ⟨a_p·a_{p+g}⟩ ± se and z vs the independence null; joint 5×5 class χ² (dof 16).
  Global verdicts via the z-ensemble (max|z|, Bonferroni at 30 tests, QQ vs
  normal) — never a single-gap headline.
- **Q2 (the contrast control):** measured pair counts N_g vs the classical
  Hardy–Littlewood singular series 2C₂·∏_{odd p|g}(p−1)/(p−2)·Li₂(X). The
  classical abelian gap structure varies by ~2.7× across the sweep; the genome
  curve is predicted flat at 0. One picture: rich classical structure, blind
  non-abelian channel.
- **Q3 (cousin verdict):** the BK_OFFDIAG +2.76σ flag at g=4 retested on the
  DISJOINT decade (10⁸, 10⁹] — independent data, pre-registered single test:
  survives only if the disjoint-decade z is again ≥ +2.8 with the same sign.

## Inputs / compute

- `../BK_OFFDIAG_2026-09-11/results/classes_100000000.npz` (exact classes ≤ 10⁸).
- Stage B: segmented sieve to 10⁹ (RAM-safe on the 16 GB machine) + flint
  Frobenius for the 45.1M new primes (~17 min, cached to
  `results/classes_1e8_1e9.npz`).

## Files

`gapmap.py` → `FINDINGS_gapmap.md`, `results/gapmap_results.json`, `results/*.png`.
