# FINDINGS — GAPMAP: independence holds uniformly; the genome is blind to ALL gap structure ≤ 60 at 10⁹; the cousin flag is retired by sign-flip

**2026-09-11 · Stenberg + Claude (Fable 5) · `gapmap.py` (84 min: 45.1M new exact
Frobenius classes in (10⁸,10⁹], checkpointed, cached), log `run_gapmap.log`,
numbers `results/gapmap_results.json`. Successor of BK_OFFDIAG (executes its menu
item 3 + its named cousin recheck).**

## §0. Scope lock (binding)

Exploration. NOT the wall (Op2/GRH(ψ)/R-209 — OPEN, claimed nowhere). Not RH/GRH,
no TPC claim. Measurements are exact counting; the null is the independence
heuristic (structurally forced: PSL(2,7) simple ⇒ no congruence channel).
Refutations/nulls at equal prominence — this note is two nulls and one retirement,
and all three are the result.

## §1. Q1 — the map: blindness is UNIFORM across all even gaps g ≤ 60

⟨a_p·a_{p+g}⟩ vs the independence null, 30 gaps × three windows (10⁸ full; 10⁹
full; disjoint decade (10⁸,10⁹]), 3.0M–6.9M pairs per gap at 10⁹:

| window | max\|z\| (gap) | \|z\|>2 | over Bonferroni (2.93) | z-ensemble std | Σχ² joint-class / dof |
|---|---|---|---|---|---|
| 10⁸ | 2.76 (g=4) | 2/30 | **0** | 1.13 | 475 / 480 |
| 10⁹ | 2.17 (g=6) | 2/30 | **0** | 0.95 | 480 / 480 |
| disjoint | 2.26 (g=4, **−**) | 1/30 | **0** | 0.96 | 462 / 480 |

Not one gap crosses the 30-test threshold in any window; the z-ensembles are
unit-variance noise (QQ plot `results/gapmap_qq.png`); the class-level χ² totals
sit dead on their dof. Quantitative blindness bound at 10⁹: **max |corr| over all
30 gaps = 0.00101**, per-gap se 0.00033–0.00054 — every gap consistent with zero
at ~2σ·0.001 precision.

> **Verdict: independence holds uniformly to ~10⁻³ per gap — the genome's
> Frobenius carries no gap information whatsoever, for any even g ≤ 60, to 10⁹.**
> First measurement of its kind for a PSL(2,7) extension, to our knowledge
> (computed; negative-space; no novelty asserted beyond the measurement).

**Grade: MEASURED — UNIFORM NULL (the result).**

## §2. Q2 — the contrast control: Hardy–Littlewood's singular series confirmed across the whole sweep

Same pairs, classical count: N_g / (2C₂·Li₂) vs ∏_{odd p|g}(p−1)/(p−2), which
swings from 1.00 (g=2,4,8,16,…) to 2.67 (g=30, 60):

> mean ratio measured/HL = **0.9989** at 10⁹; worst deviation across all 30 gaps
> **0.0015**. The full sawtooth (2, ×2 at 3|g, ×2.4 at 21|g, ×2.67 at 105-free
> 30|g…) is reproduced point by point (`results/gapmap.png`, lower panel).

The sub-0.2% uniform offset is the expected truncation of the HL asymptotic (plus
our p>7 exclusion) — not a discrepancy. The knife this provides: **the same pair
data resolves the classical abelian gap structure to 3 digits while the
non-abelian channel shows nothing at 10⁻³** — blindness measured against a loudly
structured background, not against silence.

**Grade: MEASURED — POSITIVE CONTROL PASSED (sweep-wide).**

## §3. Q3 — the cousin flag is RETIRED

Pre-registered test (BK_OFFDIAG findings, fixed before this data existed): the
g=4 flag (+2.76σ at 10⁸) survives only if the disjoint decade shows z ≥ +2.8,
same sign.

> Disjoint decade (10⁸,10⁹], N = 2,984,422 pairs: **z = −2.26. Sign flipped.**
> **DOES NOT SURVIVE — fluctuation, retired** (as the look-elsewhere arithmetic
> predicted). The full-range 10⁹ value: z(g=4) = −1.10.

Textbook behaviour of a 1-in-30 fluctuation under independent retest. The
BK_OFFDIAG caveat is now closed; g=4 must never be quoted as a signal.

**Grade: RETIRED (pre-registered, independent data).**

## §4. Assets and successors

- **Reusable asset:** exact Frobenius classes for ALL 50,847,534 primes ≤ 10⁹
  (`results/classes_1e8_1e9.npz` + the BK_OFFDIAG 10⁸ cache). Any future
  arithmetic statistic on this field to 10⁹ is now seconds, not hours.
- Chebotarev at 10⁹: 1A .00595 / 2A .12500 / 3A .33335 / 4A .25006 / 7 .28564 —
  5-digit agreement.
- Bounded successors (none urgent): odd-prime-power gaps or (p, 2p±1)
  Sophie-Germain-type pairs on the same cache; the CHEBRACE 10⁹ race extension
  now costs seconds with this cache.

## §5. Files

`README_BRIEF.md`, `gapmap.py`, `run_gapmap.log`, `results/gapmap_results.json`,
`results/gapmap.png`, `results/gapmap_qq.png`, `results/classes_1e8_1e9.npz`.
Parent lane untouched; Merkabit self-contained, no registry row.

— Stenberg + Claude (Fable 5), 2026-09-11. Two nulls and a retirement, all three
load-bearing: the blindness family gains its sharpest member.
