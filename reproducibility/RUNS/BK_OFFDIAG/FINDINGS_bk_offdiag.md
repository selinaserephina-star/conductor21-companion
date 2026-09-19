# FINDINGS — BK_OFFDIAG: the genome is twin-blind; the off-diagonal machinery verifies as an identity

**2026-09-11 · Stenberg + Claude (Fable 5) · run `bk_explore.py`, log `run_bk_explore.log`,
results `results/bk_offdiag_results.json`. Inputs: 153 certified χ₈ zeros (parent-lane
master), 599 ζ ordinates (control), exact Frobenius classes of x⁷−7x+3 for all p ≤ 10⁸
(flint, 22 µs/prime).**

## §0. Scope lock (binding)

Exploration. **NOT the wall** (Op2 / GRH(ψ) / R-209 — OPEN, claimed nowhere). **Not
RH/GRH. No twin-prime-conjecture claim** — motivating question answered structurally,
below, and the answer is a *decoupling*, not progress on TPC. Computed/heuristic
throughout; refutations and nulls at equal prominence.

---

## §1. Q1 — the resonance instrument VERIFIES AS AN IDENTITY (calibration, not discovery)

Windowed resonance F(x) = Σⱼ w(γⱼ)cos(xγⱼ) minus the smooth main term, vs the
explicit-formula prime side −(1/2π)Σₙ Λ_χ(n)/√n·[Ŵ(x−ln n)+Ŵ(x+ln n)] with
Λ_χ(pᵏ) = Re(Σαᵢᵏ)ln p from the exact class local factors:

> **χ₈: r = +1.000** (resid rms 0.904 = pred rms 0.904; curves pointwise
> indistinguishable, `results/q1_resonance.png`). **ζ control: r = +1.000.**

**Honest framing — this is an identity check, not a statistical discovery.** The
explicit formula is exact; what the r = 1.000 buys is an **end-to-end coherence
certificate**: the certified zero list ⊗ the paper's analytic data (Q=21¹⁰, Γ_ℂ⁴)
⊗ the class local-factor tables are mutually consistent through a completely
independent pipeline (extends RUN 73's Gate A from the counting function to the
full arithmetic fluctuation). The June 21-zero probe
(`new_direction_exploration/zeros_to_primes_resonance.py`) is superseded: its
qualitative troughs are now amplitude-exact.

Fingerprints confirmed inside the identity:
- **log 7 silence** (ramified, local factor 1): χ₈ measured +0.438 vs predicted
  +0.436 (pure neighbour-kernel tail); ζ at log 7: −53.95 measured vs −51.08
  predicted (sings). The genome's spectrum audibly *omits* the ramified prime.
- p = 3 sings (local factor 1−X, Λ = ln 3 at all powers).
- Class tr(Frob²) checks: 2A silent at log p, loud at log p² (Σα² = 8); 4A silent
  at both (Σα² = 0); 3A gives −1 at both. All as the LOC8 tables force.

**Grade: VERIFIED (identity / instrument calibration).**

## §2. Q2 — pair correlation: GUE alone suffices; BK off-diagonal NOT detectable at N=153

Exact θ/π unfold (RUN 73's formula; span 151.44, mean spacing 0.9963). R₂ against
a same-estimator GUE surrogate band (301×301 GUE, central 153 levels, 300 draws;
uniform-surrogate edge correction M=2000):

> χ²(GUE alone) = **23.2 over 24 bins** — fully GUE-consistent (agrees with RUN 73/74).
> Adding the BK arithmetic shape Σₙ (Λ_χ(n)²/n)cos(d·ln n) with one fitted scale:
> Δχ² = 2.0 for 1 dof — **not significant**. Fitted scale −0.0045.

**Grade: MEASURED — NULL.** The BK off-diagonal is *bounded, not seen*, at this zero
count; at N=153 the two-level statistic cannot separate BK's conjectural content
from finite-N noise. The named (bounded, wall-untouched) successor is the certified
engine's production extension of the spectrum (SCALING_MEMO spec) — more zeros is
the only lever here.

## §3. Q3 — THE RESULT: the genome is TWIN-BLIND at 10⁸ (and the pipeline sees Hardy–Littlewood perfectly)

Exact Frobenius classes for all 5,761,455 primes ≤ 10⁸ (Chebotarev fractions match
to 4 digits: 1A .00594 / 2A .12497 / 3A .33305 / 4A .25008 / 7 .28596). Over prime
pairs, correlation of the genome trace a_p (null = independence; PSL(2,7) simple ⇒
no abelian quotient ⇒ no congruence channel ⇒ predicted 0):

| pair set | N pairs | ⟨a_p·a_p′⟩ | z vs null | joint 5×5 χ²/dof 16 |
|---|---|---|---|---|
| **twin (p, p+2)** | 440,310 | **−0.00196 ± 0.00148** | **−1.32** | **11.1** (p≈0.80) |
| cousin (p, p+4) | 440,256 | +0.00425 ± 0.00154 | +2.76 ⚠ | 15.0 (p≈0.52) |
| sexy (p, p+6) | 879,906 | +0.00129 ± 0.00111 | +1.16 | 24.7 (p≈0.08) |
| consecutive p | 5.76M | +0.00043 ± 0.00042 | +1.02 | — |

**Positive control (the knife that makes the null sharp):** the same prime set
gives π₂(10⁸) = 440,312 vs the Hardy–Littlewood prediction 440,473 — **ratio
0.9996**. The pipeline resolves the classical twin singular series to 4 digits
*while seeing no genome-class coupling whatsoever in the same pairs*.

⚠ Cousin flag, recorded honestly: +2.76σ in one of ~8 tests (look-elsewhere ⇒
global p ≈ 0.05), with an unremarkable joint-class table (15.0/16). Treated as a
fluctuation; cheap named follow-up = rerun `gap_corr(4)` at 10⁹ before it is ever
mentioned again. It does not survive as a claim.

> **Verdict: TWIN-BLIND.** The genome's field-analogue Hardy–Littlewood constant is
> the **trivial (independence) one**: C^(χ₈)_twin = −0.002 ± 0.0015 ≈ 0, exactly as
> the simplicity of PSL(2,7) forces. The classical twin constant C₂ = 0.660… lives
> on **congruence (abelian) data** — the singular series is a product over
> residue-class obstructions — and the genome's channel is **maximally non-abelian**,
> so it has nothing to couple to. A new member of the blindness family (marking
> identity: the cell reads φ only; twinning is likewise invisible to it).

**Grade: MEASURED — twin-blind at 10⁸; positive control PASSED.**

## §4. What this says about the motivating question (Selina, 2026-09-11: "can our work touch the twin prime conjecture?")

The three probes close the question from three sides, all pointing the same way:

1. **One-level (Q1):** the zero↔prime bridge on our object is an *identity* —
   nothing conjectural to harvest there.
2. **Two-level (Q2):** the BK heuristic — the only known conceptual bridge from
   zero statistics to prime-pair constants — is invisible at 153 zeros.
3. **Prime-pair side (Q3):** measured directly, the genome's arithmetic is
   **decoupled from twinning** (as its group theory predicts), while the classical
   twin structure (HL, abelian) is alive and exactly measurable in the same data.

**Honest bottom line: the conductor-21 program is structurally orthogonal to the
twin prime conjecture.** Its one genuine TPC-adjacent product is negative-space
knowledge: twin correlations are an abelian/congruence phenomenon, and a maximally
non-abelian channel measurably carries none of it. No route from this lane toward
TPC is claimed or implied.

## §5. Files

- `README_BRIEF.md` (brief + scope lock), `bk_explore.py`, `run_bk_explore.log`
- `results/bk_offdiag_results.json` (all numbers), `results/q1_curves.npz`
- `results/q1_resonance.png`, `results/q1_scatter.png`, `results/q2_r2.png`
- `results/classes_100000000.npz` (cached exact classes ≤ 10⁸ — reusable)
- Inputs under `data/` (read-only copies; parent lane untouched)

— Stenberg + Claude (Fable 5), 2026-09-11. Exploration; no registry row (Merkabit
is self-contained); parent lane sealed and untouched.
