# SIVAR — short-interval variance: the genome's bulk form factor measured from primes

**2026-09-11 · Stenberg + Claude (Fable 5) · Merkabit_crystal, successor of GAPMAP
(menu item 2). Framework: Goldston–Montgomery / Montgomery–Soundararajan
(variance of ψ in short intervals ↔ pair correlation / form factor of zeros) —
no framework priority claimed. Ours: the object, the exact prime data, and the
conductor-enabled regime.**

## Scope lock (binding)

Exploration. NOT the wall (Op2/GRH(ψ)/R-209 — OPEN, claimed nowhere). Not RH/GRH
(the prediction formulas assume zeros on the ½-line; an off-line zero would show
as anomalous x-growth of the variance — a byproduct check, not a target).
Predictions use the stationary cluster expansion (GUE Y₂ or Poisson) — heuristic
grade; measurement is exact counting. Refutations/nulls at equal prominence.

## The idea (and why our object is uniquely suited)

Var over x of [ψ_χ(xe^δ) − ψ_χ(x)]/√x, computed from primes, equals (explicit
formula + cluster expansion):

  **Var(δ; x) ≈ ∫₀^∞ ρ̄(γ) · |φ_γ(δ)|²·2 · F(α(γ,x)) dγ**,
  |φ_γ|² = (e^δ − 2e^{δ/2}cos γδ + 1)/(¼+γ²),  α = ln x/(2πρ̄(γ)),

with F the Montgomery form factor: **F = min(α, 1) for GUE, F ≡ 1 for
uncorrelated (Poisson) zeros.** The correlation-sensitive regime is α < 1.
For the genome, 2πρ̄(γ) = ln Q − 8 ln 2π + 8 ln γ ≈ 15.7 + 8 ln γ with
ln Q = 30.44 — so **α < 1 for all γ ≳ 2 at x = 10⁸**: the conductor makes the
spectrum dense enough that essentially the whole variance is
correlation-sensitive, and GUE vs Poisson differ by a factor ~2–3. For ζ at the
same (x, δ), α ≈ 3–4 ⇒ F = 1 ⇒ ζ is the pipeline control with zero
discrimination — the contrast is the design.

Sweeping δ ∈ [10⁻⁴, 0.3] and windows [10⁷,10⁸], [10⁸,10⁹] maps the genome's
form factor over α ≈ 0.22–0.7 at heights γ up to ~10⁴ — **the first two-level
statistics of the genome spectrum in the bulk**, far beyond the 153 certified
zeros (γ ≤ 27.9), where BK_OFFDIAG Q2 stalled.

Complementarity with GAPMAP: the gap between the measured variance and the
"uncorrelated-prime-signs" benchmark (δ·⟨a²⟩·ln x) equals the HL-weighted SUM of
per-gap correlations over ~xδ gaps — the integrated quantity GAPMAP's per-gap
nulls (≤10⁻³ each) could not constrain.

## Inputs

- Exact Frobenius classes for all 50,847,534 primes ≤ 10⁹ (BK_OFFDIAG + GAPMAP
  caches) → exact ψ_χ and ψ_ζ at arbitrary x (prime powers included exactly via
  the closed-form tr ρ(Frobᵏ); ramified: Λ_χ(3ᵏ)=ln 3, Λ_χ(7ᵏ)=0).
- `genome_deep_zeros_raw/signal_100000000000.npy` (ψ_χ/√x, 8000-pt grid to 10¹¹)
  → coarse extension window [10⁹,10¹¹], δ ≥ 0.014 (binned; caveat recorded).
- Exact ρ̄ = θ′/π (RUN 73's validated analytic data).

## Questions

- **Q1:** does the measured Var(δ; x) follow the GUE curve or the Poisson curve?
  Money plot: R(δ) = Var_meas/Var_Poisson vs the GUE ratio curve (≈0.3–0.6)
  and 1 (Poisson) — genome panel vs ζ panel (control pinned at 1 by design).
- **Q2:** x-consistency across the three windows (10⁷–10⁸, 10⁸–10⁹, 10⁹–10¹¹) —
  α scales with ln x; the GUE prediction moves accordingly; Poisson does not.
- **Q3:** the integrated-gap-correlation reading: measured minus the
  uncorrelated-primes benchmark, vs the GUE expectation (ties SIVAR to GAPMAP).

## Files

`sivar.py` → `FINDINGS_sivar.md`, `results/sivar_results.json`, `results/*.png`.
