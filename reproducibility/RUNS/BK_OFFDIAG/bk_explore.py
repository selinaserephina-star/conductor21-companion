# bk_explore.py — BK_OFFDIAG exploratory run (2026-09-11, Stenberg + Claude/Fable 5).
# Q1: zeros -> primes resonance, amplitude-resolved (153 chi8 zeros; zeta 599 control).
# Q2: pair-correlation residual vs BK arithmetic reconstruction (shape, fitted scale).
# Q3: twin-Frobenius correlation to 1e8 (independence / blindness test).
# Exploration only. NOT the wall, not RH/GRH, no twin-prime-conjecture claim.
import json, math, sys, time
from pathlib import Path

import numpy as np
import mpmath as mp

HERE = Path(__file__).parent
RES = HERE / "results"
RES.mkdir(exist_ok=True)
sys.path.insert(0, str(HERE / "data"))
from gen_an_chi8_REFERENCE import frob_class, prime_sieve, LOC8, AP  # noqa: E402

mp.mp.dps = 30
LNQ = 10 * math.log(21.0)          # ln(21^10)
X_TWIN = 10**8                     # Q3 prime range
RNG = np.random.default_rng(21)
trap = np.trapezoid

t00 = time.time()
def log(msg):
    print(f"[{time.time()-t00:7.1f}s] {msg}", flush=True)

# ---------------------------------------------------------------- inputs
gam = np.array([float(l.split(",")[1]) for l in
                (HERE / "data/chi8_zeros_153.csv").read_text().splitlines()[1:]])
zz = np.array([float(l) for l in
               (HERE / "data/zeta_zeros_599.txt").read_text().split()])
log(f"chi8 zeros: {len(gam)} (gamma_max={gam[-1]:.3f}); zeta zeros: {len(zz)} (max={zz[-1]:.3f})")

# ------------------------------------------------- local-factor eigenvalues per class
def class_alphas(cls):
    """Frobenius eigenvalues alpha_i: local poly is prod (1 - alpha_i X)."""
    coef = LOC8[cls]                       # low-to-high in X
    if len(coef) == 1:                     # p7: local factor 1 -> no eigenvalues
        return np.array([])
    r = np.roots(list(reversed(coef)))     # roots in X are 1/alpha_i
    return 1.0 / r

ALPH = {c: class_alphas(c) for c in LOC8}
for c in ("1A", "2A", "3A", "4A", "7", "p3", "p7"):
    a = ALPH[c]
    s1 = float(np.real(a.sum())) if len(a) else 0.0
    s2 = float(np.real((a**2).sum())) if len(a) else 0.0
    assert abs(s1 - AP[c]) < 1e-8, (c, s1)
    log(f"  class {c:2}: n_alpha={len(a)}  tr={s1:+.3f}  tr(Frob^2)={s2:+.3f}")

FROB_SMALL = {}
def frob_small(p):
    if p not in FROB_SMALL:
        FROB_SMALL[p] = frob_class(p)
    return FROB_SMALL[p]

def lambda_chi(n):
    """Lambda_chi8(n) = Re(sum alpha_i^k) ln p for n = p^k, else 0."""
    p = None
    for q in range(2, int(n**0.5) + 1):
        if n % q == 0:
            p = q
            break
    if p is None:
        p, k = n, 1
    else:
        k, m = 0, n
        while m % p == 0:
            m //= p
            k += 1
        if m != 1:
            return 0.0, None
    cls = frob_small(p)
    a = ALPH[cls]
    s = float(np.real((a**k).sum())) if len(a) else 0.0
    return s * math.log(p), (p, k, cls)

def lambda_classical(n):
    for p in range(2, int(n**0.5) + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return (math.log(p), None) if m == 1 else (0.0, None)
    return math.log(n), None

# ---------------------------------------------------------------- Q1 machinery
def cos_transform(w, gg, ugrid, chunk=800):
    """Wc(u) = int w(g) cos(u g) dg, chunked to bound memory."""
    out = np.empty(len(ugrid))
    for i in range(0, len(ugrid), chunk):
        u = ugrid[i:i + chunk]
        out[i:i + chunk] = trap(w[None, :] * np.cos(np.outer(u, gg)), gg, axis=1)
    return out

def resonance(zeros, dens, lam_fn, xgrid, tag, ngrid=4000):
    """Windowed resonance F(x) minus smooth main term; explicit-formula prediction
    -(1/2pi) sum_n Lam(n)/sqrt(n) [Wc(x - ln n) + Wc(x + ln n)]."""
    L = zeros[-1] * 1.001
    gg = np.linspace(0, L, ngrid)
    wg = 0.5 * (1 - np.cos(2 * np.pi * gg / L))          # Hann on [0, L]
    wz = 0.5 * (1 - np.cos(2 * np.pi * zeros / L))
    F = (wz[None, :] * np.cos(np.outer(xgrid, zeros))).sum(axis=1)
    rho = dens(gg)
    Fsm = cos_transform(wg * rho, gg, xgrid)
    resid = F - Fsm
    ug = np.linspace(-1.0, xgrid[-1] + 6.5, 8000)
    Wc = cos_transform(wg, gg, ug)
    pred = np.zeros_like(xgrid)
    terms = []
    for n in range(2, 400):
        lam, info = lam_fn(n)
        if lam == 0.0:
            continue
        amp = lam / math.sqrt(n) / (2 * math.pi)
        pred -= amp * (np.interp(xgrid - math.log(n), ug, Wc)
                       + np.interp(xgrid + math.log(n), ug, Wc))
        if n <= 33:
            terms.append((n, lam, info))
    log(f"Q1[{tag}]: resid rms={resid.std():.3f}, pred rms={pred.std():.3f}")
    return resid, pred, terms

def dens_chi8(g):
    out = np.empty_like(g)
    for i, t in enumerate(g):
        out[i] = (0.5 * LNQ - 4 * math.log(2 * math.pi)
                  + 4 * float(mp.re(mp.digamma(0.5 + 1j * t)))) / math.pi
    return np.maximum(out, 0)

def dens_zeta(g):
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.log(np.maximum(g, 1e-9) / (2 * np.pi)) / (2 * np.pi)
    return np.maximum(d, 0)

xg = np.arange(0.30, 3.50, 0.002)
log("Q1: chi8 resonance ...")
resid8, pred8, terms8 = resonance(gam, dens_chi8, lambda_chi, xg, "chi8")
log("Q1: zeta control ...")
residz, predz, termsz = resonance(zz, dens_zeta, lambda_classical, xg, "zeta", ngrid=12000)

def scatter_stats(resid, pred, lam_fn):
    ns = list(range(2, 34))
    meas = np.array([float(np.interp(math.log(n), xg, resid)) for n in ns])
    prd = np.array([float(np.interp(math.log(n), xg, pred)) for n in ns])
    lam = [lam_fn(n)[0] for n in ns]
    hot = np.array([l != 0 for l in lam])
    return {"n": ns, "measured": meas.tolist(), "predicted": prd.tolist(),
            "Lambda": lam, "pearson_r_all": float(np.corrcoef(meas, prd)[0, 1]),
            "rms_meas_hot": float(meas[hot].std()),
            "rms_meas_cold": float(meas[~hot].std())}

sc8 = scatter_stats(resid8, pred8, lambda_chi)
scz = scatter_stats(residz, predz, lambda_classical)
log(f"Q1 amplitude correlation: chi8 r={sc8['pearson_r_all']:+.3f}  zeta r={scz['pearson_r_all']:+.3f}")
log(f"   chi8 rms at hot n: {sc8['rms_meas_hot']:.3f}  vs cold n: {sc8['rms_meas_cold']:.3f}")

def val_at(resid, x): return float(np.interp(x, xg, resid))
disc = {
    "chi8_at_log2": val_at(resid8, math.log(2)), "pred_at_log2": val_at(pred8, math.log(2)),
    "chi8_at_log3": val_at(resid8, math.log(3)), "pred_at_log3": val_at(pred8, math.log(3)),
    "chi8_at_log7": val_at(resid8, math.log(7)), "pred_at_log7": val_at(pred8, math.log(7)),
    "zeta_at_log7": val_at(residz, math.log(7)), "zeta_pred_at_log7": val_at(predz, math.log(7)),
}
log(f"Q1 log7 discriminant: chi8 {disc['chi8_at_log7']:+.3f} (pred {disc['pred_at_log7']:+.3f}, ~silent) | "
    f"zeta {disc['zeta_at_log7']:+.3f} (pred {disc['zeta_pred_at_log7']:+.3f}, sings)")

# ---------------------------------------------------------------- Q2: R2 residual
log("Q2: unfold + pair correlation ...")
def theta_chi8(t):
    return (t / 2) * LNQ - 4 * t * math.log(2 * math.pi) \
        + 4 * float(mp.im(mp.loggamma(0.5 + 1j * t)))
u = np.array([theta_chi8(t) / math.pi + 1 for t in gam])
span = u[-1] - u[0]
log(f"   unfolded span {span:.2f} for N={len(u)} (mean spacing {span/(len(u)-1):.4f})")

BINW, DMAX = 0.25, 6.0
bins = np.arange(0, DMAX + BINW, BINW)
mid = 0.5 * (bins[1:] + bins[:-1])
def pair_hist(x):
    d = np.abs(x[:, None] - x[None, :])[np.triu_indices(len(x), 1)]
    return np.histogram(d, bins=bins)[0]

M = 2000
h_unif = np.zeros((M, len(mid)))
for i in range(M):
    h_unif[i] = pair_hist(np.sort(RNG.uniform(u[0], u[-1], len(u))))
base = h_unif.mean(axis=0)                    # uniform => R2 = 1 baseline (edge-corrected)
R2 = pair_hist(u) / base

NM, NS = 301, 300                             # GUE surrogate band, same estimator
h_gue = np.zeros((NS, len(mid)))
for i in range(NS):
    A = RNG.normal(size=(NM, NM)) + 1j * RNG.normal(size=(NM, NM))
    ev = np.linalg.eigvalsh((A + A.conj().T) / 2)
    Rr = np.abs(ev).max() * 1.0001
    cdf = NM * (0.5 + (ev * np.sqrt(np.maximum(Rr**2 - ev**2, 0)) / Rr**2
                       + np.arcsin(np.clip(ev / Rr, -1, 1))) / math.pi)
    sel = np.sort(np.argsort(np.abs(ev))[:len(u)])   # central 153 levels
    h_gue[i] = pair_hist(np.sort(cdf[sel]))
R2_gue = h_gue / base
gue_lo, gue_md, gue_hi = np.percentile(R2_gue, [2.5, 50, 97.5], axis=0)
sig = np.maximum((gue_hi - gue_lo) / 3.92, 1e-9)

sinc2 = 1 - np.sinc(mid)**2
res_data = R2 - sinc2
rho_eff = (len(u) - 1) / (gam[-1] - gam[0])
bk = np.zeros_like(mid)
for n in range(2, 400):
    lam, _ = lambda_chi(n)
    if lam:
        bk += lam**2 / n * np.cos(mid / rho_eff * math.log(n))
bk_fit = float(np.dot(res_data / sig**2, bk) / np.dot(bk / sig, bk / sig))
res_pred = bk_fit * bk
chi2_flat = float(((res_data / sig)**2).sum())
chi2_bk = float((((res_data - res_pred) / sig)**2).sum())
log(f"Q2: BK fitted scale {bk_fit:+.5f}; chi2 flat {chi2_flat:.1f} vs +BK {chi2_bk:.1f} (bins {len(mid)})")

# ---------------------------------------------------------------- Q3: twin-Frobenius
log(f"Q3: sieving to {X_TWIN:.0e} + exact Frobenius classes ...")
ps = prime_sieve(X_TWIN)
log(f"   {len(ps)} primes")
CLS_ID = {"1A": 0, "2A": 1, "3A": 2, "4A": 3, "7": 4, "p3": 5, "p7": 6}
cache = RES / f"classes_{X_TWIN}.npz"
if cache.exists():
    dat = np.load(cache)
    cls_id, ap = dat["cls_id"], dat["ap"]
    log("   (classes loaded from cache)")
else:
    cls_id = np.empty(len(ps), dtype=np.int8)
    ap = np.empty(len(ps), dtype=np.int8)
    for i, p in enumerate(ps):
        c = frob_class(int(p))
        cls_id[i] = CLS_ID[c]
        ap[i] = AP[c]
        if i and i % 1000000 == 0:
            log(f"   ... {i}/{len(ps)}")
    np.savez_compressed(cache, cls_id=cls_id, ap=ap)
frac = {c: float((cls_id == i).mean()) for c, i in CLS_ID.items() if i < 5}
log(f"   Chebotarev: { {k: round(v,5) for k,v in frac.items()} }"
    " (expect 1A .00595 2A .125 3A .33333 4A .25 7 .28571)")

def gap_corr(gap):
    """pairs (p, p+gap) both prime, p > 7; correlation + joint class table."""
    q = ps + gap
    j = np.searchsorted(ps, q)
    j = np.minimum(j, len(ps) - 1)
    ok = (ps[j] == q) & (ps > 7)
    i1 = np.nonzero(ok)[0]
    j1 = j[ok]
    a1, a2 = ap[i1].astype(float), ap[j1].astype(float)
    c1, c2 = cls_id[i1], cls_id[j1]
    N = len(a1)
    prod = a1 * a2
    corr, se = float(prod.mean()), float(prod.std() / math.sqrt(N))
    null = float(a1.mean() * a2.mean())
    O = np.zeros((5, 5))
    np.add.at(O, (c1, c2), 1)
    m1, m2 = O.sum(1) / N, O.sum(0) / N
    E = np.outer(m1, m2) * N
    mask = E > 0
    chi2 = float(((O - E)[mask]**2 / E[mask]).sum())
    dof = int(mask.sum()) - 2 * (5 - 1) - 1
    return {"gap": gap, "N_pairs": N, "corr_aa": corr, "se": se,
            "null_prod_means": null, "z_vs_null": (corr - null) / se,
            "mean_a_left": float(a1.mean()), "mean_a_right": float(a2.mean()),
            "chi2_joint": chi2, "dof": dof, "table_O": O.tolist(),
            "marg1": m1.tolist(), "marg2": m2.tolist()}

results_gaps = {}
for gap, name in [(2, "twin"), (4, "cousin"), (6, "sexy")]:
    g = gap_corr(gap)
    results_gaps[name] = g
    log(f"Q3 {name} (p,p+{gap}): N={g['N_pairs']}  <a a'>={g['corr_aa']:+.5f} +- {g['se']:.5f}  "
        f"z={g['z_vs_null']:+.2f}  chi2_joint={g['chi2_joint']:.1f} / dof {g['dof']}")

i0 = np.nonzero(ps[:-1] > 7)[0]
aa = ap[:-1][i0].astype(float) * ap[1:][i0].astype(float)
cons = {"N_pairs": int(len(aa)), "corr_aa": float(aa.mean()),
        "se": float(aa.std() / math.sqrt(len(aa)))}
log(f"Q3 consecutive-prime control: <a a'>={cons['corr_aa']:+.5f} +- {cons['se']:.5f} "
    f"(z={cons['corr_aa']/cons['se']:+.2f})")

C2 = 0.6601618158468696
tt = np.linspace(5, X_TWIN, 200000)
li2 = trap(1 / np.log(tt)**2, tt)
pi2 = results_gaps["twin"]["N_pairs"] + 2      # add (3,5),(5,7) excluded by p>7
hl = {"pi2_measured": int(pi2), "pi2_HL_pred": float(2 * C2 * li2),
      "ratio": float(pi2 / (2 * C2 * li2))}
log(f"Q3 positive control: pi2(1e8) = {hl['pi2_measured']} vs HL {hl['pi2_HL_pred']:.0f} "
    f"(ratio {hl['ratio']:.4f})")

# ---------------------------------------------------------------- outputs
out = {
    "meta": {"date": "2026-09-11", "zeros": len(gam), "gamma_max": float(gam[-1]),
             "zeta_control": len(zz), "X_twin": X_TWIN,
             "scope": "exploration; not the wall; not RH/GRH; no TPC claim"},
    "Q1": {"chi8_scatter": sc8, "zeta_scatter": scz, "discriminants": disc,
           "chi8_terms_le33": [(n, l, str(i)) for n, l, i in terms8]},
    "Q2": {"binw": BINW, "mid": mid.tolist(), "R2": R2.tolist(),
           "gue_band": [gue_lo.tolist(), gue_md.tolist(), gue_hi.tolist()],
           "bk_fitted_scale": bk_fit, "chi2_flat": chi2_flat, "chi2_bk": chi2_bk,
           "rho_eff": float(rho_eff)},
    "Q3": {name: {k: v for k, v in g.items() if k != "table_O"}
           for name, g in results_gaps.items()} |
          {"twin_table": results_gaps["twin"]["table_O"], "consecutive": cons,
           "chebotarev_frac": frac, "hl_positive_control": hl},
}
(RES / "bk_offdiag_results.json").write_text(json.dumps(out, indent=1))
np.savez(RES / "q1_curves.npz", xg=xg, resid8=resid8, pred8=pred8, residz=residz, predz=predz)
log("saved results/bk_offdiag_results.json + q1_curves.npz")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pp = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32]
    fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for A, (resid, pred, ttl) in zip(ax, [(resid8, pred8, "chi8 genome (153 certified zeros)"),
                                          (residz, predz, "zeta control (599 zeros)")]):
        A.plot(xg, resid, lw=1, label="measured residual")
        A.plot(xg, pred, lw=1, ls="--", label="explicit-formula prediction")
        for n in pp:
            A.axvline(math.log(n), color="gray", alpha=.25, lw=.6)
        A.set_title(ttl); A.legend(fontsize=8); A.set_ylabel("F - smooth")
    ax[1].set_xlabel("x   (gray: log p^k;  log 7 = 1.946 — chi8 predicted silent there)")
    fig.tight_layout(); fig.savefig(RES / "q1_resonance.png", dpi=130); plt.close(fig)

    fig, A = plt.subplots(figsize=(6, 5))
    A.scatter(sc8["predicted"], sc8["measured"], c="tab:red",
              label=f"chi8  r={sc8['pearson_r_all']:+.2f}")
    A.scatter(scz["predicted"], scz["measured"], c="tab:blue", marker="s", s=18,
              label=f"zeta  r={scz['pearson_r_all']:+.2f}")
    lim = max(abs(v) for v in A.get_xlim() + A.get_ylim())
    A.plot([-lim, lim], [-lim, lim], "k:", lw=.8)
    A.set_xlabel("predicted amplitude at log n"); A.set_ylabel("measured")
    A.legend(); fig.tight_layout(); fig.savefig(RES / "q1_scatter.png", dpi=130); plt.close(fig)

    fig, A = plt.subplots(figsize=(9, 5))
    A.fill_between(mid, gue_lo, gue_hi, alpha=.25, label="GUE 95% band (N=153, same estimator)")
    A.plot(mid, gue_md, lw=.8, ls=":", color="k")
    A.plot(mid, R2, "o-", ms=3, lw=1, label="chi8 R2 (unfolded, 153 zeros)")
    A.plot(mid, sinc2, lw=1, label="GUE  1 - sinc^2")
    A.plot(mid, sinc2 + res_pred, lw=1, ls="--", label=f"+ BK shape (fitted {bk_fit:+.4f})")
    A.set_xlabel("unfolded distance s"); A.set_ylabel("R2"); A.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(RES / "q2_r2.png", dpi=130); plt.close(fig)
    log("plots saved")
except Exception as e:  # noqa: BLE001
    log(f"plotting skipped: {e}")
log("DONE")
