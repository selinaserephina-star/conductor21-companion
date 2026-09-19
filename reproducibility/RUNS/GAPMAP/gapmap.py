# gapmap.py — systematic gap-blindness map, even gaps g<=60, X=1e8 then 1e9.
# 2026-09-11, Stenberg + Claude (Fable 5). Exploration; not the wall; no TPC claim.
import json, math, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
RES = HERE / "results"
RES.mkdir(exist_ok=True)
BK = HERE.parent / "BK_OFFDIAG_2026-09-11"
sys.path.insert(0, str(BK / "data"))
from gen_an_chi8_REFERENCE import frob_class, prime_sieve, AP  # noqa: E402

GAPS = list(range(2, 61, 2))
CLS_ID = {"1A": 0, "2A": 1, "3A": 2, "4A": 3, "7": 4, "p3": 5, "p7": 6}
AP_BY_ID = np.array([8, 0, -1, 0, 1, 1, 0], dtype=np.int8)   # a_p per class id
C2 = 0.6601618158468696
t00 = time.time()
def log(m): print(f"[{time.time()-t00:7.1f}s] {m}", flush=True)

# ---------------------------------------------------------------- stage A data
dat = np.load(BK / "results" / "classes_100000000.npz")
cls8, ap8 = dat["cls_id"], dat["ap"]
ps8 = prime_sieve(10**8)
assert len(ps8) == len(cls8)
log(f"stage A: {len(ps8)} primes <= 1e8 (cached classes)")

# ---------------------------------------------------------------- stage B cache
CACHE9 = RES / "classes_1e8_1e9.npz"
def segmented_primes(lo, hi, seg=10**7):
    """primes in (lo, hi], RAM-light."""
    base = prime_sieve(int(math.isqrt(hi)) + 1)
    out = []
    start = lo + 1
    while start <= hi:
        end = min(start + seg - 1, hi)
        sieve = np.ones(end - start + 1, dtype=bool)
        for p in base:
            p = int(p)
            if p * p > end:
                break
            first = max(p * p, ((start + p - 1) // p) * p)
            sieve[first - start::p] = False
        if start <= 1:
            sieve[:2 - start] = False
        out.append(np.nonzero(sieve)[0] + start)
        start = end + 1
    return np.concatenate(out)

PART9 = RES / "classes_1e8_1e9_PARTIAL.npy"
if CACHE9.exists():
    d9 = np.load(CACHE9)
    ps_hi, cls_hi = d9["ps"], d9["cls_id"]
    log(f"stage B: loaded cache, {len(ps_hi)} primes in (1e8,1e9]")
else:
    log("stage B: segmented sieve (1e8,1e9] ...")
    ps_hi = segmented_primes(10**8, 10**9)
    log(f"   {len(ps_hi)} new primes; Frobenius classes (~17 min, checkpointed) ...")
    cls_hi = np.full(len(ps_hi), -1, dtype=np.int8)
    i0 = 0
    if PART9.exists():                       # resume from checkpoint
        part = np.load(PART9)
        if len(part) == len(ps_hi):
            cls_hi = part
            i0 = int(np.argmax(cls_hi < 0)) if (cls_hi < 0).any() else len(ps_hi)
            log(f"   resuming at {i0}/{len(ps_hi)}")
    for i in range(i0, len(ps_hi)):
        cls_hi[i] = CLS_ID[frob_class(int(ps_hi[i]))]
        if i and i % 2000000 == 0:
            np.save(PART9, cls_hi)
            log(f"   ... {i}/{len(ps_hi)} (checkpoint)")
    np.savez_compressed(CACHE9, ps=ps_hi.astype(np.int64), cls_id=cls_hi)
    if PART9.exists():
        PART9.unlink()
    log("   cached")

ps = np.concatenate([ps8.astype(np.int64), ps_hi.astype(np.int64)])
cls = np.concatenate([cls8, cls_hi])
ap = AP_BY_ID[cls]
log(f"total: {len(ps)} primes <= 1e9")
fr = {c: float((cls == i).mean()) for c, i in CLS_ID.items() if i < 5}
log(f"Chebotarev at 1e9: { {k: round(v,5) for k,v in fr.items()} }")

# ---------------------------------------------------------------- gap machinery
def li2(x, n=400000):
    t = np.linspace(2.0, x, n)
    return float(np.trapezoid(1 / np.log(t)**2, t))

def sing(g):
    c = 1.0
    m = g
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29):
        if m % p == 0:
            c *= (p - 1) / (p - 2)
            while m % p == 0:
                m //= p
    return c

def gap_stats(g, lo_p):
    """pairs (p, p+g) both prime, lo_p < p <= max; exact stats."""
    q = ps + g
    j = np.searchsorted(ps, q)
    j[j >= len(ps)] = len(ps) - 1
    ok = (ps[j] == q) & (ps > lo_p)
    i1 = np.nonzero(ok)[0]
    j1 = j[i1]
    a1, a2 = ap[i1].astype(np.float64), ap[j1].astype(np.float64)
    c1, c2 = cls[i1], cls[j1]
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
    return {"gap": g, "N": N, "corr": corr, "se": se, "null": null,
            "z": (corr - null) / se, "chi2": chi2, "dof": int(mask.sum()) - 9}

def run_map(lo_p, X, tag):
    LI2 = li2(X) - (li2(lo_p) if lo_p > 7 else 0.0)
    rows = []
    for g in GAPS:
        r = gap_stats(g, max(lo_p, 7))
        r["N_HL"] = 2 * C2 * sing(g) * LI2
        r["ratio_HL"] = r["N"] / r["N_HL"]
        rows.append(r)
    zs = np.array([r["z"] for r in rows])
    chi2s = np.array([r["chi2"] for r in rows])
    dofs = np.array([r["dof"] for r in rows])
    glob = {
        "tag": tag, "n_gaps": len(GAPS), "max_abs_z": float(np.abs(zs).max()),
        "argmax_gap": int(GAPS[int(np.abs(zs).argmax())]),
        "n_absz_gt2": int((np.abs(zs) > 2).sum()),
        "bonf_threshold_5pct": 2.93,
        "n_over_bonf": int((np.abs(zs) > 2.93).sum()),
        "z_mean": float(zs.mean()), "z_std": float(zs.std()),
        "chi2_total": float(chi2s.sum()), "dof_total": int(dofs.sum()),
        "HL_ratio_mean": float(np.mean([r["ratio_HL"] for r in rows])),
        "HL_ratio_worst": float(max(abs(r["ratio_HL"] - 1) for r in rows)),
    }
    log(f"[{tag}] max|z|={glob['max_abs_z']:.2f} @g={glob['argmax_gap']}; "
        f"|z|>2: {glob['n_absz_gt2']}/30; >Bonf(2.93): {glob['n_over_bonf']}; "
        f"z std={glob['z_std']:.2f}; chi2 {glob['chi2_total']:.0f}/{glob['dof_total']}; "
        f"HL ratio mean {glob['HL_ratio_mean']:.4f} worst dev {glob['HL_ratio_worst']:.4f}")
    return rows, glob

log("map at X=1e8 (full cache range) ...")
ps_all, cls_all, ap_all = ps, cls, ap                      # keep refs
ps, cls, ap = ps8.astype(np.int64), cls8, AP_BY_ID[cls8]
rows8, glob8 = run_map(7, 10**8, "1e8")
c8 = [r for r in rows8 if r["gap"] == 4][0]
log(f"   cousin g=4 at 1e8: z={c8['z']:+.2f} (the BK_OFFDIAG flag, same data)")

log("map at X=1e9 (full) ...")
ps, cls, ap = ps_all, cls_all, ap_all
rows9, glob9 = run_map(7, 10**9, "1e9")

log("DISJOINT decade (1e8,1e9] — pre-registered cousin retest + full map ...")
rowsD, globD = run_map(10**8, 10**9, "disjoint")
cD = [r for r in rowsD if r["gap"] == 4][0]
verdict = "SURVIVES (bug-hunt next)" if (cD["z"] >= 2.8) else "DOES NOT SURVIVE — fluctuation, retired"
log(f"   COUSIN VERDICT g=4 disjoint decade: z={cD['z']:+.2f}, N={cD['N']} -> {verdict}")

out = {
    "meta": {"date": "2026-09-11", "gaps": GAPS,
             "scope": "exploration; exact counting; independence null; "
                      "not the wall; not RH/GRH; no TPC claim"},
    "chebotarev_1e9": fr,
    "map_1e8": {"rows": rows8, "global": glob8},
    "map_1e9": {"rows": rows9, "global": glob9},
    "map_disjoint": {"rows": rowsD, "global": globD},
    "cousin": {"z_1e8": c8["z"], "z_disjoint": cD["z"], "N_disjoint": cD["N"],
               "verdict": verdict},
}
(RES / "gapmap_results.json").write_text(json.dumps(out, indent=1))
log("saved results/gapmap_results.json")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    g = np.array(GAPS)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True,
                                 gridspec_kw={"height_ratios": [1, 1]})
    for rows, colr, lab in [(rows8, "tab:blue", "X=1e8"), (rows9, "tab:red", "X=1e9")]:
        cc = [r["corr"] for r in rows]; ee = [r["se"] for r in rows]
        a1.errorbar(g + (0.3 if colr == "tab:red" else -0.3), cc, yerr=ee, fmt="o",
                    ms=4, color=colr, label=f"genome <a_p a_(p+g)> {lab}", capsize=2)
    a1.axhline(0, color="k", lw=.8)
    a1.set_ylabel("Frobenius correlation")
    a1.set_title("Gap-blindness map: the genome channel is flat while the classical "
                 "gap structure is rich")
    a1.legend(fontsize=8)
    r9 = [r["ratio_HL"] for r in rows9]
    a2.plot(g, [sing(int(x)) for x in g], "s--", ms=4, color="tab:green",
            label="Hardy-Littlewood singular series  prod (p-1)/(p-2), p|g odd")
    a2.plot(g, np.array(r9) * np.array([sing(int(x)) for x in g]), "o", ms=4,
            color="tab:purple", label="measured N_g / (2 C2 Li2)  (X=1e9)")
    a2.set_xlabel("gap g"); a2.set_ylabel("relative pair abundance")
    a2.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(RES / "gapmap.png", dpi=130); plt.close(fig)

    fig, A = plt.subplots(figsize=(6, 5))
    for rows, colr, lab in [(rows8, "tab:blue", "1e8"), (rowsD, "tab:orange", "disjoint decade")]:
        zs = np.sort([r["z"] for r in rows])
        n = len(zs)
        # normal quantiles via inverse erf (numpy has erfinv in special? use approximation)
        pgrid = (np.arange(n) + 0.5) / n
        # Acklam-style inverse normal CDF, vectorized, adequate for a QQ plot
        def invnorm(p):
            a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
            b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                 6.680131188771972e+01, -1.328068155288572e+01]
            c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                 -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
            d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                 3.754408661907416e+00]
            p = np.asarray(p, dtype=float)
            x = np.empty_like(p)
            plow = p < 0.02425
            phigh = p > 1 - 0.02425
            mid = ~(plow | phigh)
            q = np.sqrt(-2 * np.log(p[plow]))
            x[plow] = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                      ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
            q = np.sqrt(-2 * np.log(1 - p[phigh]))
            x[phigh] = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                       ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
            q = p[mid] - 0.5; r = q * q
            x[mid] = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
                     (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
            return x
        A.plot(invnorm(pgrid), zs, "o", ms=4, color=colr, label=f"gap z-scores, {lab}")
    lim = 3.2
    A.plot([-lim, lim], [-lim, lim], "k:", lw=.8)
    A.set_xlabel("normal quantile"); A.set_ylabel("observed z")
    A.legend(fontsize=8); A.set_title("Gap-map z ensemble vs pure noise")
    fig.tight_layout(); fig.savefig(RES / "gapmap_qq.png", dpi=130); plt.close(fig)
    log("plots saved")
except Exception as e:  # noqa: BLE001
    log(f"plotting skipped: {e}")
log("DONE")
