# sivar.py — short-interval variance of psi_chi from exact primes to 1e9 (+1e11 signal):
# the genome's bulk form factor, measured. 2026-09-11, Stenberg + Claude (Fable 5).
# Exploration; not the wall; not RH/GRH.
import gc, json, math, sys, time
from pathlib import Path

import numpy as np
import mpmath as mp

HERE = Path(__file__).parent
RES = HERE / "results"
RES.mkdir(exist_ok=True)
BK = HERE.parent / "BK_OFFDIAG_2026-09-11"
GM = HERE.parent / "GAPMAP_2026-09-11"
LANE = Path(r"C:\Users\selin\OneDrive\Desktop\Ilya Riemanns")
sys.path.insert(0, str(BK / "data"))
from gen_an_chi8_REFERENCE import prime_sieve  # noqa: E402

mp.mp.dps = 25
LNQ = 10 * math.log(21.0)
AP_BY_ID = np.array([8, 0, -1, 0, 1, 1, 0], dtype=np.int8)  # 1A,2A,3A,4A,7,p3,p7
ID2CLS = ["1A", "2A", "3A", "4A", "7", "p3", "p7"]
t00 = time.time()
def log(m): print(f"[{time.time()-t00:7.1f}s] {m}", flush=True)

def chi8_pow(c, k):     # tr rho(Frob^k) closed form (deep_sweep convention, verified)
    if c == "1A": return 8
    if c == "2A": return 8 if k % 2 == 0 else 0
    if c == "3A": return 8 if k % 3 == 0 else -1
    if c == "4A": return 8 if k % 4 == 0 else 0
    if c == "7":  return 8 if k % 7 == 0 else 1
    if c == "p3": return 1
    return 0            # p7

# ---------------------------------------------------------------- exact psi engines
log("loading prime/class caches to 1e9 ...")
ps8 = prime_sieve(10**8).astype(np.int64)
cls8 = np.load(BK / "results" / "classes_100000000.npz")["cls_id"]
d9 = np.load(GM / "results" / "classes_1e8_1e9.npz")
ps = np.concatenate([ps8, d9["ps"]])
cls = np.concatenate([cls8, d9["cls_id"]])
del ps8, cls8, d9; gc.collect()
log(f"   {len(ps)} primes")

lp = np.log(ps.astype(np.float64))
cum_chi = np.cumsum(AP_BY_ID[cls].astype(np.float64) * lp).astype(np.float32)
cum2_chi = np.cumsum((AP_BY_ID[cls].astype(np.float64) * lp)**2)
cum_z = np.cumsum(lp)
del lp; gc.collect()
log("   cumulative psi arrays built (chi f32, zeta/L2 f64)")

# prime powers p^k<=1e9 (k>=2): small separate event lists
pow_n, pow_wc, pow_wz = [], [], []
nsmall = int(np.searchsorted(ps, 31623, side="right"))
for i in range(nsmall):
    p = int(ps[i]); c = ID2CLS[cls[i]]; lpp = math.log(p)
    pk, k = p * p, 2
    while pk <= 10**9:
        pow_n.append(pk); pow_wc.append(chi8_pow(c, k) * lpp); pow_wz.append(lpp)
        pk *= p; k += 1
o = np.argsort(pow_n)
pow_n = np.array(pow_n, dtype=np.int64)[o]
pow_cum_c = np.cumsum(np.array(pow_wc)[o])
pow_cum_z = np.cumsum(np.array(pow_wz)[o])
log(f"   {len(pow_n)} prime-power events added")

def psi_at(x, cum, pcum):
    i = np.searchsorted(ps, x, side="right")
    j = np.searchsorted(pow_n, x, side="right")
    out = np.where(i > 0, cum[np.maximum(i - 1, 0)], 0.0)
    return out + np.where(j > 0, pcum[np.maximum(j - 1, 0)], 0.0)

# ---------------------------------------------------------------- density / predictions
gz_master = np.geomspace(0.05, 5e6, 3000)
def theta_p(T):
    if T < 300:
        return 0.5 * LNQ - 4 * math.log(2 * math.pi) + 4 * float(mp.re(mp.digamma(0.5 + 1j * T)))
    return 0.5 * LNQ - 4 * math.log(2 * math.pi) + 4 * math.log(T)
rho_chi_m = np.array([theta_p(g) / math.pi for g in gz_master])
def rho_chi(g): return np.interp(g, gz_master, rho_chi_m)
def rho_zeta(g):
    return np.maximum(np.log(np.maximum(g, 1e-6) / (2 * math.pi)) / (2 * math.pi), 0.0)

def var_pred(delta, u_lo, u_hi, obj, model):
    g = np.geomspace(0.05, max(400 / delta, 50), 6000)
    phi2 = (math.exp(delta) - 2 * math.exp(delta / 2) * np.cos(g * delta) + 1) / (0.25 + g**2)
    rho = rho_chi(g) if obj == "chi" else rho_zeta(g)
    tot = 0.0
    uu = np.linspace(u_lo, u_hi, 7)
    for u in uu:
        F = np.clip(u / np.maximum(2 * math.pi * rho, 1e-9), 0, 1) if model == "gue" \
            else np.ones_like(g)
        tot += float(np.trapezoid(rho * 2 * phi2 * F, g))
    return tot / len(uu)

# ---------------------------------------------------------------- measurement
DELTAS = np.geomspace(1e-4, 0.3, 26)
WINDOWS = [("W1", math.log(1e7), math.log(1e8)), ("W2", math.log(1e8), math.log(1e9))]
NB = 25          # bootstrap blocks

def measure(obj, u_lo, u_hi, delta):
    n = 20000
    u = np.linspace(u_lo, u_hi - delta, n)
    x1, x2 = np.exp(u), np.exp(u + delta)
    if obj == "chi":
        D = (psi_at(x2, cum_chi, pow_cum_c) - psi_at(x1, cum_chi, pow_cum_c)) / np.sqrt(x1)
    else:
        D = (psi_at(x2, cum_z, pow_cum_z) - psi_at(x1, cum_z, pow_cum_z) - (x2 - x1)) / np.sqrt(x1)
    D = D - D.mean()
    v = float((D**2).mean())
    bl = np.array_split(D, NB)
    bv = np.array([float((b**2).mean()) for b in bl])
    return v, float(bv.std() / math.sqrt(NB))

def bench_chi(u_lo, u_hi, delta):     # exact independent-signs benchmark: sum Lambda^2 in window
    n = 4000
    u = np.linspace(u_lo, u_hi - delta, n)
    x1, x2 = np.exp(u), np.exp(u + delta)
    i1 = np.searchsorted(ps, x1, side="right"); i2 = np.searchsorted(ps, x2, side="right")
    s = np.where(i2 > 0, cum2_chi[np.maximum(i2 - 1, 0)], 0) - \
        np.where(i1 > 0, cum2_chi[np.maximum(i1 - 1, 0)], 0)
    return float((s / x1).mean())

rows = []
for wname, ulo, uhi in WINDOWS:
    for obj in ("chi", "zeta"):
        for dl in DELTAS:
            v, se = measure(obj, ulo, uhi, float(dl))
            r = {"window": wname, "obj": obj, "delta": float(dl), "var": v, "se": se,
                 "pred_gue": var_pred(float(dl), ulo, uhi, obj, "gue"),
                 "pred_poisson": var_pred(float(dl), ulo, uhi, obj, "poisson")}
            if obj == "chi":
                r["bench_signs"] = bench_chi(ulo, uhi, float(dl))
            rows.append(r)
    log(f"{wname} done")

# ---- W3 from the stored 1e11 signal (coarse, genome only) ----
W3rows = []
sig_path = LANE / "genome_deep_zeros_raw" / "signal_100000000000.npy"
try:
    U, gsig = np.load(sig_path)
    du3 = float(U[1] - U[0])
    lo3, hi3 = math.log(1e9), float(U[-1])
    i_lo = int(np.searchsorted(U, lo3))
    for m in (5, 8, 13, 21, 34, 55, 89, 112):
        dl = m * du3
        seg = slice(i_lo, len(U) - m)
        D = gsig[i_lo + m:] * math.exp(dl / 2) - gsig[seg]
        D = D - D.mean()
        v = float((D**2).mean())
        bl = np.array_split(D, NB)
        se = float(np.array([(b**2).mean() for b in bl]).std() / math.sqrt(NB))
        W3rows.append({"window": "W3(signal 1e9-1e11)", "obj": "chi", "delta": dl,
                       "var": v, "se": se,
                       "pred_gue": var_pred(dl, lo3, hi3, "chi", "gue"),
                       "pred_poisson": var_pred(dl, lo3, hi3, "chi", "poisson")})
    log("W3 (stored 1e11 signal) done")
except Exception as e:  # noqa: BLE001
    log(f"W3 skipped: {e}")

# ---- gate A: two independent prime engines agree (flint cache vs 2026-07 numba signal) ----
gate = {}
try:
    U10, g10 = np.load(LANE / "genome_deep_zeros_raw" / "signal_10000000000.npy")
    du = float(U10[1] - U10[0]); m = 19; dl = m * du
    i0, i1 = int(np.searchsorted(U10, math.log(2e8))), int(np.searchsorted(U10, math.log(9e8)))
    idx = np.arange(i0, i1)
    D_sig = g10[idx + m] * math.exp(dl / 2) - g10[idx]
    x1 = np.exp(U10[idx]); x2 = np.exp(U10[idx] + dl)
    D_ex = (psi_at(x2, cum_chi, pow_cum_c) - psi_at(x1, cum_chi, pow_cum_c)) / np.sqrt(x1)
    D_sig = D_sig - D_sig.mean(); D_ex = D_ex - D_ex.mean()
    gate = {"delta": dl, "n": len(idx),
            "corr": float(np.corrcoef(D_sig, D_ex)[0, 1]),
            "rms_ratio": float(D_sig.std() / D_ex.std())}
    log(f"gate A (numba-2026-07 vs flint-cache increments, delta={dl:.4f}): "
        f"corr={gate['corr']:+.4f}, rms ratio={gate['rms_ratio']:.4f}")
except Exception as e:  # noqa: BLE001
    log(f"gate A skipped: {e}")

# ---------------------------------------------------------------- verdicts
def chi2_against(rows_sel, key):
    return sum(((r["var"] - r[key]) / r["se"])**2 for r in rows_sel), len(rows_sel)

summary = {}
for wname in ("W1", "W2"):
    sel = [r for r in rows if r["window"] == wname and r["obj"] == "chi"]
    cg, n = chi2_against(sel, "pred_gue")
    cp, _ = chi2_against(sel, "pred_poisson")
    rat = np.mean([r["var"] / r["pred_poisson"] for r in sel])
    ratg = np.mean([r["pred_gue"] / r["pred_poisson"] for r in sel])
    summary[wname] = {"chi2_gue": cg, "chi2_poisson": cp, "n": n,
                      "mean_R_measured": float(rat), "mean_R_gue": float(ratg)}
    log(f"[{wname} chi] chi2 vs GUE {cg:.0f} vs POISSON {cp:.0f} (n={n}); "
        f"mean R=Var/Poisson: measured {rat:.3f}, GUE predicts {ratg:.3f}")
    selz = [r for r in rows if r["window"] == wname and r["obj"] == "zeta"]
    cz, nz = chi2_against(selz, "pred_gue")   # for zeta gue==poisson effectively
    ratz = np.mean([r["var"] / r["pred_poisson"] for r in selz])
    summary[wname + "_zeta"] = {"chi2_pred": cz, "n": nz, "mean_R": float(ratz)}
    log(f"[{wname} zeta] chi2 vs prediction {cz:.0f} (n={nz}); mean Var/pred {ratz:.3f}")

out = {"meta": {"date": "2026-09-11", "deltas": DELTAS.tolist(),
                "scope": "exploration; cluster-expansion predictions (GUE vs Poisson); "
                         "measurement exact; not the wall; not RH/GRH"},
       "gateA": gate, "rows": rows, "W3": W3rows, "summary": summary}
(RES / "sivar_results.json").write_text(json.dumps(out, indent=1))
log("saved results/sivar_results.json")

# ---------------------------------------------------------------- plots
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    for A, obj, ttl in [(ax[0], "chi", "genome chi8 — the conductor puts ALL heights in the "
                                       "correlation regime"),
                        (ax[1], "zeta", "zeta control — same (x, delta): diagonal regime, "
                                        "no discrimination")]:
        for wname, colr in [("W1", "tab:blue"), ("W2", "tab:red")]:
            sel = [r for r in rows if r["window"] == wname and r["obj"] == obj]
            d = np.array([r["delta"] for r in sel])
            A.errorbar(d, [r["var"] / r["delta"] for r in sel],
                       yerr=[r["se"] / r["delta"] for r in sel], fmt="o", ms=4,
                       color=colr, label=f"measured {wname}")
            A.plot(d, [r["pred_gue"] / r["delta"] for r in sel], "-", lw=1.2, color=colr,
                   label=f"GUE pred {wname}")
            A.plot(d, [r["pred_poisson"] / r["delta"] for r in sel], "--", lw=1, color=colr,
                   alpha=.6, label=f"Poisson pred {wname}")
            if obj == "chi":
                A.plot(d, [r["bench_signs"] / r["delta"] for r in sel], ":", lw=1,
                       color=colr, alpha=.8)
        A.set_xscale("log"); A.set_xlabel("delta (interval length in log x)")
        A.set_ylabel("Var / delta"); A.set_title(ttl, fontsize=9); A.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(RES / "sivar_var.png", dpi=130); plt.close(fig)

    fig, A = plt.subplots(figsize=(8, 5))
    for wname, colr, rowset in [("W1", "tab:blue", rows), ("W2", "tab:red", rows),
                                ("W3(signal 1e9-1e11)", "tab:green", W3rows)]:
        sel = [r for r in rowset if r["window"] == wname and r["obj"] == "chi"]
        if not sel:
            continue
        d = np.array([r["delta"] for r in sel])
        A.errorbar(d, [r["var"] / r["pred_poisson"] for r in sel],
                   yerr=[r["se"] / r["pred_poisson"] for r in sel], fmt="o", ms=4,
                   color=colr, label=f"measured / Poisson  {wname}")
        A.plot(d, [r["pred_gue"] / r["pred_poisson"] for r in sel], "-", lw=1.2, color=colr,
               label=f"GUE ratio {wname}")
    A.axhline(1.0, color="k", lw=.8, ls="--", label="Poisson (uncorrelated zeros)")
    A.set_xscale("log"); A.set_xlabel("delta"); A.set_ylabel("R = Var / Var_Poisson")
    A.set_title("The genome's bulk form factor: measured variance ratio vs GUE ramp")
    A.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(RES / "sivar_ratio.png", dpi=130); plt.close(fig)
    log("plots saved")
except Exception as e:  # noqa: BLE001
    log(f"plotting skipped: {e}")
log("DONE")
