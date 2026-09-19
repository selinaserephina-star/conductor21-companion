# finale.py — four closing micro-probes: A races@1e9, B Gaussianity, C multiplicative
# pairs, D spectral-purity scaling. 2026-09-11, Stenberg + Claude (Fable 5).
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
CR = HERE.parent / "CHEBRACE_2026-09-11"
LANE = Path(r"C:\Users\selin\OneDrive\Desktop\Ilya Riemanns")
sys.path.insert(0, str(BK / "data"))
from gen_an_chi8_REFERENCE import prime_sieve  # noqa: E402

mp.mp.dps = 25
LNQ = 10 * math.log(21.0)
AP_BY_ID = np.array([8, 0, -1, 0, 1, 1, 0], dtype=np.int8)
ID2CLS = ["1A", "2A", "3A", "4A", "7", "p3", "p7"]
C2 = 0.6601618158468696
t00 = time.time()
def log(m): print(f"[{time.time()-t00:7.1f}s] {m}", flush=True)

def chi8_pow(c, k):
    if c == "1A": return 8
    if c == "2A": return 8 if k % 2 == 0 else 0
    if c == "3A": return 8 if k % 3 == 0 else -1
    if c == "4A": return 8 if k % 4 == 0 else 0
    if c == "7":  return 8 if k % 7 == 0 else 1
    return 1 if c == "p3" else 0

log("loading caches ...")
ps8 = prime_sieve(10**8).astype(np.int64)
cls8 = np.load(BK / "results" / "classes_100000000.npz")["cls_id"]
d9 = np.load(GM / "results" / "classes_1e8_1e9.npz")
ps = np.concatenate([ps8, d9["ps"]])
cls = np.concatenate([cls8, d9["cls_id"]])
del ps8, cls8, d9; gc.collect()
ap = AP_BY_ID[cls]
lp = np.log(ps.astype(np.float64))
w = ap.astype(np.float64) * lp
cum_chi = np.cumsum(w).astype(np.float32)
cum3 = np.cumsum(w**3)          # exact shot-noise (diagonal) 3rd/4th moments
cum4 = np.cumsum(w**4)
del lp, w; gc.collect()
pow_n, pow_wc = [], []
for i in range(int(np.searchsorted(ps, 31623, side="right"))):
    p = int(ps[i]); c = ID2CLS[cls[i]]; lpp = math.log(p)
    pk, k = p * p, 2
    while pk <= 10**9:
        pow_n.append(pk); pow_wc.append(chi8_pow(c, k) * lpp); pk *= p; k += 1
o = np.argsort(pow_n)
pow_n = np.array(pow_n, dtype=np.int64)[o]
pow_cum_c = np.cumsum(np.array(pow_wc)[o])
def psi_chi(x):
    i = np.searchsorted(ps, x, side="right")
    j = np.searchsorted(pow_n, x, side="right")
    return np.where(i > 0, cum_chi[np.maximum(i - 1, 0)], 0.0).astype(np.float64) + \
        np.where(j > 0, pow_cum_c[np.maximum(j - 1, 0)], 0.0)
log(f"engine ready: {len(ps)} primes")

OUT = {"meta": {"date": "2026-09-11",
                "scope": "exploration; not the wall; not RH/GRH; part D = "
                         "consistency-grade only (deep-sweep beta~0 framing)"}}

# ================================================================ A — races @ 1e9
log("A: Chebotarev races at 1e9 ...")
s7 = math.sqrt(7.0)
def merge(v): return np.array([v[0], v[1], v[2], v[3], (v[4] + v[5]) / 2]).real
chi3 = np.array([3, -1, 1, 0, complex(-1, s7) / 2, complex(-1, -s7) / 2])
chi6 = np.array([6, 2, 0, 0, -1, -1], dtype=complex)
chi7v = np.array([7, -1, -1, 1, 0, 0], dtype=complex)
chi8v = np.array([8, 0, 0, -1, 1, 1], dtype=complex)
W = np.stack([merge(chi3), merge(chi6), merge(chi7v), merge(chi8v)], axis=1)
MU = 1 - (1 + merge(chi6) + merge(chi7v) + merge(chi8v))
CLSN = ["1A", "2A", "4A", "3A", "7m"]; SIZE = np.array([1, 21, 42, 56, 48]); G = 168
def loadz(p): return np.array([float(l.split(",")[1]) for l in
                               (CR / "data" / p).read_text().splitlines()[1:]])
ZE = [loadz("zeros_chi3pair_extended.csv"), loadz("zeros_chi6_extended.csv"),
      loadz("zeros_chi7_extended.csv"), loadz("chi8_zeros_153.csv")]
BK2HERE = {0: 0, 1: 1, 2: 3, 3: 2, 4: 4}
xg = np.geomspace(1e3, 1e9, 3000)
idx = np.searchsorted(ps, xg, side="right") - 1
pi_un = (idx + 1) - 2
E = np.zeros((5, len(xg))); lnx = np.log(xg)
for bk, k in BK2HERE.items():
    piC = np.cumsum(cls == bk)[idx]
    E[k] = lnx / np.sqrt(xg) * (G / SIZE[k] * piC - pi_un)
def re_term(z, L):
    c = np.cos(np.outer(L, z)); s = np.sin(np.outer(L, z))
    return (2 * (0.5 * c + z[None, :] * s) / (0.25 + z[None, :]**2)).sum(axis=1)
E_pred = MU[:, None] - sum(np.outer(W[:, j], re_term(ZE[j], lnx)) for j in range(4))
sel = xg >= 1e4
meas_mean = E[:, sel].mean(axis=1)
corr = [float(np.corrcoef(E[k, sel], E_pred[k, sel])[0, 1]) for k in range(5)]
mc_std = json.load(open(CR / "results" / "chebrace_results.json"))["MC_std"]
OUT["A"] = {"classes": CLSN, "predicted_mean": MU.tolist(),
            "measured_logavg_mean_1e4_1e9": meas_mean.tolist(),
            "LI_std": mc_std, "reconstruction_corr": corr}
log("A means: " + "  ".join(f"{c}:{m:+.2f} (pred {p:+.0f}, LI sd {s:.0f})"
    for c, m, p, s in zip(CLSN, meas_mean, MU, mc_std)))
log("A recon corr: " + "  ".join(f"{c}:{r:+.3f}" for c, r in zip(CLSN, corr)))

# ================================================================ B — Gaussianity
log("B: Gaussianity of short-interval increments ...")
def increments(u_lo, u_hi, delta, n=20000):
    u = np.linspace(u_lo, u_hi - delta, n)
    x1, x2 = np.exp(u), np.exp(u + delta)
    D = (psi_chi(x2) - psi_chi(x1)) / np.sqrt(x1)
    return D - D.mean()
NB = 25
def shot_moments(u_lo, u_hi, delta, n=4000):
    """exact diagonal (shot-noise) 3rd/4th cumulant contributions to D."""
    u = np.linspace(u_lo, u_hi - delta, n)
    x1, x2 = np.exp(u), np.exp(u + delta)
    i1 = np.searchsorted(ps, x1, side="right"); i2 = np.searchsorted(ps, x2, side="right")
    d3 = cum3[np.maximum(i2 - 1, 0)] - cum3[np.maximum(i1 - 1, 0)]
    d4 = cum4[np.maximum(i2 - 1, 0)] - cum4[np.maximum(i1 - 1, 0)]
    return float((d3 / x1**1.5).mean()), float((d4 / x1**2).mean())

Brows = []
for wname, ulo, uhi in [("W1", math.log(1e7), math.log(1e8)),
                        ("W2", math.log(1e8), math.log(1e9))]:
    for dl in np.geomspace(2e-4, 0.05, 12):
        D = increments(ulo, uhi, float(dl))
        s = D.std()
        m3, m4 = shot_moments(ulo, uhi, float(dl))
        sk = float((D**3).mean() / s**3)
        ku = float((D**4).mean() / s**4 - 3)
        sk_pred = m3 / s**3                     # diagonal shot-noise skew
        ku_pred = m4 / s**4                     # diagonal shot-noise excess kurtosis
        bl = np.array_split(D, NB)
        sks = [float((b**3).mean() / b.std()**3) for b in bl]
        kus = [float((b**4).mean() / b.std()**4 - 3) for b in bl]
        Brows.append({"window": wname, "delta": float(dl),
                      "skew": sk, "skew_pred_shot": sk_pred,
                      "skew_se": float(np.std(sks) / math.sqrt(NB)),
                      "exkurt": ku, "kurt_pred_shot": ku_pred,
                      "kurt_se": float(np.std(kus) / math.sqrt(NB))})
zs = [abs((r["skew"] - r["skew_pred_shot"]) / r["skew_se"]) for r in Brows]
zk = [abs((r["exkurt"] - r["kurt_pred_shot"]) / r["kurt_se"]) for r in Brows]
zs_naive = [abs(r["skew"] / r["skew_se"]) for r in Brows]
OUT["B"] = {"rows": Brows, "max_abs_z_skew_after_shot": max(zs),
            "max_abs_z_kurt_after_shot": max(zk),
            "max_abs_z_skew_naive": max(zs_naive), "n_tests": len(Brows) * 2}
log(f"B: {len(Brows)} points; naive max|z| skew {max(zs_naive):.2f}; AFTER exact "
    f"shot-noise prediction: skew {max(zs):.2f}, kurt {max(zk):.2f} "
    f"(Bonferroni ~3.1 for 48 tests)")

# ================================================================ C — multiplicative pairs
log("C: multiplicative pairs (p,2p+1), (p,2p-1) ...")
def mult_pair(a_coef, b_coef):
    q = a_coef * ps + b_coef
    j = np.searchsorted(ps, q)
    j = np.minimum(j, len(ps) - 1)
    ok = (ps[j] == q) & (ps > 7) & (q <= ps[-1])
    i1 = np.nonzero(ok)[0]; j1 = j[i1]
    a1, a2 = ap[i1].astype(float), ap[j1].astype(float)
    c1, c2 = cls[i1], cls[j1]
    N = len(a1); prod = a1 * a2
    corr_, se = float(prod.mean()), float(prod.std() / math.sqrt(N))
    null = float(a1.mean() * a2.mean())
    O = np.zeros((5, 5)); np.add.at(O, (c1, c2), 1)
    m1, m2 = O.sum(1) / N, O.sum(0) / N
    Ex = np.outer(m1, m2) * N; mask = Ex > 0
    chi2 = float(((O - Ex)[mask]**2 / Ex[mask]).sum())
    return {"pair": f"(p,{a_coef}p{b_coef:+d})", "N": N, "corr": corr_, "se": se,
            "z": (corr_ - null) / se, "chi2_joint": chi2, "dof": int(mask.sum()) - 9}
sg = mult_pair(2, 1); sm = mult_pair(2, -1)
tt = np.linspace(11, (int(ps[-1]) - 1) // 2, 400000)   # pair needs 2p+1 <= max prime
N_hl = float(2 * C2 * np.trapezoid(1 / (np.log(tt) * np.log(2 * tt)), tt))
for r, nm in [(sg, "Sophie-Germain"), (sm, "2p-1")]:
    log(f"C {nm} {r['pair']}: N={r['N']}  <a a'>={r['corr']:+.5f} +- {r['se']:.5f}  "
        f"z={r['z']:+.2f}  chi2={r['chi2_joint']:.1f}/dof {r['dof']}")
log(f"C count control: SG measured {sg['N']} vs HL-type {N_hl:.0f} "
    f"(ratio {sg['N']/N_hl:.4f})")
OUT["C"] = {"sophie_germain": sg, "twop_minus1": sm,
            "hl_count_pred": N_hl, "sg_ratio": sg["N"] / N_hl}

# ================================================================ D — purity scaling
log("D: spectral-purity scaling across x ...")
gz = np.geomspace(0.05, 5e6, 3000)
def theta_p(T):
    if T < 300:
        return 0.5 * LNQ - 4 * math.log(2 * math.pi) + 4 * float(mp.re(mp.digamma(0.5 + 1j * T)))
    return 0.5 * LNQ - 4 * math.log(2 * math.pi) + 4 * math.log(T)
rho_m = np.array([theta_p(g) / math.pi for g in gz])
def var_gue(delta, u_lo, u_hi):
    g = np.geomspace(0.05, 400 / delta, 6000)
    phi2 = (math.exp(delta) - 2 * math.exp(delta / 2) * np.cos(g * delta) + 1) / (0.25 + g**2)
    rho = np.interp(g, gz, rho_m)
    tot = 0.0
    for u in np.linspace(u_lo, u_hi, 7):
        F = np.clip(u / np.maximum(2 * math.pi * rho, 1e-9), 0, 1)
        tot += float(np.trapezoid(rho * 2 * phi2 * F, g))
    return tot / 7
DU3 = None
U11, g11 = np.load(LANE / "genome_deep_zeros_raw" / "signal_100000000000.npy")
DU3 = float(U11[1] - U11[0])
DELTAS_D = [5 * DU3, 11 * DU3]          # 0.0134, 0.0295 — shared across all windows
Drows = []
for lo10 in np.arange(6.0, 8.5, 0.5):    # exact half-decade windows to 1e9
    ulo, uhi = lo10 * math.log(10), (lo10 + 0.5) * math.log(10)
    for dl in DELTAS_D:
        D = increments(ulo, uhi, dl, n=12000)
        v = float((D**2).mean())
        bl = np.array_split(D, NB)
        se = float(np.array([(b**2).mean() for b in bl]).std() / math.sqrt(NB))
        pg = var_gue(dl, ulo, uhi)
        Drows.append({"x_lo10": float(lo10), "delta": dl, "src": "exact",
                      "ratio": v / pg, "se": se / pg})
for lo10 in (9.0, 9.5, 10.0, 10.5):      # stored-signal windows to 1e11
    ulo, uhi = lo10 * math.log(10), (lo10 + 0.5) * math.log(10)
    i0, i1 = int(np.searchsorted(U11, ulo)), int(np.searchsorted(U11, uhi))
    for m, dl in [(5, 5 * DU3), (11, 11 * DU3)]:
        i1c = min(i1, len(U11) - m)
        Dv = g11[i0 + m:i1c + m] * math.exp(dl / 2) - g11[i0:i1c]
        Dv = Dv - Dv.mean()
        v = float((Dv**2).mean())
        bl = np.array_split(Dv, 15)
        se = float(np.array([(b**2).mean() for b in bl]).std() / math.sqrt(15))
        pg = var_gue(dl, ulo, uhi)
        Drows.append({"x_lo10": float(lo10), "delta": dl, "src": "signal",
                      "ratio": v / pg, "se": se / pg})
# weighted linear fit of ratio vs ln x: slope consistent with 0 <=> no anomalous growth
for dl in DELTAS_D:
    ss = [r for r in Drows if abs(r["delta"] - dl) < 1e-9]
    xs = np.array([(r["x_lo10"] + 0.25) * math.log(10) for r in ss])
    ys = np.array([r["ratio"] for r in ss]); ws = 1 / np.array([r["se"] for r in ss])**2
    A = np.vstack([np.ones_like(xs), xs]).T
    cov = np.linalg.inv(A.T @ (A * ws[:, None]))
    beta = cov @ (A.T @ (ws * ys))
    slope, slope_se = float(beta[1]), float(math.sqrt(cov[1, 1]))
    log(f"D delta={dl:.4f}: ratio slope d(V/V_GUE)/dlnx = {slope:+.4f} +- {slope_se:.4f} "
        f"({'consistent with 0' if abs(slope) < 2*slope_se else 'CHECK'}); "
        f"mean ratio {ys.mean():.3f}")
    OUT.setdefault("D", {})[f"delta_{dl:.4f}"] = {
        "rows": ss, "slope_per_lnx": slope, "slope_se": slope_se,
        "mean_ratio": float(ys.mean())}

(RES / "finale_results.json").write_text(json.dumps(OUT, indent=1))
log("saved results/finale_results.json")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    A0 = ax[0, 0]
    COL = {"1A": "tab:red", "2A": "tab:blue", "4A": "tab:green", "3A": "tab:orange",
           "7m": "tab:purple"}
    for k in range(1, 5):
        A0.semilogx(xg, E[k], lw=.8, color=COL[CLSN[k]], label=CLSN[k])
        A0.semilogx(xg, E_pred[k], lw=.8, ls="--", color=COL[CLSN[k]], alpha=.6)
    A0.set_title("A: races to 1e9 (solid measured, dashed zero-table)", fontsize=9)
    A0.legend(fontsize=7, ncol=4); A0.set_xlabel("x"); A0.set_ylabel("E'_C")
    A1 = ax[0, 1]
    for wname, colr in [("W1", "tab:blue"), ("W2", "tab:red")]:
        s = [r for r in Brows if r["window"] == wname]
        d = [r["delta"] for r in s]
        A1.errorbar(d, [r["skew"] for r in s], yerr=[r["skew_se"] for r in s],
                    fmt="o", ms=3, color=colr, label=f"skew {wname}")
        A1.errorbar(d, [r["exkurt"] for r in s], yerr=[r["kurt_se"] for r in s],
                    fmt="s", ms=3, color=colr, alpha=.5, label=f"ex.kurt {wname}")
    A1.axhline(0, color="k", lw=.7); A1.set_xscale("log")
    A1.set_title("B: increment Gaussianity (skew, excess kurtosis)", fontsize=9)
    A1.legend(fontsize=7); A1.set_xlabel("delta")
    A2 = ax[1, 0]
    labels = ["SG (p,2p+1)", "(p,2p-1)"]
    vals = [sg, sm]
    A2.bar(range(2), [r["corr"] for r in vals],
           yerr=[r["se"] for r in vals], color=["tab:cyan", "tab:gray"], capsize=4)
    A2.axhline(0, color="k", lw=.7); A2.set_xticks(range(2), labels)
    A2.set_title(f"C: multiplicative-pair correlation (N~{sg['N']:.0f}; "
                 f"SG count/HL = {sg['N']/N_hl:.4f})", fontsize=9)
    A3 = ax[1, 1]
    for dl, colr in zip(DELTAS_D, ["tab:blue", "tab:red"]):
        s = [r for r in Drows if abs(r["delta"] - dl) < 1e-9]
        xs = [r["x_lo10"] + 0.25 for r in s]
        A3.errorbar(xs, [r["ratio"] for r in s], yerr=[r["se"] for r in s],
                    fmt="o", ms=4, color=colr, label=f"delta={dl:.3f}")
    A3.axhline(1, color="k", lw=.8, ls="--")
    A3.set_title("D: Var/Var_GUE across x = 1e6 .. 1e11 (flat = on-line)", fontsize=9)
    A3.set_xlabel("log10 x"); A3.set_ylabel("measured / GUE"); A3.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(RES / "finale.png", dpi=130); plt.close(fig)
    log("plot saved")
except Exception as e:  # noqa: BLE001
    log(f"plotting skipped: {e}")
log("DONE")
