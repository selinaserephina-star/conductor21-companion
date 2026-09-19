# retest_D.py — FINALE part-D tail retest with the exact (1e9,1e10] sweep:
# Var/Var_GUE for windows [1e9,10^9.5], [10^9.5,1e10]; refit purity slope with exact
# points 1e6..1e10, stored signal only for [1e10,1e11].
# 2026-09-11, Stenberg + Claude (Fable 5). Consistency-grade only; not the wall.
import json, math, sys, time
from pathlib import Path

import numpy as np
import mpmath as mp

HERE = Path(__file__).parent
LANE = Path(r"C:\Users\selin\OneDrive\Desktop\Ilya Riemanns")
mp.mp.dps = 25
LNQ = 10 * math.log(21.0)
t00 = time.time()
def log(m): print(f"[{time.time()-t00:7.1f}s] {m}", flush=True)

d = np.load(HERE / "bins_1e9_1e10.npz")
S1, U0, DU = d["S1"], float(d["U0"]), float(d["DU"])
psi = np.cumsum(S1)                      # psi_chi(u) - psi_chi(ln 1e9), on the bin grid
U = U0 + DU * np.arange(len(psi))
log(f"bins loaded: {len(psi)} (du={DU:g})")

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

NB = 25
DELTAS = [0.0134, 0.0295]
rows = []
for lo10 in (9.0, 9.5):
    ulo, uhi = lo10 * math.log(10), (lo10 + 0.5) * math.log(10)
    for dl in DELTAS:
        m = int(round(dl / DU)); dle = m * DU
        i0 = max(int((ulo - U0) / DU), 0)
        i1 = min(int((uhi - U0) / DU), len(psi) - m - 1)
        seg = slice(i0, i1)
        D = (psi[i0 + m:i1 + m] - psi[seg]) / np.exp(U[seg] / 2)
        D = D - D.mean()
        v = float((D**2).mean())
        bl = np.array_split(D, NB)
        se = float(np.array([(b**2).mean() for b in bl]).std() / math.sqrt(NB))
        pg = var_gue(dle, max(ulo, U0), uhi)
        rows.append({"x_lo10": lo10, "delta": dle, "src": "exact-1e10",
                     "ratio": v / pg, "se": se / pg})
        log(f"[{lo10},{lo10+0.5}] delta={dle:.4f}: Var/GUE = {v/pg:.3f} +- {se/pg:.3f}")

# combined refit: exact points from finale JSON (1e6..1e9) + these + signal 1e10..1e11
fin = json.load(open(HERE.parent / "results" / "finale_results.json"))
out_fits = {}
for dl_key, dl in [("delta_0.0134", 0.0134), ("delta_0.0295", 0.0295)]:
    old = fin["D"][dl_key]["rows"]
    pts = [r for r in old if r["src"] == "exact"] + \
          [r for r in rows if abs(r["delta"] - dl) < 2e-3] + \
          [r for r in old if r["src"] == "signal" and r["x_lo10"] >= 10.0]
    xs = np.array([(r["x_lo10"] + 0.25) * math.log(10) for r in pts])
    ys = np.array([r["ratio"] for r in pts])
    ws = 1 / np.array([r["se"] for r in pts])**2
    A = np.vstack([np.ones_like(xs), xs]).T
    cov = np.linalg.inv(A.T @ (A * ws[:, None]))
    beta = cov @ (A.T @ (ws * ys))
    slope, sse = float(beta[1]), float(math.sqrt(cov[1, 1]))
    # exact-only fit (1e6..1e10), no signal points at all
    pe = [r for r in pts if "exact" in r["src"]]
    xe = np.array([(r["x_lo10"] + 0.25) * math.log(10) for r in pe])
    ye = np.array([r["ratio"] for r in pe]); we = 1 / np.array([r["se"] for r in pe])**2
    Ae = np.vstack([np.ones_like(xe), xe]).T
    cove = np.linalg.inv(Ae.T @ (Ae * we[:, None]))
    be = cove @ (Ae.T @ (we * ye))
    slope_e, sse_e = float(be[1]), float(math.sqrt(cove[1, 1]))
    out_fits[dl_key] = {"rows_new": [r for r in rows if abs(r["delta"] - dl) < 2e-3],
                        "slope_all": slope, "slope_all_se": sse,
                        "slope_exact_only_1e6_1e10": slope_e, "slope_exact_se": sse_e}
    log(f"{dl_key}: slope(all incl signal tail) {slope:+.4f}+-{sse:.4f}; "
        f"EXACT-ONLY 1e6..1e10: {slope_e:+.4f}+-{sse_e:.4f} "
        f"({'consistent with 0' if abs(slope_e) < 2*sse_e else 'CHECK'})")

(HERE / "retest_D_results.json").write_text(json.dumps(
    {"rows_new": rows, "fits": out_fits,
     "scope": "consistency-grade; deep-sweep beta~0 framing; not the wall"}, indent=1))
log("saved retest_D_results.json — DONE")
