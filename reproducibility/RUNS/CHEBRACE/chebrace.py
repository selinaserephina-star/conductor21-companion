# chebrace.py — Chebotarev races for the Trinks field (PSL(2,7), x^7-7x+3),
# bias predicted from our own zero tables, measured exactly to 1e8.
# 2026-09-11, Stenberg + Claude (Fable 5). Exploration; not the wall; not RH/GRH.
import json, math, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
RES = HERE / "results"
RES.mkdir(exist_ok=True)
sys.path.insert(0, str(HERE / "data"))
from gen_an_chi8_REFERENCE import prime_sieve  # noqa: E402

RNG = np.random.default_rng(168)
t00 = time.time()
def log(m): print(f"[{time.time()-t00:7.1f}s] {m}", flush=True)

# ---------------------------------------------------------------- character data
# classes (observable): 1A, 2A, 4A, 3A, 7m(=7A u 7B); sizes:
CLS = ["1A", "2A", "4A", "3A", "7m"]
SIZE = np.array([1, 21, 42, 56, 48])
G = 168
# full table on 6 true classes (1A,2A,4A,3A,7A,7B), sizes (1,21,42,56,24,24):
s7 = math.sqrt(7.0)
chi3 = np.array([3, -1, 1, 0, complex(-1, s7) / 2, complex(-1, -s7) / 2])
chi3b = chi3.conj()
chi6 = np.array([6, 2, 0, 0, -1, -1], dtype=complex)
chi7 = np.array([7, -1, -1, 1, 0, 0], dtype=complex)
chi8 = np.array([8, 0, 0, -1, 1, 1], dtype=complex)
chi1 = np.ones(6, dtype=complex)
sizes6 = np.array([1, 21, 42, 56, 24, 24])
for name, (a, b) in {"<1,1>": (chi1, chi1), "<3,3>": (chi3, chi3), "<3,3b>": (chi3, chi3b),
                     "<6,6>": (chi6, chi6), "<7,7>": (chi7, chi7), "<8,8>": (chi8, chi8),
                     "<6,7>": (chi6, chi7), "<7,8>": (chi7, chi8)}.items():
    ip = complex((sizes6 * a * b.conj()).sum()) / G
    expected = 1.0 if name in ("<1,1>", "<3,3>", "<6,6>", "<7,7>", "<8,8>") else 0.0
    assert abs(ip - expected) < 1e-12, (name, ip)
log("character table orthogonality: OK")

# w_chi(C) = class-averaged chi value (real after 7A/7B merge); order:
# columns = [chi3pair(one L6 table, coefficient = chi3 avg), chi6, chi7, chi8]
def merge(v):
    return np.array([v[0], v[1], v[2], v[3], (v[4] + v[5]) / 2]).real
W = np.stack([merge(chi3), merge(chi6), merge(chi7), merge(chi8)], axis=1)  # (5 cls,4 chars)
r2 = 1 + merge(chi6) + merge(chi7) + merge(chi8)          # squares count per class rep
MU = 1 - r2                                                # predicted E' means
assert np.allclose(r2, [22, 2, 0, 1, 1])
assert abs((SIZE * MU).sum()) < 1e-12                      # weighted means cancel
log(f"predicted bias means mu_C = {dict(zip(CLS, MU.tolist()))}")

# ---------------------------------------------------------------- zero tables
def load_csv_col1(p):
    return np.array([float(l.split(",")[1]) for l in
                     (HERE / "data" / p).read_text().splitlines()[1:]])
ZEROS = {
    "chi3pair": load_csv_col1("zeros_chi3pair_extended.csv"),
    "chi6": load_csv_col1("zeros_chi6_extended.csv"),
    "chi7": load_csv_col1("zeros_chi7_extended.csv"),
    "chi8": load_csv_col1("chi8_zeros_153.csv"),
}
CHARS = ["chi3pair", "chi6", "chi7", "chi8"]
for c in CHARS:
    z = ZEROS[c]
    log(f"  {c}: {len(z)} zeros, gamma in [{z[0]:.3f}, {z[-1]:.3f}]")

def S_in(z):                       # in-table variance factor sum 2/(1/4+g^2)
    return float((2 / (0.25 + z**2)).sum())
def V_tail(z, k_end=5):            # flat-density tail estimate beyond table end
    T = z[-1]
    rho = k_end / (z[-1] - z[-1 - k_end])
    return float(2 * rho * 2 * (math.pi / 2 - math.atan(2 * T)))
SIN = {c: S_in(ZEROS[c]) for c in CHARS}
VTL = {c: V_tail(ZEROS[c]) for c in CHARS}
log(f"variance factors in-table: { {c: round(v,3) for c,v in SIN.items()} }")
log(f"tail estimates:            { {c: round(v,3) for c,v in VTL.items()} }")

# ---------------------------------------------------------------- measured races
log("loading cached Frobenius classes (<=1e8) + sieve ...")
X = 10**8
ps = prime_sieve(X)
dat = np.load(HERE.parent / "BK_OFFDIAG_2026-09-11" / "results" / f"classes_{X}.npz")
cls_id = dat["cls_id"]             # 0..4 = 1A,2A,3A,4A,7 ; 5,6 = ramified 3,7
assert len(cls_id) == len(ps)
# map BK ids (1A,2A,3A,4A,7) -> our order (1A,2A,4A,3A,7m)
BK2HERE = {0: 0, 1: 1, 2: 3, 3: 2, 4: 4}
xg = np.geomspace(1e3, 1e8, 2500)
idx = np.searchsorted(ps, xg, side="right") - 1
pi_un = (idx + 1) - 2                                     # minus ramified 3,7 (x>=1e3)
cum = {k: np.cumsum(cls_id == bk) for bk, k in BK2HERE.items()}
E = np.zeros((5, len(xg)))
lnx = np.log(xg)
for k in range(5):
    piC = cum[k][idx]
    E[k] = lnx / np.sqrt(xg) * (G / SIZE[k] * piC - pi_un)
sel = xg >= 1e4                                            # stats window
meas_mean = E[:, sel].mean(axis=1)
log("measured log-avg means [1e4,1e8]: "
    + "  ".join(f"{c}:{m:+.2f}" for c, m in zip(CLS, meas_mean)))
log("predicted means:                  "
    + "  ".join(f"{c}:{m:+.0f}" for c, m in zip(CLS, MU)))

# ---------------------------------------------------------------- reconstruction
def re_term(z, L):                 # sum over table zeros of 2 Re[x^{i g}/(1/2+i g)]
    # Re = [0.5 cos(gL) + g sin(gL)] / (1/4+g^2)
    c = np.cos(np.outer(L, z)); s = np.sin(np.outer(L, z))
    return (2 * (0.5 * c + z[None, :] * s) / (0.25 + z[None, :]**2)).sum(axis=1)
osc = np.stack([re_term(ZEROS[c], lnx) for c in CHARS])    # (4 chars, x)
E_pred = MU[:, None] - W @ osc                             # (5, x)
corr = [float(np.corrcoef(E[k, sel], E_pred[k, sel])[0, 1]) for k in range(5)]
rmsd = [float((E[k, sel] - E_pred[k, sel]).std()) for k in range(5)]
log("reconstruction corr(measured, zero-table): "
    + "  ".join(f"{c}:{r:+.3f}" for c, r in zip(CLS, corr)))
log("reconstruction residual rms:               "
    + "  ".join(f"{c}:{r:.2f}" for c, r in zip(CLS, rmsd)))

# ---------------------------------------------------------------- LI Monte Carlo
log("LI Monte Carlo (shared phases across classes, tails as per-char Gaussians) ...")
NMC = 200_000
allz, allw = [], []                                        # per zero: gamma + w row
for j, c in enumerate(CHARS):
    for g in ZEROS[c]:
        allz.append(g); allw.append(W[:, j])
allz = np.array(allz); allw = np.array(allw)               # (Nz,), (Nz,5)
amp = 2 / np.sqrt(0.25 + allz**2)                          # oscillation amplitude
Esim = np.zeros((NMC, 5))
CH = 20_000
for i0 in range(0, NMC, CH):
    n = min(CH, NMC - i0)
    ph = RNG.uniform(0, 2 * np.pi, size=(n, len(allz)))
    Esim[i0:i0 + n] = MU[None, :] - (np.cos(ph) * amp[None, :]) @ allw
for j, c in enumerate(CHARS):                              # tail: correlated Gaussians
    gt = RNG.normal(0, math.sqrt(VTL[c]), NMC)
    Esim -= np.outer(gt, W[:, j])
P_neg = (Esim < 0).mean(axis=0)
log("LI P(E'_C < 0): " + "  ".join(f"{c}:{p:.3f}" for c, p in zip(CLS, P_neg)))
pair_delta = {}
for a in range(5):
    for b in range(a + 1, 5):
        d = float((Esim[:, a] > Esim[:, b]).mean())
        pair_delta[f"{CLS[a]}>{CLS[b]}"] = d
log("LI pairwise lead-probabilities:")
for k, v in pair_delta.items():
    log(f"   P({k}) = {v:.4f}")
# measured lead fractions in log measure, [1e4, 1e8]
meas_lead = {}
for a in range(5):
    for b in range(a + 1, 5):
        meas_lead[f"{CLS[a]}>{CLS[b]}"] = float((E[a, sel] > E[b, sel]).mean())
log("measured lead fractions (log-time in [1e4,1e8], ONE noisy sample of the race):")
for k in pair_delta:
    log(f"   {k}: measured {meas_lead[k]:.3f} vs LI {pair_delta[k]:.4f}")

# ---------------------------------------------------------------- outputs
out = {
    "meta": {"date": "2026-09-11", "X": X, "grid": [1e3, 1e8, len(xg)],
             "zeros": {c: len(ZEROS[c]) for c in CHARS},
             "scope": "exploration; prediction side assumes GRH(family)+LI; "
                      "measurement unconditional; not the wall; not RH/GRH"},
    "classes": CLS, "sizes": SIZE.tolist(), "r2": r2.tolist(),
    "predicted_mean": MU.tolist(),
    "measured_logavg_mean_1e4_1e8": meas_mean.tolist(),
    "reconstruction_corr": corr, "reconstruction_rms": rmsd,
    "variance_factors": {"in_table": SIN, "tail_est": VTL},
    "LI_P_negative": P_neg.tolist(),
    "LI_pair_lead": pair_delta, "measured_pair_lead_logfrac": meas_lead,
    "MC_mean": Esim.mean(axis=0).tolist(), "MC_std": Esim.std(axis=0).tolist(),
}
(RES / "chebrace_results.json").write_text(json.dumps(out, indent=1))
np.savez(RES / "race_curves.npz", xg=xg, E=E, E_pred=E_pred)
log("saved results/chebrace_results.json + race_curves.npz")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    COLORS = {"1A": "tab:red", "2A": "tab:blue", "4A": "tab:green",
              "3A": "tab:orange", "7m": "tab:purple"}
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True,
                                 gridspec_kw={"height_ratios": [1, 1.4]})
    a1.semilogx(xg, E[0], color=COLORS["1A"], lw=1, label="1A measured")
    a1.semilogx(xg, E_pred[0], color=COLORS["1A"], lw=1, ls="--", alpha=.7,
                label="1A zero-table prediction")
    a1.axhline(MU[0], color="k", lw=.7, ls=":")
    a1.set_ylabel("E'_1A"); a1.legend(fontsize=8)
    a1.set_title("Chebotarev races, Trinks field x^7-7x+3 (PSL(2,7)) — measured to 1e8 "
                 "vs prediction from our zero tables")
    for k in range(1, 5):
        a2.semilogx(xg, E[k], color=COLORS[CLS[k]], lw=1, label=f"{CLS[k]} measured")
        a2.semilogx(xg, E_pred[k], color=COLORS[CLS[k]], lw=1, ls="--", alpha=.6)
        a2.axhline(MU[k], color=COLORS[CLS[k]], lw=.7, ls=":")
    a2.set_ylabel("E'_C (dashed: prediction; dotted: bias mean)")
    a2.set_xlabel("x"); a2.legend(fontsize=8, ncol=4)
    fig.tight_layout(); fig.savefig(RES / "cheb_races.png", dpi=130); plt.close(fig)

    fig, A = plt.subplots(figsize=(7, 4.5))
    xpos = np.arange(5)
    A.bar(xpos - 0.18, MU, width=0.36, label="predicted mean 1 - r2(C)", color="tab:gray")
    A.bar(xpos + 0.18, meas_mean, width=0.36, label="measured log-avg [1e4,1e8]",
          color="tab:cyan")
    A.errorbar(xpos + 0.18, meas_mean, yerr=Esim.std(axis=0), fmt="none",
               ecolor="k", capsize=3, label="LI predicted spread (1 sigma)")
    A.set_xticks(xpos, CLS); A.axhline(0, color="k", lw=.6)
    A.set_ylabel("E' bias"); A.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(RES / "cheb_bias_bars.png", dpi=130); plt.close(fig)
    log("plots saved")
except Exception as e:  # noqa: BLE001
    log(f"plotting skipped: {e}")
log("DONE")
