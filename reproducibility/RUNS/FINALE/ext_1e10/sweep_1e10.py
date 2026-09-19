# sweep_1e10.py — exact-engine extension (1e9,1e10]: a_p via the validated overflow-safe
# numba kernel (gate-A cross-certified vs flint today), binned psi_chi + Lambda^2 on a
# fine u-grid. Checkpointed + resumable. FINALE part-D follow-up.
# 2026-09-11, Stenberg + Claude (Fable 5). Exploration; not the wall; not RH/GRH.
import json, math, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
HERE.mkdir(exist_ok=True)
LANE = Path(r"C:\Users\selin\OneDrive\Desktop\Ilya Riemanns")
BK = HERE.parent.parent / "BK_OFFDIAG_2026-09-11"
GM = HERE.parent.parent / "GAPMAP_2026-09-11"
sys.path.insert(0, str(LANE / "genome_deep_zeros_raw"))
sys.path.insert(0, str(BK / "data"))
from _genome_fast_big import genome_batch          # noqa: E402  overflow-safe numba
from gen_an_chi8_REFERENCE import prime_sieve      # noqa: E402

LO, HI, BLK = 10**9, 10**10, 5 * 10**7
U0, U1, DU = math.log(LO), math.log(HI), 2e-4
NBINS = int((U1 - U0) / DU) + 2
CKPT = HERE / "checkpoint_sweep.npz"
FINAL = HERE / "bins_1e9_1e10.npz"
t00 = time.time()
def log(m): print(f"[{time.time()-t00:7.1f}s] {m}", flush=True)

if FINAL.exists():
    log("final bins already exist — nothing to do"); sys.exit(0)

base = prime_sieve(int(HI**0.5) + 1).astype(np.int64)
log(f"base primes to sqrt(1e10): {len(base)}")

S1 = np.zeros(NBINS); S2 = np.zeros(NBINS)
counts = {8: 0, 0: 0, -1: 0, 1: 0}
b0 = 0
if CKPT.exists():
    d = np.load(CKPT, allow_pickle=True)
    S1, S2, b0 = d["S1"], d["S2"], int(d["next_block"])
    counts = {int(k): int(v) for k, v in zip(d["ck"], d["cv"])}
    log(f"resuming at block {b0}")

blocks = list(range(LO, HI, BLK))
for bi in range(b0, len(blocks)):
    lo = blocks[bi]; hi = min(lo + BLK, HI)
    comp = np.zeros(hi - lo, dtype=bool)
    for q in base:
        q = int(q)
        if q * q >= hi:
            break
        start = max(q * q, ((lo + q - 1) // q) * q)
        if start < hi:
            comp[start - lo::q] = True
    pr = lo + np.nonzero(~comp)[0]
    if bi == b0:
        log(f"block {bi+1}: {len(pr)} primes; first genome_batch call — numba "
            "compile happens here (can take minutes) ...")
    ap = genome_batch(pr.astype(np.int64))
    if bi == b0:
        log(f"block {bi+1}: kernel compiled + ran")
    for v in (8, 0, -1, 1):
        counts[v] += int((ap == v).sum())
    lpr = np.log(pr.astype(np.float64))
    w = ap * lpr
    idx = ((lpr - U0) / DU).astype(np.int64)
    S1 += np.bincount(idx, weights=w, minlength=NBINS)[:NBINS]
    S2 += np.bincount(idx, weights=w * w, minlength=NBINS)[:NBINS]
    if (bi + 1) % 5 == 0 or bi == len(blocks) - 1:
        np.savez(CKPT, S1=S1, S2=S2, next_block=bi + 1,
                 ck=list(counts.keys()), cv=list(counts.values()))
        tot = sum(counts.values())
        log(f"block {bi+1}/{len(blocks)} (to {hi:.2e}); {tot} primes; "
            f"fracs 8:{counts[8]/tot:.5f} 0:{counts[0]/tot:.5f} "
            f"-1:{counts[-1]/tot:.5f} 1:{counts[1]/tot:.5f}")

# prime powers p^k in (1e9,1e10]: exact classes from the <=1e9 caches (p <= 1e5)
ID2CLS = ["1A", "2A", "3A", "4A", "7", "p3", "p7"]
def chi8_pow(c, k):
    if c == "1A": return 8
    if c == "2A": return 8 if k % 2 == 0 else 0
    if c == "3A": return 8 if k % 3 == 0 else -1
    if c == "4A": return 8 if k % 4 == 0 else 0
    if c == "7":  return 8 if k % 7 == 0 else 1
    return 1 if c == "p3" else 0
ps_small = prime_sieve(10**5).astype(np.int64)
cls_small = np.load(BK / "results" / "classes_100000000.npz")["cls_id"][:len(ps_small)]
npow = 0
for p, cid in zip(ps_small, cls_small):
    p = int(p); c = ID2CLS[cid]; lpp = math.log(p)
    pk, k = p * p, 2
    while pk <= HI:
        if pk > LO:
            u = k * lpp
            i = int((u - U0) / DU)
            wv = chi8_pow(c, k) * lpp
            S1[i] += wv; S2[i] += wv * wv; npow += 1
        pk *= p; k += 1
log(f"added {npow} prime-power terms (exact classes)")

tot = sum(counts.values())
np.savez_compressed(FINAL, S1=S1, S2=S2, U0=U0, DU=DU,
                    ck=list(counts.keys()), cv=list(counts.values()))
log(f"DONE: {tot} primes in (1e9,1e10]; Chebotarev "
    f"8:{counts[8]/tot:.5f} 0:{counts[0]/tot:.5f} -1:{counts[-1]/tot:.5f} "
    f"1:{counts[1]/tot:.5f} (expect .00595/.375/.33333/.28571)")
if CKPT.exists():
    CKPT.unlink()
