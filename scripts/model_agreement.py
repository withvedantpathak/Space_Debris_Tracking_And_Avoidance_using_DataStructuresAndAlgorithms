"""C3/C4 evidence: do close pairs found with our propagation model agree with close pairs found with SGP4?
Real catalogue, near-earth objects only (period < 225 min), common time = median epoch + t. At each t we list all
pairs closer than R km under each model (scipy cKDTree) and compare the two sets.
Run after `scripts/build.sh`; needs python packages sgp4, scipy, pandas, matplotlib.
Writes results/model_agreement.csv and results/model_agreement.png."""
import io, os, subprocess, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from sgp4 import omm
from sgp4.api import Satrec
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
merged = os.path.join(root, "build", "merged_real.csv")
cat = pd.read_csv(merged)
near = (cat.MEAN_MOTION > 1440 / 225).values
with open(merged) as fh: rows = list(omm.parse_csv(fh))
sats = []
for r in rows:
    s = Satrec(); omm.initialize(s, r); sats.append(s)
ref_jd = float(np.median([s.jdsatepoch + s.jdsatepochF for s in sats]))
def sgp4_pos(t):
    jd = ref_jd + t / 86400.0
    out = np.full((len(sats), 3), np.nan)
    for i, s in enumerate(sats):
        if not near[i]: continue
        e, p, v = s.sgp4(np.floor(jd), jd - np.floor(jd))
        if e == 0: out[i] = p
    return out
def ours_pos(t, flags=3):
    txt = subprocess.run([os.path.join(root, "build", "dlscms"), "prop", merged, str(t), str(flags), "common"], capture_output=True, text=True, check=True).stdout
    return pd.read_csv(io.StringIO(txt))[["x_km", "y_km", "z_km"]].values
def pairs(P, R):
    ok = ~np.isnan(P).any(axis=1)
    idx = np.where(ok)[0]
    tr = cKDTree(P[idx])
    return {(idx[a], idx[b]) for a, b in tr.query_pairs(R)}

out, scat = [], {}
for t in [0, 3600, 6 * 3600, 24 * 3600, 72 * 3600]:
    A, B = sgp4_pos(t), ours_pos(t)
    for R in [5, 10, 20, 50]:
        ps, po = pairs(A, R), pairs(B, R)
        both = len(ps & po)
        out.append(dict(t_h=t / 3600, R_km=R, sgp4_pairs=len(ps), our_pairs=len(po), both=both,
                        recall=both / len(ps) if ps else np.nan, precision=both / len(po) if po else np.nan))
    ps = list(pairs(A, 50))
    if ps:
        ia = np.array([p[0] for p in ps]); ib = np.array([p[1] for p in ps])
        scat[t] = (np.linalg.norm(A[ia] - A[ib], axis=1), np.linalg.norm(B[ia] - B[ib], axis=1))
res = pd.DataFrame(out)
res.round(4).to_csv(os.path.join(root, "results", "model_agreement.csv"), index=False)
print(res.round(3).to_string(index=False))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
for t, c in zip([0, 24 * 3600, 72 * 3600], ["#1b7f5c", "#c2410c", "#444444"]):
    if t in scat: ax[0].scatter(*scat[t], s=4, alpha=.5, color=c, label=f"t = {t//3600} h")
ax[0].plot([0, 50], [0, 50], "k--", lw=1); ax[0].set(xlabel="pair distance, SGP4 (km)", ylabel="same pair, our model (km)", title="Distance of SGP4-close pairs (< 50 km)")
ax[0].legend(); ax[0].grid(alpha=.3)
for R, c in zip([5, 10, 20, 50], ["#c2410c", "#1b7f5c", "#2563eb", "#444444"]):
    d = res[res.R_km == R]; ax[1].plot(d.t_h, d.recall * 100, "o-", color=c, label=f"R = {R} km")
ax[1].set(xscale="symlog", xlabel="hours after reference epoch", ylabel="% of SGP4 close pairs also close in our model", title="Pair-set recall vs SGP4", ylim=(0, 105))
ax[1].legend(); ax[1].grid(alpha=.3)
fig.suptitle("Do both models agree on which objects are close? (real catalogue, near-earth objects)")
fig.tight_layout(); fig.savefig(os.path.join(root, "results", "model_agreement.png"), dpi=130)
