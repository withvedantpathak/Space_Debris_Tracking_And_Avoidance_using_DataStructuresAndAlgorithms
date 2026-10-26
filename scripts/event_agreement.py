"""C4 evidence: do conjunction EVENTS found with our model agree with events found with SGP4?
Input: results/baseline_real_1h_25km.csv from
  build/dlscms screen results/baseline_real_1h_25km.csv 3600 10 25 3 3 6.4 data/raw/*.csv
(reference epoch = median epoch, objects with epoch within 3 days and period < 225 min, J2 + drag, 25 km threshold).
SGP4 side: positions and velocities every 10 s (SatrecArray), candidate pairs by cKDTree within a proven reach
(thr + L*step), closest approach by ternary search on the cubic Hermite curve through (r0,v0,r1,v1), interior minima
only. Events with miss < 25 km are collected for both models and compared.
Writes results/event_agreement.csv (summary), results/event_agreement_pairs.csv (matched events), results/event_agreement.png."""
import os, sys, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from sgp4 import omm
from sgp4.api import Satrec, SatrecArray
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
# usage: event_agreement.py [threshold_km=25] [our_events_csv=results/baseline_real_1h_25km.csv] [tag=""]
THR = float(sys.argv[1]) if len(sys.argv) > 1 else 25.0
OURS = sys.argv[2] if len(sys.argv) > 2 else os.path.join(root, "results", "baseline_real_1h_25km.csv")
TAG = sys.argv[3] if len(sys.argv) > 3 else ""
WINDOW, STEP, FRESH, MM_MIN = 3600.0, 10.0, 3.0, 6.4
RE_KM = 6378.137       # object scope as in src/catalog.h: LEO 175-2000 km (only applied when TAG says so)
SCOPE = TAG.startswith("scope")
MU = 398600.4418
cat = pd.read_csv(os.path.join(root, "build", "merged_real.csv"))
jd = 2440587.5 + (pd.to_datetime(cat.EPOCH) - pd.Timestamp("1970-01-01")).dt.total_seconds().values / 86400.0
ref = np.sort(jd)[len(jd) // 2]                      # same rule as the C code: sorted()[n/2]
keep = (np.abs(jd - ref) <= FRESH) & (cat.MEAN_MOTION.values >= MM_MIN)
if SCOPE:
    a_ = (MU / (cat.MEAN_MOTION.values * 2 * np.pi / 86400) ** 2) ** (1 / 3)
    keep &= (a_ * (1 - cat.ECCENTRICITY.values) - RE_KM >= 175.0) & (a_ * (1 + cat.ECCENTRICITY.values) - RE_KM <= 2000.0)
with open(os.path.join(root, "build", "merged_real.csv")) as fh: rows = list(omm.parse_csv(fh))
sats, ids, vmax = [], [], []
for i, r in enumerate(rows):
    if not keep[i]: continue
    s = Satrec(); omm.initialize(s, r); sats.append(s); ids.append(int(cat.NORAD_CAT_ID.iloc[i]))
    n = cat.MEAN_MOTION.iloc[i] * 2 * np.pi / 86400; a = (MU / n**2) ** (1 / 3); e = cat.ECCENTRICITY.iloc[i]
    vmax.append(1.01 * np.sqrt(MU / a) * np.sqrt((1 + e) / (1 - e)))
ids = np.array(ids); vmax = np.array(vmax)
N = len(sats); print(f"{N} objects, reference JD {ref:.6f}")
times = np.arange(0.0, WINDOW + STEP, STEP); K = len(times)
arr = SatrecArray(sats)
P = np.empty((K, N, 3)); V = np.empty((K, N, 3)); bad = np.zeros(N, bool)
for c0 in range(0, K, 40):
    tt = times[c0:c0 + 40]; j = ref + tt / 86400.0
    e, r, v = arr.sgp4(np.floor(j), j - np.floor(j))
    P[c0:c0 + 40] = np.transpose(r, (1, 0, 2)); V[c0:c0 + 40] = np.transpose(v, (1, 0, 2)); bad |= (e != 0).any(axis=1)
P[:, bad] = 1e9 + np.arange(bad.sum())[None, :, None] * 1e4     # failed objects: park far away
L = 2 * vmax.max()

def hermite_min(p0, v0, p1, v1, h):
    """vectorised ternary search for the minimum distance of the pair curve on s in [0,1]"""
    lo = np.zeros(len(p0)); hi = np.ones(len(p0))
    def dist2(s):
        s = s[:, None]; s2 = s * s; s3 = s2 * s
        h00 = 2 * s3 - 3 * s2 + 1; h10 = s3 - 2 * s2 + s; h01 = -2 * s3 + 3 * s2; h11 = s3 - s2
        d = h00 * p0 + h10 * h * v0 + h01 * p1 + h11 * h * v1
        return (d * d).sum(axis=1)
    for _ in range(45):
        m1 = lo + (hi - lo) / 3; m2 = hi - (hi - lo) / 3
        f1, f2 = dist2(m1), dist2(m2)
        lo = np.where(f1 < f2, lo, m1); hi = np.where(f1 < f2, m2, hi)
    s = 0.5 * (lo + hi)
    return s, np.sqrt(dist2(s))

events = []
for k in range(K - 1):
    tr = cKDTree(P[k]); pr = tr.query_pairs(THR + L * STEP, output_type="ndarray")
    if len(pr) == 0: continue
    a, b = pr[:, 0], pr[:, 1]
    d0 = np.linalg.norm(P[k][a] - P[k][b], axis=1); d1 = np.linalg.norm(P[k + 1][a] - P[k + 1][b], axis=1)
    Lij = vmax[a] + vmax[b]
    ok = d0 + d1 - Lij * STEP <= 2 * THR
    a, b = a[ok], b[ok]
    if len(a) == 0: continue
    s, dm = hermite_min(P[k][a] - P[k][b], V[k][a] - V[k][b], P[k + 1][a] - P[k + 1][b], V[k + 1][a] - V[k + 1][b], STEP)
    sel = (s > 1e-4) & (s < 1 - 1e-4) & (dm < THR)
    for ia, ib, ss, dd in zip(a[sel], b[sel], s[sel], dm[sel]):
        vr = np.linalg.norm((1 - ss) * (V[k][ia] - V[k][ib]) + ss * (V[k + 1][ia] - V[k + 1][ib]))
        events.append((min(ids[ia], ids[ib]), max(ids[ia], ids[ib]), times[k] + ss * STEP, dd, vr))
sg = pd.DataFrame(events, columns=["a", "b", "tca", "miss", "vrel"])
sg = sg[(sg.tca > 0) & (sg.tca <= WINDOW)].reset_index(drop=True)
ours = pd.read_csv(OURS)
ours["a"] = ours[["norad_id_a", "norad_id_b"]].min(axis=1); ours["b"] = ours[["norad_id_a", "norad_id_b"]].max(axis=1)
ours = ours.rename(columns={"tca_s": "tca", "miss_km": "miss", "rel_speed_kms": "vrel"})[["a", "b", "tca", "miss", "vrel"]]
print(f"SGP4 events (<{THR:.0f} km): {len(sg)}, our events: {len(ours)}")

TOL = 15.0
def match(src, dst):
    """for each row of src, the nearest-in-time event of dst with the same pair within TOL seconds (or NaN)"""
    g = {k: v for k, v in dst.groupby(["a", "b"])}
    out = []
    for r in src.itertuples():
        d = g.get((r.a, r.b)); best = None
        if d is not None:
            dt = (d.tca.values - r.tca); i = int(np.abs(dt).argmin())
            if abs(dt[i]) <= TOL: best = (dt[i], d.miss.values[i])
        out.append(best if best else (np.nan, np.nan))
    return np.array(out)
m_sg = match(sg, ours); sg["dt"], sg["our_miss"] = m_sg[:, 0], m_sg[:, 1]
m_ou = match(ours, sg); ours["dt"], ours["sg_miss"] = m_ou[:, 0], m_ou[:, 1]
rows = []
for cls, mask_s, mask_o in [("all", np.ones(len(sg), bool), np.ones(len(ours), bool)),
                            ("crossing (v_rel > 1 km/s)", sg.vrel.values > 1, ours.vrel.values > 1),
                            ("slow / co-orbiting (v_rel <= 1 km/s)", sg.vrel.values <= 1, ours.vrel.values <= 1)]:
    for thr in sorted({t for t in [5, 10, THR] if t <= THR}):
        S = sg[mask_s & (sg.miss < thr).values]; O = ours[mask_o & (ours.miss < thr).values]
        for margin in sorted({thr, 10, 15, 20, 25, THR}):
            if margin < thr: continue
            rec = (S.our_miss < margin).mean() if len(S) else np.nan
            rows.append(dict(event_class=cls, sgp4_threshold_km=thr, our_threshold_km=margin, sgp4_events=len(S),
                             sgp4_events_found_by_ours=int((S.our_miss < margin).sum()), recall=rec))
        prec = (O.sg_miss < thr).mean() if len(O) else np.nan
        rows.append(dict(event_class=cls + " [precision]", sgp4_threshold_km=thr, our_threshold_km=thr, sgp4_events=len(O),
                         sgp4_events_found_by_ours=int((O.sg_miss < thr).sum()), recall=prec))
res = pd.DataFrame(rows); res.round(4).to_csv(os.path.join(root, "results", f"event_agreement{TAG}.csv"), index=False)
mm = sg.dropna(subset=["our_miss"])
mm.round(4).to_csv(os.path.join(root, "results", f"event_agreement_pairs{TAG}.csv"), index=False)
print(res.round(3).to_string(index=False))
cr = mm[mm.vrel > 1]
print(f"matched crossing events: {len(cr)}; |miss_ours - miss_sgp4|: median {np.median(np.abs(cr.our_miss - cr.miss)):.2f} km, p95 {np.percentile(np.abs(cr.our_miss - cr.miss), 95):.2f} km; |tca diff| median {np.median(np.abs(cr.dt)):.2f} s")
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
ax[0].scatter(cr.miss, cr.our_miss, s=6, alpha=.5, color="#c2410c"); ax[0].plot([0, 25], [0, 25], "k--", lw=1)
ax[0].set(xlabel="miss distance, SGP4 (km)", ylabel="same event, our model (km)", title=f"Crossing events (v_rel > 1 km/s), n={len(cr)}"); ax[0].grid(alpha=.3)
ax[1].hist(cr.dt, bins=40, color="#1b7f5c"); ax[1].set(xlabel="TCA difference, ours - SGP4 (s)", ylabel="events", title="Timing of closest approach"); ax[1].grid(alpha=.3)
s10 = sg[(sg.miss < 10) & (sg.vrel > 1)]
margins = [10, 15, 20, 25]
ax[2].bar([str(m) for m in margins], [(s10.our_miss < m).mean() * 100 for m in margins], color="#444444")
ax[2].set(xlabel="our screening threshold (km)", ylabel="% of SGP4 events (< 10 km) found", title=f"Recall of SGP4 crossing events < 10 km (n={len(s10)})", ylim=(0, 105)); ax[2].grid(alpha=.3)
fig.suptitle("Event-level agreement with SGP4, real catalogue, 1 h window, near-earth objects with fresh elements")
fig.tight_layout(); fig.savefig(os.path.join(root, "results", f"event_agreement{TAG}.png"), dpi=130)
