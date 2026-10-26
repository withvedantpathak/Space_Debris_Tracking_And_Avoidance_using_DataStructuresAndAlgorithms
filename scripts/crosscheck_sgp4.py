"""C3 evidence: compare our propagator (two-body and two-body+J2) with SGP4 on the real catalogue.
For every object, both models start from the object's own epoch and are evaluated dt later.
SGP4 is the python `sgp4` package fed with the same OMM rows. Error = |r_ours - r_sgp4| in km.
Writes results/propagation_error.csv and results/propagation_error.png.  Run from the repo root."""
import glob, io, os, subprocess, sys
import numpy as np, pandas as pd
from sgp4 import omm
from sgp4.api import Satrec
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
exe = os.path.join(root, "build", "dlscms")
files = sorted(glob.glob(os.path.join(root, "data", "raw", "*.csv")))
MU, RE = 398600.4418, 6378.137

# merged catalogue file (unique ids) so the C tool and SGP4 see identical objects
frames = [pd.read_csv(f) for f in files]
cat = pd.concat(frames).drop_duplicates("NORAD_CAT_ID").reset_index(drop=True)
merged = os.path.join(root, "build", "merged_real.csv")
cat.to_csv(merged, index=False)

n = cat.MEAN_MOTION * 2 * np.pi / 86400
a = (MU / n**2) ** (1 / 3)
apo = a * (1 + cat.ECCENTRICITY) - RE
cls = np.where(apo < 2000, "LEO", np.where(cat.ECCENTRICITY >= 0.25, "HEO",
      np.where((a > 41664) & (a < 42664) & (cat.ECCENTRICITY < 0.1), "GEO", "MEO/other")))
cat["cls"] = cls

sats = []
with open(merged) as fh:
    for row in omm.parse_csv(fh):
        s = Satrec(); omm.initialize(s, row); sats.append(s)
assert len(sats) == len(cat)

def sgp4_at(dt):
    out = np.full((len(sats), 3), np.nan)
    for i, s in enumerate(sats):
        e, r, v = s.sgp4(s.jdsatepoch, s.jdsatepochF + dt / 86400.0)
        if e == 0: out[i] = r
    return out

def ours(dt, j2):
    txt = subprocess.run([exe, "prop", merged, str(dt), str(j2)], capture_output=True, text=True, check=True).stdout
    d = pd.read_csv(io.StringIO(txt))
    assert (d.norad_id.values == cat.NORAD_CAT_ID.values).all()
    return d[["x_km", "y_km", "z_km"]].values

rows = []
for dt_h in [0, 1, 3, 6, 12, 24, 48, 72, 168, 336]:
    dt = dt_h * 3600.0
    ref = sgp4_at(dt)
    ok = ~np.isnan(ref).any(axis=1)
    for j2, model in [(0, "two-body"), (1, "two-body + J2"), (3, "two-body + J2 + drag")]:
        err = np.linalg.norm(ours(dt, j2) - ref, axis=1)
        for c in ["LEO", "MEO/other", "GEO", "HEO", "all"]:
            m = ok & ((cat.cls == c) if c != "all" else True)
            e = err[m]
            if len(e) == 0: continue
            rows.append(dict(dt_h=dt_h, model=model, orbit_class=c, objects=int(m.sum()), sgp4_failed=int((~ok & ((cat.cls == c) if c != "all" else True)).sum()),
                             median_km=np.median(e), p95_km=np.percentile(e, 95), max_km=e.max()))
res = pd.DataFrame(rows)
res.round(4).to_csv(os.path.join(root, "results", "propagation_error.csv"), index=False)

fig, ax = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
colm = {"two-body": "#888888", "two-body + J2": "#c2410c", "two-body + J2 + drag": "#1b7f5c"}
for k, c in enumerate(["LEO", "GEO", "HEO"]):
    for model, color in colm.items():
        d = res[(res.orbit_class == c) & (res.model == model)]
        ax[k].plot(d.dt_h, d.median_km, "o-", color=color, label=f"{model} median")
        ax[k].plot(d.dt_h, d.p95_km, "s--", color=color, alpha=.7, label=f"{model} 95th pct")
    ax[k].set_xscale("symlog", linthresh=1)
    ax[k].set(yscale="log", xlabel="hours after element epoch", title=f"{c} (n={int(d.objects.iloc[0])})")
    ax[k].grid(alpha=.3)
ax[0].set_ylabel("position error vs SGP4 (km)"); ax[0].legend(fontsize=7)
fig.suptitle("Propagation error against SGP4, real catalogue 2026-10-20 (short-period terms, lunisolar not modelled)")
fig.tight_layout(); fig.savefig(os.path.join(root, "results", "propagation_error.png"), dpi=130)
print(res[res.orbit_class.isin(["LEO", "GEO", "HEO", "all"])].round(2).to_string(index=False))
