"""C2 evidence: compare the real CelesTrak snapshot with the synthetic catalogue.
Reads data/raw/*.csv and data/synthetic/synth_20000*.csv; writes results/dataset_stats.csv and
results/dataset_compare.png. Altitude uses the same two-body formula as src/catalog.c."""
import glob, os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MU, RE = 398600.4418, 6378.137
root = os.path.join(os.path.dirname(__file__), "..")
def load(files):
    d = pd.concat([pd.read_csv(f) for f in files]).drop_duplicates("NORAD_CAT_ID")
    n = d.MEAN_MOTION * 2 * np.pi / 86400
    d["a"] = (MU / n**2) ** (1 / 3); d["alt"] = d.a - RE
    d["apo"] = d.a * (1 + d.ECCENTRICITY) - RE
    return d
real = load(sorted(glob.glob(f"{root}/data/raw/*.csv")))
syn = load([f"{root}/data/synthetic/synth_20000.csv"])
sets = {"real (CelesTrak 2026-10-20)": real, "synthetic (20,000, seed 1)": syn}

rows = []
for name, d in sets.items():
    leo = d[d.apo < 2000]
    geo = d[(d.a > 41664) & (d.a < 42664) & (d.ECCENTRICITY < 0.1)]
    rows.append(dict(dataset=name, objects=len(d), leo_share=len(leo) / len(d), geo_share=len(geo) / len(d),
                     high_ecc_share=(d.ECCENTRICITY >= 0.25).mean(), median_leo_alt_km=leo.alt.median(),
                     median_leo_incl_deg=leo.INCLINATION.median(), median_ecc=d.ECCENTRICITY.median()))
pd.DataFrame(rows).round(5).to_csv(f"{root}/results/dataset_stats.csv", index=False)

col = ["#444444", "#c2410c"]
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
for (name, d), c in zip(sets.items(), col):
    leo = d[d.apo < 2000]
    ax[0].hist(leo.alt, bins=np.arange(0, 2100, 100), weights=np.full(len(leo), 100 / len(d)), histtype="step", lw=2, color=c, label=name)
    ax[1].hist(d.INCLINATION, bins=np.arange(0, 190, 10), weights=np.full(len(d), 100 / len(d)), histtype="step", lw=2, color=c, label=name)
ax[0].set(xlabel="altitude of low-orbit objects (km)", ylabel="% of all objects per 100 km", title="Altitude (LEO)")
ax[1].set(xlabel="inclination (deg)", ylabel="% of all objects per 10 deg", title="Inclination")
shares = pd.DataFrame(rows).set_index("dataset")[["leo_share", "geo_share", "high_ecc_share"]] * 100
x = np.arange(3); w = 0.38
for k, (name, c) in enumerate(zip(shares.index, col)):
    ax[2].bar(x + (k - .5) * w, shares.loc[name], w, color=c, label=name)
ax[2].set_xticks(x); ax[2].set_xticklabels(["LEO", "GEO", "e >= 0.25"]); ax[2].set(ylabel="% of objects", title="Orbit class share")
for a in ax: a.grid(alpha=.3)
ax[0].legend(fontsize=8)
fig.suptitle("Synthetic filler population vs the real catalogue it is modelled on")
fig.tight_layout(); fig.savefig(f"{root}/results/dataset_compare.png", dpi=130)
print(pd.DataFrame(rows).round(4).to_string(index=False))
