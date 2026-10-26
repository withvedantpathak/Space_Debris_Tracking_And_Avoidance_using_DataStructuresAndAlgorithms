"""L34: calibrate the epoch-age growth g (km per day) of the age-aware screening distance from element-set HISTORY.
Needs a CSV in the same OMM format as data/raw/*.csv but with SEVERAL rows per NORAD_CAT_ID (different EPOCH), e.g. Space-Track
'gp_history' or a CelesTrak archive. NOT part of the delivered data: the repository only holds one element set per object.
Method (TLE-to-TLE): for every pair (older set i, newer set j) of the same object with 0 < age = epoch_j - epoch_i <= MAXAGE days,
propagate set i to the epoch of set j with SGP4 and measure the distance to set j's own position at its epoch. The newer set is
treated as the reference (its own error, ~1 km at epoch, is therefore included: this OVER-states the growth slightly).
Output: results/c8_tle_error_growth.csv with the error percentiles per age bin, and the suggested g = the slope through the origin
of the chosen percentile (default 95) of error versus age, printed as km/day.
Usage: python3 scripts/calibrate_age_from_history.py <history.csv> [max_age_days=5] [percentile=95]"""
import sys, os, numpy as np, pandas as pd
from sgp4 import omm
from sgp4.api import Satrec
path = sys.argv[1]; MAXAGE = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0; PCT = float(sys.argv[3]) if len(sys.argv) > 3 else 95.0
with open(path) as fh: rows = list(omm.parse_csv(fh))
by = {}
for r in rows:
    s = Satrec(); omm.initialize(s, r); by.setdefault(int(r["NORAD_CAT_ID"]), []).append(s)
age, err = [], []
for sats in by.values():
    sats.sort(key=lambda s: s.jdsatepoch + s.jdsatepochF)
    for a in range(len(sats)):
        for b in range(a + 1, len(sats)):
            sa, sb = sats[a], sats[b]; days = (sb.jdsatepoch - sa.jdsatepoch) + (sb.jdsatepochF - sa.jdsatepochF)
            if days <= 0 or days > MAXAGE: continue
            e1, r1, _ = sa.sgp4(sb.jdsatepoch, sb.jdsatepochF); e2, r2, _ = sb.sgp4(sb.jdsatepoch, sb.jdsatepochF)
            if e1 or e2: continue
            age.append(days); err.append(float(np.linalg.norm(np.array(r1) - np.array(r2))))
age, err = np.array(age), np.array(err)
if len(age) < 50: sys.exit(f"only {len(age)} element-set pairs: not enough history to calibrate")
edges = np.arange(0, MAXAGE + 0.5, 0.5); out = []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (age > lo) & (age <= hi)
    if m.sum() >= 10: out.append(dict(age_from_days=lo, age_to_days=hi, pairs=int(m.sum()), err_median_km=np.median(err[m]), err_p95_km=np.percentile(err[m], PCT), err_p75_km=np.percentile(err[m], 75), err_p90_km=np.percentile(err[m], 90), err_p99_km=np.percentile(err[m], 99)))
tab = pd.DataFrame(out); tab.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "c8_tle_error_growth.csv"), index=False)
mid = (tab.age_from_days + tab.age_to_days) / 2; g = float((mid * tab.err_p95_km).sum() / (mid * mid).sum())
print(tab.to_string(index=False)); print(f"\nSuggested g (slope through origin of the p{PCT:g} error): {g:.2f} km/day from {len(age)} element-set pairs")
for col, name in (('err_median_km', 'median'), ('err_p75_km', 'p75'), ('err_p90_km', 'p90'), ('err_p95_km', f'p{PCT:g}'), ('err_p99_km', 'p99')):
    print(f"slope through origin, {name}: {float((mid * tab[col]).sum() / (mid * mid).sum()):.2f} km/day")
