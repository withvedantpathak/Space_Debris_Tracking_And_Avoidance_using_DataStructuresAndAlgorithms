#!/usr/bin/env bash
# L34 (C8/C9): evaluate the epoch-age-aware screening distance for a finer set of growth rates g (km per day of element age, per object)
# on the same real-data setup as scripts/run_c7.sh (1 h, 5 km, SGP4 + LEO scope, fresh elements <= 3 days from the median epoch).
# Also records the age distribution of the objects in scope, so the effect of g can be read off.
# Writes results/c8_age_calibration.csv and results/c8_age_distribution.csv.   g near 1 km/day is the order of magnitude implied by the abstract of Levit & Marshall (arXiv:1002.2277): their method reaches about 100 m/day, a tenfold improvement over TLEs propagated with SGP4. Not calibrated here.
set -e
cd "$(dirname "$0")/.."
bash scripts/build.sh > /dev/null
R="data/raw/active.csv data/raw/cosmos-2251-debris.csv data/raw/fengyun-1c-debris.csv data/raw/iridium-33-debris.csv"
REAL="12 3 6.4"
CSV=results/c8_age_calibration.csv
echo "growth_km_per_day,events,pairs,objects,components,largest_cluster,greedy_cover,exact_cover" > $CSV
for g in 0 0.5 1 1.5 2 3; do
    DLSCMS_AGE_GROWTH=$g ./build/dlscms fastscreen kd /tmp/c8_age_ev.csv 3600 10 5 $REAL $R > /tmp/c8_a.log
    ./build/dlscms analyze /tmp/c8_age_ev.csv data/raw/active.csv 10 /tmp/c8_age_an 300000 > /tmp/c8_an.log
    c=$(grep '^C6SUMMARY' /tmp/c8_an.log | cut -d, -f3-)
    IFS=, read ev pairs nodes mov comps single ln lp agree topk gc x lb unres nt dn worse tg tx <<< "$c"
    echo "$g,$ev,$pairs,$nodes,$comps,$ln,$gc,$x" >> $CSV; echo "g=$g done"
done
python3 - <<'PY'
import pandas as pd, numpy as np, glob
cat = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("data/raw/*.csv"))]).drop_duplicates("NORAD_CAT_ID")
jd = 2440587.5 + (pd.to_datetime(cat.EPOCH) - pd.Timestamp("1970-01-01")).dt.total_seconds().values / 86400.0
ref = np.sort(jd)[len(jd) // 2]                      # same rule as the C code: sorted()[n/2]
keep = (np.abs(jd - ref) <= 3.0) & (cat.MEAN_MOTION.values >= 6.4)
age = np.abs(jd[keep] - ref)
rows = [("objects_in_scope_before_leo_filter", int(keep.sum()))] + [(f"age_percentile_{q}_days", round(float(np.percentile(age, q)), 3)) for q in (50, 75, 90, 95, 99)] + [("age_max_days", round(float(age.max()), 3))]
pd.DataFrame(rows, columns=["statistic", "value"]).to_csv("results/c8_age_distribution.csv", index=False)
print(pd.DataFrame(rows, columns=["statistic", "value"]).to_string(index=False))
PY
cat $CSV
