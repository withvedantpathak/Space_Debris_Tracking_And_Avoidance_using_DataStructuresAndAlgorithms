"""C4 evidence: check that the exact brute-force baseline recovers every planted pair of a synthetic catalogue.
For each planted pair (data/synthetic/synth_N_truth.csv) the encounter recurs at tca + k*period/2; every recurrence
inside the window must appear in results/baseline_synth_N.csv with the planted miss distance (hits, miss < threshold)
and must be absent when the planted miss is above the threshold (decoys).
Usage: python3 scripts/verify_baseline.py <N> <window_s> <threshold_km>  -> results/baseline_planted_N.csv"""
import sys, os, numpy as np, pandas as pd
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
N, window, thr = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
truth = pd.read_csv(f"{root}/data/synthetic/synth_{N}_truth.csv")
ev = pd.read_csv(f"{root}/results/baseline_synth_{N}.csv")
ev["lo"] = ev[["norad_id_a", "norad_id_b"]].min(axis=1); ev["hi"] = ev[["norad_id_a", "norad_id_b"]].max(axis=1)
rows = []
for _, p in truth.iterrows():
    lo, hi = min(p.id_a, p.id_b), max(p.id_a, p.id_b)
    e = ev[(ev.lo == lo) & (ev.hi == hi)]
    half = p.period_s / 2
    ks = [k for k in range(-400, 400) if 1.0 < p.tca_s_after_epoch + k * half <= window - 1.0]
    exp_t = np.array([p.tca_s_after_epoch + k * half for k in ks])
    is_hit = p.miss_km < thr
    found, terr, merr = 0, 0.0, 0.0
    for t in exp_t:
        d = np.abs(e.tca_s.values - t)
        if len(d) and d.min() <= 0.05:
            found += 1; r = e.iloc[int(d.argmin())]
            terr = max(terr, float(d.min())); merr = max(merr, abs(r.miss_km - p.miss_km))
    extra = len(e) - found
    ok = (found == len(exp_t) and extra == 0 and merr < 1e-3) if is_hit else (len(e) == 0)
    rows.append(dict(id_a=p.id_a, id_b=p.id_b, planted_miss_km=p.miss_km, kind="hit" if is_hit else "decoy",
                     expected_events=len(exp_t) if is_hit else 0, found_events=found, extra_events=extra,
                     max_tca_error_s=round(terr, 5), max_miss_error_km=round(merr, 7), ok=ok))
res = pd.DataFrame(rows)
res.to_csv(f"{root}/results/baseline_planted_{N}.csv", index=False)
print(f"synth_{N}: {len(res)} planted pairs, {int((res.kind=='hit').sum())} hits, {int((res.kind=='decoy').sum())} decoys; all ok = {bool(res.ok.all())}")
print(res.groupby("kind").agg(pairs=("ok", "size"), ok=("ok", "sum"), events=("found_events", "sum"), worst_tca_s=("max_tca_error_s", "max"), worst_miss_km=("max_miss_error_km", "max")).to_string())
sys.exit(0 if res.ok.all() else 1)
