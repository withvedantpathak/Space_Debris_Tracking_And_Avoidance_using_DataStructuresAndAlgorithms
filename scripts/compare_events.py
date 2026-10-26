"""Compare two event files (screen/fastscreen output). Usage: compare_events.py <reference.csv> <test.csv> [tca_tol_s=2]
Prints one CSV line: ref_events,test_events,matched,missed,extra,recall,precision,miss_abs_diff_median_km,miss_abs_diff_p95_km,tca_abs_diff_median_s
An event matches if the pair is the same and the TCAs differ by at most tca_tol_s (each reference event is used once)."""
import sys, numpy as np, pandas as pd
ref = pd.read_csv(sys.argv[1]); tst = pd.read_csv(sys.argv[2]); tol = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
def norm(d):
    d = d.copy(); d["a"] = d[["norad_id_a", "norad_id_b"]].min(axis=1); d["b"] = d[["norad_id_a", "norad_id_b"]].max(axis=1); return d
ref, tst = norm(ref), norm(tst)
g = {k: v for k, v in tst.groupby(["a", "b"])}
used = set(); dm, dt, matched = [], [], 0
for r in ref.itertuples():
    d = g.get((r.a, r.b))
    if d is None: continue
    best = None
    for j, t in zip(d.index, d.itertuples()):
        if j in used: continue
        if abs(t.tca_s - r.tca_s) <= tol and (best is None or abs(t.tca_s - r.tca_s) < best[1]): best = (j, abs(t.tca_s - r.tca_s), abs(t.miss_km - r.miss_km))
    if best: used.add(best[0]); matched += 1; dt.append(best[1]); dm.append(best[2])
n_ref, n_tst = len(ref), len(tst)
f = lambda x: f"{x:.4f}"
print(",".join(map(str, [n_ref, n_tst, matched, n_ref - matched, n_tst - matched, f(matched / max(n_ref, 1)), f(matched / max(n_tst, 1)),
      f(np.median(dm)) if dm else "nan", f(np.percentile(dm, 95)) if dm else "nan", f(np.median(dt)) if dt else "nan"])))
