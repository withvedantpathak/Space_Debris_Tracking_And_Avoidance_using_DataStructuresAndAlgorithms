"""L33: why does linear interpolation of positions (Gomez-style) give extra events, and how accurate is it?
Uses the event files of scripts/run_c7.sh (interp step): reference results/c7_events_step5.csv (bisection, 5 s step) and
results/c7_events_lin<h>.csv (linear, step h) and results/c8_events_hyb<h>.csv (hybrid: linear for relative speed >= 1 km/s, bisection below; made by scripts/run_c8_hybrid.sh). Matching rule = scripts/compare_events.py (same pair, TCA within 2 s, each reference
event used once). Metrics are reported for ALL events and for CROSSING events only (relative speed >= 1 km/s at TCA, the same
definition as the *_topk_crossing.csv list of `analyze`). The 'extra' events are also split by relative speed.
Writes results/c8_l33_interp.csv.  Usage: python3 scripts/l33_analysis.py"""
import os, numpy as np, pandas as pd
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
TOL = 2.0
def norm(d):
    d = d.copy(); d["a"] = d[["norad_id_a", "norad_id_b"]].min(axis=1); d["b"] = d[["norad_id_a", "norad_id_b"]].max(axis=1); return d
def match(ref, tst):
    g = {k: v for k, v in tst.groupby(["a", "b"])}; used = set(); matched_ref = []
    for r in ref.itertuples():
        d = g.get((r.a, r.b)); best = None
        if d is not None:
            for j, t in zip(d.index, d.itertuples()):
                if j in used: continue
                dt = abs(t.tca_s - r.tca_s)
                if dt <= TOL and (best is None or dt < best[1]): best = (j, dt)
        if best: used.add(best[0])
        matched_ref.append(best is not None)
    return np.array(matched_ref, bool), used
ref_all = norm(pd.read_csv(os.path.join(root, "results", "c7_events_step5.csv")))
rows = []
for h, method in [(h, m) for h in (5, 10, 20, 30, 60) for m in ("linear", "hybrid")]:
    fn = f"c7_events_lin{h}.csv" if method == "linear" else f"c8_events_hyb{h}.csv"
    if not os.path.exists(os.path.join(root, "results", fn)): continue
    tst_all = norm(pd.read_csv(os.path.join(root, "results", fn)))
    for scope, vmin in (("all", 0.0), ("crossing", 1.0)):
        ref = ref_all[ref_all.rel_speed_kms >= vmin].reset_index(drop=True); tst = tst_all[tst_all.rel_speed_kms >= vmin].reset_index(drop=True)
        mref, used = match(ref, tst); tp = int(mref.sum()); nref, ntst = len(ref), len(tst)
        rec = tp / max(nref, 1); pre = tp / max(ntst, 1); f1 = 2 * pre * rec / (pre + rec) if pre + rec else 0.0
        extra = tst.loc[[i for i in tst.index if i not in used]]
        rows.append(dict(step_s=h, method=method, scope=scope, ref_events=nref, lin_events=ntst, matched=tp, missed=nref - tp, extra=ntst - tp,
                         recall=round(rec, 4), precision=round(pre, 4), f1=round(f1, 4),
                         extra_with_rel_speed_below_0p5=int((extra.rel_speed_kms < 0.5).sum())))
out = pd.DataFrame(rows); out.to_csv(os.path.join(root, "results", "c8_l33_interp.csv"), index=False)
print(out.to_string(index=False))
