#!/usr/bin/env python3
"""Build results/summary.csv (metric, value, source file) from the experiment CSVs (C4-C7)."""
import csv, os, statistics as st
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
def rd(p): return list(csv.DictReader(open(p)))
out = []
def add(m, v, s): out.append((m, v, s))

s = 'results/c7_scaling.csv'
if os.path.exists(s):
    rows = rd(s)
    add('scaling: rows identical to brute force', f"{sum(r['identical_to_bruteforce']=='yes' for r in rows)}/{len(rows)}", s)
    for n in sorted({int(r['n']) for r in rows}):
        for m in ('brute', 'kd', 'grid'):
            t = [float(r['cpu_seconds']) for r in rows if int(r['n']) == n and r['method'] == m]
            if t: add(f'scaling N={n} {m} median cpu s', f"{st.median(t):.3f}", s)
s = 'results/c7_threshold.csv'
if os.path.exists(s):
    for r in rd(s):
        add(f"threshold {r['threshold_km']} km: events/largest cluster/exact cover", f"{r['events']}/{r['largest_cluster']}/{r['exact_cover']}", s)
s = 'results/c7_step.csv'
if os.path.exists(s):
    for r in rd(s):
        add(f"step {r['step_s']} s: cpu s/candidates/recall", f"{r['cpu_seconds']}/{r['candidates']}/{r['recall_vs_5s']}", s)
s = 'results/c7_interp.csv'
if os.path.exists(s):
    for r in rd(s):
        if r['method'] == 'linear':
            add(f"linear interp step {r['step_s']} s: cpu s/recall/precision", f"{r['cpu_seconds']}/{r['recall_vs_ref']}/{r['precision_vs_ref']}", s)
s = 'results/c7_cell.csv'
if os.path.exists(s):
    for r in rd(s): add(f"cell {r['method']} scale {r['cell_scale']}: cpu s/peak RSS kB", f"{r['cpu_seconds']}/{r['peak_rss_kb']}", s)
s = 'results/c7_age.csv'
if os.path.exists(s):
    for r in rd(s): add(f"age growth {r['growth_km_per_day']} km/day: events/largest cluster/exact cover", f"{r['events']}/{r['largest_cluster']}/{r['exact_cover']}", s)
s = 'results/c6_summary.csv'
if os.path.exists(s):
    for r in rd(s):
        if r['variant'] == 'reductions':
            add(f"C6 {os.path.basename(r['events_file'])}: events/greedy/exact cover", f"{r['events']}/{r['greedy_cover']}/{r['exact_cover']}", s)
with open('results/summary.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['metric', 'value', 'source']); w.writerows(out)
print(f'wrote results/summary.csv ({len(out)} rows)')
