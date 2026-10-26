#!/usr/bin/env python3
"""C7 figures from results/c7_*.csv -> results/c7_summary.png"""
import csv, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
rd = lambda p: list(csv.DictReader(open(p)))
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
sc = rd('results/c7_scaling.csv')
for m in ('brute', 'kd', 'grid'):
    ns = sorted({int(r['n']) for r in sc if r['method'] == m})
    ys = [sorted(float(r['cpu_seconds']) for r in sc if int(r['n']) == n and r['method'] == m) for n in ns]
    ax[0][0].plot(ns, [y[len(y)//2] for y in ys], 'o-', label=m)
ax[0][0].set_xscale('log'); ax[0][0].set_yscale('log'); ax[0][0].set_xlabel('objects N'); ax[0][0].set_ylabel('CPU s (median)'); ax[0][0].set_title('Scaling (synthetic, 1 h, 10 km)'); ax[0][0].legend()
th = rd('results/c7_threshold.csv'); x = [float(r['threshold_km']) for r in th]
ax[0][1].plot(x, [int(r['largest_cluster']) for r in th], 'o-', label='largest cluster (objects)')
ax[0][1].plot(x, [int(r['exact_cover']) for r in th], 's-', label='exact cover size')
ax[0][1].set_yscale('log'); ax[0][1].set_xlabel('threshold km'); ax[0][1].set_title('Threshold sensitivity (real, 1 h)'); ax[0][1].legend()
st_ = rd('results/c7_step.csv'); hs = [int(r['step_s']) for r in st_]
ax[1][0].plot(hs, [float(r['cpu_seconds']) for r in st_], 'o-', label='bisection')
for r in rd('results/c7_interp.csv'):
    pass
li = [r for r in rd('results/c7_interp.csv') if r['method'] == 'linear']
ax[1][0].plot([int(r['step_s']) for r in li], [float(r['cpu_seconds']) for r in li], 's-', label='linear interp')
ax[1][0].set_xlabel('step s'); ax[1][0].set_ylabel('CPU s'); ax[1][0].set_title('Step size: cost'); ax[1][0].legend()
ag = rd('results/c7_age.csv'); g = [float(r['growth_km_per_day']) for r in ag]
ax[1][1].plot(g, [int(r['events']) for r in ag], 'o-', label='events')
ax[1][1].plot(g, [int(r['largest_cluster']) for r in ag], 's-', label='largest cluster')
ax[1][1].set_xlabel('age growth km/day (heuristic)'); ax[1][1].set_title('Age-aware distance'); ax[1][1].legend()
plt.tight_layout(); plt.savefig('results/c7_summary.png', dpi=130); print('wrote results/c7_summary.png')
