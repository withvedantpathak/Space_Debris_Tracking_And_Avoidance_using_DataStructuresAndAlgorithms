# C7 Experiments (all numbers trace to `results/`, produced by `scripts/run_c7.sh`)

Reproduce: `bash scripts/run_c7.sh env|scaling "<N list>"|threshold|step|interp|cell|age|summary`, then `python3 scripts/plot_c7.py`. Environment: `results/environment.txt`. Timings are CPU seconds on one VM (L32). Real data: CelesTrak active + three debris files, SGP4, LEO scope, 1 h window from the median epoch unless stated.

## Results
| Experiment | Finding | File |
|---|---|---|
| Scaling (synthetic, 1 h, 10 km, 3 seeds, repeats) | All 88 runs identical to brute force. At N=10,000: brute 26.6 s, k-d tree 3.1 s, grid 2.9 s | `c7_scaling.csv` |
| Scaling, largest point (synthetic N=50,000, 1 h, 10 km, seed 1, one run) | k-d tree 45.4 s, grid 41.3 s, same 56,780 events; brute force not run (L51) | `c7_scaling_50k.csv` |
| Threshold 1/2/5/10 km | Events 48/174/2,254/12,504; largest cluster 5/5/24/7,197; exact cover 39/154/1,346/3,687 (matches MILP in every case) | `c7_threshold.csv` |
| Step 5-60 s (bisection) | Events identical (2,254) for every step; candidates grow 375,749 -> 7,078,327 and CPU 13 s -> 89 s | `c7_step.csv` |
| Linear interpolation (E9) | 2.5x to 18x less CPU at 10-60 s steps, recall 0.991-0.998, precision ~0.92: all extra events are slow co-moving pairs (L33). The hybrid (`kd:hybrid`) removes them: precision 1.0, recall 1.0 at 5-10 s, 0.998 at 60 s | `c7_interp.csv`, `c8_l33_interp.csv` |
| Grid cell scale (10 min, 5 km) | Same events and pair checks for every scale; best near 1x (2.2 s); 0.25x is 4x slower; memory ~22.5 MB for all | `c7_cell.csv` |
| Age-aware distance (heuristic) | g = 0/1/2/3 km/day: events 2,254/3,046/3,900/4,686; largest cluster 24/169/1,354/2,502 | `c7_age.csv` |

## Reading
- Indexing is exact and scales much better than brute force; the gain is limited by propagation and refinement.
- Step size does not change bisection results here; smaller steps are cheaper, so use the smallest step tested (5 s) for this window.
- A 5 km threshold keeps clusters small enough for exact mitigation; 10 km produces a giant cluster.
- Age-aware inflation can merge clusters quickly and is not validated (L34).

Figure: `results/c7_summary.png`. Table: `results/summary.csv`.
