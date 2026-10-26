# Fast screening (C5)

Code: `src/screen.[ch]`. CLI: `build/dlscms fastscreen <kd|grid|alt> <same arguments as screen>`. Tests: `tests/test_screen.c` (98 checks). Benchmark: `scripts/bench_c5.sh`, `scripts/plot_bench_c5.py` -> `results/bench_c5.csv`, `bench_c5_summary.csv`, `bench_c5.png`.

## Methods
All three reuse the brute-force per-pair test and refinement (`test_pair`), so they differ only in how candidate pairs are found at each 10 s step.
- **k-d tree** (`SCREEN_KDTREE`): rebuild a tree over the current positions every step; each object queries radius `thr + (vmax_i + vmax_max) * h`.
- **Hashed grid** (`SCREEN_GRID`): same with the uniform grid, cell size = the largest reach.
- **Altitude filter** (`SCREEN_ALT`): sort objects by perigee radius, sweep, skip pairs whose radial ranges `[rmin - pad, rmax + pad]` are further apart than the threshold.

## Correctness
- k-d tree and grid: a pair with an event in a step has `d0 <= thr + (vmax_i + vmax_j) h` (Lipschitz bound, same as the baseline), so the query of the lower-numbered object returns it. Recall is therefore 100% by construction, and checked: the event files of kd and grid are **byte-identical** to the brute-force file on every dataset in `results/bench_c5.csv` (synthetic 1k-20k, and the real catalogue with SGP4, 12,504 events). Unit tests also require identical candidate counts.
- Altitude filter: **not a proven bound** under SGP4 (the radius leaves the mean-element range [perigee, apogee]). Measured on the real catalogue over 10 minutes: excursion up to 26.7 km. The CLI pad is 50 km. Identical output on all datasets tested (with pad 30 km on those runs), but recall is only as good as the pad. See L23.

## Results (CPU seconds, one thread, 1 h window, 10 s step, 10 km; `results/bench_c5_summary.csv`)
| Objects | brute force | k-d tree | grid | altitude filter |
|---|---|---|---|---|
| 1,000 | 0.29 | 0.16 | 0.19 | 0.34 |
| 5,000 | 6.6 | 1.2 | 1.2 | 8.1 |
| 10,000 | 25.6 | 3.3 | 3.1 | 33.0 |
| 20,000 | 98.2 | 9.8 | 9.2 | 141.5 |
| Real catalogue + SGP4, 18,191 | 102.3 | 15.6 | 14.8 | (not run for 1 h: > 170 s) |

- Pairs examined fall by ~15,000x (synthetic) / ~11,500x (real), but the time falls "only" 10x at 20k objects: after the index, the floor is propagation (about 6.5 million state evaluations per hour of window) and refining the ~600k-900k intervals that survive the Lipschitz bound. With a tiny threshold the real run still takes 12.2 s.
- Grid and k-d tree are within ~7% of each other everywhere; neither is consistently faster.
- The altitude filter is slower than brute force in LEO: only 2.5x fewer pairs, because LEO radial ranges overlap heavily, and the sweep adds overhead. It is a negative result worth reporting; it would pay off for mixed LEO/MEO/GEO catalogues, which are out of scope.
- Timing uses `clock()` on a VM with 2 cores (AMD Ryzen 9 8940HX host); single run per cell; repeats and an environment capture are planned for C7 (L11).

## Known limits / next
- Cost is now dominated by propagation and refinement; ideas for C7: a cheaper pre-refinement test, caching SGP4 states between steps, a larger step with a tighter bound.
- Index reach uses the Lipschitz speed bound with a 1% margin (`conj_vmax`); under SGP4 the true speed can differ slightly from the Kepler vmax (velocity of SGP4 vs d(position)/dt up to 1.5 m/s). The identical-output check on the real catalogue is the empirical support.
