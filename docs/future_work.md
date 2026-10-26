# Out of scope, limitations and future work (single list)

Maintained as the project goes on: **add every new limitation or out-of-scope decision here** (and the detailed evidence row in `docs/known_limitations.md`). Status: OUT = deliberately excluded from this project, OPEN = a weakness we may still fix, DONE = fixed (kept for history), DOC = documented, accepted.

## A. Out of scope (deliberate)
| # | Item | Why excluded | Idea for future work |
|---|---|---|---|
| O1 | Deep-space / HEO objects (GEO, MEO, Molniya; about 800 of 19,440) | Needs SDP4 (lunisolar terms, resonance); our model error is 25 km median, up to 9,000 km | Port SDP4 from Vallado 2006, validate against the Python reference like `src/sgp4.c`; widen scope past 2000 km |
| O2 | Collision probability and covariance | TLE/OMM carry no covariance; spheres with a single screening distance | Use operator CDM covariance where available; estimate along-track growth from epoch age; Pc computation (Foster/Alfano) |
| O3 | Operational use, live data feeds, GUI | University project, snapshot data | Scheduled CelesTrak download, alerting, web UI |
| O4 | Force models beyond SGP4 (higher-order gravity, solar pressure, third bodies, maneuvers) | SGP4 is the catalogue standard; synthetic truth is two-body | Numerical propagator (Gomez et al. use C++/AutoPas with perturbations); compare with SGP4 |
| O5 | Object size / radius / mass; passive vs active pair rules from Gomez et al. | No radius data in OMM files | Join SATCAT RCS data; apply radius-sum thresholds; pair-type filters |
| O6 | Excluding ISS / mega-constellations as the paper does | We keep all objects | Add filters to separate constellation self-conjunctions from debris risk |
| O7 | Long-term cascade dynamics (Kessler), multi-year simulation | Window of hours to days | Population evolution model, decay, launch traffic |
| O8 | Realistic manoeuvre planning (fuel, timing, secondary conjunctions, operator coordination) | Vertex cover is a stylised model | Weighted cover with cost proxies; re-screen after each manoeuvre; multi-objective planning |
| O9 | Parallel / GPU screening | Single-threaded reference implementation | Thread-parallel step loop (needs a decision on platform threads), SIMD, GPU |
| O10 | Incremental data structures (updating the k-d tree or grid instead of rebuilding each step) | Rebuild is simple and fast enough | Incremental updates; Verlet-list style reuse across steps |
| O11 | Windows-specific verification of every script | Development and tests ran on Linux (device VM) | CI on Windows. `test.bat`, `build.bat` and `demo.bat` verified by hand on Windows 2026-10-25; the bench/experiment `.sh` scripts and CI are not |

## B. Open weaknesses (could still be improved)
| # | Item | Status | Plan |
|---|---|---|---|
| L11, L26 | Timings are single runs with `clock()` CPU time, one VM, no environment capture | OPEN, C7 | Repeats, median and range, `results/environment.txt` |
| L25 | After indexing, time is dominated by propagation and refinement | OPEN, C7 | Cheaper pre-refinement test, reuse endpoint states, measure each change |
| L28 | Exact vertex cover cannot finish on the giant cluster at 3 h | OPEN, optional | LP / Nemhauser-Trotter lower bound, per-node reductions; or lower thresholds (C7 sweep) |
| L5 | Epoch age: stale elements; fixed by a 3-day freshness filter, not by an uncertainty model | DOC | Age-aware screening distance (C7 experiment, unvalidated heuristic) |
| L23 | Altitude filter is empirical, not a proven bound; slower than brute force in LEO | DOC | Only useful for mixed LEO/MEO/GEO catalogues |
| L19 | SGP4 velocity is not the exact derivative of position (~1.5 m/s) | DOC | Hermite/finite-difference refinement using positions only |
| L20 | Event agreement is SGP4 vs SGP4: validates algorithm and port, not real-world accuracy | DOC | Compare with operator ephemerides (e.g. SOCRATES reports) where public |
| L9 | Synthetic truth exact only for the two-body model | DOC | Planted pairs under J2/SGP4 would need numerically solved truth |
| L6 | Pairs that stay close for the whole window (co-orbiting) have no interior minimum and are not events | DOC | Report them separately as "persistent proximity" (counted as `close_at_start`) |
| L27, L30, L31 | Mitigation model, movable set from file membership, single snapshot | DOC | See O7, O8; take object type from SATCAT |
| L13 | libm (`-lm`) and platform-dependent last digits | DOC | Hand-written sin/cos/Kepler if strict independence is needed |
| L12 | Julian date in a double resolves ~40 us | DOC | Two-part Julian date (integer + fraction) |
| L1 | Mercurio & Singla 2016 not accessible; cited from abstract only | DOC | Obtain full text, then cite details |
| L21 | Objects that decay mid-window drop out silently | DOC (warning printed at start only) | Count objects that fail during the window |

## C. Fixed (history)
Model error vs SGP4 (L3, L7, L8: hand-written SGP4), drag (L2), deep-space error (L4: scope), documentation inconsistencies (L15, L16, L17), DSA structures unused (L14), parked-object warning (L21 at start), Top-K dominated by co-moving objects (L29).

## D. Practicality ideas agreed for C7 (done in C7: threshold sweep, age-aware heuristic, repeats + environment; not done: cheaper refinement, weighted mitigation)
Threshold sweep (1, 2, 5, 10 km) with cover sizes; epoch-age-aware thresholds (heuristic); cheaper refinement; reproducible one-command run with environment capture and repeated timings; weighted mitigation (optional).

## E. New from C7
| ID | Item | Evidence | Idea |
|---|---|---|---|
| F1 | DONE (L33): the extra events of pure linear interpolation are co-moving pairs; `kd:hybrid` removes them | `results/c8_l33_interp.csv` | See F11 for what is left |
| F2 | Calibrate age growth g from TLE-to-TLE differences (L34) | `results/c7_age.csv` | Use public historical elements |
| F3 | Bisection output is identical for 5-60 s steps, so a larger step only adds cost; coarse step plus finer refinement is untested at < 5 s | `results/c7_step.csv` | Test 2 s and 1 s steps |
| F4 | Parallel propagation (threads) not tried; propagation dominates (L25) | `results/bench_c5_summary.csv` | OpenMP or pthreads study |

## F. New from C8 (verification and QA)
| ID | Item | Evidence | Idea |
|---|---|---|---|
| F5 | Windows portability: DONE for build, tests and demo (`scripts\test.bat` ALL TESTS PASSED, `scripts\demo.bat` matches Linux: 19,440 objects, 2,254 events, largest cluster 24, greedy 1,353 vs exact 1,346; 2026-10-25, Windows 10.0.26200, MSYS2 gcc 16.1.0). STILL OPEN: `long` is 32 bits there (L38), so replace `long` counters by `int64_t`; peak-memory line reads 0 | manual run 2026-10-25 | Replace `long` counters by `int64_t` and `PRId64` |
| F6 | `analyze_events` input hardening (L44): sort events on load, bound the cost of the time-index key bump, extend the time index beyond 23 days with 64-bit keys | `docs/known_limitations.md` L44 | Small change, not needed for the 1-6 h windows used |
| F7 | Mutation testing with generated (not hand-picked) mutants, e.g. a small C-aware mutator over all operators; the 39 hand-written mutants in `scripts/mutations.tsv` are a starting set | `docs/verification.md` section 4 | Larger sample, per-module score |
| F8 | Windows equivalents of `test_oom` and `cli_smoke` (GNU ld `--wrap` works with MinGW but is untested; the smoke test is bash) | `scripts/test.bat` | F5 build/test part is done; port these next |
| F9 | Measure coverage of `main.c` (CLI) and make `cli_smoke` run under the coverage build | `docs/verification.md` section 3 | Link `main.c` in `scripts/coverage.sh` |
| F10 (PARTLY DONE 2026-10-24: median 2.5 km/day measured, see L34; remaining: longer history, per-orbit-class and per-percentile policy) | True calibration of the age-aware distance (L34): obtain multi-epoch element-set history (Space-Track `gp_history` needs a free account, or a CelesTrak archive) and run `scripts/calibrate_age_from_history.py`; use a high percentile (e.g. 95) as a bound rather than the typical growth | `docs/known_limitations.md` L34 | Needs the data from the team |
| F11 | DONE as `kd:hybrid` (see L33). Remaining: re-run the scaling and threshold experiments with `kd:hybrid`, and re-check the 1 h / 3 h / 6 h cluster results, before making it the default | `results/c8_l33_interp.csv` | The hybrid was only evaluated on the 1 h, 5 km real-data setup |

