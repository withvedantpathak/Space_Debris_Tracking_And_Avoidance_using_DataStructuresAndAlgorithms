# Verification and QA (chunk C8)

Everything here is reproducible with the scripts named in each section. Numbers are copied from the output of those scripts
(measured 2026-10-22 on the development VM), not estimated.

## 1. What is run

| Layer | Script / file | What it checks |
|---|---|---|
| Unit tests | `tests/test_ds.c`, `test_catalog.c`, `test_propagate.c`, `test_sgp4.c`, `test_conjunction.c`, `test_screen.c`, `test_graph.c` | Example-based tests of every module (10,188 / 65 / 37 / 306 / 36 / 197 / 22,452 checks) |
| Property tests | `tests/test_property.c` (12,327 checks; groups `basic screen meta valid graph analyze`) | Structural invariants of heap, k-d tree, grid, DSU, AVL; CSV parser fuzzing with an accounting identity; loader range boundaries; scope filters against naive predicates; synthetic generator; Kepler solver; **120 random screening configurations** (kd, grid and altitude filter against the brute-force baseline); metamorphic relations (input permutation, threshold monotonicity, window prefix, step independence, determinism); parameter validation; 300 random graphs (exact cover = subset-enumeration optimum, reductions on and off); `analyze_events` on damaged files |
| Out-of-memory injection | `tests/test_oom.c` | The k-th allocation fails (once, then every allocation from k on): 8 scenarios, 6,550 injected failures; each must be reported as an error or give the unfailed result; no crash, no leak |
| Command line | `tests/cli_smoke.sh` (46 checks) | All commands; determinism; `-O0` = `-O2` = sanitizer build byte for byte; kd/grid screening write exactly the baseline's event file; refusals (non-zero exit, no crash) |
| Sanitizers | `scripts/test.sh` | Every test twice: plain `-O2`, and `-fsanitize=address,undefined,float-divide-by-zero,float-cast-overflow -fno-sanitize-recover=all` (AddressSanitizer includes LeakSanitizer) |
| Compiler warnings | `scripts/build.sh`, `scripts/test.sh` | `-std=c99 -Wall -Wextra -pedantic`: no warnings. All sources also clean under `-Wshadow -Wformat=2 -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes` at `-O0/-O2/-O3` |
| Static analysis | `gcc -fanalyzer` on every source file | 0 warnings (final run, after the fixes below) |
| Coverage | `scripts/coverage.sh` -> `results/c8_coverage.csv` | See section 3 |
| Mutation check | `scripts/mutation_check.sh` -> `results/c8_mutation.csv` | See section 4 |

`scripts/test.sh` runs the first five layers (about 3 minutes in total; a test name can be given to run one test:
`scripts/test.sh test_screen`). `scripts/test.bat` (Windows) runs the unit and property tests only; it was run on Windows on 2026-10-25 (Windows 10.0.26200, MSYS2 gcc 16.1.0) and
printed `ALL TESTS PASSED` (`docs/future_work.md` F5).

## 2. Defects found by this chunk (all fixed unless stated)

| ID | Defect | Found by |
|---|---|---|
| L36 | A NUL byte in a CSV row made the loader swallow the next row | parser fuzz |
| L37 | `avl_check` overflowed where `long` is 32 bits (Windows) | code review |
| L39 | No validation of screening parameters: step 0 or NaN looped forever, huge step counts wrapped | parameter test (hung for 80 s) |
| L40 | Grid screening returned an error for an empty catalogue with threshold 0 (cell size 0) | random differential sweep |
| L41 | `grid_build` left its hash uninitialised on early failure; `grid_free` then crashed | grid failure test |
| L42 | `analyze_events`: leaks on error paths, self-pair events aborted the analysis, TCA outside 0..2e6 s overflowed the 32-bit time key (undefined behaviour) | analyze robustness test + sanitizers |
| L43 | `topk_smallest` overflowed for K = INT_MAX | top-K extremes test |
| L45 | Out-of-memory: `hash_*` crashed after a failed init; `cat_add` left array and index inconsistent; `analyze_events` ignored about ten allocation results and a failed index insert | malloc-failure injection |
| L47 | CLI exit status 0 for unknown commands and unopenable files | CLI smoke test |
| L48 | Sanitizer builds did not abort on `float-cast-overflow` (flag `-fno-sanitize-recover=undefined` does not cover it) | mutation M35 |
| (unnumbered) | `main.c`: four unchecked `malloc` calls in `screen`/`fastscreen` | `gcc -fanalyzer` |

Also documented, not fixed: L38 (Windows `long` counters), L44 (analysis input assumptions), L46 (scope of OOM injection), L50 (equivalent mutants).
The existing results were re-checked after the fixes that touch the pipeline: the analysis output on the 1 h and 6 h real event files
(top-K, crossing top-K, cover, timeline, components) is byte-identical before and after.

## 3. Coverage (`results/c8_coverage.csv`)

| | Lines | Branches |
|---|---|---|
| Before C8 (7 test programs) | 97.9 % (1,478 of 1,510) | 79.5 % (1,287 of 1,618) |
| After C8 | 99.9 % (1,565 of 1,567) | 91.8 % (1,623 of 1,768) |
| Final (two `synth.c` fallbacks tested) | 100 % (1,578 of 1,578) | 92.0 % (1,635 of 1,778) |

After C8 the only unexecuted lines were two defensive `return` statements in `synth.c` (the rounding fallback of `pick` and the give-up of `make_pair`). Both are now reached by `tests/test_catalog.c` through two small hooks (`synth_test_pick`, `synth_test_plant`; the attempt limit of `make_pair` became a parameter, still 1,000 in normal use), so no library line is left unexecuted. `main.c` is not part of the figure (the
coverage build links the library modules); it is exercised by `tests/cli_smoke.sh` only. Branch coverage is lower than line
coverage mostly because of error branches (out-of-memory paths that can only be reached by injection, already counted here, and
defensive checks).

## 4. Mutation check (`results/c8_mutation.csv`)

39 hand-written mutants (`scripts/mutations.tsv`), each a single-line change: inverted comparisons, off-by-one bounds, dropped
statements, weaker or stronger pruning tests, removed error handling. The mutated source is compiled alone and linked with the
tests; a mutant is KILLED when any test fails, crashes or times out.

**36 of 39 mutants killed (92 %).** The three survivors are judged equivalent by inspection (not proven):
M03 (heap sift-down on `>` instead of `>=` changes only the order of equal keys), M19 (top-K tie rule; a later pass re-selects the
lowest indices) and M28 (an event exactly at the window end, measure zero for floating-point times).

The first run of the mutation check left six mutants alive (M04, M21, M22, M35, M37, M39) that exposed real gaps in the tests, which were then closed
(L48, L49): the union-by-size direction of the DSU, exact range boundaries of the CSV loader, catalogue consistency after a failed
insert, the result of `analyze_events` under allocation failure, and the `float-cast-overflow` abort.

Limits of this measure: the mutants were chosen by hand (not generated at random), cover the algorithmic modules and the
loader, and a score of 92 % says nothing about bugs of kinds that were not mutated.

## 5. What is NOT verified

* Windows: `long` is 32 bits there (L38 open). `build.bat`, `test.bat` and `demo.bat` were run on 2026-10-25 and pass; `test_oom`, `cli_smoke` and the sanitizers are Linux only (F8). gcc 16.1.0 prints two warnings in test code only (`tests/test_graph.c:180` unused variable `total`; `tests/test_property.c:125` `-Walloc-size-larger-than` in `grid_ok`); no failures.
* SGP4 accuracy against reality (L20): events are SGP4 against SGP4. Validation of our SGP4 port against the reference
  implementation is in `docs/propagation.md`.
* L33 is fixed by the hybrid refinement (`kd:hybrid`), which still loses 1 to 4 events at steps of 20 s and above; its CPU cost was not measured. L34 (age-aware distance) is uncalibrated; calibration from Space-Track history is pending.
* Out-of-memory injection fails only the project's own allocations and is Linux/GNU-ld only (L46).
* The altitude filter is a heuristic, not a proven bound (L23); the tests assert that it never invents events and that it is
  exact whenever its padding assumption holds, not that it finds every event.
