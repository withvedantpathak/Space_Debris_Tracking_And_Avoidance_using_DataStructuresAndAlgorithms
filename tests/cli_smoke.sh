#!/usr/bin/env bash
# CLI smoke test: builds the program three ways (-O2, -O0, ASan+UBSan), runs the real command line on a synthetic catalogue and checks
#   - every command succeeds and prints what it should; bad input is refused with a non-zero exit and a message, never a crash
#   - determinism: two runs give byte-identical event files
#   - optimisation level / sanitizer build does not change the result (-O0 == -O2 == ASan, byte for byte)
#   - kd-tree and grid screening write exactly the baseline's event file
#   - the analysis stage runs on the events and writes its files
# Usage: tests/cli_smoke.sh    (from anywhere; uses build/smoke/)
cd "$(dirname "$0")/.." || exit 2
D=build/smoke; mkdir -p $D
fail=0; n=0
ok()  { n=$((n+1)); }
bad() { n=$((n+1)); fail=$((fail+1)); echo "FAIL: $1"; }
expect_ok()   { ok; if ! "$@" >$D/out.txt 2>$D/err.txt; then bad "exit status of: $*"; fi; }
expect_fail() { ok; "$@" >$D/out.txt 2>$D/err.txt; rc=$?; if [ $rc -eq 0 ]; then bad "should have failed: $*"; elif [ $rc -ge 124 ]; then bad "crashed or timed out (rc $rc): $*"; fi; }
same() { ok; cmp -s "$1" "$2" || bad "$3 ($1 differs from $2)"; }
has()  { ok; grep -q -- "$2" "$1" || bad "$1 should contain '$2': $3"; }

SRC="src/*.c src/ds/*.c"
gcc -O2 -std=c99 -Wall -Isrc -o $D/dl_o2 $SRC -lm || exit 2
gcc -O0 -g -std=c99 -Wall -Isrc -o $D/dl_o0 $SRC -lm || exit 2
gcc -O1 -g -std=c99 -Wall -Isrc -fsanitize=address,undefined,float-divide-by-zero,float-cast-overflow -fno-sanitize-recover=all -o $D/dl_san $SRC -lm || exit 2

B=$D/dl_o2
expect_ok $B gen 300 7 $D/cat.csv 3;                 has $D/out.txt "" ""
expect_ok $B info $D/cat.csv;                          has $D/out.txt "300" "info should report the objects"
expect_ok $B prop $D/cat.csv 60 1
args="600 10 20 1 1000 1 $D/cat.csv"                    # window step threshold flags(J2) fresh_days min_mean_motion catalogue
for v in o2 o0 san; do expect_ok $D/dl_$v screen $D/ev_$v.csv $args; done
same $D/ev_o2.csv $D/ev_o0.csv  "-O0 and -O2 events differ"
same $D/ev_o2.csv $D/ev_san.csv "sanitizer build and -O2 events differ"
expect_ok $B screen $D/ev_again.csv $args;             same $D/ev_o2.csv $D/ev_again.csv "two runs differ (not deterministic)"
ok; [ "$(wc -l < $D/ev_o2.csv)" -gt 1 ] || bad "screening found no events (test is vacuous)"
for m in kd grid; do expect_ok $B fastscreen $m $D/ev_$m.csv $args; same $D/ev_o2.csv $D/ev_$m.csv "$m screening differs from the baseline"; done
expect_ok $B fastscreen grid:0.5 $D/ev_g05.csv $args;  same $D/ev_o2.csv $D/ev_g05.csv "grid:0.5 differs from the baseline"
expect_ok $B fastscreen kd:hybrid $D/ev_hyb.csv $args; ok; [ "$(wc -l < $D/ev_hyb.csv)" -gt 1 ] || bad "hybrid screening found no events"
expect_fail $B fastscreen kd:bogus $D/x.csv $args
expect_ok $B analyze $D/ev_o2.csv - 5 $D/an 100000;    has $D/out.txt "C6SUMMARY" "analyze summary"
ok; [ -s $D/an_topk.csv ] && [ -s $D/an_components.csv ] && [ -s $D/an_cover.csv ] && [ -s $D/an_timeline.csv ] || bad "analyze output files missing"
expect_ok $D/dl_san analyze $D/ev_o2.csv - 5 $D/an_san 100000;  same $D/an_topk.csv $D/an_san_topk.csv "analyze top-K differs between builds"

# refusals: non-zero exit, a message, no crash
expect_ok $B; has $D/out.txt "usage" "no arguments"                        # bare program: banner + usage, exit 0 (scripts/build.sh run relies on it)
expect_fail $B info                                                         # command without its arguments: usage on stderr, exit 2
expect_fail $B nonsense
expect_fail $B info $D/no_such_file.csv
expect_fail $B screen $D/x.csv 600 0 20 1 1000 1 $D/cat.csv;               has $D/err.txt "invalid parameters" "step 0"
expect_fail $B screen $D/x.csv 600 -5 20 1 1000 1 $D/cat.csv
expect_fail $B screen $D/x.csv -600 10 20 1 1000 1 $D/cat.csv
expect_fail $B screen $D/x.csv 600 10 -1 1 1000 1 $D/cat.csv
expect_fail $B screen $D/x.csv 600 10 20 1 1000 1 $D/no_such_file.csv
expect_fail $B screen $D/x.csv 600
expect_fail $B fastscreen bogus $D/x.csv $args
expect_fail $B fastscreen grid:0.0001 $D/x.csv $args
expect_fail $B analyze $D/no_such_events.csv - 5 $D/an 10
expect_fail $B gen abc 1 $D/x.csv
: > $D/empty.csv;           expect_fail $B info $D/empty.csv
printf 'garbage\n\x00\x01\x02\n' > $D/garbage.csv; expect_fail $B info $D/garbage.csv; has $D/out.txt "HEADER ERROR" "garbage file"   # reported, not a crash
echo "cli_smoke: $n checks, $fail failed"
[ $fail -eq 0 ]
