#!/usr/bin/env bash
# C4: exact brute-force screening on the synthetic sets (two-body model, matches the planted truth).
# Writes results/baseline_synth_N.csv (events), results/baseline_summary.csv (cost/pruning), results/baseline_planted_N.csv (check).
# Usage: scripts/run_baseline.sh [N window_s ...]   default: 1000 86400  5000 21600  20000 3600   (brute force at 20000 over 24 h ~ 40 min; longer windows scale linearly, see docs)
set -e
cd "$(dirname "$0")/.."
bash scripts/build.sh > /dev/null
[ -f results/baseline_summary.csv ] || echo "dataset,n,window_s,step_s,threshold_km,flags,pair_checks,candidates,refinements,events,pairs_with_events,close_at_start,cpu_seconds" > results/baseline_summary.csv
args=("$@"); [ ${#args[@]} -eq 0 ] && args=(1000 86400 5000 21600 20000 3600)
for ((i=0; i<${#args[@]}; i+=2)); do
    n=${args[i]}; w=${args[i+1]}
    ./build/dlscms screen results/baseline_synth_$n.csv $w 10 10 0 0 0 data/synthetic/synth_$n.csv | tee /tmp/screen_$n.log
    grep '^SUMMARY' /tmp/screen_$n.log | cut -d, -f2- >> results/baseline_summary.csv
    python3 scripts/verify_baseline.py $n $w 10
done
