#!/usr/bin/env bash
# L33: hybrid refinement (linear interpolation for crossing pairs, bisection for slow pairs) on the same real-data setup as the
# interp step of scripts/run_c7.sh (1 h, 5 km, SGP4 + LEO scope, 3-day fresh elements). Needs results/c7_events_step5.csv (reference)
# and results/c7_events_lin<h>.csv (pure linear) from run_c7.sh. Writes results/c8_events_hyb<h>.csv, results/c8_hybrid_timing.csv,
# then runs scripts/l33_analysis.py (accuracy of linear and hybrid against the reference).
set -e
cd "$(dirname "$0")/.."
bash scripts/build.sh > /dev/null
R="data/raw/active.csv data/raw/cosmos-2251-debris.csv data/raw/fengyun-1c-debris.csv data/raw/iridium-33-debris.csv"
echo "step_s,method,events,cpu_seconds" > results/c8_hybrid_timing.csv
for h in 5 10 20 30 60; do
    for m in linear hybrid; do
        out=results/c7_events_lin$h.csv; [ $m = hybrid ] && out=results/c8_events_hyb$h.csv
        ./build/dlscms fastscreen kd:$m $out 3600 $h 5 12 3 6.4 $R > /tmp/c8_h.log
        sm=$(grep '^SUMMARY' /tmp/c8_h.log)
        echo "$h,$m,$(echo $sm | awk -F, '{print $11","$14}')" >> results/c8_hybrid_timing.csv
    done
    echo "step $h done"
done
python3 scripts/l33_analysis.py
