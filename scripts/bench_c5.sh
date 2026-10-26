#!/usr/bin/env bash
# C5 benchmark: brute force vs k-d tree vs grid vs altitude filter on synthetic catalogues (two-body, 1 h window, 10 s step, 10 km).
# Every fast method's event file must be byte-identical to the brute-force one ("identical" column).
# Usage: scripts/bench_c5.sh N [N ...]   -> appends rows to results/bench_c5.csv (CPU seconds from clock(), single thread)
set -e
cd "$(dirname "$0")/.."
bash scripts/build.sh > /dev/null
CSV=results/bench_c5.csv
[ -s $CSV ] || echo "dataset,method,n,window_s,step_s,threshold_km,pair_checks,candidates,events,pairs,cpu_seconds,identical_to_bruteforce" > $CSV
for n in "$@"; do
    f=data/synthetic/synth_$n.csv
    ./build/dlscms screen /tmp/c5_brute_$n.csv 3600 10 10 0 0 0 $f > /tmp/c5_brute_$n.log
    for m in brute kd grid alt; do
        if [ $m = brute ]; then log=/tmp/c5_brute_$n.log; out=/tmp/c5_brute_$n.csv
        else ./build/dlscms fastscreen $m /tmp/c5_${m}_$n.csv 3600 10 10 0 0 0 $f > /tmp/c5_${m}_$n.log; log=/tmp/c5_${m}_$n.log; out=/tmp/c5_${m}_$n.csv; fi
        same=no; cmp -s $out /tmp/c5_brute_$n.csv && same=yes
        echo "$(grep '^SUMMARY' $log | awk -F, -v m=$m -v same=$same '{print $2","m","$3","$4","$5","$6","$8","$9","$11","$12","$14","same}')" >> $CSV
    done
    tail -4 $CSV
done
