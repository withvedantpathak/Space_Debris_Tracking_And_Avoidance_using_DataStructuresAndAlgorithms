#!/usr/bin/env bash
# C6: conjunction graph, clusters, Top-K, mitigation on real events (SGP4, LEO scope, 10 km) of different window lengths.
# Needs results/c6_events_{1h,3h,6h}.csv, made with:
#   build/dlscms fastscreen kd results/c6_events_<w>.csv <seconds> 10 10 12 3 6.4 data/raw/active.csv data/raw/cosmos-2251-debris.csv data/raw/fengyun-1c-debris.csv data/raw/iridium-33-debris.csv
# Movable objects = those in data/raw/active.csv (active satellites); debris cannot manoeuvre.
# Writes results/c6_<w>_{components,topk,cover}.csv, results/c6_summary.csv; ablation without reductions for 1h and 3h.
set -e
cd "$(dirname "$0")/.."
bash scripts/build.sh > /dev/null
S=results/c6_summary.csv
echo "variant,events_file,events,pairs,nodes,movable_nodes,components,single_pair_components,largest_nodes,largest_pairs,algorithms_agree,topk,greedy_cover,exact_cover,lower_bound,unresolvable_pairs,nontrivial_components,exact_finished,greedy_worse_components,greedy_seconds,exact_seconds" > $S
for w in 1h 3h 6h; do
    ./build/dlscms analyze results/c6_events_$w.csv data/raw/active.csv 20 results/c6_$w 300000 | tee /tmp/c6_$w.log | head -5
    grep '^C6SUMMARY' /tmp/c6_$w.log | sed 's/^C6SUMMARY,/reductions,/' >> $S
    python3 scripts/verify_cover.py results/c6_events_$w.csv data/raw/active.csv results/c6_$w
done
for w in 1h 3h; do
    ./build/dlscms analyze results/c6_events_$w.csv data/raw/active.csv 20 /tmp/c6_nr_$w 100000 noreduce > /tmp/c6_nr_$w.log
    grep '^C6SUMMARY' /tmp/c6_nr_$w.log | sed 's/^C6SUMMARY,/branch_and_bound_only,/' >> $S
done
cut -d, -f1,3-4,13-14,15,17-18,21 $S
