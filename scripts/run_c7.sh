#!/usr/bin/env bash
# C7 experiments. Each step is independent and writes its own CSV in results/. Run all in order, or one: scripts/run_c7.sh <step> [args]
#   env | scaling "<N list>" | threshold | step | cell | interp | age | summary
# Real data: CelesTrak files in data/raw, SGP4 + LEO scope (flags 12), 1 h window from the median epoch, fresh elements (3 days).
set -e
cd "$(dirname "$0")/.."
bash scripts/build.sh > /dev/null
R="data/raw/active.csv data/raw/cosmos-2251-debris.csv data/raw/fengyun-1c-debris.csv data/raw/iridium-33-debris.csv"
REAL="12 3 6.4"
val() { grep "$1" "$2" | head -1; }
cpu_of() { grep '^SUMMARY' "$1" | awk -F, '{print $14}'; }
rss_of() { grep '^PEAKRSS_KB' "$1" | cut -d, -f2; }
sum_of() { grep '^SUMMARY' "$1"; }   # SUMMARY,dataset,n,window,step,thr,flags,pair_checks,candidates,refinements,events,pairs,close_at_start,cpu
step_env() { bash scripts/capture_env.sh > /dev/null; echo "wrote results/environment.txt"; }

step_scaling() {   # NOTE: appends; to restart, truncate the file first (: > results/c7_scaling.csv)   # synthetic catalogues, 3 seeds, 3 timing repeats (brute force once for N >= 10000); 1 h, 10 s step, 10 km, two-body
    local CSV=results/c7_scaling.csv
    [ -s $CSV ] || echo "seed,n,method,repeat,cpu_seconds,pair_checks,candidates,events,peak_rss_kb,identical_to_bruteforce" > $CSV
    for n in $1; do for seed in 1 2 3; do
        f=data/synthetic/synth_$n.csv; [ $seed != 1 ] && f=data/synthetic/synth_${n}_s$seed.csv
        [ -f $f ] || ./build/dlscms gen $n $seed $f $((n/100)) > /dev/null
        if [ $n -ge 10000 ] && [ $seed != 1 ]; then continue; fi   # N >= 10000: seed 1 only
        ./build/dlscms screen /tmp/c7_ref.csv 3600 10 10 0 0 0 $f > /tmp/c7_ref.log
        for m in brute kd grid; do
            reps=3; [ $m = brute ] && [ $n -ge 10000 ] && reps=1
            for rep in $(seq 1 $reps); do
                if [ $m = brute ]; then ./build/dlscms screen /tmp/c7_out.csv 3600 10 10 0 0 0 $f > /tmp/c7_out.log
                else ./build/dlscms fastscreen $m /tmp/c7_out.csv 3600 10 10 0 0 0 $f > /tmp/c7_out.log; fi
                same=no; cmp -s /tmp/c7_out.csv /tmp/c7_ref.csv && same=yes
                sm=$(sum_of /tmp/c7_out.log)
                echo "$seed,$n,$m,$rep,$(echo $sm | awk -F, '{print $14","$8","$9","$11}'),$(rss_of /tmp/c7_out.log),$same" >> $CSV
            done
        done
        echo "scaling seed=$seed n=$n done"
    done; done
}

step_scale50k() {   # brief 11.1 largest point: N=50000, seed 1, kd and grid only (brute force at this size exceeds the 180 s per-call limit of the dev shell and was NOT run). Usage: scripts/run_c7.sh scale50k kd|grid|compare
    local CSV=results/c7_scaling_50k.csv f=data/synthetic/synth_50000.csv
    [ -f $f ] || ./build/dlscms gen 50000 1 $f 500 > /dev/null
    [ -s $CSV ] || echo "seed,n,method,repeat,cpu_seconds,pair_checks,candidates,events,peak_rss_kb" > $CSV
    case "$1" in
      kd|grid)
        ./build/dlscms fastscreen $1 results/c7_events_50k_$1.csv 3600 10 10 0 0 0 $f > /tmp/c7_50k.log
        sm=$(sum_of /tmp/c7_50k.log)
        echo "1,50000,$1,1,$(echo $sm | awk -F, '{print $14","$8","$9","$11}'),$(rss_of /tmp/c7_50k.log)" >> $CSV; tail -1 $CSV ;;
      compare) cmp results/c7_events_50k_kd.csv results/c7_events_50k_grid.csv && echo "kd and grid events identical" ;;
    esac
}

step_threshold() {   # headline sweep: screening distance vs conjunction pairs, clusters and mitigation (real data)
    local CSV=results/c7_threshold.csv
    echo "threshold_km,events,pairs,objects,components,largest_cluster,greedy_cover,exact_cover,exact_finished_all,lower_bound,milp_optimum,fastscreen_cpu_s" > $CSV
    for t in 1 2 5 10; do
        ./build/dlscms fastscreen kd results/c7_events_thr$t.csv 3600 10 $t $REAL $R > /tmp/c7_thr.log
        ./build/dlscms analyze results/c7_events_thr$t.csv data/raw/active.csv 10 results/c7_thr$t 300000 > /tmp/c7_an.log
        python3 scripts/verify_cover.py results/c7_events_thr$t.csv data/raw/active.csv results/c7_thr$t > /tmp/c7_ver.log || true
        c=$(grep '^C6SUMMARY' /tmp/c7_an.log | cut -d, -f3-)   # events,pairs,nodes,movable,components,single,largest_nodes,largest_pairs,agree,topk,greedy,exact,lb,unres,nontrivial,done,worse,tg,tx
        IFS=, read ev pairs nodes mov comps single ln lp agree topk g x lb unres nt dn worse tg tx <<< "$c"
        milp=$(awk -F, 'NR==2{print $3}' results/c7_thr${t}_verify.csv)
        all=$([ "$nt" = "$dn" ] && echo yes || echo no)
        echo "$t,$ev,$pairs,$nodes,$comps,$ln,$g,$x,$all,$lb,$milp,$(cpu_of /tmp/c7_thr.log)" >> $CSV
        echo "threshold $t km done"
    done
    cat $CSV
}

step_step() {   # sensitivity to the time step (bisection refinement), 5 km, kd
    local CSV=results/c7_step.csv
    echo "step_s,events,pairs,cpu_seconds,candidates,recall_vs_5s,precision_vs_5s,miss_abs_diff_median_km,miss_abs_diff_p95_km,tca_abs_diff_median_s" > $CSV
    for h in 5 10 20 30 60; do
        ./build/dlscms fastscreen kd results/c7_events_step$h.csv 3600 $h 5 $REAL $R > /tmp/c7_s.log
        sm=$(sum_of /tmp/c7_s.log); cmpv=$(python3 scripts/compare_events.py results/c7_events_step5.csv results/c7_events_step$h.csv 2 | cut -d, -f6-)
        echo "$h,$(echo $sm | awk -F, '{print $11","$12","$14","$9}'),$cmpv" >> $CSV; echo "step $h s done"
    done
    cat $CSV
}

step_interp() {   # E9: linear interpolation of positions (Gomez et al.) vs bisection, against the 5 s bisection reference
    local CSV=results/c7_interp.csv
    echo "step_s,method,events,pairs,cpu_seconds,recall_vs_ref,precision_vs_ref,miss_abs_diff_median_km,miss_abs_diff_p95_km,tca_abs_diff_median_s" > $CSV
    for h in 5 10 20 30 60; do
        ./build/dlscms fastscreen kd:linear results/c7_events_lin$h.csv 3600 $h 5 $REAL $R > /tmp/c7_l.log
        sm=$(sum_of /tmp/c7_l.log); cmpv=$(python3 scripts/compare_events.py results/c7_events_step5.csv results/c7_events_lin$h.csv 2 | cut -d, -f6-)
        echo "$h,linear,$(echo $sm | awk -F, '{print $11","$12","$14}'),$cmpv" >> $CSV
        b=$(awk -F, -v h=$h 'NR>1 && $1==h{print $2","$3","$4","$6","$7","$8","$9","$10}' results/c7_step.csv); echo "$h,bisection,$b" >> $CSV; echo "interp $h s done"
    done
    cat $CSV
}

step_cell() {      # grid cell size vs query reach, and memory of kd vs grid (10 min window, 5 km)
    local CSV=results/c7_cell.csv
    echo "method,cell_scale,cpu_seconds,peak_rss_kb,pair_checks,events" > $CSV
    for spec in kd grid:0.25 grid:0.5 grid:1 grid:2 grid:4; do
        ./build/dlscms fastscreen $spec /tmp/c7_c.csv 600 10 5 $REAL $R > /tmp/c7_c.log
        sm=$(sum_of /tmp/c7_c.log); sc=${spec#*:}; [ "$spec" = kd ] && sc=NA
        echo "${spec%%:*},$sc,$(echo $sm | awk -F, '{print $14","}')$(rss_of /tmp/c7_c.log),$(echo $sm | awk -F, '{print $8","$11}')" >> $CSV
    done
    cat $CSV
}

step_age() {       # epoch-age-aware screening distance (heuristic, unvalidated): +g km per day of element age per object
    local CSV=results/c7_age.csv
    echo "growth_km_per_day,events,pairs,objects,components,largest_cluster,greedy_cover,exact_cover" > $CSV
    for g in 0 1 2 3; do
        DLSCMS_AGE_GROWTH=$g ./build/dlscms fastscreen kd results/c7_events_age$g.csv 3600 10 5 $REAL $R > /tmp/c7_a.log
        ./build/dlscms analyze results/c7_events_age$g.csv data/raw/active.csv 10 results/c7_age$g 300000 > /tmp/c7_an.log
        c=$(grep '^C6SUMMARY' /tmp/c7_an.log | cut -d, -f3-)
        IFS=, read ev pairs nodes mov comps single ln lp agree topk gc x lb unres nt dn worse tg tx <<< "$c"
        echo "$g,$ev,$pairs,$nodes,$comps,$ln,$gc,$x" >> $CSV; echo "age growth $g done"
    done
    cat $CSV
}

case "$1" in
    env) step_env;; scaling) step_scaling "$2";; threshold) step_threshold;; step) step_step;; interp) step_interp;; cell) step_cell;; age) step_age;; scale50k) step_scale50k "$2";;
    summary) python3 scripts/make_summary.py;;
    *) echo "usage: $0 env|scaling \"<N list>\"|threshold|step|interp|cell|age|scale50k kd|grid|compare|summary"; exit 2;;
esac
