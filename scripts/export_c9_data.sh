#!/usr/bin/env bash
# C9 data exports for the figures that need raw numbers from the C code (all produced with `dlscms prop`):
#   results/c9_positions.csv   V07: ECI positions of every object at the TCA of the top crossing pair (reference epoch + tca_s)
#   results/c9_pair_curve.csv  V11: relative position/velocity of the top crossing pair of the 1 h run, 1 s steps around its TCA
# The reference epoch is the median epoch of the merged catalogue ("common" in `prop`), the same start the 1 h screening used;
# the script checks that the curve reproduces the miss distance stored in results/c6_1h_topk_crossing.csv.
set -e
cd "$(dirname "$0")/.."
bash scripts/build.sh > /dev/null
mkdir -p build/c9
RAW=data/raw
( cat $RAW/active.csv; for f in cosmos-2251-debris fengyun-1c-debris iridium-33-debris; do tail -n +2 $RAW/$f.csv; done ) > build/c9/merged.csv
read A B T M <<< "$(awk -F, 'NR==2{print $2,$3,$4,$5}' results/c6_1h_topk_crossing.csv)"
./build/dlscms prop build/c9/merged.csv $T 12 common > results/c9_positions.csv
echo "pair $A $B tca_s=$T miss_km=$M"
echo "pair,$A,$B,tca_s,$T,miss_km,$M" > results/c9_pair_meta.csv
echo "dt_s,dx,dy,dz,dvx,dvy,dvz" > results/c9_pair_curve.csv
T0=$(awk -v t="$T" 'BEGIN{printf "%d", t-30}')
for i in $(seq 0 60); do
  t=$((T0+i))
  ./build/dlscms prop build/c9/merged.csv $t 12 common | awk -F, -v A=$A -v B=$B -v t=$t '
    $1==A{ax=$2;ay=$3;az=$4;avx=$5;avy=$6;avz=$7} $1==B{bx=$2;by=$3;bz=$4;bvx=$5;bvy=$6;bvz=$7}
    END{printf "%d,%.6f,%.6f,%.6f,%.9f,%.9f,%.9f\n", t, bx-ax,by-ay,bz-az,bvx-avx,bvy-avy,bvz-avz}' >> results/c9_pair_curve.csv
done
wc -l results/c9_positions.csv results/c9_pair_curve.csv
