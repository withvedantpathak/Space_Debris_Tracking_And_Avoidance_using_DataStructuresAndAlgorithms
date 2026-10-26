#!/usr/bin/env bash
# Demo: the whole pipeline on the real CelesTrak catalogue (1 h window, 5 km threshold, SGP4, LEO scope). About 30 s.
# Output goes to build/demo/. Same commands as scripts\demo.bat.
set -e
cd "$(dirname "$0")/.."
bash scripts/build.sh > /dev/null
mkdir -p build/demo
D="data/raw/active.csv data/raw/cosmos-2251-debris.csv data/raw/fengyun-1c-debris.csv data/raw/iridium-33-debris.csv"
echo "=== 1. catalogue"
./build/dlscms info $D
echo "=== 2. screening with the k-d tree (1 h, 10 s step, 5 km, SGP4 + LEO scope)"
./build/dlscms fastscreen kd build/demo/events.csv 3600 10 5 12 3 6.4 $D
echo "=== 3. clusters, top-K, mitigation"
./build/dlscms analyze build/demo/events.csv data/raw/active.csv 10 build/demo/out 300000
echo "=== files written to build/demo/"; ls build/demo
