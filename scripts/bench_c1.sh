#!/usr/bin/env bash
# C1 micro-benchmark -> results/bench_c1.csv (+ results/bench_c1.png if matplotlib is available).
set -e
cd "$(dirname "$0")/.."
mkdir -p build results
gcc -O2 -std=c99 -Wall -Wextra -pedantic -Isrc -o build/bench_spatial tests/bench_spatial.c src/ds/*.c
./build/bench_spatial 50 > results/bench_c1.csv
echo "wrote results/bench_c1.csv"
python3 scripts/plot_bench_c1.py results/bench_c1.csv results/bench_c1.png
