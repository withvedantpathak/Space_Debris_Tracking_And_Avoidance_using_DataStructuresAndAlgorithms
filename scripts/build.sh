#!/usr/bin/env bash
# Build (and optionally run) the DLSCMS DSA program.  Usage: scripts/build.sh [run]
set -e
cd "$(dirname "$0")/.."
mkdir -p build
gcc -O2 -std=c99 -Wall -Isrc -o build/dlscms src/*.c src/ds/*.c -lm
echo "Built build/dlscms"
[ "$1" = "run" ] && ./build/dlscms
exit 0
