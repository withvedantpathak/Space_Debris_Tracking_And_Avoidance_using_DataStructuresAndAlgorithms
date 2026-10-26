#!/usr/bin/env bash
# Records the machine and tool versions the results were produced with -> results/environment.txt
cd "$(dirname "$0")/.."
{
echo "captured: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "kernel: $(uname -srm)"
echo "cpu: $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ //')"
echo "logical cpus: $(nproc 2>/dev/null)"
echo "memory: $(grep -m1 MemTotal /proc/meminfo 2>/dev/null | sed 's/  */ /g')"
echo "compiler: $(gcc --version | head -1)"
echo "build flags: gcc -O2 -std=c99 -Wall -Isrc ... -lm (scripts/build.sh); timing is clock() CPU time, single thread"
echo "python: $(python3 --version 2>&1)"
python3 - <<'PY'
import importlib
for m in ["numpy", "pandas", "scipy", "matplotlib", "sgp4"]:
    try: print(f"{m}: {importlib.import_module(m).__version__}")
    except Exception as e: print(f"{m}: not installed")
PY
echo "note: the machine is a virtual machine on the project owner's laptop (2 logical CPUs visible); absolute times are indicative, ratios are the result"
} > results/environment.txt
cat results/environment.txt
