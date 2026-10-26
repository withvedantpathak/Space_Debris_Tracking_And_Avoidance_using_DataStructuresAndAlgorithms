#!/usr/bin/env bash
# C8: line and branch coverage of src/ by the unit and property tests (gcc --coverage + gcov).
# Writes results/c8_coverage.csv (one row per source file plus TOTAL) and, with --list, prints the
# source lines that no test executes.  Usage: scripts/coverage.sh [--list]
# Source objects are compiled once with coverage counters; every tests/test_*.c is linked against the same
# objects, so the counters accumulate across all test programs.  Needs gcc and gcov.
set -e
cd "$(dirname "$0")/.."
ROOT=$(pwd)
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
mkdir -p "$T/obj" "$T/gcov" results
SRCS="src/catalog.c src/synth.c src/propagate.c src/sgp4.c src/conjunction.c src/screen.c src/graph.c src/mitigate.c src/rank.c src/analyze.c src/ds/*.c"
for f in $SRCS; do
    gcc -O0 -g -std=c99 -I"$ROOT/src" --coverage -c "$ROOT/$f" -o "$T/obj/$(echo "$f" | tr '/' '_' | sed 's/\.c$/.o/')"
done
for t in tests/test_*.c; do
    n=$(basename "$t" .c)
    CF=""; LF=""
    [ "$n" = "test_oom" ] && { CF="-D_POSIX_C_SOURCE=200809L"; LF="-Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc"; }   # malloc-failure injection
    gcc -O0 -g -std=c99 -Isrc $CF -c "$t" -o "$T/$n.o"
    gcc --coverage -o "$T/$n" "$T/$n.o" "$T"/obj/*.o $LF -lm
    echo "== running $n"
    "$T/$n" > "$T/$n.log" 2>&1 || { cat "$T/$n.log"; echo "test $n FAILED under coverage build"; exit 1; }
    tail -1 "$T/$n.log"
done
cd "$T/gcov"
for o in "$T"/obj/*.o; do
    base="${o%.o}"
    gcov -b -o "$T/obj" "$base.gcda" > /dev/null 2>&1 || true
done
cd - > /dev/null
python3 - "$T/gcov" results/c8_coverage.csv "$1" <<'PY'
import glob, os, re, sys
d, out, flag = sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else ""
rows, tot = [], [0, 0, 0, 0]
for g in sorted(glob.glob(os.path.join(d, "*.gcov"))):
    src = None; lines = hit = br = brhit = 0; missing = []
    for ln in open(g, errors="replace"):
        if ln.startswith("        -:    0:Source:"):
            src = ln.split("Source:", 1)[1].strip()
            continue
        m = re.match(r"\s*([^:]+):\s*(\d+):(.*)", ln)
        if m:
            cnt, no = m.group(1).strip(), int(m.group(2))
            if no == 0 or cnt == "-": continue
            lines += 1
            if cnt.startswith("#####") or cnt.startswith("====="): missing.append((no, m.group(3)))
            else: hit += 1
        elif ln.startswith("branch"):
            br += 1
            if "taken 0" not in ln and "never executed" not in ln: brhit += 1
    m2 = re.search(r"(src/.*)$", src or "")
    if not m2: continue
    src = m2.group(1)
    rows.append((src, lines, hit, br, brhit, missing))
    tot[0] += lines; tot[1] += hit; tot[2] += br; tot[3] += brhit
with open(out, "w") as f:
    f.write("file,lines,lines_executed,line_pct,branches,branches_taken,branch_pct\n")
    for src, l, h, b, bh, _ in rows:
        f.write("%s,%d,%d,%.1f,%d,%d,%.1f\n" % (src, l, h, 100.0 * h / max(l, 1), b, bh, 100.0 * bh / max(b, 1)))
    f.write("TOTAL,%d,%d,%.1f,%d,%d,%.1f\n" % (tot[0], tot[1], 100.0 * tot[1] / max(tot[0], 1), tot[2], tot[3], 100.0 * tot[3] / max(tot[2], 1)))
print(open(out).read())
if flag == "--list":
    for src, l, h, b, bh, miss in rows:
        if miss:
            print("--- %s: %d unexecuted lines" % (src, len(miss)))
            for no, txt in miss: print("%5d:%s" % (no, txt.rstrip()))
PY
