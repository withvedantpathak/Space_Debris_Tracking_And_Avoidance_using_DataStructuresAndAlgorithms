#!/usr/bin/env bash
# Mutation check: does the test suite notice a deliberate bug?
# Each line of scripts/mutations.tsv is ONE small change to ONE source line. For every mutant this script compiles only the changed
# file (everything else is reused), relinks the tests (UBSan + float checks, -O1) and runs them until one fails.
#   KILLED   = some test failed, crashed or timed out
#   SURVIVED = all tests passed (either the tests have a gap, or the mutant does not change behaviour: see docs/verification.md)
# The sources are never modified: mutants are built from a copy in build/mut/.
# Replacement text <EMPTY> means "delete the old text".
# Usage: scripts/mutation_check.sh [first_index [last_index]]   (indices into mutations.tsv, 1-based; default all)
#        The result file results/c8_mutation.csv is rewritten when first_index is 1 (or absent) and appended to otherwise.
cd "$(dirname "$0")/.." || exit 2
M=build/mut; mkdir -p $M/obj $M/src results
FIRST=${1:-1}; LAST=${2:-9999}
FLAGS="-std=c99 -O1 -g -Isrc -fsanitize=undefined,float-divide-by-zero,float-cast-overflow -fno-sanitize-recover=all"
LIBSRC=$(ls src/*.c src/ds/*.c | grep -v 'src/main.c')
# baseline objects (rebuilt only when missing)
if [ ! -f $M/obj/.ready ] || [ -n "$(find src tests -name '*.c' -newer $M/obj/.ready | head -1)" ]; then   # rebuild when any source or test is newer than the cache
  for f in $LIBSRC; do o=$M/obj/$(echo $f | tr '/' '_').o; gcc $FLAGS -c $f -o $o || exit 2; done
  for t in test_ds test_catalog test_propagate test_conjunction test_graph test_screen test_sgp4 test_property; do gcc $FLAGS -c tests/$t.c -o $M/obj/T_$t.o || exit 2; done
  gcc $FLAGS -D_POSIX_C_SOURCE=200809L -c tests/test_oom.c -o $M/obj/T_test_oom.o || exit 2
  touch $M/obj/.ready
fi
[ "$FIRST" -le 1 ] && echo "id,file,description,result,killed_by,seconds" > results/c8_mutation.csv
RUNS_ALL=("test_ds" "test_catalog" "test_propagate" "test_conjunction" "test_graph" "test_screen" "test_sgp4" "test_property basic" "test_property graph" "test_property analyze" "test_property valid" "test_property meta" "test_property screen")
idx=0
grep -v '^#' scripts/mutations.tsv | while IFS=$'\t' read -r id file old new which desc; do
  idx=$((idx+1)); [ $idx -lt $FIRST ] && continue; [ $idx -gt $LAST ] && break
  start=$(date +%s); [ "$new" = "<EMPTY>" ] && new=""   # (an empty field would be swallowed by read)
  # apply the mutation to a copy of the file (exactly one occurrence required)
  python3 - "$file" "$old" "$new" $M/src/mut.c <<'PY' || { echo "$id,$file,\"$desc\",BAD_MUTATION,," >> results/c8_mutation.csv; echo "$id BAD_MUTATION (text not found exactly once)"; continue; }
import sys
f, old, new, out = sys.argv[1:5]
s = open(f).read()
if s.count(old) != 1: sys.exit(1)
open(out, "w").write(s.replace(old, new))
PY
  obj=$M/obj/$(echo $file | tr '/' '_').o; mo=$M/mut.o
  if ! gcc $FLAGS -Isrc/ds -Isrc -c $M/src/mut.c -o $mo 2>$M/cc.err; then
    echo "$id,$file,\"$desc\",BUILD_ERROR,," >> results/c8_mutation.csv; echo "$id BUILD_ERROR"; continue; fi
  OBJS=""; for f in $LIBSRC; do [ "$f" = "$file" ] && OBJS="$OBJS $mo" || OBJS="$OBJS $M/obj/$(echo $f | tr '/' '_').o"; done
  result=SURVIVED; by=""
  if [ "$which" = "oom" ]; then
    gcc $FLAGS -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc -o $M/run_oom $M/obj/T_test_oom.o $OBJS -lm 2>/dev/null || { result=KILLED; by="link"; }
    if [ "$result" = SURVIVED ]; then timeout 120 $M/run_oom >/dev/null 2>&1 || { result=KILLED; by="test_oom"; }; fi
  else
    for entry in "${RUNS_ALL[@]}"; do
      set -- $entry; t=$1; arg=$2
      gcc $FLAGS -o $M/run_$t $M/obj/T_$t.o $OBJS -lm 2>/dev/null || { result=KILLED; by="link $t"; break; }
      timeout 90 $M/run_$t $arg >/dev/null 2>&1 || { result=KILLED; by="$entry"; break; }
    done
  fi
  echo "$id,$file,\"$desc\",$result,$by,$(( $(date +%s) - start ))" >> results/c8_mutation.csv
  echo "$id $result $by"
done
