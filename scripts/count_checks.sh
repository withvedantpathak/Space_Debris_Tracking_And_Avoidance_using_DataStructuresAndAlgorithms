#!/usr/bin/env bash
# Records how many automated checks each test program runs (plain -O2 build) in results/c8_test_counts.csv (used by V28).
# Usage: scripts/count_checks.sh reset | <test_program> ... | cli_smoke      (one program per call is fine; the dev shell has a 180 s limit)
cd "$(dirname "$0")/.."
F=results/c8_test_counts.csv
if [ "$1" = reset ]; then echo "program,checks,failed" > $F; exit 0; fi
[ -s $F ] || echo "program,checks,failed" > $F
for t in "$@"; do
  sed -i "/^$t,/d" $F
  if [ "$t" = cli_smoke ]; then line=$(bash tests/cli_smoke.sh | grep -E "^cli_smoke: [0-9]+ checks" | tail -1); n=$(echo "$line" | sed -E 's/.* ([0-9]+) checks, ([0-9]+) failed/\1,\2/')
  else line=$(bash scripts/test.sh $t 2>&1 | awk '/plain/{p=1;next} p&&/^[0-9]+ checks, [0-9]+ failed/{print;exit}'); n=$(echo "$line" | sed -E 's/^([0-9]+) checks, ([0-9]+) failed.*/\1,\2/'); fi
  echo "$t,$n" >> $F; echo "$t,$n"
done
