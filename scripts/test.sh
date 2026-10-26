#!/usr/bin/env bash
# Build and run all unit tests, plain and with AddressSanitizer + UndefinedBehaviorSanitizer
# (+ float-divide-by-zero and float-cast-overflow, which -fsanitize=undefined does not include).
# Usage: scripts/test.sh [test_name ...]   (no argument = all tests; run from anywhere)
set -e
cd "$(dirname "$0")/.."
mkdir -p build
LIB="src/catalog.c src/synth.c src/propagate.c src/sgp4.c src/conjunction.c src/screen.c src/graph.c src/mitigate.c src/rank.c src/analyze.c src/ds/*.c"
WARN="-std=c99 -Wall -Wextra -pedantic -Isrc"
SAN="-fsanitize=address,undefined,float-divide-by-zero,float-cast-overflow -fno-sanitize-recover=all"
TESTS="test_ds test_catalog test_propagate test_sgp4 test_conjunction test_screen test_graph test_property test_oom"
[ $# -gt 0 ] && TESTS="$*"
for t in $TESTS; do
    EXTRA=""
    # test_oom injects malloc failures through GNU ld --wrap (Linux / MinGW); it needs POSIX for freopen /dev/null
    [ "$t" = "test_oom" ] && EXTRA="-D_POSIX_C_SOURCE=200809L -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc"
    echo "== $t: plain -O2 build"
    gcc -O2 $WARN $EXTRA -o build/$t tests/$t.c $LIB -lm
    ./build/$t
    echo "== $t: ASan + UBSan + float checks build"
    gcc -O1 -g $WARN $SAN $EXTRA -o build/${t}_san tests/$t.c $LIB -lm
    ./build/${t}_san
done
# command-line smoke test (builds -O2, -O0 and sanitizer binaries itself); skipped when test names were given on the command line
if [ $# -eq 0 ]; then echo "== cli_smoke"; tests/cli_smoke.sh; fi
echo "ALL TESTS PASSED"
