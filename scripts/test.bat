@echo off
REM Build and run the unit and property tests (no sanitizers on Windows MinGW; use scripts/test.sh for those).
REM NOT included: tests\test_oom.c (needs GNU ld --wrap and POSIX freopen) and tests\cli_smoke.sh (bash); run those with scripts/test.sh.
REM Verified on Windows (2026-10-25, Windows 10.0.26200, MSYS2 gcc 16.1.0): ALL TESTS PASSED.
cd /d "%~dp0\.."
if not exist build mkdir build
for %%t in (test_ds test_catalog test_propagate test_sgp4 test_conjunction test_screen test_graph test_property) do (
    gcc -O2 -std=c99 -Wall -Wextra -pedantic -D__USE_MINGW_ANSI_STDIO=1 -Isrc -o build\%%t.exe tests\%%t.c src\catalog.c src\synth.c src\propagate.c src\sgp4.c src\conjunction.c src\screen.c src\graph.c src\mitigate.c src\rank.c src\analyze.c src\ds\*.c -lm || exit /b 1
    build\%%t.exe || exit /b 1
)
echo ALL TESTS PASSED
