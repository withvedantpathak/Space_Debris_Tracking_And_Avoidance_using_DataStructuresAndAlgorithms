# Running the program

There is **no graphical interface**: `dlscms` is a command-line program (C99, no libraries). The "interface" is the command line;
results are CSV files plus a printed summary. 

## Windows (no Linux needed)

The source is plain C99 and builds with the MinGW-w64 `gcc`. **Verified on Windows 2026-10-25** (Windows 10.0.26200, MSYS2 gcc 16.1.0): `scripts\test.bat` prints
`ALL TESTS PASSED` and `scripts\demo.bat` gives the expected output below. Development and the other experiments ran on Linux (`docs/future_work.md` F5).

1. **Install a C compiler** (one-off). In PowerShell: `winget install BrechtSanders.WinLibs.POSIX.UCRT` (a MinGW-w64 gcc), then close
   and reopen the terminal. Alternative: install MSYS2 and run `pacman -S mingw-w64-ucrt-x86_64-gcc` inside its UCRT64 shell.
2. **Check it works**: open *Command Prompt*, run `gcc --version`. If it says "not recognised", the install folder is not on PATH.
3. **Open the project folder**: `cd /d C:\Vedant\NMIMS\projects\DSA`.
4. **Run the tests** (proves the build on your machine): `scripts\test.bat`. Expected last line: `ALL TESTS PASSED`.
   Please send me the output if anything fails.
5. **Build**: `scripts\build.bat` (creates `build\dlscms.exe`).
6. **Run the demo** (real catalogue, about 30 s): `scripts\demo.bat`. Expected: 19,440 objects, 2,254 events on 2,204 pairs,
   largest cluster 24 objects, greedy cover 1,353 vs exact 1,346. Files appear in `build\demo\` (open the `.csv` files in Excel).

## Commands (Windows spelling; on Linux/macOS use `./build/dlscms` and `/`)

| Task | Command |
|---|---|
| Describe catalogues | `build\dlscms.exe info data\raw\active.csv` |
| Make a synthetic catalogue (n, seed, output, planted pairs) | `build\dlscms.exe gen 5000 1 build\synth5000.csv 50` |
| Brute-force screening | `build\dlscms.exe screen out.csv <window_s> <step_s> <threshold_km> <flags> <fresh_days> <min_mean_motion> <catalog.csv> [more.csv ...]` |
| Fast screening (k-d tree, grid, altitude filter) | `build\dlscms.exe fastscreen kd out.csv <same arguments>` (also `grid`, `grid:0.5`, `alt`; `kd:hybrid` = faster closest-approach refinement, equal to the default in tests, recommended for long windows; `kd:linear` = pure linear interpolation, kept for comparison only, it reports extra events for co-moving pairs) |
| Clusters, top-K, mitigation | `build\dlscms.exe analyze events.csv data\raw\active.csv 10 prefix 300000` (`-` instead of the CSV = every object movable) |

`flags`: 1 = J2, 2 = drag, 4 = LEO scope (175-2000 km), 8 = SGP4; the real-data runs use 12 (SGP4 + LEO scope) with `3 6.4`
(elements at most 3 days from the median epoch, mean motion at least 6.4 rev/day).
Invalid parameters are refused with a message and exit status 1 or 2. Run `build\dlscms.exe` with no arguments for the usage text.

## Linux / macOS

`scripts/build.sh`, `scripts/demo.sh`, `scripts/test.sh` (adds the sanitizers), `tests/cli_smoke.sh`.

## Known differences on Windows

`long` is 32 bits there (L38: counters overflow only above 2.1e9, i.e. very large brute-force runs); the line `PEAKRSS_KB` (peak
memory) reads `/proc` and is therefore 0 or missing; the shell scripts (`*.sh`, `run_c7.sh`, ...) need Git Bash or WSL; the
Python evidence scripts need `pip install sgp4 pandas matplotlib scipy`.

## Calibrating the age-aware distance (L34, needs a free Space-Track account)
1. Register at space-track.org (free) and note your username and password. Never put them in a file or in chat.
2. Windows (cmd), in the project folder: `set SPACETRACK_USER=your_email` then `set SPACETRACK_PASS=your_password` (valid for that window only). Linux/macOS: `export` the same names.
3. Check the plan first (no network): `python scripts/fetch_spacetrack_history.py --dry-run`
4. Download (about 6 requests, about 1 minute): `python scripts/fetch_spacetrack_history.py`  -> `data/history/spacetrack_gp_history.csv`
5. Calibrate: `python scripts/calibrate_age_from_history.py data/history/spacetrack_gp_history.csv 5 95`  -> `results/c8_tle_error_growth.csv` and a suggested g in km/day.
The fetch script has not been tested against the live service; if it stops with an error message, send me the message (not your password).
The downloaded history is for this project only; check Space-Track's terms before sharing it.

