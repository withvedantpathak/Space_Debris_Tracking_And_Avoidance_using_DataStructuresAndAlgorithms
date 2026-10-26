@echo off
REM Demo: the whole pipeline on the real CelesTrak catalogue (1 h window, 5 km threshold, SGP4, LEO scope). About 30 s.
REM Output goes to build\demo\. Same commands as scripts/demo.sh.  Verified on Windows (2026-10-25).
cd /d "%~dp0\.."
call scripts\build.bat || exit /b 1
if not exist build\demo mkdir build\demo
set D=data\raw\active.csv data\raw\cosmos-2251-debris.csv data\raw\fengyun-1c-debris.csv data\raw\iridium-33-debris.csv
echo === 1. catalogue
build\dlscms.exe info %D% || exit /b 1
echo === 2. screening with the k-d tree (1 h, 10 s step, 5 km, SGP4 + LEO scope)
build\dlscms.exe fastscreen kd build\demo\events.csv 3600 10 5 12 3 6.4 %D% || exit /b 1
echo === 3. clusters, top-K, mitigation
build\dlscms.exe analyze build\demo\events.csv data\raw\active.csv 10 build\demo\out 300000 || exit /b 1
echo === files written to build\demo\
dir /b build\demo
