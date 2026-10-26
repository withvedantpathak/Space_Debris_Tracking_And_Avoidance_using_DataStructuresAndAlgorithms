@echo off
REM Build (and optionally run) the DLSCMS DSA program. Usage: scripts\build.bat [run]
cd /d "%~dp0\.."
if not exist build mkdir build
gcc -O2 -std=c99 -Wall -D__USE_MINGW_ANSI_STDIO=1 -Isrc -o build\dlscms.exe src\*.c src\ds\*.c -lm || exit /b 1
echo Built build\dlscms.exe
if "%1"=="run" build\dlscms.exe
