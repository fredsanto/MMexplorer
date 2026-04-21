@echo off
:: Protein Missense Mutation Explorer — Windows build script
:: Copyright (c) 2026 fredsanto (Federico Santoni)
:: Requires: conda env with openmm/pdbfixer/biopython, PyInstaller, Inno Setup 6

setlocal enabledelayedexpansion

echo ============================================================
echo  MMexplorer — Windows installer build
echo ============================================================

:: ---- 1. Install PyInstaller if missing ----
python -m pip show pyinstaller >nul 2>&1
if errorlevel 1 (
    echo Installing PyInstaller...
    python -m pip install pyinstaller
)

:: ---- 2. Clean previous build ----
if exist dist\MMexplorer rmdir /s /q dist\MMexplorer
if exist build rmdir /s /q build

:: ---- 3. PyInstaller bundle ----
echo.
echo [1/2] Building executable with PyInstaller...
pyinstaller protein_mutation.spec --noconfirm
if errorlevel 1 (
    echo ERROR: PyInstaller failed.
    exit /b 1
)

:: ---- 4. Inno Setup ----
echo.
echo [2/2] Creating installer with Inno Setup...

set ISCC="%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe"
if not exist %ISCC% set ISCC="%ProgramFiles%\Inno Setup 6\ISCC.exe"

if not exist %ISCC% (
    echo ERROR: Inno Setup 6 not found. Download from https://jrsoftware.org/isdl.php
    exit /b 1
)

%ISCC% installer\installer.iss
if errorlevel 1 (
    echo ERROR: Inno Setup compilation failed.
    exit /b 1
)

echo.
echo ============================================================
echo  Done! Installer: installer\output\MMexplorer-1.0.0-Setup.exe
echo ============================================================
