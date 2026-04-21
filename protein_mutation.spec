# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for Protein Missense Mutation Explorer
# Copyright (c) 2026 fredsanto (Federico Santoni)

import sys
from pathlib import Path
from PyInstaller.utils.hooks import collect_data_files, collect_submodules

block_cipher = None

# Collect OpenMM XML force-field files
openmm_datas = collect_data_files("openmm")
pdbfixer_datas = collect_data_files("pdbfixer")

a = Analysis(
    ["app.py"],
    pathex=["."],
    binaries=[],
    datas=[
        ("templates", "templates"),
        *openmm_datas,
        *pdbfixer_datas,
    ],
    hiddenimports=[
        "openmm",
        "openmm.app",
        "openmm.unit",
        "pdbfixer",
        "Bio",
        "Bio.PDB",
        "Bio.PDB.Polypeptide",
        "Bio.SeqUtils",
        "flask",
        "werkzeug",
        "jinja2",
        "numpy",
        "requests",
        *collect_submodules("openmm"),
        *collect_submodules("pdbfixer"),
    ],
    hookspath=[],
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="MMexplorer",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,           # no black console window
    icon="installer\\icon.ico",
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    name="MMexplorer",
)
