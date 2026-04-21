# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec — macOS
# Copyright (c) 2026 fredsanto (Federico Santoni)

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

block_cipher = None

openmm_datas   = collect_data_files("openmm")
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
        "openmm", "openmm.app", "openmm.unit",
        "pdbfixer",
        "Bio", "Bio.PDB", "Bio.PDB.Polypeptide", "Bio.SeqUtils",
        "flask", "werkzeug", "jinja2",
        "numpy", "requests",
        *collect_submodules("openmm"),
        *collect_submodules("pdbfixer"),
    ],
    hookspath=[],
    runtime_hooks=[],
    excludes=[],
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
    strip=False,
    upx=False,          # UPX unreliable on macOS ARM
    console=False,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    name="MMexplorer",
)

app = BUNDLE(
    coll,
    name="MMexplorer.app",
    icon=None,
    bundle_identifier="com.fredsanto.mmexplorer",
    info_plist={
        "CFBundleName":               "MMexplorer",
        "CFBundleDisplayName":        "MMexplorer",
        "CFBundleVersion":            "1.0.0",
        "CFBundleShortVersionString": "1.0.0",
        "CFBundleIdentifier":         "com.fredsanto.mmexplorer",
        "NSHighResolutionCapable":    True,
        "LSMinimumSystemVersion":     "12.0",
        "NSHumanReadableCopyright":   "Copyright © 2026 fredsanto (Federico Santoni)",
    },
)
