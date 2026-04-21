#!/usr/bin/env bash
# Protein Missense Mutation Explorer — macOS build script
# Copyright (c) 2026 fredsanto (Federico Santoni)
# Requires: conda env with openmm/pdbfixer/biopython, PyInstaller, create-dmg

set -euo pipefail

APP_NAME="MMexplorer"
APP_BUNDLE="MMexplorer.app"
DMG_NAME="MMexplorer-1.0.0-mac"
VOLUME_NAME="$APP_NAME"

echo "============================================================"
echo "  $APP_NAME — macOS build"
echo "============================================================"

# ---- 1. PyInstaller ----
echo ""
echo "[1/3] Installing PyInstaller if needed..."
pip show pyinstaller &>/dev/null || pip install pyinstaller

echo "[1/3] Bundling app with PyInstaller..."
rm -rf "dist/$APP_BUNDLE" build/MMexplorer
pyinstaller protein_mutation_mac.spec --noconfirm

# ---- 2. Ad-hoc code sign (allows Gatekeeper to run without a paid cert) ----
echo ""
echo "[2/3] Code signing (ad-hoc)..."
codesign --force --deep --sign - "dist/$APP_BUNDLE"

# ---- 3. DMG ----
echo ""
echo "[3/3] Creating DMG..."

if ! command -v create-dmg &>/dev/null; then
    echo "  create-dmg not found — installing via Homebrew..."
    brew install create-dmg
fi

rm -f "dist/${DMG_NAME}.dmg"

create-dmg \
    --volname "$VOLUME_NAME" \
    --window-pos 200 120 \
    --window-size 600 400 \
    --icon-size 128 \
    --icon "$APP_BUNDLE" 150 185 \
    --hide-extension "$APP_BUNDLE" \
    --app-drop-link 450 185 \
    "dist/${DMG_NAME}.dmg" \
    "dist/$APP_BUNDLE" \
|| {
    # Fallback: plain DMG without background/positioning if assets missing
    hdiutil create \
        -volname "$VOLUME_NAME" \
        -srcfolder "dist/$APP_BUNDLE" \
        -ov -format UDZO \
        "dist/${DMG_NAME}.dmg"
}

echo ""
echo "============================================================"
echo "  Done!  dist/${DMG_NAME}.dmg"
echo "============================================================"
