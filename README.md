# MMexplorer

**Fetch a protein structure, apply a missense mutation, energy-minimise with OpenMM, and compare original vs mutant side-by-side in the browser.**

> By [@fredsanto](https://github.com/fredsanto) · MIT License

---

## Features

- **Auto structure selection** — searches RCSB PDB for crystal structures (sorted by resolution); falls back to AlphaFold v4+ if none found
- **Missense mutation** — applies the substitution via PDBFixer (`Arg175His`, `R175H`, etc.)
- **Energy minimisation** — AMBER14 force field + OBC2 implicit solvent; backbone restraints keep the global fold intact while side chains relax
- **NVT dynamics** — 8 ps short MD for per-residue Cα RMSF
- **3D viewer** — 3Dmol.js; cartoon / sticks / surface / spheres; auto-highlights mutation site
- **Biochemistry report** — Claude AI summarises structural impact, conservation, and clinical relevance
- **UniProt annotations** — active sites, binding sites, PTMs, disease variants at the mutation position
- **Session save/load** — results persisted locally; reload any previous run
- **Stop button** — cancel any in-progress calculation including energy minimisation

---

## Quick Start

### Requirements

- Python ≥ 3.10
- [conda](https://docs.conda.io/) recommended (OpenMM easiest via conda-forge)

### Install

```bash
conda create -n mmexplorer -c conda-forge python=3.11 openmm pdbfixer
conda activate mmexplorer
pip install flask biopython requests numpy
```

### Run

```bash
python app.py
# Opens http://localhost:5050
```

---

## Usage

1. Enter a **protein name or gene symbol** (e.g. `TP53`, `BRCA1`, `p53`)
2. Click **Browse structures** to pick a PDB entry or AlphaFold model
3. Enter a **missense mutation** in any format: `R175H`, `Arg175His`, `R175H`
4. Click **Analyse** — progress streams live
5. Inspect 3D viewers, RMSF plot, biochemistry report, and UniProt annotations
6. Use **Stop** at any time to cancel

---

## Windows Installer

Build a standalone `.exe` installer on Windows:

```bat
conda activate mmexplorer
pip install pyinstaller
build_windows.bat
```

Requires [Inno Setup 6](https://jrsoftware.org/isdl.php). Output: `installer\output\MMexplorer-1.0.0-Setup.exe`

---

## Project Structure

```
MMexplorer/
├── app.py                  # Flask backend — all routes + pipeline
├── templates/
│   └── index.html          # Single-page frontend (3Dmol.js)
├── sessions/               # Saved session JSON files
├── protein_mutation.spec   # PyInstaller build spec
├── installer/
│   ├── installer.iss       # Inno Setup script
│   └── license.txt         # MIT license
├── build_windows.bat       # One-click Windows build
└── pyproject.toml          # Package metadata
```

---

## Dependencies

| Package | Purpose |
|---------|---------|
| Flask | Web server |
| OpenMM | MD / energy minimisation |
| PDBFixer | Mutation application, missing atom repair |
| Biopython | PDB parsing, residue lookup |
| Requests | RCSB / UniProt / AlphaFold API calls |
| 3Dmol.js | In-browser 3D molecular viewer (CDN) |

---

## License

MIT © 2026 [fredsanto](https://github.com/fredsanto)
