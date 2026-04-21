"""
Standalone debug script for the 'object of type Atom has no len()' bug.

Run normally:
    python3.10 debug_report.py

Run with pdb:
    python3.10 -m pdb debug_report.py

The script loads the saved 4E26 PDB, applies a mutation with PDBFixer
(no full minimization so it finishes fast), then calls _analyze_biochemistry
on both the original and mutant structures — the same code path as /api/report.
"""

import io, sys, traceback, tempfile, os
sys.path.insert(0, os.path.dirname(__file__))

# ── imports from app ──────────────────────────────────────────────────────────
from app import (
    _analyze_biochemistry, _diff_biochemistry,
    _protein_only_pdb, _clean_pdb, _normalize_resnames,
    _highlight_residue, _locate_residue, _restore_numbering,
    AA1, AA3,
)
from pdbfixer import PDBFixer
from openmm import app as omm_app
import openmm as mm
from openmm import unit
from Bio.PDB import PDBParser

# ── config ────────────────────────────────────────────────────────────────────
PDB_FILE   = "/tmp/tmp9hg_ww6g.pdb"   # most recent saved 4E26 structure
MUTATION   = "E451K"                   # change to whatever mutation was tested
MINIMIZE   = True                      # set False for a ~30-second run; True for full

# ── helpers (same logic as _fix_and_minimise in app.py) ──────────────────────
def parse_mutation(s):
    import re
    s = s.strip().lstrip("p.")
    m3 = re.match(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", s)
    if m3:
        from app import AA1 as _AA1
        ref = _AA1.get(m3.group(1).upper(), m3.group(1)[0].upper())
        alt = _AA1.get(m3.group(3).upper(), m3.group(3)[0].upper())
        return ref, int(m3.group(2)), alt
    m1 = re.match(r"([A-Z])(\d+)([A-Z])", s.upper())
    if m1:
        return m1.group(1), int(m1.group(2)), m1.group(3)
    raise ValueError(f"Cannot parse mutation {s!r}")

def run_pipeline(raw_pdb, ref_aa, resnum, alt_aa, minimize=True):
    clean = _protein_only_pdb(raw_pdb)

    # Snapshot original numbering
    parser = PDBParser(QUIET=True)
    orig_seq = {}
    for ch in parser.get_structure("_o", io.StringIO(clean))[0].get_chains():
        rl = [(r.get_id()[1], r.get_resname().strip())
              for r in ch.get_residues()
              if r.get_id()[0] == " " and r.get_resname().strip() in AA1]
        if rl:
            orig_seq[ch.id] = rl

    actual_chain, actual_resnum, actual_resname = _locate_residue(clean, ref_aa, resnum)
    print(f"Mutation site: {actual_resname}{actual_resnum} chain {actual_chain}")

    alt3 = AA3.get(alt_aa.upper(), alt_aa.upper())
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w") as f:
        f.write(clean)
        tmp = f.name

    try:
        print("Running PDBFixer …")
        fixer = PDBFixer(filename=tmp)
        fixer.applyMutations([f"{actual_resname}-{actual_resnum}-{alt3}"], actual_chain)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        print(f"  atoms after fixer: {fixer.topology.getNumAtoms()}")

        if minimize:
            print("Building force field …")
            ff = omm_app.ForceField("amber14-all.xml", "implicit/obc2.xml")
            modeller = omm_app.Modeller(fixer.topology, fixer.positions)
            system = ff.createSystem(
                modeller.topology,
                nonbondedMethod=omm_app.NoCutoff,
                constraints=omm_app.HBonds,
            )
            restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
            restraint.addGlobalParameter("k", 1000.0 * unit.kilojoules_per_mole / unit.nanometer**2)
            restraint.addPerParticleParameter("x0")
            restraint.addPerParticleParameter("y0")
            restraint.addPerParticleParameter("z0")
            pos_nm = modeller.positions.value_in_unit(unit.nanometers)
            BACKBONE = {"CA", "C", "N", "O"}
            for atom in modeller.topology.atoms():
                if atom.name in BACKBONE:
                    x0, y0, z0 = pos_nm[atom.index]
                    restraint.addParticle(atom.index, [x0, y0, z0])
            system.addForce(restraint)

            integrator = mm.LangevinMiddleIntegrator(
                300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds
            )
            sim = omm_app.Simulation(modeller.topology, system, integrator)
            sim.context.setPositions(modeller.positions)
            print("Minimising (500 steps) …")
            for i in range(10):
                sim.minimizeEnergy(maxIterations=50)
                e = sim.context.getState(getEnergy=True).getPotentialEnergy()
                print(f"  step {(i+1)*50}: {e}")
            positions = sim.context.getState(getPositions=True).getPositions()
            buf = io.StringIO()
            omm_app.PDBFile.writeFile(sim.topology, positions, buf)
        else:
            # Skip minimisation — use PDBFixer positions directly
            print("Skipping minimisation (MINIMIZE=False)")
            buf = io.StringIO()
            omm_app.PDBFile.writeFile(fixer.topology, fixer.positions, buf)

        mutant_raw = buf.getvalue()
    finally:
        os.unlink(tmp)

    renumbered  = _restore_numbering(mutant_raw, orig_seq)
    mutant_pdb  = _normalize_resnames(_clean_pdb(renumbered))
    original_pdb = _clean_pdb(raw_pdb)

    mutant_pdb   = _highlight_residue(mutant_pdb,   actual_chain, actual_resnum)
    original_pdb = _highlight_residue(original_pdb, actual_chain, actual_resnum)

    return original_pdb, mutant_pdb, actual_chain, actual_resnum

# ── main ──────────────────────────────────────────────────────────────────────
def main():
    print(f"Loading {PDB_FILE}")
    with open(PDB_FILE) as f:
        raw_pdb = f.read()

    ref_aa, resnum, alt_aa = parse_mutation(MUTATION)
    print(f"Mutation: {ref_aa}{resnum}{alt_aa}")

    original_pdb, mutant_pdb, chain, actual_resnum = run_pipeline(
        raw_pdb, ref_aa, resnum, alt_aa, minimize=MINIMIZE
    )

    print(f"\nAnalysing original PDB …")
    orig_a = _analyze_biochemistry(original_pdb, chain, actual_resnum)
    print(f"  ok — hbonds={len(orig_a.get('hbonds', []))}, "
          f"charge_interactions={len(orig_a.get('charge_interactions', []))}")

    print(f"Analysing mutant PDB …")
    mut_a = _analyze_biochemistry(mutant_pdb, chain, actual_resnum)
    print(f"  ok — hbonds={len(mut_a.get('hbonds', []))}, "
          f"charge_interactions={len(mut_a.get('charge_interactions', []))}")

    print(f"Computing diff …")
    diff = _diff_biochemistry(orig_a, mut_a)
    print(f"  charge_delta={diff['charge_delta']}, "
          f"hbonds_delta={diff['hbonds_delta']}")

    print("\nAll steps completed without error.")

if __name__ == "__main__":
    main()
