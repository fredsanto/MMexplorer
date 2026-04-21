"""
Protein Missense Mutation Explorer
-----------------------------------
1. Fetch best available PDB structure (crystal > AlphaFold)
2. Apply the requested missense mutation with PDBFixer
3. Energy-minimise the mutant with OpenMM / AMBER14
4. Return original + minimised-mutant PDB strings to the browser

Copyright (c) 2026 fredsanto (Federico Santoni)
SPDX-License-Identifier: MIT
"""

__author__    = "fredsanto (Federico Santoni)"
__copyright__ = "Copyright (c) 2026 fredsanto"
__license__   = "MIT"
__version__   = "1.0.0"

import io, os, json, tempfile, re, requests, gzip, queue, threading, concurrent.futures, math
from datetime import datetime, timezone
import numpy as np
from flask import Flask, request, jsonify, render_template, Response, stream_with_context

# Biopython — structure parsing / chain/residue inspection
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Polypeptide import PPBuilder
from Bio.SeqUtils import seq1

# OpenMM ecosystem
from pdbfixer import PDBFixer
import openmm as mm
from openmm import app as omm_app, unit

app = Flask(__name__)

SESSIONS_DIR = os.path.join(os.path.dirname(__file__), "sessions")
os.makedirs(SESSIONS_DIR, exist_ok=True)

_stop_event = threading.Event()


class _Stopped(Exception):
    pass

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RCSB_SEARCH   = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA     = "https://data.rcsb.org/rest/v1/core/entry/{}"
RCSB_DL       = "https://files.rcsb.org/download/{}.pdb"
UNIPROT_SRCH  = "https://rest.uniprot.org/uniprotkb/search?query={}&format=json&fields=accession,gene_names&size=1"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{}.json"
PDBE_SIFTS    = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{}"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{}"
ALPHAFOLD_DL  = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb"

# UniProt feature types to include (exact strings returned by rest.uniprot.org)
_FEAT_KEEP = {
    "Domain", "Region", "Motif", "Coiled coil", "Repeat",
    "Zinc finger", "DNA binding",
    "Active site", "Binding site", "Site",
    "Modified residue", "Lipidation", "Glycosylation", "Cross-link",
    "Disulfide bond",
    "Signal peptide", "Propeptide", "Transit peptide",
    "Transmembrane", "Topological domain", "Intramembrane",
    "Natural variant", "Mutagenesis",
}

AA3 = {
    "A":"ALA","C":"CYS","D":"ASP","E":"GLU","F":"PHE",
    "G":"GLY","H":"HIS","I":"ILE","K":"LYS","L":"LEU",
    "M":"MET","N":"ASN","P":"PRO","Q":"GLN","R":"ARG",
    "S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR",
}
AA1 = {v: k for k, v in AA3.items()}

# ---------------------------------------------------------------------------
# Helpers — structure fetching
# ---------------------------------------------------------------------------

def _best_pdb(protein_query: str) -> tuple[str, str, str]:
    """
    Returns (pdb_id, chain_id, source) where source is 'crystal' or 'alphafold'.
    pdb_id may be a UniProt accession when source == 'alphafold'.
    """
    # --- try RCSB text search ---
    q = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": protein_query}
        },
        "request_options": {
            "paginate": {"start": 0, "rows": 20},
            "sort": [{"sort_by": "score", "direction": "desc"}]
        },
        "return_type": "entry"
    }
    r = requests.post(RCSB_SEARCH, json=q, timeout=10)
    if r.ok:
        hits = r.json().get("result_set", [])
        # Filter for X-ray crystal structures, pick highest resolution
        best_pdb, best_res, best_chain = None, 999, "A"
        for h in hits:
            pid = h["identifier"]
            meta = requests.get(RCSB_DATA.format(pid), timeout=5)
            if not meta.ok:
                continue
            info = meta.json()
            method = info.get("exptl", [{}])[0].get("method", "")
            if method.upper() not in ("X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"):
                continue
            res = info.get("rcsb_entry_info", {}).get("resolution_combined", [999])[0]
            try:
                res = float(res)
            except (TypeError, ValueError):
                res = 999
            chains = info.get("rcsb_entry_info", {}).get("polymer_entity_count_protein", 0)
            if res < best_res and chains > 0:
                best_res = res
                best_pdb = pid
                # pick first protein chain
                poly = requests.get(
                    f"https://data.rcsb.org/rest/v1/core/entry/{pid}/polymer_entities",
                    timeout=5
                )
                if poly.ok:
                    for ent in poly.json():
                        cids = ent.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_chain_ids", [])
                        if cids:
                            best_chain = cids[0]; break
        if best_pdb:
            return best_pdb, best_chain, f"crystal ({best_res:.1f} Å)"

    # --- fallback: AlphaFold via UniProt ---
    uni_r = requests.get(UNIPROT_SRCH.format(requests.utils.quote(protein_query)), timeout=8)
    if uni_r.ok:
        results = uni_r.json().get("results", [])
        if results:
            acc = results[0]["primaryAccession"]
            af_check = requests.get(ALPHAFOLD_API.format(acc), timeout=8)
            if af_check.ok and af_check.json():
                return acc, "A", "AlphaFold v4"

    raise ValueError(f"No structure found for '{protein_query}'")


def _fetch_entry_info(pid: str) -> dict | None:
    """Fetch metadata for one RCSB entry. Returns None on failure or non-protein."""
    METHOD_LABEL = {
        "X-RAY DIFFRACTION":  "X-Ray",
        "ELECTRON MICROSCOPY":"Cryo-EM",
        "SOLUTION NMR":       "NMR",
        "NEUTRON DIFFRACTION":"Neutron",
    }
    try:
        meta = requests.get(RCSB_DATA.format(pid), timeout=6)
        if not meta.ok:
            return None
        info   = meta.json()
        method = info.get("exptl", [{}])[0].get("method", "").upper()
        if method not in METHOD_LABEL:
            return None
        n_prot = info.get("rcsb_entry_info", {}).get("polymer_entity_count_protein", 0)
        if not n_prot:
            return None
        res = info.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
        try:
            res = round(float(res), 2) if res is not None else None
        except (TypeError, ValueError):
            res = None
        title = info.get("struct", {}).get("title", pid)[:100]
        # first protein chain
        chain = "A"
        poly = requests.get(
            f"https://data.rcsb.org/rest/v1/core/entry/{pid}/polymer_entities",
            timeout=6
        )
        if poly.ok:
            for ent in poly.json():
                cids = (ent.get("rcsb_polymer_entity_container_identifiers", {})
                           .get("auth_chain_ids", []))
                if cids:
                    chain = cids[0]; break
        return {
            "id":         pid,
            "chain":      chain,
            "source":     "crystal",
            "method":     METHOD_LABEL[method],
            "resolution": res,
            "title":      title,
        }
    except Exception:
        return None


def _list_structures(protein_query: str, max_crystal: int = 10) -> list[dict]:
    """
    Return up to `max_crystal` RCSB structures (sorted by resolution) plus
    an AlphaFold prediction if one exists, for the given protein query.
    RCSB metadata is fetched in parallel to keep latency low.
    """
    pdb_ids: list[str] = []
    q = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": protein_query},
        },
        "request_options": {
            "paginate": {"start": 0, "rows": 40},
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
        "return_type": "entry",
    }
    try:
        r = requests.post(RCSB_SEARCH, json=q, timeout=10)
        if r.ok:
            pdb_ids = [h["identifier"] for h in r.json().get("result_set", [])[:40]]
    except Exception:
        pass

    structures: list[dict] = []
    if pdb_ids:
        with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
            for info in ex.map(_fetch_entry_info, pdb_ids):
                if info:
                    structures.append(info)

    # Sort by resolution (NMR / no-res go last among crystal entries)
    structures.sort(key=lambda s: (s["resolution"] is None, s["resolution"] or 0))
    structures = structures[:max_crystal]

    # AlphaFold
    try:
        uni_r = requests.get(
            UNIPROT_SRCH.format(requests.utils.quote(protein_query)), timeout=8
        )
        if uni_r.ok:
            hits = uni_r.json().get("results", [])
            if hits:
                acc  = hits[0]["primaryAccession"]
                gene = (hits[0].get("genes") or [{}])[0].get("geneName", {}).get("value", "")
                prot = (hits[0].get("proteinDescription", {})
                              .get("recommendedName", {})
                              .get("fullName", {})
                              .get("value", ""))
                af = requests.get(ALPHAFOLD_API.format(acc), timeout=8)
                if af.ok and af.json():
                    structures.append({
                        "id":         acc,
                        "chain":      "A",
                        "source":     "alphafold",
                        "method":     "AlphaFold v4",
                        "resolution": None,
                        "title":      prot or gene or acc,
                    })
    except Exception:
        pass

    return structures


def _fetch_pdb_string(pdb_id: str, source: str) -> str:
    """Download PDB text for the given id."""
    if "alphafold" in source.lower():
        api = requests.get(ALPHAFOLD_API.format(pdb_id), timeout=10)
        api.raise_for_status()
        entries = api.json()
        if not entries:
            raise ValueError(f"No AlphaFold entry for {pdb_id}")
        url = entries[0]["pdbUrl"]
    else:
        url = RCSB_DL.format(pdb_id.upper())
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.text


def _locate_residue(pdb_string: str, ref_aa: str, resnum: int) -> tuple[str, int, str]:
    """
    Scan every protein chain in the PDB and return (chain_id, pdb_resnum, resname3)
    for the residue that best matches the user's request.

    Strategy (in order):
      1. Exact author-sequence-number match + ref_aa match in any chain
      2. Exact number match ignoring ref_aa (structure may have a different AA)
      3. Sequential (1-indexed) position within a chain that matches ref_aa
      4. Closest residue matching ref_aa across all chains
    Raises ValueError if nothing is found.
    """
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("X", io.StringIO(pdb_string))
    model  = struct[0]
    ref3   = AA3.get(ref_aa.upper())

    exact_aa, exact_num, seq_match = None, None, None
    closest, closest_dist = None, 10**9

    for chain in model.get_chains():
        protein_res = [
            r for r in chain.get_residues()
            if r.get_id()[0] == " " and r.get_resname().strip() in AA1
        ]
        for idx, res in enumerate(protein_res):
            rnum  = res.get_id()[1]
            rname = res.get_resname().strip()
            cid   = chain.id

            # Strategy 1
            if rnum == resnum and (ref3 is None or rname == ref3):
                exact_aa = (cid, rnum, rname)

            # Strategy 2
            if rnum == resnum and exact_num is None:
                exact_num = (cid, rnum, rname)

            # Strategy 3: sequential position (1-indexed)
            if idx + 1 == resnum and (ref3 is None or rname == ref3) and seq_match is None:
                seq_match = (cid, rnum, rname)

            # Strategy 4: closest by author-number distance among matching AA
            if ref3 and rname == ref3:
                d = abs(rnum - resnum)
                if d < closest_dist:
                    closest_dist = d
                    closest = (cid, rnum, rname)

    for candidate in (exact_aa, exact_num, seq_match, closest):
        if candidate is not None:
            return candidate

    raise ValueError(
        f"Residue {ref_aa}{resnum} not found in any chain. "
        "Check that the residue number corresponds to PDB author numbering "
        "or to the sequential position within the chain."
    )


def _protein_only_pdb(pdb_string: str) -> str:
    """Keep only ATOM / TER / END records (drop HETATM, waters, ligands)."""
    lines = [
        ln for ln in pdb_string.splitlines()
        if ln[:6].rstrip() in ("ATOM", "TER", "END", "HEADER", "REMARK")
    ]
    return "\n".join(lines)


def _parse_residue_number(mut_str: str) -> tuple[str, int, str]:
    """
    Parse mutation string such as 'A123G', 'Ala123Gly', 'p.Ala123Gly'.
    Returns (ref_aa_1letter, residue_number, alt_aa_1letter).
    """
    mut_str = mut_str.strip().lstrip("p.")
    # Three-letter form: Ala123Gly
    m3 = re.match(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", mut_str)
    if m3:
        ref = AA1.get(m3.group(1).upper(), m3.group(1)[0].upper())
        alt = AA1.get(m3.group(3).upper(), m3.group(3)[0].upper())
        return ref, int(m3.group(2)), alt
    # One-letter form: A123G
    m1 = re.match(r"([A-Z])(\d+)([A-Z])", mut_str.upper())
    if m1:
        return m1.group(1), int(m1.group(2)), int(m1.group(3)) if False else m1.group(3)
    raise ValueError(f"Cannot parse mutation '{mut_str}'. Use format A123G or Ala123Gly.")


# ---------------------------------------------------------------------------
# Core pipeline
# ---------------------------------------------------------------------------

def _restore_numbering(mutant_pdb: str, orig_seq: dict) -> str:
    """
    OpenMM's PDBFile.writeFile resets all residue numbers to 1-based sequential.
    This function maps them back to the original PDB author numbering by matching
    residues positionally (same sequence order) within each chain.

    orig_seq: {chain_id: [(orig_resnum, resname), …]}  captured before PDBFixer ran.
    """
    parser = PDBParser(QUIET=True)
    mut_model = parser.get_structure("_m", io.StringIO(mutant_pdb))[0]

    mut_seq: dict[str, list] = {}
    for ch in mut_model.get_chains():
        rl = [(r.get_id()[1], r.get_resname().strip())
              for r in ch.get_residues() if r.get_id()[0] == " "]
        if rl:
            mut_seq[ch.id] = rl

    orig_chain_ids = list(orig_seq.keys())
    mut_chain_ids  = list(mut_seq.keys())

    # remap[(mut_chain, mut_resnum)] = (orig_chain, orig_resnum)
    remap: dict[tuple, tuple] = {}
    for i, mc in enumerate(mut_chain_ids):
        if i >= len(orig_chain_ids):
            break
        oc      = orig_chain_ids[i]
        o_list  = orig_seq[oc]
        m_list  = mut_seq[mc]
        for j in range(min(len(o_list), len(m_list))):
            remap[(mc, m_list[j][0])] = (oc, o_list[j][0])

    fixed = []
    for line in mutant_pdb.splitlines():
        if line[:4] in ("ATOM", "HETA") and len(line) >= 26:
            try:
                mc, mr = line[21], int(line[22:26])
                if (mc, mr) in remap:
                    oc, ornum = remap[(mc, mr)]
                    line = line[:21] + oc + f"{ornum:4d}" + line[26:]
            except ValueError:
                pass
        elif line.startswith("TER") and len(line) >= 26:
            try:
                mc, mr = line[21], int(line[22:26])
                if (mc, mr) in remap:
                    oc, ornum = remap[(mc, mr)]
                    line = line[:21] + oc + f"{ornum:4d}" + line[26:]
            except (ValueError, IndexError):
                pass
        fixed.append(line)

    return "\n".join(fixed)


def _fix_and_minimise(pdb_string: str, ref_aa: str, resnum: int, alt_aa: str,
                      on_progress=None, stop_event: threading.Event | None = None) -> tuple[str, str, int]:
    """
    Locate the residue in the actual structure, apply the mutation, minimise.
    Returns (minimised_pdb_string, actual_chain, actual_pdb_resnum).
    on_progress(msg): optional callable invoked at each major step.
    stop_event: optional threading.Event; when set, raises _Stopped immediately.
    """
    def _check_stop():
        if stop_event and stop_event.is_set():
            raise _Stopped()

    def _prog(msg):
        if on_progress:
            on_progress(msg)

    # Strip HETATM / waters so AMBER force field doesn't choke on unknown ligands
    clean = _protein_only_pdb(pdb_string)

    # Snapshot original residue numbering BEFORE PDBFixer touches anything.
    # OpenMM's PDBFile.writeFile resets to 1-based sequential; we restore from this.
    _parser = PDBParser(QUIET=True)
    _orig_seq: dict[str, list] = {}   # chain_id → [(orig_resnum, resname), …]
    for _ch in _parser.get_structure("_o", io.StringIO(clean))[0].get_chains():
        _rl = [(r.get_id()[1], r.get_resname().strip())
               for r in _ch.get_residues()
               if r.get_id()[0] == " " and r.get_resname().strip() in AA1]
        if _rl:
            _orig_seq[_ch.id] = _rl

    _prog("Locating residue in structure…")
    actual_chain, actual_resnum, actual_resname = _locate_residue(clean, ref_aa, resnum)

    _prog(f"Applying mutation {actual_resname}{actual_resnum}→{AA3.get(alt_aa.upper(), alt_aa)} on chain {actual_chain}…")

    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w") as f:
        f.write(clean)
        tmp = f.name
    try:
        fixer = PDBFixer(filename=tmp)

        # PDBFixer mutation format: "EXISTING_RESNAME-RESNUM-NEW_RESNAME"
        # chain_id is the second argument to applyMutations (not part of the string)
        alt3 = AA3.get(alt_aa.upper(), alt_aa.upper())
        fixer.applyMutations([f"{actual_resname}-{actual_resnum}-{alt3}"], actual_chain)

        _prog("Repairing missing residues and heavy atoms…")
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        _prog("Adding hydrogens (pH 7.0)…")
        fixer.addMissingHydrogens(7.0)

        # OBC2 is ~5x faster than GBn2 and accurate enough for local side-chain moves
        _prog("Building AMBER14 force field (OBC2 implicit solvent)…")
        ff = omm_app.ForceField("amber14-all.xml", "implicit/obc2.xml")
        modeller = omm_app.Modeller(fixer.topology, fixer.positions)
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=omm_app.NoCutoff,
            constraints=omm_app.HBonds,
        )

        # Backbone positional restraints: fix CA/C/N/O so only side chains relax.
        # Reduces effective DOF by ~75% and prevents global drift.
        _prog("Applying backbone restraints…")
        restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint.addGlobalParameter(
            "k", 1000.0 * unit.kilojoules_per_mole / unit.nanometer**2
        )
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        positions_nm = modeller.positions.value_in_unit(unit.nanometers)
        BACKBONE = {"CA", "C", "N", "O"}
        for atom in modeller.topology.atoms():
            if atom.name in BACKBONE:
                x0, y0, z0 = positions_nm[atom.index]
                restraint.addParticle(atom.index, [x0, y0, z0])
        system.addForce(restraint)

        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds
        )
        simulation = omm_app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        # Chunked minimisation with live energy readout
        _prog("Energy minimisation starting…")
        chunk, n_chunks = 50, 10          # 500 iterations; enough with backbone restraints
        prev_e = None
        for i in range(n_chunks):
            _check_stop()
            simulation.minimizeEnergy(maxIterations=chunk)
            state = simulation.context.getState(getEnergy=True)
            e = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            delta = f"  Δ={e - prev_e:+.0f}" if prev_e is not None else ""
            _prog(f"Minimising… step {(i+1)*chunk}/{n_chunks*chunk} "
                  f"· E = {e:.0f} kJ/mol{delta}")
            prev_e = e

        # ── Short NVT for per-residue RMSF ───────────────────────────────────
        _prog("Running short NVT dynamics for Cα RMSF (8 ps)…")
        simulation.context.setParameter("k", 0)   # release backbone restraints
        simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

        # Index Cα atoms once
        ca_idx = {}   # (chain_id_str, res_id_str) → atom index
        for atom in simulation.topology.atoms():
            if atom.name == "CA":
                ca_idx[(atom.residue.chain.id, atom.residue.id)] = atom.index

        n_steps_per_frame, n_frames = 50, 40   # 40 × 50 × 4 fs = 8 ps
        ca_traj = {k: [] for k in ca_idx}

        for _ in range(n_frames):
            _check_stop()
            simulation.step(n_steps_per_frame)
            state = simulation.context.getState(getPositions=True)
            pos_A = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)
            for key, idx in ca_idx.items():
                ca_traj[key].append(pos_A[idx].copy())

        rmsf_dict = {}
        for (ch, rn_str), frames in ca_traj.items():
            if len(frames) < 2:
                continue
            arr  = np.array(frames)
            rmsf = float(np.sqrt(np.mean(np.sum((arr - arr.mean(axis=0))**2, axis=1))))
            rmsf_dict[f"{ch}:{rn_str}"] = round(rmsf, 3)

        _prog("Writing minimised structure…")
        positions = simulation.context.getState(getPositions=True).getPositions()
        buf = io.StringIO()
        omm_app.PDBFile.writeFile(simulation.topology, positions, buf)

        _prog("Restoring original residue numbering…")
        renumbered = _restore_numbering(buf.getvalue(), _orig_seq)
        return renumbered, actual_chain, actual_resnum, rmsf_dict
    finally:
        os.unlink(tmp)


def _normalize_resnames(pdb_string: str) -> str:
    """
    Map AMBER14/OpenMM variant residue names back to canonical 3-letter codes.

    AMBER14 assigns protonation-state-specific names (HIE/HID/HIP, LYN, ASH,
    GLH, CYM, CYX) and terminal prefixes (NALA, CALA, …) that end up in the
    OpenMM PDB output and look wrong in the viewer.  We normalise them here so
    mutant labels stay consistent with the original PDB.
    """
    VARIANTS = {
        "HIE": "HIS", "HID": "HIS", "HIP": "HIS",
        "LYN": "LYS",
        "ASH": "ASP",
        "GLH": "GLU",
        "CYM": "CYS", "CYX": "CYS",
    }
    AA3_SET = set(AA3.values())
    lines = []
    for line in pdb_string.splitlines():
        if line[:4] in ("ATOM", "HETA"):
            # PDB residue name field: columns 17-20 (0-indexed)
            resn = line[17:21].strip()
            canonical = None
            if resn in VARIANTS:
                canonical = VARIANTS[resn]
            # N-terminal AMBER names: NALA, NGLY, NVAL, …
            elif len(resn) == 4 and resn[0] == "N" and resn[1:] in AA3_SET:
                canonical = resn[1:]
            # C-terminal AMBER names: CALA, CGLY, CVAL, …
            elif len(resn) == 4 and resn[0] == "C" and resn[1:] in AA3_SET:
                canonical = resn[1:]
            if canonical:
                # Replace columns 17-20 with left-justified canonical name (+ space)
                line = line[:17] + canonical.ljust(4) + line[21:]
        lines.append(line)
    return "\n".join(lines)


def _clean_pdb(pdb_string: str) -> str:
    """Remove HETATM / water / non-standard lines for clean display."""
    lines = []
    for line in pdb_string.splitlines():
        rec = line[:6].strip()
        if rec in ("ATOM", "TER", "END", "HEADER", "REMARK"):
            lines.append(line)
    return "\n".join(lines)


def _highlight_residue(pdb_string: str, chain: str, resnum: int) -> str:
    """Insert a REMARK line noting the mutation site (used by frontend)."""
    return f"REMARK 999 MUTATION_SITE {chain} {resnum}\n" + pdb_string


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Functional site annotation (UniProt via PDBe SIFTS)
# ---------------------------------------------------------------------------

def _sifts_uniprot(pdb_id: str, chain_id: str, resnum: int):
    """
    Map a PDB chain residue to a UniProt accession + residue number.
    Returns (accession, gene_id, unp_resnum) or (None, None, None).

    For AlphaFold models (pdb_id like "AF-P04637-F1") the accession is
    embedded in the ID and residue numbering equals the UniProt position.
    """
    # AlphaFold — accession is in the ID, numbering is already UniProt-based
    m = re.match(r'^AF-([A-Z0-9]+)-F\d+$', pdb_id, re.IGNORECASE)
    if m:
        acc = m.group(1).upper()
        return acc, acc, resnum

    try:
        resp = requests.get(PDBE_SIFTS.format(pdb_id.lower()), timeout=8)
        if resp.status_code != 200:
            return None, None, None
        uniprot_block = resp.json().get(pdb_id.lower(), {}).get("UniProt", {})
        for acc, info in uniprot_block.items():
            gene_id = info.get("identifier", acc)
            for seg in info.get("mappings", []):
                if seg.get("chain_id") != chain_id and seg.get("struct_asym_id") != chain_id:
                    continue
                start_d = seg.get("start") or {}
                end_d   = seg.get("end")   or {}
                # Prefer author numbering; fall back to sequential residue_number
                pdb_s = start_d.get("author_residue_number") or start_d.get("residue_number")
                pdb_e = end_d.get("author_residue_number")   or end_d.get("residue_number")
                unp_s = seg.get("unp_start")
                if pdb_s is None or pdb_e is None or unp_s is None:
                    continue
                if pdb_s <= resnum <= pdb_e:
                    return acc, gene_id, unp_s + (resnum - pdb_s)
    except Exception:
        pass
    return None, None, None


def _annotate_residue(pdb_id: str, chain_id: str, resnum: int) -> dict:
    """
    Fetch UniProt functional annotations for (pdb_id, chain_id, resnum).
    Always returns a dict; errors are reported via the "error" key.
    """
    acc, gene_id, unp_pos = _sifts_uniprot(pdb_id, chain_id, resnum)
    if acc is None:
        return {"error": "No UniProt mapping for this chain/residue."}

    try:
        resp = requests.get(UNIPROT_ENTRY.format(acc), timeout=10,
                            headers={"Accept": "application/json"})
        if resp.status_code != 200:
            return {"error": f"UniProt API error ({resp.status_code})."}
        entry = resp.json()
    except Exception as e:
        return {"error": f"UniProt fetch failed: {e}"}

    # Protein / gene name
    rec = (entry.get("proteinDescription", {})
                .get("recommendedName", {})
                .get("fullName", {})
                .get("value", ""))
    genes = entry.get("genes", [])
    gene_name = gene_id
    if genes:
        gn = genes[0].get("geneName", {})
        if isinstance(gn, dict):
            gene_name = gn.get("value", gene_id)

    # Collect features that overlap unp_pos
    hits = []
    for feat in entry.get("features", []):
        ftype = feat.get("type", "")
        if ftype not in _FEAT_KEEP:
            continue
        loc   = feat.get("location", {})
        start = (loc.get("start") or {}).get("value")
        end   = (loc.get("end")   or {}).get("value")
        if start is None or end is None:
            continue
        if not (start <= unp_pos <= end):
            continue
        hits.append({
            "type":        ftype,
            "description": feat.get("description", ""),
            "range":       f"{start}–{end}" if start != end else str(start),
        })

    return {
        "accession":   acc,
        "gene":        gene_name,
        "protein":     rec,
        "unp_resnum":  unp_pos,
        "annotations": hits,
    }


@app.route("/api/annotations", methods=["POST"])
def api_annotations():
    d      = request.json or {}
    pdb_id = d.get("pdb_id", "").strip()
    chain  = d.get("chain",  "A").strip()
    resnum = int(d.get("resnum", 1))
    if not pdb_id:
        return jsonify({"error": "pdb_id required"}), 400
    return jsonify(_annotate_residue(pdb_id, chain, resnum))


# ---------------------------------------------------------------------------
# Biochemical analysis
# ---------------------------------------------------------------------------

_CHARGED_POS = {"ARG", "LYS", "HIS"}
_CHARGED_NEG = {"ASP", "GLU"}
_POLAR       = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "TRP"}
_NONPOLAR    = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "PRO", "GLY"}

_CHARGE_ATOMS = {
    "ARG": ["NH1", "NH2", "NE"],
    "LYS": ["NZ"],
    "HIS": ["ND1", "NE2"],
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}

_POLAR_HEAVY = {"N","O","OG","OG1","OD1","OD2","OE1","OE2",
                "ND1","NE2","NZ","NH1","NH2","NE","ND2","OH","NE1"}

def _btype(resname: str) -> str:
    if resname in _CHARGED_POS: return "positive"
    if resname in _CHARGED_NEG: return "negative"
    if resname in _POLAR:       return "polar"
    if resname in _NONPOLAR:    return "nonpolar"
    return "other"

def _charge(resname: str) -> int:
    return +1 if resname in _CHARGED_POS else (-1 if resname in _CHARGED_NEG else 0)

def _v(atom) -> np.ndarray:
    """Return (3,) float64 coords without using BioPython Vector (numpy 2.x safe)."""
    # atom.coord works for both Atom (direct attribute) and DisorderedAtom
    # (__getattr__ delegates to disordered_get().coord).  Bypassing Vector
    # avoids the numpy 2.x `len(atom)` probe that raised TypeError.
    try:
        return np.asarray(atom.coord, dtype=np.float64)
    except AttributeError:
        return atom.get_vector().get_array()

def _d(a, b) -> float:
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)))

def _ang(p1, vertex, p2) -> float:
    a = np.asarray(p1) - np.asarray(vertex)
    b = np.asarray(p2) - np.asarray(vertex)
    cos = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12)
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))

def _res_atom(residue, name):
    """Safe atom lookup; resolves DisorderedAtom to its selected child."""
    try:
        a = residue[name]
        # DisorderedAtom has disordered_get(); unwrap to a concrete Atom so
        # numpy never sees the disordered wrapper.
        if hasattr(a, 'disordered_get'):
            a = a.disordered_get()
        return a
    except KeyError:
        return None

def _charge_centroid(residue):
    names  = _CHARGE_ATOMS.get(residue.get_resname().strip(), [])
    atoms  = (_res_atom(residue, n) for n in names if n in residue)
    coords = [_v(a) for a in atoms if a is not None]
    return np.mean(coords, axis=0) if coords else None

def _find_res(model, chain_id, resnum):
    for ch in model.get_chains():
        if ch.id == chain_id:
            for r in ch.get_residues():
                if r.get_id() == (' ', resnum, ' '):
                    return r
    for ch in model.get_chains():          # fallback: ignore chain
        for r in ch.get_residues():
            if r.get_id() == (' ', resnum, ' '):
                return r
    return None

def _analyze_biochemistry(pdb_string: str, chain_id: str, resnum: int,
                           prox_cut=6.0, salt_cut=4.5, hb_cut=3.5, hydro_cut=5.0) -> dict:
    """
    Analyse the biochemical environment of (chain_id, resnum):
      • proximal residues classified by charge/polarity
      • charge interactions (salt bridges) with centroid distance and approach angle
      • H-bond count (heavy-atom N/O criterion, no explicit H required)
      • hydrophobic contacts
    """
    parser = PDBParser(QUIET=True)
    model  = parser.get_structure("S", io.StringIO(pdb_string))[0]

    mut_r  = _find_res(model, chain_id, resnum)
    if mut_r is None:
        return {"error": f"Residue {chain_id}{resnum} not found"}

    rname = mut_r.get_resname().strip()
    ca_m  = _res_atom(mut_r, "CA")
    if ca_m is None:
        return {"error": "No Cα at mutation site"}
    ca_pos = _v(ca_m)

    all_prot = [r for ch in model.get_chains()
                for r in ch.get_residues()
                if r.get_id()[0] == ' ' and r.get_resname().strip() in AA1]

    # ---------- proximal residues ----------
    prox = []
    for r in all_prot:
        if r is mut_r: continue
        ca = _res_atom(r, "CA")
        if ca is None: continue
        d = _d(_v(ca), ca_pos)
        if d <= prox_cut:
            prox.append((r, round(d, 2)))
    prox.sort(key=lambda x: x[1])

    proximal_info = [{
        "resnum":  r.get_id()[1],
        "resname": r.get_resname().strip(),
        "chain":   r.get_parent().id,
        "type":    _btype(r.get_resname().strip()),
        "dist_CA": d,
    } for r, d in prox]

    # ---------- charge interactions ----------
    c_mut = _charge_centroid(mut_r)
    charged_prox = [(r, d) for r, d in prox if _charge(r.get_resname().strip()) != 0]
    ci_list = []

    if c_mut is not None:
        for r, _ in charged_prox:
            if _charge(rname) * _charge(r.get_resname().strip()) >= 0: continue
            c_r = _charge_centroid(r)
            if c_r is None: continue
            dist = _d(c_mut, c_r)
            if dist <= salt_cut:
                angle = _ang(ca_pos, c_mut, c_r)
                ci_list.append({
                    "res1": f"{rname}{resnum}", "res2": f"{r.get_resname().strip()}{r.get_id()[1]}",
                    "dist_Å": round(dist,2), "angle_deg": round(angle,1),
                })

    seen_ci = set()
    for i, (r1, _) in enumerate(charged_prox):
        for j, (r2, _) in enumerate(charged_prox):
            if j <= i: continue
            if _charge(r1.get_resname().strip()) * _charge(r2.get_resname().strip()) >= 0: continue
            k = (id(r1), id(r2))
            if k in seen_ci: continue
            seen_ci.add(k)
            c1, c2 = _charge_centroid(r1), _charge_centroid(r2)
            if c1 is None or c2 is None: continue
            dist = _d(c1, c2)
            if dist <= salt_cut:
                ca1 = _res_atom(r1, "CA")
                angle = _ang(_v(ca1), c1, c2) if ca1 else 0.0
                ci_list.append({
                    "res1": f"{r1.get_resname().strip()}{r1.get_id()[1]}",
                    "res2": f"{r2.get_resname().strip()}{r2.get_id()[1]}",
                    "dist_Å": round(dist,2), "angle_deg": round(angle,1),
                })
    ci_list.sort(key=lambda x: x["dist_Å"])

    # ---------- hydrogen bonds (heavy-atom N/O criterion) ----------
    def _polar_atoms(res):
        result = []
        for a in res.get_atoms():
            if hasattr(a, 'disordered_get'):
                a = a.disordered_get()
            name = a.get_name().strip()
            if name in _POLAR_HEAVY and not name.startswith('H'):
                result.append(a)
        return result

    hbonds = []
    seen_hb: set = set()
    for prox_r, _ in prox:
        for a1 in _polar_atoms(mut_r):
            for a2 in _polar_atoms(prox_r):
                d = _d(_v(a1), _v(a2))
                if d <= hb_cut:
                    key = tuple(sorted([id(a1), id(a2)]))
                    if key not in seen_hb:
                        seen_hb.add(key)
                        hbonds.append({
                            "atom1":  f"{rname}{resnum}:{a1.get_name().strip()}",
                            "atom2":  f"{prox_r.get_resname().strip()}{prox_r.get_id()[1]}:{a2.get_name().strip()}",
                            "dist_Å": round(d, 2),
                        })
    hbonds.sort(key=lambda x: x["dist_Å"])

    # ---------- hydrophobic contacts ----------
    hydrophobic = []
    cb_m = _res_atom(mut_r, "CB") or _res_atom(mut_r, "CA")
    if rname in _NONPOLAR and cb_m is not None:
        for prox_r, _ in prox:
            if prox_r.get_resname().strip() not in _NONPOLAR: continue
            cb_p = _res_atom(prox_r, "CB") or _res_atom(prox_r, "CA")
            if cb_p is None: continue
            d = _d(_v(cb_m), _v(cb_p))
            if d <= hydro_cut:
                hydrophobic.append({
                    "partner": f"{prox_r.get_resname().strip()}{prox_r.get_id()[1]}",
                    "dist_Å": round(d, 2),
                })
    hydrophobic.sort(key=lambda x: x["dist_Å"])

    return {
        "mutation_site":       {"resname": rname, "type": _btype(rname), "charge": _charge(rname)},
        "proximal_residues":   proximal_info,
        "charge_interactions": ci_list,
        "hbonds":              hbonds,
        "hydrophobic_contacts":hydrophobic,
    }


def _diff_biochemistry(orig: dict, mut: dict) -> dict:
    def key(ci): return tuple(sorted([ci["res1"], ci["res2"]]))
    o_ci = {key(c): c for c in orig.get("charge_interactions", [])}
    m_ci = {key(c): c for c in mut.get("charge_interactions", [])}
    changed = []
    for k in set(o_ci) & set(m_ci):
        oc, mc = o_ci[k], m_ci[k]
        if abs(oc["dist_Å"] - mc["dist_Å"]) > 0.3 or abs(oc["angle_deg"] - mc["angle_deg"]) > 5:
            changed.append({"orig": oc, "mut": mc})
    return {
        "charge_delta":        mut["mutation_site"]["charge"] - orig["mutation_site"]["charge"],
        "type_change":         f'{orig["mutation_site"]["type"]} → {mut["mutation_site"]["type"]}',
        "charge_gained":       [m_ci[k] for k in set(m_ci) - set(o_ci)],
        "charge_lost":         [o_ci[k] for k in set(o_ci) - set(m_ci)],
        "charge_changed":      changed,
        "hbonds_orig":         len(orig.get("hbonds", [])),
        "hbonds_mut":          len(mut.get("hbonds", [])),
        "hbonds_delta":        len(mut.get("hbonds", [])) - len(orig.get("hbonds", [])),
        "hydrophobic_orig":    len(orig.get("hydrophobic_contacts", [])),
        "hydrophobic_mut":     len(mut.get("hydrophobic_contacts", [])),
        "hydrophobic_delta":   len(mut.get("hydrophobic_contacts", [])) - len(orig.get("hydrophobic_contacts", [])),
    }


# ---------------------------------------------------------------------------
# Structural analysis: Ramachandran · backbone H-bonds · RMSF · ΔΔG
# ---------------------------------------------------------------------------

# Ramachandran region boundaries (non-Gly, non-Pro general case)
def _rama_region(phi, psi, resname):
    if phi is None or psi is None:
        return "undefined"
    if resname == "GLY":
        return "allowed"          # Gly has no β-carbon; all quadrants accessible
    if resname == "PRO":
        if -90 <= phi <= -30:  return "favored"
        if -120 <= phi <= 0:   return "allowed"
        return "outlier"
    # General residues
    if -90 <= phi <= -30  and -60 <= psi <= 5:   return "core α"
    if -170 <= phi <= -55 and (100 <= psi <= 180 or -180 <= psi <= -150):
        return "core β"
    if 30 <= phi <= 90 and 10 <= psi <= 80:      return "αL"
    if phi <= 0:                                  return "allowed"
    return "outlier"


def _phi_psi_window(pdb_string: str, chain_id: str, resnum: int, window: int = 3) -> dict:
    """Return {resnum: {phi, psi, resname, region}} for residues in [resnum±window]."""
    parser  = PDBParser(QUIET=True)
    model   = parser.get_structure("S", io.StringIO(pdb_string))[0]
    builder = PPBuilder()
    target  = set(range(resnum - window, resnum + window + 1))
    result  = {}
    for pp in builder.build_peptides(model):
        phi_psi = pp.get_phi_psi_list()
        for i, res in enumerate(pp):
            if res.get_parent().id != chain_id:
                continue
            rid = res.get_id()[1]
            if rid not in target:
                continue
            phi, psi   = phi_psi[i]
            phi_d = round(math.degrees(phi), 1) if phi is not None else None
            psi_d = round(math.degrees(psi), 1) if psi is not None else None
            rn    = res.get_resname().strip()
            result[rid] = {
                "resnum":  rid,
                "resname": rn,
                "phi":     phi_d,
                "psi":     psi_d,
                "region":  _rama_region(phi_d, psi_d, rn),
            }
    return result


def _backbone_nhbonds(pdb_string: str, chain_id: str, resnum: int,
                      cut: float = 3.5) -> tuple:
    """
    Find backbone N–H…O=C H-bonds at the mutation site.

    Returns (donated, accepted, is_pro):
      donated  – bonds where the site N donates to a partner C=O
      accepted – bonds where the site C=O accepts from a partner N
      is_pro   – True when the site residue is PRO (N cannot donate)
    """
    parser  = PDBParser(QUIET=True)
    model   = parser.get_structure("S", io.StringIO(pdb_string))[0]
    site    = _find_res(model, chain_id, resnum)
    if site is None:
        return [], [], False

    site_rn = site.get_resname().strip()
    is_pro  = site_rn == "PRO"
    site_N  = _res_atom(site, "N")
    site_O  = _res_atom(site, "O")

    all_res = [r for ch in model.get_chains()
               for r in ch.get_residues()
               if r.get_id()[0] == " " and r.get_resname().strip() in AA1
               and r is not site]

    donated, accepted = [], []
    for r in all_res:
        r_rn  = r.get_resname().strip()
        r_rid = r.get_id()[1]
        r_N   = _res_atom(r, "N")
        r_O   = _res_atom(r, "O")

        if not is_pro and site_N and r_O:
            d = _d(_v(site_N), _v(r_O))
            if d <= cut:
                donated.append({"partner_resnum": r_rid, "partner_resname": r_rn,
                                 "dist_Å": round(d, 2)})
        if r_rn != "PRO" and r_N and site_O:
            d = _d(_v(r_N), _v(site_O))
            if d <= cut:
                accepted.append({"partner_resnum": r_rid, "partner_resname": r_rn,
                                  "dist_Å": round(d, 2)})

    return donated, accepted, is_pro


def _local_bfactor_rmsf(pdb_string: str, chain_id: str, resnum: int,
                         window: int = 5) -> dict:
    """
    Mean backbone B-factor per residue in [resnum±window].
    Converted to RMSF via  RMSF = sqrt(3B / 8π²)  [Å].
    """
    parser = PDBParser(QUIET=True)
    model  = parser.get_structure("S", io.StringIO(pdb_string))[0]
    target = set(range(resnum - window, resnum + window + 1))
    result = {}
    for ch in model.get_chains():
        if ch.id != chain_id:
            continue
        for r in ch.get_residues():
            if r.get_id()[0] != " ":
                continue
            rid = r.get_id()[1]
            if rid not in target:
                continue
            bvals = [_res_atom(r, nm).get_bfactor()
                     for nm in ("N", "CA", "C", "O")
                     if _res_atom(r, nm) is not None]
            if not bvals:
                continue
            b_mean = sum(bvals) / len(bvals)
            result[rid] = {
                "resnum":       rid,
                "resname":      r.get_resname().strip(),
                "bfactor":      round(b_mean, 2),
                "rmsf_bfactor": round(math.sqrt(max(3 * b_mean / (8 * math.pi**2), 0)), 3),
            }
    return result


def _local_ss_content(phi_psi_dict: dict) -> dict:
    """Count helix / strand / coil from a φ/ψ window dict."""
    counts = {"helix": 0, "strand": 0, "coil": 0, "undefined": 0}
    for info in phi_psi_dict.values():
        reg = info.get("region", "undefined")
        if "α" in reg:     counts["helix"]   += 1
        elif "β" in reg:   counts["strand"]  += 1
        elif reg == "undefined": counts["undefined"] += 1
        else:              counts["coil"]    += 1
    total = counts["helix"] + counts["strand"] + counts["coil"]
    if total == 0:
        return {**counts, "helix_frac": 0.0, "strand_frac": 0.0, "coil_frac": 0.0}
    return {**counts,
            "helix_frac":  round(counts["helix"]  / total, 3),
            "strand_frac": round(counts["strand"] / total, 3),
            "coil_frac":   round(counts["coil"]   / total, 3)}


# Empirical scales for ΔΔG
_KD = {                        # Kyte-Doolittle hydrophobicity
    "ILE":4.5,"VAL":4.2,"LEU":3.8,"PHE":2.8,"CYS":2.5,"MET":1.9,"ALA":1.8,
    "GLY":-0.4,"THR":-0.7,"SER":-0.8,"TRP":-0.9,"TYR":-1.3,"PRO":-1.6,
    "HIS":-3.2,"GLU":-3.5,"GLN":-3.5,"ASP":-3.5,"ASN":-3.5,"LYS":-3.9,"ARG":-4.5,
}
_NROT = {                      # approximate rotamer count per residue
    "ALA":1,"GLY":1,"PRO":1,"VAL":3,"THR":3,"SER":3,"CYS":3,
    "ILE":7,"LEU":9,"ASP":6,"ASN":6,"HIS":8,"PHE":4,"TYR":4,"TRP":4,
    "GLU":12,"GLN":12,"MET":9,"LYS":27,"ARG":27,
}


def _ddg_estimate(ref_aa: str, alt_aa: str,
                  bb_hbond_delta: int,
                  orig_site_region: str) -> dict:
    """
    Rough empirical ΔΔG (kcal/mol).  Positive = destabilising.

    Components:
      hydrophobic  — KD scale × 0.15 kcal/mol per unit (50% burial assumed)
      bb_hbond     — ±1.5 kcal/mol per backbone H-bond lost/gained
      sc_entropy   — RT·Δln(N_rot): more rotamers in mutant → less stable folded state
      proline      — backbone-geometry penalty/bonus specific to Pro
    """
    ref3 = AA3.get(ref_aa.upper(), ref_aa.upper())
    alt3 = AA3.get(alt_aa.upper(), alt_aa.upper())

    ddg_hydro   = -(_KD.get(alt3, 0) - _KD.get(ref3, 0)) * 0.15
    ddg_bb      = -bb_hbond_delta * 1.5
    ddg_entropy = (math.log(_NROT.get(alt3, 3)) - math.log(_NROT.get(ref3, 3))) * 0.6
    ddg_pro     = 0.0
    if alt3 == "PRO" and ref3 != "PRO":
        if "α" in orig_site_region:   ddg_pro = +3.5
        elif "β" in orig_site_region: ddg_pro = +2.0
        else:                         ddg_pro = +0.5
    elif ref3 == "PRO" and alt3 != "PRO":
        ddg_pro = -2.0 if "α" in orig_site_region else -0.5

    total = ddg_hydro + ddg_bb + ddg_entropy + ddg_pro
    if   total >  1.0: verdict = "likely destabilising"
    elif total < -1.0: verdict = "likely stabilising"
    else:              verdict = "approximately neutral"

    return {
        "ddg_kcal_mol":  round(total, 2),
        "components": {
            "hydrophobic":    round(ddg_hydro,   2),
            "backbone_hbond": round(ddg_bb,      2),
            "sc_entropy":     round(ddg_entropy, 2),
            "proline":        round(ddg_pro,     2),
        },
        "verdict":    verdict,
        "confidence": "rough (±2 kcal/mol); positive = destabilising",
    }


def _structural_analysis(original_pdb: str, mutant_pdb: str,
                          chain_id: str, resnum: int,
                          ref_aa: str, alt_aa: str,
                          rmsf_mutant: dict | None = None) -> dict:
    """
    Master function: Ramachandran angles, backbone H-bond loss, RMSF proxy,
    local secondary-structure content, and empirical ΔΔG estimate.
    """
    # ── Ramachandran ──────────────────────────────────────────────────────────
    orig_rama = _phi_psi_window(original_pdb, chain_id, resnum, window=3)
    mut_rama  = _phi_psi_window(mutant_pdb,   chain_id, resnum, window=3)
    rama = {
        "original":      sorted(orig_rama.values(), key=lambda x: x["resnum"]),
        "mutant":        sorted(mut_rama.values(),  key=lambda x: x["resnum"]),
        "outliers_orig": [r for r in orig_rama.values() if r["region"] == "outlier"],
        "outliers_mut":  [r for r in mut_rama.values()  if r["region"] == "outlier"],
    }

    # ── Backbone H-bonds ──────────────────────────────────────────────────────
    o_don, o_acc, o_pro = _backbone_nhbonds(original_pdb, chain_id, resnum)
    m_don, m_acc, m_pro = _backbone_nhbonds(mutant_pdb,   chain_id, resnum)

    def _hb_set(lst):  return {h["partner_resnum"]: h for h in lst}
    od, md = _hb_set(o_don), _hb_set(m_don)
    oa, ma = _hb_set(o_acc), _hb_set(m_acc)

    note = None
    if m_pro and not o_pro:
        note = "PRO cannot donate backbone N–H; H-bond(s) to C=O i−4 are lost."
    elif o_pro and not m_pro:
        note = "PRO replaced: site N regains ability to donate backbone H-bond."

    bb_hbonds = {
        "orig_donates":   o_don, "mut_donates":    m_don,
        "orig_accepts":   o_acc, "mut_accepts":    m_acc,
        "donated_lost":   [od[k] for k in set(od) - set(md)],
        "donated_gained": [md[k] for k in set(md) - set(od)],
        "accepted_lost":  [oa[k] for k in set(oa) - set(ma)],
        "accepted_gained":[ma[k] for k in set(ma) - set(oa)],
        "mut_is_pro":  m_pro, "orig_is_pro": o_pro,
        "note":        note,
    }
    bb_delta = (len(m_don) + len(m_acc)) - (len(o_don) + len(o_acc))

    # ── Flexibility (RMSF) ───────────────────────────────────────────────────
    orig_bf = _local_bfactor_rmsf(original_pdb, chain_id, resnum, window=5)
    per_res = []
    for rn, info in sorted(orig_bf.items()):
        key   = f"{chain_id}:{rn}"
        entry = {
            "resnum":       rn,
            "resname":      info["resname"],
            "rmsf_orig":    info["rmsf_bfactor"],
            "bfactor_orig": info["bfactor"],
        }
        if rmsf_mutant:
            entry["rmsf_mut"] = rmsf_mutant.get(key)
        per_res.append(entry)

    orig_rmsf_vals = [v["rmsf_bfactor"] for v in orig_bf.values()]
    mut_rmsf_vals  = ([rmsf_mutant.get(f"{chain_id}:{rn}")
                       for rn in orig_bf if rmsf_mutant
                       and rmsf_mutant.get(f"{chain_id}:{rn}") is not None]
                      if rmsf_mutant else [])

    flexibility = {
        "per_residue":   per_res,
        "mean_rmsf_orig": round(float(np.mean(orig_rmsf_vals)), 3) if orig_rmsf_vals else None,
        "mean_rmsf_mut":  round(float(np.mean(mut_rmsf_vals)),  3) if mut_rmsf_vals  else None,
        "source_orig":    "B-factor → RMSF = √(3B/8π²)",
        "source_mut":     "Cα RMSF from NVT dynamics" if rmsf_mutant else "not available",
    }

    # ── Secondary structure ───────────────────────────────────────────────────
    orig_ss = _local_ss_content(orig_rama)
    mut_ss  = _local_ss_content(mut_rama)
    ss = {
        "original": orig_ss, "mutant": mut_ss,
        "helix_delta":  round(mut_ss["helix_frac"]  - orig_ss["helix_frac"],  3),
        "strand_delta": round(mut_ss["strand_frac"] - orig_ss["strand_frac"], 3),
    }

    # ── ΔΔG ──────────────────────────────────────────────────────────────────
    site_region = orig_rama.get(resnum, {}).get("region", "allowed")
    ddg = _ddg_estimate(ref_aa, alt_aa, bb_delta, site_region)

    return {"rama": rama, "bb_hbonds": bb_hbonds,
            "flexibility": flexibility, "ss": ss, "ddg": ddg}


@app.route("/api/report", methods=["POST"])
def api_report():
    d      = request.json or {}
    chain  = d.get("chain",  "A")
    resnum = int(d.get("resnum", 1))
    ref_aa = d.get("ref_aa", "")
    alt_aa = d.get("alt_aa", "")
    rmsf_mutant = d.get("rmsf_mutant") or {}
    try:
        orig_a = _analyze_biochemistry(d.get("original_pdb",""), chain, resnum)
        mut_a  = _analyze_biochemistry(d.get("mutant_pdb", ""), chain, resnum)
        if "error" in orig_a: return jsonify(orig_a), 400
        if "error" in mut_a:  return jsonify(mut_a),  400
        structural = _structural_analysis(
            d.get("original_pdb",""), d.get("mutant_pdb",""),
            chain, resnum, ref_aa, alt_aa,
            rmsf_mutant=rmsf_mutant or None,
        )
        return jsonify({
            "original":   orig_a,
            "mutant":     mut_a,
            "diff":       _diff_biochemistry(orig_a, mut_a),
            "structural": structural,
            "ref_aa":     ref_aa,
            "alt_aa":     alt_aa,
            "resnum":     resnum,
        })
    except Exception as e:
        import traceback, sys
        tb = traceback.format_exc()
        print(tb, file=sys.stderr, flush=True)
        return jsonify({"error": str(e), "trace": tb}), 500


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/structures", methods=["POST"])
def api_structures():
    protein = (request.json or {}).get("protein", "").strip()
    if not protein:
        return jsonify({"error": "protein is required"}), 400
    try:
        return jsonify(_list_structures(protein))
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/mutate", methods=["POST"])
def api_mutate():
    data     = request.json or {}
    protein  = data.get("protein",  "").strip()
    mutation = data.get("mutation", "").strip()
    # optional pre-selected structure from the picker
    sel_id     = data.get("pdb_id",  "").strip()
    sel_chain  = data.get("chain",   "").strip()
    sel_source = data.get("source",  "crystal").strip()

    def _sse(obj):
        return f"data: {json.dumps(obj)}\n\n"

    def generate():
        # ---- parse ----
        try:
            ref_aa, resnum, alt_aa = _parse_residue_number(mutation)
        except ValueError as e:
            yield _sse({"error": str(e)}); return

        if not protein or not mutation:
            yield _sse({"error": "protein and mutation are required"}); return

        # ---- find / use structure ----
        if sel_id:
            pdb_id = sel_id
            chain  = sel_chain or "A"
            source = sel_source
            yield _sse({"step": f"Using selected structure {pdb_id}…"})
        else:
            yield _sse({"step": f"Searching for best structure for '{protein}'…"})
            try:
                pdb_id, chain, source = _best_pdb(protein)
            except ValueError as e:
                yield _sse({"error": str(e)}); return

        yield _sse({"step": f"Downloading {pdb_id} ({source})…"})
        try:
            raw_pdb = _fetch_pdb_string(pdb_id, source)
        except Exception as e:
            yield _sse({"error": f"Failed to download structure: {e}"}); return

        # ---- locate residue for original highlight ----
        try:
            actual_chain, actual_resnum, _ = _locate_residue(
                _protein_only_pdb(raw_pdb), ref_aa, resnum
            )
        except ValueError as e:
            yield _sse({"error": f"Residue not found: {e}"}); return

        original_pdb = _clean_pdb(raw_pdb)
        original_pdb = _highlight_residue(original_pdb, actual_chain, actual_resnum)

        # ---- mutate + minimise (with per-step progress) ----
        _stop_event.clear()
        q = queue.Queue()

        def run_pipeline():
            try:
                result = _fix_and_minimise(
                    raw_pdb, ref_aa, resnum, alt_aa,
                    on_progress=lambda msg: q.put({"step": msg}),
                    stop_event=_stop_event,
                )
                q.put({"_result": result})
            except _Stopped:
                q.put({"cancelled": True})
            except Exception as e:
                q.put({"error": f"Mutation / minimisation failed: {e}"})

        t = threading.Thread(target=run_pipeline, daemon=True)
        t.start()

        while True:
            item = q.get()
            if "error" in item:
                yield _sse(item); return
            if "cancelled" in item:
                yield _sse({"cancelled": True}); return
            if "_result" in item:
                mutant_raw, actual_chain, actual_resnum, rmsf_mutant = item["_result"]
                break
            yield _sse(item)

        mutant_pdb = _clean_pdb(mutant_raw)
        mutant_pdb = _normalize_resnames(mutant_pdb)
        mutant_pdb = _highlight_residue(mutant_pdb, actual_chain, actual_resnum)

        result_payload = {
            "done":           True,
            "pdb_id":         pdb_id,
            "chain":          actual_chain,
            "source":         source,
            "ref_aa":         ref_aa,
            "alt_aa":         alt_aa,
            "resnum":         actual_resnum,
            "original_pdb":   original_pdb,
            "mutant_pdb":     mutant_pdb,
            "mutation_label": f"{ref_aa}{resnum}{alt_aa}",
            "rmsf_mutant":    rmsf_mutant,
        }

        # ---- persist session ----
        try:
            ts   = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            sid  = f"{ts}_{pdb_id}_{ref_aa}{resnum}{alt_aa}"
            sess = dict(result_payload, protein=protein, session_id=sid, saved_at=ts)
            with open(os.path.join(SESSIONS_DIR, f"{sid}.json"), "w") as _f:
                json.dump(sess, _f)
        except Exception:
            pass   # never block the response for a save failure

        yield _sse(result_payload)

    return Response(
        stream_with_context(generate()),
        mimetype="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


@app.route("/api/sessions", methods=["GET"])
def api_sessions():
    """List saved sessions, newest first, without the large PDB strings."""
    sessions = []
    for fname in sorted(os.listdir(SESSIONS_DIR), reverse=True):
        if not fname.endswith(".json"):
            continue
        try:
            with open(os.path.join(SESSIONS_DIR, fname)) as f:
                s = json.load(f)
            sessions.append({
                "session_id":     s.get("session_id", fname[:-5]),
                "protein":        s.get("protein", ""),
                "pdb_id":         s.get("pdb_id", ""),
                "source":         s.get("source", ""),
                "chain":          s.get("chain", ""),
                "mutation_label": s.get("mutation_label", ""),
                "ref_aa":         s.get("ref_aa", ""),
                "alt_aa":         s.get("alt_aa", ""),
                "resnum":         s.get("resnum", 0),
                "saved_at":       s.get("saved_at", ""),
            })
        except Exception:
            continue
    return jsonify(sessions)


@app.route("/api/sessions/<session_id>", methods=["GET"])
def api_session_load(session_id):
    """Return a saved session including PDB strings."""
    path = os.path.join(SESSIONS_DIR, f"{session_id}.json")
    if not os.path.isfile(path):
        return jsonify({"error": "Session not found"}), 404
    with open(path) as f:
        return jsonify(json.load(f))


@app.route("/api/stop", methods=["POST"])
def api_stop():
    _stop_event.set()
    return jsonify({"ok": True})


def main():
    import webbrowser, threading as _t
    _t.Timer(1.5, lambda: webbrowser.open("http://localhost:5050")).start()
    app.run(host="127.0.0.1", port=5050)


if __name__ == "__main__":
    main()
