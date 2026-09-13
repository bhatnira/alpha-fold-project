#!/usr/bin/env python3
"""7. Medicinal-chemistry QC (Chen et al. framework) applied to panel poses.

Metrics (no experimental pose exists -> report internal consistency only):
  * ligand chemical validity (RDKit parse of panel SMILES)
  * stereochemistry (undefined/warning stereo centers in generated poses)
  * protein-ligand clashes (heavy-atom pairs < CLASH_CUTOFF)
  * pocket occupancy (ligand heavy atoms within BURIED_CUTOFF of pocket CAs)
  * pharmacophore feature conservation across models
  * whole-ligand RMSD vs pose centroid does NOT imply experimental accuracy

Output: outputs/07_medchem_qc.json
"""

import json
import warnings
from collections import defaultdict
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np

from config import RESULTS, COFOLDING_RESULTS, MODULATOR_CSV, CLASH_CUTOFF

try:
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    HAS_RDKit = True
except Exception:
    HAS_RDKit = False


def parse_name(name):
    parts = name.split("_")
    if len(parts) >= 4:
        return parts[1], parts[2], int(parts[3].replace("cpd", ""))
    return None, None, None


def main():
    from common import parse_pdb, split_ligand, parse_ligand_any_chain

    panel = COFOLDING_RESULTS / "full_panel"
    n_models = n_clash_free = n_valid_pos = 0
    n_undefined = 0
    clash_total = 0
    occupancy = []
    pharmacophore_positions = defaultdict(list)

    for cond_dir in sorted(panel.iterdir()) if panel.exists() else []:
        if not cond_dir.is_dir():
            continue
        stoich, site_str, cpd = parse_name(cond_dir.name)
        if site_str is None:
            continue
        pred_dir = cond_dir / f"boltz_results_{cond_dir.name}" / "predictions" / cond_dir.name
        if not pred_dir.exists():
            continue
        for pdb in sorted(pred_dir.glob(f"{cond_dir.name}_model_*.pdb")):
            atoms = parse_pdb(pdb)
            protein, ligand = split_ligand(atoms)
            if not ligand:
                protein, ligand = parse_ligand_any_chain(atoms)
            if not (ligand and protein):
                continue
            n_models += 1
            lig_coords = np.vstack([a.coord for a in ligand])
            prot_coords = np.vstack([a.coord for a in protein])
            d = np.linalg.norm(lig_coords[:, None, :] - prot_coords[None, :, :], axis=2)
            n_clash = int((d.min(axis=1) < CLASH_CUTOFF).sum())
            clash_total += n_clash
            n_clash_free += (n_clash == 0)
            ca = np.vstack([a.coord for a in protein if a.atomname == "CA"])
            if len(ca):
                buried = float((np.linalg.norm(lig_coords[:, None, :] - ca[None, :, :], axis=2).min(axis=1) <= 6.0).mean())
                occupancy.append(buried)
            pharm = np.array([a.coord for a in ligand if a.atomname.startswith(("O", "N"))])
            if len(pharm):
                centroids = pharm.mean(axis=0)
                pharmacophore_positions[cpd].append(centroids)

    pharm_spread = {}
    for cpd, pos in pharmacophore_positions.items():
        if len(pos) < 2:
            pharm_spread[cpd] = None
            continue
        pts = np.vstack(pos)
        cent = pts.mean(axis=0)
        # pooled RMS distance of the per-model pharmacophore centroid from the
        # within-compound mean centroid = genuine pose-consistency metric
        rms = float(np.sqrt(np.mean(((pts - cent) ** 2).sum(axis=1))))
        pharm_spread[cpd] = round(rms, 3)

    # RDKit quality + stereo on the experimental SMILES
    rdkit_qc = {"available": HAS_RDKit}
    import csv
    if HAS_RDKit:
        parsed = valid = with_stereo = 0
        bad = []
        with open(MODULATOR_CSV) as f:
            for row in csv.DictReader(f):
                parsed += 1
                mol = Chem.MolFromSmiles(row["Smiles"])
                if mol is None:
                    bad.append((row["Identifier"], "unparseable"))
                    continue
                valid += 1
                if any(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for atom in mol.GetAtoms()):
                    with_stereo += 1
        rdkit_qc.update({
            "n_smiles": parsed, "n_valid": valid,
            "n_with_stereo": with_stereo,
            "unparseable": bad,
        })

    result = {
        "framework": "Chen et al., A Medicinal Chemistry-Centered Evaluation of AF3 and Boltz-2",
        "n_models_analyzed": n_models,
        "models_clash_free": n_clash_free,
        "fraction_clash_free": round(n_clash_free / max(1, n_models), 3),
        "mean_clashing_ligand_atoms_per_model": round(clash_total / max(1, n_models), 3),
        "clash_cutoff_A": CLASH_CUTOFF,
        "mean_pocket_occupancy_fraction": round(float(np.mean(occupancy)), 3) if occupancy else None,
        "pharmacophore_centroid_spread_A_per_cpds": pharm_spread,
        "rdkit": rdkit_qc,
        "note": "No experimental a9a10 PAM-bound structure exists; pose-level numbers are internal consistency only, not accuracy.",
    }

    with open(RESULTS / "07_medchem_qc.json", "w") as f:
        json.dump(result, f, indent=2)

    print(f"models: {n_models}; clash-free {n_clash_free} ({result['fraction_clash_free']:.1%})")
    print(f"mean clashing ligand atoms/model: {result['mean_clashing_ligand_atoms_per_model']}")
    print(f"pocket occupancy: {result['mean_pocket_occupancy_fraction']}")
    print(f"pharmacophore spread (top): {list(pharm_spread.items())[:5]}")
    print(f"RDKit: {rdkit_qc}")
    print(f"wrote {RESULTS / '07_medchem_qc.json'}")


if __name__ == "__main__":
    main()