#!/usr/bin/env python3
"""Phase V — Flexible docking attempt.

Vina 1.2.7 Python SWIG binding does not support flex receptor (set_receptor
with 2 args crashes). This script does rigid docking with proper PDBQT
format and records which residues would be flexible.
"""

import csv, sys, os, math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_docking import (
    prepare_receptor_pdbqt, prepare_ligand_pdbqt, run_vina_docking, RESULTS
)

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent

FLEX_CONFIG = {
    "exhaustiveness": 16,
    "n_poses": 10,
    "energy_range": 4,
    "box_padding_A": 10,
}


def get_residues_from_site(site_row):
    contacts = site_row.get("recurrent_contacts_ge50pct", "")
    if not contacts:
        return []
    chain_map = {"alpha9": "B", "alpha10": "C"}
    residues = []
    for item in contacts.split(";"):
        item = item.strip()
        if ":" in item:
            subunit, res_info = item.split(":", 1)
            chain = chain_map.get(subunit, "A")
            if res_info.strip().isdigit():
                residues.append(f"{chain}:{res_info.strip()}")
    return residues


def main():
    print("Phase V — Docking with flex residue annotation (rigid mode)")
    print("=" * 60)

    sites_csv = ROOT / "data/deliverable/05_structures/2to3/sites_2to3.csv"
    receptor_pdb = ROOT / "data/deliverable/05_structures/2to3/a9a10_2to3_apo_reference.pdb"
    if not sites_csv.exists() or not receptor_pdb.exists():
        print("Error: required input files not found")
        return

    sites = []
    with open(sites_csv) as f:
        for row in csv.DictReader(f):
            sites.append(row)
    print(f"Loaded {len(sites)} sites")

    ligands = [
        {"name": "ascorbate", "smiles": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O"},
        {"name": "acetate", "smiles": "CC(=O)[O-]"},
        {"name": "O-ethyl_ascorbate", "smiles": "CCOC(=O)C(=O)OC[C@@H](O)[C@@H]1OC(=O)C(O)=C1O"},
    ]

    rec_pdbqt = RESULTS / "receptor_rigid.pdbqt"
    if not rec_pdbqt.exists():
        print("Preparing receptor...")
        prepare_receptor_pdbqt(str(receptor_pdb), str(rec_pdbqt))

    for lig in ligands:
        lig_pdbqt = RESULTS / f"{lig['name']}.pdbqt"
        if not lig_pdbqt.exists() or lig_pdbqt.stat().st_size < 100:
            prepare_ligand_pdbqt(lig['smiles'], str(lig_pdbqt))

    results = []
    for site_row in sites[:8]:
        site_id = site_row["site"]
        copy = site_row["copy"]
        try:
            cx = float(site_row["centroid_x"])
            cy = float(site_row["centroid_y"])
            cz = float(site_row["centroid_z"])
        except (ValueError, KeyError):
            continue

        flex_residues = get_residues_from_site(site_row)
        try:
            spread = float(site_row.get("centroid_spread_A", 3.0))
        except ValueError:
            spread = 3.0
        box_dim = max(15, spread * 2 + FLEX_CONFIG["box_padding_A"] * 2)
        box_size = (box_dim, box_dim, box_dim)

        print(f"\nSite {site_id} copy {copy}  center=({cx:.1f}, {cy:.1f}, {cz:.1f})  flex={len(flex_residues)}")

        for lig in ligands:
            lig_pdbqt = RESULTS / f"{lig['name']}.pdbqt"
            if not lig_pdbqt.exists() or lig_pdbqt.stat().st_size < 100:
                continue

            out_pdb = RESULTS / f"site{site_id}_copy{copy}_{lig['name']}_docked.pdb"
            print(f"  Docking {lig['name']}...", end=" ", flush=True)

            success, energy_val = run_vina_docking(
                str(rec_pdbqt), str(lig_pdbqt),
                (cx, cy, cz), box_size, str(out_pdb),
                exhaustiveness=FLEX_CONFIG["exhaustiveness"],
                n_poses=FLEX_CONFIG["n_poses"],
            )

            results.append({
                "site_id": site_id,
                "copy": copy,
                "ligand": lig['name'],
                "smiles": lig['smiles'],
                "center_x": cx,
                "center_y": cy,
                "center_z": cz,
                "flex_residues": ";".join(flex_residues) if flex_residues else "",
                "box_size": box_dim,
                "success": success,
                "binding_energy": energy_val,
                "output": str(out_pdb) if success else "",
            })
            print(f"{'OK' if success else 'FAILED'}" + (f"  E={energy_val:.2f}" if energy_val else ""))

    if results:
        csv_path = OUT / "flexible_docking_results.csv"
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults: {csv_path}")

    n_ok = sum(1 for r in results if r["success"])
    print(f"\n{n_ok}/{len(results)} dockings successful")


if __name__ == "__main__":
    main()
