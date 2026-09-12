#!/usr/bin/env python3
"""
Step 4: Structural Features - Extract structural features from PDB files:
contact counts, H-bond counts, residue distances between ligand and protein.
Parse PDB ATOM/HETATM records. Save to features/structural/.
"""
import csv, json, os, re, math
from pathlib import Path
import numpy as np

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
FEATURES = MECHANISTIC / "features"
STRUCTURAL = FEATURES / "structural"
STRUCTURAL.mkdir(parents=True, exist_ok=True)

CACHE = RESULTS / ".cache"
CACHE.mkdir(parents=True, exist_ok=True)


def _load_cache(name):
    p = CACHE / f"{name}.json"
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return None


def _save_cache(name, data):
    with open(CACHE / f"{name}.json", "w") as f:
        json.dump(data, f)


def parse_pdb(pdb_path):
    """
    Parse PDB file, extract protein atoms and ligand atoms.
    Returns dict with 'protein' and 'ligand' atom lists.
    """
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    record = line[:6].strip()
                    serial = int(line[6:11].strip()) if line[6:11].strip() else 0
                    name = line[12:16].strip()
                    alt_loc = line[16].strip()
                    res_name = line[17:20].strip()
                    chain = line[21].strip()
                    res_seq = int(line[22:26].strip()) if line[22:26].strip() else 0
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    occupancy = float(line[54:60].strip()) if line[54:60].strip() else 1.0
                    bfactor = float(line[60:66].strip()) if line[60:66].strip() else 0.0
                    element = line[76:78].strip() if len(line) > 76 else ""
                    atoms.append({
                        "record": record,
                        "serial": serial,
                        "name": name,
                        "alt_loc": alt_loc,
                        "res_name": res_name,
                        "chain": chain,
                        "res_seq": res_seq,
                        "x": x, "y": y, "z": z,
                        "occupancy": occupancy,
                        "bfactor": bfactor,
                        "element": element,
                        "is_ligand": record == "HETATM" and res_name not in
                            ("HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "ACT",
                             "CL", "NA", "MG", "CA", "ZN", "MN", "FE"),
                    })
                except (ValueError, IndexError):
                    continue
    return atoms


def get_ca_atoms(atoms):
    """Extract C-alpha atoms (protein backbone)."""
    return [a for a in atoms if a["name"] == "CA" and not a["is_ligand"]]


def get_ligand_atoms(atoms):
    """Extract ligand atoms (non-solvent HETATM)."""
    return [a for a in atoms if a["is_ligand"]]


def compute_contacts(protein_atoms, ligand_atoms, cutoff=5.0):
    """Compute contact counts between protein and ligand."""
    contacts = 0
    hbonds = 0
    hydrophobic = 0
    polar = 0

    # Define polar and hydrophobic elements
    polar_elements = {"N", "O", "S"}
    hydrophobic_elements = {"C"}

    for la in ligand_atoms:
        for pa in protein_atoms:
            dx = la["x"] - pa["x"]
            dy = la["y"] - pa["y"]
            dz = la["z"] - pa["z"]
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if dist <= cutoff:
                contacts += 1
                le = la.get("element", "") or la["name"][0]
                pe = pa.get("element", "") or pa["name"][0]
                if le in polar_elements and pe in polar_elements and dist <= 3.5:
                    hbonds += 1
                if le in hydrophobic_elements and pe in hydrophobic_elements:
                    hydrophobic += 1
                if le in polar_elements or pe in polar_elements:
                    polar += 1

    return {
        "contact_count": contacts,
        "hbond_count": hbonds,
        "hydrophobic_contacts": hydrophobic,
        "polar_contacts": polar,
    }


def compute_min_distances(protein_atoms, ligand_atoms):
    """Compute minimum distances between ligand and protein residue groups."""
    # Group protein atoms by residue
    res_groups = {}
    for a in protein_atoms:
        key = (a["chain"], a["res_seq"])
        if key not in res_groups:
            res_groups[key] = []
        res_groups[key].append(a)

    min_dists = []
    for la in ligand_atoms:
        for key, atoms_list in res_groups.items():
            for pa in atoms_list:
                dx = la["x"] - pa["x"]
                dy = la["y"] - pa["y"]
                dz = la["z"] - pa["z"]
                dist = math.sqrt(dx * dx + dy * dy + dz * dz)
                min_dists.append(dist)

    if not min_dists:
        return {
            "min_distance": 999.0,
            "mean_contact_distance": 999.0,
            "n_residues_in_contact": 0,
        }

    # Count residues in contact (< 5 Angstrom)
    res_dists = {}
    for la in ligand_atoms:
        for pa in protein_atoms:
            key = (pa["chain"], pa["res_seq"])
            dx = la["x"] - pa["x"]
            dy = la["y"] - pa["y"]
            dz = la["z"] - pa["z"]
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if key not in res_dists or dist < res_dists[key]:
                res_dists[key] = dist

    contact_residues = sum(1 for d in res_dists.values() if d < 5.0)
    contact_dists = [d for d in res_dists.values() if d < 5.0]

    return {
        "min_distance": float(min(min_dists)) if min_dists else 999.0,
        "mean_contact_distance": float(np.mean(contact_dists)) if contact_dists else 999.0,
        "n_residues_in_contact": contact_residues,
    }


def compute_backbone_coords(atoms):
    """Extract backbone N, CA, C, O coordinates for RMSD."""
    backbone = []
    for a in atoms:
        if a["name"] in ("N", "CA", "C", "O") and not a["is_ligand"]:
            backbone.append({
                "name": a["name"],
                "chain": a["chain"],
                "res_seq": a["res_seq"],
                "x": a["x"], "y": a["y"], "z": a["z"],
            })
    return backbone


def compute_bfactors(atoms):
    """Compute average B-factor statistics for protein atoms."""
    protein_atoms = [a for a in atoms if not a["is_ligand"]]
    if not protein_atoms:
        return {"mean_bfactor": 0, "std_bfactor": 0, "max_bfactor": 0}
    bfactors = [a["bfactor"] for a in protein_atoms]
    return {
        "mean_bfactor": float(np.mean(bfactors)),
        "std_bfactor": float(np.std(bfactors)),
        "max_bfactor": float(np.max(bfactors)),
    }


def find_pdb_files():
    """Find all relevant PDB files."""
    pdb_files = []
    # SAR structures
    sar_dir = PROJECT / "phase08_boltz_affinity" / "boltz2_sar_outputs"
    if sar_dir.exists():
        for p in sar_dir.glob("*_model_*.pdb"):
            pdb_files.append(p)

    # Inactive structures
    inactive_dir = PROJECT / "phase08_boltz_affinity" / "boltz2_inactive_outputs"
    if inactive_dir.exists():
        for p in inactive_dir.glob("*_model_*.pdb"):
            pdb_files.append(p)

    # Reference receptor
    ref = PROJECT / "09_docking" / "receptor_2to3.pdb"
    if ref.exists():
        pdb_files.append(ref)

    return pdb_files


def main():
    print("=" * 70)
    print("STRUCTURAL FEATURES - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/3] Finding PDB files...")
    pdb_files = find_pdb_files()
    print(f"  Found {len(pdb_files)} PDB files")

    print("\n[2/3] Extracting structural features...")
    all_features = {}
    reference_atoms = None
    ref_pdb = PROJECT / "09_docking" / "receptor_2to3.pdb"

    if ref_pdb.exists():
        print(f"  Loading reference: {ref_pdb.name}")
        reference_atoms = parse_pdb(str(ref_pdb))

    processed = 0
    skipped = 0
    for pdb_path in sorted(pdb_files):
        try:
            atoms = parse_pdb(str(pdb_path))
            if not atoms:
                skipped += 1
                continue

            protein_atoms = [a for a in atoms if not a["is_ligand"]]
            ligand_atoms = get_ligand_atoms(atoms)
            ca_atoms = get_ca_atoms(atoms)

            # Extract compound name from path
            name = pdb_path.stem
            # Try to extract compound ID from path
            parts = str(pdb_path).split("/")
            compound_name = name

            contacts = compute_contacts(protein_atoms, ligand_atoms, cutoff=5.0)
            min_dists = compute_min_distances(protein_atoms, ligand_atoms)
            bfactors = compute_bfactors(atoms)

            feat = {
                "pdb_file": str(pdb_path.relative_to(PROJECT)),
                "n_atoms": len(atoms),
                "n_protein_atoms": len(protein_atoms),
                "n_ligand_atoms": len(ligand_atoms),
                "n_ca_atoms": len(ca_atoms),
                **contacts,
                **min_dists,
                **bfactors,
            }

            all_features[name] = feat
            processed += 1
            if processed % 10 == 0:
                print(f"  Processed {processed} files...")

        except Exception as e:
            print(f"  WARNING: Error processing {pdb_path.name}: {e}")
            skipped += 1

    print(f"  Processed: {processed}, Skipped: {skipped}")

    print("\n[3/3] Saving structural features...")
    # Save as JSON
    with open(STRUCTURAL / "structural_features.json", "w") as f:
        json.dump(all_features, f, indent=2)

    # Save as CSV
    if all_features:
        names = sorted(all_features.keys())
        fieldnames = list(all_features[names[0]].keys())
        with open(STRUCTURAL / "structural_features.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for name in names:
                writer.writerow(all_features[name])

    # Save feature matrix as numpy
    if all_features:
        names = sorted(all_features.keys())
        feat_keys = [k for k in all_features[names[0]].keys() if k != "pdb_file"]
        matrix = np.array([[all_features[n][k] for k in feat_keys] for n in names])
        np.save(STRUCTURAL / "structural_features.npy", matrix)
        with open(STRUCTURAL / "feature_names.json", "w") as f:
            json.dump(feat_keys, f)

    print(f"\nResults saved to: {STRUCTURAL}")
    print("  structural_features.json")
    print("  structural_features.csv")
    print("  structural_features.npy")
    print("  feature_names.json")


if __name__ == "__main__":
    main()
