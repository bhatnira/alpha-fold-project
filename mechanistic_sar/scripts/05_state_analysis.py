#!/usr/bin/env python3
"""
Step 5: State Analysis - Compare ligand-bound structures against apo reference
to compute state-like signatures: backbone RMSD, pocket RMSD, residue displacements.
Save to features/structural/.
"""
import csv, json, os, math
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


def parse_pdb_atoms(pdb_path):
    """Parse PDB file and return list of atom dicts."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    record = line[:6].strip()
                    name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain = line[21].strip()
                    res_seq = int(line[22:26].strip()) if line[22:26].strip() else 0
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    bfactor = float(line[60:66].strip()) if line[60:66].strip() else 0.0
                    atoms.append({
                        "record": record,
                        "name": name,
                        "res_name": res_name,
                        "chain": chain,
                        "res_seq": res_seq,
                        "x": x, "y": y, "z": z,
                        "bfactor": bfactor,
                        "is_ligand": record == "HETATM" and res_name not in
                            ("HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "ACT"),
                    })
                except (ValueError, IndexError):
                    continue
    return atoms


def get_ca_coords(atoms):
    """Get C-alpha coordinates as array, keyed by (chain, res_seq)."""
    ca = {}
    for a in atoms:
        if a["name"] == "CA" and not a["is_ligand"]:
            key = (a["chain"], a["res_seq"])
            ca[key] = np.array([a["x"], a["y"], a["z"]])
    return ca


def get_backbone_coords(atoms):
    """Get backbone atom coordinates grouped by residue."""
    backbone = {}
    for a in atoms:
        if a["name"] in ("N", "CA", "C", "O") and not a["is_ligand"]:
            key = (a["chain"], a["res_seq"], a["name"])
            backbone[key] = np.array([a["x"], a["y"], a["z"]])
    return backbone


def reference_alignment(mobile, target):
    """
    Simple Kabsch-like alignment: find rotation + translation to align
    mobile onto target using matched atom pairs.
    """
    if len(mobile) < 3:
        return None, None

    # Center both sets
    centroid_m = np.mean(mobile, axis=0)
    centroid_t = np.mean(target, axis=0)
    mobile_c = mobile - centroid_m
    target_c = target - centroid_t

    # SVD for optimal rotation
    H = mobile_c.T @ target_c
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # Ensure proper rotation (det = +1)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    t = centroid_t - R @ centroid_m
    return R, t


def apply_transform(coords, R, t):
    """Apply rotation R and translation t to coordinates."""
    return (R @ coords.T).T + t


def compute_rmsd(coords1, coords2):
    """Compute RMSD between two coordinate sets."""
    if len(coords1) != len(coords2) or len(coords1) == 0:
        return float("inf")
    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def compute_residue_displacements(ref_ca, bound_ca):
    """Compute per-residue C-alpha displacements."""
    displacements = {}
    for key in ref_ca:
        if key in bound_ca:
            d = np.linalg.norm(ref_ca[key] - bound_ca[key])
            displacements[f"{key[0]}_{key[1]}"] = float(d)
    return displacements


def find_pocket_residues(atoms, cutoff=8.0):
    """
    Identify pocket residues: residues within cutoff of any ligand atom.
    """
    ligand_atoms = [a for a in atoms if a["is_ligand"]]
    if not ligand_atoms:
        return set()

    pocket_residues = set()
    for a in atoms:
        if a["is_ligand"]:
            continue
        key = (a["chain"], a["res_seq"])
        for la in ligand_atoms:
            dx = a["x"] - la["x"]
            dy = a["y"] - la["y"]
            dz = a["z"] - la["z"]
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if dist <= cutoff:
                pocket_residues.add(key)
                break
    return pocket_residues


def find_reference_pdb():
    """Find the apo reference PDB."""
    ref = PROJECT / "09_docking" / "receptor_2to3.pdb"
    if ref.exists():
        return ref
    return None


def find_bound_pdbs():
    """Find all bound (ligand-bound) PDB files."""
    pdbs = []
    for root_dir in [
        PROJECT / "phase08_boltz_affinity" / "boltz2_sar_outputs",
        PROJECT / "phase08_boltz_affinity" / "boltz2_inactive_outputs",
    ]:
        if root_dir.exists():
            for p in root_dir.glob("*_model_*.pdb"):
                pdbs.append(p)
    return pdbs


def main():
    print("=" * 70)
    print("STATE ANALYSIS - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/4] Finding reference structure...")
    ref_pdb = find_reference_pdb()
    if ref_pdb is None:
        print("  ERROR: No reference PDB found. Exiting.")
        return
    print(f"  Reference: {ref_pdb.name}")

    ref_atoms = parse_pdb_atoms(str(ref_pdb))
    ref_ca = get_ca_coords(ref_atoms)
    ref_backbone = get_backbone_coords(ref_atoms)
    print(f"  Reference C-alpha atoms: {len(ref_ca)}")

    print("\n[2/4] Finding bound structures...")
    bound_pdbs = find_bound_pdbs()
    print(f"  Found {len(bound_pdbs)} bound structures")

    print("\n[3/4] Computing state signatures...")
    state_features = {}
    processed = 0

    for pdb_path in sorted(bound_pdbs):
        try:
            bound_atoms = parse_pdb_atoms(str(pdb_path))
            if not bound_atoms:
                continue

            bound_ca = get_ca_coords(bound_atoms)
            bound_backbone = get_backbone_coords(bound_atoms)

            # Match C-alpha atoms
            common_keys = sorted(set(ref_ca.keys()) & set(bound_ca.keys()))
            if len(common_keys) < 10:
                continue

            ref_matched = np.array([ref_ca[k] for k in common_keys])
            bound_matched = np.array([bound_ca[k] for k in common_keys])

            # Superimpose
            R, t = reference_alignment(bound_matched, ref_matched)
            if R is None:
                continue

            # Apply transformation
            aligned = apply_transform(bound_matched, R, t)

            # Compute RMSD
            ca_rmsd = compute_rmsd(aligned, ref_matched)

            # Backbone RMSD (matched N, CA, C, O)
            bb_keys = sorted(set(ref_backbone.keys()) & set(bound_backbone.keys()))
            if len(bb_keys) > 0:
                ref_bb = np.array([ref_backbone[k] for k in bb_keys])
                bound_bb = np.array([bound_backbone[k] for k in bb_keys])
                aligned_bb = apply_transform(bound_bb, R, t)
                backbone_rmsd = compute_rmsd(aligned_bb, ref_bb)
            else:
                backbone_rmsd = ca_rmsd

            # Pocket RMSD
            pocket_residues = find_pocket_residues(bound_atoms, cutoff=8.0)
            pocket_keys = [k for k in common_keys if k in pocket_residues]
            if len(pocket_keys) > 0:
                ref_pocket = np.array([ref_ca[k] for k in pocket_keys])
                bound_pocket = np.array([bound_ca[k] for k in pocket_keys])
                aligned_pocket = apply_transform(bound_pocket, R, t)
                pocket_rmsd = compute_rmsd(aligned_pocket, ref_pocket)
            else:
                pocket_rmsd = ca_rmsd

            # Per-residue displacements
            displacements = compute_residue_displacements(ref_ca, bound_ca)

            # Compute state scores (simple heuristic)
            # Lower RMSD -> more "closed" (apo-like), higher -> more "open"
            open_score = min(ca_rmsd / 5.0, 1.0)
            closed_score = max(0, 1.0 - ca_rmsd / 5.0)
            desensitized_score = pocket_rmsd / max(ca_rmsd, 0.01)

            # Mean displacement
            disp_values = list(displacements.values()) if displacements else [0]
            mean_disp = float(np.mean(disp_values))
            max_disp = float(np.max(disp_values))

            name = pdb_path.stem
            state_features[name] = {
                "pdb_file": str(pdb_path.relative_to(PROJECT)),
                "backbone_rmsd": backbone_rmsd,
                "ca_rmsd": ca_rmsd,
                "pocket_rmsd": pocket_rmsd,
                "n_residues_compared": len(common_keys),
                "n_pocket_residues": len(pocket_residues),
                "mean_displacement": mean_disp,
                "max_displacement": max_disp,
                "open_state_score": open_score,
                "closed_state_score": closed_score,
                "desensitized_score": desensitized_score,
                "state_preference": "open" if open_score > 0.6 else ("closed" if closed_score > 0.6 else "intermediate"),
            }

            processed += 1
            if processed % 10 == 0:
                print(f"  Processed {processed} structures...")

        except Exception as e:
            print(f"  WARNING: Error processing {pdb_path.name}: {e}")

    print(f"  Processed: {processed} structures")

    print("\n[4/4] Saving state analysis results...")
    # Save JSON
    with open(STRUCTURAL / "state_analysis.json", "w") as f:
        json.dump(state_features, f, indent=2)

    # Save CSV
    if state_features:
        names = sorted(state_features.keys())
        fieldnames = list(state_features[names[0]].keys())
        with open(STRUCTURAL / "state_analysis.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for name in names:
                writer.writerow(state_features[name])

    # Summary
    if state_features:
        rmsds = [v["ca_rmsd"] for v in state_features.values()]
        print(f"\n  CA-RMSD range: [{min(rmsds):.3f}, {max(rmsds):.3f}] Angstrom")
        print(f"  Mean CA-RMSD: {np.mean(rmsds):.3f} Angstrom")

        state_counts = {}
        for v in state_features.values():
            s = v["state_preference"]
            state_counts[s] = state_counts.get(s, 0) + 1
        print(f"  State distribution: {state_counts}")

    print(f"\nResults saved to: {STRUCTURAL}")
    print("  state_analysis.json")
    print("  state_analysis.csv")


if __name__ == "__main__":
    main()
