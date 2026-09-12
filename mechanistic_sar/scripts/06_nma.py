#!/usr/bin/env python3
"""
Step 6: NMA (Normal Mode Analysis) - Implement simple Elastic Network Model
using numpy: build contact matrix from C-alpha atoms, compute Hessian,
extract lowest modes, compute fluctuations. Save to features/nma/.
"""
import csv, json, os, math
from pathlib import Path
import numpy as np

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
FEATURES = MECHANISTIC / "features"
NMA = FEATURES / "nma"
NMA.mkdir(parents=True, exist_ok=True)

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


def parse_ca_atoms(pdb_path):
    """Extract C-alpha atom coordinates from PDB."""
    ca_atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and "CA" in line[12:16]:
                try:
                    chain = line[21].strip()
                    res_seq = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    ca_atoms.append({
                        "chain": chain,
                        "res_seq": res_seq,
                        "x": x, "y": y, "z": z,
                    })
                except (ValueError, IndexError):
                    continue
    return ca_atoms


def build_distance_matrix(coords):
    """Build pairwise distance matrix."""
    n = len(coords)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(coords[i] - coords[j])
            dist[i, j] = d
            dist[j, i] = d
    return dist


def build_contact_matrix(dist_matrix, cutoff=8.0):
    """Build binary contact matrix: 1 if distance < cutoff."""
    return (dist_matrix < cutoff).astype(float)


def build_hessian(coords, contact_matrix, cutoff=8.0, spring_constant=1.0):
    """
    Build the 3N x 3N Hessian matrix for the Elastic Network Model.
    H_ij = -k * (r_i - r_j)(r_i - r_j)^T / |r_i - r_j|^2 for contacting pairs
    H_ii = sum over j of -H_ij
    """
    n = len(coords)
    H = np.zeros((3 * n, 3 * n))

    for i in range(n):
        for j in range(i + 1, n):
            if contact_matrix[i, j] == 0:
                continue
            dr = coords[i] - coords[j]
            dist_sq = np.sum(dr * dr)
            if dist_sq < 1e-10:
                continue
            # Outer product
            outer = np.outer(dr, dr)
            h_block = -spring_constant * outer / dist_sq

            # Off-diagonal blocks
            H[3 * i:3 * i + 3, 3 * j:3 * j + 3] = h_block
            H[3 * j:3 * j + 3, 3 * i:3 * i + 3] = h_block

            # Diagonal blocks
            H[3 * i:3 * i + 3, 3 * i:3 * i + 3] -= h_block
            H[3 * j:3 * j + 3, 3 * j:3 * j + 3] -= h_block

    return H


def extract_modes(H, n_modes=10):
    """
    Extract lowest non-trivial modes from Hessian.
    Uses eigenvalue decomposition.
    """
    # Eigenvalue decomposition
    eigenvalues, eigenvectors = np.linalg.eigh(H)

    # Remove trivial modes (lowest 6 are translational/rotational)
    # The first 6 eigenvalues should be ~0
    non_trivial_start = 6
    if len(eigenvalues) <= non_trivial_start:
        return eigenvalues, eigenvectors

    evals = eigenvalues[non_trivial_start:non_trivial_start + n_modes]
    evecs = eigenvectors[:, non_trivial_start:non_trivial_start + n_modes]

    return evals, evecs


def compute_fluctuations(eigenvalues, eigenvectors, n_atoms):
    """
    Compute mean-square fluctuations from modes.
    B_i = sum_m (1/lambda_m) * |e_m(i)|^2
    """
    n_modes = len(eigenvalues)
    fluctuations = np.zeros(n_atoms)

    for m in range(n_modes):
        if eigenvalues[m] < 1e-10:
            continue
        mode = eigenvectors[:, m]
        # Reshape to (n_atoms, 3)
        mode_reshaped = mode.reshape(n_atoms, 3)
        # Sum squares over x, y, z for each atom
        sq = np.sum(mode_reshaped * mode_reshaped, axis=1)
        fluctuations += sq / eigenvalues[m]

    return fluctuations


def compute_collectivity(eigenvector, n_atoms):
    """
    Compute mode collectivity: measure of how many atoms participate.
    kappa = 1/N * exp(-sum_i p_i * ln(p_i))
    where p_i = |u_i|^2 / sum_j |u_j|^2
    """
    mode = eigenvector.reshape(n_atoms, 3)
    sq = np.sum(mode * mode, axis=1)
    total = np.sum(sq)
    if total < 1e-20:
        return 0.0

    p = sq / total
    p = p[p > 1e-20]
    entropy = -np.sum(p * np.log(p))
    return float(np.exp(entropy) / n_atoms)


def compute_cross_correlation(eigenvectors, eigenvalues, n_atoms):
    """
    Compute dynamic cross-correlation matrix.
    C_ij = sum_m (1/lambda_m) * u_m(i) . u_m(j)
    """
    C = np.zeros((n_atoms, n_atoms))
    n_modes = len(eigenvalues)

    for m in range(n_modes):
        if eigenvalues[m] < 1e-10:
            continue
        mode = eigenvectors[:, m].reshape(n_atoms, 3)
        for i in range(n_atoms):
            for j in range(i, n_atoms):
                dot = np.dot(mode[i], mode[j])
                C[i, j] += dot / eigenvalues[m]
                C[j, i] = C[i, j]

    # Normalize
    diag = np.diag(C)
    diag = np.where(diag > 0, diag, 1.0)
    D = np.sqrt(np.outer(diag, diag))
    C_norm = C / D
    np.fill_diagonal(C_norm, 1.0)

    return C_norm


def main():
    print("=" * 70)
    print("NORMAL MODE ANALYSIS (ENM) - Mechanistic SAR Analysis")
    print("=" * 70)

    # Find PDB files
    pdb_files = []
    for root_dir in [
        PROJECT / "phase08_boltz_affinity" / "boltz2_sar_outputs",
        PROJECT / "phase08_boltz_affinity" / "boltz2_inactive_outputs",
    ]:
        if root_dir.exists():
            for p in root_dir.glob("*_model_*.pdb"):
                pdb_files.append(p)

    ref = PROJECT / "09_docking" / "receptor_2to3.pdb"
    if ref.exists():
        pdb_files.insert(0, ref)

    print(f"\n[1/4] Found {len(pdb_files)} PDB files")

    print("\n[2/4] Computing ENM for each structure...")
    nma_features = {}
    processed = 0

    for pdb_path in sorted(pdb_files):
        try:
            ca_atoms = parse_ca_atoms(str(pdb_path))
            if len(ca_atoms) < 10:
                print(f"  WARNING: {pdb_path.name}: too few C-alpha atoms ({len(ca_atoms)})")
                continue

            coords = np.array([[a["x"], a["y"], a["z"]] for a in ca_atoms])
            n_atoms = len(coords)

            # Build matrices
            dist_matrix = build_distance_matrix(coords)
            contact_matrix = build_contact_matrix(dist_matrix, cutoff=8.0)
            n_contacts = int(np.sum(contact_matrix) / 2)

            # Build and diagonalize Hessian
            H = build_hessian(coords, contact_matrix)
            evals, evecs = extract_modes(H, n_modes=min(20, 3 * n_atoms - 6))

            # Compute properties
            fluctuations = compute_fluctuations(evals, evecs, n_atoms)
            mean_fluctuation = float(np.mean(fluctuations))

            # Collectivity of first few modes
            collectivities = []
            for m in range(min(5, len(evals))):
                c = compute_collectivity(evecs[:, m], n_atoms)
                collectivities.append(c)

            # Cross-correlation (only for small systems)
            if n_atoms <= 200:
                cc = compute_cross_correlation(evecs, evals, n_atoms)
                # Average correlation in pocket region (first 30 residues)
                pocket_cc = float(np.mean(cc[:min(30, n_atoms), :min(30, n_atoms)]))
            else:
                pocket_cc = 0.0

            name = pdb_path.stem
            nma_features[name] = {
                "pdb_file": str(pdb_path.relative_to(PROJECT)),
                "n_atoms": n_atoms,
                "n_contacts": n_contacts,
                "lowest_eigenvalue": float(evals[0]) if len(evals) > 0 else 0.0,
                "second_eigenvalue": float(evals[1]) if len(evals) > 1 else 0.0,
                "mean_collectivity": float(np.mean(collectivities)) if collectivities else 0.0,
                "max_collectivity": float(np.max(collectivities)) if collectivities else 0.0,
                "mean_fluctuation": mean_fluctuation,
                "max_fluctuation": float(np.max(fluctuations)),
                "pocket_cross_correlation": pocket_cc,
                "eigenvalues": [float(e) for e in evals[:10]],
                "fluctuations_per_residue": {
                    f"{ca_atoms[i]['chain']}_{ca_atoms[i]['res_seq']}": float(fluctuations[i])
                    for i in range(n_atoms)
                },
            }

            processed += 1
            print(f"  Processed {pdb_path.name}: {n_atoms} atoms, "
                  f"lambda1={evals[0]:.4f}")

        except Exception as e:
            print(f"  WARNING: Error processing {pdb_path.name}: {e}")

    print(f"\n  Processed: {processed} structures")

    print("\n[3/4] Saving NMA features...")
    # Save JSON
    with open(NMA / "nma_features.json", "w") as f:
        json.dump(nma_features, f, indent=2)

    # Save CSV (without per-residue data)
    if nma_features:
        names = sorted(nma_features.keys())
        csv_keys = [k for k in nma_features[names[0]].keys()
                    if k not in ("eigenvalues", "fluctuations_per_residue")]
        with open(NMA / "nma_features.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=csv_keys)
            writer.writeheader()
            for name in names:
                writer.writerow({k: nma_features[name][k] for k in csv_keys})

    # Save eigenvalues comparison
    if nma_features:
        names = sorted(nma_features.keys())
        max_modes = max(len(nma_features[n].get("eigenvalues", [])) for n in names)
        eigenval_matrix = np.zeros((len(names), max_modes))
        for i, name in enumerate(names):
            evals = nma_features[name].get("eigenvalues", [])
            eigenval_matrix[i, :len(evals)] = evals
        np.save(NMA / "eigenvalues.npy", eigenval_matrix)

    print("\n[4/4] Summary...")
    if nma_features:
        lambdas = [v["lowest_eigenvalue"] for v in nma_features.values()]
        print(f"  Lowest eigenvalue range: [{min(lambdas):.6f}, {max(lambdas):.6f}]")
        fluc = [v["mean_fluctuation"] for v in nma_features.values()]
        print(f"  Mean fluctuation range: [{min(fluc):.4f}, {max(fluc):.4f}]")

    print(f"\nResults saved to: {NMA}")
    print("  nma_features.json")
    print("  nma_features.csv")
    print("  eigenvalues.npy")


if __name__ == "__main__":
    main()
