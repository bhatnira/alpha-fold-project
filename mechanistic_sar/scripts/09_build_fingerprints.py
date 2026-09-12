#!/usr/bin/env python3
"""
Step 9: Build Fingerprints - Build Morgan/ECFP-like fingerprints from SMILES
using basic string operations (no rdkit), compute simple physicochemical
descriptors (MW, logP, HBD, HBA from SMILES patterns). Save to features/chemical/.
"""
import csv, json, os, re, math
from pathlib import Path
import numpy as np

PROJECT = Path("/cluster/home/nbhatt04/lean_pipeline")
MECHANISTIC = PROJECT / "mechanistic_sar"
RESULTS = MECHANISTIC / "results"
FEATURES = MECHANISTIC / "features"
CHEMICAL = FEATURES / "chemical"
CHEMICAL.mkdir(parents=True, exist_ok=True)

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


# Atomic weights for MW calculation
ATOMIC_WEIGHTS = {
    "C": 12.011, "c": 12.011,
    "N": 14.007, "n": 14.007,
    "O": 15.999, "o": 15.999,
    "S": 32.065, "s": 32.065,
    "P": 30.974, "p": 30.974,
    "F": 18.998,
    "Cl": 35.453,
    "Br": 79.904,
    "I": 126.904,
    "H": 1.008,
    "B": 10.811,
    "Si": 28.086,
}


def parse_smiles_atoms(smiles):
    """
    Parse SMILES string into list of atom tokens.
    Simple tokenizer handling brackets and common atoms.
    """
    atoms = []
    i = 0
    while i < len(smiles):
        ch = smiles[i]
        if ch in "0123456789" or ch in "=-#/\\+.":
            i += 1
            continue
        if ch in "()[]":
            i += 1
            continue
        if ch == "%" and i + 2 < len(smiles):
            # Two-digit ring closure
            i += 3
            continue
        # Check for two-letter atoms
        if i + 1 < len(smiles):
            two = ch + smiles[i + 1]
            if two in ("Cl", "Br", "Si", "Se", "Te"):
                atoms.append(two)
                i += 2
                continue
        if ch in ATOMIC_WEIGHTS:
            atoms.append(ch)
        i += 1
    return atoms


def compute_molecular_weight(smiles):
    """Compute molecular weight from SMILES."""
    atoms = parse_smiles_atoms(smiles)
    mw = 0.0
    for atom in atoms:
        mw += ATOMIC_WEIGHTS.get(atom, 0.0)
    return mw


def count_hbd(smiles):
    """Count hydrogen bond donors (NH, OH groups)."""
    # Simple pattern matching
    count = 0
    # Count [NH] and [OH] patterns (explicit H in brackets)
    count += len(re.findall(r'\[N[Hh]\]', smiles))
    count += len(re.findall(r'\[O[Hh]\]', smiles))
    count += smiles.count("N") - smiles.count("[N]")
    count += smiles.count("O") - smiles.count("[O]")
    # Rough approximation: NH and OH in rings
    count += smiles.count("c1") * 0  # aromatic N doesn't always have H
    return max(0, min(count, 20))


def count_hba(smiles):
    """Count hydrogen bond acceptors (N, O atoms)."""
    count = 0
    for ch in smiles:
        if ch in "No":
            count += 1
    # Count bracketed atoms
    count += len(re.findall(r'\[N\]', smiles))
    count += len(re.findall(r'\[O\]', smiles))
    count += len(re.findall(r'\[n\]', smiles))
    count += len(re.findall(r'\[o\]', smiles))
    return count


def compute_logp_estimate(smiles):
    """
    Estimate logP using a simple fragment-based approach.
    Based on Wildman-Crippen fragment contributions.
    """
    atoms = parse_smiles_atoms(smiles)
    n_carbon = sum(1 for a in atoms if a in ("C", "c"))
    n_nitrogen = sum(1 for a in atoms if a in ("N", "n"))
    n_oxygen = sum(1 for a in atoms if a in ("O", "o"))
    n_sulfur = sum(1 for a in atoms if a in ("S", "s"))
    n_halogen = sum(1 for a in atoms if a in ("F", "Cl", "Br", "I"))

    # Simple logP estimation
    logp = 0.0
    logp += n_carbon * 0.20
    logp -= n_nitrogen * 0.70
    logp -= n_oxygen * 1.00
    logp += n_sulfur * 0.30
    logp += n_halogen * 0.50

    # Ring correction
    n_rings = smiles.count("1") // 2 + smiles.count("c1") // 2
    logp += n_rings * 0.10

    return round(logp, 3)


def count_rotatable_bonds(smiles):
    """Count rotatable bonds (simplified)."""
    count = 0
    # Single bonds not in rings
    count += smiles.count("-") - smiles.count("=") - smiles.count("#")
    # Account for ring closures
    ring_closures = len(re.findall(r'[0-9]', smiles)) // 2
    count = max(0, count - ring_closures)
    return min(count, 30)


def count_rings(smiles):
    """Count number of rings."""
    # Count aromatic rings (lowercase letters)
    aromatic = 0
    i = 0
    while i < len(smiles):
        if smiles[i] in "cnos":
            aromatic += 1
        i += 1
    aromatic = aromatic // 5  # rough: 5-6 atoms per aromatic ring

    # Count aliphatic rings (digits)
    ring_closures = len(re.findall(r'[0-9]', smiles))
    aliphatic = ring_closures // 2

    return aromatic + aliphatic


def compute_tpsa(smiles):
    """Estimate topological polar surface area."""
    # Simple approximation: sum of fragment contributions
    n_oxygen = sum(1 for ch in smiles if ch in "Oo")
    n_nitrogen = sum(1 for ch in smiles if ch in "Nn")
    n_sulfur = sum(1 for ch in smiles if ch in "Ss")

    tpsa = 0.0
    tpsa += n_oxygen * 20.0
    tpsa += n_nitrogen * 26.0
    tpsa += n_sulfur * 32.0

    return round(tpsa, 2)


def smiles_to_fingerprint(smiles, n_bits=2048, radius=2):
    """
    Build a Morgan/ECFP-like fingerprint using string hashing.
    No rdkit - use simple substring hashing.
    """
    fp = np.zeros(n_bits, dtype=int)

    if not smiles:
        return fp

    # Generate substrings of varying lengths (mimicking Morgan neighborhoods)
    substrings = set()
    # Add individual atoms
    for ch in smiles:
        if ch in "CNOSPFIBcnosClBr":
            substrings.add(ch)

    # Add 2-grams
    for i in range(len(smiles) - 1):
        sub = smiles[i:i + 2]
        if not all(c in "0123456789=-#/\\+().[]" for c in sub):
            substrings.add(sub)

    # Add 3-grams
    for i in range(len(smiles) - 2):
        sub = smiles[i:i + 3]
        if not all(c in "0123456789=-#/\\+().[]" for c in sub):
            substrings.add(sub)

    # Add 4-grams
    for i in range(len(smiles) - 3):
        sub = smiles[i:i + 4]
        if not all(c in "0123456789=-#/\\+().[]" for c in sub):
            substrings.add(sub)

    # Hash each substring to a bit position
    for sub in substrings:
        h = hash(sub) % n_bits
        fp[h] = 1

    # Add some ring-specific patterns
    ring_patterns = re.findall(r'[a-z]{3,}', smiles)
    for pat in ring_patterns:
        h = hash(pat) % n_bits
        fp[h] = 1

    # Bracket patterns (special atoms)
    bracket_atoms = re.findall(r'\[.*?\]', smiles)
    for ba in bracket_atoms:
        h = hash(ba) % n_bits
        fp[h] = 1

    return fp


def load_smiles_data():
    """Load SMILES from experimental manifest."""
    manifest = RESULTS / "data_audit" / "experimental_manifest.csv"
    if not manifest.exists():
        return {}

    data = {}
    with open(manifest) as f:
        for row in csv.DictReader(f):
            cid = row.get("compound_id", "")
            if cid:
                data[cid] = {
                    "compound_id": cid,
                    "smiles": row.get("smiles", ""),
                    "is_active": row.get("is_active", "False") == "True",
                }
    return data


def main():
    print("=" * 70)
    print("CHEMICAL FINGERPRINTS - Mechanistic SAR Analysis")
    print("=" * 70)

    print("\n[1/4] Loading SMILES data...")
    smiles_data = load_smiles_data()
    print(f"  Compounds: {len(smiles_data)}")

    print("\n[2/4] Computing fingerprints and descriptors...")
    features = {}

    for cid, sdata in sorted(smiles_data.items()):
        smiles = sdata["smiles"]
        if not smiles:
            print(f"  WARNING: No SMILES for compound {cid}")
            continue

        # Compute fingerprint
        fp = smiles_to_fingerprint(smiles, n_bits=2048, radius=2)

        # Compute descriptors
        mw = compute_molecular_weight(smiles)
        logp = compute_logp_estimate(smiles)
        hbd = count_hbd(smiles)
        hba = count_hba(smiles)
        n_rot = count_rotatable_bonds(smiles)
        n_rings = count_rings(smiles)
        tpsa = compute_tpsa(smiles)

        # Derived descriptors
        qed_approx = max(0, min(1, 1.0 - abs(mw - 300) / 500))  # rough drug-likeness

        features[cid] = {
            "compound_id": cid,
            "smiles": smiles,
            "is_active": sdata["is_active"],
            "molecular_weight": mw,
            "logP": logp,
            "hbd": hbd,
            "hba": hba,
            "n_rotatable_bonds": n_rot,
            "n_rings": n_rings,
            "tpsa": tpsa,
            "qed_approx": qed_approx,
            "fingerprint": fp.tolist(),
            "fp_density": float(np.mean(fp)),
        }

    print(f"  Processed: {len(features)} compounds")

    if features:
        # Show some statistics
        mws = [f["molecular_weight"] for f in features.values()]
        logps = [f["logP"] for f in features.values()]
        print(f"  MW range: [{min(mws):.1f}, {max(mws):.1f}]")
        print(f"  logP range: [{min(logps):.2f}, {max(logps):.2f}]")

    print("\n[3/4] Saving chemical features...")
    # Save JSON (without fingerprint arrays for readability)
    features_json = {}
    for cid, feat in features.items():
        features_json[cid] = {k: v for k, v in feat.items() if k != "fingerprint"}
    with open(CHEMICAL / "chemical_descriptors.json", "w") as f:
        json.dump(features_json, f, indent=2)

    # Save CSV
    if features:
        cids = sorted(features.keys())
        csv_keys = [k for k in features[cids[0]].keys() if k != "fingerprint"]
        with open(CHEMICAL / "chemical_descriptors.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=csv_keys)
            writer.writeheader()
            for cid in cids:
                writer.writerow({k: features[cid][k] for k in csv_keys})

    # Save fingerprints as numpy
    if features:
        cids = sorted(features.keys())
        fp_matrix = np.array([features[cid]["fingerprint"] for cid in cids])
        np.save(CHEMICAL / "fingerprints.npy", fp_matrix)

        # Save compound order
        with open(CHEMICAL / "fingerprint_compounds.json", "w") as f:
            json.dump(cids, f)

    print("\n[4/4] Summary...")
    if features:
        print(f"  Fingerprint shape: {fp_matrix.shape}")
        print(f"  Bits set per compound: [{np.min(np.sum(fp_matrix, axis=1))}, "
              f"{np.max(np.sum(fp_matrix, axis=1))}]")

    print(f"\nResults saved to: {CHEMICAL}")
    print("  chemical_descriptors.json")
    print("  chemical_descriptors.csv")
    print("  fingerprints.npy")
    print("  fingerprint_compounds.json")


if __name__ == "__main__":
    main()
