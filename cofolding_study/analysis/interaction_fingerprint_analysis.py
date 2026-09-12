#!/usr/bin/env python
"""
Interaction Fingerprint Analysis & SAR Correlation

Processes Boltz-2/AF3 results to:
1. Extract residue-level interactions for all 30 compounds
2. Build interaction fingerprints
3. Correlate with experimental potency
4. Test active vs inactive discrimination
5. Identify key pharmacophoric interactions

Input: PDB/CIF models from cofolding predictions
Output: Interaction matrices, SAR statistics, figures
"""

import csv
import json
import os
import numpy as np
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Configuration
PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
COFOLDING_DIR = PIPELINE_DIR / "cofolding_study"
RESULTS_DIR = COFOLDING_DIR / "results"
ANALYSIS_DIR = COFOLDING_DIR / "analysis"
ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

# Load compound data
CSV_PATH = PIPELINE_DIR / "modulator-dataset-a9a10.csv"

# Key residues from literature and pipeline analysis
KEY_RESIDUES = {
    "alpha9": {
        "W176": "H-bond/polar anchor",
        "S175": "H-bond",
        "Y120": "H-bond/polar",
        "Y224": "H-bond/polar",
        "Y217": "H-bond/polar",
        "W55": "Aromatic",
        "C192": "Disulfide",
        "Y93": "Aromatic",
        "Y197": "Aromatic",
    },
    "alpha10": {
        "D145": "Electrostatic",
        "R83": "Electrostatic",
        "W81": "H-bond/polar",
        "R143": "Electrostatic",
        "H153": "Unique to alpha10",
    }
}

# Interaction types
INTERACTION_TYPES = [
    "H_bond",
    "electrostatic",
    "hydrophobic",
    "aromatic",
    "cation_pi",
    "pi_pi",
    "halogen_mediated",
    "vdw",
    "steric_clash",
]


def load_compounds():
    """Load compound dataset."""
    compounds = []
    with open(CSV_PATH, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            activity = float(row["Activity (uM)"])
            potentiation = float(row["%Potentiation"])
            compounds.append({
                "id": int(row["Identifier"]),
                "smiles": row["Smiles"],
                "activity_uM": activity,
                "potentiation_pct": potentiation,
                "activity_class": classify(activity, potentiation),
                "is_active": activity > 0,
            })
    return compounds


def classify(activity, potentiation):
    if activity == 0 and potentiation == 0:
        return "inactive"
    elif activity < 1:
        return "highly_potent"
    elif activity < 10:
        return "potent"
    elif activity < 2000:
        return "moderate"
    else:
        return "weak"


def parse_pdb_contacts(pdb_path, ligand_chain="L", cutoff=5.0):
    """
    Parse PDB file to extract ligand-residue contacts.
    Returns dict of residue -> interaction list.
    """
    contacts = defaultdict(list)
    
    try:
        with open(pdb_path, "r") as f:
            lines = f.readlines()
    except:
        return contacts
    
    # Extract ligand atoms
    ligand_atoms = []
    receptor_atoms = []
    
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain = line[21].strip()
            res_name = line[17:20].strip()
            res_num = int(line[22:26].strip())
            atom_name = line[12:16].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            
            atom_info = {
                "chain": chain,
                "res_name": res_name,
                "res_num": res_num,
                "atom_name": atom_name,
                "coords": np.array([x, y, z]),
            }
            
            if chain == ligand_chain:
                ligand_atoms.append(atom_info)
            else:
                receptor_atoms.append(atom_info)
    
    # Calculate contacts
    for lig_atom in ligand_atoms:
        for rec_atom in receptor_atoms:
            dist = np.linalg.norm(lig_atom["coords"] - rec_atom["coords"])
            if dist < cutoff:
                key = f"{rec_atom['chain']}:{rec_atom['res_name']}{rec_atom['res_num']}"
                contacts[key].append({
                    "distance": dist,
                    "ligand_atom": lig_atom["atom_name"],
                    "receptor_atom": rec_atom["atom_name"],
                    "interaction_type": classify_interaction(
                        lig_atom, rec_atom, dist
                    ),
                })
    
    return contacts


def classify_interaction(lig_atom, rec_atom, distance):
    """Classify interaction type based on atom types and distance."""
    lig_name = lig_atom["atom_name"].upper()
    rec_name = rec_atom["atom_name"].upper()
    
    # H-bond: donor/acceptor atoms within 3.5 A
    h_bond_donors = {"N", "O", "S"}
    h_bond_acceptors = {"N", "O", "S", "F"}
    
    if distance < 3.5:
        if (lig_name[0] in h_bond_donors and rec_name[0] in h_bond_acceptors) or \
           (lig_name[0] in h_bond_acceptors and rec_name[0] in h_bond_donors):
            return "H_bond"
    
    # Electrostatic: charged groups
    if distance < 4.0:
        if any(x in lig_name for x in ["N+", "O-", "CA"]) or \
           any(x in rec_name for x in ["N+", "O-", "CA", "GLU", "ASP", "ARG", "LYS"]):
            return "electrostatic"
    
    # Hydrophobic: C atoms, non-polar
    if distance < 4.5:
        if lig_name[0] == "C" and rec_name[0] == "C":
            return "hydrophobic"
    
    # Aromatic: ring systems
    if distance < 5.0:
        if any(x in lig_name for x in ["CG", "CD", "CE", "CZ"]) and \
           any(x in rec_name for x in ["CG", "CD", "CE", "CZ"]):
            return "aromatic"
    
    # Van der Waals: close contacts
    if distance < 4.0:
        return "vdw"
    
    return "vdw"


def build_interaction_fingerprint(contacts, key_residues):
    """Build binary interaction fingerprint for key residues."""
    fp = {}
    
    for subunit, residues in key_residues.items():
        for res_name, res_type in residues.items():
            # Find this residue in contacts
            found = False
            min_dist = float("inf")
            interaction_type = None
            
            for key, contact_list in contacts.items():
                if res_name in key:
                    found = True
                    for c in contact_list:
                        if c["distance"] < min_dist:
                            min_dist = c["distance"]
                            interaction_type = c["interaction_type"]
            
            fp[res_name] = {
                "contact": found,
                "min_distance": min_dist if found else None,
                "interaction_type": interaction_type,
                "residue_type": res_type,
                "subunit": subunit,
            }
    
    return fp


def compute_sar_correlation(compounds, fingerprints, site_id):
    """
    Compute correlation between interaction features and experimental potency.
    Uses Spearman correlation for non-parametric analysis.
    """
    from scipy import stats
    
    results = {
        "site_id": site_id,
        "n_compounds": len(compounds),
        "n_active": sum(1 for c in compounds if c["is_active"]),
        "n_inactive": sum(1 for c in compounds if not c["is_active"]),
        "correlations": {},
        "discrimination": {},
    }
    
    # Extract activity values (use log-transform for potency)
    activities = []
    for c in compounds:
        if c["activity_uM"] > 0:
            activities.append(np.log10(c["activity_uM"]))
        else:
            activities.append(None)
    
    # For each key residue, compute contact frequency by activity class
    for res_name in KEY_RESIDUES["alpha9"]:
        if res_name not in fingerprints.get(list(fingerprints.keys())[0], {}):
            continue
        
        # Contact presence (1/0) for each compound
        contact_present = []
        valid_activities = []
        
        for i, c in enumerate(compounds):
            if c["id"] in fingerprints:
                fp = fingerprints[c["id"]]
                if res_name in fp:
                    contact_present.append(1 if fp[res_name]["contact"] else 0)
                    if activities[i] is not None:
                        valid_activities.append(activities[i])
                    else:
                        valid_activities.append(0)  # inactive = log(0) undefined, use 0
        
        if len(contact_present) > 2:
            # Point-biserial correlation (binary vs continuous)
            try:
                r, p = stats.pointbiserialr(contact_present, valid_activities)
                results["correlations"][res_name] = {
                    "r": float(r),
                    "p_value": float(p),
                    "n_contacts": sum(contact_present),
                    "contact_rate": float(np.mean(contact_present)),
                }
            except:
                pass
    
    # Active vs inactive discrimination
    active_fps = []
    inactive_fps = []
    
    for c in compounds:
        if c["id"] in fingerprints:
            fp = fingerprints[c["id"]]
            n_contacts = sum(1 for v in fp.values() if v.get("contact", False))
            if c["is_active"]:
                active_fps.append(n_contacts)
            else:
                inactive_fps.append(n_contacts)
    
    if active_fps and inactive_fps:
        try:
            t_stat, p_val = stats.ttest_ind(active_fps, inactive_fps)
            results["discrimination"] = {
                "active_mean_contacts": float(np.mean(active_fps)),
                "inactive_mean_contacts": float(np.mean(inactive_fps)),
                "t_statistic": float(t_stat),
                "p_value": float(p_val),
                "effect_size_cohens_d": float(
                    (np.mean(active_fps) - np.mean(inactive_fps)) /
                    np.sqrt((np.var(active_fps) + np.var(inactive_fps)) / 2)
                ),
            }
        except:
            pass
    
    return results


def main():
    print("=" * 60)
    print("Interaction Fingerprint Analysis")
    print("=" * 60)
    
    compounds = load_compounds()
    print(f"Loaded {len(compounds)} compounds")
    
    # Process existing PDB results if available
    pdb_dir = RESULTS_DIR / "boltz2"
    if not pdb_dir.exists():
        pdb_dir = RESULTS_DIR / "af3"
    
    all_fingerprints = {}
    site_results = {}
    
    # Process each site
    for site_id in [23, 5, 21]:
        print(f"\n--- Site {site_id} ---")
        site_dir = pdb_dir / f"site{site_id}"
        
        if not site_dir.exists():
            print(f"  No results found at {site_dir}")
            continue
        
        fingerprints = {}
        for compound in compounds:
            # Look for PDB files
            pdb_files = list(site_dir.glob(f"*cpd{compound['id']}*/*.pdb"))
            pdb_files.extend(site_dir.glob(f"*cpd{compound['id']}*.pdb"))
            
            if pdb_files:
                contacts = parse_pdb_contacts(pdb_files[0])
                fp = build_interaction_fingerprint(contacts, KEY_RESIDUES)
                fingerprints[compound["id"]] = fp
                print(f"  Compound {compound['id']}: {len(contacts)} contacts")
            else:
                print(f"  Compound {compound['id']}: No PDB found")
        
        if fingerprints:
            all_fingerprints[site_id] = fingerprints
            
            # Compute SAR correlation
            sar_results = compute_sar_correlation(compounds, fingerprints, site_id)
            site_results[site_id] = sar_results
            
            print(f"\n  SAR Analysis (Site {site_id}):")
            print(f"    Active: {sar_results['n_active']}, Inactive: {sar_results['n_inactive']}")
            if sar_results["discrimination"]:
                d = sar_results["discrimination"]
                print(f"    Active mean contacts: {d.get('active_mean_contacts', 'N/A'):.1f}")
                print(f"    Inactive mean contacts: {d.get('inactive_mean_contacts', 'N/A'):.1f}")
                print(f"    Cohen's d: {d.get('effect_size_cohens_d', 'N/A'):.3f}")
                print(f"    p-value: {d.get('p_value', 'N/A'):.4f}")
    
    # Save results
    output_path = ANALYSIS_DIR / "interaction_fingerprint_analysis.json"
    with open(output_path, "w") as f:
        json.dump({
            "site_results": site_results,
            "n_compounds": len(compounds),
            "compounds": compounds,
        }, f, indent=2, default=str)
    
    print(f"\nResults saved to: {output_path}")
    
    # Generate CSV summary
    csv_path = ANALYSIS_DIR / "interaction_fingerprints.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "site_id", "compound_id", "activity_class", "is_active",
            "activity_uM", "n_contacts", "key_residue_contacts"
        ])
        
        for site_id, fps in all_fingerprints.items():
            for compound in compounds:
                if compound["id"] in fps:
                    fp = fps[compound["id"]]
                    n_contacts = sum(1 for v in fp.values() if v.get("contact", False))
                    key_contacts = [k for k, v in fp.items() if v.get("contact", False)]
                    writer.writerow([
                        site_id, compound["id"], compound["activity_class"],
                        compound["is_active"], compound["activity_uM"],
                        n_contacts, ";".join(key_contacts)
                    ])
    
    print(f"Fingerprints CSV: {csv_path}")


if __name__ == "__main__":
    main()
