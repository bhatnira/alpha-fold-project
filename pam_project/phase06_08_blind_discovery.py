#!/usr/bin/env python
"""
PHASE 6-8: Blind PAM Site Discovery + Candidate Site Generation
Master Prompt Compliance: Sections 11-12

Extracts ligand-binding locations from 468 existing AF3 CIF models.
Clusters ligand positions to identify candidate PAM sites.
Works with what exists: 2 stoichiometries × 9 ligand conditions × multiple seeds.
"""

import json
import csv
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
AF3_DIR = PIPELINE_DIR / "data/allostery/af3_outputs"
OUTPUT_DIR = PIPELINE_DIR / "pam_project/06_af3_discovery"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_cif_contacts(cif_path):
    """
    Parse CIF file to extract ligand atoms and receptor contacts.
    Returns ligand center of mass and contact residues.
    
    AF3 CIF format: HETATM lines with ligand residues (LIG_F, etc.)
    """
    ligand_atoms = []
    receptor_atoms = []
    
    try:
        with open(cif_path, "r") as f:
            in_atom = False
            atom_start = False
            columns = []
            
            for line in f:
                line = line.rstrip()
                
                # Find atom_site loop
                if "_atom_site.group_PDB" in line:
                    in_atom = True
                    continue
                
                if in_atom and not atom_start:
                    # Collect column names
                    if line.startswith("_atom_site."):
                        columns.append(line.strip())
                        continue
                    else:
                        atom_start = True
                
                if not atom_start:
                    continue
                
                # Parse atom line
                parts = line.split()
                if len(parts) < 15:
                    continue
                
                try:
                    group = parts[0]  # HETATM or ATOM
                    atom_id = parts[1]
                    type_symbol = parts[2]
                    atom_name = parts[3]
                    comp_id = parts[5]  # residue name
                    asym_id = parts[6]  # chain
                    seq_id = parts[8]  # sequence ID
                    
                    # Coordinates
                    x = float(parts[10])
                    y = float(parts[11])
                    z = float(parts[12])
                    
                    atom_info = {
                        "group": group,
                        "chain": asym_id,
                        "res_name": comp_id,
                        "res_num": int(seq_id) if seq_id != "." else 0,
                        "atom_name": atom_name,
                        "coords": np.array([x, y, z]),
                    }
                    
                    # Ligand residues start with LIG or are HETATM
                    if group == "HETATM" or comp_id.startswith("LIG"):
                        ligand_atoms.append(atom_info)
                    else:
                        receptor_atoms.append(atom_info)
                except (ValueError, IndexError):
                    continue
    
    except Exception as e:
        return None, None, None, []
    
    if not ligand_atoms:
        return None, None, None, []
    
    # Calculate ligand center of mass
    lig_coords = np.array([a["coords"] for a in ligand_atoms])
    lig_com = np.mean(lig_coords, axis=0)
    
    # Find receptor contacts (within 8 A of ligand COM)
    contacts = []
    for rec_atom in receptor_atoms:
        dist = np.linalg.norm(rec_atom["coords"] - lig_com)
        if dist < 8.0:
            contacts.append({
                "chain": rec_atom["chain"],
                "res_name": rec_atom["res_name"],
                "res_num": rec_atom["res_num"],
                "distance": float(dist),
                "atom_name": rec_atom["atom_name"],
            })
    
    # Group by residue
    residue_contacts = defaultdict(list)
    for c in contacts:
        key = f"{c['chain']}:{c['res_name']}{c['res_num']}"
        residue_contacts[key].append(c)
    
    # Get closest residues
    closest = sorted(residue_contacts.keys(), 
                     key=lambda k: min(c["distance"] for c in residue_contacts[k]))
    
    return lig_com, lig_coords, closest[:20], list(residue_contacts.keys())


def parse_pdb_contacts(pdb_path, cutoff=8.0):
    """Parse PDB file for ligand contacts."""
    ligand_atoms = []
    receptor_atoms = []
    
    try:
        with open(pdb_path, "r") as f:
            for line in f:
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
                    
                    if line.startswith("HETATM"):
                        ligand_atoms.append(atom_info)
                    else:
                        receptor_atoms.append(atom_info)
    except:
        return None, None, []
    
    if not ligand_atoms:
        return None, None, []
    
    lig_coords = np.array([a["coords"] for a in ligand_atoms])
    lig_com = np.mean(lig_coords, axis=0)
    
    residue_contacts = defaultdict(list)
    for rec_atom in receptor_atoms:
        dist = np.linalg.norm(rec_atom["coords"] - lig_com)
        if dist < cutoff:
            key = f"{rec_atom['chain']}:{rec_atom['res_name']}{rec_atom['res_num']}"
            residue_contacts[key].append({"distance": float(dist), **rec_atom})
    
    closest = sorted(residue_contacts.keys(),
                     key=lambda k: min(c["distance"] for c in residue_contacts[k]))
    
    return lig_com, lig_coords, closest[:20]


def scan_af3_outputs():
    """Scan all AF3 output directories and extract ligand binding information."""
    results = []
    
    for condition_dir in sorted(AF3_DIR.iterdir()):
        if not condition_dir.is_dir():
            continue
        
        condition_name = condition_dir.name
        
        # Parse condition name
        parts = condition_name.split("_", 2)
        stoich = parts[0] + "_" + parts[1] if len(parts) > 1 else "unknown"
        ligand = parts[2] if len(parts) > 2 else "unknown"
        
        # Find ALL CIF files (not just ranked model)
        cif_files = list(condition_dir.rglob("*.cif"))
        
        for cif_path in cif_files:
            # Determine if this is a ranked model or seed/sample
            is_ranked = "model" in cif_path.name and "seed" not in str(cif_path)
            seed_info = None
            if "seed-" in str(cif_path):
                parts_path = str(cif_path).split("seed-")
                if len(parts_path) > 1:
                    seed_num = parts_path[1].split("_")[0]
                    seed_info = f"seed-{seed_num}"
            
            lig_com, lig_coords, top_contacts, all_contacts = parse_cif_contacts(cif_path)
            
            if lig_com is not None:
                results.append({
                    "condition": condition_name,
                    "stoichiometry": stoich,
                    "ligand_type": ligand,
                    "cif_file": str(cif_path),
                    "is_ranked_model": is_ranked,
                    "seed_info": seed_info,
                    "lig_com_x": float(lig_com[0]),
                    "lig_com_y": float(lig_com[1]),
                    "lig_com_z": float(lig_com[2]),
                    "n_contacts": len(all_contacts),
                    "top_contacts": top_contacts,
                    "all_contacts": all_contacts,
                })
    
    return results


def cluster_binding_sites(results, distance_threshold=10.0):
    """Cluster ligand positions to identify candidate binding sites."""
    
    # Collect all ligand centers
    centers = []
    for r in results:
        centers.append(np.array([r["lig_com_x"], r["lig_com_y"], r["lig_com_z"]]))
    
    if not centers:
        return []
    
    centers = np.array(centers)
    
    # Simple greedy clustering
    clusters = []
    assigned = set()
    
    for i in range(len(centers)):
        if i in assigned:
            continue
        
        cluster = [i]
        assigned.add(i)
        
        for j in range(i + 1, len(centers)):
            if j in assigned:
                continue
            
            dist = np.linalg.norm(centers[i] - centers[j])
            if dist < distance_threshold:
                cluster.append(j)
                assigned.add(j)
        
        clusters.append(cluster)
    
    # Build cluster info
    cluster_info = []
    for ci, cluster_indices in enumerate(clusters):
        cluster_centers = centers[cluster_indices]
        centroid = np.mean(cluster_centers, axis=0)
        
        # Get contact residue frequency
        all_contacts = []
        for idx in cluster_indices:
            all_contacts.extend(results[idx]["top_contacts"])
        
        contact_freq = Counter(all_contacts)
        top_residues = [r for r, _ in contact_freq.most_common(15)]
        
        # Ligand types in this cluster
        ligand_types = Counter(results[idx]["ligand_type"] for idx in cluster_indices)
        stoichiometries = Counter(results[idx]["stoichiometry"] for idx in cluster_indices)
        
        cluster_info.append({
            "cluster_id": ci,
            "n_models": len(cluster_indices),
            "centroid": centroid.tolist(),
            "ligand_types": dict(ligand_types),
            "stoichiometries": dict(stoichiometries),
            "top_residues": top_residues,
            "contact_frequency": dict(contact_freq.most_common(20)),
            "model_indices": cluster_indices,
            "recurrence_across_stoichiometries": len(stoichiometries),
            "recurrence_across_ligands": len(ligand_types),
        })
    
    # Sort by size (most recurrent first)
    cluster_info.sort(key=lambda x: -x["n_models"])
    
    return cluster_info


def identify_interfaces(cluster_info):
    """Identify which subunit interface each cluster represents."""
    
    for cluster in cluster_info:
        residues = cluster["top_residues"]
        
        # Determine which subunits are involved
        subunits = set()
        for r in residues:
            chain = r.split(":")[0]
            subunits.add(chain)
        
        # Classify interface
        subunit_list = sorted(subunits)
        if len(subunit_list) == 1:
            interface = f"intra_{subunit_list[0]}"
        elif len(subunit_list) == 2:
            interface = f"{subunit_list[0]}_{subunit_list[1]}"
        else:
            interface = "multi_subunit"
        
        cluster["interface"] = interface
        cluster["subunits"] = subunit_list
        
        # Classify pocket type based on residue composition
        # Alpha9 key residues: W176, S175, Y120, Y224, Y217, W55
        # Alpha10 key residues: D145, R83, W81, R143
        alpha9_residues = ["W176", "S175", "Y120", "Y224", "Y217", "W55", "C192", "Y93", "Y197"]
        alpha10_residues = ["D145", "R83", "W81", "R143", "H153"]
        
        a9_count = sum(1 for r in residues if any(ar in r for ar in alpha9_residues))
        a10_count = sum(1 for r in residues if any(ar in r for ar in alpha10_residues))
        
        if a9_count > 0 and a10_count > 0:
            cluster["pocket_type"] = "inter_subunit"
        elif a9_count > 0:
            cluster["pocket_type"] = "intra_alpha9"
        elif a10_count > 0:
            cluster["pocket_type"] = "intra_alpha10"
        else:
            cluster["pocket_type"] = "other"
    
    return cluster_info


def main():
    print("=" * 70)
    print("PHASE 6-8: BLIND PAM SITE DISCOVERY")
    print("=" * 70)
    
    # Scan AF3 outputs
    print("\nScanning AF3 output directories...")
    results = scan_af3_outputs()
    print(f"Found {len(results)} ranked models with ligand coordinates")
    
    if not results:
        print("No ligand coordinates found in AF3 outputs")
        print("Checking PDB files instead...")
        
        # Try PDB files
        pdb_dir = PIPELINE_DIR / "data/deliverable/05_structures"
        pdb_files = list(pdb_dir.rglob("*.pdb"))
        print(f"Found {len(pdb_files)} PDB files")
        
        for pdb in pdb_files[:5]:
            lig_com, lig_coords, contacts = parse_pdb_contacts(pdb)
            if lig_com is not None:
                results.append({
                    "condition": pdb.stem,
                    "stoichiometry": "unknown",
                    "ligand_type": "unknown",
                    "cif_file": str(pdb),
                    "is_ranked_model": True,
                    "lig_com_x": float(lig_com[0]),
                    "lig_com_y": float(lig_com[1]),
                    "lig_com_z": float(lig_com[2]),
                    "n_contacts": len(contacts),
                    "top_contacts": contacts,
                    "all_contacts": contacts,
                })
        
        print(f"Extracted {len(results)} models from PDB files")
    
    # Summary by condition
    print("\n--- AF3 Model Summary ---")
    by_condition = defaultdict(list)
    for r in results:
        by_condition[r["condition"]].append(r)
    
    for cond, models in sorted(by_condition.items()):
        print(f"  {cond}: {len(models)} models")
    
    # Cluster binding sites
    print("\n--- Clustering Binding Sites ---")
    clusters = cluster_binding_sites(results, distance_threshold=12.0)
    print(f"Identified {len(clusters)} candidate binding site clusters")
    
    # Identify interfaces
    clusters = identify_interfaces(clusters)
    
    # Report top candidates
    print("\n--- TOP CANDIDATE SITES ---")
    for i, cluster in enumerate(clusters[:10]):
        print(f"\n  Site {i+1} (Cluster {cluster['cluster_id']}):")
        print(f"    Models: {cluster['n_models']}")
        print(f"    Interface: {cluster['interface']}")
        print(f"    Pocket type: {cluster['pocket_type']}")
        print(f"    Stoichiometries: {cluster['stoichiometries']}")
        print(f"    Ligand types: {cluster['ligand_types']}")
        print(f"    Top residues: {cluster['top_residues'][:10]}")
        print(f"    Centroid: [{cluster['centroid'][0]:.1f}, {cluster['centroid'][1]:.1f}, {cluster['centroid'][2]:.1f}]")
    
    # Save results
    output = {
        "n_total_models": len(results),
        "n_clusters": len(clusters),
        "clusters": clusters,
        "all_models": results,
    }
    
    with open(OUTPUT_DIR / "af3_site_discovery.json", "w") as f:
        json.dump(output, f, indent=2, default=str)
    
    # Save cluster summary CSV
    with open(OUTPUT_DIR / "candidate_sites.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "rank", "cluster_id", "n_models", "interface", "pocket_type",
            "stoichiometries", "ligand_types", "top_residues", "centroid_x", "centroid_y", "centroid_z"
        ])
        for i, c in enumerate(clusters):
            writer.writerow([
                i+1, c["cluster_id"], c["n_models"], c["interface"], c["pocket_type"],
                str(c["stoichiometries"]), str(c["ligand_types"]),
                ";".join(c["top_residues"]),
                c["centroid"][0], c["centroid"][1], c["centroid"][2],
            ])
    
    print(f"\nResults saved to: {OUTPUT_DIR}")
    
    return clusters


if __name__ == "__main__":
    main()
