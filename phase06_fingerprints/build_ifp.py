#!/usr/bin/env python3
"""Phase VI — Build interaction fingerprints from complex PDBs.

Parses contact data from site_clusters_v3.csv and builds
per-site residue-level contact matrices (interaction fingerprints).
"""

import csv, sys, os
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def load_site_clusters():
    """Load per-model site assignments with contacts."""
    records = []
    with open(DATA / "site_clusters_v3.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            contacts = row["contacts"].split(";") if row["contacts"] else []
            # Parse contacts into (subunit, residue_number) pairs
            parsed = []
            for c in contacts:
                if ":" in c:
                    parts = c.split(":")
                    parsed.append((parts[0], int(parts[1])))
            records.append({
                "site": int(row["site"]),
                "method": row["method"],
                "stoichiometry": row["stoichiometry"],
                "ligand": row["ligand"],
                "seed": int(row["seed"]),
                "sample": int(row["sample"]),
                "n_contacts": int(row["n_contact_residues"]) if row["n_contact_residues"] else 0,
                "contacts": parsed,
                "subunits": row["subunits"],
                "domains": row["domains"],
            })
    return records

def build_contact_matrix(records):
    """Build per-site residue contact frequency matrix."""
    by_site = defaultdict(list)
    for r in records:
        by_site[r["site"]].append(r)
    
    matrices = {}
    for site_id, site_records in sorted(by_site.items()):
        # Count contact frequency per residue
        residue_counts = defaultdict(lambda: {"total": 0, "by_ligand": defaultdict(int), "by_method": defaultdict(int)})
        
        n_total = len(site_records)
        for r in site_records:
            seen = set()
            for subunit, resnum in r["contacts"]:
                key = f"{subunit}:{resnum}"
                if key not in seen:
                    residue_counts[key]["total"] += 1
                    residue_counts[key]["by_ligand"][r["ligand"]] += 1
                    residue_counts[key]["by_method"][r["method"]] += 1
                    seen.add(key)
        
        # Convert to sorted list
        matrix = []
        for reskey, counts in sorted(residue_counts.items(), 
                                      key=lambda x: x[1]["total"], 
                                      reverse=True):
            freq = counts["total"] / n_total if n_total > 0 else 0
            matrix.append({
                "site": site_id,
                "residue": reskey,
                "subunit": reskey.split(":")[0],
                "residue_number": int(reskey.split(":")[1]),
                "contact_frequency": round(freq, 4),
                "n_models_contact": counts["total"],
                "n_models_total": n_total,
                "n_ligands_contact": len(counts["by_ligand"]),
                "ligands": ";".join(sorted(counts["by_ligand"].keys())),
                "ligand_freqs": ";".join(
                    f"{k}:{round(v/n_total, 3)}" 
                    for k, v in sorted(counts["by_ligand"].items(), 
                                        key=lambda x: x[1], reverse=True)
                ),
                "af3_freq": round(counts["by_method"].get("AF3", 0) / n_total, 4),
                "boltz2_freq": round(counts["by_method"].get("Boltz2", 0) / n_total, 4),
            })
        
        matrices[site_id] = matrix
    
    return matrices

def build_ifp_summary(matrices):
    """Create a summary of interaction fingerprints per site."""
    summaries = []
    for site_id, matrix in sorted(matrices.items()):
        n_residues = len(matrix)
        core = [r for r in matrix if r["contact_frequency"] >= 0.5]
        peripheral = [r for r in matrix if 0.2 <= r["contact_frequency"] < 0.5]
        
        # Domain breakdown
        ecd = sum(1 for r in matrix if r["subunit"].startswith("alpha") and r["residue_number"] < 200)
        tm = sum(1 for r in matrix if r["residue_number"] >= 200)
        
        summaries.append({
            "site": site_id,
            "n_contact_residues": n_residues,
            "n_core": len(core),
            "n_peripheral": len(peripheral),
            "n_ECD": ecd,
            "n_TM": tm,
            "core_residues": ";".join(r["residue"] for r in core[:20]),
            "top_residue": matrix[0]["residue"] if matrix else "",
            "top_frequency": matrix[0]["contact_frequency"] if matrix else 0,
        })
    
    return summaries

def write_csv(rows, path):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

def write_matrix_csv(matrix, path):
    """Write the contact matrix for one site."""
    if not matrix:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=matrix[0].keys())
        writer.writeheader()
        writer.writerows(matrix)

def write_summary(matrices, summaries, path):
    lines = [
        "INTERACTION FINGERPRINT SUMMARY",
        "=" * 60,
        "",
    ]
    
    for site_id in sorted(matrices.keys()):
        matrix = matrices[site_id]
        summ = next(s for s in summaries if s["site"] == site_id)
        
        lines.append(f"Site {site_id}:")
        lines.append(f"  Contact residues: {summ['n_contact_residues']}")
        lines.append(f"  Core (≥50%):      {summ['n_core']}")
        lines.append(f"  Peripheral (20-50%): {summ['n_peripheral']}")
        lines.append(f"  ECD:              {summ['n_ECD']}")
        lines.append(f"  TM:               {summ['n_TM']}")
        lines.append(f"  Top residue:      {summ['top_residue']} ({summ['top_frequency']:.3f})")
        
        # Show core residues
        core = [r for r in matrix if r["contact_frequency"] >= 0.5]
        if core:
            lines.append(f"  Core residues:")
            for r in core[:10]:
                lines.append(f"    {r['residue']:15s}  freq={r['contact_frequency']:.3f}  "
                              f"ligands={r['n_ligands_contact']}")
        lines.append("")
    
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    print("Loading site clusters...")
    records = load_site_clusters()
    print(f"  {len(records)} model records")
    
    print("Building contact matrices...")
    matrices = build_contact_matrix(records)
    print(f"  {len(matrices)} sites")
    
    print("Building IFP summaries...")
    summaries = build_ifp_summary(matrices)
    
    print("Writing outputs...")
    # Write per-site contact matrices
    for site_id, matrix in matrices.items():
        write_matrix_csv(matrix, OUT / f"site{site_id:02d}_contacts.csv")
    
    # Write summaries
    write_csv(summaries, OUT / "ifp_summary.csv")
    write_summary(matrices, summaries, OUT / "ifp_summary.txt")
    
    print("Done.")

if __name__ == "__main__":
    main()
