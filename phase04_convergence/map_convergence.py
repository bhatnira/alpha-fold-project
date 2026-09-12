#!/usr/bin/env python3
"""Phase IV — AF3 + Boltz-2 convergence map.

Parses site_clusters_v3.csv to determine whether AF3 and Boltz-2
converge on the same spatial region for each ligand.
"""

import csv, sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def load_site_clusters():
    """Load per-model site assignments."""
    records = []
    with open(DATA / "site_clusters_v3.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            contacts = row["contacts"].split(";") if row["contacts"] else []
            records.append({
                "site": int(row["site"]),
                "pose_cluster": int(row["pose_cluster"]),
                "method": row["method"],
                "stoichiometry": row["stoichiometry"],
                "ligand": row["ligand"],
                "seed": int(row["seed"]),
                "sample": int(row["sample"]),
                "receptor_rmsd": float(row["receptor_rmsd"]),
                "chain_map": row["chain_map"],
                "confidence": float(row["confidence"]) if row["confidence"] else None,
                "ligand_confidence": float(row["ligand_confidence"]) if row["ligand_confidence"] else None,
                "n_contacts": int(row["n_contact_residues"]) if row["n_contact_residues"] else 0,
                "subunits": row["subunits"],
                "domains": row["domains"],
                "cx": float(row["cx"]) if row["cx"] else None,
                "cy": float(row["cy"]) if row["cy"] else None,
                "cz": float(row["cz"]) if row["cz"] else None,
                "contacts": contacts,
            })
    return records

def compute_convergence(records):
    """Compute AF3 vs Boltz-2 convergence per site per ligand."""
    # Group by site
    by_site = defaultdict(list)
    for r in records:
        by_site[r["site"]].append(r)
    
    convergence = []
    for site_id, site_records in sorted(by_site.items()):
        # Group by ligand
        by_ligand = defaultdict(list)
        for r in site_records:
            by_ligand[r["ligand"]].append(r)
        
        for ligand, lig_records in by_ligand.items():
            af3 = [r for r in lig_records if r["method"] == "AF3"]
            boltz = [r for r in lig_records if r["method"] == "Boltz2"]
            
            # Centroid overlap
            af3_centroids = [(r["cx"], r["cy"], r["cz"]) for r in af3 if r["cx"] is not None]
            boltz_centroids = [(r["cx"], r["cy"], r["cz"]) for r in boltz if r["cx"] is not None]
            
            # Contact overlap
            af3_contacts = set()
            for r in af3:
                af3_contacts.update(r["contacts"])
            boltz_contacts = set()
            for r in boltz:
                boltz_contacts.update(r["contacts"])
            
            # Jaccard similarity of contacts
            if af3_contacts or boltz_contacts:
                jaccard = len(af3_contacts & boltz_contacts) / len(af3_contacts | boltz_contacts)
            else:
                jaccard = 0.0
            
            # Subunit overlap
            af3_subunits = set()
            for r in af3:
                af3_subunits.update(r["subunits"].split("+"))
            boltz_subunits = set()
            for r in boltz:
                boltz_subunits.update(r["subunits"].split("+"))
            
            convergence.append({
                "site": site_id,
                "ligand": ligand,
                "n_af3": len(af3),
                "n_boltz2": len(boltz),
                "n_total": len(lig_records),
                "af3_contact_residues": ";".join(sorted(af3_contacts)),
                "boltz2_contact_residues": ";".join(sorted(boltz_contacts)),
                "shared_contacts": ";".join(sorted(af3_contacts & boltz_contacts)),
                "contact_jaccard": round(jaccard, 4),
                "af3_subunits": "+".join(sorted(af3_subunits)),
                "boltz2_subunits": "+".join(sorted(boltz_subunits)),
                "subunit_overlap": bool(af3_subunits & boltz_subunits),
                "af3_centroid_count": len(af3_centroids),
                "boltz2_centroid_count": len(boltz_centroids),
            })
    
    return convergence

def classify_convergence(conv):
    """Classify convergence strength."""
    for r in conv:
        if r["contact_jaccard"] >= 0.5 and r["subunit_overlap"]:
            r["convergence_class"] = "STRONG"
        elif r["contact_jaccard"] >= 0.2 or r["subunit_overlap"]:
            r["convergence_class"] = "MODERATE"
        elif r["n_af3"] > 0 and r["n_boltz2"] > 0:
            r["convergence_class"] = "WEAK"
        else:
            r["convergence_class"] = "SINGLE_METHOD"
    return conv

def write_csv(rows, path):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

def write_summary(conv, records, path):
    lines = [
        "AF3 + BOLTZ-2 CONVERGENCE MAP",
        "=" * 60,
        "",
        f"Total model records:     {len(records)}",
        f"AF3 records:             {sum(1 for r in records if r['method'] == 'AF3')}",
        f"Boltz2 records:          {sum(1 for r in records if r['method'] == 'Boltz2')}",
        f"Unique site×ligand:      {len(conv)}",
        "",
    ]
    
    # Convergence classes
    by_class = defaultdict(list)
    for r in conv:
        by_class[r["convergence_class"]].append(r)
    
    lines.append("Convergence classification:")
    for cls, items in sorted(by_class.items()):
        sites = sorted(set(r["site"] for r in items))
        lines.append(f"  {cls:15s}  {len(items):3d} pairs  sites={sites}")
    lines.append("")
    
    # Top convergent sites
    strong = [r for r in conv if r["convergence_class"] == "STRONG"]
    strong.sort(key=lambda r: r["contact_jaccard"], reverse=True)
    lines.append("Strong convergence sites:")
    for r in strong[:10]:
        lines.append(f"  site {r['site']:2d}  {r['ligand']:12s}  "
                      f"jaccard={r['contact_jaccard']:.3f}  "
                      f"shared={len(r['shared_contacts'].split(';')) if r['shared_contacts'] else 0}")
    lines.append("")
    
    # Sites with only AF3 (no Boltz-2)
    af3_only = [r for r in conv if r["n_boltz2"] == 0]
    lines.append(f"Sites with only AF3 (no Boltz-2): {len(af3_only)}")
    for r in af3_only[:5]:
        lines.append(f"  site {r['site']:2d}  {r['ligand']:12s}  n_af3={r['n_af3']}")
    
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    print("Loading site clusters...")
    records = load_site_clusters()
    print(f"  {len(records)} model records")
    
    print("Computing convergence...")
    conv = compute_convergence(records)
    print(f"  {len(conv)} site×ligand pairs")
    
    print("Classifying convergence...")
    conv = classify_convergence(conv)
    
    print("Writing outputs...")
    write_csv(conv, OUT / "convergence_map.csv")
    write_summary(conv, records, OUT / "convergence_summary.txt")
    
    print("Done.")

if __name__ == "__main__":
    main()
