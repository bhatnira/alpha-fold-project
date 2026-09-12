#!/usr/bin/env python3
"""Phase XI — Robustness without MD.

Vary AF3 seed, Boltz-2 seed, receptor conformation, stoichiometry,
ligand representation, docking seed, docking protocol, candidate pocket.
Calculate site recurrence, pose recurrence, residue recurrence, framework recurrence.
"""

import csv, sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def load_site_clusters():
    records = []
    with open(DATA / "site_clusters_v3.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            contacts = row["contacts"].split(";") if row["contacts"] else []
            records.append({
                "site": int(row["site"]),
                "method": row["method"],
                "stoichiometry": row["stoichiometry"],
                "ligand": row["ligand"],
                "seed": int(row["seed"]),
                "sample": int(row["sample"]),
                "receptor_rmsd": float(row["receptor_rmsd"]),
                "n_contacts": int(row["n_contact_residues"]) if row["n_contact_residues"] else 0,
                "contacts": contacts,
                "subunits": row["subunits"],
            })
    return records

def compute_robustness(records):
    """Compute robustness metrics across seeds, states, stoichiometries."""
    by_site = defaultdict(list)
    for r in records:
        by_site[r["site"]].append(r)
    
    robustness = []
    for site_id, site_records in sorted(by_site.items()):
        n_total = len(site_records)
        
        # Site recurrence: how many unique seed×sample combinations
        unique_configs = set((r["seed"], r["sample"]) for r in site_records)
        
        # AF3 vs Boltz2 recurrence
        af3 = [r for r in site_records if r["method"] == "AF3"]
        boltz = [r for r in site_records if r["method"] == "Boltz2"]
        
        # Stoichiometry recurrence
        stoichs = set(r["stoichiometry"] for r in site_records)
        
        # Ligand recurrence
        ligands = set(r["ligand"] for r in site_records)
        
        # Residue recurrence
        all_contacts = defaultdict(int)
        for r in site_records:
            for c in r["contacts"]:
                all_contacts[c] += 1
        
        core_residues = [res for res, count in all_contacts.items() 
                         if count / n_total >= 0.5]
        peripheral_residues = [res for res, count in all_contacts.items()
                               if 0.2 <= count / n_total < 0.5]
        
        # RMSD spread (lower = more robust)
        rmsds = [r["receptor_rmsd"] for r in site_records]
        mean_rmsd = sum(rmsds) / len(rmsds) if rmsds else 0
        rmsd_spread = max(rmsds) - min(rmsds) if len(rmsds) > 1 else 0
        
        # Framework recurrence (subunit composition)
        frameworks = set(r["subunits"] for r in site_records if r["subunits"])
        
        # Robustness score
        seed_score = len(unique_configs) / max(n_total, 1)
        method_score = min(len(af3), len(boltz)) / max(min(len(af3), 1), len(boltz), 1) if boltz else 0.5
        stoich_score = len(stoichs) / 2  # max 2 stoichiometries
        residue_score = len(core_residues) / max(len(all_contacts), 1)
        
        overall_score = (seed_score + method_score + stoich_score + residue_score) / 4
        
        robustness.append({
            "site": site_id,
            "n_models": n_total,
            "n_unique_configs": len(unique_configs),
            "n_af3": len(af3),
            "n_boltz2": len(boltz),
            "n_stoichiometries": len(stoichs),
            "stoichiometries": ";".join(sorted(stoichs)),
            "n_ligands": len(ligands),
            "ligands": ";".join(sorted(ligands)),
            "n_core_residues": len(core_residues),
            "n_peripheral_residues": len(peripheral_residues),
            "core_residues": ";".join(sorted(core_residues)[:10]),
            "mean_rmsd": round(mean_rmsd, 3),
            "rmsd_spread": round(rmsd_spread, 3),
            "n_frameworks": len(frameworks),
            "frameworks": ";".join(sorted(frameworks)),
            "seed_score": round(seed_score, 4),
            "method_score": round(method_score, 4),
            "stoich_score": round(stoich_score, 4),
            "residue_score": round(residue_score, 4),
            "overall_robustness": round(overall_score, 4),
        })
    
    return robustness

def write_csv(rows, path):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

def write_summary(robustness, path):
    lines = [
        "COMPUTATIONAL ROBUSTNESS MATRIX",
        "=" * 60,
        "",
        "Robustness = persistence across seeds, methods, stoichiometries, conformations",
        "No MD required — uses replication instead of dynamics",
        "",
    ]
    
    # Overall robustness distribution
    scores = [r["overall_robustness"] for r in robustness]
    lines.append(f"Overall robustness: {min(scores):.3f} – {max(scores):.3f} "
                  f"(median {sorted(scores)[len(scores)//2]:.3f})")
    lines.append("")
    
    # Top robust sites
    lines.append("Top 10 most robust sites:")
    top10 = sorted(robustness, key=lambda r: r["overall_robustness"], reverse=True)[:10]
    for r in top10:
        lines.append(f"  site {r['site']:2d}  robustness={r['overall_robustness']:.3f}  "
                      f"models={r['n_models']:3d}  seeds={r['n_unique_configs']:2d}  "
                      f"stoich={r['n_stoichiometries']}  core={r['n_core_residues']}")
    lines.append("")
    
    # Component breakdown
    lines.append("Component scores (site, seed, method, stoich, residue):")
    for r in sorted(robustness, key=lambda x: x["overall_robustness"], reverse=True)[:5]:
        lines.append(f"  site {r['site']:2d}  seed={r['seed_score']:.3f}  "
                      f"method={r['method_score']:.3f}  "
                      f"stoich={r['stoich_score']:.3f}  "
                      f"residue={r['residue_score']:.3f}")
    
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    print("Loading site clusters...")
    records = load_site_clusters()
    print(f"  {len(records)} records")
    
    print("Computing robustness...")
    robustness = compute_robustness(records)
    print(f"  {len(robustness)} sites")
    
    print("Writing outputs...")
    write_csv(robustness, OUT / "robustness_matrix.csv")
    write_summary(robustness, OUT / "robustness_summary.txt")
    
    print("Done.")

if __name__ == "__main__":
    main()
