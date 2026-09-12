#!/usr/bin/env python3
"""Phase IX — Specificity / falsification analysis.

Compares L-ascorbate vs acetate (charge control),
L-ascorbate vs D-ascorbate (stereochemical control),
and parent vs O-ethyl ascorbate (chemical perturbation).
"""

import csv, sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def load_ligand_comparison():
    rows = []
    with open(DATA / "LIGAND_COMPARISON.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "site": int(row["site"]),
                "ligand_class": row["ligand_class"],
                "n_models": int(row["n_models"]),
                "enrichment": float(row["enrichment"]),
                "mean_confidence": float(row["mean_confidence"]),
                "stoichiometry": row["stoichiometry"],
                "n_contacts": int(row["n_contact_residues_recurrent"]),
                "recurrent_contacts": row["recurrent_contacts"].split(";") if row["recurrent_contacts"] else [],
                "recurrent_freqs": [
                    float(x.split(":")[-1])
                    for x in row["recurrent_contact_freqs"].split(";")
                    if ":" in x
                ] if row["recurrent_contact_freqs"] else [],
            })
    return rows

def specificity_analysis(rows):
    """Compare controls vs active compounds."""
    by_site = defaultdict(list)
    for r in rows:
        by_site[r["site"]].append(r)
    
    results = []
    for site_id, site_records in sorted(by_site.items()):
        by_ligand = defaultdict(list)
        for r in site_records:
            by_ligand[r["ligand_class"]].append(r)
        
        # L-ascorbate vs acetate
        l_asc = by_ligand.get("L-ASC", [])
        acetate = by_ligand.get("ACETATE", [])
        
        l_asc_contacts = set()
        for r in l_asc:
            l_asc_contacts.update(r["recurrent_contacts"])
        acetate_contacts = set()
        for r in acetate:
            acetate_contacts.update(r["recurrent_contacts"])
        
        shared = l_asc_contacts & acetate_contacts
        total = l_asc_contacts | acetate_contacts
        jaccard = len(shared) / len(total) if total else 0
        
        # L-ascorbate vs D-ascorbate
        d_asc = by_ligand.get("D-ASC", [])
        d_asc_contacts = set()
        for r in d_asc:
            d_asc_contacts.update(r["recurrent_contacts"])
        
        ld_shared = l_asc_contacts & d_asc_contacts
        ld_total = l_asc_contacts | d_asc_contacts
        ld_jaccard = len(ld_shared) / len(ld_total) if ld_total else 0
        
        # Enrichment comparison
        l_asc_enrich = max([r["enrichment"] for r in l_asc], default=0)
        acetate_enrich = max([r["enrichment"] for r in acetate], default=0)
        d_asc_enrich = max([r["enrichment"] for r in d_asc], default=0)
        
        # Discrimination scores
        charge_discrimination = 1.0 - jaccard  # high = good (acetate ≠ L-asc)
        stereo_discrimination = 1.0 - ld_jaccard  # high = good (D ≠ L)
        
        # Classification
        if charge_discrimination >= 0.7 and l_asc_enrich > acetate_enrich:
            charge_class = "MOLECULAR_RECOGNITION"
        elif charge_discrimination >= 0.4:
            charge_class = "PARTIAL_DISCRIMINATION"
        else:
            charge_class = "CHARGE_DRIVEN"
        
        if stereo_discrimination >= 0.5 and l_asc_enrich > d_asc_enrich:
            stereo_class = "STEREOSPECIFIC"
        elif stereo_discrimination >= 0.3:
            stereo_class = "PARTIAL_STEREOSELECTIVITY"
        else:
            stereo_class = "NON_STEREOSPECIFIC"
        
        results.append({
            "site": site_id,
            "l_asc_models": len(l_asc),
            "acetate_models": len(acetate),
            "d_asc_models": len(d_asc),
            "l_asc_enrichment": round(l_asc_enrich, 3),
            "acetate_enrichment": round(acetate_enrich, 3),
            "d_asc_enrichment": round(d_asc_enrich, 3),
            "l_asc_contacts": ";".join(sorted(l_asc_contacts)[:15]),
            "acetate_contacts": ";".join(sorted(acetate_contacts)[:15]),
            "d_asc_contacts": ";".join(sorted(d_asc_contacts)[:15]),
            "shared_l_asc_acetate": ";".join(sorted(shared)[:10]),
            "charge_jaccard": round(jaccard, 4),
            "charge_discrimination": round(charge_discrimination, 4),
            "charge_class": charge_class,
            "stereo_jaccard": round(ld_jaccard, 4),
            "stereo_discrimination": round(stereo_discrimination, 4),
            "stereo_class": stereo_class,
        })
    
    return results

def write_csv(rows, path):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

def write_summary(results, path):
    lines = [
        "SPECIFICITY / FALSIFICATION ANALYSIS",
        "=" * 60,
        "",
        "Controls:",
        "  L-ascorbate vs acetate  →  charge-only vs molecular recognition",
        "  L-ascorbate vs D-ascorbate  →  stereochemical specificity",
        "",
    ]
    
    # Charge discrimination summary
    charge_classes = defaultdict(list)
    for r in results:
        charge_classes[r["charge_class"]].append(r["site"])
    
    lines.append("Charge discrimination (L-asc vs acetate):")
    for cls, sites in sorted(charge_classes.items()):
        lines.append(f"  {cls:25s}  {len(sites):3d} sites  {sites[:8]}")
    lines.append("")
    
    # Stereo discrimination summary
    stereo_classes = defaultdict(list)
    for r in results:
        stereo_classes[r["stereo_class"]].append(r["site"])
    
    lines.append("Stereochemical discrimination (L-asc vs D-asc):")
    for cls, sites in sorted(stereo_classes.items()):
        lines.append(f"  {cls:25s}  {len(sites):3d} sites  {sites[:8]}")
    lines.append("")
    
    # Top discriminating sites
    lines.append("Top sites by charge discrimination:")
    by_charge = sorted(results, key=lambda r: r["charge_discrimination"], reverse=True)
    for r in by_charge[:5]:
        lines.append(f"  site {r['site']:2d}  disc={r['charge_discrimination']:.3f}  "
                      f"jaccard={r['charge_jaccard']:.3f}  class={r['charge_class']}")
    lines.append("")
    
    lines.append("Top sites by stereo discrimination:")
    by_stereo = sorted(results, key=lambda r: r["stereo_discrimination"], reverse=True)
    for r in by_stereo[:5]:
        lines.append(f"  site {r['site']:2d}  disc={r['stereo_discrimination']:.3f}  "
                      f"jaccard={r['stereo_jaccard']:.3f}  class={r['stereo_class']}")
    
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    print("Loading ligand comparison data...")
    rows = load_ligand_comparison()
    print(f"  {len(rows)} records")
    
    print("Running specificity analysis...")
    results = specificity_analysis(rows)
    print(f"  {len(results)} sites analyzed")
    
    print("Writing outputs...")
    write_csv(results, OUT / "specificity_results.csv")
    write_summary(results, OUT / "specificity_summary.txt")
    
    print("Done.")

if __name__ == "__main__":
    main()
