#!/usr/bin/env python3
"""Phase II — Receptor ensemble inventory and selection.

Inventory all apo structures from the campaign, select a compact
representative ensemble by stoichiometry and conformation.
"""

import csv, os, glob
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

def find_apo_structures():
    """Find all apo receptor PDBs in the deliverable."""
    apo_dir = DATA / "deliverable" / "05_structures"
    pdbs = []
    if apo_dir.exists():
        for f in sorted(apo_dir.rglob("*.pdb")):
            pdbs.append(f)
    # Also check the publication set
    pub_dir = DATA / "campaign" / ".." / "alpha9alpha10_publication" / "01_receptor"
    if pub_dir.exists():
        for f in sorted(pub_dir.rglob("*.pdb")):
            pdbs.append(f)
    return pdbs

def inventory():
    """Create an inventory of all receptor structures."""
    records = []
    
    # Check deliverable structures
    deliverable = DATA / "deliverable" / "05_structures"
    if deliverable.exists():
        for pdb in sorted(deliverable.rglob("*.pdb")):
            rel = pdb.relative_to(deliverable)
            stoich = "unknown"
            if "2to3" in str(rel):
                stoich = "2to3"
            elif "3to2" in str(rel):
                stoich = "3to2"
            
            state = "unknown"
            if "apo" in str(rel).lower():
                state = "apo"
            elif "holo" in str(rel).lower():
                state = "holo"
            
            records.append({
                "path": str(pdb),
                "relative_path": str(rel),
                "source": "deliverable",
                "stoichiometry": stoich,
                "state": state,
                "size_kb": os.path.getsize(pdb) / 1024,
            })
    
    # Check publication receptor folder
    pub_receptor = ROOT.parent / "alpha9alpha10_publication" / "01_receptor"
    if pub_receptor.exists():
        for pdb in sorted(pub_receptor.rglob("*.pdb")):
            rel = pdb.relative_to(pub_receptor)
            stoich = "unknown"
            if "2to3" in str(rel):
                stoich = "2to3"
            elif "3to2" in str(rel):
                stoich = "3to2"
            
            records.append({
                "path": str(pdb),
                "relative_path": str(rel),
                "source": "publication",
                "stoichiometry": stoich,
                "state": "apo",
                "size_kb": os.path.getsize(pdb) / 1024,
            })
    
    # Check AF3 apo models
    af3_dir = DATA / "allostery" / "af3_outputs"
    if af3_dir.exists():
        apo_jobs = [d for d in af3_dir.iterdir() if d.is_dir() and "apo" in d.name.lower()]
        for job in sorted(apo_jobs)[:20]:  # sample first 20
            for cif in sorted(job.rglob("*.cif"))[:5]:
                records.append({
                    "path": str(cif),
                    "relative_path": str(cif.relative_to(af3_dir)),
                    "source": "af3",
                    "stoichiometry": "unknown",
                    "state": "apo",
                    "size_kb": os.path.getsize(cif) / 1024,
                })
    
    return records

def write_inventory(records, path):
    if not records:
        print("No structures found!")
        return
    
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)

def select_ensemble(records):
    """Select a compact representative ensemble.
    
    Strategy: pick top-ranked apo reference per stoichiometry + 
    a few representative AF3 holo models (highest confidence).
    """
    selected = []
    
    # Always include the publication apo references
    for r in records:
        if r["source"] == "publication" and r["state"] == "apo":
            selected.append(r)
    
    # Add top AF3 models by source
    af3_holo = [r for r in records if r["source"] == "af3" and r["state"] == "holo"]
    af3_holo.sort(key=lambda r: r["size_kb"], reverse=True)  # proxy for quality
    
    seen_stoich = set()
    for r in af3_holo:
        stoich = r["stoichiometry"]
        if stoich not in seen_stoich and len(selected) < 10:
            selected.append(r)
            seen_stoich.add(stoich)
    
    return selected

def write_summary(records, selected, path):
    lines = [
        "RECEPTOR ENSEMBLE INVENTORY",
        "=" * 60,
        "",
        f"Total structures found: {len(records)}",
        "",
        "By source:",
    ]
    
    by_source = defaultdict(list)
    for r in records:
        by_source[r["source"]].append(r)
    for src, recs in sorted(by_source.items()):
        lines.append(f"  {src:15s}  {len(recs):3d} files")
    
    lines.extend(["", "By stoichiometry:"])
    by_stoich = defaultdict(list)
    for r in records:
        by_stoich[r["stoichiometry"]].append(r)
    for sto, recs in sorted(by_stoich.items()):
        lines.append(f"  {sto:10s}  {len(recs):3d} files")
    
    lines.extend(["", "Selected ensemble:", ""])
    for r in selected:
        lines.append(f"  {r['source']:12s}  {r['stoichiometry']:6s}  {r['state']:6s}  {r['relative_path']}")
    
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    print("Inventorying receptor structures...")
    records = inventory()
    print(f"  Found {len(records)} structures")
    
    print("Selecting compact ensemble...")
    selected = select_ensemble(records)
    print(f"  Selected {len(selected)} structures")
    
    print("Writing outputs...")
    write_inventory(records, OUT / "receptor_inventory.csv")
    write_inventory(selected, OUT / "selected_ensemble.csv")
    write_summary(records, selected, OUT / "receptor_summary.txt")
    
    print("Done.")

if __name__ == "__main__":
    main()
