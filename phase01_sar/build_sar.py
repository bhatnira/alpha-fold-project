#!/usr/bin/env python3
"""Phase I — Build the α9α10 SAR dataset.

Parses:
  FINAL_SITE_RANKING.csv  → per-site composite scores and discriminations
  LIGAND_COMPARISON.csv   → per-site per-ligand enrichment, contacts, stoichiometry
  MUTAGENESIS_PREDICTIONS.csv → per-site residue predictions

Outputs:
  sar_dataset.csv          — full SAR matrix (ligand × site × stoichiometry)
  discovery_set.csv        — discovery split (80%)
  held_out_set.csv         — held-out split (20%)
  chemical_space.csv       — per-ligand molecular descriptors
  sar_summary.txt          — human-readable summary
"""

import csv, sys, os, json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT  = Path(__file__).resolve().parent

# ── Load site ranking ──────────────────────────────────────────────────────
def load_site_ranking():
    sites = {}
    with open(DATA / "FINAL_SITE_RANKING.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sites[int(row["site_id"])] = {
                "rank": int(row["rank_primary"]),
                "composite_score": float(row["SiteScore"]),
                "n_models": int(row["n_models"]),
                "n_af3": int(row["n_af3"]),
                "n_boltz2": int(row["n_boltz2"]),
                "ligands": row["ligands"].split(";"),
                "domains": row["domains"],
                "level": row["evidence_level"],
                "sufficient_support": row["sufficient_support"] == "True",
                "stereo_test": row["stereo_test"],
                "acetate_test": row["acetate_test"],
                "analogue_test": row["analogue_test"],
                "comparative_class": row["comparative_class"],
            }
    return sites

# ── Load ligand comparison ─────────────────────────────────────────────────
def load_ligand_comparison():
    rows = []
    with open(DATA / "LIGAND_COMPARISON.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "site": int(row["site"]),
                "ligand_class": row["ligand_class"],
                "n_models": int(row["n_models"]),
                "n_models_af3": int(row["n_models_af3"]),
                "n_models_boltz2": int(row["n_models_boltz2"]),
                "site_share": float(row["site_share"]),
                "enrichment": float(row["enrichment"]),
                "fisher_p": float(row["fisher_p_greater"]),
                "mean_confidence": float(row["mean_confidence"]),
                "stoichiometry": row["stoichiometry"],
                "n_contact_residues": int(row["n_contact_residues_recurrent"]),
                "recurrent_contacts": row["recurrent_contacts"].split(";") if row["recurrent_contacts"] else [],
                "recurrent_freqs": [
                    float(x.split(":")[-1])
                    for x in row["recurrent_contact_freqs"].split(";")
                    if ":" in x
                ] if row["recurrent_contact_freqs"] else [],
            })
    return rows

# ── Load mutagenesis predictions ──────────────────────────────────────────
def load_mutagenesis():
    preds = []
    with open(DATA / "MUTAGENESIS_PREDICTIONS.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            preds.append({
                "site": int(row["site"]),
                "subunit": row["subunit"],
                "residue_number": int(row["residue_number"]),
                "residue_identity": row["residue_identity"],
                "domain": row["domain"],
                "contact_frequency": float(row["contact_frequency"]),
                "classification": row["classification"],
                "anchor_type": row["anchor_type"],
                "priority": row["priority"],
                "hypothesis_A": row["prediction_if_A_specific_recognition"],
                "hypothesis_B": row["prediction_if_B_electrostatic"],
            })
    return preds

# ── Ligand metadata ────────────────────────────────────────────────────────
LIGAND_META = {
    "L-ASC": {
        "full_name": "L-ascorbate",
        "smiles": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
        "mw": 176.12,
        "charge": -1,
        "class": "active",
        "scaffold": "ascorbate",
        "stereochem": "L",
        "logp": -1.85,
        "activity_rank": 1,
    },
    "D-ASC": {
        "full_name": "D-ascorbate",
        "smiles": "OC[C@@H](O)[C@@H]1OC(=O)C(O)=C1O",
        "mw": 176.12,
        "charge": -1,
        "class": "weak",
        "scaffold": "ascorbate",
        "stereochem": "D",
        "logp": -1.85,
        "activity_rank": 3,
    },
    "O-ETHYL": {
        "full_name": "O-ethyl ascorbate",
        "smiles": "CCOC(=O)C(=O)OC[C@@H](O)[C@@H]1OC(=O)C(O)=C1O",
        "mw": 204.17,
        "charge": 0,
        "class": "active",
        "scaffold": "ascorbate",
        "stereochem": "L",
        "logp": -0.92,
        "activity_rank": 2,
    },
    "ACETATE": {
        "full_name": "acetate",
        "smiles": "CC(=O)[O-]",
        "mw": 59.04,
        "charge": -1,
        "class": "control",
        "scaffold": "acetate",
        "stereochem": "none",
        "logp": -0.17,
        "activity_rank": 0,
    },
    "RYANODINE": {
        "full_name": "ryanodine",
        "smiles": "CC1CCC2C(C)(C)CCC3C(C)(C)CC1C23OC4C(O)C(CO)OC4O",
        "mw": 433.52,
        "charge": 0,
        "class": "modulator",
        "scaffold": "ryanoid",
        "stereochem": "natural",
        "logp": 0.84,
        "activity_rank": 1,
    },
}

# ── Build SAR matrix ───────────────────────────────────────────────────────
def build_sar_matrix(sites, ligand_comp):
    """Ligand × site × stoichiometry matrix with enrichment, contacts, score."""
    matrix = []
    for lc in ligand_comp:
        site_id = lc["site"]
        if site_id not in sites:
            continue
        site = sites[site_id]
        lig = lc["ligand_class"]
        meta = LIGAND_META.get(lig, {})
        matrix.append({
            "site_id": site_id,
            "site_rank": site["rank"],
            "site_composite_score": site["composite_score"],
            "ligand": lig,
            "ligand_full_name": meta.get("full_name", lig),
            "ligand_class": meta.get("class", "unknown"),
            "ligand_charge": meta.get("charge", "unknown"),
            "ligand_scaffold": meta.get("scaffold", "unknown"),
            "ligand_stereochem": meta.get("stereochem", "unknown"),
            "stoichiometry": lc["stoichiometry"],
            "n_models": lc["n_models"],
            "n_models_af3": lc["n_models_af3"],
            "n_models_boltz2": lc["n_models_boltz2"],
            "enrichment": lc["enrichment"],
            "fisher_p": lc["fisher_p"],
            "mean_confidence": lc["mean_confidence"],
            "n_contacts": lc["n_contact_residues"],
            "recurrent_contacts": ";".join(lc["recurrent_contacts"][:10]),
            "site_evidence_level": site["level"],
            "site_sufficient": site["sufficient_support"],
            "stereo_test": site["stereo_test"],
            "acetate_test": site["acetate_test"],
            "analogue_test": site["analogue_test"],
        })
    return matrix

# ── Split discovery / held-out ─────────────────────────────────────────────
def split_dataset(matrix, sites):
    """80/20 split by site rank. Top sites go to discovery, tail to held-out."""
    # sort by site rank
    ranked = sorted(matrix, key=lambda r: (r["site_rank"], r["ligand"]))
    n = len(ranked)
    cutoff = int(n * 0.8)
    discovery = ranked[:cutoff]
    held_out = ranked[cutoff:]
    return discovery, held_out

# ── Write CSV ──────────────────────────────────────────────────────────────
def write_csv(rows, path):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

# ── Summary ────────────────────────────────────────────────────────────────
def write_summary(sites, matrix, discovery, held_out, path):
    lines = ["α9α10 SAR DATASET SUMMARY", "=" * 60, ""]
    lines.append(f"Total sites:          {len(sites)}")
    lines.append(f"Total SAR records:    {len(matrix)}")
    lines.append(f"Discovery records:    {len(discovery)}")
    lines.append(f"Held-out records:     {len(held_out)}")
    lines.append("")

    # Ligand breakdown
    by_lig = defaultdict(list)
    for r in matrix:
        by_lig[r["ligand"]].append(r)
    lines.append("Per-ligand breakdown:")
    for lig, rows in sorted(by_lig.items()):
        sites_visited = len(set(r["site_id"] for r in rows))
        total_models = sum(r["n_models"] for r in rows)
        lines.append(f"  {lig:12s}  {len(rows):3d} records  {sites_visited:2d} sites  {total_models:4d} models")
    lines.append("")

    # Stoichiometry breakdown
    by_sto = defaultdict(list)
    for r in matrix:
        by_sto[r["stoichiometry"]].append(r)
    lines.append("Per-stoichiometry breakdown:")
    for sto, rows in sorted(by_sto.items()):
        lines.append(f"  {sto:6s}  {len(rows):3d} records  {sum(r['n_models'] for r in rows):4d} models")
    lines.append("")

    # Enrichment distribution
    enrichments = [r["enrichment"] for r in matrix]
    lines.append(f"Enrichment range:     {min(enrichments):.2f} – {max(enrichments):.2f}")
    lines.append(f"Median enrichment:    {sorted(enrichments)[len(enrichments)//2]:.2f}")
    lines.append("")

    # Top 5 enriched ligand×site pairs
    top5 = sorted(matrix, key=lambda r: r["enrichment"], reverse=True)[:5]
    lines.append("Top 5 enriched ligand×site pairs:")
    for r in top5:
        lines.append(f"  site {r['site_id']:2d}  {r['ligand']:12s}  enrich={r['enrichment']:.2f}  "
                      f"models={r['n_models']:3d}  contacts={r['n_contacts']}")
    lines.append("")

    # Sites with sufficient support
    sufficient = [s for s in sites.values() if s["sufficient_support"]]
    lines.append(f"Sites with sufficient support: {len(sufficient)} / {len(sites)}")
    for s in sufficient:
        lines.append(f"  site {s['rank']:2d}  composite={s['composite_score']:.4f}  "
                      f"level={s['level']}  class={s['comparative_class']}")
    lines.append("")

    # Held-out sites
    held_out_sites = sorted(set(r["site_id"] for r in held_out))
    lines.append(f"Held-out sites ({len(held_out_sites)}): {held_out_sites}")

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

# ── Main ───────────────────────────────────────────────────────────────────
def main():
    print("Loading site ranking...")
    sites = load_site_ranking()
    print(f"  {len(sites)} sites")

    print("Loading ligand comparison...")
    lig_comp = load_ligand_comparison()
    print(f"  {len(lig_comp)} records")

    print("Loading mutagenesis predictions...")
    mut = load_mutagenesis()
    print(f"  {len(mut)} residue predictions")

    print("Building SAR matrix...")
    matrix = build_sar_matrix(sites, lig_comp)
    print(f"  {len(matrix)} records")

    print("Splitting discovery / held-out...")
    discovery, held_out = split_dataset(matrix, sites)
    print(f"  discovery={len(discovery)}, held_out={len(held_out)}")

    print("Writing outputs...")
    write_csv(matrix, OUT / "sar_dataset.csv")
    write_csv(discovery, OUT / "discovery_set.csv")
    write_csv(held_out, OUT / "held_out_set.csv")
    write_summary(sites, matrix, discovery, held_out, OUT / "sar_summary.txt")

    print("Done. Outputs in", OUT)

if __name__ == "__main__":
    main()
