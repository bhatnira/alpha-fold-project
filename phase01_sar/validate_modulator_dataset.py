#!/usr/bin/env python3
"""Phase I — Validate binding model against experimental modulator dataset.

Correlates:
  modulator-dataset-a9a10.csv (EC50, %Potentiation)
  with computational enrichment predictions (sar_dataset.csv)

Outputs:
  modulator_validation.csv  — per-compound validation metrics
  modulator_validation.txt  — summary report
  site_validation.csv       — per-site validation against experimental data
"""

import csv, sys, math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT  = Path(__file__).resolve().parent

# ── Load modulator dataset ─────────────────────────────────────────────────
def load_modulators():
    """Load experimental EC50 and %Potentiation data."""
    mods = []
    path = ROOT / "modulator-dataset-a9a10.csv"
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ec50 = float(row["Activity (uM)"])
            pot = float(row["%Potentiation"])
            mods.append({
                "compound_id": int(row["Identifier"]),
                "smiles": row["Smiles"],
                "ec50_um": ec50,
                "potentiation_pct": pot,
                "is_active": ec50 > 0 and pot > 0,
                "is_inactive": ec50 == 0 and pot == 0,
                "pIC50": -math.log10(ec50 * 1e-6) if ec50 > 0 else 0,
            })
    return mods

# ── Load SAR dataset ───────────────────────────────────────────────────────
def load_sar():
    """Load computational enrichment predictions per site per ligand."""
    sar = []
    path = OUT / "sar_dataset.csv"
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sar.append({
                "site_id": int(row["site_id"]),
                "site_rank": int(row["site_rank"]),
                "site_composite_score": float(row["site_composite_score"]),
                "ligand": row["ligand"],
                "ligand_class": row["ligand_class"],
                "stoichiometry": row["stoichiometry"],
                "n_models": int(row["n_models"]),
                "enrichment": float(row["enrichment"]),
                "fisher_p": float(row["fisher_p"]),
                "n_contacts": int(row["n_contacts"]),
                "recurrent_contacts": row["recurrent_contacts"],
            })
    return sar

# ── Load evidence matrix ──────────────────────────────────────────────────
def load_evidence():
    """Load integrated evidence scores."""
    ev = {}
    path = ROOT / "phase14_evidence/evidence_matrix.csv"
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ev[int(row["site_id"])] = {
                "evidence_score": float(row["evidence_score"]),
                "classification": row["classification"],
            }
    return ev

# ── SMILES similarity (Tanimoto on fragments) ─────────────────────────────
def smiles_fragments(smiles):
    """Simple fragment fingerprint from SMILES string."""
    frags = set()
    for i in range(len(smiles) - 1):
        frags.add(smiles[i:i+2])
    for i in range(len(smiles) - 2):
        frags.add(smiles[i:i+3])
    return frags

def tanimoto(smiles_a, smiles_b):
    """Tanimoto similarity between two SMILES fragment sets."""
    a = smiles_fragments(smiles_a)
    b = smiles_fragments(smiles_b)
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)

# ── Map modulator compounds to SAR ligands ────────────────────────────────
SAR_SMILES = {
    "L-ASC": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
    "D-ASC": "OC[C@@H](O)[C@@H]1OC(=O)C(O)=C1O",
    "O-ETHYL": "CCOC(=O)C(=O)OC[C@@H](O)[C@@H]1OC(=O)C(O)=C1O",
    "ACETATE": "CC(=O)[O-]",
    "RYANODINE": "CC1CCC2C(C)(C)CCC3C(C)(C)CC1C23OC4C(O)C(CO)OC4O",
}

def auto_map_compound(smiles, threshold=0.45):
    """Map a compound to the closest SAR ligand by SMILES Tanimoto similarity."""
    best_ligand = None
    best_sim = 0
    for lig, ref_smiles in SAR_SMILES.items():
        sim = tanimoto(smiles, ref_smiles)
        if sim > best_sim:
            best_sim = sim
            best_ligand = lig
    if best_sim >= threshold:
        return best_ligand, best_sim
    return None, best_sim

# ── Core validation logic ─────────────────────────────────────────────────
def validate_compounds(mods, sar, evidence):
    """For each modulator, find matching SAR data and validate."""
    results = []

    # Group SAR by ligand
    sar_by_ligand = defaultdict(list)
    for row in sar:
        sar_by_ligand[row["ligand"]].append(row)

    # Group SAR by site
    sar_by_site = defaultdict(lambda: defaultdict(list))
    for row in sar:
        sar_by_site[row["site_id"]][row["ligand"]].append(row)

    # L-ASC enrichment per site (for correlation)
    l_asc_by_site = {}
    for row in sar:
        if row["ligand"] == "L-ASC":
            key = (row["site_id"], row["stoichiometry"])
            if key not in l_asc_by_site or row["enrichment"] > l_asc_by_site[key]["enrichment"]:
                l_asc_by_site[key] = row

    for mod in mods:
        cid = mod["compound_id"]

        # Auto-map using SMILES similarity
        sar_ligand, sim_score = auto_map_compound(mod["smiles"])

        # Find best matching SAR data
        best_enrichment = 0
        best_site = None
        best_fisher = 1.0
        matched_sites = []

        if sar_ligand and sar_ligand in sar_by_ligand:
            for row in sar_by_ligand[sar_ligand]:
                if row["enrichment"] > best_enrichment:
                    best_enrichment = row["enrichment"]
                    best_site = row["site_id"]
                    best_fisher = row["fisher_p"]
                if row["enrichment"] > 1.0:
                    matched_sites.append(row["site_id"])

        # Evidence score at best site
        ev_score = 0
        ev_class = "UNKNOWN"
        if best_site and best_site in evidence:
            ev_score = evidence[best_site]["evidence_score"]
            ev_class = evidence[best_site]["classification"]

        # SMILES similarity to L-ASC
        l_asc_smiles = "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O"
        sim_to_lasc = tanimoto(mod["smiles"], l_asc_smiles)

        results.append({
            "compound_id": cid,
            "smiles": mod["smiles"],
            "ec50_um": mod["ec50_um"],
            "pIC50": mod["pIC50"],
            "potentiation_pct": mod["potentiation_pct"],
            "is_active": mod["is_active"],
            "is_inactive": mod["is_inactive"],
            "mapped_sar_ligand": sar_ligand or "UNMAPPED",
            "mapping_similarity": round(sim_score, 4),
            "best_enrichment": best_enrichment,
            "best_site": best_site,
            "best_fisher_p": best_fisher,
            "n_enriched_sites": len(matched_sites),
            "evidence_score_at_best_site": ev_score,
            "evidence_class_at_best_site": ev_class,
            "similarity_to_lasc": sim_to_lasc,
        })

    return results

# ── Site-level validation ─────────────────────────────────────────────────
def validate_sites(sar, evidence, mods):
    """For each site, check if enrichment predicts experimental activity."""
    # Group SAR by site
    site_data = defaultdict(lambda: {"l_asc": [], "d_asc": [], "o_ethyl": [], "acetate": []})
    ligand_key = {"L-ASC": "l_asc", "D-ASC": "d_asc", "O-ETHYL": "o_ethyl", "ACETATE": "acetate"}
    for row in sar:
        site_id = int(row["site_id"])
        lig = row["ligand"]
        key = ligand_key.get(lig)
        if key:
            site_data[site_id][key].append(row)

    site_results = []
    for site_id, ligands in sorted(site_data.items()):
        l_asc_enrich = max((r["enrichment"] for r in ligands["l_asc"]), default=0)
        d_asc_enrich = max((r["enrichment"] for r in ligands["d_asc"]), default=0)
        o_ethyl_enrich = max((r["enrichment"] for r in ligands["o_ethyl"]), default=0)
        acetate_enrich = max((r["enrichment"] for r in ligands["acetate"]), default=0)

        # Specificity ratio: L-ASC vs acetate
        lasc_vs_acetate = l_asc_enrich / acetate_enrich if acetate_enrich > 0 else 0

        # Stereo selectivity: L-ASC vs D-ASC
        stereo_ratio = l_asc_enrich / d_asc_enrich if d_asc_enrich > 0 else 0

        ev = evidence.get(site_id, {})

        # Molecular recognition: L-ASC enriches > acetate OR shows stereo preference
        has_lasc_enrichment = l_asc_enrich > 1.0
        has_lasc_advantage = lasc_vs_acetate > 1.2 if acetate_enrich > 0 else has_lasc_enrichment
        has_stereo_pref = stereo_ratio > 1.1 if d_asc_enrich > 0 else True

        site_results.append({
            "site_id": site_id,
            "evidence_score": ev.get("evidence_score", 0),
            "evidence_class": ev.get("classification", "UNKNOWN"),
            "l_asc_enrichment": l_asc_enrich,
            "d_asc_enrichment": d_asc_enrich,
            "o_ethyl_enrichment": o_ethyl_enrich,
            "acetate_enrichment": acetate_enrich,
            "lasc_vs_acetate_ratio": lasc_vs_acetate,
            "stereo_selectivity_ratio": stereo_ratio,
            "molecular_recognition": has_lasc_enrichment and (has_lasc_advantage or has_stereo_pref),
        })

    return site_results

# ── Statistical correlation ───────────────────────────────────────────────
def spearman_rho(x, y):
    """Spearman rank correlation coefficient."""
    n = len(x)
    if n < 3:
        return 0.0

    def rankdata(vals):
        indexed = sorted(enumerate(vals), key=lambda t: t[1])
        ranks = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j < n - 1 and indexed[j+1][1] == indexed[j][1]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                ranks[indexed[k][0]] = avg_rank
            i = j + 1
        return ranks

    rx = rankdata(x)
    ry = rankdata(y)

    d2 = sum((rx[i] - ry[i])**2 for i in range(n))
    return 1 - (6 * d2) / (n * (n**2 - 1))

# ── Write outputs ─────────────────────────────────────────────────────────
def write_csv(rows, path):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

def write_summary(compounds, sites, path):
    lines = ["MODULATOR DATASET VALIDATION", "=" * 60, ""]

    # Overall statistics
    active = [c for c in compounds if c["is_active"]]
    inactive = [c for c in compounds if c["is_inactive"]]
    mapped = [c for c in compounds if c["mapped_sar_ligand"] != "UNMAPPED"]

    lines.append(f"Total compounds:       {len(compounds)}")
    lines.append(f"Active (EC50>0, pot>0): {len(active)}")
    lines.append(f"Inactive (EC50=0):     {len(inactive)}")
    lines.append(f"Mapped to SAR ligand:  {len(mapped)}")
    lines.append("")

    # Mapping quality
    sims = [c["mapping_similarity"] for c in compounds if c["mapped_sar_ligand"] != "UNMAPPED"]
    if sims:
        lines.append(f"SMILES mapping quality:")
        lines.append(f"  Mean Tanimoto similarity: {sum(sims)/len(sims):.4f}")
        lines.append(f"  Min: {min(sims):.4f}  Max: {max(sims):.4f}")
        lines.append("")

    # Enrichment correlation for active compounds
    active_with_enrich = [c for c in active if c["best_enrichment"] > 0]
    if len(active_with_enrich) >= 3:
        pic50_vals = [c["pIC50"] for c in active_with_enrich]
        enrich_vals = [c["best_enrichment"] for c in active_with_enrich]
        rho = spearman_rho(pic50_vals, enrich_vals)
        lines.append(f"CORRELATION: pIC50 vs enrichment (active compounds with enrichment)")
        lines.append(f"  Spearman rho = {rho:.4f}")
        lines.append(f"  N = {len(active_with_enrich)}")
        lines.append(f"  Interpretation: {'STRONG' if abs(rho) > 0.6 else 'MODERATE' if abs(rho) > 0.4 else 'WEAK'}")
        lines.append("")
    elif active_with_enrich:
        lines.append(f"CORRELATION: N={len(active_with_enrich)} active compounds with enrichment (need >= 3)")
        lines.append("")

    # Enrichment correlation for potentiation
    if len(active_with_enrich) >= 3:
        pot_vals = [c["potentiation_pct"] for c in active_with_enrich]
        enrich_vals = [c["best_enrichment"] for c in active_with_enrich]
        rho = spearman_rho(pot_vals, enrich_vals)
        lines.append(f"CORRELATION: %Potentiation vs enrichment (active compounds with enrichment)")
        lines.append(f"  Spearman rho = {rho:.4f}")
        lines.append(f"  N = {len(active_with_enrich)}")
        lines.append(f"  Interpretation: {'STRONG' if abs(rho) > 0.6 else 'MODERATE' if abs(rho) > 0.4 else 'WEAK'}")
        lines.append("")

    # Active vs inactive discrimination
    active_enrichments = [c["best_enrichment"] for c in active]
    inactive_enrichments = [c["best_enrichment"] for c in inactive]
    active_ev_scores = [c["evidence_score_at_best_site"] for c in active]
    inactive_ev_scores = [c["evidence_score_at_best_site"] for c in inactive]

    if active_enrichments and inactive_enrichments:
        lines.append("ACTIVE vs INACTIVE DISCRIMINATION:")
        lines.append(f"  Active  enrichment: mean={sum(active_enrichments)/len(active_enrichments):.2f}  "
                      f"max={max(active_enrichments):.2f}")
        lines.append(f"  Inactive enrichment: mean={sum(inactive_enrichments)/len(inactive_enrichments):.2f}  "
                      f"max={max(inactive_enrichments):.2f}")
        lines.append(f"  Active  evidence:   mean={sum(active_ev_scores)/len(active_ev_scores):.4f}")
        lines.append(f"  Inactive evidence:  mean={sum(inactive_ev_scores)/len(inactive_ev_scores):.4f}")
        lines.append("")

    # Top active compounds
    lines.append("TOP ACTIVE COMPOUNDS (by enrichment):")
    lines.append(f"{'ID':>4s} {'EC50(uM)':>10s} {'pIC50':>7s} {'%Pot':>7s} {'Enrich':>8s} {'Site':>5s} {'Evidence':>8s} {'Mapped':>10s} {'Sim':>6s}")
    for c in sorted(active, key=lambda x: x["best_enrichment"], reverse=True)[:10]:
        site_str = str(c['best_site']) if c['best_site'] else 'N/A'
        lines.append(f"{c['compound_id']:4d} {c['ec50_um']:10.2f} {c['pIC50']:7.2f} {c['potentiation_pct']:7.1f} "
                      f"{c['best_enrichment']:8.2f} {site_str:>5s} {c['evidence_score_at_best_site']:8.4f} "
                      f"{c['mapped_sar_ligand']:>10s} {c['mapping_similarity']:6.4f}")
    lines.append("")

    # Inactive compounds
    unmapped_inactive = [c for c in inactive if c["mapped_sar_ligand"] == "UNMAPPED"]
    mapped_inactive = [c for c in inactive if c["mapped_sar_ligand"] != "UNMAPPED"]
    lines.append(f"INACTIVE COMPOUNDS ({len(inactive)}): {len(mapped_inactive)} mapped, {len(unmapped_inactive)} unmapped")
    for c in mapped_inactive:
        site_str = str(c['best_site']) if c['best_site'] else 'N/A'
        lines.append(f"  ID={c['compound_id']:3d}  enrichment={c['best_enrichment']:.2f}  "
                      f"site={site_str}  mapped={c['mapped_sar_ligand']}  sim={c['mapping_similarity']:.4f}")
    lines.append("")

    # Site-level validation
    mol_recog_sites = [s for s in sites if s["molecular_recognition"]]
    lines.append(f"SITES SHOWING MOLECULAR RECOGNITION (L-ASC enriches + stereo preference): {len(mol_recog_sites)}")
    for s in sorted(mol_recog_sites, key=lambda x: x["evidence_score"], reverse=True)[:15]:
        site_str = f"  site {s['site_id']:2d}  ev={s['evidence_score']:.4f}  "
        site_str += f"L-ASC={s['l_asc_enrichment']:.2f}  acet={s['acetate_enrichment']:.2f}  "
        site_str += f"L/D={s['stereo_selectivity_ratio']:.2f}  L/acet={s['lasc_vs_acetate_ratio']:.2f}"
        lines.append(site_str)
    lines.append("")

    # Summary conclusion
    lines.append("VALIDATION CONCLUSION:")
    if mol_recog_sites:
        top = sorted(mol_recog_sites, key=lambda x: x["evidence_score"], reverse=True)[0]
        lines.append(f"  Best validated site: {top['site_id']} (evidence={top['evidence_score']:.4f})")
        lines.append(f"  L-ASC enrichment: {top['l_asc_enrichment']:.2f}")
        lines.append(f"  Stereo selectivity (L/D): {top['stereo_selectivity_ratio']:.2f}")
    else:
        lines.append("  No sites pass molecular recognition criteria.")
        lines.append("  This suggests either:")
        lines.append("    1. The computational model does not distinguish active from inactive compounds")
        lines.append("    2. The enrichment metric is not sensitive enough")
        lines.append("    3. The binding site is not correctly identified")
    lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

# ── Main ───────────────────────────────────────────────────────────────────
def main():
    print("Loading modulator dataset...")
    mods = load_modulators()
    print(f"  {len(mods)} compounds")

    print("Loading SAR dataset...")
    sar = load_sar()
    print(f"  {len(sar)} records")

    print("Loading evidence matrix...")
    evidence = load_evidence()
    print(f"  {len(evidence)} sites")

    print("Validating compounds...")
    compounds = validate_compounds(mods, sar, evidence)
    print(f"  {len(compounds)} validated")

    print("Validating sites...")
    sites = validate_sites(sar, evidence, mods)
    print(f"  {len(sites)} sites analyzed")

    print("Writing outputs...")
    write_csv(compounds, OUT / "modulator_validation.csv")
    write_csv(sites, OUT / "site_validation.csv")
    write_summary(compounds, sites, OUT / "modulator_validation.txt")

    # Print key results
    active = [c for c in compounds if c["is_active"]]
    mapped_active = [c for c in active if c["mapped_sar_ligand"] != "UNMAPPED"]
    print(f"\n  Active compounds: {len(active)}")
    print(f"  Mapped to SAR:    {len(mapped_active)}")

    if mapped_active:
        pic50 = [c["pIC50"] for c in mapped_active if c["pIC50"] > 0]
        enrich = [c["best_enrichment"] for c in mapped_active if c["best_enrichment"] > 0]
        if len(pic50) == len(enrich) and len(pic50) >= 3:
            rho = spearman_rho(pic50, enrich)
            print(f"  Spearman rho (pIC50 vs enrichment): {rho:.4f}")

    mol_recog = [s for s in sites if s["molecular_recognition"]]
    print(f"  Sites with molecular recognition: {len(mol_recog)}")
    print("Done.")

if __name__ == "__main__":
    main()
