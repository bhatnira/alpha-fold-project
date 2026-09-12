#!/usr/bin/env python3
"""
Answer two questions:
  1. Do all active ascorbic acid analogs bind at the same site?
  2. Does the binding interaction correlate with activity (active vs inactive)?

Reads Boltz-2 cofolding results (confidence, affinity, contacts).
"""

import json
import csv
from pathlib import Path
from collections import defaultdict

PIPELINE = Path("/cluster/home/nbhatt04/lean_pipeline")
COFOLDING = PIPELINE / "site_directed_cofolding" / "results"
EXISTING_BOLTZ2 = PIPELINE / "data" / "allostery" / "boltz2" / "outputs"
AF3_OUTPUTS = PIPELINE / "data" / "allostery" / "af3_outputs"

# Activity data
COMPOUNDS = {
    1:  {"name": "CPD1",  "uM": 1797,  "pot": 286,  "class": "active"},
    2:  {"name": "CPD2",  "uM": 1316,  "pot": 300,  "class": "active"},
    3:  {"name": "CPD3",  "uM": 6077,  "pot": 293,  "class": "active"},
    12: {"name": "CPD12", "uM": 0.198, "pot": 150,  "class": "active"},
    18: {"name": "CPD18", "uM": 1288,  "pot": 500,  "class": "active"},
    24: {"name": "CPD24", "uM": 1202,  "pot": 190,  "class": "active"},
    25: {"name": "CPD25", "uM": 2.63,  "pot": 180,  "class": "active"},
    4:  {"name": "CPD4",  "uM": 0,     "pot": 0,    "class": "inactive"},
    5:  {"name": "CPD5",  "uM": 0,     "pot": 0,    "class": "inactive"},
    6:  {"name": "CPD6",  "uM": 0,     "pot": 0,    "class": "inactive"},
    7:  {"name": "CPD7",  "uM": 0,     "pot": 0,    "class": "inactive"},
    8:  {"name": "CPD8",  "uM": 0,     "pot": 0,    "class": "inactive"},
    9:  {"name": "CPD9",  "uM": 0,     "pot": 0,    "class": "inactive"},
    10: {"name": "CPD10", "uM": 0,     "pot": 0,    "class": "inactive"},
    11: {"name": "CPD11", "uM": 0,     "pot": 0,    "class": "inactive"},
    13: {"name": "CPD13", "uM": 0,     "pot": 0,    "class": "inactive"},
    14: {"name": "CPD14", "uM": 0,     "pot": 0,    "class": "inactive"},
    15: {"name": "CPD15", "uM": 0,     "pot": 0,    "class": "inactive"},
    16: {"name": "CPD16", "uM": 0,     "pot": 0,    "class": "inactive"},
    17: {"name": "CPD17", "uM": 0,     "pot": 0,    "class": "inactive"},
    19: {"name": "CPD19", "uM": 0,     "pot": 0,    "class": "inactive"},
    20: {"name": "CPD20", "uM": 0,     "pot": 0,    "class": "inactive"},
    21: {"name": "CPD21", "uM": 0,     "pot": 0,    "class": "inactive"},
    22: {"name": "CPD22", "uM": 0,     "pot": 0,    "class": "inactive"},
    23: {"name": "CPD23", "uM": 0,     "pot": 0,    "class": "inactive"},
    26: {"name": "CPD26", "uM": 0,     "pot": 0,    "class": "inactive"},
    27: {"name": "CPD27", "uM": 0,     "pot": 0,    "class": "inactive"},
    28: {"name": "CPD28", "uM": 0,     "pot": 0,    "class": "inactive"},
    29: {"name": "CPD29", "uM": 0,     "pot": 0,    "class": "inactive"},
    30: {"name": "CPD30", "uM": 0,     "pot": 0,    "class": "inactive"},
}


def parse_cofolding_results():
    """Parse all Boltz-2 cofolding results."""
    records = []
    for complex_type in ["binary", "ternary", "full_panel"]:
        type_dir = COFOLDING / complex_type
        if not type_dir.exists():
            continue
        for result_dir in type_dir.iterdir():
            if not result_dir.is_dir():
                continue
            name = result_dir.name

            # Find boltz_results subdirectory
            boltz_dirs = list(result_dir.glob("boltz_results_*"))
            if boltz_dirs:
                pred_dir = boltz_dirs[0] / "predictions"
            else:
                pred_dir = result_dir

            # Parse confidence
            conf = {}
            for p in pred_dir.rglob("confidence_*.json"):
                try:
                    with open(p) as f:
                        conf = json.load(f)
                    break
                except Exception:
                    continue

            # Parse affinity
            aff = {}
            for p in pred_dir.rglob("affinity_*.json"):
                try:
                    with open(p) as f:
                        aff = json.load(f)
                    break
                except Exception:
                    continue

            # Parse contacts from PDB
            contacts = []
            for pdb in pred_dir.rglob("*_model_0.pdb"):
                try:
                    with open(pdb) as f:
                        for line in f:
                            if line.startswith("ATOM"):
                                resname = line[17:20].strip()
                                resnum = int(line[22:26].strip())
                                chain = line[21].strip()
                                if resname not in ["HOH", "WAT"]:
                                    contacts.append(f"{chain}:{resname}{resnum}")
                except Exception:
                    continue

            if conf or aff:
                # Extract compound ID from name
                cpd_id = None
                for cid in COMPOUNDS:
                    if f"cpd{cid}" in name.lower() or f"_cpd{cid}_" in name:
                        cpd_id = cid
                        break

                records.append({
                    "name": name,
                    "complex_type": complex_type,
                    "compound_id": cpd_id,
                    "confidence": conf.get("confidence_score"),
                    "ptm": conf.get("ptm"),
                    "iptm": conf.get("iptm"),
                    "ligand_iptm": conf.get("ligand_iptm"),
                    "complex_plddt": conf.get("complex_plddt"),
                    "affinity_pred": aff.get("affinity_pred_value"),
                    "affinity_prob": aff.get("affinity_probability_binary"),
                    "n_contacts": len(set(contacts)),
                })
    return records


def parse_existing_boltz2():
    """Parse existing Boltz-2 blind cofolding results from data/allostery."""
    records = []
    boltz_dir = EXISTING_BOLTZ2
    if not boltz_dir.exists():
        return records

    for result_dir in boltz_dir.iterdir():
        if not result_dir.is_dir():
            continue
        name = result_dir.name

        # Parse confidence
        conf = {}
        for p in result_dir.rglob("confidence_*.json"):
            try:
                with open(p) as f:
                    conf = json.load(f)
                break
            except Exception:
                continue

        # Parse ligand contacts from CIF
        ligand_contacts = set()
        for cif in result_dir.rglob("summary_confidences.cif"):
            try:
                with open(cif) as f:
                    for line in f:
                        if not line.startswith("#") and not line.startswith("data_") and not line.startswith("_"):
                            parts = line.split()
                            if len(parts) > 3 and parts[1] == "ligand":
                                ligand_contacts.add(parts[2])
            except Exception:
                continue

        if conf:
            # Determine ligand type
            ligand = "unknown"
            for lg in ["LASC", "DASC", "ACETATE", "OETHYL", "RYANODINE", "DERIVATIVE"]:
                if lg in name.upper():
                    ligand = lg
                    break

            # Determine stoichiometry
            stoich = "unknown"
            if "2to3" in name:
                stoich = "2to3"
            elif "3to2" in name:
                stoich = "3to2"

            records.append({
                "name": name,
                "ligand": ligand,
                "stoichiometry": stoich,
                "confidence": conf.get("confidence_score"),
                "ptm": conf.get("ptm"),
                "iptm": conf.get("iptm"),
                "ligand_iptm": conf.get("ligand_iptm"),
                "complex_plddt": conf.get("complex_plddt"),
                "n_ligand_contacts": len(ligand_contacts),
            })
    return records


def question1_same_site(records):
    """
    Q1: Do all active analogs bind at the same site?

    Uses existing Boltz-2 blind cofolding results.
    If all active compounds show high ligand_iptm at the same site,
    that's evidence for a shared binding site.
    """
    print("=" * 70)
    print("Q1: DO ALL ACTIVE ANALOGS BIND AT THE SAME SITE?")
    print("=" * 70)

    # Use existing Boltz-2 results (blind cofolding)
    existing = parse_existing_boltz2()
    if not existing:
        print("\n  No existing Boltz-2 results found.")
        print("  Using cofolding study results instead...")
        existing = records

    # Map ligands to activity
    ligand_activity = {
        "LASC": "active",      # L-ascorbate
        "DASC": "active",      # D-ascorbate (weak)
        "ACETATE": "inactive", # acetate control
        "OETHYL": "active",    # O-ethyl ascorbate
        "RYANODINE": "active", # ryanodine
        "DERIVATIVE": "active", # derivative experimental
    }

    # Group by ligand and stoichiometry
    by_ligand = defaultdict(list)
    for r in existing:
        lig = r.get("ligand", "unknown")
        by_ligand[lig].append(r)

    print("\n  LIGAND BINDING SUMMARY (from existing Boltz-2 blind cofolding):")
    print(f"  {'Ligand':<18s} {'Activity':<10s} {'n_models':<10s} {'avg_conf':<10s} {'avg_iptm':<10s} {'avg_plddt':<10s}")
    print("  " + "-" * 70)

    active_iptms = []
    inactive_iptms = []

    for lig in ["LASC", "DASC", "OETHYL", "RYANODINE", "DERIVATIVE", "ACETATE"]:
        recs = by_ligand.get(lig, [])
        activity = ligand_activity.get(lig, "unknown")
        if recs:
            confs = [r["confidence"] for r in recs if r.get("confidence")]
            iptms = [r["ligand_iptm"] for r in recs if r.get("ligand_iptm")]
            plddts = [r["complex_plddt"] for r in recs if r.get("complex_plddt")]

            avg_c = sum(confs) / len(confs) if confs else 0
            avg_i = sum(iptms) / len(iptms) if iptms else 0
            avg_p = sum(plddts) / len(plddts) if plddts else 0

            print(f"  {lig:<18s} {activity:<10s} {len(recs):<10d} {avg_c:<10.3f} {avg_i:<10.3f} {avg_p:<10.3f}")

            if activity == "active":
                active_iptms.extend(iptms)
            elif activity == "inactive":
                inactive_iptms.extend(iptms)

    # Answer Q1
    print("\n  ANALYSIS:")
    if active_iptms and inactive_iptms:
        avg_active = sum(active_iptms) / len(active_iptms)
        avg_inactive = sum(inactive_iptms) / len(inactive_iptms)
        print(f"    Active compounds avg ligand_iptm:  {avg_active:.3f} (n={len(active_iptms)})")
        print(f"    Inactive compounds avg ligand_iptm: {avg_inactive:.3f} (n={len(inactive_iptms)})")
        if avg_active > avg_inactive:
            print(f"    → Active compounds show HIGHER ligand_iptm ({avg_active:.3f} vs {avg_inactive:.3f})")
            print(f"    → This suggests active analogs bind more favorably at the predicted site")
        else:
            print(f"    → No clear difference in ligand_iptm between active and inactive")
    else:
        print("    Insufficient data to compare active vs inactive")


def question2_correlation(records):
    """
    Q2: Does binding interaction correlate with activity?

    Uses cofolding study results (30 compounds).
    Checks if confidence/iptm/affinity correlate with experimental potency.
    """
    print("\n" + "=" * 70)
    print("Q2: DOES BINDING INTERACTION CORRELATE WITH ACTIVITY?")
    print("=" * 70)

    if not records:
        print("\n  No cofolding records found.")
        print("  Waiting for job 3332533 to complete...")
        return

    # Group by compound
    by_compound = defaultdict(list)
    for r in records:
        cid = r.get("compound_id")
        if cid:
            by_compound[cid].append(r)

    print("\n  COMPOUND-LEVEL RESULTS:")
    print(f"  {'Cpd':<6s} {'Activity':<10s} {'Potency':<10s} {'n':<5s} {'avg_conf':<10s} {'avg_iptm':<10s} {'avg_aff':<10s}")
    print("  " + "-" * 65)

    active_data = []
    inactive_data = []

    for cid in sorted(COMPOUNDS.keys()):
        cpd = COMPOUNDS[cid]
        recs = by_compound.get(cid, [])
        if recs:
            confs = [r["confidence"] for r in recs if r.get("confidence")]
            iptms = [r["ligand_iptm"] for r in recs if r.get("ligand_iptm")]
            affs = [r["affinity_pred"] for r in recs if r.get("affinity_pred") is not None]

            avg_c = sum(confs) / len(confs) if confs else 0
            avg_i = sum(iptms) / len(iptms) if iptms else 0
            avg_a = sum(affs) / len(affs) if affs else 0

            pot_str = f"{cpd['uM']:.1f}" if cpd['uM'] > 0 else "inactive"
            print(f"  {cid:<6d} {cpd['class']:<10s} {pot_str:<10s} {len(recs):<5d} {avg_c:<10.3f} {avg_i:<10.3f} {avg_a:<10.3f}")

            if cpd["class"] == "active" and iptms:
                active_data.append({"id": cid, "uM": cpd["uM"], "iptm": avg_i, "conf": avg_c})
            elif cpd["class"] == "inactive" and iptms:
                inactive_data.append({"id": cid, "iptm": avg_i, "conf": avg_c})

    # Statistical comparison
    print("\n  ACTIVE vs INACTIVE COMPARISON:")
    if active_data and inactive_data:
        a_iptm = [d["iptm"] for d in active_data]
        i_iptm = [d["iptm"] for d in inactive_data]
        a_conf = [d["conf"] for d in active_data]
        i_conf = [d["conf"] for d in inactive_data]

        print(f"    Active (n={len(active_data)}):   avg_iptm={sum(a_iptm)/len(a_iptm):.3f}  avg_conf={sum(a_conf)/len(a_conf):.3f}")
        print(f"    Inactive (n={len(inactive_data)}): avg_iptm={sum(i_iptm)/len(i_iptm):.3f}  avg_conf={sum(i_conf)/len(i_conf):.3f}")

        # Potency correlation for actives
        if len(active_data) > 2:
            import math
            # Use log potency
            log_uMs = [math.log10(d["uM"]) for d in active_data if d["uM"] > 0]
            iptms_for_corr = [d["iptm"] for d in active_data if d["uM"] > 0]

            if len(log_uMs) > 2:
                # Simple Pearson correlation
                n = len(log_uMs)
                mean_x = sum(log_uMs) / n
                mean_y = sum(iptms_for_corr) / n
                cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(log_uMs, iptms_for_corr)) / n
                std_x = (sum((x - mean_x) ** 2 for x in log_uMs) / n) ** 0.5
                std_y = (sum((y - mean_y) ** 2 for y in iptms_for_corr) / n) ** 0.5
                r = cov / (std_x * std_y) if std_x * std_y > 0 else 0
                print(f"\n    Potency vs ligand_iptm (actives):")
                print(f"      Pearson r = {r:.3f}")
                if r > 0.3:
                    print(f"      → POSITIVE correlation (higher potency = stronger predicted binding)")
                elif r < -0.3:
                    print(f"      → NEGATIVE correlation (higher potency = weaker predicted binding — unexpected)")
                else:
                    print(f"      → NO clear correlation")
    else:
        print("    Insufficient data (need both active and inactive results)")


def main():
    print("BOLTZ-2 COFOLDING ANALYSIS")
    print("Two questions:")
    print("  Q1: Do all active analogs bind at the same site?")
    print("  Q2: Does binding interaction correlate with activity?")
    print()

    # Parse results
    records = parse_cofolding_results()
    print(f"Parsed {len(records)} cofolding results")

    # Q1
    question1_same_site(records)

    # Q2
    question2_correlation(records)

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
