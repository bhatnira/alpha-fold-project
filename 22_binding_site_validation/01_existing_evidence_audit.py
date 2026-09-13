#!/usr/bin/env python3
"""Section 1: Existing-evidence audit.

Inventories the completed computational work and classifies each result by
evidence type, what it supports, independence, and known limitations.
"""

import json
from datetime import datetime
from pathlib import Path

from config import RESULTS, AF3_OUTPUTS, BOLTZ2_OUTPUTS, COFOLDING_RESULTS, PIPELINE

RESULTS.mkdir(parents=True, exist_ok=True)
MODULATOR = PIPELINE / "modulator-dataset-a9a10.csv"
MMP = PIPELINE / "10_sar_validation" / "mmp_analysis.json"


def count_files(root, pattern):
    root = Path(root)
    if not root.exists():
        return 0
    return len(list(root.rglob(pattern)))


def count_result_dirs(group):
    gdir = COFOLDING_RESULTS / group
    if not gdir.exists():
        return 0
    n = 0
    for d in gdir.iterdir():
        if d.is_dir() and count_files(d, "confidence_*.json") > 0:
            n += 1
    return n


def mmp_count():
    try:
        with open(MMP) as f:
            d = json.load(f)
        return len(d.get("informative_mmps", []))
    except Exception:
        return 0


def inventory():
    rows = []
    add = rows.append

    add({
        "layer": "AF3 cofolding (existing)",
        "evidence_type": "receptor_structural / ligand_pose",
        "n": count_files(AF3_OUTPUTS, "*_model.cif"),
        "supports": "Receptor architecture; recurring ligand regions across conditions",
        "does_not_support": "Affinity / potency; validated residue-level interface mapping",
        "independence": "Independent method (AlphaFold3)",
        "correlated": "Multiple models of the same run are seeds, not independent observations",
        "caveat": "No experimental bound structure exists to benchmark pose accuracy",
    })

    add({
        "layer": "Boltz-2 cofolding (existing outputs)",
        "evidence_type": "receptor_structural / ligand_pose",
        "n": count_files(BOLTZ2_OUTPUTS, "*.cif"),
        "supports": "Independent cofolding corroboration of regions identified by AF3",
        "does_not_support": "Affinity (old runs used alpha10 monomer input)",
        "independence": "Independent method (Boltz-2)",
        "correlated": "MSA derived from same sequence database",
        "caveat": "Pre-correction runs on alpha10 monomer are INVALID (FALLACIES_AUDIT C7)",
    })

    group_meta = {
        "state_models": ("receptor_structural", "Boltz-2 seed series"),
        "ach_states": ("receptor_structural", "Boltz-2 seed series"),
        "binary": ("receptor_ligand_interaction", "Boltz-2 ligand-conditioned"),
        "ternary": ("receptor_ligand_interaction", "Boltz-2 ligand-conditioned w/ ACh"),
        "full_panel": ("receptor_ligand_interaction", "Boltz-2 ligand-conditioned 30-cpd panel"),
    }
    for group, (kind, indep) in group_meta.items():
        n = count_result_dirs(group)
        caveat = ""
        if group == "full_panel":
            caveat = f"run in progress ({n}/180 complete at audit time)"
        add({
            "layer": f"Boltz-2 site-directed {group}",
            "evidence_type": kind,
            "n": n,
            "supports": "Ligand-conditioned placement at candidate sites; residue contacts; pose reproducibility",
            "does_not_support": "Absolute affinity (no affinity_*.json emitted in this run)",
            "independence": indep,
            "correlated": "Same receptor model family used across compounds",
            "caveat": caveat,
        })

    add({
        "layer": "Molecular docking (existing)",
        "evidence_type": "binding_site / ligand_pose",
        "n": count_files(PIPELINE / "09_docking", "*_docked.pdbqt"),
        "supports": "Orthogonal pocket-occupancy evidence at candidate sites; ligand pose library",
        "does_not_support": "SAR ranking alone (scores are relative, not affinity)",
        "independence": "Independent method (AutoDock family)",
        "correlated": "Docking input receptor models derive from the AF3/Boltz-2 ensemble",
        "caveat": "",
    })

    add({
        "layer": "Interaction fingerprints (existing)",
        "evidence_type": "receptor_ligand_interaction",
        "n": count_files(PIPELINE / "06_interaction_fingerprints", "*.csv"),
        "supports": "Residue contact frequencies from prior analyses",
        "does_not_support": "",
        "independence": "Same structures re-analyzed -> correlated",
        "correlated": "Derived from the same predicted complexes",
        "caveat": "Partly stale; recomputed on fresh full_panel data in step 05",
    })

    n_sar = 30 if MODULATOR.exists() else 0
    add({
        "layer": "Experimental SAR (existing)",
        "evidence_type": "experimental_sar",
        "n": n_sar,
        "supports": "Activity labels, potency, stereochemical series, % potentiation",
        "does_not_support": "Binding-site location (no spatial information in SAR alone)",
        "independence": "Independent (experimental)",
        "correlated": "n/a",
        "caveat": "30-compound panel: 7 active / 23 inactive",
    })

    add({
        "layer": "Boltz-2 affinity prediction",
        "evidence_type": "affinity",
        "n": count_files(COFOLDING_RESULTS, "affinity_*.json"),
        "supports": "Nothing in the current run (no affinity_*.json emitted)",
        "does_not_support": "Affinity ranking",
        "independence": "n/a",
        "correlated": "n/a",
        "caveat": "Earlier affinity 'proxy' matrix was circular; affinity layer is ABSENT until a corrected run exists",
    })

    add({
        "layer": "XAI",
        "evidence_type": "explainability",
        "n": 0,
        "supports": "No valid XAI feature importances exist",
        "does_not_support": "Ligand-fragment attributions",
        "independence": "n/a",
        "correlated": "n/a",
        "caveat": "Prior XAI permutation importances were computed on an overfit model; not usable (FALLACIES_AUDIT)",
    })

    add({
        "layer": "Activity-cliff / MMP (existing)",
        "evidence_type": "activity_cliff",
        "n": mmp_count(),
        "supports": "Matched-pair SAR rules: 12vs7 (alkynyl), 25vs1 (Br), 18vs19 (benzyl vs vinyl), 3vs1",
        "does_not_support": "Structural mechanism of the cliff",
        "independence": "Independent (experimental + computational)",
        "correlated": "n/a",
        "caveat": "Used as an SAR challenge layer, not as structural evidence",
    })

    return rows


def main():
    rows = inventory()
    out = {
        "audit_time": datetime.now().isoformat(),
        "summary_note": "Do not count multiple analyses of the same structure as independent evidence.",
        "layers": rows,
    }
    out_path = RESULTS / "01_evidence_audit.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print("Evidence audit rows:", len(rows))
    for r in rows:
        tag = "[CORRELATED]" if r["correlated"] else ""
        if r["caveat"]:
            tag += " [CAVEAT]"
        print(f"  {r['layer'][:52]:52s} n={r['n']:5d}  {r['evidence_type']:<38s} {tag}")
    print("Wrote:", out_path)


if __name__ == "__main__":
    main()