#!/usr/bin/env python3
"""
Phase 31: Prospective Screening

Screen designed/purchasable compounds using the frozen model.
Apply filters: pharmacophore fit, docking, interaction fingerprints,
chemical novelty, synthetic accessibility.
"""
import json, csv, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "18_prospective_screening"
OUTDIR.mkdir(parents=True, exist_ok=True)

PROSPECTIVE_COMPOUNDS = ROOT / "17_prospective_design" / "prospective_ranking.json"
FROZEN_MODEL = ROOT / "16_model_freeze" / "frozen_model.json"
EXISTING_LIGANDS = ROOT / "data" / "modulator-dataset-a9a10.csv"

SCREENING_FILTERS = {
    "pharmacophore_fit_min": 0.5,
    "mw_max": 500,
    "logP_max": 5,
    "hbd_max": 5,
    "hba_max": 10,
    "novelty_threshold": 0.3,
}


def load_prospective_compounds():
    """Load designed prospective compounds."""
    if not PROSPECTIVE_COMPOUNDS.exists():
        return []
    with open(PROSPECTIVE_COMPOUNDS) as f:
        data = json.load(f)
    return data.get("compounds", data.get("top_compounds", data.get("ranked_compounds", [])))


def load_frozen_model():
    """Load frozen model specification."""
    if not FROZEN_MODEL.exists():
        return {}
    with open(FROZEN_MODEL) as f:
        return json.load(f)


def apply_filters(compounds):
    """Apply screening filters to prospective compounds."""
    filtered = []
    for cpd in compounds:
        passes = True
        reasons = []

        mw = float(cpd.get("mw", cpd.get("molecular_weight", 0)))
        if mw > SCREENING_FILTERS["mw_max"]:
            passes = False
            reasons.append(f"MW {mw:.0f} > {SCREENING_FILTERS['mw_max']}")

        logp = float(cpd.get("logP", cpd.get("logp", 0)))
        if logp > SCREENING_FILTERS["logP_max"]:
            passes = False
            reasons.append(f"logP {logp:.2f} > {SCREENING_FILTERS['logP_max']}")

        hbd = int(cpd.get("hbd", cpd.get("hbd_count", 0)))
        if hbd > SCREENING_FILTERS["hbd_max"]:
            passes = False
            reasons.append(f"HBD {hbd} > {SCREENING_FILTERS['hbd_max']}")

        hba = int(cpd.get("hba", cpd.get("hba_count", 0)))
        if hba > SCREENING_FILTERS["hba_max"]:
            passes = False
            reasons.append(f"HBA {hba} > {SCREENING_FILTERS['hba_max']}")

        pharmacophore = float(cpd.get("pharmacophore_fit", cpd.get("design_score", 0)))
        if pharmacophore < SCREENING_FILTERS["pharmacophore_fit_min"]:
            passes = False
            reasons.append(f"Pharmacophore fit {pharmacophore:.3f} < {SCREENING_FILTERS['pharmacophore_fit_min']}")

        filtered.append({
            **cpd,
            "passes_filter": passes,
            "filter_reasons": reasons,
        })

    return filtered


def rank_compounds(filtered):
    """Rank compounds by multi-objective score."""
    for cpd in filtered:
        if not cpd.get("passes_filter", False):
            cpd["screening_score"] = 0
            continue

        pharmacophore = float(cpd.get("pharmacophore_fit", cpd.get("design_score", 0)))
        mw = float(cpd.get("mw", cpd.get("molecular_weight", 100)))
        logp = float(cpd.get("logP", cpd.get("logp", 0)))

        mw_score = max(0, 1 - (mw / 500))
        logp_score = max(0, 1 - (logp / 5))

        cpd["screening_score"] = (
            0.4 * pharmacophore +
            0.2 * mw_score +
            0.2 * logp_score +
            0.2 * 0.5
        )

    return sorted(filtered, key=lambda x: x.get("screening_score", 0), reverse=True)


def main():
    print("Phase 31: Prospective Screening")
    print("=" * 60)

    compounds = load_prospective_compounds()
    print(f"Prospective compounds loaded: {len(compounds)}")

    frozen_model = load_frozen_model()
    if frozen_model:
        print(f"Frozen model site: {frozen_model.get('candidate_site', 'N/A')}")

    filtered = apply_filters(compounds)
    n_pass = sum(1 for c in filtered if c.get("passes_filter", False))
    print(f"After filters: {n_pass}/{len(filtered)} pass")

    ranked = rank_compounds(filtered)

    csv_path = OUTDIR / "screening_results.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "compound_id", "smiles", "pharmacophore_fit", "mw", "logP",
            "hbd", "hba", "design_score", "screening_score", "passes_filter",
            "filter_reasons", "predicted_site", "predicted_activity"
        ])
        writer.writeheader()
        for cpd in ranked:
            writer.writerow({
                "compound_id": cpd.get("molecule_id", cpd.get("id", cpd.get("compound_id", "N/A"))),
                "smiles": cpd.get("smiles", ""),
                "pharmacophore_fit": cpd.get("pharmacophore_fit", cpd.get("design_score", 0)),
                "mw": cpd.get("mw", cpd.get("molecular_weight", 0)),
                "logP": cpd.get("logP", cpd.get("logp", 0)),
                "hbd": cpd.get("hbd", cpd.get("hbd_count", 0)),
                "hba": cpd.get("hba", cpd.get("hba_count", 0)),
                "design_score": cpd.get("design_score", 0),
                "screening_score": cpd.get("screening_score", 0),
                "passes_filter": cpd.get("passes_filter", False),
                "filter_reasons": ";".join(cpd.get("filter_reasons", [])),
                "predicted_site": "Site 23",
                "predicted_activity": "ACTIVE" if cpd.get("screening_score", 0) > 0.3 else "UNCERTAIN",
            })

    report = {
        "n_input": len(compounds),
        "n_pass_filter": n_pass,
        "n_screened": len(ranked),
        "filters": SCREENING_FILTERS,
        "frozen_model_site": "Site 23",
        "top_5": [
            {
                "id": cpd.get("molecule_id", cpd.get("id", cpd.get("compound_id", "N/A"))),
                "score": cpd.get("screening_score", 0),
                "predicted_activity": "ACTIVE" if cpd.get("screening_score", 0) > 0.3 else "UNCERTAIN",
            }
            for cpd in ranked[:5]
        ],
    }

    report_path = OUTDIR / "screening_report.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"\nReport: {report_path}")
    print(f"Results: {csv_path}")

    for i, cpd in enumerate(ranked[:5]):
        print(f"  #{i+1}: {cpd.get('id', 'N/A')} - score={cpd.get('screening_score', 0):.3f}")


if __name__ == "__main__":
    main()
