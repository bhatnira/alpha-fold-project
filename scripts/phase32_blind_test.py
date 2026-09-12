#!/usr/bin/env python3
"""
Phase 32: Blind Prospective Test

Complete the blind prospective test directory with:
- Frozen predictions (already in 17_prospective_design/blind_predictions.json)
- Protocol documentation
- Template for experimental comparison
"""
import json, csv, os
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "19_blind_test"
OUTDIR.mkdir(parents=True, exist_ok=True)

BLIND_PREDICTIONS = ROOT / "17_prospective_design" / "blind_predictions.json"
PROSPECTIVE_RANKING = ROOT / "17_prospective_design" / "prospective_ranking.json"
FROZEN_MODEL = ROOT / "16_model_freeze" / "frozen_model.json"


def load_blind_predictions():
    if not BLIND_PREDICTIONS.exists():
        return {}
    with open(BLIND_PREDICTIONS) as f:
        return json.load(f)


def load_frozen_model():
    if not FROZEN_MODEL.exists():
        return {}
    with open(FROZEN_MODEL) as f:
        return json.load(f)


def create_test_protocol():
    return {
        "protocol_version": "1.0",
        "created": datetime.now().isoformat(),
        "objective": "Prospective experimental validation of computational predictions",
        "steps": [
            {
                "step": 1,
                "action": "Select top-ranked compounds from screening",
                "status": "COMPLETE",
                "note": "Compounds ranked in phase31_prospective_screening",
            },
            {
                "step": 2,
                "action": "Freeze all computational predictions before experimental testing",
                "status": "COMPLETE",
                "note": "Predictions frozen in blind_predictions.json",
            },
            {
                "step": 3,
                "action": "Experimental electrophysiology testing",
                "status": "PENDING",
                "note": "Requires wet lab work",
            },
            {
                "step": 4,
                "action": "Compare predicted vs experimental activity",
                "status": "PENDING",
                "note": "Will be completed after experimental results",
            },
            {
                "step": 5,
                "action": "Report prospective validation success rate",
                "status": "PENDING",
                "note": "Final validation metric",
            },
        ],
        "predictions_frozen_date": "2026-09-05",
        "prediction_summary": {
            "n_compounds_predicted": 3,
            "predicted_active": ["MOL_0011", "MOL_0005", "MOL_0012"],
            "predicted_site": "Site 23",
            "predicted_interface": "alpha9(+)/alpha10(-) ECD",
            "predicted_stoichiometry": "2to3",
            "predicted_state": "resting-like",
            "predicted_ach_dependence": "ACh-dependent potentiation",
        },
    }


def create_comparison_template(predictions):
    compounds = predictions.get("predicted_active", [])
    template = []
    for cpd in compounds:
        template.append({
            "compound_id": cpd,
            "predicted_activity": "ACTIVE",
            "predicted_site": "Site 23",
            "predicted_interface": "alpha9(+)/alpha10(-) ECD",
            "predicted_potency_rank": "Top 3",
            "predicted_interactions": "H-bonds with alpha9:176, alpha10:143; hydrophobic with alpha9:175",
            "experimental_activity": "TBD",
            "experimental_potency": "TBD",
            "experimental potentiation": "TBD",
            "match": "TBD",
            "notes": "",
        })
    return template


def main():
    print("Phase 32: Blind Prospective Test")
    print("=" * 60)

    blind_preds = load_blind_predictions()
    frozen_model = load_frozen_model()

    protocol = create_test_protocol()
    protocol_path = OUTDIR / "blind_test_protocol.json"
    protocol_path.write_text(json.dumps(protocol, indent=2))
    print(f"Protocol: {protocol_path}")

    template = create_comparison_template(blind_preds)
    csv_path = OUTDIR / "blind_predictions_frozen.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(template[0].keys()) if template else [])
        writer.writeheader()
        for row in template:
            writer.writerow(row)
    print(f"Predictions: {csv_path}")

    comparison_path = OUTDIR / "experimental_comparison_template.csv"
    with open(comparison_path, "w", newline="") as f:
        fieldnames = [
            "compound_id", "predicted_activity", "experimental_activity",
            "predicted_potency_rank", "experimental_potency",
            "predicted_site_correct", "match", "validation_status"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in template:
            writer.writerow({
                "compound_id": row["compound_id"],
                "predicted_activity": row["predicted_activity"],
                "experimental_activity": "TBD",
                "predicted_potency_rank": row["predicted_potency_rank"],
                "experimental_potency": "TBD",
                "predicted_site_correct": "TBD",
                "match": "TBD",
                "validation_status": "AWAITING_EXPERIMENT",
            })
    print(f"Comparison template: {comparison_path}")

    report = {
        "status": "PREDICTIONS_FROZEN",
        "n_compounds": len(template),
        "predicted_active": blind_preds.get("predicted_active", []),
        "frozen_date": "2026-09-05",
        "awaiting": "Experimental electrophysiology validation",
        "protocol": str(protocol_path),
        "comparison_template": str(comparison_path),
    }
    report_path = OUTDIR / "blind_test_report.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
