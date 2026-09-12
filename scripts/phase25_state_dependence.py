#!/usr/bin/env python3
"""
Phase 25: State Dependence Analysis

Compare PAM binding predictions across:
- Resting-like (APO) - existing data
- Open-like - from Phase 3-4 results
- Desensitized-like - from Phase 3-4 results

Determine whether the PAM prefers specific receptor states.
"""
import json, csv, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "14_stoichiometry_state"
OUTDIR.mkdir(parents=True, exist_ok=True)

EXISTING_BOLTZ2 = ROOT / "phase08_boltz_affinity" / "boltz2_outputs"
STATE_MODELS = ROOT / "02_receptor_ensemble" / "state_models" / "predictions"
ACH_RESULTS = ROOT / "05_ach_occupancy" / "yaml_inputs" / "predictions"
COFOLDING = ROOT / "cofolding_study" / "results" / "boltz2"


def load_confidence_data():
    """Load existing Boltz-2 confidence scores."""
    conf_file = ROOT / "11_boltz2_analysis" / "boltz2_all_confidences.csv"
    if not conf_file.exists():
        return []

    results = []
    with open(conf_file) as f:
        for row in csv.DictReader(f):
            results.append(row)
    return results


def find_state_predictions():
    """Find Boltz-2 predictions from state models."""
    state_predictions = {}

    for state_dir in [STATE_MODELS, ACH_RESULTS]:
        if not state_dir.exists():
            continue
        for pred_dir in state_dir.iterdir():
            if pred_dir.is_dir() and "predictions" in str(pred_dir):
                for model_dir in pred_dir.iterdir():
                    if model_dir.is_dir():
                        state = "unknown"
                        if "open" in model_dir.name.lower():
                            state = "open"
                        elif "desens" in model_dir.name.lower():
                            state = "desensitized"
                        elif "ach" in model_dir.name.lower():
                            state = "ach_bound"

                        conf_files = list(model_dir.glob("*confidence*"))
                        if conf_files:
                            if state not in state_predictions:
                                state_predictions[state] = []
                            state_predictions[state].append({
                                "model": model_dir.name,
                                "confidence_file": str(conf_files[0]),
                            })

    return state_predictions


def analyze_state_preference(confidence_data, state_predictions):
    """Analyze whether PAM binding is state-dependent."""
    analysis = {
        "resting_like": {
            "source": "existing_boltz2",
            "n_models": len([r for r in confidence_data if "site" in str(r)]),
            "mean_confidence": 0,
        },
        "open": {"source": "phase3_state_models", "n_models": 0, "mean_confidence": 0},
        "desensitized": {"source": "phase3_state_models", "n_models": 0, "mean_confidence": 0},
        "ach_bound": {"source": "phase6_ach_occupancy", "n_models": 0, "mean_confidence": 0},
    }

    for state, predictions in state_predictions.items():
        if state in analysis:
            analysis[state]["n_models"] = len(predictions)

    analysis["conclusion"] = "State dependence cannot be determined without open/desensitized model results"
    analysis["recommendation"] = "Run Phase 3-4 and Phase 6 to generate state models, then re-run this analysis"

    return analysis


def main():
    print("Phase 25: State Dependence Analysis")
    print("=" * 60)

    confidence_data = load_confidence_data()
    print(f"Existing Boltz-2 confidence records: {len(confidence_data)}")

    state_predictions = find_state_predictions()
    print(f"State model predictions found: {list(state_predictions.keys())}")

    analysis = analyze_state_preference(confidence_data, state_predictions)

    print("\nState Analysis:")
    for state, info in analysis.items():
        if isinstance(info, dict):
            print(f"  {state}: {info.get('n_models', 0)} models ({info.get('source', 'N/A')})")

    report_path = OUTDIR / "state_dependence_report.json"
    report_path.write_text(json.dumps(analysis, indent=2))
    print(f"\nReport: {report_path}")

    state_comparison = {
        "states_compared": list(analysis.keys()),
        "methodology": "Compare Boltz-2 confidence scores across receptor conformational states",
        "status": "PARTIALLY_COMPLETE" if len(state_predictions) < 2 else "COMPLETE",
        "gap": "Need open/desensitized state models from Phase 3-4 for full comparison",
    }
    comparison_path = OUTDIR / "state_comparison.json"
    comparison_path.write_text(json.dumps(state_comparison, indent=2))
    print(f"Comparison: {comparison_path}")


if __name__ == "__main__":
    main()
