#!/usr/bin/env python3
"""Aggregate Boltz-2 affinity prediction results into a clean CSV + ranking."""

import csv, json, sys, yaml
from pathlib import Path

ROOT = Path(__file__).resolve().parent
YAML_DIR = ROOT / "yaml_inputs"
OUT_ROOT = ROOT / "boltz2_outputs"
RESULTS_CSV = ROOT / "boltz2_affinity_results.csv"


def read_inputs():
    parsed = {}
    for p in sorted(YAML_DIR.glob("*.yaml")):
        d = yaml.safe_load(p.read_text())
        stem = p.stem
        parsed[stem] = {
            "site_id": stem.split("_")[0].replace("site", ""),
            "ligand": stem.split("_", 1)[1],
            "smiles": d["sequences"][1]["ligand"]["smiles"],
        }
    return parsed


def main():
    inputs = read_inputs()
    row_template = ["site_id", "ligand", "smiles",
                    "affinity_pred_value", "affinity_probability_binary", "boltz2_success"]
    rows = []
    for stem, meta in sorted(inputs.items()):
        aff = OUT_ROOT.rglob(f"affinity_{stem}.json")
        aff_path = next(iter(aff), None)
        if aff_path is None or not aff_path.exists():
            rows.append({**meta, "affinity_pred_value": None,
                         "affinity_probability_binary": None,
                         "boltz2_success": False})
            continue
        data = json.loads(aff_path.read_text())
        rows.append({
            **meta,
            "affinity_pred_value": data.get("affinity_pred_value"),
            "affinity_probability_binary": data.get("affinity_probability_binary"),
            "boltz2_success": True,
        })

    with RESULTS_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=row_template)
        writer.writeheader()
        writer.writerows(rows)

    n_ok = sum(1 for r in rows if r["boltz2_success"])
    print(f"Wrote {len(rows)} rows to {RESULTS_CSV}")
    print(f"{n_ok}/{len(rows)} predictions successful\n")

    valid = [r for r in rows if r["affinity_pred_value"] is not None]
    valid.sort(key=lambda r: r["affinity_pred_value"])
    print("Ranked by Boltz-2 affinity log-ratio (lower favours binding):")
    print(f"{'site':>5}  {'ligand':20s}  {'affinity_pred':>12}  {'P(binder)':>10}")
    for r in valid:
        print(f"{r['site_id']:>5}  {r['ligand']:20s}  "
              f"{r['affinity_pred_value']:>12.4f}  {r['affinity_probability_binary']:>10.4f}")

    print("\nBest binder per site (most negative affinity_pred_value):")
    best = {}
    for r in valid:
        s = r["site_id"]
        if s not in best or r["affinity_pred_value"] < best[s]["affinity_pred_value"]:
            best[s] = r
    for s in sorted(best):
        r = best[s]
        print(f"  site {s:>2}: {r['ligand']:20s} pred={r['affinity_pred_value']:.4f} P={r['affinity_probability_binary']:.4f}")


if __name__ == "__main__":
    sys.exit(main())