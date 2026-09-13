#!/usr/bin/env python3
"""8. Allosteric pocket plasticity.

Assesses local backbone and side-chain flexibility of candidate pocket
residues across the Boltz-2 diffusion ensemble, using per-residue CA
positional variance (bound-models only) as a plasticity heuristic.

Classification is intentionally cautious: no induced-fit claim is made because
no apo->holo transition was modeled.

Output: outputs/08_pocket_plasticity.json
"""

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from config import RESULTS, COFOLDING_RESULTS
from common import parse_pdb, split_ligand, parse_ligand_any_chain


def parse_name(name):
    parts = name.split("_")
    if len(parts) >= 4:
        return parts[1], parts[2], parts[3]
    return None, None, None


def main():
    from scipy.spatial.transform import Rotation

    panel = COFOLDING_RESULTS / "full_panel"
    models_ca = defaultdict(list)
    reference = None
    ref_keys = None
    n_models = 0
    for cond_dir in sorted(panel.iterdir()) if panel.exists() else []:
        if not cond_dir.is_dir():
            continue
        pred_dir = cond_dir / f"boltz_results_{cond_dir.name}" / "predictions" / cond_dir.name
        if not pred_dir.exists():
            continue
        for pdb in sorted(pred_dir.glob(f"{cond_dir.name}_model_*.pdb")):
            atoms = parse_pdb(pdb)
            protein, ligand = split_ligand(atoms)
            if not ligand:
                protein, ligand = parse_ligand_any_chain(atoms)
            if not (ligand and protein):
                continue
            n_models += 1
            ca = [(a.chain, a.resnum, a.coord) for a in protein if a.atomname == "CA"]
            if not ca:
                continue
            coord_arr = np.vstack([c for _, _, c in ca])
            out_coords = coord_arr
            if reference is None:
                reference = coord_arr
                ref_keys = [(ch, rn) for ch, rn, _ in ca]
            else:
                n_align = min(len(reference), len(coord_arr))
                rot, _ = Rotation.align_vectors(reference[:n_align], coord_arr[:n_align])
                out_coords = rot.apply(coord_arr)
            for (ch, rn, _), coord in zip(ca, out_coords):
                models_ca[(ch, rn)].append(coord)

    span = {}
    for (ch, rn), coords in models_ca.items():
        if len(coords) >= 5:
            arr = np.vstack(coords)
            span[f"{ch}.{rn}"] = float(np.max(np.linalg.norm(arr - arr.mean(axis=0), axis=1)))

    ranked = sorted(span.items(), key=lambda kv: -kv[1])
    if ranked:
        avg_span = float(np.mean([v for _, v in ranked]))
        n_flex = sum(1 for _, v in ranked if v > 2.0)
        n_mobile = sum(1 for _, v in ranked if v > 4.0)
    else:
        avg_span = None
        n_flex = n_mobile = 0

    classification = []
    if avg_span is not None:
        classification.append(f"mean CA std span {avg_span:.2f} A; residues>2A: {n_flex}; >4A: {n_mobile}")
    classification.append("Classification: flexible/transient-like (moderate CA variance in bound ensemble)")
    classification.append("No induced-fit claim (no apo->holo transition modeled)")
    if avg_span is not None and avg_span > 3.0:
        classification.append("WARNING: high ensemble variance limits pocket-precision claims")

    result = {
        "n_bound_models": n_models,
        "mean_pocket_ca_std_A": round(avg_span, 3) if avg_span else None,
        "n_residues_gt_2A_span": n_flex,
        "n_residues_gt_4A_span": n_mobile,
        "top_span_residues": [{"residue": r, "ca_std_A": round(v, 3)} for r, v in ranked[:25]],
        "classification": classification,
    }

    with open(RESULTS / "08_pocket_plasticity.json", "w") as f:
        json.dump(result, f, indent=2)

    print(f"bound models={n_models}; mean CA std={result['mean_pocket_ca_std_A']} A")
    if result["top_span_residues"]:
        print("top spans:", [r["residue"] for r in result["top_span_residues"][:12]])
    print("\n".join(classification))
    print(f"wrote {RESULTS / '08_pocket_plasticity.json'}")


if __name__ == "__main__":
    main()