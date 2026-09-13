#!/usr/bin/env python3
"""3. SiteAF3 input preparation (site-conditioned candidate-site testing).

Reuses the validated SiteAF3 JSON template from phase2_siteaf3_validation:
sites x stoichiometries x representative compound panel. Produces a
self-contained input manifest in this package.

Output: outputs/03_siteaf3_manifest.json
"""

import json
import shutil
from pathlib import Path

from config import (
    RESULTS,
    CHAIN_MAP,
    STOICHIOMETRIES,
    SITE_INFO,
    CANDIDATE_SITES,
    REPRESENTATIVE_COMPOUNDS,
    SITE_AF3_INPUTS,
    COFOLDING,
)

RECEPTOR_PDB = {
    "2to3": Path("/cluster/home/nbhatt04/lean_pipeline/09_docking/receptor_2to3.pdb"),
    "3to2": Path("/cluster/home/nbhatt04/lean_pipeline/09_docking/receptor_3to2.pdb"),
}
HOTSPOT_PB = Path("/cluster/home/nbhatt04/lean_pipeline/phase2_siteaf3_validation/hotspot_pocket_pdbs")
PHASE2_INPUTS = Path("/cluster/home/nbhatt04/lean_pipeline/phase2_siteaf3_validation/siteaf3_inputs")

SEEDS = [42, 123, 456]


def site_fixed_chains(stoich):
    """Fixed (non-diffused) receptor chains for a stoichiometry, order preserved."""
    ordered = CHAIN_MAP[stoich]["alpha9"] + CHAIN_MAP[stoich]["alpha10"]
    return ordered


def validate_receptor(pdb, expected_chains, stoich):
    chains = set()
    with open(pdb) as f:
        for line in f:
            if line.startswith("ATOM") and line[21].strip():
                chains.add(line[21].strip())
    missing = set(expected_chains) - chains
    return {
        "path": str(pdb),
        "present_chains": sorted(chains),
        "missing_expected": sorted(missing),
        "ok": not missing,
    }


def build_config(name, stoich, site, compound_key, cmpd):
    return {
        "name": f"siteaf3_{name}",
        "receptor": [
            {
                "rec_struct_path": str(RECEPTOR_PDB[stoich]),
                "fixed_chain_id": site_fixed_chains(stoich),
                "hotspot_path": str(HOTSPOT_PB / f"site{site}_{stoich}_hotspot.pdb"),
                "pocket_path": str(HOTSPOT_PB / f"site{site}_{stoich}_pocket.pdb"),
            }
        ],
        "ligand": [
            {
                "small_molecule": {
                    "id": "LIG",
                    "smiles": REPRESENTATIVE_COMPOUNDS[compound_key]["smiles"],
                }
            }
        ],
        "modelSeeds": SEEDS,
    }


def main():
    SITE_AF3_INPUTS.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)

    manifest = {
        "method": "SiteAF3 site-conditioned modeling (Tang & Wang, PNAS 10.1073/pnas.2521048122)",
        "seeds": SEEDS,
        "candidate_sites": {},
        "configs": [],
        "warnings": [],
    }

    n = 0
    for site in CANDIDATE_SITES:
        manifest["candidate_sites"][site] = SITE_INFO[site]
        for stoich in STOICHIOMETRIES:
            rec = validate_receptor(RECEPTOR_PDB[stoich], site_fixed_chains(stoich), stoich)
            if not rec["ok"]:
                manifest["warnings"].append(f"receptor {stoich} missing chains {rec['missing_expected']}")
            for compound_key, cmpd in REPRESENTATIVE_COMPOUNDS.items():
                name = f"site{site}_{stoich}_{compound_key}"
                cfg = build_config(name, stoich, site, compound_key, cmpd)
                out = SITE_AF3_INPUTS / f"{name}.json"
                with open(out, "w") as f:
                    json.dump(cfg, f, indent=2)
                manifest["configs"].append(str(out))
                n += 1

    # Reference the original phase2 templates for traceability.
    phase2 = sorted(PHASE2_INPUTS.glob("*_CPD*.json")) if PHASE2_INPUTS.exists() else []
    manifest["phase2_template_count"] = len(phase2)
    manifest["phase2_available"] = PHASE2_INPUTS.exists()

    with open(RESULTS / "03_siteaf3_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"[ok] wrote {n} SiteAF3 configs to {SITE_AF3_INPUTS}")
    print(f"[ok] manifest -> {RESULTS / '03_siteaf3_manifest.json'}")
    if manifest["warnings"]:
        for w in manifest["warnings"]:
            print(f"[warn] {w}")


if __name__ == "__main__":
    main()