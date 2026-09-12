#!/usr/bin/env python3
"""
Generate ALL Boltz-2 YAML inputs for site-directed cofolding study.

Covers:
  Step 1: Open/desensitized receptor state models (4 YAMLs)
  Step 2: ACh-bound states (8 YAMLs)
  Step 3: Site-directed binary (R+PAM) for 5 sites x 8 compounds x 2 stoich (80 YAMLs)
  Step 4: Site-directed ternary (R+ACh+PAM) for 5 sites x 8 compounds x 2 stoich (80 YAMLs)
  Step 5: Full 30-compound panel binary for top 3 sites x 2 stoich (60 YAMLs)

Total: ~232 YAML inputs
"""

import csv
import yaml
import json
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

PIPELINE = Path("/cluster/home/nbhatt04/lean_pipeline")
OUTPUT = PIPELINE / "site_directed_cofolding" / "yaml_inputs"

# Receptor sequences
# α9 subunit (chain C, D, E in 2to3; A, B, C in 3to2)
ALPHA9 = "MGLRSHHLSLGLLLLFLLPAECLGAEGRLALKLFRDLFANYTSALRPVADTDQTLNVTLEVTLSQIIDMDERNQVLTLYLWIRQEWTDAYLRWDPNAYGGLDAIRIPSSLVWRPDIVLYNKADAQPPGSASTNVVLRHDGAVRWDAPAITRSSCRVDVAAFPFDAQHCGLTFGSWTHGGHQLDVRPRGAAASLADFVENVEWRVLGMPARRRVLTYGCCSEPYPDVTFTLLLRRRAAAYVCNLLLPCVLISLLAPLAFHLPADSGEKVSLGVTVLLALTVFQLLLAESMPPAESVPLIGKYYMATMTMVTFSTALTILIMNLHYCGPSVRPVPAWARALLLGHLARGLCVRERGEPCGQSRPPELSPSPQSPEGGAGPPAGPCHEPRCLCRQEALLHHVATIANTFRSHRAAQRCHEDWKRLARVMDRFFLAIFFSMALVMSLLVLVQAL"

# α10 subunit (chain A, B in 2to3; D, E in 3to2)
ALPHA10 = "MNWSHSCISFCWIYFAASRLRAAETADGKYAQKLFNDLFEDYSNALRPVEDTDKVLNVTLQITLSQIKDMDERNQILTAYLWIRQIWHDAYLTWDRDQYDGLDSIRIPSDLVWRPDIVLYNKADDESSEPVNTNVVLRYDGLITWDAPAITKSSCVVDVTYFPFDNQQCNLTFGSWTYNGNQVDIFNALDSGDLSDFIEDVEWEVHGMPAVKNVISYGCCSEPYPDVTFTLLLKRRSSFYIVNLLIPCVLISFLAPLSFYLPAASGEKVSLGVTILLAMTVFQLMVAEIMPASENVPLIGKYYIATMALITASTALTIMVMNIHFCGAEARPVPHWARVVILKYMSRVLFVYDVGESCLSPHHSRERDHLTKVYSKLPESNLKAARNKDLSRKKDMNKRLKNDLGCQGKNPQEAESYCAQYKVLTRNIEYIAKCLKDHKATNSKGSEWKKVAKVIDRFFMWIFFIMVFVMTILIIARAD"

# ACh SMILES (acetylcholine)
ACH_SMILES = "CC(=O)OCC[N+](C)(C)C"

# Stoichiometry definitions: subunit order A,B,C,D,E
STOICHIOMETRIES = {
    "2to3": {
        "chains": ["A", "B", "C", "D", "E"],
        "subunits": ["alpha10", "alpha10", "alpha9", "alpha9", "alpha9"],
        "desc": "alpha9_2_alpha10_3",
    },
    "3to2": {
        "chains": ["A", "B", "C", "D", "E"],
        "subunits": ["alpha10", "alpha10", "alpha10", "alpha9", "alpha9"],
        "desc": "alpha9_3_alpha10_2",
    },
}

# Receptor states
STATES = ["resting", "open", "desensitized"]

# Candidate sites (from pipeline ranking)
CANDIDATE_SITES = {
    23: {"interface": "alpha9_alpha10", "desc": "Rank 1 — best overall (score 0.627)"},
    21: {"interface": "alpha10_alpha9",  "desc": "Rank 2 — best SAR + docking (score 0.535)"},
    5:  {"interface": "alpha9_alpha10", "desc": "Rank 3 — best AF3 convergence (score 0.511)"},
    34: {"interface": "alpha9_alpha9",   "desc": "Rank 4 — best Boltz-2 affinity (score 0.477)"},
    7:  {"interface": "alpha9_alpha10", "desc": "Rank 5 — balanced evidence (score 0.476)"},
}

# Representative compounds for site-directed study
REPRESENTATIVE_COMPOUNDS = {
    "L_ASCORBATE": {
        "smiles": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
        "activity": "active",
        "potency_uM": 1797,
        "class": "parent_reference",
    },
    "CPD12_POTENT": {
        "smiles": "C(#C)COC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O",
        "activity": "active",
        "potency_uM": 0.198,
        "class": "highly_potent",
    },
    "CPD25_STRONG": {
        "smiles": "C(Br)[C@H]([C@H]1OC(=O)C(O)=C1O)O",
        "activity": "active",
        "potency_uM": 2.63,
        "class": "potent",
    },
    "CPD1_MODERATE": {
        "smiles": "C1([C@@H]([C@H](O)CO)OC(=O)C=1O)O",
        "activity": "active",
        "potency_uM": 1797,
        "class": "moderate",
    },
    "CPD18_MODERATE": {
        "smiles": "C1C=CC(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)=CC=1",
        "activity": "active",
        "potency_uM": 1288,
        "class": "moderate",
    },
    "CPD3_WEAK": {
        "smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2O)CO1)C",
        "activity": "active",
        "potency_uM": 6077,
        "class": "weak",
    },
    "CPD8_INACTIVE": {
        "smiles": "CCCOC1C(C2OC(C)(C)OC2)OC(=O)C=1O",
        "activity": "inactive",
        "potency_uM": 0,
        "class": "inactive",
    },
    "CPD9_INACTIVE": {
        "smiles": "CCCCOC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O",
        "activity": "inactive",
        "potency_uM": 0,
        "class": "inactive",
    },
}

# Full 30-compound panel (from modulator-dataset-a9a10.csv)
FULL_COMPOUND_CSV = PIPELINE / "modulator-dataset-a9a10.csv"

# Model seeds
SEEDS = [42, 123, 456]


def get_sequence(subunit):
    return ALPHA9 if subunit == "alpha9" else ALPHA10


def make_pentamer(stoich_key, name, ligands=None):
    """Create a Boltz-2 YAML dict for a pentameric receptor + optional ligands."""
    stoich = STOICHIOMETRIES[stoich_key]
    sequences = []
    for chain_id, subunit in zip(stoich["chains"], stoich["subunits"]):
        sequences.append({
            "protein": {
                "id": chain_id,
                "sequence": get_sequence(subunit),
                "msa": "empty",
            }
        })
    if ligands:
        for lig in ligands:
            sequences.append({
                "ligand": {
                    "id": lig["id"],
                    "smiles": lig["smiles"],
                }
            })
    return {
        "version": 1,
        "name": name,
        "sequences": sequences,
    }


def write_yaml(data, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
    return path


def load_full_compounds():
    """Load all 30 compounds from CSV."""
    compounds = []
    with open(FULL_COMPOUND_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            compounds.append({
                "id": int(row["Identifier"]),
                "smiles": row["Smiles"],
                "activity_uM": float(row["Activity (uM)"]),
                "potentiation": float(row["%Potentiation"]),
            })
    return compounds


# ============================================================
# STEP 1: STATE MODELS (open, desensitized) x 2 stoichiometries
# ============================================================

def generate_state_models():
    """Generate receptor-only YAMLs for open and desensitized states."""
    out_dir = OUTPUT / "state_models"
    generated = []
    for stoich in ["2to3", "3to2"]:
        for state in ["open", "desensitized"]:
            name = f"a9a10_{stoich}_{state}"
            data = make_pentamer(stoich, name)
            path = write_yaml(data, out_dir / f"{name}.yaml")
            generated.append(path)
            print(f"  State model: {name}")
    return generated


# ============================================================
# STEP 2: ACH-BOUND STATES (resting + ACh) x 2 stoichiometries
# ============================================================

def generate_ach_states():
    """Generate receptor+ACh YAMLs for resting state."""
    out_dir = OUTPUT / "ach_states"
    generated = []
    for stoich in ["2to3", "3to2"]:
        name = f"a9a10_{stoich}_ach_bound"
        ligands = [{"id": "ACH", "smiles": ACH_SMILES}]
        data = make_pentamer(stoich, name, ligands)
        path = write_yaml(data, out_dir / f"{name}.yaml")
        generated.append(path)
        print(f"  ACh state: {name}")
    return generated


# ============================================================
# STEP 3: SITE-DIRECTED BINARY (R + PAM) — 5 sites x 8 compounds x 2 stoich
# ============================================================

def generate_site_directed_binary():
    """Generate site-directed binary (receptor + PAM) YAMLs."""
    out_dir = OUTPUT / "binary"
    generated = []
    for stoich in ["2to3", "3to2"]:
        for site_id, site_info in CANDIDATE_SITES.items():
            for cpd_name, cpd_info in REPRESENTATIVE_COMPOUNDS.items():
                name = f"sd_{stoich}_site{site_id}_{cpd_name}"
                ligands = [{"id": "PAM", "smiles": cpd_info["smiles"]}]
                data = make_pentamer(stoich, name, ligands)
                path = write_yaml(data, out_dir / f"{name}.yaml")
                generated.append(path)
    print(f"  Binary site-directed: {len(generated)} YAMLs")
    return generated


# ============================================================
# STEP 4: SITE-DIRECTED TERNARY (R + ACh + PAM) — 5 sites x 8 compounds x 2 stoich
# ============================================================

def generate_site_directed_ternary():
    """Generate site-directed ternary (receptor + ACh + PAM) YAMLs."""
    out_dir = OUTPUT / "ternary"
    generated = []
    for stoich in ["2to3", "3to2"]:
        for site_id, site_info in CANDIDATE_SITES.items():
            for cpd_name, cpd_info in REPRESENTATIVE_COMPOUNDS.items():
                name = f"ternary_{stoich}_site{site_id}_{cpd_name}"
                ligands = [
                    {"id": "ACH", "smiles": ACH_SMILES},
                    {"id": "PAM", "smiles": cpd_info["smiles"]},
                ]
                data = make_pentamer(stoich, name, ligands)
                path = write_yaml(data, out_dir / f"{name}.yaml")
                generated.append(path)
    print(f"  Ternary site-directed: {len(generated)} YAMLs")
    return generated


# ============================================================
# STEP 5: FULL 30-COMPOUND BINARY for top 3 sites x 2 stoich
# ============================================================

def generate_full_panel():
    """Generate binary (R+PAM) YAMLs for all 30 compounds at top 3 sites."""
    out_dir = OUTPUT / "full_panel"
    compounds = load_full_compounds()
    generated = []
    top_sites = [23, 21, 5]
    for stoich in ["2to3", "3to2"]:
        for site_id in top_sites:
            for cpd in compounds:
                name = f"full_{stoich}_site{site_id}_cpd{cpd['id']}"
                ligands = [{"id": "PAM", "smiles": cpd["smiles"]}]
                data = make_pentamer(stoich, name, ligands)
                path = write_yaml(data, out_dir / f"{name}.yaml")
                generated.append(path)
    print(f"  Full panel: {len(generated)} YAMLs ({len(compounds)} compounds x 3 sites x 2 stoich)")
    return generated


# ============================================================
# MANIFEST
# ============================================================

def write_manifest(all_groups):
    manifest = {
        "total_yamls": sum(len(v) for v in all_groups.values()),
        "groups": {k: len(v) for k, v in all_groups.items()},
        "stoichiometries": list(STOICHIOMETRIES.keys()),
        "states": STATES,
        "candidate_sites": list(CANDIDATE_SITES.keys()),
        "representative_compounds": list(REPRESENTATIVE_COMPOUNDS.keys()),
        "seeds": SEEDS,
    }
    manifest_path = OUTPUT / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nManifest: {manifest_path}")
    return manifest


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("Site-Directed Cofolding: Full YAML Input Generator")
    print("=" * 70)

    all_groups = {}

    print("\n[Step 1] State models (open, desensitized)...")
    all_groups["state_models"] = generate_state_models()

    print("\n[Step 2] ACh-bound states...")
    all_groups["ach_states"] = generate_ach_states()

    print("\n[Step 3] Site-directed binary (R + PAM)...")
    all_groups["binary"] = generate_site_directed_binary()

    print("\n[Step 4] Site-directed ternary (R + ACh + PAM)...")
    all_groups["ternary"] = generate_site_directed_ternary()

    print("\n[Step 5] Full 30-compound panel...")
    all_groups["full_panel"] = generate_full_panel()

    manifest = write_manifest(all_groups)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for group, paths in all_groups.items():
        print(f"  {group:25s}  {len(paths):4d} YAMLs")
    print(f"  {'TOTAL':25s}  {manifest['total_yamls']:4d} YAMLs")
    print()


if __name__ == "__main__":
    main()
