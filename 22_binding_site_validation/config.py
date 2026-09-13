#!/usr/bin/env python3
"""Central configuration for the alpha9alpha10 PAM binding-site validation analysis."""

import csv
from pathlib import Path

PIPELINE = Path("/cluster/home/nbhatt04/lean_pipeline")
PACKAGE = PIPELINE / "22_binding_site_validation"
RESULTS = PACKAGE / "outputs"

PYTHON = "/cluster/scratch/nbhatt04/conda/envs/polymer_ai/bin/python"

COFOLDING = PIPELINE / "site_directed_cofolding"
COFOLDING_RESULTS = COFOLDING / "results"
AF3_OUTPUTS = PIPELINE / "data" / "allostery" / "af3_outputs"
BOLTZ2_OUTPUTS = PIPELINE / "data" / "allostery" / "boltz2" / "outputs"
SITE_AF3_INPUTS = PACKAGE / "siteaf3_inputs"
SITE_AF3_REPO = COFOLDING / "SiteAF3"

MODULATOR_CSV = PIPELINE / "modulator-dataset-a9a10.csv"
NORMALIZED_CSV = PIPELINE / "01_qc" / "normalized_dataset.csv"
MMP_JSON = PIPELINE / "10_sar_validation" / "mmp_analysis.json"
MUTAGENESIS_CSV = PIPELINE / "data" / "MUTAGENESIS_PREDICTIONS.csv"
IFP_CSV = PIPELINE / "06_interaction_fingerprints" / "full_ifp_summary.csv"
SITE_TABLE = PIPELINE / "07_candidate_sites" / "candidate_site_table.json"
COMPETITION_REPORT = PIPELINE / "03_binding_sites" / "competition" / "site_competition_report.md"

CANDIDATE_SITES = [23, 21, 5, 34, 7]

SITE_INFO = {
    23: {"interface": "alpha9(+)/alpha10(-) ECD", "pocket": "vestibular_inter_subunit"},
    21: {"interface": "alpha10(+)/alpha9(-) ECD", "pocket": "vestibular_inter_subunit"},
    5: {"interface": "alpha9(+)/alpha10(-) ECD", "pocket": "vestibular_inter_subunit"},
    34: {"interface": "TMD pore", "pocket": "transmembrane_allosteric"},
    7: {"interface": "alpha9(+)/alpha10(-) ECD", "pocket": "vestibular_inter_subunit"},
}

FULL_PANEL_SITES = [23, 21, 5]

STOICHIOMETRIES = ["2to3", "3to2"]

CHAIN_MAP = {
    "2to3": {"alpha9": ["A", "B"], "alpha10": ["C", "D", "E"]},
    "3to2": {"alpha9": ["A", "B", "C"], "alpha10": ["D", "E"]},
}

CHAIN_SUBUNIT = {
    "2to3": {"A": "alpha9", "B": "alpha9", "C": "alpha10", "D": "alpha10", "E": "alpha10"},
    "3to2": {"A": "alpha9", "B": "alpha9", "C": "alpha9", "D": "alpha10", "E": "alpha10"},
}

ALPHA9_LENGTH = 479
ALPHA10_LENGTH = 450

REPRESENTATIVE_COMPOUNDS = {
    "L_ASCORBATE": {"smiles": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", "class": "parent_reference"},
    "CPD12_POTENT": {"smiles": "C(#C)COC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O", "class": "highly_potent"},
    "CPD25_STRONG": {"smiles": "C(Br)[C@H]([C@H]1OC(=O)C(O)=C1O)O", "class": "potent"},
    "CPD1_MODERATE": {"smiles": "C1([C@@H]([C@H](O)CO)OC(=O)C=1O)O", "class": "moderate"},
    "CPD18_MODERATE": {"smiles": "C1C=CC(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)=CC=1", "class": "moderate"},
    "CPD3_WEAK": {"smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2O)CO1)C", "class": "weak"},
    "CPD8_INACTIVE": {"smiles": "CCCOC1C(C2OC(C)(C)OC2)OC(=O)C=1O", "class": "inactive"},
    "CPD9_INACTIVE": {"smiles": "CCCCOC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O", "class": "inactive"},
}

POCKET_RESIDUES = {
    23: {
        "alpha9": [176, 175, 177, 120, 121, 119, 224, 223, 225],
        "alpha10": [145, 143, 146, 83, 81, 84, 195, 62],
    },
    21: {
        "alpha9": [178, 177, 179, 47, 46, 48, 82],
        "alpha10": [103, 105, 104, 30, 31, 29],
    },
    5: {
        "alpha9": [176, 175, 177, 120, 121],
        "alpha10": [145, 143, 146, 83, 81, 84],
    },
    34: {
        "alpha9": [270, 274, 277, 271, 273],
        "alpha10": [273, 276, 274],
    },
    7: {
        "alpha9": [176, 175, 177, 120, 121],
        "alpha10": [145, 143, 146, 83, 81],
    },
}

REPRESENTATIVE_TO_CPD = {
    "CPD1_MODERATE": 1,
    "CPD12_POTENT": 12,
    "CPD25_STRONG": 25,
    "L_ASCORBATE": 1,
    "CPD18_MODERATE": 18,
    "CPD3_WEAK": 3,
    "CPD8_INACTIVE": 8,
    "CPD9_INACTIVE": 9,
}

SEEDS = [42, 123, 456]

CONTACT_CUTOFF = 4.5
CLASH_CUTOFF = 2.2
BURIED_CUTOFF = 4.5


def load_compounds(csv_path=None):
    path = Path(csv_path or MODULATOR_CSV)
    compounds = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            cid = int(row["Identifier"])
            compounds[cid] = {
                "id": cid,
                "smiles": row["Smiles"],
                "activity_uM": float(row["Activity (uM)"]),
                "potentiation": float(row["%Potentiation"]),
                "active": float(row["Activity (uM)"]) > 0.0,
            }
    return compounds


def activity_class(cid, compounds=None):
    compounds = compounds or load_compounds()
    if compounds[cid]["active"]:
        return "active"
    return "inactive"


ACTIVE_CPDS = [c["id"] for c in load_compounds().values() if c["active"]]
INACTIVE_CPDS = [c["id"] for c in load_compounds().values() if not c["active"]]

NANOMOLAR_PAM = 12
NEAR_ACTIVE_CANDIDATES = [21, 23, 26, 28, 29]

MMP_CLIFF_PAIRS = [
    {"pair": "12 vs 7", "cpd_a": 12, "cpd_b": 7, "change": "alkynyl vs vinyl at C5",
     "effect": "10000x potency (active vs inactive)"},
    {"pair": "25 vs 1", "cpd_a": 25, "cpd_b": 1, "change": "bromo vs OH at C4",
     "effect": "500x potency increase (both active)"},
    {"pair": "18 vs 19", "cpd_a": 18, "cpd_b": 19, "change": "benzyloxymethyl vs allyloxy at C5",
     "effect": "active vs inactive"},
    {"pair": "3 vs 1", "cpd_a": 3, "cpd_b": 1, "change": "gem-dimethyl acetonide vs diol",
     "effect": "3x potency decrease"},
    {"pair": "24 vs 20", "cpd_a": 24, "cpd_b": 20, "change": "3-O-propyl vs 5-O-butyl",
     "effect": "active vs inactive"},
]