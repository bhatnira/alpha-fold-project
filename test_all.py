#!/usr/bin/env python3
"""Lean Pipeline — Automated Test Suite.

Run:  python test_all.py
      python test_all.py --verbose
      python test_all.py --regenerate   (re-runs all phases, then tests)

Exit codes: 0 = all pass, 1 = failures
"""

import sys, os, time, csv, json, traceback, argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PYTHON = sys.executable

# ═══════════════════════════════════════════════════════════════════
# UTILITIES
# ═══════════════════════════════════════════════════════════════════
class Result:
    def __init__(self, name):
        self.name = name
        self.status = "PASS"
        self.error = ""
        self.time = 0.0
        self.details = []

    def fail(self, msg):
        self.status = "FAIL"
        self.error = msg

    def check(self, condition, msg):
        if not condition:
            self.fail(msg)
            return False
        return True

    def csv_rows(self, path):
        p = ROOT / path
        if not p.exists():
            self.fail(f"File not found: {path}")
            return []
        with open(p) as f:
            return list(csv.DictReader(f))

    def json_load(self, path):
        p = ROOT / path
        if not p.exists():
            self.fail(f"File not found: {path}")
            return {}
        with open(p) as f:
            return json.load(f)

    def file_exists(self, path):
        return (ROOT / path).exists()

    def dir_exists(self, path):
        return (ROOT / path).is_dir()

def run(phase_name, func):
    r = Result(phase_name)
    t0 = time.time()
    try:
        func(r)
    except Exception as e:
        r.fail(f"Exception: {e}\n{traceback.format_exc()}")
    r.time = time.time() - t0
    return r

# ═══════════════════════════════════════════════════════════════════
# PHASE TESTS
# ═══════════════════════════════════════════════════════════════════

def test_phase01(r):
    """Phase I — SAR dataset"""
    rows = r.csv_rows("phase01_sar/sar_dataset.csv")
    r.check(len(rows) >= 50, f"SAR records {len(rows)} < 50")

    disc = r.csv_rows("phase01_sar/discovery_set.csv")
    held = r.csv_rows("phase01_sar/held_out_set.csv")
    r.check(len(disc) + len(held) == len(rows), "Split doesn't add up")

    # Validate required columns
    required = {"site_id", "ligand", "enrichment", "n_models"}
    actual = set(rows[0].keys())
    r.check(required.issubset(actual), f"Missing columns: {required - actual}")

    # Check site_id range
    sites = set(int(row["site_id"]) for row in rows)
    r.check(min(sites) >= 0 and max(sites) <= 50, f"Site IDs out of range: {sites}")

    # Check ligand classes present
    ligands = set(row["ligand"] for row in rows)
    r.check("L-ASC" in ligands, "L-ASC missing from SAR")
    r.check("ACETATE" in ligands, "ACETATE missing from SAR")

def test_phase02(r):
    """Phase II — Receptor ensemble"""
    inv = r.csv_rows("phase02_receptor/receptor_inventory.csv")
    r.check(len(inv) >= 10, f"Inventory {len(inv)} < 10 structures")

    sel = r.csv_rows("phase02_receptor/selected_ensemble.csv")
    r.check(len(sel) >= 1, "No receptors selected")

    # Check stoichiometries present
    stoichs = set(row["stoichiometry"] for row in inv)
    r.check("2to3" in stoichs, "No 2to3 stoichiometry")
    r.check("3to2" in stoichs, "No 3to2 stoichiometry")

    # Validate selected ensemble points to real files
    for row in sel:
        path = Path(row["path"])
        r.check(path.exists(), f"Selected receptor missing: {path}")

def test_phase03(r):
    """Phase III — fpocket"""
    fpocket = Path.home() / ".venvs/pymol/bin/fpocket"
    r.check(fpocket.exists(), "fpocket binary not found")

    # Check fpocket outputs exist
    deliverable = ROOT / "data/deliverable/05_structures"
    fpocket_dirs = list(deliverable.rglob("*_out/"))
    r.check(len(fpocket_dirs) >= 2, f"Only {len(fpocket_dirs)} fpocket outputs")

    # Validate at least one info.txt has pocket data
    info_files = list(deliverable.rglob("*info.txt"))
    r.check(len(info_files) >= 1, "No fpocket info.txt files")
    if info_files:
        content = info_files[0].read_text()
        r.check("Pocket" in content, "fpocket info.txt has no pocket data")

def test_phase04(r):
    """Phase IV — Convergence"""
    rows = r.csv_rows("phase04_convergence/convergence_map.csv")
    r.check(len(rows) >= 20, f"Convergence pairs {len(rows)} < 20")

    # Check convergence classes
    classes = set(row["convergence_class"] for row in rows)
    r.check("STRONG" in classes or "MODERATE" in classes, "No strong/moderate convergence")

    # Check contact_jaccard values in [0, 1]
    for row in rows:
        j = float(row["contact_jaccard"])
        r.check(0 <= j <= 1, f"Jaccard out of range: {j}")

    # Check method coverage
    af3_count = sum(1 for row in rows if int(row["n_af3"]) > 0)
    boltz_count = sum(1 for row in rows if int(row["n_boltz2"]) > 0)
    r.check(af3_count >= 5, f"Only {af3_count} AF3 pairs")
    r.check(boltz_count >= 5, f"Only {boltz_count} Boltz2 pairs")

def test_phase05(r):
    """Phase V — Docking"""
    # Check scripts exist
    r.check(r.file_exists("phase05_docking/run_docking.py"), "run_docking.py missing")
    r.check(r.file_exists("phase05_docking/docking_template.txt"), "template missing")

    # Check Vina importable (try system python first, then pymol venv)
    vina_ok = False
    try:
        from vina import Vina
        v = Vina(sf_name='vina')
        vina_ok = True
    except Exception:
        # Try pymol venv
        pymol_python = Path.home() / ".venvs/pymol/bin/python"
        if pymol_python.exists():
            import subprocess
            result = subprocess.run(
                [str(pymol_python), "-c", "from vina import Vina; v = Vina(sf_name='vina'); print('OK')"],
                capture_output=True, text=True, timeout=30
            )
            if result.returncode == 0 and "OK" in result.stdout:
                vina_ok = True
    r.check(vina_ok, "Vina not available in system python or pymol venv")

    # Check docking summary has content
    summary = ROOT / "phase05_docking/docking_summary.txt"
    if summary.exists():
        content = summary.read_text()
        r.check("ENSEMBLE DOCKING" in content, "Summary header missing")

def test_phase06(r):
    """Phase VI — Fingerprints"""
    summary = r.csv_rows("phase06_fingerprints/ifp_summary.csv")
    r.check(len(summary) >= 20, f"IFP sites {len(summary)} < 20")

    # Validate required columns
    required = {"site", "n_contact_residues", "n_core", "core_residues"}
    actual = set(summary[0].keys())
    r.check(required.issubset(actual), f"Missing columns: {required - actual}")

    # Check per-site contact files
    contact_files = list((ROOT / "phase06_fingerprints").glob("site*_contacts.csv"))
    r.check(len(contact_files) >= 15, f"Contact files {len(contact_files)} < 15")

    # Validate a contact file
    if contact_files:
        with open(contact_files[0]) as f:
            contacts = list(csv.DictReader(f))
        r.check(len(contacts) >= 1, f"Contact file empty: {contact_files[0].name}")
        r.check("contact_frequency" in contacts[0], "Missing contact_frequency column")

def test_phase07(r):
    """Phase VII — SAR model"""
    rules = r.csv_rows("phase07_sar_model/sar_rules.csv")
    r.check(len(rules) >= 50, f"SAR rules {len(rules)} < 50")

    # Check rule types
    rule_types = set(row.get("rule_type", "") for row in rules)
    r.check("REQUIRED" in rule_types, "No REQUIRED rules")
    r.check("STEREOSENSITIVE" in rule_types, "No STEREOSENSITIVE rules")

    # Validate residues referenced
    residues = set(row.get("residue", "") for row in rules)
    r.check(len(residues) >= 10, f"Only {len(residues)} unique residues in rules")

def test_phase08(r):
    """Phase VIII — Boltz-2 affinity"""
    matrix = r.csv_rows("phase08_boltz_affinity/affinity_matrix.csv")
    r.check(len(matrix) >= 30, f"Affinity records {len(matrix)} < 30")

    # Check affinity scores are numeric
    for row in matrix[:5]:
        score = float(row["affinity_score"])
        r.check(score >= 0, f"Negative affinity score: {score}")

    # Check comparisons exist
    comps = r.csv_rows("phase08_boltz_affinity/affinity_comparisons.csv")
    r.check(len(comps) >= 5, f"Affinity comparisons {len(comps)} < 5")

def test_phase09(r):
    """Phase IX — Specificity"""
    results = r.csv_rows("phase09_specificity/specificity_results.csv")
    r.check(len(results) >= 15, f"Specificity records {len(results)} < 15")

    # Check discrimination classes
    charge_classes = set(row["charge_class"] for row in results)
    r.check("MOLECULAR_RECOGNITION" in charge_classes or "CHARGE_DRIVEN" in charge_classes,
            "No charge discrimination classes")

    stereo_classes = set(row["stereo_class"] for row in results)
    r.check("STEREOSPECIFIC" in stereo_classes or "NON_STEREOSPECIFIC" in stereo_classes,
            "No stereo discrimination classes")

    # Check Jaccard values in range
    for row in results:
        j = float(row["charge_jaccard"])
        r.check(0 <= j <= 1, f"Charge Jaccard out of range: {j}")

def test_phase10(r):
    """Phase X — Electrostatics"""
    pots = r.csv_rows("phase10_electrostatics/electrostatic_potentials.csv")
    r.check(len(pots) >= 10, f"Potential records {len(pots)} < 10")

    # Check PQR files exist
    pqr_files = list((ROOT / "phase10_electrostatics/results").glob("*.pqr"))
    r.check(len(pqr_files) >= 1, "No PQR files")

    # Validate potential values are numeric
    for row in pots[:3]:
        pot = float(row["potential_kcal_mol"])
        r.check(-1000 < pot < 1000, f"Potential out of range: {pot}")

def test_phase11(r):
    """Phase XI — Robustness"""
    matrix = r.csv_rows("phase11_robustness/robustness_matrix.csv")
    r.check(len(matrix) >= 20, f"Robustness records {len(matrix)} < 20")

    # Check robustness scores in [0, 1]
    for row in matrix[:5]:
        score = float(row["overall_robustness"])
        r.check(0 <= score <= 1, f"Robustness out of range: {score}")

    # Check component scores exist
    for key in ["seed_score", "method_score", "stoich_score", "residue_score"]:
        r.check(key in matrix[0], f"Missing component: {key}")

def test_phase12(r):
    """Phase XII — Subtype comparison"""
    comp = r.csv_rows("phase12_subtype/subtype_comparison.csv")
    r.check(len(comp) >= 2, f"Subtypes {len(comp)} < 2")

    # Check α7 and α4β2 present
    subtypes = set(row["subtype"] for row in comp)
    r.check("α7" in subtypes, "α7 missing")
    r.check("α4β2" in subtypes, "α4β2 missing")

    # Check conservation scores
    for row in comp:
        conserv = float(row["a9_pharma_conservation"])
        r.check(0 <= conserv <= 1, f"Conservation out of range: {conserv}")

def test_phase13(r):
    """Phase XIII — Ryanodine"""
    analysis = r.csv_rows("phase13_ryanodine/ryanodine_analysis.csv")
    r.check(len(analysis) >= 1, "No ryanodine analysis")

    # Check key metrics present
    row = analysis[0]
    r.check("ryan_n_models" in row, "Missing ryan_n_models")
    r.check("contact_jaccard" in row, "Missing contact_jaccard")
    r.check("same_site" in row, "Missing same_site")

    # Ryanodine should NOT be at the same site as L-ascorbate
    jaccard = float(row["contact_jaccard"])
    r.check(jaccard < 0.3, f"Ryanodine Jaccard {jaccard} >= 0.3 (should be distinct)")

def test_phase14(r):
    """Phase XIV — Evidence"""
    matrix = r.csv_rows("phase14_evidence/evidence_matrix.csv")
    r.check(len(matrix) >= 30, f"Evidence records {len(matrix)} < 30")

    # Check classification distribution
    classes = set(row["classification"] for row in matrix)
    r.check("HIGH_CONFIDENCE" in classes or "INTERMEDIATE" in classes,
            "No high/intermediate confidence sites")

    # Check evidence scores in [0, 1]
    for row in matrix[:5]:
        score = float(row["evidence_score"])
        r.check(0 <= score <= 1, f"Evidence score out of range: {score}")

    # Check all evidence components present
    ev_keys = [k for k in matrix[0].keys() if k.startswith("evidence_")]
    r.check(len(ev_keys) >= 10, f"Only {len(ev_keys)} evidence components")

def test_phase15(r):
    """Phase XV — Pharmacophore"""
    features = r.csv_rows("phase15_pharmacophore/pharmacophore_features.csv")
    r.check(len(features) >= 3, f"Features {len(features)} < 3")

    # Check feature types
    types = set(row["type"] for row in features)
    r.check("H-bond donor" in types, "No HBD feature")
    r.check("H-bond acceptor" in types, "No HBA feature")

    # Check PML file
    pml = ROOT / "phase15_pharmacophore/pharmacophore.pml"
    r.check(pml.exists(), "PML file missing")
    content = pml.read_text()
    r.check("load" in content, "PML has no load commands")
    r.check("set_color" in content, "PML has no colors")
    r.check("spectrum" in content.lower() or "color" in content.lower(), "PML has no visualization")

def test_phase16(r):
    """Phase XVI — Molecule generation"""
    lib = r.csv_rows("phase16_molgen/molecule_library.csv")
    r.check(len(lib) >= 5, f"Molecules {len(lib)} < 5")

    # Validate SMILES
    try:
        from rdkit import Chem
        valid = sum(1 for row in lib if Chem.MolFromSmiles(row.get("smiles", "")) is not None)
        r.check(valid >= 3, f"Only {valid} valid SMILES")
    except ImportError:
        r.check(True, "RDKit not available, skipping SMILES validation")

    # Check required columns
    required = {"molecule_id", "smiles", "source"}
    actual = set(lib[0].keys())
    r.check(required.issubset(actual), f"Missing columns: {required - actual}")

def test_phase17(r):
    """Phase XVII — Cascading filter"""
    filtered = r.csv_rows("phase17_cascade/filtered_candidates.csv")
    r.check(len(filtered) >= 1, "No filtered candidates")

    library = r.csv_rows("phase16_molgen/molecule_library.csv")
    r.check(len(filtered) <= len(library), "Filter produced more than input")

    # Check cascade summary exists
    summary = ROOT / "phase17_cascade/cascade_summary.txt"
    r.check(summary.exists(), "Cascade summary missing")
    content = summary.read_text()
    r.check("Funnel stages" in content, "Summary incomplete")

def test_phase18(r):
    """Phase XVIII — Boltz-2 design"""
    subs = r.csv_rows("phase18_boltz_design/boltz2_submissions.csv")
    r.check(len(subs) >= 1, "No Boltz-2 submissions")

    ranking = r.csv_rows("phase18_boltz_design/design_ranking.csv")
    r.check(len(ranking) >= 1, "No design ranking")

    # Check submission JSON files exist
    sub_dir = ROOT / "phase18_boltz_design/boltz2_submissions"
    r.check(sub_dir.exists(), "Submissions directory missing")
    jsons = list(sub_dir.glob("*.json"))
    r.check(len(jsons) >= 1, "No submission JSON files")

    # Validate JSON structure
    if jsons:
        with open(jsons[0]) as f:
            data = json.load(f)
        r.check("name" in data, "Submission missing 'name'")
        r.check("sequences" in data, "Submission missing 'sequences'")

def test_phase19(r):
    """Phase XIX — Selectivity"""
    scores = r.csv_rows("phase19_selectivity/selectivity_scores.csv")
    r.check(len(scores) >= 20, f"Selectivity scores {len(scores)} < 20")

    # Check scores in [0, 1]
    for row in scores[:5]:
        score = float(row["selectivity_score"])
        r.check(0 <= score <= 1, f"Selectivity out of range: {score}")

    # Check all sites have scores
    site_ids = set(row["site_id"] for row in scores)
    r.check(len(site_ids) >= 20, f"Only {len(site_ids)} sites have selectivity scores")

def test_phase20(r):
    """Phase XX — Ranking"""
    ranking = r.csv_rows("phase20_ranking/site_ranking.csv")
    r.check(len(ranking) >= 20, f"Ranked sites {len(ranking)} < 20")

    # Check ranking is sorted by weighted_score (descending)
    scores = [float(row["weighted_score"]) for row in ranking]
    r.check(scores == sorted(scores, reverse=True), "Ranking not sorted")

    # Check rank column is sequential
    ranks = [int(row["rank"]) for row in ranking]
    r.check(ranks == list(range(1, len(ranks) + 1)), "Ranks not sequential")

    # Check objective columns exist
    obj_keys = [k for k in ranking[0].keys() if k.startswith("obj_")]
    r.check(len(obj_keys) >= 5, f"Only {len(obj_keys)} objective columns")

def test_phase21(r):
    """Phase XXI — Validation"""
    plan = r.json_load("phase21_validation/validation_plan.json")
    r.check("top_sites" in plan, "No top_sites")
    r.check(len(plan["top_sites"]) >= 1, "No top sites")
    r.check("mutation_plan" in plan, "No mutation_plan")
    r.check("assay_plan" in plan, "No assay_plan")

    # Check modulator validation section
    r.check("modulator_validation" in plan, "No modulator_validation")
    mv = plan["modulator_validation"]
    r.check("spearman_rho_pic50_enrichment" in mv, "No Spearman rho")
    rho = float(mv["spearman_rho_pic50_enrichment"])
    r.check(rho > 0.5, f"Spearman rho {rho} < 0.5 (expected STRONG)")

    # Check mutation plan has required fields
    if plan["mutation_plan"]:
        mut = plan["mutation_plan"][0]
        r.check("mutation" in mut, "Mutation missing 'mutation'")
        r.check("hypothesis" in mut, "Mutation missing 'hypothesis'")

    # Check assay plan
    r.check(len(plan["assay_plan"]) >= 2, "Fewer than 2 assays planned")

    # Check modulator compounds to test
    r.check("modulator_compounds_to_test" in plan, "No modulator_compounds_to_test")
    r.check(len(plan["modulator_compounds_to_test"]) >= 1, "No compounds to test")

def test_phase22(r):
    """Phase XXII — Closed-loop"""
    template = r.json_load("phase22_closed_loop/closed_loop_template.json")
    r.check("steps" in template, "No steps")
    r.check(len(template["steps"]) >= 5, f"Only {len(template['steps'])} steps")

    # Check each step has required fields
    for step in template["steps"]:
        r.check("name" in step, f"Step missing 'name'")
        r.check("action" in step, f"Step missing 'action'")

    # Check success criteria
    r.check("success_criteria" in template, "No success_criteria")
    r.check("stop_criteria" in template, "No stop_criteria")

# ═══════════════════════════════════════════════════════════════════
# CROSS-PHASE INTEGRITY TESTS
# ═══════════════════════════════════════════════════════════════════

def test_cross_phase_data_flow(r):
    """Verify data flows correctly between dependent phases."""
    # Phase I → VI: SAR sites appear in fingerprints
    sar = r.csv_rows("phase01_sar/sar_dataset.csv")
    sar_sites = set(int(row["site_id"]) for row in sar)
    ifp = r.csv_rows("phase06_fingerprints/ifp_summary.csv")
    ifp_sites = set(int(row["site"]) for row in ifp)
    overlap = sar_sites & ifp_sites
    r.check(len(overlap) >= 10, f"SAR-IFP overlap only {len(overlap)}")

    # Phase VI → XIV: IFP feeds evidence
    ev = r.csv_rows("phase14_evidence/evidence_matrix.csv")
    ev_sites = set(int(row["site_id"]) for row in ev)
    r.check(len(ifp_sites & ev_sites) >= 10, "IFP-Evidence gap")

    # Phase XIV → XX: Evidence feeds ranking
    rank = r.csv_rows("phase20_ranking/site_ranking.csv")
    rank_sites = set(int(row["site_id"]) for row in rank)
    r.check(len(ev_sites & rank_sites) >= 10, "Evidence-Ranking gap")

    # Phase XX → XXI: Ranking feeds validation
    plan = r.json_load("phase21_validation/validation_plan.json")
    top_site_ids = set(int(s["site_id"]) for s in plan["top_sites"])
    r.check(len(top_site_ids & rank_sites) >= 1, "Top sites not in ranking")

    # Phase XV → XVI → XVII: Pharmacophore → molecules → filter
    lib = r.csv_rows("phase16_molgen/molecule_library.csv")
    filt = r.csv_rows("phase17_cascade/filtered_candidates.csv")
    r.check(len(filt) <= len(lib), "Filter inflated library")

    # Phase XVI → XVIII: Molecules → Boltz-2 submissions
    subs = r.csv_rows("phase18_boltz_design/boltz2_submissions.csv")
    r.check(len(subs) >= 1, "No submissions from filtered candidates")

def test_site_consistency(r):
    """Verify site IDs are consistent across all phases."""
    site_sets = {}
    for phase, path, col in [
        ("I",   "phase01_sar/sar_dataset.csv", "site_id"),
        ("IV",  "phase04_convergence/convergence_map.csv", "site"),
        ("VI",  "phase06_fingerprints/ifp_summary.csv", "site"),
        ("IX",  "phase09_specificity/specificity_results.csv", "site"),
        ("XI",  "phase11_robustness/robustness_matrix.csv", "site"),
        ("XIV", "phase14_evidence/evidence_matrix.csv", "site_id"),
        ("XIX", "phase19_selectivity/selectivity_scores.csv", "site_id"),
        ("XX",  "phase20_ranking/site_ranking.csv", "site_id"),
    ]:
        rows = r.csv_rows(path)
        if rows:
            site_sets[phase] = set(int(row[col]) for row in rows)

    # All phases should share a core set of sites
    if len(site_sets) >= 3:
        core = set.intersection(*site_sets.values())
        r.check(len(core) >= 5, f"Core site set only {len(core)} sites")

def test_output_completeness(r):
    """Verify every phase has at least one output file."""
    phase_dirs = sorted(ROOT.glob("phase*"))
    for d in phase_dirs:
        if d.is_dir():
            files = list(d.glob("*"))
            csvs = [f for f in files if f.suffix == ".csv"]
            txts = [f for f in files if f.suffix == ".txt"]
            others = [f for f in files if f.suffix in (".json", ".pml")]
            dirs = [f for f in files if f.is_dir()]
            total = len(csvs) + len(txts) + len(others) + len(dirs)
            r.check(total >= 1, f"{d.name} has no outputs")

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--regenerate", action="store_true",
                        help="Re-run all phase scripts before testing")
    args = parser.parse_args()

    print("=" * 70)
    print("LEAN PIPELINE — AUTOMATED TEST SUITE")
    print("=" * 70)
    print()

    # Optionally regenerate all phases
    if args.regenerate:
        print("REGENERATING ALL PHASES...")
        import subprocess
        scripts = sorted(ROOT.glob("phase*/build_*.py")) + \
                  sorted(ROOT.glob("phase*/compute_*.py")) + \
                  sorted(ROOT.glob("phase*/derive_*.py")) + \
                  sorted(ROOT.glob("phase*/map_*.py")) + \
                  sorted(ROOT.glob("phase*/analyze_*.py")) + \
                  sorted(ROOT.glob("phase*/compare_*.py")) + \
                  sorted(ROOT.glob("phase*/design_*.py")) + \
                  sorted(ROOT.glob("phase*/filter_*.py")) + \
                  sorted(ROOT.glob("phase*/run_*.py")) + \
                  sorted(ROOT.glob("phase*/create_*.py")) + \
                  sorted(ROOT.glob("phase*/generate_*.py"))
        for script in scripts:
            if "test_" in script.name:
                continue
            print(f"  Running {script.name}...")
            subprocess.run([PYTHON, str(script)], capture_output=True, timeout=120)
        print("  Regeneration complete.\n")

    tests = [
        ("Phase I — SAR dataset",            test_phase01),
        ("Phase II — Receptor ensemble",      test_phase02),
        ("Phase III — fpocket",               test_phase03),
        ("Phase IV — Convergence",            test_phase04),
        ("Phase V — Docking",                 test_phase05),
        ("Phase VI — Fingerprints",           test_phase06),
        ("Phase VII — SAR model",             test_phase07),
        ("Phase VIII — Boltz-2 affinity",     test_phase08),
        ("Phase IX — Specificity",            test_phase09),
        ("Phase X — Electrostatics",          test_phase10),
        ("Phase XI — Robustness",             test_phase11),
        ("Phase XII — Subtype comparison",    test_phase12),
        ("Phase XIII — Ryanodine",            test_phase13),
        ("Phase XIV — Evidence",              test_phase14),
        ("Phase XV — Pharmacophore",          test_phase15),
        ("Phase XVI — Molecule generation",   test_phase16),
        ("Phase XVII — Cascading filter",     test_phase17),
        ("Phase XVIII — Boltz-2 design",      test_phase18),
        ("Phase XIX — Selectivity",           test_phase19),
        ("Phase XX — Ranking",                test_phase20),
        ("Phase XXI — Validation",            test_phase21),
        ("Phase XXII — Closed-loop",          test_phase22),
        ("Cross-phase data flow",             test_cross_phase_data_flow),
        ("Site ID consistency",               test_site_consistency),
        ("Output completeness",               test_output_completeness),
    ]

    results = []
    for name, func in tests:
        icon = "⏳"
        print(f"  {icon} {name:40s}", end="", flush=True)
        result = run(name, func)
        results.append(result)
        icon = "✅" if result.status == "PASS" else "❌"
        print(f"\r  {icon} {name:40s} {result.time:>5.1f}s  {result.error[:60] if result.error else ''}")

    # Summary
    print()
    print("=" * 70)
    n_pass = sum(1 for r in results if r.status == "PASS")
    n_fail = sum(1 for r in results if r.status == "FAIL")
    total = len(results)

    if n_fail == 0:
        print(f"✅ ALL {total} TESTS PASSED")
    else:
        print(f"❌ {n_fail}/{total} TESTS FAILED")
        print()
        for r in results:
            if r.status == "FAIL":
                print(f"  FAILED: {r.name}")
                print(f"          {r.error[:200]}")
                print()

    # Stats
    total_time = sum(r.time for r in results)
    print(f"\nTotal time: {total_time:.1f}s")

    # Save results
    report = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "total": total,
        "passed": n_pass,
        "failed": n_fail,
        "total_time_s": round(total_time, 1),
        "tests": [{
            "name": r.name,
            "status": r.status,
            "time_s": round(r.time, 2),
            "error": r.error[:200] if r.error else "",
        } for r in results],
    }
    report_path = ROOT / "test_results.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Report saved to: {report_path}")

    return n_fail == 0

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
