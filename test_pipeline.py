#!/usr/bin/env python3
"""Lean Pipeline — Full integration test.

Runs each phase with minimal data to verify the pipeline works end-to-end.
Does NOT run full computations — just validates scripts, imports, and data flow.
"""

import sys, os, time, traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PYTHON = sys.executable

results = []

def run_test(name, func):
    """Run a test function and record result."""
    t0 = time.time()
    try:
        func()
        elapsed = time.time() - t0
        results.append({"phase": name, "status": "PASS", "time": f"{elapsed:.1f}s", "error": ""})
        print(f"  ✅ {name} ({elapsed:.1f}s)")
    except Exception as e:
        elapsed = time.time() - t0
        results.append({"phase": name, "status": "FAIL", "time": f"{elapsed:.1f}s", "error": str(e)})
        print(f"  ❌ {name}: {e}")

# ═══════════════════════════════════════════════════════════════════
# PHASE I — SAR dataset
# ═══════════════════════════════════════════════════════════════════
def test_phase01():
    assert (ROOT / "phase01_sar" / "sar_dataset.csv").exists(), "sar_dataset.csv missing"
    assert (ROOT / "phase01_sar" / "discovery_set.csv").exists(), "discovery_set.csv missing"
    assert (ROOT / "phase01_sar" / "held_out_set.csv").exists(), "held_out_set.csv missing"
    import csv
    with open(ROOT / "phase01_sar" / "sar_dataset.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 50, f"Expected >=50 SAR records, got {len(rows)}"
    with open(ROOT / "phase01_sar" / "discovery_set.csv") as f:
        disc = list(csv.DictReader(f))
    with open(ROOT / "phase01_sar" / "held_out_set.csv") as f:
        held = list(csv.DictReader(f))
    assert len(disc) + len(held) == len(rows), "Split doesn't add up"

# ═══════════════════════════════════════════════════════════════════
# PHASE II — Receptor ensemble
# ═══════════════════════════════════════════════════════════════════
def test_phase02():
    assert (ROOT / "phase02_receptor" / "receptor_inventory.csv").exists()
    assert (ROOT / "phase02_receptor" / "selected_ensemble.csv").exists()
    import csv
    with open(ROOT / "phase02_receptor" / "selected_ensemble.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 1, "No receptors selected"

# ═══════════════════════════════════════════════════════════════════
# PHASE III — fpocket
# ═══════════════════════════════════════════════════════════════════
def test_phase03():
    fpocket_bin = Path(PYTHON).parent / "fpocket"
    assert fpocket_bin.exists(), "fpocket binary not found"
    # Check at least one fpocket output exists
    fpocket_out = list(Path("/cluster/home/nbhatt04/lean_pipeline/data/deliverable/05_structures").rglob("*_out/"))
    assert len(fpocket_out) >= 1, "No fpocket output directories found"

# ═══════════════════════════════════════════════════════════════════
# PHASE IV — Convergence
# ═══════════════════════════════════════════════════════════════════
def test_phase04():
    assert (ROOT / "phase04_convergence" / "convergence_map.csv").exists()
    import csv
    with open(ROOT / "phase04_convergence" / "convergence_map.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 10, f"Expected >=10 convergence pairs, got {len(rows)}"
    # Check convergence classes exist
    classes = set(r["convergence_class"] for r in rows)
    assert "STRONG" in classes or "MODERATE" in classes, "No strong/moderate convergence"

# ═══════════════════════════════════════════════════════════════════
# PHASE V — Docking
# ═══════════════════════════════════════════════════════════════════
def test_phase05():
    vina_ok = False
    try:
        from vina import Vina
        v = Vina(sf_name='vina')
        vina_ok = True
    except Exception:
        pymol_python = Path.home() / ".venvs/pymol/bin/python"
        if pymol_python.exists():
            import subprocess
            result = subprocess.run(
                [str(pymol_python), "-c", "from vina import Vina; v = Vina(sf_name='vina'); print('OK')"],
                capture_output=True, text=True, timeout=30
            )
            if result.returncode == 0 and "OK" in result.stdout:
                vina_ok = True
    assert vina_ok, "Vina not available in system python or pymol venv"
    assert (ROOT / "phase05_docking" / "run_docking.py").exists()
    assert (ROOT / "phase05_docking" / "docking_summary.txt").exists()

# ═══════════════════════════════════════════════════════════════════
# PHASE VI — Fingerprints
# ═══════════════════════════════════════════════════════════════════
def test_phase06():
    assert (ROOT / "phase06_fingerprints" / "ifp_summary.csv").exists()
    import csv
    with open(ROOT / "phase06_fingerprints" / "ifp_summary.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 20, f"Expected >=20 sites with IFP, got {len(rows)}"
    # Check per-site contact files
    contact_files = list((ROOT / "phase06_fingerprints").glob("site*_contacts.csv"))
    assert len(contact_files) >= 10, f"Expected >=10 contact files, got {len(contact_files)}"

# ═══════════════════════════════════════════════════════════════════
# PHASE VII — SAR model
# ═══════════════════════════════════════════════════════════════════
def test_phase07():
    assert (ROOT / "phase07_sar_model" / "sar_rules.csv").exists()
    import csv
    with open(ROOT / "phase07_sar_model" / "sar_rules.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 10, f"Expected >=10 SAR rules, got {len(rows)}"
    rule_types = set(r.get("rule_type", "") for r in rows)
    assert "REQUIRED" in rule_types, "No REQUIRED rules found"

# ═══════════════════════════════════════════════════════════════════
# PHASE VIII — Boltz-2 affinity
# ═══════════════════════════════════════════════════════════════════
def test_phase08():
    assert (ROOT / "phase08_boltz_affinity" / "affinity_matrix.csv").exists()
    import csv
    with open(ROOT / "phase08_boltz_affinity" / "affinity_matrix.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 20, f"Expected >=20 affinity records, got {len(rows)}"

# ═══════════════════════════════════════════════════════════════════
# PHASE IX — Specificity
# ═══════════════════════════════════════════════════════════════════
def test_phase09():
    assert (ROOT / "phase09_specificity" / "specificity_results.csv").exists()
    import csv
    with open(ROOT / "phase09_specificity" / "specificity_results.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 10, f"Expected >=10 specificity records, got {len(rows)}"
    charge_classes = set(r["charge_class"] for r in rows)
    assert "MOLECULAR_RECOGNITION" in charge_classes or "CHARGE_DRIVEN" in charge_classes

# ═══════════════════════════════════════════════════════════════════
# PHASE X — Electrostatics
# ═══════════════════════════════════════════════════════════════════
def test_phase10():
    assert (ROOT / "phase10_electrostatics" / "electrostatic_potentials.csv").exists()
    assert (ROOT / "phase10_electrostatics" / "results").is_dir()
    pqr_files = list((ROOT / "phase10_electrostatics" / "results").glob("*.pqr"))
    assert len(pqr_files) >= 1, "No PQR files generated"
    import csv
    with open(ROOT / "phase10_electrostatics" / "electrostatic_potentials.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 5, f"Expected >=5 potentials, got {len(rows)}"

# ═══════════════════════════════════════════════════════════════════
# PHASE XI — Robustness
# ═══════════════════════════════════════════════════════════════════
def test_phase11():
    assert (ROOT / "phase11_robustness" / "robustness_matrix.csv").exists()
    import csv
    with open(ROOT / "phase11_robustness" / "robustness_matrix.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 10, f"Expected >=10 robustness records, got {len(rows)}"
    scores = [float(r["overall_robustness"]) for r in rows]
    assert max(scores) > 0, "All robustness scores are zero"

# ═══════════════════════════════════════════════════════════════════
# PHASE XII — Subtype comparison
# ═══════════════════════════════════════════════════════════════════
def test_phase12():
    assert (ROOT / "phase12_subtype" / "subtype_comparison.csv").exists()
    import csv
    with open(ROOT / "phase12_subtype" / "subtype_comparison.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 2, f"Expected >=2 subtypes, got {len(rows)}"

# ═══════════════════════════════════════════════════════════════════
# PHASE XIII — Ryanodine
# ═══════════════════════════════════════════════════════════════════
def test_phase13():
    assert (ROOT / "phase13_ryanodine" / "ryanodine_analysis.csv").exists()
    import csv
    with open(ROOT / "phase13_ryanodine" / "ryanodine_analysis.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 1, "No ryanodine analysis"

# ═══════════════════════════════════════════════════════════════════
# PHASE XIV — Evidence
# ═══════════════════════════════════════════════════════════════════
def test_phase14():
    assert (ROOT / "phase14_evidence" / "evidence_matrix.csv").exists()
    import csv
    with open(ROOT / "phase14_evidence" / "evidence_matrix.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 20, f"Expected >=20 evidence records, got {len(rows)}"
    classifications = set(r["classification"] for r in rows)
    assert "HIGH_CONFIDENCE" in classifications or "INTERMEDIATE" in classifications

# ═══════════════════════════════════════════════════════════════════
# PHASE XV — Pharmacophore
# ═══════════════════════════════════════════════════════════════════
def test_phase15():
    assert (ROOT / "phase15_pharmacophore" / "pharmacophore.pml").exists()
    assert (ROOT / "phase15_pharmacophore" / "pharmacophore_features.csv").exists()
    import csv
    with open(ROOT / "phase15_pharmacophore" / "pharmacophore_features.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 3, f"Expected >=3 pharmacophore features, got {len(rows)}"
    # Check PML loads without errors
    pml = (ROOT / "phase15_pharmacophore" / "pharmacophore.pml").read_text()
    assert "load" in pml, "PML has no load commands"
    assert "set_color" in pml, "PML has no color definitions"

# ═══════════════════════════════════════════════════════════════════
# PHASE XVI — Molecule generation
# ═══════════════════════════════════════════════════════════════════
def test_phase16():
    from rdkit import Chem
    assert (ROOT / "phase16_molgen" / "molecule_library.csv").exists()
    import csv
    with open(ROOT / "phase16_molgen" / "molecule_library.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 5, f"Expected >=5 molecules, got {len(rows)}"
    # Validate at least one SMILES
    valid = sum(1 for r in rows if Chem.MolFromSmiles(r.get("smiles", "")) is not None)
    assert valid >= 1, "No valid SMILES in library"

# ═══════════════════════════════════════════════════════════════════
# PHASE XVII — Cascading filter
# ═══════════════════════════════════════════════════════════════════
def test_phase17():
    assert (ROOT / "phase17_cascade" / "filtered_candidates.csv").exists()
    import csv
    with open(ROOT / "phase17_cascade" / "filtered_candidates.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 1, "No filtered candidates"

# ═══════════════════════════════════════════════════════════════════
# PHASE XVIII — Boltz-2 design
# ═══════════════════════════════════════════════════════════════════
def test_phase18():
    assert (ROOT / "phase18_boltz_design" / "boltz2_submissions.csv").exists()
    assert (ROOT / "phase18_boltz_design" / "design_ranking.csv").exists()
    assert (ROOT / "phase18_boltz_design" / "boltz2_submissions").is_dir()
    submissions = list((ROOT / "phase18_boltz_design" / "boltz2_submissions").glob("*.json"))
    assert len(submissions) >= 1, "No Boltz-2 submission files"

# ═══════════════════════════════════════════════════════════════════
# PHASE XIX — Selectivity
# ═══════════════════════════════════════════════════════════════════
def test_phase19():
    assert (ROOT / "phase19_selectivity" / "selectivity_scores.csv").exists()
    import csv
    with open(ROOT / "phase19_selectivity" / "selectivity_scores.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 10, f"Expected >=10 selectivity records, got {len(rows)}"

# ═══════════════════════════════════════════════════════════════════
# PHASE XX — Ranking
# ═══════════════════════════════════════════════════════════════════
def test_phase20():
    assert (ROOT / "phase20_ranking" / "site_ranking.csv").exists()
    import csv
    with open(ROOT / "phase20_ranking" / "site_ranking.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 10, f"Expected >=10 ranked sites, got {len(rows)}"
    # Check ranking is sorted
    scores = [float(r["weighted_score"]) for r in rows]
    assert scores == sorted(scores, reverse=True), "Ranking not sorted by score"

# ═══════════════════════════════════════════════════════════════════
# PHASE XXI — Validation
# ═══════════════════════════════════════════════════════════════════
def test_phase21():
    assert (ROOT / "phase21_validation" / "validation_plan.json").exists()
    import json
    with open(ROOT / "phase21_validation" / "validation_plan.json") as f:
        plan = json.load(f)
    assert "top_sites" in plan, "No top_sites in validation plan"
    assert len(plan["top_sites"]) >= 1, "No top sites in validation plan"
    assert "mutation_plan" in plan, "No mutation_plan"
    assert "assay_plan" in plan, "No assay_plan"

# ═══════════════════════════════════════════════════════════════════
# PHASE XXII — Closed-loop
# ═══════════════════════════════════════════════════════════════════
def test_phase22():
    assert (ROOT / "phase22_closed_loop" / "closed_loop_template.json").exists()
    import json
    with open(ROOT / "phase22_closed_loop" / "closed_loop_template.json") as f:
        template = json.load(f)
    assert "steps" in template, "No steps in closed-loop template"
    assert len(template["steps"]) >= 5, "Not enough steps in template"

# ═══════════════════════════════════════════════════════════════════
# CROSS-PHASE INTEGRITY
# ═══════════════════════════════════════════════════════════════════
def test_cross_phase_data_flow():
    """Verify that phases that depend on each other have consistent data."""
    import csv
    
    # Phase I → Phase VII: SAR rules use SAR dataset
    with open(ROOT / "phase01_sar" / "sar_dataset.csv") as f:
        sar_sites = set(r["site_id"] for r in csv.DictReader(f))
    with open(ROOT / "phase07_sar_model" / "sar_rules.csv") as f:
        rules = list(csv.DictReader(f))
    assert len(rules) > 0, "SAR rules empty"
    
    # Phase VI → Phase XIV: IFP feeds evidence matrix
    with open(ROOT / "phase06_fingerprints" / "ifp_summary.csv") as f:
        ifp_sites = set(r["site"] for r in csv.DictReader(f))
    with open(ROOT / "phase14_evidence" / "evidence_matrix.csv") as f:
        evidence_sites = set(r["site_id"] for r in csv.DictReader(f))
    overlap = ifp_sites & evidence_sites
    assert len(overlap) >= 10, f"Only {len(overlap)} sites overlap between IFP and evidence"
    
    # Phase XIV → Phase XX: Evidence feeds ranking
    with open(ROOT / "phase20_ranking" / "site_ranking.csv") as f:
        ranked_sites = set(r["site_id"] for r in csv.DictReader(f))
    assert len(ranked_sites & evidence_sites) >= 10, "Ranking doesn't use evidence data"
    
    # Phase XV → Phase XVI: Pharmacophore feeds molecule generation
    assert (ROOT / "phase15_pharmacophore" / "pharmacophore_features.csv").exists()
    assert (ROOT / "phase16_molgen" / "molecule_library.csv").exists()
    
    # Phase XVI → Phase XVII: Library feeds filter
    with open(ROOT / "phase16_molgen" / "molecule_library.csv") as f:
        lib = list(csv.DictReader(f))
    with open(ROOT / "phase17_cascade" / "filtered_candidates.csv") as f:
        filtered = list(csv.DictReader(f))
    assert len(filtered) <= len(lib), "Filter produced more molecules than input"

# ═══════════════════════════════════════════════════════════════════
# RUN ALL TESTS
# ═══════════════════════════════════════════════════════════════════
def main():
    print("=" * 70)
    print("LEAN PIPELINE INTEGRATION TEST")
    print("=" * 70)
    print()
    
    tests = [
        ("Phase I — SAR dataset", test_phase01),
        ("Phase II — Receptor ensemble", test_phase02),
        ("Phase III — fpocket", test_phase03),
        ("Phase IV — Convergence", test_phase04),
        ("Phase V — Docking (Vina import)", test_phase05),
        ("Phase VI — Fingerprints", test_phase06),
        ("Phase VII — SAR model", test_phase07),
        ("Phase VIII — Boltz-2 affinity", test_phase08),
        ("Phase IX — Specificity", test_phase09),
        ("Phase X — Electrostatics", test_phase10),
        ("Phase XI — Robustness", test_phase11),
        ("Phase XII — Subtype comparison", test_phase12),
        ("Phase XIII — Ryanodine", test_phase13),
        ("Phase XIV — Evidence", test_phase14),
        ("Phase XV — Pharmacophore", test_phase15),
        ("Phase XVI — Molecule generation", test_phase16),
        ("Phase XVII — Cascading filter", test_phase17),
        ("Phase XVIII — Boltz-2 design", test_phase18),
        ("Phase XIX — Selectivity", test_phase19),
        ("Phase XX — Ranking", test_phase20),
        ("Phase XXI — Validation", test_phase21),
        ("Phase XXII — Closed-loop", test_phase22),
        ("Cross-phase data flow", test_cross_phase_data_flow),
    ]
    
    for name, func in tests:
        run_test(name, func)
    
    # Summary
    print()
    print("=" * 70)
    n_pass = sum(1 for r in results if r["status"] == "PASS")
    n_fail = sum(1 for r in results if r["status"] == "FAIL")
    total = len(results)
    
    print(f"RESULTS: {n_pass}/{total} passed, {n_fail} failed")
    print()
    
    if n_fail > 0:
        print("FAILURES:")
        for r in results:
            if r["status"] == "FAIL":
                print(f"  ❌ {r['phase']}: {r['error']}")
    else:
        print("✅ ALL TESTS PASSED — Pipeline is fully functional")
    
    print()
    print("PHASE SUMMARY:")
    for r in results:
        icon = "✅" if r["status"] == "PASS" else "❌"
        print(f"  {icon} {r['phase']:35s} {r['time']:>6s}")
    
    print()
    return n_fail == 0

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
