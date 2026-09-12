# α9α10 nAChR Allosteric-Site Discovery — Lean Pipeline

No MD. No expensive simulations. Maximum evidence per unit of compute.

---

## Principles

1. Pocket discovery is ligand-independent.
2. No single AF3 prediction carries the argument.
3. No single Boltz-2 prediction carries the argument.
4. No single docking pose carries the argument.
5. No docking score alone explains SAR.
6. Experimental SAR is the biological anchor.
7. Held-out compounds test the model.
8. Acetate and stereochemical controls test specificity.
9. Boltz-2 affinity is relative predictive evidence, not proof.
10. Clearly distinguish prediction from validation.

## Phase map

| # | Phase | Input | Output | status |
|---|---|---|---|---|
| I | SAR dataset | experimental data + campaign CSVs | discovery/held-out split | `phase01_sar/` |
| II | Receptor ensemble | apo PDBs | compact representative set | `phase02_receptor/` |
| III | Blind pocket discovery | apo PDBs | pocket atlas | `phase03_pockets/` |
| IV | AF3 + Boltz-2 convergence | site_clusters, ligand complexes | convergence map | `phase04_convergence/` |
| V | Ensemble docking | receptor × ligand × pocket | consensus docking model | `phase05_docking/` |
| VI | Interaction fingerprints | complex PDBs | contact matrix | `phase06_fingerprints/` |
| VII | Structural SAR | activity + fingerprints | SAR model | `phase07_sar_model/` |
| VIII | Boltz-2 affinity | complex ensembles | affinity matrix | `phase08_boltz_affinity/` |
| IX | Specificity / falsification | controls | discrimination | `phase09_specificity/` |
| X | Electrostatic analysis | PDBs + APBS | E-vs-recognition | `phase10_electrostatics/` |
| XI | Robustness without MD | seeds/states | consensus metrics | `phase11_robustness/` |
| XII | α9α10 vs α7 vs α4β2 | subtype PDBs | selectivity map | `phase12_subtype/` |
| XIII | Ryanodine | ryanodine complexes | site analysis | `phase13_ryanodine/` |
| XIV | Integrated evidence | all phases | confidence score | `phase14_evidence/` |
| XV | Pharmacophore | SAR + IFP + pocket | pharmacophore | `phase15_pharmacophore/` |
| XVI | Molecule generation | pharmacophore | library | `phase16_molgen/` |
| XVII | Cascading filter | library | 10-30 candidates | `phase17_cascade/` |
| XVIII | Boltz-2 guided design | survivors | complex + affinity | `phase18_boltz_design/` |
| XIX | α9α10 selectivity | candidates | selectivity scores | `phase19_selectivity/` |
| XX | Multi-objective ranking | all metrics | top 10-30 | `phase20_ranking/` |
| XXI | Prospective validation | ranking | experiment design | `phase21_validation/` |
| XXII | Closed-loop | results | updated model | `phase22_closed_loop/` |

## File convention

* `<site>_contacts.csv` — residue-level contact data.
* `<site>_ifp.npy` — interaction fingerprint matrix (numPy binary).
* `<site>_binding.csv` — site × ligand × stoichiometry binding summary.
* `job_tracker.csv` — every computational job recorded.

## Data sources

| source | path | contents |
|---|---|---|
| Campaign | `data/campaign/` | site_clusters_v3, FINAL_SITE_RANKING, LIGAND_COMPARISON, structlib.py, mutagenesis |
| AF3 outputs | `data/allostery/af3_outputs/` | 468 CIF models |
| Boltz-2 outputs | `data/allostery/boltz2/outputs/` | 90 CIF models |
| Deliverable | `data/deliverable/05_structures/` | apo refs, site PDBs, markers |
| MD (ECD) | `data/allostery/md/` | NOT AVAILABLE — no MD output exists (all runs timed out, 0 DCD written) |
