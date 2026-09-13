#!/usr/bin/env python3
"""12. MD decision.

Decides, based on the consensus + plasticity results, whether a molecular
dynamics simulation would resolve a specific open question, or is not
necessary (per prompt section 18).

Output: outputs/12_md_decision.json
"""

import json
from pathlib import Path

from config import RESULTS


def load(o):
    p = RESULTS / o
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return {}


def main():
    plasticity = load("08_pocket_plasticity.json")
    consensus = load("10_residue_consensus.json")
    conv = load("06_three_method_convergence.json")

    open_questions = []
    if (plasticity.get("mean_pocket_ca_std_A") or 0) > 3.0:
        open_questions.append(
            "High receptor-conformation variance in the Boltz-2 ensemble "
            "(ca std %.1f A) leaves the pocket-open/closed state unresolved." % (plasticity.get("mean_pocket_ca_std_A") or 0))
    if not (consensus.get("tier_1") or []):
        open_questions.append("No Tier-1 residue reached; residue prioritization is not yet decisive.")
    if (conv or {}).get("siteaf3_status", "").startswith("NOT RUN"):
        open_questions.append("SiteAF3 layer not yet run; three-method convergence incomplete.")

    decision = {
        "necessary": False,
        "open_questions": open_questions,
        "decision": "MD not necessary at this stage"
        if not open_questions
        else "MD deferred: pending SiteAF3 results and tighter residue prioritization.",
        "rationale": [
            "No experimental structure to benchmark MD against.",
            "Primary actionable output (mutagenesis residue list) is decided by "
            "structural+SAR convergence, not by MD.",
            "MD would not alter residue prioritization unless a specific "
            "pocket-open/closed or cryptic-pocket question can be posed.",
        ],
    }
    if open_questions:
        decision["deferred_reason"] = open_questions

    with open(RESULTS / "12_md_decision.json", "w") as f:
        json.dump(decision, f, indent=2)

    for q in open_questions:
        print("[open] " + q)
    print(decision["decision"])
    print(f"wrote {RESULTS / '12_md_decision.json'}")


if __name__ == "__main__":
    main()