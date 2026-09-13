#!/usr/bin/env python3
"""11. Mutagenesis plan (final receptor-level validation step).

For each Tier-1/Tier-2 residue defines:
  * wild-type identity (from the verified alpha9 / alpha10 sequences)
  * predicted role (SAR-defining contact vs pocket-lining vs support)
  * a specific mutation (default Ala probe) + a conservative control where useful
  * expected effect if the binding site is correct vs incorrect
  * receptor-function confounders
  * which SAR cliff-pair (if any) implicates the residue

Mutagenesis occurs AFTER all computational + SAR analyses.

Output: outputs/11_mutagenesis_plan.json
"""

import json
from pathlib import Path

from config import RESULTS, POCKET_RESIDUES

# Verified sequences (AGENTS.md: alpha9 = UniProt Q9UGM1 479 aa,
# alpha10 = UniProt Q9GZZ6 450 aa). 1-indexed residue numbers.
SEQ = {
    "alpha9": ("MNWSHSCISFCWIYFAASRLRAAETADGKYAQKLFNDLFEDYSNALRPVEDTDKVLNVTL"
               "QITLSQIKDMDERNQILTAYLWIRQIWHDAYLTWDRDQYDGLDSIRIPSDLVWRPDIVLYNK"
               "ADDESSEPVNTNVVLRYDGLITWDAPAITKSSCVVDVTYFPFDNQQCNLTFGSWTYNGNQVD"
               "IFNALDSGDLSDFIEDVEWEVHGMPAVKNVISYGCCSEPYPDVTFTLLLKRRSSFYIVNLLIP"
               "CVLISFLAPLSFYLPAASGEKVSLGVTILLAMTVFQLMVAEIMPASENVPLIGKYYIATMAL"
               "ITASTALTIMVMNIHFCGAEARPVPHWARVVILKYMSRVLFVYDVGESCLSPHHSRERDHLT"
               "KVYSKLPESNLKAARNKDLSRKKDMNKRLKNDLGCQGKNPQEAESYCAQYKVLTRNIEYIAK"
               "CLKDHKATNSKGSEWKKVAKVIDRFFMWIFFIMVFVMTILIIARAD"),
    "alpha10": ("MGLRSHHLSLGLLLLFLLPAECLGAEGRLALKLFRDLFANYTSALRPVADTDQTLNVTLE"
                "VTLSQIIDMDERNQVLTLYLWIRQEWTDAYLRWDPNAYGGLDAIRIPSSLVWRPDIVLYNKA"
                "DAQPPGSASTNVVLRHDGAVRWDAPAITRSSCRVDVAAFPFDAQHCGLTFGSWTHGGHQLDV"
                "RPRGAAASLADFVENVEWRVLGMPARRRVLTYGCCSEPYPDVTFTLLLRRRAAAYVCNLLLP"
                "CVLISLLAPLAFHLPADSGEKVSLGVTVLLALTVFQLLLAESMPPAESVPLIGKYYMATMTM"
                "VTFSTALTILIMNLHYCGPSVRPVPAWARALLLGHLARGLCVRERGEPCGQSRPPELSPSPQ"
                "SPEGGAGPPAGPCHEPRCLCRQEALLHHVATIANTFRSHRAAQRCHEDWKRLARVMDRFFLA"
                "IFFSMALVMSLLVLVQAL"),
}

CONFOUNDER = ("residue may alter receptor expression/assembly, surface delivery, or "
              "ACh gating independent of PAM binding (control: confirm unchanged ACh "
              "EC50 and Imax in the mutant before interpreting PAM readouts)")

# type -> (probe, rationale)
PROBE = {
    "D": ("Ala", "removes negative charge"),
    "E": ("Ala", "removes negative charge"),
    "K": ("Ala", "removes positive charge"),
    "R": ("Ala", "removes positive charge"),
    "H": ("Ala", "removes titratable/aromatic side chain"),
    "F": ("Ala", "removes aromatic contact"),
    "Y": ("Ala", "removes aromatic/polar contact"),
    "W": ("Ala", "removes aromatic contact"),
    "S": ("Ala", "removes polar side chain"),
    "T": ("Ala", "removes polar side chain"),
    "N": ("Ala", "removes polar amide"),
    "Q": ("Ala", "removes polar amide"),
    "C": ("Ala", "removes thiol"),
    "M": ("Ala", "removes hydrophobic thioether"),
    "A": ("Gly", "removes methyl (size probe)"),
    "G": ("Ala", "adds methyl (size probe)"),
    "V": ("Ala", "reduces branching"),
    "I": ("Ala", "reduces size"),
    "L": ("Ala", "reduces size"),
    "P": ("Ala", "removes kink"),
}

CONSERVATIVE = {
    "F": "Tyr", "Y": "Phe", "W": "Phe",
    "D": "Glu", "E": "Asp",
    "K": "Arg", "R": "Lys",
    "S": "Thr", "T": "Ser",
    "N": "Asp", "Q": "Glu",
    "V": "Leu", "I": "Leu", "L": "Ile",
    "H": "Asn",
}


def load_consensus():
    p = RESULTS / "10_residue_consensus.json"
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return {"rows": [], "tier_1": [], "tier_2": [], "tier_3": []}


def load_cliffs():
    p = RESULTS / "09_sar_activity_cliff_validation.json"
    if not p.exists():
        return {}
    with open(p) as f:
        return json.load(f)


def main():
    cons = load_consensus()
    rows = cons.get("rows", [])
    cliffs = load_cliffs()

    # residue -> cliff pairs that implicate it
    residue_cliffs = {}
    for pair in cliffs.get("pairs", []):
        for res, v in pair.get("defining_residues_sorted", []):
            residue_cliffs.setdefault(res, []).append(
                f"{pair['pair']} (dCF {v:+.2f})")
    for res in cons.get("cliff_residues", []):
        residue_cliffs.setdefault(res, [])

    site_of = {}
    for site, subs in POCKET_RESIDUES.items():
        for sub, reslist in subs.items():
            for r in reslist:
                site_of[f"{sub}.{r}"] = site

    plan = []
    for r in rows:
        if r["tier"] not in (1, 2):
            continue
        sub, rn = r["residue"].split(".")
        rn = int(rn)
        wt = SEQ.get(sub, "")
        aa = wt[rn - 1] if wt and 0 < rn <= len(wt) else "?"
        probe, rationale = PROBE.get(aa, ("Ala", "generic size probe"))
        conservative = CONSERVATIVE.get(aa, None)
        site = r.get("site") or site_of.get(r["residue"])
        sar_defining = bool(r.get("sar_support"))
        cliff_pairs = residue_cliffs.get(r["residue"], [])
        if sar_defining:
            role = "SAR-defining contact residue"
        elif site is not None:
            role = "candidate pocket-lining residue"
        else:
            role = "predicted contact residue (site-agnostic)"
        mutations = [f"{sub}.{rn}:{aa}{probe}"]
        if conservative:
            mutations.append(f"{sub}.{rn}:{aa}{conservative} (conservative control)")
        plan.append({
            "residue": r["residue"],
            "subunit": sub,
            "resnum": rn,
            "wt_aa": aa,
            "probe_mutation": f"{aa}{rn}{probe}",
            "probe_rationale": f"{probe}: {rationale}",
            "proposed_mutations": mutations,
            "resident_site": site,
            "tier": r["tier"],
            "boltz2_cf": r.get("boltz2_cf"),
            "af3_n": r.get("af3_n"),
            "dcf": r.get("dcf"),
            "sar_defining": sar_defining,
            "predicted_role": role,
            "cliff_pairs": cliff_pairs,
            "expected_if_site_correct": (
                "reduced/lost PAM potentiation of ACh-evoked current with unchanged "
                "ACh EC50/Imax (site-correct prediction)"),
            "expected_if_site_incorrect": "no PAM-effect change (WT-like potentiation)",
            "receptor_function_confounders": CONFOUNDER,
        })

    n_t1 = len([r for r in rows if r["tier"] == 1])
    n_t2 = len([r for r in rows if r["tier"] == 2])

    with open(RESULTS / "11_mutagenesis_plan.json", "w") as f:
        json.dump({
            "note": ("Mutagenesis is the final receptor-level validation step, run "
                     "AFTER all computational/SAR work. Mutate, express, confirm "
                     "unchanged ACh EC50/Imax in the mutant, then run the PAM "
                     "potentiation assay; compare WT vs mutant."),
            "order": ["1. confirm receptor expression + ACh response in mutant",
                      "2. measure ACh dose-response (EC50, Imax)",
                      "3. measure PAM potentiation with the nano/PAM and reference "
                      "ascorbate",
                      "4. report delta vs WT; loss-of-PAM-effect = support; "
                      "no-change = falsification at that residue"],
            "tier1_count": n_t1,
            "tier2_count": n_t2,
            "plan": plan,
        }, f, indent=2)

    print(f"tier1={n_t1}, tier2={n_t2}")
    t1 = [p for p in plan if p["tier"] == 1]
    print(f"\nTier-1 mutation plan ({len(t1)}):")
    for p in t1:
        print(f"  {p['wt_aa']}{p['resnum']} ({p['residue']}) -> {p['probe_mutation']}  "
              f"site={p['resident_site']}  dcf={p['dcf']}  role={p['predicted_role']}"
              + (f"  cliffs={p['cliff_pairs']}" if p["cliff_pairs"] else ""))
    print(f"\nplan residues: {[p['residue'] for p in plan]}")
    print(f"wrote {RESULTS / '11_mutagenesis_plan.json'}")


if __name__ == "__main__":
    main()