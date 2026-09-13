#!/usr/bin/env python3
"""6. AF3 + Boltz-2 + SiteAF3 three-method structural convergence.

Merges:
  * AF3 cofolding contacts    (data/allostery/af3_outputs *.cif)
  * Boltz-2 panel contacts    (outputs/05_boltz2_consensus.json)
  * SiteAF3 site-conditioned  (results/siteaf3/*) -- only if already run

Distinguishes model agreement (same method, multiple seeds) from
independent evidence (different methods agreeing on the same region).

n_methods_agree counts the *independent methods* that support a residue:
AF3 contacts present, Boltz-2 cf>0, SiteAF3 contact present. A residue
listed in a candidate pocket definition is NOT by itself a method (avoid
circularity) -- candidate-site membership is reported separately.

Output: outputs/06_three_method_convergence.json
"""

import json
from collections import defaultdict

from config import RESULTS, AF3_OUTPUTS, BOLTZ2_OUTPUTS, POCKET_RESIDUES, CHAIN_SUBUNIT

RESULTS_DIR = RESULTS  # alias


def af3_contacts(limit=200):
    """Minimal residue-level contacts from AF3 CIF outputs.

    Converts AF3 chain IDs to subunits using the condition's stoichiometry so
    the residue keys match the Boltz-2 layer. Uses only root-level
    per-condition *_model.cif files (not per-seed replicates).
    """
    contacts = defaultdict(int)
    try:
        from Bio.PDB.MMCIFParser import MMCIFParser
    except Exception:
        return contacts
    parser = MMCIFParser(QUIET=True)
    cif_files = []
    if AF3_OUTPUTS.exists():
        for cif in sorted(AF3_OUTPUTS.rglob("*_model.cif")):
            if cif.parent.name not in ("seed-1", "seed-2", "seed-3", "seed-4",
                                       "seed-5", "sample-0", "sample-1"):
                cif_files.append(cif)
    cif_files = cif_files[:limit]
    for cif in cif_files:
        try:
            structure = parser.get_structure("x", str(cif))
        except Exception:
            continue
        name = str(cif)
        stoich = "2to3" if "2to3" in name else ("3to2" if "3to2" in name else None)
        chain_sub = CHAIN_SUBUNIT.get(stoich, {}) if stoich else {}
        model = next(iter(structure))
        try:
            import numpy as np
            prot = []
            lig = []
            for chain in model:
                for res in chain:
                    for atom in res:
                        if chain.id in "ABCDE":
                            prot.append((chain.id, res.id[1], atom.coord.copy()))
                        else:
                            lig.append((chain.id, res.id[1], atom.coord.copy()))
            if not lig or not prot:
                continue
            lig = np.vstack([l[2] for l in lig])
            all_c = np.vstack([p[2] for p in prot])
            dist = np.linalg.norm(all_c[:, None, :] - lig[None, :, :], axis=2).min(axis=1)
            for i, (ch, rn, coord) in enumerate(prot):
                if dist[i] <= 4.5:
                    subunit = chain_sub.get(ch, ch)
                    contacts[(subunit, rn)] += 1
        except Exception:
            continue
    return contacts


def siteaf3_contacts():
    """Parse SiteAF3 outputs if present: results/siteaf3/<name>/...model.cif."""
    contacts = defaultdict(int)
    sa_root = RESULTS / "results" / "siteaf3"
    if not sa_root.exists():
        return contacts
    cifs = list(sa_root.rglob("*.cif"))
    if not cifs:
        return contacts
    try:
        from Bio.PDB.MMCIFParser import MMCIFParser
    except Exception:
        return contacts
    parser = MMCIFParser(QUIET=True)
    for cif in cifs[:200]:
        name = str(cif)
        stoich = "2to3" if "2to3" in name else ("3to2" if "3to2" in name else None)
        chain_sub = CHAIN_SUBUNIT.get(stoich, {}) if stoich else {}
        try:
            structure = parser.get_structure("x", cif)
        except Exception:
            continue
        model = next(iter(structure))
        try:
            import numpy as np
            prot, lig = [], []
            for chain in model:
                for res in chain:
                    for atom in res:
                        if chain.id in "ABCDE":
                            prot.append((chain.id, res.id[1], atom.coord.copy()))
                        else:
                            lig.append((chain.id, res.id[1], atom.coord.copy()))
            if not lig or not prot:
                continue
            lig = np.vstack([l[2] for l in lig])
            all_c = np.vstack([p[2] for p in prot])
            dist = np.linalg.norm(all_c[:, None, :] - lig[None, :, :], axis=2).min(axis=1)
            for i, (ch, rn, coord) in enumerate(prot):
                if dist[i] <= 4.5:
                    subunit = chain_sub.get(ch, ch)
                    contacts[(subunit, rn)] += 1
        except Exception:
            continue
    return contacts


def main():
    site_res = {}
    # iterate in reverse so overlapping pocket residues prefer the favored site 23
    for site in reversed(list(POCKET_RESIDUES)):
        for sub, reslist in POCKET_RESIDUES[site].items():
            for r in reslist:
                site_res[(sub, r)] = site

    af3 = af3_contacts()
    print(f"[1/4] AF3: {len(af3)} contacting residues")

    b2 = {}
    b2_path = RESULTS / "05_boltz2_consensus.json"
    if b2_path.exists():
        with open(b2_path) as f:
            b2_data = json.load(f)
        for row in b2_data["residue_table"]:
            sub, rn = row["residue"].split(".")
            b2[(sub, int(rn))] = row["cf_all"]
    print(f"[2/4] Boltz-2: {len(b2)} contacting residues")

    sa = siteaf3_contacts()
    sa_root = RESULTS / "results" / "siteaf3"
    if len(sa):
        sa_notes = f"parsed {len(sa)} contacting residues from {sa_root}"
    else:
        sa_notes = "NOT RUN — results/siteaf3 empty; fill after 04_siteaf3_submit.sh"
    print(f"[3/4] SiteAF3: {sa_notes}")

    residues = sorted(set(af3) | set(b2) | set(sa), key=lambda k: (k[0], k[1]))
    table = []
    for sub, rn in residues:
        cf_b2 = b2.get((sub, rn), 0.0)
        n_af3 = af3.get((sub, rn), 0)
        n_sa = sa.get((sub, rn), 0)
        site = site_res.get((sub, rn))
        n_methods = (1 if n_af3 > 0 else 0) + (1 if cf_b2 > 0 else 0) + (1 if n_sa > 0 else 0)
        table.append({
            "subunit": sub, "resnum": rn,
            "site_candidate": site,
            "af3_n_contacts": n_af3,
            "boltz2_cf": round(cf_b2, 4),
            "siteaf3_n_contacts": n_sa,
            "n_methods_agree": n_methods,
        })
    table.sort(key=lambda r: (-r["n_methods_agree"], -r["boltz2_cf"], -r["af3_n_contacts"]))

    result = {
        "af3_n_contacts": len(af3),
        "boltz2_n_contacts": len(b2),
        "siteaf3_status": sa_notes,
        "method_note": ("n_methods_agree counts independent methods (AF3/Boltz-2/SiteAF3). "
                        "Candidate-site membership is reported separately and is NOT counted "
                        "as a method (avoids circularity)."),
        "residue_agreement": table,
        "site_hit_summary": {
            str(s): sum(1 for r in table if r["site_candidate"] == s and r["boltz2_cf"] > 0)
            for s in POCKET_RESIDUES
        },
    }

    with open(RESULTS / "06_three_method_convergence.json", "w") as f:
        json.dump(result, f, indent=2)

    print("Residues with AF3 + Boltz-2 agreement (top 20):")
    for r in table[:20]:
        print(f"  {r['subunit']}.{r['resnum']}  site={r['site_candidate']}  "
              f"boltz2_cf={r['boltz2_cf']:.3f}  af3_n={r['af3_n_contacts']}  "
              f"siteaf3_n={r['siteaf3_n_contacts']}  n_methods={r['n_methods_agree']}")
    print(f"wrote {RESULTS / '06_three_method_convergence.json'}")


if __name__ == "__main__":
    main()