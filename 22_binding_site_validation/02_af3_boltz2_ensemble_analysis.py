#!/usr/bin/env python3
"""Section 2: AF3/Boltz-2 ensemble convergence analysis.

Assesses global (receptor), local (pocket), pose, interaction and
stoichiometry convergence across the full-panel Boltz-2 runs and the
existing AF3/Boltz-2 cofolding outputs.
"""

import json
from collections import defaultdict
from pathlib import Path

from Bio.PDB import MMCIFParser

from config import (RESULTS, COFOLDING_RESULTS, AF3_OUTPUTS, BOLTZ2_OUTPUTS,
                    STOICHIOMETRIES, CHAIN_SUBUNIT, POCKET_RESIDUES, NANOMOLAR_PAM)
from common import (parse_pdb, split_ligand, parse_ligand_any_chain,
                    find_contacts, ligand_centroid, first_model_pdbs,
                    parse_confidence, _rmsd_align)

RESULTS.mkdir(parents=True, exist_ok=True)


def classify_site(contact_keys, stoich):
    by_subunit = defaultdict(set)
    for (chain, resnum) in contact_keys:
        sub = CHAIN_SUBUNIT[stoich].get(chain)
        if sub:
            by_subunit[sub].add(resnum)
    best = None
    best_score = -1
    for site, res_map in POCKET_RESIDUES.items():
        score = 0
        for sub, residues in res_map.items():
            score += len(by_subunit.get(sub, set()) & set(residues))
        if score > best_score:
            best_score = score
            best = site
    return best, best_score


def analyze_full_panel():
    records = []
    gdir = COFOLDING_RESULTS / "full_panel"
    if not gdir.exists():
        return records
    for cond_dir in sorted(gdir.iterdir()):
        if not cond_dir.is_dir():
            continue
        name = cond_dir.name
        conf = parse_confidence(cond_dir)
        if not conf:
            continue
        pdbs = first_model_pdbs(cond_dir)
        if not pdbs:
            continue
        parts = name.split("_")
        stoich = parts[1]
        site = int(parts[2].replace("site", ""))
        cpd = int(parts[3].replace("cpd", ""))
        models = []
        for pdb in pdbs:
            atoms = parse_pdb(pdb)
            protein, ligand = split_ligand(atoms)
            if not ligand:
                protein, ligand = parse_ligand_any_chain(atoms)
            valid = len(ligand) > 0
            contacts = find_contacts(protein, ligand) if valid else {}
            centroid = ligand_centroid(ligand) if valid else None
            models.append({
                "path": str(pdb),
                "n_ligand_atoms": len(ligand),
                "n_contacts": len(contacts),
                "contacts": [f"{c[0]}.{c[1]}" for c in contacts],
                "contact_keys": [list(c) for c in contacts],
                "centroid": centroid.tolist() if centroid is not None else None,
                "site_hit": classify_site(contacts, stoich)[0] if contacts else None,
                "site_score": classify_site(contacts, stoich)[1] if contacts else 0,
            })
        pose_centroids = [m["centroid"] for m in models if m["centroid"] is not None]
        pose_spread = _pose_spread(pose_centroids)
        records.append({
            "name": name, "stoich": stoich, "site": site, "cpd": cpd,
            "n_models": len(models), "n_valid_models": sum(1 for m in models if m["n_ligand_atoms"] > 0),
            "pose_spread_ang": pose_spread,
            "site_hits": [m["site_hit"] for m in models],
            "models": models,
            "confidence_score": conf.get("confidence_score"),
            "ligand_iptm": conf.get("ligand_iptm"),
        })
    return records


def _pose_spread(centroids):
    if len(centroids) < 2:
        return None
    import itertools
    dists = []
    for a, b in itertools.combinations(centroids, 2):
        import numpy as np
        dists.append(float(np.linalg.norm(np.array(a) - np.array(b))))
    return sum(dists) / len(dists) if dists else None


def analyze_af3_outputs(limit=40):
    records = []
    parser = MMCIFParser(QUIET=True)
    cif_files = sorted(AF3_OUTPUTS.rglob("*_model.cif")) if AF3_OUTPUTS.exists() else []
    for cif in cif_files[:limit]:
        name = cif.parent.name
        try:
            structure = parser.get_structure("x", str(cif))
        except Exception:
            continue
        model = next(iter(structure))
        atoms = []
        for chain in model:
            for res in chain:
                hetero = res.id[0] != " "
                for atom in res:
                    coord = atom.coord.copy()
                    atoms.append({"chain": chain.id, "resnum": res.id[1],
                                  "resname": res.resname.strip(), "atomname": atom.name.strip(),
                                  "coord": coord, "het": hetero})
        protein = [a for a in atoms if a["chain"] in "ABCDE"]
        ligand = [a for a in atoms if a["chain"] not in "ABCDE"]
        valid = len(ligand) > 0
        contacts = {}
        if valid:
            import numpy as np
            lig = np.vstack([a["coord"] for a in ligand])
            res_map = defaultdict(list)
            for a in protein:
                res_map[(a["chain"], a["resnum"])].append(a["coord"])
            for (chain, resnum), coords in res_map.items():
                dist = min(np.linalg.norm(c - lig, axis=1).min() for c in coords)
                if dist <= 4.5:
                    contacts[(chain, resnum)] = float(dist)
        stoich = "2to3" if "2to3" in name else ("3to2" if "3to2" in name else None)
        records.append({
            "name": name, "stoich": stoich,
            "n_ligand_atoms": len(ligand), "n_contacts": len(contacts),
            "contacts": [f"{c[0]}.{c[1]}" for c in contacts],
            "site_hit": classify_site(contacts, stoich)[0] if contacts and stoich else None,
            "centroid": [float(x) for x in np.vstack([a["coord"] for a in ligand]).mean(axis=0)] if valid else None,
        })
    return records


def analyze_existing_boltz2(limit=60):
    records = []
    cifs = sorted(BOLTZ2_OUTPUTS.rglob("*.cif")) if BOLTZ2_OUTPUTS.exists() else []
    parser = MMCIFParser(QUIET=True)
    for cif in cifs[:limit]:
        try:
            structure = parser.get_structure("x", str(cif))
        except Exception:
            continue
        name = cif.parents[1].name
        model = next(iter(structure))
        atoms = []
        for chain in model:
            for res in chain:
                for atom in res:
                    atoms.append({"chain": chain.id, "resnum": res.id[1],
                                  "resname": res.resname.strip(),
                                  "atomname": atom.name.strip(),
                                  "coord": atom.coord.copy()})
        protein = [a for a in atoms if a["chain"] in "ABCDE"]
        ligand = [a for a in atoms if a["chain"] not in "ABCDE"]
        import numpy as np
        valid = len(ligand) > 0
        contacts = {}
        if valid:
            lig = np.vstack([a["coord"] for a in ligand])
            res_map = defaultdict(list)
            for a in protein:
                res_map[(a["chain"], a["resnum"])].append(a["coord"])
            for (chain, resnum), coords in res_map.items():
                dist = min(np.linalg.norm(c - lig, axis=1).min() for c in coords)
                if dist <= 4.5:
                    contacts[(chain, resnum)] = float(dist)
        stoich = "2to3" if "2to3" in name else ("3to2" if "3to2" in name else None)
        records.append({
            "name": name, "stoich": stoich, "n_ligand_atoms": len(ligand),
            "n_contacts": len(contacts),
            "contacts": [f"{c[0]}.{c[1]}" for c in contacts],
            "site_hit": classify_site(contacts, stoich)[0] if contacts and stoich else None,
        })
    return records


def receptor_rmsd_per_condition(records):
    from collections import defaultdict
    by_name = defaultdict(list)
    for r in records:
        for m in r.get("models", []):
            by_name[r["name"]].append(m)
    results = []
    for name, models in by_name.items():
        if len(models) < 2:
            continue
        ca_by_model = []
        for m in models:
            atoms = parse_pdb(m["path"])
            ca = defaultdict(list)
            for a in atoms:
                if a.atomname == "CA" and a.chain in "ABCDE":
                    ca[a.chain].append(a.coord)
            ca_by_model.append(ca)
        chain_rmsds = []
        for chain in "ABCDE":
            ref = np_vstack(ca_by_model[0].get(chain, []))
            if len(ref) == 0:
                continue
            for ca in ca_by_model[1:]:
                mob = np_vstack(ca.get(chain, []))
                if len(mob) == len(ref) and len(ref) > 10:
                    chain_rmsds.append(_rmsd_align(mob, ref))
        if chain_rmsds:
            results.append({"name": name, "receptor_ca_rmsd_mean": sum(chain_rmsds) / len(chain_rmsds),
                            "n_chains": len(chain_rmsds)})
    return results


def np_vstack(coords):
    import numpy as np
    return np.vstack(coords) if coords else np.zeros((0, 3))


def aggregate(records, key="site"):
    agg = defaultdict(lambda: {"n": 0, "hits": 0})
    for r in records:
        val = r.get(key)
        if val is None:
            continue
        for rh, m in zip(r.get("site_hits", []), r.get("models", [])):
            agg[val]["n"] += 1
            if rh is not None and rh == val:
                agg[val]["hits"] += 1
    return {k: {"conditions": v["n"], "self_hits": v["hits"],
                "self_recurrence": round(v["hits"] / v["n"], 3) if v["n"] else 0.0}
            for k, v in agg.items()}


def main():
    out = {}
    print("[1/4] Parsing full-panel Boltz-2 results...")
    full = analyze_full_panel()
    print(f"  {len(full)} completed conditions, "
          f"{sum(r['n_valid_models'] for r in full)} valid ligand models")
    out["full_panel"] = {"n_conditions": len(full), "conditions": full}

    print("[2/4] Parsing existing AF3 outputs...")
    af3 = analyze_af3_outputs()
    print(f"  {len(af3)} AF3 models parsed")
    out["af3_existing"] = {"n_models": len(af3), "models": af3}

    print("[3/4] Parsing existing Boltz-2 outputs...")
    b2 = analyze_existing_boltz2()
    print(f"  {len(b2)} Boltz-2 models parsed")
    out["boltz2_existing"] = {"n_models": len(b2), "models": b2}

    print("[4/4] Computing convergence metrics...")
    receptor = receptor_rmsd_per_condition(full)
    pose_spreads = [r["pose_spread_ang"] for r in full if r["pose_spread_ang"] is not None]
    out["receptor_convergence"] = {
        "n_conditions": len(receptor),
        "chain_ca_rmsd_mean": round(sum(r["receptor_ca_rmsd_mean"] for r in receptor) / len(receptor), 3) if receptor else None,
        "per_condition": receptor,
    }
    out["pose_convergence"] = {
        "n_conditions": len(pose_spreads),
        "ligand_centroid_spread_mean_ang": round(sum(pose_spreads) / len(pose_spreads), 3) if pose_spreads else None,
        "spread_percentiles": _percentiles(pose_spreads),
    }

    site_ranking = defaultdict(lambda: {"conditions": 0, "models": 0,
                                        "site_hits": defaultdict(int)})
    for r in full:
        site_ranking[r["site"]]["conditions"] += 1
        for hit in r["site_hits"]:
            site_ranking[r["site"]]["models"] += 1
            if hit is not None:
                site_ranking[r["site"]]["site_hits"][hit] += 1
    local = {}
    for site, d in site_ranking.items():
        total_hits = sum(d["site_hits"].values())
        local[site] = {
            "conditions": d["conditions"],
            "models": d["models"],
            "self_hits": d["site_hits"].get(site, 0),
            "self_recurrence": round(d["site_hits"].get(site, 0) / total_hits, 3) if total_hits else 0.0,
            "hit_distribution": dict(d["site_hits"]),
        }
    out["local_convergence"] = local

    stoich_persist = defaultdict(lambda: defaultdict(lambda: {"n": 0, "hits": 0}))
    for r in full:
        if r["stoich"] not in STOICHIOMETRIES:
            continue
        for hit in r["site_hits"]:
            stoich_persist[r["site"]][r["stoich"]]["n"] += 1
            if hit == r["site"]:
                stoich_persist[r["site"]][r["stoich"]]["hits"] += 1
    out["stoichiometry_convergence"] = {
        site: {st: {"self_hits": v["hits"], "n_models": v["n"]}
               for st, v in stoich_by_stoich.items()}
        for site, stoich_by_stoich in stoich_persist.items()
    }

    out["nanomolar_pam"] = {
        "cpd": NANOMOLAR_PAM,
        "conditions": [r for r in full if r["cpd"] == NANOMOLAR_PAM],
    }

    af3_by_site = defaultdict(int)
    for r in af3:
        if r["site_hit"]:
            af3_by_site[r["site_hit"]] += 1
    b2_by_site = defaultdict(int)
    for r in b2:
        if r["site_hit"]:
            b2_by_site[r["site_hit"]] += 1
    out["af3_existing_site_recurrence"] = dict(af3_by_site)
    out["boltz2_existing_site_recurrence"] = dict(b2_by_site)

    write(out)


def _percentiles(vals):
    if not vals:
        return {}
    import numpy as np
    v = np.array(vals)
    return {"p25": round(float(np.percentile(v, 25)), 3),
            "median": round(float(np.median(v)), 3),
            "p75": round(float(np.percentile(v, 75)), 3)}


def write(out):
    path = RESULTS / "02_ensemble_convergence.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print("\nSummary")
    print("-------")
    rc = out["receptor_convergence"]
    print(f"Receptor Cα RMSD (model 0 vs 1/2, chain-level): {rc['chain_ca_rmsd_mean']} A"
          if rc["chain_ca_rmsd_mean"] else "Receptor RMSD: n/a")
    print(f"Ligand pose spread (mean centroid dist across models): "
          f"{out['pose_convergence']['ligand_centroid_spread_mean_ang']} A")
    print("Site self-recurrence pct:")
    for site, d in sorted(out["local_convergence"].items()):
        print(f"  site {site:>2}: {d['self_recurrence']*100:6.1f}%  hits={d['self_hits']}/{d['models']}")
    print("AF3 existing site hits:", dict(out["af3_existing_site_recurrence"]))
    print("Boltz-2 existing site hits:", dict(out["boltz2_existing_site_recurrence"]))
    print("Wrote:", path)


if __name__ == "__main__":
    main()