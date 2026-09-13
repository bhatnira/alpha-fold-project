#!/usr/bin/env python3
"""Shared parsing, geometry, and reporting helpers for the validation analysis."""

import json
import dataclasses
from pathlib import Path

import numpy as np

from config import CHAIN_SUBUNIT, ALPHA9_LENGTH, ALPHA10_LENGTH


@dataclasses.dataclass
class Atom:
    chain: str
    resnum: int
    resname: str
    atomname: str
    coord: np.ndarray
    element: str


def parse_pdb(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21].strip()
                resname = line[17:20].strip()
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    resnum = -1
                atomname = line[12:16].strip()
                if line.startswith("ATOM"):
                    coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                else:
                    coord = _parse_hetatm_coords(line)
                element = line[76:78].strip() or atomname[0]
                atoms.append(Atom(chain, resnum, resname, atomname, coord, element))
    return atoms


def _parse_hetatm_coords(line):
    tokens = line.split()
    for i, tok in enumerate(tokens):
        if tok.count(".") == 1 and i + 2 < len(tokens) and all(
            t.count(".") == 1 and t.lstrip("-").replace(".", "", 1).isdigit() for t in tokens[i + 1 : i + 3]
        ):
            try:
                xyz = np.array([float(t) for t in tokens[i : i + 3]])
                if np.isfinite(xyz).all():
                    return xyz
            except ValueError:
                pass
    raise ValueError(f"Could not locate coordinates in HETATM line: {line}")


def split_ligand(atoms, ligand_chain="P", ligand_resname="LIG"):
    protein, ligand = [], []
    for atom in atoms:
        if atom.chain == ligand_chain and atom.resname == ligand_resname:
            ligand.append(atom)
        elif atom.chain in "ABCDE":
            protein.append(atom)
    return protein, ligand


def parse_ligand_any_chain(atoms):
    protein, ligand = [], []
    for atom in atoms:
        if atom.chain in "ABCDE":
            protein.append(atom)
        else:
            ligand.append(atom)
    return protein, ligand


def protein_residues(protein_atoms):
    residues = {}
    for atom in protein_atoms:
        key = (atom.chain, atom.resnum)
        residues.setdefault(key, []).append(atom)
    return residues


def residue_name(chain, resnum, stoich):
    subunit = CHAIN_SUBUNIT[stoich].get(chain, "?")
    if subunit == "alpha9":
        return f"a9.{resnum}"
    if subunit == "alpha10":
        return f"a10.{resnum}"
    return f"{chain}.{resnum}"


def subunit_of(chain, stoich):
    return CHAIN_SUBUNIT[stoich].get(chain, "?")


def ligand_atom_coords(ligand):
    if not ligand:
        return np.zeros((0, 3))
    return np.vstack([a.coord for a in ligand])


def ligand_centroid(ligand):
    coords = ligand_atom_coords(ligand)
    if len(coords) == 0:
        return None
    return coords.mean(axis=0)


def find_contacts(protein_atoms, ligand, cutoff=4.5):
    contacts = {}
    if not ligand:
        return contacts
    lig_coords = ligand_atom_coords(ligand)
    residue_map = protein_residues(protein_atoms)
    for (chain, resnum), atoms in residue_map.items():
        best = np.inf
        for atom in atoms:
            d = np.linalg.norm(atom.coord - lig_coords, axis=1).min()
            best = min(best, d)
        if best <= cutoff:
            contacts[(chain, resnum)] = float(best)
    return contacts


def find_hbond_contacts(protein_atoms, ligand, cutoff=3.5):
    hbonds = {}
    heavy_ligand = [a for a in ligand if a.element != "H"]
    residue_map = protein_residues(protein_atoms)
    for (chain, resnum), atoms in residue_map.items():
        heavy = [a for a in atoms if a.element != "H"]
        if not heavy or not heavy_ligand:
            continue
        best = np.inf
        for a in heavy:
            d = np.linalg.norm(a.coord - ligand_atom_coords(heavy_ligand), axis=1).min()
            best = min(best, d)
        if best <= cutoff:
            hbonds[(chain, resnum)] = float(best)
    return hbonds


def buried_fraction(ligand, protein_atoms, cutoff=4.5):
    if not ligand:
        return 0.0, 0
    lig_coords = ligand_atom_coords(ligand)
    if len(lig_coords) == 0:
        return 0.0, 0
    residues = protein_residues(protein_atoms)
    all_prot = np.vstack([a.coord for atoms in residues.values() for a in atoms])
    if len(all_prot) == 0:
        return 0.0, 0
    buried = 0
    total = len(lig_coords)
    for coord in lig_coords:
        if np.linalg.norm(all_prot - coord, axis=1).min() <= cutoff:
            buried += 1
    return buried / total, total


def count_clashes(protein_atoms, ligand, cutoff=2.2):
    if not ligand:
        return 0
    lig_coords = ligand_atom_coords(ligand)
    residues = protein_residues(protein_atoms)
    all_prot = np.vstack([a.coord for atoms in residues.values() for a in atoms])
    if len(all_prot) == 0:
        return 0
    n = 0
    for coord in lig_coords:
        if np.linalg.norm(all_prot - coord, axis=1).min() < cutoff:
            n += 1
    return n


def parse_confidence(result_dir):
    result_dir = Path(result_dir)
    for p in sorted(result_dir.rglob("confidence_*.json")):
        try:
            with open(p) as f:
                return json.load(f)
        except Exception:
            continue
    return {}


def first_model_pdbs(result_dir):
    result_dir = Path(result_dir)
    pdbs = []
    for p in sorted(result_dir.rglob("*_model_*.pdb")):
        pdbs.append(str(p))
    return pdbs


def best_model_pdb(result_dir):
    pdbs = sorted(Path(result_dir).rglob("*_model_*.pdb"))
    if not pdbs:
        return None
    return str(pdbs[0])


def load_compound_activity():
    from config import load_compounds
    return load_compounds()


def rmsd(coords_a, coords_b):
    if len(coords_a) != len(coords_b) or len(coords_a) == 0:
        return float("nan")
    d = coords_a - coords_b
    return float(np.sqrt((d * d).sum(axis=1).mean()))


def fit_and_rmsd(mobile_coords, ref_coords):
    from scipy.spatial.transform import Rotation
    if len(mobile_coords) != len(ref_coords) or len(mobile_coords) < 3:
        return float("nan")
    return _rmsd_align(mobile_coords, ref_coords)


def _rmsd_align(mobile, ref):
    m = mobile - mobile.mean(axis=0)
    r = ref - ref.mean(axis=0)
    u, _, vt = np.linalg.svd(r.T @ m)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    diag = np.array([[1, 0, 0], [0, 1, 0], [0, 0, d]])
    rot = vt.T @ diag @ u.T
    aligned = m @ rot.T
    diff = aligned - r
    return float(np.sqrt((diff * diff).sum(axis=1).mean()))


def ca_coords_from_ligand(paths, ligand_chain="P", ligand_resname="LIG"):
    parsed = []
    for path in paths:
        atoms = parse_pdb(path)
        protein, ligand = split_ligand(atoms, ligand_chain, ligand_resname)
        if not ligand:
            protein, ligand = parse_ligand_any_chain(atoms)
        ca = {a.chain: [] for a in protein}
        for a in protein:
            if a.atomname == "CA":
                ca[a.chain].append(a.coord)
        parsed.append({"path": path, "ca": ca, "ligand": ligand, "protein": protein})
    return parsed


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    return path