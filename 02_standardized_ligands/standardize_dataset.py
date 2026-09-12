#!/usr/bin/env python3
"""Tier 1.1: Create standardized master dataset.

Standardizes all 30 compounds from modulator-dataset-a9a10.csv with:
- Original and standardized SMILES
- Isomeric SMILES
- Stereochemistry
- Molecular formula
- Molecular weight
- Experimental activity
- Potency units
- % potentiation
- Compound class
- Parent relationship
- Matched-pair relationship
"""

import csv, json
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("WARNING: RDKit not available. Using manual SMILES parsing.")

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "02_standardized_ligands"
OUT.mkdir(exist_ok=True)


# ─── Compound Definitions ───────────────────────────────────────────
# From modulator-dataset-a9a10.csv + manual curation

COMPOUNDS = [
    {
        "id": "1",
        "original_smiles": "C1([C@@H]([C@H](O)CO)OC(=O)C=1O)O",
        "name": "L-ascorbate (parent)",
        "parent": "self",
        "matched_pair": [],
        "stereochem": "L",
        "class": "active",
        "activity_um": 1797.0,
        "potentiation_pct": 286.0,
    },
    {
        "id": "2",
        "original_smiles": "C(O)C(C1OC(=O)C(O)=C1O)O",
        "name": "L-ascorbate analog (open ring)",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "active",
        "activity_um": 1316.0,
        "potentiation_pct": 300.0,
    },
    {
        "id": "3",
        "original_smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2O)CO1)C",
        "name": "Methyl ether analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "weak",
        "activity_um": 6077.0,
        "potentiation_pct": 293.0,
    },
    {
        "id": "4",
        "original_smiles": "CCCCCCCCC1O[C@H](C2OC(=O)C(O)=C2O)CO1",
        "name": "Nonanyl chain analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "5",
        "original_smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OC)CO1)C",
        "name": "Methyl ester analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "6",
        "original_smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OCC2C=CC=CC=2)CO1)C",
        "name": "Benzyl ether analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "7",
        "original_smiles": "C=CCOC1C([C@@H]2OC(C)(C)OC2)OC(=O)C=1O",
        "name": "Allyl acetonide analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "8",
        "original_smiles": "CCCOC1C(C2OC(C)(C)OC2)OC(=O)C=1O",
        "name": "Propyl acetonide analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "9",
        "original_smiles": "CCCCOC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O",
        "name": "Butyl acetonide analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "10",
        "original_smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OCC2CC2)CO1)C",
        "name": "Cyclopropylmethyl analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "11",
        "original_smiles": "CC1(O[C@H]([C@H]2OC(=O)C(O)=C2OCCC2CCCCC2)CO1)C",
        "name": "Cyclohexylpropyl analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "12",
        "original_smiles": "C(#C)COC1[C@@H]([C@@H]2OC(C)(C)OC2)OC(=O)C=1O",
        "name": "Alkyne acetonide (HIGHLY POTENT)",
        "parent": "1",
        "matched_pair": ["1", "7"],
        "stereochem": "L",
        "class": "potent",
        "activity_um": 0.1981,
        "potentiation_pct": 150.0,
    },
    {
        "id": "13",
        "original_smiles": "CC1(O[C@H]([C@H]2OC(=O)C(OC)=C2OC)CO1)C",
        "name": "Dimethyl ester analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "14",
        "original_smiles": "CC1(O[C@H]([C@H]2OC(=O)C(OCC3C=CC=CC=3)=C2OCC2C=CC=CC=2)CO1)C",
        "name": "Dibenzyl analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "15",
        "original_smiles": "C=CCOC1C(C2OC(C)(C)OC2)OC(=O)C=1OCC=C",
        "name": "Diallyl acetonide analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "16",
        "original_smiles": "CCCOC1[C@H]([C@H]2OC(C)(C)OC2)OC(=O)C=1OCCC",
        "name": "Dipropyl acetonide analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "17",
        "original_smiles": "CC1(O[C@H]([C@H]2OC(=O)C(OCC#N)=C2OCC#N)CO1)C",
        "name": "Dicyanomethyl analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "18",
        "original_smiles": "C1C=CC(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)=CC=1",
        "name": "Benzyl analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "weak",
        "activity_um": 1288.0,
        "potentiation_pct": 500.0,
    },
    {
        "id": "19",
        "original_smiles": "C=CCOC1[C@@H]([C@H](O)CO)OC(=O)C=1O",
        "name": "Allyl analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "20",
        "original_smiles": "CCCCOC1[C@@H]([C@H](O)CO)OC(=O)C=1O",
        "name": "Butyl analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "21",
        "original_smiles": "C1C(COC2[C@@H]([C@H](O)CO)OC(=O)C=2O)C1",
        "name": "Cyclopropylmethyl analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "22",
        "original_smiles": "COC1[C@@H]([C@H](O)CO)OC(=O)C=1OC",
        "name": "Dimethyl analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "23",
        "original_smiles": "C(O)[C@H]1OC(=O)[C@H](O)[C@@H](O)[C@@H]1O",
        "name": "Gluconolactone analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "24",
        "original_smiles": "CCCOC1C(=O)O[C@H]([C@H](O)CO)C=1O",
        "name": "Propyl ester analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "weak",
        "activity_um": 1202.0,
        "potentiation_pct": 190.0,
    },
    {
        "id": "25",
        "original_smiles": "C(Br)[C@H]([C@H]1OC(=O)C(O)=C1O)O",
        "name": "Brominated analog (POTENT)",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "potent",
        "activity_um": 2.63,
        "potentiation_pct": 180.0,
    },
    {
        "id": "26",
        "original_smiles": "C(O)[C@H]([C@H]1OC(=O)[C@@H](O)[C@H]1O)O",
        "name": "Gluconic acid analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "27",
        "original_smiles": "CC1(O[C@H]([C@@H]2OC(=O)[C@@H]3OC(O[C@H]23)(C)C)CO1)C",
        "name": "Bicyclic analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "28",
        "original_smiles": "C(O)[C@H]1OC(=O)[C@H](O)[C@@H]1O",
        "name": "Gulonolactone analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "29",
        "original_smiles": "C([C@@H](O)C1OC(=O)[C@@H](O)[C@H]1O)=O",
        "name": "Fructonic acid analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
    {
        "id": "30",
        "original_smiles": "C(O)[C@@H](O)[C@H]1OC2O[C@@H](O[C@@H]2[C@H]1O)C(Cl)(Cl)Cl",
        "name": "Chloralose analog",
        "parent": "1",
        "matched_pair": ["1"],
        "stereochem": "L",
        "class": "inactive",
        "activity_um": 0.0,
        "potentiation_pct": 0.0,
    },
]


def standardize_smiles(smiles):
    """Standardize SMILES using RDKit."""
    if not RDKIT_AVAILABLE:
        return smiles, smiles, smiles
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles, smiles, smiles
    
    # Standardize
    Chem.SanitizeMol(mol)
    
    # Canonical SMILES (non-isomeric)
    canonical = Chem.MolToSmiles(mol, isomericSmiles=False)
    
    # Isomeric SMILES
    isomeric = Chem.MolToSmiles(mol, isomericSmiles=True)
    
    return smiles, canonical, isomeric


def compute_properties(smiles):
    """Compute molecular properties."""
    if not RDKIT_AVAILABLE:
        return {}
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    
    return {
        "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "formal_charge": Chem.GetFormalCharge(mol),
        "qed": round(Descriptors.qed(mol), 3),
        "tpsa": round(Descriptors.TPSA(mol), 2),
    }


def classify_compound(compound):
    """Classify compound based on activity."""
    ic50 = compound["activity_um"]
    
    if ic50 == 0:
        return "inactive"
    elif ic50 < 10:
        return "potent"
    elif ic50 < 1000:
        return "moderate"
    else:
        return "weak"


def main():
    print("Tier 1.1: Creating Standardized Master Dataset")
    print("=" * 60)
    
    master_records = []
    
    for compound in COMPOUNDS:
        print(f"\nProcessing compound {compound['id']}: {compound['name']}")
        
        # Standardize SMILES
        original, canonical, isomeric = standardize_smiles(compound["original_smiles"])
        
        # Compute properties
        props = compute_properties(canonical)
        
        # Classify
        activity_class = classify_compound(compound)
        
        # Compute pIC50
        ic50 = compound["activity_um"]
        pic50 = -math.log10(ic50 * 1e-6) if ic50 > 0 else 0.0
        
        record = {
            # Identity
            "compound_id": compound["id"],
            "name": compound["name"],
            
            # SMILES
            "original_smiles": original,
            "canonical_smiles": canonical,
            "isomeric_smiles": isomeric,
            
            # Properties
            "molecular_formula": props.get("molecular_formula", ""),
            "molecular_weight": props.get("molecular_weight", 0),
            "logp": props.get("logp", 0),
            "hba": props.get("hba", 0),
            "hbd": props.get("hbd", 0),
            "rotatable_bonds": props.get("rotatable_bonds", 0),
            "aromatic_rings": props.get("aromatic_rings", 0),
            "heavy_atoms": props.get("heavy_atoms", 0),
            "formal_charge": props.get("formal_charge", 0),
            "qed": props.get("qed", 0),
            "tpsa": props.get("tpsa", 0),
            
            # Stereochemistry
            "stereochemistry": compound["stereochem"],
            
            # Experimental data
            "activity_um": compound["activity_um"],
            "potentiation_pct": compound["potentiation_pct"],
            "pic50": round(pic50, 3),
            "activity_class": activity_class,
            
            # Relationships
            "parent_compound": compound["parent"],
            "matched_pairs": ";".join(compound["matched_pair"]),
            
            # Classification
            "compound_class": compound["class"],
            
            # Computational predictions (to be filled)
            "af3_model_ids": "",
            "boltz2_model_ids": "",
            "candidate_pocket": "",
            "receptor_stoichiometry": "",
            "receptor_state": "",
            "interaction_fingerprint": "",
            "confidence": "",
            "final_interpretation": "",
        }
        
        master_records.append(record)
        
        print(f"  Formula: {props.get('molecular_formula', 'N/A')}")
        print(f"  MW: {props.get('molecular_weight', 0):.2f}")
        print(f"  Activity: {compound['activity_um']:.2f} μM")
        print(f"  Class: {activity_class}")
    
    # Write master dataset
    output_path = OUT / "master_dataset.csv"
    if master_records:
        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=master_records[0].keys())
            writer.writeheader()
            writer.writerows(master_records)
    
    print(f"\n{'=' * 60}")
    print(f"Master dataset written: {output_path}")
    print(f"Total compounds: {len(master_records)}")
    
    # Summary statistics
    active = sum(1 for r in master_records if r["activity_class"] in ["potent", "moderate"])
    weak = sum(1 for r in master_records if r["activity_class"] == "weak")
    inactive = sum(1 for r in master_records if r["activity_class"] == "inactive")
    
    print(f"\nActivity distribution:")
    print(f"  Potent (IC50 < 10 μM): {sum(1 for r in master_records if r['activity_class'] == 'potent')}")
    print(f"  Moderate (10-1000 μM): {sum(1 for r in master_records if r['activity_class'] == 'moderate')}")
    print(f"  Weak (>1000 μM): {weak}")
    print(f"  Inactive (IC50 = 0): {inactive}")
    
    # Write JSON version
    json_path = OUT / "master_dataset.json"
    with open(json_path, "w") as f:
        json.dump(master_records, f, indent=2)
    
    print(f"\nJSON dataset written: {json_path}")


if __name__ == "__main__":
    import math
    main()
