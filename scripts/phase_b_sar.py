#!/usr/bin/env python3
"""
Phase B: Experimental SAR Analysis
Analyzes active/inactive compounds WITHOUT structural modeling.
Identifies SAR anchors, MMPs, and generates SAR hypothesis.
"""
import csv
import json
import os
from pathlib import Path
from collections import defaultdict

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors, DataStructs
    from rdkit.Chem import Draw, rdRAlign
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

PIPELINE_DIR = Path("/cluster/home/nbhatt04/lean_pipeline")
INPUT_CSV = PIPELINE_DIR / "modulator-dataset-a9a10.csv"
OUTPUT_DIR = PIPELINE_DIR / "02_sar_analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

def load_dataset():
    compounds = []
    with open(INPUT_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            compounds.append(row)
    return compounds

def compute_fingerprint(smiles):
    if not HAS_RDKIT:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

def compute_similarity(fp1, fp2):
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def analyze_mmp(compounds):
    """Identify matched molecular pairs."""
    if not HAS_RDKIT:
        return []
    
    fps = {}
    for c in compounds:
        fp = compute_fingerprint(c['Smiles'])
        fps[c['Identifier']] = fp
    
    mmps = []
    ids = list(fps.keys())
    for i in range(len(ids)):
        for j in range(i+1, len(ids)):
            sim = compute_similarity(fps[ids[i]], fps[ids[j]])
            if sim > 0.5:  # Similar enough for MMP analysis
                c1 = next(c for c in compounds if c['Identifier'] == ids[i])
                c2 = next(c for c in compounds if c['Identifier'] == ids[j])
                
                act1 = float(c1.get('Activity (uM)', 0) or 0)
                act2 = float(c2.get('Activity (uM)', 0) or 0)
                pot1 = float(c1.get('%Potentiation', 0) or 0)
                pot2 = float(c2.get('%Potentiation', 0) or 0)
                
                active1 = 1 if (act1 > 0 and pot1 > 0) else 0
                active2 = 1 if (act2 > 0 and pot2 > 0) else 0
                
                mmps.append({
                    'pair': f"{ids[i]}-{ids[j]}",
                    'similarity': round(sim, 3),
                    'activity_change': active2 - active1,
                    'potentiation_change': pot2 - pot1,
                    'compound_a': ids[i],
                    'compound_b': ids[j],
                    'smiles_a': c1['Smiles'],
                    'smiles_b': c2['Smiles'],
                })
    
    return sorted(mmps, key=lambda x: x['similarity'], reverse=True)

def main():
    print("=" * 70)
    print("PHASE B: EXPERIMENTAL SAR ANALYSIS")
    print("=" * 70)
    
    compounds = load_dataset()
    print(f"\nLoaded {len(compounds)} compounds")
    
    # Separate active and inactive
    active = []
    inactive = []
    for c in compounds:
        activity = float(c.get('Activity (uM)', 0) or 0)
        potentiation = float(c.get('%Potentiation', 0) or 0)
        if activity > 0 and potentiation > 0:
            active.append(c)
        else:
            inactive.append(c)
    
    print(f"\nActive: {len(active)} | Inactive: {len(inactive)}")
    
    # Potency ranking
    print("\n--- POTENCY RANKING (active compounds) ---")
    ranked = sorted(active, key=lambda x: float(x['Activity (uM)']))
    for i, c in enumerate(ranked, 1):
        act = float(c['Activity (uM)'])
        pot = float(c['%Potentiation'])
        print(f"  {i}. ID {c['Identifier']}: {act:.2f} uM ({pot:.0f}% pot)")
    
    # Strongest active
    strongest = ranked[0]
    print(f"\nStrongest active: ID {strongest['Identifier']} ({float(strongest['Activity (uM)']):.2f} uM)")
    
    # SAR anchors
    print("\n--- SAR ANCHORS ---")
    print(f"  Strongest: ID {strongest['Identifier']} ({float(strongest['Activity (uM)']):.2f} uM)")
    if len(ranked) > 1:
        print(f"  Second: ID {ranked[1]['Identifier']} ({float(ranked[1]['Activity (uM)']):.2f} uM)")
    if len(ranked) > 2:
        print(f"  Weakest active: ID {ranked[-1]['Identifier']} ({float(ranked[-1]['Activity (uM)']):.2f} uM)")
    
    # MMP analysis
    print("\n--- MATCHED MOLECULAR PAIR ANALYSIS ---")
    mmps = analyze_mmp(compounds)
    print(f"  Found {len(mmps)} matched molecular pairs (similarity > 0.5)")
    
    # Show informative MMPs (activity changes)
    informative = [m for m in mmps if m['activity_change'] != 0]
    print(f"  {len(informative)} pairs show activity changes")
    
    for m in informative[:10]:
        print(f"    Pair {m['pair']}: sim={m['similarity']:.3f}, "
              f"activity change={m['activity_change']}, "
              f"potentiation change={m['potentiation_change']:.0f}%")
    
    # Potentiation analysis
    print("\n--- POTENTIATION ANALYSIS ---")
    for c in ranked:
        act = float(c['Activity (uM)'])
        pot = float(c['%Potentiation'])
        ratio = pot / act if act > 0 else 0
        print(f"  ID {c['Identifier']}: {pot:.0f}% potentiation, "
              f"activity={act:.2f} uM, ratio={ratio:.2f}")
    
    # Save results
    sar_report = {
        'total_compounds': len(compounds),
        'active_count': len(active),
        'inactive_count': len(inactive),
        'strongest_active': strongest['Identifier'],
        'strongest_activity_uM': float(strongest['Activity (uM)']),
        'ranked_active': [{'id': c['Identifier'], 
                          'activity_uM': float(c['Activity (uM)']),
                          'potentiation_pct': float(c['%Potentiation'])} for c in ranked],
        'mmp_count': len(mmps),
        'informative_mmp_count': len(informative),
    }
    
    with open(OUTPUT_DIR / "sar_report.json", 'w') as f:
        json.dump(sar_report, f, indent=2)
    
    # Save MMPs
    if mmps:
        with open(OUTPUT_DIR / "mmp_analysis.csv", 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=mmps[0].keys())
            writer.writeheader()
            writer.writerows(mmps)
    
    print("\n" + "=" * 70)
    print("PHASE B COMPLETE")
    print("=" * 70)

if __name__ == "__main__":
    main()
