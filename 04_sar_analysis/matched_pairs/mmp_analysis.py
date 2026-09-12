#!/usr/bin/env python
"""
Tier 5: Matched Molecular Pair Analysis
Systematic comparison of structurally similar compound pairs.
"""
import math
from pathlib import Path

BASE = Path("/cluster/home/nbhatt04/lean_pipeline")

COMPOUNDS = {
    1:  {"name": "L-ascorbate (parent)", "ic50": 1797, "pot": 286, "mw": 176.12, "features": ["lactone", "OH", "open_ring"]},
    2:  {"name": "Open ring analog", "ic50": 1316, "pot": 300, "mw": 176.12, "features": ["lactone", "OH", "open_ring"]},
    3:  {"name": "Methyl ether", "ic50": 6077, "pot": 293, "mw": 216.19, "features": ["lactone", "methyl_ether", "acetonide"]},
    12: {"name": "Alkyne acetonide", "ic50": 0.198, "pot": 150, "mw": 254.24, "features": ["lactone", "alkyne", "acetonide"]},
    18: {"name": "Benzyl analog", "ic50": 1288, "pot": 500, "mw": 266.25, "features": ["lactone", "OH", "benzyl"]},
    24: {"name": "Propyl ester", "ic50": 1202, "pot": 190, "mw": 218.20, "features": ["lactone", "OH", "propyl_ester"]},
    25: {"name": "Brominated analog", "ic50": 2.63, "pot": 180, "mw": 239.02, "features": ["lactone", "OH", "bromine"]},
}

MATCHED_PAIRS = [
    {
        "id": "MMP-1",
        "pair": (1, 18),
        "transformation": "+O-benzyl",
        "chemical_change": "Replace H with benzyl at C5-OH",
        "delta_ic50": 18.197,  # ratio
        "delta_pot": 1.748,  # ratio
        "potency_change": "1.4x weaker",
        "potentiation_change": "1.75x increase",
        "interpretation": "Benzyl provides hydrophobic contact but disrupts H-bond network",
    },
    {
        "id": "MMP-2",
        "pair": (1, 25),
        "transformation": "+Br",
        "chemical_change": "Replace OH with Br at C5",
        "delta_ic50": 683.27,
        "delta_pot": 0.629,
        "potency_change": "683x stronger",
        "potentiation_change": "1.6x decrease",
        "interpretation": "Br fills hydrophobic pocket, halogen bonding possible",
    },
    {
        "id": "MMP-3",
        "pair": (1, 12),
        "transformation": "+Alkyne+Acetonide",
        "chemical_change": "Add alkyne + protect hydroxyls as acetonide",
        "delta_ic50": 9075.76,
        "delta_pot": 0.524,
        "potency_change": "9000x stronger",
        "potentiation_change": "1.9x decrease",
        "interpretation": "Alkyne fills hydrophobic subpocket; acetonide improves shape",
    },
    {
        "id": "MMP-4",
        "pair": (1, 3),
        "transformation": "+Methyl ether",
        "chemical_change": "Add methyl ether at C6",
        "delta_ic50": 0.296,
        "delta_pot": 1.024,
        "potency_change": "3.4x weaker",
        "potentiation_change": "Negligible change",
        "interpretation": "Methyl ether adds steric bulk without favorable interactions",
    },
    {
        "id": "MMP-5",
        "pair": (1, 24),
        "transformation": "+Propyl ester",
        "chemical_change": "Replace OH with propyl ester at C3",
        "delta_ic50": 0.667,
        "delta_pot": 0.664,
        "potency_change": "1.5x weaker",
        "potentiation_change": "1.5x decrease",
        "interpretation": "Ester disrupts H-bond; propyl adds some hydrophobic contact",
    },
    {
        "id": "MMP-6",
        "pair": (25, 12),
        "transformation": "Br→Alkyne+Acetonide",
        "chemical_change": "Replace Br with alkyne; add acetonide",
        "delta_ic50": 13.26,
        "delta_pot": 0.833,
        "potency_change": "13x stronger",
        "potentiation_change": "1.2x decrease",
        "interpretation": "Alkyne + acetonide > bromine for potency",
    },
    {
        "id": "MMP-7",
        "pair": (18, 24),
        "transformation": "Benzyl→Propyl ester",
        "chemical_change": "Replace benzyl with propyl ester",
        "delta_ic50": 0.933,
        "delta_pot": 0.38,
        "potency_change": "Similar (1.07x)",
        "potentiation_change": "2.6x decrease",
        "interpretation": "Similar binding but different functional coupling",
    },
]


def analyze_mmp():
    """Perform matched molecular pair analysis."""
    print("=" * 70)
    print("TIER 5: MATCHED MOLECULAR PAIR ANALYSIS")
    print("=" * 70)
    
    print(f"\n{'Pair':<8} {'Transform':<25} {'ΔIC50':>8} {'ΔPot':>8} {'Interpretation'}")
    print("-" * 90)
    
    for mmp in MATCHED_PAIRS:
        c1, c2 = COMPOUNDS[mmp["pair"][0]], COMPOUNDS[mmp["pair"][1]]
        print(f"{mmp['id']:<8} {mmp['transformation']:<25} {mmp['potency_change']:>12} {mmp['potentiation_change']:>12}")
        print(f"         {mmp['interpretation']}")
    
    # Key insights
    print("\n" + "=" * 70)
    print("KEY MMP INSIGHTS")
    print("=" * 70)
    
    print("""
1. HYDROPHOBIC GROUPS:
   - Br: 683x potency increase (favorable)
   - Benzyl: 1.4x decrease (marginally unfavorable)
   - Alkyne: 9000x increase (highly favorable)
   - Methyl: 3.4x decrease (unfavorable)
   
   INSIGHT: Size/shape matters more than hydrophobicity alone.
   Compact groups (Br, alkyne) >> bulky groups (benzyl, methyl)

2. HYDROXYL PROTECTION:
   - Acetonide alone: inactive
   - Acetonide + alkyne: 9000x
   - Acetonide + methyl: 3.4x weaker
   
   INSIGHT: Acetonide improves shape complementarity only when 
   paired with favorable hydrophobic group

3. H-BOND MODIFICATIONS:
   - Free OH → O-benzyl: 1.4x weaker
   - Free OH → Br: 683x stronger
   - Free OH → propyl ester: 1.5x weaker
   
   INSIGHT: Replacing OH is tolerated if new group fills pocket well

4. POTENTIATION vs BINDING:
   - Strongest binder (C12): weakest PAM (150%)
   - Weakest binder (C18): strongest PAM (500%)
   
   INSIGHT: Binding affinity and PAM function are mechanistically distinct
""")
    
    # SAR rules from MMP
    print("=" * 70)
    print("SAR RULES DERIVED FROM MMP")
    print("=" * 70)
    
    print("""
RULE 1: Compact hydrophobic groups enhance potency
  SUPPORTED BY: Br (683x), Alkyne (9000x)
  VIOLATED BY: Benzyl (1.4x weaker), Methyl (3.4x weaker)
  CONSTRAINT: Group must fit pocket without steric clash

RULE 2: Acetonide improves shape complementarity
  SUPPORTED BY: C12 (9000x with acetonide)
  VIOLATED BY: C7-C11 (inactive with acetonide alone)
  CONSTRAINT: Must have favorable hydrophobic group

RULE 3: Free hydroxyl contributes to binding
  SUPPORTED BY: C1 > C3 (removing free OH weakens)
  VIOLATED BY: C12 (no free OH but most potent)
  CONSTRAINT: Alkyne can compensate for lost H-bond

RULE 4: Long alkyl chains are detrimental
  SUPPORTED BY: C4, C20 (inactive with long chains)
  CONSTRAINT: Max chain length ~2 carbons

RULE 5: Ester modifications reduce potency
  SUPPORTED BY: C5, C13 (inactive with esters)
  CONSTRAINT: Avoid ester modifications at C2/C3
""")
    
    # Write report
    outpath = BASE / "04_sar_analysis" / "matched_pairs" / "mmp_report.md"
    with open(outpath, "w") as f:
        f.write("# Matched Molecular Pair Analysis Report\n\n")
        f.write("**Generated:** Tier 5\n\n")
        f.write("## Pairs Analyzed\n\n")
        f.write("| ID | Pair | Transform | ΔIC50 | ΔPot | Interpretation |\n")
        f.write("|----|------|-----------|-------|------|----------------|\n")
        for mmp in MATCHED_PAIRS:
            f.write(f"| {mmp['id']} | C{mmp['pair'][0]}→C{mmp['pair'][1]} | {mmp['transformation']} | {mmp['potency_change']} | {mmp['potentiation_change']} | {mmp['interpretation']} |\n")
        f.write("\n## SAR Rules\n\n")
        f.write("See console output for detailed rules.\n")
    
    print(f"\nReport written: {outpath}")


if __name__ == "__main__":
    analyze_mmp()
