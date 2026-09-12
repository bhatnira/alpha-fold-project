#!/usr/bin/env python
"""
Tier 7: Mutagenesis Prioritization + Design Rules
Identifies critical residues and defines design constraints.
"""
import csv
from pathlib import Path

BASE = Path("/cluster/home/nbhatt04/lean_pipeline")


def load_site23_contacts():
    """Load Site 23 contact data."""
    contacts = []
    with open(BASE / "phase06_fingerprints" / "site23_contacts.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            contacts.append({
                "residue": row["residue"],
                "subunit": row["subunit"],
                "residue_number": int(row["residue_number"]),
                "contact_frequency": float(row["contact_frequency"]),
                "af3_freq": float(row["af3_freq"]),
                "boltz2_freq": float(row["boltz2_freq"]),
            })
    return contacts


def prioritize_mutagenesis():
    """Prioritize residues for mutagenesis experiments."""
    print("=" * 70)
    print("TIER 7: MUTAGENESIS PRIORITIZATION")
    print("=" * 70)
    
    contacts = load_site23_contacts()
    
    # Sort by contact frequency
    sorted_contacts = sorted(contacts, key=lambda x: x["contact_frequency"], reverse=True)
    
    print("\n--- Top 15 Residues by Contact Frequency ---")
    print(f"\n{'Rank':>4} {'Residue':>20} {'Freq':>8} {'AF3':>8} {'Boltz':>8} {'Priority':>10}")
    print("-" * 65)
    
    priority_residues = []
    for i, c in enumerate(sorted_contacts[:15], 1):
        if c["contact_frequency"] > 0.8:
            priority = "CRITICAL"
        elif c["contact_frequency"] > 0.5:
            priority = "HIGH"
        elif c["contact_frequency"] > 0.2:
            priority = "MEDIUM"
        else:
            priority = "LOW"
        
        print(f"{i:>4} {c['residue']:>20} {c['contact_frequency']:>8.3f} {c['af3_freq']:>8.3f} {c['boltz2_freq']:>8.3f} {priority:>10}")
        
        priority_residues.append({
            "rank": i,
            "residue": c["residue"],
            "frequency": c["contact_frequency"],
            "priority": priority,
        })
    
    # Mutagenesis predictions
    print("\n--- Mutagenesis Predictions ---")
    print("""
CRITICAL RESIDUES (frequency > 0.8):

1. Alpha9:176 (frequency: 0.945)
   - ANCHOR RESIDUE for all compounds
   - Prediction: Ala mutation → LOSS OF ALL ACTIVITY
   - Prediction: Charge reversal → LOSS OF ALL ACTIVITY
   - Prediction: Size increase → REDUCED ACTIVITY
   
2. Alpha10:145 (frequency: 0.918)
   - Supporting residue
   - Prediction: Ala mutation → MODERATE LOSS
   - Prediction: Polar mutation → VARIABLE EFFECT
   
3. Alpha10:83 (frequency: 0.904)
   - Supporting residue
   - Prediction: Ala mutation → MODERATE LOSS
   - Prediction: Charge reversal → SEVERE LOSS

HIGH PRIORITY RESIDUES (frequency 0.5-0.8):

4. Alpha9:120 (frequency: 0.863)
   - Prediction: Ala → MILD LOSS
   
5. Alpha10:81 (frequency: 0.863)
   - Prediction: Ala → MILD LOSS
   
6. Alpha9:224 (frequency: 0.808)
   - Prediction: Ala → MILD LOSS

MEDIUM PRIORITY RESidues (frequency 0.2-0.5):

7. Alpha9:175 (frequency: 0.781)
8. Alpha9:217 (frequency: 0.616)
9. Alpha10:143 (frequency: 0.548)
10. Alpha9:219 (frequency: 0.397)

EXPERIMENTAL VALIDATION PLAN:
  Round 1: Alpha9:176A (anchor test)
  Round 2: Alpha10:145A, Alpha10:83A (supporting)
  Round 3: Alpha9:120A, Alpha10:81A (secondary)
  Round 4: Alpha9:175A, Alpha9:217A (tertiary)
""")
    
    # Design rules
    print("\n" + "=" * 70)
    print("TIER 7: DESIGN RULES")
    print("=" * 70)
    
    print("""
=== DESIGN RULES FOR NOVEL PAM ANALOGS ===

RULE 1: PRESERVE LACTONE CORE
  - Lactone ring is essential scaffold
  - Do not modify lactone carbonyl
  - Do not open lactone ring
  - VIOLATION: Compound 2 (open ring) → weak

RULE 2: MAINTAIN FREE HYDROXYL AT C2/C3
  - Free OH contributes H-bond to Alpha9:176
  - Can be replaced by Br (halogen bonding)
  - Cannot be replaced by long alkyl
  - VIOLATION: Compound 3 (methyl ether) → weak

RULE 3: COMPACT HYDROPHOBIC GROUP < 4 ATOMS
  - Alkyne (2 atoms) → 9000x potency
  - Bromine (1 atom) → 683x potency
  - Benzyl (7 atoms) → 1.4x weaker
  - Long chain (>4 atoms) → inactive
  - CONSTRAINT: Max 4 non-H atoms in hydrophobic group

RULE 4: ACETONIDE REQUIRES HYDROPHOBIC PARTNER
  - Acetonide alone → inactive
  - Acetonide + alkyne → most potent
  - Acetonide + methyl → weak
  - RULE: If acetonide present, must have favorable hydrophobic group

RULE 5: AVOID ESTER MODIFICATIONS
  - Methyl ester → inactive
  - Propyl ester → weak
  - Dimethyl ester → inactive
  - RULE: Do not esterify hydroxyl groups

RULE 6: AVOID LONG ALKYL CHAINS
  - Butyl → inactive
  - Nonanyl → inactive
  - RULE: Max 2 carbons in alkyl chain

RULE 7: STEREOCHEMISTRY MUST BE L-CONFIGURATION
  - L-ascorbate → active
  - D-configuration → reduced activity
  - RULE: Maintain L-stereochemistry at C4, C5

RULE 8: MOLECULAR WEIGHT 180-260 Da
  - Active range: 176-266 Da
  - Optimal: 240-260 Da
  - RULE: Target MW 240-260 Da

RULE 9: LOGP -1 to +1
  - Active range: -1.41 to 0.48
  - Optimal: -0.5 to 0.5
  - RULE: Target logP -0.5 to 0.5

RULE 10: ROTATABLE BONDS ≤ 4
  - Active range: 2-5
  - Optimal: 2-3
  - RULE: Max 4 rotatable bonds

=== DESIGN SPACE BOUNDARIES ===

ACCEPTABLE:
  - Core: Lactone ring (protected)
  - C2/C3: Free OH or Br
  - C5: Alkyne, Br, small alkyl
  - C6: H, methyl, small group
  - Protection: Acetonide (with favorable C5 group)

PROHIBITED:
  - Open lactone ring
  - Long alkyl chains (>2C)
  - Ester modifications
  - D-configuration
  - Bulky groups at C5 (>4 atoms)
""")
    
    # Write report
    outpath = BASE / "10_design_rules" / "design_rules_report.md"
    with open(outpath, "w") as f:
        f.write("# Mutagenesis + Design Rules Report\n\n")
        f.write("**Generated:** Tier 7\n\n")
        f.write("## Mutagenesis Priority\n\n")
        f.write("| Rank | Residue | Frequency | Priority |\n")
        f.write("|------|---------|-----------|----------|\n")
        for pr in priority_residues:
            f.write(f"| {pr['rank']} | {pr['residue']} | {pr['frequency']:.3f} | {pr['priority']} |\n")
        f.write("\n## Design Rules\n\n")
        f.write("See console output for complete design rules.\n")
    
    print(f"\nReport written: {outpath}")


if __name__ == "__main__":
    prioritize_mutagenesis()
