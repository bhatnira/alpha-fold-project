#!/bin/bash
# Phase III — Run fpocket on apo structures for blind pocket discovery.
# fpocket must be installed (conda install -c conda-forge fpocket)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"
OUTDIR="$SCRIPT_DIR/results"
mkdir -p "$OUTDIR"

# Find apo PDBs
DELIVERABLE="$ROOT/data/deliverable/05_structures"
PUB_RECEPTOR="$ROOT/../alpha9alpha10_publication/01_receptor"

echo "=== Blind pocket discovery with fpocket ==="
echo ""

# Process deliverable structures
if [ -d "$DELIVERABLE" ]; then
    echo "Processing deliverable structures..."
    for pdb in "$DELIVERABLE"/*.pdb; do
        name=$(basename "$pdb" .pdb)
        echo "  Running fpocket on $name..."
        fpocket -f "$pdb" -o "$OUTDIR/${name}_fpocket" 2>/dev/null || true
    done
fi

# Process publication receptor structures
if [ -d "$PUB_RECEPTOR" ]; then
    echo "Processing publication receptor structures..."
    for pdb in "$PUB_RECEPTOR"/*.pdb; do
        name=$(basename "$pdb" .pdb)
        echo "  Running fpocket on $name..."
        fpocket -f "$pdb" -o "$OUTDIR/${name}_fpocket" 2>/dev/null || true
    done
fi

echo ""
echo "Done. Results in $OUTDIR"
echo ""
echo "fpocket output structure:"
echo "  <name>_fpocket/<name>_out/"
echo "    info.txt          — pocket rankings"
echo "    pocket<n>_atm.pdb — pocket residues"
echo "    pocket<n>_vert.pdb — pocket vertices"
