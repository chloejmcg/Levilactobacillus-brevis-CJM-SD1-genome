#!/usr/bin/env bash
# 05b_replicon_typing.sh — plasmid replicon typing with ABRicate/PlasmidFinder
# Usage: bash 05b_replicon_typing.sh <assembly.fasta>
set -euo pipefail

GENOME=$(realpath "$1")
OUT=/workdir/cjm423/DLVP55_2/annotation/replicons
SIF=/programs/abricate-1.0.0/abricate.sif
mkdir -p "$OUT"

# Default thresholds
singularity exec --bind /workdir,/local/workdir "$SIF" \
  abricate --db plasmidfinder "$GENOME" > "$OUT/plasmidfinder_default.tsv"

# Relaxed thresholds: Gram-positive/lactobacillus replicons are divergent from
# the (mostly Enterobacteriaceae + Staph/Entero) PlasmidFinder database
singularity exec --bind /workdir,/local/workdir "$SIF" \
  abricate --db plasmidfinder --minid 50 --mincov 40 "$GENOME" \
  > "$OUT/plasmidfinder_relaxed.tsv"

singularity exec --bind /workdir,/local/workdir "$SIF" abricate --version > "$OUT/versions.txt"
echo "Done. Results in $OUT"
