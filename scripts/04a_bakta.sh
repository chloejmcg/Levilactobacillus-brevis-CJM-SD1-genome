#!/usr/bin/env bash
# 04a_bakta.sh — genome annotation with Bakta 1.9.3 (BioHPC container) + light DB
# Usage: bash 04a_bakta.sh <final_assembly.fasta>
set -euo pipefail

GENOME=$(realpath "$1")
THREADS=24
OUT=/workdir/cjm423/DLVP55_2/annotation/bakta
DB=/workdir/cjm423/databases/bakta/db-light
SIF=/programs/bakta-1.9.3/bakta.sif
REPLICONS=$HOME/DLVP55_2-genome/notes/replicons.tsv
mkdir -p "$(dirname "$OUT")"

# Container PATH must include /opt/conda/bin (where bakta lives inside the image)
singularity exec --bind /workdir,/local/workdir \
  --env PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin \
  "$SIF" \
  bakta --db "$DB" \
    --genus Levilactobacillus --species brevis --strain DLVP55_2 \
    --gram '?' --complete --compliant \
    --replicons "$REPLICONS" \
    --locus-tag DLVP55 --prefix DLVP55_2 \
    --threads $THREADS --output "$OUT" --force \
    "$GENOME"

echo "Done. Summary: $OUT/DLVP55_2.txt"
