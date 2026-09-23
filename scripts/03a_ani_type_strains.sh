#!/usr/bin/env bash
# 03a_ani_type_strains.sh — ANI of DLVP55_2 vs type-strain genomes of every
# Levilactobacillus species (species boundary: ANI >= 95%)
# Usage: bash 03a_ani_type_strains.sh <assembly.fasta>
set -euo pipefail

GENOME=$(realpath "$1")
THREADS=24
OUT=/workdir/cjm423/DLVP55_2/qc_assembly/skani
DATASETS=/programs/ncbi_datasets_18.35.0/datasets
SKANI=/programs/skani-0.3.1/skani
mkdir -p "$OUT"
cd "$OUT"

# 1) Download type-strain genomes for the whole genus from NCBI
$DATASETS download genome taxon Levilactobacillus --from-type \
  --include genome --filename type_refs.zip
unzip -o -q type_refs.zip -d type_refs
find type_refs -name "*.fna" > ref_list.txt
echo "Type-strain genomes downloaded: $(wc -l < ref_list.txt)"

# 2) ANI of our genome against every reference
$SKANI dist -t $THREADS -q "$GENOME" --rl ref_list.txt -o skani_raw.tsv

# 3) Sort best hit first; keep ANI, aligned fractions, reference name
{ head -1 skani_raw.tsv; tail -n +2 skani_raw.tsv | sort -t$'\t' -k3,3gr; } > skani_vs_type_strains.tsv
$SKANI --version > versions.txt
echo "Done."
