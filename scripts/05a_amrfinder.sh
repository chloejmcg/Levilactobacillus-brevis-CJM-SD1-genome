#!/usr/bin/env bash
# 05a_amrfinder.sh — antimicrobial resistance gene screen (AMRFinderPlus)
# Note: amrfinder 4.0.23 on BioHPC fails (needs GLIBCXX_3.4.32); 3.12.8 works.
# Database reused from the Bakta light DB download (no separate download needed).
# Usage: bash 05a_amrfinder.sh <assembly.fasta> <proteins.faa>
set -euo pipefail

GENOME=$(realpath "$1")
PROTEINS=$(realpath "$2")
THREADS=24
OUT=/workdir/cjm423/DLVP55_2/annotation/amrfinder
DB=/workdir/cjm423/databases/bakta/db-light/amrfinderplus-db/latest
export PATH=/programs/amrfinder-3.12.8:/programs/hmmer-3.4/bin:/programs/bin/blast+:$PATH
mkdir -p "$OUT"

# Nucleotide mode: search the assembly directly
amrfinder --nucleotide "$GENOME" --database "$DB" --plus \
  --threads $THREADS --name DLVP55_2 > "$OUT/amrfinder_nucleotide.tsv"

# Protein mode: search the predicted proteome
amrfinder --protein "$PROTEINS" --database "$DB" --plus \
  --threads $THREADS --name DLVP55_2 > "$OUT/amrfinder_protein.tsv"

amrfinder --version > "$OUT/versions.txt"
readlink -f "$DB" >> "$OUT/versions.txt"
echo "Done. Results in $OUT"
