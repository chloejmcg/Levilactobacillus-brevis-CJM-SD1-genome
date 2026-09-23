#!/usr/bin/env bash
# 04c_pfam_domains.sh — Pfam domain assignment for selected proteins
# Usage: bash 04c_pfam_domains.sh <proteins.faa> <outdir>
set -euo pipefail

FAA=$(realpath "$1")
OUT=$(realpath "$2")
THREADS=24
PFAM=/programs/PfamScan/db/Pfam-A.hmm
export PATH=/programs/hmmer-3.4/bin:$PATH
mkdir -p "$OUT"

hmmscan --cpu $THREADS -E 1e-5 --domE 1e-5 \
  --tblout "$OUT/pfam_hits.tbl" \
  --domtblout "$OUT/pfam_domains.tbl" \
  "$PFAM" "$FAA" > "$OUT/pfam_full.txt"

hmmscan -h | head -2 > "$OUT/versions.txt"
echo "Done. Results in $OUT"
