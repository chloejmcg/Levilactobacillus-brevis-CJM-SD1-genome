#!/usr/bin/env bash
# 02a_filter_reads.sh — keep the longest, highest-quality reads for assembly
# Usage: bash 02a_filter_reads.sh <reads.fastq.gz>
set -euo pipefail

READS=$1
OUT=/workdir/cjm423/DLVP55_2/filtered
export PATH=/programs/Filtlong-0.2.1/bin:/programs/seqkit-2.13.0:$PATH
mkdir -p "$OUT"

# --min_length 1000     drop reads < 1 kb (too short to help resolve repeats)
# --keep_percent 90     discard the worst 10% of bases (scored on length + quality)
# --target_bases 250M   keep at most ~100x of a 2.5 Mb genome
filtlong --min_length 1000 --keep_percent 90 --target_bases 250000000 "$READS" \
  | gzip > "$OUT/DLVP55_2_filtered.fastq.gz"

seqkit stats -a -T "$OUT/DLVP55_2_filtered.fastq.gz" > "$OUT/filtered_seqkit_stats.tsv"
filtlong --version > "$OUT/versions.txt"
echo "Done."
