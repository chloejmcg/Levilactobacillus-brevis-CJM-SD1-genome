#!/usr/bin/env bash
# 01_read_qc.sh — raw ONT read QC for DLVP55_2 (seqkit)
# Usage: bash 01_read_qc.sh <reads.fastq.gz>
set -euo pipefail

READS=$1
THREADS=24
OUT=/workdir/cjm423/DLVP55_2/qc
export PATH=/programs/seqkit-2.13.0:$PATH
mkdir -p "$OUT"

# Summary: read count, total bases, N50, quality metrics
seqkit stats -a -T -j $THREADS "$READS" > "$OUT/raw_seqkit_stats.tsv"

# One row per read: id, length, mean quality (input for the R plots)
seqkit fx2tab -n -i -l -q -j $THREADS "$READS" > "$OUT/per_read_len_qual.tsv"

seqkit version > "$OUT/versions.txt"
echo "Done. Results in $OUT"
