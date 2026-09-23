#!/usr/bin/env bash
# 02b_assemble_flye.sh — de novo assembly of filtered ONT reads with Flye
# Usage: bash 02b_assemble_flye.sh <filtered_reads.fastq.gz>
set -euo pipefail

READS=$1
THREADS=24
OUT=/workdir/cjm423/DLVP55_2/assembly/flye
export PATH=/programs/Flye-2.9.2/bin:$PATH
mkdir -p "$(dirname "$OUT")"

# --nano-hq : mode for high-accuracy reads (R10.4.1 / Q20+ basecalling)
# Flye detects circular contigs (chromosome, plasmids) automatically
flye --nano-hq "$READS" --out-dir "$OUT" --threads $THREADS

flye --version > "$OUT/versions.txt"
echo "Done. Key outputs: assembly.fasta, assembly_info.txt, assembly_graph.gfa"
