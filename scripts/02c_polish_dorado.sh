#!/usr/bin/env bash
# 02c_polish_dorado.sh — polish the Flye assembly with dorado polish (--bacteria)
# Usage: bash 02c_polish_dorado.sh <draft.fasta> <reads.fastq.gz>
set -euo pipefail

DRAFT=$1
READS=$2
THREADS=24
OUT=/workdir/cjm423/DLVP55_2/assembly/dorado_polish
DORADO=/programs/dorado-1.1.1-linux-x64/bin/dorado
# Basecaller model, taken from the FASTQ read headers (basecall_model_version_id=...)
MODEL=dna_r10.4.1_e8.2_400bps_sup@v4.3.0
export PATH=/programs/samtools-1.20/bin:$PATH
mkdir -p "$OUT" "$OUT/models"

# 1) Align reads to the draft, then sort
$DORADO aligner -t $THREADS "$DRAFT" "$READS" \
  | samtools sort -@ 8 -o "$OUT/reads_to_draft.noRG.bam"

# 2) Add a read group recording the basecaller model (FASTQ input loses this;
#    dorado polish reads it to choose the matching polishing model)
samtools addreplacerg -@ 8 \
  -r "@RG\tID:DLVP55_2\tSM:DLVP55_2\tDS:basecall_model=${MODEL}" \
  -o "$OUT/reads_to_draft.bam" "$OUT/reads_to_draft.noRG.bam"
samtools index "$OUT/reads_to_draft.bam"
rm "$OUT/reads_to_draft.noRG.bam"

# 3) Polish. --bacteria: model tuned for bacterial genomes/plasmids (incl. methylation)
#    --vcf: also write a list of every change made to the draft
$DORADO polish --bacteria --vcf -t $THREADS --device cpu \
  --models-directory "$OUT/models" -o "$OUT" \
  "$OUT/reads_to_draft.bam" "$DRAFT"

$DORADO --version > "$OUT/versions.txt" 2>&1
echo "Done. Outputs in $OUT"
