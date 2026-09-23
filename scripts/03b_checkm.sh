#!/usr/bin/env bash
# 03b_checkm.sh — genome completeness/contamination with CheckM v1.2.4 (lineage workflow)
# Usage: bash 03b_checkm.sh <assembly.fasta>
set -euo pipefail

GENOME=$(realpath "$1")
THREADS=24
OUT=/workdir/cjm423/DLVP55_2/qc_assembly/checkm
DB=/workdir/cjm423/databases/checkm_data

# CheckM install on BioHPC (python 3.9); data path set by env var (no write to /programs)
export PATH=/programs/checkm-1.2.4/bin:/programs/hmmer-3.4/bin:/programs/pplacer-Linux-v1.1.alpha19:/programs/prodigal-2.6.3:$PATH
export PYTHONPATH=/programs/checkm-1.2.4/lib/python3.9/site-packages
export CHECKM_DATA_PATH=$DB

# Download CheckM reference data once (checkm_data_2015_01_16, ~275 MB)
if [ ! -d "$DB/genome_tree" ]; then
  mkdir -p "$DB" && cd "$DB"
  wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
  tar -xzf checkm_data_2015_01_16.tar.gz
fi

# CheckM takes a folder of genomes ("bins"), so put our assembly in its own folder
rm -rf "$OUT" && mkdir -p "$OUT/bins"
cp "$GENOME" "$OUT/bins/DLVP55_2.fasta"

checkm lineage_wf -x fasta -t $THREADS --tab_table \
  -f "$OUT/checkm_results.tsv" "$OUT/bins" "$OUT/lineage_wf"

checkm 2>&1 | grep -m1 "CheckM v" > "$OUT/versions.txt"
echo "Done."
