#!/usr/bin/env bash
# =============================================================================
# 06_prepare_ncbi_submission.sh
# Build the GenBank-ready FASTA (.fsa) and checksums for L. brevis CJM-SD1
# (Plasmidsaurus order DLVP55, sample 2; renamed 2026-10-08).
#
#   - optionally rotates the chromosome to start at dnaA, using the dnaA
#     coordinates from the Bakta annotation (seqkit restart; no dnaapler needed)
#   - verifies the rotated chromosome is a pure rotation of the original
#   - writes NCBI definition lines with source/topology/plasmid modifiers
#   - runs the pre-submission checks NCBI applies (no terminal Ns, >199 nt,
#     SeqID rules, plasmid-name rules, expected total length)
#   - records MD5s for the .fsa and the raw reads (for SRA)
#
# Run on BioHPC inside tmux:   ROTATE=1 bash scripts/06_prepare_ncbi_submission.sh
# Needs only seqkit (>=2.0, for `restart`). Optional: BAKTA_TSV=<path to Bakta .tsv>
# Set ROTATE=0 to submit the Flye orientation unchanged.
#
# Refs: NCBI genome submission guide, https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/
# =============================================================================
set -euo pipefail

ROTATE="${ROTATE:-1}"

WD=/workdir/cjm423/DLVP55_2
REPO="$HOME/DLVP55_2-genome"
ASM="$REPO/results/DLVP55_2_final_assembly.fasta"
RAW="$WD/raw/DLVP55_2_sample_2.fastq.gz"
RAW_MD5_EXPECTED="8be94936d257c7a4a52c56eaf9c8444f"
OUT="$WD/06_ncbi"
FSA="$OUT/CJM-SD1.fsa"

ORG="Levilactobacillus brevis"
STRAIN="CJM-SD1"
# Input contig names (as in the repo assembly) -> names submitted to NCBI.
# NCBI plasmid names: start with "p", <=20 chars, letters/digits/._ only (no hyphens).
declare -A PNAME=( [pDLVP55_2_1]=pCJMSD1_1 [pDLVP55_2_2]=pCJMSD1_2 )
EXPECTED_TOTAL=2514906

mkdir -p "$OUT" "$REPO/results/06_ncbi"
command -v seqkit >/dev/null || { echo "seqkit not on PATH"; exit 1; }

echo "== 0. Input replicons =="
seqkit fx2tab -n -i -l -g "$ASM" | tee "$OUT/input_replicons.tsv"
ids=$(seqkit seq -n -i "$ASM")
for id in chromosome pDLVP55_2_1 pDLVP55_2_2; do
  # exact-match count; avoid `| grep -q`, which trips pipefail via SIGPIPE
  [[ $(grep -cx "$id" <<< "$ids") -eq 1 ]] || { echo "Missing or duplicated replicon: $id"; exit 1; }
done
[[ $(wc -l <<< "$ids") -eq 3 ]] || { echo "Expected exactly 3 records in $ASM"; exit 1; }

# -----------------------------------------------------------------------------
# 1. Optional dnaA rotation of the chromosome (plasmids left as assembled)
# -----------------------------------------------------------------------------
seqkit grep -p chromosome "$ASM" | seqkit seq -w 0 > "$OUT/chromosome.orig.fasta"
seqkit grep -p pDLVP55_2_1 -p pDLVP55_2_2 "$ASM" | seqkit seq -w 0 > "$OUT/plasmids.fasta"

if [[ "$ROTATE" == "1" ]]; then
  # Rotate so the chromosome starts at the dnaA start codon, using the dnaA
  # coordinates from the existing Bakta annotation (same operation dnaapler
  # performs; avoids an extra dependency). Bakta annotated this exact sequence,
  # so its coordinates apply directly. Override the TSV path with BAKTA_TSV=...
  if [[ -z "${BAKTA_TSV:-}" ]]; then
    mapfile -t cands < <(find "$REPO/results" -name '*.tsv' -size +100k -exec grep -lP '\tdnaA\t' {} + 2>/dev/null || true)
    [[ ${#cands[@]} -eq 1 ]] || { echo "Found ${#cands[@]} Bakta TSVs with dnaA; set BAKTA_TSV=<path>. Candidates: ${cands[*]:-none}"; exit 1; }
    BAKTA_TSV="${cands[0]}"
  fi
  echo "dnaA coordinates from: $BAKTA_TSV"
  # Bakta TSV: Sequence Id, Type, Start, Stop, Strand, Locus Tag, Gene, Product, DbXrefs
  hits=$(awk -F'\t' '!/^#/ && $2=="cds" && $7=="dnaA"' "$BAKTA_TSV")
  [[ $(grep -c . <<< "$hits") -eq 1 ]] || { echo "Expected exactly one dnaA CDS in $BAKTA_TSV, found:"; echo "$hits"; exit 1; }
  read -r d_seq d_start d_stop d_strand d_tag < <(awk -F'\t' '{print $1, $3, $4, $5, $6}' <<< "$hits")
  echo "dnaA: $d_tag  $d_seq:$d_start-$d_stop ($d_strand)" | tee "$OUT/dnaA_rotation.txt"
  L=$(seqkit fx2tab -n -l "$OUT/chromosome.orig.fasta" | awk '{print $NF}')
  (( d_stop <= L )) || { echo "dnaA end ($d_stop) beyond chromosome length ($L) — wrong TSV?"; exit 1; }
  if [[ "$d_strand" == "+" ]]; then
    seqkit restart -i "$d_start" "$OUT/chromosome.orig.fasta"
  else
    # minus strand: reverse-complement, so dnaA reads forward; its start codon
    # (Bakta 'Stop' coordinate) lands at position L - stop + 1
    seqkit seq -r -p -t dna "$OUT/chromosome.orig.fasta" | seqkit restart -i $(( L - d_stop + 1 ))
  fi | seqkit seq -w 0 | seqkit replace -p '^.*$' -r 'chromosome' > "$OUT/chromosome.fasta"

  # Confirm the new start is really dnaA: an ORF from position 1 with a start
  # codon, a terminal stop and no internal stops (fails if coordinates don't match).
  glen=$(( d_stop - d_start + 1 ))
  aa=$(seqkit subseq -r 1:"$glen" "$OUT/chromosome.fasta" | seqkit translate -T 11 -w 0 | grep -v '^>')
  codon=$(seqkit subseq -r 1:3 "$OUT/chromosome.fasta" | seqkit seq -s -u)
  [[ "$codon" =~ ^(ATG|GTG|TTG)$ ]] || { echo "Position 1 is $codon, not a start codon — rotation wrong"; exit 1; }
  [[ "${aa: -1}" == "*" && "${aa%?}" != *"*"* ]] \
    || { echo "dnaA ORF check FAILED (internal stop or no terminal stop)"; exit 1; }
  echo "dnaA ORF check: ${codon} start, $(( ${#aa} - 1 )) aa, no internal stops — OK" | tee -a "$OUT/dnaA_rotation.txt"

  # Verify: same length, and rotated seq occurs in the doubled original
  # (forward or reverse complement, since rotation flips strand when dnaA is on -).
  # grep -F on files is fast for a single 2.4-Mb pattern; bash globbing is not.
  seqkit seq -s -w 0 -u "$OUT/chromosome.orig.fasta" > "$OUT/orig.txt"
  seqkit seq -s -w 0 -u "$OUT/chromosome.fasta"      > "$OUT/rot.txt"
  seqkit seq -r -p -t dna -s -w 0 -u "$OUT/chromosome.fasta" > "$OUT/rot_rc.txt"
  len_o=$(tr -d '\n' < "$OUT/orig.txt" | wc -c)
  len_r=$(tr -d '\n' < "$OUT/rot.txt"  | wc -c)
  [[ "$len_o" -eq "$len_r" ]] || { echo "Rotation changed length ($len_o -> $len_r)!"; exit 1; }
  printf '%s%s\n' "$(cat "$OUT/orig.txt")" "$(cat "$OUT/orig.txt")" > "$OUT/doubled.txt"
  if grep -qF -f "$OUT/rot.txt" "$OUT/doubled.txt"; then
    echo "Rotation check: forward rotation of original — OK"
  elif grep -qF -f "$OUT/rot_rc.txt" "$OUT/doubled.txt"; then
    echo "Rotation check: reverse-complement rotation of original — OK"
  else
    echo "Rotation check FAILED: output is not a rotation of the input"; exit 1
  fi
  rm -f "$OUT/doubled.txt"
  echo "First 60 nt (should start with the dnaA start codon): $(head -c 60 "$OUT/rot.txt")"
else
  cp "$OUT/chromosome.orig.fasta" "$OUT/chromosome.fasta"
  echo "ROTATE=0: chromosome left in Flye orientation"
fi

# -----------------------------------------------------------------------------
# 2. Write NCBI definition lines
#    SeqIDs: <50 chars, [A-Za-z0-9_.:*#-]. Plasmid names: start with 'p',
#    <=20 chars, no word 'plasmid'. Bacterial single chromosome: [location=chromosome]
#    with no chromosome name.
# -----------------------------------------------------------------------------
{
  seqkit replace -p '^.*$' \
    -r "chromosome [organism=${ORG}] [strain=${STRAIN}] [location=chromosome] [topology=circular] [completeness=complete]" \
    "$OUT/chromosome.fasta"
  for p in pDLVP55_2_1 pDLVP55_2_2; do
    np=${PNAME[$p]}
    seqkit grep -p "$p" "$OUT/plasmids.fasta" | seqkit replace -p '^.*$' \
      -r "${np} [organism=${ORG}] [strain=${STRAIN}] [plasmid-name=${np}] [topology=circular] [completeness=complete]"
  done
} | seqkit seq -w 80 > "$FSA"

echo "== 2. Definition lines =="
grep '^>' "$FSA"

# -----------------------------------------------------------------------------
# 3. Pre-submission checks
# -----------------------------------------------------------------------------
echo "== 3. Checks =="
fail=0
seqkit fx2tab -n -i -l -g -B N "$FSA" | tee "$OUT/final_replicons.tsv"

total=$(seqkit stats -T "$FSA" | awk 'NR==2{print $5}')
[[ "$total" -eq "$EXPECTED_TOTAL" ]] && echo "Total length $total — OK" \
  || { echo "Total length $total != $EXPECTED_TOTAL"; fail=1; }

# any N at all (closed genome should have none)
nN=$(seqkit fx2tab -n -B N "$FSA" | awk '{s+=$NF} END{print s+0}')
[[ "$nN" == "0" ]] && echo "No Ns — OK" || { echo "Contains Ns (${nN}% in some record)"; fail=1; }

# length > 199 nt
seqkit fx2tab -n -i -l "$FSA" | awk '$2<200{print "Too short: "$1; bad=1} END{exit bad}' \
  && echo "All records >199 nt — OK" || fail=1

# SeqID characters and length
grep '^>' "$FSA" | awk '{id=substr($1,2); if (length(id)>=50 || id !~ /^[A-Za-z0-9_.:*#-]+$/) {print "Bad SeqID: "id; bad=1}} END{exit bad}' \
  && echo "SeqIDs valid — OK" || fail=1

# plasmid-name rules
grep -o 'plasmid-name=[^]]*' "$FSA" | cut -d= -f2 | while read -r pn; do
  if [[ ! "$pn" =~ ^p[A-Za-z0-9._]+$ || ${#pn} -gt 20 || "$pn" == *plasmid* ]]; then
    echo "Bad plasmid name: $pn"; exit 1
  fi
done && echo "Plasmid names valid — OK" || fail=1

[[ "$fail" -eq 0 ]] || { echo "One or more checks failed — fix before submitting."; exit 1; }

# -----------------------------------------------------------------------------
# 4. Checksums (assembly + raw reads for SRA)
# -----------------------------------------------------------------------------
echo "== 4. Checksums =="
gzip -t "$RAW" && echo "Raw reads gzip integrity — OK"
raw_md5=$(md5sum "$RAW" | cut -d' ' -f1)
[[ "$raw_md5" == "$RAW_MD5_EXPECTED" ]] && echo "Raw MD5 matches notes — OK" \
  || { echo "Raw MD5 $raw_md5 != $RAW_MD5_EXPECTED"; exit 1; }

md5sum "$FSA" "$RAW" | tee "$OUT/md5sums.txt"

# Copy small outputs into the repo (reads stay out of git)
cp "$FSA" "$OUT/final_replicons.tsv" "$OUT/md5sums.txt" "$REPO/results/06_ncbi/"
[[ -f "$OUT/dnaA_rotation.txt" ]] && cp "$OUT/dnaA_rotation.txt" "$REPO/results/06_ncbi/"

echo
echo "Done. Upload to the Genome Submission Portal: $FSA"
echo "Upload to SRA: $RAW  (MD5 $raw_md5)"
