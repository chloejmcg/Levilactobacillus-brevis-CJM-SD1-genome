#!/usr/bin/env bash
# 06c_screen_p1_targets.sh — pick a pDLVP55_2_1-SPECIFIC target for the P1R pair.
#
# Why this exists: DLVP55_12230 (RepA) is NOT p1-specific. Primers designed in
# it also amplify pDLVP55_2_2, because p2 carries a truncated Rep that shares
# sequence with p1's functional one (the in-trans replication dependency
# predicted in the EPS doc, B5). A p1 primer that also hits p2 would score a
# culture that had LOST p1 as positive — the exact failure the assay exists to
# detect.
#
# This screens every non-transposase p1 CDS for single-replicon uniqueness by
# blasting the whole gene against the assembly and counting where it lands.
#
#   cd ~/DLVP55_2-genome
#   bash scripts/06c_screen_p1_targets.sh
#
set -euo pipefail

REPO="${HOME}/DLVP55_2-genome"
ASM="${ASM:-${REPO}/results/DLVP55_2_final_assembly.fasta}"
BAKTA_TSV="${BAKTA_TSV:-${REPO}/results/04_annotation/DLVP55_2.tsv}"
OUT="${REPO}/results/06_benchwork/primer_design"
mkdir -p "${OUT}"

add_programs_dir() {
  command -v "$2" >/dev/null 2>&1 && return 0
  local d
  for d in $(compgen -G "$1" 2>/dev/null | sort -Vr); do
    [[ -x "${d}/$2"     ]] && { export PATH="${d}:${PATH}";     return 0; }
    [[ -x "${d}/bin/$2" ]] && { export PATH="${d}/bin:${PATH}"; return 0; }
  done
  return 1
}
add_programs_dir '/programs/seqkit-*'    seqkit     || { echo "no seqkit" >&2; exit 1; }
add_programs_dir '/programs/ncbi-blast*' makeblastdb || { echo "no blast" >&2; exit 1; }

# Candidate CDS on p1: 800-2500 bp, excluding transposases/IS and Rep itself.
echo "== candidate p1 CDS =="
awk -F'\t' '$1=="pDLVP55_2_1" && $2=="cds" {
  len=$4-$3+1
  if (len>=800 && len<=2500 &&
      $8 !~ /[Tt]ransposase/ && $8 !~ /IS[0-9]/ &&
      $6 != "DLVP55_12230")
    printf "%s\t%s\t%s\t%s\t%s\n", $6, $3, $4, len, $8
}' "${BAKTA_TSV}" | tee "${OUT}/p1_candidates.tsv"

[[ -s "${OUT}/p1_candidates.tsv" ]] || { echo "no candidates" >&2; exit 1; }

# Extract each gene.
: > "${OUT}/p1_candidates.fa"
while IFS=$'\t' read -r LT START END LEN PROD; do
  seqkit subseq --chr pDLVP55_2_1 -r "${START}:${END}" "${ASM}" \
    | seqkit replace -p '.*' -r "${LT}" >> "${OUT}/p1_candidates.fa"
done < "${OUT}/p1_candidates.tsv"

# Blast whole genes against the assembly. A p1-specific gene aligns to
# pDLVP55_2_1 and nowhere else at any meaningful length.
makeblastdb -in "${ASM}" -dbtype nucl -out "${OUT}/asm_db" -title DLVP55_2 >/dev/null
blastn -query "${OUT}/p1_candidates.fa" -db "${OUT}/asm_db" \
  -evalue 1e-5 -num_threads "${THREADS:-24}" \
  -outfmt '6 qseqid sseqid pident length qlen' \
  -out "${OUT}/p1_candidates_blast.tsv"

echo
echo "== uniqueness screen =="
echo "A clean target shows ONE row: pDLVP55_2_1, near-full-length, ~100% id."
echo
awk -F'\t' '{
  # only alignments covering >=20% of the gene count as "lands here"
  if ($4 >= 0.2*$5) {
    key=$1"\t"$2
    n[key]++
    if ($4>best[key]) { best[key]=$4; pid[key]=$3 }
  }
} END {
  for (k in n) printf "%s\t%d hit(s)\tbest_len=%d\tpid=%.1f\n", k, n[k], best[k], pid[k]
}' "${OUT}/p1_candidates_blast.tsv" | sort

echo
echo "== verdict =="
awk -F'\t' '{ if ($4 >= 0.2*$5) { c[$1]++; if ($2!="pDLVP55_2_1") off[$1]++ } }
END { for (g in c) printf "%s\t%s\n", g, (off[g] ? "SHARED — reject ("off[g]" off-replicon)" : "p1-SPECIFIC — usable") }' \
  "${OUT}/p1_candidates_blast.tsv" | sort -k2,2 -k1,1

echo
echo "Pick a 'p1-SPECIFIC — usable' gene >= 800 bp, then rerun the primer"
echo "design with P1R retargeted:"
echo "  P1R_LOCUS=<locus_tag> bash scripts/06c_primer_design.sh"
