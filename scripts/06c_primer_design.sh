#!/usr/bin/env bash
# 06c_primer_design.sh — retention multiplex primer design (Step 0.3)
#
# Four pairs, designed from the closed assembly:
#   CH1  DLVP55_03050  chromosome       priming GT       strain ID + positive control   ~150 bp
#   P2M  DLVP55_12405  pDLVP55_2_2      MobA/MobL        plasmid 2 retention            ~325 bp
#   P1R  DLVP55_12230  pDLVP55_2_1      functional RepA  plasmid 1 retention            ~650 bp
#   P2C  DLVP55_12475  pDLVP55_2_2      pyruvyl transferase  p2 confirmation (singleplex) ~475 bp
#
# ⚠ COORDINATE SYSTEM. Bakta locus tags and coordinates refer to the UNROTATED
#   Flye assembly. The chromosome was rotated to dnaA on 2026-10-10 for NCBI
#   (seqkit restart -i 2357098). Using CJM-SD1.fsa with Bakta coordinates would
#   extract the wrong window and the oligos would be worthless. This script
#   reads results/DLVP55_2_final_assembly.fasta and refuses anything else.
#
#   tmux attach -t asm || tmux new -s asm
#   cd ~/DLVP55_2-genome
#   bash scripts/06c_primer_design.sh 2>&1 | tee results/06_benchwork/06c_primer_design.log
#
set -euo pipefail

THREADS="${THREADS:-24}"
REPO="${HOME}/DLVP55_2-genome"
ASM="${ASM:-${REPO}/results/DLVP55_2_final_assembly.fasta}"
OUT="${REPO}/results/06_benchwork/primer_design"
mkdir -p "${OUT}"

say() { printf '\n== %s ==\n' "$*"; }

# ---------------------------------------------------------------------------
# 0. Guard rails on the assembly
# ---------------------------------------------------------------------------
say "Assembly check"
[[ -s "${ASM}" ]] || { echo "ERROR: assembly not found: ${ASM}" >&2; exit 1; }
case "${ASM}" in
  *CJM-SD1.fsa|*06_ncbi*) echo "ERROR: that is the ROTATED NCBI FASTA. Use the Flye assembly." >&2; exit 1;;
esac

# BioHPC keeps tools in versioned directories under /programs with no bin/
# subdirectory and nothing on the default PATH. Resolve by glob, newest first,
# and prefer the version recorded in notes/VERSIONS.md (seqkit 2.13.0 was used
# for the dnaA rotation; `seqkit restart` needs >= 2.0).
add_programs_dir() {  # $1 = glob, $2 = binary name
  command -v "$2" >/dev/null 2>&1 && return 0
  local d
  for d in $(compgen -G "$1" 2>/dev/null | sort -Vr); do
    if [[ -x "${d}/$2" ]]; then export PATH="${d}:${PATH}"; return 0; fi
    if [[ -x "${d}/bin/$2" ]]; then export PATH="${d}/bin:${PATH}"; return 0; fi
  done
  return 1
}
add_programs_dir '/programs/seqkit-*' seqkit || {
  echo "ERROR: seqkit not found under /programs/seqkit-*" >&2; exit 1; }
seqkit version
seqkit stats "${ASM}" | tee "${OUT}/assembly_stats.txt"
TOTAL=$(seqkit stats -T "${ASM}" | awk 'NR>1{s+=$5} END{print s}')
[[ "${TOTAL}" -eq 2514906 ]] || { echo "ERROR: total length ${TOTAL}, expected 2514906" >&2; exit 1; }
echo "Total length 2514906 — OK"

# ---------------------------------------------------------------------------
# 1. Bakta TSV → target coordinates
# ---------------------------------------------------------------------------
say "Bakta annotation table"
if [[ -z "${BAKTA_TSV:-}" ]]; then
  mapfile -t BT < <(find "${REPO}/results" -name '*.tsv' -path '*annotation*' 2>/dev/null | sort)
  [[ ${#BT[@]} -eq 1 ]] || { echo "ERROR: set BAKTA_TSV=. Candidates:" >&2; printf '  %s\n' "${BT[@]:-none}" >&2; exit 1; }
  BAKTA_TSV="${BT[0]}"
fi
echo "Bakta TSV: ${BAKTA_TSV}"

# P1R target is overridable. The default, DLVP55_12230 (RepA), was shown on
# 2026-10-20 to be NON-SPECIFIC: primers in it amplify pDLVP55_2_2 as well,
# because p2's truncated Rep shares sequence with p1's functional one. Use
# 06c_screen_p1_targets.sh to pick a p1-specific replacement and pass it as
# P1R_LOCUS=<locus_tag>.
P1R_LOCUS="${P1R_LOCUS:-DLVP55_12230}"
if [[ "${P1R_LOCUS}" == "DLVP55_12230" ]]; then
  echo "WARNING: P1R is targeting RepA (DLVP55_12230), which is shared with" >&2
  echo "         pDLVP55_2_2. Expect it to fail the specificity filter." >&2
  echo "         Run scripts/06c_screen_p1_targets.sh and set P1R_LOCUS=." >&2
fi

TARGETS=( "CH1:DLVP55_03050:130:180" \
          "P2M:DLVP55_12405:300:360" \
          "P1R:${P1R_LOCUS}:600:700" \
          "P2C:DLVP55_12475:440:510" )

: > "${OUT}/target_windows.tsv"
printf 'pair\tlocus_tag\tcontig\tstart\tend\tstrand\tlen\n' >> "${OUT}/target_windows.tsv"

for T in "${TARGETS[@]}"; do
  IFS=: read -r PAIR LT LO HI <<< "${T}"
  # Bakta TSV: #Sequence Id, Type, Start, Stop, Strand, Locus Tag, Gene, Product
  LINE=$(awk -F'\t' -v lt="${LT}" '$0 !~ /^#/ && $6==lt {print; exit}' "${BAKTA_TSV}")
  [[ -n "${LINE}" ]] || { echo "ERROR: ${LT} not found in ${BAKTA_TSV}" >&2; exit 1; }
  CONTIG=$(cut -f1 <<< "${LINE}"); START=$(cut -f3 <<< "${LINE}")
  END=$(cut -f4 <<< "${LINE}");   STRAND=$(cut -f5 <<< "${LINE}")
  LEN=$(( END - START + 1 ))
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${PAIR}" "${LT}" "${CONTIG}" "${START}" "${END}" "${STRAND}" "${LEN}" \
    >> "${OUT}/target_windows.tsv"
  # The amplicon must fit inside the gene with room for primers at both ends.
  if (( LEN < HI + 80 )); then
    echo "WARNING: ${PAIR} (${LT}) is only ${LEN} bp — tight for a ${HI} bp amplicon." >&2
  fi
  seqkit subseq --chr "${CONTIG}" -r "${START}:${END}" "${ASM}" \
    | seqkit replace -p '.*' -r "${PAIR}" > "${OUT}/${PAIR}.fa"
done
column -t "${OUT}/target_windows.tsv"

# ---------------------------------------------------------------------------
# 2. primer3
#    Multiplex-compatible settings: tight Tm window so all four pairs share one
#    touchdown program; GC clamp; poly-X capped. Tm is primer3's default
#    SantaLucia thermodynamics, which is NOT the vendor's OligoAnalyzer value —
#    do not re-derive the annealing temperature from the vendor spec sheet.
# ---------------------------------------------------------------------------
say "primer3"
# On cbsumm13: /programs/primer3-2.6/primer3_core with primer3_config/ beside
# it. Take the newest version available.
P3=""
if command -v primer3_core >/dev/null 2>&1; then
  P3="$(command -v primer3_core)"
else
  for D in $(compgen -G '/programs/primer3-*' 2>/dev/null | sort -Vr); do
    for CAND in "${D}/primer3_core" "${D}/bin/primer3_core" "${D}/src/primer3_core"; do
      [[ -x "${CAND}" ]] && { P3="${CAND}"; break 2; }
    done
  done
fi
[[ -n "${P3}" ]] || { echo "ERROR: primer3_core not found under /programs/primer3-*" >&2; exit 1; }
echo "primer3_core: ${P3}"
"${P3}" --version 2>&1 | head -1 || true

# The thermodynamic parameter directory must be passed explicitly or primer3
# refuses the *_TH alignment model that the multiplex settings depend on.
P3CFG=""
for C in "$(dirname "${P3}")/primer3_config" \
         "$(dirname "$(dirname "${P3}")")/primer3_config" \
         "$(dirname "$(dirname "${P3}")")/share/primer3/primer3_config"; do
  [[ -d "${C}" ]] && { P3CFG="${C}/"; break; }
done
[[ -n "${P3CFG}" ]] || { echo "ERROR: primer3_config/ not found next to ${P3}" >&2; exit 1; }
echo "primer3_config: ${P3CFG}"

: > "${OUT}/primer3_input.txt"
for T in "${TARGETS[@]}"; do
  IFS=: read -r PAIR LT LO HI <<< "${T}"
  SEQ=$(seqkit fx2tab "${OUT}/${PAIR}.fa" | cut -f2)
  cat >> "${OUT}/primer3_input.txt" <<EOF
SEQUENCE_ID=${PAIR}_${LT}
SEQUENCE_TEMPLATE=${SEQ}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_NUM_RETURN=5
PRIMER_PRODUCT_SIZE_RANGE=${LO}-${HI}
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=58.0
PRIMER_MAX_TM=62.0
PRIMER_PAIR_MAX_DIFF_TM=2.0
PRIMER_MIN_GC=40.0
PRIMER_OPT_GC_PERCENT=50.0
PRIMER_MAX_GC=60.0
PRIMER_GC_CLAMP=1
PRIMER_MAX_POLY_X=4
PRIMER_MAX_SELF_ANY_TH=45.0
PRIMER_MAX_SELF_END_TH=35.0
PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0
PRIMER_PAIR_MAX_COMPL_END_TH=35.0
PRIMER_MAX_HAIRPIN_TH=40.0
PRIMER_SALT_MONOVALENT=50.0
PRIMER_SALT_DIVALENT=3.0
PRIMER_DNTP_CONC=0.2
PRIMER_DNA_CONC=250.0
PRIMER_EXPLAIN_FLAG=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=${P3CFG}
=
EOF
done

"${P3}" --output="${OUT}/primer3_output.txt" "${OUT}/primer3_input.txt"
grep -c 'PRIMER_LEFT_0_SEQUENCE' "${OUT}/primer3_output.txt" | xargs echo "pairs with >=1 candidate:"
grep 'PRIMER_LEFT_EXPLAIN\|PRIMER_RIGHT_EXPLAIN\|PRIMER_PAIR_EXPLAIN' "${OUT}/primer3_output.txt" || true

# ---------------------------------------------------------------------------
# 3. Specificity. Every candidate must hit its own replicon exactly once and
#    nothing else. Short-word blastn — the default megablast word size will
#    miss 20-mers.
# ---------------------------------------------------------------------------
say "BLAST specificity against the assembly"
add_programs_dir '/programs/ncbi-blast*' makeblastdb || {
  echo "ERROR: BLAST+ not found under /programs/ncbi-blast*" >&2; exit 1; }
makeblastdb -version | head -1
makeblastdb -in "${ASM}" -dbtype nucl -out "${OUT}/asm_db" -title DLVP55_2 >/dev/null
echo "blast db built"

awk -F= '/^PRIMER_(LEFT|RIGHT)_[0-9]+_SEQUENCE=/ {
  split($1,a,"_"); print ">" id "_" a[2] "_" a[3] "\n" $2
}
/^SEQUENCE_ID=/ {id=$2}' "${OUT}/primer3_output.txt" > "${OUT}/candidates.fa"
grep -c '^>' "${OUT}/candidates.fa" | xargs echo "candidate oligos:"

blastn -task blastn-short \
  -query "${OUT}/candidates.fa" \
  -db "${OUT}/asm_db" \
  -word_size 7 -evalue 1000 -dust no -num_threads "${THREADS}" \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  -out "${OUT}/candidates_blast.tsv"
wc -l < "${OUT}/candidates_blast.tsv" | xargs echo "blast hit rows:"

say "Done — now parse and choose"
cat <<EOF

Outputs in ${OUT}:
  target_windows.tsv     locus tag -> contig/coords used
  primer3_output.txt     5 candidate pairs per target
  candidates_blast.tsv   specificity hits

Next:
  Rscript scripts/06c_primer_design.R
  # writes notes/primers.tsv and notes/oligo_order.csv
EOF
