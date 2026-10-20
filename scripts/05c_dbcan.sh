#!/usr/bin/env bash
# 05c_dbcan.sh — CAZyme annotation of DLVP55_2 (run_dbcan / dbCAN3)
#
# Purpose: positive evidence for the GH70/GH68-negative claim. Keyword search of
# Bakta output is not sufficient evidence for a manuscript statement.
#
# Run on BioHPC (cbsumm13). Databases go to /workdir (temporary); results and
# the summary tables go to ~/DLVP55_2-genome/results (permanent, git-tracked).
#
#   tmux attach -t asm || tmux new -s asm
#   cd ~/DLVP55_2-genome
#   bash scripts/05c_dbcan.sh 2>&1 | tee results/05_tractability/05c_dbcan.log
#
set -euo pipefail

THREADS="${THREADS:-24}"                      # nproc lies on BioHPC (OMP_NUM_THREADS=1)
USER_ID="$(whoami)"
WORKDIR="/workdir/${USER_ID}/DLVP55_2"
DB_DIR="/workdir/${USER_ID}/databases/dbcan"
REPO="${HOME}/DLVP55_2-genome"
OUT_REPO="${REPO}/results/05_tractability/dbcan"
OUT_WORK="${WORKDIR}/dbcan_out"
VENV="/workdir/${USER_ID}/venv_dbcan"

mkdir -p "${DB_DIR}" "${OUT_REPO}" "${OUT_WORK}" "${REPO}/results/05_tractability"

say() { printf '\n== %s ==\n' "$*"; }

# ---------------------------------------------------------------------------
# 1. Locate the proteome. Bakta .faa, NOT the PGAP file (not back yet) and NOT
#    the rotated NCBI FASTA. Protein sequences are rotation-invariant, so the
#    GH70/GH68 absence claim is unaffected by the dnaA rotation; only the
#    locus tags will need remapping to PGAP later.
# ---------------------------------------------------------------------------
say "Locating Bakta proteome"
if [[ -n "${FAA:-}" ]]; then
  :
else
  mapfile -t CAND < <(find "${REPO}/results" -name '*.faa' -type f 2>/dev/null | sort)
  if [[ ${#CAND[@]} -eq 0 ]]; then
    echo "ERROR: no .faa found under ${REPO}/results" >&2
    echo "       Bakta protein FASTA may only exist in /workdir. Set FAA=/path/to/x.faa" >&2
    exit 1
  elif [[ ${#CAND[@]} -gt 1 ]]; then
    echo "ERROR: more than one .faa found. Rerun with FAA=<one of these>:" >&2
    printf '  %s\n' "${CAND[@]}" >&2
    exit 1
  fi
  FAA="${CAND[0]}"
fi
echo "Proteome: ${FAA}"
NPROT=$(grep -c '^>' "${FAA}")
echo "Proteins: ${NPROT}"
[[ "${NPROT}" -gt 1000 ]] || { echo "ERROR: only ${NPROT} proteins — wrong file?" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 2. Get run_dbcan. Prefer a system install; otherwise a throwaway venv in
#    /workdir (do not pollute home — it is network storage and slow).
# ---------------------------------------------------------------------------
say "Resolving run_dbcan"
DBCAN_BIN=""
if command -v run_dbcan >/dev/null 2>&1; then
  DBCAN_BIN="$(command -v run_dbcan)"
elif compgen -G "/programs/*dbcan*/bin/run_dbcan" >/dev/null 2>&1; then
  DBCAN_BIN="$(compgen -G '/programs/*dbcan*/bin/run_dbcan' | head -1)"
  export PATH="$(dirname "${DBCAN_BIN}"):${PATH}"
else
  echo "Not on PATH or in /programs — building a venv at ${VENV}"
  python3 -m venv "${VENV}"
  # shellcheck disable=SC1091
  source "${VENV}/bin/activate"
  pip install --upgrade pip >/dev/null
  pip install dbcan
  DBCAN_BIN="$(command -v run_dbcan)"
fi
echo "run_dbcan: ${DBCAN_BIN}"
DBCAN_VER="$(run_dbcan --version 2>&1 | head -1 || true)"
echo "version string: ${DBCAN_VER}"

# DIAMOND is a hard dependency of the diamond module.
command -v diamond >/dev/null 2>&1 || export PATH="/programs/diamond:${PATH}"
command -v diamond >/dev/null 2>&1 && diamond --version || echo "WARNING: diamond not found"

# CLI generation differs: v5.x uses subcommands, v4.x uses positional args.
if run_dbcan --help 2>&1 | grep -qi 'CAZyme_annotation'; then
  GEN=5
else
  GEN=4
fi
echo "CLI generation: ${GEN}"

# ---------------------------------------------------------------------------
# 3. Database. Several GB — this is the step that needs the reservation.
# ---------------------------------------------------------------------------
say "Database"
if [[ -s "${DB_DIR}/.dbcan_db_ok" ]]; then
  echo "Already staged at ${DB_DIR} (remove .dbcan_db_ok to force refetch)"
else
  if [[ "${GEN}" -eq 5 ]]; then
    run_dbcan database --db_dir "${DB_DIR}" --aws_s3 || \
    run_dbcan database --db_dir "${DB_DIR}"
  else
    echo "ERROR: v4 CLI needs manual DB staging; see dbcan.readthedocs.io" >&2
    exit 1
  fi
  touch "${DB_DIR}/.dbcan_db_ok"
fi
du -sh "${DB_DIR}" || true

# ---------------------------------------------------------------------------
# 4. Annotate. Protein mode — the proteome is already called by Bakta, so do
#    not let dbCAN re-call genes with pyrodigal and renumber everything.
# ---------------------------------------------------------------------------
say "CAZyme annotation (${NPROT} proteins, ${THREADS} threads)"
if [[ "${GEN}" -eq 5 ]]; then
  run_dbcan CAZyme_annotation \
    --input_raw_data "${FAA}" \
    --mode protein \
    --output_dir "${OUT_WORK}" \
    --db_dir "${DB_DIR}" \
    --methods hmm --methods diamond --methods dbCANsub \
    --threads "${THREADS}"
else
  run_dbcan "${FAA}" protein \
    --out_dir "${OUT_WORK}" \
    --db_dir "${DB_DIR}" \
    --tools all \
    --hmm_cpu "${THREADS}" --dia_cpu "${THREADS}" --dbcan_thread "${THREADS}"
fi

# ---------------------------------------------------------------------------
# 5. The actual question: any GH70 or GH68?
#    GH70 = glucansucrase (dextran/reuteran), GH68 = levansucrase/inulosucrase.
#    Both act on sucrose. Their absence is what predicts "sucrose should not
#    induce EPS" in Gate A, and it is a claim the manuscript makes.
# ---------------------------------------------------------------------------
say "GH70 / GH68 verdict"
OVERVIEW="$(find "${OUT_WORK}" -maxdepth 2 -iname 'overview*.tsv' -o -maxdepth 2 -iname 'overview*.txt' | head -1 || true)"
[[ -n "${OVERVIEW}" ]] || OVERVIEW="$(find "${OUT_WORK}" -maxdepth 2 -iname '*hmm_results.tsv' | head -1 || true)"

{
  echo "# dbCAN run — DLVP55_2"
  echo "date: $(date -Iseconds)"
  echo "host: $(hostname)"
  echo "run_dbcan: ${DBCAN_VER}"
  echo "proteome: ${FAA}"
  echo "proteins: ${NPROT}"
  echo "db_dir: ${DB_DIR}"
  echo "methods: hmm,diamond,dbCANsub (protein mode)"
  echo
  echo "## GH70 / GH68 hits (any method, any file)"
  if grep -rhoE 'GH(70|68)([^0-9]|$)' "${OUT_WORK}" 2>/dev/null | sort | uniq -c | grep . ; then
    echo
    echo "NON-EMPTY — the sucrose-negative prediction needs revisiting. Inspect:"
    grep -rlE 'GH(70|68)([^0-9]|$)' "${OUT_WORK}" 2>/dev/null || true
  else
    echo "none"
    echo
    echo "Negative confirmed by HMM + DIAMOND + dbCAN-sub. This is the citable"
    echo "evidence for the GH70/GH68-negative statement."
  fi
  echo
  echo "## CAZyme family counts"
  if [[ -n "${OVERVIEW}" ]]; then
    echo "source: ${OVERVIEW}"
    awk -F'\t' 'NR>1{print $0}' "${OVERVIEW}" | wc -l | xargs echo "rows:"
  fi
  echo
  echo "## Glycosyltransferase families present (GT*)"
  grep -rhoE 'GT[0-9]+' "${OUT_WORK}" 2>/dev/null | sort | uniq -c | sort -rn || echo "none"
} | tee "${OUT_REPO}/dbcan_summary.txt"

# ---------------------------------------------------------------------------
# 6. Copy results out of /workdir before the reservation ends.
# ---------------------------------------------------------------------------
say "Copying results to repo"
find "${OUT_WORK}" -maxdepth 2 -type f \( -name '*.tsv' -o -name '*.txt' -o -name '*.out' \) \
  -size -20M -exec cp -v {} "${OUT_REPO}/" \;

say "Done"
cat <<EOF

Results:  ${OUT_REPO}
Summary:  ${OUT_REPO}/dbcan_summary.txt

Next:
  cd ${REPO}
  git add results/05_tractability/dbcan scripts/05c_dbcan.sh
  git commit -m "dbCAN CAZyme annotation: GH70/GH68 status confirmed"
  git push

Caveat to record in the manuscript: these are Bakta locus tags (DLVP55_xxxxx).
Remap to PGAP tags when the GenBank annotation returns (SUB16549748).
EOF
