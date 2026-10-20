#!/usr/bin/env Rscript
# 06c_primer_design.R — parse primer3, screen specificity and cross-dimers,
# emit the order sheet for the DLVP55_2 retention multiplex.
#
#   cd ~/DLVP55_2-genome && Rscript scripts/06c_primer_design.R
#
# Writes:
#   notes/primers.tsv        full record: sequence, Tm, GC, amplicon, coords, hits
#   notes/oligo_order.csv    paste-into-vendor sheet (8 multiplex + 2 singleplex)
#   results/06_benchwork/primer_design/qc_report.txt

suppressPackageStartupMessages({
  library(tidyverse)
})

repo    <- path.expand("~/DLVP55_2-genome")
dir_in  <- file.path(repo, "results/06_benchwork/primer_design")
dir_out <- file.path(repo, "notes")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(file.path(dir_in, "primer3_output.txt")))

# --- expected replicon per pair -------------------------------------------
windows <- read_tsv(file.path(dir_in, "target_windows.tsv"), show_col_types = FALSE)

# --- 1. parse boulder-IO ---------------------------------------------------
lines <- read_lines(file.path(dir_in, "primer3_output.txt"))
recs  <- split(lines, cumsum(lag(lines == "=", default = FALSE)))

parse_rec <- function(x) {
  kv <- str_split_fixed(x[str_detect(x, "=")], "=", 2)
  v  <- setNames(kv[, 2], kv[, 1])
  sid <- v[["SEQUENCE_ID"]]
  idx <- as.integer(unique(str_match(names(v), "^PRIMER_LEFT_([0-9]+)_SEQUENCE$")[, 2]))
  idx <- idx[!is.na(idx)]
  if (!length(idx)) return(NULL)
  map_dfr(idx, function(i) {
    g <- function(k, default = NA_character_) {
      key <- sprintf(k, i); if (key %in% names(v)) unname(v[[key]]) else default
    }
    pos <- str_split_fixed(c(g("PRIMER_LEFT_%d"), g("PRIMER_RIGHT_%d")), ",", 2)
    tibble(
      seq_id    = sid,
      pair      = str_split_fixed(sid, "_", 2)[, 1],
      locus_tag = str_remove(sid, "^[A-Z0-9]+_"),
      cand      = i,
      left_seq  = g("PRIMER_LEFT_%d_SEQUENCE"),
      right_seq = g("PRIMER_RIGHT_%d_SEQUENCE"),
      left_tm   = as.numeric(g("PRIMER_LEFT_%d_TM")),
      right_tm  = as.numeric(g("PRIMER_RIGHT_%d_TM")),
      left_gc   = as.numeric(g("PRIMER_LEFT_%d_GC_PERCENT")),
      right_gc  = as.numeric(g("PRIMER_RIGHT_%d_GC_PERCENT")),
      left_pos  = as.integer(pos[1, 1]),
      right_pos = as.integer(pos[2, 1]),
      amplicon  = as.integer(g("PRIMER_PAIR_%d_PRODUCT_SIZE")),
      pair_any  = as.numeric(g("PRIMER_PAIR_%d_COMPL_ANY_TH", "0")),
      pair_end  = as.numeric(g("PRIMER_PAIR_%d_COMPL_END_TH", "0")),
      penalty   = as.numeric(g("PRIMER_PAIR_%d_PENALTY"))
    )
  })
}

cand <- map_dfr(recs, parse_rec)
if (!nrow(cand)) stop("primer3 returned no candidate pairs — check PRIMER_*_EXPLAIN in the log")
message(sprintf("parsed %d candidate pairs across %d targets",
                nrow(cand), n_distinct(cand$pair)))

# --- 2. specificity --------------------------------------------------------
# A primer is specific if exactly one hit is a full-length, zero-mismatch match.
# Partial 5' matches are tolerated; a 3'-anchored partial is not, because that
# is what primes off-target.
bl <- read_tsv(file.path(dir_in, "candidates_blast.tsv"),
               col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                             "gapopen", "qstart", "qend", "sstart", "send",
                             "evalue", "bitscore"),
               show_col_types = FALSE)

oligo_len <- cand %>%
  select(seq_id, cand, left_seq, right_seq) %>%
  pivot_longer(c(left_seq, right_seq), names_to = "side", values_to = "seq") %>%
  mutate(side   = if_else(side == "left_seq", "LEFT", "RIGHT"),
         qseqid = sprintf("%s_%s_%d", seq_id, side, cand),
         qlen   = nchar(seq))

spec <- bl %>%
  inner_join(oligo_len, by = "qseqid") %>%
  mutate(
    full_perfect = (length == qlen & mismatch == 0 & gapopen == 0),
    # Off-target priming risk is a PERFECT 3'-anchored match long enough to
    # extend. A 12-mer with a mismatch is not that: in a 2.5 Mb genome short
    # degenerate matches occur by the dozen, and scoring them as risk rejects
    # every candidate. Require >= 15 nt, zero mismatch, 3' terminus included.
    three_prime_anchored = (qend == qlen & length >= 15 &
                            mismatch == 0 & gapopen == 0)
  ) %>%
  group_by(qseqid, seq_id, cand, side) %>%
  summarise(
    n_perfect   = sum(full_perfect),
    n_anchored  = sum(three_prime_anchored),
    hit_contigs = paste(sort(unique(sseqid[full_perfect])), collapse = ","),
    .groups = "drop"
  )

spec_wide <- spec %>%
  select(seq_id, cand, side, n_perfect, n_anchored, hit_contigs) %>%
  pivot_wider(names_from = side, values_from = c(n_perfect, n_anchored, hit_contigs))

cand <- cand %>%
  left_join(spec_wide, by = c("seq_id", "cand")) %>%
  left_join(windows %>% select(pair, expect_contig = contig), by = "pair") %>%
  mutate(
    ok_unique   = n_perfect_LEFT == 1 & n_perfect_RIGHT == 1,
    ok_replicon = hit_contigs_LEFT == expect_contig &
                  hit_contigs_RIGHT == expect_contig,
    ok_anchor   = n_anchored_LEFT <= 1 & n_anchored_RIGHT <= 1,
    specific    = ok_unique & ok_replicon & ok_anchor
  )

# Diagnostics. A filter that rejects everything must say which criterion did
# it, or the next person spends an hour guessing (this one cost 20 minutes).
diag <- cand %>%
  group_by(pair) %>%
  summarise(
    n            = n(),
    pass_unique  = sum(ok_unique,   na.rm = TRUE),
    pass_replicon= sum(ok_replicon, na.rm = TRUE),
    pass_anchor  = sum(ok_anchor,   na.rm = TRUE),
    pass_all     = sum(specific,    na.rm = TRUE),
    max_anchor_L = max(n_anchored_LEFT,  na.rm = TRUE),
    max_anchor_R = max(n_anchored_RIGHT, na.rm = TRUE),
    .groups = "drop"
  )
message("\n--- specificity diagnostics (candidates passing each criterion) ---")
print(as.data.frame(diag))
message("")

# --- 3. cross-pair 3' dimer screen ----------------------------------------
# primer3 checks within a pair only. In a multiplex, eight oligos share a tube.
# Crude but effective: count reverse-complement matches of the terminal 5 nt of
# one primer anywhere in another. Anything >0 at the 3' terminus gets flagged.
revcomp <- function(s) {
  chartr("ACGT", "TGCA", s) %>% str_split("") %>% map_chr(~ paste(rev(.x), collapse = ""))
}
three_prime_clash <- function(seqs, names, k = 5) {
  out <- tibble()
  for (i in seq_along(seqs)) {
    tail_i <- str_sub(seqs[i], -k)
    rc     <- revcomp(tail_i)
    for (j in seq_along(seqs)) {
      if (i == j) next
      if (str_detect(seqs[j], fixed(rc))) {
        out <- bind_rows(out, tibble(oligo = names[i], clashes_with = names[j],
                                     motif = tail_i))
      }
    }
  }
  out
}

# --- 4. pick the best specific candidate per pair --------------------------
best <- cand %>%
  filter(specific %in% TRUE) %>%
  group_by(pair) %>%
  arrange(penalty, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

missing <- setdiff(unique(cand$pair), best$pair)
if (length(missing)) {
  warning("no strictly specific candidate for: ", paste(missing, collapse = ", "))
}

# Never emit empty order sheets. If a target has no strictly specific pair,
# fall back to the best available and flag it for manual inspection rather
# than writing a header-only CSV that looks like a clean run.
if (length(missing)) {
  fallback <- cand %>%
    filter(pair %in% missing, ok_unique %in% TRUE, ok_replicon %in% TRUE) %>%
    group_by(pair) %>%
    arrange(n_anchored_LEFT + n_anchored_RIGHT, penalty, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  if (nrow(fallback)) {
    message("fallback (unique + correct replicon, but 3'-anchored off-target ",
            "matches present) for: ", paste(fallback$pair, collapse = ", "))
    best <- bind_rows(best, fallback)
  }
}

best <- best %>% mutate(flag = if_else(specific %in% TRUE, "", "REVIEW"))

if (!nrow(best)) {
  stop("No candidate passed even the relaxed filter. Inspect the diagnostics ",
       "table above and results/06_benchwork/primer_design/candidates_blast.tsv")
}

final <- best %>%
  select(pair, locus_tag, replicon = expect_contig, amplicon, flag,
         left_seq, left_tm, left_gc, left_pos,
         right_seq, right_tm, right_gc, right_pos,
         pair_any, pair_end, penalty,
         n_perfect_LEFT, n_perfect_RIGHT,
         n_anchored_LEFT, n_anchored_RIGHT) %>%
  arrange(desc(amplicon))

oligos <- final %>%
  pivot_longer(c(left_seq, right_seq), names_to = "side", values_to = "sequence") %>%
  mutate(side = if_else(side == "left_seq", "F", "R"),
         name = paste0(pair, "_", side)) %>%
  select(name, sequence, pair, side)

clash <- three_prime_clash(oligos$sequence, oligos$name)

# --- 5. write out ----------------------------------------------------------
write_tsv(final, file.path(dir_out, "primers.tsv"))

stopifnot(nrow(oligos) > 0)
order_sheet <- oligos %>%
  transmute(Name = name,
            Sequence = sequence,
            Scale = "25nm",
            Purification = "STD",
            Note = if_else(pair == "P2C",
                           "singleplex confirmation only",
                           "retention multiplex"))
write_csv(order_sheet, file.path(dir_out, "oligo_order.csv"))

report <- c(
  "# Retention multiplex primer design — QC",
  paste("date:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste("candidates parsed:", nrow(cand)),
  paste("targets with a specific pair:", nrow(final), "of", n_distinct(cand$pair)),
  "",
  "## Amplicon ladder (must be resolvable on 2% agarose with a 100 bp ladder)",
  paste0("  ", final$pair, "  ", final$amplicon, " bp  (", final$replicon, ")"),
  "",
  if (nrow(final) >= 2) {
    sep <- min(diff(sort(final$amplicon)))
    c(paste("  minimum separation:", sep, "bp"),
      if (sep < 100)
        "  ⚠ under 100 bp apart — bands will be hard to call. Widen the size ranges."
      else "  separation OK")
  } else "  (fewer than 2 pairs — no separation to assess)",
  "",
  "## Tm spread across all oligos (one touchdown program must serve all)",
  sprintf("  range: %.1f – %.1f °C, spread %.1f °C",
          min(c(final$left_tm, final$right_tm)),
          max(c(final$left_tm, final$right_tm)),
          diff(range(c(final$left_tm, final$right_tm)))),
  if (diff(range(c(final$left_tm, final$right_tm))) > 4)
    "  ⚠ spread over 4 °C — tighten PRIMER_MIN_TM/MAX_TM and re-run"
  else "  spread OK — touchdown 68→58 °C will serve all four pairs",
  "",
  "## Specificity diagnostics (candidates passing each criterion, of n)",
  capture.output(print(as.data.frame(diag))),
  "",
  if (any(final$flag == "REVIEW"))
    paste0("  ⚠ REVIEW: ", paste(final$pair[final$flag == "REVIEW"], collapse = ", "),
           " — unique and on the right replicon, but with perfect 3'-anchored",
           " off-target matches >= 15 nt. Check before ordering.")
  else "  all reported pairs are strictly specific",
  "",
  "## Cross-pair 3' dimer flags",
  if (nrow(clash) == 0) "  none" else
    paste0("  ⚠ ", clash$oligo, " 3' end complements a site in ", clash$clashes_with,
           " (", clash$motif, ")"),
  "",
  "## Reminders",
  "  - Tm here is primer3/SantaLucia. The vendor's OligoAnalyzer value will differ;",
  "    do not re-derive the annealing temperature from the spec sheet.",
  "  - Coordinates are Bakta/unrotated. Remap to PGAP after SUB16549748 returns.",
  "  - Part C singleplex validation is mandatory before the multiplex is trusted.",
  "  - There is no plasmid-free negative control strain. NTC + extraction blank only."
)
writeLines(report, file.path(dir_in, "qc_report.txt"))
cat(paste(report, collapse = "\n"), "\n")

message("\nwrote notes/primers.tsv and notes/oligo_order.csv")
