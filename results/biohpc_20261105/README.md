Targeted genome queries resolving the Gate A Stage 1 design and several Lane B unknowns.
## Inputs

Genome pulled from the public record **GCA_062071506.1** (= BioProject PRJNA1542414,
BioSample SAMN64125442, GenBank JCFRYO000000000). Verified on download: 3 replicons,
2,514,906 bp, 2,391 CDS, 2,365 proteins (= 2,391 − 26 pseudogenes). Matches Table 1.

NCBI ships `protein.faa` headed by protein accessions, not locus tags. `protid2locus.tsv`
is the mapping extracted from the GFF; the proteome was relabelled to `ADATZ4_*` headers
before any analysis. **Anything run against NCBI's unmodified protein.faa will not match
locus tags used elsewhere in this project.**
## Files

| File | Query | Result |
|---|---|---|
| `Q_02205_swissprot.tsv` | Q1 | ADATZ4_02205 is α-phosphoglucomutase (bits 554, 100% qcov; next protein in genome 58.9) |
| `Q_02560_swissprot.tsv`, `Q4_align.txt` | Q4 | ADATZ4_02560 is a functional CAT — catalytic His conserved against E. coli CAT-III (P00484) |
| `pfam_domains.txt` | Q2/Q3/Q5/Q6/Q12 | hmmscan 3.4 vs Pfam-A, E≤1e-5, complete 2365/2365 |
| `Q9_eps_cluster_coords.txt` | Q9A | eps locus CDS coordinates and intergenic gaps |
| `Q9_codon_usage.tsv` | Q9A | codon frequencies per 1000 (frame-naive 3-mer walk; bias only, not a publication table) |
| `Q9_up_03440.fa`, `Q9_up_03445.fa` | Q9A | upstream regions scanned for promoter and `cre` site |
| `Q11_nucleases_signalp6.txt` | Q11 | neither surface nuclease is exported |
| `CJM-SD1_genome_virus_summary.tsv` | geNomad | **two complete proviruses, 42.8 + 48.7 kb** |
| `CJM-SD1_genome_plasmid_summary.tsv` | geNomad | both plasmids confirmed; p2 conjugation genes |
| `minced.crisprs`, `minced.gff` | Q8 | **empty by design — no CRISPR arrays detected.** The empty file IS the result |
| `locus_products.tsv` | — | locus_tag / replicon / coords / strand / product, from the PGAP GFF |
| `protid2locus.tsv` | — | NCBI protein accession → ADATZ4 locus tag |

## Headline results

- **Q5**: PFK, FBA (both classes), FBPase (both classes incl. `glpX`) — **all absent** on a
  complete scan. The §4f carbon-entry-topology argument holds.
- **Q2**: complete Leloir operon (GalM-PGM-permease-GalK-GalE-GalT) → **galactose arm added**
  to the Fig 4 panel.
- **Q3**: `mtlD` present (ADATZ4_11390) → **fructose arm dropped**.
- **Q4**: functional CAT → **erythromycin as the primary selection marker**; pNZ8048 demoted.
- **Q6 + Q8**: no restriction-modification system, no CRISPR arrays → no known genomic
  barrier to electroporation.
- **geNomad**: two complete, lysis-competent prophages (holin + GH25 endolysin each);
  the second carries ImmA/IrrE and an AntA/AntB antirepressor, i.e. SOS-inducible.

## Not done in this session

- Q7 (ADATZ4_12090 vs nr, replicon incompatibility) — still running at time of commit
- Proteome-wide SignalP 6 — still running at time of commit
- Q9B (*eps* operon terminator scan) — `09_eps_operon_structure.sh` is not in this repo
- ISEScan — not installed on BioHPC and no substitute found
