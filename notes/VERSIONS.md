# Software versions (Cornell BioHPC, cbsumm13, Rocky Linux 9.8)
| Step | Tool | Version | Path |
|---|---|---|---|
| Read QC | seqkit | 2.13.0 | /programs/seqkit-2.13.0 |
| Read QC plots | R / ggplot2 | 4.4.2 / 4.0.2 | /programs/R-4.4.2 |
| Read filtering | Filtlong | 0.2.1 | /programs/Filtlong-0.2.1 |

Note: BioHPC NanoPlot 1.38.1 fails on Rocky 9 (Python 3.6 build, numpy error),
so read QC plots are made in R instead.
| Assembly | Flye | 2.9.2-b1795 | /programs/Flye-2.9.2 |

Reads: R10.4.1 flow cell, Dorado basecalling model dna_r10.4.1_e8.2_400bps_sup@v4.3.0 (from FASTQ header)
| Polishing | dorado polish --bacteria | 1.1.1 (model dna_r10.4.1_e8.2_400bps_polish_bacterial_methylation_v5.0.0) | /programs/dorado-1.1.1-linux-x64 |
| Alignment/BAM | samtools | 1.20 | /programs/samtools-1.20 |
| Taxonomy (ANI) | skani | 0.3.1 | /programs/skani-0.3.1 |
| Reference download | NCBI datasets | 18.35.0 | /programs/ncbi_datasets_18.35.0 |
| Completeness/contamination | CheckM (lineage_wf) | 1.2.4, data checkm_data_2015_01_16 | /programs/checkm-1.2.4 |
| CheckM dependency | HMMER | 3.4 | /programs/hmmer-3.4 |
| CheckM dependency | prodigal | 2.6.3 | /programs/prodigal-2.6.3 |
| CheckM dependency | pplacer | v1.1.alpha19 | /programs/pplacer-Linux-v1.1.alpha19 |
| Annotation | Bakta (container bakta.sif) | 1.9.3, DB v5.1 light (doi:10.5281/zenodo.10522951) | /programs/bakta-1.9.3 |
| Domain annotation | HMMER hmmscan vs Pfam-A | 3.4 / Pfam-A Dec 2017 | /programs/hmmer-3.4, /programs/PfamScan/db |
| AMR screen | AMRFinderPlus | 3.12.8, DB 2024-07-22.1 | /programs/amrfinder-3.12.8 |
