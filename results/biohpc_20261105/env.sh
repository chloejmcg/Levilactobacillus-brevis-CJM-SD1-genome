# BioHPC environment for CJM-SD1 (GCA_062071506.1) — source before any step
export PATH=/programs/ncbi-blast-2.16.0+/bin:$PATH
export PATH=/programs/hmmer-3.4/bin:$PATH
export PATH=/programs/seqkit-2.13.0:$PATH
export PATH=/programs/EMBOSS-6.6.0/bin:$PATH
export PATH=/programs/ncbi_datasets_18.35.0:$PATH
export PATH=/programs/edirect:$PATH

export WORK=/workdir/cjm423
export ACC=GCA_062071506.1
export DATA=$WORK/ncbi_raw/ncbi_dataset/data/$ACC

export PFAM_DB=/programs/PfamScan/db/Pfam-A.hmm      # Dec 2017, 16712 families
export NR_DB=/shared_data/genome_db/BLAST_NCBI/nr
export SP_DB=$WORK/db/swissprot                       # shared copy is broken; local pull
