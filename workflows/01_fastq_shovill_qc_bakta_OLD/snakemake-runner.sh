
snakemake \
    --use-conda \
    --cores 8 \
    --snakefile fastq_shovill_qc_bakta/Snakefile \
    --configfile fastq_shovill_qc_bakta/config.yaml 
