conda activate snakemake

snakemake --use-conda --cores 4 --snakefile Download_SRA_data/snakemake.smk --configfile Download_SRA_data/config.yml



