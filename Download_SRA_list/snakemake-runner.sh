conda activate snakemake

snakemake --use-conda --cores 4 --snakefile Download_SRA_list/snakemake.smk --configfile Download_SRA_list/config.yml

