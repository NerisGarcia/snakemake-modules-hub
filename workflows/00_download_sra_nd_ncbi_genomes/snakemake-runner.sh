conda activate snakemake

snakemake --use-conda --cores 1 --snakefile code/0_download_genomes/snakemake.smk --configfile code/0_download_genomes/config.yml --conda-frontend conda --keep-going





