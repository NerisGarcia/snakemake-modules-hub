conda activate snakemake

snakemake --use-conda --cores 1 --snakefile Fetch_SRA_taxid/snakemake.smk --configfile Fetch_SRA_taxid/config.yml

