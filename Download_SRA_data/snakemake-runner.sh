conda activate snakemake

snakemake --use-conda --cores 4 --snakefile code/Module1_download_SRA_data/snakemake.smk --rerun-incomplete



datasets summary genome taxon lactobacillus --as-json-lines | dataformat tsv genome
