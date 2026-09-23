
snakemake \
    --use-conda \
    --cores 8 \
    --snakefile pangenome_tree/Snakefile \
    --configfile pangenome_tree/config.yaml 
