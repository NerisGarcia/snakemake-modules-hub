
snakemake \
-n \
    --use-conda \
    --cores 8 \
    --snakefile panaroo_tree/Snakefile \
    --configfile panaroo_tree/config.yaml 
