

module download_genomes:
    snakefile:
        github("NerisGarcia/snakemake-modules-hub", path="modules/00_data_acquisition/fetch_ncbi_genomes_taxid/Snakefile", branch="main")
    config:
        config

module download_sra:
    snakefile:
        github("NerisGarcia/snakemake-modules-hub", path="modules/00_data_acquisition/fetch_sra_taxid/Snakefile", branch="main")
    config:
        config

rule all:
    input:
        f"{config['dataset_out_dir']}/{config['taxid']}_ncbi_genomes_summary.tsv",
        f"{config['dataset_out_dir']}/{config['taxid']}_ncbi_genomes_accessions.txt",
        f"{config['dataset_out_dir']}/{config['taxid']}_sra_runexperiments_filtered.metadata.tsv",
        f"{config['dataset_out_dir']}/{config['taxid']}_sra_runexperiments_accessions.txt"

use rule * from download_genomes as ncbi_*
use rule * from download_sra as sra_*