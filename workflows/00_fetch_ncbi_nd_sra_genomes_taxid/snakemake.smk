

module download_genomes:
    snakefile:
        "../../modules/00_data_acquisition/fetch_ncbi_genomes_taxid/Snakefile"
    config:
        config

module download_sra:
    snakefile:
        "../../modules/00_data_acquisition/fetch_sra_taxid/Snakefile"
    config:
        config

rule all:
    input:
        f"{config['out_dir']}/NCBI_genomes/{config['taxid']}_ncbi_genomes_summary.tsv",
        f"{config['out_dir']}/NCBI_genomes/{config['taxid']}_ncbi_genomes_accessions.txt",
        f"{config['out_dir']}/SRA/{config['taxid']}_sra_runexperiments_filtered.metadata.tsv",
        f"{config['out_dir']}/SRA/{config['taxid']}_sra_runexperiments_accessions.txt"

use rule * from download_genomes as ncbi_*
use rule * from download_sra as sra_*