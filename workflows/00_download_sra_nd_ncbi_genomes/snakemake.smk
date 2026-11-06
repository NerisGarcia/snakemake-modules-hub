
module download_genomes:
    snakefile:
        "../../modules/00_data_acquisition/download_ncbi_genomes_list/Snakefile"
    config:
        {**config, "accessions_file": config["ncbi_accessions_file"]}

module download_sra:
    snakefile:
        "../../modules/00_data_acquisition/download_sra_list/Snakefile"
    config:
        {**config, "accessions_file": config["sra_accessions_file"]}

use rule * from download_genomes as ncbi_*
use rule * from download_sra as sra_*

rule all:
    default_target: True
    input:
        rules.ncbi_all.input,
        rules.sra_all.input
