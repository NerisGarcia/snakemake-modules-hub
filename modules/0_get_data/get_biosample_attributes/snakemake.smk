

ACCESSIONS_FILE = config.get("biosamp_accessions_file")
DATASET_DIR = config.get("dataset_out_dir")

FILE_PREFIX = config.get("file_prefix")

TMP_DIR = DATASET_DIR + "/" + FILE_PREFIX + "_biosample_attributes_tmp"

import re

with open(ACCESSIONS_FILE) as f:
    BIOSAMPLES = [line.strip() for line in f if re.match(r"^SAM[EDN][AG]?\d+$", line.strip())]

rule all:
    input:
         DATASET_DIR + "/" + FILE_PREFIX + ".biosample_attributes.tsv"

rule get_biosample_xml:
    input:
        ancient(ACCESSIONS_FILE)
    output:
        TMP_DIR + "/{biosample}.xml"
    conda:
        "ncbi_download"
    shell:
        """
        esearch -db biosample -query "{wildcards.biosample}" | \
        efetch -format xml > {output}
        """

rule merge_biosample_attributes:
    input:
        ancient(expand(TMP_DIR + "/{biosample}.xml", biosample=BIOSAMPLES))
    output:
        DATASET_DIR + "/" + FILE_PREFIX + ".biosample_attributes.tsv"
    params:
        tmp_dir = TMP_DIR
    script:
        "parse_biosamples.py"