configfile: "Fetch_SRA_taxid/config.yml"

TAXID = config.get("taxid")
DATASET_DIR = config.get("dataset_out_dir")

OUT_PREFIX = f"{DATASET_DIR}/{TAXID}_sra_runexperiments"

METADATA_FILE = f"{OUT_PREFIX}.metadata.tsv"
METADATA_FILTERED_FILE = f"{OUT_PREFIX}_filtered.metadata.tsv"
METADATA_ACCESSIONS_FILE = f"{OUT_PREFIX}_accessions.txt"


rule all:
    input:
        METADATA_FILTERED_FILE,
        METADATA_ACCESSIONS_FILE

rule get_sra_metadata:
    output:
        METADATA_FILE
    params:
        TAXID=TAXID
    conda:
        "ncbi_download"
    shell:
        """
        esearch -db sra -query "txid{params.TAXID}[Organism:exp]" | efetch -format runinfo > {output}
        """


rule get_sra_metadata_filtered:
    input:
        METADATA_FILE
    output:
        METADATA_FILTERED_FILE
    shell:
        """
        awk -F',' '
            NR==1 {{
                for (i = 1; i <= NF; i++) {{
                    gsub(/"/, "", $i)
                    k = tolower($i)
                    gsub(/_/, "", k)
                    col[k] = i
                }}
                print
                next
            }}
            $col["librarystrategy"] == "WGS" &&
            $col["libraryselection"] == "RANDOM" &&
            $col["librarylayout"] == "PAIRED" &&
            $col["librarysource"] == "GENOMIC" &&
            $col["platform"] == "ILLUMINA"
        ' {input[0]} > {output[0]}
        """


rule extract_sra_accessions:
    input:
        METADATA_FILTERED_FILE
    output:
        METADATA_ACCESSIONS_FILE
    shell:
        """
        awk -F',' 'NR>1 {{print $1}}' {input[0]} > {output[0]}
        """
