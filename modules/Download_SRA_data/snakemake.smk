
TAXON = config.get("taxon").replace(" ", "_")
DATASET_DIR = config.get("dataset_out_dir")

OUT_PREFIX = f"{DATASET_DIR}/{TAXON}_sra_runexperiments"

METADATA_FILE = f"{OUT_PREFIX}.metadata.tsv"
METADATA_FILTERED_FILE = f"{OUT_PREFIX}_filtered.metadata.tsv"
METADATA_ACCESSIONS_FILE = f"{OUT_PREFIX}_accessions.txt"
PREFETCHED_ACCESSIONS_FILE = f"{OUT_PREFIX}_prefetched_accessions.txt"
FASTQ_DIR = f"{DATASET_DIR}/SRA_fastq"


DONEFILE = f"{DATASET_DIR}/download.log"


NREAD= config.get("n_reads")

def get_accessions(_wildcards):
    accessions_file = checkpoints.extract_sra_accessions.get().output[0]
    with open(accessions_file) as handle:
        return [line.strip() for line in handle if line.strip()]


def get_prefetched_accessions(_wildcards):
    accessions_file = checkpoints.collect_prefetched_accessions.get().output[0]
    with open(accessions_file) as handle:
        return [line.strip() for line in handle if line.strip()]


def get_prefetch_status_targets(wildcards):
    accessions = get_accessions(wildcards)
    return expand(f"{FASTQ_DIR}/{{acc}}.sra.status", acc=accessions)


def get_fastq_targets(wildcards):
    accessions = get_prefetched_accessions(wildcards)
    return expand(f"{FASTQ_DIR}/{{acc}}_{{read}}.fastq.gz", acc=accessions, read=["1", "2"])


def get_all_targets(wildcards):
    targets = get_fastq_targets(wildcards)
    return [METADATA_FILTERED_FILE, METADATA_ACCESSIONS_FILE] + targets


rule all:
    input:
        get_all_targets

rule get_sra_metadata:
    output:
        METADATA_FILE
    params:
        TAXON=TAXON
    conda:
        "ncbi_download"
    shell:
        """
        esearch -db sra -query "{params.TAXON}[Organism]" | efetch -format runinfo > {output}
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


checkpoint extract_sra_accessions:
    input:
        METADATA_FILTERED_FILE
    output:
        METADATA_ACCESSIONS_FILE
    shell:
        """
        awk -F',' 'NR>1 {{print $1}}' {input[0]} > {output[0]}
        """


rule prefetch_sra:
    output:
        status = f"{FASTQ_DIR}/{{acc}}.sra.status"
    params:
        fastq_dir = f"{FASTQ_DIR}",
        sra = f"{FASTQ_DIR}/{{acc}}/{{acc}}.sra",
        donefile = DONEFILE
    conda:
        "ncbi_download"
    shell:
        """
        mkdir -p {params.fastq_dir}
        if prefetch --output-directory {params.fastq_dir} "{wildcards.acc}" >> {params.donefile} 2>&1 && [ -s {params.sra} ]; then
            echo OK > {output.status}
        else
            echo MISSING > {output.status}
        fi
        """


checkpoint collect_prefetched_accessions:
    input:
        get_prefetch_status_targets
    output:
        PREFETCHED_ACCESSIONS_FILE
    params:
        donefile = DONEFILE
    shell:
        """
        : > {output[0]}
        for status in {input}; do
            if [ "$(cat "$status")" = "OK" ]; then
                accession=$(basename "$status")
                accession=${{accession%.sra.status}}
                echo "$accession" >> {output[0]}
            else
                echo "Warning: Prefetch failed for accession $(basename "$status" .sra.status)" >> {params.donefile}
            fi
        done
        """

rule fastq_dump:
    input:
        sra = f"{FASTQ_DIR}/{{acc}}/{{acc}}.sra"
    output:
        r1_gz = f"{FASTQ_DIR}/{{acc}}_1.fastq.gz",
        r2_gz = f"{FASTQ_DIR}/{{acc}}_2.fastq.gz"
    params:
        fastq_dir = f"{FASTQ_DIR}",
        sra = f"{FASTQ_DIR}/{{acc}}/{{acc}}.sra",
        n_reads = config["n_reads"]
    conda:
        "ncbi_download"
    shell:
        """
        fastq-dump --gzip --skip-technical --split-3 -X {params.n_reads} --outdir {params.fastq_dir} {input.sra}
        """

rule cleanup_sra:
    input:
        sra = f"{FASTQ_DIR}/{{acc}}/{{acc}}.sra",
         status = f"{FASTQ_DIR}/{{acc}}.sra.status"
    shell:
        """
        rm -f {input.sra}
        rm -f {input.status}
        """

