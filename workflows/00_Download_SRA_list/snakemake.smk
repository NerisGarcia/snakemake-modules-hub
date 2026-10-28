
configfile: "Download_SRA_list/config.yml"

ACCESSIONS_FILE = config.get("accessions_file")
DATASET_DIR = config.get("dataset_out_dir")

FASTQ_DIR = f"{DATASET_DIR}/SRA_fastq"
DONEFILE = f"{DATASET_DIR}/download.log"
PREFETCHED_ACCESSIONS_FILE = f"{DATASET_DIR}/prefetched_accessions.txt"

NREAD= config.get("n_reads")

def get_accessions(_wildcards):
    with open(ACCESSIONS_FILE) as handle:
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


rule all:
    input:
        get_fastq_targets

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
