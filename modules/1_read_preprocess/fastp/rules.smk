
rule fastp:
    input:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz"
    output:
        r1=f"{CLEAN_DIR}/{{sample}}_1.clean.fastq.gz",
        r2=f"{CLEAN_DIR}/{{sample}}_2.clean.fastq.gz",
        json=f"{CLEAN_DIR}/{{sample}}.fastp.json",
        html=f"{CLEAN_DIR}/{{sample}}.fastp.html"
    threads: config["threads"]["fastp"]
    params:
        out_dir=f"{CLEAN_DIR}"
    conda:
        "genomics_env"
    shell:
        """
        mkdir -p {params.out_dir}
        fastp --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --json {output.json} --html {output.html} \
            --thread {threads} --dedup --trim_tail1 3 --trim_tail2 3 \
            --length_required 120 -q 30 -u 80 \
            --trim_front1 2 --trim_front2 2
        """"


rule fastp_report:
    input:
        json = f'{FASTP_DIR}/{{sample}}_report.json'

    output:
        tsv = f'{FASTP_DIR}/{{sample}}_report.tsv'
    params:
        delimiter = '\t'
    script:
        "scripts/fastp-json2csv.py"

