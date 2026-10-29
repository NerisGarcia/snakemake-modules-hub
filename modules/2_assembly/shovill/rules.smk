
rule shovill:
    input:
        r1=f"{CLEAN_DIR}/{{sample}}_1.clean.fastq.gz",
        r2=f"{CLEAN_DIR}/{{sample}}_2.clean.fastq.gz"
    output:
        contigs=f"{ASS_DIR}/{{sample}}_shovill/{{sample}}.contigs.fa"
    params:
        out_dir=f"{ASS_DIR}/{{sample}}_shovill"
    threads: config["threads"]["shovill"]
    conda:
        "genomics_env"
    shell:
        """
        mkdir -p {params.out_dir}
        shovill --R1 {input.r1} --R2 {input.r2} \
            --outdir {params.out_dir} --cpus {threads} --force
        
        mv {params.out_dir}/contigs.fa {output.contigs}

        """