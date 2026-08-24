
# ASSEMBLY, ANNOTATION AND QC ——————————————————————————————————————————————————————————————

rule assembly_isolate:
    input:
        r1=lambda wc: f"data/0_raw_reads_QC/2_subsampled/{wc.isolate}.sub3M.R1.fq.gz",
        r2=lambda wc: f"data/0_raw_reads_QC/2_subsampled/{wc.isolate}.sub3M.R2.fq.gz"
    output:
        "data/1_Assembly/1_Assembly/{isolate}_shovill/{isolate}.contigs.fa"
    conda:
        "envs/shovill_env.yaml"
    shell:
        r'''
        mkdir -p $(dirname {output})
        shovill --R1 {input.r1} --R2 {input.r2} --outdir data/1_Assembly/1_Assembly/{wildcards.isolate}_shovill --cpus 16 --ram 32
        '''

rule assembly:
    input:
        expand("data/1_Assembly/1_Assembly/{isolate}_shovill/{isolate}.contigs.fa", isolate=ALL_ISOLATES)

rule annotation_isolate:
    input:
        "data/1_Assembly/1_Assembly/{isolate}_shovill/{isolate}.contigs.fa"
    output:
        "data/1_Assembly/3_Annotation/{isolate}_bakta/{isolate}.gff3"
    conda:
        "envs/baktav1.11.4.yaml"
    shell:
        r'''
        mkdir -p $(dirname {output})
        bakta --db {config[bakta][database]} -o data/1_Assembly/3_Annotation/{wildcards.isolate}_bakta/ --prefix {wildcards.isolate}  --locus-tag {wildcards.isolate} --threads 16 {input} --keep-contig-headers --force
        '''

rule annotation:
    input:
        expand("data/1_Assembly/3_Annotation/{isolate}_bakta/{isolate}.gff3", isolate=ALL_ISOLATES)

rule assembly_qc_isolate:
    input:
        "data/1_Assembly/1_Assembly/{isolate}_shovill/{isolate}.contigs.fa"
    output:
        "data/1_Assembly/2_QC/1_Quast/{isolate}.quast.transposed_report.tsv"
    conda:
        "envs/quast_env.yaml"
    shell:
        r'''
        mkdir -p $(dirname {output})
        quast.py {input} -o data/1_Assembly/2_QC/1_Quast/{wildcards.isolate}_quast --threads 16 --glimmer --no-sv --no-icarus --no-plots --min-contig 200

        mv data/1_Assembly/2_QC/1_Quast/{wildcards.isolate}_quast/transposed_report.tsv {output}
        rm -rf data/1_Assembly/2_QC/1_Quast/{wildcards.isolate}_quast
        '''

rule assembly_qc:
    input:
        expand("data/1_Assembly/2_QC/1_Quast/{isolate}.quast.transposed_report.tsv", isolate=ALL_ISOLATES)


rule assembly_annotation_qc:
    input:
        rules.assembly.input,
        rules.annotation.input,
        rules.assembly_qc.input
