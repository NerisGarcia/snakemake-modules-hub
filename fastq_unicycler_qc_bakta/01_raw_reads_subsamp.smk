# RAW READS QC ——————————————————————————————————————————————————————————————

rule raw_reads_stats:
    input:
        lambda wc: expand(
            "resources/Isolate_seqs/reads/0_raw_reads/{isolate}.R1.fq.gz",
            isolate=[
                line.strip().strip('"')
                for line in open(config["datasets"][str(wc.dataset)]["samples"])
                if line.strip()
                and line.strip().lower() not in {"cultureid", "culture_id", "sampleid", "sample_id"}
            ],
        )
    output:
        "data/0_raw_reads_QC/Lcrispatus{dataset}.reads.seqkitstats.txt"
    conda:
        "envs/basics_env.yaml"
    shell:
        r'''
        mkdir -p $(dirname {output})
        rm -f {output}
        for f in {input}; do
            seqkit stats "$f" >> {output}
        done
        '''

# SUBSAMPLED READS QC ——————————————————————————————————————————————————————————————
rule subsample_read:
    input:
        R1="resources/Isolate_seqs/reads/0_raw_reads/{isolate}.R1.fq.gz",
        R2="resources/Isolate_seqs/reads/0_raw_reads/{isolate}.R2.fq.gz"
    output:
        R1_sub="data/0_raw_reads_QC/2_subsampled/{isolate}.sub3M.R1.fq.gz",
        R2_sub="data/0_raw_reads_QC/2_subsampled/{isolate}.sub3M.R2.fq.gz"
    conda:
        "envs/basics_env.yaml"
    shell:
        r'''
        mkdir -p data/0_raw_reads_QC/2_subsampled
        seqtk sample -s100 {input.R1} 1500000 | gzip > {output.R1_sub}
        seqtk sample -s100 {input.R2} 1500000 | gzip > {output.R2_sub}
        '''


rule subsampled_reads_stats:
    input:
        lambda wc: expand(
            [
                "data/0_raw_reads_QC/2_subsampled/{isolate}.sub3M.R1.fq.gz",
                "data/0_raw_reads_QC/2_subsampled/{isolate}.sub3M.R2.fq.gz",
            ],
            isolate=[
                line.strip().strip('"')
                for line in open(config["datasets"][str(wc.dataset)]["samples"])
                if line.strip()
                and line.strip().lower() not in {"cultureid", "culture_id", "sampleid", "sample_id"}
            ],
        )
    output:
        "data/0_raw_reads_QC/Lcrispatus{dataset}.reads.sub3M.seqkitstats.txt"
    conda:
        "envs/basics_env.yaml"
    shell:
        r'''
        mkdir -p $(dirname {output})
        rm -f {output}
        for f in {input}; do
            seqkit stats "$f" >> {output}
        done
        '''


#!!!!!!!!!!! FALTA FASTP
