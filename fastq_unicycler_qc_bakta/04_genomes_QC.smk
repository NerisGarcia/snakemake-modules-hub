
# QC ANALYSIS ——————————————————————————————————————————————————————————————

rule qc_analysis:
    input:
        isolates_metadata=manifest_parsed,
        qc_raw_path="data/0_raw_reads_QC/Lcrispatus{dataset}.reads.seqkitstats.txt",
        qc_sub_path="data/0_raw_reads_QC/Lcrispatus{dataset}.reads.sub3M.seqkitstats.txt",
        species_file=species_summary,
        mapping_stats_file="data/2_mapping/2_snippy_cores/Lcrispatus_{dataset}/Lcrispatus_{dataset}.txt",
        quast_stats=quast_summary
    output:
        "data/Lcrispatus{dataset}_NCBI.raw_reads_QC_metadata_combined.ods",
        "data/Lcrispatus{dataset}.metadata.QC_stats_combined.ods",
        "data/Lcrispatus{dataset}.metadata.QC_outlier_criteria_summary.ods"
    conda:
        "envs/r_env.yaml"
    shell:
        "Rscript code/02_QC.R {input.isolates_metadata} {input.qc_raw_path} {input.qc_sub_path} {input.species_file} {input.mapping_stats_file} {input.quast_stats} {output[0]} {output[1]} {output[2]}"


rule qc_plots:
    input:
        mash_dist=INPUTS_CFG["mash_distances_fallback"],
        qc_metadata_final="data/Lcrispatus{dataset}.metadata.QC_stats_combined.ods",
        qc_filters="data/Lcrispatus{dataset}.metadata.QC_outlier_criteria_summary.ods",
        tree_file=lambda wc: f"data/3_pangenomes/{wc.dataset}/core{int(round(float(config['panaroo']['core_threshold']) * 100))}/core_gene_alignment_filtered.treefile",
        species_file=species_summary
    output:
        directory("figures/1_QC_Assembly_mapping/Lcrispatus_{dataset}")
    conda:
        "envs/r_env.yaml"
    shell:
        "Rscript code/02_QC.plots.R {input.mash_dist} {input.qc_metadata_final} {input.qc_filters} {input.tree_file} {input.species_file}"


rule v2_description:
    input:
        qc_stats="data/Lcrispatus{dataset}.metadata.QC_stats_combined.ods"
    output:
        "data/derived/Lcrispatus{dataset}.CultureID.txt",
        directory("figures/2_Lcrispatus_v2/Lcrispatus_{dataset}")
    conda:
        "envs/r_env.yaml"
    shell:
        "Rscript code/02_v2_description.R {input.qc_stats} {output[0]} {output[1]}"


rule run_qc:
    input:
        "data/1_Assembly/stages/{dataset}/genomics.done",
        qc_raw="data/Lcrispatus{dataset}_NCBI.raw_reads_QC_metadata_combined.ods",
        qc_meta="data/Lcrispatus{dataset}.metadata.QC_stats_combined.ods",
        qc_filters="data/Lcrispatus{dataset}.metadata.QC_outlier_criteria_summary.ods",
        qc_figs="figures/1_QC_Assembly_mapping/Lcrispatus_{dataset}",
        filtered_samples="data/derived/Lcrispatus{dataset}.CultureID.txt"
    output:
        touch("data/0_raw_reads_QC/stages/{dataset}/qc.done")

