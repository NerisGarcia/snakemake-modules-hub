#!!!!!!!!!!! FALTA KRAKEN


# MAPPING STATS ——————————————————————————————————————————————————————————————

rule mapping_stats_dataset:
    input:
        base_mapping=BASE_MAPPING_STATS,
        sample_file=lambda wc: config["datasets"][str(wc.dataset)]["samples"]
    output:
        "data/2_mapping/2_snippy_cores/Lcrispatus_{dataset}/Lcrispatus_{dataset}.txt"
    shell:
        r'''
        mkdir -p $(dirname {output})
        if [[ "{wildcards.dataset}" == "417" ]]; then
            cp {input.base_mapping} {output}
        else
            awk 'NR==1 {print; next} {print}' {input.base_mapping} > {output}.tmp
            awk 'NR==FNR {keep[$1]=1; next} FNR==1 || keep[$1]' {input.sample_file} {output}.tmp > {output}
            rm -f {output}.tmp
        fi
        '''
