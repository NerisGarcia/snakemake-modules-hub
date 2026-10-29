rule bakta:
    input:
        contigs=assembly_input
    output:
        gff3=f"{ANN_DIR}/{{sample}}_bakta/{{sample}}.gff3",
        fna=f"{ANN_DIR}/{{sample}}_bakta/{{sample}}.fna"
    params:
        out_dir=f"{ANN_DIR}/{{sample}}_bakta",
        prefix="{sample}",
        bakta_db=BAKTA_DB_DIR
    threads: config["threads"]["bakta"]
    conda:
        "genomics_env"
    shell:
        """
        mkdir -p {params.out_dir}
        bakta --db {params.bakta_db} --output {params.out_dir} \
            --prefix {params.prefix} --locus-tag {params.prefix} \
            --threads {threads} {input.contigs} --keep-contig-headers --force
        """


rule bakta_cleanup:
    input:
        gff3=f"{ANN_DIR}/{{sample}}_bakta/{{sample}}.gff3",
        fna=f"{ANN_DIR}/{{sample}}_bakta/{{sample}}.fna"
    output:
        gff3=f"{ANN_DIR}/{{sample}}.gff3",
        fna=f"{ANN_DIR}/{{sample}}.fna"
    params:
        out_dir=f"{ANN_DIR}/{{sample}}_bakta",
        prefix="{sample}"
    shell:
        """
        bakta_dir=$(cd "{params.out_dir}" && pwd)
        gff3_path=$(find "$bakta_dir" -type f -name "{params.prefix}.gff3" -print -quit)
        fna_path=$(find "$bakta_dir" -type f -name "{params.prefix}.fna" -print -quit)
        test -n "$gff3_path" && test -n "$fna_path"
        ln -sfn "$gff3_path" "{output.gff3}"
        ln -sfn "$fna_path" "{output.fna}"

        rm -f {params.out_dir}/{params.prefix}.embl \
              {params.out_dir}/{params.prefix}.png \
              {params.out_dir}/{params.prefix}.svg \
              {params.out_dir}/{params.prefix}.hypotheticals.faa \
              {params.out_dir}/{params.prefix}.hypotheticals.tsv
        """
