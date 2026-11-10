# Bakta

`sample_file` lists sample IDs. Inputs are
`raw_dir/{sample}{input_extension}`. For Shovill assemblies, use
`input_extension: ".contigs.fa"` and point `raw_dir` to the Shovill output
directory. Configure `bakta_db` with the Bakta database path.

Results are written under `out_dir/bakta`, with full per-sample Bakta
directories and parent-level `{sample}.gff3` and `{sample}.fna` links.

```bash
snakemake --snakefile modules/03_gene_annotation/bakta/Snakefile \
  --configfile modules/03_gene_annotation/bakta/config.yaml --use-conda --cores 8
```
