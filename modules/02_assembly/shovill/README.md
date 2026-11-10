# Shovill

Provide `sample_file`, `raw_dir`, `pairend_extension_pattern` (with `{read}`
for reads 1 and 2), and `input_extension`. Each assembly is written as
`out_dir/shovill/{sample}_shovill/{sample}.contigs.fa`. A flat relative
symlink, `out_dir/shovill/{sample}.contigs.fa`, lets downstream modules use
`input_extension: ".contigs.fa"`. Existing nested assemblies are reused to
create the link without rerunning Shovill.

```bash
snakemake --snakefile modules/02_assembly/shovill/Snakefile \
  --configfile modules/02_assembly/shovill/config.yaml --use-conda --cores 8
```
