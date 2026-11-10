# QUAST

Assess each assembly listed in `sample_file` (one sample identifier per line).
The input path is `raw_dir/{sample}{input_extension}`. For example, use
`".fna"` for `raw_dir/{sample}.fna`, or `".contigs.fa"` for flat Shovill
assemblies.
Choose one extension per run; input files need not exist before upstream
workflow rules produce them.

Reports are written to `out_dir/QUAST/{sample}.quast.transposed_report.tsv`.
To run all listed samples independently:

```bash
snakemake --snakefile modules/02_assembly/quast/Snakefile \
  --configfile modules/02_assembly/quast/config.yaml --use-conda --cores 8
```
