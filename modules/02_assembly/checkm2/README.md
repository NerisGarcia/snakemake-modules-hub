# CheckM2

List one sample ID per line in `sample_file`. Inputs are selected from
`raw_dir` using `{sample}{input_extension}`. Set `checkm2_db` to the
CheckM2 database file.
Each report is written to `out_dir/CheckM2/{sample}/quality_report.tsv`.

```bash
snakemake --snakefile modules/02_assembly/checkm2/Snakefile \
  --configfile modules/02_assembly/checkm2/config.yaml --use-conda --cores 8
```
