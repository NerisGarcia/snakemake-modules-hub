# BUSCO

List one sample ID per line in `sample_file`. Inputs are selected from
`input_dir` using `{sample}{input_extension}`. Set `busco_lineage` and
`threads` in the config.
The short summary and full table for each sample are written under
`out_dir/BUSCO`.

```bash
snakemake --snakefile modules/02_assembly/busco/Snakefile \
  --configfile modules/02_assembly/busco/config.yaml --use-conda --cores 8
```
