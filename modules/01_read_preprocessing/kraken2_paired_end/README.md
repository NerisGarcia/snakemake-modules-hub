# Kraken2 paired and unpaired reads

Provide `sample_file`, `input_dir`, `pairend_extension_pattern` (with `{read}`
for reads 1 and 2), `unpaired_input_pattern` (with `{sample}`), and a shared
`input_extension`. This variant requires both paired and unpaired inputs.
It combines their Kraken reports and runs Bracken, writing results under
`out_dir/kraken2`.

```bash
snakemake --snakefile modules/01_read_preprocessing/kraken2_paired_end/Snakefile \
  --configfile modules/01_read_preprocessing/kraken2_paired_end/config.yaml --use-conda --cores 8
```
