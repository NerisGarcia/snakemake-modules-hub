# Kraken2 and Bracken

Provide `sample_file`, `input_dir`, `pairend_extension_pattern` (with `{read}`
for reads 1 and 2), and `input_extension`. `kraken_db` selects the database;
configure the Bracken read length and classification level in the config.
Reports are written under `out_dir/kraken2`.

```bash
snakemake --snakefile modules/01_read_preprocessing/kraken2/Snakefile \
  --configfile modules/01_read_preprocessing/kraken2/config.yaml --use-conda --cores 8
```
