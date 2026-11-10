# fastp

`sample_file` lists sample IDs. Paired inputs are found beneath `input_dir` using `pairend_extension_pattern` (`{read}` is 1 or 2) and `input_extension`. Filtering flags are configured in `fastp.extra_args`; there are no additional fixed filtering flags.

`out_dir/fastp` contains `{sample}.fastp.R1{input_extension}` and `{sample}.fastp.R2{input_extension}`, a combined `{sample}.fastp.unpaired{input_extension}`, and JSON, HTML, and TSV reports. The intermediate unpaired R1/R2 files are managed as temporary outputs.

``` bash
snakemake --snakefile modules/01_read_preprocessing/fastp/Snakefile \
  --configfile modules/01_read_preprocessing/fastp/config.yaml --use-conda --cores 8
```