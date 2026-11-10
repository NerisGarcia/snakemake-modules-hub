# Cleanifier

Remove human reads from paired FASTQ inputs. `sample_file` contains one sample
name per line. Each paired input is `raw_dir/{sample}` followed by
`pairend_extension_pattern` and `input_extension`; `{read}` is replaced by 1
or 2. Optionally, set `unpaired_input_pattern` to a path relative to `raw_dir`
containing `{sample}`; the same `input_extension` is appended. Use null when
there is no unpaired input.

Outputs are compressed and written under `out_dir/cleanifier` as
`{sample}.cleanifier.R1.fastq.gz`, `{sample}.cleanifier.R2.fastq.gz`, and,
when configured, `{sample}.cleanifier.unpaired.fastq.gz`. This module can run
before fastp on paired reads or after fastp with its combined unpaired reads.
Per-sample seqkit statistics are written to `out_dir/cleanifier/qc`.

```bash
snakemake --snakefile modules/01_read_preprocessing/cleanifier/Snakefile \
  --configfile modules/01_read_preprocessing/cleanifier/config.yaml --use-conda --cores 8
```
