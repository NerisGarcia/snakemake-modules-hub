# AMAS summary

Summarize a set of DNA alignments with AMAS. The sample file lists alignment
names, and each input is read from `input_dir/{sample}{input_extension}`.

## Input
- `sample_file`: one alignment name per line; blank lines and lines beginning
  with `#` are ignored.
- `input_dir`: directory containing the alignments.
- `input_extension`: extension appended to each name in `sample_file`.

## Output
- `out_dir/amas_summary/amas_alignment_summary.tsv`: AMAS alignment summary.
- `out_dir/amas_summary/logs/amas_summary.log`: command log.

## Config
See [config.yaml](config.yaml).

## Environment
- Conda environment: [amas.yaml](amas.yaml)

## Usage
```bash
snakemake --snakefile modules/06_msa/amas_summary/Snakefile \
  --configfile modules/06_msa/amas_summary/config.yaml --use-conda --cores 1
```
