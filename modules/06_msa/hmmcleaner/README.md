# HmmCleaner

Clean one alignment per name in `sample_file`. The input path is
`input_dir/{sample}{input_extension}`, with `{sample}` replaced by the alignment
name. All inputs in a run use the same extension.

Cleaned alignments are written to `out_dir/HmmCleaner/{sample}_hmm.fasta`.
HmmCleaner runs on a copy in a per-alignment work directory, so input files
are not modified. Its executable must be available on the command path.

```bash
snakemake --snakefile modules/06_msa/hmmcleaner/Snakefile \
  --configfile modules/06_msa/hmmcleaner/config.yaml --cores 8
```
