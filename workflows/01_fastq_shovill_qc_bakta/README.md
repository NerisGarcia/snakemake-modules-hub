# FASTQ assembly, annotation, and QC with modules

This workflow composes the reusable `fastp`, `kraken2`, `shovill`, `bakta`, `busco`, `quast`, and `checkm2` modules. Tool-specific execution and output processing live in their respective module directories; this workflow handles samples and the pipeline summary.

Edit `config.yaml` before running. The workflow maps its common input/output paths and tool settings into each imported module's configuration. The input list contains one SRA accession per line, and paired FASTQ files must be named `{accession}_1.fastq.gz` and `{accession}_2.fastq.gz`.

Run from the repository root, or invoke the supplied runner from any directory:

```bash
bash workflows/01_fastq_shovill_qc_bakta/run.sh
```

The runner activates the `snakemake` Conda environment before starting Snakemake. To run manually:

```bash
conda activate snakemake
snakemake --use-conda --cores 8 --keep-going \
  --snakefile workflows/01_fastq_shovill_qc_bakta/Snakefile \
  --configfile workflows/01_fastq_shovill_qc_bakta/config.yaml
```

Samples without both non-empty raw FASTQ files are skipped, marked `FAIL` in the input columns and in the downstream summary columns, and do not prevent other samples from running. The workflow writes cleaned reads, Kraken2 and Bracken reports (on the fastp-cleaned reads), Shovill assemblies, Bakta annotations, BUSCO/QUAST/CheckM2 results, and `data/results/03_qc/pipeline_summary.tsv`.
