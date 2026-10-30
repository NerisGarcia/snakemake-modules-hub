# FASTQ assembly, annotation, and QC with modules

This workflow composes the reusable `fastp`, `shovill`, `bakta`, `busco`, `quast`, and `checkm2` modules. Tool-specific execution and output processing live in their respective module directories; this workflow handles samples and the pipeline summary.

Edit `config.yaml` before running. The workflow maps its common input/output paths and tool settings into each imported module's configuration. The input list contains one SRA accession per line, and paired FASTQ files must be named `{accession}_1.fastq.gz` and `{accession}_2.fastq.gz`.

Run from the repository root, or invoke the supplied runner from any directory:

```bash
bash workflows/01_fastq_shovill_qc_bakta_modules/run.sh
```

The runner activates the `snakemake` Conda environment before starting Snakemake. To run manually:

```bash
conda activate snakemake
snakemake --use-conda --cores 8 --keep-going \
  --snakefile workflows/01_fastq_shovill_qc_bakta_modules/Snakefile \
  --configfile workflows/01_fastq_shovill_qc_bakta_modules/config.yaml
```

The workflow writes cleaned reads, Shovill assemblies, Bakta annotations, BUSCO/QUAST/CheckM2 results, and `data/results/03_qc/pipeline_summary.tsv`.
