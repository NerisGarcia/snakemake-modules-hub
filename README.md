# Snakemake genomics pipelines

This repository contains Snakemake workflows for downloading NCBI/SRA data, assembling and annotating genomes, running genome QC, and building a reference-free pangenome and tree.


## 1. Obtain SRA data

### From a taxonomic ID

Edit [Fetch_SRA_taxid/config.yml](Fetch_SRA_taxid/config.yml), then run:

```bash
snakemake --use-conda --cores 1 --snakefile Fetch_SRA_taxid/snakemake.smk --configfile Fetch_SRA_taxid/config.yml
```

The important settings are `taxid` and `dataset_out_dir`.

### From an accession list

Edit [Download_SRA_list/config.yml](Download_SRA_list/config.yml), then run:

```bash
bash Download_SRA_list/snakemake-runner.sh
```

The accession file must contain one SRA accession per line. The resulting paired reads must be available as:

```text
{accession}_1.fastq.gz
{accession}_2.fastq.gz
```

## 2. Run the FASTQ assembly and annotation pipeline

This workflow runs FASTP, Shovill, Bakta, BUSCO, and QUAST.

Edit [fastq_shovill_qc_bakta/config.yaml](fastq_shovill_qc_bakta/config.yaml):

```yaml
sra_acc_list: "path/to/accessions.txt"
sra_raw: "path/to/paired_fastq_files"
clean_dir: "data/clean"
genomics_results_dir: "data/results"
summary_table: "data/results/03_qc/pipeline_summary.tsv"
bakta_db: "/path/to/bakta/database"
busco_lineage: "lactobacillales_odb10"
```

Then run Snakemake directly:

```bash
snakemake \
  --use-conda \
  --cores 8 \
  --keep-going \
  --snakefile fastq_shovill_qc_bakta/Snakefile \
  --configfile fastq_shovill_qc_bakta/config.yaml
```

`--keep-going` allows independent samples to continue when another sample fails. The configured summary table records which expected files exist for each sample.

The important outputs are:

```text
data/results/01_assembly/{sample}_shovill/{sample}.contigs.fa
data/results/02_annotation/{sample}.gff3
data/results/02_annotation/{sample}.fna
data/results/03_qc/pipeline_summary.tsv
data/results/03_qc/BUSCO/
data/results/03_qc/QUAST/
```

Bakta writes its complete output inside `{sample}_bakta`. A separate cleanup step creates the parent-level `.gff3` and `.fna` symlinks and removes unwanted `.embl`, `.png`, `.svg`, and hypothetical files. If the nested Bakta outputs already exist, rerunning the workflow can perform only this cleanup step.

## 3. Obtain NCBI genome assemblies

### From a taxonomic ID

Edit [Fetch_ncbi_genomes_taxid/config.yml](Fetch_ncbi_genomes_taxid/config.yml), then run:

```bash
snakemake --use-conda --cores 1 --snakefile Fetch_ncbi_genomes_taxid/snakemake.smk
```
This dowloads all metadata and creates an accession list. You cna filter the accession list using the metadata and feed the filtered accession list to the dowload step

### Dowload ncbi genomes from an accession list

Edit [Download_ncbi_genomes_list/config.yml](Download_ncbi_genomes_list/config.yml), then run:

```bash
snakemake --use-conda --cores 1 --snakefile Download_ncbi_genomes_list/snakemake.smk
```

The accession list contains one NCBI genome accession per line. Set the output directory in the config before running.

## 4. Run the FASTA annotation and QC pipeline

This workflow starts from existing genome assemblies and runs Bakta, BUSCO, and QUAST.

Edit [fasta_qc_bakta/config.yaml](fasta_qc_bakta/config.yaml):

```yaml
genome_acc_list: "path/to/genome_accessions.txt"
ass_dir: "path/to/genome_fasta_files"
genomics_results_dir: "data/results"
summary_table: "data/results/03_qc/pipeline_summary.tsv"
bakta_db: "/path/to/bakta/database"
busco_lineage: "lactobacillales_odb10"
```

Then run:

```bash
snakemake \
  --use-conda \
  --cores 8 \
  --keep-going \
  --snakefile fasta_qc_bakta/Snakefile \
  --configfile fasta_qc_bakta/config.yaml
```

Expected assembly names are based on the sample IDs in `genome_acc_list`. The main outputs are the same annotation and QC layout described for the FASTQ workflow.

## 5. Run the reference-free pangenome and tree

The pangenome workflow expects one isolate ID per line in [panaroo_tree/config.yaml](panaroo_tree/config.yaml), and a parent-level GFF3 for every isolate:

```text
{annot_dir}/{isolate}.gff3
```

Edit these settings:

```yaml
sample_list: "path/to/isolates.txt"
annot_dir: "data/results/02_annotation"
genomics_results_dir: "data/results"
threads: 8
```

Run:

```bash
snakemake \
  --use-conda \
  --cores 8 \
  --keep-going \
  --snakefile panaroo_tree/Snakefile \
  --configfile panaroo_tree/config.yaml
```

The main outputs are written to:

```text
data/results/04_pangenome/gene_presence_absence.Rtab
data/results/04_pangenome/gene_presence_absence_roary.csv
data/results/04_pangenome/iqtree/core_gene_alignment_filtered.treefile
```


You can also run a panacota pipeline


snakemake \
    --use-conda \
    --cores 8 \
    --snakefile panacota_tree/Snakefile \
    --configfile panacota_tree/config.yaml 
