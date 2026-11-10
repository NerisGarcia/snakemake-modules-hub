# Snakemake genomics pipelines

This repository contains Snakemake workflows for downloading NCBI/SRA data, assembling and annotating genomes, running genome QC, and building a reference-free pangenome and tree.

## 1. Obtain SRA data

### From a taxonomic ID

Edit [the SRA taxid module config](modules/00_data_acquisition/fetch_sra_taxid/config.yaml), then run:

``` bash
snakemake --use-conda --cores 1 --snakefile modules/00_data_acquisition/fetch_sra_taxid/Snakefile --configfile modules/00_data_acquisition/fetch_sra_taxid/config.yaml
```

The important settings are `taxid` and `out_dir`; results go under `out_dir/SRA`.

### From an accession list

Edit [the SRA download module config](modules/00_data_acquisition/download_sra_list/config.yaml), then run:

``` bash
snakemake --use-conda --cores 1 --snakefile modules/00_data_acquisition/download_sra_list/Snakefile --configfile modules/00_data_acquisition/download_sra_list/config.yaml
```

The accession file must contain one SRA accession per line. Paired reads are
written under `out_dir/SRA/SRA_fastq` as:

``` text
{accession}_1.fastq.gz
{accession}_2.fastq.gz
```

## 3. Obtain NCBI genome assemblies

### From a taxonomic ID

Edit [the NCBI taxid module config](modules/00_data_acquisition/fetch_ncbi_genomes_taxid/config.yaml), then run:

``` bash
snakemake --use-conda --cores 1 --snakefile modules/00_data_acquisition/fetch_ncbi_genomes_taxid/Snakefile --configfile modules/00_data_acquisition/fetch_ncbi_genomes_taxid/config.yaml
```

This downloads metadata and creates an accession list under `out_dir/NCBI_genomes`. Filter the accession list using the metadata and supply it to the download step.

### Download NCBI genomes from an accession list

Edit [the NCBI download module config](modules/00_data_acquisition/download_ncbi_genomes_list/config.yaml), then run:

``` bash
snakemake --use-conda --cores 1 --snakefile modules/00_data_acquisition/download_ncbi_genomes_list/Snakefile --configfile modules/00_data_acquisition/download_ncbi_genomes_list/config.yaml
```

The accession list contains one NCBI genome accession per line. Results go under `out_dir/NCBI_genomes`.

## X. get Biosample info

``` bash
snakemake --use-conda --cores 1 --snakefile modules/00_data_acquisition/get_biosample_attributes/Snakefile --configfile modules/00_data_acquisition/get_biosample_attributes/config.yaml
```

-- This has to be core 1 if not it gives error

## 2. Run the FASTQ assembly and annotation pipeline

This workflow runs FASTP, Shovill, Bakta, BUSCO, and QUAST.

Edit [the FASTQ workflow config](workflows/01_fastq_shovill_qc_bakta/config.yaml):

``` yaml
sra_acc_list: "path/to/accessions.txt"
sra_raw: "path/to/paired_fastq_files"
clean_dir: "data/clean"
genomics_results_dir: "data/results"
summary_table: "data/results/03_qc/pipeline_summary.tsv"
bakta_db: "/path/to/bakta/database"
busco_lineage: "lactobacillales_odb10"
```

Then run Snakemake directly:

``` bash
snakemake \
  --use-conda \
  --cores 8 \
  --keep-going \
  --snakefile workflows/01_fastq_shovill_qc_bakta/Snakefile \
  --configfile workflows/01_fastq_shovill_qc_bakta/config.yaml
```

`--keep-going` allows independent samples to continue when another sample fails. The configured summary table records which expected files exist for each sample.

The important outputs are:

``` text
data/results/01_assembly/shovill/{sample}.contigs.fa
data/results/02_annotation/bakta/{sample}.gff3
data/results/02_annotation/bakta/{sample}.fna
data/results/03_qc/pipeline_summary.tsv
data/results/03_qc/BUSCO/
data/results/03_qc/QUAST/
```

Bakta writes its complete output inside `{sample}_bakta`. A separate cleanup step creates the parent-level `.gff3` and `.fna` symlinks and removes unwanted `.embl`, `.png`, `.svg`, and hypothetical files. If the nested Bakta outputs already exist, rerunning the workflow can perform only this cleanup step.

### Compose the same pipeline from reusable modules

The module-based workflow is in [workflows/01_fastq_shovill_qc_bakta](workflows/01_fastq_shovill_qc_bakta/README.md). It uses the configured `fastp`, `shovill`, `bakta`, `busco`, `quast`, and `checkm2` modules. Edit that workflow's `config.yaml`, then run:

``` bash
bash workflows/01_fastq_shovill_qc_bakta/run.sh
```

The runner activates the `snakemake` Conda environment before invoking Snakemake. Module directories and filenames follow the convention described in [modules/README.md](modules/README.md).

## 4. Run the FASTA annotation and QC pipeline

This workflow starts from existing genome assemblies and runs Bakta, BUSCO, and QUAST.

Edit [the FASTA workflow config](workflows/01_fasta_qc_bakta/config.yaml):

``` yaml
sample_file: "path/to/genome_accessions.txt"
raw_dir: "path/to/genome_fasta_files"
input_extension: ".fna"
genomics_results_dir: "data/results"
summary_table: "data/results/03_qc/pipeline_summary.tsv"
bakta_db: "/path/to/bakta/database"
busco_lineage: "lactobacillales_odb10"
```

Then run:

``` bash
snakemake \
  --use-conda \
  --cores 8 \
  --keep-going \
  --snakefile workflows/01_fasta_qc_bakta/Snakefile \
  --configfile workflows/01_fasta_qc_bakta/config.yaml
```

Expected assembly names are based on `sample_file`, `raw_dir`, and
`input_extension`. The main outputs are the same annotation and QC layout
described for the FASTQ workflow.

## 5. Run the reference-free pangenome and tree

The pangenome workflow expects one isolate ID per line in [panaroo_tree/config.yaml](panaroo_tree/config.yaml), and a parent-level GFF3 for every isolate:

``` text
{annot_dir}/{isolate}.gff3
```

Edit these settings:

``` yaml
sample_list: "path/to/isolates.txt"
annot_dir: "data/results/02_annotation"
genomics_results_dir: "data/results"
threads: 8
```

Run:

``` bash
snakemake \
  --use-conda \
  --cores 8 \
  --keep-going \
  --snakefile panaroo_tree/Snakefile \
  --configfile panaroo_tree/config.yaml
```

The main outputs are written to:

``` text
data/results/04_pangenome/gene_presence_absence.Rtab
data/results/04_pangenome/gene_presence_absence_roary.csv
data/results/04_pangenome/iqtree/core_gene_alignment_filtered.treefile
```

You can also run a panacota pipeline

snakemake\
--use-conda\
--cores 8\
--snakefile panacota_tree/Snakefile\
--configfile panacota_tree/config.yaml

# Panaroo changing core

## 6. Reuse a previous Panaroo run

`panaroo_editcore_previous` performs the core-gene selection, optional per-gene realignment, concatenation, and tree steps from an existing Panaroo output directory. It does not require GFF files and never runs Panaroo again.

Set `panaroo.existing_pangenome_dir` in [panaroo_editcore_previous/config.yaml](panaroo_editcore_previous/config.yaml) to the directory created by the previous run. The directory must contain `gene_presence_absence.Rtab`, `gene_presence_absence_roary.csv`, `combined_DNA_CDS.fasta`, both core alignment files, and `aligned_gene_sequences/`.

Run:

``` bash
snakemake \
  --use-conda \
  --cores 8 \
  --keep-going \
  --snakefile panaroo_editcore_previous/Snakefile \
  --configfile panaroo_editcore_previous/config.yaml
```

This pipeline creates a pangenome Edit these settings:

``` yaml
sample_list: "path/to/isolates.txt"
annot_dir: "data/results/02_annotation"
genomics_results_dir: "data/results"
threads: 8
```

Run:

``` bash
snakemake \
  --use-conda \
  --cores 8 \
  --keep-going \
  --snakefile panaroo_editcore/Snakefile \
  --configfile panaroo_editcore/config.yaml
```