# panaroo

Builds a pangenome from Bakta GFF3 files with Panaroo and writes a pangenome summary table.

## Input
- `sample_list`: one isolate ID per line.
- `annot_dir/{isolate}.gff3`: annotations.

## Output
All under `pangenome_dir/core_{core_threshold}/`:
- `gene_presence_absence.Rtab`, `gene_presence_absence_roary.csv`
- `core_gene_alignment.aln`, `core_gene_alignment_filtered.aln` and `.embl` headers
- `combined_DNA_CDS.fasta`, `aligned_gene_sequences/` (used by `core_subset` and `gene_alignment`)
- `{dataset_name}.core{core_threshold}.pangenome_summary.tsv`
- `logs/`

The edited GFFs (`;gene=` removed) are temporary and deleted after Panaroo finishes. Panaroo output files keep their original names.

## Config
See [config.yaml](config.yaml).

## Environment / databases
- Conda env: `environment.yaml`; no databases.

## Usage
```bash
snakemake --use-conda --cores 8 --configfile modules/04_pangenome/panaroo/config.yaml
```
