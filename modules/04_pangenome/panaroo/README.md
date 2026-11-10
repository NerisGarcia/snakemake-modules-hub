# panaroo

Builds a pangenome from Bakta GFF3 files with Panaroo and writes a pangenome summary table.

## Input
- `sample_file`: one isolate ID per line.
- `input_dir` and `input_extension` identify each annotation. For example,
  `".gff3"` selects `input_dir/{sample}.gff3`.

## Output
All under `out_dir/panaroo/{run_name}/`, where `run_name` defaults to `panaroo_original.core_{pct}` (`pct` = `core_threshold` × 100, e.g. 0.9 → `core_90`); set `run_name` to override:
- `gene_presence_absence.Rtab`, `gene_presence_absence_roary.csv`
- `core_gene_alignment.aln`, `core_gene_alignment_filtered.aln` and `.embl` headers
- `combined_DNA_CDS.fasta`, `aligned_gene_sequences/` (used by `core_subset`)
- `{run_name}.pangenome_summary.tsv`
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
