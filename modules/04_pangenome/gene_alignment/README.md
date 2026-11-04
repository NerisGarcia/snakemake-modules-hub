# gene_alignment

Takes the Panaroo core-gene alignments selected by `core_subset`. Genes with at least `min_variable_sites` variable columns have gappy sequences removed (`gap_seq_max`), are realigned with MAFFT and cleaned with HmmCleaner; the remaining genes keep their Panaroo alignment. All genes are padded to the analysis samples and concatenated with AMAS.

## Input
- `pangenome_dir`, `core_threshold`, `core_sample_list`, `analysis_sample_list` (same as `core_subset`; they determine `{name}`, e.g. `subcore_subsamp.{dataset_name}.{core_threshold}`): reads the `core_subset` output in `{pangenome_dir}/{name}/` (`selected_core_groups.tsv`, `selected_gene_fastas/`).
- `analysis_sample_list`: samples (and order) of the final alignment.

## Output
Under `{pangenome_dir}/{name}/`:
- `{name}.gene_alignment.aln`: concatenated alignment.
- `amas/{name}.partitions.txt`, `amas/amas_alignment_summary.tsv` (one row per gene alignment; no per-sequence summaries).
- `core_gene_aln/{name}_{gene}.fasta`: all per-gene alignments used for the concatenation (HmmCleaner output for variable genes, original Panaroo alignment for constant genes), padded to the analysis samples.
- `removed_sequences/{gene}.tsv`: sequences dropped as gappy (`gene`, `sequence`, `gap_fraction`), one file per variable gene (header only if none removed).
- `selected_genes.txt`, `variable_genes.txt` (MAFFT + HmmCleaner path), `constant_genes.txt` (original Panaroo alignment), `logs/`

Filtered, MAFFT and HmmCleaner intermediates are temporary and deleted after use.

## Config
See [config.yaml](config.yaml).

## Environment / databases
- Conda env: `environment.yaml` (mafft, AMAS, HmmCleaner); no databases.

## Usage
```bash
snakemake --use-conda --cores 8 --configfile modules/04_pangenome/gene_alignment/config.yaml
```

## Notes
Gappy sequences are removed with a small script (`scripts/remove_gappy_sequences.py`) instead of trimAl, which cannot drop sequences by gap fraction. A gene with fewer than 2 sequences left fails. HmmCleaner output is assumed to be `<stem>_hmm.fasta`.

`list_selected_genes` and `classify_genes` are checkpoints, as the number of genes is only known after `core_subset` ran.
