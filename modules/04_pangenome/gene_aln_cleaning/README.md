# gene_aln_cleaning

Cleans a set of per-gene alignments and concatenates them. Genes with at least `min_variable_sites` variable columns have gappy sequences removed (`gap_seq_max`), are realigned with MAFFT and cleaned with HmmCleaner; the remaining genes are used as given. All genes are padded to the sample list and concatenated with AMAS. It does not depend on Panaroo: any list of per-gene FASTA alignments works. `core_subset` writes the inputs for the Panaroo case.

## Input
- `alignment_list`: text file with one FASTA alignment path per line. The gene ID is the file name without its last extension (must be unique).
- `sample_file`: samples (and order) of the final alignment; genes lacking a sample are padded with gaps.
- `out_dir`, `prefix`: where results are written and how the concatenated alignment is named.

## Output
Under `out_dir/gene_aln_cleaning/{prefix}/`:
- `{prefix}.gene_alignment.aln`: concatenated alignment.
- `amas/{prefix}.partitions.txt`, `amas/amas_alignment_summary.tsv` (one row per gene alignment; no per-sequence summaries).
- `cleaned_gene_alignments/{gene}.fasta`: per-gene alignments used for the concatenation (HmmCleaner output for variable genes, the input alignment for constant genes), padded to the sample list.
- `removed_sequences/{gene}.tsv`: sequences dropped as gappy (`gene`, `sequence`, `gap_fraction`), one file per variable gene (header only if none removed).
- `gene_alignments.tsv` (gene and path), `variable_genes.txt` (MAFFT + HmmCleaner path), `constant_genes.txt`, `logs/`

Filtered, MAFFT and HmmCleaner intermediates are temporary and deleted after use.

## Config
See [config.yaml](config.yaml).

## Environment / databases
- Conda env: `environment.yaml` (mafft, AMAS, HmmCleaner); no databases.

## Usage
```bash
snakemake --use-conda --cores 8 --configfile modules/04_pangenome/gene_aln_cleaning/config.yaml
```

## Notes
Gappy sequences are removed with a small script (`scripts/remove_gappy_sequences.py`) instead of trimAl, which cannot drop sequences by gap fraction. A gene with fewer than 2 sequences left fails. HmmCleaner output is assumed to be `<stem>_hmm.fasta`.

`read_alignment_list` and `classify_genes` are checkpoints, as the genes are only known once the list exists (it may be produced by another rule, such as `core_subset`).
