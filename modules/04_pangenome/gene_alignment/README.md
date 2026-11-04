# gene_alignment

Aligns each selected core gene with MAFFT, pads/finalizes it to the analysis samples, and concatenates the genes into one core alignment with AMAS.

## Input
- `pangenome_dir`, `core_threshold`, `core_sample_list`, `analysis_sample_list` (same as `core_subset`; they determine `{name}`, e.g. `subcore_subsamp.{dataset_name}.{core_threshold}`): reads the `core_subset` output in `{pangenome_dir}/{name}/` (`selected_core_groups.tsv`, `selected_gene_fastas/`).
- `analysis_sample_list`: samples (and order) of the final alignment.

## Output
Under `{pangenome_dir}/{name}/`:
- `{name}.gene_alignment.aln`: concatenated alignment.
- `amas/{name}.partitions.txt`, `amas/amas_alignment_summary.tsv` (one row per gene alignment; no per-sequence summaries).
- `gene_alignments/`: finalized per-gene alignments.
- `selected_genes.txt`, `logs/`

Raw MAFFT alignments are temporary and removed after finalization.

## Config
See [config.yaml](config.yaml).

## Environment / databases
- Conda env: `environment.yaml` (mafft, AMAS); no databases.

## Usage
```bash
snakemake --use-conda --cores 8 --configfile modules/04_pangenome/gene_alignment/config.yaml
```

## Notes
`list_selected_genes` is a checkpoint, as the number of genes is only known after `core_subset` ran.
