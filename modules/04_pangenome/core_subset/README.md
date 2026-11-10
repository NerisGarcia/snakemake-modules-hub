# core_subset

Redefines the core genome of an existing Panaroo pangenome using a chosen set of core-definition samples and restricts the results to an analysis subsample.

## Input
- `pangenome_dir`, `panaroo.run_name`: the `panaroo` module output is read from `{pangenome_dir}/{panaroo.run_name}` (needs `alignment: "pan"`).
- `subcore`: `core_sample_list` (samples that define the core genes), `core_threshold` (0-1, for this subset) and `label` (short description, e.g. `isolates_only`; required with the list).
- `subsample`: `analysis_sample_list` (samples kept in the outputs) and `label` (e.g. `cladeA1`; required with the list).
- `sample_file`: one sample ID per line. The two lists are optional (`false`/empty = not applied) but at least one is required.

## Output
All under `out_dir/core_subset/{name}/`, where the name is fixed in the Snakefile: `subcore.{subcore.label}.core_{pct}` (core list only), `subsamp.{subsample.label}.core_{pct}` (analysis list only) or `subcore_subsamp.{subcore.label}__{subsample.label}.core_{pct}` (both), with `pct` = `subcore.core_threshold` × 100:
- `params.txt`: lists (with checksums), labels, thresholds and Panaroo run used. If the folder already exists with different settings the module stops; use a new label/threshold or delete the folder.
- `gene_presence_absence.selected_samples.csv`: Roary table restricted to the analysis samples.
- `selected_core_groups.tsv`: gene families passing `subcore.core_threshold` among the core samples (with status and duplicates).
- `sample_selection.tsv`: effective core and analysis sample sets.
- `selected_gene_fastas/`: one FASTA per selected core gene.
- `selected_alignments.txt` (paths of those FASTAs) and `analysis_samples.txt` (ordered analysis samples): the inputs of `gene_aln_cleaning`.
- `{name}.pangenome_summary.tsv`: summary of the `gene_aln_cleaning` concatenated alignment, generated when requested by a workflow that provides that alignment.
- `.panaroo_results_validated`, `logs/`

## Config
See [config.yaml](config.yaml).

## Environment / databases
- Python scripts only; no conda env or database.

## Usage
```bash
snakemake --cores 1 --configfile modules/04_pangenome/core_subset/config.yaml
```

## Notes
A sample counts as present in a family only when it has exactly one sequence in it.
