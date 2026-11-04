# core_subset

Redefines the core genome of an existing Panaroo pangenome using a chosen set of core-definition samples and restricts the results to an analysis subsample.

## Input
- `pangenome_dir`, `core_threshold`: the `panaroo` module output is read from `{pangenome_dir}/core_{core_threshold}` (needs `alignment: "pan"`).
- `sample_list`, `core_sample_list`, `analysis_sample_list`: one sample ID per line. The last two are optional (`false`/empty = not applied) but at least one is required.

## Output
All under `{pangenome_dir}/{name}/`, where the name is fixed in the Snakefile: `subcore.{dataset_name}.{core_threshold}` (core list only), `core_subsamp.{dataset_name}.{core_threshold}` (analysis list only) or `subcore_subsamp.{dataset_name}.{core_threshold}` (both):
- `gene_presence_absence.selected_samples.csv`: Roary table restricted to the analysis samples.
- `selected_core_groups.tsv`: gene families passing `core_threshold` among the core samples (with status and duplicates).
- `sample_selection.tsv`: effective core and analysis sample sets.
- `selected_gene_fastas/`: one FASTA per selected core gene (input of `gene_alignment`).
- `{name}.pangenome_summary.tsv`: summary of the `gene_alignment` concatenated alignment.
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
