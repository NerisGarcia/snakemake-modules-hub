# Panaroo followed by custom-core alignment

This workflow composes modules from this repository only, fetched from GitHub with `github()` (`REPO`/`BRANCH` at the top of the Snakefile), so it can run from any directory, e.g. a server. Modules must be pushed to that branch before running:

1.  `modules/04_pangenome/panaroo` builds the pangenome.
2.  `modules/05_phylogenetics/iqtree` builds a tree from Panaroo's filtered core alignment.
3.  `modules/04_pangenome/core_subset` selects single-copy core families and filters the Roary table to the requested samples.
4.  `modules/04_pangenome/gene_alignment` aligns the selected families and concatenates them.
5.  A second instance of `modules/05_phylogenetics/iqtree` builds a tree from the custom concatenated alignment.

The Snakemake DAG connects each module through declared outputs and inputs. `run.sh` preserves the two-stage invocation: stage one creates the Panaroo results and core tree, and stage two consumes them to create the custom core and its tree.

## Run

From the repository root:

``` bash
bash workflows/02_panaroo_then_editcore/run.sh [CORES]
```

Edit [config.yaml](config.yaml) to choose the dataset, sample subsets, core threshold and IQ-TREE options.

## Outputs

- Panaroo files and its summary under `pangenome_dir/core_{core_threshold}/`.
- Custom core group table, sample-selection log, selected gene FASTAs, alignment, AMAS summary, and pangenome summary under `pangenome_dir/<subset name>/`, named `subcore.`, `core_subsamp.` or `subcore_subsamp.` + `{dataset_name}.{core_threshold}` depending on whether `core_sample_list`, `analysis_sample_list` or both are set.
- A core tree and a custom-core tree under their respective IQ-TREE prefixes.

If neither list is set, only the Panaroo stage and core tree run.