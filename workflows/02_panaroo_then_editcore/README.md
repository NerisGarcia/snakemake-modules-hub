# Panaroo followed by custom-core alignment

This workflow composes modules from this repository only, fetched from GitHub with `github()` (`REPO`/`BRANCH` at the top of the Snakefile), so it can run from any directory, e.g. a server. Modules must be pushed to that branch before running:

1.  `modules/04_pangenome/panaroo` builds the pangenome.
2.  `modules/05_phylogenetics/iqtree` builds a tree from Panaroo's filtered core alignment.
3.  `modules/04_pangenome/core_subset` selects single-copy core families and filters the Roary table to the requested samples.
4.  `modules/04_pangenome/gene_aln_cleaning` cleans and realigns the selected families (listed by `core_subset`) and concatenates them.
5.  A second instance of `modules/05_phylogenetics/iqtree` builds a tree from the custom concatenated alignment.

The Snakemake DAG connects each module through declared outputs and inputs. `run.sh` preserves the two-stage invocation: stage one creates the Panaroo results and core tree, and stage two consumes them to create the custom core and its tree.

## Run

From the repository root:

``` bash
bash workflows/02_panaroo_then_editcore/run.sh [CORES]
```

Edit [config.yaml](config.yaml) to choose the dataset, sample subsets, the two core thresholds (`panaroo.core_threshold` and `subcore.core_threshold`, which defaults to the Panaroo one) and IQ-TREE options (`iqtree.fast: true` by default, which skips bootstrap; `iqtree.model: "TEST"` runs ModelFinder).

## Outputs

- Panaroo files and its summary under `pangenome_dir/panaroo_original.core_{pct}/` (`pct` = `panaroo.core_threshold` × 100; override with `panaroo.run_name`). It is reused, not re-run, when you add other subsets.
- Custom core group table, sample-selection log, selected gene FASTAs, alignment, AMAS summary, and pangenome summary under `pangenome_dir/<subset name>/`, named `subcore.`, `subsamp.` or `subcore_subsamp.` + `{subcore.label}[__{subsample.label}].core_{pct}` (labels are required with their list and describe it, e.g. `isolates_only`, `cladeA1`) depending on whether `subcore.core_sample_list`, `subsample.analysis_sample_list` or both are set.
- `params.txt` in each subset folder records its settings; re-using a label with different settings stops with an error.
- A core tree and a custom-core tree under their respective IQ-TREE prefixes.

If neither `subcore.core_sample_list` nor `subsample.analysis_sample_list` is set, only the Panaroo stage and core tree run.