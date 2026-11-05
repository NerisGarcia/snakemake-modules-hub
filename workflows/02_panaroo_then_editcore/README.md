# Subcore / subsample on an existing Panaroo run

This workflow does **not** run Panaroo. It reads a finished Panaroo run (made with `workflows/02_panaroo_tree`, which also gives the original core tree; the folder is `.../{dataset}/panaroo_original.core_{pct}`) and builds custom cores and trees from it. Run it as many times as you like with different lists, labels or thresholds. The folder name `02_panaroo_then_editcore` is kept for now.

It composes modules from this repository only, fetched from GitHub with `github()` (`REPO`/`BRANCH` at the top of the Snakefile), so it can run from any directory, e.g. a server. Modules must be pushed to that branch before running:

1.  `modules/04_pangenome/core_subset` selects single-copy core families and filters the Roary table to the requested samples.
2.  `modules/04_pangenome/gene_aln_cleaning` cleans and realigns the selected families (listed by `core_subset`) and concatenates them.
3.  `modules/05_phylogenetics/iqtree` builds a tree from the custom concatenated alignment.

## Run

From the repository root, after the Panaroo run exists:

``` bash
snakemake --snakefile workflows/02_panaroo_then_editcore/Snakefile \
    --configfile workflows/02_panaroo_then_editcore/config.yaml \
    --cores 8 --use-conda --conda-frontend conda
```

Edit [config.yaml](config.yaml): `panaroo.dir` (full path ending at `panaroo_original.core_{pct}`), the `subcore` / `subsample` lists with their labels, `subcore.core_threshold` and the IQ-TREE options (`iqtree.fast: true` by default, which skips bootstrap; `iqtree.model: "TEST"` runs ModelFinder). The workflow stops with a clear error if the Panaroo files are missing.

## Outputs

Next to the Panaroo folder, under `<parent of panaroo.dir>/<subset name>/`, named `subcore.`, `subsamp.` or `subcore_subsamp.` + `{subcore.label}[__{subsample.label}].core_{pct}` (`pct` = `subcore.core_threshold` × 100; labels are required with their list and describe it, e.g. `isolates_only`, `cladeA1`):

- Custom core group table, sample-selection log, selected gene FASTAs, alignment, AMAS summary, pangenome summary and the IQ-TREE tree.
- `params.txt` records the settings; re-using a label with different settings stops with an error.

At least one of `subcore.core_sample_list` and `subsample.analysis_sample_list` is required.