# Subcore / subsample on an existing Panaroo run

This workflow does **not** run Panaroo. It reads a finished Panaroo run (made with `workflows/02_panaroo_tree`, which also gives the original core tree; the folder is `.../{dataset}/panaroo_original.core_{pct}`) and builds custom cores and trees from it. Run it as many times as you like with different lists, labels or thresholds. The folder name `02_panaroo_then_editcore` is kept for now.

It loads modules from the local repository copy using paths relative to the workflow Snakefile. The same repository layout must be present on the server; changes do not need to be pushed to GitHub before running. Paths for data and outputs in the configuration are relative to the working directory unless absolute paths are provided:

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

Edit [config.yaml](config.yaml), setting:

- `panaroo.rtab`: the existing `gene_presence_absence.Rtab`, used in the pangenome summary.
- `panaroo.roary`: the existing `gene_presence_absence_roary.csv`, used to select and filter gene families.
- `panaroo.gene_sequences`: the existing `aligned_gene_sequences/` directory, created by Panaroo with `alignment: "pan"`.
- `out_dir`: where to write the new results.
- The `subcore` / `subsample` lists and their labels, `subcore.core_threshold`, and the IQ-TREE options (`iqtree.fast: true` skips bootstrap; `iqtree.model: "TEST"` runs ModelFinder).

The full sample list and its order come from the Rtab. Each sample must appear exactly once in the Roary CSV, and both selection lists must be subsets of the Rtab samples. Missing inputs or mismatched samples stop the workflow before rules run. The three Panaroo inputs must come from the same run; the original combined CDS and core alignment files are not needed.

## Outputs

Under `out_dir/core_subset/<subset name>/`, named `subcore.`, `subsamp.` or `subcore_subsamp.` + `{subcore.label}[__{subsample.label}].core_{pct}` (`pct` = `subcore.core_threshold` × 100; labels are required with their list and describe it, e.g. `isolates_only`, `cladeA1`):

- Custom core group table, sample-selection log, selected gene FASTAs, alignment, AMAS summary, pangenome summary and the IQ-TREE tree.
- `params.txt` records the settings; re-using a label with different settings stops with an error.
- Cleaned alignments go under `out_dir/gene_aln_cleaning/<subset name>/` and the tree under `out_dir/iqtree/`.

At least one of `subcore.core_sample_list` and `subsample.analysis_sample_list` is required. An omitted list uses all Panaroo samples in Rtab order. Existing `params.txt` files from the previous configuration format will not match; use a new label or rebuild the derived output folder.