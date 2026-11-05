# iqtree

Builds a maximum-likelihood tree from a core gene alignment: extracts SNP sites with snp-sites (plus a SNP distance matrix with snp-dists) and runs IQ-TREE with ascertainment-bias correction (`-fconst`).

## Input
- `alignment`: FASTA alignment (e.g. from `panaroo` or `gene_aln_cleaning`).

## Output
- `{prefix}.treefile`, `{prefix}.iqtree`: tree and IQ-TREE report.
- `{prefix}.vcf`, `{prefix}.fasta.DM`, `{prefix}.SNPs.fasta`, `{prefix}.fconst.txt`
- `{prefix}_logs/`

## Config
See [config.yaml](config.yaml): `alignment`, `prefix`, `model` (a fixed model such as `GTR+I+G`, or `TEST` to run ModelFinder and use the best model; the chosen model is in `{prefix}.iqtree`), `bootstrap` (ultrafast bootstrap replicates; 0 = none), `fast` (true adds IQ-TREE `--fast` and skips bootstrap, because IQ-TREE cannot combine them), `threads`. Instantiate the module several times with different `alignment`/`prefix` to build several trees.

## Environment / databases
- Conda env: `environment.yaml`; no databases.

## Usage
```bash
snakemake --use-conda --cores 8 --configfile modules/05_phylogenetics/iqtree/config.yaml
```
