#!/usr/bin/env bash
# Run the Panaroo stage, then the subset/alignment/tree stage, from the repository root:
#   bash workflows/02_panaroo_then_editcore/run.sh [CORES]
set -euo pipefail

CORES="${1:-8}"
W="workflows/02_panaroo_then_editcore"

cd "$(git rev-parse --show-toplevel)"

# Stage 1: Panaroo pangenome + core tree
snakemake panaroo_core_all \
    --snakefile "$W/Snakefile" --configfile "$W/config.yaml" \
    --cores "$CORES" --use-conda --conda-frontend conda

# Stage 2: select core genes, align them, summarise the pangenome and build the
# edited-core tree (reuses the Panaroo output from stage 1)
snakemake all \
    --snakefile "$W/Snakefile" --configfile "$W/config.yaml" \
    --cores "$CORES" --use-conda --conda-frontend conda
