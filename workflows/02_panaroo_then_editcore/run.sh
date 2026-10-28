#!/usr/bin/env bash
# Run panaroo core, then editcore, from the repository root:
#   bash workflows/02_panaroo_then_editcore/run.sh [CORES]
set -euo pipefail

CORES="${1:-8}"
W="workflows/02_panaroo_then_editcore"

cd "$(git rev-parse --show-toplevel)"

# Stage 1: panaroo core + tree
snakemake panaroo_core_all \
    --snakefile "$W/Snakefile" --configfile "$W/config.yaml" \
    --cores "$CORES" --use-conda --conda-frontend conda

# Stage 2: editcore (reuses the panaroo output from stage 1)
snakemake all \
    --snakefile "$W/Snakefile" --configfile "$W/config.yaml" \
    --cores "$CORES" --use-conda --conda-frontend conda
