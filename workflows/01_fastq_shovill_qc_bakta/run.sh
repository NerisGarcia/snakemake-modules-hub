#!/usr/bin/env bash
set -euo pipefail

workflow_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "$workflow_dir/../.." && pwd)"
cores="${1:-8}"

cd "$repo_root"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snakemake

snakemake \
    --use-conda \
    --cores "$cores" \
    --keep-going \
    --snakefile "$workflow_dir/Snakefile" \
    --configfile "$workflow_dir/config.yaml"
