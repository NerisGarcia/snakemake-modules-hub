#!/bin/bash
# Convenience script to run the genome characterization workflow
# Usage: ./run.sh [dry-run|run|profile]
#   dry-run: Show what would run without executing
#   run: Execute the workflow (default)
#   profile: Run with job submission profile (requires slurm config)

set -e

WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_HUB_DIR="$(cd "${WORKFLOW_DIR}/../.." && pwd)"

SNAKEFILE="${WORKFLOW_DIR}/Snakefile"
CONFIGFILE="${WORKFLOW_DIR}/config.yaml"
THREADS=${THREADS:-32}
MODE=${1:-run}

cd "${MODULE_HUB_DIR}"

echo "======================================"
echo "Genome Characterization Workflow"
echo "======================================"
echo "Workflow:  ${SNAKEFILE}"
echo "Config:    ${CONFIGFILE}"
echo "Threads:   ${THREADS}"
echo "Mode:      ${MODE}"
echo "Directory: $(pwd)"
echo "======================================"
echo ""

case "${MODE}" in
  dry-run)
    echo "Running in DRY-RUN mode (no execution)..."
    snakemake --use-conda --cores "${THREADS}" --dry-run \
      --snakefile "${SNAKEFILE}" \
      --configfile "${CONFIGFILE}"
    ;;
  
  profile)
    echo "Running with SLURM profile..."
    snakemake --use-conda --cores "${THREADS}" --profile slurm \
      --snakefile "${SNAKEFILE}" \
      --configfile "${CONFIGFILE}"
    ;;
  
  run|*)
    echo "Running workflow..."
    snakemake --use-conda --cores "${THREADS}" \
      --snakefile "${SNAKEFILE}" \
      --configfile "${CONFIGFILE}"
    ;;
esac

echo ""
echo "======================================"
echo "Workflow complete!"
echo "======================================"
