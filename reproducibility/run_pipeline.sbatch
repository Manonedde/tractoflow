#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
# One independent run (own work dir, results, trace) of one pipeline on one
# BIDS dataset, Nextflow local executor inside the allocation.
#
#   sbatch --account=def-desco --job-name=<label> run_pipeline.sbatch \
#          tractoflow|sftractomics <bids_dir> <out_dir> [extra pipeline args...]
# Extra args carry the profile, e.g. "-profile fully_reproducible" (Tractoflow)
# or "-profile full_reproducible,apptainer" (sf-tractomics).
#
# Required environment (exported by submit_test_retest.sh):
#   TRACTOFLOW_DIR, SFTRACTOMICS_DIR, IMG (scilus_2.2.2.sif), APPTAINER_IMG_DIR
set -euo pipefail

PIPELINE=$1; BIDS=$2; OUT=$3; shift 3
HERE=$TRACTOFLOW_DIR/reproducibility

module load apptainer java/17 nextflow
export NXF_OFFLINE=true
export NXF_ANSI_LOG=false
mkdir -p "$OUT" && cd "$OUT"
COMMON=(-w "$OUT/work" -with-trace "$OUT/trace.txt" -with-report "$OUT/report.html" -resume)

case $PIPELINE in
tractoflow)
    # DSL1: Nextflow >=19.04.2, <=21.12.1
    NXF_VER=21.10.6 nextflow run "$TRACTOFLOW_DIR/main.nf" \
        --bids "$BIDS" --output_dir "$OUT/results" \
        --processes "$SLURM_CPUS_PER_TASK" \
        -with-singularity "$IMG" -c "$HERE/rorqual.config" \
        "${COMMON[@]}" "$@"
    ;;
sftractomics)
    export NXF_APPTAINER_CACHEDIR=$APPTAINER_IMG_DIR
    nextflow run "$SFTRACTOMICS_DIR/main.nf" \
        --input "$BIDS" --outdir "$OUT/results" \
        -c "$HERE/rorqual_sftractomics.config" \
        "${COMMON[@]}" "$@"
    ;;
*)
    echo "Unknown pipeline: $PIPELINE"; exit 1 ;;
esac
