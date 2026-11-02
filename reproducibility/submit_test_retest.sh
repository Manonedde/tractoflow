#!/bin/bash
# Submit the Tractoflow (scilpy 2.2.2) / sf-tractomics reproducibility study on
# Rorqual, from a login node:
#
#   ./submit_test_retest.sh runretest <bids_dir> <out_root> [n_subjects=3]
#   ./submit_test_retest.sh sessions  <bids_dir> <out_root>
#
# runretest : computational reproducibility. The first n subjects (all their
#             sessions) are processed twice (A/B) by each pipeline, with the
#             reproducible profiles ("repro") and with default threading ("mt").
#             Every step is compared from the Nextflow trace + work/ (kept).
# sessions  : intra/inter-subject variability. Every subject, all sessions,
#             one job per subject and pipeline, reproducible profiles.
#
# Required environment:
#   IMG_DIR (from setup_login_node.sh), SFTRACTOMICS_DIR; ACCOUNT (default def-desco)
set -euo pipefail

MODE=$1; BIDS=$2; OUT=$3
ACCOUNT=${ACCOUNT:-def-desco}
export TRACTOFLOW_DIR=$(cd "$(dirname "$0")/.." && pwd)
export SFTRACTOMICS_DIR=${SFTRACTOMICS_DIR:?}
export IMG=$IMG_DIR/scilus_2.2.2.sif APPTAINER_IMG_DIR=$IMG_DIR
cd "$TRACTOFLOW_DIR/reproducibility"

# Same steps in both pipelines: no Gibbs correction, no T1 denoising (nlmeans),
# local tracking on (off by default in Tractoflow).
TF_ARGS=(--run_gibbs_correction false --run_t1_denoising false --run_local_tracking true)
SF_ARGS=(--run_gibbs_correction false --run_t1_denoising false)
TF_REPRO=(-profile fully_reproducible)
SF_REPRO=(-profile full_reproducible,apptainer)
SF_MT=(-profile apptainer)

subject_bids() {  # <sub-dir> -> BIDS dataset with only this subject
    local d=$OUT/inputs/$(basename "$1")
    mkdir -p "$d"; cp "$BIDS/dataset_description.json" "$d/"
    ln -sfn "$1" "$d/"
    [ -f "$BIDS/participants.tsv" ] && cp "$BIDS/participants.tsv" "$d/"
    echo "$d"
}

submit() {  # <job_name> <pipeline> <bids> <out_dir> [args...]
    local name=$1; shift
    sbatch --account="$ACCOUNT" --job-name="$name" --export=ALL run_pipeline.sbatch "$@"
}

mapfile -t SUBS < <(find "$BIDS" -mindepth 1 -maxdepth 1 -type d -name 'sub-*' | sort)
echo "${#SUBS[@]} subjects"

case $MODE in
runretest)
    N=${4:-3}
    IN=$OUT/inputs/runretest; mkdir -p "$IN"; cp "$BIDS/dataset_description.json" "$IN/"
    for s in "${SUBS[@]:0:N}"; do ln -sfn "$s" "$IN/"; done
    for run in A B; do
        submit tf_rr_repro_$run tractoflow "$IN" "$OUT/runretest/tractoflow/repro_$run" "${TF_REPRO[@]}" "${TF_ARGS[@]}"
        submit tf_rr_mt_$run tractoflow "$IN" "$OUT/runretest/tractoflow/mt_$run" "${TF_ARGS[@]}"
        submit sf_rr_repro_$run sftractomics "$IN" "$OUT/runretest/sftractomics/repro_$run" "${SF_REPRO[@]}" "${SF_ARGS[@]}"
        submit sf_rr_mt_$run sftractomics "$IN" "$OUT/runretest/sftractomics/mt_$run" "${SF_MT[@]}" "${SF_ARGS[@]}"
    done
    ;;
sessions)
    for s in "${SUBS[@]}"; do
        sub=$(basename "$s"); IN=$(subject_bids "$s")
        submit tf_$sub tractoflow "$IN" "$OUT/sessions/tractoflow/$sub" "${TF_REPRO[@]}" "${TF_ARGS[@]}"
        submit sf_$sub sftractomics "$IN" "$OUT/sessions/sftractomics/$sub" "${SF_REPRO[@]}" "${SF_ARGS[@]}"
    done
    ;;
*)
    echo "Unknown mode: $MODE"; exit 1 ;;
esac
