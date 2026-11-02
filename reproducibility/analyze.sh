#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
# Compare runs once all jobs are done. Submit from tractoflow/reproducibility/:
#   sbatch --account=def-desco --job-name=analyze analyze.sh <out_root> <image.sif>
set -euo pipefail
OUT=$1; IMG=$2
HERE=$SLURM_SUBMIT_DIR
module load apptainer
export MPLCONFIGDIR=$SLURM_TMPDIR
RUN="apptainer exec --no-home -B /project -B /scratch -B /home -B $HERE $IMG python"
A=$OUT/analysis; mkdir -p "$A"

# 1. Computational run-retest: every file of every step, run A vs run B
for p in tractoflow sftractomics; do
    for cfg in repro mt; do
        R=$OUT/runretest/$p/${cfg}
        [ -f "${R}_A/trace.txt" ] || continue
        $RUN "$HERE/compare_runs.py" "${R}_A" "${R}_B" \
            "$A/runretest_${p}_${cfg}.csv" > "$A/runretest_${p}_${cfg}.txt"
    done
done

# 2. Intra/inter-subject variability: gather the per-subject jobs per pipeline
ARGS=()
for p in tractoflow sftractomics; do
    M=$OUT/sessions/$p/merged; mkdir -p "$M"
    for s in "$OUT"/sessions/$p/sub-*/results/*/; do
        n=$(basename "$s")
        if [[ $n == sub-* ]]; then ln -sfn "$s" "$M/$n"; fi
    done
    layout=$p
    ARGS+=(--run "$p" "$layout" "$M")
done
$RUN "$HERE/variability.py" "${ARGS[@]}" --out_dir "$A/variability" \
    --nproc "$SLURM_CPUS_PER_TASK"

# 3. Tractoflow vs sf-tractomics: QC images per session + agreement
$RUN "$HERE/report.py" qc --tractoflow "$OUT/sessions/tractoflow/merged" \
    --sftractomics "$OUT/sessions/sftractomics/merged" --out_dir "$A/qc" \
    --nproc "$SLURM_CPUS_PER_TASK"

# 4. Figures + HTML report: $A/report/index.html
$RUN "$HERE/report.py" plots --analysis_dir "$A"
