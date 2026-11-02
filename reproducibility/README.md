# Reproducibility and variability study: Tractoflow (scilpy 2.2.2) vs sf-tractomics

Dataset: Myelo-inferno (BIDS, 20 subjects × 5-6 sessions), on Alliance Rorqual:
`/home/medde/projects/rrg-descotea/datasets/CHUS/myelo_inferno/`. Account: `def-desco`.

Both pipelines run the same tools from one image, `scilus/scilus:2.2.2`
(FSL 6.0.7.18, ANTs 2.6.2, MRtrix 3.0.7, scilpy 2.2.2), with equivalent
reproducibility settings:

| | Tractoflow | sf-tractomics |
|---|---|---|
| code | branch `scilpy_2.2.2` (port of 2.4.2) | `-profile full_reproducible` |
| threads | `-profile fully_reproducible` | same per step (`conf/full_reproducible.config`) |
| image | `-with-singularity scilus_2.2.2.sif` | profile forces `scilus/scilus:2.2.2` on every Tractoflow-block step |
| steps | `--run_gibbs_correction false --run_t1_denoising false --run_local_tracking true` | `--run_gibbs_correction false --run_t1_denoising false` |

Same steps in both: no Gibbs correction, no T1 denoising (nlmeans), PFT and
local tracking. Threads in both: dwidenoise 4, fODF 4, everything else 1. ANTs seed 1234, `OPENBLAS_NUM_THREADS=1`, tracking seeds 0. The T1
template (`mni_152_sym_09c/t1`) is identical in both (same md5).

## Port: behaviour changes to be aware of

| Step | Change |
|---|---|
| all b0 handling | `--force_b0_threshold` → `--skip_b0_check`, `--b0_thr` → `--b0_threshold` |
| Denoise_T1 | `scil_denoising_nlmeans --basic_sigma --number_coils 1` (same as sf-tractomics) |
| Crop_DWI / Crop_T1 | bounding boxes in `.json` instead of `.pkl` |
| Compute_FRF | `--fa/--min_fa` → `--fa_thresh/--min_fa_thresh` |
| FODF_Metrics | `--fa_t/--md_t` → `--fa_threshold/--md_threshold` |
| Local_Tracking | `--sf` → `--sfthres` |
| SH basis | `descoteaux07` is the **non-legacy** basis in scilpy 2.x (legacy in 1.x). Kept as `descoteaux07`, like sf-tractomics. `--basis descoteaux07_legacy` reproduces Tractoflow 2.4.2 coefficients. |

## 1. Setup (login node, once)

```bash
git clone -b scilpy_2.2.2 git@github.com:Manonedde/tractoflow.git
git clone -b scilpy_2.2.2 git@github.com:Manonedde/sf-tractomics.git
export SFTRACTOMICS_DIR=$PWD/sf-tractomics IMG_DIR=$PROJECT/containers
export SRC_BIDS=/home/medde/projects/rrg-descotea/datasets/CHUS/myelo_inferno
cd tractoflow/reproducibility
./setup_login_node.sh $IMG_DIR $SRC_BIDS
```

Builds `scilus_2.2.2.sif`, caches Nextflow 21.10.6 (Tractoflow is DSL1: needs
Nextflow ≤ 21.12.1 and Java 17), and pulls every sf-tractomics image with a
stub run on one subject (compute nodes have no internet).

## 2. Check the BIDS input (login node, no computation)

```bash
# Reverse b0: add IntendedFor (Tractoflow needs it). Report first, then write a
# mirror (symlinks to the original files + edited JSONs, source untouched).
python3 fix_intendedfor.py $SRC_BIDS
python3 fix_intendedfor.py $SRC_BIDS --out_dir $SCRATCH/repro_study/bids   # [--sbref_as_epi]
export BIDS=$SCRATCH/repro_study/bids

./check_bids.sh $BIDS $SCRATCH/check_bids
```

`fix_intendedfor.py` gives each DWI one reverse b0 with the opposite
`PhaseEncodingDirection` (`fmap/*_epi`, else `fmap/*_sbref`, else
`dwi/*_sbref`; same run first) and adds the DWI to its `IntendedFor`. With a
single DWI, Tractoflow rejects a reverse **sbref** (only `fmap/*_epi`):
`--sbref_as_epi` exposes it as `fmap/*_dir-XX_epi` in the mirror. Full
reverse DWIs (`dir-PA_dwi` with bval/bvec) need nothing.

One row per session: DWI, reverse b0/DWI and T1 selected by Tractoflow
(`scil_bids_validate`), and whether topup runs in each pipeline (sf-tractomics
through a stub run). Any `MISMATCH` must be fixed before submitting. Known
difference: Tractoflow only uses a reverse b0 in `fmap/` with `IntendedFor`
pointing to the DWI, or a reverse DWI (`dir-PA_dwi`); sf-tractomics also
matches `sbref`/`epi` by direction.

## 3. Runs (login node)

```bash
export ACCOUNT=def-desco      # IMG_DIR, SFTRACTOMICS_DIR, BIDS as above
OUT=$SCRATCH/repro_study

./submit_test_retest.sh runretest $BIDS $OUT 3   # computational: 3 subjects, runs A and B
./submit_test_retest.sh sessions  $BIDS $OUT     # all subjects/sessions, 1 job per subject
```

- `runretest`: 8 jobs (2 pipelines × {reproducible, default threading} × A/B).
  Every step is compared from the Nextflow trace and `work/` (keep `work/`).
- `sessions`: 40 jobs (20 subjects × 2 pipelines), reproducible profiles.
- Every job uses `-resume`: resubmit the same command after a timeout
  (72 h, 32 CPUs, 120 GB per job; eddy and tracking are single-threaded).

## 4. Analysis

```bash
sbatch --account=$ACCOUNT --job-name=analyze analyze.sh $OUT $IMG_DIR/scilus_2.2.2.sif
```

**Report: `$OUT/analysis/report/index.html`** (figures + tables; blue
Tractoflow, orange sf-tractomics; subjects coded S01…, mapping in
`qc/subject_codes.csv`):

1. run-retest: % bitwise-identical files per step, size of the differences;
2. variability: CV intra vs inter, ICC, per-subject boxplots (sessions as
   dots) for each step's measures, histograms of Dice / r with the first
   session, FA trajectories across sessions;
3. Tractoflow vs sf-tractomics: Dice / r per step, Bland-Altman of WM means,
   and one QC image per session in `qc/sessions/` (both pipelines side by
   side + difference map, mask contours overlaid).

Tables in `$OUT/analysis/`:

- `runretest_<pipeline>_<repro|mt>.csv`: every file of every step (from trace.txt + work/), A vs B:
  bitwise identical?, max |diff|, relative mean |diff|, r, Dice (masks),
  streamline counts (tractograms).
- `variability/sessions.csv`: per pipeline × subject × session, native space:
  WM means of DTI/fODF maps, mask volumes, FRF, streamline counts.
- `variability/agreement.csv` / `agreement_summary.csv`: each session rigidly
  registered (FA) on the subject's first session: Dice of masks and
  tractogram density maps, r and relative |diff| of maps in WM.
- `variability/group.csv`: per pipeline × measure: CV intra-subject,
  CV inter-subject, their ratio, ICC(1) (one-way, unbalanced: 5 or 6
  sessions per subject).
- `qc/pipelines.csv`: per session, Tractoflow vs sf-tractomics (Dice of
  masks and tractogram densities, r and WM means of maps).




## Command run on rorqual
git clone -b scilpy_2.2.2 git@github.com:Manonedde/tractoflow.git
git clone -b scilpy_2.2.2 git@github.com:Manonedde/sf-tractomics.git
export SFTRACTOMICS_DIR=$PWD/sf-tractomics IMG_DIR=$PROJECT/containers
export SRC_BIDS=/home/medde/projects/rrg-descotea/datasets/CHUS/myelo_inferno
cd tractoflow/reproducibility
./setup_login_node.sh $IMG_DIR $SRC_BIDS
python3 fix_intendedfor.py $SRC_BIDS                                    # lire le rapport
python3 fix_intendedfor.py $SRC_BIDS --out_dir $SCRATCH/repro_study/bids
export BIDS=$SCRATCH/repro_study/bids
./check_bids.sh $BIDS $SCRATCH/check_bids                               # 0 mismatch attendu
./submit_test_retest.sh runretest $BIDS $SCRATCH/repro_study 3