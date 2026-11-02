#!/bin/bash
# To run ONCE on a Rorqual login node (internet access required).
#   ./setup_login_node.sh <image_dir> <bids_dir>
# - builds scilus_2.2.2.sif (Tractoflow port)
# - caches Nextflow 21.10.6 (Tractoflow, DSL1)
# - pulls every image sf-tractomics needs, through a stub run (no computation)
set -euo pipefail

IMG_DIR=${1:?"Usage: $0 <image_dir> <bids_dir>"}; BIDS=${2:?"Usage: $0 <image_dir> <bids_dir>"}
HERE=$(cd "$(dirname "$0")" && pwd)
SFTRACTOMICS_DIR=${SFTRACTOMICS_DIR:?"export SFTRACTOMICS_DIR=/path/to/sf-tractomics"}
mkdir -p "$IMG_DIR"

module load apptainer java/17 nextflow
export APPTAINER_CACHEDIR=${SCRATCH}/.apptainer_cache
export NXF_APPTAINER_CACHEDIR=$IMG_DIR

[ -f "$IMG_DIR/scilus_2.2.2.sif" ] || \
    apptainer build "$IMG_DIR/scilus_2.2.2.sif" docker://scilus/scilus:2.2.2
NXF_VER=21.10.6 nextflow -version

# One subject is enough for the stub run to pull the images
STUB=$SCRATCH/sftractomics_stub; rm -rf "$STUB"; mkdir -p "$STUB/bids"
cp "$BIDS/dataset_description.json" "$STUB/bids/"
ln -s "$(find "$BIDS" -mindepth 1 -maxdepth 1 -type d -name 'sub-*' | sort | head -1)" "$STUB/bids/"
(cd "$STUB" && nextflow run "$SFTRACTOMICS_DIR/main.nf" -stub \
    -profile full_reproducible,apptainer --input "$STUB/bids" --outdir "$STUB/out")
ls -lh "$IMG_DIR"

# Tool versions, for the methods section
apptainer exec "$IMG_DIR/scilus_2.2.2.sif" bash -c '
echo "FSL $(cat $FSLDIR/etc/fslversion)"; mrconvert --version | head -1
antsRegistration --version | head -1; scil_get_version'
