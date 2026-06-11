#!/bin/bash

set -Eeo pipefail

# Defaults for the Kube HPC case. Override from the command line if needed:
#   CASE_DIR=tutorials/laserbeamFoam/connor316AMR REMOTE_CASE_NAME=connor316AMR ./upload_hpc_case.sh
CASE_DIR="${CASE_DIR:-tutorials/laserbeamFoam/kubeAL6061_hpc}"
REMOTE_CASE_NAME="${REMOTE_CASE_NAME:-kubeAL6061_keyhole_validations}"
ZIP_NAME="${ZIP_NAME:-$REMOTE_CASE_NAME.zip}"
REMOTE="${REMOTE:-anvil}"
REMOTE_DIR="${REMOTE_DIR:-/anvil/scratch/x-rkanak1}"
ARCHIVE_DIR="${ARCHIVE_DIR:-hpc_uploads}"
UNZIP_REMOTE="${UNZIP_REMOTE:-1}"

repo_dir="$(cd "$(dirname "$0")" && pwd)"
case_path="$repo_dir/$CASE_DIR"
archive_dir="$repo_dir/$ARCHIVE_DIR"
archive_path="$archive_dir/$ZIP_NAME"
staging_root="$archive_dir/.upload_staging"
staged_case="$staging_root/$REMOTE_CASE_NAME"

if [ ! -d "$case_path" ]; then
    echo "ERROR: case folder not found: $case_path" >&2
    exit 1
fi

command -v zip >/dev/null || { echo "ERROR: zip not found" >&2; exit 1; }
command -v scp >/dev/null || { echo "ERROR: scp not found" >&2; exit 1; }

mkdir -p "$archive_dir"
rm -f "$archive_path"
rm -rf "$staging_root"
mkdir -p "$staging_root"
cp -a "$case_path" "$staged_case"

echo "Zipping: $CASE_DIR"
echo "Remote folder name: $REMOTE_CASE_NAME"
echo "Archive: $archive_path"

(
    cd "$staging_root"
    zip -r "$archive_path" "$REMOTE_CASE_NAME" \
        -x '*/processor*/*' \
        -x '*/postProcessing/*' \
        -x '*/post-processing-data/*' \
        -x '*/VTK/*' \
        -x '*/VTKs/*' \
        -x '*/absorptivity_vs_time/*' \
        -x '*/log.*' \
        -x '*/output_*.out' \
        -x '*/output_*.err' \
        -x '*/calibrate_*.out' \
        -x '*/calibrate_*.err' \
        -x '*/meltpool_calibration/runs/*' \
        -x '*/meltpool_calibration/plots/*' \
        -x '*/keyhole_validation/runs/*' \
        -x '*/keyhole_validation/plots/*' \
        -x '*/keyhole_validation/smoke_runs/*' \
        -x '*/keyhole_validation/smoke_plots/*' \
        -x '*/__pycache__/*'
)
rm -rf "$staging_root"

echo "Uploading to: $REMOTE:$REMOTE_DIR/$ZIP_NAME"
scp "$archive_path" "$REMOTE:$REMOTE_DIR/$ZIP_NAME"

remote_case_dir="$REMOTE_DIR/$REMOTE_CASE_NAME"

if [ "$UNZIP_REMOTE" = "1" ]; then
    echo "Unzipping on remote: $remote_case_dir"
    ssh "$REMOTE" "cd '$REMOTE_DIR' && unzip -o '$ZIP_NAME'"
fi

cat <<EOF

Upload complete.

Remote case folder:

$remote_case_dir

Run these on HPC:

ssh $REMOTE
cd $REMOTE_DIR
cd $REMOTE_CASE_NAME

# single simulation:
sbatch job.sh

# calibration sweep:
sbatch calibrate_job.sh

EOF
