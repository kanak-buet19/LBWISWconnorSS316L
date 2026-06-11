#!/bin/bash

set -Eeo pipefail

# Upload a selected Tao Sun HPC case.
# Override from the command line if needed:
#   CASE_DIR=refined_mesh_convergence REMOTE=anvil REMOTE_DIR=/anvil/scratch/x-rkanak1 ./upload_hpc.sh
REMOTE="${REMOTE:-anvil}"
REMOTE_DIR="${REMOTE_DIR:-/anvil/scratch/x-rkanak1}"
ARCHIVE_DIR="${ARCHIVE_DIR:-hpc_uploads}"
UNZIP_REMOTE="${UNZIP_REMOTE:-1}"

case_root="$(cd "$(dirname "$0")" && pwd)"

case_options=(
    "mesh_convergence:TaoSun_mesh_convergence"
    "refined_mesh_convergence:TaoSun_refined_mesh_convergence"
    "reference_case:TaoSun_reference_case"
    "case_original_thermal_prop:TaoSun_original_thermal_prop"
    "case_zhengtao_thermal:TaoSun_zhengtao_thermal"
)

selectCase()
{
    local choice=""
    local i=""
    local option=""

    echo "Choose case to zip and upload:"
    for i in "${!case_options[@]}"
    do
        option="${case_options[$i]}"
        printf '  %d) %s\n' "$((i + 1))" "${option%%:*}"
    done

    while true
    do
        read -r -p "Selection [1-${#case_options[@]}]: " choice
        if [[ "$choice" =~ ^[0-9]+$ ]] \
            && [ "$choice" -ge 1 ] \
            && [ "$choice" -le "${#case_options[@]}" ]; then
            option="${case_options[$((choice - 1))]}"
            CASE_DIR="${option%%:*}"
            REMOTE_CASE_NAME="${option#*:}"
            return
        fi

        echo "Invalid selection: $choice" >&2
    done
}

if [ -z "${CASE_DIR:-}" ]; then
    selectCase
else
    REMOTE_CASE_NAME="${REMOTE_CASE_NAME:-TaoSun_${CASE_DIR}}"
fi

ZIP_NAME="${ZIP_NAME:-$REMOTE_CASE_NAME.zip}"
case_path="$case_root/$CASE_DIR"
archive_dir="$case_root/$ARCHIVE_DIR"
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
        -x '*/.mesh_convergence_complete' \
        -x '*/__pycache__/*'
)
rm -rf "$staging_root"

echo "Uploading to: $REMOTE:$REMOTE_DIR/$ZIP_NAME"
scp "$archive_path" "$REMOTE:$REMOTE_DIR/$ZIP_NAME"

remote_case_dir="$REMOTE_DIR/$REMOTE_CASE_NAME"

if [ "$UNZIP_REMOTE" = "1" ]; then
    echo "Unzipping on remote: $remote_case_dir"
    ssh "$REMOTE" "cd '$REMOTE_DIR' && unzip -o '$ZIP_NAME' && rm -f '$ZIP_NAME'"
fi

cat <<EOF

Upload complete.

Remote case folder:

$remote_case_dir

Run these on HPC:

ssh $REMOTE
cd $remote_case_dir
sbatch job.sh

EOF
