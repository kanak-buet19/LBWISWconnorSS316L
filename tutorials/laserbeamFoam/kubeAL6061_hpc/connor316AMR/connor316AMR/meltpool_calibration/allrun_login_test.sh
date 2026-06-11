#!/bin/bash

set -Eeo pipefail

fail()
{
    echo "FAIL: $*" >&2
    exit 1
}

check()
{
    echo
    echo "=== $* ==="
}

scriptDir="$(cd "$(dirname "$0")" && pwd)"
baseCase="$(cd "$scriptDir/.." && pwd)"
testCase="$scriptDir/runs/_login_allrun_test"
OF2506_IMAGE="${OF2506_IMAGE:-$HOME/openfoam-dev_2506.sif}"
OF2506_USER="${OF2506_USER:-x-rkanak1}"
LOGIN_TEST_TIMEOUT="${LOGIN_TEST_TIMEOUT:-180}"
LOGIN_TEST_SUBDOMAINS="${LOGIN_TEST_SUBDOMAINS:-8}"
LOGIN_TEST_END_TIME="${LOGIN_TEST_END_TIME:-0.008}"
LOGIN_TEST_WRITE_INTERVAL="${LOGIN_TEST_WRITE_INTERVAL:-0.00025}"
PYTHON="${PYTHON:-/home/x-rkanak1/.conda/envs/isw_env/bin/python3}"

of2506()
{
    apptainer exec --cleanenv --env USER="$OF2506_USER" "$OF2506_IMAGE" bash -lc "source /openfoam/bash.rc && $*"
}

patch_scalar()
{
    local file="$1"
    local key="$2"
    local value="$3"
    perl -0pi -e "s/(^\\s*\\Q$key\\E\\s+)[^;]+(;)/\${1}$value\${2}/m" "$file"
}

check "Login Allrun Smoke Test"
echo "base case: $baseCase"
echo "test case: $testCase"
echo "timeout: ${LOGIN_TEST_TIMEOUT}s"
echo "subdomains: $LOGIN_TEST_SUBDOMAINS"
echo "endTime: $LOGIN_TEST_END_TIME"
echo "writeInterval: $LOGIN_TEST_WRITE_INTERVAL"
echo "OF2506_IMAGE: $OF2506_IMAGE"
echo "OF2506_USER: $OF2506_USER"

command -v apptainer >/dev/null || fail "apptainer not found"
[ -f "$OF2506_IMAGE" ] || fail "missing image: $OF2506_IMAGE"
[ -x "$baseCase/Allrun" ] || fail "base Allrun not executable"
[ -x "$baseCase/Allrun_long" ] || fail "base Allrun_long not executable"

check "Create isolated mini case"
rm -rf "$testCase"
mkdir -p "$scriptDir/runs"
rsync -a \
    --exclude '0' \
    --exclude 'processor*' \
    --exclude 'VTK' \
    --exclude 'post-processing-data' \
    --exclude 'meltpool_calibration' \
    --exclude 'log.*' \
    --exclude '*.foam' \
    "$baseCase/" "$testCase/"

patch_scalar "$testCase/system/controlDict" "endTime" "$LOGIN_TEST_END_TIME"
patch_scalar "$testCase/system/controlDict" "writeInterval" "$LOGIN_TEST_WRITE_INTERVAL"
patch_scalar "$testCase/system/decomposeParDict" "numberOfSubdomains" "$LOGIN_TEST_SUBDOMAINS"

check "Environment"
of2506 "echo WM_PROJECT_VERSION=\$WM_PROJECT_VERSION; echo FOAM_USER_APPBIN=\$FOAM_USER_APPBIN; command -v laserbeamFoam; command -v setFields; command -v blockMesh; command -v decomposePar"

check "Run ./Allrun in test case"
set +e
of2506 "cd '$testCase' && export PYTHON='$PYTHON' && export MPLCONFIGDIR=/tmp && timeout '$LOGIN_TEST_TIMEOUT' ./Allrun"
status=$?
set -e

echo "Allrun exit status: $status"
if [ "$status" -eq 124 ]; then
    echo "WARN: timeout hit. Treating as smoke-test pass only if solver started."
elif [ "$status" -ne 0 ]; then
    echo "Recent logs:"
    find "$testCase" -maxdepth 1 -name 'log.*' -print -exec tail -40 {} \;
    fail "Allrun failed"
fi

check "Log checks"
[ -f "$testCase/log.blockMesh" ] || fail "missing log.blockMesh"
[ -f "$testCase/log.setFields" ] || fail "missing log.setFields"
[ -f "$testCase/log.decomposePar" ] || fail "missing log.decomposePar"
[ -f "$testCase/log.laserbeamFoam" ] || fail "missing log.laserbeamFoam"

grep -q "End" "$testCase/log.blockMesh" || fail "blockMesh did not complete"
grep -q "End" "$testCase/log.setFields" || fail "setFields did not complete"
grep -q "End" "$testCase/log.decomposePar" || fail "decomposePar did not complete"
grep -q "Starting time loop\|Time =" "$testCase/log.laserbeamFoam" || fail "laserbeamFoam did not start time loop"

check "Result"
echo "PASS: ./Allrun starts correctly inside OF2506 Apptainer on login node."
echo "Test case kept at: $testCase"
echo "Remove later with: rm -rf '$testCase'"
