#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

runner="$repo_root/scripts/run_three_city_reproducibility.sh"
checker="$repo_root/scripts/check_three_city_results.py"

echo "=========================================="
echo "THREE-CITY REPRODUCIBILITY VERIFICATION"
echo "=========================================="

if [ ! -x "$runner" ]; then
    echo "ERROR: Simulation runner not executable"
    exit 1
fi

if [ ! -f "$checker" ]; then
    echo "ERROR: Regression checker not found"
    exit 1
fi

logfile="$(mktemp /tmp/uhi-three-city-run.XXXXXX)"

if ! "$runner" > "$logfile" 2>&1; then
    cat "$logfile"
    echo "FAIL: Three-city simulation"
    exit 1
fi

cat "$logfile"

workdir="$(
    sed -n 's/^Results directory: //p' "$logfile" | tail -n 1
)"

if [ -z "$workdir" ] || [ ! -d "$workdir" ]; then
    echo "ERROR: Cannot locate simulation outputs"
    exit 1
fi

echo
echo "=========================================="
echo "CHECKING REGENERATED RESULTS"
echo "=========================================="

python3 "$checker" "$workdir"

echo
echo "=========================================="
echo "FINAL VERIFICATION"
echo "=========================================="

echo "Simulation execution: PASS"
echo "Numerical regression: PASS"
echo "Overall status: PASS"
echo "Results: $workdir"
