#!/usr/bin/env bash
set -uo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workdir="$(mktemp -d /tmp/uhi-verification.XXXXXX)"

programs=(
    verify_heat2d
    grid_convergence
    time_convergence
    verify_coupled_spatial_operator
)

echo "======================================"
echo "UHI NUMERICAL VERIFICATION SYSTEM"
echo "======================================"

echo
echo "Building Fortran programs..."

if ! make -C "$repo_root" -j2; then
    echo "FAIL: Fortran compilation"
    exit 1
fi

mkdir -p \
    "$workdir/results" \
    "$workdir/figures" \
    "$workdir/data"

passed=0
failed=0

for program in "${programs[@]}"; do

    echo
    echo "===== $program ====="

    (
        cd "$workdir"
        timeout 120s "$repo_root/build/$program"
    ) > "$workdir/${program}.log" 2>&1

    status=$?

    if [ "$status" -eq 0 ]; then
        echo "PASS: $program"
        passed=$((passed + 1))
    else
        echo "FAIL: $program (exit code $status)"
        tail -n 20 "$workdir/${program}.log"
        failed=$((failed + 1))
    fi

done

echo
echo "======================================"
echo "FORTRAN EXECUTION SUMMARY"
echo "======================================"

echo "Passed: $passed"
echo "Failed: $failed"

if [ "$failed" -ne 0 ]; then
    echo "FAIL: One or more Fortran programs failed"
    echo "Logs: $workdir"
    exit 1
fi

echo
echo "======================================"
echo "NUMERICAL REGRESSION VERIFICATION"
echo "======================================"

if ! python3 \
    "$repo_root/scripts/check_numerical_results.py" \
    "$workdir/results"; then

    echo
    echo "FAIL: Numerical regression verification"
    echo "Results: $workdir/results"
    exit 1
fi

echo
echo "======================================"
echo "FINAL VERIFICATION SUMMARY"
echo "======================================"

echo "Fortran execution: 4/4 passed"
echo "Numerical regression: PASSED"
echo "Overall status: PASS"
echo "Logs and results: $workdir"

exit 0
