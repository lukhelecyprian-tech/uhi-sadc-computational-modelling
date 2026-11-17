#!/usr/bin/env bash
set -uo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_data="$HOME/UHI_Thesis_Johannesburg/physics_model/data"

cities=(mbabane manzini johannesburg)
suffixes=(grid temperature_K ndvi albedo reflectance mask)

echo "=========================================="
echo "THREE-CITY UHI REPRODUCIBILITY TEST"
echo "=========================================="

if [ ! -d "$source_data" ]; then
    echo "ERROR: Original data directory not found"
    exit 1
fi

echo
echo "Checking original input files..."

for city in "${cities[@]}"; do
    for suffix in "${suffixes[@]}"; do
        filename="${city}_energy_${suffix}.txt"

        if [ ! -f "$source_data/$filename" ]; then
            echo "MISSING: $filename"
            exit 1
        fi
    done
done

echo "PASS: All 18 input files available"

echo
echo "Compiling mitigation programs..."

for city in "${cities[@]}"; do
    if ! make -C "$repo_root" \
        "build/${city}_mitigation_scenarios"; then
        echo "ERROR: Compilation failed for $city"
        exit 1
    fi
done

workdir="$(mktemp -d /tmp/uhi-three-city.XXXXXX)" || exit 1

mkdir -p "$workdir/data"

echo
echo "Preparing isolated input data..."

for city in "${cities[@]}"; do
    for suffix in "${suffixes[@]}"; do
        filename="${city}_energy_${suffix}.txt"

        ln -s \
            "$source_data/$filename" \
            "$workdir/data/$filename"
    done
done

(
    cd "$workdir"
    sha256sum data/* > input_checksums.sha256
)

echo "Input checksums recorded."

passed=0
failed=0

for city in "${cities[@]}"; do

    echo
    echo "=========================================="
    echo "RUNNING: $city"
    echo "=========================================="

    citydir="$workdir/$city"
    mkdir -p "$citydir/results"

    ln -s "$workdir/data" "$citydir/data"

    (
        cd "$citydir"
        "$repo_root/build/${city}_mitigation_scenarios"
    ) > "$citydir/run.log" 2>&1

    status=$?

    if [ "$status" -eq 0 ]; then

        expected=(
            "${city}_mitigation_summary.csv"
            "${city}_mitigation_S0_K.txt"
            "${city}_mitigation_S1_K.txt"
            "${city}_mitigation_S2_K.txt"
            "${city}_mitigation_S3_K.txt"
        )

        complete=1

        for filename in "${expected[@]}"; do
            if [ ! -s "$citydir/results/$filename" ]; then
                echo "MISSING OR EMPTY OUTPUT: $filename"
                complete=0
            fi
        done

        if [ "$complete" -eq 1 ]; then
            echo "PASS: $city"
            passed=$((passed + 1))
            echo "Summary:"
            head -n 5 \
                "$citydir/results/${city}_mitigation_summary.csv"
        else
            echo "FAIL: $city - incomplete outputs"
            failed=$((failed + 1))
        fi

    else
        echo "FAIL: $city (exit code $status)"
        tail -n 20 "$citydir/run.log"
        failed=$((failed + 1))
    fi

done

echo
echo "=========================================="
echo "FINAL REPRODUCIBILITY RUN SUMMARY"
echo "=========================================="

echo "Successful city runs: $passed"
echo "Failed city runs: $failed"
echo "Results directory: $workdir"
echo "Input checksums: $workdir/input_checksums.sha256"

if [ "$failed" -ne 0 ]; then
    exit 1
fi

echo "Overall execution status: PASS"
echo "Note: Comparison with Chapter 4 values is still required."
