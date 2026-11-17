#!/usr/bin/env python3
"""
Three-City UHI Mitigation Regression Test

MSc Computational Physics
University of Eswatini

Checks:
1. Scenario summary CSV structure
2. Baseline and mitigation cooling values
3. Agreement with Chapter 4 reference results
4. Scenario ordering and physical consistency

Uses Python standard library only.
"""

import csv
import math
import sys
from pathlib import Path

REFERENCE = {
    "mbabane": 0.0624,
    "manzini": 0.0188,
    "johannesburg": 0.2200,
}

TOLERANCE = 0.0001

def read_summary(path):
    if not path.is_file():
        raise ValueError(f"Missing summary: {path}")

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)

        required = {
            "scenario",
            "active_pixels",
            "eligible_pixels",
            "initial_mean_C",
            "final_mean_C",
            "mean_change_C",
        }

        if reader.fieldnames is None:
            raise ValueError(f"Missing CSV header: {path}")

        if not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"Incorrect CSV columns: {path}")

        records = {}

        for row in reader:
            scenario = int(row["scenario"])

            if scenario in records:
                raise ValueError(
                    f"Duplicate scenario {scenario}: {path}"
                )

            records[scenario] = {
                "active": int(row["active_pixels"]),
                "eligible": int(row["eligible_pixels"]),
                "initial": float(row["initial_mean_C"]),
                "final": float(row["final_mean_C"]),
                "change": float(row["mean_change_C"]),
            }

    if set(records) != {0, 1, 2, 3}:
        raise ValueError(f"Expected scenarios 0-3: {path}")

    return records

def check_city(city, directory):
    path = (
        directory
        / city
        / "results"
        / f"{city}_mitigation_summary.csv"
    )

    records = read_summary(path)

    baseline = records[0]["final"]

    cooling = {}

    for scenario in (1, 2, 3):
        cooling[scenario] = (
            baseline - records[scenario]["final"]
        )

    for scenario, record in records.items():
        if not all(
            math.isfinite(record[key])
            for key in ("initial", "final", "change")
        ):
            raise ValueError(
                f"{city}: non-finite values in S{scenario}"
            )

        if record["active"] <= 0:
            raise ValueError(f"{city}: invalid active pixels")

        if not 0 <= record["eligible"] <= record["active"]:
            raise ValueError(f"{city}: invalid eligible pixels")

        calculated_change = (
            record["final"] - record["initial"]
        )

        if abs(calculated_change - record["change"]) > 1e-8:
            raise ValueError(
                f"{city}: inconsistent mean change"
            )

    if not (
        cooling[1] >= 0
        and cooling[2] >= 0
        and cooling[3] >= cooling[1]
        and cooling[3] >= cooling[2]
    ):
        raise ValueError(
            f"{city}: unexpected mitigation ordering"
        )

    observed = cooling[3]
    expected = REFERENCE[city]

    difference = abs(observed - expected)

    print(f"\nCITY: {city.upper()}")
    print(f"Active pixels: {records[0]['active']}")
    print(f"Eligible pixels: {records[0]['eligible']}")
    print(f"S1 cooling: {cooling[1]:.8f} C")
    print(f"S2 cooling: {cooling[2]:.8f} C")
    print(f"S3 cooling: {cooling[3]:.8f} C")
    print(f"Reference S3: {expected:.8f} C")
    print(f"Absolute difference: {difference:.8e} C")

    if difference > TOLERANCE:
        raise ValueError(
            f"{city}: S3 regression tolerance exceeded"
        )

    print("STATUS: PASS")

def main():
    if len(sys.argv) != 2:
        print(
            "Usage: python3 check_three_city_results.py "
            "RESULTS_DIRECTORY"
        )
        return 2

    directory = Path(sys.argv[1])

    if not directory.is_dir():
        print(f"ERROR: Directory not found: {directory}")
        return 1

    print("======================================")
    print("THREE-CITY MITIGATION REGRESSION")
    print("======================================")

    try:
        for city in REFERENCE:
            check_city(city, directory)
    except (OSError, ValueError, TypeError) as error:
        print(f"\nFAIL: {error}", file=sys.stderr)
        return 1

    print("\n======================================")
    print("ALL THREE CITY REGRESSION TESTS PASSED")
    print("======================================")

    return 0

if __name__ == "__main__":
    sys.exit(main())
