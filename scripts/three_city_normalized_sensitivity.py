#!/usr/bin/env python3

import csv
import math
from pathlib import Path

RESULTS = (
    Path.home()
    / "UHI_Thesis_Johannesburg"
    / "physics_model"
    / "results"
)

CITIES = ["mbabane", "manzini", "johannesburg"]

PARAMETERS = [
    "solar",
    "heat_capacity",
    "emissivity",
    "heat_transfer",
    "diffusivity",
]

OUTPUT = RESULTS / "three_city_normalized_sensitivity.csv"


def read_results(city):
    path = RESULTS / f"{city}_physical_sensitivity.csv"

    if not path.is_file():
        raise FileNotFoundError(path)

    grouped = {}

    with path.open(newline="") as file:
        reader = csv.DictReader(file)

        for row in reader:
            parameter = row["parameter"].strip()

            value = float(row["parameter_value"])
            cooling = float(row["domain_mean_cooling_C"])

            if not math.isfinite(value) or not math.isfinite(cooling):
                raise ValueError(
                    f"Non-finite value in {path}"
                )

            grouped.setdefault(parameter, []).append(
                (value, cooling)
            )

    return grouped


def normalized_sensitivity(points):
    points = sorted(points, key=lambda x: x[0])

    if len(points) != 3:
        raise ValueError(
            "Expected exactly three parameter values"
        )

    (p_minus, c_minus), (p0, c0), (p_plus, c_plus) = points

    if p_plus == p_minus:
        raise ValueError("Zero parameter interval")

    if c0 == 0:
        raise ValueError("Zero baseline cooling")

    coefficient = (
        (p0 / c0)
        * (c_plus - c_minus)
        / (p_plus - p_minus)
    )

    return coefficient


def main():
    results = {
        city: read_results(city)
        for city in CITIES
    }

    rows = []

    for parameter in PARAMETERS:
        row = {"parameter": parameter}

        for city in CITIES:
            if parameter not in results[city]:
                raise ValueError(
                    f"Missing {parameter} for {city}"
                )

            row[f"{city}_sensitivity"] = (
                normalized_sensitivity(
                    results[city][parameter]
                )
            )

        rows.append(row)

    fieldnames = [
        "parameter",
        "mbabane_sensitivity",
        "manzini_sensitivity",
        "johannesburg_sensitivity",
    ]

    with OUTPUT.open("w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames
        )

        writer.writeheader()
        writer.writerows(rows)

    print("\nTHREE-CITY NORMALIZED SENSITIVITY")
    print("=" * 76)

    print(
        f"{'Parameter':<18}"
        f"{'Mbabane':>18}"
        f"{'Manzini':>18}"
        f"{'Johannesburg':>18}"
    )

    print("-" * 76)

    for row in rows:
        print(
            f"{row['parameter']:<18}"
            f"{row['mbabane_sensitivity']:>18.8f}"
            f"{row['manzini_sensitivity']:>18.8f}"
            f"{row['johannesburg_sensitivity']:>18.8f}"
        )

    print("\nOutput:", OUTPUT)


if __name__ == "__main__":
    main()
