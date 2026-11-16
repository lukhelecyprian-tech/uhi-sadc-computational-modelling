#!/usr/bin/env python3

"""
Numerical Regression Verification
MSc Computational Physics
University of Eswatini

Tests:
1. Analytical heat equation verification
2. Spatial grid convergence
3. Temporal convergence
4. Coupled spatial operator convergence

Uses only the Python standard library.
"""

import csv
import math
import sys
from pathlib import Path


def read_table(path, columns):
    """Read numerical CSV data and validate columns."""

    with path.open(newline="") as file:
        reader = csv.DictReader(file)

        if reader.fieldnames is None:
            raise ValueError(
                f"{path.name}: CSV header is missing"
            )

        # Remove whitespace from column names.
        reader.fieldnames = [
            name.strip() for name in reader.fieldnames
        ]

        missing = [
            column
            for column in columns
            if column not in reader.fieldnames
        ]

        if missing:
            raise ValueError(
                f"{path.name}: missing columns {missing}"
            )

        rows = []

        for row in reader:
            values = [
                float(row[column])
                for column in columns
            ]

            if not all(math.isfinite(v) for v in values):
                raise ValueError(
                    f"{path.name}: non-finite numerical values"
                )

            rows.append(values)

    if len(rows) < 2:
        raise ValueError(
            f"{path.name}: insufficient numerical data"
        )

    return rows


def check_orders(
    path,
    spacing_column,
    lower,
    upper,
    error_column="rmse"
):
    """Calculate and verify observed convergence orders."""

    rows = read_table(
        path,
        [spacing_column, error_column]
    )

    orders = []

    for previous, current in zip(rows, rows[1:]):

        h1, e1 = previous
        h2, e2 = current

        if not (
            h1 > h2 > 0
            and e1 > 0
            and e2 > 0
        ):
            raise ValueError(
                f"{path.name}: invalid refinement data"
            )

        order = (
            math.log(e1 / e2)
            / math.log(h1 / h2)
        )

        orders.append(order)

        if not lower <= order <= upper:
            raise ValueError(
                f"{path.name}: convergence order "
                f"{order:.6f} outside "
                f"[{lower}, {upper}]"
            )

    print(f"PASS: {path.name}")

    for index, order in enumerate(orders, start=1):
        print(
            f"  Refinement {index}: "
            f"observed order = {order:.6f}"
        )


def read_field(path):
    """Read a two-dimensional temperature field."""

    field = []

    with path.open(newline="") as file:

        reader = csv.reader(file)

        for row in reader:

            if not row:
                continue

            values = [
                float(value)
                for value in row
            ]

            if not values:
                raise ValueError(
                    f"{path.name}: empty field row"
                )

            if not all(math.isfinite(v) for v in values):
                raise ValueError(
                    f"{path.name}: non-finite field values"
                )

            field.append(values)

    if not field:
        raise ValueError(
            f"{path.name}: empty temperature field"
        )

    widths = {
        len(row)
        for row in field
    }

    if len(widths) != 1:
        raise ValueError(
            f"{path.name}: inconsistent field dimensions"
        )

    return field


def check_heat_field(directory):
    """Compare numerical and analytical temperature fields."""

    numerical = read_field(
        directory / "verification_numerical.csv"
    )

    analytical = read_field(
        directory / "verification_exact.csv"
    )

    if len(numerical) != len(analytical):
        raise ValueError(
            "Temperature field row counts differ"
        )

    if any(
        len(a) != len(b)
        for a, b in zip(numerical, analytical)
    ):
        raise ValueError(
            "Temperature field column counts differ"
        )

    errors = [
        numerical_value - analytical_value
        for numerical_row, analytical_row
        in zip(numerical, analytical)
        for numerical_value, analytical_value
        in zip(numerical_row, analytical_row)
    ]

    rmse = math.sqrt(
        math.fsum(error ** 2 for error in errors)
        / len(errors)
    )

    maximum_error = max(
        abs(error)
        for error in errors
    )

    rmse_limit = 1.0e-4
    maximum_error_limit = 1.0e-3

    print("\n===== HEAT EQUATION VERIFICATION =====")

    print(f"RMSE: {rmse:.10e}")
    print(f"Maximum error: {maximum_error:.10e}")

    if rmse >= rmse_limit:
        raise ValueError(
            "Heat equation RMSE exceeds tolerance"
        )

    if maximum_error >= maximum_error_limit:
        raise ValueError(
            "Heat equation maximum error exceeds tolerance"
        )

    print("PASS: Analytical heat equation verification")


def main():

    if len(sys.argv) != 2:
        print(
            "Usage: python3 check_numerical_results.py "
            "RESULTS_DIRECTORY"
        )
        return 2

    directory = Path(sys.argv[1])

    if not directory.is_dir():
        print(
            f"ERROR: Results directory not found: {directory}"
        )
        return 1

    try:

        # TEST 1: Analytical heat equation
        check_heat_field(directory)

        # TEST 2: Spatial convergence
        print("\n===== SPATIAL GRID CONVERGENCE =====")

        check_orders(
            path=directory / "convergence.csv",
            spacing_column="dx",
            lower=1.90,
            upper=2.10,
            error_column="rmse"
        )

        # TEST 3: Temporal convergence
        print("\n===== TEMPORAL CONVERGENCE =====")

        check_orders(
            path=directory / "time_convergence.csv",
            spacing_column="dt",
            lower=0.95,
            upper=1.05,
            error_column="rmse"
        )

        # TEST 4: Coupled spatial operator
        print("\n===== COUPLED SPATIAL OPERATOR =====")

        check_orders(
            path=(
                directory
                / "coupled_spatial_operator_verification.csv"
            ),
            spacing_column="dx_m",
            lower=1.95,
            upper=2.05,
            error_column="rmse_K_per_m2"
        )

    except (
        OSError,
        ValueError,
        ZeroDivisionError,
        OverflowError
    ) as error:

        print(
            f"\nFAIL: {error}",
            file=sys.stderr
        )

        return 1

    print("\n======================================")
    print("ALL NUMERICAL REGRESSION CHECKS PASSED")
    print("======================================")

    return 0


if __name__ == "__main__":
    sys.exit(main())
