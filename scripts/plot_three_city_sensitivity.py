#!/usr/bin/env python3
"""Generate thesis figures from three-city sensitivity CSV results."""

import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = (
    Path.home()
    / "UHI_Thesis_Johannesburg"
    / "physics_model"
    / "results"
)

FIGURES = RESULTS / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

CITIES = ["mbabane", "manzini", "johannesburg"]

LABELS = {
    "mbabane": "Mbabane",
    "manzini": "Manzini",
    "johannesburg": "Johannesburg",
}

COLORS = {
    "mbabane": "#2166AC",
    "manzini": "#D6604D",
    "johannesburg": "#4D9221",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "legend.fontsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def read_csv(filename):
    path = RESULTS / filename
    if not path.is_file():
        raise FileNotFoundError(path)

    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def save_figure(fig, name):
    for extension in ("png", "pdf"):
        path = FIGURES / f"{name}.{extension}"
        fig.savefig(
            path,
            dpi=300,
            bbox_inches="tight"
        )
        print("Saved:", path)

    plt.close(fig)


def plot_normalized():
    rows = read_csv(
        "three_city_normalized_sensitivity.csv"
    )

    parameters = [
        "heat_capacity",
        "solar",
        "heat_transfer",
        "emissivity",
        "diffusivity",
    ]

    names = [
        "Heat capacity",
        "Solar radiation",
        "Heat transfer",
        "Emissivity",
        "Diffusivity",
    ]

    lookup = {
        row["parameter"]: row
        for row in rows
    }

    y = np.arange(len(parameters))
    width = 0.24

    fig, ax = plt.subplots(figsize=(9, 5.5))

    for index, city in enumerate(CITIES):
        values = [
            float(
                lookup[p][f"{city}_sensitivity"]
            )
            for p in parameters
        ]

        ax.barh(
            y + (index - 1) * width,
            values,
            height=width,
            label=LABELS[city],
            color=COLORS[city],
        )

    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.invert_yaxis()

    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Normalized sensitivity coefficient")
    ax.set_title(
        "Physical parameter sensitivity across three cities"
    )
    ax.legend(frameon=False)
    ax.grid(axis="x", alpha=0.2)
    ax.set_axisbelow(True)

    fig.tight_layout()

    save_figure(
        fig,
        "three_city_normalized_sensitivity"
    )



def plot_mitigation():
    data = {
        city: read_csv(f"{city}_sensitivity_analysis.csv")
        for city in CITIES
    }

    experiments = [
        ("NDVI_threshold", "ndvi_threshold",
         "NDVI eligibility threshold", "(a) NDVI threshold"),
        ("albedo_increase", "delta_albedo",
         "Albedo increase", "(b) Albedo modification"),
        ("latent_increase", "delta_latent_Wm2",
         "Latent heat increase (W/m²)", "(c) Latent heat enhancement"),
    ]

    fig, axes = plt.subplots(
        1, 3, figsize=(13.5, 4.8),
        sharey=True
    )

    for ax, (parameter, xcolumn, xlabel, title) in zip(
        axes, experiments
    ):
        for city in CITIES:
            rows = [
                row for row in data[city]
                if row["parameter"].strip() == parameter
            ]

            pairs = sorted(
                (
                    float(row[xcolumn]),
                    float(row["domain_mean_cooling_C"])
                )
                for row in rows
            )

            ax.plot(
                [point[0] for point in pairs],
                [point[1] for point in pairs],
                marker="o",
                linewidth=2,
                markersize=5,
                label=LABELS[city],
                color=COLORS[city],
            )

        ax.set_title(title, fontsize=11, pad=12)
        ax.set_xlabel(xlabel)
        ax.set_ylim(0, 0.40)
        ax.grid(alpha=0.25)
        ax.set_axisbelow(True)

    fig.supylabel(
        "Domain-mean mitigation cooling (°C)",
        x=0.015,
        fontsize=11
    )

    handles, labels = axes[0].get_legend_handles_labels()

    fig.legend(
        handles, labels,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.01)
    )

    fig.suptitle(
        "Mitigation sensitivity across three SADC cities",
        fontsize=13,
        y=0.99
    )

    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        top=0.82,
        bottom=0.23,
        wspace=0.18
    )

    save_figure(
        fig,
        "three_city_mitigation_sensitivity"
    )


if __name__ == "__main__":
    plot_normalized()
    plot_mitigation()
    print("\nAll sensitivity figures generated successfully.")
