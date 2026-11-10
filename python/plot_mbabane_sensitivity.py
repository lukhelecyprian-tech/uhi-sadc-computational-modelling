
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------------------------------
# MBABANE PARAMETER SENSITIVITY ANALYSIS
# Scientific figures for MSc Computational Physics thesis
# -------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
CSV_FILE = ROOT / "results" / "mbabane_sensitivity_analysis.csv"
FIG_DIR = ROOT / "figures"

FIG_DIR.mkdir(parents=True, exist_ok=True)

if not CSV_FILE.exists():
    raise FileNotFoundError(f"Missing sensitivity data: {CSV_FILE}")

df = pd.read_csv(CSV_FILE)

required = [
    "test",
    "parameter",
    "ndvi_threshold",
    "delta_albedo",
    "delta_latent_Wm2",
    "eligible_pixels",
    "domain_mean_cooling_C",
    "eligible_mean_cooling_C",
]

missing = [col for col in required if col not in df.columns]
if missing:
    raise ValueError(f"Missing CSV columns: {missing}")

if len(df) != 9:
    raise ValueError("Expected nine sensitivity tests.")

if df[required[2:]].isna().any().any():
    raise ValueError("Missing numerical values in CSV.")

# Plot settings
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.fontsize": 10,
    "figure.dpi": 120,
    "savefig.dpi": 300,
})

def make_figure(parameter, xcolumn, xlabel, filename):
    subset = df.loc[df["parameter"] == parameter].copy()
    subset = subset.sort_values(xcolumn)

    if len(subset) != 3:
        raise ValueError(
            f"Expected three tests for {parameter}; got {len(subset)}"
        )

    x = subset[xcolumn].to_numpy(dtype=float)
    domain = subset["domain_mean_cooling_C"].to_numpy(dtype=float)
    eligible = subset["eligible_mean_cooling_C"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(8.2, 5.5))

    ax.plot(
        x, domain,
        marker="o",
        linewidth=2,
        markersize=7,
        label="Entire valid-data footprint",
    )

    ax.plot(
        x, eligible,
        marker="s",
        linewidth=2,
        markersize=7,
        label="Intervention pixels only",
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Mean simulated cooling after 6 h (°C)")
    ax.set_title("Mbabane: mitigation parameter sensitivity")

    ax.set_xticks(x)
    ax.set_ylim(bottom=0)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(loc="best")

    # Mark the original parameter setting
    baseline_x = {
        "NDVI_threshold": 0.30,
        "albedo_increase": 0.10,
        "latent_increase": 30.0,
    }[parameter]

    ax.axvline(
        baseline_x,
        linestyle=":",
        linewidth=1.5,
        label="Original S3 setting",
    )

    # Use a compact annotation rather than overcrowding
    ax.annotate(
        "Original S3",
        xy=(baseline_x, 0),
        xycoords=("data", "axes fraction"),
        xytext=(5, 8),
        textcoords="offset points",
        rotation=90,
        fontsize=9,
        va="bottom",
    )

    fig.tight_layout()

    output = FIG_DIR / filename
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {output}")

make_figure(
    parameter="NDVI_threshold",
    xcolumn="ndvi_threshold",
    xlabel="NDVI eligibility threshold (dimensionless)",
    filename="mbabane_sensitivity_ndvi_threshold.png",
)

make_figure(
    parameter="albedo_increase",
    xcolumn="delta_albedo",
    xlabel="Increase in surface albedo (dimensionless)",
    filename="mbabane_sensitivity_albedo.png",
)

make_figure(
    parameter="latent_increase",
    xcolumn="delta_latent_Wm2",
    xlabel="Additional latent heat flux (W/m²)",
    filename="mbabane_sensitivity_latent_heat.png",
)

print()
print("========================================")
print("MBABANE SENSITIVITY FIGURES COMPLETED")
print("========================================")
print("Number of figures: 3")
print("Output directory:", FIG_DIR)
print()
print("Important interpretation:")
print("NDVI threshold changes intervention coverage.")
print("Albedo and latent flux change intervention strength.")
print("Results are numerical sensitivity experiments,")
print("not physically calibrated cooling predictions.")
