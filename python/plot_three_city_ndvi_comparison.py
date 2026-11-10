from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

cities = ["mbabane", "manzini", "johannesburg"]
labels = ["Mbabane", "Manzini", "Johannesburg"]

outdir = Path("figures/comparison")
outdir.mkdir(parents=True, exist_ok=True)

fig, axes = plt.subplots(
    1, 2, figsize=(13, 5.5),
    constrained_layout=True
)

for city, label in zip(cities, labels):
    path = Path("results") / f"{city}_sensitivity_analysis.csv"
    df = pd.read_csv(path)

    subset = df[
        df["parameter"] == "NDVI_threshold"
    ].sort_values("ndvi_threshold")

    if len(subset) != 3:
        raise ValueError(f"Expected three NDVI tests for {city}")

    thresholds = subset["ndvi_threshold"].to_numpy()
    cooling = subset["domain_mean_cooling_C"].to_numpy()

    summary = pd.read_csv(
        Path("results") / f"{city}_mitigation_summary.csv"
    )
    active = int(
        summary.loc[
            summary["scenario"] == 0, "active_pixels"
        ].iloc[0]
    )

    eligible_percent = (
        100 * subset["eligible_pixels"].to_numpy() / active
    )

    axes[0].plot(
        thresholds, cooling,
        marker="o", linewidth=2, label=label
    )

    axes[1].plot(
        thresholds, eligible_percent,
        marker="s", linewidth=2, label=label
    )

    print(f"\n{label}")
    for t, c, p in zip(
        thresholds, cooling, eligible_percent
    ):
        print(
            f"NDVI={t:.2f}  "
            f"Cooling={c:.6f} °C  "
            f"Eligible={p:.2f}%"
        )

axes[0].set_title("A. Domain-wide cooling")
axes[0].set_ylabel("Mean cooling relative to S0 (°C)")

axes[1].set_title("B. Mitigation-eligible area")
axes[1].set_ylabel("Eligible fraction of valid pixels (%)")

for ax in axes:
    ax.set_xlabel("NDVI eligibility threshold")
    ax.set_xticks([0.2, 0.3, 0.4])
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

fig.suptitle(
    "NDVI threshold sensitivity across three cities",
    fontsize=15
)

for ext in ("png", "pdf"):
    output = outdir / f"three_city_ndvi_sensitivity.{ext}"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print("Saved:", output)

plt.close(fig)
