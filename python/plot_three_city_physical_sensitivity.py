from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

cities = ["mbabane", "manzini", "johannesburg"]
labels = ["Mbabane", "Manzini", "Johannesburg"]

parameters = [
    ("solar", "Solar radiation", "W/m²", 500.0),
    ("heat_capacity", "Effective heat capacity", "J/(m² K)", 2e6),
    ("emissivity", "Surface emissivity", "", 0.95),
    ("heat_transfer", "Heat-transfer coefficient", "W/(m² K)", 10.0),
    ("diffusivity", "Thermal diffusivity", "m²/s", 1e-5)
]

outdir = Path("figures/comparison")
outdir.mkdir(parents=True, exist_ok=True)

fig, axes = plt.subplots(
    2, 3, figsize=(16, 10),
    constrained_layout=True
)
axes = axes.flatten()

for city, label in zip(cities, labels):
    df = pd.read_csv(
        Path("results") / f"{city}_physical_sensitivity.csv"
    )

    reference = df[
        (df["parameter"] == "solar")
        & np.isclose(df["parameter_value"], 500.0)
    ]

    if len(reference) != 1:
        raise ValueError(f"Missing reference experiment: {city}")

    reference_cooling = float(
        reference["domain_mean_cooling_C"].iloc[0]
    )

    for ax, (key, title, unit, baseline) in zip(
        axes, parameters
    ):
        subset = df[
            df["parameter"] == key
        ].sort_values("parameter_value")

        if len(subset) != 3:
            raise ValueError(
                f"Expected three tests for {city}: {key}"
            )

        x = (
            subset["parameter_value"].to_numpy()
            / baseline
        )
        y = (
            subset["domain_mean_cooling_C"].to_numpy()
            / reference_cooling
        )

        ax.plot(
            x, y,
            marker="o",
            linewidth=2,
            label=label
        )

for ax, (_, title, unit, _) in zip(axes, parameters):
    ax.set_title(title)
    ax.set_xlabel("Parameter / reference value")
    ax.set_ylabel("Cooling / reference cooling")
    ax.axhline(1.0, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(1.0, color="gray", linestyle=":", alpha=0.5)
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=9)

axes[-1].axis("off")

fig.suptitle(
    "Three-city physical-parameter sensitivity of combined mitigation",
    fontsize=16
)

for ext in ("png", "pdf"):
    output = outdir / (
        f"three_city_physical_sensitivity.{ext}"
    )
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print("Saved:", output)

plt.close(fig)

print("\nREFERENCE S3 COOLING")
for city, label in zip(cities, labels):
    df = pd.read_csv(
        Path("results") / f"{city}_physical_sensitivity.csv"
    )
    ref = df[
        (df["parameter"] == "solar")
        & np.isclose(df["parameter_value"], 500.0)
    ]
    print(
        label,
        f"{ref['domain_mean_cooling_C'].iloc[0]:.6f} °C"
    )
