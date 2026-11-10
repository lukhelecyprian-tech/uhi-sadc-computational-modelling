from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

CITIES = {
    "Mbabane": "mbabane",
    "Manzini": "manzini",
    "Johannesburg": "johannesburg",
}

OUT = Path("figures/comparison")
OUT.mkdir(parents=True, exist_ok=True)

data = {}

for city, filename in CITIES.items():
    path = Path("results") / f"{filename}_sensitivity_analysis.csv"
    df = pd.read_csv(path)

    required = {
        "parameter",
        "ndvi_threshold",
        "delta_albedo",
        "delta_latent_Wm2",
        "domain_mean_cooling_C",
    }

    if not required.issubset(df.columns):
        raise ValueError(f"Missing columns in {path}")

    data[city] = df

fig, axes = plt.subplots(
    1, 2,
    figsize=(13, 5),
    constrained_layout=True
)

experiments = [
    (
        "albedo_increase",
        "delta_albedo",
        "Increase in surface albedo",
        "Albedo sensitivity",
    ),
    (
        "latent_increase",
        "delta_latent_Wm2",
        "Additional latent heat flux (W/m²)",
        "Latent heat sensitivity",
    ),
]

for ax, (parameter, xcol, xlabel, title) in zip(
    axes, experiments
):
    for city, df in data.items():

        subset = df[
            (df["parameter"] == parameter)
            & (df["ndvi_threshold"] == 0.3)
        ].sort_values(xcol)

        if len(subset) != 3:
            raise ValueError(
                f"{city}: expected 3 rows for {parameter}, "
                f"found {len(subset)}"
            )

        ax.plot(
            subset[xcol],
            subset["domain_mean_cooling_C"],
            marker="o",
            linewidth=2,
            markersize=6,
            label=city,
        )

    ax.set_title(title, fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Domain-average cooling relative to S0 (°C)")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    ax.legend()

fig.suptitle(
    "Three-city sensitivity to mitigation strength",
    fontsize=15,
    fontweight="bold",
)

png = OUT / "three_city_mitigation_sensitivity.png"
pdf = OUT / "three_city_mitigation_sensitivity.pdf"

fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")

plt.close(fig)

print("Generated:", png)
print("Generated:", pdf)
print("All three cities included.")
