from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

cities = ["mbabane", "manzini", "johannesburg"]
labels = ["Mbabane", "Manzini", "Johannesburg"]
scenarios = [1, 2, 3]

outdir = Path("figures/comparison")
outdir.mkdir(parents=True, exist_ok=True)

cooling = []
coverage = []

for city in cities:
    path = Path("results") / f"{city}_mitigation_summary.csv"
    df = pd.read_csv(path).set_index("scenario")

    required = {0, 1, 2, 3}
    if not required.issubset(set(df.index)):
        raise ValueError(f"Missing scenarios for {city}")

    baseline = float(df.loc[0, "final_mean_C"])

    cooling.append([
        baseline - float(df.loc[s, "final_mean_C"])
        for s in scenarios
    ])

    active = int(df.loc[0, "active_pixels"])
    eligible = int(df.loc[0, "eligible_pixels"])
    coverage.append(100 * eligible / active)

cooling = np.array(cooling)
x = np.arange(len(cities))
width = 0.24

fig, ax = plt.subplots(figsize=(10, 6))

scenario_labels = [
    "S1: Reflective surfaces",
    "S2: Increased latent heat",
    "S3: Combined mitigation"
]

for i, label in enumerate(scenario_labels):
    bars = ax.bar(
        x + (i - 1) * width,
        cooling[:, i],
        width,
        label=label
    )
    ax.bar_label(
        bars,
        fmt="%.3f",
        padding=3,
        fontsize=9
    )

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Domain-wide mean cooling relative to S0 (°C)")
ax.set_title(
    "Simulated urban heat mitigation across three Southern African cities"
)
ax.legend(frameon=False)
ax.grid(axis="y", alpha=0.2)
ax.set_axisbelow(True)
ax.set_ylim(0, max(cooling.flatten()) * 1.25)

fig.tight_layout()

for ext in ["png", "pdf"]:
    output = outdir / f"three_city_mitigation_comparison.{ext}"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print("Saved:", output)

plt.close(fig)

print("\nTHREE-CITY MITIGATION SUMMARY")
for i, city in enumerate(labels):
    print(
        f"{city:14s} "
        f"S1={cooling[i,0]:.6f}  "
        f"S2={cooling[i,1]:.6f}  "
        f"S3={cooling[i,2]:.6f}  "
        f"Eligible={coverage[i]:.2f}%"
    )
