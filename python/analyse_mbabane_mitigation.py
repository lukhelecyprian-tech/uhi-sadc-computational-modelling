
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path.home() / "UHI_Thesis_Johannesburg" / "physics_model"
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

FIGURES.mkdir(parents=True, exist_ok=True)

REFERENCE = DATA / "Mbabane_LST.tif"

SCENARIOS = {
    "S0": RESULTS / "mbabane_mitigation_S0_K.txt",
    "S1": RESULTS / "mbabane_mitigation_S1_K.txt",
    "S2": RESULTS / "mbabane_mitigation_S2_K.txt",
    "S3": RESULTS / "mbabane_mitigation_S3_K.txt",
}

with rasterio.open(REFERENCE) as src:
    profile = src.profile.copy()
    height, width = src.height, src.width
    bounds = src.bounds
    transform = src.transform

def read_grid(path, dtype=float):
    arr = np.loadtxt(path, dtype=dtype)
    if arr.shape != (height, width):
        raise ValueError(
            f"{path.name}: expected {(height, width)}, "
            f"received {arr.shape}"
        )
    return arr

mask = read_grid(
    DATA / "mbabane_energy_mask.txt", int
) == 1

ndvi = read_grid(DATA / "mbabane_energy_ndvi.txt")
eligible = mask & (ndvi < 0.30)

temperatures = {
    name: read_grid(path) - 273.15
    for name, path in SCENARIOS.items()
}

baseline = temperatures["S0"]

if not np.all(np.isfinite(baseline[mask])):
    raise ValueError("Non-finite baseline temperatures")

if not np.isclose(
    baseline[mask].mean(),
    23.056383579857027,
    atol=1e-6,
):
    raise ValueError("Baseline mean does not match Fortran")

print("\nMBABANE MITIGATION ANALYSIS")
print("=" * 55)
print(f"Active pixels: {mask.sum():,}")
print(f"Eligible pixels: {eligible.sum():,}")
print(f"Eligible fraction: {100 * eligible.sum()/mask.sum():.3f}%")

records = []
cooling_fields = {}

profile.update(
    driver="GTiff",
    dtype="float32",
    count=1,
    nodata=-9999.0,
    compress="lzw",
)

for scenario in ("S1", "S2", "S3"):
    cooling = baseline - temperatures[scenario]
    cooling_fields[scenario] = cooling

    if not np.all(np.isfinite(cooling[mask])):
        raise ValueError(f"Non-finite cooling in {scenario}")

    domain_mean = cooling[mask].mean()
    eligible_mean = cooling[eligible].mean()
    max_cooling = cooling[mask].max()
    min_cooling = cooling[mask].min()

    records.append({
        "scenario": scenario,
        "domain_mean_cooling_C": domain_mean,
        "eligible_mean_cooling_C": eligible_mean,
        "maximum_cooling_C": max_cooling,
        "minimum_cooling_C": min_cooling,
        "eligible_pixels": int(eligible.sum()),
    })

    output = RESULTS / f"Mbabane_{scenario}_Cooling.tif"
    raster = np.full((height, width), -9999.0,
                     dtype=np.float32)
    raster[mask] = cooling[mask].astype(np.float32)

    with rasterio.open(output, "w", **profile) as dst:
        dst.write(raster, 1)

    print(f"\n{scenario}")
    print(f"  Domain mean cooling:   {domain_mean:.8f} C")
    print(f"  Eligible mean cooling: {eligible_mean:.8f} C")
    print(f"  Maximum cooling:       {max_cooling:.8f} C")
    print(f"  Minimum cooling:       {min_cooling:.8f} C")
    print(f"  Saved: {output.name}")

summary = pd.DataFrame(records)
summary.to_csv(
    RESULTS / "mbabane_mitigation_detailed_summary.csv",
    index=False,
)

# Test whether combined cooling is approximately additive
interaction = (
    cooling_fields["S3"]
    - cooling_fields["S1"]
    - cooling_fields["S2"]
)

print("\nCOMBINED-SCENARIO INTERACTION")
print(
    "Domain mean interaction:",
    f"{interaction[mask].mean():.10f} C"
)
print(
    "Eligible mean interaction:",
    f"{interaction[eligible].mean():.10f} C"
)

# Consistent colour scale across the three cooling maps
maximum = max(
    np.max(cooling_fields[s][mask])
    for s in ("S1", "S2", "S3")
)

fig, axes = plt.subplots(
    1, 3, figsize=(16, 5), constrained_layout=True
)

titles = {
    "S1": "S1: Reflective surfaces",
    "S2": "S2: Increased latent heat",
    "S3": "S3: Combined mitigation",
}

last_image = None

for ax, scenario in zip(axes, ("S1", "S2", "S3")):
    field = np.ma.masked_where(
        ~mask, cooling_fields[scenario]
    )

    last_image = ax.imshow(
        field,
        cmap="viridis",
        norm=Normalize(vmin=0, vmax=maximum),
        interpolation="nearest",
        extent=[
            bounds.left, bounds.right,
            bounds.bottom, bounds.top
        ],
    )

    ax.set_title(titles[scenario])
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_aspect("equal")

fig.colorbar(
    last_image,
    ax=axes,
    shrink=0.8,
    label="Cooling relative to baseline (°C)"
)

fig.suptitle(
    "Mbabane: Illustrative cooling after 6 hours",
    fontsize=14
)

figure_path = FIGURES / "mbabane_mitigation_cooling_maps.png"
fig.savefig(figure_path, dpi=300, bbox_inches="tight")
plt.close(fig)

print("\nSaved detailed summary and cooling figure.")
print("Analysis completed successfully.")
