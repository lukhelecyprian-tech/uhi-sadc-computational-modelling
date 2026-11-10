from pathlib import Path
import numpy as np
import rasterio
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import ScalarFormatter

root = Path(".")
outdir = root / "figures" / "improved"
outdir.mkdir(parents=True, exist_ok=True)

raster_path = Path(
    "/mnt/c/Users/Cyprian Lukhele/Desktop/"
    "MSc_UHI_Thesis/Johannesburg/Raster/"
    "Johannesburg_LST.tif"
)

boundary_path = Path(
    "/mnt/c/Users/Cyprian Lukhele/Desktop/"
    "MSc_UHI_Thesis/Johannesburg/Johannesburg.gpkg"
)

with rasterio.open(raster_path) as src:
    shape = src.shape
    crs = src.crs
    transform = src.transform
    bounds = src.bounds
    extent = (
        bounds.left, bounds.right,
        bounds.bottom, bounds.top
    )

with rasterio.open(
    root / "data/Johannesburg_four_variable_mask.tif"
) as src:
    if (
        src.shape != shape
        or src.crs != crs
        or src.transform != transform
    ):
        raise ValueError("Modelling mask is not aligned")

    valid = src.read(1) == 1

boundary = gpd.read_file(boundary_path).to_crs(crs)

def load_temperature(scenario):
    path = (
        root / "results"
        / f"johannesburg_mitigation_{scenario}_K.txt"
    )

    data = np.loadtxt(path)

    if data.size != valid.size:
        raise ValueError(f"Grid mismatch: {scenario}")

    return data.reshape(shape)

baseline = load_temperature("S0")

scenarios = [
    ("S1", "Reflective surfaces"),
    ("S2", "Enhanced latent heat"),
    ("S3", "Combined mitigation")
]

fig, axes = plt.subplots(
    1, 3,
    figsize=(17, 6),
    constrained_layout=True
)

norm = Normalize(vmin=0, vmax=0.8)

for ax, (scenario, label) in zip(axes, scenarios):

    temperature = load_temperature(scenario)
    cooling = baseline - temperature

    cooling = np.ma.masked_where(
        ~valid | ~np.isfinite(cooling),
        cooling
    )

    image = ax.imshow(
        cooling,
        extent=extent,
        origin="upper",
        cmap="viridis",
        norm=norm,
        interpolation="nearest"
    )

    boundary.boundary.plot(
        ax=ax,
        color="black",
        linewidth=0.8,
        zorder=3
    )

    ax.set_title(
        f"{scenario} — {label}",
        fontsize=12,
        fontweight="bold"
    )

    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_aspect("equal")

    ax.tick_params(labelsize=8)

    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

    print(
        f"{scenario}: "
        f"mean cooling = {cooling.mean():.6f} °C; "
        f"maximum = {cooling.max():.6f} °C"
    )

fig.suptitle(
    "Johannesburg: Simulated surface cooling after six hours",
    fontsize=16,
    fontweight="bold"
)

cbar = fig.colorbar(
    image,
    ax=axes,
    shrink=0.78,
    pad=0.02
)

cbar.set_label(
    "Cooling relative to baseline S0 (°C)"
)

for extension in ["png", "pdf"]:
    output = (
        outdir /
        f"johannesburg_cooling_improved.{extension}"
    )

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight"
    )

    print("Saved:", output)

plt.close(fig)
