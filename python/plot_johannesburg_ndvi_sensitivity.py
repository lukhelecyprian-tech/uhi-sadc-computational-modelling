from pathlib import Path
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

root = Path(".")
outdir = root / "figures" / "sensitivity"
outdir.mkdir(parents=True, exist_ok=True)

raster_path = Path(
    "/mnt/c/Users/Cyprian Lukhele/Desktop/"
    "MSc_UHI_Thesis/Johannesburg/Raster/"
    "Johannesburg_LST.tif"
)

mask_path = root / "data/Johannesburg_four_variable_mask.tif"

with rasterio.open(raster_path) as src:
    shape = src.shape
    extent = (
        src.bounds.left,
        src.bounds.right,
        src.bounds.bottom,
        src.bounds.top
    )
    transform = src.transform
    crs = src.crs

with rasterio.open(mask_path) as src:
    if (
        src.shape != shape
        or src.transform != transform
        or src.crs != crs
    ):
        raise ValueError("Raster and mask are not aligned")

    valid = src.read(1) == 1

fig, axes = plt.subplots(
    1, 3,
    figsize=(17, 6),
    constrained_layout=True
)

thresholds = [
    ("020", "0.20"),
    ("030", "0.30"),
    ("040", "0.40")
]

norm = Normalize(vmin=0, vmax=0.8)

for ax, (code, threshold) in zip(axes, thresholds):
    path = root / "results" / (
        f"johannesburg_ndvi{code}_cooling.txt"
    )

    cooling = np.loadtxt(path)

    if cooling.size != valid.size:
        raise ValueError(f"Grid mismatch: {path}")

    cooling = cooling.reshape(shape)

    cooling = np.ma.masked_where(
        ~valid | ~np.isfinite(cooling),
        cooling
    )

    image = ax.imshow(
        cooling,
        origin="upper",
        extent=extent,
        cmap="viridis",
        norm=norm,
        interpolation="nearest"
    )

    ax.set_title(
        f"NDVI threshold < {threshold}",
        fontsize=12
    )
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_aspect("equal")

    print(
        f"NDVI {threshold}:",
        f"Mean cooling = {cooling.mean():.6f} °C;",
        f"Maximum cooling = {cooling.max():.6f} °C"
    )

fig.suptitle(
    "Johannesburg: Spatial sensitivity of combined mitigation",
    fontsize=15
)

cbar = fig.colorbar(
    image,
    ax=axes,
    shrink=0.78,
    pad=0.02
)
cbar.set_label("Cooling relative to S0 (°C)")

for extension in ("png", "pdf"):
    output = outdir / (
        f"johannesburg_ndvi_sensitivity.{extension}"
    )
    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight"
    )
    print("Saved:", output)

plt.close(fig)
