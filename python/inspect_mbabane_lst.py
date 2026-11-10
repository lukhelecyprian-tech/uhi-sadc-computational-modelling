
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from pathlib import Path

raster_path = "data/Mbabane_LST.tif"
output_path = "figures/mbabane_lst_quality.png"

Path("figures").mkdir(exist_ok=True)

with rasterio.open(raster_path) as src:
    raw = src.read(1)
    masked = src.read(1, masked=True)
    lst = np.ma.masked_invalid(masked)

    valid_mask = ~np.ma.getmaskarray(lst)

    bounds = src.bounds
    extent = [
        bounds.left / 1000,
        bounds.right / 1000,
        bounds.bottom / 1000,
        bounds.top / 1000
    ]

fig, axes = plt.subplots(
    1, 2, figsize=(13, 6),
    constrained_layout=True
)

im = axes[0].imshow(
    lst,
    cmap="inferno",
    extent=extent,
    origin="upper"
)

axes[0].set_title("Mbabane Landsat LST")
axes[0].set_xlabel("UTM Easting (km)")
axes[0].set_ylabel("UTM Northing (km)")

fig.colorbar(
    im, ax=axes[0],
    label="Land Surface Temperature (°C)",
    shrink=0.8
)

axes[1].imshow(
    valid_mask,
    cmap="gray_r",
    extent=extent,
    origin="upper",
    vmin=0,
    vmax=1
)

axes[1].set_title("Valid-Data Mask")
axes[1].set_xlabel("UTM Easting (km)")
axes[1].set_ylabel("UTM Northing (km)")

fig.suptitle(
    "Mbabane LST: Observational Data Quality",
    fontsize=14,
    fontweight="bold"
)

fig.savefig(output_path, dpi=300)
plt.close(fig)

print("Figure saved to:", output_path)
print("Valid finite pixels:", valid_mask.sum())
print("Excluded pixels:", valid_mask.size - valid_mask.sum())

