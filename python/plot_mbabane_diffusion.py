
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from pathlib import Path

Path("results").mkdir(exist_ok=True)
Path("figures").mkdir(exist_ok=True)

input_raster = "data/Mbabane_LST_prepared.tif"
fortran_output = "results/mbabane_diffusion_final.txt"

with rasterio.open(input_raster) as src:
    initial = src.read(1, masked=True)
    profile = src.profile.copy()
    bounds = src.bounds

    ny, nx = initial.shape
    domain = ~np.ma.getmaskarray(initial)

    final = np.loadtxt(fortran_output)

    if final.shape != (ny, nx):
        raise ValueError(
            f"Shape mismatch: {final.shape} != {(ny, nx)}"
        )

    if not np.isfinite(final[domain]).all():
        raise ValueError("Non-finite simulated temperatures.")

    final_masked = np.ma.array(
        final, mask=~domain
    )

    change = final_masked - initial

    nodata = -9999.0

    profile.update(
        dtype="float32",
        count=1,
        nodata=nodata,
        compress="deflate"
    )

    with rasterio.open(
        "results/Mbabane_LST_diffusion_final.tif",
        "w", **profile
    ) as dst:
        dst.write(
            final_masked.filled(nodata).astype("float32"),
            1
        )
        dst.set_band_description(1, "Diffusion_LST_C")

    with rasterio.open(
        "results/Mbabane_LST_diffusion_change.tif",
        "w", **profile
    ) as dst:
        dst.write(
            change.filled(nodata).astype("float32"),
            1
        )
        dst.set_band_description(1, "Diffusion_Delta_LST_C")

    extent = [
        bounds.left / 1000,
        bounds.right / 1000,
        bounds.bottom / 1000,
        bounds.top / 1000
    ]

vmin = float(initial.min())
vmax = float(initial.max())

fig, axes = plt.subplots(
    1, 3,
    figsize=(17, 6),
    constrained_layout=True
)

im0 = axes[0].imshow(
    initial,
    cmap="inferno",
    origin="upper",
    extent=extent,
    vmin=vmin,
    vmax=vmax
)

axes[0].set_title("Initial Landsat LST")

im1 = axes[1].imshow(
    final_masked,
    cmap="inferno",
    origin="upper",
    extent=extent,
    vmin=vmin,
    vmax=vmax
)

axes[1].set_title("Final Diffusion Temperature")

max_abs = float(np.max(np.abs(change.compressed())))

im2 = axes[2].imshow(
    change,
    cmap="RdBu_r",
    origin="upper",
    extent=extent,
    vmin=-max_abs,
    vmax=max_abs
)

axes[2].set_title("Temperature Change")

for ax in axes:
    ax.set_xlabel("UTM Easting (km)")
    ax.set_ylabel("UTM Northing (km)")
    ax.set_aspect("equal")

fig.colorbar(
    im0, ax=axes[:2],
    label="Temperature (°C)",
    shrink=0.75
)

fig.colorbar(
    im2, ax=axes[2],
    label="Temperature change (°C)",
    shrink=0.75
)

fig.suptitle(
    "Mbabane: Landsat-Initialised 2D Heat Diffusion\n"
    "Zero-flux boundaries | 30 m grid | 200,000 s",
    fontsize=14,
    fontweight="bold"
)

figure_path = "figures/mbabane_diffusion_comparison.png"

fig.savefig(figure_path, dpi=300)
plt.close(fig)

print("=" * 50)
print("MBABANE DIFFUSION POST-PROCESSING")
print("=" * 50)
print("Initial mean (C):", initial.mean())
print("Final mean (C):", final_masked.mean())
print("Maximum absolute change (C):", max_abs)
print()
print("Created:")
print("results/Mbabane_LST_diffusion_final.tif")
print("results/Mbabane_LST_diffusion_change.tif")
print(figure_path)

