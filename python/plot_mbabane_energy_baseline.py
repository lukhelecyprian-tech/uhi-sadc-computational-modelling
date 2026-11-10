
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# MSc COMPUTATIONAL PHYSICS - UHI THESIS
# MBABANE SPATIAL ENERGY-BALANCE POST-PROCESSING
# ============================================================

Path("results").mkdir(exist_ok=True)
Path("figures").mkdir(exist_ok=True)

reference_path = "data/Mbabane_LST.tif"
mask_path = "data/Mbabane_four_variable_mask.tif"

model_path = "results/mbabane_energy_baseline_K.txt"

# ------------------------------------------------------------
# 1. Read original Landsat temperature
# ------------------------------------------------------------

with rasterio.open(reference_path) as src:
    observed = src.read(1).astype(np.float64)
    profile = src.profile.copy()

    ny, nx = src.shape

# ------------------------------------------------------------
# 2. Read common observation mask
# ------------------------------------------------------------

with rasterio.open(mask_path) as src:
    common_mask = src.read(1) == 1

# ------------------------------------------------------------
# 3. Read Fortran simulation output
# ------------------------------------------------------------

model_kelvin = np.loadtxt(model_path)

if model_kelvin.shape != (ny, nx):
    raise ValueError(
        "Fortran output dimensions do not match Landsat."
    )

model_celsius = model_kelvin - 273.15

# ------------------------------------------------------------
# 4. Apply the same mask to all outputs
# ------------------------------------------------------------

initial = np.where(common_mask, observed, np.nan)

final = np.where(
    common_mask,
    model_celsius,
    np.nan
)

change = final - initial

# ------------------------------------------------------------
# 5. Export georeferenced GeoTIFF files
# ------------------------------------------------------------

profile.update(
    driver="GTiff",
    dtype="float32",
    count=1,
    nodata=-9999.0,
    compress="lzw"
)

outputs = {
    "results/Mbabane_Energy_Initial_LST.tif": initial,
    "results/Mbabane_Energy_Baseline_Final.tif": final,
    "results/Mbabane_Energy_Baseline_Change.tif": change
}

for filename, array in outputs.items():

    export = np.where(
        np.isfinite(array),
        array,
        -9999.0
    ).astype(np.float32)

    with rasterio.open(filename, "w", **profile) as dst:
        dst.write(export, 1)

    print("Saved:", filename)

# ------------------------------------------------------------
# 6. Create scientific comparison figure
# ------------------------------------------------------------

valid_initial = initial[common_mask]
valid_final = final[common_mask]
valid_change = change[common_mask]

vmin = min(valid_initial.min(), valid_final.min())
vmax = max(valid_initial.max(), valid_final.max())

fig, axes = plt.subplots(
    1, 3,
    figsize=(16, 6),
    constrained_layout=True
)

im0 = axes[0].imshow(
    initial,
    cmap="turbo",
    vmin=vmin,
    vmax=vmax
)

axes[0].set_title("Observed Initial LST (°C)")
axes[0].set_axis_off()

im1 = axes[1].imshow(
    final,
    cmap="turbo",
    vmin=vmin,
    vmax=vmax
)

axes[1].set_title("Modelled Baseline Final (°C)")
axes[1].set_axis_off()

fig.colorbar(
    im1,
    ax=axes[:2],
    shrink=0.7,
    label="Temperature (°C)"
)

max_abs = np.max(np.abs(valid_change))

im2 = axes[2].imshow(
    change,
    cmap="RdBu_r",
    vmin=-max_abs,
    vmax=max_abs
)

axes[2].set_title("Modelled Temperature Change (°C)")
axes[2].set_axis_off()

fig.colorbar(
    im2,
    ax=axes[2],
    shrink=0.7,
    label="Temperature change (°C)"
)

fig.suptitle(
    "Mbabane: Illustrative Six-Hour Spatial Energy Balance",
    fontsize=14
)

figure_path = "figures/mbabane_energy_baseline.png"

plt.savefig(
    figure_path,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ------------------------------------------------------------
# 7. Verify numerical results
# ------------------------------------------------------------

print()
print("=" * 55)
print("MBABANE SPATIAL ENERGY BALANCE RESULTS")
print("=" * 55)

print("Active pixels:", common_mask.sum())

print(
    "Initial mean (C):",
    valid_initial.mean()
)

print(
    "Final mean (C):",
    valid_final.mean()
)

print(
    "Mean change (C):",
    valid_change.mean()
)

print(
    "Maximum absolute change (C):",
    np.abs(valid_change).max()
)

print()
print("Figure saved:", figure_path)
print("Post-processing complete.")

