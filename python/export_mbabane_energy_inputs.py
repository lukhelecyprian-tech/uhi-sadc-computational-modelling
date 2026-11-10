
import numpy as np
import rasterio
from pathlib import Path

# ============================================================
# MSc COMPUTATIONAL PHYSICS - UHI MODELLING
#
# STEP 15D: PREPARE MBABANE SPATIAL ENERGY-BALANCE INPUTS
#
# Inputs:
#   LST, NDVI, Albedo, Mean Reflectance
#   Four-variable common observation mask
#
# Outputs:
#   Fortran-readable spatial data files
# ============================================================

data_dir = Path("data")

paths = {
    "LST": data_dir / "Mbabane_LST.tif",
    "NDVI": data_dir / "Mbabane_NDVI.tif",
    "Albedo": data_dir / "Mbabane_Albedo.tif",
    "Reflectance": data_dir / "Mbabane_Mean_Reflectance.tif"
}

mask_path = data_dir / "Mbabane_four_variable_mask.tif"

arrays = {}
valid_masks = {}

reference = None

# ------------------------------------------------------------
# 1. Read and verify satellite rasters
# ------------------------------------------------------------

for name, path in paths.items():

    with rasterio.open(path) as src:

        if reference is None:
            reference = (
                src.crs,
                src.transform,
                src.shape
            )
            ny, nx = src.shape
            dx, dy = src.res

        elif (
            src.crs,
            src.transform,
            src.shape
        ) != reference:
            raise ValueError(
                f"{name} does not match the reference grid."
            )

        band = src.read(1, masked=True)

        values = band.data.astype(np.float64)

        valid = (
            ~np.ma.getmaskarray(band)
            & np.isfinite(values)
        )

        arrays[name] = values
        valid_masks[name] = valid

# ------------------------------------------------------------
# 2. Read the verified common observation mask
# ------------------------------------------------------------

with rasterio.open(mask_path) as src:

    if (
        src.crs,
        src.transform,
        src.shape
    ) != reference:
        raise ValueError("Common mask grid mismatch.")

    common = src.read(1) == 1

# Confirm that the mask excludes all invalid values
for name in paths:

    if not np.all(valid_masks[name][common]):
        raise ValueError(
            f"Common mask contains invalid {name} pixels."
        )

# ------------------------------------------------------------
# 3. Convert observed LST from Celsius to Kelvin
# ------------------------------------------------------------

lst_kelvin = arrays["LST"] + 273.15

# ------------------------------------------------------------
# 4. Export Fortran-readable data
#
# Each output has ny rows and nx columns.
# Excluded pixels contain zero placeholders.
# The mask determines which pixels are active.
# ------------------------------------------------------------

outputs = {
    "temperature_K": lst_kelvin,
    "ndvi": arrays["NDVI"],
    "albedo": arrays["Albedo"],
    "reflectance": arrays["Reflectance"]
}

for name, values in outputs.items():

    prepared = np.where(common, values, 0.0)

    output_path = (
        data_dir / f"mbabane_energy_{name}.txt"
    )

    np.savetxt(
        output_path,
        prepared,
        fmt="%.10f"
    )

    print("Saved:", output_path)

# ------------------------------------------------------------
# 5. Export active-domain mask
# ------------------------------------------------------------

mask_output = data_dir / "mbabane_energy_mask.txt"

np.savetxt(
    mask_output,
    common.astype(np.int32),
    fmt="%d"
)

print("Saved:", mask_output)

# ------------------------------------------------------------
# 6. Export grid metadata
# ------------------------------------------------------------

grid_output = data_dir / "mbabane_energy_grid.txt"

with open(grid_output, "w") as f:

    f.write(f"{nx} {ny}\n")
    f.write(f"{dx:.6f} {dy:.6f}\n")

print("Saved:", grid_output)

# ------------------------------------------------------------
# 7. Final checks
# ------------------------------------------------------------

active = int(np.count_nonzero(common))

print()
print("=" * 60)
print("MBABANE SPATIAL ENERGY MODEL INPUT SUMMARY")
print("=" * 60)

print(f"Grid: {nx} columns x {ny} rows")
print(f"Pixel size: {dx:.1f} x {dy:.1f} metres")
print(f"Active pixels: {active}")
print(f"Excluded pixels: {common.size - active}")

print()
print("MEAN VALUES ON COMMON DOMAIN")
print(
    f"LST: {arrays['LST'][common].mean():.6f} C"
)
print(
    f"Temperature: {lst_kelvin[common].mean():.6f} K"
)
print(
    f"NDVI: {arrays['NDVI'][common].mean():.6f}"
)
print(
    f"Albedo: {arrays['Albedo'][common].mean():.6f}"
)
print(
    f"Reflectance: {arrays['Reflectance'][common].mean():.6f}"
)

print()
print("All four satellite variables exported successfully.")
print("Ready for spatial Fortran energy-balance modelling.")

