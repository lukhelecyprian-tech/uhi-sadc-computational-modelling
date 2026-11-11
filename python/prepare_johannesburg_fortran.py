import numpy as np
import rasterio
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

data_dir = PROJECT_ROOT / "data"
source_dir = data_dir / "raw" / "johannesburg"

names = {
    "temperature": "Johannesburg_LST.tif",
    "ndvi": "Johannesburg_NDVI.tif",
    "albedo": "Johannesburg_Albedo.tif",
    "reflectance": "Johannesburg_Mean_Reflectance.tif",
}

data_dir.mkdir(exist_ok=True)
arrays = {}
valid_masks = []

reference = None
dx = dy = None

for variable, filename in names.items():
    with rasterio.open(source_dir / filename) as src:
        geometry = (
            src.width, src.height, src.crs, src.transform
        )

        if reference is None:
            reference = geometry
            dx, dy = src.res
        elif geometry != reference:
            raise ValueError(
                f"Raster alignment mismatch: {filename}"
            )

        arr = src.read(1, masked=True)
        values = arr.data.astype(np.float64)

        valid = (
            ~np.ma.getmaskarray(arr)
            & np.isfinite(values)
        )

        arrays[variable] = values
        valid_masks.append(valid)

# Same common-mask principle as Manzini.
mask = np.logical_and.reduce(valid_masks)

print("Initial common valid pixels:", int(mask.sum()))

# Apply physical validity checks consistently.
lst = arrays["temperature"]
ndvi = arrays["ndvi"]
albedo = arrays["albedo"]
reflectance = arrays["reflectance"]

physical = (
    (lst > -50) & (lst < 80)
    & (ndvi >= -1) & (ndvi <= 1)
    & (albedo >= 0) & (albedo <= 1)
    & (reflectance >= 0) & (reflectance <= 1)
)

mask &= physical

print("Physically valid pixels:", int(mask.sum()))
print("Excluded by physical checks:",
      int(np.logical_and.reduce(valid_masks).sum() - mask.sum()))

if not mask.any():
    raise ValueError("No valid modelling pixels.")

# Celsius to Kelvin.
temperature_K = np.zeros(mask.shape, dtype=np.float64)
temperature_K[mask] = lst[mask] + 273.15

outputs = {
    "johannesburg_energy_temperature_K.txt": temperature_K,
    "johannesburg_energy_ndvi.txt": ndvi,
    "johannesburg_energy_albedo.txt": albedo,
    "johannesburg_energy_reflectance.txt": reflectance,
    "johannesburg_energy_mask.txt": mask.astype(np.int32),
}

for filename, arr in outputs.items():
    clean = np.where(mask, arr, 0)
    fmt = "%d" if filename.endswith("_mask.txt") else "%.12e"
    np.savetxt(data_dir / filename, clean, fmt=fmt)

with open(data_dir / "johannesburg_energy_grid.txt", "w") as f:
    f.write(f"{mask.shape[1]} {mask.shape[0]}\n")
    f.write(f"{dx:.6f} {dy:.6f}\n")

# Preserve geographical reference for mapping.
with rasterio.open(source_dir / names["temperature"]) as src:
    profile = src.profile.copy()
    profile.update(
        dtype="uint8",
        nodata=0,
        count=1,
        compress="lzw"
    )

    with rasterio.open(
        data_dir / "Johannesburg_four_variable_mask.tif",
        "w",
        **profile
    ) as dst:
        dst.write(mask.astype(np.uint8), 1)

print("\nJOHANNESBURG FORTRAN INPUT PREPARATION")
print("-------------------------------------")
print("Columns:", mask.shape[1])
print("Rows:", mask.shape[0])
print("Valid pixels:", int(mask.sum()))
print("Grid spacing:", dx, dy, "metres")
print("Mean initial LST (C):", float(lst[mask].mean()))
print("Mean initial temperature (K):",
      float(temperature_K[mask].mean()))
print("Mean NDVI:", float(ndvi[mask].mean()))
print("Mean albedo:", float(albedo[mask].mean()))
print("Mean reflectance:", float(reflectance[mask].mean()))
print("Fortran input files saved successfully.")
