
import numpy as np
import rasterio
from pathlib import Path

data_dir = Path("data")

names = {
    "temperature": "Manzini_LST.tif",
    "ndvi": "Manzini_NDVI.tif",
    "albedo": "Manzini_Albedo.tif",
    "reflectance": "Manzini_Mean_Reflectance.tif",
}

arrays = {}

# Load the common mask and its geographical properties.
with rasterio.open(data_dir / "Manzini_four_variable_mask.tif") as src:
    mask = src.read(1).astype(bool)
    reference = (src.width, src.height, src.crs, src.transform)
    dx, dy = src.res

assert mask.sum() == 68158, "Unexpected number of valid pixels."

for variable, filename in names.items():
    with rasterio.open(data_dir / filename) as src:
        geometry = (src.width, src.height, src.crs, src.transform)

        if geometry != reference:
            raise ValueError(f"Raster alignment mismatch: {filename}")

        arr = src.read(1, masked=True)

        if np.any(np.ma.getmaskarray(arr)[mask]):
            raise ValueError(f"Missing data inside common mask: {filename}")

        values = arr.data.astype(np.float64)

        if not np.all(np.isfinite(values[mask])):
            raise ValueError(f"Non-finite values: {filename}")

        arrays[variable] = values

# Check units and physically meaningful ranges.
lst = arrays["temperature"]

if not np.all((-50 < lst[mask]) & (lst[mask] < 80)):
    raise ValueError("LST does not appear to be in Celsius.")

if not np.all((-1 <= arrays["ndvi"][mask]) &
              (arrays["ndvi"][mask] <= 1)):
    raise ValueError("NDVI outside [-1, 1].")

for name in ("albedo", "reflectance"):
    if not np.all((0 <= arrays[name][mask]) &
                  (arrays[name][mask] <= 1)):
        raise ValueError(f"{name} outside [0, 1].")

# Convert Celsius to kelvin on valid cells.
temperature_K = np.zeros(mask.shape, dtype=np.float64)
temperature_K[mask] = lst[mask] + 273.15

outputs = {
    "manzini_energy_temperature_K.txt": temperature_K,
    "manzini_energy_ndvi.txt": arrays["ndvi"],
    "manzini_energy_albedo.txt": arrays["albedo"],
    "manzini_energy_reflectance.txt": arrays["reflectance"],
    "manzini_energy_mask.txt": mask.astype(np.int32),
}

for filename, arr in outputs.items():
    # Set excluded cells to zero. Fortran must use the mask
    # and must never treat those zeros as physical observations.
    clean = np.where(mask, arr, 0)

    fmt = "%d" if filename.endswith("_mask.txt") else "%.12e"

    np.savetxt(data_dir / filename, clean, fmt=fmt)

# Grid metadata: columns, rows, then dx and dy in metres.
with open(data_dir / "manzini_energy_grid.txt", "w") as f:
    f.write(f"{mask.shape[1]} {mask.shape[0]}\n")
    f.write(f"{dx:.6f} {dy:.6f}\n")

print("MANZINI FORTRAN INPUT PREPARATION")
print("--------------------------------")
print("Columns:", mask.shape[1])
print("Rows:", mask.shape[0])
print("Valid pixels:", int(mask.sum()))
print("Grid spacing:", dx, dy, "metres")
print("Mean initial LST (C):", float(lst[mask].mean()))
print("Mean initial temperature (K):",
      float(temperature_K[mask].mean()))
print("Mean NDVI:", float(arrays["ndvi"][mask].mean()))
print("Mean albedo:", float(arrays["albedo"][mask].mean()))
print("Mean reflectance:",
      float(arrays["reflectance"][mask].mean()))
print("Fortran input files saved successfully.")
