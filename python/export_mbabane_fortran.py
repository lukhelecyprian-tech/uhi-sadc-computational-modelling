
import numpy as np
import rasterio
from pathlib import Path

Path("data").mkdir(exist_ok=True)

input_path = "data/Mbabane_LST_prepared.tif"

with rasterio.open(input_path) as src:
    raster = src.read(1, masked=True)

    # True inside the prepared modelling domain
    domain = ~np.ma.getmaskarray(raster)

    # Use zero only as a storage placeholder outside
    # the domain. Fortran must never treat it as
    # a physical temperature.
    temperature = raster.filled(0.0).astype(np.float64)

    ny, nx = temperature.shape
    dx, dy = src.res

# Write metadata for the Fortran solver
with open("data/mbabane_grid.txt", "w") as f:
    f.write(f"{nx} {ny}\n")
    f.write(f"{dx:.8f} {dy:.8f}\n")

# Python rows correspond to raster rows.
# Fortran will read these values in row order.
np.savetxt(
    "data/mbabane_temperature.txt",
    temperature,
    fmt="%.8f"
)

np.savetxt(
    "data/mbabane_domain_mask.txt",
    domain.astype(np.int32),
    fmt="%d"
)

print("=" * 45)
print("MBABANE PYTHON-TO-FORTRAN EXPORT")
print("=" * 45)
print("Grid:", nx, "x", ny)
print("Resolution:", dx, "x", dy, "m")
print("Active domain cells:", domain.sum())
print("Excluded cells:", (~domain).sum())
print()
print("Files created:")
print("data/mbabane_grid.txt")
print("data/mbabane_temperature.txt")
print("data/mbabane_domain_mask.txt")

