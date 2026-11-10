
import numpy as np
import rasterio
from scipy.ndimage import binary_fill_holes, label

path = "data/Mbabane_LST.tif"

with rasterio.open(path) as src:
    data = src.read(1, masked=True)
    clean = np.ma.masked_invalid(data)
    valid = ~np.ma.getmaskarray(clean)

# Fill enclosed holes in the binary mask only.
# This does NOT modify temperature values.
filled_mask = binary_fill_holes(valid)

internal_holes = filled_mask & ~valid

# Identify connected internal missing-data regions.
labels, number = label(internal_holes)

sizes = np.bincount(labels.ravel())[1:]
sizes = np.sort(sizes)[::-1]

print("=" * 45)
print("MBABANE LST MASK ANALYSIS")
print("=" * 45)

print("Total raster pixels:", valid.size)
print("Valid pixels:", np.count_nonzero(valid))
print("Excluded pixels:", np.count_nonzero(~valid))
print("Enclosed missing pixels:",
      np.count_nonzero(internal_holes))
print("Number of enclosed regions:", number)

print()
print("Largest enclosed regions (pixels):")
print(sizes[:10])

print()
print("Largest enclosed regions (m²):")
print(sizes[:10] * 30 * 30)

print("=" * 45)

