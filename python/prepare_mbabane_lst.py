
import numpy as np
import rasterio
from scipy.ndimage import (
    binary_fill_holes,
    distance_transform_edt
)
from pathlib import Path

input_path = "data/Mbabane_LST.tif"
output_path = "data/Mbabane_LST_prepared.tif"
mask_path = "data/Mbabane_LST_observed_mask.tif"

Path("data").mkdir(exist_ok=True)

with rasterio.open(input_path) as src:
    raw = src.read(1)
    original_mask = src.read_masks(1) > 0

    profile = src.profile.copy()

    # Original valid observations:
    # respect the GeoTIFF mask and exclude NaN/Inf
    valid = original_mask & np.isfinite(raw)

    # Find enclosed missing-data gaps
    filled_domain = binary_fill_holes(valid)
    holes = filled_domain & ~valid

    # Locate nearest observed pixel for each missing cell
    nearest_indices = distance_transform_edt(
        ~valid,
        return_distances=False,
        return_indices=True
    )

    # Preserve observed values exactly
    prepared = raw.astype(np.float32).copy()

    # Fill enclosed holes only
    prepared[holes] = raw[
        nearest_indices[0][holes],
        nearest_indices[1][holes]
    ]

    # Keep excluded exterior as NoData
    prepared[~filled_domain] = -9999.0

    # Save prepared raster
    profile.update(
        dtype="float32",
        count=1,
        nodata=-9999.0,
        compress="deflate"
    )

    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(prepared, 1)
        dst.set_band_description(1, "LST_C_prepared")

    # Save original observation mask
    mask_profile = profile.copy()
    mask_profile.update(
        dtype="uint8",
        nodata=255
    )

    with rasterio.open(mask_path, "w", **mask_profile) as dst:
        dst.write(valid.astype(np.uint8), 1)
        dst.set_band_description(1, "Observed_valid_mask")

print("=" * 45)
print("MBABANE LST MODEL PREPARATION")
print("=" * 45)
print("Original observed pixels:", valid.sum())
print("Enclosed pixels filled:", holes.sum())
print("Prepared domain pixels:", filled_domain.sum())
print("Excluded exterior pixels:", (~filled_domain).sum())
print()
print("Prepared raster:", output_path)
print("Observation mask:", mask_path)

