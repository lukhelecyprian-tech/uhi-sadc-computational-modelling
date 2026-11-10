
import numpy as np
import rasterio
from pathlib import Path

# ============================================================
# MSc COMPUTATIONAL PHYSICS - URBAN HEAT ISLAND THESIS
# MBABANE: FOUR-RASTER CONSISTENCY CHECK
#
# Variables:
#   1. Land Surface Temperature (LST)
#   2. Normalized Difference Vegetation Index (NDVI)
#   3. Surface Albedo
#   4. Mean Surface Reflectance
# ============================================================

data_dir = Path("data")

paths = {
    "LST": data_dir / "Mbabane_LST.tif",
    "NDVI": data_dir / "Mbabane_NDVI.tif",
    "Albedo": data_dir / "Mbabane_Albedo.tif",
    "Reflectance": data_dir / "Mbabane_Mean_Reflectance.tif"
}

arrays = {}
masks = {}

reference_crs = None
reference_transform = None
reference_shape = None
reference_profile = None

print("=" * 60)
print("MBABANE FOUR-RASTER CONSISTENCY CHECK")
print("=" * 60)

# ------------------------------------------------------------
# STEP 1: Read and verify raster alignment
# ------------------------------------------------------------

for name, path in paths.items():

    if not path.exists():
        raise FileNotFoundError(f"Missing raster: {path}")

    with rasterio.open(path) as src:

        data = src.read(1, masked=True)

        values = data.data.astype(np.float64)

        valid = (
            ~np.ma.getmaskarray(data)
        ) & np.isfinite(values)

        arrays[name] = values
        masks[name] = valid

        if reference_crs is None:

            reference_crs = src.crs
            reference_transform = src.transform
            reference_shape = src.shape
            reference_profile = src.profile.copy()

        else:

            if src.crs != reference_crs:
                raise ValueError(f"{name}: CRS mismatch")

            if src.shape != reference_shape:
                raise ValueError(f"{name}: dimensions mismatch")

            if src.transform != reference_transform:
                raise ValueError(f"{name}: grid alignment mismatch")

        print(f"\n{name}")
        print("  Valid pixels:", np.count_nonzero(valid))
        print("  CRS:", src.crs)
        print("  Dimensions:", src.width, "x", src.height)
        print("  Resolution:", src.res)

# ------------------------------------------------------------
# STEP 2: Identify pixels valid in ALL FOUR rasters
# ------------------------------------------------------------

common_mask = np.logical_and.reduce(
    [masks[name] for name in paths]
)

common_count = np.count_nonzero(common_mask)

pixel_area = abs(
    reference_transform.a * reference_transform.e
    - reference_transform.b * reference_transform.d
)

common_area_km2 = common_count * pixel_area / 1e6

print("\n" + "=" * 60)
print("COMMON OBSERVATION DOMAIN")
print("=" * 60)

print("Common valid pixels:", common_count)
print("Excluded pixels:", common_mask.size - common_count)
print(f"Common area: {common_area_km2:.4f} km2")

if common_count == 0:
    raise ValueError("No common valid pixels found.")

# ------------------------------------------------------------
# STEP 3: Examine negative reflectance values
# ------------------------------------------------------------

reflectance = arrays["Reflectance"]
reflectance_valid = masks["Reflectance"]

negative_mask = (
    reflectance_valid
    & (reflectance < 0)
)

negative_count = np.count_nonzero(negative_mask)
reflectance_count = np.count_nonzero(reflectance_valid)

negative_percentage = (
    100 * negative_count / reflectance_count
)

print("\n" + "=" * 60)
print("REFLECTANCE QUALITY CHECK")
print("=" * 60)

print("Negative reflectance pixels:", negative_count)
print(f"Percentage negative: {negative_percentage:.4f}%")

# ------------------------------------------------------------
# STEP 4: Statistics using identical observation pixels
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("STATISTICS ON COMMON VALID PIXELS")
print("=" * 60)

for name in paths:

    values = arrays[name][common_mask]

    print(f"\n{name}")
    print(f"  Minimum: {np.min(values):.6f}")
    print(f"  Maximum: {np.max(values):.6f}")
    print(f"  Mean: {np.mean(values):.6f}")
    print(f"  Median: {np.median(values):.6f}")
    print(f"  Standard deviation: {np.std(values):.6f}")

# ------------------------------------------------------------
# STEP 5: Export the common observation mask as GeoTIFF
# ------------------------------------------------------------

output_path = data_dir / "Mbabane_four_variable_mask.tif"

profile = reference_profile.copy()

profile.update(
    driver="GTiff",
    dtype="uint8",
    count=1,
    nodata=255,
    compress="lzw"
)

with rasterio.open(output_path, "w", **profile) as dst:

    dst.write(
        common_mask.astype(np.uint8),
        1
    )

    dst.set_band_description(
        1, "Four-variable common observation mask"
    )

print("\n" + "=" * 60)
print("OUTPUT")
print("=" * 60)

print("Saved common mask:", output_path)
print("Mask values: 1 = common valid, 0 = excluded")

print("\nFour-raster consistency check completed successfully.")

