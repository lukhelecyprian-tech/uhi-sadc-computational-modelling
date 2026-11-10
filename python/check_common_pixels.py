
import numpy as np
import rasterio

files = {
    "LST": "data/Mbabane_LST.tif",
    "NDVI": "data/Mbabane_NDVI.tif",
    "Albedo": "data/Mbabane_Albedo.tif"
}

masks = {}

for name, path in files.items():
    with rasterio.open(path) as src:
        data = src.read(1, masked=True)
        clean = np.ma.masked_invalid(data)
        masks[name] = ~np.ma.getmaskarray(clean)

common = masks["LST"] & masks["NDVI"] & masks["Albedo"]

print("=" * 50)
print("MBABANE COMMON VALID-PIXEL ANALYSIS")
print("=" * 50)

for name in files:
    print(f"{name} valid pixels:", masks[name].sum())

print()
print("Common valid pixels:", common.sum())

for name in files:
    missing_from_common = masks[name] & ~common
    print(
        f"{name} pixels excluded from common mask:",
        missing_from_common.sum()
    )

print()
print("Common valid area (km²):",
      common.sum() * 30 * 30 / 1e6)

# Save a common-validity mask GeoTIFF.
# 1 = valid in all three rasters
# 0 = not valid in all three rasters
with rasterio.open(files["LST"]) as src:
    profile = src.profile.copy()

profile.update(
    dtype="uint8",
    count=1,
    nodata=255,
    compress="deflate"
)

output = "data/Mbabane_common_observation_mask.tif"

with rasterio.open(output, "w", **profile) as dst:
    dst.write(common.astype(np.uint8), 1)
    dst.set_band_description(1, "Common_valid_mask")

print()
print("Saved:", output)

