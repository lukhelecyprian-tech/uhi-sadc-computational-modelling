import numpy as np
import rasterio

names = [
    "Manzini_LST",
    "Manzini_NDVI",
    "Manzini_Albedo",
    "Manzini_Mean_Reflectance",
]

masks = []
profile = None
reference = None

for name in names:
    with rasterio.open(f"data/{name}.tif") as src:
        geometry = (src.width, src.height, src.crs, src.transform)

        if reference is None:
            reference = geometry
            profile = src.profile.copy()
        elif geometry != reference:
            raise ValueError(f"{name} is not aligned.")

        arr = src.read(1, masked=True)
        valid = ~np.ma.getmaskarray(arr) & np.isfinite(arr.data)

        if name == "Manzini_LST":
            valid &= (arr.data > -50) & (arr.data < 80)
        elif name == "Manzini_NDVI":
            valid &= (arr.data >= -1) & (arr.data <= 1)
        else:
            valid &= (arr.data >= 0) & (arr.data <= 1)

        masks.append(valid)

common = np.logical_and.reduce(masks)

profile.update(
    dtype="uint8",
    count=1,
    nodata=None,
    compress="lzw",
)

output = "data/Manzini_four_variable_mask.tif"

with rasterio.open(output, "w", **profile) as dst:
    dst.write(common.astype(np.uint8), 1)

print("Saved:", output)
print("Valid pixels:", np.count_nonzero(common))
print("Valid area (km²):", np.count_nonzero(common) * 900 / 1e6)
