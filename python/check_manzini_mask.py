import numpy as np
import rasterio

files = {
    "LST": "data/Manzini_LST.tif",
    "NDVI": "data/Manzini_NDVI.tif",
    "Albedo": "data/Manzini_Albedo.tif",
    "Reflectance": "data/Manzini_Mean_Reflectance.tif",
}

masks = []
reference = None

for name, filename in files.items():

    with rasterio.open(filename) as src:
        geometry = (
            src.width,
            src.height,
            src.crs,
            src.transform,
        )

        if reference is None:
            reference = geometry
        elif geometry != reference:
            raise ValueError(
                f"{name} is not aligned with the LST raster."
            )

        data = src.read(1, masked=True)

        # Valid according to the raster mask and finite values
        valid = (
            ~np.ma.getmaskarray(data)
            & np.isfinite(data.data)
        )

        # Basic physical plausibility checks
        if name == "LST":
            valid &= (data.data > -50) & (data.data < 80)

        elif name == "NDVI":
            valid &= (data.data >= -1) & (data.data <= 1)

        elif name in ("Albedo", "Reflectance"):
            valid &= (data.data >= 0) & (data.data <= 1)

        masks.append(valid)

        print(f"{name}: {np.count_nonzero(valid):,} usable pixels")

common = np.logical_and.reduce(masks)

pixel_area = 30 * 30
area_km2 = np.count_nonzero(common) * pixel_area / 1_000_000

print("\nMANZINI COMMON MASK")
print("-------------------")
print("Common valid pixels:", np.count_nonzero(common))
print("Common valid area (km²):", area_km2)
print("Total grid pixels:", common.size)
print("Valid fraction (%):", 100 * common.mean())
