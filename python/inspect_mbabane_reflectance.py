
import numpy as np
import rasterio

reflectance_path = "data/Mbabane_Mean_Reflectance.tif"
reference_path = "data/Mbabane_LST.tif"

print("=" * 55)
print("MBABANE MEAN SURFACE REFLECTANCE INSPECTION")
print("=" * 55)

with rasterio.open(reflectance_path) as src:

    data = src.read(1, masked=True)
    values = data.compressed()
    values = values[np.isfinite(values)]

    print("\nRASTER INFORMATION")
    print("CRS:", src.crs)
    print("Width:", src.width)
    print("Height:", src.height)
    print("Resolution:", src.res)
    print("Number of bands:", src.count)
    print("Band descriptions:", src.descriptions)
    print("NoData:", src.nodata)

    print("\nREFLECTANCE STATISTICS")
    print("Valid pixels:", len(values))

    if len(values) > 0:
        print("Minimum:", np.min(values))
        print("Maximum:", np.max(values))
        print("Mean:", np.mean(values))
        print("Standard deviation:", np.std(values))
        print("Median:", np.median(values))

    print("\nGRID COMPATIBILITY")

    with rasterio.open(reference_path) as ref:

        print("CRS matches LST:", src.crs == ref.crs)
        print(
            "Dimensions match LST:",
            src.shape == ref.shape
        )
        print(
            "Transform matches LST:",
            src.transform == ref.transform
        )

print("\nInspection complete.")

