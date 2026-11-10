
import numpy as np
import rasterio

files = {
    "LST": "data/Mbabane_LST.tif",
    "NDVI": "data/Mbabane_NDVI.tif",
    "Albedo": "data/Mbabane_Albedo.tif"
}

metadata = {}

print("=" * 55)
print("MBABANE MULTI-RASTER COMPATIBILITY CHECK")
print("=" * 55)

for name, path in files.items():

    with rasterio.open(path) as src:

        data = src.read(1, masked=True)
        data = np.ma.masked_invalid(data)
        values = data.compressed()

        metadata[name] = {
            "crs": src.crs,
            "shape": src.shape,
            "transform": src.transform,
            "resolution": src.res
        }

        print()
        print("DATASET:", name)
        print("CRS:", src.crs)
        print("Shape:", src.shape)
        print("Resolution:", src.res)
        print("NoData:", src.nodata)
        print("Band description:", src.descriptions[0])
        print("Valid finite pixels:", values.size)

        if values.size > 0:
            print("Minimum:", np.min(values))
            print("Maximum:", np.max(values))
            print("Mean:", np.mean(values))

print()
print("=" * 55)
print("SPATIAL ALIGNMENT")
print("=" * 55)

reference = metadata["LST"]

for name in ["NDVI", "Albedo"]:

    item = metadata[name]

    crs_match = item["crs"] == reference["crs"]
    shape_match = item["shape"] == reference["shape"]

    transform_match = np.allclose(
        tuple(item["transform"]),
        tuple(reference["transform"]),
        rtol=0,
        atol=1e-8
    )

    print()
    print(name)
    print("CRS matches LST:", crs_match)
    print("Shape matches LST:", shape_match)
    print("Grid transform matches LST:", transform_match)

    if crs_match and shape_match and transform_match:
        print("RESULT: ALIGNED")
    else:
        print("RESULT: REQUIRES ALIGNMENT")

