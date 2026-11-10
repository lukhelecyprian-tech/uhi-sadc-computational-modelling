import numpy as np
import rasterio

files = {
    "LST": "data/Manzini_LST.tif",
    "NDVI": "data/Manzini_NDVI.tif",
    "Albedo": "data/Manzini_Albedo.tif",
    "Reflectance": "data/Manzini_Mean_Reflectance.tif",
}

for name, filename in files.items():

    with rasterio.open(filename) as src:
        data = src.read(1, masked=True)

        values = data.compressed()
        values = values[np.isfinite(values)]

        print("\n-------------------------")
        print(name)
        print("-------------------------")

        print("Valid pixels:", len(values))
        print("Minimum:", np.min(values))
        print("Maximum:", np.max(values))
        print("Mean:", np.mean(values))
        print("Standard deviation:", np.std(values))
