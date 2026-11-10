
import numpy as np
import rasterio
from scipy.stats import pearsonr, spearmanr

paths = {
    "LST": "data/Mbabane_LST.tif",
    "NDVI": "data/Mbabane_NDVI.tif",
    "Albedo": "data/Mbabane_Albedo.tif"
}

arrays = {}
masks = {}

for name, path in paths.items():
    with rasterio.open(path) as src:
        data = src.read(1, masked=True)
        clean = np.ma.masked_invalid(data)

        arrays[name] = clean.data
        masks[name] = ~np.ma.getmaskarray(clean)

common = masks["LST"] & masks["NDVI"] & masks["Albedo"]

lst = arrays["LST"][common].astype(np.float64)
ndvi = arrays["NDVI"][common].astype(np.float64)
albedo = arrays["Albedo"][common].astype(np.float64)

print("=" * 55)
print("MBABANE LST-NDVI-ALBEDO CORRELATION ANALYSIS")
print("=" * 55)
print("Common valid pixels:", common.sum())

pairs = [
    ("LST", lst, "NDVI", ndvi),
    ("LST", lst, "Albedo", albedo),
    ("NDVI", ndvi, "Albedo", albedo)
]

for name1, x, name2, y in pairs:
    r, p = pearsonr(x, y)
    rho, ps = spearmanr(x, y)

    print()
    print(f"{name1} versus {name2}")
    print(f"Pearson r: {r:.6f}")
    print(f"Spearman rho: {rho:.6f}")
    print(f"Pearson p-value: {p:.6e}")
    print(f"Spearman p-value: {ps:.6e}")

print()
print("Mean LST (C):", lst.mean())
print("Mean NDVI:", ndvi.mean())
print("Mean albedo:", albedo.mean())
print("=" * 55)

