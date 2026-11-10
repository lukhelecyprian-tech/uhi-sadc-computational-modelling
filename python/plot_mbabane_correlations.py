
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from pathlib import Path

Path("figures").mkdir(exist_ok=True)

paths = {
    "LST": "data/Mbabane_LST.tif",
    "NDVI": "data/Mbabane_NDVI.tif",
    "Albedo": "data/Mbabane_Albedo.tif"
}

arrays = {}
masks = {}

for name, path in paths.items():
    with rasterio.open(path) as src:
        data = np.ma.masked_invalid(
            src.read(1, masked=True)
        )
        arrays[name] = data.data
        masks[name] = ~np.ma.getmaskarray(data)

common = masks["LST"] & masks["NDVI"] & masks["Albedo"]

lst = arrays["LST"][common].astype(float)
ndvi = arrays["NDVI"][common].astype(float)
albedo = arrays["Albedo"][common].astype(float)

pairs = [
    (ndvi, lst, "NDVI", "LST (°C)", "LST vs NDVI"),
    (albedo, lst, "Albedo", "LST (°C)", "LST vs Albedo"),
    (ndvi, albedo, "NDVI", "Albedo", "NDVI vs Albedo")
]

fig, axes = plt.subplots(
    1, 3,
    figsize=(16, 5),
    constrained_layout=True
)

for ax, (x, y, xlabel, ylabel, title) in zip(axes, pairs):

    # Density-aware plot using hexagonal bins
    hb = ax.hexbin(
        x, y,
        gridsize=55,
        mincnt=1,
        bins="log",
        cmap="viridis"
    )

    r, _ = pearsonr(x, y)

    # Simple least-squares trend line
    slope, intercept = np.polyfit(x, y, 1)
    xx = np.linspace(x.min(), x.max(), 200)

    ax.plot(
        xx,
        slope * xx + intercept,
        color="red",
        linewidth=1.5,
        label="Linear fit"
    )

    ax.set_title(f"{title}\nPearson r = {r:.3f}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc="best")

    fig.colorbar(
        hb, ax=ax,
        label="log10(pixel count)"
    )

fig.suptitle(
    "Mbabane: Pixelwise Relationships Between LST, "
    "NDVI and Albedo",
    fontsize=14,
    fontweight="bold"
)

output = "figures/mbabane_correlations.png"

fig.savefig(output, dpi=300)
plt.close(fig)

print("Figure saved:", output)
print("Common pixels:", len(lst))

