from pathlib import Path
import argparse
import numpy as np
import matplotlib.pyplot as plt
import rasterio

parser = argparse.ArgumentParser()
parser.add_argument("city", choices=["mbabane", "manzini", "johannesburg"])
args = parser.parse_args()

city = args.city
root = Path(".")
results = root / "results"
figures = root / "figures"

raster_path = root / "data" / f"{city.capitalize()}_LST.tif"

if city == "johannesburg":
    raster_path = Path(
        "/mnt/c/Users/Cyprian Lukhele/Desktop/"
        "MSc_UHI_Thesis/Johannesburg/Raster/"
        "Johannesburg_LST.tif"
    )

with rasterio.open(raster_path) as src:
    shape = (src.height, src.width)
    extent = (
        src.bounds.left, src.bounds.right,
        src.bounds.bottom, src.bounds.top
    )
    raster_mask = np.ma.getmaskarray(
        src.read(1, masked=True)
    )

if city == "johannesburg":
    with rasterio.open(
        "data/Johannesburg_four_variable_mask.tif"
    ) as mask_src:
        if mask_src.shape != shape:
            raise ValueError("Johannesburg mask shape mismatch")
        raster_mask = raster_mask | (mask_src.read(1) != 1)

def load_scenario(scenario):
    path = results / f"{city}_mitigation_{scenario}_K.txt"
    data = np.loadtxt(path)

    if data.size != shape[0] * shape[1]:
        raise ValueError(
            f"{path}: {data.size} values, "
            f"but raster requires {shape[0] * shape[1]}"
        )

    return data.reshape(shape)

baseline = load_scenario("S0")

scenarios = [
    ("S1", "Reflective surfaces"),
    ("S2", "Increased latent heat"),
    ("S3", "Combined mitigation"),
]

fig, axes = plt.subplots(
    1, 3, figsize=(17, 5.5), constrained_layout=True
)

for ax, (scenario, label) in zip(axes, scenarios):
    temperature = load_scenario(scenario)
    cooling = baseline - temperature

    cooling = np.ma.masked_where(
        raster_mask | ~np.isfinite(cooling),
        cooling
    )

    image = ax.imshow(
        cooling,
        origin="upper",
        extent=extent,
        cmap="viridis",
        vmin=0,
        vmax=0.8,
        interpolation="nearest"
    )

    ax.set_title(f"{scenario}: {label}", fontsize=12)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_aspect("equal")

    print(
        scenario,
        "Mean cooling:",
        float(cooling.mean()),
        "Maximum cooling:",
        float(cooling.max())
    )

fig.suptitle(
    f"{city.capitalize()}: Illustrative cooling after 6 hours",
    fontsize=16
)

cbar = fig.colorbar(
    image, ax=axes, shrink=0.8, pad=0.02
)
cbar.set_label("Cooling relative to baseline (°C)")

output = figures / f"{city}_mitigation_cooling_maps.png"
fig.savefig(output, dpi=300, bbox_inches="tight")

pdf = output.with_suffix(".pdf")
fig.savefig(pdf, bbox_inches="tight")

plt.close(fig)

print("Saved:", output)
print("Saved:", pdf)
