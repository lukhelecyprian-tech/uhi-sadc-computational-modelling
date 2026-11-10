from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import rasterio

figdir = Path("figures/baseline")
figdir.mkdir(parents=True, exist_ok=True)

cities = ["mbabane", "manzini", "johannesburg"]
datasets = []

for city in cities:
    raster_path = Path("data") / f"{city.capitalize()}_LST.tif"
    model_path = Path("results") / f"{city}_mitigation_S0_K.txt"

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
        mask = np.ma.getmaskarray(src.read(1, masked=True))

    if city == "johannesburg":
        with rasterio.open(
            "data/Johannesburg_four_variable_mask.tif"
        ) as mask_src:
            if mask_src.shape != shape:
                raise ValueError("Johannesburg mask shape mismatch")
            mask = mask | (mask_src.read(1) != 1)

    temperature = np.loadtxt(model_path)

    if temperature.size != shape[0] * shape[1]:
        raise ValueError(f"Grid mismatch: {city}")

    temperature = temperature.reshape(shape)

    # Exclude missing, zero and invalid Kelvin values.
    invalid = (
        mask
        | ~np.isfinite(temperature)
        | (temperature <= 0)
    )

    temperature = np.ma.masked_where(
        invalid,
        temperature - 273.15
    )

    datasets.append((city, temperature, extent))

vmin = min(float(t.min()) for _, t, _ in datasets)
vmax = max(float(t.max()) for _, t, _ in datasets)

for city, temperature, extent in datasets:
    fig, ax = plt.subplots(figsize=(8, 7))

    image = ax.imshow(
        temperature,
        extent=extent,
        origin="upper",
        cmap="inferno",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest"
    )

    ax.set_title(
        f"{city.capitalize()}: Baseline S0 temperature after 6 hours"
    )
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_aspect("equal")

    cbar = fig.colorbar(image, ax=ax, shrink=0.8)
    cbar.set_label("Modelled temperature (°C)")

    for extension in ["png", "pdf"]:
        output = figdir / f"{city}_baseline_S0.{extension}"
        fig.savefig(output, dpi=300, bbox_inches="tight")
        print("Saved:", output)

    print(
        city.capitalize(),
        "Mean baseline temperature:",
        round(float(temperature.mean()), 4),
        "°C"
    )

    plt.close(fig)
