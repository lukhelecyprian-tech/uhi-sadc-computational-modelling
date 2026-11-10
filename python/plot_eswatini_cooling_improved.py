from pathlib import Path

import numpy as np
import rasterio
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

ROOT = Path(
    "/mnt/c/Users/Cyprian Lukhele/Desktop/MSc_UHI_Thesis"
)

CITIES = {
    "Mbabane": ROOT / "Mbabane/Figures/Mbabane.gpkg",
    "Manzini": ROOT / "Manzini/Raster/manzini_boundary.gpkg",
}

OUT = Path("figures/improved")
OUT.mkdir(parents=True, exist_ok=True)

SCENARIOS = {
    "S1": "Reflective surfaces",
    "S2": "Enhanced latent heat",
    "S3": "Combined mitigation",
}

for city, boundary_path in CITIES.items():

    city_lower = city.lower()

    with rasterio.open(f"data/{city}_LST.tif") as src:
        shape = src.shape
        crs = src.crs
        transform = src.transform
        bounds = src.bounds

    with rasterio.open(
        f"data/{city}_four_variable_mask.tif"
    ) as src:
        assert src.shape == shape
        assert src.crs == crs
        assert src.transform == transform
        mask = src.read(1) == 1

    boundary = gpd.read_file(boundary_path).to_crs(crs)

    def load_temperature(scenario):
        path = (
            Path("results")
            / f"{city_lower}_mitigation_{scenario}_K.txt"
        )
        values = np.loadtxt(path)

        if values.size != np.prod(shape):
            raise ValueError(
                f"{path}: incorrect number of pixels"
            )

        return values.reshape(shape)

    baseline = load_temperature("S0")

    extent = [
        bounds.left,
        bounds.right,
        bounds.bottom,
        bounds.top,
    ]

    fig, axes = plt.subplots(
        1, 3, figsize=(17, 6),
        sharex=True, sharey=True,
        constrained_layout=True,
    )

    for ax, (scenario, description) in zip(
        axes, SCENARIOS.items()
    ):
        simulated = load_temperature(scenario)

        cooling = baseline - simulated
        cooling = np.where(mask, cooling, np.nan)

        im = ax.imshow(
            cooling,
            cmap="viridis",
            vmin=0,
            vmax=0.8,
            extent=extent,
            origin="upper",
            interpolation="nearest",
        )

        boundary.boundary.plot(
            ax=ax,
            color="black",
            linewidth=0.85,
            zorder=3,
        )

        ax.set_title(
            f"{scenario} – {description}",
            fontsize=11,
            fontweight="bold",
        )

        ax.set_xlabel("Easting (m)", fontsize=10)
        ax.set_ylabel("Northing (m)", fontsize=10)

        ax.set_aspect("equal")

        for axis in [ax.xaxis, ax.yaxis]:
            formatter = ScalarFormatter(
                useOffset=False
            )
            formatter.set_scientific(False)
            axis.set_major_formatter(formatter)

        ax.tick_params(
            axis="both",
            labelsize=8,
            rotation=0,
        )

    fig.suptitle(
        f"{city}: Simulated surface cooling after six hours",
        fontsize=16,
        fontweight="bold",
        y=1.08,
    )

    cbar = fig.colorbar(
        im,
        ax=axes,
        orientation="horizontal",
        fraction=0.045,
        pad=0.10,
        shrink=0.65,
    )
    cbar.set_label(
        "Cooling relative to baseline S0 (°C)",
        fontsize=11,
    )

    png_path = OUT / f"{city_lower}_cooling_improved.png"
    pdf_path = OUT / f"{city_lower}_cooling_improved.pdf"

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")

    plt.close(fig)

    print(f"\n{city}")
    print("Active model pixels:", np.count_nonzero(mask))
    print("PNG:", png_path)
    print("PDF:", pdf_path)

print("\nBoth improved figures generated successfully.")
