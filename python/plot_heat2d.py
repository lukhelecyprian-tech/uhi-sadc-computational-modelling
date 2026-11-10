
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Project directories
BASE = Path(__file__).resolve().parent.parent
RESULTS = BASE / "results"
FIGURES = BASE / "figures"

FIGURES.mkdir(exist_ok=True)

# Read temperature fields produced by Fortran
initial = np.loadtxt(RESULTS / "initial.csv", delimiter=",")
final = np.loadtxt(RESULTS / "final.csv", delimiter=",")

# Check dimensions and values
assert initial.shape == final.shape == (51, 51)
assert np.isfinite(initial).all()
assert np.isfinite(final).all()

difference = final - initial

print("Initial temperature range:",
      initial.min(), initial.max())

print("Final temperature range:",
      final.min(), final.max())

print("Final mean temperature:", final.mean())

print("Maximum temperature change:",
      np.max(np.abs(difference)))

# Create three-panel scientific figure
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

extent = [0, 50, 0, 50]

im0 = axes[0].imshow(
    initial, origin="lower", cmap="turbo",
    vmin=20, vmax=40, extent=extent
)
axes[0].set_title("Initial Temperature")

im1 = axes[1].imshow(
    final, origin="lower", cmap="turbo",
    vmin=20, vmax=40, extent=extent
)
axes[1].set_title("Final Temperature")

im2 = axes[2].imshow(
    difference, origin="lower",
    cmap="coolwarm", extent=extent,
    vmin=-2, vmax=2
)
axes[2].set_title("Temperature Change (Final - Initial)")

for ax in axes:
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")

fig.colorbar(im1, ax=axes[:2], label="Temperature (°C)",
             shrink=0.85, pad=0.02)

fig.colorbar(im2, ax=axes[2], label="Change (°C)",
             shrink=0.85)

fig.suptitle(
    "Two-Dimensional Heat Diffusion: Fortran Simulation",
    fontsize=14
)

fig.subplots_adjust(
    left=0.06, right=0.94, bottom=0.14,
    top=0.82, wspace=0.35
)

output = FIGURES / "heat2d_comparison.png"
fig.savefig(output, dpi=300, bbox_inches="tight")

print("Figure saved:", output)

plt.show()
