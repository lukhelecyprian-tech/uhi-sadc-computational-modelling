
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

Path("figures").mkdir(exist_ok=True)

sigma = 5.670374419e-8
epsilon = 0.95

C_A = 2.0e6
h = 10.0

T_air = 293.15
T_sky = 283.15

S = 500.0
albedo = 0.12

T0 = 293.15
dt = 60.0
hours = 6.0
nt = int(hours * 3600 / dt)

latent_fluxes = [50.0, 75.0, 100.0]

def simulate(LE):
    T = np.zeros(nt + 1)
    T[0] = T0

    for n in range(nt):
        solar = (1.0 - albedo) * S

        longwave = epsilon * sigma * (
            T[n]**4 - T_sky**4
        )

        sensible = h * (T[n] - T_air)

        net_flux = solar - longwave - sensible - LE

        T[n + 1] = T[n] + dt * net_flux / C_A

    return T - 273.15

time_hours = np.arange(nt + 1) * dt / 3600

results = {}

plt.figure(figsize=(9, 5))

for LE in latent_fluxes:
    temperature = simulate(LE)
    results[LE] = temperature

    plt.plot(
        time_hours,
        temperature,
        linewidth=2,
        label=f"Latent heat flux = {LE:.0f} W/m²"
    )

print("=" * 55)
print("VEGETATION-RELATED LATENT HEAT SENSITIVITY")
print("=" * 55)

baseline_final = results[50.0][-1]

for LE in latent_fluxes:
    final = results[LE][-1]
    difference = final - baseline_final

    print()
    print("Latent heat flux:", LE, "W/m²")
    print("Final temperature:", final, "°C")
    print("Difference from baseline:", difference, "°C")

plt.xlabel("Simulation time (hours)")
plt.ylabel("Surface temperature (°C)")
plt.title(
    "Illustrative Surface Energy Balance:\n"
    "Sensitivity to Latent Heat Flux"
)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

output = "figures/vegetation_energy_balance_test.png"
plt.savefig(output, dpi=300)
plt.close()

print()
print("Figure saved:", output)

