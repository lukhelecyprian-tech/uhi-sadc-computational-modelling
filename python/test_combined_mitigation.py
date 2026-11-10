
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
T0 = 293.15

dt = 60.0
hours = 6.0
nt = int(hours * 3600 / dt)

scenarios = {
    "Baseline": (0.12, 50.0),
    "Higher albedo": (0.22, 50.0),
    "Higher latent cooling": (0.12, 100.0),
    "Combined mitigation": (0.22, 100.0)
}

def simulate(albedo, latent_flux):

    T = np.zeros(nt + 1)
    T[0] = T0

    for n in range(nt):

        absorbed_solar = (1.0 - albedo) * S

        longwave_loss = epsilon * sigma * (
            T[n]**4 - T_sky**4
        )

        sensible_loss = h * (T[n] - T_air)

        net_flux = (
            absorbed_solar
            - longwave_loss
            - sensible_loss
            - latent_flux
        )

        T[n + 1] = T[n] + dt * net_flux / C_A

    return T - 273.15

time_hours = np.arange(nt + 1) * dt / 3600

results = {}

plt.figure(figsize=(10, 6))

for name, (albedo, latent_flux) in scenarios.items():

    temperatures = simulate(albedo, latent_flux)
    results[name] = temperatures

    plt.plot(
        time_hours,
        temperatures,
        linewidth=2,
        label=name
    )

baseline_final = results["Baseline"][-1]

print("=" * 55)
print("COMBINED UHI MITIGATION SENSITIVITY TEST")
print("=" * 55)

for name, temperatures in results.items():

    final = temperatures[-1]
    difference = final - baseline_final

    print()
    print("Scenario:", name)
    print(f"Final temperature: {final:.6f} C")
    print(f"Change from baseline: {difference:.6f} C")

plt.xlabel("Simulation time (hours)")
plt.ylabel("Surface temperature (°C)")
plt.title(
    "Illustrative Surface Energy Balance:\n"
    "Albedo and Latent Cooling Scenarios"
)

plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

output = "figures/combined_mitigation_test.png"

plt.savefig(output, dpi=300)
plt.close()

print()
print("Figure saved:", output)

