
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

Path("figures").mkdir(exist_ok=True)

# -------------------------------------
# Illustrative physical parameters
# -------------------------------------

sigma = 5.670374419e-8  # W m^-2 K^-4
epsilon = 0.95

C_A = 2.0e6  # J m^-2 K^-1
h = 10.0     # W m^-2 K^-1

T_air = 293.15  # K (20 C)
T_sky = 283.15  # K (10 C)

S = 500.0       # W m^-2
LE = 50.0       # W m^-2

# Compare baseline and higher albedo
albedo_baseline = 0.12
albedo_scenario = 0.22

# Initial temperature: 20 C
T0 = 293.15

dt = 60.0        # seconds
hours = 6.0
nt = int(hours * 3600 / dt)

def simulate(albedo):

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
            - LE
        )

        T[n + 1] = T[n] + dt * net_flux / C_A

    return T - 273.15

baseline = simulate(albedo_baseline)
scenario = simulate(albedo_scenario)

time_hours = np.arange(nt + 1) * dt / 3600

print("=" * 50)
print("SURFACE ENERGY BALANCE: ALBEDO EXPERIMENT")
print("=" * 50)

print("Initial temperature: 20 C")
print("Baseline albedo:", albedo_baseline)
print("Scenario albedo:", albedo_scenario)

print()
print("Final baseline temperature:", baseline[-1], "C")
print("Final scenario temperature:", scenario[-1], "C")
print(
    "Scenario minus baseline:",
    scenario[-1] - baseline[-1],
    "C"
)

plt.figure(figsize=(9, 5))

plt.plot(
    time_hours, baseline,
    label="Baseline albedo = 0.12",
    linewidth=2
)

plt.plot(
    time_hours, scenario,
    label="Higher albedo = 0.22",
    linewidth=2
)

plt.xlabel("Simulation time (hours)")
plt.ylabel("Surface temperature (°C)")
plt.title(
    "Illustrative Surface Energy Balance:\n"
    "Effect of Increasing Albedo"
)

plt.legend()
plt.grid(alpha=0.3)

plt.tight_layout()

output = "figures/albedo_energy_balance_test.png"
plt.savefig(output, dpi=300)
plt.close()

print()
print("Figure saved:", output)

