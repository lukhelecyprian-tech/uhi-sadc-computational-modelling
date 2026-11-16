# Numerical Model

## Governing equation

The reduced surface-temperature model is

dT/dt = kappa * Laplacian(T)
      + [(1-alpha)*S
      - epsilon*sigma*(T^4 - Tsky^4)
      - h*(T-Tair)
      - LE] / CA

where:

- T: surface temperature (K)
- kappa: effective lateral thermal diffusivity
- alpha: surface albedo
- S: incoming shortwave radiation
- epsilon: surface emissivity
- sigma: Stefan-Boltzmann constant
- Tsky: prescribed sky temperature
- Tair: prescribed air temperature
- h: sensible-heat exchange coefficient
- LE: prescribed latent-heat cooling flux
- CA: effective areal heat capacity

## Numerical method

The model uses a two-dimensional finite difference
approximation for lateral thermal diffusion.

Forward Euler time integration is used.

## Reference numerical settings

| Parameter | Value |
|---|---|
| Grid spacing | 30 m |
| Time step | 60 s |
| Simulation duration | 6 hours |
| Time steps | 360 |
| Effective thermal diffusivity | 1.0e-5 m2/s |
| Incoming shortwave radiation | 500 W/m2 |
| Areal heat capacity | 2.0e6 J/m2/K |
| Reference emissivity | 0.95 |
| Sensible-heat coefficient | 10 W/m2/K |
| Reference latent-heat flux | 50 W/m2 |

These are reference experimental settings, not universal
physical constants for the three cities.

## Mitigation scenarios

S0: Baseline without additional mitigation.

S1: Increased surface albedo.

S2: Enhanced latent-heat cooling.

S3: Combined albedo and latent-heat intervention.

For the reference simulations, mitigation eligibility
is defined using:

NDVI < 0.30

The same threshold is held constant for S1, S2 and S3.

Only eligible valid grid cells receive the prescribed
intervention.

Reference treatment increments:

- Albedo increase: up to 0.10
- Additional latent-heat flux: 30 W/m2

Separate sensitivity experiments investigate changes
in the NDVI eligibility threshold and physical parameters.

## Verification

The research includes spatial and temporal refinement,
known-solution tests and coupled-operator verification.

The expected orders for the tested numerical methods are:

- Second-order spatial discretisation
- First-order forward Euler temporal discretisation

These statements describe numerical verification and
should not be confused with observational validation.
