# Research Methodology

## Research project

Computational Modelling of Urban Heat Island Mitigation in SADC Cities

MSc Computational Physics — University of Eswatini

## Study areas

The research investigates three Southern African urban areas:

1. Mbabane, Eswatini
2. Manzini, Eswatini
3. Johannesburg, South Africa

## Satellite observations

Landsat 8 Collection 2 Level-2 observations provide the
environmental inputs used in the study.

Reference acquisition dates:

| City | Acquisition date |
|---|---|
| Mbabane | 11 July 2025 |
| Manzini | 8 May 2025 |
| Johannesburg | 22 May 2025 |

The principal variables are:

- Land surface temperature (LST)
- Normalised Difference Vegetation Index (NDVI)
- Surface albedo
- Surface reflectance

The computational analysis uses a nominal 30 m grid.

The native thermal resolution is coarser than 30 m.
The different acquisition dates prevent a controlled
comparison of urban heat island intensity between cities.

## Computational workflow

1. Obtain and process Landsat observations.
2. Prepare spatially varying input fields using Python.
3. Implement the reduced surface-energy and diffusion model
   using Fortran.
4. Conduct baseline and mitigation simulations.
5. Perform numerical verification and sensitivity analysis.
6. Analyse and visualise results using Python.

## Model scope

The numerical framework represents selected surface-energy
exchange processes and lateral thermal diffusion.

It is intended for controlled computational experiments.

It is not a fully coupled atmospheric or urban climate model.

## Scientific limitations

Numerical verification does not establish independent
physical validation.

The prescribed mitigation scenarios are idealised.
Their simulated cooling responses should not be interpreted
as measured effects of installed interventions.
