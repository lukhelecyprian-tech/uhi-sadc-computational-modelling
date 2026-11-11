# Computational Modelling of Urban Heat Island Mitigation in SADC Cities

MSc Computational Physics Research — University of Eswatini

## Project Overview

This repository contains Python and Fortran source codes for investigating urban heat island behaviour and mitigation in Mbabane, Manzini and Johannesburg.

The research combines Landsat-derived environmental data with numerical modelling of surface-energy processes and lateral heat diffusion.

## Research Components

- Satellite-derived LST, NDVI, albedo and surface reflectance
- Surface-energy and coupled heat-diffusion models
- Numerical verification and convergence studies
- Physical parameter sensitivity analyses
- Urban heat island mitigation scenarios
- Scientific visualisation and comparison of results

## Repository Structure

- `fortran/` — numerical models, verification and sensitivity programs
- `python/` — preprocessing, analysis and plotting scripts

## Software Requirements

Python libraries: NumPy, Pandas, Matplotlib, SciPy, Rasterio and GeoPandas.

Fortran compiler used: GNU Fortran 13.3.0.

## Reproducibility

This repository provides the Python and Fortran source
codes developed for an ongoing MSc Computational Physics
research project at the University of Eswatini.

The computational workflow consists of:

1. Processing Landsat-derived LST, NDVI, albedo and
   surface reflectance data.
2. Preparing spatial input fields using Python.
3. Running numerical models and mitigation scenarios
   using Fortran.
4. Performing numerical verification and sensitivity
   analyses.
5. Generating scientific figures using Python.

The source codes are organised into the `python/`
and `fortran/` directories.

Some programs require external satellite datasets,
prepared input files and specific directory structures.
Large raster datasets and generated simulation outputs
are not distributed in this repository.

The repository is a research-in-progress archive.
Source codes, documentation and results may be updated
following academic review.

Full reproducibility instructions will be added as
the research is finalised.
## Author

Cyprian Lukhele
MSc Computational Physics
University of Eswatini
