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

Some programs require external raster datasets and prepared numerical inputs. Large datasets and generated outputs are not included in this repository. Complete reproduction instructions and numerical verification results will be documented separately.

## Author

Cyprian Lukhele
MSc Computational Physics
University of Eswatini
