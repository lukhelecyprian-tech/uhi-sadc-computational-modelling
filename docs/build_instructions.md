# Fortran Build and Reproducibility Instructions

## Research project

Computational Modelling of Urban Heat Island Mitigation
in SADC Cities.

MSc Computational Physics, University of Eswatini.

## Software requirements

- GNU Fortran (gfortran)
- GNU Make
- Linux or WSL Ubuntu

## Building the numerical programs

From the repository root, execute:

    make -j2

This compiles the Fortran source files and places
the generated executables in the build/ directory.

## Checking the executables

    find build -maxdepth 1 -type f -executable | wc -l

The current source collection contains 33 Fortran
programs.

## Cleaning the build directory

    make clean

## Running numerical experiments

Compilation and execution are separate steps.

Some research programs require prepared satellite
input files, spatial masks, or other external data.

These datasets are not distributed in the repository.

The exact runtime requirements and execution commands
must be documented and tested separately.

## Verification status

All 33 Fortran source files successfully compiled
and linked during the November 2026 build tests.

All 34 Python source files passed syntax checking.

These checks do not establish numerical correctness,
independent physical validation, or complete
end-to-end reproducibility.
