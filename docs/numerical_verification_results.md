# Numerical Verification Results

## Research project

Computational Modelling of Urban Heat Island Mitigation
in SADC Cities

MSc Computational Physics — University of Eswatini

## Verification environment

The numerical verification programs were compiled
using GNU Fortran through the repository Makefile.

The tests were executed in a temporary working directory
with the required output folders.

All four verification programs completed with exit code 0.

## Analytical heat-equation verification

Program: verify_heat2d

Grid: 51 x 51
Final simulation time: 200000 s

Numerical maximum: 29.8433716268 degrees C
Analytical maximum: 29.8433266287 degrees C

RMSE: 2.2057917728e-5 degrees C

Maximum absolute error: 4.4998152148e-5 degrees C

## Spatial grid convergence

Program: grid_convergence

| Grid refinement | Observed order |
|---|---:|
| 26 to 51 | 1.971455 |
| 51 to 101 | 1.985653 |

The observed spatial convergence approaches
the theoretical second-order accuracy.

## Temporal convergence

Program: time_convergence

| Time-step refinement (s) | Observed order |
|---|---:|
| 1000 to 500 | 1.000038 |
| 500 to 250 | 1.000019 |
| 250 to 125 | 1.000009 |

The observed temporal convergence is consistent
with first-order time integration.

## Coupled spatial operator verification

Program: verify_coupled_spatial_operator

| Grid refinement | Observed order |
|---|---:|
| 20 to 40 | 1.999110 |
| 40 to 80 | 1.999778 |
| 80 to 160 | 1.999944 |
| 160 to 320 | 1.999986 |

The results are consistent with second-order
spatial convergence of the tested operator.

## Scientific interpretation

The numerical experiments support the correctness
of the tested discretisation implementations under
their respective verification conditions.

They do not establish independent observational
validation of the urban heat island model.

The full three-city mitigation workflow requires
separate end-to-end reproducibility testing.
