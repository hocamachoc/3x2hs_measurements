## Installation

Run the `etc/bootstrap-cosmosis.sh` script, it should install a `3x2pths` environment containing all the required software.

The script provides functionality for NERSC-CORI, SDumont and GRID-UNESP, other computational facilities should not be too dificult to add as well.
The script expects to be run from the root folder of the repository.
For example, for NERSC-CORI, simply run:

```{sh}
./etc/bootstrap-cosmosis.sh cori
```

Relevant files:

- conda environment: `environment.yml`

- Bootstraping script: `etc/bootstrap-cosmosis.sh` 

## FLASK E2E pipeline: 

AIM: Produce, process, measure and remove Lognormal catalogs 

## 3x2hs_measurements
