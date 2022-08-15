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

### How to run

We use [slurm's job array syntax](https://slurm.schedmd.com/job_array.html) to specify a set of seeds to run.
For example, to run the seeds 1 to 5,

```{sh}
cd conf
./run_2p0.sh 1-5
```

## 3x2hs_measurements
