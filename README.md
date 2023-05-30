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

We use the `run_v2p0.sh` script for the E2E pipeline.

The first argument should especify the set of seeds to run.
We use [slurm's job array syntax](https://slurm.schedmd.com/job_array.html) to specify a set of seeds to run.

The second argument sets the slurm queue for the run.

A third (optional) argument can specify the directory where the outputs will be placed.
If no third argument is given, the code will set a temporary folder.

For example, to run the seeds 1 to 5 on the `debug` queue,

```{sh}
cd conf
./run_v2p0.sh 1-5 debug
```

## 3x2hs_measurements

## Gaussian covariance

** still under development **

TO DO: apply shot-noise for clustering

In order to calculate clustering and GGL gaussian covariances, there are two scripts that you have to run. 
The first of them (`save_cov_ws.py`) saves each covariance workspace, which is the most time consuming operation in the computation. The other (`2x2gcovtest.py`) gets the covariance and save all blocks in a numpy npz file. Currently, for both codes you have to specify which sector of the covariance you want: gcl-gcl for clustering, ggl-ggl for galaxy-galaxy lensing and gcl-ggl for their cross covariance. 

To run the `save_cov_ws.py`:

```{sh}
./save_cov_ws.sh 0-25  'gcl-gcl' regular # for clustering-clustering
./save_cov_ws.sh 0-100 'gcl-ggl' regular # for clustering-ggl
./save_cov_ws.sh 0-210 'ggl-ggl' regular # for ggl-ggl
```

The first argument is the blocks of the covariance to be calculated (each zbin combination is a block), the second argument is sector of the covariance you want and the third is the type of queue you want in slurm (use debug if you are computing up to 5 blocks). 

There are: 
- 25 unique blocks for gcl-gcl (only "diagonal blocks")
- 100 unique blocks for gcl-ggl; they are the same as ggl-gcl because cov matrix is symmetric
- 210 unique blocks for ggl-ggl; these account for one of the triangles of the cov matrix

To run the `2x2gcovtest.py` symply do:

```{sh}
sbatch 2x2gcovtest.sh
```
