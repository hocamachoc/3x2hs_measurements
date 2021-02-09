#!/bin/bash

# 0) Load the conda module
# Currently supporting
# - GRID/UNESP: miniconda/3
# - NERSC/CORI: python (anaconda3)
case ${1} in
    grid)
	module load miniconda/3
	;;
    cori)
	module load python
	;;
    *)
	echo "Usage: ${0} <grid||cori>"
	exit 1
	;;
esac

# 1) Clone CosmoSIS repos
# You will need to pass de DES password for the cosmosis-des-library
mkdir -p etc/src
cd etc/src
git clone http://bitbucket.org/joezuntz/cosmosis
cd cosmosis
git clone http://bitbucket.org/joezuntz/cosmosis-standard-library
git clone https://darkenergysurvey@bitbucket.org/joezuntz/cosmosis-des-library
cd ../../../

# 2) Setup the conda env
# Note the strategy here is: if there's a `3x2pths` module already is to start
# fresh, remove it and create a new.
sed -i -e 's/^CosmoloPy/# CosmoloPy/' etc/src/cosmosis/config/requirements.txt
sed -i -e 's/^sklearn/scikit-learn/' etc/src/cosmosis/config/requirements.txt
conda remove --name 3x2pths --all -y
conda env create --name 3x2pths --file etc/env_conda.yml
source activate 3x2pths
conda install -f etc/src/cosmosis/config/requirements.txt

# 3) Install cosmosis
# We're skipping from now `polychord` as its having some compilation issues,
# sorry!
sed -i -e 's/^SUBDIRS/# SUBDIRS/' etc/src/cosmosis/cosmosis/samplers/polychord/Makefile
source etc/setup_cosmosis
cd etc/src/cosmosis
make

# 4) Clean up
# Just deactivate the conda env
conda deactivate
