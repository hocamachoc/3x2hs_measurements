#!/bin/bash

# 0) Load the conda module
# Currently supporting: GRID/UNESP: miniconda/3, CORI/NERSC: python (anaconda3) 
# & SD anaconda3.
case ${1} in
    grid)
	module load miniconda/3
	NJ=16	# Number of allowed parallel jobs
	;;
    cori)
	module load python
	NJ=16	# Number of allowed parallel jobs
	;;
    sd)
	module load anaconda3
	NJ=16	# Number of allowed parallel jobs
	;;
    *)
	echo "Usage: ${0} <grid|cori|sd>"
	exit 1
	;;
esac

BASE=${PWD}

# 1) Setup the conda env
# Strategy here: if there's a `3x2pths` module, start fresh, remove it and
# create a new one - TODO: There should be a better way...
# conda remove --name 3x2pths --all -y
# conda env create --name 3x2pths --file environment.yml
source activate 3x2pths

# 2) Clone the CosmoSIS repos
rm -rf ${CONDA_PREFIX}/cosmosis
# You will need to pass the DES password for the cosmosis-des-library
cd ${CONDA_PREFIX}
# 1.1) cosmosis (develop)
URL=https://bitbucket.org/joezuntz/cosmosis
git clone ${URL}
cd cosmosis
git checkout develop
# 1.2) cosmosis-standard-library (des-y3)
URL=https://bitbucket.org/joezuntz/cosmosis-standard-library
git clone ${URL}
cd cosmosis-standard-library
git checkout des-y3
cd ..
# 1.3) cosmosis-des-library (master)
URL=git@bitbucket.org:joezuntz/cosmosis-des-library.git
git clone ${URL}
cd cosmosis-des-library
git checkout master
cd ..
# Cosmosis requirements
pip install -r ${CONDA_PREFIX}/cosmosis/config/nersc-pip-req.txt

# 3) Install cosmosis
SETUPCOSMOSIS=${CONDA_PREFIX}/etc/setup_cosmosis
cp -v ${BASE}/etc/setup_cosmosis.template ${SETUPCOSMOSIS}
sed -i -e "/export COSMOSIS_SRC_DIR=/cexport COSMOSIS_SRC_DIR=\\${CONDA_PREFIX}/cosmosis" ${SETUPCOSMOSIS}
# sed -i -e "s#export PATH=#export PATH=${PWD}/build/bin:#" ${SETUPCOSMOSIS}
# sed -i -e "s#export LD_LIBRARY_PATH=#export LD_LIBRARY_PATH=${PWD}/build/lib:#" ${SETUPCOSMOSIS}
source ${SETUPCOSMOSIS}
cd ${CONDA_PREFIX}/cosmosis
make -j${NJ}

# 4) Install FLASK
TMP=$(mktemp --tmpdir=${CONDA_PREFIX}/tmp -d)
cd ${TMP}
# 4.1) Healpix CXX
VER=3.82
URL=https://downloads.sourceforge.net/project/healpix/Healpix_3.82/Healpix_${VER}_2022Jul28.tar.gz
wget -c --no-check-certificate $URL
tar -xzf *.tar.gz
cd Healpix*/
FITSDIR=${CONDA_PREFIX}/lib FITSINC=${CONDA_PREFIX}/include \
	./configure -L --auto=cxx
make -j${NJ}
cp -ur include lib bin data Version ${CONDA_PREFIX}/
# 4.2) flask
cd ${TMP}
URLREPO=https://github.com/ucl-cosmoparticles/flask.git
git clone $URLREPO
cd flask/src
sed -i -e "/HEALDIR  =  \/home\/jayesh\/Documents\/Euclid\/Software\/Healpix_3.50/cHEALDIR = \\${CONDA_PREFIX}" Makefile
sed -i -e '/CXXHEAL  = -I$(HEALDIR)\/src\/cxx\/generic_gcc\/include/cCXXHEAL = -I$(HEALDIR)/include/healpix_cxx' Makefile
sed -i -e '/LDHEAL   = -L$(HEALDIR)\/src\/cxx\/generic_gcc\/lib/cLDHEAL = -L$(HEALDIR)/lib' Makefile
sed -i -e "/#CXXFITS  = -I\/mnt\/c\/User\/arthu\/Work\/lib\/cfitsio\//cCXXFITS = -I\\${CONDA_PREFIX}/include" Makefile
sed -i -e "/#LDFITS   = -L\/mnt\/c\/User\/arthu\/Work\/lib\/cfitsio\//cLDFITS = -L\\${CONDA_PREFIX}/lib" Makefile
sed -i -e "/#CXXGSL   = -I\/path\/to\/gsl\/headers\/folder/cCXXGSL = -I\\${CONDA_PREFIX}/include" Makefile
sed -i -e "/#LDGSL    = -L\/path\/to\/gsl\/libraries\/folder/cLDGSL = -L\\${CONDA_PREFIX}/lib" Makefile
make -j${NJ}
cd ..
cp -ru bin ${CONDA_PREFIX}/
# Finalize
rm -rf ${TMP}

# 4) Clean up
conda deactivate
