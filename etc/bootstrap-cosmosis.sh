#!/bin/bash
# This script installs cosmosis and flask inside the 3x2pths conda environment
# Note the script do not handle the creation of the environment, you will need 
# to do create it yourself.
#
# Usage:
# bash etc/bootstrap-cosmosis.sh

BASE=${PWD}
TMP=$(mktemp -d)
NJ=$(echo $(lscpu | grep ^CPU\(s\) | awk '{print $2}') - 1 | bc -l)
echo "Running with $NJ parallel jobs"

# 0) Check if we are in a conda env
if [ -z "$CONDA_PREFIX" ]; then
  echo "Oops, you are not inside a conda environment."
  echo "You may need to run first `conda create -y --file environment.yml`"
  echo "to create it, then run `conda activate 3x2pths` before trying to run"
  echo "this script again."
  exit 1
else
  echo "Running inside conda env: ${CONDA_PREFIX}"
fi

# 1) Clone the CosmoSIS repos
rm -rf ${CONDA_PREFIX}/cosmosis
# You will need to pass the DES password for the cosmosis-des-library
cd ${CONDA_PREFIX}
# 1.1) cosmosis (develop)
URL=https://bitbucket.org/joezuntz/cosmosis
git clone ${URL}
cd cosmosis
git checkout tags/des-y3 -b develop
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

# 2) Install cosmosis
SETUPCOSMOSIS=${CONDA_PREFIX}/etc/setup_cosmosis
cp -v ${BASE}/etc/setup_cosmosis.template ${SETUPCOSMOSIS}
sed -i -e "/export COSMOSIS_SRC_DIR=/cexport COSMOSIS_SRC_DIR=\\${CONDA_PREFIX}/cosmosis" ${SETUPCOSMOSIS}
source ${SETUPCOSMOSIS}
cd ${CONDA_PREFIX}/cosmosis
pip install -r ${CONDA_PREFIX}/cosmosis/config/nersc-pip-req.txt
make -j${NJ}

# 3) Install FLASK
cd ${TMP}
# 3.1) Healpix CXX
VER=3.82
URL=https://downloads.sourceforge.net/project/healpix/Healpix_${VER}/Healpix_${VER}_2022Jul28.tar.gz
curl -LO $URL
tar -xzf Healpix_${VER}_2022Jul28.tar.gz
cd Healpix*/
FITSDIR=${CONDA_PREFIX}/lib FITSINC=${CONDA_PREFIX}/include ./configure -L --auto=sharp
cd src/cxx
CFITSIO_CFLAGS=-I$CONDA_PREFIX/include CFITSIO_LIBS=-L$CONDA_PREFIX/lib SHARP_CFLAGS=-I$PWD/../../include SHARP_LIBS=-L$PWD/../../ ./configure --prefix=$PWD/../../
make -j${NJ} install
cd ../../
cp -ur include lib bin data Version ${CONDA_PREFIX}/
rm /tmp/Healpix_autolist.txt
# 3.2) flask
cd ${TMP}
# URLREPO=https://github.com/ucl-cosmoparticles/flask.git # Does not read CATALOG_COLS(?)
URLREPO=https://github.com/hsxavier/flask.git
git clone $URLREPO
cd flask/src
sed -i -e "/HEALDIR  = \/home\/skems\/prog\/Healpix_3.60_2019Dec18\/Healpix_3.60/cHEALDIR = \\${CONDA_PREFIX}" Makefile
sed -i -e '/CXXHEAL  = -I$(HEALDIR)\/src\/cxx\/generic_gcc\/include/cCXXHEAL = -I$(HEALDIR)/include/healpix_cxx' Makefile
sed -i -e '/LDHEAL   = -L$(HEALDIR)\/src\/cxx\/generic_gcc\/lib/cLDHEAL = -L$(HEALDIR)/lib' Makefile
sed -i -e "/#CXXFITS  = -I\/mnt\/c\/User\/arthu\/Work\/lib\/cfitsio\//cCXXFITS = -I\\${CONDA_PREFIX}/include" Makefile
sed -i -e "/#LDFITS   = -L\/mnt\/c\/User\/arthu\/Work\/lib\/cfitsio\//cLDFITS = -L\\${CONDA_PREFIX}/lib" Makefile
sed -i -e "/#CXXGSL   = -I\/path\/to\/gsl\/headers\/folder/cCXXGSL = -I\\${CONDA_PREFIX}/include" Makefile
sed -i -e "/#LDGSL    = -L\/path\/to\/gsl\/libraries\/folder/cLDGSL = -L\\${CONDA_PREFIX}/lib" Makefile
make -j${NJ}
cd ..
cp -ru bin ${CONDA_PREFIX}/

# 4) Clean up
rm -rf ${TMP}
