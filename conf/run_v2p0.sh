#!/bin/bash
#
SEED=${1}
TMP=$(mktemp --tmpdir=${SCRATCH}/tmp -d)
DIROUT=${TMP}/seed${SEED}

conda activate 3x2pths
source ${CONDA_PREFIX}/etc/setup_cosmosis

mkdir -p ${DIROUT}

# 1) Cosmosis -> fiducial Cls
cd ${PWD}
export SCALE_CUT_DIR=${PWD}
export SCALE_CUTS="scales_all.ini"
export DATAFILE="sim_3x2_fiducial_nla.fits"
cosmosis ${PWD}/params.ini

# 2) Flask config file
echo "DIST:      LOGNORMAL
RNDSEED:   ${SEED}
POISSON:   1" >> ${DIROUT}/tmpfile
cat ${DIROUT}/tmpfile $PWD/template_v2p0.config > ${DIROUT}/run.config
rm ${DIROUT}/tmpfile
sed -i 's|output|'$DIROUT'|g' ${DIROUT}/run.config

# 3) Flask submission file
echo "#!/bin/bash
#SBATCH -q debug #debug or regular
#SBATCH --nodes=1
#SBATCH -t 00:29:50
#SBATCH -o ${DIROUT}/outputfile
#SBATCH -e ${DIROUT}/errorfile
#SBATCH -L SCRATCH
#SBATCH --constraint=haswell
# #SBATCH --account=des
#SBATCH -J seed${SEED}
#SBATCH --mail-type=ALL

module load python
conda activate 3x2pths
source ${CONDA_PREFIX}/etc/setup_cosmosis
export PMI_NO_FORK=1
export PMI_NO_PREINITIALIZE=1

cd ${PWD}
${CONDA_PREFIX}/bin/flask ${DIROUT}/run.config" >> ${DIROUT}/submit_job${SEED}

# 4) Run Flask
echo "* Output run dir: ${DIROUT}"
sbatch ${DIROUT}/submit_job${SEED}

# Clean up
# rm -rf ${TMP}
