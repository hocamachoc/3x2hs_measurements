#!/bin/bash
#
SEEDS=${1}
QUEUE=${2:-debug}
TMP=${3:-$(mktemp --tmpdir=${SCRATCH}/tmp -d)}
mkdir -p ${SCRATCH}/tmp

# module load python
# conda activate 3x2pths
source ${CONDA_PREFIX}/etc/setup_cosmosis

# 1) Cosmosis -> fiducial Cls
mkdir -p Cl_flaskv2p0_nolimber_emu_Nsource4
cd ${PWD}
export SCALE_CUT_DIR=${PWD}
export SCALE_CUTS="scales_all.ini"
export DATAFILE="3x2_cls_y3_maglim_metacal_nocov_v0.4.fits"
cosmosis ${PWD}/params_maglim.ini

# 2) Measurements config
cp -r cookies ${TMP} 
cp ../etc/binCDFid.txt ${TMP}

echo "type: 'flask'
nz_src: 4
nz_lns: 6
nck: 2
nside: 1024
flaskdir: '${TMP}'
maskedcatsdir: '${TMP}/maskedcats'
odir: '${TMP}/cls'
elledges: '${TMP}/binCDFid.txt'
neff: [1.476, 1.479, 1.484, 1.461]
sigma_e: [0.247, 0.266, 0.263, 0.314]
compute_cross: False
dolens: True
nonoise: False
save_maps: False
pixwin: True" >> ${TMP}/flask.yml

# 3) Flask + measurements submission file
mkdir -p ${TMP}/4096
cat <<EOF > ${TMP}/4096/submit_job${SEEDS}
#!/bin/bash
#SBATCH -q ${QUEUE}
#SBATCH --nodes=1
#SBATCH --ntasks 1 
#SBATCH --cpus-per-task 8
# #SBATCH -t 00:29:50
#SBATCH -o ${TMP}/4096/outputfile-${SEEDS}_%a
#SBATCH -e ${TMP}/4096/errorfile-${SEEDS}_%a
# #SBATCH -L SCRATCH
# #SBATCH --constraint=cpu
# #SBATCH --account=des
#SBATCH -J seed${SEEDS}
#SBATCH --array=${SEEDS}
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK}

# module load python
conda activate 3x2pths
source ${CONDA_PREFIX}/etc/setup_cosmosis
# export PMI_NO_FORK=1
# export PMI_NO_PREINITIALIZE=1

SEED=\${SLURM_ARRAY_TASK_ID}
DIROUT=${TMP}/4096/seed\${SEED}
mkdir -p \${DIROUT}

echo "DIST:      LOGNORMAL
RNDSEED:   \${SEED}
POISSON:   1" >> \${DIROUT}/tmpfile
cat \${DIROUT}/tmpfile ${PWD}/template_v2p0_maglim.config > \${DIROUT}/run.config
rm \${DIROUT}/tmpfile
sed -i 's|output|'\$DIROUT'|g' \${DIROUT}/run.config

cd ${PWD}
${CONDA_PREFIX}/bin/flask \${DIROUT}/run.config

time python3 ../flask.py ${TMP}/flask.yml --iseed \${SEED} --des_release y3 --processes \${SLURM_CPUS_PER_TASK} # $(grep -c processor /proc/cpuinfo)
for CK in 1 2 ; do
	time python3 ../3x2test.py ${TMP}/flask.yml \${SEED} \${CK}
done
EOF

# 4) Run Flask + Measurements
echo "* Output run dir: ${TMP}"
sbatch ${TMP}/4096/submit_job${SEEDS}

# Clean up
# rm -rf ${TMP}
