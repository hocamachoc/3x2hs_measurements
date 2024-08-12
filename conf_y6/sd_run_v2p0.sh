#!/bin/bash
#
SEEDS=${1}
#QUEUE=${2:-debug}
TMP=${3:-$(mktemp --tmpdir=${SCRATCH}/tmp -d)}
mkdir -p ${SCRATCH}/tmp

#module load python
#conda activate 3x2pths
#source ${CONDA_PREFIX}/etc/setup_cosmosis

## 1) Cosmosis -> fiducial Cls
#mkdir -p Cl_flaskv2p0_nolimber_emu_Nsource4
#cd ${PWD}
#export SCALE_CUT_DIR=${PWD}
#export SCALE_CUTS="scales_all.ini"
#export DATAFILE="sim_3x2_fiducial_nla.fits"
#cosmosis ${PWD}/params.ini

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
neff: [2.1402, 2.14455, 2.1518, 2.11845]
sigma_e: [0.27216, 0.29344, 0.29008, 0.33712]
compute_cross: False
dolens: True
nonoise: False
save_maps: False
pixwin: True" >> ${TMP}/flask.yml

# 3) Flask + measurements submission file
mkdir -p ${TMP}/4096
cat <<EOF > ${TMP}/4096/submit_job${SEEDS}
#!/bin/bash
#SBATCH -N 1
#SBATCH --tasks-per-node=24
#SBATCH -p cpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user lucas.faga@usp.br
#SBATCH -t 01:00:00
#SBATCH -o ${TMP}/4096/outputfile-${SEEDS}_%a
#SBATCH -e ${TMP}/4096/errorfile-${SEEDS}_%a
#SBATCH -J seed${SEEDS}
#SBATCH --array=${SEEDS}

echo $SHELL
echo ${SLURM_NTASKS}
# module load python
source ${HOME}/.bashrc
conda activate /prj/eubd/hugo.chavez2/micromamba/envs/3x2pths
source ${CONDA_PREFIX}/etc/setup_cosmosis
export OMP_NUM_THREADS=2
echo "Number of tasks:" $SLURM_NTASKS

SEED=\${SLURM_ARRAY_TASK_ID}
DIROUT=${TMP}/4096/seed\${SEED}
mkdir -p \${DIROUT}

echo "DIST:      LOGNORMAL
RNDSEED:   \${SEED}
POISSON:   1" >> \${DIROUT}/tmpfile
cat \${DIROUT}/tmpfile ${PWD}/template_v2p0.config > \${DIROUT}/run.config
rm \${DIROUT}/tmpfile
sed -i 's|output|'\$DIROUT'|g' \${DIROUT}/run.config

cd ${PWD}
${CONDA_PREFIX}/bin/flask \${DIROUT}/run.config

# LJF - measurement part starts below. I'm commenting it for now.
#time python3 ../flask.py ${TMP}/flask.yml --iseed \${SEED} --des_release y3 --processes 10 	# $(grep -c processor /proc/cpuinfo)
#for CK in 1 2 ; do
#	time python3 ../3x2test.py ${TMP}/flask.yml \${SEED} \${CK}
#done
EOF

# 4) Run Flask + Measurements
cd /scratch/eubd/lucas.faga/y3-3x2pt_harmonic/cosmosis/
echo "* Output run dir: ${TMP}"
sbatch ${TMP}/4096/submit_job${SEEDS}

# Clean up
# rm -rf ${TMP}
