#!/bin/bash
# last = 614,   debug, ate 659
# regular = 660 679
for i in $(seq 701 900 )
do
dir_out="/global/cscratch1/sd/faoli/flask_desy3/4096/seed${i}/"
# dir_out=$SCRATCH/..._seed${i}/
mkdir ${dir_out}

echo '#!/bin/bash' > submit_job${i}

echo "#SBATCH -q regular #debug or regular" >> submit_job${i}
echo "#SBATCH --nodes=1" >> submit_job${i}

echo "#SBATCH -t 00:29:50">> submit_job${i}
echo "#SBATCH -o ${dir_out}outputfile" >> submit_job${i}
echo "#SBATCH -e ${dir_out}errorfile" >> submit_job${i}
echo "#SBATCH -L SCRATCH" >> submit_job${i}
echo "#SBATCH --constraint=haswell" >> submit_job${i}
echo "#SBATCH --account=des" >> submit_job${i}
echo "#SBATCH -J seed${i}" >> submit_job${i}
echo "#SBATCH --mail-user=felipeaoli@gmail.com" >> submit_job${i}
echo "#SBATCH --mail-type=END" >> submit_job${i}

#echo "source ~/.profile.ext" >> submit_job${i}
echo "module load gsl" >> submit_job${i}
echo "export PMI_NO_FORK=1" >> submit_job${i}
echo "export PMI_NO_PREINITIALIZE=1" >> submit_job${i}
echo "cd $PWD" >> submit_job${i}

echo "DIST:      LOGNORMAL" > ${dir_out}/tmpfile
echo "RNDSEED:   $i" >> ${dir_out}/tmpfile
echo "POISSON:   1"  >> ${dir_out}/tmpfile

cat ${dir_out}/tmpfile $PWD/template_v2p0_yo.config > ${dir_out}/run.config
rm ${dir_out}/tmpfile

sed -i 's|output|'$dir_out'|g' ${dir_out}run.config
echo "/global/homes/f/faoli/flask-root/flask/bin/flask ${dir_out}/run.config">> submit_job${i}
## echo "flask ${dir_out}/run.config">> submit_job${i}
cp submit_job${i} ${dir_out}
sbatch submit_job${i}
done

