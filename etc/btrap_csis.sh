#!/bin/bash

# 0) Load the conda (on grid we use miniconda/3)
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
conda remove --name 3x2pths --all -y

# 1) Clone CosmoSIS repos
# You will need to pass de DES password for the cosmosis-des-library
mkdir -p etc/src
cd etc/src
git clone http://bitbucket.org/joezuntz/cosmosis
cd cosmosis
git clone http://bitbucket.org/joezuntz/cosmosis-standard-library
git clone https://darkenergysurvey@bitbucket.org/joezuntz/cosmosis-des-library
cd ../../../

# # 2) Setup conda env
case ${1} in
    grid)
	# >>> conda initialize >>>
	# !! Contents within this block are managed by 'conda init' !!
	__conda_setup="$('/opt/gridunesp/dist/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
	if [ $? -eq 0 ]; then
	    eval "$__conda_setup"
	else
	    if [ -f "/opt/gridunesp/dist/miniconda3/etc/profile.d/conda.sh" ]; then
		. "/opt/gridunesp/dist/miniconda3/etc/profile.d/conda.sh"
	    else
		export PATH="/opt/gridunesp/dist/miniconda3/bin:$PATH"
	    fi
	fi
	unset __conda_setup
	# <<< conda initialize <<<
	;;
    cori)
	# >>> conda initialize >>>
	# !! Contents within this block are managed by 'conda init' !!
	__conda_setup="$('/usr/common/software/python/3.8-anaconda-2020.11/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
	if [ $? -eq 0 ]; then
	    eval "$__conda_setup"
	else
	    if [ -f "/usr/common/software/python/3.8-anaconda-2020.11/etc/profile.d/conda.sh" ]; then
		. "/usr/common/software/python/3.8-anaconda-2020.11/etc/profile.d/conda.sh"
	    else
		export PATH="/usr/common/software/python/3.8-anaconda-2020.11/bin:$PATH"
	    fi
	fi
	unset __conda_setup
	# <<< conda initialize <<<
	;;
    *)
	echo "Usage: ${0} <grid||cori>"
	exit 1
	;;
esac
conda deactivate

sed -i -e 's/^CosmoloPy/# CosmoloPy/' etc/src/cosmosis/config/requirements.txt
sed -i -e 's/^sklearn/scikit-learn/' etc/src/cosmosis/config/requirements.txt
conda env create -n 3x2pths -f etc/env_conda3.yml
conda activate 3x2pths
conda install -f etc/src/cosmosis/config/requirements.txt

# 3) Install cosmosis
sed -i -e 's/^SUBDIRS/# SUBDIRS/' etc/src/cosmosis/cosmosis/samplers/polychord/Makefile
source etc/setup_cosmosis
cd etc/src/cosmosis
make
conda deactivate
