#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# This job run the code hyper_optimz_probmod.py on Orion

ulimit -s unlimited
ulimit -c 0

DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/4postproc/AI_ML/jobs"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/4postproc/AI_ML"
export DIRJOUT=${DIRJOUT}
export DIRSCRIPTS=${DIRSCRIPTS}
mkdir -p "${DIRJOUT}"
cd "${DIRSCRIPTS}" || exit 1

export rst=${rst}
export hl=${hl}
export alp=${alp}
export lri=${lri}

echo " Starting at "$(date +"%T")

# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh

# run
/work/noaa/marine/ricardo.campos/progs/python/anaconda3/bin/python3 ${DIRSCRIPTS}/optimize_NNMLP.py ${rst} ${hl} ${alp} ${lri}

wait $!
sleep 1

echo " Complete at "$(date +"%T")

