#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128G
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# This job script run the python code col_altimeter.py on Orion
#
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/modeling/PModel/reanalysis"
# export YEAR=2022
# export MONTH=6
# export CTAG="Default"
# sbatch --output=${DIRJOUT}/jreanalysis_pmodel.sh_${YEAR}_${MONTH}_${CTAG}.out jreanalysis_pmodel.sh

echo " Starting at "$(date +"%T")

ulimit -s unlimited
ulimit -c 0

DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/modeling/PModel/reanalysis"
DIRO="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/modeling/PModel/reanalysis"

export YEAR=${YEAR}
export MONTH=${MONTH}
export CTAG=${CTAG}

# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh

echo "  "
echo " reanalysis_pmodel.py ${YEAR} ${MONTH} ${CTAG} "
echo "  "

# work dir
cd ${DIRSCRIPTS}

 python3 ${DIRSCRIPTS}/reanalysis_pmodel.py ${YEAR} ${MONTH} ${CTAG}
 wait $!

echo " Complete reanalysis_pmodel.py ${YEAR} ${MONTH} ${CTAG} at  "$(date +"%T")

