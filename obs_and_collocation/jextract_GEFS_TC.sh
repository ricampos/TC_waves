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
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/2collocation/mextract"
# export MBR=0
# sbatch --output=${DIRJOUT}/jextract_GEFS_TC_${MBR}.out jextract_GEFS_TC.sh

echo " Starting at "$(date +"%T")

ulimit -s unlimited
ulimit -c 0

DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/2collocation/mextract"
DIRO="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/2collocation/mextract"

export MBR=${MBR}

# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh

echo "  "
echo " Starting extract_GEFS_TC.py ${MBR} "
echo "  "

# work dir
cd ${DIRSCRIPTS}

 python3 ${DIRSCRIPTS}/extract_GEFS_TC.py ${MBR}
 wait $!

echo " Complete extract_GEFS_TC.py ${MBR} at  "$(date +"%T")

