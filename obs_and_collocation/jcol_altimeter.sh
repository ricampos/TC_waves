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
# DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/2collocation/CDIP"
# export NSEG=10
# export SEG=1
# export SAT=0
# sbatch --output=${DIRJOUT}/jcol_altimeter_${SAT}_${SEG}_${NSEG}.out jcol_altimeter.sh

echo " Starting at "$(date +"%T")

ulimit -s unlimited
ulimit -c 0

DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/2collocation/CDIP"
DIRO="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/2collocation/CDIP"

export SAT=${SAT}
export NSEG=${NSEG}
export SEG=${SEG}

# python env
source /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh
sh /work/noaa/marine/ricardo.campos/progs/python/setanaconda3.sh

echo "  "
echo " Starting col_altimeter.py ${SAT} ${NSEG} ${SEG} "
echo "  "

# work dir
cd ${DIRSCRIPTS}

 python3 ${DIRSCRIPTS}/col_altimeter.py ${SAT} ${NSEG} ${SEG}
 wait $!

echo " Complete col_altimeter.py ${SAT} ${SEG} ${NSEG} at  "$(date +"%T")

