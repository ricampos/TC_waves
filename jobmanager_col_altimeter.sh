#!/bin/bash

DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/2collocation/CDIP"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/2collocation/CDIP"

cd ${DIRSCRIPTS}

NSEG=10
for SAT in 0 2 5 6 7 14 15 16; do
  for SEG in `seq 1 10`; do
    export NSEG=${NSEG}
    export SAT=${SAT}
    export SEG=${SEG}
    sbatch --output=${DIRJOUT}"/jcol_altimeter_"${SAT}"_"${SEG}"_"${NSEG}".out" ${DIRSCRIPTS}"/jcol_altimeter.sh"
    echo " job jcol_altimeter_"${SAT}"_"${SEG}"_"${NSEG}" submitted OK at "$(date +"%T")
    sleep 2
  done
done


