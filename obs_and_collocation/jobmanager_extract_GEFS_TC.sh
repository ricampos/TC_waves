#!/bin/bash

DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/2collocation/mextract"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/2collocation/mextract"

cd ${DIRSCRIPTS}

for MBR in `seq 0 30`; do
    export MBR=${MBR}
    sbatch --output=${DIRJOUT}"/jextract_GEFS_TC_"${MBR}".out" ${DIRSCRIPTS}"/jextract_GEFS_TC.sh"
    echo " job jextract_GEFS_TC_"${MBR}" submitted OK at "$(date +"%T")
    sleep 2
done


