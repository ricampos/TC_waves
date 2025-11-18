#!/bin/bash

DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/modeling/PModel/reanalysis"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/modeling/PModel/reanalysis"

cd ${DIRSCRIPTS}

CTAG="Default"

for YEAR in `seq 2022 2024`; do
  for MONTH in `seq 6 11`; do
    export YEAR=${YEAR}
    export MONTH=${MONTH}
    export CTAG=${CTAG}
    sbatch --output=${DIRJOUT}"/jreanalysis_pmodel.sh_"${YEAR}"_"${MONTH}"_"${CTAG}".out" ${DIRSCRIPTS}"/jreanalysis_pmodel.sh"
    echo " job jreanalysis_pmodel.sh_"${YEAR}"_"${MONTH}"_"${CTAG}" submitted OK at "$(date +"%T")
    sleep 2
  done
done


