#!/bin/bash

DIRJOUT="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GEFSv12Waves_AWS/jobs"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GEFSv12Waves_AWS"

cd ${DIRSCRIPTS}

YEAR="2022"
for MONTH in `seq 8 9`; do
  for DAY in `seq 1 31`; do
    export YEAR=${YEAR}
    export MONTH=${MONTH}
    export DAY=${DAY}
    sbatch --output="${DIRJOUT}/postproc_AWSarchive_GEFS_Field_${YEAR}${MONTH}${DAY}.out" "${DIRSCRIPTS}/postproc_AWSarchive_GEFS_Field.
sh"
      
    echo " Job postproc_AWSarchive_GEFS_Field_${YEAR}${MONTH}${DAY} submitted OK at $(date +"%T")"
    wait $!
    sleep 1
  done
done

