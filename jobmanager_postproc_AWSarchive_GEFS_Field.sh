#!/bin/bash

DIRJOUT="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GEFSv12Waves_AWS/jobs"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/data/archiveOPruns/GEFSv12Waves_AWS"

cd ${DIRSCRIPTS}

YEAR="2020"
for MONTH in $(seq 12 12); do
  for DAY in $(seq 1 31); do
    MM=$(printf "%02d" $MONTH)
    DD=$(printf "%02d" $DAY)
    # Check if date is valid
    if date -d "${YEAR}-${MM}-${DD}" >/dev/null 2>&1; then
      export YEAR
      export MONTH=$MM
      export DAY=$DD

      sbatch --output="${DIRJOUT}/postproc_AWSarchive_GEFS_Field_${YEAR}${MONTH}${DAY}.out" "${DIRSCRIPTS}/postproc_AWSarchive_GEFS_Field.sh"
      
      echo " Job postproc_AWSarchive_GEFS_Field_${YEAR}${MONTH}${DAY} submitted OK at $(date +"%T")"
      wait $!
      sleep 1
    fi
  done
done


