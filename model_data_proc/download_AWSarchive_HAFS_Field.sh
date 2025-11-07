#!/bin/bash

# ---------------
# Download HAFS wave (fields) from AWS
# https://wpo.noaa.gov/the-hurricane-analysis-and-forecast-system-hafs/
# It requires aws s3 to be installed
# ---------------

# 3 input arguments
# cycle date and time
CTIME="$1" # YYYYMMDD
CHOUR="$2" # HH 0,6,12,18
# destination path
DIRW="$3"

# https://noaa-nws-hafs-pds.s3.amazonaws.com/index.html#hfsa/20241010/06/
# server address
SERVER=https://noaa-nws-hafs-pds.s3.amazonaws.com
# model
model="hfsa" 

cd ${DIRW}

mapfile -t fnames < <(aws s3 ls s3://noaa-nws-hafs-pds/hfsa/${CTIME}/"$(printf "%02.f" $CHOUR)"/ --no-sign-request | awk '{print $4}' | grep "\.${model}\.ww3\.grb2$")
if [ ${#fnames[@]} -eq 0 ]; then
  echo "No matching files found. "${CTIME}$(printf "%02.f" $CHOUR)
else
  for fname in "${fnames[@]}"; do
    echo "Downloading $fname..."
    wget -l1 -H -t1 -nd -N -np -erobots=off --tries=3 ${SERVER}/${model}/${CTIME}/"$(printf "%02.f" $CHOUR)"/${fname} -O $DIRW/${fname} 2>&1
  done
fi

find $DIR -empty -type f -delete

