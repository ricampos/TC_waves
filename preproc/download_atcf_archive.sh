#!/bin/bash

########################################################################
# download_atcf_archive.sh
#
# VERSION AND LAST UPDATE:
#   v1.0  06/05/2025
#
# PURPOSE:
#  Script to download cyclone tracks, Automated Tropical Cyclone Forecast (ATCF)
#    Archive Data Files.
#
#  https://web.uwm.edu/hurricane-models/models/models.html
#  https://ftp.nhc.noaa.gov/atcf/archive/2024/
#  https://ftp.nhc.noaa.gov/atcf/archive/README
#
# USAGE:
#  Two input arguments are required: Year and destination path.
#
# OUTPUT:
#  Several atcf .dat files in the given destination path.
#
# Example:
#  This program can be run with: 
#    bash download_atcf_archive.sh 2024 /work/noaa/marine/ricardo.campos/data/
#
# DEPENDENCIES:
#  wget, gunzip
#
# AUTHOR and DATE:
#  06/05/2025: Ricardo M. Campos, first version
#
# PERSON OF CONTACT:
#  Ricardo M. Campos: ricardo.campos@noaa.gov
#
########################################################################

# Two input arguments
# Year
YEAR="$1"
# destination path
DIRW="$2"

SERVER=https://ftp.nhc.noaa.gov/atcf/archive
mdomain="aal" # (guidance information, Atlantic Ocean)

# Download index of files
wget -q ${SERVER}/${YEAR} -O index.html
# Extract list of .gz files starting with mdomain
grep -oP ${mdomain}"\d{6}\.dat\.gz" index.html | sort -u > file_list.txt

# Download each file
while read -r file; do
  echo "Downloading $file ..."
  wget -l1 -H -nd -N -np -erobots=off --tries=3 ${SERVER}/${YEAR}/$file -O $DIRW/$file
done < file_list.txt

echo "Unzipping files ..."
gunzip -f *.gz

echo "Done. Files are in: $DIRW"

