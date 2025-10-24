#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -p orion

# DIRJOUT="/work/noaa/marine/ricardo.campos/data/archiveOPruns/HAFS_AWS/jobs"
# sbatch --output=${DIRJOUT}/postproc_AWSarchive_HAFS_Field.out postproc_AWSarchive_HAFS_Field.sh

echo " Starting at "$(date +"%T")

ulimit -s unlimited
ulimit -c 0

module load nco
module load cdo
module load wgrib2

# Paths and list of filest to be converted
DIRS="/work/noaa/marine/ricardo.campos/data/archiveOPruns/HAFS_AWS/grib2" # archive path
DIRO="/work/noaa/marine/ricardo.campos/data/archiveOPruns/HAFS_AWS/netcdf" # final netcdf4 output path
file_list="/work/noaa/marine/ricardo.campos/data/archiveOPruns/HAFS_AWS/list.txt"

# cutoff decimals to reduce file size
dp=2

# Post-processing: compress, select variables, and reduce decimals resolution, to save disk space. ------------------
echo " "
echo " Post-Processing. Netcdf4 compression "
cd $DIRS

# Read each line (i.e., each file name) and process it
while IFS= read -r file; do
  echo "Processing file: $file"
  arqn="${file%.grb2}"

  test -f ${file}
  TE=$?
  if [ ${TE} -eq 1 ]; then
    echo " File ${file} does not exist. Download failed."
    echo " File ${file} does not exist. Download failed." >> "${DIRS}/file_doesnt_extist.txt"
    exit 1
  else
    wgrib2 ${file} -netcdf ${arqn}.saux.nc 2>&1
    wait $!
    cdo delvar,WDIR_surface ${arqn}.saux.nc ${arqn}.saux2.nc
    wait $!
    ncks -4 -L 1 ${arqn}.saux2.nc ${arqn}.saux3.nc 2>&1
    wait $!
    ncks --ppc default=.$dp ${arqn}.saux3.nc ${arqn}.nc 2>&1
    wait $!
    ncatted -a _FillValue,,o,f,NaN ${arqn}.nc 2>&1
    wait $!
    rm -f ${arqn}.saux*
    rm -f ${arqn}.*idx*
    echo " File ${arqn} converted to netcdf and compressed with success. "
    sleep 1
    mv ${arqn}.nc ${DIRO}
  fi

done < "$file_list"


echo "Completed postproc_AWSarchive_HAFS_Field.sh at "$(date +"%T")

