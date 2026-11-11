#!/bin/bash

#  simulations
DIRJOUT="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/4postproc/AI_ML/jobs"
DIRSCRIPTS="/work/noaa/marine/ricardo.campos/work/analysis/TC_waves/4postproc/AI_ML"
export DIRJOUT=${DIRJOUT}
export DIRSCRIPTS=${DIRSCRIPTS}
mkdir -p "${DIRJOUT}"
cd "${DIRSCRIPTS}" || exit 1

rst_values="0 1"
hl_values=$(seq 0 10)     # 11 architectures (indices into the list in the python script)
alp_values=$(seq 0 3)     # 4 alpha values
lri_values=$(seq 0 4)     # 5 learning rates

njobs=0
for rst in ${rst_values}; do
  for hl in ${hl_values}; do
    for alp in ${alp_values}; do
      for lri in ${lri_values}; do
        tag="${rst}_${hl}_${alp}_${lri}"
        export rst=${rst}
        export hl=${hl}
        export alp=${alp}
        export lri=${lri}
        sbatch --output="${DIRJOUT}/joptimize_NNMLP_${tag}.out" \
               --export=ALL,rst=${rst},hl=${hl},alp=${alp},lri=${lri} \
               "${DIRSCRIPTS}/joptimize_NNMLP.sh"

        njobs=$((njobs+1))
        echo " job joptimize_NNMLP_${tag} submitted OK at $(date +"%T")"
        sleep 1
      done
    done
  done
done
echo "${njobs} jobs submitted"

