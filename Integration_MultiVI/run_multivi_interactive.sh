#!/usr/bin/env bash

set -euxo pipefail

# JOB CONFIG
FARM_GROUP="teichlab"
FARM_QUEUE="gpu-basement"
RAM=32000
CPU=6
GPU_RAM=6000

# COMMAND CONFIG
NOTEBOOK_ROOT_FOLDER="/nfs/team205/kk18/data/6region_v2/MultiVI/notebooks"
PASSWORD="123"
PORT="1555"
CMD="jupyter notebook --notebook-dir=${NOTEBOOK_ROOT_FOLDER} --NotebookApp.token=${PASSWORD} --ip=0.0.0.0 --port=${PORT} --no-browser --allow-root"

bsub -G $FARM_GROUP \
    -q $FARM_QUEUE \
    "-n${CPU}" \
    "-M${RAM}" -R"select[mem>${RAM}] rusage[mem=${RAM}] span[hosts=1]" \
    -gpu"mode=shared:j_exclusive=no:gmem=${GPU_RAM}:num=1" -Is \
    /software/singularity-v3.6.4/bin/singularity exec --nv --bind /nfs,/lustre \
      /nfs/cellgeni/singularity/images/multivi0.13.0_cuda10.2_pytorch1.8.1_py3.8.sif \
      $CMD
