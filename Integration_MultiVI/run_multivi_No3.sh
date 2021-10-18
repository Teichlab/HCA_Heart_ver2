#!/usr/bin/env bash

set -euxo pipefail

# JOB CONFIG
FARM_GROUP="teichlab"
FARM_QUEUE="gpu-basement"
RAM=32000
CPU=24
GPU_RAM=32000

SCRIPT_PATH=/lustre/scratch117/cellgen/team205/dm19/kazu/notebooks/Integration_MultiVI/24CPUs_full
SCRIPT=No3_run_MultiVI.py


bsub -G $FARM_GROUP \
     -q $FARM_QUEUE \
     -o job_%J_output.txt \
     -e job_%J_error.txt \
    "-n${CPU}" \
    "-M${RAM}" -R"select[mem>${RAM}] rusage[mem=${RAM}] span[hosts=1]" \
    -gpu"mode=shared:j_exclusive=yes:gmem=${GPU_RAM}:num=1" \
    /software/singularity-v3.6.4/bin/singularity exec --nv --bind /nfs,/lustre \
    /nfs/cellgeni/singularity/images/multivi0.13.0_cuda10.2_pytorch1.8.1_py3.8.sif \
    python ${SCRIPT_PATH}/${SCRIPT}
