#!/bin/bash
#
# Move this file to the directory of RunPUMA.m
# Then run this file using 'nohup bash RunPUMA.sh &'
#
#$ -S /bin/bash
#$ -cwd
#$ -N RunPUMA
#$ -j y
#$ -o RunPUMA.qlog

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "Task index number : $SGE_TASK_ID"
echo "=========================================================="

matlab -nodisplay -nodesktop -nojvm -nosplash < RunPUMA.m > RunPUMA.qlog

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="
