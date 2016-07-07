#!/bin/bash
#
# Run this file using 'nohup bash RunPUMA.sh &'
#
#$ -S /bin/bash
#$ -cwd
#$ -N RunPUMALIONESS
#$ -j y
#$ -o RunPUMALIONESS.qlog

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
