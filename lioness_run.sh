#!/bin/bash
# Working script for running LIONESS in different scenarios.
# Must run PANDA first to get preprocessed middle files and PANDA network.
# Edit 'lioness_config.m' to set parameters.

# Console running
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('lioness_config.m'); run('lioness_run.m');"

# Foreground running
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('lioness_config.m'); run('lioness_run.m'); quit;"

# Background running (e.g., ./lioness_run.sh &)
#nohup matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('lioness_config.m'); run('lioness_run.m'); quit;" >& lioness.`hostname`.log

# Email notification when done
#echo "LIONESS run on `hostname` has just finished: `date`." | mail -s "Task finished on `hostname`" -a lioness.`hostname`.log `whoami`

# Example: Batch running on different machines
# Usage: Run this script by specifying the sample range, e.g., ./lioness_run.sh 1 100
if [[ $# -eq 0 ]] ; then
    echo 'Must assign a range of samples to run! For example, ./lioness_run.sh 1 100'
    exit 1
fi
start=`date`
host=`hostname`
nohup matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('lioness_config.m'); START=str2num('$1'); END=str2num('$2'); run('lioness_run.m'); quit;" >& lioness.$host.$1-$2.log
echo "LIONESS run on $host for sample $1 - $2 starting $start ends `date`." | mail -s "Task finished on $host" -a lioness.$host.$1-$2.log `whoami`
