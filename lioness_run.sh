#!/bin/bash
# Working script for running LIONESS in different scenarios.
# Must run PANDA first to get preprocessed middle files and PANDA network.
# Edit 'lioness_config.m' to set parameters.

# Console running
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('lioness_config.m'); run('lioness_run.m');"

# Foreground running
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('lioness_config.m'); run('lioness_run.m'); quit;"

# Background running
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('lioness_config.m'); run('lioness_run.m'); quit;" > lioness.`hostname`.log &

# Email notification when done
#echo "LIONESS run on `hostname` has just finished: `date`." | mail -s "Task finished on `hostname`" `whoami`

# Example: Batch running on different machines
# Usage: 1) Set up a range of samples in the corresponding config file.
#        2) Run this script by specifying the range, e.g., ./lioness_run.sh 1 100
if [[ $# -eq 0 ]] ; then
    echo 'Must assign a range! See the configuration m-file.'
    exit 1
fi
start=`date`
host=`hostname`
echo "LIONESS starts! Logging: lioness.$host.$1-$2.log. Date: $start"
matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('lioness_config_$host.m'); run('lioness_run.m'); quit;" > "lioness.$host.$1-$2.log" &
echo "LIONESS run on $host for sample $1 - $2 starting $start ends `date`." | mail -s "Task finished on $host" `whoami`
