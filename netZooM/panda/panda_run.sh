#!/bin/bash
# Working script for running PANDA in different scenarios.
# Edit 'panda_config.m' to set program parameters.

# Console running
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('panda_config.m'); run('panda_run.m');"

# Foreground running
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('panda_config.m'); run('panda_run.m'); quit;"

# Background running (./panda_run.sh &)
nohup matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('panda_config.m'); run('panda_run.m'); quit;" >& panda.`hostname`.log

# Email notification when done
echo "PANDA run on `hostname` has just finished: `date`." | mail -s "Task finished on `hostname`" -a panda.`hostname`.log `whoami`

# Testing
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('panda_config_test.m'); run('panda_run.m'); quit;"
#diff /tmp/panda.test.txt test_data/panda.txt

# Benchmark (5 repeats)
#matlab -nodisplay -nosplash -nodesktop -nojvm -r "run('panda_config_test.m'); run('panda_run.m'); run('panda_run.m'); run('panda_run.m'); run('panda_run.m'); run('panda_run.m'); quit;" >& panda.test.`hostname`.log
