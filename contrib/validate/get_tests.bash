#!/bin/bash
#
# get_tests.bash
# --------------
#
# Go through the output of the QA run and extract the names of all the test
# cases that were run. The list of tests run is echoed to standard output.
#
# Arguments:
#
# - $0: The name of this script
# - $1: The log filename from the QA test suite run
#
job_list=`grep -v 'QM: Running' $1 | grep Running`
for job in $job_list; do
   if [ $job != "Running" ] ; then
     testname=`basename $job`
     echo $testname
   fi
done
