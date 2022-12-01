#!/bin/bash

logpath="log"
computeModes=("msswall" "mssfast" "mssfine" "msssplit")
extend=512
runthreads=(1 2 4 8 12 16 20 24)

for m in ${computeModes[@]}; do
    for t in ${runthreads[@]}; do
        logfilepath="${logpath}/ssm_${m}_${t}T.log"
        echo $logfilepath
        echo miranda.py ${m} ${extend} ${t} | tee $logfilepath
        pvpython miranda.py ${m} ${extend} ${t} | tee -a $logfilepath
        sleep 1
    done
done

exit 0