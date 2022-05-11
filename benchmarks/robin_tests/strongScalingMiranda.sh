#!/bin/bash

logpath="/data/log"
computeModes=("msswall" "mssfast" "mssfine" "msssplit")
extend=256
runthreads=(1 2 4 6)

for m in ${computeModes[@]}; do
    for t in ${runthreads[@]}; do
        logfilepath="${logpath}/ssm_${m}_${t}.log"
        echo $logfilepath
        pvpython miranda.py ${m} ${extend} ${t} >> $logfilepath
        sleep 10
    done
done

exit 0