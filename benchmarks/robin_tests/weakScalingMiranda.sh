#!/bin/bash

logpath="/data/log"
computeModes=("msswall" "mssfast" "mssfine" "msssplit")
maxextend=256
runthreads=(1 2 4 6)
maxthreads=0

for t in ${runthreads[@]}; do
    if [[ $maxthreads -lt ${runthreads[$t]} ]]
    then
        maxthreads=${runthreads[$t]}
    fi
done

for m in ${computeModes[@]}; do
    for t in ${runthreads[@]}; do
        logfilepath="${logpath}/wsm_${m}_${t}.log"
        echo $logfilepath
        pvpython weakScalingMiranda.py ${m} ${maxextend} ${t} ${maxthreads} >> $logfilepath
        sleep 10
    done
done

exit 0