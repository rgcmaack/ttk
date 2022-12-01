#!/bin/bash

logpath="log"
computeModes=("mssfast" "mssfine" "msssplit" "msswall")
maxextend=512
runthreads=(24 20 16 12 8 4 1)
maxthreads=24

for t in ${runthreads[@]}; do
    if [[ $maxthreads -lt ${runthreads[$t]} ]]
    then
        maxthreads=${runthreads[$t]}
    fi
done

echo $maxthreads

for m in ${computeModes[@]}; do
    for t in ${runthreads[@]}; do
        logfilepath="${logpath}/wsm_${maxextend}_${m}_${t}T.log"
        echo $logfilepath
        echo weakScalingMiranda.py ${m} ${maxextend} ${t} ${maxthreads} | tee $logfilepath
        pvpython weakScalingMiranda.py ${m} ${maxextend} ${t} ${maxthreads} | tee -a $logfilepath
        grep processed $logfilepath
    done
done

exit 0