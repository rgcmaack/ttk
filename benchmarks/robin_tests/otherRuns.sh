#!/bin/bash

logpath="/data/log"
computeModes3D=("msswall" "mssfast" "mssfine" "msssplit")
computeModes2D=("mss" "msc")
threads=6

for m in ${computeModes2D[@]}; do
    logfilepath="${logpath}/hugeNoisyTerrain_${m}.log"
    echo $logfilepath
    pvpython hugeNoisyTerrain.py ${m} ${threads} >> $logfilepath
    sleep 10
done

for m in ${computeModes2D[@]}; do
    logfilepath="${logpath}/noisyTerrain_${m}.log"
    echo $logfilepath
    pvpython noisyTerrain.py ${m} ${threads} >> $logfilepath
    sleep 10
done

for m in ${computeModes3D[@]}; do
    logfilepath="${logpath}/viscous_${m}_${threads}.log"
    echo $logfilepath
    pvpython viscousFingering.py ${m} ${threads} >> $logfilepath
    sleep 10
done

exit 0