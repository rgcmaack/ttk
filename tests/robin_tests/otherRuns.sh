#!/bin/bash

logpath="log"
computeModes3D=("msswall" "mssfast" "mssfine" "msssplit", "msc")
computeModes2D=("mss" "msc")
# threads=24

for m in ${computeModes2D[@]}; do
    for threads in 1 4 8 12 16 20 24; do
        logfilepath="${logpath}/hugeNoisyTerrain_${m}_${threads}T.log"
        echo $logfilepath
        echo hugeNoisyTerrain ${m} ${threads} > $logfilepath
        pvpython hugeNoisyTerrain.py ${m} ${threads} >> $logfilepath
        grep "processed" $logfilepath
        sleep 1
    done
done

for m in ${computeModes2D[@]}; do
    for threads in 1 4 8 12 16 20 24; do
        logfilepath="${logpath}/noisyTerrain_${m}_${threads}T.log"
        echo $logfilepath
        echo noisyTerrain ${m} ${threads} > $logfilepath
        pvpython noisyTerrain.py ${m} ${threads} >> $logfilepath
        grep "processed" $logfilepath
        sleep 1
    done
done

for m in ${computeModes3D[@]}; do
    for threads in 1 4 8 12 16 20 24; do
        logfilepath="${logpath}/viscous_${m}_${threads}T.log"
        echo $logfilepath
        echo viscousFingering ${m} ${threads} > $logfilepath
        pvpython viscousFingering.py ${m} ${threads} >> $logfilepath
        grep "processed" $logfilepath
        sleep 1
    done
done

exit 0