#!/bin/bash

mkdir -p output

source /cvmfs/cepc.ihep.ac.cn/software/cepcenv/setup.sh
cepcenv -r /cvmfs/cepc.ihep.ac.cn/software/cepcsoft use 0.1.0-rc9

echo "time before simu ", `date`
rm -f simu_SB1_SB2.log; (time Mokka -U simu_SB1_SB2.macro) | cat > simu_SB1_SB2.log 
echo "time after running", `date`

echo "time before recon ", `date`
rm -f reco_SB1_SB2.log; (time Marlin reco_SB1_SB2.xml) | cat > reco_SB1_SB2.log;
 #rm -f Gear*.xml
echo "time after running", `date`
