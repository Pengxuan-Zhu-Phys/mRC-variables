#!/bin/bash

# get procid from command line
proc=$1
procid=$2
JobId=$_CONDOR_IHEP_JOB_ID

# map 0,1,2,...,30 to 1,2,3,...,31
sub_name_number=`expr $procid + 1`

# format 1,2,3,...,31 to 01,02,03,...,31
sub_name=`printf "%05d\n" $sub_name_number` 

# run the real job script by the formatted file name
bdir=/publicfs/atlas/atlasnew/SUSY/users/yuanjiarong/CEPC240/bkg
cd ${bdir}/"$proc"/
mkdir -p output
mkdir -p Log/
mkdir -p Log/$JobId/
cd Log/$JobId/
rm -rf *
cp ${bdir}/"$proc"/run_"$proc"_"${sub_name}".sh .
cp ${bdir}/"$proc"/simu_"$proc"_"${sub_name}".macro .
cp ${bdir}/"$proc"/event_"$proc"_"${sub_name}".macro .
cp ${bdir}/"$proc"/reco_"$proc"_"${sub_name}".xml .

source run_"$proc"_"${sub_name}".sh

mv *.log ../
mv output/* ../../output/
cd ../
rm -rf $JobId
