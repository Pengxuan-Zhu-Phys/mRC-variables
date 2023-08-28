#!/bin/bash

path=`pwd`
cd $path

#pars=("ww_l" "e3e3" "nnh_e3e3" "zzorww_l")

#par="sznu_l"
par="sw_l"

export SimuWorkDir=/publicfs/atlas/atlasnew/SUSY/users/yuanjiarong/CEPC240/bkg
export InputFilesDir=/cefs/data/stdhep/CEPC240/4fermions/E240.P${par}.e0.p0.whizard195
#export InputFilesDir=/cefs/data/stdhep/CEPC240/2fermions/E240.P${par}.e0.p0.whizard195
#export InputFilesDir=/cefs/data/stdhep/CEPC240/higgs/exclusive/E240.P${par}.e0.p0.whizard195

for file in ${InputFilesDir}/*.e0.p0.*.stdhep
do
	

	InputFiles_stdhep=${file}
#	echo $InputFiles_stdhep
	s=${file##*/}
	sjob=${s%%.e0*}
	ss=${s%%.stdhep}
	job=${ss##*p0.}
	
	OUTPUTDATA=${sjob}
#	OUTPUT_sim_slcio=$SimuWorkDir/$OUTPUTDATA/${sjob}_${job}.sim.slcio
#	OUTPUT_rec_slcio=$SimuWorkDir/$OUTPUTDATA/${sjob}_${job}.rec.slcio
#	OUTPUT_ana_root=$SimuWorkDir/$OUTPUTDATA/${sjob}_${job}.root
	OUTPUT_sim_slcio=${sjob}_${job}.sim.slcio
	OUTPUT_rec_slcio=${sjob}_${job}.rec.slcio
	OUTPUT_ana_root=${sjob}_${job}.root

	mkdir -p $SimuWorkDir/$OUTPUTDATA/
	cd $SimuWorkDir/$OUTPUTDATA/
#	echo "in dir:"
#	echo `pwd`

#	if [ -f "$SimuWorkDir/$OUTPUTDATA/filter_${sjob}_${job}.xml" ]; then
#	rm -f $SimuWorkDir/$OUTPUTDATA/filter_${sjob}_${job}.xml
#	continue
#	fi
	if [ -f "$OUTPUT_sim_slcio" ]; then
	rm -f $OUTPUT_sim_slcio
	continue
	fi
	if [ -f "$OUTPUT_rec_slcio" ]; then
	rm -f $OUTPUT_rec_slcio
	continue
	fi
	cp /workfs/bes/yuanjiarong/HepSub/CEPCstdhep/simu.macro $SimuWorkDir/$OUTPUTDATA/simu_${sjob}_${job}.macro
	cp /workfs/bes/yuanjiarong/HepSub/CEPCstdhep/event.macro $SimuWorkDir/$OUTPUTDATA/event_${sjob}_${job}.macro
	cp /workfs/bes/yuanjiarong/HepSub/CEPCstdhep/reco.xml $SimuWorkDir/$OUTPUTDATA/reco_${sjob}_${job}.xml
	cp /workfs/bes/yuanjiarong/HepSub/CEPCstdhep/run_all.sh $SimuWorkDir/$OUTPUTDATA/run_${sjob}_${job}.sh
       
	sed -i "s#STDHEPFILE#$file#g" $SimuWorkDir/$OUTPUTDATA/event_${sjob}_${job}.macro
	sed -i "s#SB1#${sjob}#g" $SimuWorkDir/$OUTPUTDATA/simu_${sjob}_${job}.macro
	sed -i "s#SB2#${job}#g" $SimuWorkDir/$OUTPUTDATA/simu_${sjob}_${job}.macro
	sed -i "s#SIMUFILE#$OUTPUT_sim_slcio#g" $SimuWorkDir/$OUTPUTDATA/simu_${sjob}_${job}.macro
	sed -i "s#SIMUFILE#$OUTPUT_sim_slcio#g" $SimuWorkDir/$OUTPUTDATA/reco_${sjob}_${job}.xml
	sed -i "s#RECOFILE#$OUTPUT_rec_slcio#g" $SimuWorkDir/$OUTPUTDATA/reco_${sjob}_${job}.xml
	sed -i "s#ROOTFILE#$OUTPUT_ana_root#g" $SimuWorkDir/$OUTPUTDATA/reco_${sjob}_${job}.xml
	sed -i "s#SB1#${sjob}#g" $SimuWorkDir/$OUTPUTDATA/run_${sjob}_${job}.sh
	sed -i "s#SB2#${job}#g" $SimuWorkDir/$OUTPUTDATA/run_${sjob}_${job}.sh
	
#	echo "$SimuWorkDir/$OUTPUTDATA/event_${sjob}_${job}.macro"
#	echo "$SimuWorkDir/$OUTPUTDATA/simu_${sjob}_${job}.macro"
#	echo "$SimuWorkDir/$OUTPUTDATA/reco_${sjob}_${job}.xml"
#	echo "$SimuWorkDir/$OUTPUTDATA/run_${sjob}_${job}.sh"
	chmod +x $SimuWorkDir/$OUTPUTDATA/run_${sjob}_${job}.sh
	chmod +x $SimuWorkDir/$OUTPUTDATA/event_${sjob}_${job}.macro
	chmod +x $SimuWorkDir/$OUTPUTDATA/simu_${sjob}_${job}.macro

#	heps $SimuWorkDir/$OUTPUTDATA/run_${sjob}_${job}.sh
#	sh $SimuWorkDir/$OUTPUTDATA/run_${sjob}_${job}.sh


done

cd /workfs/bes/yuanjiarong/HepSub/CEPCstdhep;

