for j in {1..1000}
do
    file=/data2/pidcott/HKOD/2021merge/Collimators/SDO/barrel_default/reduced_barrel_e10000_n10000_${j}.root
    if [ -f "$file" ]; then
	echo "$file exists"
    else
	cp /home/pidcott/OD_calib/calib_devel/WCSim/submitWaterParams.sh /home/pidcott/OD_calib/calib_devel/WCSim/submitWaterParams_${j}.sh
	cp /home/pidcott/OD_calib/calib_devel/WCSim/water_properties_testing.mac /home/pidcott/OD_calib/calib_devel/WCSim/wpt_barrel_${j}.mac

	sleep 4

	sed -i "s/seed 3254/seed 3254$j/" /home/pidcott/OD_calib/calib_devel/WCSim/wpt_barrel_${j}.mac
	sed -i "s/barrel_e10000_n10000.root/barrel_e10000_n10000_${j}.root/" /home/pidcott/OD_calib/calib_devel/WCSim/wpt_barrel_${j}.mac

	#sed -i "s/barrel_rayff/barrel_default/" /home/pidcott/OD_calib/calib_devel/WCSim/submitWaterParams_${j}.sh


	sed -i "s/barrel_e10000_n10000.root/barrel_e10000_n10000_${j}.root/" /home/pidcott/OD_calib/calib_devel/WCSim/submitWaterParams_${j}.sh
	sed -i "s/reduced_barrel_e10000_n10000.root/reduced_barrel_e10000_n10000_${j}.root/" /home/pidcott/OD_calib/calib_devel/WCSim/submitWaterParams_${j}.sh
	sed -i "s/wpt_barrel.mac/wpt_barrel_${j}.mac/" /home/pidcott/OD_calib/calib_devel/WCSim/submitWaterParams_${j}.sh
	sed -i "s/default_barrel.err/default_barrel_${j}.err/" /home/pidcott/OD_calib/calib_devel/WCSim/submitWaterParams_${j}.sh
	sed -i "s/default_barrel.out/default_barrel_${j}.out/" /home/pidcott/OD_calib/calib_devel/WCSim/submitWaterParams_${j}.sh
	#sed -i "s/local_water_params.mac/local_water_params_$count.mac/" /home/pidcott/OD_calib/calib_devel/WCSim/submitWaterParams_${count}_${j}.sh

	condor_qsub /home/pidcott/OD_calib/calib_devel/WCSim/submitWaterParams_${j}.sh -e /scratch/pidcott/default_barrel_${j}.err -o /scratch/pidcott/default_barrel_${j}.out
	fi
done;
