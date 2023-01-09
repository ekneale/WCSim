count=0

for i in 0.8955 0.891 0.8865 0.882 0.8775 0.873 0.8685 0.864 0.8595 0.855 0.846 0.837 0.828 0.819 0.81 0.765 0.72 0.63 #### refl

do
    count=$((count+1))
    echo $count

	cp /path/to/your/installation/of/WCSim/local_water_params.mac /path/to/your/installation/of/WCSim/local_water_params_$count.mac

	sleep 2

	sed -i "s/WCODTyvekReflectivity\s0.90/WCODTyvekReflectivity $i/" /path/to/your/installation/of/WCSim/local_water_params_$count.mac

	for j in {1..1000}
	do
	    file=/data/user/HK_OD_work/collimators/barrel_refl/reduced_barrel_e10000_n10000_refl_${i}_${j}.root
	    if [ -f "$file" ]; then
		echo "$file exists"
	    else
		cp /path/to/your/installation/of/WCSim/submitWaterParams.sh /path/to/your/installation/of/WCSim/submitWaterParamsre_${count}_${j}.sh
		cp /path/to/your/installation/of/WCSim/water_properties_testing.mac /path/to/your/installation/of/WCSim/wpt_barrel_re_${count}_${j}.mac

		sleep 2

		sed -i "s/barrel_default/barrel_refl/" /path/to/your/installation/of/WCSim/submitWaterParamsre_${count}_${j}.sh

		sed -i "s/seed 3254/seed 3254$j/" /path/to/your/installation/of/WCSim/wpt_barrel_re_${count}_${j}.mac
		sed -i "s/barrel_e10000_n10000.root/barrel_e10000_n10000_refl_${i}_${j}.root/" /path/to/your/installation/of/WCSim/wpt_barrel_re_${count}_${j}.mac

		sed -i "s/barrel_e10000_n10000.root/barrel_e10000_n10000_refl_${i}_${j}.root/" /path/to/your/installation/of/WCSim/submitWaterParamsre_${count}_${j}.sh
		sed -i "s/reduced_barrel_e10000_n10000.root/reduced_barrel_e10000_n10000_refl_${i}_${j}.root/" /path/to/your/installation/of/WCSim/submitWaterParamsre_${count}_${j}.sh
		sed -i "s/wpt_barrel.mac/wpt_barrel_re_${count}_${j}.mac/" /path/to/your/installation/of/WCSim/submitWaterParamsre_${count}_${j}.sh
		sed -i "s/local_water_params.mac/local_water_params_$count.mac/" /path/to/your/installation/of/WCSim/submitWaterParamsre_${count}_${j}.sh

		sed -i "s/default_barrel.err/refl_barrel_${count}_${j}.err/" /path/to/your/installation/of/WCSim/submitWaterParamsre_${count}_${j}.sh
		sed -i "s/default_barrel.out/refl_barrel_${count}_${j}.out/" /path/to/your/installation/of/WCSim/submitWaterParamsre_${count}_${j}.sh

		condor_qsub /path/to/your/installation/of/WCSim/submitWaterParamsre_${count}_${j}.sh  -e /scratch/user/refl_barrel_${count}_${j}.err -o /scratch/user/refl_barrel_${count}_${j}.out
	    fi
	done;
 #   fi
done;
