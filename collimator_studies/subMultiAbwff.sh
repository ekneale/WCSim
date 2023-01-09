count=0

for i in 1.69 1.625 1.56 1.495 1.43 1.365 1.352 1.339 1.326 1.313 1.287 1.274 1.261 1.248 1.235 1.17 1.105 1.04 0.975 0.91 #### abwff

do
    count=$((count+1))
    echo $count

	cp /path/to/your/installation/of/WCSim/local_water_params.mac /path/to/your/installation/of/WCSim/local_water_params_$count.mac

	sleep 2

	sed -i "s/abwff\s1.3/abwff $i/" /path/to/your/installation/of/WCSim/local_water_params_$count.mac

	for j in {1..1000}
	do
	    file=/data/user/HK_OD_work/collimators/barrel_abwff/reduced_barrel_e10000_n10000_abwff_${i}_${j}.root
	    if [ -f "$file" ]; then
		echo "$file exists"
	    else
		cp /path/to/your/installation/of/WCSim/submitWaterParams.sh /path/to/your/installation/of/WCSim/submitWaterParamsa_${count}_${j}.sh
		cp /path/to/your/installation/of/WCSim/water_properties_testing.mac /path/to/your/installation/of/WCSim/wpt_barrel_a_${count}_${j}.mac

		sleep 2

		sed -i "s/barrel_default/barrel_abwff/" /path/to/your/installation/of/WCSim/submitWaterParamsa_${count}_${j}.sh

		sed -i "s/seed 3254/seed 3254$j/" /path/to/your/installation/of/WCSim/wpt_barrel_a_${count}_${j}.mac
		sed -i "s/barrel_e10000_n10000.root/barrel_e10000_n10000_abwff_${i}_${j}.root/" /path/to/your/installation/of/WCSim/wpt_barrel_a_${count}_${j}.mac

		sed -i "s/barrel_e10000_n10000.root/barrel_e10000_n10000_abwff_${i}_${j}.root/" /path/to/your/installation/of/WCSim/submitWaterParamsa_${count}_${j}.sh
		sed -i "s/reduced_barrel_e10000_n10000.root/reduced_barrel_e10000_n10000_abwff_${i}_${j}.root/" /path/to/your/installation/of/WCSim/submitWaterParamsa_${count}_${j}.sh
		sed -i "s/wpt_barrel.mac/wpt_barrel_a_${count}_${j}.mac/" /path/to/your/installation/of/WCSim/submitWaterParamsa_${count}_${j}.sh
		sed -i "s/local_water_params.mac/local_water_params_$count.mac/" /path/to/your/installation/of/WCSim/submitWaterParamsa_${count}_${j}.sh

		sed -i "s/default_barrel.err/abwff_barrel_${count}_${j}.err/" /path/to/your/installation/of/WCSim/submitWaterParamsa_${count}_${j}.sh
		sed -i "s/default_barrel.out/abwff_barrel_${count}_${j}.out/" /path/to/your/installation/of/WCSim/submitWaterParamsa_${count}_${j}.sh

		condor_qsub /path/to/your/installation/of/WCSim/submitWaterParamsa_${count}_${j}.sh  -e /scratch/user/abwff_barrel_${count}_${j}.err -o /scratch/user/abwff_barrel_${count}_${j}.out
	    fi
	done;
 #   fi
done;
