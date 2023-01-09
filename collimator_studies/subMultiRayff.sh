count=0

for i in 0.975 0.9375 0.9 0.8625 0.825 0.7875 0.78 0.7725 0.765 0.7575 0.7425 0.735 0.7275 0.72 0.7125 0.675 0.6375 0.6 0.5625 0.525 #### rayff

do
    count=$((count+1))
    echo $count

	cp /path/to/your/installation/of/WCSim/local_water_params.mac /path/to/your/installation/of/WCSim/local_water_params_$count.mac

	sleep 2

	sed -i "s/rayff\s0.75/rayff $i/" /path/to/your/installation/of/WCSim/local_water_params_$count.mac

	for j in {1..1000}
	do
	    file=/data/user/HK_OD_work/collimators/barrel_rayff/reduced_barrel_e10000_n10000_rayff_${i}_${j}.root
	    if [ -f "$file" ]; then
		echo "$file exists"
	    else
		cp /path/to/your/installation/of/WCSim/submitWaterParams.sh /path/to/your/installation/of/WCSim/submitWaterParamsra_${count}_${j}.sh
		cp /path/to/your/installation/of/WCSim/water_properties_testing.mac /path/to/your/installation/of/WCSim/wpt_barrel_ra_${count}_${j}.mac

		sleep 2

		sed -i "s/barrel_default/barrel_rayff/" /path/to/your/installation/of/WCSim/submitWaterParamsra_${count}_${j}.sh

		sed -i "s/seed 3254/seed 3254$j/" /path/to/your/installation/of/WCSim/wpt_barrel_ra_${count}_${j}.mac
		sed -i "s/barrel_e10000_n10000.root/barrel_e10000_n10000_rayff_${i}_${j}.root/" /path/to/your/installation/of/WCSim/wpt_barrel_ra_${count}_${j}.mac

		sed -i "s/barrel_e10000_n10000.root/barrel_e10000_n10000_rayff_${i}_${j}.root/" /path/to/your/installation/of/WCSim/submitWaterParamsra_${count}_${j}.sh
		sed -i "s/reduced_barrel_e10000_n10000.root/reduced_barrel_e10000_n10000_rayff_${i}_${j}.root/" /path/to/your/installation/of/WCSim/submitWaterParamsra_${count}_${j}.sh
		sed -i "s/wpt_barrel.mac/wpt_barrel_ra_${count}_${j}.mac/" /path/to/your/installation/of/WCSim/submitWaterParamsra_${count}_${j}.sh
		sed -i "s/local_water_params.mac/local_water_params_$count.mac/" /path/to/your/installation/of/WCSim/submitWaterParamsra_${count}_${j}.sh

		sed -i "s/default_barrel.err/rayff_barrel_${count}_${j}.err/" /path/to/your/installation/of/WCSim/submitWaterParamsra_${count}_${j}.sh
		sed -i "s/default_barrel.out/rayff_barrel_${count}_${j}.out/" /path/to/your/installation/of/WCSim/submitWaterParamsra_${count}_${j}.sh

		condor_qsub /path/to/your/installation/of/WCSim/submitWaterParamsra_${count}_${j}.sh  -e /scratch/user/rayff_barrel_${count}_${j}.err -o /scratch/user/rayff_barrel_${count}_${j}.out
	    fi
	done;
 #   fi
done;
