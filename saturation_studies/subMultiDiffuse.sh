for j in {1..1000}
do
    # if some jobs fail and you need to run again, this basically checks which jobs need to be re-run
    file=/data/user/HK_OD_work/saturation_studies/reduced_one_diffuse_barrel_9000k_1ev_${j}.root
    if [ -f "$file" ]; then
	echo "$file exists"
    else 
	cp /path/to/your/installation/of/WCSim/submitWCSim.sh /path/to/your/installation/of/WCSim/submitWCSim_${j}.sh
	cp /path/to/your/installation/of/WCSim/one_diffuse_fibre.mac /path/to/your/installation/of/WCSim/one_diffuse_fibre_${j}.mac

	sleep 4 # sometimes the Sheffield cluster would complain if I didn't let it rest between commands ¯\_(ツ)_/¯

	sed -i "s/seed 3254/seed 32545$j/" /path/to/your/installation/of/WCSim/one_diffuse_fibre_${j}.mac
	sed -i "s/one_diffuse_barrel_9000k_1ev.root/one_diffuse_barrel_9000k_1ev_${j}.root/" /path/to/your/installation/of/WCSim/one_diffuse_fibre_${j}.mac

	sed -i "s/one_diffuse_barrel_9000k_1ev.root/one_diffuse_barrel_9000k_1ev_${j}.root/" /path/to/your/installation/of/WCSim/submitWCSim_${j}.sh
	sed -i "s/reduced_one_diffuse_barrel_9000k_1ev.root/reduced_one_diffuse_barrel_9000k_1ev_${j}.root/" /path/to/your/installation/of/WCSim/submitWCSim_${j}.sh
	sed -i "s/one_diffuse_fibre.mac/one_diffuse_fibre_${j}.mac/" /path/to/your/installation/of/WCSim/submitWCSim_${j}.sh
	sed -i "s/default_diffuse.err/default_diffuse_${j}.err/" /path/to/your/installation/of/WCSim/submitWCSim_${j}.sh
	sed -i "s/default_diffuse.out/default_diffuse_${j}.out/" /path/to/your/installation/of/WCSim/submitWCSim_${j}.sh

	condor_qsub /path/to/your/installation/of/WCSim/submitWCSim_${j}.sh -e /scratch/user/default_diffuse_${j}.err -o /scratch/user/default_diffuse_${j}.out
	fi
done;
