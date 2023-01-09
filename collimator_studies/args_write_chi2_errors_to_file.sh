barrelORbottom=$1 # wants barrel, bottom, or top
rayffORabwff=$2   # wants rayff, abwff, or refl

cp reduced_chi_squared_errors.C ${1}_${2}_reduced_chi_squared_errors.C

sleep 4

sed -i "s/reduced_chi_squared_errors/${1}_${2}_reduced_chi_squared_errors/" ${1}_${2}_reduced_chi_squared_errors.C

if [ "$rayffORabwff" == "rayff" ]
    then
    for i in 0.525 0.5625 0.6 0.6375 0.675 0.7125 0.72 0.7275 0.735 0.7425 0.7575 0.765 0.7725 0.78 0.7875 0.825 0.8625 0.9 0.9375 0.975
    #### rayff
    do
	if [ "$barrelORbottom" == "barrel" ]
	    then
	    sed -i "s/cylLoc4_rayff_chi2_errors.txt/cylLoc4_rayff_chi2_errors.txt/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/rayff_cylLoc4_reduced_chi_sq.txt/rayff_cylLoc4_reduced_chi_sq.txt/" ${1}_${2}_reduced_chi_squared_errors.C
            root -l -x -b -q "${1}_${2}_reduced_chi_squared_errors.C(\"/data/user/HK_OD_work/collimators/barrel_rayff/reduced_barrel_e1000000_n10000_rayff_$i.root\", \"/data/user/HK_OD_work/collimators/barrel_default/reduced_barrel_e1000000_n10000.root\")"
	fi
	if [ "$barrelORbottom" == "bottom" ]
	    then
	    sed -i "s/6481/1527/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc==4/cylLoc==3/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc4_rayff_chi2_errors.txt/cylLoc3_rayff_chi2_errors.txt/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/rayff_cylLoc4_reduced_chi_sq.txt/rayff_cylLoc3_reduced_chi_sq.txt/" ${1}_${2}_reduced_chi_squared_errors.C
            root -l -x -b -q "${1}_${2}_reduced_chi_squared_errors.C(\"/data/user/HK_OD_work/collimators/bottom_rayff/reduced_bottom_e100000_n10000_rayff_$i.root\", \"/data/user/HK_OD_work/collimators/bottom_default/bottom_default.root\")"
	fi 
	if [ "$barrelORbottom" == "top" ]
	    then
	    sed -i "s/6481/1527/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc==4/cylLoc==5/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc4_rayff_chi2_errors.txt/cylLoc5_rayff_chi2_errors.txt/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/rayff_cylLoc4_reduced_chi_sq.txt/rayff_cylLoc5_reduced_chi_sq.txt/" ${1}_${2}_reduced_chi_squared_errors.C
            root -l -x -b -q "${1}_${2}_reduced_chi_squared_errors.C(\"/data/user/HK_OD_work/collimators/top_rayff/reduced_top_e100000_n10000_rayff_$i.root\", \"/data/user/HK_OD_work/collimators/top_default/top_default.root\")"
	fi 
    done;
fi


if [ "$rayffORabwff" == "abwff" ]
    then
    for i in 0.91 0.975 1.04 1.105 1.17 1.235 1.248 1.261 1.274 1.287 1.313 1.326 1.339 1.352 1.365 1.43 1.495 1.56 1.625 1.69
    #### abwff
    do
	if [ "$barrelORbottom" == "barrel" ]
	    then
	    sed -i "s/cylLoc4_rayff_chi2_errors.txt/cylLoc4_abwff_chi2_errors.txt/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/rayff_cylLoc4_reduced_chi_sq.txt/abwff_cylLoc4_reduced_chi_sq.txt/" ${1}_${2}_reduced_chi_squared_errors.C
            root -l -x -b -q "${1}_${2}_reduced_chi_squared_errors.C(\"/data/user/HK_OD_work/collimators/barrel_abwff/reduced_barrel_e1000000_n10000_abwff_$i.root\", \"/data/user/HK_OD_work/collimators/barrel_default/reduced_barrel_e1000000_n10000.root\")"
	fi
	if [ "$barrelORbottom" == "bottom" ]
	    then
	    sed -i "s/6481/1527/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc==4/cylLoc==3/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc4_rayff_chi2_errors.txt/cylLoc3_abwff_chi2_errors.txt/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/rayff_cylLoc4_reduced_chi_sq.txt/abwff_cylLoc3_reduced_chi_sq.txt/" ${1}_${2}_reduced_chi_squared_errors.C
            root -l -x -b -q "${1}_${2}_reduced_chi_squared_errors.C(\"/data/user/HK_OD_work/collimators/bottom_abwff/reduced_bottom_e100000_n10000_abwff_$i.root\", \"/data/user/HK_OD_work/collimators/bottom_default/bottom_default.root\")"
	fi
	if [ "$barrelORbottom" == "top" ]
	then
	    sed -i "s/6481/1527/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc==4/cylLoc==5/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc4_rayff_chi2_errors.txt/cylLoc5_abwff_chi2_errors.txt/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/rayff_cylLoc4_reduced_chi_sq.txt/abwff_cylLoc5_reduced_chi_sq.txt/" ${1}_${2}_reduced_chi_squared_errors.C
            root -l -x -b -q "${1}_${2}_reduced_chi_squared_errors.C(\"/data/user/HK_OD_work/collimators/top_abwff/reduced_top_e100000_n10000_abwff_$i.root\", \"/data/user/HK_OD_work/collimators/top_default/top_default.root\")"
	fi 

    done;
fi

if [ "$rayffORabwff" == "refl" ]
    then
    for i in 0.63 0.72 0.765 0.81 0.819 0.828 0.837 0.846 0.855 0.8595 0.864 0.8685 0.873 0.8775 0.882 0.8865 0.891 0.8955
    #### refl
    do
	if [ "$barrelORbottom" == "barrel" ]
	    then
	    sed -i "s/cylLoc4_rayff_chi2_errors.txt/cylLoc4_refl_chi2_errors.txt/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/rayff_cylLoc4_reduced_chi_sq.txt/refl_cylLoc4_reduced_chi_sq.txt/" ${1}_${2}_reduced_chi_squared_errors.C
            root -l -x -b -q "${1}_${2}_reduced_chi_squared_errors.C(\"/data/user/HK_OD_work/collimators/barrel_refl/reduced_barrel_e1000000_n10000_refl_$i.root\", \"/data2/pidcott/HKOD/2021merge/Collimators/barrel_default/reduced_barrel_e1000000_n10000.root\")"
	fi
	if [ "$barrelORbottom" == "bottom" ]
	    then
	    sed -i "s/6481/1527/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc==4/cylLoc==3/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc4_rayff_chi2_errors.txt/cylLoc3_refl_chi2_errors.txt/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/rayff_cylLoc4_reduced_chi_sq.txt/refl_cylLoc3_reduced_chi_sq.txt/" ${1}_${2}_reduced_chi_squared_errors.C
            root -l -x -b -q "${1}_${2}_reduced_chi_squared_errors.C(\"/data/user/HK_OD_work/collimators/bottom_refl/reduced_bottom_e100000_n10000_refl_$i.root\", \"/data/user/HK_OD_work/collimators/bottom_default/bottom_default.root\")"
	fi
	if [ "$barrelORbottom" == "top" ]
	then
	    sed -i "s/6481/1527/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc==4/cylLoc==5/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/cylLoc4_rayff_chi2_errors.txt/cylLoc5_refl_chi2_errors.txt/" ${1}_${2}_reduced_chi_squared_errors.C
	    sed -i "s/rayff_cylLoc4_reduced_chi_sq.txt/refl_cylLoc5_reduced_chi_sq.txt/" ${1}_${2}_reduced_chi_squared_errors.C
            root -l -x -b -q "${1}_${2}_reduced_chi_squared_errors.C(\"/data/user/HK_OD_work/collimators/top_refl/reduced_top_e100000_n10000_refl_$i.root\", \"/data/user/HK_OD_work/collimators/top_default/top_default.root\")"
	fi 

    done;
fi


rm ${1}_${2}_reduced_chi_squared_errors.C
