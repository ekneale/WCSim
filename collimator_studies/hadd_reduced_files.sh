for i in 0.975 0.9375 0.9 0.8625 0.825 0.7875 0.78 0.7725 0.765 0.7575 0.7425 0.735 0.7275 0.72 0.7125 0.675 0.6375 0.6 0.5625 0.525 #### rayff
do

    hadd -f /data/user/HK_OD_work/collimators/barrel_rayff/reduced_barrel_e100000_n10000_rayff_${i}.root /data/user/HK_OD_work/collimators/barrel_rayff/reduced_barrel_e10000_n10000_rayff_${i}_*.root
done;

for i in 1.69 1.625 1.56 1.495 1.43 1.365 1.352 1.339 1.326 1.313 1.287 1.274 1.261 1.248 1.235 1.17 1.105 1.04 0.975 0.91 #### abwff
do

    hadd /data/user/HK_OD_work/collimators/barrel_abwff/reduced_barrel_e100000_n10000_abwff_${i}.root /data/user/HK_OD_work/collimators/barrel_abwff/reduced_barrel_e10000_n10000_abwff_${i}_*.root
done;

for i in 0.63 0.72 0.765 0.81 0.819 0.828 0.837 0.846 0.855 0.8595 0.864 0.8685 0.873 0.8775 0.882 0.8865 0.891 0.8955 #### refl
do

    hadd -f /data/user/HK_OD_work/collimators/barrel_refl/reduced_barrel_e100000_n10000_refl_${i}.root /data/user/HK_OD_work/collimators/barrel_refl/reduced_barrel_e10000_n10000_refl_${i}_*.root
done;
