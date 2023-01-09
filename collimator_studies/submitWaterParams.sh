#!/bin/sh

source /path/to/your/WCSim/environment_setup_script.sh

SCRATCH_DIR=/scratch/user/ # essentially whatever temp/scratch area you want your jobs to run in
if [ ! -e $SCRATCH_DIR ]; then mkdir $SCRATCH_DIR; fi

BASE_DATA_DIR=/data/user/HK_OD_work/collimators/ # wherever you want to save your OD work
if [ ! -e $BASE_DATA_DIR ]; then mkdir $BASE_DATA_DIR; fi

DATA_DIR=$BASE_DATA_DIR/barrel_default/
if [ ! -e $DATA_DIR ]; then mkdir $DATA_DIR; fi

LOG_DIR=$DATA_DIR/logs/
if [ ! -e $LOG_DIR ]; then mkdir $LOG_DIR; fi

JOB_DIR=$SCRATCH_DIR
if [ ! -e $JOB_DIR ]; then mkdir $JOB_DIR; fi

(>&2 echo `hostname`)
cd /path/to/your/installation/of/WCSim/

WCSim wpt_barrel.mac local_water_params.mac

root -l -x -b -q 'ReduceODOutput.C("/scratch/user/barrel_e10000_n10000.root","/scratch/user/reduced_barrel_e10000_n10000.root")'

rm /scratch/user/barrel_e10000_n10000.root
mv /scratch/user/reduced_barrel_e10000_n10000.root $DATA_DIR
mv /scratch/user/default_barrel.err $LOG_DIR
mv /scratch/user/default_barrel.out $LOG_DIR

rm wpt_barrel.mac
#rm local_water_params.mac # this is used by multiple jobs so you'll just need to remove manually once all have run
rm -- "$0"

