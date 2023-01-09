#!/bin/sh

source /path/to/your/WCSim/environment_setup_script.sh

SCRATCH_DIR=/scratch/user/ # essentially whatever temp/scratch area you want your jobs to run in
if [ ! -e $SCRATCH_DIR ]; then mkdir $SCRATCH_DIR; fi

BASE_DATA_DIR=/data/user/HK_OD_work/ # wherever you want to save your OD work
if [ ! -e $BASE_DATA_DIR ]; then mkdir $BASE_DATA_DIR; fi

DATA_DIR=$BASE_DATA_DIR/saturation_studies/
if [ ! -e $DATA_DIR ]; then mkdir $DATA_DIR; fi

LOG_DIR=$DATA_DIR/logs/
if [ ! -e $LOG_DIR ]; then mkdir $LOG_DIR; fi

JOB_DIR=$SCRATCH_DIR
if [ ! -e $JOB_DIR ]; then mkdir $JOB_DIR; fi

(>&2 echo `hostname`)
cd /path/to/your/installation/of/WCSim/

WCSim one_diffuse_fibre.mac

root -l -x -b -q 'ReduceODOutput.C("/scratch/user/one_diffuse_barrel_9000k_1ev.root","/scratch/user/reduced_one_diffuse_barrel_9000k_1ev.root")'

rm /scratch/user/one_diffuse_barrel_9000k_1ev.root
mv /scratch/user/reduced_one_diffuse_barrel_9000k_1ev.root $DATA_DIR
mv /scratch/user/default_diffuse.err $LOG_DIR
mv /scratch/user/default_diffuse.out $LOG_DIR


