#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --mem-per-cpu=1500

echo 'working directory ='
echo $PWD

runnum=$1
adccut=$2

# analyzer -b -q '/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/xtalk_replay/replay_xtalk.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments')' 


#echo 'swif_job_work_dir='$SWIF_JOB_WORK_DIR

# source login stuff since swif2 completely nukes any sensible default software environment
#source /home/puckett/.cshrc
#source /home/puckett/.login

#echo 'before sourcing environment, env = '
#env 
# source /site/12gev_phys/softenv.sh 2.4

# ldd /work/halla/sbs/ANALYZER/install/bin/analyzer |& grep not

# MODULES=/etc/profile.d/modules.sh 

# if [[ $(type -t module) != function && -r ${MODULES} ]]; then 
# source ${MODULES} 
# fi 

# if [ -d /apps/modulefiles ]; then y
# module use /apps/modulefiles 
# fi 

# module load gcc/9.2.0 

# source /site/12gev_phys/softenv.sh 2.4

# cd /w/halla-scshelf2102/sbs/jboyd/SBS_OFFLINE/install/run_replay_here

# ldd /work/halla/sbs/ANALYZER/install/bin/analyzer |& grep not
# ldd /work/halla/sbs/jboyd/analyzer/bin |& grep not

# export ROOTSYS=/site/12gev_phys/2.5/Linux_CentOS7.7.1908-gcc9.2.0/root/6.24.06

# source $ROOTSYS/bin/thisroot.sh

source /site/12gev_phys/softenv.sh 2.4

ldd /work/halla/sbs/ANALYZER/install/bin/analyzer |& grep not

echo 'working directory = '$PWD

export ANALYZER=/work/halla/sbs/ANALYZER/install
# export ANALYZER=/w/halla-scshelf2102/sbs/jboyd/ANALYZER/install
source $ANALYZER/bin/setup.sh

#source $ANALYZER/build/cmake/setup.sh

export SBS_ONLINE=/w/halla-scshelf2102/sbs/SBS_OFFLINE/SBS-offline
# export	SBS_ONLINE=/w/halla-scshelf2102/sbs/jboyd/SBS_OFFLINE/SBS-offline
# source /w/halla-scshelf2102/sbs/jboyd/SBS_OFFLINE/install/bin/sbsenv.sh

# source /w/halla-scshelf2102/sbs/jboyd/SBS_OFFLINE/build/sbsenv.sh
# source /w/halla-scshelf2102/sbs/jboyd/SBS_OFFLINE/install/bin/jb_sbsenv.sh
source /w/halla-scshelf2102/sbs/SBS_OFFLINE/install/bin/sbsenv.sh
# source /w/halla-scshelf2102/sbs/jboyd/analysis/jb_sbsenv.sh

export SBS_REPLAY=/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay
# export SBS_REPLAY=/w/halla-scshelf2102/sbs/jboyd/SBS_REPLAY/SBS-replay/
export DB_DIR=/w/halla-scshelf2102/sbs/jboyd/SBS_REPLAY/SBS-replay/DB
# export DATA_DIR=/lustre19/enp/cache/mss/halla/sbs/raw
# export DATA_DIR=/w/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd
export DATA_DIR=/lustre19/expphy/cache/halla/sbs/raw

export OUT_DIR=/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/xtalk_by_events/$runnum
export LOG_DIR=$SWIF_JOB_WORK_DIR

echo 'OUT_DIR='$OUT_DIR
echo 'LOG_DIR='$LOG_DIR

export ANALYZER_CONFIGPATH=$SBS_REPLAY/replay
export SBS=/w/halla-scshelf2102/sbs/SBS_OFFLINE/install
# export SBS=/w/halla-scshelf2102/sbs/jboyd/SBS_OFFLINE/install
cp $SBS/run_replay_here/.rootrc ./

#cmd="analyzer -b -q 'replay/replay_BBGEM.C+("$runnum","$firstsegment","$maxsegments")'"

#echo $cmd

echo "Using this DB directory: " $DB_DIR

echo "Showing contents of SWIF_JOB_WORK_DIR: "
ls -la $SWIF_JOB_WORK_DIR

# movepdfdir=/work/halla/sbs/jboyd/xtalk/xtalk_ratio_outputs/plots/$runnum
movepdfdir=/w/halla-scshelf2102/sbs/jboyd/xtalk/xtalk_ratio_outputs/plots/$runnum/ADCcut$adccut
echo "PDF output dir: " $movepdfdir

if [ ! -d $movepdfdir ]
then
	echo "movepdfdir doesn't exist... creating it. " 
	mkdir -p $movepdfdir
fi

movedatdir=/w/halla-scshelf2102/sbs/jboyd/xtalk/xtalk_ratio_outputs/txt/$runnum/ADCcut$adccut
echo "Dat-file output dir: " $movedatdir

if [ ! -d $movedatdir ]
then
	echo "movedatdir doesn't exist... creating it. " 
	mkdir -p $movedatdir
fi

# analyzer -b -q '/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/xtalk_replay/replay_xtalk.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments')' 

echo "Running the following script: "
# echo "/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/test/xtalk_farm.C"
echo "/work/halla/sbs/jboyd/xtalk/fit_xtalk/xtalk_fits_farm.C+('$runnum','$adccut')"
# echo "/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/test/xtalk_histograms_farm.C+('$runnum', '$lastevent', '$segment')"

# analyzer -b -q /w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/test/xtalk_farm.C+
# analyzer -b -q '/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/test/xtalk_all_farm.C+('$runnum', '$lastevent', '$segment')'
# analyzer -b -q '/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/test/xtalk_histograms_farm.C+('$runnum', '$lastevent', '$segment')'
# analyzer -b -q '/work/halla/sbs/jboyd/xtalk/fit_xtalk/xtalk_fits.C+('$runnum','$adccut')' 
analyzer -b -q '/work/halla/sbs/jboyd/xtalk/fit_xtalk/xtalk_fits_farm.C+('$runnum','$adccut')' 

# moverootdir=/lustre19/expphy/volatile/halla/sbs/jboyd/swif_output/xtalk/histograms/$runnum/
# movelogdir=/lustre19/expphy/volatile/halla/sbs/jboyd/logs

echo "Analysis finished. Looking for created output files...." 

echo "PDF directory: " $movepdfdir
echo "Contents: "
ls $movepdfdir

echo ""
echo "TXT directory: " $movedatdir
echo "Contents: " 
ls $movedatdir

# echo "After the analyzer OUT_DIR is: " $OUT_DIR
#mv $outfilename $newfilename
#mv $logfilename $newlogname

# echo "Moving: " $SWIF_JOB_WORK_DIR/*.pdf " from " $PWD " to " $movepdfdir
# mv -v $SWIF_JOB_WORK_DIR/*.pdf $movepdfdir

# echo "Moving: " $SWIF_JOB_WORK_DIR/*.dat " from " $PWD " to " $movedatdir
# mv -v $SWIF_JOB_WORK_DIR/*.pdf $movedatdir

# echo "Moving: " $SWIF_JOB_WORK_DIR/*.log " from " $PWD " to " $movelogdir
# mv -v $SWIF_JOB_WORK_DIR/*.log $movelogdir


