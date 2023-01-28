#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --mem-per-cpu=1500

echo 'working directory ='
echo $PWD

runnum=$1
lastevent=$2
lastmodule=$3
segment=$4


#echo 'swif_job_work_dir='$SWIF_JOB_WORK_DIR

# source login stuff since swif2 completely nukes any sensible default software environment
#source /home/puckett/.cshrc
#source /home/puckett/.login

#echo 'before sourcing environment, env = '
#env 

#module load gcc/9.2.0

#!/bin/bash 

# MODULES=/etc/profile.d/modules.sh 

# if [[ $(type -t module) != function && -r ${MODULES} ]]; then 
# source ${MODULES} 
# fi 

# if [ -d /apps/modulefiles ]; then y
# module use /apps/modulefiles 
# fi 

# module load gcc/9.2.0 

source /site/12gev_phys/softenv.sh 2.5

ldd /work/halla/sbs/ANALYZER/install/bin/analyzer |& grep not


echo 'working directory = '$PWD

export ANALYZER=/work/halla/sbs/ANALYZER/install
source $ANALYZER/bin/setup.sh
source /w/halla-scshelf2102/sbs/SBS_OFFLINE/install/bin/sbsenv.sh

export SBS_REPLAY=/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay
export DB_DIR=/w/halla-scshelf2102/sbs/jboyd/SBS_REPLAY/SBS-replay/DB
# export DATA_DIR=/lustre19/enp/cache/mss/halla/sbs/raw
export DATA_DIR=/w/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/$runnum

export OUT_DIR=$SWIF_JOB_WORK_DIR
export LOG_DIR=$SWIF_JOB_WORK_DIR

echo 'OUT_DIR='$OUT_DIR
echo 'LOG_DIR='$LOG_DIR

export ANALYZER_CONFIGPATH=$SBS_REPLAY/replay
cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR

#cmd="analyzer -b -q 'replay/replay_BBGEM.C+("$runnum","$firstsegment","$maxsegments")'"

#echo $cmd

echo "Using this DB directory: " $DB_DIR

echo "Showing contents of SWIF_JOB_WORK_DIR: "
ls $SWIF_JOB_WORK_DIR

# analyzer -b -q '/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/xtalk_replay/replay_xtalk.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments')' 

echo "Running the following script: "
# echo "/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/test/xtalk_farm.C"
echo "/w/halla-scshelf2102/sbs/jboyd/xtalk/xtalk_by_events_farm.C+('$runnum', '$lastevent',  '$segment', '$lastmodule')"
# echo "/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/test/xtalk_histograms_farm.C+('$runnum', '$lastevent', '$segment')"

# analyzer -b -q /w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/test/xtalk_farm.C+
# analyzer -b -q '/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/test/xtalk_all_farm.C+('$runnum', '$lastevent', '$segment')'
# analyzer -b -q '/w/halla-scshelf2102/sbs/jboyd/analysis/xtalk/test/xtalk_histograms_farm.C+('$runnum', '$lastevent', '$segment')'
analyzer -b -q '/w/halla-scshelf2102/sbs/jboyd/xtalk/xtalk_by_events_farm.C+('$runnum', '$lastevent',  '$segment', '$lastmodule')'

moverootdir=/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/xtalk_by_events/$runnum/
# moverootdir=/lustre19/expphy/volatile/halla/sbs/jboyd/swif_output/xtalk/histograms/$runnum/
movelogdir=/lustre19/expphy/volatile/halla/sbs/jboyd/logs

echo "After the analyzer OUT_DIR is: " $OUT_DIR
#mv $outfilename $newfilename
#mv $logfilename $newlogname

echo "Moving: " *.log " to " $movelogdir
mv -v *.log $movelogdir

if [ $lastmodule -eq 1 ]
	then

		if [ ! -d $moverootdir ]
		then
			echo "movedatdir doesn't exist... creating it. " 
			mkdir -p $moverootdir
		fi

		echo "Moving: " *_last_module_only.root " to " $moverootdir
		mv -v *_last_module_only.root $moverootdir
	else

		if [ ! -d $moverootdir ]
		then
			echo "movedatdir doesn't exist... creating it. " 
			mkdir -p $moverootdir
		fi		

		echo "Moving: " *_ALL.root " to " $moverootdir
		mv -v *_ALL.root $moverootdir
fi

