#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --mem-per-cpu=1500
source /site/12gev_phys/softenv.sh 2.4

export ANALYZER=/w/halla-scshelf2102/sbs/ANALYZER
echo "Sourcing ANALYZER env setup: /w/halla-scshelf2102/sbs/ANALYZER/install/bin/setup.sh"
source /w/halla-scshelf2102/sbs/ANALYZER/install/bin/setup.sh

echo "Sourcing SBS_OFFLINE en setup: /w/halla-scshelf2102/sbs/SBS_OFFLINE/install/bin/sbsenv.sh"
source /w/halla-scshelf2102/sbs/SBS_OFFLINE/install/bin/sbsenv.sh

echo 'working directory ='
echo $PWD

echo 'swif_job_work_dir='$SWIF_JOB_WORK_DIR

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

# if [ -d /apps/modulefiles ]; then 
# module use /apps/modulefiles 
# fi 
# # source $jboyd/jbenv.sh
# module load gcc/9.2.0 

# ldd /work/halla/sbs/ANALYZER/install/bin/analyzer |& grep not

echo 'working directory = '$PWD

export SBS_REPLAY=/work/halla/sbs/SBS_REPLAY/SBS-replay
export SBS=/work/halla/sbs/SBS_OFFLINE/install
cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR
export DB_DIR=/work/halla/sbs/jboyd/SBS_REPLAY/SBS-replay/DB
# ------
# export JBOYD_REPLAY=/work/halla/sbs/jboyd/SBS_OFFLINE/SBS_REPLAY/SBS-replay
# export DB_DIR=$JBOYD_REPLAY/DB
# ------
export DATA_DIR=/cache/mss/halla/sbs/raw

export OUT_DIR=$SWIF_JOB_WORK_DIR
#export OUT_DIR=/volatile/halla/sbs/jboyd/swif_output
export LOG_DIR=$SWIF_JOB_WORK_DIR
#export LOG_DIR=/volatile/halla/sbs/jboyd/swif_output

echo 'OUT_DIR='$OUT_DIR
echo 'LOG_DIR='$LOG_DIR

export ANALYZER_CONFIGPATH=$SBS_REPLAY/replay

# ##-+-+-+-+-+-+-+-+-+
# #SBATCH --partition=production
# #SBATCH --account=halla
# #SBATCH --mem-per-cpu=1500

# #cd /work/halla/sbs/puckett/GMN_ANALYSIS

# echo 'working directory ='
# echo $PWD

# echo 'swif_job_work_dir='$SWIF_JOB_WORK_DIR

# # source login stuff since swif2 completely nukes any sensible default software environment
# #source /home/puckett/.cshrc
# #source /home/puckett/.login

# #echo 'before sourcing environment, env = '
# #env 

# #module load gcc/9.2.0

# #!/bin/bash 

# # MODULES=/etc/profile.d/modules.sh 

# # if [[ $(type -t module) != function && -r ${MODULES} ]]; then 
# # source ${MODULES} 
# # fi 

# # if [ -d /apps/modulefiles ]; then 
# # module use /apps/modulefiles 
# # fi 

# # module load gcc/9.2.0 

# source /site/12gev_phys/softenv.sh 2.4

# ldd /work/halla/sbs/ANALYZER/install/bin/analyzer |& grep not

# #source /etc/profile.d/modules.sh

# #module load gcc/9.2.0


# # setup environment for ANALYZER and SBS-offline:

# echo 'working directory = '$PWD

# export ANALYZER=/work/halla/sbs/ANALYZER/install
# source $ANALYZER/bin/setup.sh
# # source /work/halla/sbs/jboyd/SBS_OFFLINE/install/bin/sbsenv.sh
# source /w/halla-scshelf2102/sbs/SBS_OFFLINE/install/bin/sbsenv.sh
# ##cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR


# export SBS_REPLAY=/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay
# export SBS=/w/halla-scshelf2102/sbs/SBS_OFFLINE/install

# export DB_DIR=/w/halla-scshelf2102/sbs/jboyd/SBS_REPLAY/SBS-replay/DB
# # export DB_DIR=/work/halla/sbs/jboyd/SBS_OFFLINE/SBS_REPLAY/SBS-replay/DB_xtalk
# # export DB_DIR=/w/halla-scshelf2102/sbs/SBS_REPLAY/SBS-replay/DB
# #export DB_DIR=/work/halla/sbs/jboyd/SBS_OFFLINE/SBS_REPLAY/SBS-replay/DB
# export DATA_DIR=/lustre19/enp/cache/mss/halla/sbs/raw

# export OUT_DIR=$SWIF_JOB_WORK_DIR
# export LOG_DIR=$SWIF_JOB_WORK_DIR

# #export OUT_DIR=/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/GainMatching
# #export OUT_DIR=/work/halla/sbs/jboyd/Rootfiles/GainMatch
# #export LOG_DIR=/work/halla/sbs/jboyd/logs
# #export LOG_DIR=/lustre19/expphy/volatile/halla/sbs/jboyd/logs/
# echo 'OUT_DIR='$OUT_DIR
# echo 'LOG_DIR='$LOG_DIR
# # try this under swif2:
# #export DATA_DIR=$PWD
# #export OUT_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS1_OPTICS/rootfiles
# #export LOG_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS1_OPTICS/logs
# #mkdir -p /volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles
# #mkdir -p /volatile/halla/sbs/puckett/GMN_REPLAYS/logs

# #export OUT_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS4/rootfiles
# #export LOG_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS4/rootfiles
# export ANALYZER_CONFIGPATH=$SBS_REPLAY/replay

# ##-+-+-+-+-+-+-+-+



runnum=$1
maxevents=$2
firstevent=$3

prefix=$4
firstsegment=$5
maxsegments=$6

# cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR

#cmd="analyzer -b -q 'replay/replay_BBGEM.C+("$runnum","$firstsegment","$maxsegments")'"

#echo $cmd

echo "Using this DB directory: " $DB_DIR
echo "Using this SBS_REPLAY director: " $SBS_REPLAY

analyzer -b -q '/w/halla-scshelf2102/sbs/jboyd/SBS_REPLAY/SBS-replay/replay/replay_xtalk_farm_jobs.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments')' 
# analyzer -b -q 'replay_xtalk.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments')' 

#analyzer -b -q 'replay_gmn.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments')'
#analyzer -b -q 'replay_BBGEM.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments')'

# outfilename=$OUT_DIR'/e1209019_*'$runnum'*.root'
# logfilename=$LOG_DIR'/replay_gmn_'$runnum'*.log'

moverootdir=/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/$runnum/
movelogdir=/lustre19/expphy/volatile/halla/sbs/jboyd/logs


# newfilename=$moverootdir'/GainMatch_1209019_*'$runnum'_'$maxevents'events.root'
# newlogname=$movelogdir'/GainMatch_1209019_*'$runnum'_'$maxevents'events.log'

echo "After the analyzer this is the DB dir: " $DB_DIR
#mv $outfilename $newfilename
#mv $logfilename $newlogname

if [ ! -d $moverootdir ]
	then
		echo "moverootdir doesn't exist for $runnum in xtalk_out dir... creating it. " 
		mkdir -p $moverootdir
fi

echo "Moving: " *.root " from " $PWD " to " $moverootdir
mv -vf *.root $moverootdir

echo "Moving: " *.log " from " $PWD " to " $movelogdir
mv -vf *.log $movelogdir




