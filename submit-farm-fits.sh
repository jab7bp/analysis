#!/bin/bash

# single-run version for swif2:

runnum=$1
adccut=$2

echo "Runnum: " $runnum ", ADCcut: " $adccut 

disk='115GB'
memory='8GB'

script='/w/halla-scshelf2102/sbs/jboyd/SBS_OFFLINE/install/run_replay_here/farm_fits.sh'
jobname='xtalk_fits'

outfilename='match:*.root'
logfilename='match:*.log'
# echo ls /mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770 | wc -l
# DIR="/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/${runnum}"
SWIF_DIR=$SWIF_JOB_WORK_DIR
# FILES="/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/${runnum}/*"
# FILES="/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/${runnum}/e1209019*"
# FILES="/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/'$runnum'/e1209019*"
# FILES="/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/histo_with_corr/'$runnum'/e1209019*"

cd /lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/xtalk_by_events/$runnum

echo '---------------------------------------------------'
echo 'Adding new swif2 job, runnum='$runnum', segment='$j
echo

##I had a hard time figuring out properly load a bunch of root files into a single root/analyzer run so....
##I opted to just merge all the root files into one and call just that one file. 

##Define local files to use:
local_file_all='$SWIF_JOB_WORK_DIR/'$runnum'/'$runnum'_xtalk_ratios_ALL_merged.root'
local_file_last='$SWIF_JOB_WORK_DIR/'$runnum'/'$runnum'_xtalk_ratios_LAST_merged.root'


##Terminal line format for the filenames:
remote_uri_all='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/xtalk_by_events/'$runnum'/merged/'$runnum'_xtalk_ratios_ALL_merged.root'
remote_uri_last='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/xtalk_by_events/'$runnum'/merged/'$runnum'_xtalk_ratios_LAST_merged.root'

##submit script
swif2 add-job -workflow jboyd_xtalk_fits -partition production -name $jobname -cores 1 -disk $disk -ram $memory -input $local_file_all $remote_uri_all -input $local_file_last $remote_uri_last $script $runnum $adccut

