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

local_file_all='$SWIF_JOB_WORK_DIR/'$runnum'/'$runnum'_xtalk_ratios_ALL_merged.root'
local_file_last='$SWIF_JOB_WORK_DIR/'$runnum'/'$runnum'_xtalk_ratios_LAST_merged.root'

# echo "Local_file: " $local_file
# remote_uri='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_Histogramming_WITH_corr.root'
# remote_uri='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_xtalk_replay_Histogramming_WITH_corr.root'
remote_uri_all='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/xtalk_by_events/'$runnum'/merged/'$runnum'_xtalk_ratios_ALL_merged.root'
remote_uri_last='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/xtalk_by_events/'$runnum'/merged/'$runnum'_xtalk_ratios_LAST_merged.root'
# remote_uri='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_Histogramming_WITH_corr.root'
# echo "Remote file: " $remote_uri

swif2 add-job -workflow jboyd_xtalk_fits -partition production -name $jobname -cores 1 -disk $disk -ram $memory -input $local_file_all $remote_uri_all -input $local_file_last $remote_uri_last $script $runnum $adccut

