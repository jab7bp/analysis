#!/bin/bash

# single-run version for swif2:
## allsegs is -1 for all segments
runnum=$1
allsegs=$2
lastevent=$3
firstsegment=$4
lastmodule=$5
lastsegment=$6

echo "Last module is: " $lastmodule

disk='40GB'
memory='40GB'
workflow='jboyd_xtalk_analysis'

script='/w/halla-scshelf2102/sbs/jboyd/SBS_OFFLINE/install/run_replay_here/xtalk-analysis.sh'
jobname='xtalk_analysis'

outfilename='match:*.root'
logfilename='match:*.log'
# echo ls /mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770 | wc -l
DIR="/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/${runnum}"
SWIF_DIR=$SWIF_JOB_WORK_DIR
# FILES="/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/${runnum}/*"
FILES="/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/${runnum}/e1209019*"
# FILES="/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/'$runnum'/e1209019*"
# FILES="/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/histo_with_corr/'$runnum'/e1209019*"

if [ "$5" == "0" ]; then
    echo "Last Module set to 0 (Analyze all Mods) "
    lastmodule=0;
elif [ "$5" == "1" ]; then
    echo "Set to only analyze Last lastmodule. " 
    lastmodule=1;
else
    echo "No input for Last Module. Setting all modules by default"
    lastmodule=0;
fi

echo "Last module set to: " $lastmodule

if [ -n "$6" ]
then 
    echo "Analyzing subset of segments from Segment " $4 " to Segment " $2
    i=$(($2 + 1))

    for ((j=$4; j<$i; j++))
    do
        echo '---------------------------------------------------'
        echo 'Adding new swif2 job, runnum='$runnum', segment='$j
        echo
        local_file='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_'$runnum'_stream0_seg'${j}_${j}'_xtalk_replay_Histogramming_WITH_corr.root'
        echo "Local_file: " $local_file
        # remote_uri='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_Histogramming_WITH_corr.root'
        remote_uri='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_xtalk_replay_Histogramming_WITH_corr.root'
        # remote_uri='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/histo_with_corr/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_Histogramming_WITH_corr.root'
        # remote_uri='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_Histogramming_WITH_corr.root'
        echo "Remote file: " $remote_uri

        swif2 add-job -workflow $workflow -partition production -name $jobname -cores 1 -disk $disk -ram $memory -input $local_file $remote_uri $script $runnum $lastevent $lastmodule $j

    done

else
    
    if [ "$2" == "-1" ] 
    then 
        echo "Analyzing all events in " $runnum " directory."
        i=0
        for f in ${FILES}
        do
        #     declare local_file${i}=""
        #     declare remote_uri${i}=""

        #     local_file$i="$SWIF_JOB_WORK_DIR/e1209019_fullreplay_${runnum}_stream0_seg${i}_${i}_Histogramming_WITH_corr.root"
        #     remote_uri$i="mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/${runnum}/e1209019_fullreplay_${runnum}_stream0_seg${i}_${i}_Histogramming_WITH_corr.root"
            i=$((i+1))
        done
        echo "Number of files in File Dir: " $i
        
        echo "Last segment: " "$(($i-1))"

        for ((j=0; j<$i; j++))
        do
            echo 'Adding new swif2 job, runnum='$runnum', segment='$j
            local_file='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_'$runnum'_stream0_seg'${j}_${j}'_xtalk_replay_Histogramming_WITH_corr.root'
            echo "Local_file: " $local_file
            # remote_uri='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_Histogramming_WITH_corr.root'
            remote_uri='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_xtalk_replay_Histogramming_WITH_corr.root'
            # remote_uri='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/histo_with_corr/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_Histogramming_WITH_corr.root'
            # remote_uri='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_Histogramming_WITH_corr.root'
            echo "Remote file: " $remote_uri

            swif2 add-job -workflow $workflow -partition production -name $jobname -cores 1 -disk $disk -ram $memory -input $local_file $remote_uri $script $runnum $lastevent $lastmodule $j

        done

    else
            num_segs=$2

            echo "Number of segments to analyze: " $2

            for ((j=0; j<$num_segs; j++))
            do
                echo 'Adding new swif2 job, runnum='$runnum', segment='$j
                local_file='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_'$runnum'_stream0_seg'${j}_${j}'_xtalk_replay_Histogramming_WITH_corr.root'
                echo "Local_file: " $local_file
                # remote_uri='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_Histogramming_WITH_corr.root'
                remote_uri='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_xtalk_replay_Histogramming_WITH_corr.root'
                # remote_uri='file:/lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/histo_with_corr/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}'_'${j}'_Histogramming_WITH_corr.root'
                # remote_uri='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/xtalk/'$runnum'/e1209019_fullreplay_'$runnum'_stream0_seg'${j}_${j}'_Histogramming_WITH_corr.root'
                echo "Remote file: " $remote_file

                swif2 add-job -workflow $workflow -partition production -name $jobname -cores 1 -disk $disk -ram $memory -input $local_file $remote_uri $script $runnum $lastevent $lastmodule $j

            done
    fi
fi

echo "Workflow: " $workflow 


# # echo "local_file3: " ${local_file3}
# # echo "remote_uri4: " ${remote_uri4}

# # swif2 add-job -workflow jboyd_xtalk_analysis -partition production -name $jobname -cores 1 -disk 500GB -ram 500GB -input $local_file0 $remote_uri0

# # # j=0
# # # swif2 add-job -workflow jboyd_xtalk_analysis -partition production -name $jobname -cores 1 -disk 300GB -ram 10000MB -input $local_file0 $remote_uri0 -input $local_file1 $remote_uri1 $script

# # # local_file='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg*'
# # # local_file0='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_$runnum_stream0_seg0_0_Histogramming_WITH_corr.root'
# # # local_file1='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_$runnum_stream0_seg1_1_Histogramming_WITH_corr.root'
# # # local_file2='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_$runnum_stream0_seg2_2_Histogramming_WITH_corr.root'
# # # local_file3='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_$runnum_stream0_seg3_3_Histogramming_WITH_corr.root'
# # # local_file4='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_$runnum_stream0_seg4_4_Histogramming_WITH_corr.root'
# # # # remote_uri='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/e1209019_fullreplay_13770_stream0_seg*root'
# # # remote_uri0='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/$runnum/e1209019_fullreplay_$runnum_stream0_seg0_0_Histogramming_WITH_corr.root'
# # # remote_uri1='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/$runnum/e1209019_fullreplay_$runnum_stream0_seg1_1_Histogramming_WITH_corr.root'
# # # remote_uri2='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/$runnum/e1209019_fullreplay_$runnum_stream0_seg2_2_Histogramming_WITH_corr.root'
# # # remote_uri3='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/$runnum/e1209019_fullreplay_$runnum_stream0_seg3_3_Histogramming_WITH_corr.root'
# # # remote_uri4='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/$runnum/e1209019_fullreplay_$runnum_stream0_seg4_4_Histogramming_WITH_corr.root'

# # local_file0='$SWIF_JOB_WORK_DIR/13770_Histogramming_WITH_corr_MERGED.root'
# # local_file0='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg0_0_Histogramming_WITH_corr.root'
# # local_file1='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg1_1_Histogramming_WITH_corr.root'
# # local_file2='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg2_2_Histogramming_WITH_corr.root'
# # local_file3='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg3_3_Histogramming_WITH_corr.root'
# # local_file4='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg4_4_Histogramming_WITH_corr.root'
# # local_file5='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg5_5_Histogramming_WITH_corr.root'
# # local_file6='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg6_6_Histogramming_WITH_corr.root'
# # local_file7='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg7_7_Histogramming_WITH_corr.root'
# # local_file8='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg8_8_Histogramming_WITH_corr.root'
# # local_file9='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg9_9_Histogramming_WITH_corr.root'
# # local_file10='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg10_10_Histogramming_WITH_corr.root'
# # local_file11='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg11_11_Histogramming_WITH_corr.root'
# # local_file12='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg12_12_Histogramming_WITH_corr.root'
# # local_file13='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg13_13_Histogramming_WITH_corr.root'
# # local_file14='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg14_14_Histogramming_WITH_corr.root'
# # local_file15='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg15_15_Histogramming_WITH_corr.root'
# # local_file16='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg16_16_Histogramming_WITH_corr.root'
# # local_file17='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg17_17_Histogramming_WITH_corr.root'
# # local_file18='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg18_18_Histogramming_WITH_corr.root'
# # local_file19='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg19_19_Histogramming_WITH_corr.root'
# # local_file20='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg20_20_Histogramming_WITH_corr.root'
# # local_file21='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg21_21_Histogramming_WITH_corr.root'
# # local_file22='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg22_22_Histogramming_WITH_corr.root'
# # local_file23='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg23_23_Histogramming_WITH_corr.root'
# # local_file24='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg24_24_Histogramming_WITH_corr.root'
# # local_file25='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg25_25_Histogramming_WITH_corr.root'
# # local_file26='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg26_26_Histogramming_WITH_corr.root'
# # local_file27='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg27_27_Histogramming_WITH_corr.root'
# # local_file28='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg28_28_Histogramming_WITH_corr.root'
# # local_file29='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg29_29_Histogramming_WITH_corr.root'
# # local_file30='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg30_30_Histogramming_WITH_corr.root'
# # local_file31='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg31_31_Histogramming_WITH_corr.root'
# # local_file32='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg32_32_Histogramming_WITH_corr.root'
# # local_file33='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg33_33_Histogramming_WITH_corr.root'
# # local_file34='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg34_34_Histogramming_WITH_corr.root'
# # local_file35='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg35_35_Histogramming_WITH_corr.root'
# # local_file36='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg36_36_Histogramming_WITH_corr.root'
# # local_file37='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg37_37_Histogramming_WITH_corr.root'
# # local_file38='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg38_38_Histogramming_WITH_corr.root'
# # local_file39='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg39_39_Histogramming_WITH_corr.root'
# # local_file40='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg40_40_Histogramming_WITH_corr.root'
# # local_file41='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg41_41_Histogramming_WITH_corr.root'
# # local_file42='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg42_42_Histogramming_WITH_corr.root'
# # local_file43='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg43_43_Histogramming_WITH_corr.root'
# # local_file44='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg44_44_Histogramming_WITH_corr.root'
# # local_file45='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg45_45_Histogramming_WITH_corr.root'
# # local_file46='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg46_46_Histogramming_WITH_corr.root'
# # local_file47='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg47_47_Histogramming_WITH_corr.root'
# # local_file48='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg48_48_Histogramming_WITH_corr.root'
# # local_file49='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg49_49_Histogramming_WITH_corr.root'
# # local_file50='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg50_50_Histogramming_WITH_corr.root'
# # local_file51='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg51_51_Histogramming_WITH_corr.root'
# # local_file52='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg52_52_Histogramming_WITH_corr.root'
# # local_file53='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg53_53_Histogramming_WITH_corr.root'
# # local_file54='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg54_54_Histogramming_WITH_corr.root'
# # local_file55='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg55_55_Histogramming_WITH_corr.root'
# # local_file56='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg56_56_Histogramming_WITH_corr.root'
# # local_file57='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg57_57_Histogramming_WITH_corr.root'
# # local_file58='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg58_58_Histogramming_WITH_corr.root'
# # local_file59='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg59_59_Histogramming_WITH_corr.root'
# # local_file60='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg60_60_Histogramming_WITH_corr.root'
# # local_file61='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg61_61_Histogramming_WITH_corr.root'
# # local_file62='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg62_62_Histogramming_WITH_corr.root'
# # local_file63='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg63_63_Histogramming_WITH_corr.root'
# # local_file64='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg64_64_Histogramming_WITH_corr.root'
# # local_file65='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg65_65_Histogramming_WITH_corr.root'
# # local_file66='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg66_66_Histogramming_WITH_corr.root'
# # local_file67='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg67_67_Histogramming_WITH_corr.root'
# # local_file68='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg68_68_Histogramming_WITH_corr.root'
# # local_file69='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg69_69_Histogramming_WITH_corr.root'
# # local_file70='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg70_70_Histogramming_WITH_corr.root'
# # local_file71='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg71_71_Histogramming_WITH_corr.root'
# # local_file72='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg72_72_Histogramming_WITH_corr.root'
# # local_file73='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg73_73_Histogramming_WITH_corr.root'
# # local_file74='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg74_74_Histogramming_WITH_corr.root'
# # local_file75='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg75_75_Histogramming_WITH_corr.root'
# # local_file76='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg76_76_Histogramming_WITH_corr.root'
# # local_file77='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg77_77_Histogramming_WITH_corr.root'
# # local_file78='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg78_78_Histogramming_WITH_corr.root'
# # local_file79='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg79_79_Histogramming_WITH_corr.root'
# # local_file80='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg80_80_Histogramming_WITH_corr.root'
# # local_file81='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg81_81_Histogramming_WITH_corr.root'
# # local_file82='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg82_82_Histogramming_WITH_corr.root'
# # local_file83='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg83_83_Histogramming_WITH_corr.root'
# # local_file84='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg84_84_Histogramming_WITH_corr.root'
# # local_file85='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg85_85_Histogramming_WITH_corr.root'
# # local_file86='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg86_86_Histogramming_WITH_corr.root'
# # local_file87='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg87_87_Histogramming_WITH_corr.root'
# # local_file88='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg88_88_Histogramming_WITH_corr.root'
# # local_file89='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg89_89_Histogramming_WITH_corr.root'
# # local_file90='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg90_90_Histogramming_WITH_corr.root'
# # local_file91='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg91_91_Histogramming_WITH_corr.root'
# # local_file92='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg92_92_Histogramming_WITH_corr.root'
# # local_file93='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg93_93_Histogramming_WITH_corr.root'
# # local_file94='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg94_94_Histogramming_WITH_corr.root'
# # local_file95='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg95_95_Histogramming_WITH_corr.root'
# # local_file96='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg96_96_Histogramming_WITH_corr.root'
# # local_file97='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg97_97_Histogramming_WITH_corr.root'
# # local_file98='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg98_98_Histogramming_WITH_corr.root'
# # local_file99='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg99_99_Histogramming_WITH_corr.root'
# # local_file100='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg100_100_Histogramming_WITH_corr.root'
# # local_file101='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg101_101_Histogramming_WITH_corr.root'
# # local_file102='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg102_102_Histogramming_WITH_corr.root'
# # local_file103='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg103_103_Histogramming_WITH_corr.root'
# # local_file104='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg104_104_Histogramming_WITH_corr.root'
# # local_file105='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg105_105_Histogramming_WITH_corr.root'
# # local_file106='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg106_106_Histogramming_WITH_corr.root'
# # local_file107='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg107_107_Histogramming_WITH_corr.root'
# # local_file108='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg108_108_Histogramming_WITH_corr.root'
# # local_file109='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg109_109_Histogramming_WITH_corr.root'
# # local_file110='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg110_110_Histogramming_WITH_corr.root'
# # local_file111='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg111_111_Histogramming_WITH_corr.root'
# # local_file112='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg112_112_Histogramming_WITH_corr.root'
# # local_file113='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg113_113_Histogramming_WITH_corr.root'
# # local_file114='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg114_114_Histogramming_WITH_corr.root'
# # local_file115='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg115_115_Histogramming_WITH_corr.root'
# # local_file116='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg116_116_Histogramming_WITH_corr.root'
# # local_file117='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg117_117_Histogramming_WITH_corr.root'
# # local_file118='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg118_118_Histogramming_WITH_corr.root'
# # local_file119='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg119_119_Histogramming_WITH_corr.root'
# # local_file120='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg120_120_Histogramming_WITH_corr.root'
# # local_file121='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg121_121_Histogramming_WITH_corr.root'
# # local_file122='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg122_122_Histogramming_WITH_corr.root'
# # local_file123='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg123_123_Histogramming_WITH_corr.root'
# # local_file124='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg124_124_Histogramming_WITH_corr.root'
# # local_file125='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg125_125_Histogramming_WITH_corr.root'
# # local_file126='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg126_126_Histogramming_WITH_corr.root'
# # local_file127='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg127_127_Histogramming_WITH_corr.root'
# # local_file128='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg128_128_Histogramming_WITH_corr.root'
# # local_file129='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg129_129_Histogramming_WITH_corr.root'
# # local_file130='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg130_130_Histogramming_WITH_corr.root'
# # local_file131='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg131_131_Histogramming_WITH_corr.root'
# # local_file132='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg132_132_Histogramming_WITH_corr.root'
# # local_file133='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg133_133_Histogramming_WITH_corr.root'
# # local_file134='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg134_134_Histogramming_WITH_corr.root'
# # local_file135='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg135_135_Histogramming_WITH_corr.root'
# # local_file136='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg136_136_Histogramming_WITH_corr.root'
# # local_file137='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg137_137_Histogramming_WITH_corr.root'
# # local_file138='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg138_138_Histogramming_WITH_corr.root'
# # local_file139='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg139_139_Histogramming_WITH_corr.root'
# # local_file140='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg140_140_Histogramming_WITH_corr.root'
# # local_file141='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg141_141_Histogramming_WITH_corr.root'
# # local_file142='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg142_142_Histogramming_WITH_corr.root'
# # local_file143='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg143_143_Histogramming_WITH_corr.root'
# # local_file144='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg144_144_Histogramming_WITH_corr.root'
# # local_file145='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg145_145_Histogramming_WITH_corr.root'
# # local_file146='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg146_146_Histogramming_WITH_corr.root'
# # local_file147='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg147_147_Histogramming_WITH_corr.root'
# # local_file148='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg148_148_Histogramming_WITH_corr.root'
# # local_file149='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg149_149_Histogramming_WITH_corr.root'
# # local_file150='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg150_150_Histogramming_WITH_corr.root'
# # local_file151='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg151_151_Histogramming_WITH_corr.root'
# # local_file152='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg152_152_Histogramming_WITH_corr.root'
# # local_file153='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg153_153_Histogramming_WITH_corr.root'
# # local_file154='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg154_154_Histogramming_WITH_corr.root'
# # local_file155='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg155_155_Histogramming_WITH_corr.root'
# # local_file156='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg156_156_Histogramming_WITH_corr.root'
# # local_file157='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg157_157_Histogramming_WITH_corr.root'
# # local_file158='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg158_158_Histogramming_WITH_corr.root'
# # local_file159='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg159_159_Histogramming_WITH_corr.root'
# # local_file160='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg160_160_Histogramming_WITH_corr.root'
# # local_file161='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg161_161_Histogramming_WITH_corr.root'
# # local_file162='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg162_162_Histogramming_WITH_corr.root'
# # local_file163='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg163_163_Histogramming_WITH_corr.root'
# # local_file164='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg164_164_Histogramming_WITH_corr.root'
# # local_file165='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg165_165_Histogramming_WITH_corr.root'
# # local_file166='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg166_166_Histogramming_WITH_corr.root'
# # local_file167='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg167_167_Histogramming_WITH_corr.root'
# # local_file168='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg168_168_Histogramming_WITH_corr.root'
# # local_file169='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg169_169_Histogramming_WITH_corr.root'
# # local_file170='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg170_170_Histogramming_WITH_corr.root'
# # local_file171='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg171_171_Histogramming_WITH_corr.root'
# # local_file172='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg172_172_Histogramming_WITH_corr.root'
# # local_file173='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg173_173_Histogramming_WITH_corr.root'
# # local_file174='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg174_174_Histogramming_WITH_corr.root'
# # local_file175='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg175_175_Histogramming_WITH_corr.root'
# # local_file176='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg176_176_Histogramming_WITH_corr.root'
# # local_file177='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg177_177_Histogramming_WITH_corr.root'
# # local_file178='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg178_178_Histogramming_WITH_corr.root'
# # local_file179='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg179_179_Histogramming_WITH_corr.root'
# # local_file180='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg180_180_Histogramming_WITH_corr.root'
# # local_file181='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg181_181_Histogramming_WITH_corr.root'
# # local_file182='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg182_182_Histogramming_WITH_corr.root'
# # local_file183='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg183_183_Histogramming_WITH_corr.root'
# # local_file184='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg184_184_Histogramming_WITH_corr.root'
# # local_file185='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg185_185_Histogramming_WITH_corr.root'
# # local_file186='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg186_186_Histogramming_WITH_corr.root'
# # local_file187='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg187_187_Histogramming_WITH_corr.root'
# # local_file188='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg188_188_Histogramming_WITH_corr.root'
# # local_file189='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg189_189_Histogramming_WITH_corr.root'
# # local_file190='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg190_190_Histogramming_WITH_corr.root'
# # local_file191='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg191_191_Histogramming_WITH_corr.root'
# # local_file192='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg192_192_Histogramming_WITH_corr.root'
# # local_file193='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg193_193_Histogramming_WITH_corr.root'
# # local_file194='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg194_194_Histogramming_WITH_corr.root'
# # local_file195='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg195_195_Histogramming_WITH_corr.root'
# # local_file196='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg196_196_Histogramming_WITH_corr.root'
# # local_file197='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg197_197_Histogramming_WITH_corr.root'
# # local_file198='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg198_198_Histogramming_WITH_corr.root'
# # local_file199='$SWIF_JOB_WORK_DIR/e1209019_fullreplay_13770_stream0_seg199_199_Histogramming_WITH_corr.root'

# # remote_uri0='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/merged/13770_Histogramming_WITH_corr_MERGED.root'
# # remote_uri0='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg0_0_Histogramming_WITH_corr.root'
# # remote_uri1='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg1_1_Histogramming_WITH_corr.root'
# # remote_uri2='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg2_2_Histogramming_WITH_corr.root'
# # remote_uri3='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg3_3_Histogramming_WITH_corr.root'
# # remote_uri4='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg4_4_Histogramming_WITH_corr.root'
# # remote_uri5='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg5_5_Histogramming_WITH_corr.root'
# # remote_uri6='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg6_6_Histogramming_WITH_corr.root'
# # remote_uri7='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg7_7_Histogramming_WITH_corr.root'
# # remote_uri8='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg8_8_Histogramming_WITH_corr.root'
# # remote_uri9='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg9_9_Histogramming_WITH_corr.root'
# # remote_uri10='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg10_10_Histogramming_WITH_corr.root'
# # remote_uri11='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg11_11_Histogramming_WITH_corr.root'
# # remote_uri12='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg12_12_Histogramming_WITH_corr.root'
# # remote_uri13='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg13_13_Histogramming_WITH_corr.root'
# # remote_uri14='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg14_14_Histogramming_WITH_corr.root'
# # remote_uri15='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg15_15_Histogramming_WITH_corr.root'
# # remote_uri16='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg16_16_Histogramming_WITH_corr.root'
# # remote_uri17='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg17_17_Histogramming_WITH_corr.root'
# # remote_uri18='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg18_18_Histogramming_WITH_corr.root'
# # remote_uri19='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg19_19_Histogramming_WITH_corr.root'
# # remote_uri20='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg20_20_Histogramming_WITH_corr.root'
# # remote_uri21='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg21_21_Histogramming_WITH_corr.root'
# # remote_uri22='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg22_22_Histogramming_WITH_corr.root'
# # remote_uri23='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg23_23_Histogramming_WITH_corr.root'
# # remote_uri24='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg24_24_Histogramming_WITH_corr.root'
# # remote_uri25='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg25_25_Histogramming_WITH_corr.root'
# # remote_uri26='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg26_26_Histogramming_WITH_corr.root'
# # remote_uri27='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg27_27_Histogramming_WITH_corr.root'
# # remote_uri28='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg28_28_Histogramming_WITH_corr.root'
# # remote_uri29='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg29_29_Histogramming_WITH_corr.root'
# # remote_uri30='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg30_30_Histogramming_WITH_corr.root'
# # remote_uri31='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg31_31_Histogramming_WITH_corr.root'
# # remote_uri32='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg32_32_Histogramming_WITH_corr.root'
# # remote_uri33='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg33_33_Histogramming_WITH_corr.root'
# # remote_uri34='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg34_34_Histogramming_WITH_corr.root'
# # remote_uri35='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg35_35_Histogramming_WITH_corr.root'
# # remote_uri36='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg36_36_Histogramming_WITH_corr.root'
# # remote_uri37='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg37_37_Histogramming_WITH_corr.root'
# # remote_uri38='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg38_38_Histogramming_WITH_corr.root'
# # remote_uri39='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg39_39_Histogramming_WITH_corr.root'
# # remote_uri40='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg40_40_Histogramming_WITH_corr.root'
# # remote_uri41='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg41_41_Histogramming_WITH_corr.root'
# # remote_uri42='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg42_42_Histogramming_WITH_corr.root'
# # remote_uri43='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg43_43_Histogramming_WITH_corr.root'
# # remote_uri44='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg44_44_Histogramming_WITH_corr.root'
# # remote_uri45='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg45_45_Histogramming_WITH_corr.root'
# # remote_uri46='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg46_46_Histogramming_WITH_corr.root'
# # remote_uri47='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg47_47_Histogramming_WITH_corr.root'
# # remote_uri48='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg48_48_Histogramming_WITH_corr.root'
# # remote_uri49='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg49_49_Histogramming_WITH_corr.root'
# # remote_uri50='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg50_50_Histogramming_WITH_corr.root'
# # remote_uri51='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg51_51_Histogramming_WITH_corr.root'
# # remote_uri52='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg52_52_Histogramming_WITH_corr.root'
# # remote_uri53='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg53_53_Histogramming_WITH_corr.root'
# # remote_uri54='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg54_54_Histogramming_WITH_corr.root'
# # remote_uri55='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg55_55_Histogramming_WITH_corr.root'
# # remote_uri56='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg56_56_Histogramming_WITH_corr.root'
# # remote_uri57='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg57_57_Histogramming_WITH_corr.root'
# # remote_uri58='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg58_58_Histogramming_WITH_corr.root'
# # remote_uri59='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg59_59_Histogramming_WITH_corr.root'
# # remote_uri60='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg60_60_Histogramming_WITH_corr.root'
# # remote_uri61='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg61_61_Histogramming_WITH_corr.root'
# # remote_uri62='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg62_62_Histogramming_WITH_corr.root'
# # remote_uri63='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg63_63_Histogramming_WITH_corr.root'
# # remote_uri64='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg64_64_Histogramming_WITH_corr.root'
# # remote_uri65='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg65_65_Histogramming_WITH_corr.root'
# # remote_uri66='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg66_66_Histogramming_WITH_corr.root'
# # remote_uri67='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg67_67_Histogramming_WITH_corr.root'
# # remote_uri68='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg68_68_Histogramming_WITH_corr.root'
# # remote_uri69='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg69_69_Histogramming_WITH_corr.root'
# # remote_uri70='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg70_60_Histogramming_WITH_corr.root'
# # remote_uri71='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg71_71_Histogramming_WITH_corr.root'
# # remote_uri72='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg72_72_Histogramming_WITH_corr.root'
# # remote_uri73='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg73_73_Histogramming_WITH_corr.root'
# # remote_uri74='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg74_74_Histogramming_WITH_corr.root'
# # remote_uri75='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg75_75_Histogramming_WITH_corr.root'
# # remote_uri76='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg76_76_Histogramming_WITH_corr.root'
# # remote_uri77='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg77_77_Histogramming_WITH_corr.root'
# # remote_uri78='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg78_78_Histogramming_WITH_corr.root'
# # remote_uri79='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg79_79_Histogramming_WITH_corr.root'
# # remote_uri80='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg80_80_Histogramming_WITH_corr.root'
# # remote_uri81='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg81_81_Histogramming_WITH_corr.root'
# # remote_uri82='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg82_82_Histogramming_WITH_corr.root'
# # remote_uri83='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg83_83_Histogramming_WITH_corr.root'
# # remote_uri84='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg84_84_Histogramming_WITH_corr.root'
# # remote_uri85='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg85_85_Histogramming_WITH_corr.root'
# # remote_uri86='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg86_86_Histogramming_WITH_corr.root'
# # remote_uri87='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg87_87_Histogramming_WITH_corr.root'
# # remote_uri88='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg88_88_Histogramming_WITH_corr.root'
# # remote_uri89='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg89_89_Histogramming_WITH_corr.root'
# # remote_uri90='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg90_90_Histogramming_WITH_corr.root'
# # remote_uri91='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg91_91_Histogramming_WITH_corr.root'
# # remote_uri92='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg92_92_Histogramming_WITH_corr.root'
# # remote_uri93='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg93_93_Histogramming_WITH_corr.root'
# # remote_uri94='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg94_94_Histogramming_WITH_corr.root'
# # remote_uri95='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg95_95_Histogramming_WITH_corr.root'
# # remote_uri96='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg96_96_Histogramming_WITH_corr.root'
# # remote_uri97='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg97_97_Histogramming_WITH_corr.root'
# # remote_uri98='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg98_98_Histogramming_WITH_corr.root'
# # remote_uri99='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg99_99_Histogramming_WITH_corr.root'
# # remote_uri100='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg100_100_Histogramming_WITH_corr.root'
# # remote_uri101='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg101_101_Histogramming_WITH_corr.root'
# # remote_uri102='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg102_102_Histogramming_WITH_corr.root'
# # remote_uri103='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg103_103_Histogramming_WITH_corr.root'
# # remote_uri104='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg104_104_Histogramming_WITH_corr.root'
# # remote_uri105='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg105_105_Histogramming_WITH_corr.root'
# # remote_uri106='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg106_106_Histogramming_WITH_corr.root'
# # remote_uri107='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg107_107_Histogramming_WITH_corr.root'
# # remote_uri108='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg108_108_Histogramming_WITH_corr.root'
# # remote_uri109='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg109_109_Histogramming_WITH_corr.root'
# # remote_uri110='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg110_110_Histogramming_WITH_corr.root'
# # remote_uri111='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg111_111_Histogramming_WITH_corr.root'
# # remote_uri112='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg112_112_Histogramming_WITH_corr.root'
# # remote_uri113='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg113_113_Histogramming_WITH_corr.root'
# # remote_uri114='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg114_114_Histogramming_WITH_corr.root'
# # remote_uri115='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg115_115_Histogramming_WITH_corr.root'
# # remote_uri116='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg116_116_Histogramming_WITH_corr.root'
# # remote_uri117='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg117_117_Histogramming_WITH_corr.root'
# # remote_uri118='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg118_118_Histogramming_WITH_corr.root'
# # remote_uri119='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg119_119_Histogramming_WITH_corr.root'
# # remote_uri120='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg120_120_Histogramming_WITH_corr.root'
# # remote_uri121='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg121_121_Histogramming_WITH_corr.root'
# # remote_uri122='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg122_122_Histogramming_WITH_corr.root'
# # remote_uri123='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg123_123_Histogramming_WITH_corr.root'
# # remote_uri124='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg124_124_Histogramming_WITH_corr.root'
# # remote_uri125='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg125_125_Histogramming_WITH_corr.root'
# # remote_uri126='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg126_126_Histogramming_WITH_corr.root'
# # remote_uri127='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg127_127_Histogramming_WITH_corr.root'
# # remote_uri128='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg128_128_Histogramming_WITH_corr.root'
# # remote_uri129='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg129_129_Histogramming_WITH_corr.root'
# # remote_uri130='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg130_130_Histogramming_WITH_corr.root'
# # remote_uri131='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg131_131_Histogramming_WITH_corr.root'
# # remote_uri132='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg132_132_Histogramming_WITH_corr.root'
# # remote_uri133='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg133_133_Histogramming_WITH_corr.root'
# # remote_uri134='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg134_134_Histogramming_WITH_corr.root'
# # remote_uri135='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg135_135_Histogramming_WITH_corr.root'
# # remote_uri136='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg136_136_Histogramming_WITH_corr.root'
# # remote_uri137='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg137_137_Histogramming_WITH_corr.root'
# # remote_uri138='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg138_138_Histogramming_WITH_corr.root'
# # remote_uri139='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg139_139_Histogramming_WITH_corr.root'
# # remote_uri140='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg140_140_Histogramming_WITH_corr.root'
# # remote_uri141='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg141_141_Histogramming_WITH_corr.root'
# # remote_uri142='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg142_142_Histogramming_WITH_corr.root'
# # remote_uri143='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg143_143_Histogramming_WITH_corr.root'
# # remote_uri144='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg144_144_Histogramming_WITH_corr.root'
# # remote_uri145='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg145_145_Histogramming_WITH_corr.root'
# # remote_uri146='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg146_146_Histogramming_WITH_corr.root'
# # remote_uri147='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg147_147_Histogramming_WITH_corr.root'
# # remote_uri148='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg148_148_Histogramming_WITH_corr.root'
# # remote_uri149='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg149_149_Histogramming_WITH_corr.root'
# # remote_uri150='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg150_150_Histogramming_WITH_corr.root'
# # remote_uri151='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg151_151_Histogramming_WITH_corr.root'
# # remote_uri152='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg152_152_Histogramming_WITH_corr.root'
# # remote_uri153='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg153_153_Histogramming_WITH_corr.root'
# # remote_uri154='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg154_154_Histogramming_WITH_corr.root'
# # remote_uri155='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg155_155_Histogramming_WITH_corr.root'
# # remote_uri156='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg156_156_Histogramming_WITH_corr.root'
# # remote_uri157='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg157_157_Histogramming_WITH_corr.root'
# # remote_uri158='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg158_158_Histogramming_WITH_corr.root'
# # remote_uri159='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg159_159_Histogramming_WITH_corr.root'
# # remote_uri160='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg160_160_Histogramming_WITH_corr.root'
# # remote_uri161='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg161_161_Histogramming_WITH_corr.root'
# # remote_uri162='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg162_162_Histogramming_WITH_corr.root'
# # remote_uri163='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg163_163_Histogramming_WITH_corr.root'
# # remote_uri164='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg164_164_Histogramming_WITH_corr.root'
# # remote_uri165='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg165_165_Histogramming_WITH_corr.root'
# # remote_uri166='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg166_166_Histogramming_WITH_corr.root'
# # remote_uri167='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg167_167_Histogramming_WITH_corr.root'
# # remote_uri168='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg168_168_Histogramming_WITH_corr.root'
# # remote_uri169='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg169_169_Histogramming_WITH_corr.root'
# # remote_uri170='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg170_160_Histogramming_WITH_corr.root'
# # remote_uri171='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg171_171_Histogramming_WITH_corr.root'
# # remote_uri172='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg172_172_Histogramming_WITH_corr.root'
# # remote_uri173='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg173_173_Histogramming_WITH_corr.root'
# # remote_uri174='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg174_174_Histogramming_WITH_corr.root'
# # remote_uri175='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg175_175_Histogramming_WITH_corr.root'
# # remote_uri176='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg176_176_Histogramming_WITH_corr.root'
# # remote_uri177='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg177_177_Histogramming_WITH_corr.root'
# # remote_uri178='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg178_178_Histogramming_WITH_corr.root'
# # remote_uri179='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg179_179_Histogramming_WITH_corr.root'
# # remote_uri180='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg180_180_Histogramming_WITH_corr.root'
# # remote_uri181='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg181_181_Histogramming_WITH_corr.root'
# # remote_uri182='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg182_182_Histogramming_WITH_corr.root'
# # remote_uri183='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg183_183_Histogramming_WITH_corr.root'
# # remote_uri184='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg184_184_Histogramming_WITH_corr.root'
# # remote_uri185='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg185_185_Histogramming_WITH_corr.root'
# # remote_uri186='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg186_186_Histogramming_WITH_corr.root'
# # remote_uri187='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg187_187_Histogramming_WITH_corr.root'
# # remote_uri188='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg188_188_Histogramming_WITH_corr.root'
# # remote_uri189='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg189_189_Histogramming_WITH_corr.root'
# # remote_uri190='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg190_190_Histogramming_WITH_corr.root'
# # remote_uri191='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg191_191_Histogramming_WITH_corr.root'
# # remote_uri192='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg192_192_Histogramming_WITH_corr.root'
# # remote_uri193='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg193_193_Histogramming_WITH_corr.root'
# # remote_uri194='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg194_194_Histogramming_WITH_corr.root'
# # remote_uri195='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg195_195_Histogramming_WITH_corr.root'
# # remote_uri196='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg196_196_Histogramming_WITH_corr.root'
# # remote_uri197='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg197_197_Histogramming_WITH_corr.root'
# # remote_uri198='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg198_198_Histogramming_WITH_corr.root'
# # remote_uri199='mss:/mss/halla/sbs/prod/GEM_test/UVA_GEM/jboyd/13770/e1209019_fullreplay_13770_stream0_seg199_199_Histogramming_WITH_corr.root'

# # swif2 add-job -workflow jboyd_xtalk_analysis -partition production -name $jobname -cores 1 -disk 500GB -ram 500GB -input $local_file0 $remote_uri0 -input $local_file1 $remote_uri1 -input $local_file2 $remote_uri2 -input $local_file3 $remote_uri3 -input $local_file4 $remote_uri4 -input $local_file5 $remote_uri5 -input $local_file6 $remote_uri6 -input $local_file7 $remote_uri7 -input $local_file8 $remote_uri8 -input $local_file9 $remote_uri9 -input $local_file10 $remote_uri0 -input $local_file11 $remote_uri1 -input $local_file12 $remote_uri2 -input $local_file13 $remote_uri3 -input $local_file14 $remote_uri4 -input $local_file15 $remote_uri5 -input $local_file16 $remote_uri6 -input $local_file17 $remote_uri7 -input $local_file18 $remote_uri8 -input $local_file19 $remote_uri9 -input $local_file20 $remote_uri20 -input $local_file21 $remote_uri21 -input $local_file22 $remote_uri22 -input $local_file23 $remote_uri23 -input $local_file24 $remote_uri24 -input $local_file25 $remote_uri25 -input $local_file26 $remote_uri26 -input $local_file27 $remote_uri27 -input $local_file28 $remote_uri28 -input $local_file29 $remote_uri29 -input $local_file30 $remote_uri30 -input $local_file31 $remote_uri31 -input $local_file32 $remote_uri32 -input $local_file33 $remote_uri33 -input $local_file34 $remote_uri34 -input $local_file35 $remote_uri35 -input $local_file36 $remote_uri36 -input $local_file37 $remote_uri37 -input $local_file38 $remote_uri38 -input $local_file39 $remote_uri39 -input $local_file40 $remote_uri40 -input $local_file41 $remote_uri41 -input $local_file42 $remote_uri42 -input $local_file43 $remote_uri43 -input $local_file44 $remote_uri44 -input $local_file45 $remote_uri45 -input $local_file46 $remote_uri46 -input $local_file47 $remote_uri47 -input $local_file48 $remote_uri48 -input $local_file49 $remote_uri49 -input $local_file50 $remote_uri50 -input $local_file51 $remote_uri51 -input $local_file52 $remote_uri52 -input $local_file53 $remote_uri53 -input $local_file54 $remote_uri54 -input $local_file55 $remote_uri55 -input $local_file56 $remote_uri56 -input $local_file57 $remote_uri57 -input $local_file58 $remote_uri58 -input $local_file59 $remote_uri59 -input $local_file60 $remote_uri60 -input $local_file61 $remote_uri61 -input $local_file62 $remote_uri62 -input $local_file63 $remote_uri63 -input $local_file64 $remote_uri64 -input $local_file65 $remote_uri65 -input $local_file66 $remote_uri66 -input $local_file67 $remote_uri67 -input $local_file68 $remote_uri68 -input $local_file69 $remote_uri69 -input $local_file70 $remote_uri70 -input $local_file71 $remote_uri71 -input $local_file72 $remote_uri72 -input $local_file73 $remote_uri73 -input $local_file74 $remote_uri74 -input $local_file75 $remote_uri75 -input $local_file76 $remote_uri76 -input $local_file77 $remote_uri77 -input $local_file78 $remote_uri78 -input $local_file79 $remote_uri79 -input $local_file80 $remote_uri80 -input $local_file81 $remote_uri81 -input $local_file82 $remote_uri82 -input $local_file83 $remote_uri83 -input $local_file84 $remote_uri84 -input $local_file85 $remote_uri85 -input $local_file86 $remote_uri86 -input $local_file87 $remote_uri87 -input $local_file88 $remote_uri88 -input $local_file89 $remote_uri89 -input $local_file90 $remote_uri90 -input $local_file91 $remote_uri91 -input $local_file92 $remote_uri92 -input $local_file93 $remote_uri93 -input $local_file94 $remote_uri94 -input $local_file95 $remote_uri95 -input $local_file96 $remote_uri96 -input $local_file97 $remote_uri97 -input $local_file98 $remote_uri98 -input $local_file99 $remote_uri99 -input $local_file100 $remote_uri100 -input $local_file101 $remote_uri101 -input $local_file102 $remote_uri102 -input $local_file103 $remote_uri103 -input $local_file104 $remote_uri104 -input $local_file105 $remote_uri105 -input $local_file106 $remote_uri106 -input $local_file107 $remote_uri107 -input $local_file108 $remote_uri108 -input $local_file109 $remote_uri109 -input $local_file110 $remote_uri100 -input $local_file111 $remote_uri111 -input $local_file112 $remote_uri112 -input $local_file113 $remote_uri113 -input $local_file114 $remote_uri114 -input $local_file115 $remote_uri115 -input $local_file116 $remote_uri116 -input $local_file117 $remote_uri117 -input $local_file118 $remote_uri118 -input $local_file119 $remote_uri119 -input $local_file120 $remote_uri120 -input $local_file121 $remote_uri121 -input $local_file122 $remote_uri122 -input $local_file123 $remote_uri123 -input $local_file124 $remote_uri124 -input $local_file125 $remote_uri125 -input $local_file126 $remote_uri126 -input $local_file127 $remote_uri127 -input $local_file128 $remote_uri128 -input $local_file129 $remote_uri129 -input $local_file130 $remote_uri130 -input $local_file131 $remote_uri131 -input $local_file132 $remote_uri132 -input $local_file133 $remote_uri133 -input $local_file134 $remote_uri134 -input $local_file135 $remote_uri135 -input $local_file136 $remote_uri136 -input $local_file137 $remote_uri137 -input $local_file138 $remote_uri138 -input $local_file139 $remote_uri139 -input $local_file140 $remote_uri140 -input $local_file141 $remote_uri141 -input $local_file142 $remote_uri142 -input $local_file143 $remote_uri143 -input $local_file144 $remote_uri144 -input $local_file145 $remote_uri145 -input $local_file146 $remote_uri146 -input $local_file147 $remote_uri147 -input $local_file148 $remote_uri148 -input $local_file149 $remote_uri149 -input $local_file150 $remote_uri150 -input $local_file151 $remote_uri151 -input $local_file152 $remote_uri152 -input $local_file153 $remote_uri153 -input $local_file154 $remote_uri154 -input $local_file155 $remote_uri155 -input $local_file156 $remote_uri156 -input $local_file157 $remote_uri157 -input $local_file158 $remote_uri158 -input $local_file159 $remote_uri159 -input $local_file160 $remote_uri160 -input $local_file161 $remote_uri161 -input $local_file162 $remote_uri162 -input $local_file163 $remote_uri163 -input $local_file164 $remote_uri164 -input $local_file165 $remote_uri165 $script $runnum

# # -input $local_file166 $remote_uri166 -input $local_file167 $remote_uri167 -input $local_file168 $remote_uri168 -input $local_file169 $remote_uri169 -input $local_file170 $remote_uri170 -input $local_file171 $remote_uri171 -input $local_file172 $remote_uri172 -input $local_file173 $remote_uri173 -input $local_file174 $remote_uri174 -input $local_file175 $remote_uri175 -input $local_file176 $remote_uri176 -input $local_file177 $remote_uri177 -input $local_file178 $remote_uri178 -input $local_file179 $remote_uri179 -input $local_file180 $remote_uri180 -input $local_file181 $remote_uri181 -input $local_file182 $remote_uri182 -input $local_file183 $remote_uri183 -input $local_file184 $remote_uri184 -input $local_file185 $remote_uri185 -input $local_file186 $remote_uri186 -input $local_file187 $remote_uri187 -input $local_file188 $remote_uri188 -input $local_file189 $remote_uri189 -input $local_file190 $remote_uri190 -input $local_file191 $remote_uri191 -input $local_file192 $remote_uri192 -input $local_file193 $remote_uri193 -input $local_file194 $remote_uri194 -input $local_file195 $remote_uri195 -input $local_file196 $remote_uri196 -input $local_file197 $remote_uri197 -input $local_file198 $remote_uri198 -input $local_file199 $remote_uri199 
     

   


# # for ((i=0; i<=$2; i++))
# # do
# #     fnameout_pattern='/farm_out/jboyd/bb_gmn_xtalk_'$runnum'_segment'$i'.out'
# #     #    sbatch --output=$fnameout_pattern run_GMN_sbatch_nohodo.sh $runnum -1 0 e1209019 $i 1
# #     jobname='bb_gmn_xtalk_'$runnum'_segment'$
    
# #     # look for first segment on cache disk:
# #     firstsegname='e1209019_'$runnum'.evio.0.0'
# #     mssfirst='mss:/mss/halla/sbs/raw/'$firstsegname
# #     cachefirst='/cache/mss/halla/sbs/raw/'$firstsegname
    
# #     eviofilename='e1209019_'$runnum'.evio.0.'$i
# #     mssfilename='mss:/mss/halla/sbs/raw/'$eviofilename
# #     cachefile='/cache/mss/halla/sbs/raw/'$eviofilename
    
# #     script='/w/halla-scshelf2102/sbs/jboyd/analysis/run_replay_here/run-xtalk.sh'
# #     #script='/work/halla/sbs/jboyd/SBS_OFFLINE/install/run_replay_here/run-gmn-replay.sh'
    
# #     testfilename='/mss/halla/sbs/raw/'$eviofilename
    
# #     #outfilename='match:e1209019_fullreplay_'$runnum'*seg'$i'*.root'
# #     outfilename='match:*.root'
# #     #logfilename='match:replay_gmn_'$runnum'*seg'$i'*.log'
# #     logfilename='match:*.log'

# #     if [ -f "$testfilename" ]; 
# #     then
# # 	echo 'Adding new swif2 job, runnum='$runnum', segment='$i 
    
# # 	if [ $i -gt 0 ]
# # 	then
# # 	    echo 'segment '$i' also requires first segment'
# # 	    swif2 add-job -workflow jboyd_GMN_analysis -partition production -name $jobname -cores 1 -disk 25GB -ram 1500MB -input $cachefile $mssfilename -input $cachefirst $mssfirst $script $runnum -1 0 e1209019 $i 1
# # 	else
# # 	    echo 'segment '$i' IS first segment'
# # 	    swif2 add-job -workflow jboyd_GMN_analysis -partition production -name $jobname -cores 1 -disk 25GB -ram 1500MB -input $cachefile $mssfilename $script $runnum -1 0 e1209019 $i 1
# # 	fi
# #     fi
# # done
 


# ##-output 'match:${SWIF_JOB_WORK_DIR}/*.root' '/volatile/halla/sbs/jboyd/Rootfiles/' -output 'match:${SWIF_JOB_WORK_DIR}/*.log' '/volatile/halla/sbs/jboyd/logs/'