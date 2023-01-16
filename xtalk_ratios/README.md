This directory is for the scripts that will generate a Neighbor Ratio plot for all APVs in a run.

The script relies on the prime GMn replay. No strip variables are need in output file. The strip comparison is performed on the replay-level.
The only criteria is that there be two sets of files: XTALK_ON and XTALK_OFF.

Another interesting and annoying aspect of this is that due to a bug in one of the codes, a middle step requires analyzing the replay files into ALL_MODULES AND LAST_MODULE. For some reason, in the ALL_MODULES version the very last module is not properly analyzed. So, we must run a separate analysis where we only look at the final, LAST, module.

# 1. MAIN REPLAY

  ## submit-xtalk-jobs.sh

  Description: This produces a controlled set of replays to be used down the line for the Crosstalk Ratio analysis. These replays contains a certain set of strip variables needed by the crosstalk analysis. 
  These variables are: bb.gem.m#.strip.--> istrip, ADCsum, isampmax, IsU, IsV.
  This also pulls from DB_XTALK and therefore pulls some certain flags and criteria from that db file.

  ## run-xtalk.sh

   Called by submit-xtalk-jobs.sh
   
  ## replay_xtalk_farm_jobs.C
   
  -------------------------------------------------------------------
  ## Inputs:
  
  ### submit-xtalk-jobs.sh $runnum $maxsegments
  
  ### Pre-requisites:
  
    Must have files in /cache/halla/sbs/ 

  ## Outputs:
  
  ### "$xtalk_out/RUNNUM/e1209019_fullreplay_$RUNNUM$_stream0_seg#_#_xtalk_replay_Histogramming_WITH_corr.root"
  
  Directory:  $xtalk_out/RUNNUM = /lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/RUNNUM

--------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------
   
# 2. XTALK ANALYSIS

  ## submit-xtalk-analysis.sh
  
  Description: This takes as input the outputs from "submit-xtalk-jobs.sh". This will run an event-by-event crosstalk analysis. This analysis is the neighboring ADC channel ratio analysis.
  
  This one produces two outputs: **ALL_MODULES** and **LAST_MODULE**
  Due to a bug in the script, when analyzing multiple consecutive modules, the last module is not read correctly. So, as a quick fix, we run on all modules and expect the last one to be faulty. In parallel, we also run just the last module. Then we use both of these as input. 
  
  ## xtalk-analysis.sh
    Called by submit-xtalk-analysis.sh  
    
  ## xtalk_by_events_farm.C
  
 -------------------------------------------------------------------
 
  ## Inputs:
  
  ### submit-xtalk-analysis.sh $runnum $allsegs $lastevent  $firstsegment $lastmodule $lastsegment
  
  $runnum: The run number
  
  $allsegs: This defines which segments to analyze. 
    -1 --> Analyze all segments
    0->inf --> Analyze to that segment
   
  $lastevent: This is the last event to analyze:
    -1 --> All events
    0->inf --> Analyze to that event
    
  $firstsegment: Which segment to start on, typically 0.
  
  $last module:
    0 --> ALL MODULES. You are saying you DO NOT want to analyze only the last module
    1 --> ONLY analyze the last module.
    
  $lastsegment
    ___ LEAVE BLANK TO ANALYZE ALL SEGMENTS
    
  This analysis will require to steps. They are typically:
  
      ./submit-xtalk-analysis.sh $runnum -1 -1 0 0
      ./submit-xtalk-analysis.sh $runnum -1 -1 0 1
      
  
  ### Pre-requisites:
  
     Requires outputs from ./submit-xtalk-jobs: e1209019_fullreplay_$RUNNUM$_stream0_seg#_#_xtalk_replay_Histogramming_WITH_corr.root
     
  ## Outputs:
  
  ### "$xtalk_out/xtalk_by_events/RUNNUM/$RUNNUM$_xtalk_ratios_0_thru0_events_U1_V1_seg_#_ALL.root"
  ### "$xtalk_out/xtalk_by_events/RUNNUM/$RUNNUM$_xtalk_ratios_0_thru0_events_U1_V1_seg_#_last_module_only.root"
  
  Directory:  $xtalk_out/RUNNUM = /lustre19/expphy/volatile/halla/sbs/jboyd/Rootfiles/xtalk/RUNNUM
  
  --------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------
   
# 3. Fits (farm_fits)

  ## submit-farm-fits.sh
  
  Description: Takes the crosstalk analysis outputs (ALL and last module only) and plots/fits them for each APV. Outputs a GEM xtalk DB file along with corresponding plots.
  
  This one produces two outputs: DB file that ends up getting called up like a pedestal file. PDF plots of each module.
  
  ## farm_fits.sh
    Called by submit-xtalk-analysis.sh  
    
  ## Main code: /work/halla/sbs/jboyd/xtalk/fit_xtalk/xtalk_fits_farm.C
  
 -------------------------------------------------------------------
 
  ## Inputs:
  
  ### REQUIRES that the outputs for ALL and last_mod_only be merged and placed into a /merged/ subfolder. File naming nomenclature is:
    --> $RUNNUM$_xtalk_ratios_ALL_merged.root and $RUNNUM$_xtalk_ratios_LAST_merged.root
  
  ### submit-farm-fits.sh $runnum $adccut
  
  $runnum: The run number
  
  $adccut The threshold for the numerator (larger) ADC for it to be considered for Crosstalk Correction.
  
  For ADCsum: 1500 
  
  This value will get propagated to the output DB file and should be recorded therefore through.
      
  ### Pre-requisites:
  
     Requires outputs from ./submit-xtalk-analysis.sh: ..._ALL.root & ..._last_module_only.root
     
  ## Outputs:
  
  ### DB File: db_xtalk-bb_gem_run$RUNNUM$.dat
    Contains the Crosstalk ratios for all APVs in the configuration. 
    **Must cross-reference with the PDF plots**
  ### PDF Files, 1 for U, and 1 for V: xtalk_ratios_$RUNNUM$_$UV$_ADCcut$ADCCUT$_ALL.pdf"
  
  Directory: 
  **DB FILE**: "/w/halla-scshelf2102/sbs/jboyd/xtalk/xtalk_ratio_outputs/txt/$RUNNUM$/ADCcut$ADCCUT/"
  **Plots**: "/w/halla-scshelf2102/sbs/jboyd/xtalk/xtalk_ratio_outputs/plots/$RUNNUM$/ADCcut$ADCCUT/"
  
