# XTALK PERFORMANCE

The primary function of these files is to run through the metrics and performance of the Crosstalk output files.
The primary scripts are:

**submit-farm-strips.sh**, **farm_analyzer_strips.sh**, **xtalk_perf_multirun.C**, **xtalk_perf_histos.C**, and **xtalk_recon_efficiency.C**

# 1. submit-farm-strips.sh & farm_analyzer_strips.sh

  Description: This runs a replay of CACHE files that has Crosstalk ON or OFF and with some STRIP variables enabled. This also contains **histograms** needed by xtalk_perf_multirun.C.  These histograms are things like counts, number of tracks, etc. 
  
  ## Inputs: 
  
    Terminal line input: ./submit-farm-strips.sh $runnum $last_segment $first_segment

    Typical input for ALL segments: ./submit-farm-strips.sh 11449 33 0

  ## Pre-Requisites
  
    The .evio files need to be in the cache drive.
  
  ## Outputs:
    Creates: 
      e1209019_fullreplay_11449_stream0_seg5_5_xtalk_replay_Histogramming_WITH_corr_XTALK_OFF_STRIPS.root
      e1209019_fullreplay_11449_stream0_seg5_5_xtalk_replay_Histogramming_WITH_corr_XTALK_ON_STRIPS.root

    These are in **$farm_jobs/farm_replays/strips/RUNNUM/**
--------------------------------------------------------------------------------------------------------------------------------------

# 2. xtalk_perf_multirun.C

  This is a version of xtalk_perf_histos.C that lets you (not yet... still working on that) run multiple runs instead of just one.
  The idea will be to compile all of the strip, tracks, xtalk loops, etc, into one single histogram and table.
  
  **As a note, xtalk_perf_histos.C runs basically the same way it just can only take a single run number for input**
  
## Inputs:
  
  **Can select whether or not to "match segment counts". If this is true the number of files in the directory MUST match. This is improtant to ensure that we are comparing equal numbe of run & entries to each other. 

  ### Input: ./xtalk_perf_multirun.C 
  
  ** As of yet, there is no terminal line input. All of that is input in the file. Right now you enter the run numbers into a vector. Stuff like SBS config and magfield are automatically looked up from the beam_variables.h header.
  
## PRE-REQUISITES:  
  
  ### Requires STRIP variable crosstalk (Xtalk) files with XTALK_ON and XTALK_ON from **submit-farm-strips.sh**.
  
## Outputs: 

  Text file with histogram strips, tracking, crosstalk loop counts, etc. 
  PDF files with plots
  
  
--------------------------------------------------------------------------------------------------------------------------------------

# 3. calibrate.C

  Runs a vector of runnums and ouputs a file with lots of histograms that can be used to plot the various deltaplots for crosstalk recon efficiency. 
  
  This creates outputs that will be pasted into the beam_variables.h header for the cut/fit lookups
  
  
## Inputs:
  The running configurations depend on whether or not all cut/fit parameters DO or DO NOT exist in beam_variables.h. If you are running it to fill in the table then you need to run a multi-step process:
  
  ### **Full Process (No previous fit/cut parameters)**
  
  ### a) Crosstalk ON/OFF For each of the following:
  
    ### b) calc_W = false
        This provides cuts and fits for all of the variables that ARE NOT W. So, E/p, PS, PS + SH, HCal_clust, etc.

    ### c) calc_W = true
        This provides cuts for the W and W^2 including all previous parameters from previous !calc_W
      
  ### Input: 
  
  Simply: ./calibrate.C --> No terminal line inputs. Everthing is specified inside. 
  
## PRE-REQUISITES:  
  
  ### Three selections include:
    
    -- No Crosstalk (!crosstalk): Work/Volatile drive for Pass 0 and Pass replays
    
    -- Crosstalk On:
      --- Strips: /lustre19/expphy/volatile/halla/sbs/jboyd/swif_output/xtalk/farm_replays/strips/$RUNNUM$/$XTALK_ONOFF$
      
      --- No Strips: /lustre19/expphy/volatile/halla/sbs/jboyd/swif_output/xtalk/farm_replays
      
  ### /w/halla-scshelf2102/sbs/jboyd/analysis/gmn/deltaplots/rootfiles/*_XTALK_$ON/OFF$_full_calibration.root
  
  
## Outputs: 
  
  
--------------------------------------------------------------------------------------------------------------------------------------

# 4. xtalk_recon_efficiency.C

  **ONLY RUNS ON A SINGLE RUN RIGHT NOW** 
  This script creates the dxdy plots for runs. It will create it 
  
  
## Inputs:

  ### Input: 
  
## PRE-REQUISITES:  
  
  ### /w/halla-scshelf2102/sbs/jboyd/analysis/gmn/deltaplots/rootfiles/*_XTALK_$ON/OFF$_full_calibration.root
  
## Outputs: 
  
  
--------------------------------------------------------------------------------------------------------------------------------------


# 5. count_and_organize_files.C

  
 
## Inputs:

  ### Input: 
  
## PRE-REQUISITES:  
  
  ### 
  
## Outputs: 
  
  
--------------------------------------------------------------------------------------------------------------------------------------
